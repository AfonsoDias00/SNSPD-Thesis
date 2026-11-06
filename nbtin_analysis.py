#!/usr/bin/env python3
"""
nbtin_analysis.py -- the gamma = 1 kinetic results against the "nbtin" rerun (gamma = 10, NbTiN hotspot times,
settled start; run_nbtin.sh), for every device and both criteria, from the stored maps (no simulation).

Everything is in units of I_dep of the narrowest free cross-section, which does not depend on gamma or the hotspot.
Per device and criterion: dark I_sw, the efficiency quantiles i_10 / i_50 / i_95 / i_max of the interpolated map,
the site-by-site change, the gap between the two bounds, the latency of the switching run at threshold, flags, and
the selectivity (probe sites at f = 0.6 and 0.35).  Then the go / no-go tests of the hole study under each model.

  python nbtin_analysis.py            -> results/nbtin_analysis.txt, results/nbtin_efficiency.png,
                                         results/nbtin_sites.png, results/nbtin_latency.png
Works on partial results: devices or criteria without rows are skipped.
"""
from __future__ import annotations

import numpy as np

from snspd import analysis as an
from snspd import geometry as geo
from snspd import plotting as pl
from snspd import store
from snspd.units import FILM

MODELS = {"gamma 1": "kinetic", "nbtin": "kinetic-nbtin"}
CRITERIA = ("trigger", "sustained")
T0 = 20.0
GO_MARGIN = 0.03


def devices():
    """(label, geometry, csv stem path builder, isw json path builder)."""
    out = []
    for d, keys in (("results/hitmap", ("Plain_Wire", "Long_Zipper")),
                    ("results/slots", ("Plain_Wire_2um", "Slot_110", "Slot_200", "Slot_400"))):
        for k in keys:
            out.append((k, geo.get(k), lambda tag, d=d: f"{d}/hitmap_{tag}.csv", lambda tag, d=d: f"{d}/isw_{tag}.json"))
    for w in (30, 56, 100):
        d = "results/narrow_baseline"
        out.append((f"{w} nm", geo.plain_strip(w), lambda tag, w=w: f"{d}/narrow{w}nm_hitmap_{tag}.csv",
                    lambda tag, w=w: f"{d}/narrow{w}nm_isw_{tag}.json"))
    return out


def tag(model, criterion):
    return f"{model}_T400_m1.4" + ("_trigger" if criterion == "trigger" else "")


def load(label, g, csv_of, isw_of, model, criterion):
    """dict(sites={(x, y): row}, vals=pixel thresholds [I_dep], isw=I_sw/I_dep, probes={(x, y, f): i}) or None."""
    t = tag(model, criterion)
    rows = [r for r in store.load_rows(csv_of(t)) if r["Geometry"] == g.key]
    if not rows:
        return None
    i_dep = FILM.I_dep_uA(g.free_width_min_nm())
    isw_j = store.load_json(isw_of(t), {}).get(g.key)
    isw = float(isw_j["I_sw"]) if isw_j else float(rows[0]["I_sw_uA"])
    sites = {(float(r["x_nm"]), float(r["y_nm"])): r for r in rows if abs(float(r["energy_f"]) - 1) < 1e-6}
    samples = an.photon_rows(rows, g.key)
    vals = an.pixel_values(an.unit_map(g, samples)[2]) * isw / i_dep if samples else None
    probes = {(float(r["x_nm"]), float(r["y_nm"]), float(r["energy_f"])):
              (isw if r["flag"] == "at_wall" else float(r["I_det_uA"])) / i_dep for r in rows}
    return dict(sites=sites, vals=vals, isw=isw / i_dep, i_dep=i_dep, probes=probes, n=len(sites))


def thr(r, i_dep):
    """Threshold of one row [I_dep]; a site that never switched counts as I_sw."""
    return float(r["I_sw_uA" if r["flag"] == "at_wall" else "I_det_uA"]) / i_dep


def q(vals, p):
    return float(np.quantile(vals, p)) if vals is not None and len(vals) else float("nan")


def selectivity(D):
    """Median over probe sites of i(f) - i(1) at f = 0.6 and 0.35 [I_dep]."""
    out = {}
    for f in (0.6, 0.35):
        d = [D["probes"][(x, y, f)] - D["probes"][(x, y, 1.0)] for (x, y, ff) in D["probes"]
             if ff == f and (x, y, 1.0) in D["probes"]]
        out[f] = float(np.median(d)) if d else float("nan")
    return out


def main():
    data = {}
    for label, g, csv_of, isw_of in devices():
        for m, model in MODELS.items():
            for c in CRITERIA:
                D = load(label, g, csv_of, isw_of, model, c)
                if D is not None:
                    data[(label, m, c)] = D
    labels = [d[0] for d in devices()]
    geoms = {d[0]: d[1] for d in devices()}
    L = ["GAMMA 1 vs NBTIN (gamma 10, NbTiN hotspot times D 7.1 nm^2/tau0 and escape 530 tau0, settled start)",
         "kinetic hotspot f = 1, 400 tau0, 1.4 xi, pinned meshes; thresholds in I_dep of the narrowest free width", ""]

    # --- quantiles -------------------------------------------------------------------------------------------
    for c in CRITERIA:
        L += [f"=== {c} {'(lower bound)' if c == 'trigger' else '(upper bound)'} ===",
              f"{'device':<16}{'sites':>6}{'I_sw':>15}{'i_10':>15}{'i_50':>15}{'i_95':>15}{'i_max':>15}",
              f"{'':<16}{'':>6}" + "   g1 -> nbtin " * 5]
        for lab in labels:
            a, b = data.get((lab, "gamma 1", c)), data.get((lab, "nbtin", c))
            if a is None and b is None:
                continue
            cells = []
            for key in ("isw", 0.10, 0.50, 0.95, 1.00):
                va = (a["isw"] if key == "isw" else q(a["vals"], key)) if a else float("nan")
                vb = (b["isw"] if key == "isw" else q(b["vals"], key)) if b else float("nan")
                cells.append(f"  {va:5.3f}->{vb:5.3f}" if b else f"  {va:5.3f}->  -  ")
            n = f"{b['n'] if b else 0}/{a['n'] if a else 0}"
            L.append(f"{lab:<16}{n:>6}" + "".join(cells))
        L.append("")

    # --- site by site ----------------------------------------------------------------------------------------
    L += ["=== site by site: nbtin - gamma 1 [I_dep] at the sites present in both (f = 1) ===",
          f"{'device':<16}{'criterion':<11}{'n':>4}{'mean':>8}{'min':>8}{'max':>8}   largest changes"]
    site_pairs = {}
    for lab in labels:
        for c in CRITERIA:
            a, b = data.get((lab, "gamma 1", c)), data.get((lab, "nbtin", c))
            if not a or not b:
                continue
            common = sorted(set(a["sites"]) & set(b["sites"]))
            if not common:
                continue
            d = {s: thr(b["sites"][s], b["i_dep"]) - thr(a["sites"][s], a["i_dep"]) for s in common}
            site_pairs[(lab, c)] = [(thr(a["sites"][s], a["i_dep"]), thr(b["sites"][s], b["i_dep"])) for s in common]
            v = np.array(list(d.values()))
            ext = sorted(d, key=lambda s: -abs(d[s]))[:3]
            L.append(f"{lab:<16}{c:<11}{len(v):4d}{v.mean():+8.3f}{v.min():+8.3f}{v.max():+8.3f}   "
                     + ", ".join(f"({s[0]:g},{s[1]:g}) {d[s]:+.3f}" for s in ext))
    L.append("")

    # --- gap between the bounds --------------------------------------------------------------------------------
    L += ["=== gap between the bounds, sustained - trigger, per site [I_dep]: median (max) ===",
          f"{'device':<16}{'gamma 1':>18}{'nbtin':>18}"]
    for lab in labels:
        cells = []
        for m in MODELS:
            lo, hi = data.get((lab, m, "trigger")), data.get((lab, m, "sustained"))
            if not lo or not hi:
                cells.append(f"{'-':>18}")
                continue
            common = sorted(set(lo["sites"]) & set(hi["sites"]))
            gap = np.array([thr(hi["sites"][s], hi["i_dep"]) - thr(lo["sites"][s], lo["i_dep"]) for s in common])
            cells.append(f"{np.median(gap):10.3f} ({gap.max():.3f})" if len(gap) else f"{'-':>18}")
        L.append(f"{lab:<16}" + "".join(cells))
    L.append("")

    # --- latency, branches, flags ------------------------------------------------------------------------------
    L += ["=== the switching run at threshold: latency after the photon [tau0] (median, 90 %, max), voltage branch, "
          "flags ===", f"{'device':<16}{'criterion':<11}{'model':<9}{'latency median / 90% / max':>30}"
          f"{'late (>300)':>12}{'V_peak/V_trail':>16}{'flags':>18}"]
    lat = {}
    for lab in labels:
        for c in CRITERIA:
            for m in MODELS:
                D = data.get((lab, m, c))
                if not D:
                    continue
                rs = list(D["sites"].values())
                t = np.array([float(r["t_trigger"]) - T0 for r in rs if r["t_trigger"] not in ("", "None")])
                lat[(lab, c, m)] = t
                br = [r["branch"] for r in rs]
                fl = [r["flag"] for r in rs if r["flag"]]
                stat = f"{np.median(t):8.0f} / {np.quantile(t, 0.9):5.0f} / {t.max():5.0f}" if len(t) else "-"
                L.append(f"{lab:<16}{c:<11}{m:<9}{stat:>30}{int((t > 300 - T0).sum()):>12}"
                         f"{br.count('V_peak'):>8}/{br.count('V_trail'):<7}"
                         f"{(', '.join(f'{x} {fl.count(x)}' for x in sorted(set(fl))) or '-'):>18}")
    L.append("")

    # --- selectivity -------------------------------------------------------------------------------------------
    L += ["=== selectivity: median over probe sites of i(f) - i(1) [I_dep]; weaker photons need more bias ===",
          f"{'device':<16}{'criterion':<11}{'f = 0.6  g1 -> nbtin':>24}{'f = 0.35  g1 -> nbtin':>24}"]
    for lab in labels:
        for c in CRITERIA:
            a, b = data.get((lab, "gamma 1", c)), data.get((lab, "nbtin", c))
            if not a or not b:
                continue
            sa, sb = selectivity(a), selectivity(b)
            L.append(f"{lab:<16}{c:<11}{sa[0.6]:>12.3f} -> {sb[0.6]:6.3f}{sa[0.35]:>14.3f} -> {sb[0.35]:6.3f}")
    L.append("")

    # --- go / no-go ----------------------------------------------------------------------------------------------
    L += [f"=== go / no-go: does the hole device beat its plain reference by {GO_MARGIN} I_dep? (i_50 onset, "
          "i_95 plateau, i_max) ==="]
    for c in CRITERIA:
        for m in MODELS:
            for ref, holes in (("Plain_Wire", ("Long_Zipper",)), ("Plain_Wire_2um", ("Slot_110", "Slot_200", "Slot_400")),
                               ("56 nm", ("Long_Zipper",))):
                R = data.get((ref, m, c))
                for h in holes:
                    H = data.get((h, m, c))
                    if not R or not H:
                        continue
                    d = {k: q(R["vals"], k) - q(H["vals"], k) for k in (0.5, 0.95, 1.0)}
                    verdict = "GO" if d[0.95] >= GO_MARGIN else "no-go"
                    L.append(f"  {c:<10}{m:<9}{h:<12} vs {ref:<15} gain i_50 {d[0.5]:+.3f}  i_95 {d[0.95]:+.3f}  "
                             f"i_max {d[1.0]:+.3f}   -> {verdict} (plateau)")
    text = "\n".join(L)
    print(text)
    with open("results/nbtin_analysis.txt", "w") as fh:
        fh.write(text + "\n")
    figures(data, geoms, site_pairs, lat)


def figures(data, geoms, site_pairs, lat):
    import matplotlib.pyplot as plt
    pl.apply_style()
    # efficiency curves: rows = 0.8 um family / 2 um family, columns = criteria
    fam = [("0.8 um strips", ["Plain_Wire", "Long_Zipper", "30 nm", "56 nm", "100 nm"]),
           ("2 um strips", ["Plain_Wire_2um", "Slot_110", "Slot_200", "Slot_400"])]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4), sharex=True, sharey=True, layout="constrained")
    ib = np.linspace(0.2, 1.02, 1200)
    for i, (fname, labs) in enumerate(fam):
        for j, c in enumerate(CRITERIA):
            ax = axes[i, j]
            for lab in labs:
                col = geoms[lab].color
                for m, ls, lw in (("gamma 1", (0, (3, 2)), 1.2), ("nbtin", "-", 2.0)):
                    D = data.get((lab, m, c))
                    if D is None or D["vals"] is None:
                        continue
                    ax.plot(ib, an.eta_curve(D["vals"], ib), color=col, ls=ls, lw=lw,
                            label=f"{lab}" if m == "nbtin" else None)
            ax.set_title(f"{fname}: {'voltage trigger (lower bound)' if c == 'trigger' else 'sustained (upper bound)'}",
                         fontsize=10)
            if i == 1:
                ax.set_xlabel("I_bias / I_dep (narrowest cross-section)")
            if j == 0:
                ax.set_ylabel("share of the unit cell that detects [%]")
            ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.suptitle("Detection efficiency: gamma 1 (dashed) vs gamma 10 + NbTiN hotspot times (solid)", fontsize=11)
    fig.savefig("results/nbtin_efficiency.png", dpi=200)
    plt.close(fig)

    # site by site
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharex=True, sharey=True, layout="constrained")
    for ax, c in zip(axes, CRITERIA):
        for (lab, cc), pairs in site_pairs.items():
            if cc != c:
                continue
            p = np.array(pairs)
            ax.plot(p[:, 0], p[:, 1], "o", ms=5, color=geoms[lab].color, mec=pl.SURFACE, mew=0.6, label=lab, alpha=0.9)
        ax.plot([0.2, 1.02], [0.2, 1.02], color=pl.INK2, lw=1, ls=":")
        ax.set_title("voltage trigger" if c == "trigger" else "sustained", fontsize=10)
        ax.set_xlabel("threshold, gamma 1 [I_dep]")
        ax.legend(loc="upper left", fontsize=8, frameon=False)
    axes[0].set_ylabel("threshold, gamma 10 + NbTiN times [I_dep]")
    fig.suptitle("Site by site: every photon threshold before and after (dotted: unchanged)", fontsize=11)
    fig.savefig("results/nbtin_sites.png", dpi=200)
    plt.close(fig)

    # latency at threshold
    labs = [lab for lab in geoms if any((lab, c, m) in lat for c in CRITERIA for m in MODELS)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.45 * len(labs) + 1.8), sharey=True, layout="constrained")
    for ax, c in zip(axes, CRITERIA):
        for k, lab in enumerate(labs):
            y = len(labs) - 1 - k
            for m, dy, col in (("gamma 1", 0.15, "#b5b3ac"), ("nbtin", -0.15, "#2a78d6")):
                t = lat.get((lab, c, m))
                if t is None or not len(t):
                    continue
                ax.plot(t, np.full(len(t), y + dy), "|", ms=8, color=col, mew=1.4)
        ax.set_yticks(range(len(labs)), labs[::-1], fontsize=9)
        ax.set_xlabel("latency of the switching run at threshold, tau0 after the photon")
        ax.set_title("voltage trigger" if c == "trigger" else "sustained", fontsize=10)
        ax.axvline(400 - T0, color=pl.INK2, lw=0.8, ls=":")
    fig.legend(handles=[plt.Line2D([], [], color="#b5b3ac", marker="|", ls="", ms=8, mew=1.4, label="gamma 1"),
                        plt.Line2D([], [], color="#2a78d6", marker="|", ls="", ms=8, mew=1.4, label="nbtin")],
               loc="outside lower center", ncol=2, frameon=False, fontsize=8.5)
    fig.savefig("results/nbtin_latency.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
