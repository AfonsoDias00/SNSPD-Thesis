#!/usr/bin/env python3
"""
energy_analysis.py -- efficiency against photon energy, from the thresholds of energy_scan.py (11 Nov 2026).

For each geometry and wavelength:
  * the hard-step efficiency over the unit cell (thresholds interpolated as in snspd/analysis.py) and its quantiles,
    with i_90 = I_det as Simon et al. define it (90 % of the plateau of the count-rate curve);
  * the same with the energy fluctuations of the deposited energy added analytically: at each pixel the threshold
    i_det(E) is interpolated in log E between the simulated energies, a deposited energy E ~ N(E0, (s E0)^2) makes
    the pixel count with probability P(E > E_th(i)), and the curve is the area average.  s = 0 is the hard step.

Only sites simulated at every wavelength enter the energy comparison (for the double column: the medium grid).  At
1550 nm the double column's full grid, when complete, is compared with its medium-grid subset.

  python energy_analysis.py                      writes results/energy/summary.txt and figures
"""
from __future__ import annotations

import os

import numpy as np
from scipy.special import ndtr

from snspd import analysis as an
from snspd import geometry as geo
from snspd import store

WAVELENGTHS = (1550.0, 2500.0, 4000.0)
HC = 1239.84198
GEOMS = ("Plain_Wire", "Double_400")
SIGMAS = (0.0, 0.1, 0.2)              # relative spread of the deposited energy
OUT = "results/energy"
COLORS = {"Plain_Wire": "#1f5fa6", "Double_400": "#c2410c"}
LS = {1550.0: "-", 2500.0: "--", 4000.0: ":"}


def path(w):
    return f"results/thresholds/thr_n750_w{w:g}_m2.csv"


def dark_of(rows, key):
    for r in rows:
        if r["Geometry"] == key and r["x_nm"] == "":
            return float(r["i_det"])
    return None


def samples(rows, key, dark):
    """(x, y) -> threshold; 'at_wall' and 'dark' rows become the dark switching current (the device counts there
    anyway, photon or not)."""
    out = {}
    for r in rows:
        if r["Geometry"] != key or r["x_nm"] == "":
            continue
        v = float(r["i_det"])
        if r["flag"] in ("at_wall", "dark") and dark is not None:
            v = dark
        out[(float(r["x_nm"]), float(r["y_nm"]))] = v
    return out


def maps_for(g, by_w, sites):
    """Threshold map per wavelength on the same pixels, from the given sites only."""
    Z = {}
    for w, s in by_w.items():
        pts = [(x, y, s[(x, y)]) for x, y in sites]
        X, Y, z = an.unit_map(g, pts)
        Z[w] = z
    return X, Y, Z


def smooth_eta(Z, ib, sigma):
    """Area-averaged count probability with a Gaussian spread sigma (relative) of the deposited energy."""
    ws = sorted(Z)
    E = np.array([HC / w for w in ws])
    mask = np.isfinite(Z[ws[0]])
    I = np.stack([Z[w][mask] for w in ws])               # (n_energy, n_pix): i_det at each energy
    if sigma == 0:
        return {w: an.eta_curve(np.sort(Z[w][mask]), ib) for w in ws}
    lnE = np.log(E)
    o = np.argsort(lnE)
    lnE, I = lnE[o], I[o]
    # the threshold cannot rise with the photon energy: enforce it (removes bisection noise), then make it strictly
    # decreasing so a flat stretch becomes a near-vertical one in ln E
    I = np.minimum.accumulate(I, axis=0) - 1e-6 * np.arange(len(lnE))[:, None]
    out = {}
    for w in ws:
        E0 = HC / w
        res = np.zeros(len(ib))
        for k, i in enumerate(ib):
            # ln E_th(i) per pixel, piecewise linear in i between the simulated energies, extrapolated with the end
            # segments: segment j covers I[j+1] <= i (j = 0 also everything above I[0]; the last also everything
            # below the lowest threshold)
            lth = np.empty(I.shape[1])
            seg = np.zeros(I.shape[1], int)
            for j in range(1, len(lnE) - 1):
                seg[i < I[j]] = j
            a, b = np.take_along_axis(I, seg[None], 0)[0], np.take_along_axis(I, seg[None] + 1, 0)[0]
            la, lb = lnE[seg], lnE[seg + 1]
            lth = la + (a - i) / (a - b) * (lb - la)
            Eth = np.exp(np.clip(lth, -50, 50))
            res[k] = 100.0 * np.mean(1.0 - ndtr((Eth - E0) / (sigma * E0)))
        out[w] = res
    return out


def quant(ib, eta, q):
    k = np.argmax(eta >= q)
    return float(ib[k]) if eta[k] >= q else float("nan")


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = {w: store.load_rows(path(w)) for w in WAVELENGTHS}
    dark = {k: dark_of(rows[1550.0], k) for k in GEOMS}
    L = ["ENERGY SCAN  (2.0 xi, circuit 750 squares, bath 0.3 T_c; bias in I_dep of the narrowest free width)", ""]
    ib = np.linspace(0.5, 1.0, 501)
    curves, sitetab = {}, {}
    for key in GEOMS:
        g = geo.get(key)
        by_w = {w: samples(rows[w], key, dark[key]) for w in WAVELENGTHS}
        common = sorted(set.intersection(*(set(s) for s in by_w.values())), key=lambda p: (p[1], p[0]))
        L.append(f"{g.label}: dark switching {an.fmt(dark[key])}; sites at every wavelength: {len(common)} "
                 f"(simulated: " + ", ".join(f"{w:g} nm {len(by_w[w])}" for w in WAVELENGTHS) + ")")
        sitetab[key] = (common, by_w)
        L.append(f"  {'site':>12}  " + "  ".join(f"{w:>7g}" for w in WAVELENGTHS))
        for s in common:
            L.append(f"  {'(%g, %g)' % s:>12}  " + "  ".join(f"{by_w[w][s]:7.3f}" for w in WAVELENGTHS))
        if len(common) < 3:
            L.append("  (too few common sites for curves)\n")
            continue
        X, Y, Z = maps_for(g, by_w, common)
        curves[key] = {}
        for sig in SIGMAS:
            eta = smooth_eta(Z, ib, sig)
            curves[key][sig] = eta
            for w in WAVELENGTHS:
                e = eta[w]
                L.append(f"  sigma {sig:.1f}  {w:>5g} nm ({HC / w:.2f} eV): i_10 {quant(ib, e, 10):.3f}  i_50 "
                         f"{quant(ib, e, 50):.3f}  i_90 {quant(ib, e, 90):.3f}  i_99 {quant(ib, e, 99):.3f}  "
                         f"eta(0.80) {np.interp(0.80, ib, e):5.1f} %  eta(0.90) {np.interp(0.90, ib, e):5.1f} %")
        L.append("")
    # medium vs full grid at 1550
    g = geo.get("Double_400")
    full = samples(rows[1550.0], "Double_400", dark["Double_400"])
    fg = g.hit_grid("full")
    if all(s in full for s in fg):
        from energy_scan import MEDIUM
        med = [(float(x), float(y)) for x, y in MEDIUM]
        for name, sites in (("full 63", fg), ("medium 23", med)):
            _, _, z = an.unit_map(g, [(x, y, full[(x, y)]) for x, y in sites])
            m = an.metrics(an.pixel_values(z))
            L.append(f"Double_400 at 1550 nm, {name} sites: i_10 {m['i_10']:.3f} i_50 {m['i_50']:.3f} "
                     f"i_90 {m['i_90']:.3f} i_99 {m['i_99']:.3f} i_max {m['i_max']:.3f}")
    else:
        L.append(f"Double_400 full grid at 1550 nm: {sum(s in full for s in fg)} of {len(fg)} sites so far")
    text = "\n".join(L)
    print(text)
    with open(f"{OUT}/summary.txt", "w") as fh:
        fh.write(text + "\n")
    plot(ib, curves, sitetab, dark)
    absorbed_check(rows, dark, ib)


ABS_FILE = "results/thresholds/thr_n750_w1550_m2_abs0.5.csv"


def absorbed_check(rows, dark, ib):
    """1550 nm with half of h nu deposited (0.40 eV): measured thresholds against the ln-E interpolation of the
    0.31 / 0.50 / 0.80 eV runs, and the efficiency curves next to the full-deposit ones."""
    if not os.path.exists(ABS_FILE):
        return
    va = store.load_rows(ABS_FILE)
    E = 0.5 * HC / 1550.0
    lnE = np.log([HC / w for w in (4000.0, 2500.0, 1550.0)])
    L = [f"DEPOSITED FRACTION 0.5 AT 1550 nm ({E:.3f} eV)  vs  full deposit, and vs ln-E interpolation of the "
         f"other runs", ""]
    cur = {}
    for key in GEOMS:
        g = geo.get(key)
        sa = samples(va, key, dark[key])
        by_w = {w: samples(rows[w], key, dark[key]) for w in WAVELENGTHS}
        sites = sorted((s for s in sa if all(s in by_w[w] for w in WAVELENGTHS)), key=lambda p: (p[1], p[0]))
        L.append(f"{g.label}: {len(sites)} sites")
        L.append(f"  {'site':>12}  {'0.80 eV':>7}  {'0.40 eV':>7}  {'interp':>7}  {'diff':>6}")
        diffs = []
        for s in sites:
            pred = float(np.interp(np.log(E), lnE, [by_w[w][s] for w in (4000.0, 2500.0, 1550.0)]))
            diffs.append(sa[s] - pred)
            L.append(f"  {'(%g, %g)' % s:>12}  {by_w[1550.0][s]:7.3f}  {sa[s]:7.3f}  {pred:7.3f}  {diffs[-1]:+6.3f}")
        if diffs:
            L.append(f"  measured - interpolated: mean {np.mean(diffs):+.3f}, max |.| {np.max(np.abs(diffs)):.3f}")
        if len(sites) >= 3:
            c = {}
            for name, src in (("absorbed 1.0", by_w[1550.0]), ("absorbed 0.5", sa)):
                _, _, z = an.unit_map(g, [(x, y, src[(x, y)]) for x, y in sites])
                e = an.eta_curve(np.sort(z[np.isfinite(z)]), ib)
                c[name] = e
                L.append(f"  {name}: i_10 {quant(ib, e, 10):.3f}  i_50 {quant(ib, e, 50):.3f}  i_90 "
                         f"{quant(ib, e, 90):.3f}  i_99 {quant(ib, e, 99):.3f}  eta(0.80) "
                         f"{np.interp(0.80, ib, e):5.1f} %  eta(0.90) {np.interp(0.90, ib, e):5.1f} %")
            cur[key] = c
        L.append("")
    text = "\n".join(L)
    print(text)
    with open(f"{OUT}/summary_abs0.5.txt", "w") as fh:
        fh.write(text + "\n")
    if cur:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        for key, c in cur.items():
            for name, ls in (("absorbed 1.0", "-"), ("absorbed 0.5", "--")):
                ax.plot(ib, c[name], color=COLORS[key], ls=ls, lw=2, label=f"{geo.get(key).label}, {name}")
            ax.axvline(dark[key], color=COLORS[key], lw=1, alpha=0.4)
        ax.set_xlim(0.6, 1.0)
        ax.set_xlabel("bias  [I_dep of the narrowest free width]")
        ax.set_ylabel("share of the cell that counts  [%]")
        ax.set_title("1550 nm, deposited fraction of h nu (hard step, common sites)", fontsize=10)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(f"{OUT}/efficiency_abs0.5.png", dpi=150)
        plt.close(fig)


def plot(ib, curves, sitetab, dark):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if curves:
        fig, axes = plt.subplots(1, len(SIGMAS), figsize=(4.6 * len(SIGMAS), 4.8), sharey=True)
        for ax, sig in zip(np.atleast_1d(axes), SIGMAS):
            for key, c in curves.items():
                for w in WAVELENGTHS:
                    ax.plot(ib, c[sig][w], color=COLORS[key], ls=LS[w], lw=2,
                            label=f"{geo.get(key).label}, {w:g} nm")
                if dark.get(key):
                    ax.axvline(dark[key], color=COLORS[key], lw=1, alpha=0.4)
            ax.set_title("hard step" if sig == 0 else f"energy spread {sig:.0%}", fontsize=10)
            ax.set_xlabel("bias  [I_dep of the narrowest free width]")
            ax.grid(alpha=0.25)
            ax.set_xlim(0.6, 1.0)
        np.atleast_1d(axes)[0].set_ylabel("share of the cell that counts  [%]")
        h, lab = np.atleast_1d(axes)[0].get_legend_handles_labels()
        fig.legend(h, lab, fontsize=8, loc="lower center", ncol=2, frameon=False)
        fig.tight_layout(rect=(0, 0.16, 1, 1))
        fig.savefig(f"{OUT}/efficiency_vs_energy.png", dpi=150)
        plt.close(fig)
    # threshold vs photon energy per site
    fig, axes = plt.subplots(1, len(sitetab), figsize=(5.2 * len(sitetab), 4.0), sharey=True)
    for ax, (key, (common, by_w)) in zip(np.atleast_1d(axes), sitetab.items()):
        E = [HC / w for w in WAVELENGTHS]
        for s in common:
            ax.plot(E, [by_w[w][s] for w in WAVELENGTHS], "o-", color=COLORS[key], alpha=0.5, lw=1, ms=4)
        ax.set_title(geo.get(key).label, fontsize=10)
        ax.set_xlabel("photon energy  [eV]")
        ax.grid(alpha=0.25)
    np.atleast_1d(axes)[0].set_ylabel("threshold  [I_dep]")
    fig.tight_layout()
    fig.savefig(f"{OUT}/threshold_vs_energy.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
