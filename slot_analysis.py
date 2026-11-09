#!/usr/bin/env python3
"""
slot_analysis.py -- the slot layouts against their plain references, from the stored maps (no simulation).  7 Nov 2026.

Everything in units of I_dep of each device's narrowest free cross-section (120 nm for the slot layouts, the strip
width for plain strips), preset nbtin (gamma 10), kinetic hotspot f = 1, 400 tau0, voltage trigger.

  python slot_analysis.py        -> results/slot_analysis.txt, results/slot_efficiency.png, results/slot_maps.png

Sections: dark I_sw; efficiency quantiles and eta at fixed biases; the cascade picture (channel thresholds against
1 - w_channel / w_free and the plain strip of the channel width); gains over the plain references; where the
insensitive area is; latency; selectivity; slot length (150 vs 400 nm, critical sites); mesh (2.0 vs 1.4 xi, critical
sites); flux scans if present.
"""
from __future__ import annotations

import os

import numpy as np

from snspd import analysis as an
from snspd import geometry as geo
from snspd import plotting as pl
from snspd import sim, store

FILM = sim.PRESETS["nbtin"].film
TAG = "kinetic-nbtin_T400_m{m}_trigger"
MAPS, MAPS150, FLUX = "results/maps", "results/maps150", "results/flux"
REFS = ("Plain_Wire", "Plain_120", "Plain_60", "Plain_40")
LAYOUTS = ("Single_400", "Double_400", "Stagger_400")
CHANNEL = {"Single_400": 60.0, "Double_400": 40.0, "Stagger_400": 40.0}
GO_MARGIN = 0.03
HEAT_CSV = "results/heat/check_dt0.025.csv"   # heated thresholds with dt_max 0.025 (check.csv: dt_max 0.2, superseded)
T0 = 20.0


def i_dep(key):
    return FILM.I_dep_uA(geo.get(key).free_width_min_nm())


def load(d, m):
    tag = TAG.format(m=m)
    rows = store.load_rows(f"{d}/hitmap_{tag}.csv")
    isw = {k: v["I_sw"] for k, v in store.load_json(f"{d}/isw_{tag}.json", {}).items()}
    return rows, isw


def site_rows(rows, key, f=1.0):
    """{(x, y): row} at energy factor f."""
    return {(float(r["x_nm"]), float(r["y_nm"])): r for r in rows
            if r["Geometry"] == key and abs(float(r["energy_f"]) - f) < 1e-6}


def i_of(r, key):
    """Threshold in I_dep units; no detection below I_sw ('at_wall') counts as I_sw."""
    return float(r["I_sw_uA" if r["flag"] == "at_wall" else "I_det_uA"]) / i_dep(key)


def censored(r):
    """A threshold run that fired at t = 55.0: the settled start failed (near dark switching) and the bias-on
    transient of the uniform start fired the trigger.  The site does not detect below ~ that bias."""
    return r["t_trigger"] not in ("", None) and abs(float(r["t_trigger"]) - 55.0) < 0.01


def eta_map(key, rows, isw):
    g = geo.get(key)
    samples = [(x, y, i_of(r, key)) for (x, y), r in site_rows(rows, key).items()]
    X, Y, Z = an.unit_map(g, samples)
    return g, X, Y, Z, an.pixel_values(Z)


def region_of(g, x, y):
    """channel / slot end / bridge for a slot layout site (rough, for the insensitive-area breakdown)."""
    b = g.cell_half_y_nm - geo.BRIDGE_NM / 2.0
    if abs(y) > b:
        return "bridge row"
    if abs(y) > b - 40.0:
        return "slot end"
    return "channel"


def main():
    L = []
    P = L.append
    rows, isw = load(MAPS, 2)
    keys = [k for k in REFS + LAYOUTS if k in isw and site_rows(rows, k)]
    P("SLOT LAYOUTS vs PLAIN REFERENCES  (preset nbtin: gamma 10, NbTiN hotspot times; kinetic f = 1; 400 tau0; "
      "voltage trigger; dark I_sw sustained; 2.0 xi, shortened strips)")
    P("units: I_dep of the narrowest free cross-section (slot layouts: 120 nm, I_dep = %.1f uA; plain strips: their "
      "width)" % i_dep("Double_400"))
    P("")
    # ---------------------------------------------------------------- quantiles
    maps, met = {}, {}
    for k in keys:
        g, X, Y, Z, vals = eta_map(k, rows, isw)
        maps[k] = (g, X, Y, Z)
        met[k] = an.metrics(vals)
        met[k]["vals"] = vals
    P(f"{'device':<13}{'free nm':>8}{'sites':>6}{'I_sw':>7}{'i_10':>7}{'i_50':>7}{'i_90':>7}{'i_95':>7}{'i_max':>7}"
      f"{'eta(0.70)':>10}{'eta(0.80)':>10}{'eta(0.90)':>10}{'eta(0.95)':>10}")
    for k in keys:
        m, v = met[k], met[k]["vals"]
        e = lambda b: an.eta_curve(v, b)
        P(f"{k:<13}{geo.get(k).free_width_min_nm():8.0f}{len(site_rows(rows, k)):6d}{isw[k] / i_dep(k):7.3f}"
          f"{m['i_10']:7.3f}{m['i_50']:7.3f}{m['i_90']:7.3f}{m['i_95']:7.3f}{m['i_max']:7.3f}"
          f"{e(0.70):9.0f}%{e(0.80):9.0f}%{e(0.90):9.0f}%{e(0.95):9.0f}%")
    # ---------------------------------------------------------------- cascade
    P("")
    P("CASCADE: a photon that blocks one channel sends its current into the others, which switch if the bias exceeds "
      "1 - w_channel / w_free.  The hit channel must also detect on its own: compare the plain strip of the channel width.")
    P(f"{'device':<13}{'channels':>14}{'cascade limit':>15}{'plain strip of w_ch':>21}{'channel sites (median, min-max)':>34}")
    for k in LAYOUTS:
        if k not in met:
            continue
        g, wc = geo.get(k), CHANNEL[k]
        lim = 1.0 - wc / g.free_width_min_nm()
        ref = f"Plain_{wc:g}"
        refv = np.median([i_of(r, ref) for r in site_rows(rows, ref).values()]) if ref in keys else float("nan")
        ch = [i_of(r, k) for (x, y), r in site_rows(rows, k).items() if region_of(g, x, y) == "channel"]
        nch = "2 x 60 nm" if k.startswith("Single") else "3 x 40 nm"
        P(f"{k:<13}{nch:>14}{lim:15.3f}{refv:21.3f}{np.median(ch):18.3f} ({min(ch):.3f}-{max(ch):.3f})")
    P("  -> the channel threshold is about max(cascade limit, local threshold of a strip of the channel width).")
    # ---------------------------------------------------------------- gains
    P("")
    P(f"GAINS (reference - device, positive = the slot device needs less bias; go if > {GO_MARGIN} at the plateau i_95)")
    for ref in ("Plain_Wire", "Plain_120"):
        if ref not in met:
            continue
        for k in LAYOUTS:
            if k not in met:
                continue
            d = {q: met[ref][q] - met[k][q] for q in ("i_10", "i_50", "i_95", "i_max")}
            go = "GO" if d["i_95"] > GO_MARGIN else "no-go"
            P(f"  {k:<12} vs {ref:<11}" + "".join(f"  {q} {v:+.3f}" for q, v in d.items()) + f"   -> {go} (plateau)")
    # ---------------------------------------------------------------- where the insensitive area is
    P("")
    P("WHERE IS THE INSENSITIVE AREA (f = 1 sites by region; '*' = site censored at the dark-switching edge, t = 55)")
    for k in LAYOUTS:
        if k not in met:
            continue
        g = geo.get(k)
        by = {}
        for (x, y), r in site_rows(rows, k).items():
            by.setdefault(region_of(g, x, y), []).append((i_of(r, k), censored(r)))
        P(f"  {k}: " + ";  ".join(f"{reg} {len(v)} sites median {np.median([a for a, _ in v]):.3f} "
                                  f"max {max(a for a, _ in v):.3f}" + (f" ({sum(c for _, c in v)}*)" if any(c for _, c in v) else "")
                                  for reg, v in sorted(by.items())))
        Z = maps[k][3]
        film_px = np.isfinite(Z)
        P(f"      area fraction above 0.85: {np.mean(Z[film_px] > 0.85) * 100:.0f} %, above 0.90: "
          f"{np.mean(Z[film_px] > 0.90) * 100:.0f} %   (bridge rows are {2 * geo.BRIDGE_NM / 2 / (2 * g.cell_half_y_nm) * 100:.0f} % of the cell length)")
    # ---------------------------------------------------------------- latency
    P("")
    P("LATENCY of the switching run at threshold, after the photon [tau0]: median / 90 % / max")
    for k in keys:
        ts = [float(r["t_trigger"]) - T0 for r in site_rows(rows, k).values() if r["t_trigger"] and not censored(r)]
        if ts:
            P(f"  {k:<13}{np.median(ts):6.0f} /{np.quantile(ts, 0.9):5.0f} /{max(ts):5.0f}")
    # ---------------------------------------------------------------- selectivity
    P("")
    P("SELECTIVITY: i_det(f) - i_det(1) at the probe sites (weak photons need more bias), median")
    for k in keys:
        tab = an.probe_table(rows, geo.get(k))
        sel = an.selectivity(tab)
        if sel:
            P(f"  {k:<13}" + "  ".join(f"f = {f:g}: {v:+.3f}" for f, v in sorted(sel.items(), reverse=True)))
    # ---------------------------------------------------------------- slot length
    r150, isw150 = load(MAPS150, 2)
    if r150:
        P("")
        P("SLOT LENGTH: critical sites, 150 nm vs 400 nm slots (2.0 xi)")
        for lay in ("Single", "Double", "Stagger"):
            a, b = f"{lay}_150", f"{lay}_400"
            sa, sb = site_rows(r150, a), site_rows(rows, b)
            if not sa or not sb:
                continue
            ga, gb = geo.get(a), geo.get(b)
            ba, bb = ga.cell_half_y_nm - 20.0, gb.cell_half_y_nm - 20.0
            # match sites by region: same x, y scaled from the slot centre / slot end / bridge
            def label(g, b_, x, y):
                if y == 0:
                    return f"({x:g}, 0)"
                return f"({x:g}, end{y - b_:+g})"
            la = {label(ga, ba, x, y): i_of(r, a) for (x, y), r in sa.items()}
            lb = {label(gb, bb, x, y): i_of(r, b) for (x, y), r in sb.items()}
            common = [s for s in la if s in lb]
            P(f"  {lay}: dark I_sw {isw150.get(a, np.nan) / i_dep(a):.3f} vs {isw[b] / i_dep(b):.3f};  "
              + ",  ".join(f"{s} {la[s]:.3f} vs {lb[s]:.3f}" for s in sorted(common)))
    # ---------------------------------------------------------------- mesh
    r14, isw14 = load(MAPS, 1.4)
    if r14:
        P("")
        P("MESH: critical sites at 1.4 xi vs the 2.0 xi maps")
        for k in ("Plain_Wire",) + LAYOUTS:
            s14, s20 = site_rows(r14, k), site_rows(rows, k)
            common = [s for s in s14 if s in s20]
            if not common:
                continue
            d = [i_of(s14[s], k) - i_of(s20[s], k) for s in common]
            P(f"  {k:<12} dark I_sw {isw14.get(k, np.nan) / i_dep(k):.3f} vs {isw[k] / i_dep(k):.3f};  "
              f"{len(common)} sites, 1.4 - 2.0: mean {np.mean(d):+.3f}, range {min(d):+.3f}..{max(d):+.3f};  "
              + ", ".join(f"({x:g},{y:g}) {i_of(s14[(x, y)], k):.3f}" for x, y in sorted(common, key=lambda s: (s[1], s[0]))))
    # ---------------------------------------------------------------- flux
    for k in LAYOUTS:
        p = f"{FLUX}/scan_{k}_kinetic-nbtin_T400_m2.csv"
        sc = store.load_rows(p)
        if not sc:
            continue
        P("")
        P(f"FLUX SCAN {k}: outcome counts per bias (critical sites)")
        biases = sorted({float(r["bias"]) for r in sc})
        outs = ("clean", "trapped", "transient", "direct", "via_trap")
        P(f"  {'bias':>6}" + "".join(f"{o:>11}" for o in outs))
        for b in biases:
            rr = [r for r in sc if float(r["bias"]) == b]
            P(f"  {b:6.2f}" + "".join(f"{sum(r['outcome'] == o for r in rr):11d}" for o in outs))
    # ---------------------------------------------------------------- heat
    heat = [r for r in store.load_rows(HEAT_CSV) if r["variant"] == "NbTiN 2T"]
    hcurves = {}
    if heat:
        P("")
        P("WITH HEAT (NbTiN two-temperature model, heat_check.py): critical sites only, the same sites for every "
          "criterion,")
        P("interpolated over the unit cell.  'sustained' from the flux scans (0.03 steps); plain wire: the bounds coincide.")
        P(f"{'device':<13}{'criterion':<11}{'sites':>6}{'i_10':>7}{'i_50':>7}{'i_90':>7}{'i_95':>7}{'i_max':>7}")
        for k in ("Plain_Wire",) + LAYOUTS:
            g = geo.get(k)
            crit = {(float(x), float(y)) for x, y in g.grid_critical}
            hv = {(float(r["x_nm"]), float(r["y_nm"])): float(r["i_heat"]) for r in heat if r["geom"] == k}
            tv = {s: i_of(r, k) for s, r in site_rows(rows, k).items() if s in crit}
            sv = {}
            for r in store.load_rows(f"{FLUX}/scan_{k}_kinetic-nbtin_T400_m2.csv"):
                if r["outcome"] in ("direct", "via_trap"):
                    s = (float(r["x_nm"]), float(r["y_nm"]))
                    sv[s] = min(sv.get(s, 9.0), float(r["bias"]))
            if k == "Plain_Wire":
                sv = dict(tv)
            for name, vals in (("trigger", tv), ("sustained", sv), ("heat", hv)):
                common = [s for s in crit if s in vals]
                if len(common) < len(crit):
                    continue
                _, _, Z = an.unit_map(g, [(x, y, vals[(x, y)]) for x, y in common])
                v = an.pixel_values(Z)
                m = an.metrics(v)
                hcurves[(k, name)] = v
                P(f"{k:<13}{name:<11}{len(common):6d}{m['i_10']:7.3f}{m['i_50']:7.3f}{m['i_90']:7.3f}{m['i_95']:7.3f}"
                  f"{m['i_max']:7.3f}")
        # full heated maps where every full-grid site has a heated threshold
        for k in ("Plain_Wire",) + LAYOUTS:
            g = geo.get(k)
            hv = {(float(r["x_nm"]), float(r["y_nm"])): float(r["i_heat"]) for r in heat if r["geom"] == k}
            full = [(float(x), float(y)) for x, y in g.hit_grid("full")]
            if all(s in hv for s in full):
                _, _, Z = an.unit_map(g, [(x, y, hv[(x, y)]) for x, y in full])
                v = an.pixel_values(Z)
                m = an.metrics(v)
                hcurves[(k, "heat, full map")] = v
                P(f"{k:<13}{'heat, FULL':<11}{len(full):6d}{m['i_10']:7.3f}{m['i_50']:7.3f}{m['i_90']:7.3f}"
                  f"{m['i_95']:7.3f}{m['i_max']:7.3f}   (trigger full map: i_50 {met[k]['i_50']:.3f}, "
                  f"i_95 {met[k]['i_95']:.3f})")
        P("  site by site (heat / trigger / sustained):")
        for k in ("Plain_Wire",) + LAYOUTS:
            hv = sorted(((float(r["x_nm"]), float(r["y_nm"])), float(r["i_heat"]), r["normal_nm2"]) for r in heat
                        if r["geom"] == k)
            sr = site_rows(rows, k)
            P(f"    {k:<12}" + ",  ".join(f"({x:g},{y:g}) {h:.3f}/{i_of(sr[(x, y)], k):.3f}" if (x, y) in sr
                                         else f"({x:g},{y:g}) {h:.3f}" for (x, y), h, _ in hv))
    # ---------------------------------------------------------------- optimisation candidates (results/opt)
    ro, iswo = load("results/opt", 2)
    if ro:
        heat_all = {(r["geom"], float(r["x_nm"]), float(r["y_nm"])): float(r["i_heat"])
                    for r in store.load_rows(HEAT_CSV) if r["variant"] == "NbTiN 2T"}
        P("")
        P("OPTIMISATION CANDIDATES (every feature >= 40 nm): critical sites, 2.0 xi; trigger (no heat) and with heat "
          "(NbTiN 2T),")
        P("interpolated over the unit cell; the reference rows use the same critical sites")
        P(f"{'device':<17}{'free nm':>8}{'I_sw':>7}{'criterion':>11}{'i_10':>7}{'i_50':>7}{'i_90':>7}{'i_95':>7}{'i_max':>7}"
          f"{'channel med':>13}{'bridge max':>12}")
        cand = ("Plain_Wire", "Plain_80", "Double_400", "Double_400_n20", "Double_600", "Double_800", "Double_1200",
                "Stagger_400", "Single_400", "Single_400_w120", "Single_800_w120")
        opt_curves = {}
        for k in cand:
            if k not in geo.GEOMETRIES:
                continue
            g = geo.get(k)
            src_rows, src_isw = (rows, isw) if k == "Double_400" else (ro, iswo)
            sr = site_rows(src_rows, k)
            crit = [(float(x), float(y)) for x, y in g.grid_critical]
            for name in ("trigger", "heat"):
                if name == "trigger":
                    vals = {s: i_of(sr[s], k) for s in crit if s in sr}
                else:          # every heated site of this layout: the critical grid plus any slot-end sites
                    vals = {s[1:]: v for s, v in heat_all.items() if s[0] == k}
                if not all(s_ in vals for s_ in crit):
                    continue
                _, _, Z = an.unit_map(g, [(x, y, v) for (x, y), v in vals.items()])
                if name == "heat":
                    opt_curves[k] = an.pixel_values(Z)
                m = an.metrics(an.pixel_values(Z))
                ch = [v for (x, y), v in vals.items() if not g.holes or region_of(g, x, y) == "channel"]
                br = [v for (x, y), v in vals.items() if g.holes and region_of(g, x, y) == "bridge row"]
                P(f"{k:<17}{g.free_width_min_nm():8.0f}{src_isw.get(k, np.nan) / i_dep(k):7.3f}{name:>11}"
                  f"{m['i_10']:7.3f}{m['i_50']:7.3f}{m['i_90']:7.3f}{m['i_95']:7.3f}{m['i_max']:7.3f}"
                  f"{np.median(ch):13.3f}{(max(br) if br else float('nan')):12.3f}")
        if opt_curves:
            P("  heat rows: every heated site of the layout (critical grid + slot-end sites where run), dt_max 0.025")
            fig_opt(opt_curves)
    text = "\n".join(L)
    print(text)
    with open("results/slot_analysis.txt", "w") as fh:
        fh.write(text + "\n")
    figures(keys, maps, met, isw)
    if hcurves:
        if "Double_400" in met:
            hcurves[("Double_400", "trigger, full map")] = met["Double_400"]["vals"]
        fig_heat(hcurves)


def fig_opt(curves):
    """Heated efficiency of every layout from its heated sites (dt_max 0.025)."""
    import matplotlib.pyplot as plt
    pl.apply_style()
    order = ["Plain_Wire", "Plain_80", "Double_400", "Double_600", "Double_800", "Double_1200", "Double_400_n20",
             "Stagger_400", "Single_400", "Single_400_w120", "Single_800_w120"]
    colors = {"Plain_Wire": "#2a78d6", "Plain_80": "#52514e", "Double_400": "#eb6834", "Double_600": "#f0a35e",
              "Double_800": "#c4511a", "Double_1200": "#8a3410", "Double_400_n20": "#e8b21c", "Stagger_400": "#d55181",
              "Single_400": "#1baf7a", "Single_400_w120": "#0f7a55", "Single_800_w120": "#7ad3b2"}
    fig, ax = plt.subplots(figsize=(8.8, 5.4), layout="constrained")
    ib = np.linspace(0.6, 1.0, 401)
    for k in order:
        if k in curves:
            dash = "--" if k.startswith("Plain") else "-"
            ax.plot(ib, an.eta_curve(curves[k], ib), dash, color=colors[k], lw=2.0, label=geo.get(k).label)
    ax.set_xlabel("bias  I / I_dep(narrowest free width)")
    ax.set_ylabel("hard-step efficiency η with heat  [% of the unit cell]")
    ax.set_xlim(0.6, 1.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=7.8, loc="upper left", frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID)
    ax.set_title("With heat (NbTiN 2T, dt_max 0.025, 2.0 ξ): critical sites + slot-end sites, interpolated", fontsize=10)
    fig.savefig("results/slot_efficiency_heat_all.png", dpi=170)
    plt.close(fig)
    print("wrote results/slot_efficiency_heat_all.png")


def fig_heat(hc):
    import matplotlib.pyplot as plt
    pl.apply_style()
    colors = {"Plain_Wire": "#2a78d6", "Single_400": "#1baf7a", "Double_400": "#eb6834", "Stagger_400": "#d55181"}
    styles = {"trigger": ":", "sustained": "--", "heat": "-", "heat, full map": "-"}
    fig, ax = plt.subplots(figsize=(8.4, 5.2), layout="constrained")
    ib = np.linspace(0.6, 1.0, 401)
    full = {k for k, name in hc if name == "heat, full map"}
    styles["trigger, full map"] = ":"
    for (k, name), v in hc.items():
        if k == "Plain_Wire" and name not in ("heat", "heat, full map"):
            continue
        if name == "heat" and k in full:          # the full heated map replaces the critical-site estimate
            continue
        if k in ("Single_400", "Stagger_400") and not name.startswith("heat"):
            continue                              # without heat: shown for the double column only
        if k == "Double_400" and name == "trigger":
            continue                              # the full trigger map instead
        heated = name.startswith("heat")
        what = ("with heat, full map" if name == "heat, full map" else "with heat, critical sites") if heated \
            else ("trigger, no heat, full map" if name == "trigger, full map" else name + ", no heat, critical sites")
        ax.plot(ib, an.eta_curve(v, ib), styles[name], color=colors[k], lw=2.4 if heated else 1.3,
                label=f"{geo.get(k).label}: {what}")
    ax.set_xlabel("bias  I / I_dep(narrowest free width)")
    ax.set_ylabel("hard-step efficiency η  [% of the unit cell]")
    ax.set_xlim(0.6, 1.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=7.8, loc="lower right", frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID, ncol=1)
    ax.set_title("Without heat (trigger / sustained) and with heat (NbTiN two-temperature model), γ = 10, 2.0 ξ", fontsize=10)
    fig.savefig("results/slot_efficiency_heat.png", dpi=170)
    plt.close(fig)
    print("wrote results/slot_efficiency_heat.png")


def figures(keys, maps, met, isw):
    import matplotlib.pyplot as plt
    pl.apply_style()
    colors = {"Plain_Wire": "#2a78d6", "Plain_120": "#52514e", "Plain_60": "#9a9893", "Plain_40": "#c3c1bb",
              "Single_400": "#1baf7a", "Double_400": "#eb6834", "Stagger_400": "#d55181"}
    fig, ax = plt.subplots(figsize=(8.4, 5.2), layout="constrained")
    ib = np.linspace(0.45, 1.0, 551)
    for k in keys:
        dash = "--" if k in ("Plain_60", "Plain_40", "Plain_120") else "-"
        ax.plot(ib, an.eta_curve(met[k]["vals"], ib), dash, color=colors[k], lw=2.0 if dash == "-" else 1.4,
                label=f"{geo.get(k).label}  (I_sw {isw[k] / i_dep(k):.2f})")
    for lim, txt in ((2 / 3, "3-channel cascade 0.67"), (0.5, "2-channel cascade 0.50")):
        ax.axvline(lim, color=pl.INK2, lw=0.8, ls=":")
        ax.text(lim + 0.004, 4, txt, rotation=90, fontsize=8, color=pl.INK2, va="bottom")
    ax.set_xlabel("bias  I / I_dep(narrowest free width)")
    ax.set_ylabel("hard-step efficiency η  [% of the unit cell]")
    ax.set_xlim(0.45, 1.0)
    ax.set_ylim(0, 102)
    ax.legend(fontsize=8.5, loc="lower right", frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID)
    ax.set_title("Slot layouts vs plain strips (γ = 10, NbTiN hotspot, f = 1, voltage trigger, 2.0 ξ)", fontsize=10.5)
    fig.savefig("results/slot_efficiency.png", dpi=170)
    plt.close(fig)

    lays = [k for k in LAYOUTS if k in maps]
    if not lays:
        return
    fig, axes = plt.subplots(1, len(lays), figsize=(2.6 * len(lays) + 1.4, 6.4), layout="constrained", sharey=True)
    axes = np.atleast_1d(axes)
    im = None
    for ax, k in zip(axes, lays):
        g, X, Y, Z = maps[k]
        im = ax.pcolormesh(X, Y, Z, cmap="viridis", vmin=0.6, vmax=1.0, shading="auto")
        pl.draw_holes(ax, g)
        ax.set_xlim(-g.width_nm / 2, g.width_nm / 2)
        ax.set_ylim(-g.cell_half_y_nm, g.cell_half_y_nm)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_title(f"{g.label}\ni_50 {met[k]['i_50']:.3f}, i_95 {met[k]['i_95']:.3f}", fontsize=9.5)
        ax.set_xlabel("x [nm]")
    axes[0].set_ylabel("y [nm]  (one unit cell)")
    cb = fig.colorbar(im, ax=axes, shrink=0.85)
    cb.set_label("photon threshold  I_det / I_dep(120 nm)")
    fig.savefig("results/slot_maps.png", dpi=170)
    plt.close(fig)
    print("wrote results/slot_efficiency.png, results/slot_maps.png")


if __name__ == "__main__":
    main()
