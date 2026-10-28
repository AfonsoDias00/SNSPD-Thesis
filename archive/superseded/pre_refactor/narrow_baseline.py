#!/usr/bin/env python3
"""
narrow_baseline.py -- plain *narrow* wire control for photon_hitmap_v2.py
==========================================================================

Why this exists
---------------
The zipper's channels are 56 nm wide.  With the energy-conserving hotspot the zipper's photon thresholds drop
by 0.01-0.04 relative to the legacy hotspot (the plain 200 nm wire barely moves, 0.002-0.007).  That is what
hotspot-energy confinement by narrow walls looks like -- and a plain 56 nm wire should show it too, without the
12 nm bridges that create the dead zone.  Whether the zipper beats a narrow plain wire on *normalised* bias is
therefore the test that decides whether this geometry is worth iterating on.

It re-uses every numerical ingredient of photon_hitmap_v2.py (same hotspot model, same detection criterion,
same bisection, same mesh density) and only changes the film width.  Keep photon_hitmap_v2.py in the same folder.

Usage
-----
  python narrow_baseline.py --hotspot conserving --tmax 150 --compare hitmap_conserving_T150_m1.4.csv
  python narrow_baseline.py --hotspot legacy     --tmax 150 --compare hitmap_legacy_T150_m1.4.csv
  python narrow_baseline.py --report-only --compare ...        # redo the table / figure from saved files

--hotspot, --tmax and --mesh-xi must match the run you compare against (a warning is printed if they do not).
--compare may be a path or a basename from results_photon_hitmap_v2/.
The run is resumable (rows are appended to a CSV as they finish) and serial (about 8 + 5 x 8 + 6 x 8 solver
runs; a 56 nm film is small, so each run is much faster than a 200 nm one).

Settings that scale with the width (so the numerical resolution matches the 200 nm runs)
  * bisection tolerance: 4 uA * width / 200 nm  (override with --tol)
  * collapse threshold of the monitored supercurrent slice: 50 uA * width / 200 nm
  * bracket for the dark switching current: 4.2 ... 7.0 uA per nm of width
Lower photon-threshold bracket: --lo-frac x I_sw (default 0.30; photon_hitmap_v2 uses 0.70, which the narrow
wire's f = 1 thresholds fall below).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

import photon_hitmap_v2 as p

GEOM = "Narrow_Wire"
W_ORIG = p.W   # the 200 nm film width, needed for the zipper's free-width figure after we patch p.W
LABEL = {"Plain_Wire": "Plain 200 nm", "Long_Zipper": "Long zipper", GEOM: "Plain narrow"}


def setup(width_nm: float, lo_frac: float):
    """Patch photon_hitmap_v2's film width and register the narrow wire as a geometry."""
    wn = width_nm * 1e-3
    half = wn / 2.0
    p.W = wn
    p.COLLAPSE_UA = 50.0 * width_nm / 200.0
    p.HIT_LO_FRAC = lo_frac
    p.GEOMETRIES[GEOM] = dict(holes=False, slices=[(0.002, half - 0.002, 0.0)],
                              isw_bracket=(4.2 * width_nm, 7.0 * width_nm))
    p.PROBES[GEOM] = [(int(round(x)), 0) for x in probe_positions(width_nm)]


def hit_positions(width_nm: float):
    """Sites 2, 8, 14, 21 nm in from the edge (those that fit the half-width), then the quarter points and the
    centre; a site within 3 nm of one already chosen is skipped.  56 nm gives the original five sites."""
    h = width_nm / 2.0
    xs = [-(h - o) for o in (2.0, 8.0, 14.0, 21.0) if o <= h - 3.0]
    for x in (-float(round(h / 2.0)), -float(round(h / 4.0)), 0.0):
        if all(abs(x - s) >= 3.0 for s in xs):
            xs.append(x)
    return sorted(xs)


def probe_positions(width_nm: float):
    """Edge site, the site nearest the half-way point to the edge, and the centre."""
    xs = hit_positions(width_nm)
    return [xs[0], min(xs, key=lambda x: abs(x + width_nm / 4.0)), 0.0]


def j_dep():
    """GL depairing sheet current d * Phi0 / (3 sqrt(3) pi mu0 lambda^2 xi) of the uniform film [uA/nm]."""
    phi0, mu0 = 2.067833848e-15, 4e-7 * np.pi
    lam, xi_, d_ = p.london_lambda * 1e-6, p.xi * 1e-6, p.d * 1e-6
    return d_ * phi0 / (3.0 * np.sqrt(3.0) * np.pi * mu0 * lam ** 2 * xi_) * 1e-3   # A/m -> uA/nm


def row_exists(have, x, y, f):
    return p.row_key(GEOM, x, y, f) in have


def narrow_csv_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def run(cfg, width_nm, csv_path, isw_path):
    os.makedirs(cfg.out_dir, exist_ok=True)
    p._need_tdgl()
    isw = p.load_json(isw_path, {})
    if GEOM not in isw:
        p.log(f"-- dark switching current of the {width_nm:g} nm wire")
        r = p.task_fn(dict(kind="isw", geom=GEOM), cfg)
        if not r["ok"]:
            sys.exit(f"I_sw failed: {r['error']}")
        isw[GEOM] = r["res"]
        p.save_json(isw_path, isw)
        p.log(f"   I_sw = {r['res']['I_sw']:.1f} uA  ({r['res']['n_sims']} runs, {r['wall_s']:.0f} s)"
              + (f"  [{r['res']['flag']}]" if r["res"]["flag"] else ""))
    i_sw = isw[GEOM]["I_sw"]
    have = {p.row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]) for r in narrow_csv_rows(csv_path)}
    probes = {(float(x), float(y)) for x, y in p.PROBES[GEOM]}
    tasks = []
    for x in hit_positions(width_nm):
        x = float(round(x, 3))
        tasks.append((x, 1.0))
        if (x, 0.0) in probes or (round(x), 0.0) in probes:
            tasks += [(x, f) for f in cfg.energies]
    tasks = [(x, f) for x, f in tasks if not row_exists(have, x, 0.0, f)]
    p.log(f"-- {len(tasks)} threshold searches to run")
    new_file = not os.path.exists(csv_path)
    for k, (x, f) in enumerate(tasks, 1):
        r = p.task_fn(dict(kind="hit", geom=GEOM, x_nm=x, y_nm=0.0, f=f, I_sw=i_sw), cfg)
        if not r["ok"]:
            p.log(f"   FAILED x={x} f={f}: {r['error']}")
            continue
        res = r["res"]
        row = dict(Geometry=GEOM, x_nm=x, y_nm=0.0, energy_f=f, tmax=cfg.tmax, hotspot=cfg.hotspot,
                   mesh_xi=cfg.mesh_xi, I_sw_uA=i_sw, I_det_uA=res["I_det"], i_det=res["I_det"] / i_sw,
                   flag=res["flag"], n_sims=res["n_sims"], wall_s=round(r["wall_s"], 1))
        with open(csv_path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=p.CSV_FIELDS)
            if new_file:
                w.writeheader()
                new_file = False
            w.writerow(row)
        p.log(f"   [{k}/{len(tasks)}] x = {x:6.1f} nm  f = {f:.2f}  i_det = {row['i_det']:.3f}"
              f"{'  [' + res['flag'] + ']' if res['flag'] else ''}   ({res['n_sims']} runs, {r['wall_s']:.0f} s)")


def narrow_vals(samples, half_nm):
    xs = np.array([abs(s[0]) for s in samples])
    vs = np.array([s[2] for s in samples])
    o = np.argsort(xs)
    grid = np.linspace(-half_nm, half_nm, 1121)
    return np.sort(np.interp(np.abs(grid), xs[o], vs[o]))


def compare_rows(path):
    candidates = [path]
    if not os.path.isabs(path):
        project_dir = os.path.dirname(os.path.abspath(__file__))
        candidates.extend([
            os.path.join(project_dir, path),
            os.path.join(project_dir, "results_photon_hitmap_v2", path),
        ])
    compare_path = next((candidate for candidate in candidates if os.path.isfile(candidate)), None)
    if compare_path is None:
        searched = ", ".join(candidates)
        raise SystemExit(
            f"comparison CSV not found: {path!r}. "
            f"Generate the matching photon_hitmap_v2.py run or pass its path "
            f"(searched: {searched})"
        )
    with open(compare_path, newline="") as fh:
        return list(csv.DictReader(fh))


def isw_from_rows(rows, geom):
    for r in rows:
        if r["Geometry"] == geom:
            return float(r["I_sw_uA"])
    return None


def fig_eff(cfg, vals, met, isw, path):
    import matplotlib.pyplot as plt
    p.apply_style()
    colors = {"Plain_Wire": p.C_PLAIN, "Long_Zipper": p.C_ZIP, GEOM: p.INK2}
    styles = {"Plain_Wire": "-", "Long_Zipper": "-", GEOM: "--"}
    lo = min(0.75, np.floor(min(v[0] for v in vals.values()) * 20) / 20 - 0.05)
    ib = np.linspace(lo, 1.0, 800)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), layout="constrained")
    for g, v in vals.items():
        ax.plot(ib, p.eta_curve(v, ib), color=colors[g], ls=styles[g],
                label=f"{LABEL[g]} (I_sw = {isw[g]:.0f} uA):  i_50 = {met[g]['i_50']:.3f},  i_95 = {met[g]['i_95']:.3f}")
        ax.axvline(met[g]["i_95"], color=colors[g], ls=":", lw=1.1)
    ax.set_xlabel("normalised bias  i_b = I_bias / I_sw   (dotted: i_95)")
    ax.set_ylabel("area-integrated efficiency [%]")
    ax.set_ylim(-3, 103)
    ax.legend(loc="lower right", fontsize=8.5)
    ax.set_title(f"Zipper vs plain wires, normalised bias  [{cfg.hotspot}, T = {cfg.tmax:g} tau_0]")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def report(cfg, width_nm, csv_path, isw_path, compare):
    isw_j = p.load_json(isw_path, {})
    nrows = narrow_csv_rows(csv_path)
    samples = p.photon_rows(nrows, GEOM)
    if not samples or GEOM not in isw_j:
        p.log("no narrow-wire photon rows yet.")
        return
    rows = list(nrows)
    geoms = [GEOM]
    if compare:
        crow = compare_rows(compare)
        bad = {(r["hotspot"], float(r["tmax"]), float(r["mesh_xi"])) for r in crow} - {(cfg.hotspot, cfg.tmax, cfg.mesh_xi)}
        if bad:
            p.log(f"WARNING: --compare file uses {sorted(bad)}, this run uses {(cfg.hotspot, cfg.tmax, cfg.mesh_xi)}")
        rows += crow
        geoms = ["Plain_Wire", "Long_Zipper", GEOM]
    isw = {GEOM: isw_j[GEOM]["I_sw"]}
    if compare:
        for g in ("Plain_Wire", "Long_Zipper"):
            isw[g] = isw_from_rows(crow, g)
    vals, met, tabs, sel = {}, {}, {}, {}
    for g in geoms:
        s = p.photon_rows(rows, g)
        if not s:
            continue
        if g == "Plain_Wire":
            vals[g] = p.pixel_values(p.plain_map(s)[2])
        elif g == "Long_Zipper":
            vals[g] = p.pixel_values(p.zipper_map(s)[2])
        else:
            vals[g] = narrow_vals(s, width_nm / 2.0)
        met[g] = p.metrics(vals[g])
        tabs[g] = p.probe_table(rows, g)
        sel[g] = p.selectivity(tabs[g])
    geoms = [g for g in geoms if g in vals]
    fig_eff(cfg, vals, met, isw, cfg.path(f"narrow{width_nm:g}nm_efficiency", "png"))

    L = []
    P = L.append
    P("=" * 78)
    P(f"NARROW-WIRE BASELINE   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau_0   mesh={cfg.mesh_xi:g} xi   "
      f"narrow width={width_nm:g} nm")
    P("=" * 78)
    P(f"{'':<30}" + "".join(f"{LABEL[g]:>15}" for g in geoms))
    P(f"{'I_sw [uA]':<30}" + "".join(f"{isw[g]:>15.1f}" for g in geoms))
    wid = {"Plain_Wire": 200.0, "Long_Zipper": None, GEOM: width_nm}
    if "Long_Zipper" in geoms:
        w_patched, p.W = p.W, W_ORIG
        wid["Long_Zipper"] = p.free_width_min() * 1e3
        p.W = w_patched
    P(f"{'I_sw / free width [uA/nm]':<30}" + "".join(f"{isw[g] / wid[g]:>15.2f}" for g in geoms))
    P(f"{'I_sw / I_dep (narrowest)':<30}" + "".join(f"{isw[g] / (j_dep() * wid[g]):>15.3f}" for g in geoms))
    for key, name in [("i_10", "i_10"), ("i_50", "i_50"), ("i_90", "i_90"), ("i_95", "i_95"), ("i_99", "i_99"),
                      ("i_max", "i_max (100 % point)")]:
        P(f"{name:<30}" + "".join(f"{met[g][key]:>15.3f}" for g in geoms))
    for key in ("i_50", "i_95"):
        P(f"{key + '  in units of I_dep':<30}"
          + "".join(f"{met[g][key] * isw[g] / (j_dep() * wid[g]):>15.3f}" for g in geoms))
    P(f"{'window below the wall: 1 - i_95':<30}" + "".join(f"{1 - met[g]['i_95']:>15.3f}" for g in geoms))
    P(f"{'  (photon-level sites used)':<30}" + "".join(f"{len(p.photon_rows(rows, g)):>15d}" for g in geoms))
    if any(sel.values()):
        P("-" * 78)
        P("selectivity  median[i_det(f) - i_det(1)] over probe sites")
        for f in sorted({f for s in sel.values() for f in s}):
            P(f"  f = {f:<5g}".ljust(30) + "".join(f"{p.fmt(sel[g].get(f)):>15}" for g in geoms))
    if "Long_Zipper" in geoms:
        z, n = met["Long_Zipper"], met[GEOM]
        P("-" * 78)
        P("zipper minus narrow plain wire  (positive = the narrow wire turns on / saturates at the lower bias)")
        P(f"  onset    i_50: zipper - narrow = {z['i_50'] - n['i_50']:+.3f}      i_10: {z['i_10'] - n['i_10']:+.3f}")
        P(f"  plateau  i_95: zipper - narrow = {z['i_95'] - n['i_95']:+.3f}      i_99: {z['i_99'] - n['i_99']:+.3f}")
        dominated = z["i_50"] >= n["i_50"] and z["i_95"] >= n["i_95"]
        if dominated:
            P("  RESULT: the narrow plain wire is at least as good as the zipper on onset AND plateau in normalised bias.")
            P("          -> on the stated goal (lower normalised bias, wider window below I_sw) the zipper is dominated;")
            P("             what is left to argue for it is total signal current, see note.")
        else:
            P("  RESULT: the narrow wire does NOT dominate the zipper on both measures -- inspect the table above.")
        P(f"  note: the zipper carries I_sw = {isw['Long_Zipper']:.0f} uA against {isw[GEOM]:.0f} uA for one narrow wire; this")
        P("        baseline only decides the normalised-bias claim, not a signal-amplitude (SNAP-style) argument.")
    P("=" * 78)
    txt = "\n".join(L)
    print(txt)
    with open(cfg.path(f"narrow{width_nm:g}nm_summary", "txt"), "w") as fh:
        fh.write(txt + "\n")


def scan_report(cfg, widths, lo_frac, compare):
    """One table over every plain-wire width already simulated for this tag, plus the zipper from --compare."""
    cols = []                                                # (label, I_sw, I_dep, metrics, sites, selectivity)
    for w in widths:
        setup(w, lo_frac)
        isw_j = p.load_json(cfg.path(f"narrow{w:g}nm_isw", "json"), {})
        rows = narrow_csv_rows(cfg.path(f"narrow{w:g}nm_hitmap", "csv"))
        s = p.photon_rows(rows, GEOM)
        if not s or GEOM not in isw_j:
            p.log(f"  {w:g} nm: no results for {cfg.tag}, skipped")
            continue
        flags = sorted({r["flag"] for r in rows if r["flag"]})
        if flags:
            p.log(f"  {w:g} nm: flagged rows present {flags} -- these thresholds are bounds, not values")
        cols.append((f"{w:g} nm", isw_j[GEOM]["I_sw"], j_dep() * w, p.metrics(narrow_vals(s, w / 2.0)), len(s),
                     p.selectivity(p.probe_table(rows, GEOM))))
    if compare:
        crow = compare_rows(compare)
        w_patched, p.W = p.W, W_ORIG
        free = p.free_width_min() * 1e3
        p.W = w_patched
        for g, label, wid, mapper in (("Plain_Wire", "200 nm", 200.0, p.plain_map), ("Long_Zipper", "zipper", free, p.zipper_map)):
            s = p.photon_rows(crow, g)
            if s:
                cols.append((label, isw_from_rows(crow, g), j_dep() * wid, p.metrics(p.pixel_values(mapper(s)[2])),
                             len(s), p.selectivity(p.probe_table(crow, g))))
    L = []
    P = L.append
    P("=" * 96)
    P(f"WIDTH SCAN   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau_0   mesh={cfg.mesh_xi:g} xi   "
      f"I_dep = {j_dep():.3f} uA/nm x narrowest width (GL)")
    P("=" * 96)
    P(f"{'':<28}" + "".join(f"{c[0]:>11}" for c in cols))
    P(f"{'I_sw [uA]':<28}" + "".join(f"{c[1]:>11.1f}" for c in cols))
    P(f"{'I_sw / I_dep':<28}" + "".join(f"{c[1] / c[2]:>11.3f}" for c in cols))
    for key in ("i_10", "i_50", "i_95", "i_max"):
        P(f"{key + '  [I_sw]':<28}" + "".join(f"{c[3][key]:>11.3f}" for c in cols))
    for key in ("i_50", "i_95"):
        P(f"{key + '  [I_dep]':<28}" + "".join(f"{c[3][key] * c[1] / c[2]:>11.3f}" for c in cols))
    P(f"{'window 1 - i_95':<28}" + "".join(f"{1 - c[3]['i_95']:>11.3f}" for c in cols))
    P(f"{'  (photon-level sites)':<28}" + "".join(f"{c[4]:>11d}" for c in cols))
    fs = sorted({f for c in cols for f in c[5]})
    if fs:
        P("selectivity  median[i_det(f) - i_det(1)] over probe sites")
        for f in fs:
            P(f"  f = {f:<5g}".ljust(28) + "".join(f"{p.fmt(c[5].get(f)):>11}" for c in cols))
    P("=" * 96)
    txt = "\n".join(L)
    print(txt)
    with open(cfg.path("width_scan_summary", "txt"), "w") as fh:
        fh.write(txt + "\n")


def parse():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--width-nm", type=float, default=56.0)
    ap.add_argument("--hotspot", choices=["conserving", "legacy"], default="conserving")
    ap.add_argument("--tmax", type=float, default=150.0)
    ap.add_argument("--mesh-xi", type=float, default=1.4)
    ap.add_argument("--tol", type=float, default=None, help="bisection tolerance [uA]; default 4 uA * width / 200 nm")
    ap.add_argument("--energies", default="0.6,0.35")
    ap.add_argument("--lo-frac", type=float, default=0.30, help="lower photon-threshold bracket [fraction of I_sw]")
    ap.add_argument("--compare", default=None, help="hitmap_<tag>.csv written by photon_hitmap_v2.py")
    ap.add_argument("--out", default="results_narrow_baseline")
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--scan", default=None, help="comma-separated widths [nm]: print the width-scan table only")
    return ap.parse_args()


def main():
    a = parse()
    setup(a.width_nm, a.lo_frac)
    cfg = p.Cfg(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi,
                tol=a.tol if a.tol is not None else 4.0 * a.width_nm / 200.0,
                energies=tuple(float(s) for s in a.energies.split(",") if s), out_dir=a.out, workers=1)
    if a.scan:
        scan_report(cfg, [float(s) for s in a.scan.split(",") if s], a.lo_frac, a.compare)
        return
    csv_path = cfg.path(f"narrow{a.width_nm:g}nm_hitmap", "csv")
    isw_path = cfg.path(f"narrow{a.width_nm:g}nm_isw", "json")
    if not a.report_only:
        run(cfg, a.width_nm, csv_path, isw_path)
    report(cfg, a.width_nm, csv_path, isw_path, a.compare)


if __name__ == "__main__":
    main()