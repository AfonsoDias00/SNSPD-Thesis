#!/usr/bin/env python3
"""
narrow_baseline.py -- plain strips of any width as the control for the internal-hole geometries.

The long zipper's channels are 56 nm wide, and narrow strips are known to saturate at lower normalised bias
(Marsili et al., Nano Lett. 11, 2048 (2011); Lusche et al., J. Appl. Phys. 116, 043906 (2014)).  A hole geometry is
only interesting if it beats a plain strip as narrow as its channels.  This script runs plain strips with exactly the
model of hitmap.py (same hotspot, criterion, bisection and mesh; only the width changes), compares one width with the
hitmap.py results, and tabulates a width scan.

Usage
-----
  python narrow_baseline.py run    --width-nm 56 --mesh-xi 0.8 --compare hitmap_conserving_T150_m0.8.csv
  python narrow_baseline.py report --width-nm 56 --mesh-xi 0.8 --compare hitmap_conserving_T150_m0.8.csv
  python narrow_baseline.py scan   --widths 30,56,100 --mesh-xi 0.8 --compare hitmap_conserving_T150_m0.8.csv

--compare is a hitmap.py CSV (a path, or a name inside results/hitmap/); its hotspot, window and mesh must match.
Photon sites: 2, 8, 14, 21 nm in from the edge (those that fit), the quarter points and the centre.  Probe sites
(edge, half-way, centre) also get the weaker hotspots of --energies.
"""
from __future__ import annotations

import argparse
import os

import numpy as np

from snspd import analysis as an
from snspd import geometry as geo
from snspd import hotspot as hsm
from snspd import plotting as pl
from snspd import sim, store
from snspd.units import Film

OUT_DIR = "results/narrow_baseline"
HITMAP_DIR = "results/hitmap"


def paths(cfg, width_nm):
    return cfg.path(f"narrow{width_nm:g}nm_hitmap", "csv"), cfg.path(f"narrow{width_nm:g}nm_isw", "json")


def cmd_run(cfg, a):
    sim.need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    g = geo.plain_strip(a.width_nm)
    csv_path, isw_path = paths(cfg, a.width_nm)
    isw = store.load_json(isw_path, {})
    if g.key not in isw:
        sim.log(f"-- dark switching current of the {a.width_nm:g} nm strip")
        r = sim.task_fn(dict(kind="isw", geom=g, cfg=cfg, film=a.film))
        if not r["ok"]:
            raise SystemExit(f"I_sw failed: {r['error']}")
        isw[g.key] = r["res"]
        store.save_json(isw_path, isw)
        sim.log(f"   I_sw = {r['res']['I_sw']:.1f} uA  ({r['res']['n_sims']} runs, {r['wall_s']:.0f} s)"
                + (f"  [{r['res']['flag']}]" if r["res"]["flag"] else ""))
    i_sw = isw[g.key]["I_sw"]
    have = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]) for r in store.load_rows(csv_path)}
    tasks = []
    for x, y in g.hit_grid("full"):
        for f in (1.0,) + (a.energies if (x, y) in set(g.probes) else ()):
            if store.row_key(g.key, x, y, f) not in have:
                tasks.append(dict(kind="hit", geom=g, cfg=cfg, film=a.film, x_nm=x, y_nm=y, f=f, I_sw=i_sw))
    sim.log(f"-- {len(tasks)} threshold searches to run")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED x={t['x_nm']} f={t['f']}: {r['error']}")
            return
        row = store.hit_row(g.key, r["res"], t, cfg, r["wall_s"])
        store.append_row(csv_path, row)
        sim.log(f"   x = {t['x_nm']:6.1f} nm  f = {t['f']:.2f}  i_det = {row['i_det']:.3f}"
                + (f"  [{row['flag']}]" if row["flag"] else "") + f"   ({row['n_sims']} runs, {r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result)
    cmd_report(cfg, a)


def compare_rows(name):
    for cand in (name, os.path.join(HITMAP_DIR, name)):
        if os.path.isfile(cand):
            return store.load_rows(cand)
    raise SystemExit(f"comparison CSV not found: {name!r} (looked in . and {HITMAP_DIR}/)")


def isw_from_rows(rows, key):
    return next((float(r["I_sw_uA"]) for r in rows if r["Geometry"] == key), None)


def column(g, rows, i_sw, film):
    """Everything one table column needs for geometry g, from its rows."""
    s = an.photon_rows(rows, g.key)
    if not s or i_sw is None:
        return None
    vals = an.pixel_values(an.unit_map(g, s)[2])
    w_free = g.free_width_min_nm()
    return dict(g=g, I_sw=i_sw, I_dep=film.I_dep_uA(w_free), w_free=w_free, vals=vals, met=an.metrics(vals),
                n=len(s), sel=an.selectivity(an.probe_table(rows, g)),
                flags=sorted({r["flag"] for r in rows if r["Geometry"] == g.key and r["flag"]}))


def check_compare(cfg, crow):
    bad = {(r["hotspot"], float(r["tmax"]), float(r["mesh_xi"])) for r in crow} - {(cfg.hotspot, cfg.tmax, cfg.mesh_xi)}
    if bad:
        sim.log(f"WARNING: --compare file uses {sorted(bad)}, this run uses {(cfg.hotspot, cfg.tmax, cfg.mesh_xi)}")


def cmd_report(cfg, a):
    g = geo.plain_strip(a.width_nm)
    csv_path, isw_path = paths(cfg, a.width_nm)
    rows = store.load_rows(csv_path)
    isw = store.load_json(isw_path, {})
    cols = [column(g, rows, isw.get(g.key, {}).get("I_sw"), a.film)]
    if cols[0] is None:
        sim.log("no narrow-strip photon rows yet.")
        return
    if a.compare:
        crow = compare_rows(a.compare)
        check_compare(cfg, crow)
        cols = [column(h, crow, isw_from_rows(crow, h.key), a.film) for h in (geo.PLAIN_WIRE, geo.LONG_ZIPPER)] + cols
        cols = [c for c in cols if c]
    pl.apply_style()
    geoms = [c["g"] for c in cols]
    pl.fig_efficiency(geoms, {c["g"].key: c["vals"] for c in cols}, {c["g"].key: c["met"] for c in cols},
                      {c["g"].key: c["I_sw"] for c in cols},
                      f"Hole geometries vs plain strips, normalised bias  [{cfg.hotspot}, T = {cfg.tmax:g} tau0]",
                      cfg.path(f"narrow{a.width_nm:g}nm_efficiency", "png"), styles={g.key: "--"})

    L = []
    P = L.append
    P("=" * 78)
    P(f"NARROW-STRIP BASELINE   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau0   mesh={cfg.mesh_xi:g} xi   "
      f"narrow width={a.width_nm:g} nm")
    for line in a.film.describe():
        P(line)
    P("=" * 78)
    P(f"{'':<30}" + "".join(f"{c['g'].label:>15}" for c in cols))
    P(f"{'I_sw [uA]':<30}" + "".join(f"{c['I_sw']:>15.1f}" for c in cols))
    P(f"{'narrowest free width [nm]':<30}" + "".join(f"{c['w_free']:>15.0f}" for c in cols))
    P(f"{'I_sw / I_dep(free width)':<30}" + "".join(f"{c['I_sw'] / c['I_dep']:>15.3f}" for c in cols))
    for key, name in [("i_10", "i_10"), ("i_50", "i_50"), ("i_90", "i_90"), ("i_95", "i_95"), ("i_99", "i_99"),
                      ("i_max", "i_max (100 % point)")]:
        P(f"{name:<30}" + "".join(f"{c['met'][key]:>15.3f}" for c in cols))
    for key in ("i_50", "i_95"):
        P(f"{key + '  in units of I_dep':<30}" + "".join(f"{c['met'][key] * c['I_sw'] / c['I_dep']:>15.3f}" for c in cols))
    P(f"{'window below the wall: 1 - i_95':<30}" + "".join(f"{1 - c['met']['i_95']:>15.3f}" for c in cols))
    P(f"{'  (photon-level sites used)':<30}" + "".join(f"{c['n']:>15d}" for c in cols))
    if any(c["sel"] for c in cols):
        P("-" * 78)
        P("selectivity  median[i_det(f) - i_det(1)] over probe sites")
        for f in sorted({f for c in cols for f in c["sel"]}):
            P(f"  f = {f:<5g}".ljust(30) + "".join(f"{an.fmt(c['sel'].get(f)):>15}" for c in cols))
    by = {c["g"].key: c for c in cols}
    if "Long_Zipper" in by:
        z, n = by["Long_Zipper"]["met"], by[g.key]["met"]
        P("-" * 78)
        P("zipper minus narrow strip  (positive = the narrow strip turns on / saturates at the lower bias)")
        P(f"  onset    i_50: {z['i_50'] - n['i_50']:+.3f}      i_10: {z['i_10'] - n['i_10']:+.3f}")
        P(f"  plateau  i_95: {z['i_95'] - n['i_95']:+.3f}      i_99: {z['i_99'] - n['i_99']:+.3f}")
        if z["i_50"] >= n["i_50"] and z["i_95"] >= n["i_95"]:
            P("  RESULT: the narrow strip is at least as good as the zipper on onset AND plateau in normalised bias.")
        else:
            P("  RESULT: the narrow strip does NOT dominate the zipper on both measures.")
    for c in cols:
        if c["flags"]:
            P(f"  {c['g'].label}: flagged rows {c['flags']} -- those thresholds are bounds, not values")
    P("=" * 78)
    text = "\n".join(L)
    sim.log(text)
    with open(cfg.path(f"narrow{a.width_nm:g}nm_summary", "txt"), "w") as fh:
        fh.write(text + "\n")


def cmd_scan(cfg, a):
    cols = []
    for w in a.widths:
        g = geo.plain_strip(w)
        csv_path, isw_path = paths(cfg, w)
        c = column(g, store.load_rows(csv_path), store.load_json(isw_path, {}).get(g.key, {}).get("I_sw"), a.film)
        if c is None:
            sim.log(f"  {w:g} nm: no results for {cfg.tag}, skipped")
            continue
        c["label"] = f"{w:g} nm"
        cols.append(c)
    if a.compare:
        crow = compare_rows(a.compare)
        check_compare(cfg, crow)
        for h, lab in ((geo.PLAIN_WIRE, "200 nm"), (geo.LONG_ZIPPER, "zipper")):
            c = column(h, crow, isw_from_rows(crow, h.key), a.film)
            if c:
                c["label"] = lab
                cols.append(c)
    for c in cols:
        if c["flags"]:
            sim.log(f"  {c['label']}: flagged rows {c['flags']} -- these thresholds are bounds, not values")
    L = []
    P = L.append
    P("=" * 96)
    P(f"WIDTH SCAN   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau0   mesh={cfg.mesh_xi:g} xi   "
      f"I_dep = {a.film.j_dep_uA_per_nm:.3f} uA/nm x narrowest free width (GL)")
    P("=" * 96)
    P(f"{'':<28}" + "".join(f"{c['label']:>11}" for c in cols))
    P(f"{'I_sw [uA]':<28}" + "".join(f"{c['I_sw']:>11.1f}" for c in cols))
    P(f"{'I_sw / I_dep':<28}" + "".join(f"{c['I_sw'] / c['I_dep']:>11.3f}" for c in cols))
    for key in ("i_10", "i_50", "i_95", "i_max"):
        P(f"{key + '  [I_sw]':<28}" + "".join(f"{c['met'][key]:>11.3f}" for c in cols))
    for key in ("i_50", "i_95"):
        P(f"{key + '  [I_dep]':<28}" + "".join(f"{c['met'][key] * c['I_sw'] / c['I_dep']:>11.3f}" for c in cols))
    P(f"{'window 1 - i_95':<28}" + "".join(f"{1 - c['met']['i_95']:>11.3f}" for c in cols))
    P(f"{'(1 - i_50) x width [nm]':<28}" + "".join(f"{(1 - c['met']['i_50']) * c['w_free']:>11.1f}" for c in cols))
    P(f"{'  (photon-level sites)':<28}" + "".join(f"{c['n']:>11d}" for c in cols))
    fs = sorted({f for c in cols for f in c["sel"]})
    if fs:
        P("selectivity  median[i_det(f) - i_det(1)] over probe sites")
        for f in fs:
            P(f"  f = {f:<5g}".ljust(28) + "".join(f"{an.fmt(c['sel'].get(f)):>11}" for c in cols))
    P("=" * 96)
    text = "\n".join(L)
    sim.log(text)
    with open(cfg.path("width_scan_summary", "txt"), "w") as fh:
        fh.write(text + "\n")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report", "scan"])
    ap.add_argument("--width-nm", type=float, default=56.0)
    ap.add_argument("--widths", default="30,56,100", help="scan: comma-separated widths [nm]")
    ap.add_argument("--hotspot", choices=["conserving", "legacy"], default=sim.SimConfig.hotspot)
    ap.add_argument("--tmax", type=float, default=sim.SimConfig.tmax)
    ap.add_argument("--mesh-xi", type=float, default=sim.SimConfig.mesh_xi)
    ap.add_argument("--hit-lo-frac", type=float, default=sim.SimConfig.hit_lo_frac)
    ap.add_argument("--energies", default="0.6,0.35")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--compare", default=None, help="hitmap.py CSV with the same hotspot, window and mesh")
    ap.add_argument("--rho-uohm-cm", type=float, default=None, help="film resistivity; only used to quote times in ps")
    ap.add_argument("--out", default=OUT_DIR)
    a = ap.parse_args(argv)
    a.energies = tuple(float(s) for s in a.energies.split(",") if s)
    a.widths = [float(s) for s in a.widths.split(",") if s]
    a.film = Film(rho_uohm_cm=a.rho_uohm_cm)
    cfg = sim.SimConfig(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi, hit_lo_frac=a.hit_lo_frac, out_dir=a.out)
    return a, cfg


def main(argv=None):
    a, cfg = parse_args(argv)
    {"run": cmd_run, "report": cmd_report, "scan": cmd_scan}[a.command](cfg, a)


if __name__ == "__main__":
    main()
