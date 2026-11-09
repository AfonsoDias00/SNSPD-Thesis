#!/usr/bin/env python3
"""
heat_check.py -- does Joule heat turn the slot layouts' cascade bursts into detections?  (7 Nov 2026)

Without heat the two criteria are far apart in the slot layouts.  In the double column, a channel photon fires the
voltage trigger from 0.67 I_dep, but the strip only stays resistive from 0.80-0.83: in between, the hit channel dumps
its current into the others with 4-5 phase slips and the device heals with trapped flux (docs/slots.md, sec. 6-7).
Each slip dissipates I Phi0, about 0.5 eV here, so heat may hold a normal domain.

This bisects the threshold with the electro-thermal TDGL (snspd/thermal.py).  Detection is one criterion: the strip
is still resistive at the end of the 400 tau0 window.  Preset nbtin (gamma 10, settled start), the kinetic hotspot
at f = 1, 2.0 xi, shortened strips.  Heat variants: the NbTiN two-temperature model with the measured parameters
(thermal.NBTIN_2T) everywhere; the one-temperature model (eta 0.59) at two sites as a check.

  python heat_check.py run --workers 10
  python heat_check.py run --workers 11 --grid-geoms Plain_Wire,Double_400,Stagger_400,Single_400
  python heat_check.py run --workers 11 --full-geoms Double_400,Plain_Wire
  python heat_check.py report

Writes results/heat/check.csv.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import time

os.environ["SNSPD_MESH_DIR"] = ""
os.environ["SNSPD_NO_PROGRESS"] = "1"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")         # before numpy is imported: one BLAS thread per worker process

from snspd import geometry as geo
from snspd import sim, store

OUT = "results/heat"
TOL = 0.005                         # [I_dep of the free width]
BRACKET = (0.55, 1.00)
SITES = [("Plain_Wire", -58.0, 0.0), ("Plain_Wire", 0.0, 0.0),
         ("Double_400", -80.0, 0.0), ("Double_400", 0.0, 0.0), ("Double_400", -40.0, 220.0),
         ("Stagger_400", -80.0, 0.0), ("Stagger_400", -40.0, 220.0),
         ("Single_400", -70.0, 0.0), ("Single_400", -44.0, 0.0)]
ONE_T_SITES = (("Plain_Wire", -58.0, 0.0), ("Double_400", -80.0, 0.0))


def variants():
    from snspd import thermal
    return {"NbTiN 2T": thermal.NBTIN_2T, "1T eta 0.59": thermal.ThermalParams(eta=0.59, D_nm2_per_tau0=7.1,
                                                                              esc_tau0=530.0)}


def task_fn(t):
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    from snspd import hotspot as hsm
    from snspd import thermal
    t_start = time.time()
    try:
        g = geo.get(t["geom"])
        from dataclasses import replace
        cfg = replace(sim.SimConfig(tmax=400.0, mesh_xi=t.get("mesh_xi", 2.0), hotspot="kinetic", preset="nbtin",
                                    sustained=False, out_dir=OUT), dt_max=t.get("dt_max", sim.SimConfig.dt_max))
        film, params = cfg.pset.film, cfg.pset.params
        ctx = sim.context(g, cfg.mesh_xi, film, params)
        dev = ctx["device"]
        idep = film.I_dep_uA(g.free_width_min_nm())
        eps = hsm.build_eps(ctx["hs"], "kinetic", t["x"], t["y"], 1.0, cfg.tmax, params)
        tp = variants()[t["variant"]]
        lo, hi, n, last = BRACKET[0], BRACKET[1], 0, None
        while hi - lo > TOL:
            mid = 0.5 * (lo + hi)
            r = thermal.run_thermal(dev, mid * idep, eps, cfg, tp, params.t0_tau0)
            n += 1
            if r["latched"]:
                hi, last = mid, r
            else:
                lo = mid
        flag = "at_wall" if hi >= BRACKET[1] else ("below_bracket" if lo <= BRACKET[0] else "")
        last = last or {}
        return dict(geom=g.key, x_nm=t["x"], y_nm=t["y"], variant=t["variant"], i_heat=round(hi, 4), flag=flag,
                    n_runs=n, branch=last.get("branch", ""), t_trigger=last.get("t_trigger", ""),
                    theta_peak=round(last.get("theta_peak", 0.0), 3), normal_nm2=round(last.get("normal_nm2", 0.0)),
                    wall_s=round(time.time() - t_start))
    except Exception as exc:
        import traceback
        return dict(error=f"{t}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")


def check_csv(mesh_xi, dt_max=None):
    """results/heat/check[_m<mesh>][_dt<dt>].csv; the files without _dt are dt_max 0.2 (before 7 Nov 2026 15:00,
    affected by the local instability of docs/methods.md)."""
    dt = sim.SimConfig.dt_max if dt_max is None else dt_max
    return (f"{OUT}/check" + ("" if mesh_xi == 2.0 else f"_m{mesh_xi:g}") + ("" if dt == 0.2 else f"_dt{dt:g}") + ".csv")


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    csv = check_csv(a.mesh_xi, a.dt_max)
    done = {(r["geom"], float(r["x_nm"]), float(r["y_nm"]), r["variant"]) for r in store.load_rows(csv)}
    sites = [] if (a.only_grids or a.extra) else list(SITES)
    if a.grid_geoms:             # every critical-grid site of these geometries as well (NbTiN 2T)
        for k in a.grid_geoms.split(","):
            sites += [(k, float(x), float(y)) for x, y in geo.get(k).grid_critical]
    if a.full_geoms:             # every full-grid site (NbTiN 2T): a heated map
        for k in a.full_geoms.split(","):
            sites += [(k, float(x), float(y)) for x, y in geo.get(k).hit_grid("full")
]
    for spec in filter(None, a.extra.split(";")):     # GEOM:x,y
        k, xy = spec.split(":")
        sites.append((k, *(float(v) for v in xy.split(","))))
    sites = list(dict.fromkeys(sites))
    tasks = [dict(geom=k, x=x, y=y, variant="NbTiN 2T") for k, x, y in sites] + \
            [dict(geom=k, x=x, y=y, variant="1T eta 0.59") for k, x, y in ONE_T_SITES]
    if a.only_grids or a.extra:
        tasks = [t for t in tasks if t["variant"] == "NbTiN 2T"]
    tasks = [dict(t, mesh_xi=a.mesh_xi, dt_max=a.dt_max) for t in tasks if (t["geom"], t["x"], t["y"], t["variant"]) not in done]
    sim.log(f"-- heat check: {len(tasks)} bisections")
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for r in pool.imap_unordered(task_fn, tasks):
            if "error" in r:
                sim.log("   FAILED " + r["error"])
                continue
            store.append_row(csv, r, list(r))
            sim.log("   " + "  ".join(f"{k}={v}" for k, v in r.items()))
    cmd_report(a)


def cmd_report(a):
    rows = store.load_rows(check_csv(2.0, a.dt_max))
    film = sim.PRESETS["nbtin"].film
    hit = store.load_rows("results/maps/hitmap_kinetic-nbtin_T400_m2_trigger.csv")
    trig = {(r["Geometry"], float(r["x_nm"]), float(r["y_nm"])): float(r["I_det_uA"]) /
            film.I_dep_uA(geo.get(r["Geometry"]).free_width_min_nm()) for r in hit if float(r["energy_f"]) == 1.0}
    sust = {}
    for k in ("Double_400", "Stagger_400", "Single_400"):
        for r in store.load_rows(f"results/flux/scan_{k}_kinetic-nbtin_T400_m2.csv"):
            if r["outcome"] in ("direct", "via_trap"):
                key = (k, float(r["x_nm"]), float(r["y_nm"]))
                sust[key] = min(sust.get(key, 9.0), float(r["bias"]))
    L = ["HEAT CHECK: thresholds [I_dep of the free width] with Joule heat (still resistive at 400 tau0) against the "
         "two bounds without heat", "(trigger: hitmap 2.0 xi; sustained: flux scan grid, 0.03 steps; plain wire: "
         "the bounds coincide)", "",
         f"{'site':<26}{'variant':<13}{'with heat':>10}{'trigger':>9}{'sustained':>11}  peak theta, normal area at end"]
    for r in sorted(rows, key=lambda r: (r["geom"], float(r["y_nm"]), float(r["x_nm"]), r["variant"])):
        key = (r["geom"], float(r["x_nm"]), float(r["y_nm"]))
        s = sust.get(key)
        L.append(f"{r['geom'] + ' (' + format(float(r['x_nm']), 'g') + ',' + format(float(r['y_nm']), 'g') + ')':<26}"
                 f"{r['variant']:<13}{float(r['i_heat']):10.3f}{trig.get(key, float('nan')):9.3f}"
                 f"{('' if s is None else format(s, '.2f')):>11}  {r['theta_peak']}, {r['normal_nm2']} nm^2"
                 + (f"  [{r['flag']}]" if r["flag"] else ""))
    text = "\n".join(L)
    print(text)
    with open(check_csv(2.0, a.dt_max).replace("check", "summary").replace(".csv", ".txt"), "w") as fh:
        fh.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--grid-geoms", default="", help="also every critical-grid site of these geometries (comma separated)")
    ap.add_argument("--full-geoms", default="", help="also every full-grid site of these geometries (a heated map)")
    ap.add_argument("--mesh-xi", type=float, default=2.0, help="1.4: writes results/heat/check_m1.4.csv")
    ap.add_argument("--extra", default="", help="more sites, GEOM:x,y;GEOM:x,y (NbTiN 2T)")
    ap.add_argument("--dt-max", type=float, default=sim.SimConfig.dt_max, help="solver dt cap (0.2 = the old runs)")
    ap.add_argument("--only-grids", action="store_true", help="only the --grid-geoms / --full-geoms sites, NbTiN 2T")
    a = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
