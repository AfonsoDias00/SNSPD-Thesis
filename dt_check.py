#!/usr/bin/env python3
"""
dt_check.py -- the solver time step: what dt_max = 0.2 did to the results.  (7 Nov 2026)

With dt_max = 0.2 (the production value since the first scripts) pyTDGL's GTDGL update develops a local instability at
single small mesh cells: a site whose |psi| sits above its neighbours and whose potential mu is offset by ~0.2, with a
stationary |grad mu|^2 50x the 99.9th percentile of the rest of the film (dark, settled, no photon).  Found on 7 Nov
through the electro-thermal model, which turns that |grad mu|^2 into a permanent hot spot: the heated snapshot of the
double column latched there, not at the photon.  dt_max 0.05 removes it on the 2.0 xi meshes, 0.025 on the 1.4 xi
double column.

This re-bisects a few thresholds with dt_max = 0.025, without heat (voltage trigger) and with heat (NbTiN 2T), on the
2.0 xi shortened strips, against the stored dt_max = 0.2 values, and times the runs.

  python dt_check.py run --workers 11
  python dt_check.py report

Writes results/dtcheck/check.csv.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import time

os.environ["SNSPD_MESH_DIR"] = ""
os.environ["SNSPD_NO_PROGRESS"] = "1"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from snspd import geometry as geo
from snspd import sim, store

OUT = "results/dtcheck"
DT = 0.025
SITES = [("Plain_Wire", -58.0, 0.0), ("Plain_Wire", 0.0, 0.0), ("Double_400", -80.0, 0.0), ("Double_400", 0.0, 0.0),
         ("Double_400", -40.0, 220.0), ("Double_400", -80.0, 220.0)]
TOL = 0.005


def task_fn(t):
    from dataclasses import replace
    from snspd import hotspot as hsm
    from snspd import thermal
    t_start = time.time()
    try:
        g = geo.get(t["geom"])
        cfg = replace(sim.SimConfig(tmax=400.0, mesh_xi=2.0, hotspot="kinetic", preset="nbtin",
                                    sustained=t["kind"] == "dark", out_dir=OUT), dt_max=t["dt"])
        film, params = cfg.pset.film, cfg.pset.params
        ctx = sim.context(g, cfg.mesh_xi, film, params)
        dev = ctx["device"]
        idep = film.I_dep_uA(g.free_width_min_nm())
        if t["kind"] == "dark":
            I, flag, n, trig = sim.bisect(dev, 0.85 * idep, 1.15 * idep, None, cfg, TOL * idep)
            val = I / idep
        else:
            eps = hsm.build_eps(ctx["hs"], "kinetic", t["x"], t["y"], 1.0, cfg.tmax, params)
            lo, hi, n, flag = 0.55, 1.0, 0, ""
            while hi - lo > TOL:
                mid = 0.5 * (lo + hi)
                n += 1
                if t["kind"] == "heat":
                    hit = thermal.run_thermal(dev, mid * idep, eps, cfg, thermal.NBTIN_2T, params.t0_tau0)["latched"]
                else:
                    hit = sim.run_simulation(dev, mid * idep, eps, cfg, t0=params.t0_tau0).switched
                hi, lo = (mid, lo) if hit else (hi, mid)
            val = hi
        return dict(geom=g.key, kind=t["kind"], x_nm=t.get("x", ""), y_nm=t.get("y", ""), dt_max=t["dt"],
                    i=round(val, 4), n_runs=n, wall_s=round(time.time() - t_start))
    except Exception as exc:
        import traceback
        return dict(error=f"{t}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    tasks = [dict(kind="dark", geom=k, dt=DT) for k in ("Plain_Wire", "Double_400")]
    tasks += [dict(kind=kind, geom=k, x=x, y=y, dt=DT) for kind in ("heat", "trigger") for k, x, y in SITES]
    done = {(r["geom"], r["kind"], r["x_nm"], r["y_nm"]) for r in store.load_rows(f"{OUT}/check.csv")}
    tasks = [t for t in tasks if (t["geom"], t["kind"], str(t.get("x", "")), str(t.get("y", ""))) not in done]
    sim.log(f"-- dt check: {len(tasks)} bisections at dt_max {DT}")
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for r in pool.imap_unordered(task_fn, tasks):
            if "error" in r:
                sim.log("   FAILED " + r["error"])
                continue
            store.append_row(f"{OUT}/check.csv", r, list(r))
            sim.log("   " + "  ".join(f"{k}={v}" for k, v in r.items()))
    cmd_report(a)


def cmd_report(a):
    film = sim.PRESETS["nbtin"].film
    hit = store.load_rows("results/maps/hitmap_kinetic-nbtin_T400_m2_trigger.csv")
    isw = store.load_json("results/maps/isw_kinetic-nbtin_T400_m2_trigger.json", {})
    heat = {(r["geom"], float(r["x_nm"]), float(r["y_nm"])): float(r["i_heat"])
            for r in store.load_rows("results/heat/check.csv") if r["variant"] == "NbTiN 2T"}
    trig = {(r["Geometry"], float(r["x_nm"]), float(r["y_nm"])): float(r["I_det_uA"]) /
            film.I_dep_uA(geo.get(r["Geometry"]).free_width_min_nm()) for r in hit if float(r["energy_f"]) == 1.0}
    L = [f"TIME STEP: dt_max {DT} against the stored dt_max 0.2 (2.0 xi, nbtin, kinetic f = 1, 400 tau0) [I_dep]", "",
         f"{'what':<34}{'dt 0.2':>9}{f'dt {DT}':>9}{'change':>9}{'s/run':>8}"]
    for r in store.load_rows(f"{OUT}/check.csv"):
        k = r["geom"]
        if r["kind"] == "dark":
            old = isw.get(k, {}).get("I_sw", float("nan")) / film.I_dep_uA(geo.get(k).free_width_min_nm())
            what = f"{k} dark I_sw (sustained)"
        else:
            key = (k, float(r["x_nm"]), float(r["y_nm"]))
            old = (heat if r["kind"] == "heat" else trig).get(key, float("nan"))
            what = f"{k} ({float(r['x_nm']):g},{float(r['y_nm']):g}) {r['kind']}"
        new = float(r["i"])
        L.append(f"{what:<34}{old:9.3f}{new:9.3f}{new - old:+9.3f}{float(r['wall_s']) / max(int(r['n_runs']), 1):8.0f}")
    text = "\n".join(L)
    print(text)
    with open(f"{OUT}/summary.txt", "w") as fh:
        fh.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--workers", type=int, default=11)
    a = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
