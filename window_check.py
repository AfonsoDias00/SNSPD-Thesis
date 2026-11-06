#!/usr/bin/env python3
"""
window_check.py -- are the gamma-10 ("nbtin") thresholds set by the 400 tau0 window?  (6 Nov 2026)

At gamma 10 most threshold runs of the sustained maps fire at t = 330-395 tau0, against 85-210 tau0 at gamma 1, and
the NbTiN hotspot (escape 530 tau0) is still there at the end of the window.  This bisects a few sites again with a
longer window (production mesh and preset, voltage trigger: the bounds coincide at gamma 10), in a narrow bracket
[LO, 1] x the 400 tau0 sustained map value.

  python window_check.py --tmax 800 --workers 4      -> results/window/check.csv, summary on stdout
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import time

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")
os.environ["SNSPD_NO_PROGRESS"] = "1"

from snspd import sim, store
from snspd.geometry import get

OUT = "results/window"
LO = 0.80
SITES = [("Plain_Wire", 0.0, 0.0), ("Plain_Wire", -96.0, 0.0), ("Long_Zipper", 0.0, 0.0), ("Long_Zipper", -36.0, 61.0)]
MAP = "results/hitmap/hitmap_kinetic-nbtin_T400_m1.4.csv"


def ref_uA(gk, x, y):
    for r in store.load_rows(MAP):
        if r["Geometry"] == gk and float(r["x_nm"]) == x and float(r["y_nm"]) == y and float(r["energy_f"]) == 1.0:
            return float(r["I_det_uA"]), float(r["t_trigger"])
    raise KeyError((gk, x, y))


def run(task):
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    from snspd import hotspot as hsm
    gk, x, y, tmax, criterion = task
    g = get(gk)
    cfg = sim.SimConfig(tmax=tmax, mesh_xi=1.4, hotspot="kinetic", preset="nbtin",
                        sustained=criterion == "sustained", out_dir=OUT)
    film, params = cfg.pset.film, cfg.pset.params
    ctx = sim.context(g, 1.4, film, params)
    I400, t400 = ref_uA(gk, x, y)
    tol, _ = sim.scales(g, film)
    eps = hsm.build_eps(ctx["hs"], "kinetic", x, y, 1.0, tmax, params)
    t = time.time()
    I, flag, n, trig = sim.bisect(ctx["device"], LO * I400, I400 * 1.0001, eps, cfg, tol)
    idep = film.I_dep_uA(g.free_width_min_nm())
    return dict(geom=gk, x_nm=x, y_nm=y, tmax=tmax, criterion=criterion, i_400=round(I400 / idep, 4),
                t_400=t400, i_det=round(I / idep, 4), change=round((I - I400) / idep, 4), flag=flag, n_runs=n,
                branch=trig.branch, t_trigger="" if trig.t is None else round(trig.t, 1), wall_s=round(time.time() - t))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tmax", type=float, default=800.0)
    ap.add_argument("--criterion", default="trigger", choices=["trigger", "sustained"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    tasks = [(gk, x, y, a.tmax, a.criterion) for gk, x, y in SITES]
    with mp.get_context("spawn").Pool(a.workers) as pool:
        for r in pool.imap_unordered(run, tasks):
            store.append_row(os.path.join(OUT, "check.csv"), r, list(r))
            print(r, flush=True)
