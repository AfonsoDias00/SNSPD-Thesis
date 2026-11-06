#!/usr/bin/env python3
"""
mesh_test.py -- how coarse can the mesh be before the photon thresholds stop making sense?  (6 Nov 2026)

Photon thresholds at the sites of thermal_check.py (kinetic hotspot f = 1, 400 tau0, gamma 1: the stored production
parameter set, so the 1.4 xi values are the reference) on other meshes, both criteria, each mesh with its own dark
I_sw (bisected here; photon bracket [0.30, 1] x I_sw as in production):
  U2.0, U2.5   uniform, pyTDGL's mesher with max edge 2.0 / 2.5 xi (Simon et al. use 2-2.5 xi)
  G1.4/4       graded (snspd/meshgen.py): production density within 150 nm of the photon cell and 60 nm of the contacts,
               4 xi elsewhere -- only the far field is coarser
  G2.0/4       graded, 2.0 xi near the photon cell
The 0.8 um devices get the uniform meshes (grading saves only ~15 % there), Slot 400 (2 um) the graded ones.

  python mesh_test.py run --workers 4
  python mesh_test.py report
  python mesh_test.py run --preset nbtin --meshes U2.0 --criteria sustained     # the gamma-10 spot check

Writes results/meshtest[_<preset>]/check.csv (thresholds and wall times), isw.csv, summary.txt.  With a preset the
reference is the production 1.4 xi map of that preset (factorial_check.map_value) instead of thermal_check.CASES.
Same-load wall times per mesh: mesh_bench.py.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")

from snspd import sim, store
from thermal_check import CASES, geometry

OUT = "results/meshtest"
CRITERIA = ("trigger", "sustained")
PLAN = {"Plain_Wire": ("U2.0", "U2.5"), "Long_Zipper": ("U2.0", "U2.5"), "Slot_400": ("G1.4/4", "G2.0/4")}
SITES = ("plain (0,0)", "plain (-96,0)", "zipper (0,0)", "zipper (-68,0)", "zipper (-42,58)", "zipper (-36,61)",
         "Slot 400 (-58,0)", "Slot 400 (-28,200)")
FIELDS = ["mesh", "label", "geom", "criterion", "nodes", "i_det", "flag", "n_runs", "t_trigger", "I_sw_uA", "wall_s"]
ISW_FIELDS = ["mesh", "geom", "nodes", "I_sw_uA", "i_sw", "flag", "n_runs", "wall_s"]

_DEV = {}
PRESET = ""          # sim preset of the run (set by main; worker processes get it through the task)


def device(gk, mesh):
    """(device, hotspot solver) of geometry gk on mesh variant `mesh`, cached per process."""
    key = (gk, mesh)
    if key not in _DEV:
        from snspd import hotspot as hsm
        from snspd import meshgen
        g = geometry(gk)
        film, params = sim.PRESETS[PRESET].film, sim.PRESETS[PRESET].params
        if mesh.startswith("U"):
            os.environ["SNSPD_MESH_DIR"] = ""               # a fresh pyTDGL mesh (no pinned file at this edge)
            dev = sim.make_device(g, float(mesh[1:]), film)
            os.environ["SNSPD_MESH_DIR"] = "meshes"
        else:
            fine, coarse = (float(v) for v in mesh[1:].split("/"))
            dev = sim.make_device(g, 1.4, film)
            pts, tri = meshgen.graded_mesh(dev, g, film.xi_nm, fine, coarse)
            dev._create_dimensionless_mesh(pts, tri)
        _DEV[key] = (dev, hsm.HotspotSolver(dev.points, dev.triangles, params))
    return _DEV[key]


def task_fn(task):
    global PRESET
    from snspd import hotspot as hsm
    PRESET = task["preset"]
    t_start = time.time()
    try:
        g = geometry(task["geom"])
        dev, hs = device(task["geom"], task["mesh"])
        cfg = sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", sustained=task["criterion"] == "sustained",
                            preset=task["preset"], out_dir=OUT)
        film = cfg.pset.film
        tol, (lo, hi) = sim.scales(g, film)
        eps = None
        if task["kind"] == "hit":
            lo, hi = cfg.hit_lo_frac * task["I_sw_uA"], task["I_sw_uA"]
            eps = hsm.build_eps(hs, "kinetic", task["x"], task["y"], 1.0, cfg.tmax, cfg.pset.params)
        I, flag, n, trig = sim.bisect(dev, lo, hi, eps, cfg, tol)
        return dict(ok=True, task=task, I=I, flag=flag, n=n, t=trig.t, nodes=len(dev.points),
                    i_dep=film.I_dep_uA(g.free_width_min_nm()), wall_s=time.time() - t_start)
    except Exception as exc:
        import traceback
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}")


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    done = {(r["mesh"], r["label"], r["criterion"]) for r in store.load_rows(os.path.join(OUT, "check.csv"))}
    isw = {(r["mesh"], r["geom"]): float(r["I_sw_uA"]) for r in store.load_rows(os.path.join(OUT, "isw.csv"))}
    cases = {c[0]: c for c in CASES if c[0] in SITES}
    waiting, first = {}, []
    for label, (_, gk, x, y, _, _) in cases.items():
        for mesh in PLAN[gk]:
            if a.meshes and mesh not in a.meshes:
                continue
            for crit in a.criteria:
                if (mesh, label, crit) in done:
                    continue
                t = dict(kind="hit", mesh=mesh, label=label, geom=gk, x=float(x), y=float(y), criterion=crit,
                         preset=a.preset)
                if (mesh, gk) in isw:
                    first.append(dict(t, I_sw_uA=isw[(mesh, gk)]))
                else:
                    if (mesh, gk) not in waiting:
                        first.append(dict(kind="isw", mesh=mesh, geom=gk, criterion="sustained", label="dark",
                                          preset=a.preset))
                    waiting.setdefault((mesh, gk), []).append(t)
    first.sort(key=lambda t: t["kind"] != "isw")
    sim.log(f"-- mesh test: {sum(t['kind'] == 'isw' for t in first)} dark I_sw, "
            f"{sum(t['kind'] == 'hit' for t in first) + sum(map(len, waiting.values()))} thresholds")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t}: {r['error']}")
            return []
        if t["kind"] == "isw":
            store.append_row(os.path.join(OUT, "isw.csv"), dict(
                mesh=t["mesh"], geom=t["geom"], nodes=r["nodes"], I_sw_uA=round(r["I"], 3),
                i_sw=round(r["I"] / r["i_dep"], 4), flag=r["flag"], n_runs=r["n"], wall_s=round(r["wall_s"])), ISW_FIELDS)
            sim.log(f"   dark {t['mesh']} {t['geom']}: {r['I'] / r['i_dep']:.3f} I_dep ({r['nodes']} nodes)")
            return [dict(s, I_sw_uA=r["I"]) for s in waiting.pop((t["mesh"], t["geom"]), [])]
        store.append_row(os.path.join(OUT, "check.csv"), dict(
            mesh=t["mesh"], label=t["label"], geom=t["geom"], criterion=t["criterion"], nodes=r["nodes"],
            i_det=round(r["I"] / r["i_dep"], 4), flag=r["flag"], n_runs=r["n"],
            t_trigger="" if r["t"] is None else round(r["t"], 1), I_sw_uA=round(t["I_sw_uA"], 3),
            wall_s=round(r["wall_s"])), FIELDS)
        sim.log(f"   {t['mesh']:<7} {t['label']:<20} {t['criterion']:<9} {r['I'] / r['i_dep']:.3f}  {r['wall_s']:.0f} s")
        return []

    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn")) as pool:
        pending = {pool.submit(task_fn, t) for t in first}
        while pending:
            fin, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in fin:
                for nt in on_result(fut.result()):
                    pending.add(pool.submit(task_fn, nt))
    cmd_report(a)


def cmd_report(a):
    rows = store.load_rows(os.path.join(OUT, "check.csv"))
    isw = store.load_rows(os.path.join(OUT, "isw.csv"))
    if a.preset:
        from factorial_check import map_value
        ref = {(c[0], crit): map_value(c[1], c[2], c[3], crit) for c in CASES for crit in CRITERIA}
        ref = {k: (float("nan") if v is None else v) for k, v in ref.items()}
    else:
        ref = {(c[0], "trigger"): c[4] for c in CASES} | {(c[0], "sustained"): c[5] for c in CASES}
    fac = {(r["label"], r["criterion"]): float(r["wall_s"]) / max(int(r["n_runs"]), 1)
           for r in store.load_rows("results/factorial/check.csv") if r["corner"] == "base"}
    L = [f"MESH TEST (preset '{a.preset or 'gamma 1'}', kinetic f = 1, 400 tau0): thresholds [I_dep] against the "
         "1.4 xi production values",
         "s/run: wall time per bisection run under the same machine load (1.4 xi: factorial_check base corner)", "",
         "dark I_sw / I_dep (production 1.4 xi: " + ("plain 1.007, zipper 0.995, Slot 400 ~1.00):" if a.preset == "nbtin"
                                                     else "plain 0.970, zipper 0.968, Slot 400 0.973):")]
    for r in isw:
        L.append(f"  {r['mesh']:<8}{r['geom']:<14}{float(r['i_sw']):.3f}   {r['nodes']} nodes")
    L += ["", f"{'site':<20}{'criterion':<11}{'1.4 xi':>8}" + "".join(f"{m:>16}" for m in ("U2.0", "U2.5", "G1.4/4", "G2.0/4"))
          + "     s/run at 1.4 xi"]
    by = {(r["mesh"], r["label"], r["criterion"]): r for r in rows}
    for label in SITES:
        for crit in CRITERIA:
            cells = []
            for m in ("U2.0", "U2.5", "G1.4/4", "G2.0/4"):
                r = by.get((m, label, crit))
                cells.append(f"{float(r['i_det']):7.3f} ({float(r['i_det']) - ref[(label, crit)]:+.3f})"
                             if r else f"{'':>16}")
            spr = [f"{float(by[(m, label, crit)]['wall_s']) / int(by[(m, label, crit)]['n_runs']):.0f}"
                   for m in ("U2.0", "U2.5", "G1.4/4", "G2.0/4") if (m, label, crit) in by]
            L.append(f"{label:<20}{crit:<11}{ref[(label, crit)]:8.3f}" + "".join(cells)
                     + f"     {fac.get((label, crit), float('nan')):.0f}  vs  {' / '.join(spr)}")
    # the comparisons the study actually makes
    L += ["", "differences that carry the conclusions (zipper - plain at matching sites, I_dep):"]
    for crit in CRITERIA:
        for m in ("1.4", "U2.0", "U2.5"):
            def v(lab):
                return ref[(lab, crit)] if m == "1.4" else (float(by[(m, lab, crit)]["i_det"]) if (m, lab, crit) in by else None)
            pairs = [("zipper (0,0)", "plain (0,0)"), ("zipper (-68,0)", "plain (-96,0)"), ("zipper (-36,61)", "plain (0,0)")]
            vals = [(z, p, v(z), v(p)) for z, p in pairs]
            L.append(f"  {crit:<10}{m:<6}" + "   ".join(f"{z.split()[1]} - {p.split()[1]}: "
                     + (f"{vz - vp:+.3f}" if vz is not None and vp is not None else "  -  ") for z, p, vz, vp in vals))
    text = "\n".join(L)
    print(text)
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--preset", default="", choices=list(sim.PRESETS))
    ap.add_argument("--meshes", nargs="*", default=None, help="subset of the meshes in PLAN")
    ap.add_argument("--criteria", nargs="*", default=list(CRITERIA), choices=list(CRITERIA))
    a = ap.parse_args()
    global OUT, PRESET
    PRESET = a.preset
    if a.preset:
        OUT = f"results/meshtest_{a.preset}"
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
