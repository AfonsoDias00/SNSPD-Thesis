#!/usr/bin/env python3
"""
mesh_check.py -- how short and how coarse can the slot devices be?  (6 Nov 2026)

Variants of a slot layout (snspd.geometry.slot_layout), production preset nbtin (gamma 10, settled start):
  2um/U1.4   2 um strip, pyTDGL mesh at max edge 1.4 xi              (the reference)
  2um/G1.4   2 um strip, graded: 1.4 xi over the slotted section, 4 xi in the plain leads (snspd.meshgen)
  short/U1.4 shortest strip: 150 nm of plain lead beyond the outermost slot end
  short/U2.0 shortest strip, 2.0 xi (screening)

  bench      one sub-threshold photon run per variant, all at the same time (same machine load): wall time
  accuracy   dark I_sw (sustained) and voltage-trigger photon thresholds at a few sites per variant, in a fixed bracket

  python mesh_check.py bench    --layout double --slots 150,400
  python mesh_check.py accuracy --layout double --slots 400 --workers 8
  python mesh_check.py report

Writes results/meshcheck/bench.csv, accuracy.csv.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import time

os.environ["SNSPD_NO_PROGRESS"] = "1"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")         # before numpy is imported: one BLAS thread per worker process
os.environ["SNSPD_MESH_DIR"] = ""            # fresh meshes: nothing is pinned for the new layouts yet

from snspd import geometry as geo
from snspd import sim, store

OUT = "results/meshcheck"
VARIANTS = ("2um/U1.4", "2um/G1.4", "short/U1.4", "short/U2.0")
BENCH_BIAS = 0.75                            # [I_dep of the free width], below threshold
DARK_BRACKET = (0.85, 1.15)                  # [I_dep of the free width]
HIT_BRACKET = (0.45, 1.02)
TOL_UA = 0.20                                # 0.0033 I_dep of the 120 nm free width


def sites(g):
    b = g.cell_half_y_nm - geo.BRIDGE_NM / 2.0      # half the slot length
    if g.key.startswith("Single"):
        return [(-70.0, 0.0), (-44.0, 0.0), (0.0, b + 20.0)]
    return [(-80.0, 0.0), (0.0, 0.0), (-40.0, b + 20.0)]


def build(layout, slot, variant):
    from snspd import hotspot as hsm
    from snspd import meshgen
    length, mesh = variant.split("/")
    g = geo.slot_layout(layout, slot, 2000.0 if length == "2um" else None)
    cfg = sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", preset="nbtin", sustained=False, out_dir=OUT)
    film, params = cfg.pset.film, cfg.pset.params
    if mesh == "G1.4":
        dev = sim.make_device(g, 1.4, film)
        dev._create_dimensionless_mesh(*meshgen.graded_mesh(dev, g, film.xi_nm, 1.4, 4.0))
    else:
        dev = sim.make_device(g, float(mesh[1:]), film)
    return g, dev, hsm.HotspotSolver(dev.points, dev.triangles, params), cfg


_CACHE = {}


def task_fn(t):
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    from dataclasses import replace
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        key = (t["layout"], t["slot"], t["variant"])
        if key not in _CACHE:
            _CACHE[key] = build(*key)
        g, dev, hs, cfg = _CACHE[key]
        idep = cfg.pset.film.I_dep_uA(g.free_width_min_nm())
        base = dict(layout=t["layout"], slot=t["slot"], variant=t["variant"], length_nm=g.length_nm,
                    nodes=len(dev.points))
        if t["kind"] == "bench":
            x, y = sites(g)[0]
            eps = hsm.build_eps(hs, "kinetic", x, y, 1.0, cfg.tmax, cfg.pset.params)
            t1 = time.time()
            trig = sim.run_simulation(dev, BENCH_BIAS * idep, eps, cfg, t0=cfg.pset.params.t0_tau0)
            return dict(base, x_nm=x, y_nm=y, bias=BENCH_BIAS, switched=trig.switched, wall_s=round(time.time() - t1))
        if t["kind"] == "dark":
            c = replace(cfg, sustained=True)
            I, flag, n, trig = sim.bisect(dev, DARK_BRACKET[0] * idep, DARK_BRACKET[1] * idep, None, c, TOL_UA)
            x = y = ""
        else:
            x, y = t["x"], t["y"]
            eps = hsm.build_eps(hs, "kinetic", x, y, 1.0, cfg.tmax, cfg.pset.params)
            I, flag, n, trig = sim.bisect(dev, HIT_BRACKET[0] * idep, HIT_BRACKET[1] * idep, eps, cfg, TOL_UA)
        return dict(base, kind=t["kind"], x_nm=x, y_nm=y, i=round(I / idep, 4), flag=flag, n_runs=n,
                    branch=trig.branch, t_trigger="" if trig.t is None else round(trig.t, 1),
                    wall_s=round(time.time() - t_start))
    except Exception as exc:
        import traceback
        return dict(error=f"{t}: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")


def run_pool(tasks, workers, csv):
    os.makedirs(OUT, exist_ok=True)
    with mp.get_context("spawn").Pool(workers, maxtasksperchild=None) as pool:
        for r in pool.imap_unordered(task_fn, tasks):
            if "error" in r:
                sim.log("   FAILED " + r["error"])
                continue
            store.append_row(os.path.join(OUT, csv), r, list(r))
            sim.log("   " + "  ".join(f"{k}={v}" for k, v in r.items()))


def cmd_bench(a):
    tasks = [dict(kind="bench", layout=a.layout, slot=s, variant=v) for s in a.slots for v in VARIANTS]
    run_pool(tasks, len(tasks), "bench.csv")


def cmd_accuracy(a):
    done = {(r["layout"], float(r["slot"]), r["variant"], r["kind"], r["x_nm"], r["y_nm"])
            for r in store.load_rows(os.path.join(OUT, "accuracy.csv"))}
    tasks = []
    for s in a.slots:
        g = geo.slot_layout(a.layout, s)
        for v in VARIANTS:
            todo = [dict(kind="dark")] + [dict(kind="hit", x=x, y=y) for x, y in sites(g)]
            for t in todo:
                k = (a.layout, s, v, t["kind"], str(t.get("x", "")), str(t.get("y", "")))
                if k not in done:
                    tasks.append(dict(t, layout=a.layout, slot=s, variant=v))
    tasks.sort(key=lambda t: t["variant"].startswith("short"))       # the slow 2 um ones first
    sim.log(f"-- {len(tasks)} bisections")
    run_pool(tasks, a.workers, "accuracy.csv")


def cmd_report(a):
    L = []
    rows = store.load_rows(os.path.join(OUT, "bench.csv"))
    if rows:
        L.append("BENCH: one sub-threshold photon run (0.75 I_dep, 400 tau0, nbtin), all variants at once")
        L.append(f"{'device':<16}{'variant':<12}{'length':>8}{'nodes':>8}{'wall s':>8}{'speed-up':>10}")
        for s in sorted({float(r["slot"]) for r in rows}):
            rs = [r for r in rows if float(r["slot"]) == s]
            ref = next((float(r["wall_s"]) for r in rs if r["variant"] == "2um/U1.4"), None)
            for v in VARIANTS:
                for r in rs:
                    if r["variant"] == v:
                        sp = f"{ref / float(r['wall_s']):.1f}x" if ref else ""
                        L.append(f"{r['layout'] + '_' + format(s, 'g'):<16}{v:<12}{float(r['length_nm']):8.0f}"
                                 f"{int(r['nodes']):8d}{float(r['wall_s']):8.0f}{sp:>10}")
    rows = store.load_rows(os.path.join(OUT, "accuracy.csv"))
    if rows:
        L += ["", "ACCURACY [I_dep of the 120 nm free width]: dark I_sw (sustained) and trigger thresholds, "
              "difference to 2um/U1.4"]
        keys = sorted({(r["layout"], float(r["slot"]), r["kind"], r["x_nm"], r["y_nm"]) for r in rows},
                      key=lambda k: (k[0], k[1], k[2] != "dark", k[3], k[4]))
        L.append(f"{'device':<16}{'what':<18}" + "".join(f"{v:>18}" for v in VARIANTS))
        by = {(r["layout"], float(r["slot"]), r["kind"], r["x_nm"], r["y_nm"], r["variant"]): r for r in rows}
        for k in keys:
            ref = by.get(k + ("2um/U1.4",))
            cells = []
            for v in VARIANTS:
                r = by.get(k + (v,))
                if r is None:
                    cells.append(f"{'':>18}")
                    continue
                d = f" ({float(r['i']) - float(ref['i']):+.3f})" if ref and v != "2um/U1.4" else ""
                cells.append(f"{float(r['i']):.3f}{d}{('*' if r['flag'] else ''):1}".rjust(18))
            what = "dark I_sw" if k[2] == "dark" else f"({float(k[3]):g},{float(k[4]):g})"
            L.append(f"{k[0] + '_' + format(k[1], 'g'):<16}{what:<18}" + "".join(cells))
        L.append("* flag: at_wall / below_bracket")
    text = "\n".join(L)
    print(text)
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["bench", "accuracy", "report"])
    ap.add_argument("--layout", default="double", choices=geo.LAYOUTS)
    ap.add_argument("--slots", default="150,400", type=lambda s: [float(v) for v in s.split(",")])
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    {"bench": cmd_bench, "accuracy": cmd_accuracy, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
