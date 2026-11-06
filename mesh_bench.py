#!/usr/bin/env python3
"""
mesh_bench.py -- wall time of one photon run per mesh variant, all variants at the same time (same machine load).
(6 Nov 2026, companion of mesh_test.py, whose timings come from different moments of the day.)

One run each: production preset nbtin (gamma 10, settled start), kinetic photon f = 1 at the cell centre site, a
current below threshold (no switching, so every run lasts the full 400 tau0 and no settle runs are added).

  python mesh_bench.py            # writes results/meshtest/bench.csv
"""
from __future__ import annotations

import multiprocessing as mp
import os
import time

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")
os.environ["SNSPD_NO_PROGRESS"] = "1"

from snspd import sim, store
from thermal_check import geometry

OUT = "results/meshtest"
# (geometry, mesh, site, current [uA]): currents ~0.75 I_dep, below every gamma-10 threshold at these sites
VARIANTS = [("Plain_Wire", m, (0.0, 0.0), 75.0) for m in ("P1.4", "U2.0", "U2.5")] + \
           [("Slot_400", m, (-58.0, 0.0), 63.0) for m in ("P1.4", "G1.4/4", "G1.4/2", "U2.0")]


def run(v):
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    from snspd import hotspot as hsm
    from snspd import meshgen
    gk, mesh, (x, y), I = v
    g = geometry(gk)
    cfg = sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", preset="nbtin", sustained=False,
                        out_dir=OUT)
    film, params = cfg.pset.film, cfg.pset.params
    if mesh == "P1.4":
        dev = sim.make_device(g, 1.4, film)
    elif mesh.startswith("U"):
        os.environ["SNSPD_MESH_DIR"] = ""
        dev = sim.make_device(g, float(mesh[1:]), film)
    else:
        fine, coarse = (float(s) for s in mesh[1:].split("/"))
        dev = sim.make_device(g, 1.4, film)
        dev._create_dimensionless_mesh(*meshgen.graded_mesh(dev, g, film.xi_nm, fine, coarse))
    hs = hsm.HotspotSolver(dev.points, dev.triangles, params)
    eps = hsm.build_eps(hs, "kinetic", x, y, 1.0, cfg.tmax, params)
    t = time.time()
    trig = sim.run_simulation(dev, I, eps, cfg, t0=params.t0_tau0)
    return dict(geom=gk, mesh=mesh, nodes=len(dev.points), I_uA=I, switched=trig.switched,
                wall_s=round(time.time() - t))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with mp.get_context("spawn").Pool(len(VARIANTS)) as pool:
        for r in pool.imap_unordered(run, VARIANTS):
            store.append_row(os.path.join(OUT, "bench.csv"), r, list(r))
            print(r, flush=True)
