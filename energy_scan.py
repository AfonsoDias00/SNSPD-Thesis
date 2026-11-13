#!/usr/bin/env python3
"""
energy_scan.py -- thresholds at several photon energies, for the plain wire and the double column (11 Nov 2026).

Question: does a lower-energy photon lose efficiency faster in the plain 200 nm wire than in the slotted wire?
The same model as thresholds.py (2.0 xi screening mesh, circuit with 750 squares, bath 0.3 T_c, count = 20 % of
I_b in the load).  Rows go to the same files as thresholds.py (results/thresholds/thr_n750_w<lambda>_m2.csv), so
`python thresholds.py report --wavelength <lambda> --mesh-xi 2.0` works on them, and a rerun skips finished rows.

All tasks go into one pool, in priority order:
  1. the plain wire, full grid (8 sites across the half width), at every wavelength;
  2. the double column on the MEDIUM grid (23 of the 63 full-grid sites: the critical sites plus rows at
     y = 100, 180, 200, 204, 210, 220 nm, so the interpolation no longer spreads bridge values along the channel),
     at every wavelength;
  3. the rest of the double column's full grid at 1550 nm (the full map; also tells how good the medium grid is).
Dark switching does not depend on the photon and is not rerun.

  python energy_scan.py run --workers 11
  python energy_scan.py plan            (list the tasks, run nothing)
  python energy_scan.py run --wavelengths 1550 --photon absorbed=0.5 --label abs0.5 --no-full
                                        (sensitivity to the deposited fraction: own file, ..._m2_abs0.5.csv)
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from snspd import cli, sim, store
from snspd import geometry as geo

WAVELENGTHS = (1550.0, 2500.0, 4000.0)     # 0.80, 0.50, 0.31 eV
MEDIUM = ((-96, 0), (-80, 0), (-64, 0), (-16, 0), (0, 0),
          (-80, 100), (0, 100),
          (-80, 180), (0, 180),
          (-80, 200), (-16, 200), (0, 200),
          (-56, 204), (-40, 204), (-24, 204),
          (-80, 210), (0, 210),
          (-96, 220), (-80, 220), (-56, 220), (-40, 220), (-16, 220), (0, 220))


def model_args(wavelength, mesh_xi, photon="", label="", thermal=""):
    ap = argparse.ArgumentParser()
    cli.add_model_args(ap)
    argv = ["--wavelength", f"{wavelength:g}", "--mesh-xi", f"{mesh_xi:g}"]
    if photon:
        argv += ["--photon", photon]
    if thermal:
        argv += ["--thermal", thermal]
    if label:
        argv += ["--label", label]
    return ap.parse_args(argv)


def build_tasks(mesh_xi, wavelengths=WAVELENGTHS, photon="", label="", rest_1550=True, thermal=""):
    plain, double = geo.get("Plain_Wire"), geo.get("Double_400")
    full_d = double.hit_grid("full")
    medium = [(float(x), float(y)) for x, y in MEDIUM]
    assert all(s in full_d for s in medium), "medium grid must be a subset of the full grid"
    rest = [s for s in full_d if s not in medium]
    stages = [(plain, w, plain.hit_grid("full")) for w in wavelengths]
    stages += [(double, w, medium) for w in wavelengths]
    if rest_1550:
        stages += [(double, 1550.0, rest)]
    tasks, skipped = [], 0
    for g, w, sites in stages:
        a = model_args(w, mesh_xi, photon, label, thermal)
        cfg, tag = cli.config(a), cli.tag(a)
        path = os.path.join("results/thresholds", f"thr_{tag}.csv")
        done = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"]) for r in store.load_rows(path)}
        for s in sites:
            if store.row_key(g.key, *s) in done:
                skipped += 1
                continue
            tasks.append(dict(kind="threshold", geom=g, cfg=cfg, site=s, tag=tag, path=path, wavelength=w))
    return tasks, skipped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "plan"])
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--mesh-xi", type=float, default=2.0)
    ap.add_argument("--wavelengths", default=",".join(f"{w:g}" for w in WAVELENGTHS), help="comma-separated [nm]")
    ap.add_argument("--photon", default="", help="photon overrides, e.g. absorbed=0.5 (needs --label)")
    ap.add_argument("--thermal", default="", help="heat-model overrides, e.g. latent=False (needs --label)")
    ap.add_argument("--label", default="", help="suffix of the result file for a model variant, e.g. abs0.5")
    ap.add_argument("--no-full", action="store_true", help="skip the rest of the double column's full grid")
    a = ap.parse_args(argv)
    if (a.photon or a.thermal) and not a.label:
        ap.error("--photon/--thermal need --label, or the rows would land in the default model's file")
    wl = tuple(float(w) for w in a.wavelengths.split(","))
    tasks, skipped = build_tasks(a.mesh_xi, wl, a.photon, a.label, rest_1550=not a.no_full, thermal=a.thermal)
    sim.log(f"{len(tasks)} thresholds to run, {skipped} already done")
    for t in tasks:
        sim.log(f"  {t['wavelength']:>6g} nm  {t['geom'].key:<11} {t['site']}")
    if a.command == "plan":
        return
    for line in tasks[0]["cfg"].describe():
        sim.log(line)

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"FAILED {t['wavelength']:g} nm {t['geom'].key} {t['site']}: {r['error']}\n{r.get('tb', '')}")
            return
        row = store.threshold_row(t["geom"].key, t["site"], t["tag"], r["res"], r["wall_s"])
        store.append_row(t["path"], row)
        sim.log(f"  {t['wavelength']:>6g} nm  {t['geom'].key:<11} {t['site']}: {row['i_det']:.3f} {row['flag']}"
                f"  ({row['n_runs']} runs, {row['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result)


if __name__ == "__main__":
    main()
