#!/usr/bin/env python3
"""
thresholds.py -- the lowest bias at which a photon at each site makes a count, and the dark switching current.

  python thresholds.py run --geoms Plain_Wire,Double_400 --grid critical --mesh-xi 2.0 --workers 11
  python thresholds.py run --geoms Double_400 --grid full --workers 11          (63 sites)
  python thresholds.py run ... --source                    current source instead of the circuit
  python thresholds.py run ... --thermal latent=False --label nolatent          a model variant, kept in its own file
  python thresholds.py report --geoms Plain_Wire,Double_400 --mesh-xi 2.0

Every threshold is a bisection in [0.50, 1.00] I_dep of the free width to 0.005 (7 events).  Rows are appended to
results/thresholds/thr_<tag>.csv as they finish, and a rerun skips them.  The row without a site is the dark run:
the lowest bias with no steady superconducting state.  `report` prints the table and, where a geometry has enough
sites, draws the threshold map and the hard-step efficiency.
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")         # before numpy is imported: one BLAS thread per worker process

from snspd import analysis as an
from snspd import cli, sim, store
from snspd import geometry as geo


def csv_path(a):
    return os.path.join(a.out, f"thr_{cli.tag(a)}.csv")


def cmd_run(a):
    cfg, tag, path = cli.config(a), cli.tag(a), csv_path(a)
    for line in cfg.describe():
        sim.log(line)
    done = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"]) for r in store.load_rows(path)}
    tasks = []
    for key in a.geoms.split(","):
        g = geo.get(key)
        sites = ([] if a.no_dark else [None]) + g.hit_grid(a.grid)
        for site in sites:
            k = store.row_key(g.key, *(site if site is not None else ("", "")))
            if k not in done:
                tasks.append(dict(kind="threshold", geom=g, cfg=cfg, site=site))
    sim.log(f"{len(tasks)} thresholds to run ({len(done)} already in {path})")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"FAILED {t['geom'].key} {t['site']}: {r['error']}\n{r.get('tb', '')}")
            return
        row = store.threshold_row(t["geom"].key, t["site"], tag, r["res"], r["wall_s"])
        store.append_row(path, row)
        sim.log(f"  {t['geom'].key:<11} {'dark' if t['site'] is None else t['site']}: {row['i_det']:.3f} {row['flag']}"
                f"  ({row['n_runs']} runs, {row['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result)
    cmd_report(a)


def cmd_report(a):
    path = csv_path(a)
    rows = store.load_rows(path)
    L = [f"THRESHOLDS  {path}   (bias in I_dep of the narrowest free cross-section)", ""]
    geoms, maps, samples, vals, met, dark = [], {}, {}, {}, {}, {}
    for key in a.geoms.split(","):
        g = geo.get(key)
        mine = [r for r in rows if r["Geometry"] == g.key]
        if not mine:
            continue
        d = [r for r in mine if r["x_nm"] == ""]
        dark[g.key] = float(d[0]["i_det"]) if d else None
        L.append(f"{g.label}   ({g.key}; free width {g.free_width_min_nm():.0f} nm; dark switching "
                 f"{an.fmt(dark[g.key])})")
        L.append(f"  {'site [nm]':>14}  {'i_det':>6}  {'flag':<8} {'t_count':>8} {'T_e peak':>9} {'load max':>9}")
        for r in sorted((r for r in mine if r["x_nm"] != ""), key=lambda r: (float(r["y_nm"]), float(r["x_nm"]))):
            L.append(f"  {'(%g, %g)' % (float(r['x_nm']), float(r['y_nm'])):>14}  {float(r['i_det']):.3f}  {r['flag']:<8} "
                     f"{r['t_count']:>8} {r['te_peak']:>9} {r['load_max']:>9}")
        s = an.photon_rows(mine, g.key)
        if len(s) >= 3:
            samples[g.key] = s
            maps[g.key] = an.unit_map(g, s)
            vals[g.key] = an.pixel_values(maps[g.key][2])
            met[g.key] = m = an.metrics(vals[g.key])
            geoms.append(g)
            L.append(f"  hard-step efficiency over the unit cell ({len(s)} sites, interpolated): i_10 {m['i_10']:.3f}, "
                     f"i_50 {m['i_50']:.3f}, i_95 {m['i_95']:.3f}, i_max {m['i_max']:.3f}")
        L.append("")
    text = "\n".join(L)
    sim.log(text)
    with open(path.replace(".csv", ".txt"), "w") as fh:
        fh.write(text + "\n")
    if geoms:
        from snspd import plotting as pl
        pl.fig_hitmaps(geoms, maps, samples, f"Threshold bias for a count  ({cli.tag(a)})", path.replace(".csv", "_map.png"))
        pl.fig_efficiency(geoms, vals, met, dark, f"Share of the cell that counts  ({cli.tag(a)})",
                          path.replace(".csv", "_efficiency.png"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--geoms", default="Plain_Wire,Double_400")
    ap.add_argument("--grid", default="critical", choices=["critical", "full", "probes"])
    ap.add_argument("--no-dark", action="store_true", help="skip the dark run of each geometry")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="results/thresholds")
    cli.add_model_args(ap)
    a = ap.parse_args(argv)
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
