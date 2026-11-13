#!/usr/bin/env python3
"""
event.py -- one photon, start to finish: what the model does at one site and one bias.

Each case is "<geometry>:<x>,<y>:<bias>" with the bias in I_dep of the free width; "<geometry>:dark:<bias>" has no
photon.  The whole window is recorded (the run does not stop at the count).

  python event.py --cases "Double_400:-80,0:0.80;Plain_Wire:-58,0:0.90" --mesh-xi 2.0 --tmax 600 --workers 2
  python event.py --cases "Double_400:-80,0:0.80" --source              current source instead of the circuit
  python event.py --cases "Double_400:-80,0:0.80" --thermal latent=False --label nolatent

Writes, per case, results/events/<case>_<tag>.json (the outcome), .npz (traces, snapshots, mesh) and .png.
"""
from __future__ import annotations

import argparse
import os
from dataclasses import replace

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")         # before numpy is imported: one BLAS thread per worker process

import numpy as np

from snspd import cli, sim, store
from snspd import geometry as geo


def parse_cases(text):
    out = []
    for item in filter(None, (s.strip() for s in text.split(";"))):
        g, site, bias = item.split(":")
        out.append((g, None if site == "dark" else tuple(float(v) for v in site.split(",")), float(bias)))
    return out


def stem(g, site, bias, tag):
    where = "dark" if site is None else f"{site[0]:g}_{site[1]:g}"
    return f"{g}_{where}_i{bias:.3f}_{tag}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--snaps", type=int, default=6, help="number of snapshots, spread over the window after the photon")
    ap.add_argument("--snap-times", default="", help="explicit snapshot times [tau0], comma separated")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default="results/events")
    ap.add_argument("--replot", action="store_true", help="redraw the figures of these cases from their saved files")
    cli.add_model_args(ap)
    a = ap.parse_args(argv)
    cfg = replace(cli.config(a, stop_on_count=False), tmax_ext=a.tmax)   # a fixed window: this script records only
    tag = cli.tag(a)
    for line in cfg.describe():
        sim.log(line)
    ta = cfg.photon.t_arrival
    times = ([float(v) for v in a.snap_times.split(",")] if a.snap_times
             else [round(ta + (cfg.tmax - ta) * f, 1) for f in np.linspace(0.04, 1.0, a.snaps) ** 1.6])
    os.makedirs(a.out, exist_ok=True)
    if a.replot:
        import json
        from snspd import plotting as pl
        for g, site, bias in parse_cases(a.cases):
            name = os.path.join(a.out, stem(g, site, bias, tag))
            d, ev = np.load(name + ".npz"), json.load(open(name + ".json"))
            ev["trace"] = {k[6:]: d[k] for k in d.files if k.startswith("trace_")}
            ev["snaps"] = [dict(t=float(t), n=n, te=te) for t, n, te in zip(d["snap_t"], d["snap_n"], d["snap_te"])]
            pl.fig_event(ev, geo.get(g), cfg, d["points"], d["triangles"], name + ".png", title=ev.get("title", ""))
            sim.log(name + ".png")
        return
    tasks = [dict(kind="event", geom=geo.get(g), cfg=cfg, site=site, i_bias=bias, trace=True, snap_times=tuple(times))
             for g, site, bias in parse_cases(a.cases)]

    def on_result(r):
        t = r["task"]
        g, site, bias = t["geom"], t["site"], t["i_bias"]
        name = stem(g.key, site, bias, tag)
        if not r["ok"]:
            sim.log(f"FAILED {name}: {r['error']}\n{r.get('tb', '')}")
            return
        ev = r["res"]
        ev["site"] = site
        if ev["dark"]:
            sim.log(f"{name}: no steady dark state at this bias (dark switching)")
            return
        tr, snaps = ev["trace"], ev.get("snaps", [])
        what = "no photon" if site is None else f"photon at ({site[0]:g}, {site[1]:g}) nm"
        outcome = (f"count at t = {ev['t_count']:.0f}" if ev["t_count"] is not None else "no count") if cfg.circuit \
            else ("resistive at the end" if ev["count"] else "superconducting at the end")
        title = f"{g.label}, {what}, bias {bias:.3f} I_dep: {outcome}"
        summary = {k: v for k, v in ev.items() if k not in ("trace", "snaps", "final", "mesh")}
        store.save_json(os.path.join(a.out, name + ".json"),
                        dict(summary, geometry=g.key, tag=tag, title=title, config=cfg.describe()))
        pts, tri = ev["mesh"]["points"], ev["mesh"]["triangles"]
        np.savez_compressed(os.path.join(a.out, name + ".npz"), points=pts, triangles=tri,
                            snap_t=np.array([s["t"] for s in snaps]), snap_n=np.array([s["n"] for s in snaps]),
                            snap_te=np.array([s["te"] for s in snaps]), snap_tph=np.array([s["tph"] for s in snaps]),
                            **{f"trace_{k}": v for k, v in tr.items()})
        from snspd import plotting as pl
        pl.fig_event(ev, g, cfg, pts, tri, os.path.join(a.out, name + ".png"), title=title)
        au = ev["audit"]
        sim.log(f"{name}: {outcome}; load max {ev['load_max'] if ev['load_max'] is None else round(ev['load_max'], 3)}, "
                f"I_end/I_b {ev['I_end']:.3f}, R_end {ev['R_end']:.2f} sq, slips {ev['slips']:.1f}, T_e peak "
                f"{ev['te_peak']:.2f}, T_ph peak {ev['tph_peak']:.2f} T_c, normal area max {ev['normal_nm2']:.0f} nm^2, "
                f"heat [h nu]: Joule {au['joule'] / au['photon'] if au['photon'] else float('nan'):.1f}, relaxation "
                f"{au['relax'] / au['photon'] if au['photon'] else float('nan'):.1f}; {r['wall_s']:.0f} s")

    sim.execute(tasks, a.workers, on_result)


if __name__ == "__main__":
    main()
