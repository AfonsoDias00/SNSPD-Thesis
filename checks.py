#!/usr/bin/env python3
"""
checks.py -- photon thresholds at extra sites that test assumptions behind the hitmap.py maps.  The results go to
results/checks/ and never into the map data.

  symmetry   hitmap.py simulates only x <= 0 and fills the rest by symmetry: mirror x -> -x for the plain wire, the
             y-mirror and the glide (x, y) -> (-x, 61 - y) for the zipper.  The hole array is finite, so the device
             is not exactly symmetric.  Re-simulate image sites and compare with the stored originals.
  edge       the outer 4 nm of the strip are filled with the x = -96 nm value.  Simulate x = -98 nm.
  delay      single runs at fixed bias with a long window: when (and by which branch) does a photon at one site
             trigger?  A sharp threshold shows as a delay that diverges just above it; a slow, window-limited
             "threshold" shows as triggers spread over the whole window.

Usage
-----
  python checks.py run    symmetry --mesh-xi 0.8 --workers 4
  python checks.py run    edge     --mesh-xi 0.8
  python checks.py report symmetry --mesh-xi 0.8
  python checks.py delay  --geom Long_Zipper --site 0,0 --biases 0.76,0.78,0.80,0.82,0.84 --tmax 500

The I_sw and the reference thresholds come from results/hitmap/ with the same hotspot, window and mesh.
"""
from __future__ import annotations

import argparse
import os

from snspd import geometry as geo
from snspd import sim, store
from snspd.units import Film

CHECKS = {   # (geometry, test site, stored reference site) [nm]
    "symmetry": [
        ("Plain_Wire", (14.0, 0.0), (-14.0, 0.0)),
        ("Plain_Wire", (58.0, 0.0), (-58.0, 0.0)),
        ("Plain_Wire", (96.0, 0.0), (-96.0, 0.0)),
        ("Long_Zipper", (68.0, 61.0), (-68.0, 0.0)),       # glide image, right channel
        ("Long_Zipper", (96.0, 61.0), (-96.0, 0.0)),       # glide image, right edge
        ("Long_Zipper", (18.0, 31.0), (-18.0, 30.0)),      # glide image, centre channel
        ("Long_Zipper", (36.0, 0.0), (-36.0, 61.0)),       # glide image, right-column bridge
        ("Long_Zipper", (-50.0, -30.0), (-50.0, 30.0)),    # y-mirror
        ("Long_Zipper", (-36.0, -61.0), (-36.0, 61.0)),    # y-mirror, left-column bridge below
    ],
    "edge": [
        ("Plain_Wire", (-98.0, 0.0), (-96.0, 0.0)),
        ("Long_Zipper", (-98.0, 0.0), (-96.0, 0.0)),
        ("Long_Zipper", (-98.0, 61.0), (-96.0, 61.0)),
    ],
}
HITMAP_DIR = "results/hitmap"
OUT_DIR = "results/checks"


def stored(cfg):
    hcfg = sim.with_out_dir(cfg, HITMAP_DIR)
    isw = {k: v["I_sw"] for k, v in store.load_json(hcfg.path("isw", "json"), {}).items()}
    ref = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]): r
           for r in store.load_rows(hcfg.path("hitmap", "csv"))}
    return isw, ref


def cmd_run(cfg, a):
    sim.need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    isw, _ = stored(cfg)
    path = cfg.path(a.check, "csv")
    have = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]) for r in store.load_rows(path)}
    tasks = []
    for key, (x, y), _ in CHECKS[a.check]:
        if key not in isw:
            raise SystemExit(f"no I_sw for {key} at {cfg.tag} in {HITMAP_DIR}/; run hitmap.py first")
        if store.row_key(key, x, y, 1.0) not in have:
            tasks.append(dict(kind="hit", geom=geo.get(key), cfg=cfg, film=a.film, x_nm=x, y_nm=y, f=1.0, I_sw=isw[key]))
    sim.log(f"-- {a.check}: {len(tasks)} threshold searches at {cfg.tag}")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['geom'].key} ({t['x_nm']:g},{t['y_nm']:g}): {r['error']}")
            return
        row = store.hit_row(t["geom"].key, r["res"], t, cfg, r["wall_s"])
        store.append_row(path, row)
        sim.log(f"   {t['geom'].key:<11} ({t['x_nm']:6.1f},{t['y_nm']:6.1f})  i_det = {row['i_det']:.4f}  "
                f"[{row['branch']} at t = {row['t_trigger']}]  ({row['n_sims']} runs, {r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result)
    cmd_report(cfg, a)


def cmd_report(cfg, a):
    _, ref = stored(cfg)
    rows = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]): r for r in store.load_rows(cfg.path(a.check, "csv"))}
    L = [f"CHECK {a.check}   {cfg.tag}   (tol = bisection tolerance in units of I_sw)",
         f"{'geometry':<12}{'test site':>16}{'i_det':>9}{'branch':>10}{'t_trig':>8}{'ref site':>16}{'ref i_det':>11}{'diff':>9}{'tol':>8}"]
    for key, (x, y), (xr, yr) in CHECKS[a.check]:
        r, rr = rows.get(store.row_key(key, x, y, 1.0)), ref.get(store.row_key(key, xr, yr, 1.0))
        if r is None or rr is None:
            L.append(f"{key:<12}{f'({x:g},{y:g})':>16}   {'missing' if r is None else 'no reference'}")
            continue
        g = geo.get(key)
        tol = sim.scales(g, a.film)[0] / float(r["I_sw_uA"])
        d = float(r["i_det"]) - float(rr["i_det"])
        L.append(f"{key:<12}{f'({x:g},{y:g})':>16}{float(r['i_det']):>9.4f}{r['branch']:>10}{r['t_trigger']:>8}"
                 f"{f'({xr:g},{yr:g})':>16}{float(rr['i_det']):>11.4f}{d:>+9.4f}{tol:>8.4f}"
                 + ("" if abs(d) <= 1.01 * tol else "   <-- outside tolerance"))
    text = "\n".join(L)
    sim.log(text)
    with open(cfg.path(f"{a.check}_summary", "txt"), "w") as fh:
        fh.write(text + "\n")


def delay_task(task):
    """One fixed-bias run with a photon (worker process)."""
    import time
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        g, cfg = task["geom"], task["cfg"]
        ctx = sim.context(g, cfg.mesh_xi)
        eps = hsm.build_eps(ctx["hs"], cfg.hotspot, task["x_nm"], task["y_nm"], 1.0, cfg.tmax)
        trig = sim.run_simulation(ctx["device"], task["I"], eps, cfg)
        return dict(ok=True, task=task, res=trig._asdict(), wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def cmd_delay(cfg, a):
    """Trigger time vs bias at one site; bias is in units of the stored I_sw at --isw-tmax (same mesh)."""
    sim.need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    g = geo.get(a.geom)
    ref_cfg = sim.with_out_dir(sim.SimConfig(tmax=a.isw_tmax, hotspot=cfg.hotspot, mesh_xi=cfg.mesh_xi), HITMAP_DIR)
    i_sw = store.load_json(ref_cfg.path("isw", "json"), {})[g.key]["I_sw"]
    x, y = (float(v) for v in a.site.split(","))
    path = cfg.path(f"delay_{g.key}_{x:g}_{y:g}", "csv")
    fields = ["Geometry", "x_nm", "y_nm", "bias", "I_uA", "I_sw_ref_uA", "tmax", "switched", "branch", "t_trigger", "wall_s"]
    have = {round(float(r["bias"]), 4) for r in store.load_rows(path)}
    tasks = [dict(geom=g, cfg=cfg, x_nm=x, y_nm=y, bias=b, I=b * i_sw) for b in a.biases if round(b, 4) not in have]

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED bias {t['bias']}: {r['error']}")
            return
        res = r["res"]
        store.append_row(path, dict(Geometry=g.key, x_nm=x, y_nm=y, bias=t["bias"], I_uA=t["I"], I_sw_ref_uA=i_sw,
                                    tmax=cfg.tmax, switched=res["switched"], branch=res["branch"],
                                    t_trigger="" if res["t"] is None else round(res["t"], 2), wall_s=round(r["wall_s"], 1)), fields)

    sim.execute(tasks, a.workers, on_result, fn=delay_task)
    rows = sorted(store.load_rows(path), key=lambda r: float(r["bias"]))
    sim.log(f"DELAY {g.key} ({x:g},{y:g})  {cfg.tag}  bias in units of I_sw(T = {a.isw_tmax:g}) = {i_sw:.2f} uA")
    for r in rows:
        sim.log(f"  i = {float(r['bias']):.3f}   {'switched' if r['switched'] == 'True' else 'no switch':<9}  "
                f"{r['branch']:<9} t = {r['t_trigger'] or '-'}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report", "delay"])
    ap.add_argument("check", nargs="?", choices=list(CHECKS))
    ap.add_argument("--geom", default="Long_Zipper", help="delay: geometry key")
    ap.add_argument("--site", default="0,0", help="delay: photon site x,y [nm]")
    ap.add_argument("--biases", default="0.76,0.78,0.80,0.82,0.84", help="delay: biases in units of I_sw")
    ap.add_argument("--isw-tmax", type=float, default=150.0, help="delay: window of the stored I_sw used as unit")
    ap.add_argument("--hotspot", choices=["conserving", "legacy"], default=sim.SimConfig.hotspot)
    ap.add_argument("--tmax", type=float, default=sim.SimConfig.tmax)
    ap.add_argument("--mesh-xi", type=float, default=sim.SimConfig.mesh_xi)
    ap.add_argument("--hit-lo-frac", type=float, default=sim.SimConfig.hit_lo_frac)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default=OUT_DIR)
    a = ap.parse_args(argv)
    a.film = Film()
    a.biases = [float(s) for s in a.biases.split(",") if s]
    if a.command != "delay" and a.check is None:
        ap.error("run / report need a check name")
    cfg = sim.SimConfig(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi, hit_lo_frac=a.hit_lo_frac, out_dir=a.out)
    {"run": cmd_run, "report": cmd_report, "delay": cmd_delay}[a.command](cfg, a)


if __name__ == "__main__":
    main()
