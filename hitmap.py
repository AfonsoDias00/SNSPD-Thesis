#!/usr/bin/env python3
"""
hitmap.py -- photon-threshold maps for internal-hole geometries vs the plain wire.

For every geometry: the dark switching current I_sw, the photon threshold i_det = I_det / I_sw on a grid of absorption
sites over one unit cell, thresholds for weaker hotspots at probe sites (selectivity), the steady-state current
crowding, and a go / no-go comparison of each hole geometry against the plain wire.  The model, units and detection
criterion are documented in snspd/ (units.py, hotspot.py, sim.py).

Usage
-----
  python hitmap.py selftest                     # hotspot-solver checks, legacy energy-loss table, hotspot figure
  python hitmap.py plan   [options]             # how many solver runs are still missing
  python hitmap.py run    [options]             # simulate (resumable) and report
  python hitmap.py report [options]             # rebuild figures and summary from the saved files

Options
  --geoms Plain_Wire,Long_Zipper   first = reference for the go / no-go comparison
  --grid {full,critical}           critical = ~10 sites per geometry, full = the reference grid
  --mesh-xi 1.4                    max mesh edge in xi (0.8 for the final numbers)
  --tmax 150                       solver window in tau0
  --hotspot {conserving,legacy}
  --energies 0.6,0.35              hotspot energy factors at the probe sites
  --workers 4                      parallel processes (~1.8 GB each at 0.8 xi)
  --rho-uohm-cm 150                film resistivity, only to quote times in ps

What the dark-count proxies mean: deterministic TDGL has no fluctuations, so I_sw is an instability limit, not a
dark-count rate.  Selectivity (weak seeds need more bias) and K_peak at the 95 % operating point are cheap proxies that
can falsify a geometry, not prove one.
"""
from __future__ import annotations

import argparse
import os
from dataclasses import replace

import numpy as np

from snspd import analysis as an
from snspd import geometry as geo
from snspd import hotspot as hsm
from snspd import plotting as pl
from snspd import sim, store
from snspd.units import Film

OUT_DIR = "results/hitmap"


def load_isw(cfg):
    """{geom key: I_sw [uA]} (the JSON also keeps flag and run count)."""
    return {k: v["I_sw"] for k, v in store.load_json(cfg.path("isw", "json"), {}).items()}


def plan_tasks(cfg, a, geoms, isw):
    """All hit and crowding tasks still missing, most informative first."""
    have = {store.row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]) for r in store.load_rows(cfg.path("hitmap", "csv"))}
    have_crowd = {(c["geom"], round(c["bias"], 2)) for c in store.load_json(cfg.path("crowd", "json"), [])}
    tasks = []
    base = dict(cfg=cfg, film=a.film)
    for g in geoms:
        if g.key not in isw:
            continue
        probes = set(g.probes)
        for x, y in g.hit_grid(a.grid):
            edge = x <= -(g.width_nm / 2 - 6)
            bridge = any(abs(x - h.cx_nm) <= 7 and abs(y - h.cy_nm) >= h.b_nm - 9 for h in g.holes)
            prio = 0 if ((x, y) in probes or edge or bridge) else 1
            if store.row_key(g.key, x, y, 1.0) not in have:
                tasks.append(dict(base, kind="hit", geom=g, x_nm=x, y_nm=y, f=1.0, I_sw=isw[g.key], prio=prio))
            if not a.no_probes and (x, y) in probes:
                for f in a.energies:
                    if store.row_key(g.key, x, y, f) not in have:
                        tasks.append(dict(base, kind="hit", geom=g, x_nm=x, y_nm=y, f=f, I_sw=isw[g.key], prio=1))
        if not a.no_crowding:
            for b in a.crowd_biases:
                if (g.key, round(b, 2)) not in have_crowd:
                    tasks.append(dict(base, kind="crowd", geom=g, bias=b, I_sw=isw[g.key], prio=0))
    tasks.sort(key=lambda t: t["prio"])
    return tasks


def cmd_plan(cfg, a, geoms):
    n_bis = int(np.ceil(np.log2((1 - cfg.hit_lo_frac) / sim.TOL_REL)))
    tasks = plan_tasks(cfg, a, geoms, {g.key: 1.0 for g in geoms})
    n_hit = sum(t["kind"] == "hit" for t in tasks)
    n_crowd = sum(t["kind"] == "crowd" for t in tasks)
    n_isw = sum(g.key not in load_isw(cfg) for g in geoms)
    sim.log(f"config {cfg.tag}, grid={a.grid}: {n_isw} I_sw + {n_hit} threshold searches + {n_crowd} crowding runs"
            f"  ~ {n_isw * 8 + n_hit * n_bis + n_crowd} solver runs of {cfg.tmax:g} tau0")


def cmd_run(cfg, a, geoms):
    sim.need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    sim.log(f"== run: {cfg.tag}  grid={a.grid}  workers={a.workers}  out={cfg.out_dir}")
    isw_path = cfg.path("isw", "json")
    isw_j = store.load_json(isw_path, {})
    todo = [dict(kind="isw", geom=g, cfg=cfg, film=a.film) for g in geoms if g.key not in isw_j]
    if todo:
        sim.log(f"-- dark switching current for {[t['geom'].key for t in todo]}")

        def on_isw(r):
            if not r["ok"]:
                sim.log(f"   FAILED {r['task']['geom'].key}: {r['error']}")
                return
            isw_j[r["task"]["geom"].key] = r["res"]
            store.save_json(isw_path, isw_j)
            sim.log(f"   {r['task']['geom'].key}: I_sw = {r['res']['I_sw']:.1f} uA  ({r['res']['n_sims']} runs, "
                    f"{r['wall_s']:.0f} s)" + (f"  [{r['res']['flag']}]" if r["res"]["flag"] else ""))

        sim.execute(todo, a.workers, on_isw)
    isw = load_isw(cfg)
    missing = [g.key for g in geoms if g.key not in isw]
    if missing:
        raise SystemExit(f"I_sw missing for {missing}; cannot continue.")

    tasks = plan_tasks(cfg, a, geoms, isw)
    sim.log(f"-- {len(tasks)} threshold / crowding tasks")
    crowd_path = cfg.path("crowd", "json")
    crowd = store.load_json(crowd_path, [])

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['kind']} {t['geom'].key}: {r['error']}")
            return
        if t["kind"] == "hit":
            row = store.hit_row(t["geom"].key, r["res"], t, cfg, r["wall_s"])
            store.append_row(cfg.path("hitmap", "csv"), row)
            sim.log(f"   {t['geom'].key:<11} ({t['x_nm']:6.1f},{t['y_nm']:5.1f}) f={t['f']:.2f}  i_det = {row['i_det']:.3f}"
                    + (f"  [{row['flag']}]" if row["flag"] else "") + f"   ({row['n_sims']} runs, {r['wall_s']:.0f} s)")
        else:
            res = dict(r["res"], geom=t["geom"].key)
            crowd.append(res)
            store.save_json(crowd_path, crowd)
            sim.log(f"   crowding {t['geom'].key:<11} bias {t['bias']:.2f}: K_peak/K_uniform = {res['Ksm'] / res['K_uniform']:.2f}")

    sim.execute(tasks, a.workers, on_result)
    cmd_report(cfg, a, geoms)


def enrich_crowd(cfg, crowd, film, geoms):
    """(Re)compute the interior disc-averaged peak Ksm from the saved K maps, so every entry uses the current
    definition (older runs stored a boundary-dominated value)."""
    by_key = {g.key: g for g in geoms}
    for c in crowd:
        f = cfg.path(f"crowd_{c['geom']}_{c['bias']:.2f}", "npz")
        if c["geom"] in by_key and os.path.exists(f):
            d_ = np.load(f)
            c["Ksm"], _ = sim.smoothed_peak(d_["points"], d_["K"], by_key[c["geom"]], 1.5 * film.xi_nm * 1e-3)
    return crowd


def cmd_report(cfg, a, geoms):
    pl.apply_style()
    rows = store.load_rows(cfg.path("hitmap", "csv"))
    isw = load_isw(cfg)
    crowd = enrich_crowd(cfg, store.load_json(cfg.path("crowd", "json"), []), a.film, geoms)
    geoms = [g for g in geoms if g.key in isw and an.photon_rows(rows, g.key)]
    if not geoms:
        sim.log("nothing to report yet (no photon-level rows).")
        return
    samples = {g.key: an.photon_rows(rows, g.key) for g in geoms}
    maps = {g.key: an.unit_map(g, samples[g.key]) for g in geoms}
    vals = {k: an.pixel_values(m[2]) for k, m in maps.items()}
    met = {k: an.metrics(v) for k, v in vals.items()}
    tabs = {g.key: an.probe_table(rows, g) for g in geoms}
    sel = {k: an.selectivity(t) for k, t in tabs.items()}
    label = f"[{cfg.hotspot} hotspot, T = {cfg.tmax:g} tau0, mesh {cfg.mesh_xi:g} xi]"
    pl.fig_hitmaps(geoms, maps, samples, isw, f"Photon-threshold maps {label}, f = 1", cfg.path("hitmap", "png"))
    pl.fig_efficiency(geoms, vals, met, isw, f"Hard-step efficiency from the maps {label}", cfg.path("efficiency", "png"))
    pl.fig_selectivity(geoms, tabs, cfg.path("selectivity", "png"))
    for b in a.crowd_biases[-2:]:
        pl.fig_crowding(geoms, {g.key: cfg.path(f"crowd_{g.key}_{b:.2f}", "npz") for g in geoms}, b,
                        cfg.path(f"crowding_{b:.2f}", "png"), a.film)

    L = []
    P = L.append
    P("=" * 78)
    P(f"SUMMARY   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau0   mesh={cfg.mesh_xi:g} xi   grid={a.grid}")
    if a.grid == "critical":
        P("NOTE: critical grid -> the maps are interpolated from few points; quantiles are indicative only.")
    for line in a.film.describe() + hsm.HOTSPOT.describe(a.film):
        P(line)
    P("=" * 78)
    P(f"{'':<30}" + "".join(f"{g.label:>14}" for g in geoms))
    P(f"{'I_sw [uA]':<30}" + "".join(f"{isw[g.key]:>14.1f}" for g in geoms))
    P(f"{'narrowest free width [nm]':<30}" + "".join(f"{g.free_width_min_nm():>14.0f}" for g in geoms))
    P(f"{'I_sw / I_dep(free width)':<30}" + "".join(f"{isw[g.key] / a.film.I_dep_uA(g.free_width_min_nm()):>14.3f}" for g in geoms))
    for key, name in [("i_10", "i_10"), ("i_50", "i_50"), ("i_90", "i_90"), ("i_95", "i_95"), ("i_99", "i_99"),
                      ("i_max", "i_max (100 % point)"), ("width_10_90", "10-90 % width")]:
        P(f"{name:<30}" + "".join(f"{met[g.key][key]:>14.3f}" for g in geoms))
    for key in ("i_10", "i_50", "i_95", "i_max"):
        P(f"{key + '  in units of I_dep':<30}"
          + "".join(f"{met[g.key][key] * isw[g.key] / a.film.I_dep_uA(g.free_width_min_nm()):>14.3f}" for g in geoms))
    P(f"{'  (samples used)':<30}" + "".join(f"{len(samples[g.key]):>14d}" for g in geoms))

    ref, cands = geoms[0], geoms[1:]
    for c in cands:
        pw, zp = met[ref.key], met[c.key]
        P(f"{'onset gain  ref - ' + c.key:<30}" + f"  i_10: {pw['i_10'] - zp['i_10']:+.3f}   i_50: {pw['i_50'] - zp['i_50']:+.3f}"
          "   (positive = holes turn on earlier)")
        P(f"{'plateau gain ref - ' + c.key:<30}" + f"  i_90: {pw['i_90'] - zp['i_90']:+.3f}   i_95: {pw['i_95'] - zp['i_95']:+.3f}"
          f"   i_99: {pw['i_99'] - zp['i_99']:+.3f}")
        vz = vals[c.key]
        P(f"{c.key + ' area above 0.87/0.90':<30}" + f"  {100 * np.mean(vz > 0.87):.1f} % / {100 * np.mean(vz > 0.90):.1f} %"
          "   (the tail that sets i_95 / i_99)")

    K95, Kfail, Kcrowd, Kraw = {}, {}, {}, {}
    for g in geoms:
        ent = [c for c in crowd if c["geom"] == g.key]
        if len(ent) >= 2:
            K95[g.key] = an.K_of_I(ent, met[g.key]["i_95"] * isw[g.key])
            Kfail[g.key] = an.K_of_I(ent, isw[g.key])
            top = max(ent, key=lambda c: c["bias"])
            Kcrowd[g.key] = top["Ksm"] / top["K_uniform"]
            Kraw[g.key] = top["Kmax"] / top["K_uniform"]
    if K95:
        P("-" * 78)
        P("current crowding (steady state; K_peak = max of K averaged over discs of 1.5 xi, >= 1.5 xi from all walls)")
        P(f"{'K_peak / K_uniform @ top bias':<30}" + "".join(f"{an.fmt(Kcrowd.get(g.key), 2):>14}" for g in geoms))
        P(f"{'  (single-node max, noisy)':<30}" + "".join(f"{an.fmt(Kraw.get(g.key), 2):>14}" for g in geoms))
        P(f"{'K_peak at own i_95 [uA/um]':<30}" + "".join(f"{an.fmt(K95.get(g.key), 0):>14}" for g in geoms))
        P(f"{'K_peak at own I_sw [uA/um]':<30}" + "".join(f"{an.fmt(Kfail.get(g.key), 0):>14}" for g in geoms))
        for g in geoms:
            ent = sorted([c for c in crowd if c["geom"] == g.key], key=lambda c: c["bias"])
            if len(ent) >= 2:
                lin = [c["Ksm"] / c["I"] for c in ent]
                P(f"  linearity {g.label}: K/I changes by {100 * (lin[-1] / lin[0] - 1):+.1f} % from bias "
                  f"{ent[0]['bias']:.2f} to {ent[-1]['bias']:.2f}")
    if any(sel.values()):
        P("-" * 78)
        P("selectivity  median[i_det(f) - i_det(1)] over probe sites (larger = weak seeds need more bias)")
        for f in sorted({f for s in sel.values() for f in s}):
            P(f"  f = {f:<5g}".ljust(30) + "".join(f"{an.fmt(sel[g.key].get(f)):>14}" for g in geoms))
        for g in geoms:
            nw = [(s, f) for s, d_ in tabs[g.key].items() for f, (v, fl) in d_.items() if fl == "at_wall"]
            if nw:
                P(f"  {g.label}: no detection below I_sw at {len(nw)} (site, f) pairs -> counted as i_det = 1")

    P("-" * 78)
    P(f"GO / NO-GO vs {ref.label}   (margins: --go-margin {a.go_margin:g}, --sel-margin {a.sel_margin:g})")
    for c in cands:
        verdict = []
        dz = met[ref.key]["i_95"] - met[c.key]["i_95"]
        verdict.append(dz >= a.go_margin)
        P(f"  {c.label}  C1 detection : ref i_95 - i_95 = {dz:+.3f}   need >= {a.go_margin:.3f}   -> {'PASS' if verdict[-1] else 'FAIL'}")
        if ref.key in K95 and c.key in K95:
            verdict.append(K95[c.key] <= K95[ref.key])
            P(f"  {c.label}  C2 stress    : K_peak@i_95 {K95[c.key]:.0f} vs ref {K95[ref.key]:.0f} uA/um   need <= ref"
              f"   -> {'PASS' if verdict[-1] else 'FAIL'}")
        common = sorted(set(sel.get(ref.key, {})) & set(sel.get(c.key, {})))
        if common:
            f = common[0]
            ds = sel[c.key][f] - sel[ref.key][f]
            verdict.append(ds >= a.sel_margin)
            P(f"  {c.label}  C3 selectivity: - ref at f = {f:g}: {ds:+.3f}   need >= {a.sel_margin:.3f}"
              f"   -> {'PASS' if verdict[-1] else 'FAIL'}")
        P(f"  {c.label}  OVERALL: {'GO' if all(verdict) else 'NO-GO / not shown'}   [{sum(verdict)}/{len(verdict)} criteria passed]")
    P("-" * 78)
    P("caveats: deterministic TDGL; one photon energy (plus probes); sparse sampling near the edge is extrapolated;"
      " K_peak at hole tips depends on the mesh.")
    text = "\n".join(L)
    sim.log(text)
    with open(cfg.path("summary", "txt"), "w") as fh:
        fh.write(text + "\n")
    store.save_json(cfg.path("summary", "json"), dict(
        metrics=met, selectivity={k: {str(f): v for f, v in s.items()} for k, s in sel.items()},
        K95=K95, Kfail=Kfail, I_sw=isw, free_width_nm={g.key: g.free_width_min_nm() for g in geoms}))
    sim.log(f"figures and summary written to {cfg.out_dir}/ (tag {cfg.tag})")


def cmd_selftest(cfg, a, geoms):
    import time
    pl.apply_style()
    os.makedirs(cfg.out_dir, exist_ok=True)
    sim.log("== 1. Energy kept by the legacy (free-space) hotspot inside the film and outside the holes")
    dts = (0.0, 10.0, 30.0, 60.0, 100.0)
    sim.log(f"   {'site':<26}" + "".join(f"  t0+{int(t):<4d}" for t in dts))
    for g in geoms:
        for x, y in g.probes:
            _, ret = hsm.legacy_energy_retention(g, x, y, dts)
            sim.log(f"   {g.key:<12} ({x:4.0f},{y:3.0f}) nm  " + "".join(f"  {r:6.2f} " for r in ret))
    sim.log("   (the conserving model keeps 1.00 everywhere by construction)")
    sim.need_tdgl()
    sim.log("\n== 2. Hotspot solver on the real meshes")
    p = hsm.HOTSPOT
    for g in geoms:
        t0 = time.time()
        hs = sim.context(g, cfg.mesh_xi, a.film)["hs"]
        snaps = hs.evolve(hs.initial_condition(0.0, 0.0), 20.0)
        mass = snaps.astype(float) @ hs.m
        drift = float(abs(mass.max() - mass.min()) / mass[0])
        lines = [f"  [{'PASS' if drift < 1e-5 else 'FAIL'}] mass conserved: max drift {drift:.1e}"]
        if not g.holes:
            P = hs.P
            r2 = (snaps.astype(float) * (hs.m * (P[:, 0] ** 2 + P[:, 1] ** 2))).sum(1) / mass
            slope = np.polyfit(np.arange(len(snaps)) * hs.dt, r2, 1)[0] / (4 * hs.D)
            lines.append(f"  [{'PASS' if abs(slope - 1) < 0.08 else 'FAIL'}] free-space <r^2> growth / 4D: {slope:.3f}")
            k = int(round(10.0 / hs.dt))
            r0 = p.r0_nm * 1e-3
            rt2 = r0 ** 2 + 4 * hs.D * 10.0
            ana = (r0 ** 2 / rt2) * np.exp(-(P[:, 0] ** 2 + P[:, 1] ** 2) / rt2)
            s = np.abs(P[:, 1]) < 0.15
            err = float(np.linalg.norm(snaps[k][s] - ana[s]) / np.linalg.norm(ana[s]))
            lines.append(f"  [{'PASS' if err < 0.08 else 'FAIL'}] matches the legacy profile in free space (dt = 10): rel L2 err {err:.3f}")
        sim.log(f"{g.key}: {len(hs.P)} sites, {time.time() - t0:.1f} s\n" + "\n".join(lines))
    sim.log("\n== 3. Figure: legacy vs conserving epsilon for an edge hit and a bridge hit")
    import matplotlib.pyplot as plt
    cases = [(g, *g.probes[0]) for g in geoms[:1]] + [(g, *g.probes[-1]) for g in geoms[1:2]]
    fig, axes = plt.subplots(len(cases), 2, figsize=(9.5, 4.25 * len(cases)), squeeze=False)
    im = None
    for r_i, (g, x, y) in enumerate(cases):
        hs = sim.context(g, cfg.mesh_xi, a.film)["hs"]
        for c_i, mode in enumerate(("legacy", "conserving")):
            e = hsm.build_eps(hs, mode, x, y, 1.0, cfg.tmax)(hs.P, t=p.t0_tau0 + 10.0)
            ax = axes[r_i, c_i]
            im = ax.tripcolor(hs.P[:, 0] * 1e3, hs.P[:, 1] * 1e3, hs.T, 1.0 - e, shading="gouraud", cmap=pl.SEQ_CMAP, vmin=0, vmax=2.6)
            pl.draw_holes(ax, g)
            ax.plot(x, y, "*", ms=11, mfc="white", mec=pl.INK, zorder=5)
            ax.set_xlim(-g.width_nm / 2, g.width_nm / 2)
            ax.set_ylim(y - 90, y + 90)
            ax.set_aspect("equal")
            ax.set_title(f"{g.label} hit ({x:.0f},{y:.0f})  -  {mode}", fontsize=10)
            ax.grid(False)
    fig.colorbar(im, ax=axes, shrink=0.8, label="1 - epsilon  at t0 + 10 tau0")
    path = cfg.path("selftest_hotspot", "png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    sim.log(f"   wrote {path}")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "plan", "run", "report"])
    ap.add_argument("--geoms", default="Plain_Wire,Long_Zipper", help=f"comma separated, from {list(geo.GEOMETRIES)}")
    ap.add_argument("--tmax", type=float, default=sim.SimConfig.tmax)
    ap.add_argument("--hotspot", choices=["conserving", "legacy"], default=sim.SimConfig.hotspot)
    ap.add_argument("--mesh-xi", type=float, default=sim.SimConfig.mesh_xi)
    ap.add_argument("--save-every", type=int, default=sim.SimConfig.save_every)
    ap.add_argument("--hit-lo-frac", type=float, default=sim.SimConfig.hit_lo_frac)
    ap.add_argument("--grid", choices=["full", "critical"], default="full")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--energies", default="0.6,0.35", help="probe-site energy factors, comma separated")
    ap.add_argument("--crowd-biases", default="0.80,0.90,0.95")
    ap.add_argument("--no-crowding", action="store_true")
    ap.add_argument("--no-probes", action="store_true")
    ap.add_argument("--go-margin", type=float, default=0.03, help="holes' i_95 must beat the reference by this much")
    ap.add_argument("--sel-margin", type=float, default=0.01, help="holes' selectivity must beat the reference by this much")
    ap.add_argument("--rho-uohm-cm", type=float, default=None, help="film resistivity; only used to quote times in ps")
    ap.add_argument("--out", default=OUT_DIR)
    a = ap.parse_args(argv)
    a.energies = tuple(float(x) for x in a.energies.split(",") if x)
    a.crowd_biases = tuple(float(x) for x in a.crowd_biases.split(",") if x)
    a.film = Film(rho_uohm_cm=a.rho_uohm_cm)
    cfg = sim.SimConfig(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi, save_every=a.save_every,
                        hit_lo_frac=a.hit_lo_frac, out_dir=a.out)
    geoms = [geo.get(k) for k in a.geoms.split(",") if k]
    return a, cfg, geoms


def main(argv=None):
    a, cfg, geoms = parse_args(argv)
    if cfg.tmax <= hsm.HOTSPOT.t0_tau0 + 10:
        raise SystemExit("--tmax must exceed the photon time by a safe margin")
    {"selftest": cmd_selftest, "plan": cmd_plan, "run": cmd_run, "report": cmd_report}[a.command](cfg, a, geoms)


if __name__ == "__main__":
    main()
