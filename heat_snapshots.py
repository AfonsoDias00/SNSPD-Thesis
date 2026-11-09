#!/usr/bin/env python3
"""
heat_snapshots.py -- snapshots.py with Joule heat: what a detection looks like in the electro-thermal TDGL.  (7 Nov 2026)

One photon per case, 0.01 I_dep above and below the heated threshold at that site (results/heat/check.csv, NbTiN
two-temperature model).  Preset nbtin, settled start, 2.0 xi (the mesh of the heated thresholds), frames every ~2 tau0.
The temperature field is not in pyTDGL's output file, so it is recorded during the run.

  python heat_snapshots.py --cases "Plain_Wire:-58,0:0;Double_400:-80,0:0;Double_400:-40,220:100" --workers 6

Writes docs/fig/heat_cascade.png (channel currents, voltage, peak temperature) and
docs/fig/heat_snapshots_<geometry>_<x>_<y>.png (|psi|^2 with the normal region theta >= 1 outlined).
"""
from __future__ import annotations

import argparse
import os
from dataclasses import replace

os.environ["SNSPD_MESH_DIR"] = ""
os.environ["SNSPD_NO_PROGRESS"] = "1"
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import h5py
import numpy as np

from snspd import geometry as geo
from snspd import hotspot as hsm
from snspd import sim, store
from snspd import thermal
from snspd.units import FILM

import snapshots as snap

OFFSET = 0.01
T0 = 20.0


def run_case(task):
    from validation.common import polyline
    g, (x, y), y_cut, I_uA, cfg = task["geom"], task["site"], task["y_cut"], task["I"], task["cfg"]
    film, params = cfg.pset.film, cfg.pset.params
    ctx = sim.context(g, cfg.mesh_xi, film, params)
    dev = ctx["device"]
    eps = hsm.build_eps(ctx["hs"], "kinetic", x, y, 1.0, cfg.tmax, params)
    init = sim.settled_state(dev, I_uA, cfg)
    meters = [snap.ChannelCurrent(dev, polyline([(a + 0.3, y_cut), (b - 0.3, y_cut)]))
              for a, b in snap.channels(g, y_cut)]
    h5 = sim.temp_h5(cfg)
    sol = None
    try:
        solver = thermal.thermal_solver_class()(dev, sim.solver_options(cfg, cfg.tmax, h5), thermal=thermal.NBTIN_2T,
                                                eps_photon=eps, terminal_currents=dict(source=I_uA, drain=-I_uA))
        if init is not None:
            solver.psi_init, solver.mu_init = np.array(init[0]), np.array(init[1])
        heat = solver.heat
        rec_t, rec_th, clock = [], [], [0.0]
        step = heat.step

        def recording_step(dt, mu):              # keep theta every ~2 tau0 (solver time = sum of the steps)
            step(dt, mu)
            clock[0] += dt
            if not rec_t or clock[0] - rec_t[-1] >= 2.0:
                rec_t.append(clock[0])
                rec_th.append(heat.theta.astype(np.float32).copy())
        heat.step = recording_step
        sol = solver.solve()
        trig = sim.detect(sol, T0)
        t_v, v = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
        th_t = np.array(rec_t)
        th_peak = np.array([float(a.max()) for a in rec_th])
        want = [0.0, 10.0, 25.0]
        if trig.t is not None:
            want += [trig.t - T0 + d for d in (-10.0, 0.0, 15.0, 40.0, 100.0)]
        else:
            want += [80.0, 160.0, 300.0]
        want = sorted({round(w, 1) for w in want if 0 <= w <= cfg.tmax - T0})
        times, cur, frames = [], [], {}
        with h5py.File(sol.path, "r") as f:
            for gk in sorted(f["data"].keys(), key=int):
                t = float(f["data"][gk].attrs["time"])
                psi = np.asarray(f["data"][gk]["psi"])
                times.append(t)
                cur.append([m(psi) for m in meters])
                for r in want:
                    d = abs(t - (T0 + r))
                    if r not in frames or d < frames[r][0]:
                        frames[r] = (d, t, np.abs(psi) ** 2)
        out_frames = []
        for r, (_, t, p2) in sorted(frames.items()):
            k = int(np.argmin(np.abs(th_t - t)))
            out_frames.append((r, t, p2, rec_th[k]))
        latched = trig.slips_last >= sim.SUSTAIN_SLIPS
        return dict(task=task, points=np.asarray(dev.points), triangles=np.asarray(dev.triangles), t=np.array(times),
                    I_ch=np.array(cur), channels=snap.channels(g, y_cut), t_v=t_v, v=v, trig=trig, latched=latched,
                    th_t=th_t, th_peak=th_peak, frames=out_frames)
    finally:
        path = getattr(sol, "path", None) or h5
        if path and os.path.exists(path):
            os.remove(path)


def plot_cascade(results, cases, path):
    import matplotlib.pyplot as plt
    from snspd import plotting as pl
    pl.apply_style()
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]
    fig, axes = plt.subplots(3, len(cases), figsize=(5.4 * len(cases), 9.0), sharex=True, layout="constrained",
                             gridspec_kw=dict(height_ratios=[2.0, 1.0, 1.0]))
    axes = np.atleast_2d(axes).reshape(3, len(cases))
    for col, (gk, site, ycut) in enumerate(cases):
        runs = sorted([r for r in results if (r["task"]["geom"].key, r["task"]["site"]) == (gk, site)],
                      key=lambda r: -r["task"]["I"])
        ax, axv, axt = axes[:, col]
        for r in runs:
            above = r["task"]["I"] > r["task"]["I_det"]
            ls, lw = ("-", 2.0) if above else ("--", 1.3)
            for k, nm in enumerate(snap.channel_names(r["channels"], gk)):
                ax.plot(r["t"], r["I_ch"][:, k], color=colors[k % 3], ls=ls, lw=lw, label=nm if above else None)
            state = "latched (detects)" if r["latched"] else ("pulse, heals" if r["trig"].t is not None else "nothing")
            axv.plot(r["t_v"], np.abs(r["v"]), color=pl.INK if above else pl.INK2, ls=ls, lw=1.1,
                     label=f"{r['task']['I_dep_frac']:.3f} I_dep: {state}")
            axt.plot(r["th_t"], r["th_peak"], color=pl.INK if above else pl.INK2, ls=ls, lw=1.2)
        for a in (ax, axv, axt):
            a.axvline(T0, color=pl.INK, lw=0.9, ls="--")
        axt.axhline(1.0, color="#eb6834", lw=0.8, ls=":")
        g = geo.get(gk)
        ax.set_title(f"{g.label}: photon at ({site[0]:g}, {site[1]:g}) nm\nheated threshold "
                     f"{runs[0]['task']['I_det_frac']:.3f} I_dep;  cut at y = {ycut:g} nm", fontsize=10)
        ax.legend(fontsize=8, loc="lower left", title="channel (solid: above, dashed: below)", title_fontsize=7.5,
                  frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID)
        axv.legend(fontsize=7.5, loc="upper left", frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID)
        axv.set_ylim(0, 3.0)
        axt.set_xlabel("time t [τ0]")
        axt.set_xlim(0, 400)
    axes[0, 0].set_ylabel("channel supercurrent [µA]")
    axes[1, 0].set_ylabel("|probe voltage| [V0]")
    axes[2, 0].set_ylabel("peak θ  (1 = T_c)")
    fig.suptitle("After one photon, with Joule heat (NbTiN two-temperature model, γ = 10, 2.0 ξ): 0.01 I_dep above "
                 "(solid) and below (dashed) the heated threshold", fontsize=11)
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print("wrote", path)


def plot_snapshots(r, path, half_y=260.0):
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from snspd import plotting as pl
    pl.apply_style()
    g, (x0, y0) = r["task"]["geom"], r["task"]["site"]
    P = r["points"] * 1e3
    tri = mtri.Triangulation(P[:, 0], P[:, 1], r["triangles"])
    frames = r["frames"]
    fig, axes = plt.subplots(1, len(frames), figsize=(1.55 * len(frames) + 1.6, 5.0), sharey=True, layout="constrained")
    im = None
    for ax, (rel, t, p2, th) in zip(axes, frames):
        im = ax.tripcolor(tri, p2, shading="gouraud", cmap="Blues_r", vmin=0, vmax=1)
        if (th > 0.5).any():
            ax.tricontour(tri, th, levels=[0.5, 1.0], colors=["#eb6834", "#c0392b"], linewidths=[0.8, 1.4],
                          linestyles=["--", "-"])
        pl.draw_holes(ax, g)
        ax.plot(x0, y0, "*", ms=9, mfc="white", mec=pl.INK, zorder=5)
        ax.set_xlim(-g.width_nm / 2, g.width_nm / 2)
        ax.set_ylim(y0 - half_y, y0 + half_y)
        ax.set_aspect("equal")
        ax.grid(False)
        ax.set_xticks([-100, 0, 100])
        ax.set_title(f"t0 + {max(t - T0, 0.0):.0f} τ0", fontsize=9.5)
    axes[0].set_ylabel("y [nm]")
    cb = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.01)
    cb.set_label("|ψ|²   (red: θ = 1, normal;  orange dashed: θ = 0.5)")
    state = "latched: detects" if r["latched"] else "heals"
    fig.suptitle(f"{g.label}, photon at ({x0:g}, {y0:g}) nm, {r['task']['I_dep_frac']:.3f} I_dep, with heat: {state}",
                 fontsize=11)
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True, help="GEOM:X,Y:YCUT;...")
    ap.add_argument("--heat-csv", default="results/heat/check_dt0.025.csv")
    ap.add_argument("--mesh-xi", type=float, default=2.0)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="docs/fig")
    a = ap.parse_args()
    cases = []
    for c in a.cases.split(";"):
        gk, site, ycut = c.split(":")
        cases.append((gk, tuple(float(v) for v in site.split(",")), float(ycut)))
    cfg = replace(sim.SimConfig(tmax=400.0, mesh_xi=a.mesh_xi, hotspot="kinetic", sustained=False,
                                out_dir="results/tmp_heat_snap", preset="nbtin"), save_every=10)
    rows = store.load_rows(a.heat_csv)
    tasks = []
    for gk, site, ycut in cases:
        g = geo.get(gk)
        i_dep = cfg.pset.film.I_dep_uA(g.free_width_min_nm())
        r = [r for r in rows if r["geom"] == gk and r["variant"] == "NbTiN 2T"
             and (float(r["x_nm"]), float(r["y_nm"])) == site][0]
        i_det = float(r["i_heat"])
        for s in (+1, -1):
            i = i_det + s * OFFSET
            tasks.append(dict(geom=g, site=site, y_cut=ycut, I=i * i_dep, I_det=i_det * i_dep, I_det_frac=i_det,
                              I_dep_frac=i, cfg=cfg))
    results = []
    sim.execute(tasks, a.workers, results.append, fn=run_case)
    os.makedirs(a.out, exist_ok=True)
    plot_cascade(results, cases, os.path.join(a.out, "heat_cascade.png"))
    for r in results:
        if r["task"]["I"] > r["task"]["I_det"]:
            x, y = r["task"]["site"]
            plot_snapshots(r, os.path.join(a.out, f"heat_snapshots_{r['task']['geom'].key}_{x:g}_{y:g}.png"))
    for r in results:
        t = r["task"]
        print(f"{t['geom'].key:<12} {t['site']} {t['I_dep_frac']:.3f} I_dep: latched {r['latched']}, trigger {r['trig'].t}")


if __name__ == "__main__":
    main()
