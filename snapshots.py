"""
What happens after a photon: channel supercurrents, probe voltage and |psi| snapshots for one photon per device, just
above and just below the photon threshold at that site (read from the device's hitmap CSV).

  python snapshots.py --cases "Plain_Wire_X:-58,0:0;Double_400:-80,0:0" --hitmap-dir results/slots2 \
                      --preset nbtin --criterion trigger --mesh-xi 2.0

A case is GEOMETRY:X,Y:YCUT (photon site and the y of the cut across the strip for the channel currents).
Writes <out>/cascade_<tag>.png and <out>/snapshots_<tag>_<geometry>.png, frames every 10 solver steps (~2 tau0).
Photon runs of a preset with a settled start (nbtin) start from sim.settled_state, as in hitmap.py.
"""
import argparse
import os
from dataclasses import replace

import h5py
import numpy as np

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")
os.environ.setdefault("SNSPD_NO_PROGRESS", "1")

from snspd import geometry as geo  # noqa: E402
from snspd import hotspot as hsm  # noqa: E402
from snspd import sim, store  # noqa: E402
from snspd.units import FILM  # noqa: E402

T0 = hsm.HOTSPOT.t0_tau0
CASES = []         # (geometry key, photon site [nm], y of the channel cut [nm]); set from --cases
OFFSET = 0.01      # bias = threshold +/- OFFSET [I_dep]
TMAX_PLOT = 400.0  # time axis of the cascade figure [tau0] (gamma 10 triggers come late, up to ~390 tau0)


def channels(geom, y_nm, step=0.5):
    """Segments of the cut y = y_nm that are film, left to right: [(x_start, x_end)] in nm.  A strip without holes is
    split into its left and right halves."""
    if not geom.holes:
        h = geom.width_nm / 2 - 0.25
        return [(-h, 0.0), (0.0, h)]
    xs = np.arange(-geom.width_nm / 2 + 0.25, geom.width_nm / 2, step)
    film = ~geom.inside_hole(xs, np.full_like(xs, y_nm))
    segs, start = [], None
    for x, f in zip(xs, film):
        if f and start is None:
            start = x
        if not f and start is not None:
            segs.append((start, x - step))
            start = None
    if start is not None:
        segs.append((start, xs[-1]))
    return segs


class ChannelCurrent:
    """Supercurrent [uA] through a straight cut, from psi alone: K_s = (K0/4) Im(psi* grad psi) with psi linear on each
    triangle (as validation.common.supercurrent_through), the geometry precomputed once so each frame is cheap.
    Positive along +y (the bias direction)."""

    def __init__(self, device, xy_um):
        from matplotlib.tri import Triangulation
        pts, tri = np.asarray(device.points) * 1e3, np.asarray(device.triangles)
        xy = np.asarray(xy_um, float) * 1e3
        mid, d = 0.5 * (xy[1:] + xy[:-1]), np.diff(xy, axis=0)
        k = Triangulation(pts[:, 0], pts[:, 1], tri).get_trifinder()(mid[:, 0], mid[:, 1])
        if (k < 0).any():
            raise ValueError("cut leaves the mesh")
        v = pts[tri[k]]
        T = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]], axis=1)
        self.Tinv = np.linalg.inv(T)
        self.lam = np.linalg.solve(np.transpose(T, (0, 2, 1)), (mid - v[:, 0])[..., None])[..., 0]
        self.idx = tri[k]
        self.nrm = np.stack([-d[:, 1], d[:, 0]], axis=1)          # left normal x length
        self.scale = -FILM.K0_uA_per_nm / 4 * FILM.xi_nm         # minus: the cut runs left to right

    def __call__(self, psi):
        z = psi[self.idx]
        dz = np.stack([z[:, 1] - z[:, 0], z[:, 2] - z[:, 0]], axis=1)
        grad = np.einsum("nij,nj->ni", self.Tinv, dz)
        zm = z[:, 0] + self.lam[:, 0] * dz[:, 0] + self.lam[:, 1] * dz[:, 1]
        K = self.scale * np.imag(np.conj(zm)[:, None] * grad)
        return float(np.sum(K * self.nrm))


def run_case(task):
    """One TDGL run; returns time series and the snapshot frames (all small)."""
    from validation.common import polyline
    geom, (x, y), y_cut, I_uA, cfg = task["geom"], task["site"], task["y_cut"], task["I"], task["cfg"]
    ctx = sim.context(geom, cfg.mesh_xi, cfg.pset.film, cfg.pset.params)
    dev = ctx["device"]
    eps = hsm.build_eps(ctx["hs"], cfg.hotspot, x, y, 1.0, cfg.tmax, cfg.pset.params)
    init = sim.settled_state(dev, I_uA, cfg) if cfg.pset.settled_start else None
    meters = [ChannelCurrent(dev, polyline([(a + 0.3, y_cut), (b - 0.3, y_cut)])) for a, b in channels(geom, y_cut)]
    out = dict(task=task, points=np.asarray(dev.points), triangles=np.asarray(dev.triangles), channels=[])
    with sim.solved(dev, I_uA, eps, cfg, init=init) as sol:
        trig = sim.detect(sol)                                     # as hitmap.py: the trigger, then (sustained) settle
        trig = sim.settle(dev, I_uA, sol, trig, cfg) if cfg.sustained else trig._replace(switched=bool(trig.branch))
        out.update(t_v=np.asarray(sol.dynamics.time), v=np.asarray(sol.dynamics.voltage()), trig=trig)
        t_tr = trig.t if trig.t is not None else None
        rel = [0.0, 5.0, 12.0, 25.0]                                   # build-up of the kinetic hotspot
        if t_tr is not None:
            rel += [t_tr - T0 + d for d in (0.0, 15.0, 40.0, 80.0)]     # the trigger and what follows
        else:
            rel += [60.0, 120.0, 250.0]
        want = []
        for r in sorted(round(r, 1) for r in rel if 0.0 <= r <= cfg.tmax - T0):
            if not want or r - want[-1] >= 4.0:                         # drop frames closer than 4 tau0
                want.append(r)
        times, cur, frames = [], [], {}
        with h5py.File(sol.path, "r") as f:
            for g in sorted(f["data"].keys(), key=int):
                t = float(f["data"][g].attrs["time"])
                psi = np.asarray(f["data"][g]["psi"])
                times.append(t)
                cur.append([m(psi) for m in meters])
                for r in want:                               # keep the frame closest to each wanted time
                    d = abs(t - (T0 + r))
                    if r not in frames or d < frames[r][0]:
                        frames[r] = (d, t, np.abs(psi) ** 2)
        out.update(t=np.array(times), I_ch=np.array(cur), channels=channels(geom, y_cut),
                   frames=[(r, t, p2, np.asarray(eps(dev.points, t=t))) for r, (_, t, p2) in sorted(frames.items())])
    return out


# ----------------------------------------------------------------------------------------------------------------------
def outcome(trig):
    """Short description of a run under the full sustained criterion (window + continuation)."""
    if trig.t is None:
        return "no trigger"
    if trig.switched:
        return "detects"
    if trig.slips_settled is not None:
        return "resistive to 400 τ0, heals in continuation"
    return "triggers, then heals"


def plot_cascade(results, hotspot, path):
    import matplotlib.pyplot as plt
    from snspd import plotting as pl
    pl.apply_style()
    colors = ["#2a78d6", "#eb6834", "#1baf7a"]                 # categorical slots 1-3: left, middle, right channel
    fig, axes = plt.subplots(2, len(CASES), figsize=(5.4 * len(CASES), 7.4), sharex=True, layout="constrained",
                             gridspec_kw=dict(height_ratios=[2.2, 1]))
    for col, (gkey, site, y_cut) in enumerate(CASES):
        runs = sorted([r for r in results if r["task"]["geom"].key == gkey], key=lambda r: -r["task"]["I"])
        ax, axv = axes[0, col], axes[1, col]
        for r in runs:
            above = r["task"]["I"] > r["task"]["I_det"]
            ls, lw = ("-", 2.0) if above else ("--", 1.3)
            names = channel_names(r["channels"], gkey)
            for k, nm in enumerate(names):
                ax.plot(r["t"], r["I_ch"][:, k], color=colors[k], ls=ls, lw=lw,
                        label=nm if above else None)
            axv.plot(r["t_v"], np.abs(r["v"]), color=pl.INK if above else pl.INK2, ls=ls, lw=1.2 if above else 1.0,
                     label=f"{r['task']['I_dep_frac']:.3f} I_dep: {outcome(r['trig'])}")
            if r["trig"].t is not None and above:
                for a in (ax, axv):
                    a.axvline(r["trig"].t, color=pl.INK2, lw=0.8, ls=":")
        for a in (ax, axv):
            a.axvline(T0, color=pl.INK, lw=0.9, ls="--")
        g = geo.get(gkey)
        I_det = runs[0]["task"]["I_det"]
        ax.set_title(f"{g.label}: photon at ({site[0]:g}, {site[1]:g}) nm\n"
                     f"threshold {I_det:.1f} µA = {I_det / runs[0]['task']['I_dep']:.3f} I_dep;  cut at y = {y_cut:g} nm",
                     fontsize=10)
        ax.text(T0 + 3, 0.98, "photon", fontsize=8, color=pl.INK2, va="top", transform=ax.get_xaxis_transform())
        ax.legend(fontsize=8, loc="lower left", title="channel (solid: above, dashed: below)", title_fontsize=7.5,
                  frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID, framealpha=0.95)
        ax.set_xlim(0, TMAX_PLOT)
        axv.set_ylim(0, 3.0)
        axv.legend(fontsize=7.5, loc="upper left", frameon=True, facecolor=pl.SURFACE, edgecolor=pl.GRID,
                   framealpha=0.95)
        axv.set_xlabel("time t [τ0]")
    axes[0, 0].set_ylabel("channel supercurrent [µA]")
    axes[1, 0].set_ylabel("|probe voltage| [V0]")
    fig.suptitle(f"After one photon ({hotspot}, f = 1, 400 τ0): bias 0.01 I_dep above (solid) and "
                 f"below (dashed) the threshold at that site; dotted = trigger; first {TMAX_PLOT:g} of 400 τ0", fontsize=11)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print("wrote", path)


def channel_names(segs, gkey):
    if not geo.get(gkey).holes:
        names = ["left half (hit)", "right half"]
    elif len(segs) == 2:
        names = ["left (hit)", "right"]
    elif len(segs) == 3:
        names = ["left (hit)", "middle", "right"]
    else:
        names = [f"channel {i + 1}" for i in range(len(segs))]
    return [f"{n}  {b - a:.0f} nm" for n, (a, b) in zip(names, segs)]


def plot_snapshots(r, hotspot, path, half_y=260.0):
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
    for ax, (rel, t, p2, eps) in zip(axes, frames):
        im = ax.tripcolor(tri, p2, shading="gouraud", cmap="Blues_r", vmin=0, vmax=1)
        if (eps < 0.9).any():
            ax.tricontour(tri, eps, levels=[0.5, 0.9], colors=["#eb6834", "#eb6834"], linewidths=[1.2, 0.7],
                          linestyles=["-", "--"])
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
    cb.set_label("|ψ|²   (orange: hotspot ε = 0.5 solid, 0.9 dashed)")
    trig = r["trig"]
    state = outcome(trig) if trig.t is None else f"triggers at t0 + {trig.t - T0:.0f} τ0; {outcome(trig)}"
    fig.suptitle(f"{g.label}, photon at ({x0:g}, {y0:g}) nm, {r['task']['I_dep_frac']:.3f} I_dep "
                 f"({hotspot} hotspot): {state}", fontsize=11)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True, help="GEOM:X,Y:YCUT;GEOM:X,Y:YCUT;...")
    ap.add_argument("--hitmap-dir", required=True)
    ap.add_argument("--hotspot", choices=list(hsm.MODELS), default="kinetic")
    ap.add_argument("--preset", choices=list(sim.PRESETS), default="nbtin")
    ap.add_argument("--criterion", choices=["sustained", "trigger"], default="trigger")
    ap.add_argument("--mesh-xi", type=float, default=1.4)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="docs/fig")
    a = ap.parse_args()
    for c in a.cases.split(";"):
        gk, site, ycut = c.split(":")
        CASES.append((gk, tuple(float(v) for v in site.split(",")), float(ycut)))
    scratch = os.path.join("results", "tmp_snapshots")
    cfg = replace(sim.SimConfig(tmax=400.0, mesh_xi=a.mesh_xi, hotspot=a.hotspot, sustained=a.criterion == "sustained",
                                out_dir=scratch, preset=a.preset), save_every=10)
    hit_cfg = replace(cfg, save_every=sim.SimConfig.save_every)
    tag = cfg.tag
    tasks = []
    for gkey, site, y_cut in CASES:
        g = geo.get(gkey)
        I_dep = cfg.pset.film.I_dep_uA(g.free_width_min_nm())
        row = [r for r in store.load_rows(os.path.join(a.hitmap_dir, f"hitmap_{hit_cfg.tag}.csv")) if r["Geometry"] == gkey
               and (float(r["x_nm"]), float(r["y_nm"])) == site and float(r["energy_f"]) == 1.0][0]
        I_det = float(row["I_det_uA"])
        for s in (+1, -1):
            I = I_det + s * OFFSET * I_dep
            tasks.append(dict(geom=g, site=site, y_cut=y_cut, I=I, I_det=I_det, I_dep=I_dep, I_dep_frac=I / I_dep,
                              cfg=cfg))
    results = []
    sim.execute(tasks, a.workers, results.append, fn=run_case)
    os.makedirs(a.out, exist_ok=True)
    plot_cascade(results, tag, os.path.join(a.out, f"cascade_{tag}.png"))
    for r in results:
        if r["task"]["I"] > r["task"]["I_det"]:
            plot_snapshots(r, tag, os.path.join(a.out, f"snapshots_{tag}_{r['task']['geom'].key}.png"))
    for r in results:
        t = r["task"]
        print(f"{t['geom'].key:<12} {t['I_dep_frac']:.3f} I_dep: trigger {r['trig']}")


if __name__ == "__main__":
    main()
