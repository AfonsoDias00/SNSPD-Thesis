"""Figure style and the standard figures: one event (traces and snapshots), threshold maps, efficiency curves."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.tri import Triangulation

from .analysis import eta_curve, pixel_values

INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SURFACE = "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"                      # categorical slots 1 and 2
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", ["#eaf2fc", "#9cc3ee", "#2a78d6", "#0e3a73"])
HEAT_CMAP = LinearSegmentedColormap.from_list("seq_orange", ["#fdf0e9", "#f6b596", "#eb6834", "#7a2b0b"])


def apply_style():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11,
        "legend.frameon": False, "lines.linewidth": 2.0,
    })


def draw_holes(ax, geom):
    for h in geom.holes:
        ax.fill(*h.outline(120).T, facecolor=SURFACE, edgecolor=INK, lw=0.9, zorder=3)


# ==============================================================================
# One event
# ==============================================================================
def fig_event(ev, geom, cfg, points_um, triangles, path, title=""):
    """Time traces of one event (left) and snapshots of |psi|^2 and of the electron temperature (right).  The strip
    is drawn as it is defined: current flows down the page; a ring marks the photon site."""
    apply_style()
    tr, snaps, film = ev["trace"], ev.get("snaps", []), cfg.film
    t, I, v = tr["t"], tr["I"], tr["v"]
    fig = plt.figure(figsize=(16.0, 8.2), layout="constrained")
    left, right = fig.subfigures(1, 2, width_ratios=[1.0, 1.25])
    ax = left.subplots(4, 1, sharex=True)

    ax[0].plot(t, I / I[0], color=BLUE)
    ax[0].set_ylabel("detector current\nI / I_bias")
    if cfg.circuit is not None:
        ax[0].axhline(1 - cfg.count_fraction, color=INK2, lw=1.0, ls="--")
        ax[0].annotate("count level", (t[0], 1 - cfg.count_fraction), textcoords="offset points", xytext=(4, 4),
                       ha="left", fontsize=8.5, color=INK2)
        ax[0].set_ylim(min(0.78, float((I / I[0]).min()) - 0.03), 1.02)
    ax[1].plot(t, v * film.V_unit_uV * 1e-3, color=BLUE)
    ax[1].set_ylabel("voltage across\nthe strip [mV]")
    ax[2].plot(t, tr["te"], color=BLUE, label="electrons")
    ax[2].plot(t, tr["tph"], color=ORANGE, label="phonons")
    ax[2].axhline(1.0, color=INK2, lw=1.0, ls="--")
    ax[2].set_ylabel("hottest site\nT / T_c")
    ax[2].legend(loc="upper right", ncols=2, fontsize=8.5)
    ax[3].plot(t, v * film.I_unit_uA / I, color=BLUE)
    ax[3].set_ylabel("resistance\n[squares]")
    ax[3].set_xlabel("time [tau0]   (photon at %g; 1 tau0 = %.3f ps; grey lines: the snapshots)"
                     % (cfg.photon.t_arrival, film.tau0_ps))
    for a in ax:
        for s in snaps:
            a.axvline(s["t"], color=GRID, lw=0.9, zorder=0)
    if title:
        fig.suptitle(title, fontsize=11, x=0.01, ha="left")

    if snaps:
        tri = Triangulation(points_um[:, 0] * 1e3, points_um[:, 1] * 1e3, triangles)
        te_max = max(1.5, float(np.ceil(max(s["te"].max() for s in snaps) * 4) / 4))
        groups = right.subfigures(1, 2)
        specs = (("n", "|psi|^2", SEQ_CMAP, 0.0, 1.0), ("te", "T_e / T_c   (line: T_c)", HEAT_CMAP, cfg.thermal.t0, te_max))
        for sub, (key, label, cmap, vmin, vmax) in zip(groups, specs):
            axs = sub.subplots(1, len(snaps), squeeze=False)[0]
            im = None
            for a, s in zip(axs, snaps):
                im = a.tripcolor(tri, s[key], cmap=cmap, vmin=vmin, vmax=vmax, shading="gouraud")
                if key == "te" and s["te"].max() >= 1.0:
                    a.tricontour(tri, s["te"], levels=[1.0], colors=INK, linewidths=0.7)
                if ev.get("site") is not None:
                    a.plot(ev["site"][0], ev["site"][1], "o", ms=6, mfc="none", mec=INK, mew=0.8)
                a.set_xlim(-geom.width_nm / 2, geom.width_nm / 2)
                a.set_ylim(-geom.length_nm / 2, geom.length_nm / 2)
                a.set_aspect("equal")
                a.set_anchor("N")
                a.grid(False)
                a.set_xticks([]), a.set_yticks([])
                for sp in a.spines.values():
                    sp.set_visible(False)
                a.set_title(f"{s['t']:.0f}", fontsize=9)
            axs[0].set_title(f"t = {snaps[0]['t']:.0f}", fontsize=9)
            sub.colorbar(im, ax=list(axs), location="bottom", shrink=0.8, pad=0.02, label=label)
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ==============================================================================
# Threshold maps and efficiency
# ==============================================================================
def fig_hitmaps(geoms, maps, samples, title, path):
    apply_style()
    fig, axes = plt.subplots(1, len(geoms), squeeze=False, figsize=(4.6 * len(geoms) + 1.2, 6.4), layout="constrained")
    axes = list(axes[0])
    allv = np.concatenate([pixel_values(maps[g.key][2]) for g in geoms])
    vmin, vmax = np.floor(allv.min() * 50) / 50, np.ceil(allv.max() * 50) / 50
    im = None
    for ax, g in zip(axes, geoms):
        X, Y, Z = maps[g.key]
        im = ax.pcolormesh(X, Y, Z, cmap=SEQ_CMAP, vmin=vmin, vmax=vmax, shading="auto")
        draw_holes(ax, g)
        for x, y, _ in samples[g.key]:
            ax.plot(x, y, "o", ms=4.5, mfc="white", mec=INK, mew=0.8, zorder=4)
        ax.set_xlim(-g.width_nm / 2, g.width_nm / 2)
        ax.set_ylim(-g.cell_half_y_nm, g.cell_half_y_nm)
        ax.set_aspect("equal")
        ax.set_title(g.label, fontsize=10)
        ax.set_xlabel("x [nm]")
        ax.grid(False)
    axes[0].set_ylabel("y [nm]  (one unit cell)")
    cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    cb.set_label("threshold bias  I / I_dep(free width)   (darker = needs more bias)")
    fig.suptitle(title)
    fig.savefig(path, dpi=170)
    plt.close(fig)


def fig_efficiency(geoms, vals, met, dark, title, path):
    """Hard-step efficiency against bias; the dashed vertical line of each colour is that device's dark switching."""
    apply_style()
    fig, ax = plt.subplots(figsize=(7.6, 4.8), layout="constrained")
    lo = np.floor(min(vals[g.key][0] for g in geoms) * 20) / 20 - 0.05
    ib = np.linspace(lo, 1.0, 800)
    for g in geoms:
        m = met[g.key]
        ax.plot(ib, eta_curve(vals[g.key], ib), color=g.color,
                label=f"{g.label}:  i_50 = {m['i_50']:.3f},  i_95 = {m['i_95']:.3f}")
        if dark.get(g.key):
            ax.axvline(dark[g.key], color=g.color, ls="--", lw=1.0)
    ax.set_xlabel("bias  I / I_dep(free width)   (dashed: dark switching)")
    ax.set_ylabel("share of the unit cell that counts [%]")
    ax.set_ylim(-3, 103)
    ax.legend(loc="upper left", fontsize=8.5)
    ax.set_title(title)
    fig.savefig(path, dpi=170)
    plt.close(fig)
