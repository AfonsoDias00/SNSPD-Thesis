"""Figure style and the standard figures (threshold maps, efficiency curves, selectivity, current crowding)."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from .analysis import eta_curve, pixel_values

INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SURFACE = "#fcfcfb"
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", ["#eaf2fc", "#9cc3ee", "#2a78d6", "#0e3a73"])


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


def _panels(n, **kw):
    fig, axes = plt.subplots(1, n, squeeze=False, **kw)
    return fig, list(axes[0])


def fig_hitmaps(geoms, maps, samples, isw, title, path):
    fig, axes = _panels(len(geoms), figsize=(6.25 * len(geoms), 5.2), sharey=True, layout="constrained")
    allv = np.concatenate([pixel_values(maps[g.key][2]) for g in geoms])
    vmin, vmax = np.floor(allv.min() * 200) / 200, np.ceil(allv.max() * 200) / 200
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
        ax.set_title(f"{g.label}\nI_sw = {isw[g.key]:.0f} uA,  mean i_det = {np.nanmean(Z):.3f}", fontsize=10)
        ax.set_xlabel("x0 [nm]")
        ax.grid(False)
    axes[0].set_ylabel("y0 [nm]  (one unit cell)")
    if im is not None:
        cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
        cb.set_label("i_det = I_det / I_sw   (darker = less sensitive)")
    fig.suptitle(title)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_efficiency(geoms, vals, met, isw, title, path, styles=None):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    lo = min(0.75, np.floor(min(vals[g.key][0] for g in geoms) * 20) / 20 - 0.05)
    ib = np.linspace(lo, 1.0, 800)
    for g in geoms:
        ls = (styles or {}).get(g.key, "-")
        v, m = vals[g.key], met[g.key]
        a1.plot(ib, eta_curve(v, ib), color=g.color, ls=ls,
                label=f"{g.label}:  i_50 = {m['i_50']:.3f},  i_95 = {m['i_95']:.3f}")
        a1.axvline(m["i_95"], color=g.color, ls=":", lw=1.2)
        a2.plot(ib * isw[g.key], eta_curve(v, ib), color=g.color, ls=ls, label=f"{g.label}  (I_sw = {isw[g.key]:.0f} uA)")
        a2.axvline(isw[g.key], color=g.color, ls="--", lw=1.0)
    a1.set_xlabel("normalised bias  i_b = I_bias / I_sw   (dotted: i_95)")
    a2.set_xlabel("absolute bias  I_bias [uA]   (dashed: I_sw)")
    for a in (a1, a2):
        a.set_ylabel("area-integrated efficiency [%]")
        a.set_ylim(-3, 103)
        a.legend(loc="upper left", fontsize=8.5)
    a1.set_title("Normalised: window below the wall")
    a2.set_title("Absolute")
    fig.suptitle(title)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_selectivity(geoms, tabs, path):
    geoms = [g for g in geoms if tabs.get(g.key)]
    if not geoms:
        return
    fig, axes = _panels(len(geoms), figsize=(5.25 * len(geoms), 4.4), sharey=True)
    for ax, g in zip(axes, geoms):
        tab = tabs[g.key]
        fs = sorted({f for site in tab.values() for f in site})
        for site, d_ in sorted(tab.items()):
            xs = [f for f in fs if f in d_]
            ys = [d_[f][0] for f in xs]
            ax.plot(xs, ys, "-", color=INK2, lw=1.0, alpha=0.7)
            for f in xs:
                if d_[f][1] == "at_wall":
                    ax.plot(f, d_[f][0], "^", color=INK2, ms=7)
            ax.annotate(f"({site[0]},{site[1]})", (xs[0], ys[0]), textcoords="offset points", xytext=(-4, 0),
                        ha="right", fontsize=8, color=INK2)
        med = [np.median([d_[f][0] for d_ in tab.values() if f in d_]) for f in fs]
        ax.plot(fs, med, "-o", color=g.color, lw=2.6, ms=6, label="median over probe sites")
        ax.set_xlabel("hotspot energy factor f  (1 = photon)")
        ax.set_title(g.label)
        ax.legend(loc="lower right")
        ax.set_xlim(min(fs) - 0.12, 1.05)
    axes[0].set_ylabel("i_det   (triangle = no detection below I_sw)")
    fig.suptitle("Threshold vs hotspot energy: steeper rise towards 1 for weak seeds = more selective", y=0.99)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_crowding(geoms, npz_paths, bias, path, film):
    from .sim import smoothed_peak
    data = {g.key: np.load(npz_paths[g.key]) for g in geoms if os.path.exists(npz_paths[g.key])}
    geoms = [g for g in geoms if g.key in data]
    if not geoms:
        return
    fig, axes = _panels(len(geoms), figsize=(5.5 * len(geoms), 5.4), sharey=True, layout="constrained")
    ims = []
    for ax, g in zip(axes, geoms):
        d_ = data[g.key]
        pts, tri, K = d_["points"], d_["triangles"], d_["K"] / float(d_["K_uniform"])
        ims.append(ax.tripcolor(pts[:, 0] * 1e3, pts[:, 1] * 1e3, tri, K, shading="gouraud", cmap=SEQ_CMAP,
                                vmin=0.0, vmax=1.6))
        draw_holes(ax, g)
        sel = np.abs(pts[:, 1]) <= g.crowd_half_y_nm * 1e-3
        raw = float(K[sel].max())
        sm, loc = smoothed_peak(pts, K, g, 1.5 * film.xi_nm * 1e-3)
        ax.plot(loc[0] * 1e3, loc[1] * 1e3, "*", ms=13, mfc="white", mec=INK, zorder=5)
        ax.set_title(f"{g.label}\nsmoothed peak {sm:.2f}x uniform  (single-node max {raw:.2f}x)", fontsize=10)
        ax.set_xlim(-g.width_nm / 2, g.width_nm / 2)
        ax.set_ylim(-140, 140)
        ax.set_aspect("equal")
        ax.set_xlabel("x [nm]")
        ax.grid(False)
    axes[0].set_ylabel("y [nm]")
    fig.colorbar(ims[0], ax=axes, shrink=0.8, pad=0.02, label="|K| / (I / W)   (steady state before the photon)")
    fig.suptitle(f"Current crowding at {bias:.2f} I_sw  (star = interior disc-averaged peak, radius 1.5 xi)")
    fig.savefig(path, dpi=200)
    plt.close(fig)
