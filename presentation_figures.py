"""
Figures for the presentation (docs/presentation_outline.md) that the analysis scripts do not make.  Reads only stored
results; no simulation.  Writes to docs/fig/.

  python presentation_figures.py            # all four
  python presentation_figures.py devices width validation sensitivity
"""
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from snspd import geometry as geo
from snspd import plotting as pl
from snspd.units import FILM

OUT = "docs/fig"
# categorical slots 1-3 (reference palette, validated all-pairs in both modes) and the reserved status colours
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GOOD, WARNING, SERIOUS, CRITICAL, NEUTRAL = "#0ca30c", "#fab219", "#ec835a", "#d03b3b", "#b9b8b2"


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


# ----------------------------------------------------------------------------------------------------------------------
def fig_devices():
    """All devices to one nm scale, centre 0.8 um of each strip (the whole 0.8 um devices)."""
    half = 400.0
    devices = [geo.PLAIN_WIRE, geo.get("Long_Zipper"), geo.get("Slot_110"), geo.get("Slot_200"), geo.get("Slot_400"),
               geo.plain_strip(100), geo.plain_strip(56), geo.plain_strip(30)]
    gap = 140.0
    narrow_gap = 175.0
    fig, ax = plt.subplots(figsize=(13, 6.2))
    x = 0.0
    for i, g in enumerate(devices):
        if i == 5:
            x += 120.0                                          # separate the narrow-strip controls
        if i >= 5:
            gap = narrow_gap
        w = g.width_nm
        cx = x + w / 2
        ax.add_patch(Rectangle((x, -half), w, 2 * half, facecolor="#cfe0f5", edgecolor=pl.INK, lw=0.9))
        for h in g.holes:
            pts = h.outline(120)
            if np.abs(pts[:, 1]).min() > half:
                continue
            pts = pts.copy()
            pts[:, 1] = np.clip(pts[:, 1], -half, half)
            ax.fill(pts[:, 0] + cx, pts[:, 1], facecolor=pl.SURFACE, edgecolor=pl.INK, lw=0.8)
        wf = g.free_width_min_nm()
        length = "2 µm" if g.length_nm > 1000 else "0.8 µm"
        name = {"Long_Zipper": "long zipper", "Plain_Wire": "plain wire"}.get(
            g.key, g.label.replace("Slot column, ", "slots ").replace("Plain ", "plain\n"))
        ax.text(cx, half + 25, name, ha="center", va="bottom", fontsize=9.5, color=pl.INK)
        ax.text(cx, -half - 25, f"w = {w:g} nm\nfree {wf:.0f} nm\nI_dep {FILM.I_dep_uA(wf):.1f} µA\n{length} long",
                ha="center", va="top", fontsize=8, color=pl.INK2, linespacing=1.35)
        x += w + gap
    ax.annotate("", xy=(-60, half), xytext=(-60, -half), arrowprops=dict(arrowstyle="->", color=pl.INK2, lw=1.2))
    ax.text(-75, 0, "current", rotation=90, ha="right", va="center", fontsize=9, color=pl.INK2)
    ax.plot([x - gap - 100, x - gap], [-half - 200] * 2, color=pl.INK, lw=2)
    ax.text(x - gap - 50, -half - 215, "100 nm", ha="center", va="top", fontsize=8.5)
    ax.set_xlim(-110, x - gap + 20)
    ax.set_ylim(-half - 260, half + 210)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Devices to scale (centre 0.8 µm of each strip shown; film 9 nm thick, ξ = 5 nm)", fontsize=11)
    ax.text(x - gap - 250, half + 175, "narrow-strip controls", ha="center", fontsize=9, color=pl.INK2)
    save(fig, "devices.png")


# ----------------------------------------------------------------------------------------------------------------------
def read_width_scan(path):
    rows = {}
    for line in open(path):
        m = re.match(r"(i_50|i_95)\s+\[I_dep\]\s+(.*)", line.strip())
        if m:
            rows[m.group(1)] = [float(v) for v in m.group(2).split()]
    return rows      # values for 30, 56, 100, 200 nm, zipper


def fig_width():
    NB, V2 = "results/narrow_baseline", "results/v2_voltage_trigger/narrow_baseline"
    runs = [  # (label, file, colour, linestyle)
        ("lower bound (voltage trigger), conserving", f"{V2}/width_scan_summary_conserving_T400_m1.4.txt", BLUE, "-"),
        ("lower bound (voltage trigger), kinetic", f"{NB}/width_scan_summary_kinetic_T400_m1.4_trigger.txt", BLUE, "--"),
        ("upper bound (sustained), conserving", f"{NB}/width_scan_summary_conserving_T400_m1.4.txt", ORANGE, "-"),
        ("upper bound (sustained), kinetic", f"{NB}/width_scan_summary_kinetic_T400_m1.4.txt", ORANGE, "--"),
    ]
    widths = [30, 56, 100, 200]
    wz = geo.get("Long_Zipper").free_width_min_nm()
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True, layout="constrained")
    for ax, q, title in zip(axes, ("i_50", "i_95"), ("onset  i_50", "plateau  i_95")):
        for label, path, c, ls in runs:
            v = read_width_scan(path)[q]
            mk = "o" if ls == "-" else "s"
            ax.plot(widths, v[:4], color=c, ls=ls, marker=mk, ms=7, mfc=c if ls == "-" else pl.SURFACE, mew=1.6,
                    label=label)
            ax.plot([wz * (0.975 if ls == "-" else 1.025)], [v[4]], marker="*", ms=14, color=c, mfc=c if ls == "-" else pl.SURFACE, mew=1.4, ls="none")
        ax.set_xscale("log")
        ax.set_xticks(widths + [wz])
        ax.set_xticklabels(["30", "56", "100", "200", f"{wz:.0f}\nzipper"])
        ax.minorticks_off()
        ax.set_xlabel("strip width (zipper: narrowest free width)  [nm]")
        ax.set_title(title, fontsize=11)
    axes[0].set_ylabel("threshold  [I_dep]")
    axes[0].annotate("narrow strips heal without heat:\nthe upper bound is a floor", xy=(56, 0.81), xytext=(60, 0.55),
                     fontsize=8.5, color=pl.INK2, arrowprops=dict(arrowstyle="->", color=pl.INK2, lw=0.8))
    handles = [Line2D([], [], color=BLUE, lw=2, label="lower bound (voltage trigger)"),
               Line2D([], [], color=ORANGE, lw=2, label="upper bound (sustained)"),
               Line2D([], [], color=pl.INK2, lw=2, marker="o", label="conserving hotspot"),
               Line2D([], [], color=pl.INK2, lw=2, ls="--", marker="s", mfc=pl.SURFACE, label="kinetic hotspot"),
               Line2D([], [], color=pl.INK2, marker="*", ms=12, ls="none", label="long zipper")]
    axes[1].legend(handles=handles, fontsize=8.5, loc="lower right")
    fig.suptitle("Photon threshold vs strip width, f = 1, 400 τ0, 1.4 ξ", fontsize=12)
    save(fig, "width_scan.png")


# ----------------------------------------------------------------------------------------------------------------------
FAIL_CLASS = {  # failed check -> (class, short reason); classes: revised, limit, issue
    "1d.2": ("revised", "i_det(f) not linear at the edge"),
    "2a.5": ("limit", "gate misses cores < 0.4 nm from loop"),
    "2a.8": ("limit", "error if ψ = 0 exactly on loop"),
    "2a.9": ("limit", "bridge vortex counted twice, 1 frame"),
    "2b.5": ("revised", "ring overshoots to n = 2-4 (hysteresis)"),
    "2b.6": ("revised", "several slips between frames in bursts"),
    "2d.2": ("limit", "11 % multi-quantum steps (frame rate)"),
    "2d.4": ("limit", "same for Slot 400"),
    "2d.8": ("issue", "zipper meshes differently elsewhere"),
    "2d.9": ("issue", "flux sites move 0.04 I_dep (now pinned)"),
}
STATES = {  # key -> (colour, symbol, legend label)
    "PASS": (GOOD, "✓", "pass"),
    "INFO": (NEUTRAL, "i", "information only"),
    "revised": (WARNING, "~", "fail: the expectation was wrong"),
    "limit": (SERIOUS, "!", "fail: known limit of the fluxoid meter"),
    "issue": (CRITICAL, "✗", "fail: real issue (fixed by mesh pinning)"),
}


def fig_validation():
    groups, order = {}, []
    for line in open("results/validation/summary.txt"):
        m = re.match(r"^(\d[a-z])\s{2}(.+?)\s+\(\d+ s\)", line)
        if m:
            order.append(m.group(1))
            groups[m.group(1)] = (m.group(2), [])
            continue
        m = re.match(r"^\s+\[(PASS|FAIL|INFO)\]\s+(\d[a-z]\.\d+)", line)
        if m:
            st, cid = m.group(1), m.group(2)
            groups[cid.split(".")[0]][1].append((cid, FAIL_CLASS[cid][0] if st == "FAIL" else st))
    for gid in order:                                  # numeric order; ids used twice by the suite get a prime
        checks = sorted(groups[gid][1], key=lambda t: int(t[0].split('.')[1]))
        seen, out = set(), []
        for cid, st in checks:
            num = cid.split('.')[1]
            out.append((cid, num + ('′' if num in seen else ''), st))
            seen.add(num)
        groups[gid] = (groups[gid][0], out)
    ncol = max(len(v[1]) for v in groups.values())
    fig, ax = plt.subplots(figsize=(13, 0.62 * len(order) + 2.6))
    for r, gid in enumerate(order):
        title, checks = groups[gid]
        y = -r
        ax.text(-0.25, y, f"{gid}  {title}", ha="right", va="center", fontsize=9.5)
        for c, (cid, num, st) in enumerate(checks):
            col, sym, _ = STATES[st]
            ink = pl.INK if st in ('revised', 'INFO') else 'white'
            ax.add_patch(Rectangle((c + 0.04, y - 0.4), 0.92, 0.8, facecolor=col, edgecolor=pl.SURFACE, lw=2))
            ax.text(c + 0.5, y + 0.08, sym, ha="center", va="center", fontsize=11, color=ink, fontweight="bold")
            ax.text(c + 0.5, y - 0.26, num, ha="center", va="center", fontsize=6.5, color=ink)
    # failure notes under the grid
    notes = [f"{cid}  {why}" for cid, (_, why) in FAIL_CLASS.items()]
    y0 = -len(order) - 0.3
    for i, n in enumerate(notes):
        ax.text((i % 2) * ncol / 2, y0 - 0.45 * (i // 2), n, fontsize=8.2, color=pl.INK2, va="top")
    sts = [st for v in groups.values() for _, _, st in v[1]]
    n_pass, n_info = sts.count("PASS"), sts.count("INFO")
    n_fail = len(sts) - n_pass - n_info
    ax.legend(handles=[Patch(facecolor=STATES[k][0], label=f"{STATES[k][1]}  {STATES[k][2]}") for k in STATES],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=8.5, handlelength=1.2)
    ax.set_title(f"Validation suite (3 Nov 2026): {n_pass} pass, {n_fail} fail, {n_info} information only",
                 fontsize=11, pad=52)
    ax.set_xlim(-0.1, ncol + 0.1)
    ax.set_ylim(y0 - 0.45 * ((len(notes) + 1) // 2) - 0.2, 0.6)
    ax.set_aspect("equal")
    ax.axis("off")
    save(fig, "validation_grid.png")


# ----------------------------------------------------------------------------------------------------------------------
def fig_sensitivity():
    rows = []
    for line in open("results/kinetic_sensitivity.txt"):
        if re.match(r"^(plain|zipper|\d+ nm)", line) and "kinetic" not in line[:40]:
            vals = [float(v) for v in re.findall(r"\b\d\.\d{3}\b", line[33:])]
            if len(vals) == 5:
                crit = line[22:33].strip()
                rows.append((line[:22].strip() + ("  · upper" if crit == "sustained" else "  · lower"), vals))
    names = [("kinetic (stored)", pl.INK, "o"), ("7 / 75 τ0", BLUE, "D"), ("7 / 50 τ0", ORANGE, "s"),
             ("no ε floor", AQUA, "^")]
    fig, ax = plt.subplots(figsize=(10.5, 6.2), layout="constrained")
    for i, (label, v) in enumerate(rows):
        y = -i
        cons, kin = v[0], v[1:]
        d = np.array(kin) - cons
        ax.plot([d.min(), d.max()], [y, y], color=pl.GRID, lw=6, solid_capstyle="round", zorder=1)
        for (nm, c, mk), dv, dy in zip(names, d, (0.21, 0.07, -0.07, -0.21)):
            ax.plot(dv, y + dy, marker=mk, ms=8, color=c, mec=pl.SURFACE, mew=1.2, ls="none", zorder=3)
        ax.text(0.112, y, f"{d.max() - d.min():.3f}", va="center", fontsize=8.5, color=pl.INK2)
    ax.axvline(0, color=pl.INK2, lw=1)
    ax.axvspan(-0.015, 0.015, color=pl.GRID, alpha=0.35, zorder=0)
    ax.text(0, 0.55, "band 0.03 wide\n= the go margin", ha="center", va="bottom", fontsize=8, color=pl.INK2)
    ax.text(0.112, 0.55, "spread", fontsize=8.5, color=pl.INK2, va="bottom")
    ax.set_yticks([-i for i in range(len(rows))])
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_ylim(-len(rows) + 0.5, 1.2)
    ax.set_xlim(-0.05, 0.13)
    ax.set_xlabel("threshold minus the conserving-hotspot threshold  [I_dep]")
    ax.grid(axis="y", visible=False)
    ax.legend(handles=[Line2D([], [], marker=mk, color=c, ls="none", ms=8, label=nm) for nm, c, mk in names],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=8.5,
              title="kinetic variant (t_rise / t_esc); 0 = conserving hotspot; rows: site · bound", title_fontsize=8.5)
    ax.set_title("Kinetic hotspot: how much the choices Simon et al. leave open move the thresholds (f = 1)",
                 fontsize=11, pad=48)
    save(fig, "sensitivity.png")


# ----------------------------------------------------------------------------------------------------------------------
def _eps_at_site(model, ts, geom_key="Plain_Wire", site=(0.0, 0.0)):
    """epsilon at the hit site vs time after arrival, f = 1 (the hotspot field only, no TDGL)."""
    from snspd import hotspot as hsm
    from snspd import sim
    hs, p = sim.context(geo.get(geom_key), 1.4)["hs"], hsm.HOTSPOT
    k = int(np.argmin((hs.P[:, 0] * 1e3 - site[0]) ** 2 + (hs.P[:, 1] * 1e3 - site[1]) ** 2))
    eps = hsm.build_eps(hs, model, site[0], site[1], 1.0, p.t0_tau0 + ts[-1] + 1.0)
    return np.array([eps(hs.P, t=p.t0_tau0 + t)[k] for t in ts])


def fig_gl():
    """Slide: the GL equation pyTDGL solves, with epsilon(r, t) marked, and how the two hotspot models set it."""
    fig = plt.figure(figsize=(13.33, 7.5))
    INK, INK2, HL = pl.INK, pl.INK2, "#fde9b5"

    # --- the equation, epsilon highlighted ---------------------------------------------------------------------
    pre = (r"$\frac{u}{\sqrt{1+\gamma^2|\psi|^2}}\left(\frac{\partial}{\partial t}+i\mu"
           r"+\frac{\gamma^2}{2}\frac{\partial |\psi|^2}{\partial t}\right)\psi\ =\ ($")
    eps_s = r"$\epsilon(\mathbf{r},t)$"
    post = r"$\ -\ |\psi|^2)\,\psi\ +\ (\nabla - i\mathbf{A})^2\psi$"
    fs, y = 25, 0.865
    r = fig.canvas.get_renderer()
    parts = []
    for s in (pre, eps_s, post):
        t = fig.text(0, y, s, fontsize=fs, va="baseline", color=INK)
        parts.append((t, t.get_window_extent(r).width / fig.bbox.width))
    x = 0.5 - sum(w for _, w in parts) / 2
    for i, (t, w) in enumerate(parts):
        t.set_x(x)
        if i == 1:
            bb = t.get_window_extent(r).transformed(fig.transFigure.inverted())
            fig.patches.append(Rectangle((x - 0.004, bb.y0 - 0.012), w + 0.008, bb.height + 0.024,
                                         transform=fig.transFigure, facecolor=HL, edgecolor="none", zorder=-1))
        x += w
    fig.text(0.5, 0.965, "Generalised TDGL equation solved by pyTDGL (dimensionless)", ha="center", fontsize=15,
             color=INK, weight="bold")
    fig.text(0.5, 0.755, r"$\epsilon(\mathbf{r},t)$ sets how superconducting each point is: at rest $|\psi|^2 = \epsilon$."
             r"   The photon enters only here, as a dip in $\epsilon$.", ha="center", fontsize=13, color=INK2)

    # --- epsilon at the hit site vs time, with the three regimes ------------------------------------------------
    ax = fig.add_axes([0.07, 0.10, 0.42, 0.55])
    ts = np.arange(0.0, 150.01, 0.5)
    ax.axhspan(-1.7, 0, color="#ebe9e4", zorder=0)
    ax.axhline(0, color=INK2, lw=0.8)
    for y_, s in [(1.04, r"$\epsilon = 1$: fully superconducting"), (0.27, r"$0 < \epsilon < 1$: weakened, still no resistance"),
                  (-1.1, r"$\epsilon \leq 0$: normal ($\psi = 0$, resistive)")]:
        ax.text(148, y_, s, ha="right", va="bottom", fontsize=10.5, color=INK2)
    for model, c, lab in [("conserving", BLUE, "conserving"), ("kinetic", ORANGE, "kinetic (default)")]:
        ax.plot(ts, _eps_at_site(model, ts), color=c, lw=2.4, label=lab)
    ax.set_xlim(0, 150)
    ax.set_ylim(-1.7, 1.25)
    ax.set_xlabel(r"time after the photon arrives  [$\tau_0 \approx 0.15$ ps]", fontsize=11)
    ax.set_ylabel(r"$\epsilon$ at the photon's position", fontsize=11)
    ax.set_title("Centre hit of a 200 nm strip, f = 1", fontsize=11.5, color=INK)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.10), fontsize=10.5, frameon=False)

    # --- the two models ------------------------------------------------------------------------------------------
    def card(y0, h, color, title, lines):
        fig.patches.append(Rectangle((0.555, y0), 0.405, h, transform=fig.transFigure, facecolor="none",
                                     edgecolor=color, lw=2))
        fig.text(0.572, y0 + h - 0.045, title, fontsize=13.5, weight="bold", color=INK)
        for i, s in enumerate(lines):
            fig.text(0.572, y0 + h - 0.105 - 0.058 * i, s, fontsize=12.5, color=INK if s.startswith("$") else INK2)

    card(0.395, 0.27, BLUE, "Conserving hotspot (until 4 Nov)", [
        r"$\epsilon = 1 + A\,f\,e^{-t/\tau_r}\,n(\mathbf{r},t),\quad \partial_t n = D\nabla^2 n$",
        r"$A = -2.5,\ \ \tau_r = 35\,\tau_0,\ \ D = 10\ \mathrm{nm^2}/\tau_0$",
        r"Instant normal core ($\epsilon = -1.5$), heals in $\approx 6\,\tau_0$"])
    card(0.06, 0.315, ORANGE, "Kinetic hotspot (Simon et al. 2025)", [
        r"$\epsilon = 1 - s,\quad \partial_t s = [\,s_{eq}(e) - s\,]\,/\,\tau_{th}$",
        r"$e = f\,n(\mathbf{r},t)\,e^{-t/\tau_{esc}},\quad s_{eq} = s_{max}\tanh(|A|\,e/s_{max})$",
        r"$\tau_{th} = 10\,\tau_0,\ \ \tau_{esc} = 100\,\tau_0$  (same $n$, same $D$)",
        r"Gradual, never normal ($\epsilon_{min} = 0.43$), lasts $\approx 100\,\tau_0$"])
    save(fig, "gl_hotspot.png")


FIGS = {"devices": fig_devices, "width": fig_width, "validation": fig_validation, "sensitivity": fig_sensitivity,
        "gl": fig_gl}

if __name__ == "__main__":
    pl.apply_style()
    for name in sys.argv[1:] or list(FIGS):
        FIGS[name]()
