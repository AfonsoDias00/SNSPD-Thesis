"""Zipper and slot columns against both plain wires on one axis, in units of I_dep (no simulation).

The per-figure efficiency plots of hitmap.py use i_b = I_bias / I_sw. The dark I_sw is a slow instability whose value
moves with strip length, window and criterion (plain 200 nm: 97.7 uA at 0.8 um, 98.7 / 95.8 uA at 2 um, sustained /
trigger). That moves every curve sideways. The photon thresholds themselves are bisected in uA, so here they are
converted back to I_bias / I_dep(narrowest free width), which does not depend on I_sw.

  python efficiency_compare.py [--hotspot kinetic] [--preset nbtin]
                                  -> results/efficiency_compare_<hotspot>[-<preset>]_T400_m1.4.png
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from snspd import analysis as an
from snspd import geometry as geo
from snspd import plotting as pl
from snspd import sim, store
from snspd.units import Film

RUNS = [("results/hitmap", ["Plain_Wire", "Long_Zipper"]),
        ("results/slots", ["Plain_Wire_2um", "Slot_110", "Slot_200", "Slot_400"])]
LABELS = {"Plain_Wire": "Plain wire, 0.8 um", "Plain_Wire_2um": "Plain wire, 2 um", "Long_Zipper": "Zipper",
          "Slot_110": "Slots 110 nm", "Slot_200": "Slots 200 nm", "Slot_400": "Slots 400 nm"}
STYLES = {"Plain_Wire_2um": "--"}


def curves(hotspot, sustained, film, preset=""):
    """{key: sorted thresholds in units of I_dep(free width)}, {key: I_sw / I_dep}."""
    out, isw_out = {}, {}
    for out_dir, keys in RUNS:
        cfg = sim.SimConfig(tmax=400, mesh_xi=1.4, hotspot=hotspot, sustained=sustained, out_dir=out_dir,
                            preset=preset)
        rows = store.load_rows(cfg.path("hitmap", "csv"))
        isw = {k: v["I_sw"] for k, v in store.load_json(cfg.path("isw", "json"), {}).items()}
        for k in keys:
            g = geo.get(k)
            samples = an.photon_rows(rows, k)
            if k not in isw or not samples:
                continue
            idep = film.I_dep_uA(g.free_width_min_nm())
            out[k] = an.pixel_values(an.unit_map(g, samples)[2]) * isw[k] / idep
            isw_out[k] = isw[k] / idep
    return out, isw_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hotspot", default="kinetic")
    ap.add_argument("--preset", default="", choices=list(sim.PRESETS))
    a = ap.parse_args()
    pl.apply_style()
    film = Film()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True, layout="constrained")
    ib = np.linspace(0.70, 1.0, 900)
    for ax, sustained, title in [(axes[0], False, "Voltage trigger (lower bound)"),
                                 (axes[1], True, "Sustained switching (upper bound)")]:
        vals, isw = curves(a.hotspot, sustained, film, a.preset)
        for k, v in vals.items():
            ax.plot(ib, an.eta_curve(v, ib), color=geo.get(k).color, ls=STYLES.get(k, "-"), lw=1.8,
                    label=f"{LABELS[k]}:  i_50 = {np.quantile(v, 0.5):.3f}")
        plain = ", ".join(f"{isw[k]:.3f}" for k in ("Plain_Wire", "Plain_Wire_2um") if k in isw)
        ax.set_title(title)
        ax.set_xlabel("I_bias / I_dep (narrowest cross-section)")
        ax.text(0.02, 0.97, f"plain-wire dark I_sw / I_dep: {plain}", transform=ax.transAxes, va="top",
                fontsize=8, color=pl.INK2)
        ax.legend(loc="lower right", fontsize=8.5, frameon=False)
    axes[0].set_ylabel("share of the unit cell that detects [%]")
    axes[0].set_ylim(-3, 103)
    model = a.hotspot + (f"-{a.preset}" if a.preset else "")
    fig.suptitle(f"Detection efficiency vs bias, all devices on one axis [{model}, 400 tau0, 1.4 xi]")
    path = f"results/efficiency_compare_{model}_T400_m1.4.png"
    fig.savefig(path, dpi=200)
    print(path)


if __name__ == "__main__":
    main()
