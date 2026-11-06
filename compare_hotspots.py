"""
Compare photon thresholds between two hotspot models (default: conserving vs kinetic) site by site, from the stored
hitmap and narrow-strip CSVs, in units of I_dep of the narrowest free width.

  python compare_hotspots.py                                   # sustained criterion, 400 tau0, 1.4 xi
  python compare_hotspots.py --criterion trigger
  python compare_hotspots.py --dir results/slots --geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400
  python compare_hotspots.py --profiles results/hotspot_profiles.png   # the two epsilon(r, t) fields, no TDGL
"""
import argparse
import os

import numpy as np

from snspd import geometry as geo
from snspd import store
from snspd.units import FILM

PROFILE_TIMES = (5.0, 20.0, 50.0, 100.0, 200.0)     # tau0 after arrival, for the printed table


def profiles(a, path, geom_key="Plain_Wire", site=(0.0, 0.0), cuts=(10.0, 50.0)):
    """Figure and table of the hotspot fields of models a and b for one photon (f = 1): epsilon at the hit site,
    the integrated suppression  integral (1 - epsilon) dA  and the cut across the strip at two times."""
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from snspd import hotspot as hsm
    from snspd import plotting as pl
    from snspd import sim

    pl.apply_style()
    colors = {a.a: "#2a78d6", a.b: "#eb6834"}          # categorical slots 1 and 2 (validated pair)
    g, p = geo.get(geom_key), hsm.HOTSPOT
    hs = sim.context(g, a.mesh_xi)["hs"]
    x0, y0 = site
    k = int(np.argmin((hs.P[:, 0] * 1e3 - x0) ** 2 + (hs.P[:, 1] * 1e3 - y0) ** 2))
    ts = np.arange(0.0, 200.01, 0.5)
    tri = mtri.Triangulation(hs.P[:, 0] * 1e3, hs.P[:, 1] * 1e3, hs.T)
    xs = np.linspace(-g.width_nm / 2, g.width_nm / 2, 401)
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), layout="constrained")
    print(f"hotspot fields, f = 1, {g.label} hit ({x0:g}, {y0:g}) nm, times after arrival [tau0]")
    for model in (a.a, a.b):
        eps = hsm.build_eps(hs, model, x0, y0, 1.0, p.t0_tau0 + ts[-1] + 1.0)
        E = np.array([eps(hs.P, t=p.t0_tau0 + t) for t in ts])
        S = ((1.0 - E) * hs.m).sum(1) * 1e6                                  # um^2 -> nm^2
        c = colors[model]
        ax[0].plot(ts, E[:, k], color=c, label=model)
        ax[1].plot(ts, S, color=c, label=model)
        for t, ls in zip(cuts, ("-", "--")):
            cut = mtri.LinearTriInterpolator(tri, E[int(round(t / 0.5))])(xs, np.full_like(xs, y0))
            ax[2].plot(xs, cut, color=c, ls=ls, label=f"{model}, t0 + {t:g}")
        i = int(np.argmin(E[:, k]))
        normal = ts[(E < 0).any(1)]
        print(f"  {model:<11} min epsilon at the site {E[i, k]:+.3f} at {ts[i]:g};  epsilon < 0 somewhere until "
              + (f"{normal[-1]:g}" if len(normal) else "never")
              + ";  integral(1 - eps) dA [nm^2] at " + "/".join(f"{t:g}" for t in PROFILE_TIMES) + ": "
              + "  ".join(f"{S[int(round(t / 0.5))]:.0f}" for t in PROFILE_TIMES))
    for axis in ax[:2]:
        axis.set_xlabel("time after photon arrival  t - t0  [tau0]")
        axis.set_xlim(0, ts[-1])
    ax[0].axhline(0.0, color=pl.INK2, lw=0.8, ls=":")
    ax[0].text(ts[-1], 0.03, "normal (epsilon < 0) below", ha="right", va="bottom", fontsize=8, color=pl.INK2)
    ax[0].set_ylabel("epsilon at the hit site")
    ax[0].set_title("depth: instant core vs delayed, saturating build-up", fontsize=10)
    ax[1].set_ylabel("integral of (1 - epsilon) over the film  [nm^2]")
    ax[1].set_title("weight: relaxes over t_relax vs escapes over t_esc", fontsize=10)
    ax[2].set_xlabel("x across the strip  [nm]")
    ax[2].set_ylabel(f"epsilon along y = {y0:g} nm")
    ax[2].set_xlim(xs[0], xs[-1])
    ax[2].set_title("cut across the strip", fontsize=10)
    for axis, loc in zip(ax, ("lower right", "upper right", "lower right")):
        axis.legend(fontsize=8, loc=loc)
    fig.suptitle(f"Photon hotspot, f = 1, {g.label} hit at ({x0:g}, {y0:g}) nm:  {a.a} vs {a.b}", fontsize=11)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"wrote {path}")


def tag(model, a):
    return f"{model}_T{a.tmax:g}_m{a.mesh_xi:g}" + ("_trigger" if a.criterion == "trigger" else "")


def rows_by_site(path, energy=1.0):
    out = {}
    for r in store.load_rows(path):
        if abs(float(r["energy_f"]) - energy) < 1e-9:
            out[(r["Geometry"], float(r["x_nm"]), float(r["y_nm"]))] = r
    return out


def i_dep(geom_key, width_by_key):
    return FILM.I_dep_uA(width_by_key[geom_key])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="conserving")
    ap.add_argument("--b", default="kinetic")
    ap.add_argument("--criterion", choices=["sustained", "trigger"], default="sustained")
    ap.add_argument("--tmax", type=float, default=400.0)
    ap.add_argument("--mesh-xi", type=float, default=1.4)
    ap.add_argument("--dir", default="results/hitmap")
    ap.add_argument("--geoms", default="Plain_Wire,Long_Zipper")
    ap.add_argument("--a-dir", default=None, help="directory of model a if different (e.g. results/v2_voltage_trigger/hitmap)")
    ap.add_argument("--narrow", default="30,56,100", help="narrow-strip widths to include ('' = none)")
    ap.add_argument("--profiles", default=None, metavar="PNG",
                    help="instead: plot the two models' epsilon(r, t) for a centre hit of the plain wire")
    a = ap.parse_args()
    if a.profiles:
        return profiles(a, a.profiles)

    widths = {k: geo.get(k).free_width_min_nm() for k in a.geoms.split(",") if k}
    pairs = [(os.path.join(a.a_dir or a.dir, f"hitmap_{tag(a.a, a)}.csv"), os.path.join(a.dir, f"hitmap_{tag(a.b, a)}.csv"))]
    for w in [s for s in a.narrow.split(",") if s]:
        widths["Narrow_Wire_" + w] = float(w)
        nb_a = "results/narrow_baseline" if a.criterion == "sustained" or a.a == "kinetic" else "results/v2_voltage_trigger/narrow_baseline"
        ta = tag(a.a, a).replace("_trigger", "") if nb_a.startswith("results/v2") else tag(a.a, a)
        pairs.append((os.path.join(nb_a, f"narrow{w}nm_hitmap_{ta}.csv"),
                      os.path.join("results/narrow_baseline", f"narrow{w}nm_hitmap_{tag(a.b, a)}.csv")))

    print(f"photon thresholds at f = 1 in units of I_dep (narrowest free width)   criterion = {a.criterion}")
    print(f"{'geometry':<18}{'x_nm':>7}{'y_nm':>7}{a.a:>12}{a.b:>10}{'change':>9}   flags / trigger time [tau0] {a.a} -> {a.b}")
    summary = {}
    for pa, pb in pairs:
        if not (os.path.exists(pa) and os.path.exists(pb)):
            continue
        ra, rb = rows_by_site(pa), rows_by_site(pb)
        for key in sorted(set(ra) & set(rb)):
            g, x, y = key
            gk = g if g in widths else None
            if g == "Narrow_Wire":
                w = os.path.basename(pb).split("nm_")[0].replace("narrow", "")
                gk = "Narrow_Wire_" + w
            if gk is None:
                continue
            Id = i_dep(gk, widths)
            va, vb = float(ra[key]["I_det_uA"]) / Id, float(rb[key]["I_det_uA"]) / Id
            fa, fb = ra[key]["flag"] or "-", rb[key]["flag"] or "-"
            ta_, tb_ = ra[key].get("t_trigger") or "", rb[key].get("t_trigger") or ""
            print(f"{gk:<18}{x:7.0f}{y:7.0f}{va:12.3f}{vb:10.3f}{vb - va:+9.3f}   {fa}/{fb}  {ta_} -> {tb_}")
            summary.setdefault(gk, []).append((va, vb))
    print("\nper geometry: median and range of the change (b - a) over common sites")
    for gk, v in summary.items():
        d = np.array([b - a_ for a_, b in v])
        print(f"  {gk:<18} n = {len(d):3d}   median {np.median(d):+.3f}   min {d.min():+.3f}   max {d.max():+.3f}")


if __name__ == "__main__":
    main()
