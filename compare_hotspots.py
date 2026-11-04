"""
Compare photon thresholds between two hotspot models (default: conserving vs kinetic) site by site, from the stored
hitmap and narrow-strip CSVs, in units of I_dep of the narrowest free width.

  python compare_hotspots.py                                   # sustained criterion, 400 tau0, 1.4 xi
  python compare_hotspots.py --criterion trigger
  python compare_hotspots.py --dir results/slots --geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400
"""
import argparse
import os

import numpy as np

from snspd import geometry as geo
from snspd import store
from snspd.units import FILM


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
    a = ap.parse_args()

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
