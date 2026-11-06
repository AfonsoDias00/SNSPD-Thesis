#!/usr/bin/env python3
"""
hotspot_sensitivity.py -- photon thresholds of the kinetic hotspot under other choices of its open parameters.

Writes results/kinetic_sensitivity.csv (one row per variant and case) and results/kinetic_sensitivity.txt (table).
f = 1, 400 tau0, 1.4 xi, pinned production meshes.  Rows computed on 4 Nov 2026 by scratch scripts were imported
into the CSV unchanged (variants D-anchor, esc10ps, no-sat).

  python hotspot_sensitivity.py run --variants NbTiN --workers 4
  python hotspot_sensitivity.py table

Variants (stored kinetic results: t_rise 10, t_esc 100 tau0, D 10 nm^2/tau0, epsilon >= -0.25; tau0 = 0.12-0.17 ps):
  D-anchor  t_rise 7, t_esc 75 tau0: tau0 = 0.200 ps, where D = 10 nm^2/tau0 = 0.5 cm^2/s (tau_th 1.4, tau_esc 15 ps)
  esc10ps   t_rise 7, t_esc 50 tau0: same anchor, tau_esc = 10 ps (appendix VII.1 of Simon et al.)
  no-sat    no floor on epsilon (s_eq = |A| e, as the paper's unbounded alpha)
  NbTiN     the times measured on 9 nm NbTiN (KTH films) by Sidorova et al., PRB 104, 184514 (2021) and
            arXiv:2204.01423: tau0 = hbar/(pi Delta0) = 0.15 ps (T_c 8.6-9.5 K; with R_sq = 381 Ohm it implies
            lambda ~ 660 nm, close to the 600 nm used), D = 0.472 cm^2/s = 7.1 nm^2/tau0, energy escape 80-90 ps
            (phonon escape tau_esc = 80 ps, electron energy relaxation tau_eps(T_c) = 90 ps) = 530 tau0; tau_th is not
            measured for NbTiN, so t_rise stays 10 tau0 (1.5 ps, as NbN in Simon et al.).
"""
from __future__ import annotations

import argparse
import os
from dataclasses import replace

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")

from snspd import geometry as geo
from snspd import hotspot as hsm
from snspd import sim, store
from snspd.units import FILM

OUT_CSV, OUT_TXT = "results/kinetic_sensitivity.csv", "results/kinetic_sensitivity.txt"
FIELDS = ["variant", "case", "criterion", "geometry", "width_nm", "x_nm", "y_nm", "I_det_uA", "flag", "t_trigger",
          "branch"]
VARIANTS = {
    "D-anchor": dict(kin_rise_tau0=7.0, kin_esc_tau0=75.0),
    "esc10ps": dict(kin_rise_tau0=7.0, kin_esc_tau0=50.0),
    "no-sat": dict(kin_s_max=1e3),
    "NbTiN": dict(kin_rise_tau0=10.0, kin_esc_tau0=530.0, D_nm2_per_tau0=7.1),
}
NB, V2, HM = "results/narrow_baseline", "results/v2_voltage_trigger", "results/hitmap"


def cases():
    """(case label, criterion, geometry, I_sw json, conserving csv, kinetic csv, x, y)."""
    out = []
    for gk, lab in (("Plain_Wire", "plain 200 nm"), ("Long_Zipper", "zipper")):
        sites = ((0.0, 0.0), (-96.0, 0.0)) if gk == "Plain_Wire" else ((0.0, 0.0), (-68.0, 0.0), (-96.0, 0.0))
        for x, y in sites:
            out.append((f"{lab} ({x:g},{y:g})", "sustained", geo.get(gk), f"{HM}/isw_kinetic_T400_m1.4.json",
                        f"{HM}/hitmap_conserving_T400_m1.4.csv", f"{HM}/hitmap_kinetic_T400_m1.4.csv", x, y))
    for w in (30, 56):
        for crit in ("sustained", "trigger"):
            suf = "" if crit == "sustained" else "_trigger"
            cons = (f"{NB}/narrow{w}nm_hitmap_conserving_T400_m1.4.csv" if crit == "sustained" else
                    f"{V2}/narrow_baseline/narrow{w}nm_hitmap_conserving_T400_m1.4.csv")
            out.append((f"{w} nm strip (0,0)", crit, geo.plain_strip(w), f"{NB}/narrow{w}nm_isw_kinetic_T400_m1.4{suf}.json",
                        cons, f"{NB}/narrow{w}nm_hitmap_kinetic_T400_m1.4{suf}.csv", 0.0, 0.0))
    for x, y in ((-36.0, 61.0), (-42.0, 58.0)):
        out.append((f"zipper ({x:g},{y:g})", "sustained", geo.LONG_ZIPPER, f"{HM}/isw_kinetic_T400_m1.4.json",
                    f"{HM}/hitmap_conserving_T400_m1.4.csv", f"{HM}/hitmap_kinetic_T400_m1.4.csv", x, y))
    return out


def stored(path, gkey, x, y):
    for r in store.load_rows(path):
        if r["Geometry"] == gkey and float(r["x_nm"]) == x and float(r["y_nm"]) == y and float(r["energy_f"]) == 1.0:
            return float(r["I_det_uA"])
    return float("nan")


def cmd_run(a):
    done = {(r["variant"], r["case"], r["criterion"]) for r in store.load_rows(OUT_CSV)}
    names = [v for v in a.variants.split(",") if v] or list(VARIANTS)
    tasks = []
    for name in names:
        params = replace(hsm.HOTSPOT, **VARIANTS[name])
        for lab, crit, g, isw, _, _, x, y in cases():
            if (name, lab, crit) in done:
                continue
            cfg = sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", sustained=crit == "sustained",
                                out_dir="results")
            tasks.append(dict(kind="hit", geom=g, cfg=cfg, params=params, x_nm=x, y_nm=y, f=1.0,
                              I_sw=store.load_json(isw, {})[g.key]["I_sw"], variant=name, case=lab, criterion=crit))
    sim.log(f"-- hotspot sensitivity: {len(tasks)} thresholds")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['variant']} {t['case']} {t['criterion']}: {r['error']}")
            return
        res = r["res"]
        store.append_row(OUT_CSV, dict(variant=t["variant"], case=t["case"], criterion=t["criterion"],
                                       geometry=t["geom"].key, width_nm=t["geom"].width_nm, x_nm=t["x_nm"],
                                       y_nm=t["y_nm"], I_det_uA=res["I_det"], flag=res["flag"] or "",
                                       t_trigger=res["t_trigger"], branch=res["branch"]), FIELDS)
        sim.log(f"   {t['variant']} {t['case']} {t['criterion']}: {res['I_det']:.3f} uA  ({r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result)
    cmd_table(a)


def cmd_table(a):
    rows = store.load_rows(OUT_CSV)
    V = list(VARIANTS)
    L = ["Kinetic hotspot: sensitivity of photon thresholds to the choices Simon et al. (arXiv:2501.13791) leave open,",
         "and to the time scales measured on 9 nm NbTiN (variant NbTiN).  f = 1, 400 tau0 window, 1.4 xi, production",
         "(pinned) meshes.  Thresholds in units of I_dep of the narrowest free width; [t] = trigger time minus t0 [tau0].",
         "Written by hotspot_sensitivity.py (variants: see its docstring).", "",
         "  kinetic   stored results: t_rise 10, t_esc 100 tau0, D 10 nm^2/tau0, epsilon >= -0.25 (tau0 = 0.12-0.17 ps)",
         "  D-anchor  t_rise  7, t_esc  75 tau0 (tau0 = 0.200 ps, where D = 10 nm^2/tau0 = 0.5 cm^2/s)",
         "  esc10ps   t_rise  7, t_esc  50 tau0 (same anchor, tau_esc = 10 ps of the paper's appendix VII.1)",
         "  no-sat    t_rise 10, t_esc 100 tau0, no floor on epsilon (s_eq = |A| e, as the paper's unbounded alpha)",
         "  NbTiN     t_rise 10, t_esc 530 tau0, D 7.1 nm^2/tau0 (measured on 9 nm NbTiN, tau0 = 0.15 ps; 4-5 Nov 2026)",
         "",
         f"{'site':<22}{'criterion':<11}{'conserving':>11}{'kinetic':>9}" + "".join(f"{v:>16}" for v in V)]
    spread = []
    for lab, crit, g, _, cons, kin, x, y in cases():
        Id = FILM.I_dep_uA(g.free_width_min_nm())
        c, k = stored(cons, g.key, x, y) / Id, stored(kin, g.key, x, y) / Id
        row, allv = f"{lab:<22}{crit:<11}{c:11.3f}{k:9.3f}", [k]
        for v in V:
            hit = [r for r in rows if (r["variant"], r["case"], r["criterion"]) == (v, lab, crit)]
            if hit:
                i = float(hit[-1]["I_det_uA"]) / Id
                t = hit[-1]["t_trigger"]
                allv.append(i)
                row += f"{i:10.3f} [{float(t) - 20:3.0f}]" if t not in ("", "None") else f"{i:10.3f}      "
            else:
                row += f"{'n/a':>16}"
        spread.append((lab, crit, c, k, min(allv), max(allv)))
        L.append(row)
    L += ["", "range over the kinetic variants (stored kinetic included) vs the conserving value:"]
    for lab, crit, c, k, lo, hi in spread:
        L.append(f"  {lab:<22}{crit:<11} kinetic {lo:.3f}-{hi:.3f}   conserving {c:.3f}   "
                 f"change {lo - c:+.3f} to {hi - c:+.3f}")
    text = "\n".join(L)
    sim.log(text)
    with open(OUT_TXT, "w") as fh:
        fh.write(text + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Variants")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "table"])
    ap.add_argument("--variants", default="", help="comma-separated (default: all)")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    {"run": cmd_run, "table": cmd_table}[a.command](a)


if __name__ == "__main__":
    main()
