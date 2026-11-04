"""
1d  Photon thresholds of the plain 200 nm wire against known behaviour (production model: conserving hotspot,
400 tau0, 1.4 xi, sustained switching; I_sw and the f = 0.35 / 0.6 / 1 thresholds are the stored hitmap values).

  energy     i_det at the centre and the edge site for f = 0.2 ... 1.  Measured detection current falls linearly with
             photon energy over a wide range (Renema et al. PRL 112, 117604 (2014); Engel et al. SuST 28, 114003
             (2015)).  In this model f scales the epsilon suppression, i.e. the integrated depairing, so a linear
             i_det(f) is the analogue.
  position   weak photons: the edge is the most sensitive line; strong photons: the profile across the width
             flattens (Zotova & Vodolazov SuST 27, 125001 (2014); Renema et al. Nano Lett. 15, 4541 (2015)).
  width      i_det ~ 1 - d_hs / w across 30-200 nm strips (stored results; voltage-trigger bound only: under the
             sustained criterion every narrow strip sits at its own floor, README "Two criteria as bounds").
  arrival    the photon at t0 = 40 instead of 20 tau0 (same window after it): the bias transient has settled.
  mirror     sites at +x against the stored -x: mesh asymmetry.
"""
from __future__ import annotations

import time
from dataclasses import replace

import numpy as np

from snspd import geometry as geo
from snspd import sim, store
from snspd.hotspot import HOTSPOT

from . import common as C

ID, TITLE = "1d", "plain-wire photon thresholds against known behaviour"
G = geo.PLAIN_WIRE
HITMAP = "results/hitmap/hitmap_conserving_T400_m1.4.csv"
ISW = "results/hitmap/isw_conserving_T400_m1.4.json"
TRIGGER_DIR = "results/v2_voltage_trigger"
SITES = {"centre": (0.0, 0.0), "edge": (-96.0, 0.0)}
ENERGIES = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
T0_LATE = 40.0
FIELDS = ["case", "x_nm", "y_nm", "energy_f", "t0", "tmax", "I_sw_uA", "I_det_uA", "i_det", "flag", "n_sims",
          "branch", "t_trigger", "wall_s"]


def csv_path():
    return C.out_path(f"photon_{C.CFG.tag}.csv")


def task_fn(task):
    t_start = time.time()
    try:
        params = replace(HOTSPOT, t0_tau0=task["t0"])
        cfg = replace(C.CFG, tmax=C.CFG.tmax + task["t0"] - HOTSPOT.t0_tau0)     # same window after the photon
        I, flag, n, trig = C.photon_threshold(G, task["x"], task["y"], task["f"], task["I_sw"], cfg, params)
        return dict(ok=True, task=task, res=dict(I=I, flag=flag, n=n, branch=trig.branch, t=trig.t, tmax=cfg.tmax),
                    wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def stored():
    """{(x, y, f): (i_det in I_sw, flag)} from the hitmap run, f = 1 grid and probe energies."""
    return {(float(r["x_nm"]), float(r["y_nm"]), round(float(r["energy_f"]), 3)): (float(r["i_det"]), r["flag"])
            for r in store.load_rows(HITMAP) if r["Geometry"] == G.key}


def run(a):
    I_sw = store.load_json(ISW, {})[G.key]["I_sw"]
    ref = stored()
    have = {(r["case"], float(r["x_nm"]), float(r["y_nm"]), round(float(r["energy_f"]), 3), float(r["t0"]))
            for r in store.load_rows(csv_path())}
    want = [("energy", *SITES[s], f, HOTSPOT.t0_tau0) for s in SITES for f in ENERGIES
            if (*SITES[s], round(f, 3)) not in ref]
    want += [("arrival", *SITES[s], 1.0, T0_LATE) for s in SITES]
    want += [("mirror", -x, 0.0, 1.0, HOTSPOT.t0_tau0) for x in (-58.0, -96.0)]
    tasks = [dict(case=c, x=x, y=y, f=f, t0=t0, I_sw=I_sw) for c, x, y, f, t0 in want
             if (c, x, y, round(f, 3), t0) not in have]
    sim.log(f"-- 1d: {len(tasks)} threshold searches to do ({len(want) - len(tasks)} done), I_sw = {I_sw:.2f} uA")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['case']} ({t['x']:g},{t['y']:g}) f={t['f']}: {r['error']}")
            return
        q = r["res"]
        store.append_row(csv_path(), dict(case=t["case"], x_nm=t["x"], y_nm=t["y"], energy_f=t["f"], t0=t["t0"],
                                          tmax=q["tmax"], I_sw_uA=I_sw, I_det_uA=q["I"], i_det=q["I"] / I_sw,
                                          flag=q["flag"], n_sims=q["n"], branch=q["branch"],
                                          t_trigger="" if q["t"] is None else round(q["t"], 2),
                                          wall_s=round(r["wall_s"], 1)), FIELDS)
        sim.log(f"   {t['case']:<8} ({t['x']:6.1f},{t['y']:4.1f}) f={t['f']:.2f} t0={t['t0']:g}: "
                f"i_det = {q['I'] / I_sw:.4f} {q['flag']}  ({r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result, fn=task_fn)
    return analyse(I_sw, ref)


def analyse(I_sw, ref):
    tol = sim.scales(G)[0] / I_sw                      # bisection tolerance in units of I_sw
    i_dep = sim.FILM.I_dep_uA(G.width_nm)
    rows = store.load_rows(csv_path())
    val = dict(ref)                                    # (x, y, f) -> (i_det, flag), t0 = 20
    for r in rows:
        if r["case"] in ("energy", "mirror"):
            val[(float(r["x_nm"]), float(r["y_nm"]), round(float(r["energy_f"]), 3))] = (float(r["i_det"]), r["flag"])
    late = {(float(r["x_nm"]), float(r["y_nm"])): float(r["i_det"]) for r in rows if r["case"] == "arrival"}
    checks, extra = [], {}

    # energy dependence
    for k, (name, (x, y)) in enumerate(SITES.items()):
        fs = sorted({f for (xx, yy, f) in val if (xx, yy) == (x, y)})
        pts = [(f, *val[(x, y, f)]) for f in fs]
        extra[f"energy_{name}"] = [(f, i, fl) for f, i, fl in pts]
        det = [(f, i) for f, i, fl in pts if fl != "at_wall"]
        if len(det) >= 4:
            F, Iv = np.array(det).T
            p = np.polyfit(F, Iv, 1)
            res = Iv - np.polyval(p, F)
            f_sat = min(f for f in F if all(abs(i - Iv[-1]) <= tol for ff, i in zip(F, Iv) if ff >= f))
            checks.append(C.check(
                f"1d.{k + 1}", f"{name} site ({x:g},{y:g}): linearity of i_det(f) over the detected range "
                f"f = {F.min():g}-{F.max():g}", f"max |residual| {np.abs(res).max():.4f}",
                f"<= 3 x bisection tol = {3 * tol:.4f} (linear in photon energy)", np.abs(res).max() <= 3 * tol,
                f"slope {p[0] * I_sw / i_dep:+.3f} I_dep per unit f; i_det stays within tol of its f = 1 value for "
                f"f >= {f_sat:g}; i_det(f) = "
                + ", ".join(f"{f:g}:{i:.3f}{'*' if fl == 'at_wall' else ''}" for f, i, fl in pts)
                + "   (* = no detection below I_sw)"))
            mono = all(b <= a_ + tol for a_, b in zip(Iv, Iv[1:]))
            checks.append(C.check(f"1d.{k + 3}", f"{name} site: i_det falls monotonically with f", "yes" if mono else "no",
                                  "monotonic within tol", mono))
        else:
            checks.append(C.check(f"1d.{k + 1}", f"{name} site: linearity of i_det(f)", f"{len(det)} detected points",
                                  ">= 4 detected energies needed", False))

    # position dependence: edge (-96) vs centre at each energy
    spread = {f: val[(0.0, 0.0, f)][0] - val[(-96.0, 0.0, f)][0] for f in sorted({f for (_, _, f) in val})
              if (0.0, 0.0, f) in val and (-96.0, 0.0, f) in val}
    extra["spread_centre_minus_edge"] = spread
    weak = [f for f in spread if f <= 0.6 and val[(-96.0, 0.0, f)][1] != "at_wall"]
    ok = all(spread[f] > tol for f in weak)
    checks.append(C.check("1d.5", "weak photons (f <= 0.6, where the edge detects): the edge is more sensitive than "
                          "the centre", ", ".join(f"f={f:g}: centre-edge {spread[f]:+.3f}" for f in weak),
                          "centre - edge > tol at every such f", ok,
                          "f = 0.2: neither site detects below I_sw" if 0.2 in spread else ""))
    if 0.6 in spread and 1.0 in spread:
        checks.append(C.check("1d.6", "strong photons: the profile across the width flattens",
                              f"|centre - edge| {abs(spread[1.0]):.4f} at f = 1 vs {abs(spread[0.6]):.4f} at f = 0.6",
                              "smaller at f = 1", abs(spread[1.0]) < abs(spread[0.6]),
                              "all f: " + ", ".join(f"{f:g}:{s:+.3f}" for f, s in spread.items())))

    # width law, voltage-trigger bound (stored)
    widths, ivals = [], []
    for w in (30, 56, 100):
        rr = [r for r in store.load_rows(f"{TRIGGER_DIR}/narrow_baseline/narrow{w}nm_hitmap_conserving_T400_m1.4.csv")
              if float(r["x_nm"]) == 0.0 and float(r["energy_f"]) == 1.0]
        if rr:
            widths.append(w), ivals.append(float(rr[0]["i_det"]))
    rr = [r for r in store.load_rows(f"{TRIGGER_DIR}/hitmap/hitmap_conserving_T400_m1.4.csv")
          if r["Geometry"] == G.key and float(r["x_nm"]) == 0.0 and float(r["energy_f"]) == 1.0]
    if rr:
        widths.append(200), ivals.append(float(rr[0]["i_det"]))
    if len(widths) >= 3 and 200 in widths:
        w, iv = np.array(widths, float), np.array(ivals)
        nar = w <= 100
        d_hs = float(np.sum((1 - iv[nar]) / w[nar]) / np.sum(1 / w[nar] ** 2))
        res = iv - (1 - d_hs / w)
        checks.append(C.check("1d.7", "width law i_det = 1 - d_hs/w at the centre, 30-100 nm strips, voltage-trigger "
                              "bound", f"d_hs = {d_hs:.1f} nm, max |residual| {np.abs(res[nar]).max():.4f}",
                              f"max |residual| <= 3 x tol = {3 * tol:.4f}", np.abs(res[nar]).max() <= 3 * tol,
                              "(w, i_det, (1 - i) w): " + ", ".join(f"({a_:g}, {b:.3f}, {(1 - b) * a_:.1f} nm)"
                                                                   for a_, b in zip(w, iv))))
        checks.append(C.check("1d.8", "the 200 nm wire against that law", f"{iv[~nar][0]:.3f} vs {1 - d_hs / 200:.3f} "
                              f"({res[~nar][0]:+.3f})", "on the law if detection is a normal belt of width d_hs",
                              note="(1 - i) w = 26 nm instead of 20: the wide wire detects at lower bias than a d_hs "
                              "belt would allow"))

    # arrival time and mirror symmetry
    for k, (name, (x, y)) in enumerate(SITES.items()):
        if (x, y) in late:
            d = late[(x, y)] - ref[(x, y, 1.0)][0]
            checks.append(C.check(f"1d.{9 + k}", f"{name} site: photon at t0 = {T0_LATE:g} vs 20 tau0",
                                  f"{late[(x, y)]:.4f} vs {ref[(x, y, 1.0)][0]:.4f} ({d:+.4f})",
                                  f"|diff| <= tol = {tol:.4f}", abs(d) <= 1.01 * tol))
    for k, x in enumerate((58.0, 96.0)):
        if (x, 0.0, 1.0) in val:
            d = val[(x, 0.0, 1.0)][0] - ref[(-x, 0.0, 1.0)][0]
            checks.append(C.check(f"1d.{11 + k}", f"mirror: ({x:g},0) vs ({-x:g},0)",
                                  f"{val[(x, 0.0, 1.0)][0]:.4f} vs {ref[(-x, 0.0, 1.0)][0]:.4f} ({d:+.4f})",
                                  f"|diff| <= tol = {tol:.4f}", abs(d) <= 1.01 * tol))
    return checks, extra
