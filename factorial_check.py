#!/usr/bin/env python3
"""
factorial_check.py -- what gamma and the NbTiN hotspot times each do to the photon thresholds, at the 14 sites of
thermal_check.py, under both criteria (voltage trigger = lower bound, sustained = upper bound).

The "nbtin" parameter set (snspd.sim.PRESETS, run_nbtin.sh) changes two things at once: gamma 1 -> 10 and the
hotspot times (D 10 -> 7.1 nm^2/tau0, escape 100 -> 530 tau0).  The four corners of that 2 x 2:

  base         gamma 1,  NbN times    the stored production results (re-run here for the per-run details)
  nbtin-times  gamma 1,  NbTiN times
  g10          gamma 10, NbN times    photon into the settled state
  nbtin        gamma 10, NbTiN times  photon into the settled state (also in the full maps of run_nbtin.sh)

Same model as the maps: kinetic hotspot f = 1, 400 tau0, 1.4 xi, pinned meshes, bisection in [0.30, 1] x I_sw to
sim.TOL_REL.  The dark I_sw does not depend on the hotspot, so gamma 1 corners use the stored production values and
gamma 10 corners the gamma 10 values of run_nbtin.sh (results/*/..._kinetic-nbtin_...json) or results/gamma/isw.csv;
missing ones are bisected here first.  Every run is stored (results/factorial/runs.csv) with the phase it slipped
after t0, so the highest run below threshold shows how large a pulse that healed was.

  python factorial_check.py run --corners base,nbtin-times,g10 --workers 3      # resumable
  python factorial_check.py report

Writes results/factorial/check.csv, runs.csv, isw.csv, summary.txt, check.png.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")

from snspd import sim, store
from thermal_check import CASES, geometry

OUT = "results/factorial"
CRITERIA = ("trigger", "sustained")
CORNERS = {"base": "", "nbtin-times": "nbtin-times", "g10": "g10", "nbtin": "nbtin"}     # corner -> sim preset
CHECK_FIELDS = ["corner", "label", "geom", "x_nm", "y_nm", "criterion", "I_det_uA", "i_det", "flag", "n_runs", "branch",
                "t_trigger", "slips_at", "quiet_bias", "quiet_branch", "quiet_slips", "I_sw_uA", "wall_s"]
RUN_FIELDS = ["corner", "label", "geom", "criterion", "bias", "switched", "branch", "t_trigger", "slips", "slips_last",
              "slips_settled", "settled_start"]
ISW_FIELDS = ["gamma", "geom", "criterion", "I_sw_uA", "i_sw", "flag", "n_runs", "wall_s"]


def gkey(g):
    return f"{g} nm" if isinstance(g, (int, float)) else g


def config(corner, criterion):
    return sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", sustained=criterion == "sustained", out_dir=OUT,
                         preset=CORNERS[corner])


def isw_json(g, preset, criterion):
    """Path and key of the production dark I_sw of geometry g (thermal_check CASES key) for a preset."""
    suf = "_trigger" if criterion == "trigger" else ""
    model = "kinetic" + (f"-{preset}" if preset else "")
    if isinstance(g, (int, float)):
        return f"results/narrow_baseline/narrow{g:g}nm_isw_{model}_T400_m1.4{suf}.json", "Narrow_Wire"
    d = "results/slots" if g.startswith("Slot") or g == "Plain_Wire_2um" else "results/hitmap"
    return f"{d}/isw_{model}_T400_m1.4{suf}.json", g


def known_isw(g, gamma, criterion):
    """Dark I_sw [uA] at this gamma from stored results, or None."""
    for preset in (("",) if gamma == 1.0 else ("nbtin", "g10")):
        path, key = isw_json(g, preset, criterion)
        v = store.load_json(path, {}).get(key)
        if v:
            return float(v["I_sw"])
    for src in ((os.path.join(OUT, "isw.csv"), "gamma"), ("results/gamma/isw.csv", "gamma")):
        for r in store.load_rows(src[0]):
            if float(r[src[1]]) == gamma and r["geom"] == gkey(g) and r["criterion"] == criterion:
                return float(r["I_sw_uA"])
    return None


def run_once(ctx, current_uA, eps, cfg):
    """One run: the criterion's verdict, the phase slipped between the probes after t0, and whether it started from
    the settled state (sim.settled_state; only for presets that ask for it, photon runs only)."""
    t0 = cfg.pset.params.t0_tau0
    init = sim.settled_state(ctx["device"], current_uA, cfg) if (eps is not None and cfg.pset.settled_start) else None
    with sim.solved(ctx["device"], current_uA, eps, cfg, init=init) as sol:
        t, vs = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
        slips = float(sim.phase_slips(t, vs, t0)[-1])
        trig = sim.detect(sol, t0)
        trig = (sim.settle(ctx["device"], current_uA, sol, trig, cfg) if cfg.sustained
                else trig._replace(switched=bool(trig.branch)))
    return trig, slips, init is not None


def task_fn(task):
    """kind 'isw': dark switching current; kind 'hit': photon threshold at (x, y).  Never raises."""
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        cfg = config(task["corner"], task["criterion"])
        film, params, g = cfg.pset.film, cfg.pset.params, geometry(task["geom"])
        ctx = sim.context(g, cfg.mesh_xi, film, params)
        tol, (lo, hi) = sim.scales(g, film)
        eps = None
        if task["kind"] == "hit":
            lo, hi = cfg.hit_lo_frac * task["I_sw_uA"], task["I_sw_uA"]
            eps = hsm.build_eps(ctx["hs"], "kinetic", task["x"], task["y"], 1.0, cfg.tmax, params)
        lo0, hi0, runs, at, quiet = lo, hi, [], None, None
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            trig, slips, settled = run_once(ctx, mid, eps, cfg)
            r = dict(bias_uA=mid, switched=trig.switched, branch=trig.branch, t=trig.t, slips=slips,
                     slips_last=trig.slips_last, slips_settled=trig.slips_settled, settled=settled)
            runs.append(r)
            if trig.switched:
                hi, at = mid, r
            else:
                lo, quiet = mid, r                  # lo only rises: the last one is the highest non-switching run
        flag = "at_wall" if hi >= hi0 else ("below_bracket" if lo <= lo0 else "")
        return dict(ok=True, task=task, I=hi, flag=flag, runs=runs, at=at, quiet=quiet,
                    i_dep=film.I_dep_uA(g.free_width_min_nm()), wall_s=time.time() - t_start)
    except Exception as exc:
        import traceback
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
                    wall_s=time.time() - t_start)


def _f(v, nd):
    return "" if v is None else round(float(v), nd)


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    corners = [c for c in a.corners.split(",") if c]
    done = {(r["corner"], r["label"], r["criterion"]) for r in store.load_rows(os.path.join(OUT, "check.csv"))}
    first, waiting = [], {}                  # waiting: (gamma, geometry, criterion) -> photon tasks needing its I_sw
    for corner in corners:
        gamma = config(corner, "trigger").pset.film.gamma
        for label, g, x, y, _, _ in CASES:
            for crit in CRITERIA:
                if (corner, label, crit) in done:
                    continue
                t = dict(kind="hit", corner=corner, label=label, geom=g, x=float(x), y=float(y), criterion=crit)
                isw = known_isw(g, gamma, crit)
                if isw is not None:
                    first.append(dict(t, I_sw_uA=isw))
                else:
                    key = (gamma, gkey(g), crit)
                    if key not in waiting:
                        first.append(dict(kind="isw", corner=corner, label="dark", geom=g, criterion=crit))
                    waiting.setdefault(key, []).append(t)
    first.sort(key=lambda t: t["kind"] != "isw")
    sim.log(f"-- factorial check {corners}: {sum(t['kind'] == 'isw' for t in first)} dark I_sw, "
            f"{sum(t['kind'] == 'hit' for t in first) + sum(len(v) for v in waiting.values())} photon thresholds")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['corner']} {t['label']} {t['criterion']}: {r['error']}")
            return []
        for q in r["runs"]:
            store.append_row(os.path.join(OUT, "runs.csv"), dict(
                corner=t["corner"], label=t["label"], geom=gkey(t["geom"]), criterion=t["criterion"],
                bias=round(q["bias_uA"] / r["i_dep"], 4), switched=int(q["switched"]), branch=q["branch"],
                t_trigger=_f(q["t"], 1), slips=round(q["slips"], 2), slips_last=round(q["slips_last"], 2),
                slips_settled=_f(q["slips_settled"], 2), settled_start=int(q["settled"])), RUN_FIELDS)
        at, qu, gk = r["at"] or {}, r["quiet"] or {}, gkey(t["geom"])
        if t["kind"] == "isw":
            gamma = config(t["corner"], t["criterion"]).pset.film.gamma
            store.append_row(os.path.join(OUT, "isw.csv"), dict(
                gamma=gamma, geom=gk, criterion=t["criterion"], I_sw_uA=round(r["I"], 3),
                i_sw=round(r["I"] / r["i_dep"], 4), flag="above_bracket" if r["flag"] == "at_wall" else r["flag"],
                n_runs=len(r["runs"]), wall_s=round(r["wall_s"])), ISW_FIELDS)
            sim.log(f"   dark I_sw gamma {gamma:g} {gk} {t['criterion']}: {r['I'] / r['i_dep']:.3f} I_dep")
            return [dict(s, I_sw_uA=r["I"]) for s in waiting.pop((gamma, gk, t["criterion"]), [])]
        store.append_row(os.path.join(OUT, "check.csv"), dict(
            corner=t["corner"], label=t["label"], geom=gk, x_nm=t["x"], y_nm=t["y"], criterion=t["criterion"],
            I_det_uA=round(r["I"], 3), i_det=round(r["I"] / r["i_dep"], 4), flag=r["flag"], n_runs=len(r["runs"]),
            branch=at.get("branch", ""), t_trigger=_f(at.get("t"), 1), slips_at=_f(at.get("slips"), 2),
            quiet_bias="" if not qu else round(qu["bias_uA"] / r["i_dep"], 4), quiet_branch=qu.get("branch", ""),
            quiet_slips=_f(qu.get("slips"), 2), I_sw_uA=round(t["I_sw_uA"], 3), wall_s=round(r["wall_s"])), CHECK_FIELDS)
        sim.log(f"   {t['corner']:<12} {t['label']:<20} {t['criterion']:<9} {r['I'] / r['i_dep']:.3f} I_dep "
                f"{r['flag']}  {r['wall_s']:.0f} s")
        return []

    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    n_done = 0
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn")) as pool:
        pending = {pool.submit(task_fn, t) for t in first}
        while pending:
            fin, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in fin:
                for nt in on_result(fut.result()):
                    pending.add(pool.submit(task_fn, nt))
                n_done += 1
                sim.log(f"    [{n_done} done, {len(pending)} running or queued]")
    cmd_report(a)


def map_value(g, x, y, criterion):
    """The nbtin threshold [I_dep] of a site from the full maps of run_nbtin.sh, or None."""
    suf = "_trigger" if criterion == "trigger" else ""
    if isinstance(g, (int, float)):
        path, key = f"results/narrow_baseline/narrow{g:g}nm_hitmap_kinetic-nbtin_T400_m1.4{suf}.csv", "Narrow_Wire"
    else:
        d = "results/slots" if g.startswith("Slot") else "results/hitmap"
        path, key = f"{d}/hitmap_kinetic-nbtin_T400_m1.4{suf}.csv", g
    for r in store.load_rows(path):
        if (r["Geometry"] == key and float(r["x_nm"]) == x and float(r["y_nm"]) == y and float(r["energy_f"]) == 1.0):
            return float(r["I_det_uA"]) / sim.FILM.I_dep_uA(geometry(g).free_width_min_nm())
    return None


def table(rows):
    """{(corner, label, criterion): i_det}, the nbtin corner filled from the maps where this script has no row."""
    val = {(r["corner"], r["label"], r["criterion"]): float(r["i_det"]) for r in rows}
    for label, g, x, y, lo1, hi1 in CASES:
        for crit, ref in zip(CRITERIA, (lo1, hi1)):
            val.setdefault(("base", label, crit), ref)
            if ("nbtin", label, crit) not in val:
                m = map_value(g, x, y, crit)
                if m is not None:
                    val[("nbtin", label, crit)] = m
    return val


def cmd_report(a):
    rows = store.load_rows(os.path.join(OUT, "check.csv"))
    val = table(rows)
    by = {(r["corner"], r["label"], r["criterion"]): r for r in rows}
    fm = lambda v: f"{v:7.3f}" if v is not None else f"{'-':>7}"
    sg = lambda v: f"{v:+7.3f}" if v is not None else f"{'-':>7}"
    L = ["FACTORIAL CHECK: gamma (1 -> 10) x hotspot times (NbN -> NbTiN), kinetic hotspot f = 1, 400 tau0, 1.4 xi",
         "thresholds in I_dep of the narrowest free width; corners: base (gamma 1, NbN times), times (gamma 1, NbTiN "
         "times), g10 (gamma 10, NbN times), nbtin (gamma 10, NbTiN times)",
         "effects: d_times = times - base, d_gamma = g10 - base, both = nbtin - base, interaction = nbtin - times - g10 "
         "+ base", ""]
    for crit in CRITERIA:
        L += [f"--- {crit} ---",
              f"{'site':<20}{'base':>7}{'times':>7}{'g10':>7}{'nbtin':>7}   {'d_times':>7}{'d_gamma':>8}{'both':>7}"
              f"{'interact':>9}"]
        for label, *_ in CASES:
            b, t, g, n = (val.get((c, label, crit)) for c in CORNERS)
            dt = None if t is None or b is None else t - b
            dg = None if g is None or b is None else g - b
            dn = None if n is None or b is None else n - b
            ia = None if None in (b, t, g, n) else n - t - g + b
            L.append(f"{label:<20}{fm(b)}{fm(t)}{fm(g)}{fm(n)}   {sg(dt)}{sg(dg)} {sg(dn)}{sg(ia)}")
        L.append("")
    L += ["--- gap between the bounds (sustained - trigger) ---", f"{'site':<20}" + "".join(f"{c:>13}" for c in CORNERS)]
    for label, *_ in CASES:
        cells = []
        for c in CORNERS:
            lo, hi = val.get((c, label, "trigger")), val.get((c, label, "sustained"))
            cells.append(f"{hi - lo:13.3f}" if lo is not None and hi is not None else f"{'-':>13}")
        L.append(f"{label:<20}" + "".join(cells))
    L += ["", "--- run at threshold: trigger time after t0 [tau0] and phase slipped [2 pi]; highest run below threshold: "
          "bias, branch ('-' = no voltage trigger), slips ---"]
    for c in CORNERS:
        for label, *_ in CASES:
            for crit in CRITERIA:
                r = by.get((c, label, crit))
                if not r:
                    continue
                tt = f"{float(r['t_trigger']) - 20:6.1f}" if r["t_trigger"] else f"{'-':>6}"
                q = (f"{float(r['quiet_bias']):.3f} {r['quiet_branch'] or '-':<7} {r['quiet_slips']:>7}"
                     if r["quiet_bias"] else "(none)")
                L.append(f"  {c:<12}{label:<20}{crit:<10} at {fm(float(r['i_det']))} t {tt} slips {r['slips_at']:>7}"
                         f"   below: {q}" + (f"   [{r['flag']}]" if r["flag"] else ""))
    text = "\n".join(L)
    sim.log(text)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    figure(val)


def figure(val):
    """One row per site, one bar per corner from the trigger to the sustained threshold."""
    import matplotlib.pyplot as plt
    from snspd import plotting as pl
    pl.apply_style()
    labels = [c[0] for c in CASES]
    fig, ax = plt.subplots(figsize=(9.0, 0.78 * len(labels) + 1.6), layout="constrained")
    style = {"base": ("#b5b3ac", "gamma 1, NbN times (stored)"), "nbtin-times": ("#eda100", "gamma 1, NbTiN times"),
             "g10": ("#9fc4ee", "gamma 10, NbN times"), "nbtin": ("#2a78d6", "gamma 10, NbTiN times (nbtin)")}
    offs = dict(zip(CORNERS, (0.27, 0.09, -0.09, -0.27)))
    for k, label in enumerate(labels):
        yy = len(labels) - 1 - k
        for c in CORNERS:
            lo, hi = val.get((c, label, "trigger")), val.get((c, label, "sustained"))
            if lo is None or hi is None:
                continue
            ax.plot([lo, max(hi, lo + 0.003)], [yy + offs[c]] * 2, color=style[c][0], lw=4.5, solid_capstyle="butt")
    ax.set_yticks(range(len(labels)), labels[::-1], fontsize=9)
    ax.set_xlabel("photon threshold  I / I_dep   (bar: voltage trigger to sustained; kinetic hotspot f = 1, 400 tau0)")
    ax.grid(axis="y", visible=False)
    fig.legend(handles=[plt.Line2D([], [], color=s[0], lw=4.5, label=s[1]) for s in style.values()],
               loc="outside lower center", ncol=2, frameon=False, fontsize=8.5)
    ax.set_title("What gamma and the NbTiN hotspot times each do to the two bounds", fontsize=11)
    fig.savefig(os.path.join(OUT, "check.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--corners", default="base,nbtin-times,g10", help=f"comma-separated, from {list(CORNERS)}")
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
