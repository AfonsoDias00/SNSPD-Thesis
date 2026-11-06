#!/usr/bin/env python3
"""
gamma_check.py -- does the GTDGL inelastic parameter gamma change the photon thresholds and the gap between the two
bounds?  Every result so far uses gamma = 1 (snspd/units.py); pyTDGL's default, and Simon et al. (PRB 112, 174512
(2025)), use gamma = 10, and gamma = 2 tau_E Delta0 / hbar with the electron-phonon time measured on 9 nm NbTiN
(tau_EP(T_c) 2.9-4.0 ps, Sidorova et al.) gives ~12-16.

gamma sets how fast |psi| follows a change: in the GTDGL the amplitude relaxes with u sqrt(1 + gamma^2 |psi|^2), i.e.
~8 tau0 at gamma = 1 and ~58 tau0 at gamma = 10 for |psi| ~ 1.  A larger gamma makes a phase-slip region heal more
slowly.

At the 14 sites of thermal_check.py: kinetic hotspot f = 1, 400 tau0, 1.4 xi, pinned meshes, no heat, both criteria
(voltage trigger = lower bound, sustained = upper bound), exactly as the stored gamma = 1 results, with Film(gamma).
The dark I_sw is bisected again at the new gamma first (bracket 0.70-1.20 I_dep), and each photon threshold is
bisected in [0.30, 1] x that I_sw, as in the production runs.  The voltage thresholds of the criteria are unchanged.
Every run is stored with the phase it slipped after t0 (slips, in 2 pi), so a run below threshold that triggered and
healed shows how large its pulse was.

Two protocols.  "fresh" is the production protocol: the bias is switched on at t = 0 on the uniform state and the
photon arrives at t0 = 20 tau0.  At gamma = 1 the strip has settled by then (plain wire at 0.80 I_dep: <|psi|^2> 0.880
at t0, 0.874 settled); at gamma = 10 it has not (0.953 at t0, 0.912 at 50 tau0, 0.875 at 200 tau0), so the photon
meets a stronger condensate.  "settled": every photon run starts from the final state of a DARK_SETTLE tau0 dark run
at the same bias (pyTDGL seed; the clock restarts at 0, so the photon still arrives at t0).  gamma = 1 runs use the
stored dark I_sw; with --protocols settled they are the control for the seeding itself.

  python gamma_check.py run --workers 11 [--gamma 10] [--protocols fresh,settled]     # resumable
  python gamma_check.py run --workers 3 --gamma 1 --protocols settled --sites control   # or "30 nm (0,0);plain (0,0)"
  python gamma_check.py report

Writes results/gamma/isw.csv (dark I_sw), check.csv (one row per site, criterion and gamma), runs.csv (every run),
summary.txt and check.png.
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
from snspd.units import Film
from thermal_check import CASES, geometry

OUT = "results/gamma"
CRITERIA = ("trigger", "sustained")
PROTOCOLS = ("fresh", "settled")
DARK_SETTLE = 250.0                      # dark run before each "settled" photon run [tau0]
CONTROL_SITES = ("30 nm (0,0)", "plain (0,0)", "zipper (0,0)", "zipper (-68,0)", "Slot 400 (-58,0)")
ISW_FIELDS = ["gamma", "geom", "criterion", "I_sw_uA", "i_sw", "flag", "n_runs", "branch", "t_trigger", "wall_s"]
CHECK_FIELDS = ["gamma", "protocol", "label", "geom", "x_nm", "y_nm", "criterion", "I_det_uA", "i_det", "flag", "n_runs", "branch",
                "t_trigger", "slips_at", "quiet_bias", "quiet_branch", "quiet_slips", "I_sw_uA", "wall_s"]
RUN_FIELDS = ["gamma", "protocol", "label", "geom", "criterion", "bias", "switched", "branch", "t_trigger", "slips", "slips_last",
              "slips_settled"]

# gamma = 1 dark I_sw of the stored kinetic runs, (directory, file stem prefix, key) per geometry
ISW_G1 = {30: ("results/narrow_baseline", "narrow30nm_", "Narrow_Wire"),
          56: ("results/narrow_baseline", "narrow56nm_", "Narrow_Wire"),
          100: ("results/narrow_baseline", "narrow100nm_", "Narrow_Wire"),
          "Plain_Wire": ("results/hitmap", "", "Plain_Wire"), "Long_Zipper": ("results/hitmap", "", "Long_Zipper"),
          "Slot_400": ("results/slots", "", "Slot_400")}


def gkey(g):
    return f"{g} nm" if isinstance(g, (int, float)) else g


def config(criterion):
    return sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", sustained=criterion == "sustained", out_dir=OUT)


def run_once(ctx, current_uA, eps, cfg, settled=False):
    """One run: the criterion's verdict plus the phase slipped between the probes after t0.  settled: start from the
    end of a DARK_SETTLE tau0 dark run at the same current instead of the uniform state."""
    from contextlib import nullcontext
    from dataclasses import replace
    from snspd import hotspot as hsm
    t0 = hsm.HOTSPOT.t0_tau0
    dark = (sim.solved(ctx["device"], current_uA, None, replace(cfg, tmax=DARK_SETTLE)) if settled
            else nullcontext(None))
    with dark as seed, sim.solved(ctx["device"], current_uA, eps, cfg, seed=seed) as sol:
        t, vs = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
        slips = float(sim.phase_slips(t, vs, t0)[-1])
        trig = sim.detect(sol, t0)
        trig = (sim.settle(ctx["device"], current_uA, sol, trig, cfg) if cfg.sustained
                else trig._replace(switched=bool(trig.branch)))
    return trig, slips


def task_fn(task):
    """kind 'isw': dark switching current; kind 'hit': photon threshold at (x, y).  Never raises."""
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        film, g, cfg = Film(gamma=task["gamma"]), geometry(task["geom"]), config(task["criterion"])
        ctx = sim.context(g, cfg.mesh_xi, film)
        tol, (lo, hi) = sim.scales(g, film)
        eps = None
        if task["kind"] == "hit":
            lo, hi = cfg.hit_lo_frac * task["I_sw_uA"], task["I_sw_uA"]
            eps = hsm.build_eps(ctx["hs"], "kinetic", task["x"], task["y"], 1.0, cfg.tmax)
        lo0, hi0, runs, at, quiet = lo, hi, [], None, None
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            trig, slips = run_once(ctx, mid, eps, cfg, settled=task.get("protocol") == "settled")
            r = dict(bias_uA=mid, switched=trig.switched, branch=trig.branch, t=trig.t, slips=slips,
                     slips_last=trig.slips_last, slips_settled=trig.slips_settled)
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
    gam = float(a.gamma)
    protocols = [p for p in a.protocols.split(",") if p]
    labels = (list(CONTROL_SITES) if a.sites == "control" else
              [v.strip() for v in a.sites.split(";") if v.strip()] if a.sites else [c[0] for c in CASES])
    if gam == 1.0:          # the stored dark I_sw of the production runs
        isw = {(gkey(c[1]), cr): isw_gamma1(c[1], cr) * Film().I_dep_uA(geometry(c[1]).free_width_min_nm())
               for c in CASES for cr in CRITERIA}
    else:
        isw = {(r["geom"], r["criterion"]): float(r["I_sw_uA"])
               for r in store.load_rows(os.path.join(OUT, "isw.csv")) if float(r["gamma"]) == gam}
    done = {(r["protocol"], r["label"], r["criterion"]) for r in store.load_rows(os.path.join(OUT, "check.csv"))
            if float(r["gamma"]) == gam}
    sites = {}                               # (geometry, criterion) -> photon tasks still to run
    for label, g, x, y, _, _ in CASES:
        for crit in CRITERIA:
            for prot in protocols:
                if label in labels and (prot, label, crit) not in done:
                    sites.setdefault((gkey(g), crit), []).append(
                        dict(kind="hit", gamma=gam, protocol=prot, label=label, geom=g, x=float(x), y=float(y),
                             criterion=crit))
    geoms = {gkey(c[1]): c[1] for c in CASES}
    first = []
    for (gk, crit), ts in sites.items():
        if (gk, crit) in isw:
            first += [dict(t, I_sw_uA=isw[(gk, crit)]) for t in ts]
        else:
            first.append(dict(kind="isw", gamma=gam, protocol="fresh", label="dark", geom=geoms[gk], criterion=crit))
    first.sort(key=lambda t: t["kind"] != "isw")
    sim.log(f"-- gamma = {gam:g}: {sum(t['kind'] == 'isw' for t in first)} dark I_sw, "
            f"{sum(len(v) for v in sites.values())} photon thresholds")

    def store_runs(t, runs, i_dep):
        for q in runs:
            store.append_row(os.path.join(OUT, "runs.csv"), dict(
                gamma=gam, protocol=t["protocol"], label=t["label"], geom=gkey(t["geom"]), criterion=t["criterion"],
                bias=round(q["bias_uA"] / i_dep, 4), switched=int(q["switched"]), branch=q["branch"],
                t_trigger=_f(q["t"], 1), slips=round(q["slips"], 2), slips_last=round(q["slips_last"], 2),
                slips_settled=_f(q["slips_settled"], 2)), RUN_FIELDS)

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['label']} {gkey(t['geom'])} {t['criterion']}: {r['error']}")
            return []
        store_runs(t, r["runs"], r["i_dep"])
        at, qu, gk = r["at"] or {}, r["quiet"] or {}, gkey(t["geom"])
        if t["kind"] == "isw":
            flag = "above_bracket" if r["flag"] == "at_wall" else r["flag"]
            store.append_row(os.path.join(OUT, "isw.csv"), dict(
                gamma=gam, geom=gk, criterion=t["criterion"], I_sw_uA=round(r["I"], 3), i_sw=round(r["I"] / r["i_dep"], 4),
                flag=flag, n_runs=len(r["runs"]), branch=at.get("branch", ""), t_trigger=_f(at.get("t"), 1),
                wall_s=round(r["wall_s"])), ISW_FIELDS)
            sim.log(f"   dark I_sw {gk:<12} {t['criterion']:<9} {r['I'] / r['i_dep']:.3f} I_dep {flag}  "
                    f"{r['wall_s']:.0f} s")
            return [dict(s, I_sw_uA=r["I"]) for s in sites.get((gk, t["criterion"]), [])]
        store.append_row(os.path.join(OUT, "check.csv"), dict(
            gamma=gam, protocol=t["protocol"], label=t["label"], geom=gk, x_nm=t["x"], y_nm=t["y"], criterion=t["criterion"],
            I_det_uA=round(r["I"], 3), i_det=round(r["I"] / r["i_dep"], 4), flag=r["flag"], n_runs=len(r["runs"]),
            branch=at.get("branch", ""), t_trigger=_f(at.get("t"), 1), slips_at=_f(at.get("slips"), 2),
            quiet_bias="" if not qu else round(qu["bias_uA"] / r["i_dep"], 4), quiet_branch=qu.get("branch", ""),
            quiet_slips=_f(qu.get("slips"), 2), I_sw_uA=round(t["I_sw_uA"], 3), wall_s=round(r["wall_s"])), CHECK_FIELDS)
        sim.log(f"   {t['label']:<20} {t['criterion']:<9} {t['protocol']:<7} {r['I'] / r['i_dep']:.3f} I_dep {r['flag']}  "
                f"{r['wall_s']:.0f} s")
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


def isw_gamma1(g, crit):
    d, prefix, key = ISW_G1[g]
    suf = "_trigger" if crit == "trigger" else ""
    v = store.load_json(f"{d}/{prefix}isw_kinetic_T400_m1.4{suf}.json", {}).get(key)
    return None if v is None else v["I_sw"] / Film().I_dep_uA(geometry(g).free_width_min_nm())


def heat_values():
    """{label: (one temperature eta 0.59, NbTiN two temperatures)} from thermal_check.py."""
    rows = {}
    for r in store.load_rows("results/thermal/check.csv"):
        rows[(r.get("variant") or f"1T eta {float(r['eta']):.2f}", r["label"])] = float(r["i_thermal"])
    return {c[0]: (rows.get(("1T eta 0.59", c[0])), rows.get(("NbTiN 2T", c[0]))) for c in CASES}


def cmd_report(a):
    gam = float(a.gamma)
    rows = {(float(r["gamma"]), r["protocol"], r["label"], r["criterion"]): r
            for r in store.load_rows(os.path.join(OUT, "check.csv"))}
    isw = {(r["geom"], r["criterion"]): r for r in store.load_rows(os.path.join(OUT, "isw.csv"))
           if float(r["gamma"]) == gam}
    heat = heat_values()
    fmt = lambda v, w=7: f"{v:{w}.3f}" if v is not None else f"{'-':>{w}}"

    def val(g_, prot, label, crit):
        r = rows.get((g_, prot, label, crit))
        return None if r is None else float(r["i_det"])

    L = [f"GAMMA CHECK  gamma = 1 vs {gam:g}: kinetic hotspot f = 1, 400 tau0, 1.4 xi, pinned meshes, no heat "
         "(heat: thermal_check.py, gamma = 1)",
         "thresholds in I_dep (narrowest free width); trig = voltage trigger (lower bound), sust = sustained (upper)",
         "fresh = production protocol (bias on at t = 0, photon at 20 tau0); settled = photon into the state reached "
         f"after a {DARK_SETTLE:g} tau0 dark run", "",
         f"dark I_sw / I_dep (fresh)  gamma 1: trigger  sustained    gamma {gam:g}: trigger  sustained"]
    for g in dict.fromkeys(c[1] for c in CASES):
        v = [isw_gamma1(g, c) for c in CRITERIA]
        w = [float(isw[(gkey(g), c)]["i_sw"]) if (gkey(g), c) in isw else None for c in CRITERIA]
        fl = " ".join(isw[(gkey(g), c)]["flag"] for c in CRITERIA if (gkey(g), c) in isw and isw[(gkey(g), c)]["flag"])
        L.append(f"  {gkey(g):<24}{fmt(v[0], 16)}{fmt(v[1], 11)}{fmt(w[0], 20)}{fmt(w[1], 11)}  {fl}")
    hdr = (f"{'site':<20}{'g1 stored':>16}{'g1 settled':>16}{f'g{gam:g} fresh':>16}{f'g{gam:g} settled':>16}"
           f"{'gap g1':>8}{f'gap g{gam:g}':>9}{'heat 1T':>9}")
    L += ["", "photon thresholds, each column: trig  sust", hdr]
    repro = []
    for label, g, x, y, lo1, hi1 in CASES:
        c1 = [val(1.0, "settled", label, c) for c in CRITERIA]
        cf = [val(gam, "fresh", label, c) for c in CRITERIA]
        cs = [val(gam, "settled", label, c) for c in CRITERIA]
        gap = f"{cs[1] - cs[0]:9.3f}" if None not in cs else f"{'-':>9}"
        flags = " ".join(f"[{p} {c} {r['flag']}]" for (g_, p, lab, c), r in rows.items()
                         if lab == label and r["flag"] and g_ in (1.0, gam))
        L.append(f"{label:<20}{lo1:9.3f}{hi1:7.3f}{fmt(c1[0], 9)}{fmt(c1[1])}{fmt(cf[0], 9)}{fmt(cf[1])}"
                 f"{fmt(cs[0], 9)}{fmt(cs[1])}{hi1 - lo1:8.3f}{gap}{fmt(heat[label][0], 9)}  {flags}")
        for c, ref in zip(CRITERIA, (lo1, hi1)):
            r = val(1.0, "fresh", label, c)
            if r is not None:
                repro.append(f"{label} {c}: {r:.4f} vs stored {ref:.3f}")
    if repro:
        L += ["", "reproduction check, gamma 1 fresh (this script) vs the stored production thresholds:"]
        L += ["  " + q for q in repro]
    L += ["", "runs just below threshold (the highest bias that did not count) and the run at threshold:",
          "  bias, the voltage branch that fired ('-' = none), phase slipped after t0 [2 pi] | at threshold: trigger "
          "time after t0 [tau0], slips"]
    for g_, prot in ((1.0, "settled"), (gam, "settled"), (gam, "fresh")):
        for label, *_ in CASES:
            for c in CRITERIA:
                r = rows.get((g_, prot, label, c))
                if not r:
                    continue
                q = (f"{float(r['quiet_bias']):.3f} {r['quiet_branch'] or '-':<7} {r['quiet_slips']:>7}"
                     if r["quiet_bias"] else f"{'(none)':<23}")
                tt = f"{float(r['t_trigger']) - 20:6.1f}" if r["t_trigger"] else f"{'-':>6}"
                L.append(f"  g{g_:g} {prot:<8}{label:<20}{c:<10}{q}  | {tt} {r['slips_at']:>7}")
    text = "\n".join(L)
    sim.log(text)
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    figure(val, heat, gam)


def figure(val, heat, gam):
    """One row per site: the two bounds as a bar at gamma = 1 (stored) and at the new gamma (fresh, settled), and the
    threshold with heat."""
    import matplotlib.pyplot as plt
    from snspd import plotting as pl
    pl.apply_style()
    labels = [c[0] for c in CASES]
    fig, ax = plt.subplots(figsize=(8.8, 0.62 * len(labels) + 1.6), layout="constrained")
    grey, light, blue, orange = "#b5b3ac", "#9fc4ee", "#2a78d6", "#eb6834"
    bars = [(-0.22, grey, None), (0.0, light, "fresh"), (0.22, blue, "settled")]
    for k, (label, g, x, y, lo1, hi1) in enumerate(CASES):
        yy = len(labels) - 1 - k
        for dy, col, prot in bars:
            lo, hi = (lo1, hi1) if prot is None else (val(gam, prot, label, "trigger"), val(gam, prot, label, "sustained"))
            if lo is None or hi is None:
                continue
            ax.plot([lo, max(hi, lo + 0.003)], [yy - dy] * 2, color=col, lw=5, solid_capstyle="butt", zorder=1)
        if heat[label][0] is not None:
            ax.plot(heat[label][0], yy + 0.22, "D", ms=6, color=orange, mec=pl.SURFACE, mew=1.2, zorder=3)
    ax.set_yticks(range(len(labels)), labels[::-1], fontsize=9)
    ax.set_xlabel("photon threshold  I / I_dep   (kinetic hotspot, f = 1, 400 tau0, no heat unless marked)")
    ax.grid(axis="y", visible=False)
    handles = [plt.Line2D([], [], color=grey, lw=5, label="gamma = 1 (stored): trigger to sustained"),
               plt.Line2D([], [], color=light, lw=5, label=f"gamma = {gam:g}, photon at 20 tau0 after the bias"),
               plt.Line2D([], [], color=blue, lw=5, label=f"gamma = {gam:g}, photon into the settled state"),
               plt.Line2D([], [], marker="D", ls="", ms=6, color=orange, label="gamma = 1 with heat (one temperature)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=8.5)
    ax.set_title(f"The two bounds at gamma = 1 and gamma = {gam:g}", fontsize=11)
    fig.savefig(os.path.join(OUT, "check.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--gamma", type=float, default=10.0)
    ap.add_argument("--protocols", default="fresh,settled", help="comma-separated: fresh, settled")
    ap.add_argument("--sites", default="", help='";"-separated site labels, or "control" (default: all 14)')
    a = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
