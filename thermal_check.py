#!/usr/bin/env python3
"""
thermal_check.py -- first test of the electro-thermal TDGL (snspd/thermal.py): does adding Joule heat merge the two
bounds on the photon threshold?

For a set of sites whose kinetic-hotspot thresholds are stored under both criteria (voltage trigger = lower bound,
sustained = upper bound), bisect the threshold with heat, where detection is "still resistive at the end of the
400 tau0 window" (a latched normal domain does not heal).  Same photon (kinetic hotspot, f = 1), window, mesh
(pinned production meshes) and film as the stored runs.  The bracket is [lower - 0.10, upper + 0.01] I_dep (the first
one-temperature runs used lower - 0.02).  Heat can take the threshold below the lower bound: a site that switches at
the bottom of its bracket is re-run on the next `run` with lower - 0.10 (the report keeps the latest row per site and
variant).

  python thermal_check.py run    --workers 3            # every variant below that is not done yet
  python thermal_check.py run    --variants "NbTiN 2T"
  python thermal_check.py report

Variants of the heat model (snspd/thermal.py): one temperature with eta = 0.59 (nominal), 0.30 and 1.20 (C_v, T_c and
T0 not measured for this film; factor ~2 either way), and "NbTiN 2T": two temperatures with the parameters measured
on 9 nm NbTiN by Sidorova et al. (thermal.NBTIN_2T).  The photon is the nominal kinetic hotspot in every variant.

Writes results/thermal/check.csv (one row per site and variant) and runs.csv (every run).
"""
from __future__ import annotations

import argparse
import os
import time

import numpy as np

os.environ.setdefault("SNSPD_MESH_DIR", "meshes")

from snspd import geometry as geo
from snspd import sim, store
from snspd.units import Film

OUT = "results/thermal"
TOL = 0.004                        # bisection tolerance [I_dep]


def variants():
    from dataclasses import replace
    from snspd import thermal
    p = thermal.ThermalParams()
    return {"1T eta 0.59": p, "1T eta 0.30": replace(p, eta=0.30), "1T eta 1.20": replace(p, eta=1.20),
            "NbTiN 2T": thermal.NBTIN_2T}


ALL_SITES = ("1T eta 0.59", "NbTiN 2T")    # the other variants run only at VARIANT_CASES

# (label, geometry key or width in nm, x, y, lower bound, upper bound), kinetic hotspot, T400, 1.4 xi (stored results)
CASES = [
    ("30 nm (0,0)", 30, 0, 0, 0.376, 0.775),
    ("30 nm (-13,0)", 30, -13, 0, 0.390, 0.778),
    ("56 nm (0,0)", 56, 0, 0, 0.620, 0.801),
    ("100 nm (0,0)", 100, 0, 0, 0.774, 0.819),
    ("100 nm (-36,0)", 100, -36, 0, 0.732, 0.817),
    ("plain (0,0)", "Plain_Wire", 0, 0, 0.848, 0.848),
    ("plain (-96,0)", "Plain_Wire", -96, 0, 0.832, 0.832),
    ("zipper (0,0)", "Long_Zipper", 0, 0, 0.756, 0.804),
    ("zipper (-68,0)", "Long_Zipper", -68, 0, 0.772, 0.817),
    ("zipper (-96,0)", "Long_Zipper", -96, 0, 0.774, 0.817),
    ("zipper (-42,58)", "Long_Zipper", -42, 58, 0.856, 0.856),
    ("zipper (-36,61)", "Long_Zipper", -36, 61, 0.891, 0.891),
    ("Slot 400 (-58,0)", "Slot_400", -58, 0, 0.762, 0.808),
    ("Slot 400 (-28,200)", "Slot_400", -28, 200, 0.807, 0.830),
]
VARIANT_CASES = ("30 nm (0,0)", "100 nm (0,0)", "plain (0,0)", "zipper (0,0)", "zipper (-68,0)", "Slot 400 (-58,0)")
CHECK_FIELDS = ["variant", "label", "geometry", "x_nm", "y_nm", "eta", "lower", "upper", "i_thermal", "flag", "n_runs",
                "t_trigger", "normal_nm2", "theta_peak", "highest_healed", "healed_slips", "healed_theta_peak", "wall_s"]
RUN_FIELDS = ["variant", "label", "eta", "bias", "latched", "branch", "t_trigger", "slips", "slips_last", "theta_peak",
              "theta_end", "normal_nm2", "joule_xi2"]


def geometry(g):
    return geo.plain_strip(float(g)) if isinstance(g, (int, float)) else geo.get(g)


def bisect_task(task):
    from dataclasses import replace
    from snspd import hotspot as hsm
    from snspd import thermal
    t_start = time.time()
    try:
        film = Film()
        g = geometry(task["geom"])
        cfg = sim.SimConfig(tmax=400.0, mesh_xi=1.4, hotspot="kinetic", out_dir=OUT)
        ctx = sim.context(g, cfg.mesh_xi, film)
        i_dep = film.I_dep_uA(g.free_width_min_nm())
        tp = task["params"]
        lo, hi = task["lower"] - task.get("lo_margin", 0.10), task["upper"] + 0.01
        lo0, hi0, runs, at_thr, healed = lo, hi, [], None, None
        while hi - lo > TOL:
            mid = 0.5 * (lo + hi)
            eps = hsm.build_eps(ctx["hs"], "kinetic", task["x"], task["y"], 1.0, cfg.tmax)
            r = thermal.run_thermal(ctx["device"], mid * i_dep, eps, cfg, tp, hsm.HOTSPOT.t0_tau0)
            runs.append(dict(bias=mid, **r))
            if r["latched"]:
                hi, at_thr = mid, r
            else:
                lo = mid
                healed = dict(bias=mid, **r)          # lo only rises: the last healed run is the highest
        flag = "at_wall" if hi >= hi0 else ("below_bracket" if lo <= lo0 else "")
        return dict(ok=True, task=task, i=hi, flag=flag, runs=runs, at=at_thr, healed=healed,
                    wall_s=time.time() - t_start)
    except Exception as exc:
        import traceback
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}",
                    wall_s=time.time() - t_start)


def variant_of(r):
    return r.get("variant") or f"1T eta {float(r['eta']):.2f}"     # rows of 5 Nov have no variant column


def cmd_run(a):
    os.makedirs(OUT, exist_ok=True)
    latest = {(variant_of(r), r["label"]): r for r in store.load_rows(os.path.join(OUT, "check.csv"))}
    done = {k for k, r in latest.items() if r["flag"] != "below_bracket"}
    names = [v for v in a.variants.split(",") if v] if a.variants else list(variants())
    tasks = []
    for name in names:
        params = variants()[name]
        for label, g, x, y, lo, hi in CASES:
            if name not in ALL_SITES and label not in VARIANT_CASES:
                continue
            if (name, label) in done:
                continue
            tasks.append(dict(variant=name, params=params, label=label, geom=g, x=float(x), y=float(y), lower=lo,
                              upper=hi, eta=params.eta, lo_margin=0.10))
    sim.log(f"-- thermal check: {len(tasks)} bisections")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['label']} eta {t['eta']}: {r['error']}")
            return
        for q in r["runs"]:
            store.append_row(os.path.join(OUT, "runs.csv"), dict(
                variant=t["variant"], label=t["label"], eta=t["eta"], bias=round(q["bias"], 4), latched=int(q["latched"]),
                branch=q["branch"], t_trigger="" if q["t_trigger"] is None else round(q["t_trigger"], 1),
                slips=round(q["slips"], 2), slips_last=round(q["slips_last"], 2), theta_peak=round(q["theta_peak"], 3),
                theta_end=round(q["theta_end"], 3), normal_nm2=round(q["normal_nm2"]), joule_xi2=round(q["joule_xi2"], 1)),
                RUN_FIELDS)
        at, he = r["at"] or {}, r["healed"] or {}
        store.append_row(os.path.join(OUT, "check.csv"), dict(
            variant=t["variant"], label=t["label"], geometry=geometry(t["geom"]).key, x_nm=t["x"], y_nm=t["y"], eta=t["eta"],
            lower=t["lower"], upper=t["upper"], i_thermal=round(r["i"], 4), flag=r["flag"], n_runs=len(r["runs"]),
            t_trigger="" if at.get("t_trigger") is None else round(at["t_trigger"], 1),
            normal_nm2=round(at.get("normal_nm2", 0.0)), theta_peak=round(at.get("theta_peak", 0.0), 3),
            highest_healed="" if not he else round(he["bias"], 4),
            healed_slips="" if not he else round(he["slips"], 2),
            healed_theta_peak="" if not he else round(he["theta_peak"], 3), wall_s=round(r["wall_s"])), CHECK_FIELDS)
        sim.log(f"   {t['label']:<20} {t['variant']}: {r['i']:.3f}  (bounds {t['lower']:.3f}-{t['upper']:.3f}) "
                f"{r['flag']}  {r['wall_s']:.0f} s")

    sim.execute(tasks, a.workers, on_result, fn=bisect_task)
    cmd_report(a)


def cmd_report(a):
    rows = list({(variant_of(r), r["label"]): r for r in store.load_rows(os.path.join(OUT, "check.csv"))}.values())
    if not rows:
        raise SystemExit("no results yet")
    order = {c[0]: k for k, c in enumerate(CASES)}
    vorder = {v: k for k, v in enumerate(variants())}
    rows.sort(key=lambda r: (vorder.get(variant_of(r), 99), order.get(r["label"], 99)))
    L = ["ELECTRO-THERMAL CHECK (snspd/thermal.py, prototype)  kinetic hotspot f = 1, 400 tau0, 1.4 xi, pinned meshes",
         "threshold with heat = lowest bias at which the run is still resistive at 400 tau0 (bisection, tol "
         f"{TOL} I_dep); bounds from the stored kinetic runs: lower = voltage trigger, upper = sustained",
         "position = (i_thermal - lower) / (upper - lower): 0 = at the lower bound, 1 = at the upper bound",
         f"{'variant':<13}{'site':<20}{'lower':>8}{'upper':>8}{'thermal':>9}{'position':>10}{'t_trig':>8}"
         f"{'normal nm2':>12}{'normal domain':>16}   highest healed run (bias, slips, peak theta)",
         "  t_trig: trigger time of the run at the threshold [tau0]; normal nm2: area with theta >= 1 at its end; normal "
         "domain: (highest bias without, lowest with) a normal domain at 400 tau0 among the bisection runs"]
    runs = store.load_rows(os.path.join(OUT, "runs.csv"))

    def normal_bracket(r):
        """(highest run without a normal domain at the end, lowest run with one) for this site and variant."""
        rs = [(float(q["bias"]), float(q["normal_nm2"]) > 0) for q in runs
              if variant_of(q) == variant_of(r) and q["label"] == r["label"]]
        with_nd = [b for b, nd in rs if nd]
        without = [b for b, nd in rs if not nd and (not with_nd or b < min(with_nd))]
        return (max(without) if without else None), (min(with_nd) if with_nd else None)

    for r in rows:
        lo, hi, it = float(r["lower"]), float(r["upper"]), float(r["i_thermal"])
        nb = normal_bracket(r)
        pos = f"{(it - lo) / (hi - lo):.2f}" if hi - lo > 0.005 else "-"
        heal = "" if not r["highest_healed"] else (f"{float(r['highest_healed']):.3f}, {r['healed_slips']} slips, "
                                                   f"theta {r['healed_theta_peak']}")
        L.append(f"{variant_of(r):<13}{r['label']:<20}{lo:>8.3f}{hi:>8.3f}{it:>9.3f}{pos:>10}{r['t_trigger']:>8}"
                 f"{r['normal_nm2']:>12}{(f'{nb[0]:.3f}' if nb[0] else '-'):>8}-{(f'{nb[1]:.3f}' if nb[1] else '-'):<7}   {heal}" + (f"   [{r['flag']}]" if r["flag"] else ""))
    text = "\n".join(L)
    sim.log(text)
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    figure(rows)


def figure(rows):
    """One row per site: the no-heat bounds as a bar, the threshold with heat as markers (nominal one-temperature
    model, two-temperature NbTiN model) and the eta range as a whisker."""
    from snspd import plotting as pl
    import matplotlib.pyplot as plt
    pl.apply_style()
    by = {}
    for r in rows:
        by.setdefault(r["label"], {})[variant_of(r)] = r
    labels = [c[0] for c in CASES if c[0] in by]
    fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(labels) + 1.5), layout="constrained")
    blue, orange = "#2a78d6", "#eb6834"
    for k, lab in enumerate(labels):
        y = len(labels) - 1 - k
        d = by[lab]
        r0 = next(iter(d.values()))
        lo, hi = float(r0["lower"]), float(r0["upper"])
        ax.plot([lo, max(hi, lo + 0.002)], [y, y], color="#c9c7c0", lw=9, solid_capstyle="butt", zorder=1)
        etas = [float(d[v]["i_thermal"]) for v in ("1T eta 0.30", "1T eta 1.20") if v in d]
        if etas and "1T eta 0.59" in d:
            vals = etas + [float(d["1T eta 0.59"]["i_thermal"])]
            ax.plot([min(vals), max(vals)], [y, y], color=pl.INK2, lw=1.2, zorder=2)
            for v in etas:
                ax.plot([v, v], [y - 0.16, y + 0.16], color=pl.INK2, lw=1.2, zorder=2)
        if "1T eta 0.59" in d:
            ax.plot(float(d["1T eta 0.59"]["i_thermal"]), y + 0.12, "o", ms=8, color=blue, mec=pl.SURFACE, mew=1.5,
                    zorder=4)
        if "NbTiN 2T" in d:
            ax.plot(float(d["NbTiN 2T"]["i_thermal"]), y - 0.12, "D", ms=7, color=orange, mec=pl.SURFACE, mew=1.5,
                    zorder=3)
    ax.set_yticks(range(len(labels)), labels[::-1], fontsize=9)
    ax.set_xlabel("photon threshold  I / I_dep   (kinetic hotspot, f = 1, 400 tau0)")
    ax.grid(axis="y", visible=False)
    handles = [plt.Line2D([], [], color="#c9c7c0", lw=9, label="no heat: voltage trigger to sustained (the two bounds)"),
               plt.Line2D([], [], marker="o", ls="", ms=8, color=blue, label="with heat, one temperature, eta = 0.59"),
               plt.Line2D([], [], color=pl.INK2, lw=1.2, label="... eta = 0.30 to 1.20"),
               plt.Line2D([], [], marker="D", ls="", ms=7, color=orange, label="with heat, two temperatures (NbTiN)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=8.5)
    ax.set_title("Electro-thermal TDGL (prototype): where heat puts the threshold", fontsize=11)
    fig.savefig(os.path.join(OUT, "check.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("Writes")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "report"])
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--variants", default="", help="comma-separated variant names (default: all)")
    a = ap.parse_args()
    {"run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
