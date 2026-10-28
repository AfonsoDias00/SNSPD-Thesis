#!/usr/bin/env python3
"""
edge_bank.py -- PARKED: a straight strip with nominal edge banks around a weakened centre
=========================================================================================

Parked on 28 Oct 2026: the thesis focuses on internal holes first.  The 27 Oct sweep (results/ next to this file)
found no plateau gain from edge banks in 16 configurations with either hotspot scaling.  The script runs on the
shared snspd library; run it from the repository root:  python parked/edge_bank/edge_bank.py report

What it simulates
-----------------
py-tdgl caps epsilon at 1, so "strong edges" are built the way bank devices are made: nominal film (epsilon = 1) in
bands of width b along both edges, around a centre weakened to epsilon_c < 1 (a pristine rim around an irradiated
centre).  Every bank strip is compared with two references of the same width, mesh and hotspot:
  * the uniform-epsilon_c strip (b = 0)   -- the fair reference: same absorber, no banks
  * the uniform nominal strip (eps = 1)   -- the plain wire of hitmap.py

The static profile is added to the photon hotspot,  eps(r, t) = eps_s(x) + [eps_hotspot(r, t) - 1] * eps_s(x)^-alpha,
floored at the hotspot's eps_floor.  alpha = 0 (default) puts the same suppression into nominal and weakened film;
--hs-alpha 1 lets the same photon remove a larger share of a weaker condensate (GL condensation energy ~ eps^2, so a
fixed deposited energy suppresses eps by ~ 1 / eps_s) -- the "sensitised absorber" half of the bank idea.  Detection
criterion, bisection, mesh and hotspot model are those of snspd (hitmap.py).

Per configuration
  * I_sw: dark bisection with the static profile
  * the weakest interior line: one dark run at I_sw, nodes > 150 nm from the terminals.  In every strip the switching
    itself starts in the 20-50 nm next to the psi = 0 terminals while the interior sits at the GL depairing state
    (|psi| ~ sqrt(2/3)), so the interior minimum of |psi| / sqrt(eps_s) before that disturbance arrives is what
    tells whether the weakened centre or the bank is closest to failing (x_nuc, psi_min; t_nuc only if it collapses)
  * i_det(x) = I_det / I_sw at photon level f = 1 over the half-width, including sites 3 nm either side of the
    bank/centre interface, which can become a new weak line
  * I_dep of the profile: GL depairing in the uniform-phase-gradient approximation (w << Pearl length),
    I(q) = j_dep * (3 sqrt(3) / 2) * q * integral (eps_s(x) - q^2)_+ dx  at its first local maximum in q.
    Uniform eps gives eps^(3/2) * j_dep * w.

Usage
-----
  python edge_bank.py profile  --bank-nm 30 --eps-c 0.7        # profile on the mesh + I_dep; no solver runs
  python edge_bank.py validate --workers 4                      # short checks; run before any sweep
  python edge_bank.py run      --banks 10,20,30,50 --eps-cs 0.5,0.7,0.9 --workers 8
  python edge_bank.py report
  python edge_bank.py run      --hs-alpha 1 --banks 30,50 --eps-cs 0.5,0.7 --workers 6   # separate files (_a1);
                                                     # dark rows and eps = 1 strips are copied from the alpha = 0 run

The run is resumable (rows are appended to a CSV as they finish).  Defaults: 200 nm strip, mesh 0.8 xi, conserving
hotspot, T = 150 tau0, 2 nm tanh ramp at the bank/centre interface.  The barrier criterion of the ranked plan
(dark-count barrier at I_plateau must not fall) is NOT computed here.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))   # repository root
from snspd import analysis as an          # noqa: E402
from snspd import geometry as geo         # noqa: E402
from snspd import hotspot as hsm          # noqa: E402
from snspd import sim, store              # noqa: E402
from snspd.units import FILM              # noqa: E402

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
NUC_PSI = 0.3           # |psi| / sqrt(eps_s) below this marks where the dark switching starts
TERMINAL_GUARD = 0.15   # [um] nodes this close to the source / drain are ignored for the nucleation site
HS_ALPHA = 0.0          # hotspot scaling exponent, set from --hs-alpha in main()
FIELDS = ["config", "width_nm", "bank_nm", "eps_c", "ramp_nm", "kind", "x_nm", "energy_f", "tmax", "hotspot",
          "mesh_xi", "I_dep_uA", "I_uA", "i_rel", "flag", "n_sims", "t_nuc", "x_nuc_nm", "y_nuc_nm", "psi_min", "psi_centre", "psi_bank", "wall_s"]


# ==============================================================================
# 1. Profile (pure numpy)
# ==============================================================================
def conf_label(bank_nm, eps_c):
    return f"b{bank_nm:g}_e{eps_c:g}"


def eps_profile(x_um, width_nm, bank_nm, eps_c, ramp_nm):
    """Static epsilon across the strip: 1 within bank_nm of either edge, eps_c in the centre."""
    x_um = np.asarray(x_um, float)
    if bank_nm <= 0 or eps_c >= 1.0:
        return np.full(x_um.shape, 1.0 if eps_c >= 1.0 else float(eps_c))
    edge_dist = width_nm / 2.0 - np.abs(x_um * 1e3)                 # nm to the nearest edge
    if ramp_nm > 0:
        s = 0.5 * (1.0 - np.tanh((edge_dist - bank_nm) / ramp_nm))
    else:
        s = (edge_dist <= bank_nm).astype(float)
    return eps_c + (1.0 - eps_c) * s


def i_dep_profile(width_nm, bank_nm, eps_c, ramp_nm):
    """GL depairing current [uA] of the profile (uniform phase gradient q in 1/xi, first local maximum of I(q))."""
    x = np.linspace(-width_nm / 2.0, width_nm / 2.0, 1001)
    eps = eps_profile(x * 1e-3, width_nm, bank_nm, eps_c, ramp_nm)
    qs = np.linspace(0.0, 1.0, 4001)
    dens = np.clip(eps[None, :] - qs[:, None] ** 2, 0.0, None)
    I = qs * np.trapezoid(dens, x, axis=1) if hasattr(np, "trapezoid") else qs * np.trapz(dens, x, axis=1)
    k = int(np.argmax(np.diff(I) < 0)) if np.any(np.diff(I) < 0) else len(I) - 1
    return FILM.j_dep_uA_per_nm * 1.5 * np.sqrt(3.0) * float(I[k])


def hit_positions(width_nm, bank_nm):
    """Sites over the left half-width.  Banks: 2 nm from the edge, mid-bank, 3 nm either side of the interface,
    12 nm into the centre; uniform strips: 2, 8, 14, 21, 30 nm.  Both: the quarter points and the centre."""
    h = width_nm / 2.0
    offs = [2.0, bank_nm / 2.0, bank_nm - 3.0, bank_nm + 3.0, bank_nm + 12.0] if bank_nm > 0 else [2.0, 8.0, 14.0, 21.0, 30.0]
    offs += [h / 2.0, 0.75 * h, h]
    kept = []
    for o in offs:
        if 0.0 < o <= h and all(abs(o - k) >= 2.5 for k in kept):
            kept.append(o)
    return sorted(round(-(h - o), 2) + 0.0 for o in kept)


# ==============================================================================
# 2. Solver tasks (run in worker processes; everything they need travels in the task dict)
# ==============================================================================
def context(width_nm, mesh_xi):
    """Device, monitored right-half slice and hotspot solver of a plain strip (shared cache in snspd.sim)."""
    return sim.context(geo.plain_strip(width_nm), mesh_xi)


def static_eps(conf):
    def eps_func(r, *, vectorized=True):
        return eps_profile(np.asarray(r)[:, 0], conf["width"], conf["bank"], conf["eps_c"], conf["ramp"])
    return eps_func


def with_hotspot(static_fn, hot, alpha=0.0):
    cache = {}

    def eps_func(r, *, t, vectorized=True):
        r = np.asarray(r)
        if cache.get("n") != len(r):
            s = static_fn(r)
            cache["s"], cache["w"], cache["n"] = s, np.maximum(s, 0.05) ** -alpha, len(r)
        return np.maximum(cache["s"] + (hot(r, t=t) - 1.0) * cache["w"], hsm.HOTSPOT.eps_floor)

    return eps_func


def nucleation(ctx, conf, current, cfg):
    """One dark run at `current`, interior only (the 20-50 nm next to the psi = 0 terminals deplete first in every
    strip and would otherwise always be reported).  If the interior collapses (|psi| / sqrt(eps_s) < NUC_PSI):
    time and place.  If not: t_nuc is blank and x/y/psi_min give the weakest interior line between t = 30 and
    0.6 tmax, i.e. in the steady state before the terminal disturbance arrives."""
    import tdgl
    dev = ctx["device"]
    pts = np.asarray(dev.points)
    ok = np.abs(pts[:, 1]) < geo.plain_strip(conf["width"]).length_nm * 5e-4 - TERMINAL_GUARD
    ref = np.sqrt(np.maximum(eps_profile(pts[:, 0], conf["width"], conf["bank"], conf["eps_c"], conf["ramp"]), 1e-3))
    edge_dist = (conf["width"] / 2.0 - np.abs(pts[ok, 0] * 1e3))
    regions = {"psi_centre": edge_dist >= conf["bank"] + 2.0 * conf["ramp"],
               "psi_bank": (edge_dist <= conf["bank"] - 2.0 * conf["ramp"]) if conf["bank"] > 0 else np.zeros(ok.sum(), bool)}
    reg_min = {k: np.inf for k in regions}
    h5 = sim.temp_h5(cfg)
    sol = None
    try:
        sol = tdgl.solve(dev, sim.solver_options(cfg, cfg.tmax, h5), terminal_currents=dict(source=current, drain=-current),
                         disorder_epsilon=static_eps(conf))
        steps, times = sim.saved_steps_and_times(sol)
        weakest = (np.inf, None)
        for s, tm in zip(steps, times):
            sol.solve_step = int(s)
            a = (np.abs(sol.tdgl_data.psi) / ref)[ok]
            k = int(np.argmin(a))
            x, y = pts[ok][k] * 1e3
            if a[k] < NUC_PSI:
                return dict(t_nuc=round(float(tm), 2), x_nuc_nm=round(float(x), 1), y_nuc_nm=round(float(y), 1),
                            psi_min=round(float(a[k]), 3))
            if 30.0 <= tm <= 0.6 * cfg.tmax:
                if a[k] < weakest[0]:
                    weakest = (float(a[k]), (x, y))
                for name, m in regions.items():
                    if m.any():
                        reg_min[name] = min(reg_min[name], float(a[m].min()))
        x, y = weakest[1] if weakest[1] is not None else ("", "")
        out = dict(t_nuc="", x_nuc_nm=x if x == "" else round(float(x), 1), y_nuc_nm=y if y == "" else round(float(y), 1),
                   psi_min=round(weakest[0], 3) if np.isfinite(weakest[0]) else "")
        out.update({k: round(v, 3) if np.isfinite(v) else "" for k, v in reg_min.items()})
        return out
    finally:
        path = getattr(sol, "path", None) or h5
        if path and os.path.exists(path):
            os.remove(path)


def task_fn(task):
    """Run one task. Never raises (a failed task is simply not recorded and is retried on resume)."""
    t_start = time.time()
    try:
        conf, kind = task["conf"], task["kind"]
        i_dep = i_dep_profile(conf["width"], conf["bank"], conf["eps_c"], conf["ramp"])
        cfg = sim.SimConfig(**task["cfg"])
        tol, collapse = sim.TOL_REL * i_dep, sim.COLLAPSE_REL * i_dep       # both scale with the profile's I_dep
        ctx = context(conf["width"], cfg.mesh_xi)
        res = dict(I_dep_uA=round(i_dep, 2))
        if kind == "isw":
            val, flag, n = sim.bisect(ctx["device"], ctx["slices"], 0.55 * i_dep, 1.15 * i_dep, static_eps(conf), cfg,
                                      tol, collapse)
            res.update(I_uA=val, i_rel=val / i_dep, flag="above_bracket" if flag == "at_wall" else flag, n_sims=n)
        elif kind == "nucleation":
            res.update(nucleation(ctx, conf, task["I_sw"], dataclasses.replace(cfg, save_every=20)), I_uA=task["I_sw"], n_sims=1)
        elif kind == "hit":
            hot = hsm.build_eps(ctx["hs"], cfg.hotspot, task["x_nm"], 0.0, task["f"], cfg.tmax)
            i_sw = task["I_sw"]
            val, flag, n = sim.bisect(ctx["device"], ctx["slices"], cfg.hit_lo_frac * i_sw, i_sw,
                                      with_hotspot(static_eps(conf), hot, task.get("alpha", 0.0)), cfg, tol, collapse)
            res.update(I_uA=val, i_rel=val / i_sw, flag=flag, n_sims=n)
        else:
            raise ValueError(kind)
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:  # keep the pool alive
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


# ==============================================================================
# 3. Orchestration with resumable CSV
# ==============================================================================
def stem():
    return "edge_bank" if HS_ALPHA == 0 else f"edge_bank_a{HS_ALPHA:g}"


def csv_path(cfg):
    return cfg.path(stem(), "csv")


def load_rows(cfg):
    return store.load_rows(csv_path(cfg))


def row_key(config, kind, x="", f=""):
    return (config, kind, "" if x == "" else round(float(x), 2), "" if f == "" else round(float(f), 3))


def append_row(cfg, result):
    t, res = result["task"], result["res"]
    c = t["conf"]
    row = dict(config=conf_label(c["bank"], c["eps_c"]), width_nm=c["width"], bank_nm=c["bank"], eps_c=c["eps_c"],
               ramp_nm=c["ramp"], kind=t["kind"], x_nm=t.get("x_nm", ""), energy_f=t.get("f", ""), tmax=cfg.tmax,
               hotspot=cfg.hotspot, mesh_xi=cfg.mesh_xi, wall_s=round(result["wall_s"], 1))
    row.update({k: res.get(k, "") for k in FIELDS if k in res})
    store.append_row(csv_path(cfg), row, FIELDS)
    return row


def execute(cfg, tasks, workers):
    if not tasks:
        return
    sim.log(f"-- {len(tasks)} tasks, {workers} worker(s)")

    def handle(r):
        t = r["task"]
        tag = f"{conf_label(t['conf']['bank'], t['conf']['eps_c'])} {t['kind']}" + (f" x={t['x_nm']:g} f={t['f']:g}" if t["kind"] == "hit" else "")
        if not r["ok"]:
            sim.log(f"   FAILED {tag}: {r['error']}")
            return
        row = append_row(cfg, r)
        extra = (f"x_nuc={row['x_nuc_nm']} nm t={row['t_nuc']}" if t["kind"] == "nucleation"
                 else f"I={float(row['I_uA']):.1f} uA  i={float(row['i_rel']):.3f}" + (f" [{row['flag']}]" if row["flag"] else ""))
        sim.log(f"   {tag:<34} {extra}   ({row['n_sims']} runs, {r['wall_s']:.0f} s)")

    sim.execute(tasks, workers, handle, fn=task_fn)


def isw_of(rows):
    return {r["config"]: float(r["I_uA"]) for r in rows if r["kind"] == "isw"}


def run_configs(cfg, a, confs, sites=None):
    """Phase 1: I_sw of every configuration.  Phase 2: nucleation run + photon thresholds (f = 1)."""
    base = dict(cfg=dataclasses.asdict(cfg), alpha=HS_ALPHA)
    seed_from_alpha0(cfg, confs)
    have = {row_key(r["config"], r["kind"], r["x_nm"], r["energy_f"]) for r in load_rows(cfg)}
    execute(cfg, [dict(base, conf=c, kind="isw") for c in confs
                  if row_key(conf_label(c["bank"], c["eps_c"]), "isw") not in have], a.workers)
    isw = isw_of(load_rows(cfg))
    tasks = []
    for c in confs:
        lab = conf_label(c["bank"], c["eps_c"])
        if lab not in isw:
            continue
        if row_key(lab, "nucleation") not in have:
            tasks.append(dict(base, conf=c, kind="nucleation", I_sw=isw[lab]))
        for x in (sites or {}).get(lab, hit_positions(c["width"], c["bank"])):
            if row_key(lab, "hit", x, 1.0) not in have:
                tasks.append(dict(base, conf=c, kind="hit", x_nm=float(x), f=1.0, I_sw=isw[lab]))
    execute(cfg, tasks, a.workers)


def seed_from_alpha0(cfg, confs):
    """Dark rows (I_sw, interior diagnostic) do not involve the hotspot, and alpha does nothing where eps_s = 1:
    copy those rows from the alpha = 0 file of the same tag instead of re-running them."""
    src = cfg.path("edge_bank", "csv")
    if HS_ALPHA == 0 or not os.path.exists(src):
        return
    labels = {conf_label(c["bank"], c["eps_c"]) for c in confs}
    have = {row_key(r["config"], r["kind"], r["x_nm"], r["energy_f"]) for r in load_rows(cfg)}
    with open(src, newline="") as fh:
        take = [r for r in csv.DictReader(fh) if r["config"] in labels
                and (r["kind"] in ("isw", "nucleation") or float(r["eps_c"]) >= 1.0)
                and row_key(r["config"], r["kind"], r["x_nm"], r["energy_f"]) not in have]
    if take:
        new = not os.path.exists(csv_path(cfg))
        with open(csv_path(cfg), "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            if new:
                w.writeheader()
            w.writerows({k: r.get(k, "") for k in FIELDS} for r in take)
        sim.log(f"-- reused {len(take)} hotspot-independent rows from {src}")


def make_conf(a, bank, eps_c):
    return dict(width=a.width_nm, bank=float(bank), eps_c=float(eps_c), ramp=a.ramp_nm)


def sweep_confs(a):
    eps_cs = [float(s) for s in a.eps_cs.split(",") if s]
    confs = [make_conf(a, 0, 1.0)] + [make_conf(a, 0, e) for e in eps_cs]
    confs += [make_conf(a, b, e) for e in eps_cs for b in (float(s) for s in a.banks.split(",") if s)]
    return confs


# ==============================================================================
# 4. Commands
# ==============================================================================
def mesh_bank_width(ctx, conf):
    """Bank width as resolved on the mesh: area with eps > (1 + eps_c) / 2, per edge, divided by the strip length."""
    hs = ctx["hs"]
    eps = eps_profile(hs.P[:, 0], conf["width"], conf["bank"], conf["eps_c"], conf["ramp"])
    return float(hs.m[eps > 0.5 * (1.0 + conf["eps_c"])].sum() / (geo.plain_strip(conf["width"]).length_nm * 1e-3) / 2.0 * 1e3)


def cmd_profile(cfg, a):
    conf = make_conf(a, a.bank_nm, a.eps_c)
    ctx = context(conf["width"], cfg.mesh_xi)
    xs = np.array(hit_positions(conf["width"], conf["bank"]))
    sim.log(f"profile {conf_label(a.bank_nm, a.eps_c)}: width {a.width_nm:g} nm, ramp {a.ramp_nm:g} nm, "
          f"mesh {cfg.mesh_xi:g} xi = {cfg.mesh_xi * FILM.xi_nm:.1f} nm, {len(ctx['device'].points)} nodes")
    sim.log(f"  bank width on the mesh   {mesh_bank_width(ctx, conf):.1f} nm  (set {a.bank_nm:g} nm)")
    for lab, b, e in (("this profile", a.bank_nm, a.eps_c), ("uniform eps_c", 0, a.eps_c), ("uniform nominal", 0, 1.0)):
        sim.log(f"  I_dep {lab:<16} {i_dep_profile(a.width_nm, b, e, a.ramp_nm):8.1f} uA")
    sim.log("  hit sites x [nm]: " + ", ".join(f"{x:g}" for x in xs))
    sim.log("  eps_s at those sites: " + ", ".join(f"{v:.3f}" for v in eps_profile(xs * 1e-3, a.width_nm, a.bank_nm, a.eps_c, a.ramp_nm)))


def reference(cfg):
    """Plain-wire I_sw and the centre f = 1 threshold from hitmap.py at the same tag, if that run exists."""
    tag_cfg = dataclasses.replace(cfg, out_dir="results/hitmap")
    isw = store.load_json(tag_cfg.path("isw", "json"), {}).get("Plain_Wire", {}).get("I_sw")
    i0 = None
    if os.path.exists(tag_cfg.path("hitmap", "csv")):
        with open(tag_cfg.path("hitmap", "csv"), newline="") as fh:
            for r in csv.DictReader(fh):
                if r["Geometry"] == "Plain_Wire" and float(r["x_nm"]) == 0.0 and float(r["energy_f"]) == 1.0:
                    i0 = float(r["i_det"])
    return isw, i0


def cmd_validate(cfg, a):
    vb, ve = a.val_bank, a.val_eps
    nominal, uniform, bank = make_conf(a, 0, 1.0), make_conf(a, 0, ve), make_conf(a, vb, ve)
    sim.log("== V0  profile on the mesh (no solver)")
    ctx = context(a.width_nm, cfg.mesh_xi)
    hs_eps = eps_profile(ctx["hs"].P[:, 0], a.width_nm, vb, ve, a.ramp_nm)
    checks = [("V0 eps <= 1 everywhere", hs_eps.max() <= 1.0, f"max {hs_eps.max():.4f}, min {hs_eps.min():.4f}"),
              ("V0 bank width resolved", abs(mesh_bank_width(ctx, bank) - vb) <= 0.5 * cfg.mesh_xi * FILM.xi_nm,
               f"{mesh_bank_width(ctx, bank):.1f} nm on the mesh vs {vb:g} nm set")]
    sites ={conf_label(0, 1.0): [0.0], conf_label(0, ve): [0.0],
             conf_label(vb, ve): [0.0, -(a.width_nm / 2.0 - vb / 2.0)]}
    sim.log("== V1-V3  short solver runs: nominal, uniform eps_c and one bank strip")
    run_configs(cfg, a, [nominal, uniform, bank], sites)
    rows = load_rows(cfg)
    isw = isw_of(rows)
    get = lambda c, k, x="": next((r for r in rows if r["config"] == c and r["kind"] == k
                                   and (x == "" or float(r["x_nm"]) == x)), None)
    ln, lu, lb = conf_label(0, 1.0), conf_label(0, ve), conf_label(vb, ve)
    if any(c not in isw for c in (ln, lu, lb)):
        sim.log(f"  FAIL  I_sw missing for {[c for c in (ln, lu, lb) if c not in isw]} -- see the FAILED lines above")
        return
    ref_isw, ref_i0 = reference(cfg)
    tol_n = sim.TOL_REL * i_dep_profile(a.width_nm, 0, 1.0, a.ramp_nm)
    if ref_isw:
        checks.append(("V1 static-eps path reproduces plain I_sw", abs(isw[ln] - ref_isw) <= 2 * tol_n,
                       f"{isw[ln]:.1f} vs {ref_isw:.1f} uA (hitmap.py, same tag), tol {tol_n:.1f}"))
    if ref_i0 is not None and get(ln, "hit", 0.0):
        v = float(get(ln, "hit", 0.0)["i_rel"])
        checks.append(("V1 centre threshold reproduces plain wire", abs(v - ref_i0) <= 0.007, f"{v:.3f} vs {ref_i0:.3f}"))
    r_n = isw[ln] / i_dep_profile(a.width_nm, 0, 1.0, a.ramp_nm)
    r_u = isw[lu] / i_dep_profile(a.width_nm, 0, ve, a.ramp_nm)
    r_b = isw[lb] / i_dep_profile(a.width_nm, vb, ve, a.ramp_nm)
    checks.append(("V2 uniform eps_c follows GL eps^(3/2)", abs(r_u / r_n - 1) <= 0.03,
                   f"I_sw/I_dep {r_u:.3f} vs nominal {r_n:.3f}"))
    checks.append(("V3 banks sit between the two references", isw[lu] < isw[lb] < isw[ln],
                   f"{isw[lu]:.1f} < {isw[lb]:.1f} < {isw[ln]:.1f} uA"))
    checks.append(("V3 bank I_sw / I_dep(profile) near nominal", abs(r_b / r_n - 1) <= 0.05,
                   f"{r_b:.3f} vs {r_n:.3f}  (depairing model check)"))
    flagged = [f"{r['config']} {r['kind']} {r['x_nm']}: {r['flag']}" for r in rows
               if r["config"] in (ln, lu, lb) and r["flag"]]
    checks.append(("no bracket flags", not flagged, "; ".join(flagged) or "none"))
    sim.log("-" * 96)
    for name, ok, info in checks:
        sim.log(f"  {'PASS' if ok else 'FAIL'}  {name:<46} {info}")
    for c in (ln, lu, lb):
        n = get(c, "nucleation")
        if n:
            sim.log(f"  weakest interior line  {c:<10} x = {n['x_nuc_nm']} nm, y = {n['y_nuc_nm']} nm, "
                  f"|psi|/sqrt(eps) = {n['psi_min']}  (centre {n['psi_centre']}, bank {n['psi_bank'] or '-'})"
                  + (f", collapses at t = {n['t_nuc']} tau_0" if n["t_nuc"] else ""))
    for c in (lu, lb):
        for r in rows:
            if r["config"] == c and r["kind"] == "hit":
                sim.log(f"  photon f=1  {c:<10} x = {float(r['x_nm']):7.1f} nm  i_det = {float(r['i_rel']):.3f}")
    hits = [r for r in rows if r["kind"] == "hit" and r["config"] in (ln, lu, lb)]
    if hits:
        per_sim = np.mean([float(r["wall_s"]) / float(r["n_sims"]) for r in hits])
        n_conf = 1 + len(a.eps_cs.split(",")) * (1 + len(a.banks.split(",")))
        n_sims = n_conf * (9 + 1 + 8 * 9)
        sim.log(f"  cost: {per_sim:.0f} s per solver run -> full sweep ({n_conf} configs) ~ {n_sims} runs ~ "
              f"{n_sims * per_sim / 3600:.1f} h serial, ~{n_sims * per_sim / 3600 / max(a.workers, 1):.1f} h on {a.workers} workers")
    sim.log("-" * 96)
    sim.log(f"  VALIDATION {'PASSED' if all(ok for _, ok, _ in checks) else 'FAILED -- fix before sweeping'}")


def collect(rows):
    confs = {}
    for r in rows:
        c = confs.setdefault(r["config"], dict(bank=float(r["bank_nm"]), eps_c=float(r["eps_c"]), width=float(r["width_nm"]),
                                               ramp=float(r["ramp_nm"]), hits=[]))
        if r["kind"] == "isw":
            c.update(I_sw=float(r["I_uA"]), I_dep=float(r["I_dep_uA"]), isw_flag=r["flag"])
        elif r["kind"] == "nucleation":
            c["weak"] = (f"{float(r['psi_centre']):.2f}/{float(r['psi_bank']):.2f}" if r["psi_bank"]
                         else f"{float(r['psi_centre']):.2f}/-" if r["psi_centre"] else "")
        elif r["kind"] == "hit" and float(r["energy_f"]) == 1.0:
            c["hits"].append((float(r["x_nm"]), 0.0, 1.0 if r["flag"] == "at_wall" else float(r["i_rel"]), r["flag"]))
    for c in confs.values():
        if c["hits"] and "I_sw" in c:
            g = geo.plain_strip(c["width"])
            c["met"] = an.metrics(an.pixel_values(an.unit_map(g, [h[:3] for h in c["hits"]])[2]))
    return confs


def cmd_report(cfg, a):
    rows = load_rows(cfg)
    if not rows:
        sim.log("no rows yet.")
        return
    confs = collect(rows)
    base = {}
    if HS_ALPHA != 0 and os.path.exists(cfg.path("edge_bank", "csv")):
        with open(cfg.path("edge_bank", "csv"), newline="") as fh:
            base = collect(list(csv.DictReader(fh)))
    i_dep_nom = i_dep_profile(a.width_nm, 0, 1.0, a.ramp_nm)
    order = sorted(confs, key=lambda k: (-confs[k]["eps_c"], confs[k]["bank"]))
    L = []
    P = L.append
    P("=" * 118)
    P(f"EDGE-BANK STRIP   width={a.width_nm:g} nm   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau0   mesh={cfg.mesh_xi:g} xi   "
      f"ramp={a.ramp_nm:g} nm   photon f = 1   hotspot scaling alpha={HS_ALPHA:g}")
    P("=" * 118)
    P(f"{'config':<11}{'b':>5}{'eps_c':>7}{'I_sw':>9}{'I_dep':>9}{'Isw/Idep':>9}{'psi c/b':>11}{'i_10':>8}{'i_50':>8}"
      f"{'i_max':>8}{'I_plat':>9}{'Iplat/':>8}{'d i_max':>9}{'sites':>6}  verdict")
    P(f"{'':<11}{'[nm]':>5}{'':>7}{'[uA]':>9}{'[uA]':>9}{'':>9}{'':>11}{'[Isw]':>8}{'[Isw]':>8}{'[Isw]':>8}{'[uA]':>9}"
      f"{'Idep1':>8}{'vs unif':>9}")
    for k in order:
        c = confs[k]
        if "I_sw" not in c:
            continue
        m = c.get("met")
        uni = confs.get(conf_label(0, c["eps_c"]))
        line = (f"{k:<11}{c['bank']:>5g}{c['eps_c']:>7g}{c['I_sw']:>9.1f}{c['I_dep']:>9.1f}{c['I_sw'] / c['I_dep']:>9.3f}"
                f"{c.get('weak', ''):>11}")
        if not m:
            P(line + "   (no photon thresholds yet)")
            continue
        plat = m["i_max"] * c["I_sw"]
        line += (f"{m['i_10']:>8.3f}{m['i_50']:>8.3f}{m['i_max']:>8.3f}{plat:>9.1f}"
                 f"{plat / i_dep_nom if i_dep_nom else float('nan'):>8.3f}")
        verdict = ""
        n_expected = len(hit_positions(c["width"], c["bank"]))
        if len(c["hits"]) < n_expected or (uni and len(uni["hits"]) < len(hit_positions(uni["width"], 0))):
            line += f"{'':>9}"
            verdict = f"partial ({len(c['hits'])}/{n_expected} sites), no verdict"
        elif c["bank"] > 0 and uni and uni.get("met"):
            d = m["i_max"] - uni["met"]["i_max"]
            line += f"{d:>+9.3f}"
            better = d < -0.007 and c["I_sw"] / c["I_dep"] >= uni["I_sw"] / uni["I_dep"] - 0.01
            verdict = "candidate (barrier not checked)" if better else "no gain"
        else:
            line += f"{'':>9}"
        flags = sorted({h[3] for h in c["hits"] if h[3]} | ({c["isw_flag"]} if c.get("isw_flag") else set()))
        b0 = base.get(k, {}).get("met")
        P(line + f"{len(c['hits']):>6}  {verdict}" + (f"  FLAGS {flags}" if flags else "")
          + (f"   [alpha=0: i_50 {b0['i_50']:.3f}, i_max {b0['i_max']:.3f}]" if b0 else ""))
    P("-" * 118)
    P("i_max = I_plateau / I_sw (100 % point); d i_max vs unif = bank minus uniform-eps_c strip (negative = plateau at lower")
    P("normalised bias).  'candidate' needs d i_max < -0.007 at equal or higher I_sw / I_dep; the ranked plan also requires the")
    P("dark-count barrier at I_plateau not to fall, which this script does not compute.  psi c/b: lowest |psi|/sqrt(eps_s) in")
    P("the centre / the banks at I_sw, > 150 nm from the terminals, before switching (GL depairing state = sqrt(2/3) = 0.816;")
    P("the lower of the two is the line closest to failing).")
    P("=" * 118)
    txt = "\n".join(L)
    print(txt)
    with open(cfg.path(stem() + "_summary", "txt"), "w") as fh:
        fh.write(txt + "\n")


def parse():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["profile", "validate", "run", "report"])
    ap.add_argument("--width-nm", type=float, default=200.0)
    ap.add_argument("--hotspot", choices=["conserving", "legacy"], default="conserving")
    ap.add_argument("--tmax", type=float, default=150.0)
    ap.add_argument("--mesh-xi", type=float, default=0.8)
    ap.add_argument("--ramp-nm", type=float, default=2.0, help="tanh half-width of the bank/centre interface [nm]")
    ap.add_argument("--hit-lo-frac", type=float, default=0.30, help="lower photon-threshold bracket [fraction of I_sw]")
    ap.add_argument("--bank-nm", type=float, default=30.0, help="profile: bank width")
    ap.add_argument("--eps-c", type=float, default=0.7, help="profile: centre epsilon")
    ap.add_argument("--val-bank", type=float, default=30.0, help="validate: bank width of the test strip")
    ap.add_argument("--val-eps", type=float, default=0.7, help="validate: centre epsilon of the test strip")
    ap.add_argument("--banks", default="10,20,30,50")
    ap.add_argument("--eps-cs", default="0.5,0.7,0.9")
    ap.add_argument("--hs-alpha", type=float, default=0.0, help="hotspot suppression scales as eps_s^-alpha")
    ap.add_argument("--workers", type=int, default=1, help="each worker needs ~1.8 GB at 200 nm / 0.8 xi")
    ap.add_argument("--out", default=OUT_DIR)
    return ap.parse_args()


def main():
    global HS_ALPHA
    a = parse()
    HS_ALPHA = a.hs_alpha
    sim.need_tdgl()
    cfg = sim.SimConfig(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi, hit_lo_frac=a.hit_lo_frac, out_dir=a.out)
    os.makedirs(cfg.out_dir, exist_ok=True)
    if a.command == "profile":
        cmd_profile(cfg, a)
    elif a.command == "validate":
        cmd_validate(cfg, a)
    elif a.command == "run":
        run_configs(cfg, a, sweep_confs(a))
        cmd_report(cfg, a)
    else:
        cmd_report(cfg, a)


if __name__ == "__main__":
    main()
