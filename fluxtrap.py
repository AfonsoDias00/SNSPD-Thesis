#!/usr/bin/env python3
"""
fluxtrap.py -- magnetic flux that a photon leaves trapped in the holes.

A vortex that a photon lets into the film and that reaches a hole is held there: leaving again means crossing more
film.  The hole then carries a fluxoid (an integer winding of the order-parameter phase around it, see
snspd/fluxoid.py) and a circulating current that adds to the transport current on one side of the hole, so after
such a photon the device is no longer the one the photon found.

  scan     every hit site of a geometry (full grid) at fixed biases, in units of I_dep of
           the narrowest cross-section, over the full window: did the voltage trigger fire and when, did the device
           stay resistive (both together: the detection criterion of snspd/sim.py), the phase slipped between the probes since the photon (|int V dt| / 2 pi: one vortex crossing the full width
           slips 1; a vortex that stops in a hole slips only the fraction of the width it crossed, ~0.3-0.4 here) and
           the fluxoid of every hole on every saved frame.  Each run is one of
             clean      no trigger, no flux in the holes at the end of the window
             trapped    no trigger, flux left in the holes
             transient  the voltage trigger fired but the device did not stay resistive (< sim.SUSTAIN_SLIPS slips
                        in the last sim.SUSTAIN_WIN tau0 of the window or of the sim.SETTLE continuation), e.g. on flux
                        moving between holes: not a detection
             direct     switched (detected), without first sitting in a zero-voltage trapped state
             via_trap   switched (detected) after >= DWELL_MIN tau0 in a zero-voltage trapped state
                        (max |V| < V_QUIET) before the resistive onset (ONSET_SLIPS slips since the photon)
  trapped  the most common trapped states from `scan` (or states given with --from, or placed directly with --place:
           the phase is wound around the holes and relaxed for PLACE_RELAX tau0 at zero bias, no photon),
           regenerated and used as the initial state of new runs: their
           own switching current (no photon) and what happens at it (a sustained switch, or flux that only moves to
           another hole), whether they survive at zero bias, and the photon threshold at a few sites compared with the
           pristine device
  report   per-site thresholds, area fractions, figures
  trace    one site at a few biases: probe voltage and the fluxoid of every hole vs time (text and figure)

Usage
-----
  python fluxtrap.py scan    --workers 9
  python fluxtrap.py trapped --workers 6
  python fluxtrap.py report
  python fluxtrap.py trapped --mesh-xi 0.8 --from=-68,0@0.77 --no-second      (mesh check of one trapped state)
  python fluxtrap.py trapped --geom Slot_400 --hitmap-dir results/slots --place="0 -1 0" --no-second
  python fluxtrap.py trace   --site=-68,0 --biases 0.70,0.77,0.80 --workers 3
  python fluxtrap.py scan    --geom Slot_400 --hitmap-dir results/slots --biases 0.76:0.86:0.02 --sites="-74,0;-42,190"

Biases are "start:stop:step" ranges and single values, comma separated, in units of I_dep.  The window, mesh and
hotspot default to the current hitmap production settings (400 tau0, 1.4 xi, conserving).
"""
from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter

import numpy as np

from snspd import analysis as an
from snspd import geometry as geo
from snspd import sim, store
from snspd.units import Film

OUT_DIR = "results/fluxtrap"
HITMAP_DIR = "results/hitmap"
V_QUIET = 0.05           # max |V| between two frames for the device to count as sitting still [V0]
DWELL_MIN = 20.0         # zero-voltage trapped time before a switch for it to count as via_trap [tau0]
PROMPT_T = 150.0         # "prompt" by latency: resistive onset by this time (130 tau0 after the photon) [tau0]
ONSET_SLIPS = 2.0        # resistive onset: this many phase slips since the photon (beyond any trapping)
PLACE_RELAX = 50.0       # relaxation of a placed fluxoid state at zero bias before it is used [tau0]
SAVE_EVERY = 40          # solver steps between saved frames (<= 8 tau0 at dt_max = 0.2)
PROBE_SITES = ((-68.0, 0.0), (-18.0, 0.0), (0.0, 0.0), (68.0, 61.0))   # second-photon sites [nm]
SCAN_FIELDS = ["Geometry", "x_nm", "y_nm", "bias", "I_uA", "tmax", "hotspot", "mesh_xi", "outcome", "branch",
               "t_trigger", "slips", "slips_last", "slips_settled", "t_onset", "n_end", "n_abs", "t_trap", "dwell", "history", "wall_s"]
OUTCOMES = ("clean", "trapped", "transient", "direct", "via_trap")
FIRED = ("transient", "direct", "via_trap")       # the voltage trigger fired
SUSTAINED = ("direct", "via_trap")                # ... and the device stayed resistive: detected


def parse_biases(text):
    out = []
    for part in text.split(","):
        if ":" in part:
            a, b, s = (float(v) for v in part.split(":"))
            out += list(np.round(np.arange(a, b + s / 2, s), 4))
        elif part:
            out.append(round(float(part), 4))
    return sorted(set(float(v) for v in out))


def fmt_n(n):
    return " ".join(str(int(v)) for v in n)


def parse_n(text):
    return tuple(int(v) for v in text.split()) if text else ()


# ==============================================================================
# One run with the fluxoid record
# ==============================================================================
_METERS = {}


def meter_for(device):
    from snspd.fluxoid import FluxoidMeter
    if id(device) not in _METERS:
        _METERS[id(device)] = FluxoidMeter(device)
    return _METERS[id(device)]


def classify(trig, frames, t, v, t0):
    """Outcome, phase slips, resistive onset, flux state before the onset (or at the end), first trapping time,
    zero-voltage trapped time before the onset, n history."""
    from snspd.fluxoid import PSI_LOOP_MIN
    sl = sim.phase_slips(t, v, t0)
    cum = sl[-1]
    k = np.flatnonzero(sl >= ONSET_SLIPS)
    t_onset = float(t[k[0]]) if k.size else None
    sustained = trig.switched
    t_end = t_onset if sustained and t_onset is not None else np.inf
    valid = [(t, n, vmax) for t, n, floor, vmax in frames if floor >= PSI_LOOP_MIN and t <= t_end]
    n_end, t_trap, dwell, hist, prev = (), None, 0.0, [], None
    t_prev = 0.0
    for t, n, vmax in valid:
        n = tuple(n)
        if any(n):
            t_trap = t if t_trap is None else t_trap
            if vmax < V_QUIET:
                dwell += t - t_prev
        if n != prev:
            hist.append(f"{t:.0f}:{','.join(map(str, n))}")
            prev = n
        n_end, t_prev = n, t
    if sustained:
        outcome = "via_trap" if dwell >= DWELL_MIN else "direct"
    elif trig.branch:
        outcome = "transient"
    else:
        outcome = "trapped" if any(n_end) else "clean"
    return dict(outcome=outcome, slips=float(cum), slips_last=float(trig.slips_last),
                slips_settled=trig.slips_settled, t_onset=t_onset, n_end=n_end,
                t_trap=t_trap, dwell=dwell, history=" ".join(hist), frames=frames)


def observed_run(device, I, eps, cfg, seed=None, phase=None):
    """One run: (Trigger, classification, final Solution with its last frame loaded -- usable as a seed).  Phase
    slips count from the photon, or from t = 0 without one.  phase: see sim.solved."""
    from snspd.fluxoid import frames
    t0 = cfg.pset.params.t0_tau0
    init = None
    if seed is None and phase is None and eps is not None and cfg.pset.settled_start:
        init = sim.settled_state(device, I, cfg)          # a photon run of a preset that starts settled (sim.py)
    with sim.solved(device, I, eps, cfg, seed, phase, init=init) as sol:
        trig = sim.settle(device, I, sol, sim.detect(sol, t0), cfg)
        rec = classify(trig, frames(sol, meter_for(device)), sol.dynamics.time, sol.dynamics.voltage(),
                       t0 if eps is not None else 0.0)
    return trig, rec, sol


# ==============================================================================
# scan
# ==============================================================================
def scan_task(task):
    from dataclasses import replace
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        g, cfg = task["geom"], replace(task["cfg"], save_every=SAVE_EVERY)
        ctx = sim.context(g, cfg.mesh_xi, task["film"], cfg.pset.params)
        eps = hsm.build_eps(ctx["hs"], cfg.hotspot, task["x_nm"], task["y_nm"], 1.0, cfg.tmax, cfg.pset.params)
        trig, rec, _ = observed_run(ctx["device"], task["I"], eps, cfg)
        return dict(ok=True, task=task, res=dict(trig=trig._asdict(), **rec), wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def scan_path(cfg, g):
    return cfg.path(f"scan_{g.key}", "csv")


def cmd_scan(cfg, a):
    sim.need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    g = geo.get(a.geom)
    i_dep = a.film.I_dep_uA(g.free_width_min_nm())
    path = scan_path(cfg, g)
    have = {(round(float(r["x_nm"]), 2), round(float(r["y_nm"]), 2), round(float(r["bias"]), 4))
            for r in store.load_rows(path)}
    sites = a.sites or g.hit_grid(a.grid)
    tasks = [dict(geom=g, cfg=cfg, film=a.film, x_nm=x, y_nm=y, bias=b, I=b * i_dep)
             for b in a.biases for x, y in sites if (round(x, 2), round(y, 2), round(b, 4)) not in have]
    sim.log(f"-- scan {g.key}: {len(sites)} sites x {len(a.biases)} biases, {len(tasks)} runs to do  "
            f"({cfg.tag}, I_dep = {i_dep:.2f} uA)")

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED ({t['x_nm']:g},{t['y_nm']:g}) bias {t['bias']}: {r['error']}")
            return
        res, trig = r["res"], r["res"]["trig"]
        store.append_row(path, dict(
            Geometry=g.key, x_nm=t["x_nm"], y_nm=t["y_nm"], bias=t["bias"], I_uA=round(t["I"], 4), tmax=cfg.tmax,
            hotspot=cfg.hotspot, mesh_xi=cfg.mesh_xi, outcome=res["outcome"], branch=trig["branch"],
            t_trigger="" if trig["t"] is None else round(trig["t"], 2), slips=round(res["slips"], 3),
            slips_last=round(res["slips_last"], 3),
            slips_settled="" if res["slips_settled"] is None else round(res["slips_settled"], 3), t_onset="" if res["t_onset"] is None else round(res["t_onset"], 1),
            n_end=fmt_n(res["n_end"]),
            n_abs=sum(abs(v) for v in res["n_end"]), t_trap="" if res["t_trap"] is None else round(res["t_trap"], 1),
            dwell=round(res["dwell"], 1), history=res["history"], wall_s=round(r["wall_s"], 1)), SCAN_FIELDS)

    sim.execute(tasks, a.workers, on_result, fn=scan_task)
    cmd_report(cfg, a)


# ==============================================================================
# Per-site summary of the scan
# ==============================================================================
def load_scan(cfg, g):
    by_site = {}
    for r in store.load_rows(scan_path(cfg, g)):
        by_site.setdefault((float(r["x_nm"]), float(r["y_nm"])), []).append(r)
    return {s: sorted(rs, key=lambda r: float(r["bias"])) for s, rs in by_site.items()}


def site_summary(rows, b_max):
    """Lowest bias for each kind of event; b_max + step stands for "not seen in the scanned range"."""
    bias = [float(r["bias"]) for r in rows]
    step = min(np.diff(bias)) if len(bias) > 1 else 0.01
    none = b_max + step

    def first(pred):
        return next((float(r["bias"]) for r in rows if pred(r)), none)

    fired = lambda r: r["outcome"] in FIRED
    sus = lambda r: r["outcome"] in SUSTAINED
    s = dict(
        i_trap=first(lambda r: r["outcome"] != "clean"),
        i_fire=first(fired),
        i_sust=first(sus),
        i_direct=first(lambda r: r["outcome"] == "direct"),
        i_prompt=first(lambda r: sus(r) and r["t_onset"] != "" and float(r["t_onset"]) <= PROMPT_T),
        none=none)
    below = [r for r in rows if float(r["bias"]) < s["i_sust"]]
    s["n_below"] = parse_n(below[-1]["n_end"]) if below else ()
    s["outcome_below"] = below[-1]["outcome"] if below else ""
    # non-monotonic: an event that is seen at some bias and not at a higher one
    seq = [r["outcome"] for r in rows]
    s["nonmono"] = any(seq[i] != "clean" and seq[j] == "clean" for i in range(len(seq)) for j in range(i + 1, len(seq))) \
        or any(sus(rows[i]) and not sus(rows[j]) for i in range(len(rows)) for j in range(i + 1, len(rows)))
    return s


def hitmap_ref(cfg, g, film, hitmap_dir=HITMAP_DIR):
    """{site: i_det in units of I_dep} from the hitmap results at the same window, mesh and hotspot."""
    hcfg = sim.with_out_dir(cfg, hitmap_dir)
    i_dep = film.I_dep_uA(g.free_width_min_nm())
    return {(float(r["x_nm"]), float(r["y_nm"])): float(r["I_sw_uA" if r["flag"] == "at_wall" else "I_det_uA"]) / i_dep
            for r in store.load_rows(hcfg.path("hitmap", "csv"))
            if r["Geometry"] == g.key and abs(float(r["energy_f"]) - 1.0) < 1e-6}


def trapped_states(scan, g, max_states):
    """The most common trapped states (by number of sub-threshold runs that end in them) and, for each, the run at
    the highest bias that leaves it -- the photon just below detection -- plus, for the most common state, the run at
    the lowest bias."""
    runs = [(tuple(parse_n(r["n_end"])), float(r["bias"]), s) for s, rows in scan.items() for r in rows
            if r["outcome"] == "trapped"]
    count = Counter(n for n, _, _ in runs)
    picks = []
    for k, (n, c) in enumerate(count.most_common(max_states)):
        rs = sorted((b, s) for m, b, s in runs if m == n)
        picks.append(dict(n=n, count=c, bias=rs[-1][0], site=rs[-1][1], role="highest"))
        if k == 0 and rs[0][0] < rs[-1][0]:
            picks.append(dict(n=n, count=c, bias=rs[0][0], site=rs[0][1], role="lowest"))
    return picks


# ==============================================================================
# trapped
# ==============================================================================
def trapped_task(task):
    """kind "state": regenerate one trapped state, then its switching current, a zero-bias hold, a hold just below
    its switching current and the photon threshold at PROBE_SITES with the trapped flux.  kind "pristine": the same
    photon thresholds in the device without flux."""
    from dataclasses import replace
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        g, film = task["geom"], task["film"]
        cfg = replace(task["cfg"], save_every=SAVE_EVERY)
        ctx = sim.context(g, cfg.mesh_xi, film, cfg.pset.params)
        dev, i_dep = ctx["device"], film.I_dep_uA(g.free_width_min_nm())
        tol = sim.scales(g, film)[0]
        if task["kind"] == "pristine":
            return dict(ok=True, task=task, res=second_photon(ctx, cfg, None, task["I_sw_pristine"], tol, i_dep),
                        wall_s=time.time() - t_start)
        if task["site"] is None:
            from snspd.fluxoid import winding_psi
            trig, rec, seed = observed_run(dev, 0.0, None, replace(cfg, tmax=PLACE_RELAX),
                                           phase=winding_psi(dev, task["n"]))
        else:
            x, y = task["site"]
            eps = hsm.build_eps(ctx["hs"], cfg.hotspot, x, y, 1.0, cfg.tmax, cfg.pset.params)
            trig, rec, seed = observed_run(dev, task["bias"] * i_dep, eps, cfg)
        res = dict(regenerated=fmt_n(rec["n_end"]), regen_outcome=rec["outcome"], regen_history=rec["history"])
        if rec["outcome"] != "trapped" or (task["site"] is None and tuple(rec["n_end"]) != tuple(task["n"])):
            return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
        lo, hi = 0.30 * i_dep, 1.05 * i_dep
        I_sw, flag, n_runs, tr = sim.bisect(dev, lo, hi, None, cfg, tol, seed=seed)
        res.update(i_sw=I_sw / i_dep, isw_flag=flag, isw_t_trigger=tr.t, isw_branch=tr.branch)
        for label, I in (("hold_zero", 0.0), ("hold_below", I_sw - tol), ("at_isw", I_sw)):
            _, h, _ = observed_run(dev, I, None, cfg, seed)
            res[label] = dict(i=I / i_dep, outcome=h["outcome"], n_end=fmt_n(h["n_end"]), history=h["history"],
                              slips=round(h["slips"], 2), slips_last=round(h["slips_last"], 2),
                              slips_settled=None if h["slips_settled"] is None else round(h["slips_settled"], 2))
        res["second"] = second_photon(ctx, cfg, seed, I_sw, tol, i_dep) if task["second"] else {}
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def second_photon(ctx, cfg, seed, I_top, tol, i_dep):
    """Photon threshold at each of PROBE_SITES, starting from `seed` (None: the pristine device), in units of I_dep;
    a site that never switches below I_top is given I_top."""
    from snspd import hotspot as hsm
    out = {}
    for xs, ys in PROBE_SITES:
        eps = hsm.build_eps(ctx["hs"], cfg.hotspot, xs, ys, 1.0, cfg.tmax, cfg.pset.params)
        I_det, fl, _, tr = sim.bisect(ctx["device"], cfg.hit_lo_frac * I_top, I_top, eps, cfg, tol, seed=seed)
        out[f"{xs:g},{ys:g}"] = dict(i_det=(I_top if fl == "at_wall" else I_det) / i_dep, flag=fl, t_trigger=tr.t)
    return out


def trapped_path(cfg, g):
    return cfg.path(f"trapped_{g.key}", "json")


def cmd_trapped(cfg, a):
    """States from the scan, or the ones given with --from (no scan needed, e.g. for a mesh check)."""
    sim.need_tdgl()
    g = geo.get(a.geom)
    if a.from_ or a.place:
        picks = [dict(n=(), count=0, site=tuple(float(v) for v in s.split("@")[0].split(",")),
                      bias=float(s.split("@")[1]), role="given") for s in a.from_]
        picks += [dict(n=parse_n(s), count=0, site=None, bias=0.0, role="placed") for s in a.place]
    else:
        scan = load_scan(cfg, g)
        if not scan:
            raise SystemExit(f"no scan for {g.key} at {cfg.tag}; run `fluxtrap.py scan` first")
        picks = trapped_states(scan, g, a.max_states)
    isw = store.load_json(sim.with_out_dir(cfg, a.hitmap_dir).path("isw", "json"), {}).get(g.key, {}).get("I_sw")
    if isw is None and not a.no_second:
        raise SystemExit(f"no I_sw for {g.key} at {cfg.tag} in {a.hitmap_dir}/ (needed for the second photon)")
    path = trapped_path(cfg, g)
    done = store.load_json(path, {})
    key = lambda p: (f"{fmt_n(p['n'])} @ placed" if p["site"] is None
                     else f"{fmt_n(p['n'])} @ ({p['site'][0]:g},{p['site'][1]:g}) {p['bias']:.3f}".strip())
    tasks = [dict(kind="state", geom=g, cfg=cfg, film=a.film, I_sw_pristine=isw, second=not a.no_second, **p)
             for p in picks if key(p) not in done]
    if "pristine" not in done and not a.no_second:
        tasks.insert(0, dict(kind="pristine", geom=g, cfg=cfg, film=a.film, I_sw_pristine=isw))
    sim.log(f"-- trapped {g.key}: {len(picks)} states, {len(tasks)} to do ({cfg.tag})")

    def on_result(r):
        t = r["task"]
        name = "pristine" if t["kind"] == "pristine" else key(t)
        if not r["ok"]:
            sim.log(f"   FAILED {name}: {r['error']}")
            return
        if t["kind"] == "pristine":
            done["pristine"] = r["res"]
            store.save_json(path, done)
            sim.log(f"   pristine second-photon thresholds done ({r['wall_s']:.0f} s)")
            return
        done[key(t)] = dict(n=fmt_n(t["n"]), count=t["count"], site=t["site"] and list(t["site"]), bias=t["bias"],
                            role=t["role"],
                            wall_s=round(r["wall_s"], 1), **r["res"])
        store.save_json(path, done)
        sim.log(f"   {key(t)}: i_sw = {r['res'].get('i_sw', float('nan')):.4f}  ({r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result, fn=trapped_task)
    cmd_report(cfg, a)


# ==============================================================================
# trace
# ==============================================================================
def trace_task(task):
    """One fixed-bias run with a photon: classification, frames, V(t) every 0.25 tau0 and the phase slipped between
    the probes since the photon, |int V dt| / 2 pi (one vortex crossing the full width = 1) (worker)."""
    from dataclasses import replace
    from snspd import hotspot as hsm
    t_start = time.time()
    try:
        g, cfg = task["geom"], replace(task["cfg"], save_every=SAVE_EVERY)
        ctx = sim.context(g, cfg.mesh_xi, task["film"], cfg.pset.params)
        eps = hsm.build_eps(ctx["hs"], cfg.hotspot, task["x_nm"], task["y_nm"], 1.0, cfg.tmax, cfg.pset.params)
        trig, rec, sol = observed_run(ctx["device"], task["I"], eps, cfg)
        t, v = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
        tg = np.arange(0.0, t[-1], 0.25)
        slips = np.interp(tg, t, sim.phase_slips(t, v, task["t0"]))
        return dict(ok=True, task=task, res=dict(trig=trig._asdict(), t=tg, v=np.abs(np.interp(tg, t, v)), slip_series=slips, **rec),
                    wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def cmd_trace(cfg, a):
    sim.need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    from snspd.hotspot import HOTSPOT
    g = geo.get(a.geom)
    x, y = (float(v) for v in a.site.split(","))
    i_dep = a.film.I_dep_uA(g.free_width_min_nm())
    tasks = [dict(geom=g, cfg=cfg, film=a.film, x_nm=x, y_nm=y, bias=b, I=b * i_dep, t0=HOTSPOT.t0_tau0)
             for b in a.biases]
    res = {}
    sim.execute(tasks, a.workers, lambda r: res.__setitem__(r["task"]["bias"], r), fn=trace_task)
    L = [f"TRACE {g.key} photon at ({x:g},{y:g}) nm, t0 = 20 tau0   {cfg.tag}   bias in units of I_dep = {i_dep:.2f} uA",
         "holes: " + ", ".join(f"{i}:({h.cx_nm:g},{h.cy_nm:g})" for i, h in enumerate(g.holes)),
         "per frame: t, fluxoid per hole, min |psi| on the fluxoid loops, max |V| since the previous frame, "
         "phase slips since the photon"]
    for b in sorted(res):
        r = res[b]
        if not r["ok"]:
            L.append(f"bias {b:.3f}  FAILED {r['error']}")
            continue
        q, tr = r["res"], r["res"]["trig"]
        L.append(f"bias {b:.3f}: {q['outcome']}" + (f", {tr['branch']} at t = {tr['t']:.1f}" if tr["branch"] else "")
                 + f", ends [{fmt_n(q['n_end'])}], {q['slips']:.2f} phase slips by t = {q['t'][-1]:.0f}")
        L += [f"  t = {tm:6.1f}   [{fmt_n(n)}]   {fl:.2f}   {vm:.3f}   {np.interp(tm, q['t'], q['slip_series']):.2f}"
              for tm, n, fl, vm in q["frames"]]
    text = "\n".join(L)
    sim.log(text)
    stem = f"trace_{g.key}_{x:g}_{y:g}"
    with open(cfg.path(stem, "txt"), "w") as fh:
        fh.write(text + "\n")
    ok = {b: r["res"] for b, r in res.items() if r["ok"]}
    if a.detailed:
        fig_trace(g, ok, cfg.path(stem + "_detailed", "png"), f"{g.label}: photon at ({x:g}, {y:g}) nm  ({cfg.tag})")
    fig_trace_simple(g, ok, cfg.path(stem, "png"), f"{g.label}, photon at ({x:g}, {y:g}) nm, {cfg.hotspot} hotspot")


# ==============================================================================
# compare: the same scan with two hotspot models
# ==============================================================================
def cmd_compare(cfg, a):
    """Conserving vs kinetic scan of one geometry: per-site thresholds, area shares, trapped states, timing."""
    from dataclasses import replace
    g = geo.get(a.geom)
    cfgs = {m: replace(cfg, hotspot=m) for m in ("conserving", "kinetic")}
    scans = {m: load_scan(c, g) for m, c in cfgs.items()}
    if not all(scans.values()):
        raise SystemExit(f"need both scans of {g.key} at {cfg.tmax:g} tau0, {cfg.mesh_xi:g} xi")
    common = sorted(set(scans["conserving"]) & set(scans["kinetic"]), key=lambda s: (s[1], s[0]))
    bias = {m: sorted({float(r["bias"]) for s in common for r in scans[m][s]}) for m in scans}
    # compare on the biases scanned in both
    shared = sorted(set(bias["conserving"]) & set(bias["kinetic"]))
    sub = {m: {s: [r for r in scans[m][s] if float(r["bias"]) in shared] for s in common} for m in scans}
    summ = {m: {s: site_summary(sub[m][s], shared[-1]) for s in common} for m in scans}
    fm = lambda v, s: f">{shared[-1]:.2f}" if v >= s["none"] - 1e-9 else f"{v:.2f}"
    L = [f"FLUX TRAPPING, CONSERVING vs KINETIC HOTSPOT  {g.label}  T{cfg.tmax:g} m{cfg.mesh_xi:g}  "
         f"({len(common)} sites, biases {shared[0]:.2f}..{shared[-1]:.2f} scanned in both, units of I_dep)", "",
         "Per site, lowest bias at which the photon leaves flux or fires (trap), fires the voltage trigger (fire), "
         "is detected (sust), is detected without a trapped dwell (direct); C = conserving, K = kinetic",
         f"{'site':>12}{'trap C':>8}{'trap K':>8}{'fire C':>8}{'fire K':>8}{'sust C':>8}{'sust K':>8}"
         f"{'direct C':>10}{'direct K':>10}   K: state just below sust"]
    for s in common:
        c, k = summ["conserving"][s], summ["kinetic"][s]
        L.append(f"{f'({s[0]:g},{s[1]:g})':>12}{fm(c['i_trap'], c):>8}{fm(k['i_trap'], k):>8}{fm(c['i_fire'], c):>8}"
                 f"{fm(k['i_fire'], k):>8}{fm(c['i_sust'], c):>8}{fm(k['i_sust'], k):>8}{fm(c['i_direct'], c):>10}"
                 f"{fm(k['i_direct'], k):>10}   {k['outcome_below']:<9} [{fmt_n(k['n_below'])}]")
    full = len(common) == len(g.hit_grid("full"))
    L += ["", "Share of the unit cell (thresholds interpolated as in hitmap.py" + ("" if full else
          "; SUBSET of the sites, indicative only") + "), C / K:",
          f"{'bias':>6}{'flux, silent':>16}{'pulse, heals':>16}{'via trap':>14}{'direct':>14}{'detected':>14}"]
    shares = {}
    for m in scans:
        maps = {k: an.unit_map(g, [(x, y, summ[m][(x, y)][k]) for x, y in common])[2]
                for k in ("i_trap", "i_fire", "i_sust", "i_direct")}
        fin = np.isfinite(maps["i_fire"])
        shares[m] = lambda k, b, maps=maps, fin=fin: (maps[k][fin] <= b + 1e-9).mean() * 100
    for b in [v for v in (0.66, 0.70, 0.74, 0.76, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.88, 0.90)
              if shared[0] <= v <= shared[-1] and v in shared]:
        cells = []
        for k1, k2 in (("i_trap", "i_fire"), ("i_fire", "i_sust"), ("i_sust", "i_direct"), ("i_direct", None),
                       ("i_sust", None)):
            v = [shares[m](k1, b) - (shares[m](k2, b) if k2 else 0.0) for m in ("conserving", "kinetic")]
            cells.append(f"{v[0]:5.0f} /{v[1]:4.0f}")
        L.append(f"{b:>6.2f}" + "".join(f"{c:>14}" if i else f"{c:>16}" for i, c in enumerate(cells)))
    L += ["", "Trapped states left by sub-threshold photons (fluxoid per hole: runs C / K, over the shared biases):"]
    cnt = {m: Counter(r["n_end"] for s in common for r in sub[m][s] if r["outcome"] == "trapped") for m in scans}
    for n in sorted(set(cnt["conserving"]) | set(cnt["kinetic"]), key=lambda n: -(cnt["conserving"][n] + cnt["kinetic"][n])):
        L.append(f"  [{n}]  {cnt['conserving'][n]} / {cnt['kinetic'][n]}")
    L += ["", "Runs by outcome over the shared biases (C / K):  " + ",  ".join(
        f"{o} {sum(r['outcome'] == o for s in common for r in sub['conserving'][s])} / "
        f"{sum(r['outcome'] == o for s in common for r in sub['kinetic'][s])}" for o in OUTCOMES)]

    def stats(vals, nd=0):
        v = np.array(vals, float)
        return "-" if not v.size else (f"median {np.median(v):.{nd}f}, range {v.min():.{nd}f}-{v.max():.{nd}f} "
                                       f"(n = {v.size})")

    t0 = 20.0
    L += ["", "Timing after the photon [tau0], C | K:"]
    for lab, pick in (("flux reaches a hole (trapped runs)", lambda r: r["outcome"] == "trapped" and r["t_trap"]),
                      ("resistive onset, direct", lambda r: r["outcome"] == "direct" and r["t_onset"]),
                      ("resistive onset, via trap", lambda r: r["outcome"] == "via_trap" and r["t_onset"])):
        col = "t_trap" if "hole" in lab else "t_onset"
        v = {m: [float(r[col]) - t0 for s in common for r in sub[m][s] if pick(r)] for m in scans}
        L.append(f"  {lab:<36} {stats(v['conserving'])}  |  {stats(v['kinetic'])}")
    sl = {m: [float(r["slips"]) for s in common for r in sub[m][s] if r["outcome"] == "transient"] for m in scans}
    L.append(f"  {'phase slips in a pulse that heals':<36} {stats(sl['conserving'], 1)}  |  {stats(sl['kinetic'], 1)}")
    text = "\n".join(L)
    sim.log(text)
    with open(os.path.join(cfg.out_dir, f"compare_{g.key}_T{cfg.tmax:g}_m{cfg.mesh_xi:g}.txt"), "w") as fh:
        fh.write(text + "\n")


# ==============================================================================
# report
# ==============================================================================
def cmd_report(cfg, a):
    g = geo.get(a.geom)
    i_dep = a.film.I_dep_uA(g.free_width_min_nm())
    isw = store.load_json(sim.with_out_dir(cfg, a.hitmap_dir).path("isw", "json"), {}).get(g.key, {}).get("I_sw")
    L = [f"FLUX TRAPPING  {g.label}  {cfg.tag}   biases in units of I_dep = {i_dep:.2f} uA "
         f"(narrowest cross-section {g.free_width_min_nm():.0f} nm)" + ("" if isw is None else f",  dark I_sw = {isw / i_dep:.3f}"),
         "holes: " + ", ".join(f"{i}:({h.cx_nm:g},{h.cy_nm:g})" for i, h in enumerate(g.holes))]
    scan = load_scan(cfg, g)
    if scan:
        L += report_scan(cfg, a, g, scan, i_dep)
        fig_outcomes(cfg, a, g)
    L += report_trapped(cfg, g)
    text = "\n".join(L)
    sim.log(text)
    os.makedirs(cfg.out_dir, exist_ok=True)
    with open(cfg.path(f"summary_{g.key}", "txt"), "w") as fh:
        fh.write(text + "\n")


def report_scan(cfg, a, g, scan, i_dep):
    """Scan tables (and the state-diagram and map figures)."""
    biases = sorted({float(r["bias"]) for rows in scan.values() for r in rows})
    b_max = biases[-1]
    summ = {s: site_summary(rows, b_max) for s, rows in scan.items()}
    ref = hitmap_ref(cfg, g, a.film, a.hitmap_dir)
    complete = all(len(rows) == len(biases) for rows in scan.values())
    L = [f"{len(scan)} sites x {len(biases)} biases {biases[0]:.2f}..{b_max:.2f}" + ("" if complete else "   (INCOMPLETE)"),
         f"outcome per run: clean / trapped (no trigger, flux left) / transient (voltage trigger, < {sim.SUSTAIN_SLIPS:g} "
         f"slip in the last {sim.SUSTAIN_WIN:g} tau0) / direct (sustained switch) / via_trap (sustained, after >= "
         f"{DWELL_MIN:g} tau0 trapped at |V| < {V_QUIET:g} before {ONSET_SLIPS:g} slips)",
         "",
         "Per site, lowest scanned bias at which the photon: leaves flux or fires the voltage trigger (i_trap); fires the "
         "voltage trigger (i_fire: the detection criterion until 29 Oct 2026); switches and stays resistive (i_sust: "
         "the detection criterion, compare 'hitmap'); does so without a trapped dwell (i_direct); does so with the "
         f"resistive onset by t = {PROMPT_T:g} tau0 (i_prompt).  '>' = not in the range.",
         f"{'site':>12}{'i_trap':>8}{'i_fire':>8}{'i_sust':>8}{'hitmap':>8}{'i_direct':>10}{'i_prompt':>10}   state just below i_sust"]
    fm = lambda v, s: f">{b_max:.2f}" if v >= s["none"] - 1e-9 else f"{v:.2f}"
    for site in sorted(summ, key=lambda s: (s[1], s[0])):
        s = summ[site]
        h = ref.get(site)
        L.append(f"{f'({site[0]:g},{site[1]:g})':>12}{fm(s['i_trap'], s):>8}{fm(s['i_fire'], s):>8}"
                 f"{fm(s['i_sust'], s):>8}{'' if h is None else f'{h:.3f}':>8}{fm(s['i_direct'], s):>10}"
                 f"{fm(s['i_prompt'], s):>10}   "
                 f"{s['outcome_below']:<8} [{fmt_n(s['n_below'])}]" + ("   non-monotonic" if s["nonmono"] else ""))

    # consistency with the hitmap: the scan's i_sust must be the first grid bias at or above the hitmap threshold
    # (the hitmap threshold h is the top of a bisection bracket: the true one lies in (h - tol, h])
    tol_i = sim.scales(g, a.film)[0] / i_dep
    off = []
    for site, s in summ.items():
        h = ref.get(site)
        if h is None or s["i_sust"] >= s["none"] - 1e-9:
            continue
        k = biases.index(s["i_sust"])
        prev = biases[k - 1] if k else -np.inf
        if not (s["i_sust"] >= h - tol_i - 1e-9 and prev < h + 1e-9):
            off.append(f"({site[0]:g},{site[1]:g}) scan {s['i_sust']:.2f} vs hitmap {h:.3f}")
    L += ["", f"Check against the hitmap thresholds: {len(summ) - len(off)} of {len(summ)} sites switch first at the "
              f"grid bias just above the hitmap threshold" + ("" if not off else ";  differ: " + ", ".join(off))]

    # area fractions from the threshold maps
    maps = {}
    for k in ("i_trap", "i_fire", "i_sust", "i_direct", "i_prompt"):
        X, Y, Z = an.unit_map(g, [(x, y, summ[(x, y)][k]) for x, y in summ])
        maps[k] = Z
    fin = np.isfinite(maps["i_fire"])
    subset = len(scan) < len(g.hit_grid("full"))
    L += ["", "Area fractions of the unit cell (thresholds interpolated as in hitmap.py; grid resolution "
              f"{min(np.diff(biases)):.2f}" + ("; SUBSET of the sites, indicative only" if subset else "") + "):",
          "  fired = voltage trigger; sustained = stays resistive (detected); silent trap = flux left, no trigger",
          f"{'bias':>6}{'fired':>10}{'sustained':>11}{'direct':>9}{'prompt':>9}{'silent trap':>13}"]
    frac = lambda k, b: (maps[k][fin] <= b + 1e-9).mean() * 100
    for b in [v for v in (0.70, 0.75, 0.78, 0.80, 0.82, 0.85, 0.88, 0.90, 0.95) if biases[0] <= v <= b_max]:
        sil = ((maps["i_trap"][fin] <= b + 1e-9) & (maps["i_fire"][fin] > b + 1e-9)).mean() * 100
        L.append(f"{b:>6.2f}{frac('i_fire', b):>9.0f}%{frac('i_sust', b):>10.0f}%{frac('i_direct', b):>8.0f}%"
                 f"{frac('i_prompt', b):>8.0f}%{sil:>12.0f}%")
    L += ["", "Efficiency-curve points for each definition of a detection (grid biases: up to one step above the "
              "bisected thresholds):", f"{'':>11}{'i_10':>7}{'i_50':>7}{'i_95':>7}{'i_max':>7}"]
    for k, lab in (("i_fire", "fired"), ("i_sust", "sustained"), ("i_prompt", "prompt"), ("i_direct", "direct")):
        m = an.metrics(an.pixel_values(maps[k]))
        L.append(f"{lab:>11}{m['i_10']:>7.3f}{m['i_50']:>7.3f}{m['i_95']:>7.3f}{m['i_max']:>7.3f}")

    states = Counter(r["n_end"] for rows in scan.values() for r in rows if r["outcome"] == "trapped")
    L += ["", "Trapped states left by sub-threshold photons (fluxoid per hole, runs):"]
    L += [f"  [{n}]  {c}" for n, c in states.most_common()]
    if a.detailed:
        fig_state(scan, summ, ref, biases, g, cfg.path(f"states_{g.key}", "png"),
                  f"{g.label}: outcome of one photon  ({cfg.tag})")
        fig_maps(g, maps, summ, cfg.path(f"maps_{g.key}", "png"), f"{g.label}: thresholds in units of I_dep  ({cfg.tag})")
    return L


def report_trapped(cfg, g):
    L = []
    tr = store.load_json(trapped_path(cfg, g), {})
    pristine = tr.pop("pristine", {})
    if tr:
        L += ["", "Trapped states as initial conditions (no photon unless stated); i in units of I_dep:"]
        for k, v in tr.items():
            if "i_sw" not in v:
                L.append(f"  {k}: regenerated as {v['regen_outcome']} [{v['regenerated']}] -- skipped")
                continue
            src = (f"placed, relaxed {PLACE_RELAX:g} tau0 at zero bias" if v["site"] is None else
                   f"from a photon at ({v['site'][0]:g},{v['site'][1]:g}) at {v['bias']:.3f}")
            L.append(f"  [{v['n'] or v['regenerated']}] {src}  (regenerated [{v['regenerated']}], {v['count']} scan "
                     f"runs end in it)")
            L.append(f"      switching current of the trapped state  i_sw = {v['i_sw']:.4f}  {v['isw_flag']}  "
                     f"[{v['isw_branch']} at t = {v['isw_t_trigger']}]")
            for lab in ("at_isw", "hold_below", "hold_zero"):
                h = v.get(lab)
                if h:
                    # (total slips are not shown: from a seed they include the phase relaxing to the new current)
                    L.append(f"      held at i = {h['i']:.3f} for {cfg.tmax:g} tau0: {h['outcome']}, "
                             f"{h.get('slips_last', '?')} slips in the last {sim.SUSTAIN_WIN:g}"
                             + ("" if h.get("slips_settled") is None else f" ({h['slips_settled']} after settling)")
                             + f", ends [{h['n_end']}]   "
                             f"{h['history']}")
            if v["second"]:
                L.append(f"      second photon:   {'site':>10}{'with trap':>11}{'pristine':>10}{'change':>9}")
            for site, d in v["second"].items():
                a_, p_ = d["i_det"], pristine.get(site, {}).get("i_det")
                L.append(f"                       {f'({site})':>10}{a_:>11.4f}"
                         + ("" if p_ is None else f"{p_:>10.4f}{a_ - p_:>+9.4f}"))
    return L


# ==============================================================================
# Figures
# ==============================================================================
def fig_state(scan, summ, ref, biases, g, path, title):
    from matplotlib.colors import ListedColormap
    from snspd import plotting as pl
    import matplotlib.pyplot as plt
    pl.apply_style()
    sites = sorted(scan, key=lambda s: (summ[s]["i_sust"], summ[s]["i_trap"], s[1], s[0]))
    colors = ["#eeede9", "#9cc3ee", "#e8c547", "#eb6834", "#7a2d0d"]
    code = {o: i for i, o in enumerate(OUTCOMES)}
    Z = np.full((len(sites), len(biases)), np.nan)
    for i, s in enumerate(sites):
        for r in scan[s]:
            Z[i, biases.index(float(r["bias"]))] = code[r["outcome"]]
    edges = np.concatenate([[biases[0] - (biases[1] - biases[0]) / 2],
                            (np.array(biases[1:]) + np.array(biases[:-1])) / 2,
                            [biases[-1] + (biases[-1] - biases[-2]) / 2]])
    fig, ax = plt.subplots(figsize=(10, 0.22 * len(sites) + 1.8), layout="constrained")
    ax.pcolormesh(edges, np.arange(len(sites) + 1), Z, cmap=ListedColormap(colors), vmin=-0.5, vmax=len(colors) - 0.5)
    for i, s in enumerate(sites):
        h = ref.get(s)
        if h is not None:
            ax.plot([h, h], [i, i + 1], color=pl.INK, lw=1.4)
        for r in scan[s]:
            if r["outcome"] == "trapped":
                ax.text(float(r["bias"]), i + 0.5, str(r["n_abs"]), ha="center", va="center", fontsize=5.5, color=pl.INK)
    ax.set_yticks(np.arange(len(sites)) + 0.5, [f"({x:g},{y:g})" for x, y in sites], fontsize=6.5)
    ax.set_xlabel("bias  I / I_dep   (black tick: hitmap threshold; digits: flux quanta left in the holes)")
    ax.grid(False)
    for c, o in zip(colors, OUTCOMES):
        ax.plot([], [], "s", color=c, ms=8, label=o.replace("_", " "))
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    ax.set_title(title)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_maps(g, maps, summ, path, title):
    from snspd import plotting as pl
    import matplotlib.pyplot as plt
    pl.apply_style()
    X, Y = an.map_grid(g)
    labels = {"i_trap": "photon leaves flux or switches", "i_sust": "switches and stays resistive",
              "i_direct": "same, without a trapped dwell"}
    allv = np.concatenate([maps[k][np.isfinite(maps[k])] for k in labels])
    vmin, vmax = np.floor(allv.min() * 50) / 50, np.ceil(allv.max() * 50) / 50
    fig, axes = plt.subplots(1, len(labels), figsize=(5.2 * len(labels), 4.6), sharey=True, layout="constrained")
    im = None
    for ax, (k, lab) in zip(axes, labels.items()):
        im = ax.pcolormesh(X, Y, maps[k], cmap=pl.SEQ_CMAP, vmin=vmin, vmax=vmax, shading="auto")
        pl.draw_holes(ax, g)
        for (x, y), s in summ.items():
            ax.plot(x, y, "o", ms=3.5, mfc="white", mec=pl.INK, mew=0.7, zorder=4)
        ax.set_xlim(-g.width_nm / 2, g.width_nm / 2)
        ax.set_ylim(-g.cell_half_y_nm, g.cell_half_y_nm)
        ax.set_aspect("equal")
        ax.set_title(f"{k}: {lab}", fontsize=9.5)
        ax.set_xlabel("x0 [nm]")
        ax.grid(False)
    axes[0].set_ylabel("y0 [nm]  (one unit cell)")
    fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02).set_label("bias  I / I_dep")
    fig.suptitle(title)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_trace(g, res, path, title):
    from snspd import plotting as pl
    import matplotlib.pyplot as plt
    from snspd.fluxoid import PSI_LOOP_MIN
    pl.apply_style()
    bs = sorted(res)
    fig, axes = plt.subplots(3, len(bs), figsize=(4.2 * len(bs), 7.6), sharex=True, squeeze=False, layout="constrained")
    cols = ["#2a78d6", "#eb6834", "#1f9e74", "#a347c2", "#c49a0e", "#52514e", "#d6336c"]
    for k, b in enumerate(bs):
        q = res[b]
        av, an_, asl = axes[0][k], axes[1][k], axes[2][k]
        asl.plot(q["t"], q["slip_series"], color=pl.INK, lw=1.2)
        av.plot(q["t"], q["v"], color=pl.INK, lw=1.0)
        av.axhline(sim.V_POST_MAX, color=pl.INK2, ls=":", lw=1.0)
        av.set_yscale("symlog", linthresh=0.01)
        av.set_ylim(0, 60)
        tr = q["trig"]
        av.set_title(f"I = {b:.3f} I_dep: {q['outcome'].replace('_', ' ')}"
                     + (f" ({tr['branch']} at t = {tr['t']:.0f})" if tr["branch"] else ""), fontsize=9.5)
        fr = [f for f in q["frames"] if f[2] >= PSI_LOOP_MIN]
        t = np.array([f[0] for f in fr])
        n = np.array([f[1] for f in fr])
        for h in range(n.shape[1] if n.size else 0):
            if n[:, h].any():
                c = g.holes[h]
                an_.step(t, n[:, h], where="post", color=cols[h % len(cols)], label=f"hole {h} ({c.cx_nm:g},{c.cy_nm:g})")
        if tr["branch"]:
            for ax in (av, an_, asl):
                ax.axvline(tr["t"], color=pl.INK2, ls="--", lw=1.0)
        asl.set_xlabel("t [tau0]   (photon at 20)")
        if an_.lines:
            an_.legend(fontsize=7.5, loc="best")
    axes[0][0].set_ylabel("|V| [V0]   (dotted: V_peak threshold)")
    axes[1][0].set_ylabel("fluxoid [flux quanta]")
    axes[2][0].set_ylabel("phase slips since the photon\n|int V dt| / 2 pi")
    fig.suptitle(title)
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ---- simple figures (the default since 5 Nov 2026) ---------------------------------------------------------------
# Bands from the bottom.  Colours checked for adjacent-pair separation (OKLab, normal and simulated protan / deutan /
# tritan vision: >= 16 for every pair).
OUTCOME_BANDS = (("detected directly", "#1d4e89", "white"),
                 ("detected via trapped flux", "#c84b17", "white"),
                 ("pulse, then heals", "#e3a832", None),
                 ("flux left, no signal", "#99bfeb", None),
                 ("nothing", "#f1f0eb", None))
OUTCOME_WORDS = {"clean": "nothing", "trapped": "flux left", "transient": "pulse, heals", "direct": "detected",
                 "via_trap": "detected via trapped flux"}


def outcome_bands(cfg, g):
    """(grid biases, share of the unit cell in each OUTCOME_BANDS band, subset?) from the scan at cfg, or None.  Same
    threshold maps and area rule as the report tables."""
    scan = load_scan(cfg, g)
    if not scan:
        return None
    biases = sorted({float(r["bias"]) for rows in scan.values() for r in rows})
    summ = {s: site_summary(rows, biases[-1]) for s, rows in scan.items()}
    maps = {k: an.unit_map(g, [(x, y, summ[(x, y)][k]) for x, y in summ])[2]
            for k in ("i_trap", "i_fire", "i_sust", "i_direct")}
    fin = np.isfinite(maps["i_fire"])
    F = {k: np.array([(m[fin] <= b + 1e-9).mean() for b in biases]) for k, m in maps.items()}
    bands = np.array([F["i_direct"], F["i_sust"] - F["i_direct"], F["i_fire"] - F["i_sust"],
                      F["i_trap"] - F["i_fire"], 1.0 - F["i_trap"]])
    return np.array(biases), np.clip(bands, 0.0, 1.0), len(scan) < len(g.hit_grid("full"))


def plain_detected(cfg, g, a):
    """(bias, detected share) of the plain reference wire from the hitmap at the same hotspot, or None."""
    if not g.holes:
        return None
    h = geo.get("Plain_Wire_2um" if g.length_nm == geo.LONG_LENGTH_NM else "Plain_Wire")
    ref = hitmap_ref(cfg, h, a.film, a.hitmap_dir)
    if not ref:
        return None
    Z = an.unit_map(h, [(x, y, v) for (x, y), v in ref.items()])[2]
    Z = Z[np.isfinite(Z)]
    b = np.linspace(0.5, 1.0, 1001)
    return b, np.array([(Z <= v).mean() for v in b]), h


def fig_outcomes(cfg, a, g):
    """What one photon does, as a share of the unit cell, vs bias: one panel per hotspot model with a scan."""
    from dataclasses import replace
    from snspd import plotting as pl
    import matplotlib.pyplot as plt
    pl.apply_style()
    models = [(m, replace(cfg, hotspot=m)) for m in ("conserving", "kinetic")]
    data = [(m, c, outcome_bands(c, g)) for m, c in models]
    data = [(m, c, d) for m, c, d in data if d is not None]
    if not data:
        return
    # common x range: from the first bias with any event to the first at which everything is detected
    first = min(b[0] for _, _, (b, _, _) in data)
    lo = max(first, min(b[np.argmax(bands[4] < 0.995)] for _, _, (b, bands, _) in data) - 0.03)
    hi = max(b[min(np.argmax(bands[0] + bands[1] > 0.995), len(b) - 1)] for _, _, (b, bands, _) in data) + 0.02
    hi = min(hi, min(b[-1] for _, _, (b, _, _) in data))
    width = max(8.0, 5.4 * len(data))
    fig, axes = plt.subplots(1, len(data), figsize=(width, 4.0), sharey=True, layout="constrained", squeeze=False)
    ax_in = 0.82 * width / len(data)                                 # approximate axes width [inch]
    for ax, (m, c, (b, bands, subset)) in zip(axes[0], data):
        ax.stackplot(b, bands * 100, colors=[col for _, col, _ in OUTCOME_BANDS], edgecolor=pl.SURFACE, linewidth=1.0)
        if b[0] > lo + 1e-9:                                         # this model's scan starts later
            ax.axvspan(lo, b[0], facecolor="none", hatch="///", edgecolor=pl.INK2, lw=0, alpha=0.3)
            ax.text((lo + b[0]) / 2, 50, "not scanned", ha="center", va="center", fontsize=8, color=pl.INK2,
                    rotation=90)
        base = np.cumsum(bands, axis=0) - bands / 2
        ref = plain_detected(c, g, a)
        cum = np.cumsum(bands, axis=0)
        xg = np.linspace(lo, hi, 241)
        for k, (lab, _, tc) in enumerate(OUTCOME_BANDS[:4]):
            half = 0.5 * len(lab) * 0.062 / ax_in * (hi - lo)       # half the label width in bias units (8 pt)
            bot, top = cum[k] - bands[k], cum[k]
            th = np.interp(xg, b, bands[k])
            for x in xg[np.argsort(-th)][:60]:                       # thickest places first; label only if it fits
                if th[np.argmin(np.abs(xg - x))] < 0.10 or x - half < lo or x + half > hi:
                    continue
                y = np.interp(x, b, base[k])
                xs = np.linspace(x - half, x + half, 9)
                if ref is not None and abs(np.interp(y, ref[1], ref[0]) - x) < half + 0.004:
                    continue                                         # keep labels off the reference line
                if np.all(np.interp(xs, b, bot) < y - 0.035) and np.all(np.interp(xs, b, top) > y + 0.035):
                    ax.text(x, y * 100, lab, ha="center", va="center", fontsize=8, color=tc or pl.INK)
                    break
        if ref is not None:
            rb, rf, h = ref
            ax.plot(rb, rf * 100, color=pl.INK, lw=1.6)
            j = int(np.argmax(rf >= 0.97))
            ax.annotate(f"{h.label}: detected", (rb[j], 97), xytext=(-6, 0), textcoords="offset points",
                        ha="right", va="center", fontsize=8, color=pl.INK)
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, 100)
        ax.set_xlabel("bias  I / I_dep")
        ax.set_title(f"{m} hotspot" + ("  (subset of sites: indicative)" if subset else ""), fontsize=10)
        ax.grid(False)
    axes[0][0].set_ylabel("share of the unit cell [%]")
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=col, edgecolor=pl.INK2, lw=0.4) for _, col, _ in OUTCOME_BANDS]
    fig.legend(handles, [lab for lab, _, _ in OUTCOME_BANDS], loc="outside lower center", ncol=5 if len(data) > 1 else 3, frameon=False,
               fontsize=8.5)
    fig.suptitle(f"{g.label}: what one photon does ({cfg.tmax:g} tau0 window, {cfg.mesh_xi:g} xi mesh)", fontsize=11)
    fig.savefig(os.path.join(cfg.out_dir, f"outcomes_{g.key}_T{cfg.tmax:g}_m{cfg.mesh_xi:g}.png"), dpi=200)
    plt.close(fig)


def fig_trace_simple(g, res, path, title):
    """Two panels: flux quanta held by the holes and phase slips since the photon, one line per bias."""
    from snspd import plotting as pl
    import matplotlib.pyplot as plt
    from snspd.fluxoid import PSI_LOOP_MIN
    pl.apply_style()
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#7a5bc2"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.7), sharex=True, layout="constrained")
    for b, col in zip(sorted(res), cols):
        q = res[b]
        fr = [f for f in q["frames"] if f[2] >= PSI_LOOP_MIN]
        a1.step([f[0] for f in fr], [sum(abs(v) for v in f[1]) for f in fr], where="post", color=col, lw=1.8,
                label=f"{b:.2f} I_dep: {OUTCOME_WORDS[q['outcome']]}, {q['slips']:.1f} slips by t = {q['t'][-1]:.0f}")
        a2.plot(q["t"], q["slip_series"], color=col, lw=1.8)
    for ax in (a1, a2):
        ax.axvline(20.0, color=pl.INK2, ls=":", lw=1.0)
        ax.set_xlabel("t [tau0]   (photon at 20, dotted)")
    a1.set_title("flux quanta held by the holes (a line ends where the loops go normal)", fontsize=9)
    a2.set_title("phase slips since the photon (|int V dt| / 2 pi)", fontsize=9)
    a2.set_yscale("symlog", linthresh=1.0)
    a2.set_ylim(0, None)
    fig.legend(*a1.get_legend_handles_labels(), loc="outside lower center", ncol=2, fontsize=8.5, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["scan", "trapped", "report", "trace", "compare"])
    ap.add_argument("--geom", required=True, help="a registered holed geometry (snspd.geometry.GEOMETRIES)")
    ap.add_argument("--grid", choices=["full", "critical"], default="full")
    ap.add_argument("--biases", default="0.56:0.70:0.02,0.71:0.96:0.01", help="units of I_dep (nothing happens below ~0.64)")
    ap.add_argument("--sites", default="", help="scan: only these sites, x,y;x,y [nm] (default: the whole --grid)")
    ap.add_argument("--hitmap-dir", default=HITMAP_DIR, help="hitmap results to compare with")
    ap.add_argument("--site", default="-68,0", help="trace: photon site x,y [nm] (use --site=-68,0)")
    ap.add_argument("--max-states", type=int, default=5, help="trapped: most common trapped states to study")
    ap.add_argument("--from", dest="from_", action="append", default=[], metavar="X,Y@BIAS",
                    help="trapped: generate the state with a photon at (X, Y) nm at BIAS instead of taking it from "
                         "the scan (repeatable; use --from=-68,0@0.77)")
    ap.add_argument("--place", action="append", default=[], metavar="N0 N1 ...",
                    help="trapped: place this fluxoid state (one integer per hole) instead of generating it with a "
                         "photon (repeatable; use --place=\"0 -1 0\")")
    ap.add_argument("--no-second", action="store_true", help="trapped: skip the second-photon thresholds")
    ap.add_argument("--hotspot", choices=["conserving", "legacy", "kinetic"], default=sim.SimConfig.hotspot)
    ap.add_argument("--tmax", type=float, default=400.0)
    ap.add_argument("--mesh-xi", type=float, default=sim.SimConfig.mesh_xi)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--preset", choices=list(sim.PRESETS), default="",
                    help="parameter set (snspd.sim.PRESETS): '' = gamma 1 results, 'nbtin' = gamma 10 + NbTiN times")
    ap.add_argument("--detailed", action="store_true",
                    help="report/trace: also the detailed figures (per-site state diagram, threshold maps, 3-row trace)")
    a = ap.parse_args(argv)
    a.film = sim.PRESETS[a.preset].film
    a.biases = parse_biases(a.biases)
    a.sites = [tuple(float(v) for v in s.split(",")) for s in a.sites.split(";") if s]
    cfg = sim.SimConfig(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi, out_dir=a.out, preset=a.preset)
    {"scan": cmd_scan, "trapped": cmd_trapped, "report": cmd_report, "trace": cmd_trace,
     "compare": cmd_compare}[a.command](cfg, a)


if __name__ == "__main__":
    main()
