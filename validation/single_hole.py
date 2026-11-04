"""
2c  One hole in the plain 200 nm wire (production mesh, window and criterion): the step between the plain wire and the
zipper.

  phase     A vortex that enters from the left edge and stops in a hole shifts the phase between the probes by
            Delta theta measured along a path on the far (right) side of the hole, which the vortex never crosses.
            The plan's estimate is x/w (a point vortex at distance x from the entry edge, w << Lambda): 0.32 for the
            zipper's left-centre hole against ~0.4 measured in the flux scan.  For a hole of finite size the exact
            London-limit statement is a reciprocity identity: Delta theta / 2 pi = I_L / I, the fraction of a
            transport current that passes between the entry edge and the hole (Green's identity for the two harmonic
            phase fields; it holds with other holes present).  Both sides are measured at zero bias for 20 nm circles
            at five x, the zipper's 16 x 110 nm ellipse alone and the zipper's hole 1 in the zipper; the phase per
            quantum is also measured at 0.70 I_dep, the bias of the scan runs, and compared with the scan.
            Channel currents are computed from psi (common.supercurrent_through): pyTDGL's interpolated current
            density reads up to ~10 % high in channels a few xi wide (see 2b).
  hold      one quantum in a centred 40 nm hole, zero bias, 400 tau0: it stays, with a circulating current
  I_sw      the same hole with and without the quantum.  Plan: "earlier by roughly the circulating current".  A linear
            two-channel picture gives 2 I_c; with the GL current-phase relation in each channel and the fluxoid
            fixing their phase difference, the drop is quadratic in I_c (common.two_channel_isw).
"""
from __future__ import annotations

import time
from dataclasses import replace

import numpy as np

from snspd import geometry as geo
from snspd import sim, store
from snspd.fluxoid import PSI_LOOP_MIN, FluxoidMeter, winding_psi

from . import common as C
from .disc_vs_hole import BRACKET, isw_values

ID, TITLE = "2c", "one hole in the plain strip"
W, L = 200.0, 800.0
RELAX = 50.0                       # zero-bias relaxation of a placed quantum [tau0] (fluxtrap.PLACE_RELAX)
I_LOW = 0.1                        # transport bias for the current split [I_dep of the free width]
I_SCAN = 0.70                      # bias for the phase per quantum under current [I_dep of the free width]
T_BIAS = 100.0                     # settling time at that bias [tau0]
R_SW = 20.0                        # the switching-current hole (shared with 1c)
SCAN = "results/fluxtrap/scan_Long_Zipper_conserving_T400_m1.4.csv"
TRAPPED = "results/fluxtrap/trapped_Long_Zipper_conserving_T400_m1.4.json"
ZIP = "zipper hole 1 (-36, 0) in the zipper"


def cases():
    out = [(f"circle d = 20 nm at x = {x:g}", C.holed_strip(x, 10.0), 0) for x in (-72.0, -36.0, 0.0, 36.0, 72.0)]
    out.append(("zipper ellipse 16 x 110 nm alone at x = -36", C.holed_strip(-36.0, 8.0, 55.0), 0))
    out.append((ZIP, geo.LONG_ZIPPER, 1))
    return out


def cuts(h):
    """Left and right cross-section cuts at the hole's centre height, and the probe-to-probe paths on either side
    of it (nm)."""
    yp, e = L / 3, 0.05                                                      # e: keep the cuts off the walls
    xl, xr = (h.cx_nm - h.a_nm - W / 2) / 2, (h.cx_nm + h.a_nm + W / 2) / 2
    xr = 0.0 if h.cx_nm + h.a_nm < 0 else xr                                 # the probes' own line if it is clear
    return (C.polyline([(-W / 2 + e, h.cy_nm), (h.cx_nm - h.a_nm - e, h.cy_nm)]),
            C.polyline([(h.cx_nm + h.a_nm + e, h.cy_nm), (W / 2 - e, h.cy_nm)]),
            C.polyline([(0, -yp), (xl, -yp), (xl, yp), (0, yp)]),
            C.polyline([(0, -yp), (xr, -yp), (xr, yp), (0, yp)]))


def phase_task(task):
    t_start = time.time()
    try:
        g, k = task["geom"], task["hole"]
        cfg = replace(C.CFG, save_every=40)
        dev = sim.context(g, cfg.mesh_xi)["device"]
        cut_l, cut_r, path_l, path_r = (np.asarray(p) for p in cuts(g.holes[k]))
        op_l, op_r = C.interp_op(dev, path_l), C.interp_op(dev, path_r)
        i_free = sim.FILM.I_dep_uA(g.free_width_min_nm())
        Is = lambda psi, cut: C.supercurrent_through(dev, psi, cut)
        res = {}
        for sign in (1, -1):                                   # +1 at zero bias; -1 = entry from the left at I > 0
            n = [0] * len(g.holes)
            n[k] = sign
            with sim.solved(dev, 0.0, None, replace(cfg, tmax=RELAX), phase=winding_psi(dev, n)) as sol:
                psi = np.asarray(sol.tdgl_data.psi)
                if sign == 1:
                    res["n_relaxed"] = FluxoidMeter(dev)(psi)[0]
                    (th_l, fl_l), (th_r, fl_r) = C.phase_change(op_l, psi), C.phase_change(op_r, psi)
                    res.update(th_l=th_l, th_r=th_r, floor=min(fl_l, fl_r), Ic_l=Is(psi, cut_l), Ic_r=Is(psi, cut_r),
                               Ic_l_tdgl=C.current_through(sol, cut_l), Ic_r_tdgl=C.current_through(sol, cut_r))
                else:
                    seed = sol
        for label, sd in (("with", seed), ("without", None)):
            with sim.solved(dev, I_SCAN * i_free, None, replace(cfg, tmax=T_BIAS), seed=sd) as sol:
                psi = np.asarray(sol.tdgl_data.psi)
                res[f"th_r_bias_{label}"] = C.phase_change(op_r, psi)[0]
                res[f"n_bias_{label}"] = FluxoidMeter(dev)(psi)[0][k]
        with sim.solved(dev, I_LOW * i_free, None, replace(cfg, tmax=sim.CROWD_TIME)) as sol:
            psi = np.asarray(sol.tdgl_data.psi)
            res.update(It_l=Is(psi, cut_l), It_r=Is(psi, cut_r), I=I_LOW * i_free,
                       It_l_tdgl=C.current_through(sol, cut_l), It_r_tdgl=C.current_through(sol, cut_r))
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def trapped_task(task):
    """The centred R_SW hole with one quantum: zero-bias hold with the fluxoid on every frame, circulating current,
    the switching current seeded from the relaxed state, and the fluxoid just below it."""
    t_start = time.time()
    try:
        import fluxtrap
        g = C.holed_strip(0.0, R_SW)
        cfg = C.CFG
        dev = sim.context(g, cfg.mesh_xi)["device"]
        cut_l, cut_r, _, _ = (np.asarray(p) for p in cuts(g.holes[0]))
        _, _, seed = fluxtrap.observed_run(dev, 0.0, None, replace(cfg, tmax=RELAX), phase=winding_psi(dev, [1]))
        _, hold, sol = fluxtrap.observed_run(dev, 0.0, None, replace(cfg, save_every=fluxtrap.SAVE_EVERY), seed)
        ns = [n[0] for _, n, fl, _ in hold["frames"] if fl >= PSI_LOOP_MIN]
        psi = np.asarray(sol.tdgl_data.psi)
        res = dict(hold_n=sorted(set(ns)), hold_frames=len(ns), hold_rejected=len(hold["frames"]) - len(ns),
                   Ic_l=C.supercurrent_through(dev, psi, cut_l), Ic_r=C.supercurrent_through(dev, psi, cut_r))
        i_free = sim.FILM.I_dep_uA(W - 2 * R_SW)
        tol = sim.scales(geo.PLAIN_WIRE)[0]
        I, flag, nr, trig = C.bisect(lambda I: sim.run_simulation(dev, I, None, cfg, seed=seed),
                                     BRACKET[0] * i_free, BRACKET[1] * i_free, tol)
        _, below, _ = fluxtrap.observed_run(dev, I - tol, None, replace(cfg, save_every=fluxtrap.SAVE_EVERY), seed)
        res.update(I_sw=I, flag=flag, n_runs=nr, below_outcome=below["outcome"], below_n_end=list(below["n_end"]),
                   below_history=below["history"], below_slips=below["slips"])
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def scan_slips(bias):
    """Phase slips of the zipper scan runs at this bias that end in [0 -1 0 0 0 0 0] (one quantum in hole 1)."""
    s = [float(r["slips"]) for r in store.load_rows(SCAN) if r["outcome"] == "trapped" and r["n_end"] == "0 -1 0 0 0 0 0"
         and abs(float(r["bias"]) - bias) < 1e-6]
    return (float(np.median(s)), float(np.percentile(s, 10)), float(np.percentile(s, 90)), len(s)) if s else None


def baseline_task(task):
    """Phase slipped between the probes after t0 in a dark run of the zipper at I_SCAN: what every scan run counts
    on top of any flux (the bias transient is still settling at t0)."""
    t_start = time.time()
    try:
        g = geo.LONG_ZIPPER
        dev = sim.context(g, C.CFG.mesh_xi)["device"]
        with sim.solved(dev, I_SCAN * sim.FILM.I_dep_uA(g.free_width_min_nm()), None, C.CFG) as sol:
            t, v = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
        s = sim.phase_slips(t, v, sim.hs_mod.HOTSPOT.t0_tau0)
        return dict(ok=True, task=task, res=dict(slips=float(s[-1]), slips_by_100=float(np.interp(100.0, t, s))),
                    wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def task_fn(task):
    return {"phase": phase_task, "trapped": trapped_task, "baseline": baseline_task}[task["fn"]](task)


def run(a):
    prev = C.load_groups().get(ID, {}).get("extra", {})                  # finished solver results are reused
    out = dict(prev.get("phase", {}))
    out.update({k: prev[k] for k in ("trapped", "baseline") if prev.get(k)})

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['name']}: {r['error']}")
            return
        out[t["name"]] = r["res"]
        sim.log(f"   {t['name']}: done ({r['wall_s']:.0f} s)")

    tasks = [dict(name=n, geom=g, hole=k, fn="phase") for n, g, k in cases()]
    tasks += [dict(name="trapped", fn="trapped"), dict(name="baseline", fn="baseline")]
    tasks = [t for t in tasks if t["name"] not in out]
    isw0 = isw_values((("hole", R_SW, None),), a.workers).get(("hole", R_SW, None))      # 1c's run, if done
    sim.execute(tasks, a.workers, on_result, fn=task_fn)
    checks = []
    extra = dict(phase={k: v for k, v in out.items() if k not in ("trapped", "baseline")}, trapped=out.get("trapped"),
                 baseline=out.get("baseline"))

    rows = []
    for name, g, k in cases():
        q = out.get(name)
        if q is None:
            continue
        h = g.holes[k]
        f_cur = q["It_l"] / (q["It_l"] + q["It_r"])
        f_tdgl = q["It_l_tdgl"] / (q["It_l_tdgl"] + q["It_r_tdgl"])
        f_bias = abs(q["th_r_bias_with"] - q["th_r_bias_without"]) / (2 * np.pi)
        w_ch = min(h.cx_nm - h.a_nm + W / 2, W / 2 - h.cx_nm - h.a_nm) if g.key.startswith("Hole_") \
            else h.cx_nm - h.a_nm + W / 2                                  # the zipper: its 56 nm left channel
        j_circ = 0.5 * (abs(q["Ic_l"]) + abs(q["Ic_r"])) / w_ch / sim.FILM.j_dep_uA_per_nm
        rows.append(dict(name=name, k=k, f_r=q["th_r"] / (2 * np.pi), f_l=-q["th_l"] / (2 * np.pi), f_cur=f_cur,
                         f_tdgl=f_tdgl, xw=(h.cx_nm + W / 2) / W, f_bias=f_bias, q=q, j_circ=j_circ, w_ch=w_ch))
    ok_wind = all(r["q"]["n_relaxed"][r["k"]] == 1 and abs(r["f_r"] + r["f_l"] - 1) < 1e-6
                  and r["q"]["floor"] >= PSI_LOOP_MIN for r in rows)
    checks.append(C.check("2c.1", "placed quantum after 50 tau0 at zero bias: still in its hole, and the two probe "
                          "paths differ by exactly 2 pi", "yes" if ok_wind else "no", "fluxoid kept; f_right + f_left = 1",
                          ok_wind, f"still held after 100 tau0 at {I_SCAN:g} I_dep: "
                          + ", ".join(f"{r['name']}: {'yes' if r['q']['n_bias_with'] == -1 else 'no'}" for r in rows)))
    lin = [r for r in rows if r["j_circ"] < 0.5]
    dev_max = max(abs(r["f_r"] - r["f_cur"]) for r in lin)
    dev_tdgl = max(abs(r["f_r"] - r["f_tdgl"]) for r in lin)
    checks.append(C.check("2c.2", "zero bias, entry from the left: phase per quantum Delta theta(right path) / 2 pi "
                          "against the transport share I_L / I passing left of the hole (reciprocity), holes whose "
                          "circulating current stays below 0.5 j_dep in the narrower channel (London regime)",
                          f"max |difference| {dev_max:.4f} over {len(lin)} holes", "< 0.005", dev_max < 0.005,
                          "; ".join(f"{r['name']}: {r['f_r']:.4f} vs {r['f_cur']:.4f} (circulation {r['j_circ']:.2f} j_dep "
                                    f"in {r['w_ch']:.0f} nm)" for r in rows)
                          + f".  With pyTDGL's interpolated current density: max |difference| {dev_tdgl:.4f}"))
    circ = [r for r in rows if r["name"].startswith("circle")]
    dev_xw = max(abs(r["f_r"] - r["xw"]) for r in circ)
    checks.append(C.check("2c.3", "20 nm circles: phase per quantum against the point-vortex value x/w",
                          f"max |difference| {dev_xw:.3f}", "< 0.03 (small hole ~ point vortex)", dev_xw < 0.03,
                          ", ".join(f"x/w {r['xw']:.2f}: {r['f_r']:.3f}" for r in circ)))
    by = {r["name"]: r for r in rows}
    sc = scan_slips(I_SCAN)
    base = out.get("baseline")
    if ZIP in by and sc is not None and base is not None:
        z, alone = by[ZIP], by["zipper ellipse 16 x 110 nm alone at x = -36"]
        corr = sc[0] - base["slips"]
        checks.append(C.check("2c.4", f"zipper hole 1: phase per quantum at {I_SCAN:g} I_dep vs the flux-scan runs at "
                              f"{I_SCAN:g} that end in [0 -1 0 0 0 0 0], less the dark-run baseline",
                              f"{z['f_bias']:.3f} vs {sc[0]:.3f} - {base['slips']:.3f} = {corr:.3f} ({sc[3]} runs)",
                              "equal within 0.01", abs(z["f_bias"] - corr) < 0.01,
                              f"the plan's 0.32 vs ~0.4: ellipse alone at zero bias {alone['f_r']:.3f} (x/w 0.32); other "
                              f"zipper holes {z['f_r'] - alone['f_r']:+.3f}; bias {I_SCAN:g} I_dep {z['f_bias'] - z['f_r']:+.3f}; "
                              f"slip baseline {base['slips']:+.3f}: a dark run slips that much after t0 = 20 tau0 "
                              f"({base['slips_by_100']:.3f} of it by 100 tau0), the bias transient still settling"))

    t = out.get("trapped")
    if t is not None:
        I_c = 0.5 * (abs(t["Ic_l"]) + abs(t["Ic_r"]))
        net = t["Ic_l"] + t["Ic_r"]
        checks.append(C.check("2c.5", f"one quantum in a centred {2 * R_SW:g} nm hole, zero bias, 400 tau0",
                              f"fluxoid {t['hold_n']} on {t['hold_frames']} frames ({t['hold_rejected']} rejected); "
                              f"circulating current {I_c:.2f} uA per side, net {net:+.3f} uA", "stays 1; net ~ 0",
                              t["hold_n"] == [1] and abs(net) < 0.02 * I_c))
        if isw0 is not None:
            dI = isw0["I_sw"] - t["I_sw"]
            model, D = C.two_channel_isw(isw0["I_sw"] / 2, isw0["I_sw"] / 2, I_c)
            dI_m = isw0["I_sw"] - model
            checks.append(C.check("2c.6", f"centred {2 * R_SW:g} nm hole: I_sw drop with one quantum, against the GL "
                                  "two-channel estimate", f"{t['I_sw']:.2f} vs {isw0['I_sw']:.2f} uA: {dI:.2f} uA "
                                  f"(estimate {dI_m:.2f} uA)", "estimate within 25 %", abs(dI / dI_m - 1) <= 0.25,
                                  f"the plan's linear estimate 2 I_c = {2 * I_c:.1f} uA is {2 * I_c / dI:.1f}x too large; "
                                  f"phase offset between the channels Delta = {D:.3f} / xi.  Just below I_sw "
                                  f"({t['I_sw'] - sim.scales(geo.PLAIN_WIRE)[0]:.2f} uA): {t['below_outcome']}, "
                                  f"{t['below_slips']:.2f} slips, fluxoid history {t['below_history']}"))
    z = out.get(ZIP)
    st = next((v for v in store.load_json(TRAPPED, {}).values() if isinstance(v, dict)
               and v.get("regenerated") == "0 -1 0 0 0 0 0" and "i_sw" in v), None)
    if z is not None and st is not None:
        i_dep = sim.FILM.I_dep_uA(geo.LONG_ZIPPER.free_width_min_nm())
        isw_empty = store.load_json("results/hitmap/isw_conserving_T400_m1.4.json", {})["Long_Zipper"]["I_sw"]
        I_c = 0.5 * (abs(z["Ic_l"]) + abs(z["Ic_r"]))
        f_l = z["It_l"] / (z["It_l"] + z["It_r"])
        n_st = int(st["regenerated"].split()[1])                 # the stored state's quantum (placed here: +1)
        adds_left = np.sign(n_st * z["Ic_l"]) == np.sign(z["It_l"])
        model, _ = C.two_channel_isw(f_l * isw_empty, (1 - f_l) * isw_empty, I_c, adds_left)
        drop_m, drop = (isw_empty - model) / i_dep, isw_empty / i_dep - st["i_sw"]
        checks.append(C.check("2c.7", "zipper [0 -1 0 0 0 0 0]: GL two-channel estimate of the I_sw drop vs the stored "
                              "trapped-state result", f"estimate {drop_m:.3f} I_dep, stored {drop:.3f} I_dep",
                              "same order (two channels standing in for a seven-hole device)",
                              note=f"I_c = {I_c:.2f} uA at zero bias, circulation along the transport on the "
                              f"{'left' if adds_left else 'right'}, transport share left of the hole {f_l:.3f}; "
                              f"linear estimate I_c / share = {I_c / (f_l if adds_left else 1 - f_l) / i_dep:.3f} I_dep"))
    return checks, extra
