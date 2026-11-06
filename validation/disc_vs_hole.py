"""
1c  A frozen hotspot against a real hole: a permanent disc of epsilon = -1 (and -3), radius R, at the centre of the
plain 200 nm wire (no diffusion, no decay, kept in the sustained-switching continuation), against a vacuum hole of
the same radius.  Dark switching current, production window and criterion, R = 20 and 40 nm.

  agreement   the disc and the hole should differ only through the proximity effect over ~xi: the disc depletes
              psi around itself (a flat N/S interface with epsilon_N = -1 costs ~0.5 xi of superfluid width per
              side, -3 ~0.9 xi; a hole's insulating wall costs nothing).  Measured as the extra radius the hole
              would need to switch at the disc's current, assuming I_sw follows the free width.
  free width  the plan expects i_sw = I_sw / I_sw(plain) ~ (w - 2R) / w "with a small crowding penalty".  A circular
              hole crowds the current at its sides (London: up to 2x for a small hole in a wide sheet), unlike the
              zipper's long ellipses, so (w - 2R) / w is an upper bound; how far below it the holes switch is
              measured.
"""
from __future__ import annotations

import time

import numpy as np

from snspd import geometry as geo
from snspd import sim, store

from . import common as C

ID, TITLE = "1c", "frozen normal disc vs a real hole"
CASES = (("plain", 0.0, None), ("hole", 20.0, None), ("disc", 20.0, -1.0), ("disc", 20.0, -3.0), ("hole", 40.0, None),
         ("disc", 40.0, -1.0))
BRACKET = (0.4, 1.2)          # [I_dep of the free width]; wider below than sim.ISW_BRACKET in case of strong crowding
W = geo.PLAIN_WIRE.width_nm


def cache_key(kind, R, eps):
    return f"isw_{kind}_R{R:g}" + ("" if eps is None else f"_eps{eps:g}") + f"_{C.CFG.tag}"


def eps_disc(R_nm, value):
    R2 = (R_nm * 1e-3) ** 2

    def eps(r, *, vectorized=True):
        r = np.asarray(r)
        return np.where(r[:, 0] ** 2 + r[:, 1] ** 2 < R2, value, 1.0)

    return eps


def isw_task(task):
    """Dark switching current of a centred hole (production run_simulation) or frozen disc (eps kept throughout)."""
    t_start = time.time()
    try:
        kind, R, e = task["kind"], task["R"], task["eps"]
        cfg = C.CFG
        if kind in ("hole", "plain"):
            g = C.holed_strip(0.0, R) if kind == "hole" else geo.PLAIN_WIRE
            dev = sim.context(g, cfg.mesh_xi)["device"]
            run = lambda I: sim.run_simulation(dev, I, None, cfg)
        else:
            dev = sim.context(geo.PLAIN_WIRE, cfg.mesh_xi)["device"]
            eps = eps_disc(R, e)
            run = lambda I: C.run_frozen(dev, I, eps, cfg)
        i_free = sim.FILM.I_dep_uA(W - 2 * R)
        I, flag, n, trig = C.bisect(run, BRACKET[0] * i_free, BRACKET[1] * i_free, sim.scales(geo.PLAIN_WIRE)[0])
        return dict(ok=True, task=task, res=dict(I_sw=I, flag=flag, n=n, branch=trig.branch, t=trig.t),
                    wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def isw_values(cases, workers):
    """{(kind, R, eps): result}, computed where missing and kept in the cache."""
    cache = C.load_cache()
    todo = [dict(kind=k, R=R, eps=e) for k, R, e in cases if cache_key(k, R, e) not in cache]

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t}: {r['error']}")
            return
        cache[cache_key(t["kind"], t["R"], t["eps"])] = dict(r["res"], wall_s=round(r["wall_s"], 1))
        C.save_cache(cache)
        sim.log(f"   {t['kind']} R = {t['R']:g} nm eps = {t['eps']}: I_sw = {r['res']['I_sw']:.2f} uA "
                f"{r['res']['flag']}  ({r['wall_s']:.0f} s)")

    sim.log(f"-- {len(todo)} switching-current searches to do")
    sim.execute(todo, workers, on_result, fn=isw_task)
    return {(k, R, e): cache[cache_key(k, R, e)] for k, R, e in cases if cache_key(k, R, e) in cache}


def interface_loss(eps_n):
    """Superfluid width lost at a flat N/S interface (N: epsilon_n for x < 0) against an insulating wall, in xi:
    int_0^inf (1 - psi^2) dx - int_-inf^0 psi^2 dx, from the 1-D GL equation psi'' + epsilon psi - psi^3 = 0."""
    from scipy.integrate import solve_bvp
    L = 12.0
    x = np.linspace(-L, L, 801)
    e = lambda x: np.where(x < 0, eps_n, 1.0)
    sol = solve_bvp(lambda x, y: np.vstack([y[1], -(e(x) - y[0] ** 2) * y[0]]),
                    lambda ya, yb: np.array([ya[0], yb[0] - 1.0]), x,
                    np.vstack([0.5 * (1 + np.tanh(x)), 0.5 / np.cosh(x) ** 2]), tol=1e-8, max_nodes=100000)
    xs = np.linspace(-L, L, 4001)
    p2 = sol.sol(xs)[0] ** 2
    return float(np.trapezoid(np.where(xs > 0, 1 - p2, -p2), xs))


def run(a):
    vals = isw_values(CASES, a.workers)
    isw_plain = vals[("plain", 0.0, None)]["I_sw"]
    checks, extra = [], dict(isw={f"{k}_{R:g}_{e}": v for (k, R, e), v in vals.items()}, I_sw_plain=isw_plain)
    for R in sorted({R for _, R, _ in CASES}):
        hole = vals.get(("hole", R, None))
        if hole is None:
            continue
        for k, R_, e in CASES:
            if k != "disc" or R_ != R or (k, R_, e) not in vals:
                continue
            disc = vals[(k, R_, e)]
            dR = (1 - disc["I_sw"] / hole["I_sw"]) * (W - 2 * R) / 2
            delta = interface_loss(e) * C.XI
            checks.append(C.check(
                f"1c.{len(checks) + 1}", f"R = {R:g} nm: disc (epsilon = {e:g}) vs hole, dark I_sw",
                f"{disc['I_sw']:.2f} vs {hole['I_sw']:.2f} uA: the disc acts as a hole {dR:+.1f} nm larger in radius",
                "0 to 2 xi = 10 nm larger (proximity effect only)", 0 <= dR <= 2 * C.XI,
                f"GL flat-interface estimate: {delta:.1f} nm; flags: disc '{disc['flag']}', hole '{hole['flag']}'"))
        ratio = hole["I_sw"] / isw_plain
        free = (W - 2 * R) / W
        checks.append(C.check(f"1c.{len(checks) + 1}", f"R = {R:g} nm hole: i_sw = I_sw / I_sw(plain) against the "
                              "free-width fraction (w - 2R)/w", f"{ratio:.3f} vs {free:.3f} ({ratio / free:.3f} of it)",
                              "within 10 % of (w - 2R)/w (plan: free width, small crowding penalty)",
                              abs(ratio / free - 1) <= 0.10,
                              f"I_sw / I_dep(free width) = {hole['I_sw'] / sim.FILM.I_dep_uA(W - 2 * R):.3f} vs "
                              f"{isw_plain / sim.FILM.I_dep_uA(W):.3f} for the plain wire"))
    return checks, extra
