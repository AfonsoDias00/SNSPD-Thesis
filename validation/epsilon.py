"""
1b  How pyTDGL 0.9 treats epsilon < -1, and which epsilon values the production hotspot reaches.

With n(0) = 1 at the arrival site, the hotspot centre starts at epsilon = 1 + A f = -1.5 for f = 1.  The tdgl.solve
docstring gives epsilon in [-1, 1].  The solver code (tdgl/solver/solver.py) only rejects epsilon > 1, once, at
t = 0; the update uses (epsilon - |psi|^2) psi as given, with no clipping.  So values below -1 are used as given, and
the only clip is the hotspot's own guard eps_floor (snspd/hotspot.py).

  accepted   the solver takes epsilon = -3 and rejects 1.2
  rate       a uniform film held at constant epsilon < 0, started at psi = 1, relaxes as |psi| ~ exp(epsilon t / u)
             once |psi| is small (linearised GTDGL).  pyTDGL's update is explicit in this term, so at the late-time
             step dt the discrete rate is exactly -ln(1 + epsilon dt / u) / dt (1.7-5.6 % above |epsilon| / u at
             dt = 0.2).  Used as given: the rates for -1, -1.5 and -3 are in the ratio 1 : 1.52 : 3.12 at dt = 0.2.
             Clipped to -1, they would all be equal.  This also checks the time unit through u.
  range      epsilon at the hotspot centre at t0, for every production geometry and full-grid hit site: the
             conserving model renormalises the truncated Gaussian near walls, so edge and bridge sites start deeper
             than -1.5.  It is the minimum over the run (n and exp(-t/t_relax) only decay), checked on one run.
"""
from __future__ import annotations

import numpy as np

from snspd import geometry as geo
from snspd import hotspot as hsm
from snspd import sim
from snspd.hotspot import HOTSPOT
from snspd.units import FILM

from . import common as C

ID, TITLE = "1b", "epsilon mapping"
EPS_RATE = (-1.0, -1.5, -3.0)


def decay_rate(dev, eps):
    """(late-time decay rate of the mean |psi| of a uniform film held at epsilon, fitted where 1e-5 < |psi| < 1e-2;
    the solver step there)."""
    with C.solved_free(dev, 100.0, eps=eps, save_every=5) as sol:
        fr = C.psi_frames(sol)
        t_step, dt = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.dt)
    t = np.array([tm for tm, _ in fr])
    a = np.array([np.abs(p).mean() for _, p in fr])
    s = (a > 1e-5) & (a < 1e-2)
    late = (t_step >= t[s].min()) & (t_step <= t[s].max())
    return -np.polyfit(t[s], np.log(a[s]), 1)[0], float(np.median(dt[late]))


def run(a):
    import tdgl
    checks = []
    dev = C.square(40.0, 1.4)

    # 1. what the solver accepts
    opts = tdgl.SolverOptions(solve_time=1.0, output_file=None)
    acc = {}
    for e in (-3.0, 1.2):
        try:
            tdgl.TDGLSolver(dev, opts, disorder_epsilon=e)
            acc[e] = "accepted"
        except ValueError:
            acc[e] = "rejected"
    checks.append(C.check("1b.1", "pyTDGL 0.9 TDGLSolver given epsilon = -3 / 1.2", f"{acc[-3.0]} / {acc[1.2]}",
                          "accepted / rejected (only epsilon > 1 is checked, at t = 0)",
                          acc == {-3.0: "accepted", 1.2: "rejected"}))

    # 2. relaxation rate at constant epsilon < 0
    rates, steps = {}, {}
    for e in EPS_RATE:
        rates[e], steps[e] = decay_rate(dev, e)
    disc = {e: -np.log(1.0 + e * steps[e] / FILM.u) / steps[e] for e in EPS_RATE}
    rel = {e: rates[e] / disc[e] for e in EPS_RATE}
    checks.append(C.check("1b.2", "uniform film at epsilon = -1 / -1.5 / -3: decay rate of |psi| over the discrete "
                          "prediction -ln(1 + epsilon dt / u) / dt", " / ".join(f"{rel[e]:.4f}" for e in EPS_RATE),
                          "1 within 1 % each (u = 5.79)", all(abs(v - 1) < 0.01 for v in rel.values()),
                          "rates [1/tau0]: " + ", ".join(f"{e:g}: {rates[e]:.4f} (|eps|/u = {abs(e) / FILM.u:.4f}, "
                                                         f"dt = {steps[e]:.3f})" for e in EPS_RATE)))
    ratio = rates[-1.5] / rates[-1.0], rates[-3.0] / rates[-1.0]
    want = disc[-1.5] / disc[-1.0], disc[-3.0] / disc[-1.0]
    checks.append(C.check("1b.3", "rate ratio epsilon = -1.5 : -1 and -3 : -1", f"{ratio[0]:.3f}, {ratio[1]:.3f}",
                          f"{want[0]:.3f}, {want[1]:.3f} if used as given (1, 1 if clipped to -1)",
                          abs(ratio[0] / want[0] - 1) < 0.02 and abs(ratio[1] / want[1] - 1) < 0.02))

    # 3. epsilon at the hotspot centre at t0 for every production site
    geoms = list(geo.GEOMETRIES.values()) + [geo.plain_strip(w) for w in (30.0, 56.0, 100.0)]
    lines, eps_min, worst = [], np.inf, None
    for g in geoms:
        hs = sim.context(g, 1.4)["hs"]
        e0 = {(x, y): 1.0 + HOTSPOT.amp * float(hs.initial_condition(x * 1e-3, y * 1e-3).max())
              for x, y in g.hit_grid("full")}
        (lo_site, lo), (hi_site, hi) = min(e0.items(), key=lambda kv: kv[1]), max(e0.items(), key=lambda kv: kv[1])
        lines.append(f"{g.label} {hi:.2f} at {hi_site} .. {lo:.2f} at {lo_site}")
        if lo < eps_min:
            eps_min, worst = lo, (g, lo_site)
    checks.append(C.check("1b.4", "epsilon at the hotspot centre at t0, f = 1, over every production hit site",
                          f"{eps_min:.2f} at {worst[0].label} {worst[1]} (deepest)",
                          f"> eps_floor = {HOTSPOT.eps_floor:g}: the floor never acts, epsilon is linear in f everywhere",
                          eps_min > HOTSPOT.eps_floor, "per geometry, shallowest .. deepest: " + "; ".join(lines)))

    # 4. the minimum over a run is the value at t0 (deepest site, the full window)
    g, site = worst
    hs = sim.context(g, 1.4)["hs"]
    eps = hsm.build_eps(hs, "conserving", *site, 1.0, C.CFG.tmax)
    ts = HOTSPOT.t0_tau0 + np.arange(0.0, C.CFG.tmax - HOTSPOT.t0_tau0, 0.5)
    mins = np.array([eps(hs.P, t=t).min() for t in ts])
    checks.append(C.check("1b.5", f"min epsilon over 20-400 tau0 at {g.label} {site}",
                          f"{mins.min():.3f} at t = {ts[mins.argmin()]:g} tau0", f"at t0 = {HOTSPOT.t0_tau0:g}",
                          ts[mins.argmin()] == HOTSPOT.t0_tau0))
    f_normal = 1.0 / abs(HOTSPOT.amp)
    checks.append(C.check("1b.6", "photon energy factor below which the centre of an interior hotspot never goes "
                          "normal (epsilon > 0)", f"f < {f_normal:.2f}", "1 / |A|",
                          note="epsilon < -1 is outside the documented range but is physical in a two-temperature "
                          "picture (T_e > 2 T_c - T at the core); it sets how fast psi collapses (rate |epsilon| / u) and "
                          "the proximity length xi / sqrt|epsilon| at the hotspot edge, not the size of the normal core"))
    return checks, dict(rates=rates, eps_by_geometry=lines)
