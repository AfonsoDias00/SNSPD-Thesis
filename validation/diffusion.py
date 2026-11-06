"""
1a  The conserving hotspot's diffusion solver alone (snspd/hotspot.HotspotSolver, no TDGL): dn/dt = D lap n with
insulating walls, P1 elements on the device mesh, lumped mass, Crank-Nicolson.  Exact solutions, with
s^2 = r0^2 + 4 D t:

  conservation  integral n dA = pi r0^2 at every step, plain wire, 400 nm and 110 nm slot columns, 380 tau0
  free space    hotspot at the centre of a 600 nm square (300 nm from every wall): n = (r0^2/s^2) exp(-r^2/s^2),
                <r^2> = r0^2 + 4 D t
  wall          source a nm from one straight insulating wall.  The code truncates the Gaussian at the wall and
                renormalises it to pi r0^2, so the exact solution is the free evolution of that truncated profile
                plus its mirror image: per axis (r0/s) exp(-(x-a)^2/s^2) erfc(-(4Dt a + r0^2 x)/(r0 s sqrt(4Dt)))/2,
                plus the same at -x, times 2/erfc(-a/r0).  (A bare image source is exact only for a >> r0.)
  strip modes   integrated over a whole strip (all walls insulating) the problem is 1-D in x, and the first mode
                across the width, a_1 = int n cos(pi (x + w/2)/w) dA, decays exactly as exp(-D (pi/w)^2 t)
  1-D regime    once the width modes have decayed, the peak is sqrt(pi) r0^2 / (w s), i.e. ~ t^(-1/2)
  convergence   Crank-Nicolson: error ~ dt^2 at fixed mesh; P1: error ~ h^2 at small dt
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.special import erfc

from snspd import geometry as geo
from snspd import sim
from snspd.hotspot import HOTSPOT, HotspotSolver

from . import common as C

ID, TITLE = "1a", "hotspot diffusion solver"
R0, D = HOTSPOT.r0_nm, HOTSPOT.D_nm2_per_tau0
MESH = 1.4


def free(P, t, x0=0.0, y0=0.0):
    s2 = R0 ** 2 + 4 * D * t
    return (R0 ** 2 / s2) * np.exp(-((P[:, 0] - x0) ** 2 + (P[:, 1] - y0) ** 2) / s2)


def wall(P, t, a, x_wall):
    """Exact Neumann solution for the code's initial condition: source at x_wall + a, wall at x = x_wall, y0 = 0."""
    s, L = np.sqrt(R0 ** 2 + 4 * D * t), np.sqrt(4 * D * t)
    U = lambda x: (R0 / s) * np.exp(-(x - a) ** 2 / s ** 2) * 0.5 * erfc(-(L ** 2 * a + R0 ** 2 * x) / (R0 * s * L))
    X = P[:, 0] - x_wall
    return 2.0 / erfc(-a / R0) * (U(X) + U(-X)) * (R0 / s) * np.exp(-P[:, 1] ** 2 / s ** 2)


def rel_err(m, n, ref):
    return float(np.sqrt(m @ (n - ref) ** 2 / (m @ ref ** 2)))


def solver(dev, dt=None):
    return HotspotSolver(dev.points, dev.triangles, HOTSPOT if dt is None else replace(HOTSPOT, dt_tau0=dt))


def evolve64(hs, n, t_end):
    """The production Crank-Nicolson step in float64 (HotspotSolver.evolve stores float32 snapshots)."""
    for _ in range(int(round(t_end / hs.dt))):
        n = hs._lu.solve(hs._B @ n)
    return n


def mean_edge_nm(dev):
    T, P = np.asarray(dev.triangles), np.asarray(dev.points) * 1e3
    e = np.concatenate([np.hypot(*(P[T[:, i]] - P[T[:, (i + 1) % 3]]).T) for i in range(3)])
    return float(e.mean())


def run(a):
    checks, extra = [], {}
    window = C.CFG.tmax - HOTSPOT.t0_tau0

    # conservation on the production meshes
    for k, (g, sites) in enumerate(((geo.PLAIN_WIRE, ((-96.0, 0.0), (0.0, 0.0))),
                                    (C.SLOT_400, ((0.0, 220.0), (-28.0, 200.0), (-96.0, 0.0))),
                                    (C.SLOT_110, ((0.0, 75.0), (-96.0, 0.0))))):
        hs = sim.context(g, MESH)["hs"]
        drift, neg = 0.0, 0.0
        for x, y in sites:
            snaps = hs.evolve(hs.initial_condition(x * 1e-3, y * 1e-3), window)
            mass = snaps.astype(float) @ hs.m / (np.pi * (R0 * 1e-3) ** 2)
            drift = max(drift, float(np.abs(mass - 1).max()))
            neg = min(neg, float(snaps.min() / snaps[0].max()))
        checks.append(C.check(f"1a.{k + 1}", f"{g.label}: integral n dA / pi r0^2 over {window:g} tau0, sites {sites}",
                              f"max |deviation| {drift:.1e}", "< 1e-5 (float32 snapshots)", drift < 1e-5,
                              f"most negative n / n(0) peak: {neg:.1e}"))

    # free space: centre of a 600 nm square
    dev = C.square(600.0, MESH)
    hs = solver(dev)
    P = hs.P * 1e3
    times = (10.0, 50.0, 100.0, 200.0)
    snaps = hs.evolve(hs.initial_condition(0.0, 0.0), max(times)).astype(float)
    e0 = rel_err(hs.m, snaps[0], np.exp(-(P[:, 0] ** 2 + P[:, 1] ** 2) / R0 ** 2))
    errs = {t: rel_err(hs.m, snaps[int(round(t / hs.dt))], free(P, t)) for t in times}
    peaks = {t: snaps[int(round(t / hs.dt))].max() / (R0 ** 2 / (R0 ** 2 + 4 * D * t)) for t in times}
    checks.append(C.check("1a.4", f"free space (600 nm square, {MESH:g} xi): profile vs (r0^2/s^2) exp(-r^2/s^2)",
                          f"rel. L2 error {max(errs.values()):.4f} (max over t = {', '.join(f'{t:g}' for t in times)})",
                          "< 0.03", max(errs.values()) < 0.03,
                          f"t = 0 (interpolated Gaussian): {e0:.4f}; per t: " + ", ".join(f"{t:g}: {e:.4f}" for t, e in errs.items())
                          + "; peak / exact: " + ", ".join(f"{t:g}: {p:.3f}" for t, p in peaks.items())))
    tt = np.arange(len(snaps)) * hs.dt
    r2 = (snaps * (hs.m * (P[:, 0] ** 2 + P[:, 1] ** 2))).sum(1) / (snaps @ hs.m)
    slope, icpt = np.polyfit(tt, r2, 1)
    checks.append(C.check("1a.5", "free space: d<r^2>/dt over 4D (a wrong D or tau0 unit shows here)",
                          f"{slope / (4 * D):.4f}", "1 within 1 %", abs(slope / (4 * D) - 1) < 0.01,
                          f"intercept {icpt:.1f} nm^2 vs r0^2 = {R0 ** 2:g}"))

    # one straight wall (the square's left side)
    for k, a_nm in enumerate((0.0, 4.0, 3 * R0)):
        snaps = hs.evolve(hs.initial_condition((-300.0 + a_nm) * 1e-3, 0.0), 100.0).astype(float)
        ew = {t: rel_err(hs.m, snaps[int(round(t / hs.dt))], wall(P, t, a_nm, -300.0)) for t in (10.0, 50.0, 100.0)}
        note = ", ".join(f"t = {t:g}: {e:.4f}" for t, e in ew.items())
        if a_nm == 3 * R0:
            nf = rel_err(hs.m, snaps[-1], free(P, 100.0, -300.0 + a_nm) * (P[:, 0] >= -300.0))
            note += f";  without the image term at t = 100: {nf:.3f}"
        checks.append(C.check(f"1a.{6 + k}", f"source {a_nm:g} nm from an insulating wall vs the exact truncated-"
                              "Gaussian + image solution", f"rel. L2 error {max(ew.values()):.4f}", "< 0.03",
                              max(ew.values()) < 0.03, note))

    # strip modes and the 1-D regime
    for k, w in enumerate((200.0, 56.0)):
        g = geo.PLAIN_WIRE if w == 200.0 else geo.plain_strip(w)
        hs = sim.context(g, MESH)["hs"]
        Pw = hs.P * 1e3
        t_fit = 380.0 if w == 200.0 else 150.0
        snaps = hs.evolve(hs.initial_condition((-w / 2 + 4) * 1e-3, 0.0), t_fit).astype(float)
        a1 = snaps @ (hs.m * np.cos(np.pi * (Pw[:, 0] + w / 2) / w))
        tt = np.arange(len(snaps)) * hs.dt
        rate = -np.polyfit(tt, np.log(np.abs(a1)), 1)[0]
        want = D * (np.pi / w) ** 2
        checks.append(C.check(f"1a.{9 + k}", f"{w:g} nm strip, edge source: decay rate of the first width mode over "
                              f"D (pi/w)^2 = {want:.2e} / tau0 (fit 0-{t_fit:g} tau0)", f"{rate / want:.4f}",
                              "1 within 2 %", abs(rate / want - 1) < 0.02, f"1/rate = {1 / rate:.0f} tau0"))
    ratios = {}
    for w in (56.0, 200.0):
        g = geo.PLAIN_WIRE if w == 200.0 else geo.plain_strip(w)
        hs = sim.context(g, MESH)["hs"]
        snaps = hs.evolve(hs.initial_condition(0.0, 0.0), window).astype(float)
        ratios[w] = {t: snaps[int(round(t / hs.dt))].max() / (np.sqrt(np.pi) * R0 ** 2 / (w * np.sqrt(R0 ** 2 + 4 * D * t)))
                     for t in (50.0, 100.0, 200.0, 380.0)}
    late = [v for t, v in ratios[56.0].items() if t >= 100]
    checks.append(C.check("1a.11", "56 nm strip, centre source: peak over the 1-D value sqrt(pi) r0^2 / (w s), "
                          "t >= 100 tau0", ", ".join(f"{v:.4f}" for v in late), "1 within 1 % (peak ~ t^-1/2)",
                          all(abs(v - 1) < 0.01 for v in late),
                          "200 nm wire, same ratio at t = 50/100/200/380 tau0: "
                          + ", ".join(f"{v:.3f}" for v in ratios[200.0].values()) + " (still 2-D early in the window)"))

    # convergence: time step at the production mesh, mesh at a small step (300 nm square, centre source, t = 20)
    t_c = 20.0
    dev = C.square(300.0, MESH)
    ref = evolve64(solver(dev, 1 / 128), solver(dev, 1 / 128).initial_condition(0.0, 0.0), t_c)
    m = solver(dev).m
    dts = (1.0, 0.5, 0.25, 0.125)
    et = [rel_err(m, evolve64(solver(dev, dt), solver(dev, dt).initial_condition(0.0, 0.0), t_c), ref) for dt in dts]
    ot = [np.log2(et[i] / et[i + 1]) for i in range(len(dts) - 1)]
    checks.append(C.check("1a.12", "Crank-Nicolson: error at t = 20 tau0 vs dt = 1/128, order for dt = 1 -> 0.125",
                          ", ".join(f"{o:.2f}" for o in ot), "2 (error falls ~4x per halving), within 1.8-2.2",
                          all(1.8 <= o <= 2.2 for o in ot),
                          ", ".join(f"dt {d:g}: {e:.1e}" for d, e in zip(dts, et)) + f" (production dt = {HOTSPOT.dt_tau0:g})"))
    meshes = (2.8, 1.4, 0.7, 0.35)
    es, hs_nm = [], []
    for mx in meshes:
        dv = C.square(300.0, mx)
        h = solver(dv, 1 / 32)
        n = evolve64(h, h.initial_condition(0.0, 0.0), t_c)
        es.append(rel_err(h.m, n, free(h.P * 1e3, t_c)))
        hs_nm.append(mean_edge_nm(dv))
    os_ = [np.log(es[i] / es[i + 1]) / np.log(hs_nm[i] / hs_nm[i + 1]) for i in range(len(meshes) - 1)]
    checks.append(C.check("1a.13", "P1 mesh convergence: error vs the exact solution at t = 20 tau0, order in the "
                          "mean edge for 2.8 -> 0.35 xi", ", ".join(f"{o:.2f}" for o in os_), "2, within 1.5-2.5",
                          all(1.5 <= o <= 2.5 for o in os_),
                          ", ".join(f"{mx:g} xi ({h_:.1f} nm): {e:.1e}" for mx, h_, e in zip(meshes, hs_nm, es))))
    extra.update(free_errors=errs, strip_1d_ratio={str(k): v for k, v in ratios.items()}, dt_errors=et, mesh_errors=es)
    return checks, extra
