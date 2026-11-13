"""
3a  The heat model alone (snspd/thermal.py), with no TDGL: the order parameter is a fixed array.

  energy      one photon into a film at rest: electrons + phonons + escaped + sunk - latent = h nu at every time.
  relaxation  a uniformly hot film cools as the two coupled ODEs of the two-temperature model say (reference:
              scipy's stiff integrator on the same equations), for a normal film and for |psi|^2 = 0.8.
  diffusion   in the normal state the electron energy U = t_e^2 spreads with the diffusion constant D:
              <y^2> of a small Gaussian grows as 2 D t.
  walls       energy put in at the film edge stays in the film (no flux through edges).
  k_s         the superconducting-state heat conduction factor at its limits.
  functional  what the two-fluid free energy implies: entropy at the bath, specific-heat jump, condensation energy
              (the model's values against BCS: characterises the approximation, no pass mark).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from snspd import geometry as geo
from snspd import sim, thermal
from snspd.units import FILM

from . import common as C

ID, TITLE = "3a", "heat model alone: conservation, two-temperature relaxation, diffusion"
MESH, DT = 2.0, 0.025


def advance(h, psi, t0, t1, dt=DT):
    mu = np.zeros(h.N)
    t = t0
    while t < t1 - 1e-9:
        h.step(t, dt, psi, mu, psi, mu)
        t += dt
    return t


def run(a):
    checks = []
    dev = sim.make_device(geo.PLAIN_WIRE, MESH)
    N = len(dev.points)
    P0 = thermal.NBTIN
    ph = thermal.Photon()

    # ---- energy conservation ----------------------------------------------------------------------------------
    worst = 0.0
    for label, psi in (("psi = 1", np.ones(N, complex)), ("psi = 0", np.zeros(N, complex))):
        h = thermal.HeatField(dev, P0, FILM, ph, (-58.0, 0.0))
        h.reset(psi)
        t = 0.0
        for t1 in (25.0, 40.0, 100.0, 300.0):
            t = advance(h, psi, t, t1)
            e, p = h.energy()
            A = h.audit
            worst = max(worst, abs(e + p + A["escape"] + A["sink"] - A["latent"] - A["photon"]) / h.E_photon)
    checks.append(C.check("3a.1", "energy audit of one photon, worst residual over 4 times x 2 states [fraction of h nu]",
                          f"{worst:.1e}", "<= 1e-9", worst <= 1e-9))
    checks.append(C.check("3a.2", "photon energy deposited by t = 300 tau0 [fraction of h nu]",
                          f"{h.audit['photon'] / h.E_photon:.6f}", "1 (to 1e-6)", abs(h.audit["photon"] / h.E_photon - 1) < 1e-6))

    # ---- uniform two-temperature relaxation ---------------------------------------------------------------------
    from scipy.integrate import solve_ivp
    for cid, n_s in (("3a.3", 0.0), ("3a.4", 0.8)):
        p = replace(P0, sink_nm=0.0)
        h = thermal.HeatField(dev, p, FILM)
        psi = np.full(N, np.sqrt(n_s), complex)
        h.reset(psi)
        h.U[:] = 1.6 ** 2
        g = 1.0 - p.b * n_s

        def rhs(t, y):
            U, P = y
            X = (U ** (p.q / 2) - P ** (p.q / 4)) / (p.q * p.tau_ep)
            return [-2 * X / g, 4 * p.r * X - (P - p.t0 ** 4) / p.tau_esc]

        times = (5.0, 30.0, 200.0, 1000.0)
        ref = solve_ivp(rhs, (0, times[-1]), [1.6 ** 2, p.t0 ** 4], method="Radau", t_eval=times, rtol=1e-10, atol=1e-14)
        t, err, vals = 0.0, 0.0, []
        for k, t1 in enumerate(times):
            t = advance(h, psi, t, t1)
            te, tp = float(np.sqrt(h.U.mean())), float(h.P.mean() ** 0.25)
            err = max(err, abs(te / np.sqrt(ref.y[0, k]) - 1), abs(tp / ref.y[1, k] ** 0.25 - 1))
            vals.append(f"{te:.3f}/{tp:.3f}")
        checks.append(C.check(cid, f"uniform cooling from 1.6 T_c, |psi|^2 = {n_s:g}: T_e/T_ph [T_c] at t = 5, 30, 200, "
                                   f"1000 tau0; worst deviation from the ODE reference",
                              f"{', '.join(vals)}; {err:.1e}", "<= 5e-3 (first-order steps of 0.025 tau0)", err <= 5e-3,
                              note="spread over the mesh: " + f"{np.ptp(h.U):.1e}"))

    # ---- diffusion in the normal state --------------------------------------------------------------------------
    p = replace(P0, sink_nm=0.0, tau_ep=1e12, tau_esc=1e12)
    h = thermal.HeatField(dev, p, FILM)
    psi = np.zeros(N, complex)
    h.reset(psi)
    x, y = h.xy_nm[:, 0], h.xy_nm[:, 1]
    R0 = 15.0
    h.U += 0.5 * np.exp(-(x ** 2 + y ** 2) / R0 ** 2)
    m0 = float(h.area @ (h.U - h.U0))

    def y2():
        w = h.area * (h.U - h.U0)
        return float(w @ y ** 2 / w.sum())

    a0, t = y2(), 0.0
    t = advance(h, psi, t, 20.0)
    D_meas = (y2() - a0) / (2 * t)
    checks.append(C.check("3a.5", "normal state: D from the growth of <y^2> of a 15 nm Gaussian over 20 tau0 [nm^2/tau0]",
                          f"{D_meas:.3f}", f"{p.D_nm2:.3f} (within 2 %)", abs(D_meas / p.D_nm2 - 1) < 0.02))
    checks.append(C.check("3a.6", "... and the energy is conserved [relative change]",
                          f"{abs(float(h.area @ (h.U - h.U0)) / m0 - 1):.1e}", "<= 1e-10",
                          abs(float(h.area @ (h.U - h.U0)) / m0 - 1) <= 1e-10))

    # ---- reflecting walls -----------------------------------------------------------------------------------------
    h = thermal.HeatField(dev, p, FILM)
    h.reset(psi)
    h.U += 0.5 * np.exp(-((x + 100.0) ** 2 + y ** 2) / R0 ** 2)
    m0 = float(h.area @ (h.U - h.U0))
    advance(h, psi, 0.0, 60.0)
    checks.append(C.check("3a.7", "edge deposit after 60 tau0 (it has spread 40 nm): energy left in the film [fraction]",
                          f"{float(h.area @ (h.U - h.U0)) / m0:.12f}", "1 (to 1e-10)",
                          abs(float(h.area @ (h.U - h.U0)) / m0 - 1) <= 1e-10))

    # ---- k_s ------------------------------------------------------------------------------------------------------
    f = thermal.ks_over_kn(np.array([0.0, 1.0, thermal.gap_at_bath(P0) / P0.t0, 40.0]))
    checks.append(C.check("3a.8", "k_s / k_n at Delta / k_B T = 0, 1, the bath value, 40",
                          ", ".join(f"{v:.4f}" for v in f), "1, between, small, 0",
                          abs(f[0] - 1) < 1e-12 and f[3] < 1e-6 and f[0] > f[1] > f[2] > f[3]))

    # ---- the two-fluid functional ---------------------------------------------------------------------------------
    a0 = P0.eta0 / (1 - P0.t0 ** 2) ** 2
    checks.append(C.check("3a.9", "condensation energy 4 E_c(0) / (gamma T_c^2): model; from lambda, xi and c_e",
                          f"{a0:.3f}; {thermal.eta0_from_film() / (1 - P0.t0 ** 2) ** 2:.3f}",
                          "BCS (3/pi^2)(1.764)^2 = 0.946; Gorter-Casimir 1"))
    checks.append(C.check("3a.10", "entropy of the superconducting state at the bath, s_s / (gamma T) = 1 - a0 (1 - t0^2)",
                          f"{1 - a0 * (1 - P0.t0 ** 2):.3f}", "> 0; BCS c_es / (gamma T) at 0.3 T_c is about 0.2",
                          1 - a0 * (1 - P0.t0 ** 2) > 0))
    checks.append(C.check("3a.11", "specific-heat jump at T_c, Delta C / (gamma T_c) = 2 a0", f"{2 * a0:.2f}",
                          "BCS 1.43; Gorter-Casimir 2"))
    return checks, {}
