"""
3c  The readout circuit (snspd/circuit.py) and the unit identities it rests on.

  LR          with no voltage on the strip the detector current returns to the bias as exp(-t R_L / L_k).
  squares     a strip held normal (epsilon = -1) shows mu-drop / (I / I_u) = its resistance in squares between the
              probes: V_u / I_u = R_sq.  Plain strip: exact geometry.  Double column: against the series sum of
              cross-sections, which ignores the current spreading at the slot ends (a lower bound).
  divider     the normal strip in the circuit: the current settles at I_b k / (k + n), k = R_L / R_sq, n squares,
              with the time constant N_sq / (k + n) tau0.
"""
from __future__ import annotations

import numpy as np
import tdgl

from snspd import geometry as geo
from snspd import sim
from snspd.circuit import Circuit, CircuitState
from snspd.units import FILM

from . import common as C

ID, TITLE = "3c", "readout circuit: LR relaxation, resistance in squares, current divider"
MESH = 2.0


def normal_strip(geom, circuit=None, tmax=60.0, i_b=0.5):
    """A strip with epsilon = -1 (normal), biased through the circuit or a current source.  Returns
    (v_end / (I_end / I_u), I_end / I_b, t at which I has made 1 - 1/e of its way)."""
    cfg = sim.RunConfig(mesh_xi=MESH)
    ctx = sim.context(geom, cfg)
    I_b = i_b * ctx["I_dep"]
    circ = CircuitState(circuit, FILM, I_b) if circuit else None
    cur = circ.currents if circ else dict(source=I_b, drain=-I_b)
    solver = tdgl.TDGLSolver(ctx["device"], sim.solver_options(cfg), terminal_currents=cur, disorder_epsilon=-1.0)
    p0, p1 = ctx["probes"]
    T, V, I = [], [], []

    def hook(t, dt, res):
        v = float(res.mu[p0] - res.mu[p1])
        if circ:
            circ.step(dt, v)
        T.append(t), V.append(v), I.append(circ.I if circ else I_b)

    sim.march(solver, solver.psi_init * 1e-3, solver.mu_init, tmax, hook)
    T, I = np.array(T), np.array(I)
    t63 = float(T[np.argmax(I <= I[-1] + (I_b - I[-1]) / np.e)]) if circ else None
    return V[-1] / (I[-1] / FILM.I_unit_uA), I[-1] / I_b, t63


def run(a):
    checks = []
    # LR relaxation
    c = Circuit(n_sq=200.0)
    s = CircuitState(c, FILM, 50.0)
    s.I = 10.0
    tau, t, dt = c.tau_reset(FILM), 0.0, 0.025
    while t < 2 * tau:
        s.step(dt, 0.0)
        t += dt
    exact = 50.0 - 40.0 * np.exp(-t / tau)
    checks.append(C.check("3c.1", f"LR recovery after 2 L_k/R_L ({2 * tau:.0f} tau0): I [uA]", f"{s.I:.4f}",
                          f"{exact:.4f} (to 1e-3)", abs(s.I - exact) < 1e-3))

    # resistance in squares
    r, _, _ = normal_strip(geo.PLAIN_WIRE)
    sq = geo.PLAIN_WIRE.squares()
    checks.append(C.check("3c.2", "normal plain strip: mu-drop / (I / I_u) between the probes [squares]", f"{r:.3f}",
                          f"{sq:.3f} (within 1 %)", abs(r / sq - 1) < 0.01))
    r2, _, _ = normal_strip(geo.DOUBLE_400)
    sq2 = geo.DOUBLE_400.squares()
    checks.append(C.check("3c.3", "normal double column: the same [squares]", f"{r2:.2f}",
                          f">= {sq2:.2f} (series sum of cross-sections), by the spreading resistance at the slot ends",
                          sq2 <= r2 <= 1.15 * sq2))

    # current divider
    c = Circuit(n_sq=40.0)
    k = c.RL_ohm / FILM.R_sq_ohm
    r3, frac, t63 = normal_strip(geo.PLAIN_WIRE, c, tmax=120.0)
    checks.append(C.check("3c.4", "normal plain strip in the circuit (40 squares in series): I_end / I_b",
                          f"{frac:.4f}", f"k / (k + n) = {k / (k + r3):.4f} (within 1 %)", abs(frac / (k / (k + r3)) - 1) < 0.01))
    checks.append(C.check("3c.5", "... time for 63 % of the current change [tau0]", f"{t63:.1f}",
                          f"N_sq / (k + n) = {c.n_sq / (k + r3):.1f} (within 10 %)",
                          abs(t63 / (c.n_sq / (k + r3)) - 1) < 0.10))
    return checks, {}
