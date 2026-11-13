"""
3b  Energy balance of the generalised TDGL: is the heat that snspd/thermal.py takes from the solver all of it?

For pyTDGL's equations with a fixed epsilon, the free energy F = sum_sites A (-eps n + n^2/2) + sum_edges w |d psi|^2
(n = |psi|^2, units 2 E_c xi^2 d) obeys

    dF/dt = 2 I V - 2 sum_sites A ( |grad mu|^2 + W_relax ),
    W_relax = u / sqrt(1 + gamma^2 n) ( |(d_t + i mu) psi|^2 + gamma^2/4 (d_t n)^2 ),

with I the bias current [I_u] and V the potential difference between the terminals [V_u].  So the work of the source
that does not stay in F is dissipated as Joule heat of the normal current plus order-parameter relaxation, and the
heat model's source term eta0 (|grad mu|^2 + W_relax) is complete.  Checked on a plain strip that is driven hard:
the bias is switched on at t = 0 above the dark switching current, so the strip goes through the bias transient,
phase slips and a growing resistive state.  The same run without the relaxation term shows how much heat a
Joule-only model (Zotova & Vodolazov 2012) leaves out.
"""
from __future__ import annotations

import numpy as np
import tdgl

from snspd import geometry as geo
from snspd import sim, thermal
from snspd.units import FILM

from . import common as C

ID, TITLE = "3b", "energy balance of the TDGL: source work = change of free energy + Joule + relaxation"
MESH = 2.0


def run(a):
    checks = []
    cfg = sim.RunConfig(mesh_xi=MESH)
    ctx = sim.context(geo.PLAIN_WIRE, cfg)
    dev = ctx["device"]
    h = thermal.HeatField(dev, cfg.thermal, FILM)          # used for its mesh operators only
    A, w, e0, e1 = h.area, h.w, h.e0, h.e1
    info = {t.name: t for t in dev.terminal_info()}
    src, drn = info["source"].site_indices, info["drain"].site_indices

    def free_energy(psi):
        n = np.abs(psi) ** 2
        return float(A @ (-n + 0.5 * n ** 2) + w @ np.abs(psi[e1] - psi[e0]) ** 2)

    for cid, i_b, tmax, label in (("3b.1", 0.80, 60.0, "bias 0.80 I_dep switched on at t = 0 (transient, then steady)"),
                                  ("3b.2", 1.15, 60.0, "bias 1.15 I_dep switched on at t = 0 (phase slips, resistive)")):
        I = i_b * ctx["I_dep"]
        i_dimless = I / FILM.I_unit_uA
        solver = tdgl.TDGLSolver(dev, sim.solver_options(cfg), terminal_currents=dict(source=I, drain=-I))
        acc = dict(work=0.0, joule=0.0, relax=0.0)
        prev = dict(psi=np.array(solver.psi_init), mu=np.array(solver.mu_init))
        F0 = free_energy(prev["psi"])

        def hook(t, dt, res):
            psi, mu = np.asarray(res.psi), np.asarray(res.mu).real
            acc["work"] += 2 * i_dimless * (mu[src].mean() - mu[drn].mean()) * dt
            acc["joule"] += 2 * float(A @ h.joule(mu)) * dt
            acc["relax"] += 2 * float(A @ h.relaxation(dt, prev["psi"], prev["mu"], psi)) * dt
            prev.update(psi=psi.copy(), mu=mu.copy())

        psi, mu, t = sim.march(solver, prev["psi"], prev["mu"], tmax, hook)
        dF = free_energy(np.asarray(psi)) - F0
        res = acc["work"] - dF - acc["joule"] - acc["relax"]
        scale = acc["joule"] + acc["relax"]
        checks.append(C.check(cid, f"{label}: residual of work - dF - Joule - relaxation over {tmax:g} tau0 "
                                   f"[fraction of the dissipated energy]",
                              f"{res / scale:+.4f}  (work {acc['work']:.1f}, dF {dF:+.1f}, Joule {acc['joule']:.1f}, "
                              f"relaxation {acc['relax']:.1f})", "|residual| <= 0.02", abs(res / scale) <= 0.02))
        checks.append(C.check(cid + "b", "... share of the dissipation that is order-parameter relaxation",
                              f"{acc['relax'] / scale:.2f}", "what a Joule-only heat source leaves out"))
    return checks, {}
