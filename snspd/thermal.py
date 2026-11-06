"""
Electro-thermal TDGL (prototype, 5 Nov 2026): pyTDGL coupled to a local electron temperature with Joule heating.

Model (Zotova & Vodolazov, PRB 85, 024509 (2012), Eq. (2); arXiv v1 read 5 Nov 2026):

    dT/dt = D_T lap T + rho_n j_n^2 / C_v - (T - T0) / tau_e-ph,      j_n = -grad(phi) / rho_n,

with C_v constant (its normal-state value at T_c, as in the paper) and no latent-heat term.  The temperature enters
the TDGL through its linear coefficient.  pyTDGL normalises lengths and |psi| at the bath temperature T0, so
epsilon = (T_c - T) / (T_c - T0) = 1 - theta with theta = (T - T0) / (T_c - T0).  In pyTDGL units (time tau0, length
xi, normal current -grad mu):

    d theta / dt = D lap theta + eta |grad mu|^2 - theta / t_esc,       epsilon = max(epsilon_photon - theta, floor)

    eta = tau0 rho_n (K0/4)^2 / (C_v d^2 (T_c - T0)) = 4 E_cond / (C_v (T_c - T0)),
    E_cond = Phi0^2 / (16 pi^2 mu0 lambda^2 xi^2)  (the GL condensation energy density at T0).

rho_n cancels (tau0 = mu0 lambda^2 / rho_n), so eta needs only lambda, xi, C_v, T_c and T0.  For the film here
(lambda = 600 nm, xi = 5 nm: E_cond = 2.4 kJ/m^3) and C_v = 2.4 mJ cm^-3 K^-1 (NbN, the paper's value), eta = 0.59 at
T_c - T0 = 6.75 K (T_c = 9 K, T0 = 0.25 T_c).  It is uncertain by about a factor 2 either way.

The photon is the kinetic hotspot of snspd/hotspot.py, unchanged; Joule heat adds to it.  D and t_esc default to the
same values as that hotspot (10 nm^2/tau0, 100 tau0 = tau_esc of Simon et al.).

Our choices, not the paper's:
  * the photon is the kinetic hotspot (the paper uses an initial hot disc in the same heat equation);
  * heat sinks at the contacts: theta = 0 within sink_nm of the two terminals (the paper holds T = T0 at the contacts
    and raises T_c over 5 xi next to them).  pyTDGL injects the bias as normal current through psi = 0 terminals, so
    without a sink the current conversion next to them heats steadily;
  * epsilon >= floor (-1).  pyTDGL evaluates the psi update under np.errstate(all="raise"), and in a normal domain
    |psi| decays as exp(-|epsilon| t / u): with epsilon = -5, |psi|^3 underflows after ~270 tau0 and the solver stops
    (dt -> 0).  With -1 that takes > 1300 tau0.  psi is ~0 in such a domain either way; at the zipper's (0,0)
    threshold the floors -5 and -1 give the same outcomes (peak theta 0.195 / 0.194 in the healed run).
  * operator splitting: Joule source and escape are applied every solver step, diffusion (backward Euler, reflecting
    walls, the hotspot's FEM matrices) every dt_tau0 of accumulated time.

Detection becomes a single criterion: the run is still resistive at the end of the window (>= sim.SUSTAIN_SLIPS phase
slips in the last sim.SUSTAIN_WIN tau0).  Without heat this needs the extra continuation of sim.settle; with heat a
latched normal domain does not heal, which is the point of the model.
"""
from __future__ import annotations

import contextlib
import os
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .units import MU0, PHI0, Film


@dataclass(frozen=True)
class ThermalParams:
    eta: float = 0.59                 # Joule heating: d theta/dt per unit |grad mu|^2 (pyTDGL units), see module doc
    D_nm2_per_tau0: float = 10.0      # heat diffusion (= the hotspot's D; D_T = D for electrons, Wiedemann-Franz)
    esc_tau0: float = 100.0           # escape to the bath (= the kinetic hotspot's t_esc)
    sink_nm: float = 25.0             # theta = 0 within this distance of the terminals
    dt_tau0: float = 0.25             # diffusion step
    eps_floor: float = -1.0           # see module doc (runs of 5 Nov before 02:00 used -5: same outcomes at threshold)
    # two-temperature option (0 = one temperature, as Zotova & Vodolazov): the electrons (heated by Joule heat, theta
    # drives epsilon) exchange energy with the film phonons over ep_tau0; the phonons hold 1/ce_over_cph times more
    # heat per kelvin and alone escape to the substrate, over esc_tau0
    ep_tau0: float = 0.0              # electron -> phonon energy relaxation time tau_EP [tau0]
    ce_over_cph: float = 0.16         # c_e / c_ph at T_c


# 9 nm NbTiN on SiO2/Si (KTH films), measured by Sidorova et al., arXiv:2108.04717 (PRB 104, 184514 (2021)) and
# arXiv:2204.01423 (strips, Table 1): T_c 8.6 K (strip) / 9.51 K (film), R_sq 381 Ohm, D 0.472 cm^2/s, tau_EP(T_c)
# 2.9-4.0 ps, tau_esc 80 ps, c_e 1303 and c_ph 8173 J/(K m^3) at T_c.  With tau0 = hbar / (pi Delta0) = 0.15 ps
# (dirty BCS, T_c 8.6-9.5 K) and T0 = 0.3 T_c: eta = 4 E_cond / (c_e (T_c - T0)) = 1.22 (electrons only), D = 7.1
# nm^2/tau0, tau_EP = 27 tau0, tau_esc = 530 tau0.  All taken at T_c (c_ph ~ T^3, tau_EP ~ T^-3.4): a linearised,
# order-of-magnitude calibration.
NBTIN_2T = ThermalParams(eta=1.22, D_nm2_per_tau0=7.1, esc_tau0=530.0, ep_tau0=27.0, ce_over_cph=0.16)


def eta_estimate(film: Film = Film(), C_v_J_m3K: float = 2400.0, Tc_K: float = 9.0, T0_over_Tc: float = 0.25) -> float:
    """eta = 4 E_cond / (C_v (T_c - T0)), E_cond = Phi0^2 / (16 pi^2 mu0 lambda^2 xi^2)."""
    lam, xi = film.lambda_nm * 1e-9, film.xi_nm * 1e-9
    e_cond = PHI0 ** 2 / (16 * np.pi ** 2 * MU0 * lam ** 2 * xi ** 2)
    return 4.0 * e_cond / (C_v_J_m3K * Tc_K * (1.0 - T0_over_Tc))


class HeatField:
    """theta on the mesh sites of a pyTDGL device, advanced from the scalar potential mu after each solver step."""

    def __init__(self, device, params: ThermalParams):
        self.p = params
        mesh = device.mesh                                   # dimensionless (xi)
        xi_um = device.coherence_length.magnitude
        em = mesh.edge_mesh
        self.e0, self.e1 = em.edges[:, 0], em.edges[:, 1]
        self.w = em.dual_edge_lengths / em.edge_lengths      # FV weight: sum_e w (dmu)^2 = int |grad mu|^2 dA
        self.area = np.asarray(mesh.areas, float)
        n = len(self.area)
        # reflecting-wall diffusion on the triangles (same FEM as snspd.hotspot.HotspotSolver), lengths in xi
        P, T = np.asarray(mesh.sites, float), np.asarray(mesh.elements, int)
        p = P[T]
        x, y = p[..., 0], p[..., 1]
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
        c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
        a = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
        Kloc = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) / (4.0 * a)[:, None, None]
        rows = np.broadcast_to(T[:, :, None], Kloc.shape).ravel()
        cols = np.broadcast_to(T[:, None, :], Kloc.shape).ravel()
        K = sp.coo_matrix((Kloc.ravel(), (rows, cols)), shape=(n, n)).tocsc()
        self.m = np.bincount(T.ravel(), weights=np.repeat(a / 3.0, 3), minlength=n)
        D = params.D_nm2_per_tau0 / (xi_um * 1e3) ** 2       # xi^2 / tau0
        self._lu = spla.splu((sp.diags(self.m) + params.dt_tau0 * D * K).tocsc())
        # contact sinks
        y_um = P[:, 1] * xi_um
        half = np.abs(y_um).max()
        self.sink = np.abs(y_um) >= half - params.sink_nm * 1e-3
        self.theta = np.zeros(n)
        self.theta_ph = np.zeros(n)                          # phonon temperature (two-temperature option)
        self._tacc = 0.0
        self.peak = 0.0                                      # max theta outside the sinks so far
        self.joule = 0.0                                     # integral of eta |grad mu|^2 dA dt outside the sinks [xi^2]

    def __getstate__(self):
        # pyTDGL pickles disorder_epsilon (a closure over this object) into the solution file; keep that small
        return {"theta": self.theta, "theta_ph": self.theta_ph, "peak": self.peak, "joule": self.joule}

    def __setstate__(self, state):
        self.__dict__.update(state)

    def power(self, mu):
        """|grad mu|^2 per site (finite-volume, consistent with pyTDGL's edge currents)."""
        mu = np.asarray(mu).real
        pe = 0.5 * self.w * (mu[self.e0] - mu[self.e1]) ** 2
        return (np.bincount(self.e0, pe, len(self.area)) + np.bincount(self.e1, pe, len(self.area))) / self.area

    def step(self, dt, mu):
        q = self.p.eta * self.power(mu)
        q[self.sink] = 0.0
        self.joule += float(self.m @ q) * dt
        if self.p.ep_tau0 > 0:
            ex = (self.theta - self.theta_ph) / self.p.ep_tau0
            self.theta += dt * (q - ex)
            self.theta_ph += dt * (self.p.ce_over_cph * ex - self.theta_ph / self.p.esc_tau0)
        else:
            self.theta += dt * (q - self.theta / self.p.esc_tau0)
        self._tacc += dt
        while self._tacc >= self.p.dt_tau0:
            self.theta = self._lu.solve(self.m * self.theta)
            self._tacc -= self.p.dt_tau0
        self.theta[self.sink] = 0.0
        self.theta_ph[self.sink] = 0.0
        self.peak = max(self.peak, float(self.theta.max()))


def thermal_solver_class():
    import tdgl

    class ThermalTDGLSolver(tdgl.TDGLSolver):
        """tdgl.TDGLSolver with epsilon = max(eps_photon - theta, floor) and theta advanced after every step."""

        def __init__(self, device, options, *, thermal: ThermalParams, eps_photon=None, **kw):
            self.heat = HeatField(device, thermal)
            heat, floor = self.heat, thermal.eps_floor

            def eps(r, *, t, vectorized=True):
                base = np.ones(len(r)) if eps_photon is None else eps_photon(r, t=t)
                return np.maximum(base - heat.theta, floor)

            super().__init__(device, options, disorder_epsilon=eps, **kw)

        def update(self, state, running_state, dt, **kw):
            res = super().update(state, running_state, dt, **kw)
            self.heat.step(float(res.dt), res.mu)
            return res

    return ThermalTDGLSolver


@contextlib.contextmanager
def solved_thermal(device, current_uA, eps_photon, cfg, thermal: ThermalParams, seed=None):
    """Like snspd.sim.solved, with heat.  Yields (Solution, HeatField)."""
    from . import sim
    h5_path = sim.temp_h5(cfg)
    sol = None
    try:
        solver = thermal_solver_class()(device, sim.solver_options(cfg, cfg.tmax, h5_path), thermal=thermal,
                                        eps_photon=eps_photon, seed_solution=seed,
                                        terminal_currents=dict(source=current_uA, drain=-current_uA))
        sol = solver.solve()
        yield sol, solver.heat
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


def run_thermal(device, current_uA, eps_photon, cfg, thermal: ThermalParams, t0: float):
    """One run: dict(latched, trigger branch and time, slips in the last window, peak theta, area of theta >= 1 at
    the end [nm^2], Joule energy)."""
    from . import sim
    with solved_thermal(device, current_uA, eps_photon, cfg, thermal) as (sol, heat):
        trig = sim.detect(sol, t0)
        t = np.asarray(sol.dynamics.time)
        v = np.asarray(sol.dynamics.voltage())
        slips = float(sim.phase_slips(t, v, t0)[-1])
        xi_nm = device.coherence_length.magnitude * 1e3
        hot = float(heat.m[heat.theta >= 1.0].sum()) * xi_nm ** 2
        return dict(latched=trig.slips_last >= sim.SUSTAIN_SLIPS, branch=trig.branch, t_trigger=trig.t,
                    slips_last=trig.slips_last, slips=slips, theta_peak=heat.peak,
                    theta_end=float(heat.theta.max()), normal_nm2=hot, joule_xi2=heat.joule)
