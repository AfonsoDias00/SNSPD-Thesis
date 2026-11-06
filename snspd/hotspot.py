"""
Photon hotspot as a time-dependent suppression of the pyTDGL disorder parameter epsilon(r, t).

kinetic (default since 4 Nov 2026; after Simon et al., PRB 112, 174512 (2025), arXiv:2501.13791)
    The conserving model's diffusion of n(r, t) (below) carries the energy, but the suppression s = 1 - epsilon follows
    the local energy density e = f n exp(-(t - t0)/t_esc) with a delay, and saturates:
        ds/dt = (s_eq(e) - s) / t_rise,   s = 0 at t0,   s_eq(e) = s_max tanh(|A| e / s_max).
    From the paper:
      * pyTDGL's epsilon is the paper's alpha = (1 - T/Tc - eps_qp)/(1 - T/Tc) (Eq. 10), and the quasiparticle
        suppression builds up as eps_qp(t) = eps_inf (1 - exp(-t/tau_th)) (Sec. III).  For a constant drive this model
        gives exactly that curve with t_rise = tau_th: 1.4 ps (D = 0.5 cm^2/s) to 1.8 ps (1.5 cm^2/s), NbN, T = 0.2 Tc.
      * energy escapes to the substrate with tau_esc = 15 ps (Table 1; appendix VII.1 uses 10 ps).
    Not in the paper, which gives no spatial profile, no eps_inf, no bound on alpha and no wall treatment; the choices
    made here:
      * transport: the paper grows a disc of radius sqrt(D t); here the reflecting-wall diffusion of the conserving
        model carries the energy, which the hole geometries need.
      * amplitude: for weak density s_eq ~ |A| e, the conserving model's calibration, so f = 1 is the same photon.
      * s_max = 1.25, i.e. epsilon >= -0.25: eps_qp is the relative gap suppression in linear response, capped here at
        1 (|Delta| -> 0) at T = 0.2 Tc.  For f = 1 the suppression peaks at s = 0.6 (centre of a 200 nm strip) to 0.8
        (edge hit, energy reflected), so the cap is never reached; tanh only bends the response.
    Time unit: tau0 ~ hbar / (pi Delta0) = 0.12-0.17 ps for Tc = 8-12 K (README, "Time unit without a measured
    resistivity"), so t_rise = 10 tau0 = 1.2-1.7 ps and t_esc = 100 tau0 = 12-17 ps.  Anchoring tau0 on D = 10 nm^2/tau0
    = 0.5 cm^2/s instead (tau0 = 0.200 ps) gives 7 and 75 tau0.  These choices (and the cap) move photon
    thresholds by 0.01-0.04 I_dep, more for narrow strips under the trigger criterion (results/kinetic_sensitivity.txt).
    Result (compare_hotspots.py --profiles): for f = 1 no normal core forms (epsilon >= +0.43 at the centre of a
    200 nm strip; the conserving core starts at -1.5), but the suppression is broader and outlasts the conserving one:
    integral (1 - epsilon) dA is 1.3x the conserving value at t0 + 20, 2.8x at t0 + 50 and 7x at t0 + 100 tau0.
conserving (every result before 4 Nov 2026)
    epsilon = 1 + A f exp(-(t - t0)/t_relax) n(r, t - t0),  dn/dt = D laplacian(n), no flux through film edges or
    hole walls, n(r, 0) = Gaussian of radius r0 rescaled so that  integral(n dA) = pi r0^2.  P1 finite elements on the
    device mesh, lumped mass, Crank-Nicolson.  In free space this is exactly the legacy profile; at walls the energy
    is reflected instead of lost.
legacy
    epsilon = 1 + A f (r0^2/r_t^2) exp(-(t - t0)/t_relax) exp(-|r - r_hit|^2 / r_t^2),  r_t^2 = r0^2 + 4 D (t - t0).
    Free-space Gaussian: energy that would sit outside the film or inside a hole is lost.

f scales the amplitude as a proxy for photon energy (f = 1 is the reference photon); in the kinetic model it scales
the energy density, so the response saturates.  Hotspot picture: Semenov et al., Physica C 351, 349 (2001);
epsilon-suppression in TDGL: Zotova & Vodolazov, PRB 85, 024509 (2012).  No model here solves a heat equation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree

from .units import Film


@dataclass(frozen=True)
class HotspotParams:
    """Hotspot parameters in pyTDGL units (lengths in nm, times in tau0).  Defaults = every result in results/."""
    r0_nm: float = 14.0               # initial radius (2.8 xi)
    amp: float = -2.5                 # epsilon suppression at the core for f = 1 (epsilon_core = -1.5)
    D_nm2_per_tau0: float = 10.0      # diffusion constant
    t_relax_tau0: float = 35.0        # quasiparticle relaxation time
    t0_tau0: float = 20.0             # photon arrival, after the bias transient has settled
    dt_tau0: float = 0.25             # diffusion time step
    eps_floor: float = -6.0           # numerical guard: epsilon is never pushed below this
    # kinetic model only (see module docstring); t_relax is not used there
    kin_s_max: float = 1.25           # saturated suppression: epsilon >= 1 - s_max = -0.25 (a choice, see docstring)
    kin_rise_tau0: float = 10.0       # build-up time = tau_th (1.2-1.7 ps at tau0 = 0.12-0.17 ps)
    kin_esc_tau0: float = 100.0       # energy escape time = tau_esc (12-17 ps); np.inf = no loss

    def describe(self, film: Film, model: str = "conserving") -> list[str]:
        D_si, tr, t0 = film.cm2_per_s(self.D_nm2_per_tau0), film.ps(self.t_relax_tau0), film.ps(self.t0_tau0)
        if model == "kinetic":
            phys = "" if D_si is None else (f"  (= {D_si:.2f} cm^2/s, rise {film.ps(self.kin_rise_tau0):.2f} ps, "
                                            f"escape {film.ps(self.kin_esc_tau0):.2f} ps)")
            return [f"hotspot (kinetic): r0 = {self.r0_nm:g} nm, slope |A| = {abs(self.amp):g}, "
                    f"s_max = {self.kin_s_max:g}, D = {self.D_nm2_per_tau0:g} nm^2/tau0, rise = "
                    f"{self.kin_rise_tau0:g} tau0, escape = {self.kin_esc_tau0:g} tau0, t0 = {self.t0_tau0:g} tau0{phys}"]
        phys = "" if D_si is None else f"  (= {D_si:.2f} cm^2/s, {tr:.2f} ps, arrival {t0:.2f} ps)"
        return [f"hotspot: r0 = {self.r0_nm:g} nm, A = {self.amp:g}, D = {self.D_nm2_per_tau0:g} nm^2/tau0, "
                f"t_relax = {self.t_relax_tau0:g} tau0, t0 = {self.t0_tau0:g} tau0{phys}"]


HOTSPOT = HotspotParams()


class HotspotSolver:
    """Diffusion of the hotspot density n(r, t) on a triangle mesh (coordinates in um), reflecting boundaries."""

    def __init__(self, points, triangles, params: HotspotParams = HOTSPOT):
        self.params = params
        self.P = np.ascontiguousarray(points, dtype=float)
        self.T = np.ascontiguousarray(triangles, dtype=int)
        self.D, self.dt = params.D_nm2_per_tau0 * 1e-6, params.dt_tau0          # um^2 / tau0
        N = len(self.P)
        p = self.P[self.T]                                   # (M, 3, 2)
        x, y = p[..., 0], p[..., 1]
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1)
        c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1)
        area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
        Kloc = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) / (4.0 * area)[:, None, None]
        rows = np.broadcast_to(self.T[:, :, None], Kloc.shape).ravel()
        cols = np.broadcast_to(self.T[:, None, :], Kloc.shape).ravel()
        self.K = sp.coo_matrix((Kloc.ravel(), (rows, cols)), shape=(N, N)).tocsc()
        self.m = np.bincount(self.T.ravel(), weights=np.repeat(area / 3.0, 3), minlength=N)  # lumped mass
        M = sp.diags(self.m)
        h = 0.5 * self.dt * self.D
        self._lu = spla.splu((M + h * self.K).tocsc())
        self._B = (M - h * self.K).tocsr()

    def initial_condition(self, x0_um, y0_um):
        r0 = self.params.r0_nm * 1e-3
        n0 = np.exp(-((self.P[:, 0] - x0_um) ** 2 + (self.P[:, 1] - y0_um) ** 2) / r0 ** 2)
        got = float(self.m @ n0)
        if got < 1e-12:
            raise ValueError(f"hotspot centre ({x0_um}, {y0_um}) um is not on the mesh")
        return n0 * (np.pi * r0 ** 2 / got)

    def evolve(self, n0, t_end):
        steps = max(2, int(np.ceil(t_end / self.dt)))
        out = np.empty((steps + 1, len(n0)), dtype=np.float32)
        n = np.array(n0, float)
        out[0] = n
        for k in range(steps):
            n = self._lu.solve(self._B @ n)
            out[k + 1] = n
        return out


def _snapshot_interpolator(mesh_points, snaps, params: HotspotParams, transform):
    """pyTDGL epsilon(r, t) = transform(value interpolated linearly in time from snaps, t - t0); 1 before t0."""
    mesh_points = np.asarray(mesh_points)
    K = snaps.shape[0]
    state = {"done": False, "idx": None}
    t0, dtp = params.t0_tau0, params.dt_tau0

    def eps_func(r, *, t, vectorized=True):
        r = np.asarray(r)
        dt = t - t0
        if dt < 0.0:
            return np.ones(len(r))
        if not state["done"]:                      # the solver passes device.points; otherwise map to nearest node
            same = r.shape == mesh_points.shape and np.allclose(r, mesh_points, atol=1e-9)
            state["idx"] = None if same else cKDTree(mesh_points).query(r)[1]
            state["done"] = True
        s = min(dt / dtp, K - 1.0)
        k = min(int(s), K - 2)
        f = s - k
        v = (1.0 - f) * snaps[k].astype(float) + f * snaps[k + 1].astype(float)
        if state["idx"] is not None:
            v = v[state["idx"]]
        return transform(v, dt)

    return eps_func


def eps_conserving(mesh_points, snaps, amp, params: HotspotParams = HOTSPOT):
    """Vectorized pyTDGL epsilon(r, t) from precomputed diffusion snapshots."""
    tr, floor = params.t_relax_tau0, params.eps_floor
    return _snapshot_interpolator(
        mesh_points, snaps, params,
        lambda n, dt: np.maximum(1.0 + amp * np.exp(-dt / tr) * np.maximum(n, 0.0), floor))


def kinetic_suppression(snaps, f, params: HotspotParams = HOTSPOT):
    """Convert diffusion snapshots n_k (in place, float32) into epsilon_k of the kinetic model.

    s_eq = s_max tanh(|A| e / s_max) with e = f n exp(-t/t_esc) at the snapshots, and ds/dt = (s_eq - s)/t_rise
    integrated exactly for s_eq varying linearly over each step (first-order-hold exponential integrator; for
    t_rise -> 0 it returns s = s_eq exactly).  Returns snaps (now epsilon)."""
    dt, s_max, tau = params.dt_tau0, params.kin_s_max, params.kin_rise_tau0
    gain = abs(params.amp) * f / s_max
    q = np.exp(-dt / tau)
    c = (tau / dt) * (1.0 - q)       # weight of the change of s_eq within the step
    esc = params.kin_esc_tau0

    def s_eq(k):
        loss = np.exp(-k * dt / esc) if np.isfinite(esc) else 1.0
        return s_max * np.tanh(gain * np.maximum(snaps[k].astype(float), 0.0) * loss)

    s = np.zeros(snaps.shape[1])     # no suppression at t0: it builds up from there
    eq_prev = s_eq(0)
    snaps[0] = 1.0
    for k in range(1, snaps.shape[0]):
        eq = s_eq(k)
        s = eq + (s - eq_prev) * q - (eq - eq_prev) * c
        snaps[k] = 1.0 - s
        eq_prev = eq
    return snaps


def eps_kinetic(mesh_points, snaps, f, params: HotspotParams = HOTSPOT):
    """Vectorized pyTDGL epsilon(r, t) of the kinetic model from precomputed diffusion snapshots (consumed)."""
    return _snapshot_interpolator(mesh_points, kinetic_suppression(snaps, f, params), params, lambda e, dt: e)


def eps_legacy(x0_um, y0_um, amp, params: HotspotParams = HOTSPOT):
    """Vectorized free-space Gaussian (energy outside the film or inside holes is lost)."""
    t0, tr, r0_sq, D = params.t0_tau0, params.t_relax_tau0, (params.r0_nm * 1e-3) ** 2, params.D_nm2_per_tau0 * 1e-6

    def eps_func(r, *, t, vectorized=True):
        r = np.asarray(r)
        dt = t - t0
        if dt < 0.0:
            return np.ones(len(r))
        dist_sq = (r[:, 0] - x0_um) ** 2 + (r[:, 1] - y0_um) ** 2
        r_t_sq = r0_sq + 4.0 * D * dt
        amp_t = amp * (r0_sq / r_t_sq) * np.exp(-dt / tr)
        return 1.0 + amp_t * np.exp(-dist_sq / r_t_sq)

    return eps_func


def build_eps(hs: HotspotSolver, model: str, x_nm, y_nm, f, tmax, params: HotspotParams = HOTSPOT):
    """epsilon(r, t) for a photon of energy factor f absorbed at (x_nm, y_nm)."""
    x0, y0, amp = x_nm * 1e-3, y_nm * 1e-3, params.amp * f
    if model == "legacy":
        return eps_legacy(x0, y0, amp, params)
    if model not in MODELS:
        raise ValueError(model)
    snaps = hs.evolve(hs.initial_condition(x0, y0), tmax - params.t0_tau0)
    if model == "kinetic":
        return eps_kinetic(hs.P, snaps, f, params)
    return eps_conserving(hs.P, snaps, amp, params)


MODELS = ("conserving", "legacy", "kinetic")


def legacy_energy_retention(geom, x_nm, y_nm, dts=(0.0, 10.0, 30.0, 60.0, 100.0), params: HotspotParams = HOTSPOT):
    """Fraction of the legacy hotspot's integrated density inside the film and outside the holes (1 = nothing lost),
    on a 1 nm grid, no mesh needed."""
    h = 1.0
    xs = np.arange(-geom.width_nm / 2 + h / 2, geom.width_nm / 2, h)
    out = []
    for dt in dts:
        rt2 = params.r0_nm ** 2 + 4.0 * params.D_nm2_per_tau0 * dt
        half = 4.0 * np.sqrt(rt2)
        X, Y = np.meshgrid(xs, np.arange(y_nm - half, y_nm + half, h))
        n = (params.r0_nm ** 2 / rt2) * np.exp(-((X - x_nm) ** 2 + (Y - y_nm) ** 2) / rt2)
        n = np.where(geom.inside_hole(X, Y), 0.0, n)
        out.append(float(n.sum() * h * h / (np.pi * params.r0_nm ** 2)))
    return list(dts), out
