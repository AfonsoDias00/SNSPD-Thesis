"""
Photon hotspot as a time-dependent suppression of the pyTDGL disorder parameter epsilon(r, t).

conserving (default)
    epsilon = 1 + A f exp(-(t - t0)/t_relax) n(r, t - t0),  dn/dt = D laplacian(n), no flux through film edges or
    hole walls, n(r, 0) = Gaussian of radius r0 rescaled so that  integral(n dA) = pi r0^2.  P1 finite elements on the
    device mesh, lumped mass, Crank-Nicolson.  In free space this is exactly the legacy profile; at walls the energy
    is reflected instead of lost.
legacy
    epsilon = 1 + A f (r0^2/r_t^2) exp(-(t - t0)/t_relax) exp(-|r - r_hit|^2 / r_t^2),  r_t^2 = r0^2 + 4 D (t - t0).
    Free-space Gaussian: energy that would sit outside the film or inside a hole is lost.

f scales the amplitude as a proxy for photon energy (f = 1 is the reference photon).  Hotspot picture: Semenov et
al., Physica C 351, 349 (2001); epsilon-suppression in TDGL: Zotova & Vodolazov, PRB 85, 024509 (2012).
Neither model has an energy sink other than exp(-t/t_relax) and neither solves a heat equation.
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

    def describe(self, film: Film) -> list[str]:
        D_si, tr, t0 = film.cm2_per_s(self.D_nm2_per_tau0), film.ps(self.t_relax_tau0), film.ps(self.t0_tau0)
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


def eps_conserving(mesh_points, snaps, amp, params: HotspotParams = HOTSPOT):
    """Vectorized pyTDGL epsilon(r, t) from precomputed diffusion snapshots."""
    mesh_points = np.asarray(mesh_points)
    K = snaps.shape[0]
    state = {"done": False, "idx": None}
    t0, tr, dtp, floor = params.t0_tau0, params.t_relax_tau0, params.dt_tau0, params.eps_floor

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
        n = (1.0 - f) * snaps[k].astype(float) + f * snaps[k + 1].astype(float)
        if state["idx"] is not None:
            n = n[state["idx"]]
        return np.maximum(1.0 + amp * np.exp(-dt / tr) * np.maximum(n, 0.0), floor)

    return eps_func


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
    if model != "conserving":
        raise ValueError(model)
    return eps_conserving(hs.P, hs.evolve(hs.initial_condition(x0, y0), tmax - params.t0_tau0), amp, params)


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
