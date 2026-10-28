#!/usr/bin/env python3
"""
photon_hitmap_v2.py -- plain wire vs. "long zipper", with an energy-conserving hotspot,
a dark-count proxy, a photon-energy sweep and a go/no-go report.

What is different from photon_hitmap.py (v1)
--------------------------------------------
1. Hotspot model (``--hotspot conserving``, default).
   v1 prescribes a free-space Gaussian and lets the film edges / hole walls simply cut it
   off, so energy is lost exactly where the comparison is decided (edge hits, hole walls,
   the 12 nm bridge).  v2 instead puts the same amount of energy (pi*r0^2, same exp(-dt/t_relax)
   decay, same amplitude) on the mesh and lets it diffuse with reflecting boundaries
   (P1 finite elements, Crank-Nicolson, natural Neumann BC at film edges AND hole walls).
   In free space this reproduces the v1 formula exactly; ``--hotspot legacy`` is the v1 formula
   (vectorized) so old results can be regenerated.
2. Window.  Default t_max = 150 tau_0 and the "sustained voltage" criterion is a *trailing*
   25 tau_0 mean instead of "mean over t>55", so lengthening the window can only add detections
   (v1's average was diluted by longer windows).  Normalisation uses I_sw measured with the same window.
3. Speed.  The epsilon callback is vectorized (``vectorized=True``).  py-tdgl otherwise calls it
   once per mesh site at every time step.  Optional process-level parallelism (``--workers``).
4. Dark-count proxies (see "What the proxies mean" below):
     a. photon-energy sweep: threshold i_det for weaker hotspots at probe sites;
     b. current crowding: peak local sheet current K_peak in the steady state before the photon.
5. Resumable CSV output, sampling extended to the film edge (x = -96 nm) and the bridge,
   hard-step efficiency curves (no invented tanh broadening), quantile metrics (i_50, i_90, i_95, i_99)
   and an explicit go / no-go table.

What the proxies mean (read this before quoting numbers)
--------------------------------------------------------
Deterministic TDGL has no fluctuations, so I_sw is an instability limit and not a dark-count rate.
  * Selectivity: a perturbation much weaker than a photon is a crude stand-in for a fluctuation.
    A geometry that meets your goal triggers on photons well below the wall but needs i -> 1 for weak seeds.
    selectivity = median[i_det(weak) - i_det(photon)] over the probe sites (larger = better).
  * K_peak at the 95 % operating point: if, at the bias where it detects 95 % of photons, the zipper's worst
    spot carries more sheet current than the plain wire's, it is expected to have the higher dark-count rate.
Neither is a dark-count rate.  They are cheap, and they can falsify a geometry, not prove one.

Usage
-----
  python photon_hitmap_v2.py selftest                    # numerics checks, energy-leak table, hotspot figure
  python photon_hitmap_v2.py plan   [options]            # how many solver runs, no simulation
  python photon_hitmap_v2.py run    [options]            # simulate (resumable) and then report
  python photon_hitmap_v2.py report [options]            # rebuild figures / summary from the CSV files

Useful options
  --grid {full,critical}  critical = ~10 well-chosen points per geometry (start here)
  --tmax 150              solver window in tau_0 (try 250 once to test convergence)
  --hotspot {conserving,legacy}
  --mesh-xi 1.4           max mesh edge in units of xi (try 0.8 on a few bridge points)
  --workers 4             parallel processes (each runs its own solver; set threads to 1)
  --energies 0.6,0.35     hotspot-energy factors for the probe sites (amp_spot is scaled)
  --only {both,plain,zipper}

Suggested order: selftest -> run --grid critical --hotspot legacy --tmax 150  (CONTROL: v1 hotspot with the v2
criterion; must reproduce your earlier T=150 sweep, otherwise the criterion rather than the hotspot moves the
thresholds) -> run --grid critical -> run --grid full -> run --grid critical --mesh-xi 0.8 (bridge points are only
~2 mesh cells wide at 1.4 xi) -> run --grid critical --tmax 250 (window convergence).
"""
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

try:  # the numerics, the report and the selftest of the numerics work without the solver
    import h5py
    import tdgl
    from tdgl.geometry import box, ellipse
    _TDGL_ERR = None
except ImportError as _exc:  # pragma: no cover
    h5py = tdgl = box = ellipse = None
    _TDGL_ERR = _exc

# ==============================================================================
# 1. Parameters (identical to v1 unless stated)
# ==============================================================================
W = 0.200              # wire width [um]
L = 0.800              # wire length [um]
d = 0.009              # film thickness [um]
xi = 0.005             # coherence length [um]
london_lambda = 0.175  # London penetration depth [um]

t0_photon = 20.0       # photon arrival [tau_0]
r_spot = 0.014         # initial hotspot radius [um]
amp_spot = -2.5        # reference core suppression (energy factor f = 1)
D_tdgl = 1.0e-5        # diffusion coefficient [um^2 / tau_0]  (= 10 nm^2 / tau_0)
t_relax = 35.0         # quasiparticle relaxation time [tau_0]
EPS_FLOOR = -6.0       # numerical guard: epsilon is never pushed below this

# detection criterion (same three branches as v1; branch 2 is now a trailing mean)
V_POST_MAX = 0.45      # peak |voltage| after t0 + 4
V_TRAIL_MEAN = 0.18    # trailing-mean |voltage| ...
V_TRAIL_WIN = 25.0     # ... over this many tau_0 ...
V_TRAIL_START = 55.0   # ... evaluated for window ends >= this time
COLLAPSE_UA = 50.0     # monitored channel supercurrent collapses below this [uA]
HIT_LO_FRAC = 0.70     # lower bisection bracket for photon thresholds [fraction of I_sw]

CROWD_TIME = 40.0      # steady-state run length for the current-crowding proxy [tau_0]
ANALYSIS_HALF_Y = 0.122  # crowding is evaluated in |y| <= this (central unit cell, away from the array ends)

# zipper geometry: (x_center, y_center); every hole has semi-axes (HOLE_A, HOLE_B)
HOLE_A, HOLE_B = 0.008, 0.055
HOLES = [(-0.036, y) for y in (-0.122, 0.0, 0.122)] + [(0.036, y) for y in (-0.183, -0.061, 0.061, 0.183)]
PERIOD_Y = 0.122

GEOMETRIES = {
    "Plain_Wire": dict(holes=False, slices=[(0.002, 0.098, 0.000)], isw_bracket=(1000.0, 1260.0)),
    "Long_Zipper": dict(holes=True, slices=[(0.046, 0.098, 0.061)], isw_bracket=(880.0, 1120.0)),
}
PROBES = {
    "Plain_Wire": [(-96, 0), (-58, 0), (0, 0)],
    "Long_Zipper": [(-96, 0), (-68, 0), (0, 0), (-18, 30), (-36, 61)],
}

# colours: categorical slots 1 and 2 of the reference palette (validated), neutral ink
C_PLAIN, C_ZIP = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SEQ_CMAP = LinearSegmentedColormap.from_list("seq_blue", ["#eaf2fc", "#9cc3ee", "#2a78d6", "#0e3a73"])
GEOM_COLOR = {"Plain_Wire": C_PLAIN, "Long_Zipper": C_ZIP}
GEOM_LABEL = {"Plain_Wire": "Plain wire", "Long_Zipper": "Long zipper"}


@dataclass
class Cfg:
    tmax: float = 150.0
    hotspot: str = "conserving"
    mesh_xi: float = 1.4
    tol: float = 4.0                 # bisection tolerance [uA]
    save_every: int = 120
    dt_max: float = 0.20
    energies: tuple = (0.6, 0.35)    # probe-site energy factors (photon level f = 1 is always run)
    grid: str = "full"
    workers: int = 1
    out_dir: str = "results_photon_hitmap_v2"
    do_crowding: bool = True
    do_probes: bool = True
    crowd_biases: tuple = (0.80, 0.90, 0.95)
    go_margin: float = 0.03          # zipper i_95 must beat plain by at least this
    sel_margin: float = 0.01         # zipper selectivity must beat plain by at least this
    only: str = "both"

    @property
    def tag(self) -> str:
        return f"{self.hotspot}_T{self.tmax:g}_m{self.mesh_xi:g}"

    def path(self, stem: str, ext: str) -> str:
        return os.path.join(self.out_dir, f"{stem}_{self.tag}.{ext}")


def log(*a):
    print(*a, flush=True)


# ==============================================================================
# 2. Geometry helpers (pure numpy, no solver needed)
# ==============================================================================
def inside_any_hole(x, y, margin=0.0):
    """True where (x, y) [um] lies inside (or within `margin` of) any zipper ellipse."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    out = np.zeros(np.broadcast(x, y).shape, bool)
    for cx, cy in HOLES:
        out |= ((x - cx) / (HOLE_A + margin)) ** 2 + ((y - cy) / (HOLE_B + margin)) ** 2 < 1.0
    return out


def free_width_min():
    """Smallest free width [um] of the zipper over one period (x-extent of film minus hole chords)."""
    ys = np.linspace(-PERIOD_Y / 2, PERIOD_Y / 2, 241)
    xs = np.linspace(-W / 2, W / 2, 2001)
    dx = xs[1] - xs[0]
    best = W
    for y in ys:
        best = min(best, float((~inside_any_hole(xs, y)).sum() * dx))
    return best


def hit_grid(geom: str, grid: str):
    """Photon-hit positions (x_nm, y_nm) simulated for one geometry."""
    if geom == "Plain_Wire":
        xs = [-96, -58, 0] if grid == "critical" else [-96, -88, -74, -58, -42, -28, -14, 0]
        return [(float(x), 0.0) for x in xs]
    if grid == "critical":
        pts = [(-96, 0), (-68, 0), (0, 0), (-50, 30), (-18, 30), (-36, 58), (-36, 61), (-42, 58), (-30, 58)]
    else:
        base = [(x, y) for x in (-96, -88, -68, -50, -18, 0) for y in (0, 15, 30, 45, 61)]
        bridge = [(x, y) for x in (-42, -36, -30) for y in (46, 50, 54, 58, 61)]
        pts = base + bridge
    keep = {(float(x), float(y)) for x, y in pts if not inside_any_hole(x * 1e-3, y * 1e-3, 0.002)}
    keep |= {(float(x), float(y)) for x, y in PROBES[geom]}
    return sorted(keep, key=lambda p: (p[1], p[0]))


# ==============================================================================
# 3. Energy-conserving hotspot: P1-FEM diffusion with reflecting boundaries
# ==============================================================================
class HotspotSolver:
    """Diffusion of the hotspot density n(r, t) on the device mesh.

        dn/dt = D * laplacian(n),   no flux through film edges or hole walls (natural BC),
        n(r, 0) = exp(-|r - r_hit|^2 / r0^2), rescaled so that  integral(n dA) = pi * r0^2.

    In free space n = (r0^2 / r_t^2) exp(-d^2 / r_t^2) with r_t^2 = r0^2 + 4 D t, i.e. the v1 profile.
    The integral is conserved exactly (lumped mass matrix, 1^T K = 0).  Crank-Nicolson, dt = 0.25 tau_0.
    """

    def __init__(self, points, triangles, D=D_tdgl, dt=0.25):
        self.P = np.ascontiguousarray(points, dtype=float)
        self.T = np.ascontiguousarray(triangles, dtype=int)
        self.D, self.dt = float(D), float(dt)
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
        self._tree = None

    def initial_condition(self, x0, y0, r0=r_spot):
        n0 = np.exp(-((self.P[:, 0] - x0) ** 2 + (self.P[:, 1] - y0) ** 2) / r0 ** 2)
        got = float(self.m @ n0)
        if got < 1e-12:
            raise ValueError(f"hotspot centre ({x0}, {y0}) is not on the mesh")
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

    def nearest(self, r):
        if self._tree is None:
            self._tree = cKDTree(self.P)
        return self._tree.query(np.asarray(r))[1]


def make_eps_conserving(mesh_points, snaps, amp, dtp):
    """epsilon(r, t) = 1 + amp * exp(-(t-t0)/t_relax) * n(r, t - t0), vectorized for py-tdgl."""
    mesh_points = np.asarray(mesh_points)
    K = snaps.shape[0]
    state = {"done": False, "idx": None}

    def eps_func(r, *, t, vectorized=True):
        r = np.asarray(r)
        dt = t - t0_photon
        if dt < 0.0:
            return np.ones(len(r))
        if not state["done"]:                      # solver passes device.points; otherwise map to nearest node
            same = r.shape == mesh_points.shape and np.allclose(r, mesh_points, atol=1e-9)
            if same:
                state["idx"] = None
            else:
                state["idx"] = cKDTree(mesh_points).query(r)[1]
            state["done"] = True
        s = min(dt / dtp, K - 1.0)
        k = min(int(s), K - 2)
        f = s - k
        n = (1.0 - f) * snaps[k].astype(float) + f * snaps[k + 1].astype(float)
        if state["idx"] is not None:
            n = n[state["idx"]]
        eps = 1.0 + amp * np.exp(-dt / t_relax) * np.maximum(n, 0.0)
        return np.maximum(eps, EPS_FLOOR)

    return eps_func


def make_eps_legacy(x0, y0, amp):
    """v1 free-space Gaussian (vectorized). Energy that falls outside the film / in holes is lost."""

    def eps_func(r, *, t, vectorized=True):
        r = np.asarray(r)
        dt = t - t0_photon
        if dt < 0.0:
            return np.ones(len(r))
        dist_sq = (r[:, 0] - x0) ** 2 + (r[:, 1] - y0) ** 2
        r0_sq = r_spot ** 2
        r_t_sq = r0_sq + 4.0 * D_tdgl * dt
        amp_t = amp * (r0_sq / r_t_sq) * np.exp(-dt / t_relax)
        return 1.0 + amp_t * np.exp(-dist_sq / r_t_sq)

    return eps_func


def legacy_energy_retention(x0, y0, with_holes, dts=(0.0, 10.0, 30.0, 60.0, 100.0)):
    """Fraction of the legacy hotspot's integrated density that lies inside the film (not outside the
    edges, not inside a hole).  1.0 means no energy lost.  Computed on a 1 nm grid, no mesh needed."""
    h = 0.001
    xs = np.arange(-W / 2 + h / 2, W / 2, h)
    out = []
    for dt in dts:
        rt2 = r_spot ** 2 + 4.0 * D_tdgl * dt
        half = 4.0 * np.sqrt(rt2)
        ys = np.arange(y0 - half, y0 + half, h)
        X, Y = np.meshgrid(xs, ys)
        n = (r_spot ** 2 / rt2) * np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / rt2)
        if with_holes:
            n = np.where(inside_any_hole(X, Y), 0.0, n)
        out.append(float(n.sum() * h * h / (np.pi * r_spot ** 2)))
    return list(dts), out


# ==============================================================================
# 4. Solver wrappers (need py-tdgl)
# ==============================================================================
def _need_tdgl():
    if tdgl is None:
        raise SystemExit(f"py-tdgl / h5py could not be imported ({_TDGL_ERR}). "
                         "Install with `pip install tdgl h5py`.")


def build_holes():
    _need_tdgl()
    return [tdgl.Polygon(f"hole_{i}", points=ellipse(HOLE_A, HOLE_B, center=(cx, cy), points=38)).resample(38)
            for i, (cx, cy) in enumerate(HOLES)]


def make_device(holes, mesh_xi):
    _need_tdgl()
    layer = tdgl.Layer(coherence_length=xi, london_lambda=london_lambda, thickness=d, gamma=1.0)
    film = tdgl.Polygon("film", points=box(W, L))
    terminals = [
        tdgl.Polygon("source", points=box(W, L / 100, center=(0, L / 2))),
        tdgl.Polygon("drain", points=box(W, L / 100, center=(0, -L / 2))),
    ]
    device = tdgl.Device(
        "snspd", layer=layer, film=film, holes=holes if len(holes) > 0 else None,
        terminals=terminals, probe_points=[(0, L / 3), (0, -L / 3)], length_units="um",
    )
    device.make_mesh(max_edge_length=xi * mesh_xi, smooth=50)
    return device


_CTX = {}  # per-process cache: (geometry, mesh_xi) -> device, slices, hotspot solver


def get_context(geom, mesh_xi):
    key = (geom, mesh_xi)
    if key not in _CTX:
        spec = GEOMETRIES[geom]
        dev = make_device(build_holes() if spec["holes"] else [], mesh_xi)
        _CTX[key] = dict(device=dev, slices=spec["slices"], hs=HotspotSolver(dev.points, dev.triangles))
    return _CTX[key]


def _trapz(y, dx):
    return (np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(y, dx=dx)


def saved_steps_and_times(sol):
    with h5py.File(sol.path, "r") as f:
        groups = sorted(f["data"].keys(), key=int)
        times = np.array([f["data"][g].attrs["time"] for g in groups])
    return np.array([int(g) for g in groups]), times


def trailing_mean_exceeds(t, v):
    """True if any trailing V_TRAIL_WIN-long mean of |V| ending at t >= V_TRAIL_START exceeds V_TRAIL_MEAN."""
    dtg = 0.25
    if t[-1] < V_TRAIL_START:
        return False
    tg = np.arange(0.0, t[-1], dtg)
    vg = np.interp(tg, t, v)
    n = int(round(V_TRAIL_WIN / dtg))
    if len(vg) < n:
        return False
    cs = np.concatenate([[0.0], np.cumsum(vg)])
    ends = np.arange(n, len(vg) + 1)
    means = (cs[ends] - cs[ends - n]) / n
    sel = tg[ends - 1] >= V_TRAIL_START
    return bool(sel.any() and means[sel].max() > V_TRAIL_MEAN)


def slice_collapsed(sol, slices):
    steps, times = saved_steps_and_times(sol)
    frames = [int(s) for s, tm in zip(steps, times) if tm > t0_photon + 4.0]
    for x_min, x_max, y_sl in slices:
        x_pts = np.linspace(x_min, x_max, 30)
        pts = np.column_stack([x_pts, np.full_like(x_pts, y_sl)])
        dx = x_pts[1] - x_pts[0]
        for s in frames:
            sol.solve_step = s
            J = sol.interp_current_density(pts, dataset="supercurrent", units="uA/um")
            if _trapz(np.abs(J[:, 1]), dx) < COLLAPSE_UA:
                return True
    return False


def _options(cfg, solve_time, h5_path):
    return tdgl.SolverOptions(solve_time=solve_time, dt_init=1e-3, dt_max=cfg.dt_max, save_every=cfg.save_every,
                              field_units="mT", current_units="uA", output_file=h5_path)


def _temp_h5(cfg):
    fd, h5_path = tempfile.mkstemp(prefix="tmp_v2-", suffix=".h5", dir=cfg.out_dir)
    os.close(fd)
    os.unlink(h5_path)
    return h5_path


def run_simulation(device, current_uA, slices, eps_func, cfg):
    """One TDGL run. Returns True if the run counts as a switching / detection event."""
    h5_path = _temp_h5(cfg)
    kw = {} if eps_func is None else {"disorder_epsilon": eps_func}
    sol = None
    try:
        sol = tdgl.solve(device, _options(cfg, cfg.tmax, h5_path),
                         terminal_currents=dict(source=current_uA, drain=-current_uA), **kw)
        t = np.asarray(sol.dynamics.time)
        v = np.abs(np.asarray(sol.dynamics.voltage()))
        post = t > t0_photon + 4.0
        if post.any() and v[post].max() > V_POST_MAX:
            return True
        if trailing_mean_exceeds(t, v):
            return True
        return slice_collapsed(sol, slices)
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


def smoothed_peak(points, K, radius=1.5 * xi):
    """Peak of K averaged over a disc of `radius` (default 1.5 xi) inside the analysis region.
    A single boundary node can carry a spurious spike (node currents are built from the adjacent edges only);
    averaging over a coherence-length-scale disc removes that without hiding real crowding at hole tips."""
    points = np.asarray(points)
    tree = cKDTree(points)
    sel = np.flatnonzero(np.abs(points[:, 1]) <= ANALYSIS_HALF_Y)
    nb = tree.query_ball_point(points[sel], r=radius)
    Ks = np.array([K[ix].mean() for ix in nb])
    j = int(np.argmax(Ks))
    return float(Ks[j]), points[sel[j]]


def run_crowding(device, current_uA, cfg):
    """Steady state without a photon: peak total sheet current K [uA/um] in the central unit cell."""
    h5_path = _temp_h5(cfg)
    sol = None
    try:
        sol = tdgl.solve(device, _options(cfg, CROWD_TIME, h5_path),
                         terminal_currents=dict(source=current_uA, drain=-current_uA))
        J = np.asarray(sol.current_density.to("uA/um").magnitude)          # (n_sites, 2)
        K = np.linalg.norm(J, axis=1)
        pts = np.asarray(device.points)
        sel = np.abs(pts[:, 1]) <= ANALYSIS_HALF_Y
        i = np.flatnonzero(sel)[int(np.argmax(K[sel]))]
        Ksm, loc = smoothed_peak(pts, K)
        return dict(Kmax=float(K[i]), Ksm=Ksm, K99=float(np.quantile(K[sel], 0.99)), x_nm=float(pts[i, 0] * 1e3),
                    y_nm=float(pts[i, 1] * 1e3), xsm_nm=float(loc[0] * 1e3), ysm_nm=float(loc[1] * 1e3),
                    K_uniform=float(current_uA / W)), (pts, np.asarray(device.triangles), K)
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


def bisect_threshold(device, slices, lo, hi, eps_func, cfg):
    lo0, hi0, n = lo, hi, 0
    while hi - lo > cfg.tol:
        mid = 0.5 * (lo + hi)
        n += 1
        if run_simulation(device, mid, slices, eps_func, cfg):
            hi = mid
        else:
            lo = mid
    flag = ""
    if hi >= hi0:
        flag = "at_wall"           # never switched below the upper bracket
    elif lo <= lo0:
        flag = "below_bracket"     # switched even at the lower bracket
    return hi, flag, n


def build_eps(ctx, x_nm, y_nm, amp, cfg):
    x0, y0 = x_nm * 1e-3, y_nm * 1e-3
    if cfg.hotspot == "legacy":
        return make_eps_legacy(x0, y0, amp)
    hs = ctx["hs"]
    return make_eps_conserving(
        hs.P,
        hs.evolve(hs.initial_condition(x0, y0), cfg.tmax - t0_photon),
        amp,
        hs.dt,
    )


def task_fn(task, cfg):
    """Run one task (module-level so it can be sent to worker processes). Never raises."""
    t_start = time.time()
    try:
        ctx = get_context(task["geom"], cfg.mesh_xi)
        kind = task["kind"]
        if kind == "isw":
            lo, hi = GEOMETRIES[task["geom"]]["isw_bracket"]
            val, flag, n = bisect_threshold(ctx["device"], ctx["slices"], lo, hi, None, cfg)
            res = dict(I_sw=val, flag=("above_bracket" if flag == "at_wall" else flag), n_sims=n)
        elif kind == "hit":
            i_sw = task["I_sw"]
            eps = build_eps(ctx, task["x_nm"], task["y_nm"], amp_spot * task["f"], cfg)
            val, flag, n = bisect_threshold(ctx["device"], ctx["slices"], HIT_LO_FRAC * i_sw, i_sw, eps, cfg)
            res = dict(I_det=val, flag=flag, n_sims=n, I_sw=i_sw)
        elif kind == "crowd":
            I = task["bias"] * task["I_sw"]
            res, arrays = run_crowding(ctx["device"], I, cfg)
            res.update(I=I, bias=task["bias"], I_sw=task["I_sw"])
            pts, tri, K = arrays
            np.savez_compressed(cfg.path(f"crowd_{task['geom']}_{task['bias']:.2f}", "npz"),
                                points=pts, triangles=tri, K=K, K_uniform=res["K_uniform"])
        else:
            raise ValueError(kind)
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:  # keep the pool alive; the task is simply not recorded and will be retried on resume
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


# ==============================================================================
# 5. Orchestration with resumable CSV
# ==============================================================================
CSV_FIELDS = ["Geometry", "x_nm", "y_nm", "energy_f", "tmax", "hotspot", "mesh_xi",
              "I_sw_uA", "I_det_uA", "i_det", "flag", "n_sims", "wall_s"]


def load_rows(cfg):
    path = cfg.path("hitmap", "csv")
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def row_key(geom, x, y, f):
    return (geom, round(float(x), 2), round(float(y), 2), round(float(f), 3))


def load_json(path, default):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return default


def save_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)


def selected_geoms(cfg):
    return [g for g in GEOMETRIES if cfg.only == "both" or g.lower().startswith(cfg.only)]


def execute(tasks, cfg, on_result):
    """Run tasks serially or in a spawn pool; call on_result(result) in the parent as they finish."""
    total = len(tasks)
    if total == 0:
        return
    done = 0
    if cfg.workers <= 1:
        for task in tasks:
            on_result(task_fn(task, cfg))
            done += 1
            log(f"    [{done}/{total}]")
        return
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    with ProcessPoolExecutor(max_workers=cfg.workers, mp_context=mp.get_context("spawn")) as pool:
        futs = [pool.submit(task_fn, t, cfg) for t in tasks]
        for fut in as_completed(futs):
            on_result(fut.result())
            done += 1
            log(f"    [{done}/{total}]")


def plan_tasks(cfg, isw):
    """All hit and crowding tasks still missing, most informative first."""
    have = {row_key(r["Geometry"], r["x_nm"], r["y_nm"], r["energy_f"]) for r in load_rows(cfg)}
    crowd = load_json(cfg.path("crowd", "json"), [])
    have_crowd = {(c["geom"], round(c["bias"], 2)) for c in crowd}
    tasks = []
    for g in selected_geoms(cfg):
        if g not in isw:
            continue
        probes = set(PROBES[g])
        for x, y in hit_grid(g, cfg.grid):
            prio = 0 if ((x, y) in probes or x <= -94 or (-43 <= x <= -29 and y >= 46)) else 1
            if row_key(g, x, y, 1.0) not in have:
                tasks.append(dict(kind="hit", geom=g, x_nm=x, y_nm=y, f=1.0, I_sw=isw[g]["I_sw"], prio=prio))
            if cfg.do_probes and (x, y) in probes:
                for f in cfg.energies:
                    if row_key(g, x, y, f) not in have:
                        tasks.append(dict(kind="hit", geom=g, x_nm=x, y_nm=y, f=f, I_sw=isw[g]["I_sw"], prio=1))
        if cfg.do_crowding:
            for b in cfg.crowd_biases:
                if (g, round(b, 2)) not in have_crowd:
                    tasks.append(dict(kind="crowd", geom=g, bias=b, I_sw=isw[g]["I_sw"], prio=0))
    tasks.sort(key=lambda t: t["prio"])
    return tasks


def cmd_plan(cfg):
    n_bis = int(np.ceil(np.log2(0.30 * 1000.0 / cfg.tol)))
    log(f"config: {cfg.tag}, grid={cfg.grid}, tol={cfg.tol} uA  ->  ~{n_bis} solver runs per bisection")
    fake_isw = {g: dict(I_sw=1000.0) for g in selected_geoms(cfg)}
    tasks = plan_tasks(cfg, fake_isw)
    n_hit = sum(t["kind"] == "hit" for t in tasks)
    n_crowd = sum(t["kind"] == "crowd" for t in tasks)
    n_sim = 2 * 7 + n_hit * n_bis + n_crowd
    log(f"missing: {n_hit} threshold searches + {n_crowd} crowding runs + I_sw  ~ {n_sim} solver runs of {cfg.tmax:g} tau_0")
    log("time one run first (e.g. `run --grid critical` prints per-task wall time) and divide by --workers.")


def cmd_run(cfg):
    _need_tdgl()
    os.makedirs(cfg.out_dir, exist_ok=True)
    log(f"== run: {cfg.tag}  grid={cfg.grid}  workers={cfg.workers}  out={cfg.out_dir}")
    isw_path = cfg.path("isw", "json")
    isw = load_json(isw_path, {})

    # phase A: dark switching current with this window
    todo = [dict(kind="isw", geom=g) for g in selected_geoms(cfg) if g not in isw]
    if todo:
        log(f"-- phase A: spontaneous I_sw for {[t['geom'] for t in todo]}")

        def on_isw(r):
            if r["ok"]:
                isw[r["task"]["geom"]] = r["res"]
                save_json(isw_path, isw)
                log(f"   {r['task']['geom']}: I_sw = {r['res']['I_sw']:.1f} uA  ({r['res']['n_sims']} runs, {r['wall_s']:.0f} s)"
                    + (f"  [{r['res']['flag']}]" if r["res"]["flag"] else ""))
            else:
                log(f"   FAILED {r['task']['geom']}: {r['error']}")

        execute(todo, cfg, on_isw)
    for g in selected_geoms(cfg):
        if g not in isw:
            raise SystemExit(f"I_sw for {g} is missing; cannot continue.")

    # phase B: hits, probes, crowding
    tasks = plan_tasks(cfg, isw)
    log(f"-- phase B: {len(tasks)} tasks")
    csv_path = cfg.path("hitmap", "csv")
    new_file = not os.path.exists(csv_path)
    crowd_path = cfg.path("crowd", "json")
    crowd = load_json(crowd_path, [])

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            log(f"   FAILED {t}: {r['error']}")
            return
        if t["kind"] == "hit":
            res = r["res"]
            row = dict(Geometry=t["geom"], x_nm=t["x_nm"], y_nm=t["y_nm"], energy_f=t["f"], tmax=cfg.tmax,
                       hotspot=cfg.hotspot, mesh_xi=cfg.mesh_xi, I_sw_uA=res["I_sw"], I_det_uA=res["I_det"],
                       i_det=res["I_det"] / res["I_sw"], flag=res["flag"], n_sims=res["n_sims"], wall_s=round(r["wall_s"], 1))
            nonlocal new_file
            with open(csv_path, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
                if new_file:
                    w.writeheader()
                    new_file = False
                w.writerow(row)
            log(f"   {t['geom']:<11} ({t['x_nm']:6.1f},{t['y_nm']:5.1f}) f={t['f']:.2f}  i_det = {row['i_det']:.3f}"
                f"{'  [' + res['flag'] + ']' if res['flag'] else ''}   ({res['n_sims']} runs, {r['wall_s']:.0f} s)")
        else:
            res = dict(r["res"], geom=t["geom"])
            crowd.append(res)
            save_json(crowd_path, crowd)
            log(f"   crowding {t['geom']:<11} bias {t['bias']:.2f}: K_peak = {res['Kmax']:.2f} uA/um at "
                f"({res['x_nm']:.0f},{res['y_nm']:.0f}) nm, K_peak/K_uniform = {res['Kmax'] / res['K_uniform']:.2f}")

    execute(tasks, cfg, on_result)
    cmd_report(cfg)


# ==============================================================================
# 6. Analysis
# ==============================================================================
def photon_rows(rows, geom, f=1.0):
    out = []
    for r in rows:
        if r["Geometry"] == geom and abs(float(r["energy_f"]) - f) < 1e-6:
            v = float(r["i_det"])
            if r["flag"] == "at_wall":
                v = 1.0
            out.append((float(r["x_nm"]), float(r["y_nm"]), v))
    return out


def symmetrize_zipper(samples):
    """y-mirror and glide symmetry (x, y) -> (-x, 61 - y) [mod period] of the zipper unit cell."""
    seen, pts, vals = set(), [], []

    def add(x, y, v):
        k = (round(x, 3), round(y, 3))
        if k not in seen:
            seen.add(k)
            pts.append((x, y))
            vals.append(v)

    for x, y, v in samples:
        for yy in (y, -y):
            add(x, yy, v)
        if x <= 0:
            yg = 61.0 - y
            for yy in (yg, -yg):
                add(-x, yy, v)
    return np.array(pts), np.array(vals)


XI_NM = np.linspace(-100, 100, 201)
YI_NM = np.linspace(-61, 61, 123)


def zipper_map(samples):
    pts, vals = symmetrize_zipper(samples)
    X, Y = np.meshgrid(XI_NM, YI_NM)
    Z = griddata(pts, vals, (X, Y), method="linear")
    Zn = griddata(pts, vals, (X, Y), method="nearest")
    Z = np.where(np.isnan(Z), Zn, Z)
    Z[inside_any_hole(X * 1e-3, Y * 1e-3)] = np.nan
    return X, Y, Z


def plain_map(samples):
    xs = np.array([abs(s[0]) for s in samples])
    vs = np.array([s[2] for s in samples])
    order = np.argsort(xs)
    xs, vs = xs[order], vs[order]
    X, Y = np.meshgrid(XI_NM, YI_NM)
    return X, Y, np.interp(np.abs(X), xs, vs)


def pixel_values(Z):
    return np.sort(Z[np.isfinite(Z)].ravel())


def eta_curve(vals, ib):
    return np.searchsorted(vals, ib, side="right") / len(vals) * 100.0


def metrics(vals):
    q = lambda p: float(np.quantile(vals, p))
    return dict(i_10=q(0.10), i_50=q(0.50), i_90=q(0.90), i_95=q(0.95), i_99=q(0.99), i_max=float(vals[-1]),
                width_10_90=q(0.90) - q(0.10))


def enrich_crowd(cfg, crowd):
    """Add the smoothed peak `Ksm` to entries written by an older run, recomputing it from the saved K maps."""
    for c in crowd:
        if "Ksm" in c:
            continue
        f = cfg.path(f"crowd_{c['geom']}_{c['bias']:.2f}", "npz")
        if os.path.exists(f):
            d_ = np.load(f)
            c["Ksm"], loc = smoothed_peak(d_["points"], d_["K"])
            c["xsm_nm"], c["ysm_nm"] = float(loc[0] * 1e3), float(loc[1] * 1e3)
        else:
            c["Ksm"] = c["Kmax"]
    return crowd


def K_of_I(crowd_entries, I, key="Ksm"):
    """Peak sheet current at bias I from the (bias -> K_peak) points, linear inter/extrapolation."""
    pts = sorted((c["I"], c[key]) for c in crowd_entries)
    if len(pts) < 2:
        return float("nan")
    xs, ys = zip(*pts)
    if xs[0] <= I <= xs[-1]:
        return float(np.interp(I, xs, ys))
    (x0, y0), (x1, y1) = (pts[0], pts[1]) if I < xs[0] else (pts[-2], pts[-1])
    return float(y0 + (y1 - y0) * (I - x0) / (x1 - x0))


def probe_table(rows, geom):
    tab = {}
    for r in rows:
        if r["Geometry"] == geom and (int(float(r["x_nm"])), int(float(r["y_nm"]))) in set(PROBES[geom]):
            v = 1.0 if r["flag"] == "at_wall" else float(r["i_det"])
            tab.setdefault((int(float(r["x_nm"])), int(float(r["y_nm"]))), {})[round(float(r["energy_f"]), 3)] = (v, r["flag"])
    return tab


def selectivity(tab):
    """median over probe sites of i_det(f) - i_det(1.0) for every weak factor f present."""
    fs = sorted({f for site in tab.values() for f in site if f < 0.999})
    out = {}
    for f in fs:
        diffs = [site[f][0] - site[1.0][0] for site in tab.values() if f in site and 1.0 in site]
        if diffs:
            out[f] = float(np.median(diffs))
    return out


# ==============================================================================
# 7. Figures
# ==============================================================================
def apply_style():
    plt.rcParams.update({
        "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11,
        "legend.frameon": False, "lines.linewidth": 2.0,
    })


def draw_holes(ax):
    th = np.linspace(0, 2 * np.pi, 120)
    for cx, cy in HOLES:
        ax.fill(1e3 * (cx + HOLE_A * np.cos(th)), 1e3 * (cy + HOLE_B * np.sin(th)), facecolor="#fcfcfb",
                edgecolor=INK, lw=0.9, zorder=3)


def fig_hitmaps(cfg, maps, samples, isw, path):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True, layout="constrained")
    allv = np.concatenate([pixel_values(m[2]) for m in maps.values()])
    vmin, vmax = np.floor(allv.min() * 200) / 200, np.ceil(allv.max() * 200) / 200
    im = None
    for ax, g in zip(axes, ["Plain_Wire", "Long_Zipper"]):
        if g not in maps:
            ax.axis("off")
            continue
        X, Y, Z = maps[g]
        im = ax.pcolormesh(X, Y, Z, cmap=SEQ_CMAP, vmin=vmin, vmax=vmax, shading="auto")
        if g == "Long_Zipper":
            draw_holes(ax)
        for x, y, _ in samples[g]:
            ax.plot(x, y, "o", ms=4.5, mfc="white", mec=INK, mew=0.8, zorder=4)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-61, 61)
        ax.set_aspect("equal")
        ax.set_title(f"{GEOM_LABEL[g]}\nI_sw = {isw[g]['I_sw']:.0f} uA,  mean i_det = {np.nanmean(Z):.3f}", fontsize=10)
        ax.set_xlabel("x0 [nm]")
        ax.grid(False)
    axes[0].set_ylabel("y0 [nm]  (one unit cell)")
    if im is not None:
        cb = fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
        cb.set_label("i_det = I_det / I_sw   (darker = less sensitive)")
    fig.suptitle(f"Photon-threshold hitmaps  [{cfg.hotspot} hotspot, T = {cfg.tmax:g} tau_0, f = 1]")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_efficiency(cfg, vals, met, isw, path):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    ib = np.linspace(0.75, 1.0, 600)
    for g, v in vals.items():
        c = GEOM_COLOR[g]
        a1.plot(ib, eta_curve(v, ib), color=c, label=f"{GEOM_LABEL[g]}:  i_50 = {met[g]['i_50']:.3f},  i_95 = {met[g]['i_95']:.3f}")
        a1.axvline(met[g]["i_95"], color=c, ls=":", lw=1.2)
        I = ib * isw[g]["I_sw"]
        a2.plot(I, eta_curve(v, ib), color=c, label=f"{GEOM_LABEL[g]}  (I_sw = {isw[g]['I_sw']:.0f} uA)")
        a2.axvline(isw[g]["I_sw"], color=c, ls="--", lw=1.0)
    a1.set_xlabel("normalised bias  i_b = I_bias / I_sw   (dotted: i_95)")
    a2.set_xlabel("absolute bias  I_bias [uA]   (dashed: I_sw)")
    for a in (a1, a2):
        a.set_ylabel("area-integrated efficiency [%]")
        a.set_ylim(-3, 103)
        a.legend(loc="lower right")
    a1.set_title("Normalised: window below the wall")
    a2.set_title("Absolute")
    fig.suptitle(f"Hard-step efficiency from the hitmap  [{cfg.hotspot}, T = {cfg.tmax:g} tau_0]")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def fig_selectivity(cfg, tabs, path):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), sharey=True)
    for ax, g in zip(axes, ["Plain_Wire", "Long_Zipper"]):
        tab = tabs.get(g, {})
        if not tab:
            ax.axis("off")
            continue
        fs = sorted({f for site in tab.values() for f in site})
        for site, d_ in sorted(tab.items()):
            xs = [f for f in fs if f in d_]
            ys = [d_[f][0] for f in xs]
            ax.plot(xs, ys, "-", color=INK2, lw=1.0, alpha=0.7)
            for f in xs:
                if d_[f][1] == "at_wall":
                    ax.plot(f, d_[f][0], "^", color=INK2, ms=7)
            ax.annotate(f"({site[0]},{site[1]})", (xs[0], ys[0]), textcoords="offset points", xytext=(-4, 0),
                        ha="right", fontsize=8, color=INK2)
        med = [np.median([d_[f][0] for d_ in tab.values() if f in d_]) for f in fs]
        ax.plot(fs, med, "-o", color=GEOM_COLOR[g], lw=2.6, ms=6, label="median over probe sites")
        ax.set_xlabel("hotspot energy factor f  (1 = photon)")
        ax.set_title(GEOM_LABEL[g])
        ax.legend(loc="lower right")
        ax.set_xlim(min(fs) - 0.12, 1.05)
    axes[0].set_ylabel("i_det   (triangle = no detection below I_sw)")
    fig.suptitle("Threshold vs hotspot energy: steeper rise towards 1 for weak seeds = more selective", y=0.99)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_crowding(cfg, geoms, path, bias):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), sharey=True, layout="constrained")
    ims = []
    data = {}
    for g in geoms:
        f = cfg.path(f"crowd_{g}_{bias:.2f}", "npz")
        if os.path.exists(f):
            data[g] = np.load(f)
    if not data:
        plt.close(fig)
        return
    for ax, g in zip(axes, ["Plain_Wire", "Long_Zipper"]):
        if g not in data:
            ax.axis("off")
            continue
        d_ = data[g]
        pts, tri, K = d_["points"], d_["triangles"], d_["K"] / float(d_["K_uniform"])
        ims.append(ax.tripcolor(pts[:, 0] * 1e3, pts[:, 1] * 1e3, tri, K, shading="gouraud", cmap=SEQ_CMAP, vmin=0.0, vmax=1.6))
        if g == "Long_Zipper":
            draw_holes(ax)
        sel = np.abs(pts[:, 1]) <= ANALYSIS_HALF_Y
        raw = float(K[sel].max())
        sm, loc = smoothed_peak(pts, K)
        ax.plot(loc[0] * 1e3, loc[1] * 1e3, "*", ms=13, mfc="white", mec=INK, zorder=5)
        ax.set_title(f"{GEOM_LABEL[g]}\nsmoothed peak {sm:.2f}x uniform  (single-node max {raw:.2f}x)", fontsize=10)
        ax.set_xlim(-100, 100)
        ax.set_ylim(-140, 140)
        ax.set_aspect("equal")
        ax.set_xlabel("x [nm]")
        ax.grid(False)
    axes[0].set_ylabel("y [nm]")
    fig.colorbar(ims[0], ax=axes, shrink=0.8, pad=0.02, label="|K| / (I / W)   (steady state before the photon)")
    fig.suptitle(f"Current crowding at {bias:.2f} I_sw  (star = smoothed peak, disc radius 1.5 xi)")
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ==============================================================================
# 8. Report and go / no-go
# ==============================================================================
def fmt(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"


def cmd_report(cfg):
    apply_style()
    os.makedirs(cfg.out_dir, exist_ok=True)
    rows = load_rows(cfg)
    isw = load_json(cfg.path("isw", "json"), {})
    crowd = enrich_crowd(cfg, load_json(cfg.path("crowd", "json"), []))
    geoms = [g for g in selected_geoms(cfg) if g in isw and photon_rows(rows, g)]
    if not geoms:
        log("nothing to report yet (no photon-level rows).")
        return
    samples = {g: photon_rows(rows, g) for g in geoms}
    maps = {"Plain_Wire": plain_map(samples["Plain_Wire"]) if "Plain_Wire" in geoms else None,
            "Long_Zipper": zipper_map(samples["Long_Zipper"]) if "Long_Zipper" in geoms else None}
    maps = {g: m for g, m in maps.items() if m is not None}
    vals = {g: pixel_values(m[2]) for g, m in maps.items()}
    met = {g: metrics(v) for g, v in vals.items()}
    tabs = {g: probe_table(rows, g) for g in geoms}
    sel = {g: selectivity(t) for g, t in tabs.items()}

    fig_hitmaps(cfg, maps, samples, isw, cfg.path("hitmap", "png"))
    fig_efficiency(cfg, vals, met, isw, cfg.path("efficiency", "png"))
    if any(t for t in tabs.values()):
        fig_selectivity(cfg, tabs, cfg.path("selectivity", "png"))
    for b in cfg.crowd_biases[-2:]:
        fig_crowding(cfg, geoms, cfg.path(f"crowding_{b:.2f}", "png"), b)

    wfree = free_width_min()
    lines = []
    P = lines.append
    P("=" * 78)
    P(f"SUMMARY   hotspot={cfg.hotspot}   T={cfg.tmax:g} tau_0   mesh={cfg.mesh_xi:g} xi   grid={cfg.grid}")
    if cfg.grid == "critical":
        P("NOTE: critical grid -> the hitmaps are interpolated from few points; quantiles are indicative only.")
    P("=" * 78)
    hdr = f"{'':<30}" + "".join(f"{GEOM_LABEL[g]:>14}" for g in geoms)
    P(hdr)
    P(f"{'I_sw [uA]':<30}" + "".join(f"{isw[g]['I_sw']:>14.1f}" for g in geoms))
    for key, name in [("i_10", "i_10"), ("i_50", "i_50"), ("i_90", "i_90"), ("i_95", "i_95"), ("i_99", "i_99"),
                      ("i_max", "i_max (100 % point)"), ("width_10_90", "10-90 % width")]:
        P(f"{name:<30}" + "".join(f"{met[g][key]:>14.3f}" for g in geoms))
    P(f"{'  (samples used)':<30}" + "".join(f"{len(samples[g]):>14d}" for g in geoms))

    if set(geoms) >= {"Plain_Wire", "Long_Zipper"}:
        pw, zp = met["Plain_Wire"], met["Long_Zipper"]
        P(f"{'onset gain  plain - zipper':<30}" + f"  i_10: {pw['i_10'] - zp['i_10']:+.3f}   i_50: {pw['i_50'] - zp['i_50']:+.3f}"
          f"   (positive = zipper turns on earlier)")
        P(f"{'plateau gain plain - zipper':<30}" + f"  i_90: {pw['i_90'] - zp['i_90']:+.3f}   i_95: {pw['i_95'] - zp['i_95']:+.3f}"
          f"   i_99: {pw['i_99'] - zp['i_99']:+.3f}")
        vz = vals["Long_Zipper"]
        P(f"{'zipper area above 0.87 / 0.90':<30}" + f"  {100 * np.mean(vz > 0.87):.1f} % / {100 * np.mean(vz > 0.90):.1f} %"
          f"   (the tail that sets i_95 / i_99)")

    # crowding / stress at the working point
    K95, Kfail, Kcrowd, Kraw = {}, {}, {}, {}
    for g in geoms:
        ent = [c for c in crowd if c["geom"] == g]
        if len(ent) >= 2:
            K95[g] = K_of_I(ent, met[g]["i_95"] * isw[g]["I_sw"])
            Kfail[g] = K_of_I(ent, isw[g]["I_sw"])
            top = max(ent, key=lambda c: c["bias"])
            Kcrowd[g] = top["Ksm"] / top["K_uniform"]
            Kraw[g] = top["Kmax"] / top["K_uniform"]
    if K95:
        P("-" * 78)
        P(f"current crowding (steady state; K_peak = max of K averaged over 1.5 xi; zipper free width in the busiest cell = {wfree * 1e3:.0f} nm)")
        P(f"{'K_peak / K_uniform @ top bias':<30}" + "".join(f"{fmt(Kcrowd.get(g), 2):>14}" for g in geoms))
        P(f"{'  (single-node max, noisy)':<30}" + "".join(f"{fmt(Kraw.get(g), 2):>14}" for g in geoms))
        P(f"{'K_peak at own i_95 [uA/um]':<30}" + "".join(f"{fmt(K95.get(g), 0):>14}" for g in geoms))
        P(f"{'K_peak at own I_sw [uA/um]':<30}" + "".join(f"{fmt(Kfail.get(g), 0):>14}" for g in geoms))
        P(f"{'I_sw / free width [uA/nm]':<30}" + "".join(
            f"{isw[g]['I_sw'] / (wfree * 1e3 if g == 'Long_Zipper' else W * 1e3):>14.2f}" for g in geoms))
        biases = sorted({c["bias"] for c in crowd})
        if len(biases) >= 2:
            for g in geoms:
                ent = sorted([c for c in crowd if c["geom"] == g], key=lambda c: c["bias"])
                g_lin = [c["Ksm"] / c["I"] for c in ent]
                P(f"  linearity {GEOM_LABEL[g]}: K/I changes by {100 * (g_lin[-1] / g_lin[0] - 1):+.1f} % from bias "
                  f"{ent[0]['bias']:.2f} to {ent[-1]['bias']:.2f}")
    # selectivity
    if any(sel.values()):
        P("-" * 78)
        P("selectivity  median[i_det(f) - i_det(1)] over probe sites (larger = weak seeds need more bias)")
        for f in sorted({f for s in sel.values() for f in s}):
            P(f"  f = {f:<5g}".ljust(30) + "".join(f"{fmt(sel[g].get(f)):>14}" for g in geoms))
        for g in geoms:
            nw = [(s, f) for s, d_ in tabs[g].items() for f, (v, fl) in d_.items() if fl == "at_wall"]
            if nw:
                P(f"  {GEOM_LABEL[g]}: no detection below I_sw at {len(nw)} (site, f) pairs -> counted as i_det = 1")

    # go / no-go
    P("-" * 78)
    P("GO / NO-GO   (thresholds are suggestions: --go-margin, --sel-margin)")
    verdict = []
    if set(geoms) >= {"Plain_Wire", "Long_Zipper"}:
        dz = met["Plain_Wire"]["i_95"] - met["Long_Zipper"]["i_95"]
        ok1 = dz >= cfg.go_margin
        P(f"  C1 detection : plain i_95 - zipper i_95 = {dz:+.3f}   need >= {cfg.go_margin:.3f}   -> {'PASS' if ok1 else 'FAIL'}")
        verdict.append(ok1)
        if "Plain_Wire" in K95 and "Long_Zipper" in K95:
            ok2 = K95["Long_Zipper"] <= K95["Plain_Wire"]
            P(f"  C2 stress    : K_peak@i_95  zipper {K95['Long_Zipper']:.0f} vs plain {K95['Plain_Wire']:.0f} uA/um"
              f"   need zipper <= plain   -> {'PASS' if ok2 else 'FAIL'}")
            verdict.append(ok2)
            P(f"                 K_peak@I_sw  zipper {Kfail['Long_Zipper']:.0f} vs plain {Kfail['Plain_Wire']:.0f} uA/um"
              "   (lower for the zipper = locally weaker spot)")
        else:
            P("  C2 stress    : n/a (no crowding data)")
        common = sorted(set(sel["Plain_Wire"]) & set(sel["Long_Zipper"]))
        if common:
            f = common[0]
            ds = sel["Long_Zipper"][f] - sel["Plain_Wire"][f]
            ok3 = ds >= cfg.sel_margin
            P(f"  C3 selectivity: zipper - plain at f = {f:g}: {ds:+.3f}   need >= {cfg.sel_margin:.3f}   -> {'PASS' if ok3 else 'FAIL'}")
            verdict.append(ok3)
        else:
            P("  C3 selectivity: n/a (no probe data for both geometries)")
        P(f"  OVERALL: {'GO (all evaluated criteria pass)' if verdict and all(verdict) else 'NO-GO / not shown'}"
          f"   [{sum(verdict)}/{len(verdict)} evaluated criteria passed]")
    else:
        P("  need both geometries for the comparison")
    P("-" * 78)
    P("caveats: deterministic TDGL; one photon energy (plus probes); sparse sampling; edge band |x|>96 nm is extrapolated;"
      " K_peak at hole tips depends on the mesh (check --mesh-xi 0.8).")
    text = "\n".join(lines)
    log(text)
    with open(cfg.path("summary", "txt"), "w") as f:
        f.write(text + "\n")
    save_json(cfg.path("summary", "json"), dict(metrics=met, selectivity={g: {str(k): v for k, v in s.items()} for g, s in sel.items()},
                                                K95=K95, Kfail=Kfail, I_sw={g: isw[g]["I_sw"] for g in geoms}))
    log(f"figures and summary written to {cfg.out_dir}/ (tag {cfg.tag})")


# ==============================================================================
# 9. Selftest
# ==============================================================================
def selftest_numerics(P, T, label, free_space=True):
    """Checks of the hotspot solver on a given mesh. Returns a dict of numbers and a list of PASS/FAIL lines."""
    t0 = time.time()
    hs = HotspotSolver(P, T)
    out, lines = {}, []

    def check(name, value, ok, note=""):
        lines.append(f"  [{'PASS' if ok else 'FAIL'}] {name}: {value}{('  ' + note) if note else ''}")

    n0 = hs.initial_condition(0.0, 0.0)
    snaps = hs.evolve(n0, 20.0)
    mass = snaps.astype(float) @ hs.m
    drift = float(abs(mass.max() - mass.min()) / mass[0])
    check("mass conserved", f"max drift {drift:.1e} (float32 storage)", drift < 1e-5)
    if free_space:
        r2 = (snaps.astype(float) * (hs.m * (P[:, 0] ** 2 + P[:, 1] ** 2))).sum(1) / mass
        t = np.arange(len(snaps)) * hs.dt
        slope = np.polyfit(t, r2, 1)[0] / (4 * D_tdgl)
        check("free-space <r^2> growth / 4D", f"{slope:.3f}", abs(slope - 1) < 0.08, "(1 = exact diffusion)")
        # compare with the analytic v1 profile at dt = 10
        k = int(round(10.0 / hs.dt))
        rt2 = r_spot ** 2 + 4 * D_tdgl * 10.0
        ana = (r_spot ** 2 / rt2) * np.exp(-(P[:, 0] ** 2 + P[:, 1] ** 2) / rt2)
        sel = np.abs(P[:, 1]) < 0.15
        err = float(np.linalg.norm(snaps[k][sel] - ana[sel]) / np.linalg.norm(ana[sel]))
        check("matches v1 profile in free space (dt=10)", f"rel L2 err {err:.3f}", err < 0.08)
        # reflecting edge: image solution for a hit 10 nm from the left edge
        x0 = -W / 2 + 0.010
        n0e = hs.initial_condition(x0, 0.0)
        se = hs.evolve(n0e, 10.0)
        ke = int(round(10.0 / hs.dt))
        xim = -W - x0
        g = lambda cx: (r_spot ** 2 / rt2) * np.exp(-((P[:, 0] - cx) ** 2 + P[:, 1] ** 2) / rt2)
        # initial truncated-and-rescaled Gaussian = Gaussian + its image (up to normalisation); compare shapes
        ana_e = g(x0) + g(xim)
        sel = np.abs(P[:, 1]) < 0.12
        a = se[ke][sel].astype(float)
        b = ana_e[sel]
        scale = float(a @ b / (b @ b))
        err_e = float(np.linalg.norm(a - scale * b) / np.linalg.norm(scale * b))
        check("reflecting edge = image solution (shape, dt=10)", f"rel L2 err {err_e:.3f}, rescale factor {scale:.2f}", err_e < 0.10,
              "(factor >1 : truncated initial Gaussian is renormalised to the full energy)")
        out.update(slope=slope, err_free=err, err_edge=err_e)
    out.update(mass_drift=drift, wall_s=time.time() - t0, n_sites=len(P))
    return out, [f"{label}: {len(P)} sites, {len(T)} triangles, {time.time() - t0:.1f} s"] + lines


def cmd_selftest(cfg):
    apply_style()
    os.makedirs(cfg.out_dir, exist_ok=True)
    log("== 1. How much energy does the v1 (legacy) hotspot lose?  (mesh-free, 1 nm grid)")
    log("   fraction of the hotspot's integrated density that is inside the film and outside the holes")
    sites = [("Plain_Wire", False, -96, 0), ("Plain_Wire", False, -88, 0), ("Plain_Wire", False, 0, 0),
             ("Long_Zipper", True, -96, 0), ("Long_Zipper", True, -68, 0), ("Long_Zipper", True, 0, 0),
             ("Long_Zipper", True, -36, 61), ("Long_Zipper", True, -42, 58)]
    dts = (0.0, 10.0, 30.0, 60.0, 100.0)
    log(f"   {'site':<26}" + "".join(f"  t0+{int(t):<4d}" for t in dts))
    for g, holes, x, y in sites:
        _, ret = legacy_energy_retention(x * 1e-3, y * 1e-3, holes, dts)
        log(f"   {g:<12} ({x:4d},{y:3d}) nm  " + "".join(f"  {r:6.2f} " for r in ret))
    log("   (the conserving model keeps 1.00 everywhere by construction)")

    if tdgl is None:
        log(f"\n== py-tdgl not importable ({_TDGL_ERR}); skipping the mesh-based checks.")
        return
    log("\n== 2. Hotspot solver on the real meshes")
    ctxs = {}
    for g in GEOMETRIES:
        ctxs[g] = get_context(g, cfg.mesh_xi)
        P, T = ctxs[g]["hs"].P, ctxs[g]["hs"].T
        _, lines = selftest_numerics(P, T, g, free_space=(g == "Plain_Wire"))
        log("\n".join(lines))
    log("\n== 3. Figure: legacy vs conserving epsilon for an edge hit and a bridge hit")
    cases = [("Plain_Wire", -96.0, 0.0), ("Long_Zipper", -36.0, 61.0)]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.5))
    for r_i, (g, x, y) in enumerate(cases):
        ctx = ctxs[g]
        hs = ctx["hs"]
        for c_i, mode in enumerate(("legacy", "conserving")):
            c2 = Cfg(**{**asdict(cfg), "hotspot": mode})
            eps = build_eps(ctx, x, y, amp_spot, c2)
            e = eps(hs.P, t=t0_photon + 10.0)
            ax = axes[r_i, c_i]
            im = ax.tripcolor(hs.P[:, 0] * 1e3, hs.P[:, 1] * 1e3, hs.T, 1.0 - e, shading="gouraud", cmap=SEQ_CMAP, vmin=0, vmax=2.6)
            if g == "Long_Zipper":
                draw_holes(ax)
            ax.plot(x, y, "*", ms=11, mfc="white", mec=INK, zorder=5)
            ax.set_xlim(-100, 100)
            ax.set_ylim(y - 90, y + 90)
            ax.set_aspect("equal")
            ax.set_title(f"{GEOM_LABEL[g]} hit ({x:.0f},{y:.0f})  -  {mode}", fontsize=10)
            ax.grid(False)
    fig.colorbar(im, ax=axes, shrink=0.8, label="1 - epsilon  at t0 + 10 tau_0")
    path = cfg.path("selftest_hotspot", "png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    log(f"   wrote {path}")


# ==============================================================================
# 10. CLI
# ==============================================================================
def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("Usage")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["selftest", "plan", "run", "report"])
    ap.add_argument("--tmax", type=float, default=Cfg.tmax)
    ap.add_argument("--hotspot", choices=["conserving", "legacy"], default=Cfg.hotspot)
    ap.add_argument("--mesh-xi", type=float, default=Cfg.mesh_xi)
    ap.add_argument("--tol", type=float, default=Cfg.tol, help="bisection tolerance [uA]")
    ap.add_argument("--save-every", type=int, default=Cfg.save_every)
    ap.add_argument("--grid", choices=["full", "critical"], default=Cfg.grid)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--energies", default="0.6,0.35", help="probe-site energy factors, comma separated")
    ap.add_argument("--crowd-biases", default="0.80,0.90,0.95")
    ap.add_argument("--no-crowding", action="store_true")
    ap.add_argument("--no-probes", action="store_true")
    ap.add_argument("--only", choices=["both", "plain", "zipper"], default="both")
    ap.add_argument("--go-margin", type=float, default=Cfg.go_margin)
    ap.add_argument("--sel-margin", type=float, default=Cfg.sel_margin)
    ap.add_argument("--out", default=Cfg.out_dir)
    a = ap.parse_args(argv)
    cfg = Cfg(tmax=a.tmax, hotspot=a.hotspot, mesh_xi=a.mesh_xi, tol=a.tol, save_every=a.save_every, grid=a.grid,
              workers=a.workers, energies=tuple(float(x) for x in a.energies.split(",") if x),
              crowd_biases=tuple(float(x) for x in a.crowd_biases.split(",") if x), do_crowding=not a.no_crowding,
              do_probes=not a.no_probes, only=a.only, go_margin=a.go_margin, sel_margin=a.sel_margin, out_dir=a.out)
    return a.command, cfg


def main(argv=None):
    cmd, cfg = parse_args(argv)
    if cfg.tmax <= t0_photon + 10:
        raise SystemExit("--tmax must exceed the photon time by a safe margin")
    {"selftest": cmd_selftest, "plan": cmd_plan, "run": cmd_run, "report": cmd_report}[cmd](cfg)


if __name__ == "__main__":
    main()