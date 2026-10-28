"""
Device construction, one TDGL run with the detection criterion, threshold bisection, the current-crowding run and a
process-pool task runner.

Detection criterion (a run "switched" if ANY holds):
  1. peak |V| after t0 + 4 tau0 exceeds V_POST_MAX
  2. a V_TRAIL_WIN-long trailing mean of |V| ending at t >= V_TRAIL_START exceeds V_TRAIL_MEAN
  3. the supercurrent through a monitored cross-section of the geometry drops below COLLAPSE_REL * I_dep(strip width)
     at any saved frame after t0 + 4 tau0   (50 uA for the 200 nm strip with the default film)
Voltages are pyTDGL's normalised probe voltage (units of V0); the probes sit at y = +/- L/3.

Thresholds come from bisection to TOL_REL * I_dep(strip width) (4 uA for 200 nm, 0.0034 in i):
  I_sw   in ISW_BRACKET * I_dep(narrowest free width), no photon
  I_det  in [hit_lo_frac, 1] * I_sw
"""
from __future__ import annotations

import multiprocessing as mp
import os
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, replace

import numpy as np

from . import hotspot as hs_mod
from .geometry import Geometry
from .units import FILM, Film

try:
    import h5py
    import tdgl
    from tdgl.geometry import box, ellipse
    _TDGL_ERR = None
except ImportError as _exc:  # pragma: no cover  (analysis and reports work without the solver)
    h5py = tdgl = box = ellipse = None
    _TDGL_ERR = _exc

V_POST_MAX = 0.45        # peak |V| after t0 + 4 tau0
V_TRAIL_MEAN = 0.18      # trailing-mean |V| ...
V_TRAIL_WIN = 25.0       # ... over this many tau0 ...
V_TRAIL_START = 55.0     # ... for windows ending at t >= this
COLLAPSE_REL = 0.0422    # monitored supercurrent collapse threshold [fraction of I_dep(strip width)]
TOL_REL = 0.003376       # bisection tolerance [fraction of I_dep(strip width)]
ISW_BRACKET = (0.70, 1.20)   # dark-switching bracket [fraction of I_dep(narrowest free width)]
CROWD_TIME = 40.0        # steady-state run length for the current-crowding proxy [tau0]


@dataclass(frozen=True)
class SimConfig:
    tmax: float = 150.0          # solver window [tau0]
    mesh_xi: float = 1.4         # max mesh edge [xi]
    hotspot: str = "conserving"  # "conserving" or "legacy"
    save_every: int = 120
    dt_max: float = 0.20
    hit_lo_frac: float = 0.30    # lower photon-threshold bracket [fraction of I_sw]
    out_dir: str = "results/hitmap"

    @property
    def tag(self) -> str:
        return f"{self.hotspot}_T{self.tmax:g}_m{self.mesh_xi:g}"

    def path(self, stem: str, ext: str) -> str:
        return os.path.join(self.out_dir, f"{stem}_{self.tag}.{ext}")


def need_tdgl():
    if tdgl is None:
        raise SystemExit(f"py-tdgl / h5py could not be imported ({_TDGL_ERR}). Install with `pip install tdgl h5py`.")


def log(*a):
    print(*a, flush=True)


# ==============================================================================
# Device
# ==============================================================================
def make_device(geom: Geometry, mesh_xi: float, film: Film = FILM):
    need_tdgl()
    W, L = geom.width_nm * 1e-3, geom.length_nm * 1e-3
    holes = [tdgl.Polygon(f"hole_{i}", points=ellipse(h.a_nm * 1e-3, h.b_nm * 1e-3, center=(h.cx_nm * 1e-3, h.cy_nm * 1e-3),
                                                      points=geom.hole_points)).resample(geom.hole_points)
             for i, h in enumerate(geom.holes)]
    terminals = [tdgl.Polygon("source", points=box(W, L / 100, center=(0, L / 2))),
                 tdgl.Polygon("drain", points=box(W, L / 100, center=(0, -L / 2)))]
    device = tdgl.Device("snspd", layer=film.layer(), film=tdgl.Polygon("film", points=box(W, L)),
                         holes=holes or None, terminals=terminals, probe_points=[(0, L / 3), (0, -L / 3)],
                         length_units="um")
    device.make_mesh(max_edge_length=film.xi_nm * 1e-3 * mesh_xi, smooth=50)
    return device


_CTX = {}   # per-process cache: (geometry, mesh, film, hotspot params) -> device, monitored slices [um], hotspot solver


def context(geom: Geometry, mesh_xi: float, film: Film = FILM, params: hs_mod.HotspotParams = hs_mod.HOTSPOT):
    key = (geom, mesh_xi, film, params)
    if key not in _CTX:
        dev = make_device(geom, mesh_xi, film)
        _CTX[key] = dict(device=dev, slices=[(a * 1e-3, b * 1e-3, y * 1e-3) for a, b, y in geom.slices],
                         hs=hs_mod.HotspotSolver(dev.points, dev.triangles, params))
    return _CTX[key]


# ==============================================================================
# One run and the bisection
# ==============================================================================
def temp_h5(cfg: SimConfig) -> str:
    os.makedirs(cfg.out_dir, exist_ok=True)
    fd, path = tempfile.mkstemp(prefix="tmp_run-", suffix=".h5", dir=cfg.out_dir)
    os.close(fd)
    os.unlink(path)
    return path


def solver_options(cfg: SimConfig, solve_time: float, h5_path: str):
    return tdgl.SolverOptions(solve_time=solve_time, dt_init=1e-3, dt_max=cfg.dt_max, save_every=cfg.save_every,
                              field_units="mT", current_units="uA", output_file=h5_path)


def saved_steps_and_times(sol):
    with h5py.File(sol.path, "r") as f:
        groups = sorted(f["data"].keys(), key=int)
        times = np.array([f["data"][g].attrs["time"] for g in groups])
    return np.array([int(g) for g in groups]), times


def _trapz(y, dx):
    return (np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(y, dx=dx)


def trailing_mean_exceeds(t, v):
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


def slice_collapsed(sol, slices_um, collapse_uA, t_after):
    steps, times = saved_steps_and_times(sol)
    frames = [int(s) for s, tm in zip(steps, times) if tm > t_after]
    for x_min, x_max, y_sl in slices_um:
        x_pts = np.linspace(x_min, x_max, 30)
        pts = np.column_stack([x_pts, np.full_like(x_pts, y_sl)])
        dx = x_pts[1] - x_pts[0]
        for s in frames:
            sol.solve_step = s
            J = sol.interp_current_density(pts, dataset="supercurrent", units="uA/um")
            if _trapz(np.abs(J[:, 1]), dx) < collapse_uA:
                return True
    return False


def run_simulation(device, current_uA, slices_um, eps_func, cfg: SimConfig, collapse_uA: float,
                   t0: float = hs_mod.HOTSPOT.t0_tau0) -> bool:
    """One TDGL run at a fixed current. True if it counts as switching / detection."""
    h5_path = temp_h5(cfg)
    kw = {} if eps_func is None else {"disorder_epsilon": eps_func}
    sol = None
    try:
        sol = tdgl.solve(device, solver_options(cfg, cfg.tmax, h5_path),
                         terminal_currents=dict(source=current_uA, drain=-current_uA), **kw)
        t = np.asarray(sol.dynamics.time)
        v = np.abs(np.asarray(sol.dynamics.voltage()))
        post = t > t0 + 4.0
        if post.any() and v[post].max() > V_POST_MAX:
            return True
        if trailing_mean_exceeds(t, v):
            return True
        return slice_collapsed(sol, slices_um, collapse_uA, t0 + 4.0)
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


def bisect(device, slices_um, lo, hi, eps_func, cfg: SimConfig, tol_uA: float, collapse_uA: float):
    """Lowest switching current in [lo, hi] to tol_uA. Returns (I, flag, n_runs); flag 'at_wall' = never switched
    below hi, 'below_bracket' = switched even at lo."""
    lo0, hi0, n = lo, hi, 0
    while hi - lo > tol_uA:
        mid = 0.5 * (lo + hi)
        n += 1
        if run_simulation(device, mid, slices_um, eps_func, cfg, collapse_uA):
            hi = mid
        else:
            lo = mid
    flag = ""
    if hi >= hi0:
        flag = "at_wall"
    elif lo <= lo0:
        flag = "below_bracket"
    return hi, flag, n


# ==============================================================================
# Current crowding (steady state, no photon)
# ==============================================================================
def smoothed_peak(points_um, K, half_y_um, radius_um):
    """Peak of K averaged over a disc of radius_um inside |y| <= half_y_um.  Removes single-node spikes at boundaries
    without hiding real crowding at hole tips."""
    from scipy.spatial import cKDTree
    points_um = np.asarray(points_um)
    tree = cKDTree(points_um)
    sel = np.flatnonzero(np.abs(points_um[:, 1]) <= half_y_um)
    nb = tree.query_ball_point(points_um[sel], r=radius_um)
    Ks = np.array([K[ix].mean() for ix in nb])
    j = int(np.argmax(Ks))
    return float(Ks[j]), points_um[sel[j]]


def run_crowding(device, geom: Geometry, current_uA, cfg: SimConfig, film: Film = FILM):
    """Peak total sheet current K [uA/um] in the central unit cell after CROWD_TIME tau0 at a fixed current."""
    h5_path = temp_h5(cfg)
    sol = None
    half_y, radius = geom.crowd_half_y_nm * 1e-3, 1.5 * film.xi_nm * 1e-3
    try:
        sol = tdgl.solve(device, solver_options(cfg, CROWD_TIME, h5_path),
                         terminal_currents=dict(source=current_uA, drain=-current_uA))
        J = np.asarray(sol.current_density.to("uA/um").magnitude)
        K = np.linalg.norm(J, axis=1)
        pts = np.asarray(device.points)
        sel = np.abs(pts[:, 1]) <= half_y
        i = np.flatnonzero(sel)[int(np.argmax(K[sel]))]
        Ksm, loc = smoothed_peak(pts, K, half_y, radius)
        return dict(Kmax=float(K[i]), Ksm=Ksm, K99=float(np.quantile(K[sel], 0.99)), x_nm=float(pts[i, 0] * 1e3),
                    y_nm=float(pts[i, 1] * 1e3), xsm_nm=float(loc[0] * 1e3), ysm_nm=float(loc[1] * 1e3),
                    K_uniform=float(current_uA / (geom.width_nm * 1e-3))), (pts, np.asarray(device.triangles), K)
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


# ==============================================================================
# Tasks
# ==============================================================================
def scales(geom: Geometry, film: Film = FILM):
    """(tol_uA, collapse_uA, I_sw bracket) for a geometry."""
    I_strip = film.I_dep_uA(geom.width_nm)
    I_free = film.I_dep_uA(geom.free_width_min_nm())
    return TOL_REL * I_strip, COLLAPSE_REL * I_strip, (ISW_BRACKET[0] * I_free, ISW_BRACKET[1] * I_free)


def task_fn(task):
    """Run one task in this or a worker process. Never raises: failures come back with ok=False.

    task: dict(kind in {"isw", "hit", "crowd"}, geom=Geometry, cfg=SimConfig, film=Film, params=HotspotParams, ...)
      hit:   x_nm, y_nm, f, I_sw        crowd: bias, I_sw
    """
    t_start = time.time()
    try:
        geom, cfg = task["geom"], task["cfg"]
        film, params = task.get("film", FILM), task.get("params", hs_mod.HOTSPOT)
        ctx = context(geom, cfg.mesh_xi, film, params)
        tol, collapse, (lo, hi) = scales(geom, film)
        kind = task["kind"]
        if kind == "isw":
            val, flag, n = bisect(ctx["device"], ctx["slices"], lo, hi, None, cfg, tol, collapse)
            res = dict(I_sw=val, flag=("above_bracket" if flag == "at_wall" else flag), n_sims=n)
        elif kind == "hit":
            i_sw = task["I_sw"]
            eps = hs_mod.build_eps(ctx["hs"], cfg.hotspot, task["x_nm"], task["y_nm"], task["f"], cfg.tmax, params)
            val, flag, n = bisect(ctx["device"], ctx["slices"], cfg.hit_lo_frac * i_sw, i_sw, eps, cfg, tol, collapse)
            res = dict(I_det=val, flag=flag, n_sims=n, I_sw=i_sw)
        elif kind == "crowd":
            I = task["bias"] * task["I_sw"]
            res, (pts, tri, K) = run_crowding(ctx["device"], geom, I, cfg, film)
            res.update(I=I, bias=task["bias"], I_sw=task["I_sw"])
            np.savez_compressed(cfg.path(f"crowd_{geom.key}_{task['bias']:.2f}", "npz"),
                                points=pts, triangles=tri, K=K, K_uniform=res["K_uniform"])
        else:
            raise ValueError(kind)
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:  # keep the pool alive; the task is retried on resume
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def execute(tasks, workers, on_result, fn=task_fn):
    """Run tasks serially or in a spawn pool (each worker needs ~1.8 GB at 200 nm / 0.8 xi); on_result(result) is
    called in the parent as tasks finish."""
    total = len(tasks)
    if total == 0:
        return
    done = 0
    if workers <= 1:
        for task in tasks:
            on_result(fn(task))
            done += 1
            log(f"    [{done}/{total}]")
        return
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = "1"
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
        for fut in as_completed([pool.submit(fn, t) for t in tasks]):
            on_result(fut.result())
            done += 1
            log(f"    [{done}/{total}]")


def with_out_dir(cfg: SimConfig, out_dir: str) -> SimConfig:
    return replace(cfg, out_dir=out_dir)
