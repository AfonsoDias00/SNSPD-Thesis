"""
Device construction, one TDGL run with the detection criterion, threshold bisection, the current-crowding run and a
process-pool task runner.

Detection criterion: sustained switching.  A run "switched" if the probe voltage triggers, by either
  1. peak |V| after t0 + 4 tau0 exceeds V_POST_MAX
  2. a V_TRAIL_WIN-long trailing mean of |V| ending at t >= V_TRAIL_START exceeds V_TRAIL_MEAN,
and the device is still resistive at the end of the window: at least SUSTAIN_SLIPS phase slips between the probes in
the last SUSTAIN_WIN tau0 (|int V dt| / 2 pi; one vortex crossing the full width slips 1), and again at the end of a
SETTLE tau0 continuation from the final state (no photon; the hotspot has relaxed by then).  The continuation is needed
because near threshold a photon can trigger a flux rearrangement between holes (1-6 slips over up to ~100 tau0, then
quiet again) so late that it falls in the last SUSTAIN_WIN tau0 of the window.
Voltages are pyTDGL's normalised probe voltage (units of V0); the probes sit at y = +/- L/3, so any resistive belt
between them counts, wherever it forms.

History.  Until 29 Oct 2026 a third branch also counted a collapse of the supercurrent through one monitored
cross-section; it fired on local channel blocking without any voltage (a photon inside the zipper's monitored channel
"detected" at 0.724 vs 0.801 at its glide image) and was removed: results in results/v1_voltage_or_collapse/.  Until
29 Oct the voltage trigger alone counted; in the hole geometries it also fires on flux that a photon pushes from hole
to hole without the device staying resistive (fluxtrap.py): results in results/v2_voltage_trigger/.

Thresholds come from bisection to TOL_REL * I_dep(strip width) (0.34 uA for 200 nm, 0.0034 in i):
  I_sw   in ISW_BRACKET * I_dep(narrowest free width), no photon
  I_det  in [hit_lo_frac, 1] * I_sw
Each threshold also records which voltage branch fired earliest in the switching run at the threshold ("V_peak" or
"V_trail") and when (t_trigger, tau0 after t = 0; the photon arrives at t0): a latency figure of merit.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, replace
from typing import NamedTuple, Optional

import numpy as np

from . import hotspot as hs_mod
from .geometry import Geometry
from .units import FILM, Film

try:
    import h5py
    import tdgl
    from tdgl.geometry import box
    _TDGL_ERR = None
except ImportError as _exc:  # pragma: no cover  (analysis and reports work without the solver)
    h5py = tdgl = box = None
    _TDGL_ERR = _exc

V_POST_MAX = 0.45        # peak |V| after t0 + 4 tau0
V_TRAIL_MEAN = 0.18      # trailing-mean |V| ...
V_TRAIL_WIN = 25.0       # ... over this many tau0 ...
V_TRAIL_START = 55.0     # ... for windows ending at t >= this
SUSTAIN_WIN = 50.0       # sustained: at least SUSTAIN_SLIPS phase slips between the probes ...
SUSTAIN_SLIPS = 1.0      # ... in the last SUSTAIN_WIN tau0 of the window ...
SETTLE = 150.0           # ... and of a continuation this long [tau0]
TOL_REL = 0.003376       # bisection tolerance [fraction of I_dep(strip width)]
ISW_BRACKET = (0.70, 1.20)   # dark-switching bracket [fraction of I_dep(narrowest free width)]
CROWD_TIME = 40.0        # steady-state run length for the current-crowding proxy [tau0]
START_SETTLE = 100.0     # gamma = 1 dark run that prepares the settled start of a photon run [tau0]
SETTLE_TRIES = 3         # ... repeated (continuing) up to this many times until stationary


# ==============================================================================
# Parameter sets
# ==============================================================================
@dataclass(frozen=True)
class ParamSet:
    """Film, hotspot parameters and start protocol of one family of results (SimConfig.preset)."""
    film: Film
    params: hs_mod.HotspotParams
    settled_start: bool          # photon runs start from the steady current-carrying state (settled_state)
    note: str


PRESETS = {
    # every result up to 5 Nov 2026: gamma = 1, the kinetic hotspot with the NbN times of Simon et al.
    "": ParamSet(FILM, hs_mod.HOTSPOT, False,
                 "gamma = 1; hotspot times of NbN (Simon et al.): D 10 nm^2/tau0, escape 100 tau0; bias on at t = 0"),
    # 6 Nov 2026: the GTDGL inelastic parameter at pyTDGL's default and Simon et al.'s value (gamma = 2 tau_E Delta0 /
    # hbar = 12-18 for the measured tau_EP(T_c) = 2.9-4.0 ps of 9 nm NbTiN), and the hotspot times measured on 9 nm
    # NbTiN (Sidorova et al. 2021/22): D = 0.472 cm^2/s = 7.1 nm^2/tau0, energy escape 80 ps = 530 tau0 (tau0 = 0.15 ps).
    # xi 5 nm, lambda 600 nm, d 9 nm, u 5.79, r0, |A|, s_max and t_rise (unmeasured for NbTiN) are unchanged.
    "nbtin": ParamSet(replace(FILM, gamma=10.0),
                      replace(hs_mod.HOTSPOT, D_nm2_per_tau0=7.1, kin_esc_tau0=530.0), True,
                      "gamma = 10; NbTiN hotspot times: D 7.1 nm^2/tau0, escape 530 tau0; photon into the settled state"),
}


@dataclass(frozen=True)
class SimConfig:
    tmax: float = 150.0          # solver window [tau0]
    mesh_xi: float = 1.4         # max mesh edge [xi]
    hotspot: str = "kinetic"     # snspd.hotspot.MODELS; every result before 4 Nov 2026 is "conserving"
    save_every: int = 120
    dt_max: float = 0.025        # solver dt cap [tau0]: 0.2 (until 7 Nov 2026) left a local instability at small cells, see docs/methods.md
    hit_lo_frac: float = 0.30    # lower photon-threshold bracket [fraction of I_sw]
    sustained: bool = True       # False: the voltage trigger alone counts (the lower bound on the threshold)
    out_dir: str = "results/hitmap"
    preset: str = ""             # PRESETS key: film, hotspot parameters and start protocol ("" = the gamma = 1 results)

    @property
    def tag(self) -> str:
        model = self.hotspot + (f"-{self.preset}" if self.preset else "")
        return f"{model}_T{self.tmax:g}_m{self.mesh_xi:g}" + ("" if self.sustained else "_trigger")

    @property
    def pset(self) -> ParamSet:
        return PRESETS[self.preset]

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
    probe_y = L / 3 if geom.probe_y_nm is None else geom.probe_y_nm * 1e-3
    holes = [tdgl.Polygon(f"hole_{i}", points=h.outline(geom.hole_points) * 1e-3).resample(geom.hole_points)
             for i, h in enumerate(geom.holes)]
    terminals = [tdgl.Polygon("source", points=box(W, L / 100, center=(0, L / 2))),
                 tdgl.Polygon("drain", points=box(W, L / 100, center=(0, -L / 2)))]
    film_pts = box(W, L)
    if geom.notches:                       # edge notches: subtract them from the film outline
        from shapely.geometry import Polygon as SPoly
        outline = SPoly(film_pts)
        for n_ in geom.notches:
            outline = outline.difference(SPoly(n_.outline(geom.hole_points) * 1e-3))
        film_pts = np.asarray(outline.exterior.coords)[:-1]
    device = tdgl.Device("snspd", layer=film.layer(), film=tdgl.Polygon("film", points=film_pts),
                         holes=holes or None, terminals=terminals, probe_points=[(0, probe_y), (0, -probe_y)],
                         length_units="um")
    pinned = pinned_mesh_path(geom, mesh_xi)
    if pinned:
        d = np.load(pinned)
        device._create_dimensionless_mesh(np.asarray(d["points"], float), np.asarray(d["triangles"], int))
    else:
        device.make_mesh(max_edge_length=film.xi_nm * 1e-3 * mesh_xi, smooth=50)
    return device


def pinned_mesh_path(geom: Geometry, mesh_xi: float) -> Optional[str]:
    """A stored mesh to use instead of generating one, for reproducing results on another machine.

    Meshes of geometries with holes depend on the shapely / GEOS version (Polygon.resample uses GEOS segmentize and
    interpolate, and Triangle's refinement amplifies last-bit differences in the hole outlines): on 4 Nov 2026 a
    source-built GEOS 3.13.1 gave the long zipper 10 543 nodes instead of 10 710, and the (-68, 0) threshold moved
    from 68.80 to 72.03 uA.  Plain strips were identical.  Set SNSPD_MESH_DIR to a folder of
    mesh_<geometry key>_m<mesh_xi>.npz files (points in um, triangles), e.g. extracted from the crowding .npz files
    with `python -m snspd.sim pin-meshes`; geometries without a file there are meshed as usual."""
    d = os.environ.get("SNSPD_MESH_DIR")
    if not d:
        return None
    p = os.path.join(d, f"mesh_{geom.key}_m{mesh_xi:g}.npz")
    return p if os.path.exists(p) else None


_CTX = {}   # per-process cache: (geometry, mesh, film, hotspot params) -> device, hotspot solver


def context(geom: Geometry, mesh_xi: float, film: Film = FILM, params: hs_mod.HotspotParams = hs_mod.HOTSPOT):
    key = (geom, mesh_xi, film, params)
    if key not in _CTX:
        dev = make_device(geom, mesh_xi, film)
        _CTX[key] = dict(device=dev, hs=hs_mod.HotspotSolver(dev.points, dev.triangles, params))
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
    """SNSPD_NO_PROGRESS=1 turns pyTDGL's progress bar off (long background runs: the bars fill the logs)."""
    quiet = {"progress_interval": 10 ** 9} if os.environ.get("SNSPD_NO_PROGRESS") else {}
    return tdgl.SolverOptions(solve_time=solve_time, dt_init=1e-3, dt_max=cfg.dt_max, save_every=cfg.save_every,
                              field_units="mT", current_units="uA", output_file=h5_path, **quiet)


def saved_steps_and_times(sol):
    with h5py.File(sol.path, "r") as f:
        groups = sorted(f["data"].keys(), key=int)
        times = np.array([f["data"][g].attrs["time"] for g in groups])
    return np.array([int(g) for g in groups]), times


def trailing_mean_time(t, v) -> Optional[float]:
    """End time of the first V_TRAIL_WIN-long trailing mean of |V| above V_TRAIL_MEAN (ending at t >= V_TRAIL_START)."""
    dtg = 0.25
    if t[-1] < V_TRAIL_START:
        return None
    tg = np.arange(0.0, t[-1], dtg)
    vg = np.interp(tg, t, v)
    n = int(round(V_TRAIL_WIN / dtg))
    if len(vg) < n:
        return None
    cs = np.concatenate([[0.0], np.cumsum(vg)])
    ends = np.arange(n, len(vg) + 1)
    means = (cs[ends] - cs[ends - n]) / n
    hit = np.flatnonzero((tg[ends - 1] >= V_TRAIL_START) & (means > V_TRAIL_MEAN))
    return float(tg[ends[hit[0]] - 1]) if hit.size else None


def phase_slips(t, v, t0):
    """|int_t0^t V dt| / 2 pi on the solver time grid (signed V in V0, t in tau0): the phase slipped between the
    probes since t0, in units of 2 pi."""
    t, v = np.asarray(t), np.asarray(v)
    dv = 0.5 * (v[1:] + v[:-1]) * np.diff(t) * (t[1:] > t0)
    return np.abs(np.concatenate([[0.0], np.cumsum(dv)])) / (2 * np.pi)


class Trigger(NamedTuple):
    switched: bool                   # the criterion: a voltage trigger and sustained phase slips
    branch: str = ""                 # "V_peak" or "V_trail": the voltage branch that fired earliest ("" = none)
    t: Optional[float] = None        # when it fired [tau0]
    slips_last: float = 0.0          # phase slips in the last SUSTAIN_WIN tau0 of the window
    slips_settled: Optional[float] = None   # ... and of the SETTLE continuation (None: not run)


@contextmanager
def solved(device, current_uA, eps_func, cfg: SimConfig, seed=None, phase=None, init=None):
    """tdgl.solve at a fixed current for cfg.tmax; yields the Solution and deletes its file afterwards.  seed: a
    Solution of the same device to start from (default: the uniform superconducting state).  phase: a unit complex
    array over the mesh sites multiplying the uniform start instead (e.g. fluxoid.winding_psi).  init: (psi, mu) over
    the mesh sites to start from (settled_state).  At most one of the three."""
    if sum(v is not None for v in (seed, phase, init)) > 1:
        raise ValueError("give at most one of seed, phase and init")
    h5_path = temp_h5(cfg)
    kw = {} if eps_func is None else {"disorder_epsilon": eps_func}
    sol = None
    try:
        solver = tdgl.TDGLSolver(device, solver_options(cfg, cfg.tmax, h5_path), seed_solution=seed,
                                 terminal_currents=dict(source=current_uA, drain=-current_uA), **kw)
        if phase is not None:
            solver.psi_init = solver.psi_init * phase      # keeps psi = 0 on the terminals
        if init is not None:
            solver.psi_init, solver.mu_init = np.array(init[0]), np.array(init[1])
        sol = solver.solve()
        yield sol
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


def detect(sol, t0: float = hs_mod.HOTSPOT.t0_tau0) -> Trigger:
    """Apply the in-window part of the detection criterion to a finished run: which voltage branch fired first and
    when (also reported when the device did not stay resistive), and the phase slips at the end of the window.  A
    switch still needs settle()."""
    t = np.asarray(sol.dynamics.time)
    vs = np.asarray(sol.dynamics.voltage())
    v = np.abs(vs)
    last = float(phase_slips(t, vs, t[-1] - SUSTAIN_WIN)[-1])
    post = t > t0 + 4.0
    times = {}
    if post.any() and v[post].max() > V_POST_MAX:
        times["V_peak"] = float(t[post][np.argmax(v[post] > V_POST_MAX)])
    t_tr = trailing_mean_time(t, v)
    if t_tr is not None:
        times["V_trail"] = t_tr
    if not times:
        return Trigger(False, slips_last=last)
    branch = min(times, key=times.get)
    return Trigger(last >= SUSTAIN_SLIPS, branch, times[branch], last)


def settle(device, current_uA, sol, trig: Trigger, cfg: SimConfig) -> Trigger:
    """Confirm a switch found by detect(): continue from the final state of `sol` for SETTLE tau0 at the same current,
    without the photon, and require SUSTAIN_SLIPS phase slips in its last SUSTAIN_WIN tau0 as well."""
    if not trig.switched:
        return trig
    with solved(device, current_uA, None, replace(cfg, tmax=SETTLE), seed=sol) as ext:
        t = np.asarray(ext.dynamics.time)
        later = float(phase_slips(t, np.asarray(ext.dynamics.voltage()), t[-1] - SUSTAIN_WIN)[-1])
    return trig._replace(switched=later >= SUSTAIN_SLIPS, slips_settled=later)


_GAMMA1 = {}   # per-process cache: id(device) -> a gamma = 1 copy of it (same mesh) for settled_state


def settled_state(device, current_uA, cfg: SimConfig):
    """(psi, mu) of the steady current-carrying state at current_uA, or None if there is none (the device switches).

    The bias is switched on at t = 0 on the uniform state.  At gamma = 1 the strip has settled by the photon's t0 =
    20 tau0, but |psi| relaxes as u sqrt(1 + gamma^2 |psi|^2), so at gamma = 10 it takes ~200 tau0 (plain wire at 0.80
    I_dep: <|psi|^2> 0.953 at t0, 0.874 settled) and the photon would meet a condensate that is still stronger than
    the steady state.  The steady state is a stationary GL solution and does not depend on gamma, so it is prepared
    by a START_SETTLE tau0 dark run of a gamma = 1 copy of the device (checked on 6 Nov 2026: under gamma = 10 dynamics
    it stays put, <|psi|^2> unchanged to 1e-4 and |V| <= 0.011 over 200 tau0, for the plain wire and the long zipper
    at 0.80 I_dep).  Returns None if the state is still slipping phase (near or above the gamma = 1 dark I_sw) after
    SETTLE_TRIES x START_SETTLE tau0, first with gamma = 1 and then with the device's own gamma; the caller then
    starts from the uniform state.

    Stationarity is judged on the last quarter of the latest stretch.  Until 6 Nov 2026 it was judged on the second
    half of a single 100 tau0 run.  On the 2 um devices above ~0.92 I_dep the supercurrent phase between the probes
    still builds up there (Slot 400 at 0.93 I_dep: 0.07 x 2 pi over 50-100 tau0, stationary from ~120 tau0), so
    those runs fell back to the uniform start, whose bias-on transient fires V_trail at 55 tau0.  That affected the
    slot bridge-centre sites (thresholds >= 0.92) of the first nbtin maps.  The 0.8 um devices settle to 0.94 I_dep."""
    d1 = _GAMMA1.get(id(device))
    if d1 is None:
        L = device.layer
        d1 = device.copy(with_mesh=True)
        d1.layer = tdgl.Layer(coherence_length=L.coherence_length, london_lambda=L.london_lambda,
                              thickness=L.thickness, conductivity=L.conductivity, u=L.u, gamma=1.0)
        _GAMMA1[id(device)] = d1
    sub = replace(cfg, tmax=START_SETTLE)
    # gamma = 1 first (fast); if its bias step overshoots into switching (between the gamma = 1 and the gamma = 10
    # dark I_sw, ~0.96-0.98 I_dep: 7 Nov 2026, Single_400 bridge sites), the device's own gamma = 10 dynamics, which
    # settle more slowly but do not overshoot.  The stationary state is the same GL solution either way.
    for dev in (d1, device) if device.layer.gamma != 1.0 else (device,):
        seed = None
        with ExitStack() as stack:
            for _ in range(SETTLE_TRIES):
                sol = stack.enter_context(solved(dev, current_uA, None, sub, seed=seed))
                t, v = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
                if phase_slips(t, v, t[-1] - 0.25 * START_SETTLE)[-1] <= 0.02 and abs(v[-1]) <= 0.01:
                    return np.array(sol.tdgl_data.psi), np.array(sol.tdgl_data.mu)
                seed = sol
    return None


def run_simulation(device, current_uA, eps_func, cfg: SimConfig, t0: float = hs_mod.HOTSPOT.t0_tau0,
                   seed=None) -> Trigger:
    """One TDGL run at a fixed current: did it count as switching / detection, by which branch, and when.  With a
    preset that asks for it, a photon run starts from the settled state at that current (settled_state)."""
    init = None
    if seed is None and eps_func is not None and cfg.pset.settled_start:
        init = settled_state(device, current_uA, cfg)
    with solved(device, current_uA, eps_func, cfg, seed, init=init) as sol:
        trig = detect(sol, t0)
        return settle(device, current_uA, sol, trig, cfg) if cfg.sustained else trig._replace(switched=bool(trig.branch))


def bisect(device, lo, hi, eps_func, cfg: SimConfig, tol_uA: float, seed=None):
    """Lowest switching current in [lo, hi] to tol_uA. Returns (I, flag, n_runs, trigger), where trigger describes
    the switching run at I (empty if it never switched); flag 'at_wall' = never switched below hi, 'below_bracket' =
    switched even at lo.  seed: initial state of every run (see solved)."""
    lo0, hi0, n, trig = lo, hi, 0, Trigger(False)
    while hi - lo > tol_uA:
        mid = 0.5 * (lo + hi)
        n += 1
        r = run_simulation(device, mid, eps_func, cfg, seed=seed)
        if r.switched:
            hi, trig = mid, r
        else:
            lo = mid
    flag = ""
    if hi >= hi0:
        flag = "at_wall"
    elif lo <= lo0:
        flag = "below_bracket"
    return hi, flag, n, trig


# ==============================================================================
# Current crowding (steady state, no photon)
# ==============================================================================
def smoothed_peak(points_um, K, geom: Geometry, radius_um):
    """Peak of K averaged over discs of radius_um whose centres are at least radius_um from the film edges and the
    hole walls, inside |y| <= geom.crowd_half_y_nm.

    pyTDGL builds node currents from the adjacent edges, so single nodes, and boundary nodes in particular, carry
    +15-40 % spikes: a uniform strip shows a node maximum of 1.2-1.4 K_uniform at the film edge where the true value
    is 1.00.  The interior disc average gives 1.04 (0.8 xi) / 1.08 (1.4 xi) for the plain strip and 1.235 / 1.239 for
    the long zipper, i.e. it is mesh-converged and keeps the real crowding beside the holes."""
    from scipy.spatial import cKDTree
    points_um = np.asarray(points_um)
    r_nm = radius_um * 1e3
    x_nm, y_nm = points_um[:, 0] * 1e3, points_um[:, 1] * 1e3
    sel = np.flatnonzero((np.abs(y_nm) <= geom.crowd_half_y_nm) & (np.abs(x_nm) <= geom.width_nm / 2 - r_nm)
                         & ~geom.inside_hole(x_nm, y_nm, r_nm))
    nb = cKDTree(points_um).query_ball_point(points_um[sel], r=radius_um)
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
        Ksm, loc = smoothed_peak(pts, K, geom, radius)
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
    """(tol_uA, I_sw bracket) for a geometry."""
    I_free = film.I_dep_uA(geom.free_width_min_nm())
    return TOL_REL * film.I_dep_uA(geom.width_nm), (ISW_BRACKET[0] * I_free, ISW_BRACKET[1] * I_free)


def task_fn(task):
    """Run one task in this or a worker process. Never raises: failures come back with ok=False.

    task: dict(kind in {"isw", "hit", "crowd"}, geom=Geometry, cfg=SimConfig, film=Film, params=HotspotParams, ...)
      hit:   x_nm, y_nm, f, I_sw        crowd: bias, I_sw
    """
    t_start = time.time()
    try:
        geom, cfg = task["geom"], task["cfg"]
        film, params = task.get("film", cfg.pset.film), task.get("params", cfg.pset.params)
        if replace(film, rho_uohm_cm=None) != replace(cfg.pset.film, rho_uohm_cm=None):
            raise ValueError(f"task film {film} does not match preset {cfg.preset!r}")
        ctx = context(geom, cfg.mesh_xi, film, params)
        tol, (lo, hi) = scales(geom, film)
        kind = task["kind"]
        if kind == "isw":
            val, flag, n, trig = bisect(ctx["device"], lo, hi, None, cfg, tol)
            res = dict(I_sw=val, flag=("above_bracket" if flag == "at_wall" else flag), n_sims=n, branch=trig.branch,
                       t_trigger=trig.t)
        elif kind == "hit":
            i_sw = task["I_sw"]
            eps = hs_mod.build_eps(ctx["hs"], cfg.hotspot, task["x_nm"], task["y_nm"], task["f"], cfg.tmax, params)
            val, flag, n, trig = bisect(ctx["device"], cfg.hit_lo_frac * i_sw, i_sw, eps, cfg, tol)
            res = dict(I_det=val, flag=flag, n_sims=n, I_sw=i_sw, branch=trig.branch, t_trigger=trig.t)
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


# ==============================================================================
# Pinned meshes (see pinned_mesh_path)
# ==============================================================================
def pin_meshes(out_dir: str = "meshes", result_dirs=("results/hitmap", "results/slots")):
    """Write mesh_<key>_m<mesh_xi>.npz for every geometry whose crowding run stored its mesh, from the conserving
    files (the crowding run has no photon, so its mesh is the production mesh of that tag)."""
    import glob
    import re
    os.makedirs(out_dir, exist_ok=True)
    for d in result_dirs:
        for f in sorted(glob.glob(os.path.join(d, "crowd_*_0.80_conserving_T*_m*.npz"))):
            m = re.match(r"crowd_(.+)_0\.80_conserving_T[0-9.]+_m([0-9.]+)(_trigger)?\.npz", os.path.basename(f))
            if not m or m.group(3):
                continue
            dst = os.path.join(out_dir, f"mesh_{m.group(1)}_m{m.group(2)}.npz")
            if os.path.exists(dst):
                continue
            src = np.load(f)
            np.savez_compressed(dst, points=src["points"], triangles=src["triangles"])
            log(f"{f} -> {dst}  ({len(src['points'])} nodes)")


if __name__ == "__main__":
    import sys
    if sys.argv[1:2] == ["pin-meshes"]:
        pin_meshes(*sys.argv[2:3])
