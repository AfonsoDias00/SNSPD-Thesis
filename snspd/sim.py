"""
The device, one electro-thermal TDGL event, the detection criterion and the threshold bisection.

One event (run_event):
  1. the strip is brought to its steady current-carrying state at the bias, at the bath temperature (settled_state);
  2. the photon's energy enters the electrons at one site (snspd/thermal.py);
  3. pyTDGL advances psi and mu; after every solver step the heat field is advanced with the dissipation of that
     step and returns epsilon(T_e) for the next one, and the circuit (snspd/circuit.py) is advanced with the voltage
     of that step and returns the detector current for the next one.

Detection criterion.  With the circuit, an event is a count when the load carries at least count_fraction (0.2) of
the bias current at some time: a pulse of at least 20 % of the full amplitude I_b R_L.  This is what the readout
sees.  It needs a resistive domain: a burst of phase slips that heals moves 2 pi I_u / N_sq per slip (0.05 uA for
750 squares) and stays far below the level.  An event that is still resistive but below the level at tmax runs on
until it counts, becomes quiet, or reaches tmax_ext.
Without the circuit (circuit=None: a current source, the limit N_sq -> infinity) a count is the criterion of the
earlier heated runs: still resistive at the end of the window, at least SUSTAIN_SLIPS phase slips between the
probes in the last SUSTAIN_WIN tau0.

Biases are in units of I_dep of the narrowest free cross-section (i = I / I_dep).  A bias with no steady dark state
(above the dark switching current) counts without a photon and is flagged "dark".

Voltages are differences of pyTDGL's mu between probes at (0, +/-probe_y), in V_u = Phi0 / (2 pi tau0); one vortex
crossing the full width, or one phase slip, is an integral of 2 pi.
"""
from __future__ import annotations

import bisect as _bisect
import itertools
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Optional

import numpy as np
import tdgl
from tdgl.geometry import box

from .circuit import Circuit, CircuitState
from .geometry import Geometry
from .thermal import NBTIN, HeatField, Photon, ThermalParams
from .units import FILM, Film

SUSTAIN_WIN = 50.0       # current-source criterion: at least SUSTAIN_SLIPS phase slips between the probes ...
SUSTAIN_SLIPS = 1.0      # ... in the last SUSTAIN_WIN tau0 of the window
QUIET_WIN = 25.0         # an event is quiet when less than QUIET_SLIPS slipped in the last QUIET_WIN tau0 ...
QUIET_SLIPS = 0.05       # ... and no site is above T_c
START_SETTLE = 100.0     # dark run that prepares the settled start [tau0]
SETTLE_TRIES = 3         # ... repeated (continuing) up to this many times until stationary
BRACKET = (0.50, 1.00)   # threshold bisection bracket [I_dep of the free width]
TOL = 0.005              # bisection tolerance [I_dep]


@dataclass(frozen=True)
class RunConfig:
    mesh_xi: float = 1.4             # max mesh edge [xi]; 2.0 for screening only (docs/methods.md)
    dt_max: float = 0.025            # solver dt cap [tau0]; larger leaves a local instability that heat amplifies
    tmax: float = 400.0              # window [tau0]; the photon arrives at photon.t_arrival
    tmax_ext: float = 800.0          # an undecided event (resistive, not yet a count) runs on until here
    film: Film = FILM
    thermal: ThermalParams = NBTIN
    circuit: Optional[Circuit] = Circuit()       # None: current source
    photon: Photon = Photon()
    count_fraction: float = 0.2      # load-current fraction that makes a count
    stop_on_count: bool = True       # threshold runs stop at the count; False records the whole window

    def describe(self) -> list[str]:
        L = self.film.describe() + self.thermal.describe(self.film)
        L += self.circuit.describe(self.film) if self.circuit else ["circuit: none (current source)"]
        ph = self.photon
        L.append(f"photon: {ph.wavelength_nm:g} nm ({ph.energy_eV:.3f} eV absorbed = "
                 f"{self.thermal.photon_energy(ph.energy_eV):.1f} energy units), r0 = {ph.r0_nm:g} nm, tau_th = "
                 f"{ph.tau_th_ps:g} ps ({self.film.tau0(ph.tau_th_ps):.1f} tau0), arrival {ph.t_arrival:g} tau0")
        L.append(f"run: mesh {self.mesh_xi:g} xi, dt_max {self.dt_max:g}, window {self.tmax:g} (to {self.tmax_ext:g}) "
                 f"tau0, count at {self.count_fraction:g} of I_b in the load")
        return L


def log(*a):
    print(*a, flush=True)


# ==============================================================================
# Device
# ==============================================================================
def make_device(geom: Geometry, mesh_xi: float, film: Film = FILM):
    W, L = geom.width_nm * 1e-3, geom.length_nm * 1e-3
    probe_y = L / 3 if geom.probe_y_nm is None else geom.probe_y_nm * 1e-3
    holes = [tdgl.Polygon(f"hole_{i}", points=h.outline(geom.hole_points) * 1e-3).resample(geom.hole_points)
             for i, h in enumerate(geom.holes)]
    terminals = [tdgl.Polygon("source", points=box(W, L / 100, center=(0, L / 2))),
                 tdgl.Polygon("drain", points=box(W, L / 100, center=(0, -L / 2)))]
    device = tdgl.Device("snspd", layer=film.layer(), film=tdgl.Polygon("film", points=box(W, L)),
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

    Meshes of geometries with holes depend on the shapely / GEOS version (a different GEOS build gave a holed strip
    10 543 nodes instead of 10 710 and moved one threshold by 0.04 I_dep; plain strips were identical).  Set
    SNSPD_MESH_DIR to a folder of mesh_<geometry key>_m<mesh_xi>.npz files, written by
    `python -m snspd.sim pin-meshes <dir> <mesh_xi>`; geometries without a file there are meshed as usual."""
    d = os.environ.get("SNSPD_MESH_DIR")
    if not d:
        return None
    p = os.path.join(d, f"mesh_{geom.key}_m{mesh_xi:g}.npz")
    return p if os.path.exists(p) else None


_CTX = {}   # per-process cache: (geometry, mesh, film) -> device and its scales


def context(geom: Geometry, cfg: RunConfig):
    key = (geom, cfg.mesh_xi, cfg.film)
    if key not in _CTX:
        dev = make_device(geom, cfg.mesh_xi, cfg.film)
        _CTX[key] = dict(device=dev, geom=geom, I_dep=cfg.film.I_dep_uA(geom.free_width_min_nm()),
                         probes=tuple(int(i) for i in dev.probe_point_indices))
    return _CTX[key]


# ==============================================================================
# Time stepping without pyTDGL's Runner
# ==============================================================================
class _NoRunningState:
    """Stands in for pyTDGL's RunningState: TDGLSolver.update appends dt and the probe values to it."""

    def append(self, name, value):
        pass


def march(solver, psi, mu, tmax, hook=None):
    """Advance a pyTDGL solver from (psi, mu) by tmax [tau0], calling TDGLSolver.update directly: the same steps as
    tdgl.solve (adaptive Euler, dt <= dt_max), but no output file, the final state comes back as arrays, and the run
    can stop early.  hook(t, dt, result) is called after every step with the time reached; returning True stops.
    Returns (psi, mu, t)."""
    ne = solver.num_edges
    values = {"psi": np.array(psi, dtype=complex), "mu": np.array(mu, dtype=float), "supercurrent": np.zeros(ne),
              "normal_current": np.zeros(ne), "induced_vector_potential": np.zeros((ne, 2))}
    if solver.dynamic_epsilon:
        values["epsilon"] = solver.epsilon
    names, rs, state = list(values), _NoRunningState(), {}
    t, dt = 0.0, solver.options.dt_init
    for i in itertools.count():
        state.update(step=i, time=t, dt=dt)
        res = solver.update(state, rs, dt, **values)
        dt = float(res[0])
        values = dict(zip(names, res[1:]))
        t += dt
        if (hook is not None and hook(t, dt, res)) or t >= tmax:
            break
    return values["psi"], values["mu"], t


def solver_options(cfg: RunConfig):
    return tdgl.SolverOptions(solve_time=1.0, dt_init=1e-3, dt_max=cfg.dt_max, save_every=10 ** 9,
                              field_units="mT", current_units="uA")


def phase_slips(t, v, t_from):
    """|integral of V dt from t_from| / 2 pi: the phase slipped between the probes, in units of 2 pi."""
    t, v = np.asarray(t), np.asarray(v)
    dv = 0.5 * (v[1:] + v[:-1]) * np.diff(t) * (t[1:] > t_from)
    return float(abs(dv.sum())) / (2 * np.pi)


_GAMMA1 = {}   # per-process cache: id(device) -> a gamma = 1 copy of it (same mesh)


def settled_state(ctx, current_uA, cfg: RunConfig):
    """(psi, mu) of the steady current-carrying state at current_uA, or None if there is none (dark switching).

    At gamma = 10, |psi| relaxes as u sqrt(1 + gamma^2 |psi|^2): after the bias step a strip needs ~200 tau0 to
    settle.  The steady state is a stationary GL solution and does not depend on gamma, so it is prepared by a
    START_SETTLE tau0 dark run of a gamma = 1 copy of the device (continued up to SETTLE_TRIES times), judged on the
    last quarter: at most 0.02 phase slips and |V_end| <= 0.01.  Between the gamma = 1 and gamma = 10 dark switching
    currents the gamma = 1 bias step overshoots into switching; then the device's own dynamics are tried.  Checked
    on the branch slot-geometries: the gamma = 10 dynamics leave the settled state in place (<|psi|^2> to 1e-4)."""
    device, film = ctx["device"], cfg.film
    d1 = _GAMMA1.get(id(device))
    if d1 is None:
        d1 = device.copy(with_mesh=True)
        d1.layer = film.layer(gamma=1.0)
        _GAMMA1[id(device)] = d1
    p0, p1 = ctx["probes"]
    for dev in (d1, device) if film.gamma != 1.0 else (device,):
        solver = tdgl.TDGLSolver(dev, solver_options(cfg), terminal_currents=dict(source=current_uA, drain=-current_uA))
        psi, mu = solver.psi_init, solver.mu_init
        for _ in range(SETTLE_TRIES):
            T, V = [], []

            def hook(t, dt, res):
                T.append(t)
                V.append(float(res.mu[p0] - res.mu[p1]))

            psi, mu, _t = march(solver, psi, mu, START_SETTLE, hook)
            if phase_slips(T, V, T[-1] - 0.25 * START_SETTLE) <= 0.02 and abs(V[-1]) <= 0.01:
                return np.array(psi), np.array(mu)
    return None


# ==============================================================================
# One event
# ==============================================================================
class ETSolver(tdgl.TDGLSolver):
    """tdgl.TDGLSolver with epsilon from the heat field and the terminal current from the circuit."""

    def __init__(self, device, options, heat: HeatField, circuit: Optional[CircuitState], current_uA: float, probes):
        self.heat, self.circuit, self.p0, self.p1 = heat, circuit, probes[0], probes[1]
        self.v = 0.0

        def eps(r, *, t, vectorized=True):
            return heat.eps

        cur = circuit.currents if circuit is not None else dict(source=current_uA, drain=-current_uA)
        super().__init__(device, options, terminal_currents=cur, disorder_epsilon=eps)

    def update(self, state, running_state, dt, **kw):
        psi_old, mu_old = kw["psi"], kw["mu"]
        res = super().update(state, running_state, dt, **kw)
        self.heat.step(state["time"], float(res.dt), psi_old, mu_old, res.psi, res.mu)
        self.v = float(res.mu[self.p0] - res.mu[self.p1])
        if self.circuit is not None:
            self.circuit.step(float(res.dt), self.v)
        return res


def run_event(ctx, i_bias: float, site, cfg: RunConfig, trace: bool = False, snap_times=()):
    """One photon at site (x, y) [nm] (None: no photon) with the strip biased at i_bias [I_dep of the free width].

    Returns a dict: count, t_count, load_max (peak load fraction), R_end (resistance at the end, in squares),
    resistive_end, slips, te_peak and tph_peak [T_c], normal_nm2 (largest area above T_c), the energy audit, and with
    trace=True the time series; snap_times [tau0] asks for |psi|^2, T_e and T_ph on the mesh at those times."""
    t_wall = time.time()
    device, film, th = ctx["device"], cfg.film, cfg.thermal
    I_b = i_bias * ctx["I_dep"]
    init = settled_state(ctx, I_b, cfg)
    if init is None:
        return dict(count=True, dark=True, i_bias=i_bias, t_count=None, load_max=None, R_end=None, resistive_end=None,
                    slips=None, te_peak=None, tph_peak=None, normal_nm2=None, t_end=0.0, wall_s=time.time() - t_wall)
    heat = HeatField(device, th, film, cfg.photon if site is not None else None, site)
    heat.reset(init[0])
    circ = CircuitState(cfg.circuit, film, I_b) if cfg.circuit is not None else None
    solver = ETSolver(device, solver_options(cfg), heat, circ, I_b, ctx["probes"])

    T, V, I, PHI, TE, TP, AN = [0.0], [0.0], [I_b], [0.0], [th.t0], [th.t0], [0.0]
    out = dict(t_count=None, load_max=0.0, normal_max=0.0)
    snaps, todo = [], sorted(snap_times)
    t_a = cfg.photon.t_arrival

    def slipped(t_from):
        k = max(_bisect.bisect_left(T, t_from) - 1, 0)
        return abs(PHI[-1] - PHI[k]) / (2 * np.pi)

    def hook(t, dt, res):
        cur = circ.I if circ is not None else I_b
        PHI.append(PHI[-1] + 0.5 * (V[-1] + solver.v) * dt)
        T.append(t), V.append(solver.v), I.append(cur)
        te, an = float(np.sqrt(heat.U[heat.live].max())), heat.normal_area_nm2()
        TE.append(te), TP.append(float(heat.P[heat.live].max() ** 0.25)), AN.append(an)
        out["normal_max"] = max(out["normal_max"], an)
        while todo and t >= todo[0]:
            todo.pop(0)
            snaps.append(dict(t=t, n=(np.abs(res.psi) ** 2).astype(np.float32), te=heat.te.astype(np.float32),
                              tph=heat.tph.astype(np.float32)))
        if circ is None:
            return False
        f = 1.0 - cur / I_b
        out["load_max"] = max(out["load_max"], f)
        if out["t_count"] is None and f >= cfg.count_fraction:
            out["t_count"] = t
        if out["t_count"] is not None and cfg.stop_on_count:
            return True
        if t >= cfg.tmax:
            if out["t_count"] is not None or t >= cfg.tmax_ext:
                return True
            return slipped(t - QUIET_WIN) < QUIET_SLIPS and te < 1.0      # quiet: stop; resistive: run on
        return False

    psi, mu, t_end = march(solver, init[0], init[1], cfg.tmax_ext if circ is not None else cfg.tmax, hook)
    T_, V_, I_ = np.array(T), np.array(V), np.array(I)
    last = T_ >= t_end - 5.0
    resistive = slipped(t_end - SUSTAIN_WIN) >= SUSTAIN_SLIPS
    count = out["t_count"] is not None if circ is not None else resistive
    e_el, e_ph = heat.energy()
    ev = dict(count=bool(count), dark=False, i_bias=i_bias, t_count=out["t_count"],
              load_max=out["load_max"] if circ is not None else None,
              R_end=float(np.mean(V_[last] * film.I_unit_uA / I_[last])), resistive_end=bool(resistive),
              slips=slipped(t_a), te_peak=heat.te_peak, tph_peak=heat.tp_peak, normal_nm2=out["normal_max"],
              I_end=float(I_[-1] / I_b), t_end=t_end, audit=dict(heat.audit, electrons=e_el, phonons=e_ph),
              wall_s=time.time() - t_wall)
    if trace:
        ev["trace"] = dict(t=T_, v=V_, I=I_, te=np.array(TE), tph=np.array(TP), normal_nm2=np.array(AN))
        ev["final"] = dict(psi=np.array(psi), mu=np.array(mu))
        ev["mesh"] = dict(points=np.asarray(device.points), triangles=np.asarray(device.triangles))
    if snaps:
        ev["snaps"] = snaps
    return ev


# ==============================================================================
# Thresholds
# ==============================================================================
def threshold(ctx, site, cfg: RunConfig, lo: float = BRACKET[0], hi: float = BRACKET[1], tol: float = TOL):
    """Lowest bias in [lo, hi] at which the event counts, by bisection to tol [I_dep].  site None: no photon (the
    dark switching current of the model).  Returns (i, flag, n_runs, event at i): flag "at_wall" = never counted,
    "below_bracket" = counted even at the lowest bias tried, "dark" = the count at i was dark switching."""
    lo0, hi0, n, best = lo, hi, 0, None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        n += 1
        ev = run_event(ctx, mid, site, cfg)
        if ev["count"]:
            hi, best = mid, ev
        else:
            lo = mid
    flag = "at_wall" if best is None else "below_bracket" if lo <= lo0 else "dark" if best["dark"] else ""
    return hi, flag, n, best


def task_fn(task):
    """Run one task in this or a worker process; never raises (failures come back with ok=False).

    task: dict(kind="threshold", geom=Geometry, cfg=RunConfig, site=(x, y) or None, [lo, hi, tol])
          dict(kind="event", geom, cfg, site, i_bias, [trace, snap_times])"""
    t_start = time.time()
    try:
        ctx = context(task["geom"], task["cfg"])
        if task["kind"] == "threshold":
            i, flag, n, ev = threshold(ctx, task["site"], task["cfg"], task.get("lo", BRACKET[0]),
                                       task.get("hi", BRACKET[1]), task.get("tol", TOL))
            res = dict(i_det=i, flag=flag, n_runs=n, event=ev)
        elif task["kind"] == "event":
            res = run_event(ctx, task["i_bias"], task["site"], task["cfg"], task.get("trace", False),
                            task.get("snap_times", ()))
        else:
            raise ValueError(task["kind"])
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start)
    except Exception as exc:  # keep the pool alive; the task is retried on resume
        import traceback
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", tb=traceback.format_exc(),
                    wall_s=time.time() - t_start)


def execute(tasks, workers, on_result, fn=task_fn):
    """Run tasks serially or in a spawn pool; on_result(result) is called in the parent as tasks finish."""
    total, done = len(tasks), 0
    if total == 0:
        return
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


# ==============================================================================
# Pinned meshes (see pinned_mesh_path)
# ==============================================================================
def pin_meshes(out_dir: str = "meshes", mesh_xi: float = 1.4):
    """Write mesh_<key>_m<mesh_xi>.npz for every registered geometry, from the mesh this machine generates."""
    from .geometry import GEOMETRIES
    os.makedirs(out_dir, exist_ok=True)
    for g in GEOMETRIES.values():
        dst = os.path.join(out_dir, f"mesh_{g.key}_m{mesh_xi:g}.npz")
        if os.path.exists(dst):
            continue
        dev = make_device(g, mesh_xi)
        np.savez_compressed(dst, points=np.asarray(dev.points), triangles=np.asarray(dev.triangles))
        log(f"{dst}  ({len(dev.points)} nodes)")


if __name__ == "__main__":
    import sys
    if sys.argv[1:2] == ["pin-meshes"]:
        pin_meshes(sys.argv[2] if len(sys.argv) > 2 else "meshes", float(sys.argv[3]) if len(sys.argv) > 3 else 1.4)
