"""
Shared pieces of the validation suite: the result record, devices that are not hit-map geometries (square films,
an annulus, a strip with one hole), runs without terminals, bisection with a frozen epsilon, and line integrals of
the phase and of the sheet current.

Every test reports Checks: one number each, the value it is expected to have, and PASS / FAIL where the expectation
is sharp, INFO where the number characterises the model without a pass mark.
"""
from __future__ import annotations

import json
import os
import tempfile
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace

import numpy as np

from snspd import geometry as geo
from snspd import sim
from snspd.hotspot import HOTSPOT
from snspd.units import FILM

OUT_DIR = "results/validation"
XI = FILM.xi_nm
CFG = sim.SimConfig(tmax=400.0, mesh_xi=1.4, out_dir=OUT_DIR)      # production window, mesh and criterion


@dataclass
class Check:
    id: str            # e.g. "1a.3"
    name: str          # what is measured
    value: str         # the measured number, formatted
    expected: str      # expected value or rule
    verdict: str       # PASS, FAIL or INFO
    note: str = ""


def check(id_, name, value, expected, ok=None, note=""):
    return Check(id_, name, value, expected, "INFO" if ok is None else ("PASS" if ok else "FAIL"), note)


def out_path(name):
    os.makedirs(OUT_DIR, exist_ok=True)
    return os.path.join(OUT_DIR, name)


def save_group(gid, title, checks, wall_s, extra=None):
    with open(out_path(f"{gid}.json"), "w") as fh:
        json.dump(dict(id=gid, title=title, wall_s=round(wall_s, 1), checks=[asdict(c) for c in checks],
                       extra=extra or {}), fh, indent=1)


def load_groups():
    out = {}
    if os.path.isdir(OUT_DIR):
        for f in sorted(os.listdir(OUT_DIR)):
            if f.endswith(".json") and f[:1].isdigit() and len(f.split(".")[0]) == 2:
                with open(os.path.join(OUT_DIR, f)) as fh:
                    out[f[:2]] = json.load(fh)
    return out


def load_cache():
    p = out_path("cache.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def save_cache(c):
    with open(out_path("cache.json"), "w") as fh:
        json.dump(c, fh, indent=1)


# ==============================================================================
# Devices
# ==============================================================================
def film_device(film_nm, holes_nm=(), mesh_xi=1.4, name="film"):
    """tdgl.Device without terminals or probes from polygons in nm (film outline, hole outlines)."""
    import tdgl
    holes = [tdgl.Polygon(f"hole_{i}", points=np.asarray(h) * 1e-3) for i, h in enumerate(holes_nm)]
    dev = tdgl.Device(name, layer=FILM.layer(), film=tdgl.Polygon("film", points=np.asarray(film_nm) * 1e-3),
                      holes=holes or None, length_units="um")
    dev.make_mesh(max_edge_length=XI * 1e-3 * mesh_xi, smooth=50)
    return dev


def square(side_nm, mesh_xi):
    from tdgl.geometry import box
    return film_device(box(side_nm, points=int(4 * side_nm / (XI * mesh_xi)) + 4), mesh_xi=mesh_xi, name="square")


def annulus(r1_nm, r2_nm, mesh_xi):
    from tdgl.geometry import circle
    n = lambda r: max(48, int(np.ceil(2 * np.pi * r / (0.8 * XI * mesh_xi))))
    return film_device(circle(r2_nm, points=n(r2_nm)), [circle(r1_nm, points=n(r1_nm))], mesh_xi, name="ring")


def holed_strip(x_nm, a_nm, b_nm=None, y_nm=0.0):
    """The 200 nm x 0.8 um plain wire with one elliptical (b = a: circular) hole."""
    b_nm = a_nm if b_nm is None else b_nm
    return geo.Geometry(key=f"Hole_{a_nm:g}x{b_nm:g}_at_{x_nm:g}", label=f"one hole {2 * a_nm:g} x {2 * b_nm:g} nm",
                        width_nm=200.0, holes=(geo.Ellipse(x_nm, y_nm, a_nm, b_nm),))


# ==============================================================================
# Runs
# ==============================================================================
@contextmanager
def solved_free(device, tmax, eps=1.0, phase=None, seed=None, save_every=20, dt_max=0.2):
    """tdgl.solve with no terminals (no transport current): a ring, or a film relaxing on its own.  phase multiplies
    the uniform start (see fluxoid.winding_psi).  Yields the Solution and deletes its file afterwards."""
    import tdgl
    os.makedirs(OUT_DIR, exist_ok=True)
    fd, h5 = tempfile.mkstemp(prefix="tmp_val-", suffix=".h5", dir=OUT_DIR)
    os.close(fd)
    os.unlink(h5)
    sol = None
    try:
        opts = tdgl.SolverOptions(solve_time=tmax, dt_init=1e-3, dt_max=dt_max, save_every=save_every,
                                  field_units="mT", current_units="uA", output_file=h5)
        solver = tdgl.TDGLSolver(device, opts, disorder_epsilon=eps, seed_solution=seed)
        if phase is not None:
            solver.psi_init = solver.psi_init * phase
        sol = solver.solve()
        yield sol
    finally:
        path = getattr(sol, "path", None) or h5
        if path and os.path.exists(path):
            os.remove(path)


def psi_frames(sol):
    """[(t, psi)] for every saved frame."""
    import h5py
    with h5py.File(sol.path, "r") as f:
        return [(float(f["data"][g].attrs["time"]), np.asarray(f["data"][g]["psi"]))
                for g in sorted(f["data"].keys(), key=int)]


def run_frozen(device, I, eps, cfg, seed=None):
    """sim.run_simulation for a time-independent epsilon (a frozen normal disc): unlike a photon, the disc is kept in
    the sustained-switching continuation."""
    with sim.solved(device, I, eps, cfg, seed) as sol:
        trig = sim.detect(sol)
        if not cfg.sustained:
            return trig._replace(switched=bool(trig.branch))
        if not trig.switched:
            return trig
        with sim.solved(device, I, eps, replace(cfg, tmax=sim.SETTLE), seed=sol) as ext:
            t = np.asarray(ext.dynamics.time)
            later = float(sim.phase_slips(t, np.asarray(ext.dynamics.voltage()), t[-1] - sim.SUSTAIN_WIN)[-1])
    return trig._replace(switched=later >= sim.SUSTAIN_SLIPS, slips_settled=later)


def bisect(run, lo, hi, tol):
    """sim.bisect for any run(I) -> sim.Trigger.  Returns (I, flag, n_runs, trigger at I)."""
    lo0, hi0, n, trig = lo, hi, 0, sim.Trigger(False)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        n += 1
        r = run(mid)
        if r.switched:
            hi, trig = mid, r
        else:
            lo = mid
    return hi, ("at_wall" if hi >= hi0 else "below_bracket" if lo <= lo0 else ""), n, trig


def photon_threshold(geom, x_nm, y_nm, f, I_sw, cfg=CFG, params=HOTSPOT):
    """i_det = I_det / I_sw of one site, as hitmap.py computes it, for any hotspot parameters (e.g. a later t0)."""
    from snspd import hotspot as hsm
    ctx = sim.context(geom, cfg.mesh_xi, FILM, params)
    eps = hsm.build_eps(ctx["hs"], cfg.hotspot, x_nm, y_nm, f, cfg.tmax, params)
    tol = sim.scales(geom)[0]
    run = lambda I: sim.run_simulation(ctx["device"], I, eps, cfg, t0=params.t0_tau0)
    return bisect(run, cfg.hit_lo_frac * I_sw, I_sw, tol)


# ==============================================================================
# Line integrals on the mesh
# ==============================================================================
def interp_op(device, xy_um):
    """Sparse operator: mesh-site values -> linear interpolation at the points xy_um (n, 2)."""
    from matplotlib.tri import Triangulation
    from scipy.sparse import csr_matrix
    pts, tri = np.asarray(device.points), np.asarray(device.triangles)
    xy = np.asarray(xy_um, float)
    k = Triangulation(pts[:, 0], pts[:, 1], tri).get_trifinder()(xy[:, 0], xy[:, 1])
    if (k < 0).any():
        raise ValueError(f"{int((k < 0).sum())} path points are off the mesh")
    v = pts[tri[k]]
    d = (v[:, 1, 1] - v[:, 2, 1]) * (v[:, 0, 0] - v[:, 2, 0]) + (v[:, 2, 0] - v[:, 1, 0]) * (v[:, 0, 1] - v[:, 2, 1])
    l0 = ((v[:, 1, 1] - v[:, 2, 1]) * (xy[:, 0] - v[:, 2, 0]) + (v[:, 2, 0] - v[:, 1, 0]) * (xy[:, 1] - v[:, 2, 1])) / d
    l1 = ((v[:, 2, 1] - v[:, 0, 1]) * (xy[:, 0] - v[:, 2, 0]) + (v[:, 0, 0] - v[:, 2, 0]) * (xy[:, 1] - v[:, 2, 1])) / d
    w = np.stack([l0, l1, 1.0 - l0 - l1], axis=1)
    return csr_matrix((w.ravel(), (np.repeat(np.arange(len(xy)), 3), tri[k].ravel())), shape=(len(xy), len(pts)))


def polyline(corners_nm, step_nm=0.5):
    """Points along a polyline (nm in, um out), at most step_nm apart."""
    c = np.asarray(corners_nm, float)
    out = [c[:1]]
    for p, q in zip(c[:-1], c[1:]):
        n = max(1, int(np.ceil(np.hypot(*(q - p)) / step_nm)))
        out.append(p + (q - p) * np.linspace(0, 1, n + 1)[1:, None])
    return np.vstack(out) * 1e-3


def phase_change(op, psi):
    """(change of arg psi along the path [rad], smallest |psi| on it)."""
    z = op @ psi
    return float(np.angle(z[1:] / z[:-1]).sum()), float(np.abs(z).min())


def current_through(sol, xy_um):
    """Total current [uA] through a path (n, 2) [um], positive to the left of the walking direction, from pyTDGL's
    current density (site values rebuilt from the edge currents, then interpolated).  pyTDGL's own
    Solution.current_through_path fails with numpy 2 (2-vector np.cross)."""
    K = np.asarray(sol.interp_current_density(xy_um, units="uA/um", with_units=False))
    d = np.diff(xy_um, axis=0)
    Km = 0.5 * (K[1:] + K[:-1])
    return float(np.sum(Km[:, 0] * -d[:, 1] + Km[:, 1] * d[:, 0]))


def supercurrent_through(device, psi, xy_um):
    """Supercurrent [uA] through a path, from psi alone: K_s = (K0/4) Im(psi* grad psi), grad in 1/xi, with psi and
    its gradient linear on each triangle (no boundary-site reconstruction).  Steady states only (no normal current)."""
    from matplotlib.tri import Triangulation
    pts, tri = np.asarray(device.points) * 1e3, np.asarray(device.triangles)
    xy = np.asarray(xy_um, float) * 1e3
    mid, d = 0.5 * (xy[1:] + xy[:-1]), np.diff(xy, axis=0)
    k = Triangulation(pts[:, 0], pts[:, 1], tri).get_trifinder()(mid[:, 0], mid[:, 1])
    if (k < 0).any():
        raise ValueError("path leaves the mesh")
    v = pts[tri[k]]
    T = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]], axis=1)               # (n, 2, 2): rows = edge vectors
    z = psi[tri[k]]
    grad = np.linalg.solve(T, np.stack([z[:, 1] - z[:, 0], z[:, 2] - z[:, 0]], axis=1)[..., None])[..., 0]
    lam = np.linalg.solve(np.transpose(T, (0, 2, 1)), (mid - v[:, 0])[..., None])[..., 0]
    zm = z[:, 0] + lam[:, 0] * (z[:, 1] - z[:, 0]) + lam[:, 1] * (z[:, 2] - z[:, 0])
    K = FILM.K0_uA_per_nm / 4 * XI * np.imag(np.conj(zm)[:, None] * grad)         # uA/nm
    return float(np.sum(K[:, 0] * -d[:, 1] + K[:, 1] * d[:, 0]))


def two_channel_isw(cap_l, cap_r, I_c, adds_left=True):
    """Switching current of two parallel channels (capacities cap_l, cap_r) around a hole holding one quantum, each
    with the GL current-phase relation I = cap g(q), g(q) = (3 sqrt 3 / 2)(1 - q^2) q, max 1 at q = 1/sqrt 3, and
    the fluxoid fixing the difference of their phase gradients, q_l - q_r = Delta.  Delta is set by the circulating
    current I_c at zero bias (I_l = -I_r = I_c); I_sw is the largest total current over q.  Equal channels give
    I_sw = I_sw(0) (1 - 3 delta^2)^(3/2), delta = Delta / 2: a drop quadratic in I_c, not the linear 2 I_c."""
    from scipy.optimize import brentq, minimize_scalar
    g = lambda q: 1.5 * np.sqrt(3) * (1 - q ** 2) * q
    qs = 1 / np.sqrt(3)

    def zero_bias(D):            # q_r with cap_l g(q_r + D) + cap_r g(q_r) = 0, and the circulation it carries
        q_r = brentq(lambda q: cap_l * g(q + D) + cap_r * g(q), -qs, 0.0)
        return q_r, cap_l * g(q_r + D)

    D = brentq(lambda D: zero_bias(D)[1] - abs(I_c), 1e-9, qs)
    D = D if adds_left else -D
    best = minimize_scalar(lambda q: -(cap_l * g(q + D) + cap_r * g(q)), bounds=(-qs, qs), method="bounded")
    return float(-best.fun), float(D)
