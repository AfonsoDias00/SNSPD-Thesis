"""
Fluxoid of every hole: how many flux quanta a hole holds.

With screening off and no applied field the vector potential is zero, so the fluxoid of a hole is the winding number
of arg(psi) along a loop around it -- the same quantity pyTDGL's Solution.hole_fluxoid integrates from the
interpolated supercurrent, but exact (an integer) and fast enough to evaluate on every saved frame.  The loop is
pyTDGL's fluxoid polygon: the hole outline pushed out by half the distance to the nearest other boundary (20 nm, the
middle of the 40 nm channels and bridges, in the double column), sampled at n_points and interpolated linearly from
the mesh.  The
winding is only defined while |psi| stays finite on the loop; frames where a vortex core or a normal belt crosses
it are flagged by the loop minimum of |psi|.
"""
from __future__ import annotations

import numpy as np

PSI_LOOP_MIN = 0.05      # frames with min |psi| on a loop below this have no well-defined winding


class FluxoidMeter:
    def __init__(self, device, n_points: int = 400):
        import tdgl
        from matplotlib.tri import Triangulation
        from scipy.sparse import csr_matrix
        pts, tri = np.asarray(device.points), np.asarray(device.triangles)
        finder = Triangulation(pts[:, 0], pts[:, 1], tri).get_trifinder()
        loops = tdgl.make_fluxoid_polygons(device, interp_points=n_points)
        self.names = [h.name for h in device.holes]
        self.ops = []
        for name in self.names:
            xy = np.asarray(loops[name])
            if np.allclose(xy[0], xy[-1]):
                xy = xy[:-1]
            if _signed_area(xy) < 0:            # counter-clockwise: positive winding = positive fluxoid
                xy = xy[::-1]
            k = finder(xy[:, 0], xy[:, 1])
            if (k < 0).any():
                raise ValueError(f"fluxoid loop of {name} leaves the mesh")
            v = pts[tri[k]]                     # (n, 3, 2) triangle corners
            d = (v[:, 1, 1] - v[:, 2, 1]) * (v[:, 0, 0] - v[:, 2, 0]) + (v[:, 2, 0] - v[:, 1, 0]) * (v[:, 0, 1] - v[:, 2, 1])
            l0 = ((v[:, 1, 1] - v[:, 2, 1]) * (xy[:, 0] - v[:, 2, 0]) + (v[:, 2, 0] - v[:, 1, 0]) * (xy[:, 1] - v[:, 2, 1])) / d
            l1 = ((v[:, 2, 1] - v[:, 0, 1]) * (xy[:, 0] - v[:, 2, 0]) + (v[:, 0, 0] - v[:, 2, 0]) * (xy[:, 1] - v[:, 2, 1])) / d
            w = np.stack([l0, l1, 1.0 - l0 - l1], axis=1)
            rows = np.repeat(np.arange(len(xy)), 3)
            self.ops.append(csr_matrix((w.ravel(), (rows, tri[k].ravel())), shape=(len(xy), len(pts))))

    def __call__(self, psi):
        """(fluxoid per hole as ints, smallest |psi| on any loop)."""
        n, floor = [], np.inf
        for op in self.ops:
            z = op @ psi
            n.append(int(round(np.angle(np.roll(z, -1) / z).sum() / (2 * np.pi))))
            floor = min(floor, float(np.abs(z).min()))
        return n, floor


def winding_psi(device, n):
    """exp(i theta) on the mesh sites, theta winding n[k] times 2 pi around hole k (positive counter-clockwise, as
    FluxoidMeter counts): a starting order parameter that holds fluxoid n[k] in each hole, with no vortex in the
    film.  The phase is centred on each hole's outline points; the solver relaxes the current distribution."""
    xy = np.asarray(device.points)
    theta = np.zeros(len(xy))
    for k, hole in zip(n, device.holes):
        cx, cy = np.asarray(hole.points).mean(axis=0)
        theta += k * np.arctan2(xy[:, 1] - cy, xy[:, 0] - cx)
    return np.exp(1j * theta)


def _signed_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def frames(sol, meter: FluxoidMeter):
    """[(t, fluxoid per hole, min |psi| on the loops, max |V| since the previous frame)] for every saved frame."""
    import h5py
    t = np.asarray(sol.dynamics.time)
    v = np.abs(np.asarray(sol.dynamics.voltage()))
    out, t_prev = [], -np.inf
    with h5py.File(sol.path, "r") as f:
        for g in sorted(f["data"].keys(), key=int):
            tm = float(f["data"][g].attrs["time"])
            n, floor = meter(np.asarray(f["data"][g]["psi"]))
            win = (t > t_prev) & (t <= tm)
            out.append((tm, n, floor, float(v[win].max()) if win.any() else 0.0))
            t_prev = tm
    return out
