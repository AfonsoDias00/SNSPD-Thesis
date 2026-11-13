"""
From sampled thresholds to maps, efficiency curves and quantiles.

The hard-step efficiency eta(i_b) is the area fraction of the unit cell where i_det(r) < i_b: every photon absorbed
where the local threshold is below the bias counts, none elsewhere.  Sampled thresholds are interpolated over the
unit cell using the geometry's symmetry.  Biases are in I_dep of the narrowest free cross-section.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import griddata

from .geometry import Geometry


def photon_rows(rows, geom_key):
    """(x, y, i_det) of every photon row of this geometry; 'at_wall' (no count up to the bracket top) counts as 1."""
    out = []
    for r in rows:
        if r["Geometry"] == geom_key and r["x_nm"] not in ("", None):
            v = 1.0 if r["flag"] == "at_wall" else float(r["i_det"])
            out.append((float(r["x_nm"]), float(r["y_nm"]), v))
    return out


def _symmetrize_column(samples):
    """x-mirror and y-mirror of a single centred hole column."""
    pts = {(sx * x, sy * y): v for x, y, v in samples for sx in (1, -1) for sy in (1, -1)}
    return np.array(list(pts)), np.array(list(pts.values()))


def map_grid(geom: Geometry):
    xs = np.linspace(-geom.width_nm / 2, geom.width_nm / 2, 201)
    ys = np.linspace(-geom.cell_half_y_nm, geom.cell_half_y_nm, int(round(2 * geom.cell_half_y_nm)) + 1)
    return np.meshgrid(xs, ys)


def unit_map(geom: Geometry, samples):
    """Threshold map i_det(x, y) over one unit cell (NaN inside holes)."""
    X, Y = map_grid(geom)
    if geom.symmetry == "strip":
        xs = np.array([abs(s[0]) for s in samples])
        vs = np.array([s[2] for s in samples])
        o = np.argsort(xs)
        return X, Y, np.interp(np.abs(X), xs[o], vs[o])
    if geom.symmetry == "column":
        pts, vals = _symmetrize_column(samples)
        Z = griddata(pts, vals, (X, Y), method="linear")
        Z = np.where(np.isnan(Z), griddata(pts, vals, (X, Y), method="nearest"), Z)
        Z[geom.inside_hole(X, Y)] = np.nan
        return X, Y, Z
    raise ValueError(f"no map rule for symmetry {geom.symmetry!r}")


def pixel_values(Z):
    return np.sort(Z[np.isfinite(Z)].ravel())


def eta_curve(vals, ib):
    return np.searchsorted(vals, ib, side="right") / len(vals) * 100.0


def metrics(vals):
    q = lambda p: float(np.quantile(vals, p))
    return dict(i_10=q(0.10), i_50=q(0.50), i_90=q(0.90), i_95=q(0.95), i_99=q(0.99), i_max=float(vals[-1]),
                width_10_90=q(0.90) - q(0.10))


def fmt(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"
