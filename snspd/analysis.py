"""
From sampled thresholds to maps, efficiency curves and quantiles.

The hard-step efficiency eta(i_b) is the area fraction of the unit cell where i_det(r) < i_b, the deterministic
definition used for intrinsic detection efficiency in hotspot/TDGL models (Zotova & Vodolazov, SUST 27, 125001
(2014)).  Sampled thresholds are interpolated over the unit cell using the geometry's symmetry.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import griddata

from .geometry import Geometry


def photon_rows(rows, geom_key, f=1.0):
    """(x, y, i_det) of every row of this geometry at energy factor f; 'at_wall' (no detection) counts as 1."""
    out = []
    for r in rows:
        if r["Geometry"] == geom_key and abs(float(r["energy_f"]) - f) < 1e-6:
            v = 1.0 if r["flag"] == "at_wall" else float(r["i_det"])
            out.append((float(r["x_nm"]), float(r["y_nm"]), v))
    return out


def _symmetrize_zipper(samples, half_period):
    """y-mirror and the glide (x, y) -> (-x, half_period - y) of two staggered hole columns."""
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
            yg = half_period - y
            for yy in (yg, -yg):
                add(-x, yy, v)
    return np.array(pts), np.array(vals)


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
    if geom.symmetry == "zipper":
        pts, vals = _symmetrize_zipper(samples, geom.cell_half_y_nm)
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


def probe_table(rows, geom: Geometry):
    """{probe site: {f: (i_det, flag)}} for the geometry's probe sites (matched to whole nm)."""
    probes = {(int(x), int(y)) for x, y in geom.probes}
    tab = {}
    for r in rows:
        site = (int(float(r["x_nm"])), int(float(r["y_nm"])))
        if r["Geometry"] == geom.key and site in probes:
            v = 1.0 if r["flag"] == "at_wall" else float(r["i_det"])
            tab.setdefault(site, {})[round(float(r["energy_f"]), 3)] = (v, r["flag"])
    return tab


def selectivity(tab):
    """Median over probe sites of i_det(f) - i_det(1) for every weak factor f present (larger = weak seeds need
    more bias)."""
    fs = sorted({f for site in tab.values() for f in site if f < 0.999})
    out = {}
    for f in fs:
        diffs = [site[f][0] - site[1.0][0] for site in tab.values() if f in site and 1.0 in site]
        if diffs:
            out[f] = float(np.median(diffs))
    return out


def fmt(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"
