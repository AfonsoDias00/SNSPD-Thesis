"""
Strip and hole geometries, photon-hit grids and unit cells.  All coordinates in nm; the strip is centred on the
origin, current flows along y, x runs across the width.

A Geometry carries everything the simulation and the analysis need to know about one layout:
  * holes (ellipses), strip width and length
  * the monitored cross-sections whose supercurrent collapse counts as a detection
  * the photon-hit grids ("full" and "critical") and the probe sites used for the photon-energy sweep
  * the unit cell and symmetry used to turn the sampled thresholds into a map:
      "strip"  -- threshold depends on |x| only (plain strips)
      "zipper" -- y-mirror plus the glide (x, y) -> (-x, half_period - y) of two staggered hole columns

To study a new internal-hole layout, add a builder here and register it in GEOMETRIES.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

C_PLAIN, C_ZIP, C_NARROW = "#2a78d6", "#eb6834", "#52514e"


@dataclass(frozen=True)
class Ellipse:
    cx_nm: float
    cy_nm: float
    a_nm: float        # semi-axis along x
    b_nm: float        # semi-axis along y

    def contains(self, x_nm, y_nm, margin_nm=0.0):
        return ((np.asarray(x_nm, float) - self.cx_nm) / (self.a_nm + margin_nm)) ** 2 + \
               ((np.asarray(y_nm, float) - self.cy_nm) / (self.b_nm + margin_nm)) ** 2 < 1.0


Site = tuple[float, float]


@dataclass(frozen=True)
class Geometry:
    key: str                                   # name used in CSV rows and file names
    label: str
    width_nm: float
    length_nm: float = 800.0
    holes: tuple[Ellipse, ...] = ()
    slices: tuple[tuple[float, float, float], ...] = ()   # monitored cross-sections (x_min, x_max, y) [nm]
    probes: tuple[Site, ...] = ()
    grid_full: tuple[Site, ...] = ()
    grid_critical: tuple[Site, ...] = ()
    symmetry: str = "strip"                    # "strip" or "zipper"
    cell_half_y_nm: float = 61.0               # unit cell for threshold maps: |y| <= this
    crowd_half_y_nm: float = 122.0             # current crowding is evaluated in |y| <= this, away from the terminals
    hole_points: int = 38                      # boundary points per hole polygon
    color: str = field(default=C_NARROW, compare=False)

    def inside_hole(self, x_nm, y_nm, margin_nm=0.0):
        x, y = np.asarray(x_nm, float), np.asarray(y_nm, float)
        out = np.zeros(np.broadcast(x, y).shape, bool)
        for h in self.holes:
            out |= h.contains(x, y, margin_nm)
        return out

    def free_width_min_nm(self) -> float:
        """Narrowest free cross-section: strip width minus the hole chords, minimised over one unit cell."""
        if not self.holes:
            return self.width_nm
        ys = np.linspace(-self.cell_half_y_nm, self.cell_half_y_nm, 241)
        xs = np.linspace(-self.width_nm / 2, self.width_nm / 2, 2001)
        dx = xs[1] - xs[0]
        return float(min((~self.inside_hole(xs, y)).sum() * dx for y in ys))

    def hit_grid(self, grid: str) -> list[Site]:
        """Photon-hit sites (x, y) [nm]: the chosen grid minus sites within 2 nm of a hole, plus the probe sites."""
        pts = self.grid_critical if grid == "critical" else self.grid_full
        keep = {(float(x), float(y)) for x, y in pts if not self.inside_hole(x, y, 2.0)}
        keep |= {(float(x), float(y)) for x, y in self.probes}
        return sorted(keep, key=lambda p: (p[1], p[0]))


# ==============================================================================
# Plain strips
# ==============================================================================
def strip_sites(width_nm: float) -> list[float]:
    """x positions across the left half: 2, 8, 14, 21 nm in from the edge (those that fit), the quarter points and
    the centre; a site within 3 nm of one already chosen is skipped.  56 nm gives -26, -20, -14, -7, 0."""
    h = width_nm / 2.0
    xs = [-(h - o) for o in (2.0, 8.0, 14.0, 21.0) if o <= h - 3.0]
    for x in (-float(round(h / 2.0)), -float(round(h / 4.0)), 0.0):
        if all(abs(x - s) >= 3.0 for s in xs):
            xs.append(x)
    return sorted(xs)


def strip_probes(width_nm: float) -> list[Site]:
    """Edge site, the site nearest the half-way point to the edge, and the centre."""
    xs = strip_sites(width_nm)
    return [(float(x), 0.0) for x in (xs[0], min(xs, key=lambda x: abs(x + width_nm / 4.0)), 0.0)]


def plain_strip(width_nm: float, key: str = "Narrow_Wire", label: str | None = None, xs_full=None, xs_critical=None,
                probes=None, color: str = C_NARROW) -> Geometry:
    """Hole-free strip.  The monitored cross-section is the right half at y = 0 (the photon hits the left half)."""
    h = width_nm / 2.0
    xs_full = list(xs_full) if xs_full is not None else strip_sites(width_nm)
    xs_crit = list(xs_critical) if xs_critical is not None else xs_full
    return Geometry(
        key=key, label=label or f"Plain {width_nm:g} nm", width_nm=float(width_nm),
        slices=((2.0, h - 2.0, 0.0),),
        probes=tuple(probes) if probes is not None else tuple(strip_probes(width_nm)),
        grid_full=tuple((float(x), 0.0) for x in xs_full),
        grid_critical=tuple((float(x), 0.0) for x in xs_crit),
        symmetry="strip", color=color)


PLAIN_WIRE = plain_strip(
    200.0, key="Plain_Wire", label="Plain wire",
    xs_full=(-96, -88, -74, -58, -42, -28, -14, 0), xs_critical=(-96, -58, 0),
    probes=((-96.0, 0.0), (-58.0, 0.0), (0.0, 0.0)), color=C_PLAIN)


# ==============================================================================
# Internal holes
# ==============================================================================
def long_zipper() -> Geometry:
    """Two staggered columns of 16 x 110 nm ellipses at x = +/-36 nm, period 122 nm, 12 nm bridges between the holes
    of a column, three 56 nm channels.  Narrowest free width 173 nm.  The right channel at y = 61 nm (beside the
    left column's bridge) is the monitored cross-section."""
    a, b = 8.0, 55.0
    holes = [Ellipse(-36.0, y, a, b) for y in (-122.0, 0.0, 122.0)] + \
            [Ellipse(36.0, y, a, b) for y in (-183.0, -61.0, 61.0, 183.0)]
    base = [(x, y) for x in (-96, -88, -68, -50, -18, 0) for y in (0, 15, 30, 45, 61)]
    bridge = [(x, y) for x in (-42, -36, -30) for y in (46, 50, 54, 58, 61)]
    critical = [(-96, 0), (-68, 0), (0, 0), (-50, 30), (-18, 30), (-36, 58), (-36, 61), (-42, 58), (-30, 58)]
    return Geometry(
        key="Long_Zipper", label="Long zipper", width_nm=200.0, holes=tuple(holes),
        slices=((46.0, 98.0, 61.0),),
        probes=((-96.0, 0.0), (-68.0, 0.0), (0.0, 0.0), (-18.0, 30.0), (-36.0, 61.0)),
        grid_full=tuple((float(x), float(y)) for x, y in base + bridge),
        grid_critical=tuple((float(x), float(y)) for x, y in critical),
        symmetry="zipper", cell_half_y_nm=61.0, color=C_ZIP)


LONG_ZIPPER = long_zipper()

GEOMETRIES = {g.key: g for g in (PLAIN_WIRE, LONG_ZIPPER)}


def get(key: str, width_nm: float | None = None) -> Geometry:
    """Registered geometry by key; "Narrow_Wire" needs width_nm."""
    if key == "Narrow_Wire":
        if width_nm is None:
            raise ValueError("Narrow_Wire needs a width")
        return plain_strip(width_nm)
    return GEOMETRIES[key]
