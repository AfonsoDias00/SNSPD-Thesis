"""
The device and its references.  All coordinates in nm; the strip is centred on the origin, current flows along y, x
runs across the width.

  Double_400   a 200 nm strip with two aligned columns of slots, 40 nm everywhere: 40 | slot 40 | 40 | slot 40 | 40,
               slots 400 nm long with semicircular ends, 40 nm bridges between consecutive slots.  Three 40 nm
               channels in parallel; the narrowest free cross-section is 120 nm.
  Plain_Wire   the plain 200 nm strip: same footprint.
  Plain_120    plain strip of the free width: same depairing current.
  Plain_40     plain strip of the channel width.

A Geometry carries what the simulation and the analysis need: the holes, the photon-hit grids ("full", "critical"),
the probe sites, and the symmetry that completes a threshold map over the unit cell ("strip": depends on |x| only;
"column": x-mirror and y-mirror).

The other layouts of the earlier phases (single and staggered columns, notches, zipper, ellipses) are on the branch
slot-geometries.  double_column() keeps the slot length as a parameter: longer slots were the one lever that helped.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

C_PLAIN, C_HOLES, C_NARROW = "#2a78d6", "#eb6834", "#52514e"


@dataclass(frozen=True)
class Slot:
    """Straight-sided hole with semicircular ends, long along y."""
    cx_nm: float
    cy_nm: float
    a_nm: float        # half-width along x, also the end radius
    b_nm: float        # half-length along y, ends included

    def contains(self, x_nm, y_nm, margin_nm=0.0):
        dx = np.asarray(x_nm, float) - self.cx_nm
        dy = np.maximum(np.abs(np.asarray(y_nm, float) - self.cy_nm) - (self.b_nm - self.a_nm), 0.0)
        return dx ** 2 + dy ** 2 < (self.a_nm + margin_nm) ** 2

    def perimeter(self):
        return 4.0 * (self.b_nm - self.a_nm) + 2.0 * np.pi * self.a_nm

    def outline(self, n):
        """n points equally spaced along the boundary, counter-clockwise from (cx + a, cy)."""
        a, s = self.a_nm, self.b_nm - self.a_nm
        arc, per = np.pi * a, self.perimeter()
        pts = []
        for d in np.linspace(0.0, per, n, endpoint=False):
            if d < 2 * s:
                pts.append((a, -s + d))
            elif d < 2 * s + arc:
                th = (d - 2 * s) / a
                pts.append((a * np.cos(th), s + a * np.sin(th)))
            elif d < 4 * s + arc:
                pts.append((-a, s - (d - 2 * s - arc)))
            else:
                th = np.pi + (d - 4 * s - arc) / a
                pts.append((a * np.cos(th), -s + a * np.sin(th)))
        return np.array(pts) + (self.cx_nm, self.cy_nm)


Site = tuple[float, float]


@dataclass(frozen=True)
class Geometry:
    key: str                                   # name used in CSV rows and file names
    label: str
    width_nm: float
    length_nm: float = 800.0
    holes: tuple = ()                          # objects with contains(x, y, margin) and outline(n), e.g. Slot
    probes: tuple[Site, ...] = ()
    grid_full: tuple[Site, ...] = ()
    grid_critical: tuple[Site, ...] = ()
    symmetry: str = "strip"                    # "strip" or "column"
    cell_half_y_nm: float = 61.0               # unit cell for threshold maps: |y| <= this
    hole_points: int = 38                      # boundary points per hole polygon
    probe_y_nm: float | None = None            # voltage probes at (0, +/-this); None = length / 3
    channels: tuple[tuple[float, float], ...] = ()   # x-ranges of the parallel channels at y = 0
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

    def squares(self) -> float:
        """Normal-state resistance of the strip between the voltage probes, in squares (series of cross-sections)."""
        py = self.length_nm / 3 if self.probe_y_nm is None else self.probe_y_nm
        ys = np.linspace(-py, py, 4001)
        xs = np.linspace(-self.width_nm / 2, self.width_nm / 2, 801)
        dx = xs[1] - xs[0]
        w = np.array([(~self.inside_hole(xs, y)).sum() * dx for y in ys])
        return float(np.sum((ys[1] - ys[0]) / w[:-1]))

    def hit_grid(self, grid: str) -> list[Site]:
        """Photon-hit sites (x, y) [nm]: the chosen grid minus sites within 2 nm of a hole."""
        pts = self.grid_critical if grid == "critical" else self.probes if grid == "probes" else self.grid_full
        keep = {(float(x), float(y)) for x, y in pts if not self.inside_hole(x, y, 2.0)}
        return sorted(keep, key=lambda p: (p[1], p[0]))


# ==============================================================================
# Plain strips
# ==============================================================================
def strip_sites(width_nm: float) -> list[float]:
    """x positions across the left half: 2, 8, 14, 21 nm in from the edge (those that fit), the quarter points and
    the centre; a site within 3 nm of one already chosen is skipped."""
    h = width_nm / 2.0
    xs = [-(h - o) for o in (2.0, 8.0, 14.0, 21.0) if o <= h - 3.0]
    for x in (-float(round(h / 2.0)), -float(round(h / 4.0)), 0.0):
        if all(abs(x - s) >= 3.0 for s in xs):
            xs.append(x)
    return sorted(xs)


def plain_strip(width_nm: float, key: str, label: str | None = None, xs_full=None, xs_critical=None,
                color: str = C_NARROW) -> Geometry:
    """Hole-free strip, 0.8 um long; photons are simulated on the left half and mirrored."""
    xs_full = list(xs_full) if xs_full is not None else strip_sites(width_nm)
    xs_crit = list(xs_critical) if xs_critical is not None else [xs_full[0], xs_full[len(xs_full) // 2], 0.0]
    crit = tuple((float(x), 0.0) for x in xs_crit)
    return Geometry(key=key, label=label or f"Plain {width_nm:g} nm", width_nm=float(width_nm), probes=crit,
                    grid_full=tuple((float(x), 0.0) for x in xs_full), grid_critical=crit, symmetry="strip",
                    channels=((-width_nm / 2.0, width_nm / 2.0),), color=color)


PLAIN_WIRE = plain_strip(200.0, "Plain_Wire", "Plain wire, 200 nm", xs_full=(-96, -88, -74, -58, -42, -28, -14, 0),
                         xs_critical=(-96, -58, 0), color=C_PLAIN)
PLAIN_120 = plain_strip(120.0, "Plain_120", "Plain 120 nm (the free width)")
PLAIN_40 = plain_strip(40.0, "Plain_40", "Plain 40 nm (one channel)")


# ==============================================================================
# The double slot column
# ==============================================================================
WIDTH_NM = 200.0
FEATURE_NM = 40.0         # slot width, channel width and bridge length: the fabrication limit
LEAD_NM = 150.0           # plain strip between the outermost slot end and each terminal
HOLE_POINT_SPACING_NM = 5.0


def double_column(slot_nm: float = 400.0, neighbours: int = 1, length_nm: float | None = None) -> Geometry:
    """Two aligned columns of 40 nm wide slots at x = +/-40 nm in a 200 nm strip; slots slot_nm long (semicircular
    ends included), 40 nm bridges, period slot_nm + 40.

    A slot pair is centred on y = 0 with `neighbours` pairs on each side.  The photon unit cell is |x| <= 100,
    |y| <= half a period; x-mirror and y-mirror complete it.  length_nm None: the shortest device, with LEAD_NM of
    plain strip beyond the outermost slot end on each side (1.58 um for 400 nm slots).  On the branch slot-geometries
    this short device reproduced a 2 um strip within 0.011 I_dep at 1.4 xi."""
    period, b, a = slot_nm + FEATURE_NM, slot_nm / 2.0, FEATURE_NM / 2.0
    holes = tuple(Slot(xc, k * period, a, b) for xc in (-FEATURE_NM, FEATURE_NM) for k in range(-neighbours, neighbours + 1))
    half = neighbours * period + b
    if length_nm is None:
        length_nm = 2.0 * (half + LEAD_NM)
    if length_nm / 2.0 < half + 50.0:
        raise ValueError(f"length {length_nm:g} nm leaves no lead beyond the slots (half-extent {half:g} nm)")
    ys = []
    for y in sorted((0.0, b / 2.0, b - a, b, b + 10.0, b + 20.0)):     # rows closer than 5 nm merge
        if not ys or y - ys[-1] >= 5.0:
            ys.append(y)
    channel = [(x, y) for x in (-96.0, -88.0, -80.0, -72.0, -64.0, -16.0, -8.0, 0.0) for y in ys]
    bridge = [(x, y) for x in (-56.0, -48.0, -40.0, -32.0, -24.0) for y in (b + 4.0, b + 12.0, b + 20.0)]
    # edge channel (edge, middle, slot side), centre channel (slot side, middle), bridge row (slot end, mid-bridge)
    critical = [(-96.0, 0.0), (-80.0, 0.0), (-64.0, 0.0), (-16.0, 0.0), (0.0, 0.0), (-40.0, b + 4.0),
                (-80.0, b + 20.0), (-40.0, b + 20.0), (0.0, b + 20.0)]
    probes = ((-80.0, 0.0), (0.0, 0.0), (-40.0, b + 20.0))            # edge channel, centre channel, bridge centre
    tag = "" if abs(length_nm - 2.0 * (half + LEAD_NM)) < 1e-6 else f"_{length_nm / 1000:g}um"
    return Geometry(
        key=f"Double_{slot_nm:g}{tag}", label=f"Double column, {slot_nm:g} nm slots", width_nm=WIDTH_NM,
        length_nm=float(length_nm), holes=holes, probes=probes, grid_full=tuple(channel + bridge),
        grid_critical=tuple(critical), symmetry="column", cell_half_y_nm=period / 2.0,
        hole_points=int(np.ceil(holes[0].perimeter() / HOLE_POINT_SPACING_NM)),
        probe_y_nm=length_nm / 2.0 - LEAD_NM / 2.0, channels=((-100.0, -60.0), (-20.0, 20.0), (60.0, 100.0)),
        color=C_HOLES)


DOUBLE_400 = double_column(400.0)

GEOMETRIES = {g.key: g for g in (DOUBLE_400, PLAIN_WIRE, PLAIN_120, PLAIN_40)}


def get(key: str) -> Geometry:
    """Registered geometry by key; "Double_<L>" builds the double column with L nm slots."""
    if key in GEOMETRIES:
        return GEOMETRIES[key]
    if key.startswith("Double_"):
        return double_column(float(key.split("_")[1]))
    raise KeyError(f"unknown geometry {key!r}; registered: {', '.join(GEOMETRIES)}")
