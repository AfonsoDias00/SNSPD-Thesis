"""
Strip and hole geometries, photon-hit grids and unit cells.  All coordinates in nm; the strip is centred on the
origin, current flows along y, x runs across the width.

A Geometry carries everything the simulation and the analysis need to know about one layout:
  * holes (ellipses or slots), strip width and length
  * the photon-hit grids ("full" and "critical") and the probe sites used for the photon-energy sweep
  * the unit cell and symmetry used to turn the sampled thresholds into a map:
      "strip"  -- threshold depends on |x| only (plain strips)
      "zipper" -- y-mirror plus the glide (x, y) -> (-x, half_period - y) of two staggered hole columns
      "column" -- x-mirror and y-mirror of a single centred hole column

To study a new internal-hole layout, add a builder here and register it in GEOMETRIES.  Registered now: the plain
reference strips only.  slot_column() is the single-column builder of the first phase (30 nm slots, 40 nm bridges),
kept unregistered for validation/ until the new fabrication-friendly layouts replace it.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

C_PLAIN, C_HOLES, C_NARROW = "#2a78d6", "#eb6834", "#52514e"


@dataclass(frozen=True)
class Ellipse:
    cx_nm: float
    cy_nm: float
    a_nm: float        # semi-axis along x
    b_nm: float        # semi-axis along y

    def contains(self, x_nm, y_nm, margin_nm=0.0):
        return ((np.asarray(x_nm, float) - self.cx_nm) / (self.a_nm + margin_nm)) ** 2 + \
               ((np.asarray(y_nm, float) - self.cy_nm) / (self.b_nm + margin_nm)) ** 2 < 1.0

    def outline(self, n):
        th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        return np.column_stack([self.cx_nm + self.a_nm * np.cos(th), self.cy_nm + self.b_nm * np.sin(th)])


@dataclass(frozen=True)
class Slot:
    """Straight-sided hole with semicircular ends, long along y.  A long ellipse would need a ~1 nm tip radius, which
    cannot be fabricated; here the tip radius is a_nm."""
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
    holes: tuple[Ellipse | Slot, ...] = ()
    probes: tuple[Site, ...] = ()
    grid_full: tuple[Site, ...] = ()
    grid_critical: tuple[Site, ...] = ()
    symmetry: str = "strip"                    # "strip", "zipper" or "column"
    cell_half_y_nm: float = 61.0               # unit cell for threshold maps: |y| <= this
    crowd_half_y_nm: float = 122.0             # current crowding is evaluated in |y| <= this, away from the terminals
    hole_points: int = 38                      # boundary points per hole polygon
    probe_y_nm: float | None = None            # probes at (0, +/-this); None = length / 3
    notches: tuple[Ellipse, ...] = ()          # cuts into the film edge (circles centred on x = +/-width/2)
    color: str = field(default=C_NARROW, compare=False)

    def inside_hole(self, x_nm, y_nm, margin_nm=0.0):
        x, y = np.asarray(x_nm, float), np.asarray(y_nm, float)
        out = np.zeros(np.broadcast(x, y).shape, bool)
        for h in self.holes + self.notches:
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
    """Hole-free strip; photons are simulated on the left half and mirrored."""
    xs_full = list(xs_full) if xs_full is not None else strip_sites(width_nm)
    xs_crit = list(xs_critical) if xs_critical is not None else xs_full
    return Geometry(
        key=key, label=label or f"Plain {width_nm:g} nm", width_nm=float(width_nm),
        probes=tuple(probes) if probes is not None else tuple(strip_probes(width_nm)),
        grid_full=tuple((float(x), 0.0) for x in xs_full),
        grid_critical=tuple((float(x), 0.0) for x in xs_crit),
        symmetry="strip", color=color)


PLAIN_WIRE = plain_strip(
    200.0, key="Plain_Wire", label="Plain wire",
    xs_full=(-96, -88, -74, -58, -42, -28, -14, 0), xs_critical=(-96, -58, 0),
    probes=((-96.0, 0.0), (-58.0, 0.0), (0.0, 0.0)), color=C_PLAIN)


# Plain reference strips: 0.8 um long (the first phase) and 2 um long, the length of the slot devices
LONG_LENGTH_NM = 2000.0
LONG_PROBE_Y_NM = 750.0
PLAIN_WIRE_2UM = replace(
    PLAIN_WIRE, key="Plain_Wire_2um", label="Plain wire, 2 um", length_nm=LONG_LENGTH_NM, probe_y_nm=LONG_PROBE_Y_NM)


# ==============================================================================
# Internal holes
# ==============================================================================
SLOT_WIDTH_NM = 30.0
SLOT_BRIDGE_NM = 40.0
HOLE_POINT_SPACING_NM = 5.0


def slot_column(slot_nm: float, n_holes: int) -> Geometry:
    """One centred column of 30 nm wide slots, slot_nm long, 40 nm bridges, two 85 nm channels, narrowest free width
    170 nm.  Hole number (n_holes - 1) / 2 is centred on y = 0; the unit cell runs from its centre to a bridge centre."""
    period = slot_nm + SLOT_BRIDGE_NM
    a, b = SLOT_WIDTH_NM / 2.0, slot_nm / 2.0
    holes = tuple(Slot(0.0, (i - (n_holes - 1) / 2.0) * period, a, b) for i in range(n_holes))
    ys = (0.0, b / 2.0, b - 10.0, b, b + 10.0, b + 20.0)
    channel = [(x, y) for x in (-96.0, -88.0, -74.0, -58.0, -42.0, -28.0, -19.0) for y in ys]
    bridge = [(x, y) for x in (-12.0, -6.0, 0.0) for y in (b + 4.0, b + 12.0, b + 20.0)]
    critical = [(-96.0, 0.0), (-58.0, 0.0), (-19.0, 0.0), (-58.0, b + 20.0), (-28.0, b), (0.0, b + 4.0), (0.0, b + 20.0)]
    return Geometry(
        key=f"Slot_{slot_nm:g}", label=f"Slot column, {slot_nm:g} nm", width_nm=200.0, length_nm=LONG_LENGTH_NM,
        holes=holes, probes=((-96.0, 0.0), (-58.0, 0.0), (-19.0, 0.0), (-58.0, b + 20.0), (0.0, b + 20.0)),
        grid_full=tuple(channel + bridge), grid_critical=tuple(critical), symmetry="column",
        cell_half_y_nm=period / 2.0, crowd_half_y_nm=period / 2.0,
        hole_points=int(np.ceil(holes[0].perimeter() / HOLE_POINT_SPACING_NM)), probe_y_nm=LONG_PROBE_Y_NM, color=C_HOLES)


# ==============================================================================
# Second phase (from 6 Nov 2026): fabrication-friendly slot layouts in a 200 nm strip
# ==============================================================================
WIDTH_NM = 200.0
BRIDGE_NM = 40.0          # bridge between consecutive slots of a column, along y
LEAD_NM = 150.0           # plain strip between the outermost slot end and each terminal (short devices)
LAYOUTS = ("single", "double", "stagger")


def slotted_half_nm(geom: Geometry) -> float:
    """Half-extent along y of the slotted section: the outermost slot end."""
    return max((abs(h.cy_nm) + h.b_nm for h in geom.holes), default=0.0)


def slot_layout(layout: str, slot_nm: float, length_nm: float | None = None, neighbours: int = 1,
                single_w_nm: float = 80.0, notch_nm: float = 0.0) -> Geometry:
    """A 200 nm strip with slots slot_nm long (semicircular ends included) and 40 nm bridges, period slot_nm + 40.

      single   one centred column of single_w_nm wide slots (default 80: channels 60 | slot 80 | 60; 120 gives
               40 | 120 | 40 and a free width of 80 nm, key suffix _w120)
      double   two aligned columns of 40 nm slots at x = +/-40: 40 | 40 | 40 | 40 | 40
      stagger  as double, the right column shifted by half a period (glide symmetry)

    With the defaults all three have a narrowest free width of 120 nm.  A slot is centred on y = 0 (left column for stagger) with
    `neighbours` slots on each side of it; the shifted column of stagger has neighbours + 1 slots on each side of
    y = 0, so that its slots at +/- half a period have neighbours too.  The photon unit cell is |y| <= half a period.  length_nm
    None: the shortest device with LEAD_NM of plain strip beyond the outermost slot end on each side."""
    if layout not in LAYOUTS:
        raise ValueError(f"layout must be one of {LAYOUTS}")
    period = slot_nm + BRIDGE_NM
    b = slot_nm / 2.0
    ks = range(-neighbours, neighbours + 1)
    if layout == "single":
        a, cols = single_w_nm / 2.0, ((0.0, 0.0),)
    elif layout == "double":
        a, cols = 20.0, ((-40.0, 0.0), (40.0, 0.0))
    else:
        a, cols = 20.0, ((-40.0, 0.0), (40.0, period / 2.0))
    holes = []
    for xc, y0 in cols:
        if y0:          # the shifted column, at +/-(k + 1/2) period: every slot keeps a neighbour on both sides
            ys = [s * (k + 0.5) * period for k in range(neighbours + 1) for s in (-1, 1)]
        else:
            ys = [k * period for k in ks]
        holes += [Slot(xc, y, a, b) for y in ys]
    holes = tuple(sorted(holes, key=lambda h: (h.cx_nm, h.cy_nm)))
    # notch_nm > 0: a circular notch of that radius cut into each outer edge at every bridge row of the (left) column,
    # so the bridge-row cross-section (otherwise the full strip) narrows towards the free width at the slots
    # (only between two slots of a column; notch_nm <= BRIDGE_NM / 2 keeps it out of the slot ends)
    notches = ()
    if notch_nm > 0:
        if notch_nm > BRIDGE_NM / 2.0:
            raise ValueError(f"notch radius {notch_nm:g} > half the bridge: it would cut into the slot ends")
        y_max = max(abs(h.cy_nm) for h in holes)
        rows = sorted({h.cy_nm + s_ * (b + BRIDGE_NM / 2.0) for h in holes if h.cx_nm <= 0 for s_ in (-1, 1)})
        notches = tuple(Ellipse(sx * WIDTH_NM / 2.0, y, notch_nm, notch_nm) for y in rows if abs(y) < y_max
                        for sx in (-1, 1))
    half = max(abs(h.cy_nm) + h.b_nm for h in holes)
    if length_nm is None:
        length_nm = 2.0 * (half + LEAD_NM)
    if length_nm / 2.0 < half + 50.0:
        raise ValueError(f"length {length_nm:g} nm leaves no lead beyond the slots (half-extent {half:g} nm)")
    ys = []
    for y in sorted((0.0, b / 2.0, b - a, b, b + 10.0, b + 20.0)):     # rows closer than 5 nm merge
        if not ys or y - ys[-1] >= 5.0:
            ys.append(y)
    bridge_y = (b + 4.0, b + 12.0, b + 20.0)
    if layout == "single" and single_w_nm == 80.0:
        xs = (-96.0, -88.0, -78.0, -70.0, -60.0, -50.0, -44.0)
        bridge = [(x, y) for x in (-30.0, -15.0, 0.0) for y in bridge_y]
        critical = [(-96.0, 0.0), (-70.0, 0.0), (-44.0, 0.0), (-44.0, b - a), (-70.0, b + 20.0), (0.0, b + 4.0),
                    (0.0, b + 20.0)]
        symmetry = "column"
    elif layout == "single":
        edge = -(a + 4.0)                                   # 4 nm from the slot wall
        mid = round(-(100.0 + a) / 2.0)                     # channel centre
        xs = tuple(sorted({-96.0, -88.0, float(mid), round((mid + edge) / 2.0), edge}))
        bridge = [(x, y) for x in (round(-a / 2.0), round(-a / 4.0), 0.0) for y in bridge_y]
        critical = [(-96.0, 0.0), (float(mid), 0.0), (edge, 0.0), (edge, b - a), (float(mid), b + 20.0),
                    (0.0, b + 4.0), (0.0, b + 20.0)]
        symmetry = "column"
    else:
        xs = (-96.0, -88.0, -80.0, -72.0, -64.0, -16.0, -8.0, 0.0)
        bridge = [(x, y) for x in (-56.0, -48.0, -40.0, -32.0, -24.0) for y in bridge_y]
        critical = [(-96.0, 0.0), (-80.0, 0.0), (-64.0, 0.0), (-16.0, 0.0), (0.0, 0.0), (-40.0, b + 4.0),
                    (-80.0, b + 20.0), (-40.0, b + 20.0), (0.0, b + 20.0)]
        symmetry = "column" if layout == "double" else "zipper"
    channel = [(x, y) for x in xs for y in ys]
    if notch_nm > 0:          # sites in or within 4 nm of a notch move inward to 4 nm from it
        x_in = -(WIDTH_NM / 2.0 - notch_nm - 4.0)

        def off_notch(site):
            x, y = site
            return (x_in, y) if abs(abs(y) - (b + BRIDGE_NM / 2.0)) <= notch_nm + 4.0 and x < x_in else (x, y)
        channel = list(dict.fromkeys(off_notch(p_) for p_ in channel))
        bridge = list(dict.fromkeys(off_notch(p_) for p_ in bridge))
        critical = list(dict.fromkeys(off_notch(p_) for p_ in critical))
    probes = tuple(critical[i] for i in (0, 1, 2, len(critical) - 2))
    tag = "" if abs(length_nm - 2.0 * (half + LEAD_NM)) < 1e-6 else f"_{length_nm / 1000:g}um"
    return Geometry(
        key=f"{layout.capitalize()}_{slot_nm:g}" + (f"_w{single_w_nm:g}" if layout == "single" and single_w_nm != 80.0
                                                     else "") + (f"_n{notch_nm:g}" if notch_nm else "") + tag,
        label=f"{layout} column{'s' if layout != 'single' else ''}, {slot_nm:g} nm slots"
              + (f" {single_w_nm:g} nm wide" if layout == "single" and single_w_nm != 80.0 else "")
              + (f", {notch_nm:g} nm edge notches" if notch_nm else "")
              + (f", {length_nm / 1000:g} um" if tag else ""),
        width_nm=WIDTH_NM, length_nm=float(length_nm), holes=holes, probes=probes,
        grid_full=tuple(channel + bridge), grid_critical=tuple(critical), symmetry=symmetry,
        cell_half_y_nm=period / 2.0, crowd_half_y_nm=period / 2.0,
        hole_points=int(np.ceil(holes[0].perimeter() / HOLE_POINT_SPACING_NM)),
        probe_y_nm=length_nm / 2.0 - LEAD_NM / 2.0, color=C_HOLES, notches=notches)


# Plain references besides the 200 nm wire: the free width of the slot layouts (120 nm) and their channels (40, 60 nm)
PLAIN_120 = plain_strip(120.0, key="Plain_120", label="Plain 120 nm")
PLAIN_60 = plain_strip(60.0, key="Plain_60", label="Plain 60 nm")
PLAIN_40 = plain_strip(40.0, key="Plain_40", label="Plain 40 nm")
PLAIN_80 = plain_strip(80.0, key="Plain_80", label="Plain 80 nm")
SLOT_LAYOUTS = tuple(slot_layout(lay, L) for L in (150.0, 400.0) for lay in LAYOUTS)
# optimisation candidates of 7 Nov (fabrication limit: every feature >= 40 nm): longer slots (fewer bridge rows) and
# two 40 nm channels (a single column of 120 nm slots, free width 80 nm)
OPT_LAYOUTS = (slot_layout("double", 400.0, notch_nm=20.0), slot_layout("double", 600.0), slot_layout("double", 800.0), slot_layout("double", 1200.0), slot_layout("single", 400.0, single_w_nm=120.0),
               slot_layout("single", 800.0, single_w_nm=120.0))

GEOMETRIES = {g.key: g for g in (PLAIN_WIRE, PLAIN_WIRE_2UM, PLAIN_120, PLAIN_60, PLAIN_40, PLAIN_80) + SLOT_LAYOUTS + OPT_LAYOUTS}


def get(key: str, width_nm: float | None = None) -> Geometry:
    """Registered geometry by key; "Narrow_Wire" needs width_nm."""
    if key == "Narrow_Wire":
        if width_nm is None:
            raise ValueError("Narrow_Wire needs a width")
        return plain_strip(width_nm)
    return GEOMETRIES[key]
