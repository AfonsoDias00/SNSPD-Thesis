"""
2a  The fluxoid meter (snspd/fluxoid.py) on hand-built order parameters, on the real meshes around real holes, with
no TDGL.  This separates bugs in the meter from the physics.

  winding     psi = exp(i sum_k n_k phi_k) around every hole of the zipper, the 400 nm slot column and a single
              60 nm hole, n_k = -3..3, plus a uniform-current phase gradient: exactly n_k.
  noise       the same with independent random phase noise in [-0.3, 0.3] rad on every mesh site: still n_k.
  vortex      a vortex psi_v = (z - z_v) / sqrt(|z - z_v|^2 + 2 xi^2) between the hole wall and the loop counts +1
              (an antivortex -1); just outside the loop it counts 0.
  crossing    a vortex moved through the loop: the plan expects frames with the core on the loop to be rejected by
              the min |psi| >= PSI_LOOP_MIN gate.  Measured: the rejection band and the band in which the count
              disagrees with the core's side of the loop.
  belt        a normal hotspot on the loop is rejected; the same hotspot 11 nm off the loop is not, and counts right.
  bridge      a vortex in the middle of a 12 nm zipper bridge, where the loops of the two holes meet.
"""
from __future__ import annotations

import numpy as np

from snspd import geometry as geo
from snspd import sim
from snspd.fluxoid import PSI_LOOP_MIN, FluxoidMeter, winding_psi

from . import common as C

ID, TITLE = "2a", "fluxoid meter on synthetic fields"
R_HOLE = 30.0
MESH = 1.4


def vortex(xy_nm, z_v, sign=1):
    z = xy_nm[:, 0] + 1j * xy_nm[:, 1] - z_v
    return (z if sign > 0 else np.conj(z)) / np.sqrt(np.abs(z) ** 2 + 2 * C.XI ** 2)


def run(a):
    import tdgl
    from shapely.geometry import Point, Polygon
    rng = np.random.default_rng(1)
    checks = []
    devs = {"zipper": sim.make_device(geo.LONG_ZIPPER, MESH), "slot400": sim.make_device(geo.SLOT_400, MESH),
            "one_hole": sim.make_device(C.holed_strip(0.0, R_HOLE), MESH)}
    meters = {k: FluxoidMeter(d) for k, d in devs.items()}

    # winding, with and without noise and a transport phase gradient
    bad, bad_noise, n_cases = [], [], 0
    for k, dev in devs.items():
        nh = len(dev.holes)
        xy = np.asarray(dev.points) * 1e3
        for n in [tuple(((np.arange(nh) + s) % 7) - 3) for s in range(7)]:
            for q in (0.0, 0.3):                      # phase gradient along the strip [1/xi]: a transport current
                psi = winding_psi(dev, n) * np.exp(1j * q * xy[:, 1] / C.XI)
                n_cases += 1
                got, floor = meters[k](psi)
                if tuple(got) != n or floor < PSI_LOOP_MIN:
                    bad.append((k, n, q, got, floor))
                for _ in range(3):
                    got, _ = meters[k](psi * np.exp(1j * rng.uniform(-0.3, 0.3, len(psi))))
                    if tuple(got) != n:
                        bad_noise.append((k, n, q, got))
    checks.append(C.check("2a.1", "winding n = -3..3 around every hole (zipper, slot column, one hole), with and "
                          "without a transport phase gradient", f"{n_cases - len(bad)} of {n_cases} exact",
                          f"all exact and accepted (min |psi| >= {PSI_LOOP_MIN})", not bad, "; ".join(map(str, bad[:3]))))
    checks.append(C.check("2a.2", "the same with phase noise of +-0.3 rad on every site (3 draws each)",
                          f"{3 * n_cases - len(bad_noise)} of {3 * n_cases} exact", "all exact", not bad_noise,
                          "; ".join(map(str, bad_noise[:3]))))

    # one hole: vortices inside and outside the loop
    dev, meter = devs["one_hole"], meters["one_hole"]
    xy = np.asarray(dev.points) * 1e3
    loop = np.asarray(tdgl.make_fluxoid_polygons(dev, interp_points=400)[dev.holes[0].name]) * 1e3
    r_loop = float(np.hypot(*loop.T).mean())
    poly = Polygon(loop)
    inside_ok = outside_ok = True
    for n in (-1, 0, 1):
        base = winding_psi(dev, [n])
        for ang in np.linspace(0, 2 * np.pi, 8, endpoint=False):
            for s in (1, -1):
                got, floor = meter(base * vortex(xy, (R_HOLE + r_loop) / 2 * np.exp(1j * ang), s))
                inside_ok &= got[0] == n + s and floor >= PSI_LOOP_MIN
                got, floor = meter(base * vortex(xy, (r_loop + 2 * C.XI) * np.exp(1j * ang), s))
                outside_ok &= got[0] == n and floor >= PSI_LOOP_MIN
    checks.append(C.check("2a.3", f"(anti)vortex between the hole wall (r = {R_HOLE:g} nm) and the loop "
                          f"(r = {r_loop:.1f} nm), at r = {(R_HOLE + r_loop) / 2:.1f} nm, 8 angles, n_hole = -1, 0, 1",
                          "all +-1" if inside_ok else "some wrong", "n_hole +- 1, accepted", inside_ok))
    checks.append(C.check("2a.4", f"(anti)vortex 2 xi outside the loop (r = {r_loop + 2 * C.XI:.1f} nm)",
                          "all unchanged" if outside_ok else "some wrong", "n_hole, accepted", outside_ok))

    # a vortex moved through the loop along 8 radial lines
    edge = C.XI * MESH
    rows = []
    for ang in np.linspace(0, 2 * np.pi, 8, endpoint=False) + 0.1:
        for d in np.arange(-10.0, 10.01, 0.25):
            z_v = (r_loop + d) * np.exp(1j * ang)
            p = Point(z_v.real, z_v.imag)
            dist = poly.exterior.distance(p) * (1 if poly.contains(p) else -1)     # > 0 inside the loop
            got, floor = meter(vortex(xy, z_v))
            rows.append((dist, got[0], floor))
    rows = np.array(rows)
    on = np.abs(rows[:, 0]) <= 0.5
    rej = rows[:, 2] < PSI_LOOP_MIN
    expect = (rows[:, 0] > 0).astype(int)
    wrong = ~rej & (rows[:, 1] != expect)
    band_rej = float(np.abs(rows[rej, 0]).max()) if rej.any() else 0.0
    band_wrong = float(np.abs(rows[wrong, 0]).max()) if wrong.any() else 0.0
    checks.append(C.check("2a.5", "vortex core within 0.5 nm of the loop: frame rejected (expectation in the plan)",
                          f"{int((on & rej).sum())} of {int(on.sum())} rejected",
                          f"all rejected (min |psi| < {PSI_LOOP_MIN})", bool(rej[on].all()),
                          f"rejected frames have the core at most {band_rej:.2f} nm from the loop; loop points are "
                          f"{2 * np.pi * r_loop / 400:.2f} nm apart, mesh edges ~{edge:.1f} nm"))
    out_of_set = ~rej & ~np.isin(rows[:, 1], (0, 1))
    checks.append(C.check("2a.6", "vortex moved through the loop: accepted frames count the core's side, except "
                          "within half a mesh edge of the loop", f"wrong side up to {band_wrong:.2f} nm from the loop "
                          f"({int(wrong.sum())} of {int((~rej).sum())} accepted); counts outside {{0, 1}}: "
                          f"{int(out_of_set.sum())}", f"wrong only within {edge / 2:.1f} nm, never outside {{0, 1}}",
                          band_wrong <= edge / 2 and not out_of_set.any(),
                          "a vortex crossing the loop moves the count by exactly 1; the gate catches it only if a "
                          "frame falls while the core is within ~0.5 nm of a loop point"))

    # a normal hotspot on the loop, and next to it
    def belt(cx, cy, r=14.0, floor=1e-6):
        d = np.hypot(xy[:, 0] - cx, xy[:, 1] - cy)
        return np.clip((d - r) / C.XI, floor, 1.0)

    base = winding_psi(dev, [1])
    got_on, floor_on = meter(base * belt(r_loop, 0.0))
    got_off, floor_off = meter(base * belt(r_loop + 25.0, 0.0))
    checks.append(C.check("2a.7", "normal disc (r = 14 nm, |psi| = 1e-6 inside) centred on the loop / 11 nm clear "
                          "of it", f"min |psi| {floor_on:.1e} / {floor_off:.3f}, count off-loop {got_off[0]}",
                          f"rejected / accepted with n = 1", floor_on < PSI_LOOP_MIN <= floor_off and got_off[0] == 1))
    try:
        got, floor = meter(base * belt(r_loop, 0.0, floor=0.0))
        res, ok = f"returns n = {got[0]}, min |psi| = {floor:.1e}", floor < PSI_LOOP_MIN
    except ValueError as exc:
        res, ok = f"raises ValueError ({exc})", False
    checks.append(C.check("2a.8", "the same disc with psi exactly 0 on the loop", res,
                          "a rejected frame (min |psi| = 0), not an exception", ok,
                          "FluxoidMeter.__call__ divides by the loop samples (np.roll(z, -1) / z) and rounds the NaN; "
                          "TDGL frames are never exactly 0 off the terminals, so production runs are not affected"))

    # a vortex in the middle of a zipper bridge: on the loops of both holes of the column
    dev, meter = devs["zipper"], meters["zipper"]
    xyz = np.asarray(dev.points) * 1e3
    step = 0.25
    both = []
    for xv in np.arange(-44.0, -27.99, step):
        for yv in np.arange(55.0, 67.01, step):
            got, floor = meter(vortex(xyz, complex(xv, yv)))
            if got[1] == 1 and got[2] == 1 and floor >= PSI_LOOP_MIN:
                both.append((xv, yv))
    both = np.array(both).reshape(-1, 2)
    ext = "" if not len(both) else (f"x {both[:, 0].min():g}..{both[:, 0].max():g}, "
                                    f"y {both[:, 1].min():g}..{both[:, 1].max():g} nm")
    checks.append(C.check("2a.9", "vortex anywhere in the 12 nm bridge between zipper holes 1 and 2 "
                          "(x = -44..-28, y = 55..67 nm, 0.25 nm grid): counted in both holes in an accepted frame",
                          f"{len(both) * step ** 2:.1f} nm^2 of double counting" + (f" ({ext})" if ext else ""),
                          "never (each vortex counted at most once)", not len(both),
                          "the two fluxoid loops (each pushed 6 nm out, half the bridge) overlap at the bridge middle; "
                          "a frame taken while a vortex crosses there adds a spurious quantum for one frame"))
    return checks, dict(r_loop_nm=r_loop, crossing=rows.tolist(), bridge_double=both.tolist())
