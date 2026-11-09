"""
Graded meshes: the production edge length over the whole slotted section, a coarser one only in the plain leads
between the outermost slot end and the terminals, where the film just carries the transport current.

Lesson of 6 Nov 2026 (first version): a fine band around the photon cell only left the outer slots of a 2 um slot
column in the coarse zone, and the dark I_sw fell from 0.97 to 0.83 I_dep, because dark switching nucleates at the
hole rims.  So the fine region is |y| <= slotted_half_nm(geom) + band_margin_nm (every hole rim and bridge), plus
contact_nm at each terminal (the psi = 0 terminals also nucleate dark switching).  In between, the target edge ramps
linearly to coarse_xi over ramp_nm.  Triangle (meshpy, as pyTDGL) refines every triangle larger than the local target
area (sqrt(3)/4 h^2) with pyTDGL's minimum angle, and the result gets pyTDGL's Laplacian smoothing.

The finest cells set the time step (run time ~ h^-4 for uniform meshes), so grading only removes nodes: expect at
most the node ratio as speed-up.
"""
from __future__ import annotations

import numpy as np

SQ3_4 = np.sqrt(3.0) / 4.0
# pyTDGL's max_edge_length caps the longest edge; Triangle with an area target sets the typical one.  The production
# mesh at max edge 1.4 xi = 7 nm has a mean triangle of 7.1 nm^2 (plain wire: 160 000 nm^2 / (2 x 11 257 nodes)), an
# equilateral edge of 4.05 nm = 0.58 x 7 nm.  Triangle refines every triangle above the target, so its triangles end up
# ~0.6 x the target area: 0.58 / sqrt(0.6) = 0.75 matches the production density in the photon cell (checked 6 Nov).
EDGE_PER_MAX = 0.75


def target_edge_nm(y_nm, geom, xi_nm, fine_xi, coarse_xi, band_margin_nm=40.0, ramp_nm=100.0, contact_nm=60.0):
    """Target edge length [nm] at height y_nm (array), in the production convention (max edge = mesh_xi * xi)."""
    from snspd.geometry import slotted_half_nm
    y = np.abs(np.asarray(y_nm, dtype=float))
    band = max(slotted_half_nm(geom), geom.cell_half_y_nm) + band_margin_nm
    to_contact = geom.length_nm / 2.0 - y
    # distance (>= 0) out of the fine band, and out of the fine contact strip
    out_band = np.maximum(y - band, 0.0)
    out_contact = np.maximum(to_contact - contact_nm, 0.0)
    d = np.minimum(out_band, out_contact)
    w = np.clip(d / ramp_nm, 0.0, 1.0)
    return xi_nm * (fine_xi + (coarse_xi - fine_xi) * w)


def graded_mesh(device, geom, xi_nm, fine_xi, coarse_xi, smooth=50, min_angle=32.5, **kw):
    """(points [um], triangles) of a graded mesh for a tdgl.Device built by snspd.sim.make_device (its film and hole
    polygons are used as given)."""
    import meshpy.triangle as triangle
    from shapely.geometry import Polygon
    from tdgl.finite_volume.mesh import Mesh

    from tdgl.geometry import ensure_unique                  # duplicate vertices segfault triangle.build()
    film = ensure_unique(np.asarray(device.film.points, dtype=float))
    holes = [ensure_unique(np.asarray(h.points, dtype=float)) for h in device.holes]
    coords = np.concatenate([film] + holes)
    facets, start = [], 0
    for ring in [film] + holes:
        idx = np.arange(start, start + len(ring))
        facets += list(zip(idx, np.roll(idx, -1)))
        start += len(ring)
    info = triangle.MeshInfo()
    info.set_points(coords)
    info.set_facets(facets)
    if holes:
        info.set_holes([tuple(Polygon(h).representative_point().coords[0]) for h in holes])

    def refine(vertices, area):
        y_nm = 1e3 * (vertices[0][1] + vertices[1][1] + vertices[2][1]) / 3.0
        h_um = 1e-3 * EDGE_PER_MAX * float(target_edge_nm(y_nm, geom, xi_nm, fine_xi, coarse_xi, **kw))
        return bool(area > SQ3_4 * h_um ** 2)

    mesh = triangle.build(info, refinement_func=refine, min_angle=min_angle)
    points, triangles = np.array(mesh.points), np.array(mesh.elements)
    if smooth:
        m = Mesh.from_triangulation(points, triangles, create_submesh=False).smooth(smooth, create_submesh=False)
        points, triangles = m.sites, m.elements
    return points, triangles
