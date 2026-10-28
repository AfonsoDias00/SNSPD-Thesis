import os
import csv
import tempfile
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import tdgl
from tdgl.geometry import box, ellipse

# ==============================================================================
# 1. Material & Device Parameters (Exact NbTiN setup)
# ==============================================================================
W = 0.200              # Wire width: 200 nm (in um)
L = 0.800              # Wire length: 800 nm (in um)
d = 0.009              # Film thickness: 9 nm
xi = 0.005             # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm

# 2D Diffusive Photon Hotspot Parameters
t0_photon = 20.0       # Photon arrival time [tau_0]
r_spot = 0.014         # Initial hotspot radius: 14 nm
amp_spot = -2.5        # Initial core suppression
D_tdgl = 1.0e-5        # 2D diffusion coefficient: 10 nm^2 / tau_0
t_relax = 35.0         # Quasiparticle relaxation time [tau_0]
t_max = 80.0           # Total simulation window [tau_0]
BISECT_TOL = 5.5       # Bisection tolerance [uA] (~5 steps per point for fast 2D grid)

OUT_DIR = "results_photon_hitmap_2d"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Build Geometries & 2D Irreducible Hitmap Grids
# ==============================================================================
def make_device(holes):
    layer = tdgl.Layer(coherence_length=xi, london_lambda=london_lambda, thickness=d, gamma=1.0)
    film = tdgl.Polygon("film", points=box(W, L))
    terminals = [
        tdgl.Polygon("source", points=box(W, L / 100, center=(0,  L / 2))),
        tdgl.Polygon("drain",  points=box(W, L / 100, center=(0, -L / 2))),
    ]
    device = tdgl.Device(
        "snspd",
        layer=layer,
        film=film,
        holes=holes if len(holes) > 0 else None,
        terminals=terminals,
        probe_points=[(0, L / 3), (0, -L / 3)],
        length_units="um",
    )
    device.make_mesh(max_edge_length=xi * 1.4, smooth=50)
    return device


def build_long_zipper_holes():
    holes = []
    # Left column: x = -36 nm, centers at y = -122, 0, +122 nm (by = 55 nm, ax = 8 nm)
    for i, y0 in enumerate([-0.122, 0.000, 0.122]):
        holes.append(
            tdgl.Polygon(f"l_l_{i}", points=ellipse(0.008, 0.055, center=(-0.036, y0), points=38)).resample(38)
        )
    # Right column: x = +36 nm, offset by dy/2 = 61 nm -> centers at y = -183, -61, +61, +183 nm
    for i, y0 in enumerate([-0.183, -0.061, 0.061, 0.183]):
        holes.append(
            tdgl.Polygon(f"l_r_{i}", points=ellipse(0.008, 0.055, center=( 0.036, y0), points=38)).resample(38)
        )
    return holes


def point_inside_any_hole(x, y, holes, margin=0.003):
    """Check if (x, y) is inside or within `margin` (3 nm) of any ellipse hole."""
    pt = np.array([[x, y]])
    for h in holes:
        if h.contains_points(pt)[0]:
            return True
        # Also check distance to hole boundary vertices
        dists = np.linalg.norm(h.points - pt, axis=1)
        if np.min(dists) < margin:
            return True
    return False


# ==============================================================================
# 3. Solver & Bisection Functions
# ==============================================================================
def make_hotspot_epsilon(x0, y0, t_start=t0_photon):
    def eps_func(r, *, t):
        r = np.asarray(r)
        if t < t_start:
            return np.ones(r.shape[:-1], dtype=float)
        dt = t - t_start
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2

        r0_sq = r_spot ** 2
        r_t_sq = r0_sq + 4.0 * D_tdgl * dt
        amp_t = amp_spot * (r0_sq / r_t_sq) * np.exp(-dt / t_relax)
        return 1.0 + amp_t * np.exp(-dist_sq / r_t_sq)
    return eps_func


def saved_steps_and_times(sol):
    with h5py.File(sol.path, "r") as h5_file:
        groups = sorted(h5_file["data"].keys(), key=int)
        times = np.array([h5_file["data"][group].attrs["time"] for group in groups])
    return np.array([int(group) for group in groups]), times


def run_simulation(device, current_uA, check_slices, photon_xy=None, save_every=120):
    fd, h5_path = tempfile.mkstemp(prefix="temp_hm-", suffix=".h5", dir=OUT_DIR)
    os.close(fd)
    os.unlink(h5_path)

    options = tdgl.SolverOptions(
        solve_time=t_max,
        dt_init=1e-3,
        dt_max=0.20,
        save_every=save_every,
        field_units="mT",
        current_units="uA",
        output_file=h5_path,
    )
    solve_kwargs = {}
    if photon_xy is not None:
        solve_kwargs["disorder_epsilon"] = make_hotspot_epsilon(*photon_xy, t_start=t0_photon)

    sol = None
    try:
        sol = tdgl.solve(
            device,
            options,
            terminal_currents=dict(source=current_uA, drain=-current_uA),
            **solve_kwargs,
        )
        t_arr = sol.dynamics.time
        v_arr = np.abs(sol.dynamics.voltage())

        post_mask = t_arr > (t0_photon + 4.0)
        max_v_post = np.max(v_arr[post_mask]) if np.any(post_mask) else 0.0
        mean_v_late = np.mean(v_arr[t_arr > 55.0]) if np.any(t_arr > 55.0) else 0.0

        # Check supercurrent collapse on the un-hit side of the wire
        saved_steps, saved_times = saved_steps_and_times(sol)
        post_saved_mask = saved_times > (t0_photon + 4.0)
        collapsed = False
        for x_min, x_max, y_sl in check_slices:
            x_pts = np.linspace(x_min, x_max, 30)
            pts = np.column_stack([x_pts, np.full_like(x_pts, y_sl)])
            dx = x_pts[1] - x_pts[0]
            curr_trace = []
            for s_step, s_time in zip(saved_steps, saved_times):
                if s_time > (t0_photon + 4.0):
                    sol.solve_step = int(s_step)
                    J_s = sol.interp_current_density(pts, dataset="supercurrent", units="uA/um")
                    curr_trace.append(np.trapezoid(np.abs(J_s[:, 1]), dx=dx))
            if len(curr_trace) > 0 and np.min(curr_trace) < 50.0:
                collapsed = True
                break

        switched = (max_v_post > 0.45) or (mean_v_late > 0.18) or collapsed
        return switched
    finally:
        output_path = getattr(sol, "path", None) or h5_path
        if output_path and os.path.exists(output_path):
            os.remove(output_path)


def find_threshold_bisection(device, check_slices, i_low, i_high, photon_xy=None, tol=BISECT_TOL):
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        if run_simulation(device, i_mid, check_slices, photon_xy=photon_xy):
            i_high = i_mid
        else:
            i_low = i_mid
    return i_high


# ==============================================================================
# 4. Main Execution: 1D Plain Wire Scan + 2D Long Zipper Hitmap
# ==============================================================================
if __name__ == "__main__":
    csv_path = os.path.join(OUT_DIR, "hitmap_results.csv")
    csv_rows = []

    # --- PART 1: Plain Wire 1D Transverse Scan (7 points from x = 0 to -88 nm) ---
    print("=" * 88)
    print("[1/2] Running Plain Wire (No Holes) Transverse Scan...")
    dev_plain = make_device([])
    slices_plain = [(0.002, 0.098, 0.000)]  # Monitor Right Half when hitting Left Half
    i_sw_plain = find_threshold_bisection(dev_plain, slices_plain, 1080.0, 1220.0, photon_xy=None)
    print(f" -> Plain Wire Spontaneous I_sw = {i_sw_plain:.1f} uA")

    x_plain_scan = [-0.088, -0.074, -0.058, -0.042, -0.028, -0.014, 0.000]
    plain_results = []
    for x0 in x_plain_scan:
        i_det = find_threshold_bisection(
            dev_plain, slices_plain, 0.78 * i_sw_plain, i_sw_plain, photon_xy=(x0, 0.0)
        )
        ratio = i_det / i_sw_plain
        plain_results.append((x0, 0.0, i_det, ratio))
        csv_rows.append(["Plain_Wire", x0 * 1e3, 0.0, i_sw_plain, i_det, ratio])
        print(f"    Plain Wire | x0 = {x0*1e3:6.1f} nm | I_det = {i_det:6.1f} uA | i_det = {ratio:.3f}")

    # --- PART 2: Long Zipper 2D Unit-Cell Hitmap ---
    print("\n" + "=" * 88)
    print("[2/2] Running Long Zipper (by=55nm) 2D Unit-Cell Hitmap...")
    holes_zip = build_long_zipper_holes()
    dev_zip = make_device(holes_zip)
    # Monitor Right Channel (Ch 3) at y = +61 nm (or Left Channel if hit is at x=0)
    slices_zip = [(0.046, 0.098, 0.061)]
    i_sw_zip = find_threshold_bisection(dev_zip, slices_zip, 960.0, 1080.0, photon_xy=None)
    print(f" -> Long Zipper Spontaneous I_sw = {i_sw_zip:.1f} uA")

    # 2D Grid over the left half of the unit cell:
    # x0 samples Ch 1 (-88, -68, -50 nm), the Left Bridge/Wall (-36 nm), and Ch 2 (-18, 0 nm)
    # y0 samples from y = 0 up to the half-period y = +61 nm
    x_grid_vals = [-0.088, -0.068, -0.050, -0.036, -0.018, 0.000]
    y_grid_vals = [0.000, 0.015, 0.030, 0.045, 0.061]

    zip_results = []
    for y0 in y_grid_vals:
        for x0 in x_grid_vals:
            if point_inside_any_hole(x0, y0, holes_zip, margin=0.003):
                print(f"    Long Zipper | ({x0*1e3:5.1f}, {y0*1e3:4.1f}) nm -> [Inside Hole - Skipped]")
                continue
            i_det = find_threshold_bisection(
                dev_zip, slices_zip, 0.78 * i_sw_zip, i_sw_zip, photon_xy=(x0, y0)
            )
            ratio = i_det / i_sw_zip
            zip_results.append((x0, y0, i_det, ratio))
            csv_rows.append(["Long_Zipper", x0 * 1e3, y0 * 1e3, i_sw_zip, i_det, ratio])
            print(f"    Long Zipper | ({x0*1e3:5.1f}, {y0*1e3:4.1f}) nm | I_det = {i_det:6.1f} uA | i_det = {ratio:.3f}")

    # Save CSV immediately
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Geometry", "x_nm", "y_nm", "I_sw_uA", "I_det_uA", "i_det"])
        writer.writerows(csv_rows)
    print(f"\nSaved raw hitmap table to: {csv_path}")

    # ==========================================================================
    # 5. Symmetrize 2D Hitmap Over Full Unit Cell (x in [-100, 100], y in [-61, 61] nm)
    # ==========================================================================
    # By Zipper symmetry:
    # 1. Reflection in y: (x, -y) has the exact same i_det as (x, +y)
    # 2. Glide reflection in x: (+x, y) is equivalent to (-x, |61 nm - y|)
    sym_pts = []
    sym_vals_norm = []
    sym_vals_abs = []

    for x0, y0, i_det, ratio in zip_results:
        # Left half: (x0, +y0) and (x0, -y0)
        for s_y in ([y0] if y0 == 0 else [y0, -y0]):
            sym_pts.append((x0 * 1e3, s_y * 1e3))
            sym_vals_norm.append(ratio)
            sym_vals_abs.append(i_det)
        # Right half via Zipper glide symmetry (x -> -x, y -> 0.061 - y)
        if x0 < 0:
            y_glide = 0.061 - y0
            for s_yg in ([y_glide] if y_glide == 0 else [y_glide, -y_glide]):
                sym_pts.append((-x0 * 1e3, s_yg * 1e3))
                sym_vals_norm.append(ratio)
                sym_vals_abs.append(i_det)

    sym_pts = np.array(sym_pts)
    sym_vals_norm = np.array(sym_vals_norm)
    sym_vals_abs = np.array(sym_vals_abs)

    # Interpolate onto a fine 2D unit-cell pixel grid for exact area histogram & heatmap
    xi_nm = np.linspace(-95, 95, 160)
    yi_nm = np.linspace(-61, 61, 120)
    X_nm, Y_nm = np.meshgrid(xi_nm, yi_nm)

    Z_zip_norm = griddata(sym_pts, sym_vals_norm, (X_nm, Y_nm), method="linear")
    Z_zip_norm_nn = griddata(sym_pts, sym_vals_norm, (X_nm, Y_nm), method="nearest")
    Z_zip_norm = np.where(np.isnan(Z_zip_norm), Z_zip_norm_nn, Z_zip_norm)

    Z_zip_abs = Z_zip_norm * i_sw_zip

    # Mask out the white ellipse holes on the fine 2D grid
    pts_um = np.column_stack([X_nm.ravel() * 1e-3, Y_nm.ravel() * 1e-3])
    hole_mask = np.zeros(len(pts_um), dtype=bool)
    for h in holes_zip:
        hole_mask |= h.contains_points(pts_um)
    hole_mask = hole_mask.reshape(X_nm.shape)
    Z_zip_norm[hole_mask] = np.nan
    Z_zip_abs[hole_mask] = np.nan

    # Plain wire 2D grid (depends only on |x|)
    x_pl_nm = np.array([abs(r[0]) * 1e3 for r in plain_results])
    idx_sort = np.argsort(x_pl_nm)
    x_pl_nm = x_pl_nm[idx_sort]
    norm_pl = np.array([r[3] for r in plain_results])[idx_sort]
    Z_plain_norm = np.interp(np.abs(X_nm), x_pl_nm, norm_pl)
    Z_plain_abs = Z_plain_norm * i_sw_plain

    # ==========================================================================
    # 6. Plot Figure 1: 2D Spatial Photon Hitmaps (Plain Wire vs. Long Zipper)
    # ==========================================================================
    fig_hm, (ax_p, ax_z) = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    vmin_hm, vmax_hm = 0.82, 0.96

    im_p = ax_p.pcolormesh(X_nm, Y_nm, Z_plain_norm, cmap="plasma_r", vmin=vmin_hm, vmax=vmax_hm, shading="auto")
    for x0, y0, _, _ in plain_results:
        ax_p.plot(x0 * 1e3, y0 * 1e3, "wo", markersize=6, markeredgecolor="k")
    ax_p.set_title(f"Plain Wire 2D Hitmap ($I_{{sw}} = {i_sw_plain:.0f}\\,\\mu$A)\nMean $i_{{det}} = {np.nanmean(Z_plain_norm):.3f}$",
                   fontweight="bold")
    ax_p.set_xlabel("$x_0$ [nm]")
    ax_p.set_ylabel("$y_0$ [nm] (One Unit Cell $P_y = 122$ nm)")
    ax_p.set_aspect("equal")

    im_z = ax_z.pcolormesh(X_nm, Y_nm, Z_zip_norm, cmap="plasma_r", vmin=vmin_hm, vmax=vmax_hm, shading="auto")
    for h in holes_zip:
        ax_z.fill(h.points[:, 0] * 1e3, h.points[:, 1] * 1e3, color="white", edgecolor="black", lw=1.2)
    for x0, y0, _, _ in zip_results:
        ax_z.plot(x0 * 1e3, y0 * 1e3, "wo", markersize=6, markeredgecolor="k")
    ax_z.set_title(f"Long Zipper 2D Hitmap ($I_{{sw}} = {i_sw_zip:.0f}\\,\\mu$A)\nMean $i_{{det}} = {np.nanmean(Z_zip_norm):.3f}$",
                   fontweight="bold")
    ax_z.set_xlabel("$x_0$ [nm]")
    ax_z.set_xlim(-100, 100)
    ax_z.set_ylim(-61, 61)
    ax_z.set_aspect("equal")

    cbar = fig_hm.colorbar(im_z, ax=[ax_p, ax_z], shrink=0.88, pad=0.03)
    cbar.set_label("Normalized Photon Threshold $i_{\\mathrm{det}}(x_0, y_0) = I_{\\mathrm{det}} / I_{\\mathrm{sw}}$\n(Brighter / Lower = More Sensitive)", fontsize=10.5)

    hm_path = os.path.join(OUT_DIR, "photon_hitmap_2d_comparison.png")
    fig_hm.savefig(hm_path, dpi=220, bbox_inches="tight")
    plt.close(fig_hm)
    print(f"Saved 2D photon hitmap comparison to: {hm_path}")

    # ==========================================================================
    # 7. Plot Figure 2: Exact Pixel-Integrated Efficiency Curves from 2D Hitmap
    # ==========================================================================
    fig_eff, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.5))
    i_norm_axis = np.linspace(0.75, 1.005, 300)
    i_abs_axis = np.linspace(750, 1210, 300)

    valid_plain_norm = Z_plain_norm[~np.isnan(Z_plain_norm)]
    valid_plain_abs = Z_plain_abs[~np.isnan(Z_plain_abs)]
    valid_zip_norm = Z_zip_norm[~np.isnan(Z_zip_norm)]
    valid_zip_abs = Z_zip_abs[~np.isnan(Z_zip_abs)]

    # Exact pixel fraction with i_det <= i_b (smoothed with 0.005 thermal broadening)
    eta_p_norm = [np.mean(0.5 * (1.0 + np.tanh((ib - valid_plain_norm) / 0.005))) * 100 for ib in i_norm_axis]
    eta_z_norm = [np.mean(0.5 * (1.0 + np.tanh((ib - valid_zip_norm) / 0.005))) * 100 for ib in i_norm_axis]

    eta_p_abs = [np.mean(0.5 * (1.0 + np.tanh((I - valid_plain_abs) / 5.0))) * 100 if I <= i_sw_plain else np.nan for I in i_abs_axis]
    eta_z_abs = [np.mean(0.5 * (1.0 + np.tanh((I - valid_zip_abs) / 5.0))) * 100 if I <= i_sw_zip else np.nan for I in i_abs_axis]

    ax1.plot(i_norm_axis, eta_p_norm, "k--", lw=2.5, label=f"Plain Wire (2D Integrated, $I_{{sw}}={i_sw_plain:.0f}\\,\\mu$A)")
    ax1.plot(i_norm_axis, eta_z_norm, "r-",  lw=2.6, label=f"Long Zipper (2D Hitmap Integrated, $I_{{sw}}={i_sw_zip:.0f}\\,\\mu$A)")
    ax1.axvline(1.0, color="gray", linestyle=":", label="Dark-count wall ($I_{sw}$)")
    ax1.set_xlabel("Normalized Bias Current $i_b = I_{\\mathrm{bias}} / I_{\\mathrm{sw}}$")
    ax1.set_ylabel("Active Area Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax1.set_title("True 2D Area-Integrated Normalized Efficiency")
    ax1.set_ylim(-5, 105)
    ax1.grid(True, alpha=0.4)
    ax1.legend(loc="lower right", fontsize=9)

    ax2.plot(i_abs_axis, eta_p_abs, "k--", lw=2.5, label="Plain Wire (No Holes)")
    ax2.plot(i_abs_axis, eta_z_abs, "r-",  lw=2.6, label="Long Zipper (by=55nm)")
    ax2.axvline(i_sw_plain, color="k", linestyle=":", alpha=0.75)
    ax2.axvline(i_sw_zip,   color="r", linestyle=":", alpha=0.75)
    ax2.set_xlabel("Absolute Bias Current $I_{\\mathrm{bias}}$ [$\\mu$A]")
    ax2.set_ylabel("Active Area Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax2.set_title("True 2D Area-Integrated Absolute Operating Window")
    ax2.set_ylim(-5, 105)
    ax2.grid(True, alpha=0.4)
    ax2.legend(loc="lower right", fontsize=9)

    plt.tight_layout()
    eff_path = os.path.join(OUT_DIR, "hitmap_efficiency_curves.png")
    fig_eff.savefig(eff_path, dpi=220)
    plt.close(fig_eff)
    print(f"Saved 2D-integrated efficiency curves to: {eff_path}")