import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
import tdgl
from tdgl.geometry import box, ellipse

# ==============================================================================
# 1. Material & Device Parameters (Adapted to your setup)
# ==============================================================================
W = 0.200              # Wire width: 200 nm (in um)
L = 0.800              # Wire length: 800 nm (in um)
d = 0.009              # Film thickness: 9 nm (set to 0.005 if using 5 nm)
xi = 0.005             # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm


I_BIAS_STATIC = 850.0  # ~1530 uA for d=9nm (850 uA for d=5nm)
I_SW_LOW = 880.0       # Lower bisection bound for spontaneous I_sw
I_SW_HIGH = 965.0       # Upper bisection bound for spontaneous I_sw
BISECT_TOL = 4.0       # Bisection convergence tolerance [uA]

T_STATIC = 25.0        # Relaxation time for static maps [tau_0]

# Hotspot parameters for dynamic sweep
r_spot = 0.018         # Hotspot radius: 18 nm
amp_spot = -1.0        # Core suppression
t_relax = 15.0         # Thermal relaxation time [tau_0]
t_max = 60.0           # Dynamic simulation window [tau_0]

OUT_DIR = "results_staggered_verification"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Define Geometries & 6-Point Complete Unit-Cell Grid
# ==============================================================================
def build_geometries():
    geoms = {}

    # 1. Baseline 2x2 Ellipses (40 nm vertical gap)
    holes_base = [
        tdgl.Polygon(f"h_{i}", points=ellipse(0.012, 0.025, center=(x0, y0), points=28)).resample(28)
        for i, (x0, y0) in enumerate([
            (-0.042,  0.045), (0.042,  0.045),
            (-0.042, -0.045), (0.042, -0.045),
        ])
    ]
    pts_base = [
        (0.000, 0.045, 0.25, "1_row_mid_(0,45nm)"),
        (0.077, 0.045, 0.25, "2_row_edge_(77,45nm)"),
        (0.020, 0.022, 0.20, "3_saddle_(20,22nm)"),
        (0.000, 0.000, 0.15, "4_center_(0,0nm)"),
        (0.042, 0.012, 0.08, "5_hole_tip_(42,12nm)"),
        (0.085, 0.000, 0.07, "6_far_edge_(85,0nm)"),
    ]
    geoms["Baseline 2x2\n(40 nm Gap)"] = (holes_base, pts_base, "k--o", 0.045)

    # 2. Original Width-Balanced Staggered (ax_c = 18 nm, dy = 38 nm)
    holes_stag18 = [
        tdgl.Polygon("h_ul", points=ellipse(0.010, 0.024, center=(-0.046,  0.038), points=28)).resample(28),
        tdgl.Polygon("h_ur", points=ellipse(0.010, 0.024, center=( 0.046,  0.038), points=28)).resample(28),
        tdgl.Polygon("h_c",  points=ellipse(0.018, 0.024, center=( 0.000,  0.000), points=32)).resample(32),
        tdgl.Polygon("h_ll", points=ellipse(0.010, 0.024, center=(-0.046, -0.038), points=28)).resample(28),
        tdgl.Polygon("h_lr", points=ellipse(0.010, 0.024, center=( 0.046, -0.038), points=28)).resample(28),
    ]
    pts_stag18 = [
        (0.000, 0.038, 0.20, "1_outer_mid_(0,38nm)"),
        (0.078, 0.038, 0.22, "2_outer_edge_(78,38nm)"),
        (0.028, 0.019, 0.22, "3_diag_saddle_(28,19nm)"),
        (0.059, 0.000, 0.18, "4_center_ch_(59,0nm)"),
        (0.000, 0.029, 0.08, "5_tip_wake_(0,29nm)"),
        (0.085, 0.000, 0.10, "6_far_edge_(85,0nm)"),
    ]
    geoms["Staggered 2-1-2\n(ax_c = 18 nm)"] = (holes_stag18, pts_stag18, "r-s", 0.038)

    # 3. Tuned Staggered 2-1-2 (ax_c = 15 nm, dy = 40 nm -> Recovers I_sw)
    holes_stag15 = [
        tdgl.Polygon("h_ul", points=ellipse(0.010, 0.024, center=(-0.046,  0.040), points=28)).resample(28),
        tdgl.Polygon("h_ur", points=ellipse(0.010, 0.024, center=( 0.046,  0.040), points=28)).resample(28),
        tdgl.Polygon("h_c",  points=ellipse(0.015, 0.024, center=( 0.000,  0.000), points=30)).resample(30),
        tdgl.Polygon("h_ll", points=ellipse(0.010, 0.024, center=(-0.046, -0.040), points=28)).resample(28),
        tdgl.Polygon("h_lr", points=ellipse(0.010, 0.024, center=( 0.046, -0.040), points=28)).resample(28),
    ]
    pts_stag15 = [
        (0.000, 0.040, 0.20, "1_outer_mid_(0,40nm)"),
        (0.078, 0.040, 0.22, "2_outer_edge_(78,40nm)"),
        (0.027, 0.020, 0.22, "3_diag_saddle_(27,20nm)"),
        (0.057, 0.000, 0.18, "4_center_ch_(57,0nm)"),
        (0.000, 0.029, 0.08, "5_tip_wake_(0,29nm)"),
        (0.085, 0.000, 0.10, "6_far_edge_(85,0nm)"),
    ]
    geoms["Tuned Staggered 2-1-2\n(ax_c = 15 nm)"] = (holes_stag15, pts_stag15, "g-^", 0.040)

    return geoms


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
        holes=holes,
        terminals=terminals,
        probe_points=[(0, L / 3), (0, -L / 3)],
        length_units="um",
    )
    device.make_mesh(max_edge_length=xi * 1.4, smooth=50)
    return device


# ==============================================================================
# 3. Helper Functions for Dynamic Photon Bisection (Adapted to your pytdgl API)
# ==============================================================================
def make_hotspot_epsilon(x0, y0):
    def eps_func(r, *, t):
        r = np.asarray(r)
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2
        return 1.0 + amp_spot * np.exp(-dist_sq / (r_spot ** 2)) * np.exp(-t / t_relax)
    return eps_func


def did_switch(device, current_uA, photon_xy=None):
    options = tdgl.SolverOptions(
        solve_time=t_max,
        dt_init=1e-3,
        dt_max=0.25,
        save_every=200,
        field_units="mT",
        current_units="uA",
        output_file=None,
    )
    solve_kwargs = {}
    if photon_xy is not None:
        solve_kwargs["disorder_epsilon"] = make_hotspot_epsilon(*photon_xy)

    sol = tdgl.solve(
        device,
        options,
        terminal_currents=dict(source=current_uA, drain=-current_uA),
        **solve_kwargs,
    )
    if sol is None or sol.dynamics is None:
        raise RuntimeError("TDGL returned no dynamics data.")
    t_arr = sol.dynamics.time
    v_arr = np.abs(sol.dynamics.voltage())
    late_v = np.mean(v_arr[t_arr > 0.6 * t_max])
    return late_v > 0.05, sol


def find_threshold_bisection(device, i_low, i_high, photon_xy=None, tol=BISECT_TOL):
    best_sol = None
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        switched, sol = did_switch(device, i_mid, photon_xy=photon_xy)
        if switched:
            i_high = i_mid
            best_sol = sol
        else:
            i_low = i_mid
    return i_high, best_sol


# ==============================================================================
# 4. Main Execution: Part A (Static 2D Maps & 1D Slices) -> Part B (6-Pt Sweep)
# ==============================================================================
if __name__ == "__main__":
    geoms = build_geometries()
    devices = {}
    static_sols = {}

    print(f"\n[PART A] Solving Static |psi| & |J| Maps at I = {I_BIAS_STATIC:.1f} uA...")
    
    for g_name, (holes, pts, style, y_row) in geoms.items():
        clean_name = g_name.replace("\n", " ")
        slug = clean_name.split("(")[0].strip().replace(" ", "_").replace("-", "_")
        print(f" -> Meshing & solving static state for: {clean_name}")

        dev = make_device(holes)
        devices[g_name] = dev

        static_opts = tdgl.SolverOptions(
            solve_time=T_STATIC,
            dt_init=1e-3,
            dt_max=0.25,
            save_every=100,
            field_units="mT",
            current_units="uA",
            output_file=os.path.join(OUT_DIR, f"{slug}_static.h5"),
        )
        sol = tdgl.solve(
            dev,
            static_opts,
            terminal_currents=dict(source=I_BIAS_STATIC, drain=-I_BIAS_STATIC),
        )
        static_sols[g_name] = sol

        # Save individual built-in plots
        fig_psi, _ = sol.plot_order_parameter(figsize=(7, 3.5))
        fig_psi.suptitle(f"{clean_name} at I = {I_BIAS_STATIC:.0f} uA", fontsize=11, y=1.02)
        fig_psi.savefig(os.path.join(OUT_DIR, f"{slug}_order_parameter.png"), dpi=200, bbox_inches="tight")
        plt.close(fig_psi)

        fig_j, _ = sol.plot_currents(figsize=(7, 3.5), streamplot=True)
        fig_j.suptitle(f"{clean_name} Supercurrent J at I = {I_BIAS_STATIC:.0f} uA", fontsize=11, y=1.02)
        fig_j.savefig(os.path.join(OUT_DIR, f"{slug}_currents.png"), dpi=200, bbox_inches="tight")
        plt.close(fig_j)

    # --- Plot 2x3 Combined |psi| and |J| + Streamlines Figure ---
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9), sharex=True, sharey=True)

    for col, (g_name, (holes, pts, style, y_row)) in enumerate(geoms.items()):
        dev = devices[g_name]
        sol = static_sols[g_name]
        sol.solve_step = -1
        points = dev.points
        triangles = dev.triangles

        psi_mag = np.abs(sol.tdgl_data.psi)
        J_vec = sol.interp_current_density(
            points, dataset="supercurrent", units="uA/um"
        )
        J_mag = np.linalg.norm(J_vec, axis=1)

        # Top row: |psi| + 6 Impact Points
        ax_psi = axes[0, col]
        im_psi = ax_psi.tripcolor(
            points[:, 0], points[:, 1], psi_mag,
            triangles=triangles, shading="gouraud", cmap="viridis", vmin=0.55, vmax=1.0
        )
        for hole in dev.holes:
            ax_psi.fill(hole.points[:, 0], hole.points[:, 1], color="white", edgecolor="black", lw=0.8)

        for idx, (x0, y0, _, lbl) in enumerate(pts, 1):
            ax_psi.plot(x0, y0, "r*", markersize=9, markeredgecolor="k", markeredgewidth=0.5)
            ax_psi.annotate(
                f"P{idx}", (x0, y0), xytext=(3, 3), textcoords="offset points",
                fontsize=8, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.7, lw=0)
            )

        ax_psi.set_title(g_name, fontsize=11, fontweight="bold")
        ax_psi.set_xlim(-W / 2, W / 2)
        ax_psi.set_ylim(-0.12, 0.12)
        ax_psi.set_aspect("equal")
        if col == 0:
            ax_psi.set_ylabel("$y$ [$\\mu$m]", fontsize=11)

        # Bottom row: |J| + Cyan Streamlines
        ax_j = axes[1, col]
        im_j = ax_j.tripcolor(
            points[:, 0], points[:, 1], J_mag,
            triangles=triangles, shading="gouraud", cmap="inferno", vmin=0, vmax=np.percentile(J_mag, 99)
        )
        for hole in dev.holes:
            ax_j.fill(hole.points[:, 0], hole.points[:, 1], color="white", edgecolor="black", lw=0.8)

        # Explicit meshgrid for streamplot (compatible with all pytdgl versions)
        x_values = np.linspace(-W / 2, W / 2, 160)
        y_values = np.linspace(-0.12, 0.12, 160)
        x_grid, y_grid = np.meshgrid(x_values, y_values)
        Jx_i = LinearNDInterpolator(points, J_vec[:, 0])(x_grid, y_grid)
        Jy_i = LinearNDInterpolator(points, J_vec[:, 1])(x_grid, y_grid)

        grid_pts = np.column_stack([x_grid.ravel(), y_grid.ravel()])
        for hole in dev.holes:
            mask = hole.contains_points(grid_pts).reshape(x_grid.shape)
            Jx_i[mask] = np.nan
            Jy_i[mask] = np.nan

        ax_j.streamplot(
            x_grid[0, :], y_grid[:, 0], Jx_i, Jy_i,
            color="cyan", density=1.1, linewidth=0.65, arrowsize=0.7
        )
        ax_j.set_xlim(-W / 2, W / 2)
        ax_j.set_ylim(-0.12, 0.12)
        ax_j.set_aspect("equal")
        ax_j.set_xlabel("$x$ [$\\mu$m]", fontsize=11)
        if col == 0:
            ax_j.set_ylabel("$y$ [$\\mu$m]", fontsize=11)

    cbar_psi = fig.colorbar(im_psi, ax=axes[0, :], shrink=0.85, pad=0.02)
    cbar_psi.set_label("Order Parameter $|\\psi|$", fontsize=11)
    cbar_j = fig.colorbar(im_j, ax=axes[1, :], shrink=0.85, pad=0.02)
    cbar_j.set_label("Sheet Current Density $|\\mathbf{J}|$ [$\\mu$A/$\\mu$m]", fontsize=11)

    fig.suptitle(f"Order Parameter $|\\psi|$ and Supercurrent $|\\mathbf{{J}}|$ at $I_{{\\mathrm{{bias}}}} = {I_BIAS_STATIC:.0f}\\,\\mu$A", fontsize=13, y=0.98)
    map_path = os.path.join(OUT_DIR, "combined_psi_and_J_maps.png")
    fig.savefig(map_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f" -> Saved 2x3 thesis figure to: {map_path}")

    # --- Plot 1D Transverse Slices at y = 0 and y = y_row ---
    fig1d, (ax_y0, ax_yr) = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    x_line = np.linspace(-W / 2 + 0.002, W / 2 - 0.002, 250)

    for g_name, (holes, pts, style, y_row) in geoms.items():
        clean_name = g_name.replace("\n", " ")
        dev = devices[g_name]
        sol = static_sols[g_name]
        points = dev.points

        J_vec = sol.interp_current_density(
            points, dataset="supercurrent", units="uA/um"
        )
        J_mag = np.linalg.norm(J_vec, axis=1)
        J_interp = LinearNDInterpolator(points, J_mag)

        pts_y0 = np.column_stack([x_line, np.zeros_like(x_line)])
        j_y0 = J_interp(pts_y0)
        for hole in dev.holes:
            j_y0[hole.contains_points(pts_y0)] = np.nan
        ax_y0.plot(x_line * 1e3, j_y0, style[:2], lw=2.2, label=clean_name)

        pts_yr = np.column_stack([x_line, np.full_like(x_line, y_row)])
        j_yr = J_interp(pts_yr)
        for hole in dev.holes:
            j_yr[hole.contains_points(pts_yr)] = np.nan
        ax_yr.plot(x_line * 1e3, j_yr, style[:2], lw=2.2, label=f"{clean_name} ($y={y_row*1e3:.0f}$ nm)")

    ax_y0.set_title("Transverse Current Profile at $y = 0$ nm", fontweight="bold")
    ax_y0.set_xlabel("$x$ [nm]")
    ax_y0.set_ylabel("Sheet Current Density $|\\mathbf{J}|$ [$\\mu$A/$\\mu$m]")
    ax_y0.grid(True, alpha=0.4)
    ax_y0.legend(fontsize=8.5)

    ax_yr.set_title("Transverse Current Profile at Outer Hole Row ($y = y_{\\mathrm{row}}$)", fontweight="bold")
    ax_yr.set_xlabel("$x$ [nm]")
    ax_yr.grid(True, alpha=0.4)
    ax_yr.legend(fontsize=8.5)

    plt.tight_layout()
    slice_path = os.path.join(OUT_DIR, "transverse_J_profiles.png")
    fig1d.savefig(slice_path, dpi=200)
    plt.close(fig1d)
    print(f" -> Saved 1D transverse slices to: {slice_path}\n")

    # ==========================================================================
    # PART B: 6-Point Dynamic Photon Bisection Sweep (Checking Points 5 & 6)
    # ==========================================================================
    print("=" * 88)
    print(f"{'Geometry':<34} | {'Point':<26} | {'I_det [uA]':<10} | {'i_det'}")
    print("-" * 88)

    summary = {}
    for g_name, (holes, pts, style, _) in geoms.items():
        clean_name = g_name.replace("\n", " ")
        dev = devices[g_name]
        i_sw, _ = find_threshold_bisection(dev, I_SW_LOW, I_SW_HIGH, photon_xy=None, tol=BISECT_TOL)
        print(f"{clean_name:<34} | {'Spontaneous I_sw (No Photon)':<26} | {i_sw:<10.1f} | 1.000")

        pt_results = []
        for x0, y0, weight, label in pts:
            i_det, _ = find_threshold_bisection(dev, 0.78 * i_sw, i_sw, photon_xy=(x0, y0), tol=BISECT_TOL)
            ratio = i_det / i_sw
            pt_results.append((i_det, ratio, weight, label))
            print(f"{'':<34} | {label:<26} | {i_det:<10.1f} | {ratio:.3f}")
        print("-" * 88)
        summary[clean_name] = {"I_sw": i_sw, "points": pt_results, "style": style}

    # --- Plot 6-Point Efficiency Plateau Comparison ---
    fig_eff, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    i_norm_grid = np.linspace(0.80, 1.01, 300)
    i_abs_grid = np.linspace(0.78 * I_SW_LOW, 1.02 * I_SW_HIGH, 300)

    for clean_name, data in summary.items():
        i_sw = data["I_sw"]
        pts = data["points"]
        style = data["style"]

        eta_norm = np.zeros_like(i_norm_grid)
        eta_abs = np.zeros_like(i_abs_grid)
        for i_det, ratio, w, _ in pts:
            eta_norm += w * 0.5 * (1.0 + np.tanh((i_norm_grid - ratio) / 0.008))
            eta_abs  += w * 0.5 * (1.0 + np.tanh((i_abs_grid - i_det) / 7.5 ))

        eta_norm[i_norm_grid > 1.0] = np.nan
        eta_abs[i_abs_grid > i_sw] = np.nan

        ax1.plot(i_norm_grid, eta_norm * 100, style, lw=2.4, markevery=35,
                 label=f"{clean_name} ($I_{{sw}}={i_sw:.0f}\\,\\mu$A)")
        ax2.plot(i_abs_grid, eta_abs * 100, style, lw=2.4, markevery=35,
                 label=f"{clean_name}")
        ax2.axvline(i_sw, color=style[0], linestyle=":", alpha=0.75)

    ax1.axvline(1.0, color="gray", linestyle="--", label="Dark Count Limit ($I_{sw}$)")
    ax1.set_xlabel("Normalized Bias Current $i_b = I_{\\mathrm{bias}} / I_{\\mathrm{sw}}$")
    ax1.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax1.set_title("6-Point Normalized Efficiency Plateau ($\\eta_{\\mathrm{int}}$ vs. $I_b/I_{sw}$)")
    ax1.set_ylim(-5, 105)
    ax1.grid(True, alpha=0.4)
    ax1.legend(loc="lower right", fontsize=8.5)

    ax2.set_xlabel("Absolute Bias Current $I_{\\mathrm{bias}}$ [$\\mu$A]")
    ax2.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax2.set_title("6-Point Absolute Operating Window (Plateau vs. $I_{sw}$ Wall)")
    ax2.set_ylim(-5, 105)
    ax2.grid(True, alpha=0.4)
    ax2.legend(loc="lower right", fontsize=8.5)

    plt.tight_layout()
    eff_path = os.path.join(OUT_DIR, "staggered_6pt_efficiency_plateau.png")
    fig_eff.savefig(eff_path, dpi=200)
    plt.close(fig_eff)
    print(f"Saved 6-point efficiency plateau plot to: {eff_path}")