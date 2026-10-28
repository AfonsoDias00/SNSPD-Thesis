import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
import tdgl
from tdgl.geometry import box, ellipse

# Use non-interactive backend to avoid FigureCanvasAgg warnings
import matplotlib
matplotlib.use("Agg")

# ==============================================================================
# 1. Material & Device Parameters (Matching efficiency_curves.py)
# ==============================================================================
W = 0.200              # Wire width: 200 nm (in um)
L = 0.800              # Wire length: 400 nm (in um)
d = 0.009              # Film thickness: 9 nm
xi = 0.005             # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm

I_BIAS_STATIC = 850.0  # Sub-critical bias [uA] (~0.94 I_sw) to inspect static crowding
T_RELAX = 25.0         # Relaxation time to reach steady state (no vortices)

OUT_DIR = "results_plateau_maps"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Build the 3 Geometries (with points=28 to prevent Voronoi mesh warnings)
# ==============================================================================
def get_geometries():
    geoms = {}

    # 1. Baseline 2x2 (40 nm vertical gap at y = 0)
    holes_base = [
        tdgl.Polygon(f"h_{i}", points=ellipse(0.012, 0.025, center=(x0, y0), points=28)).resample(28)
        for i, (x0, y0) in enumerate([
            (-0.042,  0.045), (0.042,  0.045),
            (-0.042, -0.045), (0.042, -0.045),
        ])
    ]
    pts_base = [
        (0.000, 0.045, "1 (907 uA)"),
        (0.077, 0.045, "2 (911 uA)"),
        (0.020, 0.022, "3 (923 uA)"),
        (0.000, 0.000, "4 (923 uA)"),
    ]
    geoms["Baseline 2x2\n(40 nm Gap)"] = (holes_base, pts_base)

    # 2. Close-Packed 2x2 (10 nm vertical bridge at y = 0)
    holes_cp = [
        tdgl.Polygon(f"h_{i}", points=ellipse(0.012, 0.026, center=(x0, y0), points=28)).resample(28)
        for i, (x0, y0) in enumerate([
            (-0.042,  0.031), (0.042,  0.031),
            (-0.042, -0.031), (0.042, -0.031),
        ])
    ]
    pts_cp = [
        (0.000, 0.031, "1 (907 uA)"),
        (0.077, 0.031, "2 (911 uA)"),
        (0.020, 0.015, "3 (919 uA)"),
        (0.000, 0.000, "4 (923 uA)"),
    ]
    geoms["Close-Packed 2x2\n(10 nm Bridge)"] = (holes_cp, pts_cp)

    # 3. Width-Balanced Staggered (2-1-2)
    holes_stag = [
        tdgl.Polygon("h_ul", points=ellipse(0.010, 0.024, center=(-0.046,  0.038), points=28)).resample(28),
        tdgl.Polygon("h_ur", points=ellipse(0.010, 0.024, center=( 0.046,  0.038), points=28)).resample(28),
        tdgl.Polygon("h_c",  points=ellipse(0.018, 0.024, center=( 0.000,  0.000), points=32)).resample(32),
        tdgl.Polygon("h_ll", points=ellipse(0.010, 0.024, center=(-0.046, -0.038), points=28)).resample(28),
        tdgl.Polygon("h_lr", points=ellipse(0.010, 0.024, center=( 0.046, -0.038), points=28)).resample(28),
    ]
    pts_stag = [
        (0.000, 0.038, "1 (833 uA)"),
        (0.078, 0.038, "2 (826 uA)"),
        (0.028, 0.019, "3 (786 uA)"),
        (0.059, 0.000, "4 (830 uA)"),
    ]
    geoms["Width-Balanced Staggered\n(2-1-2 Lattice)"] = (holes_stag, pts_stag)

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
# 3. Solve Static Steady-State & Plot 2D Maps + 1D Profiles
# ==============================================================================
if __name__ == "__main__":
    geoms = get_geometries()
    solutions = {}

    for name, (holes, pts) in geoms.items():
        clean_name = name.replace("\n", " ")
        slug = clean_name.split("(")[0].strip().replace(" ", "_").replace("-", "_")
        print(f"Solving steady-state at I = {I_BIAS_STATIC} uA for: {clean_name}...")
        device = make_device(holes)
        options = tdgl.SolverOptions(
            solve_time=T_RELAX,
            dt_init=1e-3,
            dt_max=0.25,
            save_every=100,
            field_units="mT",
            current_units="uA",
            output_file=os.path.join(OUT_DIR, f"{slug}.h5"),
        )
        sol = tdgl.solve(
            device,
            options,
            terminal_currents=dict(source=I_BIAS_STATIC, drain=-I_BIAS_STATIC),
        )
        solutions[name] = (device, sol, pts)

        # Save individual full order parameter (|psi| + phase) & current plots using tdgl built-ins
        fig_psi, _ = sol.plot_order_parameter(figsize=(7, 3.5))
        fig_psi.suptitle(f"{clean_name} at I = {I_BIAS_STATIC:.0f} uA", fontsize=11, y=1.02)
        fig_psi.savefig(os.path.join(OUT_DIR, f"{slug}_order_parameter.png"), dpi=200, bbox_inches="tight")
        plt.close(fig_psi)

        fig_j, _ = sol.plot_currents(figsize=(7, 3.5), streamplot=True)
        fig_j.suptitle(f"{clean_name} Supercurrent J at I = {I_BIAS_STATIC:.0f} uA", fontsize=11, y=1.02)
        fig_j.savefig(os.path.join(OUT_DIR, f"{slug}_currents.png"), dpi=200, bbox_inches="tight")
        plt.close(fig_j)

    # ==========================================================================
    # 4. Build Side-by-Side 2x3 Thesis Figure (Top: |psi|, Bottom: |J| + Streamlines)
    # ==========================================================================
    fig, axes = plt.subplots(2, 3, figsize=(13, 9), sharex=True, sharey=True)

    for col, (name, (device, sol, pts)) in enumerate(solutions.items()):
        sol.solve_step = -1  # Final relaxed frame
        points = device.points
        triangles = device.triangles

        # Extract |psi| and sheet current density J [uA/um]
        psi_mag = np.abs(sol.tdgl_data.psi)
        J_vec = sol.interp_current_density(
            points, dataset="supercurrent", units="uA/um"
        )
        J_mag = np.linalg.norm(J_vec, axis=1)

        # --- Top Row: |psi(x, y)| ---
        ax_psi = axes[0, col]
        im_psi = ax_psi.tripcolor(
            points[:, 0], points[:, 1], psi_mag,
            triangles=triangles, shading="gouraud", cmap="viridis", vmin=0.55, vmax=1.0
        )
        for hole in device.holes:
            ax_psi.fill(hole.points[:, 0], hole.points[:, 1], color="white", edgecolor="black", lw=0.8)

        # Overlay tested photon impact points
        for x0, y0, lbl in pts:
            ax_psi.plot(x0, y0, "r*", markersize=9, markeredgecolor="k", markeredgewidth=0.5)
            ax_psi.annotate(
                lbl, (x0, y0), xytext=(4, 3), textcoords="offset points",
                fontsize=7.5, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.65, lw=0)
            )

        ax_psi.set_title(name, fontsize=11, fontweight="bold")
        ax_psi.set_xlim(-W / 2, W / 2)
        ax_psi.set_ylim(-0.12, 0.12)
        ax_psi.set_aspect("equal")
        if col == 0:
            ax_psi.set_ylabel("$y$ [$\\mu$m]", fontsize=11)

        # --- Bottom Row: |J(x, y)| + Streamlines ---
        ax_j = axes[1, col]
        im_j = ax_j.tripcolor(
            points[:, 0], points[:, 1], J_mag,
            triangles=triangles, shading="gouraud", cmap="inferno", vmin=0, vmax=np.percentile(J_mag, 99)
        )
        for hole in device.holes:
            ax_j.fill(hole.points[:, 0], hole.points[:, 1], color="white", edgecolor="black", lw=0.8)

        # Interpolate J onto regular grid for streamlines
        x_values = np.linspace(-W / 2, W / 2, 160)
        y_values = np.linspace(-0.12, 0.12, 160)
        x_grid, y_grid = np.meshgrid(x_values, y_values)
        Jx_interp = LinearNDInterpolator(points, J_vec[:, 0])(x_grid, y_grid)
        Jy_interp = LinearNDInterpolator(points, J_vec[:, 1])(x_grid, y_grid)

        # Mask points inside holes so streamlines don't cross holes
        grid_pts = np.column_stack([x_grid.ravel(), y_grid.ravel()])
        for hole in device.holes:
            in_hole = hole.contains_points(grid_pts).reshape(x_grid.shape)
            Jx_interp[in_hole] = np.nan
            Jy_interp[in_hole] = np.nan

        ax_j.streamplot(
            x_grid[0, :], y_grid[:, 0], Jx_interp, Jy_interp,
            color="cyan", density=1.1, linewidth=0.65, arrowsize=0.7
        )

        ax_j.set_xlim(-W / 2, W / 2)
        ax_j.set_ylim(-0.12, 0.12)
        ax_j.set_aspect("equal")
        ax_j.set_xlabel("$x$ [$\\mu$m]", fontsize=11)
        if col == 0:
            ax_j.set_ylabel("$y$ [$\\mu$m]", fontsize=11)

    # Add colorbars
    cbar_psi = fig.colorbar(im_psi, ax=axes[0, :], shrink=0.85, pad=0.02)
    cbar_psi.set_label("Order Parameter $|\\psi|$", fontsize=11)

    cbar_j = fig.colorbar(im_j, ax=axes[1, :], shrink=0.85, pad=0.02)
    cbar_j.set_label("Sheet Current Density $|\\mathbf{J}|$ [$\\mu$A/$\\mu$m]", fontsize=11)

    fig.suptitle(f"Comparison of Order Parameter $|\\psi|$ and Supercurrent $|\\mathbf{{J}}|$ at $I_{{\\mathrm{{bias}}}} = {I_BIAS_STATIC:.0f}\\,\\mu$A", fontsize=13, y=0.98)
    comb_path = os.path.join(OUT_DIR, "combined_psi_and_J_maps.png")
    fig.savefig(comb_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined 2x3 thesis map to: {comb_path}")

    # ==========================================================================
    # 5. Plot 1D Transverse Slices at y = 0 and y = y_row (Quantitative Proof)
    # ==========================================================================
    fig1d, (ax_y0, ax_yrow) = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    x_line = np.linspace(-W / 2 + 0.002, W / 2 - 0.002, 250)
    styles = ["k--", "b-", "r-"]

    for (name, (device, sol, _)), style in zip(solutions.items(), styles):
        clean_name = name.replace("\n", " ")
        points = device.points
        J_vec = sol.interp_current_density(
            points, dataset="supercurrent", units="uA/um"
        )
        J_mag = np.linalg.norm(J_vec, axis=1)
        J_interp = LinearNDInterpolator(points, J_mag)

        # Slice at y = 0 (Inter-row gap / Center staggered row)
        pts_y0 = np.column_stack([x_line, np.zeros_like(x_line)])
        j_y0 = J_interp(pts_y0)
        for hole in device.holes:
            j_y0[hole.contains_points(pts_y0)] = np.nan
        ax_y0.plot(x_line * 1e3, j_y0, style, lw=2.2, label=clean_name)

        # Slice at outer hole row (y = 45nm, 31nm, or 38nm)
        y_r = 0.045 if "Baseline" in name else (0.031 if "Close" in name else 0.038)
        pts_yr = np.column_stack([x_line, np.full_like(x_line, y_r)])
        j_yr = J_interp(pts_yr)
        for hole in device.holes:
            j_yr[hole.contains_points(pts_yr)] = np.nan
        ax_yrow.plot(x_line * 1e3, j_yr, style, lw=2.2, label=f"{clean_name} ($y={y_r*1e3:.0f}$ nm)")

    ax_y0.set_title("Transverse Current Profile at $y = 0$ nm", fontsize=11, fontweight="bold")
    ax_y0.set_xlabel("$x$ [nm]")
    ax_y0.set_ylabel("Sheet Current Density $|\\mathbf{J}(x, y)|$ [$\\mu$A/$\\mu$m]")
    ax_y0.grid(True, alpha=0.4)
    ax_y0.legend(fontsize=8.5)

    ax_yrow.set_title("Transverse Current Profile at Outer Hole Row ($y = y_{\\mathrm{row}}$)", fontsize=11, fontweight="bold")
    ax_yrow.set_xlabel("$x$ [nm]")
    ax_yrow.grid(True, alpha=0.4)
    ax_yrow.legend(fontsize=8.5)

    plt.tight_layout()
    slice_path = os.path.join(OUT_DIR, "transverse_J_profiles_y0_vs_yrow.png")
    fig1d.savefig(slice_path, dpi=200)
    plt.close(fig1d)
    print(f"Saved 1D transverse current profiles to: {slice_path}")