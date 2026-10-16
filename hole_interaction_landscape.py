from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from scipy.interpolate import griddata
import tdgl
from tdgl.geometry import box, circle

# ============================================================
# DIRECTORY SETUP
# ============================================================
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_interaction_landscape"
FIG_DIR = RESULTS_DIR / "figures"

for directory in (RESULTS_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 130

# ============================================================
# PHYSICAL & BASELINE PARAMETERS
# ============================================================
LAYER = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,     # 175 nm
    thickness=0.009,         # 9 nm
)

WIDTH = 0.200          # 200 nm wire width
LENGTH = 2.000         # 2.0 um wire length
I_BIAS_UA = 400.0      # Subcritical static current
K0_BG = I_BIAS_UA / WIDTH  

T_RAMP = 10.0          
SOLVE_TIME = 60.0      

#Hole radius
R0 = 0.018


def terminal_currents(time: float) -> dict[str, float]:
    """Smooth linear ramp to static subcritical bias."""
    if time < T_RAMP:
        curr = I_BIAS_UA * (time / T_RAMP)
    else:
        curr = I_BIAS_UA
    return {"source": curr, "drain": -curr}


# ============================================================
# EQUAL-AREA SHAPE GENERATORS
# ============================================================
def make_ellipse_points(a_x: float, b_y: float, num_pts: int = 200) -> np.ndarray:
    """Returns points for an ellipse with semi-axes a_x and b_y."""
    angles = np.linspace(0, 2 * np.pi, num_pts, endpoint=False)
    return np.column_stack([a_x * np.cos(angles), b_y * np.sin(angles)])


def make_rounded_triangle_points(
    target_area: float, corner_radius: float = 0.004, num_pts: int = 240
) -> np.ndarray:

    # Start with a unit equilateral-like triangle and scale to match target_area
    angles = np.deg2rad([90.0, 210.0, 330.0])
    verts = np.column_stack([np.cos(angles), np.sin(angles)])
    
    # Interpolate perimeter and apply Gaussian corner smoothing
    t = np.linspace(0, 3, num_pts, endpoint=False)
    pts = np.zeros((num_pts, 2))
    for i in range(3):
        v_start = verts[i]
        v_end = verts[(i + 1) % 3]
        mask = (t >= i) & (t < i + 1)
        frac = t[mask] - i
        pts[mask] = (1 - frac)[:, None] * v_start + frac[:, None] * v_end

    # Smooth sharp corners using a periodic moving-average kernel
    kernel_size = int(num_pts * 0.08)
    kernel = np.ones(kernel_size) / kernel_size
    x_smooth = np.convolve(np.pad(pts[:, 0], kernel_size, mode="wrap"), kernel, mode="valid")[:num_pts]
    y_smooth = np.convolve(np.pad(pts[:, 1], kernel_size, mode="wrap"), kernel, mode="valid")[:num_pts]
    pts_smooth = np.column_stack([x_smooth, y_smooth])

    # Compute polygon area via Shoelace formula and scale to exact target_area
    x_p, y_p = pts_smooth[:, 0], pts_smooth[:, 1]
    current_area = 0.5 * np.abs(np.dot(x_p, np.roll(y_p, 1)) - np.dot(y_p, np.roll(x_p, 1)))
    scale_factor = np.sqrt(target_area / current_area)
    return pts_smooth * scale_factor


def build_configurations() -> dict[str, list[tdgl.Polygon]]:
    """Constructs 5 equal-area 4-hole configurations."""
    target_area = np.pi * (R0 ** 2)

    # Standard symmetric 2x2 centers: dx = +/- 42 nm, dy = +/- 45 nm
    sym_centers = [
        (-0.042, -0.045),
        (0.042, -0.045),
        (-0.042, 0.045),
        (0.042, 0.045),
    ]

    # Sheared/staggered asymmetric centers (breaks x and y mirror symmetries)
    asym_centers = [
        (-0.048, -0.055),
        (0.036, -0.030),
        (-0.030, 0.030),
        (0.048, 0.055),
    ]

    configs = {}

    # 1. Symmetric Circles (R = 18 nm)
    configs["1_circles_symmetric"] = [
        tdgl.Polygon(f"hole_{i}", points=circle(R0)).translate(dx=cx, dy=cy).resample(180)
        for i, (cx, cy) in enumerate(sym_centers)
    ]

    # 2. Asymmetric Sheared Circles (R = 18 nm)
    configs["2_circles_asymmetric"] = [
        tdgl.Polygon(f"hole_{i}", points=circle(R0)).translate(dx=cx, dy=cy).resample(180)
        for i, (cx, cy) in enumerate(asym_centers)
    ]

    # 3. Longitudinal Ellipses (a_x = 12 nm, b_y = 27 nm -> Area = pi * 12 * 27 = pi * 18^2)
    ell_long_pts = make_ellipse_points(a_x=0.012, b_y=0.027)
    configs["3_ellipses_longitudinal"] = [
        tdgl.Polygon(f"hole_{i}", points=ell_long_pts).translate(dx=cx, dy=cy).resample(180)
        for i, (cx, cy) in enumerate(sym_centers)
    ]

    # 4. Transverse Ellipses (a_x = 24 nm, b_y = 13.5 nm -> Area = pi * 24 * 13.5 = pi * 18^2)
    ell_trans_pts = make_ellipse_points(a_x=0.024, b_y=0.0135)
    configs["4_ellipses_transverse"] = [
        tdgl.Polygon(f"hole_{i}", points=ell_trans_pts).translate(dx=cx, dy=cy).resample(180)
        for i, (cx, cy) in enumerate(sym_centers)
    ]

    # 5. Rounded Triangles (Exact same area, pointing along +y)
    tri_pts = make_rounded_triangle_points(target_area=target_area)
    configs["5_triangles"] = [
        tdgl.Polygon(f"hole_{i}", points=tri_pts).translate(dx=cx, dy=cy).resample(180)
        for i, (cx, cy) in enumerate(sym_centers)
    ]

    return configs


def make_device(name: str, holes: list[tdgl.Polygon]) -> tdgl.Device:
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(600)
    source = (
        tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.2))
        .translate(dy=LENGTH / 2)
        .resample(200)
    )
    drain = source.scale(yfact=-1).set_name("drain")

    device = tdgl.Device(
        name,
        layer=LAYER,
        film=film,
        holes=holes,
        terminals=[source, drain],
        probe_points=[(0.0, 0.3 * LENGTH), (0.0, -0.3 * LENGTH)],
    )
    device.make_mesh(max_edge_length=0.1, smooth=10)
    return device


# ============================================================
# INTERACTION FIELD EXTRACTION & PLOTTING
# ============================================================
def plot_interaction_landscape(solution: tdgl.Solution, name: str) -> dict[str, float]:
    """Interpolates and plots the 3 static interaction maps in the active zone."""
    device = solution.device
    pts = device.points
    x, y = pts[:, 0], pts[:, 1]

    # Extract raw quantities without pint unit warnings
    K_vec = solution.current_density.magnitude
    K_x = K_vec[:, 0]
    K_y = K_vec[:, 1]
    K_mag = np.linalg.norm(K_vec, axis=1)
    delta_K = K_mag - K0_BG  # Excess current above uniform background

    psi_mag_sq = np.abs(solution.tdgl_data.psi) ** 2
    energy_suppression = 1.0 - psi_mag_sq  # 1 - |Psi|^2

    # Build a high-resolution 2D grid over the zoomed interaction zone
    gx = np.linspace(-WIDTH / 2, WIDTH / 2, 240)
    gy = np.linspace(-0.16, 0.16, 380)
    GX, GY = np.meshgrid(gx, gy)
    grid_pts = np.column_stack([GX.ravel(), GY.ravel()])

    # Create a mask for points inside any of the 4 holes
    hole_mask = np.zeros(len(grid_pts), dtype=bool)
    for hole in device.holes:
        poly_path = MplPath(hole.points)
        hole_mask |= poly_path.contains_points(grid_pts, radius=0.001)
    hole_mask_2d = hole_mask.reshape(GX.shape)

    # Interpolate fields onto the 2D grid using linear interpolation
    dK_grid = griddata((x, y), delta_K, (GX, GY), method="linear")
    Kx_grid = griddata((x, y), K_x, (GX, GY), method="linear")
    Ky_grid = griddata((x, y), K_y, (GX, GY), method="linear")
    supp_grid = griddata((x, y), energy_suppression, (GX, GY), method="linear")

    for arr in (dK_grid, Kx_grid, Ky_grid, supp_grid):
        arr[hole_mask_2d] = np.nan

    # Create the 3-panel figure
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 5.2), sharey=True)
    fig.suptitle(f" {name} ($I_b = {I_BIAS_UA:.0f}\\,\\mu$A)", fontsize=13)

    # Panel 1: Current Perturbation (Delta |K| = |K| - K_0) + Streamlines
    vmax_dk = 2200.0
    im0 = axes[0].pcolormesh(
        GX, GY, dK_grid, cmap="RdBu_r", vmin=-vmax_dk, vmax=vmax_dk, shading="auto"
    )
    axes[0].streamplot(
        gx, gy, Kx_grid, Ky_grid, color="black", density=1.1, linewidth=0.6, arrowsize=0.7
    )
    axes[0].set_title("$\\Delta |K| = |K| - K_0$ [$\\mu$A/$\\mu$m]")
    fig.colorbar(im0, ax=axes[0], fraction=0.05, pad=0.03)

    # Panel 2: Transverse Steering Current (K_x)
    vmax_kx = 1200.0
    im1 = axes[1].pcolormesh(
        GX, GY, Kx_grid, cmap="coolwarm", vmin=-vmax_kx, vmax=vmax_kx, shading="auto"
    )
    axes[1].set_title("$K_x(x, y)$ [$\\mu$A/$\\mu$m]")
    fig.colorbar(im1, ax=axes[1], fraction=0.05, pad=0.03)

    # Panel 3: Condensation Energy Suppression (1 - |Psi|^2)
    vmax_supp = float(np.nanpercentile(supp_grid, 99.5))
    im2 = axes[2].pcolormesh(
        GX, GY, supp_grid, cmap="inferno", vmin=0.0, vmax=max(vmax_supp, 0.05), shading="auto"
    )
    axes[2].set_title("\n$1 - |\\Psi(x,y)|^2$")
    fig.colorbar(im2, ax=axes[2], fraction=0.05, pad=0.03)

    # Draw crisp hole outlines on all panels
    for ax in axes:
        for hole in device.holes:
            hp = np.vstack([hole.points, hole.points[0]])
            ax.fill(hp[:, 0], hp[:, 1], color="white", zorder=3)
            ax.plot(hp[:, 0], hp[:, 1], color="black", linewidth=1.0, zorder=4)
        ax.set_xlim(-WIDTH / 2, WIDTH / 2)
        ax.set_ylim(-0.16, 0.16)
        ax.set_aspect("equal")
        ax.set_xlabel("$x$ [$\\mu$m]")

    axes[0].set_ylabel("$y$ [$\\mu$m]")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"{name}_landscape.png", bbox_inches="tight")
    plt.close(fig)

    # Compute summary metrics for quantitative comparison
    interior_zone = (np.abs(GX) < 0.035) & (np.abs(GY) < 0.070)
    edge_zone = (np.abs(GX) > 0.065) & (np.abs(GY) < 0.070)

    metrics = {
        "max_K_interior": float(np.nanmax(dK_grid[interior_zone]) + K0_BG),
        "max_K_edge": float(np.nanmax(dK_grid[edge_zone]) + K0_BG),
        "max_abs_Kx": float(np.nanmax(np.abs(Kx_grid))),
        "net_center_Kx": float(np.nanmean(Kx_grid[np.abs(GX) < 0.020])),
        "max_psi_supp": float(np.nanmax(supp_grid)),
    }
    return metrics


def main():
    configs = build_configurations()
    summary = {}

    for name, holes in configs.items():
        print(f"\n=======================================================")
        print(f"Simulating static interaction landscape: {name}")
        print(f"=======================================================")
        device = make_device(name, holes)


        options = tdgl.SolverOptions(
            solve_time=SOLVE_TIME,
            dt_init=0.01,
            save_every=60,
            current_units="uA",
        )

        solution = tdgl.solve(
            device,
            options,
            applied_vector_potential=0.0,
            terminal_currents=terminal_currents,
        )

        summary[name] = plot_interaction_landscape(solution, name)

    # Print quantitative comparison table
    print("\n" + "=" * 88)
    print(
        f"{'Configuration':26s} | {'K_int_max':>10s} | {'K_edge_max':>10s} | "
        f"{'Max |Kx|':>9s} | {'<Kx>_cent':>10s} | {'Max(1-|Psi|^2)':>14s}"
    )
    print("-" * 88)
    for name, m in summary.items():
        print(
            f"{name:26s} | {m['max_K_interior']:10.1f} | {m['max_K_edge']:10.1f} | "
            f"{m['max_abs_Kx']:9.1f} | {m['net_center_Kx']:+10.1f} | {m['max_psi_supp']:14.4f}"
        )
    print("=" * 88)
    print(f"\nAll 3-panel interaction figures saved to: {FIG_DIR}")


if __name__ == "__main__":
    main()