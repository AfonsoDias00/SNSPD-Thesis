"""Systematic 2-hole staggering sweep (Delta y) for current redistribution.

Sweeps the longitudinal offset of two 30 nm circular holes placed near the 
wire center (x = +/- 20 nm) to show how diagonal symmetry breaking relieves 
the Dayem-bridge bottleneck while keeping the central channel sensitized.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import tdgl
from tdgl.geometry import box, circle

# ============================================================
# DIRECTORY SETUP
# ============================================================
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_two_holes_stagger"
H5_DIR = RESULTS_DIR / "h5"
FIG_DIR = RESULTS_DIR / "figures"

for directory in (RESULTS_DIR, H5_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 120

# ============================================================
# PHYSICAL AND GEOMETRIC PARAMETERS
# ============================================================
LAYER = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,     # 175 nm
    thickness=0.009,         # 9 nm
)

WIDTH = 0.200          # 200 nm wire width
LENGTH = 2.000         # 2.0 um wire length
HOLE_DIAMETER = 0.030  # 30 nm fab-friendly diameter
X_OFFSET = 0.020       # +/- 20 nm from center (10 nm gap when dy = 0)

TEST_CURRENT_UA = 400.0  # Subcritical bias for static current mapping
SOLVE_TIME = 80.0
T_RAMP = 5.0


def make_two_hole_device(dy_um: float) -> tdgl.Device:
    """Builds a 200 nm wire with two 30 nm holes at (-x0, -dy) and (+x0, +dy)."""
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(500)

    hole_left = (
        tdgl.Polygon("hole_left", points=circle(HOLE_DIAMETER / 2))
        .translate(dx=-X_OFFSET, dy=-dy_um)
        .resample(120)
    )
    hole_right = (
        tdgl.Polygon("hole_right", points=circle(HOLE_DIAMETER / 2))
        .translate(dx=X_OFFSET, dy=dy_um)
        .resample(120)
    )

    source = (
        tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.2))
        .translate(dy=LENGTH / 2)
        .resample(120)
    )
    drain = source.scale(yfact=-1).set_name("drain")

    name = f"two_holes_dy_{dy_um*1000:.0f}nm"
    device = tdgl.Device(
        name,
        layer=LAYER,
        film=film,
        holes=[hole_left, hole_right],
        terminals=[source, drain],
        probe_points=[(0.0, 0.3 * LENGTH), (0.0, -0.3 * LENGTH)],
    )
    device.make_mesh(max_edge_length=0.010, smooth=10)
    return device



def terminal_currents(time: float) -> dict[str, float]:
    """Smoothly ramp the static test current to prevent vortex shock at t=0."""
    if time < T_RAMP:
        curr = TEST_CURRENT_UA * (time / T_RAMP)
    else:
        curr = TEST_CURRENT_UA
    return {"source": curr, "drain": -curr}


def analyze_staggered_currents(
    solution: tdgl.Solution, dy_um: float
) -> tuple[float, float, float]:
    """Extracts true superconducting current peaks while masking out hole voids."""
    pts = solution.device.points
    x = pts[:, 0]
    y = pts[:, 1]
    J_mag = np.linalg.norm(solution.current_density.magnitude, axis=1)

    line_x = np.linspace(-WIDTH / 2, WIDTH / 2, 500)
    hole_radius = HOLE_DIAMETER / 2
    margin = 0.002  # 2 nm safety margin outside the hole boundary

    # 1. Midplane slice (y = 0)
    y_mid = np.zeros_like(line_x)
    J_mid = griddata((x, y), J_mag, (line_x, y_mid), method="linear")
    
    # Mask out any point on the y=0 line that falls inside either hole
    in_left_mid = np.hypot(line_x - (-X_OFFSET), y_mid - (-dy_um)) <= (hole_radius + margin)
    in_right_mid = np.hypot(line_x - X_OFFSET, y_mid - dy_um) <= (hole_radius + margin)
    J_mid[in_left_mid | in_right_mid] = np.nan

    center_idx = np.argmin(np.abs(line_x))
    j_bridge_center = float(J_mid[center_idx])

    # 2. Upper hole plane slice (y = +dy_um)
    y_hp = np.full_like(line_x, dy_um)
    J_hole_plane = griddata((x, y), J_mag, (line_x, y_hp), method="linear")
    
    in_left_hp = np.hypot(line_x - (-X_OFFSET), y_hp - (-dy_um)) <= (hole_radius + margin)
    in_right_hp = np.hypot(line_x - X_OFFSET, y_hp - dy_um) <= (hole_radius + margin)
    J_hole_plane[in_left_hp | in_right_hp] = np.nan

    # Outer edge channels (|x| > X_OFFSET + hole_radius)
    outer_mask = np.abs(line_x) > (X_OFFSET + hole_radius + margin)
    j_edge_max = float(np.nanmax(J_hole_plane[outer_mask]))

    # True interior peak along the diagonal bridge (|x| < X_OFFSET) outside hole voids
    interior_mask = np.abs(line_x) < X_OFFSET
    j_interior_max = float(
        np.nanmax([np.nanmax(J_mid[interior_mask]), np.nanmax(J_hole_plane[interior_mask])])
    )

    return j_bridge_center, j_interior_max, j_edge_max


def save_zoomed_streamplot(solution: tdgl.Solution, dy_nm: float, gap_nm: float):
    """Saves a zoomed-in 2D current density plot around the active hole zone."""
    fig, ax = solution.plot_currents(streamplot=True)
    ax.set_xlim(-WIDTH / 2, WIDTH / 2)
    ax.set_ylim(-0.25, 0.25)  # Zoom into the +/- 250 nm active region
    ax.set_aspect("equal")
    ax.set_title(f"Δy = {dy_nm:.0f} nm (Diag. gap = {gap_nm:.1f} nm)")
    fig.savefig(FIG_DIR / f"streamplot_dy_{dy_nm:.0f}nm.png", bbox_inches="tight")
    plt.close(fig)


def main():
    # Sweep longitudinal offset dy from 0 nm (side-by-side) to 50 nm
    dy_values_nm = np.array([0.0, 10.0, 15.0, 20.0, 25.0, 35.0, 50.0])
    dy_values_um = dy_values_nm / 1000.0

    diag_gaps_nm = []
    j_centers = []
    j_interiors = []
    j_edges = []

    for dy_nm, dy_um in zip(dy_values_nm, dy_values_um):
        # Calculate edge-to-edge diagonal distance between the two holes
        center_dist_nm = np.hypot(2 * X_OFFSET * 1000.0, 2 * dy_nm)
        gap_nm = center_dist_nm - (HOLE_DIAMETER * 1000.0)
        diag_gaps_nm.append(gap_nm)

        print(f"\n--- Simulating dy = {dy_nm:.0f} nm (Diagonal bridge = {gap_nm:.1f} nm) ---")
        device = make_two_hole_device(dy_um)

        out_h5 = H5_DIR / f"stagger_dy_{dy_nm:.0f}nm.h5"
        out_h5.unlink(missing_ok=True)

        options = tdgl.SolverOptions(
            solve_time=SOLVE_TIME,
            save_every=80,
            current_units="uA",
            output_file=str(out_h5),
        )

        solution = tdgl.solve(
            device,
            options,
            applied_vector_potential=0.0,
            terminal_currents=terminal_currents,
        )

        j_cen, j_int, j_edg = analyze_staggered_currents(solution, dy_um)
        j_centers.append(j_cen)
        j_interiors.append(j_int)
        j_edges.append(j_edg)

        save_zoomed_streamplot(solution, dy_nm, gap_nm)
        print(
            f"J_origin(0,0): {j_cen:.1f} | J_interior_peak: {j_int:.1f} | J_edge_peak: {j_edg:.1f}"
        )

    # Plot summary: Current densities & Diagonal gap vs. Longitudinal Stagger
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    ax1.plot(dy_values_nm, j_interiors, "o-", color="tab:orange", linewidth=2, label="Interior Peak ($J_{interior}$)")
    ax1.plot(dy_values_nm, j_centers, "d--", color="tab:red", linewidth=1.8, label="Midpoint Bridge $J(0,0)$")
    ax1.plot(dy_values_nm, j_edges, "s-", color="tab:blue", linewidth=2, label="Outer Edge Peak ($J_{edge}$)")
    ax1.set_xlabel("Longitudinal Half-Stagger $\Delta y$ [nm]")
    ax1.set_ylabel("Sheet Current Density $|K|$ [$\\mu$A/$\\mu$m]")
    ax1.set_title("Current Steering vs. Longitudinal Stagger")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.plot(dy_values_nm, diag_gaps_nm, "^-", color="tab:green", linewidth=2)
    ax2.axhline(25.0, color="gray", linestyle="--", label="Safe EBL Fab Limit (~25 nm)")
    ax2.axhline(10.0, color="red", linestyle=":", label="Dayem Bridge Limit (2ξ = 10 nm)")
    ax2.set_xlabel("Longitudinal Half-Stagger $\Delta y$ [nm]")
    ax2.set_ylabel("Edge-to-Edge Diagonal Gap [nm]")
    ax2.set_title("Bridge Width Recovery via Staggering")
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    fig.tight_layout()
    fig.savefig(FIG_DIR / "stagger_optimization_summary.png", bbox_inches="tight")
    print(f"\nSaved summary plots and zoomed streamplots to {FIG_DIR}")


if __name__ == "__main__":
    main()