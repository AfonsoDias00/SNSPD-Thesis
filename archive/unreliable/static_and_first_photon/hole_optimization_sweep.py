"""Automated geometric optimization to achieve uniform current density.

This script sweeps the transverse offset of two 30 nm holes and plots 
the difference between edge current density and central current density.
The zero-crossing of this curve indicates perfect uniform sensitivity.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import tdgl
from tdgl.geometry import box, circle
from scipy.interpolate import griddata
from pathlib import Path

# ============================================================
# DIRECTORY SETUP
# ============================================================
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_optimization"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================
LAYER = tdgl.Layer(
    coherence_length=0.005,
    london_lambda=0.175,
    thickness=0.009,
)

WIDTH = 0.2
LENGTH = 2.0
HOLE_DIAMETER = 0.03  # 30 nm holes (easier for e-beam lithography)
TEST_CURRENT_UA = 400.0  # Safe static current well below Ic


def make_device(x_offset: float) -> tdgl.Device:
    """Builds a device with two holes at +/- x_offset."""
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(400)
    
    hole_left = (
        tdgl.Polygon("hole_left", points=circle(HOLE_DIAMETER / 2))
        .translate(dx=-x_offset)
        .resample(100)
    )
    hole_right = (
        tdgl.Polygon("hole_right", points=circle(HOLE_DIAMETER / 2))
        .translate(dx=x_offset)
        .resample(100)
    )
    
    source = (
        tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.2))
        .translate(dy=LENGTH / 2)
        .resample(100)
    )
    drain = source.scale(yfact=-1).set_name("drain")
    
    device = tdgl.Device(
        f"opt_{x_offset*1000:.0f}nm",
        layer=LAYER,
        film=film,
        holes=[hole_left, hole_right],
        terminals=[source, drain],
        probe_points=[(0.0, 0.3 * LENGTH), (0.0, -0.3 * LENGTH)]
    )
    
    device.make_mesh(max_edge_length=0.1, smooth=10)
    return device


def get_j_profile(solution: tdgl.Solution, x_offset: float) -> tuple[float, float]:
    """Extracts J_y along the x-axis and returns the center and edge peaks."""
    x = solution.device.points[:, 0]
    y = solution.device.points[:, 1]
    
    # Extract the longitudinal supercurrent density
    J_y = np.abs(solution.current_density[:, 1])
    
    # Interpolate onto a high-resolution 1D line across the holes (y=0)
    line_x = np.linspace(-WIDTH / 2, WIDTH / 2, 500)
    line_y = np.zeros_like(line_x)
    J_y_interp = griddata((x, y), J_y, (line_x, line_y), method="cubic")
    
    # 1. Find the current density exactly in the center gap
    center_idx = np.argmin(np.abs(line_x))
    j_center = J_y_interp[center_idx]
    
    # 2. Find the maximum current density in the outer edge channels
    outer_mask = np.abs(line_x) > (x_offset + HOLE_DIAMETER / 2)
    if np.any(outer_mask):
        j_edge = np.nanmax(J_y_interp[outer_mask])
    else:
        j_edge = np.nan
        
    return j_center, j_edge


def main():
    # Sweep hole offset from 30 nm to 70 nm
    offsets_nm = np.linspace(30, 70, 9)
    offsets_um = offsets_nm / 1000.0
    
    j_centers = []
    j_edges = []
    delta_js = []
    
    for offset in offsets_um:
        print(f"\n--- Simulating offset: {offset*1000:.1f} nm ---")
        device = make_device(offset)
        
        # Short solve time: we only need the static Meissner distribution, no phase slips
        options = tdgl.SolverOptions(
            solve_time=50.0,
            save_every=50,
            output_file=str(RESULTS_DIR / f"sol_{offset*1000:.0f}nm.h5")
        )
        
        solution = tdgl.solve(
            device,
            options,
            applied_vector_potential=0.0,
            terminal_currents={"source": TEST_CURRENT_UA, "drain": -TEST_CURRENT_UA}
        )
        
        j_cen, j_edg = get_j_profile(solution, offset)
        delta_j = j_edg - j_cen
        
        j_centers.append(j_cen)
        j_edges.append(j_edg)
        delta_js.append(delta_j)
        
        print(f"J_center: {j_cen:.1f}, J_edge: {j_edg:.1f} | Delta: {delta_j:.1f}")
        
    # Plot 1: Absolute Current Densities
    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(offsets_nm, j_centers, 'o-', color="tab:orange", label="Center peak ($J_{center}$)")
    ax1.plot(offsets_nm, j_edges, 's-', color="tab:blue", label="Edge peak ($J_{edge}$)")
    ax1.set_xlabel("Hole Offset from Center [nm]")
    ax1.set_ylabel("Current Density $|J_y|$")
    ax1.set_title("Current Density Peaks vs. Transverse Hole Placement")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    fig1.savefig(RESULTS_DIR / "optimization_peaks.png", bbox_inches="tight")
    
    # Plot 2: The Optimization Curve (Zero-Crossing)
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(offsets_nm, delta_js, '^-', color='tab:red', markersize=8)
    ax2.axhline(0, color='black', linestyle='--', linewidth=2, label="Perfect Uniformity (ΔJ = 0)")
    ax2.set_xlabel("Hole Offset from Center [nm]")
    ax2.set_ylabel("ΔJ ($J_{edge} - J_{center}$)")
    ax2.set_title("Finding the Zero-Crossing")
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    fig2.savefig(RESULTS_DIR / "optimization_delta.png", bbox_inches="tight")
    
    print(f"\nOptimization plots saved to {RESULTS_DIR}")
    plt.show()

if __name__ == "__main__":
    main()