"""Two-hole TDGL simulations
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tdgl
from tdgl.geometry import box, circle


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_two_holes"
H5_DIR = RESULTS_DIR / "h5"
MESH_DIR = RESULTS_DIR / "meshes"
FIGURE_DIR = RESULTS_DIR / "figures"

for directory in (H5_DIR, MESH_DIR, FIGURE_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 120


# ============================================================
# PHYSICAL AND DEVICE PARAMETERS
# ============================================================

layer = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,     # 175 nm
    thickness=0.009,         # 9 nm
)

width = 0.2                # 200 nm
length = 1.0               # 1 um
hole_diameter = 0.050      # 50 nm

# The current is along y, so a transverse pair is separated in x.
# These positions leave 20 nm from each hole to the side boundary and
# 60 nm between the hole edges.
HOLE_CONFIGURATIONS = {
    "no_hole": [],
    "pair_symmetric": [(-0.055, 0.0), (0.055, 0.0)],
    "pair_shifted": [(-0.055, -0.015), (0.055, 0.015)],
}


# ============================================================
# CURRENT AND FIELD SETTINGS
# ============================================================

I_MAX = 500.0             # uA
T_START = 0.0
T_END_RAMP = 2.0
SOLVE_TIME = 50.0

# APPLIED_FIELDS_MT = (0.0, 20.0, 50.0, 100.0, 250.0)
APPLIED_FIELDS_MT = (-100.0,)

MAX_EDGE_LENGTH = 0.005  
MESH_SMOOTHING = 10


def terminal_currents(time):
    """Ramp equal and opposite currents through the two terminals."""
    if time < T_START:
        current = 0.0
    elif time < T_END_RAMP:
        current = I_MAX * (time - T_START) / (T_END_RAMP - T_START)
    else:
        current = I_MAX

    return {"source": current, "drain": -current}


def make_device(name, hole_centers):
    """Construct and mesh one device for a specified hole configuration."""
    film = (
        tdgl.Polygon("film", points=box(width, length))
        .resample(200)
    )

    holes = []
    for index, (hole_x, hole_y) in enumerate(hole_centers):
        hole = (
            tdgl.Polygon(
                f"hole_{index}",
                points=circle(hole_diameter / 2),
            )
            .translate(dx=hole_x, dy=hole_y)
            .resample(128)
        )
        holes.append(hole)

    source = (
        tdgl.Polygon("source", points=box(1.1 * width, 0.2))
        .resample(20)
        .translate(dy=length / 2)
    )
    drain = source.scale(yfact=-1).set_name("drain")

    device = tdgl.Device(
        name,
        layer=layer,
        film=film,
        holes=holes,
        terminals=[source, drain],
        probe_points=[(0.0, 0.4 * length), (0.0, -0.4 * length)],
    )
    device.make_mesh(
        max_edge_length=MAX_EDGE_LENGTH,
        smooth=MESH_SMOOTHING,
    )
    return device


def make_options(output_file):
    return tdgl.SolverOptions(
        solve_time=SOLVE_TIME,
        dt_init=0.01,
        save_every=10,
        field_units="mT",
        current_units="uA",
        output_file=str(output_file),
    )


def safe_field_name(field_mT):
    return f"{field_mT:g}".replace("-", "minus").replace(".", "p")


def run_simulation(name, hole_centers, field_mT):
    """Build, save, solve, and return one simulation."""
    device = make_device(name, hole_centers)
    field_name = safe_field_name(field_mT)
    stem = f"{name}_B{field_name}mT"

    mesh_file = MESH_DIR / f"{stem}_mesh.h5"
    output_file = H5_DIR / f"{stem}.h5"
    mesh_file.unlink(missing_ok=True)
    output_file.unlink(missing_ok=True)

    device.to_hdf5(str(mesh_file))

    fig, ax = device.plot(mesh=True)
    ax.set_aspect("equal")
    ax.set_title(f"{name}: mesh, hole diameter = {hole_diameter * 1000:g} nm")
    fig.savefig(FIGURE_DIR / f"{stem}_mesh.png", bbox_inches="tight")
    plt.close(fig)

    print("\n" + "=" * 70)
    print(f"Running {stem}")
    print("=" * 70)
    print(device.mesh_stats())

    solution = tdgl.solve(
        device,
        make_options(output_file),
        applied_vector_potential=field_mT,
        terminal_currents=terminal_currents,
    )

    return solution, stem


def save_solution_figures(solution, stem):
    """Save current, order-parameter, and vorticity diagnostics."""
    fig, ax = solution.plot_currents(streamplot=True)
    ax.set_title(f"{stem}: current density")
    fig.savefig(FIGURE_DIR / f"{stem}_currents.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = solution.plot_order_parameter()
    fig.suptitle(f"{stem}: order parameter")
    fig.savefig(FIGURE_DIR / f"{stem}_order_parameter.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = solution.plot_vorticity()
    fig.suptitle(f"{stem}: vorticity")
    fig.savefig(FIGURE_DIR / f"{stem}_vorticity.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = solution.dynamics.plot()
    fig.suptitle(f"{stem}: dynamics")
    fig.savefig(FIGURE_DIR / f"{stem}_dynamics.png", bbox_inches="tight")
    plt.close(fig)


def main():
    for field_mT in APPLIED_FIELDS_MT:
        for name, hole_centers in HOLE_CONFIGURATIONS.items():
            solution, stem = run_simulation(name, hole_centers, field_mT)
            save_solution_figures(solution, stem)


if __name__ == "__main__":
    main()
