import os
import tempfile
import numpy as np
import matplotlib.pyplot as plt

import tdgl
from tdgl.geometry import box, circle

from shapely.geometry import Polygon as ShapelyPolygon
from shapely.ops import unary_union

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
H5_DIR = RESULTS_DIR / "h5"
MESH_DIR = RESULTS_DIR / "meshes"
FIGURE_DIR = RESULTS_DIR / "figures"

H5_DIR.mkdir(parents=True, exist_ok=True)
MESH_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

#plot sizes
plt.rcParams["figure.figsize"] = (8, 10)
plt.rcParams["figure.dpi"] = 120

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================

layer = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,       # 175 nm
    thickness=0.009,           # 9 nm
)

# Device geometry
width = 0.2          # 200 nm
length = 2.0         # 2 um

# Hole parameters
hole_diameter = 0.100       # 100 nm
hole_x = 0.0
hole_y = 0.0

# Current ramp parameters
I_max = 500.0
t_start = 0.0
t_end_ramp = 2.0
solve_time = 20.0


# ============================================================
# GEOMETRY FUNCTION
# ============================================================

def make_device(with_hole=True):

    # Outer superconducting strip
    film = tdgl.Polygon(
        "film",
        points=box(width, length)
    ).resample(600)

    holes = []

    if with_hole:
        hole = (
            tdgl.Polygon(
                "hole",
                points=circle(hole_diameter / 2),
            )
            .translate(dx=hole_x, dy=hole_y)
            .resample(200)
        )

        holes.append(hole)

    # Current terminals at top and bottom
    source = (
        tdgl.Polygon(
            "source",
            points=box(1.1 * width, 0.2),
        )
        .resample(200)
        .translate(dy=length / 2)
    )

    drain = (
        source
        .scale(yfact=-1)
        .set_name("drain")
    )

    # Voltage probes inside the device
    probe_points = [
        (0.0, 0.4 * length),
        (0.0, -0.4 * length),
    ]

    device = tdgl.Device(
        "with_hole" if with_hole else "no_hole",
        layer=layer,
        film=film,
        holes=holes,
        terminals=[source, drain],
        probe_points=probe_points,
    )

    # The mesh must be significantly finer than xi.
    device.make_mesh(
        max_edge_length=0.1,
        smooth=10,
    )
    
    return device


# ============================================================
# 4. CURRENT RAMP
# ============================================================

def current_ramp(t):
    """
    Ramp the current from 0 to I_max.

    pyTDGL expects terminal_currents(time) -> dictionary.
    The sum of all terminal currents must be zero.
    """

    if t < t_start:
        I = 0.0

    elif t < t_end_ramp:
        I = I_max * (t - t_start) / (t_end_ramp - t_start)

    else:
        I = I_max

    return {
        "source": I,
        "drain": -I,
    }


# ============================================================
# 5. SOLVER OPTIONS
# ============================================================

def make_options(output_file):
    return tdgl.SolverOptions(
        solve_time=20.0,
        dt_init=0.01,
        save_every=10,
        field_units="mT",
        current_units="uA",
        output_file=output_file,
    )


# ============================================================
# 6. RUN ONE DEVICE
# ============================================================

def run_simulation(with_hole):
    device = make_device(with_hole=with_hole)

    name = "hole" if with_hole else "no_hole"
    mesh_file = MESH_DIR / f"mesh_{name}.h5"
    if mesh_file.exists():
        mesh_file.unlink()
    device.to_hdf5(str(mesh_file))

    fig, ax = device.plot(mesh=True)
    ax.set_aspect("equal")
    ax.set_title("Graded mesh near hole")
    plt.savefig(FIGURE_DIR / "graded_mesh.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    
    output_file = H5_DIR / f"transport_{name}.h5"

    if output_file.exists():
        output_file.unlink()

    options = make_options(str(output_file))

    print("\n" + "=" * 70)
    print(f"Running: {name}")
    print("=" * 70)

    print(device.mesh_stats())

    return tdgl.solve(
        device,
        options,
        applied_vector_potential=20.0,
        terminal_currents=current_ramp,
    )


# ============================================================
# 7. RUN BOTH CASES
# ============================================================

solution_no_hole = run_simulation(with_hole=False)
solution_hole = run_simulation(with_hole=True)


# ============================================================
# 8. PLOT CURRENT DENSITY
# ============================================================

fig, axes = plt.subplots(
    1, 2,
    figsize=(16, 7),
    constrained_layout=True,
)

solution_no_hole.plot_currents(
    ax=axes[0],
    streamplot=True,
)

axes[0].set_title("No hole")

solution_hole.plot_currents(
    ax=axes[1],
    streamplot=True,
)

axes[1].set_title("100 nm hole")

fig.savefig(
    FIGURE_DIR / "current_density_comparison.png",
    dpi=150,
    bbox_inches="tight",
)
plt.close(fig)


# ============================================================
# 9. PLOT ORDER PARAMETER
# ============================================================

fig1, axes1 = solution_no_hole.plot_order_parameter()
fig1.suptitle("No hole: |psi|")

fig2, axes2 = solution_hole.plot_order_parameter()
fig2.suptitle("100 nm hole: |psi|")

fig1.savefig(FIGURE_DIR / "order_parameter_no_hole.png", dpi=150)
fig2.savefig(FIGURE_DIR / "order_parameter_hole.png", dpi=150)


# ============================================================
# 10. PLOT VORTICITY
# ============================================================

fig1, axes1 = solution_no_hole.plot_vorticity()
fig1.suptitle("No hole: vorticity")

fig2, axes2 = solution_hole.plot_vorticity()
fig2.suptitle("100 nm hole: vorticity")

fig1.savefig(FIGURE_DIR / "vorticity_no_hole.png", dpi=150)
fig2.savefig(FIGURE_DIR / "vorticity_hole.png", dpi=150)



# ============================================================
# 11. PLOT VOLTAGE / PHASE DYNAMICS
# ============================================================

fig1, axes1 = solution_no_hole.dynamics.plot()
fig1.suptitle("No hole: dynamics")
fig1.savefig(FIGURE_DIR / "dynamics_no_hole.png", dpi=150, bbox_inches="tight")

fig2, axes2 = solution_hole.dynamics.plot()
fig2.suptitle("100 nm hole: dynamics")
fig2.savefig(FIGURE_DIR / "dynamics_hole.png", dpi=150, bbox_inches="tight")