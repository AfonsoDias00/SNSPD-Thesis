"""Static-current sweeps for estimating the critical current.

The magnetic field is zero. Each current value is solved independently with
a constant terminal current; there is no current ramp and no continuation
from the previous current. The reported critical current is operational:
the first current whose relaxed mean voltage exceeds VOLTAGE_THRESHOLD.
"""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np
import tdgl
from tdgl.geometry import box, circle


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_critical_current"
H5_DIR = RESULTS_DIR / "h5"
MESH_DIR = RESULTS_DIR / "meshes"
FIGURE_DIR = RESULTS_DIR / "figures"

for directory in (H5_DIR, MESH_DIR, FIGURE_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 120


# ============================================================
# PHYSICAL PARAMETERS AND CONFIGURATIONS
# ============================================================

layer = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,      # 175 nm
    thickness=0.009,          # 9 nm
)

width = 0.2                 # 200 nm
length = 2.0                # 2 um
hole_diameter = 0.050       # 50 nm

CONFIGURATIONS = {
    "no_hole": [],
    "single_hole": [(0.0, 0.0)],
    "four_symmetric": [
        (-0.055, -0.100),
        (0.055, -0.100),
        (-0.055, 0.100),
        (0.055, 0.100),
    ],
}


# ============================================================
# SWEEP AND SOLVER SETTINGS
# ============================================================

CURRENTS_UA = {
    "no_hole": np.arange(404.0, 408.0, 2.0), #406 304
    "single_hole": np.arange(392.0, 396.0, 2.0), #394 300
    "four_symmetric": np.arange(338.0, 342.0, 2.0), #340 286
}

APPLIED_FIELD_MT = 0.0
SOLVE_TIME = 20.0
DT_INIT = 0.01
SAVE_EVERY = 10

MAX_EDGE_LENGTH = 0.005  # 5 nm; repeat final points with 0.0025 for convergence
MESH_SMOOTHING = 10

VOLTAGE_THRESHOLD = 0.01
RELAXATION_START_FRACTION = 0.5


def constant_terminal_currents(current_uA):
    """Return the static, equal-and-opposite terminal currents."""
    return {"source": float(current_uA), "drain": -float(current_uA)}


def make_device(name, hole_centers):
    film = tdgl.Polygon("film", points=box(width, length)).resample(200)

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
        probe_points=[(0.0, 0.3 * length), (0.0, -0.3 * length)],
    )
    device.make_mesh(
        max_edge_length=MAX_EDGE_LENGTH,
        smooth=MESH_SMOOTHING,
    )
    return device


def make_options(output_file):
    return tdgl.SolverOptions(
        solve_time=SOLVE_TIME,
        dt_init=DT_INIT,
        save_every=SAVE_EVERY,
        field_units="mT",
        current_units="uA",
        output_file=str(output_file),
    )


def current_name(current_uA):
    return f"{current_uA:g}".replace(".", "p")


def save_mesh(device, name):
    mesh_file = MESH_DIR / f"{name}_mesh.h5"
    mesh_file.unlink(missing_ok=True)
    device.to_hdf5(str(mesh_file))

    fig, ax = device.plot(mesh=True)
    ax.set_aspect("equal")
    ax.set_title(f"{name}: mesh")
    fig.savefig(FIGURE_DIR / f"{name}_mesh.png", bbox_inches="tight")
    plt.close(fig)


def solve_one(device, name, current_uA):
    """Solve one constant-current point and return diagnostics."""
    current_label = current_name(current_uA)
    output_file = H5_DIR / f"{name}_I{current_label}uA.h5"
    output_file.unlink(missing_ok=True)

    solution = tdgl.solve(
        device,
        make_options(output_file),
        applied_vector_potential=APPLIED_FIELD_MT,
        terminal_currents=constant_terminal_currents(current_uA),
    )

    dynamics = solution.dynamics
    if dynamics is None:
        raise RuntimeError(f"No dynamics data returned for {name}, I={current_uA:g} uA")

    t_start = RELAXATION_START_FRACTION * SOLVE_TIME
    relaxed_voltage = abs(float(dynamics.mean_voltage(tmin=t_start)))
    final_voltage = abs(float(dynamics.voltage()[-1]))

    psi = solution.tdgl_data.psi
    min_abs_psi = float(np.min(np.abs(psi)))
    final_abs_psi_mean = float(np.mean(np.abs(psi)))

    return {
        "configuration": name,
        "current_uA": float(current_uA),
        "relaxed_mean_voltage": relaxed_voltage,
        "final_voltage": final_voltage,
        "min_abs_psi": min_abs_psi,
        "mean_abs_psi": final_abs_psi_mean,
        "switched": relaxed_voltage >= VOLTAGE_THRESHOLD,
    }


def save_iv_plot(rows):
    fig, ax = plt.subplots()
    for name in CONFIGURATIONS:
        selected = [row for row in rows if row["configuration"] == name]
        currents = [row["current_uA"] for row in selected]
        voltages = [row["relaxed_mean_voltage"] for row in selected]
        ax.plot(currents, voltages, "o-", label=name)
    ax.axhline(
        VOLTAGE_THRESHOLD,
        color="black",
        linestyle="--",
        label=f"threshold = {VOLTAGE_THRESHOLD:g}",
    )
    ax.set_xlabel("Applied current [uA]")
    ax.set_ylabel("Relaxed mean voltage [V0]")
    ax.set_title("Static-current sweep at B = 0")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(FIGURE_DIR / "critical_current_sweep.png", bbox_inches="tight")
    plt.close(fig)


def save_summary(rows):
    summary_file = RESULTS_DIR / "critical_current_summary.csv"
    fieldnames = list(rows[0])
    with summary_file.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print("\nOperational critical-current estimates:")
    for name in CONFIGURATIONS:
        switched_rows = [
            row for row in rows
            if row["configuration"] == name and row["switched"]
        ]
        estimate = (
            switched_rows[0]["current_uA"] if switched_rows else "not reached"
        )
        print(f"  {name}: {estimate} uA")
    print(f"\nSummary written to {summary_file}")


def main():
    rows = []

    devices = {}
    for name, hole_centers in CONFIGURATIONS.items():
        device = make_device(name, hole_centers)
        devices[name] = device
        print(f"\n{name} mesh statistics:")
        print(device.mesh_stats())
        save_mesh(device, name)

    for name, device in devices.items():
        for current_uA in CURRENTS_UA[name]:
            print(f"Solving {name}, I={current_uA:g} uA")
            rows.append(solve_one(device, name, current_uA))

    save_summary(rows)
    save_iv_plot(rows)


if __name__ == "__main__":
    main()
