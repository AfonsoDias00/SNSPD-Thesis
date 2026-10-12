"""No-hole mesh-convergence study for the operational critical current.

All lengths are in micrometres, so 5 nm corresponds to max_edge_length=0.005
and 100 nm corresponds to max_edge_length=0.1.  The script solves the same
static-current linear sweep on every mesh and writes enough diagnostics to
separate mesh effects from a changed geometry or solver setup. The sweep
records every current in the requested interval instead of assuming a
monotonic threshold crossing.

Example:
    .venv/bin/python mesh_convergence_no_hole.py
    .venv/bin/python mesh_convergence_no_hole.py --meshes-nm 5 10 25 50 100
    .venv/bin/python mesh_convergence_no_hole.py --current-lower-ua 390 \
        --current-upper-ua 420 --current-step-ua 1 \
        --thermalization-time 50 --measurement-time 50 \
        --steady-relative-tolerance 0.01
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tdgl
from tdgl.geometry import box


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_mesh_convergence_no_hole"
H5_DIR = RESULTS_DIR / "h5"
MESH_DIR = RESULTS_DIR / "meshes"
FIGURE_DIR = RESULTS_DIR / "figures"

for directory in (H5_DIR, MESH_DIR, FIGURE_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 120


# Physical parameters, matching critical_current_sweep.py.
LAYER = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,      # 175 nm
    thickness=0.009,          # 9 nm
)
WIDTH = 0.2
LENGTH = 2.0

# Defaults bracket the transition seen in the existing no-hole sweep.
DEFAULT_MESHES_NM = (5.0, 10.0, 25.0, 50.0, 100.0)
DEFAULT_CURRENT_LOWER_UA = 746.0
DEFAULT_CURRENT_UPPER_UA = 767.0
DEFAULT_CURRENT_STEP_UA = 2.0
MESH_SMOOTHING = 10
THERMALIZATION_TIME = 50.0
MEASUREMENT_TIME = 50.0
DT_INIT = 0.01
SAVE_EVERY = 10
VOLTAGE_THRESHOLD = 0.01
STEADY_RELATIVE_TOLERANCE = 0.01
RELATIVE_VOLTAGE_EPSILON = 1e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a no-hole TDGL mesh-convergence study."
    )
    parser.add_argument(
        "--meshes-nm",
        nargs="+",
        type=float,
        default=DEFAULT_MESHES_NM,
        help="Maximum triangle edge lengths in nm (default: %(default)s).",
    )
    parser.add_argument(
        "--current-lower-ua",
        type=float,
        default=DEFAULT_CURRENT_LOWER_UA,
        help="Known non-switched lower bracket in uA (default: %(default)g).",
    )
    parser.add_argument(
        "--current-upper-ua",
        type=float,
        default=DEFAULT_CURRENT_UPPER_UA,
        help="Known switched upper bracket in uA (default: %(default)g).",
    )
    parser.add_argument(
        "--current-step-ua",
        type=float,
        default=DEFAULT_CURRENT_STEP_UA,
        help="Current spacing for the linear sweep in uA (default: %(default)g).",
    )
    parser.add_argument(
        "--thermalization-time",
        type=float,
        default=THERMALIZATION_TIME,
        help="Unrecorded TDGL relaxation time before measurements (default: %(default)g).",
    )
    parser.add_argument(
        "--measurement-time",
        type=float,
        default=MEASUREMENT_TIME,
        help="Recorded time used for the late steady-state windows (default: %(default)g).",
    )
    parser.add_argument(
        "--steady-relative-tolerance",
        type=float,
        default=STEADY_RELATIVE_TOLERANCE,
        help=(
            "Maximum relative difference between the two late-window means "
            "for steady state (default: %(default)g)."
        ),
    )
    return parser.parse_args()


def safe_number(value: float) -> str:
    return f"{value:g}".replace("-", "minus").replace(".", "p")


def constant_terminal_currents(current_uA: float) -> dict[str, float]:
    return {"source": float(current_uA), "drain": -float(current_uA)}


def make_device(mesh_nm: float) -> tdgl.Device:
    """Build the same no-hole geometry for one requested mesh size."""
    device = tdgl.Device(
        f"no_hole_{safe_number(mesh_nm)}nm",
        layer=LAYER,
        film=tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(200),
        holes=[],
        terminals=[
            tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.2))
            .resample(20)
            .translate(dy=LENGTH / 2),
            tdgl.Polygon("drain", points=box(1.1 * WIDTH, 0.2))
            .resample(20)
            .translate(dy=-LENGTH / 2),
        ],
        probe_points=[(0.0, 0.3 * LENGTH), (0.0, -0.3 * LENGTH)],
    )
    device.make_mesh(
        max_edge_length=mesh_nm / 1000.0,
        smooth=MESH_SMOOTHING,
    )
    return device


def make_options(
    output_file: Path,
    measurement_time: float,
    thermalization_time: float,
) -> tdgl.SolverOptions:
    return tdgl.SolverOptions(
        solve_time=measurement_time,
        skip_time=thermalization_time,
        dt_init=DT_INIT,
        save_every=SAVE_EVERY,
        field_units="mT",
        current_units="uA",
        output_file=str(output_file),
    )


def save_mesh_plot(device: tdgl.Device, mesh_nm: float) -> None:
    label = safe_number(mesh_nm)
    mesh_file = MESH_DIR / f"no_hole_{label}nm.h5"
    mesh_file.unlink(missing_ok=True)
    device.to_hdf5(str(mesh_file))

    fig, ax = device.plot(mesh=True)
    ax.set_aspect("equal")
    ax.set_title(f"No hole, requested maximum edge = {mesh_nm:g} nm")
    fig.savefig(FIGURE_DIR / f"no_hole_{label}nm_mesh.png", bbox_inches="tight")
    plt.close(fig)


def solve_one(
    device: tdgl.Device,
    mesh_nm: float,
    current_uA: float,
    measurement_time: float,
    thermalization_time: float,
    steady_relative_tolerance: float,
    iteration: int,
) -> tuple[dict[str, object], tdgl.Solution]:
    mesh_label = safe_number(mesh_nm)
    current_label = safe_number(current_uA)
    output_file = H5_DIR / f"no_hole_{mesh_label}nm_I{current_label}uA.h5"
    output_file.unlink(missing_ok=True)

    solution = tdgl.solve(
        device,
        make_options(output_file, measurement_time, thermalization_time),
        applied_vector_potential=0.0,
        terminal_currents=constant_terminal_currents(current_uA),
    )
    dynamics = solution.dynamics
    if dynamics is None:
        raise RuntimeError(
            f"No dynamics data returned for mesh={mesh_nm:g} nm, "
            f"I={current_uA:g} uA"
        )

    time = np.asarray(dynamics.time)
    voltage = np.asarray(dynamics.voltage())
    dt = np.asarray(dynamics.dt)
    if len(time) < 4:
        raise RuntimeError(
            f"Too few recorded voltage samples for mesh={mesh_nm:g} nm, "
            f"I={current_uA:g} uA; increase measurement time."
        )

    # Use the actual recorded time range. Adaptive stepping and save_every
    # mean that requested window boundaries need not be sample times.
    measurement_start = float(time[0])
    measurement_end = float(time[-1])
    measurement_range = measurement_end - measurement_start
    late_start = measurement_start + 0.5 * measurement_range
    late_midpoint = measurement_start + 0.75 * measurement_range
    late_first_mask = (time >= late_start) & (time < late_midpoint)
    late_second_mask = time >= late_midpoint

    def weighted_mean(mask: np.ndarray, window_name: str) -> float:
        weights = dt[mask]
        if not np.any(mask) or np.sum(weights) <= 0:
            raise RuntimeError(
                f"No positive time-step weights in {window_name} window for "
                f"mesh={mesh_nm:g} nm, I={current_uA:g} uA; "
                "increase measurement time."
            )
        return float(np.average(voltage[mask], weights=weights))

    late_first_window_mean = weighted_mean(late_first_mask, "late first")
    late_second_window_mean = weighted_mean(late_second_mask, "late second")
    window_difference = abs(late_second_window_mean - late_first_window_mean)
    relative_window_difference = window_difference / max(
        abs(late_second_window_mean),
        RELATIVE_VOLTAGE_EPSILON,
    )
    relaxed_voltage = abs(late_second_window_mean)
    final_voltage = abs(float(dynamics.voltage()[-1]))
    psi = solution.tdgl_data.psi
    row: dict[str, object] = {
        "mesh_nm": mesh_nm,
        "current_uA": current_uA,
        "sweep_iteration": iteration,
        "relaxed_mean_voltage": relaxed_voltage,
        "late_first_window_mean_voltage": late_first_window_mean,
        "late_second_window_mean_voltage": late_second_window_mean,
        "window_mean_difference": window_difference,
        "relative_window_mean_difference": relative_window_difference,
        "measurement_start": measurement_start,
        "late_window_start": late_start,
        "late_window_midpoint": late_midpoint,
        "measurement_end": measurement_end,
        "steady_relative_tolerance": steady_relative_tolerance,
        "final_voltage": final_voltage,
        "min_abs_psi": float(np.min(np.abs(psi))),
        "mean_abs_psi": float(np.mean(np.abs(psi))),
        "steady": relative_window_difference <= steady_relative_tolerance,
        # The sweep uses the late-window voltage. Steadiness is reported
        # separately because a large resistive voltage can still be relaxing.
        "switched": relaxed_voltage >= VOLTAGE_THRESHOLD,
        "steady_switched": (
            relaxed_voltage >= VOLTAGE_THRESHOLD
            and relative_window_difference <= steady_relative_tolerance
        ),
    }
    return row, solution


def sweep_critical_current(
    device: tdgl.Device,
    mesh_nm: float,
    lower_uA: float,
    upper_uA: float,
    step_uA: float,
    measurement_time: float,
    thermalization_time: float,
    steady_relative_tolerance: float,
) -> tuple[list[dict[str, object]], tdgl.Solution, float, float]:
    """Find the first threshold crossing in an increasing current sweep."""
    rows: list[dict[str, object]] = []
    solutions: dict[float, tdgl.Solution] = {}

    currents = list(np.arange(lower_uA, upper_uA, step_uA))
    if not currents or currents[-1] != upper_uA:
        currents.append(upper_uA)

    for iteration, current_uA in enumerate(currents):
        row, solution = solve_one(
            device,
            mesh_nm,
            current_uA,
            measurement_time,
            thermalization_time,
            steady_relative_tolerance,
            iteration,
        )
        rows.append(row)
        solutions[current_uA] = solution
        print(
            f"Solving no_hole, mesh={mesh_nm:g} nm, "
            f"I={current_uA:g} uA, switched={row['switched']}, "
            f"steady={row['steady']}"
        )
    if bool(rows[0]["switched"]):
        raise ValueError(
            f"Lower bracket {lower_uA:g} uA already switched for "
            f"{mesh_nm:g} nm mesh; choose a lower current."
        )
    switched_indices = [
        index for index, row in enumerate(rows) if bool(row["switched"])
    ]
    if not switched_indices:
        raise ValueError(
            f"Upper bracket {upper_uA:g} uA did not switch for "
            f"{mesh_nm:g} nm mesh; choose a higher current."
        )

    upper_index = switched_indices[0]
    lower_index = upper_index - 1
    bracket_lower = float(rows[lower_index]["current_uA"])
    bracket_upper = float(rows[upper_index]["current_uA"])

    for row in rows:
        row["bracket_lower_uA"] = bracket_lower
        row["bracket_upper_uA"] = bracket_upper
        row["bracket_width_uA"] = bracket_upper - bracket_lower
        row["critical_current_upper_bound_uA"] = bracket_upper

    return rows, solutions[bracket_upper], bracket_lower, bracket_upper


def save_transition_plots(
    solution: tdgl.Solution,
    mesh_nm: float,
    current_uA: float,
) -> None:
    """Save the threshold-relevant voltage trace and order parameter."""
    label = safe_number(mesh_nm)
    current_label = safe_number(current_uA)
    stem = f"no_hole_{label}nm_I{current_label}uA"

    dynamics = solution.dynamics
    if dynamics is None:
        raise RuntimeError(f"Cannot plot missing dynamics for {stem}")

    time = np.asarray(dynamics.time)
    voltage = np.asarray(dynamics.voltage())
    measurement_start = float(time[0])

    above = np.abs(voltage) >= VOLTAGE_THRESHOLD
    upward = np.flatnonzero(~above[:-1] & above[1:]) + 1
    downward = np.flatnonzero(above[:-1] & ~above[1:]) + 1

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(time, np.abs(voltage), color="tab:blue", label="|voltage|")
    ax.axhline(
        VOLTAGE_THRESHOLD,
        color="black",
        linestyle="--",
        label=f"threshold = {VOLTAGE_THRESHOLD:g}",
    )
    ax.axvline(
        measurement_start,
        color="tab:orange",
        linestyle=":",
        label=f"after thermalization (t={measurement_start:g})",
    )
    if upward.size:
        ax.scatter(
            time[upward],
            np.abs(voltage[upward]),
            color="tab:green",
            zorder=3,
            label="crossing upward",
        )
    if downward.size:
        ax.scatter(
            time[downward],
            np.abs(voltage[downward]),
            color="tab:red",
            zorder=3,
            label="crossing downward",
        )
    direction = (
        f"upward at t={time[upward[0]]:g}"
        if upward.size
        else f"downward at t={time[downward[0]]:g}"
        if downward.size
        else "no crossing"
    )
    ax.set_title(
        f"{stem}: voltage relaxation\n"
        f"first threshold crossing: {direction}"
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("|Voltage| [V0]")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.savefig(FIGURE_DIR / f"{stem}_voltage.png", bbox_inches="tight")
    plt.close(fig)
    fig.savefig(FIGURE_DIR / f"{stem}_voltage.png", bbox_inches="tight")
    plt.close(fig)

    fig = solution.plot_order_parameter()[0]
    fig.suptitle(f"{stem}: order parameter")
    fig.savefig(FIGURE_DIR / f"{stem}_order_parameter.png", bbox_inches="tight")
    plt.close(fig)


def save_summary(rows: list[dict[str, object]], mesh_stats: list[dict[str, object]]) -> None:
    summary_file = RESULTS_DIR / "no_hole_mesh_convergence.csv"
    fieldnames = list(rows[0])
    with summary_file.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    stats_file = RESULTS_DIR / "no_hole_mesh_stats.csv"
    stats_fields = ["mesh_nm", "requested_max_edge_um", "mesh_stats"]
    with stats_file.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=stats_fields)
        writer.writeheader()
        writer.writerows(mesh_stats)

    print("\nOperational critical-current brackets:")
    for mesh_nm in sorted({float(row["mesh_nm"]) for row in rows}):
        mesh_rows = [row for row in rows if row["mesh_nm"] == mesh_nm]
        upper = mesh_rows[0]["critical_current_upper_bound_uA"]
        lower = mesh_rows[0]["bracket_lower_uA"]
        print(f"  {mesh_nm:g} nm: {lower:g} < Ic <= {upper:g} uA")
    print(f"\nSummary written to {summary_file}")
    print(f"Mesh statistics written to {stats_file}")


def main() -> None:
    args = parse_args()
    if any(mesh_nm <= 0 for mesh_nm in args.meshes_nm):
        raise ValueError("All mesh sizes must be positive.")
    if args.current_lower_ua < 0 or args.current_upper_ua <= args.current_lower_ua:
        raise ValueError("The upper current must be greater than the lower current.")
    if args.current_step_ua <= 0:
        raise ValueError("The current step must be positive.")
    if args.thermalization_time < 0:
        raise ValueError("The thermalization time must be non-negative.")
    if args.measurement_time <= 0:
        raise ValueError("The measurement time must be positive.")
    if args.steady_relative_tolerance < 0:
        raise ValueError("The steady-state relative tolerance must be non-negative.")

    rows: list[dict[str, object]] = []
    mesh_stats: list[dict[str, object]] = []

    for mesh_nm in args.meshes_nm:
        device = make_device(mesh_nm)
        stats = device.mesh_stats()
        stats_text = getattr(stats, "data", str(stats))
        mesh_stats.append(
            {
                "mesh_nm": mesh_nm,
                "requested_max_edge_um": mesh_nm / 1000.0,
                "mesh_stats": stats_text,
            }
        )
        print(f"\n{mesh_nm:g} nm mesh statistics:")
        print(stats_text)
        save_mesh_plot(device, mesh_nm)

        mesh_rows, transition_solution, _, upper_uA = (
            sweep_critical_current(
                device,
                mesh_nm,
                args.current_lower_ua,
                args.current_upper_ua,
                args.current_step_ua,
                args.measurement_time,
                args.thermalization_time,
                args.steady_relative_tolerance,
            )
        )
        rows.extend(mesh_rows)
        save_transition_plots(
            transition_solution,
            mesh_nm,
            upper_uA,
        )

    save_summary(rows, mesh_stats)


if __name__ == "__main__":
    main()
