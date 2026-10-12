"""Inspect the no-hole order parameter and voltage dynamics.

Each current is solved independently on a 10 nm mesh. The script saves the
TDGL result, an order-parameter plot, and CSV/summary plots containing the
minimum and mean final |psi|. It also plots the voltage dynamics for selected
currents. In pyTDGL 0.9.0, ``dynamics.voltage()`` is the normalized electric
potential difference between probe points, not a value in volts.

Example:
    .venv/bin/python order_parameter_current_sweep.py
    .venv/bin/python order_parameter_current_sweep.py \
        --currents-ua 0 100 200 300 400 500 600 700 750 800 900
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
RESULTS_DIR = BASE_DIR / "results_order_parameter_sweep_10nm"
H5_DIR = RESULTS_DIR / "h5"
FIGURE_DIR = RESULTS_DIR / "figures"

for directory in (H5_DIR, FIGURE_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 120

LAYER = tdgl.Layer(
    coherence_length=0.005,
    london_lambda=0.175,
    thickness=0.009,
)
WIDTH = 0.2
LENGTH = 2.0
MESH_NM = 10.0
MESH_SMOOTHING = 10
DEFAULT_CURRENTS_UA = tuple(float(value) for value in np.arange(1080.0, 1120.0, 10.0))
THERMALIZATION_TIME = 200.0
MEASUREMENT_TIME = 50.0
DT_INIT = 0.01
SAVE_EVERY = 10
VOLTAGE_THRESHOLD = 0.01
VOLTAGE_PLOT_CURRENTS_UA = (1080.0, 1100.0, 1110.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep final order-parameter magnitude on a 10 nm mesh."
    )
    parser.add_argument(
        "--currents-ua",
        nargs="+",
        type=float,
        default=DEFAULT_CURRENTS_UA,
        help="Independent currents to solve in uA (default: %(default)s).",
    )
    parser.add_argument(
        "--thermalization-time",
        type=float,
        default=THERMALIZATION_TIME,
        help="Unrecorded TDGL relaxation time (default: %(default)g).",
    )
    parser.add_argument(
        "--measurement-time",
        type=float,
        default=MEASUREMENT_TIME,
        help="Recorded TDGL time (default: %(default)g).",
    )
    parser.add_argument(
        "--analyze-existing",
        action="store_true",
        help="Analyze existing 10 nm HDF5 files without rerunning simulations.",
    )
    return parser.parse_args()


def safe_number(value: float) -> str:
    return f"{value:g}".replace("-", "minus").replace(".", "p")


def constant_terminal_currents(current_uA: float) -> dict[str, float]:
    return {"source": float(current_uA), "drain": -float(current_uA)}


def make_device() -> tdgl.Device:
    device = tdgl.Device(
        "no_hole_10nm",
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
        max_edge_length=MESH_NM / 1000.0,
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


def solve_one(
    device: tdgl.Device,
    current_uA: float,
    measurement_time: float,
    thermalization_time: float,
) -> tuple[dict[str, object], tdgl.Solution]:
    current_label = safe_number(current_uA)
    output_file = H5_DIR / f"no_hole_10nm_I{current_label}uA.h5"
    output_file.unlink(missing_ok=True)

    solution = tdgl.solve(
        device,
        make_options(output_file, measurement_time, thermalization_time),
        applied_vector_potential=0.0,
        terminal_currents=constant_terminal_currents(current_uA),
    )
    return summarize_solution(solution, current_uA)


def summarize_solution(
    solution: tdgl.Solution,
    current_uA: float,
) -> tuple[dict[str, object], tdgl.Solution]:
    dynamics = solution.dynamics
    if dynamics is None:
        raise RuntimeError(f"No dynamics data returned for I={current_uA:g} uA")

    psi = np.asarray(solution.tdgl_data.psi)
    if psi.ndim != 1:
        raise RuntimeError(
            f"Expected final psi snapshot to be one-dimensional, got shape {psi.shape}"
        )
    time = np.asarray(dynamics.time)
    voltage = np.asarray(dynamics.voltage())
    late_start = time[0] + 0.5 * (time[-1] - time[0])
    late_mask = time >= late_start
    late_voltage = np.abs(voltage[late_mask])
    late_mean = float(np.mean(late_voltage))
    late_std = float(np.std(late_voltage))
    relaxed_voltage = abs(float(dynamics.mean_voltage(tmin=late_start)))
    abs_psi = np.abs(psi)
    row: dict[str, object] = {
        "mesh_nm": MESH_NM,
        "current_uA": current_uA,
        "relaxed_mean_voltage": relaxed_voltage,
        "late_abs_voltage_mean": late_mean,
        "late_abs_voltage_std": late_std,
        "final_voltage": abs(float(voltage[-1])),
        "psi_shape": str(psi.shape),
        "min_abs_psi": float(np.min(abs_psi)),
        "mean_abs_psi": float(np.mean(abs_psi)),
        "switched": relaxed_voltage >= VOLTAGE_THRESHOLD,
    }
    return row, solution


def save_voltage_dynamics_plot(
    solutions: dict[float, tdgl.Solution],
) -> None:
    figure, axis = plt.subplots(figsize=(8, 5))
    for current_uA in VOLTAGE_PLOT_CURRENTS_UA:
        solution = solutions.get(current_uA)
        if solution is None or solution.dynamics is None:
            continue
        time = np.asarray(solution.dynamics.time)
        voltage = np.abs(np.asarray(solution.dynamics.voltage()))
        axis.plot(time, voltage, label=f"{current_uA:g} uA")

    axis.axhline(
        VOLTAGE_THRESHOLD,
        color="black",
        linestyle="--",
        label=f"threshold = {VOLTAGE_THRESHOLD:g}",
    )
    axis.set_xlabel("Dimensionless TDGL time")
    axis.set_ylabel(r"$|\Delta\mu|$ [normalized $V_0$]")
    axis.set_title("No-hole voltage dynamics: 10 nm mesh")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.savefig(
        FIGURE_DIR / "voltage_dynamics_600_650_700_750uA_10nm.png",
        bbox_inches="tight",
    )
    plt.close(figure)


def save_order_parameter_plot(
    solution: tdgl.Solution,
    current_uA: float,
) -> None:
    figure = solution.plot_order_parameter()[0]
    figure.suptitle(f"No hole, 10 nm mesh, I = {current_uA:g} uA")
    figure.savefig(
        FIGURE_DIR / f"no_hole_10nm_I{safe_number(current_uA)}uA_order_parameter.png",
        bbox_inches="tight",
    )
    plt.close(figure)


def save_summary(rows: list[dict[str, object]]) -> None:
    summary_file = RESULTS_DIR / "order_parameter_current_sweep_10nm.csv"
    with summary_file.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    currents = [float(row["current_uA"]) for row in rows]
    minimum = [float(row["min_abs_psi"]) for row in rows]
    mean = [float(row["mean_abs_psi"]) for row in rows]
    voltage = [float(row["relaxed_mean_voltage"]) for row in rows]

    figure, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(currents, minimum, "o-", label="min |psi|")
    axes[0].plot(currents, mean, "o-", label="mean |psi|")
    axes[0].set_ylabel("Final order parameter magnitude")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(currents, voltage, "o-", color="tab:red", label="|mean voltage|")
    axes[1].axhline(
        VOLTAGE_THRESHOLD,
        color="black",
        linestyle="--",
        label=f"threshold = {VOLTAGE_THRESHOLD:g}",
    )
    axes[1].set_xlabel("Applied current [uA]")
    axes[1].set_ylabel(r"$|\Delta\mu|$ [normalized $V_0$]")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    figure.suptitle("No-hole order parameter and voltage: 10 nm mesh")
    figure.savefig(
        FIGURE_DIR / "order_parameter_current_sweep_10nm.png",
        bbox_inches="tight",
    )
    plt.close(figure)

    print(f"\nSummary written to {summary_file}")


def load_existing_solution(current_uA: float) -> tdgl.Solution:
    path = H5_DIR / f"no_hole_10nm_I{safe_number(current_uA)}uA.h5"
    if not path.exists():
        raise FileNotFoundError(
            f"Existing result not found for I={current_uA:g} uA: {path}"
        )
    return tdgl.Solution.from_hdf5(str(path))


def main() -> None:
    args = parse_args()
    if not args.currents_ua:
        raise ValueError("At least one current is required.")
    if any(current < 0 for current in args.currents_ua):
        raise ValueError("Currents must be non-negative.")
    if args.thermalization_time < 0:
        raise ValueError("The thermalization time must be non-negative.")
    if args.measurement_time <= 0:
        raise ValueError("The measurement time must be positive.")

    rows: list[dict[str, object]] = []
    solutions: dict[float, tdgl.Solution] = {}
    device = None if args.analyze_existing else make_device()
    if device is not None:
        print("10 nm mesh statistics:")
        print(device.mesh_stats())

    for current_uA in args.currents_ua:
        if args.analyze_existing:
            print(f"Analyzing existing no_hole_10nm, I={current_uA:g} uA")
            solution = load_existing_solution(current_uA)
            row, _ = summarize_solution(solution, current_uA)
        else:
            print(f"Solving no_hole_10nm, I={current_uA:g} uA")
            row, solution = solve_one(
                device,
                current_uA,
                args.measurement_time,
                args.thermalization_time,
            )
        rows.append(row)
        solutions[float(current_uA)] = solution
        save_order_parameter_plot(solution, current_uA)
        print(
            f"  psi shape={row['psi_shape']}, "
            f"  min|psi|={row['min_abs_psi']:.5g}, "
            f"mean|psi|={row['mean_abs_psi']:.5g}, "
            f"switched={row['switched']}"
        )

    save_summary(rows)
    save_voltage_dynamics_plot(solutions)
    missing = sorted(set(VOLTAGE_PLOT_CURRENTS_UA) - set(solutions))
    if missing:
        print(
            "Voltage dynamics plot omitted currents not included in the sweep: "
            + ", ".join(f"{current:g} uA" for current in missing)
        )


if __name__ == "__main__":
    main()
