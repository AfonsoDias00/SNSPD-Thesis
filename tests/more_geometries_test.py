"""Compare switching and photon-detection thresholds for several hole geometries."""

from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import tdgl
from tdgl.geometry import box


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_more_geometries"

FIG_DIR = RESULTS_DIR / "figures"
for directory in (RESULTS_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 130

LAYER = tdgl.Layer(
    coherence_length=0.005,
    london_lambda=0.175,
    thickness=0.009,
)

WIDTH = 0.200
LENGTH = 0.800
MAX_EDGE = 0.01

T_RAMP = 20.0
T_PHOTON_ON = 30.0
T_PHOTON_OFF = 45.0
SOLVE_TIME = 95.0
R_HOTSPOT = 0.022
A_HOTSPOT = 2.5
VOLTAGE_THRESHOLD = 0.35


def ellipse_points(a_x: float, b_y: float, num_pts: int = 28) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, num_pts, endpoint=False)
    return np.column_stack((a_x * np.cos(angles), b_y * np.sin(angles)))


def stadium_points(
    half_straight: float,
    radius: float,
    num_pts: int = 12,
) -> np.ndarray:
    """Return a capsule with vertical straight sides and semicircular ends."""
    angles_top = np.linspace(0, np.pi, num_pts, endpoint=False)
    angles_bottom = np.linspace(np.pi, 2 * np.pi, num_pts, endpoint=False)
    top = np.column_stack(
        (radius * np.cos(angles_top), half_straight + radius * np.sin(angles_top))
    )
    bottom = np.column_stack(
        (radius * np.cos(angles_bottom), -half_straight + radius * np.sin(angles_bottom))
    )
    return np.vstack((top, bottom))


def superellipse_points(
    a_x: float,
    b_y: float,
    exponent: float = 3.0,
    num_pts: int = 40,
) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, num_pts, endpoint=False)
    cos_term = np.sign(np.cos(angles)) * np.abs(np.cos(angles)) ** (2.0 / exponent)
    sin_term = np.sign(np.sin(angles)) * np.abs(np.sin(angles)) ** (2.0 / exponent)
    return np.column_stack((a_x * cos_term, b_y * sin_term))


def translated_holes(
    points: np.ndarray,
    centers: list[tuple[float, float]],
) -> list[tdgl.Polygon]:
    return [
        tdgl.Polygon(f"h_{index}", points=points).translate(dx=x, dy=y)
        for index, (x, y) in enumerate(centers)
    ]


def geometry_holes(name: str) -> list[tdgl.Polygon]:
    if name == "baseline":
        points = ellipse_points(0.012, 0.027)
        centers = [(-0.042, -0.045), (0.042, -0.045), (-0.042, 0.045), (0.042, 0.045)]
    elif name == "stadium":
        points = stadium_points(0.015, 0.012)
        centers = [(-0.042, -0.045), (0.042, -0.045), (-0.042, 0.045), (0.042, 0.045)]
    elif name == "superellipse":
        points = superellipse_points(0.012, 0.027)
        centers = [(-0.042, -0.045), (0.042, -0.045), (-0.042, 0.045), (0.042, 0.045)]
    elif name == "staggered":
        points = ellipse_points(0.012, 0.022)
        centers = [
            (-0.042, -0.060),
            (0.042, -0.060),
            (0.0, 0.0),
            (-0.042, 0.060),
            (0.042, 0.060),
        ]
    else:
        raise ValueError(f"Unknown geometry: {name}")
    return translated_holes(points, centers)


def make_device(name: str) -> tdgl.Device:
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(240)
    source = (
        tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.10))
        .translate(dy=LENGTH / 2)
        .resample(60)
    )
    drain = source.scale(yfact=-1).set_name("drain")
    device = tdgl.Device(
        f"more_geometries_{name}",
        layer=LAYER,
        film=film,
        holes=geometry_holes(name),
        terminals=[source, drain],
        probe_points=[(0.0, 0.22), (0.0, -0.22)],
    )
    device.make_mesh(max_edge_length=MAX_EDGE, smooth=10)
    return device


def make_epsilon_func(
    photon_active: bool,
    r_photon: tuple[float, float],
):
    x0, y0 = r_photon

    def epsilon_func(r, *, t=0.0):
        r = np.asarray(r)
        x, y = r[..., 0], r[..., 1]
        epsilon = np.ones_like(x, dtype=float)
        if photon_active and T_PHOTON_ON <= t <= T_PHOTON_OFF:
            distance_sq = (x - x0) ** 2 + (y - y0) ** 2
            epsilon -= A_HOTSPOT * np.exp(
                -((distance_sq / R_HOTSPOT**2) ** 2)
            )
        return epsilon

    return epsilon_func


def safe_label(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label).strip("_")


def run_trial(
    device: tdgl.Device,
    bias_uA: float,
    photon_active: bool,
    r_photon: tuple[float, float],
    stem: str,
) -> tuple[tdgl.Solution, bool, float]:


    def terminal_currents(time: float) -> dict[str, float]:
        current = bias_uA * min(time / T_RAMP, 1.0)
        return {"source": current, "drain": -current}

    options = tdgl.SolverOptions(
        solve_time=SOLVE_TIME,
        dt_init=0.01,
        save_every=10,
        current_units="uA",
    )
    solution = tdgl.solve(
        device,
        options,
        applied_vector_potential=0.0,
        terminal_currents=terminal_currents,
        disorder_epsilon=make_epsilon_func(photon_active, r_photon),
    )
    if solution is None or solution.dynamics is None:
        raise RuntimeError(f"TDGL returned no dynamics for {stem}")

    times = solution.dynamics.time
    voltages = np.abs(solution.dynamics.mu[0] - solution.dynamics.mu[1])
    late_mask = times > 65.0
    if not np.any(late_mask):
        raise RuntimeError(f"No late-time samples available for {stem}")
    mean_voltage = float(np.mean(voltages[late_mask]))
    return solution, mean_voltage > VOLTAGE_THRESHOLD, mean_voltage


def find_threshold(
    device: tdgl.Device,
    geometry_name: str,
    label: str,
    photon_active: bool,
    r_photon: tuple[float, float],
    i_low: float,
    i_high: float,
    tolerance: float = 15.0,
) -> float:
    print(
        f"\n--- Bisection [{geometry_name}: {label}] "
        f"at r_ph={r_photon} (Nodes: {len(device.points)}) ---"
    )
    last_solution: tdgl.Solution | None = None
    last_switched_current = i_high

    while i_high - i_low > tolerance:
        i_mid = 0.5 * (i_low + i_high)
        stem = f"{geometry_name}_{label}_{i_mid:.0f}uA"
        solution, switched, mean_voltage = run_trial(
            device, i_mid, photon_active, r_photon, stem
        )
        status = "SWITCHED" if switched else "RECOVERED"
        print(
            f"  I = {i_mid:6.1f} uA -> <|V|> = "
            f"{mean_voltage:.4f} V0 [{status}]"
        )
        if switched:
            i_high = i_mid
            last_solution = solution
            last_switched_current = i_mid
        else:
            i_low = i_mid

    threshold = 0.5 * (i_low + i_high)
    if last_solution is not None:
        figure, axes = last_solution.plot_order_parameter()
        for axis in axes:
            axis.set_xlim(-WIDTH / 2, WIDTH / 2)
            axis.set_ylim(-0.20, 0.20)
            if photon_active:
                axis.plot(*r_photon, "r*", markersize=10)
        figure.suptitle(
            f"{geometry_name} {label}: switched at "
            f"{last_switched_current:.0f} uA"
        )
        figure.savefig(
            FIG_DIR / f"{safe_label(geometry_name)}_{safe_label(label)}_psi.png",
            bbox_inches="tight",
        )
        plt.close(figure)
    return threshold


def main() -> None:
    geometry_labels = {
        "baseline": "Baseline Ellipse (2x2)",
        "stadium": "Stadium Slots",
        "superellipse": "Superellipse (n=3)",
        "staggered": "Staggered Lattice (2-1-2)",
    }
    impact_site = (0.0, 0.0)
    results: dict[str, tuple[float, float, float]] = {}

    for geometry_name, display_name in geometry_labels.items():
        device = make_device(geometry_name)
        spontaneous = find_threshold(
            device,
            display_name,
            "spontaneous_Isw",
            False,
            impact_site,
            i_low=850.0,
            i_high=1180.0,
        )
        detection = find_threshold(
            device,
            display_name,
            "photon_center",
            True,
            impact_site,
            i_low=650.0,
            i_high=spontaneous,
        )
        results[geometry_name] = (
            spontaneous,
            detection,
            detection / spontaneous,
        )

    print("\n" + "=" * 78)
    print(f"{'Geometry':30s} | {'I_sw [uA]':>10s} | {'I_det [uA]':>10s} | {'I_det/I_sw':>11s}")
    print("-" * 78)
    for geometry_name, display_name in geometry_labels.items():
        spontaneous, detection, normalized = results[geometry_name]
        print(
            f"{display_name:30s} | {spontaneous:10.1f} | "
            f"{detection:10.1f} | {normalized:11.3f}"
        )
    print("=" * 78)


if __name__ == "__main__":
    main()
