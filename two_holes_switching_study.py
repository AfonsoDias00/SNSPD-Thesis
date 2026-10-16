"""Dynamic switching current (I_sw) and polarity asymmetry study for 2 holes.

Uses a fast ramped bisection to find the switching current I_sw of each 2-hole 
topology and saves zoomed-in order-parameter and vorticity plots at switching.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tdgl
from tdgl.geometry import box, circle

# ============================================================
# DIRECTORY SETUP
# ============================================================
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_two_holes_switching"

FIG_DIR = RESULTS_DIR / "figures"

for directory in (RESULTS_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 120

# ============================================================
# DEVICE PARAMETERS
# ============================================================
LAYER = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,     # 175 nm
    thickness=0.009,         # 9 nm
)

WIDTH = 0.200
LENGTH = 2.000
HOLE_DIAMETER = 0.030  # 30 nm

T_START = 0.0
T_END_RAMP = 25.0
SOLVE_TIME = 120.0
VOLTAGE_THRESHOLD = 0.05  # Mean late-time voltage threshold [V0] for switching

# Coordinates in micrometres: [(x1, y1), (x2, y2)]
TWO_HOLE_CONFIGS = {
    # 0. Solid control wire (no holes) to measure the true mesh baseline I_sw
    "no_hole_reference": [],
    # 1. Side-by-side: 10 nm central gap (Dayem bridge bottleneck)
    "aligned_symmetric": [
        (-0.020, 0.000),
        (0.020, 0.000),
    ],
    # 2. Staggered along y: ~34 nm diagonal bridge, no simultaneous x-choke
    "diagonal_staggered": [
        (-0.020, -0.025),
        (0.020, 0.025),
    ],
    # 3. Shifted toward left edge + staggered: breaks x and y reflection
    "transverse_asymmetric": [
        (-0.045, -0.020),
        (0.015, 0.020),
    ],
}


def make_device(name: str, hole_centers: list[tuple[float, float]]) -> tdgl.Device:
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(600)

    holes = []
    for idx, (hx, hy) in enumerate(hole_centers):
        hole = (
            tdgl.Polygon(f"hole_{idx}", points=circle(HOLE_DIAMETER / 2))
            .translate(dx=hx, dy=hy)
            .resample(200)
        )
        holes.append(hole)

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
        holes=holes if holes else None,
        terminals=[source, drain],
        probe_points=[(0.0, 0.3 * LENGTH), (0.0, -0.3 * LENGTH)],
    )
    device.make_mesh(max_edge_length=0.1, smooth=10)
    return device


def did_device_switch(solution: tdgl.Solution) -> tuple[bool, float]:
    """Checks mean voltage over the second half of the simulation (after ramp)."""
    times = solution.dynamics.time
    voltages = np.abs(solution.dynamics.mu[0] - solution.dynamics.mu[1])

    # Ignore the initial inductive ramp spike by masking t > 0.5 * SOLVE_TIME
    steady_mask = times > (0.5 * SOLVE_TIME)
    mean_v = float(np.mean(voltages[steady_mask]))
    return mean_v > VOLTAGE_THRESHOLD, mean_v


def run_single_bias(
    device: tdgl.Device, signed_current_uA: float, stem: str
) -> tuple[tdgl.Solution, bool, float]:
    """Runs a single time-dependent simulation at a given signed bias current."""


    def current_ramp(time):
        if time < T_START:
            val = 0.0
        elif time < T_END_RAMP:
            val = signed_current_uA * (time - T_START) / (T_END_RAMP - T_START)
        else:
            val = signed_current_uA
        return {"source": val, "drain": -val}

    options = tdgl.SolverOptions(
        solve_time=SOLVE_TIME,
        dt_init=0.01,
        save_every=100,
        current_units="uA",
    )

    solution = tdgl.solve(
        device,
        options,
        applied_vector_potential=0.0,
        terminal_currents=current_ramp,
    )
    if solution is None:
        raise RuntimeError(
            f"TDGL solve was cancelled before producing a solution for {stem!r}."
        )
    switched, mean_v = did_device_switch(solution)
    return solution, switched, mean_v


def _iter_axes(axes):
    """Yield plotting axes regardless of whether TDGL returns one or many."""
    return np.atleast_1d(axes).flat


def save_switching_diagnostics(solution: tdgl.Solution, stem: str):
    """Saves zoomed-in order parameter, vorticity, and voltage dynamics."""
    fig, axes = solution.plot_order_parameter()
    for ax in _iter_axes(axes):
        ax.set_xlim(-WIDTH / 2, WIDTH / 2)
        ax.set_ylim(-0.30, 0.30)
    fig.suptitle(f"{stem}: Order Parameter (Zoomed)")
    fig.savefig(FIG_DIR / f"{stem}_psi.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = solution.plot_vorticity()
    for ax in _iter_axes(axes):
        ax.set_xlim(-WIDTH / 2, WIDTH / 2)
        ax.set_ylim(-0.30, 0.30)
    fig.suptitle(f"{stem}: Vorticity (Zoomed)")
    fig.savefig(FIG_DIR / f"{stem}_vorticity.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = solution.dynamics.plot()
    fig.suptitle(f"{stem}: Voltage & Phase Dynamics")
    fig.savefig(FIG_DIR / f"{stem}_dynamics.png", bbox_inches="tight")
    plt.close(fig)


def find_switching_current(
    name: str,
    hole_centers: list[tuple[float, float]],
    polarity: int = 1,
    i_low: float = 450.0,
    i_high: float = 950.0,
    tolerance: float = 12.5,
) -> float:
    """Bisection search for I_sw under +1 or -1 current polarity."""
    device = make_device(name, hole_centers)
    pol_tag = "pos" if polarity > 0 else "neg"
    last_switched_sol = None
    last_switched_i = i_high

    print(f"\n=======================================================")
    print(f"Searching I_sw for '{name}' (Polarity: {pol_tag})")
    print(f"=======================================================")

    while (i_high - i_low) > tolerance:
        i_mid = 0.5 * (i_low + i_high)
        signed_i = polarity * i_mid
        stem = f"{name}_{pol_tag}_{i_mid:.1f}uA"

        sol, switched, mean_v = run_single_bias(device, signed_i, stem)
        status = "SWITCHED (Resistive)" if switched else "SUPERCONDUCTING"
        print(f"  I = {signed_i:+.1f} uA -> <|V|> = {mean_v:.4f} V0 [{status}]")

        if switched:
            i_high = i_mid
            last_switched_sol = sol
            last_switched_i = i_mid
        else:
            i_low = i_mid

    i_sw_est = 0.5 * (i_low + i_high)
    print(f"--> Estimated I_sw ({name}, {pol_tag}): {i_sw_est:.1f} uA")

    if last_switched_sol is not None:
        save_switching_diagnostics(
            last_switched_sol, f"{name}_{pol_tag}_switched_{last_switched_i:.0f}uA"
        )

    return i_sw_est


def main():
    results = {}
    # 0. Measure the baseline switching current of a solid wire (no holes)
    results["no_hole_reference"] = find_switching_current(
        "no_hole_reference", TWO_HOLE_CONFIGS["no_hole_reference"], polarity=+1
    )
    # 1. Compare symmetric aligned vs. diagonally staggered (+I)
    results["aligned_symmetric (+)"] = find_switching_current(
        "aligned_symmetric", TWO_HOLE_CONFIGS["aligned_symmetric"], polarity=+1
    )
    results["diagonal_staggered (+)"] = find_switching_current(
        "diagonal_staggered", TWO_HOLE_CONFIGS["diagonal_staggered"], polarity=+1
    )

    # 2. Check polarity dependence (+I vs -I) on the transversely asymmetric pair
    results["transverse_asymmetric (+)"] = find_switching_current(
        "transverse_asymmetric", TWO_HOLE_CONFIGS["transverse_asymmetric"], polarity=+1
    )
    results["transverse_asymmetric (-)"] = find_switching_current(
        "transverse_asymmetric", TWO_HOLE_CONFIGS["transverse_asymmetric"], polarity=-1
    )

    print("\n================ FINAL I_sw SUMMARY ================")
    for label, val in results.items():
        print(f"{label:28s}: {val:.1f} uA")


if __name__ == "__main__":
    main()