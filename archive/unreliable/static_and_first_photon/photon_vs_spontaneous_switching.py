"""Direct comparison of spontaneous switching (I_sw) vs. photon-induced 
switching (I_det) across solid and nanohole-engineered nanowires.

Computes the normalized photon detection threshold:
    i_det = I_det / I_sw
and the normalized single-photon operating margin:
    Delta_i_op = (I_sw - I_det) / I_sw = 1 - i_det
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
RESULTS_DIR = BASE_DIR / "results_photon_margin"

FIG_DIR = RESULTS_DIR / "figures"

for directory in (RESULTS_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 130

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================
LAYER = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,     # 175 nm
    thickness=0.009,         # 9 nm
)

WIDTH = 0.200   # 200 nm
LENGTH = 1.00  # 1.6 um (shorter length speeds up meshing while keeping active zone isolated)
R0 = 0.018      # 18 nm radius (36 nm diameter)

# Time settings (in units of tau_0)
T_RAMP = 25.0        # Smooth current ramp to prevent inductive shock
T_PHOTON_ON = 30.0   # Inject photon after current has fully stabilized
T_PHOTON_OFF = 45.0  # Hotspot duration = 12 tau_0
SOLVE_TIME = 110.0   # Evaluate steady-state voltage over t in [75, 110] tau_0

# Photon hotspot parameters (Gaussian suppression of epsilon)
R_HOTSPOT = 0.022    # 18 nm Gaussian radius (~36 nm effective hotspot diameter)
A_HOTSPOT = 2.5      # Depresses epsilon from 1.0 down to -0.6 at hotspot center

VOLTAGE_THRESHOLD = 0.35  # Threshold [V0] for sustained vortex-flow switching


# ============================================================
# GEOMETRY BUILDERS (Equal Single-Hole Area = pi * 18^2 nm^2)
# ============================================================
def make_ellipse_points(a_x: float, b_y: float, num_pts: int = 180) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, num_pts, endpoint=False)
    return np.column_stack([a_x * np.cos(angles), b_y * np.sin(angles)])


def build_geometries() -> dict[str, list[tdgl.Polygon]]:
    """Constructs the solid reference and 3 hole-engineered candidates."""
    geoms = {}

    # 0. Solid unperforated nanowire reference
    #geoms["0_solid_wire"] = []

    # 1. Symmetric 2x2 Circles (dx = +/- 42 nm, dy = +/- 45 nm -> pitch_y = 90 nm)
    sym_2x2 = [(-0.042, -0.045), (0.042, -0.045), (-0.042, 0.045), (0.042, 0.045)]
    #geoms["1_circles_2x2"] = [
    #    tdgl.Polygon(f"h_{i}", points=circle(R0)).translate(dx=x, dy=y).resample(160)
    #    for i, (x, y) in enumerate(sym_2x2)
    #]

    # 2. Closely-spaced 2x4 Circular Chain (Zhu et al. wake-shielding: pitch_y = 52 nm)
    y_chain = [-0.078, -0.026, 0.026, 0.078]
    chain_centers = [(x, y) for y in y_chain for x in (-0.042, 0.042)]
    geoms["2_circles_close_chain"] = [
        tdgl.Polygon(f"h_{i}", points=circle(R0)).translate(dx=x, dy=y).resample(160)
        for i, (x, y) in enumerate(chain_centers)
    ]

    # 3. Longitudinal Ellipses 2x2 (a_x = 12 nm, b_y = 27 nm -> lowest rim crowding)
    ell_pts = make_ellipse_points(a_x=0.012, b_y=0.027)
    geoms["3_ellipses_longitudinal"] = [
        tdgl.Polygon(f"h_{i}", points=ell_pts).translate(dx=x, dy=y).resample(160)
        for i, (x, y) in enumerate(sym_2x2)
    ]

    return geoms


def make_device(name: str, holes: list[tdgl.Polygon]) -> tdgl.Device:
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(550)
    source = (
        tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.16))
        .translate(dy=LENGTH / 2)
        .resample(160)
    )
    drain = source.scale(yfact=-1).set_name("drain")

    device = tdgl.Device(
        name,
        layer=LAYER,
        film=film,
        holes=holes if holes else None,
        terminals=[source, drain],
        probe_points=[(0.0, 0.25 * LENGTH), (0.0, -0.25 * LENGTH)],
    )
    device.make_mesh(max_edge_length=0.01, smooth=10)
    return device


# ============================================================
# CONTACT REINFORCEMENT & PHOTON HOTSPOT PARAMETER (epsilon)
# ============================================================
def make_epsilon_func(photon_active: bool, r_photon: tuple[float, float] = (0.0, 0.0)):
    """Creates the space- and time-dependent disorder parameter epsilon(r, t).
    
    1. Contact banks remain at the bulk value. TDGL requires epsilon <= 1, so
       they cannot be reinforced by assigning a value above 1.
    2. Photon hotspot: if photon_active is True, a Gaussian depression of radius
       R_HOTSPOT turns on between T_PHOTON_ON and T_PHOTON_OFF at r_photon.
    """
    x0, y0 = r_photon

    def epsilon_func(r, *, t=0):
        r = np.asarray(r)
        x, y = r[..., 0], r[..., 1]
        eps = np.ones_like(x)

        # Inject transient Gaussian photon hotspot in the active zone
        if photon_active and (T_PHOTON_ON <= t <= T_PHOTON_OFF):
            dist_sq = (x - x0) ** 2 + (y - y0) ** 2
            eps -= A_HOTSPOT * np.exp(-dist_sq / (2.0 * R_HOTSPOT ** 2))

        return eps

    return epsilon_func


# ============================================================
# SOLVER & BISECTION ENGINE
# ============================================================
def run_trial(
    device: tdgl.Device,
    bias_uA: float,
    photon_active: bool,
    r_photon: tuple[float, float],
    stem: str,
) -> tuple[tdgl.Solution, bool, float]:
    """Runs one simulation and checks if sustained vortex flow persists at late times."""

    def terminal_currents(time):
        if time < T_RAMP:
            curr = bias_uA * (time / T_RAMP)
        else:
            curr = bias_uA
        return {"source": curr, "drain": -curr}

    options = tdgl.SolverOptions(
        solve_time=SOLVE_TIME,
        dt_init=0.01,
        save_every=10,
        current_units="uA",
    )

    eps_fn = make_epsilon_func(photon_active=photon_active, r_photon=r_photon)

    solution = tdgl.solve(
        device,
        options,
        applied_vector_potential=0.0,
        terminal_currents=terminal_currents,
        disorder_epsilon=eps_fn,
    )

    # Check mean voltage in the post-hotspot late window (t > 75 tau_0)
    times = solution.dynamics.time
    voltages = np.abs(solution.dynamics.mu[0] - solution.dynamics.mu[1])
    late_mask = times > 75.0
    mean_v = float(np.mean(voltages[late_mask]))
    switched = mean_v > VOLTAGE_THRESHOLD

    return solution, switched, mean_v


def find_threshold_current(
    device: tdgl.Device,
    mode_label: str,
    photon_active: bool,
    r_photon: tuple[float, float],
    i_low: float,
    i_high: float,
    tolerance: float = 15.0,
) -> float:
    """Bisection search for either I_sw (no photon) or I_det (with photon)."""
    print(f"\n--- Bisection: {device.name} [{mode_label}] ---")
    last_switched_sol = None
    last_switched_i = i_high

    while (i_high - i_low) > tolerance:
        i_mid = 0.5 * (i_low + i_high)
        stem = f"{device.name}_{mode_label}_{i_mid:.0f}uA"
        sol, switched, mean_v = run_trial(device, i_mid, photon_active, r_photon, stem)

        state_str = "SWITCHED" if switched else "RECOVERED / SUPERCONDUCTING"
        print(f"  I = {i_mid:6.1f} uA -> Late <|V|> = {mean_v:.4f} V0 [{state_str}]")

        if switched:
            i_high = i_mid
            last_switched_sol = sol
            last_switched_i = i_mid
        else:
            i_low = i_mid

    threshold_i = 0.5 * (i_low + i_high)
    print(f"==> {device.name} [{mode_label}] Threshold = {threshold_i:.1f} uA")

    if last_switched_sol is not None:
        fig, axes = last_switched_sol.plot_order_parameter()
        for ax in axes:
            ax.set_xlim(-WIDTH / 2, WIDTH / 2)
            ax.set_ylim(-0.30, 0.30)
        fig.suptitle(f"{device.name} ({mode_label}) switched at {last_switched_i:.0f} uA")
        fig.savefig(
            FIG_DIR / f"{device.name}_{mode_label}_psi_{last_switched_i:.0f}uA.png",
            bbox_inches="tight",
        )
        plt.close(fig)

    return threshold_i


def main():
    geoms = build_geometries()
    results = {}

    for name, holes in geoms.items():
        device = make_device(name, holes)

        # Photon lands in the central channel at (x = 0, y = 0.045 um)
        # (For solid wire, (0, 0) is equivalent; for hole grids, (0, 0.045) is right between the top row holes)
        r_ph = (0.0, 0.0) if name == "0_solid_wire" else (0.0, 0.045)

        # 1. Measure Spontaneous Switching Current I_sw (No Photon)
        i_sw = find_threshold_current(
            device,
            mode_label="spontaneous_Isw",
            photon_active=False,
            r_photon=r_ph,
            i_low=500.0,
            i_high=1180.0,
            tolerance=15.0,
        )

        # 2. Measure Photon-Induced Detection Current I_det (With Photon Hotspot)
        i_det = find_threshold_current(
            device,
            mode_label="photon_Idet",
            photon_active=True,
            r_photon=r_ph,
            i_low=300.0,
            i_high=i_sw,
            tolerance=15.0,
        )

        i_norm = i_det / i_sw
        margin = 1.0 - i_norm
        results[name] = (i_sw, i_det, i_norm, margin)

    # Print final summary table
    print("\n" + "=" * 82)
    print(
        f"{'Geometry':26s} | {'I_sw [uA]':>10s} | {'I_det [uA]':>10s} | "
        f"{'i_det = I_det/I_sw':>18s} | {'Margin (1 - i_det)':>18s}"
    )
    print("-" * 82)
    for name, (i_sw, i_det, i_norm, margin) in results.items():
        print(
            f"{name:26s} | {i_sw:10.1f} | {i_det:10.1f} | "
            f"{i_norm:18.3f} | {margin*100:17.1f}%"
        )
    print("=" * 82)


if __name__ == "__main__":
    main()