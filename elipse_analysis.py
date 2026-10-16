"""Fast Spatial Photon Sensitivity Map i_det(x, y) for 3_ellipses_longitudinal.

Tests photon strikes across the central sub-channel, near the hole rim, in the 
outer edge sub-channel, and vertically between rows to verify full-width 
detection efficiency.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import tdgl
from tdgl.geometry import box

# ============================================================
# DIRECTORY SETUP
# ============================================================
BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_spatial_sensitivity"
H5_DIR = RESULTS_DIR / "h5"
FIG_DIR = RESULTS_DIR / "figures"

for directory in (RESULTS_DIR, H5_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams["figure.dpi"] = 130

# ============================================================
# PHYSICAL & FAST MESH PARAMETERS (~2100 nodes -> ~6 sec/step)
# ============================================================
LAYER = tdgl.Layer(
    coherence_length=0.005,  # 5 nm
    london_lambda=0.175,     # 175 nm
    thickness=0.009,         # 9 nm
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


def make_ellipse_points(a_x: float, b_y: float, num_pts: int = 28) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, num_pts, endpoint=False)
    return np.column_stack([a_x * np.cos(angles), b_y * np.sin(angles)])


def make_elliptical_device() -> tdgl.Device:
    """Builds the 3_ellipses_longitudinal device with a fast, uniform ~2100-node mesh."""
    film = tdgl.Polygon("film", points=box(WIDTH, LENGTH)).resample(240)

    sym_2x2 = [(-0.042, -0.045), (0.042, -0.045), (-0.042, 0.045), (0.042, 0.045)]
    ell_pts = make_ellipse_points(a_x=0.012, b_y=0.027, num_pts=28)

    # IMPORTANT: Do NOT call .resample(160) on the holes so rim edges stay ~5 nm!
    holes = [
        tdgl.Polygon(f"h_{i}", points=ell_pts).translate(dx=x, dy=y)
        for i, (x, y) in enumerate(sym_2x2)
    ]

    source = (
        tdgl.Polygon("source", points=box(1.1 * WIDTH, 0.10))
        .translate(dy=LENGTH / 2)
        .resample(60)
    )
    drain = source.scale(yfact=-1).set_name("drain")

    device = tdgl.Device(
        "3_ellipses_longitudinal",
        layer=LAYER,
        film=film,
        holes=holes,
        terminals=[source, drain],
        probe_points=[(0.0, 0.22), (0.0, -0.22)],
    )
    device.make_mesh(max_edge_length=MAX_EDGE, smooth=10)
    return device


def make_epsilon_func(photon_active: bool, r_photon: tuple[float, float]):
    x0, y0 = r_photon

    def epsilon_func(r, *, t=0.0):
        r = np.asarray(r)
        x, y = r[..., 0], r[..., 1]
        eps = np.ones_like(x, dtype=float)

        if photon_active and (T_PHOTON_ON <= t <= T_PHOTON_OFF):
            dist_sq = (x - x0) ** 2 + (y - y0) ** 2
            eps -= A_HOTSPOT * np.exp(-((dist_sq / (R_HOTSPOT ** 2)) ** 2))

        return eps

    return epsilon_func


def run_trial(
    device: tdgl.Device,
    bias_uA: float,
    photon_active: bool,
    r_photon: tuple[float, float],
    stem: str,
) -> tuple[tdgl.Solution, bool, float]:
    out_h5 = H5_DIR / f"{stem}.h5"
    out_h5.unlink(missing_ok=True)

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
        output_file=str(out_h5),
    )

    solution = tdgl.solve(
        device,
        options,
        applied_vector_potential=0.0,
        terminal_currents=terminal_currents,
        disorder_epsilon=make_epsilon_func(photon_active, r_photon),
    )

    times = solution.dynamics.time
    voltages = np.abs(solution.dynamics.mu[0] - solution.dynamics.mu[1])
    mean_v = float(np.mean(voltages[times > 65.0]))
    return solution, (mean_v > VOLTAGE_THRESHOLD), mean_v


def find_threshold(
    device: tdgl.Device,
    label: str,
    photon_active: bool,
    r_photon: tuple[float, float],
    i_low: float,
    i_high: float,
    tolerance: float = 15.0,
) -> float:
    print(f"\n--- Bisection [{label}] at r_ph={r_photon} (Nodes: {len(device.points)}) ---")
    last_sol = None
    last_i = i_high

    while (i_high - i_low) > tolerance:
        i_mid = 0.5 * (i_low + i_high)
        stem = f"ellipses_{label}_{i_mid:.0f}uA"
        sol, switched, mean_v = run_trial(device, i_mid, photon_active, r_photon, stem)
        status = "SWITCHED" if switched else "RECOVERED"
        print(f"  I = {i_mid:6.1f} uA -> <|V|> = {mean_v:.4f} V0 [{status}]")

        if switched:
            i_high = i_mid
            last_sol = sol
            last_i = i_mid
        else:
            i_low = i_mid

    thresh = 0.5 * (i_low + i_high)
    if last_sol is not None:
        fig, axes = last_sol.plot_order_parameter()
        for ax in axes:
            ax.set_xlim(-WIDTH / 2, WIDTH / 2)
            ax.set_ylim(-0.20, 0.20)
            if photon_active:
                ax.plot(r_photon[0], r_photon[1], "r*", markersize=10)
        fig.suptitle(f"{label}: switched at {last_i:.0f} uA (red star = photon)")
        fig.savefig(FIG_DIR / f"ellipses_{label}_psi.png", bbox_inches="tight")
        plt.close(fig)

    return thresh


def main():
    device = make_elliptical_device()

    # 1. Spontaneous I_sw on this mesh
    i_sw = find_threshold(
        device, "0_spontaneous_Isw", False, (0.0, 0.0), i_low=850.0, i_high=1180.0
    )

    # 2. Test photon impact across different positions (x_ph, y_ph) in the unit cell
    impact_sites = {
        "1_center_bridge_(0,45nm)": (0.000, 0.045),
        "2_near_hole_rim_(20,45nm)": (0.020, 0.045),
        "3_outer_edge_ch_(72,45nm)": (0.072, 0.045),
        "4_inter_row_mid_(0,0nm)": (0.000, 0.000),
    }

    results = {}
    for site_name, r_ph in impact_sites.items():
        i_det = find_threshold(
            device, site_name, True, r_ph, i_low=650.0, i_high=i_sw
        )
        i_norm = i_det / i_sw
        results[site_name] = (r_ph, i_det, i_norm, 1.0 - i_norm)

    print("\n" + "=" * 86)
    print(f"Reference Spontaneous Switching Current I_sw = {i_sw:.1f} uA")
    print("-" * 86)
    print(
        f"{'Photon Impact Site':28s} | {'(x, y) [nm]':>14s} | "
        f"{'I_det [uA]':>10s} | {'i_det = I_det/I_sw':>18s} | {'Margin':>8s}"
    )
    print("-" * 86)
    for site_name, (r_ph, i_det, i_norm, margin) in results.items():
        coord_str = f"({r_ph[0]*1000:.0f}, {r_ph[1]*1000:.0f})"
        print(
            f"{site_name:28s} | {coord_str:>14s} | "
            f"{i_det:10.1f} | {i_norm:18.3f} | {margin*100:7.1f}%"
        )
    print("=" * 86)


if __name__ == "__main__":
    main()