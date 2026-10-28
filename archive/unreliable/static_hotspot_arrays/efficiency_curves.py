from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tdgl
from tdgl.geometry import box


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results_efficiency_curves"

FIG_DIR = RESULTS_DIR / "figures"
for directory in (RESULTS_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)

# ==============================================================================
# 1. Material & Mesh Parameters (Matching your 200nm wide, xi=5nm, lambda=175nm setup)
# ==============================================================================
W = 0.200          # Wire width: 200 nm (in um)
L = 0.800          # Simulated length: 800 nm (in um)
d = 0.009          # Thickness: 9 nm
xi = 0.005         # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm

# Photon hotspot parameters
r_spot = 0.018     # Effective hotspot radius (18 nm -> 36 nm diameter)
amp_spot = -1.0    # Depress epsilon to -1 (normal core) at impact point
t_relax = 15.0     # Hotspot thermal relaxation time (in units of tau_0)
t_max = 60.0       # Total simulation time per current step
T_RAMP = 20.0
VOLTAGE_THRESHOLD = 0.05

# ==============================================================================
# 2. Define Discrete-Hole Geometries & Unit-Cell Impact Points
# ==============================================================================
def build_geometries():
    geoms = {}

    def ellipse_points(a_x, b_y, num_pts=40):
        angles = np.linspace(0, 2 * np.pi, num_pts, endpoint=False)
        return np.column_stack((a_x * np.cos(angles), b_y * np.sin(angles)))

    def make_holes(points, centers):
        return [
            tdgl.Polygon(f"h_{index}", points=points).translate(dx=x, dy=y)
            for index, (x, y) in enumerate(centers)
        ]

    # 1. Baseline 2x2 Ellipses (40 nm vertical gap at y = 0)
    centers_base = [
        (-0.042, 0.045), (0.042, 0.045),
        (-0.042, -0.045), (0.042, -0.045),
    ]
    holes_base = make_holes(ellipse_points(0.012, 0.025), centers_base)
    # 4 representative impact points in the irreducible quadrant + (x0, y0, area_weight, label)
    pts_base = [
        (0.000, 0.045, 0.30, "1_row_center_(0,45nm)"),
        (0.077, 0.045, 0.35, "2_row_outer_(77,45nm)"),
        (0.020, 0.022, 0.20, "3_saddle_(20,22nm)"),
        (0.000, 0.000, 0.15, "4_inter_row_(0,0nm)"),
    ]
    geoms["Baseline 2x2 (40nm Gap)"] = (holes_base, pts_base, "k--")

    # 2. Close-Packed Discrete 2x2 (10 nm = 2*xi vertical bridge at y = 0)
    # y_center = +/- 31 nm, b_y = 26 nm -> hole tips at +/- 5 nm -> 10 nm bridge!
    centers_cp = [
        (-0.042, 0.031), (0.042, 0.031),
        (-0.042, -0.031), (0.042, -0.031),
    ]
    holes_cp = make_holes(ellipse_points(0.012, 0.026), centers_cp)
    pts_cp = [
        (0.000, 0.031, 0.35, "1_row_center_(0,31nm)"),
        (0.077, 0.031, 0.35, "2_row_outer_(77,31nm)"),
        (0.020, 0.015, 0.15, "3_saddle_(20,15nm)"),
        (0.000, 0.000, 0.15, "4_bridge_mid_(0,0nm)"),
    ]
    geoms["Close-Packed 2x2 (10nm Bridge)"] = (holes_cp, pts_cp, "b-o")

    # 3. Width-Balanced Staggered Lattice (2-1-2 with overlapping y-extents)
    # Outer rows: 2 holes of ax=10nm (40nm total blocked). Center row: 1 hole of ax=18nm (36nm blocked).
    centers_stag = [
        (-0.046, 0.038), (0.046, 0.038), (0.0, 0.0),
        (-0.046, -0.038), (0.046, -0.038),
    ]
    points_stag = ellipse_points(0.010, 0.024)
    holes_stag = make_holes(points_stag, centers_stag)
    holes_stag[2] = tdgl.Polygon("h_c", points=ellipse_points(0.018, 0.024)).translate(
        dx=0.0, dy=0.0
    )
    # Note: at y=0, the center hole occupies |x| <= 18nm, so the active channel at y=0 is at x = 59nm!
    pts_stag = [
        (0.000, 0.038, 0.25, "1_outer_row_mid_(0,38nm)"),
        (0.078, 0.038, 0.30, "2_outer_row_edge_(78,38nm)"),
        (0.028, 0.019, 0.25, "3_stagger_saddle_(28,19nm)"),
        (0.059, 0.000, 0.20, "4_center_row_ch_(59,0nm)"),
    ]
    geoms["Width-Balanced Staggered (2-1-2)"] = (holes_stag, pts_stag, "r-s")

    return geoms


def make_device(holes):
    layer = tdgl.Layer(coherence_length=xi, london_lambda=london_lambda, thickness=d, gamma=1.0)
    film = tdgl.Polygon("film", points=box(W, L)).resample(240)
    terminals = [
        tdgl.Polygon("source", points=box(W, L / 100))
        .translate(dy=L / 2)
        .resample(60),
        tdgl.Polygon("drain", points=box(W, L / 100))
        .translate(dy=-L / 2)
        .resample(60),
    ]
    probe_points = [(0, L / 3), (0, -L / 3)]
    device = tdgl.Device(
        "snspd",
        layer=layer,
        film=film,
        holes=holes,
        terminals=terminals,
        probe_points=probe_points,
        length_units="um",
    )
    device.make_mesh(max_edge_length=xi * 1.4, smooth=10)
    return device


def make_hotspot_epsilon(x0: float, y0: float, photon_active: bool):
    def eps_func(r, *, t=0.0):
        r = np.asarray(r)
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2
        time_decay = np.exp(-t / t_relax)
        if photon_active:
            return 1.0 + amp_spot * np.exp(-dist_sq / r_spot**2) * time_decay
        return np.ones_like(x, dtype=float)
    return eps_func


def did_switch(device, current_uA, photon_xy=None, stem="trial"):
    options = tdgl.SolverOptions(
        solve_time=t_max,
        dt_init=0.01,
        save_every=10,
        current_units="uA",
    )

    def terminal_currents(time):
        current = current_uA * min(time / T_RAMP, 1.0)
        return {"source": current, "drain": -current}

    sol = tdgl.solve(
        device,
        options,
        applied_vector_potential=0.0,
        terminal_currents=terminal_currents,
        disorder_epsilon=make_hotspot_epsilon(
            *(photon_xy or (0.0, 0.0)), photon_xy is not None
        ),
    )
    if sol is None or sol.dynamics is None:
        raise RuntimeError(f"TDGL returned no dynamics for {stem}")
    t_arr = sol.dynamics.time
    v_arr = np.abs(sol.dynamics.mu[0] - sol.dynamics.mu[1])
    late_mask = t_arr > 0.6 * t_max
    if not np.any(late_mask):
        raise RuntimeError(f"No late-time samples available for {stem}")
    late_v = float(np.mean(v_arr[late_mask]))
    return late_v > VOLTAGE_THRESHOLD, sol


def find_threshold_bisection(
    device, i_low, i_high, photon_xy=None, geometry_label="geometry", tol=6.0
):
    best_sol = None
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        photon_label = "photon" if photon_xy is not None else "dark"
        stem = f"{geometry_label}_{photon_label}_{i_mid:.0f}uA".replace(" ", "_")
        switched, sol = did_switch(
            device, i_mid, photon_xy=photon_xy, stem=stem
        )
        if switched:
            i_high = i_mid
            best_sol = sol
        else:
            i_low = i_mid
    return i_high, best_sol


# ==============================================================================
# 3. Run Multi-Point Sweep & Reconstruct the Efficiency Plateau Curve
# ==============================================================================
if __name__ == "__main__":
    geoms = build_geometries()
    summary = {}

    print("=" * 86)
    print(f"{'Geometry':<34} | {'Point':<26} | {'I_det [uA]':<10} | {'i_det'}")
    print("-" * 86)

    for g_name, (holes, pts, style) in geoms.items():
        device = make_device(holes)
        # 1. Find spontaneous switching current I_sw (dark count ceiling)
        i_sw, _ = find_threshold_bisection(
            device,
            900.0,
            1200.0,
            photon_xy=None,
            geometry_label=g_name,
            tol=6.0,
        )
        print(f"{g_name:<34} | {'Spontaneous I_sw (No Photon)':<26} | {i_sw:<10.1f} | 1.000")

        pt_results = []
        for x0, y0, weight, label in pts:
            i_det, _ = find_threshold_bisection(
                device,
                0.72 * i_sw,
                i_sw,
                photon_xy=(x0, y0),
                geometry_label=f"{g_name}_{label}",
                tol=6.0,
            )
            ratio = i_det / i_sw
            pt_results.append((i_det, ratio, weight, label))
            print(f"{'':<34} | {label:<26} | {i_det:<10.1f} | {ratio:.3f}")
        print("-" * 86)
        summary[g_name] = {"I_sw": i_sw, "points": pt_results, "style": style}

    # ==========================================================================
    # 4. Plot Internal Detection Efficiency Plateau Curves (eta_int vs Bias)
    # ==========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    i_norm_grid = np.linspace(0.75, 1.02, 300)
    i_abs_grid = np.linspace(800, 1130, 300)

    for g_name, data in summary.items():
        i_sw = data["I_sw"]
        pts = data["points"]
        style = data["style"]

        # Reconstruct smoothed cumulative efficiency curve (broadened by ~1% thermal/spatial过渡)
        eta_norm = np.zeros_like(i_norm_grid)
        eta_abs = np.zeros_like(i_abs_grid)
        for i_det, ratio, w, _ in pts:
            eta_norm += w * 0.5 * (1.0 + np.tanh((i_norm_grid - ratio) / 0.008))
            eta_abs  += w * 0.5 * (1.0 + np.tanh((i_abs_grid - i_det) / 8.0))

        # Mask out currents above I_sw (where dark counts saturate)
        eta_norm[i_norm_grid > 1.0] = np.nan
        eta_abs[i_abs_grid > i_sw] = np.nan

        ax1.plot(i_norm_grid, eta_norm * 100, style, lw=2.5, markevery=35,
                 label=f"{g_name} ($I_{{sw}}={i_sw:.0f}\\,\\mu$A)")
        ax2.plot(i_abs_grid, eta_abs * 100, style, lw=2.5, markevery=35,
                 label=f"{g_name}")
        ax2.axvline(i_sw, color=style[0], linestyle=":", alpha=0.7)

    ax1.axvline(1.0, color="gray", linestyle="--", label="Dark Count Limit ($I_{sw}$)")
    ax1.set_xlabel("Normalized Bias Current $i_b = I_{\\mathrm{bias}} / I_{\\mathrm{sw}}$")
    ax1.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax1.set_title("Normalized Efficiency Plateau ($\\eta_{\\mathrm{int}}$ vs. $I_b/I_{sw}$)")
    ax1.set_ylim(-5, 105)
    ax1.grid(True, alpha=0.4)
    ax1.legend(loc="lower right", fontsize=9)

    ax2.set_xlabel("Absolute Bias Current $I_{\\mathrm{bias}}$ [$\\mu$A]")
    ax2.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax2.set_title("Absolute Operating Window (Plateau vs. $I_{sw}$ Wall)")
    ax2.set_ylim(-5, 105)
    ax2.grid(True, alpha=0.4)
    ax2.legend(loc="lower right", fontsize=9)

    plt.tight_layout()
    plt.savefig(FIG_DIR / "efficiency_plateau_comparison.png", dpi=200)
    