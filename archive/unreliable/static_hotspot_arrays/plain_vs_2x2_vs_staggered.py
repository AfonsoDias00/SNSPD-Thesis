import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tdgl
from tdgl.geometry import box, ellipse

# ==============================================================================
# 1. Material & Device Parameters (Exact match to your setup)
# ==============================================================================
W = 0.200              # Wire width: 200 nm (in um)
L = 0.800              # Wire length: 800 nm (in um)
d = 0.009              # Film thickness: 9 nm
xi = 0.005             # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm

# True normal-core photon hotspot parameters (eps_core = 1 + amp_spot = -2.0 < 0)
r_spot = 0.018         # Hotspot radius: 18 nm (36 nm diameter)
amp_spot = -3.0        # Strong core suppression so |psi| -> 0 inside hotspot
t_relax = 15.0         # Thermal relaxation time [tau_0]
t_max = 60.0           # Dynamic simulation window [tau_0]
BISECT_TOL = 4.0       # Bisection convergence tolerance [uA]

OUT_DIR = "results_plain_vs_holes"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Define Geometries: Plain Wire vs. Baseline 2x2 vs. Staggered 2-1-2
# ==============================================================================
def build_geometries():
    geoms = {}

    # 1. Plain Wire (No Holes) - 4 transverse points from center (x=0) to edge (x=85nm)
    holes_plain = []
    pts_plain = [
        (0.000, 0.000, 0.30, "1_center_(0,0nm)"),
        (0.030, 0.000, 0.30, "2_mid_ch_(30,0nm)"),
        (0.060, 0.000, 0.25, "3_outer_ch_(60,0nm)"),
        (0.085, 0.000, 0.15, "4_near_edge_(85,0nm)"),
    ]
    # Tuple: (holes, pts, plot_style, (i_sw_low, i_sw_high))
    geoms["Plain Wire (No Holes)"] = (holes_plain, pts_plain, "k--o", (950.0, 1350.0))

    # 2. Baseline 2x2 Ellipses (40 nm vertical gap)
    holes_base = [
        tdgl.Polygon(f"h_{i}", points=ellipse(0.012, 0.025, center=(x0, y0), points=28)).resample(28)
        for i, (x0, y0) in enumerate([
            (-0.042,  0.045), (0.042,  0.045),
            (-0.042, -0.045), (0.042, -0.045),
        ])
    ]
    pts_base = [
        (0.000, 0.045, 0.25, "1_row_mid_(0,45nm)"),
        (0.077, 0.045, 0.25, "2_row_edge_(77,45nm)"),
        (0.020, 0.022, 0.20, "3_saddle_(20,22nm)"),
        (0.000, 0.000, 0.15, "4_center_(0,0nm)"),
        (0.042, 0.012, 0.08, "5_hole_tip_(42,12nm)"),
        (0.085, 0.000, 0.07, "6_far_edge_(85,0nm)"),
    ]
    geoms["Baseline 2x2 (40 nm Gap)"] = (holes_base, pts_base, "b-^", (780.0, 1000.0))

    # 3. Staggered 2-1-2 (ax_c = 18 nm, dy = 38 nm)
    holes_stag18 = [
        tdgl.Polygon("h_ul", points=ellipse(0.010, 0.024, center=(-0.046,  0.038), points=28)).resample(28),
        tdgl.Polygon("h_ur", points=ellipse(0.010, 0.024, center=( 0.046,  0.038), points=28)).resample(28),
        tdgl.Polygon("h_c",  points=ellipse(0.018, 0.024, center=( 0.000,  0.000), points=32)).resample(32),
        tdgl.Polygon("h_ll", points=ellipse(0.010, 0.024, center=(-0.046, -0.038), points=28)).resample(28),
        tdgl.Polygon("h_lr", points=ellipse(0.010, 0.024, center=( 0.046, -0.038), points=28)).resample(28),
    ]
    pts_stag18 = [
        (0.000, 0.038, 0.20, "1_outer_mid_(0,38nm)"),
        (0.078, 0.038, 0.22, "2_outer_edge_(78,38nm)"),
        (0.028, 0.019, 0.22, "3_diag_saddle_(28,19nm)"),
        (0.059, 0.000, 0.18, "4_center_ch_(59,0nm)"),
        (0.000, 0.029, 0.08, "5_tip_wake_(0,29nm)"),
        (0.085, 0.000, 0.10, "6_far_edge_(85,0nm)"),
    ]
    geoms["Staggered 2-1-2 (ax_c = 18 nm)"] = (holes_stag18, pts_stag18, "r-s", (780.0, 1000.0))

    return geoms


def make_device(holes):
    layer = tdgl.Layer(coherence_length=xi, london_lambda=london_lambda, thickness=d, gamma=1.0)
    film = tdgl.Polygon("film", points=box(W, L))
    terminals = [
        tdgl.Polygon("source", points=box(W, L / 100, center=(0,  L / 2))),
        tdgl.Polygon("drain",  points=box(W, L / 100, center=(0, -L / 2))),
    ]
    device = tdgl.Device(
        "snspd",
        layer=layer,
        film=film,
        holes=holes if len(holes) > 0 else None,
        terminals=terminals,
        probe_points=[(0, L / 4), (0, -L / 4)],
        length_units="um",
    )
    device.make_mesh(max_edge_length=xi * 1.4, smooth=50)
    return device


# ==============================================================================
# 3. Dynamic Solver & Bisection Functions
# ==============================================================================
def make_hotspot_epsilon(x0, y0):
    def eps_func(r, *, t):
        r = np.asarray(r)
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2
        return 1.0 + amp_spot * np.exp(-dist_sq / (r_spot ** 2)) * np.exp(-t / t_relax)
    return eps_func


def did_switch(device, current_uA, photon_xy=None):
    options = tdgl.SolverOptions(
        solve_time=t_max,
        dt_init=1e-3,
        dt_max=0.25,
        save_every=200,
        field_units="mT",
        current_units="uA",
        output_file=None,
    )
    solve_kwargs = {}
    if photon_xy is not None:
        solve_kwargs["disorder_epsilon"] = make_hotspot_epsilon(*photon_xy)

    sol = tdgl.solve(
        device,
        options,
        terminal_currents=dict(source=current_uA, drain=-current_uA),
        **solve_kwargs,
    )
    if sol is None or sol.dynamics is None:
        raise RuntimeError("TDGL returned no dynamics data.")
    t_arr = sol.dynamics.time
    v_arr = np.abs(sol.dynamics.voltage())
    # Evaluate persistent voltage in the final 35% of the simulation window
    late_v = np.mean(v_arr[t_arr > 0.65 * t_max])
    return late_v > 0.08, sol


def find_threshold_bisection(device, i_low, i_high, photon_xy=None, tol=BISECT_TOL):
    orig_low, orig_high = i_low, i_high
    best_sol = None
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        switched, sol = did_switch(device, i_mid, photon_xy=photon_xy)
        if switched:
            i_high = i_mid
            best_sol = sol
        else:
            i_low = i_mid

    if abs(i_high - orig_low) <= tol:
        print(f"  [WARNING] Bisection hit lower bound ({orig_low:.1f} uA)!")
    elif abs(i_low - orig_high) <= tol:
        print(f"  [WARNING] Bisection hit upper bound ({orig_high:.1f} uA)!")

    return i_high, best_sol


# ==============================================================================
# 4. Run Diagnostic Sweep & Plot Efficiency Curves
# ==============================================================================
if __name__ == "__main__":
    geoms = build_geometries()
    summary = {}

    print("=" * 90)
    print(f"{'Geometry':<32} | {'Point':<26} | {'I_det [uA]':<10} | {'i_det'}")
    print("-" * 90)

    for g_name, (holes, pts, style, (sw_low, sw_high)) in geoms.items():
        dev = make_device(holes)
        i_sw, _ = find_threshold_bisection(dev, sw_low, sw_high, photon_xy=None, tol=BISECT_TOL)
        print(f"{g_name:<32} | {'Spontaneous I_sw (No Photon)':<26} | {i_sw:<10.1f} | 1.000")

        pt_results = []
        for x0, y0, weight, label in pts:
            # Search wide photon window from 0.55 * I_sw up to I_sw
            i_det, _ = find_threshold_bisection(dev, 0.55 * i_sw, i_sw, photon_xy=(x0, y0), tol=BISECT_TOL)
            ratio = i_det / i_sw
            pt_results.append((i_det, ratio, weight, label))
            print(f"{'':<32} | {label:<26} | {i_det:<10.1f} | {ratio:.3f}")
        print("-" * 90)
        summary[g_name] = {"I_sw": i_sw, "points": pt_results, "style": style}

    # --- Plot Normalized & Absolute Efficiency Curves ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.5))
    i_norm_grid = np.linspace(0.55, 1.01, 350)
    i_abs_grid = np.linspace(480, 1300, 350)

    for g_name, data in summary.items():
        i_sw = data["I_sw"]
        pts = data["points"]
        style = data["style"]

        eta_norm = np.zeros_like(i_norm_grid)
        eta_abs = np.zeros_like(i_abs_grid)
        for i_det, ratio, w, _ in pts:
            eta_norm += w * 0.5 * (1.0 + np.tanh((i_norm_grid - ratio) / 0.008))
            eta_abs  += w * 0.5 * (1.0 + np.tanh((i_abs_grid - i_det) / 7.5))

        eta_norm[i_norm_grid > 1.0] = np.nan
        eta_abs[i_abs_grid > i_sw] = np.nan

        ax1.plot(
            i_norm_grid, eta_norm * 100, style, lw=2.4, markevery=35,
            label=f"{g_name} ($I_{{sw}}={i_sw:.0f}\\,\\mu$A)"
        )
        ax2.plot(
            i_abs_grid, eta_abs * 100, style, lw=2.4, markevery=35,
            label=f"{g_name}"
        )
        ax2.axvline(i_sw, color=style[0], linestyle=":", alpha=0.75)

    ax1.axvline(1.0, color="gray", linestyle="--", label="Dark Count Limit ($I_{sw}$)")
    ax1.set_xlabel("Normalized Bias Current $i_b = I_{\\mathrm{bias}} / I_{\\mathrm{sw}}$")
    ax1.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax1.set_title("Normalized Efficiency Plateau: Plain Wire vs. Perforated")
    ax1.set_ylim(-5, 105)
    ax1.grid(True, alpha=0.4)
    ax1.legend(loc="lower right", fontsize=8.5)

    ax2.set_xlabel("Absolute Bias Current $I_{\\mathrm{bias}}$ [$\\mu$A]")
    ax2.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax2.set_title("Absolute Operating Window ($\\eta_{\\mathrm{int}}$ vs. $I_{\\mathrm{bias}}$)")
    ax2.set_ylim(-5, 105)
    ax2.grid(True, alpha=0.4)
    ax2.legend(loc="lower right", fontsize=8.5)

    plt.tight_layout()
    out_plot = os.path.join(OUT_DIR, "plain_vs_2x2_vs_staggered_plateau.png")
    fig.savefig(out_plot, dpi=200)
    plt.close(fig)
    print(f"Saved diagnostic comparison plot to: {out_plot}")