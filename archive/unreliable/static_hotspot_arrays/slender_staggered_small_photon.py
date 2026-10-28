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

# Lower-energy / smaller photon hotspot (r_spot = 11 nm, eps_core = -1.0)
# In this regime, a 200 nm plain wire struggles to detect center impacts!
r_spot = 0.011         # Hotspot radius: 11 nm (22 nm diameter)
amp_spot = -2.0        # Core suppression: eps(0) = 1 - 2 = -1.0
t_relax = 12.0         # Thermal relaxation time [tau_0]
t_max = 60.0           # Dynamic simulation window [tau_0]
BISECT_TOL = 4.0       # Bisection convergence tolerance [uA]

OUT_DIR = "results_slender_small_photon"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Geometries: Plain Wire vs. Slender 2x2 vs. Slender Staggered 2-1-2
# ==============================================================================
def build_geometries():
    geoms = {}

    # 1. Plain Wire (No Holes)
    holes_plain = []
    pts_plain = [
        (0.000, 0.000, 0.30, "1_center_(0,0nm)"),
        (0.030, 0.000, 0.30, "2_mid_ch_(30,0nm)"),
        (0.060, 0.000, 0.25, "3_outer_ch_(60,0nm)"),
        (0.085, 0.000, 0.15, "4_near_edge_(85,0nm)"),
    ]
    geoms["Plain Wire (No Holes)"] = (holes_plain, pts_plain, "k--o", (950.0, 1250.0))

    # 2. Slender 2x2 Ellipses (ax = 8 nm, by = 26 nm, y = +/- 50 nm)
    # Minimal cross-section loss (only 16% blocked) + low crowding (ax/by = 0.31)
    holes_base = [
        tdgl.Polygon(f"h_{i}", points=ellipse(0.008, 0.026, center=(x0, y0), points=28)).resample(28)
        for i, (x0, y0) in enumerate([
            (-0.046,  0.050), (0.046,  0.050),
            (-0.046, -0.050), (0.046, -0.050),
        ])
    ]
    pts_base = [
        (0.000, 0.050, 0.25, "1_row_mid_(0,50nm)"),
        (0.077, 0.050, 0.25, "2_row_edge_(77,50nm)"),
        (0.030, 0.050, 0.20, "3_near_hole_(30,50nm)"),
        (0.000, 0.000, 0.15, "4_center_gap_(0,0nm)"),
        (0.046, 0.012, 0.08, "5_hole_tip_(46,12nm)"),
        (0.085, 0.000, 0.07, "6_far_edge_(85,0nm)"),
    ]
    geoms["Slender 2x2 (ax=8nm)"] = (holes_base, pts_base, "b-^", (880.0, 1120.0))

    # 3. Slender Staggered 2-1-2 (ax_outer = 8 nm, ax_c = 10 nm, dy = 50 nm)
    # Wide 52 nm diagonal throats prevent I_sw choking while eliminating the (0,0) blind spot!
    holes_stag = [
        tdgl.Polygon("h_ul", points=ellipse(0.008, 0.026, center=(-0.048,  0.050), points=28)).resample(28),
        tdgl.Polygon("h_ur", points=ellipse(0.008, 0.026, center=( 0.048,  0.050), points=28)).resample(28),
        tdgl.Polygon("h_c",  points=ellipse(0.010, 0.026, center=( 0.000,  0.000), points=28)).resample(28),
        tdgl.Polygon("h_ll", points=ellipse(0.008, 0.026, center=(-0.048, -0.050), points=28)).resample(28),
        tdgl.Polygon("h_lr", points=ellipse(0.008, 0.026, center=( 0.048, -0.050), points=28)).resample(28),
    ]
    pts_stag = [
        (0.000, 0.050, 0.20, "1_outer_mid_(0,50nm)"),
        (0.078, 0.050, 0.22, "2_outer_edge_(78,50nm)"),
        (0.024, 0.025, 0.22, "3_diag_throat_(24,25nm)"),
        (0.022, 0.000, 0.18, "4_center_hole_edge_(22,0nm)"),
        (0.055, 0.000, 0.10, "5_center_ch_mid_(55,0nm)"),
        (0.085, 0.000, 0.08, "6_far_edge_(85,0nm)"),
    ]
    geoms["Slender Staggered 2-1-2 (ax=8/10nm)"] = (holes_stag, pts_stag, "r-s", (880.0, 1120.0))

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
# 4. Run Sweep & Plot Efficiency Curves
# ==============================================================================
if __name__ == "__main__":
    geoms = build_geometries()
    summary = {}

    print("=" * 92)
    print(f"{'Geometry':<36} | {'Point':<28} | {'I_det [uA]':<10} | {'i_det'}")
    print("-" * 92)

    for g_name, (holes, pts, style, (sw_low, sw_high)) in geoms.items():
        dev = make_device(holes)
        i_sw, _ = find_threshold_bisection(dev, sw_low, sw_high, photon_xy=None, tol=BISECT_TOL)
        print(f"{g_name:<36} | {'Spontaneous I_sw (No Photon)':<28} | {i_sw:<10.1f} | 1.000")

        pt_results = []
        for x0, y0, weight, label in pts:
            i_det, _ = find_threshold_bisection(dev, 0.65 * i_sw, i_sw, photon_xy=(x0, y0), tol=BISECT_TOL)
            ratio = i_det / i_sw
            pt_results.append((i_det, ratio, weight, label))
            print(f"{'':<36} | {label:<28} | {i_det:<10.1f} | {ratio:.3f}")
        print("-" * 92)
        summary[g_name] = {"I_sw": i_sw, "points": pt_results, "style": style}

    # --- Plot Normalized & Absolute Efficiency Curves ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.5))
    i_norm_grid = np.linspace(0.65, 1.01, 350)
    i_abs_grid = np.linspace(600, 1250, 350)

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
    ax1.set_title("Small Hotspot ($r_{\\mathrm{spot}}=11$ nm): Plain vs. Slender Perforated")
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
    out_plot = os.path.join(OUT_DIR, "slender_small_photon_plateau.png")
    fig.savefig(out_plot, dpi=200)
    plt.close(fig)
    print(f"Saved small-photon comparison plot to: {out_plot}")