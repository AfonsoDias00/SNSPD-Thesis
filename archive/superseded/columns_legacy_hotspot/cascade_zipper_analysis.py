import os
import tempfile
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tdgl
from tdgl.geometry import box, ellipse

# ==============================================================================
# 1. Material and Device Parameters (NbTiN setup)
# ==============================================================================
W = 0.200              # Wire width: 200 nm (in um)
L = 0.800              # Wire length: 800 nm (in um)
d = 0.009              # Film thickness: 9 nm
xi = 0.005             # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm

# 2D diffusive photon-hotspot parameters (impact at t0 = 20 tau_0)
t0_photon = 20.0       # Photon impact time [tau_0]
r_spot = 0.014         # Initial hotspot radius: 14 nm (28 nm diameter)
amp_spot = -2.5         # Initial core suppression: epsilon_core = -1.5
D_tdgl = 1.0e-5        # 2D electronic diffusion: 10 nm^2 / tau_0
t_relax = 35.0          # Quasiparticle relaxation time [tau_0]
t_max = 85.0            # Total simulation window [tau_0]
BISECT_TOL = 4.0        # Bisection tolerance [uA]

OUT_DIR = "results_cascade_zipper"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Geometries: Plain Wire vs. Aligned Three-Channel vs. Three-Channel Zipper
# ==============================================================================
def build_geometries():
    geoms = {}
    y_left = [-0.064, 0.000, 0.064]           # Column period dy = 64 nm; gaps at y = +/-32 nm

    # --- 1. Plain Wire (No Holes - One 200 nm-Wide Channel) ---
    holes_plain = []
    channels_plain = [
        (-0.098, -0.002, "Left Half (Hit)",      "tab:red"),
        ( 0.002,  0.098, "Right Half (Neighbor)", "tab:blue"),
    ]
    pts_plain = [
        (0.000, 0.000, 0.30, "1_center_(0nm)"),
        (0.030, 0.000, 0.30, "2_intermediate_(30nm)"),
        (0.060, 0.000, 0.25, "3_outer_(60nm)"),
        (0.082, 0.000, 0.15, "4_edge_(82nm)"),
    ]
    geoms["Plain Wire (No Holes)"] = (
        holes_plain, (-0.055, 0.0), channels_plain, pts_plain, "k--o", (1050.0, 1250.0)
    )

    # --- 2. Aligned Three-Channel Geometry (Two Symmetric Rows at x = +/-36 nm) ---
    holes_3ch_aligned = []
    for i, y0 in enumerate(y_left):
        holes_3ch_aligned.append(
            tdgl.Polygon(f"h_l_{i}", points=ellipse(0.008, 0.026, center=(-0.036, y0), points=28)).resample(28)
        )
        holes_3ch_aligned.append(
            tdgl.Polygon(f"h_r_{i}", points=ellipse(0.008, 0.026, center=( 0.036, y0), points=28)).resample(28)
        )
    channels_3ch = [
        (-0.098, -0.046, "Ch 1 (Left - Hit)",       "tab:red"),
        (-0.026,  0.026, "Ch 2 (Center - Neighbor)", "tab:orange"),
        ( 0.046,  0.098, "Ch 3 (Right - Neighbor 2)", "tab:blue"),
    ]
    pts_3ch_aligned = [
        ( 0.000, 0.000, 0.30, "1_Ch2_center_(0,0nm)"),
        (-0.072, 0.000, 0.35, "2_Ch1_middle_(-72,0nm)"),
        (-0.085, 0.000, 0.20, "3_wire_edge_(-85,0nm)"),
        (-0.036, 0.032, 0.15, "4_vertical_gap_(-36,32nm)"),
    ]
    geoms["3-Channel Aligned (2 Lines)"] = (
        holes_3ch_aligned, (-0.072, 0.0), channels_3ch, pts_3ch_aligned, "b-^", (920.0, 1100.0)
    )

    # --- 3. Staggered Three-Channel Zipper (Right Column Offset by dy/2 = +32 nm) ---
    # The right column is shifted so its holes bridge the left-column gaps.
    holes_3ch_zipper = []
    for i, y0 in enumerate([-0.064, 0.000, 0.064]):
        holes_3ch_zipper.append(
            tdgl.Polygon(f"z_l_{i}", points=ellipse(0.008, 0.026, center=(-0.036, y0), points=28)).resample(28)
        )
    for i, y0 in enumerate([-0.096, -0.032, 0.032, 0.096]):
        holes_3ch_zipper.append(
            tdgl.Polygon(f"z_r_{i}", points=ellipse(0.008, 0.026, center=( 0.036, y0), points=28)).resample(28)
        )

    # In the zipper geometry, the two columns alternate their vertical gaps.
    pts_3ch_zipper = [
        ( 0.000, 0.016, 0.30, "1_Ch2_diagonal_middle_(0,16nm)"),
        (-0.072, 0.000, 0.35, "2_Ch1_middle_(-72,0nm)"),
        (-0.085, 0.000, 0.20, "3_wire_edge_(-85,0nm)"),
        (-0.036, 0.032, 0.15, "4_left_gap_(-36,32nm)"),
    ]
    geoms["3-Channel Zipper (Offset dy/2)"] = (
        holes_3ch_zipper, (-0.072, 0.0), channels_3ch, pts_3ch_zipper, "r-s", (920.0, 1120.0)
    )

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
# 3. Energy-Conserving 2D Diffusion Model and Cascade Criterion
# ==============================================================================
def make_hotspot_epsilon(x0, y0, t_start=t0_photon):
    def eps_func(r, *, t):
        r = np.asarray(r)
        if t < t_start:
            return np.ones(r.shape[:-1], dtype=float)
        dt = t - t_start
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2

        # 2D diffusive expansion: r^2(t) = r_0^2 + 4 * D * dt
        r0_sq = r_spot ** 2
        r_t_sq = r0_sq + 4.0 * D_tdgl * dt

        # The 2D Green-function prefactor conserves integrated hotspot energy during diffusion.
        amp_t = amp_spot * (r0_sq / r_t_sq) * np.exp(-dt / t_relax)
        return 1.0 + amp_t * np.exp(-dist_sq / r_t_sq)
    return eps_func


def saved_steps_and_times(sol):
    if not sol.path or not os.path.exists(sol.path):
        raise RuntimeError("The solution HDF5 file is unavailable.")
    with h5py.File(sol.path, "r") as h5_file:
        groups = sorted(h5_file["data"].keys(), key=int)
        times = np.array([h5_file["data"][group].attrs["time"] for group in groups])
    return np.array([int(group) for group in groups]), times


def extract_channel_currents_vs_time(sol, channels, y_slice=0.0):
    saved_steps, saved_times = saved_steps_and_times(sol)
    n_steps = len(saved_steps)
    ch_grids = []
    for x_min, x_max, lbl, col in channels:
        x_pts = np.linspace(x_min, x_max, 35)
        pts = np.column_stack([x_pts, np.full_like(x_pts, y_slice)])
        dx = x_pts[1] - x_pts[0]
        ch_grids.append((pts, dx, lbl, col, np.zeros(n_steps)))

    for frame_idx, solve_step in enumerate(saved_steps):
        sol.solve_step = int(solve_step)
        for _, (pts, dx, _, _, arr) in enumerate(ch_grids):
            J_s = sol.interp_current_density(pts, dataset="supercurrent", units="uA/um")
            arr[frame_idx] = np.trapezoid(np.abs(J_s[:, 1]), dx=dx)

    return saved_times, ch_grids


def run_simulation(device, current_uA, channels, photon_xy=None, save_every=100, h5_name="temp.h5"):
    temporary_output = h5_name == "temp_bisect.h5"
    if temporary_output:
        fd, h5_path = tempfile.mkstemp(
            prefix="temp_bisect-", suffix=".h5", dir=OUT_DIR
        )
        os.close(fd)
        os.unlink(h5_path)
    else:
        h5_path = None if h5_name is None else os.path.join(OUT_DIR, h5_name)

    options = tdgl.SolverOptions(
        solve_time=t_max,
        dt_init=1e-3,
        dt_max=0.20,
        save_every=save_every,
        field_units="mT",
        current_units="uA",
        output_file=h5_path,
    )
    solve_kwargs = {}
    if photon_xy is not None:
        solve_kwargs["disorder_epsilon"] = make_hotspot_epsilon(*photon_xy, t_start=t0_photon)

    sol = None
    try:
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

        # Evaluate only after the photon impact (t > t0 + 2 tau_0).
        post_mask = t_arr > (t0_photon + 2.0)
        max_v_post = np.max(v_arr[post_mask]) if np.any(post_mask) else 0.0
        mean_v_late = np.mean(v_arr[t_arr > 55.0]) if np.any(t_arr > 55.0) else 0.0

        # Integrate channels at the photon y-coordinate, or y = 0 for dark runs.
        y_eval = photon_xy[1] if photon_xy is not None else 0.0
        saved_times, ch_grids = extract_channel_currents_vs_time(sol, channels, y_slice=y_eval)
        post_saved_mask = saved_times > (t0_photon + 2.0)
        last_ch_curr = ch_grids[-1][4][post_saved_mask]
        last_ch_collapsed = (np.min(last_ch_curr) < 50.0) if len(last_ch_curr) > 0 else False

        switched = (max_v_post > 0.45) or (mean_v_late > 0.18) or last_ch_collapsed
        return switched, sol, ch_grids
    finally:
        if temporary_output:
            output_path = getattr(sol, "path", None) or h5_path
            if output_path and os.path.exists(output_path):
                os.remove(output_path)


def find_threshold_bisection(device, channels, i_low, i_high, photon_xy=None, tol=BISECT_TOL):
    orig_low, orig_high = i_low, i_high
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        switched, _, _ = run_simulation(
            device, i_mid, channels, photon_xy=photon_xy, save_every=120, h5_name="temp_bisect.h5"
        )
        if switched:
            i_high = i_mid
        else:
            i_low = i_mid

    if abs(i_high - orig_low) <= tol:
        print(f"  [WARNING] Bisection reached the lower bound ({orig_low:.1f} uA)!")
    elif abs(i_low - orig_high) <= tol:
        print(f"  [WARNING] Bisection reached the upper bound ({orig_high:.1f} uA)!")
    return i_high


# ==============================================================================
# 4. Main Execution: Time Cascade, Synchronized Snapshots, and Curves
# ==============================================================================
if __name__ == "__main__":
    geoms = build_geometries()
    devices = {}
    summary = {}

    print("=" * 94)
    print(f"{'Geometry':<36} | {'Impact Point':<26} | {'I_det [uA]':<10} | {'i_det'}")
    print("-" * 94)

    fig_cascade, axes = plt.subplots(2, 3, figsize=(16, 8.5), sharex=True)

    for col_idx, (g_name, (holes, cascade_xy, channels, pts, style, (sw_low, sw_high))) in enumerate(geoms.items()):
        slug = g_name.split("(")[0].strip().replace(" ", "_").replace("-", "_")
        dev = make_device(holes)
        devices[g_name] = dev

        # 1. Find the spontaneous switching current without a photon.
        i_sw = find_threshold_bisection(dev, channels, sw_low, sw_high, photon_xy=None)
        print(f"{g_name:<36} | {'Spontaneous I_sw (No Photon)':<26} | {i_sw:<10.1f} | 1.000")

        # 2. Run a high-time-resolution cascade simulation at I_bias = 0.86 * I_sw.
        i_test = 0.86 * i_sw
        _, sol_cas, ch_grids = run_simulation(
            dev, i_test, channels, photon_xy=cascade_xy, save_every=15, h5_name=f"{slug}_cascade.h5"
        )
        cascade_steps, cascade_times = saved_steps_and_times(sol_cas)
        dynamics_indices = np.array(
            [np.argmin(np.abs(sol_cas.dynamics.time - time)) for time in cascade_times]
        )
        t_arr = cascade_times
        v_arr = np.abs(sol_cas.dynamics.voltage())[dynamics_indices]

        # --- Top panel: supercurrent in each channel I_ch(t) ---
        ax_i = axes[0, col_idx]
        for _, _, lbl, color, curr_arr in ch_grids:
            ax_i.plot(t_arr, curr_arr, lw=2.3, color=color, label=lbl)
        ax_i.axvline(t0_photon, color="k", linestyle="--", alpha=0.7, label=f"Photon ($t_0={t0_photon:.0f}\\,\\tau_0$)")
        ax_i.set_title(
            f"{g_name}\n($I_{{sw}}={i_sw:.0f}\\,\\mu$A, Testado a $I_b={i_test:.0f}\\,\\mu$A = $0.86\\,I_{{sw}}$)",
            fontsize=10.5, fontweight="bold"
        )
        if col_idx == 0:
            ax_i.set_ylabel("Channel Supercurrent $I_{s,\\mathrm{ch}}(t)$ [$\\mu$A]")
        ax_i.grid(True, alpha=0.4)
        ax_i.legend(loc="upper right", fontsize=8)

        # --- Bottom panel: terminal voltage V(t) ---
        ax_v = axes[1, col_idx]
        ax_v.plot(t_arr, v_arr, color="purple", lw=2.2, label="Voltage $|\\Delta \\mu(t)|$")
        ax_v.axvline(t0_photon, color="k", linestyle="--", alpha=0.7)
        ax_v.set_ylim(-0.1, 3.0)
        ax_v.set_xlabel("Time $t$ [$\\tau_0$]")
        if col_idx == 0:
            ax_v.set_ylabel("Normalized Voltage $|\\Delta \\mu|$ [$V_0$]")
        ax_v.grid(True, alpha=0.4)
        ax_v.legend(loc="upper right", fontsize=8)

        # --- Save four snapshots synchronized with t0 = 20 tau_0 ---
        fig_snap, snap_axes = plt.subplots(1, 4, figsize=(13.5, 4.2), sharex=True, sharey=True)
        target_times = [23.0, 39.0, 56.0, 72.0]
        for s_ax, t_target in zip(snap_axes, target_times):
            frame_i = int(np.argmin(np.abs(cascade_times - t_target)))
            sol_cas.solve_step = int(cascade_steps[frame_i])
            psi_mag = np.abs(sol_cas.tdgl_data.psi)
            im = s_ax.tripcolor(
                dev.points[:, 0], dev.points[:, 1], psi_mag,
                triangles=dev.triangles, shading="gouraud", cmap="viridis", vmin=0.0, vmax=1.0
            )
            for hole in dev.holes:
                s_ax.fill(hole.points[:, 0], hole.points[:, 1], color="white", edgecolor="black", lw=0.8)
            s_ax.plot(cascade_xy[0], cascade_xy[1], "r*", markersize=10, markeredgecolor="w")
            s_ax.set_title(f"$t = {cascade_times[frame_i]:.1f}\\,\\tau_0$", fontweight="bold")
            s_ax.set_xlim(-W / 2, W / 2)
            s_ax.set_ylim(-0.14, 0.14)
            s_ax.set_aspect("equal")
            s_ax.set_xlabel("$x$ [$\\mu$m]")
        snap_axes[0].set_ylabel("$y$ [$\\mu$m]")
        fig_snap.colorbar(im, ax=snap_axes, shrink=0.85, label="Order Parameter $|\\psi|$")
        fig_snap.suptitle(f"Vortex Cascade Evolution: {g_name} ($I_b = {i_test:.0f}\\,\\mu$A)", fontsize=12, y=0.98)
        snap_path = os.path.join(OUT_DIR, f"{slug}_snapshots_t23_39_56_72.png")
        fig_snap.savefig(snap_path, dpi=200, bbox_inches="tight")
        plt.close(fig_snap)

        # 3. Bisection at the four unit-cell impact points.
        pt_results = []
        for x0, y0, weight, label in pts:
            i_det = find_threshold_bisection(dev, channels, 0.65 * i_sw, i_sw, photon_xy=(x0, y0))
            ratio = i_det / i_sw
            pt_results.append((i_det, ratio, weight, label))
            print(f"{'':<36} | {label:<26} | {i_det:<10.1f} | {ratio:.3f}")
        print("-" * 94)
        summary[g_name] = {"I_sw": i_sw, "points": pt_results, "style": style}

    plt.tight_layout()
    cas_path = os.path.join(OUT_DIR, "cascade_comparison_aligned_vs_zipper.png")
    fig_cascade.savefig(cas_path, dpi=200)
    plt.close(fig_cascade)
    print(f"Saved cascade time trace to: {cas_path}")

    # ==========================================================================
    # 5. Final Plot: Internal Detection Efficiency Curves
    # ==========================================================================
    fig_eff, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.5))
    i_norm_grid = np.linspace(0.65, 1.01, 350)
    i_abs_grid = np.linspace(650, 1220, 350)

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

        ax1.plot(i_norm_grid, eta_norm * 100, style, lw=2.4, markevery=35,
                 label=f"{g_name} ($I_{{sw}}={i_sw:.0f}\\,\\mu$A)")
        ax2.plot(i_abs_grid, eta_abs * 100, style, lw=2.4, markevery=35,
                 label=f"{g_name}")
        ax2.axvline(i_sw, color=style[0], linestyle=":", alpha=0.75)

    ax1.axvline(1.0, color="gray", linestyle="--", label="Dark-count limit ($I_{sw}$)")
    ax1.set_xlabel("Normalized Bias Current $i_b = I_{\\mathrm{bias}} / I_{\\mathrm{sw}}$")
    ax1.set_ylabel("Internal Detection Efficiency $\\eta_{\\mathrm{int}}$ [%]")
    ax1.set_title("Normalized Efficiency: Plain Wire vs. Aligned vs. Zipper")
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
    eff_path = os.path.join(OUT_DIR, "efficiency_plateau_aligned_vs_zipper.png")
    fig_eff.savefig(eff_path, dpi=200)
    plt.close(fig_eff)
    print(f"Saved efficiency plateau plot to: {eff_path}")