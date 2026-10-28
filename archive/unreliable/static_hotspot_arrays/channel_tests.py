import os
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
import tdgl
from tdgl.geometry import box, ellipse

# ==============================================================================
# 1. Material & Device Parameters (Matched to your setup)
# ==============================================================================
W = 0.200              # Wire width: 200 nm (in um)
L = 0.800              # Wire length: 800 nm (in um)
d = 0.009              # Film thickness: 9 nm
xi = 0.005             # Coherence length: 5 nm
london_lambda = 0.175  # London penetration depth: 175 nm

# Photon hotspot parameters (fired at t0 = 5.0 tau_0 so we see pre-photon steady state!)
t0_photon = 5.0        # Time when photon hits [tau_0]
r_spot = 0.015         # Hotspot radius: 15 nm (30 nm diameter)
amp_spot = -2.5        # Core suppression: eps_core = 1 - 2.5 = -1.5 < 0
t_relax = 25.0         # Quasiparticle relaxation time [tau_0] (~6-10 ps)
t_max = 65.0           # Total simulation window [tau_0]
BISECT_TOL = 4.0       # Bisection tolerance [uA]

OUT_DIR = "results_channel_cascade"
os.makedirs(OUT_DIR, exist_ok=True)


# ==============================================================================
# 2. Build Simple 2-Channel (1 Line) and 3-Channel (2 Lines) Geometries
# ==============================================================================
def build_channel_geometries():
    geoms = {}

    # Vertical positions for a column of 3 closely-spaced ellipses:
    # b_y = 26 nm (length = 52 nm), centers at y = -64, 0, +64 nm -> vertical gap g_y = 12 nm!
    y_centers = [-0.064, 0.000, 0.064]

    # --- Geometry A: 1 Line of Ellipses at x = 0 -> Forms 2 Sub-Channels ---
    # Ellipse ax = 10 nm, by = 26 nm.
    # Left Channel (Ch 1): x in [-100, -10] nm (90 nm wide)
    # Right Channel (Ch 2): x in [+10, +100] nm (90 nm wide)
    holes_2ch = [
        tdgl.Polygon(f"h_c_{i}", points=ellipse(0.010, 0.026, center=(0.000, y0), points=28)).resample(28)
        for i, y0 in enumerate(y_centers)
    ]
    # Photon hits the center of Left Channel (Ch 1) at (x = -55 nm, y = 0 nm)
    photon_2ch = (-0.055, 0.000)
    # Channel x-integration bounds at y = 0: [(x_min, x_max, label, color)]
    channels_2ch = [
        (-0.098, -0.012, "Ch 1 (Left - HIT BY PHOTON)", "tab:red"),
        ( 0.012,  0.098, "Ch 2 (Right - Neighbor)",     "tab:blue"),
    ]
    geoms["2-Channel (1 Line of Ellipses)"] = (holes_2ch, photon_2ch, channels_2ch, (920.0, 1150.0))

    # --- Geometry B: 2 Lines of Ellipses at x = +/- 36 nm -> Forms 3 Sub-Channels ---
    # Ellipse ax = 8 nm, by = 26 nm at x = -36 nm and +36 nm.
    # Left Channel (Ch 1):   x in [-100, -44] nm (56 nm wide)
    # Center Channel (Ch 2): x in [-28,  +28] nm (56 nm wide)
    # Right Channel (Ch 3):  x in [+44, +100] nm (56 nm wide)
    holes_3ch = []
    for i, y0 in enumerate(y_centers):
        holes_3ch.append(
            tdgl.Polygon(f"h_l_{i}", points=ellipse(0.008, 0.026, center=(-0.036, y0), points=28)).resample(28)
        )
        holes_3ch.append(
            tdgl.Polygon(f"h_r_{i}", points=ellipse(0.008, 0.026, center=( 0.036, y0), points=28)).resample(28)
        )
    # Photon hits the center of Left Channel (Ch 1) at (x = -72 nm, y = 0 nm)
    photon_3ch = (-0.072, 0.000)
    channels_3ch = [
        (-0.098, -0.046, "Ch 1 (Left - HIT BY PHOTON)", "tab:red"),
        (-0.026,  0.026, "Ch 2 (Center - Neighbor 1)",  "tab:orange"),
        ( 0.046,  0.098, "Ch 3 (Right - Neighbor 2)",   "tab:blue"),
    ]
    geoms["3-Channel (2 Lines of Ellipses)"] = (holes_3ch, photon_3ch, channels_3ch, (880.0, 1120.0))

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
        holes=holes,
        terminals=terminals,
        probe_points=[(0, L / 4), (0, -L / 4)],
        length_units="um",
    )
    device.make_mesh(max_edge_length=xi * 1.4, smooth=50)
    return device


# ==============================================================================
# 3. Hotspot & Solver Functions
# ==============================================================================
def make_hotspot_epsilon(x0, y0, t_start=0.0):
    def eps_func(r, *, t):
        r = np.asarray(r)
        if t < t_start:
            return np.ones(r.shape[:-1], dtype=float)
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2
        return 1.0 + amp_spot * np.exp(-dist_sq / (r_spot ** 2)) * np.exp(-(t - t_start) / t_relax)
    return eps_func


def run_simulation(device, current_uA, photon_xy=None, t_start=0.0, save_every=25, h5_name=None):
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
        solve_kwargs["disorder_epsilon"] = make_hotspot_epsilon(*photon_xy, t_start=t_start)

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
    # Check if a sustained voltage / phase-slip cascade occurred during active window t in [15, 50]
    active_mask = (t_arr > 15.0) & (t_arr < 55.0)
    mean_v = np.mean(v_arr[active_mask])
    return mean_v > 0.08, sol


def find_threshold_bisection(device, i_low, i_high, photon_xy=None, tol=BISECT_TOL):
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        switched, _ = run_simulation(
            device, i_mid, photon_xy=photon_xy, t_start=0.0, save_every=200
        )
        if switched:
            i_high = i_mid
        else:
            i_low = i_mid
    return i_high


# ==============================================================================
# 4. Extract Time-Resolved Channel Currents I_ch(t) Across y = 0
# ==============================================================================
def saved_steps_and_times(sol):
    if not sol.path or not os.path.exists(sol.path):
        raise RuntimeError("The solution HDF5 file is unavailable.")
    with h5py.File(sol.path, "r") as h5_file:
        groups = sorted(h5_file["data"].keys(), key=int)
        times = np.array([h5_file["data"][group].attrs["time"] for group in groups])
    return np.array([int(group) for group in groups]), times


def extract_channel_currents_vs_time(device, sol, channels):
    saved_steps, saved_times = saved_steps_and_times(sol)
    n_steps = len(saved_steps)
    # Build integration points for each channel along y = 0
    ch_grids = []
    for x_min, x_max, lbl, col in channels:
        x_pts = np.linspace(x_min, x_max, 40)
        pts = np.column_stack([x_pts, np.zeros_like(x_pts)])
        dx = x_pts[1] - x_pts[0]
        ch_grids.append((pts, dx, lbl, col, np.zeros(n_steps)))

    # Sample every frame
    for frame_idx, solve_step in enumerate(saved_steps):
        sol.solve_step = int(solve_step)
        for ch_idx, (pts, dx, _, _, arr) in enumerate(ch_grids):
            # Extract supercurrent density J_s [uA/um] at channel cross-section
            J_s = sol.interp_current_density(pts, dataset="supercurrent", units="uA/um")
            # Current flows from +y to -y, so magnitude along wire is -J_y
            I_ch = np.trapezoid(np.abs(J_s[:, 1]), dx=dx)
            arr[frame_idx] = I_ch

    return saved_times, ch_grids


# ==============================================================================
# 5. Main Execution: Thresholds + Time-Resolved Channel Cascade Plots
# ==============================================================================
if __name__ == "__main__":
    geoms = build_channel_geometries()

    print("=" * 88)
    print(f"{'Geometry':<34} | {'I_sw [uA]':<12} | {'I_det (Ch1 Hit) [uA]':<22} | {'i_det = I_det/I_sw'}")
    print("-" * 88)

    fig_cascade, axes = plt.subplots(2, 2, figsize=(13.5, 8.5), sharex=True)

    for col_idx, (g_name, (holes, photon_xy, channels, (sw_low, sw_high))) in enumerate(geoms.items()):
        slug = g_name.split("(")[0].strip().replace(" ", "_").replace("-", "_")
        dev = make_device(holes)

        # 1. Find Spontaneous I_sw (No Photon)
        i_sw = find_threshold_bisection(dev, sw_low, sw_high, photon_xy=None)
        # 2. Find Photon Threshold I_det when hitting Channel 1 (Left Channel)
        i_det = find_threshold_bisection(dev, 0.60 * i_sw, i_sw, photon_xy=photon_xy)
        ratio = i_det / i_sw
        print(f"{g_name:<34} | {i_sw:<12.1f} | {i_det:<22.1f} | {ratio:.3f}")

        # 3. Run high-time-resolution cascade trace at I_bias = 0.86 * I_sw (or slightly above i_det)
        i_test = max(0.86 * i_sw, i_det + 8.0)
        print(f" -> Recording time-resolved channel currents at I_bias = {i_test:.1f} uA ({i_test/i_sw:.2f} I_sw)...")
        _, sol = run_simulation(
            dev, i_test, photon_xy=photon_xy, t_start=t0_photon, save_every=15, h5_name=f"{slug}_cascade.h5"
        )

        t_arr, ch_grids = extract_channel_currents_vs_time(dev, sol, channels)
        _, saved_times = saved_steps_and_times(sol)
        saved_indices = np.array(
            [np.argmin(np.abs(sol.dynamics.time - time)) for time in saved_times]
        )
        v_arr = np.abs(sol.dynamics.voltage())[saved_indices]

        # --- Top Panel: Supercurrent in Each Sub-Channel vs. Time ---
        ax_i = axes[0, col_idx]
        for _, _, lbl, color, curr_arr in ch_grids:
            ax_i.plot(t_arr, curr_arr, lw=2.3, color=color, label=lbl)

        ax_i.axvline(t0_photon, color="k", linestyle="--", alpha=0.7, label="Photon Impact ($t_0 = 5\\,\\tau_0$)")
        ax_i.set_title(
            f"{g_name}\n($I_{{sw}}={i_sw:.0f}\\,\\mu$A, Ch1 Hit $i_{{det}}={ratio:.3f}$, Tested at $I_b={i_test:.0f}\\,\\mu$A)",
            fontsize=10.5, fontweight="bold"
        )
        ax_i.set_ylabel("Channel Supercurrent $I_{s,\\mathrm{ch}}(t)$ [$\\mu$A]")
        ax_i.grid(True, alpha=0.4)
        ax_i.legend(loc="upper right", fontsize=8.5)

        # --- Bottom Panel: Terminal Voltage V(t) Across Wire ---
        ax_v = axes[1, col_idx]
        ax_v.plot(t_arr, v_arr, color="purple", lw=2.2, label="Probe Voltage $|\\Delta \\mu(t)|$")
        ax_v.axvline(t0_photon, color="k", linestyle="--", alpha=0.7)
        ax_v.set_xlabel("Time $t$ [$\\tau_0$]")
        ax_v.set_ylabel("Normalized Voltage $|\\Delta \\mu|$ [$V_0$]")
        ax_v.grid(True, alpha=0.4)
        ax_v.legend(loc="upper right", fontsize=8.5)

        # --- Save 3 Snapshots of |psi| During the Cascade (t = 4, 10, 22 tau_0) ---
        fig_snap, snap_axes = plt.subplots(1, 3, figsize=(11, 4.2), sharex=True, sharey=True)
        target_times = [4.0, 9.0, 20.0]
        for s_ax, t_target in zip(snap_axes, target_times):
            step_i = int(np.argmin(np.abs(t_arr - t_target)))
            saved_steps, _ = saved_steps_and_times(sol)
            sol.solve_step = int(saved_steps[step_i])
            psi_mag = np.abs(sol.tdgl_data.psi)
            im = s_ax.tripcolor(
                dev.points[:, 0], dev.points[:, 1], psi_mag,
                triangles=dev.triangles, shading="gouraud", cmap="viridis", vmin=0.0, vmax=1.0
            )
            for hole in dev.holes:
                s_ax.fill(hole.points[:, 0], hole.points[:, 1], color="white", edgecolor="black", lw=0.8)
            s_ax.plot(photon_xy[0], photon_xy[1], "r*", markersize=10, markeredgecolor="w")
            s_ax.set_title(f"$t = {t_arr[step_i]:.1f}\\,\\tau_0$", fontweight="bold")
            s_ax.set_xlim(-W / 2, W / 2)
            s_ax.set_ylim(-0.14, 0.14)
            s_ax.set_aspect("equal")
            s_ax.set_xlabel("$x$ [$\\mu$m]")
        snap_axes[0].set_ylabel("$y$ [$\\mu$m]")
        fig_snap.colorbar(im, ax=snap_axes, shrink=0.85, label="Order Parameter $|\\psi|$")
        fig_snap.suptitle(f"Cascade Evolution in {g_name} at $I_b = {i_test:.0f}\\,\\mu$A", fontsize=12, y=0.98)
        snap_path = os.path.join(OUT_DIR, f"{slug}_cascade_snapshots.png")
        fig_snap.savefig(snap_path, dpi=200, bbox_inches="tight")
        plt.close(fig_snap)
        print(f" -> Saved cascade snapshots to: {snap_path}")

    print("=" * 88)
    plt.tight_layout()
    cascade_plot_path = os.path.join(OUT_DIR, "channel_cascade_2ch_vs_3ch.png")
    fig_cascade.savefig(cascade_plot_path, dpi=200)
    plt.close(fig_cascade)
    print(f"Saved time-resolved channel current cascade plot to: {cascade_plot_path}")