"""
tmax_sensitivity.py

Does the photon detection threshold i_det = I_det / I_sw depend on the length of
the simulation window t_max?

Method
------
Every solver run goes out to the LONGEST window (default 150 tau_0). The detection
criterion from photon_hitmap.py is then re-evaluated on the truncated data for each
shorter window T (the TDGL trajectory up to T is identical whatever t_max is, so
this is exact, not an approximation). Bisections for different windows share a
common bracket, so most solver runs are reused between windows via a JSON cache.
The spontaneous I_sw is also re-determined per window, because a longer window can
lower it, and i_det is normalised by I_sw(T).

Criterion per window T (same constants as photon_hitmap.py):
  switched = (max |V| for t in (t0+4, T] > 0.45)
          or (mean |V| for t in (55, T] > 0.18)
          or (monitored-channel supercurrent < 50 uA at any saved step in (t0+4, T])

Usage
-----
  python tmax_sensitivity.py --quick          # fewer points, windows 80 and 150
  python tmax_sensitivity.py                  # full run
  python tmax_sensitivity.py --mesh-factor 0.7 --out results_tmax_fine
      # second check: finer mesh (max_edge_length = 0.7 xi instead of 1.4 xi)

Results are cached in <out>/cache.json, so an interrupted run can be restarted.
Outputs: <out>/tmax_results.csv and <out>/tmax_sensitivity.png
"""
import argparse
import csv
import json
import os
import tempfile
import time

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tdgl
from tdgl.geometry import box, ellipse

trapz = getattr(np, "trapezoid", None) or np.trapz

# ==============================================================================
# 1. Parameters (identical to photon_hitmap.py)
# ==============================================================================
W = 0.200
L = 0.800
d = 0.009
xi = 0.005
london_lambda = 0.175

t0_photon = 20.0
r_spot = 0.014
amp_spot = -2.5
D_tdgl = 1.0e-5
t_relax = 35.0

# Detection criterion constants
T_SETTLE = t0_photon + 4.0   # ignore the initial transient
V_POST_MAX = 0.45
T_LATE_START = 55.0
V_LATE_MAX = 0.18
I_COLLAPSE = 50.0            # uA, monitored-channel supercurrent collapse


# ==============================================================================
# 2. Geometry
# ==============================================================================
def make_device(holes, mesh_factor):
    layer = tdgl.Layer(coherence_length=xi, london_lambda=london_lambda, thickness=d, gamma=1.0)
    film = tdgl.Polygon("film", points=box(W, L))
    terminals = [
        tdgl.Polygon("source", points=box(W, L / 100, center=(0, L / 2))),
        tdgl.Polygon("drain", points=box(W, L / 100, center=(0, -L / 2))),
    ]
    device = tdgl.Device(
        "snspd",
        layer=layer,
        film=film,
        holes=holes if len(holes) > 0 else None,
        terminals=terminals,
        probe_points=[(0, L / 3), (0, -L / 3)],
        length_units="um",
    )
    device.make_mesh(max_edge_length=xi * mesh_factor, smooth=50)
    return device


def build_long_zipper_holes():
    holes = []
    for i, y0 in enumerate([-0.122, 0.000, 0.122]):
        holes.append(
            tdgl.Polygon(f"l_l_{i}", points=ellipse(0.008, 0.055, center=(-0.036, y0), points=38)).resample(38)
        )
    for i, y0 in enumerate([-0.183, -0.061, 0.061, 0.183]):
        holes.append(
            tdgl.Polygon(f"l_r_{i}", points=ellipse(0.008, 0.055, center=(0.036, y0), points=38)).resample(38)
        )
    return holes


def inside_any_hole(x, y, holes, margin=0.003):
    pt = np.array([[x, y]])
    for h in holes:
        if h.contains_points(pt)[0]:
            return True
        if np.min(np.linalg.norm(h.points - pt, axis=1)) < margin:
            return True
    return False


# ==============================================================================
# 3. Hotspot and single solver run
# ==============================================================================
def make_hotspot_epsilon(x0, y0, t_start=t0_photon):
    def eps_func(r, *, t):
        r = np.asarray(r)
        if t < t_start:
            return np.ones(r.shape[:-1], dtype=float)
        dt = t - t_start
        x, y = r[..., 0], r[..., 1]
        dist_sq = (x - x0) ** 2 + (y - y0) ** 2
        r0_sq = r_spot ** 2
        r_t_sq = r0_sq + 4.0 * D_tdgl * dt
        amp_t = amp_spot * (r0_sq / r_t_sq) * np.exp(-dt / t_relax)
        return 1.0 + amp_t * np.exp(-dist_sq / r_t_sq)

    return eps_func


def saved_steps_and_times(sol):
    with h5py.File(sol.path, "r") as h5_file:
        groups = sorted(h5_file["data"].keys(), key=int)
        times = np.array([h5_file["data"][g].attrs["time"] for g in groups])
    return np.array([int(g) for g in groups]), times


def simulate(device, current_uA, photon_xy, slc, windows, t_long, save_every, out_dir):
    """One solver run to t_long. Returns the detection outcome for every window."""
    fd, h5_path = tempfile.mkstemp(prefix="tmax_", suffix=".h5", dir=out_dir)
    os.close(fd)
    os.unlink(h5_path)

    options = tdgl.SolverOptions(
        solve_time=t_long,
        dt_init=1e-3,
        dt_max=0.20,
        save_every=save_every,
        field_units="mT",
        current_units="uA",
        output_file=h5_path,
    )
    kwargs = {}
    if photon_xy is not None:
        kwargs["disorder_epsilon"] = make_hotspot_epsilon(*photon_xy, t_start=t0_photon)

    sol = None
    try:
        sol = tdgl.solve(
            device,
            options,
            terminal_currents=dict(source=current_uA, drain=-current_uA),
            **kwargs,
        )
        t = np.asarray(sol.dynamics.time)
        v = np.abs(sol.dynamics.voltage())

        # first time the voltage criterion trips
        over = np.where((t > T_SETTLE) & (v > V_POST_MAX))[0]
        t_v = float(t[over[0]]) if len(over) else None

        # first saved time the monitored channel has collapsed
        x_min, x_max, y_sl = slc
        x_pts = np.linspace(x_min, x_max, 30)
        pts = np.column_stack([x_pts, np.full_like(x_pts, y_sl)])
        dx = x_pts[1] - x_pts[0]
        t_c = None
        steps, times = saved_steps_and_times(sol)
        for s_step, s_time in zip(steps, times):
            if s_time <= T_SETTLE:
                continue
            sol.solve_step = int(s_step)
            J_s = sol.interp_current_density(pts, dataset="supercurrent", units="uA/um")
            if trapz(np.abs(J_s[:, 1]), dx=dx) < I_COLLAPSE:
                t_c = float(s_time)
                break

        switched = {}
        for T in windows:
            late = (t > T_LATE_START) & (t <= T)
            mean_late = float(np.mean(v[late])) if np.any(late) else 0.0
            switched[str(T)] = bool(
                (t_v is not None and t_v <= T)
                or (t_c is not None and t_c <= T)
                or (mean_late > V_LATE_MAX)
            )
        triggers = [x for x in (t_v, t_c) if x is not None]
        return {"switched": switched, "t_trig": (min(triggers) if triggers else None)}
    finally:
        path = getattr(sol, "path", None) or h5_path
        if path and os.path.exists(path):
            os.remove(path)


# ==============================================================================
# 4. Cache + bisection
# ==============================================================================
class Cache:
    def __init__(self, path):
        self.path = path
        self.data = {}
        if os.path.exists(path):
            with open(path) as f:
                self.data = json.load(f)

    def save(self):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.data, f)
        os.replace(tmp, self.path)


class Runner:
    def __init__(self, tag, device, slc, windows, t_long, save_every, mesh_factor, out_dir, cache):
        self.tag, self.device, self.slc = tag, device, slc
        self.windows, self.t_long, self.save_every = windows, t_long, save_every
        self.mesh_factor, self.out_dir, self.cache = mesh_factor, out_dir, cache
        self.n_new = 0

    def __call__(self, current, photon_xy):
        who = "dark" if photon_xy is None else f"{photon_xy[0] * 1e3:.1f},{photon_xy[1] * 1e3:.1f}"
        key = f"{self.tag}|m{self.mesh_factor}|se{self.save_every}|tl{self.t_long}|{who}|{current:.2f}"
        hit = self.cache.data.get(key)
        if hit is not None and all(str(T) in hit["switched"] for T in self.windows):
            return hit
        t_start = time.time()
        res = simulate(self.device, current, photon_xy, self.slc, self.windows,
                       self.t_long, self.save_every, self.out_dir)
        self.cache.data[key] = res
        self.cache.save()
        self.n_new += 1
        print(f"      run #{self.n_new}: {who:>12s} nm @ {current:7.1f} uA "
              f"({time.time() - t_start:5.0f} s, first trigger t = {res['t_trig']})", flush=True)
        return res


def bisect(runner, i_low, i_high, tol, T, photon_xy):
    orig_low, orig_high = i_low, i_high
    while (i_high - i_low) > tol:
        i_mid = 0.5 * (i_low + i_high)
        if runner(i_mid, photon_xy)["switched"][str(T)]:
            i_high = i_mid
        else:
            i_low = i_mid
    flag = ""
    if i_high - orig_low <= tol:
        flag = "LOWER_BRACKET"   # switched even at the lowest current tried
    elif orig_high - i_low <= tol:
        flag = "UPPER_BRACKET"   # never switched inside the bracket
    return i_high, flag


# ==============================================================================
# 5. Main
# ==============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", type=float, nargs="+", default=[60.0, 80.0, 100.0, 125.0, 150.0])
    ap.add_argument("--tol", type=float, default=4.0, help="bisection tolerance [uA]")
    ap.add_argument("--mesh-factor", type=float, default=1.4, help="max_edge_length in units of xi")
    ap.add_argument("--save-every", type=int, default=40,
                    help="snapshot interval in solver steps (photon_hitmap.py used 120, which is coarse)")
    ap.add_argument("--out", default="results_tmax_sensitivity")
    ap.add_argument("--quick", action="store_true", help="fewer points, windows 80 and 150 only")
    args = ap.parse_args()

    windows = sorted(args.windows)
    if args.quick:
        windows = [80.0, 150.0]
    t_long = windows[-1]
    ref_T = 80.0 if 80.0 in windows else windows[0]
    os.makedirs(args.out, exist_ok=True)
    cache = Cache(os.path.join(args.out, "cache.json"))

    plain_pts = [(-0.088, 0.0), (-0.058, 0.0), (-0.028, 0.0), (0.0, 0.0)]
    zip_pts = [(-0.088, 0.0), (-0.068, 0.0), (-0.050, 0.030), (-0.018, 0.030), (0.0, 0.0), (-0.036, 0.061)]
    if args.quick:
        plain_pts = [(-0.088, 0.0), (0.0, 0.0)]
        zip_pts = [(-0.088, 0.0), (0.0, 0.0), (-0.036, 0.061)]

    holes_zip = build_long_zipper_holes()
    geoms = [
        # tag, holes, monitored slice, dark-bracket, points
        ("Plain_Wire", [], (0.002, 0.098, 0.000), (1000.0, 1250.0), plain_pts),
        ("Long_Zipper", holes_zip, (0.046, 0.098, 0.061), (880.0, 1120.0), zip_pts),
    ]

    rows = []
    csv_path = os.path.join(args.out, "tmax_results.csv")

    def write_csv():
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Geometry", "x_nm", "y_nm", "t_max", "I_sw_uA", "I_det_uA", "i_det", "flag"])
            w.writerows(rows)

    for tag, holes, slc, (sw_lo, sw_hi), pts in geoms:
        print("=" * 90)
        print(f"{tag}: building mesh (max_edge_length = {args.mesh_factor} xi)", flush=True)
        dev = make_device(holes, args.mesh_factor)
        runner = Runner(tag, dev, slc, windows, t_long, args.save_every, args.mesh_factor, args.out, cache)

        # --- spontaneous I_sw for every window ---------------------------------
        print("  Spontaneous I_sw per window:", flush=True)
        i_sw = {}
        for T in windows:
            i_sw[T], flag = bisect(runner, sw_lo, sw_hi, args.tol, T, None)
            print(f"    T = {T:5.0f} tau_0 : I_sw = {i_sw[T]:7.1f} uA {flag}", flush=True)
            rows.append([tag, "", "", T, i_sw[T], i_sw[T], 1.0, flag or "dark"])
        write_csv()

        # common bracket for all windows so that solver runs are shared via the cache
        lo = 0.72 * min(i_sw.values())
        hi = max(i_sw.values())

        # --- photon thresholds -------------------------------------------------
        for x0, y0 in pts:
            if holes and inside_any_hole(x0, y0, holes):
                print(f"  ({x0 * 1e3:6.1f}, {y0 * 1e3:5.1f}) nm  [inside hole, skipped]")
                continue
            print(f"  Hit at ({x0 * 1e3:6.1f}, {y0 * 1e3:5.1f}) nm", flush=True)
            line = []
            for T in windows:
                i_det, flag = bisect(runner, lo, hi, args.tol, T, (x0, y0))
                ratio = i_det / i_sw[T]
                if not flag and ratio >= 0.99:
                    flag = "DARK_LIMITED"
                rows.append([tag, x0 * 1e3, y0 * 1e3, T, i_sw[T], i_det, ratio, flag])
                line.append(f"T={T:.0f}: {ratio:.3f}{'*' if flag else ''}")
            print("    " + "   ".join(line), flush=True)
            write_csv()

    write_csv()

    # ==========================================================================
    # 6. Summary table + plot
    # ==========================================================================
    print("\n" + "=" * 90)
    print(f"i_det vs window (reference window T = {ref_T:.0f}, longest T = {t_long:.0f}); '*' = flagged")
    print("-" * 90)
    header = f"{'Geometry':<12}{'x,y [nm]':<14}" + "".join(f"{'T=' + str(int(T)):>9}" for T in windows) + f"{'shift':>9}"
    print(header)
    series = {}
    for r in rows:
        tag, x, y, T, isw, idet, ratio, flag = r
        if x == "":
            continue
        series.setdefault((tag, round(x, 1), round(y, 1)), {})[T] = (ratio, flag)
    for (tag, x, y), d_ in series.items():
        vals = "".join(f"{d_[T][0]:>8.3f}{'*' if d_[T][1] else ' '}" for T in windows)
        shift = d_[t_long][0] - d_[ref_T][0]
        print(f"{tag:<12}{f'{x:.0f},{y:.0f}':<14}{vals}{shift:>+9.3f}")
    print("-" * 90)
    print(f"Bisection resolution is about {args.tol / 1000:.3f} in i_det (tol / I_sw). "
          "Shifts well above ~0.01 mean the window matters.")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    for tag, color in (("Plain_Wire", "k"), ("Long_Zipper", "r")):
        sw = [(r[3], r[4]) for r in rows if r[0] == tag and r[1] == ""]
        ax1.plot([a for a, _ in sw], [b for _, b in sw], "o-", color=color, label=tag.replace("_", " "))
    ax1.set_xlabel("Simulation window $t_{max}$ [$\\tau_0$]")
    ax1.set_ylabel("Spontaneous $I_{sw}$ [$\\mu$A]")
    ax1.set_title("Dark switching current vs window")
    ax1.grid(True, alpha=0.4)
    ax1.legend()

    for (tag, x, y), d_ in series.items():
        Ts = sorted(d_)
        ls = "--" if tag == "Plain_Wire" else "-"
        ax2.plot(Ts, [d_[T][0] for T in Ts], ls, marker="o", label=f"{tag.replace('_', ' ')} ({x:.0f},{y:.0f})")
    ax2.axvline(80.0, color="gray", ls=":", label="window used in photon_hitmap.py")
    ax2.set_xlabel("Simulation window $t_{max}$ [$\\tau_0$]")
    ax2.set_ylabel("$i_{det} = I_{det} / I_{sw}(t_{max})$")
    ax2.set_title("Photon threshold vs window")
    ax2.grid(True, alpha=0.4)
    ax2.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    png_path = os.path.join(args.out, "tmax_sensitivity.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved {csv_path}\nSaved {png_path}")


if __name__ == "__main__":
    main()
