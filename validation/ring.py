"""
2b  A closed ring: an annulus R1 = 13.5 xi, R2 = 16.5 xi (3 xi wide, about as narrow as the zipper's 12 nm bridges),
no terminals, no field.  psi starts as exp(i n phi) (times 1 % random noise to seed instabilities) and relaxes for
600 tau0, on three meshes.

  winding   a stable state keeps its winding: the meter reports n at every frame
  current   psi = f(r) exp(i n phi) with f'' + f'/r + (1 - n^2/r^2 - f^2) f = 0, f' = 0 at both edges (solved here
            as a BVP).  The circulating current is I = (K0/4) xi int f^2 n/r dr: pyTDGL's sheet current is
            (K0/4) Im(psi* grad psi) with grad in 1/xi (K0 = device.K0; the plan's "(1 - q^2) q in pyTDGL units"
            needs this factor 1/4, as j_dep = K0 / (6 sqrt 3) in snspd/units.py).  Compared with the current
            through radial cuts of the TDGL state.
  depairing a thin ring of radius R holds winding n only while q = n/R satisfies q^2 < 1/3 + 1/(6 R^2) (Eckhaus /
            depairing limit; the second term is the finite-ring correction): n <= 8 at R = 15 xi, n <= 7 if the
            inner radius sets it.  Above it, phase slips take n down by one at a time.  The plan expected n = 12 to
            end at 8 or slightly below; in this model (no heat) a phase-slip centre, once formed, keeps slipping
            until the current is far below depairing, so the end state is measured against that instead.
  slips     the same n = 12 run with a frame every 2 solver steps over the slip cascade: each slip a step of -1?
  mesh      the stable set and the end states from n = 9..12 on 1.4, 0.8 and 0.5 xi meshes (2, 4, 6 edges across)
"""
from __future__ import annotations

import time

import numpy as np

from snspd import sim
from snspd.fluxoid import PSI_LOOP_MIN, FluxoidMeter
from snspd.units import FILM

from . import common as C

ID, TITLE = "2b", "closed ring: winding, persistent current, depairing limit"
R1, R2 = 13.5, 16.5                    # [xi]
MESHES = (1.4, 0.8, 0.5)
N_INIT = tuple(range(13))
SEEDS_12 = (0, 1, 2)
T_RUN = 600.0
_DEV = {}


def device(mesh_xi):
    if mesh_xi not in _DEV:
        dev = C.annulus(R1 * C.XI, R2 * C.XI, mesh_xi)
        _DEV[mesh_xi] = (dev, FluxoidMeter(dev))
    return _DEV[mesh_xi]


def ring_current(dev, sol, n_lines=16, n_pts=61):
    """Circulating current [uA], mean over radial cuts: from psi (common.supercurrent_through) and from pyTDGL's
    interpolated current density."""
    r = np.linspace(R1 * C.XI + 0.01, R2 * C.XI - 0.01, n_pts) * 1e-3
    psi = np.asarray(sol.tdgl_data.psi)
    from_psi, from_tdgl = [], []
    for ph in np.linspace(0, 2 * np.pi, n_lines, endpoint=False):
        xy = np.column_stack([r * np.cos(ph), r * np.sin(ph)])
        from_psi.append(C.supercurrent_through(dev, psi, xy))
        from_tdgl.append(C.current_through(sol, xy))
    return float(np.mean(from_psi)), float(np.mean(from_tdgl))


def task_fn(task):
    t_start = time.time()
    try:
        dev, meter = device(task["mesh"])
        xy = np.asarray(dev.points)
        rng = np.random.default_rng(task["seed"] * 1000 + task["n"])
        phase = np.exp(1j * task["n"] * np.arctan2(xy[:, 1], xy[:, 0]))
        phase *= 1 + 0.01 * (rng.standard_normal(len(xy)) + 1j * rng.standard_normal(len(xy)))
        with C.solved_free(dev, task.get("T", T_RUN), phase=phase, save_every=task.get("every", 10)) as sol:
            frames = [(t, *meter(p)) for t, p in C.psi_frames(sol)]
            frames = [(t, n[0], fl) for t, n, fl in frames]
            res = dict(frames=frames, n_end=frames[-1][1])
            ok = [(t, n) for t, n, fl in frames if fl >= PSI_LOOP_MIN]
            res["steps"] = [(t2, n2 - n1) for (t1, n1), (t2, n2) in zip(ok, ok[1:]) if n2 != n1]
            res["rejected"] = sum(fl < PSI_LOOP_MIN for _, _, fl in frames)
            if not res["steps"] and res["n_end"] == task["n"]:
                res["I"], res["I_tdgl"] = ring_current(dev, sol)
                mid = R1 + (R2 - R1) / 2
                ph = np.linspace(0, 2 * np.pi, 200, endpoint=False)
                op = C.interp_op(dev, np.column_stack([np.cos(ph), np.sin(ph)]) * mid * C.XI * 1e-3)
                res["f2_mid"] = float(np.mean(np.abs(op @ frames_psi_last(sol)) ** 2))
        return dict(ok=True, task=task, res=res, wall_s=time.time() - t_start, nodes=len(xy))
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def frames_psi_last(sol):
    return np.asarray(sol.tdgl_data.psi)


def gl_profile(n):
    """f(r) of the annulus with winding n (r in xi), by solve_bvp; None if it does not converge."""
    from scipy.integrate import solve_bvp
    r = np.linspace(R1, R2, 101)
    f0 = np.sqrt(np.clip(1 - n ** 2 / r ** 2, 0.05, 1))
    sol = solve_bvp(lambda r, y: np.vstack([y[1], -y[1] / r - (1 - n ** 2 / r ** 2 - y[0] ** 2) * y[0]]),
                    lambda ya, yb: np.array([ya[1], yb[1]]), r, np.vstack([f0, np.zeros_like(r)]), tol=1e-8)
    return sol.sol if sol.success else None


def run(a):
    tasks = [dict(mesh=m, n=n, seed=0) for m in MESHES for n in N_INIT]
    tasks += [dict(mesh=m, n=12, seed=s) for m in MESHES for s in SEEDS_12 if s]
    tasks += [dict(mesh=MESHES[-1], n=12, seed=0, T=100.0, every=2, dense=True)]     # every=1 breaks pyTDGL 0.9 loading
    out = {}

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED mesh {t['mesh']} n {t['n']} seed {t['seed']}: {r['error']}")
            return
        out[(t["mesh"], t["n"], t["seed"]) + (("dense",) if t.get("dense") else ())] = r
        q = r["res"]
        sim.log(f"   mesh {t['mesh']:g} ({r['nodes']} sites) n = {t['n']:2d} seed {t['seed']}: ends {q['n_end']}, "
                f"steps {[d for _, d in q['steps']]}, {q['rejected']} rejected frames  ({r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result, fn=task_fn)
    checks, extra = [], {}
    K0, mid = FILM.K0_uA_per_nm, R1 + (R2 - R1) / 2
    n_ecks = lambda R: R * np.sqrt(1 / 3 + 1 / (6 * R ** 2))

    stable = {m: [n for n in N_INIT if not out[(m, n, 0)]["res"]["steps"] and out[(m, n, 0)]["res"]["n_end"] == n]
              for m in MESHES}
    n_max = {m: max(s) for m, s in stable.items()}
    glitch = [(m, n) for m in MESHES for n in stable[m] if out[(m, n, 0)]["res"]["rejected"]
              or any(k != n for _, k, _ in out[(m, n, 0)]["res"]["frames"])]
    checks.append(C.check("2b.1", "stable rings: winding read at every frame (all meshes)",
                          f"{sum(map(len, stable.values())) - len(glitch)} of {sum(map(len, stable.values()))} runs "
                          "clean", "every frame accepted and = n", not glitch, f"glitches: {glitch}" if glitch else ""))

    # current and amplitude against the radial GL solution, finest mesh first
    rows, worst = [], {}
    for m in MESHES:
        errs = []
        for n in range(1, n_max[m] + 1):
            prof = gl_profile(n)
            if prof is None:
                continue
            rr = np.linspace(R1, R2, 401)
            f2 = prof(rr)[0] ** 2
            I_gl = K0 * C.XI / 4 * np.trapezoid(f2 * n / rr, rr)
            I_loc = K0 * C.XI / 4 * np.trapezoid((1 - n ** 2 / rr ** 2) * n / rr, rr)
            q = out[(m, n, 0)]["res"]
            errs.append((n, q["I"] / I_gl - 1, q["f2_mid"] / prof(mid)[0] ** 2 - 1, I_loc / I_gl - 1,
                         q["I_tdgl"] / I_gl - 1))
            rows.append(dict(mesh=m, n=n, I_psi=q["I"], I_tdgl=q["I_tdgl"], I_gl=I_gl, I_local=I_loc,
                             f2_mid=q["f2_mid"], f2_gl=float(prof(mid)[0] ** 2)))
        worst[m] = errs
    fine = MESHES[-1]
    eI = max(abs(e[1]) for e in worst[fine])
    ef = max(abs(e[2]) for e in worst[fine])
    checks.append(C.check("2b.2", f"circulating current from psi vs (K0/4) xi int f^2 n/r dr (radial GL), n = 1..8, "
                          f"{fine:g} xi", f"max |I / I_GL - 1| = {eI:.4f}", "< 0.01 (coarser meshes in the note)",
                          eI < 0.01,
                          "; ".join(f"{m:g} xi: " + ", ".join(f"{n}:{e:+.3f}" for n, e, _, _, _ in worst[m]) for m in MESHES)
                          + ";  uniform-state formula (1 - q^2) q at the local q = n/r vs GL: "
                          + ", ".join(f"{n}:{e:+.3f}" for n, _, _, e, _ in worst[fine])))
    checks.append(C.check("2b.2b", "the same current from pyTDGL's interpolated current density (what hitmap.py's "
                          "crowding analysis reads)", " / ".join(f"{m:g} xi: {max(e[4] for e in worst[m]):+.3f}"
                                                                for m in MESHES),
                          "I / I_GL - 1 (largest over n)", note="a ring 3 xi wide is mostly boundary sites at 1.4 xi; "
                          "the excess falls with the mesh"))
    checks.append(C.check("2b.3", f"|psi|^2 at mid-radius vs the radial GL solution, {fine:g} xi",
                          f"max |rel. diff| = {ef:.4f}", "< 0.01", ef < 0.01,
                          "; ".join(f"{m:g} xi: " + ", ".join(f"{n}:{e:+.4f}" for n, _, e, _, _ in worst[m]) for m in MESHES)))

    # depairing limit
    contiguous = all(stable[m] == list(range(n_max[m] + 1)) for m in MESHES)
    checks.append(C.check("2b.4", "largest winding that survives 600 tau0 (1.4 / 0.8 / 0.5 xi)",
                          " / ".join(str(n_max[m]) for m in MESHES),
                          f"8 (q^2 < 1/3 + 1/6R^2 at R = 15 xi: n < {n_ecks(15.0):.2f}); 7 if set at R1 (n < "
                          f"{n_ecks(R1):.2f}); all lower n stable", contiguous and all(n_max[m] in (7, 8) for m in MESHES),
                          f"q = n_max / R_mid = " + ", ".join(f"{n_max[m] / mid:.3f}" for m in MESHES)
                          + f" vs 1/sqrt 3 = {1 / np.sqrt(3):.3f}"))
    unstable = [n for n in N_INIT if n > max(n_max.values())]
    ends = {m: {n: [out[(m, n, s)]["res"]["n_end"] for s in (SEEDS_12 if n == 12 else (0,))] for n in unstable}
            for m in MESHES}
    all_ends = [e for m in MESHES for v in ends[m].values() for e in v]
    j_frac = lambda n: (1 - (n / mid) ** 2) * (n / mid) / (2 / (3 * np.sqrt(3)))     # uniform-state j / j_dep
    checks.append(C.check("2b.5", f"start above the limit, n = {unstable[0]}..12: end state, all meshes and seeds",
                          f"n = {min(all_ends)}..{max(all_ends)}, i.e. current {j_frac(min(all_ends)):.2f}-"
                          f"{j_frac(max(all_ends)):.2f} j_dep", "8 or slightly below (plan)",
                          min(all_ends) >= max(n_max.values()) - 2,
                          "; ".join(f"{m:g} xi: " + ", ".join(f"{n}->{v}" for n, v in ends[m].items()) for m in MESHES)
                          + ".  Each slip leaves the slip point depressed and the next follows ~5 tau0 later, down to "
                          "about half the depairing current: TDGL switching/retrapping hysteresis, no heat needed"))
    dense = out[(fine, 12, 0, "dense")]["res"]
    sizes = [d for _, d in dense["steps"]]
    coarse = [d for t, d in out[(fine, 12, 0)]["res"]["steps"] if t <= 100]
    checks.append(C.check("2b.6", f"n = 12, {fine:g} xi, a frame every 2 solver steps for 100 tau0: size of every "
                          "change between accepted frames", f"{sizes}  ({dense['rejected']} of {len(dense['frames'])} "
                          "frames rejected)", "all -1 (one slip at a time, plan)", all(d == -1 for d in sizes),
                          f"same run with a frame every 10 steps: {coarse}.  In a fast burst |psi| at the slip point "
                          "stays below the gate through several slips, so only the net change is defined"))
    unbridged = []
    for k, r in out.items():
        fr, last = r["res"]["frames"], None
        for i, (t, n, fl) in enumerate(fr):
            if fl < PSI_LOOP_MIN:
                continue
            if last is not None and abs(n - fr[last][1]) > 1 and i == last + 1:
                unbridged.append((k, t, n - fr[last][1]))
            last = i
    n_multi = sum(abs(d) > 1 for r in out.values() for _, d in r["res"]["steps"])
    checks.append(C.check("2b.6b", "every change of more than one quantum spans rejected frames (all runs, both "
                          "frame rates)", f"{n_multi - len(unbridged)} of {n_multi} bridged",
                          "all (a multi-slip jump between two consecutive accepted frames would be a meter fault)",
                          not unbridged, f"unbridged: {unbridged}" if unbridged else ""))
    seed_spread = max(max(v) - min(v) for m in MESHES for v in [ends[m][12]])
    checks.append(C.check("2b.7", "mesh convergence: largest stable winding on 1.4 / 0.8 / 0.5 xi",
                          " / ".join(str(n_max[m]) for m in MESHES), "identical", len(set(n_max.values())) == 1,
                          "end states from n = 12 (seeds 0-2): " + "; ".join(f"{m:g} xi {sorted(ends[m][12])}" for m in MESHES)
                          + f" -- they differ by up to {seed_spread} between seeds at one mesh, so they are not a mesh test"))
    extra.update(stable={str(m): s for m, s in stable.items()}, current=rows,
                 runs={"_".join(map(str, k)): dict(n_end=r["res"]["n_end"], steps=r["res"]["steps"],
                                                   rejected=r["res"]["rejected"]) for k, r in out.items()})
    return checks, extra
