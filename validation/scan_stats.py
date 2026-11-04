"""
2d  Fluxoid steps in the existing flux scans (results/fluxtrap/, no new physics), and a dense-frame re-run of a few
busy runs.

  steps     between consecutive accepted frames the number of quanta in a hole should change by 1 at a time.  The
            scan CSV keeps each run's history of accepted states, so this is checked per outcome class: the
            trapped states the thesis quotes, and the resistive episodes (transient / direct / via_trap).
  gate      the plan expects each change to coincide with rejected frames (a vortex core crossing the loop).
            2a shows the gate catches a core only within ~0.5 nm of a loop point, so this is measured, from the
            saved traces (results/fluxtrap/trace_*.txt), not required.
  dense     busy runs re-simulated with 5x denser frames (save_every 8 instead of 40; the solver does not depend
            on it): the production-rate frames are every 5th dense frame, so the stored history must come back
            exactly, and the dense record shows whether its jumps of 2+ are frames too far apart (they split into
            single steps) or bursts during which the loop's |psi| stays below the gate (rejected frames between).
  mesh      the stored zipper results were made with the committed sim.make_device, which takes the hole outline
            from tdgl.geometry.ellipse.  The working tree takes it from geometry.Ellipse.outline: the same points to
            ~1e-14 nm, but the mesher returns a different triangulation.  The zipper re-runs use the committed
            construction (head_device); the size of the change, and its effect on the zipper's I_sw and three
            thresholds, are measured on the working-tree mesh.
"""
from __future__ import annotations

import re
import time
from collections import Counter
from dataclasses import replace

import numpy as np

from snspd import geometry as geo
from snspd import sim, store
from snspd.fluxoid import PSI_LOOP_MIN

from . import common as C

ID, TITLE = "2d", "fluxoid steps in the flux scans"
GEOMS = ("Long_Zipper", "Slot_400")
SCAN = "results/fluxtrap/scan_{}_conserving_T400_m1.4.csv"
TRACES = ("results/fluxtrap/trace_Long_Zipper_-68_0_conserving_T400_m1.4.txt",)
DENSE = (("Long_Zipper", 0.0, 0.0, 0.81), ("Long_Zipper", -96.0, 0.0, 0.78), ("Slot_400", -42.0, 0.0, 0.82))
EVERY, PROD_EVERY = 8, 40          # dense and production (fluxtrap.SAVE_EVERY) solver steps per frame
RESISTIVE = ("transient", "direct", "via_trap")
HITMAP = "results/hitmap/hitmap_conserving_T400_m1.4.csv"
ISW = "results/hitmap/isw_conserving_T400_m1.4.json"
MESH_SITES = ((-96.0, 0.0), (-68.0, 0.0), (0.0, 0.0))     # zipper thresholds re-run on the working-tree mesh


def head_device(g, mesh_xi):
    """sim.make_device as committed (HEAD): hole outlines from tdgl.geometry.ellipse.  Ellipse holes only."""
    import tdgl
    from tdgl.geometry import box, ellipse
    W, L = g.width_nm * 1e-3, g.length_nm * 1e-3
    holes = [tdgl.Polygon(f"hole_{i}", points=ellipse(h.a_nm * 1e-3, h.b_nm * 1e-3, center=(h.cx_nm * 1e-3, h.cy_nm * 1e-3),
                                                      points=g.hole_points)).resample(g.hole_points)
             for i, h in enumerate(g.holes)]
    terminals = [tdgl.Polygon("source", points=box(W, L / 100, center=(0, L / 2))),
                 tdgl.Polygon("drain", points=box(W, L / 100, center=(0, -L / 2)))]
    dev = tdgl.Device("snspd", layer=sim.FILM.layer(), film=tdgl.Polygon("film", points=box(W, L)), holes=holes or None,
                      terminals=terminals, probe_points=[(0, L / 3), (0, -L / 3)], length_units="um")
    dev.make_mesh(max_edge_length=sim.FILM.xi_nm * 1e-3 * mesh_xi, smooth=50)
    return dev


def frame_steps(sol):
    """Solver step of every saved frame, in the order fluxoid.frames reads them (h5 groups are save counters)."""
    import h5py
    with h5py.File(sol.path, "r") as f:
        return [int(f["data"][g].attrs["step"]) for g in sorted(f["data"].keys(), key=int)]


def history(text):
    return [(float(t), tuple(int(v) for v in s.split(",")) if s else ()) for t, s in (e.split(":") for e in text.split())]


def max_step(a, b):
    return max((abs(x - y) for x, y in zip(a, b)), default=0)


def steps(frames):
    """[(t, max per-hole |dn|, bridged by rejected frames)] for every change between accepted frames."""
    out, last = [], None
    for i, (t, n, fl, _) in enumerate(frames):
        if fl < PSI_LOOP_MIN:
            continue
        if last is not None and tuple(n) != tuple(frames[last][1]):
            out.append((t, max_step(n, frames[last][1]), i > last + 1))
        last = i
    return out


def dense_task(task):
    """fluxtrap.observed_run with frames every EVERY steps; the production frames are those at steps divisible by
    PROD_EVERY, plus the final step (pyTDGL always saves it)."""
    t_start = time.time()
    try:
        import fluxtrap
        from snspd import hotspot as hsm
        from snspd.fluxoid import frames
        g = geo.get(task["geom"])
        cfg = replace(C.CFG, save_every=EVERY)
        if all(isinstance(h, geo.Ellipse) for h in g.holes):           # made with the committed make_device
            dev = head_device(g, cfg.mesh_xi)
            hs = hsm.HotspotSolver(dev.points, dev.triangles)
        else:                                                           # slots: made with the working tree
            ctx = sim.context(g, cfg.mesh_xi)
            dev, hs = ctx["device"], ctx["hs"]
        I = task["bias"] * sim.FILM.I_dep_uA(g.free_width_min_nm())
        eps = hsm.build_eps(hs, cfg.hotspot, task["x"], task["y"], 1.0, cfg.tmax)
        with sim.solved(dev, I, eps, cfg) as sol:
            trig = sim.settle(dev, I, sol, sim.detect(sol), cfg)
            fr = frames(sol, fluxtrap.meter_for(dev))
            idx = frame_steps(sol)
            t, v = np.asarray(sol.dynamics.time), np.asarray(sol.dynamics.voltage())
        keep = [k for k, s in enumerate(idx) if s % PROD_EVERY == 0 or k == len(idx) - 1]
        sub = [(fr[k][0], fr[k][1], fr[k][2], max(f[3] for f in fr[(p + 1 if j else 0):k + 1]))
               for j, (p, k) in enumerate(zip([0] + keep[:-1], keep))]        # max |V| over the coarse interval
        rec = fluxtrap.classify(trig, sub, t, v, hsm.HOTSPOT.t0_tau0)
        return dict(ok=True, task=task, res=dict(dense=steps(fr), prod=steps(sub), history_sub=rec["history"],
                                                 outcome=rec["outcome"], n_frames=len(fr)),
                    wall_s=time.time() - t_start)
    except Exception as exc:
        return dict(ok=False, task=task, error=f"{type(exc).__name__}: {exc}", wall_s=time.time() - t_start)


def trace_frames(path):
    """{bias: [(t, n, floor, vmax)]} from a fluxtrap.py trace text file."""
    out, cur = {}, None
    for line in open(path):
        m = re.match(r"bias ([\d.]+):", line)
        if m:
            cur = out.setdefault(float(m.group(1)), [])
            continue
        m = re.match(r"\s+t =\s+([\d.]+)\s+\[([-\d ]*)\]\s+([\d.]+)\s+([\d.]+)", line)
        if m and cur is not None:
            cur.append((float(m.group(1)), tuple(int(v) for v in m.group(2).split()), float(m.group(3)),
                        float(m.group(4))))
    return out


def run(a):
    checks, extra = [], {}
    rows = {g: store.load_rows(SCAN.format(g)) for g in GEOMS}

    # steps between accepted states in the stored histories
    for k, g in enumerate(GEOMS):
        cnt = {c: Counter() for c in ("trapped",) + RESISTIVE}
        for r in rows[g]:
            if r["outcome"] in cnt:
                h = history(r["history"])
                for (_, a_), (_, b) in zip(h, h[1:]):
                    cnt[r["outcome"]][min(max_step(a_, b), 3)] += 1
        tr = cnt["trapped"]
        res = Counter()
        for c in RESISTIVE:
            res.update(cnt[c])
        extra[g] = {c: dict(v) for c, v in cnt.items()}
        checks.append(C.check(f"2d.{2 * k + 1}", f"{g} scan, runs that end trapped (no trigger): changes of one quantum "
                              "per hole between accepted states", f"{tr[1]} of {sum(tr.values())}",
                              "all", sum(tr.values()) > 0 and tr[1] == sum(tr.values())))
        checks.append(C.check(f"2d.{2 * k + 2}", f"{g} scan, resistive runs ({'/'.join(RESISTIVE)}): changes of "
                              "2 or more quanta in one hole", f"{res[2] + res[3]} of {sum(res.values())} "
                              f"({100 * (res[2] + res[3]) / max(1, sum(res.values())):.1f} %)",
                              "0 if frames resolve every vortex", res[2] + res[3] == 0,
                              "per outcome: " + "; ".join(f"{c} {dict(cnt[c])}" for c in RESISTIVE)))

    # rejected frames at the changes, from the saved traces
    n_chg = n_rej = 0
    floors = []
    for path in TRACES:
        for b, fr in trace_frames(path).items():
            for t, size, bridged in steps(fr):
                n_chg += 1
                n_rej += bridged
                i = next(j for j, f in enumerate(fr) if f[0] == t)
                floors.append(min(f[2] for f in fr[max(0, i - 1):i + 1]))
    checks.append(C.check("2d.5", "trace of the zipper photon at (-68, 0), 4 biases: changes with a rejected frame "
                          "between the two accepted states", f"{n_rej} of {n_chg}",
                          "all (plan); not required, see 2a.5-2a.6", note=f"min |psi| on the loops at the frame before "
                          f"or at each change: median {np.median(floors):.2f}, min {min(floors):.2f}" if floors else ""))

    # dense re-runs, and the zipper's I_sw and thresholds on the working-tree mesh (production task code)
    zg = geo.LONG_ZIPPER
    isw_z = store.load_json(ISW, {})[zg.key]["I_sw"]
    tasks = [dict(kind="dense", geom=g, x=x, y=y, bias=b) for g, x, y, b in DENSE]
    tasks += [dict(kind="isw", geom=zg, cfg=C.CFG)]
    tasks += [dict(kind="hit", geom=zg, cfg=C.CFG, x_nm=x, y_nm=y, f=1.0, I_sw=isw_z) for x, y in MESH_SITES]
    out, now = {}, {}

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"   FAILED {t['kind']}: {r['error']}")
            return
        if t["kind"] == "dense":
            out[(t["geom"], t["x"], t["y"], t["bias"])] = r["res"]
            sim.log(f"   {t['geom']} ({t['x']:g},{t['y']:g}) @ {t['bias']}: {r['res']['n_frames']} frames "
                    f"({r['wall_s']:.0f} s)")
        else:
            now["isw" if t["kind"] == "isw" else (t["x_nm"], t["y_nm"])] = r["res"]
            sim.log(f"   working-tree zipper mesh: {t['kind']} done ({r['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result, fn=task_fn)
    same, lines, unbr, multi_d, multi_p, split = True, [], 0, 0, 0, 0
    for key, q in out.items():
        g, x, y, b = key
        stored = next((r for r in rows[g] if (float(r["x_nm"]), float(r["y_nm"]), round(float(r["bias"]), 4)) ==
                       (x, y, b)), None)
        same &= stored is not None and stored["history"] == q["history_sub"]
        md = [s for s in q["dense"] if s[1] > 1]
        mp = [s for s in q["prod"] if s[1] > 1]
        multi_d += len(md)
        multi_p += len(mp)
        unbr += sum(not s[2] for s in md)
        split += len(mp) - len(md)
        lines.append(f"{g} ({x:g},{y:g}) @ {b}: {len(mp)} jumps of 2+ at the production rate, {len(md)} at 5x "
                     f"({sum(s[2] for s in md)} bridged by rejected frames)")
    checks.append(C.check("2d.6", "dense re-runs: every 5th frame reproduces the stored history", "yes" if same else "no",
                          "identical (deterministic solver; frames do not change the dynamics)", same))
    checks.append(C.check("2d.7", "dense re-runs: jumps of 2+ quanta between accepted frames, production rate vs 5x",
                          f"{multi_p} -> {multi_d}; at 5x, {multi_d - unbr} of {multi_d} span rejected frames",
                          "every remaining jump spans rejected frames (a burst, winding undefined in between)",
                          unbr == 0, "; ".join(lines)))

    n_head, n_now = len(head_device(zg, C.CFG.mesh_xi).points), len(sim.context(zg, C.CFG.mesh_xi)["device"].points)
    checks.append(C.check("2d.8", "zipper mesh (1.4 xi): committed sim.make_device vs the working tree",
                          f"{n_head} vs {n_now} sites", "identical (the outlines agree to ~1e-14 nm)", n_head == n_now,
                          "stored zipper results (hitmap, flux scan, trapped states) reproduce bit for bit only with "
                          "the committed construction (2d.6); zipper runs made now use another triangulation"))
    tol = sim.scales(zg)[0]
    ref = {(float(r["x_nm"]), float(r["y_nm"])): float(r["I_det_uA"]) for r in store.load_rows(HITMAP)
           if r["Geometry"] == zg.key and float(r["energy_f"]) == 1.0}
    diffs = []
    if "isw" in now:
        diffs.append(("I_sw", now["isw"]["I_sw"], isw_z))
    diffs += [(f"({x:g},{y:g})", now[(x, y)]["I_det"], ref[(x, y)]) for x, y in MESH_SITES if (x, y) in now]
    if diffs:
        worst = max(abs(a_ - b) for _, a_, b in diffs) / tol
        checks.append(C.check("2d.9", "zipper on the working-tree mesh vs stored: I_sw and the f = 1 thresholds at "
                              + ", ".join(f"({x:g},{y:g})" for x, y in MESH_SITES),
                              f"largest change {worst:.1f} x the bisection tolerance ({tol:.2f} uA)",
                              "<= 2 x tol (same resolution, only the triangulation differs)", worst <= 2.0,
                              "; ".join(f"{n}: {a_:.2f} vs {b:.2f} uA ({(a_ - b) / sim.FILM.I_dep_uA(zg.free_width_min_nm()):+.4f}"
                                        " I_dep)" for n, a_, b in diffs)))
    extra.update(mesh_sites=dict(committed=n_head, working_tree=n_now),
                 working_tree_mesh={str(k): v for k, v in now.items()})
    return checks, extra


def task_fn(task):
    return dense_task(task) if task["kind"] == "dense" else sim.task_fn(task)
