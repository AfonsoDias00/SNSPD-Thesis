#!/usr/bin/env python3
"""
sensitivity_scan.py -- how much do the hotspot model's weakest assumptions move the thresholds?  (13 Nov 2026)

At 1550 nm, 2.0 xi, circuit 750 squares: the plain wire (8 sites) and the double column (23-site medium grid of
energy_scan.py), once per model variant.  Each variant writes its own file thr_n750_w1550_m2_<label>.csv; the
reference is the default model (thr_n750_w1550_m2.csv).  All tasks share one pool, variant by variant, so a cut
run leaves complete variants.

  variant    override                              question
  2fNL       latent=False                          reference for linNL: the two-fluid closure without latent heat
  linNL      eps_model=linear, latent=False        closure: linear GL eps(T) instead of two-fluid (latent heat
                                                    needs the two-fluid entropy, so both are run without it)
  tth2.8     tau_th_ps=2.8                         first picosecond: slower thermalisation
  tth0.7     tau_th_ps=0.7                         first picosecond: faster thermalisation
  r015       r0_nm=15                              first picosecond: wider initial deposit

The deposited fraction is not rerun: the model depends only on the deposited energy (abs0.5 check, 12 Nov), so
absorbed = 0.75 (Simon et al.) is the ln-E interpolation of the 1550 / 2500 nm runs.

  python sensitivity_scan.py run --workers 11
  python sensitivity_scan.py plan
  python sensitivity_scan.py report     (results/energy/summary_sensitivity.txt, sensitivity.png)
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from energy_scan import build_tasks
from snspd import sim, store

VARIANTS = (  # label, photon overrides, thermal overrides
    ("2fNL", "", "latent=False"),
    ("linNL", "", "eps_model=linear,latent=False"),
    ("tth2.8", "tau_th_ps=2.8", ""),
    ("tth0.7", "tau_th_ps=0.7", ""),
    ("r015", "r0_nm=15", ""),
)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["run", "plan", "report"])
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--mesh-xi", type=float, default=2.0)
    a = ap.parse_args(argv)
    if a.command == "report":
        return report()
    tasks, skipped = [], 0
    for label, photon, thermal in VARIANTS:
        t, s = build_tasks(a.mesh_xi, (1550.0,), photon, label, rest_1550=False, thermal=thermal)
        for x in t:
            x["label"] = label
        tasks += t
        skipped += s
        sim.log(f"{label}: {len(t)} to run, {s} done")
        if t:
            for line in t[0]["cfg"].describe():
                sim.log("    " + line)
    sim.log(f"{len(tasks)} thresholds to run, {skipped} already done")
    if a.command == "plan":
        return

    def on_result(r):
        t = r["task"]
        if not r["ok"]:
            sim.log(f"FAILED {t['label']} {t['geom'].key} {t['site']}: {r['error']}\n{r.get('tb', '')}")
            return
        row = store.threshold_row(t["geom"].key, t["site"], t["tag"], r["res"], r["wall_s"])
        store.append_row(t["path"], row)
        sim.log(f"  {t['label']:<7} {t['geom'].key:<11} {t['site']}: {row['i_det']:.3f} {row['flag']}"
                f"  ({row['n_runs']} runs, {row['wall_s']:.0f} s)")

    sim.execute(tasks, a.workers, on_result)


def report():
    import numpy as np
    from energy_analysis import dark_of, quant, samples
    from energy_scan import MEDIUM
    from snspd import analysis as an
    from snspd import geometry as geo
    ib = np.linspace(0.5, 1.0, 501)
    base = "results/thresholds/thr_n750_w1550_m2"
    ref = store.load_rows(base + ".csv")
    labels = ("default",) + tuple(v[0] for v in VARIANTS)
    L = ["SENSITIVITY OF THE THRESHOLDS AT 1550 nm (2.0 xi, circuit 750 squares; shifts against the default model)", ""]
    curves = {}
    for key in ("Plain_Wire", "Double_400"):
        g = geo.get(key)
        d = dark_of(ref, key)
        R = samples(ref, key, d)
        sites = [tuple(map(float, s)) for s in (g.hit_grid("full") if key == "Plain_Wire" else MEDIUM)]
        L.append(f"{g.label} (dark switching {d:.3f}, {len(sites)} sites; channels y <= 180 nm, bridge rows y > 180)")
        L.append(f"  {'variant':<8} {'mean':>7} {'chan.':>7} {'bridge':>7} {'max|.|':>7}   i_10   i_50   i_90   i_99"
                 f"  eta(0.80)")
        curves[key] = {}
        for lab in labels:
            V = R if lab == "default" else samples(store.load_rows(f"{base}_{lab}.csv"), key, d)
            ss = [s for s in sites if s in V]
            if len(ss) < len(sites):
                L.append(f"  {lab:<8} ({len(ss)} of {len(sites)} sites)")
                continue
            dif = np.array([V[s] - R[s] for s in ss])
            ch = [V[s] - R[s] for s in ss if s[1] <= 180]
            br = [V[s] - R[s] for s in ss if s[1] > 180]
            _, _, z = an.unit_map(g, [(x, y, V[(x, y)]) for x, y in ss])
            e = an.eta_curve(np.sort(z[np.isfinite(z)]), ib)
            curves[key][lab] = e
            f = lambda v: f"{np.mean(v):+7.3f}" if v else " " * 7
            L.append(f"  {lab:<8} {dif.mean():+7.3f} {f(ch)} {f(br)} {np.abs(dif).max():7.3f}   "
                     + "  ".join(f"{quant(ib, e, q):.3f}" for q in (10, 50, 90, 99))
                     + f"  {np.interp(0.80, ib, e):6.1f} %")
        L.append("")
    L.append("plain minus double column (negative: the plain wire is ahead)")
    L.append(f"  {'variant':<8} {'d i_50':>7} {'d i_90':>7} {'d i_99':>7}")
    for lab in labels:
        if lab in curves["Plain_Wire"] and lab in curves["Double_400"]:
            p, q = curves["Plain_Wire"][lab], curves["Double_400"][lab]
            L.append(f"  {lab:<8} " + " ".join(f"{quant(ib, p, x) - quant(ib, q, x):+7.3f}" for x in (50, 90, 99)))
    text = "\n".join(L)
    print(text)
    with open("results/energy/summary_sensitivity.txt", "w") as fh:
        fh.write(text + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"default": "#111827", "2fNL": "#9ca3af", "linNL": "#1f5fa6", "tth2.8": "#c2410c", "tth0.7": "#f59e0b",
           "r015": "#15803d"}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    for ax, key in zip(axes, ("Plain_Wire", "Double_400")):
        for lab, e in curves[key].items():
            ax.plot(ib, e, color=col[lab], lw=2.4 if lab == "default" else 1.6, label=lab)
        ax.set_title(geo.get(key).label + ", 1550 nm", fontsize=10)
        ax.set_xlim(0.55, 1.0)
        ax.set_xlabel("bias  [I_dep of the narrowest free width]")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("share of the cell that counts  [%]")
    axes[1].legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig("results/energy/sensitivity.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
