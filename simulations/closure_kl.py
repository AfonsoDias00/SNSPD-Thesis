#!/usr/bin/env python3
"""
closure_kl.py -- the depairing current of a dirty superconductor (Kupriyanov & Lukichev 1980) against the eps(T)
closures of the heat model (13 Nov 2026).

Usadel equation for a uniform current state (phase gradient q), theta parametrisation, units k_B = hbar = T_c = 1:

    omega sin(theta) + Gamma sin(theta) cos(theta) = Delta cos(theta),      Gamma = D q^2 / 2  (pair breaking)
    gap:      ln(1/T) = 2 pi T sum_{omega > 0} ( 1/omega - sin(theta) / Delta )
    current:  j ~ sqrt(Gamma) * 2 pi T sum_{omega > 0} sin^2(theta)                (same prefactor at every T)

j_dep(T) is the maximum of j over Gamma.  The closures give j_dep(T) / j_dep(T0) = eps(T)^(3/2) with
    two-fluid  eps = (1 - t^2) / (1 - t0^2)      (snspd/thermal.py default)
    linear GL  eps = (1 - t)   / (1 - t0)
and the KL-matched closure is eps_KL(t) = (j_KL(t) / j_KL(t0))^(2/3).

  python closure_kl.py        writes results/closure/kl.csv, summary.txt, kl_closure.png
"""
from __future__ import annotations

import os

import numpy as np
from scipy.optimize import brentq, minimize_scalar

OUT = "results/closure"
T0 = 0.3
OMEGA_MAX = 400.0      # Matsubara cut-off [T_c]; the summands fall as 1/omega^2, the remainder is added analytically


def theta(omega, delta, gam):
    """Solve omega s + gam s c = delta c for theta in (0, pi/2), vectorised bisection on the monotone residual."""
    lo, hi = np.zeros_like(omega), np.full_like(omega, np.pi / 2)
    for _ in range(60):
        m = 0.5 * (lo + hi)
        f = omega * np.sin(m) + gam * np.sin(m) * np.cos(m) - delta * np.cos(m)
        hi = np.where(f > 0, m, hi)
        lo = np.where(f > 0, lo, m)
    return 0.5 * (lo + hi)


def matsubara(T):
    n = np.arange(int(OMEGA_MAX / (2 * np.pi * T)) + 1)
    return np.pi * T * (2 * n + 1)


def gap(T, gam):
    """Self-consistent Delta at pair breaking gam (0 if superconductivity is gone)."""
    w = matsubara(T)
    wmax = w[-1] + np.pi * T

    def F(d):
        s = np.sin(theta(w, d, gam))
        # remainder for omega > wmax: 1/w - s/d ~ gam / w^2 + O(d^2 / w^3) -> 2 pi T sum ~ gam / wmax
        return 2 * np.pi * T * np.sum(1 / w - s / d) + gam / wmax - np.log(1 / T)

    if F(1e-6) >= 0:      # F rises with Delta; no root above 0 means no superconductivity
        return 0.0
    return brentq(F, 1e-6, 3.0, xtol=1e-10)


def current(T, gam):
    d = gap(T, gam)
    if d == 0.0:
        return 0.0, 0.0
    w = matsubara(T)
    s = np.sin(theta(w, d, gam))
    return np.sqrt(gam) * 2 * np.pi * T * np.sum(s ** 2), d


def jdep(T):
    """max over Gamma of j; returns (j_dep, Gamma*, Delta at j_dep, Delta at zero current)."""
    d0 = gap(T, 0.0)
    # j(Gamma) is zero beyond the pair-breaking limit, which traps a bounded minimiser: scan, then refine
    gs = np.linspace(1e-4, 1.0, 81) * d0
    js = np.array([current(T, g)[0] for g in gs])
    k = int(np.argmax(js))
    lo, hi = gs[max(k - 1, 0)], gs[min(k + 1, len(gs) - 1)]
    r = minimize_scalar(lambda g: -current(T, g)[0], bounds=(lo, hi), method="bounded", options=dict(xatol=1e-8))
    j, d = current(T, r.x)
    return j, r.x, d, d0


def main():
    os.makedirs(OUT, exist_ok=True)
    ts = np.concatenate([np.linspace(0.05, 0.9, 18), np.linspace(0.92, 0.99, 8)])
    ts = np.unique(np.round(np.append(ts, T0), 4))
    rows = [(t, *jdep(t)) for t in ts]
    t = np.array([r[0] for r in rows])
    j = np.array([r[1] for r in rows])
    d0 = np.array([r[4] for r in rows])
    j0 = j[np.argmin(abs(t - T0))]
    jr = j / j0
    two = ((1 - t ** 2) / (1 - T0 ** 2)) ** 1.5
    lin = ((1 - t) / (1 - T0)) ** 1.5
    eps_kl = jr ** (2 / 3)
    eps_2f = (1 - t ** 2) / (1 - T0 ** 2)
    eps_lin = (1 - t) / (1 - T0)
    with open(f"{OUT}/kl.csv", "w") as fh:
        fh.write("t,jdep_KL,Gamma_star,Delta_at_jdep,Delta_BCS,jKL_rel,j_twofluid_rel,j_linear_rel,eps_KL,eps_twofluid,"
                 "eps_linear\n")
        for k, r in enumerate(rows):
            fh.write(",".join(f"{v:.6g}" for v in (*r, jr[k], two[k], lin[k], eps_kl[k], eps_2f[k], eps_lin[k])) + "\n")
    L = ["KUPRIYANOV-LUKICHEV DEPAIRING CURRENT vs the eps(T) closures (all relative to the bath, t0 = 0.3)", "",
         f"checks: Delta_BCS(0.05) = {d0[0]:.4f} (BCS 1.764); Delta_BCS(0.99) = {d0[-1]:.4f} (GL limit 3.063 sqrt(1-t) "
         f"= {3.063 * np.sqrt(1 - t[-1]):.4f}); j_KL(0.05)/j_KL(t0) = {jr[0]:.4f}", "",
         f"  {'t':>5}  {'j_KL':>7} {'two-fl':>7} {'linear':>7}   {'eps_KL':>7} {'eps_2f':>7} {'eps_lin':>7}   "
         f"{'2f/KL':>6} {'lin/KL':>6}"]
    for k in range(len(t)):
        L.append(f"  {t[k]:5.2f}  {jr[k]:7.4f} {two[k]:7.4f} {lin[k]:7.4f}   {eps_kl[k]:7.4f} {eps_2f[k]:7.4f} "
                 f"{eps_lin[k]:7.4f}   {two[k] / jr[k]:6.3f} {lin[k] / jr[k]:6.3f}")
    # near T_c: GL limit j ~ (1 - t)^(3/2); the KL ratio to that should tend to a constant
    nt = t > 0.9
    L += ["", "near T_c, j_KL / (1 - t)^1.5: " + "  ".join(f"{a:.2f}:{b:.4f}" for a, b in
                                                          zip(t[nt], jr[nt] / (1 - t[nt]) ** 1.5))]
    text = "\n".join(L)
    print(text)
    with open(f"{OUT}/summary.txt", "w") as fh:
        fh.write(text + "\n")
    plot(t, jr, two, lin, eps_kl, eps_2f, eps_lin)


def plot(t, jr, two, lin, eps_kl, eps_2f, eps_lin):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tt = np.linspace(0.0, 1.0, 201)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    a = ax[0]
    a.plot(t, jr, "o-", color="#111827", lw=2, ms=4, label="Kupriyanov-Lukichev (Usadel)")
    a.plot(tt, np.clip((1 - tt ** 2) / (1 - T0 ** 2), 0, None) ** 1.5, color="#c2410c", lw=2,
           label="two-fluid  (1-t²)^1.5  [default]")
    a.plot(tt, np.clip((1 - tt) / (1 - T0), 0, None) ** 1.5, color="#1f5fa6", lw=2, ls="--",
           label="linear GL  (1-t)^1.5")
    a.axvspan(0.5, 0.9, color="0.5", alpha=0.12, lw=0)
    a.text(0.7, 1.12, "hotspot halo", ha="center", fontsize=8, color="0.35")
    a.set_xlabel("electron temperature  T_e / T_c")
    a.set_ylabel("j_dep(T_e) / j_dep(bath 0.3 T_c)")
    a.set_xlim(0, 1)
    a.set_ylim(0, 1.2)
    a.grid(alpha=0.25)
    a.legend(fontsize=8, frameon=False, loc="lower left")
    a = ax[1]
    a.plot(t, eps_2f / eps_kl, color="#c2410c", lw=2, label="two-fluid / KL")
    a.plot(t, eps_lin / eps_kl, color="#1f5fa6", lw=2, ls="--", label="linear GL / KL")
    a.axhline(1, color="0.4", lw=1)
    a.axvspan(0.5, 0.9, color="0.5", alpha=0.12, lw=0)
    a.set_xlabel("electron temperature  T_e / T_c")
    a.set_ylabel("eps of the closure / eps matched to KL")
    a.set_xlim(0, 1)
    a.grid(alpha=0.25)
    a.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/kl_closure.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
