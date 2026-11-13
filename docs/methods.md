# Method: numerics, protocol and the checks behind them

What the simulation solves is in [model.md](model.md). This note keeps the numerical and protocol decisions and the
evidence for them. Most of the evidence is from the earlier phases (16 Oct – 9 Nov 2026: hole arrays, then slot
layouts), whose code and results are on the branches `repo-cleanup` (commit `00f9b0d`) and `slot-geometries`.
Numbers from those phases are quoted only as evidence for a method choice; the scripts named with them
(`mesh_check.py`, `dt_check.py`) are on `slot-geometries`.

## 1 · Solver and film

**Solver.** pyTDGL 0.9 [1]: generalised TDGL (Kramer–Watts-Tobin [2, 3]) for ψ and the potential μ on a 2D
finite-volume mesh.

    u/√(1+γ²|ψ|²) (∂ₜ + iμ + γ²/2 ∂ₜ|ψ|²) ψ = (ε(r,t) − |ψ|²) ψ + ∇²ψ,   ∇²μ = ∇·Jₛ

- Screening off, no applied field (w ≪ Pearl length Λ = 2λ²/d ≈ 80 µm).
- Current enters through two film-width terminals (ψ = 0). Film edges and slot walls are insulating.
- Temperature enters only through ε(r, t) = (1 − t_e²)/(1 − t0²) ([model.md](model.md) §4).
- The voltage is measured between probes at (0, ±probe_y): 75 nm inside each terminal for the double column.

**Film: 9 nm NbTiN** (KTH process; measured values from Sidorova et al. [39, 40]).

| Quantity | Value | Basis |
| --- | --- | --- |
| ξ | 5 nm | pyTDGL's length unit; measured ξ_GL(0) = 4.6 nm [39] |
| λ | 600 nm | fit to 9 nm NbTiN hairpins [13] |
| d, Λ | 9 nm, 80 µm | w/Λ ≤ 0.003 |
| R_sq | 381.4 Ω | measured [39, 40] |
| τ0 | 0.132 ps | L_sq/R_sq with L_sq = μ0λ²/d = 50.3 pH. (ħ/πΔ0 gives 0.15–0.16 ps.) |
| u | 5.79 | dirty limit, pyTDGL default |
| γ | 10 | 2τ_EΔ0/ħ ≈ 12–18 for the measured electron-phonon time; pyTDGL's default and Simon et al.'s value [27] |
| j_dep, I_dep | 0.504 µA/nm × narrowest free width | GL depairing; **the bias unit throughout** |

- **λ only relabels units.** With screening off it enters only the current, voltage and time units. Normalised
  results do not depend on it. The heat model's coupling is fixed by thermodynamics, not by λ ([model.md](model.md) §4).
- **GL at T ≪ T_c is formally outside its range.** Geometry, current distribution, vortex paths and thresholds in
  I_dep units are trustworthy. Absolute currents and times are not.
- **Quote everything in I_dep of the narrowest free cross-section.** The simulated dark switching current drifts
  with the window, γ and the mesh, so ratios to it would mix those effects into geometry comparisons.

## 2 · Measuring detection

- **One event** (`sim.run_event`): settled start, one photon, electro-thermal TDGL with the circuit.
- **A count:** the load carries ≥ 20 % of the bias current at some time. An event that is resistive but below that
  level at 400 τ0 runs on, to 800 τ0 at most, until it counts or is quiet (< 0.05 phase slips in 25 τ0 and no site
  above T_c).
- **Photon threshold i_det(x, y):** bisection at each site in [0.50, 1.00] I_dep, to 0.005 (7 events).
- **Dark switching:** the same bisection without a photon. A bias with no steady dark state counts and is flagged
  `dark`, so a photon threshold can never exceed it.
- **Hard-step efficiency η(I):** area fraction of the unit cell with i_det < I, interpolated between sites and
  completed by the geometry's symmetry. i_10, i_50, i_95 are where η reaches 10, 50, 95 %; i_max is 100 %.
- **Go / no-go:** the slotted device must beat its plain reference of the same width by ≥ 0.03 I_dep at the plateau
  (i_95). An onset (i_50) gain alone was never robust in the earlier phases: it depended on the criterion, window,
  mesh and hotspot.

**Trapped flux** (`snspd/fluxoid.py`): the winding of arg ψ along a loop around each slot, counted only while
min |ψ| ≥ 0.05 on the loop.

## 3 · Protocol at γ = 10: the settled start

At γ = 10, |ψ| relaxes ~6× slower (u√(1+γ²|ψ|²)). After the bias step a strip needs ~200 τ0 to settle, so a photon
at 20 τ0 would meet a condensate that is still too strong. The steady state is a stationary GL solution and does not
depend on γ, so each run starts from a state prepared by a γ = 1 dark run (`sim.settled_state`):
- 100 τ0 at the run's current, continued up to 3 × 100 τ0 until the last 25 τ0 has ≤ 0.02 phase slips and
  |V_end| ≤ 0.01.
- Checked: under γ = 10 dynamics the settled state stays put (<|ψ|²> unchanged to 1e-4, |V| ≤ 0.011 over 200 τ0).
- **Between the γ = 1 and γ = 10 dark switching currents** (≈ 0.96–0.98 I_dep) the γ = 1 bias step overshoots into
  switching. The settling then retries with the device's own γ = 10 dynamics.
- **If neither settles, there is no steady state** and the event is a dark count. (The earlier code started such
  runs from the uniform state; their bias-on transient then faked a detection at t = 55 τ0.)
- Judging stationarity on the second half of a single 100 τ0 run was too early for 2 µm devices above ~0.92 I_dep
  (their probe phase is still building up then). The last quarter of the latest stretch is used.

## 4 · Numerics

| Choice | Evidence |
| --- | --- |
| Mesh: max edge 1.4 ξ = 7 nm (mean edge ≈ 4 nm) | 1.4 vs 0.8 ξ: \|Δi_det\| ≤ 0.007 for photon thresholds in strips. **Not converged for flux crossing narrow bridges:** a trapped-flux switching current moved 0.06 between meshes in 12 nm bridges. The 40 nm features here are 8 ξ wide |
| Uniform 2.0 ξ for screening only | double column, earlier model: 1.4 − 2.0 ξ at the critical sites +0.011 on average, at most +0.016; channel sites ≤ 0.011; 6–8× faster. **Not yet repeated with the new model** |
| No graded meshes | run time scales ~h⁻⁴ (stable time step ∝ h_min²). Coarsening near a hole rim dropped dark I_sw from 0.97 to 0.83: **slot rims, bridges and contacts must stay fine** |
| Pinned meshes (`SNSPD_MESH_DIR`, `python -m snspd.sim pin-meshes meshes 1.4`) | holed geometries triangulate differently on other GEOS/shapely versions; 0.04 I_dep of triangulation noise was seen for flux-mediated quantities |
| **Solver time step: dt_max = 0.025** | dt_max 0.2 leaves a local instability at single small cells (a site with \|ψ\| above its neighbours and a stationary \|∇μ\|² 50× the 99.9th percentile, in the dark). Any heat model turns it into a permanent hot spot that seeds switching away from the photon: heated thresholds moved by up to +0.07. 0.05 removes it on 2.0 ξ meshes, 0.025 on the 1.4 ξ double column |
| Heat on pyTDGL's own Voronoi cells | energy sums are exact: audit of one photon closes to 4e-14 (validation 3a.1) |
| Heat: conduction explicit, exchange linearly implicit | two-temperature relaxation within 0.07 % of a stiff ODE reference at dt = 0.025 (3a.3–4); D recovered to 0.02 % (3a.5) |
| Contact heat sinks 50 nm | 25 nm left the electrons at 0.41 T_c beside the terminals in the dark ([model.md](model.md) §6) |
| Window 400 τ0 (to 800 if undecided) | earlier model, γ = 10: 800 τ0 changed thresholds by ≤ 0.003 at 4 sites, although threshold runs fire late (330–395 τ0). **Not yet repeated with the new model**; one bridge-row event at 0.90 I_dep started switching only at 430 τ0 |
| Shortened strip | 150 nm of plain lead beyond the outermost slot reproduces a 2 µm strip within 0.011 I_dep at 1.4 ξ (earlier model) |
| Photon sites away from the contacts | ψ = 0 terminals nucleate dark switching 20–50 nm from the contacts |

## 5 · What the model leaves out

- **No thermal noise.** The dark switching current is an instability current, not a dark-count rate, so the model
  cannot show dark counts or their suppression. Barrier estimates need GL saddle points [21].
- **Ideal strips.** Real 9 nm NbTiN strips switch at 0.63 I_dep (3 K) [40]; 0.43–0.74 for NbN/WSi devices [23].
  Every slot adds rim length, where fabrication defects and constrictions occur.
- **Near-T_c dynamics at 0.3 T_c**, and the other assumptions listed in [model.md](model.md) §7.

## 6 · Lessons from the earlier geometries

- Every hole layout of the first phase failed the plateau test; holes add an insensitive tail at bridges and hole
  tips. Gains appeared only at onset and depended on definitions.
- Narrow bridges (≲ 4–5 ξ) are under-resolved and behave as 1D links.
- Trapped flux below threshold makes a slotted device history-dependent.
- Always run the plain strip of the same width as the reference, with the same model and mesh.

## References

Numbers follow [references.md](references.md).

1. L. Bishop-Van Horn, pyTDGL, *Comput. Phys. Commun.* **291**, 108799 (2023), arXiv:2302.03812.
2. L. Kramer and R. J. Watts-Tobin, *Phys. Rev. Lett.* **40**, 1041 (1978).
3. R. J. Watts-Tobin, Y. Krähenbühl and L. Kramer, *J. Low Temp. Phys.* **42**, 459 (1981).
8. A. N. Zotova and D. Y. Vodolazov, *Phys. Rev. B* **85**, 024509 (2012), arXiv:1112.3790.
13. M. Jönsson *et al.*, *Phys. Rev. Applied* **17**, 064046 (2022), arXiv:2112.05443.
21. D. Y. Vodolazov, *Phys. Rev. B* **85**, 174507 (2012), arXiv:1204.0365.
23. S. Frasca *et al.*, *Phys. Rev. B* **100**, 054520 (2019), arXiv:1904.08596.
27. Simon *et al.*, *Phys. Rev. B* **112**, 174512 (2025), arXiv:2501.13791.
39. M. Sidorova *et al.*, *Phys. Rev. B* **104**, 184514 (2021), arXiv:2108.04717.
40. M. Sidorova *et al.*, *Supercond. Sci. Technol.* (2022), DOI 10.1088/1361-6668/ac8454, arXiv:2204.01423.
