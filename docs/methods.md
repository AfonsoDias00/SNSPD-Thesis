# Method: what the simulation does and why

This note keeps the method decisions of the first phase of the study (16 Oct – 6 Nov 2026) and the evidence behind
them. The first phase studied hole arrays in 200 nm strips (ellipse "zipper", 30 nm slot columns, narrow strips). Its
geometries, results and full report were removed on the `slot-geometries` branch. They are all in commit `00f9b0d`
("Snapshot before cleanup") on `repo-cleanup`; `git show 00f9b0d:docs/report.md` prints the old report.

Numbers below that come from old geometries are quoted only as evidence for a method choice.

## 1 · Model

**Solver.** pyTDGL 0.9 [1]: generalized TDGL (Kramer–Watts-Tobin [2, 3]) for ψ and the potential μ on a 2D
finite-volume mesh.

    u/√(1+γ²|ψ|²) (∂ₜ + iμ + γ²/2 ∂ₜ|ψ|²) ψ = (ε(r,t) − |ψ|²) ψ + (∇ − iA)² ψ,   ∇²μ = ∇·Jₛ

- Screening off, no applied field (w ≪ Pearl length Λ = 2λ²/d ≈ 80 µm).
- Current enters through two film-width terminals (ψ = 0). Film edges and hole walls are insulating.
- The photon enters only through ε(r, t) ≤ 1.
- The voltage is measured between probes at (0, ±L/3), or at the geometry's `probe_y_nm`.

**Film: 9 nm NbTiN** (KTH process; Sidorova et al. [39, 40] for the measured values).

| Quantity | Value | Basis |
| --- | --- | --- |
| ξ | 5 nm | sets the mesh and hotspot scale |
| λ | 600 nm | fit to 9 nm NbTiN hairpins [13] |
| d, Λ | 9 nm, 80 µm | w/Λ ≤ 0.003 |
| u | 5.79 | dirty limit, pyTDGL default |
| γ | **10** (preset `nbtin`) | physical estimate γ = 2τ_EΔ0/ħ ≈ 12–18 (τ_EP(T_c) = 2.9–4.0 ps [39, 40]); pyTDGL's default and Simon et al.'s value [27]. γ = 1 was used until 6 Nov (preset `""`) |
| τ0 | ≈ 0.13–0.16 ps | τ0 = μ0λ²/ρ_n ≈ ħ/(πΔ0); both routes agree for these films |
| D | 7.1 nm²/τ0 (0.47 cm²/s) | measured [39] |
| j_dep, I_dep | 0.504 µA/nm × narrowest free width | GL depairing; **the bias unit throughout** |

- **λ and ρ_n only relabel units.** With screening off they enter only through K0 and τ0. Normalised results do not
  depend on them; absolute currents scale as 1/λ².
- **GL at T ≪ T_c is formally outside its range.** Geometry, current distribution, vortex paths and thresholds in
  I_dep units are trustworthy. Absolute currents and times are not.
- **Units: quote everything in I_dep of the narrowest free cross-section.** I_sw drifts with the window and with γ, so
  I_det/I_sw would mix those effects into geometry comparisons.

## 2 · Photon: the kinetic hotspot

The energy diffuses inside the film with reflecting walls (P1 finite elements, Crank–Nicolson; energy conserved to
1e-8 on every mesh tested). This replaced a free-space Gaussian, which leaked energy at exactly the edges and hole
walls where geometries differ.

    ∂n/∂t = D∇²n,   ∂n/∂ν = 0 on edges and hole walls,   n(r, t0) = Gaussian, r0 = 14 nm (2.8 ξ), ∫n dA = π r0²
    e = f n e^{−(t−t0)/t_esc},   ds/dt = (s_eq(e) − s)/t_rise,   s_eq = s_max tanh(|A| e / s_max),   ε = 1 − s

After Simon et al. [27] (ab initio cascade → gap suppression → pyTDGL):

| Parameter | Value (preset `nbtin`) | Source |
| --- | --- | --- |
| t_rise (thermalisation τ_th) | 10 τ0 | NbN value of [27]; not measured for NbTiN |
| t_esc (phonon escape) | 530 τ0 (80 ps) | measured for 9 nm NbTiN [39]; [27] uses 15 ps for NbN |
| s_max | 1.25 (ε ≥ −0.25) | our assumption; never reached at f = 1 |
| \|A\| | 2.5 | calibration from the earlier conserving model, so f = 1 is the same photon |
| photon arrival t0 | 20 τ0 | |
| f | 1 (0.6, 0.35 for selectivity probes) | photon-energy proxy, not calibrated |

What [27] does not state, and we chose: the wall treatment, the amplitude, the lower bound on ε. The spread over
reasonable choices moves thresholds by 0.01–0.03 I_dep; the measured escape time (530 vs 100 τ0) moves them by up to
0.10. **The absolute photon scale is not calibrated.** Compare geometries, not absolute thresholds.

## 3 · Measuring detection

- **Dark switching current I_sw:** bisection without a photon in [0.70, 1.20] × I_dep.
- **Photon threshold I_det(x, y):** bisection at each site in [0.30, 1] × I_sw, to 0.34 % of I_dep(strip width).
- **Hard-step efficiency η(I):** area fraction of the unit cell with I_det < I, interpolated between sites and
  completed by the geometry's symmetry. i_10, i_50, i_95 are where η reaches 10, 50, 95 %; i_max is 100 %.
- **Go / no-go:** a hole device must beat its plain reference of the same width by ≥ 0.03 I_dep at the plateau (i_95).
  An onset (i_50) gain alone has never been robust: it depended on the criterion, window, mesh and hotspot.

**Two criteria, which are bounds, not alternatives:**
- **Voltage trigger (lower bound):** peak |V| > 0.45 V0 after t0 + 4 τ0, or a 25 τ0 trailing mean of |V| > 0.18 V0.
  It counts every resistive episode, including flux moving between holes and episodes that heal.
- **Sustained (upper bound):** the trigger, plus ≥ 1 phase slip in the last 50 τ0 of the window, plus ≥ 1 slip in the
  last 50 τ0 of a 150 τ0 photon-free continuation. Without heat, the TDGL must stay resistive on its own.
- **At γ = 10 the bounds coincide in plain strips and nearly did in the first phase's hole layouts** (median gap
  0.000 I_dep, up to 0.03 at a few slot sites). At γ = 1 they were 0.02 (200 nm) to 0.40 (30 nm strip) apart.
- **They do not coincide in the 40–60 nm channel layouts** of the second phase (docs/slots.md). A channel photon
  dumps the channel's current into the other channels with a few phase slips and the device heals with trapped
  flux. The gap is 0.15 I_dep in the double column.
- **With heat** (`snspd/thermal.py`, Zotova–Vodolazov heat equation [8]; 14-site check on 5 Nov) thresholds fell at or
  below the trigger bound at 11 of 14 sites. The trigger is the closer bound.

**Trapped flux** (`snspd/fluxoid.py`): the winding of arg ψ along a loop around each hole, counted only while
min |ψ| ≥ 0.05 on the loop. A trapped state can seed a new run.

## 4 · Protocol at γ = 10: the settled start

At γ = 10, |ψ| relaxes ~6× slower (u√(1+γ²|ψ|²)). After the bias step a strip needs ~200 τ0 to settle, so a photon
at t0 = 20 τ0 meets a condensate that is still too strong. The steady state is a stationary GL solution and does not
depend on γ, so each photon run starts from a state prepared by a γ = 1 dark run (`sim.settled_state`):
- 100 τ0 at the photon run's current, continued up to 3 × 100 τ0 until the last 25 τ0 has ≤ 0.02 phase slips and
  |V_end| ≤ 0.01; otherwise the run starts from the uniform state.
- Checked: under γ = 10 dynamics the settled state stays put (<|ψ|²> unchanged to 1e-4, |V| ≤ 0.011 over 200 τ0).
- **Pitfall (fixed 6 Nov):** judging stationarity on the second half of a single 100 τ0 run rejected 2 µm devices above
  ~0.92 I_dep. Their probe phase is still building up then (0.07 × 2π over 50–100 τ0), though they are stationary from
  ~120 τ0. The fallback start's bias-on transient then fires the trailing-mean trigger at t = 55 τ0. **A threshold run
  with t_trigger = 55.0 is a protocol artefact.**
- **Fallback to γ = 10 (7 Nov).** Between the γ = 1 and the γ = 10 dark I_sw (≈ 0.96–0.98 I_dep), the γ = 1 bias
  step overshoots into switching. The settling then retries with the device's own γ = 10 dynamics, again up to
  3 × 100 τ0. That only helps where a stationary state exists. The single slot column at 0.965 I_dep creeps
  (|V| ≈ 0.02) and switches after ~800 τ0, so it has none, and those runs still start uniform.
- Dark I_sw is bisected from the uniform start; at γ = 10 it comes out at ≈ 0.98–1.00 I_dep. That is partly a window
  effect: the single slot column switches at 0.965 after ~800 τ0. On 1.5–2 µm devices the bias-on transient alone
  can fire the trigger, so take the dark I_sw from the sustained criterion (`hitmap.py --isw-sustained`).

## 5 · Numerics and the checks behind them

| Choice | Evidence |
| --- | --- |
| Mesh: max edge 1.4 ξ = 7 nm (mean edge ≈ 4 nm) | 1.4 vs 0.8 ξ: \|Δi_det\| ≤ 0.007 for photon thresholds in strips. **Not converged for flux crossing narrow bridges:** a trapped-flux switching current moved 0.06 between meshes in 12 nm (2.4 ξ) bridges. Keep bridges ≥ ~8 ξ wide or check at 0.8 ξ |
| Uniform 2.0 ξ for screening only | γ = 1: shifts −0.006 to −0.013, differences kept within 0.007. γ = 10: per-site shifts ≤ 0.012, but one hole-vs-plain difference moved 0.02. 2.5 ξ: shifts up to −0.05, a key difference shrank 30 % |
| No graded meshes | run time scales ~h⁻⁴ (stable time step ∝ h_min²): 1.4 → 2.0 ξ is 3.7× faster, a graded mesh at most 2.4×. Coarsening near a hole rim dropped dark I_sw from 0.97 to 0.83: **hole rims, bridges and contacts must stay fine** |
| Pinned meshes (`SNSPD_MESH_DIR`, `python -m snspd.sim pin-meshes`) | holed geometries triangulate differently on other GEOS/shapely versions; 0.04 I_dep of triangulation noise was seen for flux-mediated quantities |
| **Solver time step: dt_max ≤ 0.025 (from 7 Nov 2026; 0.2 before)** | dt_max 0.2 leaves a local instability at single small cells: a site with |ψ| above its neighbours and μ offset by ~0.2, with a stationary |∇μ|² 50× the 99.9th percentile (dark, no photon). dt_max 0.05 removes it on 2.0 ξ meshes, 0.025 on the 1.4 ξ double column. Without heat it changes nothing measurable (dark I_sw and trigger thresholds within 0.003, `dt_check.py`). With heat it becomes a permanent hot spot that seeds the latch: heated thresholds moved by up to +0.07. Every heated result before 7 Nov 15:00 (`results/heat/check.csv`, `check_m1.4.csv`) is superseded by `check_dt0.025.csv`. Not the cause: the singular Neumann Laplacian of μ (pinning μ did not help) |
| Window 400 τ0 | γ = 10: 800 τ0 changes thresholds by ≤ 0.003 at 4 sites, although threshold runs fire late (330–395 τ0) |
| Shortened strips (second phase) | 150 nm of plain lead beyond the outermost slot reproduces a 2 µm strip within 0.011 I_dep at 1.4 ξ (double column, 3 sites + dark I_sw, `mesh_check.py`) |
| 2.0 ξ maps (second phase) | 1.4 − 2.0 ξ at the critical sites: mean +0.005 to +0.011, at most +0.019 (next to dark switching); channel sites ≤ 0.011 |
| Current crowding: disc average (1.5 ξ radius, ≥ 1.5 ξ from walls) | single-node maxima are mesh noise (a straight strip gave 1.2–1.4) |
| Strip length | 0.8 µm was short enough that flux exited at the far edge; ψ = 0 terminals nucleate dark switching 20–50 nm from the contacts; keep photon sites ≥ 0.27 µm from them |

## 6 · What the model leaves out

- **No thermal noise.** I_sw is an instability current, not a dark-count rate, so the model cannot show dark-count
  suppression. Barrier estimates need GL saddle points [21].
- **No circuit.** A current source keeps the bias; "latched" is not a readout pulse through 50 Ω.
- **Ideal strips.** Real 9 nm NbTiN strips switch at 0.63 I_dep (3 K) [40]; 0.43–0.74 for NbN/WSi devices [23].
  Every 200 nm f = 1 threshold so far was above 0.63, so the absolute photon scale is suspect. Every hole adds rims where
  fabrication defects and constrictions occur [31]; a hole close to the edge gave the lowest I_c in experiment.
- **Near-T_c GTDGL at 0.2–0.3 T_c** (see §1).

## 7 · Lessons for new geometries

- Every hole layout of the first phase failed the plateau test; holes add an insensitive tail at bridges and hole
  tips. Gains appeared only at onset and depended on definitions.
- Narrow bridges (≲ 4–5 ξ) are under-resolved and behave as 1D links; vortices may not form below ~25 nm [Incalza 2025].
- Trapped flux below threshold makes a holed device history-dependent.
- Always run the plain strip of the same width and length as the reference, with the same preset and mesh.

## References

1. L. Bishop-Van Horn, pyTDGL, *Comput. Phys. Commun.* **291**, 108799 (2023), arXiv:2302.03812.
2. L. Kramer and R. J. Watts-Tobin, *Phys. Rev. Lett.* **40**, 1041 (1978).
3. R. J. Watts-Tobin, Y. Krähenbühl and L. Kramer, *J. Low Temp. Phys.* **42**, 459 (1981).
8. A. N. Zotova and D. Y. Vodolazov, *Phys. Rev. B* **85**, 024509 (2012), arXiv:1112.3790.
13. M. Jönsson *et al.*, *Phys. Rev. Applied* **17**, 064046 (2022), arXiv:2112.05443.
21. D. Y. Vodolazov, *Phys. Rev. B* **85**, 174507 (2012), arXiv:1204.0365.
23. S. Frasca *et al.*, *Phys. Rev. B* **100**, 054520 (2019), arXiv:1904.08596.
27. Simon *et al.*, *Phys. Rev. B* **112**, 174512 (2025), arXiv:2501.13791.
31. Zhang *et al.*, constriction-limited dark counts in NbN nanowires (2014), arXiv:1406.2043.
39. M. Sidorova *et al.*, *Phys. Rev. B* **104**, 184514 (2021), arXiv:2108.04717.
40. M. Sidorova *et al.*, *Supercond. Sci. Technol.* (2022), DOI 10.1088/1361-6668/ac8454, arXiv:2204.01423.

Numbering follows the old report, so citations in old notes still match.
