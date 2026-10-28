# Nanohole SNSPD simulations

TDGL simulations (pyTDGL 0.9) of superconducting-nanowire single-photon detector strips with **internal holes**: does
perforating the strip lower the bias needed to detect photons (the plateau) relative to the dark switching current?

## Layout

```
snspd/                  shared library, imported by every current script
  units.py              film parameters and the unit system (K0, j_dep, tau0, Pearl length)
  geometry.py           strips and hole layouts, photon-hit grids, unit cells  <- add new hole layouts here
  hotspot.py            photon hotspot models (energy-conserving FEM diffusion; legacy free-space Gaussian)
  sim.py                device, detection criterion, bisection, current crowding, process-pool runner
  analysis.py           threshold maps, efficiency curves, quantiles, selectivity
  store.py              resumable CSV / JSON results
  plotting.py           figure style and standard figures
hitmap.py               photon-threshold maps: hole geometries vs the plain 200 nm wire, go / no-go
narrow_baseline.py      plain strips of any width (control) and the width scan
checks.py               tests of the map assumptions: symmetry, edge band, trigger delay vs bias
results/hitmap/         plain wire vs long zipper            (hitmap.py)
results/narrow_baseline/ plain strips 30-200 nm               (narrow_baseline.py)
results/checks/         symmetry, edge-band and delay checks (checks.py)
parked/edge_bank/       edge-bank strip: code + results, parked until the hole study is done
archive/                unreliable and superseded early work, see archive/README.md
docs/                   literature review
tests/                  early single/two-hole tests (not maintained)
```

## Running

Use the project environment (`.venv/bin/python`; the system Python has no scipy/tdgl). Every run is resumable:
finished thresholds are appended to a CSV and skipped next time, so an interrupted run just continues.

```bash
.venv/bin/python hitmap.py plan   --mesh-xi 0.8                  # what is still missing
.venv/bin/python hitmap.py run    --mesh-xi 0.8 --workers 6       # ~1.8 GB per worker at 0.8 xi; 15 GB WSL -> <= 6
.venv/bin/python hitmap.py report --mesh-xi 0.8                  # rebuild figures + summary from the CSV
.venv/bin/python narrow_baseline.py run  --width-nm 56 --mesh-xi 0.8 --compare hitmap_conserving_T150_m0.8.csv
.venv/bin/python narrow_baseline.py scan --widths 30,56,100 --mesh-xi 0.8 --compare hitmap_conserving_T150_m0.8.csv
.venv/bin/python hitmap.py selftest                              # hotspot-solver checks
```

Result files are tagged `<hotspot>_T<tmax>_m<mesh_xi>`, e.g. `summary_conserving_T150_m0.8.txt`.

## Model and units

* **Solver:** generalized TDGL (Kramer & Watts-Tobin, PRL 40, 1041 (1978)) in pyTDGL (Bishop-Van Horn, Comput.
  Phys. Commun. 291, 108799 (2023)); no screening, no applied field, γ = 1, u = 5.79.
* **Film** (`snspd/units.py`): ξ = 5 nm, λ = 600 nm, d = 9 nm; strip 200 nm × 0.8 µm with film-width terminals.
* **Units:** lengths in ξ, sheet current in K0 = 4ξB_c2/(μ0Λ) = 5.238 µA/nm, time in τ0 = μ0σλ². The GL depairing
  current is j_dep = K0/(6√3) = 0.504 µA/nm, so I_dep = j_dep × narrowest free width (100.8 µA for 200 nm). These are
  computed from pyTDGL's own definitions (`Film.K0_uA_per_nm` equals `device.K0`).
* **Bisection tolerance and the collapse threshold** are fractions of I_dep (0.34 % and 4.2 %), not fixed µA, so they
  stay correct if λ changes.
* **Hotspot** (`snspd/hotspot.py`): r0 = 14 nm, A = −2.5, D = 10 nm²/τ0, t_relax = 35 τ0, photon at 20 τ0, window
  150 τ0; energy-conserving diffusion with reflecting walls.

**λ and ρ only relabel the units.** With screening off, λ enters only K0 and τ0, so every normalised result
(i_det = I_det/I_sw, I_sw/I_dep) is independent of λ, absolute currents scale as 1/λ² and times as λ²/ρ_n. ρ_n has
never been set (pass `--rho-uohm-cm` to print times in ps). What the hotspot parameters mean physically depends on
both:

| λ | ρ_n [µΩ cm] | τ0 | D = 10 nm²/τ0 | t_relax = 35 τ0 | I_dep(200 nm) |
|---|---|---|---|---|---|
| 175 nm (until 28 Oct) | 100–250 | 0.015–0.038 ps | 2.6–6.5 cm²/s | 0.5–1.4 ps | 1185 µA |
| **600 nm (current)** | 100–250 | 0.18–0.45 ps | 0.22–0.55 cm²/s | 6–16 ps | 101 µA |

λ = 600 nm follows the fit to 9 nm NbTiN hairpins (Jönsson et al., PR Applied 17, 064046 (2022)). With it the hotspot
parameters match literature values for NbN-like films, D ≈ 0.5 cm²/s and τ_e-ph ≈ 17 ps (Zotova & Vodolazov, PRB 85,
024509 (2012)), and I_sw(200 nm) ≈ 100 µA is in the range of real NbTiN strips. ρ_n is still open: set it to the
film's measured value to pin τ0.

**Change of λ on 28 Oct 2026.** All results were simulated with λ = 175 nm. Their stored currents (CSV, JSON and
crowding maps) were multiplied once by (175/600)² = 0.08507; every normalised value is unchanged. Fresh λ = 600 nm
runs reproduce the rescaled values within the bisection tolerance (plain I_sw 99.50 vs 99.41 µA, 56 nm I_sw 27.92 vs
27.93 µA, plain centre i_det 0.852 vs 0.852). Absolute currents in older documents (e.g. the 28 Oct meeting report)
are in λ = 175 nm units: multiply by 0.0851 to compare.

## Current results

| Result | File | Key numbers |
|---|---|---|
| Plain wire vs long zipper (reference) | `results/hitmap/summary_conserving_T150_m0.8.txt` | i_50 0.857 vs 0.831, i_95 0.861 vs 0.872, i_max 0.862 vs 0.925; go / no-go 0/3 |
| Window check, critical grid, 1.4 ξ | `results/hitmap/summary_conserving_T{150,250,400}_m1.4.txt` | in I_dep units: plain stable (i_50 0.843–0.845); zipper i_50 0.817 → 0.803 → 0.800 |
| Same at 1.4 ξ, critical grid; legacy hotspot | `results/hitmap/summary_{conserving,legacy}_T150_m1.4.txt` | mesh convergence ≤ 0.007 at shared sites |
| Narrow-strip control | `results/narrow_baseline/narrow56nm_summary_conserving_T150_m0.8.txt` | 56 nm: i_50 0.607, i_95 0.618 (beats the zipper by 0.22–0.25) |
| Width scan | `results/narrow_baseline/width_scan_summary_conserving_T150_m0.8.txt` | i_50 = 0.33 / 0.61 / 0.78 / 0.86 at 30 / 56 / 100 / 200 nm |
| Edge banks (parked) | `parked/edge_bank/results/edge_bank*_summary_*.txt` | no plateau gain in 16 configurations |

Every geometry switches at the depairing current of its narrowest cross-section (I_sw/I_dep = 0.97–1.03).

### Checks of the method (29 Oct 2026, `checks.py`, `results/checks/`)

* **Mesh:** thresholds at 0.8 ξ and 1.4 ξ agree within 0.007.
* **Edge band:** x = −98 nm agrees with x = −96 nm within 0.003, so filling the outer 4 nm is fine.
* **Symmetry:** plain-wire mirror images (including hits on the monitored slice), the zipper's y-mirror and its glide
  images in the centre channel and at the bridges all agree within the bisection tolerance. **Exception:** a photon
  inside the zipper's monitored channel, (68, 61), "detects" at 0.724 vs 0.801 at its glide image, by supercurrent
  collapse with no voltage: the per-channel collapse branch fires on local channel blocking. The maps are not affected
  (only x ≤ 0 is simulated), but the criterion is not a device-level one for channel geometries.
* **Window:** the plain wire's photon thresholds are window-independent (delay ≤ 95 τ0 just above threshold), but
  several zipper channel and edge sites switch slowly (150–370 τ0 after the photon) and their thresholds keep falling
  with the window (site (−96, 0): 0.803 → 0.787 → 0.777 I_dep at 150 / 250 / 400 τ0). The dark I_sw of the plain wire
  also drifts (0.986 → 0.970 I_dep), so i_det = I_det/I_sw mixes in a window effect; quote thresholds in I_dep units.
  The conclusions survive: zipper onset earlier by 0.03–0.045 I_dep, plateau later by 0.025–0.03 I_dep, and the 56 nm
  strip (i_50 0.60, i_95 0.61 I_dep at 400 τ0) beats it on both.
* **Current crowding:** the old K_peak/K_uniform (1.21 plain, 1.45 zipper at 1.4 ξ) was dominated by single-node and
  boundary-node noise; a straight strip must give 1.00. The crowding peak is now the disc average over centres at
  least 1.5 ξ from all walls: 1.04 / 1.08 (plain, 0.8 / 1.4 ξ) and 1.24 / 1.24 (zipper). At I_sw both reach 1.04–1.06
  j_dep locally, i.e. both switch at the depairing current.

## Adding a hole geometry

Write a builder in `snspd/geometry.py` that returns a `Geometry` (holes as ellipses, monitored cross-section, hit
grids, probe sites, symmetry), register it in `GEOMETRIES`, and run `hitmap.py run --geoms Plain_Wire,<Key>`. If the
layout has a symmetry other than `strip` or `zipper`, add its map rule to `snspd/analysis.unit_map`.

## Next steps for the hole study

0. **Fix the detection definition before new production runs.** (a) Replace the per-channel collapse branch by the
   total supercurrent through a full-width cross-section, or drop it and rely on the probe voltage. (b) Choose the
   window as a physical detection time once ρ_n is known (150 / 400 τ0 = 27–68 / 72–180 ps for ρ_n = 100–250 µΩ cm)
   and keep the recorded trigger time as a latency figure of merit: the zipper's channel sites are slow.
   (c) Re-run the plain wire, long zipper and 56 nm strip under that definition (1.4 ξ full grid, 0.8 ξ spot checks).
   Also worth a check: 0.5 ξ mesh at the 12 nm bridges, and a longer strip with more hole periods (end effects).
1. Hole-length sweep of the zipper (channel length / width from 2 to ~6): separates "current bypasses the short
   channel" from "detection needs the whole width to switch" as the reason the 56 nm channels do not act as narrow
   strips.
2. Re-run the column geometries of the archived cascade scripts (2-channel, aligned 3-channel, short zipper) with the
   current protocol (conserving hotspot, 0.8 ξ, full grid).
3. Hotspot with an energy sink (loss time 5–20 τ0 and ∞) on the zipper and the 56 nm strip.
4. Detection criterion: record "vortex crossed" and "sustained collapse" separately; repeat key thresholds at 250 τ0.
