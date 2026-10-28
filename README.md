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
results/hitmap/         plain wire vs long zipper            (hitmap.py)
results/narrow_baseline/ plain strips 30-200 nm               (narrow_baseline.py)
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
* **Film** (`snspd/units.py`): ξ = 5 nm, λ = 175 nm, d = 9 nm; strip 200 nm × 0.8 µm with film-width terminals.
* **Units:** lengths in ξ, sheet current in K0 = 4ξB_c2/(μ0Λ) = 61.57 µA/nm, time in τ0 = μ0σλ². The GL depairing
  current is j_dep = K0/(6√3) = 5.925 µA/nm, so I_dep = j_dep × narrowest free width. These are computed from pyTDGL's
  own definitions (`Film.K0_uA_per_nm` equals `device.K0`).
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
| 175 nm (current) | 100–250 | 0.015–0.038 ps | 2.6–6.5 cm²/s | 0.5–1.4 ps | 1185 µA |
| 600 nm | 100–250 | 0.18–0.45 ps | 0.22–0.55 cm²/s | 6–16 ps | 101 µA |

Literature values for NbN-like films are D ≈ 0.5 cm²/s and τ_e-ph ≈ 17 ps (Zotova & Vodolazov, PRB 85, 024509
(2012)), and a fit to 9 nm NbTiN gives λ ≈ 600 nm (Jönsson et al., PR Applied 17, 064046 (2022)). So the simulations
already run are physically consistent if the film has λ ≈ 600 nm. Choosing λ and ρ_n for the target film is open; it
needs no new simulations, only a rescaling of the stored µA values.

## Current results

| Result | File | Key numbers |
|---|---|---|
| Plain wire vs long zipper (reference) | `results/hitmap/summary_conserving_T150_m0.8.txt` | i_50 0.857 vs 0.831, i_95 0.861 vs 0.872, i_max 0.862 vs 0.925; go / no-go 0/3 |
| Same at 1.4 ξ, critical grid; legacy hotspot | `results/hitmap/summary_{conserving,legacy}_T150_m1.4.txt` | mesh convergence ≤ 0.007 at shared sites |
| Narrow-strip control | `results/narrow_baseline/narrow56nm_summary_conserving_T150_m0.8.txt` | 56 nm: i_50 0.607, i_95 0.618 (beats the zipper by 0.22–0.25) |
| Width scan | `results/narrow_baseline/width_scan_summary_conserving_T150_m0.8.txt` | i_50 = 0.33 / 0.61 / 0.78 / 0.86 at 30 / 56 / 100 / 200 nm |
| Edge banks (parked) | `parked/edge_bank/results/edge_bank*_summary_*.txt` | no plateau gain in 16 configurations |

Every geometry switches at the depairing current of its narrowest cross-section (I_sw/I_dep = 0.99–1.03).

## Adding a hole geometry

Write a builder in `snspd/geometry.py` that returns a `Geometry` (holes as ellipses, monitored cross-section, hit
grids, probe sites, symmetry), register it in `GEOMETRIES`, and run `hitmap.py run --geoms Plain_Wire,<Key>`. If the
layout has a symmetry other than `strip` or `zipper`, add its map rule to `snspd/analysis.unit_map`.

## Next steps for the hole study

1. Hole-length sweep of the zipper (channel length / width from 2 to ~6): separates "current bypasses the short
   channel" from "detection needs the whole width to switch" as the reason the 56 nm channels do not act as narrow
   strips.
2. Re-run the column geometries of the archived cascade scripts (2-channel, aligned 3-channel, short zipper) with the
   current protocol (conserving hotspot, 0.8 ξ, full grid).
3. Hotspot with an energy sink (loss time 5–20 τ0 and ∞) on the zipper and the 56 nm strip.
4. Detection criterion: record "vortex crossed" and "sustained collapse" separately; repeat key thresholds at 250 τ0.
