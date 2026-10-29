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
  fluxoid.py            fluxoid (trapped flux quanta) of every hole from the phase winding of psi
  analysis.py           threshold maps, efficiency curves, quantiles, selectivity
  store.py              resumable CSV / JSON results
  plotting.py           figure style and standard figures
hitmap.py               photon-threshold maps: hole geometries vs the plain 200 nm wire, go / no-go
narrow_baseline.py      plain strips of any width (control) and the width scan
checks.py               tests of the map assumptions: symmetry, edge band, trigger delay
fluxtrap.py             flux a photon leaves trapped in the holes: state scan, trapped-state I_sw, traces
results/hitmap/         plain wire vs long zipper            (hitmap.py)        } voltage-only criterion,
results/narrow_baseline/ plain strips                         (narrow_baseline.py)} from 29 Oct 2026
results/checks/         symmetry, edge-band and delay checks (checks.py)          }
results/fluxtrap/       flux-trapping scan, trapped states, traces (fluxtrap.py)  }
results/v1_voltage_or_collapse/  all results with the first criterion (23-29 Oct), see its README
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
.venv/bin/python fluxtrap.py scan    --workers 9                 # every zipper site x 0.40-0.96 I_dep, 400 tau0
.venv/bin/python fluxtrap.py scan    --geom Plain_Wire --biases 0.80:0.92:0.01   # control
.venv/bin/python fluxtrap.py trapped --workers 6                 # I_sw of the trapped states, second photon
.venv/bin/python fluxtrap.py trapped --mesh-xi 0.8 --from=-68,0@0.77 --no-second   # one state, finer mesh
.venv/bin/python fluxtrap.py trace   --site=-68,0 --biases 0.70,0.77,0.80,0.85 --workers 4
```

At 1.4 xi a 400 tau0 run takes ~16 s on one core. The laptop CPU (Core Ultra 7 155H) has 6 performance cores and 10
efficiency cores; more than ~9 workers, or other jobs alongside, slow every run down sharply.

Result files are tagged `<hotspot>_T<tmax>_m<mesh_xi>`, e.g. `summary_conserving_T150_m0.8.txt`.

## Model and units

* **Solver:** generalized TDGL (Kramer & Watts-Tobin, PRL 40, 1041 (1978)) in pyTDGL (Bishop-Van Horn, Comput.
  Phys. Commun. 291, 108799 (2023)); no screening, no applied field, γ = 1, u = 5.79.
* **Film** (`snspd/units.py`): ξ = 5 nm, λ = 600 nm, d = 9 nm; strip 200 nm × 0.8 µm with film-width terminals.
* **Units:** lengths in ξ, sheet current in K0 = 4ξB_c2/(μ0Λ) = 5.238 µA/nm, time in τ0 = μ0σλ². The GL depairing
  current is j_dep = K0/(6√3) = 0.504 µA/nm, so I_dep = j_dep × narrowest free width (100.8 µA for 200 nm). These are
  computed from pyTDGL's own definitions (`Film.K0_uA_per_nm` equals `device.K0`).
* **Detection** (`snspd/sim.py`): voltage only. A run counts as a detection if the probe voltage (probes at
  y = ±L/3) exceeds 0.45 V0 at any time after t0 + 4 τ0, or its 25 τ0 trailing mean exceeds 0.18 V0. Each threshold
  records the branch and the trigger time (a latency figure of merit). Until 29 Oct a collapse of the supercurrent
  through one monitored cross-section also counted; it fired on local channel blocking with no voltage and was removed.
* **Bisection tolerance** is a fraction of I_dep (0.34 %), not fixed µA, so it stays correct if λ changes.
* **Hotspot** (`snspd/hotspot.py`): r0 = 14 nm, A = −2.5, D = 10 nm²/τ0, t_relax = 35 τ0, photon at 20 τ0, window
  150 τ0; energy-conserving diffusion with reflecting walls.

**λ and ρ only relabel the units.** With screening off, λ enters only K0 and τ0, so every normalised result
(i_det = I_det/I_sw, I_sw/I_dep) is independent of λ, absolute currents scale as 1/λ² and times as λ²/ρ_n. No result
needs ρ_n; it only converts τ0 to picoseconds (pass `--rho-uohm-cm` to print times in ps). What the hotspot parameters
mean physically depends on both:

| λ | ρ_n [µΩ cm] | τ0 | D = 10 nm²/τ0 | t_relax = 35 τ0 | I_dep(200 nm) |
|---|---|---|---|---|---|
| 175 nm (until 28 Oct) | 100–250 | 0.015–0.038 ps | 2.6–6.5 cm²/s | 0.5–1.4 ps | 1185 µA |
| **600 nm (current)** | 100–250 | 0.18–0.45 ps | 0.22–0.55 cm²/s | 6–16 ps | 101 µA |

λ = 600 nm follows the fit to 9 nm NbTiN hairpins (Jönsson et al., PR Applied 17, 064046 (2022)). With it the hotspot
parameters match literature values for NbN-like films, D ≈ 0.5 cm²/s and τ_e-ph ≈ 17 ps (Zotova & Vodolazov, PRB 85,
024509 (2012)), and I_sw(200 nm) ≈ 100 µA is in the range of real NbTiN strips.

**Time unit without a measured resistivity.** ρ_n depends on fabrication, but in a dirty superconductor it is tied to
λ: λ² ≈ ħρ_n/(πμ0Δ0) (dirty-limit BCS, Δ0 = 1.764 k_B T_c). Then τ0 = μ0λ²/ρ_n ≈ ħ/(πΔ0), which depends on T_c only.
For T_c = 8–12 K this gives τ0 ≈ 0.17–0.12 ps (and implies ρ_n ≈ 260–390 µΩ cm for λ = 600 nm), so 150 τ0 ≈ 17–26 ps,
400 τ0 ≈ 46–69 ps, D = 10 nm²/τ0 ≈ 0.6–0.9 cm²/s and t_relax = 35 τ0 ≈ 4–6 ps. This uses zero-temperature BCS
quantities in a model valid near T_c, so treat it as an order-of-magnitude scale for choosing the detection window.

**Change of λ on 28 Oct 2026.** All results were simulated with λ = 175 nm. Their stored currents (CSV, JSON and
crowding maps) were multiplied once by (175/600)² = 0.08507; every normalised value is unchanged. Fresh λ = 600 nm
runs reproduce the rescaled values within the bisection tolerance (plain I_sw 99.50 vs 99.41 µA, 56 nm I_sw 27.92 vs
27.93 µA, plain centre i_det 0.852 vs 0.852). Absolute currents in older documents (e.g. the 28 Oct meeting report)
are in λ = 175 nm units: multiply by 0.0851 to compare.

## Current results

### Voltage-only criterion, 400 τ0 window, 1.4 ξ mesh (29 Oct 2026)

| Result | File | Key numbers (units of I_dep of the narrowest cross-section) |
|---|---|---|
| Plain wire vs long zipper, full grid | `results/hitmap/summary_conserving_T400_m1.4.txt` | plain i_50 0.845, i_95 0.850; zipper i_50 0.795, i_95 0.853, i_max 0.907; go / no-go 1/3 (C3 passes by 0.001, with one non-detecting probe per geometry) |
| Width scan + zipper | `results/narrow_baseline/width_scan_summary_conserving_T400_m1.4.txt` | i_50 = 0.327 / 0.603 / 0.766 / 0.845 at 30 / 56 / 100 / 200 nm; zipper 0.795 |
| Narrow-strip comparison | `results/narrow_baseline/narrow56nm_summary_conserving_T400_m1.4.txt` | 56 nm: i_50 0.603, i_95 0.614: beats the zipper on onset and plateau |
| Flux trapping: zipper scan (42 sites × 0.40–0.96 I_dep), plain-wire control, trapped states | `results/fluxtrap/summary_{Long_Zipper,Plain_Wire}_conserving_T400_m1.4.txt` | zipper i_50 0.800 / 0.806 / 0.821 counting any trigger / sustained / prompt switches (plain wire 0.850 for all); sub-threshold photons leave 1–2 Φ0 from 0.64 I_dep; trapped flux lowers I_sw from 0.968 to 0.87–0.91 |

All strips and the zipper switch at I_sw/I_dep = 0.968–0.971 with this window. The 30 nm strip agrees with a 0.8 ξ
run within 0.006. Detection latency at threshold (t_trigger − t0): plain wire median 93 τ0, max 129 τ0; zipper median
175 τ0, max 366 τ0, with 18 of 42 sites above 200 τ0.

### Flux trapping (29 Oct 2026, `fluxtrap.py`, `results/fluxtrap/`)

Every scan run records two things:
- the fluxoid of each hole, read from the phase winding of ψ (`snspd/fluxoid.py`) at least every 8 τ0;
- the phase slipped between the probes, |∫V dt|/2π. One vortex crossing the full width counts as 1.

Each run gets one of five outcomes:
- *clean*;
- *trapped*: no trigger, flux left in a hole;
- *transient*: the hitmap criterion fires, but the last 50 τ0 have < 1 phase slip, so the device is quiet again;
- *direct*: the device stays resistive;
- *via trap*: the device stays resistive, after ≥ 20 τ0 at zero voltage with flux trapped.

A switch is *prompt* if the device stays resistive and has slipped 2 phases by t = 150 τ0, i.e. 130 τ0 after the photon.

* **Sub-threshold photons leave flux behind.**
  - Near the edges and in the left channel this starts at 0.64 I_dep, ~0.14 I_dep below the detection threshold. The
    trapped area is 59 % of the unit cell at 0.70 I_dep and 85 % at 0.75.
  - There is no voltage signal: a quantum entering a hole costs ~0.4 phase slip.
  - The flux arrives 18–102 τ0 after the photon.
  - The 319 trapped runs end in four states (holes numbered as in the summary), never more than two quanta:

    | State | Runs | Where the photon lands |
    |---|---|---|
    | one quantum in the central left hole, [0 −1 0 0 0 0 0] | 180 | near the left edge |
    | a vortex–antivortex pair split between the left and right columns, [0 1 0 0 0 −1 0] | 103 | centre channel (the current tears the pair apart) |
    | one quantum in each of two left holes, [0 −1 −1 0 0 0 0] | 23 | at the top bridge |
    | two quanta in the central left hole, [0 −2 0 0 0 0 0] | 10 | left channel, just below threshold |

* **Trapped flux stays and lowers the switching current.**
  - Each state was used as the initial condition, without a photon. It keeps its flux for 400 τ0 at zero bias, and
    also just below its own switching current.
  - That switching current is 0.909 I_dep for one quantum and 0.871–0.874 for two, against 0.968 for the empty device.
    Producing the same state from different photons gives the same value.
  - The exception is [0 −2]. At 0.775 I_dep its flux leaves the left hole (2 phase slips, no lasting voltage), and
    two antivortices enter a right-column hole from the right edge, giving [0 0 0 0 0 2 0]. So it depins but does not
    switch.
* **The zipper's onset is flux-mediated.**
  - At 0.78 I_dep, photons in the left channel and on the centre line (24 % of the area) fire the hitmap criterion
    through a rearrangement like this one. The device is left quiet with new trapped flux: the hitmap counts these as detections, but they are
    not sustained.
  - From 0.79 the same sequence ends in a resistive state, after a zero-voltage dwell with the flux trapped. The
    resistive onset comes 180–210 τ0 after the photon at 0.79, and 80–250 τ0 over all via-trap runs.
  - At 0.80 I_dep, 45 % of the area switches and stays resistive, all of it via trapped flux.
  - Prompt switching reaches 48 % of the area only at 0.82.
  - Direct switching starts at 0.81–0.85 in the channels.
  - The plain wire never traps flux, and every detection is prompt and sustained. Counting only prompt detections,
    the zipper's onset gain shrinks from 0.05 to 0.03 I_dep, and its plateau (i_95) comes 0.008 I_dep after the plain
    wire's.

  Efficiency-curve points in units of I_dep. The scan grid is 0.01, so values can be up to 0.01 above the bisected
  hitmap thresholds.

  | | i_10 | i_50 | i_95 | i_max |
  |---|---|---|---|---|
  | zipper, any trigger (hitmap criterion) | 0.780 | 0.800 | 0.855 | 0.910 |
  | zipper, sustained | 0.790 | 0.806 | 0.856 | 0.910 |
  | zipper, prompt | 0.810 | 0.821 | 0.858 | 0.910 |
  | zipper, direct (no trapped dwell) | 0.825 | 0.834 | 0.859 | 0.910 |
  | plain wire, all four | 0.840 | 0.850 | 0.850 | 0.850 |

* **A trapped state changes the response to the next photon.** I measured the photon threshold at four sites,
  starting from a trapped state instead of the empty device:
  - [0 −2] lowers it by 0.09–0.15 I_dep away from the trap, down to 0.63 at (68, 61).
  - One or two quanta shift it by −0.07 to +0.09. The sign depends on the side of the hole: one quantum in the
    central left hole raises the left-channel threshold by 0.045 and lowers the others by 0.01–0.04.
  - So the zipper's response depends on what it absorbed before.
* **Mesh.** At 0.8 ξ, the (−68, 0) photon at 0.77 I_dep leaves the same [0 −2] state.
  - Its depinning current is 0.795 I_dep, against 0.775 at 1.4 ξ. The rearrangement is the same, mirrored to the
    lower right-column hole.
  - The 12 nm bridges are only 2.4 ξ wide, so the trapping thresholds depend more on the mesh than the photon
    thresholds do (≤ 0.007).
  - The zipper's long-window onset may therefore sit ~0.02 I_dep higher at 0.8 ξ.
* **Limits.**
  - The TDGL here is deterministic, without thermal noise, so a trapped vortex never escapes by activation.
  - Rough estimate for a real device: the pinning barrier of a hole is of order ε0 d = Φ0²/(2πμ0Λ) ≈ 0.04 eV
    ≈ 500 K × k_B (Λ = 80 µm). So at zero bias the flux should also stay in a real device at 1–2 K. Near the
    depinning current the barrier vanishes.
  - The hotspot's heat and the readout circuit are not modelled. "Sustained" means the TDGL state keeps slipping
    phase over the last 50 τ0 of 400; it does not mean a latched hotspot.

### First criterion (voltage or channel collapse), 150 τ0 window

These live in `results/v1_voltage_or_collapse/`. A voltage-only re-run of the 1.4 ξ critical grid at 150 τ0
reproduces them within 0.003.

| Result | File (under `results/v1_voltage_or_collapse/`) | Key numbers |
|---|---|---|
| Plain wire vs long zipper (reference) | `hitmap/summary_conserving_T150_m0.8.txt` | i_50 0.857 vs 0.831, i_95 0.861 vs 0.872, i_max 0.862 vs 0.925; go / no-go 0/3 |
| Window check, critical grid, 1.4 ξ | `hitmap/summary_conserving_T{150,250,400}_m1.4.txt` | in I_dep units: plain stable (i_50 0.843–0.845); zipper i_50 0.817 → 0.803 → 0.800 |
| Same at 1.4 ξ, critical grid; legacy hotspot | `hitmap/summary_{conserving,legacy}_T150_m1.4.txt` | mesh convergence ≤ 0.007 at shared sites |
| Narrow-strip control | `narrow_baseline/narrow56nm_summary_conserving_T150_m0.8.txt` | 56 nm: i_50 0.607, i_95 0.618 (beats the zipper by 0.22–0.25) |
| Width scan | `narrow_baseline/width_scan_summary_conserving_T150_m0.8.txt` | i_50 = 0.33 / 0.61 / 0.78 / 0.86 at 30 / 56 / 100 / 200 nm |
| Edge banks (parked) | `parked/edge_bank/results/edge_bank*_summary_*.txt` | no plateau gain in 16 configurations |

Every geometry switches at the depairing current of its narrowest cross-section (I_sw/I_dep = 0.97–1.03).

### Checks of the method (29 Oct 2026, `checks.py`, `results/v1_voltage_or_collapse/checks/`)

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

Write a builder in `snspd/geometry.py` that returns a `Geometry` (holes as ellipses, hit grids, probe sites,
symmetry), register it in `GEOMETRIES`, and run `hitmap.py run --geoms Plain_Wire,<Key>`. If the
layout has a symmetry other than `strip` or `zipper`, add its map rule to `snspd/analysis.unit_map`.

## Next steps for the hole study

0. **Flux trapping, follow-up.** The first map is done (see above).
   - Decide which detection counts for the thesis: any trigger, sustained, or prompt. Then make `sim.detect`
     require it (e.g. ≥ 1 phase slip in the last 50 τ0), so the hitmap no longer counts flux rearrangements.
   - Check the mesh on the trapping thresholds themselves: a 0.8 ξ scan of a few left-channel and centre sites at
     0.62–0.84 I_dep, and 0.5 ξ at the 12 nm bridges.
   - Photon sequences at sub-threshold bias: does accumulated flux reach a state that switches? Does a switching
     event clear the trapped flux?
   - A longer strip with more hole periods (end effects).
1. Hole-length sweep of the zipper (channel length / width from 2 to ~6): separates "current bypasses the short
   channel" from "detection needs the whole width to switch" as the reason the 56 nm channels do not act as narrow
   strips.
2. Re-run the column geometries of the archived cascade scripts (2-channel, aligned 3-channel, short zipper) with the
   current protocol (conserving hotspot, 0.8 ξ, full grid).
3. Hotspot with an energy sink (loss time 5–20 τ0 and ∞) on the zipper and the 56 nm strip.
4. Detection criterion: `fluxtrap.py` now separates sustained from transient switching (item 0). Still to do: repeat
   the key thresholds at 250 τ0.
