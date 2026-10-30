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
results/hitmap/         plain wire vs long zipper            (hitmap.py)        } sustained-switching criterion,
results/narrow_baseline/ plain strips                         (narrow_baseline.py)} from 29 Oct 2026
results/checks/         symmetry, edge-band and delay checks (checks.py)          }
results/fluxtrap/       flux-trapping scan, trapped states, traces (fluxtrap.py)  }
results/v2_voltage_trigger/      results with the voltage trigger alone (28-29 Oct), see its README
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
.venv/bin/python fluxtrap.py scan    --workers 9                 # every zipper site x 0.56-0.96 I_dep, 400 tau0
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
* **Detection** (`snspd/sim.py`): sustained switching. A run counts as a detection if two conditions hold:
  - the probe voltage triggers: it exceeds 0.45 V0 at any time after t0 + 4 τ0, or its 25 τ0 trailing mean exceeds
    0.18 V0 (probes at y = ±L/3);
  - the device is still resistive at the end of the window: at least 1 phase slip between the probes (|∫V dt|/2π) in
    the last 50 τ0;
  - it is still resistive after another 150 τ0. The run is continued from its final state without the photon (the
    hotspot has relaxed by then), and the same test is applied to the last 50 τ0.

  The continuation is needed because near threshold a flux rearrangement between holes (1–6 slips over up to
  ~100 τ0) can be delayed into the last 50 τ0 of the window. Without the continuation, six zipper sites "detected" at
  0.777 I_dep on such a rearrangement.

  Each threshold records the voltage branch and the trigger time (a latency figure of merit).

  History:
  - Until 29 Oct a collapse of the supercurrent through one monitored cross-section also counted. It fired on local
    channel blocking with no voltage, so it was removed (`results/v1_voltage_or_collapse/`).
  - Until 29 Oct the voltage trigger alone counted. In the zipper it also fires on trapped flux moving between holes
    with no lasting voltage (`results/v2_voltage_trigger/`).
  Reported as bounds: the voltage trigger is a lower bound on the detection threshold and sustained switching an
  upper bound, because the model has no Joule heating to latch the resistive state (see Current results).
* **Temperature** is not a parameter. It enters only through ξ(T) and λ(T), which are fixed at 5 nm and 600 nm,
  and through the GL equation itself, which is an expansion near T_c. The values fit 9 nm NbTiN at its operating
  temperature (~0.2–0.3 T_c). GL is therefore used outside its formal range: currents in units of I_dep (GL
  depairing current) are the more trustworthy output, absolute values less so.
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

### Two criteria as bounds, 400 τ0 window, 1.4 ξ mesh (29 Oct 2026)

Neither criterion is the real detector:
- The **voltage trigger** (`results/v2_voltage_trigger/`) counts any resistive episode, including flux moving
  between holes and episodes that heal on their own. It is the **lower bound** on the detection threshold.
- **Sustained switching** (current `sim.detect`, `results/hitmap/`, `results/narrow_baseline/`) asks the TDGL state
  to stay resistive by itself. The model has no Joule heating, so a strip that is resistive for 100–200 τ0 can heal
  once the hotspot has relaxed. It is the **upper bound**.

A real strip lies between the two. Joule heating keeps a strip resistive at lower currents than this model does.
Whether a short episode is counted also depends on the readout.

The narrow strips show this most clearly. The 30 nm strip at 0.75 I_dep fires the trigger 4 τ0 after the photon and
slips 65 phases, about 10 eV at that current (12 photons' worth at 1550 nm), before it heals at ~200 τ0. At 0.40 it slips 2 phases and heals. Under
the sustained criterion every strip from 30 to 100 nm switches at 0.82 ± 0.01 I_dep, whatever photon it absorbs
(selectivity 0.008 at 30 nm). That floor is where this model can stay resistive by itself. It does not measure photon
detection.

Efficiency-curve points in units of I_dep (narrowest cross-section). Each cell gives the voltage-trigger value, then
the sustained value.

| | i_10 | i_50 | i_95 | i_max | selectivity f = 0.35 |
|---|---|---|---|---|---|
| plain wire, 200 nm | 0.840 / 0.840 | 0.845 / 0.845 | 0.850 / 0.850 | 0.850 / 0.850 | 0.118 / 0.118 |
| long zipper | 0.777 / 0.788 | 0.795 / 0.805 | 0.853 / 0.853 | 0.907 / 0.907 | 0.129 / 0.118 |
| 100 nm strip | — | 0.766 / 0.825 | 0.775 / 0.825 | — | 0.170 / 0.109 |
| 56 nm strip | — | 0.603 / 0.823 | 0.614 / 0.823 | — | 0.282 / 0.060 |
| 30 nm strip | — | 0.327 / 0.815 | 0.334 / 0.815 | — | 0.451 / 0.008 |

Files: `results/hitmap/summary_conserving_T400_m1.4.txt`,
`results/narrow_baseline/width_scan_summary_conserving_T400_m1.4.txt`, and the same names under
`results/v2_voltage_trigger/`.

* **The plain wire is the same under both criteria.** Every photon that fires the trigger in the plain wire also
  leaves it resistive.
* **The zipper never beats the plain wire on the plateau.**
  - i_95 is 0.853 against 0.850 under both criteria, and the zipper has a tail up to 0.907.
  - It turns on 0.040–0.050 I_dep earlier at i_50.
  - Go / no-go is 0/3 with the sustained criterion and 1/3 with the trigger. Criterion C3 passed by 0.001 there.
* **The narrow-strip control is not decisive.**
  - With the trigger, a 56 nm strip beats the zipper by about 0.2 I_dep on both onset and plateau.
  - With the sustained criterion the zipper's onset is lower than every strip's (0.805 vs ≥ 0.815) and its plateau is
    0.03 higher.
  - Neither bound is the real device. A comparison with the narrow strips needs a heat model.
* **Some of the zipper's onset gain may come from the bridges.** Under the sustained criterion the zipper's onset
  lies below the floor of the plain strips. Its 12 nm bridges and trapped flux may keep a resistive state going at
  lower current than a plain strip of any width. That part of the gain says less about detection than about how
  easily the zipper stays resistive. It is also the part that depends most on the mesh (see Flux trapping).
* All strips and the zipper switch dark at I_sw/I_dep = 0.968–0.971. The 30 nm strip at 0.8 ξ agrees with 1.4 ξ under
  both criteria (sustained i_95 0.839 I_sw at both).
* **Latency** at threshold (t_trigger − t0, sustained criterion):
  - plain wire: median 92 τ0, max 129 τ0;
  - zipper: median 118 τ0, max 342 τ0, with 15 of 42 sites above 200 τ0.
  - With the trigger alone the zipper's median was 174 τ0.

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
  - With the sustained criterion that switching current is 0.909 I_dep for one quantum and 0.871–0.874 for two,
    against 0.968 for the empty device. Producing the same state from different photons gives the same value.
  - [0 −2] switches at 0.786 I_dep (1.4 ξ).
    - At 0.782 its flux leaves the left hole, moves through the right column for ~250 τ0 and comes back as
      [0 −2]: 1.3 slips in the last 50 τ0, 0.7 after settling, so it is not counted.
    - At 0.786 the same rearrangement ends resistive: 10 slips in the last 50 τ0 of the continuation.
    - With the trigger alone, this state "switched" at 0.775 on the rearrangement itself.
* **The zipper's onset is flux-mediated.**
  - At 0.78 I_dep, photons in the left channel and on the centre line (24 % of the area) fire the voltage trigger
    through a rearrangement like this one. The device is left quiet with new trapped flux, so none of it counts as
    sustained (0 %).
  - From 0.79 the same sequence ends in a resistive state, after a zero-voltage dwell with the flux trapped. The
    resistive onset comes 180–210 τ0 after the photon at 0.79, and 80–250 τ0 over all via-trap runs.
  - At 0.80 I_dep, 45 % of the area switches and stays resistive (52 % fires), all of it via trapped flux.
  - At 0.82: 57 % sustained, 44 % prompt, 1 % direct. At 0.85 all four are 92–94 %.
  - Direct switching starts at 0.81–0.85 in the channels.
  - The plain wire never traps flux, and every detection is prompt and sustained. Counting only prompt detections,
    the zipper's onset gain shrinks from 0.05 to 0.027 I_dep, and its plateau (i_95) comes 0.008 I_dep after the plain
    wire's.
  - The scan and the hitmap agree: every site first switches at the grid bias just above its hitmap threshold
    (42 of 42 zipper sites, 8 of 8 plain-wire sites).

  Efficiency-curve points in units of I_dep. The scan grid is 0.01, so values can be up to 0.01 above the bisected
  hitmap thresholds.

  | | i_10 | i_50 | i_95 | i_max |
  |---|---|---|---|---|
  | zipper, any trigger (lower bound) | 0.780 | 0.800 | 0.855 | 0.910 |
  | zipper, sustained (upper bound, current criterion) | 0.790 | 0.810 | 0.856 | 0.910 |
  | zipper, prompt | 0.811 | 0.823 | 0.858 | 0.910 |
  | zipper, direct (no trapped dwell) | 0.827 | 0.834 | 0.859 | 0.910 |
  | plain wire, all four | 0.840 | 0.850 | 0.850 | 0.850 |

* **A trapped state changes the response to the next photon.** I measured the photon threshold at four sites,
  starting from a trapped state instead of the empty device:
  - With the sustained criterion, from [0 −2] every site switches at 0.786. That is the state's own switching
    current, so the photon no longer matters. Against the empty device this is 0.001–0.028 lower.
  - One or two quanta shift the threshold by −0.005 to +0.057, mostly upwards. The largest rise, +0.057, is in the
    left channel next to [0 −1 −1]. Trapped flux makes the zipper less sensitive at the sites it blocks.
  - With the trigger alone, [0 −2] seemed to lower the threshold to 0.63. Those were rearrangements of the trapped
    flux, not switches.
  - So the zipper's response depends on what it absorbed before.
* **Mesh.** At 0.8 ξ, the (−68, 0) photon at 0.77 I_dep leaves the same [0 −2] state.
  - It switches at 0.845 I_dep, against 0.786 at 1.4 ξ. At 0.841 its flux already moves to the right column
    ([0 0 0 0 1 −1 0]) without switching. With the trigger alone the two meshes gave 0.795 and 0.775.
  - The 12 nm bridges are only 2.4 ξ wide: about 2 mesh cells at 1.4 ξ and 3 at 0.8 ξ. The thresholds that depend on
    flux crossing a bridge are therefore not converged. They move by 0.02 (trigger) to 0.06 (sustained) between the
    two meshes. The photon thresholds of plain strips agree within 0.007.
  - The zipper's via-trap onset (0.79–0.81) is exactly this kind of quantity, so it may sit higher on a finer mesh.
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

0. **Close the gap between the two bounds: add heat.** Couple the TDGL to a local electron-temperature balance:
   Joule heating, in-plane diffusion and escape to the substrate (Vodolazov, PR Applied 7, 034014 (2017); Zotova &
   Vodolazov, PRB 85, 024509 (2012)). T(r, t) then sets the local ε = 1 − T/T_c each step. Detection becomes "the
   resistive state latches", and the voltage-trigger / sustained bounds reduce to one threshold. This decides the
   narrow-strip comparison.
1. **Mesh at the bridges.** 0.5 ξ (2.5 nm, ~5 cells across a 12 nm bridge) for the [0 −2] switching current and a
   few left-channel and centre sites at 0.76–0.84 I_dep. If 0.5 ξ is too costly on the full strip, use one hole
   period with periodic-like long leads.
2. **Flux trapping, follow-up.** The first map is done (see above).
   - Photon sequences at sub-threshold bias: does accumulated flux reach a state that switches? Does a switching
     event clear the trapped flux?
   - A longer strip with more hole periods (end effects).
3. Hole-length sweep of the zipper (channel length / width from 2 to ~6): separates "current bypasses the short
   channel" from "detection needs the whole width to switch" as the reason the 56 nm channels do not act as narrow
   strips.
4. Re-run the column geometries of the archived cascade scripts (2-channel, aligned 3-channel, short zipper) with the
   current protocol (conserving hotspot, 0.8 ξ, full grid).
5. Hotspot with an energy sink (loss time 5–20 τ0 and ∞) on the zipper and the 56 nm strip.
6. Detection criterion: both bounds are now computed (voltage trigger in `results/v2_voltage_trigger/`, sustained
   in `results/`). Still to do: repeat the key thresholds at 250 τ0.
