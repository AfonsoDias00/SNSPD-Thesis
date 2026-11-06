# Nanohole SNSPD simulations

TDGL simulations (pyTDGL 0.9) of superconducting-nanowire single-photon detector strips with **internal holes**: does
perforating the strip lower the bias needed to detect photons (the plateau) relative to the dark switching current?

This file is the working manual: layout, commands, model settings, current numbers. For the full explanation and how
the pieces fit together, start with `docs/report.md` (section 0).

## Layout

```
snspd/                  shared library, imported by every current script
  units.py              film parameters and the unit system (K0, j_dep, tau0, Pearl length)
  geometry.py           strips and hole layouts, photon-hit grids, unit cells  <- add new hole layouts here
  hotspot.py            photon hotspot models: kinetic (default, after Simon et al. 2025), energy-conserving FEM
                        diffusion (all results before 4 Nov 2026), legacy free-space Gaussian
  sim.py                device, detection criterion, bisection, current crowding, process-pool runner
  fluxoid.py            fluxoid (trapped flux quanta) of every hole from the phase winding of psi
  analysis.py           threshold maps, efficiency curves, quantiles, selectivity
  store.py              resumable CSV / JSON results
  plotting.py           figure style and standard figures
  thermal.py            electro-thermal TDGL (prototype, 5 Nov 2026): local temperature with Joule heating, diffusion,
                        escape and optional electron-phonon split, fed back into epsilon each step
hitmap.py               photon-threshold maps: hole geometries vs the plain 200 nm wire, go / no-go
narrow_baseline.py      plain strips of any width (control) and the width scan
checks.py               tests of the map assumptions: symmetry, edge band, trigger delay
fluxtrap.py             flux a photon leaves trapped in the holes: state scan, trapped-state I_sw, traces, and
                        `compare` of the two hotspot models; simple figures by default (--detailed for the old ones)
compare_hotspots.py     two hotspot models side by side: thresholds site by site, or the epsilon(r, t) fields
hotspot_sensitivity.py  kinetic-hotspot thresholds under other time constants and the times measured on 9 nm NbTiN
thermal_check.py        thresholds with Joule heat (snspd/thermal.py) at 14 sites, against the two bounds
presentation_figures.py figures for docs/presentation_outline.md from stored results (no simulation)
run_kinetic.sh          the kinetic-hotspot rerun of the main comparisons (resumable)
run_flux_kinetic.sh     the kinetic-hotspot flux scans: zipper, Slot 400, plain wire (resumable)
meshes/                 pinned production meshes (not in git; `python -m snspd.sim pin-meshes` rebuilds them)
results/slots/          2 um strips with one column of 30 nm slots, 40 nm bridges (hitmap.py), both criteria
results/hitmap/         plain wire vs long zipper            (hitmap.py)        } sustained-switching criterion,
results/narrow_baseline/ plain strips                         (narrow_baseline.py)} from 29 Oct 2026
results/fluxtrap/       flux-trapping scan, trapped states, traces (fluxtrap.py)  }  both hotspot models
results/thermal/        electro-thermal check (thermal_check.py)
results/v2_voltage_trigger/      results with the voltage trigger alone (28-29 Oct), see its README
results/v1_voltage_or_collapse/  all results with the first criterion (23-29 Oct), see its README
parked/edge_bank/       edge-bank strip: code + results, parked until the hole study is done
archive/                unreliable and superseded early work, see archive/README.md
docs/                   report.md (the full story, read this first), literature review
validation/             tests against known answers (hotspot solver, epsilon, disc vs hole, photon thresholds, fluxoid
                        meter, ring, single hole, scan statistics); results/validation/summary.txt
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
.venv/bin/python hitmap.py run --geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400 --tmax 400 --workers 9 --out results/slots
.venv/bin/python hitmap.py run --geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400 --tmax 400 --workers 9 --out results/slots --criterion trigger
.venv/bin/python hitmap.py selftest                              # hotspot-solver checks
.venv/bin/python -m validation run quick                         # validation tests without long TDGL runs (1a 1b 2a), ~6 min
.venv/bin/python -m validation run all --workers 6               # every group, ~40 min; results/validation/summary.txt
.venv/bin/python fluxtrap.py scan    --workers 9                 # every zipper site x 0.56-0.96 I_dep, 400 tau0
.venv/bin/python fluxtrap.py scan    --geom Plain_Wire --biases 0.80:0.92:0.01   # control
.venv/bin/python fluxtrap.py trapped --workers 6                 # I_sw of the trapped states, second photon
.venv/bin/python fluxtrap.py trapped --mesh-xi 0.8 --from=-68,0@0.77 --no-second   # one state, finer mesh
.venv/bin/python fluxtrap.py scan    --geom Slot_400 --hitmap-dir results/slots --biases 0.76:0.86:0.02 --sites="-96,0;-74,0;-58,0;-42,0;-19,0;-74,190;-42,190;-74,200;-42,200;-74,210;-58,210;-74,220" --workers 9
.venv/bin/python fluxtrap.py trapped --geom Slot_400 --hitmap-dir results/slots --place="0 -1 0" --no-second
.venv/bin/python fluxtrap.py trace   --site=-68,0 --biases 0.70,0.77,0.80,0.85 --workers 4
.venv/bin/python fluxtrap.py compare [--geom Slot_400 --hitmap-dir results/slots]   # conserving vs kinetic scans
WORKERS=9 bash run_flux_kinetic.sh                               # the kinetic flux scans (sets SNSPD_MESH_DIR)
.venv/bin/python thermal_check.py run --workers 3                # thresholds with heat; then `report`
.venv/bin/python hotspot_sensitivity.py run --variants NbTiN --workers 4   # measured NbTiN times; `table` to print
```

At 1.4 xi a 400 tau0 run takes ~16 s on one core. The laptop CPU (Core Ultra 7 155H) has 6 performance cores and 10
efficiency cores; more than ~9 workers, or other jobs alongside, slow every run down sharply.

Result files are tagged `<hotspot>_T<tmax>_m<mesh_xi>`, e.g. `summary_conserving_T150_m0.8.txt`.

**Hotspot model.** Since 4 Nov 2026 the default is `--hotspot kinetic`. Every result before that date, and the
commands above as they were run then, used `--hotspot conserving`; add that flag to reproduce or extend them. The
validation suite checks the stored conserving results and keeps that model (`validation/common.py`).

**Pinned meshes.** The long zipper meshes differently on some machines (10 543 instead of the production 10 710
nodes at 1.4 ξ; plain strips and slot columns are identical), and flux-mediated zipper thresholds move by ~0.04
I_dep with the mesh. To extend stored results, use the production meshes:

```bash
.venv/bin/python -m snspd.sim pin-meshes meshes                  # extract them from the stored crowding files
SNSPD_MESH_DIR=meshes .venv/bin/python hitmap.py run --tmax 400 --mesh-xi 1.4 --workers 9
WORKERS=9 bash run_kinetic.sh                                    # the kinetic rerun; sets SNSPD_MESH_DIR itself
.venv/bin/python compare_hotspots.py [--criterion trigger]      # conserving vs kinetic thresholds, site by site
.venv/bin/python compare_hotspots.py --profiles results/hotspot_profiles_conserving_vs_kinetic.png
```

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
* **Hotspot** (`snspd/hotspot.py`). Both models share the energy transport: the photon's energy starts as a Gaussian
  of radius r0 = 14 nm and diffuses with D = 10 nm²/τ0 inside the film, reflected at edges and hole walls. Photon
  at t0 = 20 τ0.
  - **kinetic** (default, after Simon et al., PRB 112, 174512 (2025), arXiv:2501.13791): the suppression
    s = 1 − ε builds up over t_rise = 10 τ0 (the paper's thermalisation time τ_th ≈ 1.4–1.8 ps), energy escapes to
    the substrate over t_esc = 100 τ0 (τ_esc = 15 ps), and ε ≥ −0.25. For f = 1 no normal core forms (ε ≥ 0.43 at
    the centre of a 200 nm strip) but the suppression is broader and lasts longer.
  - **conserving** (all results before 4 Nov 2026): ε = 1 + A f e^(−t/t_relax) n with A = −2.5, t_relax = 35 τ0.
    The core is normal (ε = −1.5) at once and relaxes within ~50 τ0.
  - What the paper fixes and what it leaves open is in the `snspd/hotspot.py` docstring; the comparison is under
    Current results.

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

### Kinetic hotspot after Simon et al. (4 Nov 2026, `--hotspot kinetic`, 400 τ0, 1.4 ξ, production meshes)

The hotspot now follows Simon et al. (PRB 112, 174512 (2025)). The suppression builds up over the thermalisation
time instead of appearing at once, energy escapes to the substrate, and ε stays ≥ −0.25 (`snspd/hotspot.py`, which
also lists what the paper leaves open). Dark I_sw and crowding maps involve no photon and are reused from the
conserving runs (`results/kinetic_reused_dark_results.sh`); only the photon thresholds are new. Figure:
`results/hotspot_profiles_conserving_vs_kinetic.png`; site-by-site table: `python compare_hotspots.py`.

**The hotspot itself** (f = 1, centre of the 200 nm wire):

| | conserving | kinetic |
|---|---|---|
| ε at the hit site | −1.45 at t0 | minimum +0.43 at t0 + 12.5 τ0 |
| normal core (ε < 0) | 580 nm², gone by t0 + 5.5 τ0 | never forms |
| ∫(1 − ε) dA at t0 + 5 / 20 / 50 / 100 τ0 | 1334 / 869 / 369 / 88 nm² | 508 / 1124 / 1021 / 629 nm² |

The kinetic photon is a shallow, broad, long-lived suppression rather than a short normal spot. The weight it puts
on the strip is 1.3× the conserving one at t0 + 20 τ0, 2.8× at t0 + 50 and 7× at t0 + 100.

**Thresholds**, in units of I_dep (narrowest cross-section), conserving → kinetic. Each cell gives the voltage
trigger (lower bound), then sustained switching (upper bound).

| | i_50 | i_95 | i_max |
|---|---|---|---|
| plain wire, 200 nm | 0.845 → 0.833 / 0.845 → 0.833 | 0.850 → 0.848 / 0.850 → 0.848 | 0.850 → 0.848 / same |
| long zipper | 0.795 → 0.779 / 0.805 → 0.814 | 0.853 → 0.835 / 0.853 → 0.835 | 0.907 → 0.891 / same |
| 100 nm strip | 0.766 → 0.737 / 0.825 → 0.818 | 0.775 → 0.772 / 0.825 → 0.819 | — |
| 56 nm strip | 0.603 → 0.617 / 0.823 → 0.801 | 0.614 → 0.622 / 0.823 → 0.804 | — |
| 30 nm strip | 0.327 → 0.382 / 0.815 → 0.775 | 0.334 → 0.390 / 0.815 → 0.778 | — |
| plain wire, 2 µm | 0.849 → 0.833 / 0.848 → 0.832 | 0.851 → 0.847 / 0.851 → 0.846 | 0.851 → 0.849 / 0.851 → 0.848 |
| slot column, 110 nm | 0.826 → 0.820 / 0.837 → 0.833 | 0.914 → 0.893 / 0.914 → 0.896 | 0.973 → 0.952 / 0.973 → 0.952 |
| slot column, 200 nm | 0.812 → 0.811 / 0.837 → 0.835 | 0.872 → 0.872 / 0.875 → 0.879 | 0.963 → 0.952 / 0.971 → 0.952 |
| slot column, 400 nm | 0.787 → 0.790 / 0.835 → 0.821 | 0.849 → 0.835 / 0.860 → 0.846 | 0.954 → 0.951 / 0.973 → 0.952 |

Files: `results/{hitmap,slots}/summary_kinetic_T400_m1.4{,_trigger}.txt`,
`results/narrow_baseline/width_scan_summary_kinetic_T400_m1.4{,_trigger}.txt`. The slot columns use the critical grid.

* **The plain wire barely moves**: i_50 −0.012, i_95 −0.002, the same under both criteria. Hits near the edges get
  0.016–0.019 easier and centre hits 0.005–0.008 harder: the wall reflects the energy, so an edge hit at −96 nm
  reaches ε = 0.21 against 0.43 for a centre hit. Latency at threshold grows from a median of 92 to 134 τ0.
* **The zipper's plateau now edges past the plain wire's, inside the hotspot uncertainty.** i_95 is 0.835 against
  0.848 under both criteria
  (0.853 against 0.850 before) and the tail is shorter (i_max 0.891, was 0.907). The margin is 0.011–0.012, well
  short of the 0.03 go margin; go / no-go is 1/3 under both criteria (C2 stress passes by 2 µA/µm).
* **The zipper's onset gain shrinks under the sustained criterion**, from 0.040 to 0.019. Its lowest sites
  (0.788–0.798 I_dep at 15 sites left of the centre, 13 of them flux-mediated and triggering 195–220 τ0 after the
  photon) detect at 0.806–0.819 instead, 38–54 τ0 after the photon. Sites triggering more than 200 τ0 after the photon drop
  from 15 of 42 to 1. With the trigger alone the onset gain stays
  (0.054). These flux-mediated thresholds are the ones that also move with the mesh (Flux trapping) and with the
  hotspot's time constants (sensitivity below), so treat them as uncertain by ~0.03 I_dep.
* **The sustained floor of the narrow strips is not a strip property.** The ~0.82 I_dep floor of the conserving runs
  becomes 0.775 (30 nm), 0.801 (56 nm) and 0.818 (100 nm): it depends on how long the hotspot keeps the strip
  suppressed. With it, the 56 nm strip beats the zipper on onset and plateau under both criteria.
* **Slot columns**: every slot site gets 0.003–0.029 easier (medians −0.016 to −0.027 under both criteria), about as
  much as the 2 µm plain wire off its centre (−0.019), so their standing against it barely moves. The trigger onset
  gain shrinks (400 nm column: 0.075 → 0.056 at i_50), the 400 nm column's plateau edge grows (i_95 +0.013 → +0.023
  with the trigger, −0.014 → −0.004 sustained) and the tails shorten (i_max 0.97 → 0.95). Every column is still
  no-go under both criteria.
* **Selectivity** (f = 0.35) of the zipper falls below the plain wire's under the sustained criterion (−0.038; it
  was +0.000) and stays near it with the trigger (+0.005, was +0.011).

**How much of this comes from the paper** (`results/kinetic_sensitivity.txt`). The paper fixes the shape of the
build-up and the two time constants in picoseconds. It does not fix τ0 (it needs the film's ρ_n), the amplitude, the
profile or a floor on ε. Eleven sites were rerun with t_rise / t_esc = 7 / 75 τ0 (τ0 anchored on D = 0.5 cm²/s),
7 / 50 τ0 (τ_esc = 10 ps, the paper's other value) and with no floor on ε:
* **Robust in sign:** no normal core forms for a centre hit; sustained thresholds of the narrow strips fall (30 nm by
  0.011–0.040, 56 nm by 0.005–0.021); edge hits of the plain wire get easier (0.005–0.024).
* **Not robust:**
  - the change at the plain-wire centre (−0.013 to +0.021);
  - the zipper's plateau-tail sites (−0.026 to +0.008), so its 0.012 plateau edge over the plain wire is inside the
    spread;
  - narrow-strip trigger thresholds (30 nm: −0.02 to +0.10);
  - the zipper's flux-mediated sites, which flip between 0.788 and 0.817 (a 7 τ0 rise brings the delayed channel
    back at (−96, 0)).
* **Consequence:** differences up to ~0.03 I_dep between geometries depend on hotspot details the paper leaves
  open. That is the size of the go margin. These include the zipper's onset gain and plateau edge, and the
  56 nm strip vs zipper order under the sustained criterion. Settling them needs the film's τ_th and ρ_n (or T_c),
  not more geometry runs. Larger gaps survive every variant, e.g. the narrow strips under the trigger.

### Measured NbTiN parameters and the hotspot (5 Nov 2026, `results/kinetic_sensitivity.txt`)

Sidorova et al. measured 9 nm NbTiN films from KTH, the group whose hairpins set λ (PRB 104, 184514 (2021),
arXiv:2108.04717; arXiv:2204.01423): T_c 9.51 K (film) / 8.60 K (200 nm strip), R_sq 381 Ω (ρ_n 343 µΩ cm),
D 0.472 cm²/s, τ_EP(T_c) 2.9–4.0 ps, phonon escape τ_esc 80 ps, electron energy relaxation 90 ps, c_e/c_ph 0.13–0.16,
I_C/I_dep 0.63 at 3 K. τ_th is not measured.
* τ0 = ħ/(πΔ0) = 0.145–0.16 ps, or μ0λ²/ρ_n = 0.13 ps: the 0.12–0.17 ps used so far stands.
* The escape time does not: 80 ps = 530 τ0, against the kinetic hotspot's 100 τ0 (15 ps, NbN).
* Variant `NbTiN` (t_esc 530 τ0, D 7.1 nm²/τ0, t_rise 10 τ0, no heat): every threshold drops by 0.013–0.10 I_dep
  (plain 0.810–0.813, zipper onset 0.790–0.804, zipper tail 0.809–0.835, 30 nm v3 0.677, 56 nm v3 0.738). The order
  of the geometries stays. It overstates the effect: most of the energy goes to the phonons within τ_EP.

### Flux trapping with the kinetic hotspot (5 Nov 2026, `run_flux_kinetic.sh`, `results/fluxtrap/`)

Same scans as below, with `--hotspot kinetic`: zipper 42 sites × 34 biases, Slot 400 12 sites × 17 biases
(0.70–0.88), plain wire 8 × 13 (0.78–0.90). Comparison: `results/fluxtrap/compare_<geom>_T400_m1.4.txt`; figure
`outcomes_<geom>_T400_m1.4.png` (share of the unit cell by outcome vs bias, one panel per hotspot).

| Zipper, share of the cell, conserving / kinetic | 0.70 | 0.78 | 0.80 | 0.82 | 0.84 |
| --- | --- | --- | --- | --- | --- |
| flux left, no signal | 59 / 25 | 68 / 47 | 43 / 26 | 30 / 6 | 7 / 2 |
| pulse, then heals | 0 / 0 | 24 / 43 | 7 / 57 | 10 / 1 | 1 / 0 |
| detected via trapped flux | 0 / 0 | 0 / 0 | 45 / 10 | 56 / 1 | 6 / 0 |
| detected directly | 0 / 0 | 0 / 0 | 0 / 0 | 1 / 87 | 84 / 95 |

* Less flux, later: first flux at 0.68 instead of 0.64 I_dep, 54 instead of 34 τ0 after the photon (median); 211
  instead of 319 runs end with silent flux. Same states, more double quanta ([0 −2]: 29 vs 10 runs).
* Switching via trapped flux: 11 runs instead of 68. Pulses that heal: 76 instead of 30, median 15 slips (was 3),
  up to 88 (≈ 80 eV at 0.80 I_dep).
* Zipper efficiency points from the scan, conserving → kinetic: trigger i_50 0.800 → 0.784, sustained i_50 0.810 →
  0.815 and i_95 0.856 → 0.840; direct-only i_50 0.834 → 0.820, so the direct-only onset gain over the plain wire is
  0.016 → 0.015.
* Six left-channel sites are resistive at 0.79, heal at 0.80–0.81 and are resistive again from 0.82. The hitmap
  bisection, which assumes a monotone outcome, found the upper crossing (0.817).
* Slot 400: still no detection via trapped flux; flux entry and pulses that heal cover 86–92 % of the sampled cell
  at 0.78–0.82; sustained onset 0.81–0.84. Plain wire: no flux, no pulses, every detection direct (0.83–0.85).
* Trapped states (`summary_Long_Zipper_kinetic_T400_m1.4.txt`): the same switching currents as with the conserving
  hotspot within 0.003 ([0 −1] 0.909, [0 +1 0 0 0 −1 0] 0.874, [0 −1 −1] 0.874, [0 −2] 0.786; new [0 −1 0 0 0 −1 0]
  0.824). All hold at zero bias, and a second photon shifts thresholds by −0.03 to +0.04 vs pristine.

### Electro-thermal check (5 Nov 2026, `snspd/thermal.py`, `thermal_check.py`, `results/thermal/`)

dθ/dt = D∇²θ + η|∇μ|² − θ/t_esc, ε = max(ε_photon − θ, −1), θ = (T − T0)/(T_c − T0), after Zotova & Vodolazov,
PRB 85, 024509 (2012). η = 4E_cond/(C_v(T_c − T0)) = 0.59 (one temperature, nominal; ρ_n cancels), varied ×0.5 and
×2; "NbTiN 2T" adds the phonons with the measured values above (η 1.22, D 7.1 nm²/τ0, τ_EP 27 τ0, c_e/c_ph 0.16,
escape 530 τ0). θ = 0 within 25 nm of the contacts. Detection: still resistive at 400 τ0 (bisection, 0.004 I_dep).

| Site | trigger (lower) | sustained (upper) | heat, 1T | heat, NbTiN 2T |
| --- | --- | --- | --- | --- |
| 30 nm (0, 0) | 0.376 | 0.775 | 0.457 | 0.439 |
| 56 nm (0, 0) | 0.620 | 0.801 | 0.613 | 0.609 |
| 100 nm (0, 0) / (−36, 0) | 0.774 / 0.732 | 0.819 / 0.817 | 0.756 / 0.723 | 0.759 / 0.717 |
| plain (0, 0) / (−96, 0) | 0.848 / 0.832 | same | 0.827 / 0.808 | 0.827 / 0.801 |
| zipper (0, 0) / (−68, 0) / (−96, 0) | 0.756 / 0.772 / 0.774 | 0.804 / 0.817 / 0.817 | 0.768 / 0.764 / 0.765 | 0.760 / 0.762 / 0.762 |
| zipper tail (−42, 58) / (−36, 61) | 0.856 / 0.891 | same | 0.835 / 0.870 | 0.835 / 0.870 |
| Slot 400 (−58, 0) / (−28, 200) | 0.762 / 0.807 | 0.808 / 0.830 | 0.754 / 0.788 | 0.747 / 0.782 |

* Heat merges the bounds at or below the lower one (30 nm and the zipper centre are the only sites above it).
* Narrow strips detect below every zipper site: the narrow-strip comparison goes to the narrow strips.
* Zipper onset sites gain ~0.05 over the plain wire; the tail stays 0.04 above it.
* η ×0.5 to ×2 moves thresholds by ≤ 0.03 (≤ 0.06 at 30 nm); the 2T model is within 0.008 of 1T (0.018 at 30 nm).
* Prototype: 14 sites, estimated heat parameters, no circuit.

### Slot column, 2 µm strips (30 Oct 2026, `results/slots/`, `_trigger` files are the lower bound)

The zipper's 16 nm holes and 12 nm bridges cannot be fabricated, and its 12 nm bridges are not mesh-converged. The new
layouts (`snspd/geometry.py`) are one centred column of slots in a 200 nm × 2 µm strip:
- slots 30 nm wide, with semicircular ends (a 400 nm long ellipse would have a ~1 nm tip radius);
- 40 nm bridges, two 85 nm channels, narrowest free width 170 nm (I_dep = 85.7 µA);
- lengths 110, 200 and 400 nm, with 7, 5 and 3 slots, so the channels are 1.3, 2.4 and 4.7 times as long as wide.

The plain reference is rerun at 2 µm. Probes sit at y = ±750 nm, outside the hole arrays.

| Device | I_sw [µA], lower / upper | i_50 lower / upper | i_95 lower / upper | i_max lower / upper | Go / no-go, lower / upper |
|---|---|---|---|---|---|
| Plain wire, 2 µm | 95.8 / 98.7 | 0.849 / 0.848 | 0.851 / 0.851 | 0.851 / 0.851 | reference |
| Slot 110 nm | 83.5 / 83.5 | 0.826 / 0.837 | 0.914 / 0.914 | 0.973 / 0.973 | 0/3 / 0/3 |
| Slot 200 nm | 83.5 / 83.5 | 0.812 / 0.837 | 0.872 / 0.875 | 0.963 / 0.971 | 0/3 / 0/3 |
| Slot 400 nm | 82.5 / 83.5 | 0.787 / 0.835 | 0.849 / 0.860 | 0.954 / 0.973 | 1/3 / 0/3 |

Thresholds are in units of I_dep (100.8 µA for the wire, 85.7 µA for the slots). The lower bound is the voltage trigger, the upper bound sustained switching.

- **The plateau is not lowered.** Even the best case, Slot 400 at the lower bound, is 0.002 I_dep below the wire, within the bisection tolerance. The zipper's i_95 is 0.853.
- **The tail sits in the bridges.** The highest thresholds, 0.98–1.0 I_sw, are at photon sites in the 40 nm bridges between slots, under both criteria. Shorter slots have more bridge area, which fits i_95 improving with length. Why bridges are insensitive (probably little current flows through them) was not checked directly.
- **Onset depends on the bound.** The i_50 gain over the wire is 0.011–0.013 I_dep at the upper bound and 0.023 / 0.037 / 0.062 at the lower bound for 110 / 200 / 400 nm. The gap between the bounds grows with length, the opposite of what was expected. A reduced flux scan explains the extra lower-bound events (below).
- **Peak current is higher than the wire's.** At i_95 it is 475–501 µA/µm against 452, so C2 fails everywhere. Slot 400 at the lower bound passes C3 only (selectivity +0.022).
- **Caveats.** The go / no-go tests use thresholds in units of I_sw. The wire's I_sw falls to 95.8 µA under the trigger, which inflates the lower-bound gains in those units: Slot 400's i_95 gain is +0.013 in I_sw but +0.002 in I_dep. The slots have no 0.8 ξ mesh check yet.

**Flux scan of Slot 400 (`results/fluxtrap/`, 12 sites × 0.76–0.86 I_dep, 72 runs; a subset, so its area fractions are indicative).** The onset gain at the lower bound is flux entering the slots, not a resistive switch.

| Bias [I_dep] | What most sites do | Phase slips | Ends |
|---|---|---|---|
| 0.76–0.80 | Flux enters a slot silently, then the second or third quantum fires the trigger | 0.6–2.5 | Quiet, 1–4 quanta left in the slots |
| 0.80–0.84 | Burst: the slot fills to 4–5 quanta, empties and refills with the opposite sign | 5–14 | Quiet, 2–3 quanta left |
| 0.84–0.86 | Sustained switching, with no zero-voltage dwell first | 80–1100 | Resistive |

All 24 transient runs end with flux in a slot and no voltage. All 19 sustained runs are direct, with none via a trapped state, unlike the zipper. Sustained onset agrees with the hitmap at 12 of 12 sites. At the 0.84–0.86 onset, 4 runs near the slot ends sit in a slower resistive state (13–26 slips per 50 τ0, against about 190), whose stability was not studied. A [0 −1 0] state placed directly (`fluxtrap.py trapped --place`, phase wound around the slot and relaxed 50 τ0) matches the one a photon at (−96, 0) leaves at 0.76: same fluxoid, kept for 400 τ0 at zero bias, and the same switching current, 0.962 I_dep against 0.973 for the empty device. So the upper bound, i_50 = 0.835 and a gain of 0.013, is the better reading unless heat makes the bursts latch. A burst dissipates about 5–13 eV, 6–16 times a 1550 nm photon, so a heat model is what decides this.

**Mesh coarser towards the terminals: tested, not used.** On the plain 2 µm wire, with the mesh edge growing from 1.4 ξ
at |y| = 700 nm to the value below at the terminals:

| Edge at the terminals | Nodes | I_sw [µA] | i_det at x = −96 / 0 nm [I_dep] | Run time |
|---|---|---|---|---|
| 1.4 ξ (uniform) | 27 757 | 98.7 | 0.853 / 0.843 | 100 % |
| 2 ξ | 25 277 (−9 %) | 97.9 (−0.8 %) | 0.855 / 0.841 | ~95 % |
| 3 ξ | 23 302 (−16 %) | 96.2 (−2.6 %) | 0.857 / 0.844 | ~91 % |

Photon thresholds stay within about the bisection tolerance (0.003), but the dark switching current falls, probably
because a coarse lattice carries less current than the continuum, which would make the leads the weak point (the
location was not checked). The time saved does not justify that. The 2 µm thresholds also agree with the 0.8 µm strip (0.850 / 0.842).

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
  - The 12 nm bridges are only 2.4 ξ wide: about 3 mean mesh edges at 1.4 ξ (4.1 nm) and 5 at 0.8 ξ (2.3 nm). The thresholds that depend on
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

Status of the 4 Nov plan: electro-thermal TDGL (prototype done, `thermal_check.py`); ρ_n / T_c / τ_th (literature
values for the same film class: τ0 confirmed, escape 5× longer, τ_th unknown); kinetic flux study (done); defect
baseline and barrier dark counts (assessed, not run: docs/report.md §12).

1. **Full maps with heat** (plain, zipper, 56 nm, Slot 400) with `snspd/thermal.py`, NbTiN 2T as default: i_50 and
   i_95 under one criterion.
2. **Calibrate η** on the measured retrapping current of 9 nm NbTiN strips (arXiv:2204.01423).
3. **Two-temperature photon:** the photon's energy in the same electron/phonon equations (τ_EP, c_e/c_ph measured).
4. **Our film:** R(T) of a bar from the same chip for T_c and R_sq; τ_th needs a time-resolved measurement.
5. **Circuit:** kinetic inductance and a 50 Ω load, to tell a latch from a pulse that resets.
6. **Mesh at the bridges:** 0.5 ξ for the [0 −2] switching current, if the zipper's flux numbers are quoted.
7. **Defect-limited baseline:** only meaningful once the absolute photon scale is fixed. Every f = 1 threshold of
   the 200 nm devices (≥ 0.75 I_dep even with heat) lies above the measured I_C/I_dep = 0.63 of this film.
8. **Barrier dark counts** (Vodolazov, PRB 85, 174507 (2012)): F0 = Φ0²d/(4πμ0λ²) ≈ 490 K; a straight film reaches
   ΔF ≈ 25 k_BT (2.5 K) only at ~0.96 I_dep. Needs a pinned-phase relaxation in the solver; matters only once a
   geometry wins on the plateau.
9. **Photon sequences** at sub-threshold bias (does accumulated flux switch? does a switch clear it?).
10. Hole-length sweep of the zipper (channel length / width 2 to ~6), and a longer strip with more hole periods.
11. Re-run the column geometries of the archived cascade scripts (2-channel, aligned 3-channel, short zipper) with
    the current protocol.
12. Repeat the key thresholds at 250 τ0 (window dependence of both bounds).
