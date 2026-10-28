# Archive

Frozen scripts and results that are **not** part of the current analysis. Nothing here is imported by the current
code. The scripts are kept as they were; they write into folders relative to the directory they are run from, so run
one only from inside its own folder. Large solver files (`*.h5`, about 1.3 GB in total, mostly in
`superseded/columns_legacy_hotspot/`) are ignored by git and can be deleted.

## unreliable/ — do not quote these numbers

### static_and_first_photon/  (8–16 Oct)

| Script | Results | Problem |
|---|---|---|
| `critical_current_sweep.py` | `results_critical_current/` | Every simulated point exceeded the 0.01 V0 threshold, so no sub-critical point was reached; the 340–406 µA "I_c" values are about a third of the depairing current and reflect the bias turn-on transient in a 20 τ0 window. |
| `order_parameter_current_sweep.py`, `four_hole_sweep.py` | `results_order_parameter_sweep_10nm/` | `four_hole_sweep.py` is a copy that writes the same `no_hole_10nm_*` files and CSV; the CSV (switching at 600–610 µA) is the four-hole device, not the plain wire. |
| `four_hole_test.py` | `results_four_holes/` | `max_edge_length=0.1` is 100 nm (the comment says 10 nm). |
| `current_analysis.py`, `hole_optimization_sweep.py` | `results_optimization/` | 100 nm mesh; the edge-minus-centre current never crossed zero and is noisy. |
| `hole_interaction_landscape.py` | `results_interaction_landscape/` | 100 nm mesh; qualitative crowding maps only. An annotated duplicate of the script sat in `results_four_holes/`. |
| `photon_vs_spontaneous_switching.py`, `elipse_analysis.py` | `results_photon_margin/`, `results_spatial_sensitivity/` | Strict late-mean voltage criterion (> 0.35 V0), 15 µA tolerance, γ = 10, 1.1 w terminals: the solid wire's photon threshold came out above its dark I_sw, and the hole devices' I_sw exceeded their own depairing current by 22 %. |

### static_hotspot_arrays/  (18–20 Oct)

`efficiency_curves.py`, `plot_staggered_vs_2x2_maps.py`, `staggered_complete_verification.py`,
`plain_vs_2x2_vs_staggered.py`, `slender_staggered_small_photon.py`, `channel_tests.py` and their results.
2×2, staggered 2-1-2 and first column-of-holes layouts with a static Gaussian hotspot that was switched on at t = 0
together with the bias step, with a different hotspot size and strength in almost every script (r = 11–18 nm,
A = −1 to −3). The plain-vs-holes ranking flips with the hotspot parameters, so these photon thresholds cannot be
compared with each other or with later work. Their I_sw values are close to depairing and are fine.

### narrow_baseline_floor070/

The first 56 nm narrow-wire rows (26 Oct) before the photon-threshold bracket was lowered: the f = 1 thresholds sat
on the 0.70 I_sw floor (`below_bracket`), so they were bounds, not values. Replaced by the rows in
`results/narrow_baseline/`.

## superseded/ — valid, but replaced by better-resolved work

### columns_legacy_hotspot/  (20–22 Oct)

`cascade_full_analysis.py` (plain vs 2 vs 3 channels), `cascade_zipper_analysis.py` (aligned vs zipper),
`cascade_long_ellipse_zipper.py` (short vs long zipper) share ~90 % of their code; `photon_hitmap.py` is the first 2D
threshold map of the long zipper; `tmax_sensitivity.py` checks the window length. All use the standardised legacy
hotspot (r0 = 14 nm, A = −2.5, D = 10 nm²/τ0, t_relax = 35 τ0, photon at 20 τ0) and a 1.4 ξ mesh, so the numbers are
valid, but the efficiency curves come from 4 hand-weighted sites. The current-redistribution traces
(`cascade_comparison_*.png`) are the one result not repeated in the current code. The column geometries (2-channel,
aligned 3-channel, short zipper) are defined in these scripts and can be moved into `snspd/geometry.py` when they are
re-studied.

Window check (`tmax_sensitivity.py`): from 60 to 150 τ0, I_sw drops by 1 % (plain) and 2.6 % (zipper); between 125 and
150 τ0 thresholds still move by ≤ 0.009 (plain) and ≤ 0.018 (zipper edge site).

### pre_refactor/

`photon_hitmap_v2.py`, `narrow_baseline.py`, `edge_bank.py` exactly as they produced `results/` and
`parked/edge_bank/results/` (23–27 Oct). They were replaced on 28 Oct by the `snspd/` library, `hitmap.py`,
`narrow_baseline.py` and `parked/edge_bank/edge_bank.py`, which regenerate every summary from the stored data with
identical numbers (narrow-strip and edge-bank quantiles within 0.001 because the threshold profile is now sampled at
201 instead of 1121 points) and reproduce live thresholds within the bisection tolerance.

## exploration/

`baseline.ipynb` (first pyTDGL run, ξ = 11 nm, λ = 125 nm), `mesh_interactive.ipynb` (mesh viewer) and the 8 Oct
single-hole output it reads (`results_oct08_single_hole/`, written by `tests/single_hole_test.py`).
