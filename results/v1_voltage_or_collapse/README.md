# Results with the first detection criterion (23-29 Oct 2026)

Everything here was simulated with the detection criterion "probe voltage (peak or sustained mean) **or** collapse of
the supercurrent through one monitored cross-section" (the right half at y = 0 for plain strips, the zipper's right
channel at y = 61 nm). The collapse branch was removed on 29 Oct 2026: a photon inside the zipper's monitored channel,
site (68, 61), triggered it at i_det = 0.724 with no voltage at all, against 0.801 at its glide image (see
`checks/symmetry_summary_conserving_T150_m0.8.txt`). The maps themselves only simulate x <= 0, away from that channel,
so the numbers below are expected to change little, but they are not comparable row by row with voltage-only runs.

The conclusions they support (to be confirmed with the voltage-only re-run):

* long zipper vs plain wire: earlier onset (i_50 0.025-0.045 I_dep lower), later plateau (i_95 0.025-0.03 I_dep
  higher); go / no-go 0/3
* plain strips: i_50 = 0.33 / 0.61 / 0.78 / 0.86 at 30 / 56 / 100 / 200 nm; the 56 nm strip beats the zipper
* window: plain-wire photon thresholds are window-independent; several zipper channel and edge sites switch late
  (150-370 tau0) and fall with the window (T150 / T250 / T400 runs in `hitmap/`, delay curves in `checks/`)

Reports regenerate with the current code, e.g.
`python hitmap.py report --mesh-xi 0.8 --out results/v1_voltage_or_collapse/hitmap`.
