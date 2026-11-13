# Baseline: thresholds of the previous model

Copied on 10 Nov 2026 from the branch `slot-geometries` (commit `b41b37d`), rows of `Double_400`, `Plain_Wire`,
`Plain_120` and `Plain_40` only. Model: the kinetic ε-hotspot at f = 1 (not tied to a photon energy), preset
`nbtin` (γ = 10, settled start), a current source, 400 τ0 window. Thresholds are in units of the dark switching
current I_sw in the `i_det` column and in µA in `I_det_uA`; divide `I_det_uA` by I_dep (60.6 µA for the double
column, 100.8 µA for the 200 nm wire) to compare with the new results.

| File | Content |
| --- | --- |
| `thresholds_trigger_m2.csv`, `thresholds_trigger_m1.4.csv` | voltage-trigger thresholds without heat (full grid at 2.0 ξ, critical grid at 1.4 ξ) |
| `thresholds_heat_m2.csv` | with the linearised two-temperature Joule heat, dt_max 0.025, column `i_heat` in I_dep |
| `isw_m2.json`, `isw_m1.4.json` | dark switching currents (sustained criterion) |
| `fluxscan_m2.csv`, `fluxscan_m2_summary.txt` | fixed-bias scan of the double column: clean / trapped / transient / switched, with the flux left in the slots |

The scripts that produced them (`hitmap.py`, `heat_check.py`, `fluxtrap.py`) are on `slot-geometries`.
