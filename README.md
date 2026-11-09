# Perforated SNSPD simulations (pyTDGL)

TDGL simulations of photon detection in superconducting nanowire single-photon detectors (SNSPDs) with internal holes,
for fabrication-friendly slot layouts in 9 nm NbTiN.

This branch (`slot-geometries`, from 6 Nov 2026) starts the second phase. The first phase (ellipse zipper, 30 nm slot
columns, narrow strips; 16 Oct – 6 Nov) was removed. Its code, results and report are in commit `00f9b0d` on
`repo-cleanup`. What it established about the method is in [docs/methods.md](docs/methods.md).

## Layout

```
snspd/          shared code
  units.py      film parameters and units (K0, j_dep, tau0, Pearl length)
  geometry.py   strips and hole layouts, photon-hit grids, unit cells and symmetries; register new layouts here
  hotspot.py    photon hotspot models (kinetic after Simon et al.; conserving; legacy)
  sim.py        device and mesh, parameter presets, settled start, detection criteria, bisection, task runner
  analysis.py   threshold maps, efficiency curves, quantiles, selectivity
  fluxoid.py    trapped-flux meter
  thermal.py    electro-thermal TDGL prototype
  store.py, plotting.py
  meshgen.py    graded meshes (fine over the slotted section, coarse in the leads); measured not worth it
hitmap.py       photon-threshold maps and efficiency curves for a set of geometries (the first is the reference)
fluxtrap.py     flux trapping: fixed-bias scans, trapped states, traces
snapshots.py    channel currents, voltage and |psi| snapshots just above / below a site's threshold
heat_check.py   thresholds with Joule heat (snspd/thermal.py) at a few sites
mesh_check.py   strip length and mesh: benchmark and accuracy against the 2 um, 1.4 xi reference
slot_analysis.py  tables and figures of the slot layouts against the plain references
run_slots.sh    the pipeline of 7 Nov 2026 (maps, slot length, mesh check, snapshots, flux scans)
validation/     solver checks against known answers (python -m validation list)
docs/methods.md the method, its parameters and the checks behind them
docs/slots.md   the slot layouts: results and what they mean
docs/SNSPD geometries beyond hole arrays.md   literature review on hole arrays vs edge protection (27 Oct)
```

## Running

Use the project environment (`.venv/bin/python`, pyTDGL 0.9.0). Every run is resumable: rows are appended to CSVs as
they finish, and a rerun skips them.

```bash
.venv/bin/python hitmap.py run --geoms Plain_Wire,Double_400 --preset nbtin --tmax 400 --criterion trigger --isw-sustained --mesh-xi 2.0 --workers 10 --out results/maps
.venv/bin/python hitmap.py report --geoms Plain_Wire,Double_400 --preset nbtin --tmax 400 --criterion trigger --mesh-xi 2.0 --out results/maps
.venv/bin/python slot_analysis.py
.venv/bin/python -m validation run quick
```

- **Use `--preset nbtin`** (γ = 10, measured NbTiN hotspot times, settled start). The default preset `""` is the
  γ = 1 model of the first phase.
- **`--criterion trigger`** gives the lower bound and the default the upper (sustained) bound. In plain strips at
  γ = 10 they coincide. In the slot layouts they do not: a channel photon can dump the channel's current with a few
  phase slips and heal, leaving flux trapped (docs/slots.md §6–7). With the trigger criterion add
  `--isw-sustained`: at γ = 10 the bias-on transient fakes a trigger in dark runs.
- **Pin the meshes** of holed geometries once results exist: `python -m snspd.sim pin-meshes meshes`, then run with
  `SNSPD_MESH_DIR=meshes`. Triangulations of holed devices differ between machines.
- **Memory:** about 0.55 GB per worker for 0.8 µm strips and 1 GB for 2 µm strips at 1.4 ξ. On 11 physical cores,
  more than ~16 workers in total slows every task.
- **Run time** scales roughly as (mesh edge)⁻⁴. Use `--mesh-xi 2.0` for screening only (docs/methods.md §5).

## Registered geometries

- **Slot layouts** (`geometry.slot_layout`): `Single_400`, `Double_400`, `Stagger_400`, and the same with 150 nm
  slots. 200 nm strips with 40 nm bridges and a 120 nm free width: one column of 80 nm slots (60 | 80 | 60), two
  aligned columns of 40 nm slots (40 nm everywhere), or two staggered by half a period. Each strip is as short as
  possible: 150 nm of plain lead beyond the outermost slot.
- **Plain references**: `Plain_Wire` (200 nm, 0.8 µm), `Plain_Wire_2um`, `Plain_120` (the free width), `Plain_60` and
  `Plain_40` (the channel widths).

## Adding a layout

Add a builder to `snspd/geometry.py` that returns a `Geometry`. It needs the holes, a photon-hit grid over one unit
cell, the probe sites and the symmetry used to complete the map (`strip`, `column` or `zipper`). Then register it in
`GEOMETRIES`. Always compare against a plain strip of the same width and length, with the same preset and mesh.
