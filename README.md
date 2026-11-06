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
hitmap.py       photon-threshold maps and efficiency curves for a set of geometries (the first is the reference)
fluxtrap.py     flux trapping: fixed-bias scans, trapped states, traces
validation/     solver checks against known answers (python -m validation list)
docs/methods.md the method, its parameters and the checks behind them
docs/SNSPD geometries beyond hole arrays.md   literature review on hole arrays vs edge protection (27 Oct)
```

## Running

Use the project environment (`.venv/bin/python`, pyTDGL 0.9.0). Every run is resumable: rows are appended to CSVs as
they finish, and a rerun skips them.

```bash
.venv/bin/python hitmap.py run --geoms Plain_Wire_2um,<new layout> --preset nbtin --tmax 400 --workers 6
.venv/bin/python hitmap.py report --geoms Plain_Wire_2um,<new layout> --preset nbtin --tmax 400
.venv/bin/python -m validation run quick
```

- **Use `--preset nbtin`** (γ = 10, measured NbTiN hotspot times, settled start). The default preset `""` is the
  γ = 1 model of the first phase.
- **`--criterion trigger`** gives the lower bound and the default the upper (sustained) bound. At γ = 10 they nearly
  coincide (docs/methods.md §3).
- **Pin the meshes** of holed geometries once results exist: `python -m snspd.sim pin-meshes meshes`, then run with
  `SNSPD_MESH_DIR=meshes`. Triangulations of holed devices differ between machines.
- **Memory:** about 0.55 GB per worker for 0.8 µm strips and 1 GB for 2 µm strips at 1.4 ξ. On 11 physical cores,
  more than ~16 workers in total slows every task.
- **Run time** scales roughly as (mesh edge)⁻⁴. Use `--mesh-xi 2.0` for screening only (docs/methods.md §5).

## Adding a layout

Add a builder to `snspd/geometry.py` that returns a `Geometry`. It needs the holes, a photon-hit grid over one unit
cell, the probe sites and the symmetry used to complete the map (`strip`, `column` or `zipper`). Then register it in
`GEOMETRIES`. Always compare against a plain strip of the same width and length, with the same preset and mesh.
