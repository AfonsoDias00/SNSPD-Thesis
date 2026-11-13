# SNSPD with a double slot column (pyTDGL)

Electro-thermal TDGL simulations of single-photon detection in a 200 nm strip of 9 nm NbTiN with two columns of
400 nm slots (three 40 nm channels in parallel), against plain strips.

Branch `double400-models` (from 10 Nov 2026): one device, and a photon, heat and circuit model built from the
literature. Earlier work is on other branches: `slot-geometries` (single, double and staggered slot columns, the
previous hotspot and heat models) and `repo-cleanup` (hole arrays, commit `00f9b0d`).

## Read first

- [docs/model.md](docs/model.md): what is solved, where each piece comes from, what is still assumed
- [docs/double400.md](docs/double400.md): the device and what is known about it
- [docs/methods.md](docs/methods.md): mesh, time step, settled start, detection criterion, with the evidence
- [docs/references.md](docs/references.md): the papers, and how much of each was read

## Layout

```
snspd/          shared code
  units.py      the film and pyTDGL's units (tau0, I_u, V_u, L_sq, j_dep)
  geometry.py   the double column (double_column(slot_nm)) and the plain reference strips
  thermal.py    electron and phonon temperatures, the photon as an energy deposit, epsilon(T_e)
  circuit.py    series kinetic inductance and 50 Ohm load
  sim.py        device, settled start, one event, the count criterion, threshold bisection
  cli.py        model options shared by the scripts
  analysis.py, plotting.py, store.py, fluxoid.py
event.py        one photon: traces of current, voltage, temperatures, and snapshots
thresholds.py   threshold bias for a count at each site, dark switching, efficiency
validation/     checks against known answers (python -m validation list)
results/
  baseline_kinetic/   thresholds of the previous model for Double_400 and the plain strips
  events/, thresholds/, validation/
```

## Running

Use the project environment (`.venv/bin/python`, pyTDGL 0.9.0).

```bash
.venv/bin/python -m validation run quick
```

```bash
.venv/bin/python event.py --cases "Double_400:-80,0:0.80;Plain_Wire:-58,0:0.90" --mesh-xi 2.0 --tmax 500 --workers 2
```

```bash
.venv/bin/python thresholds.py run --geoms Plain_Wire,Double_400 --grid critical --mesh-xi 2.0 --workers 11
```

- A case is `<geometry>:<x>,<y>:<bias>`, with x, y in nm and the bias in I_dep of the narrowest free
  cross-section (120 nm for the double column); `<geometry>:dark:<bias>` has no photon.
- **`--nsq`** is the series kinetic inductance in squares of film (default 750 = a 200 nm × 150 µm strip). It
  describes your device: set it to your wire length. `--source` replaces the circuit by a current source.
- Model variants go in their own files: `--thermal latent=False --label nolatent`,
  `--photon tau_th_ps=0.6 --label fast`, `--wavelength 1064`, `--t0 0.2`.
- `thresholds.py` is resumable: rows are appended as they finish and a rerun skips them. `--grid full` is the
  63-site map of the double column.
- **Cost.** At 2.0 ξ one 500 τ0 event of the double column takes about 2 min on a free core; a threshold is 7
  events. Run time scales roughly as (mesh edge)⁻⁴; 2.0 ξ is for screening only. More than about 11 workers on this
  machine slows every task.
- **Pin the meshes** before comparing across machines: `python -m snspd.sim pin-meshes meshes 1.4`, then run with
  `SNSPD_MESH_DIR=meshes`.

## Geometries

`Double_400` (the device), `Plain_Wire` (200 nm), `Plain_120` (the free width), `Plain_40` (one channel).
`Double_<L>` builds the double column with L nm slots on the fly (e.g. `Double_800`).
