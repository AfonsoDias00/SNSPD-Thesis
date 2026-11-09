#!/usr/bin/env bash
# Slot layouts, 7 Nov 2026 (overnight).  Preset nbtin (gamma 10, NbTiN hotspot times, settled start), kinetic hotspot,
# 400 tau0, voltage trigger (at gamma 10 the two bounds nearly coincide), dark I_sw from the sustained criterion.
# Shortened devices (150 nm of plain lead beyond the outermost slot): mesh_check.py showed they reproduce the 2 um
# strips within 0.011 I_dep at 1.4 xi and 0.019 at 2.0 xi.  Resumable.  Run from the repository root.
set -u
PY=${PY:-.venv/bin/python}
W=${W:-10}
export SNSPD_MESH_DIR= SNSPD_NO_PROGRESS=1
C="--hotspot kinetic --preset nbtin --tmax 400 --criterion trigger --isw-sustained"
REF=Plain_Wire,Plain_120,Plain_60,Plain_40
L400=Single_400,Double_400,Stagger_400
L150=Single_150,Double_150,Stagger_150
step() { echo; echo "######## $(date '+%d %b %H:%M:%S')  $*"; }

step "1. full maps at 2.0 xi: plain references and the 400 nm layouts (with the f = 0.6 / 0.35 probes)"
$PY hitmap.py run $C --mesh-xi 2.0 --geoms $REF,$L400 --workers $W --out results/maps
$PY hitmap.py report $C --mesh-xi 2.0 --geoms $REF,$L400 --out results/maps
step "2. critical grids at 2.0 xi: the 150 nm layouts (slot length)"
$PY hitmap.py run $C --mesh-xi 2.0 --geoms Plain_Wire,$L150 --grid critical --no-probes --workers $W --out results/maps150
$PY hitmap.py report $C --mesh-xi 2.0 --geoms Plain_Wire,$L150 --grid critical --out results/maps150
step "3. critical grids at 1.4 xi: the 400 nm layouts (mesh confirmation)"
$PY hitmap.py run $C --mesh-xi 1.4 --geoms Plain_Wire,$L400 --grid critical --no-probes --workers $W --out results/maps
$PY hitmap.py report $C --mesh-xi 1.4 --geoms Plain_Wire,$L400 --grid critical --out results/maps
step "4. snapshots at 1.4 xi, 0.01 I_dep above and below the threshold"
$PY snapshots.py --cases "Plain_Wire:-58,0:0;Single_400:-70,0:0;Double_400:-80,0:0;Stagger_400:-80,0:100" \
    --hitmap-dir results/maps --preset nbtin --criterion trigger --mesh-xi 1.4 --workers 8 --out docs/fig
step "5. flux scans at 2.0 xi, critical grids, 0.50-0.95 I_dep"
for g in Double_400 Single_400 Stagger_400; do
  $PY fluxtrap.py scan --geom $g --grid critical --biases 0.50:0.95:0.03 --hotspot kinetic --preset nbtin \
      --tmax 400 --mesh-xi 2.0 --hitmap-dir results/maps --hitmap-criterion trigger --out results/flux --workers $W
done
step "done"
