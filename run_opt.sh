#!/usr/bin/env bash
# Optimisation candidates, 7 Nov 2026 (fabrication limit: every feature >= 40 nm), against the double column with
# 400 nm slots.  Double_800, Double_1200: fewer bridge rows (the bridge rows set the heated plateau).  Single_400_w120 / Single_800_w120: two 40 nm channels (40 | 120 | 40,
# free width 80 nm, cascade limit 0.50).  Plain_80: their same-I_c reference.  Then the heated snapshots.
# Preset nbtin, 2.0 xi, shortened strips.  Resumable.  Run from the repository root.
set -u
PY=${PY:-.venv/bin/python}
W=${W:-11}
export SNSPD_MESH_DIR= SNSPD_NO_PROGRESS=1
C="--hotspot kinetic --preset nbtin --tmax 400 --criterion trigger --isw-sustained --mesh-xi 2.0"
OPT=Double_800,Double_1200,Single_400_w120,Single_800_w120
step() { echo; echo "######## $(date '+%d %b %H:%M:%S')  $*"; }

step "1. trigger maps, critical grids: candidates and references"
$PY hitmap.py run $C --geoms Plain_Wire,Plain_80,$OPT --grid critical --no-probes --workers $W --out results/opt
$PY hitmap.py report $C --geoms Plain_Wire,Plain_80,$OPT --grid critical --out results/opt
step "2. heat, critical grids"
$PY heat_check.py run --workers $W --only-grids --grid-geoms Plain_80,$OPT
step "3. heated snapshots"
$PY heat_snapshots.py --cases "Plain_Wire:-58,0:0;Double_400:-80,0:0;Double_400:-40,220:100" --workers 6
step "done"

step "4. Double_600 (between the best heated slot length, 400 nm, and 800 nm): trigger and heat, critical grid"
$PY hitmap.py run $C --geoms Plain_Wire,Double_600 --grid critical --no-probes --workers $W --out results/opt
$PY heat_check.py run --workers $W --only-grids --grid-geoms Double_600
step "done (step 4)"

step "5. slot-end zone of the long double columns, with heat (y = b - 40 and b, edge and centre channel)"
$PY heat_check.py run --workers $W --extra "Double_600:-80,260;Double_600:-80,300;Double_600:0,260;Double_600:0,300;Double_800:-80,360;Double_800:-80,400;Double_800:0,360;Double_800:0,400;Double_1200:-80,560;Double_1200:-80,600;Double_1200:0,560;Double_1200:0,600"
step "done (step 5)"
