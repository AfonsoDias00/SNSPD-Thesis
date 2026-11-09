#!/usr/bin/env bash
# Heated thresholds again, with dt_max 0.025 (7 Nov 2026): the dt_max 0.2 runs had a local instability at small cells
# that the heat model turned into a hot spot seeding the latch (docs/methods.md).  NbTiN 2T, preset nbtin, 2.0 xi.
set -u
PY=${PY:-.venv/bin/python}
W=${W:-11}
export SNSPD_MESH_DIR= SNSPD_NO_PROGRESS=1
step() { echo; echo "######## $(date '+%d %b %H:%M:%S')  $*"; }
step "1. heat, critical grids: plain wire and double column (400 nm)"
$PY heat_check.py run --workers $W --only-grids --grid-geoms Plain_Wire,Double_400
step "2. heated snapshots"
$PY heat_snapshots.py --cases "Plain_Wire:-58,0:0;Double_400:-80,0:0;Double_400:-40,220:100" --workers 6
step "3. heat, critical grids: the other layouts and the optimisation candidates"
$PY heat_check.py run --workers $W --only-grids --grid-geoms Stagger_400,Single_400,Plain_80,Single_400_w120,Single_800_w120,Double_600,Double_800,Double_1200
step "4. slot-end zone of the long double columns"
EX=""; for g in "Double_600:300" "Double_800:400" "Double_1200:600"; do k=${g%%:*}; b=${g##*:}; for x in -80 0; do EX="$EX;$k:$x,$((b-40));$k:$x,$b"; done; done
$PY heat_check.py run --workers $W --extra "${EX#;}"
step "done"

step "5. edge notches at the bridge rows (Double_400_n20): dark I_sw, trigger and heat, critical grid"
$PY hitmap.py run --hotspot kinetic --preset nbtin --tmax 400 --criterion trigger --isw-sustained --mesh-xi 2.0 \
    --geoms Plain_Wire,Double_400_n20 --grid critical --no-probes --workers $W --out results/opt
$PY heat_check.py run --workers $W --only-grids --grid-geoms Double_400_n20
step "done (step 5)"

step "6. the same slot-end sites for the 400 nm double column (fair area weighting)"
$PY heat_check.py run --workers $W --extra "Double_400:-80,160;Double_400:-80,200;Double_400:0,160;Double_400:0,200"
step "done (step 6)"

step "7. full heated maps of the passing candidates: Double_800 and Double_1200"
$PY heat_check.py run --workers $W --only-grids --full-geoms Double_800,Double_1200
step "done (step 7)"
