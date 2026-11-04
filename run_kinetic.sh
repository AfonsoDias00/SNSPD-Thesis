#!/usr/bin/env bash
# Kinetic-hotspot rerun of the main comparisons (4 Nov 2026).  Every step is resumable: rerunning this script skips
# finished thresholds.  Run from the repository root; the dark I_sw and crowding results are reused from the
# conserving runs first (results/kinetic_reused_dark_results.sh).  On the laptop use --workers 6-9; here 2.
set -u
PY=${PY:-.venv/bin/python}
W=${WORKERS:-2}
COMMON="--hotspot kinetic --tmax 400 --mesh-xi 1.4 --workers $W"
step() { echo; echo "######## $(date +%H:%M:%S)  $*"; }

bash results/kinetic_reused_dark_results.sh > /dev/null

step "1. plain 200 nm wire vs long zipper, full grid, sustained (upper bound)"
$PY hitmap.py run $COMMON --criterion sustained
step "2. narrow strips 30/56/100 nm, sustained"
for w in 30 56 100; do $PY narrow_baseline.py run $COMMON --criterion sustained --width-nm $w; done
$PY narrow_baseline.py scan $COMMON --criterion sustained --compare hitmap_kinetic_T400_m1.4.csv
step "3. plain wire vs long zipper, full grid, voltage trigger (lower bound)"
$PY hitmap.py run $COMMON --criterion trigger
step "4. narrow strips, voltage trigger"
for w in 30 56 100; do $PY narrow_baseline.py run $COMMON --criterion trigger --width-nm $w; done
$PY narrow_baseline.py scan $COMMON --criterion trigger --compare hitmap_kinetic_T400_m1.4_trigger.csv
step "5. 2 um plain wire vs slot columns, critical grid, sustained"
$PY hitmap.py run $COMMON --criterion sustained --geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400 --grid critical --out results/slots
step "6. 2 um plain wire vs slot columns, critical grid, voltage trigger"
$PY hitmap.py run $COMMON --criterion trigger --geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400 --grid critical --out results/slots
step "done"
