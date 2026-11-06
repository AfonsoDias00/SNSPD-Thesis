#!/usr/bin/env bash
# Rerun of the main comparisons with the "nbtin" parameter set (6 Nov 2026, snspd.sim.PRESETS): gamma = 10 and the
# hotspot times measured on 9 nm NbTiN (D 7.1 nm^2/tau0, escape 530 tau0), photon runs started from the settled
# current-carrying state.  Results sit next to the gamma = 1 ones under the tag kinetic-nbtin_T400_m1.4[_trigger].
# Every step is resumable: rerunning skips finished thresholds.  Run from the repository root.
#
# Two lanes run at once, sized to memory (~0.55 GB per worker for the 0.8 um devices, ~1 GB for the 2 um slot devices):
#   lane A (WA workers): plain wire + long zipper, narrow strips and width scans (both criteria), then slots, trigger
#   lane B (WB workers): slots, sustained
# The current-crowding maps are steady states without a photon; the steady state is a stationary GL solution and does
# not depend on gamma or the hotspot, so they are copied from the gamma = 1 kinetic runs.  The dark I_sw is recomputed:
# it is a dynamical instability within the window, and gamma changes it.
set -u
PY=${PY:-.venv/bin/python}
WA=${WA:-5}
WB=${WB:-6}
export SNSPD_MESH_DIR=${SNSPD_MESH_DIR:-meshes}
export SNSPD_NO_PROGRESS=1
COMMON="--hotspot kinetic --preset nbtin --tmax 400 --mesh-xi 1.4"
SLOTS="--geoms Plain_Wire_2um,Slot_110,Slot_200,Slot_400 --grid critical --out results/slots"
step() { echo; echo "######## $(date '+%d %b %H:%M:%S')  $*"; }

copy_crowding() {   # dir: crowding maps of the gamma = 1 kinetic runs -> the nbtin tag (both criteria)
  local d=$1 f
  for f in "$d"/crowd_kinetic_T400_m1.4*.json "$d"/crowd_*_kinetic_T400_m1.4*.npz; do
    [ -e "$f" ] || continue
    local g=${f/_kinetic_T400/_kinetic-nbtin_T400}
    [ -e "$g" ] || cp "$f" "$g"
  done
}
copy_crowding results/hitmap
copy_crowding results/slots

lane_a() {
  step "A1. plain 200 nm wire vs long zipper, full grid, sustained (upper bound)"
  $PY hitmap.py run $COMMON --workers $WA --criterion sustained
  step "A2. plain wire vs long zipper, full grid, voltage trigger (lower bound)"
  $PY hitmap.py run $COMMON --workers $WA --criterion trigger
  for crit in sustained trigger; do
    step "A3. narrow strips 30/56/100 nm, $crit"
    for w in 30 56 100; do $PY narrow_baseline.py run $COMMON --workers $WA --criterion $crit --width-nm $w; done
    suf=$([ $crit = trigger ] && echo _trigger || echo "")
    $PY narrow_baseline.py scan $COMMON --criterion $crit --compare hitmap_kinetic-nbtin_T400_m1.4$suf.csv
  done
  step "A4. 2 um plain wire vs slot columns, critical grid, voltage trigger"
  $PY hitmap.py run $COMMON --workers $WA --criterion trigger $SLOTS
  step "lane A done"
}

lane_b() {
  step "B1. 2 um plain wire vs slot columns, critical grid, sustained"
  $PY hitmap.py run $COMMON --workers $WB --criterion sustained $SLOTS
  step "lane B done"
}

lane_a > results/nbtin_lane_A.log 2>&1 &
lane_b > results/nbtin_lane_B.log 2>&1 &
wait
echo "all done $(date)"
