#!/usr/bin/env bash
# The rest of lane A of run_nbtin.sh, started on 6 Nov 2026 in its place: the plain wire vs long zipper trigger maps
# WITHOUT the f = 0.6 / 0.35 probe sites (--no-probes; rerun without the flag to add them), then the narrow strips
# (with their probes, they are cheap) and the width scans.  Resumable.  Run from the repository root.
set -u
PY=${PY:-.venv/bin/python}
WA=${WA:-5}
export SNSPD_MESH_DIR=${SNSPD_MESH_DIR:-meshes}
export SNSPD_NO_PROGRESS=1
COMMON="--hotspot kinetic --preset nbtin --tmax 400 --mesh-xi 1.4"
step() { echo; echo "######## $(date '+%d %b %H:%M:%S')  $*"; }

step "A2'. plain wire vs long zipper, full grid, voltage trigger (lower bound), f = 1 only (--no-probes)"
$PY hitmap.py run $COMMON --workers $WA --criterion trigger --no-probes
for crit in sustained trigger; do
  step "A3. narrow strips 30/56/100 nm, $crit"
  for w in 30 56 100; do $PY narrow_baseline.py run $COMMON --workers $WA --criterion $crit --width-nm $w; done
  suf=$([ $crit = trigger ] && echo _trigger || echo "")
  $PY narrow_baseline.py scan $COMMON --criterion $crit --compare hitmap_kinetic-nbtin_T400_m1.4$suf.csv
done
step "lane A done (A4, the slots with the trigger, runs separately with --no-probes)"
