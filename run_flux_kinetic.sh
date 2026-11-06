#!/usr/bin/env bash
# Flux-trapping study with the kinetic hotspot (5 Nov 2026).  Resumable: finished runs are skipped.
# Same sites as the conserving scans of 29 Oct / 2 Nov; Slot 400 and the plain wire get a wider, finer bias grid.
set -u
PY=${PY:-.venv/bin/python}
W=${WORKERS:-9}
export SNSPD_MESH_DIR=${SNSPD_MESH_DIR:-meshes}
$PY -m snspd.sim pin-meshes "$SNSPD_MESH_DIR"
COMMON="--hotspot kinetic --tmax 400 --mesh-xi 1.4 --workers $W"
step() { echo; echo "######## $(date +%H:%M:%S)  $*"; }
step "1. long zipper, full grid"
$PY fluxtrap.py scan $COMMON --geom Long_Zipper
step "2. Slot 400, the 12 sites of the conserving scan"
$PY fluxtrap.py scan $COMMON --geom Slot_400 --hitmap-dir results/slots --biases 0.70:0.74:0.02,0.75:0.88:0.01 \
    --sites="-96,0;-74,0;-58,0;-42,0;-19,0;-74,190;-42,190;-74,200;-42,200;-74,210;-58,210;-74,220"
step "3. plain wire (control: nothing to trap)"
$PY fluxtrap.py scan $COMMON --geom Plain_Wire --biases 0.78:0.90:0.01
step "done"
