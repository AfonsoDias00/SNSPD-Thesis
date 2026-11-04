#!/usr/bin/env bash
# The kinetic-hotspot runs (4 Nov 2026) reuse every photon-independent result of the conserving runs: the dark
# switching current I_sw (bisected with no hotspot) and the steady-state current-crowding maps (no photon).  Neither
# depends on the hotspot model, so they are copied under the kinetic tag instead of being recomputed.  Only the photon
# thresholds are new.  Run from the repository root.
set -eu
cp_tag() {  # src dst
  [ -e "$2" ] || { cp "$1" "$2"; echo "copied $1 -> $2"; }
}
reuse() {  # dir src_suffix dst_suffix  (suffix = tag after the stem, e.g. conserving_T400_m1.4)
  local d=$1 s=$2 t=$3
  for f in "$d"/isw_"$s".json "$d"/crowd_"$s".json "$d"/crowd_*_"$s".npz; do
    [ -e "$f" ] || continue
    cp_tag "$f" "${f/_$s./_$t.}"
  done
}
# plain wire and long zipper: sustained, and the trigger criterion (I_sw identical in results/v2_voltage_trigger/)
reuse results/hitmap conserving_T400_m1.4 kinetic_T400_m1.4
reuse results/v2_voltage_trigger/hitmap conserving_T400_m1.4 kinetic_T400_m1.4_trigger
for f in results/v2_voltage_trigger/hitmap/*_kinetic_T400_m1.4_trigger.*; do mv "$f" results/hitmap/; done
# slot column: both criteria live in results/slots/
reuse results/slots conserving_T400_m1.4 kinetic_T400_m1.4
reuse results/slots conserving_T400_m1.4_trigger kinetic_T400_m1.4_trigger
# narrow strips
for w in 30 56 100; do
  cp_tag results/narrow_baseline/narrow${w}nm_isw_conserving_T400_m1.4.json \
         results/narrow_baseline/narrow${w}nm_isw_kinetic_T400_m1.4.json
  cp_tag results/v2_voltage_trigger/narrow_baseline/narrow${w}nm_isw_conserving_T400_m1.4.json \
         results/narrow_baseline/narrow${w}nm_isw_kinetic_T400_m1.4_trigger.json
done
