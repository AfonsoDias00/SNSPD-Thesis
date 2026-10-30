# Results with the voltage-trigger criterion (28-29 Oct 2026)

The results here used the criterion "the probe voltage peaks above 0.45 V0 after t0 + 4 tau0, or its 25 tau0 trailing
mean exceeds 0.18 V0". It is the first criterion (`../v1_voltage_or_collapse/`) without the channel-collapse branch.

On 29 Oct 2026 the criterion changed to **sustained switching** (`snspd/sim.py`). The voltage trigger must fire, and
there must also be at least 1 phase slip between the probes in the last 50 tau0 of the window. The same must hold again
after the run is continued for another 150 tau0. The reason is `fluxtrap.py` (results in `../fluxtrap/`):
- In the long zipper, a photon can push trapped flux from one hole to another.
- That costs 1-6 phase slips, enough to fire the trailing-mean branch, and then the device goes quiet again with new
  trapped flux.
- At 0.78 I_dep this happens over 24 % of the zipper's unit cell, and the old criterion counted all of it as detection.
- The plain wire never does this: its thresholds are the same under both criteria.

Contents:
- `hitmap/`: plain wire vs long zipper at 400 tau0 and 1.4 xi (full grid), plus the 150 tau0 critical-grid check.
  Zipper i_50 0.795 / i_95 0.853 I_dep; plain wire 0.845 / 0.850.
- `narrow_baseline/`: 30 / 56 / 100 nm strips at 400 tau0 (30 nm also at 0.8 xi).
- `checks/`: the first fluxoid trace at (-68, 0); superseded by `../fluxtrap/trace_*`.
- `fluxtrap/`: the outcome scans of the zipper and the plain wire. Their transient / sustained split used the
  in-window test only, without the continuation.
- `fluxtrap/` also holds the trapped states' switching currents and second-photon thresholds, bisected with this
  criterion.
  - The [0 -2 0 0 0 0 0] state "switches" at 0.775 I_dep (0.795 at 0.8 xi), but its flux only moves to another hole.
  - The other states switch for real at 0.87-0.91.

Reports regenerate with the current code, e.g.
`python hitmap.py report --tmax 400 --mesh-xi 1.4 --out results/v2_voltage_trigger/hitmap`.
