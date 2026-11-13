# The device: a double column of 400 nm slots

A 200 nm strip of 9 nm NbTiN with two aligned columns of slots, 40 nm everywhere:

    edge | 40 nm channel | 40 nm slot | 40 nm channel | 40 nm slot | 40 nm channel | edge

- Slots 400 nm long with semicircular ends; 40 nm bridges between consecutive slots (period 440 nm).
- Three 40 nm channels in parallel. The narrowest free cross-section is 120 nm (I_dep = 60.6 µA); at the bridge
  rows the full 200 nm carries the current.
- Simulated device (`geometry.DOUBLE_400`): the central slot pair with one neighbour on each side and 150 nm of plain
  lead to each terminal, 1.58 µm long, 11 squares between the voltage probes.
- References: `Plain_Wire` (200 nm, the same footprint), `Plain_120` (the free width), `Plain_40` (one channel).

The model is in [model.md](model.md), numerics in [methods.md](methods.md). Biases are in I_dep of the narrowest
free cross-section.

## 1 · What the earlier model found (branch `slot-geometries`, 7–9 Nov 2026)

These numbers used the previous photon and heat models (a calibrated ε-hotspot, linearised Joule heat, a current
source). They are kept as the baseline in `results/baseline_kinetic/`; the full record, with the single and
staggered columns, notches and longer slots, is `git show slot-geometries:docs/slots.md`.

| Site | voltage trigger, no heat | sustained, no heat | with heat |
| --- | --- | --- | --- |
| double column, channels | 0.665–0.681 | 0.80–0.83 | 0.78 |
| double column, slot-end zone | | | 0.78–0.83 |
| double column, bridge rows | 0.87–0.93 | 0.92–0.95 | 0.87–0.93 |
| plain 200 nm wire | 0.87–0.90 | 0.87–0.90 | 0.87–0.90 |
| plain 40 nm strip | 0.57 | | |

What that established, and what the new model has to confirm or correct:

1. **The channels work as a cascade**, like parallel nanowires (a SNAP): a photon blocks one 40 nm channel, its
   current moves to the other two, and they switch if they cannot carry it. The voltage trigger sat at the
   three-channel limit 1 − 40/120 = 0.67.
2. **Without heat the cascade healed**, leaving 2–5 flux quanta in the slots. The device stayed resistive only from
   0.80–0.83, so the two criteria were 0.15 apart. With heat the channels detected at 0.78, between the two.
3. **The bridge rows set the plateau.** They are 9 % of the cell length and needed 0.87–0.93, above the plain wire's
   worst site (0.90). With 400 nm slots the plateau test failed (i_95 0.901 against 0.899).
4. **Longer slots were the one lever that helped**: thresholds did not move, the bridge share shrank. Full heated
   maps gave i_95 0.842 (800 nm) and 0.821 (1200 nm). `geometry.double_column(800)` builds them.
5. **Why 400 nm and not 150:** with 150 nm slots the diverted current re-enters the same channel at the next bridge
   and the channels never reach the cascade (0.78–0.82 instead of 0.67).
6. **Method checks that carry over** (methods.md): the short device reproduces a 2 µm strip within 0.011; a 2.0 ξ
   mesh is within 0.019 of 1.4 ξ (channel sites 0.005) and 6–8× faster; dark switching 0.98–1.00.
7. **Limits that carry over:** no thermal noise (the dark switching current is an instability, not a dark-count
   rate; every slot end crowds the current to 1.7× the uniform value); real 9 nm NbTiN strips switch at 0.63 I_dep
   [40]; below threshold the slots keep trapped flux, so the device has memory.

## 2 · With the new model: first look (10 Nov 2026)

Screening runs only: 2.0 ξ mesh, the critical sites, 1550 nm, bath 0.3 T_c, 750 squares in series, a count when
the load carries 20 % of the bias. `results/thresholds/thr_n750_w1550_m2.csv`, `results/events/`,
`results/variants/`. Treat every number as ±0.02 until the 1.4 ξ check and a full map exist.

| Site | new model | previous, with heat | previous, trigger |
| --- | --- | --- | --- |
| plain wire, edge (−96, 0) | **0.664** | 0.866 | 0.872 |
| plain wire, (−58, 0) | **0.727** | 0.877 | 0.877 |
| plain wire, centre (0, 0) | **0.762** | 0.902 | 0.901 |
| double, edge channel (−96 / −80 / −64, 0) | 0.746 / 0.742 / 0.746 | 0.782 | 0.67–0.68 |
| double, centre channel (−16 / 0, 0) | 0.758 / 0.754 | 0.78 | 0.67 |
| double, edge channel beside a bridge row (−80, 220) | 0.750 | 0.870 | 0.869 |
| double, centre channel beside a bridge row (0, 220) | 0.840 | 0.902 | 0.901 |
| double, slot end (−40, 204) | 0.875 | 0.919 | 0.923 |
| double, bridge centre (−40, 220) | 0.906 | 0.933 | 0.933 |
| dark switching: plain wire / double | 0.949 / 0.957 | | 0.997 / 0.981 |

**1. At 1550 nm the slots no longer help.** The energy-calibrated photon is stronger than the old "f = 1": it makes
a normal core about 40 nm across for 20–30 τ0. In the plain 200 nm wire that is enough to start a belt at
0.66–0.76 I_dep, 0.14–0.20 lower than before. The channels of the double column moved by only 0.03 (0.74–0.76),
because they were never limited by the photon: a 40 nm channel is blocked either way, and the threshold is set by
whether the other two channels and the heat carry the cascade through. So the plain wire is now as sensitive as
the channels at its centre and more sensitive towards its edges, and the bridges (0.84–0.91) are an insensitive
tail that the plain wire does not have.

**2. For a lower-energy photon the advantage comes back.** With 4000 nm (0.31 eV) at 0.80 I_dep the double
column's edge channel counts and the plain wire at (−58, 0) does not; both count at 0.90. The channel threshold
barely depends on the photon's energy as long as the photon blocks 40 nm; the plain wire's threshold rises as the
photon weakens. The crossover wavelength is somewhere between 1550 and 4000 nm and has not been located.

**3. What a count looks like** ([event at 0.80](../results/events/Double_400_-80_0_i0.800_n750_w1550_m2.png)).
The hit channel is blocked within 30 τ0 and its current moves to the other two (a small voltage step). About
150 τ0 later the other channels switch, the slotted section between two bridge rows goes normal (3.3 squares,
its geometric resistance), and the circuit takes 65 % of the current out of the detector by 500 τ0. The normal
domain stops at the bridge rows, where the full 200 nm carries the current. At 0.70 the same photon gives 2.5
phase slips and the device heals. A bridge-row photon at 0.90 does nothing visible for 410 τ0 and then starts to
switch.

**4. Area-weighted efficiency is not yet known.** Interpolating nine sites gives i_50 0.78 and i_95 0.87 for the
double column against 0.73 and 0.76 for the plain wire, but the interpolation spreads the bridge-row values along
the whole centre channel (`thr_n750_w1550_m2_map.png`), so the double column's numbers are pessimistic. The sign of
the comparison does not depend on it: every double-column site is at or above the plain wire's centre. A full map
(63 sites) is the run that settles the numbers.

**5. Dark switching in the model is set by the contacts.** With heat, the plain wire's dark switching (0.949) starts
at both normal terminals and runs inward
([event](../results/events/Plain_Wire_dark_i0.955_n750_w1550_m2.png)). The double column is quiet for 800 τ0 at
0.90 and 0.95, because its 200 nm leads carry 0.6 of the plain wire's current density at the same i. Thresholds
above about 0.93 should not be trusted until the contacts are made harmless (wider leads, or a raised T_c next to
them as Zotova & Vodolazov do).

**6. How firm is the comparison?** Events just below and above the thresholds, with one model ingredient changed
at a time ([model.md](model.md) §8), show that the channel threshold is robust to the photon's details (radius,
deposit time) and that the plain wire's is not: with τ_th = 0.6 ps instead of 1.4 ps the plain wire at (−58, 0)
no longer counts at 0.75 while the channel still counts at 0.76. So at 1550 nm the two devices are within the
model's uncertainty of each other in the channels. The bridges are not: they are 0.08–0.15 above the plain wire in
every variant tried.

## 3 · What this means for the design

- The earlier "onset gain of 0.10" came from a photon too weak to open a belt in a 200 nm strip. With hν tied to
  the measured heat capacities, a 1550 nm photon is strong enough, and the gain disappears.
- The slots are a way to make a wide strip respond like a 40 nm wire. That pays when the photon is too weak for
  the wide strip: longer wavelengths, or a film in which the same photon makes a smaller hot spot (larger heat
  capacity, thicker film, higher T_c).
- The lever found earlier still applies: longer slots shrink the bridge share (`Double_800`, `Double_1200`).
- **Before changing the design**, three things decide whether point 1 stands: the 1.4 ξ mesh, the photon's
  thermalisation time for NbTiN, and above all a measured anchor. Real 200 nm NbTiN strips either show a saturated
  1550 nm plateau below their switching current (0.63 I_dep) or they do not; the model says a plain strip needs
  0.66–0.76 I_dep, so it predicts no full plateau at 1550 nm for a strip that switches at 0.63. If your own or
  KTH's data disagree, the absolute photon scale is still off and that is the number to fix first.

