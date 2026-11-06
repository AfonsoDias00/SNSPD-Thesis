# Presentation outline: perforated SNSPD nanowires in TDGL

About 16 slides for a 20–25 minute talk. Each slide gives the message (the one sentence the audience should leave with), what to put on it, the figure, and the references. Reference numbers are those of `docs/report.md` §13. Figures not made by the analysis scripts come from `python presentation_figures.py` (written to `docs/fig/`).

Running thread: *Do holes widen the detection window? No, not by more than the model can resolve. The work shows why, and how far the answer can be trusted.*

---

## Part 1 · Motivation and question (3 slides)

### 1. Title
- Title, name, supervisor, date.
- One device image as background: the zipper threshold map (`results/hitmap/hitmap_conserving_T400_m1.4.png`, right panel only).

### 2. How an SNSPD detects a photon
- **Message:** a photon makes a hotspot; above a threshold bias the strip switches; efficiency saturates at a "plateau".
- Sketch: strip, bias current, hotspot, current diverted around it, resistive belt or vortex crossing.
- A schematic efficiency curve with i_det (onset), the plateau and I_sw marked, and the detection window 1 − i_95.
- **Figure:** still to draw by hand: a schematic. Or reuse one efficiency curve (`results/hitmap/efficiency_conserving_T400_m1.4.png`, plain wire only) with the window shaded.
- **Refs:** Gol'tsman 2001 [7] (first SNSPD); Semenov 2001 [6] (hotspot picture); Engel 2015 [14] (detection mechanism review).

### 3. The question: can internal holes widen the window?
- **Message:** narrow strips saturate at lower bias. Holes split a wide strip into narrow channels, so maybe they give the narrow-strip gain without losing area.
- Experimental motivation: narrower wires saturate earlier (Marsili 2011 [16], Lusche 2014 [17]). Parallel or cascade designs exist (SNAP: Ejrnaes 2007 [18], Marsili [19]).
- Current crowding at holes and bends (Clem & Berggren 2011 [11], Hortensius 2012 [12]): the possible downside.
- State the success test: **i_95 lower than the plain wire's by more than 0.03 I_dep** (the pre-registered go / no-go).

---

## Part 2 · Method (3 slides)

### 4. Model: TDGL in pyTDGL
- **Message:** generalized TDGL on a 2D mesh; the photon enters only through ε(r, t).
- One equation line (GTDGL), the parameter table in brief: ξ = 5 nm, λ = 600 nm, d = 9 nm (NbTiN), w ≪ Λ so no screening, I_dep = j_dep × narrowest width.
- Units: thresholds in I_dep, so λ drops out of every normalised result.
- **Refs:** Kramer & Watts-Tobin 1978 [2]; pyTDGL, Bishop-Van Horn 2023 [1]; Jönsson 2022 [13] (NbTiN λ fit); Tinkham [4] (depairing current); Zotova & Vodolazov 2012 [8] (TDGL for photon detection).

### 5. Devices
- **Message:** plain 200 nm wire vs long zipper vs fabrication-feasible slot columns, with narrow plain strips as the control.
- **Figure:** `docs/fig/devices.png`: every device to one scale, with width, narrowest free width, I_dep and length.
- Table from report §1 (devices).

### 6. How a threshold is measured
- **Message:** bisect the bias at each photon site → threshold map → efficiency curve → i_50, i_95.
- Flow: device → mesh → dark I_sw → photon → 400 τ0 run → detected? → bisect → map → curve (report §0 table, as a diagram).
- **Two detection criteria = two bounds:** the voltage trigger (lower bound) and sustained switching (upper bound), because there is no Joule heating (report §4.4).
- **Figure:** `results/hitmap/hitmap_conserving_T400_m1.4.png` (map) next to `results/hitmap/efficiency_conserving_T400_m1.4.png` (curve).
- **Refs:** Zotova & Vodolazov 2014 [9] (efficiency from local thresholds); Vodolazov 2017 [10] (resistive-state criteria).

---

## Part 3 · Results (6 slides)

### 7. Main result: holes do not lower the plateau
- **Message:** the zipper turns on earlier but saturates no earlier than the plain wire.
- Numbers (conserving, 400 τ0): i_50 0.805 vs 0.845; i_95 0.853 vs 0.850; i_max 0.907 vs 0.850 I_dep. Go / no-go 0/3.
- Point at the map: channels and edges are more sensitive; bridges and hole tips are much less sensitive and form the tail.
- **Figures:** `results/hitmap/efficiency_conserving_T400_m1.4.png` (main) and `results/hitmap/hitmap_conserving_T400_m1.4.png`.

### 8. Narrow-strip control and the two bounds
- **Message:** with the voltage trigger a 56 nm strip beats the zipper by about 0.2 I_dep. With sustained switching all narrow strips collapse to a floor, because without heat a narrow strip heals. Adding heat settles it for the narrow strips (slide 14b).
- Width-scan numbers: i_50 0.33 / 0.60 / 0.77 / 0.85 at 30 / 56 / 100 / 200 nm (lower bound) vs ~0.82 for all (upper bound).
- The 30 nm example: 65 phase slips (~10 eV) and then it recovers.
- **Figure:** `docs/fig/width_scan.png`: i_50 and i_95 vs width, both bounds × both hotspots, zipper as stars at its 173 nm free width.
- **Refs:** Marsili 2011 [16], Lusche 2014 [17] (experiment favours narrow strips → the lower bound is closer); Zotova & Vodolazov [8], Vodolazov [10] (electro-thermal retrapping).

### 8b. What a photon does: the cascade (kinetic hotspot)
- **Message:** the hit channel loses its current; the neighbours take it over and are pushed past their limit in turn. That happens in the plain wire too, but holes split it into separate channels that switch one after another.
- What to point at:
  - **Zipper:** the hit channel drops from 23 to 11 µA. The middle channel overloads to 31 µA and collapses about 65 τ0 after the photon, then the belts cross the right column.
  - **Slot 400:** the hit channel collapses within about 40 τ0. The other channel takes 43 µA and the slip belts appear at the slot ends.
  - **Plain wire:** the belt forms only about 90 τ0 after the photon.
- **Figures:** `docs/fig/cascade_kinetic.png` (channel supercurrents and probe voltage, 0.01 I_dep above and below each site's threshold); `docs/fig/snapshots_kinetic_Long_Zipper.png` and `docs/fig/snapshots_kinetic_Slot_400.png` (|ψ|² with the hotspot contours); `docs/fig/snapshots_kinetic_Plain_Wire.png` for reference. Old version for contrast: `docs/fig/cascade_traces.png` (21 Oct, legacy hotspot).
- **Caveat to say aloud:** Slot 400 at 0.818 I_dep is resistive for all 400 τ0 but heals in the 150 τ0 continuation. At 0.8180 instead of 0.8186 I_dep it detects, so its resistive state near threshold is marginal (report §9b).
- **Refs:** SNAP / cascade switching [18, 19]; Zotova & Vodolazov [8] (belt formation in TDGL).

### 9. Where the zipper's early turn-on comes from: trapped flux
- **Message:** below threshold a photon leaves 1–2 flux quanta in a hole with no voltage. What happens next depends on the hotspot. With the conserving one, the zipper switches through that flux 80–250 τ0 later. With the kinetic one, there is less flux, and the delayed switches become pulses that heal (report §9d).
- Bullets:
  - silent flux at 0.70 I_dep: 59 % of the unit cell (conserving), 25 % (kinetic); at 0.82: 30 % vs 6 %;
  - switches via trapped flux: 68 runs vs 11; pulses that heal: 30 vs 76, median 3 vs 15 phase slips;
  - counting only direct switches, the onset gain is 0.016 vs 0.015 I_dep in both models.
- **Figures:** `results/fluxtrap/outcomes_Long_Zipper_T400_m1.4.png` (main: share of the unit cell by outcome vs bias, conserving and kinetic side by side, plain wire as a line). For one site: `results/fluxtrap/trace_Long_Zipper_-68_0_conserving_T400_m1.4.png` and `..._kinetic_T400_m1.4.png` (flux held and phase slips vs time at four biases).
- **Refs:** Mkrtchyan & Shmidt 1972 [32], Buzdin 1993 [33] (hole saturation number); Bulaevskii 2012 [20] (vortex-assisted detection).

### 10. Trapped flux gives the device a memory (optional, or merge with slide 9)
- **Message:** a trapped state lowers the zipper's own I_sw from 0.968 to 0.79–0.91 I_dep and shifts the next photon's threshold by −0.03 to +0.06.
- Table: trapped state → own I_sw (report §9.3).
- An experimental test to propose: the I_sw distribution after dark cooling vs after sub-threshold light vs after warming.
- **Refs:** Casaburi 2015 [35] (photo-induced vortex crossing); Onen 2019 [36] (single-flux memory in a loop); Kirtley 2001 [38] (fluxoid dynamics in rings).

### 11. A fabrication-feasible layout: slot columns
- **Message:** 30 nm slots with 40 nm bridges in a 2 µm strip; no slot length lowers the plateau. The tail sits in the bridges.
- Numbers: Slot 400 i_95 0.849–0.860 vs 0.851 for the plain 2 µm wire (conserving); every column is no-go.
- **Figures:** `results/slots/efficiency_conserving_T400_m1.4.png` and `results/slots/hitmap_conserving_T400_m1.4.png`.

### 12. A better hotspot: Simon et al. kinetic model
- **Message:** with a delayed build-up and loss to the substrate, the photon no longer makes a normal core. It becomes a shallow, longer-lived suppression. The verdict does not change.
- Left: what the paper gives (α mapping, τ_th = 1.4–1.8 ps, τ_esc = 15 ps) vs what we chose (walls, amplitude, ε floor), from report §3.3.
- Right: the threshold changes. Plain i_50 0.845 → 0.833; zipper i_95 0.853 → 0.835 vs plain 0.848; narrow-strip floor 0.82 → 0.775–0.818; slots shift with the wire.
- **Figures:** `results/hotspot_profiles_conserving_vs_kinetic.png` (main); efficiency curves with the kinetic hotspot beside the conserving ones: `results/hitmap/efficiency_kinetic_T400_m1.4.png` (zipper vs plain, sustained; `_trigger` for the lower bound) and `results/slots/efficiency_kinetic_T400_m1.4.png` (slots vs 2 µm plain, same naming).
- **Refs:** Simon et al. 2025 [27] (PRB 112, 174512); Zotova & Vodolazov [8] (τ_e-ph scale); Masaki & Matsueda 2025 [28] (sensitivity to hotspot size).

---

## Part 4 · How far to trust it (3 slides)

### 13. Validation and checks
- **Message:** the pieces reproduce known answers; the weak spot is anything involving flux in narrow bridges.
- Three groups:
  1. The hotspot solver matches exact diffusion solutions (second order in time and space).
  2. The ε mapping is used as given by pyTDGL.
  3. Switching happens at depairing for every geometry (I_sw/I_dep = 0.968–0.971).
- Then the fragile points: mesh 1.4 vs 0.8 ξ moves the [0 −2] state by 0.06; the zipper re-meshed on another machine moves flux-mediated sites by 0.04.
- **Figure:** `docs/fig/validation_grid.png`: every check by group, with failures classed as wrong expectation / meter limit / real issue.
- **Refs:** Masaki & Matsueda [28] (mesh ≤ ξ/2 at rims).

### 14. Sensitivity: the hotspot decides differences of 0.03
- **Message:** the choices Simon et al. leave open move thresholds by 0.01–0.03 I_dep. That is the size of the go margin and of every hole-vs-plain difference.
- What is robust: no normal core forms; narrow strips get easier under the upper bound; edge hits get easier.
- What is not: the zipper's plateau edge; its flux-mediated sites flip between 0.788 and 0.817.
- **Figure:** `docs/fig/sensitivity.png`: change from the conserving threshold per site for the four kinetic variants, with the spread and the 0.03 go-margin band.
- **Ref:** Simon et al. [27].

### 14b. Adding heat merges the two bounds
- **Message:** a heat equation coupled to the TDGL (Joule heating, diffusion, escape; one- and two-temperature) puts every threshold at or below the lower bound. Narrow strips then beat every zipper site. The zipper's tail stays above the plain wire.
- Numbers (report §9e): 30 / 56 / 100 nm at 0.44–0.46 / 0.61 / 0.72–0.76 I_dep; zipper onset 0.760–0.768 vs plain 0.801–0.827; zipper tail 0.870 vs plain ≤ 0.827. Halving or doubling the heating moves thresholds by ≤ 0.03 (0.06 at 30 nm).
- **Figure:** `results/thermal/check.png` (each site's bracket between the bounds, with the threshold with heat marked).
- **Caveat to say aloud:** a prototype at 14 sites, with estimated heat parameters and no readout circuit.
- **Refs:** Zotova & Vodolazov 2012 [8] (the heat equation); Sidorova 2021 / 2022 [39, 40] (measured NbTiN parameters).

### 15. What the model leaves out
- No readout circuit, and heat only as a prototype → latch vs reset is not decided (Zotova & Vodolazov [8], Vodolazov [10]).
- No thermal noise → no dark-count rates. Vodolazov 2012 [21] gives the GL barriers; in an ideal straight film they reach a few tens of k_BT only at ≈ 0.96 I_dep.
- Ideal strips: the seven NbN/WSi devices of Frasca 2019 [23] switch at 0.43–0.74 I_dep, and 9 nm NbTiN strips from the same KTH process at 0.63 I_dep [40], limited by defects (Zhang 2014 [31] on constrictions). Every f = 1 threshold of the 200 nm devices is above 0.63, so the absolute photon scale is not calibrated.
- Measured times for this film class [39, 40]: τ0 ≈ 0.15 ps is confirmed, but energy leaves the film in 80–90 ps (530 τ0, not 100). With those times every threshold drops by 0.01–0.10 I_dep, and the order of the geometries stays (report §9c.4).
- GL is used well below T_c: trust ratios in I_dep, not absolute times or currents (Kramer & Watts-Tobin [2], Watts-Tobin 1981 [3]).

---

## Part 5 · Conclusion (1–2 slides)

### 16. Conclusions and next steps
- **Conclusions:**
  1. Holes do not widen the window by the 0.03 I_dep margin, under either bound or hotspot, for the zipper or for any slot length.
  2. Their earlier turn-on comes with trapped flux, which gives the device a memory. With heat it is about 0.05 I_dep at the sites checked; it is not a wider window.
  3. With heat, narrow plain strips beat every zipper site.
- **Next steps:**
  1. Full threshold maps with heat (the prototype merges the two bounds at 14 sites) [8, 10].
  2. Calibrate the heat and the photon on this film: retrapping current, a two-temperature photon, R(T) of our own film; τ_th is still unmeasured [39, 40].
  3. A readout circuit (kinetic inductance + 50 Ω), to tell a latch from a pulse that resets.
  4. Defect-limited baseline [23, 40] and barrier-based dark counts [21], once the absolute photon scale is fixed.
- Optional closing image: the zipper's threshold map, kinetic (`results/hitmap/hitmap_kinetic_T400_m1.4.png`).

---

## Backup slides
- **Method history** (the criterion fixed three times, report §4.2, §5): `docs/fig/plain_vs_holes.png` vs `docs/fig/slender_small.png` show that the early ranking flipped with the hotspot size. Good for "why one fixed hotspot".
- **Column cascade** (21 Oct): `docs/fig/cascade_traces.png`, with the SNAP refs [18, 19].
- **Current crowding:** `results/hitmap/crowding_0.95_conserving_T400_m1.4.png` and `results/slots/crowding_0.95_conserving_T400_m1.4.png` [11, 12, 13].
- **Selectivity vs photon energy:** `results/hitmap/selectivity_conserving_T400_m1.4.png`; weak photons are detected at the edge first [14, 15, 20].
- **Slot 400 flux scan** (pulses that end quiet, both hotspots): `results/fluxtrap/outcomes_Slot_400_T400_m1.4.png`.
- **Measured NbTiN times:** report §3.2 and §9c.4 (thresholds drop by 0.01–0.10 with the measured escape time) [39, 40].
- **Edge banks (parked):** one line plus refs [24, 25, 26] if asked about edge designs.
- **Full numbers table:** report §9c.1 (conserving → kinetic for every device).

## Practical notes
- Use I_dep units everywhere, and say once why (λ-independent; I_sw drifts with the window).
- Always say which bound a number is. A safe default for talks is "upper bound, lower bound in brackets".
- The [R]-marked references (15, 19, 22, 24–31) were not re-checked against the publisher; verify the ones you cite (report §13). References 8, 10, 21, 23 and 27 were checked on 4–5 Nov; 10, 21 and 23 had wrong titles or details and are corrected in the report.
