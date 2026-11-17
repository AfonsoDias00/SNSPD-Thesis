# The model: photon, heat and circuit

Branch `double400-models`, 10 Nov 2026. This note explains what the simulation now solves, where each piece comes
from, what was checked, and what is still an assumption. Numerics and protocol (mesh, time step, settled start) are
in [methods.md](methods.md); the device and the earlier results are in [double400.md](double400.md); the papers
are listed in [references.md](references.md).

## 1 · What changed, in one table

| | Before (branch `slot-geometries`) | Now |
| --- | --- | --- |
| Photon | a suppression ε(r, t) of the TDGL coefficient with amplitude \|A\| = 2.5, cap s_max = 1.25, radius 14 nm, all chosen by us; "f = 1" was not tied to a photon energy | the photon's energy hν (0.80 eV at 1550 nm) put into the electrons; what it does to ψ follows from the measured heat capacities |
| Heat | a second field θ added to the photon's ε: constant heat capacities taken at T_c, linear ε = 1 − θ, Joule heat only, no condensation energy | one field for everything: electron and phonon temperatures with c_e ∝ T, c_ph ∝ T³ and the measured T^5.4 exchange; Joule **and** order-parameter relaxation heat; condensation energy (latent heat); reduced heat conduction in the superconductor |
| ε(T) | linear, (T_c − T)/(T_c − T0) | two-fluid, (1 − t²)/(1 − t0²): gives the depairing current's measured temperature law |
| Circuit | none (current source) | series kinetic inductance and 50 Ω load |
| Detection | two criteria that were bounds (voltage trigger; sustained), 0.15 I_dep apart in the double column | one: the load carries ≥ 20 % of the bias current |
| Free parameters of ours | \|A\|, s_max, r0, the wall treatment, η (factor 2) | none in the heat model's coupling; see §7 for what is still assumed |

## 2 · What the papers do

Read this session from the source (how much of each: [references.md](references.md)).

| Paper | Photon | Heat | Circuit |
| --- | --- | --- | --- |
| Simon et al. 2025 [27] | ab initio quasiparticle cascade → ε(t) = ε_∞(1 − e^{−t/τ_th}) in α = (1 − T/T_c − ε)/(1 − T/T_c); initial volume πξ_c²d, radius grows as √(Dt). **ε_∞, profile, dilution and edge treatment are not given** | none; ~25 % of hν lost to the substrate in 5 nm NbN (τ_esc 10 ps) | none. Outlook asks for electron and phonon temperature equations and a circuit |
| Zotova & Vodolazov 2012 [8] | energy in the heat equation: hν = ΔT · πR²d · C_v, R = 9 nm (Eq. 1) | one temperature, constant C_v at T_c, Joule heat ρ_n j_n² only, phonons at the bath (Eq. 2) | L_k dI_s/dt = (I − I_s) R_shunt − V_s (Eq. 5); without the shunt the film latches |
| Vodolazov 2017 | hot spot 2ξ_c × 2ξ_c with T_e = T_ph = T_init from energy conservation, **including the condensation energy** | electrons and phonons, T_e⁵ − T_ph⁵ exchange, k_s(T_e, Δ), electron energy with the condensate term (Eqs. 30–32); TDGL with modified coefficients for low T (Eq. 36) | none (fixed bias) |
| Yasukawa et al. 2024 | — (an nTron) | one temperature, constant C; **Joule heat plus order-parameter relaxation heat** (Eq. 3) | L_k dI/dt + V = R_L (I_bias − I) (Eq. 4) |
| Sidorova et al. 2021, 2022 [39, 40] | — | the measured two-temperature model of our film: c_e ∝ T, c_ph ∝ T³, τ_EP ∝ T^−3.4, τ_esc (Eq. 6 and Table 1 of [40]) | — |
| Kerman et al. 2009 | — | lumped hot-spot growth | same loop equation; latching when the electrical time L/R_L is too short against the thermal time |
| Annunziata et al. 2010 | hν into one grid cell of the electron temperature (Eq. 8) | two temperatures with T-dependent capacities, GL I_c(T) | same loop equation (Eq. 3); ½L_kI_b² is dissipated in the hot spot |

No paper found couples all three to a TDGL for a single-photon detector. Every piece used below is in at least one
of them.

## 3 · The photon

### How the previous model differed from Simon et al.

Simon et al. compute the suppression from microscopic kinetics and feed pyTDGL a coefficient α(r, t). What the paper
states and what it leaves open:

- **Stated:** the definition of ε (their Eq. 3, a weighted integral of the excess quasiparticle distribution), the
  link α = (1 − T/T_c − ε)/(1 − T/T_c) (Eq. 10), the build-up time τ_th = 1.4–1.8 ps (NbN, 1064 nm), the initial
  volume πξ_c²d, growth as √(Dt), τ_esc, γ = 10, strips of 20–85 nm, bias 0.5 I_dep, runs of 5–25 ps.
- **Not stated:** any value of ε_∞, the spatial profile, how the suppression dilutes as it spreads, what happens at
  edges, and what limits α from below.

The previous model took the build-up time from the paper and filled the gaps with its own choices: a Gaussian of
14 nm, an amplitude calibrated to an older model so that "f = 1 is the same photon", a tanh cap, and diffusion with
reflecting walls. The transport choice was sound (the paper itself says a better treatment of diffusion would
"significantly improve the accuracy for wires with larger w"). The amplitude was the weak point: nothing tied it
to a photon energy, which is why [methods.md](methods.md) had to say "the absolute photon scale is not calibrated".

### The decision: deposit the energy, not the suppression

For a quasiparticle distribution that has thermalised to a temperature T_e, Simon's ε reduces to the ordinary
temperature dependence of the GL coefficient. With the BCS spectral function R/|Δ| = 1/√(E² − Δ²) and Δ → 0,

    ε = 2 ∫ dE [f(E, T_e) − f(E, T)] / E = ln(T_e / T),      so      α = ln(T_c/T_e) / ln(T_c/T) ≈ (T_c − T_e)/(T_c − T).

(This step is my derivation from their Eq. 3 with the standard BCS integral; the paper does not write it.) So an
electron temperature is the thermalised limit of their model, and the missing amplitude becomes a question with a
measured answer: how hot do hν = 0.80 eV make the electrons and phonons of this film? That is what Zotova &
Vodolazov, Vodolazov 2017 and Annunziata do, and it is what the code does now:

- **Energy:** hν × (absorbed fraction, 1 by default) goes into the electron system.
- **Where:** a Gaussian of 1/e radius r0 = 5 nm ≈ ξ (Simon: πξ_c²d; Vodolazov: 2ξ_c × 2ξ_c). The discrete sum is
  normalised, so the energy is exact on any mesh.
- **When:** at the rate e^{−(t − t_a)/τ_th}/τ_th with τ_th = 1.4 ps (10.6 τ0), Simon's build-up time. The deposited
  fraction is then 1 − e^{−t/τ_th}, the same curve as their ε(t).
- **Afterwards:** nothing special. The energy diffuses, passes to the phonons and escapes by the heat model of §4.

What this predicts for one 1550 nm photon in 9 nm NbTiN (heat model alone, no current): the electrons reach 2 T_c
in the core within a few τ0; within ~30 τ0 the phonons take three quarters of the energy and a region of about
1200 nm² (radius 20 nm) is above T_c for 20–30 τ0; after 80 τ0 the hottest electrons are at 0.76 T_c (ε = 0.47)
while the phonons under them stay near T_c and keep feeding them for the escape time (607 τ0). A quick count gives
the same scale: bringing film from the bath to T_c costs 0.91 (electrons) + 3.11 (phonons) + about 0.5
(condensate) energy units per ξ², and the photon is 101.7 units, enough for 22 ξ², a disc of 13 nm radius. So a
telecom photon makes a short-lived normal core of about the old model's size and then a long-lived warm spot. The
old kinetic hotspot never made a normal core at f = 1 (ε ≥ +0.43).

**Caveat from the same papers.** During the first τ_th the distribution is not thermal: Simon et al. warn that
"caution must be exercised when assuming a well-defined electron and phonon temperature", and Vodolazov 2017
finds that NbN is not thermalised by the time the spot has grown to ξ. Spreading the deposit over τ_th is the
mitigation; it is not a kinetic calculation.

## 4 · Heat

### Equations

State on every mesh site, with temperatures in units of T_c (t = T/T_c) and n = |ψ|²:

    U = t_e²      electron thermal energy in units of γT_c²/2        (c_e = γT)
    P = t_ph⁴     phonon energy in units of c_ph(T_c) T_c/4          (c_ph ∝ T³)

In pyTDGL units (lengths ξ, time τ0):

    T_e ∂s_e/∂t = ∇·(k_s ∇T_e) + W − P_ep + photon,        s_e = γT_e (1 − b n)                        (1)
    ∂P/∂t       = 4r X − (P − P0)/τ_esc                                                                (2)
    X           = (t_e^q − t_ph^q)/(q τ_ep),   q = 2 + 3.4 = 5.4                                       (3)
    W           = η0 [ |∇μ|² + u/√(1+γ²n) · ( |(∂_t + iμ)ψ|² + (γ²/4)(∂_t n)² ) ]                       (4)
    ε           = (1 − U)/(1 − U0)                                                                     (5)

| Piece | Source | Why it is there |
| --- | --- | --- |
| (1)–(3): two temperatures, c_e ∝ T, c_ph ∝ T³, exchange ∝ T_e^q − T_ph^q, escape | Sidorova 2022 Eq. (6), with q = 2 + n and s = 4. For q = 5 it is Vodolazov 2017 Eqs. (30)–(31) with the same prefactors (checked by algebra: 4r/(5τ_ep) on both sides) | the capacities at the bath are 0.3× (electrons) and 0.03× (phonons) of their T_c values. Holding them at T_c, as before, overstated the energy needed to reach T_c by 1.5× and 2.8× |
| electron conduction k = c_e D | Wiedemann–Franz; in U it is plain diffusion, ∂U/∂t = D∇²U | — |
| k_s/k_n = 1 − (6/π²)∫₀^{Δ/k_BT} x²eˣ/(eˣ+1)² dx | Vodolazov 2017 Eq. (32) | cold superconductor conducts 9 % of the normal value: heat is held in the hot spot |
| Joule term \|∇μ\|² | Zotova & Vodolazov Eq. (2) | — |
| relaxation term in (4) | Yasukawa 2024 Eq. (3) (after Schmid 1966), extended to the generalised TDGL that pyTDGL solves | **it is 28–30 % of the dissipation** in a strip with phase slips (validation 3b). The old Joule-only source lost it |
| s_e = γT(1 − bn): condensate entropy | the entropy −∂f/∂T of the same free energy whose TDGL equation is solved; the role of the term E₀𝓔_s(T_e, \|Δ\|) in Vodolazov 2017 Eq. (30) | destroying superconductivity costs energy. The photon has to pay it, and a recovering condensate gives it back |
| ε(T) = (1 − t²)/(1 − t0²) | Bardeen's (1 − t²)^{3/2} depairing law, which is the temperature factor of the Kupriyanov–Lukichev formula that Sidorova 2022 (Eq. 5) use for these very strips and of Vodolazov 2017 Eq. (29) | pyTDGL has one knob, ε, and j_dep ∝ ε^{3/2}. Linear GL gives (1 − t)^{3/2}: heating from 0.3 to 0.5 T_c would cut the local critical current to 0.60 instead of 0.75 |
| bath temperature at the contacts | Zotova & Vodolazov; Vodolazov 2017 | pyTDGL's normal terminals dissipate steadily; a real strip between superconducting leads does not |

### Why W is complete

For pyTDGL's equations at fixed ε, the free energy F = Σ A(−εn + n²/2) + Σ_edges w|Δψ|² obeys

    dF/dt = 2 I V − 2 Σ A ( |∇μ|² + W_relax ),

with I and V the bias current and terminal voltage in pyTDGL units. The derivation is three lines (multiply the TDGL
by ∂_tψ*, take the real part, use the imaginary part of the equation for ∇·J_s and integrate by parts). So every
joule the source delivers either stays in F or appears as one of the two terms in (4). `validation 3b` checks this
identity on a strip driven through phase slips: the balance closes to 0.1 %.

### Why the two-fluid ε, and what it buys thermodynamically

The free energy density behind the TDGL is f = f_n(T) + 2E_c[−ε(T)n + n²/2 + ξ²|∇ψ|²]. Its entropy is
s = s_n + 2E_c ε′(T) n, and the heat equation is T ∂s/∂t = (heat in). With the two-fluid ε this entropy is
γT(1 − bn), b = η0/(1 − t0²): proportional to T, as the third law wants, and the free energy difference is
E_c(0)(1 − t²)², the parabolic critical-field law. With the linear GL ε the condensate's entropy deficit is a
constant; at 0.3 T_c it is 795 J/(K m³) against a normal-state electron entropy of 391, so the superconducting state
would have negative entropy. That is the concrete form of "GL does not hold far below T_c", and it is why a latent
heat cannot be added to the linear model. The code refuses that combination.

What the two-fluid functional implies (validation 3a.9–11): entropy of the superconductor at the bath 0.14 γT (BCS
about 0.2), specific-heat jump 1.89 γT_c (BCS 1.43, Gorter–Casimir 2). It is the Gorter–Casimir two-fluid model
with the BCS condensation energy: right at both ends, about 30 % high in the jump.

In the variable U, (1) becomes (1 − bn) ∂U/∂t = 2(η0 w + …) + 2bU ∂n/∂t. The last term is the latent heat.

### The coupling constant η0 has no free parameter

η0 = 4E_c/(γT_c²) converts the TDGL's dissipation into electron energy. The old code computed it from λ and ξ
(E_c = Φ0²/(16π²μ0λ²ξ²)) and gave it a factor-2 uncertainty. But the same number fixes the condensate entropy, and
that must stay positive: 4E_c(0)/(γT_c²) ≤ 1. BCS gives (3/π²)(Δ0/k_BT_c)² = 0.946. So the code takes

    η0 = 0.946 (1 − t0²)² = 0.783     (t0 = 0.3)

from thermodynamics. λ = 600 nm, ξ = 5 nm and the measured c_e give 0.855, that is 1.03 instead of 0.946: the two
routes agree to 9 %, well inside the uncertainty of λ²ξ². It is equivalent to λ = 629 nm, between the hairpin fit
(600 nm) and the dirty-limit BCS value from the measured R_sq and T_c (630–660 nm). The heat model therefore no
longer depends on λ or ξ except as units.

### Numbers (9 nm NbTiN strip, bath 0.3 T_c = 2.6 K)

| Quantity | Value | Source |
| --- | --- | --- |
| T_c | 8.60 K (200 nm strip) | [40] Table 1 |
| R_sq | 381.4 Ω → τ0 = L_sq/R_sq = 0.132 ps | [39, 40]; with λ = 600 nm [13] |
| D | 0.472 cm²/s = 6.22 nm²/τ0 | [39] |
| τ_EP(T_c), n | 4.0 ps = 30.4 τ0, n = 3.4 (energy relaxation time) | [40] Table 1 |
| τ_esc | 80 ps = 607 τ0 | [39, 40] |
| c_e(T_c), c_ph(T_c) | 1303, 8173 J/(K m³); r = 0.159 | [40] Table 1 |
| Δ0/k_BT_c | 1.764 (BCS) | **assumed**; [40] also uses 1.76 and says the NbTiN value is not known |
| energy unit (γT_c²/2) d ξ² | 7.87 meV per ξ² | from c_e, T_c, d |
| photon, 1550 nm | 0.800 eV = 101.7 units | — |
| τ_th | 1.4 ps = 10.6 τ0 | Simon et al., for NbN. **Not measured for NbTiN** |

τ0 was 0.15 ps in the old heat preset (from ħ/πΔ0). It is now fixed by the measured sheet resistance and λ, so D and
the times in τ0 moved by 12 %.

## 5 · Circuit

    L_k dI/dt + V = R_L (I_b − I)

A current source feeds the detector (the simulated strip in series with the rest of the nanowire, a kinetic
inductance L_k) in parallel with the load. Zotova & Vodolazov Eq. (5), Yasukawa Eq. (4), Kerman's loop equation
and Annunziata Eq. (3) are this equation.

In pyTDGL's units two identities hold: V_u/I_u = R_sq and τ0 = L_sq/R_sq (L_sq = μ0λ²/d, the sheet kinetic
inductance). The equation becomes

    dI/dt = [ (R_L/R_sq)(I_b − I) − I_u v ] / N_sq,        N_sq = L_k / L_sq,

with v pyTDGL's voltage between the probes. **The only parameter is N_sq, the number of squares of nanowire in
series**; R_L/R_sq = 50/381.4 is measured. Checked in `validation 3c`: a strip held normal shows exactly its
geometric number of squares, and divides the current with the load as k/(k + n).

- Default N_sq = 750: a 200 nm × 150 µm strip, the devices of [40]. L_k = 38 nH, L_k/R_L = 0.75 ns = 5700 τ0.
  **This is a property of your device, not of the model. Set it to your wire length** (`--nsq`).
- A normal region of n squares empties the detector in N_sq/n τ0 (250 τ0 for three squares). One phase slip moves
  Φ0/L_k = 0.055 µA.
- N_sq → ∞ is the current source (`--source`).

### The detection criterion that follows

A count is a pulse: **the load carries at least 20 % of the bias current at some time.** A burst of phase slips that
heals moves 0.1 % of the current per slip and cannot reach that; a resistive domain can. This replaces the two
bounds of the earlier work. An event that is resistive but below the level at 400 τ0 runs on (to 800 τ0) until it
counts or goes quiet. With `--source` the criterion is the old heated one (still resistive at the end of the
window).

## 6 · Implementation

| File | Content |
| --- | --- |
| `snspd/units.py` | the film, unit identities |
| `snspd/thermal.py` | `Material`, `thermal_params` (measured values → pyTDGL units), `Photon`, `HeatField` |
| `snspd/circuit.py` | `Circuit`, `CircuitState` |
| `snspd/sim.py` | `march` (time stepping), `settled_state`, `ETSolver`, `run_event`, `threshold` |
| `event.py`, `thresholds.py` | one event with traces and snapshots; threshold bisection over sites |

**Coupling.** `ETSolver` subclasses `tdgl.TDGLSolver`. After each solver step it hands the fields before and after
the step to `HeatField.step`, which returns ε for the next step, and the probe voltage to `CircuitState.step`,
which returns the terminal current for the next step. Both couplings lag by one step of ≤ 0.025 τ0.

**Inside `HeatField.step`** (finite volumes on pyTDGL's own Voronoi cells, so sums of energy are exact):
1. latent heat: the entropy γ t (1 − bn) is conserved while n changes, so U is rescaled by (g_old/g_new)²;
2. sources: dissipation (4) and the photon, as energy added at fixed n;
3. conduction: explicit, with the harmonic mean of k_s on each edge, sub-stepped if the step would be unstable;
4. exchange and escape: linearly implicit in the exchange power X. It is unconditionally stable for the stiff T^5.4
   term and moves exactly the same energy out of the electrons as into the phonons;
5. contacts: sites within 50 nm of a terminal are reset to the bath. (25 nm, the old value, left the electrons at
   0.41 T_c just outside it in the dark, because at γ = 10 the current-conversion zone is about 40 nm deep.)

**`sim.march`** calls `TDGLSolver.update` directly instead of `tdgl.solve`: the same adaptive steps, but no output
file, the final state returned as arrays, and a run that can stop as soon as its outcome is known.

## 7 · What is still assumed

Ranked by how much I expect each to matter.

1. **GTDGL at 0.3 T_c.** The two-fluid ε fixes the temperature law of the depairing current, not the dynamics.
   Vodolazov 2017 changes two more coefficients (ξ_mod, Δ_mod) and adds a Usadel correction to the supercurrent;
   pyTDGL cannot vary those in space.
2. **τ_th and r0 of the deposit.** τ_th is an NbN number; Vodolazov estimates 0.57 ps for NbN. Both are options
   (`--photon tau_th_ps=0.6,r0_nm=10`).
3. **Δ0/k_BT_c = 1.764.** Strong coupling (2.0–2.1 is quoted for NbN) would need a different thermodynamic closure:
   with the two-fluid functional the entropy at a 0.3 T_c bath goes negative above about 1.9.
4. **Phonons do not diffuse and the exchange uses the normal-state rate.** Both as in Sidorova and Vodolazov.
   Annunziata includes a phonon diffusivity (0.1 cm²/s); over the escape time that would be about 30 nm.
5. **Absorbed fraction 1.** Simon's 25 % loss is for 5 nm NbN with a 10 ps escape; here the escape is modelled and
   takes 80 ps.
6. **Constant L_k**, dc-coupled load, no amplifier.
7. **The count level (20 %) and the window (400 τ0).**
8. Bath temperature 0.3 T_c and photon wavelength are settings, not assumptions (`--t0`, `--wavelength`).

Items 2–5 and 7 are one flag each in `thresholds.py`. §8 has a first test of some of them.

## 8 · What each ingredient does: a first test

One event just below and one just above the full model's threshold, at one site of each device, with one thing
changed at a time (2.0 ξ, 600 τ0 window, `results/variants/`). The full model's thresholds are 0.742 (double
column, edge channel) and 0.727 (plain wire, (−58, 0)). "count" / "no" is the outcome at the lower / upper bias.

| Variant | double channel (−80, 0) at 0.72 / 0.76 | plain wire (−58, 0) at 0.705 / 0.75 | Reading |
| --- | --- | --- | --- |
| full model | no / count | no / count | as the bisection found |
| current source (N_sq → ∞) | quiet / resistive | quiet / resistive | same within the bracket: at 750 squares the circuit shapes the pulse, it does not move the threshold |
| 3000 squares | no / no (load 0.13, rising) | no / no (load 0.14, rising) | the same resistive domain, but moving 20 % of the current takes four times longer. **The count level is tied to the window**; see below |
| no latent heat | count / count | count / count | thresholds lower by more than 0.02 in both devices |
| no relaxation heat (Joule only) | no / no | no / no | thresholds higher by more than 0.02 in both |
| no k_s (normal-state conduction everywhere) | no / count | no / count, 75 τ0 later | inside the bracket |
| linear GL ε(T), no latent heat | count / count | count / count | lower by more than 0.02 in both |
| τ_th = 0.6 ps | no / count | no / **no** | channel unchanged; plain wire higher by more than 0.02 |
| r0 = 10 nm | no / count | no / count | inside the bracket |

- **Latent heat and relaxation heat both matter, with opposite signs.** The previous model had neither. Their sizes
  are not yet measured: this test only says each exceeds 0.02 I_dep. A threshold run per variant would give them.
- **The two-fluid ε(T) matters** (together with the latent heat it brings): the linear GL form detects at lower
  bias, as expected from its steeper loss of critical current on heating.
- **The channel threshold does not care about the photon's details; the plain wire's does.** That is the cascade:
  once 40 nm is blocked, what follows is set by the geometry and the heat. It also means the comparison between the
  two devices at 1550 nm inherits the uncertainty of τ_th ([double400.md](double400.md) §2).

**The count level and the window.** At threshold the count arrives late, 520–780 τ0 after the start, close to the
800 τ0 limit, and with 3000 squares a resistive domain has moved only 13 % of the current by 600 τ0. So for long
wires the present limit would miss counts that are certain to come. Two ways out, not yet chosen: let an undecided
event run while it is resistive with a limit that scales with N_sq, or locate thresholds with the current-source
criterion (resistive at the end of the window, independent of N_sq) and use the circuit for the pulse shape and
for latching. The table says the two agree at 750 squares.
