# References: model basis, heat and circuit (10 Nov 2026)

This file collects the papers the simulations are built on. How each is used is in [model.md](model.md).

**Status column:**
- *read*: the text was checked in this project. Where only part was read verbatim, the part is named.
- *abstract*: only the abstract or search excerpt was checked.
- *memory*: a standard reference cited from memory. Check volume and page before citing it in the thesis.

On 10 Nov 2026 the papers marked **(re-read)** were read again from arXiv to build the model in
`snspd/thermal.py` and `snspd/circuit.py`: equations were extracted with their numbers, and anything a paper does
not state is recorded as such in [model.md](model.md) §2–3.

Numbers in [brackets] are the reference numbers used in docs/methods.md.

## 1 · Model basis (used in the code)

| Ref | Paper | Used for | Status |
|---|---|---|---|
| [1] | L. Bishop-Van Horn, "pyTDGL: Time-dependent Ginzburg-Landau in Python", *Comput. Phys. Commun.* **291**, 108799 (2023), [arXiv:2302.03812](https://arxiv.org/abs/2302.03812) | the solver. The heat equation is listed there only as a "possible extension" and is not implemented | read |
| [2] | L. Kramer and R. J. Watts-Tobin, *Phys. Rev. Lett.* **40**, 1041 (1978) | the generalised TDGL equation (γ, u = 5.79) | memory |
| [3] | R. J. Watts-Tobin, Y. Krähenbühl and L. Kramer, *J. Low Temp. Phys.* **42**, 459 (1981) | the generalised TDGL equation | memory |
| [27] | Simon *et al.*, "Ab initio modeling of nonequilibrium dynamics in superconducting detectors and qubits", *Phys. Rev. B* **112**, 174512 (2025), [arXiv:2501.13791](https://arxiv.org/abs/2501.13791) | the photon's build-up time τ_th and initial volume πξ_c²d; γ = 10; pyTDGL. Eq. (3) defines ε, Eq. (10) links it to α. **No value of ε_∞, no profile, no heat or circuit.** The outlook asks for "thermal equations for the electron and phonon temperatures" and "a circuit model" | read **(re-read:** main text, appendix and Table 1 from arXiv v3) |
| [8] | A. N. Zotova and D. Y. Vodolazov, "Photon detection by current-carrying superconducting film: A time-dependent Ginzburg-Landau approach", *Phys. Rev. B* **85**, 024509 (2012), [arXiv:1112.3790](https://arxiv.org/abs/1112.3790) | the photon as an energy in the heat equation (Eq. 1: hν = ΔT·πR²d·C_v), the Joule term (Eq. 2), bath temperature at the contacts, the circuit (Eq. 5). Its heat equation has one temperature, constant C_v and no latent heat; `snspd/thermal.py` goes beyond it | read **(re-read:** through Sec. IV; the conclusion was not read verbatim) |
| [39] | M. Sidorova *et al.*, *Phys. Rev. B* **104**, 184514 (2021), [arXiv:2108.04717](https://arxiv.org/abs/2108.04717) | measured 9 nm NbTiN films (KTH): R_sq, D, N(0), ξ(0), τ_e-ph and its exponent n, τ_esc, c_e/c_ph. τ_EP is the energy relaxation time, shorter than the scattering time τ_e-ph by a factor that depends on n | read **(re-read:** through Sec. III; appendices A–C only as a summary) |
| [40] | M. Sidorova *et al.*, *Supercond. Sci. Technol.* (2022), [DOI 10.1088/1361-6668/ac8454](https://doi.org/10.1088/1361-6668/ac8454), [arXiv:2204.01423](https://arxiv.org/abs/2204.01423) | **the heat model's numbers and form**: Table 1 (200 nm strips: T_c, D, τ_esc, τ_EP(T_c), n, c_e, c_ph) and the two-temperature model with T-dependent capacities, Eq. (6); Eq. (5) for I_dep(T); I_C/I_dep = 0.63 | read **(re-read:** through Eq. (6) in App. B.1; the rest of App. B only as a summary) |
| [13] | M. Jönsson *et al.*, *Phys. Rev. Applied* **17**, 064046 (2022), [arXiv:2112.05443](https://arxiv.org/abs/2112.05443) | λ = 600 nm (KTH NbTiN hairpins); crowding in GTDGL | read |
| [21] | D. Y. Vodolazov, *Phys. Rev. B* **85**, 174507 (2012), [arXiv:1204.0365](https://arxiv.org/abs/1204.0365) | the vortex-entry barrier (planned dark-count estimate) | read |
| [23] | S. Frasca *et al.*, *Phys. Rev. B* **100**, 054520 (2019), [arXiv:1904.08596](https://arxiv.org/abs/1904.08596) | measured I_sw / I_dep | read |

## 2 · TDGL with a heat equation and/or a circuit

| Paper | What it does | Status |
|---|---|---|
| Zotova & Vodolazov 2012 [8] | **The circuit:** the film's kinetic inductance L_k in series, with a 50 Ω shunt: (L_k/c²) dI_s/dt = (I − I_s) R_shunt − V_s (Eq. 5), with L_k = 4πλ²l/(wd) (Eq. 18). "When biased with the fixed current, the film latches in the normal state"; with the shunt the result is a voltage pulse of about 900 τ0 for a 200 µm film | read (re-read) |
| N. Yasukawa, T. Nishio and Y. Mawatari, "Ginzburg-Landau simulations of three-terminal operation of a superconducting nanowire cryotron", *Supercond. Sci. Technol.* (2024), [DOI 10.1088/1361-6668/ad44e4](https://doi.org/10.1088/1361-6668/ad44e4), [arXiv:2405.07466](https://arxiv.org/abs/2405.07466) | TDGL + heat equation (constant C, K, τ_el-ph; one temperature; linear T/T_c in the TDGL). **Heat sources: Joule heating and the dissipation from order-parameter relaxation** (Eq. 3, after Schmid 1966). Circuit: L_k dI/dt + V = R_L (I_bias − I) (Eq. 4), L_k = 1 nH, loads 50 Ω to 10 kΩ. Solved in COMSOL. An nTron, not a photon detector. The closest published TDGL + heat + circuit model | read (re-read: the full arXiv PDF) |
| G. R. Berdiyorov, M. V. Milošević and F. M. Peeters, "Spatially dependent sensitivity of superconducting meanders as single-photon detectors", *Appl. Phys. Lett.* **100**, 262603 (2012), [arXiv:1206.4298](https://arxiv.org/abs/1206.4298) | TDGL + heat diffusion at a meander turn; photons at the crowded corner do not switch | abstract |
| D. Y. Vodolazov *et al.*, *Phys. Rev. B* **71**, 184502 (2005), [arXiv:cond-mat/0502457](https://arxiv.org/abs/cond-mat/0502457) | TDGL + heat; hysteretic I-V from self-heating | abstract |
| Z. Jing, H. Yong and Y. Zhou, "Thermal coupling effect on the vortex dynamics of superconducting thin films: time-dependent Ginzburg-Landau simulations", *Supercond. Sci. Technol.* (2018), [DOI 10.1088/1361-6668/aab3be](https://doi.org/10.1088/1361-6668/aab3be) | TDGL + heat; vortex paths in films | abstract |
| L. R. Cadorim *et al.*, "Describing heat dissipation in the resistive state of three-dimensional superconductors" (2024), [arXiv:2405.10415](https://arxiv.org/abs/2405.10415) | GTDGL + Ampère + heat, in 3D | abstract |
| "Intrinsic timing jitter and latency in superconducting single photon nanowire detectors" (2018), [arXiv:1805.00130](https://arxiv.org/abs/1805.00130) | calls electrothermal TDGL the most advanced microscopic SNSPD model; latency and jitter | abstract |

No published SNSPD model that couples heat to pyTDGL was found (search of 10 Nov 2026).

## 3 · Latching, reset and lumped electrothermal models

| Paper | What it is for | Status |
|---|---|---|
| A. J. Kerman *et al.*, "Electrothermal feedback in superconducting nanowire single-photon detectors", *Phys. Rev. B* **79**, 100509(R) (2009), [arXiv:0812.0290](https://arxiv.org/abs/0812.0290) | the loop equation I_d R_n + L dI_d/dt = R_L(I_0 − I_d); hot-spot growth velocity (Eq. 1); latching when the electrical time τ_e = L/R_L is too short against the thermal time; the hot spot is always unstable when its resistance is below about R_L | read (10 Nov 2026, arXiv text; the journal reference was not on that page) |
| A. J. Annunziata *et al.*, "Reset dynamics and latching in niobium superconducting nanowire single-photon detectors", *J. Appl. Phys.* **108**, 084507 (2010), [arXiv:1008.0895](https://arxiv.org/abs/1008.0895) | two-temperature model with T-dependent heat capacities (Eqs. 1–2), the loop equation (Eq. 3), hν into one grid cell of the electron temperature (Eq. 8), GL I_c(T) (Eq. 7); ½L_kI_b² is dissipated in the hot spot; L_k varies ≈ 5 % with current in Nb | read (10 Nov 2026: pp. 1–12 of the arXiv PDF, through the model section) |
| A. J. Kerman *et al.*, "Kinetic-inductance-limited reset time of superconducting nanowire photon counters", *Appl. Phys. Lett.* **88**, 111116 (2006) | reset time τ = L_k / R_L | memory |
| J. K. W. Yang *et al.*, "Modeling the electrical and thermal response of superconducting nanowire single-photon detectors", *IEEE Trans. Appl. Supercond.* **17**, 581 (2007) | 1D electrothermal hotspot model + circuit, with T-dependent parameters | memory |
| K. K. Berggren *et al.*, "A superconducting nanowire can be modeled by using SPICE", *Supercond. Sci. Technol.* **31**, 055010 (2018) | SPICE model for the circuit after the TDGL window | memory |
| W. J. Skocpol, M. R. Beasley and M. Tinkham, "Self-heating hotspots in superconducting thin-film microbridges", *J. Appl. Phys.* **45**, 4054 (1974) | the self-heating hotspot and its minimum (hotspot) current; the physics of a latched domain | memory |

## 4 · Parallel nanowires (SNAPs): the long-slot limit

| Paper | What it is for | Status |
|---|---|---|
| M. Ejrnaes *et al.*, "A cascade switching superconducting single photon detector", *Appl. Phys. Lett.* **91**, 262509 (2007) | cascade switching of parallel wires: our cascade limit 1 − w_ch / w_free | memory |
| F. Marsili *et al.*, "Single-photon detectors based on ultranarrow superconducting nanowires", *Nano Lett.* **11**, 2048 (2011) | SNAPs with 20–30 nm wires: avalanche vs arm-trigger regime, series inductor against latching | memory |

## 5 · Beyond the model (limits of validity)

| Paper | What it is for | Status |
|---|---|---|
| D. Y. Vodolazov, "Single-photon detection by a dirty current-carrying superconducting strip based on the kinetic-equation approach", *Phys. Rev. Applied* **7**, 034014 (2017), [arXiv:1611.06060](https://arxiv.org/abs/1611.06060) (the arXiv title begins "Theory of single photon detection by 'dirty' …") | **the two-temperature model with energy conservation** (Eqs. 30–31: T_e⁵ − T_ph⁵ exchange, condensate term E₀𝓔_s in the electron energy), k_s (Eq. 32), I_det = I_dep(0)(1 − (T_e/T_c)²)^{3/2} (Eq. 29), the TDGL with modified coefficients for low T (Eq. 36), initial hot spot 2ξ_c × 2ξ_c. `snspd/thermal.py` uses its heat model with Sidorova's measured exponent; the modified TDGL coefficients are beyond pyTDGL | read (10 Nov 2026: through Sec. V from arXiv; Sec. VI and the conclusion only as a summary. The journal reference was not on that page) |
| Simon *et al.* [27], App. VII.2 | Usadel corrections to GTDGL at low T, and caution about defining electron and phonon temperatures | read |

## Status of the model extensions (10 Nov 2026)

Done on the branch `double400-models` ([model.md](model.md)):
1. **Circuit:** series inductance + 50 Ω load, in pyTDGL units with one parameter (squares of film in series).
2. **Temperature-dependent c_e, c_ph, τ_EP and ε(T):** Sidorova 2022 Eq. (6); two-fluid ε.
3. **Latent heat:** as the condensate entropy of the two-fluid free energy; and the order-parameter relaxation heat
   of Yasukawa *et al.*
4. **The photon as an energy** (Zotova & Vodolazov Eq. 1) instead of a calibrated ε.

Still beyond the model:
5. **Modified TDGL coefficients and the Usadel supercurrent** (Vodolazov 2017 Eq. 36; Simon *et al.* App. VII.2):
   pyTDGL cannot vary them in space. Quote as the limit of validity.
6. **Kinetic equations for the first picosecond** (Vodolazov 2017 Sec. III; Simon *et al.*).
7. **A lumped model for the reset** after the TDGL window (Berggren 2018, Yang 2007).
