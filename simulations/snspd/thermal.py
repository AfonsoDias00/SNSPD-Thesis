"""
Electron and phonon temperatures on the pyTDGL mesh: the heat model and the photon.

One model replaces the two of the earlier phases (an epsilon(r, t) "hotspot" with a calibrated amplitude, plus a
linearised Joule-heat field added to it).  The photon is now an amount of energy h nu put into the electrons; what
it does to the order parameter follows from the measured heat capacities.  docs/model.md derives everything below
and says which paper each piece comes from.

State on every mesh site (temperatures in units of T_c, t = T / T_c):

    U = t_e^2      electron thermal energy in units of gamma T_c^2 / 2        (c_e = gamma T, Drude)
    P = t_ph^4     phonon energy in units of c_ph(T_c) T_c / 4                (c_ph ~ T^3, Debye)

Equations (pyTDGL units: lengths xi, time tau0; n = |psi|^2 from the TDGL):

    T_e ds_e/dt = div(k_s grad T_e) + W - P_ep + photon,     s_e = gamma T_e (1 - b n)                       (1)
    d(P)/dt     = 4 r X - (P - P_0) / tau_esc                                                                (2)
    X     = (t_e^q - t_ph^q) / (q tau_ep)            electron-phonon exchange, q = 2 + n_ep                 (3)
    W     = eta0 [ |grad mu|^2 + u / sqrt(1 + gamma^2 n) ( |(d_t + i mu) psi|^2 + gamma^2/4 (d_t n)^2 ) ]   (4)
    eps   = (1 - U) / (1 - U_0)                      two-fluid: j_dep(T) ~ (1 - t^2)^(3/2), H_c ~ 1 - t^2    (5)

  (1)-(3) two-temperature model with T-dependent heat capacities and exchange: Sidorova et al., Supercond. Sci.
          Technol. (2022), DOI 10.1088/1361-6668/ac8454, Eq. (6), with their measured 9 nm NbTiN values; for q = 5 it is Eqs. (30)-(31) of Vodolazov,
          PR Applied 7, 034014 (2017), prefactors included.  Electron heat conduction k = c_e D (Wiedemann-Franz),
          reduced in the superconducting state by Vodolazov's Eq. (32).
  (4)     Joule heat of the normal current (Zotova & Vodolazov, PRB 85, 024509 (2012), Eq. (2)) plus the dissipation
          of order-parameter relaxation (Yasukawa et al., Supercond. Sci. Technol. (2024), arXiv:2405.07466, Eq. (3), after Schmid 1966), written for
          the generalised TDGL that pyTDGL solves.  Together they are exactly the rate at which the TDGL loses free
          energy, so no work done by the bias source goes missing (validation 3b).
  s_e     the entropy that belongs to the same free energy as the TDGL equation with (5): the condensate lowers it
          by b n.  When superconductivity is destroyed the electrons cool (latent heat), and they warm again when it
          recovers.  This is the condensation-energy term E_s(T_e, |Delta|) of Vodolazov 2017 Eq. (30) in
          Ginzburg-Landau form.
  photon  energy h nu into the electrons, Gaussian of radius r0 ~ xi, at the rate exp(-(t - t_a)/tau_th)/tau_th
          (Zotova & Vodolazov Eq. (1): the photon is an energy in the heat equation; Simon et al., PRB 112, 174512
          (2025): initial volume pi xi^2 d, thermalisation time tau_th).

Numerics: finite volumes on pyTDGL's own Voronoi cells (so energy sums are exact), operator splitting inside each
solver step: isentropic response to the change of n, heat sources, explicit conduction (sub-stepped when needed),
linearly implicit exchange and escape (stable for the stiff T^q exchange, and it moves exactly the same energy out of
the electrons as into the phonons).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .units import EV, FILM, HC_EV_NM, Film


# ==============================================================================
# Material and parameters
# ==============================================================================
@dataclass(frozen=True)
class Material:
    """Measured film properties that the heat model needs (SI-style units, as tabulated)."""
    name: str
    Tc_K: float
    D_cm2_s: float            # electron diffusion constant
    tau_ep_ps: float          # electron-phonon ENERGY relaxation time at T_c (c_e / tau_ep = dP_ep/dT_e)
    n_ep: float               # tau_ep ~ T^-n
    tau_esc_ps: float         # phonon escape time to the substrate
    ce_Tc: float              # electron heat capacity at T_c [J/(K m^3)]
    cph_Tc: float             # phonon heat capacity at T_c [J/(K m^3)]
    gap_ratio: float = 1.764  # Delta(0) / (k_B T_c); BCS (not measured for these films)
    tau_th_ps: float = 1.4    # thermalisation time of the photon's energy


# 9 nm NbTiN, 200 nm wide strip on SiO2/Si (KTH): Sidorova et al., Supercond. Sci. Technol. (2022), DOI
# 10.1088/1361-6668/ac8454, arXiv:2204.01423, Table 1.  c_e, c_ph and tau_EP are the film values of Sidorova et al., PRB 104, 184514 (2021)
# extrapolated to the strip's T_c with c_e ~ T, c_ph ~ T^3, tau_EP ~ T^-n.  tau_th is not measured for NbTiN: 1.4 ps
# is the ab initio value for NbN with D = 0.5 cm^2/s (Simon et al. 2025, Sec. III).
NBTIN_9NM = Material("9 nm NbTiN strip (Sidorova 2022, Table 1)", Tc_K=8.60, D_cm2_s=0.472, tau_ep_ps=4.0, n_ep=3.4,
                     tau_esc_ps=80.0, ce_Tc=1303.0, cph_Tc=8173.0)


@dataclass(frozen=True)
class ThermalParams:
    """The heat model in pyTDGL units (lengths nm or xi as named, times tau0, temperatures T_c)."""
    t0: float                 # bath temperature / T_c
    D_nm2: float              # electron heat diffusion [nm^2 / tau0]
    tau_ep: float             # electron-phonon energy relaxation time at T_c [tau0]
    q: float                  # exchange exponent, 2 + n_ep
    tau_esc: float            # phonon escape [tau0]
    r: float                  # c_e(T_c) / c_ph(T_c)
    eta0: float               # 4 E_c / (gamma T_c^2): TDGL dissipation -> electron energy
    e_unit_eV: float          # (gamma T_c^2 / 2) d xi^2 [eV]: the energy unit of U per xi^2 of film
    gap_ratio: float = 1.764
    eps_model: str = "two_fluid"      # "two_fluid": eps = (1 - t^2)/(1 - t0^2);  "linear": (1 - t)/(1 - t0) (GL)
    latent: bool = True               # condensate entropy (latent heat); needs the two-fluid eps
    ks: bool = True                   # reduced heat conduction in the superconducting state
    relax_heat: bool = True           # order-parameter relaxation heats the electrons, as well as Joule heat
    sink_nm: float = 50.0             # bath temperature within this distance of the terminals (see HeatField)
    eps_floor: float = -1.0           # numerical guard (pyTDGL underflows |psi|^3 in long-lived normal domains)

    @property
    def b(self) -> float:
        """Condensate entropy coefficient: s_e = gamma T (1 - b |psi|^2)."""
        return self.eta0 / (1.0 - self.t0 ** 2) if self.latent else 0.0

    def photon_energy(self, eV: float) -> float:
        """A photon energy in the units of U (per xi^2)."""
        return eV / self.e_unit_eV

    def describe(self, film: Film = FILM) -> list[str]:
        ps = film.tau0_ps
        return [
            f"heat: bath {self.t0:g} T_c; D = {self.D_nm2:.2f} nm^2/tau0; tau_ep(T_c) = {self.tau_ep:.1f} tau0 "
            f"({self.tau_ep * ps:.2f} ps), exchange ~ T_e^{self.q:g} - T_ph^{self.q:g}; tau_esc = {self.tau_esc:.0f} "
            f"tau0 ({self.tau_esc * ps:.0f} ps); c_e/c_ph(T_c) = {self.r:.3f}",
            f"      eta0 = {self.eta0:.3f}; energy unit {self.e_unit_eV * 1e3:.2f} meV per xi^2; eps(T): "
            f"{self.eps_model}; latent heat {'on' if self.latent else 'off'} (b = {self.b:.3f}); k_s "
            f"{'on' if self.ks else 'off'}; relaxation heat {'on' if self.relax_heat else 'off'}",
        ]


def thermal_params(material: Material = NBTIN_9NM, film: Film = FILM, t0: float = 0.3, **overrides) -> ThermalParams:
    """Dimensionless heat model of a material in the film's pyTDGL units.

    eta0 is not taken from lambda and xi but from thermodynamics: the two-fluid free energy must give a positive
    entropy, which ties the condensation energy to the heat capacity, 4 E_c(0) / (gamma T_c^2) = a0 <= 1, and BCS
    gives a0 = (3 / pi^2) (Delta0 / k_B T_c)^2 = 0.946.  eta0 = a0 (1 - t0^2)^2 is the same quantity at the bath
    temperature.  lambda = 600 nm, xi = 5 nm and the measured c_e give 1.03 instead: the two agree to 9 %, within
    the uncertainty of lambda^2 xi^2.  (eta0_from_film returns that value.)"""
    a0 = 3.0 / np.pi ** 2 * material.gap_ratio ** 2
    gamma_Tc2 = material.ce_Tc * material.Tc_K                                  # gamma T_c^2 [J/m^3]
    e_unit = 0.5 * gamma_Tc2 * (film.d_nm * 1e-9) * (film.xi_nm * 1e-9) ** 2 / EV
    p = ThermalParams(t0=t0, D_nm2=film.nm2_per_tau0(material.D_cm2_s), tau_ep=film.tau0(material.tau_ep_ps),
                      q=2.0 + material.n_ep, tau_esc=film.tau0(material.tau_esc_ps), r=material.ce_Tc / material.cph_Tc,
                      eta0=a0 * (1.0 - t0 ** 2) ** 2, e_unit_eV=e_unit, gap_ratio=material.gap_ratio)
    p = replace(p, **overrides)
    if p.eps_model not in ("two_fluid", "linear"):
        raise ValueError(p.eps_model)
    if p.latent and p.eps_model != "two_fluid":
        raise ValueError("latent heat needs eps_model='two_fluid': the linear GL entropy is negative far below T_c")
    if p.b >= 1.0:
        raise ValueError(f"b = {p.b:.3f} >= 1: the superconducting state would have negative entropy")
    return p


def eta0_from_film(material: Material = NBTIN_9NM, film: Film = FILM) -> float:
    """4 E_c / (gamma T_c^2) with E_c from the film's lambda and xi (the cross-check of thermal_params)."""
    return 4.0 * film.E_c_J_m3 / (material.ce_Tc * material.Tc_K)


NBTIN = thermal_params()


@dataclass(frozen=True)
class Photon:
    """One absorbed photon."""
    wavelength_nm: float = 1550.0
    absorbed: float = 1.0         # fraction of h nu that stays in the film while the energy thermalises
    r0_nm: float = 5.0            # 1/e radius of the Gaussian in which the energy appears (~ xi)
    tau_th_ps: float = NBTIN_9NM.tau_th_ps
    t_arrival: float = 20.0       # [tau0] after the start of the run

    @property
    def energy_eV(self) -> float:
        return self.absorbed * HC_EV_NM / self.wavelength_nm


# ==============================================================================
# Thermal conductivity in the superconducting state
# ==============================================================================
def _ks_table(y_max=40.0, n=4001):
    """k_s / k_n = 1 - (6 / pi^2) int_0^y x^2 e^x / (e^x + 1)^2 dx,  y = |Delta| / (k_B T_e)
    (Vodolazov 2017, Eq. (32))."""
    y = np.linspace(0.0, y_max, n)
    f = y ** 2 / (4.0 * np.cosh(np.minimum(y, 700.0) / 2.0) ** 2)
    integral = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(y))])
    return y, np.clip(1.0 - 6.0 / np.pi ** 2 * integral, 0.0, 1.0)


_KS_Y, _KS_F = _ks_table()


def ks_over_kn(y):
    return np.interp(y, _KS_Y, _KS_F)


def gap_at_bath(p: ThermalParams) -> float:
    """Delta(T_0) / (k_B T_c), BCS interpolation (Vodolazov 2017, Eq. (28))."""
    return p.gap_ratio * np.tanh(1.74 * np.sqrt(1.0 / p.t0 - 1.0))


# ==============================================================================
# The heat field
# ==============================================================================
AUDIT_KEYS = ("joule", "relax", "photon", "escape", "sink", "latent")


class HeatField:
    """U = t_e^2 and P = t_ph^4 on the mesh sites of a pyTDGL device, advanced after every solver step.

    Contacts.  pyTDGL injects the bias as normal current through psi = 0 terminals.  The conversion to supercurrent
    next to them dissipates steadily, which a real strip between superconducting leads does not do.  Sites within
    sink_nm of a terminal are therefore kept at the bath temperature and receive no heat (Zotova & Vodolazov hold
    T = T0 at the contacts and raise T_c over 5 xi next to them; Vodolazov 2017 holds T_e = T at the ends).  With
    gamma = 10 the conversion zone reaches about 40 nm into the strip: a 25 nm sink left the electrons at 0.41 T_c
    just beyond it in the dark, 50 nm leaves them at the bath temperature."""

    def __init__(self, device, p: ThermalParams, film: Film = FILM, photon: Photon | None = None, site_nm=None):
        self.p, self.film = p, film
        mesh, em = device.mesh, device.mesh.edge_mesh               # dimensionless (xi)
        self.xi_nm = film.xi_nm
        self.e0, self.e1 = np.asarray(em.edges[:, 0]), np.asarray(em.edges[:, 1])
        self.w = np.asarray(em.dual_edge_lengths, float) / np.asarray(em.edge_lengths, float)
        self.area = np.asarray(mesh.areas, float)                    # Voronoi cell areas [xi^2]
        self.N = len(self.area)
        self.xy_nm = np.asarray(mesh.sites, float) * film.xi_nm
        ay = np.abs(self.xy_nm[:, 1])
        self.sink = ay >= ay.max() - p.sink_nm if p.sink_nm > 0 else np.zeros(self.N, bool)
        self.live = ~self.sink
        self.wD = self.w * p.D_nm2 / film.xi_nm ** 2                 # edge conductance for U [xi^2 / tau0 scale]
        self.gap0 = gap_at_bath(p)
        self.U0, self.P0 = p.t0 ** 2, p.t0 ** 4
        self.U = np.full(self.N, self.U0)
        self.P = np.full(self.N, self.P0)
        self.n = np.ones(self.N)                                     # |psi|^2 that the entropy was last evaluated with
        self.g = 1.0 - p.b * self.n
        self.eps = np.ones(self.N)
        self.audit = dict.fromkeys(AUDIT_KEYS, 0.0)                  # cumulative energies [units of U per xi^2]
        self.te_peak = self.tp_peak = p.t0
        self._e0 = float(self.area[self.live] @ (self.g * self.U)[self.live])
        self.G = None
        if photon is not None and site_nm is not None:
            self.set_photon(photon, site_nm)

    def __getstate__(self):                       # keep pickles small if a closure over this object is ever stored
        return {"U": self.U, "P": self.P}

    def __setstate__(self, state):
        self.__dict__.update(state)

    # ---- set-up --------------------------------------------------------------------------------------------------
    def set_photon(self, photon: Photon, site_nm):
        r2 = (self.xy_nm[:, 0] - site_nm[0]) ** 2 + (self.xy_nm[:, 1] - site_nm[1]) ** 2
        if r2.min() > (2.0 * self.xi_nm) ** 2:
            raise ValueError(f"photon site {site_nm} nm is not on the film")
        G = np.exp(-r2 / photon.r0_nm ** 2)
        G[self.sink] = 0.0
        self.G = G / float(self.area @ G)                            # sum(area * G) = 1: the energy is exact
        self.E_photon = self.p.photon_energy(photon.energy_eV)
        self.t_a, self.tau_th = photon.t_arrival, self.film.tau0(photon.tau_th_ps)

    def reset(self, psi):
        """Bath temperature everywhere, in equilibrium with the order parameter psi."""
        self.U[:], self.P[:] = self.U0, self.P0
        self.n = np.minimum(np.abs(np.asarray(psi)) ** 2, 1.0)
        self.g = 1.0 - self.p.b * self.n
        self.eps = np.ones(self.N)
        self.audit = dict.fromkeys(AUDIT_KEYS, 0.0)
        self.te_peak = self.tp_peak = self.p.t0
        self._e0 = float(self.area[self.live] @ (self.g * self.U)[self.live])

    # ---- pieces --------------------------------------------------------------------------------------------------
    def _scatter(self, edge_values):
        """Sum of an antisymmetric edge quantity into the sites: +value at e0, -value at e1."""
        return np.bincount(self.e0, edge_values, self.N) - np.bincount(self.e1, edge_values, self.N)

    def joule(self, mu):
        """|grad mu|^2 per site: sum_e w (d mu)^2 is the integral of |grad mu|^2, half an edge to each end."""
        pe = 0.5 * self.w * (mu[self.e0] - mu[self.e1]) ** 2
        return (np.bincount(self.e0, pe, self.N) + np.bincount(self.e1, pe, self.N)) / self.area

    def relaxation(self, dt, psi_old, mu_old, psi_new):
        """u / sqrt(1 + gamma^2 n) ( |(d_t + i mu) psi|^2 + gamma^2/4 (d_t n)^2 ) per site, with the gauge-invariant
        time derivative that pyTDGL's own step uses, psi_new - exp(-i mu dt) psi_old."""
        f = self.film
        n_old = np.abs(psi_old) ** 2
        d = (psi_new - np.exp(-1j * mu_old * dt) * psi_old) / dt
        dn = (np.abs(psi_new) ** 2 - n_old) / dt
        return f.u / np.sqrt(1.0 + f.gamma ** 2 * n_old) * (np.abs(d) ** 2 + 0.25 * f.gamma ** 2 * dn ** 2)

    def conduct(self, dt):
        """Explicit heat conduction, g dU/dt = D div(F grad U) with F = k_s / k_n (harmonic mean on each edge)."""
        if self.p.ks:
            F = ks_over_kn(self.gap0 * np.sqrt(self.n / self.U))
            F = np.maximum(F, 1e-6)
            k = self.wD * 2.0 * F[self.e0] * F[self.e1] / (F[self.e0] + F[self.e1])
        else:
            k = self.wD
        cap = self.area * self.g
        rate = (np.bincount(self.e0, k, self.N) + np.bincount(self.e1, k, self.N)) / cap
        m = max(1, int(np.ceil(dt * rate.max() / 0.7)))
        h = dt / m
        for _ in range(m):
            self.U += h * self._scatter(k * (self.U[self.e1] - self.U[self.e0])) / cap
        return m

    def exchange(self, dt):
        """Electron-phonon exchange and phonon escape, linearly implicit in the exchange power X."""
        p, U, P = self.p, self.U, self.P
        te_q, tp_q = U ** (0.5 * p.q), P ** (0.25 * p.q)
        X = (te_q - tp_q) / (p.q * p.tau_ep)
        XU, XP = te_q / (2.0 * p.tau_ep * U), -tp_q / (4.0 * p.tau_ep * P)
        L = (P - self.P0) / p.tau_esc
        Xn = (X - XP * dt * L) / (1.0 + 2.0 * XU * dt / self.g - 4.0 * p.r * XP * dt)
        self.U = U - 2.0 * Xn * dt / self.g
        self.P = P + dt * (4.0 * p.r * Xn - L)
        self.audit["escape"] += float(self.area @ L) * dt / (2.0 * p.r)

    # ---- one solver step -----------------------------------------------------------------------------------------
    def step(self, t, dt, psi_old, mu_old, psi_new, mu_new):
        """Advance from t to t + dt, given the TDGL fields before and after the solver step."""
        p, A = self.p, self.area
        mu_new, mu_old = np.asarray(mu_new).real, np.asarray(mu_old).real
        # 1. latent heat: the entropy gamma t g is conserved while |psi|^2 changes, so t^2 scales as 1 / g^2
        if p.latent:
            n_new = np.minimum(np.abs(psi_new) ** 2, 1.0)
            g_new = 1.0 - p.b * n_new
            before = self.g * self.U
            self.U *= (self.g / g_new) ** 2
            self.U[self.sink] = self.U0
            self.audit["latent"] += float(A[self.live] @ (g_new * self.U - before)[self.live])
            self.n, self.g = n_new, g_new
        else:
            self.n = np.minimum(np.abs(psi_new) ** 2, 1.0)
        # 2. heat sources: dissipation of the TDGL and the photon
        wj = self.joule(mu_new)
        wj[self.sink] = 0.0
        src = 2.0 * p.eta0 * dt * wj
        self.audit["joule"] += float(A @ src)
        if p.relax_heat:
            wr = self.relaxation(dt, psi_old, mu_old, psi_new)
            wr[self.sink] = 0.0
            self.audit["relax"] += 2.0 * p.eta0 * dt * float(A @ wr)
            src += 2.0 * p.eta0 * dt * wr
        if self.G is not None and t + dt > self.t_a:
            a, b = max(t - self.t_a, 0.0), t + dt - self.t_a
            frac = np.exp(-a / self.tau_th) - np.exp(-b / self.tau_th)
            src += self.E_photon * frac * self.G
            self.audit["photon"] += self.E_photon * frac
        self.U += src / self.g
        # 3. conduction, 4. exchange with the phonons and escape
        self.conduct(dt)
        self.exchange(dt)
        # 5. contacts at the bath temperature
        if self.sink.any():
            s = self.sink
            self.audit["sink"] += float(A[s] @ (self.g[s] * (self.U[s] - self.U0) + (self.P[s] - self.P0) / (2.0 * p.r)))
            self.U[s], self.P[s] = self.U0, self.P0
        np.maximum(self.U, 1e-6, out=self.U)
        np.maximum(self.P, 1e-12, out=self.P)
        # 6. the TDGL's linear coefficient
        if p.eps_model == "two_fluid":
            eps = (1.0 - self.U) / (1.0 - self.U0)
        else:
            eps = (1.0 - np.sqrt(self.U)) / (1.0 - p.t0)
        self.eps = np.maximum(eps, p.eps_floor)
        self.te_peak = max(self.te_peak, float(np.sqrt(self.U.max())))
        self.tp_peak = max(self.tp_peak, float(self.P.max() ** 0.25))

    # ---- read-outs -----------------------------------------------------------------------------------------------
    @property
    def te(self):
        return np.sqrt(self.U)

    @property
    def tph(self):
        return self.P ** 0.25

    def energy(self):
        """(electron thermal energy, phonon energy) gained since reset(), summed over the film outside the sinks
        [units of U per xi^2].  The electron thermal energy density is g U = t_e^2 (1 - b |psi|^2).  The audit
        closes as  electrons + phonons = joule + relax + photon + latent - escape - sink."""
        A = self.area[self.live]
        return (float(A @ (self.g * self.U)[self.live]) - self._e0,
                float(A @ (self.P - self.P0)[self.live]) / (2.0 * self.p.r))

    def normal_area_nm2(self):
        """Area with T_e >= T_c."""
        return float(self.area[self.U >= 1.0].sum()) * self.xi_nm ** 2
