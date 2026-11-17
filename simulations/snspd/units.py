"""
Film parameters and the unit system shared by every script.

pyTDGL solves the generalized TDGL equations in dimensionless form [1]:
  lengths        in xi
  time           in tau0 = mu0 sigma lambda^2 = L_sq / R_sq      (sheet kinetic inductance over sheet resistance)
  sheet current  K = (K0 / 4) Im(psi* grad psi),  K0 = 4 xi B_c2 / (mu0 Lambda),  Lambda = lambda^2 / d
  potential      mu in V_u = Phi0 / (2 pi tau0) = hbar / (2 e tau0)   (one 2 pi phase slip = integral of mu dt = 2 pi)

pyTDGL's own "K0" and "V0" are 4 times these internal units (verified on pyTDGL 0.9.0: terminal currents are divided
by K0 / 4, and the equation carries i mu psi with mu = phi / V_u).  Two identities follow and are used by the circuit
(snspd/circuit.py) and the heat model (snspd/thermal.py):

  V_u / I_u = R_sq                 with I_u = (K0 / 4) xi: a normal region of n squares has mu-drop = n * I / I_u
  rho_n j_n^2 = (4 E_c / tau0) |grad mu|^2,   E_c = Phi0^2 / (16 pi^2 mu0 lambda^2 xi^2) = B_c^2 / (2 mu0)

Screening is off (pyTDGL's default, used throughout), so lambda enters ONLY the current, voltage and time units: every
normalised result is independent of lambda.  The GL depairing sheet current [2] in the same units is

  j_dep = K0 / (6 sqrt 3) = d Phi0 / (3 sqrt 3 pi mu0 lambda^2 xi)        (0.504 uA/nm for the default film)

and a uniform strip of width w has I_dep = j_dep * w.  Biases are quoted in I_dep of the narrowest free cross-section.

[1] L. Bishop-Van Horn, Comput. Phys. Commun. 291, 108799 (2023).
[2] M. Tinkham, Introduction to Superconductivity, 2nd ed. (1996), ch. 4.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PHI0 = 2.067833848e-15     # flux quantum [Wb]
MU0 = 1.25663706212e-6     # vacuum permeability [H/m]
KB = 1.380649e-23          # Boltzmann constant [J/K]
EV = 1.602176634e-19       # electron volt [J]
HC_EV_NM = 1239.841984     # h c [eV nm]: photon energy = HC_EV_NM / wavelength


@dataclass(frozen=True)
class Film:
    """The superconducting film: 9 nm NbTiN (KTH process).

    xi      5 nm.  Measured xi_GL(0) = 4.6 nm from B_c2 (Sidorova et al., PRB 104, 184514 (2021), Table 1).
    lambda  600 nm, fit to 9 nm NbTiN hairpins (Jonsson et al., PR Applied 17, 064046 (2022)).  The dirty-limit BCS
            value from the measured R_sq and T_c, sqrt(hbar rho_n / (pi mu0 Delta0)), is 630-660 nm.
    R_sq    381.4 Ohm, measured (Sidorova 2021 Table 1; 2022 Table 1).  With lambda it fixes tau0 = 0.132 ps.
    gamma   10: 2 tau_E Delta0 / hbar = 12-18 for the measured electron-phonon time; pyTDGL's default and the value
            of Simon et al. (PRB 112, 174512 (2025)).
    u       5.79, dirty limit.
    """
    xi_nm: float = 5.0
    lambda_nm: float = 600.0
    d_nm: float = 9.0
    gamma: float = 10.0
    u: float = 5.79
    R_sq_ohm: float = 381.4

    # ---- current -------------------------------------------------------------------------------------------------
    @property
    def K0_uA_per_nm(self) -> float:
        """pyTDGL's K0 [uA/nm] (4 times the internal sheet-current unit)."""
        xi, lam, d = self.xi_nm * 1e-9, self.lambda_nm * 1e-9, self.d_nm * 1e-9
        return 2.0 * PHI0 * d / (np.pi * xi * MU0 * lam ** 2) * 1e-3      # A/m -> uA/nm

    @property
    def j_dep_uA_per_nm(self) -> float:
        """GL depairing sheet current of the uniform film [uA/nm]."""
        return self.K0_uA_per_nm / (6.0 * np.sqrt(3.0))

    def I_dep_uA(self, width_nm: float) -> float:
        """Depairing current of a uniform strip (or of the narrowest free cross-section) of the given width."""
        return self.j_dep_uA_per_nm * width_nm

    @property
    def I_unit_uA(self) -> float:
        """Internal current unit I_u = (K0 / 4) xi [uA]: dimensionless current = I / I_u."""
        return self.K0_uA_per_nm / 4.0 * self.xi_nm

    @property
    def pearl_um(self) -> float:
        """Pearl length 2 lambda^2 / d [um]."""
        return 2.0 * self.lambda_nm ** 2 / self.d_nm * 1e-3

    # ---- time, voltage, inductance -------------------------------------------------------------------------------
    @property
    def L_sq_pH(self) -> float:
        """Sheet kinetic inductance mu0 lambda^2 / d [pH per square]."""
        return MU0 * (self.lambda_nm * 1e-9) ** 2 / (self.d_nm * 1e-9) * 1e12

    @property
    def tau0_ps(self) -> float:
        """pyTDGL time unit tau0 = mu0 sigma lambda^2 = L_sq / R_sq [ps]."""
        return self.L_sq_pH / self.R_sq_ohm

    @property
    def V_unit_uV(self) -> float:
        """Internal voltage unit Phi0 / (2 pi tau0) = R_sq I_u [uV]."""
        return self.R_sq_ohm * self.I_unit_uA

    @property
    def rho_uohm_cm(self) -> float:
        return self.R_sq_ohm * self.d_nm * 0.1

    def tau0(self, t_ps: float) -> float:
        """A time in ps expressed in tau0."""
        return t_ps / self.tau0_ps

    def nm2_per_tau0(self, D_cm2_s: float) -> float:
        """A diffusion constant in cm^2/s expressed in nm^2/tau0."""
        return D_cm2_s * 1e14 * self.tau0_ps * 1e-12

    # ---- energy --------------------------------------------------------------------------------------------------
    @property
    def E_c_J_m3(self) -> float:
        """GL condensation energy density B_c^2 / (2 mu0) implied by lambda and xi [J/m^3]."""
        lam, xi = self.lambda_nm * 1e-9, self.xi_nm * 1e-9
        return PHI0 ** 2 / (16 * np.pi ** 2 * MU0 * lam ** 2 * xi ** 2)

    # ---- pyTDGL objects ------------------------------------------------------------------------------------------
    def layer(self, gamma: float | None = None):
        import tdgl
        return tdgl.Layer(coherence_length=self.xi_nm * 1e-3, london_lambda=self.lambda_nm * 1e-3,
                          thickness=self.d_nm * 1e-3, gamma=self.gamma if gamma is None else gamma, u=self.u)

    def describe(self) -> list[str]:
        return [
            f"film: xi = {self.xi_nm:g} nm, lambda = {self.lambda_nm:g} nm, d = {self.d_nm:g} nm, R_sq = "
            f"{self.R_sq_ohm:g} Ohm ({self.rho_uohm_cm:.0f} uOhm cm), gamma = {self.gamma:g}, u = {self.u:g}",
            f"units: tau0 = {self.tau0_ps:.4f} ps, I_u = {self.I_unit_uA:.3f} uA, V_u = {self.V_unit_uV:.1f} uV, "
            f"L_sq = {self.L_sq_pH:.1f} pH, j_dep = {self.j_dep_uA_per_nm:.3f} uA/nm, E_c = {self.E_c_J_m3:.0f} J/m^3, "
            f"Pearl length = {self.pearl_um:.0f} um",
        ]


FILM = Film()
