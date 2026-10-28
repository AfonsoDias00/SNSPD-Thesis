"""
Film parameters and the unit system shared by every script.

pyTDGL solves the generalized TDGL equations in dimensionless form [1]:
  lengths      in xi
  sheet current in K0 = 4 xi B_c2 / (mu0 Lambda) = 2 Phi0 d / (pi xi mu0 lambda^2),   Lambda = lambda^2 / d
  time         in tau0 = mu0 sigma lambda^2  (needs the normal-state conductivity sigma = 1 / rho_n)

Screening is off (pyTDGL's default, used throughout), so lambda enters ONLY these two scales: every normalised result
(i_det = I_det / I_sw, I_sw / I_dep) is independent of lambda, absolute currents scale as 1 / lambda^2 and absolute
times as lambda^2 / rho_n.  The GL depairing sheet current [2] in the same units is

  j_dep = K0 / (6 sqrt 3) = d Phi0 / (3 sqrt 3 pi mu0 lambda^2 xi)        (5.925 uA/nm for the default film)

and a uniform strip of width w has I_dep = j_dep * w.  Bisection tolerances and detection thresholds are expressed
as fractions of I_dep so they stay meaningful if lambda changes.

[1] L. Bishop-Van Horn, Comput. Phys. Commun. 291, 108799 (2023).
[2] M. Tinkham, Introduction to Superconductivity, 2nd ed. (1996), ch. 4.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

PHI0 = 2.067833848e-15     # flux quantum [Wb]
MU0 = 1.25663706212e-6     # vacuum permeability [H/m]


@dataclass(frozen=True)
class Film:
    """Material and model parameters of the superconducting film.

    The defaults are the values used for every result in results/.  rho_uohm_cm is only needed to quote times in
    picoseconds; it does not change any simulation.  Set it to the measured normal-state resistivity of the film.
    """
    xi_nm: float = 5.0             # coherence length
    lambda_nm: float = 175.0       # London penetration depth (a fit to 9 nm NbTiN hairpins gives ~600 nm; see README)
    d_nm: float = 9.0              # thickness
    gamma: float = 1.0             # GTDGL inelastic-scattering parameter (pyTDGL default is 10)
    u: float = 5.79                # GTDGL ratio of relaxation times (dirty limit, pyTDGL default)
    rho_uohm_cm: Optional[float] = None   # normal-state resistivity [uOhm cm]; None = times stay in tau0

    # ---- current scale ------------------------------------------------------------------------------------------
    @property
    def K0_uA_per_nm(self) -> float:
        """pyTDGL sheet-current unit K0 [uA/nm]."""
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
    def pearl_um(self) -> float:
        """Pearl length 2 lambda^2 / d [um]."""
        return 2.0 * self.lambda_nm ** 2 / self.d_nm * 1e-3

    # ---- time scale ---------------------------------------------------------------------------------------------
    @property
    def tau0_ps(self) -> Optional[float]:
        """pyTDGL time unit tau0 = mu0 sigma lambda^2 [ps], or None if rho_n is not set."""
        if self.rho_uohm_cm is None:
            return None
        rho = self.rho_uohm_cm * 1e-8                                       # uOhm cm -> Ohm m
        return MU0 * (self.lambda_nm * 1e-9) ** 2 / rho * 1e12

    def ps(self, t_tau0: float) -> Optional[float]:
        return None if self.tau0_ps is None else t_tau0 * self.tau0_ps

    def cm2_per_s(self, D_nm2_per_tau0: float) -> Optional[float]:
        return None if self.tau0_ps is None else D_nm2_per_tau0 * 1e-14 / (self.tau0_ps * 1e-12)

    # ---- pyTDGL objects -----------------------------------------------------------------------------------------
    def layer(self):
        import tdgl
        kw = {} if self.rho_uohm_cm is None else {"conductivity": 1.0 / (self.rho_uohm_cm * 1e-2)}  # S/um
        return tdgl.Layer(coherence_length=self.xi_nm * 1e-3, london_lambda=self.lambda_nm * 1e-3,
                          thickness=self.d_nm * 1e-3, gamma=self.gamma, u=self.u, **kw)

    def describe(self) -> list[str]:
        t = self.tau0_ps
        return [
            f"film: xi = {self.xi_nm:g} nm, lambda = {self.lambda_nm:g} nm, d = {self.d_nm:g} nm, gamma = {self.gamma:g}, "
            f"u = {self.u:g}, rho_n = {'not set' if self.rho_uohm_cm is None else f'{self.rho_uohm_cm:g} uOhm cm'}",
            f"units: K0 = {self.K0_uA_per_nm:.3f} uA/nm, j_dep = K0/(6 sqrt 3) = {self.j_dep_uA_per_nm:.3f} uA/nm, "
            f"Pearl length = {self.pearl_um:.2f} um, tau0 = {'n/a (set rho_n)' if t is None else f'{t:.4f} ps'}",
        ]


FILM = Film()
