"""
The readout circuit: the rest of the detector as a series kinetic inductance, and a load resistor.

                 I_b ->--+-----------+
                         |           |
                        L_k         R_L      a current source I_b feeds the detector (L_k in series with the
                         |           |       simulated strip, voltage V across it) in parallel with the load
                     [ TDGL: V ]     |
                         |           |
                         +-----------+

    L_k dI/dt + V = R_L (I_b - I)

This one equation is Eq. (5) of Zotova & Vodolazov (PRB 85, 024509 (2012)), Eq. (4) of Yasukawa et al. (2024,
arXiv:2405.07466), the loop equation of Kerman et al. (PRB 79, 100509 (2009)) and Eq. (3) of Annunziata et al.
(J. Appl. Phys. 108, 084507 (2010)).  The load is dc-coupled, as in those four.

In pyTDGL's units (snspd/units.py: V_u / I_u = R_sq and tau0 = L_sq / R_sq) it has no free material constant:

    dI/dt = [ (R_L / R_sq) (I_b - I) - I_u v ] / N_sq          t in tau0, v = V / V_u (pyTDGL's mu difference)

with N_sq = L_k / L_sq, the kinetic inductance counted in squares of film, and R_L / R_sq = 50 / 381.4 measured.
So the circuit is set by one number that describes the device: how many squares of nanowire are in series.

  * 750 squares (the default) is a 200 nm x 150 um strip, the devices of Sidorova et al. 2022: L_k = 38 nH,
    L_k / R_L = 0.75 ns = 5700 tau0.  A 16 x 16 um meander of 200 nm strip at 50 % fill is 3200 squares (160 nH).
  * The simulated strip's own kinetic inductance is already in the TDGL; N_sq is the rest.
  * N_sq -> infinity is the current source of the earlier phases.

What the circuit decides: a normal region of n squares takes the current out of the detector in a time N_sq / n
(in tau0), so with N_sq in the thousands the bias is constant while a hot spot forms and the circuit only shapes
the pulse and the reset (Kerman 2009; Annunziata 2010).  One 2 pi phase slip moves Phi0 / L_k = 2 pi I_u / N_sq out
of the detector: 0.05 uA for 750 squares.  A burst of a few slips is therefore invisible at the load; only a
resistive domain makes a pulse.  That is the detection criterion used in snspd/sim.py.

Not modelled: the dependence of L_k on current and temperature (Annunziata 2010 measured 5 % in Nb up to I_c), the
ac coupling of a real bias tee, amplifier bandwidth and reflections.
"""
from __future__ import annotations

from dataclasses import dataclass

from .units import FILM, Film


@dataclass(frozen=True)
class Circuit:
    n_sq: float = 750.0        # series kinetic inductance [squares of film]
    RL_ohm: float = 50.0       # load

    def L_k_nH(self, film: Film = FILM) -> float:
        return self.n_sq * film.L_sq_pH * 1e-3

    def tau_reset(self, film: Film = FILM) -> float:
        """L_k / R_L [tau0]."""
        return self.n_sq * film.R_sq_ohm / self.RL_ohm

    def describe(self, film: Film = FILM) -> list[str]:
        t = self.tau_reset(film)
        return [f"circuit: {self.n_sq:g} squares in series (L_k = {self.L_k_nH(film):.1f} nH), R_L = {self.RL_ohm:g} Ohm; "
                f"L_k/R_L = {t:.0f} tau0 ({t * film.tau0_ps * 1e-3:.2f} ns); one phase slip moves "
                f"{2 * 3.141592653589793 * film.I_unit_uA / self.n_sq:.3f} uA"]


class CircuitState:
    """The detector current I [uA], advanced with the strip's voltage after every solver step."""

    def __init__(self, circuit: Circuit, film: Film, I_bias_uA: float):
        self.k = circuit.RL_ohm / film.R_sq_ohm
        self.n_sq = circuit.n_sq
        self.I_u = film.I_unit_uA
        self.I_bias = self.I = float(I_bias_uA)

    def step(self, dt: float, v: float):
        """Backward Euler in the load term, the strip voltage v [V_u] taken from the step just made."""
        a = dt / self.n_sq
        self.I = (self.I + a * (self.k * self.I_bias - self.I_u * v)) / (1.0 + a * self.k)

    @property
    def load_fraction(self) -> float:
        """Share of the bias current flowing in the load, (I_b - I) / I_b: the pulse, as a fraction of I_b R_L."""
        return 1.0 - self.I / self.I_bias

    def currents(self, t=None):
        return {"source": self.I, "drain": -self.I}
