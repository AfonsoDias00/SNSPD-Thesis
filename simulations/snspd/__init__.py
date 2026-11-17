"""
snspd -- shared simulation code for the double-slot-column SNSPD study (pyTDGL 0.9).

  units     film parameters and the unit system (tau0, I_u, V_u, L_sq, j_dep)
  geometry  the double column of 40 nm slots and the plain reference strips
  thermal   electron and phonon temperatures, the photon as an energy deposit, epsilon(T_e)
  circuit   series kinetic inductance and load resistor
  sim       device, settled start, one electro-thermal event, detection criterion, threshold bisection
  analysis  threshold maps, efficiency curves, quantiles
  store     resumable CSV results
  plotting  figure style and the standard figures
  fluxoid   trapped-flux meter (winding of arg psi around each hole)
"""
