"""
snspd -- shared simulation code for the nanohole SNSPD study (pyTDGL 0.9).

  units     film parameters and the unit system (K0, j_dep, tau0, Pearl length)
  geometry  strips and hole layouts, photon-hit grids, unit cells
  hotspot   photon hotspot models (kinetic, after Simon et al.; energy-conserving FEM diffusion; legacy)
  sim       device construction, detection criterion, bisection, current crowding, task runner
  analysis  threshold maps, efficiency curves, quantiles, selectivity
  store     resumable CSV / JSON results
  plotting  figure style and the standard figures
  fluxoid   trapped-flux meter (winding of arg psi around each hole)
  thermal   electro-thermal TDGL prototype (local temperature, Joule heat)

Scripts at the repository root (hitmap.py, fluxtrap.py) are thin command-line layers on top of this.
"""
