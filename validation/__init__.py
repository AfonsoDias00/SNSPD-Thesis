"""
validation -- tests of the hotspot model and the fluxoid analysis against cases with known answers.

Each test group reports Checks (one number, its expected value, PASS / FAIL / INFO) to results/validation/<id>.json;
`python -m validation report` collects them in results/validation/summary.txt.

  1a  diffusion.py      the hotspot diffusion solver alone: conservation, free space, a wall, strip modes, convergence
  1b  epsilon.py        how pyTDGL treats epsilon < -1, and the epsilon range the production hotspot reaches
  1c  disc_vs_hole.py   a frozen normal disc against a real hole of the same radius: switching current
  1d  photon.py         plain-wire thresholds vs photon energy, position, width, arrival time, mirror symmetry
  2a  meter.py          the fluxoid meter on hand-built order parameters on real meshes
  2b  ring.py           a closed ring: winding conservation, persistent current, depairing limit, phase slips
  2c  single_hole.py    one hole in the plain strip: phase per trapped quantum, persistence, switching current
  2d  scan_stats.py     fluxoid steps in the existing flux scans, and a dense-frame re-run of busy runs

Not included: 1e, a re-run of a published pyTDGL hotspot simulation.  The candidate (Simon et al., PRB 112, 174512
(2025), arXiv:2501.13791) does not state the TDGL xi and lambda, the epsilon amplitude and profile, the strip length,
or the voltage threshold, and publishes no code.
"""
