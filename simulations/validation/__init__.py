"""
validation -- tests of the model's parts against cases with known answers.

Each test group reports Checks (one number, its expected value, PASS / FAIL / INFO) to results/validation/<id>.json;
`python -m validation report` collects them in results/validation/summary.txt.

  2a  meter.py      the fluxoid meter on hand-built order parameters on real meshes
  2b  ring.py       a closed ring: winding conservation, persistent current, depairing limit, phase slips
  3a  heat.py       the heat model alone: energy audit, two-temperature relaxation, diffusion, walls
  3b  balance.py    energy balance of the TDGL: source work = change of free energy + Joule + relaxation heat
  3c  circuit.py    the readout circuit: LR recovery, resistance in squares, current divider

Groups 1a-1d of the earlier phases tested the epsilon-hotspot models, which this branch replaces by an energy
deposit in the heat equation; they are on the branch slot-geometries.
"""
