# MSc thesis: SNSPDs with slot geometries (TDGL)

The repository has two halves that do not depend on each other:

| Folder | What | Start with |
| --- | --- | --- |
| [`simulations/`](simulations/) | Python code (`snspd/`, `validation/`, run scripts), results, and the working notes on the model (`docs/`) | [simulations/README.md](simulations/README.md) |
| [`report/`](report/) | The thesis in LaTeX (TU Delft template, XeLaTeX + biber) | [report/README.md](report/README.md) |

The Python environment `.venv/` (pyTDGL 0.9.0) stays at the root. Run the simulation scripts from inside
`simulations/`, as `../.venv/bin/python event.py ...`; all their paths (`results/`, `docs/`) are relative to it.

Build the report from inside `report/` with `latexmk report.tex`, or open the folder in VS Code and save a `.tex` file.

## From results to the report

Figures for the thesis are copied, not linked, from `simulations/results/` (or `simulations/docs/fig/`) to
`report/figures/`, so the report always builds on its own. Note the script and the result file a figure came
from in a `%` comment next to its `\includegraphics`.
