"""Command-line options shared by event.py and thresholds.py: they build one RunConfig and a tag for file names."""
from __future__ import annotations

import ast
from dataclasses import replace

from . import sim
from .circuit import Circuit
from .thermal import NBTIN_9NM, Photon, thermal_params


def add_model_args(ap):
    ap.add_argument("--mesh-xi", type=float, default=1.4, help="max mesh edge [xi]; 2.0 for screening")
    ap.add_argument("--tmax", type=float, default=400.0, help="window [tau0]")
    ap.add_argument("--tmax-ext", type=float, default=None, help="undecided events run on until here (default 2 tmax)")
    ap.add_argument("--nsq", type=float, default=750.0, help="series kinetic inductance [squares of film]")
    ap.add_argument("--source", action="store_true", help="current source instead of the circuit")
    ap.add_argument("--wavelength", type=float, default=1550.0, help="photon wavelength [nm]")
    ap.add_argument("--t0", type=float, default=0.3, help="bath temperature [T_c]")
    ap.add_argument("--thermal", default="", help="heat-model overrides, e.g. latent=False,ks=False,eps_model=linear")
    ap.add_argument("--photon", default="", help="photon overrides, e.g. r0_nm=10,tau_th_ps=0.6,absorbed=0.75")
    ap.add_argument("--label", default="", help="suffix of the result files for a model variant")


def _pairs(text):
    out = {}
    for item in filter(None, (s.strip() for s in text.split(","))):
        k, v = item.split("=", 1)
        try:
            out[k.strip()] = ast.literal_eval(v.strip())
        except (ValueError, SyntaxError):
            out[k.strip()] = v.strip()
    return out


def config(a, **kw) -> sim.RunConfig:
    th = thermal_params(NBTIN_9NM, t0=a.t0, **_pairs(a.thermal))
    ph = replace(Photon(wavelength_nm=a.wavelength), **_pairs(a.photon))
    return sim.RunConfig(mesh_xi=a.mesh_xi, tmax=a.tmax, tmax_ext=a.tmax_ext or 2 * a.tmax, thermal=th,
                         circuit=None if a.source else Circuit(n_sq=a.nsq), photon=ph, **kw)


def tag(a) -> str:
    """Identifies the model variant in file names and CSV rows."""
    parts = ["src" if a.source else f"n{a.nsq:g}", f"w{a.wavelength:g}", f"m{a.mesh_xi:g}"]
    if a.label:
        parts.append(a.label)
    return "_".join(parts)
