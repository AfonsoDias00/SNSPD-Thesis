"""Resumable result files: one CSV row per finished threshold search."""
from __future__ import annotations

import csv
import json
import os

CSV_FIELDS = ["Geometry", "x_nm", "y_nm", "model", "i_det", "flag", "n_runs", "t_count", "te_peak", "tph_peak",
              "load_max", "normal_nm2", "wall_s"]


def row_key(geom, x, y):
    """Key of a row; x = y = "" (or None) is the dark run of that geometry."""
    blank = x in ("", None)
    return (geom, None if blank else round(float(x), 2), None if blank else round(float(y), 2))


def load_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def append_row(path, row, fields=CSV_FIELDS):
    new = not os.path.exists(path)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in fields})


def threshold_row(geom_key, site, model, res, wall_s):
    ev = res.get("event") or {}
    r3 = lambda v: "" if v is None else round(float(v), 3)
    return dict(Geometry=geom_key, x_nm="" if site is None else site[0], y_nm="" if site is None else site[1],
                model=model, i_det=round(res["i_det"], 4), flag=res["flag"], n_runs=res["n_runs"],
                t_count=r3(ev.get("t_count")), te_peak=r3(ev.get("te_peak")), tph_peak=r3(ev.get("tph_peak")),
                load_max=r3(ev.get("load_max")), normal_nm2=r3(ev.get("normal_nm2")), wall_s=round(wall_s, 1))


def save_json(path, obj):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=float)
