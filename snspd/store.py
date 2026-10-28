"""Resumable result files: one CSV row per finished threshold search, small JSON files for I_sw and crowding."""
from __future__ import annotations

import csv
import json
import os

CSV_FIELDS = ["Geometry", "x_nm", "y_nm", "energy_f", "tmax", "hotspot", "mesh_xi",
              "I_sw_uA", "I_det_uA", "i_det", "flag", "n_sims", "wall_s", "branch", "t_trigger"]


def row_key(geom, x, y, f):
    return (geom, round(float(x), 2), round(float(y), 2), round(float(f), 3))


def load_rows(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def append_row(path, row, fields=CSV_FIELDS):
    """Append one row.  A file written with an older column set is upgraded in place (new columns left blank)."""
    new = not os.path.exists(path)
    if not new:
        with open(path, newline="") as fh:
            rd = csv.DictReader(fh)
            header, old = rd.fieldnames or [], list(rd)
        if header != list(fields):
            missing = [k for k in header if k not in fields]
            if missing:
                raise ValueError(f"{path} has columns {missing} that are not in the expected set")
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields)
                w.writeheader()
                w.writerows({k: r.get(k, "") for k in fields} for r in old)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in fields})


def hit_row(geom_key, res, task, cfg, wall_s):
    return dict(Geometry=geom_key, x_nm=task["x_nm"], y_nm=task["y_nm"], energy_f=task["f"], tmax=cfg.tmax,
                hotspot=cfg.hotspot, mesh_xi=cfg.mesh_xi, I_sw_uA=res["I_sw"], I_det_uA=res["I_det"],
                i_det=res["I_det"] / res["I_sw"], flag=res["flag"], n_sims=res["n_sims"], wall_s=round(wall_s, 1),
                branch=res.get("branch", ""), t_trigger="" if res.get("t_trigger") is None else round(res["t_trigger"], 2))


def load_json(path, default):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return default


def save_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1)
