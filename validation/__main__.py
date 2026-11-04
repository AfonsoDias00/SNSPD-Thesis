"""
python -m validation list                              the tests
python -m validation run 1a 1b 2a [--workers 6]       run test groups ("all" = every group, "quick" = 1a 1b 2a)
python -m validation report                            collect results/validation/<id>.json into summary.txt
"""
from __future__ import annotations

import argparse
import importlib
import time

from snspd import sim

from . import common as C

MODULES = {"1a": "diffusion", "1b": "epsilon", "1c": "disc_vs_hole", "1d": "photon",
           "2a": "meter", "2b": "ring", "2c": "single_hole", "2d": "scan_stats"}
QUICK = ("1a", "1b", "2a")            # no long TDGL runs, ~6 min


def group(gid):
    return importlib.import_module(f"validation.{MODULES[gid]}")


def cmd_list(a):
    for gid in MODULES:
        m = group(gid)
        sim.log(f"{gid}  {m.TITLE}")


def cmd_run(a):
    ids = list(MODULES) if "all" in a.ids else list(QUICK) if "quick" in a.ids else a.ids
    for gid in ids:
        m = group(gid)
        sim.log(f"== {gid}: {m.TITLE}")
        t0 = time.time()
        checks, extra = m.run(a)
        C.save_group(gid, m.TITLE, checks, time.time() - t0, extra)
        for c in checks:
            sim.log(f"   [{c.verdict}] {c.id} {c.name}: {c.value}   (expected {c.expected})")
        sim.log(f"   {time.time() - t0:.0f} s")
    cmd_report(a)


def cmd_report(a):
    groups = C.load_groups()
    L = ["VALIDATION SUMMARY   (results/validation/<id>.json; PASS / FAIL against the stated expectation, INFO = no "
         "pass mark)", ""]
    n = {"PASS": 0, "FAIL": 0, "INFO": 0}
    for gid, g in groups.items():
        L.append(f"{gid}  {g['title']}   ({g['wall_s']:.0f} s)")
        for c in g["checks"]:
            n[c["verdict"]] += 1
            L.append(f"  [{c['verdict']}] {c['id']:<6} {c['name']}")
            L.append(f"         measured: {c['value']}")
            L.append(f"         expected: {c['expected']}")
            if c["note"]:
                L.append(f"         note:     {c['note']}")
        L.append("")
    L.append(f"{n['PASS']} PASS, {n['FAIL']} FAIL, {n['INFO']} INFO")
    text = "\n".join(L)
    with open(C.out_path("summary.txt"), "w") as fh:
        fh.write(text + "\n")
    sim.log(f"summary written to {C.out_path('summary.txt')}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["list", "run", "report"])
    ap.add_argument("ids", nargs="*", default=["quick"])
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args(argv)
    {"list": cmd_list, "run": cmd_run, "report": cmd_report}[a.command](a)


if __name__ == "__main__":
    main()
