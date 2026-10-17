#!/usr/bin/env python3
"""Compare matching final snapshots from two Marée runs."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

FIELDS = ("C", "R", "rho", "Ci", "Ri", "rhoi", "A", "F_total", "B_total", "P_total")


def final_snapshot(run: Path) -> tuple[int, dict[tuple[int, int], dict[str, float]]]:
    path = run / "snapshots.csv"
    final = -1
    data: dict[tuple[int, int], dict[str, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            step = int(row["step"])
            if step > final:
                final, data = step, {}
            if step == final:
                data[(int(row["q"]), int(row["r"]))] = {name: float(row[name]) for name in FIELDS}
    return final, data


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    ref_step, ref = final_snapshot(args.reference)
    cand_step, cand = final_snapshot(args.candidate)
    common = sorted(ref.keys() & cand.keys())
    if not common:
        raise SystemExit("final snapshots have no sites in common")
    result: dict[str, object] = {
        "reference_step": ref_step, "candidate_step": cand_step,
        "reference_sites": len(ref), "candidate_sites": len(cand), "common_sites": len(common),
        "fields": {},
    }
    fields = result["fields"]
    assert isinstance(fields, dict)
    for name in FIELDS:
        a = [ref[key][name] for key in common]
        b = [cand[key][name] for key in common]
        diff2 = sum((x-y)**2 for x,y in zip(a,b))
        norm2 = sum(x*x for x in a)
        ma = sum(a)/len(a); mb = sum(b)/len(b)
        covariance = sum((x-ma)*(y-mb) for x,y in zip(a,b))
        va = sum((x-ma)**2 for x in a); vb = sum((y-mb)**2 for y in b)
        fields[name] = {
            "relative_l2": math.sqrt(diff2/max(norm2,1e-300)),
            "max_absolute": max(abs(x-y) for x,y in zip(a,b)),
            "correlation": covariance/math.sqrt(va*vb) if va and vb else (1.0 if diff2 == 0 else None),
        }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
