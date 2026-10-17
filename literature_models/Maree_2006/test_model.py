#!/usr/bin/env python3
"""Fast integration and analytic checks for the Marée reconstruction."""
from __future__ import annotations

import csv
import math
import pathlib
import shutil
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "test-output"


def thermal_ratchet_error() -> float:
    """Compare the target ratchet and implicit CPM calibration in spec §16."""
    errors = []
    for j in range(21):
        b = 0.05 + j * 0.01
        target = 500.0 * math.exp(-0.05 / b)
        lo, hi = 0.0, 500.0
        for _ in range(100):
            v = (lo + hi) / 2
            rhs = 100.0 / 0.1 * math.exp(-0.046 / 0.008) * (
                math.exp((1-v/500.0)*b/0.008)-math.exp(-(1-v/500.0)*b/0.008)
            )
            if rhs > v:
                lo = v
            else:
                hi = v
        errors.append(abs((lo + hi) / 2 - target) / target)
    return max(errors)


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    subprocess.run(
        [str(HERE / "maree2006"), "--nx", "64", "--ny", "64",
         "--initial-sites", "700", "--steps", "2", "--seed", "17",
         "--output", str(OUT), "--no-pulse", "--no-snapshots", "--initial-p-total", "0.2"],
        check=True,
    )
    with (OUT / "metrics.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    assert all(int(r["sites"]) > 600 for r in rows)
    assert int(rows[-1]["sites"]) != int(rows[0]["sites"]), "forced moving-boundary check made no copies"
    # Reactions exchange active/inactive pools, while geometry correction is conservative.
    for active, inactive in (("C_mass", "Ci_mass"), ("R_mass", "Ri_mass"), ("rho_mass", "rhoi_mass")):
        total0 = float(rows[0][active]) + float(rows[0][inactive])
        for row in rows[1:]:
            total = float(row[active]) + float(row[inactive])
            assert abs(total-total0) / total0 < 2e-10
    err = thermal_ratchet_error()
    # This checks the paper's calibration claim without pretending the match is exact.
    assert err < 0.09, err
    print(f"integration test passed; thermal-ratchet max relative discrepancy={err:.3%}")
    shutil.rmtree(OUT)


if __name__ == "__main__":
    main()
