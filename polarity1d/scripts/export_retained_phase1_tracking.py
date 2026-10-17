#!/usr/bin/env python3
"""Export the compact tracking table consumed by retained Figure 1."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def number(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    return None if value in {"", "None", "null", "not_applicable"} else float(value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    arguments = parser.parse_args()
    with (arguments.output_root / "task_summary.csv").open(newline="") as handle:
        source = list(csv.DictReader(handle))
    rows = []
    for row in source:
        period = number(row, "tracking_period_s")
        lag = number(row, "phase_lag_seconds")
        rows.append({
            "task_id": row["task_id"],
            "model": row["model"],
            "coupling_variant": row.get("coupling_variant", "not_applicable"),
            "length_um": number(row, "length_um"),
            "period_s": period,
            "response_amplitude": number(row, "response_amplitude"),
            "phase_lag_seconds": lag,
            "phase_lag_fraction_of_period": None if lag is None or period is None else lag / period,
            "tracking_correlation": number(row, "tracking_correlation"),
            "tracking_reliable": row.get("tracking_reliable", "").lower() == "true",
        })
    rows.sort(key=lambda row: (row["model"], row["coupling_variant"], row["length_um"], row["period_s"]))
    destination = arguments.output_root / "jilkine_dynamic_tracking.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(destination)


if __name__ == "__main__":
    main()
