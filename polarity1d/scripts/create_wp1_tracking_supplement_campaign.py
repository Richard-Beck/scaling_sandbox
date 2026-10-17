#!/usr/bin/env python3
"""Create retained wave-pinning trajectories for the tracking supplement."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


def build_manifest() -> dict:
    source = load_manifest(ROOT / "campaigns" / "phase1_production_v1.json")
    tasks = []
    for period in (100.0, 400.0):
        match = next(
            task for task in source["tasks"]
            if task["assay"] == "dynamic_tracking" and task["model"] == "wave_pinning"
            and float(task["simulation"]["length"]) == 15.0
            and float(task["protocol"]["period"]) == period
        )
        task = copy.deepcopy(match)
        task["tier"] = "figure_supplement"
        task["save_trajectory"] = True
        task["simulation"]["output_interval"] = period / 80.0
        tasks.append(task)
    return make_manifest(
        "wp1_tracking_supplement_v1", tasks,
        description="Retained L=15 um wave-pinning tracking trajectories at T=100 and 400 s.",
        metadata={"length_um": 15.0, "periods_s": [100.0, 400.0]},
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "wp1_tracking_supplement_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
