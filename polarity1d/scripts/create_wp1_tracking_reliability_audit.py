#!/usr/bin/env python3
"""Create retained trajectories for a complete audit of panel-I tracking conditions."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


MODELS = ("wave_pinning", "otsuji", "goryachev", "legi")


def build_manifest() -> dict:
    source = load_manifest(ROOT / "campaigns" / "phase1_production_v1.json")
    tasks = []
    for original in source["tasks"]:
        if original["assay"] != "dynamic_tracking" or original["model"] not in MODELS:
            continue
        task = copy.deepcopy(original)
        task["tier"] = "tracking_reliability_audit"
        task["save_trajectory"] = True
        period = float(task["protocol"]["period"])
        task["simulation"]["output_interval"] = period / 80.0
        tasks.append(task)
    return make_manifest(
        "wp1_tracking_reliability_audit_v1", tasks,
        description=(
            "Retained trajectories for all 64 Jilkine model-by-length-by-period "
            "conditions shown in WP1 Figure 1I."
        ),
        metadata={
            "models": list(MODELS), "condition_count": 64,
            "purpose": "audit harmonic phase and reliability classification",
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "wp1_tracking_reliability_audit_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
