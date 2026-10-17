#!/usr/bin/env python3
"""Create fixed-20,000-s competition and reversal tasks for WP1 Figure 1."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


def build_manifest() -> dict:
    competition_manifest = load_manifest(ROOT / "campaigns" / "wp1_revision_audit_v1.json")
    reversal_manifest = load_manifest(ROOT / "campaigns" / "jilkine_reversal_long_conditioning_v1.json")
    tasks = []
    for source in competition_manifest["tasks"]:
        if source["assay"] != "competing_fronts_extended_4x":
            continue
        task = copy.deepcopy(source)
        task["tier"] = "figure1_common_horizon"
        task["simulation"]["end_time"] = 20_000.0
        task["save_trajectory"] = False
        tasks.append(task)
    for source in reversal_manifest["tasks"]:
        task = copy.deepcopy(source)
        task["tier"] = "figure1_common_horizon"
        reversal_time = float(task["protocol"]["reversal_time"])
        task["simulation"]["end_time"] = reversal_time + 20_000.0
        task["save_trajectory"] = False
        tasks.append(task)
    return make_manifest(
        "wp1_figure1_common_horizon_v1", tasks,
        description=(
            "Cross-model competition and reversal assays with a common 20,000-s "
            "observation horizon for WP1 Figure 1."
        ),
        metadata={
            "competition_end_time_s": 20_000.0,
            "reversal_followup_time_s": 20_000.0,
            "source_campaigns": [
                "wp1_revision_audit_v1", "jilkine_reversal_long_conditioning_v1"
            ],
        },
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "wp1_figure1_common_horizon_v1.json",
    )
    arguments = parser.parse_args(argv)
    manifest = build_manifest()
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
