#!/usr/bin/env python3
"""Create the minimal Phase-1 subset needed by retained Figure 1."""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / "retained_phase1_tracking_v1.json",
    )
    arguments = parser.parse_args()
    source = load_manifest(ROOT / "campaigns" / "phase1_production_v1.json")
    tasks = [
        copy.deepcopy(task) for task in source["tasks"]
        if task["tier"] == "production"
        and task["assay"] == "dynamic_tracking"
        and task["model"] in {"wave_pinning", "otsuji", "goryachev", "legi"}
    ]
    manifest = make_manifest(
        "retained_phase1_tracking_v1", tasks,
        description="Minimal 64-condition Phase-1 tracking subset used by retained Figure 1.",
        metadata={"source_campaign": "phase1_production_v1", "condition_count": len(tasks)},
    )
    destination = write_manifest(arguments.output, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
