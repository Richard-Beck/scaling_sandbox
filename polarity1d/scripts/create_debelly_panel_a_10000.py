#!/usr/bin/env python3
"""Create the 10,000-s panel-A campaign from the Julia-activation-default base."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


SOURCE_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_base"
CAMPAIGN_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000"


def main() -> None:
    source_path = ROOT / "campaigns" / f"{SOURCE_NAME}.json"
    source = json.loads(source_path.read_text())
    tasks = []
    for source_task in source["tasks"]:
        if source_task["assay"] != "single_stimulus_length_test":
            continue
        task = copy.deepcopy(source_task)
        for key in ("task_id", "task_index", "task_digest"):
            task.pop(key, None)
        task.update(
            implementation_revision="julia_activation_defaults_panel_a_10000_v1",
            source_campaign=SOURCE_NAME,
            source_assay="single_stimulus_length_test",
            save_trajectory=True,
        )
        task["simulation"]["end_time"] = 10_000.0
        tasks.append(task)
    if len(tasks) != 80:
        raise RuntimeError(f"expected 80 panel-A tasks, found {len(tasks)}")
    manifest = make_manifest(
        CAMPAIGN_NAME,
        tasks,
        description=(
            "Figure 2 panel-A single-stimulus trajectories extended to an absolute "
            "simulation end time of 10,000 s."
        ),
        metadata={
            "source_campaign": SOURCE_NAME,
            "simulation_end_time_s": 10_000.0,
            "changed_setting_only": "simulation.end_time",
        },
    )
    destination = ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json"
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
