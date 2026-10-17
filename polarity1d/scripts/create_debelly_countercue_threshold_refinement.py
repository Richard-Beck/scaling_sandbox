#!/usr/bin/env python3
"""Create one adaptive bisection round for short-counter-cue thresholds."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


SOURCE_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_short_cues_response"
REFINEMENT_PREFIX = "debelly_2026_standalone_figure2_julia_activation_defaults_threshold_refinement_r"
SOURCE_MANIFEST = ROOT / "campaigns" / f"{SOURCE_NAME}.json"
RESULTS = ROOT / "results/campaigns"
RESOLUTION_S = 0.5
FIELDS = ("rac", "mca")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _success(row: dict[str, str], field: str) -> bool:
    return row[f"countercue_{field}_reversal_success"].lower() == "true"


def _all_observations(round_number: int) -> list[dict[str, str]]:
    observations = [
        row for row in _rows(RESULTS / SOURCE_NAME / "task_summary.csv")
        if row["assay"] == "countercue_threshold_sweep"
    ]
    for prior in range(1, round_number):
        path = RESULTS / f"{REFINEMENT_PREFIX}{prior}" / "task_summary.csv"
        if not path.exists():
            raise RuntimeError(f"missing collected refinement round: {path}")
        observations.extend(_rows(path))
    return observations


def _candidate_durations(observations: list[dict[str, str]]) -> dict[tuple[str, float], set[float]]:
    grouped: dict[tuple[str, float], list[dict[str, str]]] = {}
    for row in observations:
        key = (row["coupling_variant"], float(row["length_um"]))
        grouped.setdefault(key, []).append(row)

    candidates: dict[tuple[str, float], set[float]] = {}
    for key, rows in grouped.items():
        for field in FIELDS:
            values = sorted((float(row["countercue_duration_s"]), _success(row, field)) for row in rows)
            successes = [duration for duration, passed in values if passed]
            failures = [duration for duration, passed in values if not passed]
            if not successes:
                continue
            upper = min(successes)
            lower_values = [duration for duration in failures if duration < upper]
            if not lower_values:
                continue
            lower = max(lower_values)
            if upper - lower <= RESOLUTION_S:
                continue
            midpoint = round(((lower + upper) / 2.0) / RESOLUTION_S) * RESOLUTION_S
            midpoint = min(upper - RESOLUTION_S, max(lower + RESOLUTION_S, midpoint))
            candidates.setdefault(key, set()).add(midpoint)
    return candidates


def _manifest(round_number: int) -> tuple[dict, dict[tuple[str, float], set[float]]]:
    source = json.loads(SOURCE_MANIFEST.read_text())
    templates = {}
    for task in source["tasks"]:
        if task["assay"] == "countercue_threshold_sweep":
            templates.setdefault((task["coupling_variant"], float(task["simulation"]["length"])), task)

    candidates = _candidate_durations(_all_observations(round_number))
    tasks = []
    for key in sorted(candidates):
        for duration in sorted(candidates[key]):
            task = copy.deepcopy(templates[key])
            for field in ("task_id", "task_index", "task_digest"):
                task.pop(field, None)
            start = float(task["opposing_cue_start_s"])
            task.update(
                assay="countercue_threshold_refinement",
                countercue_duration_s=duration,
                opposing_cue_duration_s=duration,
                implementation_revision="julia_activation_defaults_threshold_bisection_v1",
                source_campaign=SOURCE_NAME,
                source_assay="countercue_threshold_sweep",
                save_trajectory=False,
            )
            task["protocol"]["pulses"][1]["stop"] = start + duration
            task["polarity_reversal_audit"]["opposite_stop"] = start + duration
            task["simulation"]["end_time"] = start + duration + 10_000.0
            tasks.append(task)

    name = f"{REFINEMENT_PREFIX}{round_number}"
    manifest = make_manifest(
        name,
        tasks,
        description=(
            f"Adaptive counter-cue threshold bisection round {round_number}; "
            f"target duration resolution {RESOLUTION_S:g} s."
        ),
        metadata={
            "round": round_number,
            "resolution_s": RESOLUTION_S,
            "source_campaign": SOURCE_NAME,
            "method": "monotone bisection of independently scored Rac and MCA success brackets",
        },
    )
    return manifest, candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("round", type=int)
    arguments = parser.parse_args()
    if arguments.round < 1:
        raise SystemExit("round must be positive")
    manifest, candidates = _manifest(arguments.round)
    destination = ROOT / "campaigns" / f"{REFINEMENT_PREFIX}{arguments.round}.json"
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks\t{len(candidates)} conditions")


if __name__ == "__main__":
    main()
