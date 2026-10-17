#!/usr/bin/env python3
"""Create the 10-min conditioning/reversal campaign for standalone Figure 2 v2."""

from __future__ import annotations

import copy
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


NAME = "debelly_2026_standalone_figure2_v2_reversal"
SOURCE = ROOT / "campaigns/debelly_2026_continuous_polarity_alternatives_v1.json"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
BETA_VARIANTS = ("beta_0x", "beta_1x", "beta_2x")
CONDITIONING_STOP = 600.0
REVERSAL_STOP = 1_200.0
END_TIME = 20_600.0


def campaign_tasks() -> list[dict]:
    """Reuse the approved full-model beta arms while replacing only the protocol."""
    source = load_manifest(SOURCE)
    templates = {
        (task["coupling_variant"], float(task["simulation"]["length"])): task
        for task in source["tasks"]
        if task["assay"] == "figure2_B_beta_sweep"
        and task["coupling_variant"] in BETA_VARIANTS
    }
    expected = {(variant, length) for variant in BETA_VARIANTS for length in LENGTHS}
    if templates.keys() != expected:
        raise RuntimeError("source campaign does not contain the required beta/length grid")

    tasks = []
    for variant in BETA_VARIANTS:
        for length in LENGTHS:
            task = copy.deepcopy(templates[(variant, length)])
            for key in ("task_id", "task_index", "task_digest"):
                task.pop(key, None)
            task.update(
                tier="debelly_2026_standalone_figure2_v2",
                implementation_revision="standalone_figure2_v2_reversal_v1",
                assay="figure2_v2_reversal",
                input_peak=4.0,
                protocol={
                    "name": "finite_pulse_sequence",
                    "pulses": [
                        {
                            "amplitude": 3.6,
                            "start": 0.0,
                            "stop": CONDITIONING_STOP,
                            "taper": 0.0,
                            "width_fraction": 0.2,
                            "sides": ["left"],
                        },
                        {
                            "amplitude": 7.2,
                            "start": CONDITIONING_STOP,
                            "stop": REVERSAL_STOP,
                            "taper": 0.0,
                            "width_fraction": 0.4,
                            "sides": ["right"],
                        },
                    ],
                },
                reversal_design={
                    "conditioning_start_s": 0.0,
                    "conditioning_stop_s": CONDITIONING_STOP,
                    "reversal_start_s": CONDITIONING_STOP,
                    "reversal_stop_s": REVERSAL_STOP,
                    "crossing_position_fraction": 0.5,
                },
                source_campaign=source["campaign_name"],
                source_assay="figure2_B_beta_sweep",
                save_trajectory=True,
            )
            task["simulation"].update(
                end_time=END_TIME,
                output_interval=5.0,
                step=0.5,
            )
            task["metrics"].update(
                analysis_start=CONDITIONING_STOP,
                sustained_duration=50.0,
            )
            tasks.append(task)
    return tasks


def build_manifest() -> dict:
    return make_manifest(
        NAME,
        campaign_tasks(),
        description=(
            "Standalone Figure 2 v2 reversal assay: a 600-s left conditioning pulse "
            "followed immediately by a 600-s right reversal pulse in the intact model "
            "at beta=0x, 1x, and 2x, with a 20,000-s post-reversal-start horizon."
        ),
        metadata={
            "lengths_um": list(LENGTHS),
            "beta_variants": list(BETA_VARIANTS),
            "default_beta": 0.54,
            "conditioning_protocol": (
                "left outer-20% rectangular cue at 2*alpha0 from 0 through 600 s"
            ),
            "reversal_protocol": (
                "right outer-40% rectangular cue at 4*alpha0 from 600 through 1200 s"
            ),
            "observation_horizon": "20,000 s after reversal onset",
            "trajectory_output_interval_s": 5.0,
        },
    )


def main() -> None:
    destination = ROOT / "campaigns" / f"{NAME}.json"
    manifest = build_manifest()
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
