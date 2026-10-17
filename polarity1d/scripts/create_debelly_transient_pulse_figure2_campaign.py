#!/usr/bin/env python3
"""Create the transient-pulse backup campaign for de Belly Figure 2."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import make_manifest, write_manifest  # noqa: E402


CAMPAIGN_NAME = "debelly_2026_transient_pulse_figure2_v1"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
ALPHA0 = 3.5
BETA0 = 4.0
ALPHA = 1.05
BETA = 1.20
REVERSAL_AMPLITUDE_RATIOS = (1.0, 2.0, 4.0, 8.0)
REVERSAL_DURATIONS = (80.0, 160.0, 320.0, 640.0, 1_280.0)
REPLICATE_MULTIPLIERS = (0.90, 0.95, 1.00, 1.05, 1.10)
OPPOSITE_WIDTH = 0.4
REPEATED_RATIO = 4.0
REPEATED_DURATION = 320.0


def _cells(length: float) -> int:
    return max(5, round(length / 0.2))


def _parameters(local: bool = False) -> dict:
    values = {
        "assay_mode": "reflection_symmetric_extension",
        "initialization_mode": "algebraic_fixed_point",
        "input_sampling_mode": "cell_centred",
        "rac_basal_activation": ALPHA0,
        "rho_basal_activation": BETA0,
        "rac_mechanical_feedback": ALPHA,
        "rho_mechanical_feedback": BETA,
    }
    if local:
        values.update(
            protrusion_speed=0.0,
            rac_mechanical_feedback=0.0,
            rho_mechanical_feedback=0.0,
            active_stress=0.0,
        )
    return values


def _pulse(
    *, side: str, start: float, duration: float, amplitude: float,
    width: float = 0.1,
) -> dict:
    return {
        "amplitude": amplitude,
        "start": start,
        "stop": start + duration,
        "taper": duration / 2.0,
        "width_fraction": width,
        "sides": [side],
    }


def _base(assay: str, length: float, *, local: bool = False) -> dict:
    return {
        "tier": "debelly_transient_pulse_backup",
        "implementation_revision": "selected_regime_transient_pulses_v1",
        "model": "debelly_2026",
        "model_label": (
            "Reflection-symmetric extension of the 2026 de Belly "
            "mechanochemical model"
        ),
        "assay": assay,
        "coupling_variant": "local_inhibition" if local else "full",
        "model_parameters": _parameters(local),
        "simulation": {
            "length": length,
            "cells": _cells(length),
            "step": 0.5,
            "output_interval": 5.0,
        },
        "metrics": {
            "analysis_start": 40.0,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 100.0,
        },
        "mechanochemical_metrics": {"stimulus_start": 0.0},
        "field_audit": True,
        "replicate": 0,
        "save_trajectory": False,
    }


def campaign_tasks() -> list[dict]:
    tasks: list[dict] = []
    # Induction, 5,000-s cue-free stabilization, one opposite pulse, then
    # 5,000 s of cue-free observation.
    opposite_start = 5_040.0
    for length in LENGTHS:
        for ratio in REVERSAL_AMPLITUDE_RATIOS:
            for duration in REVERSAL_DURATIONS:
                events = [
                    {"side": "left", "start": 0.0, "stop": 40.0},
                    {
                        "side": "right",
                        "start": opposite_start,
                        "stop": opposite_start + duration,
                    },
                ]
                task = _base("transient_pulse_reversal", length)
                task.update(
                    input_peak=ratio,
                    pulse_amplitude_ratio=ratio,
                    pulse_duration_s=duration,
                    protocol={
                        "name": "finite_pulse_sequence",
                        "pulses": [
                            _pulse(side="left", start=0.0, duration=40.0, amplitude=ALPHA0),
                            _pulse(
                                side="right", start=opposite_start,
                                duration=duration, amplitude=ratio * ALPHA0,
                                width=OPPOSITE_WIDTH,
                            ),
                        ],
                    },
                    pulse_switching={
                        "events": events,
                        "orientation_threshold": 0.05,
                        "sustain_s": 100.0,
                        "minimum_persistence_s": 1_000.0,
                    },
                )
                task["simulation"]["end_time"] = opposite_start + duration + 5_000.0
                task["save_trajectory"] = bool(
                    length == 15.0 and ratio == REPEATED_RATIO
                    and duration == REPEATED_DURATION
                )
                tasks.append(task)

        # One mechanics-removed reference at the central pulse condition.
        duration = REPEATED_DURATION
        events = [
            {"side": "left", "start": 0.0, "stop": 40.0},
            {"side": "right", "start": opposite_start, "stop": opposite_start + duration},
        ]
        local = _base("transient_pulse_reversal", length, local=True)
        local.update(
            input_peak=REPEATED_RATIO,
            pulse_amplitude_ratio=REPEATED_RATIO,
            pulse_duration_s=duration,
            protocol={
                "name": "finite_pulse_sequence",
                "pulses": [
                    _pulse(side="left", start=0.0, duration=40.0, amplitude=ALPHA0),
                    _pulse(
                        side="right", start=opposite_start, duration=duration,
                        amplitude=REPEATED_RATIO * ALPHA0, width=OPPOSITE_WIDTH,
                    ),
                ],
            },
            pulse_switching={
                "events": events,
                "orientation_threshold": 0.05,
                "sustain_s": 100.0,
                "minimum_persistence_s": 1_000.0,
            },
        )
        local["simulation"]["end_time"] = opposite_start + duration + 5_000.0
        local["save_trajectory"] = length == 15.0
        tasks.append(local)

    # Five alternating pulses; each is followed by approximately 5,000 s
    # without input. Four post-induction events are scored as switches.
    starts = (0.0, 5_040.0, 10_360.0, 15_680.0, 21_000.0)
    sides = ("left", "right", "left", "right", "left")
    for length in LENGTHS:
        for replicate, multiplier in enumerate(REPLICATE_MULTIPLIERS):
            pulses = [_pulse(side="left", start=0.0, duration=40.0, amplitude=ALPHA0)]
            pulses.extend(
                _pulse(
                    side=side, start=start, duration=REPEATED_DURATION,
                    amplitude=REPEATED_RATIO * ALPHA0 * multiplier,
                    width=OPPOSITE_WIDTH,
                )
                for side, start in zip(sides[1:], starts[1:])
            )
            events = [{"side": "left", "start": 0.0, "stop": 40.0}]
            events.extend(
                {
                    "side": side, "start": start,
                    "stop": start + REPEATED_DURATION,
                }
                for side, start in zip(sides[1:], starts[1:])
            )
            task = _base("repeated_transient_switching", length)
            task.update(
                replicate=replicate,
                perturbation_multiplier=multiplier,
                input_peak=REPEATED_RATIO * multiplier,
                protocol={"name": "finite_pulse_sequence", "pulses": pulses},
                pulse_switching={
                    "events": events,
                    "orientation_threshold": 0.05,
                    "sustain_s": 100.0,
                    "minimum_persistence_s": 1_000.0,
                },
            )
            task["simulation"]["end_time"] = 26_320.0
            task["save_trajectory"] = length == 15.0 and replicate == 2
            tasks.append(task)
    return tasks


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output", nargs="?", type=Path,
        default=ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json",
    )
    args = parser.parse_args(argv)
    manifest = make_manifest(
        CAMPAIGN_NAME,
        campaign_tasks(),
        description=(
            "Transient-pulse backup for Figure 2: autonomous reversal after "
            "withdrawal and repeated finite-pulse switching."
        ),
        metadata={
            "selected_parameters": {
                "alpha0": ALPHA0, "beta0": BETA0,
                "alpha": ALPHA, "beta": BETA,
            },
            "selection_independent_of_length_sweep": True,
            "lengths_um": list(LENGTHS),
            "reversal_amplitude_ratios": list(REVERSAL_AMPLITUDE_RATIOS),
            "reversal_durations_s": list(REVERSAL_DURATIONS),
            "repeated_switching_multipliers": list(REPLICATE_MULTIPLIERS),
            "opposite_pulse_width_fraction": OPPOSITE_WIDTH,
            "repeated_switching_nominal_ratio": REPEATED_RATIO,
            "repeated_switching_duration_s": REPEATED_DURATION,
            "switch_definition": (
                "target Rac/Rho orientations first established after pulse "
                "withdrawal for >=100 s and retained cue-free for >=1000 s"
            ),
            "competition_source": "debelly_2026_jilkine_common_v1",
        },
    )
    write_manifest(args.output, manifest)
    print(f"{args.output}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
