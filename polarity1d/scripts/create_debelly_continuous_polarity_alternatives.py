#!/usr/bin/env python3
"""Create the standalone continuous-polarity alternatives requested for Figure 2."""

from __future__ import annotations

import copy
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.campaign import load_manifest, make_manifest, write_manifest  # noqa: E402


NAME = "debelly_2026_continuous_polarity_alternatives_v1"
LENGTHS = tuple(float(value) for value in range(1, 40, 2))
VARIANTS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
BETA_SWEEP = (
    ("beta_0x", 0.0, 0.0),
    ("beta_0p5x", 0.5, 0.27),
    ("beta_1x", 1.0, 0.54),
    ("beta_1p5x", 1.5, 0.81),
    ("beta_2x", 2.0, 1.08),
)
OVERRIDES = {
    "full": {},
    "beta_zero": {"rho_mechanical_feedback": 0.0},
    "alpha_zero": {"rac_mechanical_feedback": 0.0},
    "local_inhibition": {
        "protrusion_speed": 0.0,
        "rac_mechanical_feedback": 0.0,
        "rho_mechanical_feedback": 0.0,
        "active_stress": 0.0,
    },
}


def _single_pulse_task(
    length: float,
    variant: str,
    *,
    duration: float = 160.0,
    assay: str = "single_stimulus_length_test",
) -> dict:
    cue_start = 200.0
    cue_stop = cue_start + duration
    return {
        "tier": "debelly_2026_continuous_polarity_alternatives",
        "implementation_revision": "single_stimulus_length_test_v1",
        "model": "debelly_2026",
        "model_label": "Reflection-symmetric extension of the 2026 de Belly mechanochemical model",
        "assay": assay,
        "coupling_variant": variant,
        "input_peak": 2.0,
        "model_parameters": {
            "assay_mode": "reflection_symmetric_extension",
            "initialization_mode": "algebraic_fixed_point",
            "input_sampling_mode": "cell_centred",
            "rac_basal_activation": 1.8,
            "rho_basal_activation": 1.8,
            "rac_mechanical_feedback": 0.54,
            "rho_mechanical_feedback": 0.54,
            **OVERRIDES[variant],
        },
        "protocol": {
            "name": "localized_pulse",
            "amplitude": 3.6,
            "start": cue_start,
            "stop": cue_stop,
            "taper": duration / 2.0,
            "width_fraction": 0.2,
            "sides": ["left"],
            "asymmetry": 0.0,
        },
        "simulation": {
            "length": length,
            "cells": max(5, round(length / 0.2)),
            "end_time": cue_stop + 6_000.0,
            "output_interval": 2.0,
            "step": 0.5,
        },
        "metrics": {
            "analysis_start": cue_start,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.25,
            "sustained_duration": 10.0,
        },
        "single_cue_kinetics": {
            "cue_start": cue_start,
            "cue_stop": cue_stop,
            "cue_width_fraction": 0.2,
            "cue_amplitude": 2.0,
            "orientation_threshold": 0.25,
            "polarization_sustain_s": 10.0,
            "depolarization_sustain_s": 50.0,
            "cue_side": "left",
        },
        "field_audit": True,
        "save_trajectory": True,
        "replicate": 0,
    }


def _brief_conditioning_reversal_task(length: float, interval: float) -> dict:
    conditioning_start = 200.0
    conditioning_stop = conditioning_start + 40.0
    counter_start = conditioning_stop + interval
    counter_stop = counter_start + 40.0
    return {
        "tier": "debelly_2026_continuous_polarity_alternatives",
        "implementation_revision": "brief_conditioning_interval_reversal_v1",
        "model": "debelly_2026",
        "model_label": "Reflection-symmetric extension of the 2026 de Belly mechanochemical model",
        "assay": f"brief_conditioning_reversal_interval_{interval:g}s",
        "coupling_variant": "full",
        "input_peak": 8.0,
        "model_parameters": {
            "assay_mode": "reflection_symmetric_extension",
            "initialization_mode": "algebraic_fixed_point",
            "input_sampling_mode": "cell_centred",
            "rac_basal_activation": 1.8,
            "rho_basal_activation": 1.8,
            "rac_mechanical_feedback": 0.54,
            "rho_mechanical_feedback": 0.54,
        },
        "protocol": {
            "name": "finite_pulse_sequence",
            "pulses": [
                {
                    "amplitude": 3.6,
                    "start": conditioning_start,
                    "stop": conditioning_stop,
                    "taper": 0.0,
                    "width_fraction": 0.4,
                    "sides": ["left"],
                },
                {
                    "amplitude": 14.4,
                    "start": counter_start,
                    "stop": counter_stop,
                    "taper": 20.0,
                    "width_fraction": 0.4,
                    "sides": ["right"],
                },
            ],
        },
        "simulation": {
            "length": length,
            "cells": max(5, round(length / 0.2)),
            "end_time": counter_stop + 6_000.0,
            "output_interval": 2.0,
            "step": 0.5,
        },
        "metrics": {
            "analysis_start": conditioning_start,
            "polarization_threshold": 0.1,
            "orientation_threshold": 0.05,
            "sustained_duration": 10.0,
        },
        "field_audit": True,
        "save_trajectory": True,
        "replicate": 0,
        "conditioning_duration_s": 40.0,
        "counterstimulus_duration_s": 40.0,
        "interstimulus_interval_s": interval,
    }


def _copied_length_tasks() -> list[dict]:
    sources = (
        (
            ROOT / "campaigns/debelly_2026_responsive_figure2_v2.json",
            "responsive_competing_fronts",
            "figure2_B_continuous",
        ),
        (
            ROOT / "campaigns/debelly_2026_conditioned_reversal_candidate_v1.json",
            "conditioned_finite_reversal_candidate",
            "figure2_D_continuous",
        ),
    )
    output = []
    for path, source_assay, assay in sources:
        manifest = load_manifest(path)
        selected = [
            task for task in manifest["tasks"]
            if task["assay"] == source_assay
            and task["coupling_variant"] in VARIANTS
        ]
        if len(selected) != len(VARIANTS) * len(LENGTHS):
            raise RuntimeError(f"expected 80 length/arm tasks in {path}, found {len(selected)}")
        for source in selected:
            task = copy.deepcopy(source)
            for key in ("task_id", "task_index", "task_digest"):
                task.pop(key, None)
            task["tier"] = "debelly_2026_continuous_polarity_alternatives"
            task["assay"] = assay
            task["source_campaign"] = manifest["campaign_name"]
            task["source_assay"] = source_assay
            task["save_trajectory"] = True
            output.append(task)
    return output


def _beta_sweep_tasks() -> list[dict]:
    """Copy the full-model B protocol while varying only tension feedback beta."""
    path = ROOT / "campaigns/debelly_2026_responsive_figure2_v2.json"
    manifest = load_manifest(path)
    selected = sorted(
        (
            task for task in manifest["tasks"]
            if task["assay"] == "responsive_competing_fronts"
            and task["coupling_variant"] == "full"
        ),
        key=lambda task: float(task["simulation"]["length"]),
    )
    if len(selected) != len(LENGTHS):
        raise RuntimeError(f"expected 20 full-model B tasks in {path}, found {len(selected)}")
    output = []
    for variant, multiplier, beta in BETA_SWEEP:
        for source in selected:
            task = copy.deepcopy(source)
            for key in ("task_id", "task_index", "task_digest"):
                task.pop(key, None)
            task["tier"] = "debelly_2026_continuous_polarity_alternatives"
            task["implementation_revision"] = "figure2_B_beta_sweep_v1"
            task["assay"] = "figure2_B_beta_sweep"
            task["coupling_variant"] = variant
            task["model_parameters"]["rho_mechanical_feedback"] = beta
            task["beta_multiplier"] = multiplier
            task["beta_tension_feedback"] = beta
            task["source_campaign"] = manifest["campaign_name"]
            task["source_assay"] = "responsive_competing_fronts"
            task["save_trajectory"] = True
            output.append(task)
    return output


def build_manifest() -> dict:
    tasks = [
        _single_pulse_task(length, variant)
        for length in LENGTHS
        for variant in VARIANTS
    ]
    tasks.extend(_copied_length_tasks())
    for duration in (40.0, 80.0):
        assay = f"single_stimulus_length_test_{duration:g}s"
        tasks.extend(
            _single_pulse_task(length, "full", duration=duration, assay=assay)
            for length in LENGTHS
        )
    for interval in (40.0, 80.0, 160.0):
        tasks.extend(
            _brief_conditioning_reversal_task(length, interval)
            for length in LENGTHS
        )
    for variant in VARIANTS[1:]:
        for duration in (40.0, 80.0):
            assay = f"single_stimulus_length_test_{duration:g}s"
            tasks.extend(
                _single_pulse_task(length, variant, duration=duration, assay=assay)
                for length in LENGTHS
            )
    tasks.extend(_beta_sweep_tasks())
    return make_manifest(
        NAME,
        tasks,
        description=(
            "Standalone Figure 2 alternatives: one left stimulus followed by withdrawal "
            "across length, plus raw continuous polarity trajectories across length for panels B and D, "
            "including a five-level full-model beta sweep for B."
        ),
        metadata={
            "lengths_um": list(LENGTHS),
            "variants": list(VARIANTS),
            "figure2_B_beta_sweep": [
                {"variant": variant, "multiplier": multiplier, "beta": beta}
                for variant, multiplier, beta in BETA_SWEEP
            ],
            "polarity_score": (
                "Relative half-mass position x50/L: the interpolated point along the cell at which "
                "50% of the relevant field's content lies on either side"
            ),
            "single_stimulus_protocol": (
                "unstimulated 0-200 s; left outer-20% cue at 2*alpha0 for 40, 80, or "
                "160 s, full for the first half then cosine tapered for the second half; "
                "observe 6000 s after withdrawal"
            ),
            "threshold_latencies": (
                "polarize at cue-aligned Rac score >=0.25 for 10 s; depolarize at absolute "
                "Rac score <0.25 for 50 s; 0.25 matches the strong-orientation gate in Figure 2"
            ),
            "brief_conditioning_reversal_protocol": (
                "unstimulated 0-200 s; left outer-40% cue at 2*alpha0 for 40 s; "
                "cue-free interval of 40, 80, or 160 s; right outer-40% cue at "
                "8*alpha0 for 40 s with a 20-s cosine taper; observe 6000 s after withdrawal"
            ),
        },
    )


def main() -> None:
    destination = ROOT / "campaigns" / f"{NAME}.json"
    manifest = build_manifest()
    write_manifest(destination, manifest)
    print(f"{destination}\t{manifest['task_count']} tasks")


if __name__ == "__main__":
    main()
