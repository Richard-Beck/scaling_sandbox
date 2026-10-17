#!/usr/bin/env python3
"""Freeze the constant-per-cell Taniguchi firing sensitivity design."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=32)
    parser.add_argument("--first-seed", type=int, default=3001)
    parser.add_argument("--end-time", type=float, default=400.0)
    parser.add_argument(
        "--exclude-deformable-seed", type=int, default=3027,
        help="Numerically failed paired seed excluded from the frozen deformable analysis",
    )
    parser.add_argument(
        "--replacement-deformable-seed", type=int, default=3033,
        help="Next unused seed replacing the excluded deformable pair",
    )
    parser.add_argument(
        "--directory", type=Path,
        default=ROOT / "experiments/wgd_v1/designs/taniguchi_firing_sensitivity_v1",
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=Path("results/wgd_v1/taniguchi_firing_sensitivity_v1"),
    )
    args = parser.parse_args()
    if args.replicates < 30 or args.replicates > 100:
        raise SystemExit("--replicates must be between 30 and 100")
    directory = args.directory if args.directory.is_absolute() else ROOT / args.directory
    directory.mkdir(parents=True, exist_ok=True)
    reference_area = 25.0 * np.pi
    common = {
        "parameter_set": "4F",
        "dt": 5e-4,
        "dx": 0.1,
        "end_time": args.end_time,
        "output_interval": 0.1,
        "activity_threshold": 0.5,
        "minimum_domain_area": 0.2,
        "radial_bins": 72,
        "snapshot_bins": 8,
        "baseline_box_half_width": 10.0,
        "firing_rule": "per_cell_reference",
        "reference_event_area": reference_area,
    }
    configs: list[Path] = []
    base_seeds = list(range(args.first_seed, args.first_seed + args.replicates))
    if args.exclude_deformable_seed not in base_seeds:
        raise SystemExit("--exclude-deformable-seed must lie in the base seed range")
    if args.replacement_deformable_seed in base_seeds:
        raise SystemExit("--replacement-deformable-seed must be unused")
    arm_seeds = {
        "fixed_spontaneous": base_seeds,
        "deformable_spontaneous": [
            seed for seed in base_seeds if seed != args.exclude_deformable_seed
        ] + [args.replacement_deformable_seed],
    }
    for arm in ("fixed_spontaneous", "deformable_spontaneous"):
        for seed in arm_seeds[arm]:
            for volume in (1, 2):
                label = f"{arm}_{volume}x_seed{seed}"
                path = directory / f"{label}.json"
                path.write_text(json.dumps({
                    "campaign": "taniguchi_firing_sensitivity_v1",
                    "arm": arm,
                    "volume_ratio": volume,
                    "seed": seed,
                    "output_path": str(args.output_root / label),
                    "model_config": common,
                }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                configs.append(path)

    relative = [path.relative_to(ROOT) for path in configs]
    (directory / "configs.txt").write_text(
        "".join(f"{path}\n" for path in relative), encoding="utf-8"
    )
    design = {
        "campaign": "taniguchi_firing_sensitivity_v1",
        "purpose": (
            "sensitivity analysis holding stochastic firing attempts per cell "
            "constant across 1x and 2x geometry"
        ),
        "parameter_set": "4F",
        "arms": ["fixed_spontaneous", "deformable_spontaneous"],
        "replicates_per_condition": args.replicates,
        "fixed_paired_seed_range": [args.first_seed,
                                    args.first_seed + args.replicates - 1],
        "deformable_paired_seeds": arm_seeds["deformable_spontaneous"],
        "numerical_failure_handling": {
            "excluded_seed": args.exclude_deformable_seed,
            "replacement_seed": args.replacement_deformable_seed,
            "reason": (
                "The 1x deformable run collapsed to zero phase-field area after "
                "an extreme late firing amplitude (55.346 at t=378.6895). The "
                "entire paired seed was excluded and replaced before analysis."
            ),
        },
        "volume_ratios": [1, 2],
        "linear_scale_2x": 2 ** (1 / 3),
        "area_scale_2x": 2 ** (2 / 3),
        "end_time": args.end_time,
        "configuration_count": len(configs),
        "output_format": "taniguchi_geometry_compact_v1",
        "firing_rule": "per_cell_reference",
        "reference_event_area": reference_area,
        "prescribed_event_rate_per_cell_time": 1.0e-3 * reference_area,
        "pairing": (
            "Within each arm and seed, 1x and 2x use the same Poisson-arrival "
            "and amplitude stream and a separate same-seed spatial-placement stream."
        ),
        "unchanged_from_primary": (
            "Fig. 4F chemistry, concentrations, diffusion and reaction constants, "
            "event width and amplitude law, dx, dt, duration, and mechanics"
        ),
        "changed_from_primary": (
            "Firing intensity is per cell rather than per instantaneous physical area."
        ),
    }
    (directory / "design.json").write_text(
        json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(directory / "configs.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
