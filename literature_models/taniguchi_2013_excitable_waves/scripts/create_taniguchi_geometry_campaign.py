#!/usr/bin/env python3
"""Freeze the three-arm Taniguchi pure-geometry production design."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=32)
    parser.add_argument("--first-seed", type=int, default=3001)
    parser.add_argument("--end-time", type=float, default=400.0)
    parser.add_argument("--controlled-end-time", type=float, default=80.0)
    parser.add_argument(
        "--directory", type=Path,
        default=ROOT / "experiments/wgd_v1/designs/taniguchi_geometry_v1",
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=Path("results/wgd_v1/taniguchi_geometry_v1"),
    )
    args = parser.parse_args()
    if args.replicates < 30 or args.replicates > 100:
        raise SystemExit("--replicates must be between 30 and 100")
    directory = args.directory if args.directory.is_absolute() else ROOT / args.directory
    directory.mkdir(parents=True, exist_ok=True)
    common = {
        "parameter_set": "4F",
        "dt": 5e-4,
        "dx": 0.1,
        "output_interval": 0.1,
        "activity_threshold": 0.5,
        "minimum_domain_area": 0.2,
        "radial_bins": 72,
        "snapshot_bins": 8,
        "baseline_box_half_width": 10.0,
    }
    configs: list[Path] = []

    def write_config(name: str, payload: dict) -> None:
        path = directory / f"{name}.json"
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
        configs.append(path)

    for arm in ("fixed_spontaneous", "deformable_spontaneous"):
        for seed in range(args.first_seed, args.first_seed + args.replicates):
            for volume in (1, 2):
                label = f"{arm}_{volume}x_seed{seed}"
                write_config(label, {
                    "campaign": "taniguchi_geometry_v1",
                    "arm": arm,
                    "volume_ratio": volume,
                    "seed": seed,
                    "output_path": str(args.output_root / label),
                    "model_config": {**common, "end_time": args.end_time},
                })
    for volume in (1, 2):
        label = f"fixed_controlled_{volume}x"
        write_config(label, {
            "campaign": "taniguchi_geometry_v1",
            "arm": "fixed_controlled",
            "volume_ratio": volume,
            "seed": None,
            "output_path": str(args.output_root / label),
            "model_config": {
                **common,
                "end_time": args.controlled_end_time,
                "output_interval": 0.05,
                "trigger_relative_x": -0.5,
                "trigger_relative_y": 0.0,
                "trigger_amplitude": 5.7,
                "trigger_width": 0.07,
                "snapshot_times": [0.0, 1.0, 2.0, 3.5, 5.0, 7.0, 9.0, 12.0],
            },
        })
    relative = []
    for path in configs:
        try:
            relative.append(path.relative_to(ROOT))
        except ValueError:
            relative.append(path)
    (directory / "configs.txt").write_text(
        "".join(f"{path}\n" for path in relative), encoding="utf-8"
    )
    design = {
        "campaign": "taniguchi_geometry_v1",
        "purpose": "pure geometric scaling across fixed spontaneous, controlled-wave, and deformable arms",
        "parameter_set": "4F in all arms",
        "replicates_per_stochastic_condition": args.replicates,
        "paired_seed_range": [args.first_seed, args.first_seed + args.replicates - 1],
        "volume_ratios": [1, 2],
        "linear_scale_2x": 2 ** (1 / 3),
        "area_scale_2x": 2 ** (2 / 3),
        "spontaneous_end_time": args.end_time,
        "controlled_end_time": args.controlled_end_time,
        "stochastic_configuration_count": 4 * args.replicates,
        "controlled_configuration_count": 2,
        "total_configuration_count": len(configs),
        "output_format": "taniguchi_geometry_compact_v1",
        "fixed_boundary_definition": (
            "Hold the initialized diffuse phi circle fixed while advancing the identical "
            "Fig. 4F conservative U/V equations, including interface loss."
        ),
        "event_process": "fixed intensity per physical area-time; fixed physical width and amplitude distribution",
        "controlled_trigger": {
            "relative_position": [-0.5, 0.0],
            "physical_amplitude": 5.7,
            "physical_width": 0.07,
            "stochasticity": False,
        },
    }
    (directory / "design.json").write_text(
        json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(directory / "configs.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
