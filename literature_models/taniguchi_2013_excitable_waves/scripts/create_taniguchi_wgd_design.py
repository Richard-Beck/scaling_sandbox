#!/usr/bin/env python3
"""Create a predeclared paired-seed Taniguchi WGD design."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("pilot", "production"), required=True)
    parser.add_argument(
        "--replicates", type=int, required=True,
        help="number of paired seeds declared before any run is launched",
    )
    parser.add_argument("--first-seed", type=int)
    parser.add_argument("--directory", type=Path)
    parser.add_argument("--output-root")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.replicates < 2:
        raise SystemExit("--replicates must be at least 2")
    first_seed = args.first_seed
    if first_seed is None:
        first_seed = 1 if args.stage == "pilot" else 1001
    if first_seed < 0:
        raise SystemExit("--first-seed must be nonnegative")
    directory = args.directory or (
        ROOT / "experiments" / "wgd_v1" / "designs" /
        f"taniguchi_{args.stage}"
    )
    if not directory.is_absolute():
        directory = ROOT / directory
    output_root = args.output_root or f"results/wgd_v1/taniguchi_{args.stage}"
    directory.mkdir(parents=True, exist_ok=True)
    config_paths: list[Path] = []
    for seed in range(first_seed, first_seed + args.replicates):
        for volume in (1, 2):
            label = f"{volume}x_seed{seed}"
            config = {
                "protocol_version": "wgd_v1",
                "model": "taniguchi_phase_field",
                "volume_ratio": volume,
                "scaling_condition": "wgd_homeostatic",
                "stimulus_condition": "relative_default",
                "seed": seed,
                "output_path": f"{output_root}/{label}",
                "model_config": {
                    "backend": "fast",
                    "checkpoint_interval": 2.0,
                    "parameter_set": "4F",
                },
            }
            path = directory / f"taniguchi_phase_field_{label}.json"
            path.write_text(
                json.dumps(config, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            config_paths.append(path)

    def manifest_path(path: Path) -> Path:
        try:
            return path.relative_to(ROOT)
        except ValueError:
            return path

    manifest = directory / "configs.txt"
    manifest.write_text(
        "".join(f"{manifest_path(path)}\n" for path in config_paths),
        encoding="utf-8",
    )
    seeds = list(range(first_seed, first_seed + args.replicates))
    design = {
        "model": "taniguchi_phase_field",
        "stage": args.stage,
        "purpose": (
            "paired-seed variability pilot; not the final inferential design"
            if args.stage == "pilot" else
            "frozen paired-seed WGD production design"
        ),
        "replicates": args.replicates,
        "first_seed": first_seed,
        "seeds": seeds,
        "backend": "fast",
        "parameter_set": "4F",
        "raw_format": "taniguchi_wgd_compact_v1",
        "checkpoint_interval": 2.0,
        "volume_ratios": [1, 2],
    }
    if args.stage == "production":
        design["replicate_count_basis"] = (
            "Eight-pair pilot paired SD 0.08459 for the fraction of frames "
            "with an active domain; n=30.54 for a normal-approximation 95% "
            "half-width of 0.03, rounded up to 32 paired seeds."
        )
    (directory / "design.json").write_text(
        json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
