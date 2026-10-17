#!/usr/bin/env python3
"""Plot the constant-per-cell Taniguchi firing sensitivity analysis."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from plot_taniguchi_geometry_campaign import (  # noqa: E402
    LABELS,
    groups,
    migration_dynamics_figure,
    migration_figure,
    paired_arrays,
    write_csv,
    write_track_csv,
)


RULE_COLOURS = {"per_area": "#4C78A8", "per_cell_reference": "#E45756"}
RULE_LABELS = {
    "per_area": "Primary: rate per area",
    "per_cell_reference": "Sensitivity: rate per cell",
}


def load_rows(root: Path, expected: int, *, deformable_seeds: set[int] | None = None) -> list[dict]:
    rows = []
    for path in sorted(root.glob("*/summary.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        if (deformable_seeds is not None
                and row.get("arm") == "deformable_spontaneous"
                and int(row["seed"]) not in deformable_seeds):
            continue
        row["run_dir"] = str(path.parent)
        rows.append(row)
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} completed runs under {root}, found {len(rows)}")
    return rows


def interval(values: np.ndarray) -> tuple[float, float]:
    mean = float(np.mean(values))
    half = 1.96 * float(np.std(values, ddof=1)) / math.sqrt(values.size)
    return mean, half


def interaction_figure(
    primary: list[dict], sensitivity: list[dict], arm: str,
    specs: tuple[tuple[str, str], ...], output: Path, title: str,
) -> dict:
    datasets = {"per_area": primary, "per_cell_reference": sensitivity}
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), constrained_layout=True)
    summaries: dict[str, dict] = {}
    for ax, (metric, label) in zip(axes.ravel(), specs):
        summaries[metric] = {}
        for rule, rows in datasets.items():
            group = groups(rows, arm)
            means, halves = [], []
            for volume in (1, 2):
                values = np.asarray([float(row[metric]) for row in group[volume]
                                     if row.get(metric) is not None])
                mean, half = interval(values)
                means.append(mean)
                halves.append(half)
            first, second = paired_arrays(group, metric)
            difference = second - first
            summaries[metric][rule] = {
                "mean_1x": means[0], "mean_2x": means[1],
                "paired_difference": float(np.mean(difference)),
                "paired_difference_ci95_half_width": (
                    1.96 * float(np.std(difference, ddof=1)) /
                    math.sqrt(difference.size)
                ),
            }
            ax.errorbar([1, 2], means, yerr=halves, marker="o", linewidth=2,
                        capsize=4, color=RULE_COLOURS[rule], label=RULE_LABELS[rule])
        ax.set_xticks([1, 2], ["1×", "2×"])
        ax.set_ylabel(label)
        ax.grid(color="#dddddd", lw=0.5)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.supxlabel("Volume")
    fig.suptitle(title, fontsize=14)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return summaries


def paired_effect_figure(
    primary: list[dict], sensitivity: list[dict], arm: str,
    specs: tuple[tuple[str, str], ...], output: Path, title: str,
) -> None:
    datasets = {"per_area": primary, "per_cell_reference": sensitivity}
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), constrained_layout=True)
    for ax, (metric, label) in zip(axes.ravel(), specs):
        for index, rule in enumerate(("per_area", "per_cell_reference"), start=1):
            first, second = paired_arrays(groups(datasets[rule], arm), metric)
            difference = second - first
            mean, half = interval(difference)
            ax.errorbar(index, mean, yerr=half, marker="o", markersize=6,
                        capsize=5, color=RULE_COLOURS[rule])
        ax.axhline(0, color="black", lw=0.8, ls="--")
        ax.set_xticks([1, 2], ["Per area", "Per cell"])
        ax.set_ylabel(f"2× − 1×\n{label}")
        ax.grid(axis="y", color="#dddddd", lw=0.5)
    fig.suptitle(title, fontsize=14)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def verify_paired_schedules(rows: list[dict]) -> dict:
    checks = {}
    for arm in ("fixed_spontaneous", "deformable_spontaneous"):
        group = groups(rows, arm)
        matching = 0
        for left, right in zip(group[1], group[2], strict=True):
            if int(left["seed"]) != int(right["seed"]):
                raise RuntimeError("paired schedule seed order mismatch")
            with np.load(Path(left["run_dir"]) / "events_and_tracks.npz") as first, \
                 np.load(Path(right["run_dir"]) / "events_and_tracks.npz") as second:
                if (np.array_equal(first["event_time"], second["event_time"])
                        and np.array_equal(first["event_amplitude"],
                                           second["event_amplitude"])):
                    matching += 1
        if matching != len(group[1]):
            raise RuntimeError(
                f"only {matching}/{len(group[1])} {arm} pairs have matched schedules"
            )
        checks[arm] = {"matched_pairs": matching}
    return checks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primary-root", type=Path,
        default=ROOT / "results/wgd_v1/taniguchi_geometry_v1",
    )
    parser.add_argument(
        "--sensitivity-root", type=Path,
        default=ROOT / "results/wgd_v1/taniguchi_firing_sensitivity_v1",
    )
    args = parser.parse_args()
    primary_root = args.primary_root if args.primary_root.is_absolute() else ROOT / args.primary_root
    sensitivity_root = (args.sensitivity_root if args.sensitivity_root.is_absolute()
                        else ROOT / args.sensitivity_root)
    primary = load_rows(primary_root, 130)
    sensitivity_design_path = (ROOT / "experiments" / "wgd_v1" / "designs" /
                               "taniguchi_firing_sensitivity_v1" / "design.json")
    sensitivity_design = json.loads(
        sensitivity_design_path.read_text(encoding="utf-8")
    )
    sensitivity = load_rows(
        sensitivity_root, 128,
        deformable_seeds=set(sensitivity_design["deformable_paired_seeds"]),
    )
    figure_dir = sensitivity_root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(figure_dir / "run_summaries.csv", sensitivity)
    write_track_csv(figure_dir / "wave_tracks.csv", sensitivity)
    schedule_checks = verify_paired_schedules(sensitivity)

    fixed_specs = (
        ("event_rate_per_cell_time", "Firing attempts / cell / time"),
        ("successful_wave_initiation_rate_per_cell_time", "Wave initiations / cell / time"),
        ("fraction_time_active", "Fraction of time active"),
        ("mean_simultaneous_active_patches", "Mean simultaneous patches"),
        ("median_wave_lifetime", "Median patch lifetime"),
        ("mean_wave_travel_per_diameter", "Travel / cell diameter"),
    )
    deformable_specs = (
        ("event_rate_per_cell_time", "Firing attempts / cell / time"),
        ("fraction_time_active", "Fraction of time active"),
        ("net_displacement", "Net displacement"),
        ("path_length", "Path length"),
        ("mean_protrusion_count", "Mean protrusion count"),
        ("pip3_to_future_boundary_extension_correlation", "PIP3 → extension corr."),
    )
    fixed_summary = interaction_figure(
        primary, sensitivity, "fixed_spontaneous", fixed_specs,
        figure_dir / "01_fixed_firing_rule_interaction.png",
        "Fixed domains — sensitivity to the firing-rate scaling rule",
    )
    paired_effect_figure(
        primary, sensitivity, "fixed_spontaneous", fixed_specs,
        figure_dir / "02_fixed_firing_rule_paired_effects.png",
        "Fixed domains — paired geometry effect under each firing rule",
    )
    deformable_summary = interaction_figure(
        primary, sensitivity, "deformable_spontaneous", deformable_specs,
        figure_dir / "03_deformable_firing_rule_interaction.png",
        "Deformable cells — sensitivity to the firing-rate scaling rule",
    )
    paired_effect_figure(
        primary, sensitivity, "deformable_spontaneous", deformable_specs,
        figure_dir / "04_deformable_firing_rule_paired_effects.png",
        "Deformable cells — paired geometry effect under each firing rule",
    )
    migration_figure(sensitivity, figure_dir)
    migration_dynamics_figure(sensitivity, figure_dir)

    manifest = {
        "campaign": "taniguchi_firing_sensitivity_v1",
        "primary_campaign": "taniguchi_geometry_v1",
        "figures": sorted(path.name for path in figure_dir.glob("*.png")),
        "derived_tables": ["run_summaries.csv", "wave_tracks.csv"],
        "schedule_checks": schedule_checks,
        "fixed_summaries": fixed_summary,
        "deformable_summaries": deformable_summary,
        "migration_units": "native model length units; no initial-radius normalization",
        "turning_angle_definition": (
            "Angle between successive non-overlapping 2-model-time-unit centroid displacements"
        ),
    }
    (figure_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(figure_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
