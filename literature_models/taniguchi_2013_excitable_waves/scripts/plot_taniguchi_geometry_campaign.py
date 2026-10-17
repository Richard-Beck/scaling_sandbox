#!/usr/bin/env python3
"""Create the standalone figure and table suite for Taniguchi geometry v1."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
COLOURS = {1: "#2878B5", 2: "#D95319"}
LABELS = {1: "1× volume", 2: "2× volume"}


def load_runs(root: Path) -> list[dict]:
    rows = []
    for path in sorted(root.glob("*/summary.json")):
        row = json.loads(path.read_text(encoding="utf-8"))
        row["run_dir"] = str(path.parent)
        rows.append(row)
    expected = 130
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} completed runs, found {len(rows)}")
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    keys = sorted({key for row in rows for key in row if key not in {"trigger"}})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in keys})


def write_track_csv(path: Path, rows: list[dict]) -> None:
    track_rows: list[dict] = []
    for row in rows:
        with np.load(Path(row["run_dir"]) / "events_and_tracks.npz") as data:
            tracks = json.loads(str(data["track_json"]))
        for track in tracks:
            track_rows.append({
                "arm": row["arm"],
                "volume_ratio": row["volume_ratio"],
                "seed": row.get("seed"),
                **track,
            })
    if not track_rows:
        return
    keys = sorted({key for row in track_rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(track_rows)


def groups(rows: list[dict], arm: str) -> dict[int, list[dict]]:
    return {
        volume: sorted(
            [row for row in rows if row["arm"] == arm
             and int(row["volume_ratio"]) == volume],
            key=lambda row: (-1 if row.get("seed") is None else int(row["seed"])),
        )
        for volume in (1, 2)
    }


def paired_arrays(group: dict[int, list[dict]], metric: str) -> tuple[np.ndarray, np.ndarray]:
    maps = {
        volume: {row.get("seed"): row.get(metric) for row in group[volume]}
        for volume in (1, 2)
    }
    seeds = sorted(set(maps[1]) & set(maps[2]), key=lambda value: -1 if value is None else value)
    first, second = [], []
    for seed in seeds:
        if maps[1][seed] is not None and maps[2][seed] is not None:
            first.append(float(maps[1][seed]))
            second.append(float(maps[2][seed]))
    return np.asarray(first), np.asarray(second)


def paired_axis(ax: plt.Axes, first: np.ndarray, second: np.ndarray,
                ylabel: str, *, normalize_difference: bool = False) -> dict:
    for left, right in zip(first, second):
        ax.plot([1, 2], [left, right], color="#888888", alpha=0.28, lw=0.7)
    ax.scatter(np.ones(first.size), first, color=COLOURS[1], s=18, alpha=0.75, zorder=2)
    ax.scatter(np.full(second.size, 2), second, color=COLOURS[2], s=18, alpha=0.75, zorder=2)
    means = [float(np.mean(first)), float(np.mean(second))]
    ax.plot([1, 2], means, color="black", lw=2.2, marker="o", ms=5, zorder=3)
    ax.set_xticks([1, 2], ["1×", "2×"])
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#dddddd", lw=0.6)
    difference = second - first
    if normalize_difference:
        difference = difference / np.maximum(np.abs(first), 1e-12)
    return {
        "n_pairs": int(first.size),
        "mean_1x": means[0],
        "mean_2x": means[1],
        "mean_paired_difference": float(np.mean(difference)),
        "ci95_half_width_paired_difference": (
            float(1.96 * np.std(difference, ddof=1) / math.sqrt(difference.size))
            if difference.size > 1 else 0.0
        ),
    }


def design_figure(rows: list[dict], figure_dir: Path) -> None:
    fixed = groups(rows, "fixed_spontaneous")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    for volume, radius in ((1, 1.0), (2, 2 ** (1 / 3))):
        circle = plt.Circle((0, 0), radius, fill=False, lw=3, color=COLOURS[volume])
        axes[0].add_patch(circle)
    axes[0].set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), aspect="equal")
    axes[0].axis("off")
    axes[0].set_title("Pure geometry change")
    axes[0].text(0, -1.43, "radius ×1.260  |  area ×1.587",
                 ha="center", fontsize=10)
    axes[1].axis("off")
    axes[1].set_title("Shared Fig. 4F model")
    axes[1].text(
        0.02, 0.94,
        "Unchanged\n"
        "• U and V concentrations\n"
        "• reaction and diffusion constants\n"
        "• physical mesh spacing\n"
        "• event intensity per area-time\n"
        "• event width and amplitude law\n"
        "• interface/mechanical coefficients\n\n"
        "Changed\n"
        "• target area and initial radius\n"
        "• computational-box extent\n"
        "• boundary mobility between arms",
        va="top", fontsize=10.2,
    )
    for volume in (1, 2):
        event_rate = np.asarray([row["event_rate_per_cell_time"] for row in fixed[volume]])
        area_rate = np.asarray([row["event_rate_per_area_time"] for row in fixed[volume]])
        axes[2].scatter(np.full(event_rate.size, volume), event_rate,
                        color=COLOURS[volume], alpha=0.5, s=18)
        mean = float(np.mean(event_rate))
        half = 1.96 * float(np.std(event_rate, ddof=1)) / math.sqrt(event_rate.size)
        axes[2].errorbar(volume, mean, yerr=half, color="black", capsize=4,
                         marker="o", zorder=3)
        axes[2].text(volume, mean, f"  area-normalized\n  {np.mean(area_rate):.4g}",
                     va="bottom", fontsize=8)
    axes[2].set_xticks([1, 2], ["1×", "2×"])
    axes[2].set_ylabel("Stochastic attempts / cell / time")
    axes[2].set_title("Prescribed per-area firing process")
    axes[2].grid(axis="y", color="#dddddd", lw=0.6)
    fig.suptitle("Taniguchi geometry campaign — design and realized firing", fontsize=14)
    fig.savefig(figure_dir / "01_campaign_design_and_firing.png", dpi=180)
    plt.close(fig)


def fixed_summary_figure(rows: list[dict], figure_dir: Path) -> dict:
    group = groups(rows, "fixed_spontaneous")
    specs = (
        ("successful_wave_initiation_rate_per_area_time", "Wave initiations / area / time"),
        ("mean_simultaneous_active_patches", "Mean simultaneous patches"),
        ("maximum_simultaneous_active_patches", "Maximum simultaneous patches"),
        ("fraction_time_active", "Fraction of time active"),
        ("mean_active_area_fraction", "Mean PIP3-high area fraction"),
        ("reexcited_area_fraction", "Area excited at least twice"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.5), constrained_layout=True)
    summaries = {}
    for ax, (metric, label) in zip(axes.ravel(), specs):
        first, second = paired_arrays(group, metric)
        summaries[metric] = paired_axis(ax, first, second, label)
    fig.suptitle("Fixed circular domains — spontaneous-wave activity (32 paired seeds)",
                 fontsize=14)
    fig.savefig(figure_dir / "02_fixed_spontaneous_activity.png", dpi=180)
    plt.close(fig)
    return summaries


def wave_dynamics_figure(rows: list[dict], figure_dir: Path) -> dict:
    group = groups(rows, "fixed_spontaneous")
    specs = (
        ("median_wave_lifetime", "Median patch lifetime"),
        ("median_wave_speed", "Median leading-edge speed"),
        ("mean_wave_travel_per_diameter", "Travel / cell diameter"),
        ("visited_area_fraction", "Fraction of cell ever active"),
    )
    fig, axes = plt.subplots(1, 4, figsize=(14.5, 4.1), constrained_layout=True)
    summaries = {}
    for ax, (metric, label) in zip(axes, specs):
        first, second = paired_arrays(group, metric)
        summaries[metric] = paired_axis(ax, first, second, label)
    fig.suptitle("Fixed circular domains — tracked wave dynamics", fontsize=14)
    fig.savefig(figure_dir / "03_fixed_wave_dynamics.png", dpi=180)
    plt.close(fig)
    return summaries


def regime_figure(rows: list[dict], figure_dir: Path) -> None:
    group = groups(rows, "fixed_spontaneous")
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    maximum = max(int(row["maximum_simultaneous_active_patches"])
                  for volume in (1, 2) for row in group[volume])
    bins = np.arange(-0.5, maximum + 1.5, 1)
    for volume in (1, 2):
        values = [row["maximum_simultaneous_active_patches"] for row in group[volume]]
        axes[0].hist(values, bins=bins, color=COLOURS[volume], alpha=0.55,
                     label=LABELS[volume])
    axes[0].set_xlabel("Maximum simultaneous active patches")
    axes[0].set_ylabel("Replicates")
    axes[0].legend(frameon=False)
    outcome_names = (
        ("extinction_candidate_count", "Extinction"),
        ("reflection_candidate_count", "Reflection"),
        ("collision_merge_count", "Collision/merge"),
    )
    x = np.arange(len(outcome_names))
    width = 0.34
    for offset, volume in ((-width / 2, 1), (width / 2, 2)):
        rates = []
        errors = []
        for metric, _ in outcome_names:
            values = np.asarray([
                row[metric] / max(row["wave_track_count"], 1)
                for row in group[volume]
            ])
            rates.append(float(np.mean(values)))
            errors.append(1.96 * float(np.std(values, ddof=1)) / math.sqrt(values.size))
        axes[1].bar(x + offset, rates, width, yerr=errors, color=COLOURS[volume],
                    alpha=0.8, capsize=3, label=LABELS[volume])
    axes[1].set_xticks(x, [label for _, label in outcome_names], rotation=20)
    axes[1].set_ylabel("Candidates per tracked patch")
    axes[1].legend(frameon=False)
    for volume in (1, 2):
        values = np.asarray([
            row["fraction_frames_with_phase_singularity_candidate"]
            for row in group[volume]
        ])
        axes[2].scatter(np.full(values.size, volume), values,
                        color=COLOURS[volume], alpha=0.65, s=20)
        axes[2].plot(volume, np.mean(values), "ko", ms=5)
    axes[2].set_xticks([1, 2], ["1×", "2×"])
    axes[2].set_ylabel("Fraction frames with phase winding")
    axes[2].set_title("Operational spiral-core candidates")
    fig.suptitle("Fixed-domain dynamical-regime diagnostics", fontsize=14)
    fig.savefig(figure_dir / "04_fixed_regime_outcomes.png", dpi=180)
    plt.close(fig)


def representative_seed(group: dict[int, list[dict]]) -> int:
    second = sorted(
        group[2],
        key=lambda row: (
            row["maximum_simultaneous_active_patches"],
            row["fraction_time_active"],
            row["phase_singularity_candidate_frames"],
        ),
        reverse=True,
    )
    return int(second[0]["seed"])


def snapshot_figure(group: dict[int, list[dict]], seed: int, figure_dir: Path,
                    filename: str, title: str) -> None:
    fig, axes = plt.subplots(2, 5, figsize=(14, 5.8), constrained_layout=True)
    for row_index, volume in enumerate((1, 2)):
        row = next(item for item in group[volume] if int(item["seed"]) == seed)
        with np.load(Path(row["run_dir"]) / "snapshots.npz") as data:
            times = data["times"]
            phi = data["phi"]
            V = data["V"]
            x, y = data["x"], data["y"]
        indices = np.linspace(0, len(times) - 1, 5).round().astype(int)
        vmax = max(float(np.percentile(V[phi >= 0.5], 99.5)), 0.5)
        for column, index in enumerate(indices):
            ax = axes[row_index, column]
            masked = np.ma.masked_where(phi[index] < 0.5, V[index])
            ax.imshow(masked, origin="lower", extent=(x[0], x[-1], y[0], y[-1]),
                      cmap="magma", vmin=0, vmax=vmax)
            ax.contour(x, y, phi[index], levels=[0.5], colors="white", linewidths=0.7)
            ax.set_title(f"t={times[index]:.1f}", fontsize=9)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            if column == 0:
                ax.set_ylabel(LABELS[volume], color=COLOURS[volume], fontweight="bold")
    fig.suptitle(f"{title} — activity-rich paired seed {seed}", fontsize=14)
    fig.savefig(figure_dir / filename, dpi=180)
    plt.close(fig)


def controlled_figure(rows: list[dict], figure_dir: Path) -> None:
    group = groups(rows, "fixed_controlled")
    fig, axes = plt.subplots(2, 5, figsize=(14, 6.2), constrained_layout=True)
    for row_index, volume in enumerate((1, 2)):
        row = group[volume][0]
        run = Path(row["run_dir"])
        with np.load(run / "snapshots.npz") as data:
            times_s, phi, V, x, y = (data[key] for key in ("times", "phi", "V", "x", "y"))
            indices = np.linspace(0, len(times_s) - 1, 4).round().astype(int)
            vmax = max(float(np.percentile(V[phi >= 0.5], 99.5)), 0.5)
            for column, index in enumerate(indices):
                ax = axes[row_index, column]
                ax.imshow(np.ma.masked_where(phi[index] < 0.5, V[index]), origin="lower",
                          extent=(x[0], x[-1], y[0], y[-1]), cmap="magma",
                          vmin=0, vmax=vmax)
                ax.contour(x, y, phi[index], levels=[0.5], colors="white", linewidths=0.7)
                ax.set_title(f"t={times_s[index]:.2f}", fontsize=9)
                ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        with np.load(run / "trajectory.npz") as data:
            ax = axes[row_index, 4]
            ax.plot(data["times"], data["active_area_fraction"], color=COLOURS[volume],
                    label="currently high")
            ax.plot(data["times"], data["visited_area_fraction"], color="black",
                    ls="--", label="ever reached")
            ax.set_ylim(0, 1)
            ax.set_xlabel("Time")
            ax.set_ylabel("Cell fraction")
            ax.grid(color="#dddddd", lw=0.5)
            if row_index == 0:
                ax.legend(frameon=False, fontsize=8)
        axes[row_index, 0].set_ylabel(LABELS[volume], color=COLOURS[volume],
                                      fontweight="bold")
    fig.suptitle("Controlled off-centre excitation — identical physical trigger", fontsize=14)
    fig.savefig(figure_dir / "06_controlled_single_wave.png", dpi=180)
    plt.close(fig)


def load_trajectory(row: dict) -> dict[str, np.ndarray]:
    with np.load(Path(row["run_dir"]) / "trajectory.npz") as data:
        return {key: data[key].copy() for key in data.files}


def migration_figure(rows: list[dict], figure_dir: Path) -> None:
    group = groups(rows, "deformable_spontaneous")
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.8), constrained_layout=True,
                             sharex=True, sharey=True)
    extent = 0.0
    paths: dict[int, list[np.ndarray]] = {1: [], 2: []}
    for volume in (1, 2):
        for row in group[volume]:
            data = load_trajectory(row)
            path = np.column_stack((data["centroid_x"] - data["centroid_x"][0],
                                    data["centroid_y"] - data["centroid_y"][0]))
            paths[volume].append(path)
            extent = max(extent, float(np.max(np.abs(path))))
    extent = max(0.05, 1.08 * extent)
    for ax, volume in zip(axes, (1, 2)):
        for path in paths[volume]:
            ax.plot(path[:, 0], path[:, 1], color=COLOURS[volume], alpha=0.28, lw=0.75)
            ax.scatter(path[-1, 0], path[-1, 1], color=COLOURS[volume], alpha=0.45, s=8)
        ax.scatter(0, 0, marker="+", color="black", s=60)
        ax.set(xlim=(-extent, extent), ylim=(-extent, extent), aspect="equal")
        ax.set_title(LABELS[volume], color=COLOURS[volume], fontweight="bold")
        ax.set_xlabel("x displacement (model length units)")
        ax.grid(color="#dddddd", lw=0.5)
    axes[0].set_ylabel("y displacement (model length units)")
    fig.suptitle("Deformable cells — centroid trajectories over 400 model-time units",
                 fontsize=14)
    fig.savefig(figure_dir / "07_deformable_centroid_trajectories.png", dpi=180)
    plt.close(fig)


def migration_dynamics_figure(rows: list[dict], figure_dir: Path) -> None:
    """Compare migration dynamics in native, unnormalized length units."""
    group = groups(rows, "deformable_spontaneous")
    trajectories = {
        volume: [load_trajectory(row) for row in group[volume]]
        for volume in (1, 2)
    }
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.4), constrained_layout=True)

    # Ensemble displacement squared from the common starting point.
    for volume in (1, 2):
        times = trajectories[volume][0]["times"]
        squared = []
        for data in trajectories[volume]:
            dx = data["centroid_x"] - data["centroid_x"][0]
            dy = data["centroid_y"] - data["centroid_y"][0]
            squared.append(dx * dx + dy * dy)
        values = np.vstack(squared)
        mean = np.mean(values, axis=0)
        half = 1.96 * np.std(values, axis=0, ddof=1) / math.sqrt(values.shape[0])
        axes[0].plot(times, mean, color=COLOURS[volume], label=LABELS[volume])
        axes[0].fill_between(times, np.maximum(0, mean - half), mean + half,
                             color=COLOURS[volume], alpha=0.16, linewidth=0)
    axes[0].set(xlabel="Time", ylabel="Mean squared displacement\n(model length units²)",
                title="Displacement from start")

    # Time-averaged MSD uses all valid time origins within each trajectory;
    # uncertainty is across the 32 replicate-level estimates.
    base_dt = float(np.median(np.diff(trajectories[1][0]["times"])))
    lag_steps = np.unique(np.geomspace(
        max(1, int(round(0.5 / base_dt))),
        min(int(round(100.0 / base_dt)), len(trajectories[1][0]["times"]) // 2),
        32,
    ).round().astype(int))
    for volume in (1, 2):
        per_run = []
        for data in trajectories[volume]:
            coordinates = np.column_stack((data["centroid_x"], data["centroid_y"]))
            per_run.append([
                float(np.mean(np.sum((coordinates[lag:] - coordinates[:-lag]) ** 2,
                                     axis=1)))
                for lag in lag_steps
            ])
        values = np.asarray(per_run)
        mean = np.mean(values, axis=0)
        half = 1.96 * np.std(values, axis=0, ddof=1) / math.sqrt(values.shape[0])
        lag_times = lag_steps * base_dt
        axes[1].plot(lag_times, mean, color=COLOURS[volume], label=LABELS[volume])
        axes[1].fill_between(lag_times, np.maximum(0, mean - half), mean + half,
                             color=COLOURS[volume], alpha=0.16, linewidth=0)
    axes[1].set(xlabel="Lag time", ylabel="Time-averaged MSD\n(model length units²)",
                title="All time origins")

    # Turning angles are formed from successive, non-overlapping 2-time-unit
    # displacement vectors.  Pooling is used only for the distributional plot.
    interval = 2.0
    mean_angles: dict[int, list[float]] = {1: [], 2: []}
    for volume in (1, 2):
        for data in trajectories[volume]:
            angles: list[float] = []
            stride = max(1, int(round(interval / np.median(np.diff(data["times"])))))
            coordinates = np.column_stack((data["centroid_x"], data["centroid_y"]))[::stride]
            vectors = np.diff(coordinates, axis=0)
            norms = np.linalg.norm(vectors, axis=1)
            good = norms > 1e-5
            for first, second, first_good, second_good in zip(
                vectors[:-1], vectors[1:], good[:-1], good[1:]
            ):
                if first_good and second_good:
                    cosine = np.dot(first, second) / (
                        np.linalg.norm(first) * np.linalg.norm(second)
                    )
                    angles.append(float(np.arccos(np.clip(cosine, -1.0, 1.0))))
            mean_angles[volume].append(float(np.degrees(np.mean(angles))))
    paired_axis(
        axes[2], np.asarray(mean_angles[1]), np.asarray(mean_angles[2]),
        "Mean turning angle (degrees)",
    )
    axes[2].set(xlabel="Volume", title="Successive 2-time-unit steps")
    for ax in axes:
        ax.grid(color="#dddddd", lw=0.5)
    axes[0].legend(frameon=False)
    fig.suptitle(
        "Deformable cells — migration dynamics in native model length units",
        fontsize=14,
    )
    fig.savefig(figure_dir / "12_deformable_msd_turning.png", dpi=180)
    plt.close(fig)


def deformable_metrics_figure(rows: list[dict], figure_dir: Path) -> dict:
    group = groups(rows, "deformable_spontaneous")
    specs = (
        ("net_displacement", "Net displacement"),
        ("path_length", "Path length"),
        ("mean_speed", "Mean speed"),
        ("straightness", "Net / path length"),
        ("mean_aspect_ratio", "Mean aspect ratio"),
        ("sd_circularity", "Circularity fluctuations (SD)"),
        ("mean_protrusion_count", "Mean protrusion count"),
        ("pip3_to_future_boundary_extension_correlation", "PIP3 → boundary extension corr."),
    )
    fig, axes = plt.subplots(2, 4, figsize=(14.8, 7.3), constrained_layout=True)
    summaries = {}
    for ax, (metric, label) in zip(axes.ravel(), specs):
        first, second = paired_arrays(group, metric)
        if first.size:
            summaries[metric] = paired_axis(ax, first, second, label)
        else:
            ax.text(0.5, 0.5, "Not estimable", ha="center", va="center")
            ax.set_axis_off()
    fig.suptitle("Deformable cells — migration, morphology, and wave–boundary coupling",
                 fontsize=14)
    fig.savefig(figure_dir / "08_deformable_cell_metrics.png", dpi=180)
    plt.close(fig)
    return summaries


def deformable_timecourses(rows: list[dict], figure_dir: Path) -> None:
    group = groups(rows, "deformable_spontaneous")
    specs = (
        ("active_area_fraction", "PIP3-high area fraction"),
        ("active_patch_count", "Active patches"),
        ("aspect_ratio", "Aspect ratio"),
        ("circularity", "Circularity"),
        ("protrusion_count", "Protrusion count"),
        ("boundary_radial_roughness", "Boundary roughness"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), constrained_layout=True)
    for ax, (metric, label) in zip(axes.ravel(), specs):
        for volume in (1, 2):
            trajectories = [load_trajectory(row) for row in group[volume]]
            times = trajectories[0]["times"]
            values = np.vstack([trajectory[metric] for trajectory in trajectories])
            mean = np.mean(values, axis=0)
            half = 1.96 * np.std(values, axis=0, ddof=1) / math.sqrt(values.shape[0])
            # A 2-time-unit running mean keeps the ensemble figure legible.
            width = max(1, int(round(2.0 / np.median(np.diff(times)))))
            kernel = np.ones(width) / width
            left = width // 2
            right = width - 1 - left
            smooth_mean = np.convolve(
                np.pad(mean, (left, right), mode="edge"), kernel, mode="valid"
            )
            smooth_half = np.convolve(
                np.pad(half, (left, right), mode="edge"), kernel, mode="valid"
            )
            ax.plot(times, smooth_mean, color=COLOURS[volume], label=LABELS[volume])
            ax.fill_between(times, smooth_mean - smooth_half, smooth_mean + smooth_half,
                            color=COLOURS[volume], alpha=0.15, linewidth=0)
        ax.set_xlabel("Time")
        ax.set_ylabel(label)
        ax.grid(color="#dddddd", lw=0.5)
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Deformable cells — ensemble time courses (mean ± 95% CI, smoothed)",
                 fontsize=14)
    fig.savefig(figure_dir / "09_deformable_timecourses.png", dpi=180)
    plt.close(fig)


def coupling_heatmap(rows: list[dict], seed: int, figure_dir: Path) -> None:
    group = groups(rows, "deformable_spontaneous")
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), constrained_layout=True)
    for row_index, volume in enumerate((1, 2)):
        row = next(item for item in group[volume] if int(item["seed"]) == seed)
        data = load_trajectory(row)
        extent = (0, 360, float(data["times"][-1]), float(data["times"][0]))
        first = axes[row_index, 0].imshow(data["boundary_pip3_profiles"], aspect="auto",
                                          extent=extent, cmap="magma")
        second = axes[row_index, 1].imshow(
            data["radial_profiles"] / row["initial_radius"], aspect="auto",
            extent=extent, cmap="viridis"
        )
        axes[row_index, 0].set_ylabel(f"{LABELS[volume]}\nTime",
                                      color=COLOURS[volume], fontweight="bold")
        fig.colorbar(first, ax=axes[row_index, 0], label="Boundary PIP3")
        fig.colorbar(second, ax=axes[row_index, 1], label="Radius / initial radius")
    axes[0, 0].set_title("Chemical activity around boundary")
    axes[0, 1].set_title("Subsequent boundary geometry")
    axes[1, 0].set_xlabel("Boundary angle (degrees)")
    axes[1, 1].set_xlabel("Boundary angle (degrees)")
    fig.suptitle(f"Wave–boundary coupling map — activity-rich paired seed {seed}", fontsize=14)
    fig.savefig(figure_dir / "10_wave_boundary_coupling_maps.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root", type=Path,
        default=ROOT / "results/wgd_v1/taniguchi_geometry_v1",
    )
    args = parser.parse_args()
    output = args.root if args.root.is_absolute() else ROOT / args.root
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    rows = load_runs(output)
    write_csv(figure_dir / "run_summaries.csv", rows)
    write_track_csv(figure_dir / "wave_tracks.csv", rows)
    design_figure(rows, figure_dir)
    fixed_summaries = fixed_summary_figure(rows, figure_dir)
    wave_summaries = wave_dynamics_figure(rows, figure_dir)
    regime_figure(rows, figure_dir)
    fixed_group = groups(rows, "fixed_spontaneous")
    fixed_seed = representative_seed(fixed_group)
    snapshot_figure(fixed_group, fixed_seed, figure_dir,
                    "05_fixed_representative_waves.png", "Fixed-domain PIP3 waves")
    controlled_figure(rows, figure_dir)
    migration_figure(rows, figure_dir)
    migration_dynamics_figure(rows, figure_dir)
    deformable_summaries = deformable_metrics_figure(rows, figure_dir)
    deformable_timecourses(rows, figure_dir)
    deformable_group = groups(rows, "deformable_spontaneous")
    deformable_seed = representative_seed(deformable_group)
    snapshot_figure(deformable_group, deformable_seed, figure_dir,
                    "11_deformable_representative_waves.png",
                    "Deformable-cell PIP3 waves and shapes")
    coupling_heatmap(rows, deformable_seed, figure_dir)
    controlled = {
        str(volume): groups(rows, "fixed_controlled")[volume][0]
        for volume in (1, 2)
    }
    manifest = {
        "campaign": "taniguchi_geometry_v1",
        "source": str(output),
        "report_integration": "none",
        "figure_count": len(list(figure_dir.glob("*.png"))),
        "figures": sorted(path.name for path in figure_dir.glob("*.png")),
        "derived_tables": ["run_summaries.csv", "wave_tracks.csv"],
        "fixed_representative_seed": fixed_seed,
        "deformable_representative_seed": deformable_seed,
        "representative_selection": (
            "Paired seed ranked by 2x maximum simultaneous patches, active-time fraction, "
            "then phase-singularity-candidate frames; illustrative rather than typical."
        ),
        "fixed_paired_summaries": fixed_summaries,
        "fixed_wave_paired_summaries": wave_summaries,
        "deformable_paired_summaries": deformable_summaries,
        "controlled_summaries": controlled,
        "operational_cautions": {
            "event": "A stochastic U-to-V firing attempt; distinct from a successful tracked wave birth.",
            "patch": "An eight-connected V>=0.5 region of area at least 0.2 inside phi>=0.5.",
            "reflection": "A tracked patch contacts the 0.5-unit boundary zone and later detaches while persisting.",
            "spiral": "A state-space phase-winding candidate, requiring persistence for biological interpretation.",
            "collision": "A merge of overlapping tracked active components; annihilating near-misses are not counted.",
        },
    }
    (figure_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(figure_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
