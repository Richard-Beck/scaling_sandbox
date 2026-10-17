#!/usr/bin/env python3
"""Create standalone diagnostic figures for the Taniguchi WGD production run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
COLOURS = {1: "#0072B2", 2: "#D55E00"}
LABELS = {1: "1× volume", 2: "2× volume"}


def style() -> None:
    plt.rcParams.update({
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.2,
        "figure.dpi": 150,
        "savefig.dpi": 240,
    })


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float:
    return float(value) if value not in {"", "NA", "nan"} else np.nan


def mean_ci(values: np.ndarray, axis: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Return the mean and normal-approximation 95% half-width."""
    values = np.asarray(values, dtype=float)
    count = np.sum(np.isfinite(values), axis=axis)
    mean = np.nanmean(values, axis=axis)
    sd = np.nanstd(values, axis=axis, ddof=1)
    half_width = 1.96 * sd / np.sqrt(np.maximum(count, 1))
    return mean, half_width


def write_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_rows(output: Path) -> list[dict[str, str]]:
    return read_csv(output / "comparison" / "run_summary.csv")


def paired_values(rows: list[dict[str, str]], statistic: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    by_seed: dict[int, dict[int, float]] = {}
    for row in rows:
        seed = int(row["seed"])
        volume = int(float(row["volume_ratio"]))
        by_seed.setdefault(seed, {})[volume] = as_float(row[statistic])
    seeds = np.asarray(sorted(seed for seed, values in by_seed.items()
                              if 1 in values and 2 in values), dtype=int)
    one = np.asarray([by_seed[int(seed)][1] for seed in seeds], dtype=float)
    two = np.asarray([by_seed[int(seed)][2] for seed in seeds], dtype=float)
    return seeds, one, two


def choose_representative_seed(rows: list[dict[str, str]]) -> tuple[int, str]:
    """Choose an active pair closest to the multimetric production center."""
    statistics = (
        "fraction_frames_with_active_domain",
        "maximum_active_area_fraction",
        "directional_persistence",
        "stochastic_event_rate_per_time",
    )
    matrices = []
    seeds = None
    for statistic in statistics:
        current_seeds, one, two = paired_values(rows, statistic)
        if seeds is None:
            seeds = current_seeds
        elif not np.array_equal(seeds, current_seeds):
            raise ValueError("paired seed sets differ across statistics")
        matrices.extend((one, two))
    assert seeds is not None
    matrix = np.column_stack(matrices)
    center = np.nanmedian(matrix, axis=0)
    scale = np.nanmedian(np.abs(matrix - center), axis=0)
    fallback = np.nanstd(matrix, axis=0, ddof=1)
    scale = np.where(scale > 0, scale, np.where(fallback > 0, fallback, 1.0))
    distance = np.nanmean(((matrix - center) / scale) ** 2, axis=1)
    # Field panels are useful only when both paired runs exhibit some activity.
    active_both = (matrix[:, 0] > 0) & (matrix[:, 1] > 0)
    distance[~active_both] = np.inf
    index = int(np.argmin(distance))
    if not np.isfinite(distance[index]):
        raise ValueError("no paired seed has activity in both conditions")
    rule = (
        "Active in both conditions; smallest robust standardized distance to "
        "the production medians across four trajectory-wide metrics."
    )
    return int(seeds[index]), rule


def allocated_bytes(paths: list[Path]) -> int:
    return sum(path.stat().st_blocks * 512 for path in paths if path.is_file())


def output_format_figure(output: Path, figure_dir: Path) -> dict[str, float]:
    trajectories = sorted(output.glob("*x_seed*/trajectory.npz"))
    checkpoints = sorted(output.glob("*x_seed*/checkpoints.npz"))
    all_files = [
        path for path in output.rglob("*")
        if path.is_file() and figure_dir not in path.parents
    ]
    raw_set = set(trajectories + checkpoints)
    other = [path for path in all_files if path not in raw_set]
    storage = {
        "compact_observables_mb": allocated_bytes(trajectories) / 1e6,
        "full_field_checkpoints_mb": allocated_bytes(checkpoints) / 1e6,
        "configs_summaries_tables_mb": allocated_bytes(other) / 1e6,
    }
    storage["total_mb"] = sum(storage.values())
    storage["legacy_projection_mb"] = 4700.0

    with np.load(trajectories[0], allow_pickle=False) as trajectory:
        fine_times = trajectory["times"]
    with np.load(checkpoints[0], allow_pickle=False) as checkpoint:
        coarse_times = checkpoint["times"]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)
    ax = axes[0]
    ax.hlines(1, fine_times[0], fine_times[-1], color=COLOURS[1], lw=2)
    ax.vlines(fine_times, 0.91, 1.09, color=COLOURS[1], lw=0.35, alpha=0.65)
    ax.hlines(0, coarse_times[0], coarse_times[-1], color=COLOURS[2], lw=2)
    ax.vlines(coarse_times, -0.13, 0.13, color=COLOURS[2], lw=2)
    ax.set_yticks((0, 1), ("Full fields\n$\\phi$, PIP2, PIP3", "Compact\nobservables"))
    ax.set_xlabel("Model time")
    ax.set_title("Two-resolution retained output")
    ax.text(0.5, 1.20, "401 frames  (Δt = 0.05)", color=COLOURS[1],
            transform=ax.get_yaxis_transform(), va="center")
    ax.text(0.5, -0.25, "11 checkpoints  (Δt = 2)", color=COLOURS[2],
            transform=ax.get_yaxis_transform(), va="center")
    ax.set_ylim(-0.55, 1.5)

    ax = axes[1]
    labels = ("Compact\nproduction", "Legacy full-field\nprojection")
    actual_parts = (
        storage["compact_observables_mb"],
        storage["full_field_checkpoints_mb"],
        storage["configs_summaries_tables_mb"],
    )
    bottom = 0.0
    part_labels = ("Compact observables", "Full-field checkpoints", "Other")
    part_colours = (COLOURS[1], COLOURS[2], "#999999")
    for value, label, colour in zip(actual_parts, part_labels, part_colours):
        ax.bar(0, value, bottom=bottom, width=0.58, color=colour, label=label)
        bottom += value
    ax.bar(1, storage["legacy_projection_mb"], width=0.58, color="#BBBBBB",
           hatch="///", edgecolor="#777777", label="Original projection")
    ax.set_yscale("log")
    ax.set_xticks((0, 1), labels)
    ax.set_ylabel("Allocated storage (MB, log scale)")
    ax.set_title("Production footprint: 64 runs")
    ax.text(0, storage["total_mb"] * 1.08, f"{storage['total_mb']:.0f} MB",
            ha="center", va="bottom", fontweight="bold")
    ax.text(1, storage["legacy_projection_mb"] * 1.08, "~4.7 GB",
            ha="center", va="bottom", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Taniguchi WGD output design", fontsize=14)
    fig.savefig(figure_dir / "output_format_and_footprint.png", bbox_inches="tight")
    plt.close(fig)
    return storage


def paired_fields_figure(output: Path, figure_dir: Path, seed: int) -> list[float]:
    loaded = {}
    pooled_activity = None
    checkpoint_times = None
    for volume in (1, 2):
        run = output / f"{volume}x_seed{seed}"
        with np.load(run / "trajectory.npz", allow_pickle=False) as trajectory:
            activity = trajectory["active_fractions"].copy()
        with np.load(run / "checkpoints.npz", allow_pickle=False) as checkpoint:
            data = {key: checkpoint[key].copy() for key in checkpoint.files}
        checkpoint_activity = activity[data["source_indices"]]
        loaded[volume] = data
        if checkpoint_times is None:
            checkpoint_times = data["times"]
            pooled_activity = checkpoint_activity.copy()
        else:
            if not np.allclose(checkpoint_times, data["times"]):
                raise ValueError("paired checkpoint times differ")
            pooled_activity += checkpoint_activity
    assert checkpoint_times is not None and pooled_activity is not None
    candidates = np.arange(1, checkpoint_times.size)
    selected = [0]
    for index in candidates[np.argsort(pooled_activity[candidates])[::-1]]:
        if all(abs(int(index) - old) >= 1 for old in selected):
            selected.append(int(index))
        if len(selected) == 4:
            break
    selected = sorted(selected)

    vmax = max(float(loaded[volume]["V"][selected].max()) for volume in (1, 2))
    norm = Normalize(vmin=0, vmax=max(vmax, np.finfo(float).eps))
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.1), constrained_layout=True)
    artist = None
    for row, volume in enumerate((1, 2)):
        data = loaded[volume]
        extent = (data["x"][0], data["x"][-1], data["y"][0], data["y"][-1])
        for column, index in enumerate(selected):
            ax = axes[row, column]
            support = data["phi"][index] >= 0.5
            field = np.ma.masked_where(~support, data["V"][index])
            artist = ax.imshow(field, origin="lower", extent=extent, cmap="magma",
                               norm=norm, interpolation="bilinear")
            ax.contour(data["x"], data["y"], data["phi"][index], levels=(0.5,),
                       colors="white", linewidths=0.8)
            ax.set_aspect("equal")
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"t = {data['times'][index]:g}")
            if column == 0:
                ax.set_ylabel(LABELS[volume], color=COLOURS[volume],
                              fontweight="bold")
    assert artist is not None
    fig.colorbar(artist, ax=axes, shrink=0.77, label="PIP3 field V")
    fig.suptitle(f"Paired full-field checkpoints — representative seed {seed}", fontsize=14)
    fig.savefig(figure_dir / "representative_paired_fields.png", bbox_inches="tight")
    plt.close(fig)
    return [float(checkpoint_times[index]) for index in selected]


def time_course_figure(output: Path, figure_dir: Path, seed: int) -> None:
    trajectories = {}
    for volume in (1, 2):
        path = output / f"{volume}x_seed{seed}" / "trajectory.npz"
        with np.load(path, allow_pickle=False) as data:
            trajectories[volume] = {key: data[key].copy() for key in data.files}
    specifications = (
        ("active_fractions", "Active area / cell area"),
        ("domain_counts", "Number of active domains"),
        ("directional_moments", "Directional moment"),
        ("shape_anisotropy", "Cell-shape anisotropy"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), sharex=True,
                             constrained_layout=True)
    for ax, (key, ylabel) in zip(axes.flat, specifications):
        for volume in (1, 2):
            data = trajectories[volume]
            ax.plot(data["times"], data[key], color=COLOURS[volume], lw=1.5,
                    label=LABELS[volume])
        ax.set_ylabel(ylabel)
    for ax in axes[-1]:
        ax.set_xlabel("Model time")
    axes[0, 0].legend(frameon=False)
    fig.suptitle(f"Compact trajectory observables — representative seed {seed}", fontsize=14)
    fig.savefig(figure_dir / "representative_compact_trajectories.png",
                bbox_inches="tight")
    plt.close(fig)


def load_migration_trajectories(
    output: Path,
) -> dict[int, dict[int, dict[str, np.ndarray | float]]]:
    trajectories: dict[int, dict[int, dict[str, np.ndarray | float]]] = {1: {}, 2: {}}
    for path in sorted(output.glob("*x_seed*/trajectory.npz")):
        condition, seed_text = path.parent.name.split("_seed")
        volume = int(condition.removesuffix("x"))
        seed = int(seed_text)
        with np.load(path, allow_pickle=False) as data:
            trajectories[volume][seed] = {
                "times": data["times"].copy(),
                "x": data["centroid_x"].copy(),
                "y": data["centroid_y"].copy(),
                "area0": float(data["area0"]),
            }
    seeds_1x = set(trajectories[1])
    seeds_2x = set(trajectories[2])
    if seeds_1x != seeds_2x or not seeds_1x:
        raise ValueError("migration trajectories do not form a nonempty paired design")
    return trajectories


def migration_metrics(data: dict[str, np.ndarray | float]) -> dict[str, float]:
    times = np.asarray(data["times"], dtype=float)
    x = np.asarray(data["x"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    radius = float(np.sqrt(float(data["area0"]) / np.pi))
    coordinates = np.column_stack((x - x[0], y - y[0]))
    steps = np.diff(coordinates, axis=0)
    step_lengths = np.linalg.norm(steps, axis=1)
    displacement = float(np.linalg.norm(coordinates[-1]))
    path_length = float(np.sum(step_lengths))
    duration = float(times[-1] - times[0])
    return {
        "net_displacement": displacement,
        "path_length": path_length,
        "maximum_excursion": float(np.max(np.linalg.norm(coordinates, axis=1))),
        "straightness": displacement / path_length if path_length > 0 else 0.0,
        "net_displacement_per_radius": displacement / radius,
        "path_length_per_radius": path_length / radius,
        "maximum_excursion_per_radius": (
            float(np.max(np.linalg.norm(coordinates, axis=1))) / radius
        ),
        "mean_speed_per_radius": path_length / radius / duration if duration > 0 else 0.0,
    }


def centroid_trajectory_figure(
    trajectories: dict[int, dict[int, dict[str, np.ndarray | float]]],
    figure_dir: Path,
    representative_seed: int,
) -> None:
    normalized: dict[int, list[np.ndarray]] = {1: [], 2: []}
    for volume in (1, 2):
        for data in trajectories[volume].values():
            radius = float(np.sqrt(float(data["area0"]) / np.pi))
            x = np.asarray(data["x"], dtype=float)
            y = np.asarray(data["y"], dtype=float)
            normalized[volume].append(
                np.column_stack(((x - x[0]) / radius, (y - y[0]) / radius))
            )
    extent = max(float(np.max(np.abs(path))) for paths in normalized.values()
                 for path in paths)
    extent = max(0.05, 1.08 * extent)

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.8), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, volume in zip(axes, (1, 2)):
        for seed, data in trajectories[volume].items():
            radius = float(np.sqrt(float(data["area0"]) / np.pi))
            x = (np.asarray(data["x"], dtype=float) - float(np.asarray(data["x"])[0])) / radius
            y = (np.asarray(data["y"], dtype=float) - float(np.asarray(data["y"])[0])) / radius
            representative = seed == representative_seed
            ax.plot(x, y, color=COLOURS[volume],
                    lw=2.2 if representative else 0.85,
                    alpha=1.0 if representative else 0.35,
                    zorder=2 if representative else 1)
            ax.scatter(x[-1], y[-1], s=24 if representative else 8,
                       color=COLOURS[volume], alpha=1.0 if representative else 0.45,
                       edgecolor="black" if representative else "none",
                       linewidth=0.7, zorder=3)
        ax.scatter(0, 0, marker="+", s=55, color="black", lw=1.2, zorder=4)
        ax.axhline(0, color="#777777", lw=0.6, alpha=0.45)
        ax.axvline(0, color="#777777", lw=0.6, alpha=0.45)
        ax.set_aspect("equal")
        ax.set_xlim(-extent, extent)
        ax.set_ylim(-extent, extent)
        ax.set_title(LABELS[volume], color=COLOURS[volume], fontweight="bold")
        ax.set_xlabel("Centroid displacement x / initial radius")
    axes[0].set_ylabel("Centroid displacement y / initial radius")
    axes[1].plot([], [], color=COLOURS[2], lw=2.2,
                 label=f"Representative seed {representative_seed}")
    axes[1].scatter([], [], s=18, color="#777777", label="Final position")
    axes[1].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Centroid migration trajectories — 32 stochastic runs per condition",
                 fontsize=14)
    fig.savefig(figure_dir / "centroid_migration_trajectories.png", bbox_inches="tight")
    plt.close(fig)


def paired_migration_figure(
    trajectories: dict[int, dict[int, dict[str, np.ndarray | float]]],
    figure_dir: Path,
) -> tuple[list[dict[str, float | int | str]], dict[str, dict[str, float]]]:
    seeds = sorted(trajectories[1])
    per_run: dict[int, dict[int, dict[str, float]]] = {1: {}, 2: {}}
    table: list[dict[str, float | int | str]] = []
    for volume in (1, 2):
        for seed in seeds:
            metrics = migration_metrics(trajectories[volume][seed])
            per_run[volume][seed] = metrics
            table.append({"seed": seed, "volume_ratio": volume, **metrics})

    specifications = (
        ("net_displacement_per_radius", "Net displacement\n/ initial radius"),
        ("path_length_per_radius", "Path length\n/ initial radius"),
        ("maximum_excursion_per_radius", "Maximum excursion\n/ initial radius"),
        ("straightness", "Path straightness\n(net / path length)"),
    )
    summaries: dict[str, dict[str, float]] = {}
    rng = np.random.default_rng(20260823)
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.0), constrained_layout=True)
    for ax, (key, ylabel) in zip(axes.flat, specifications):
        one = np.asarray([per_run[1][seed][key] for seed in seeds])
        two = np.asarray([per_run[2][seed][key] for seed in seeds])
        jitter = rng.normal(0, 0.018, size=len(seeds))
        for a, b, j in zip(one, two, jitter):
            ax.plot((1 + j, 2 + j), (a, b), color="#A8A8A8", lw=0.7,
                    alpha=0.62, zorder=1)
        ax.scatter(1 + jitter, one, s=18, color=COLOURS[1], alpha=0.82, zorder=2)
        ax.scatter(2 + jitter, two, s=18, color=COLOURS[2], alpha=0.82, zorder=2)
        means = (float(np.mean(one)), float(np.mean(two)))
        ax.plot((1, 2), means, "o-", color="black", lw=2.2, markersize=6, zorder=3)
        ax.set_xticks((1, 2), ("1×", "2×"))
        ax.set_xlim(0.72, 2.28)
        ax.set_ylabel(ylabel)
        difference = two - one
        summaries[key] = {
            "n_pairs": len(seeds),
            "mean_1x": means[0],
            "mean_2x": means[1],
            "mean_paired_difference_2x_minus_1x": float(np.mean(difference)),
            "sd_paired_difference": float(np.std(difference, ddof=1)),
        }
        ax.text(0.5, 0.98, f"mean Δ = {np.mean(difference):.3g}",
                transform=ax.transAxes, ha="center", va="top", fontsize=9)
    fig.suptitle("Paired centroid-migration metrics — 32 seed pairs", fontsize=14)
    fig.savefig(figure_dir / "paired_migration_metrics.png", bbox_inches="tight")
    plt.close(fig)
    return table, summaries


def checkpoint_shape_metrics(
    phi: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    area0: float,
) -> dict[str, np.ndarray]:
    """Calculate scale-aware morphology from retained diffuse-cell fields."""
    dx = float(x[1] - x[0])
    pixel_area = dx * dx
    X, Y = np.meshgrid(x, y)
    target_perimeter = 2.0 * np.sqrt(np.pi * area0)
    result = {
        "relative_area": np.empty(phi.shape[0]),
        "relative_perimeter": np.empty(phi.shape[0]),
        "shape_anisotropy": np.empty(phi.shape[0]),
        "boundary_radial_roughness": np.empty(phi.shape[0]),
    }
    for index, frame in enumerate(np.asarray(phi, dtype=float)):
        weight = np.maximum(frame, 0.0)
        total = float(np.sum(weight))
        cx = float(np.sum(weight * X) / max(total, 1e-12))
        cy = float(np.sum(weight * Y) / max(total, 1e-12))
        xx = float(np.sum(weight * (X - cx) ** 2) / max(total, 1e-12))
        yy = float(np.sum(weight * (Y - cy) ** 2) / max(total, 1e-12))
        xy = float(np.sum(weight * (X - cx) * (Y - cy)) / max(total, 1e-12))
        eigenvalues = np.linalg.eigvalsh(np.asarray([[xx, xy], [xy, yy]]))
        anisotropy = ((eigenvalues[-1] - eigenvalues[0]) /
                      max(eigenvalues[-1] + eigenvalues[0], 1e-12))
        gradient_y, gradient_x = np.gradient(frame, dx, dx)
        perimeter = float(np.sum(np.hypot(gradient_x, gradient_y)) * pixel_area)

        support = frame >= 0.5
        boundary = support & ~ndimage.binary_erosion(
            support, structure=np.ones((3, 3), dtype=bool), border_value=0
        )
        radii = np.hypot(X[boundary] - cx, Y[boundary] - cy)
        radial_roughness = (
            float(np.std(radii) / np.mean(radii)) if radii.size else np.nan
        )
        result["relative_area"][index] = total * pixel_area / area0
        result["relative_perimeter"][index] = perimeter / target_perimeter
        result["shape_anisotropy"][index] = float(anisotropy)
        result["boundary_radial_roughness"][index] = radial_roughness
    return result


def checkpoint_morphology_figure(
    output: Path,
    figure_dir: Path,
) -> tuple[list[dict[str, float | int | str]], dict[str, dict[str, list[float]]]]:
    seeds = sorted(int(path.parent.name.split("_seed")[1])
                   for path in output.glob("1x_seed*/trajectory.npz"))
    profiles: dict[int, dict[str, list[np.ndarray]]] = {
        1: {}, 2: {},
    }
    table: list[dict[str, float | int | str]] = []
    times = None
    metric_names = (
        "relative_area", "relative_perimeter", "shape_anisotropy",
        "boundary_radial_roughness",
    )
    for volume in (1, 2):
        profiles[volume] = {name: [] for name in metric_names}
        for seed in seeds:
            run = output / f"{volume}x_seed{seed}"
            with np.load(run / "trajectory.npz", allow_pickle=False) as trajectory:
                area0 = float(trajectory["area0"])
            with np.load(run / "checkpoints.npz", allow_pickle=False) as checkpoint:
                current_times = checkpoint["times"].copy()
                metrics = checkpoint_shape_metrics(
                    checkpoint["phi"], checkpoint["x"], checkpoint["y"], area0
                )
            if times is None:
                times = current_times
            elif not np.allclose(times, current_times):
                raise ValueError("checkpoint time grids differ among production runs")
            for name in metric_names:
                profiles[volume][name].append(metrics[name])
            for index, time in enumerate(current_times):
                table.append({
                    "seed": seed,
                    "volume_ratio": volume,
                    "time": float(time),
                    **{name: float(metrics[name][index]) for name in metric_names},
                })
    assert times is not None
    arrays = {
        volume: {name: np.vstack(profiles[volume][name]) for name in metric_names}
        for volume in (1, 2)
    }
    specifications = (
        ("relative_area", "Area / target area"),
        ("relative_perimeter", "Diffuse perimeter /\ntarget-circle perimeter"),
        ("shape_anisotropy", "Second-moment\nshape anisotropy"),
        ("boundary_radial_roughness", "Boundary radial roughness\n(SD radius / mean radius)"),
    )
    summaries: dict[str, dict[str, list[float]]] = {}
    fig, axes = plt.subplots(2, 4, figsize=(14.2, 6.8), sharex="col",
                             constrained_layout=True)
    for column, (name, ylabel) in enumerate(specifications):
        summaries[name] = {}
        for volume in (1, 2):
            mean, half_width = mean_ci(arrays[volume][name], axis=0)
            axes[0, column].plot(times, mean, color=COLOURS[volume], lw=2,
                                 marker="o", ms=3.5, label=LABELS[volume])
            axes[0, column].fill_between(
                times, mean - half_width, mean + half_width,
                color=COLOURS[volume], alpha=0.16, linewidth=0,
            )
            summaries[name][f"mean_{volume}x"] = mean.tolist()
            summaries[name][f"ci95_half_width_{volume}x"] = half_width.tolist()
        paired_difference = arrays[2][name] - arrays[1][name]
        difference_mean, difference_half_width = mean_ci(paired_difference, axis=0)
        axes[1, column].axhline(0, color="#777777", lw=0.9, ls="--")
        axes[1, column].errorbar(
            times, difference_mean, yerr=difference_half_width,
            color="black", marker="o", ms=3.5, lw=1.5, capsize=2,
        )
        axes[0, column].set_ylabel(ylabel)
        axes[1, column].set_ylabel("Paired difference\n2× − 1×")
        axes[1, column].set_xlabel("Model time")
        summaries[name]["mean_paired_difference_2x_minus_1x"] = difference_mean.tolist()
        summaries[name]["ci95_half_width_paired_difference"] = difference_half_width.tolist()
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(
        "Checkpoint-resolved cell morphology — means and 95% CIs over 32 paired seeds\n"
        "Scale-normalized features; t = 0 offsets reflect the fixed diffuse-interface width",
        fontsize=13,
    )
    fig.savefig(figure_dir / "checkpoint_morphology_comparison.png",
                bbox_inches="tight")
    plt.close(fig)
    summaries["checkpoint_times"] = {"values": times.tolist()}
    return table, summaries


def paired_results_figure(rows: list[dict[str, str]], figure_dir: Path) -> dict[str, dict[str, float]]:
    specifications = (
        ("fraction_frames_with_active_domain", "Fraction of frames\nwith active domain"),
        ("maximum_active_area_fraction", "Maximum active\narea fraction"),
        ("maximum_n_domains", "Maximum simultaneous\nactive domains"),
        ("directional_persistence", "Directional persistence"),
        ("stochastic_event_rate_per_time", "Event rate\nper model time"),
        ("stochastic_event_rate_per_area_time", "Event rate\nper area×time"),
    )
    summaries: dict[str, dict[str, float]] = {}
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.1), constrained_layout=True)
    rng = np.random.default_rng(4473331)
    for ax, (statistic, ylabel) in zip(axes.flat, specifications):
        seeds, one, two = paired_values(rows, statistic)
        finite = np.isfinite(one) & np.isfinite(two)
        one, two = one[finite], two[finite]
        jitter = rng.normal(0, 0.018, size=one.size)
        for a, b, j in zip(one, two, jitter):
            ax.plot((1 + j, 2 + j), (a, b), color="#A8A8A8", lw=0.7, alpha=0.62,
                    zorder=1)
        ax.scatter(1 + jitter, one, s=18, color=COLOURS[1], alpha=0.82, zorder=2)
        ax.scatter(2 + jitter, two, s=18, color=COLOURS[2], alpha=0.82, zorder=2)
        means = (float(np.mean(one)), float(np.mean(two)))
        ax.plot((1, 2), means, "o-", color="black", lw=2.2, markersize=6, zorder=3)
        ax.set_xticks((1, 2), ("1×", "2×"))
        ax.set_xlim(0.72, 2.28)
        ax.set_ylabel(ylabel)
        difference = two - one
        summaries[statistic] = {
            "n_pairs": int(one.size),
            "mean_1x": means[0],
            "mean_2x": means[1],
            "mean_paired_difference_2x_minus_1x": float(np.mean(difference)),
            "sd_paired_difference": float(np.std(difference, ddof=1)),
        }
        ax.text(0.5, 0.98, f"mean Δ = {np.mean(difference):.3g}",
                transform=ax.transAxes, ha="center", va="top", fontsize=9)
    fig.suptitle("Taniguchi WGD production results — 32 paired stochastic seeds",
                 fontsize=14)
    fig.savefig(figure_dir / "paired_production_results.png", bbox_inches="tight")
    plt.close(fig)
    return summaries


def state_outcomes_figure(rows: list[dict[str, str]], figure_dir: Path) -> dict[str, dict[str, int]]:
    states = sorted({row["state_class"] for row in rows})
    counts = {
        volume: {state: sum(int(float(row["volume_ratio"])) == volume
                                  and row["state_class"] == state for row in rows)
                 for state in states}
        for volume in (1, 2)
    }
    x = np.arange(len(states))
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    width = 0.36
    ax.bar(x - width / 2, [counts[1][state] for state in states], width,
           color=COLOURS[1], label=LABELS[1])
    ax.bar(x + width / 2, [counts[2][state] for state in states], width,
           color=COLOURS[2], label=LABELS[2])
    ax.set_xticks(x, [state.replace("_", " ").title() for state in states])
    ax.set_ylabel("Number of stochastic runs")
    ax.set_ylim(0, 32)
    ax.legend(frameon=False)
    ax.set_title("Terminal state classification (32 runs per condition)")
    fig.savefig(figure_dir / "terminal_state_counts.png", bbox_inches="tight")
    plt.close(fig)
    return {LABELS[volume]: counts[volume] for volume in (1, 2)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results/wgd_v1/taniguchi_production",
    )
    parser.add_argument("--figure-dir", type=Path, default=None)
    args = parser.parse_args()
    figure_dir = args.figure_dir or args.output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    style()

    rows = run_rows(args.output)
    seed, selection_rule = choose_representative_seed(rows)
    storage = output_format_figure(args.output, figure_dir)
    checkpoint_times = paired_fields_figure(args.output, figure_dir, seed)
    time_course_figure(args.output, figure_dir, seed)
    trajectories = load_migration_trajectories(args.output)
    centroid_trajectory_figure(trajectories, figure_dir, seed)
    migration_table, migration_summaries = paired_migration_figure(
        trajectories, figure_dir
    )
    morphology_table, morphology_summaries = checkpoint_morphology_figure(
        args.output, figure_dir
    )
    write_rows(figure_dir / "migration_metrics.csv", migration_table)
    write_rows(figure_dir / "checkpoint_morphology.csv", morphology_table)
    paired_summaries = paired_results_figure(rows, figure_dir)
    state_counts = state_outcomes_figure(rows, figure_dir)
    manifest = {
        "source": str(args.output),
        "figures": sorted(path.name for path in figure_dir.glob("*.png")),
        "derived_tables": [
            "checkpoint_morphology.csv",
            "migration_metrics.csv",
        ],
        "representative_seed": seed,
        "representative_seed_selection_rule": selection_rule,
        "representative_checkpoint_times": checkpoint_times,
        "storage": storage,
        "paired_metric_summaries": paired_summaries,
        "paired_migration_summaries": migration_summaries,
        "migration_metric_definitions": {
            "net_displacement_per_radius": "Endpoint distance from the initial centroid divided by initial target radius.",
            "path_length_per_radius": "Sum of consecutive centroid step lengths divided by initial target radius.",
            "maximum_excursion_per_radius": "Largest distance from the initial centroid divided by initial target radius.",
            "straightness": "Net displacement divided by path length.",
        },
        "checkpoint_morphology_summaries": morphology_summaries,
        "morphology_metric_definitions": {
            "relative_area": "Integral of diffuse phi divided by target area.",
            "relative_perimeter": "Integral of gradient magnitude of phi divided by the target-circle perimeter.",
            "shape_anisotropy": "Difference over sum of eigenvalues of the phi-weighted second-moment tensor.",
            "boundary_radial_roughness": "Coefficient of variation of distances from the phi>=0.5 boundary to the phi-weighted centroid.",
        },
        "terminal_state_counts": state_counts,
        "report_integration": "none",
    }
    (figure_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
