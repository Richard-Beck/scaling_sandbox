"""Outputs and figures for the explicit Holmes bulk--surface model."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .core import SimulationResult
from .holmes import FAMILIES


COLORS = ("#1b9e77", "#7570b3", "#d95f02")


def extract_holmes_profiles(result: SimulationResult,
                            stimulus: np.ndarray) -> dict[str, np.ndarray]:
    traces = np.empty((result.times.size, 3, result.geometry.n_surface))
    for index in range(result.times.size):
        traces[index] = result.geometry.trace_bulk(result.bulk[index])
    load = result.geometry.area / result.geometry.perimeter
    projected_inactive = result.surface[:, 3:, :] + load * traces
    return {
        "times": result.times,
        "arc_length": result.geometry.membrane_fraction * result.geometry.perimeter,
        "membrane_points": result.geometry.membrane_points,
        "active_membrane": result.surface[:, :3, :],
        "inactive_membrane": result.surface[:, 3:, :],
        "inactive_cytosol_trace": traces,
        "inactive_projected": projected_inactive,
        "stimulus": np.asarray(stimulus),
        "families": np.asarray(FAMILIES),
    }


def _front_back(values: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    front = int(np.argmax(points[:, 0]))
    back = int(np.argmin(points[:, 0]))
    return values[..., front], values[..., back]


def _diagnostics(result: SimulationResult, profiles: dict[str, np.ndarray],
                 condition: str, equilibrium_fraction: float,
                 association_rate: float) -> dict[str, object]:
    active = profiles["active_membrane"]
    front, back = _front_back(active, result.geometry.membrane_points)
    denominator = np.maximum(np.abs(front) + np.abs(back), 1e-12)
    normalized = (front - back) / denominator
    partition = []
    for family in range(3):
        membrane_amount = (
            result.surface[:, family + 3, :] @ result.geometry.surface_weights
        )
        bulk_amount = result.bulk[:, family, :] @ result.geometry.bulk_weights
        partition.append(membrane_amount / np.maximum(membrane_amount + bulk_amount, 1e-12))
    partition_array = np.asarray(partition).T
    conservation = {}
    for name, values in result.pool_totals.items():
        conservation[name] = float(np.max(np.abs(values / values[0] - 1.0)))
    rac_contrast = normalized[:, 1]
    threshold_indices = np.flatnonzero(np.abs(rac_contrast) >= 0.2)
    time_to_polarity = (
        float(result.times[threshold_indices[0]]) if threshold_indices.size else None
    )
    return {
        "condition": condition,
        "area": result.geometry.area,
        "perimeter": result.geometry.perimeter,
        "bulk_shape": list(result.geometry.bulk_shape),
        "n_surface": result.geometry.n_surface,
        "association_rate": association_rate,
        "equilibrium_inactive_membrane_fraction": equilibrium_fraction,
        "maximum_relative_pool_errors": conservation,
        "minimum_bulk_value": float(np.min(result.bulk)),
        "minimum_surface_value": float(np.min(result.surface)),
        "maximum_bulk_value": float(np.max(result.bulk)),
        "maximum_surface_value": float(np.max(result.surface)),
        "terminal_front_back_contrast": {
            family: float(normalized[-1, index])
            for index, family in enumerate(FAMILIES)
        },
        "terminal_active_peak_x": {
            family: float(result.geometry.membrane_points[
                int(np.argmax(active[-1, index])), 0
            ]) for index, family in enumerate(FAMILIES)
        },
        "time_to_abs_rac_contrast_0.2": time_to_polarity,
        "initial_inactive_membrane_fractions": {
            family: float(partition_array[0, index])
            for index, family in enumerate(FAMILIES)
        },
        "terminal_inactive_membrane_fractions": {
            family: float(partition_array[-1, index])
            for index, family in enumerate(FAMILIES)
        },
        "state_class": (
            "front_polarized" if rac_contrast[-1] >= 0.2 else
            "back_polarized" if rac_contrast[-1] <= -0.2 else
            "weak_or_unpolarized"
        ),
    }


def _plot_kymographs(result: SimulationResult, profiles: dict[str, np.ndarray],
                     output: Path, condition: str) -> None:
    active = profiles["active_membrane"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True,
                             sharey=True)
    for index, (family, color) in enumerate(zip(FAMILIES, COLORS)):
        image = axes[index].imshow(
            active[:, index], origin="lower", aspect="auto", cmap="viridis",
            extent=(0.0, result.geometry.perimeter,
                    result.times[0], result.times[-1]),
        )
        axes[index].set(title=family.capitalize(), xlabel="membrane arc length")
        axes[index].title.set_color(color)
        fig.colorbar(image, ax=axes[index], label="active membrane density")
    axes[0].set_ylabel("time")
    fig.suptitle(f"{condition}: active GTPase membrane dynamics")
    fig.savefig(output / "active_membrane_kymographs.png", dpi=180)
    plt.close(fig)


def _plot_terminal_profiles(result: SimulationResult,
                            profiles: dict[str, np.ndarray], output: Path,
                            condition: str) -> None:
    arc = profiles["arc_length"]
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.0), constrained_layout=True,
                             sharex=True)
    for index, (family, color) in enumerate(zip(FAMILIES, COLORS)):
        axes[0, index].plot(arc, profiles["active_membrane"][-1, index],
                            color=color, lw=2)
        axes[0, index].plot(arc, profiles["inactive_membrane"][-1, index],
                            color="#555555", lw=1.5, ls="--")
        axes[0, index].set_title(family.capitalize(), color=color)
        axes[0, index].set_ylabel("membrane density")
        axes[0, index].legend(("active", "inactive"), frameon=False)
        axes[1, index].plot(arc, profiles["inactive_cytosol_trace"][-1, index],
                            color=color, lw=2)
        axes[1, index].set(xlabel="membrane arc length",
                           ylabel="cytosolic boundary trace")
    fig.suptitle(f"{condition}: terminal Holmes fields")
    fig.savefig(output / "terminal_family_profiles.png", dpi=180)
    plt.close(fig)


def _plot_bulk(result: SimulationResult, output: Path, condition: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0), constrained_layout=True)
    extent = (result.geometry.x[0], result.geometry.x[-1],
              result.geometry.y[0], result.geometry.y[-1])
    for index, (family, color) in enumerate(zip(FAMILIES, COLORS)):
        values = result.geometry.reshape_bulk(result.bulk[-1, index])
        values = np.ma.masked_where(result.geometry.phi < 0.05, values)
        image = axes[index].imshow(values, origin="lower", extent=extent,
                                   cmap="viridis", interpolation="nearest")
        axes[index].set(title=family.capitalize(), xlabel="x", ylabel="y",
                        aspect="equal")
        axes[index].title.set_color(color)
        fig.colorbar(image, ax=axes[index], label="inactive cytosolic concentration")
    fig.suptitle(f"{condition}: terminal inactive cytosolic fields")
    fig.savefig(output / "bulk_inactive_terminal.png", dpi=180)
    plt.close(fig)


def _plot_conservation(result: SimulationResult, output: Path,
                       condition: str) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    for family, color in zip(FAMILIES, COLORS):
        values = result.pool_totals[f"{family}_total"]
        ax.plot(result.times, values / values[0] - 1.0, color=color,
                label=family.capitalize())
    ax.set(xlabel="time", ylabel="relative total error",
           title=f"{condition}: family-wise conservation")
    ax.legend(frameon=False)
    fig.savefig(output / "pool_conservation.png", dpi=180)
    plt.close(fig)


def _plot_partition(result: SimulationResult, output: Path,
                    condition: str) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4), constrained_layout=True)
    for index, (family, color) in enumerate(zip(FAMILIES, COLORS)):
        membrane = result.surface[:, index + 3] @ result.geometry.surface_weights
        bulk = result.bulk[:, index] @ result.geometry.bulk_weights
        ax.plot(result.times, membrane / (membrane + bulk), color=color,
                label=family.capitalize())
    ax.set(xlabel="time", ylabel="inactive fraction on membrane", ylim=(0, 1),
           title=f"{condition}: inactive partition")
    ax.legend(frameon=False)
    fig.savefig(output / "membrane_partition.png", dpi=180)
    plt.close(fig)


def render_holmes_run(result: SimulationResult, output: str | Path, *,
                      condition: str, stimulus: np.ndarray,
                      equilibrium_fraction: float,
                      association_rate: float) -> dict[str, object]:
    directory = Path(output)
    figures = directory / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    profiles = extract_holmes_profiles(result, stimulus)
    np.savez_compressed(directory / "profiles.npz", **profiles)
    diagnostics = _diagnostics(
        result, profiles, condition, equilibrium_fraction, association_rate
    )
    (directory / "diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _plot_kymographs(result, profiles, figures, condition)
    _plot_terminal_profiles(result, profiles, figures, condition)
    _plot_bulk(result, figures, condition)
    _plot_partition(result, figures, condition)
    _plot_conservation(result, figures, condition)
    return diagnostics


def render_holmes_comparison(results: dict[str, SimulationResult],
                             output: str | Path) -> Path:
    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    line_colors = {"1x_direct": "#2166ac", "2x_direct": "#b2182b"}
    for condition, result in results.items():
        points = result.geometry.membrane_points
        order = np.argsort(points[:, 0])
        for index, family in enumerate(FAMILIES):
            axes[index].plot(
                points[order, 0], result.surface[-1, index, order], lw=2,
                color=line_colors.get(condition), label=condition,
            )
            axes[index].set(title=family.capitalize(), xlabel="x position",
                            ylabel="terminal active membrane density")
    for ax in axes:
        ax.legend(frameon=False)
    path = directory / "scaling_response.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path
