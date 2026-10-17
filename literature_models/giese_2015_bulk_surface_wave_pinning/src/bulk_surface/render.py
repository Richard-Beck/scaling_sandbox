"""Plot and profile-output helpers for bulk--surface trajectories."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .core import SimulationResult


def extract_surface_profiles(result: SimulationResult) -> dict[str, np.ndarray]:
    """Extract every saved surface field and every bulk membrane trace."""
    traces = np.empty((result.times.size, result.bulk.shape[1],
                       result.geometry.n_surface))
    for time_index in range(result.times.size):
        traces[time_index] = result.geometry.trace_bulk(result.bulk[time_index])
    return {
        "times": result.times,
        "arc_length": result.geometry.membrane_fraction * result.geometry.perimeter,
        "membrane_points": result.geometry.membrane_points,
        "surface": result.surface,
        "bulk_trace": traces,
        "surface_species": np.asarray(result.surface_species),
        "bulk_species": np.asarray(result.bulk_species),
    }


def _time_index(times: np.ndarray, requested: float) -> int:
    return int(np.argmin(np.abs(times - requested)))


def render_snapshot(result: SimulationResult, profiles: dict[str, np.ndarray],
                    index: int, output: str | Path, *, title: str,
                    bulk_limits: tuple[float, float],
                    surface_limits: tuple[float, float]) -> Path:
    """Render the full bulk, membrane field, and corresponding line profiles."""
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    g = result.geometry
    bulk = g.reshape_bulk(result.bulk[index, 0])
    bulk = np.ma.masked_where(g.phi < 0.05, bulk)
    active = result.surface[index, 0]
    trace = profiles["bulk_trace"][index, 0]
    arc = profiles["arc_length"]

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.2), constrained_layout=True)
    extent = (g.x[0], g.x[-1], g.y[0], g.y[-1])
    image = axes[0].imshow(
        bulk, origin="lower", extent=extent, cmap="viridis",
        vmin=bulk_limits[0], vmax=bulk_limits[1], interpolation="nearest",
    )
    axes[0].plot(g.membrane_points[:, 0], g.membrane_points[:, 1],
                 color="white", lw=0.6, alpha=0.8)
    axes[0].set(title="Inactive bulk", xlabel="x", ylabel="y", aspect="equal")
    fig.colorbar(image, ax=axes[0], shrink=0.82, label="concentration")

    ring = axes[1].scatter(
        g.membrane_points[:, 0], g.membrane_points[:, 1], c=active,
        cmap="magma", vmin=surface_limits[0], vmax=surface_limits[1],
        s=28, linewidths=0,
    )
    axes[1].set(title="Active surface", xlabel="x", ylabel="y", aspect="equal")
    margin = 0.12 * np.max(np.abs(g.membrane_points))
    axes[1].set_xlim(g.membrane_points[:, 0].min()-margin,
                     g.membrane_points[:, 0].max()+margin)
    axes[1].set_ylim(g.membrane_points[:, 1].min()-margin,
                     g.membrane_points[:, 1].max()+margin)
    fig.colorbar(ring, ax=axes[1], shrink=0.82, label="surface density")

    axes[2].plot(arc, active, color="#b2182b", lw=2, label="active surface")
    axes[2].plot(arc, trace, color="#2166ac", lw=2, label="inactive bulk trace")
    axes[2].set(title="Membrane profiles", xlabel="arc length s",
                ylabel="field value", xlim=(0.0, g.perimeter))
    axes[2].legend(frameon=False)
    fig.suptitle(f"{title} — t = {result.times[index]:g}", fontsize=14)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def render_run(result: SimulationResult, output: str | Path, *, condition: str,
               snapshot_times: tuple[float, ...]) -> dict[str, object]:
    """Save profile data and the standard figure suite for one run."""
    directory = Path(output)
    snapshots = directory / "snapshots"
    figures = directory / "figures"
    snapshots.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    profiles = extract_surface_profiles(result)
    np.savez_compressed(directory / "profiles.npz", **profiles)

    bulk_limits = (float(np.min(result.bulk)), float(np.max(result.bulk)))
    surface_limits = (float(np.min(result.surface)), float(np.max(result.surface)))
    rendered_indices: list[int] = []
    for requested in snapshot_times:
        index = _time_index(result.times, requested)
        if index in rendered_indices:
            continue
        rendered_indices.append(index)
        render_snapshot(
            result, profiles, index,
            snapshots / f"t{result.times[index]:010.3f}.png",
            title=f"{condition.capitalize()} shared-chemistry condition",
            bulk_limits=bulk_limits, surface_limits=surface_limits,
        )

    arc = profiles["arc_length"]
    fig, ax = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(rendered_indices)))
    for color, index in zip(colors, rendered_indices):
        ax.plot(arc, result.surface[index, 0], color=color, lw=1.8,
                label=f"t={result.times[index]:g}")
    ax.set(xlabel="arc length s", ylabel="active surface density",
           title=f"{condition.capitalize()} selected surface profiles",
           xlim=(0.0, result.geometry.perimeter))
    ax.legend(ncol=2, frameon=False, fontsize=8)
    fig.savefig(figures / "surface_profiles_selected.png", dpi=180)
    plt.close(fig)

    for values, filename, label, cmap in (
        (result.surface[:, 0], "surface_kymograph.png", "active surface", "magma"),
        (profiles["bulk_trace"][:, 0], "bulk_trace_kymograph.png",
         "inactive bulk trace", "viridis"),
    ):
        fig, ax = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
        image = ax.imshow(
            values, origin="lower", aspect="auto", cmap=cmap,
            extent=(0.0, result.geometry.perimeter,
                    result.times[0], result.times[-1]),
        )
        ax.set(xlabel="arc length s", ylabel="time",
               title=f"{condition.capitalize()} {label}")
        fig.colorbar(image, ax=ax, label=label)
        fig.savefig(figures / filename, dpi=180)
        plt.close(fig)

    pool = result.pool_totals["total_protein"]
    relative_error = pool / pool[0] - 1.0
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    ax.plot(result.times, relative_error, color="#333333", lw=1.5)
    ax.set(xlabel="time", ylabel="relative pool error",
           title=f"{condition.capitalize()} conservation diagnostic")
    fig.savefig(figures / "pool_conservation.png", dpi=180)
    plt.close(fig)

    diagnostics = {
        "condition": condition,
        "snapshot_times": [float(result.times[i]) for i in rendered_indices],
        "initial_total_protein": float(pool[0]),
        "maximum_relative_pool_error": float(np.max(np.abs(relative_error))),
        "bulk_minimum": float(np.min(result.bulk)),
        "bulk_maximum": float(np.max(result.bulk)),
        "surface_minimum": float(np.min(result.surface)),
        "surface_maximum": float(np.max(result.surface)),
        "initial_surface_contrast": float(np.ptp(result.surface[0, 0])),
        "final_surface_contrast": float(np.ptp(result.surface[-1, 0])),
    }
    (directory / "diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return diagnostics


def render_comparison(results: dict[str, SimulationResult], output: str | Path) -> Path:
    """Render terminal surface profiles and time-dependent contrasts."""
    directory = Path(output)
    directory.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), constrained_layout=True)
    colors = {"mori": "#2166ac", "chiou": "#b2182b"}
    for name, result in results.items():
        arc = result.geometry.membrane_fraction * result.geometry.perimeter
        axes[0].plot(arc, result.surface[-1, 0], lw=2,
                     color=colors.get(name), label=name.capitalize())
        contrast = np.ptp(result.surface[:, 0], axis=1)
        axes[1].plot(result.times, contrast, lw=2,
                     color=colors.get(name), label=name.capitalize())
    axes[0].set(xlabel="arc length s", ylabel="active surface density",
                title="Terminal surface profiles")
    axes[1].set(xlabel="time", ylabel="surface max − min",
                title="Surface contrast")
    for ax in axes:
        ax.legend(frameon=False)
    path = directory / "terminal_profiles_and_contrast.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path
