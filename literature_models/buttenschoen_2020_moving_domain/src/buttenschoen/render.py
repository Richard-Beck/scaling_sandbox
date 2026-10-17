"""Kymograph rendering for the moving-domain Buttenschön model."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .model import SimulationResult


def _time_edges(times: np.ndarray) -> np.ndarray:
    if times.size == 1:
        return np.asarray([times[0]-0.5, times[0]+0.5])
    edges = np.empty(times.size + 1)
    edges[1:-1] = 0.5*(times[:-1] + times[1:])
    edges[0] = times[0] - 0.5*(times[1]-times[0])
    edges[-1] = times[-1] + 0.5*(times[-1]-times[-2])
    return edges


def _extrapolated_edges(values: np.ndarray) -> np.ndarray:
    if values.size == 1:
        return np.asarray([values[0], values[0]])
    out = np.empty(values.size + 1)
    out[1:-1] = 0.5*(values[:-1] + values[1:])
    out[0] = values[0] - 0.5*(values[1]-values[0])
    out[-1] = values[-1] + 0.5*(values[-1]-values[-2])
    return out


def kymograph_mesh(result: SimulationResult) -> tuple[np.ndarray, np.ndarray]:
    """Return time and moving physical-position cell edges for pcolormesh."""
    time_edges = _time_edges(result.times)
    centroid_edges = _extrapolated_edges(result.centroid)
    radius_edges = _extrapolated_edges(result.radius)
    xbar_edges = np.linspace(-1.0, 1.0, result.xbar.size + 1)
    positions = centroid_edges[:, None] + radius_edges[:, None]*xbar_edges
    return np.broadcast_to(time_edges[:, None], positions.shape), positions


def plot_kymograph(result: SimulationResult, output: str | Path,
                   title: str | None = None, vmax: float | None = None) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    time_mesh, position_mesh = kymograph_mesh(result)
    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    image = ax.pcolormesh(time_mesh, position_mesh, result.active,
                          shading="flat", cmap="viridis", vmin=0, vmax=vmax)
    ax.plot(result.times, result.left_edge, color="white", lw=0.5, alpha=0.8)
    ax.plot(result.times, result.right_edge, color="white", lw=0.5, alpha=0.8)
    ax.set(xlabel="Time (s)", ylabel="Physical position (µm)",
           title=title or f"{result.mechanics.coupling.title()} GTPase")
    fig.colorbar(image, ax=ax, label="Active GTPase concentration")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_four_regimes(results: dict[str, SimulationResult], output: str | Path) -> Path:
    order = ("contractile_stiff", "contractile_soft",
             "protrusive_stiff", "protrusive_soft")
    missing = set(order)-set(results)
    if missing:
        raise ValueError("missing regimes: " + ", ".join(sorted(missing)))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True,
                             sharex=True)
    images = []
    for ax, name in zip(axes.flat, order):
        result = results[name]
        time_mesh, position_mesh = kymograph_mesh(result)
        image = ax.pcolormesh(time_mesh, position_mesh, result.active,
                              shading="flat", cmap="viridis", vmin=0)
        images.append(image)
        kind, stiffness = name.split("_")
        ax.set_title(
            f"{kind.title()}, {stiffness}\n"
            f"E={result.mechanics.elastic_modulus:g}, "
            f"total={result.initial_condition.total_amount:g}"
        )
        ax.set_ylabel("Physical position (µm)")
        ax.set_xlabel("Time (s)")
        fig.colorbar(image, ax=ax, label="Active GTPase")
    fig.suptitle("Buttenschön spatial mechanochemical regimes", fontsize=15)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_wgd_companion(
    results: dict[str, dict[str, SimulationResult]],
    output: str | Path,
) -> Path:
    """Render Fig. 9 regimes as matched 1x/2x-volume kymograph pairs."""
    order = ("contractile_stiff", "contractile_soft",
             "protrusive_stiff", "protrusive_soft")
    missing = set(order)-set(results)
    if missing:
        raise ValueError("missing regimes: " + ", ".join(sorted(missing)))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 2, figsize=(12.5, 15.0), constrained_layout=True,
                             sharex=True)
    for row, name in enumerate(order):
        pair = results[name]
        if set(pair) != {"1x", "2x"}:
            raise ValueError(f"{name} must contain exactly 1x and 2x results")
        vmax = max(float(np.max(pair[scale].active)) for scale in ("1x", "2x"))
        position_min = min(float(np.min(pair[scale].left_edge))
                           for scale in ("1x", "2x"))
        position_max = max(float(np.max(pair[scale].right_edge))
                           for scale in ("1x", "2x"))
        padding = 0.025*max(position_max-position_min, 1.0)
        images = []
        for column, scale in enumerate(("1x", "2x")):
            result = pair[scale]
            ax = axes[row, column]
            time_mesh, position_mesh = kymograph_mesh(result)
            image = ax.pcolormesh(time_mesh, position_mesh, result.active,
                                  shading="flat", cmap="viridis", vmin=0,
                                  vmax=vmax)
            images.append(image)
            kind, stiffness = name.split("_")
            scale_text = "baseline" if scale == "1x" else "2× volume"
            ax.set_title(
                f"{kind.title()}, {stiffness} · {scale_text}\n"
                f"L₀={result.mechanics.rest_length:.2f} µm, "
                f"E={result.mechanics.elastic_modulus:.1f}, "
                f"Gᵀ={result.initial_condition.total_amount:.2f}"
            )
            ax.set_ylabel("Physical position (µm)")
            ax.set_xlabel("Time (s)")
            ax.set_ylim(position_min-padding, position_max+padding)
        fig.colorbar(images[-1], ax=axes[row, :], label="Active GTPase",
                     shrink=0.92)
    fig.suptitle(
        "Buttenschön Fig. 9 regimes under matched 1× and 2×-volume geometry",
        fontsize=15,
    )
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_soft_protrusive_mechanism(
    pair: dict[str, SimulationResult], output: str | Path,
) -> Path:
    """Show both chemical fields and the global radius dynamics."""
    if set(pair) != {"1x", "2x"}:
        raise ValueError("pair must contain exactly 1x and 2x results")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 2, figsize=(12.5, 11.5), constrained_layout=True)
    position_min = min(float(np.min(pair[scale].left_edge)) for scale in pair)
    position_max = max(float(np.max(pair[scale].right_edge)) for scale in pair)
    padding = 0.025*max(position_max-position_min, 1.0)
    for row, (field_name, label) in enumerate((
        ("active", "Active GTPase, G"),
        ("inactive", "Inactive GTPase, Gᵢ"),
    )):
        vmax = max(float(np.max(getattr(pair[scale], field_name))) for scale in pair)
        images = []
        for column, scale in enumerate(("1x", "2x")):
            result = pair[scale]
            time_mesh, position_mesh = kymograph_mesh(result)
            image = axes[row, column].pcolormesh(
                time_mesh, position_mesh, getattr(result, field_name),
                shading="flat", cmap="viridis", vmin=0, vmax=vmax,
            )
            images.append(image)
            axes[row, column].set_title(f"{label} · {scale}")
            axes[row, column].set_ylabel("Physical position (µm)")
            axes[row, column].set_xlabel("Time (s)")
            axes[row, column].set_ylim(position_min-padding, position_max+padding)
        fig.colorbar(images[-1], ax=axes[row, :], label="Concentration",
                     shrink=0.92)

    colours = {"1x": "#0072B2", "2x": "#D55E00"}
    for scale in ("1x", "2x"):
        result = pair[scale]
        dradius = np.gradient(result.radius, result.times)
        axes[2, 0].plot(result.times, result.radius, color=colours[scale],
                        label=scale, lw=1.5)
        axes[2, 1].plot(result.radius, dradius, color=colours[scale],
                        label=scale, lw=1.3)
        axes[2, 1].scatter(result.radius[::20], dradius[::20],
                           c=result.times[::20], cmap="Greys", s=7, alpha=0.5)
    axes[2, 0].set(title="Global half-length", xlabel="Time (s)",
                   ylabel="R(t) (µm)")
    axes[2, 1].axhline(0, color="grey", lw=0.7, ls="--")
    axes[2, 1].set(title="Mechanical expansion–recoil trajectory",
                   xlabel="R(t) (µm)", ylabel="dR/dt (µm/s)")
    axes[2, 0].legend(frameon=False)
    axes[2, 1].legend(frameon=False)
    fig.suptitle("How chemistry moves the soft protrusive domain", fontsize=15)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output


def plot_length_centroid_comparison(
    results: dict[str, dict[str, SimulationResult]], output: str | Path,
) -> Path:
    """Compare length and centroid for each mechanochemical regime."""
    order = ("contractile_stiff", "contractile_soft",
             "protrusive_stiff", "protrusive_soft")
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 2, figsize=(12.0, 12.5), constrained_layout=True,
                             sharex=True)
    colours = {"1x": "#0072B2", "2x": "#D55E00"}
    for row, name in enumerate(order):
        if name not in results or set(results[name]) != {"1x", "2x"}:
            raise ValueError(f"missing matched pair for {name}")
        for scale in ("1x", "2x"):
            result = results[name][scale]
            axes[row, 0].plot(result.times, result.length,
                              color=colours[scale], label=scale, lw=1.4)
            axes[row, 1].plot(result.times, result.centroid,
                              color=colours[scale], label=scale, lw=1.4)
        kind, stiffness = name.split("_")
        title = f"{kind.title()}, {stiffness}"
        axes[row, 0].set_title(title+" · cell length")
        axes[row, 1].set_title(title+" · centroid")
        axes[row, 0].set_ylabel("L(t) (µm)")
        axes[row, 1].set_ylabel("xᶜ(t) (µm)")
        axes[row, 0].set_xlabel("Time (s)")
        axes[row, 1].set_xlabel("Time (s)")
        axes[row, 0].legend(frameon=False)
        axes[row, 1].legend(frameon=False)
    fig.suptitle("Size and displacement across Buttenschön regimes", fontsize=15)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output
