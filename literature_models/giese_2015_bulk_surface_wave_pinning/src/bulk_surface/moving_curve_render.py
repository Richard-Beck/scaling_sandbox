"""Rendering for a phase-field bulk with an explicit 1D membrane."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

from .moving_curve import MovingCurveResult


def _surface_limits(result: MovingCurveResult) -> tuple[float, float]:
    values = result.surface[:, 0].ravel()
    return float(np.percentile(values, 1)), float(np.percentile(values, 99))


def _colored_curve(axis, points: np.ndarray, values: np.ndarray,
                   limits: tuple[float, float], linewidth: float = 3.0):
    following = np.roll(points, -1, axis=0)
    segments = np.stack((points, following), axis=1)
    edge_values = 0.5 * (values + np.roll(values, -1))
    collection = LineCollection(
        segments, cmap="magma", linewidths=linewidth, capstyle="round",
        norm=plt.Normalize(*limits),
    )
    collection.set_array(edge_values)
    axis.add_collection(collection)
    return collection


def render_moving_curve_snapshots(
    result: MovingCurveResult, output: str | Path
) -> list[Path]:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    extent = (result.x[0], result.x[-1], result.y[0], result.y[-1])
    bulk_values = result.bulk[:, 0][result.phi > 0.1]
    bulk_limits = (
        float(np.percentile(bulk_values, 1)), float(np.percentile(bulk_values, 99))
    )
    surface_limits = _surface_limits(result)
    paths: list[Path] = []
    for index, time in enumerate(result.times):
        phi = result.phi[index]
        bulk = result.bulk[index, 0]
        points = result.membrane_points[index]
        surface = result.surface[index, 0]
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.0), constrained_layout=True)
        image = axes[0, 0].imshow(
            phi, origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=1
        )
        axes[0, 0].plot(*np.vstack((points, points[0])).T, color="tab:cyan", lw=1.3)
        axes[0, 0].set_title(r"Phase field $\phi$ and 1D membrane")
        fig.colorbar(image, ax=axes[0, 0], shrink=0.82)

        masked_bulk = np.ma.masked_where(phi < 0.05, bulk)
        image = axes[0, 1].imshow(
            masked_bulk, origin="lower", extent=extent, cmap="viridis",
            vmin=bulk_limits[0], vmax=bulk_limits[1],
        )
        axes[0, 1].plot(*np.vstack((points, points[0])).T, color="white", lw=0.8)
        axes[0, 1].set_title("Inactive bulk concentration")
        fig.colorbar(image, ax=axes[0, 1], shrink=0.82)

        axes[1, 0].imshow(phi, origin="lower", extent=extent, cmap="Greys",
                          vmin=0, vmax=1, alpha=0.22)
        curve = _colored_curve(axes[1, 0], points, surface, surface_limits)
        axes[1, 0].set_title("Active concentration on explicit 1D membrane")
        fig.colorbar(curve, ax=axes[1, 0], shrink=0.82)

        axes[1, 1].contourf(result.x, result.y, phi, levels=np.linspace(0, 1, 11),
                           cmap="Greys", alpha=0.45)
        axes[1, 1].plot(*np.vstack((points, points[0])).T, color="black", lw=1.0)
        axes[1, 1].plot(result.centroid[:index + 1, 0],
                        result.centroid[:index + 1, 1], color="tab:blue", lw=2)
        axes[1, 1].scatter(*result.centroid[index], color="tab:red", s=35, zorder=3)
        axes[1, 1].set_title(
            f"Centroid track\narea={result.area[index]:.4f}, "
            f"perimeter={result.perimeter[index]:.4f}"
        )
        for axis in axes.ravel():
            axis.set_aspect("equal")
            axis.set_xlim(result.x[0], result.x[-1])
            axis.set_ylim(result.y[0], result.y[-1])
            axis.set_xlabel("x")
            axis.set_ylabel("y")
        fig.suptitle(f"Frozen Mori polarity released at t = 0; t = {time:.1f} s")
        path = output / f"t{time:07.1f}.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        paths.append(path)
    return paths


def render_moving_curve_overview(
    result: MovingCurveResult, output: str | Path,
    times: tuple[float, ...] = (0, 10, 20, 30, 40, 50),
) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    indices = [int(np.argmin(np.abs(result.times - time))) for time in times]
    limits = _surface_limits(result)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.2), constrained_layout=True)
    curve = None
    for axis, index in zip(axes.ravel(), indices):
        phi = result.phi[index]
        points = result.membrane_points[index]
        surface = result.surface[index, 0]
        axis.imshow(
            phi, origin="lower",
            extent=(result.x[0], result.x[-1], result.y[0], result.y[-1]),
            cmap="Greys", vmin=0, vmax=1, alpha=0.22,
        )
        curve = _colored_curve(axis, points, surface, limits, linewidth=3.5)
        axis.plot(result.centroid[:index + 1, 0], result.centroid[:index + 1, 1],
                  color="tab:cyan", lw=1.4)
        axis.set_title(f"t = {result.times[index]:.0f} s")
        axis.set_aspect("equal")
        axis.set_xlim(result.x[0], result.x[-1])
        axis.set_ylim(result.y[0], result.y[-1])
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    if curve is not None:
        fig.colorbar(curve, ax=axes.ravel().tolist(), shrink=0.8,
                     label="active membrane concentration")
    fig.suptitle("Released Mori polarity on an explicit moving 1D membrane")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output
