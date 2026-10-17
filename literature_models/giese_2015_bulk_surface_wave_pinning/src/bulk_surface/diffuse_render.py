"""Rendering utilities for dynamic diffuse-domain trajectories."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .diffuse import DiffuseSimulationResult, geometry_fields


def _limits(result: DiffuseSimulationResult) -> tuple[tuple[float, float],
                                                       tuple[float, float]]:
    bulk_values: list[np.ndarray] = []
    surface_values: list[np.ndarray] = []
    dx = float(result.x[1] - result.x[0])
    for phi, bulk, surface in zip(result.phi, result.bulk[:, 0], result.surface[:, 0]):
        fields = geometry_fields(phi, result.x, result.y, dx)
        bulk_values.append(bulk[phi > 0.1])
        interface = ((fields.delta > 0.15 * fields.delta.max())
                     & (phi > 0.1) & (phi < 0.9))
        surface_values.append(surface[interface])
    bulk_all = np.concatenate(bulk_values)
    surface_all = np.concatenate(surface_values)
    return (
        (float(np.percentile(bulk_all, 1)), float(np.percentile(bulk_all, 99))),
        (float(np.percentile(surface_all, 1)), float(np.percentile(surface_all, 99))),
    )


def render_diffuse_snapshots(
    result: DiffuseSimulationResult, output: str | Path
) -> list[Path]:
    """Render one four-panel PNG for every saved trajectory frame."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    extent = (result.x[0], result.x[-1], result.y[0], result.y[-1])
    dx = float(result.x[1] - result.x[0])
    bulk_limits, surface_limits = _limits(result)
    paths: list[Path] = []
    for index, time in enumerate(result.times):
        phi = result.phi[index]
        bulk = result.bulk[index, 0]
        surface = result.surface[index, 0]
        fields = geometry_fields(phi, result.x, result.y, dx)
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 9.0), constrained_layout=True)
        phase_image = axes[0, 0].imshow(
            phi, origin="lower", extent=extent, cmap="Greys", vmin=0, vmax=1
        )
        axes[0, 0].contour(result.x, result.y, phi, levels=[0.5], colors="tab:cyan",
                           linewidths=1.5)
        axes[0, 0].set_title(r"Phase field $\phi$")
        fig.colorbar(phase_image, ax=axes[0, 0], shrink=0.82)

        bulk_masked = np.ma.masked_where(phi < 0.05, bulk)
        bulk_image = axes[0, 1].imshow(
            bulk_masked, origin="lower", extent=extent, cmap="viridis",
            vmin=bulk_limits[0], vmax=bulk_limits[1],
        )
        axes[0, 1].contour(result.x, result.y, phi, levels=[0.5], colors="white",
                           linewidths=0.8)
        axes[0, 1].set_title("Inactive bulk concentration")
        fig.colorbar(bulk_image, ax=axes[0, 1], shrink=0.82)

        interface = ((fields.delta > 0.12 * fields.delta.max())
                     & (phi > 0.1) & (phi < 0.9))
        surface_masked = np.ma.masked_where(~interface, surface)
        surface_image = axes[1, 0].imshow(
            surface_masked, origin="lower", extent=extent, cmap="magma",
            vmin=surface_limits[0], vmax=surface_limits[1],
        )
        axes[1, 0].contour(result.x, result.y, phi, levels=[0.5], colors="white",
                           linewidths=0.7)
        axes[1, 0].set_title("Active membrane concentration")
        fig.colorbar(surface_image, ax=axes[1, 0], shrink=0.82)

        axes[1, 1].contourf(result.x, result.y, phi, levels=np.linspace(0, 1, 11),
                           cmap="Greys", alpha=0.45)
        axes[1, 1].contour(result.x, result.y, phi, levels=[0.5], colors="black",
                           linewidths=1.2)
        axes[1, 1].plot(result.centroid[:index + 1, 0],
                        result.centroid[:index + 1, 1], color="tab:blue", lw=2)
        axes[1, 1].scatter(result.centroid[index, 0], result.centroid[index, 1],
                           color="tab:red", s=35, zorder=3)
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


def render_release_overview(
    result: DiffuseSimulationResult,
    output: str | Path,
    times: tuple[float, ...] = (0, 10, 20, 30, 40, 50),
) -> Path:
    """Render selected active-membrane snapshots in one comparison figure."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    dx = float(result.x[1] - result.x[0])
    _, surface_limits = _limits(result)
    indices = [int(np.argmin(np.abs(result.times - time))) for time in times]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.2), constrained_layout=True)
    image = None
    for axis, index in zip(axes.ravel(), indices):
        phi = result.phi[index]
        surface = result.surface[index, 0]
        fields = geometry_fields(phi, result.x, result.y, dx)
        interface = ((fields.delta > 0.12 * fields.delta.max())
                     & (phi > 0.1) & (phi < 0.9))
        masked = np.ma.masked_where(~interface, surface)
        image = axis.imshow(
            masked, origin="lower",
            extent=(result.x[0], result.x[-1], result.y[0], result.y[-1]),
            cmap="magma", vmin=surface_limits[0], vmax=surface_limits[1],
        )
        axis.contour(result.x, result.y, phi, levels=[0.5], colors="black",
                     linewidths=1.0)
        axis.plot(result.centroid[:index + 1, 0], result.centroid[:index + 1, 1],
                  color="tab:cyan", lw=1.4)
        axis.set_title(f"t = {result.times[index]:.0f} s")
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.8,
                     label="active membrane concentration")
    fig.suptitle("Release of induced Mori polarization into a deformable cell")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    return output
