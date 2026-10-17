#!/usr/bin/env python3
"""Plot scalar and directional-vector kymographs from an axisymmetric run."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCALARS = (
    ("rac", "Rac"),
    ("rho", "Rho"),
    ("mcab", "Bound MCA (mcab)"),
    ("mcau", "Unbound MCA (mcau)"),
    ("tension", "Tension"),
)


def field(frame: pd.DataFrame, name: str, times: np.ndarray, arclength: np.ndarray) -> np.ndarray:
    values = frame.pivot(index="time", columns="arclength", values=name)
    return values.reindex(index=times, columns=arclength).to_numpy()


def add_colorbar(fig: plt.Figure, image, ax: plt.Axes, label: str) -> None:
    bar = fig.colorbar(image, ax=ax, pad=0.012, aspect=28)
    bar.set_label(label)


def vector_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    time: np.ndarray,
    arclength: np.ndarray,
    x_component: np.ndarray,
    y_component: np.ndarray,
    title: str,
    radius: float,
) -> None:
    magnitude = np.hypot(x_component, y_component)
    image = ax.pcolormesh(time, arclength, magnitude.T, shading="auto", cmap="inferno")
    add_colorbar(fig, image, ax, "magnitude")

    # Express vectors in the local arc basis.  Increasing arclength runs from
    # the rear toward the front, and the chosen normal points outward.
    theta = arclength / radius
    tangent_x = np.sin(theta)
    tangent_y = np.cos(theta)
    normal_x = -np.cos(theta)
    normal_y = np.sin(theta)
    tangential = x_component * tangent_x[None, :] + y_component * tangent_y[None, :]
    normal = x_component * normal_x[None, :] + y_component * normal_y[None, :]

    # Sparse, normalized arrows show local direction only; magnitude is in color.
    ti = np.unique(np.linspace(0, len(time) - 1, min(14, len(time)), dtype=int))
    si = np.unique(np.linspace(0, len(arclength) - 1, min(15, len(arclength)), dtype=int))
    sampled_magnitude = magnitude[np.ix_(ti, si)]
    cutoff = max(float(np.nanmax(magnitude)) * 1e-4, np.finfo(float).eps)
    with np.errstate(divide="ignore", invalid="ignore"):
        u = tangential[np.ix_(ti, si)] / sampled_magnitude
        v = normal[np.ix_(ti, si)] / sampled_magnitude
    u[sampled_magnitude <= cutoff] = np.nan
    v[sampled_magnitude <= cutoff] = np.nan
    xx, yy = np.meshgrid(time[ti], arclength[si], indexing="ij")
    ax.quiver(
        xx,
        yy,
        u,
        v,
        angles="uv",
        pivot="middle",
        color="#19F9FF",
        edgecolor="#111827",
        linewidth=0.25,
        width=0.0022,
        headwidth=3.4,
        headlength=4.2,
        headaxislength=3.8,
        scale=24,
        zorder=3,
    )
    ax.set_title(f"{title} magnitude + local direction")


def direction_key(ax: plt.Axes) -> None:
    """Use the spare panel to explain arrow-glyph coordinates."""
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.35, 1.15)
    ax.set_aspect("equal")
    ax.set_facecolor("#111827")
    arrow = dict(arrowstyle="-|>", color="#19F9FF", lw=1.2, mutation_scale=13)
    ax.annotate("", xy=(0.85, 0), xytext=(0, 0), arrowprops=arrow)
    ax.annotate("", xy=(-0.85, 0), xytext=(0, 0), arrowprops=arrow)
    ax.annotate("", xy=(0, 0.85), xytext=(0, 0), arrowprops=arrow)
    ax.annotate("", xy=(0, -0.85), xytext=(0, 0), arrowprops=arrow)
    ax.text(0.9, 0, "+t  front", color="#19F9FF", va="center", ha="left")
    ax.text(-0.9, 0, "−t  rear", color="#19F9FF", va="center", ha="right")
    ax.text(0, 0.9, "+n  outward", color="#19F9FF", va="bottom", ha="center")
    ax.text(0, -0.9, "−n  inward", color="#19F9FF", va="top", ha="center")
    ax.text(
        0,
        -1.12,
        "+t follows increasing meridional arclength\n"
        "Arrow glyphs use this local basis, not the plot axes",
        color="white",
        va="top",
        ha="center",
        fontsize=9,
    )
    ax.set_title("Vector direction key: local tangent/normal basis")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path, help="Tidy axisymmetric trajectory CSV")
    parser.add_argument("output", type=Path, help="Output PNG path")
    parser.add_argument("--radius", type=float, default=6.0, help="Reference-arc radius")
    args = parser.parse_args()

    frame = pd.read_csv(args.csv)
    times = np.sort(frame["time"].unique())
    arclength = np.sort(frame["arclength"].unique())
    expected_rows = len(times) * len(arclength)
    if len(frame) != expected_rows:
        raise ValueError(f"Expected a complete rectangular trajectory ({expected_rows} rows), got {len(frame)}")

    fig, axes_grid = plt.subplots(
        2, 4, figsize=(22, 10), constrained_layout=True
    )
    axes = axes_grid.ravel()
    for ax, (name, label) in zip(axes[:5], SCALARS):
        values = field(frame, name, times, arclength)
        image = ax.pcolormesh(times, arclength, values.T, shading="auto", cmap="viridis")
        add_colorbar(fig, image, ax, label)
        ax.set_title(label)

    vector_panel(
        fig,
        axes[5],
        times,
        arclength,
        field(frame, "displacement_x", times, arclength),
        field(frame, "displacement_y", times, arclength),
        "Membrane displacement",
        args.radius,
    )
    vector_panel(
        fig,
        axes[6],
        times,
        arclength,
        field(frame, "velocity_x", times, arclength),
        field(frame, "velocity_y", times, arclength),
        "Cortical velocity",
        args.radius,
    )

    direction_key(axes[7])
    for ax in axes[:7]:
        ax.set_xlabel("time")
        ax.set_ylabel("meridional arclength")
        ax.axvline(50, color="white", linewidth=0.6, linestyle="--", alpha=0.65)
        ax.axvline(150, color="white", linewidth=0.6, linestyle="--", alpha=0.65)
    fig.suptitle(
        "2D axisymmetric mechanochemical model — published example\n"
        "Optogenetic activation steps 50–150 (dashed lines); vector color = magnitude, arrows = local direction",
        fontsize=15,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
