#!/usr/bin/env python3
"""Render the actual discrete De Belly mechanics as scientific dashboards."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import animation, colors
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from debelly_1d.discrete import (  # noqa: E402
    DiscreteParameters,
    DiscreteProtocol,
    DiscreteTrajectory,
    simulate_discrete,
    smooth_step,
    zero_protocol,
)


OUT = ROOT / "reports" / "assets" / "debelly_1d" / "v1"
BLUE = "#2166ac"
RED = "#b2182b"
GOLD = "#d89400"
PURPLE = "#762a83"
GREY = "#4d4d4d"
CMAP = plt.get_cmap("coolwarm")


def setup_style():
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.dpi": 180,
            "savefig.bbox": "tight",
        }
    )


def pull_protocol(time):
    return smooth_step(time, amplitude=2.0, ramp_time=1.0)


def cortex_drive_protocol(time):
    value, _ = smooth_step(time, amplitude=0.35, ramp_time=1.0)
    return value, 0.0


def active_contraction(time, links):
    amplitude, _ = smooth_step(time, amplitude=0.6, ramp_time=1.0)
    return amplitude * np.exp(-0.5 * ((links - 10.0) / 1.3) ** 2)


def active_contraction_near_clutch(time, links):
    amplitude, _ = smooth_step(time, amplitude=0.6, ramp_time=1.0)
    return amplitude * np.exp(-0.5 * ((links - 16.0) / 1.3) ** 2)


def experiments(times):
    base = DiscreteParameters(nodes=21, cortical_viscosity=5.0, linker_friction=1.0)
    release = np.zeros(base.nodes)
    release[base.nodes // 2] = 0.8
    return {
        "membrane_pull": simulate_discrete(
            base,
            DiscreteProtocol(
                membrane_left=zero_protocol,
                membrane_right=pull_protocol,
            ),
            times,
        ),
        "node_release": simulate_discrete(
            base,
            DiscreteProtocol(),
            times,
            initial_membrane=release,
        ),
        "cortex_drive": simulate_discrete(
            base,
            DiscreteProtocol(cortex_right=cortex_drive_protocol),
            times,
        ),
        "active_contraction": simulate_discrete(
            base,
            DiscreteProtocol(active_stress=active_contraction),
            times,
        ),
    }


def _spring_polyline(x0, x1, y=1.0, coils=5, amplitude=0.055):
    if x1 <= x0:
        return np.array([[x0, y], [x1, y]])
    lead = 0.13 * (x1 - x0)
    core = np.linspace(x0 + lead, x1 - lead, 2 * coils + 1)
    offsets = np.zeros_like(core)
    offsets[1:-1] = amplitude * np.where(np.arange(1, len(core) - 1) % 2, 1, -1)
    xs = np.concatenate([[x0], core, [x1]])
    ys = np.concatenate([[y], y + offsets, [y]])
    return np.column_stack([xs, ys])


def _dashpot(ax, membrane_x, cortex_x, strength, force, force_scale):
    """Draw a standard dashpot symbol along the membrane--cortex connector."""
    p0 = np.array([cortex_x, 0.0])
    p1 = np.array([membrane_x, 1.0])
    direction = p1 - p0
    normal = np.array([-direction[1], direction[0]])
    normal /= np.linalg.norm(normal)
    color = CMAP(0.5 + 0.5 * np.clip(force / force_scale, -1, 1))
    linewidth = 0.5 + 1.7 * np.sqrt(strength)
    a, b, c, d = 0.0, 0.34, 0.66, 1.0
    ax.plot(*np.vstack([p0 + a * direction, p0 + b * direction]).T, color=color, lw=linewidth)
    half = 0.055 + 0.045 * np.sqrt(strength)
    box = np.vstack(
        [
            p0 + b * direction - half * normal,
            p0 + c * direction - half * normal,
            p0 + c * direction + half * normal,
            p0 + b * direction + half * normal,
        ]
    )
    ax.plot(*np.vstack([box, box[0]]).T, color=color, lw=linewidth)
    ax.plot(*np.vstack([p0 + c * direction, p0 + d * direction]).T, color=color, lw=linewidth)


def draw_mechanics(
    ax,
    trajectory,
    frame,
    *,
    title=None,
    show_colorbar=False,
    tension_scale=None,
    force_scale=None,
    friction_scale=None,
):
    reference = trajectory.reference
    membrane_x = reference + trajectory.membrane_displacement[frame]
    cortex_x = reference + trajectory.cortex_displacement[frame]
    tension = trajectory.spring_tension[frame]
    friction_force = trajectory.linker_force[frame]
    tension_scale = max(
        np.max(np.abs(trajectory.spring_tension)) if tension_scale is None else tension_scale,
        1e-8,
    )
    force_scale = max(
        np.max(np.abs(trajectory.linker_force)) if force_scale is None else force_scale,
        1e-8,
    )
    friction_scale = np.max(trajectory.friction) if friction_scale is None else friction_scale
    relative_strength = np.clip(trajectory.friction / friction_scale, 0, 1)

    for index, value in enumerate(tension):
        color = CMAP(0.5 + 0.5 * np.clip(value / tension_scale, -1, 1))
        width = 0.8 + 3.2 * abs(value) / tension_scale
        spring = _spring_polyline(membrane_x[index], membrane_x[index + 1])
        ax.plot(spring[:, 0], spring[:, 1], color=color, lw=width, solid_capstyle="round")
    ax.scatter(membrane_x, np.ones_like(membrane_x), s=32, facecolor="white", edgecolor=GREY, zorder=5)
    ax.plot(cortex_x, np.zeros_like(cortex_x), color="0.35", lw=3, alpha=0.65)
    ax.scatter(cortex_x, np.zeros_like(cortex_x), s=28, marker="s", color="0.25", zorder=5)
    for xm, xc, strength, force in zip(membrane_x, cortex_x, relative_strength, friction_force):
        _dashpot(ax, xm, xc, strength, force, force_scale)

    velocity = trajectory.cortex_velocity[frame]
    velocity_scale = max(np.max(np.abs(trajectory.cortex_velocity)), 1e-8)
    ax.quiver(
        cortex_x,
        np.full_like(cortex_x, -0.16),
        0.65 * velocity / velocity_scale,
        np.zeros_like(velocity),
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.004,
        color=PURPLE,
    )
    ax.text(-0.7, 1.0, "membrane", ha="right", va="center", color=GREY)
    ax.text(-0.7, 0.0, "cortex", ha="right", va="center", color="0.25")
    ax.set_xlim(-2.5, reference[-1] + 3.0)
    ax.set_ylim(-0.4, 1.3)
    ax.axis("off")
    if title:
        ax.set_title(title, loc="left", weight="bold", fontsize=11)
    if show_colorbar:
        sm = plt.cm.ScalarMappable(norm=colors.Normalize(-tension_scale, tension_scale), cmap=CMAP)
        plt.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.03, label="spring tension")


def _symmetric_limit(values, padding=1.12):
    maximum = max(np.max(np.abs(values)), 1e-7)
    return (-padding * maximum, padding * maximum)


def render_dashboard(trajectory, filename, title, subtitle):
    node = np.arange(len(trajectory.reference))
    link = node[:-1] + 0.5
    fig = plt.figure(figsize=(9.0, 7.3))
    grid = fig.add_gridspec(4, 1, height_ratios=[2.15, 1, 1, 1], hspace=0.44)
    axes = [fig.add_subplot(grid[index]) for index in range(4)]
    displacement_limit = _symmetric_limit(trajectory.membrane_displacement)
    tension_limit = _symmetric_limit(trajectory.spring_tension)
    friction_limit = _symmetric_limit(trajectory.linker_force)

    def update(frame):
        for ax in axes:
            ax.clear()
        draw_mechanics(axes[0], trajectory, frame, title=f"{title}   t/τ = {trajectory.time[frame] / 5.0:.2f}")
        axes[0].text(0.01, 0.02, subtitle, transform=axes[0].transAxes, color="0.3", fontsize=8)
        axes[1].plot(node, trajectory.membrane_displacement[frame], "o-", color=BLUE, ms=3)
        axes[1].set(ylabel=r"membrane $X_i-X_i^0$", ylim=displacement_limit)
        axes[2].plot(link, trajectory.spring_tension[frame], "o-", color=RED, ms=3)
        axes[2].set(ylabel=r"spring tension $T_i$", ylim=tension_limit)
        axes[3].plot(node, trajectory.linker_force[frame], "o-", color=GOLD, ms=3)
        axes[3].set(ylabel=r"linker force $\bar\mu(V_i-\dot X_i)$", xlabel="node index", ylim=friction_limit)
        for ax in axes[1:]:
            ax.axhline(0, color="0.75", lw=0.7)
            ax.set_xlim(-0.5, node[-1] + 0.5)
            ax.grid(alpha=0.18)
        return []

    movie = animation.FuncAnimation(fig, update, frames=len(trajectory.time), interval=150, blit=False)
    movie.save(OUT / filename, writer=animation.PillowWriter(fps=7))
    plt.close(fig)


def render_snapshot_sequence(all_experiments):
    selected = [0, 5, 10, 20, 40]
    rows = [
        ("membrane_pull", "A  Pull membrane node and hold"),
        ("node_release", "B  Displace one node and release"),
        ("cortex_drive", "C  Drive the cortex; do not prescribe membrane motion"),
        ("active_contraction", "D  Local active cortical stress"),
    ]
    fig, axes = plt.subplots(len(rows), len(selected), figsize=(15.0, 8.2))
    for row, (key, label) in enumerate(rows):
        trajectory = all_experiments[key]
        for column, frame in enumerate(selected):
            draw_mechanics(axes[row, column], trajectory, frame)
            time_title = f"t/τ = {trajectory.time[frame] / 5.0:g}"
            axes[row, column].set_title(
                f"{label}\n{time_title}" if column == 0 else time_title,
                fontsize=9,
                weight="bold" if column == 0 else "normal",
            )
    fig.suptitle("Snapshots are rendered states of the discrete mechanical equations", weight="bold", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "08_discrete_experiments.png", facecolor="white")
    plt.close(fig)


def render_comparison_animation(cases, filename, title, row_order, column_order):
    tension_scale = max(np.max(np.abs(value.spring_tension[-1])) for value in cases.values())
    force_scale = max(np.max(np.abs(value.linker_force[-1])) for value in cases.values())
    friction_scale = max(np.max(value.friction) for value in cases.values())
    fig, axes = plt.subplots(
        len(row_order),
        len(column_order),
        figsize=(11.2, 6.2),
        facecolor="white",
    )
    fig.subplots_adjust(top=0.88, bottom=0.10, hspace=0.28, wspace=0.12)
    first = next(iter(cases.values()))
    heading = fig.suptitle("", weight="bold", y=0.98)
    fig.text(
        0.5,
        0.025,
        "spring: red extension, blue compression  •  dashpot: thickness = friction, color = linker force  •  purple arrows = cortical velocity",
        ha="center",
        color="0.35",
        fontsize=8,
    )

    def update(frame):
        for row, row_label in enumerate(row_order):
            for column, column_label in enumerate(column_order):
                ax = axes[row, column]
                ax.clear()
                trajectory = cases[(column_label, row_label)]
                draw_mechanics(
                    ax,
                    trajectory,
                    frame,
                    tension_scale=tension_scale,
                    force_scale=force_scale,
                    friction_scale=friction_scale,
                )
                ax.set_title(f"{column_label} — {row_label}", fontsize=10, weight="bold")
        heading.set_text(f"{title}   t/τ = {first.time[frame] / 5.0:.2f}")

    # Matplotlib's PillowWriter may encode multi-axes frames as transparent
    # deltas without a disposal instruction. Some viewers then retain the old
    # axes and create a diagonal smear. Store complete opaque RGB frames and
    # explicitly restore the background before drawing the next frame.
    frames = []
    for frame in range(len(first.time)):
        update(frame)
        fig.canvas.draw()
        rgba = np.asarray(fig.canvas.buffer_rgba())
        frames.append(Image.fromarray(rgba[:, :, :3].copy(), mode="RGB"))
    frames[0].save(
        OUT / filename,
        save_all=True,
        append_images=frames[1:],
        duration=143,
        loop=0,
        disposal=2,
        optimize=True,
    )
    plt.close(fig)


def render_friction_reversal(times, rows, summary):
    cases = {}
    for friction, label in [(0.05, "weak"), (5.0, "strong")]:
        params = DiscreteParameters(nodes=21, cortical_viscosity=5.0, linker_friction=friction)
        cases[("pull", label)] = simulate_discrete(
            params,
            DiscreteProtocol(membrane_left=zero_protocol, membrane_right=pull_protocol),
            times,
        )
        cases[("drive", label)] = simulate_discrete(
            params,
            DiscreteProtocol(cortex_right=cortex_drive_protocol),
            times,
        )

    tension_scale = max(np.max(np.abs(value.spring_tension[-1])) for value in cases.values())
    force_scale = max(np.max(np.abs(value.linker_force[-1])) for value in cases.values())
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.2))
    for column, perturbation in enumerate(["pull", "drive"]):
        for row, label in enumerate(["weak", "strong"]):
            trajectory = cases[(perturbation, label)]
            draw_mechanics(
                axes[row, column],
                trajectory,
                -1,
                tension_scale=tension_scale,
                force_scale=force_scale,
                friction_scale=5.0,
            )
            title = ("membrane pull" if perturbation == "pull" else "cortical drive") + f" — {label} friction"
            axes[row, column].set_title(title, fontsize=10, weight="bold")
    fig.suptitle("Same discrete rules, different loading layer: the friction dependence reverses", weight="bold", y=1.01)
    fig.text(0.5, 0.01, "red/blue springs: extension/compression  •  dashpot thickness: friction  •  dashpot color: linker force  •  arrows: cortical velocity", ha="center", color="0.35", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "09_discrete_friction_reversal.png", facecolor="white")
    plt.close(fig)
    render_comparison_animation(
        cases,
        "dashpot_05_friction_reversal.gif",
        "Friction reversal from the simulated discrete rules",
        ["weak", "strong"],
        ["pull", "drive"],
    )

    weak_pull = cases[("pull", "weak")]
    strong_pull = cases[("pull", "strong")]
    weak_drive = cases[("drive", "weak")]
    strong_drive = cases[("drive", "strong")]
    summary.extend(
        [
            {
                "quantity": "discrete pull distal/proximal weak friction",
                "value": abs(weak_pull.spring_tension[-1, 0] / weak_pull.spring_tension[-1, -1]),
            },
            {
                "quantity": "discrete pull distal/proximal strong friction",
                "value": abs(strong_pull.spring_tension[-1, 0] / strong_pull.spring_tension[-1, -1]),
            },
            {
                "quantity": "discrete cortical-drive tension norm weak friction",
                "value": np.linalg.norm(weak_drive.spring_tension[-1]),
            },
            {
                "quantity": "discrete cortical-drive tension norm strong friction",
                "value": np.linalg.norm(strong_drive.spring_tension[-1]),
            },
        ]
    )
    for (perturbation, label), trajectory in cases.items():
        for frame, time in enumerate(times):
            rows.append(
                {
                    "experiment": "friction_reversal",
                    "case": f"{perturbation}_{label}",
                    "time": time,
                    "tension_norm": np.linalg.norm(trajectory.spring_tension[frame]),
                    "distal_tension": trajectory.spring_tension[frame, 0],
                    "proximal_tension": trajectory.spring_tension[frame, -1],
                }
            )


def render_mca_clutch(times, rows, summary):
    uniform = np.ones(21)
    depleted = uniform.copy()
    depleted[14:19] = 0.1
    cases = {}
    for profile, label in [(uniform, "uniform MCA"), (depleted, "MCA-low clutch")]:
        params = DiscreteParameters(nodes=21, cortical_viscosity=5.0, linker_friction=profile)
        cases[("pull", label)] = simulate_discrete(
            params,
            DiscreteProtocol(membrane_left=zero_protocol, membrane_right=pull_protocol),
            times,
        )
        cases[("active stress", label)] = simulate_discrete(
            params,
            DiscreteProtocol(active_stress=active_contraction_near_clutch),
            times,
        )
    tension_scale = max(np.max(np.abs(value.spring_tension[-1])) for value in cases.values())
    force_scale = max(np.max(np.abs(value.linker_force[-1])) for value in cases.values())
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.2))
    for column, perturbation in enumerate(["pull", "active stress"]):
        for row, label in enumerate(["uniform MCA", "MCA-low clutch"]):
            trajectory = cases[(perturbation, label)]
            draw_mechanics(
                axes[row, column],
                trajectory,
                -1,
                tension_scale=tension_scale,
                force_scale=force_scale,
                friction_scale=1.0,
            )
            axes[row, column].set_title(f"{perturbation} — {label}", fontsize=10, weight="bold")
            if label == "MCA-low clutch":
                axes[row, column].axvspan(13.5, 18.5, color=GOLD, alpha=0.12)
    fig.suptitle(r"Cut the linkers: a prescribed local reduction in $\bar\mu_i=\chi\rho_{b,i}$", weight="bold", y=1.01)
    fig.text(0.5, 0.01, "red/blue springs: extension/compression  •  dashpot thickness: local friction  •  shaded band: MCA-low region", ha="center", color="0.35", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "10_discrete_mca_clutch.png", facecolor="white")
    plt.close(fig)
    render_comparison_animation(
        cases,
        "dashpot_06_mca_clutch.gif",
        "MCA-low linkers behave as a local mechanical clutch",
        ["uniform MCA", "MCA-low clutch"],
        ["pull", "active stress"],
    )

    for (perturbation, label), trajectory in cases.items():
        for frame, time in enumerate(times):
            rows.append(
                {
                    "experiment": "mca_clutch",
                    "case": f"{perturbation}_{label}",
                    "time": time,
                    "tension_norm": np.linalg.norm(trajectory.spring_tension[frame]),
                    "linker_force_norm": np.linalg.norm(trajectory.linker_force[frame]),
                }
            )
    summary.extend(
        [
            {
                "quantity": "MCA clutch active-stress tension ratio depleted/uniform",
                "value": np.linalg.norm(cases[("active stress", "MCA-low clutch")].spring_tension[-1])
                / np.linalg.norm(cases[("active stress", "uniform MCA")].spring_tension[-1]),
            },
            {
                "quantity": "MCA clutch pull distal/proximal transmission depleted",
                "value": abs(
                    cases[("pull", "MCA-low clutch")].spring_tension[-1, 0]
                    / cases[("pull", "MCA-low clutch")].spring_tension[-1, -1]
                ),
            },
            {
                "quantity": "MCA clutch pull distal/proximal transmission uniform",
                "value": abs(
                    cases[("pull", "uniform MCA")].spring_tension[-1, 0]
                    / cases[("pull", "uniform MCA")].spring_tension[-1, -1]
                ),
            },
        ]
    )


def trajectory_rows(name, trajectory):
    rows = []
    for frame, time in enumerate(trajectory.time):
        for node in range(len(trajectory.reference)):
            rows.append(
                {
                    "experiment": name,
                    "time": time,
                    "node": node,
                    "membrane_displacement": trajectory.membrane_displacement[frame, node],
                    "cortex_displacement": trajectory.cortex_displacement[frame, node],
                    "membrane_velocity": trajectory.membrane_velocity[frame, node],
                    "cortex_velocity": trajectory.cortex_velocity[frame, node],
                    "linker_force": trajectory.linker_force[frame, node],
                    "spring_tension_right": trajectory.spring_tension[frame, node]
                    if node < len(trajectory.reference) - 1
                    else np.nan,
                    "friction": trajectory.friction[node],
                }
            )
    return rows


def main():
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    times = np.linspace(0, 10, 41)
    all_experiments = experiments(times)
    dashboard_specs = {
        "membrane_pull": (
            "dashpot_01_membrane_pull.gif",
            "Exogenous: prescribe the end membrane node",
            "Opposite membrane end fixed; cortex endpoints stationary",
        ),
        "node_release": (
            "dashpot_02_node_release.gif",
            "Initial-value experiment: displace one membrane node and release",
            "No displacement or velocity is prescribed after t = 0",
        ),
        "cortex_drive": (
            "dashpot_03_cortex_drive.gif",
            "Endogenous: prescribe cortical velocity at the right end",
            "Membrane motion is entirely simulated",
        ),
        "active_contraction": (
            "dashpot_04_active_contraction.gif",
            "Endogenous: localized active cortical stress",
            "Cortex endpoints are stationary; membrane nodes are not prescribed",
        ),
    }
    for key, (filename, title, subtitle) in dashboard_specs.items():
        render_dashboard(all_experiments[key], filename, title, subtitle)
    render_snapshot_sequence(all_experiments)

    rows = []
    for name, trajectory in all_experiments.items():
        rows.extend(trajectory_rows(name, trajectory))
    summary = []
    render_friction_reversal(times, rows, summary)
    render_mca_clutch(times, rows, summary)
    pd.DataFrame(rows).to_csv(OUT / "dashpot_plot_data.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "dashpot_summary.csv", index=False)
    print(f"Rendered {len(dashboard_specs) + 2} animations and 3 static figures in {OUT}")


if __name__ == "__main__":
    main()
