#!/usr/bin/env python3
"""Repeat the 1x/WGD Chiou comparison with mechanics released after stimulation."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1] / "migration1d" / "src"))

from migration1d import (  # noqa: E402
    ChiouChemistry,
    ChiouParameters,
    FieldForce,
    MechanicsParameters,
    MigrationModel,
    SimulationConfig,
    TwoEdgeStimulus,
    run_simulation,
)
from run_chiou_endpoint_wgd_kymographs import (  # noqa: E402
    BEST,
    LINEAR_WGD_SCALE,
    competitor_series,
    resolution_time,
)

OUTPUT = ROOT / "results" / "migration1d_chiou_endpoint_wgd_motile"
FORCE_SCALE = 0.01
FORCE_HALF = 10.0
FORCE_HILL = 4.0


def force_response(value: float) -> float:
    positive = max(value, 0.0)
    power = positive**FORCE_HILL
    return power / (FORCE_HALF**FORCE_HILL + power)


def make_model(length: float, linear_scale: float) -> MigrationModel:
    diffusion_scale = BEST["diffusion_scale"]
    chemistry = ChiouChemistry(
        parameters=ChiouParameters(
            activator_diffusivity=0.01 * diffusion_scale,
            substrate_diffusivity=diffusion_scale,
            initial_substrate=BEST["initial_substrate"],
        ),
        stimulus=TwoEdgeStimulus(
            left_amplitude=BEST["left_stimulus"],
            right_amplitude=BEST["left_stimulus"] * BEST["right_left_ratio"],
            width_fraction=BEST["stimulus_width_fraction"],
            start=50.0,
            stop=70.0,
        ),
    )
    endpoint_force = FieldForce(
        "chiou_activator",
        response=force_response,
        gain=FORCE_SCALE,
        direction="outward",
    )
    return MigrationModel(
        mechanics=MechanicsParameters(
            reference_length=length,
            # Preserve E=kL and mu=eta*L for the unit-area 1D element.
            spring_stiffness=1.0 / linear_scale,
            dashpot_viscosity=2.0 / linear_scale,
            left_friction=1.0,
            right_friction=1.0,
        ),
        chemistry=chemistry,
        left_force=endpoint_force,
        right_force=endpoint_force,
    )


def run_condition(length: float, nodes: int, linear_scale: float):
    return run_simulation(
        make_model(length, linear_scale),
        SimulationConfig(
            nodes=nodes,
            duration=1000.0,
            output_interval=2.0,
            max_step=0.5,
            noise_amplitude=0.0,
            mechanics_activation_time=70.0,
        ),
    )


def render(results: dict, metrics: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"1x": "#0072B2", "wgd": "#D55E00"}
    titles = {
        "1x": "1× motile",
        "wgd": "WGD-null motile",
    }
    vmax = max(float(result.membrane[:, 0].max()) for result in results.values())
    xmin = min(float(result.rear.min()) for result in results.values())
    xmax = max(float(result.front.max()) for result in results.values())
    margin = 0.025 * (xmax - xmin)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.5), sharex=True, sharey=True,
                             constrained_layout=True)
    image = None
    for ax, key in zip(axes, ("1x", "wgd")):
        result = results[key]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        image = ax.pcolormesh(
            result.position, time_grid, result.membrane[:, 0],
            shading="auto", cmap="magma", vmin=0.0, vmax=vmax, rasterized=True,
        )
        ax.plot(result.rear, result.times, color="white", linewidth=0.7)
        ax.plot(result.front, result.times, color="white", linewidth=0.7)
        ax.axhline(50.0, color="#54d6ff", linestyle="--", linewidth=0.8)
        ax.axhline(70.0, color="#54d6ff", linestyle="--", linewidth=0.8)
        resolved = resolution_time(result, metrics[key])
        resolution_label = (
            "not resolved by 1000 s" if resolved is None
            else f"resolved at t = {resolved:g}"
        )
        ax.text(
            0.98, 0.98,
            f"{resolution_label}\nΔx = {result.center[-1] - result.center[0]:+.2f}",
            transform=ax.transAxes, ha="right", va="top", color="white", fontsize=9,
            bbox={"facecolor": "black", "alpha": 0.4, "edgecolor": "none"},
        )
        ax.set(xlabel="laboratory position", title=titles[key],
               xlim=(xmin - margin, xmax + margin))
    axes[0].set_ylabel("time")
    colorbar = fig.colorbar(image, ax=axes, pad=0.02)
    colorbar.set_label("Chiou activator")
    fig.suptitle(
        "Endpoint competition with protrusive mechanics released at 70 s\n"
        "white curves: material boundaries · cyan dashes: cue window",
        fontsize=14,
    )
    fig.savefig(OUTPUT / "00_chiou_1x_vs_wgd_motile_kymographs.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 9.0), sharex=True,
                             constrained_layout=True)
    for key in ("1x", "wgd"):
        result = results[key]
        axes[0].plot(result.times, result.center - result.center[0],
                     color=colors[key], label=key.upper())
        axes[1].plot(result.times, 100.0 * (result.length / result.length[0] - 1.0),
                     color=colors[key], label=key.upper())
        axes[2].plot(result.times, metrics[key]["mass_ratio"],
                     color=colors[key], label=key.upper())
    axes[0].set(ylabel="center displacement")
    axes[0].legend(frameon=False)
    axes[1].set(ylabel="length change (%)")
    axes[2].axhline(0.1, color="0.35", linestyle=":",
                    label="resolution threshold")
    axes[2].set(xlabel="time", ylabel="loser/winner peak mass", ylim=(-0.02, 1.02))
    axes[2].legend(frameon=False)
    fig.suptitle("Motility, deformation, and chemical competition")
    fig.savefig(OUTPUT / "01_motility_and_competition_metrics.png", dpi=220)
    plt.close(fig)

    selected_times = (70.0, 200.0, 500.0, 1000.0)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, time in zip(axes.flat, selected_times):
        for key in ("1x", "wgd"):
            result = results[key]
            index = int(np.argmin(np.abs(result.times - time)))
            relative = (
                (result.position[index] - result.rear[index]) / result.length[index]
            )
            ax.plot(relative, result.membrane[index, 0], color=colors[key],
                    label=key.upper())
        ax.set(title=f"t = {time:g}", xlabel="relative position", ylabel="activator")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Motile profiles in relative cell coordinates")
    fig.savefig(OUTPUT / "02_relative_coordinate_profiles.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    conditions = {
        "1x": {"length": 10.0, "nodes": 101, "linear_scale": 1.0},
        "wgd": {
            "length": 10.0 * LINEAR_WGD_SCALE,
            "nodes": 127,
            "linear_scale": LINEAR_WGD_SCALE,
        },
    }
    results = {
        key: run_condition(**values) for key, values in conditions.items()
    }
    metrics = {key: competitor_series(result) for key, result in results.items()}
    for key, result in results.items():
        np.savez_compressed(
            OUTPUT / f"{key}_motile_trajectory.npz",
            times=result.times,
            position=result.position,
            velocity=result.velocity,
            membrane=result.membrane,
            cytoplasm=result.cytoplasm,
            applied_edge_force=result.applied_edge_force,
            spring_force=result.spring_force,
            dashpot_force=result.dashpot_force,
            **metrics[key],
        )
    render(results, metrics)

    summary = {
        "experiment": "best Chiou 1x/WGD comparison with mechanics released at 70 s",
        "mechanical_coupling": {
            "kind": "protrusive/outward endpoint force",
            "response": "a^4 / (10^4 + a^4)",
            "force_scale": FORCE_SCALE,
            "mechanics_release_time": 70.0,
            "baseline_spring_stiffness": 1.0,
            "baseline_dashpot_viscosity": 2.0,
            "endpoint_friction": 1.0,
        },
        "wgd_mechanical_null": {
            "interpretation": "preserve unit-area 1D material coefficients E=kL and mu=eta*L",
            "spring_stiffness": "x 2^(-1/3)",
            "dashpot_viscosity": "x 2^(-1/3)",
            "endpoint_friction": "preserved",
            "local_field_to_force_mapping": "preserved",
        },
        "chemical_wgd_null": {
            "length": "x 2^(1/3)",
            "concentrations": "preserved",
            "native_line_mass": "x 2^(1/3)",
            "reaction_and_diffusion_parameters": "preserved",
            "stimulus": "amplitude/duration fixed; relative width fixed",
        },
        "conditions": {},
    }
    for key, result in results.items():
        total_membrane, total_cytoplasm = result.compartment_totals()
        total = total_membrane + total_cytoplasm
        summary["conditions"][key] = {
            **conditions[key],
            "spring_stiffness": result.model.mechanics.spring_stiffness,
            "dashpot_viscosity": result.model.mechanics.dashpot_viscosity,
            "net_displacement": float(result.center[-1] - result.center[0]),
            "minimum_length": float(result.length.min()),
            "maximum_length": float(result.length.max()),
            "final_length": float(result.length[-1]),
            "final_loser_winner_mass_ratio": float(metrics[key]["mass_ratio"][-1]),
            "final_left_peak": float(metrics[key]["left_peak"][-1]),
            "final_right_peak": float(metrics[key]["right_peak"][-1]),
            "resolution_time": resolution_time(result, metrics[key]),
            "maximum_relative_material_drift": float(
                np.max(np.abs(total / total[0] - 1.0))
            ),
        }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
