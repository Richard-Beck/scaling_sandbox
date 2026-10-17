#!/usr/bin/env python3
"""Compare the best Chiou endpoint competition at 1x and WGD null scaling."""

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
    MechanicsParameters,
    MigrationModel,
    SimulationConfig,
    TwoEdgeStimulus,
    run_simulation,
)
from migration1d.model import control_volumes  # noqa: E402

OUTPUT = ROOT / "results" / "migration1d_chiou_endpoint_wgd"
LINEAR_WGD_SCALE = 2.0 ** (1.0 / 3.0)
BEST = {
    "diffusion_scale": 5.059808062389493,
    "initial_substrate": 2.146940352860838,
    "stimulus_width_fraction": 0.06418176913633943,
    "left_stimulus": 0.12696002777665855,
    "right_left_ratio": 0.7909728942997754,
}


def make_model(length: float) -> MigrationModel:
    scale = BEST["diffusion_scale"]
    chemistry = ChiouChemistry(
        parameters=ChiouParameters(
            activator_diffusivity=0.01 * scale,
            substrate_diffusivity=scale,
            initial_substrate=BEST["initial_substrate"],
        ),
        stimulus=TwoEdgeStimulus(
            left_amplitude=BEST["left_stimulus"],
            right_amplitude=(BEST["left_stimulus"] * BEST["right_left_ratio"]),
            width_fraction=BEST["stimulus_width_fraction"],
            start=50.0,
            stop=70.0,
        ),
    )
    return MigrationModel(
        mechanics=MechanicsParameters(reference_length=length),
        chemistry=chemistry,
    )


def run_condition(length: float, nodes: int):
    return run_simulation(
        make_model(length),
        SimulationConfig(
            nodes=nodes,
            duration=1000.0,
            output_interval=2.0,
            max_step=0.5,
            noise_amplitude=0.0,
            mechanics_activation_time=None,
        ),
    )


def competitor_series(result) -> dict[str, np.ndarray]:
    left_mass = np.empty(len(result.times))
    right_mass = np.empty_like(left_mass)
    left_peak = np.empty_like(left_mass)
    right_peak = np.empty_like(left_mass)
    for index, (position, raw_profile) in enumerate(
        zip(result.position, result.membrane[:, 0])
    ):
        profile = np.maximum(raw_profile, 0.0)
        weights = control_volumes(position)
        midpoint = len(weights) // 2
        left_weights = weights.copy()
        right_weights = weights.copy()
        left_weights[midpoint + 1 :] = 0.0
        right_weights[:midpoint] = 0.0
        left_weights[midpoint] *= 0.5
        right_weights[midpoint] *= 0.5
        excess = np.maximum(profile - np.min(profile), 0.0)
        left_mass[index] = excess @ left_weights
        right_mass[index] = excess @ right_weights
        band = max(2, len(profile) // 5)
        left_peak[index] = np.max(profile[:band])
        right_peak[index] = np.max(profile[-band:])
    mass_min = np.minimum(left_mass, right_mass)
    mass_max = np.maximum(left_mass, right_mass)
    peak_min = np.minimum(left_peak, right_peak)
    peak_max = np.maximum(left_peak, right_peak)
    return {
        "left_mass": left_mass,
        "right_mass": right_mass,
        "left_peak": left_peak,
        "right_peak": right_peak,
        "mass_ratio": np.divide(
            mass_min, mass_max, out=np.zeros_like(mass_min), where=mass_max > 1.0e-14
        ),
        "peak_ratio": np.divide(
            peak_min, peak_max, out=np.zeros_like(peak_min), where=peak_max > 1.0e-14
        ),
        "left_mass_fraction": np.divide(
            left_mass, left_mass + right_mass, out=np.full_like(left_mass, 0.5),
            where=(left_mass + right_mass) > 1.0e-14,
        ),
    }


def resolution_time(result, metrics) -> float | None:
    release = int(np.argmin(np.abs(result.times - 70.0)))
    release_loser_fraction = min(
        metrics["left_mass_fraction"][release],
        1.0 - metrics["left_mass_fraction"][release],
    )
    robust = (
        min(metrics["left_peak"][release], metrics["right_peak"][release]) >= 3.0
        and release_loser_fraction >= 0.15
    )
    if not robust:
        return None
    resolved = (
        (result.times >= 70.0)
        & (metrics["mass_ratio"] <= 0.10)
        & (metrics["peak_ratio"] <= 0.20)
        & (np.maximum(metrics["left_peak"], metrics["right_peak"]) >= 3.0)
    )
    indices = np.flatnonzero(resolved)
    return float(result.times[indices[0]]) if len(indices) else None


def render(results: dict, metrics: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vmax = max(float(result.membrane[:, 0].max()) for result in results.values())
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), sharey=True,
                             constrained_layout=True)
    image = None
    titles = {
        "1x": "1× baseline · L = 10.00",
        "wgd": f"WGD null · L = {10 * LINEAR_WGD_SCALE:.2f}",
    }
    for ax, key in zip(axes, ("1x", "wgd")):
        result = results[key]
        time_grid = np.broadcast_to(result.times[:, None], result.position.shape)
        image = ax.pcolormesh(
            result.position, time_grid, result.membrane[:, 0],
            shading="auto", cmap="magma", vmin=0.0, vmax=vmax, rasterized=True,
        )
        ax.axhline(50.0, color="#54d6ff", linestyle="--", linewidth=0.8)
        ax.axhline(70.0, color="#54d6ff", linestyle="--", linewidth=0.8)
        resolved = resolution_time(result, metrics[key])
        annotation = "not resolved by 1000 s" if resolved is None else f"resolved at t = {resolved:g}"
        ax.text(
            0.98, 0.98, annotation, transform=ax.transAxes, ha="right", va="top",
            color="white", fontsize=9,
            bbox={"facecolor": "black", "alpha": 0.4, "edgecolor": "none"},
        )
        ax.set(xlabel="physical position", title=titles[key])
    axes[0].set_ylabel("time")
    colorbar = fig.colorbar(image, ax=axes, pad=0.02)
    colorbar.set_label("Chiou activator")
    fig.suptitle(
        "Endpoint competition to 1000 s\n"
        "cyan dashes: 50–70 s cues · mechanics off",
        fontsize=14,
    )
    fig.savefig(OUTPUT / "00_chiou_1x_vs_wgd_physical_kymographs.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(9, 7.2), sharex=True,
                             constrained_layout=True)
    colors = {"1x": "#0072B2", "wgd": "#D55E00"}
    for key in ("1x", "wgd"):
        axes[0].plot(results[key].times, metrics[key]["mass_ratio"],
                     color=colors[key], label=key.upper())
        axes[1].plot(results[key].times, metrics[key]["left_mass_fraction"],
                     color=colors[key], label=key.upper())
    axes[0].axhline(0.1, color="0.3", linestyle=":", label="resolution threshold")
    axes[0].set(ylabel="loser/winner peak mass", ylim=(-0.02, 1.02))
    axes[0].legend(frameon=False, ncol=3)
    axes[1].axhline(0.5, color="0.5", linestyle=":")
    axes[1].set(xlabel="time", ylabel="left peak mass fraction", ylim=(0.45, 1.02))
    axes[1].legend(frameon=False)
    fig.suptitle("WGD preserves two competing endpoint peaks under the default null")
    fig.savefig(OUTPUT / "01_competition_metrics.png", dpi=220)
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
    fig.suptitle("Profiles compared in relative cell coordinates")
    fig.savefig(OUTPUT / "02_relative_coordinate_profiles.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    conditions = {
        "1x": {"length": 10.0, "nodes": 101},
        "wgd": {"length": 10.0 * LINEAR_WGD_SCALE, "nodes": 127},
    }
    results = {
        key: run_condition(values["length"], values["nodes"])
        for key, values in conditions.items()
    }
    metrics = {key: competitor_series(result) for key, result in results.items()}
    for key, result in results.items():
        np.savez_compressed(
            OUTPUT / f"{key}_trajectory.npz",
            times=result.times,
            position=result.position,
            membrane=result.membrane,
            cytoplasm=result.cytoplasm,
            **metrics[key],
        )
    render(results, metrics)

    summary = {
        "experiment": "best validated Chiou endpoint competition to 1000 s",
        "mechanics": "off",
        "wgd_null_interpretation": {
            "volume_ratio": 2.0,
            "linear_scale": LINEAR_WGD_SCALE,
            "length": "x 2^(1/3)",
            "initial_concentrations": "preserved",
            "native_line_mass": "x 2^(1/3)",
            "reaction_rates": "preserved",
            "physical_diffusion_coefficients": "preserved",
            "stimulus_amplitude_and_duration": "preserved",
            "stimulus_width": "preserved as fraction of cell length",
            "spatial_resolution": "approximately preserved",
        },
        "best_parameters": BEST,
        "conditions": {},
    }
    for key, result in results.items():
        total_membrane, total_cytoplasm = result.compartment_totals()
        total = total_membrane + total_cytoplasm
        release = int(np.argmin(np.abs(result.times - 70.0)))
        final = -1
        summary["conditions"][key] = {
            **conditions[key],
            "dx": float(result.position[0, 1] - result.position[0, 0]),
            "initial_native_line_mass": float(total[0, 0]),
            "maximum_relative_material_drift": float(
                np.max(np.abs(total / total[0] - 1.0))
            ),
            "release_loser_winner_mass_ratio": float(metrics[key]["mass_ratio"][release]),
            "final_loser_winner_mass_ratio": float(metrics[key]["mass_ratio"][final]),
            "release_left_mass_fraction": float(metrics[key]["left_mass_fraction"][release]),
            "final_left_mass_fraction": float(metrics[key]["left_mass_fraction"][final]),
            "final_left_peak": float(metrics[key]["left_peak"][final]),
            "final_right_peak": float(metrics[key]["right_peak"][final]),
            "resolution_time": resolution_time(result, metrics[key]),
        }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
