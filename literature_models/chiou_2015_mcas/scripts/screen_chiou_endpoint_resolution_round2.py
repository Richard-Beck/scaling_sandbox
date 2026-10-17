#!/usr/bin/env python3
"""Adaptive second-round search for genuine Chiou endpoint peak resolution."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parents[1] / "migration1d" / "src"))

from migration1d import (  # noqa: E402
    ChiouChemistry,
    ChiouParameters,
    MigrationModel,
    SimulationConfig,
    TwoEdgeStimulus,
    run_simulation,
)
from migration1d.model import control_volumes  # noqa: E402

OUTPUT = ROOT / "results" / "migration1d_chiou_endpoint_screen_round2"
RANGES = {
    "diffusion_scale": (2.0, 7.0),
    "initial_substrate": (2.05, 2.60),
    "stimulus_width": (0.04, 0.09),
    "left_stimulus": (0.07, 0.13),
    "right_left_ratio": (0.65, 0.98),
}


def parameters_from_unit(sample: np.ndarray) -> dict[str, float]:
    return {
        name: lower + float(value) * (upper - lower)
        for value, (name, (lower, upper)) in zip(sample, RANGES.items())
    }


def make_model(parameters: dict[str, float]) -> MigrationModel:
    scale = parameters["diffusion_scale"]
    chemistry = ChiouChemistry(
        parameters=ChiouParameters(
            activator_diffusivity=0.01 * scale,
            substrate_diffusivity=scale,
            initial_substrate=parameters["initial_substrate"],
        ),
        stimulus=TwoEdgeStimulus(
            left_amplitude=parameters["left_stimulus"],
            right_amplitude=(parameters["left_stimulus"]
                             * parameters["right_left_ratio"]),
            width_fraction=parameters["stimulus_width"],
            start=50.0,
            stop=70.0,
        ),
    )
    return MigrationModel(chemistry=chemistry)


def competitor_metrics(profile: np.ndarray, position: np.ndarray) -> dict[str, float]:
    profile = np.maximum(np.asarray(profile), 0.0)
    weights = control_volumes(position)
    midpoint = len(weights) // 2
    left_weights = weights.copy()
    right_weights = weights.copy()
    left_weights[midpoint + 1 :] = 0.0
    right_weights[:midpoint] = 0.0
    left_weights[midpoint] *= 0.5
    right_weights[midpoint] *= 0.5
    background = float(np.min(profile))
    excess = np.maximum(profile - background, 0.0)
    masses = np.asarray([excess @ left_weights, excess @ right_weights])
    fractions = masses / max(float(masses.sum()), 1.0e-14)
    band = max(2, len(profile) // 5)
    peaks = np.asarray([np.max(profile[:band]), np.max(profile[-band:])])
    return {
        "left_peak": float(peaks[0]),
        "right_peak": float(peaks[1]),
        "loser_winner_peak_ratio": float(
            np.min(peaks) / max(np.max(peaks), 1.0e-14)
        ),
        "left_mass_fraction": float(fractions[0]),
        "loser_winner_mass_ratio": float(
            np.min(masses) / max(np.max(masses), 1.0e-14)
        ),
    }


def evaluate(parameters: dict[str, float], *, nodes: int, max_step: float,
             output_interval: float = 10.0):
    result = run_simulation(
        make_model(parameters),
        SimulationConfig(
            nodes=nodes,
            duration=500.0,
            output_interval=output_interval,
            max_step=max_step,
            noise_amplitude=0.0,
            mechanics_activation_time=None,
        ),
    )
    release_index = int(np.argmin(np.abs(result.times - 70.0)))
    release = competitor_metrics(
        result.membrane[release_index, 0], result.position[release_index]
    )
    final = competitor_metrics(result.membrane[-1, 0], result.position[-1])
    release_min_peak = min(release["left_peak"], release["right_peak"])
    release_loser_fraction = min(
        release["left_mass_fraction"], 1.0 - release["left_mass_fraction"]
    )
    final_max_peak = max(final["left_peak"], final["right_peak"])
    robust_two_peaks = release_min_peak >= 3.0 and release_loser_fraction >= 0.15
    one_peak_survives = final_max_peak >= 3.0
    resolved = (
        robust_two_peaks
        and one_peak_survives
        and final["loser_winner_mass_ratio"] <= 0.10
        and final["loser_winner_peak_ratio"] <= 0.20
    )
    penalty = 0.0
    penalty += 8.0 * max(0.0, 3.0 - release_min_peak) / 3.0
    penalty += 8.0 * max(0.0, 0.15 - release_loser_fraction) / 0.15
    penalty += 8.0 * max(0.0, 3.0 - final_max_peak) / 3.0
    score = (
        final["loser_winner_mass_ratio"]
        + 0.35 * final["loser_winner_peak_ratio"]
        + penalty
    )
    row = dict(parameters)
    row.update({
        "right_stimulus": parameters["left_stimulus"] * parameters["right_left_ratio"],
        **{f"release_{key}": value for key, value in release.items()},
        **{f"final_{key}": value for key, value in final.items()},
        "release_min_peak": release_min_peak,
        "release_loser_mass_fraction": release_loser_fraction,
        "final_max_peak": final_max_peak,
        "robust_two_peaks_at_release": robust_two_peaks,
        "one_peak_survives": one_peak_survives,
        "resolved": resolved,
        "score": score,
        "nodes": nodes,
    })
    return row, result


def write_rows(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def render(search_rows: list[dict], validation_rows: list[dict], best_result) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    eligible = np.asarray([row["robust_two_peaks_at_release"] for row in search_rows])
    color = np.asarray([row["final_loser_winner_mass_ratio"] for row in search_rows])
    scatter = axes[0].scatter(
        [row["diffusion_scale"] for row in search_rows],
        [row["initial_substrate"] for row in search_rows],
        c=color, cmap="viridis", vmin=0.0, vmax=1.0,
        marker="o", s=np.where(eligible, 42, 16),
        alpha=np.where(eligible, 0.9, 0.25),
    )
    axes[0].set(
        xlabel="common diffusion multiplier",
        ylabel="initial substrate",
        title="128-point targeted search\nlarge markers: robust peaks at release",
    )
    fig.colorbar(scatter, ax=axes[0], label="final loser/winner peak mass")

    top = sorted(validation_rows, key=lambda row: row["score"])
    labels = [f"V{index + 1}" for index in range(len(top))]
    x = np.arange(len(top))
    axes[1].bar(
        x - 0.18,
        [row["release_loser_winner_mass_ratio"] for row in top],
        width=0.36, label="release", color="0.65",
    )
    axes[1].bar(
        x + 0.18,
        [row["final_loser_winner_mass_ratio"] for row in top],
        width=0.36, label="500 s", color="#0072B2",
    )
    axes[1].axhline(0.1, color="#D55E00", linestyle=":", linewidth=1.0,
                    label="resolution threshold")
    axes[1].set(
        xticks=x, xticklabels=labels,
        ylabel="loser/winner peak mass",
        title="High-resolution validation of top candidates",
        ylim=(0.0, 1.02),
    )
    axes[1].legend(frameon=False)
    fig.suptitle("Second-round Chiou endpoint competition search", fontsize=14)
    fig.savefig(OUTPUT / "00_round2_search_and_validation.png", dpi=220)
    plt.close(fig)

    selected_times = (0.0, 50.0, 70.0, 100.0, 200.0, 500.0)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    time_grid = np.broadcast_to(best_result.times[:, None], best_result.position.shape)
    image = axes[0].pcolormesh(
        best_result.position, time_grid, best_result.membrane[:, 0],
        shading="auto", cmap="magma", rasterized=True,
    )
    axes[0].axhline(50.0, color="#54d6ff", linestyle="--", linewidth=0.8)
    axes[0].axhline(70.0, color="#54d6ff", linestyle="--", linewidth=0.8)
    axes[0].set(xlabel="position", ylabel="time", title="Best candidate kymograph")
    fig.colorbar(image, ax=axes[0], label="Chiou activator")
    for time in selected_times:
        index = int(np.argmin(np.abs(best_result.times - time)))
        axes[1].plot(
            best_result.position[index], best_result.membrane[index, 0],
            label=f"t={best_result.times[index]:g}",
        )
    axes[1].set(xlabel="position", ylabel="Chiou activator",
                title="Best candidate profiles")
    axes[1].legend(ncol=2, fontsize=8, frameon=False)
    fig.savefig(OUTPUT / "01_best_candidate_dynamics.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    samples = qmc.Sobol(d=len(RANGES), scramble=True, seed=17).random_base2(m=7)
    search_rows = []
    for index, sample in enumerate(samples):
        parameters = parameters_from_unit(sample)
        row, _ = evaluate(parameters, nodes=51, max_step=1.0)
        row["search_index"] = index
        search_rows.append(row)
    write_rows(OUTPUT / "search_128.csv", search_rows)

    eligible = [row for row in search_rows if row["robust_two_peaks_at_release"]]
    top_search = sorted(eligible, key=lambda row: row["score"])[:12]
    validation_rows = []
    validation_results = []
    for rank, candidate in enumerate(top_search, start=1):
        parameters = {name: float(candidate[name]) for name in RANGES}
        row, result = evaluate(parameters, nodes=101, max_step=0.5, output_interval=2.0)
        row["validation_rank"] = rank
        row["source_search_index"] = candidate["search_index"]
        validation_rows.append(row)
        validation_results.append(result)
    write_rows(OUTPUT / "validation_top12.csv", validation_rows)

    order = np.argsort([row["score"] for row in validation_rows])
    best_index = int(order[0])
    best = validation_rows[best_index]
    best_result = validation_results[best_index]
    np.savez_compressed(
        OUTPUT / "best_candidate_trajectory.npz",
        times=best_result.times,
        position=best_result.position,
        membrane=best_result.membrane,
        cytoplasm=best_result.cytoplasm,
    )
    render(search_rows, validation_rows, best_result)
    summary = {
        "search_runs": len(search_rows),
        "high_resolution_validations": len(validation_rows),
        "mechanics": "off",
        "endpoint_time": 500.0,
        "robust_release_requirement": {
            "minimum_loser_peak": 3.0,
            "minimum_loser_mass_fraction": 0.15,
        },
        "resolution_requirement": {
            "final_loser_winner_mass_ratio": 0.10,
            "final_loser_winner_peak_ratio": 0.20,
            "minimum_surviving_peak": 3.0,
        },
        "resolved_search_runs": int(sum(row["resolved"] for row in search_rows)),
        "resolved_validations": int(sum(row["resolved"] for row in validation_rows)),
        "best_validated_candidate": best,
        "validated_candidates": sorted(validation_rows, key=lambda row: row["score"]),
    }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
