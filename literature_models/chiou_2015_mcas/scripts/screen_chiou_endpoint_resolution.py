#!/usr/bin/env python3
"""Targeted fixed-domain screen for resolving two Chiou endpoint peaks by t=500."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np

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

OUTPUT = ROOT / "results" / "migration1d_chiou_endpoint_screen"
DIFFUSION_SCALES = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
SUBSTRATE_LEVELS = (1.6, 1.8, 2.0, 2.2, 2.4)


def simulate(
    *,
    family: str,
    value: float,
    diffusion_scale: float = 1.0,
    substrate: float = 2.0,
    activation: float = 1.0,
    inactivation: float = 1.0,
    saturation: float = 0.0,
    width: float = 0.06,
) -> dict[str, float | str | bool]:
    chemistry = ChiouChemistry(
        parameters=ChiouParameters(
            activation=activation,
            inactivation=inactivation,
            saturation=saturation,
            activator_diffusivity=0.01 * diffusion_scale,
            substrate_diffusivity=1.0 * diffusion_scale,
            initial_substrate=substrate,
        ),
        stimulus=TwoEdgeStimulus(
            left_amplitude=0.085,
            right_amplitude=0.078,
            width_fraction=width,
            start=50.0,
            stop=70.0,
        ),
    )
    result = run_simulation(
        MigrationModel(chemistry=chemistry),
        SimulationConfig(
            nodes=61,
            duration=500.0,
            output_interval=10.0,
            max_step=1.0,
            noise_amplitude=0.0,
            mechanics_activation_time=None,
        ),
    )
    release_index = int(np.argmin(np.abs(result.times - 70.0)))
    weights = control_volumes(result.position[0])
    midpoint = len(weights) // 2
    left_weights = weights.copy()
    right_weights = weights.copy()
    left_weights[midpoint + 1 :] = 0.0
    right_weights[:midpoint] = 0.0
    left_weights[midpoint] *= 0.5
    right_weights[midpoint] *= 0.5

    def metrics(profile):
        profile = np.maximum(profile, 0.0)
        background = float(np.min(profile))
        excess = np.maximum(profile - background, 0.0)
        masses = np.asarray([excess @ left_weights, excess @ right_weights])
        fraction = masses / max(float(masses.sum()), 1.0e-14)
        band = max(2, len(profile) // 5)
        peaks = np.asarray([np.max(profile[:band]), np.max(profile[-band:])])
        return masses, fraction, peaks

    release_mass, release_fraction, release_peaks = metrics(
        result.membrane[release_index, 0]
    )
    final_mass, final_fraction, final_peaks = metrics(result.membrane[-1, 0])
    both_established = bool(np.min(release_peaks) >= 1.0)
    final_peak_ratio = float(np.min(final_peaks) / max(np.max(final_peaks), 1.0e-14))
    final_mass_ratio = float(np.min(final_mass) / max(np.max(final_mass), 1.0e-14))
    if np.max(final_peaks) < 1.0:
        outcome = "both_lost"
    elif both_established and final_peak_ratio <= 0.1:
        outcome = "resolved"
    elif not both_established:
        outcome = "failed_to_establish_two"
    else:
        outcome = "coexist"
    total_membrane, total_cytoplasm = result.compartment_totals()
    total = total_membrane + total_cytoplasm
    return {
        "family": family,
        "value": value,
        "diffusion_scale": diffusion_scale,
        "D_u": 0.01 * diffusion_scale,
        "D_v": diffusion_scale,
        "initial_substrate": substrate,
        "activation": activation,
        "inactivation": inactivation,
        "saturation": saturation,
        "stimulus_width": width,
        "release_left_peak": float(release_peaks[0]),
        "release_right_peak": float(release_peaks[1]),
        "release_left_mass_fraction": float(release_fraction[0]),
        "final_left_peak": float(final_peaks[0]),
        "final_right_peak": float(final_peaks[1]),
        "final_left_mass_fraction": float(final_fraction[0]),
        "final_peak_loser_winner_ratio": final_peak_ratio,
        "final_mass_loser_winner_ratio": final_mass_ratio,
        "both_established_at_release": both_established,
        "outcome": outcome,
        "relative_material_drift": float(
            np.max(np.abs(total / total[0] - 1.0))
        ),
    }


def one_at_a_time() -> list[dict]:
    specifications = {
        "diffusion_scale": DIFFUSION_SCALES,
        "initial_substrate": (1.6, 1.8, 2.0, 2.2, 2.4),
        "activation": (0.6, 0.8, 1.0, 1.2, 1.4),
        "inactivation": (0.6, 0.8, 1.0, 1.2, 1.4),
        "saturation": (0.0, 0.001, 0.005, 0.01, 0.02),
        "stimulus_width": (0.03, 0.04, 0.06, 0.08, 0.10),
    }
    keyword = {
        "diffusion_scale": "diffusion_scale",
        "initial_substrate": "substrate",
        "activation": "activation",
        "inactivation": "inactivation",
        "saturation": "saturation",
        "stimulus_width": "width",
    }
    rows = []
    for family, values in specifications.items():
        for value in values:
            rows.append(simulate(family=family, value=value, **{keyword[family]: value}))
    return rows


def diffusion_mass_grid() -> list[dict]:
    return [
        simulate(
            family="diffusion_x_mass_grid",
            value=diffusion_scale,
            diffusion_scale=diffusion_scale,
            substrate=substrate,
        )
        for substrate in SUBSTRATE_LEVELS
        for diffusion_scale in DIFFUSION_SCALES
    ]


def render(rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    families = (
        "diffusion_scale", "initial_substrate", "activation",
        "inactivation", "saturation", "stimulus_width",
    )
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.0))
    colors = {
        "resolved": "#0072B2",
        "coexist": "#E69F00",
        "failed_to_establish_two": "#CC79A7",
        "both_lost": "#666666",
    }
    for ax, family in zip(axes.flat, families):
        subset = [row for row in rows if row["family"] == family]
        x = np.asarray([row["value"] for row in subset])
        release_ratio = np.asarray([
            min(row["release_left_peak"], row["release_right_peak"])
            / max(row["release_left_peak"], row["release_right_peak"], 1e-14)
            for row in subset
        ])
        final_ratio = np.asarray([
            row["final_peak_loser_winner_ratio"] for row in subset
        ])
        ax.plot(x, release_ratio, "--", color="0.55", label="at release")
        ax.plot(x, final_ratio, "-", color="0.2", label="at 500 s")
        for row in subset:
            ax.scatter(
                row["value"], row["final_peak_loser_winner_ratio"],
                color=colors[row["outcome"]], s=42, zorder=3,
            )
        ax.axhline(0.1, color="#0072B2", linewidth=0.8, linestyle=":")
        ax.set(title=family.replace("_", " "), xlabel="parameter value",
               ylabel="loser/winner peak height")
        ax.set_ylim(-0.03, 1.03)
        if family == "diffusion_scale":
            ax.set_xscale("log", base=2)
        elif family == "saturation":
            ax.set_xscale("symlog", linthresh=0.001)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=color,
                   label=label.replace("_", " "))
        for label, color in colors.items()
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4,
               bbox_to_anchor=(0.5, 0.01))
    fig.suptitle(
        "Fixed-domain Chiou endpoint competition at 500 s\n"
        "Dᵤ and Dᵥ co-scaled; blue threshold denotes operational resolution",
        fontsize=13,
    )
    fig.subplots_adjust(top=0.88, bottom=0.12, hspace=0.36, wspace=0.24)
    fig.savefig(OUTPUT / "00_parameter_family_screen.png", dpi=210)
    plt.close(fig)

    grid = [row for row in rows if row["family"] == "diffusion_x_mass_grid"]
    ratio = np.empty((len(SUBSTRATE_LEVELS), len(DIFFUSION_SCALES)))
    established = np.empty_like(ratio, dtype=bool)
    for i, mass in enumerate(SUBSTRATE_LEVELS):
        for j, scale in enumerate(DIFFUSION_SCALES):
            row = next(
                item for item in grid
                if item["initial_substrate"] == mass and item["diffusion_scale"] == scale
            )
            ratio[i, j] = row["final_peak_loser_winner_ratio"]
            established[i, j] = row["both_established_at_release"]
    fig, ax = plt.subplots(figsize=(9, 5.2), constrained_layout=True)
    image = ax.imshow(ratio, origin="lower", aspect="auto", vmin=0.0, vmax=1.0,
                      cmap="viridis")
    for i in range(ratio.shape[0]):
        for j in range(ratio.shape[1]):
            marker = "" if established[i, j] else "*"
            ax.text(j, i, f"{ratio[i, j]:.2f}{marker}", ha="center", va="center",
                    color="white" if ratio[i, j] < 0.55 else "black", fontsize=8)
    ax.set_xticks(range(len(DIFFUSION_SCALES)), DIFFUSION_SCALES)
    ax.set_yticks(range(len(SUBSTRATE_LEVELS)), SUBSTRATE_LEVELS)
    ax.set(xlabel="common diffusion multiplier", ylabel="initial substrate concentration",
           title="Loser/winner peak height at 500 s (* = two peaks not established at 70 s)")
    fig.colorbar(image, ax=ax, label="loser/winner peak height")
    fig.savefig(OUTPUT / "01_diffusion_mass_targeted_grid.png", dpi=210)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = one_at_a_time() + diffusion_mass_grid()
    fields = list(rows[0])
    with (OUTPUT / "screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    render(rows)

    eligible = [row for row in rows if row["both_established_at_release"]]
    ranked = sorted(eligible, key=lambda row: row["final_peak_loser_winner_ratio"])
    outcomes = {}
    for row in rows:
        outcomes[row["outcome"]] = outcomes.get(row["outcome"], 0) + 1
    summary = {
        "endpoint_time": 500.0,
        "mechanics": "off",
        "baseline": {
            "D_u": 0.01,
            "D_v": 1.0,
            "initial_substrate": 2.0,
            "activation": 1.0,
            "inactivation": 1.0,
            "saturation": 0.0,
            "stimulus_width": 0.06,
            "left_stimulus": 0.085,
            "right_stimulus": 0.078,
            "stimulus_window": [50.0, 70.0],
        },
        "classification": {
            "both_established": "both end-band peak heights >= 1 at 70 s",
            "resolved": "both established and loser/winner peak height <= 0.1 at 500 s",
        },
        "runs": len(rows),
        "outcomes": outcomes,
        "best_established_candidates": ranked[:10],
    }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
