#!/usr/bin/env python3
"""Run a compact WP1 qualification panel; deliberately not a parameter screen."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly,
    DeBellyParameters,
    GaussianPulse,
    Goryachev,
    Grid1D,
    LEGI,
    LocalizedPulse,
    Otsuji,
    WavePinning,
    noisy_initial_state,
    phenotype_summary,
    simulate,
)


OUTPUT = ROOT / "results" / "qualification"
JILKINE = (WavePinning(), Goryachev(), Otsuji(), LEGI())
SHORT_NAMES = ("Wave-pinning", "Goryachev", "Otsuji", "LEGI")
SIZE_LENGTHS = (1.0, 3.0, 17.0, 25.0, 39.0)
DUAL_PROFILE_TIMES = {
    "jilkine_wave_pinning": (0.0, 200.0, 2000.0),
    "jilkine_goryachev": (0.0, 2000.0, 4000.0),
    "jilkine_otsuji": (0.0, 2000.0, 4000.0),
    "jilkine_legi": (0.0, 20.0, 200.0),
}


def _length_profiles(
    model,
    protocol: LocalizedPulse,
    end_time: float,
    *,
    lengths: tuple[float, ...] = SIZE_LENGTHS,
    noise_amplitude: float | None = None,
    noise_seed: int = 1,
):
    trajectories = []
    for length in lengths:
        cells = max(101, round(length / 0.1))
        initial_state = None
        if noise_amplitude is not None:
            initial_state = noisy_initial_state(
                model,
                Grid1D(length=length, cells=cells),
                amplitude=noise_amplitude,
                seed=noise_seed,
            )
        trajectories.append(
            simulate(
                model,
                protocol=protocol,
                length=length,
                cells=cells,
                times=np.asarray((0.0, end_time)),
                backend="bdf",
                max_step=5.0,
                initial_state=initial_state,
            )
        )
    return trajectories


def _length_profile_summary(trajectories) -> list[dict[str, object]]:
    return [
        {
            "length_um": trajectory.metadata["length_um"],
            "final_time_s": float(trajectory.time[-1]),
            "cells": trajectory.metadata["cells"],
            "dx_um": trajectory.metadata["dx_um"],
            "left_edge_observable": float(trajectory.observable[-1, 0]),
            "right_edge_observable": float(trajectory.observable[-1, -1]),
            "maximum_mirror_error": float(
                np.max(
                    np.abs(
                        trajectory.observable[-1]
                        - trajectory.observable[-1, ::-1]
                    )
                )
            ),
            "conservation_errors": trajectory.metadata.get("conservation_errors"),
        }
        for trajectory in trajectories
    ]


def _figure_guided_length_profiles(model):
    """Recreate the earlier, explicitly model-specific Figure-2f assumptions."""
    trajectories = []
    for length in SIZE_LENGTHS:
        cells = max(101, round(length / 0.1))
        initial_state = None
        if isinstance(model, WavePinning):
            run_model = WavePinning(u0=0.49, v0=1.8801794752073646)
            protocol = LocalizedPulse(
                amplitude=0.04,
                stop=40.0,
                taper=20.0,
                sides=("left", "right"),
                asymmetry=0.02,
            )
        elif isinstance(model, Goryachev):
            run_model = model
            protocol = LocalizedPulse(
                amplitude=1e-4,
                stop=40.0,
                taper=20.0,
                sides=("left", "right"),
            )
            initial_state = noisy_initial_state(
                run_model,
                Grid1D(length=length, cells=cells),
                amplitude=0.1,
                seed=11,
            )
        elif isinstance(model, Otsuji):
            run_model = model
            protocol = LocalizedPulse(
                amplitude=0.04,
                stop=40.0,
                taper=20.0,
                sides=("left", "right"),
                asymmetry=0.02,
            )
        else:
            run_model = LEGI(da=0.1, dr=0.0, r0=0.1)
            protocol = LocalizedPulse(
                amplitude=2e-4,
                stop=2001.0,
                taper=0.0,
                sides=("left", "right"),
            )
        trajectories.append(
            simulate(
                run_model,
                protocol=protocol,
                length=length,
                cells=cells,
                times=np.asarray((0.0, 2000.0)),
                backend="bdf",
                max_step=10.0,
                initial_state=initial_state,
            )
        )
    return trajectories


def qualify_jilkine() -> list[dict[str, object]]:
    """Reproduce the literal two-identical-transient-stimuli experiment."""
    symmetric_protocol = LocalizedPulse(
        amplitude=0.2,
        stop=40.0,
        taper=20.0,
        width_fraction=0.1,
        sides=("left", "right"),
        asymmetry=0.0,
    )
    rows: list[dict[str, object]] = []
    time_profile_runs = []
    size_profile_runs = {
        "common_stimulus_at_2000_s": [],
        "figure_guided_at_2000_s": [],
    }
    for model in JILKINE:
        sample_times = np.asarray(DUAL_PROFILE_TIMES[model.name])
        trajectory = simulate(
            model,
            protocol=symmetric_protocol,
            length=10.0,
            cells=101,
            times=sample_times,
            backend="bdf",
            max_step=5.0,
        )
        time_profile_runs.append(trajectory)
        common_lengths = _length_profiles(model, symmetric_protocol, 2000.0)
        figure_guided_lengths = _figure_guided_length_profiles(model)
        model_size_runs = {
            "common_stimulus_at_2000_s": common_lengths,
            "figure_guided_at_2000_s": figure_guided_lengths,
        }
        for experiment, experiment_runs in model_size_runs.items():
            size_profile_runs[experiment].append(experiment_runs)
        mirror_errors = [
            float(np.max(np.abs(profile - profile[::-1])))
            for profile in trajectory.observable
        ]
        rows.append(
            {
                "model": model.name,
                "protocol": "two_identical_transient_localized_stimuli",
                "length_um": 10.0,
                "sample_times_s": sample_times.tolist(),
                "stimulus": {
                    "S": symmetric_protocol.amplitude,
                    "t1_s": symmetric_protocol.stop - symmetric_protocol.taper,
                    "t2_s": symmetric_protocol.stop,
                    "width_fraction_each_end": symmetric_protocol.width_fraction,
                    "asymmetry": symmetric_protocol.asymmetry,
                },
                "initial_condition": "model homogeneous state",
                "imposed_symmetry_breaker": None,
                "maximum_mirror_error_by_time": mirror_errors,
                "left_edge_observable_by_time": trajectory.observable[:, 0].tolist(),
                "right_edge_observable_by_time": trajectory.observable[:, -1].tolist(),
                "conservation_errors": trajectory.metadata.get("conservation_errors"),
                "length_scaling_experiments": {
                    experiment: _length_profile_summary(experiment_runs)
                    for experiment, experiment_runs in model_size_runs.items()
                },
                "figure_guided_assumption": {
                    "jilkine_wave_pinning": (
                        "u0=0.49, v0=1.8801794752073646; transient S=0.04; "
                        "2% cue asymmetry"
                    ),
                    "jilkine_goryachev": (
                        "equal transient S=1e-4; conservative initial noise "
                        "A=0.1, seed=11"
                    ),
                    "jilkine_otsuji": "transient S=0.04; 2% cue asymmetry",
                    "jilkine_legi": (
                        "DA=0.1, DR=0, R0=0.1; persistent equal S=2e-4"
                    ),
                }[model.name],
                "solver": {
                    key: trajectory.metadata[key]
                    for key in (
                        "backend",
                        "wall_seconds",
                        "rhs_evaluations",
                        "jacobian_evaluations",
                        "linear_decompositions",
                        "maximum_step",
                    )
                },
            }
        )

    figure, axes = plt.subplots(3, 4, figsize=(14, 9.5))
    colors = ("#4D4D4D", "#D97706", "#2563A6")
    size_colors = {
        length: color
        for length, color in zip(
            SIZE_LENGTHS,
            plt.cm.viridis(np.linspace(0.08, 0.92, len(SIZE_LENGTHS))),
        )
    }
    size_experiments = tuple(size_profile_runs)
    for column, (model, label, trajectory) in enumerate(
        zip(JILKINE, SHORT_NAMES, time_profile_runs)
    ):
        axis = axes[0, column]
        for time, profile, color in zip(trajectory.time, trajectory.observable, colors):
            axis.plot(
                trajectory.x,
                profile,
                color=color,
                lw=2.0,
                label=f"t={time:g} s",
            )
        axis.axvspan(0.0, 1.0, color="#B8C4CE", alpha=0.18, lw=0)
        axis.axvspan(9.0, 10.0, color="#B8C4CE", alpha=0.18, lw=0)
        axis.set_title(label)
        axis.set_xlim(0.0, 10.0)
        axis.set_xlabel("Position x (µm)")
        axis.set_ylabel("R" if isinstance(model, LEGI) else "Active u (µM)")
        axis.legend(frameon=False, fontsize=8)

        size_titles = (
            "T=2000 s; common transient cues",
            "T=2000 s; figure-guided assumptions",
        )
        for row, (experiment, title) in enumerate(
            zip(size_experiments, size_titles), start=1
        ):
            axis = axes[row, column]
            length_trajectories = size_profile_runs[experiment][column]
            for length_trajectory in length_trajectories:
                length = length_trajectory.metadata["length_um"]
                axis.plot(
                    length_trajectory.x / length,
                    length_trajectory.observable[-1],
                    color=size_colors[length],
                    lw=1.8,
                    label=f"L={length:g} µm",
                )
            axis.axvspan(0.0, 0.1, color="#B8C4CE", alpha=0.18, lw=0)
            axis.axvspan(0.9, 1.0, color="#B8C4CE", alpha=0.18, lw=0)
            axis.set_title(title, fontsize=9)
            axis.set_xlim(0.0, 1.0)
            axis.set_ylabel("R" if isinstance(model, LEGI) else "Active u (µM)")
            if row == 2:
                axis.set_xlabel("Relative position x/L")
            else:
                axis.tick_params(labelbottom=False)
    size_handles, size_labels = axes[1, 0].get_legend_handles_labels()
    figure.legend(
        size_handles,
        size_labels,
        loc="lower center",
        ncol=len(SIZE_LENGTHS),
        frameon=False,
        fontsize=9,
    )
    figure.suptitle(
        "Jilkine two-ended stimulus qualification\n"
        "Top: L=10 µm common-stimulus time profiles; "
        "rows 2–3: L=1, 3, 17, 25, 39 µm",
        fontsize=13,
    )
    figure.tight_layout(rect=(0.0, 0.05, 1.0, 0.955))
    figure.savefig(OUTPUT / "jilkine_qualification.png", dpi=180)
    plt.close(figure)
    return rows


def qualify_debelly() -> list[dict[str, object]]:
    protocol = GaussianPulse()
    variants = {
        "intact": DeBellyParameters(),
        "rac_to_rho_off": DeBellyParameters(rac_to_rho_coupling=0.0),
        "rho_to_rac_off": DeBellyParameters(rho_to_rac_coupling=0.0),
        "both_off": DeBellyParameters(
            rac_to_rho_coupling=0.0,
            rho_to_rac_coupling=0.0,
        ),
    }
    runs = {
        name: simulate(
            DeBelly(parameters),
            protocol=protocol,
            cells=100,
            end_time=350.0,
            output_interval=2.0,
            step=0.25,
        )
        for name, parameters in variants.items()
    }
    rows = []
    for name, trajectory in runs.items():
        relative = trajectory.x / 20.0
        distal = relative >= 0.8
        summary = phenotype_summary(trajectory, analysis_start=50.0)
        rows.append(
            {
                "model": trajectory.model,
                "protocol": "debelly_single_opto",
                "coupling_variant": name,
                "distal_rac_peak_change": float(
                    np.max(np.mean(trajectory.observable[:, distal], axis=1) - np.mean(trajectory.observable[0, distal]))
                ),
                "maximum_absolute_tension": float(np.max(np.abs(trajectory.fields["tension"]))),
                "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
                **summary,
            }
        )

    intact = runs["intact"]
    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, field, title in (
        (axes[0, 0], "rac", "Intact: ΔRac"),
        (axes[0, 1], "rho", "Intact: ΔRho"),
    ):
        change = intact.fields[field] - intact.fields[field][0]
        limit = max(float(np.max(np.abs(change))), 1e-12)
        image = axis.imshow(
            change,
            origin="lower",
            aspect="auto",
            cmap="RdBu_r",
            vmin=-limit,
            vmax=limit,
            extent=(0, 1, intact.time[0], intact.time[-1]),
        )
        axis.axhline(50, color="black", ls=":", lw=0.7)
        axis.axhline(150, color="black", ls=":", lw=0.7)
        axis.set_title(title)
        axis.set_xlabel("Relative position")
        axis.set_ylabel("Time (s)")
        figure.colorbar(image, ax=axis, shrink=0.75)

    colors = dict(
        intact="#222222",
        rac_to_rho_off="#D55E00",
        rho_to_rac_off="#0072B2",
        both_off="#7A5195",
    )
    for name, trajectory in runs.items():
        axes[1, 0].plot(
            trajectory.x / 20.0,
            trajectory.observable[-1] - trajectory.observable[0],
            label=name.replace("_", " "),
            color=colors[name],
        )
        distal = trajectory.x >= 16.0
        distal_rac = np.mean(trajectory.observable[:, distal], axis=1)
        axes[1, 1].plot(
            trajectory.time,
            distal_rac - distal_rac[0],
            label=name.replace("_", " "),
            color=colors[name],
        )
    axes[1, 0].set_xlabel("Relative position")
    axes[1, 0].set_ylabel("Final ΔRac")
    axes[1, 1].set_xlabel("Time (s)")
    axes[1, 1].set_ylabel("Distal ΔRac")
    axes[1, 1].axvspan(50, 150, color="#BBBBBB", alpha=0.25)
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 1].legend(frameon=False, fontsize=8)
    figure.suptitle("De Belly 1D mechanochemical qualification", fontsize=15)
    figure.tight_layout()
    figure.savefig(OUTPUT / "debelly_qualification.png", dpi=180)
    plt.close(figure)
    return rows


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "scope": "compact model qualification; not a parameter screen",
        "jilkine": qualify_jilkine(),
        "debelly": qualify_debelly(),
    }
    destination = OUTPUT / "qualification_summary.json"
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(destination)


if __name__ == "__main__":
    main()
