"""Exploratory figures for selected single-species trajectories.

These outputs are deliberately not given a figure manifest, so the aggregate
report will not include them automatically.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from migration1d.model import (
    MechanicsParameters,
    MigrationModel,
    SimulationConfig,
    run_simulation,
)
from .model import ActivationPulse, SingleSpeciesParameters, StripWavePinningChemistry
from .observables import polarization_amplitude, polarized_domain_count


def domain_intervals(
    position: np.ndarray, profile: np.ndarray
) -> tuple[float | None, list[tuple[float, float]]]:
    """Return the plotting threshold and qualifying connected intervals."""
    if polarized_domain_count(profile) == 0:
        return None, []
    threshold = float(np.min(profile) + 0.5 * (np.max(profile) - np.min(profile)))
    mask = profile > threshold
    starts = np.flatnonzero(mask & np.r_[True, ~mask[:-1]])
    ends = np.flatnonzero(mask & np.r_[~mask[1:], True])
    return threshold, [
        (float(position[start]), float(position[end]))
        for start, end in zip(starts, ends)
    ]


def simulate_example(args):
    parameters = SingleSpeciesParameters(
        total_abundance=args.mu,
        dimensionless_length=args.dimensionless_length,
        cytoplasm_depth=args.h,
        diffusion_ratio=args.rho,
        basal_ratio=args.kappa,
        feedback_ratio=args.gamma,
        hill_coefficient=args.hill,
        domain_length=args.length,
    )
    chemistry = StripWavePinningChemistry(
        parameters,
        activation=ActivationPulse(
            amplitude=args.stimulus_amplitude,
            center_fraction=args.stimulus_center,
            width_fraction=args.stimulus_width,
            start=args.stimulus_start,
            stop=args.stimulus_stop,
        ),
    )
    result = run_simulation(
        MigrationModel(
            mechanics=MechanicsParameters(reference_length=args.length),
            chemistry=chemistry,
        ),
        SimulationConfig(
            nodes=args.nodes,
            duration=args.duration,
            output_interval=args.output_interval,
            max_step=args.max_step,
            noise_amplitude=0.0,
            mechanics_activation_time=None,
        ),
    )
    return parameters, result


def render_dual_single_pulse_example(args) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    parameters, result = simulate_example(args)
    u = result.membrane[:, 0]
    v = result.cytoplasm[:, 0]
    amplitudes = np.asarray([polarization_amplitude(profile) for profile in u])
    counts = np.asarray([polarized_domain_count(profile) for profile in u])
    dual_indices = np.flatnonzero(counts == 2)
    first_dual = float(result.times[dual_indices[0]]) if len(dual_indices) else None
    membrane_total, cytoplasm_total = result.compartment_totals()
    conserved = membrane_total[:, 0] + parameters.cytoplasm_depth * cytoplasm_total[:, 0]
    conservation_error = float(
        np.max(np.abs(conserved - conserved[0])) / max(abs(conserved[0]), 1.0e-14)
    )

    fig, ax = plt.subplots(figsize=(9.2, 5.4), constrained_layout=True)
    image = ax.pcolormesh(
        result.position,
        np.broadcast_to(result.times[:, None], result.position.shape),
        u,
        shading="auto",
        cmap="magma",
    )
    ax.axhspan(
        args.stimulus_start, args.stimulus_stop,
        color="#54d6ff", alpha=0.18, label="left activation cue",
    )
    if first_dual is not None:
        ax.axhline(first_dual, color="#55efc4", linestyle="--", linewidth=1.8,
                   label=f"first two-domain profile: {first_dual:g} s")
    ax.set(
        xlabel="position along strip x",
        ylabel="time (s)",
        title="A second active domain emerges at the unstimulated right boundary",
    )
    fig.colorbar(image, ax=ax, label="active membrane GTPase u")
    ax.legend(loc="upper center", frameon=True)
    fig.savefig(output / "active_gtpase_kymograph.png", dpi=220)
    plt.close(fig)

    selected_times = [0.0, 20.0, 30.0, 40.0, 60.0, 100.0, 200.0, 300.0]
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.2), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, selected in zip(axes.ravel(), selected_times):
        index = int(np.argmin(np.abs(result.times - selected)))
        position = result.position[index]
        profile = u[index]
        threshold, intervals = domain_intervals(position, profile)
        ax.plot(position, profile, color="#175a8a", linewidth=2.2)
        if threshold is not None:
            ax.axhline(threshold, color="#D55E00", linestyle="--", linewidth=1.2,
                       label="within-profile threshold")
            for start, stop in intervals:
                ax.axvspan(start, stop, color="#56B4E9", alpha=0.22)
        ax.set_title(
            f"t = {result.times[index]:g} s   P = {amplitudes[index]:.2f}   "
            f"domains = {counts[index]}"
        )
        ax.set_xlabel("position x")
        ax.set_ylabel("active GTPase u")
    fig.suptitle(
        "Single left cue: one stimulated domain followed by a second right domain",
        fontsize=16,
    )
    fig.savefig(output / "thresholded_spatial_profiles.png", dpi=220)
    plt.close(fig)

    midpoint = u.shape[1] // 2
    background = np.min(u, axis=1)
    left_peak = np.max(u[:, :midpoint], axis=1) - background
    right_peak = np.max(u[:, midpoint:], axis=1) - background
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 8.2), sharex=True,
                             constrained_layout=True)
    axes[0].plot(result.times, amplitudes, color="#175a8a", linewidth=2)
    axes[0].axhline(0.5, color="#D55E00", linestyle="--", label="P=0.5 cutoff")
    axes[0].set(ylabel="polarization amplitude P", title="Trajectory diagnostics")
    axes[0].legend(frameon=False)
    axes[1].step(result.times, counts, where="post", color="#009E73", linewidth=2)
    axes[1].set(ylabel="polarized domains", yticks=(0, 1, 2, 3), ylim=(-0.15, 3.15))
    axes[2].plot(result.times, left_peak, color="#0072B2", linewidth=2,
                 label="strongest left-half activity above background")
    axes[2].plot(result.times, right_peak, color="#CC79A7", linewidth=2,
                 label="strongest right-half activity above background")
    axes[2].set(xlabel="time (s)", ylabel="peak excess")
    axes[2].legend(frameon=False)
    for ax in axes:
        ax.axvspan(args.stimulus_start, args.stimulus_stop, color="#54d6ff", alpha=0.15)
        if first_dual is not None:
            ax.axvline(first_dual, color="#009E73", linestyle=":", linewidth=1.2)
    fig.savefig(output / "domain_dynamics.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    active_image = axes[0].pcolormesh(
        result.position,
        np.broadcast_to(result.times[:, None], result.position.shape),
        u,
        shading="auto",
        cmap="magma",
    )
    inactive_image = axes[1].pcolormesh(
        result.position,
        np.broadcast_to(result.times[:, None], result.position.shape),
        v,
        shading="auto",
        cmap="viridis",
    )
    axes[0].set(xlabel="position x", ylabel="time (s)", title="Active membrane pool u")
    axes[1].set(xlabel="position x", ylabel="time (s)", title="Inactive cytoplasmic pool v")
    fig.colorbar(active_image, ax=axes[0], label="u")
    fig.colorbar(inactive_image, ax=axes[1], label="v")
    fig.suptitle("Redistribution between the conserved active and inactive pools")
    fig.savefig(output / "active_inactive_kymographs.png", dpi=220)
    plt.close(fig)

    final_threshold, final_intervals = domain_intervals(result.position[-1], u[-1])
    summary = {
        "analysis": "exploratory_single_pulse_dual_example",
        "included_in_aggregate_report": False,
        "parameters": vars(parameters),
        "stimulus": {
            "amplitude": args.stimulus_amplitude,
            "center_fraction": args.stimulus_center,
            "width_fraction": args.stimulus_width,
            "start": args.stimulus_start,
            "stop": args.stimulus_stop,
        },
        "first_two_domain_time": first_dual,
        "terminal_polarization_amplitude": float(amplitudes[-1]),
        "terminal_polarized_domain_count": int(counts[-1]),
        "terminal_within_profile_threshold": final_threshold,
        "terminal_domain_intervals": final_intervals,
        "maximum_conservation_error": conservation_error,
        "figures": [
            "active_gtpase_kymograph.png",
            "thresholded_spatial_profiles.png",
            "domain_dynamics.png",
            "active_inactive_kymographs.png",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument(
        "--output", type=Path,
        default=Path("migration1d/results/single_species/exploratory/single_pulse_dual_example"),
    )
    result.add_argument("--mu", type=float, default=2.46875)
    result.add_argument("--dimensionless-length", type=float, default=83.125)
    result.add_argument("--h", type=float, default=1.0)
    result.add_argument("--length", type=float, default=10.0)
    result.add_argument("--rho", type=float, default=100.0)
    result.add_argument("--kappa", type=float, default=0.067)
    result.add_argument("--gamma", type=float, default=1.0)
    result.add_argument("--hill", type=float, default=2.0)
    result.add_argument("--stimulus-amplitude", type=float, default=1.0)
    result.add_argument("--stimulus-center", type=float, default=0.1)
    result.add_argument("--stimulus-width", type=float, default=0.05)
    result.add_argument("--stimulus-start", type=float, default=10.0)
    result.add_argument("--stimulus-stop", type=float, default=30.0)
    result.add_argument("--duration", type=float, default=300.0)
    result.add_argument("--nodes", type=int, default=81)
    result.add_argument("--output-interval", type=float, default=2.0)
    result.add_argument("--max-step", type=float, default=0.5)
    return result


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    render_dual_single_pulse_example(args)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
