"""Local bistability and single-pulse regime maps for strip wave-pinning."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import os
from pathlib import Path

import numpy as np

from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from .model import (
    ActivationPulse,
    SingleSpeciesParameters,
    StripWavePinningChemistry,
    local_reaction_roots,
)
from .observables import classify_trajectory
from .reporting import Figure, write_figure_manifest


STATE_COLORS = {
    "homogeneous_low": "#9E9E9E",
    "homogeneous_high": "#F0E442",
    "single_polarized": "#0072B2",
    "dual_polarized": "#009E73",
    "multi_peak": "#56B4E9",
    "transient_coarsening": "#D55E00",
    "oscillatory_dynamic": "#CC79A7",
}


def output_path(default: str) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else Path.cwd() / "results" / default


def parse_values(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def common_parameters(args, *, mu: float, dimensionless_length: float, h: float):
    return SingleSpeciesParameters(
        total_abundance=mu,
        dimensionless_length=dimensionless_length,
        cytoplasm_depth=h,
        diffusion_ratio=args.rho,
        basal_ratio=args.kappa,
        feedback_ratio=args.gamma,
        hill_coefficient=args.hill,
        half_saturation=1.0,
        inactivation_rate=1.0,
        domain_length=args.length,
    )


def run_single_pulse(args, mu: float, dimensionless_length: float, h: float) -> dict[str, object]:
    parameters = common_parameters(args, mu=mu, dimensionless_length=dimensionless_length, h=h)
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
    profile = result.membrane[:, 0]
    classified = classify_trajectory(
        result.times, profile, half_saturation=parameters.half_saturation
    )
    membrane_total, cytoplasm_total = result.compartment_totals()
    conserved = membrane_total[:, 0] + h * cytoplasm_total[:, 0]
    conservation_error = float(
        np.max(np.abs(conserved - conserved[0])) / max(abs(conserved[0]), 1.0e-14)
    )
    endpoint_contrast = (
        profile[:, -1] - profile[:, 0]
    ) / np.maximum(np.abs(profile[:, -1]) + np.abs(profile[:, 0]), 1.0e-12)
    return {
        "mu": mu,
        "Lambda": dimensionless_length,
        "h": h,
        "terminal_polarity_status": classified["state"],
        "time_to_polarization": classified["time_to_polarization"],
        "terminal_polarization_amplitude": classified["terminal_polarization_amplitude"],
        "terminal_pole_count": classified["terminal_pole_count"],
        "terminal_high_width_fraction": classified["terminal_high_width_fraction"],
        "terminal_polarity_contrast": float(endpoint_contrast[-1]),
        "late_maximum_normalized_change": classified["late_maximum_normalized_change"],
        "maximum_conservation_error": conservation_error,
    }


def run_single_pulse_condition(task: tuple) -> dict[str, object]:
    args, mu, dimensionless_length, h = task
    return run_single_pulse(args, mu, dimensionless_length, h)


def local_bistability(args) -> None:
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    parameters = common_parameters(
        args, mu=args.baseline_mu, dimensionless_length=args.baseline_lambda,
        h=args.baseline_h,
    )
    inactive_values = np.linspace(args.inactive_min, args.inactive_max, args.inactive_count)
    records = []
    for inactive in inactive_values:
        for root in local_reaction_roots(inactive, parameters):
            records.append({
                "inactive": float(inactive),
                "active": root["active"],
                "stable": root["stable"],
            })
    triple = [
        float(value) for value in inactive_values
        if len(local_reaction_roots(float(value), parameters)) == 3
    ]
    summary = {
        "analysis": "local_bistability",
        "parameters": vars(parameters),
        "three_root_inactive_interval": (
            [min(triple), max(triple)] if triple else None
        ),
        "has_local_bistability": bool(triple),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "roots.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("inactive", "active", "stable"))
        writer.writeheader()
        writer.writerows(records)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    for stable, color, label in ((True, "#0072B2", "stable"), (False, "#D55E00", "unstable")):
        subset = [row for row in records if row["stable"] is stable]
        ax.scatter([row["inactive"] for row in subset], [row["active"] for row in subset],
                   s=8, color=color, label=label)
    ax.set(xlabel="inactive cytoplasmic concentration v",
           ylabel="active membrane concentration u",
           title="Local roots of J(u,v)=0")
    ax.legend(frameon=False)
    fig.savefig(output / "local_bistability.png", dpi=200)
    plt.close(fig)

    counts = [sum(row["inactive"] == float(value) for row in records) for value in inactive_values]
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    ax.step(inactive_values, counts, where="mid", color="#0072B2")
    ax.fill_between(inactive_values, 0, counts, where=np.asarray(counts) == 3,
                    step="mid", color="#009E73", alpha=0.25,
                    label="three-root interval")
    ax.set(xlabel="inactive cytoplasmic concentration v", ylabel="number of roots",
           yticks=(1, 2, 3), title="Local bistability window")
    ax.legend(frameon=False)
    fig.savefig(output / "local_root_count.png", dpi=200)
    plt.close(fig)
    write_figure_manifest(
        output,
        title="Local reaction bistability",
        introduction=(
            "Roots of J(u,v)=0 are evaluated before PDE simulation. A range with "
            "two stable roots separated by one unstable root is the necessary "
            "local-kinetic filter for classical wave-pinning."
        ),
        figures=[
            Figure("local_bistability.png", "Stable and unstable reaction roots",
                   "Root branches as cytoplasmic inactive concentration is varied."),
            Figure("local_root_count.png", "Three-root bistability window",
                   "The shaded interval identifies cytoplasmic levels with three roots."),
        ],
        data_files=("summary.json", "roots.csv", "figure_manifest.json"),
    )


def regime_map(args) -> None:
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    if args.render_existing:
        with (output / "summary.json").open() as handle:
            existing_summary = json.load(handle)
        with (output / "regime_map.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            row["mu"] = float(row["mu"])
            row["Lambda"] = float(row["Lambda"])
            row["h"] = float(row["h"])
        mus = np.asarray(existing_summary["axes"]["mu"], dtype=float)
        lambdas = np.asarray(existing_summary["axes"]["Lambda"], dtype=float)
        h_values = [float(value) for value in existing_summary["axes"]["h"]]
        figures = render_regime_maps(output, rows, mus, lambdas, h_values)
        write_figure_manifest(
            output,
            title="Single-pulse wave-pinning regime map",
            introduction=(
                "A matched transient end stimulus is applied across total "
                "abundance mu and dimensionless length Lambda. Columns compare "
                "cytoplasmic-depth h values while kinetic controls remain fixed."
            ),
            figures=figures,
            data_files=("summary.json", "regime_map.csv", "figure_manifest.json"),
        )
        return
    mus = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    h_values = parse_values(args.h_values)
    tasks = [
        (args, float(mu), float(length), float(h))
        for h in h_values for mu in mus for length in lambdas
    ]
    if args.workers <= 1:
        rows = [run_single_pulse_condition(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            rows = list(executor.map(run_single_pulse_condition, tasks, chunksize=1))
    fields = tuple(rows[0])
    with (output / "regime_map.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    counts = {}
    for row in rows:
        status = row["terminal_polarity_status"]
        counts[status] = counts.get(status, 0) + 1
    summary = {
        "analysis": "single_pulse_regime_map",
        "axes": {"mu": mus.tolist(), "Lambda": lambdas.tolist(), "h": h_values},
        "fixed": {
            "rho": args.rho, "kappa": args.kappa, "Gamma": args.gamma,
            "n": args.hill, "duration_seconds": args.duration,
        },
        "numerics": {
            "nodes": args.nodes, "output_interval": args.output_interval,
            "max_step": args.max_step, "workers": args.workers,
        },
        "stimulus": {
            "amplitude": args.stimulus_amplitude, "center_fraction": args.stimulus_center,
            "width_fraction": args.stimulus_width,
            "start": args.stimulus_start, "stop": args.stimulus_stop,
        },
        "state_counts": counts,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = render_regime_maps(output, rows, mus, lambdas, h_values)
    write_figure_manifest(
        output,
        title="Single-pulse wave-pinning regime map",
        introduction=(
            "A matched transient end stimulus is applied across the (mu, Lambda) "
            "grid. Columns compare geometry slices h while all kinetic controls "
            "remain fixed. Classification uses the late trajectory."
        ),
        figures=figures,
        data_files=("summary.json", "regime_map.csv", "figure_manifest.json"),
    )


def plot_extent(values: np.ndarray) -> tuple[float, float]:
    if len(values) > 1:
        return float(values[0]), float(values[-1])
    delta = max(0.5, 0.05 * abs(float(values[0])))
    return float(values[0] - delta), float(values[0] + delta)


def render_regime_maps(output, rows, mus, lambdas, h_values) -> list[Figure]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    states = list(STATE_COLORS)
    cmap = ListedColormap([STATE_COLORS[state] for state in states])
    fig, axes = plt.subplots(1, len(h_values), figsize=(5.2 * len(h_values), 4.8),
                             squeeze=False, constrained_layout=True)
    for ax, h in zip(axes[0], h_values):
        subset = [row for row in rows if row["h"] == h]
        grid = np.asarray([
            states.index(row["terminal_polarity_status"]) for row in subset
        ]).reshape(
            len(mus), len(lambdas)
        )
        ax.imshow(grid, origin="lower", aspect="auto", cmap=cmap,
                  vmin=-0.5, vmax=len(states)-0.5,
                  extent=(*plot_extent(lambdas), *plot_extent(mus)))
        ax.set(
            title=f"cytoplasmic depth h = {h:g}",
            xlabel="dimensionless length Λ",
            ylabel="total abundance μ",
        )
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=color,
                          label=state.replace("_", " "))
               for state, color in STATE_COLORS.items()]
    fig.legend(handles=handles, loc="outside lower center",
               ncol=min(4, len(h_values) + 1), frameon=False)
    fig.suptitle("Single-pulse wave-pinning outcome")
    fig.savefig(output / "single_pulse_phase_diagram.png", dpi=200)
    plt.close(fig)
    figures = [Figure(
        "single_pulse_phase_diagram.png", "Terminal-state phase diagram",
        "Categorical late-time outcome across abundance, domain size, and geometry.",
    )]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    bottoms = np.zeros(len(h_values))
    for state in states:
        counts = np.asarray([
            sum(row["h"] == h and row["terminal_polarity_status"] == state for row in rows)
            for h in h_values
        ])
        if np.any(counts):
            ax.bar(
                [f"h = {h:g}" for h in h_values], counts, bottom=bottoms,
                color=STATE_COLORS[state], label=state.replace("_", " "),
            )
            bottoms += counts
    ax.set(
        xlabel="cytoplasmic depth",
        ylabel="number of sampled (μ, Λ) conditions",
        title="Single-pulse outcome counts by geometry",
    )
    ax.legend(frameon=False, ncol=2)
    fig.savefig(output / "single_pulse_outcome_counts.png", dpi=200)
    plt.close(fig)
    figures.append(Figure(
        "single_pulse_outcome_counts.png",
        "Outcome counts by cytoplasmic depth",
        "A direct count of the 289 sampled (μ, Λ) conditions in each h panel. Counts summarize this grid and are not biological probabilities.",
    ))
    return figures


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--analysis", choices=("local-bistability", "regime-map"), required=True)
    result.add_argument("--output", type=Path)
    result.add_argument(
        "--render-existing", action="store_true",
        help="regenerate regime-map figures from its existing CSV without simulations",
    )
    result.add_argument("--length", type=float, default=10.0)
    result.add_argument("--rho", type=float, default=100.0)
    result.add_argument("--kappa", type=float, default=0.067)
    result.add_argument("--gamma", type=float, default=1.0)
    result.add_argument("--hill", type=float, default=2.0)
    result.add_argument("--baseline-mu", type=float, default=2.2683312)
    result.add_argument("--baseline-lambda", type=float, default=31.6227766017)
    result.add_argument("--baseline-h", type=float, default=1.0)
    result.add_argument("--inactive-min", type=float, default=0.0)
    result.add_argument("--inactive-max", type=float, default=3.0)
    result.add_argument("--inactive-count", type=int, default=301)
    result.add_argument("--mu-min", type=float, default=0.5)
    result.add_argument("--mu-max", type=float, default=5.0)
    result.add_argument("--mu-count", type=int, default=17)
    result.add_argument("--lambda-min", type=float, default=10.0)
    result.add_argument("--lambda-max", type=float, default=100.0)
    result.add_argument("--lambda-count", type=int, default=17)
    result.add_argument("--h-values", default="1,2,4")
    result.add_argument("--duration", type=float, default=300.0)
    result.add_argument("--nodes", type=int, default=81)
    result.add_argument("--output-interval", type=float, default=2.0)
    result.add_argument("--max-step", type=float, default=0.5)
    result.add_argument(
        "--workers", type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
        help="parallel simulation processes; defaults to the allocated SLURM CPUs",
    )
    result.add_argument("--stimulus-amplitude", type=float, default=1.0)
    result.add_argument("--stimulus-center", type=float, default=0.1)
    result.add_argument("--stimulus-width", type=float, default=0.05)
    result.add_argument("--stimulus-start", type=float, default=10.0)
    result.add_argument("--stimulus-stop", type=float, default=30.0)
    return result


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    if args.output is None:
        default = "single_species/calibration/local_bistability" if args.analysis == "local-bistability" else "single_species/calibration/regime_map"
        args.output = output_path(default)
    if args.analysis == "local-bistability":
        local_bistability(args)
    else:
        regime_map(args)


if __name__ == "__main__":
    main()
