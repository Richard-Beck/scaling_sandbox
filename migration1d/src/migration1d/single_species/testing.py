"""Directly interpretable dual-domain and polarity-reversal tests."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import os
from pathlib import Path

import numpy as np

from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from .calibration import plot_extent
from .model import (
    ActivationPulse,
    DualPulseInitialCondition,
    PulseTrain,
    SingleSpeciesParameters,
    StripWavePinningChemistry,
)
from .observables import polarized_domain_count
from .reporting import Figure, write_figure_manifest


DOMAIN_COLORS = ("#d9d9d9", "#0072B2", "#009E73", "#CC79A7")
REVERSAL_COLORS = {
    "not_initially_left_polarized": "#9E9E9E",
    "reversed": "#009E73",
    "retained_original_pole": "#0072B2",
    "lost_polarity": "#F0E442",
    "multiple_polarized_domains": "#CC79A7",
    "ambiguous": "#D55E00",
}


def output_path(default: str) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else Path.cwd() / "results" / default


def values(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def parameters(args, mu: float, dimensionless_length: float, h: float):
    return SingleSpeciesParameters(
        total_abundance=mu,
        dimensionless_length=dimensionless_length,
        cytoplasm_depth=h,
        diffusion_ratio=args.rho,
        basal_ratio=args.kappa,
        feedback_ratio=args.gamma,
        hill_coefficient=args.hill,
        domain_length=args.length,
    )


def simulate(args, chemistry: StripWavePinningChemistry):
    return run_simulation(
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


def nearest_index(times: np.ndarray, selected: float) -> int:
    return int(np.argmin(np.abs(times - selected)))


def count_at_times(result, selected_times: list[float]) -> dict[str, int]:
    profiles = result.membrane[:, 0]
    return {
        time_key(selected): polarized_domain_count(profiles[nearest_index(result.times, selected)])
        for selected in selected_times
    }


def time_key(selected: float) -> str:
    return f"polarized_domain_count_t{selected:g}"


def endpoint_contrast(profiles: np.ndarray) -> np.ndarray:
    return (profiles[:, -1] - profiles[:, 0]) / np.maximum(
        np.abs(profiles[:, -1]) + np.abs(profiles[:, 0]), 1.0e-12
    )


def conservation_error(result, h: float) -> float:
    membrane_total, cytoplasm_total = result.compartment_totals()
    conserved = membrane_total[:, 0] + h * cytoplasm_total[:, 0]
    return float(
        np.max(np.abs(conserved - conserved[0]))
        / max(abs(conserved[0]), 1.0e-14)
    )


def first_sustained_time(
    times: np.ndarray,
    condition: np.ndarray,
    *,
    start: float,
    duration: float,
) -> float | None:
    for index in np.flatnonzero((times >= start) & condition):
        end = int(np.searchsorted(times, times[index] + duration, side="left"))
        if end < len(times) and np.all(condition[index : end + 1]):
            return float(times[index] - start)
    return None


def parallel_map(function, tasks: list[tuple], workers: int) -> list[dict[str, object]]:
    """Map independent parameter conditions onto the allocated CPU workers."""
    if workers <= 1:
        return [function(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(function, tasks, chunksize=1))


def run_dual_condition(task: tuple) -> dict[str, object]:
    args, mu, dimensionless_length, h, asymmetry, selected_times = task
    p = parameters(args, mu, dimensionless_length, h)
    chemistry = StripWavePinningChemistry(
        p,
        initial_condition=DualPulseInitialCondition(
            amplitude=args.initial_amplitude,
            width_fraction=args.initial_width,
            asymmetry=asymmetry,
        ),
    )
    result = simulate(args, chemistry)
    counts = count_at_times(result, selected_times)
    return {
        "mu": mu,
        "Lambda": dimensionless_length,
        "h": h,
        "epsilon": asymmetry,
        **counts,
        "maximum_conservation_error": conservation_error(result, p.cytoplasm_depth),
    }


def dual_pole_map(args) -> None:
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    mus = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    h_values = values(args.h_values)
    asymmetries = values(args.asymmetries)
    selected_times = values(args.domain_count_times)
    if any(time < 0.0 or time > args.duration for time in selected_times):
        raise ValueError("domain-count times must fall within the simulation duration")

    tasks = [
        (
            args, float(mu), float(dimensionless_length), float(h),
            float(asymmetry), selected_times,
        )
        for h in h_values
        for asymmetry in asymmetries
        for mu in mus
        for dimensionless_length in lambdas
    ]
    rows = parallel_map(run_dual_condition, tasks, args.workers)

    with (output / "dual_pole_map.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    transition_counts: dict[str, dict[str, int]] = {}
    for asymmetry in asymmetries:
        transitions: dict[str, int] = {}
        for row in rows:
            if row["epsilon"] != asymmetry:
                continue
            transition = "→".join(str(row[time_key(time)]) for time in selected_times)
            transitions[transition] = transitions.get(transition, 0) + 1
        transition_counts[f"epsilon={asymmetry:g}"] = transitions
    summary = {
        "analysis": "dual_pole_domain_counts",
        "question": "How many distinct polarized domains remain as two seeded sites evolve?",
        "axes": {
            "mu": mus.tolist(),
            "Lambda": lambdas.tolist(),
            "h": h_values,
            "epsilon": asymmetries,
        },
        "domain_count_times_seconds": selected_times,
        "fixed": {
            "rho": args.rho,
            "kappa": args.kappa,
            "Gamma": args.gamma,
            "n": args.hill,
            "duration_seconds": args.duration,
        },
        "numerics": {
            "nodes": args.nodes,
            "output_interval": args.output_interval,
            "max_step": args.max_step,
            "workers": args.workers,
        },
        "initial_condition": {
            "centers": [0.2, 0.8],
            "amplitude": args.initial_amplitude,
            "width_fraction": args.initial_width,
        },
        "domain_definition": {
            "minimum_polarization_amplitude": 0.5,
            "within_profile_threshold": "min(u) + 0.5 * (max(u) - min(u))",
        },
        "transition_counts": transition_counts,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = render_domain_count_maps(
        output, rows, mus, lambdas, h_values, asymmetries, selected_times
    )
    write_figure_manifest(
        output,
        title="Two seeded domains: persistence and loss",
        introduction=(
            "The directly observed output is the number of spatially distinct "
            "polarized domains at 0, 100, and 300 seconds. Each figure holds the "
            "initial left-right imbalance epsilon fixed."
        ),
        figures=figures,
        data_files=("summary.json", "dual_pole_map.csv", "figure_manifest.json"),
    )


def render_domain_count_maps(
    output, rows, mus, lambdas, h_values, asymmetries, selected_times
) -> list[Figure]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    cmap = ListedColormap(DOMAIN_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(DOMAIN_COLORS) + 0.5), cmap.N)
    figures = []
    for asymmetry in asymmetries:
        fig, axes = plt.subplots(
            len(selected_times), len(h_values),
            figsize=(5.0 * len(h_values), 3.8 * len(selected_times)),
            squeeze=False,
            constrained_layout=True,
        )
        image = None
        for row_index, selected in enumerate(selected_times):
            for column_index, h in enumerate(h_values):
                ax = axes[row_index, column_index]
                subset = [
                    row for row in rows
                    if row["h"] == h and row["epsilon"] == asymmetry
                ]
                grid = np.asarray([
                    int(row[time_key(selected)]) for row in subset
                ]).reshape(len(mus), len(lambdas))
                image = ax.imshow(
                    grid,
                    origin="lower",
                    aspect="auto",
                    cmap=cmap,
                    norm=norm,
                    extent=(*plot_extent(lambdas), *plot_extent(mus)),
                )
                if row_index == 0:
                    ax.set_title(f"cytoplasmic depth h = {h:g}")
                ax.set_xlabel("dimensionless length Λ")
                ax.set_ylabel(f"t = {selected:g} s\ntotal abundance μ")
        colorbar = fig.colorbar(
            image,
            ax=axes.ravel().tolist(),
            ticks=range(len(DOMAIN_COLORS)),
            label="number of polarized domains",
        )
        colorbar.ax.set_yticklabels(["0", "1", "2", "3+"])
        fig.suptitle(f"Polarized-domain count through time; ε = {asymmetry:g}")
        filename = f"polarized_domain_counts_epsilon_{asymmetry:g}.png"
        fig.savefig(output / filename, dpi=200)
        plt.close(fig)
        figures.append(Figure(
            filename,
            f"Polarized domains at 0, 100, and 300 s; ε={asymmetry:g}",
            "Rows are observation times, columns are cytoplasmic-depth values, and each pixel is one (μ, Λ) condition. Gray means no sufficiently contrasted domain, blue one, and green two.",
        ))
    return figures


def run_reversal_condition(task: tuple) -> dict[str, object]:
    args, mu, dimensionless_length, h = task
    p = parameters(args, mu, dimensionless_length, h)
    activation = PulseTrain((
        ActivationPulse(
            amplitude=args.stimulus_amplitude,
            center_fraction=0.1,
            width_fraction=args.stimulus_width,
            start=10.0,
            stop=30.0,
        ),
        ActivationPulse(
            amplitude=args.reversal_amplitude,
            center_fraction=0.9,
            width_fraction=args.stimulus_width,
            start=args.reversal_start,
            stop=args.reversal_start + args.reversal_duration,
        ),
    ))
    result = simulate(args, StripWavePinningChemistry(p, activation=activation))
    profiles = result.membrane[:, 0]
    contrast = endpoint_contrast(profiles)
    before_indices = np.flatnonzero(result.times < args.reversal_start)
    before = int(before_indices[-1])
    before_count = polarized_domain_count(profiles[before])
    final_count = polarized_domain_count(profiles[-1])
    initially_left = before_count == 1 and contrast[before] <= -args.reversal_threshold
    if not initially_left:
        outcome = "not_initially_left_polarized"
    elif final_count == 1 and contrast[-1] >= args.reversal_threshold:
        outcome = "reversed"
    elif final_count == 1 and contrast[-1] <= -args.reversal_threshold:
        outcome = "retained_original_pole"
    elif final_count == 0:
        outcome = "lost_polarity"
    elif final_count >= 2:
        outcome = "multiple_polarized_domains"
    else:
        outcome = "ambiguous"
    counts = np.asarray([polarized_domain_count(profile) for profile in profiles])
    reversal_time = first_sustained_time(
        result.times,
        (counts == 1) & (contrast >= args.reversal_threshold),
        start=args.reversal_start,
        duration=10.0,
    ) if initially_left else None
    return {
        "mu": mu,
        "Lambda": dimensionless_length,
        "h": h,
        "prechallenge_time": float(result.times[before]),
        "prechallenge_polarized_domain_count": before_count,
        "prechallenge_endpoint_contrast": float(contrast[before]),
        "terminal_polarized_domain_count": final_count,
        "terminal_endpoint_contrast": float(contrast[-1]),
        "reversal_outcome": outcome,
        "reversal_time_after_challenge": reversal_time,
        "maximum_conservation_error": conservation_error(result, p.cytoplasm_depth),
    }


def reversal_map(args) -> None:
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    mus = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    h_values = values(args.h_values)
    tasks = [
        (args, float(mu), float(dimensionless_length), float(h))
        for h in h_values
        for mu in mus
        for dimensionless_length in lambdas
    ]
    rows = parallel_map(run_reversal_condition, tasks, args.workers)

    with (output / "polarity_reversal_map.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    outcome_counts: dict[str, int] = {}
    for row in rows:
        outcome = str(row["reversal_outcome"])
        outcome_counts[outcome] = outcome_counts.get(outcome, 0) + 1
    summary = {
        "analysis": "polarity_reversal_map",
        "question": "Which initially left-polarized states reverse after a matched right-hand cue?",
        "axes": {"mu": mus.tolist(), "Lambda": lambdas.tolist(), "h": h_values},
        "fixed": {
            "rho": args.rho,
            "kappa": args.kappa,
            "Gamma": args.gamma,
            "n": args.hill,
            "duration_seconds": args.duration,
        },
        "numerics": {
            "nodes": args.nodes,
            "output_interval": args.output_interval,
            "max_step": args.max_step,
            "workers": args.workers,
        },
        "initial_left_pulse": {
            "center_fraction": 0.1,
            "start": 10.0,
            "stop": 30.0,
            "amplitude": args.stimulus_amplitude,
        },
        "right_reversal_pulse": {
            "center_fraction": 0.9,
            "start": args.reversal_start,
            "stop": args.reversal_start + args.reversal_duration,
            "amplitude": args.reversal_amplitude,
        },
        "direction_threshold": args.reversal_threshold,
        "outcome_counts": outcome_counts,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = render_reversal_maps(output, rows, mus, lambdas, h_values)
    write_figure_manifest(
        output,
        title="Polarity-reversal regime map",
        introduction=(
            "A left cue is followed by a matched right cue. Only conditions with "
            "one clearly left-biased domain immediately before the second cue "
            "are eligible to be called reversed."
        ),
        figures=figures,
        data_files=("summary.json", "polarity_reversal_map.csv", "figure_manifest.json"),
    )


def render_reversal_maps(output, rows, mus, lambdas, h_values) -> list[Figure]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    outcomes = list(REVERSAL_COLORS)
    cmap = ListedColormap([REVERSAL_COLORS[outcome] for outcome in outcomes])
    norm = BoundaryNorm(np.arange(-0.5, len(outcomes) + 0.5), cmap.N)
    fig, axes = plt.subplots(
        1, len(h_values), figsize=(5.2 * len(h_values), 4.8),
        squeeze=False, constrained_layout=True,
    )
    for ax, h in zip(axes[0], h_values):
        subset = [row for row in rows if row["h"] == h]
        grid = np.asarray([
            outcomes.index(str(row["reversal_outcome"])) for row in subset
        ]).reshape(len(mus), len(lambdas))
        ax.imshow(
            grid, origin="lower", aspect="auto", cmap=cmap, norm=norm,
            extent=(*plot_extent(lambdas), *plot_extent(mus)),
        )
        ax.set(
            title=f"cytoplasmic depth h = {h:g}",
            xlabel="dimensionless length Λ",
            ylabel="total abundance μ",
        )
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=REVERSAL_COLORS[item],
                   label=item.replace("_", " "))
        for item in outcomes
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Response to the opposing right-hand cue")
    fig.savefig(output / "polarity_reversal_outcome.png", dpi=200)
    plt.close(fig)
    figures = [Figure(
        "polarity_reversal_outcome.png",
        "Reversal outcome across the parameter grid",
        "Green conditions formed a left pole before the challenge and a right pole afterward. Gray conditions never formed the required initial left pole and therefore do not test reversal.",
    )]

    domain_cmap = ListedColormap(DOMAIN_COLORS)
    domain_norm = BoundaryNorm(np.arange(-0.5, len(DOMAIN_COLORS) + 0.5), domain_cmap.N)
    fig, axes = plt.subplots(
        2, len(h_values), figsize=(5.0 * len(h_values), 7.6),
        squeeze=False, constrained_layout=True,
    )
    image = None
    for row_index, (metric, label) in enumerate((
        ("prechallenge_polarized_domain_count", "before right cue (t=148 s)"),
        ("terminal_polarized_domain_count", "endpoint (t=300 s)"),
    )):
        for column_index, h in enumerate(h_values):
            ax = axes[row_index, column_index]
            subset = [row for row in rows if row["h"] == h]
            grid = np.asarray([int(row[metric]) for row in subset]).reshape(
                len(mus), len(lambdas)
            )
            image = ax.imshow(
                grid, origin="lower", aspect="auto", cmap=domain_cmap,
                norm=domain_norm,
                extent=(*plot_extent(lambdas), *plot_extent(mus)),
            )
            if row_index == 0:
                ax.set_title(f"cytoplasmic depth h = {h:g}")
            ax.set(xlabel="dimensionless length Λ", ylabel=f"{label}\ntotal abundance μ")
    colorbar = fig.colorbar(
        image, ax=axes.ravel().tolist(), ticks=range(len(DOMAIN_COLORS)),
        label="number of polarized domains",
    )
    colorbar.ax.set_yticklabels(["0", "1", "2", "3+"])
    fig.suptitle("Polarized domains before and after the opposing cue")
    fig.savefig(output / "polarity_reversal_domain_counts.png", dpi=200)
    plt.close(fig)
    figures.append(Figure(
        "polarity_reversal_domain_counts.png",
        "Domain count before the challenge and at 300 s",
        "This separates failure to establish initial polarity from loss, retention, or duplication after the right-hand cue.",
    ))

    fig, axes = plt.subplots(
        1, len(h_values), figsize=(5.2 * len(h_values), 4.5),
        squeeze=False, constrained_layout=True,
    )
    arrays = []
    for h in h_values:
        subset = [row for row in rows if row["h"] == h]
        arrays.append(np.asarray([
            np.nan if row["reversal_outcome"] != "reversed"
            or row["reversal_time_after_challenge"] is None
            else float(row["reversal_time_after_challenge"])
            for row in subset
        ]).reshape(len(mus), len(lambdas)))
    finite = np.concatenate([array[np.isfinite(array)] for array in arrays])
    if len(finite):
        vmin, vmax = float(np.min(finite)), float(np.max(finite))
        if np.isclose(vmin, vmax):
            vmax = vmin + 1.0
        image = None
        for ax, h, array in zip(axes[0], h_values, arrays):
            image = ax.imshow(
                array, origin="lower", aspect="auto", cmap="viridis",
                vmin=vmin, vmax=vmax,
                extent=(*plot_extent(lambdas), *plot_extent(mus)),
            )
            ax.set(
                title=f"cytoplasmic depth h = {h:g}",
                xlabel="dimensionless length Λ",
                ylabel="total abundance μ",
            )
        fig.colorbar(image, ax=axes.ravel().tolist(), label="seconds after right cue")
        fig.suptitle("Reversal time among successful reversals")
        fig.savefig(output / "polarity_reversal_time.png", dpi=200)
        plt.close(fig)
        figures.append(Figure(
            "polarity_reversal_time.png",
            "Time required for successful reversal",
            "Only green 'reversed' conditions from the outcome map are colored; blank cells were ineligible or did not end with a right-hand pole.",
        ))
    else:
        plt.close(fig)
    return figures


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--analysis", choices=("dual-pole-map", "reversal-map"), required=True)
    result.add_argument("--output", type=Path)
    result.add_argument("--length", type=float, default=10.0)
    result.add_argument("--rho", type=float, default=100.0)
    result.add_argument("--kappa", type=float, default=0.067)
    result.add_argument("--gamma", type=float, default=1.0)
    result.add_argument("--hill", type=float, default=2.0)
    result.add_argument("--mu-min", type=float, default=0.5)
    result.add_argument("--mu-max", type=float, default=5.0)
    result.add_argument("--mu-count", type=int, default=13)
    result.add_argument("--lambda-min", type=float, default=10.0)
    result.add_argument("--lambda-max", type=float, default=100.0)
    result.add_argument("--lambda-count", type=int, default=13)
    result.add_argument("--h-values", default="1,2,4")
    result.add_argument("--asymmetries", default="0,0.001,0.01,0.05")
    result.add_argument("--domain-count-times", default="0,100,300")
    result.add_argument("--initial-amplitude", type=float, default=2.0)
    result.add_argument("--initial-width", type=float, default=0.04)
    result.add_argument("--duration", type=float, default=300.0)
    result.add_argument("--nodes", type=int, default=81)
    result.add_argument("--output-interval", type=float, default=2.0)
    result.add_argument("--max-step", type=float, default=0.5)
    result.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
        help="parallel simulation processes; defaults to the allocated SLURM CPUs",
    )
    result.add_argument("--stimulus-amplitude", type=float, default=1.0)
    result.add_argument("--stimulus-width", type=float, default=0.05)
    result.add_argument("--reversal-start", type=float, default=150.0)
    result.add_argument("--reversal-duration", type=float, default=20.0)
    result.add_argument("--reversal-amplitude", type=float, default=1.0)
    result.add_argument("--reversal-threshold", type=float, default=0.5)
    return result


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    if args.output is None:
        default = (
            "single_species/testing/dual_pole_map"
            if args.analysis == "dual-pole-map"
            else "single_species/testing/polarity_reversal"
        )
        args.output = output_path(default)
    if args.analysis == "dual-pole-map":
        dual_pole_map(args)
    else:
        reversal_map(args)


if __name__ == "__main__":
    main()
