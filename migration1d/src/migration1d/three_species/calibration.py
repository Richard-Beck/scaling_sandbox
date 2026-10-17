"""Tier-1 Holmes screen over size, basal Rac activation, and cooperativity."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import replace
import json
import os
from pathlib import Path

import numpy as np

from migration1d.holmes import HolmesChemistry
from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from migration1d.single_species.calibration import STATE_COLORS
from migration1d.single_species.observables import classify_trajectory
from migration1d.single_species.reporting import Figure, write_figure_manifest
from .model import (
    GaussianRacPulse,
    canonical_holmes_parameters,
    homogeneous_equilibrium,
    homogeneous_stability,
)
from .observables import classify_three_species_trajectory


def output_path() -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    if configured:
        return Path(configured).resolve()
    return Path.cwd() / "results" / "three_species" / "calibration" / "tier1_regime_map"


def parse_values(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def endpoint_contrast(profile: np.ndarray) -> float:
    return float(
        (profile[-1] - profile[0]) / max(abs(profile[-1]) + abs(profile[0]), 1.0e-12)
    )


def plot_extent(values: np.ndarray) -> tuple[float, float]:
    if len(values) > 1:
        return float(values[0]), float(values[-1])
    delta = max(1.0e-3, 0.05 * max(abs(float(values[0])), 1.0))
    return float(values[0] - delta), float(values[0] + delta)


def run_condition(args, dimensionless_length: float, kappa_r: float, hill: float) -> dict[str, object]:
    parameters = canonical_holmes_parameters(kappa_r=kappa_r, hill_exponent=hill)
    active = homogeneous_equilibrium(parameters)
    parameters = replace(parameters, initial_active=tuple(float(value) for value in active))
    # d=1/s and D_m=0.1 um2/s are fixed, so Lambda is changed through physical L.
    length = dimensionless_length * np.sqrt(parameters.active_diffusivity)
    stability = homogeneous_stability(parameters, active, length)
    chemistry = HolmesChemistry(
        parameters,
        cue=GaussianRacPulse(
            amplitude=args.stimulus_amplitude,
            center_fraction=args.stimulus_center,
            width_fraction=args.stimulus_width,
            start=args.stimulus_start,
            stop=args.stimulus_stop,
        ),
    )
    result = run_simulation(
        MigrationModel(
            mechanics=MechanicsParameters(reference_length=length),
            chemistry=chemistry,
        ),
        SimulationConfig(
            nodes=args.nodes,
            duration=args.duration,
            output_interval=args.output_interval,
            max_step=args.max_step,
            noise_amplitude=0.0,
            mechanics_activation_time=None,
            jacobian_sparsity=True,
        ),
    )
    composite = classify_three_species_trajectory(
        result.times,
        result.membrane,
        parameters.totals,
        homogeneous_growth_rate=float(stability["homogeneous_growth_rate"]),
        maximum_spatial_growth_rate=float(stability["maximum_spatial_growth_rate"]),
        polarity_threshold=args.polarity_threshold,
    )
    composite_status = composite.pop("terminal_polarity_status")
    composite_change = composite.pop("late_maximum_normalized_change")
    rac_classified = classify_trajectory(
        result.times, result.membrane[:, 1],
        half_saturation=parameters.rac_inhibits_rho_half_level,
    )
    active_total, inactive_total = result.compartment_totals()
    conserved = active_total + inactive_total
    relative_error = np.max(
        np.abs(conserved - conserved[0]) / np.maximum(np.abs(conserved[0]), 1.0e-12),
        axis=0,
    )
    terminal = result.membrane[-1]
    return {
        "Lambda": dimensionless_length,
        "physical_length_um": length,
        "kappa_R": kappa_r,
        "I_R1_uM_per_s": parameters.rac_basal_activation,
        "n": hill,
        **stability,
        "terminal_polarity_status": rac_classified["state"],
        "time_to_polarization": rac_classified["time_to_polarization"],
        "terminal_polarization_amplitude": rac_classified["terminal_polarization_amplitude"],
        "terminal_pole_count": rac_classified["terminal_pole_count"],
        "terminal_high_width_fraction": rac_classified["terminal_high_width_fraction"],
        "late_maximum_normalized_change": rac_classified["late_maximum_normalized_change"],
        "terminal_composite_response_status": composite_status,
        "composite_late_maximum_normalized_change": composite_change,
        **composite,
        "terminal_cdc42_endpoint_contrast": endpoint_contrast(terminal[0]),
        "terminal_rac_endpoint_contrast": endpoint_contrast(terminal[1]),
        "terminal_rho_endpoint_contrast": endpoint_contrast(terminal[2]),
        "maximum_cdc42_conservation_error": float(relative_error[0]),
        "maximum_rac_conservation_error": float(relative_error[1]),
        "maximum_rho_conservation_error": float(relative_error[2]),
    }


def run_task(task: tuple) -> dict[str, object]:
    return run_condition(*task)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict[str, object]]:
    with path.open(newline="") as handle:
        raw = list(csv.DictReader(handle))
    numeric = set(raw[0]) - {
        "terminal_polarity_status", "terminal_composite_response_status",
        "coherent_polarity",
    }
    for row in raw:
        for key in numeric:
            row[key] = None if row[key] in {"", "None"} else float(row[key])
        row["coherent_polarity"] = row["coherent_polarity"] == "True"
    return raw


def render_maps(
    output: Path,
    rows: list[dict[str, object]],
    lambdas: np.ndarray,
    kappas: np.ndarray,
    hills: list[float],
) -> list[Figure]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, TwoSlopeNorm

    states = list(STATE_COLORS)
    cmap = ListedColormap([STATE_COLORS[state] for state in states])
    fig, axes = plt.subplots(1, len(hills), figsize=(5.3 * len(hills), 4.9),
                             squeeze=False, constrained_layout=True)
    for ax, hill in zip(axes[0], hills):
        subset = [row for row in rows if float(row["n"]) == hill]
        grid = np.asarray([
            states.index(str(row["terminal_polarity_status"])) for row in subset
        ]).reshape(len(kappas), len(lambdas))
        ax.imshow(
            grid, origin="lower", aspect="auto", cmap=cmap,
            vmin=-0.5, vmax=len(states) - 0.5,
            extent=(*plot_extent(lambdas), *plot_extent(kappas)),
        )
        ax.scatter([20.0 / np.sqrt(0.1)], [0.2 / 7.5], marker="*", s=90,
                   facecolor="white", edgecolor="black", linewidth=0.7)
        ax.set(title=f"Hill coefficient n = {hill:g}", xlabel="dimensionless length Λ",
               ylabel="basal Rac activation κ_R")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=color,
                          label=state.replace("_", " "))
               for state, color in STATE_COLORS.items()]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("Holmes six-field response to a transient local Rac cue")
    fig.savefig(output / "tier1_outcome_map.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(hills), figsize=(5.3 * len(hills), 4.5),
                             squeeze=False, constrained_layout=True)
    all_growth = np.asarray([
        max(float(row["homogeneous_growth_rate"]),
            float(row["maximum_spatial_growth_rate"])) for row in rows
    ])
    bound = max(0.02, float(np.nanpercentile(np.abs(all_growth), 95)))
    image = None
    for ax, hill in zip(axes[0], hills):
        subset = [row for row in rows if float(row["n"]) == hill]
        grid = np.asarray([
            max(float(row["homogeneous_growth_rate"]),
                float(row["maximum_spatial_growth_rate"])) for row in subset
        ]).reshape(len(kappas), len(lambdas))
        image = ax.imshow(
            grid, origin="lower", aspect="auto", cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound),
            extent=(*plot_extent(lambdas), *plot_extent(kappas)),
        )
        if len(lambdas) > 1 and len(kappas) > 1 and np.min(grid) <= 0.0 <= np.max(grid):
            ax.contour(lambdas, kappas, grid, levels=[0.0], colors="black", linewidths=1.0)
        ax.set(title=f"n = {hill:g}", xlabel="dimensionless length Λ",
               ylabel="basal Rac activation κ_R")
    fig.colorbar(image, ax=axes, label="largest linear growth rate (1/s)", shrink=0.85)
    fig.suptitle("Homogeneous and spatial linear stability")
    fig.savefig(output / "tier1_linear_stability.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.7), constrained_layout=True)
    bottoms = np.zeros(len(hills))
    for state in states:
        counts = np.asarray([
            sum(float(row["n"]) == hill and row["terminal_polarity_status"] == state
                for row in rows) for hill in hills
        ])
        if np.any(counts):
            ax.bar([f"n = {hill:g}" for hill in hills], counts, bottom=bottoms,
                   color=STATE_COLORS[state], label=state.replace("_", " "))
            bottoms += counts
    ax.set(ylabel="number of sampled (Λ, κ_R) conditions",
           title="Tier-1 outcome counts by cooperativity")
    ax.legend(frameon=False, ncol=2)
    fig.savefig(output / "tier1_outcome_counts.png", dpi=200)
    plt.close(fig)
    return [
        Figure("tier1_outcome_map.png", "Stimulated-response regime map",
               "Late outcome across physical size, basal Rac activation, and Hill cooperativity; the star marks the Holmes baseline."),
        Figure("tier1_linear_stability.png", "Linear-stability map",
               "Largest homogeneous or nonzero spatial-mode growth rate. Positive values identify spontaneous/hypersensitive regimes."),
        Figure("tier1_outcome_counts.png", "Outcome counts by cooperativity",
               "Counts summarize the deliberately uniform sampled grid and are not biological probabilities."),
    ]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--render-existing", action="store_true")
    parser.add_argument("--lambda-min", type=float, default=31.6227766017)
    parser.add_argument("--lambda-max", type=float, default=94.8683298051)
    parser.add_argument("--lambda-count", type=int, default=13)
    parser.add_argument("--kappa-min", type=float, default=0.0)
    parser.add_argument("--kappa-max", type=float, default=0.08)
    parser.add_argument("--kappa-count", type=int, default=13)
    parser.add_argument("--hill-values", default="2,3,4")
    parser.add_argument("--duration", type=float, default=300.0)
    parser.add_argument("--nodes", type=int, default=81)
    parser.add_argument("--output-interval", type=float, default=2.0)
    parser.add_argument("--max-step", type=float, default=0.5)
    parser.add_argument("--polarity-threshold", type=float, default=0.2)
    parser.add_argument("--stimulus-amplitude", type=float, default=1.0)
    parser.add_argument("--stimulus-center", type=float, default=0.1)
    parser.add_argument("--stimulus-width", type=float, default=0.05)
    parser.add_argument("--stimulus-start", type=float, default=10.0)
    parser.add_argument("--stimulus-stop", type=float, default=30.0)
    parser.add_argument(
        "--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    )
    args = parser.parse_args(argv)
    output = (args.output or output_path()).resolve()
    output.mkdir(parents=True, exist_ok=True)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    kappas = np.linspace(args.kappa_min, args.kappa_max, args.kappa_count)
    hills = parse_values(args.hill_values)
    if args.render_existing:
        rows = _read_csv(output / "tier1_regime_map.csv")
        summary = json.loads((output / "summary.json").read_text())
        lambdas = np.asarray(summary["axes"]["Lambda"], dtype=float)
        kappas = np.asarray(summary["axes"]["kappa_R"], dtype=float)
        hills = [float(value) for value in summary["axes"]["n"]]
    else:
        tasks = [
            (args, float(length), float(kappa), float(hill))
            for hill in hills for kappa in kappas for length in lambdas
        ]
        if args.workers <= 1:
            rows = [run_task(task) for task in tasks]
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                rows = list(executor.map(run_task, tasks, chunksize=1))
        _write_csv(output / "tier1_regime_map.csv", rows)
        counts: dict[str, int] = {}
        for row in rows:
            state = str(row["terminal_polarity_status"])
            counts[state] = counts.get(state, 0) + 1
        baseline_lambda = 20.0 / np.sqrt(0.1)
        wgd_lambda = baseline_lambda * 2.0 ** (1.0 / 3.0)
        matched = [
            run_condition(args, baseline_lambda, 0.2 / 7.5, 3.0),
            run_condition(args, wgd_lambda, 0.2 / 7.5, 3.0),
        ]
        matched[0]["geometry"] = "1x_volume"
        matched[1]["geometry"] = "2x_volume_isometric"
        _write_csv(output / "matched_wgd_null.csv", matched)
        summary = {
            "analysis": "three_species_tier1_regime_map",
            "model": "Holmes six-field GTPase-only reduction",
            "axes": {"Lambda": lambdas.tolist(), "kappa_R": kappas.tolist(), "n": hills},
            "fixed": {
                "D_m_um2_per_s": 0.1, "D_i_um2_per_s": 50.0,
                "diffusion_ratio": 500.0, "inactivation_per_s": 1.0,
                "C_total_uM": 2.4, "R_total_uM": 7.5, "rho_total_uM": 3.1,
                "I_C_uM_per_s": 2.95, "I_rho_uM_per_s": 6.6,
                "alpha_per_s": 0.65, "a1_uM": 1.25, "a2_uM": 1.0,
            },
            "stimulus": {
                "target": "Rac activation", "shape": "Gaussian in relative coordinate",
                "amplitude_uM_per_s": args.stimulus_amplitude,
                "center_fraction": args.stimulus_center,
                "width_fraction": args.stimulus_width,
                "start_seconds": args.stimulus_start, "stop_seconds": args.stimulus_stop,
            },
            "classification": {
                "Cdc42_and_Rac_range_contrast_threshold": args.polarity_threshold,
                "minimum_Cdc42_Rac_correlation": 0.5,
                "maximum_Cdc42_Rho_correlation": 0.0,
                "minimum_sustained_seconds": 10.0,
                "late_window_seconds": 50.0,
                "stationary_change_threshold": 0.01,
                "linear_instability_threshold_per_s": 1.0e-5,
            },
            "numerics": {
                "duration_seconds": args.duration, "nodes": args.nodes,
                "output_interval_seconds": args.output_interval,
                "max_step_seconds": args.max_step, "workers": args.workers,
            },
            "state_counts": counts,
            "matched_wgd_null": {
                "linear_scale": 2.0 ** (1.0 / 3.0),
                "Lambda_1x": baseline_lambda, "Lambda_2x_volume": wgd_lambda,
                "kinetics_and_concentrations": "fixed",
            },
        }
        (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = render_maps(output, rows, lambdas, kappas, hills)
    write_figure_manifest(
        output,
        title="Three-species Holmes Tier-1 parameter screen",
        introduction=(
            "A transient local Rac cue is screened over physical domain size, basal "
            "Rac activation, and inhibitory cooperativity. Linear stability is "
            "evaluated beside each PDE response so finite polarizability is separated "
            "from spontaneous symmetry breaking."
        ),
        figures=figures,
        data_files=("summary.json", "tier1_regime_map.csv", "matched_wgd_null.csv", "figure_manifest.json"),
    )


if __name__ == "__main__":
    main()
