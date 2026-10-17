"""Paired isometric cell-volume scaling screen for single-species wave-pinning."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from dataclasses import replace
from html import escape
import json
import os
from pathlib import Path

import numpy as np

from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from .calibration import STATE_COLORS, plot_extent
from .model import (
    ActivationPulse,
    SingleSpeciesParameters,
    StripWavePinningChemistry,
    homogeneous_equilibria,
)
from .observables import classify_trajectory, polarized_domain_count
from .reporting import STYLE


OUTCOME_ORDER = tuple(STATE_COLORS)
POLARIZED = {"single_polarized", "dual_polarized", "multi_peak"}
TRANSITION_COLORS = {
    "unchanged": "#d9d9d9",
    "single_to_dual": "#009E73",
    "dual_to_single": "#0072B2",
    "gained_polarity": "#56B4E9",
    "lost_polarity": "#D55E00",
    "homogeneous_low_to_high": "#F0E442",
    "homogeneous_high_to_low": "#8C6BB1",
    "other_change": "#CC79A7",
}
TRANSITION_LABELS = {
    "unchanged": "unchanged",
    "single_to_dual": "single polarized → dual polarized",
    "dual_to_single": "dual polarized → single polarized",
    "gained_polarity": "gained polarity",
    "lost_polarity": "lost polarity",
    "homogeneous_low_to_high": "homogeneous low → high",
    "homogeneous_high_to_low": "homogeneous high → low",
    "other_change": "other changed outcome",
}


def values(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def output_root(args) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else args.output.resolve()


def resting_state(parameters: SingleSpeciesParameters) -> tuple[float, float]:
    stable = [item for item in homogeneous_equilibria(parameters) if item["stable"]]
    if not stable:
        raise ValueError("baseline condition has no stable homogeneous equilibrium")
    resting = min(stable, key=lambda item: float(item["active"]))
    return float(resting["active"]), float(resting["inactive"])


def scaled_parameters(
    args, baseline_mu: float, baseline_lambda: float, baseline_h: float,
    volume_multiplier: float,
) -> tuple[SingleSpeciesParameters, float, float, float]:
    baseline = SingleSpeciesParameters(
        total_abundance=baseline_mu,
        dimensionless_length=baseline_lambda,
        cytoplasm_depth=baseline_h,
        diffusion_ratio=args.rho,
        basal_ratio=args.kappa,
        feedback_ratio=args.gamma,
        hill_coefficient=args.hill,
        half_saturation=1.0,
        inactivation_rate=1.0,
        domain_length=args.length,
    )
    active, inactive = resting_state(baseline)
    linear_scale = volume_multiplier ** (1.0 / 3.0)
    scaled_h = linear_scale * baseline_h
    scaled_mu = (active + scaled_h * inactive) / baseline.half_saturation
    scaled = replace(
        baseline,
        total_abundance=scaled_mu,
        dimensionless_length=linear_scale * baseline_lambda,
        cytoplasm_depth=scaled_h,
        domain_length=linear_scale * baseline.domain_length,
    )
    return scaled, active, inactive, linear_scale


def run_condition(task: tuple) -> dict[str, object]:
    args, baseline_mu, baseline_lambda, baseline_h, volume_multiplier = task
    parameters, resting_active, resting_inactive, linear_scale = scaled_parameters(
        args, baseline_mu, baseline_lambda, baseline_h, volume_multiplier
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
            mechanics=MechanicsParameters(reference_length=parameters.domain_length),
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
    profiles = result.membrane[:, 0]
    classified = classify_trajectory(
        result.times, profiles, half_saturation=parameters.half_saturation
    )
    index_100 = int(np.argmin(np.abs(result.times - min(100.0, args.duration))))
    membrane_total, cytoplasm_total = result.compartment_totals()
    conserved = membrane_total[:, 0] + parameters.cytoplasm_depth * cytoplasm_total[:, 0]
    error = float(
        np.max(np.abs(conserved - conserved[0])) / max(abs(conserved[0]), 1.0e-14)
    )
    return {
        "volume_multiplier": volume_multiplier,
        "linear_scale": linear_scale,
        "baseline_mu": baseline_mu,
        "baseline_Lambda": baseline_lambda,
        "baseline_h": baseline_h,
        "resting_active_u": resting_active,
        "resting_inactive_v": resting_inactive,
        "scaled_mu": parameters.total_abundance,
        "scaled_Lambda": parameters.dimensionless_length,
        "scaled_h": parameters.cytoplasm_depth,
        "scaled_length": parameters.domain_length,
        "terminal_outcome": classified["state"],
        "polarized_domain_count_t100": polarized_domain_count(profiles[index_100]),
        "polarized_domain_count_t300": polarized_domain_count(profiles[-1]),
        "terminal_polarization_amplitude": classified["terminal_polarization_amplitude"],
        "late_maximum_normalized_change": classified["late_maximum_normalized_change"],
        "maximum_conservation_error": error,
    }


def all_tasks(args) -> list[tuple]:
    mus = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    return [
        (args, float(mu), float(dimensionless_length), float(h), float(volume))
        for volume in values(args.volume_multipliers)
        for h in values(args.h_values)
        for mu in mus
        for dimensionless_length in lambdas
    ]


def run_shard(args) -> Path:
    root = output_root(args)
    shard_index = int(os.environ.get("SLURM_ARRAY_TASK_ID", args.shard_index))
    if not 0 <= shard_index < args.shard_count:
        raise ValueError("shard index must be smaller than shard count")
    tasks = [
        task for index, task in enumerate(all_tasks(args))
        if index % args.shard_count == shard_index
    ]
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", args.workers))
    if workers <= 1:
        rows = [run_condition(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(run_condition, tasks, chunksize=1))
    shard_dir = root / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)
    destination = shard_dir / f"shard_{shard_index:03d}.csv"
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "shard_index": shard_index,
        "shard_count": args.shard_count,
        "conditions": len(rows),
        "workers": workers,
        "slurm_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
    }
    (shard_dir / f"shard_{shard_index:03d}.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return destination


def cast_row(row: dict[str, str]) -> dict[str, object]:
    numeric = {
        "volume_multiplier", "linear_scale", "baseline_mu", "baseline_Lambda",
        "baseline_h", "resting_active_u", "resting_inactive_v", "scaled_mu",
        "scaled_Lambda", "scaled_h", "scaled_length",
        "terminal_polarization_amplitude", "late_maximum_normalized_change",
        "maximum_conservation_error",
    }
    integer = {"polarized_domain_count_t100", "polarized_domain_count_t300"}
    return {
        key: int(value) if key in integer else float(value) if key in numeric else value
        for key, value in row.items()
    }


def transition_class(source: str, destination: str) -> str:
    if source == destination:
        return "unchanged"
    if source == "single_polarized" and destination == "dual_polarized":
        return "single_to_dual"
    if source == "dual_polarized" and destination == "single_polarized":
        return "dual_to_single"
    if source not in POLARIZED and destination in POLARIZED:
        return "gained_polarity"
    if source in POLARIZED and destination not in POLARIZED:
        return "lost_polarity"
    if source == "homogeneous_low" and destination == "homogeneous_high":
        return "homogeneous_low_to_high"
    if source == "homogeneous_high" and destination == "homogeneous_low":
        return "homogeneous_high_to_low"
    return "other_change"


def paired_lookup(rows: list[dict[str, object]]) -> dict[tuple[float, float, float, float], dict[str, object]]:
    return {
        (
            float(row["volume_multiplier"]), float(row["baseline_h"]),
            float(row["baseline_mu"]), float(row["baseline_Lambda"]),
        ): row
        for row in rows
    }


def render_outcomes(output: Path, rows, mus, lambdas, h_values, volumes) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    lookup = paired_lookup(rows)
    states = list(OUTCOME_ORDER)
    cmap = ListedColormap([STATE_COLORS[state] for state in states])
    fig, axes = plt.subplots(
        len(volumes), len(h_values), figsize=(5.0 * len(h_values), 3.8 * len(volumes)),
        squeeze=False, constrained_layout=True,
    )
    for row_index, volume in enumerate(volumes):
        for column_index, h in enumerate(h_values):
            ax = axes[row_index, column_index]
            grid = np.asarray([
                states.index(str(lookup[(volume, h, float(mu), float(length))]["terminal_outcome"]))
                for mu in mus for length in lambdas
            ]).reshape(len(mus), len(lambdas))
            ax.imshow(
                grid, origin="lower", aspect="auto", cmap=cmap,
                vmin=-0.5, vmax=len(states) - 0.5,
                extent=(*plot_extent(lambdas), *plot_extent(mus)),
            )
            if row_index == 0:
                ax.set_title(f"baseline h = {h:g}")
            ax.set(xlabel="baseline dimensionless length Λ",
                   ylabel=f"{volume:g}× volume\nbaseline abundance μ")
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=STATE_COLORS[state],
                   label=state.replace("_", " ")) for state in states
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("Matched single-pulse outcomes under the concentration-preserving size null")
    path = output / "matched_outcome_maps.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def render_transition_maps(output: Path, rows, mus, lambdas, h_values, destinations) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    lookup = paired_lookup(rows)
    classes = list(TRANSITION_COLORS)
    cmap = ListedColormap([TRANSITION_COLORS[item] for item in classes])
    norm = BoundaryNorm(np.arange(-0.5, len(classes) + 0.5), cmap.N)
    fig, axes = plt.subplots(
        len(destinations), len(h_values),
        figsize=(5.0 * len(h_values), 4.0 * len(destinations)),
        squeeze=False, constrained_layout=True,
    )
    for row_index, destination_volume in enumerate(destinations):
        for column_index, h in enumerate(h_values):
            ax = axes[row_index, column_index]
            grid = []
            for mu in mus:
                for length in lambdas:
                    source = str(lookup[(1.0, h, float(mu), float(length))]["terminal_outcome"])
                    destination = str(lookup[(destination_volume, h, float(mu), float(length))]["terminal_outcome"])
                    grid.append(classes.index(transition_class(source, destination)))
            ax.imshow(
                np.asarray(grid).reshape(len(mus), len(lambdas)),
                origin="lower", aspect="auto", cmap=cmap, norm=norm,
                extent=(*plot_extent(lambdas), *plot_extent(mus)),
            )
            if row_index == 0:
                ax.set_title(f"baseline h = {h:g}")
            ax.set(
                xlabel="baseline dimensionless length Λ",
                ylabel=f"1× → {destination_volume:g}×\nbaseline abundance μ",
            )
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=TRANSITION_COLORS[item],
                   label=TRANSITION_LABELS[item])
        for item in classes
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("Where cell-volume scaling changes the terminal regime")
    path = output / "transition_maps.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def render_transition_matrices(output: Path, rows, destinations) -> Path:
    import matplotlib.pyplot as plt

    lookup = paired_lookup(rows)
    states = list(OUTCOME_ORDER)
    baseline_keys = sorted({key[1:] for key in lookup if key[0] == 1.0})
    fig, axes = plt.subplots(1, len(destinations), figsize=(7.0 * len(destinations), 6.2),
                             squeeze=False, constrained_layout=True)
    for ax, destination_volume in zip(axes[0], destinations):
        matrix = np.zeros((len(states), len(states)), dtype=int)
        for h, mu, length in baseline_keys:
            source = str(lookup[(1.0, h, mu, length)]["terminal_outcome"])
            destination = str(lookup[(destination_volume, h, mu, length)]["terminal_outcome"])
            matrix[states.index(source), states.index(destination)] += 1
        image = ax.imshow(matrix, cmap="Blues")
        for source_index in range(len(states)):
            for destination_index in range(len(states)):
                count = matrix[source_index, destination_index]
                if count:
                    ax.text(destination_index, source_index, str(count), ha="center", va="center",
                            color="white" if count > 0.55 * matrix.max() else "black")
        labels = [state.replace("_", " ") for state in states]
        ax.set(
            xticks=range(len(states)), yticks=range(len(states)),
            xticklabels=labels, yticklabels=labels,
            xlabel=f"outcome at {destination_volume:g}× volume",
            ylabel="outcome at 1× volume",
            title=f"1× → {destination_volume:g}×",
        )
        ax.tick_params(axis="x", rotation=45)
        fig.colorbar(image, ax=ax, label="paired grid conditions")
    fig.suptitle("Exact terminal-outcome transition counts")
    path = output / "transition_matrices.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def render_domain_changes(output: Path, rows, mus, lambdas, h_values, destinations) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    lookup = paired_lookup(rows)
    colors = ["#762a83", "#af8dc3", "#f7f7f7", "#7fbf7b", "#1b7837"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-2.5, 3.0, 1.0), cmap.N)
    fig, axes = plt.subplots(
        len(destinations), len(h_values),
        figsize=(5.0 * len(h_values), 4.0 * len(destinations)),
        squeeze=False, constrained_layout=True,
    )
    image = None
    for row_index, destination_volume in enumerate(destinations):
        for column_index, h in enumerate(h_values):
            ax = axes[row_index, column_index]
            grid = np.asarray([
                int(lookup[(destination_volume, h, float(mu), float(length))]["polarized_domain_count_t300"])
                - int(lookup[(1.0, h, float(mu), float(length))]["polarized_domain_count_t300"])
                for mu in mus for length in lambdas
            ]).reshape(len(mus), len(lambdas))
            image = ax.imshow(
                np.clip(grid, -2, 2), origin="lower", aspect="auto", cmap=cmap, norm=norm,
                extent=(*plot_extent(lambdas), *plot_extent(mus)),
            )
            if row_index == 0:
                ax.set_title(f"baseline h = {h:g}")
            ax.set(xlabel="baseline dimensionless length Λ",
                   ylabel=f"1× → {destination_volume:g}×\nbaseline abundance μ")
    fig.colorbar(image, ax=axes.ravel().tolist(), ticks=(-2, -1, 0, 1, 2),
                 label="change in polarized-domain count at 300 s")
    fig.suptitle("Cell-size-induced gain or loss of terminal polarized domains")
    path = output / "terminal_domain_count_changes.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def render_outcome_counts(output: Path, rows, h_values, volumes) -> Path:
    import matplotlib.pyplot as plt

    labels = [(volume, h) for volume in volumes for h in h_values]
    bottoms = np.zeros(len(labels))
    fig, ax = plt.subplots(figsize=(11.5, 5.2), constrained_layout=True)
    for state in OUTCOME_ORDER:
        counts = np.asarray([
            sum(
                row["volume_multiplier"] == volume and row["baseline_h"] == h
                and row["terminal_outcome"] == state for row in rows
            ) for volume, h in labels
        ])
        if np.any(counts):
            ax.bar(
                range(len(labels)), counts, bottom=bottoms,
                color=STATE_COLORS[state], label=state.replace("_", " "),
            )
            bottoms += counts
    ax.set(
        xticks=range(len(labels)),
        xticklabels=[f"{volume:g}×\nh={h:g}" for volume, h in labels],
        xlabel="volume multiplier and baseline h",
        ylabel="number of (μ, Λ) conditions",
        title="Outcome counts across size and geometry",
    )
    ax.legend(frameon=False, ncol=3)
    path = output / "outcome_counts_by_scale.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def summarize_transitions(rows, destinations) -> dict[str, object]:
    lookup = paired_lookup(rows)
    baseline_keys = sorted({key[1:] for key in lookup if key[0] == 1.0})
    result: dict[str, object] = {}
    for destination_volume in destinations:
        exact: dict[str, int] = {}
        grouped: dict[str, int] = {}
        changed_locations: list[dict[str, object]] = []
        for h, mu, length in baseline_keys:
            source = str(lookup[(1.0, h, mu, length)]["terminal_outcome"])
            destination = str(lookup[(destination_volume, h, mu, length)]["terminal_outcome"])
            exact_key = f"{source} -> {destination}"
            exact[exact_key] = exact.get(exact_key, 0) + 1
            group = transition_class(source, destination)
            grouped[group] = grouped.get(group, 0) + 1
            if source != destination:
                changed_locations.append({
                    "baseline_h": h,
                    "baseline_mu": mu,
                    "baseline_Lambda": length,
                    "source": source,
                    "destination": destination,
                })
        result[f"1x_to_{destination_volume:g}x"] = {
            "grouped": grouped,
            "exact": exact,
            "changed_locations": changed_locations,
        }
    return result


def write_report(output: Path, summary: dict[str, object], figures: list[Path]) -> Path:
    transitions = summary["transitions"]
    transition_blocks = []
    highlights = []
    single_to_dual_locations = []
    for comparison, records in transitions.items():
        exact = records["exact"]
        grouped = records["grouped"]
        changed_total = sum(grouped.values()) - grouped.get("unchanged", 0)
        highlights.append(
            f"<strong>{escape(comparison.replace('_', ' '))}</strong>: "
            f"{changed_total} of 867 paired conditions changed terminal class; "
            f"{grouped.get('single_to_dual', 0)} changed single polarized → dual polarized."
        )
        for location in records["changed_locations"]:
            if (location["source"], location["destination"]) == (
                "single_polarized", "dual_polarized"
            ):
                single_to_dual_locations.append(
                    f"{escape(comparison.replace('_', ' '))} at baseline "
                    f"h={location['baseline_h']:g}, μ={location['baseline_mu']:g}, "
                    f"Λ={location['baseline_Lambda']:g}"
                )
        changed = [(name, count) for name, count in exact.items() if name.split(" -> ")[0] != name.split(" -> ")[1]]
        changed.sort(key=lambda item: item[1], reverse=True)
        rows = "".join(
            f"<tr><td>{escape(name.replace('_', ' '))}</td><td>{count}</td></tr>"
            for name, count in changed
        ) or "<tr><td>No changed outcomes</td><td>0</td></tr>"
        transition_blocks.append(
            f"<div><h3>{escape(comparison.replace('_', ' '))}</h3>"
            f"<table><thead><tr><th>Exact changed transition</th><th>Conditions</th></tr>"
            f"</thead><tbody>{rows}</tbody></table></div>"
        )
    captions = {
        "transition_maps": "Primary result. Each colored pixel says how the 300 s endpoint changes from the same 1× baseline condition after scaling cell volume. Gray pixels do not change class.",
        "transition_matrices": "Counts behind the transition maps. Rows are 1× outcomes and columns are scaled-cell outcomes; off-diagonal cells are changes.",
        "terminal_domain_count_changes": "Difference in the number of disconnected polarized regions at 300 s: scaled cell minus its matched 1× cell. Positive values gain regions; negative values lose them.",
        "matched_outcome_maps": "The underlying endpoint classes at all three volumes. Axes retain the matched 1× baseline coordinates so the same pixel can be compared vertically.",
        "outcome_counts_by_scale": "Aggregate endpoint-class counts. This summarizes prevalence but does not preserve which individual conditions switched class.",
    }
    figure_html = "".join(
        f"<figure><a href='{path.name}'><img src='{path.name}' alt='{escape(path.stem)}' "
        f"loading='lazy'></a><figcaption><strong>{escape(path.stem.replace('_', ' ').title())}</strong>"
        f"<span>{escape(captions[path.stem])}</span></figcaption></figure>"
        for path in figures
    )
    data_file = "wgd_null_screen.csv"
    summary_file = "summary.json"
    extended_path = output / "extended_report_fragment.html"
    extended_html = extended_path.read_text() if extended_path.exists() else ""
    larger_grouped = transitions.get("1x_to_2x", {}).get("grouped", {})
    larger_gained = larger_grouped.get("gained_polarity", 0)
    larger_lost = larger_grouped.get("lost_polarity", 0)
    report = output / "report.html"
    report.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width"><title>Single-species WGD-null scaling report</title>
<style>{STYLE}</style></head><body><main><h1>Single-species WGD-null scaling report</h1>
<p class="lead">A paired single-pulse screen comparing 0.5×, 1×, and 2× cell
volume under an isometric, concentration-preserving null.</p>
<section class="guide"><h2>Model and symbols</h2><p>The model follows one
conserved GTPase in two forms: slowly diffusing active membrane-bound
<code>u(x,t)</code> and rapidly diffusing inactive cytoplasmic <code>v(x,t)</code>.
Position <code>x</code> runs from one end of the one-dimensional cell strip to the
other.</p><div class="table-wrap"><table><thead><tr><th>Symbol</th><th>Meaning</th>
<th>How to read it here</th></tr></thead><tbody>
<tr><td><code>μ</code> (mu)</td><td>Total normalized GTPase abundance</td><td><code>μ=mean(u+h v)/K</code>. Moving upward in a map means more total GTPase in the baseline cell.</td></tr>
<tr><td><code>Λ</code> (capital lambda)</td><td>Length relative to the active-GTPase diffusion distance</td><td><code>Λ=L√(δ/D<sub>u</sub>)</code>. Moving right means a longer effective domain. During size scaling, physical <code>L</code> and therefore <code>Λ</code> both scale.</td></tr>
<tr><td><code>h</code></td><td>Cytoplasmic reservoir per membrane area</td><td>The strip-model conversion between cytoplasmic and membrane amount. Larger <code>h</code> gives each membrane segment more adjacent cytoplasmic volume.</td></tr>
<tr><td><code>K</code></td><td>Half-feedback concentration</td><td>The active concentration at which positive feedback is half maximal; <code>K=1</code>.</td></tr>
<tr><td><code>V</code>, <code>s</code></td><td>Volume and linear multipliers</td><td><code>V</code> is 0.5, 1, or 2; isometric scaling gives <code>s=V<sup>1/3</sup></code>.</td></tr>
</tbody></table></div></section>
<section class="guide"><h2>Experiment</h2><p>Each of the 867 baseline
<code>(μ, Λ, h)</code> grid points is followed through three matched cell volumes.
For volume multiplier <code>V</code>, the linear scale is
<code>s=V<sup>1/3</sup></code>. The transformation is
<code>L→sL</code>, <code>Λ→sΛ</code>, and <code>h→sh</code>. Baseline resting
concentrations <code>u₀</code> and <code>v₀</code> are preserved, so
<code>μ→(u₀+shv₀)/K</code>. Diffusivities, reaction kinetics, cue amplitude and
timing, and the cue's fractional position and width are unchanged.</p>
<p><strong>Common protocol:</strong> all runs use <code>ρ=100</code>,
<code>κ=0.067</code>, <code>Γ=1</code>, and Hill coefficient <code>n=2</code>.
The starting state is the lowest stable homogeneous equilibrium. A Gaussian
activation cue of amplitude 1, centered 10% along the strip with width 5%, acts
from 10 to 30 s. The no-flux strip uses 81 finite-volume nodes and is followed to
300 s. Thus material cannot leave through either end.</p>
<p>The axes in all maps are the <strong>baseline 1× coordinates</strong>. Thus the
same pixel always refers to the same paired molecular condition before geometric
scaling.</p></section>
<section class="guide"><h2>What the endpoint labels mean</h2>
<p>At 300 s, polarization contrast is <code>P=(max(u)−min(u))/mean(u)</code>. A
<strong>polarized region</strong> is one continuous interval along the strip where
<code>u</code> is above the midpoint between that profile's minimum and maximum.
Separate above-midpoint intervals are separate regions. Tiny ripples are rejected:
the profile must also have <code>P≥0.5</code>. The final 50 s must have normalized
profile change ≤0.01 to receive a stationary label.</p>
<div class="table-wrap"><table><thead><tr><th>Label</th><th>Operational meaning</th><th>Interpretation</th></tr></thead><tbody>
<tr><td><strong>homogeneous low</strong></td><td><code>P&lt;0.5</code> and mean <code>u&lt;K</code></td><td>No substantial spatial polarity; membrane activity is broadly low.</td></tr>
<tr><td><strong>homogeneous high</strong></td><td><code>P&lt;0.5</code> and mean <code>u≥K</code></td><td>Whole-strip activation, not a pole.</td></tr>
<tr><td><strong>single polarized</strong></td><td><code>P≥0.5</code> and one connected high interval</td><td>One persistent active region.</td></tr>
<tr><td><strong>dual polarized</strong></td><td><code>P≥0.5</code> and two disconnected high intervals</td><td>Two regions persist to 300 s; this is finite-time persistence, not proof of mathematical stability.</td></tr>
<tr><td><strong>multi peak</strong></td><td><code>P≥0.5</code> and more than two high intervals</td><td>Several persistent activity sites.</td></tr>
<tr><td><strong>transient coarsening</strong></td><td>Still changing, with fewer regions during the last 50 s</td><td>Sites are merging or competing.</td></tr>
<tr><td><strong>oscillatory/dynamic</strong></td><td>Still changing without late loss of regions</td><td>No stationary endpoint by 300 s.</td></tr>
</tbody></table></div>
<div class="profile-example"><h3>One region versus two</h3><p>The orange dashed
line is the profile-specific midpoint. Blue shading marks disconnected intervals
above it; these are what the analysis counts.</p>
<svg viewBox="0 0 760 300" role="img" aria-label="Profiles illustrating one and two polarized regions">
<g transform="translate(25 10)"><text x="150" y="20" text-anchor="middle" class="svg-title">one connected polarized region</text><line x1="25" y1="245" x2="300" y2="245" class="axis"/><line x1="25" y1="45" x2="25" y2="245" class="axis"/><rect x="55" y="55" width="85" height="190" class="domain"/><line x1="25" y1="150" x2="300" y2="150" class="threshold"/><path d="M25 225 C45 225 48 190 58 135 C70 70 88 52 101 52 C118 52 132 90 142 155 C150 210 170 225 300 225" class="profile"/><text x="98" y="275" text-anchor="middle">region 1</text></g>
<g transform="translate(405 10)"><text x="150" y="20" text-anchor="middle" class="svg-title">two disconnected polarized regions</text><line x1="25" y1="245" x2="300" y2="245" class="axis"/><line x1="25" y1="45" x2="25" y2="245" class="axis"/><rect x="45" y="55" width="65" height="190" class="domain"/><rect x="210" y="55" width="65" height="190" class="domain"/><line x1="25" y1="150" x2="300" y2="150" class="threshold"/><path d="M25 225 C38 220 42 160 50 120 C60 66 75 52 85 52 C100 52 108 100 114 165 C120 218 145 225 170 225 C195 225 202 165 212 120 C222 66 238 52 250 52 C266 52 275 105 282 170 C288 215 294 225 300 225" class="profile"/><text x="78" y="275" text-anchor="middle">region 1</text><text x="243" y="275" text-anchor="middle">region 2</text></g>
</svg></div></section>
<section class="guide"><h2>How to read the transition maps</h2><p>Gray means the
terminal outcome is unchanged. Green specifically identifies
<strong>single polarized → dual polarized</strong>. Blue is the reverse. Cyan and
orange indicate gain and loss of polarity, respectively. Yellow and purple are
changes between homogeneous-low and homogeneous-high states. Pink collects less
common changes involving other endpoint classes.</p></section>
<aside class="result">{'<br>'.join(highlights)}
{'<br><strong>Single-to-dual location:</strong> ' + '; '.join(single_to_dual_locations) if single_to_dual_locations else ''}</aside>
<aside class="note"><strong>Interpretation for pole robustness.</strong> At 2×
volume, {larger_gained} previously homogeneous conditions gain polarity and
{larger_lost} polarized conditions lose it. This supports more robust
establishment or maintenance of spatial polarity in the larger-cell null. It
does not by itself imply more reliable monopolarity, because some transitions
change the number of persistent regions.</aside>
<section class="analysis"><h2>Exact changed transitions</h2><div class="figures">
{''.join(transition_blocks)}</div></section>
<div class="figures">{figure_html}</div>
{extended_html}
<section class="guide"><h2>Scope</h2><p>This is a geometric null, not a fitted
model of a particular WGD cell line. A terminal state means the classification at
300 s. The paired design isolates changes caused by isometric cell-volume scaling
while retaining local concentration-level chemistry.</p></section>
<h2>Downloads</h2><ul class="files"><li><a href="{data_file}">{data_file}</a></li>
<li><a href="{summary_file}">{summary_file}</a></li></ul>
</main></body></html>""")
    return report


def merge_shards(args) -> Path:
    root = output_root(args)
    shard_paths = [root / "shards" / f"shard_{index:03d}.csv" for index in range(args.shard_count)]
    missing = [str(path) for path in shard_paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing WGD shards: {missing}")
    rows: list[dict[str, object]] = []
    for path in shard_paths:
        with path.open(newline="") as handle:
            rows.extend(cast_row(row) for row in csv.DictReader(handle))
    expected = args.mu_count * args.lambda_count * len(values(args.h_values)) * len(values(args.volume_multipliers))
    keys = {
        (row["volume_multiplier"], row["baseline_h"], row["baseline_mu"], row["baseline_Lambda"])
        for row in rows
    }
    if len(rows) != expected or len(keys) != expected:
        raise ValueError(f"expected {expected} unique rows, found {len(rows)} rows and {len(keys)} keys")
    rows.sort(key=lambda row: (
        row["volume_multiplier"], row["baseline_h"], row["baseline_mu"], row["baseline_Lambda"]
    ))
    root.mkdir(parents=True, exist_ok=True)
    with (root / "wgd_null_screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    mus = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    h_values = values(args.h_values)
    volumes = values(args.volume_multipliers)
    destinations = [volume for volume in volumes if not np.isclose(volume, 1.0)]
    figures = [
        render_transition_maps(root, rows, mus, lambdas, h_values, destinations),
        render_transition_matrices(root, rows, destinations),
        render_domain_changes(root, rows, mus, lambdas, h_values, destinations),
        render_outcomes(root, rows, mus, lambdas, h_values, volumes),
        render_outcome_counts(root, rows, h_values, volumes),
    ]
    summary = {
        "analysis": "single_species_wgd_null_screen",
        "volume_multipliers": volumes,
        "linear_scales": {f"{volume:g}x": volume ** (1.0 / 3.0) for volume in volumes},
        "baseline_grid": {
            "mu": mus.tolist(), "Lambda": lambdas.tolist(), "h": h_values,
            "conditions": args.mu_count * args.lambda_count * len(h_values),
        },
        "total_simulations": len(rows),
        "null": {
            "fixed": ["resting u concentration", "resting v concentration", "kinetics",
                      "diffusivities", "stimulus amplitude", "stimulus timing",
                      "stimulus fractional position and width"],
            "scaled": ["L", "Lambda", "h"],
            "mu_rule": "(u0 + scaled_h * v0) / K",
        },
        "transitions": summarize_transitions(rows, destinations),
        "maximum_conservation_error": max(float(row["maximum_conservation_error"]) for row in rows),
        "figures": [path.name for path in figures],
    }
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return write_report(root, summary, figures)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--mode", choices=("shard", "merge"), required=True)
    result.add_argument("--output", type=Path, default=Path("migration1d/results/single_species/testing/wgd_null"))
    result.add_argument("--volume-multipliers", default="0.5,1,2")
    result.add_argument("--shard-count", type=int, default=32)
    result.add_argument("--shard-index", type=int, default=0)
    result.add_argument("--workers", type=int, default=1)
    result.add_argument("--mu-min", type=float, default=0.5)
    result.add_argument("--mu-max", type=float, default=5.0)
    result.add_argument("--mu-count", type=int, default=17)
    result.add_argument("--lambda-min", type=float, default=10.0)
    result.add_argument("--lambda-max", type=float, default=100.0)
    result.add_argument("--lambda-count", type=int, default=17)
    result.add_argument("--h-values", default="1,2,4")
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
    result = run_shard(args) if args.mode == "shard" else merge_shards(args)
    print(result)


if __name__ == "__main__":
    main()
