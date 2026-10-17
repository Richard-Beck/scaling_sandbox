"""Paired isometric volume-scaling screen for the six-field Holmes model."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from html import escape
import json
import os
from pathlib import Path

import numpy as np

from migration1d.single_species.reporting import STYLE, write_figure_manifest
from migration1d.single_species.wgd import (
    TRANSITION_COLORS,
    transition_class as transition_group,
)
from .calibration import (
    STATE_COLORS, parse_values, plot_extent, render_maps, run_condition,
)


OUTCOMES = tuple(STATE_COLORS)


def output_root(args) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else args.output.resolve()


def all_tasks(args) -> list[tuple]:
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    kappas = np.linspace(args.kappa_min, args.kappa_max, args.kappa_count)
    return [
        (args, float(length), float(kappa), float(hill), float(volume))
        for volume in parse_values(args.volume_multipliers)
        for hill in parse_values(args.hill_values)
        for kappa in kappas
        for length in lambdas
    ]


def task_key(task: tuple) -> tuple[float, float, float, float]:
    _, length, kappa, hill, volume = task
    return volume, hill, kappa, length


def run_scaled_condition(task: tuple) -> dict[str, object]:
    args, baseline_lambda, kappa, hill, volume = task
    scale = volume ** (1.0 / 3.0)
    row = run_condition(args, baseline_lambda * scale, kappa, hill)
    return {
        "volume_multiplier": volume,
        "linear_scale": scale,
        "baseline_Lambda": baseline_lambda,
        "scaled_Lambda": row.pop("Lambda"),
        "baseline_physical_length_um": baseline_lambda * np.sqrt(0.1),
        "scaled_physical_length_um": row.pop("physical_length_um"),
        **row,
    }


def run_array_task(args) -> Path:
    tasks = all_tasks(args)
    index = int(os.environ.get("SLURM_ARRAY_TASK_ID", args.task_index))
    if not 0 <= index < len(tasks):
        raise ValueError(f"task index {index} outside 0..{len(tasks) - 1}")
    task = tasks[index]
    destination_dir = output_root(args) / "tasks"
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / f"task_{index:05d}.json"
    if destination.exists():
        existing = json.loads(destination.read_text())
        if existing.get("task_index") == index and tuple(existing.get("task_key", ())) == task_key(task):
            return destination
    row = run_scaled_condition(task)
    payload = {"task_index": index, "task_key": task_key(task), "row": row}
    destination.write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    return destination


def run_shard(args) -> Path:
    root = output_root(args)
    shard_index = int(os.environ.get("SLURM_ARRAY_TASK_ID", args.shard_index))
    if not 0 <= shard_index < args.shard_count:
        raise ValueError("shard index must be smaller than shard count")
    indexed = [
        (index, task) for index, task in enumerate(all_tasks(args))
        if index % args.shard_count == shard_index
    ]
    if not indexed:
        raise ValueError(f"shard {shard_index} has no assigned conditions")
    shard_dir = root / "shards"
    destination = shard_dir / f"shard_{shard_index:03d}.csv"
    metadata_path = shard_dir / f"shard_{shard_index:03d}.json"
    if destination.exists() and metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        with destination.open(newline="") as handle:
            existing = [cast_row(row) for row in csv.DictReader(handle)]
        expected_indices = [index for index, _ in indexed]
        expected_keys = {task_key(task) for _, task in indexed}
        if (
            metadata.get("shard_index") == shard_index
            and metadata.get("shard_count") == args.shard_count
            and metadata.get("global_indices") == expected_indices
            and metadata.get("conditions") == len(indexed)
            and {row_key(row) for row in existing} == expected_keys
        ):
            return destination
    workers = int(os.environ.get("SLURM_CPUS_PER_TASK", args.workers))
    tasks = [task for _, task in indexed]
    if workers <= 1:
        rows = [run_scaled_condition(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            rows = list(executor.map(run_scaled_condition, tasks, chunksize=1))
    shard_dir.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "shard_index": shard_index,
        "shard_count": args.shard_count,
        "conditions": len(rows),
        "global_indices": [index for index, _ in indexed],
        "workers": workers,
        "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
    }
    metadata_path.write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    return destination


def cast_row(row: dict[str, str]) -> dict[str, object]:
    integers = {"fastest_spatial_mode", "rac_pole_count"}
    strings = {"terminal_polarity_status", "terminal_composite_response_status"}
    result: dict[str, object] = {}
    for key, value in row.items():
        if key in strings:
            result[key] = value
        elif key == "coherent_polarity":
            if value not in {"True", "False"}:
                raise ValueError(f"invalid boolean {value!r} in {key}")
            result[key] = value == "True"
        elif value in {"", "None"}:
            result[key] = None
        elif key in integers:
            result[key] = int(float(value))
        else:
            result[key] = float(value)
    return result


def row_key(row: dict[str, object]) -> tuple[float, float, float, float]:
    return (
        float(row["volume_multiplier"]), float(row["n"]),
        float(row["kappa_R"]), float(row["baseline_Lambda"]),
    )


def paired_lookup(rows: list[dict[str, object]]) -> dict[tuple[float, float, float, float], dict[str, object]]:
    return {row_key(row): row for row in rows}


def exact_transitions(rows, destinations) -> list[str]:
    lookup = paired_lookup(rows)
    baseline_keys = sorted(key[1:] for key in lookup if np.isclose(key[0], 1.0))
    observed = set()
    for destination_volume in destinations:
        for hill, kappa, length in baseline_keys:
            source = lookup[(1.0, hill, kappa, length)]["terminal_polarity_status"]
            destination = lookup[(destination_volume, hill, kappa, length)]["terminal_polarity_status"]
            observed.add(f"{source} -> {destination}")
    return sorted(observed, key=lambda value: (value.split(" -> ")[0] != value.split(" -> ")[1], value))


def render_exact_transition_maps(output, rows, lambdas, kappas, hills, destinations) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    lookup = paired_lookup(rows)
    transitions = exact_transitions(rows, destinations)
    colors = plt.get_cmap("tab20", max(len(transitions), 1))
    cmap = ListedColormap([colors(index) for index in range(len(transitions))])
    fig, axes = plt.subplots(
        len(destinations), len(hills),
        figsize=(5.2 * len(hills), 4.2 * len(destinations)),
        squeeze=False, constrained_layout=True,
    )
    for row_index, volume in enumerate(destinations):
        for column_index, hill in enumerate(hills):
            ax = axes[row_index, column_index]
            grid = []
            for kappa in kappas:
                for length in lambdas:
                    source = lookup[(1.0, hill, float(kappa), float(length))]["terminal_polarity_status"]
                    destination = lookup[(volume, hill, float(kappa), float(length))]["terminal_polarity_status"]
                    grid.append(transitions.index(f"{source} -> {destination}"))
            ax.imshow(
                np.asarray(grid).reshape(len(kappas), len(lambdas)),
                origin="lower", aspect="auto", cmap=cmap,
                vmin=-0.5, vmax=len(transitions) - 0.5,
                extent=(*plot_extent(lambdas), *plot_extent(kappas)),
            )
            if row_index == 0:
                ax.set_title(f"n = {hill:g}")
            ax.set(
                xlabel="baseline dimensionless length Λ",
                ylabel=f"1× → {volume:g}× volume\nbaseline κ_R",
            )
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=cmap(index),
                   label=value.replace("_", " "))
        for index, value in enumerate(transitions)
    ]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Exact Holmes terminal-outcome transitions under isometric scaling")
    path = output / "exact_transition_maps.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def render_transition_matrices(output, rows, hills, destinations) -> Path:
    import matplotlib.pyplot as plt

    lookup = paired_lookup(rows)
    fig, axes = plt.subplots(
        len(destinations), len(hills),
        figsize=(5.8 * len(hills), 5.4 * len(destinations)),
        squeeze=False, constrained_layout=True,
    )
    for row_index, volume in enumerate(destinations):
        for column_index, hill in enumerate(hills):
            ax = axes[row_index, column_index]
            keys = sorted(key[2:] for key in lookup if np.isclose(key[0], 1.0) and key[1] == hill)
            matrix = np.zeros((len(OUTCOMES), len(OUTCOMES)), dtype=int)
            for kappa, length in keys:
                source = str(lookup[(1.0, hill, kappa, length)]["terminal_polarity_status"])
                destination = str(lookup[(volume, hill, kappa, length)]["terminal_polarity_status"])
                matrix[OUTCOMES.index(source), OUTCOMES.index(destination)] += 1
            image = ax.imshow(matrix, cmap="Blues", vmin=0)
            limit = max(int(matrix.max()), 1)
            for source_index, destination_index in np.argwhere(matrix > 0):
                count = matrix[source_index, destination_index]
                ax.text(destination_index, source_index, str(count), ha="center", va="center",
                        color="white" if count > 0.55 * limit else "black")
            labels = [state.replace("_", " ") for state in OUTCOMES]
            ax.set(
                xticks=range(len(OUTCOMES)), yticks=range(len(OUTCOMES)),
                xticklabels=labels, yticklabels=labels,
                xlabel=f"outcome at {volume:g}× volume",
                ylabel="outcome at 1× volume",
                title=f"n={hill:g}; 1× → {volume:g}×",
            )
            ax.tick_params(axis="x", rotation=55, labelsize=8)
            ax.tick_params(axis="y", labelsize=8)
            fig.colorbar(image, ax=ax, label="paired conditions", shrink=0.75)
    fig.suptitle("Exact three-species transition-count matrices")
    path = output / "exact_transition_count_matrices.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def render_grouped_transition_maps(output, rows, lambdas, kappas, hills, destinations) -> Path:
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    lookup = paired_lookup(rows)
    groups = list(TRANSITION_COLORS)
    cmap = ListedColormap([TRANSITION_COLORS[group] for group in groups])
    fig, axes = plt.subplots(
        len(destinations), len(hills),
        figsize=(5.2 * len(hills), 4.2 * len(destinations)),
        squeeze=False, constrained_layout=True,
    )
    for row_index, volume in enumerate(destinations):
        for column_index, hill in enumerate(hills):
            grid = []
            for kappa in kappas:
                for length in lambdas:
                    source = str(lookup[(1.0, hill, float(kappa), float(length))]["terminal_polarity_status"])
                    destination = str(lookup[(volume, hill, float(kappa), float(length))]["terminal_polarity_status"])
                    grid.append(groups.index(transition_group(source, destination)))
            ax = axes[row_index, column_index]
            ax.imshow(
                np.asarray(grid).reshape(len(kappas), len(lambdas)),
                origin="lower", aspect="auto", cmap=cmap,
                vmin=-0.5, vmax=len(groups) - 0.5,
                extent=(*plot_extent(lambdas), *plot_extent(kappas)),
            )
            if row_index == 0:
                ax.set_title(f"n = {hill:g}")
            ax.set(xlabel="baseline dimensionless length Λ",
                   ylabel=f"1× → {volume:g}× volume\nbaseline κ_R")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=color,
                          label=group.replace("_", " "))
               for group, color in TRANSITION_COLORS.items()]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Grouped polarity-regime transitions")
    path = output / "grouped_transition_maps.png"
    fig.savefig(path, dpi=210)
    plt.close(fig)
    return path


def summarize_transitions(rows, hills, destinations) -> dict[str, object]:
    lookup = paired_lookup(rows)
    summary: dict[str, object] = {}
    for volume in destinations:
        by_hill = {}
        for hill in hills:
            keys = sorted(key[2:] for key in lookup if np.isclose(key[0], 1.0) and key[1] == hill)
            exact: dict[str, int] = {}
            grouped: dict[str, int] = {}
            for kappa, length in keys:
                source = str(lookup[(1.0, hill, kappa, length)]["terminal_polarity_status"])
                destination = str(lookup[(volume, hill, kappa, length)]["terminal_polarity_status"])
                label = f"{source} -> {destination}"
                exact[label] = exact.get(label, 0) + 1
                group = transition_group(source, destination)
                grouped[group] = grouped.get(group, 0) + 1
            by_hill[f"n_{hill:g}"] = {"exact": exact, "grouped": grouped}
        summary[f"1x_to_{volume:g}x"] = by_hill
    return summary


def export_tier1_baseline(root: Path, rows, lambdas, kappas, hills) -> Path:
    """Materialize the 1x slice as the corresponding calibration analysis."""
    output = root.parents[1] / "calibration" / "tier1_regime_map"
    output.mkdir(parents=True, exist_ok=True)
    baseline = []
    for source in rows:
        if not np.isclose(float(source["volume_multiplier"]), 1.0):
            continue
        row = dict(source)
        row["Lambda"] = row["baseline_Lambda"]
        row["physical_length_um"] = row["baseline_physical_length_um"]
        baseline.append(row)
    baseline.sort(key=lambda row: (float(row["n"]), float(row["kappa_R"]), float(row["Lambda"])))
    with (output / "tier1_regime_map.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(baseline[0]))
        writer.writeheader()
        writer.writerows(baseline)
    counts: dict[str, int] = {}
    for row in baseline:
        outcome = str(row["terminal_polarity_status"])
        counts[outcome] = counts.get(outcome, 0) + 1
    summary = {
        "analysis": "three_species_tier1_regime_map",
        "model": "Holmes six-field GTPase-only reduction",
        "source": "1x slice of paired WGD-null array",
        "axes": {"Lambda": lambdas.tolist(), "kappa_R": kappas.tolist(), "n": hills},
        "fixed": {
            "D_m_um2_per_s": 0.1, "D_i_um2_per_s": 50.0,
            "diffusion_ratio": 500.0, "alpha_per_s": 0.65,
        },
        "state_counts": counts,
        "conditions": len(baseline),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    figures = render_maps(output, baseline, lambdas, kappas, hills)
    write_figure_manifest(
        output,
        title="Three-species Holmes Tier-1 parameter screen",
        introduction=(
            "The 1x-volume slice of the paired array maps the transient Rac-cue "
            "response over physical size, basal Rac activation, and cooperativity."
        ),
        figures=figures,
        data_files=("summary.json", "tier1_regime_map.csv", "figure_manifest.json"),
    )
    return output


def write_report(output: Path, summary: dict[str, object], figures: list[Path]) -> Path:
    blocks = []
    for comparison, by_hill in summary["transitions"].items():
        rows = []
        for hill, records in by_hill.items():
            for transition, count in sorted(records["exact"].items()):
                source, destination = transition.split(" -> ")
                if source != destination:
                    rows.append(
                        f"<tr><td>{escape(hill.replace('_', '='))}</td>"
                        f"<td>{escape(transition.replace('_', ' '))}</td><td>{count}</td></tr>"
                    )
        blocks.append(
            f"<section><h3>{escape(comparison.replace('_', ' '))}</h3>"
            "<div class='table-wrap'><table><thead><tr><th>Cooperativity</th>"
            f"<th>Exact changed transition</th><th>Conditions</th></tr></thead><tbody>"
            f"{''.join(rows) or '<tr><td colspan=3>No changed outcomes</td></tr>'}"
            "</tbody></table></div></section>"
        )
    figure_html = "".join(
        f"<figure><a href='{path.name}'><img src='{path.name}' alt='{escape(path.stem)}' "
        f"loading='lazy'></a><figcaption><strong>{escape(path.stem.replace('_', ' ').title())}</strong></figcaption></figure>"
        for path in figures
    )
    extended_path = output / "extended_report_fragment.html"
    extended_html = extended_path.read_text() if extended_path.exists() else ""
    report = output / "report.html"
    report.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width"><title>Three-species WGD-null report</title>
<style>{STYLE}</style></head><body><main><h1>Three-species Holmes WGD-null report</h1>
<p class="lead">A paired Cdc42–Rac–Rho screen showing exactly where isometric
cell-volume scaling changes the response to a transient local Rac cue.</p>
<section class="guide"><h2>Paired geometric null</h2><p>Each baseline
<code>(Λ, κ<sub>R</sub>, n)</code> point is evaluated at 0.5×, 1×, and 2× volume.
For volume multiplier <code>V</code>, physical length and <code>Λ</code> scale by
<code>V<sup>1/3</sup></code>. Cdc42, Rac, and Rho concentrations, molecular
kinetics, diffusivities, cue amplitude and timing, and the cue's fractional
position and width remain fixed. The native six-field reduction has no explicit
<code>h=V/A<sub>mem</sub></code>; that requires the later nine-field exchange
model.</p><p>All map axes are baseline 1× coordinates, so one pixel is one exact
paired molecular condition. The primary figures retain exact source→destination
labels and exact count matrices; grouped maps are a secondary visual summary.</p></section>
<section class="analysis"><h2>Exact changed transitions</h2>{''.join(blocks)}</section>
<div class="figures">{figure_html}</div>
{extended_html}
<section class="guide"><h2>Scope</h2><p>This is a concentration-preserving
geometric null, not a fitted model of a specific WGD cell line. Terminal outcomes
are evaluated at 300 s. Active Rac uses the same normalized-amplitude and domain
classification as the single-species analysis. Composite Cdc42–Rac–Rho coherence
and linear stability are reported separately as mechanistic annotations.</p></section>
<h2>Downloads</h2><ul class="files"><li><a href="wgd_null_screen.csv">wgd_null_screen.csv</a></li>
<li><a href="summary.json">summary.json</a></li></ul></main></body></html>""")
    return report


def merge_results(args) -> Path:
    root = output_root(args)
    tasks = all_tasks(args)
    rows: list[dict[str, object]] = []
    for index, task in enumerate(tasks):
        path = root / "tasks" / f"task_{index:05d}.json"
        if not path.exists():
            raise FileNotFoundError(f"missing single-cue task result {index}: {path}")
        payload = json.loads(path.read_text())
        if payload.get("task_index") != index or tuple(payload.get("task_key", ())) != task_key(task):
            raise ValueError(f"single-cue task {index} does not match its assignment")
        row = payload.get("row")
        if not isinstance(row, dict) or row_key(row) != task_key(task):
            raise ValueError(f"single-cue task {index} has an invalid result row")
        rows.append(row)
    keys = {row_key(row) for row in rows}
    if len(rows) != len(tasks) or len(keys) != len(tasks):
        raise ValueError(f"expected {len(tasks)} unique rows, found {len(rows)} rows and {len(keys)} keys")
    invalid_outcomes = {str(row["terminal_polarity_status"]) for row in rows} - set(OUTCOMES)
    if invalid_outcomes:
        raise ValueError(f"invalid terminal outcomes: {sorted(invalid_outcomes)}")
    maximum_error = max(
        max(float(row[f"maximum_{species}_conservation_error"]) for species in ("cdc42", "rac", "rho"))
        for row in rows
    )
    if maximum_error > args.maximum_conservation_error:
        raise ValueError(
            f"maximum conservation error {maximum_error:g} exceeds "
            f"{args.maximum_conservation_error:g}"
        )
    rows.sort(key=row_key)
    root.mkdir(parents=True, exist_ok=True)
    with (root / "wgd_null_screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lambdas = np.linspace(args.lambda_min, args.lambda_max, args.lambda_count)
    kappas = np.linspace(args.kappa_min, args.kappa_max, args.kappa_count)
    hills = parse_values(args.hill_values)
    volumes = parse_values(args.volume_multipliers)
    if not any(np.isclose(volume, 1.0) for volume in volumes):
        raise ValueError("volume multipliers must include the 1x baseline")
    destinations = [volume for volume in volumes if not np.isclose(volume, 1.0)]
    figures = [
        render_exact_transition_maps(root, rows, lambdas, kappas, hills, destinations),
        render_transition_matrices(root, rows, hills, destinations),
        render_grouped_transition_maps(root, rows, lambdas, kappas, hills, destinations),
    ]
    tier1_output = export_tier1_baseline(root, rows, lambdas, kappas, hills)
    summary = {
        "analysis": "three_species_holmes_wgd_null_screen",
        "volume_multipliers": volumes,
        "linear_scales": {f"{volume:g}x": volume ** (1.0 / 3.0) for volume in volumes},
        "baseline_grid": {
            "Lambda": lambdas.tolist(), "kappa_R": kappas.tolist(), "n": hills,
            "conditions": len(lambdas) * len(kappas) * len(hills),
        },
        "total_simulations": len(rows),
        "tasks_validated": len(tasks),
        "null": {
            "fixed": ["Cdc42/Rac/Rho concentrations", "kinetics", "diffusivities",
                      "stimulus amplitude", "stimulus timing",
                      "stimulus fractional position and width"],
            "scaled": ["physical length L", "dimensionless length Lambda"],
            "linear_scale": "volume_multiplier^(1/3)",
        },
        "transitions": summarize_transitions(rows, hills, destinations),
        "maximum_conservation_error": maximum_error,
        "tier1_calibration_export": str(tier1_output),
        "figures": [path.name for path in figures],
    }
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return write_report(root, summary, figures)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--mode", choices=("task", "shard", "merge"), required=True)
    result.add_argument("--output", type=Path, default=Path("migration1d/results/three_species/testing/wgd_null"))
    result.add_argument("--volume-multipliers", default="0.5,1,2")
    result.add_argument("--shard-count", type=int, default=32)
    result.add_argument("--shard-index", type=int, default=0)
    result.add_argument("--task-index", type=int, default=0)
    result.add_argument("--workers", type=int, default=1)
    result.add_argument("--lambda-min", type=float, default=31.6227766017)
    result.add_argument("--lambda-max", type=float, default=94.8683298051)
    result.add_argument("--lambda-count", type=int, default=13)
    result.add_argument("--kappa-min", type=float, default=0.0)
    result.add_argument("--kappa-max", type=float, default=0.08)
    result.add_argument("--kappa-count", type=int, default=13)
    result.add_argument("--hill-values", default="2,3,4")
    result.add_argument("--duration", type=float, default=300.0)
    result.add_argument("--nodes", type=int, default=81)
    result.add_argument("--output-interval", type=float, default=2.0)
    result.add_argument("--max-step", type=float, default=0.5)
    result.add_argument("--polarity-threshold", type=float, default=0.2)
    result.add_argument("--stimulus-amplitude", type=float, default=1.0)
    result.add_argument("--stimulus-center", type=float, default=0.1)
    result.add_argument("--stimulus-width", type=float, default=0.05)
    result.add_argument("--stimulus-start", type=float, default=10.0)
    result.add_argument("--stimulus-stop", type=float, default=30.0)
    result.add_argument("--maximum-conservation-error", type=float, default=1.0e-6)
    return result


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    path = (
        run_array_task(args) if args.mode == "task"
        else run_shard(args) if args.mode == "shard"
        else merge_results(args)
    )
    print(path)


if __name__ == "__main__":
    main()
