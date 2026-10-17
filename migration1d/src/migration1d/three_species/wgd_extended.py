"""Full Holmes WGD-null suite: two seeded Rac domains and polarity reversal."""

from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from html import escape
import json
import os
from pathlib import Path

import numpy as np

from migration1d.holmes import HolmesChemistry
from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from migration1d.single_species.observables import polarized_domain_count
from migration1d.single_species.testing import endpoint_contrast, first_sustained_time
from migration1d.single_species.wgd_extended import (
    DOMAIN_TRANSITION_COLORS,
    DOMAIN_TRANSITION_LABELS,
    REVERSAL_COLORS,
    REVERSAL_OUTCOMES,
    REVERSAL_TRANSITION_COLORS,
    REVERSAL_TRANSITION_LABELS,
    domain_transition,
    reversal_transition,
)
from .calibration import plot_extent
from .model import (
    GaussianRacPulse,
    RacPulseTrain,
    canonical_holmes_parameters,
    homogeneous_equilibrium,
)
from .observables import profile_metrics
from .wgd import parse_values, write_report


class SeededRacHolmesChemistry(HolmesChemistry):
    """Holmes chemistry with two mass-conserving active-Rac enrichments."""

    def __init__(self, parameters, amplitude: float, width: float, asymmetry: float):
        super().__init__(parameters)
        self.amplitude = amplitude
        self.width = width
        self.asymmetry = asymmetry

    def initial_profiles(self, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        active = np.repeat(self.initial_membrane[:, None], len(position), axis=1)
        inactive = np.repeat(self.initial_cytoplasm[:, None], len(position), axis=1)
        fraction = (position - position[0]) / (position[-1] - position[0])
        left = np.exp(-0.5 * ((fraction - 0.2) / self.width) ** 2)
        right = np.exp(-0.5 * ((fraction - 0.8) / self.width) ** 2)
        transfer = self.amplitude * (
            (1.0 + self.asymmetry) * left + (1.0 - self.asymmetry) * right
        )
        transfer = np.minimum(transfer, 0.9 * inactive[1])
        active[1] += transfer
        inactive[1] -= transfer
        return active, inactive


def output_root(args) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else args.output.resolve()


def report_root(args) -> Path:
    if args.report_root is None:
        return output_root(args)
    project = Path(__file__).resolve().parents[3]
    return (project / args.report_root).resolve()


def task_key(task: tuple) -> tuple:
    return tuple(float(value) if isinstance(value, (float, np.floating)) else value for value in task)


def parameters(kappa: float, hill: float, baseline_lambda: float, volume: float):
    p = canonical_holmes_parameters(kappa_r=kappa, hill_exponent=hill)
    rest = homogeneous_equilibrium(p)
    p = replace(p, initial_active=tuple(float(value) for value in rest))
    scale = volume ** (1.0 / 3.0)
    scaled_lambda = baseline_lambda * scale
    length = scaled_lambda * np.sqrt(p.active_diffusivity)
    return p, rest, scale, scaled_lambda, length


def simulate(args, chemistry, length: float):
    return run_simulation(
        MigrationModel(MechanicsParameters(reference_length=length), chemistry),
        SimulationConfig(
            nodes=args.nodes, duration=args.duration,
            output_interval=args.output_interval, max_step=args.max_step,
            noise_amplitude=0.0, mechanics_activation_time=None,
            jacobian_sparsity=True,
        ),
    )


def conservation_errors(result) -> dict[str, float]:
    active, inactive = result.compartment_totals()
    total = active + inactive
    errors = np.max(
        np.abs(total - total[0]) / np.maximum(np.abs(total[0]), 1.0e-14), axis=0
    )
    return {
        f"maximum_{name}_conservation_error": float(value)
        for name, value in zip(("cdc42", "rac", "rho"), errors)
    }


def common_row(kind, volume, hill, kappa, baseline_lambda, rest, scale,
               scaled_lambda, length):
    return {
        "experiment": kind,
        "volume_multiplier": volume,
        "linear_scale": scale,
        "n": hill,
        "kappa_R": kappa,
        "baseline_Lambda": baseline_lambda,
        "scaled_Lambda": scaled_lambda,
        "scaled_length_um": length,
        "resting_cdc42": float(rest[0]),
        "resting_rac": float(rest[1]),
        "resting_rho": float(rest[2]),
    }


def run_dual(args, task) -> dict[str, object]:
    _, volume, hill, epsilon, kappa, baseline_lambda = task
    p, rest, scale, scaled_lambda, length = parameters(kappa, hill, baseline_lambda, volume)
    chemistry = SeededRacHolmesChemistry(
        p, args.initial_amplitude, args.initial_width, epsilon
    )
    result = simulate(args, chemistry, length)
    profiles = result.membrane[:, 1]
    counts = {}
    for selected in parse_values(args.domain_count_times):
        index = int(np.argmin(np.abs(result.times - selected)))
        counts[f"polarized_domain_count_t{selected:g}"] = polarized_domain_count(profiles[index])
    metrics = profile_metrics(result.membrane[-1])
    return {
        **common_row("two_seeded_domains", volume, hill, kappa, baseline_lambda,
                     rest, scale, scaled_lambda, length),
        "epsilon": epsilon,
        **counts,
        "terminal_coherent_polarity": bool(metrics["coherent_polarity"]),
        "terminal_cdc42_rac_correlation": float(metrics["cdc42_rac_correlation"]),
        "terminal_cdc42_rho_correlation": float(metrics["cdc42_rho_correlation"]),
        **conservation_errors(result),
    }


def run_reversal(args, task) -> dict[str, object]:
    _, volume, hill, kappa, baseline_lambda = task
    p, rest, scale, scaled_lambda, length = parameters(kappa, hill, baseline_lambda, volume)
    cue = RacPulseTrain((
        GaussianRacPulse(args.stimulus_amplitude, 0.1, args.stimulus_width, 10.0, 30.0),
        GaussianRacPulse(args.reversal_amplitude, 0.9, args.stimulus_width,
                         args.reversal_start, args.reversal_start + args.reversal_duration),
    ))
    result = simulate(args, HolmesChemistry(p, cue=cue), length)
    profiles = result.membrane[:, 1]
    contrast = endpoint_contrast(profiles)
    before = int(np.flatnonzero(result.times < args.reversal_start)[-1])
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
        result.times, (counts == 1) & (contrast >= args.reversal_threshold),
        start=args.reversal_start, duration=10.0,
    ) if initially_left else None
    metrics = profile_metrics(result.membrane[-1])
    return {
        **common_row("polarity_reversal", volume, hill, kappa, baseline_lambda,
                     rest, scale, scaled_lambda, length),
        "prechallenge_time": float(result.times[before]),
        "prechallenge_polarized_domain_count": before_count,
        "prechallenge_endpoint_contrast": float(contrast[before]),
        "terminal_polarized_domain_count": final_count,
        "terminal_endpoint_contrast": float(contrast[-1]),
        "reversal_outcome": outcome,
        "reversal_time_after_challenge": reversal_time,
        "terminal_coherent_polarity": bool(metrics["coherent_polarity"]),
        **conservation_errors(result),
    }


def all_tasks(args) -> list[tuple]:
    volumes = parse_values(args.volume_multipliers)
    hills = parse_values(args.hill_values)
    dual_kappas = np.linspace(args.dual_kappa_min, args.dual_kappa_max, args.dual_kappa_count)
    dual_lambdas = np.linspace(args.dual_lambda_min, args.dual_lambda_max, args.dual_lambda_count)
    reversal_kappas = np.linspace(
        args.reversal_kappa_min, args.reversal_kappa_max, args.reversal_kappa_count
    )
    reversal_lambdas = np.linspace(
        args.reversal_lambda_min, args.reversal_lambda_max, args.reversal_lambda_count
    )
    dual = [
        ("dual", float(volume), float(hill), float(epsilon), float(kappa), float(length))
        for volume in volumes for hill in hills for epsilon in parse_values(args.asymmetries)
        for kappa in dual_kappas for length in dual_lambdas
    ]
    reversal = [
        ("reversal", float(volume), float(hill), float(kappa), float(length))
        for volume in volumes for hill in hills
        for kappa in reversal_kappas for length in reversal_lambdas
    ]
    return dual + reversal


def run_task(args) -> Path:
    tasks = all_tasks(args)
    index = int(os.environ.get("SLURM_ARRAY_TASK_ID", args.task_index))
    if not 0 <= index < len(tasks):
        raise ValueError(f"task index {index} outside 0..{len(tasks) - 1}")
    task = tasks[index]
    destination = output_root(args) / "extended_tasks" / f"task_{index:05d}.json"
    if destination.exists():
        payload = json.loads(destination.read_text())
        if payload.get("task_index") == index and tuple(payload.get("task_key", ())) == task_key(task):
            print(destination)
            return destination
    row = run_dual(args, task) if task[0] == "dual" else run_reversal(args, task)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps({
        "task_index": index, "task_key": task_key(task), "row": row,
    }, separators=(",", ":")) + "\n")
    print(destination)
    return destination


def paired(rows, extra=()):
    return {
        (float(row["volume_multiplier"]), float(row["n"]),
         *(float(row[item]) for item in extra), float(row["kappa_R"]),
         float(row["baseline_Lambda"])): row
        for row in rows
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def render_dual_maps(output: Path, rows, args) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    kappas = np.linspace(args.dual_kappa_min, args.dual_kappa_max, args.dual_kappa_count)
    lambdas = np.linspace(args.dual_lambda_min, args.dual_lambda_max, args.dual_lambda_count)
    hills = parse_values(args.hill_values)
    destinations = [v for v in parse_values(args.volume_multipliers) if not np.isclose(v, 1)]
    lookup = paired(rows, ("epsilon",))
    classes = list(DOMAIN_TRANSITION_COLORS)
    cmap = ListedColormap([DOMAIN_TRANSITION_COLORS[item] for item in classes])
    norm = BoundaryNorm(np.arange(-0.5, len(classes) + 0.5), cmap.N)
    paths = []
    for epsilon in parse_values(args.asymmetries):
        fig, axes = plt.subplots(2, len(hills), figsize=(5 * len(hills), 8), squeeze=False,
                                 constrained_layout=True)
        for i, volume in enumerate(destinations):
            for j, hill in enumerate(hills):
                grid = []
                for kappa in kappas:
                    for length in lambdas:
                        source = int(lookup[(1.0, hill, epsilon, float(kappa), float(length))]["polarized_domain_count_t300"])
                        target = int(lookup[(volume, hill, epsilon, float(kappa), float(length))]["polarized_domain_count_t300"])
                        grid.append(classes.index(domain_transition(source, target)))
                axes[i, j].imshow(np.asarray(grid).reshape(len(kappas), len(lambdas)),
                                  origin="lower", aspect="auto", cmap=cmap, norm=norm,
                                  extent=(*plot_extent(lambdas), *plot_extent(kappas)))
                axes[i, j].set(title=f"n = {hill:g}", xlabel="baseline Λ",
                               ylabel=f"1× → {volume:g}×\nbaseline κ_R")
        handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=DOMAIN_TRANSITION_COLORS[c],
                              label=DOMAIN_TRANSITION_LABELS[c]) for c in classes]
        fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
        fig.suptitle(f"Active-Rac polarized-region transition at 300 s; ε={epsilon:g}")
        path = output / f"double_stimulus_transitions_epsilon_{epsilon:g}.png"
        fig.savefig(path, dpi=210); plt.close(fig); paths.append(path)
    return paths


def render_reversal_maps(output: Path, rows, args) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    kappas = np.linspace(args.reversal_kappa_min, args.reversal_kappa_max, args.reversal_kappa_count)
    lambdas = np.linspace(args.reversal_lambda_min, args.reversal_lambda_max, args.reversal_lambda_count)
    hills = parse_values(args.hill_values)
    volumes = parse_values(args.volume_multipliers)
    lookup = paired(rows)
    cmap = ListedColormap([REVERSAL_COLORS[o] for o in REVERSAL_OUTCOMES])
    norm = BoundaryNorm(np.arange(-0.5, len(REVERSAL_OUTCOMES) + 0.5), cmap.N)
    fig, axes = plt.subplots(len(volumes), len(hills), figsize=(5 * len(hills), 12),
                             squeeze=False, constrained_layout=True)
    for i, volume in enumerate(volumes):
        for j, hill in enumerate(hills):
            grid = [REVERSAL_OUTCOMES.index(str(lookup[(volume, hill, float(k), float(l))]["reversal_outcome"]))
                    for k in kappas for l in lambdas]
            axes[i, j].imshow(np.asarray(grid).reshape(len(kappas), len(lambdas)), origin="lower",
                              aspect="auto", cmap=cmap, norm=norm,
                              extent=(*plot_extent(lambdas), *plot_extent(kappas)))
            axes[i, j].set(title=f"n = {hill:g}", xlabel="baseline Λ",
                           ylabel=f"{volume:g}× volume\nbaseline κ_R")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=REVERSAL_COLORS[o],
                          label=o.replace("_", " ")) for o in REVERSAL_OUTCOMES]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Active-Rac polarity-reversal outcome under the WGD null")
    outcome_path = output / "reversal_outcomes_by_volume.png"
    fig.savefig(outcome_path, dpi=210); plt.close(fig)
    classes = list(REVERSAL_TRANSITION_COLORS)
    cmap = ListedColormap([REVERSAL_TRANSITION_COLORS[c] for c in classes])
    norm = BoundaryNorm(np.arange(-0.5, len(classes) + 0.5), cmap.N)
    destinations = [v for v in volumes if not np.isclose(v, 1)]
    fig, axes = plt.subplots(2, len(hills), figsize=(5 * len(hills), 8), squeeze=False,
                             constrained_layout=True)
    for i, volume in enumerate(destinations):
        for j, hill in enumerate(hills):
            grid = []
            for kappa in kappas:
                for length in lambdas:
                    source = str(lookup[(1.0, hill, float(kappa), float(length))]["reversal_outcome"])
                    target = str(lookup[(volume, hill, float(kappa), float(length))]["reversal_outcome"])
                    grid.append(classes.index(reversal_transition(source, target)))
            axes[i, j].imshow(np.asarray(grid).reshape(len(kappas), len(lambdas)), origin="lower",
                              aspect="auto", cmap=cmap, norm=norm,
                              extent=(*plot_extent(lambdas), *plot_extent(kappas)))
            axes[i, j].set(title=f"n = {hill:g}", xlabel="baseline Λ",
                           ylabel=f"1× → {volume:g}×\nbaseline κ_R")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=REVERSAL_TRANSITION_COLORS[c],
                          label=REVERSAL_TRANSITION_LABELS[c]) for c in classes]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Where volume scaling changes reversal behavior")
    transition_path = output / "reversal_transition_maps.png"
    fig.savefig(transition_path, dpi=210); plt.close(fig)
    return [outcome_path, transition_path]


def summarize_dual(rows, args):
    lookup = paired(rows, ("epsilon",))
    keys = sorted({key[1:] for key in lookup if key[0] == 1.0})
    summary = {"analysis": "three_species_wgd_null_two_seeded_domains",
               "simulations": len(rows), "primary_field": "active Rac",
               "domain_count_times_seconds": parse_values(args.domain_count_times),
               "transitions": {}, "time_sequences": {}}
    for volume in parse_values(args.volume_multipliers):
        for epsilon in parse_values(args.asymmetries):
            sequences = {}
            for hill, eps, kappa, length in keys:
                if not np.isclose(eps, epsilon): continue
                row = lookup[(volume, hill, eps, kappa, length)]
                sequence = "→".join(str(row[f"polarized_domain_count_t{t:g}"])
                                    for t in parse_values(args.domain_count_times))
                sequences[sequence] = sequences.get(sequence, 0) + 1
            summary["time_sequences"][f"volume={volume:g},epsilon={epsilon:g}"] = sequences
    for volume in [v for v in parse_values(args.volume_multipliers) if not np.isclose(v, 1)]:
        for epsilon in parse_values(args.asymmetries):
            exact, grouped = {}, {}
            for hill, eps, kappa, length in keys:
                if not np.isclose(eps, epsilon): continue
                source = int(lookup[(1.0, hill, eps, kappa, length)]["polarized_domain_count_t300"])
                target = int(lookup[(volume, hill, eps, kappa, length)]["polarized_domain_count_t300"])
                exact[f"{source}→{target}"] = exact.get(f"{source}→{target}", 0) + 1
                group = domain_transition(source, target)
                grouped[group] = grouped.get(group, 0) + 1
            summary["transitions"][f"1x_to_{volume:g}x,epsilon={epsilon:g}"] = {
                "exact": exact, "grouped": grouped,
            }
    return summary


def summarize_reversal(rows, args):
    lookup = paired(rows)
    keys = sorted({key[1:] for key in lookup if key[0] == 1.0})
    counts, transitions = {}, {}
    for volume in parse_values(args.volume_multipliers):
        record = {}
        for hill, kappa, length in keys:
            outcome = str(lookup[(volume, hill, kappa, length)]["reversal_outcome"])
            record[outcome] = record.get(outcome, 0) + 1
        counts[f"volume={volume:g}"] = record
    for volume in [v for v in parse_values(args.volume_multipliers) if not np.isclose(v, 1)]:
        exact, grouped = {}, {}
        for hill, kappa, length in keys:
            source = str(lookup[(1.0, hill, kappa, length)]["reversal_outcome"])
            target = str(lookup[(volume, hill, kappa, length)]["reversal_outcome"])
            exact[f"{source} -> {target}"] = exact.get(f"{source} -> {target}", 0) + 1
            group = reversal_transition(source, target)
            grouped[group] = grouped.get(group, 0) + 1
        transitions[f"1x_to_{volume:g}x"] = {"exact": exact, "grouped": grouped}
    return {"analysis": "three_species_wgd_null_polarity_reversal",
            "simulations": len(rows), "primary_field": "active Rac",
            "outcome_counts": counts, "transitions": transitions}


def report_fragment(dual, reversal, dual_figures, reversal_figures) -> str:
    dual_rows = "".join(
        f"<tr><td>{escape(key)}</td><td>{sum(value['grouped'].values()) - value['grouped'].get('unchanged', 0)}</td>"
        f"<td>{value['grouped'].get('one_to_two', 0)}</td><td>{value['grouped'].get('two_to_one', 0)}</td></tr>"
        for key, value in dual["transitions"].items()
    )
    reversal_rows = "".join(
        f"<tr><td>{escape(key)}</td><td>{sum(value['grouped'].values()) - value['grouped'].get('unchanged', 0)}</td>"
        f"<td>{value['grouped'].get('became_reversible', 0)}</td><td>{value['grouped'].get('lost_reversibility', 0)}</td></tr>"
        for key, value in reversal["transitions"].items()
    )
    dual_html = "".join(f"<figure><a href='double_stimulus/{p.name}'><img src='double_stimulus/{p.name}' loading='lazy'></a></figure>" for p in dual_figures)
    reversal_html = "".join(f"<figure><a href='polarity_reversal/{p.name}'><img src='polarity_reversal/{p.name}' loading='lazy'></a></figure>" for p in reversal_figures)
    return f"""
<section class='analysis'><h2>Two seeded competing Rac domains</h2><p>Two
mass-conserving active-Rac enrichments are seeded at 20% and 80% at t=0. Counts
at 0, 100, and 300 s use exactly the single-species normalized-amplitude and
connected-domain definition. Cdc42/Rac/Rho coherence is retained as a secondary
annotation.</p><div class='table-wrap'><table><thead><tr><th>Comparison</th><th>Changed</th><th>1 → 2</th><th>2 → 1</th></tr></thead><tbody>{dual_rows}</tbody></table></div><div class='figures'>{dual_html}</div></section>
<section class='analysis'><h2>Polarity reversal</h2><p>An active-Rac cue acts
on the left from 10–30 s and an equal cue on the right from 150–170 s. Eligibility,
endpoint contrast, and sustained reversal use the same definitions as the
single-species screen.</p><div class='table-wrap'><table><thead><tr><th>Comparison</th><th>Changed</th><th>Became reversible</th><th>Lost reversibility</th></tr></thead><tbody>{reversal_rows}</tbody></table></div><div class='figures'>{reversal_html}</div></section>"""


def merge(args) -> Path:
    work, root, tasks = output_root(args), report_root(args), all_tasks(args)
    rows = []
    for index, task in enumerate(tasks):
        path = work / "extended_tasks" / f"task_{index:05d}.json"
        if not path.exists():
            raise FileNotFoundError(f"missing extended task {index}: {path}")
        payload = json.loads(path.read_text())
        if payload.get("task_index") != index or tuple(payload.get("task_key", ())) != task_key(task):
            raise ValueError(f"extended task {index} does not match its assignment")
        rows.append(payload["row"])
    dual_rows = [row for row in rows if row.get("experiment") == "two_seeded_domains"]
    reversal_rows = [row for row in rows if row.get("experiment") == "polarity_reversal"]
    expected_dual = (len(parse_values(args.volume_multipliers)) * len(parse_values(args.hill_values))
                     * len(parse_values(args.asymmetries)) * args.dual_kappa_count * args.dual_lambda_count)
    expected_reversal = (len(parse_values(args.volume_multipliers)) * len(parse_values(args.hill_values))
                         * args.reversal_kappa_count * args.reversal_lambda_count)
    if len(dual_rows) != expected_dual or len(reversal_rows) != expected_reversal:
        raise ValueError(f"unexpected row counts: dual={len(dual_rows)}, reversal={len(reversal_rows)}")
    maximum_error = max(float(row[f"maximum_{species}_conservation_error"])
                        for row in rows for species in ("cdc42", "rac", "rho"))
    if maximum_error > args.maximum_conservation_error:
        raise ValueError(f"maximum conservation error {maximum_error:g} exceeds limit")
    dual_output, reversal_output = root / "double_stimulus", root / "polarity_reversal"
    write_csv(dual_output / "two_seeded_domains.csv", dual_rows)
    write_csv(reversal_output / "polarity_reversal.csv", reversal_rows)
    dual_summary, reversal_summary = summarize_dual(dual_rows, args), summarize_reversal(reversal_rows, args)
    dual_summary["maximum_conservation_error"] = maximum_error
    reversal_summary["maximum_conservation_error"] = maximum_error
    dual_figures = render_dual_maps(dual_output, dual_rows, args)
    reversal_figures = render_reversal_maps(reversal_output, reversal_rows, args)
    (dual_output / "summary.json").write_text(json.dumps(dual_summary, indent=2) + "\n")
    (reversal_output / "summary.json").write_text(json.dumps(reversal_summary, indent=2) + "\n")
    (root / "extended_report_fragment.html").write_text(
        report_fragment(dual_summary, reversal_summary, dual_figures, reversal_figures)
    )
    base_summary = json.loads((root / "summary.json").read_text())
    return write_report(root, base_summary, [root / name for name in base_summary["figures"]])


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--mode", choices=("task", "merge"), required=True)
    result.add_argument("--task-index", type=int, default=0)
    result.add_argument("--output", type=Path, default=Path("migration1d/results/three_species/testing/wgd_null_extended_work"))
    result.add_argument("--report-root", type=Path)
    result.add_argument("--volume-multipliers", default="0.5,1,2")
    result.add_argument("--hill-values", default="2,3,4")
    result.add_argument("--dual-kappa-min", type=float, default=0.0)
    result.add_argument("--dual-kappa-max", type=float, default=0.08)
    result.add_argument("--dual-kappa-count", type=int, default=13)
    result.add_argument("--dual-lambda-min", type=float, default=31.6227766017)
    result.add_argument("--dual-lambda-max", type=float, default=94.8683298051)
    result.add_argument("--dual-lambda-count", type=int, default=13)
    result.add_argument("--reversal-kappa-min", type=float, default=0.0)
    result.add_argument("--reversal-kappa-max", type=float, default=0.08)
    result.add_argument("--reversal-kappa-count", type=int, default=17)
    result.add_argument("--reversal-lambda-min", type=float, default=31.6227766017)
    result.add_argument("--reversal-lambda-max", type=float, default=94.8683298051)
    result.add_argument("--reversal-lambda-count", type=int, default=17)
    result.add_argument("--asymmetries", default="0,0.001,0.01,0.05")
    result.add_argument("--domain-count-times", default="0,100,300")
    result.add_argument("--initial-amplitude", type=float, default=2.0)
    result.add_argument("--initial-width", type=float, default=0.04)
    result.add_argument("--duration", type=float, default=300.0)
    result.add_argument("--nodes", type=int, default=81)
    result.add_argument("--output-interval", type=float, default=2.0)
    result.add_argument("--max-step", type=float, default=0.5)
    result.add_argument("--stimulus-amplitude", type=float, default=1.0)
    result.add_argument("--stimulus-width", type=float, default=0.05)
    result.add_argument("--reversal-start", type=float, default=150.0)
    result.add_argument("--reversal-duration", type=float, default=20.0)
    result.add_argument("--reversal-amplitude", type=float, default=1.0)
    result.add_argument("--reversal-threshold", type=float, default=0.5)
    result.add_argument("--maximum-conservation-error", type=float, default=1.0e-6)
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    path = run_task(args) if args.mode == "task" else merge(args)
    print(path)


if __name__ == "__main__":
    main()
