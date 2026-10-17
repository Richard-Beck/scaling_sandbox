"""WGD-null extensions for two seeded domains and polarity reversal."""

from __future__ import annotations

import argparse
import csv
from html import escape
import json
import os
from pathlib import Path

import numpy as np

from migration1d.model import MechanicsParameters, MigrationModel, SimulationConfig, run_simulation
from .calibration import plot_extent
from .model import ActivationPulse, DualPulseInitialCondition, PulseTrain, StripWavePinningChemistry
from .observables import polarized_domain_count
from .testing import endpoint_contrast, first_sustained_time
from .wgd import scaled_parameters, values, write_report


DOMAIN_TRANSITION_COLORS = {
    "unchanged": "#d9d9d9",
    "one_to_two": "#009E73",
    "two_to_one": "#0072B2",
    "gained_regions": "#56B4E9",
    "lost_regions": "#D55E00",
    "other_change": "#CC79A7",
}
DOMAIN_TRANSITION_LABELS = {
    "unchanged": "unchanged count",
    "one_to_two": "1 → 2 regions",
    "two_to_one": "2 → 1 regions",
    "gained_regions": "gained region(s)",
    "lost_regions": "lost region(s)",
    "other_change": "other count change",
}
REVERSAL_OUTCOMES = (
    "not_initially_left_polarized", "reversed", "retained_original_pole",
    "lost_polarity", "multiple_polarized_domains", "ambiguous",
)
REVERSAL_COLORS = {
    "not_initially_left_polarized": "#9E9E9E",
    "reversed": "#009E73",
    "retained_original_pole": "#0072B2",
    "lost_polarity": "#F0E442",
    "multiple_polarized_domains": "#CC79A7",
    "ambiguous": "#D55E00",
}
REVERSAL_TRANSITION_COLORS = {
    "unchanged": "#d9d9d9",
    "became_reversible": "#009E73",
    "lost_reversibility": "#D55E00",
    "became_eligible": "#56B4E9",
    "became_ineligible": "#8C6BB1",
    "other_change": "#CC79A7",
}
REVERSAL_TRANSITION_LABELS = {
    "unchanged": "unchanged outcome",
    "became_reversible": "became reversible",
    "lost_reversibility": "lost reversibility",
    "became_eligible": "became eligible",
    "became_ineligible": "became ineligible",
    "other_change": "other changed outcome",
}


def output_root(args) -> Path:
    configured = os.environ.get("MIGRATION1D_OUTPUT_DIR")
    return Path(configured).resolve() if configured else args.output.resolve()


def report_root(args) -> Path:
    if args.report_root is None:
        return output_root(args)
    project = Path(__file__).resolve().parents[3]
    return (project / args.report_root).resolve()


def simulate(args, parameters, chemistry):
    return run_simulation(
        MigrationModel(
            mechanics=MechanicsParameters(reference_length=parameters.domain_length),
            chemistry=chemistry,
        ),
        SimulationConfig(
            nodes=args.nodes, duration=args.duration,
            output_interval=args.output_interval, max_step=args.max_step,
            noise_amplitude=0.0, mechanics_activation_time=None,
        ),
    )


def conservation_error(result, h: float) -> float:
    membrane, cytoplasm = result.compartment_totals()
    total = membrane[:, 0] + h * cytoplasm[:, 0]
    return float(np.max(np.abs(total - total[0])) / max(abs(total[0]), 1.0e-14))


def common_row(kind, volume, baseline_mu, baseline_lambda, baseline_h, parameters,
               resting_u, resting_v, linear_scale):
    return {
        "experiment": kind,
        "volume_multiplier": volume,
        "linear_scale": linear_scale,
        "baseline_mu": baseline_mu,
        "baseline_Lambda": baseline_lambda,
        "baseline_h": baseline_h,
        "resting_active_u": resting_u,
        "resting_inactive_v": resting_v,
        "scaled_mu": parameters.total_abundance,
        "scaled_Lambda": parameters.dimensionless_length,
        "scaled_h": parameters.cytoplasm_depth,
        "scaled_length": parameters.domain_length,
    }


def run_dual(args, task) -> dict[str, object]:
    _, volume, h, epsilon, mu, dimensionless_length = task
    p, u0, v0, scale = scaled_parameters(args, mu, dimensionless_length, h, volume)
    chemistry = StripWavePinningChemistry(
        p,
        initial_condition=DualPulseInitialCondition(
            amplitude=args.initial_amplitude,
            width_fraction=args.initial_width,
            asymmetry=epsilon,
        ),
    )
    result = simulate(args, p, chemistry)
    profiles = result.membrane[:, 0]
    counts = {}
    for selected in values(args.domain_count_times):
        index = int(np.argmin(np.abs(result.times - selected)))
        counts[f"polarized_domain_count_t{selected:g}"] = polarized_domain_count(profiles[index])
    return {
        **common_row("two_seeded_domains", volume, mu, dimensionless_length, h,
                     p, u0, v0, scale),
        "epsilon": epsilon,
        **counts,
        "maximum_conservation_error": conservation_error(result, p.cytoplasm_depth),
    }


def run_reversal(args, task) -> dict[str, object]:
    _, volume, h, mu, dimensionless_length = task
    p, u0, v0, scale = scaled_parameters(args, mu, dimensionless_length, h, volume)
    activation = PulseTrain((
        ActivationPulse(
            amplitude=args.stimulus_amplitude, center_fraction=0.1,
            width_fraction=args.stimulus_width, start=10.0, stop=30.0,
        ),
        ActivationPulse(
            amplitude=args.reversal_amplitude, center_fraction=0.9,
            width_fraction=args.stimulus_width, start=args.reversal_start,
            stop=args.reversal_start + args.reversal_duration,
        ),
    ))
    result = simulate(args, p, StripWavePinningChemistry(p, activation=activation))
    profiles = result.membrane[:, 0]
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
    return {
        **common_row("polarity_reversal", volume, mu, dimensionless_length, h,
                     p, u0, v0, scale),
        "prechallenge_time": float(result.times[before]),
        "prechallenge_polarized_domain_count": before_count,
        "prechallenge_endpoint_contrast": float(contrast[before]),
        "terminal_polarized_domain_count": final_count,
        "terminal_endpoint_contrast": float(contrast[-1]),
        "reversal_outcome": outcome,
        "reversal_time_after_challenge": reversal_time,
        "maximum_conservation_error": conservation_error(result, p.cytoplasm_depth),
    }


def all_tasks(args) -> list[tuple]:
    volumes = values(args.volume_multipliers)
    h_values = values(args.h_values)
    dual_mus = np.linspace(args.dual_mu_min, args.dual_mu_max, args.dual_mu_count)
    dual_lambdas = np.linspace(
        args.dual_lambda_min, args.dual_lambda_max, args.dual_lambda_count
    )
    reversal_mus = np.linspace(
        args.reversal_mu_min, args.reversal_mu_max, args.reversal_mu_count
    )
    reversal_lambdas = np.linspace(
        args.reversal_lambda_min, args.reversal_lambda_max, args.reversal_lambda_count
    )
    dual = [
        ("dual", float(volume), float(h), float(epsilon), float(mu), float(length))
        for volume in volumes for h in h_values for epsilon in values(args.asymmetries)
        for mu in dual_mus for length in dual_lambdas
    ]
    reversal = [
        ("reversal", float(volume), float(h), float(mu), float(length))
        for volume in volumes for h in h_values
        for mu in reversal_mus for length in reversal_lambdas
    ]
    return dual + reversal


def run_task(args) -> Path:
    tasks = all_tasks(args)
    index = int(os.environ.get("SLURM_ARRAY_TASK_ID", args.task_index))
    if not 0 <= index < len(tasks):
        raise ValueError(f"task index {index} outside 0..{len(tasks) - 1}")
    task = tasks[index]
    row = run_dual(args, task) if task[0] == "dual" else run_reversal(args, task)
    destination_dir = output_root(args) / "extended_tasks"
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / f"task_{index:05d}.json"
    destination.write_text(json.dumps(row, separators=(",", ":")) + "\n")
    print(destination)
    return destination


def paired(rows, extra=()):
    return {
        (float(row["volume_multiplier"]), float(row["baseline_h"]),
         *(float(row[item]) for item in extra), float(row["baseline_mu"]),
         float(row["baseline_Lambda"])): row
        for row in rows
    }


def domain_transition(source: int, destination: int) -> str:
    if source == destination:
        return "unchanged"
    if (source, destination) == (1, 2):
        return "one_to_two"
    if (source, destination) == (2, 1):
        return "two_to_one"
    if destination > source:
        return "gained_regions"
    if destination < source:
        return "lost_regions"
    return "other_change"


def reversal_transition(source: str, destination: str) -> str:
    if source == destination:
        return "unchanged"
    if source != "reversed" and destination == "reversed":
        return "became_reversible"
    if source == "reversed" and destination != "reversed":
        return "lost_reversibility"
    ineligible = "not_initially_left_polarized"
    if source == ineligible and destination != ineligible:
        return "became_eligible"
    if source != ineligible and destination == ineligible:
        return "became_ineligible"
    return "other_change"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def render_dual_transition_maps(output, rows, args) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    mus = np.linspace(args.dual_mu_min, args.dual_mu_max, args.dual_mu_count)
    lambdas = np.linspace(args.dual_lambda_min, args.dual_lambda_max, args.dual_lambda_count)
    h_values = values(args.h_values)
    destinations = [v for v in values(args.volume_multipliers) if not np.isclose(v, 1)]
    lookup = paired(rows, ("epsilon",))
    classes = list(DOMAIN_TRANSITION_COLORS)
    cmap = ListedColormap([DOMAIN_TRANSITION_COLORS[item] for item in classes])
    norm = BoundaryNorm(np.arange(-0.5, len(classes) + 0.5), cmap.N)
    paths = []
    for epsilon in values(args.asymmetries):
        fig, axes = plt.subplots(2, len(h_values), figsize=(5 * len(h_values), 8),
                                 squeeze=False, constrained_layout=True)
        for row_index, volume in enumerate(destinations):
            for column_index, h in enumerate(h_values):
                grid = []
                for mu in mus:
                    for length in lambdas:
                        source = int(lookup[(1.0, h, epsilon, float(mu), float(length))]["polarized_domain_count_t300"])
                        destination = int(lookup[(volume, h, epsilon, float(mu), float(length))]["polarized_domain_count_t300"])
                        grid.append(classes.index(domain_transition(source, destination)))
                ax = axes[row_index, column_index]
                ax.imshow(np.asarray(grid).reshape(len(mus), len(lambdas)), origin="lower",
                          aspect="auto", cmap=cmap, norm=norm,
                          extent=(*plot_extent(lambdas), *plot_extent(mus)))
                if row_index == 0:
                    ax.set_title(f"baseline h = {h:g}")
                ax.set(xlabel="baseline dimensionless length Λ",
                       ylabel=f"1× → {volume:g}×\nbaseline abundance μ")
        handles = [plt.Line2D([0], [0], marker="s", linestyle="",
                              color=DOMAIN_TRANSITION_COLORS[item],
                              label=DOMAIN_TRANSITION_LABELS[item]) for item in classes]
        fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
        fig.suptitle(f"Change in polarized-region count at 300 s; ε={epsilon:g}")
        path = output / f"double_stimulus_transitions_epsilon_{epsilon:g}.png"
        fig.savefig(path, dpi=210)
        plt.close(fig)
        paths.append(path)
    return paths


def render_reversal_maps(output, rows, args) -> list[Path]:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    mus = np.linspace(args.reversal_mu_min, args.reversal_mu_max, args.reversal_mu_count)
    lambdas = np.linspace(
        args.reversal_lambda_min, args.reversal_lambda_max, args.reversal_lambda_count
    )
    h_values = values(args.h_values)
    volumes = values(args.volume_multipliers)
    destinations = [v for v in volumes if not np.isclose(v, 1)]
    lookup = paired(rows)
    outcome_cmap = ListedColormap([REVERSAL_COLORS[item] for item in REVERSAL_OUTCOMES])
    outcome_norm = BoundaryNorm(np.arange(-0.5, len(REVERSAL_OUTCOMES) + 0.5), outcome_cmap.N)
    fig, axes = plt.subplots(len(volumes), len(h_values), figsize=(5 * len(h_values), 12),
                             squeeze=False, constrained_layout=True)
    for row_index, volume in enumerate(volumes):
        for column_index, h in enumerate(h_values):
            grid = [REVERSAL_OUTCOMES.index(str(lookup[(volume, h, float(mu), float(length))]["reversal_outcome"]))
                    for mu in mus for length in lambdas]
            ax = axes[row_index, column_index]
            ax.imshow(np.asarray(grid).reshape(len(mus), len(lambdas)), origin="lower",
                      aspect="auto", cmap=outcome_cmap, norm=outcome_norm,
                      extent=(*plot_extent(lambdas), *plot_extent(mus)))
            if row_index == 0:
                ax.set_title(f"baseline h = {h:g}")
            ax.set(xlabel="baseline dimensionless length Λ",
                   ylabel=f"{volume:g}× volume\nbaseline abundance μ")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="", color=REVERSAL_COLORS[item],
                          label=item.replace("_", " ")) for item in REVERSAL_OUTCOMES]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Polarity-reversal outcome under the WGD null")
    outcome_path = output / "reversal_outcomes_by_volume.png"
    fig.savefig(outcome_path, dpi=210)
    plt.close(fig)

    classes = list(REVERSAL_TRANSITION_COLORS)
    transition_cmap = ListedColormap([REVERSAL_TRANSITION_COLORS[item] for item in classes])
    transition_norm = BoundaryNorm(np.arange(-0.5, len(classes) + 0.5), transition_cmap.N)
    fig, axes = plt.subplots(2, len(h_values), figsize=(5 * len(h_values), 8),
                             squeeze=False, constrained_layout=True)
    for row_index, volume in enumerate(destinations):
        for column_index, h in enumerate(h_values):
            grid = []
            for mu in mus:
                for length in lambdas:
                    source = str(lookup[(1.0, h, float(mu), float(length))]["reversal_outcome"])
                    destination = str(lookup[(volume, h, float(mu), float(length))]["reversal_outcome"])
                    grid.append(classes.index(reversal_transition(source, destination)))
            ax = axes[row_index, column_index]
            ax.imshow(np.asarray(grid).reshape(len(mus), len(lambdas)), origin="lower",
                      aspect="auto", cmap=transition_cmap, norm=transition_norm,
                      extent=(*plot_extent(lambdas), *plot_extent(mus)))
            if row_index == 0:
                ax.set_title(f"baseline h = {h:g}")
            ax.set(xlabel="baseline dimensionless length Λ",
                   ylabel=f"1× → {volume:g}×\nbaseline abundance μ")
    handles = [plt.Line2D([0], [0], marker="s", linestyle="",
                          color=REVERSAL_TRANSITION_COLORS[item],
                          label=REVERSAL_TRANSITION_LABELS[item]) for item in classes]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle("Where cell-volume scaling changes reversal behavior")
    transition_path = output / "reversal_transition_maps.png"
    fig.savefig(transition_path, dpi=210)
    plt.close(fig)
    return [outcome_path, transition_path]


def summarize_dual(rows, args):
    lookup = paired(rows, ("epsilon",))
    keys = sorted({key[1:] for key in lookup if key[0] == 1.0})
    summary = {"analysis": "wgd_null_two_seeded_domains", "simulations": len(rows),
               "domain_count_times_seconds": values(args.domain_count_times),
               "transitions": {}, "time_sequences": {}}
    for volume in values(args.volume_multipliers):
        for epsilon in values(args.asymmetries):
            sequences = {}
            for h, eps, mu, length in keys:
                if not np.isclose(eps, epsilon):
                    continue
                row = lookup[(volume, h, eps, mu, length)]
                sequence = "→".join(str(row[f"polarized_domain_count_t{time:g}"])
                                    for time in values(args.domain_count_times))
                sequences[sequence] = sequences.get(sequence, 0) + 1
            summary["time_sequences"][f"volume={volume:g},epsilon={epsilon:g}"] = sequences
    for destination_volume in [v for v in values(args.volume_multipliers) if not np.isclose(v, 1)]:
        for epsilon in values(args.asymmetries):
            exact, grouped = {}, {}
            for h, eps, mu, length in keys:
                if not np.isclose(eps, epsilon):
                    continue
                source = int(lookup[(1.0, h, eps, mu, length)]["polarized_domain_count_t300"])
                destination = int(lookup[(destination_volume, h, eps, mu, length)]["polarized_domain_count_t300"])
                exact[f"{source}→{destination}"] = exact.get(f"{source}→{destination}", 0) + 1
                group = domain_transition(source, destination)
                grouped[group] = grouped.get(group, 0) + 1
            summary["transitions"][f"1x_to_{destination_volume:g}x,epsilon={epsilon:g}"] = {
                "grouped": grouped, "exact": exact,
            }
    summary["maximum_conservation_error"] = max(float(row["maximum_conservation_error"]) for row in rows)
    return summary


def summarize_reversal(rows, args):
    lookup = paired(rows)
    keys = sorted({key[1:] for key in lookup if key[0] == 1.0})
    outcome_counts = {}
    for volume in values(args.volume_multipliers):
        counts = {}
        for h, mu, length in keys:
            outcome = str(lookup[(volume, h, mu, length)]["reversal_outcome"])
            counts[outcome] = counts.get(outcome, 0) + 1
        outcome_counts[f"volume={volume:g}"] = counts
    transitions = {}
    for destination_volume in [v for v in values(args.volume_multipliers) if not np.isclose(v, 1)]:
        exact, grouped = {}, {}
        for h, mu, length in keys:
            source = str(lookup[(1.0, h, mu, length)]["reversal_outcome"])
            destination = str(lookup[(destination_volume, h, mu, length)]["reversal_outcome"])
            exact[f"{source} -> {destination}"] = exact.get(f"{source} -> {destination}", 0) + 1
            group = reversal_transition(source, destination)
            grouped[group] = grouped.get(group, 0) + 1
        transitions[f"1x_to_{destination_volume:g}x"] = {"grouped": grouped, "exact": exact}
    return {
        "analysis": "wgd_null_polarity_reversal", "simulations": len(rows),
        "outcome_counts": outcome_counts, "transitions": transitions,
        "maximum_conservation_error": max(float(row["maximum_conservation_error"]) for row in rows),
    }


def report_fragment(dual_summary, reversal_summary, dual_figures, reversal_figures):
    dual_rows = []
    for comparison, records in dual_summary["transitions"].items():
        changed = sum(records["grouped"].values()) - records["grouped"].get("unchanged", 0)
        dual_rows.append(
            f"<tr><td>{escape(comparison)}</td><td>{changed}</td>"
            f"<td>{records['grouped'].get('one_to_two', 0)}</td>"
            f"<td>{records['grouped'].get('two_to_one', 0)}</td></tr>"
        )
    reversal_rows = []
    for comparison, records in reversal_summary["transitions"].items():
        changed = sum(records["grouped"].values()) - records["grouped"].get("unchanged", 0)
        reversal_rows.append(
            f"<tr><td>{escape(comparison)}</td><td>{changed}</td>"
            f"<td>{records['grouped'].get('became_reversible', 0)}</td>"
            f"<td>{records['grouped'].get('lost_reversibility', 0)}</td></tr>"
        )
    dual_html = "".join(
        f"<figure><a href='double_stimulus/{path.name}'><img src='double_stimulus/{path.name}' loading='lazy'></a>"
        f"<figcaption><strong>{escape(path.stem.replace('_', ' ').title())}</strong>At 300 s, each pixel compares the scaled cell with its matched 1× condition.</figcaption></figure>"
        for path in dual_figures
    )
    reversal_html = "".join(
        f"<figure><a href='polarity_reversal/{path.name}'><img src='polarity_reversal/{path.name}' loading='lazy'></a>"
        f"<figcaption><strong>{escape(path.stem.replace('_', ' ').title())}</strong></figcaption></figure>"
        for path in reversal_figures
    )
    return f"""
<aside class="result"><strong>Common interpretation across the three
experiments.</strong> Increasing cell size makes polarity domains more robust in
this model: larger cells more readily establish or retain spatial polarity, and
established patterns are less readily replaced by competition or an opposing
cue. “More robust” applies to polarity domains generally, not exclusively to a
single winning pole—larger cells can also preserve multiple competing regions.</aside>
<section class="analysis"><h2>WGD null: two seeded competing domains</h2>
<p class="section-intro">This repeats the standardized simultaneous-competition
test at 0.5×, 1×, and 2× volume. Two mass-conserving active enrichments are seeded
at 20% and 80% of the strip at <code>t=0</code>; they are initial conditions, not
time-limited external cues. Their amplitude is 2, fractional width is 0.04, and
left/right imbalance is <code>ε=0, 0.001, 0.01, 0.05</code>. The direct output is
the number of disconnected polarized regions at 0, 100, and 300 s.</p>
<aside class="note"><strong>Interpretation.</strong> Scaling down to 0.5× causes
2 → 1 region transitions at every tested <code>ε</code>, increasingly often as
the initial imbalance grows. Scaling up to 2× produces no 2 → 1 transitions and
instead produces 1 → 2 transitions at <code>ε=0.01</code> and <code>0.05</code>.
Thus larger cells are better able to preserve competing polarity regions rather
than forcing rapid winner-take-all elimination.</aside>
<div class="table-wrap"><table><thead><tr><th>Comparison</th><th>Changed terminal count</th><th>1 → 2</th><th>2 → 1</th></tr></thead><tbody>{''.join(dual_rows)}</tbody></table></div>
<div class="figures">{dual_html}</div>
<p><a class="report-link" href="double_stimulus/two_seeded_domains.csv">Download two-domain results</a>
<a class="report-link" href="double_stimulus/summary.json">Download two-domain summary</a></p></section>
<section class="analysis"><h2>WGD null: polarity reversal</h2>
<p class="section-intro">This repeats the sequential-cue experiment at all three
volumes. A left cue acts from 10–30 s and an equal right cue from 150–170 s. A
condition is eligible only if it has one clearly left-biased region immediately
before the right cue. “Reversed” means it ends with one clearly right-biased
region; gray ineligibility is not a failed reversal.</p>
<aside class="note"><strong>Interpretation.</strong> At 2× volume, clean
reversals fall from 41 to 35 while multiple-region endpoints rise from 14 to 37.
Thirteen conditions that reverse at 1× instead retain multiple regions at 2×.
The opposing cue therefore less often replaces the established polarity pattern
with one clean new pole, consistent with stronger persistence of polarity
domains in larger cells.</aside>
<div class="table-wrap"><table><thead><tr><th>Comparison</th><th>Changed reversal outcome</th><th>Became reversible</th><th>Lost reversibility</th></tr></thead><tbody>{''.join(reversal_rows)}</tbody></table></div>
<div class="figures">{reversal_html}</div>
<p><a class="report-link" href="polarity_reversal/polarity_reversal.csv">Download reversal results</a>
<a class="report-link" href="polarity_reversal/summary.json">Download reversal summary</a></p></section>
"""


def merge(args) -> Path:
    work = output_root(args)
    root = report_root(args)
    tasks = all_tasks(args)
    paths = [work / "extended_tasks" / f"task_{index:05d}.json" for index in range(len(tasks))]
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing {len(missing)} extended task results; first is {missing[0]}")
    rows = [json.loads(path.read_text()) for path in paths]
    dual_rows = [row for row in rows if row["experiment"] == "two_seeded_domains"]
    reversal_rows = [row for row in rows if row["experiment"] == "polarity_reversal"]
    expected_dual = (len(values(args.volume_multipliers)) * len(values(args.h_values))
                     * len(values(args.asymmetries)) * args.dual_mu_count
                     * args.dual_lambda_count)
    expected_reversal = (len(values(args.volume_multipliers)) * len(values(args.h_values))
                         * args.reversal_mu_count * args.reversal_lambda_count)
    if len(dual_rows) != expected_dual or len(reversal_rows) != expected_reversal:
        raise ValueError(f"unexpected row counts: dual={len(dual_rows)}, reversal={len(reversal_rows)}")
    dual_output = root / "double_stimulus"
    reversal_output = root / "polarity_reversal"
    write_csv(dual_output / "two_seeded_domains.csv", dual_rows)
    write_csv(reversal_output / "polarity_reversal.csv", reversal_rows)
    dual_summary = summarize_dual(dual_rows, args)
    reversal_summary = summarize_reversal(reversal_rows, args)
    dual_output.mkdir(parents=True, exist_ok=True)
    reversal_output.mkdir(parents=True, exist_ok=True)
    dual_figures = render_dual_transition_maps(dual_output, dual_rows, args)
    reversal_figures = render_reversal_maps(reversal_output, reversal_rows, args)
    (dual_output / "summary.json").write_text(json.dumps(dual_summary, indent=2) + "\n")
    (reversal_output / "summary.json").write_text(json.dumps(reversal_summary, indent=2) + "\n")
    (root / "extended_report_fragment.html").write_text(
        report_fragment(dual_summary, reversal_summary, dual_figures, reversal_figures)
    )
    single_summary = json.loads((root / "summary.json").read_text())
    single_figures = [root / name for name in single_summary["figures"]]
    report = write_report(root, single_summary, single_figures)
    print(report)
    return report


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--mode", choices=("task", "merge"), required=True)
    result.add_argument("--task-index", type=int, default=0)
    result.add_argument("--output", type=Path, default=Path("migration1d/results/single_species/testing/wgd_null"))
    result.add_argument("--report-root", type=Path)
    result.add_argument("--volume-multipliers", default="0.5,1,2")
    result.add_argument("--h-values", default="1,2,4")
    result.add_argument("--dual-mu-min", type=float, default=0.5)
    result.add_argument("--dual-mu-max", type=float, default=5.0)
    result.add_argument("--dual-mu-count", type=int, default=13)
    result.add_argument("--dual-lambda-min", type=float, default=10.0)
    result.add_argument("--dual-lambda-max", type=float, default=100.0)
    result.add_argument("--dual-lambda-count", type=int, default=13)
    result.add_argument("--reversal-mu-min", type=float, default=0.5)
    result.add_argument("--reversal-mu-max", type=float, default=5.0)
    result.add_argument("--reversal-mu-count", type=int, default=17)
    result.add_argument("--reversal-lambda-min", type=float, default=10.0)
    result.add_argument("--reversal-lambda-max", type=float, default=100.0)
    result.add_argument("--reversal-lambda-count", type=int, default=17)
    result.add_argument("--asymmetries", default="0,0.001,0.01,0.05")
    result.add_argument("--domain-count-times", default="0,100,300")
    result.add_argument("--initial-amplitude", type=float, default=2.0)
    result.add_argument("--initial-width", type=float, default=0.04)
    result.add_argument("--length", type=float, default=10.0)
    result.add_argument("--rho", type=float, default=100.0)
    result.add_argument("--kappa", type=float, default=0.067)
    result.add_argument("--gamma", type=float, default=1.0)
    result.add_argument("--hill", type=float, default=2.0)
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
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    run_task(args) if args.mode == "task" else merge(args)


if __name__ == "__main__":
    main()
