#!/usr/bin/env python3
"""Analyze the decisive WP1 campaign and produce manuscript-oriented outputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_NAME = "wp1_decisive_v1"
VARIANTS = (
    "intact",
    "rac_to_rho_half",
    "rho_to_rac_half",
    "rac_to_rho_off",
    "rho_to_rac_off",
    "both_off",
)
VARIANT_LABELS = {
    "intact": "Both arms intact",
    "rac_to_rho_half": "Rac→Rho 50%",
    "rho_to_rac_half": "Rho→Rac 50%",
    "rac_to_rho_off": "Rac→Rho off",
    "rho_to_rac_off": "Rho→Rac off",
    "both_off": "Both arms off",
}
VARIANT_COLORS = dict(zip(VARIANTS, plt.cm.viridis(np.linspace(0.05, 0.95, 6))))
MODEL_LABELS = {
    "jilkine_wave_pinning": "Wave-pinning",
    "jilkine_otsuji": "Otsuji",
    "jilkine_goryachev": "Goryachev",
    "jilkine_legi": "LEGI",
}
MODEL_COLORS = dict(zip(MODEL_LABELS, plt.cm.Set2(np.linspace(0.05, 0.85, 4))))


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _number(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    if value in {"", "None", "null", "not_applicable"}:
        return None
    try:
        result = float(value)
    except ValueError:
        return None
    return result if np.isfinite(result) else None


def _boolean(row: dict[str, str], key: str) -> bool | None:
    value = row.get(key, "").lower()
    if value == "true":
        return True
    if value == "false":
        return False
    return None


def _format_length(value: float) -> str:
    return f"{value:g}" if float(value).is_integer() else f"{value:.1f}"


def _clean(rows: list[dict[str, str]]) -> tuple[list[dict], list[dict], list[dict]]:
    competition, reversal, propagation = [], [], []
    for row in rows:
        length = _number(row, "length_um")
        variant = row["coupling_variant"]
        assay = row["assay"]
        if assay == "dual_site_fixed_physical_width_extended":
            establishment = _number(row, "competition_establishment_time_s")
            end_time = _number(row, "observation_end_time_s")
            latency = _number(row, "competition_resolution_latency_s")
            censor = (
                None
                if establishment is None or end_time is None
                else end_time - establishment
            )
            competition.append(
                {
                    "task_id": row["task_id"],
                    "coupling_variant": variant,
                    "length_um": length,
                    "competition_status": row["competition_status"],
                    "competition_effective_status": row["competition_effective_status"],
                    "competition_establishment_time_s": establishment,
                    "competition_resolution_latency_s": latency,
                    "competition_resolution_latency_over_L2": (
                        None if latency is None else latency / length**2
                    ),
                    "resolution_censor_limit_s": censor,
                    "resolution_censor_limit_over_L2": (
                        None if censor is None else censor / length**2
                    ),
                    "persistent_resolution": _boolean(
                        row, "competition_resolution_persistent_at_end"
                    ),
                    "transient_resolution_only": _boolean(
                        row, "competition_resolution_transient_only"
                    ),
                    "endpoint_state": row["endpoint_state"],
                    "endpoint_state_sustained": _boolean(row, "endpoint_state_sustained"),
                    "observation_end_time_s": end_time,
                }
            )
        elif assay == "gradient_reversal_length_scaled":
            reversal_latency = _number(row, "orientation_reversal_latency_s")
            pre_reversal_time = _number(row, "pre_reversal_time_s")
            observation_end = _number(row, "observation_end_time_s")
            reversal.append(
                {
                    "task_id": row["task_id"],
                    "coupling_variant": variant,
                    "input_peak": _number(row, "input_peak"),
                    "length_um": length,
                    "orientation_reversal_status": row["orientation_reversal_status"],
                    "orientation_reversal_latency_s": reversal_latency,
                    "orientation_reversal_latency_over_L2": (
                        None if reversal_latency is None else reversal_latency / length**2
                    ),
                    "reversal_censor_limit_s": (
                        None
                        if pre_reversal_time is None or observation_end is None
                        else observation_end - (pre_reversal_time + 2.0)
                    ),
                    "old_front_decay_latency_s": _number(
                        row, "old_front_decay_latency_s"
                    ),
                    "pre_reversal_state": row.get("pre_reversal_state", ""),
                    "pre_reversal_domain_count": _number(
                        row, "pre_reversal_domain_count"
                    ),
                    "pre_reversal_spatial_contrast": _number(
                        row, "pre_reversal_spatial_contrast"
                    ),
                    "pre_reversal_orientation": _number(
                        row, "pre_reversal_orientation"
                    ),
                    "observation_end_time_s": observation_end,
                }
            )
        elif assay == "single_site_fixed_physical_width":
            distal_peak_time = _number(row, "rac_distal_peak_time_s")
            observation_end = _number(row, "observation_end_time_s")
            propagation.append(
                {
                    "task_id": row["task_id"],
                    "coupling_variant": variant,
                    "length_um": length,
                    "rac_proximal_peak_change": _number(
                        row, "rac_proximal_peak_change"
                    ),
                    "rac_distal_peak_change": _number(row, "rac_distal_peak_change"),
                    "rac_distal_peak_time_s": distal_peak_time,
                    "rac_distal_peak_at_observation_end": bool(
                        distal_peak_time is not None
                        and observation_end is not None
                        and np.isclose(distal_peak_time, observation_end)
                    ),
                    "observation_end_time_s": observation_end,
                    "rac_distal_to_proximal_abs_ratio": _number(
                        row, "rac_distal_to_proximal_abs_ratio"
                    ),
                    "rho_distal_peak_change": _number(row, "rho_distal_peak_change"),
                    "maximum_absolute_tension": _number(
                        row, "maximum_absolute_tension"
                    ),
                }
            )
    competition.sort(key=lambda row: (VARIANTS.index(row["coupling_variant"]), row["length_um"]))
    reversal.sort(
        key=lambda row: (
            row["input_peak"], VARIANTS.index(row["coupling_variant"]), row["length_um"]
        )
    )
    propagation.sort(key=lambda row: (VARIANTS.index(row["coupling_variant"]), row["length_um"]))
    return competition, reversal, propagation


def _plot_competition(rows: list[dict], destination: Path) -> None:
    lengths = sorted({row["length_um"] for row in rows})
    statuses = sorted({row["competition_effective_status"] for row in rows})
    preferred = [
        "not_established",
        "resolved_persistent",
        "resolved_transient_to_single_front",
        "resolved_transient_to_polarized_no_resolved_peak",
        "resolved_transient_to_unpolarized",
        "unresolved_persistent_multifront",
        "unresolved_to_polarized_no_resolved_peak",
        "unresolved_to_unpolarized",
    ]
    statuses = [status for status in preferred if status in statuses] + [
        status for status in statuses if status not in preferred
    ]
    palette = plt.cm.tab10(np.linspace(0.05, 0.95, max(len(statuses), 2)))
    matrix = np.zeros((len(VARIANTS), len(lengths)), dtype=int)
    for row in rows:
        matrix[VARIANTS.index(row["coupling_variant"]), lengths.index(row["length_um"])] = statuses.index(
            row["competition_effective_status"]
        )
    figure, axes = plt.subplots(2, 1, figsize=(13, 8.5), gridspec_kw={"height_ratios": [0.85, 1.2]})
    axes[0].imshow(
        matrix,
        aspect="auto",
        cmap=ListedColormap(palette[: len(statuses)]),
        norm=BoundaryNorm(np.arange(-0.5, len(statuses) + 0.5), len(statuses)),
    )
    axes[0].set_yticks(range(len(VARIANTS)), [VARIANT_LABELS[value] for value in VARIANTS])
    axes[0].set_xticks(range(len(lengths)), [_format_length(value) for value in lengths])
    axes[0].set_title("A  Long-time competition outcome", loc="left", fontweight="bold")
    for i, variant in enumerate(VARIANTS):
        for row in [item for item in rows if item["coupling_variant"] == variant]:
            axes[0].text(
                lengths.index(row["length_um"]), i,
                {"single_front": "S", "multifront": "M", "unpolarized": "U", "polarized_no_resolved_peak": "P"}[row["endpoint_state"]],
                ha="center", va="center", fontsize=8, fontweight="bold",
            )
    for variant in VARIANTS:
        selected = [row for row in rows if row["coupling_variant"] == variant]
        resolved = [row for row in selected if row["competition_resolution_latency_s"] is not None]
        axes[1].plot(
            [row["length_um"] for row in resolved],
            [row["competition_resolution_latency_s"] for row in resolved],
            "-", color=VARIANT_COLORS[variant], lw=1.8, label=VARIANT_LABELS[variant],
        )
        for row in resolved:
            axes[1].scatter(
                row["length_um"], row["competition_resolution_latency_s"],
                s=42, marker="o", edgecolor=VARIANT_COLORS[variant],
                facecolor=VARIANT_COLORS[variant] if row["persistent_resolution"] else "white",
                zorder=3,
            )
        censored = [row for row in selected if row["competition_resolution_latency_s"] is None and row["resolution_censor_limit_s"] is not None]
        axes[1].scatter(
            [row["length_um"] for row in censored],
            [row["resolution_censor_limit_s"] for row in censored],
            marker="v", facecolor="none", edgecolor=VARIANT_COLORS[variant], s=58,
        )
    axes[1].set_xlabel("Cell length (µm)")
    axes[1].set_ylabel("Two-front resolution latency (s)")
    axes[1].set_title(
        "B  Filled circle: persistent at run end; open circle: transient; triangle: censored",
        loc="left", fontweight="bold",
    )
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False, ncol=2, fontsize=9)
    legend_text = "; ".join(f"{index}={status}" for index, status in enumerate(statuses))
    figure.text(0.01, 0.01, legend_text, fontsize=7)
    figure.suptitle("Fixed-width De Belly competition across mechanochemical coupling states")
    figure.tight_layout(rect=(0, 0.05, 1, 0.95))
    figure.savefig(destination, dpi=220)
    plt.close(figure)


def _plot_reversal(rows: list[dict], destination: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(14, 8.5), sharex="col")
    for column, input_peak in enumerate((0.2, 0.5)):
        selected_input = [row for row in rows if np.isclose(row["input_peak"], input_peak)]
        intact = {
            row["length_um"]: row["orientation_reversal_latency_s"]
            for row in selected_input
            if row["coupling_variant"] == "intact"
        }
        for variant in VARIANTS:
            selected = [
                row for row in selected_input
                if row["coupling_variant"] == variant
                and row["orientation_reversal_latency_s"] is not None
            ]
            axes[0, column].plot(
                [row["length_um"] for row in selected],
                [row["orientation_reversal_latency_s"] for row in selected],
                "o-", color=VARIANT_COLORS[variant], lw=1.6, ms=4,
                label=VARIANT_LABELS[variant],
            )
            censored = [
                row for row in selected_input
                if row["coupling_variant"] == variant
                and row["orientation_reversal_status"] == "right_censored"
            ]
            if censored:
                plotted_values = [
                    item["orientation_reversal_latency_s"]
                    for item in selected_input
                    if item["orientation_reversal_latency_s"] is not None
                ]
                marker_height = 1.08 * max(plotted_values)
                axes[0, column].scatter(
                    [row["length_um"] for row in censored],
                    [marker_height] * len(censored),
                    marker="x", color=VARIANT_COLORS[variant], s=42, zorder=4,
                )
            delta_rows = [row for row in selected if intact.get(row["length_um"]) is not None]
            axes[1, column].plot(
                [row["length_um"] for row in delta_rows],
                [row["orientation_reversal_latency_s"] - intact[row["length_um"]] for row in delta_rows],
                "o-", color=VARIANT_COLORS[variant], lw=1.6, ms=4,
            )
        axes[0, column].set_title(f"Input peak {input_peak:g}")
        axes[0, column].set_ylabel("Side-switch latency (s)")
        axes[1, column].set_ylabel("Latency difference from intact (s)")
        axes[1, column].set_xlabel("Cell length (µm)")
        axes[1, column].axhline(0, color="black", lw=0.8)
        for axis in axes[:, column]:
            axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    figure.suptitle(
        "Length-scaled De Belly reversal: size × arm-state interaction\n"
        "Conditioning and follow-up scale with L²; side-orientation endpoint"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    figure.savefig(destination, dpi=220)
    plt.close(figure)


def _median_by(rows: list[dict], keys: tuple[str, ...], metric: str) -> list[dict]:
    groups: dict[tuple, list[float]] = {}
    for row in rows:
        value = _number(row, metric)
        if value is not None:
            groups.setdefault(tuple(row[key] for key in keys), []).append(value)
    return [
        {**dict(zip(keys, key)), metric: median(values)}
        for key, values in groups.items()
    ]


def _plot_diffusive_rescaling(
    competition: list[dict], reversal: list[dict], phase1_rows: list[dict[str, str]], destination: Path
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(14, 8.5))
    old_comp = [
        row for row in phase1_rows
        if row["assay"] == "competing_fronts" and _number(row, "competition_resolution_latency_s") is not None
    ]
    old_rev = [row for row in phase1_rows if row["assay"] == "gradient_reversal"]
    old_comp_medians = _median_by(old_comp, ("model", "length_um"), "competition_resolution_latency_s")
    for model, label in MODEL_LABELS.items():
        selected = sorted(
            [row for row in old_comp_medians if row["model"] == model],
            key=lambda row: float(row["length_um"]),
        )
        lengths = np.asarray([float(row["length_um"]) for row in selected])
        values = np.asarray([row["competition_resolution_latency_s"] for row in selected])
        axes[0, 0].plot(lengths, values / lengths**2, "o-", ms=3.5, color=MODEL_COLORS[model], label=label)
    for model, label in MODEL_LABELS.items():
        selected = sorted(
            [
                row for row in old_rev
                if row["model"] == model
                and (
                    _number(row, "orientation_reversal_latency_s") is not None
                    or _number(row, "reversal_latency_s") is not None
                )
            ],
            key=lambda row: _number(row, "length_um"),
        )
        lengths = np.asarray([_number(row, "length_um") for row in selected])
        values = np.asarray([
            _number(row, "orientation_reversal_latency_s")
            if _number(row, "orientation_reversal_latency_s") is not None
            else _number(row, "reversal_latency_s")
            for row in selected
        ])
        axes[1, 0].plot(lengths, values / lengths**2, "o-", ms=3.5, color=MODEL_COLORS[model], label=label)
    for variant in VARIANTS:
        selected = [row for row in competition if row["competition_resolution_latency_over_L2"] is not None and row["coupling_variant"] == variant]
        axes[0, 1].plot(
            [row["length_um"] for row in selected],
            [row["competition_resolution_latency_over_L2"] for row in selected],
            "o-", ms=3.5, color=VARIANT_COLORS[variant], label=VARIANT_LABELS[variant],
        )
        selected = [
            row for row in reversal
            if np.isclose(row["input_peak"], 0.2)
            and row["coupling_variant"] == variant
            and row["orientation_reversal_latency_over_L2"] is not None
        ]
        axes[1, 1].plot(
            [row["length_um"] for row in selected],
            [row["orientation_reversal_latency_over_L2"] for row in selected],
            "o-", ms=3.5, color=VARIANT_COLORS[variant], label=VARIANT_LABELS[variant],
        )
    titles = (
        "A  Jilkine competition",
        "B  De Belly competition",
        "C  Jilkine reversal",
        "D  De Belly reversal, input 0.2",
    )
    for axis, title in zip(axes.flat, titles):
        axis.set_title(title, loc="left", fontweight="bold")
        axis.set_xlabel("Cell length (µm)")
        axis.set_ylabel("Latency / L² (s/µm²)")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 1].legend(frameon=False, fontsize=7, ncol=2)
    figure.suptitle("Diffusive-time rescaling diagnostic")
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    figure.savefig(destination, dpi=220)
    plt.close(figure)


def _plot_propagation(rows: list[dict], destination: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    for variant in VARIANTS:
        selected = [row for row in rows if row["coupling_variant"] == variant]
        axes[0].plot(
            [row["length_um"] for row in selected],
            [abs(row["rac_distal_peak_change"]) for row in selected],
            "o-", color=VARIANT_COLORS[variant], lw=1.7, ms=4, label=VARIANT_LABELS[variant],
        )
        axes[1].plot(
            [row["length_um"] for row in selected],
            [row["rac_distal_to_proximal_abs_ratio"] for row in selected],
            "o-", color=VARIANT_COLORS[variant], lw=1.7, ms=4,
        )
        endpoint_peaks = [
            row for row in selected if row["rac_distal_peak_at_observation_end"]
        ]
        axes[0].scatter(
            [row["length_um"] for row in endpoint_peaks],
            [abs(row["rac_distal_peak_change"]) for row in endpoint_peaks],
            marker="v", facecolor="white", edgecolor=VARIANT_COLORS[variant],
            s=48, zorder=4,
        )
    axes[0].set_ylabel("Absolute distal Rac peak change")
    axes[1].set_ylabel("Distal/proximal absolute Rac response")
    for axis in axes:
        axis.set_xlabel("Cell length (µm)")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False, fontsize=8, ncol=2)
    figure.suptitle(
        "Fixed-physical-width single-site Rac propagation (σ = 1.6 µm)\n"
        "Open triangles mark responses whose largest observed magnitude occurs at run end"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    figure.savefig(destination, dpi=220)
    plt.close(figure)


def _analysis_summary(competition: list[dict], reversal: list[dict], propagation: list[dict]) -> dict:
    competition_counts = {}
    for variant in VARIANTS:
        selected = [row for row in competition if row["coupling_variant"] == variant]
        competition_counts[variant] = {
            "persistent_resolution_count": sum(bool(row["persistent_resolution"]) for row in selected),
            "transient_resolution_count": sum(bool(row["transient_resolution_only"]) for row in selected),
            "right_censored_count": sum(row["competition_status"] == "right_censored" for row in selected),
            "condition_count": len(selected),
            "persistent_resolution_lengths_um": [
                row["length_um"] for row in selected if row["persistent_resolution"]
            ],
            "right_censored_lengths_um": [
                row["length_um"]
                for row in selected if row["competition_status"] == "right_censored"
            ],
        }
    reversal_counts = {}
    for input_peak in (0.2, 0.5):
        reversal_counts[str(input_peak)] = {}
        for variant in VARIANTS:
            selected = [row for row in reversal if np.isclose(row["input_peak"], input_peak) and row["coupling_variant"] == variant]
            reversal_counts[str(input_peak)][variant] = {
                "success_count": sum(row["orientation_reversal_status"] == "success" for row in selected),
                "right_censored_count": sum(
                    row["orientation_reversal_status"] == "right_censored"
                    for row in selected
                ),
                "condition_count": len(selected),
                "right_censored_lengths_um": [
                    row["length_um"]
                    for row in selected
                    if row["orientation_reversal_status"] == "right_censored"
                ],
            }
    propagation_contrasts = []
    for length in sorted({row["length_um"] for row in propagation}):
        lookup = {
            row["coupling_variant"]: abs(row["rac_distal_peak_change"])
            for row in propagation if row["length_um"] == length
        }
        propagation_contrasts.append(
            {
                "length_um": length,
                "intact_abs_distal_rac": lookup["intact"],
                "both_off_abs_distal_rac": lookup["both_off"],
                "both_off_below_1e-10_numerical_floor": lookup["both_off"] < 1e-10,
                "intact_to_both_off_fold": (
                    None
                    if lookup["both_off"] < 1e-10
                    else lookup["intact"] / lookup["both_off"]
                ),
            }
        )
    return {
        "competition": competition_counts,
        "reversal": reversal_counts,
        "fixed_width_propagation_peak_at_observation_end_count": sum(
            bool(row["rac_distal_peak_at_observation_end"]) for row in propagation
        ),
        "fixed_width_propagation_intact_vs_both_off": propagation_contrasts,
    }


def _report(output_root: Path, collection: dict, summary: dict) -> None:
    text = [
        "# Decisive WP1 campaign",
        "",
        "## Design",
        "",
        f"All {collection['success_count']}/{collection['task_count']} deterministic tasks completed successfully.",
        "The campaign crosses 14 lengths with six reciprocal-arm states for long-time",
        "fixed-width competition, length-scaled reversal at input peaks 0.2 and 0.5,",
        "and fixed-width single-site propagation. Gaussian inputs use σ = 1.6 µm.",
        "",
        "Competition outcomes explicitly distinguish persistent one-front resolution at",
        "run end from transient passage through a one-front state and from right-censoring.",
        "Counts describe a deterministic design grid, not replicate probabilities.",
        "",
        "## Competition outcome counts",
        "",
        "| Coupling state | Persistent | Transient | Right-censored | Total |",
        "|---|---:|---:|---:|---:|",
    ]
    for variant in VARIANTS:
        item = summary["competition"][variant]
        text.append(
            f"| {VARIANT_LABELS[variant]} | {item['persistent_resolution_count']} | "
            f"{item['transient_resolution_count']} | {item['right_censored_count']} | "
            f"{item['condition_count']} |"
        )
    rac_off_reversal = summary["reversal"]["0.2"]["rac_to_rho_off"]
    text.extend(
        [
            "",
            "## Main findings",
            "",
            "- Intact coupling produced a sustained one-front episode at every length,",
            "  but it remained a persistent final one-front state only at L=20 and 25.2 µm.",
            "  Smaller systems later depolarized and larger systems ended polarized without",
            "  a detector-resolved peak. Thus event resolution and long-time organization",
            "  are distinct endpoints.",
            "- With the Rac→Rho arm off, L=12–36 µm systems depolarized without resolving,",
            "  L=40–44 µm systems resolved persistently after long delays, and L=48–64 µm",
            "  systems retained persistent competing fronts. The interaction is arm-specific",
            "  and non-monotone rather than a scalar effect of total coupling strength.",
            "- At reversal input 0.2, Rac→Rho-off systems switched only through L=25.2 µm;",
            f"  {rac_off_reversal['right_censored_count']}/14 were right-censored at "
            f"{rac_off_reversal['right_censored_lengths_um']}. All arm states switched at",
            "  input 0.5, making the large-cell dependency conditional on cue strength.",
            "- Dividing latency by L² does not collapse the competition curves. Reversal",
            "  latencies generally grow much more slowly than L², while endpoint classes",
            "  change qualitatively in the competition assay.",
            f"- In the fixed-width propagation control, {summary['fixed_width_propagation_peak_at_observation_end_count']}/84",
            "  distal-response maxima occurred at the observation boundary even after",
            "  length-scaled extension. These values are flagged as lower bounds and should",
            "  not be presented as resolved peak amplitudes.",
            "",
            "## Outputs",
            "",
            "- `competition_hardened.csv`: persistent/transient/censored competition endpoints.",
            "- `reversal_length_scaled.csv`: side-switch latency and conditioning morphology.",
            "- `propagation_fixed_width.csv`: fixed-width distal Rac/Rho response metrics.",
            "- `wp1_competition_interaction.png`: long-time size × coupling outcomes.",
            "- `wp1_reversal_interaction.png`: length-scaled reversal interaction.",
            "- `wp1_diffusive_rescaling.png`: absolute-size effects expressed on the L² clock.",
            "- `wp1_fixed_width_propagation.png`: fixed-width mechanochemical propagation control.",
            "- `analysis_summary.json`: machine-readable outcome counts and matched contrasts.",
            "",
        ]
    )
    (output_root / "README.md").write_text("\n".join(text))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_root", nargs="?", type=Path,
        default=ROOT / "results" / "campaigns" / CAMPAIGN_NAME,
    )
    parser.add_argument(
        "--phase1-root", type=Path,
        default=ROOT / "results" / "campaigns" / "phase1_production_v1",
    )
    arguments = parser.parse_args(argv)
    collection = json.loads((arguments.output_root / "collection_summary.json").read_text())
    if not collection.get("complete"):
        raise SystemExit("campaign collection is incomplete")
    competition, reversal, propagation = _clean(_read(arguments.output_root / "task_summary.csv"))
    if (len(competition), len(reversal), len(propagation)) != (84, 168, 84):
        raise SystemExit(
            f"unexpected task split: {len(competition)} competition, "
            f"{len(reversal)} reversal, {len(propagation)} propagation"
        )
    _write(arguments.output_root / "competition_hardened.csv", competition)
    _write(arguments.output_root / "reversal_length_scaled.csv", reversal)
    _write(arguments.output_root / "propagation_fixed_width.csv", propagation)
    _plot_competition(competition, arguments.output_root / "wp1_competition_interaction.png")
    _plot_reversal(reversal, arguments.output_root / "wp1_reversal_interaction.png")
    _plot_diffusive_rescaling(
        competition,
        reversal,
        _read(arguments.phase1_root / "task_summary.csv"),
        arguments.output_root / "wp1_diffusive_rescaling.png",
    )
    _plot_propagation(propagation, arguments.output_root / "wp1_fixed_width_propagation.png")
    summary = _analysis_summary(competition, reversal, propagation)
    (arguments.output_root / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    _report(arguments.output_root, collection, summary)
    print(arguments.output_root)


if __name__ == "__main__":
    main()
