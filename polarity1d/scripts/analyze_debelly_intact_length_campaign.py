#!/usr/bin/env python3
"""Analyze and plot the intact-first De Belly length-effects campaign."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_NAME = "debelly_intact_length_v2"
WIDTH_ASSAYS = {
    "dual_site_fractional_width": "σ = 0.08L",
    "dual_site_fixed_physical_width": "σ = 1.6 µm",
}
WIDTH_COLORS = {
    "dual_site_fractional_width": "#6f4c9b",
    "dual_site_fixed_physical_width": "#1984a0",
}
INPUT_COLORS = {0.2: "#a0a0a0", 0.5: "#d95f02", 1.0: "#1b9e77"}
COMPETITION_ORDER = ("not_established", "resolved", "right_censored")
COMPETITION_COLORS = ("#d8d8d8", "#2878b5", "#e68a2e")
COMPETITION_LABELS = {
    "not_established": "Two fronts never sustained",
    "resolved": "Two fronts → one sustained front",
    "right_censored": "No observed one-front resolution",
}
REVERSAL_ORDER = ("not_conditioned", "success", "right_censored")
REVERSAL_COLORS = ("#d8d8d8", "#2a9d5b", "#d95f45")
REVERSAL_LABELS = {
    "not_conditioned": "No qualified old-side orientation",
    "success": "Sustained side switch",
    "right_censored": "No sustained side switch by run end",
}


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


def _clean(rows: list[dict[str, str]]) -> tuple[list[dict], list[dict]]:
    dual, reversal = [], []
    for row in rows:
        assay = row["assay"]
        length = _number(row, "length_um")
        if assay in WIDTH_ASSAYS:
            establishment = _number(row, "competition_establishment_time_s")
            observation_end = _number(row, "observation_end_time_s")
            censor_limit = None
            if establishment is not None and observation_end is not None:
                censor_limit = observation_end - establishment
            dual.append(
                {
                    "task_id": row["task_id"],
                    "assay": assay,
                    "width_mode": WIDTH_ASSAYS[assay],
                    "length_um": length,
                    "competition_status": row["competition_status"],
                    "competition_establishment_time_s": establishment,
                    "competition_resolution_latency_s": _number(
                        row, "competition_resolution_latency_s"
                    ),
                    "resolution_censor_limit_s": censor_limit,
                    "endpoint_state": row["endpoint_state"],
                    "observation_end_time_s": observation_end,
                }
            )
        elif assay == "gradient_reversal_length_scaled":
            conditioned = _boolean(row, "orientation_reversal_conditioned")
            success = _boolean(row, "orientation_reversal_success")
            censored = _boolean(row, "orientation_reversal_censored")
            if not conditioned:
                status = "not_conditioned"
            elif success:
                status = "success"
            elif censored:
                status = "right_censored"
            else:
                status = row.get("orientation_reversal_status") or "right_censored"
            reversal_time = _number(row, "pre_reversal_time_s")
            observation_end = _number(row, "observation_end_time_s")
            follow_up = None
            if reversal_time is not None and observation_end is not None:
                # Saved prestate is one output interval before the actual switch.
                follow_up = observation_end - (reversal_time + 2.0)
            reversal.append(
                {
                    "task_id": row["task_id"],
                    "input_peak": _number(row, "input_peak"),
                    "length_um": length,
                    "orientation_reversal_status": status,
                    "orientation_reversal_latency_s": _number(
                        row, "orientation_reversal_latency_s"
                    ),
                    "post_reversal_follow_up_s": follow_up,
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
                }
            )
    dual.sort(key=lambda row: (row["assay"], row["length_um"]))
    reversal.sort(key=lambda row: (row["input_peak"], row["length_um"]))
    return dual, reversal


def _format_lengths(lengths: list[float]) -> list[str]:
    return [f"{value:g}" if float(value).is_integer() else f"{value:.1f}" for value in lengths]


def _status_matrix(
    rows: list[dict], groups: list, lengths: list[float], group_key: str,
    status_key: str, order: tuple[str, ...],
) -> np.ndarray:
    matrix = np.zeros((len(groups), len(lengths)), dtype=int)
    for row in rows:
        matrix[groups.index(row[group_key]), lengths.index(row["length_um"])] = order.index(
            row[status_key]
        )
    return matrix


def _plot_story(dual: list[dict], reversal: list[dict], destination: Path) -> None:
    lengths = sorted({row["length_um"] for row in dual})
    assays = list(WIDTH_ASSAYS)
    inputs = sorted({row["input_peak"] for row in reversal})
    figure, axes = plt.subplots(
        2, 2, figsize=(14, 7.6), gridspec_kw={"height_ratios": [0.42, 1.0]}
    )

    dual_status = _status_matrix(
        dual, assays, lengths, "assay", "competition_status", COMPETITION_ORDER
    )
    axes[0, 0].imshow(
        dual_status,
        aspect="auto",
        cmap=ListedColormap(COMPETITION_COLORS),
        norm=BoundaryNorm(np.arange(-0.5, 3.5), 3),
    )
    axes[0, 0].set_yticks(range(len(assays)), [WIDTH_ASSAYS[value] for value in assays])
    axes[0, 0].set_xticks(range(len(lengths)), _format_lengths(lengths), rotation=45)
    axes[0, 0].set_title("A  Intact two-front competition status", loc="left", fontweight="bold")

    for assay in assays:
        selected = [row for row in dual if row["assay"] == assay]
        resolved = [row for row in selected if row["competition_status"] == "resolved"]
        censored = [row for row in selected if row["competition_status"] == "right_censored"]
        color = WIDTH_COLORS[assay]
        axes[1, 0].plot(
            [row["length_um"] for row in resolved],
            [row["competition_resolution_latency_s"] for row in resolved],
            "o-", color=color, lw=2, label=WIDTH_ASSAYS[assay],
        )
        axes[1, 0].scatter(
            [row["length_um"] for row in censored],
            [row["resolution_censor_limit_s"] for row in censored],
            marker="v", facecolors="none", edgecolors=color, s=60,
        )
    axes[1, 0].set_xlabel("Cell length (µm)")
    axes[1, 0].set_ylabel("Two-front resolution latency (s)")
    axes[1, 0].grid(alpha=0.2)
    axes[1, 0].legend(frameon=False, title="End-cue width")

    reversal_status = _status_matrix(
        reversal, inputs, lengths, "input_peak", "orientation_reversal_status", REVERSAL_ORDER
    )
    axes[0, 1].imshow(
        reversal_status,
        aspect="auto",
        cmap=ListedColormap(REVERSAL_COLORS),
        norm=BoundaryNorm(np.arange(-0.5, 3.5), 3),
    )
    axes[0, 1].set_yticks(range(len(inputs)), [f"input {value:g}" for value in inputs])
    axes[0, 1].set_xticks(range(len(lengths)), _format_lengths(lengths), rotation=45)
    axes[0, 1].set_title("B  Intact side-orientation reversal status", loc="left", fontweight="bold")

    for input_peak in inputs:
        selected = [row for row in reversal if row["input_peak"] == input_peak]
        success = [row for row in selected if row["orientation_reversal_status"] == "success"]
        censored = [row for row in selected if row["orientation_reversal_status"] == "right_censored"]
        axes[1, 1].plot(
            [row["length_um"] for row in success],
            [row["orientation_reversal_latency_s"] for row in success],
            "o-", color=INPUT_COLORS[input_peak], lw=2 if input_peak > 0.2 else 1.3,
            alpha=0.75 if input_peak == 0.2 else 1.0,
            label=f"input {input_peak:g}" + (" (near threshold)" if input_peak == 0.2 else ""),
        )
        axes[1, 1].scatter(
            [row["length_um"] for row in censored],
            [row["post_reversal_follow_up_s"] for row in censored],
            marker="v", facecolors="none", edgecolors=INPUT_COLORS[input_peak], s=60,
        )
    axes[1, 1].set_xlabel("Cell length (µm)")
    axes[1, 1].set_ylabel("Reversal latency after cue switch (s)")
    axes[1, 1].grid(alpha=0.2)
    axes[1, 1].legend(frameon=False)

    figure.legend(
        handles=[
            Patch(color=color, label=COMPETITION_LABELS[status])
            for status, color in zip(COMPETITION_ORDER, COMPETITION_COLORS)
        ] + [Line2D([], [], marker="v", markerfacecolor="none", color="black", ls="", label="Censored lower bound")],
        loc="lower center", ncol=4, frameon=False,
    )
    figure.suptitle(
        "Intact De Belly model: strong competition slowdown; input-dependent reversal scaling\n"
        "One deterministic run per condition; reversal uses length-scaled conditioning and follow-up"
    )
    figure.tight_layout(rect=(0, 0.08, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_two_site(dual: list[dict], destination: Path) -> None:
    lengths = sorted({row["length_um"] for row in dual})
    assays = list(WIDTH_ASSAYS)
    figure, axes = plt.subplots(2, 1, figsize=(11, 7.2), gridspec_kw={"height_ratios": [0.55, 1]})
    status = _status_matrix(
        dual, assays, lengths, "assay", "competition_status", COMPETITION_ORDER
    )
    axes[0].imshow(
        status, aspect="auto", cmap=ListedColormap(COMPETITION_COLORS),
        norm=BoundaryNorm(np.arange(-0.5, 3.5), 3),
    )
    axes[0].set_yticks(range(2), [WIDTH_ASSAYS[value] for value in assays])
    axes[0].set_xticks(range(len(lengths)), _format_lengths(lengths))
    axes[0].set_title("Competition event sequence")
    endpoint_symbols = {"unpolarized": "U", "single_front": "S", "multifront": "M", "polarized_no_resolved_peak": "P"}
    for i, assay in enumerate(assays):
        for row in [item for item in dual if item["assay"] == assay]:
            axes[0].text(
                lengths.index(row["length_um"]), i, endpoint_symbols[row["endpoint_state"]],
                ha="center", va="center", fontweight="bold", fontsize=8,
            )
    for assay in assays:
        selected = [row for row in dual if row["assay"] == assay]
        color = WIDTH_COLORS[assay]
        resolved = [row for row in selected if row["competition_status"] == "resolved"]
        censored = [row for row in selected if row["competition_status"] == "right_censored"]
        axes[1].plot(
            [row["length_um"] for row in resolved],
            [row["competition_resolution_latency_s"] for row in resolved],
            "o-", color=color, lw=2.2, label=WIDTH_ASSAYS[assay],
        )
        axes[1].scatter(
            [row["length_um"] for row in censored],
            [row["resolution_censor_limit_s"] for row in censored],
            marker="v", facecolors="none", edgecolors=color, s=65,
        )
    axes[1].set_xlabel("Cell length (µm)")
    axes[1].set_ylabel("Resolution latency (s)")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False, title="End-cue width")
    figure.suptitle(
        "Intact De Belly two-site challenge: length effect and cue-width sensitivity\n"
        "Triangles give the observed lower bound when resolution remains censored"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_reversal(reversal: list[dict], destination: Path) -> None:
    inputs = sorted({row["input_peak"] for row in reversal})
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for input_peak in inputs:
        selected = [row for row in reversal if row["input_peak"] == input_peak]
        success = [row for row in selected if row["orientation_reversal_status"] == "success"]
        lengths = np.asarray([row["length_um"] for row in success])
        latency = np.asarray([row["orientation_reversal_latency_s"] for row in success])
        color = INPUT_COLORS[input_peak]
        label = f"input {input_peak:g}" + (" (near threshold)" if input_peak == 0.2 else "")
        axes[0].plot(lengths, latency, "o-", color=color, lw=2, label=label)
        axes[1].plot(lengths, latency / lengths**2, "o-", color=color, lw=2, label=label)
        axes[2].plot(
            [row["length_um"] for row in selected],
            [row["pre_reversal_domain_count"] for row in selected],
            "o-", color=color, lw=2, label=label,
        )
    axes[0].set_ylabel("Reversal latency (s)")
    axes[1].set_ylabel("Reversal latency / L² (s/µm²)")
    axes[2].set_ylabel("Detected Rac maxima before reversal")
    for axis in axes:
        axis.set_xlabel("Cell length (µm)")
        axis.grid(alpha=0.2)
    axes[0].set_title("A  Absolute switching time", loc="left", fontweight="bold")
    axes[1].set_title("B  Diffusive scaling diagnostic", loc="left", fontweight="bold")
    axes[2].set_title("C  Conditioning morphology", loc="left", fontweight="bold")
    axes[0].legend(frameon=False)
    figure.suptitle(
        "Intact De Belly reversal after length-scaled conditioning\n"
        "Primary endpoint: sustained switch in Rac side orientation"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.88))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_arm_context(old_rows: list[dict[str, str]], destination: Path) -> None:
    variants = (
        "rac_to_rho_half", "rho_to_rac_half", "rac_to_rho_off",
        "rho_to_rac_off", "both_off",
    )
    labels = {
        "rac_to_rho_half": "Rac→Rho 50%", "rho_to_rac_half": "Rho→Rac 50%",
        "rac_to_rho_off": "Rac→Rho off", "rho_to_rac_off": "Rho→Rac off",
        "both_off": "Both off",
    }
    colors = dict(zip(variants, plt.cm.cividis(np.linspace(0.05, 0.90, len(variants)))))
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.3), sharey=False)
    for axis, input_peak in zip(axes, (0.2, 0.5, 1.0)):
        selected = [row for row in old_rows if np.isclose(_number(row, "input_peak"), input_peak)]
        intact = {
            _number(row, "length_um"): _number(row, "orientation_reversal_latency_s")
            for row in selected if row["coupling_variant"] == "intact"
        }
        for variant in variants:
            rows = [row for row in selected if row["coupling_variant"] == variant]
            x, delta = [], []
            for row in rows:
                length = _number(row, "length_um")
                latency = _number(row, "orientation_reversal_latency_s")
                if latency is not None and intact.get(length) is not None:
                    x.append(length)
                    delta.append(latency - intact[length])
            axis.plot(x, delta, "o-", lw=1.3, ms=4, color=colors[variant], label=labels[variant])
        axis.axhline(0, color="black", lw=0.8)
        axis.set_title(f"Input {input_peak:g}")
        axis.set_xlabel("Cell length (µm)")
        axis.set_ylabel("Latency change from intact (s)")
        axis.grid(alpha=0.2)
    handles, labels_out = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels_out, loc="lower center", ncol=5, frameon=False)
    figure.suptitle(
        "Secondary mechanistic context: arm perturbations relative to intact\n"
        "Original titration with fixed 250-s conditioning; deterministic differences"
    )
    figure.tight_layout(rect=(0, 0.12, 1, 0.87))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_cross_model(
    dual: list[dict], reversal: list[dict], phase1_root: Path, destination: Path
) -> None:
    competing = _read(phase1_root / "jilkine_competing_fronts.csv")
    jreversal = _read(phase1_root / "jilkine_gradient_reversal.csv")
    model_labels = {
        "jilkine_wave_pinning": "Wave-pinning", "jilkine_otsuji": "Otsuji",
        "jilkine_goryachev": "Goryachev", "jilkine_legi": "LEGI",
    }
    model_colors = dict(zip(model_labels, plt.cm.Set2(np.linspace(0.05, 0.85, 4))))
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 8))
    for model, label in model_labels.items():
        rows = [row for row in competing if row["model"] == model]
        resolved = [row for row in rows if _number(row, "resolution_latency_median_s") is not None]
        axes[0, 0].plot(
            [_number(row, "length_um") for row in resolved],
            [_number(row, "resolution_latency_median_s") for row in resolved],
            "o-", color=model_colors[model], ms=3.5, label=label,
        )
    axes[0, 0].set_title("Jilkine-reviewed models")
    axes[0, 0].set_ylabel("Two-front resolution latency (s)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    for assay in WIDTH_ASSAYS:
        rows = [row for row in dual if row["assay"] == assay and row["competition_status"] == "resolved"]
        axes[0, 1].plot(
            [row["length_um"] for row in rows],
            [row["competition_resolution_latency_s"] for row in rows],
            "o-", color=WIDTH_COLORS[assay], label=WIDTH_ASSAYS[assay],
        )
    axes[0, 1].set_title("Intact De Belly")
    axes[0, 1].legend(frameon=False, fontsize=8)
    for model, label in model_labels.items():
        rows = [row for row in jreversal if row["model"] == model and row["reversal_status"] == "success"]
        axes[1, 0].plot(
            [_number(row, "length_um") for row in rows],
            [_number(row, "reversal_latency_s") for row in rows],
            "o-", color=model_colors[model], ms=3.5, label=label,
        )
    axes[1, 0].set_ylabel("Reversal latency (s)")
    axes[1, 0].set_title("Jilkine-reviewed models")
    for input_peak in (0.2, 0.5, 1.0):
        rows = [
            row for row in reversal
            if row["input_peak"] == input_peak and row["orientation_reversal_status"] == "success"
        ]
        axes[1, 1].plot(
            [row["length_um"] for row in rows],
            [row["orientation_reversal_latency_s"] for row in rows],
            "o-", color=INPUT_COLORS[input_peak], label=f"input {input_peak:g}",
        )
    axes[1, 1].set_title("Intact De Belly")
    axes[1, 1].legend(frameon=False, fontsize=8)
    for axis in axes.flat:
        axis.set_xlabel("Cell length (µm)")
        axis.grid(alpha=0.2)
    figure.suptitle(
        "Cross-model synthesis: robust competition slowdown and qualified reversal scaling\n"
        "Absolute timescales and protocols differ; panels intentionally use independent y-scales"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _summary_rows(dual: list[dict], reversal: list[dict]) -> list[dict]:
    output = []
    for assay, label in WIDTH_ASSAYS.items():
        rows = [row for row in dual if row["assay"] == assay]
        resolved = [row for row in rows if row["competition_status"] == "resolved"]
        rho = None
        if len(resolved) >= 2:
            rho = float(spearmanr(
                [row["length_um"] for row in resolved],
                [row["competition_resolution_latency_s"] for row in resolved],
            ).statistic)
        output.append(
            {
                "assay": "two_site_competition",
                "condition": label,
                "condition_count": len(rows),
                "success_count": len(resolved),
                "censored_count": sum(row["competition_status"] == "right_censored" for row in rows),
                "spearman_rho_length_vs_latency": rho,
                "latency_at_smallest_resolved_length_s": resolved[0]["competition_resolution_latency_s"] if resolved else None,
                "latency_at_largest_resolved_length_s": resolved[-1]["competition_resolution_latency_s"] if resolved else None,
            }
        )
    for input_peak in sorted({row["input_peak"] for row in reversal}):
        rows = [row for row in reversal if row["input_peak"] == input_peak]
        success = [row for row in rows if row["orientation_reversal_status"] == "success"]
        rho = None
        if len(success) >= 2:
            rho = float(spearmanr(
                [row["length_um"] for row in success],
                [row["orientation_reversal_latency_s"] for row in success],
            ).statistic)
        output.append(
            {
                "assay": "gradient_reversal",
                "condition": f"input {input_peak:g}",
                "condition_count": len(rows),
                "success_count": len(success),
                "censored_count": sum(row["orientation_reversal_status"] == "right_censored" for row in rows),
                "spearman_rho_length_vs_latency": rho,
                "latency_at_smallest_resolved_length_s": success[0]["orientation_reversal_latency_s"] if success else None,
                "latency_at_largest_resolved_length_s": success[-1]["orientation_reversal_latency_s"] if success else None,
            }
        )
    return output


def _report(output_root: Path, collection: dict, summary: list[dict]) -> None:
    rows = {row["condition"]: row for row in summary}
    text = [
        "# Intact De Belly length-effects campaign",
        "",
        "## Design",
        "",
        f"All {collection['success_count']}/{collection['task_count']} tasks completed successfully.",
        "The campaign contains only the intact model. It uses 14 lengths from 12–64 µm,",
        "pairs the original two-site cue width (σ=0.08L) with a fixed σ=1.6 µm cue,",
        "extends two-front follow-up through at least 2,000 s at L=50 and more than",
        "3,000 s at L=64, and scales reversal conditioning/follow-up with L² above L=20.",
        "Each condition is one deterministic run.",
        "",
        "## Results",
        "",
        "The clearest result is the intact two-front competition assay. With the original",
        "fractional-width cue, latency rose monotonically from 230 s at L=12 to 1,200 s",
        "at L=60; L=64 remained unresolved after 3,219 s of post-establishment observation.",
        "With fixed σ=1.6 µm, every run resolved and latency rose monotonically from 240 s",
        "at L=12 to 2,956 s at L=64. Thus the length effect survives removal of the",
        "fractional-cue-width confound; the fractional L=64 case is slower still and censored.",
        "",
        "All 42 length-scaled reversal runs made a sustained side switch. The length effect",
        "depends on input: input 0.5 shows the clearest ordering (52 s at L=12 and 88 s",
        "at L=64, ρ=0.911); input 1 increases through intermediate lengths and then",
        "plateaus near 60–62 s (ρ=0.460); input 0.2 is U-shaped (ρ=0.187). Reversal",
        "therefore supports a qualified, regime-dependent size effect rather than the",
        "near-universal competition result.",
        "",
    ]
    for label in ("σ = 0.08L", "σ = 1.6 µm", "input 0.2", "input 0.5", "input 1"):
        row = rows[label]
        text.append(
            f"- {label}: {row['success_count']}/{row['condition_count']} resolved/successful; "
            f"{row['censored_count']} censored; Spearman ρ(length, latency) = "
            f"{row['spearman_rho_length_vs_latency']:.3f}."
        )
    text.extend(
        [
            "",
            "Spearman values describe ordering across a deterministic design grid; they are not",
            "replicate-based inferential statistics.",
            "",
            "## Figures",
            "",
            "![Primary intact-model result](debelly_intact_length_story.png)",
            "",
            "- `debelly_intact_length_story.png`: primary intact-model two-site and reversal result.",
            "- `debelly_intact_two_site_width_sensitivity.png`: extended competition and cue-width control.",
            "- `debelly_intact_reversal_length_scaled.png`: raw reversal latency, L² diagnostic, and prestate morphology.",
            "- `jilkine_debelly_length_synthesis.png`: Jilkine/deBelly cross-model synthesis.",
            "- `debelly_arm_perturbations_secondary.png`: arm perturbations retained as secondary context.",
            "",
            "## Interpretation",
            "",
            "> In the intact De Belly model, larger domains resolve competing fronts much",
            "> more slowly. Reversal also slows with length in the intermediate-input",
            "> regime, but this effect saturates at strong input and becomes non-monotone",
            "> near the conditioning/morphology transition.",
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
    parser.add_argument(
        "--titration-root", type=Path,
        default=ROOT / "results" / "campaigns" / "debelly_input_titration_v1",
    )
    arguments = parser.parse_args(argv)
    collection = json.loads((arguments.output_root / "collection_summary.json").read_text())
    if not collection.get("complete"):
        raise SystemExit("campaign collection is incomplete")
    dual, reversal = _clean(_read(arguments.output_root / "task_summary.csv"))
    if len(dual) != 28 or len(reversal) != 42:
        raise SystemExit(f"unexpected task split: {len(dual)} dual, {len(reversal)} reversal")
    _write(arguments.output_root / "debelly_intact_dual_site.csv", dual)
    _write(arguments.output_root / "debelly_intact_reversal.csv", reversal)
    summary = _summary_rows(dual, reversal)
    _write(arguments.output_root / "length_effect_summary.csv", summary)
    _plot_story(dual, reversal, arguments.output_root / "debelly_intact_length_story.png")
    _plot_two_site(dual, arguments.output_root / "debelly_intact_two_site_width_sensitivity.png")
    _plot_reversal(reversal, arguments.output_root / "debelly_intact_reversal_length_scaled.png")
    _plot_arm_context(
        _read(arguments.titration_root / "reversal_titration.csv"),
        arguments.output_root / "debelly_arm_perturbations_secondary.png",
    )
    _plot_cross_model(
        dual, reversal, arguments.phase1_root,
        arguments.output_root / "jilkine_debelly_length_synthesis.png",
    )
    _report(arguments.output_root, collection, summary)
    print(arguments.output_root)


if __name__ == "__main__":
    main()
