#!/usr/bin/env python3
"""Create a self-explaining Phase 1 campaign report bundle.

The figures deliberately keep the four assay types separate. Summary panels
are paired with representative saved trajectories so categorical phenotype
labels can be checked against the underlying space--time behavior.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import sys
from typing import Callable

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.contract import Trajectory  # noqa: E402


MODEL_ORDER = (
    "jilkine_wave_pinning",
    "jilkine_otsuji",
    "jilkine_goryachev",
    "jilkine_legi",
)
MODEL_LABELS = {
    "jilkine_wave_pinning": "Wave-pinning",
    "jilkine_otsuji": "Otsuji",
    "jilkine_goryachev": "Goryachev",
    "jilkine_legi": "LEGI",
}
VARIANT_ORDER = (
    "intact",
    "rac_to_rho_half",
    "rho_to_rac_half",
    "rac_to_rho_off",
    "rho_to_rac_off",
    "both_off",
)
VARIANT_LABELS = {
    "intact": "Both arms intact",
    "rac_to_rho_half": "Rac→Rho arm 50%",
    "rho_to_rac_half": "Rho→Rac arm 50%",
    "rac_to_rho_off": "Rac→Rho arm off",
    "rho_to_rac_off": "Rho→Rac arm off",
    "both_off": "Both arms off",
}
VARIANT_COLORS = dict(
    zip(VARIANT_ORDER, plt.cm.viridis(np.linspace(0.05, 0.95, len(VARIANT_ORDER))))
)
COMPETITION_STATUS = ("not_established", "resolved", "right_censored")
COMPETITION_LABELS = {
    "not_established": "Two fronts never sustained",
    "resolved": "Two fronts → one sustained front",
    "right_censored": "No sustained one-front resolution before run end",
}
COMPETITION_SHORT_LABELS = {
    "not_established": "No two-front phase",
    "resolved": "Resolved to one front",
    "right_censored": "No one-front resolution",
}
COMPETITION_COLORS = ("#d8d8d8", "#2878b5", "#e68a2e")
REVERSAL_STATUS = ("not_conditioned", "success", "right_censored")
REVERSAL_LABELS = {
    "not_conditioned": "No qualifying pre-reversal front",
    "success": "Opposite front sustained",
    "right_censored": "No opposite front by run end",
}
REVERSAL_COLORS = ("#d8d8d8", "#2a9d5b", "#d95f45")


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    if value in {"", "None", "null", "not_applicable"}:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if np.isfinite(number) else None


def _boolean(row: dict[str, str], key: str) -> bool | None:
    value = row.get(key, "").strip().lower()
    if value == "true":
        return True
    if value == "false":
        return False
    return None


def _write(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _manifest_index(manifest: dict) -> dict[str, dict]:
    return {task["task_id"]: task for task in manifest["tasks"]}


def _production(
    rows: list[dict[str, str]],
    *,
    model: str | None = None,
    assay: str | None = None,
) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row.get("tier") == "production"
        and row.get("status") == "success"
        and (model is None or row.get("model") == model)
        and (assay is None or row.get("assay") == assay)
    ]


def _jilkine_competition(rows: list[dict[str, str]]) -> list[dict]:
    result = []
    selected = _production(rows, assay="competing_fronts")
    for model in MODEL_ORDER:
        lengths = sorted(
            {_number(row, "length_um") for row in selected if row.get("model") == model}
            - {None}
        )
        for length in lengths:
            members = [
                row
                for row in selected
                if row.get("model") == model and _number(row, "length_um") == length
            ]
            statuses = [row.get("competition_status") for row in members]
            endpoints = [row.get("endpoint_state") for row in members]
            latencies = [
                value
                for row in members
                if (value := _number(row, "competition_resolution_latency_s"))
                is not None
            ]
            result.append(
                {
                    "model": model,
                    "length_um": length,
                    "n": len(members),
                    **{
                        f"competition_{status}_fraction": statuses.count(status)
                        / len(members)
                        for status in COMPETITION_STATUS
                    },
                    "competition_established_fraction": float(
                        np.mean([_boolean(row, "competition_established") for row in members])
                    ),
                    "final_single_front_fraction": endpoints.count("single_front")
                    / len(members),
                    "final_multifront_fraction": endpoints.count("multifront")
                    / len(members),
                    "final_unpolarized_fraction": endpoints.count("unpolarized")
                    / len(members),
                    "resolution_latency_median_s": (
                        "" if not latencies else float(np.median(latencies))
                    ),
                }
            )
    return result


def _reversal_rows(
    rows: list[dict[str, str]],
    task_index: dict[str, dict],
    *,
    models: tuple[str, ...] | None = None,
) -> list[dict]:
    selected = _production(rows, assay="gradient_reversal")
    if models is not None:
        selected = [row for row in selected if row.get("model") in models]
    result = []
    for row in selected:
        task = task_index[row["task_id"]]
        input_reversal = float(task["protocol"]["reversal_time"])
        end_time = float(task["simulation"]["end_time"])
        result.append(
            {
                "task_id": row["task_id"],
                "model": row["model"],
                "coupling_variant": row.get("coupling_variant", "not_applicable"),
                "length_um": _number(row, "length_um"),
                "input_reversal_time_s": input_reversal,
                "post_reversal_observation_s": end_time - input_reversal,
                "reversal_status": row.get("reversal_status"),
                "reversal_conditioned": _boolean(row, "reversal_conditioned"),
                "reversal_success": _boolean(row, "reversal_success"),
                "reversal_latency_s": _number(row, "reversal_latency_s"),
                "old_front_decay_latency_s": _number(
                    row, "old_front_decay_latency_s"
                ),
                "endpoint_state": row.get("endpoint_state"),
            }
        )
    return sorted(
        result,
        key=lambda item: (
            item["model"],
            item["coupling_variant"],
            item["length_um"],
        ),
    )


def _tracking_rows(
    rows: list[dict[str, str]],
    *,
    models: tuple[str, ...] | None = None,
) -> list[dict]:
    selected = _production(rows, assay="dynamic_tracking")
    if models is not None:
        selected = [row for row in selected if row.get("model") in models]
    result = []
    for row in selected:
        period = _number(row, "tracking_period_s")
        lag = _number(row, "phase_lag_seconds")
        result.append(
            {
                "task_id": row["task_id"],
                "model": row["model"],
                "coupling_variant": row.get("coupling_variant", "not_applicable"),
                "length_um": _number(row, "length_um"),
                "period_s": period,
                "response_amplitude": _number(row, "response_amplitude"),
                "phase_lag_seconds": lag,
                "phase_lag_fraction_of_period": (
                    None if lag is None or period is None else lag / period
                ),
                "tracking_correlation": _number(row, "tracking_correlation"),
                "tracking_reliable": _boolean(row, "tracking_reliable"),
            }
        )
    return sorted(
        result,
        key=lambda item: (
            item["model"],
            item["coupling_variant"],
            item["length_um"],
            item["period_s"],
        ),
    )


def _debelly_propagation(rows: list[dict[str, str]]) -> list[dict]:
    result = []
    for row in _production(
        rows,
        model="debelly_1d_mechanochemical",
        assay="single_site_propagation",
    ):
        result.append(
            {
                "task_id": row["task_id"],
                "coupling_variant": row["coupling_variant"],
                "length_um": _number(row, "length_um"),
                "rac_proximal_peak_change": _number(
                    row, "rac_proximal_peak_change"
                ),
                "rac_distal_peak_change": _number(row, "rac_distal_peak_change"),
                "rac_distal_peak_time_s": _number(row, "rac_distal_peak_time_s"),
                "rac_distal_to_proximal_abs_ratio": _number(
                    row, "rac_distal_to_proximal_abs_ratio"
                ),
                "rho_proximal_peak_change": _number(
                    row, "rho_proximal_peak_change"
                ),
                "rho_distal_peak_change": _number(row, "rho_distal_peak_change"),
                "rho_distal_peak_time_s": _number(row, "rho_distal_peak_time_s"),
                "maximum_absolute_tension": _number(
                    row, "maximum_absolute_tension"
                ),
            }
        )
    return sorted(
        result,
        key=lambda row: (
            VARIANT_ORDER.index(row["coupling_variant"]),
            row["length_um"],
        ),
    )


def _debelly_dual_site(rows: list[dict[str, str]]) -> list[dict]:
    result = []
    for row in _production(
        rows,
        model="debelly_1d_mechanochemical",
        assay="dual_site_challenge",
    ):
        result.append(
            {
                "task_id": row["task_id"],
                "coupling_variant": row["coupling_variant"],
                "length_um": _number(row, "length_um"),
                "competition_status": row.get("competition_status"),
                "competition_established": _boolean(
                    row, "competition_established"
                ),
                "competition_resolution_latency_s": _number(
                    row, "competition_resolution_latency_s"
                ),
                "endpoint_state": row.get("endpoint_state"),
                "final_spatial_contrast": _number(row, "final_spatial_contrast"),
            }
        )
    return sorted(
        result,
        key=lambda row: (
            VARIANT_ORDER.index(row["coupling_variant"]),
            row["length_um"],
        ),
    )


def _dependency_summary(rows: list[dict]) -> dict:
    intact = {
        row["length_um"]: row
        for row in rows
        if row["coupling_variant"] == "intact"
    }
    ablated = {
        row["length_um"]: row
        for row in rows
        if row["coupling_variant"] == "both_off"
    }
    paired = []
    for length in sorted(set(intact) & set(ablated)):
        intact_value = abs(intact[length]["rac_distal_peak_change"])
        ablated_value = abs(ablated[length]["rac_distal_peak_change"])
        paired.append(
            {
                "length_um": length,
                "intact_abs_distal_rac": intact_value,
                "both_off_abs_distal_rac": ablated_value,
                "absolute_coupling_effect": intact_value - ablated_value,
                "fraction_removed_by_ablation": 1.0
                - ablated_value / max(intact_value, 1e-12),
                "intact_to_both_off_fold_separation": intact_value
                / max(ablated_value, 1e-12),
            }
        )
    lengths = np.asarray([row["length_um"] for row in paired])
    effects = np.asarray([row["absolute_coupling_effect"] for row in paired])
    fractions = np.asarray([row["fraction_removed_by_ablation"] for row in paired])
    fold_separation = np.asarray(
        [row["intact_to_both_off_fold_separation"] for row in paired]
    )
    wgd_anchor = min(
        paired, key=lambda row: abs(row["length_um"] - 2 ** (1 / 3) * 20.0)
    )
    return {
        "pairs": paired,
        "absolute_effect_slope_per_um": float(np.polyfit(lengths, effects, 1)[0]),
        "fraction_removed_slope_per_um": float(
            np.polyfit(lengths, fractions, 1)[0]
        ),
        "log10_fold_separation_slope_per_um": float(
            np.polyfit(lengths, np.log10(fold_separation), 1)[0]
        ),
        "smallest_length": paired[0],
        "reference_length": next(
            row for row in paired if row["length_um"] == 20.0
        ),
        "wgd_linear_scale_anchor": wgd_anchor,
        "largest_length": paired[-1],
    }


def _trajectory_row(
    rows: list[dict[str, str]],
    output_root: Path,
    predicate: Callable[[dict[str, str]], bool],
) -> tuple[dict[str, str], Trajectory]:
    candidates = []
    for row in rows:
        path = output_root / "tasks" / row["task_id"] / "trajectory.npz"
        if predicate(row) and path.is_file():
            candidates.append((row, path))
    if len(candidates) != 1:
        raise RuntimeError(f"expected one saved trajectory, found {len(candidates)}")
    row, path = candidates[0]
    return row, Trajectory.load(path)


def _near(row: dict[str, str], key: str, value: float) -> bool:
    candidate = _number(row, key)
    return candidate is not None and np.isclose(candidate, value)


def _orientation(trajectory: Trajectory) -> np.ndarray:
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    left = np.mean(trajectory.observable[:, relative <= 0.2], axis=1)
    right = np.mean(trajectory.observable[:, relative >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _regional_change(
    trajectory: Trajectory, field: str, baseline_time: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    baseline_candidates = np.flatnonzero(trajectory.time < baseline_time)
    baseline = int(baseline_candidates[-1]) if baseline_candidates.size else 0
    change = trajectory.fields[field] - trajectory.fields[field][baseline]
    proximal = np.mean(change[:, relative <= 0.2], axis=1)
    distal = np.mean(change[:, relative >= 0.8], axis=1)
    return change, proximal, distal


def _plot_protocols(destination: Path) -> None:
    figure, axes = plt.subplots(1, 4, figsize=(13.5, 3.4), sharey=True)
    x = np.linspace(0, 1, 300)
    left = np.zeros_like(x)
    selected = x <= 0.1
    left[selected] = 0.5 * (1 + np.cos(np.pi * x[selected] / 0.1))
    right = left[::-1]
    axes[0].plot(x, 1.01 * left + 0.99 * right, color="#7030a0", lw=2.5)
    axes[0].annotate("left 1.01×", (0.03, 1.01), xytext=(0.14, 0.82),
                     arrowprops={"arrowstyle": "-", "color": "0.4"})
    axes[0].annotate("right 0.99×", (0.97, 0.99), xytext=(0.54, 0.67),
                     arrowprops={"arrowstyle": "-", "color": "0.4"})
    axes[0].set_title("Two-site challenge")
    axes[0].text(0.5, 0.12, "Transient end cues", ha="center", fontsize=9)
    gaussian = np.exp(-0.5 * (x / 0.2) ** 2)
    axes[1].plot(x, gaussian, color="#c44e52", lw=2.5)
    axes[1].axvspan(0, 0.2, color="#c44e52", alpha=0.12)
    axes[1].set_title("Single-site propagation")
    axes[1].text(0.5, 0.12, "De Belly only", ha="center", fontsize=9)
    axes[2].plot(x, 1 - x, lw=2.5, label="before reversal")
    axes[2].plot(x, x, lw=2.5, label="after reversal")
    axes[2].set_title("Gradient reversal")
    axes[2].legend(frameon=False, fontsize=8, loc="lower center")
    for phase, color, label in (
        (1.0, "#2878b5", "phase 0"),
        (0.0, "#777777", "quarter-cycle"),
        (-1.0, "#e68a2e", "half-cycle"),
    ):
        axes[3].plot(x, 0.5 * (1 - phase * (2 * x - 1)), lw=2.3,
                     color=color, label=label)
    axes[3].set_title("Alternating gradient")
    axes[3].legend(frameon=False, fontsize=8, loc="lower center")
    for axis in axes:
        axis.set_xlim(0, 1)
        axis.set_ylim(-0.03, 1.12)
        axis.set_xlabel("Relative position x/L")
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Input / peak input")
    figure.suptitle("The campaign contains four distinct stimulation assays")
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_jilkine_competition(rows: list[dict], destination: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for axis, model in zip(axes.flat, MODEL_ORDER):
        selected = [row for row in rows if row["model"] == model]
        x = np.asarray([row["length_um"] for row in selected])
        bottom = np.zeros(len(selected))
        for status, color in zip(COMPETITION_STATUS, COMPETITION_COLORS):
            values = np.asarray(
                [row[f"competition_{status}_fraction"] for row in selected]
            )
            axis.bar(x, values, width=1.55, bottom=bottom, color=color,
                     edgecolor="white", linewidth=0.4)
            bottom += values
        axis.set_title(f"{MODEL_LABELS[model]} (n={selected[0]['n']} per length)")
        axis.set_ylim(0, 1.04)
        axis.set_ylabel("Fraction of runs")
        axis.set_xlabel("Cell length (µm)")
        latency = np.asarray(
            [np.nan if row["resolution_latency_median_s"] == ""
             else row["resolution_latency_median_s"] for row in selected],
            dtype=float,
        )
        if np.any(np.isfinite(latency)):
            twin = axis.twinx()
            twin.plot(x, latency, "ko--", ms=3, lw=1)
            twin.set_ylabel("Resolution latency (s)", color="black")
            twin.tick_params(axis="y", labelsize=8)
    figure.legend(
        handles=[Patch(color=color, label=COMPETITION_LABELS[status])
                 for status, color in zip(COMPETITION_STATUS, COMPETITION_COLORS)]
        + [Line2D([], [], color="black", marker="o", ls="--",
                  label="Median latency among resolved runs")],
        loc="lower center", ncol=2, frameon=False,
    )
    figure.suptitle(
        "Jilkine-model assay 1: response to two transient end stimuli\n"
        "2% left/right cue difference; five seeded 0.01-amplitude perturbations except LEGI"
    )
    figure.tight_layout(rect=(0, 0.12, 1, 0.92))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_jilkine_competition_examples(
    rows: list[dict[str, str]], output_root: Path, destination: Path
) -> None:
    figure, axes = plt.subplots(1, 4, figsize=(13.5, 5.0), sharex=True, sharey=True)
    for axis, model in zip(axes, MODEL_ORDER):
        row, trajectory = _trajectory_row(
            rows, output_root,
            lambda candidate, model=model: candidate.get("model") == model
            and candidate.get("assay") == "competing_fronts"
            and _near(candidate, "length_um", 21.0)
            and candidate.get("replicate") == "0",
        )
        values = trajectory.observable
        low, high = float(np.min(values)), float(np.max(values))
        normalized = (values - low) / max(high - low, 1e-12)
        mesh = axis.pcolormesh(
            trajectory.x / float(trajectory.metadata["length_um"]), trajectory.time,
            normalized, shading="auto", cmap="magma", vmin=0, vmax=1,
        )
        axis.axhline(40, color="cyan", lw=1.2, ls="--")
        axis.set_title(
            f"{MODEL_LABELS[model]}\n"
            f"{COMPETITION_SHORT_LABELS[row['competition_status']]}\n"
            f"final {row['endpoint_state'].replace('_', ' ')}",
            fontsize=9,
        )
        axis.set_xlabel("x/L")
        axis.text(0.98, 0.03, f"raw range {low:.3g}–{high:.3g}",
                  transform=axis.transAxes, ha="right", va="bottom", fontsize=7,
                  color="white", bbox={"facecolor": "black", "alpha": 0.35, "pad": 2})
    axes[0].set_ylabel("Time (s)")
    colorbar = figure.colorbar(mesh, ax=axes, location="bottom", fraction=0.09, pad=0.13)
    colorbar.set_label("Observable normalized within each run (0=min, 1=max)")
    figure.suptitle(
        "Representative two-site trajectories at L=21 µm (seed 1)\n"
        "Cyan line: end of stimulus; time runs from bottom to top",
        y=0.99,
    )
    figure.subplots_adjust(left=0.06, right=0.98, top=0.69, bottom=0.22, wspace=0.16)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _status_code(status: str, order: tuple[str, ...]) -> int:
    return order.index(status)


def _plot_jilkine_reversal(rows: list[dict], destination: Path) -> None:
    figure, axes = plt.subplots(
        2, 4, figsize=(15, 6.5), sharex="col", gridspec_kw={"height_ratios": [0.6, 2.2]}
    )
    cmap = ListedColormap(REVERSAL_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, 3.5), cmap.N)
    for column, model in enumerate(MODEL_ORDER):
        selected = [row for row in rows if row["model"] == model]
        lengths = np.asarray([row["length_um"] for row in selected])
        codes = np.asarray([_status_code(row["reversal_status"], REVERSAL_STATUS)
                            for row in selected])[None, :]
        axes[0, column].imshow(codes, aspect="auto", cmap=cmap, norm=norm,
                               extent=(lengths[0] - 1, lengths[-1] + 1, 0, 1))
        axes[0, column].set_yticks([])
        axes[0, column].set_title(MODEL_LABELS[model])
        success = [row for row in selected if row["reversal_status"] == "success"]
        censored = [row for row in selected if row["reversal_status"] == "right_censored"]
        axes[1, column].plot([row["length_um"] for row in success],
                             [row["reversal_latency_s"] for row in success],
                             "o-", color="#2a9d5b")
        axes[1, column].scatter([row["length_um"] for row in censored],
                                [row["post_reversal_observation_s"] for row in censored],
                                marker="v", facecolors="none", edgecolors="#d95f45", s=45)
        axes[1, column].set_xlabel("Cell length (µm)")
        axes[1, column].set_ylabel("Seconds after cue reversal")
        axes[1, column].grid(alpha=0.2)
    figure.legend(
        handles=[Patch(color=color, label=REVERSAL_LABELS[status])
                 for status, color in zip(REVERSAL_STATUS, REVERSAL_COLORS)]
        + [Line2D([], [], color="#2a9d5b", marker="o", label="Latency if successful"),
           Line2D([], [], color="#d95f45", marker="v", markerfacecolor="none", ls="",
                  label="Follow-up limit if censored")],
        loc="lower center", ncol=3, frameon=False,
    )
    figure.suptitle(
        "Jilkine-model assay 2: reverse a fixed-peak linear gradient after conditioning\n"
        "One deterministic run per length; status strip above, latency/follow-up below"
    )
    figure.tight_layout(rect=(0, 0.16, 1, 0.90))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_reversal_examples(
    rows: list[dict[str, str]], task_index: dict[str, dict], output_root: Path,
    destination: Path, *, debelly: bool,
) -> None:
    if debelly:
        specifications = [("debelly_1d_mechanochemical", variant, VARIANT_LABELS[variant])
                          for variant in ("intact", "both_off")]
        length = 20.0
        title = "De Belly reversal examples at L=20 µm"
    else:
        specifications = [(model, "not_applicable", MODEL_LABELS[model]) for model in MODEL_ORDER]
        length = 21.0
        title = "Jilkine-model reversal examples at L=21 µm"
    figure, axes = plt.subplots(1, len(specifications), figsize=(4.1 * len(specifications), 3.5), sharey=True)
    axes = np.atleast_1d(axes)
    for axis, (model, variant, label) in zip(axes, specifications):
        row, trajectory = _trajectory_row(
            rows, output_root,
            lambda candidate, model=model, variant=variant: candidate.get("model") == model
            and candidate.get("assay") == "gradient_reversal"
            and _near(candidate, "length_um", length)
            and candidate.get("coupling_variant", "not_applicable") == variant,
        )
        task = task_index[row["task_id"]]
        reversal_time = float(task["protocol"]["reversal_time"])
        transition = float(task["protocol"]["transition"])
        cue = np.full_like(trajectory.time, -1.0)
        cue[trajectory.time >= reversal_time + transition] = 1.0
        in_transition = (trajectory.time >= reversal_time) & (trajectory.time < reversal_time + transition)
        cue[in_transition] = -1 + 2 * (trajectory.time[in_transition] - reversal_time) / transition
        axis.plot(trajectory.time, cue, "k--", lw=1.2, label="cue direction")
        axis.plot(trajectory.time, _orientation(trajectory), color="#4c72b0", lw=2,
                  label="model orientation")
        axis.axvspan(reversal_time, reversal_time + transition, color="gold", alpha=0.2)
        axis.axvline(reversal_time, color="0.4", lw=0.8)
        axis.set_title(f"{label}\n{REVERSAL_LABELS[row['reversal_status']]}")
        axis.set_xlabel("Time (s)")
        axis.set_ylim(-1.08, 1.08)
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Right-minus-left orientation")
    axes[0].legend(frameon=False, fontsize=8, loc="best")
    figure.suptitle(f"{title}\nYellow band: linear cue cross-fade")
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_tracking_grid(rows: list[dict], destination: Path, *, debelly: bool) -> None:
    if debelly:
        groups = [(variant, VARIANT_LABELS[variant]) for variant in VARIANT_ORDER]
        figure, axes = plt.subplots(2, 3, figsize=(12, 7.5), sharex=True, sharey=True)
        selection_key = "coupling_variant"
        title = ("De Belly assay 4: response to a periodically alternating gradient\n"
                 "Cell text: orientation amplitude / lag as fraction of period; one run per cell")
    else:
        groups = [(model, MODEL_LABELS[model]) for model in MODEL_ORDER]
        figure, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
        selection_key = "model"
        title = ("Jilkine-model assay 3: response to a periodically alternating gradient\n"
                 "Cell text: orientation amplitude / lag as fraction of period; one run per cell")
    amplitudes = [row["response_amplitude"] for row in rows]
    vmax = max(amplitudes)
    for axis, (group, label) in zip(axes.flat, groups):
        selected = [row for row in rows if row[selection_key] == group]
        lengths = sorted({row["length_um"] for row in selected})
        periods = sorted({row["period_s"] for row in selected})
        matrix = np.full((len(lengths), len(periods)), np.nan)
        reliable = np.zeros_like(matrix, dtype=bool)
        lag = np.full_like(matrix, np.nan)
        for row in selected:
            i = lengths.index(row["length_um"])
            j = periods.index(row["period_s"])
            matrix[i, j] = row["response_amplitude"]
            reliable[i, j] = row["tracking_reliable"]
            lag[i, j] = row["phase_lag_fraction_of_period"]
        mesh = axis.imshow(matrix, origin="lower", aspect="auto", cmap="viridis", vmin=0, vmax=vmax)
        for i in range(len(lengths)):
            for j in range(len(periods)):
                color = "white" if matrix[i, j] < 0.55 * vmax else "black"
                axis.text(j, i, f"{matrix[i,j]:.2f}\n{lag[i,j]:+.2f}P",
                          ha="center", va="center", fontsize=7.5, color=color)
                if not reliable[i, j]:
                    axis.plot(j, i, marker="x", color="#ff3b30", ms=13, mew=2)
        axis.set_title(label)
        axis.set_xticks(range(len(periods)), [f"{value:g}" for value in periods])
        axis.set_yticks(range(len(lengths)), [f"{value:g}" for value in lengths])
        axis.set_xlabel("Period (s)")
        axis.set_ylabel("Length (µm)")
    colorbar = figure.colorbar(mesh, ax=axes, location="bottom", fraction=0.05, pad=0.09)
    colorbar.set_label("Normalized orientation amplitude; red × fails declared tracking rule")
    figure.suptitle(title)
    figure.subplots_adjust(left=0.08, right=0.98, top=0.87, bottom=0.18, wspace=0.18, hspace=0.25)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_tracking_examples(
    rows: list[dict[str, str]], output_root: Path, destination: Path, *, debelly: bool
) -> None:
    if debelly:
        specifications = [("intact", 20.0, 100.0), ("intact", 64.0, 400.0),
                          ("both_off", 20.0, 100.0), ("both_off", 64.0, 400.0)]
        figure, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharey=True)
    else:
        specifications = [(model, 39.0, period) for model in MODEL_ORDER for period in (50.0, 400.0)]
        figure, axes = plt.subplots(4, 2, figsize=(11, 10), sharey=True)
    for axis, specification in zip(axes.flat, specifications):
        if debelly:
            variant, length, period = specification
            model = "debelly_1d_mechanochemical"
            label = VARIANT_LABELS[variant]
        else:
            model, length, period = specification
            variant = "not_applicable"
            label = MODEL_LABELS[model]
        row, trajectory = _trajectory_row(
            rows, output_root,
            lambda candidate, model=model, variant=variant, length=length, period=period: candidate.get("model") == model
            and candidate.get("assay") == "dynamic_tracking"
            and candidate.get("coupling_variant", "not_applicable") == variant
            and _near(candidate, "length_um", length)
            and _near(candidate, "tracking_period_s", period),
        )
        cue = -np.cos(2 * np.pi * trajectory.time / period)
        axis.plot(trajectory.time, cue, "k--", lw=1, label="cue direction")
        axis.plot(trajectory.time, _orientation(trajectory), color="#4c72b0", lw=1.8,
                  label="model orientation")
        axis.axvspan(0, period, color="0.5", alpha=0.08)
        axis.set_title(
            f"{label}; L={length:g} µm, P={period:g} s\n"
            f"A={float(row['response_amplitude']):.2f}, "
            f"lag={float(row['phase_lag_seconds']) / period:+.2f}P, "
            f"reliable={row['tracking_reliable']}"
        )
        axis.set_xlabel("Time (s)")
        axis.grid(alpha=0.15)
    for axis in axes[:, 0]:
        axis.set_ylabel("Right-minus-left orientation")
    axes.flat[0].legend(frameon=False, fontsize=8)
    title = "De Belly tracking examples" if debelly else "Jilkine-model tracking examples at L=39 µm"
    figure.suptitle(f"{title}\nGray band: first cycle excluded from metric")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_debelly_propagation(rows: list[dict], destination: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    metrics = (
        ("rac_distal_peak_change", "Signed distal ΔRac peak"),
        ("rho_distal_peak_change", "Signed distal ΔRho peak"),
        ("rac_distal_to_proximal_abs_ratio", "|Distal/proximal ΔRac peak|"),
        ("maximum_absolute_tension", "Maximum |tension| (pN)"),
    )
    for variant in VARIANT_ORDER:
        selected = [row for row in rows if row["coupling_variant"] == variant]
        for axis, (metric, ylabel) in zip(axes.flat, metrics):
            axis.plot([row["length_um"] for row in selected], [row[metric] for row in selected],
                      "o-", color=VARIANT_COLORS[variant], ms=4, label=VARIANT_LABELS[variant])
            axis.set_ylabel(ylabel)
            axis.set_xlabel("Cell length (µm)")
            axis.axhline(0, color="0.7", lw=0.7)
            axis.grid(alpha=0.15)
    figure.legend(*axes[0, 0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False)
    figure.suptitle(
        "De Belly assay 1: response after a left-edge Rac activation pulse\n"
        "Peak regional changes are relative to the pre-pulse frame; one run per point"
    )
    figure.tight_layout(rect=(0, 0.13, 1, 0.92))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_debelly_propagation_examples(
    rows: list[dict[str, str]], output_root: Path, destination: Path
) -> None:
    examples = []
    for variant in ("intact", "both_off"):
        row, trajectory = _trajectory_row(
            rows, output_root,
            lambda candidate, variant=variant: candidate.get("model") == "debelly_1d_mechanochemical"
            and candidate.get("assay") == "single_site_propagation"
            and candidate.get("coupling_variant") == variant
            and _near(candidate, "length_um", 20.0),
        )
        examples.append((variant, row, trajectory))
    changes = {(variant, field): _regional_change(trajectory, field, 50.0)
               for variant, _, trajectory in examples for field in ("rac", "rho")}
    limits = {field: max(np.max(np.abs(changes[(variant, field)][0]))
                         for variant in ("intact", "both_off")) for field in ("rac", "rho")}
    figure, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex="col", sharey="col")
    meshes = {}
    for row_index, (variant, _, trajectory) in enumerate(examples):
        for column, field in enumerate(("rac", "rho")):
            change, _, _ = changes[(variant, field)]
            limit = limits[field]
            meshes[field] = axes[row_index, column].pcolormesh(
                trajectory.x / 20.0, trajectory.time, change, shading="auto", cmap="coolwarm",
                norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
            )
            axes[row_index, column].axhline(50, color="lime", ls="--", lw=1)
            axes[row_index, column].axhline(150, color="lime", ls=":", lw=1)
            axes[row_index, column].set_title(f"{VARIANT_LABELS[variant]}: Δ{field.capitalize()}")
            axes[row_index, column].set_xlabel("x/L")
        _, proximal, distal = changes[(variant, "rac")]
        axes[row_index, 2].plot(trajectory.time, proximal, label="proximal 20%", lw=2)
        axes[row_index, 2].plot(trajectory.time, distal, label="distal 20%", lw=2)
        axes[row_index, 2].axvspan(50, 150, color="lime", alpha=0.12)
        axes[row_index, 2].axhline(0, color="0.6", lw=0.7)
        axes[row_index, 2].set_title(f"{VARIANT_LABELS[variant]}: regional ΔRac")
        axes[row_index, 2].set_xlabel("Time (s)")
        axes[row_index, 2].set_ylabel("Change from pre-pulse frame")
        axes[row_index, 2].legend(frameon=False, fontsize=8)
    axes[0, 0].set_ylabel("Time (s)")
    axes[1, 0].set_ylabel("Time (s)")
    figure.subplots_adjust(
        left=0.07, right=0.98, top=0.86, bottom=0.20, wspace=0.24, hspace=0.25
    )
    for field, position in (
        ("rac", [0.08, 0.07, 0.25, 0.025]),
        ("rho", [0.40, 0.07, 0.25, 0.025]),
    ):
        colorbar_axis = figure.add_axes(position)
        colorbar = figure.colorbar(
            meshes[field], cax=colorbar_axis, orientation="horizontal"
        )
        colorbar.set_label(f"Δ{field.capitalize()} (shared scale across rows)")
    figure.suptitle(
        "De Belly single-site examples at L=20 µm\n"
        "Green dashed/dotted lines: pulse on/off; distal means the far-right 20%"
    )
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_debelly_dual(rows: list[dict], destination: Path) -> None:
    lengths = sorted({row["length_um"] for row in rows})
    status = np.zeros((len(VARIANT_ORDER), len(lengths)), dtype=int)
    endpoints = np.empty(status.shape, dtype=object)
    latency = np.full(status.shape, np.nan)
    for row in rows:
        i = VARIANT_ORDER.index(row["coupling_variant"])
        j = lengths.index(row["length_um"])
        status[i, j] = _status_code(row["competition_status"], COMPETITION_STATUS)
        endpoints[i, j] = row["endpoint_state"]
        if row["competition_resolution_latency_s"] is not None:
            latency[i, j] = row["competition_resolution_latency_s"]
    cmap = ListedColormap(COMPETITION_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, 3.5), cmap.N)
    figure, axes = plt.subplots(
        2,
        1,
        figsize=(11, 8.4),
        gridspec_kw={"height_ratios": [1.4, 1]},
    )
    axes[0].imshow(status, aspect="auto", cmap=cmap, norm=norm)
    endpoint_letters = {"unpolarized": "U", "single_front": "S", "multifront": "M",
                        "polarized_no_resolved_peak": "P"}
    for i in range(status.shape[0]):
        for j in range(status.shape[1]):
            axes[0].text(j, i, endpoint_letters[endpoints[i, j]], ha="center", va="center",
                         fontsize=9, fontweight="bold")
    axes[0].set_xticks(range(len(lengths)), [f"{value:g}" for value in lengths])
    axes[0].set_yticks(range(len(VARIANT_ORDER)), [VARIANT_LABELS[value] for value in VARIANT_ORDER])
    axes[0].set_xlabel("Cell length (µm)")
    axes[0].set_title(
        "Cell background = event sequence; bold letter = final Rac state"
    )
    for i, variant in enumerate(VARIANT_ORDER):
        axes[1].plot(lengths, latency[i], "o-", color=VARIANT_COLORS[variant],
                     label=VARIANT_LABELS[variant])
    axes[1].set_ylabel("Resolution latency (s)")
    axes[1].set_xlabel("Cell length (µm)")
    axes[1].grid(alpha=0.2)
    resolved_variants = [
        variant
        for i, variant in enumerate(VARIANT_ORDER)
        if np.any(np.isfinite(latency[i]))
    ]
    axes[1].legend(
        handles=[
            Line2D(
                [],
                [],
                color=VARIANT_COLORS[variant],
                marker="o",
                label=VARIANT_LABELS[variant],
            )
            for variant in resolved_variants
        ],
        title="Lines: coupling state (resolved runs only)",
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    figure.legend(
        handles=[
            Patch(color=color, label=COMPETITION_LABELS[status])
            for status, color in zip(COMPETITION_STATUS, COMPETITION_COLORS)
        ],
        title="Cell background: competition event sequence",
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        bbox_to_anchor=(0.5, 0.070),
    )
    figure.text(
        0.5,
        0.025,
        "Final-state letters: U = unpolarized; S = single front; "
        "M = multiple fronts; P = polarized, no resolved peak.\n"
        "Orange is right-censored: two fronts occurred, but no sustained "
        "single-front state was observed by t = 700 s.",
        ha="center",
        va="center",
        fontsize=8,
    )
    figure.suptitle(
        "De Belly assay 2: response to two simultaneous end-localized Rac pulses\n"
        "2% left/right cue difference; one deterministic run per cell"
    )
    figure.tight_layout(rect=(0, 0.18, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_debelly_dual_examples(
    rows: list[dict[str, str]], output_root: Path, destination: Path
) -> None:
    examples = []
    for variant in ("intact", "both_off"):
        row, trajectory = _trajectory_row(
            rows, output_root,
            lambda candidate, variant=variant: candidate.get("model") == "debelly_1d_mechanochemical"
            and candidate.get("assay") == "dual_site_challenge"
            and candidate.get("coupling_variant") == variant
            and _near(candidate, "length_um", 20.0),
        )
        change, _, _ = _regional_change(trajectory, "rac", 50.0)
        examples.append((variant, row, trajectory, change))
    limit = max(np.max(np.abs(example[3])) for example in examples)
    figure, axes = plt.subplots(1, 2, figsize=(10, 4.3), sharex=True, sharey=True)
    for axis, (variant, row, trajectory, change) in zip(axes, examples):
        mesh = axis.pcolormesh(trajectory.x / 20.0, trajectory.time, change, shading="auto",
                               cmap="coolwarm", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
        axis.axhline(50, color="lime", ls="--", lw=1)
        axis.axhline(150, color="lime", ls=":", lw=1)
        axis.set_title(
            f"{VARIANT_LABELS[variant]}\n"
            f"{COMPETITION_SHORT_LABELS[row['competition_status']]}; "
            f"final {row['endpoint_state'].replace('_', ' ')}"
        )
        axis.set_xlabel("x/L")
    axes[0].set_ylabel("Time (s)")
    colorbar = figure.colorbar(mesh, ax=axes, location="bottom", fraction=0.08, pad=0.16)
    colorbar.set_label("ΔRac from pre-pulse frame (shared scale)")
    figure.suptitle("De Belly two-site examples at L=20 µm")
    figure.subplots_adjust(left=0.08, right=0.98, top=0.82, bottom=0.25, wspace=0.15)
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_debelly_reversal(rows: list[dict], destination: Path) -> None:
    lengths = sorted({row["length_um"] for row in rows})
    status = np.zeros((len(VARIANT_ORDER), len(lengths)), dtype=int)
    latency = np.full(status.shape, np.nan)
    censor = np.full(status.shape, np.nan)
    for row in rows:
        i = VARIANT_ORDER.index(row["coupling_variant"])
        j = lengths.index(row["length_um"])
        status[i, j] = _status_code(row["reversal_status"], REVERSAL_STATUS)
        if row["reversal_latency_s"] is not None:
            latency[i, j] = row["reversal_latency_s"]
        if row["reversal_status"] == "right_censored":
            censor[i, j] = row["post_reversal_observation_s"]
    cmap = ListedColormap(REVERSAL_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, 3.5), cmap.N)
    figure, axes = plt.subplots(2, 1, figsize=(11, 7), gridspec_kw={"height_ratios": [1.3, 1]})
    axes[0].imshow(status, aspect="auto", cmap=cmap, norm=norm)
    axes[0].set_xticks(range(len(lengths)), [f"{value:g}" for value in lengths])
    axes[0].set_yticks(range(len(VARIANT_ORDER)), [VARIANT_LABELS[value] for value in VARIANT_ORDER])
    axes[0].set_xlabel("Cell length (µm)")
    axes[0].set_title("Reversal classification")
    for i, variant in enumerate(VARIANT_ORDER):
        axes[1].plot(lengths, latency[i], "o-", color=VARIANT_COLORS[variant])
        valid = np.isfinite(censor[i])
        axes[1].scatter(np.asarray(lengths)[valid], censor[i, valid], marker="v",
                        facecolors="none", edgecolors=[VARIANT_COLORS[variant]], s=45)
    axes[1].set_ylabel("Seconds after cue reversal")
    axes[1].set_xlabel("Cell length (µm)")
    axes[1].grid(alpha=0.2)
    figure.legend(handles=[Patch(color=color, label=REVERSAL_LABELS[status])
                           for status, color in zip(REVERSAL_STATUS, REVERSAL_COLORS)]
                  + [Line2D([], [], color="black", marker="v", markerfacecolor="none", ls="",
                            label="Censoring follow-up limit")],
                  loc="lower center", ncol=2, frameon=False)
    figure.suptitle(
        "De Belly assay 3: reverse a linear Rac-activation gradient at t=250 s\n"
        "One deterministic run per cell"
    )
    figure.tight_layout(rect=(0, 0.12, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _assay_counts(rows: list[dict[str, str]]) -> list[dict]:
    groups: dict[tuple[str, str, str], int] = {}
    for row in rows:
        key = (row.get("tier", ""), row.get("model", ""), row.get("assay", ""))
        groups[key] = groups.get(key, 0) + 1
    return [{"tier": key[0], "model": key[1], "assay": key[2], "task_count": count}
            for key, count in sorted(groups.items())]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--manifest", type=Path,
                        default=ROOT / "campaigns" / "phase1_production_v1.json")
    arguments = parser.parse_args()
    task_table = arguments.output_root / "task_summary.csv"
    if not task_table.is_file():
        raise SystemExit(f"collect the campaign first: missing {task_table}")
    manifest = json.loads(arguments.manifest.read_text())
    rows = _read(task_table)
    task_index = _manifest_index(manifest)

    jilkine_competition = _jilkine_competition(rows)
    all_reversal = _reversal_rows(rows, task_index)
    jilkine_reversal = [row for row in all_reversal if row["model"] in MODEL_ORDER]
    debelly_reversal = [row for row in all_reversal if row["model"] == "debelly_1d_mechanochemical"]
    all_tracking = _tracking_rows(rows)
    jilkine_tracking = [row for row in all_tracking if row["model"] in MODEL_ORDER]
    debelly_tracking = [row for row in all_tracking if row["model"] == "debelly_1d_mechanochemical"]
    debelly_propagation = _debelly_propagation(rows)
    debelly_dual = _debelly_dual_site(rows)
    dependency = _dependency_summary(debelly_propagation)

    outputs = {
        "campaign_assay_counts.csv": _assay_counts(rows),
        "jilkine_competing_fronts.csv": jilkine_competition,
        "jilkine_size_response.csv": jilkine_competition,
        "jilkine_gradient_reversal.csv": jilkine_reversal,
        "jilkine_dynamic_tracking.csv": jilkine_tracking,
        "debelly_single_site_propagation.csv": debelly_propagation,
        "debelly_size_coupling.csv": debelly_propagation,
        "debelly_dual_site_challenge.csv": debelly_dual,
        "debelly_gradient_reversal.csv": debelly_reversal,
        "debelly_dynamic_tracking.csv": debelly_tracking,
        "debelly_intact_vs_both_off.csv": dependency["pairs"],
    }
    for filename, table in outputs.items():
        _write(arguments.output_root / filename, table)

    _plot_protocols(arguments.output_root / "campaign_protocols.png")
    _plot_jilkine_competition(jilkine_competition,
                              arguments.output_root / "jilkine_competing_fronts_size_scan.png")
    _plot_jilkine_competition_examples(rows, arguments.output_root,
                                       arguments.output_root / "jilkine_competing_fronts_examples.png")
    _plot_jilkine_reversal(jilkine_reversal,
                           arguments.output_root / "jilkine_gradient_reversal.png")
    _plot_reversal_examples(rows, task_index, arguments.output_root,
                            arguments.output_root / "jilkine_gradient_reversal_examples.png",
                            debelly=False)
    _plot_tracking_grid(jilkine_tracking,
                        arguments.output_root / "jilkine_dynamic_tracking.png", debelly=False)
    _plot_tracking_examples(rows, arguments.output_root,
                            arguments.output_root / "jilkine_dynamic_tracking_examples.png",
                            debelly=False)
    _plot_debelly_propagation(debelly_propagation,
                              arguments.output_root / "debelly_single_site_size_scan.png")
    _plot_debelly_propagation_examples(rows, arguments.output_root,
                                       arguments.output_root / "debelly_single_site_examples.png")
    _plot_debelly_dual(debelly_dual,
                       arguments.output_root / "debelly_dual_site_challenge.png")
    _plot_debelly_dual_examples(rows, arguments.output_root,
                                arguments.output_root / "debelly_dual_site_examples.png")
    _plot_debelly_reversal(debelly_reversal,
                           arguments.output_root / "debelly_gradient_reversal.png")
    _plot_reversal_examples(rows, task_index, arguments.output_root,
                            arguments.output_root / "debelly_gradient_reversal_examples.png",
                            debelly=True)
    _plot_tracking_grid(debelly_tracking,
                        arguments.output_root / "debelly_dynamic_tracking.png", debelly=True)
    _plot_tracking_examples(rows, arguments.output_root,
                            arguments.output_root / "debelly_dynamic_tracking_examples.png",
                            debelly=True)

    # Preserve the two original filenames as aliases, but replace their
    # ambiguous contents with explicitly named assay figures.
    shutil.copyfile(arguments.output_root / "jilkine_competing_fronts_size_scan.png",
                    arguments.output_root / "jilkine_size_response.png")
    shutil.copyfile(arguments.output_root / "debelly_single_site_size_scan.png",
                    arguments.output_root / "debelly_size_coupling.png")

    summary = {key: value for key, value in dependency.items() if key != "pairs"}
    summary.update(
        {
            "campaign_task_count": len(rows),
            "jilkine_competition_seeded_replicates": 5,
            "single_run_assays": [
                "Jilkine LEGI competing fronts",
                "all gradient reversal conditions",
                "all dynamic tracking conditions",
                "all De Belly conditions",
            ],
            "size_scan_design_limit": (
                "Localized cue widths are fixed fractions of length, so physical cue "
                "width increases with cell length. Size effects therefore combine domain "
                "size and cue-width effects."
            ),
            "jilkine_scope": (
                "Results describe the current equation-level implementations under the "
                "declared protocols. They are not claimed to reproduce transition lengths "
                "in the Jilkine publication."
            ),
            "candidate_prediction": (
                "Under the fractional-width single-site protocol, intact mechanics gives "
                "a much larger distal Rac response than dual-arm inhibition. A WGD-specific "
                "amplification is not established."
            ),
        }
    )
    (arguments.output_root / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(f"wrote {len(outputs)} tables and 15 figures to {arguments.output_root}")


if __name__ == "__main__":
    main()
