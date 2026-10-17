#!/usr/bin/env python3
"""Build the near-manuscript-ready WP1 Figures 1 and 2 package."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PHASE1 = ROOT / "results" / "campaigns" / "phase1_production_v1"
DECISIVE = ROOT / "results" / "campaigns" / "wp1_decisive_v1"
MANIFEST = ROOT / "campaigns" / "wp1_decisive_v1.json"
PHASE1_MANIFEST = ROOT / "campaigns" / "phase1_production_v1.json"
OTSUJI_SINGLE = ROOT / "results" / "campaigns" / "otsuji_single_cue_v1"
OTSUJI_GEOMETRY = ROOT / "results" / "campaigns" / "otsuji_reversal_geometry_sweep_l15_v1"
TRACKING_SUPPLEMENT = ROOT / "results" / "campaigns" / "wp1_tracking_supplement_v1"
TRACKING_AUDIT = ROOT / "results" / "campaigns" / "wp1_tracking_reliability_audit_v1"
REVISION_AUDIT = ROOT / "results" / "campaigns" / "wp1_revision_audit_v1"
LONG_REVERSAL = ROOT / "results" / "campaigns" / "jilkine_reversal_long_conditioning_v1"
FIGURE1_EXAMPLES = ROOT / "results" / "campaigns" / "wp1_figure1_examples_v1"
FIGURE1_COMMON = ROOT / "results" / "campaigns" / "wp1_figure1_common_horizon_v1"
DEFAULT_OUTPUT = ROOT / "results" / "manuscript_wp1_v3"

MODELS = (
    "jilkine_wave_pinning",
    "jilkine_otsuji",
    "jilkine_goryachev",
    "jilkine_legi",
)
MODEL_LABELS = {
    "jilkine_wave_pinning": "Wave-pinning",
    "jilkine_otsuji": "Otsuji",
    "jilkine_goryachev": "Goryachev",
    "jilkine_legi": "LEGI†",
}
MODEL_COLORS = {
    "jilkine_wave_pinning": "#2a9d8f",
    "jilkine_otsuji": "#577590",
    "jilkine_goryachev": "#90be6d",
    "jilkine_legi": "#e9c46a",
}
VARIANTS = (
    "intact",
    "rac_to_rho_half",
    "rho_to_rac_half",
    "rac_to_rho_off",
    "rho_to_rac_off",
    "both_off",
)
VARIANT_LABELS = {
    "intact": "Intact",
    "rac_to_rho_half": "Rac→mechanics→Rho 50%",
    "rho_to_rac_half": "Rho→mechanics→Rac 50%",
    "rac_to_rho_off": "Rac→mechanics→Rho off",
    "rho_to_rac_off": "Rho→mechanics→Rac off",
    "both_off": "Both mechanical arms off",
}
VARIANT_COLORS = dict(zip(VARIANTS, plt.cm.viridis(np.linspace(0.05, 0.95, 6))))
MAIN_VARIANTS = ("intact", "rac_to_rho_off", "rho_to_rac_off", "both_off")
VARIANT_COLORS.update(
    {
        "intact": "#4b176f",
        "rac_to_rho_off": "#d95f45",
        "rho_to_rac_off": "#2f7ebc",
        "both_off": "#2ca58d",
    }
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0]) if rows else []
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _boolean(row: dict[str, str], key: str) -> bool:
    return row.get(key, "").lower() == "true"


def _panel(axis, label: str) -> None:
    axis.text(
        -0.12, 1.06, label, transform=axis.transAxes, fontsize=13,
        fontweight="bold", va="top", ha="left",
    )


def _style() -> None:
    mpl.rcParams.update(
        {
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


def _box(axis, xy, width, height, text, color, fontsize=8.5) -> None:
    patch = FancyBboxPatch(
        xy, width, height, boxstyle="round,pad=0.02,rounding_size=0.02",
        transform=axis.transAxes, facecolor=color, edgecolor="#333333", lw=0.8,
    )
    axis.add_patch(patch)
    axis.text(
        xy[0] + width / 2, xy[1] + height / 2, text,
        transform=axis.transAxes, ha="center", va="center", fontsize=fontsize,
    )


def _arrow(axis, start, end, color="#333333", connectionstyle="arc3", style="-|>") -> None:
    axis.add_patch(
        FancyArrowPatch(
            start, end, transform=axis.transAxes, arrowstyle=style,
            mutation_scale=10, lw=1.1, color=color, connectionstyle=connectionstyle,
        )
    )


def _load_campaign_trajectory(campaign: Path, task_id: str) -> dict[str, np.ndarray]:
    path = campaign / "tasks" / task_id / "trajectory.npz"
    with np.load(path, allow_pickle=False) as payload:
        return {key: payload[key] for key in payload.files}


def _observable(payload: dict[str, np.ndarray]) -> np.ndarray:
    name = str(payload["observable_name"])
    return payload[f"field__{name}"]


def _normalized_profile(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    span = float(np.max(values) - np.min(values))
    return (values - np.min(values)) / span if span > 1e-12 else np.zeros_like(values)


def _localized_cue(relative_x: np.ndarray, sides=("left",), asymmetry=0.0) -> np.ndarray:
    cue = np.zeros_like(relative_x)
    if "left" in sides:
        mask = relative_x <= 0.1
        cue[mask] += 0.5 * (1.0 + np.cos(np.pi * relative_x[mask] / 0.1))
    if "right" in sides:
        mask = relative_x >= 0.9
        cue[mask] += (1.0 - asymmetry) * 0.5 * (
            1.0 + np.cos(np.pi * (1.0 - relative_x[mask]) / 0.1)
        )
    return cue


def _figure1_schematic(axis) -> None:
    axis.axis("off")
    _panel(axis, "A")
    axis.text(0.0, 0.97, "Common comparison framework", transform=axis.transAxes,
              fontweight="bold", va="top", fontsize=10)
    labels = ("Wave-pinning", "Otsuji", "Goryachev", "LEGI sensing")
    colors = ("#d8f0eb", "#dbe5f1", "#e5f2d6", "#f8eac3")
    for index, (label, color) in enumerate(zip(labels, colors)):
        _box(axis, (0.02 + 0.245 * index, 0.62), 0.205, 0.18, label, color, fontsize=8)
    _arrow(axis, (0.50, 0.60), (0.50, 0.44))
    _box(axis, (0.08, 0.22), 0.24, 0.15, "Two competing\nend cues", "#f1f1f1")
    _box(axis, (0.38, 0.22), 0.24, 0.15, "Gradient\nreversal", "#f1f1f1")
    _box(axis, (0.68, 0.22), 0.24, 0.15, "Alternating-gradient\ntracking", "#f1f1f1")
    axis.text(
        0.5, 0.06,
        "Common outputs: front number • resolution latency • final state • reversal • response bandwidth",
        transform=axis.transAxes, ha="center", color="#444444", fontsize=7.5,
    )


def _competition_grid(axis, rows: list[dict[str, str]]) -> None:
    lengths = sorted({_number(row, "length_um") for row in rows})
    codes = np.zeros((len(MODELS), len(lengths)), dtype=int)
    letters = np.empty_like(codes, dtype=object)
    for row in rows:
        i = MODELS.index(row["model"])
        j = lengths.index(_number(row, "length_um"))
        fractions = {
            0: _number(row, "competition_not_established_fraction") or 0.0,
            1: _number(row, "competition_resolved_fraction") or 0.0,
            2: _number(row, "competition_right_censored_fraction") or 0.0,
        }
        order = sorted(fractions, key=fractions.get, reverse=True)
        codes[i, j] = 3 if np.isclose(fractions[order[0]], fractions[order[1]]) else order[0]
        endpoint = {
            "U": _number(row, "final_unpolarized_fraction") or 0.0,
            "S": _number(row, "final_single_front_fraction") or 0.0,
            "M": _number(row, "final_multifront_fraction") or 0.0,
        }
        letters[i, j] = max(endpoint, key=endpoint.get)
    colors = ("#d9d9d9", "#3b82b8", "#e07a32", "#9b8bb8")
    axis.imshow(
        codes, aspect="auto", cmap=ListedColormap(colors),
        norm=BoundaryNorm(np.arange(-0.5, 4.5), 4),
    )
    for i in range(len(MODELS)):
        for j in range(len(lengths)):
            axis.text(j, i, letters[i, j], ha="center", va="center", fontsize=5.6,
                      fontweight="bold")
    axis.set_yticks(range(len(MODELS)), [MODEL_LABELS[value] for value in MODELS])
    shown = list(range(0, len(lengths), 2))
    axis.set_xticks(shown, [f"{lengths[index]:g}" for index in shown])
    axis.set_xlabel("Cell length (µm)")
    axis.set_title("Competing-front outcome and final state", loc="left")
    _panel(axis, "B")
    axis.legend(
        handles=[
            Patch(color=colors[0], label="Two fronts not sustained"),
            Patch(color=colors[1], label="Resolved"),
            Patch(color=colors[2], label="Not resolved by run end"),
            Patch(color=colors[3], label="Mixed seeds"),
        ],
        frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, -0.28),
    )


def _competition_latency(axis, rows: list[dict[str, str]]) -> None:
    for model in MODELS:
        selected = [
            row for row in rows
            if row["model"] == model and _number(row, "resolution_latency_median_s") is not None
        ]
        axis.plot(
            [_number(row, "length_um") for row in selected],
            [_number(row, "resolution_latency_median_s") for row in selected],
            "o-", lw=1.8, ms=3.5, color=MODEL_COLORS[model], label=MODEL_LABELS[model],
        )
    axis.set_xlabel("Cell length (µm)")
    axis.set_ylabel("Resolution latency (s)")
    axis.set_yscale("log")
    axis.grid(alpha=0.18)
    axis.set_title("Resolved systems slow with size", loc="left")
    axis.legend(frameon=False)
    _panel(axis, "C")


def _reversal_grid(axis, rows: list[dict[str, str]]) -> None:
    lengths = sorted({_number(row, "length_um") for row in rows})
    status_order = {"not_conditioned": 0, "success": 1, "right_censored": 2}
    matrix = np.zeros((len(MODELS), len(lengths)), dtype=int)
    for row in rows:
        matrix[MODELS.index(row["model"]), lengths.index(_number(row, "length_um"))] = status_order[row["reversal_status"]]
    colors = ("#d9d9d9", "#3aa76d", "#d95d4f")
    axis.imshow(
        matrix, aspect="auto", cmap=ListedColormap(colors),
        norm=BoundaryNorm(np.arange(-0.5, 3.5), 3),
    )
    axis.set_yticks(range(len(MODELS)), [MODEL_LABELS[value] for value in MODELS])
    shown = list(range(0, len(lengths), 2))
    axis.set_xticks(shown, [f"{lengths[index]:g}" for index in shown])
    axis.set_xlabel("Cell length (µm)")
    axis.set_title("Polarity replacement after cue reversal", loc="left")
    axis.legend(
        handles=[
            Patch(color=colors[0], label="Not conditioned"),
            Patch(color=colors[1], label="Reversed"),
            Patch(color=colors[2], label="Not reversed by run end"),
        ],
        frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.28),
    )
    _panel(axis, "D")


def _reversal_latency(axis, rows: list[dict[str, str]]) -> None:
    for model in MODELS:
        selected = [
            row for row in rows
            if row["model"] == model and row["reversal_status"] == "success"
        ]
        axis.plot(
            [_number(row, "length_um") for row in selected],
            [_number(row, "reversal_latency_s") for row in selected],
            "o-", lw=1.7, ms=3.2, color=MODEL_COLORS[model], label=MODEL_LABELS[model],
        )
        censored = [
            row for row in rows
            if row["model"] == model and row["reversal_status"] == "right_censored"
        ]
        axis.scatter(
            [_number(row, "length_um") for row in censored],
            [_number(row, "post_reversal_observation_s") for row in censored],
            marker="v", facecolor="none", edgecolor=MODEL_COLORS[model], s=28,
        )
    axis.set_xlabel("Cell length (µm)")
    axis.set_ylabel("Reversal latency (s)")
    axis.set_yscale("log")
    axis.grid(alpha=0.18)
    axis.set_title("Reversal latency", loc="left")
    _panel(axis, "E")


def _tracking_heatmap(axis, rows: list[dict[str, str]]) -> None:
    periods = sorted({_number(row, "period_s") for row in rows})
    lengths = sorted({_number(row, "length_um") for row in rows})
    row_keys = [(model, length) for model in MODELS for length in lengths]
    matrix = np.full((len(row_keys), len(periods)), np.nan)
    reliable = np.zeros_like(matrix, dtype=bool)
    for row in rows:
        i = row_keys.index((row["model"], _number(row, "length_um")))
        j = periods.index(_number(row, "period_s"))
        matrix[i, j] = _number(row, "response_amplitude")
        reliable[i, j] = _boolean(row, "tracking_reliable")
    image = axis.imshow(matrix, aspect="auto", cmap="magma", vmin=0, vmax=0.85)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if not reliable[i, j]:
                axis.text(j, i, "×", ha="center", va="center", color="white", fontsize=7)
    axis.set_xticks(range(len(periods)), [f"{value:g}" for value in periods])
    axis.set_xlabel("Forcing period (s)")
    axis.set_yticks(
        [index * len(lengths) + (len(lengths) - 1) / 2 for index in range(len(MODELS))],
        [MODEL_LABELS[value] for value in MODELS],
    )
    for boundary in range(1, len(MODELS)):
        axis.axhline(boundary * len(lengths) - 0.5, color="white", lw=0.8)
    axis.set_title(
        "Tracking amplitude\n(rows within model: 5, 15, 25, 39 µm)",
        loc="left", fontsize=8.5,
    )
    plt.colorbar(image, ax=axis, fraction=0.045, pad=0.03, label="Orientation amplitude")
    _panel(axis, "F")


def _otsuji_single_cue_result(
    metric_axis, profile_axis, rows: list[dict[str, str]], source: Path
) -> None:
    colors = {0.04: "#577590", 0.2: "#d95f45"}
    styles = {
        0.04: {"marker": "o", "mfc": "none", "ms": 5.0, "lw": 2.0, "zorder": 3},
        0.2: {"marker": "o", "mfc": colors[0.2], "ms": 2.7, "lw": 1.2, "zorder": 4},
    }
    for amplitude in (0.04, 0.2):
        selected = [row for row in rows if np.isclose(_number(row, "input_peak"), amplitude)]
        metric_axis.plot(
            [_number(row, "length_um") for row in selected],
            [abs(_number(row, "final_orientation") or 0.0) for row in selected],
            color=colors[amplitude], label=f"S = {amplitude:g}",
            markeredgecolor=colors[amplitude], **styles[amplitude],
        )
    metric_axis.axhline(0.05, color="#777777", ls="--", lw=0.8, label="orientation threshold")
    metric_axis.set(xlabel="Cell length (µm)", ylabel="Final cue-aligned orientation", ylim=(-0.03, 1.03))
    metric_axis.set_title("Final Otsuji orientation", loc="left")
    metric_axis.grid(alpha=0.16)
    metric_axis.legend(frameon=False, fontsize=6.8, loc="lower right")
    _panel(metric_axis, "A")

    manifest = json.loads((ROOT / "campaigns" / "otsuji_single_cue_v1.json").read_text())
    task_lookup = {task["task_id"]: task for task in manifest["tasks"]}
    selected_lengths = (5.0, 13.0, 21.0, 31.0, 39.0)
    exported = []
    for length, color in zip(selected_lengths, plt.cm.Blues(np.linspace(0.38, 0.95, 5))):
        row = next(
            item for item in rows
            if np.isclose(_number(item, "input_peak"), 0.04)
            and np.isclose(_number(item, "length_um"), length)
        )
        payload = _load_campaign_trajectory(OTSUJI_SINGLE, row["task_id"])
        relative_x = payload["x"] / payload["x"][-1]
        response = _normalized_profile(_observable(payload)[-1])
        profile_axis.plot(relative_x, response, lw=1.4, color=color, label=f"{length:g} µm")
        cue = _localized_cue(relative_x)
        for x, value, stimulus in zip(relative_x, response, cue):
            exported.append(
                {
                    "task_id": row["task_id"], "length_um": length,
                    "input_peak": task_lookup[row["task_id"]]["protocol"]["amplitude"],
                    "relative_position": float(x), "normalized_final_active": float(value),
                    "normalized_initial_cue": float(stimulus),
                }
            )
    relative_x = np.linspace(0, 1, 401)
    profile_axis.fill_between(
        relative_x, 0, _localized_cue(relative_x), color="#bbbbbb", alpha=0.28,
        label="Transient cue",
    )
    profile_axis.set(xlabel="Relative position", ylabel="Normalized active Otsuji")
    profile_axis.set_title("Final profiles (S = 0.04)", loc="left")
    profile_axis.legend(frameon=False, ncol=2, fontsize=6.5)
    _write(source / "figure1_otsuji_single_cue_profiles.csv", exported)


def _competition_examples(axes, rows: list[dict[str, str]], source: Path) -> None:
    exported = []
    task_rows = _read(PHASE1 / "task_summary.csv")
    for index, (axis, model) in enumerate(zip(axes, MODELS)):
        for length, ls in ((5.0, "-"), (39.0, "--")):
            task_row = next(
                item for item in task_rows
                if item["assay"] == "competing_fronts" and item["model"] == model
                and np.isclose(_number(item, "length_um"), length)
                and item.get("replicate") == "0"
            )
            payload = _load_campaign_trajectory(PHASE1, task_row["task_id"])
            relative_x = payload["x"] / payload["x"][-1]
            response = _normalized_profile(_observable(payload)[-1])
            axis.plot(relative_x, response, ls=ls, lw=1.45, color=MODEL_COLORS[model], label=f"{length:g} µm")
            for x, value in zip(relative_x, response):
                exported.append(
                    {"model": model, "task_id": task_row["task_id"], "length_um": length,
                     "relative_position": float(x), "normalized_final_response": float(value)}
                )
        relative_x = np.linspace(0, 1, 401)
        cue = _localized_cue(relative_x, sides=("left", "right"), asymmetry=0.02)
        axis.fill_between(relative_x, 0, cue, color="#bbbbbb", alpha=0.25)
        axis.set_ylim(-0.03, 1.05)
        axis.set_title(MODEL_LABELS[model], fontsize=8.2)
        axis.set_xticks((0, 0.5, 1))
        if index == 0:
            axis.set_ylabel("Normalized response")
            axis.legend(frameon=False, fontsize=6.5)
            _panel(axis, "B")
        else:
            axis.set_yticklabels([])
        if index == len(MODELS) - 1:
            axis.text(0.98, 0.94, "gray: transient cues", transform=axis.transAxes,
                      ha="right", va="top", fontsize=6.2, color="#555555")
    axes[0].text(
        0.0, 1.22, "Two competing cues: model-specific large-cell outcomes",
        transform=axes[0].transAxes, fontsize=10, fontweight="bold", ha="left",
    )
    _write(source / "figure1_competing_front_examples.csv", exported)


def _reversal_examples(axes, rows: list[dict[str, str]], source: Path) -> None:
    exported = []
    for index, (axis, model) in enumerate(zip(axes, MODELS)):
        for length, ls in ((5.0, "-"), (39.0, "--")):
            row = next(
                item for item in rows
                if item["model"] == model and np.isclose(_number(item, "length_um"), length)
            )
            payload = _load_campaign_trajectory(PHASE1, row["task_id"])
            observable = _observable(payload)
            relative_x = payload["x"] / payload["x"][-1]
            left = np.mean(observable[:, relative_x <= 0.2], axis=1)
            right = np.mean(observable[:, relative_x >= 0.8], axis=1)
            orientation = (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)
            manifest = json.loads(PHASE1_MANIFEST.read_text())
            task = next(task for task in manifest["tasks"] if task["task_id"] == row["task_id"])
            reversal_time = float(task["protocol"]["reversal_time"])
            time = payload["time"] - reversal_time
            axis.plot(time, orientation, ls=ls, lw=1.35, color=MODEL_COLORS[model], label=f"{length:g} µm")
            for t, value in zip(time, orientation):
                exported.append(
                    {"model": model, "task_id": row["task_id"], "length_um": length,
                     "time_from_reversal_s": float(t), "orientation": float(value)}
                )
        axis.axvline(0, color="black", lw=0.75, ls=":")
        axis.axhline(0, color="#777777", lw=0.6)
        axis.set_title(MODEL_LABELS[model], fontsize=8.2)
        if index == 0:
            axis.set_ylabel("Response orientation")
            axis.legend(frameon=False, fontsize=6.5)
            _panel(axis, "D")
        else:
            axis.set_yticklabels([])
        axis.set_xlabel("")
    axes[0].text(
        0.0, 1.22, "Cue reverses at t = 0: representative responses",
        transform=axes[0].transAxes, fontsize=10, fontweight="bold", ha="left",
    )
    axes[-1].annotate("cue: left → right", xy=(0.98, 0.92), xycoords="axes fraction",
                      ha="right", fontsize=6.5, color="#444444")
    axes[1].text(1.08, -0.18, "Time from cue reversal (s)", transform=axes[1].transAxes,
                 ha="center", va="top", fontsize=8.2)
    _write(source / "figure1_reversal_examples.csv", exported)


def _tracking_facets(axes, rows: list[dict[str, str]]) -> None:
    periods = sorted({_number(row, "period_s") for row in rows})
    lengths = sorted({_number(row, "length_um") for row in rows})
    images = [None, None]
    metric_specs = (
        ("response_amplitude", "magma", 0.0, 0.85, "Amp.", (0, 0.4, 0.8)),
        ("phase_lag_fraction_of_period", "coolwarm", -0.25, 0.25,
         "Phase / T", (-0.25, 0, 0.25)),
    )
    for row_index, model in enumerate(MODELS):
        matrices = [np.full((len(lengths), len(periods)), np.nan) for _ in metric_specs]
        reliable = np.zeros_like(matrices[0], dtype=bool)
        for row in rows:
            if row["model"] != model:
                continue
            i = lengths.index(_number(row, "length_um"))
            j = periods.index(_number(row, "period_s"))
            for metric_index, (metric, *_unused) in enumerate(metric_specs):
                matrices[metric_index][i, j] = _number(row, metric)
            reliable[i, j] = _boolean(row, "tracking_reliable_audited")
        for metric_index, (metric, cmap, vmin, vmax, legend_title, ticks) in enumerate(metric_specs):
            del metric, legend_title, ticks
            axis = axes[row_index][metric_index]
            plot_cmap = mpl.colormaps[cmap].copy()
            plot_cmap.set_bad("#d9d9d9")
            display_matrix = np.ma.masked_where(~reliable, matrices[metric_index])
            images[metric_index] = axis.imshow(
                display_matrix, aspect="auto", cmap=plot_cmap, vmin=vmin, vmax=vmax
            )
            for i in range(len(lengths)):
                for j in range(len(periods)):
                    if not reliable[i, j]:
                        axis.text(j, i, "×", ha="center", va="center", color="#555555", fontsize=6)
            sample_i = lengths.index(15.0)
            sample_j = periods.index(200.0)
            axis.add_patch(
                mpl.patches.Rectangle(
                    (sample_j - 0.5, sample_i - 0.5), 1, 1,
                    fill=False, edgecolor="#d62728", linewidth=1.35,
                )
            )
            axis.set_yticks(range(len(lengths)))
            if metric_index == 0:
                axis.set_yticklabels([f"{value:g}" for value in lengths])
                axis.set_ylabel(
                    MODEL_LABELS[model], rotation=90, ha="center", va="bottom",
                    labelpad=12, fontsize=6.1,
                )
            else:
                axis.set_yticklabels([])
            axis.set_xticks(range(len(periods)))
            if row_index == len(MODELS) - 1:
                axis.set_xticklabels([f"{value:g}" for value in periods])
                axis.set_xlabel("Period (s)", fontsize=7)
            else:
                axis.set_xticklabels([])
    for metric_index, (_metric, _cmap, _vmin, _vmax, legend_title, ticks) in enumerate(metric_specs):
        color_axis = axes[0][metric_index].inset_axes((0.87, 0.10, 0.035, 0.72))
        colorbar = plt.colorbar(images[metric_index], cax=color_axis, ticks=ticks)
        colorbar.ax.tick_params(labelsize=4.8, length=1.5, pad=1)
        colorbar.ax.set_title(legend_title, fontsize=5.3, pad=1.5)


def _cue_kymograph(axis, kind: str, *, period: float = 100.0) -> None:
    position = np.linspace(0, 1, 201)
    if kind == "competition":
        times = np.linspace(0, 80, 241)
        envelope = np.where(
            times <= 20, 1.0,
            np.where(times < 40, 0.5 * (1 + np.cos(np.pi * (times - 20) / 20)), 0.0),
        )
        matrix = envelope[:, None] * _localized_cue(
            position, sides=("left", "right"), asymmetry=0.02
        )[None, :]
        note = "Transient cues: plateau 0–20 s, taper to zero at 40 s\nObservation continues to 8,000–30,420 s"
        xlabel = "Time (s)"
    elif kind == "reversal":
        times = np.linspace(-200, 400, 301)
        left = 1 - position
        right = position
        matrix = np.empty((len(times), len(position)))
        for index, time in enumerate(times):
            if time < 0:
                matrix[index] = left
            elif time >= 20:
                matrix[index] = right
            else:
                weight = time / 20
                matrix[index] = (1 - weight) * left + weight * right
        note = "Persistent gradient; 20-s cross-fade at t = 0\nConditioning 2,000–6,084 s; follow-up 4,000–12,168 s"
        xlabel = "Time (s)"
    elif kind == "alternating":
        times = np.linspace(0, 0.5 * period, 301)
        centered = 2 * position - 1
        matrix = np.asarray(
            [0.5 * (1 - np.cos(2 * np.pi * time / period) * centered) for time in times]
        )
        note = ""
        xlabel = "Time (s)"
    else:
        raise ValueError(kind)
    image = axis.imshow(
        matrix.T, origin="lower", aspect="auto",
        extent=(times[0], times[-1], 0, 1), cmap="Greys", vmin=0, vmax=1,
    )
    axis.set(xlabel=xlabel, ylabel="Rel. pos.")
    axis.set_ylim(-0.03, 1.03)
    if kind == "competition":
        axis.text(times[0] + 0.02 * (times[-1] - times[0]), 0.06, "left cue",
                  fontsize=6.2, color="#333333", va="bottom")
        axis.text(times[0] + 0.02 * (times[-1] - times[0]), 0.94, "right cue",
                  fontsize=6.2, color="#333333", va="top")
    if kind == "competition":
        axis.axvline(20, color="#d62728", ls="--", lw=1.0)
    elif kind == "reversal":
        for sample_time in (-1e-6, 20, 100):
            axis.axvline(sample_time, color="#d62728", ls="--", lw=1.0)
    elif kind == "alternating":
        for phase_time in (period / 8, 2 * period / 8, 3 * period / 8):
            axis.axvline(phase_time, color="#d62728", ls="--", lw=1.05)
    return image


def _profile_snapshots(
    axes, task_rows: list[dict[str, str]], campaign: Path, assay: str,
    length: float, time_selector, source_path: Path,
) -> None:
    exported = []
    for index, (axis, model) in enumerate(zip(axes, MODELS)):
        row = next(
            item for item in task_rows
            if item["assay"] == assay and item["model"] == model
            and np.isclose(_number(item, "length_um"), length)
            and item.get("replicate", "0") == "0"
        )
        payload = _load_campaign_trajectory(campaign, row["task_id"])
        requested = time_selector(row, payload)
        indices = [int(np.argmin(np.abs(payload["time"] - time))) for time in requested]
        observable = _observable(payload)
        selected_values = observable[indices]
        minimum, maximum = float(np.min(selected_values)), float(np.max(selected_values))
        scale = max(maximum - minimum, 1e-12)
        relative_x = payload["x"] / payload["x"][-1]
        for color, time_index in zip(plt.cm.viridis(np.linspace(0.08, 0.90, len(indices))), indices):
            profile = (observable[time_index] - minimum) / scale
            actual_time = float(payload["time"][time_index])
            axis.plot(relative_x, profile, color=color, lw=1.25, label=time_selector.label(actual_time, row))
            for x, value in zip(relative_x, observable[time_index]):
                exported.append(
                    {"model": model, "task_id": row["task_id"], "length_um": length,
                     "time_s": actual_time, "relative_position": float(x),
                     "observable": str(payload["observable_name"]), "value": float(value)}
                )
        axis.set(xlim=(0, 1), ylim=(-0.03, 1.03), title=MODEL_LABELS[model])
        axis.set_xticks((0, 0.5, 1))
        if index == 0:
            axis.set_ylabel("Scaled response")
        else:
            axis.set_yticklabels([])
        axis.set_xlabel("Relative position")
    axes[-1].legend(frameon=False, fontsize=5.8, loc="upper left", bbox_to_anchor=(1.02, 1.02))
    _write(source_path, exported)


class _CompetitionTimes:
    def __call__(self, row, payload):
        del row
        return (0.0, 40.0, min(2000.0, float(payload["time"][-1])), float(payload["time"][-1]))

    @staticmethod
    def label(time, row):
        del row
        return f"t={time:g} s"


class _ReversalTimes:
    def __init__(self, manifest: dict[str, dict]):
        self.manifest = manifest

    def __call__(self, row, payload):
        reversal = float(self.manifest[row["task_id"]]["protocol"]["reversal_time"])
        return (reversal - 20.0, reversal + 100.0, reversal + 500.0, float(payload["time"][-1]))

    def label(self, time, row):
        reversal = float(self.manifest[row["task_id"]]["protocol"]["reversal_time"])
        offset = time - reversal
        return "run end" if np.isclose(time, float(self.manifest[row["task_id"]]["simulation"]["end_time"])) else f"{offset:+g} s"


def _focused_profile_snapshots(
    axes, assay: str, offsets: tuple[float, ...], source_path: Path,
    *, origin: str = "absolute", zero_from_left: bool = False,
) -> None:
    """Plot retained trajectories from the focused, common-time Figure 1 runs."""
    exported = []
    colors = plt.cm.viridis(np.linspace(0.08, 0.90, len(offsets)))
    short_model = {model: model.removeprefix("jilkine_") for model in MODELS}
    for index, (axis, model) in enumerate(zip(axes, MODELS)):
        directory = FIGURE1_EXAMPLES / "tasks" / f"{assay}_{short_model[model]}"
        with np.load(directory / "trajectory.npz", allow_pickle=False) as stored:
            payload = {key: stored[key] for key in stored.files}
        config = json.loads(
            next((FIGURE1_EXAMPLES / "configs").glob(f"*_{assay}_{short_model[model]}.json")).read_text()
        )
        if origin == "reversal":
            time_zero = float(config["protocol"]["reversal_time"])
        elif origin == "tracking_cycle":
            time_zero = float(config["simulation"]["end_time"]) - float(config["protocol"]["period"])
        else:
            time_zero = 0.0
        requested = []
        for value in offsets:
            if zero_from_left and value == 0:
                requested.append(float(np.nextafter(time_zero, -np.inf)))
            else:
                requested.append(time_zero + value)
        observable = _observable(payload)
        selected = np.asarray([
            [np.interp(time, payload["time"], observable[:, position_index])
             for position_index in range(observable.shape[1])]
            for time in requested
        ])
        minimum, maximum = float(np.min(selected)), float(np.max(selected))
        scale = max(maximum - minimum, 1e-12)
        relative_x = payload["x"] / payload["x"][-1]
        for color, offset, actual_time, snapshot in zip(colors, offsets, requested, selected):
            values = (snapshot - minimum) / scale
            if assay == "tracking":
                label = {25.0: "π/8", 50.0: "2π/8", 75.0: "3π/8"}[offset]
            elif zero_from_left and offset == 0:
                label = "0⁻ s"
            else:
                label = f"{offset:g} s"
            axis.plot(relative_x, values, color=color, lw=1.25, label=label)
            for x, value in zip(relative_x, snapshot):
                exported.append(
                    {"model": model, "assay": assay, "length_um": 15.0,
                     "time_s": actual_time, "time_from_reference_s": actual_time - time_zero,
                     "relative_position": float(x), "observable": str(payload["observable_name"]),
                     "value": float(value)}
                )
        axis.set(xlim=(0, 1), ylim=(-0.03, 1.03))
        axis.set_xticks((0, 0.5, 1))
        axis.text(0.04, 0.92, MODEL_LABELS[model], transform=axis.transAxes,
                  ha="left", va="top", fontsize=7.2)
        if index % 2 == 0:
            axis.set_ylabel("Scaled resp.")
        else:
            axis.set_yticklabels([])
        if index >= 2:
            axis.set_xlabel("Rel. pos.")
        else:
            axis.set_xticklabels([])
    axes[1].legend(frameon=False, fontsize=6.1, loc="lower right", title="Time")
    _write(source_path, exported)


def _resolved_segments(x: list[float], y: list[float]) -> list[tuple[np.ndarray, np.ndarray]]:
    """Split resolved points at gaps so shaded regions do not bridge failures."""
    if not x:
        return []
    order = np.argsort(x)
    sx, sy = np.asarray(x)[order], np.asarray(y)[order]
    cuts = np.where(np.diff(sx) > 2.1)[0] + 1
    return [(a, b) for a, b in zip(np.split(sx, cuts), np.split(sy, cuts)) if len(a)]


def _extended_competition_latency(axis, task_rows: list[dict[str, str]], source: Path) -> None:
    exported = []
    models = MODELS[:3]
    for model in models:
        model_rows = [row for row in task_rows if row["assay"] == "competing_fronts_extended_4x" and row["model"] == model]
        lengths = sorted({_number(row, "length_um") for row in model_rows})
        resolved_x, resolved_y = [], []
        for length in lengths:
            members = [row for row in model_rows if np.isclose(_number(row, "length_um"), length)]
            latencies = [_number(row, "competition_resolution_latency_s") for row in members]
            latencies = [value for value in latencies if value is not None]
            latencies = [value for value in latencies if value <= 20_000]
            if latencies:
                resolved_x.append(length); resolved_y.append(float(np.median(latencies)))
            exported.append(
                {"model": model, "length_um": length, "n": len(members),
                 "resolved_count_by_20000_s": len(latencies),
                 "resolution_latency_median_s": "" if not latencies else float(np.median(latencies)),
                 "common_observation_time_s": 20_000.0}
            )
        for segment_index, (segment_x, segment_y) in enumerate(_resolved_segments(resolved_x, resolved_y)):
            axis.plot(segment_x, segment_y, "o-", color=MODEL_COLORS[model], lw=1.7, ms=3.5,
                      label=MODEL_LABELS[model] if segment_index == 0 else None)
            axis.fill_between(segment_x, segment_y, 20_000, color=MODEL_COLORS[model], alpha=0.12)
    axis.axhline(200, color="#d62728", ls="--", lw=0.9)
    axis.axhline(20_000, color="#d62728", ls="--", lw=0.9)
    axis.set(xlabel="L (µm)", ylabel="Resolution time (s)", yscale="log", ylim=(20, 30_000))
    axis.grid(alpha=0.18)
    axis.legend(frameon=False, fontsize=6.5, loc="lower right")
    _write(source / "figure1_extended_competition.csv", exported)


def _long_reversal_latency(axis, rows: list[dict[str, str]], source: Path) -> None:
    exported = []
    for model in MODELS:
        selected = [row for row in rows if row["model"] == model]
        success = [row for row in selected if row["reversal_status"] == "success"
                   and (_number(row, "reversal_latency_s") or np.inf) <= 20_000]
        x = [_number(row, "length_um") for row in success]
        y = [_number(row, "reversal_latency_s") for row in success]
        for segment_index, (segment_x, segment_y) in enumerate(_resolved_segments(x, y)):
            axis.plot(segment_x, segment_y, "o-", color=MODEL_COLORS[model], lw=1.6, ms=3.2,
                      label=MODEL_LABELS[model] if segment_index == 0 else None)
            axis.fill_between(segment_x, segment_y, 20_000, color=MODEL_COLORS[model], alpha=0.12)
        for row in selected:
            exported.append(dict(row))
    axis.axhline(20, color="#d62728", ls="--", lw=0.9)
    axis.axhline(100, color="#d62728", ls="--", lw=0.9)
    axis.set(xlabel="L (µm)", ylabel="Reversal time (s)", yscale="log", ylim=(10, 30_000))
    axis.grid(alpha=0.18)
    axis.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9,
                ncol=2, fontsize=6.2, loc="upper left")
    _write(source / "figure1_long_reversal.csv", exported)


def build_figure1(output: Path, source: Path) -> None:
    common_tasks = _read(FIGURE1_COMMON / "task_summary.csv")
    competition_tasks = common_tasks
    competition_tasks = [row for row in competition_tasks if row["assay"] == "competing_fronts_extended_4x"]
    reversal = [row for row in common_tasks if row["assay"] == "gradient_reversal_long_conditioning"]
    tracking_rows = _read(PHASE1 / "jilkine_dynamic_tracking.csv")
    tracking_audit = _read(TRACKING_AUDIT / "tracking_reliability_audit.csv")
    audit_lookup = {
        (row["model"], _number(row, "length_um"), _number(row, "period_s")): row
        for row in tracking_audit
    }
    for row in tracking_rows:
        audit = audit_lookup[(row["model"], _number(row, "length_um"), _number(row, "period_s"))]
        row["tracking_reliable_published"] = row["tracking_reliable"]
        row["tracking_reliable_audited"] = audit["audited_reliable"]
        row["audit_failure_reasons"] = audit["audit_failure_reasons"]
        row["audit_harmonic_r_squared"] = audit["harmonic_r_squared"]
        row["audit_cycle_amplitude_cv"] = audit["cycle_amplitude_cv"]
        row["audit_cycle_phase_dispersion"] = audit["cycle_phase_dispersion"]
    _write(source / "figure1_dynamic_tracking.csv", tracking_rows)

    figure = plt.figure(figsize=(13.2, 11.0))
    grid = figure.add_gridspec(3, 3, hspace=0.42, wspace=0.36)

    competition_cue = figure.add_subplot(grid[0, 0])
    _cue_kymograph(competition_cue, "competition")
    _panel(competition_cue, "A")
    competition_profile_grid = grid[0, 1].subgridspec(2, 2, hspace=0.10, wspace=0.10)
    competition_axes = [
        figure.add_subplot(competition_profile_grid[i, j])
        for i in range(2) for j in range(2)
    ]
    _focused_profile_snapshots(
        competition_axes, "competition", (20.0, 200.0, 20_000.0),
        source / "figure1_competition_profiles.csv",
    )
    _panel(competition_axes[0], "B")
    competition_latency = figure.add_subplot(grid[0, 2])
    _extended_competition_latency(competition_latency, competition_tasks, source)
    competition_latency.axvline(15, color="black", ls=":", lw=1.1)
    competition_latency.text(15, 24, "L=15", rotation=90, ha="right", va="bottom", fontsize=6.4)
    _panel(competition_latency, "C")

    reversal_cue = figure.add_subplot(grid[1, 0])
    _cue_kymograph(reversal_cue, "reversal")
    _panel(reversal_cue, "D")
    reversal_profile_grid = grid[1, 1].subgridspec(2, 2, hspace=0.10, wspace=0.10)
    reversal_axes = [
        figure.add_subplot(reversal_profile_grid[i, j])
        for i in range(2) for j in range(2)
    ]
    _focused_profile_snapshots(
        reversal_axes, "reversal", (0.0, 20.0, 100.0),
        source / "figure1_reversal_profiles.csv",
        origin="reversal", zero_from_left=True,
    )
    _panel(reversal_axes[0], "E")
    reversal_latency = figure.add_subplot(grid[1, 2])
    _long_reversal_latency(reversal_latency, reversal, source)
    reversal_latency.axvline(15, color="black", ls=":", lw=1.1)
    reversal_latency.text(15, 12, "L=15", rotation=90, ha="right", va="bottom", fontsize=6.4)
    _panel(reversal_latency, "F")

    tracking_protocol = figure.add_subplot(grid[2, 0])
    _cue_kymograph(tracking_protocol, "alternating", period=200.0)
    _panel(tracking_protocol, "G")
    tracking_profile_grid = grid[2, 1].subgridspec(2, 2, hspace=0.10, wspace=0.10)
    tracking_profile_axes = [
        figure.add_subplot(tracking_profile_grid[i, j])
        for i in range(2) for j in range(2)
    ]
    _focused_profile_snapshots(
        tracking_profile_axes, "tracking", (25.0, 50.0, 75.0),
        source / "figure1_tracking_profiles.csv",
        origin="tracking_cycle",
    )
    _panel(tracking_profile_axes[0], "H")
    tracking_grid = grid[2, 2].subgridspec(4, 2, hspace=0.04, wspace=0.08)
    tracking_axes = [
        [figure.add_subplot(tracking_grid[row, column]) for column in range(2)]
        for row in range(4)
    ]
    _tracking_facets(tracking_axes, tracking_rows)
    _panel(tracking_axes[0][0], "I")
    figure.subplots_adjust(top=0.985, bottom=0.06, left=0.055, right=0.96)
    figure.savefig(output / "Figure1_WP1_cross_model.png", dpi=300)
    figure.savefig(output / "Figure1_WP1_cross_model.pdf")
    plt.close(figure)


def build_supplement_otsuji(output: Path, source: Path) -> None:
    rows = _read(OTSUJI_GEOMETRY / "reversal_geometry_sweep.csv")
    _write(source / "figureS1_otsuji_reversal_geometry_sweep.csv", rows)
    amplitudes = sorted({_number(row, "input_peak_per_s") for row in rows})
    widths = sorted({_number(row, "stimulus_width_fraction") for row in rows})
    success = np.zeros((len(widths), len(amplitudes)))
    centroid = np.full_like(success, np.nan, dtype=float)
    centroid_shift = np.full_like(success, np.nan, dtype=float)
    for row in rows:
        i = widths.index(_number(row, "stimulus_width_fraction"))
        j = amplitudes.index(_number(row, "input_peak_per_s"))
        success[i, j] = row["strict_reversal_status"] == "success"
        centroid[i, j] = _number(row, "final_excess_centroid")
        centroid_shift[i, j] = _number(row, "excess_centroid_shift")

    figure, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), gridspec_kw={"wspace": 0.52})
    relative_x = np.linspace(0, 1, 301)
    shown_widths = (0.10, 0.30, 0.55, 1.0)
    for color, width in zip(plt.cm.viridis(np.linspace(0.08, 0.90, len(shown_widths))), shown_widths):
        left = 0.2 * np.clip(1 - relative_x / width, 0, 1)
        right = left[::-1]
        axes[0].plot(relative_x, left, color=color, lw=1.7, label=f"{width:g}")
        axes[0].plot(relative_x, right, color=color, lw=1.0, ls="--")
    axes[0].set(xlabel="Rel. pos.", ylabel="Stimulus (s⁻¹)", ylim=(-0.01, 0.22))
    axes[0].legend(frameon=False, title="Width / L", ncol=2, loc="upper center")
    axes[0].text(0.03, 0.08, "solid: before   dashed: after", transform=axes[0].transAxes,
                 fontsize=6.6)
    _panel(axes[0], "A")

    if np.any(success):
        binary_map = ListedColormap(("#eeeeee", "#2a9d8f"))
        middle_image = axes[1].imshow(
            success, origin="lower", aspect="auto", cmap=binary_map, vmin=0, vmax=1
        )
    else:
        middle_image = axes[1].imshow(
            centroid_shift, origin="lower", aspect="auto", cmap="YlGn",
            vmin=0, vmax=max(float(np.nanmax(centroid_shift)), 1e-6),
        )
    axes[1].set(
        xlabel="Peak stimulus (s⁻¹)", ylabel="Width / L",
        xticks=range(len(amplitudes)), yticks=range(len(widths)),
        xticklabels=[f"{value:g}" for value in amplitudes],
        yticklabels=[f"{value:g}" for value in widths],
    )
    axes[1].tick_params(axis="x", rotation=45)
    if np.any(success):
        axes[1].legend(
            handles=(Patch(facecolor="#eeeeee", label="No reversal"),
                     Patch(facecolor="#2a9d8f", label="Reversal")),
            frameon=False, loc="upper left", fontsize=6.5,
        )
    else:
        plt.colorbar(middle_image, ax=axes[1], fraction=0.046, pad=0.03,
                     label="Response-centroid shift (x/L)")
        axes[1].text(0.03, 0.96, "0/64 strict reversals", transform=axes[1].transAxes,
                     ha="left", va="top", fontsize=7, color="black")
    _panel(axes[1], "B")

    centroid_min, centroid_max = float(np.nanmin(centroid)), float(np.nanmax(centroid))
    continuous = axes[2].imshow(
        centroid, origin="lower", aspect="auto", cmap="viridis",
        vmin=centroid_min, vmax=centroid_max,
    )
    axes[2].set(
        xlabel="Peak stimulus (s⁻¹)", ylabel="",
        xticks=range(len(amplitudes)), yticks=range(len(widths)),
        xticklabels=[f"{value:g}" for value in amplitudes],
        yticklabels=[f"{value:g}" for value in widths],
    )
    axes[2].tick_params(axis="x", rotation=45)
    plt.colorbar(continuous, ax=axes[2], fraction=0.046, pad=0.03,
                 label="Final response centroid (x/L)")
    _panel(axes[2], "C")
    figure.subplots_adjust(top=0.96, bottom=0.18, left=0.07, right=0.97)
    figure.savefig(output / "FigureS1_Otsuji_reversal_sweep.png", dpi=300)
    figure.savefig(output / "FigureS1_Otsuji_reversal_sweep.pdf")
    plt.close(figure)


def build_supplement_tracking(output: Path, source: Path) -> None:
    manifest = json.loads((ROOT / "campaigns" / "wp1_tracking_supplement_v1.json").read_text())
    tasks = sorted(manifest["tasks"], key=lambda task: task["protocol"]["period"])
    figure, axes = plt.subplots(3, 2, figsize=(10.5, 9.0), gridspec_kw={"hspace": 0.42, "wspace": 0.28})
    snapshot_rows, waveform_rows = [], []
    phase_fractions = np.linspace(0, 0.5, 6)
    colors = plt.cm.plasma(np.linspace(0.05, 0.90, len(phase_fractions)))
    for column, task in enumerate(tasks):
        period = float(task["protocol"]["period"])
        directory = TRACKING_SUPPLEMENT / "tasks" / task["task_id"]
        payload = _load_campaign_trajectory(TRACKING_SUPPLEMENT, task["task_id"])
        observable = _observable(payload)
        relative_x = payload["x"] / payload["x"][-1]
        cycle_start = float(payload["time"][-1]) - period
        sample_times = cycle_start + phase_fractions * period
        response_snapshots = np.asarray([
            [np.interp(time, payload["time"], observable[:, j]) for j in range(observable.shape[1])]
            for time in sample_times
        ])
        response_min, response_max = float(np.min(response_snapshots)), float(np.max(response_snapshots))
        response_scale = max(response_max - response_min, 1e-12)
        centered = 2 * relative_x - 1
        slope = float(task["protocol"]["slope"])
        length = float(task["simulation"]["length"])
        for color, fraction, time, response in zip(colors, phase_fractions, sample_times, response_snapshots):
            stimulus = 0.5 * slope * length * (1 - np.cos(2 * np.pi * fraction) * centered)
            label = f"{fraction:g}T"
            axes[0, column].plot(relative_x, stimulus, color=color, lw=1.4, label=label)
            axes[1, column].plot(relative_x, (response - response_min) / response_scale,
                                 color=color, lw=1.4, label=label)
            for x, stimulus_value, response_value in zip(relative_x, stimulus, response):
                snapshot_rows.append(
                    {"period_s": period, "phase_fraction": fraction, "time_s": time,
                     "relative_position": float(x), "stimulus": float(stimulus_value),
                     "response": float(response_value)}
                )
        axes[0, column].set(xlabel="Rel. pos.", ylabel="Stimulus (s⁻¹)" if column == 0 else "")
        axes[1, column].set(xlabel="Rel. pos.", ylabel="Scaled response" if column == 0 else "")
        axes[0, column].text(0.03, 0.92, f"T={period:g} s", transform=axes[0, column].transAxes,
                             va="top", fontsize=8)
        axes[0, column].legend(frameon=False, ncol=2, fontsize=6.2, loc="upper center")

        left = np.mean(observable[:, relative_x <= 0.2], axis=1)
        right = np.mean(observable[:, relative_x >= 0.8], axis=1)
        orientation = (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)
        scoring = payload["time"] >= float(task["tracking"]["start"])
        phase = 2 * np.pi * (payload["time"][scoring] - float(task["tracking"]["start"])) / period
        centered_response = orientation[scoring] - np.mean(orientation[scoring])
        cosine = 2 * np.mean(centered_response * np.cos(phase))
        sine = 2 * np.mean(centered_response * np.sin(phase))
        amplitude = float(np.hypot(cosine, sine))
        fit_all = np.mean(orientation[scoring]) + cosine * np.cos(phase) + sine * np.sin(phase)
        show = payload["time"][scoring] >= payload["time"][scoring][-1] - 2 * period
        relative_time = (payload["time"][scoring][show] - payload["time"][scoring][show][0]) / period
        target = -np.cos(phase[show])
        axes[2, column].plot(relative_time, target, color="#777777", ls="--", lw=1.0, label="Cue orientation")
        axes[2, column].plot(relative_time, orientation[scoring][show], color="#2a9d8f", lw=1.3,
                             label="Response orientation")
        axes[2, column].plot(relative_time, fit_all[show], color="#d62728", lw=1.2,
                             label="First harmonic")
        axes[2, column].set(xlabel="Time / T", ylabel="Orientation" if column == 0 else "", ylim=(-1.08, 1.08))
        axes[2, column].text(0.03, 0.08, f"Amp.={amplitude:.2f}", transform=axes[2, column].transAxes,
                             fontsize=7.5)
        axes[2, column].legend(frameon=False, fontsize=6.2, loc="upper right")
        for time_value, orientation_value, target_value, fit_value in zip(
            payload["time"][scoring][show], orientation[scoring][show], target, fit_all[show]
        ):
            waveform_rows.append(
                {"period_s": period, "time_s": float(time_value),
                 "time_over_period": float((time_value - payload["time"][scoring][show][0]) / period),
                 "response_orientation": float(orientation_value),
                 "cue_orientation": float(target_value), "first_harmonic": float(fit_value),
                 "response_amplitude": amplitude}
            )
    for label, axis in zip("ABCDEF", axes.flat):
        _panel(axis, label)
    _write(source / "figureS2_tracking_snapshots.csv", snapshot_rows)
    _write(source / "figureS2_tracking_waveforms.csv", waveform_rows)
    figure.subplots_adjust(top=0.97, bottom=0.07, left=0.08, right=0.97)
    figure.savefig(output / "FigureS2_tracking_examples.png", dpi=300)
    figure.savefig(output / "FigureS2_tracking_examples.pdf")
    plt.close(figure)


def _figure2_schematic(axis) -> None:
    axis.axis("off")
    _panel(axis, "A")
    axis.text(0.0, 0.98, "De Belly mechanochemical coordination", transform=axis.transAxes,
              fontweight="bold", va="top", fontsize=10)
    _box(axis, (0.03, 0.60), 0.18, 0.16, "Local Rac", "#cdece6")
    _box(axis, (0.29, 0.60), 0.21, 0.16, "Protrusion and\nmembrane tension", "#f9d7a8")
    _box(axis, (0.59, 0.60), 0.20, 0.16, "Distal Rho", "#d8ddf2")
    _arrow(axis, (0.21, 0.68), (0.29, 0.68))
    _arrow(axis, (0.50, 0.68), (0.59, 0.68))
    _box(axis, (0.59, 0.26), 0.20, 0.16, "Contractile flow", "#d8ddf2")
    _box(axis, (0.29, 0.26), 0.21, 0.16, "MCA redistribution", "#eadcf2")
    _box(axis, (0.03, 0.26), 0.18, 0.16, "Distal Rac", "#cdece6")
    _arrow(axis, (0.69, 0.60), (0.69, 0.42))
    _arrow(axis, (0.59, 0.34), (0.50, 0.34))
    _arrow(axis, (0.29, 0.34), (0.21, 0.34))
    axis.text(0.88, 0.66, "Rac→mechanics→Rho", transform=axis.transAxes,
              ha="center", color="#8c5b17", fontweight="bold", fontsize=8)
    axis.text(0.88, 0.34, "Rho→mechanics→Rac", transform=axis.transAxes,
              ha="center", color="#55457d", fontweight="bold", fontsize=8)
    axis.text(
        0.50, 0.49, "Local Rac ⇄ Rho mutual inhibition and diffusion remain when mechanical arms are off",
        transform=axis.transAxes, ha="center", fontsize=7.2, color="#333333",
        bbox={"facecolor": "white", "edgecolor": "#999999", "alpha": 0.9, "pad": 2},
    )
    axis.text(
        0.5, 0.06,
        "Question: does the size-dependent coordination burden remain when long-range mechanical signaling is present?",
        transform=axis.transAxes, ha="center", fontweight="bold", fontsize=8.2,
    )


def _intact_size_burden(axis, competition: list[dict[str, str]]) -> None:
    selected = [row for row in competition if row["coupling_variant"] == "intact"]
    lengths = [_number(row, "length_um") for row in selected]
    latency = [_number(row, "competition_resolution_latency_s") for row in selected]
    axis.plot(lengths, latency, "o-", color="#4b176f", lw=2.2, ms=4.5)
    for row in selected:
        persistent = _boolean(row, "persistent_resolution")
        axis.scatter(
            _number(row, "length_um"), _number(row, "competition_resolution_latency_s"),
            s=50, edgecolor="#4b176f", facecolor="#4b176f" if persistent else "white",
            zorder=3,
        )
    axis.set_yscale("log")
    axis.set_xlabel("Cell length (µm)")
    axis.set_ylabel("Resolution latency (s)")
    axis.grid(alpha=0.18)
    axis.set_title("The burden remains with intact mechanics", loc="left")
    axis.text(
        0.04, 0.95, "240 s → 2,956 s\n12.3-fold increase",
        transform=axis.transAxes, va="top", fontsize=8.5, fontweight="bold",
    )
    axis.text(
        0.04, 0.69, "Filled: persistent final front\nOpen: transient resolution",
        transform=axis.transAxes, va="top", fontsize=7.3, color="#444444",
    )
    inset = axis.inset_axes([0.56, 0.18, 0.40, 0.18])
    relative = np.linspace(0, 1, 301)
    sigma = 1.6 / 40.0
    cue = 1.01 * np.exp(-0.5 * (relative / sigma) ** 2) + 0.99 * np.exp(
        -0.5 * ((relative - 1) / sigma) ** 2
    )
    inset.fill_between(relative, 0, cue, color="#999999", alpha=0.45)
    inset.set(xlim=(0, 1), ylim=(0, 1.08), xticks=(0, 1), yticks=())
    inset.set_title("End cues active 50–150 s", fontsize=6.2, pad=1)
    _panel(axis, "B")


def _coupling_outcome_grid(axis, rows: list[dict[str, str]]) -> None:
    lengths = sorted({_number(row, "length_um") for row in rows})
    categories = (
        "resolved_persistent",
        "resolved_transient_to_polarized_no_resolved_peak",
        "resolved_transient_to_unpolarized",
        "unresolved_persistent_multifront",
        "unresolved_to_unpolarized",
    )
    colors = ("#2f7ebc", "#68a9cf", "#b9d8ea", "#d95f45", "#d9d9d9")
    matrix = np.zeros((len(MAIN_VARIANTS), len(lengths)), dtype=int)
    letters = np.empty_like(matrix, dtype=object)
    for row in rows:
        if row["coupling_variant"] not in MAIN_VARIANTS:
            continue
        i = MAIN_VARIANTS.index(row["coupling_variant"])
        j = lengths.index(_number(row, "length_um"))
        matrix[i, j] = categories.index(row["competition_effective_status"])
        letters[i, j] = {
            "unpolarized": "U", "single_front": "S", "multifront": "M",
            "polarized_no_resolved_peak": "P",
        }[row["endpoint_state"]]
    axis.imshow(
        matrix, aspect="auto", cmap=ListedColormap(colors),
        norm=BoundaryNorm(np.arange(-0.5, len(categories) + 0.5), len(categories)),
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            axis.text(j, i, letters[i, j], ha="center", va="center", fontsize=5.8,
                      fontweight="bold")
    axis.set_yticks(range(len(MAIN_VARIANTS)), [VARIANT_LABELS[value] for value in MAIN_VARIANTS])
    shown = list(range(0, len(lengths), 2))
    axis.set_xticks(shown, [f"{lengths[index]:g}" for index in shown])
    axis.set_xlabel("Cell length (µm)")
    axis.set_title("Mechanical-arm ablation changes the failure mode", loc="left")
    axis.legend(
        handles=[Patch(color=color, label=label) for color, label in zip(
            colors,
            ("Persistent resolution", "Transient → broad polarized", "Transient → unpolarized",
             "Persistent multifront", "Unresolved → unpolarized"),
        )],
        frameon=False, fontsize=6.5, ncol=2, loc="upper left", bbox_to_anchor=(0, -0.24),
    )
    _panel(axis, "C")


def _intact_reversal(axis, rows: list[dict[str, str]]) -> None:
    for input_peak, color, marker in ((0.2, "#4b176f", "o"), (0.5, "#d19a00", "s")):
        selected = [
            row for row in rows
            if row["coupling_variant"] == "intact"
            and np.isclose(_number(row, "input_peak"), input_peak)
        ]
        axis.plot(
            [_number(row, "length_um") for row in selected],
            [_number(row, "orientation_reversal_latency_s") for row in selected],
            marker + "-", color=color, lw=2.0, ms=3.8, label=f"Cue peak {input_peak:g}",
        )
    axis.set(xlabel="Cell length (µm)", ylabel="Side-switch latency (s)")
    axis.set_title("Intact model: weak-cue size effect is non-monotonic", loc="left")
    axis.grid(alpha=0.18)
    axis.legend(frameon=False)
    axis.text(
        0.04, 0.08, "All intact cells switch; weak-cue latency\nrises again from 48–64 µm",
        transform=axis.transAxes, fontsize=7.2,
    )
    _panel(axis, "D")


def _full_ablation_reversal(axis, rows: list[dict[str, str]]) -> None:
    selected_input = [row for row in rows if np.isclose(_number(row, "input_peak"), 0.2)]
    marker_height = 590.0
    for variant in MAIN_VARIANTS:
        selected = [row for row in selected_input if row["coupling_variant"] == variant]
        successful = [row for row in selected if row["orientation_reversal_status"] == "success"]
        axis.plot(
            [_number(row, "length_um") for row in successful],
            [_number(row, "orientation_reversal_latency_s") for row in successful],
            "o-", color=VARIANT_COLORS[variant], lw=1.7, ms=3.3, label=VARIANT_LABELS[variant],
        )
        censored = [row for row in selected if row["orientation_reversal_status"] == "right_censored"]
        axis.scatter(
            [_number(row, "length_um") for row in censored], [marker_height] * len(censored),
            marker="x", color=VARIANT_COLORS[variant], s=27,
        )
    axis.set(xlabel="Cell length (µm)", ylabel="Side-switch latency (s)", ylim=(65, 625))
    axis.set_title("Full-arm ablations under the persistent weak cue", loc="left")
    axis.grid(alpha=0.18)
    axis.legend(frameon=False, fontsize=6.8, ncol=2)
    axis.text(
        0.98, 0.43, "Both off follows the reversed cue at every length\n(140 → 126 s)",
        transform=axis.transAxes, ha="right", fontsize=7.2, color=VARIANT_COLORS["both_off"],
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )
    axis.text(0.98, 0.93, "× not switched by run end", transform=axis.transAxes,
              ha="right", va="top", fontsize=6.6, color="#555555")
    _panel(axis, "E")


def _reversal_interaction(axis, rows: list[dict[str, str]], input_peak: float) -> None:
    selected_input = [row for row in rows if np.isclose(_number(row, "input_peak"), input_peak)]
    successful_values = [
        _number(row, "orientation_reversal_latency_s") for row in selected_input
        if _number(row, "orientation_reversal_latency_s") is not None
    ]
    marker_height = 1.08 * max(successful_values)
    for variant in VARIANTS:
        selected = [row for row in selected_input if row["coupling_variant"] == variant]
        success = [row for row in selected if row["orientation_reversal_status"] == "success"]
        axis.plot(
            [_number(row, "length_um") for row in success],
            [_number(row, "orientation_reversal_latency_s") for row in success],
            "o-", color=VARIANT_COLORS[variant], lw=1.5, ms=3.2,
            label=VARIANT_LABELS[variant],
        )
        censored = [row for row in selected if row["orientation_reversal_status"] == "right_censored"]
        axis.scatter(
            [_number(row, "length_um") for row in censored],
            [marker_height] * len(censored), marker="x", color=VARIANT_COLORS[variant], s=28,
        )
    axis.set_xlabel("Cell length (µm)")
    axis.set_ylabel("Side-switch latency (s)")
    axis.grid(alpha=0.18)
    axis.set_title(f"Cue peak {input_peak:g}" + (" reveals size-dependent failure" if input_peak == 0.2 else " rescues switching"), loc="left")
    if input_peak == 0.2:
        axis.legend(frameon=False, ncol=2, fontsize=6.7)
        _panel(axis, "D")
    else:
        _panel(axis, "E")


def _manifest_lookup() -> dict[str, dict]:
    manifest = json.loads(MANIFEST.read_text())
    return {task["task_id"]: task for task in manifest["tasks"]}


def _load_trajectory(task_id: str) -> dict[str, np.ndarray]:
    path = DECISIVE / "tasks" / task_id / "trajectory.npz"
    with np.load(path, allow_pickle=False) as payload:
        return {key: payload[key] for key in payload.files}


def _orientation_trace(payload: dict[str, np.ndarray]) -> np.ndarray:
    x = payload["x"]
    rac = payload["field__rac"]
    relative = (x - x[0]) / (x[-1] - x[0])
    left = np.mean(rac[:, relative <= 0.2], axis=1)
    right = np.mean(rac[:, relative >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _representative_reversal(axis, rows: list[dict[str, str]], source: Path) -> None:
    requested = (
        ("intact", 40.0, "Intact, 40 µm", "#4b176f"),
        ("rac_to_rho_off", 25.198420997897465, "Rac→Rho off, 25.2 µm", "#2ca58d"),
        ("rac_to_rho_off", 40.0, "Rac→Rho off, 40 µm", "#d95f45"),
    )
    lookup = _manifest_lookup()
    exported = []
    for variant, length, label, color in requested:
        row = next(
            item for item in rows
            if item["coupling_variant"] == variant
            and np.isclose(_number(item, "input_peak"), 0.2)
            and np.isclose(_number(item, "length_um"), length)
        )
        task_id = row["task_id"]
        payload = _load_trajectory(task_id)
        task = lookup[task_id]
        reversal_time = float(task["protocol"]["reversal_time"])
        scaled_time = (payload["time"] - reversal_time) * (20.0 / length) ** 2
        orientation = _orientation_trace(payload)
        mask = (scaled_time >= -250) & (scaled_time <= 750)
        axis.plot(scaled_time[mask], orientation[mask], color=color, lw=1.8, label=label)
        for time, value in zip(scaled_time[mask], orientation[mask]):
            exported.append(
                {
                    "condition": label,
                    "task_id": task_id,
                    "length_um": length,
                    "scaled_time_s_at_L20": float(time),
                    "orientation": float(value),
                }
            )
    _write(source / "figure2_representative_reversal_traces.csv", exported)
    axis.axvline(0, color="black", lw=0.8, ls="--")
    axis.axhline(0, color="#777777", lw=0.7)
    axis.axhspan(0.05, 1, color="#e6f3eb", alpha=0.45)
    axis.set_xlabel("Time from cue reversal, scaled to L=20 (s)")
    axis.set_ylabel("Rac side orientation")
    axis.set_ylim(-0.55, 0.55)
    axis.grid(alpha=0.15)
    axis.set_title("Representative weak-cue polarity replacement", loc="left")
    axis.legend(frameon=False, fontsize=6.8)
    _panel(axis, "F")


def _chemical_field_examples(axes, rows: list[dict[str, str]], source: Path) -> None:
    conditions = (
        ("intact", "Intact: switches"),
        ("rac_to_rho_off", "Rac→Rho off: fails"),
        ("both_off", "Both off: cue-following"),
    )
    fields = (
        ("rac", "Rac"), ("rho", "Rho"),
        ("mca_bound", "MCA, bound"), ("mca_unbound", "MCA, unbound"),
    )
    offsets = (-20.0, 50.0, 150.0, 500.0)
    colors = plt.cm.plasma(np.linspace(0.08, 0.88, len(offsets)))
    lookup = _manifest_lookup()
    exported = []
    for row_index, (variant, condition_label) in enumerate(conditions):
        row = next(
            item for item in rows
            if item["coupling_variant"] == variant
            and np.isclose(_number(item, "input_peak"), 0.2)
            and np.isclose(_number(item, "length_um"), 40.0)
        )
        payload = _load_trajectory(row["task_id"])
        reversal_time = float(lookup[row["task_id"]]["protocol"]["reversal_time"])
        relative_x = payload["x"] / payload["x"][-1]
        for column_index, (field, field_label) in enumerate(fields):
            axis = axes[row_index][column_index]
            values = payload[f"field__{field}"]
            for offset, color in zip(offsets, colors):
                time_index = int(np.argmin(np.abs(payload["time"] - (reversal_time + offset))))
                actual_offset = float(payload["time"][time_index] - reversal_time)
                axis.plot(relative_x, values[time_index], color=color, lw=1.25,
                          label=f"{actual_offset:+g} s")
                for x, value in zip(relative_x, values[time_index]):
                    exported.append(
                        {"condition": condition_label, "coupling_variant": variant,
                         "task_id": row["task_id"], "field": field,
                         "time_from_reversal_s": actual_offset,
                         "relative_position": float(x), "value": float(value)}
                    )
            axis.set_xlim(0, 1)
            axis.set_xticks((0, 0.5, 1))
            if row_index == 0:
                axis.set_title(field_label, fontsize=8.2)
            if column_index == 0:
                axis.set_ylabel(condition_label, fontsize=7.1, rotation=0, ha="right", va="center",
                                labelpad=42)
            else:
                axis.set_yticklabels([])
            if row_index < len(conditions) - 1:
                axis.set_xticklabels([])
            else:
                axis.set_xlabel("Relative position")
    axes[0][0].text(
        0.0, 1.34, "Chemical fields around weak-cue reversal (L = 40 µm)",
        transform=axes[0][0].transAxes, fontsize=10, fontweight="bold", ha="left",
    )
    _panel(axes[0][0], "F")
    axes[0][-1].legend(frameon=False, fontsize=6.4, title="Time from switch", title_fontsize=6.5,
                       loc="upper left", bbox_to_anchor=(1.02, 1.04))
    _write(source / "figure2_chemical_field_examples.csv", exported)


def _debelly_reversal_field_summary(axis, source: Path) -> None:
    rows = _read(REVISION_AUDIT / "analysis" / "debelly_field_orientation.csv")
    selected = [row for row in rows if np.isclose(_number(row, "input_peak"), 0.5)]
    styles = {"rac": ("Rac", "#6a1b9a"), "rho": ("Rho", "#2f7ebc")}
    for field, (label, color) in styles.items():
        field_rows = [row for row in selected if row["field"] == field]
        axis.plot(
            [_number(row, "length_um") for row in field_rows],
            [_number(row, "post_switch_event_time_s") for row in field_rows],
            "o-", lw=1.9, ms=3.8, color=color, label=f"{label} reorientation",
        )
    axis.axvspan(40, 50, color="#efc050", alpha=0.15, label="latency transition")
    axis.set(xlabel="Cell length (µm)", ylabel="Field-reorientation latency (s)")
    axis.set_title("Cue peak 0.5: chemical reorientation slows near 40–50 µm", loc="left")
    axis.grid(alpha=0.18)
    axis.legend(frameon=False)
    axis.text(
        0.98, 0.93, "MCA never reverses at any tested length", transform=axis.transAxes,
        ha="right", va="top", fontsize=7.6, fontweight="bold", color="#b23a48",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.5},
    )
    axis.text(-0.16, 1.08, "D", transform=axis.transAxes, fontsize=13,
              fontweight="bold", va="top", ha="left")
    _write(source / "figure2_debelly_field_orientation.csv", rows)


def _debelly_reversal_cue(axis) -> None:
    position = np.linspace(0, 1, 201)
    times = np.linspace(-200, 500, 351)
    matrix = np.empty((len(times), len(position)))
    for index, time in enumerate(times):
        if time < 0:
            matrix[index] = 0.5 * (1 - position)
        elif time >= 25:
            matrix[index] = 0.5 * position
        else:
            weight = time / 25
            matrix[index] = 0.5 * ((1 - weight) * (1 - position) + weight * position)
    axis.imshow(matrix.T, origin="lower", aspect="auto", extent=(-200, 500, 0, 1),
                cmap="Greys", vmin=0, vmax=0.5)
    axis.set(xlabel="Time from cue reversal (s)", ylabel="Relative position")
    axis.set_title("Persistent Rac-activation gradient (peak 0.5)", loc="left")
    axis.text(
        0.03, 0.95, "25-s cross-fade; cue remains on throughout\n"
        "L40–52: condition 1,000–1,690 s; follow 3,000–5,070 s\n"
        "Fixed mechanical leading edge remains at x=0",
        transform=axis.transAxes, va="top", fontsize=7,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 1.5},
    )
    _panel(axis, "E")


def _debelly_transition_profiles(axes, source: Path) -> None:
    summary = _read(REVISION_AUDIT / "task_summary.csv")
    manifest = json.loads((ROOT / "campaigns" / "wp1_revision_audit_v1.json").read_text())
    lookup = {task["task_id"]: task for task in manifest["tasks"]}
    lengths = (40.0, 48.0, 52.0)
    fields = (("rac", "Rac"), ("rho", "Rho"), ("mca_bound", "MCA, bound"), ("mca_unbound", "MCA, unbound"))
    offsets = (-20.0, 50.0, 150.0, 500.0)
    colors = plt.cm.plasma(np.linspace(0.08, 0.88, len(offsets)))
    exported = []
    for row_index, length in enumerate(lengths):
        row = next(
            item for item in summary
            if item["assay"] == "debelly_intact_reversal_field_audit"
            and np.isclose(_number(item, "input_peak"), 0.5)
            and np.isclose(_number(item, "length_um"), length)
        )
        payload = _load_campaign_trajectory(REVISION_AUDIT, row["task_id"])
        reversal = float(lookup[row["task_id"]]["protocol"]["reversal_time"])
        relative_x = payload["x"] / payload["x"][-1]
        for column_index, (field, label) in enumerate(fields):
            axis = axes[row_index][column_index]
            for offset, color in zip(offsets, colors):
                index = int(np.argmin(np.abs(payload["time"] - (reversal + offset))))
                actual = float(payload["time"][index] - reversal)
                values = payload[f"field__{field}"][index]
                axis.plot(relative_x, values, color=color, lw=1.25, label=f"{actual:+g} s")
                for x, value in zip(relative_x, values):
                    exported.append(
                        {"length_um": length, "input_peak": 0.5, "task_id": row["task_id"],
                         "field": field, "time_from_reversal_s": actual,
                         "relative_position": float(x), "value": float(value)}
                    )
            axis.set_xlim(0, 1)
            axis.set_xticks((0, 0.5, 1))
            if row_index == 0:
                axis.set_title(label, fontsize=8.2)
            if column_index == 0:
                axis.set_ylabel(f"L = {length:g} µm", rotation=0, ha="right", va="center", labelpad=32)
            else:
                axis.set_yticklabels([])
            if row_index < len(lengths) - 1:
                axis.set_xticklabels([])
            else:
                axis.set_xlabel("Relative position")
    axes[0][0].text(
        0.0, 1.34, "Fields across the 0.5-cue latency transition", transform=axes[0][0].transAxes,
        fontsize=10, fontweight="bold", ha="left",
    )
    axes[0][0].text(-0.18, 1.06, "F", transform=axes[0][0].transAxes, fontsize=13,
                    fontweight="bold", va="top", ha="left")
    axes[0][-1].legend(frameon=False, fontsize=6.4, title="Time from reversal", title_fontsize=6.5,
                       loc="upper left", bbox_to_anchor=(1.02, 1.04))
    _write(source / "figure2_debelly_transition_profiles.csv", exported)


def build_figure2(output: Path, source: Path) -> None:
    competition = _read(DECISIVE / "competition_hardened.csv")
    reversal = _read(DECISIVE / "reversal_length_scaled.csv")
    propagation = _read(DECISIVE / "propagation_fixed_width.csv")
    _write(source / "figure2_competition_hardened.csv", competition)
    _write(source / "figure2_reversal_length_scaled.csv", reversal)
    _write(source / "figure2_fixed_width_propagation.csv", propagation)

    figure = plt.figure(figsize=(13.2, 13.2))
    grid = figure.add_gridspec(4, 2, height_ratios=(0.67, 1.0, 0.9, 1.35), hspace=0.67, wspace=0.32)
    schematic = figure.add_subplot(grid[0, :])
    intact = figure.add_subplot(grid[1, 0])
    outcomes = figure.add_subplot(grid[1, 1])
    field_summary = figure.add_subplot(grid[2, 0])
    reversal_cue = figure.add_subplot(grid[2, 1])
    field_grid = grid[3, :].subgridspec(3, 4, hspace=0.12, wspace=0.12)
    field_axes = [[figure.add_subplot(field_grid[i, j]) for j in range(4)] for i in range(3)]
    _figure2_schematic(schematic)
    _intact_size_burden(intact, competition)
    _coupling_outcome_grid(outcomes, competition)
    _debelly_reversal_field_summary(field_summary, source)
    _debelly_reversal_cue(reversal_cue)
    _debelly_transition_profiles(field_axes, source)
    figure.suptitle(
        "Figure 2. Mechanochemical signaling modifies—but does not eliminate—the size-dependent coordination burden",
        fontsize=14, fontweight="bold", y=0.995,
    )
    figure.text(
        0.01, 0.006,
        "C: letters give the final Rac state (U, unpolarized; S, single front; M, multifront; P, polarized without a resolved peak). "
        "Mechanical-arm ablations retain local Rac–Rho mutual inhibition and diffusion. D–F diagnose chemical redistribution on the released model's fixed-front geometry; they are not scored as whole-cell polarity reversal.",
        fontsize=7.2, color="#444444",
    )
    figure.savefig(output / "Figure2_WP1_mechanochemical.png", dpi=300)
    figure.savefig(output / "Figure2_WP1_mechanochemical.pdf")
    plt.close(figure)


def _write_text_package(output: Path) -> None:
    captions = """# WP1 manuscript figure captions

## Figure 1. Increasing cell size creates a polarity-coordination burden across established model classes

(A) Spatiotemporal input for the two-cue competition assay. Slightly asymmetric cues stimulate the outer 10% of each end, remain at full strength through 20 s, and taper to zero at 40 s; the red dashed line marks 20 s. (B) Representative L=15 µm response profiles at 20, 200 and 20,000 s. (C) Competition-resolution time across cell length. The vertical dotted line identifies L=15 µm. All models use a common 20,000-s observation horizon; red horizontal dashed lines mark the 200- and 20,000-s profile times. Low-alpha model-colored regions show the length–time area in which competition has resolved, and unresolved runs are not plotted as censoring symbols. LEGI is an input-sensing circuit and is not interpreted as intrinsic front competition. (D) Persistent-gradient reversal protocol; red dashed lines identify 0⁻, 20 and 100 s relative to reversal. (E) Representative L=15 µm responses immediately before reversal (0⁻) and at 20 and 100 s after reversal. (F) Reversal time across cell length, with red horizontal guides at 20 and 100 s and a black vertical L=15 µm guide; otherwise it uses the same common-horizon and resolved-region convention as C. (G) One half-period of the alternating-gradient input. Red dashed lines mark π/8, 2π/8 and 3π/8. (H) Representative L=15 µm, T=200 s response profiles at those phases during the final simulated cycle. (I) Orientation-response amplitude and phase lag (fraction of a period), shown as separate metric columns and faceted by model with cell length on the ordinate. Gray cells with crosses fail the waveform audit and are excluded from both color scales; the audit retains the amplitude and lag criteria and additionally requires first-harmonic R² ≥ 0.8, cycle-amplitude CV ≤ 0.25 and circular phase dispersion ≤ 0.1. Red edges identify the L=15 µm, T=200 s conditions shown in H. The assays distinguish basic cue responsiveness from the cell-wide burden of selecting, replacing and dynamically coordinating polarity. Numerical transition lengths remain protocol-specific properties of qualitative equation-level reconstructions.

## Figure S1. Otsuji reversal across stimulus amplitude and spatial width

(A) Finite-width linear edge gradients used before (solid) and after (dashed) reversal. Width/L=1 recovers the full-cell linear gradient used in Figure 1; smaller values truncate the ramp nearer the stimulated edge. (B) Change in the centroid of active-field mass above the profile minimum at L=15 µm after 2,000 s conditioning and 20,000 s follow-up. None of the 64 amplitude–width combinations met the strict reversal criterion, so the continuous redistribution endpoint is shown. (C) Final response centroid, where zero and one denote left- and right-edge localization.

## Figure S2. Wave-pinning examples expose the tracking-amplitude calculation

(A,B) Six stimulus profiles spanning one half-period at L=15 µm for T=100 and 400 s. (C,D) Corresponding wave-pinning response profiles during the final simulated cycle. (E,F) Normalized left–right response orientation over the final two cycles, the signed cue orientation, and the fitted first harmonic. Panel I reports the magnitude of the fitted cosine and sine coefficients, shown here as Amp.

## Figure 2. Mechanochemical signaling modifies—but does not eliminate—the size-dependent coordination burden

(A) De Belly model architecture. Rac-dependent protrusion and tension-dependent Rho form one long-range mechanical arm; Rho-driven cortical flow, MCA redistribution and MCA-dependent Rac form the reciprocal arm. The perturbation multipliers remove these mechanical paths only: local Rac–Rho mutual inhibition, basal activation and Rac/Rho diffusion remain intact. The released 1D geometry has a fixed Rac-responsive leading edge at x=0 and an anchored boundary at x=L. (B) Under a fixed 1.6-µm two-cue width, intact-model resolution latency increases from 240 s at 12 µm to 2,956 s at 64 µm. Filled symbols end as a sustained single front; open symbols only pass transiently through a one-front state. Thus the size-dependent coordination burden remains present with intact mechanochemical signaling. (C) Long-time Rac outcome in the intact model and three complete mechanical-arm ablations; 50% conditions are omitted for clarity. Mechanical signaling changes whether the burden manifests as depolarization, broad polarization or persistent competing fronts, but the perturbations are not an ordered dose response. (D) Under a persistent reversing Rac-activation gradient of peak 0.5, Rac and Rho reorientation slow sharply around 40–50 µm. MCA does not reverse at any tested length. (E) The gradient cross-fades over 25 s and then remains on. Because the model's mechanical leading edge remains fixed at x=0, this assay tests chemical redistribution on a fixed-front reference cell rather than whole-cell polarity reversal. (F) Rac, Rho, bound MCA and unbound MCA profiles at 40, 48 and 52 µm around reversal. Rac and Rho follow the reversed input, while the MCA distribution remains locked toward the original mechanical front. The former Rac-only “switch” endpoint is therefore insufficient as a mechanochemical polarity criterion and is not used to support the principal WP1 conclusion.
"""
    readme = """# Near-manuscript-ready WP1 figure package

This package is generated from the frozen Jilkine reconstruction, fixed-horizon competition and reversal audits, the Otsuji amplitude-by-width reversal sweep, and focused tracking examples. It contains raster previews, vector PDFs, manuscript captions, and panel-level source data.

The intended result hierarchy is:

1. Increasing cell size creates a polarity-coordination burden across established biochemical model classes.
2. The same burden remains present in the intact De Belly mechanochemical model.
3. Mechanical-arm perturbations change the severity and failure mode of the burden rather than producing a simple ordered buffering response.
4. The fixed-front De Belly reversal assay supports chemical Rac/Rho redistribution but not full mechanochemical polarity reversal; it is treated as a diagnostic rather than headline evidence.

The figures are near-manuscript-ready drafts. They deliberately label the Jilkine implementations as qualitative equation-level reconstructions and do not interpret deterministic condition counts as biological probabilities.
"""
    (output / "CAPTIONS.md").write_text(captions)
    (output / "README.md").write_text(readme)
    cue_audit = """# Otsuji cue audit

The source protocol defines a transient localized cue over the outer 10% of a 10-µm domain and reports a broad tested amplitude range (`S=0.01–0.5`), with a plateau ending at 10–25 s and decay ending at 15–50 s. It does not identify a unique Otsuji amplitude for the size comparison.

The original production campaign selected `S=0.2`, a plateau through 20 s, decay to zero at 40 s, and a 2% left-right asymmetry as one shared in-range protocol for all four reconstructed models. This was a cross-model normalization choice, not an Otsuji-specific value mandated by the source. The earlier figure-guided Otsuji qualification instead used `S=0.04` with the same timing and 2% two-cue asymmetry.

`otsuji_single_cue_v1` therefore tested one left-edge cue at both `S=0.04` and `S=0.2`, without noise, across odd lengths from 1–39 µm. Both amplitudes produced the same qualitative result: the 1-µm domain remained unpolarized, while every domain from 3–39 µm formed a sustained cue-aligned single front. The apparent >13-µm failure in the reversal assay and the large-cell two-cue outcomes must consequently be described as failures of polarity coordination or replacement under those protocols, not failures of cue-induced polarization.
"""
    (output / "OTSUJI_CUE_AUDIT.md").write_text(cue_audit)
    debelly_audit = """# De Belly reversal and ablation audit

The prior `orientation_reversal_success` endpoint was based only on Rac: a polarized Rac profile had to change its left-right orientation sign for a sustained interval. It did not require opposite Rho orientation, MCA redistribution, tension reversal, or relocation of the protrusive boundary.

The released 1D model and its optimized port use an asymmetric fixed-reference geometry. Rac at the left boundary (`rac[0]`) controls protrusive boundary velocity, while the right boundary is anchored. Reversing a spatial Rac-activation gradient therefore does not exchange the model's mechanical front and rear. In the field audit, Rac and Rho reoriented under cue peak 0.5, but bound and unbound MCA remained biased toward the original mechanical front at every tested length. These runs are now described as chemical cue-following, not whole-cell mechanochemical polarity reversal.

The campaign multipliers are mechanical-arm perturbations. `rac_to_rho_coupling=0` removes Rac-driven protrusion and tension-dependent Rho activation. `rho_to_rac_coupling=0` removes Rho active-stress forcing and MCA-dependent Rac activation. Setting both to zero leaves basal Rac/Rho activation, local mutual inhibition through the reaction denominators, Rac/Rho diffusion, and direct optogenetic Rac activation intact. Consequently, a “both mechanical arms off” cell can follow a persistent Rac cue; this does not represent complete Rac–Rho uncoupling.
"""
    (output / "DEBELLY_REVERSAL_ABLATION_AUDIT.md").write_text(debelly_audit)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args(argv)
    output = arguments.output
    source = output / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    for stale in (
        output / "FigureS1_Otsuji_single_cue.png",
        output / "FigureS1_Otsuji_single_cue.pdf",
        source / "figureS1_otsuji_single_cue.csv",
        source / "figure1_otsuji_single_cue_profiles.csv",
        source / "otsuji_single_cue_collection_summary.json",
    ):
        if stale.exists():
            stale.unlink()
    stale_trace = source / "figure2_representative_reversal_traces.csv"
    if stale_trace.exists():
        stale_trace.unlink()
    _style()
    build_figure1(output, source)
    build_supplement_otsuji(output, source)
    build_supplement_tracking(output, source)
    _write_text_package(output)
    shutil.copy2(REVISION_AUDIT / "collection_summary.json", source / "wp1_revision_audit_collection_summary.json")
    shutil.copy2(LONG_REVERSAL / "collection_summary.json", source / "jilkine_long_reversal_collection_summary.json")
    shutil.copy2(FIGURE1_COMMON / "collection_summary.json", source / "wp1_figure1_common_horizon_collection_summary.json")
    shutil.copy2(OTSUJI_GEOMETRY / "collection_summary.json", source / "otsuji_reversal_geometry_collection_summary.json")
    shutil.copy2(TRACKING_SUPPLEMENT / "collection_summary.json", source / "wp1_tracking_supplement_collection_summary.json")
    shutil.copy2(TRACKING_AUDIT / "collection_summary.json", source / "wp1_tracking_reliability_audit_collection_summary.json")
    shutil.copy2(TRACKING_AUDIT / "tracking_reliability_audit_summary.json", source / "wp1_tracking_reliability_audit_summary.json")
    shutil.copy2(TRACKING_AUDIT / "TRACKING_RELIABILITY_AUDIT.md", output / "TRACKING_RELIABILITY_AUDIT.md")
    # Figure 2 is intentionally rebuilt last from the 2026 de Belly campaign
    # under the Jilkine common protocol.  The retired fixed-width campaign must
    # not overwrite or contribute source data to the manuscript package.
    from build_debelly_2026_common_results import main as build_debelly_2026
    build_debelly_2026()
    print(output)


if __name__ == "__main__":
    main()
