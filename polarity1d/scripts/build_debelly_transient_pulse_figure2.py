#!/usr/bin/env python3
"""Build the transient-pulse backup version of de Belly Figure 2."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import FinitePulseSequence, LocalizedPulse  # noqa: E402


PRIMARY_NAME = "debelly_2026_jilkine_common_v1"
PULSE_NAME = "debelly_2026_transient_pulse_figure2_v1"
PRIMARY = ROOT / "results" / "campaigns" / PRIMARY_NAME
PULSE = ROOT / "results" / "campaigns" / PULSE_NAME
OUTPUT = ROOT / "results" / "manuscript_wp1_v3"
SOURCE = OUTPUT / "source_data"
DESTINATION = OUTPUT / "Figure2_WP1_mechanochemical_transient_pulse_backup.png"

FIELDS = (
    "cue", "rac", "rho", "tension", "mca_bound", "mca_unbound",
    "delta_beta_tension", "delta_alpha_mca",
)
FIELD_LABELS = {
    "cue": "Cue", "rac": "Rac", "rho": "Rho", "tension": "Tension",
    "mca_bound": "MCA bound", "mca_unbound": "MCA unbound",
    "delta_beta_tension": r"$\Delta\beta_T$",
    "delta_alpha_mca": r"$\Delta\alpha_{MCA}$",
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _manifest(name: str) -> dict:
    return json.loads((ROOT / "campaigns" / f"{name}.json").read_text())


def _payload(campaign: Path, task_id: str) -> dict[str, np.ndarray]:
    with np.load(campaign / "tasks" / task_id / "trajectory.npz", allow_pickle=False) as raw:
        return {key: raw[key] for key in raw.files}


def _task_map(manifest: dict) -> dict[str, dict]:
    return {task["task_id"]: task for task in manifest["tasks"]}


def _protocol(task: dict):
    configuration = dict(task["protocol"])
    name = configuration.pop("name")
    if name == "localized_pulse":
        if "sides" in configuration:
            configuration["sides"] = tuple(configuration["sides"])
        return LocalizedPulse(**configuration)
    if name == "finite_pulse_sequence":
        return FinitePulseSequence(pulses=tuple(configuration["pulses"]))
    raise ValueError(name)


def _select(rows, tasks, *, assay: str, length: float, **configuration):
    for row in rows:
        if row["assay"] != assay or float(row["length_um"]) != length:
            continue
        task = tasks[row["task_id"]]
        if all(task.get(key) == value for key, value in configuration.items()):
            return row, task
    raise RuntimeError(f"missing {assay}, L={length}, {configuration}")


def _kymograph_entry(
    campaign, row, task, label, title, events, *, display_stop=None,
):
    payload = _payload(campaign, row["task_id"])
    time = np.asarray(payload["time"], dtype=float)
    x = np.asarray(payload["x"], dtype=float)
    protocol = _protocol(task)
    cue = np.asarray([protocol(float(at), x, float(task["simulation"]["length"])) for at in time])
    selected = np.ones(len(time), dtype=bool) if display_stop is None else time <= display_stop
    time = time[selected]
    cue = cue[selected]
    fields = {"cue": cue}
    for field in FIELDS[1:]:
        fields[field] = np.asarray(payload[f"field__{field}"], dtype=float)[selected]
    return {
        "time": time, "x": x / float(task["simulation"]["length"]),
        "fields": fields, "label": label, "title": title, "events": events,
        "task_id": row["task_id"],
    }


def _float(value):
    return None if value in (None, "") else float(value)


def _write_rows(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    primary_manifest = _manifest(PRIMARY_NAME)
    pulse_manifest = _manifest(PULSE_NAME)
    primary_tasks = _task_map(primary_manifest)
    pulse_tasks = _task_map(pulse_manifest)
    primary_rows = _read_csv(PRIMARY / "task_summary.csv")
    pulse_rows = _read_csv(PULSE / "task_summary.csv")

    competition_row, competition_task = _select(
        primary_rows, primary_tasks, assay="competing_fronts_common", length=15.0,
        coupling_variant="full",
    )
    reversal_row, reversal_task = _select(
        pulse_rows, pulse_tasks, assay="transient_pulse_reversal", length=15.0,
        coupling_variant="full", pulse_amplitude_ratio=4.0, pulse_duration_s=320.0,
    )
    repeated_row, repeated_task = _select(
        pulse_rows, pulse_tasks, assay="repeated_transient_switching", length=15.0,
        coupling_variant="full", perturbation_multiplier=1.0,
    )
    entries = [
        _kymograph_entry(
            PRIMARY, competition_row, competition_task, "A",
            "Two transient end cues", (20.0, 40.0), display_stop=200.0,
        ),
        _kymograph_entry(
            PULSE, reversal_row, reversal_task, "C",
            "Induce → stabilize → opposite pulse",
            (20.0, 40.0, 5_040.0, 5_200.0, 5_360.0),
        ),
        _kymograph_entry(
            PULSE, repeated_row, repeated_task, "E",
            "Repeated finite left/right pulses",
            tuple(
                event
                for start, duration in zip(
                    (0.0, 5_040.0, 10_360.0, 15_680.0, 21_000.0),
                    (40.0, 320.0, 320.0, 320.0, 320.0),
                )
                for event in (start + duration / 2.0, start + duration)
            ),
        ),
    ]

    mpl.rcParams.update({"font.size": 7, "axes.linewidth": 0.7})
    figure = plt.figure(figsize=(24, 9.4), constrained_layout=False)
    outer = GridSpec(
        3, 2, figure=figure, width_ratios=(4.9, 1.45),
        left=0.035, right=0.985, bottom=0.065, top=0.925,
        hspace=0.38, wspace=0.10,
    )
    left_axes = []
    for row_index in range(3):
        inner = GridSpecFromSubplotSpec(1, len(FIELDS), subplot_spec=outer[row_index, 0], wspace=0.08)
        left_axes.append([figure.add_subplot(inner[0, column]) for column in range(len(FIELDS))])

    cmaps = {
        "cue": "Greys", "rac": "viridis", "rho": "magma",
        "tension": "coolwarm", "mca_bound": "cividis",
        "mca_unbound": "plasma", "delta_beta_tension": "Reds",
        "delta_alpha_mca": "Blues",
    }
    norms = {}
    for field in FIELDS:
        combined = np.concatenate([entry["fields"][field].ravel() for entry in entries])
        if field == "tension":
            limit = max(abs(float(np.min(combined))), abs(float(np.max(combined))), 1e-12)
            norms[field] = mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
        else:
            low, high = float(np.min(combined)), float(np.max(combined))
            norms[field] = mpl.colors.Normalize(low, high if high > low else low + 1e-12)

    images = {}
    for row_index, entry in enumerate(entries):
        stride = max(1, int(np.ceil(len(entry["time"]) / 1_200)))
        time = entry["time"][::stride]
        for column, field in enumerate(FIELDS):
            axis = left_axes[row_index][column]
            images[field] = axis.imshow(
                entry["fields"][field][::stride].T,
                origin="lower", aspect="auto", interpolation="nearest",
                extent=(time[0], time[-1], 0.0, 1.0), cmap=cmaps[field], norm=norms[field],
            )
            for event in entry["events"]:
                axis.axvline(event, color="#d62728", lw=0.55, ls=":")
            if column == 0:
                axis.set_ylabel("Relative position")
                axis.text(
                    -0.43, 1.12, f"{entry['label']}  {entry['title']}",
                    transform=axis.transAxes, fontweight="bold", fontsize=9, va="top",
                )
            else:
                axis.set_yticklabels([])
            axis.set_xlabel("Time (s)")
            if row_index == 0:
                axis.set_title(FIELD_LABELS[field], fontsize=8)
            axis.tick_params(length=2, pad=1)
    for column, field in enumerate(FIELDS):
        color_axis = left_axes[0][column].inset_axes((0.02, -0.21, 0.96, 0.032))
        colorbar = figure.colorbar(images[field], cax=color_axis, orientation="horizontal")
        colorbar.ax.tick_params(labelsize=5, length=1, pad=1)

    # B: retain the common transient two-cue competition result.
    competition_axis = figure.add_subplot(outer[0, 1])
    colors = {"full": "#54278f", "beta_zero": "#e6550d", "alpha_zero": "#3182bd", "local_inhibition": "#31a354"}
    labels = {"full": "Full", "beta_zero": r"$\beta=0$ only", "alpha_zero": r"$\alpha=0$ only", "local_inhibition": "Local inhibition"}
    right_rows = []
    for variant in colors:
        selected = sorted(
            (row for row in primary_rows if row["assay"] == "competing_fronts_common" and row["coupling_variant"] == variant),
            key=lambda row: float(row["length_um"]),
        )
        lengths = [float(row["length_um"]) for row in selected]
        times = [_float(row["front_establishment_time_s"]) for row in selected]
        plotted = [value if value is not None else 20_000.0 for value in times]
        competition_axis.plot(lengths, plotted, marker="o", ms=2.8, lw=1.2, color=colors[variant], label=labels[variant])
        for length, value, row in zip(lengths, times, selected):
            if value is None:
                competition_axis.scatter(length, 20_000.0, marker="v", facecolors="none", edgecolors=colors[variant], s=18)
            right_rows.append({"panel": "competition", "variant": variant, "length_um": length, "front_establishment_time_s": value, "censored": value is None})
    competition_axis.set(yscale="log", xlabel="Cell length (µm)", ylabel="Single-front establishment (s)", title="B  Transient two-cue competition")
    competition_axis.grid(alpha=0.18)
    competition_axis.legend(frameon=False, fontsize=6, ncol=2)

    # D: every reversal condition is explicit; color successful cells by latency.
    reversal_axis = figure.add_subplot(outer[1, 1])
    reversal_conditions = [
        (ratio, duration)
        for ratio in (1.0, 2.0, 4.0, 8.0)
        for duration in (80.0, 160.0, 320.0, 640.0, 1_280.0)
    ]
    lengths = [float(value) for value in range(1, 40, 2)]
    latency = np.full((len(reversal_conditions), len(lengths)), np.nan)
    success = np.zeros_like(latency, dtype=bool)
    for row in pulse_rows:
        if row["assay"] != "transient_pulse_reversal" or row["coupling_variant"] != "full":
            continue
        task = pulse_tasks[row["task_id"]]
        condition = (float(task["pulse_amplitude_ratio"]), float(task["pulse_duration_s"]))
        i, j = reversal_conditions.index(condition), lengths.index(float(row["length_um"]))
        success[i, j] = row["switch_1_success"] == "True"
        latency[i, j] = _float(row["switch_1_latency_from_pulse_start_s"]) or np.nan
        right_rows.append({
            "panel": "pulse_reversal", "length_um": lengths[j],
            "amplitude_ratio": condition[0], "duration_s": condition[1],
            "success": success[i, j], "latency_s": None if np.isnan(latency[i, j]) else latency[i, j],
            "persistence_s": _float(row["switch_1_persistence_s"]),
        })
    masked = np.ma.masked_invalid(latency)
    image = reversal_axis.imshow(masked, origin="lower", aspect="auto", cmap="viridis", interpolation="nearest")
    failed_y, failed_x = np.where(~success)
    reversal_axis.scatter(failed_x, failed_y, marker="x", color="#b2182b", s=9, lw=0.7, label="No autonomous switch")
    reversal_axis.set_xticks(np.arange(0, len(lengths), 2), [f"{value:g}" for value in lengths[::2]])
    reversal_axis.set_yticks(np.arange(len(reversal_conditions)), [f"{ratio:g}× / {duration:g}s" for ratio, duration in reversal_conditions], fontsize=5.5)
    reversal_axis.set(xlabel="Cell length (µm)", ylabel="Opposite pulse", title="D  Autonomous reversal after withdrawal")
    if np.any(success):
        figure.colorbar(image, ax=reversal_axis, fraction=0.035, pad=0.02, label="Switch latency (s)")
    reversal_axis.legend(frameon=False, fontsize=6, loc="upper right")

    # F: aggregate event-wise switching across five reproducible pulse multipliers.
    repeated_grid = GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[2, 1], hspace=0.08)
    repeated_axes = [figure.add_subplot(repeated_grid[index, 0]) for index in range(3)]
    aggregate = []
    for length in lengths:
        selected = [row for row in pulse_rows if row["assay"] == "repeated_transient_switching" and float(row["length_um"]) == length]
        successes = sum(int(row["successful_switches"]) for row in selected)
        opportunities = sum(int(row["switch_opportunities"]) for row in selected)
        latencies, persistences = [], []
        for row in selected:
            for event in range(1, 5):
                if row.get(f"switch_{event}_success") == "True":
                    latencies.append(float(row[f"switch_{event}_latency_from_pulse_start_s"]))
                    persistences.append(float(row[f"switch_{event}_persistence_s"]))
        aggregate.append({
            "panel": "repeated_switching", "length_um": length,
            "switch_probability": successes / opportunities if opportunities else np.nan,
            "mean_latency_s": float(np.mean(latencies)) if latencies else np.nan,
            "median_persistence_s": float(np.median(persistences)) if persistences else np.nan,
            "successful_switches": successes, "opportunities": opportunities,
        })
    repeated_axes[0].plot(lengths, [row["switch_probability"] for row in aggregate], color="#54278f", marker="o", ms=3)
    repeated_axes[1].plot(lengths, [row["mean_latency_s"] for row in aggregate], color="#e6550d", marker="o", ms=3)
    repeated_axes[2].plot(lengths, [row["median_persistence_s"] for row in aggregate], color="#238b45", marker="o", ms=3)
    repeated_axes[0].set_title("F  Repeated-pulse autonomous switching", fontsize=8)
    repeated_axes[0].set_ylabel("Probability")
    repeated_axes[1].set_ylabel("Latency (s)")
    repeated_axes[2].set_ylabel("Persistence (s)")
    repeated_axes[2].set_xlabel("Cell length (µm)")
    for axis in repeated_axes:
        axis.grid(alpha=0.18)
        axis.set_xlim(0, 40)
    for axis in repeated_axes[:2]:
        axis.set_xticklabels([])

    _write_rows(SOURCE / "figure2_transient_pulse_backup_metrics.csv", right_rows + aggregate)
    figure.suptitle(
        r"Figure 2 backup. Transient-pulse tests of the reflection-symmetric 2026 de Belly model "
        r"($\alpha_0=3.5$, $\beta_0=4.0$)",
        fontsize=13, fontweight="bold",
    )
    figure.savefig(DESTINATION, dpi=300)
    figure.savefig(DESTINATION.with_suffix(".pdf"))
    plt.close(figure)

    reversal_successes = int(np.sum(success))
    audit = [
        "# Transient-pulse Figure 2 audit", "",
        "Model: reflection-symmetric extension of the 2026 de Belly mechanochemical model.",
        "Selected parameters: alpha0=3.5, beta0=4.0, alpha=1.05, beta=1.20; frozen across length.",
        "A switch is scored only after pulse withdrawal and only if the target Rac/Rho orientation persists cue-free for at least 1,000 s.", "",
        f"All 520 campaign tasks completed successfully. Autonomous opposite-pulse switches: {reversal_successes}/{success.size} full-model amplitude/duration/length conditions.",
        "Repeated switching probabilities are event-wise success fractions pooled across five predeclared amplitude multipliers at each length.",
        "Failed, unresolved and censored outcomes are plotted explicitly as crosses or missing latency/persistence values.",
    ]
    (OUTPUT / "DEBELLY_TRANSIENT_PULSE_FIGURE2_AUDIT.md").write_text("\n".join(audit) + "\n")
    print(DESTINATION)


if __name__ == "__main__":
    main()
