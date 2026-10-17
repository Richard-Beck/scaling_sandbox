#!/usr/bin/env python3
"""Plot the beta=1x opposing-cue duration trial as a standalone subpanel."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CONDITIONING = ROOT / "results/campaigns/debelly_2026_opposing_cue_duration_conditioning_trial"
FIELDS = (("rac", "Rac"), ("mca_bound", "MCA bound"))


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _half_mass(payload: dict[str, np.ndarray], field: str) -> np.ndarray:
    values = np.asarray(payload[f"field__{field}"], dtype=float)
    cumulative = np.cumsum(values, axis=1)
    half = 0.5 * cumulative[:, -1]
    crossing = np.argmax(cumulative >= half[:, None], axis=1)
    indices = np.arange(len(values))
    before = np.where(
        crossing > 0,
        cumulative[indices, np.maximum(crossing - 1, 0)],
        0.0,
    )
    fraction = np.divide(
        half - before,
        values[indices, crossing],
        out=np.full_like(half, 0.5),
        where=values[indices, crossing] > 0.0,
    )
    return np.clip((crossing + np.clip(fraction, 0.0, 1.0)) / values.shape[1], 0.0, 1.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--double-strength", action="store_true")
    arguments = parser.parse_args()
    if arguments.double_strength:
        name = "debelly_2026_opposing_cue_duration_double_strength_trial"
        durations = (60.0, 300.0, 1_200.0)
        output = ROOT / "results/manuscript_wp1_v3/standalone_alternatives/trials/opposing_cue_duration_beta1_double_strength"
        title = r"Opposing-cue duration trial, full model ($\beta=1\times$, opposing cue $4\times\alpha_0$)"
        stem = "Figure2_trial_opposing_cue_duration_beta1_double_strength"
        strength_text = "4*alpha0"
    else:
        name = "debelly_2026_opposing_cue_duration_response_trial"
        durations = (60.0, 300.0, 1_200.0, 7_200.0, 86_400.0)
        output = ROOT / "results/manuscript_wp1_v3/standalone_alternatives/trials/opposing_cue_duration_beta1"
        title = r"Opposing-cue duration trial, full model ($\beta=1\times$)"
        stem = "Figure2_trial_opposing_cue_duration_beta1"
        strength_text = "2*alpha0"
    campaign = ROOT / "results/campaigns" / name
    mpl.rcParams.update({
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.bbox": "tight",
    })
    output.mkdir(parents=True, exist_ok=True)
    rows = _read(campaign / "task_summary.csv")
    manifest = json.loads((ROOT / "campaigns" / f"{name}.json").read_text())
    tasks = {task["task_id"]: task for task in manifest["tasks"]}
    for row in rows:
        task = tasks[row["task_id"]]
        row["opposing_cue_duration_s"] = task["opposing_cue_duration_s"]
        row["opposing_cue_start_s"] = task["opposing_cue_start_s"]
    audit = _read(CONDITIONING / "conditioning_steady_state_audit.csv")
    accepted = sorted(float(row["length_um"]) for row in audit if row["both_reached_before_10000_s"].lower() == "true")
    discarded = sorted(float(row["length_um"]) for row in audit if row["both_reached_before_10000_s"].lower() != "true")

    figure, axes = plt.subplots(
        len(durations), len(FIELDS), figsize=(9.0, 3.0 + 2.2 * len(durations)),
        sharex=True, sharey=True, constrained_layout=True,
    )
    norm = mpl.colors.Normalize(1.0, 39.0)
    cmap = mpl.colormaps["viridis"]
    export_fields = (
        "task_id", "length_um", "opposing_cue_duration_s", "field",
        "opposing_cue_start_s", "time_from_opposing_cue_s", "x50_over_L",
    )
    with gzip.open(output / "opposing_cue_duration_polarity_trajectories.csv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=export_fields)
        writer.writeheader()
        for row_index, duration in enumerate(durations):
            selected_rows = sorted(
                (row for row in rows if np.isclose(float(row["opposing_cue_duration_s"]), duration)),
                key=lambda row: float(row["length_um"]),
            )
            if len(selected_rows) != len(accepted):
                raise RuntimeError(f"duration {duration:g}: expected {len(accepted)} rows, found {len(selected_rows)}")
            for row in selected_rows:
                length = float(row["length_um"])
                start = float(row["opposing_cue_start_s"])
                path = campaign / "tasks" / row["task_id"] / "trajectory.npz"
                with np.load(path, allow_pickle=False) as raw:
                    payload = {key: raw[key] for key in raw.files}
                absolute_time = np.asarray(payload["time"], dtype=float)
                relative_time = absolute_time - start
                display = relative_time >= -1_000.0
                stride = max(1, int(np.ceil(np.count_nonzero(display) / 1_500)))
                plot_indices = np.flatnonzero(display)[::stride]
                for column, (field, label) in enumerate(FIELDS):
                    score = _half_mass(payload, field)
                    axis = axes[row_index, column]
                    axis.plot(
                        relative_time[plot_indices], score[plot_indices],
                        color=cmap(norm(length)), lw=0.9, alpha=0.9,
                    )
                    writer.writerows({
                        "task_id": row["task_id"],
                        "length_um": length,
                        "opposing_cue_duration_s": duration,
                        "field": field,
                        "opposing_cue_start_s": start,
                        "time_from_opposing_cue_s": float(t),
                        "x50_over_L": float(value),
                    } for t, value in zip(relative_time, score))
                    axis.axvspan(0.0, duration, color="#f2b134", alpha=0.16, lw=0)
                    axis.axvline(0.0, color="#9a6b00", ls="--", lw=0.7)
                    axis.axvline(duration, color="#9a6b00", ls="--", lw=0.7)
                    axis.axhline(0.5, color="#777777", ls=":", lw=0.8)
                    axis.set_xscale("symlog", linthresh=300.0)
                    axis.set_xlim(-1_000.0, duration + 10_000.0)
                    axis.set_ylim(0.0, 1.0)
                    axis.grid(alpha=0.16)
                    if row_index == 0:
                        axis.set_title(label)
                    if row_index == len(durations) - 1:
                        axis.set_xlabel("Time from opposing-cue onset (s)")
                    if column == 0:
                        axis.set_ylabel(f"{duration:g}-s opposing cue\n$x_{{50}}/L$")

    colorbar = figure.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
        label="Cell length (µm)", shrink=0.78, pad=0.025,
    )
    colorbar.ax.tick_params(labelsize=7)
    figure.suptitle(title, fontsize=12, fontweight="bold")
    figure.savefig(output / f"{stem}.png", dpi=300)
    figure.savefig(output / f"{stem}.pdf")
    plt.close(figure)

    with (output / "conditioning_steady_state_audit.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit[0]))
        writer.writeheader()
        writer.writerows(audit)
    duration_text = ", ".join(f"{duration:g}" for duration in durations)
    (output / "README.md").write_text(
        "# Opposing-cue duration trial\n\n"
        "The intact beta=1x model receives the panel-E conditioning cue: a 600-s "
        "rectangular left outer-20% cue at 2*alpha0. Each length then waits until both "
        "Rac and bound-MCA x50/L have entered and remain within 0.002L of their final "
        "1,000-s medians. Lengths not satisfying both endpoints within 10,000 s after "
        "conditioning withdrawal are discarded.\n\n"
        "At the later of the two settling times, a right outer-20% cue at "
        f"{strength_text} is applied for {duration_text} s. Curves are aligned to that "
        "opposing-cue onset and colored by cell length.\n\n"
        f"Accepted lengths: {accepted}.\n\nDiscarded lengths: {discarded}.\n"
    )
    print(output)


if __name__ == "__main__":
    main()
