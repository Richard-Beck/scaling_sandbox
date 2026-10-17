#!/usr/bin/env python3
"""Build the Figure 2 full-model kymograph supplement."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000"
CAMPAIGN = ROOT / "results/campaigns" / CAMPAIGN_NAME
OUTPUT = (
    ROOT
    / "results/manuscript_wp1_v3/standalone_alternatives/"
    "v2_julia_activation_defaults_linear50/supplement"
)
LENGTHS = (15.0, 21.0, 31.0)
FIELDS = (
    ("rac", "Rac", "viridis"),
    ("rho", "Rho", "magma"),
    ("tension", "Tension", "coolwarm"),
    ("mca_bound", "MCA bound", "cividis"),
)
ORIENTATION_FIELDS = {"rac", "mca_bound"}
STIMULUS_START = 200.0
STIMULUS_STOP = 360.0
MAXIMUM_TIME_SLICES = 1_500


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _trajectory(row: dict[str, str]) -> dict[str, np.ndarray]:
    path = CAMPAIGN / "tasks" / row["task_id"] / "trajectory.npz"
    with np.load(path, allow_pickle=False) as raw:
        return {key: raw[key] for key in raw.files}


def _half_mass_position(payload: dict[str, np.ndarray], field: str) -> np.ndarray:
    x = np.asarray(payload["x"], dtype=float)
    values = np.asarray(payload[f"field__{field}"], dtype=float)
    dx = float(np.median(np.diff(x)))
    left_edge = float(x[0] - dx / 2.0)
    length = float(x[-1] - x[0] + dx)
    cumulative = np.cumsum(values, axis=1)
    half = 0.5 * cumulative[:, -1]
    crossing = np.argmax(cumulative >= half[:, None], axis=1)
    rows = np.arange(len(values))
    before = np.where(
        crossing > 0,
        cumulative[rows, np.maximum(crossing - 1, 0)],
        0.0,
    )
    fraction = np.divide(
        half - before,
        values[rows, crossing],
        out=np.full_like(half, 0.5),
        where=values[rows, crossing] > 0.0,
    )
    position = left_edge + (crossing + np.clip(fraction, 0.0, 1.0)) * dx
    return np.clip((position - left_edge) / length, 0.0, 1.0)


def _edges(centers: np.ndarray) -> np.ndarray:
    centers = np.asarray(centers, dtype=float)
    midpoints = 0.5 * (centers[:-1] + centers[1:])
    first = centers[0] - 0.5 * (centers[1] - centers[0])
    last = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return np.concatenate(([first], midpoints, [last]))


def main() -> None:
    mpl.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    rows = _rows(CAMPAIGN / "task_summary.csv")
    selected = {}
    for length in LENGTHS:
        matches = [
            row for row in rows
            if row["assay"] == "single_stimulus_length_test"
            and row["coupling_variant"] == "full"
            and np.isclose(float(row["length_um"]), length)
        ]
        if len(matches) != 1:
            raise RuntimeError(f"expected one full-model L={length:g} trajectory, found {len(matches)}")
        selected[length] = _trajectory(matches[0])

    norms = {}
    for field, _, _ in FIELDS:
        arrays = [np.asarray(selected[length][f"field__{field}"], dtype=float) for length in LENGTHS]
        low = min(float(np.min(values)) for values in arrays)
        high = max(float(np.max(values)) for values in arrays)
        if field == "tension":
            limit = max(abs(low), abs(high), 1e-12)
            norms[field] = mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
        else:
            norms[field] = mpl.colors.Normalize(vmin=low, vmax=high)

    figure, axes = plt.subplots(
        len(LENGTHS), len(FIELDS), figsize=(16.0, 10.5),
        sharex=True, sharey=True, constrained_layout=True,
    )
    mappables = {}
    orientation_rows = []
    exported = {}
    for row_index, length in enumerate(LENGTHS):
        payload = selected[length]
        absolute_time = np.asarray(payload["time"], dtype=float)
        relative_time = absolute_time - STIMULUS_START
        stride_step = max(1, int(np.ceil(len(relative_time) / MAXIMUM_TIME_SLICES)))
        indices = np.arange(0, len(relative_time), stride_step)
        if indices[-1] != len(relative_time) - 1:
            indices = np.append(indices, len(relative_time) - 1)
        displayed_time = relative_time[indices]
        position = np.asarray(payload["x"], dtype=float) / length
        exported[f"L{length:g}_time_from_stimulus_s"] = relative_time
        exported[f"L{length:g}_position_over_L"] = position

        for column, (field, label, cmap) in enumerate(FIELDS):
            axis = axes[row_index, column]
            values = np.asarray(payload[f"field__{field}"], dtype=float)
            mesh = axis.pcolormesh(
                _edges(displayed_time), _edges(position), values[indices].T,
                shading="flat", cmap=cmap, norm=norms[field], rasterized=True,
            )
            mappables[field] = mesh
            axis.set_xscale("symlog", linthresh=200.0, linscale=1.0)
            axis.set_xlim(relative_time[0], relative_time[-1])
            axis.set_ylim(0.0, 1.0)
            axis.axvspan(0.0, STIMULUS_STOP - STIMULUS_START, color="#f2b134", alpha=0.16, lw=0)
            axis.axvline(0.0, color="white", ls="--", lw=0.7, alpha=0.9)
            axis.axvline(STIMULUS_STOP - STIMULUS_START, color="white", ls="--", lw=0.7, alpha=0.9)
            if field in ORIENTATION_FIELDS:
                orientation = _half_mass_position(payload, field)
                axis.plot(relative_time, orientation, color="white", lw=2.2, alpha=0.9)
                axis.plot(relative_time, orientation, color="black", lw=0.75, alpha=0.95)
                orientation_rows.extend({
                    "length_um": length,
                    "field": field,
                    "time_s": float(time),
                    "time_from_stimulus_s": float(time - STIMULUS_START),
                    "x50_over_L": float(value),
                } for time, value in zip(absolute_time, orientation))
                exported[f"L{length:g}_{field}_x50_over_L"] = orientation
            exported[f"L{length:g}_{field}"] = values
            if row_index == 0:
                axis.set_title(label)
            if column == 0:
                axis.set_ylabel(f"L={length:g} µm\nPosition, x/L")
            if row_index == len(LENGTHS) - 1:
                axis.set_xlabel("Time from stimulus onset (s; symlog)")

    for column, (field, label, _) in enumerate(FIELDS):
        colorbar = figure.colorbar(
            mappables[field], ax=axes[:, column].tolist(),
            orientation="horizontal", pad=0.055, fraction=0.04,
        )
        colorbar.set_label(label)
        colorbar.ax.tick_params(labelsize=7)

    figure.suptitle(
        "Full-model spatiotemporal responses across cell length",
        fontsize=13, fontweight="bold",
    )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT / "Figure2_supplement_full_model_kymographs.png", dpi=300)
    figure.savefig(OUTPUT / "Figure2_supplement_full_model_kymographs.pdf")
    plt.close(figure)

    with gzip.open(OUTPUT / "figure2_supplement_orientation_trajectories.csv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=(
            "length_um", "field", "time_s", "time_from_stimulus_s", "x50_over_L",
        ))
        writer.writeheader()
        writer.writerows(orientation_rows)
    np.savez_compressed(OUTPUT / "figure2_supplement_kymograph_data.npz", **exported)
    (OUTPUT / "README.md").write_text(
        "# Figure 2 kymograph supplement\n\n"
        "Rows show full-model cells of length 15, 21, and 31 µm; columns show Rac, Rho, "
        "tension, and bound MCA. Time is plotted relative to the 160-s stimulus onset on a "
        "symmetric-log axis with a 200-s linear region. Color normalization is shared across "
        "lengths within each field. Rac and bound-MCA panels overlay the interpolated half-mass "
        "orientation trajectory x50/L as a black line with a white outline.\n"
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
