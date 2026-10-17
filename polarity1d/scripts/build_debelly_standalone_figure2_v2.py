#!/usr/bin/env python3
"""Build standalone Figure 2 v2 from the completed alternatives and reversal runs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_NAME = "debelly_2026_continuous_polarity_alternatives_v1"
A_NAME = BASE_NAME
TRAJECTORY_NAME = "debelly_2026_opposing_cue_duration_response_trial"
SWEEP_NAME = "debelly_2026_countercue_threshold_sweep"
BASE = ROOT / "results/campaigns" / BASE_NAME
A_CAMPAIGN = BASE
TRAJECTORY_CAMPAIGN = ROOT / "results/campaigns" / TRAJECTORY_NAME
SWEEP_CAMPAIGN = ROOT / "results/campaigns" / SWEEP_NAME
OUTPUT = ROOT / "results/manuscript_wp1_v3/standalone_alternatives/v2"
OUTPUT_STEM = "Figure2_standalone_v2"
ACTIVATION_DESCRIPTION = "beta variants at 0x, 1x, and 2x the default beta=0.54"
FIGURE_NOTE = ""

MODEL_ARMS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
MODEL_LABELS = {
    "full": "Full model",
    "beta_zero": r"$\beta=0$ only",
    "alpha_zero": r"$\alpha=0$ only",
    "local_inhibition": "Local inhibition",
}
BETA_ARMS = ("beta_0x", "beta_1x", "beta_2x")
BETA_LABELS = {
    "beta_0x": r"$\beta=0\times$",
    "beta_1x": r"$\beta=1\times$",
    "beta_2x": r"$\beta=2\times$",
}
BETA_COLORS = {"beta_0x": "#d95f0e", "beta_1x": "#54278f", "beta_2x": "#238b45"}
FIELDS = (("rac", "Rac"), ("mca_bound", "MCA bound"))
SETTLING_TOLERANCE = 0.002
SETTLING_CONFIRMATION_S = 1_000.0
COMPETITION_CUE_STOP = 80.0
TRAJECTORY_DURATIONS = (300.0, 1_200.0, 7_200.0)
MINIMUM_POLARITY_DISPLACEMENT = 0.02


def _scale_polarity_axes(axes: np.ndarray) -> None:
    """Use compact, shared y limits for each field column in a panel."""
    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes[None, :]
    for column in range(axes.shape[1]):
        values = []
        for axis in axes[:, column]:
            for line in axis.lines:
                data = np.asarray(line.get_ydata(), dtype=float)
                if data.size > 2:
                    values.extend(data[np.isfinite(data)])
        if not values:
            continue
        low, high = float(np.min(values)), float(np.max(values))
        padding = max(0.02, 0.08 * (high - low))
        low = max(0.0, np.floor((low - padding) / 0.05) * 0.05)
        high = min(1.0, np.ceil((high + padding) / 0.05) * 0.05)
        if high - low < 0.1:
            midpoint = 0.5 * (low + high)
            low, high = max(0.0, midpoint - 0.05), min(1.0, midpoint + 0.05)
        for axis in axes[:, column]:
            axis.set_ylim(low, high)


def _style() -> None:
    mpl.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
    })


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], key: str) -> float:
    try:
        value = float(row.get(key, "nan"))
    except (TypeError, ValueError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def _trajectory(campaign: Path, row: dict[str, str]) -> dict[str, np.ndarray]:
    path = campaign / "tasks" / row["task_id"] / "trajectory.npz"
    with np.load(path, allow_pickle=False) as raw:
        return {key: raw[key] for key in raw.files}


def _half_mass_position(payload: dict[str, np.ndarray], field: str) -> np.ndarray:
    """Interpolate the normalized position at which cumulative field mass is 50%."""
    x = np.asarray(payload["x"], dtype=float)
    values = np.asarray(payload[f"field__{field}"], dtype=float)
    dx = float(np.median(np.diff(x)))
    left_edge = float(x[0] - dx / 2.0)
    length = float(x[-1] - x[0] + dx)
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
    position = left_edge + (crossing + np.clip(fraction, 0.0, 1.0)) * dx
    return np.clip((position - left_edge) / length, 0.0, 1.0)


def _settling_time(
    time: np.ndarray,
    score: np.ndarray,
    *,
    after: float,
) -> tuple[float | None, float]:
    """Return first sustained entry into a final-tail tolerance band."""
    tail = time >= time[-1] - SETTLING_CONFIRMATION_S
    reference = float(np.median(score[tail]))
    if float(np.max(np.abs(score[tail] - reference))) > SETTLING_TOLERANCE:
        return None, reference
    within = np.abs(score - reference) <= SETTLING_TOLERANCE
    remains_within = np.logical_and.accumulate(within[::-1])[::-1]
    eligible = (
        remains_within
        & (time >= after)
        & (time <= time[-1] - SETTLING_CONFIRMATION_S)
    )
    indices = np.flatnonzero(eligible)
    return (None if len(indices) == 0 else float(time[indices[0]])), reference


def _center_crossing_time(
    time: np.ndarray,
    score: np.ndarray,
    *,
    after: float,
    final_reference: float,
) -> float | None:
    """Find the first interpolated 0.5L crossing toward the final polarity side."""
    direction = np.sign(final_reference - 0.5)
    if direction == 0:
        return None
    indices = np.flatnonzero(time >= after)
    if len(indices) < 2:
        return None
    directed = direction * (score - 0.5)
    candidates = indices[1:][
        (directed[indices[1:]] >= 0.0) & (directed[indices[:-1]] < 0.0)
    ]
    if len(candidates) == 0:
        return None
    index = int(candidates[0])
    y0, y1 = float(score[index - 1]), float(score[index])
    if np.isclose(y0, y1):
        return float(time[index])
    weight = (0.5 - y0) / (y1 - y0)
    return float(time[index - 1] + weight * (time[index] - time[index - 1]))


def _plot_stride(time: np.ndarray, maximum: int = 1_200) -> slice:
    return slice(None, None, max(1, int(np.ceil(len(time) / maximum))))


def _format_time_axis(axis: plt.Axes, *, reversal: bool = False) -> None:
    axis.set_xscale("symlog", linthresh=200.0 if not reversal else 600.0)
    axis.grid(alpha=0.16)
    axis.axhline(0.5, color="#777777", ls=":", lw=0.7)


def _plot_length_trajectories(
    axes: np.ndarray,
    rows: list[dict[str, str]],
    campaign: Path,
    *,
    assay: str,
    facets: tuple[str, ...],
    labels: dict[str, str],
    time_offset: float,
    event_spans: tuple[tuple[float, float, str], ...],
    mark_settling: bool,
    mark_only_polarized: bool = False,
    settling_after: float,
    trajectory_writer: csv.DictWriter | None,
    panel_name: str,
) -> list[dict]:
    norm = mpl.colors.Normalize(1.0, 39.0)
    cmap = mpl.colormaps["viridis"]
    summaries = []
    for row_index, facet in enumerate(facets):
        selected = sorted(
            (row for row in rows if row["assay"] == assay and row["coupling_variant"] == facet),
            key=lambda row: _number(row, "length_um"),
        )
        if len(selected) != 20:
            raise RuntimeError(f"expected 20 rows for {assay}/{facet}, found {len(selected)}")
        for row in selected:
            length = _number(row, "length_um")
            payload = _trajectory(campaign, row)
            absolute_time = np.asarray(payload["time"], dtype=float)
            display_time = absolute_time - time_offset
            stride = _plot_stride(display_time)
            for column, (field, field_label) in enumerate(FIELDS):
                axis = axes[row_index, column]
                score = _half_mass_position(payload, field)
                settling, reference = _settling_time(absolute_time, score, after=settling_after)
                axis.plot(
                    display_time[stride], score[stride],
                    color=cmap(norm(length)), lw=0.75, alpha=0.86,
                )
                polarized = abs(reference - 0.5) >= MINIMUM_POLARITY_DISPLACEMENT
                if (
                    mark_settling
                    and settling is not None
                    and (polarized or not mark_only_polarized)
                ):
                    at = int(np.searchsorted(absolute_time, settling))
                    axis.scatter(
                        settling - time_offset, score[at], s=9,
                        marker="o" if polarized else "x",
                        facecolor=cmap(norm(length)) if polarized else "none",
                        edgecolor="black" if polarized else cmap(norm(length)),
                        linewidth=0.25 if polarized else 0.7,
                        zorder=4,
                    )
                summaries.append({
                    "panel": panel_name,
                    "task_id": row["task_id"],
                    "variant": facet,
                    "beta_multiplier": row.get("beta_multiplier", ""),
                    "field": field,
                    "length_um": length,
                    "steady_state_reached": settling is not None,
                    "steady_state_time_s": "" if settling is None else settling,
                    "polarized_steady_state": (
                        settling is not None
                        and abs(reference - 0.5) >= MINIMUM_POLARITY_DISPLACEMENT
                    ),
                    "steady_state_time_from_analysis_start_s": (
                        "" if settling is None else settling - settling_after
                    ),
                    "final_tail_median_x50_over_L": reference,
                    "settling_tolerance_fraction_L": SETTLING_TOLERANCE,
                    "settling_confirmation_s": SETTLING_CONFIRMATION_S,
                })
                if trajectory_writer is not None:
                    trajectory_writer.writerows({
                        "panel": panel_name,
                        "task_id": row["task_id"],
                        "variant": facet,
                        "field": field,
                        "length_um": length,
                        "time_s": float(t),
                        "display_time_s": float(dt),
                        "x50_over_L": float(value),
                    } for t, dt, value in zip(absolute_time, display_time, score))
        for column, (_, field_label) in enumerate(FIELDS):
            axis = axes[row_index, column]
            _format_time_axis(axis, reversal=panel_name == "E")
            for start, stop, color in event_spans:
                axis.axvspan(start - time_offset, stop - time_offset, color=color, alpha=0.12, lw=0)
            axis.set_ylim(0.0, 1.0)
            if row_index == 0:
                axis.set_title(field_label)
            if row_index == len(facets) - 1:
                axis.set_xlabel("Time from stimulus onset (s)" if panel_name != "E" else "Time from reversal onset (s)")
            if column == 0:
                axis.set_ylabel(f"{labels[facet]}\n$x_{50}/L$")
    return summaries


def _panel_a(rows, axes, writer) -> list[dict]:
    summaries = _plot_length_trajectories(
        axes, rows, A_CAMPAIGN,
        assay="single_stimulus_length_test",
        facets=MODEL_ARMS,
        labels=MODEL_LABELS,
        time_offset=200.0,
        event_spans=((200.0, 360.0, "#f1c232"),),
        mark_settling=True,
        mark_only_polarized=True,
        settling_after=360.0,
        trajectory_writer=writer,
        panel_name="A",
    )
    _scale_polarity_axes(axes)
    return summaries


def _panel_a_summary(axis: plt.Axes, steady_rows: list[dict]) -> None:
    colors = {"rac": "#1b9e77", "mca_bound": "#7570b3"}
    for field, label in FIELDS:
        selected = sorted(
            (
                row for row in steady_rows
                if row["variant"] == "full"
                and row["field"] == field
                and row["steady_state_reached"]
                and row["polarized_steady_state"]
            ),
            key=lambda row: row["length_um"],
        )
        axis.plot(
            [row["length_um"] for row in selected],
            [row["steady_state_time_s"] - 200.0 for row in selected],
            "o-", color=colors[field], ms=3.5, lw=1.2, label=label,
        )
    axis.set(
        xlabel="Cell length (µm)",
        ylabel="Polarization time from stimulus onset (s)",
        title="Full model",
    )
    axis.set_yscale("log")
    axis.grid(alpha=0.18)
    axis.legend(frameon=False)


def _panel_b(rows, axes) -> None:
    matches = [
        row for row in rows
        if row["assay"] == "single_stimulus_length_test"
        and row["coupling_variant"] == "full"
        and np.isclose(_number(row, "length_um"), 15.0)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one L=15 full-model trajectory, found {len(matches)}")
    payload = _trajectory(A_CAMPAIGN, matches[0])
    time = np.asarray(payload["time"], dtype=float) - 200.0
    x = np.asarray(payload["x"], dtype=float) / 15.0
    field_specs = (
        ("rac", "Rac", "viridis"),
        ("rho", "Rho", "magma"),
        ("tension", "Tension", "coolwarm"),
        ("mca_bound", "MCA bound", "cividis"),
    )
    exported = {"time_from_stimulus_s": time, "position_over_L": x}
    for axis, (field, label, cmap) in zip(axes, field_specs):
        values = np.asarray(payload[f"field__{field}"], dtype=float)
        stride = _plot_stride(time, maximum=1_500)
        norm = None
        if field == "tension":
            limit = max(abs(float(np.min(values))), abs(float(np.max(values))), 1e-12)
            norm = mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
        image = axis.imshow(
            values[stride].T,
            origin="lower", aspect="auto", interpolation="nearest",
            extent=(time[stride][0], time[stride][-1], x[0], x[-1]),
            cmap=cmap, norm=norm,
        )
        axis.axvspan(0.0, 160.0, color="white", alpha=0.08, lw=0)
        axis.axvline(160.0, color="white", ls="--", lw=0.7)
        axis.set_title(label)
        axis.set_xlabel("Time from stimulus onset (s)")
        if axis is axes[0]:
            axis.set_ylabel("Position, x/L")
        else:
            axis.set_yticklabels([])
        colorbar = axis.figure.colorbar(image, ax=axis, orientation="horizontal", pad=0.18, fraction=0.05)
        colorbar.ax.tick_params(labelsize=5, length=1)
        exported[field] = values
    np.savez_compressed(OUTPUT / "figure2_v2_panel_B_L15_kymographs.npz", **exported)


def _panel_c(rows, axes, writer) -> list[dict]:
    summaries = _plot_length_trajectories(
        axes, rows, BASE,
        assay="figure2_B_beta_sweep",
        facets=BETA_ARMS,
        labels=BETA_LABELS,
        time_offset=0.0,
        event_spans=((0.0, COMPETITION_CUE_STOP, "#f1c232"),),
        mark_settling=True,
        settling_after=COMPETITION_CUE_STOP,
        trajectory_writer=writer,
        panel_name="C",
    )
    _scale_polarity_axes(axes)
    axes[0, 1].legend(
        handles=(
            Line2D([], [], color="black", marker="o", ls="none", ms=4, label="Polarized steady state"),
            Line2D([], [], color="black", marker="x", ls="none", ms=4, label="Unpolarized steady state"),
        ),
        frameon=False,
        loc="upper left",
    )
    return summaries


def _panel_d(axis_pair, steady_rows: list[dict]) -> None:
    for axis, (field, label) in zip(axis_pair, FIELDS):
        for variant in BETA_ARMS:
            selected = sorted(
                (row for row in steady_rows if row["panel"] == "C" and row["field"] == field and row["variant"] == variant),
                key=lambda row: row["length_um"],
            )
            passed = [
                row for row in selected
                if row["steady_state_reached"] and row["polarized_steady_state"]
            ]
            axis.plot(
                [row["length_um"] for row in passed],
                [row["steady_state_time_s"] for row in passed],
                "o-", color=BETA_COLORS[variant], ms=3, lw=1.1, label=BETA_LABELS[variant],
            )
        axis.set(xlabel="Cell length (µm)", ylabel="Time to steady state (s)", title=label)
        axis.set_yscale("log")
        axis.grid(alpha=0.18)
    axis_pair[0].legend(frameon=False)


def _panel_e_duration_trajectories(rows, tasks, axes, writer) -> None:
    norm = mpl.colors.Normalize(1.0, 39.0)
    cmap = mpl.colormaps["viridis"]
    for row_index, duration in enumerate(TRAJECTORY_DURATIONS):
        selected = sorted(
            (
                row for row in rows
                if np.isclose(float(tasks[row["task_id"]]["opposing_cue_duration_s"]), duration)
            ),
            key=lambda row: _number(row, "length_um"),
        )
        for row in selected:
            task = tasks[row["task_id"]]
            start = float(task["opposing_cue_start_s"])
            length = _number(row, "length_um")
            payload = _trajectory(TRAJECTORY_CAMPAIGN, row)
            absolute_time = np.asarray(payload["time"], dtype=float)
            relative_time = absolute_time - start
            display = relative_time >= -1_000.0
            indices = np.flatnonzero(display)[_plot_stride(relative_time[display])]
            for column, (field, label) in enumerate(FIELDS):
                score = _half_mass_position(payload, field)
                axis = axes[row_index, column]
                axis.plot(relative_time[indices], score[indices], color=cmap(norm(length)), lw=0.75, alpha=0.88)
                writer.writerows({
                    "panel": "E", "task_id": row["task_id"], "variant": "beta_1x",
                    "field": field, "length_um": length, "time_s": float(t),
                    "display_time_s": float(t - start), "x50_over_L": float(value),
                } for t, value in zip(absolute_time, score))
                axis.axvspan(0.0, duration, color="#f2b134", alpha=0.16, lw=0)
                axis.axvline(0.0, color="#9a6b00", ls="--", lw=0.7)
                axis.axvline(duration, color="#9a6b00", ls="--", lw=0.7)
                axis.axhline(0.5, color="#777777", ls=":", lw=0.7)
                axis.set_xscale("symlog", linthresh=300.0)
                axis.set_xlim(-1_000.0, duration + 10_000.0)
                axis.set_ylim(0.0, 1.0)
                axis.grid(alpha=0.16)
                if row_index == 0:
                    axis.set_title(label)
                if row_index == len(TRAJECTORY_DURATIONS) - 1:
                    axis.set_xlabel("Time from counter-cue onset (s)")
                if column == 0:
                    axis.set_ylabel(f"{duration:g}-s cue\n$x_{{50}}/L$")
    _scale_polarity_axes(axes)


def _strict_reversal(row: dict[str, str], field: str) -> bool:
    return bool(
        row[f"countercue_{field}_reversal_success"].lower() == "true"
        and abs(_number(row, f"countercue_{field}_conditioning_reference_x50_over_L") - 0.5)
        >= MINIMUM_POLARITY_DISPLACEMENT
        and abs(_number(row, f"countercue_{field}_final_x50_over_L") - 0.5)
        >= MINIMUM_POLARITY_DISPLACEMENT
    )


def _minimum_duration_rows(rows: list[dict[str, str]]) -> list[dict]:
    result = []
    for variant in BETA_ARMS:
        lengths = sorted({_number(row, "length_um") for row in rows if row["coupling_variant"] == variant})
        for field in ("rac", "mca"):
            for length in lengths:
                candidates = [
                    round(_number(row, "countercue_duration_s"), 9) for row in rows
                    if row["coupling_variant"] == variant
                    and np.isclose(_number(row, "length_um"), length)
                    and _strict_reversal(row, field)
                ]
                if candidates:
                    result.append({
                        "variant": variant, "field": field, "length_um": length,
                        "minimum_successful_countercue_duration_s": min(candidates),
                    })
    return result


def _countercue_audit_rows(rows: list[dict[str, str]]) -> list[dict]:
    exported = []
    for row in rows:
        exported.append({
            "task_id": row["task_id"],
            "variant": row["coupling_variant"],
            "length_um": _number(row, "length_um"),
            "countercue_duration_s": round(_number(row, "countercue_duration_s"), 9),
            "rac_reversal_success": _strict_reversal(row, "rac"),
            "mca_reversal_success": _strict_reversal(row, "mca"),
            "rac_conditioning_x50_over_L": _number(row, "countercue_rac_conditioning_reference_x50_over_L"),
            "rac_final_x50_over_L": _number(row, "countercue_rac_final_x50_over_L"),
            "mca_conditioning_x50_over_L": _number(row, "countercue_mca_conditioning_reference_x50_over_L"),
            "mca_final_x50_over_L": _number(row, "countercue_mca_final_x50_over_L"),
            "rac_final_steady": row["countercue_rac_final_steady"],
            "mca_final_steady": row["countercue_mca_final_steady"],
            "target_tolerance_fraction_L": _number(row, "countercue_flip_tolerance_fraction_L"),
            "minimum_polarity_displacement_fraction_L": MINIMUM_POLARITY_DISPLACEMENT,
        })
    return exported


def _panel_f(axis_pair, rows: list[dict]) -> None:
    for axis, (field, label) in zip(axis_pair, (("rac", "Rac-defined reversal"), ("mca", "MCA-defined reversal"))):
        for variant in BETA_ARMS:
            selected = sorted(
                (row for row in rows if row["field"] == field and row["variant"] == variant),
                key=lambda row: row["length_um"],
            )
            lookup = {row["length_um"]: row["minimum_successful_countercue_duration_s"] for row in selected}
            lengths = np.asarray(tuple(float(value) for value in range(1, 40, 2)))
            durations = np.asarray([lookup.get(length, np.nan) for length in lengths])
            axis.plot(
                lengths, durations,
                "o-", color=BETA_COLORS[variant], ms=3, lw=1.1, label=BETA_LABELS[variant],
            )
        axis.set(xlabel="Cell length (µm)", ylabel="Minimum successful cue duration (s)", title=label)
        axis.grid(alpha=0.18)
    axis_pair[0].legend(frameon=False)


def _write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    global BASE_NAME, A_NAME, TRAJECTORY_NAME, SWEEP_NAME
    global BASE, A_CAMPAIGN, TRAJECTORY_CAMPAIGN, SWEEP_CAMPAIGN, OUTPUT
    global OUTPUT_STEM, ACTIVATION_DESCRIPTION, FIGURE_NOTE, TRAJECTORY_DURATIONS

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--julia-activation-defaults",
        action="store_true",
        help="Use the parallel campaigns with alpha0=beta0=1 and alpha=beta=0.2.",
    )
    parser.add_argument(
        "--julia-activation-defaults-short-cues",
        action="store_true",
        help="Use Julia activation defaults with 150/300/600-s trajectories and a 10-1000-s sweep.",
    )
    parser.add_argument(
        "--julia-activation-defaults-linear50",
        action="store_true",
        help="Use Julia activation defaults with a brute-force 50-duration linear sweep.",
    )
    arguments = parser.parse_args()
    if (
        arguments.julia_activation_defaults
        or arguments.julia_activation_defaults_short_cues
        or arguments.julia_activation_defaults_linear50
    ):
        BASE_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_base"
        A_NAME = "debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000"
        TRAJECTORY_NAME = (
            "debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response"
            if arguments.julia_activation_defaults_linear50
            else "debelly_2026_standalone_figure2_julia_activation_defaults_short_cues_response"
            if arguments.julia_activation_defaults_short_cues
            else "debelly_2026_standalone_figure2_julia_activation_defaults_response"
        )
        SWEEP_NAME = TRAJECTORY_NAME
        BASE = ROOT / "results/campaigns" / BASE_NAME
        A_CAMPAIGN = ROOT / "results/campaigns" / A_NAME
        TRAJECTORY_CAMPAIGN = ROOT / "results/campaigns" / TRAJECTORY_NAME
        SWEEP_CAMPAIGN = TRAJECTORY_CAMPAIGN
        OUTPUT = ROOT / "results/manuscript_wp1_v3/standalone_alternatives/v2_julia_activation_defaults"
        OUTPUT_STEM = "Figure2_standalone_v2_julia_activation_defaults"
        ACTIVATION_DESCRIPTION = (
            "beta variants at 0x, 1x, and 2x the Julia-default beta=0.2; "
            "alpha0=beta0=1 and alpha=0.2"
        )
        FIGURE_NOTE = r"Julia activation defaults: $\alpha_0=\beta_0=1$, $\alpha=0.2$, $\beta=0/0.2/0.4$"
        if arguments.julia_activation_defaults_short_cues or arguments.julia_activation_defaults_linear50:
            suffix = "linear50" if arguments.julia_activation_defaults_linear50 else "short_cues"
            OUTPUT = ROOT / f"results/manuscript_wp1_v3/standalone_alternatives/v2_julia_activation_defaults_{suffix}"
            OUTPUT_STEM = f"Figure2_standalone_v2_julia_activation_defaults_{suffix}"
            TRAJECTORY_DURATIONS = (150.0, 300.0, 600.0)

    _style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    base_rows = _read(BASE / "task_summary.csv")
    a_rows = _read(A_CAMPAIGN / "task_summary.csv")
    trajectory_rows = _read(TRAJECTORY_CAMPAIGN / "task_summary.csv")
    trajectory_manifest = json.loads((ROOT / "campaigns" / f"{TRAJECTORY_NAME}.json").read_text())
    trajectory_tasks = {task["task_id"]: task for task in trajectory_manifest["tasks"]}
    sweep_rows = _read(SWEEP_CAMPAIGN / "task_summary.csv")
    if (
        arguments.julia_activation_defaults
        or arguments.julia_activation_defaults_short_cues
        or arguments.julia_activation_defaults_linear50
    ):
        trajectory_rows = [row for row in trajectory_rows if row["assay"] == "countercue_duration_trajectories"]
        sweep_rows = [row for row in sweep_rows if row["assay"] == "countercue_threshold_sweep"]
    minimum_durations = _minimum_duration_rows(sweep_rows)

    figure = plt.figure(figsize=(30.0, 15.0))
    outer = figure.add_gridspec(
        2, 3,
        height_ratios=(3.2, 1.25),
        left=0.035, right=0.965, top=0.955, bottom=0.07,
        hspace=0.34, wspace=0.22,
    )
    a_grid = outer[0, 0].subgridspec(4, 2, hspace=0.30, wspace=0.20)
    c_grid = outer[0, 1].subgridspec(3, 2, hspace=0.30, wspace=0.20)
    e_grid = outer[0, 2].subgridspec(3, 2, hspace=0.30, wspace=0.20)
    d_grid = outer[1, 1].subgridspec(1, 2, wspace=0.28)
    e_summary_grid = outer[1, 2].subgridspec(1, 2, wspace=0.28)
    a_axes = np.asarray([[figure.add_subplot(a_grid[i, j]) for j in range(2)] for i in range(4)])
    b_axis = figure.add_subplot(outer[1, 0])
    c_axes = np.asarray([[figure.add_subplot(c_grid[i, j]) for j in range(2)] for i in range(3)])
    d_axes = np.asarray([figure.add_subplot(d_grid[0, j]) for j in range(2)])
    e_axes = np.asarray([[figure.add_subplot(e_grid[i, j]) for j in range(2)] for i in range(3)])
    e_summary_axes = np.asarray([figure.add_subplot(e_summary_grid[0, j]) for j in range(2)])

    trajectory_fields = ("panel", "task_id", "variant", "field", "length_um", "time_s", "display_time_s", "x50_over_L")
    with gzip.open(OUTPUT / "figure2_v2_polarity_trajectories.csv.gz", "wt", newline="") as handle:
        trajectory_writer = csv.DictWriter(handle, fieldnames=trajectory_fields)
        trajectory_writer.writeheader()
        steady_a = _panel_a(a_rows, a_axes, trajectory_writer)
        _panel_a_summary(b_axis, steady_a)
        steady_c = _panel_c(base_rows, c_axes, trajectory_writer)
        _panel_d(d_axes, steady_c)
        _panel_e_duration_trajectories(trajectory_rows, trajectory_tasks, e_axes, trajectory_writer)
        _panel_f(e_summary_axes, minimum_durations)

    for axis, label in (
        (a_axes[0, 0], "A"), (b_axis, "B"),
        (c_axes[0, 0], "C"), (d_axes[0], "D"),
        (e_axes[0, 0], "E"), (e_summary_axes[0], "F"),
    ):
        axis.text(-0.13, 1.10, label, transform=axis.transAxes, fontsize=14, fontweight="bold", va="top")
    a_axes[0, 0].text(0.0, 1.24, "Single 160-s stimulus", transform=a_axes[0, 0].transAxes, fontsize=11, fontweight="bold")
    b_axis.text(0.0, 1.24, "Full-model polarization time", transform=b_axis.transAxes, fontsize=11, fontweight="bold")
    c_axes[0, 0].text(0.0, 1.24, "Full-model beta variants: polarity and first steady state", transform=c_axes[0, 0].transAxes, fontsize=11, fontweight="bold")
    d_axes[0].text(0.0, 1.24, "Steady-state time by length", transform=d_axes[0].transAxes, fontsize=11, fontweight="bold")
    e_axes[0, 0].text(0.0, 1.24, "Counter-cue duration trajectories, full model β=1×", transform=e_axes[0, 0].transAxes, fontsize=11, fontweight="bold")
    e_summary_axes[0].text(0.0, 1.24, "Shortest cue producing steady sign-flipped repolarization", transform=e_summary_axes[0].transAxes, fontsize=10, fontweight="bold")

    color_axis = figure.add_axes((0.975, 0.54, 0.008, 0.30))
    colorbar = figure.colorbar(
        mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(1.0, 39.0), cmap="viridis"),
        cax=color_axis, orientation="vertical",
    )
    colorbar.set_label("Cell length (µm)")
    if FIGURE_NOTE:
        figure.text(0.5, 0.985, FIGURE_NOTE, ha="center", va="top", fontsize=11, fontweight="bold")
    figure.savefig(OUTPUT / f"{OUTPUT_STEM}.png", dpi=300)
    figure.savefig(OUTPUT / f"{OUTPUT_STEM}.pdf")
    plt.close(figure)

    _write_csv(OUTPUT / "figure2_v2_steady_state_metrics.csv", steady_a + steady_c)
    _write_csv(OUTPUT / "figure2_v2_minimum_countercue_duration.csv", minimum_durations)
    _write_csv(OUTPUT / "figure2_v2_countercue_repolarization_audit.csv", _countercue_audit_rows(sweep_rows))
    (OUTPUT / "README.md").write_text(
        "# Standalone Figure 2 v2\n\n"
        "Panel A retains only the 160-s single-stimulus experiment and facets the four "
        "historical model arms by row, now simulated through 10,000 s. Markers show the first "
        "sustained steady-state point only for polarized trajectories. Panel B replaces the former "
        "kymographs with polarization time versus length and includes only polarized full-model "
        "Rac and MCA results. "
        f"Panels C-D use intact-model {ACTIVATION_DESCRIPTION}.\n\n"
        "The polarity coordinate is the interpolated relative half-mass position x50/L. "
        "A steady state is called at the first time after the relevant stimulus at which the "
        "coordinate enters and remains within 0.002L of the median over the final 1,000 s. "
        "The event is censored unless at least 1,000 s of confirmation remain. Panel C points "
        "mark this event: circles denote polarized steady states and crosses denote unpolarized "
        "steady states, using a 0.02L displacement from center as the cutoff. Panel D includes "
        "only polarized, uncensored steady states. Panels A, C, and E use compact y-axis ranges "
        "shared by field within each panel.\n\n"
        + (
            "Panel E shows beta=1x trajectories for 150-, 300-, and 600-s mirrored counter-cues "
            + (
                "applied only after Rac and bound MCA have settled. Panel F summarizes a brute-force "
                "50-duration linearly spaced sweep from 10 to 1000 s for beta=0x, 1x, and 2x. "
                if arguments.julia_activation_defaults_linear50 else
                "applied only after Rac and bound MCA have settled. Panel F summarizes a 20-duration "
                "log-spaced sweep from 10 to 1000 s for beta=0x, 1x, and 2x. "
            )
            if arguments.julia_activation_defaults_short_cues or arguments.julia_activation_defaults_linear50 else
            "Panel E shows beta=1x trajectories for 300-, 1200-, and 7200-s mirrored counter-cues "
            "applied only after Rac and bound MCA have settled. Panel F summarizes a 20-duration "
            "log-spaced sweep from 300 to 10000 s for beta=0x, 1x, and 2x. "
        )
        + "A reversal requires the "
        "final field to settle within 0.002L and lie within 0.05L of the reflected pre-cue "
        "polarity on the opposite side of 0.5L. Both pre- and post-cue states must also lie at "
        "least 0.02L from the center so numerical sign changes of unpolarized states are not "
        "counted. Points without a successful duration are omitted and lines break across them.\n"
        + (
            "\nThis parallel rerun changes only the four activation parameters requested: "
            "alpha0=1, beta0=1, alpha=0.2, and beta=0/0.2/0.4 for the 0x/1x/2x arms. "
            "All other model parameters, cue amplitudes, protocols, numerical settings, and "
            "analysis criteria are inherited unchanged from the existing v2 workflow.\n"
            if (
                arguments.julia_activation_defaults
                or arguments.julia_activation_defaults_short_cues
                or arguments.julia_activation_defaults_linear50
            ) else ""
        )
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
