#!/usr/bin/env python3
"""Analyze the De Belly reversal/tracking input-strength titration."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import statistics
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle
import numpy as np
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d.contract import Trajectory  # noqa: E402


CAMPAIGN_NAME = "debelly_input_titration_v1"
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
PRESTATE_ORDER = (
    "unpolarized",
    "single_front",
    "multifront",
    "polarized_no_resolved_peak",
)
PRESTATE_LABELS = {
    "unpolarized": "U = unpolarized",
    "single_front": "S = one detected Rac maximum",
    "multifront": "M = multiple detected Rac maxima",
    "polarized_no_resolved_peak": "P = polarized, no resolved peak",
}
PRESTATE_COLORS = ("#d8d8d8", "#2878b5", "#e68a2e", "#8f6bb3")
PRESTATE_LETTERS = dict(zip(PRESTATE_ORDER, ("U", "S", "M", "P")))


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


def _median(values: list[float]) -> float | None:
    return None if not values else float(statistics.median(values))


def _clean_rows(rows: list[dict[str, str]]) -> tuple[list[dict], list[dict]]:
    reversal = []
    tracking = []
    for row in rows:
        common = {
            "task_id": row["task_id"],
            "input_peak": _number(row, "input_peak"),
            "coupling_variant": row["coupling_variant"],
            "length_um": _number(row, "length_um"),
        }
        if row["assay"] == "gradient_reversal_titration":
            reversal.append(
                {
                    **common,
                    "reversal_conditioned": _boolean(row, "reversal_conditioned"),
                    "reversal_success": _boolean(row, "reversal_success"),
                    "reversal_status": row["reversal_status"],
                    "orientation_reversal_conditioned": _boolean(
                        row, "orientation_reversal_conditioned"
                    ),
                    "orientation_reversal_success": _boolean(
                        row, "orientation_reversal_success"
                    ),
                    "orientation_reversal_status": row.get(
                        "orientation_reversal_status", ""
                    ),
                    "orientation_reversal_latency_s": _number(
                        row, "orientation_reversal_latency_s"
                    ),
                    "pre_reversal_state": row.get("pre_reversal_state") or None,
                    "pre_reversal_polarized": _boolean(
                        row, "pre_reversal_polarized"
                    ),
                    "pre_reversal_domain_count": _number(
                        row, "pre_reversal_domain_count"
                    ),
                    "pre_reversal_spatial_contrast": _number(
                        row, "pre_reversal_spatial_contrast"
                    ),
                    "pre_reversal_orientation": _number(
                        row, "pre_reversal_orientation"
                    ),
                    "reversal_latency_s": _number(row, "reversal_latency_s"),
                    "old_front_decay_latency_s": _number(
                        row, "old_front_decay_latency_s"
                    ),
                    "final_orientation": _number(row, "final_orientation"),
                    "final_spatial_contrast": _number(
                        row, "final_spatial_contrast"
                    ),
                }
            )
        elif row["assay"] == "dynamic_tracking_titration":
            period = _number(row, "tracking_period_s")
            lag = _number(row, "phase_lag_seconds")
            tracking.append(
                {
                    **common,
                    "period_s": period,
                    "response_amplitude": _number(row, "response_amplitude"),
                    "phase_lag_seconds": lag,
                    "phase_lag_fraction": (
                        None if lag is None or period in {None, 0} else lag / period
                    ),
                    "tracking_correlation": _number(row, "tracking_correlation"),
                    "tracking_reliable": _boolean(row, "tracking_reliable"),
                }
            )
    return reversal, tracking


def _reversal_summary(rows: list[dict]) -> list[dict]:
    groups: dict[tuple[float, str], list[dict]] = {}
    for row in rows:
        groups.setdefault((row["input_peak"], row["coupling_variant"]), []).append(
            row
        )
    summary = []
    for (peak, variant), members in sorted(
        groups.items(), key=lambda item: (item[0][0], VARIANT_ORDER.index(item[0][1]))
    ):
        conditioned = [row for row in members if row["reversal_conditioned"]]
        successful = [row for row in members if row["reversal_success"]]
        latencies = [
            row["reversal_latency_s"]
            for row in successful
            if row["reversal_latency_s"] is not None
        ]
        orientation_conditioned = [
            row for row in members if row["orientation_reversal_conditioned"]
        ]
        orientation_successful = [
            row for row in members if row["orientation_reversal_success"]
        ]
        orientation_latencies = [
            row["orientation_reversal_latency_s"]
            for row in orientation_successful
            if row["orientation_reversal_latency_s"] is not None
        ]
        summary.append(
            {
                "input_peak": peak,
                "coupling_variant": variant,
                "length_count": len(members),
                "conditioned_count": len(conditioned),
                "success_count": len(successful),
                "right_censored_count": sum(
                    row["reversal_status"] == "right_censored" for row in members
                ),
                "not_conditioned_count": sum(
                    row["reversal_status"] == "not_conditioned" for row in members
                ),
                "pre_reversal_unpolarized_count": sum(
                    row["pre_reversal_state"] == "unpolarized" for row in members
                ),
                "pre_reversal_single_front_count": sum(
                    row["pre_reversal_state"] == "single_front" for row in members
                ),
                "pre_reversal_multifront_count": sum(
                    row["pre_reversal_state"] == "multifront" for row in members
                ),
                "pre_reversal_no_resolved_peak_count": sum(
                    row["pre_reversal_state"] == "polarized_no_resolved_peak"
                    for row in members
                ),
                "success_fraction_all_lengths": len(successful) / len(members),
                "success_fraction_conditioned": (
                    None if not conditioned else len(successful) / len(conditioned)
                ),
                "median_reversal_latency_s": _median(latencies),
                "orientation_conditioned_count": len(orientation_conditioned),
                "orientation_success_count": len(orientation_successful),
                "orientation_right_censored_count": sum(
                    row["orientation_reversal_status"] == "right_censored"
                    for row in members
                ),
                "orientation_not_conditioned_count": sum(
                    row["orientation_reversal_status"] == "not_conditioned"
                    for row in members
                ),
                "orientation_success_fraction_all_lengths": (
                    len(orientation_successful) / len(members)
                ),
                "orientation_success_fraction_conditioned": (
                    None
                    if not orientation_conditioned
                    else len(orientation_successful) / len(orientation_conditioned)
                ),
                "median_orientation_reversal_latency_s": _median(
                    orientation_latencies
                ),
            }
        )
    return summary


def _tracking_summary(rows: list[dict]) -> list[dict]:
    groups: dict[tuple[float, float, str], list[dict]] = {}
    for row in rows:
        key = (row["input_peak"], row["period_s"], row["coupling_variant"])
        groups.setdefault(key, []).append(row)
    summary = []
    for (peak, period, variant), members in sorted(
        groups.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            VARIANT_ORDER.index(item[0][2]),
        ),
    ):
        amplitudes = [
            row["response_amplitude"]
            for row in members
            if row["response_amplitude"] is not None
        ]
        lags = [
            row["phase_lag_fraction"]
            for row in members
            if row["phase_lag_fraction"] is not None
            and row["response_amplitude"] is not None
            and row["response_amplitude"] > 1e-8
        ]
        correlations = [
            row["tracking_correlation"]
            for row in members
            if row["tracking_correlation"] is not None
        ]
        summary.append(
            {
                "input_peak": peak,
                "period_s": period,
                "coupling_variant": variant,
                "length_count": len(members),
                "median_response_amplitude": _median(amplitudes),
                "minimum_response_amplitude": min(amplitudes),
                "maximum_response_amplitude": max(amplitudes),
                "median_phase_lag_fraction": _median(lags),
                "median_tracking_correlation": _median(correlations),
                "reliable_count": sum(
                    row["tracking_reliable"] is True for row in members
                ),
            }
        )
    return summary


def _discrimination_summary(
    reversal: list[dict], tracking: list[dict]
) -> list[dict]:
    output = []
    reversal_peaks = sorted({row["input_peak"] for row in reversal})
    for peak in reversal_peaks:
        selected = [row for row in reversal if row["input_peak"] == peak]
        success_counts = [row["orientation_success_count"] for row in selected]
        latencies = [
            row["median_orientation_reversal_latency_s"]
            for row in selected
            if row["median_orientation_reversal_latency_s"] is not None
        ]
        strict_success_counts = [row["success_count"] for row in selected]
        strict_latencies = [
            row["median_reversal_latency_s"]
            for row in selected
            if row["median_reversal_latency_s"] is not None
        ]
        output.append(
            {
                "assay": "gradient_reversal",
                "input_peak": peak,
                "period_s": None,
                "variant_success_count_range": max(success_counts) - min(success_counts),
                "variant_median_latency_range_s": (
                    None if len(latencies) < 2 else max(latencies) - min(latencies)
                ),
                "strict_variant_success_count_range": (
                    max(strict_success_counts) - min(strict_success_counts)
                ),
                "strict_variant_median_latency_range_s": (
                    None
                    if len(strict_latencies) < 2
                    else max(strict_latencies) - min(strict_latencies)
                ),
                "variant_median_amplitude_range": None,
                "relative_median_amplitude_range": None,
                "intact_minus_both_off_median_amplitude": None,
            }
        )
    keys = sorted({(row["input_peak"], row["period_s"]) for row in tracking})
    for peak, period in keys:
        selected = [
            row
            for row in tracking
            if row["input_peak"] == peak and row["period_s"] == period
        ]
        by_variant = {row["coupling_variant"]: row for row in selected}
        amplitudes = [row["median_response_amplitude"] for row in selected]
        spread = max(amplitudes) - min(amplitudes)
        output.append(
            {
                "assay": "dynamic_tracking",
                "input_peak": peak,
                "period_s": period,
                "variant_success_count_range": None,
                "variant_median_latency_range_s": None,
                "strict_variant_success_count_range": None,
                "strict_variant_median_latency_range_s": None,
                "variant_median_amplitude_range": spread,
                "relative_median_amplitude_range": (
                    0.0 if max(amplitudes) <= 1e-12 else spread / max(amplitudes)
                ),
                "intact_minus_both_off_median_amplitude": (
                    by_variant["intact"]["median_response_amplitude"]
                    - by_variant["both_off"]["median_response_amplitude"]
                ),
            }
        )
    return output


def _plot_reversal(summary: list[dict], destination: Path) -> None:
    peaks = sorted({row["input_peak"] for row in summary})
    total = max(row["length_count"] for row in summary)
    conditioned = np.zeros((len(VARIANT_ORDER), len(peaks)))
    success = np.zeros_like(conditioned)
    latency = np.full_like(conditioned, np.nan)
    for row in summary:
        i = VARIANT_ORDER.index(row["coupling_variant"])
        j = peaks.index(row["input_peak"])
        conditioned[i, j] = row["orientation_conditioned_count"]
        success[i, j] = row["orientation_success_count"]
        if row["median_orientation_reversal_latency_s"] is not None:
            latency[i, j] = row["median_orientation_reversal_latency_s"]

    figure, axes = plt.subplots(1, 3, figsize=(15, 5.6), sharey=True)
    panels = (
        (conditioned, "Side-oriented before reversal", "viridis", 0, total),
        (success, "Successful sustained side switches", "Blues", 0, total),
        (latency, "Median latency among successes (s)", "magma", None, None),
    )
    for panel, (matrix, title, cmap, vmin, vmax) in zip(axes, panels):
        image = panel.imshow(matrix, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        color_min = float(np.nanmin(matrix)) if vmin is None else float(vmin)
        color_max = float(np.nanmax(matrix)) if vmax is None else float(vmax)
        mapper = plt.get_cmap(cmap)
        normalizer = plt.Normalize(color_min, color_max)
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if np.isfinite(matrix[i, j]):
                    text = (
                        f"{int(matrix[i, j])}/{total}"
                        if panel is not axes[2]
                        else f"{matrix[i, j]:.0f}"
                    )
                    red, green, blue, _ = mapper(normalizer(matrix[i, j]))
                    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                    panel.text(
                        j,
                        i,
                        text,
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if luminance < 0.45 else "black",
                    )
        panel.set_xticks(range(len(peaks)), [f"{peak:g}" for peak in peaks])
        panel.set_xlabel("Peak Rac-gradient input")
        panel.set_title(title)
        figure.colorbar(image, ax=panel, shrink=0.76)
    axes[0].set_yticks(
        range(len(VARIANT_ORDER)), [VARIANT_LABELS[value] for value in VARIANT_ORDER]
    )
    figure.suptitle(
        "De Belly reversal input titration\n"
        "Reversal is a sustained switch in Rac side orientation; counts span eight lengths"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_tracking(summary: list[dict], destination: Path) -> None:
    peaks = sorted({row["input_peak"] for row in summary})
    periods = sorted({row["period_s"] for row in summary})
    x = np.arange(len(peaks))
    figure, axes = plt.subplots(2, len(periods), figsize=(15, 8), sharex=True)
    for column, period in enumerate(periods):
        for variant in VARIANT_ORDER:
            selected = sorted(
                (
                    row
                    for row in summary
                    if row["period_s"] == period
                    and row["coupling_variant"] == variant
                ),
                key=lambda row: row["input_peak"],
            )
            axes[0, column].plot(
                x,
                [row["median_response_amplitude"] for row in selected],
                "o-",
                color=VARIANT_COLORS[variant],
                label=VARIANT_LABELS[variant],
            )
            axes[1, column].plot(
                x,
                [
                    np.nan
                    if row["median_phase_lag_fraction"] is None
                    else row["median_phase_lag_fraction"]
                    for row in selected
                ],
                "o-",
                color=VARIANT_COLORS[variant],
            )
        axes[0, column].set_title(f"Period {period:g} s")
        axes[0, column].grid(alpha=0.2)
        axes[1, column].grid(alpha=0.2)
        axes[1, column].axhline(0.25, color="0.5", ls=":", lw=0.8)
        axes[1, column].axhline(-0.25, color="0.5", ls=":", lw=0.8)
        axes[1, column].set_xticks(x, [f"{peak:g}" for peak in peaks])
        axes[1, column].set_xlabel("Peak Rac-gradient input")
    axes[0, 0].set_ylabel("Median orientation amplitude")
    axes[1, 0].set_ylabel("Median phase lag / period")
    figure.legend(
        *axes[0, 0].get_legend_handles_labels(),
        loc="lower center",
        ncol=3,
        frameon=False,
    )
    figure.suptitle(
        "De Belly alternating-gradient input titration\n"
        "Medians are across four cell lengths; one deterministic run per condition"
    )
    figure.tight_layout(rect=(0, 0.09, 1, 0.93))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_reversal_prestates(rows: list[dict], destination: Path) -> None:
    peaks = sorted({row["input_peak"] for row in rows})
    lengths = sorted({row["length_um"] for row in rows})
    cmap = ListedColormap(PRESTATE_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(PRESTATE_ORDER) + 0.5), cmap.N)
    figure, axes = plt.subplots(2, 3, figsize=(14, 8.5), sharex=True, sharey=True)
    for axis, variant in zip(axes.flat, VARIANT_ORDER):
        matrix = np.zeros((len(lengths), len(peaks)), dtype=int)
        conditioned = np.zeros_like(matrix, dtype=bool)
        orientation_conditioned = np.zeros_like(matrix, dtype=bool)
        for row in rows:
            if row["coupling_variant"] != variant:
                continue
            i = lengths.index(row["length_um"])
            j = peaks.index(row["input_peak"])
            matrix[i, j] = PRESTATE_ORDER.index(row["pre_reversal_state"])
            conditioned[i, j] = bool(row["reversal_conditioned"])
            orientation_conditioned[i, j] = bool(
                row["orientation_reversal_conditioned"]
            )
        axis.imshow(matrix, aspect="auto", cmap=cmap, norm=norm)
        for i in range(len(lengths)):
            for j in range(len(peaks)):
                state = PRESTATE_ORDER[matrix[i, j]]
                axis.text(
                    j,
                    i,
                    PRESTATE_LETTERS[state],
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                )
                if conditioned[i, j]:
                    axis.add_patch(
                        Rectangle(
                            (j - 0.45, i - 0.45),
                            0.9,
                            0.9,
                            fill=False,
                            edgecolor="black",
                            linewidth=1.4,
                        )
                    )
                elif orientation_conditioned[i, j]:
                    axis.add_patch(
                        Rectangle(
                            (j - 0.45, i - 0.45),
                            0.9,
                            0.9,
                            fill=False,
                            edgecolor="#e66101",
                            linewidth=1.8,
                        )
                    )
        axis.set_title(VARIANT_LABELS[variant])
        axis.set_xticks(range(len(peaks)), [f"{peak:g}" for peak in peaks])
        axis.set_yticks(
            range(len(lengths)),
            [
                f"{length:g}" if float(length).is_integer() else f"{length:.1f}"
                for length in lengths
            ],
        )
        axis.set_xlabel("Peak Rac-gradient input")
    axes[0, 0].set_ylabel("Cell length (µm)")
    axes[1, 0].set_ylabel("Cell length (µm)")
    handles = [
        Patch(color=color, label=PRESTATE_LABELS[state])
        for state, color in zip(PRESTATE_ORDER, PRESTATE_COLORS)
    ]
    handles.append(
        Patch(
            facecolor="none",
            edgecolor="black",
            linewidth=1.4,
            label="Black outline = qualified conditioning state",
        )
    )
    handles.append(
        Patch(
            facecolor="none",
            edgecolor="#e66101",
            linewidth=1.8,
            label="Orange = side-oriented but excluded by one-maximum rule",
        )
    )
    figure.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    figure.suptitle(
        "Rac state immediately before the De Belly gradient reversal\n"
        "Letters report the local-maximum detector; outlines report reversal eligibility"
    )
    figure.tight_layout(rect=(0, 0.10, 1, 0.92))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _profile_peak_indices(
    profile: np.ndarray,
    *,
    relative_prominence: float = 0.08,
    min_distance_fraction: float = 0.12,
) -> np.ndarray:
    """Return the exact maxima used by the common domain-count metric."""
    values = np.asarray(profile, dtype=float)
    amplitude = float(np.max(values) - np.min(values))
    scale = max(float(np.max(np.abs(values))), 1.0)
    if amplitude <= 1e-8 * scale:
        return np.asarray([], dtype=int)
    prominence = max(relative_prominence * amplitude, 1e-12)
    distance = max(1, round(min_distance_fraction * len(values)))
    peaks, _ = find_peaks(values, prominence=prominence, distance=distance)
    selected = list(peaks)
    height_threshold = float(np.min(values) + relative_prominence * amplitude)
    if values[0] > values[1] and values[0] >= height_threshold:
        selected.append(0)
    if values[-1] > values[-2] and values[-1] >= height_threshold:
        selected.append(len(values) - 1)
    return np.asarray(sorted(set(selected)), dtype=int)


def _trajectory(output_root: Path, row: dict) -> Trajectory:
    path = output_root / "tasks" / row["task_id"] / "trajectory.npz"
    if not path.is_file():
        raise FileNotFoundError(f"saved trajectory required for figure: {path}")
    return Trajectory.load(path)


def _selected_reversal_row(
    rows: list[dict], *, input_peak: float, length_um: float, variant: str = "intact"
) -> dict:
    matches = [
        row
        for row in rows
        if row["coupling_variant"] == variant
        and np.isclose(row["input_peak"], input_peak)
        and np.isclose(row["length_um"], length_um)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one reversal row for peak={input_peak}, length={length_um}, "
            f"variant={variant}; found {len(matches)}"
        )
    return matches[0]


def _plot_reversal_profile_examples(
    rows: list[dict], output_root: Path, destination: Path
) -> None:
    peaks = (0.2, 0.5, 1.0)
    lengths = (20.0, 2 ** (1 / 3) * 20.0, 32.0, 40.0, 64.0)
    figure, axes = plt.subplots(
        len(peaks), len(lengths), figsize=(16, 8.2), sharex=True, sharey="row"
    )
    for i, input_peak in enumerate(peaks):
        for j, length in enumerate(lengths):
            row = _selected_reversal_row(
                rows, input_peak=input_peak, length_um=length
            )
            trajectory = _trajectory(output_root, row)
            index = int(np.flatnonzero(trajectory.time < 250.0)[-1])
            profile = trajectory.observable[index]
            x_fraction = trajectory.x / length
            detected = _profile_peak_indices(profile)
            axis = axes[i, j]
            axis.axvspan(0.0, 0.2, color="#4c78a8", alpha=0.08)
            axis.axvspan(0.8, 1.0, color="#e45756", alpha=0.08)
            axis.plot(x_fraction, profile, color="#2f2f2f", lw=1.6)
            if detected.size:
                axis.plot(
                    x_fraction[detected],
                    profile[detected],
                    "o",
                    ms=4.5,
                    color="#e66101",
                )
            axis.grid(alpha=0.15)
            axis.set_title(
                f"L={length:.1f} µm; {len(detected)} max.; "
                f"O={row['pre_reversal_orientation']:+.2f}",
                fontsize=9,
            )
            if j == 0:
                axis.set_ylabel(f"input {input_peak:g}\nactive Rac")
            if i == len(peaks) - 1:
                axis.set_xlabel("position x/L")
    figure.suptitle(
        "Intact De Belly Rac profiles immediately before reversal (t = 248 s)\n"
        "Orange points are maxima counted by the strict domain detector; "
        "shading marks the outer fifths"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_reversal_profile_timecourses(
    rows: list[dict], output_root: Path, destination: Path
) -> None:
    examples = (
        ("intact", 0.2, 32.0),
        ("rac_to_rho_off", 0.2, 32.0),
        ("rho_to_rac_off", 0.2, 32.0),
        ("both_off", 0.2, 32.0),
    )
    requested_times = (248.0, 300.0, 400.0, 550.0, 750.0)
    colors = plt.cm.plasma(np.linspace(0.05, 0.9, len(requested_times)))
    figure, axes = plt.subplots(1, 4, figsize=(17, 4.7), sharey=False)
    for axis, (variant, input_peak, length) in zip(axes, examples):
        row = _selected_reversal_row(
            rows, input_peak=input_peak, length_um=length, variant=variant
        )
        trajectory = _trajectory(output_root, row)
        x_fraction = trajectory.x / length
        for requested, color in zip(requested_times, colors):
            index = int(np.argmin(np.abs(trajectory.time - requested)))
            axis.plot(
                x_fraction,
                trajectory.observable[index],
                lw=1.7,
                color=color,
                label=f"{trajectory.time[index]:g} s",
            )
        latency = row["orientation_reversal_latency_s"]
        latency_label = "unavailable" if latency is None else f"{latency:.0f} s"
        axis.set_title(
            f"{VARIANT_LABELS[variant]}\n"
            f"orientation-reversal latency = {latency_label}"
        )
        axis.set_xlabel("position x/L")
        axis.set_ylabel("active Rac")
        axis.grid(alpha=0.18)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False)
    figure.suptitle(
        "Rac profile evolution at input 0.2 and L = 32 µm across the t = 250 s reversal\n"
        "The conditioning cue favors the left; after reversal it favors the right"
    )
    figure.tight_layout(rect=(0, 0.12, 1, 0.91))
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _plot_reversal_latency_heatmaps(rows: list[dict], destination: Path) -> None:
    peaks = sorted({row["input_peak"] for row in rows})
    lengths = sorted({row["length_um"] for row in rows})
    available = [
        row["orientation_reversal_latency_s"]
        for row in rows
        if row["orientation_reversal_latency_s"] is not None
    ]
    vmin, vmax = min(available), max(available)
    cmap = plt.colormaps["viridis"].copy()
    cmap.set_bad("#c8c8c8")
    figure, axes = plt.subplots(2, 3, figsize=(14.5, 9), sharex=True, sharey=True)
    image = None
    for axis, variant in zip(axes.flat, VARIANT_ORDER):
        matrix = np.full((len(lengths), len(peaks)), np.nan)
        strict = np.zeros_like(matrix, dtype=bool)
        orientation_conditioned = np.zeros_like(matrix, dtype=bool)
        orientation_success = np.zeros_like(matrix, dtype=bool)
        for row in rows:
            if row["coupling_variant"] != variant:
                continue
            i = lengths.index(row["length_um"])
            j = peaks.index(row["input_peak"])
            strict[i, j] = bool(row["reversal_conditioned"])
            orientation_conditioned[i, j] = bool(
                row["orientation_reversal_conditioned"]
            )
            orientation_success[i, j] = bool(row["orientation_reversal_success"])
            if row["orientation_reversal_latency_s"] is not None:
                matrix[i, j] = row["orientation_reversal_latency_s"]
        image = axis.imshow(
            matrix,
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        midpoint = 0.5 * (vmin + vmax)
        for i in range(len(lengths)):
            for j in range(len(peaks)):
                if np.isfinite(matrix[i, j]):
                    color = "white" if matrix[i, j] < midpoint else "black"
                    label = f"{matrix[i, j]:.0f}"
                elif orientation_conditioned[i, j] and not orientation_success[i, j]:
                    color, label = "#444444", "C"
                else:
                    color, label = "#666666", "—"
                axis.text(j, i, label, ha="center", va="center", color=color, fontsize=8)
                if orientation_conditioned[i, j] and not strict[i, j]:
                    axis.add_patch(
                        Rectangle(
                            (j - 0.45, i - 0.45),
                            0.9,
                            0.9,
                            fill=False,
                            edgecolor="#e66101",
                            linewidth=1.8,
                        )
                    )
        axis.set_title(VARIANT_LABELS[variant])
        axis.set_xticks(range(len(peaks)), [f"{peak:g}" for peak in peaks])
        axis.set_yticks(
            range(len(lengths)),
            [
                f"{length:g}" if float(length).is_integer() else f"{length:.1f}"
                for length in lengths
            ],
        )
        axis.set_xlabel("Peak Rac-gradient input")
    axes[0, 0].set_ylabel("Cell length (µm)")
    axes[1, 0].set_ylabel("Cell length (µm)")
    assert image is not None
    colorbar_axis = figure.add_axes((0.90, 0.24, 0.018, 0.52))
    colorbar = figure.colorbar(image, cax=colorbar_axis)
    colorbar.set_label("Sustained orientation-reversal latency (s)")
    handles = [
        Patch(color="#c8c8c8", label="Gray / — = no qualified side orientation"),
        Patch(color="#c8c8c8", label="Gray / C = qualified but right-censored"),
        Patch(
            facecolor="none",
            edgecolor="#e66101",
            linewidth=1.8,
            label="Orange outline = multiple maxima, but a qualified old side",
        ),
    ]
    figure.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    figure.suptitle(
        "De Belly gradient-reversal latency by input, length, and inhibition-arm state\n"
        "Numbers are seconds after the reversal begins; one deterministic run per cell"
    )
    figure.subplots_adjust(
        left=0.08,
        right=0.86,
        bottom=0.13,
        top=0.88,
        wspace=0.15,
        hspace=0.22,
    )
    figure.savefig(destination, dpi=180)
    plt.close(figure)


def _report(
    output_root: Path,
    manifest: dict,
    reversal: list[dict],
    tracking: list[dict],
    discrimination: list[dict],
) -> None:
    collection = json.loads((output_root / "collection_summary.json").read_text())
    reversal_discrimination = [
        row for row in discrimination if row["assay"] == "gradient_reversal"
    ]
    tracking_discrimination = [
        row for row in discrimination if row["assay"] == "dynamic_tracking"
    ]
    best_reversal = max(
        reversal_discrimination,
        key=lambda row: (
            row["variant_success_count_range"],
            row["variant_median_latency_range_s"] or -1,
        ),
    )
    positive_tracking = [row for row in tracking_discrimination if row["input_peak"] > 0]
    best_tracking = max(
        positive_tracking, key=lambda row: row["variant_median_amplitude_range"]
    )
    strict_conditioned_total = sum(row["conditioned_count"] for row in reversal)
    strict_successful_total = sum(row["success_count"] for row in reversal)
    oriented_total = sum(row["orientation_conditioned_count"] for row in reversal)
    orientation_successful_total = sum(
        row["orientation_success_count"] for row in reversal
    )
    best_reversal_rows = {
        row["coupling_variant"]: row
        for row in reversal
        if row["input_peak"] == best_reversal["input_peak"]
    }
    best_tracking_rows = [
        row
        for row in tracking
        if row["input_peak"] == best_tracking["input_peak"]
        and row["period_s"] == best_tracking["period_s"]
    ]
    tracking_low = min(
        best_tracking_rows, key=lambda row: row["median_response_amplitude"]
    )
    tracking_high = max(
        best_tracking_rows, key=lambda row: row["median_response_amplitude"]
    )
    latency_values = [
        row["median_orientation_reversal_latency_s"]
        for row in best_reversal_rows.values()
        if row["median_orientation_reversal_latency_s"] is not None
    ]
    latency_range = None if len(latency_values) < 2 else max(latency_values) - min(latency_values)

    def latency_text(variant: str) -> str:
        value = best_reversal_rows[variant]["median_orientation_reversal_latency_s"]
        return "NA" if value is None else f"{value:.0f} s"

    lines = [
        "# De Belly input-strength titration",
        "",
        "## Design",
        "",
        f"All {collection['success_count']}/{collection['task_count']} tasks completed "
        "successfully; there were no failed, missing, or digest-invalid outputs.",
        "",
        f"All {manifest['task_count']} conditions are deterministic. Peak Rac-gradient "
        "inputs were 0, 0.02, 0.05, 0.1, 0.2, 0.5, and 1.0. Basal Rac activation "
        "is 1.0 and the maximum nominal mechanochemical increment is 0.2.",
        "",
        "Reversal retained the production timing: 250 s of conditioning, a 25-s "
        "cross-fade, and observation through 750 s at eight lengths. Tracking retained "
        "periods 100, 200, and 400 s for six cycles at four lengths.",
        "",
        "Counts across lengths are a deterministic design grid, not replicate-based "
        "probability estimates.",
        "",
        "## Reversal definition and result",
        "",
        "The primary reversal endpoint is now a sustained change in the sign of the "
        "Rac side-orientation score. It requires a polarized, directionally biased "
        "state immediately before t = 250 s, but does not require the profile to have "
        "exactly one detected local maximum.",
        "",
        "The orientation score is `(mean Rac in rightmost 20% - mean Rac in leftmost "
        "20%) / (|right mean| + |left mean|)`. A prestate qualifies when spatial "
        "contrast is at least 0.1 and `|orientation| ≥ 0.05`. Success is the first "
        "post-reversal time that begins at least 25 s of polarized orientation at "
        "least 0.05 toward the opposite side; latency is measured from t = 250 s.",
        "",
        f"- {orientation_successful_total}/{oriented_total} side-oriented runs made a "
        "sustained switch to the opposite side. Gray cells in the latency heatmap are "
        "runs without a qualified old-side orientation (or, if marked C, a censored "
        "switch).",
        f"- The older, stricter one-maximum diagnostic reports "
        f"{strict_successful_total}/{strict_conditioned_total}. It excludes "
        f"{oriented_total - strict_conditioned_total} otherwise side-oriented states. "
        "The profile panels show why this matters: several profiles have two detected "
        "local maxima clustered on the same side, so peak count alone does not make "
        "the old side ambiguous.",
        f"- The largest across-arm median-latency spread occurs at input "
        f"{best_reversal['input_peak']:g}"
        + ("." if latency_range is None else f" ({latency_range:.0f} s): ")
        + ", ".join(
            f"{VARIANT_LABELS[variant]} {latency_text(variant)}"
            for variant in VARIANT_ORDER
        )
        + ". These are deterministic differences across the length grid, not "
        "uncertainty estimates.",
        "",
        "## Tracking and inhibition-arm separation",
        "",
        f"- Largest absolute tracking-amplitude spread: input "
        f"{best_tracking['input_peak']:g}, period {best_tracking['period_s']:g} s, "
        f"with a median-amplitude range of "
        f"{best_tracking['variant_median_amplitude_range']:.4f} across variants "
        f"({VARIANT_LABELS[tracking_low['coupling_variant']]} "
        f"{tracking_low['median_response_amplitude']:.4f}; "
        f"{VARIANT_LABELS[tracking_high['coupling_variant']]} "
        f"{tracking_high['median_response_amplitude']:.4f}). Titrating downward did "
        "not improve absolute tracking separation.",
        "",
        "At the original tracking input of 1.0, all conditions should be interpreted "
        "from continuous amplitude and lag rather than only the thresholded pass label.",
        "",
        "## What moves in this model?",
        "",
        "The De Belly model is mechanochemically coupled to internal deformation, but "
        "it is not a whole-cell migration model. It evolves membrane displacement, "
        "cortical velocity, and tension on a fixed reference interval. Rac drives a "
        "protrusive left-boundary velocity, Rho drives cortical active stress, and MCA "
        "is transported by those mechanical velocities. The right boundary is anchored "
        "and there is no center-of-mass translation or freely moving cell body. Rac and "
        "Rho themselves use diffusion and reaction in this 1D port rather than whole-cell "
        "advection.",
        "",
        "## Outputs",
        "",
        "- `debelly_reversal_latency_heatmaps.png`: side-switch latency by input, length, and arm state.",
        "- `debelly_reversal_profile_examples.png`: pre-reversal Rac profiles and detected maxima.",
        "- `debelly_reversal_profile_timecourses.png`: representative profiles before and after reversal.",
        "- `debelly_reversal_input_titration.png`: aggregate side-orientation counts and median latency.",
        "- `debelly_reversal_prestate_titration.png`: local-maximum class and metric eligibility.",
        "- `debelly_tracking_input_titration.png`: response amplitude and lag versus input.",
        "- `reversal_titration.csv` and `tracking_titration.csv`: task-level results.",
        "- `reversal_by_peak_variant.csv` and `tracking_by_peak_period_variant.csv`: summaries.",
        "- `arm_discrimination_summary.csv`: between-variant separation metrics.",
        "",
        "The zero-input controls establish the no-forcing floor; they are not reversal "
        "or tracking experiments in a biological sense.",
        "",
    ]
    (output_root / "README.md").write_text("\n".join(lines))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_root",
        type=Path,
        nargs="?",
        default=ROOT / "results" / "campaigns" / CAMPAIGN_NAME,
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json",
    )
    arguments = parser.parse_args(argv)
    table = arguments.output_root / "task_summary.csv"
    if not table.is_file():
        raise SystemExit(f"collect the campaign first: missing {table}")
    manifest = json.loads(arguments.manifest.read_text())
    reversal, tracking = _clean_rows(_read(table))
    reversal_summary = _reversal_summary(reversal)
    tracking_summary = _tracking_summary(tracking)
    discrimination = _discrimination_summary(reversal_summary, tracking_summary)
    _write(arguments.output_root / "reversal_titration.csv", reversal)
    _write(arguments.output_root / "tracking_titration.csv", tracking)
    _write(
        arguments.output_root / "reversal_by_peak_variant.csv", reversal_summary
    )
    _write(
        arguments.output_root / "tracking_by_peak_period_variant.csv",
        tracking_summary,
    )
    _write(
        arguments.output_root / "arm_discrimination_summary.csv", discrimination
    )
    _plot_reversal(
        reversal_summary,
        arguments.output_root / "debelly_reversal_input_titration.png",
    )
    _plot_reversal_prestates(
        reversal,
        arguments.output_root / "debelly_reversal_prestate_titration.png",
    )
    _plot_reversal_profile_examples(
        reversal,
        arguments.output_root,
        arguments.output_root / "debelly_reversal_profile_examples.png",
    )
    _plot_reversal_profile_timecourses(
        reversal,
        arguments.output_root,
        arguments.output_root / "debelly_reversal_profile_timecourses.png",
    )
    _plot_reversal_latency_heatmaps(
        reversal,
        arguments.output_root / "debelly_reversal_latency_heatmaps.png",
    )
    _plot_tracking(
        tracking_summary,
        arguments.output_root / "debelly_tracking_input_titration.png",
    )
    _report(
        arguments.output_root,
        manifest,
        reversal_summary,
        tracking_summary,
        discrimination,
    )
    print(
        f"wrote titration tables, figures, and report to {arguments.output_root}"
    )


if __name__ == "__main__":
    main()
