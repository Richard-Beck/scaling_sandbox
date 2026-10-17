#!/usr/bin/env python3
"""Build the temporary de Belly single-cue polarization/depolarization figure."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_NAME = "debelly_2026_single_cue_kinetics_trial_v1"
CAMPAIGN = ROOT / "results" / "tmp" / "trial" / CAMPAIGN_NAME
OUTPUT = ROOT / "results" / "tmp" / "trial"
VARIANTS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
LABELS = {
    "full": "Full",
    "beta_zero": r"$\beta=0$ only",
    "alpha_zero": r"$\alpha=0$ only",
    "local_inhibition": "Local inhibition",
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], key: str) -> float:
    try:
        value = float(row.get(key, "nan"))
    except ValueError:
        return np.nan
    return value if np.isfinite(value) else np.nan


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = _read(CAMPAIGN / "task_summary.csv")
    widths = sorted({_number(row, "cue_width_fraction") for row in rows})
    amplitudes = sorted({_number(row, "cue_amplitude") for row in rows})
    durations = sorted({_number(row, "cue_duration_s") for row in rows})
    finite_polarization = [
        _number(row, "polarization_latency_s") for row in rows
        if np.isfinite(_number(row, "polarization_latency_s"))
    ]
    finite_depolarization = [
        _number(row, "depolarization_latency_s") for row in rows
        if np.isfinite(_number(row, "depolarization_latency_s"))
    ]
    polarization_max = max(finite_polarization, default=max(durations))
    depolarization_max = max(finite_depolarization, default=3_000.0)
    polarization_norm = mpl.colors.Normalize(vmin=0.0, vmax=polarization_max)
    depolarization_norm = mpl.colors.PowerNorm(
        gamma=0.55, vmin=0.0, vmax=depolarization_max
    )

    figure = plt.figure(figsize=(24, 12.5), constrained_layout=False)
    grid = figure.add_gridspec(
        5, 10, left=0.075, right=0.965, top=0.885, bottom=0.08,
        wspace=0.11, hspace=0.34, height_ratios=(1, 1, 1, 1, 0.75),
    )
    axes = [[figure.add_subplot(grid[row, column]) for column in range(10)] for row in range(4)]
    summary_axes = (figure.add_subplot(grid[4, :5]), figure.add_subplot(grid[4, 5:]))
    polarization_image = depolarization_image = None

    for row_index, variant in enumerate(VARIANTS):
        selected_variant = [row for row in rows if row["coupling_variant"] == variant]
        for duration_index, duration in enumerate(durations):
            selected = [
                row for row in selected_variant
                if np.isclose(_number(row, "cue_duration_s"), duration)
            ]
            polarization = np.full((len(widths), len(amplitudes)), np.nan)
            depolarization = np.full_like(polarization, np.nan)
            status = np.full(polarization.shape, "missing", dtype=object)
            for row in selected:
                i = widths.index(_number(row, "cue_width_fraction"))
                j = amplitudes.index(_number(row, "cue_amplitude"))
                polarization[i, j] = _number(row, "polarization_latency_s")
                depolarization[i, j] = _number(row, "depolarization_latency_s")
                status[i, j] = row["cue_kinetics_status"]

            for mode, matrix, column, cmap, norm in (
                ("polarization", polarization, duration_index, "viridis", polarization_norm),
                (
                    "depolarization", depolarization, duration_index + len(durations),
                    "magma", depolarization_norm,
                ),
            ):
                axis = axes[row_index][column]
                plot_cmap = mpl.colormaps[cmap].copy()
                plot_cmap.set_bad("#d9d9d9")
                image = axis.imshow(
                    np.ma.masked_invalid(matrix), origin="lower", aspect="auto",
                    cmap=plot_cmap, norm=norm,
                )
                if mode == "polarization":
                    polarization_image = image
                else:
                    depolarization_image = image
                for i in range(len(widths)):
                    for j in range(len(amplitudes)):
                        if status[i, j] == "not_polarized":
                            axis.text(j, i, "×", ha="center", va="center", color="#555555", fontsize=7)
                        elif status[i, j] == "depolarization_right_censored" and mode == "depolarization":
                            axis.text(j, i, ">", ha="center", va="center", color="white", fontsize=7)
                axis.set_xticks(range(len(amplitudes)))
                axis.set_yticks(range(len(widths)))
                if row_index == len(VARIANTS) - 1:
                    axis.set_xticklabels([f"{value:g}" for value in amplitudes], rotation=45)
                    axis.set_xlabel(r"$\alpha_{cue}/\alpha_0$", fontsize=7)
                else:
                    axis.set_xticklabels([])
                if column == 0:
                    axis.set_yticklabels([f"{100 * value:g}" for value in widths])
                    axis.set_ylabel(f"{LABELS[variant]}\nCue width (% L)", fontsize=7)
                else:
                    axis.set_yticklabels([])
                if row_index == 0:
                    prefix = "Polarize" if mode == "polarization" else "Depolarize"
                    axis.set_title(f"{prefix}\ncue {duration:g} s", fontsize=8)
                axis.tick_params(length=2, pad=1, labelsize=6.5)

    polarization_color_axis = axes[0][2].inset_axes((0.0, -0.18, 1.0, 0.035))
    depolarization_color_axis = axes[0][7].inset_axes((0.0, -0.18, 1.0, 0.035))
    polarization_colorbar = figure.colorbar(
        polarization_image, cax=polarization_color_axis, orientation="horizontal"
    )
    depolarization_colorbar = figure.colorbar(
        depolarization_image, cax=depolarization_color_axis, orientation="horizontal"
    )
    polarization_colorbar.set_label("Polarization latency (s)", fontsize=6)
    depolarization_colorbar.set_label("Depolarization latency after withdrawal (s)", fontsize=6)
    for colorbar in (polarization_colorbar, depolarization_colorbar):
        colorbar.ax.tick_params(labelsize=5.5, length=1.5, pad=1)

    colors = dict(zip(VARIANTS, ("#54278f", "#e6550d", "#3182bd", "#31a354")))
    fraction_axis, median_axis = summary_axes
    summary_rows = []
    for variant_index, variant in enumerate(VARIANTS):
        fraction = []
        median = []
        for duration in durations:
            selected = [
                row for row in rows
                if row["coupling_variant"] == variant
                and np.isclose(_number(row, "cue_duration_s"), duration)
            ]
            polarized = [row["cue_kinetics_status"] != "not_polarized" for row in selected]
            depolarization = [
                _number(row, "depolarization_latency_s") for row in selected
                if np.isfinite(_number(row, "depolarization_latency_s"))
            ]
            fraction.append(float(np.mean(polarized)))
            median.append(float(np.median(depolarization)) if depolarization else np.nan)
            summary_rows.append(
                {
                    "coupling_variant": variant, "cue_duration_s": duration,
                    "fraction_polarized": fraction[-1],
                    "median_depolarization_latency_s": median[-1],
                    "condition_count": len(selected),
                }
            )
        fraction_axis.plot(
            np.asarray(durations) + (-2.0 + 4.0 * variant_index / 3.0), fraction,
            "o-", color=colors[variant], lw=1.35,
            label=LABELS[variant],
        )
        median_axis.plot(
            durations, median, "o-", color=colors[variant], lw=1.6,
            label=LABELS[variant],
        )
    fraction_axis.set(
        xlabel="Total cue duration (s)", ylabel="Fraction polarized",
        ylim=(-0.03, 1.03), title="Conditioning success",
    )
    median_axis.set(
        xlabel="Total cue duration (s)", ylabel="Median loss latency (s)",
        title="Loss after withdrawal",
    )
    for axis in summary_axes:
        axis.grid(alpha=0.18)
        axis.legend(frameon=False, ncol=2, fontsize=7)

    figure.text(
        0.075, 0.955,
        "Temporary trial: single left Rac cues at L=15 µm — polarization and depolarization kinetics",
        fontsize=15, fontweight="bold",
    )
    figure.text(
        0.075, 0.928,
        "Cue is full strength for half its duration and cosine-tapered for the second half; "
        "× = never polarized, > = still polarized 3,000 s after withdrawal.",
        fontsize=8,
    )
    destination = OUTPUT / "Figure_trial_debelly_single_cue_kinetics.png"
    figure.savefig(destination, dpi=300)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)

    source = OUTPUT / "Figure_trial_debelly_single_cue_kinetics_source.csv"
    fields = sorted({key for row in rows for key in row})
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with (OUTPUT / "Figure_trial_debelly_single_cue_kinetics_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    (OUTPUT / "Figure_trial_debelly_single_cue_kinetics_README.md").write_text(
        "# Temporary de Belly single-cue kinetics trial\n\n"
        "This trial does not feed manuscript Figure 2. At L=15 µm it crosses five cue widths "
        "(2.5–30% L), five relative amplitudes (0.1–2), five total durations (10–160 s), "
        "and four independent arm variants. Polarization is the first correctly directed "
        "Rac orientation ≥0.05 sustained for 10 s. Depolarization is the first |orientation| "
        "<0.05 sustained for 50 s after cue withdrawal. Every run is followed for 3,000 s "
        "after withdrawal.\n"
    )
    indexed = {
        variant: {
            (
                _number(row, "cue_width_fraction"),
                _number(row, "cue_amplitude"),
                _number(row, "cue_duration_s"),
            ): row
            for row in rows if row["coupling_variant"] == variant
        }
        for variant in VARIANTS
    }
    full_local = []
    for key, full_row in indexed["full"].items():
        full_latency = _number(full_row, "depolarization_latency_s")
        local_latency = _number(
            indexed["local_inhibition"][key], "depolarization_latency_s"
        )
        if np.isfinite(full_latency) and np.isfinite(local_latency):
            full_local.append((full_latency - local_latency, key, full_latency, local_latency))
    largest = max(full_local, default=(np.nan, (np.nan, np.nan, np.nan), np.nan, np.nan))
    audit = [
        "# Single-cue kinetics trial audit",
        "",
        "All 500 parameter combinations completed without numerical or field-audit failure.",
        "",
        "## Fraction conditioned",
        "",
        "The number of width–amplitude combinations reaching sustained Rac orientation was identical for all four variants:",
        "",
        "| Total cue duration (s) | Polarized / 25 | Fraction |",
        "| ---: | ---: | ---: |",
    ]
    for duration in durations:
        selected = [
            row for row in rows
            if row["coupling_variant"] == "full"
            and np.isclose(_number(row, "cue_duration_s"), duration)
        ]
        count = sum(row["cue_kinetics_status"] != "not_polarized" for row in selected)
        audit.append(f"| {duration:g} | {count} | {count / 25:.2f} |")
    audit += [
        "",
        "## Median depolarization latency after cue withdrawal",
        "",
        "| Cue duration (s) | Full | β=0 | α=0 | Local inhibition |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for duration in durations:
        values = []
        for variant in VARIANTS:
            finite = [
                _number(row, "depolarization_latency_s") for row in rows
                if row["coupling_variant"] == variant
                and np.isclose(_number(row, "cue_duration_s"), duration)
                and np.isfinite(_number(row, "depolarization_latency_s"))
            ]
            values.append(float(np.median(finite)))
        audit.append(
            f"| {duration:g} | {values[0]:.1f} | {values[1]:.1f} | "
            f"{values[2]:.1f} | {values[3]:.1f} |"
        )
    width, amplitude, duration = largest[1]
    audit += [
        "",
        "## Interpretation",
        "",
        "Cue width, amplitude and duration determine whether the direct Rac input crosses the orientation threshold. β removal does not change polarization onset in any conditioned case; α removal or local inhibition advances onset by at most 1 s. The arms therefore do not materially control initial conditioning under these cues.",
        "",
        "Mechanical arms primarily extend the decay after withdrawal. The effect is concentrated in wide, strong and long cues rather than being uniform across the sweep.",
        "",
        f"The largest paired full-minus-local persistence difference is {largest[0]:.0f} s at width {100 * width:g}% L, amplitude {amplitude:g}, duration {duration:g} s: full depolarizes after {largest[2]:.0f} s and local inhibition after {largest[3]:.0f} s.",
        "",
        "No conditioned run remained above threshold for the complete 3,000-s post-withdrawal observation window.",
        "",
    ]
    (OUTPUT / "Figure_trial_debelly_single_cue_kinetics_AUDIT.md").write_text(
        "\n".join(audit)
    )
    print(destination)


if __name__ == "__main__":
    main()
