#!/usr/bin/env python3
"""Plot standalone single-stimulus and continuous-score Figure 2 alternatives."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
NAME = "debelly_2026_continuous_polarity_alternatives_v1"
CAMPAIGN = ROOT / "results/campaigns" / NAME
OUTPUT = ROOT / "results/manuscript_wp1_v3/standalone_alternatives"
VARIANTS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
LABELS = {
    "full": "Full",
    "beta_zero": r"$\beta=0$ only",
    "alpha_zero": r"$\alpha=0$ only",
    "local_inhibition": "Local inhibition",
}
COLORS = {
    "full": "#54278f",
    "beta_zero": "#e6550d",
    "alpha_zero": "#3182bd",
    "local_inhibition": "#31a354",
}
BETA_SWEEP = (
    ("beta_0x", 0.0, 0.0),
    ("beta_0p5x", 0.5, 0.27),
    ("beta_1x", 1.0, 0.54),
    ("beta_1p5x", 1.5, 0.81),
    ("beta_2x", 2.0, 1.08),
)


def _style() -> None:
    mpl.rcParams.update({
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 7.2,
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


def _trajectory(row: dict[str, str]) -> dict[str, np.ndarray]:
    path = CAMPAIGN / "tasks" / row["task_id"] / "trajectory.npz"
    with np.load(path, allow_pickle=False) as raw:
        return {key: raw[key] for key in raw.files}


def _field_half_mass_position(payload: dict[str, np.ndarray], field: str) -> np.ndarray:
    """Return x/L where cumulative cell-centred field mass reaches one half."""
    x = np.asarray(payload["x"], dtype=float)
    values = np.asarray(payload[f"field__{field}"], dtype=float)
    dx = float(np.median(np.diff(x)))
    left_edge = float(x[0] - dx / 2.0)
    length = float(x[-1] - x[0] + dx)
    cumulative = np.cumsum(values, axis=1)
    half = 0.5 * cumulative[:, -1]
    crossing = np.argmax(cumulative >= half[:, None], axis=1)
    rows = np.arange(len(values))
    before = np.where(crossing > 0, cumulative[rows, np.maximum(crossing - 1, 0)], 0.0)
    crossing_mass = values[rows, crossing]
    fraction = np.divide(
        half - before,
        crossing_mass,
        out=np.full_like(half, 0.5),
        where=crossing_mass > 0.0,
    )
    position = left_edge + (crossing + np.clip(fraction, 0.0, 1.0)) * dx
    return np.clip((position - left_edge) / length, 0.0, 1.0)


def _save(figure: plt.Figure, stem: str) -> None:
    figure.savefig(OUTPUT / f"{stem}.png", dpi=300)
    figure.savefig(OUTPUT / f"{stem}.pdf")
    plt.close(figure)


def _latency_figure(rows: list[dict[str, str]]) -> None:
    selected = [row for row in rows if row["assay"] == "single_stimulus_length_test"]
    figure, axes = plt.subplots(1, 2, figsize=(8.2, 3.2), constrained_layout=True)
    specs = (
        ("polarization_latency_s", "Time from stimulus onset (s)", "Polarization"),
        ("depolarization_latency_s", "Time from stimulus withdrawal (s)", "Depolarization"),
    )
    for axis, (metric, ylabel, title) in zip(axes, specs):
        for variant_index, variant in enumerate(VARIANTS):
            group = sorted(
                (row for row in selected if row["coupling_variant"] == variant),
                key=lambda row: _number(row, "length_um"),
            )
            x = np.asarray([_number(row, "length_um") for row in group])
            y = np.asarray([_number(row, metric) for row in group])
            finite = np.isfinite(y)
            axis.plot(x[finite], y[finite], "o-", ms=3.2, lw=1.4,
                      color=COLORS[variant], label=LABELS[variant])
            if np.any(~finite):
                limit = 160.0 if metric.startswith("polarization") else 6_000.0
                missing_rows = [row for row, keep in zip(group, ~finite) if keep]
                for missing_row in missing_rows:
                    not_polarized = missing_row["cue_kinetics_status"] == "not_polarized"
                    marker = "x" if not_polarized else "v"
                    jitter = (-0.18 + 0.12 * variant_index) if not_polarized else 0.0
                    marker_style = (
                        {"color": COLORS[variant]}
                        if not_polarized else
                        {"facecolors": "none", "edgecolors": COLORS[variant]}
                    )
                    axis.scatter(
                        _number(missing_row, "length_um") + jitter, limit,
                        marker=marker, s=24, linewidths=1.0, **marker_style,
                    )
        axis.axvline(15.0, color="#777777", lw=0.8, ls=":")
        axis.set(xlabel="Cell length (µm)", ylabel=ylabel, title=title)
        axis.grid(alpha=0.18)
    axes[0].legend(frameon=False, ncol=2)
    figure.suptitle("Single 160-s stimulus: polarization and loss after withdrawal", fontweight="bold")
    figure.text(0.5, -0.015, "× never reached strong polarity   ▽ still polarized at 6,000 s",
                ha="center", fontsize=7)
    _save(figure, "Figure2_single_stimulus_latencies_vs_length")
    fields = (
        "task_id", "coupling_variant", "length_um", "cue_kinetics_status",
        "polarization_latency_s", "depolarization_latency_s",
        "orientation_at_withdrawal", "peak_directed_orientation",
    )
    with (OUTPUT / "single_stimulus_latencies_vs_length.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in selected)


def _single_score_figure(rows: list[dict[str, str]]) -> None:
    protocols = (
        ("single_stimulus_length_test_40s", 40.0),
        ("single_stimulus_length_test_80s", 80.0),
        ("single_stimulus_length_test", 160.0),
    )
    column_count = len(VARIANTS) * 2
    figure, axes = plt.subplots(len(protocols), column_count, figsize=(24.0, 8.5),
                                constrained_layout=True, sharex=True, sharey=True)
    norm = mpl.colors.Normalize(vmin=1.0, vmax=39.0)
    cmap = mpl.colormaps["viridis"]
    fields = (("rac", "Rac polarity"), ("mca_bound", "Bound-MCA polarity"))
    export_fields = (
        "assay", "variant", "stimulus_duration_s", "field", "length_um",
        "time_from_stimulus_s", "half_mass_position",
    )
    export_path = OUTPUT / "single_stimulus_half_mass_trajectories.csv"
    with export_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=export_fields)
        writer.writeheader()
        for protocol_index, (assay, duration) in enumerate(protocols):
            for variant_index, variant in enumerate(VARIANTS):
                group = sorted(
                    (
                        row for row in rows
                        if row["assay"] == assay and row["coupling_variant"] == variant
                    ),
                    key=lambda item: _number(item, "length_um"),
                )
                for field_index, (field, field_label) in enumerate(fields):
                    column = variant_index * len(fields) + field_index
                    axis = axes[protocol_index, column]
                    for row in group:
                        length = _number(row, "length_um")
                        payload = _trajectory(row)
                        score = _field_half_mass_position(payload, field)
                        time = np.asarray(payload["time"], dtype=float) - 200.0
                        axis.plot(score, time, color=cmap(norm(length)), lw=1.0, alpha=0.9)
                        writer.writerows(
                            {
                                "assay": assay,
                                "variant": variant,
                                "stimulus_duration_s": duration,
                                "field": field,
                                "length_um": length,
                                "time_from_stimulus_s": float(t),
                                "half_mass_position": float(p),
                            }
                            for t, p in zip(time, score)
                        )
                    axis.axhspan(0, duration, color="#f3d66b", alpha=0.18, lw=0)
                    axis.axhline(0, color="#9a6b00", ls="--", lw=0.8)
                    axis.axhline(duration, color="#9a6b00", ls="--", lw=0.8)
                    axis.axvline(0.5, color="#777777", ls=":", lw=0.8)
                    axis.set(xlim=(0.0, 1.0), ylim=(-200.0, 6_500.0))
                    axis.set_yscale("symlog", linthresh=200.0, linscale=1.0)
                    axis.grid(alpha=0.16)
                    if protocol_index == 0:
                        axis.set_title(f"{LABELS[variant]}\n{field_label}")
                    if protocol_index == len(protocols) - 1:
                        axis.set_xlabel(r"50%-content position, $x_{50}/L$")
            axes[protocol_index, 0].set_ylabel(
                f"{duration:g}-s stimulus\nTime from onset (s)"
            )
    figure.suptitle("Single-stimulus duration comparison across model arms and length",
                    fontweight="bold")
    figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
                    label="Cell length (µm)", shrink=0.82)
    _save(figure, "Figure2_single_stimulus_duration_comparison_all_arms_Rac_MCA_polarity_by_length")


def _competition_all_arms_figure(rows: list[dict[str, str]]) -> None:
    assay = "figure2_B_continuous"
    figure, axes = plt.subplots(
        4, 2, figsize=(8.0, 11.5), constrained_layout=True,
        sharex=True, sharey=True,
    )
    norm = mpl.colors.Normalize(vmin=1.0, vmax=39.0)
    cmap = mpl.colormaps["viridis"]
    exported = []
    fields = (("rac", "Rac polarity"), ("mca_bound", "Bound-MCA polarity"))
    ticks = (-0.5, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5)
    tick_labels = ("−0.5", "−0.1", "−0.01", "0", "0.01", "0.1", "0.5")
    for row_index, variant in enumerate(VARIANTS):
        variant_rows = sorted(
            (
                item for item in rows
                if item["assay"] == assay and item["coupling_variant"] == variant
            ),
            key=lambda row: _number(row, "length_um"),
        )
        for column, (field, field_label) in enumerate(fields):
            axis = axes[row_index, column]
            for row in variant_rows:
                length = _number(row, "length_um")
                payload = _trajectory(row)
                time = np.asarray(payload["time"], dtype=float)
                position = _field_half_mass_position(payload, field)
                displacement = position - 0.5
                axis.plot(displacement, time, color=cmap(norm(length)), lw=1.0, alpha=0.9)
                for t, p, delta in zip(time, position, displacement):
                    exported.append({
                        "assay": assay,
                        "variant": variant,
                        "field": field,
                        "length_um": length,
                        "time_s": float(t),
                        "half_mass_position": float(p),
                        "centered_half_mass_displacement": float(delta),
                    })
            axis.axvline(0.0, color="#777777", ls=":", lw=0.8)
            axis.axhspan(0, 80, color="#f3d66b", alpha=0.18, lw=0)
            axis.axhline(80, color="#9a6b00", ls="--", lw=0.8)
            axis.set_xscale("symlog", linthresh=0.01, linscale=1.0, base=10)
            axis.set_yscale("symlog", linthresh=200.0, linscale=1.0)
            axis.set(xlim=(-0.5, 0.5), ylim=(0.0, 21_000.0))
            axis.set_xticks(ticks, tick_labels)
            axis.grid(alpha=0.16)
            if row_index == 0:
                axis.set_title(field_label)
            if row_index == len(VARIANTS) - 1:
                axis.set_xlabel(r"Centered 50%-content position, $x_{50}/L-0.5$")
        axes[row_index, 0].set_ylabel(f"{LABELS[variant]}\nTime from cue onset (s)")
    figure.suptitle(
        "B alternative: continuous two-cue response across model arms",
        fontweight="bold",
    )
    figure.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
        label="Cell length (µm)", shrink=0.82,
    )
    _save(figure, "Figure2_B_all_arms_Rac_MCA_polarity_by_length")
    _write_half_mass_csv(assay, exported)


def _competition_beta_sweep_figure(rows: list[dict[str, str]]) -> None:
    assay = "figure2_B_beta_sweep"
    figure, axes = plt.subplots(
        len(BETA_SWEEP), 2, figsize=(8.0, 14.0), constrained_layout=True,
        sharex=True, sharey=True,
    )
    norm = mpl.colors.Normalize(vmin=1.0, vmax=39.0)
    cmap = mpl.colormaps["viridis"]
    fields = (("rac", "Rac polarity"), ("mca_bound", "Bound-MCA polarity"))
    ticks = (-0.5, -0.1, -0.01, 0.0, 0.01, 0.1, 0.5)
    tick_labels = ("−0.5", "−0.1", "−0.01", "0", "0.01", "0.1", "0.5")
    export_fields = (
        "assay", "beta_multiplier", "beta_tension_feedback", "field", "length_um",
        "time_s", "half_mass_position", "centered_half_mass_displacement",
    )
    export_path = OUTPUT / "figure2_B_beta_sweep_half_mass_trajectories.csv"
    with export_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=export_fields)
        writer.writeheader()
        for row_index, (variant, multiplier, beta) in enumerate(BETA_SWEEP):
            beta_rows = sorted(
                (
                    item for item in rows
                    if item["assay"] == assay and item["coupling_variant"] == variant
                ),
                key=lambda row: _number(row, "length_um"),
            )
            for column, (field, field_label) in enumerate(fields):
                axis = axes[row_index, column]
                for row in beta_rows:
                    length = _number(row, "length_um")
                    payload = _trajectory(row)
                    time = np.asarray(payload["time"], dtype=float)
                    position = _field_half_mass_position(payload, field)
                    displacement = position - 0.5
                    axis.plot(displacement, time, color=cmap(norm(length)), lw=1.0, alpha=0.9)
                    writer.writerows(
                        {
                            "assay": assay,
                            "beta_multiplier": multiplier,
                            "beta_tension_feedback": beta,
                            "field": field,
                            "length_um": length,
                            "time_s": float(t),
                            "half_mass_position": float(p),
                            "centered_half_mass_displacement": float(delta),
                        }
                        for t, p, delta in zip(time, position, displacement)
                    )
                axis.axvline(0.0, color="#777777", ls=":", lw=0.8)
                axis.axhspan(0, 80, color="#f3d66b", alpha=0.18, lw=0)
                axis.axhline(80, color="#9a6b00", ls="--", lw=0.8)
                axis.set_xscale("symlog", linthresh=0.01, linscale=1.0, base=10)
                axis.set_yscale("symlog", linthresh=200.0, linscale=1.0)
                axis.set(xlim=(-0.5, 0.5), ylim=(0.0, 21_000.0))
                axis.set_xticks(ticks, tick_labels)
                axis.grid(alpha=0.16)
                if row_index == 0:
                    axis.set_title(field_label)
                if row_index == len(BETA_SWEEP) - 1:
                    axis.set_xlabel(r"Centered 50%-content position, $x_{50}/L-0.5$")
            axes[row_index, 0].set_ylabel(
                rf"$\beta={beta:.2f}$ ({multiplier:g}$\times$)" + "\nTime from cue onset (s)"
            )
    figure.suptitle(
        r"B alternative: tension-feedback sweep in the full model ($\beta_\mathrm{default}=0.54$)",
        fontweight="bold",
    )
    figure.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
        label="Cell length (µm)", shrink=0.82,
    )
    _save(figure, "Figure2_B_full_beta_sweep_Rac_MCA_polarity_by_length")


def _brief_conditioning_reversal_figure(rows: list[dict[str, str]]) -> None:
    intervals = (40.0, 80.0, 160.0)
    fields = (("rac", "Rac polarity"), ("mca_bound", "Bound-MCA polarity"))
    figure, axes = plt.subplots(
        3, 2, figsize=(8.0, 9.4), constrained_layout=True,
        sharex=True, sharey=True,
    )
    norm = mpl.colors.Normalize(vmin=1.0, vmax=39.0)
    cmap = mpl.colormaps["viridis"]
    exported = []
    for row_index, interval in enumerate(intervals):
        assay = f"brief_conditioning_reversal_interval_{interval:g}s"
        group = sorted(
            (
                row for row in rows
                if row["assay"] == assay and row["coupling_variant"] == "full"
            ),
            key=lambda item: _number(item, "length_um"),
        )
        counter_start = 40.0 + interval
        counter_stop = counter_start + 40.0
        for column, (field, field_label) in enumerate(fields):
            axis = axes[row_index, column]
            for row in group:
                length = _number(row, "length_um")
                payload = _trajectory(row)
                position = _field_half_mass_position(payload, field)
                time = np.asarray(payload["time"], dtype=float) - 200.0
                axis.plot(position, time, color=cmap(norm(length)), lw=1.0, alpha=0.9)
                for t, p in zip(time, position):
                    exported.append({
                        "assay": assay,
                        "variant": "full",
                        "interstimulus_interval_s": interval,
                        "conditioning_duration_s": 40.0,
                        "counterstimulus_duration_s": 40.0,
                        "field": field,
                        "length_um": length,
                        "time_from_conditioning_onset_s": float(t),
                        "half_mass_position": float(p),
                    })
            axis.axhspan(0.0, 40.0, color="#4c78a8", alpha=0.13, lw=0)
            axis.axhspan(counter_start, counter_stop, color="#f2b134", alpha=0.18, lw=0)
            for event, color in (
                (0.0, "#4c78a8"), (40.0, "#4c78a8"),
                (counter_start, "#a66f00"), (counter_stop, "#a66f00"),
            ):
                axis.axhline(event, color=color, ls="--", lw=0.75)
            axis.axvline(0.5, color="#777777", ls=":", lw=0.8)
            axis.set(xlim=(0.0, 1.0), ylim=(-200.0, 6_500.0))
            axis.set_yscale("symlog", linthresh=200.0, linscale=1.0)
            axis.grid(alpha=0.16)
            if row_index == 0:
                axis.set_title(field_label)
            if row_index == len(intervals) - 1:
                axis.set_xlabel(r"50%-content position, $x_{50}/L$")
        axes[row_index, 0].set_ylabel(
            f"{interval:g}-s interval\nTime from first cue (s)"
        )
    figure.suptitle(
        "D alternative: 40-s conditioning cue → 40-s opposing cue",
        fontweight="bold",
    )
    figure.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes,
        label="Cell length (µm)", shrink=0.82,
    )
    _save(figure, "Figure2_D_full_Rac_MCA_polarity_by_length")
    _write_half_mass_csv("figure2_D_continuous", exported)


def _write_half_mass_csv(name: str, rows: list[dict]) -> None:
    path = OUTPUT / f"{name}_half_mass_trajectories.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    _style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = _read(CAMPAIGN / "task_summary.csv")
    _latency_figure(rows)
    _single_score_figure(rows)
    _competition_all_arms_figure(rows)
    _competition_beta_sweep_figure(rows)
    _brief_conditioning_reversal_figure(rows)
    (OUTPUT / "README.md").write_text(
        "# Standalone continuous-polarity alternatives for Figure 2\n\n"
        "The continuous polarity coordinate is the relative position x50/L at which "
        "half of the relevant field's total content lies on either side. It is computed "
        "from conservative cell-centred field masses with interpolation inside the grid "
        "cell containing the 50% crossing. A symmetric field has x50/L = 0.5; a left- "
        "or right-shifted field has a value below or above 0.5, respectively. No event "
        "classification is used in the continuous panels.\n\n"
        "The single-stimulus comparison starts with 200 s without input, applies a left "
        "outer-20% cue at 2 alpha0 for 40, 80, or 160 s (full for the first half and "
        "cosine-tapered for the second half), and follows withdrawal for 6,000 s at "
        "every odd length from 1 to 39 um. "
        "The separate latency summary retains the existing strong end-orientation gate: "
        "cue-aligned orientation >= 0.25 for 10 s to polarize and absolute orientation "
        "< 0.25 for 50 s to depolarize. Open downward triangles indicate that "
        "the latter event did not occur within the observation window.\n\n"
        "The single-stimulus figure includes all four model arms as horizontally "
        "concatenated facet blocks, while B stacks them vertically and D focuses on "
        "the intact/full model. Each has paired "
        "Rac and bound-MCA panels, "
        "with one color-coded raw trajectory per odd cell length from 1 to 39 um. Both "
        "fields use the same half-mass-position definition. B plots x50/L - 0.5 on a "
        "symmetric-log horizontal axis to expand small displacements around the cell "
        "center; time remains on the vertical symmetric-log axis.\n"
        "\nThe separate B beta-sweep figure keeps all other full-model parameters fixed "
        "and varies rho_mechanical_feedback (tension feedback) over beta = 0, 0.27, "
        "0.54, 0.81, and 1.08, corresponding to 0, 0.5, 1, 1.5, and 2 times the "
        "default beta of 0.54.\n"
        "\nThe revised D protocol uses a 40-s left conditioning cue followed by a "
        "40-s right counterstimulus after a cue-free interval of 40, 80, or 160 s. "
        "It retains the D amplitudes and outer-40% stimulus width: 2 alpha0 for "
        "conditioning and 8 alpha0 for the counterstimulus.\n"
    )
    print(OUTPUT)


if __name__ == "__main__":
    main()
