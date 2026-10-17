#!/usr/bin/env python3
"""Build revised Figure 2 from the approved responsive de Belly campaign."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import FinitePulseSequence, LocalizedPulse  # noqa: E402


NAME = "debelly_2026_responsive_figure2_v2"
CAMPAIGN = ROOT / "results/campaigns" / NAME
MANIFEST = ROOT / "campaigns" / f"{NAME}.json"
REVERSAL_NAME = "debelly_2026_conditioned_reversal_candidate_v1"
REVERSAL_CAMPAIGN = ROOT / "results/campaigns" / REVERSAL_NAME
REVERSAL_MANIFEST = ROOT / "campaigns" / f"{REVERSAL_NAME}.json"
OUTPUT = ROOT / "results/manuscript_wp1_v3"
SOURCE = OUTPUT / "source_data"
DESTINATION = OUTPUT / "Figure2_WP1_mechanochemical.png"
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
VARIANTS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
VARIANT_LABELS = {
    "full": "Full", "beta_zero": r"$\beta=0$ only",
    "alpha_zero": r"$\alpha=0$ only", "local_inhibition": "Local inhibition",
}
COLORS = {
    "full": "#54278f", "beta_zero": "#e6550d",
    "alpha_zero": "#3182bd", "local_inhibition": "#31a354",
}


def _read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def _boolean(value):
    return str(value).lower() == "true"


def _number(row, key):
    value = row.get(key)
    return None if value in (None, "") else float(value)


def _select(rows, tasks, *, assay, length=15.0, variant="full", period=None):
    for row in rows:
        if row["assay"] != assay or row["coupling_variant"] != variant:
            continue
        if not np.isclose(float(row["length_um"]), length):
            continue
        task = tasks[row["task_id"]]
        if period is not None and not np.isclose(float(task["tracking_period_s"]), period):
            continue
        return row, task
    raise RuntimeError((assay, length, variant, period))


def _protocol(task):
    configuration = dict(task["protocol"])
    name = configuration.pop("name")
    if name == "localized_pulse":
        configuration["sides"] = tuple(configuration["sides"])
        return LocalizedPulse(**configuration)
    if name == "finite_pulse_sequence":
        return FinitePulseSequence(pulses=tuple(configuration["pulses"]))
    raise ValueError(name)


def _load_entry(
    rows, tasks, *, assay, label, title, display_stop, period=None,
    campaign=CAMPAIGN, display_start=0.0, time_offset=0.0,
):
    row, task = _select(rows, tasks, assay=assay, period=period)
    path = campaign / "tasks" / row["task_id"] / "trajectory.npz"
    with np.load(path, allow_pickle=False) as raw:
        payload = {key: raw[key] for key in raw.files}
    selected = (payload["time"] >= display_start) & (payload["time"] <= display_stop)
    absolute_time = payload["time"][selected]
    time = absolute_time - time_offset
    x = payload["x"]
    protocol = _protocol(task)
    fields = {
        "cue": np.asarray([protocol(float(at), x, 15.0) for at in absolute_time]),
    }
    fields.update({
        field: np.asarray(payload[f"field__{field}"])[selected]
        for field in FIELDS[1:]
    })
    return {
        "row": row, "task": task, "time": time, "x": x / 15.0,
        "fields": fields, "label": label, "title": title,
    }


def _write_csv(path, rows):
    fields = sorted({key for row in rows for key in row})
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _draw_kymographs(figure, outer_slot, entries):
    inner = [
        GridSpecFromSubplotSpec(1, len(FIELDS), subplot_spec=outer_slot[index], wspace=0.08)
        for index in range(3)
    ]
    axes = [
        [figure.add_subplot(inner[row][0, column]) for column in range(len(FIELDS))]
        for row in range(3)
    ]
    cmaps = {
        "cue": "Greys", "rac": "viridis", "rho": "magma",
        "tension": "coolwarm", "mca_bound": "cividis",
        "mca_unbound": "plasma", "delta_beta_tension": "Reds",
        "delta_alpha_mca": "Blues",
    }
    norms = {}
    for field in FIELDS:
        values = np.concatenate([entry["fields"][field].ravel() for entry in entries])
        if field == "tension":
            limit = max(abs(float(np.min(values))), abs(float(np.max(values))), 1e-12)
            norms[field] = mpl.colors.TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
        else:
            low, high = float(np.min(values)), float(np.max(values))
            norms[field] = mpl.colors.Normalize(low, high if high > low else low + 1e-12)
    images = {}
    for row_index, entry in enumerate(entries):
        stride = max(1, int(np.ceil(len(entry["time"]) / 1_300)))
        time = entry["time"][::stride]
        for column, field in enumerate(FIELDS):
            axis = axes[row_index][column]
            images[field] = axis.imshow(
                entry["fields"][field][::stride].T,
                origin="lower", aspect="auto", interpolation="nearest",
                extent=(time[0], time[-1], 0.0, 1.0),
                cmap=cmaps[field], norm=norms[field],
            )
            axis.set_xlabel("Time (s)")
            if column == 0:
                axis.set_ylabel("Relative position")
                axis.text(
                    -0.43, 1.13, f"{entry['label']}  {entry['title']}",
                    transform=axis.transAxes, fontsize=9, fontweight="bold", va="top",
                )
            else:
                axis.set_yticklabels([])
            if row_index == 0:
                axis.set_title(FIELD_LABELS[field], fontsize=8)
            axis.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4, prune="upper"))
            axis.tick_params(length=2, pad=1)
    for column, field in enumerate(FIELDS):
        color_axis = axes[0][column].inset_axes((0.03, -0.22, 0.94, 0.032))
        colorbar = figure.colorbar(images[field], cax=color_axis, orientation="horizontal")
        colorbar.ax.tick_params(labelsize=5, length=1, pad=1)
    return axes


def _competition_panel(axis, rows, source_rows):
    horizon = 20_000.0
    for variant in VARIANTS:
        selected = sorted(
            [row for row in rows if row["assay"] == "responsive_competing_fronts" and row["coupling_variant"] == variant],
            key=lambda row: float(row["length_um"]),
        )
        successful = [row for row in selected if _boolean(row["strong_competition_success"])]
        axis.plot(
            [float(row["length_um"]) for row in successful],
            [float(row["strong_resolution_time_s"]) for row in successful],
            "o-", color=COLORS[variant], ms=3, lw=1.2, label=VARIANT_LABELS[variant],
        )
        for row in selected:
            status = row["strong_competition_status"]
            if status == "two_poles_not_established":
                marker = "x"
            elif status == "two_poles_unresolved":
                marker = "v"
            elif status == "resolution_not_persistent":
                marker = "D"
            else:
                marker = None
            if marker:
                axis.scatter(
                    float(row["length_um"]), horizon, marker=marker,
                    facecolors="none" if marker != "x" else COLORS[variant],
                    edgecolors=COLORS[variant], color=COLORS[variant], s=21, lw=0.8,
                )
            source_rows.append({
                "panel": "B", "variant": variant, "length_um": float(row["length_um"]),
                "status": status, "success": _boolean(row["strong_competition_success"]),
                "two_pole_time_s": _number(row, "strong_two_pole_establishment_time_s"),
                "two_pole_peak_strength": _number(row, "strong_two_pole_peak_strength"),
                "resolution_time_s": _number(row, "strong_resolution_time_s"),
                "final_rac_contrast": _number(row, "strong_final_rac_contrast"),
                "final_rac_orientation": _number(row, "strong_final_rac_orientation"),
            })
    axis.axvline(15.0, color="#555555", ls=":", lw=0.8)
    axis.set(
        xlabel="Cell length (µm)", ylabel="Strong single-front resolution (s)",
        title="B  Bona fide two-pole resolution",
    )
    axis.set_yscale("log")
    axis.set_ylim(450, 27_000)
    axis.grid(alpha=0.18)
    axis.legend(frameon=False, fontsize=6, ncol=2, loc="upper left")
    axis.text(
        0.99, 0.03, "× no two poles   ▽ unresolved   ◇ nonpersistent",
        transform=axis.transAxes, ha="right", va="bottom", fontsize=5.7,
    )


def _reversal_panel(figure, slot, rows, source_rows):
    lengths = [float(value) for value in range(1, 40, 2)]
    grid = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=slot, height_ratios=(0.95, 1.35), hspace=0.52,
    )
    status_axis = figure.add_subplot(grid[0, 0])
    latency_axis = figure.add_subplot(grid[1, 0])
    status_order = {
        "not_preconditioned": 0,
        "no_causal_reversal": 1,
        "reversal_not_persistent": 2,
        "reversed_persistent": 3,
    }
    status_labels = ("Not conditioned", "Locked", "Transient", "Persistent")
    status_colors = ("#d9d9d9", "#4d4d4d", "#fdae6b", "#31a354")
    matrix = np.full((len(VARIANTS), len(lengths)), np.nan)
    for variant_index, variant in enumerate(VARIANTS):
        selected = sorted(
            [row for row in rows if row["assay"] == "conditioned_finite_reversal_candidate" and row["coupling_variant"] == variant],
            key=lambda row: float(row["length_um"]),
        )
        successful = [row for row in selected if _boolean(row["strong_reversal_success"])]
        for row in successful:
            latency = float(row["strong_reversal_latency_after_withdrawal_s"])
            shown_latency = max(latency, 5.0)
            latency_axis.scatter(
                float(row["length_um"]), shown_latency,
                marker="v" if latency == 0.0 else "o",
                color=COLORS[variant], s=22, linewidths=0.5,
                label=VARIANT_LABELS[variant] if row is successful[0] else None,
                zorder=3,
            )
        for row in selected:
            length = float(row["length_um"])
            status = row["strong_reversal_status"]
            matrix[variant_index, lengths.index(length)] = status_order[status]
            source_rows.append({
                "panel": "D", "variant": variant, "length_um": length,
                "status": status, "success": _boolean(row["strong_reversal_success"]),
                "preconditioned": _boolean(row["strong_reversal_preconditioned"]),
                "latency_after_withdrawal_s": _number(row, "strong_reversal_latency_after_withdrawal_s"),
                "persistence_s": _number(row, "strong_reversal_persistence_s"),
            })
    cmap = mpl.colors.ListedColormap(status_colors)
    norm = mpl.colors.BoundaryNorm(np.arange(-0.5, 4.5, 1.0), cmap.N)
    status_axis.imshow(matrix, origin="upper", aspect="auto", cmap=cmap, norm=norm)
    status_axis.set_xticks(range(0, len(lengths), 2))
    status_axis.set_xticklabels([f"{lengths[index]:g}" for index in range(0, len(lengths), 2)])
    status_axis.set_yticks(range(len(VARIANTS)))
    status_axis.set_yticklabels([VARIANT_LABELS[variant] for variant in VARIANTS], fontsize=6)
    status_axis.set_title(
        "D  Conditioned finite-pulse reversal", loc="left",
        fontweight="bold", y=1.38,
    )
    status_axis.add_patch(mpl.patches.Rectangle(
        (lengths.index(15.0) - 0.5, -0.5), 1, len(VARIANTS),
        fill=False, edgecolor="#d62728", linewidth=1.1,
    ))
    status_axis.tick_params(length=1.5, pad=1)
    legend_handles = [
        mpl.patches.Patch(facecolor=color, edgecolor="none", label=label)
        for color, label in zip(status_colors, status_labels)
    ]
    status_axis.legend(
        handles=legend_handles, frameon=False, fontsize=4.8, ncol=4,
        loc="lower right", bbox_to_anchor=(1.0, 1.03), borderaxespad=0.0,
        handlelength=1.0, columnspacing=0.8,
    )

    latency_axis.axvline(15.0, color="#555555", ls=":", lw=0.8)
    latency_axis.set_yscale("log")
    latency_axis.set_xlim(0, 40)
    latency_axis.set_ylim(4, 12_000)
    latency_axis.set_xticks(range(1, 40, 4))
    latency_axis.set(
        xlabel="Cell length (µm)", ylabel="Latency after withdrawal (s)",
    )
    latency_axis.set_title("Persistent switches only (▼ ≤5 s)", fontsize=6.2, pad=2)
    latency_axis.grid(alpha=0.18, which="both")
    latency_axis.legend(frameon=False, fontsize=5.4, ncol=2, loc="upper left")
    latency_axis.annotate(
        "delayed L=39\n5.97 ks", xy=(39.0, 5_970.0), xytext=(29.5, 2_600.0),
        fontsize=5.2, color=COLORS["full"],
        arrowprops={"arrowstyle": "-", "lw": 0.6, "color": COLORS["full"]},
    )
    return status_axis, latency_axis


def _frequency_panel(figure, slot, rows, source_rows):
    periods = [100.0, 200.0, 400.0, 800.0, 1_600.0]
    lengths = [float(value) for value in range(1, 40, 2)]
    grid = GridSpecFromSubplotSpec(4, 2, subplot_spec=slot, hspace=0.08, wspace=0.07)
    axes = [[figure.add_subplot(grid[i, j]) for j in range(2)] for i in range(4)]
    specs = (
        ("pulse_train_response_amplitude", "magma", 0.0, 1.0, "Amplitude"),
        ("pulse_train_phase_lag_fraction", "coolwarm", -0.2, 0.2, "Phase / T"),
    )
    images = [None, None]
    for row_index, variant in enumerate(VARIANTS):
        selected = [row for row in rows if row["assay"] == "responsive_pulse_train" and row["coupling_variant"] == variant]
        reliable = np.zeros((len(lengths), len(periods)), dtype=bool)
        matrices = [np.full_like(reliable, np.nan, dtype=float) for _ in specs]
        for row in selected:
            i = lengths.index(float(row["length_um"]))
            j = periods.index(float(row["tracking_period_s"]))
            reliable[i, j] = _boolean(row["pulse_train_reliable"])
            for metric_index, (metric, *_unused) in enumerate(specs):
                matrices[metric_index][i, j] = float(row[metric])
            source_rows.append({
                "panel": "F", "variant": variant, "length_um": lengths[i],
                "period_s": periods[j], "reliable": reliable[i, j],
                "tracked_fraction": float(row["pulse_train_tracking_fraction"]),
                "response_amplitude": float(row["pulse_train_response_amplitude"]),
                "phase_lag_fraction": float(row["pulse_train_phase_lag_fraction"]),
                "harmonic_r2": float(row["pulse_train_harmonic_r2"]),
            })
        for metric_index, (_metric, cmap, vmin, vmax, _title) in enumerate(specs):
            axis = axes[row_index][metric_index]
            plot_cmap = mpl.colormaps[cmap].copy()
            plot_cmap.set_bad("#d9d9d9")
            shown = np.ma.masked_where(~reliable, matrices[metric_index])
            images[metric_index] = axis.imshow(
                shown, origin="lower", aspect="auto", interpolation="nearest",
                cmap=plot_cmap, vmin=vmin, vmax=vmax,
            )
            for i in range(len(lengths)):
                for j in range(len(periods)):
                    if not reliable[i, j]:
                        axis.text(j, i, "×", ha="center", va="center", color="#555555", fontsize=5)
            axis.add_patch(mpl.patches.Rectangle(
                (periods.index(400.0) - 0.5, lengths.index(15.0) - 0.5),
                1, 1, fill=False, edgecolor="#d62728", linewidth=1.2,
            ))
            axis.set_xticks(range(len(periods)))
            if row_index == len(VARIANTS) - 1:
                axis.set_xticklabels([f"{value:g}" for value in periods], fontsize=5.5)
                axis.set_xlabel("Cycle period (s)", fontsize=6.5)
            else:
                axis.set_xticklabels([])
            axis.set_yticks(range(0, len(lengths), 2))
            if metric_index == 0:
                axis.set_yticklabels([])
                axis.set_ylabel(
                    VARIANT_LABELS[variant], fontsize=5.7, rotation=0,
                    ha="right", va="center", labelpad=18,
                )
                axis.text(
                    -0.035, 0.98, "39", transform=axis.transAxes,
                    ha="right", va="top", fontsize=4.8,
                )
                axis.text(
                    -0.035, 0.02, "1", transform=axis.transAxes,
                    ha="right", va="bottom", fontsize=4.8,
                )
            else:
                axis.set_yticklabels([])
            axis.tick_params(length=1, pad=1)
    axes[0][0].text(
        0.0, 1.38, "F  Deterministic frequency response",
        transform=axes[0][0].transAxes, fontsize=7.5, fontweight="bold", ha="left",
    )
    axes[0][0].set_title("Orientation amplitude", fontsize=6.8)
    axes[0][1].set_title("Phase lag / cycle", fontsize=6.8)
    for metric_index, (_metric, _cmap, _vmin, _vmax, title) in enumerate(specs):
        color_axis = axes[0][metric_index].inset_axes((1.04, 0.05, 0.035, 0.85))
        colorbar = figure.colorbar(images[metric_index], cax=color_axis)
        colorbar.ax.tick_params(labelsize=4.5, length=1, pad=1)
        colorbar.ax.set_title(title, fontsize=4.5, pad=1)
    return axes


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    collection = json.loads((CAMPAIGN / "collection_summary.json").read_text())
    if not collection.get("complete") or collection.get("invalid_count") or collection.get("failure_count"):
        raise SystemExit(f"campaign is not clean and complete: {collection}")
    manifest = json.loads(MANIFEST.read_text())
    tasks = {task["task_id"]: task for task in manifest["tasks"]}
    rows = _read_csv(CAMPAIGN / "task_summary.csv")
    reversal_collection = json.loads((REVERSAL_CAMPAIGN / "collection_summary.json").read_text())
    if (
        not reversal_collection.get("complete")
        or reversal_collection.get("invalid_count")
        or reversal_collection.get("failure_count")
    ):
        raise SystemExit(f"reversal campaign is not clean and complete: {reversal_collection}")
    reversal_manifest = json.loads(REVERSAL_MANIFEST.read_text())
    reversal_tasks = {task["task_id"]: task for task in reversal_manifest["tasks"]}
    reversal_rows = _read_csv(REVERSAL_CAMPAIGN / "task_summary.csv")
    entries = [
        _load_entry(
            rows, tasks, assay="responsive_competing_fronts", label="A",
            title="Two competing end cues", display_stop=1_500.0,
        ),
        _load_entry(
            reversal_rows, reversal_tasks,
            assay="conditioned_finite_reversal_candidate", label="C",
            title="Condition → finite opposite pulse (switch at 0)",
            campaign=REVERSAL_CAMPAIGN, display_start=9_000.0,
            display_stop=12_000.0, time_offset=10_000.0,
        ),
        _load_entry(
            rows, tasks, assay="responsive_pulse_train", label="E",
            title="Alternating left/right train (cycle 400 s)",
            display_stop=2_600.0, period=400.0,
        ),
    ]

    mpl.rcParams.update({"font.size": 7, "axes.linewidth": 0.7})
    figure = plt.figure(figsize=(24, 9.5), constrained_layout=False)
    outer = GridSpec(
        3, 2, figure=figure, width_ratios=(4.7, 1.75),
        left=0.035, right=0.985, bottom=0.065, top=0.925,
        hspace=0.40, wspace=0.10,
    )
    _draw_kymographs(figure, [outer[i, 0] for i in range(3)], entries)
    source_rows = []
    _competition_panel(figure.add_subplot(outer[0, 1]), rows, source_rows)
    _reversal_panel(figure, outer[1, 1], reversal_rows, source_rows)
    _frequency_panel(figure, outer[2, 1], rows, source_rows)
    figure.suptitle(
        r"Figure 2. Responsive reflection-symmetric extension of the 2026 de Belly model "
        r"($\alpha_0=\beta_0=1.8$, $\alpha=\beta=0.54$)",
        fontsize=13, fontweight="bold",
    )
    figure.savefig(DESTINATION, dpi=300)
    figure.savefig(DESTINATION.with_suffix(".pdf"))
    plt.close(figure)

    _write_csv(SOURCE / "figure2_responsive_metrics.csv", source_rows)
    shutil.copy2(CAMPAIGN / "collection_summary.json", SOURCE / "debelly_2026_responsive_collection_summary.json")
    shutil.copy2(MANIFEST, SOURCE / "debelly_2026_responsive_figure2_manifest.json")
    shutil.copy2(
        REVERSAL_CAMPAIGN / "collection_summary.json",
        SOURCE / "debelly_2026_conditioned_reversal_collection_summary.json",
    )
    shutil.copy2(
        REVERSAL_MANIFEST,
        SOURCE / "debelly_2026_conditioned_reversal_manifest.json",
    )
    np.savez_compressed(
        SOURCE / "figure2_responsive_L15_kymographs.npz",
        **{
            f"{entry['label']}__{key}": value
            for entry in entries
            for key, value in {"time": entry["time"], "x": entry["x"], **entry["fields"]}.items()
        },
    )

    status_counts = {}
    for assay, field in (
        ("responsive_competing_fronts", "strong_competition_status"),
        ("conditioned_finite_reversal_candidate", "strong_reversal_status"),
    ):
        assay_rows = reversal_rows if assay == "conditioned_finite_reversal_candidate" else rows
        status_counts[assay] = {
            variant: {
                status: sum(
                    row["assay"] == assay
                    and row["coupling_variant"] == variant
                    and row[field] == status
                    for row in assay_rows
                )
                for status in sorted({row[field] for row in assay_rows if row["assay"] == assay})
            }
            for variant in VARIANTS
        }
    l15_comp, _ = _select(rows, tasks, assay="responsive_competing_fronts")
    l15_rev, _ = _select(
        reversal_rows, reversal_tasks, assay="conditioned_finite_reversal_candidate",
    )
    audit = [
        "# Revised de Belly Figure 2 audit", "",
        "Model: reflection-symmetric extension of the 2026 de Belly mechanochemical model.",
        "Selected parameters: alpha0=beta0=1.8; alpha=beta=0.54; every other physical and kinetic parameter unchanged.",
        "The parameterization and L=15 protocols were reviewed before the production length sweep.", "",
        f"Base campaign completion: {collection['success_count']}/{collection['task_count']} successful; "
        f"{collection['failure_count']} failed; {collection['invalid_count']} invalid.",
        f"Replacement reversal campaign: {reversal_collection['success_count']}/{reversal_collection['task_count']} successful; "
        f"{reversal_collection['failure_count']} failed; {reversal_collection['invalid_count']} invalid.",
        f"All field audits passed: "
        f"{sum(_boolean(row['field_audit_passed']) for row in rows + reversal_rows)}/"
        f"{len(rows) + len(reversal_rows)}.", "",
        "## L=15 checks", "",
        f"Competition: {l15_comp['strong_competition_status']}; two poles at "
        f"{float(l15_comp['strong_two_pole_establishment_time_s']):g} s; strong resolution at "
        f"{float(l15_comp['strong_resolution_time_s']):g} s; final Rac contrast "
        f"{float(l15_comp['strong_final_rac_contrast']):.3f}.",
        f"Reversal: {l15_rev['strong_reversal_status']}; pre-cue Rac directed orientation "
        f"{float(l15_rev['strong_pre_opposite_rac_directed']):.3f}; post-withdrawal latency "
        f"{float(l15_rev['strong_reversal_latency_after_withdrawal_s']):g} s; persistence "
        f"{float(l15_rev['strong_reversal_persistence_s']):g} s.", "",
        "## Status counts", "",
        "```json", json.dumps(status_counts, indent=2), "```", "",
        "B displays all unresolved and nonpersistent outcomes at the 20,000-s censoring line.",
        "C-D use 2*alpha0 outer-40% left conditioning for 10 ks, followed by an 8*alpha0 outer-40% right pulse for 320 s.",
        "D displays every outcome categorically and uses a log latency axis for persistent switches only; the isolated delayed L=39 event is labeled explicitly.",
        "F is deterministic frequency response. Gray crossed cells fail the predeclared event-tracking/amplitude audit; no probability is reported.",
    ]
    (OUTPUT / "DEBELLY_RESPONSIVE_FIGURE2_AUDIT.md").write_text("\n".join(audit) + "\n")
    print(DESTINATION)
    print(json.dumps(collection, indent=2))


if __name__ == "__main__":
    main()
