#!/usr/bin/env python3
"""Analyze the common-protocol de Belly campaign and rebuild manuscript Figure 2."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_NAME = "debelly_2026_jilkine_common_v1"
CAMPAIGN = ROOT / "results" / "campaigns" / CAMPAIGN_NAME
MANIFEST = ROOT / "campaigns" / f"{CAMPAIGN_NAME}.json"
OUTPUT = ROOT / "results" / "manuscript_wp1_v3"
SOURCE = OUTPUT / "source_data"
XVALIDATION = ROOT / "results" / "cross_validation" / "debelly_2026_v1"
SYMMETRY_VALIDATION = (
    ROOT / "results" / "validation" / "debelly_2026_reflection_symmetric"
    / "validation.json"
)
MECHANISM_DIAGNOSTICS = (
    ROOT / "results" / "validation" / "debelly_2026_mechanism_activity"
)
VARIANTS = ("full", "beta_zero", "alpha_zero", "local_inhibition")
VARIANT_LABELS = {
    "full": "Full",
    "beta_zero": r"$\beta=0$ only",
    "alpha_zero": r"$\alpha=0$ only",
    "local_inhibition": "Local inhibition",
}
FIELDS = ("rac", "rho", "tension", "mca_bound", "mca_unbound")
FIELD_LABELS = {
    "rac": "Rac", "rho": "Rho", "tension": "Tension",
    "mca_bound": "MCA bound", "mca_unbound": "MCA unbound",
}
STATE_CODES = {
    "unpolarized": 0,
    "single_front": 1,
    "multifront": 2,
    "polarized_no_resolved_peak": 3,
}
STATE_LETTERS = {0: "U", 1: "S", 2: "M", 3: "P"}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _load_payload(task_id: str) -> dict[str, np.ndarray]:
    path = CAMPAIGN / "tasks" / task_id / "trajectory.npz"
    with np.load(path, allow_pickle=False) as payload:
        return {key: payload[key] for key in payload.files}


def _task_protocol(task_id: str) -> dict:
    return json.loads((CAMPAIGN / "tasks" / task_id / "task.json").read_text())["protocol"]


def _cue(protocol: dict, time: float, x: np.ndarray, length: float) -> np.ndarray:
    name = protocol["name"]
    if name == "localized_pulse":
        if time < protocol["start"] or time >= protocol["stop"]:
            return np.zeros_like(x)
        taper = protocol.get("taper", 0.0)
        envelope = 1.0
        if taper and time > protocol["stop"] - taper:
            phase = np.pi * (time - (protocol["stop"] - taper)) / taper
            envelope = 0.5 * (1.0 + np.cos(phase))
        width = protocol["width_fraction"] * length
        values = np.zeros_like(x)
        for side in protocol["sides"]:
            distance = x if side == "left" else length - x
            selected = distance <= width
            shape = np.zeros_like(x)
            shape[selected] = 0.5 * (1 + np.cos(np.pi * distance[selected] / width))
            asymmetry = protocol.get("asymmetry", 0.0)
            scale = 1 + (asymmetry / 2 if side == "left" else -asymmetry / 2)
            values += scale * shape
        return protocol["amplitude"] * envelope * values
    if name == "reversing_gradient":
        left = protocol["slope"] * (length - x)
        right = protocol["slope"] * x
        reversal = protocol["reversal_time"]
        transition = protocol["transition"]
        if time < reversal:
            return left
        if time >= reversal + transition:
            return right
        weight = (time - reversal) / transition
        return (1 - weight) * left + weight * right
    if name == "alternating_gradient":
        phase = np.cos(2 * np.pi * (time - protocol["start"]) / protocol["period"])
        return 0.5 * protocol["slope"] * length * (1 - phase * (2 * x / length - 1))
    raise ValueError(name)


def _orientation(values: np.ndarray, relative_x: np.ndarray) -> np.ndarray:
    front = np.mean(values[:, relative_x <= 0.2], axis=1)
    back = np.mean(values[:, relative_x >= 0.8], axis=1)
    return (back - front) / np.maximum(np.abs(back) + np.abs(front), 1e-12)


def _status_rows(
    manifest: dict,
    summaries: list[dict[str, str]],
    port_validation_passed: bool,
) -> list[dict]:
    by_id = {row["task_id"]: row for row in summaries}
    rows = []
    for task in manifest["tasks"]:
        row = by_id.get(task["task_id"])
        if row is None:
            directory = CAMPAIGN / "tasks" / task["task_id"]
            failure = directory / "failure.json"
            status = "failed" if failure.is_file() else "missing"
            rows.append(
                {
                    "task_id": task["task_id"], "status": status,
                    "assay": task["assay"], "coupling_variant": task["coupling_variant"],
                    "length_um": task["simulation"]["length"],
                    "field_audit_passed": False,
                }
            )
        else:
            item = dict(row)
            item["port_cross_validation_passed"] = port_validation_passed
            if row.get("field_audit_passed") != "True":
                item["display_status"] = "failed_numerical_audit"
            elif not port_validation_passed:
                item["display_status"] = "failed_port_validation"
            else:
                item["display_status"] = "valid"
            rows.append(item)
    return rows


def _cross_validation_results() -> list[dict]:
    reports = []
    for case in ("competition_l5", "reversal_l15", "tracking_l39"):
        path = XVALIDATION / case / "cross_validation.json"
        if not path.is_file():
            reports.append({"case": case, "status": "missing", "passed": False})
            continue
        payload = json.loads(path.read_text())
        reports.append(
            {
                "case": case,
                "status": "complete",
                "passed": bool(payload["acceptance"]["passed"]),
                "initialization_mode": payload.get("initialization_mode", "unreported"),
                "julia_reference_reused": payload.get("julia_reference_reused", False),
                "maximum_trajectory_relative_l2": payload["acceptance"][
                    "maximum_trajectory_relative_l2"
                ],
                "maximum_orientation_rmse": payload["acceptance"][
                    "maximum_orientation_rmse"
                ],
                "criteria": payload["acceptance"]["criteria"],
            }
        )
    return reports


def _diagnostic_window(assay: str, protocol: dict, times: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if assay == "competing_fronts_common":
        selected = (times >= 0) & (times <= 500)
        display = times
    elif assay == "gradient_reversal_common":
        reversal = protocol["reversal_time"]
        selected = (times >= reversal - 200) & (times <= reversal + 500)
        display = times - reversal
    else:
        period = protocol["period"]
        selected = times >= times[-1] - 2 * period
        display = times - times[-1]
    return selected, display


def _diagnostic_figure(
    assay: str,
    rows: list[dict[str, str]],
    destination: Path,
    source_rows: list[dict],
    port_validation_passed: bool,
) -> None:
    candidates = [
        row for row in rows
        if row["assay"] == assay and float(row["length_um"]) == 15.0
        and (assay != "dynamic_tracking_common" or float(row["tracking_period_s"]) == 200.0)
    ]
    candidates.sort(key=lambda row: VARIANTS.index(row["coupling_variant"]))
    figure, axes = plt.subplots(4, 8, figsize=(17.0, 9.0), constrained_layout=True)
    column_labels = ("Cue", "Rac", "Rho", "Tension", "MCA bound", "MCA unbound",
                     "Signed orientation", "Front/back traces")
    for row_index, row in enumerate(candidates):
        payload = _load_payload(row["task_id"])
        protocol = _task_protocol(row["task_id"])
        times = payload["time"]
        x = payload["x"]
        relative_x = (x - x[0]) / (x[-1] - x[0])
        selected, display_time = _diagnostic_window(assay, protocol, times)
        selected_times = display_time[selected]
        cue = np.asarray([_cue(protocol, time, x, 15.0) for time in times])
        image_fields = (cue,) + tuple(payload[f"field__{field}"] for field in FIELDS)
        for column, values in enumerate(image_fields):
            axis = axes[row_index, column]
            cmap = "coolwarm" if column == 3 else "viridis"
            axis.imshow(
                values[selected], aspect="auto", origin="lower", cmap=cmap,
                extent=(0, 1, selected_times[0], selected_times[-1]),
            )
            if row_index == 0:
                axis.set_title(column_labels[column])
            if column == 0:
                axis.set_ylabel(VARIANT_LABELS[row["coupling_variant"]] + "\nTime (s)")
            else:
                axis.set_yticklabels([])
            axis.set_xlabel("x/L" if row_index == 3 else "")
        orientation_axis = axes[row_index, 6]
        trace_axis = axes[row_index, 7]
        for field, color in zip(("rac", "rho", "tension", "mca_bound"),
                                ("#6a3d9a", "#1f78b4", "#e31a1c", "#33a02c")):
            values = payload[f"field__{field}"]
            orientation = _orientation(values, relative_x)
            orientation_axis.plot(selected_times, orientation[selected], color=color,
                                  lw=1.0, label=FIELD_LABELS[field])
            front = np.mean(values[:, relative_x <= 0.2], axis=1)
            back = np.mean(values[:, relative_x >= 0.8], axis=1)
            scale = max(float(np.ptp(np.r_[front[selected], back[selected]])), 1e-12)
            center = float(np.mean(np.r_[front[selected], back[selected]]))
            trace_axis.plot(selected_times, (front[selected] - center) / scale,
                            color=color, lw=0.9)
            trace_axis.plot(selected_times, (back[selected] - center) / scale,
                            color=color, lw=0.9, ls="--")
            for time, front_value, back_value, orientation_value in zip(
                selected_times, front[selected], back[selected], orientation[selected]
            ):
                source_rows.append(
                    {
                        "assay": assay, "variant": row["coupling_variant"],
                        "length_um": 15.0, "tracking_period_s": row.get("tracking_period_s", ""),
                        "task_id": row["task_id"], "time_display_s": float(time),
                        "field": field, "front_mean": float(front_value),
                        "back_mean": float(back_value),
                        "signed_orientation": float(orientation_value),
                    }
                )
        orientation_axis.axhline(0, color="#888888", lw=0.6)
        trace_axis.axhline(0, color="#888888", lw=0.6)
        if row_index == 0:
            orientation_axis.set_title(column_labels[6])
            trace_axis.set_title(column_labels[7])
            orientation_axis.legend(frameon=False, fontsize=6, ncol=2)
        orientation_axis.set_ylim(-1.05, 1.05)
        orientation_axis.set_xlabel("Time (s)" if row_index == 3 else "")
        trace_axis.set_xlabel("Time (s)" if row_index == 3 else "")
        trace_axis.text(0.02, 0.02, "solid: front; dashed: back\nnormalized per field",
                        transform=trace_axis.transAxes, fontsize=6, va="bottom")
    figure.suptitle(
        f"2026 de Belly common-protocol diagnostics at L=15 µm: {assay.replace('_common', '').replace('_', ' ')}"
        + ("" if port_validation_passed else "  [FAILED NATIVE-PORT AUDIT]"),
        fontsize=14, fontweight="bold",
        color="black" if port_validation_passed else "#b2182b",
    )
    figure.savefig(destination, dpi=250)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)


def _main_figure(
    rows: list[dict[str, str]],
    status_rows: list[dict],
    validation_rows: list[dict],
    destination: Path,
) -> None:
    lengths = sorted({float(row["length_um"]) for row in rows if row["assay"] == "competing_fronts_common"})
    figure = plt.figure(figsize=(14.0, 11.0), constrained_layout=True)
    grid = figure.add_gridspec(3, 2, height_ratios=(0.65, 1.0, 1.0))
    architecture = figure.add_subplot(grid[0, :])
    competition_axis = figure.add_subplot(grid[1, 0])
    reversal_axis = figure.add_subplot(grid[1, 1])
    tracking_axis = figure.add_subplot(grid[2, 0])
    audit_axis = figure.add_subplot(grid[2, 1])

    architecture.axis("off")
    architecture.text(0.0, 0.95, "A  2026 de Belly mechanochemical model under the Jilkine comparison protocol",
                      fontsize=13, fontweight="bold", va="top")
    architecture.text(
        0.01, 0.55,
        "Rac → protrusion → tension → Rho                    Rho → active stress → cortical flow → MCA redistribution → Rac",
        fontsize=11,
    )
    architecture.text(
        0.01, 0.15,
        "L=1,3,…,39 µm; Δx=0.2 µm; Δt=0.5 s; MCA density=0.1; outer-10% cue; 0–20 s full + 20–40 s taper; αcue/α0=1",
        fontsize=10,
    )

    competition = [row for row in rows if row["assay"] == "competing_fronts_common"]
    matrix = np.full((len(VARIANTS), len(lengths)), np.nan)
    for row in competition:
        i = VARIANTS.index(row["coupling_variant"])
        j = lengths.index(float(row["length_um"]))
        matrix[i, j] = STATE_CODES[row["endpoint_state"]]
    competition_axis.imshow(matrix, aspect="auto", cmap=ListedColormap(
        ["#d9d9d9", "#3182bd", "#e6550d", "#9ecae1"]), vmin=-0.5, vmax=3.5)
    competition_axis.set_title("B  Two-cue endpoint; competition never establishes")
    competition_axis.set_yticks(range(4), [VARIANT_LABELS[v] for v in VARIANTS])
    competition_axis.set_xticks(range(len(lengths)), [f"{x:g}" for x in lengths], rotation=90)
    competition_axis.set_xlabel("Cell length (µm)")
    for i in range(4):
        for j in range(len(lengths)):
            competition_axis.text(j, i, STATE_LETTERS[int(matrix[i, j])], ha="center", va="center", fontsize=7)
    competition_axis.text(0.0, -0.24, "U unpolarized; S single front; M multifront; P polarized/no resolved peak.\n"
                          "All 80 runs: competition_status=not_established (not zero-latency resolution).",
                          transform=competition_axis.transAxes, fontsize=8, va="top")

    reversal = [row for row in rows if row["assay"] == "gradient_reversal_common"]
    colors = dict(zip(VARIANTS, ("#54278f", "#e6550d", "#3182bd", "#31a354")))
    for variant in VARIANTS:
        selected = sorted((row for row in reversal if row["coupling_variant"] == variant),
                          key=lambda row: float(row["length_um"]))
        x = [float(row["length_um"]) for row in selected]
        y = [float(row["orientation_reversal_latency_s"]) if row["orientation_reversal_status"] == "success" else np.nan for row in selected]
        reversal_axis.plot(x, y, marker="o", ms=4, lw=1.2, color=colors[variant], label=VARIANT_LABELS[variant])
        for xx, row in zip(x, selected):
            if row["orientation_reversal_status"] != "success":
                reversal_axis.scatter(xx, 20_000, marker="x", s=35, color=colors[variant])
    reversal_axis.set_title("C  Rac orientation reversal; failures retained")
    reversal_axis.set_xlabel("Cell length (µm)")
    reversal_axis.set_ylabel("Latency after reversal (s)")
    reversal_axis.set_yscale("log")
    reversal_axis.legend(frameon=False, fontsize=8)
    reversal_axis.text(0.02, 0.02, "× at 20,000 s: not conditioned; no suppressed runs",
                       transform=reversal_axis.transAxes, fontsize=8)

    tracking = [row for row in rows if row["assay"] == "dynamic_tracking_common" and row["coupling_variant"] == "full"]
    periods = sorted({float(row["tracking_period_s"]) for row in tracking})
    reliability = np.zeros((len(periods), len(lengths)))
    amplitude = np.zeros_like(reliability)
    for row in tracking:
        i = periods.index(float(row["tracking_period_s"]))
        j = lengths.index(float(row["length_um"]))
        reliability[i, j] = row["tracking_reliable"] == "True"
        amplitude[i, j] = float(row["response_amplitude"])
    image = tracking_axis.imshow(amplitude, aspect="auto", cmap="magma", origin="lower")
    for i in range(len(periods)):
        for j in range(len(lengths)):
            if not reliability[i, j]:
                tracking_axis.text(j, i, "×", color="white", ha="center", va="center", fontsize=11)
    tracking_axis.set_title("D  Full-model tracking amplitude; audit failures crossed")
    tracking_axis.set_yticks(range(len(periods)), [f"{p:g}" for p in periods])
    tracking_axis.set_ylabel("Period (s)")
    tracking_axis.set_xticks(range(len(lengths)), [f"{x:g}" for x in lengths], rotation=90)
    tracking_axis.set_xlabel("Cell length (µm)")
    figure.colorbar(image, ax=tracking_axis, label="Orientation amplitude", shrink=0.8)

    audit_axis.axis("off")
    audit_axis.set_title("E  Run and field audit", loc="left")
    counts = {
        "valid": sum(row.get("display_status") == "valid" for row in status_rows),
        "failed_audit": sum(str(row.get("display_status", "")).startswith("failed_") for row in status_rows),
        "failed": sum(row.get("status") == "failed" for row in status_rows),
        "missing": sum(row.get("status") == "missing" for row in status_rows),
    }
    mass_errors = [float(row["total_mca_max_relative_error"]) for row in rows if row.get("total_mca_max_relative_error")]
    minima = {field: min(float(row[f"{field}_minimum"]) for row in rows if row.get(f"{field}_minimum")) for field in FIELDS}
    calibration = next(row for row in rows if row["assay"] == "cue_response_calibration")
    calibration_passed = (
        float(calibration["rac_proximal_peak_change"]) >= 0.05
        and float(calibration["maximum_spatial_contrast"]) >= 0.1
    )
    lines = [
        f"Valid: {counts['valid']} / {len(status_rows)}",
        f"Failed audit: {counts['failed_audit']}; numerical failures: {counts['failed']}; missing: {counts['missing']}",
        f"Maximum total-MCA relative error: {max(mass_errors):.3e}",
        "Global field minima:",
    ] + [f"  {FIELD_LABELS[field]}: {minima[field]:.5g}" for field in FIELDS]
    lines += [
        "", "Independent interventions:",
        "  β=0 changes only tension-dependent Rho activation",
        "  α=0 changes only MCA-dependent Rac activation",
        "  local inhibition sets vCTE=α=β=σa0=0",
        "", (
            f"L=15 cue gate: proximal Rac Δpeak = "
            f"{float(calibration['rac_proximal_peak_change']):.3f}; "
            f"{'PASS' if calibration_passed else 'FAIL'}"
        ),
        "", "Julia/Gridap validation (Julia warm-up IC):",
    ] + [
        f"  {row['case']}: {'PASS' if row['passed'] else 'FAIL'}"
        for row in validation_rows
    ] + [
        "Production IC: homogeneous algebraic fixed point.",
    ]
    audit_axis.text(0.02, 0.95, "\n".join(lines), va="top", family="monospace", fontsize=9)
    port_validation_passed = all(row["passed"] for row in validation_rows)
    figure.suptitle(
        "Figure 2. 2026 de Belly model under the exact Jilkine comparison protocol"
        + ("" if port_validation_passed else "  [FAILED NATIVE-PORT AUDIT]"),
        fontsize=15,
        fontweight="bold",
        color="black" if port_validation_passed else "#b2182b",
    )
    figure.savefig(destination, dpi=300)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)


def _main_figure_symmetric(
    rows: list[dict[str, str]],
    validation_rows: list[dict],
    extension_validation_passed: bool,
    destination: Path,
) -> None:
    """Figure-1-styled summary plus intact L=15 response kymographs."""
    mpl.rcParams.update(
        {
            "font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 9,
            "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
            "legend.fontsize": 7.0, "axes.spines.top": False,
            "axes.spines.right": False, "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
        }
    )
    lengths = sorted(
        {float(row["length_um"]) for row in rows if row["assay"] != "cue_response_calibration"}
    )
    colors = dict(zip(VARIANTS, ("#54278f", "#e6550d", "#3182bd", "#31a354")))
    figure = plt.figure(figsize=(22.0, 16.5), constrained_layout=False)
    outer = figure.add_gridspec(
        4, 2, width_ratios=(5.8, 1.7), height_ratios=(1, 1, 1, 1.55),
        left=0.045, right=0.985, bottom=0.055, top=0.945,
        wspace=0.11, hspace=0.46,
    )
    kymograph_axes = []
    for row_index in range(3):
        field_grid = outer[row_index, 0].subgridspec(1, 8, wspace=0.05)
        kymograph_axes.append(
            [figure.add_subplot(field_grid[0, column]) for column in range(8)]
        )
    competition_axis = figure.add_subplot(outer[0, 1])
    reversal_axis = figure.add_subplot(outer[1, 1])
    tracking_grid = outer[2, 1].subgridspec(4, 2, hspace=0.04, wspace=0.08)
    tracking_axes = [
        [figure.add_subplot(tracking_grid[row, column]) for column in range(2)]
        for row in range(4)
    ]
    activity_grid = outer[3, 0].subgridspec(2, 3, hspace=0.12, wspace=0.16)
    activity_axes = [
        [figure.add_subplot(activity_grid[row, column]) for column in range(3)]
        for row in range(2)
    ]
    causal_grid = outer[3, 1].subgridspec(3, 1, hspace=0.55)
    causal_axes = [figure.add_subplot(causal_grid[row, 0]) for row in range(3)]

    def panel(axis, label: str) -> None:
        axis.text(-0.12, 1.06, label, transform=axis.transAxes, fontsize=13,
                  fontweight="bold", va="top")

    def optional_number(row: dict, key: str) -> float:
        try:
            value = float(row.get(key, "nan"))
        except (TypeError, ValueError):
            return np.nan
        return value if np.isfinite(value) else np.nan

    competition = [row for row in rows if row["assay"] == "competing_fronts_common"]
    for variant in VARIANTS:
        selected = sorted(
            (row for row in competition if row["coupling_variant"] == variant),
            key=lambda row: float(row["length_um"]),
        )
        x = np.asarray([float(row["length_um"]) for row in selected])
        y = np.asarray([optional_number(row, "front_establishment_time_s") for row in selected])
        competition_axis.plot(x, y, "o-", lw=1.7, ms=3.2,
                              color=colors[variant], label=VARIANT_LABELS[variant])
        censored = ~np.isfinite(y)
        competition_axis.scatter(
            x[censored], np.full(np.sum(censored), 20_000.0), marker="v",
            facecolor="none", edgecolor=colors[variant], s=24,
        )
    competition_axis.set(
        xlabel="Cell length (µm)", ylabel="Single-front establishment time (s)",
        title="Two transient end cues",
    )
    competition_axis.set_yscale("log")
    competition_axis.grid(alpha=0.18)
    competition_axis.axvline(15, color="black", ls=":", lw=1.0)
    competition_axis.axhline(20_000, color="#d62728", ls="--", lw=0.9)
    competition_axis.legend(frameon=False, fontsize=6.2, loc="upper left")
    competition_axis.text(
        0.02, 0.02,
        "▽ no sustained front by 20,000 s\nNo run established sustained two-front competition",
        transform=competition_axis.transAxes, fontsize=6.2,
    )
    panel(competition_axis, "B")

    reversal = [row for row in rows if row["assay"] == "gradient_reversal_common"]
    for variant in VARIANTS:
        selected = sorted(
            (row for row in reversal if row["coupling_variant"] == variant),
            key=lambda row: float(row["length_um"]),
        )
        x = np.asarray([float(row["length_um"]) for row in selected])
        success = np.asarray([row["orientation_reversal_status"] == "success" for row in selected])
        y = np.asarray([
            optional_number(row, "orientation_reversal_latency_s") if passed else np.nan
            for row, passed in zip(selected, success)
        ])
        reversal_axis.plot(x, y, "o-", lw=1.7, ms=3.2,
                           color=colors[variant], label=VARIANT_LABELS[variant])
        reversal_axis.scatter(x[~success], np.full(np.sum(~success), 20_000.0),
                              marker="x", color=colors[variant], s=26)
    reversal_axis.set(
        xlabel="Cell length (µm)", ylabel="Rac reversal latency (s)",
        title="Persistent-gradient reversal",
    )
    reversal_axis.set_yscale("log")
    reversal_axis.grid(alpha=0.18)
    reversal_axis.axvline(15, color="black", ls=":", lw=1.0)
    for guide in (20, 100):
        reversal_axis.axhline(guide, color="#d62728", ls="--", lw=0.9)
    reversal_axis.legend(frameon=False, fontsize=6.2, loc="upper left")
    panel(reversal_axis, "D")

    tracking = [row for row in rows if row["assay"] == "dynamic_tracking_common"]
    periods = sorted({float(row["tracking_period_s"]) for row in tracking})
    images = [None, None]
    metric_specs = (("magma", 0.0, 0.35, "Amp."),
                    ("coolwarm", 0.0, 0.25, "Phase/T"))
    shown_lengths = (1.0, 15.0, 29.0, 39.0)
    shown_y = [lengths.index(value) for value in shown_lengths if value in lengths]
    for variant_index, variant in enumerate(VARIANTS):
        selected = [row for row in tracking if row["coupling_variant"] == variant]
        amplitude = np.full((len(lengths), len(periods)), np.nan)
        phase = np.full_like(amplitude, np.nan)
        reliable = np.zeros_like(amplitude, dtype=bool)
        for row in selected:
            i = lengths.index(float(row["length_um"]))
            j = periods.index(float(row["tracking_period_s"]))
            amplitude[i, j] = float(row["response_amplitude"])
            phase[i, j] = float(row["phase_lag_seconds"]) / float(row["tracking_period_s"])
            reliable[i, j] = row["tracking_reliable"] == "True"
        for column, (matrix, spec) in enumerate(zip((amplitude, phase), metric_specs)):
            cmap, vmin, vmax, title = spec
            axis = tracking_axes[variant_index][column]
            plot_cmap = mpl.colormaps[cmap].copy()
            plot_cmap.set_bad("#d9d9d9")
            images[column] = axis.imshow(
                np.ma.masked_where(~reliable, matrix), aspect="auto",
                cmap=plot_cmap, vmin=vmin, vmax=vmax,
            )
            for i in range(len(lengths)):
                for j in range(len(periods)):
                    if not reliable[i, j]:
                        axis.text(j, i, "×", ha="center", va="center",
                                  color="#555555", fontsize=5.5)
            axis.add_patch(
                mpl.patches.Rectangle(
                    (periods.index(200.0) - 0.5, lengths.index(15.0) - 0.5),
                    1, 1, fill=False, edgecolor="#d62728", linewidth=1.25,
                )
            )
            axis.set_yticks(shown_y)
            if column == 0:
                axis.set_yticklabels([f"{lengths[i]:g}" for i in shown_y])
                axis.set_ylabel(VARIANT_LABELS[variant], fontsize=5.7, labelpad=6)
            else:
                axis.set_yticklabels([])
            axis.set_xticks(range(len(periods)))
            if variant_index == len(VARIANTS) - 1:
                axis.set_xticklabels([f"{period:g}" for period in periods])
                axis.set_xlabel("Period (s)", fontsize=7)
            else:
                axis.set_xticklabels([])
            if variant_index == 0:
                axis.set_title(title, fontsize=7.2)
    for column, ticks in enumerate(((0, 0.15, 0.30), (0, 0.125, 0.25))):
        color_axis = tracking_axes[0][column].inset_axes((0.87, 0.10, 0.035, 0.72))
        colorbar = figure.colorbar(images[column], cax=color_axis, ticks=ticks)
        colorbar.ax.tick_params(labelsize=4.5, length=1.3, pad=1)
    tracking_axes[0][0].text(
        -0.28, 1.25, "F  Alternating-gradient tracking",
        transform=tracking_axes[0][0].transAxes, fontsize=10, fontweight="bold",
    )

    exported: list[dict] = []
    specifications = (
        ("competing_fronts_common", None, "A", "Two transient end cues"),
        ("gradient_reversal_common", None, "C", "Persistent-gradient reversal"),
        ("dynamic_tracking_common", 200.0, "E", "Alternating-gradient tracking, T=200 s"),
    )
    kymograph_data: list[dict] = []
    for assay, period, label, title in specifications:
        candidates = [
            row for row in rows
            if row["assay"] == assay and row["coupling_variant"] == "full"
            and float(row["length_um"]) == 15.0
            and (period is None or np.isclose(float(row["tracking_period_s"]), period))
        ]
        row = candidates[0]
        payload = _load_payload(row["task_id"])
        protocol = _task_protocol(row["task_id"])
        times, x = payload["time"], payload["x"]
        cue = np.asarray([_cue(protocol, time, x, 15.0) for time in times])
        if assay == "competing_fronts_common":
            selected = times <= 200.0
            display = times
            events = (20.0, 40.0)
        elif assay == "gradient_reversal_common":
            reversal_time = float(protocol["reversal_time"])
            selected = (times >= reversal_time - 200.0) & (times <= reversal_time + 500.0)
            display = times - reversal_time
            events = (0.0, 20.0)
        else:
            selected = times >= times[-1] - float(period)
            display = times - (times[-1] - float(period))
            events = (0.25 * period, 0.5 * period, 0.75 * period)
        selected_time = display[selected]
        selected_cue = cue[selected]
        relative_x = x / 15.0
        fields = {"cue": selected_cue}
        fields.update({field: payload[f"field__{field}"][selected] for field in FIELDS})
        fields.update(
            {
                field: payload[f"field__{field}"][selected]
                for field in ("delta_beta_tension", "delta_alpha_mca")
            }
        )
        metadata = json.loads(str(payload["metadata_json"]))
        parameters = metadata["parameters"]
        diagnostics = {
            "maximum_delta_beta_tension": np.max(
                fields["delta_beta_tension"], axis=1
            ),
            "maximum_delta_alpha_mca": np.max(
                fields["delta_alpha_mca"], axis=1
            ),
            "fraction_tension_above_threshold": np.mean(
                fields["tension"] > parameters["tension_rho_threshold"], axis=1
            ),
            "fraction_mca_below_threshold": np.mean(
                fields["mca_bound"] < parameters["mca_rac_threshold"], axis=1
            ),
            "boundary_separation": payload["field__boundary_separation"][selected, 0],
            "mean_strain": payload["field__mean_strain"][selected, 0],
            "tension_range": payload["field__tension_range"][selected, 0],
        }
        kymograph_data.append(
            {
                "task_id": row["task_id"], "assay": assay, "label": label,
                "title": title, "time": selected_time, "x": relative_x,
                "events": events, "fields": fields, "diagnostics": diagnostics,
            }
        )

    field_order = (
        "cue", "rac", "rho", "tension", "mca_bound", "mca_unbound",
        "delta_beta_tension", "delta_alpha_mca",
    )
    field_titles = {
        "cue": "Cue", **FIELD_LABELS,
        "delta_beta_tension": r"$\Delta\beta_T$",
        "delta_alpha_mca": r"$\Delta\alpha_{\rm MCA}$",
    }
    field_cmaps = {
        "cue": "Greys", "rac": "viridis", "rho": "magma",
        "tension": "coolwarm", "mca_bound": "cividis", "mca_unbound": "plasma",
        "delta_beta_tension": "Reds", "delta_alpha_mca": "Blues",
    }
    field_norms: dict[str, mpl.colors.Normalize] = {}
    for field in field_order:
        combined = np.concatenate(
            [entry["fields"][field].ravel() for entry in kymograph_data]
        )
        if field == "cue":
            field_norms[field] = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
        elif field == "tension":
            limit = max(float(np.max(np.abs(combined))), 1e-12)
            field_norms[field] = mpl.colors.TwoSlopeNorm(
                vmin=-limit, vcenter=0.0, vmax=limit
            )
        else:
            lower, upper = float(np.min(combined)), float(np.max(combined))
            if np.isclose(lower, upper):
                upper = lower + max(abs(lower), 1.0) * 1e-12
            field_norms[field] = mpl.colors.Normalize(vmin=lower, vmax=upper)

    field_images: dict[str, mpl.image.AxesImage] = {}
    for row_index, entry in enumerate(kymograph_data):
        for column, field in enumerate(field_order):
            axis = kymograph_axes[row_index][column]
            values = entry["fields"][field]
            field_images[field] = axis.imshow(
                values.T, origin="lower", aspect="auto", cmap=field_cmaps[field],
                norm=field_norms[field],
                extent=(entry["time"][0], entry["time"][-1], entry["x"][0], entry["x"][-1]),
                interpolation="nearest",
            )
            for event in entry["events"]:
                axis.axvline(event, color="#d62728", ls="--", lw=0.75)
            axis.set_xlabel("Time (s)", fontsize=7)
            if column == 0:
                axis.set_ylabel("Relative position", fontsize=7)
                axis.text(
                    -0.34, 1.10, f"{entry['label']}  {entry['title']}",
                    transform=axis.transAxes, fontsize=10, fontweight="bold", va="top",
                )
            else:
                axis.set_yticklabels([])
            if row_index == 0:
                axis.set_title(field_titles[field], fontsize=8.5)
            axis.tick_params(length=2, pad=1.5)

            for time_index, time in enumerate(entry["time"]):
                for position, value in zip(entry["x"], values[time_index]):
                    exported.append(
                        {
                            "task_id": entry["task_id"], "assay": entry["assay"],
                            "time_display_s": float(time),
                            "relative_position": float(position),
                            "field": field, "value": float(value),
                        }
                    )

    for column, field in enumerate(field_order):
        color_axis = kymograph_axes[0][column].inset_axes((0.02, -0.23, 0.96, 0.035))
        colorbar = figure.colorbar(
            field_images[field], cax=color_axis, orientation="horizontal",
        )
        colorbar.ax.tick_params(labelsize=5.5, length=1.5, pad=1)

    mechanism_rows: list[dict] = []
    for column, entry in enumerate(kymograph_data):
        feedback_axis = activity_axes[0][column]
        fraction_axis = feedback_axis.twinx()
        diagnostics = entry["diagnostics"]
        time = entry["time"]
        feedback_axis.plot(
            time, diagnostics["maximum_delta_beta_tension"],
            color="#cb181d", lw=1.4, label=r"max $\Delta\beta_T$",
        )
        feedback_axis.plot(
            time, diagnostics["maximum_delta_alpha_mca"],
            color="#2171b5", lw=1.4, label=r"max $\Delta\alpha_{MCA}$",
        )
        fraction_axis.plot(
            time, diagnostics["fraction_tension_above_threshold"],
            color="#cb181d", lw=1.0, ls="--", label=r"fraction $T>T_{th}$",
        )
        fraction_axis.plot(
            time, diagnostics["fraction_mca_below_threshold"],
            color="#2171b5", lw=1.0, ls="--", label=r"fraction $M_b<M_{th}$",
        )
        feedback_axis.set_title(entry["title"], fontsize=8)
        feedback_axis.set_ylabel("Maximum feedback increment" if column == 0 else "")
        fraction_axis.set_ylabel("Domain fraction" if column == 2 else "")
        fraction_axis.set_ylim(-0.03, 1.03)
        feedback_axis.grid(alpha=0.16)

        mechanics_axis = activity_axes[1][column]
        tension_axis = mechanics_axis.twinx()
        mechanics_axis.plot(
            time, diagnostics["boundary_separation"], color="#54278f", lw=1.35,
            label=r"$u_R-u_L$",
        )
        tension_axis.plot(
            time, diagnostics["tension_range"], color="#238b45", lw=1.25,
            label="Tension range",
        )
        mechanics_axis.set_xlabel("Time (s)")
        mechanics_axis.set_ylabel("Boundary separation (µm)" if column == 0 else "")
        mechanics_axis.ticklabel_format(axis="y", style="plain", useOffset=False)
        tension_axis.set_ylabel("Tension range" if column == 2 else "")
        mechanics_axis.grid(alpha=0.16)
        maximum_strain = float(np.max(np.abs(diagnostics["mean_strain"])))
        mechanics_axis.text(
            0.02, 0.93, f"max |mean strain|={maximum_strain:.3g}",
            transform=mechanics_axis.transAxes, va="top", fontsize=6.2,
        )
        for event in entry["events"]:
            feedback_axis.axvline(event, color="#d62728", ls=":", lw=0.7)
            mechanics_axis.axvline(event, color="#d62728", ls=":", lw=0.7)
        if column == 0:
            feedback_axis.legend(frameon=False, fontsize=5.5, loc="upper left")
            fraction_axis.legend(frameon=False, fontsize=5.5, loc="upper right")
            mechanics_axis.legend(frameon=False, fontsize=5.5, loc="upper left")
            tension_axis.legend(frameon=False, fontsize=5.5, loc="upper right")
            feedback_axis.text(
                -0.18, 1.30, "G  Mechanical-feedback activity",
                transform=feedback_axis.transAxes, fontsize=10, fontweight="bold",
            )
        for index, at_time in enumerate(time):
            mechanism_rows.append(
                {
                    "case": entry["assay"], "condition": "full",
                    "time_s": float(at_time),
                    **{
                        name: float(values[index])
                        for name, values in diagnostics.items()
                    },
                }
            )

    causal_specs = (
        (
            "rac_pulse_full", "rac_pulse_beta_zero", "rho",
            "One-sided Rac pulse → distal Rho", r"$\beta=0$ only",
        ),
        (
            "rho_pulse_full", "rho_pulse_alpha_zero", "rac",
            "One-sided Rho pulse → distal Rac", r"$\alpha=0$ only",
        ),
    )
    for index, (full_name, ablated_name, field, title, ablated_label) in enumerate(causal_specs):
        axis = causal_axes[index]
        with np.load(MECHANISM_DIAGNOSTICS / f"{full_name}.npz", allow_pickle=False) as raw:
            full_time = np.asarray(raw["time"])
            full_x = np.asarray(raw["x"])
            full_values = np.asarray(raw[f"field__{field}"])
        with np.load(MECHANISM_DIAGNOSTICS / f"{ablated_name}.npz", allow_pickle=False) as raw:
            ablated_values = np.asarray(raw[f"field__{field}"])
        distal = full_x >= 0.8 * 15.0
        full_trace = np.mean(full_values[:, distal], axis=1)
        ablated_trace = np.mean(ablated_values[:, distal], axis=1)
        full_change = full_trace - full_trace[0]
        ablated_change = ablated_trace - ablated_trace[0]
        axis.plot(full_time, full_change, color="#54278f", lw=1.5, label="Full")
        axis.plot(full_time, ablated_change, color="#e6550d", lw=1.4, label=ablated_label)
        axis.axvspan(0, 20, color="#969696", alpha=0.18, lw=0)
        axis.axvspan(20, 40, color="#969696", alpha=0.09, lw=0)
        axis.axhline(0, color="black", lw=0.6)
        axis.set(xlabel="", ylabel=f"Distal Δ{field.title()}", title=title)
        axis.grid(alpha=0.16)
        axis.legend(frameon=False, fontsize=6.2)
        if index == 0:
            axis.text(
                -0.16, 1.38, "H  Causal and persistence tests",
                transform=axis.transAxes, fontsize=10, fontweight="bold",
            )
        for time_index, at_time in enumerate(full_time):
            mechanism_rows.extend(
                (
                    {
                        "case": title, "condition": "full", "time_s": float(at_time),
                        "distal_change": float(full_change[time_index]), "field": field,
                    },
                    {
                        "case": title, "condition": ablated_name, "time_s": float(at_time),
                        "distal_change": float(ablated_change[time_index]), "field": field,
                    },
                )
            )

    weak_axis = causal_axes[2]
    mechanism_summary = json.loads(
        (MECHANISM_DIAGNOSTICS / "summary.json").read_text()
    )
    weak_rows = mechanism_summary["weak_cue_withdrawal"]
    for condition, color, marker in (
        ("full", "#54278f", "o"),
        ("local_inhibition", "#31a354", "s"),
    ):
        selected = sorted(
            (row for row in weak_rows if row["condition"] == condition),
            key=lambda row: row["amplitude"],
        )
        amplitude = [row["amplitude"] for row in selected]
        at_500 = [abs(row["orientation_at_500_s"]) for row in selected]
        at_5000 = [abs(row["orientation_at_5000_s"]) for row in selected]
        label = "Full" if condition == "full" else "Local inhibition"
        weak_axis.plot(
            amplitude, at_500, color=color, marker=marker, lw=1.35,
            label=f"{label}, 500 s",
        )
        weak_axis.plot(
            amplitude, at_5000, color=color, marker=marker, lw=1.0, ls="--",
            alpha=0.8, label=f"{label}, 5,000 s",
        )
    weak_axis.axhline(0.05, color="#d62728", lw=0.8, ls=":", label="orientation threshold")
    weak_axis.set(
        xlabel=r"Rac cue amplitude $\alpha_{cue}/\alpha_0$",
        ylabel="|Rac orientation|",
        title="One-sided weak cue: persistence after withdrawal",
        yscale="log",
    )
    weak_axis.set_ylim(1e-6, 0.3)
    weak_axis.grid(alpha=0.16)
    weak_axis.legend(frameon=False, fontsize=5.4, ncol=2)
    _write_csv(
        SOURCE / "figure2_debelly_2026_mechanism_activity.csv", mechanism_rows
    )
    _write_csv(SOURCE / "figure2_debelly_2026_l15_kymographs.csv", exported)

    validation_passed = (
        all(row["passed"] for row in validation_rows)
        and extension_validation_passed
    )
    figure.suptitle(
        r"Figure 2. Reflection-symmetric extension of the 2026 de Belly model  "
        r"($\alpha_0=3.5$, $\beta_0=4.0$)"
        + ("" if validation_passed else "  [FAILED MODEL AUDIT]"),
        fontsize=14, fontweight="bold",
        color="black" if validation_passed else "#b2182b",
    )
    figure.savefig(destination, dpi=300)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)


def _write_documentation() -> None:
    caption_path = OUTPUT / "CAPTIONS.md"
    existing = caption_path.read_text() if caption_path.exists() else "# WP1 manuscript figure captions\n"
    marker = "\n## Figure 2."
    if marker in existing:
        existing = existing.split(marker, 1)[0].rstrip() + "\n"
    figure2_caption = r"""

## Figure 2. Reflection-symmetric extension of the 2026 de Belly model

For the Jilkine assays, the published Rac–Rho–MCA circuit is retained but the one-sided boundary mechanics are extended by applying the same tension-inhibited Rac protrusion law at both ends. Mean boundary displacement and velocity are subtracted to remove rigid translation without changing strain or tension. This reflection-symmetric extension is not claimed to be the literally released model; the original one-sided mode is retained for Julia/Gridap validation. All panels use odd lengths from 1–39 µm, approximately 0.2-µm spacing, 0.5-s steps, constant MCA density and the common stimulus definitions. Each assay occupies one row. (A, C, E) Full-model L=15 cue, Rac, Rho, tension, bound-MCA, unbound-MCA, tension-to-Rho feedback increment and MCA-to-Rac feedback increment kymographs for the two-cue, reversal and 200-s alternating-gradient protocols, respectively. Values are not panel-normalized; every field uses one common color scale across all three assays, with tension centered at zero. Red dashed lines mark protocol events. (B) Time to sustained single-front establishment after two transient end cues. Open triangles are explicit 20,000-s censored outcomes; no run first established sustained two-front competition. (D) Rac reversal latency under the persistent-gradient protocol. Crosses are unconditioned runs. (F) Alternating-gradient response amplitude and phase lag, styled as in Figure 1I; gray crossed cells fail the tracking audit and red boxes identify L=15 µm, T=200 s. (G) For each L=15 assay, maximum feedback increments and fractions of the domain above the tension or below the bound-MCA threshold, followed by boundary separation, maximum mean strain and tension range. (H) Directional causal tests showing distal Rho after a one-sided Rac pulse with and without β, distal Rac after a one-sided Rho pulse with and without α, and Rac-orientation persistence 500 and 2,000 s after withdrawal of one-sided cues spanning one tenth to full strength. Production uses the algebraic homogeneous initialization. Published-mode Julia validation uses its finite warm-up initialization and the released element-wise input sampling.

## Figures S3–S5. L=15 mechanochemical pathway diagnostics

For the two-cue, reversal, and 200-s alternating-gradient assays, respectively, rows show the full extension, β=0 only, α=0 only, and the authors' local-inhibition control. Columns show the applied Rac cue and Rac, Rho, tension, bound-MCA, and unbound-MCA kymographs, followed by signed field contrasts and normalized left/right traces. These diagnostics localize responses along cue→Rac, Rac→tension, tension→Rho, and Rho→MCA→Rac.
"""
    caption_path.write_text(existing + figure2_caption.lstrip("\n"))
    (OUTPUT / "README.md").write_text(
        """# Near-manuscript-ready WP1 figure package

Figure 2 and Figures S3–S5 use a clearly labeled reflection-symmetric extension of the 2026 de Belly Rac–Rho–MCA mechanochemical model under the exact Jilkine comparison geometry, timing, length set, and phenotype definitions. The released one-sided assay remains available as `published_one_sided` and is used only for native Julia/Gridap cross-validation. The former fixed-physical-width campaign is not part of this package.

All conditions, including non-established competition, failed reversal conditioning, failed tracking audits and numerical audit states, are represented explicitly. Production uses the homogeneous algebraic Rac–Rho fixed point; published-mode cross-validation alone uses the released Julia finite warm-up initialization. Reflection and pathway validation of the extension and Julia parity of the published mode are both required. No validation waiver is applied.
"""
    )
    (OUTPUT / "DEBELLY_REVERSAL_ABLATION_AUDIT.md").write_text(
        """# Reflection-symmetric 2026 de Belly extension: protocol and ablation audit

The current campaign uses direct parameter interventions: `rho_mechanical_feedback=0` for β=0 only; `rac_mechanical_feedback=0` for α=0 only; and `protrusion_speed=rac_mechanical_feedback=rho_mechanical_feedback=active_stress=0` for local inhibition. No composite arm multipliers or 50% variants feed Figure 2.

The extension applies identical Rac-dependent, tension-inhibited protrusion at both endpoints, centers boundary displacement/velocity to remove translation, and adds no Rho-driven boundary contraction. The reversal endpoint reports Rac orientation and separately exports Rho, tension, and bound-MCA contrasts. The released anchored one-sided geometry is not used for Jilkine production assays.
"""
    )


def _write_mechanism_activity_audit(rows: list[dict[str, str]]) -> None:
    activity = _read_csv(SOURCE / "figure2_debelly_2026_mechanism_activity.csv")
    assay_labels = {
        "competing_fronts_common": "Two-cue competition",
        "gradient_reversal_common": "Persistent-gradient reversal",
        "dynamic_tracking_common": "Alternating-gradient tracking (T=200 s)",
    }
    lines = [
        "# Why the de Belly comparison outcomes are weak",
        "",
        "These diagnostics distinguish an inactive mechanical circuit from a circuit that is active but does not dominate the common assays.",
        "",
        "## Mechanical-feedback activity at L=15 µm",
        "",
        "| Assay | max Δβ_T | max Δα_MCA | max fraction T>T_th | max fraction M_b<M_th | max boundary separation (µm) | max |mean strain| | max tension range |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for assay, label in assay_labels.items():
        selected = [row for row in activity if row["case"] == assay]
        maximum = lambda key: max(float(row[key]) for row in selected if row.get(key, "") != "")
        lines.append(
            f"| {label} | {maximum('maximum_delta_beta_tension'):.4g} | "
            f"{maximum('maximum_delta_alpha_mca'):.4g} | "
            f"{maximum('fraction_tension_above_threshold'):.3g} | "
            f"{maximum('fraction_mca_below_threshold'):.3g} | "
            f"{maximum('boundary_separation'):.4g} | "
            f"{maximum('mean_strain'):.4g} | {maximum('tension_range'):.4g} |"
        )

    mechanism = json.loads((MECHANISM_DIAGNOSTICS / "summary.json").read_text())
    rac_path = mechanism["left_rac_to_distal_rho"]
    rho_path = mechanism["left_rho_to_distal_rac"]
    lines += [
        "",
        "## Directional causal tests",
        "",
        f"A one-sided Rac pulse raises distal Rho by {rac_path['full_peak_change']:.4g} in the full model but only {rac_path['beta_zero_peak_change']:.4g} when β=0; the maximum full-minus-β=0 difference is {rac_path['maximum_full_minus_beta_zero']:.4g}.",
        f"A one-sided Rho pulse raises distal Rac by {rho_path['full_peak_change']:.4g} in the full model and not detectably when α=0; the maximum full-minus-α=0 difference is {rho_path['maximum_full_minus_alpha_zero']:.4g}.",
        "Thus both directed mechanical paths are operational.",
        "",
        "## Common-assay outcomes at L=15 µm",
        "",
    ]
    for assay in assay_labels:
        selected = [
            row for row in rows
            if row["assay"] == assay and float(row["length_um"]) == 15.0
            and (assay != "dynamic_tracking_common" or float(row["tracking_period_s"]) == 200.0)
        ]
        lines.append(f"### {assay_labels[assay]}")
        lines.append("")
        if assay == "competing_fronts_common":
            values = {row["coupling_variant"]: row for row in selected}
            lines.append(
                f"The full model establishes a persistent single front at {float(values['full']['front_establishment_time_s']):.0f} s and β=0 at {float(values['beta_zero']['front_establishment_time_s']):.0f} s. α=0 and local inhibition remain unpolarized at 20,000 s. No variant first establishes a sustained two-front state, so competition resolution remains undefined under the strict common definition."
            )
        elif assay == "gradient_reversal_common":
            values = {row["coupling_variant"]: float(row["orientation_reversal_latency_s"]) for row in selected}
            lines.append(
                f"Every variant reverses under the persistent direct cue: full {values['full']:.0f} s, β=0 {values['beta_zero']:.0f} s, α=0 {values['alpha_zero']:.0f} s and local inhibition {values['local_inhibition']:.0f} s. Mechanics changes latency but is not required for the flip."
            )
        else:
            values = {row["coupling_variant"]: row for row in selected}
            lines.append(
                f"Full and local inhibition are close: amplitudes {float(values['full']['response_amplitude']):.4g} versus {float(values['local_inhibition']['response_amplitude']):.4g}, and phase lags {float(values['full']['phase_lag_seconds']):.3g} versus {float(values['local_inhibition']['phase_lag_seconds']):.3g} s. The nearly domain-wide tension feedback is largely common-mode and contributes little orientation information."
            )
        lines.append("")

    weak = mechanism["weak_cue_withdrawal"]
    full_500 = max(abs(row["orientation_at_500_s"]) for row in weak if row["condition"] == "full")
    full_2000 = max(abs(row["orientation_at_2000_s"]) for row in weak if row["condition"] == "full")
    full_5000 = max(abs(row["orientation_at_5000_s"]) for row in weak if row["condition"] == "full")
    local_500 = max(abs(row["orientation_at_500_s"]) for row in weak if row["condition"] == "local_inhibition")
    lines += [
        "## Interpretation",
        "",
        f"After withdrawal of one-sided cues spanning α_cue/α0=0.1–1, the largest full-model Rac orientation grows from {full_500:.4g} at 500 s to {full_2000:.4g} at 2,000 s and {full_5000:.4g} at 5,000 s. Local inhibition instead collapses rapidly ({local_500:.3g} at 500 s). The selected full model therefore converts a sufficiently strong transient cue into slowly developing autonomous polarity.",
        "",
        "Both directed mechanical paths remain operational. In the selected regime, mechanics is required for long-lived polarity after withdrawal and for delayed single-front establishment after the two transient end cues, while strong persistent and alternating cues still drive all variants directly. The separation between those assays is the intended distinction between autonomous memory and cue-following responsiveness.",
        "",
    ]
    (OUTPUT / "DEBELLY_MECHANISM_ACTIVITY_AUDIT.md").write_text("\n".join(lines))


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    for stale in (
        "figure2_competition_hardened.csv",
        "figure2_debelly_field_orientation.csv",
        "figure2_debelly_transition_profiles.csv",
        "figure2_fixed_width_propagation.csv",
        "figure2_reversal_length_scaled.csv",
        "wp1_decisive_collection_summary.json",
    ):
        path = SOURCE / stale
        if path.exists():
            path.unlink()
    manifest = json.loads(MANIFEST.read_text())
    summaries = _read_csv(CAMPAIGN / "task_summary.csv")
    validation_rows = _cross_validation_results()
    port_validation_passed = all(row["passed"] for row in validation_rows)
    extension_validation = json.loads(SYMMETRY_VALIDATION.read_text())
    extension_validation_passed = bool(extension_validation["passed"])
    overall_validation_passed = port_validation_passed and extension_validation_passed
    status_rows = _status_rows(manifest, summaries, overall_validation_passed)
    _write_csv(SOURCE / "figure2_debelly_2026_run_audit.csv", status_rows)
    _write_csv(SOURCE / "figure2_debelly_2026_cross_validation.csv", validation_rows)
    validation_field_rows: list[dict] = []
    for validation in validation_rows:
        path = XVALIDATION / validation["case"] / "cross_validation.json"
        if not path.is_file():
            continue
        payload = json.loads(path.read_text())
        for field, metrics in payload["fields"].items():
            validation_field_rows.append(
                {
                    "case": validation["case"], "field": field,
                    **metrics,
                }
            )
    _write_csv(
        SOURCE / "figure2_debelly_2026_cross_validation_fields.csv",
        validation_field_rows,
    )
    orientation_rows = []
    for row in summaries:
        if row["assay"] == "cue_response_calibration":
            continue
        item = {
            key: row.get(key, "") for key in (
                "task_id", "assay", "coupling_variant", "length_um", "tracking_period_s",
                "status", "field_audit_passed", "endpoint_state", "competition_status",
                "orientation_reversal_status", "orientation_reversal_latency_s",
                "reversal_status", "tracking_reliable", "response_amplitude",
                "phase_lag_seconds", "total_mca_max_relative_error",
            )
        }
        for field in FIELDS:
            item[f"{field}_final_orientation"] = row.get(f"{field}_final_orientation", "")
            item[f"{field}_minimum"] = row.get(f"{field}_minimum", "")
            item[f"{field}_front_final"] = row.get(f"{field}_front_final", "")
            item[f"{field}_back_final"] = row.get(f"{field}_back_final", "")
        orientation_rows.append(item)
    _write_csv(SOURCE / "figure2_debelly_2026_field_summary.csv", orientation_rows)

    _main_figure_symmetric(
        summaries,
        validation_rows,
        extension_validation_passed,
        OUTPUT / "Figure2_WP1_mechanochemical.png",
    )
    diagnostic_source: list[dict] = []
    for assay, name in (
        ("competing_fronts_common", "FigureS3_debelly_2026_competition_diagnostics.png"),
        ("gradient_reversal_common", "FigureS4_debelly_2026_reversal_diagnostics.png"),
        ("dynamic_tracking_common", "FigureS5_debelly_2026_tracking_diagnostics.png"),
    ):
        _diagnostic_figure(
            assay,
            summaries,
            OUTPUT / name,
            diagnostic_source,
            overall_validation_passed,
        )
    _write_csv(SOURCE / "figure2_debelly_2026_l15_front_back.csv", diagnostic_source)
    shutil.copyfile(CAMPAIGN / "collection_summary.json", SOURCE / "debelly_2026_collection_summary.json")
    shutil.copyfile(
        SYMMETRY_VALIDATION,
        SOURCE / "debelly_2026_reflection_symmetric_validation.json",
    )
    _write_documentation()
    _write_mechanism_activity_audit(summaries)
    lines = [
        "# 2026 de Belly common-protocol audit",
        "",
        "Production initialization: `algebraic_fixed_point`.",
        "Julia/Python cross-validation initialization: `julia_warmup` in both implementations.",
        "Julia/Python cross-validation input sampling: Julia element values at `j*dx`; production remains cell-centred.",
        "The existing native Julia reference trajectories were reused; Julia was not rerun.",
        "No validation waiver is used.",
        "Julia-exported tension is compared using its backward-difference interior-node representation; the direct centred-FV tension discrepancy is retained in the field-level source table.",
        "",
        "| Case | Status | max relative L2 | max orientation RMSE |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in validation_rows:
        relative = row.get("maximum_trajectory_relative_l2")
        orientation = row.get("maximum_orientation_rmse")
        lines.append(
            f"| {row['case']} | {'PASS' if row['passed'] else 'FAIL'} | "
            f"{relative:.5g} | {orientation:.5g} |"
            if relative is not None and orientation is not None
            else f"| {row['case']} | MISSING | — | — |"
        )
    native = "PASS" if port_validation_passed else "FAIL"
    symmetric = "PASS" if extension_validation_passed else "FAIL"
    overall = "PASS" if overall_validation_passed else "FAIL"
    lines += [
        "", f"Published one-sided native-port validation: **{native}**.",
        f"Reflection and pathway validation of the symmetric extension: **{symmetric}**.",
        f"Overall production audit: **{overall}**.", "",
    ]
    (OUTPUT / "DEBELLY_2026_COMMON_AUDIT.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
