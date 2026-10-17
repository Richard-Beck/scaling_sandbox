#!/usr/bin/env python3
"""Create standalone diagnostic figures for the De Belly 1D WGD pair."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "wgd_v1" / "debelly_homeostatic"
RUNS = {
    "1x": RESULTS / "debelly_one_d_1x",
    "WGD": RESULTS / "debelly_one_d_2x",
}
FINE_RUNS = {
    "1x": RESULTS / "convergence" / "debelly_one_d_1x_fine",
    "WGD": RESULTS / "convergence" / "debelly_one_d_2x_fine",
}
COLORS = {"1x": "#0072B2", "WGD": "#D55E00"}
PULSE = (50.0, 150.0)


def load_table(path: Path, fields: tuple[str, ...]):
    raw = np.genfromtxt(path, delimiter=",", names=True, dtype=float)
    times = np.unique(raw["time_s"])
    x = np.unique(raw["position_um"])
    arrays = {
        field: np.asarray(raw[field]).reshape(times.size, x.size)
        for field in fields
    }
    return times, x, arrays


def load_run(path: Path):
    times, x, chemistry = load_table(
        path / "trajectory.csv", ("rac", "rho", "mca_bound", "mca_unbound")
    )
    _, tension_x, tension = load_table(path / "tension.csv", ("tension",))
    summary = json.loads((path / "summary.json").read_text(encoding="utf-8"))
    return {
        "times": times, "x": x, **chemistry,
        "tension_x": tension_x, "tension": tension["tension"],
        "summary": summary,
    }


def save(fig, directory: Path, stem: str):
    fig.savefig(directory / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(directory / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def heatmaps(runs, directory):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    for row, field in enumerate(("rac", "rho")):
        vmax = max(np.max(np.abs(run[field] - run[field][0])) for run in runs.values())
        for column, (label, run) in enumerate(runs.items()):
            change = run[field] - run[field][0]
            image = axes[row, column].imshow(
                change, origin="lower", aspect="auto", cmap="RdBu_r",
                vmin=-vmax, vmax=vmax,
                extent=(0, 1, run["times"][0], run["times"][-1]),
            )
            axes[row, column].axhline(PULSE[0], color="black", lw=0.7, ls=":")
            axes[row, column].axhline(PULSE[1], color="black", lw=0.7, ls=":")
            axes[row, column].set_title(f"{label}: {field.capitalize()} change")
            fig.colorbar(image, ax=axes[row, column], shrink=0.82)
    axes[0, 0].set_ylabel("Time (s)")
    axes[1, 0].set_ylabel("Time (s)")
    for axis in axes[-1]:
        axis.set_xlabel("Relative position (stimulated → distal)")
    fig.suptitle("De Belly 1D WGD: propagation of Rac/Rho responses")
    fig.tight_layout()
    save(fig, directory, "01_rac_rho_response_heatmaps")


def communication_timeseries(runs, directory):
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    specifications = (
        ("rac", "x", "Rac change"),
        ("rho", "x", "Rho change"),
        ("tension", "tension_x", "Tension change"),
    )
    for axis, (field, x_key, ylabel) in zip(axes, specifications):
        for label, run in runs.items():
            relative = run[x_key] / run[x_key][-1]
            proximal = np.mean(run[field][:, relative <= 0.2], axis=1)
            distal = np.mean(run[field][:, relative >= 0.8], axis=1)
            axis.plot(run["times"], proximal - proximal[0], color=COLORS[label],
                      lw=2, label=f"{label} proximal")
            axis.plot(run["times"], distal - distal[0], color=COLORS[label],
                      lw=2, ls="--", label=f"{label} distal")
        axis.axvspan(*PULSE, color="#BBBBBB", alpha=0.2)
        axis.axhline(0, color="black", lw=0.6)
        axis.set_ylabel(ylabel)
        axis.legend(ncol=2, fontsize=8, frameon=False)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Matched proximal and distal communication responses")
    fig.tight_layout()
    save(fig, directory, "02_proximal_distal_timecourses")


def selected_profiles(runs, directory):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    selected = (50.0, 150.0, 250.0, 350.0)
    time_colors = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, len(selected)))
    specifications = (
        ("rac", "x", "Rac"), ("rho", "x", "Rho"),
        ("mca_bound", "x", "Bound MCA"),
        ("tension", "tension_x", "Tension"),
    )
    for axis, (field, x_key, title) in zip(axes.flat, specifications):
        for label, run in runs.items():
            relative = run[x_key] / run[x_key][-1]
            for time, color in zip(selected, time_colors):
                index = int(np.argmin(np.abs(run["times"] - time)))
                axis.plot(relative, run[field][index], color=color,
                          ls="-" if label == "1x" else "--", lw=1.6)
        axis.set_title(title)
        axis.set_ylabel("Native model units")
    for axis in axes[-1]:
        axis.set_xlabel("Relative position (stimulated → distal)")
    time_handles = [plt.Line2D([], [], color=color, label=f"{time:g} s")
                    for time, color in zip(selected, time_colors)]
    size_handles = [
        plt.Line2D([], [], color="black", ls="-", label="1x"),
        plt.Line2D([], [], color="black", ls="--", label="WGD"),
    ]
    fig.legend(handles=time_handles + size_handles, loc="lower center", ncol=6,
               frameon=False)
    fig.suptitle("Spatial profiles at matched physical times")
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, directory, "03_selected_spatial_profiles")


def metric_comparison(runs, directory):
    summaries = {label: run["summary"]["model_metrics"] for label, run in runs.items()}
    shared = sorted(set(summaries["1x"]) & set(summaries["WGD"]))
    rows = []
    for metric in shared:
        baseline = summaries["1x"][metric]
        wgd = summaries["WGD"][metric]
        if isinstance(baseline, (int, float)) and isinstance(wgd, (int, float)):
            rows.append({
                "metric": metric, "one_x": baseline, "wgd": wgd,
                "wgd_to_one_x_ratio": wgd / baseline if baseline != 0 else np.nan,
                "wgd_minus_one_x": wgd - baseline,
            })
    with (directory / "debelly_wgd_metric_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    selected = (
        ("spatial_contrast", "Terminal Rac\ncontrast"),
        ("max_spatial_contrast", "Maximum Rac\ncontrast"),
        ("rac_distal_peak_time_s", "Distal Rac\npeak time"),
        ("rho_distal_peak_time_s", "Distal Rho\npeak time"),
        ("tension_distal_peak_time_s", "Distal tension\npeak time"),
        ("maximum_absolute_displacement", "Maximum absolute\ndisplacement"),
    )
    ratios = [summaries["WGD"][key] / summaries["1x"][key] for key, _ in selected]
    fig, axis = plt.subplots(figsize=(10, 4.8))
    bars = axis.bar(range(len(selected)), ratios, color="#D55E00")
    axis.axhline(1, color="black", lw=1)
    axis.set_xticks(range(len(selected)), [label for _, label in selected])
    axis.set_ylabel("WGD / 1x")
    axis.set_title("Effect sizes under homeostatic MCA density")
    for bar, ratio in zip(bars, ratios):
        axis.text(bar.get_x() + bar.get_width()/2, ratio, f"{ratio:.2f}",
                  ha="center", va="bottom", fontsize=9)
    axis.set_ylim(0, max(ratios) * 1.16)
    fig.tight_layout()
    save(fig, directory, "04_wgd_effect_ratios")


def convergence_comparison(runs, fine_runs, directory):
    metrics = {label: run["summary"]["model_metrics"] for label, run in runs.items()}
    fine = {label: run["summary"]["model_metrics"] for label, run in fine_runs.items()}
    selected = (
        ("spatial_contrast", "Terminal Rac\ncontrast"),
        ("max_spatial_contrast", "Maximum Rac\ncontrast"),
        ("rac_distal_peak_delta", "Distal Rac\npeak magnitude"),
        ("rho_distal_peak_delta", "Distal Rho\npeak magnitude"),
        ("tension_distal_peak_delta", "Distal tension\npeak magnitude"),
        ("maximum_absolute_displacement", "Maximum absolute\ndisplacement"),
    )
    rows = []
    for key, _ in selected:
        rows.append({
            "metric": key,
            "primary_1x": metrics["1x"][key],
            "primary_wgd": metrics["WGD"][key],
            "primary_wgd_to_1x": metrics["WGD"][key] / metrics["1x"][key],
            "fine_1x": fine["1x"][key],
            "fine_wgd": fine["WGD"][key],
            "fine_wgd_to_1x": fine["WGD"][key] / fine["1x"][key],
        })
    with (directory / "debelly_wgd_convergence_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    x = np.arange(len(selected))
    width = 0.36
    fig, axis = plt.subplots(figsize=(11, 5))
    axis.bar(x - width/2, [row["primary_wgd_to_1x"] for row in rows], width,
             color="#D55E00", label="Primary mesh")
    axis.bar(x + width/2, [row["fine_wgd_to_1x"] for row in rows], width,
             color="#E69F00", label="Doubled resolution")
    axis.axhline(1, color="black", lw=1)
    axis.set_xticks(x, [label for _, label in selected])
    axis.set_ylabel("WGD / 1x")
    axis.set_title("WGD effect ratios are stable after doubling spatial resolution")
    axis.legend(frameon=False)
    axis.set_ylim(0.9, 1.24)
    fig.tight_layout()
    save(fig, directory, "05_spatial_convergence")


def main():
    figures = RESULTS / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    runs = {label: load_run(path) for label, path in RUNS.items()}
    fine_runs = {label: load_run(path) for label, path in FINE_RUNS.items()}
    heatmaps(runs, figures)
    communication_timeseries(runs, figures)
    selected_profiles(runs, figures)
    metric_comparison(runs, figures)
    convergence_comparison(runs, fine_runs, figures)
    manifest = {
        "source_runs": {label: str(path.relative_to(ROOT)) for label, path in RUNS.items()},
        "convergence_runs": {
            label: str(path.relative_to(ROOT)) for label, path in FINE_RUNS.items()
        },
        "figures": sorted(path.name for path in figures.glob("*.png")),
        "note": "Standalone diagnostics; no report source was created or modified.",
    }
    (figures / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote De Belly WGD figures to {figures}")


if __name__ == "__main__":
    main()
