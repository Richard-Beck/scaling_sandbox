#!/usr/bin/env python3
"""Publish report-ready figures for the conditional Giese WGD crossover."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "reports/assets/wgd_v1/giese_crossover"
TUNED_DIR = ROOT / "reports/output/giese_wp/boundary_confirmation/checkpoints"
DEFAULT_DIR = ROOT / "reports/output/giese_wp/diffusion_screen/checkpoints"
COMPETITION_DIR = ROOT / "reports/output/giese_wp/competition_screen/checkpoints"
PILOT_CSV = ROOT / "reports/output/giese_wp/boundary_pilot/boundary_pilot_summary.csv"

COLOURS = {"1x": "#0072B2", "WGD": "#D55E00"}
RUNS = {
    "crossover": {"1x": TUNED_DIR / "task_012.npz",
                  "WGD": TUNED_DIR / "task_013.npz"},
    "default": {"1x": DEFAULT_DIR / "task_004.npz",
                "WGD": DEFAULT_DIR / "task_005.npz"},
    "matched_control": {"1x": COMPETITION_DIR / "task_070.npz",
                        "WGD": COMPETITION_DIR / "task_071.npz"},
}


def style() -> None:
    plt.rcParams.update({
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "figure.dpi": 150,
        "savefig.dpi": 240,
    })


def load(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as raw:
        result = {key: raw[key] for key in raw.files}
    result["summary"] = json.loads(str(result["summary_json"]))
    order = np.argsort(result["boundary_angles"])
    result["angles"] = result["boundary_angles"][order]
    result["membrane_sorted"] = result["membrane"][:, order]
    if "metric_domains" not in result:
        result["metric_domains"] = result["metric_n_above_mean_clusters"]
    return result


def nearest_time(data: dict, time: float) -> int:
    return int(np.argmin(np.abs(data["times"] - time)))


def membrane_segments(angles: np.ndarray, radius: float) -> np.ndarray:
    extended = np.r_[angles, angles[0] + 2 * math.pi]
    xy = np.column_stack((radius * np.cos(extended), radius * np.sin(extended)))
    return np.stack((xy[:-1], xy[1:]), axis=1)


def draw_cell(ax: plt.Axes, data: dict, time: float, norm: Normalize,
              cmap, limit: float, show_cues: bool = True) -> None:
    index = nearest_time(data, time)
    radius = data["summary"]["diameter_um"] / 2.0
    values = data["membrane_sorted"][index]
    angles = data["angles"]
    segments = membrane_segments(angles, radius)
    segment_values = 0.5 * (values + np.roll(values, -1))
    ax.add_patch(plt.Circle((0, 0), radius, facecolor="#F4F6F8",
                            edgecolor="#C4CBD2", linewidth=1.0))
    ax.add_collection(LineCollection(segments, colors="white", linewidth=8.0,
                                     capstyle="round"))
    ax.add_collection(LineCollection(segments, array=segment_values, cmap=cmap,
                                     norm=norm, linewidth=6.0, capstyle="round"))
    if show_cues:
        ax.scatter((radius + 0.55, -radius - 0.55), (0, 0), marker="<", s=30,
                   c=("#6A3D9A", "#8C8C8C"), edgecolors="none", zorder=5)
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect="equal")
    ax.axis("off")


def scenario_norm(runs: dict[str, dict]) -> Normalize:
    values = np.concatenate([run["membrane"].ravel() for run in runs.values()])
    return Normalize(vmin=float(np.quantile(values, 0.01)),
                     vmax=float(np.quantile(values, 0.995)))


def snapshot_figure(runs: dict[str, dict], scenario: str, times: tuple[int, ...],
                    filename: str, subtitle: str) -> None:
    cmap = plt.get_cmap("magma")
    norm = scenario_norm(runs)
    limit = max(run["summary"]["diameter_um"] for run in runs.values()) / 2 + 1.0
    fig, axes = plt.subplots(2, len(times), figsize=(2.05 * len(times), 4.7),
                             constrained_layout=True)
    for row, label in enumerate(("1x", "WGD")):
        for column, time in enumerate(times):
            draw_cell(axes[row, column], runs[label], time, norm, cmap, limit)
            if row == 0:
                axes[row, column].set_title(f"{time} s", fontsize=11)
            if column == 0:
                diameter = runs[label]["summary"]["diameter_um"]
                axes[row, column].text(-limit - 0.15, 0, f"{label}\n{diameter:.1f} µm",
                                       ha="right", va="center", fontweight="bold",
                                       color=COLOURS[label])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(sm, ax=axes, shrink=0.76, pad=0.01,
                 label="Active membrane GTPase, u")
    fig.suptitle(subtitle, fontsize=14, fontweight="bold")
    fig.savefig(ASSETS / filename, bbox_inches="tight")
    plt.close(fig)


def kinetics_figure(all_runs: dict[str, dict[str, dict]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2), sharex=True,
                             constrained_layout=True)
    for column, scenario in enumerate(("crossover", "default")):
        for label in ("1x", "WGD"):
            run = all_runs[scenario][label]
            times = run["times"]
            bias = np.abs(run["competition_mass_bias"])
            domains = run["metric_domains"]
            axes[0, column].plot(times, bias, color=COLOURS[label],
                                 linewidth=1.7, label=label)
            axes[1, column].step(times, domains, where="post", color=COLOURS[label],
                                 linewidth=1.7, label=label)
        title = ("Crossover condition\n0.94× abundance, diffusion 0.1/20"
                 if scenario == "crossover" else
                 "Default condition\n1.00× abundance, diffusion 0.015/3")
        axes[0, column].set_title(title, fontweight="bold")
        axes[0, column].set_yscale("symlog", linthresh=1e-3)
        axes[0, column].axvspan(0, 10, color="#6A3D9A", alpha=0.10)
        axes[1, column].axvspan(0, 10, color="#6A3D9A", alpha=0.10)
        axes[1, column].set_yticks((0, 1, 2))
        axes[1, column].set_xlabel("Time (s)")
        axes[0, column].legend(frameon=False)
    axes[0, 0].set_ylabel("|Strong − weak membrane-mass bias|")
    axes[1, 0].set_ylabel("Detected domains")
    fig.suptitle("WGD changes winner selection only near the competition crossover",
                 fontsize=14, fontweight="bold")
    fig.savefig(ASSETS / "03_competition_kinetics.png", bbox_inches="tight")
    plt.close(fig)


def kymograph_figure(all_runs: dict[str, dict[str, dict]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True, sharey=True,
                             constrained_layout=True)
    for row, scenario in enumerate(("crossover", "default")):
        runs = all_runs[scenario]
        norm = scenario_norm(runs)
        for column, label in enumerate(("1x", "WGD")):
            run = runs[label]
            stride = 4
            artist = axes[row, column].imshow(
                run["membrane_sorted"][::stride], origin="lower", aspect="auto",
                extent=(-180, 180, run["times"][0], run["times"][-1]),
                cmap="magma", norm=norm, interpolation="nearest",
            )
            axes[row, column].axvline(0, color="white", linewidth=0.7, alpha=0.55)
            axes[row, column].set_title(
                f'{"Crossover" if scenario == "crossover" else "Default"}; {label}',
                fontweight="bold")
            fig.colorbar(artist, ax=axes[row, column], shrink=0.82, pad=0.02,
                         label="u")
    for ax in axes[-1]:
        ax.set_xlabel("Membrane angle (degrees)")
        ax.set_xticks((-180, -90, 0, 90, 180))
    for ax in axes[:, 0]:
        ax.set_ylabel("Time (s)")
    fig.suptitle("Space–time view of competing membrane domains",
                 fontsize=14, fontweight="bold")
    fig.savefig(ASSETS / "04_membrane_kymographs.png", bbox_inches="tight")
    plt.close(fig)


def boundary_figure() -> None:
    with PILOT_CSV.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    lookup = {(float(row["gamma"]), float(row["mass_scale"]),
               round(float(row["diameter_um"]), 1)): row for row in rows}
    masses = (0.92, 0.94, 0.96, 0.98, 1.00)
    gammas = (0.75, 1.00, 1.25)
    values = np.zeros((len(gammas), len(masses)))
    text = np.empty(values.shape, dtype=object)
    for i, gamma in enumerate(gammas):
        for j, mass in enumerate(masses):
            one = lookup[(gamma, mass, 10.0)]
            wgd = lookup[(gamma, mass, 12.6)]
            t1 = float(one["winner_selection_time_s"])
            t2 = float(wgd["winner_selection_time_s"])
            p1 = float(one["final_relative_amplitude"]) >= 0.05
            p2 = float(wgd["final_relative_amplitude"]) >= 0.05
            if np.isfinite(t1) and not np.isfinite(t2) and p2:
                values[i, j], text[i, j] = 3, "WGD\nsplit"
            elif np.isfinite(t1) and np.isfinite(t2):
                values[i, j], text[i, j] = 2, "both\nselect"
            elif p1 and p2:
                values[i, j], text[i, j] = 1, "both\nslow"
            else:
                values[i, j], text[i, j] = 0, "no stable\npolarity"
    cmap = ListedColormap(("#BDBDBD", "#E69F00", "#009E73", "#CC79A7"))
    norm = BoundaryNorm((-0.5, 0.5, 1.5, 2.5, 3.5), cmap.N)
    fig, ax = plt.subplots(figsize=(8.4, 4.4), constrained_layout=True)
    artist = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(masses)), masses)
    ax.set_yticks(range(len(gammas)), gammas)
    ax.set_xlabel("Total GTPase abundance relative to default")
    ax.set_ylabel(r"Positive-feedback strength, $\gamma$")
    ax.grid(False)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, text[i, j], ha="center", va="center", fontsize=9,
                    fontweight="bold" if values[i, j] == 3 else "normal")
    colourbar = fig.colorbar(artist, ax=ax, ticks=(0, 1, 2, 3), shrink=0.82)
    colourbar.ax.set_yticklabels(("No stable polarity", "Both slow",
                                  "Both select", "WGD split"))
    ax.set_title("The WGD effect occupies a narrow crossover regime",
                 fontweight="bold", fontsize=14)
    fig.savefig(ASSETS / "05_parameter_crossover.png", bbox_inches="tight")
    plt.close(fig)


def publish_summary(all_runs: dict[str, dict[str, dict]]) -> None:
    rows = []
    for scenario, runs in all_runs.items():
        for label, run in runs.items():
            summary = run["summary"]
            rows.append({
                "scenario": scenario,
                "geometry": label,
                "diameter_um": summary["diameter_um"],
                "abundance_scale": 0.94 if scenario == "crossover" else 1.0,
                "gamma": 1.0,
                "Dm_um2_s": 0.015 if scenario == "default" else 0.1,
                "Dc_um2_s": 3.0 if scenario == "default" else 20.0,
                "strong_stimulus_um_s": 0.44,
                "weak_stimulus_um_s": 0.40 if scenario == "default" else 0.30,
                "winner_selection_time_s": summary["winner_selection_time_s"],
                "right_censored_at_2000_s": summary[
                    "winner_selection_right_censored"],
                "domains_at_2000_s": summary["final_domains"],
                "relative_amplitude_at_2000_s": summary["final_relative_amplitude"],
                "mass_bias_at_2000_s": summary["mass_bias_final"],
                "maximum_relative_mass_error": summary["maximum_relative_mass_error"],
            })
    with (ASSETS / "giese_crossover_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    style()
    all_runs = {
        scenario: {label: load(path) for label, path in paths.items()}
        for scenario, paths in RUNS.items()
    }
    snapshot_figure(
        all_runs["crossover"], "crossover",
        (0, 10, 100, 300, 500, 1000, 2000),
        "01_crossover_cell_snapshots.png",
        "Crossover condition: WGD prevents resolution of two competing sites",
    )
    snapshot_figure(
        all_runs["default"], "default",
        (0, 10, 100, 300, 1000, 2000),
        "02_default_cell_snapshots.png",
        "Default condition: both sizes retain two domains",
    )
    snapshot_figure(
        all_runs["matched_control"], "matched_control",
        (0, 10, 100, 300, 500, 1000, 2000),
        "06_matched_abundance_control.png",
        "Matched control at 1.00× abundance: both sizes retain two domains",
    )
    kinetics_figure(all_runs)
    kymograph_figure(all_runs)
    boundary_figure()
    publish_summary(all_runs)
    print(ASSETS)


if __name__ == "__main__":
    main()
