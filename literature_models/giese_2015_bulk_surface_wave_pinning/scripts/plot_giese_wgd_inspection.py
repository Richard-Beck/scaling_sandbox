#!/usr/bin/env python3
"""Generate inspection figures for the long Giese size and WGD experiments."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NATIVE_DIAMETERS = (1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 15.0)
WGD_LABELS = ("1x: d = 10 µm", "2x volume: d = 12.60 µm")
WGD_COLOURS = ("#0072B2", "#D55E00")


def style() -> None:
    plt.rcParams.update({
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "figure.dpi": 150,
        "savefig.dpi": 220,
    })


def native_path(native: Path, mode: str, diameter: float) -> Path:
    return native / "checkpoints" / f"{mode}_diameter_{diameter:04.1f}.npz"


def wgd_path(wgd: Path, index: int) -> Path:
    return wgd / "checkpoints" / f"two_site_wgd_{('1x', '2x')[index]}.npz"


def load(path: Path) -> dict[str, np.ndarray | dict]:
    with np.load(path, allow_pickle=False) as raw:
        result = {key: raw[key].copy() for key in raw.files}
    metadata_key = "metadata_json" if "metadata_json" in result else "summary_json"
    result["metadata"] = json.loads(str(result[metadata_key]))
    return result


def vertex_measure(edge_lengths: np.ndarray) -> np.ndarray:
    return 0.5 * (edge_lengths + np.roll(edge_lengths, 1))


def angular_distance(angle: np.ndarray, center: float) -> np.ndarray:
    return np.angle(np.exp(1j * (angle - center)))


def competition_metrics(data: dict) -> dict[str, np.ndarray]:
    angle = data["boundary_angles"]
    membrane = data["membrane"]
    weights = vertex_measure(data["edge_lengths"])
    strong_half = np.abs(angular_distance(angle, 0.0)) <= math.pi / 2.0
    weak_half = ~strong_half
    strong_height = membrane[:, strong_half].max(axis=1)
    weak_height = membrane[:, weak_half].max(axis=1)
    strong_mass = membrane[:, strong_half] @ weights[strong_half]
    weak_mass = membrane[:, weak_half] @ weights[weak_half]
    return {
        "strong_height": strong_height,
        "weak_height": weak_height,
        "height_bias": ((strong_height - weak_height)
                        / np.maximum(strong_height + weak_height, 1e-15)),
        "strong_mass": strong_mass,
        "weak_mass": weak_mass,
        "mass_bias": ((strong_mass - weak_mass)
                      / np.maximum(strong_mass + weak_mass, 1e-15)),
    }


def periodic_profile(data: dict, time: float, grid: np.ndarray) -> np.ndarray:
    index = int(np.argmin(np.abs(data["times"] - time)))
    angle = data["boundary_angles"]
    values = data["membrane"][index]
    order = np.argsort(angle)
    return np.interp(grid, angle[order], values[order], period=2.0 * math.pi)


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def protocol_figure(figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    for ax, radius, title in zip(
        axes, (1.0, 2.0 ** (1.0 / 3.0)), WGD_LABELS
    ):
        theta = np.linspace(0, 2 * math.pi, 500)
        ax.fill(radius * np.cos(theta), radius * np.sin(theta), "#DCEAF7")
        ax.plot(radius * np.cos(theta), radius * np.sin(theta), color="#333333", lw=2)
        for center, colour, label in (
            (0.0, "#009E73", "S1 = 0.44"),
            (math.pi, "#CC79A7", "S2 = 0.40"),
        ):
            arc = center + np.linspace(-0.05 * math.pi, 0.05 * math.pi, 50)
            ax.plot(radius * np.cos(arc), radius * np.sin(arc), color=colour,
                    lw=8, solid_capstyle="round", label=label)
        ax.annotate("active u on membrane", (0, -1.12 * radius), ha="center")
        ax.text(0, 0, "inactive v\nin 2D bulk", ha="center", va="center")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.set_xlim(-1.45 * radius, 1.45 * radius)
        ax.set_ylim(-1.35 * radius, 1.35 * radius)
        ax.axis("off")
    axes[0].legend(frameon=False, loc="upper center", ncol=2)
    fig.suptitle("Giese two-site protocol: geometry alone changes", fontsize=14)
    fig.text(0.5, 0.01,
             "Both patches cover 5% of the perimeter and act for 0–10 s; "
             "concentrations, diffusion, and reaction parameters are fixed.",
             ha="center")
    save(fig, figure_dir / "01_protocol_schematic.png")


def native_size_metrics(native: Path, figure_dir: Path) -> None:
    rows = []
    with (native / "native_size_summary.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append({key: (float(value) if key not in {
                "mode", "polarized_final", "unique_site_final",
                "winner_selection_right_censored"
            } else value) for key, value in row.items()})
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), constrained_layout=True)
    specs = (
        ("final_relative_amplitude", "Relative membrane amplitude at endpoint"),
        ("final_POL", "POL at endpoint"),
        ("final_half_mass_PF", "Half-mass PF at endpoint"),
    )
    colours = {"single": "#0072B2", "two_site": "#D55E00"}
    for ax, (key, ylabel) in zip(axes.flat[:3], specs):
        for mode, label in (("single", "One stimulus"),
                            ("two_site", "Two opposing stimuli")):
            chosen = [row for row in rows if row["mode"] == mode]
            ax.plot([row["diameter_um"] for row in chosen],
                    [row[key] for row in chosen], "o-", lw=2,
                    color=colours[mode], label=label)
        ax.set(xlabel="Diameter (µm)", ylabel=ylabel, xticks=NATIVE_DIAMETERS)
    axes[0, 0].legend(frameon=False)
    ax = axes[1, 1]
    two = [row for row in rows if row["mode"] == "two_site"]
    x = np.asarray([row["diameter_um"] for row in two])
    y = np.asarray([row["winner_selection_time_s"] for row in two])
    finite = np.isfinite(y)
    ax.plot(x[finite], y[finite], "o-", color=colours["two_site"], lw=2,
            label="Current solver")
    censored = (~finite) & (x >= 3.0)
    ax.scatter(x[censored], np.full(censored.sum(), 2000), marker="^", s=75,
               facecolors="none", edgecolors=colours["two_site"],
               label="No winner by 2000 s")
    ax.plot((3.0, 6.0), (43.0, 953.0), "s--", color="#444444",
            label="Paper reference")
    ax.set(xlabel="Diameter (µm)", ylabel="Winner-selection time (s)",
           xticks=NATIVE_DIAMETERS, yscale="log")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Native Giese diameter sweep", fontsize=14)
    save(fig, figure_dir / "02_native_size_sweep.png")


def native_competition_trajectories(native: Path, figure_dir: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.0), sharex=True,
                             constrained_layout=True)
    colours = plt.cm.viridis(np.linspace(0.05, 0.95, len(NATIVE_DIAMETERS)))
    for diameter, colour in zip(NATIVE_DIAMETERS, colours):
        data = load(native_path(native, "two_site", diameter))
        comp = competition_metrics(data)
        axes[0].plot(data["times"], comp["height_bias"], color=colour,
                     label=f"{diameter:g} µm")
        axes[1].plot(data["times"], comp["mass_bias"], color=colour)
    for ax in axes:
        ax.axvspan(0, 10, color="#BBBBBB", alpha=0.25)
        ax.axhline(0, color="#555555", lw=0.8)
    axes[0].set_ylabel("Peak-height bias\n(S1 − S2)/(S1 + S2)")
    axes[1].set_ylabel("Half-membrane mass bias\n(S1 − S2)/(S1 + S2)")
    axes[1].set_xlabel("Time (s)")
    axes[0].legend(frameon=False, ncol=4, fontsize=8)
    fig.suptitle("Competitive divergence slows strongly with diameter", fontsize=14)
    save(fig, figure_dir / "03_native_competition_bias.png")


def native_spacetime(native: Path, figure_dir: Path) -> None:
    chosen = (3.0, 4.5, 6.0, 9.0, 15.0)
    grid = np.linspace(0, 2 * math.pi, 256, endpoint=False)
    fig, axes = plt.subplots(len(chosen), 1, figsize=(10.5, 10), sharex=True,
                             constrained_layout=True)
    for ax, diameter in zip(axes, chosen):
        data = load(native_path(native, "two_site", diameter))
        matrix = np.asarray([periodic_profile(data, time, grid)
                             for time in data["times"]])
        relative = matrix / np.maximum(matrix.mean(axis=1, keepdims=True), 1e-15)
        artist = ax.imshow(relative.T, origin="lower", aspect="auto",
                           extent=(0, data["times"][-1], 0, 2 * math.pi),
                           cmap="magma", vmin=0.25, vmax=2.5)
        ax.set_ylabel(f"{diameter:g} µm\nangle")
        ax.set_yticks((0, math.pi, 2 * math.pi), ("0", "π", "2π"))
    axes[-1].set_xlabel("Time (s)")
    fig.colorbar(artist, ax=axes, shrink=0.72, label="Membrane u / perimeter mean")
    fig.suptitle("Two-site membrane dynamics across cell sizes", fontsize=14)
    save(fig, figure_dir / "04_native_two_site_spacetime.png")


def wgd_dashboard(wgd_data: list[dict], figure_dir: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(11, 9.2), sharex=True,
                             constrained_layout=True)
    for data, label, colour in zip(wgd_data, WGD_LABELS, WGD_COLOURS):
        comp = competition_metrics(data)
        t = data["times"]
        axes[0, 0].plot(t, comp["strong_height"], color=colour, label=label)
        axes[0, 0].plot(t, comp["weak_height"], color=colour, ls="--")
        axes[0, 1].plot(t, comp["height_bias"], color=colour, label=label)
        axes[1, 0].plot(t, comp["mass_bias"], color=colour)
        axes[1, 1].plot(t, data["metric_relative_amplitude"], color=colour)
        axes[2, 0].plot(t, data["metric_PF"], color=colour)
        axes[2, 1].step(t, data["metric_n_above_mean_clusters"], where="post",
                        color=colour)
    for ax in axes.flat:
        ax.axvspan(0, 10, color="#BBBBBB", alpha=0.25)
    axes[0, 0].set_ylabel("Local peak height")
    axes[0, 0].legend(frameon=False)
    axes[0, 0].text(0.98, 0.93, "solid: S1; dashed: S2",
                    transform=axes[0, 0].transAxes, ha="right", va="top")
    axes[0, 1].set_ylabel("Peak-height bias")
    axes[1, 0].set_ylabel("Half-membrane mass bias")
    axes[1, 1].set_ylabel("Relative amplitude")
    axes[2, 0].set_ylabel("PF")
    axes[2, 1].set_ylabel("Above-mean domains")
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    fig.suptitle("Long two-site WGD comparison", fontsize=14)
    save(fig, figure_dir / "05_wgd_long_run_dashboard.png")


def wgd_spacetime(wgd_data: list[dict], figure_dir: Path) -> None:
    grid = np.linspace(0, 2 * math.pi, 320, endpoint=False)
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.4), sharex=True,
                             constrained_layout=True)
    for ax, data, label in zip(axes, wgd_data, WGD_LABELS):
        matrix = np.asarray([periodic_profile(data, time, grid)
                             for time in data["times"]])
        relative = matrix / np.maximum(matrix.mean(axis=1, keepdims=True), 1e-15)
        artist = ax.imshow(relative.T, origin="lower", aspect="auto",
                           extent=(0, data["times"][-1], 0, 2 * math.pi),
                           cmap="magma", vmin=0.25, vmax=2.3)
        ax.set_title(label)
        ax.set_ylabel("Membrane angle")
        ax.set_yticks((0, math.pi, 2 * math.pi), ("S1: 0", "S2: π", "2π"))
    axes[-1].set_xlabel("Time (s)")
    fig.colorbar(artist, ax=axes, shrink=0.78, label="Membrane u / perimeter mean")
    fig.suptitle("WGD pair: both opposing domains persist through 2000 s", fontsize=14)
    save(fig, figure_dir / "06_wgd_spacetime.png")


def wgd_profiles(wgd_data: list[dict], figure_dir: Path) -> None:
    times = (10.0, 100.0, 300.0, 1000.0, 2000.0)
    grid = np.linspace(0, 2 * math.pi, 500)
    fig, axes = plt.subplots(2, len(times), figsize=(15, 6.2), sharex=True,
                             constrained_layout=True)
    for row, (data, label, colour) in enumerate(zip(wgd_data, WGD_LABELS, WGD_COLOURS)):
        for ax, time in zip(axes[row], times):
            profile = periodic_profile(data, time, grid)
            ax.plot(grid, profile, color=colour, lw=2)
            ax.axvline(0, color="#009E73", alpha=0.5, lw=1)
            ax.axvline(math.pi, color="#CC79A7", alpha=0.5, lw=1)
            ax.set_title(f"t = {time:g} s")
            ax.set_xticks((0, math.pi, 2 * math.pi), ("0", "π", "2π"))
            if ax is axes[row, 0]:
                ax.set_ylabel(label + "\nmembrane u")
    fig.suptitle("Membrane profiles during two-site competition", fontsize=14)
    save(fig, figure_dir / "07_wgd_selected_profiles.png")


def wgd_bulk_fields(wgd_data: list[dict], figure_dir: Path) -> None:
    desired = (10.0, 300.0, 2000.0)
    fig, axes = plt.subplots(2, len(desired), figsize=(11.2, 7.0),
                             constrained_layout=True)
    for row, (data, label) in enumerate(zip(wgd_data, WGD_LABELS)):
        points = data["mesh_points"]
        triangulation = mtri.Triangulation(
            points[0], points[1], data["mesh_triangles"].T
        )
        indices = [int(np.argmin(np.abs(data["selected_times"] - time)))
                   for time in desired]
        fields = data["selected_bulk"]
        vmin = min(float(fields[index].min()) for index in indices)
        vmax = max(float(fields[index].max()) for index in indices)
        for column, (time, index) in enumerate(zip(desired, indices)):
            ax = axes[row, column]
            artist = ax.tripcolor(triangulation, fields[index], shading="gouraud",
                                  cmap="viridis", vmin=vmin, vmax=vmax)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_title(f"t = {time:g} s")
            if column == 0:
                ax.text(-0.10, 0.5, label, rotation=90, transform=ax.transAxes,
                        ha="center", va="center")
        fig.colorbar(artist, ax=axes[row], shrink=0.70,
                     label="Cytoplasmic v (µM)")
    fig.suptitle("Selected cytoplasmic fields in the WGD pair", fontsize=14)
    save(fig, figure_dir / "08_wgd_bulk_fields.png")


def endpoint_comparison(wgd_data: list[dict], figure_dir: Path) -> None:
    grid = np.linspace(0, 2 * math.pi, 600)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), constrained_layout=True)
    for data, label, colour in zip(wgd_data, WGD_LABELS, WGD_COLOURS):
        profile = periodic_profile(data, 2000.0, grid)
        axes[0].plot(grid, profile, color=colour, lw=2, label=label)
        axes[1].plot(grid, profile / profile.mean(), color=colour, lw=2)
        comp = competition_metrics(data)
        axes[2].plot(data["times"],
                     comp["strong_height"] / np.maximum(comp["weak_height"], 1e-15),
                     color=colour, lw=2)
    axes[0].set(xlabel="Membrane angle", ylabel="u at 2000 s",
                xticks=(0, math.pi, 2 * math.pi))
    axes[0].set_xticklabels(("0", "π", "2π"))
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].set(xlabel="Membrane angle", ylabel="u / perimeter mean",
                xticks=(0, math.pi, 2 * math.pi))
    axes[1].set_xticklabels(("0", "π", "2π"))
    axes[2].set(xlabel="Time (s)", ylabel="S1 peak / S2 peak")
    axes[2].axhline(1, color="#555555", lw=0.8)
    fig.suptitle("Endpoint pattern and slow competitive asymmetry", fontsize=14)
    save(fig, figure_dir / "09_wgd_endpoint_comparison.png")


def numerical_quality(wgd_data: list[dict], figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8), constrained_layout=True)
    for data, label, colour in zip(wgd_data, WGD_LABELS, WGD_COLOURS):
        relative = data["total_mass"] / data["total_mass"][0] - 1.0
        axes[0].plot(data["times"], relative, color=colour, label=label)
    axes[0].set(xlabel="Time (s)", ylabel="Relative native-mass error")
    axes[0].legend(frameon=False, fontsize=8)
    labels = [label.split(":")[0] for label in WGD_LABELS]
    meta = [data["metadata"] for data in wgd_data]
    axes[1].bar(labels, [row["actual_max_h_um"] for row in meta],
                color=WGD_COLOURS, width=0.55)
    axes[1].axhline(0.25, color="#555555", ls="--", label="target h")
    axes[1].set(ylabel="Maximum mesh edge (µm)")
    axes[1].legend(frameon=False)
    fig.suptitle("Numerical quality of the long WGD runs", fontsize=14)
    save(fig, figure_dir / "10_wgd_numerical_quality.png")


def transition_region(native: Path, wgd: Path, figure_dir: Path) -> list[dict]:
    pairs = (
        (4.5, "4p5_to_2x"),
        (6.0, "6_to_2x"),
    )
    rows = []
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.2), sharex=True,
                             constrained_layout=True)
    for row_index, (baseline, enlarged_label) in enumerate(pairs):
        baseline_data = load(native_path(native, "two_site", baseline))
        enlarged_data = load(
            wgd / "checkpoints" / f"two_site_wgd_{enlarged_label}.npz"
        )
        for data, scale_label, colour in (
            (baseline_data, f"1x: {baseline:g} µm", WGD_COLOURS[0]),
            (enlarged_data,
             f"2x volume: {baseline * 2 ** (1 / 3):.2f} µm", WGD_COLOURS[1]),
        ):
            comp = competition_metrics(data)
            time = data["times"]
            axes[row_index, 0].plot(time, comp["mass_bias"], color=colour,
                                    label=scale_label)
            axes[row_index, 1].plot(time, data["metric_relative_amplitude"],
                                    color=colour)
            axes[row_index, 2].step(time, data["metric_n_above_mean_clusters"],
                                    where="post", color=colour)
            rows.append({
                "pair_baseline_um": baseline,
                "scale": scale_label.split(":")[0],
                "diameter_um": float(data["metadata"]["diameter_um"]),
                "domains_at_300_s": int(data["metric_n_above_mean_clusters"][300]),
                "domains_at_2000_s": int(data["metric_n_above_mean_clusters"][-1]),
                "mass_bias_at_300_s": float(comp["mass_bias"][300]),
                "mass_bias_at_2000_s": float(comp["mass_bias"][-1]),
                "relative_amplitude_at_2000_s": float(
                    data["metric_relative_amplitude"][-1]
                ),
            })
        axes[row_index, 0].set_ylabel(
            f"Baseline {baseline:g} µm\nS1–S2 mass bias"
        )
        axes[row_index, 0].legend(frameon=False, fontsize=8)
        axes[row_index, 1].set_ylabel("Relative amplitude")
        axes[row_index, 2].set_ylabel("Above-mean domains")
        for ax in axes[row_index]:
            ax.axvspan(0, 10, color="#BBBBBB", alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    fig.suptitle("Transition-region WGD pairs", fontsize=14)
    save(fig, figure_dir / "11_transition_region_wgd_pairs.png")
    path = wgd / "transition_region_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_summary(wgd_data: list[dict], output: Path) -> list[dict]:
    rows = []
    for data, label in zip(wgd_data, ("1x", "2x")):
        comp = competition_metrics(data)
        row = {
            "scale": label,
            "diameter_um": data["metadata"]["diameter_um"],
            "volume_scale": data["metadata"]["volume_scale"],
            "initial_native_mass": data["metadata"]["initial_native_mass"],
            "domains_at_300_s": int(data["metric_n_above_mean_clusters"][300]),
            "domains_at_2000_s": int(data["metric_n_above_mean_clusters"][-1]),
            "relative_amplitude_at_300_s": float(data["metric_relative_amplitude"][300]),
            "relative_amplitude_at_2000_s": float(data["metric_relative_amplitude"][-1]),
            "peak_height_bias_at_300_s": float(comp["height_bias"][300]),
            "peak_height_bias_at_2000_s": float(comp["height_bias"][-1]),
            "half_mass_bias_at_300_s": float(comp["mass_bias"][300]),
            "half_mass_bias_at_2000_s": float(comp["mass_bias"][-1]),
            "final_PF": float(data["metric_PF"][-1]),
            "maximum_relative_mass_error": data["metadata"]["maximum_relative_mass_error"],
        }
        rows.append(row)
    path = output / "wgd_long_run_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--native", type=Path,
        default=ROOT / "reports/output/giese_wp/native",
    )
    parser.add_argument(
        "--wgd", type=Path,
        default=ROOT / "reports/output/giese_wp/wgd_inspection",
    )
    args = parser.parse_args()
    figure_dir = args.wgd / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    style()
    wgd_data = [load(wgd_path(args.wgd, index)) for index in range(2)]
    protocol_figure(figure_dir)
    native_size_metrics(args.native, figure_dir)
    native_competition_trajectories(args.native, figure_dir)
    native_spacetime(args.native, figure_dir)
    wgd_dashboard(wgd_data, figure_dir)
    wgd_spacetime(wgd_data, figure_dir)
    wgd_profiles(wgd_data, figure_dir)
    wgd_bulk_fields(wgd_data, figure_dir)
    endpoint_comparison(wgd_data, figure_dir)
    numerical_quality(wgd_data, figure_dir)
    rows = write_summary(wgd_data, args.wgd)
    transition_rows = transition_region(args.native, args.wgd, figure_dir)
    manifest = {
        "figures": sorted(path.name for path in figure_dir.glob("*.png")),
        "native_source": str(args.native),
        "wgd_source": str(args.wgd),
        "summary": rows,
        "transition_region_summary": transition_rows,
    }
    (args.wgd / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
