#!/usr/bin/env python3
"""Plot the Giese Dm x Dc x diameter competition screen."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DM_VALUES = (0.015, 0.03, 0.05, 0.075, 0.10)
DC_VALUES = (3.0, 5.5, 10.0)
DIAMETERS = (3.0, 4.5, 6.0, 8.0, 10.0, 10.0 * 2.0 ** (1.0 / 3.0), 15.0)
END_TIME = 2000.0


def diameter_label(value: float) -> str:
    return "12.6" if math.isclose(value, 10.0 * 2.0 ** (1.0 / 3.0)) else f"{value:g}"


def style() -> None:
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.20,
        "figure.dpi": 150,
        "savefig.dpi": 220,
    })


def read_rows(output: Path) -> list[dict]:
    with (output / "diffusion_screen_summary.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    booleans = {"winner_selection_right_censored"}
    for row in rows:
        for key, value in tuple(row.items()):
            if key in booleans:
                row[key] = value.lower() == "true"
            else:
                row[key] = float(value)
    return rows


def row_map(rows: list[dict]) -> dict[tuple[float, float, float], dict]:
    return {
        (row["Dm_um2_s"], row["Dc_um2_s"], row["diameter_um"]): row
        for row in rows
    }


def matrix(mapping: dict, dc: float, key: str,
           transform=lambda value: value) -> np.ndarray:
    return np.asarray([
        [transform(mapping[(dm, dc, diameter)][key]) for diameter in DIAMETERS]
        for dm in DM_VALUES
    ], dtype=float)


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def annotated_heatmap(ax: plt.Axes, values: np.ndarray, title: str,
                      cmap: str, vmin: float, vmax: float,
                      formatter, colourbar_label: str):
    artist = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap,
                       vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xticks(range(len(DIAMETERS)),
                  [diameter_label(diameter) for diameter in DIAMETERS], rotation=30)
    ax.set_yticks(range(len(DM_VALUES)), [f"{dm:g}" for dm in DM_VALUES])
    ax.set_xlabel("Diameter (µm)")
    ax.set_ylabel("$D_m$ (µm²/s)")
    ax.grid(False)
    midpoint = 0.5 * (vmin + vmax)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            value = values[i, j]
            colour = "white" if value < midpoint else "black"
            ax.text(j, i, formatter(value), ha="center", va="center",
                    fontsize=7.5, color=colour)
    return artist, colourbar_label


def winner_phase_maps(mapping: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    for ax, dc in zip(axes, DC_VALUES):
        values = np.asarray([
            [
                (-200.0 if (
                    not np.isfinite(mapping[(dm, dc, diameter)]["winner_selection_time_s"])
                    and not mapping[(dm, dc, diameter)]["winner_selection_right_censored"]
                ) else END_TIME if not np.isfinite(
                    mapping[(dm, dc, diameter)]["winner_selection_time_s"]
                ) else mapping[(dm, dc, diameter)]["winner_selection_time_s"])
                for diameter in DIAMETERS
            ]
            for dm in DM_VALUES
        ])
        artist, label = annotated_heatmap(
            ax, values, f"$D_c$ = {dc:g} µm²/s", "viridis", -200, END_TIME,
            lambda value: ("no pol." if value < 0 else
                           ">2000" if value >= END_TIME else f"{value:.0f}"),
            "Stable winner-selection time (s)",
        )
    fig.colorbar(artist, ax=axes, shrink=0.82, label=label)
    fig.suptitle("Two-site competition time across diffusion and size", fontsize=14)
    save(fig, figure_dir / "01_winner_time_phase_maps.png")


def domain_phase_maps(mapping: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    for ax, dc in zip(axes, DC_VALUES):
        values = matrix(mapping, dc, "final_domains")
        artist, label = annotated_heatmap(
            ax, values, f"$D_c$ = {dc:g} µm²/s", "cividis", 0, 2,
            lambda value: f"{value:.0f}", "Domains at 2000 s",
        )
    fig.colorbar(artist, ax=axes, shrink=0.82, label=label, ticks=(0, 1, 2))
    fig.suptitle("Final domain count", fontsize=14)
    save(fig, figure_dir / "02_final_domain_phase_maps.png")


def bias_phase_maps(mapping: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    for ax, dc in zip(axes, DC_VALUES):
        values = matrix(
            mapping, dc, "mass_bias_final",
            lambda value: math.log10(max(abs(value), 1e-6)),
        )
        artist, label = annotated_heatmap(
            ax, values, f"$D_c$ = {dc:g} µm²/s", "magma", -6, 0,
            lambda value: f"{value:.1f}", "log10 |S1–S2 membrane-mass bias|",
        )
    fig.colorbar(artist, ax=axes, shrink=0.82, label=label)
    fig.suptitle("Residual competitive asymmetry at 2000 s", fontsize=14)
    save(fig, figure_dir / "03_endpoint_bias_phase_maps.png")


def crossover_curves(mapping: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), sharey=True,
                             constrained_layout=True)
    colours = plt.cm.viridis(np.linspace(0.05, 0.95, len(DM_VALUES)))
    for ax, dc in zip(axes, DC_VALUES):
        for dm, colour in zip(DM_VALUES, colours):
            chosen = [mapping[(dm, dc, diameter)] for diameter in DIAMETERS]
            times = np.asarray([row["winner_selection_time_s"] for row in chosen])
            censored = np.asarray([
                row["winner_selection_right_censored"] for row in chosen
            ])
            no_polarity = ~np.isfinite(times) & ~censored
            plotted = np.where(np.isfinite(times), times,
                               np.where(censored, END_TIME, np.nan))
            ax.plot(DIAMETERS, plotted, "o-", color=colour, label=f"{dm:g}")
            ax.scatter(np.asarray(DIAMETERS)[censored], plotted[censored],
                       marker="^", facecolors="none", edgecolors=colour, s=48)
            ax.scatter(np.asarray(DIAMETERS)[no_polarity],
                       np.full(no_polarity.sum(), 25.0), marker="x", color=colour,
                       s=42)
        ax.axvspan(10, 10 * 2 ** (1 / 3), color="#D55E00", alpha=0.10)
        ax.set(title=f"$D_c$ = {dc:g}", xlabel="Diameter (µm)", yscale="log")
    axes[0].set_ylabel("Winner-selection time (s)")
    axes[-1].legend(title="$D_m$", frameon=False, fontsize=8, ncol=2)
    fig.suptitle("Competition crossover curves; triangles are censored, crosses are nonpolarized",
                 fontsize=14)
    save(fig, figure_dir / "04_crossover_curves.png")


def endpoint_metric_maps(mapping: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.2), constrained_layout=True)
    for column, dc in enumerate(DC_VALUES):
        amplitude = matrix(mapping, dc, "final_relative_amplitude")
        half_pf = matrix(mapping, dc, "final_PF_half_mass")
        a1 = axes[0, column].imshow(amplitude, origin="lower", aspect="auto",
                                    cmap="plasma", vmin=0, vmax=2.5)
        a2 = axes[1, column].imshow(half_pf, origin="lower", aspect="auto",
                                    cmap="viridis", vmin=0, vmax=0.8)
        for row, label in ((0, "Relative amplitude"), (1, "Half-mass PF")):
            ax = axes[row, column]
            ax.set_title(f"{label}; $D_c$ = {dc:g}")
            ax.set_xticks(range(len(DIAMETERS)),
                          [diameter_label(d) for d in DIAMETERS], rotation=30)
            ax.set_yticks(range(len(DM_VALUES)), [f"{dm:g}" for dm in DM_VALUES])
            ax.set_xlabel("Diameter (µm)")
            ax.set_ylabel("$D_m$ (µm²/s)")
            ax.grid(False)
    fig.colorbar(a1, ax=axes[0], shrink=0.75, label="Relative amplitude")
    fig.colorbar(a2, ax=axes[1], shrink=0.75, label="Half-mass PF")
    fig.suptitle("Endpoint polarization strength and localization", fontsize=14)
    save(fig, figure_dir / "05_endpoint_polarization_maps.png")


def wgd_scores(mapping: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray,
                                            tuple[float, float]]:
    status = np.zeros((len(DC_VALUES), len(DM_VALUES)))
    bias_ratio = np.zeros_like(status)
    amplitude_difference = np.zeros_like(status)
    best = None
    for i, dc in enumerate(DC_VALUES):
        for j, dm in enumerate(DM_VALUES):
            small = mapping[(dm, dc, 10.0)]
            large = mapping[(dm, dc, DIAMETERS[5])]
            t1 = small["winner_selection_time_s"]
            t2 = large["winner_selection_time_s"]
            if np.isfinite(t1) and not np.isfinite(t2):
                status[i, j] = 1.0
            elif np.isfinite(t1) and np.isfinite(t2):
                status[i, j] = 2.0
            elif not np.isfinite(t1) and np.isfinite(t2):
                status[i, j] = -1.0
            b1 = max(abs(small["mass_bias_final"]), 1e-6)
            b2 = max(abs(large["mass_bias_final"]), 1e-6)
            bias_ratio[i, j] = b1 / b2
            amplitude_difference[i, j] = (
                small["final_relative_amplitude"]
                - large["final_relative_amplitude"]
            )
            qualitative = np.isfinite(t1) and not np.isfinite(t2)
            score = ((1000.0 if qualitative else 0.0)
                     + math.log10(max(bias_ratio[i, j], 1e-6)))
            candidate = (score, dm, dc)
            if best is None or candidate > best:
                best = candidate
    return status, bias_ratio, amplitude_difference, (best[1], best[2])


def wgd_contrast(mapping: dict, figure_dir: Path) -> tuple[float, float]:
    status, bias_ratio, amplitude_difference, candidate = wgd_scores(mapping)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), constrained_layout=True)
    specifications = (
        (axes[0], status, "Winner-selection outcome", "cividis", "status"),
        (axes[1], bias_ratio, "Endpoint bias ratio\n10 / 12.6 µm", "magma", "number"),
        (axes[2], amplitude_difference,
         "Relative-amplitude difference\n10 − 12.6 µm", "viridis", "number"),
    )
    for ax, values, title, cmap, kind in specifications:
        artist = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(DM_VALUES)), [f"{dm:g}" for dm in DM_VALUES])
        ax.set_yticks(range(len(DC_VALUES)), [f"{dc:g}" for dc in DC_VALUES])
        ax.set_xlabel("$D_m$ (µm²/s)")
        ax.set_ylabel("$D_c$ (µm²/s)")
        ax.set_title(title)
        ax.grid(False)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                text = ("both\n>2000" if kind == "status" and values[i, j] == 0
                        else "10 wins;\n12.6 censored" if kind == "status" and values[i, j] == 1
                        else "both win" if kind == "status" and values[i, j] == 2
                        else f"{values[i, j]:.2g}")
                ax.text(j, i, text, ha="center", va="center", fontsize=8,
                        color="white" if values[i, j] <= np.nanmedian(values) else "black")
        fig.colorbar(artist, ax=ax, shrink=0.78)
    fig.suptitle("Geometry-only WGD contrast across diffusion regimes", fontsize=14)
    save(fig, figure_dir / "06_wgd_contrast_maps.png")
    return candidate


def checkpoint_for(mapping: dict, output: Path, dm: float, dc: float,
                   diameter: float) -> dict[str, np.ndarray]:
    task_index = int(mapping[(dm, dc, diameter)]["task_index"])
    with np.load(output / "checkpoints" / f"task_{task_index:03d}.npz",
                 allow_pickle=False) as raw:
        return {key: raw[key].copy() for key in raw.files}


def representative_trajectories(mapping: dict, output: Path, figure_dir: Path,
                                candidate: tuple[float, float]) -> None:
    regimes = ((0.015, 3.0, "Published consensus"),
               (candidate[0], candidate[1], "Strongest screened WGD contrast"))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.3), sharex=True,
                             constrained_layout=True)
    for row_index, (dm, dc, regime_label) in enumerate(regimes):
        for diameter, scale, colour in (
            (10.0, "1x: 10 µm", "#0072B2"),
            (DIAMETERS[5], "2x volume: 12.60 µm", "#D55E00"),
        ):
            data = checkpoint_for(mapping, output, dm, dc, diameter)
            time = data["times"]
            axes[row_index, 0].plot(time, data["competition_mass_bias"],
                                    color=colour, label=scale)
            axes[row_index, 1].step(time, data["metric_n_above_mean_clusters"],
                                    where="post", color=colour)
            axes[row_index, 2].plot(time, data["metric_relative_amplitude"],
                                    color=colour)
        axes[row_index, 0].set_ylabel(
            f"{regime_label}\n$D_m$={dm:g}, $D_c$={dc:g}\nS1–S2 mass bias"
        )
        axes[row_index, 0].legend(frameon=False, fontsize=8)
        axes[row_index, 1].set_ylabel("Domains")
        axes[row_index, 2].set_ylabel("Relative amplitude")
        for ax in axes[row_index]:
            ax.axvspan(0, 10, color="#BBBBBB", alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("Time (s)")
    fig.suptitle("Representative WGD-pair trajectories", fontsize=14)
    save(fig, figure_dir / "07_representative_wgd_trajectories.png")


def numerical_quality(rows: list[dict], figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9), constrained_layout=True)
    axes[0].scatter([row["diameter_um"] for row in rows],
                    [row["actual_max_h_um"] for row in rows], s=12, alpha=0.5)
    axes[0].axhline(0.25, color="#555555", ls="--", label="target h")
    axes[0].set(xlabel="Diameter (µm)", ylabel="Maximum mesh edge (µm)")
    axes[0].legend(frameon=False)
    axes[1].scatter([row["diameter_um"] for row in rows],
                    [row["maximum_relative_mass_error"] for row in rows],
                    s=12, alpha=0.5, color="#009E73")
    axes[1].set(xlabel="Diameter (µm)", ylabel="Maximum relative mass error",
                yscale="log")
    fig.suptitle("Diffusion-screen numerical quality", fontsize=14)
    save(fig, figure_dir / "08_numerical_quality.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "reports/output/giese_wp/diffusion_screen",
    )
    args = parser.parse_args()
    rows = read_rows(args.output)
    mapping = row_map(rows)
    figure_dir = args.output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    style()
    winner_phase_maps(mapping, figure_dir)
    domain_phase_maps(mapping, figure_dir)
    bias_phase_maps(mapping, figure_dir)
    crossover_curves(mapping, figure_dir)
    endpoint_metric_maps(mapping, figure_dir)
    candidate = wgd_contrast(mapping, figure_dir)
    representative_trajectories(mapping, args.output, figure_dir, candidate)
    numerical_quality(rows, figure_dir)
    manifest = {
        "figures": sorted(path.name for path in figure_dir.glob("*.png")),
        "candidate_Dm_um2_s": candidate[0],
        "candidate_Dc_um2_s": candidate[1],
        "source_summary": str(args.output / "diffusion_screen_summary.csv"),
    }
    (args.output / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
