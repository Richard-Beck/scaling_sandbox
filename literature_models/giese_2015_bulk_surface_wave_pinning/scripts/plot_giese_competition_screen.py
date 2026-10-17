#!/usr/bin/env python3
"""Plot the coupled-diffusion, mass, stimulus, and single-site Giese screen."""

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
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DIFFUSION_PAIRS = ((0.05, 10.0), (0.075, 15.0), (0.10, 20.0))
MASS_SCALES = (0.7, 0.8, 0.9, 1.0)
WEAK_STIMULI = (0.40, 0.35, 0.30)
DIAMETERS = (10.0, 10.0 * 2.0 ** (1.0 / 3.0))
END_TIME = 2000.0


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
    with (output / "competition_screen_summary.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    text_fields = {"mode"}
    boolean_fields = {"winner_selection_right_censored"}
    for row in rows:
        for key, value in tuple(row.items()):
            if key in text_fields:
                continue
            if key in boolean_fields:
                row[key] = value.lower() == "true"
            else:
                row[key] = float(value)
    return rows


def index_rows(rows: list[dict]) -> dict[tuple, dict]:
    result = {}
    for row in rows:
        weak = None if row["mode"] == "single" else row["weak_stimulus_um_s"]
        result[(row["mode"], row["Dm_um2_s"], row["Dc_um2_s"],
                row["mass_scale"], weak, row["diameter_um"])] = row
    return result


def save(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def json_safe(value):
    """Convert NumPy scalars and non-finite floats to strict JSON values."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def setup_heatmap(ax: plt.Axes, title: str) -> None:
    ax.set_title(title)
    ax.set_xticks(range(len(DIFFUSION_PAIRS)), [
        f"{dm:g}/{dc:g}" for dm, dc in DIFFUSION_PAIRS
    ])
    ax.set_yticks(range(len(MASS_SCALES)), [f"{mass:g}" for mass in MASS_SCALES])
    ax.set_xlabel("$D_m/D_c$ (µm²/s)")
    ax.set_ylabel("Equilibrium mass scale")
    ax.grid(False)


def condition_matrix(index: dict, mode: str, diameter: float, key: str,
                     weak: float | None = None) -> np.ndarray:
    return np.asarray([
        [index[(mode, dm, dc, mass, weak, diameter)][key]
         for dm, dc in DIFFUSION_PAIRS]
        for mass in MASS_SCALES
    ], dtype=float)


def single_time_maps(index: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    time_matrices = []
    for ax, diameter, label in zip(
        axes[:2], DIAMETERS, ("1x: 10 µm", "2x volume: 12.6 µm")
    ):
        values = condition_matrix(
            index, "single", diameter,
            "polarization_time_90_max_half_mass_PF_s",
        )
        final_amplitude = condition_matrix(
            index, "single", diameter, "final_relative_amplitude"
        )
        values[final_amplitude < 0.05] = np.nan
        time_matrices.append(values)
        artist = ax.imshow(values, origin="lower", aspect="auto", cmap="viridis")
        setup_heatmap(ax, label)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                ax.text(j, i, "no pol." if not np.isfinite(values[i, j])
                        else f"{values[i, j]:.0f}", ha="center", va="center",
                        fontsize=8, color="white")
        fig.colorbar(artist, ax=ax, shrink=0.78, label="90%-max half-mass PF time (s)")
    ratio = time_matrices[1] / time_matrices[0]
    artist = axes[2].imshow(ratio, origin="lower", aspect="auto", cmap="magma")
    setup_heatmap(axes[2], "WGD polarization-time ratio\n12.6 / 10 µm")
    for i in range(ratio.shape[0]):
        for j in range(ratio.shape[1]):
            axes[2].text(j, i, "—" if not np.isfinite(ratio[i, j])
                         else f"{ratio[i, j]:.2f}", ha="center", va="center",
                         color="white", fontsize=8)
    fig.colorbar(artist, ax=axes[2], shrink=0.78,
                 label="Persistent-polarization time ratio")
    fig.suptitle("Single-stimulus polarization kinetics; transient-only responses marked no pol.",
                 fontsize=14)
    save(fig, figure_dir / "01_single_stimulus_polarization_times.png")


def single_strength_maps(index: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.7), constrained_layout=True)
    for column, (diameter, label) in enumerate(zip(
        DIAMETERS, ("1x: 10 µm", "2x volume: 12.6 µm")
    )):
        max_pf = condition_matrix(index, "single", diameter, "maximal_half_mass_PF")
        final_amp = condition_matrix(index, "single", diameter, "final_relative_amplitude")
        a1 = axes[0, column].imshow(max_pf, origin="lower", aspect="auto",
                                    cmap="viridis", vmin=0, vmax=0.8)
        a2 = axes[1, column].imshow(final_amp, origin="lower", aspect="auto",
                                    cmap="plasma", vmin=0, vmax=2.5)
        setup_heatmap(axes[0, column], f"Maximal half-mass PF; {label}")
        setup_heatmap(axes[1, column], f"Final relative amplitude; {label}")
    fig.colorbar(a1, ax=axes[0], shrink=0.75, label="Maximal half-mass PF")
    fig.colorbar(a2, ax=axes[1], shrink=0.75, label="Final relative amplitude")
    fig.suptitle("Single-site polarization strength", fontsize=14)
    save(fig, figure_dir / "02_single_stimulus_strength.png")


def outcome_code(row: dict) -> int:
    if np.isfinite(row["winner_selection_time_s"]):
        return 2
    if row["winner_selection_right_censored"]:
        return 1
    return 0


def two_site_outcomes(index: dict, figure_dir: Path) -> None:
    cmap = ListedColormap(("#666666", "#E69F00", "#009E73"))
    norm = BoundaryNorm((-0.5, 0.5, 1.5, 2.5), cmap.N)
    fig, axes = plt.subplots(3, 2, figsize=(10, 11), constrained_layout=True)
    for row_index, weak in enumerate(WEAK_STIMULI):
        for column, (diameter, label) in enumerate(zip(
            DIAMETERS, ("1x: 10 µm", "2x volume: 12.6 µm")
        )):
            values = np.asarray([
                [outcome_code(index[("two_site", dm, dc, mass, weak, diameter)])
                 for dm, dc in DIFFUSION_PAIRS]
                for mass in MASS_SCALES
            ])
            ax = axes[row_index, column]
            artist = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap, norm=norm)
            setup_heatmap(ax, f"S1/S2 = 0.44/{weak:.2f}; {label}")
            for i, mass in enumerate(MASS_SCALES):
                for j, (dm, dc) in enumerate(DIFFUSION_PAIRS):
                    condition = index[("two_site", dm, dc, mass, weak, diameter)]
                    code = values[i, j]
                    text = (f"{condition['winner_selection_time_s']:.0f} s" if code == 2
                            else ">2000" if code == 1 else "no pol.")
                    ax.text(j, i, text, ha="center", va="center", fontsize=8,
                            color="white" if code != 1 else "black")
    colourbar = fig.colorbar(artist, ax=axes, shrink=0.65, ticks=(0, 1, 2))
    colourbar.ax.set_yticklabels(("Nonpolarized", "Censored", "Winner"))
    fig.suptitle("Two-site winner-selection outcomes", fontsize=14)
    save(fig, figure_dir / "03_two_site_outcome_maps.png")


def desired_split_maps(index: dict, figure_dir: Path) -> None:
    cmap = ListedColormap(("#BBBBBB", "#56B4E9", "#D55E00"))
    norm = BoundaryNorm((-0.5, 0.5, 1.5, 2.5), cmap.N)
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), constrained_layout=True)
    for ax, weak in zip(axes, WEAK_STIMULI):
        values = np.zeros((len(MASS_SCALES), len(DIFFUSION_PAIRS)))
        for i, mass in enumerate(MASS_SCALES):
            for j, (dm, dc) in enumerate(DIFFUSION_PAIRS):
                one = index[("two_site", dm, dc, mass, weak, DIAMETERS[0])]
                two = index[("two_site", dm, dc, mass, weak, DIAMETERS[1])]
                if np.isfinite(one["winner_selection_time_s"]) and two[
                    "winner_selection_right_censored"
                ]:
                    values[i, j] = 2
                elif np.isfinite(one["winner_selection_time_s"]) and np.isfinite(
                    two["winner_selection_time_s"]
                ):
                    values[i, j] = 1
        artist = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap, norm=norm)
        setup_heatmap(ax, f"S1/S2 = 0.44/{weak:.2f}")
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                text = ("desired\nsplit" if values[i, j] == 2 else
                        "both win" if values[i, j] == 1 else "no split")
                ax.text(j, i, text, ha="center", va="center", fontsize=8)
    colourbar = fig.colorbar(artist, ax=axes, shrink=0.75, ticks=(0, 1, 2))
    colourbar.ax.set_yticklabels(("No split", "Both win", "Desired split"))
    fig.suptitle("Does geometry-only WGD cross the competition boundary?", fontsize=14)
    save(fig, figure_dir / "04_desired_wgd_split_maps.png")


def select_candidate(index: dict) -> dict:
    best = None
    for weak in WEAK_STIMULI:
        for mass in MASS_SCALES:
            for dm, dc in DIFFUSION_PAIRS:
                single_one = index[("single", dm, dc, mass, None, DIAMETERS[0])]
                single_two = index[("single", dm, dc, mass, None, DIAMETERS[1])]
                if min(single_one["final_relative_amplitude"],
                       single_two["final_relative_amplitude"]) < 0.05:
                    continue
                one = index[("two_site", dm, dc, mass, weak, DIAMETERS[0])]
                two = index[("two_site", dm, dc, mass, weak, DIAMETERS[1])]
                qualitative = (np.isfinite(one["winner_selection_time_s"])
                               and two["winner_selection_right_censored"])
                both_win = (np.isfinite(one["winner_selection_time_s"])
                            and np.isfinite(two["winner_selection_time_s"]))
                time_score = (END_TIME - one["winner_selection_time_s"]
                              if np.isfinite(one["winner_selection_time_s"]) else 0.0)
                bias_ratio = (max(abs(one["mass_bias_final"]), 1e-6)
                              / max(abs(two["mass_bias_final"]), 1e-6))
                score = ((1e6 if qualitative else 1e4 if both_win else 0.0)
                         + time_score + 10.0 * math.log(max(bias_ratio, 1e-6)))
                candidate = {
                    "score": float(score),
                    "qualitative_split": bool(qualitative),
                    "Dm_um2_s": float(dm), "Dc_um2_s": float(dc),
                    "mass_scale": float(mass),
                    "weak_stimulus_um_s": float(weak),
                    "winner_time_1x_s": float(one["winner_selection_time_s"]),
                    "winner_time_2x_s": float(two["winner_selection_time_s"]),
                    "winner_2x_censored": bool(two["winner_selection_right_censored"]),
                    "single_time_1x_s": float(single_one[
                        "polarization_time_90_max_half_mass_PF_s"
                    ]),
                    "single_time_2x_s": float(single_two[
                        "polarization_time_90_max_half_mass_PF_s"
                    ]),
                    "bias_ratio_1x_to_2x": float(bias_ratio),
                }
                if best is None or candidate["score"] > best["score"]:
                    best = candidate
    if best is None:
        raise RuntimeError("No candidate retained single-site polarization")
    return best


def checkpoint(index: dict, output: Path, mode: str, dm: float, dc: float,
               mass: float, weak: float | None, diameter: float) -> dict:
    task = int(index[(mode, dm, dc, mass, weak, diameter)]["task_index"])
    with np.load(output / "checkpoints" / f"task_{task:03d}.npz",
                 allow_pickle=False) as raw:
        return {key: raw[key].copy() for key in raw.files}


def representative_trajectories(index: dict, output: Path, figure_dir: Path,
                                candidate: dict) -> None:
    dm, dc = candidate["Dm_um2_s"], candidate["Dc_um2_s"]
    mass, weak = candidate["mass_scale"], candidate["weak_stimulus_um_s"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    for diameter, label, colour in (
        (DIAMETERS[0], "1x: 10 µm", "#0072B2"),
        (DIAMETERS[1], "2x volume: 12.6 µm", "#D55E00"),
    ):
        single = checkpoint(index, output, "single", dm, dc, mass, None, diameter)
        two = checkpoint(index, output, "two_site", dm, dc, mass, weak, diameter)
        axes[0].plot(single["times"], single["metric_half_mass_PF"],
                     color=colour, label=label)
        axes[1].plot(two["times"], two["competition_mass_bias"], color=colour)
        axes[2].step(two["times"], two["metric_n_above_mean_clusters"],
                     where="post", color=colour)
    axes[0].set(xlabel="Time (s)", ylabel="Single-site half-mass PF")
    axes[0].legend(frameon=False)
    axes[1].set(xlabel="Time (s)", ylabel="Two-site S1–S2 mass bias")
    axes[2].set(xlabel="Time (s)", ylabel="Two-site domains")
    fig.suptitle(
        f"Selected regime: Dm/Dc={dm:g}/{dc:g}, mass={mass:g}, "
        f"S1/S2=0.44/{weak:.2f}", fontsize=13
    )
    save(fig, figure_dir / "05_selected_regime_trajectories.png")


def candidate_profiles(index: dict, output: Path, figure_dir: Path,
                       candidate: dict) -> None:
    dm, dc = candidate["Dm_um2_s"], candidate["Dc_um2_s"]
    mass, weak = candidate["mass_scale"], candidate["weak_stimulus_um_s"]
    times = (10.0, 100.0, 500.0, 2000.0)
    grid = np.linspace(0, 2 * math.pi, 500)
    fig, axes = plt.subplots(2, len(times), figsize=(13.5, 6.0), sharex=True,
                             constrained_layout=True)
    for row, (diameter, label, colour) in enumerate((
        (DIAMETERS[0], "1x: 10 µm", "#0072B2"),
        (DIAMETERS[1], "2x volume: 12.6 µm", "#D55E00"),
    )):
        data = checkpoint(index, output, "two_site", dm, dc, mass, weak, diameter)
        order = np.argsort(data["boundary_angles"])
        angle = data["boundary_angles"][order]
        for ax, time in zip(axes[row], times):
            time_index = int(np.argmin(np.abs(data["times"] - time)))
            profile = np.interp(grid, angle, data["membrane"][time_index, order],
                                period=2 * math.pi)
            ax.plot(grid, profile, color=colour, lw=2)
            ax.set_title(f"t = {time:g} s")
            ax.set_xticks((0, math.pi, 2 * math.pi), ("0", "π", "2π"))
            if ax is axes[row, 0]:
                ax.set_ylabel(label + "\nmembrane u")
    fig.suptitle("Selected-regime two-site membrane profiles", fontsize=14)
    save(fig, figure_dir / "06_selected_regime_profiles.png")


def abundance_states(rows: list[dict], figure_dir: Path) -> None:
    chosen = {}
    for row in rows:
        key = (row["mass_scale"], row["diameter_um"])
        chosen.setdefault(key, row)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), constrained_layout=True)
    for diameter, label, colour in (
        (DIAMETERS[0], "1x: 10 µm", "#0072B2"),
        (DIAMETERS[1], "2x volume: 12.6 µm", "#D55E00"),
    ):
        selected = [chosen[(mass, diameter)] for mass in MASS_SCALES]
        axes[0].plot(MASS_SCALES, [row["u0"] for row in selected], "o-",
                     color=colour, label=label)
        axes[1].plot(MASS_SCALES, [row["v0"] for row in selected], "o-",
                     color=colour)
        axes[2].plot(MASS_SCALES, [row["target_native_mass"] for row in selected],
                     "o-", color=colour)
    axes[0].set(xlabel="Mass scale", ylabel="Equilibrium membrane u0")
    axes[0].legend(frameon=False)
    axes[1].set(xlabel="Mass scale", ylabel="Equilibrium bulk v0")
    axes[2].set(xlabel="Mass scale", ylabel="Target native mass")
    fig.suptitle("Reaction-balanced initial states", fontsize=14)
    save(fig, figure_dir / "07_equilibrium_abundance_states.png")


def numerical_quality(rows: list[dict], figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), constrained_layout=True)
    axes[0].scatter([row["target_native_mass"] for row in rows],
                    [row["actual_initial_mass"] / row["target_native_mass"] - 1
                     for row in rows], s=12, alpha=0.5)
    axes[0].set(xlabel="Target native mass", ylabel="Initial mass relative error")
    axes[1].scatter([row["diameter_um"] for row in rows],
                    [row["maximum_relative_mass_error"] for row in rows],
                    s=12, alpha=0.5, color="#009E73")
    axes[1].set(xlabel="Diameter (µm)", ylabel="Maximum dynamic mass error",
                yscale="log")
    fig.suptitle("Competition-screen numerical quality", fontsize=14)
    save(fig, figure_dir / "08_numerical_quality.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "reports/output/giese_wp/competition_screen",
    )
    args = parser.parse_args()
    rows = read_rows(args.output)
    index = index_rows(rows)
    figure_dir = args.output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    style()
    single_time_maps(index, figure_dir)
    single_strength_maps(index, figure_dir)
    two_site_outcomes(index, figure_dir)
    desired_split_maps(index, figure_dir)
    candidate = select_candidate(index)
    representative_trajectories(index, args.output, figure_dir, candidate)
    candidate_profiles(index, args.output, figure_dir, candidate)
    abundance_states(rows, figure_dir)
    numerical_quality(rows, figure_dir)
    candidate = json_safe(candidate)
    (args.output / "selected_candidate.json").write_text(
        json.dumps(candidate, indent=2, allow_nan=False) + "\n"
    )
    manifest = {
        "figures": sorted(path.name for path in figure_dir.glob("*.png")),
        "selected_candidate": candidate,
        "source_summary": str(args.output / "competition_screen_summary.csv"),
    }
    (args.output / "plot_manifest.json").write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(manifest, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
