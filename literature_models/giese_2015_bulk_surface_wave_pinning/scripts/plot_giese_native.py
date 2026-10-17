#!/usr/bin/env python3
"""Create compact production figures from collected native Giese checkpoints."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DIAMETERS = (1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 15.0)
COLOURS = {"single": "#0072B2", "two_site": "#D55E00"}
LABELS = {"single": "Single stimulus", "two_site": "Two opposing stimuli"}


def checkpoint_path(output: Path, mode: str, diameter: float) -> Path:
    return output / "checkpoints" / f"{mode}_diameter_{diameter:04.1f}.npz"


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


def read_summaries(output: Path) -> list[dict]:
    with (output / "native_size_summary.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    numeric = {
        "diameter_um", "final_PF", "final_POL", "final_relative_amplitude",
        "polarization_time_90_max_PF_s", "maximum_relative_mass_error",
        "winner_selection_time_s",
        "final_half_mass_PF", "maximal_half_mass_PF",
        "polarization_time_90_final_half_mass_PF_s",
    }
    for row in rows:
        for key in numeric:
            row[key] = float(row[key])
    return rows


def size_summary_figure(output: Path, figure_dir: Path) -> None:
    rows = read_summaries(output)
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 6.8), constrained_layout=True)
    specifications = (
        ("final_POL", "Final POL", None),
        ("final_PF", "Final PF", (0, 1)),
        ("final_relative_amplitude", "Final relative membrane range", None),
        ("winner_selection_time_s", "Two-site winner-selection time (s)", None),
    )
    for ax, (key, ylabel, ylim) in zip(axes.flat, specifications):
        modes = ("two_site",) if key == "winner_selection_time_s" else ("single", "two_site")
        for mode in modes:
            chosen = sorted((row for row in rows if row["mode"] == mode),
                            key=lambda row: row["diameter_um"])
            x = np.asarray([row["diameter_um"] for row in chosen])
            y = np.asarray([row[key] for row in chosen])
            finite = np.isfinite(y)
            ax.plot(x[finite], y[finite], marker="o", lw=2,
                    color=COLOURS[mode], label=LABELS[mode])
            if key == "winner_selection_time_s":
                # The 1.5-um condition never polarizes; it is not a censored
                # winner-selection event.  Remaining non-finite cases retained
                # two sites through the 2000-s observation window.
                censored = (~finite) & (x >= 3.0)
                ax.scatter(x[censored], np.full(censored.sum(), 2000.0), marker="^",
                           s=65, facecolors="none", edgecolors=COLOURS[mode],
                           label="No unique site by 2000 s")
                ax.plot((3.0, 6.0), (43.0, 953.0), "o--", color="#555555",
                        lw=1.2, label="Published approximate targets")
                ax.set_yscale("log")
        ax.set_xlabel("Cell diameter (µm)")
        ax.set_ylabel(ylabel)
        ax.set_xticks(DIAMETERS)
        if ylim is not None:
            ax.set_ylim(*ylim)
    axes[0, 0].legend(frameon=False)
    axes[1, 1].legend(frameon=False, fontsize=8)
    fig.suptitle("Native Giese bulk–surface wave-pinning size sweep", fontsize=14)
    fig.savefig(figure_dir / "native_size_metrics.png", bbox_inches="tight")
    plt.close(fig)


def membrane_profile_figure(output: Path, figure_dir: Path) -> None:
    conditions = (("single", 3.0), ("single", 15.0),
                  ("two_site", 3.0), ("two_site", 6.0))
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), sharex=True,
                             constrained_layout=True)
    for ax, (mode, diameter) in zip(axes.flat, conditions):
        with np.load(checkpoint_path(output, mode, diameter), allow_pickle=False) as data:
            times = data["times"]
            desired = (0.0, 10.0, 50.0, times[-1])
            colours = plt.cm.viridis(np.linspace(0.1, 0.9, len(desired)))
            order = np.argsort(data["boundary_angles"])
            angle = data["boundary_angles"][order]
            for target, colour in zip(desired, colours):
                index = int(np.argmin(np.abs(times - target)))
                ax.plot(angle, data["membrane"][index, order], color=colour,
                        lw=1.8, label=f"{times[index]:g} s")
        ax.set_title(f"{LABELS[mode]}, d = {diameter:g} µm")
        ax.set_ylabel("Membrane u (µm·µM)")
        ax.set_xlim(0, 2 * np.pi)
        ax.set_xticks((0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi),
                      ("0", "π/2", "π", "3π/2", "2π"))
    for ax in axes[-1]:
        ax.set_xlabel("Membrane angle")
    axes[0, 0].legend(frameon=False, ncol=2)
    fig.suptitle("Representative membrane profiles", fontsize=14)
    fig.savefig(figure_dir / "representative_membrane_profiles.png",
                bbox_inches="tight")
    plt.close(fig)


def bulk_field_figure(output: Path, figure_dir: Path) -> None:
    rows = (
        ("single", 3.0, (0.0, 10.0, 50.0)),
        ("two_site", 3.0, (10.0, 50.0, 200.0)),
        ("two_site", 6.0, (10.0, 500.0, 2000.0)),
    )
    fig, axes = plt.subplots(3, 3, figsize=(9.3, 8.5), constrained_layout=True)
    for row_index, (mode, diameter, desired) in enumerate(rows):
        with np.load(checkpoint_path(output, mode, diameter), allow_pickle=False) as data:
            points = data["mesh_points"]
            triangulation = mtri.Triangulation(
                points[0], points[1], data["mesh_triangles"].T
            )
            selected_times = data["selected_times"]
            fields = data["selected_bulk"]
            indices = [int(np.argmin(np.abs(selected_times - target)))
                       for target in desired]
            vmin = min(float(fields[index].min()) for index in indices)
            vmax = max(float(fields[index].max()) for index in indices)
            for column, (target, index) in enumerate(zip(desired, indices)):
                ax = axes[row_index, column]
                artist = ax.tripcolor(triangulation, fields[index], shading="gouraud",
                                      cmap="magma", vmin=vmin, vmax=vmax)
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.grid(False)
                ax.set_title(f"t = {selected_times[index]:g} s")
                if column == 0:
                    ax.set_ylabel(f"{LABELS[mode]}\nd = {diameter:g} µm")
            fig.colorbar(artist, ax=axes[row_index, :], shrink=0.72,
                         label="Cytoplasmic v (µM)")
    fig.suptitle("Selected cytoplasmic concentration fields", fontsize=14)
    fig.savefig(figure_dir / "representative_bulk_fields.png", bbox_inches="tight")
    plt.close(fig)


def convergence_figure(output: Path, figure_dir: Path) -> None:
    validation = output.parent / "validation"
    with (validation / "mesh_convergence.csv").open(newline="") as handle:
        mesh = list(csv.DictReader(handle))
    with (validation / "time_convergence.csv").open(newline="") as handle:
        time = list(csv.DictReader(handle))
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.7), constrained_layout=True)
    axes[0].plot([float(row["actual_max_h_um"]) for row in mesh],
                 [float(row["final_POL"]) for row in mesh], "o-", color="#009E73")
    axes[0].invert_xaxis()
    axes[0].set_xlabel("Actual maximum mesh edge h (µm)")
    axes[0].set_ylabel("Final POL")
    axes[0].set_title("Mesh convergence")
    axes[1].plot([float(row["dt_s"]) for row in time],
                 [float(row["final_POL"]) for row in time], "o-", color="#CC79A7")
    axes[1].invert_xaxis()
    axes[1].set_xlabel("Time step Δt (s)")
    axes[1].set_ylabel("Final POL")
    axes[1].set_title("Time-step convergence")
    fig.suptitle("Numerical convergence at d = 3 µm", fontsize=13)
    fig.savefig(figure_dir / "convergence.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=ROOT / "reports/output/giese_wp/native")
    parser.add_argument("--figure-dir", type=Path,
                        default=ROOT / "reports/assets/giese_wp/v1")
    args = parser.parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    style()
    size_summary_figure(args.output, args.figure_dir)
    membrane_profile_figure(args.output, args.figure_dir)
    bulk_field_figure(args.output, args.figure_dir)
    convergence_figure(args.output, args.figure_dir)
    copied = {
        args.output / "native_size_summary.csv": args.figure_dir / "native_size_summary.csv",
        args.output.parent / "validation/mesh_convergence.csv": args.figure_dir / "mesh_convergence.csv",
        args.output.parent / "validation/time_convergence.csv": args.figure_dir / "time_convergence.csv",
        args.output.parent / "validation/validation_summary.json": args.figure_dir / "validation_summary.json",
    }
    diagnostic = args.output.parent / "fix_diagnostics.csv"
    if diagnostic.exists():
        copied[diagnostic] = args.figure_dir / "fix_diagnostics.csv"
    for source, destination in copied.items():
        shutil.copyfile(source, destination)
    manifest = {
        "figures": sorted(path.name for path in args.figure_dir.glob("*.png")),
        "source": str(args.output),
        "tables": sorted(path.name for path in copied.values()),
    }
    (args.figure_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
