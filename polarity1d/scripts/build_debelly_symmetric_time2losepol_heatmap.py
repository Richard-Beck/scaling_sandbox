#!/usr/bin/env python3
"""Reproduce the Julia Time2LosePol heatmap with the symmetric Python port.

The sweep and loss criterion follow Examples/Heatmap_Time2LosePol.jl.  The
solver intentionally remains the production reflection-symmetric extension:
both endpoints are free, Rac/Rho diffusion is active, and fields are sampled
at finite-volume cell centres.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import DeBelly, DeBellyParameters, NoStimulus, simulate  # noqa: E402


VALUES = np.asarray([0.1 + 0.8 * index for index in range(8)], dtype=float)
LENGTH = 20.0
CELLS = 100
STEP = 1.0
END_TIME = 600.0
CHECK_AFTER = 150.0
RATIO_THRESHOLD = 1.1
WARMUP_STEPS = 251


def opposed_gaussians(alpha0: float, beta0: float, x: np.ndarray) -> dict[str, np.ndarray]:
    width = LENGTH / 2.5
    return {
        "rac": 2.0 * alpha0 * np.exp(-(x**2) / width**2),
        "rho": 2.0 * beta0 * np.exp(-((x - LENGTH) ** 2) / width**2),
    }


def model_parameters(alpha0: float, beta0: float, *, warmup: bool = False) -> DeBellyParameters:
    return DeBellyParameters(
        rac_basal_activation=alpha0,
        rho_basal_activation=beta0,
        rac_mechanical_feedback=0.3 * alpha0,
        rho_mechanical_feedback=0.3 * beta0,
        assay_mode="reflection_symmetric_extension",
        initialization_mode="algebraic_fixed_point",
        input_sampling_mode="cell_centred",
        protrusion_speed=0.0 if warmup else -0.15,
        active_stress=0.0 if warmup else 100.0,
    )


def run_condition(alpha0: float, beta0: float) -> tuple[dict, object]:
    x = (np.arange(CELLS, dtype=float) + 0.5) * LENGTH / CELLS
    warmup = simulate(
        DeBelly(model_parameters(alpha0, beta0, warmup=True)),
        protocol=NoStimulus(),
        length=LENGTH,
        cells=CELLS,
        times=np.asarray([0.0, WARMUP_STEPS * STEP]),
        step=STEP,
        initial_fields=opposed_gaussians(alpha0, beta0, x),
    )
    initial_fields = {
        "rac": warmup.fields["rac"][-1].copy(),
        "rho": warmup.fields["rho"][-1].copy(),
    }
    trajectory = simulate(
        DeBelly(model_parameters(alpha0, beta0)),
        protocol=NoStimulus(),
        length=LENGTH,
        cells=CELLS,
        end_time=END_TIME,
        output_interval=1.0,
        step=STEP,
        initial_fields=initial_fields,
    )

    # Julia uses its second and penultimate nodal values.  The closest FV
    # analogue is the second cell centre from each boundary.
    rac_ratio = trajectory.fields["rac"][:, 1] / trajectory.fields["rac"][:, -2]
    rho_ratio = trajectory.fields["rho"][:, -2] / trajectory.fields["rho"][:, 1]
    eligible = trajectory.time > CHECK_AFTER
    rac_failed = rac_ratio < RATIO_THRESHOLD
    rho_failed = rho_ratio < RATIO_THRESHOLD
    failed = eligible & (rac_failed | rho_failed)
    indices = np.flatnonzero(failed)
    censored = len(indices) == 0
    loss_index = int(indices[0]) if not censored else len(trajectory.time) - 1
    initially_polarized = bool(
        rac_ratio[0] >= RATIO_THRESHOLD and rho_ratio[0] >= RATIO_THRESHOLD
    )
    if not initially_polarized:
        audit_outcome = "not_polarized_at_t0"
    elif censored:
        audit_outcome = "censored"
    elif rac_failed[loss_index] and rho_failed[loss_index]:
        audit_outcome = "both"
    elif rac_failed[loss_index]:
        audit_outcome = "rac"
    else:
        audit_outcome = "rho"
    if censored:
        failure_mode = "censored"
    elif rac_failed[loss_index] and rho_failed[loss_index]:
        failure_mode = "both"
    elif rac_failed[loss_index]:
        failure_mode = "rac"
    else:
        failure_mode = "rho"
    row = {
        "alpha0": alpha0,
        "beta0": beta0,
        "alpha_mechanical": 0.3 * alpha0,
        "beta_mechanical": 0.3 * beta0,
        "loss_time_s": float(trajectory.time[loss_index]),
        "julia_plotted_time_s": float(trajectory.time[loss_index] - 50.0),
        "censored_at_600_s": censored,
        "initially_polarized": initially_polarized,
        "audit_outcome": audit_outcome,
        "failure_mode": failure_mode,
        "rac_ratio_t0": float(rac_ratio[0]),
        "rho_ratio_t0": float(rho_ratio[0]),
        "rac_ratio_t150": float(rac_ratio[int(CHECK_AFTER)]),
        "rho_ratio_t150": float(rho_ratio[int(CHECK_AFTER)]),
        "rac_ratio_final": float(rac_ratio[-1]),
        "rho_ratio_final": float(rho_ratio[-1]),
        "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
    }
    return row, trajectory


def draw(rows: list[dict], destination: Path) -> None:
    loss = np.empty((len(VALUES), len(VALUES)))
    status = np.empty_like(loss, dtype=int)
    modes = {
        "not_polarized_at_t0": 0,
        "censored": 1,
        "rac": 2,
        "rho": 3,
        "both": 4,
    }
    for row in rows:
        i = int(np.flatnonzero(np.isclose(VALUES, row["alpha0"]))[0])
        j = int(np.flatnonzero(np.isclose(VALUES, row["beta0"]))[0])
        loss[i, j] = row["julia_plotted_time_s"]
        status[i, j] = modes[row["audit_outcome"]]

    mpl.rcParams.update({"font.size": 8, "axes.linewidth": 0.8})
    figure, axes = plt.subplots(
        1, 2, figsize=(7.0, 3.15), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.18, 1.0]},
    )
    extent = (VALUES[0] - 0.4, VALUES[-1] + 0.4, VALUES[0] - 0.4, VALUES[-1] + 0.4)
    heat = axes[0].imshow(
        loss, origin="lower", extent=extent, interpolation="nearest",
        cmap="viridis", vmin=0.0, vmax=END_TIME - 50.0, aspect="equal",
    )
    axes[0].set_title("A  Julia-style persistence time")
    colorbar = figure.colorbar(heat, ax=axes[0], fraction=0.047, pad=0.03)
    colorbar.set_label(r"reported time, $t_{loss}-50$ (s)")

    status_colors = ["#bdbdbd", "#4d4d4d", "#377eb8", "#e41a1c", "#984ea3"]
    status_cmap = mpl.colors.ListedColormap(status_colors)
    status_norm = mpl.colors.BoundaryNorm(np.arange(-0.5, 5.5), status_cmap.N)
    axes[1].imshow(
        status, origin="lower", extent=extent, interpolation="nearest",
        cmap=status_cmap, norm=status_norm, aspect="equal",
    )
    axes[1].set_title("B  Audited assay outcome")
    handles = [
        mpl.patches.Patch(color=color, label=label)
        for color, label in zip(
            status_colors,
            (
                "not polarized at t=0",
                "censored at 600 s",
                "lost: Rac criterion",
                "lost: Rho criterion",
                "lost: both criteria",
            ),
        )
    ]
    axes[1].legend(handles=handles, loc="upper left", fontsize=7, frameon=False)

    for axis in axes:
        axis.set_xlabel(r"Rho basal activation, $\beta_0$")
        axis.set_ylabel(r"Rac basal activation, $\alpha_0$")
        axis.set_xticks(VALUES)
        axis.set_yticks(VALUES)
        axis.tick_params(labelsize=7)
    figure.suptitle(
        "2026 de Belly reflection-symmetric Python extension (both endpoints free)\n"
        r"$L=20\,\mu$m; $\alpha=0.3\alpha_0$; $\beta=0.3\beta_0$; no cue",
        fontsize=9,
    )
    figure.savefig(destination, dpi=300)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_dir", nargs="?", type=Path,
        default=ROOT / "results" / "tmp" / "trial" / "debelly_symmetric_time2losepol",
    )
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    arrays = {}
    for alpha0 in VALUES:
        for beta0 in VALUES:
            row, trajectory = run_condition(float(alpha0), float(beta0))
            rows.append(row)
            key = f"alpha{alpha0:.1f}_beta{beta0:.1f}"
            arrays[f"{key}_time"] = trajectory.time
            arrays[f"{key}_rac_orientation"] = (
                trajectory.fields["rac"][:, 1] / trajectory.fields["rac"][:, -2]
            )
            arrays[f"{key}_rho_orientation"] = (
                trajectory.fields["rho"][:, -2] / trajectory.fields["rho"][:, 1]
            )
            print(
                f"alpha0={alpha0:.1f} beta0={beta0:.1f} "
                f"loss={row['loss_time_s']:.0f} mode={row['failure_mode']}"
            )

    csv_path = args.output_dir / "time2losepol_python_symmetric.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(args.output_dir / "time2losepol_endpoint_ratios.npz", **arrays)
    metadata = {
        "model": "reflection-symmetric extension of the 2026 de Belly mechanochemical model",
        "solver": "Python finite-volume port; both endpoints free",
        "upstream_example": "Examples/Heatmap_Time2LosePol.jl",
        "literal_matches": {
            "length_um": LENGTH,
            "cells": CELLS,
            "dt_s": STEP,
            "end_time_s": END_TIME,
            "basal_activation_values": VALUES.tolist(),
            "mechanical_feedback_fraction": 0.3,
            "opposed_gaussian_width": "L/2.5",
            "chemical_warmup_steps": WARMUP_STEPS,
            "loss_ratio_threshold": RATIO_THRESHOLD,
            "loss_check_after_s": CHECK_AFTER,
            "plotted_offset_s": -50.0,
        },
        "intentional_python_port_differences": [
            "reflection-symmetric free boundary mechanics instead of Julia's left-only protrusion",
            "cell-centred finite volumes instead of nodal finite elements",
            "Rac and Rho diffusion 0.5 is active; the Julia heatmap declares but does not use Drac/Drho",
            "protrusion threshold follows the Python port's homogeneous-reference rule",
        ],
        "censoring": "loss_time=600 means right-censored, retained explicitly in panel B",
        "outcome_counts": {
            outcome: sum(row["audit_outcome"] == outcome for row in rows)
            for outcome in (
                "not_polarized_at_t0", "censored", "rac", "rho", "both"
            )
        },
    }
    (args.output_dir / "PROVENANCE.json").write_text(json.dumps(metadata, indent=2) + "\n")
    draw(rows, args.output_dir / "Heatmap_Time2LosePol_python_symmetric.png")


if __name__ == "__main__":
    main()
