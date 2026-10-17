#!/usr/bin/env python3
"""Find a finite opposite Rac pulse that breaks the selected locked state."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly, DeBellyParameters, FinitePulseSequence, simulate,
)
from polarity1d.campaign import _pulse_switching_summary  # noqa: E402


OUTPUT = ROOT / "results" / "tmp" / "trial" / "debelly_opposite_pulse_breaking"
ALPHA0, BETA0 = 3.5, 4.0
WIDTHS = (0.1, 0.2, 0.4, 0.5)
RATIOS = (1.0, 2.0, 4.0, 8.0)
DURATIONS = (80.0, 160.0, 320.0, 640.0, 1_280.0)
OPPOSITE_START = 5_040.0


def _pulse(side, start, duration, amplitude, width):
    return {
        "amplitude": amplitude, "start": start, "stop": start + duration,
        "taper": duration / 2.0, "width_fraction": width, "sides": (side,),
    }


def _run(condition):
    width, ratio, duration = condition
    model = DeBelly(
        DeBellyParameters(
            assay_mode="reflection_symmetric_extension",
            initialization_mode="algebraic_fixed_point",
            input_sampling_mode="cell_centred",
            rac_basal_activation=ALPHA0,
            rho_basal_activation=BETA0,
            rac_mechanical_feedback=0.3 * ALPHA0,
            rho_mechanical_feedback=0.3 * BETA0,
        )
    )
    protocol = FinitePulseSequence(
        pulses=(
            _pulse("left", 0.0, 40.0, ALPHA0, 0.1),
            _pulse("right", OPPOSITE_START, duration, ratio * ALPHA0, width),
        )
    )
    end_time = OPPOSITE_START + duration + 3_000.0
    trajectory = simulate(
        model, protocol=protocol, length=15.0, cells=75,
        end_time=end_time, output_interval=5.0, step=0.5,
    )
    metrics = _pulse_switching_summary(
        trajectory,
        events=[
            {"side": "left", "start": 0.0, "stop": 40.0},
            {"side": "right", "start": OPPOSITE_START, "stop": OPPOSITE_START + duration},
        ],
        orientation_threshold=0.05,
        sustain_s=100.0,
        minimum_persistence_s=1_000.0,
    )
    return {
        "opposite_width_fraction": width,
        "opposite_amplitude_ratio": ratio,
        "opposite_duration_s": duration,
        "induction_success": metrics["induction_success"],
        "switch_success": metrics["switch_1_success"],
        "switch_latency_from_start_s": metrics["switch_1_latency_from_pulse_start_s"],
        "switch_latency_after_withdrawal_s": metrics["switch_1_latency_after_withdrawal_s"],
        "post_pulse_persistence_s": metrics["switch_1_persistence_s"],
        "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
    }


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    conditions = [
        (width, ratio, duration)
        for width in WIDTHS for ratio in RATIOS for duration in DURATIONS
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(_run, condition): condition for condition in conditions}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(
                row["opposite_width_fraction"], row["opposite_amplitude_ratio"],
                row["opposite_duration_s"], row["switch_success"], flush=True,
            )
    rows.sort(key=lambda row: (
        row["opposite_width_fraction"], row["opposite_amplitude_ratio"],
        row["opposite_duration_s"],
    ))
    with (OUTPUT / "screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    figure, axes = plt.subplots(1, len(WIDTHS), figsize=(10.2, 2.8), constrained_layout=True)
    for axis, width in zip(axes, WIDTHS):
        selected = [row for row in rows if row["opposite_width_fraction"] == width]
        matrix = np.zeros((len(RATIOS), len(DURATIONS)))
        for row in selected:
            i = RATIOS.index(row["opposite_amplitude_ratio"])
            j = DURATIONS.index(row["opposite_duration_s"])
            matrix[i, j] = bool(row["switch_success"])
        axis.imshow(matrix, origin="lower", cmap="Greys", vmin=0, vmax=1, aspect="auto")
        axis.set_xticks(range(len(DURATIONS)), [f"{value:g}" for value in DURATIONS], rotation=45)
        axis.set_yticks(range(len(RATIOS)), [f"{value:g}" for value in RATIOS])
        axis.set_title(f"Opposite width={width:g}L")
        axis.set_xlabel("Duration (s)")
        axis.set_ylabel(r"Amplitude $\alpha_{cue}/\alpha_0$")
    figure.suptitle(
        "Autonomous switching after withdrawal at L=15 µm\n"
        "white=switch retained ≥1,000 s; black=no switch",
        fontsize=10,
    )
    figure.savefig(OUTPUT / "opposite_pulse_breaking_screen.png", dpi=300)
    figure.savefig(OUTPUT / "opposite_pulse_breaking_screen.pdf")
    plt.close(figure)
    passing = [row for row in rows if row["switch_success"]]
    print(f"PASSING={len(passing)}")
    if passing:
        print(passing[:10])


if __name__ == "__main__":
    main()
