#!/usr/bin/env python3
"""Find an L=15 regime that is both strongly polarizable and pulse-switchable."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import DeBelly, DeBellyParameters, FinitePulseSequence, simulate  # noqa: E402


LENGTH = 15.0
CELLS = 75
STEP = 0.5
REVERSE_START = 600.0
PARAMETERS = (
    (1.5, 1.5), (1.8, 1.8), (2.0, 1.9), (2.0, 2.0),
    (2.0, 2.1), (2.2, 2.2), (3.0, 3.5),
)


def _orientation(values, x):
    left = np.mean(values[:, x <= 0.2], axis=1)
    right = np.mean(values[:, x >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _contrast(values):
    high, low = np.max(values, axis=1), np.min(values, axis=1)
    return (high - low) / np.maximum(np.abs(high) + np.abs(low), 1e-12)


def _first_sustained(time, condition, *, after, duration=100.0):
    start = int(np.searchsorted(time, after, side="left"))
    for index in range(start, len(time)):
        stop = int(np.searchsorted(time, time[index] + duration, side="left"))
        if stop < len(time) and np.all(condition[index : stop + 1]):
            return float(time[index])
    return None


def _pulse(side, start, duration, amplitude, width):
    return {
        "amplitude": amplitude,
        "start": start,
        "stop": start + duration,
        "taper": duration / 2.0,
        "width_fraction": width,
        "sides": [side],
    }


def _run(spec):
    alpha0, beta0, opposite_ratio, opposite_width, opposite_duration = spec
    model = DeBelly(DeBellyParameters(
        rac_basal_activation=alpha0,
        rho_basal_activation=beta0,
        rac_mechanical_feedback=0.3 * alpha0,
        rho_mechanical_feedback=0.3 * beta0,
        assay_mode="reflection_symmetric_extension",
        initialization_mode="algebraic_fixed_point",
        input_sampling_mode="cell_centred",
    ))
    stop = REVERSE_START + opposite_duration
    trajectory = simulate(
        model,
        protocol=FinitePulseSequence(pulses=(
            _pulse("left", 0.0, 160.0, 2.0 * alpha0, 0.2),
            _pulse(
                "right", REVERSE_START, opposite_duration,
                opposite_ratio * alpha0, opposite_width,
            ),
        )),
        length=LENGTH, cells=CELLS, end_time=stop + 800.0,
        output_interval=5.0, step=STEP,
    )
    x = trajectory.x / LENGTH
    rac = _orientation(trajectory.fields["rac"], x)
    rho = _orientation(trajectory.fields["rho"], x)
    contrast = _contrast(trajectory.fields["rac"])
    pre = int(np.argmin(np.abs(trajectory.time - REVERSE_START)))
    pre_strong = bool(-rac[pre] >= 0.25 and rho[pre] >= 0.25 and contrast[pre] >= 0.3)
    target = (rac >= 0.25) & (rho <= -0.25) & (contrast >= 0.3)
    reversal = _first_sustained(trajectory.time, target, after=stop, duration=100.0)
    tail = trajectory.time >= trajectory.time[-1] - 300.0
    persistent = bool(np.all(target[tail]))
    return {
        "alpha0": alpha0,
        "beta0": beta0,
        "opposite_ratio": opposite_ratio,
        "opposite_width_fraction": opposite_width,
        "opposite_duration_s": opposite_duration,
        "pre_strong": pre_strong,
        "pre_rac_directed": float(-rac[pre]),
        "pre_rho_directed": float(rho[pre]),
        "pre_rac_contrast": float(contrast[pre]),
        "reversal_time_s": reversal,
        "latency_after_start_s": None if reversal is None else reversal - REVERSE_START,
        "latency_after_withdrawal_s": None if reversal is None else reversal - stop,
        "persistent_reversed_tail": persistent,
        "final_rac_directed": float(rac[-1]),
        "final_rho_directed": float(-rho[-1]),
        "final_rac_contrast": float(contrast[-1]),
        "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
    }


def main():
    output = ROOT / "results/tmp/trial/debelly_l15_protocol_design"
    output.mkdir(parents=True, exist_ok=True)
    specs = [
        (*parameters, ratio, width, duration)
        for parameters in PARAMETERS
        for ratio in (2.0, 4.0, 8.0)
        for width in (0.2, 0.4)
        for duration in (160.0, 320.0, 640.0)
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(_run, spec) for spec in specs]
        for index, future in enumerate(as_completed(futures), start=1):
            rows.append(future.result())
            if index % 21 == 0:
                print(f"completed {index}/{len(futures)}", flush=True)
    rows.sort(key=lambda row: (
        not (row["pre_strong"] and row["persistent_reversed_tail"]),
        row["reversal_time_s"] is None,
        np.inf if row["reversal_time_s"] is None else row["reversal_time_s"],
        row["opposite_ratio"] * row["opposite_width_fraction"] * row["opposite_duration_s"],
    ))
    with (output / "switchability_screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows[:30]:
        print(row)


if __name__ == "__main__":
    main()
