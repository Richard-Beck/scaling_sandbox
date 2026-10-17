#!/usr/bin/env python3
"""Design a compact L=15 induction/reversal assay for the responsive regime."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly, DeBellyParameters, FinitePulseSequence, LocalizedPulse, simulate,
)


LENGTH = 15.0
CELLS = 75
STEP = 0.5
ALPHA0 = 2.0
BETA0 = 2.0
REVERSE_START = 1_000.0


def _model():
    return DeBelly(DeBellyParameters(
        rac_basal_activation=ALPHA0,
        rho_basal_activation=BETA0,
        rac_mechanical_feedback=0.3 * ALPHA0,
        rho_mechanical_feedback=0.3 * BETA0,
        assay_mode="reflection_symmetric_extension",
        initialization_mode="algebraic_fixed_point",
        input_sampling_mode="cell_centred",
    ))


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


def _pulse(side, start, duration, ratio, width):
    return {
        "amplitude": ratio * ALPHA0,
        "start": start,
        "stop": start + duration,
        "taper": duration / 2.0,
        "width_fraction": width,
        "sides": [side],
    }


def _induction(spec):
    ratio, width, duration = spec
    trajectory = simulate(
        _model(),
        protocol=LocalizedPulse(
            amplitude=ratio * ALPHA0, start=0.0, stop=duration,
            taper=duration / 2.0, width_fraction=width, sides=("left",),
        ),
        length=LENGTH, cells=CELLS, end_time=1_000.0,
        output_interval=5.0, step=STEP,
    )
    x = trajectory.x / LENGTH
    rac = -_orientation(trajectory.fields["rac"], x)
    rho = _orientation(trajectory.fields["rho"], x)
    contrast = _contrast(trajectory.fields["rac"])
    strong = (rac >= 0.3) & (rho >= 0.3) & (contrast >= 0.3)
    establishment = _first_sustained(
        trajectory.time, strong, after=0.0, duration=100.0
    )
    return {
        "induction_ratio": ratio,
        "induction_width_fraction": width,
        "induction_duration_s": duration,
        "strong_establishment_s": establishment,
        "rac_orientation_1000": float(rac[-1]),
        "rho_orientation_1000": float(rho[-1]),
        "rac_contrast_1000": float(contrast[-1]),
    }


def _reversal(spec):
    induction_ratio, induction_width, induction_duration, ratio, width, duration = spec
    stop = REVERSE_START + duration
    trajectory = simulate(
        _model(),
        protocol=FinitePulseSequence(pulses=(
            _pulse("left", 0.0, induction_duration, induction_ratio, induction_width),
            _pulse("right", REVERSE_START, duration, ratio, width),
        )),
        length=LENGTH, cells=CELLS, end_time=stop + 1_200.0,
        output_interval=5.0, step=STEP,
    )
    x = trajectory.x / LENGTH
    rac = _orientation(trajectory.fields["rac"], x)
    rho = _orientation(trajectory.fields["rho"], x)
    contrast = _contrast(trajectory.fields["rac"])
    pre = int(np.argmin(np.abs(trajectory.time - REVERSE_START)))
    reversed_strong = (rac >= 0.3) & (rho <= -0.3) & (contrast >= 0.3)
    reversal = _first_sustained(
        trajectory.time, reversed_strong, after=stop, duration=100.0
    )
    tail = trajectory.time >= trajectory.time[-1] - 500.0
    persistent = bool(np.all(reversed_strong[tail]))
    return {
        "induction_ratio": induction_ratio,
        "induction_width_fraction": induction_width,
        "induction_duration_s": induction_duration,
        "opposite_ratio": ratio,
        "opposite_width_fraction": width,
        "opposite_duration_s": duration,
        "pre_rac_directed": float(-rac[pre]),
        "pre_rho_directed": float(rho[pre]),
        "pre_rac_contrast": float(contrast[pre]),
        "strong_reversal_time_s": reversal,
        "latency_after_opposite_start_s": (
            None if reversal is None else reversal - REVERSE_START
        ),
        "latency_after_withdrawal_s": None if reversal is None else reversal - stop,
        "persistent_reversed_tail": persistent,
        "final_rac_directed": float(rac[-1]),
        "final_rho_directed": float(-rho[-1]),
        "final_rac_contrast": float(contrast[-1]),
        "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
    }


def _parallel(function, specs):
    rows = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(function, spec) for spec in specs]
        for future in as_completed(futures):
            rows.append(future.result())
    return rows


def _write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    output = ROOT / "results/tmp/trial/debelly_l15_protocol_design"
    output.mkdir(parents=True, exist_ok=True)
    induction_specs = [
        (ratio, width, duration)
        for ratio in (1.0, 2.0, 4.0)
        for width in (0.1, 0.2, 0.3)
        for duration in (40.0, 80.0, 160.0)
    ]
    induction = _parallel(_induction, induction_specs)
    induction.sort(key=lambda row: (
        row["strong_establishment_s"] is None,
        np.inf if row["strong_establishment_s"] is None else row["strong_establishment_s"],
        -row["rac_orientation_1000"],
    ))
    _write(output / "induction_screen.csv", induction)
    print("Top induction conditions")
    for row in induction[:10]:
        print(row)

    # Use the fastest condition that remains compact rather than automatically
    # selecting the widest/strongest cue.
    selected = next(
        row for row in induction
        if row["induction_ratio"] <= 2.0
        and row["induction_width_fraction"] <= 0.2
        and row["strong_establishment_s"] is not None
    )
    base = (
        selected["induction_ratio"], selected["induction_width_fraction"],
        selected["induction_duration_s"],
    )
    reversal_specs = [
        (*base, ratio, width, duration)
        for ratio in (1.0, 2.0, 4.0)
        for width in (0.1, 0.2, 0.3, 0.4)
        for duration in (80.0, 160.0, 320.0)
    ]
    reversal = _parallel(_reversal, reversal_specs)
    reversal.sort(key=lambda row: (
        not row["persistent_reversed_tail"],
        row["strong_reversal_time_s"] is None,
        np.inf if row["strong_reversal_time_s"] is None else row["strong_reversal_time_s"],
        row["opposite_ratio"] * row["opposite_width_fraction"] * row["opposite_duration_s"],
    ))
    _write(output / "reversal_screen.csv", reversal)
    print("Selected induction", selected)
    print("Top reversal conditions")
    for row in reversal[:15]:
        print(row)


if __name__ == "__main__":
    main()
