#!/usr/bin/env python3
"""Small L=15 design screen for strong two-pole establishment and resolution."""

from __future__ import annotations

import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import DeBelly, DeBellyParameters, LocalizedPulse, simulate  # noqa: E402
from polarity1d.metrics import profile_domain_count  # noqa: E402


LENGTH = 15.0
CELLS = 75
STEP = 0.5
PARAMETERS = (
    (1.8, 1.8),
    (2.0, 1.9),
    (2.0, 2.0),
    (2.0, 2.1),
    (2.2, 2.2),
    (3.0, 3.5),
)
CUES = tuple(
    (ratio, width, duration, asymmetry)
    for ratio in (1.0, 2.0)
    for width in (0.1, 0.2)
    for duration in (40.0, 80.0)
    for asymmetry in (0.05, 0.15)
)


def _first_sustained(time, condition, *, after, duration=100.0):
    start = int(np.searchsorted(time, after, side="left"))
    for index in range(start, len(time)):
        stop = int(np.searchsorted(time, time[index] + duration, side="left"))
        if stop < len(time) and np.all(condition[index : stop + 1]):
            return float(time[index])
    return None


def _run(spec):
    alpha0, beta0, ratio, width, duration, asymmetry = spec
    parameters = DeBellyParameters(
        rac_basal_activation=alpha0,
        rho_basal_activation=beta0,
        rac_mechanical_feedback=0.3 * alpha0,
        rho_mechanical_feedback=0.3 * beta0,
        assay_mode="reflection_symmetric_extension",
        initialization_mode="algebraic_fixed_point",
        input_sampling_mode="cell_centred",
    )
    trajectory = simulate(
        DeBelly(parameters),
        protocol=LocalizedPulse(
            amplitude=ratio * alpha0,
            start=0.0,
            stop=duration,
            taper=duration / 2.0,
            width_fraction=width,
            sides=("left", "right"),
            asymmetry=asymmetry,
        ),
        length=LENGTH,
        cells=CELLS,
        end_time=3_000.0,
        output_interval=5.0,
        step=STEP,
    )
    time = trajectory.time
    x = trajectory.x / LENGTH
    rac = trajectory.fields["rac"]
    left = np.mean(rac[:, x <= 0.15], axis=1)
    right = np.mean(rac[:, x >= 0.85], axis=1)
    centre = np.mean(rac[:, (x >= 0.4) & (x <= 0.6)], axis=1)
    left_strength = (left - centre) / np.maximum(np.abs(left) + np.abs(centre), 1e-12)
    right_strength = (right - centre) / np.maximum(np.abs(right) + np.abs(centre), 1e-12)
    two_pole_strength = np.minimum(left_strength, right_strength)
    domains = np.asarray([profile_domain_count(profile) for profile in rac])
    contrast = (np.max(rac, axis=1) - np.min(rac, axis=1)) / np.maximum(
        np.abs(np.max(rac, axis=1)) + np.abs(np.min(rac, axis=1)), 1e-12
    )
    orientation = (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)

    establishment_window = time <= duration + 150.0
    qualifying_two = establishment_window & (domains >= 2) & (two_pole_strength >= 0.05)
    two_time = float(time[np.flatnonzero(qualifying_two)[0]]) if np.any(qualifying_two) else None
    peak_two = float(np.max(two_pole_strength[establishment_window]))
    resolved = (domains == 1) & (contrast >= 0.25) & (np.abs(orientation) >= 0.2)
    resolution_time = None if two_time is None else _first_sustained(
        time, resolved, after=max(duration, two_time), duration=100.0
    )
    tail = time >= time[-1] - 500.0
    persistent = bool(np.all(resolved[tail]))
    score = (
        (0.0 if two_time is None else 3.0)
        + (0.0 if resolution_time is None else 3.0)
        + (2.0 if persistent else 0.0)
        + min(peak_two, 0.5)
        + min(float(contrast[-1]), 1.0)
        + min(abs(float(orientation[-1])), 1.0)
        - (0.0 if resolution_time is None else resolution_time / 3_000.0)
    )
    return {
        "alpha0": alpha0,
        "beta0": beta0,
        "cue_ratio": ratio,
        "cue_width_fraction": width,
        "cue_duration_s": duration,
        "cue_asymmetry": asymmetry,
        "two_pole_time_s": two_time,
        "peak_two_pole_strength": peak_two,
        "resolution_time_s": resolution_time,
        "resolution_latency_after_cue_s": (
            None if resolution_time is None else resolution_time - duration
        ),
        "persistent_resolved_tail": persistent,
        "final_domain_count": int(domains[-1]),
        "final_contrast": float(contrast[-1]),
        "final_orientation": float(orientation[-1]),
        "mca_conservation_error": trajectory.metadata["mca_conservation_error"],
        "score": score,
    }


def main():
    output = ROOT / "results/tmp/trial/debelly_l15_protocol_design"
    output.mkdir(parents=True, exist_ok=True)
    specs = [(*parameters, *cue) for parameters in PARAMETERS for cue in CUES]
    rows = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(_run, spec) for spec in specs]
        for index, future in enumerate(as_completed(futures), start=1):
            rows.append(future.result())
            if index % 16 == 0:
                print(f"completed {index}/{len(futures)}", flush=True)
    rows.sort(key=lambda row: row["score"], reverse=True)
    with (output / "competition_screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows[:20]:
        print(row)


if __name__ == "__main__":
    main()
