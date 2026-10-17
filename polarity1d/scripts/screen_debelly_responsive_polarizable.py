#!/usr/bin/env python3
"""Select a responsive/polarizable de Belly regime before length sweeps."""

from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    DeBelly,
    DeBellyParameters,
    LocalizedPulse,
    NoStimulus,
    ReversingGradient,
    simulate,
)


LENGTH = 15.0
CELLS = 75
STEP = 0.5
VALUES = tuple(float(value) for value in np.arange(1.0, 5.51, 0.5))
ORIENTATION_THRESHOLD = 0.05


def _orientation(trajectory, field: str) -> np.ndarray:
    relative = (trajectory.x - trajectory.x[0]) / (
        trajectory.x[-1] - trajectory.x[0]
    )
    values = trajectory.fields[field]
    left = np.mean(values[:, relative <= 0.2], axis=1)
    right = np.mean(values[:, relative >= 0.8], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _contrast(values: np.ndarray) -> np.ndarray:
    maxima = np.max(values, axis=1)
    minima = np.min(values, axis=1)
    return (maxima - minima) / np.maximum(np.abs(maxima) + np.abs(minima), 1e-12)


def _at(trajectory, values: np.ndarray, time: float) -> float:
    return float(values[int(np.argmin(np.abs(trajectory.time - time)))])


def _first_sustained(
    time: np.ndarray, condition: np.ndarray, *, after: float, duration: float
) -> float | None:
    start = int(np.searchsorted(time, after, side="left"))
    for index in range(start, len(time)):
        stop = int(np.searchsorted(time, time[index] + duration, side="left"))
        if stop < len(time) and np.all(condition[index : stop + 1]):
            return float(time[index])
    return None


def _parameters(alpha0: float, beta0: float, *, mechanics: bool) -> DeBellyParameters:
    common = dict(
        rac_basal_activation=alpha0,
        rho_basal_activation=beta0,
        rac_mechanical_feedback=0.3 * alpha0,
        rho_mechanical_feedback=0.3 * beta0,
        assay_mode="reflection_symmetric_extension",
        initialization_mode="algebraic_fixed_point",
        input_sampling_mode="cell_centred",
    )
    if not mechanics:
        common.update(
            protrusion_speed=0.0,
            rac_mechanical_feedback=0.0,
            rho_mechanical_feedback=0.0,
            active_stress=0.0,
        )
    return DeBellyParameters(**common)


def _initial_fields(model: DeBelly) -> dict[str, np.ndarray]:
    seed = simulate(
        model, protocol=NoStimulus(), length=LENGTH, cells=CELLS,
        end_time=STEP, output_interval=STEP, step=STEP,
    )
    relative = 2.0 * seed.x / LENGTH - 1.0
    # One-percent opposed perturbation tests local stability without using a
    # perturbation large enough to cross a nearby polarized basin boundary.
    return {
        "rac": seed.fields["rac"][0] * (1.0 - 0.01 * relative),
        "rho": seed.fields["rho"][0] * (1.0 + 0.01 * relative),
    }


def _pulse_metrics(alpha0: float, beta0: float, *, mechanics: bool) -> dict:
    model = DeBelly(_parameters(alpha0, beta0, mechanics=mechanics))
    pulse = LocalizedPulse(
        amplitude=alpha0, start=0.0, stop=40.0, taper=20.0,
        width_fraction=0.1, sides=("left",),
    )
    trajectory = simulate(
        model, protocol=pulse, length=LENGTH, cells=CELLS,
        end_time=5_000.0, output_interval=5.0, step=STEP,
    )
    rac = _orientation(trajectory, "rac")
    rho = _orientation(trajectory, "rho")
    directed_rac = -rac
    directed_rho = rho
    polarized = (directed_rac >= ORIENTATION_THRESHOLD) & (
        directed_rho >= ORIENTATION_THRESHOLD
    )
    loss = _first_sustained(
        trajectory.time,
        ~polarized,
        after=40.0,
        duration=100.0,
    )
    return {
        "peak_rac_directed": float(np.max(directed_rac)),
        "peak_rho_directed": float(np.max(directed_rho)),
        "rac_500": _at(trajectory, directed_rac, 500.0),
        "rho_500": _at(trajectory, directed_rho, 500.0),
        "rac_2000": _at(trajectory, directed_rac, 2_000.0),
        "rho_2000": _at(trajectory, directed_rho, 2_000.0),
        "rac_5000": float(directed_rac[-1]),
        "rho_5000": float(directed_rho[-1]),
        "loss_time": loss,
        "mass_error": trajectory.metadata["mca_conservation_error"],
    }


def _stage1(pair: tuple[float, float]) -> dict:
    alpha0, beta0 = pair
    full_model = DeBelly(_parameters(alpha0, beta0, mechanics=True))
    no_cue = simulate(
        full_model, protocol=NoStimulus(), length=LENGTH, cells=CELLS,
        end_time=3_000.0, output_interval=5.0, step=STEP,
        initial_fields=_initial_fields(full_model),
    )
    rac0 = _orientation(no_cue, "rac")
    rho0 = _orientation(no_cue, "rho")
    final_contrast = max(
        float(_contrast(no_cue.fields["rac"])[-1]),
        float(_contrast(no_cue.fields["rho"])[-1]),
    )
    tail = no_cue.time >= 2_000.0
    no_cue_pass = bool(
        max(abs(rac0[-1]), abs(rho0[-1])) < 0.01
        and max(np.max(np.abs(rac0[tail])), np.max(np.abs(rho0[tail]))) < 0.02
        and final_contrast < 0.05
    )

    full = _pulse_metrics(alpha0, beta0, mechanics=True)
    local = _pulse_metrics(alpha0, beta0, mechanics=False)
    responsive = bool(
        full["peak_rac_directed"] >= ORIENTATION_THRESHOLD
        and full["peak_rho_directed"] >= ORIENTATION_THRESHOLD
    )
    persistent = bool(
        min(full["rac_2000"], full["rho_2000"], full["rac_5000"], full["rho_5000"])
        >= ORIENTATION_THRESHOLD
    )
    mechanics_required = bool(
        min(local["rac_2000"], local["rho_2000"]) < ORIENTATION_THRESHOLD
        or min(local["rac_5000"], local["rho_5000"])
        < 0.5 * min(full["rac_5000"], full["rho_5000"])
    )
    return {
        "alpha0": alpha0,
        "beta0": beta0,
        "no_cue_rac_final": float(rac0[-1]),
        "no_cue_rho_final": float(rho0[-1]),
        "no_cue_final_max_contrast": final_contrast,
        "no_cue_pass": no_cue_pass,
        **{f"full_{key}": value for key, value in full.items()},
        **{f"local_{key}": value for key, value in local.items()},
        "responsive_pass": responsive,
        "persistent_pass": persistent,
        "mechanics_required_pass": mechanics_required,
        "stage1_pass": no_cue_pass and responsive and persistent and mechanics_required,
    }


def _reversal_metrics(pair: tuple[float, float]) -> dict:
    alpha0, beta0 = pair
    reversal_time = 2_000.0
    trajectory = simulate(
        DeBelly(_parameters(alpha0, beta0, mechanics=True)),
        protocol=ReversingGradient(
            slope=alpha0 / LENGTH, reversal_time=reversal_time, transition=20.0
        ),
        length=LENGTH, cells=CELLS, end_time=6_000.0,
        output_interval=5.0, step=STEP,
    )
    rac = _orientation(trajectory, "rac")
    rho = _orientation(trajectory, "rho")
    preconditioned = bool(
        -_at(trajectory, rac, 1_900.0) >= ORIENTATION_THRESHOLD
        and _at(trajectory, rho, 1_900.0) >= ORIENTATION_THRESHOLD
    )
    reversed_state = (rac >= ORIENTATION_THRESHOLD) & (
        rho <= -ORIENTATION_THRESHOLD
    )
    reversal = _first_sustained(
        trajectory.time, reversed_state, after=reversal_time, duration=100.0
    )
    reversible = bool(
        preconditioned
        and reversal is not None
        and reversal - reversal_time <= 2_000.0
        and rac[-1] >= ORIENTATION_THRESHOLD
        and rho[-1] <= -ORIENTATION_THRESHOLD
    )
    return {
        "alpha0": alpha0,
        "beta0": beta0,
        "reversal_pre_rac_directed": -_at(trajectory, rac, 1_900.0),
        "reversal_pre_rho_directed": _at(trajectory, rho, 1_900.0),
        "reversal_time_s": reversal,
        "reversal_latency_s": None if reversal is None else reversal - reversal_time,
        "reversal_final_rac_directed": float(rac[-1]),
        "reversal_final_rho_directed": float(-rho[-1]),
        "reversal_pass": reversible,
        "reversal_mass_error": trajectory.metadata["mca_conservation_error"],
    }


def _long_memory_metrics(pair: tuple[float, float]) -> dict:
    """Resolve slow instability versus autonomous cue-induced memory."""
    alpha0, beta0 = pair
    model = DeBelly(_parameters(alpha0, beta0, mechanics=True))
    no_cue = simulate(
        model, protocol=NoStimulus(), length=LENGTH, cells=CELLS,
        end_time=20_000.0, output_interval=20.0, step=STEP,
        initial_fields=_initial_fields(model),
    )
    pulse = simulate(
        model,
        protocol=LocalizedPulse(
            amplitude=alpha0, start=0.0, stop=40.0, taper=20.0,
            width_fraction=0.1, sides=("left",),
        ),
        length=LENGTH, cells=CELLS, end_time=20_000.0,
        output_interval=20.0, step=STEP,
    )
    no_rac, no_rho = _orientation(no_cue, "rac"), _orientation(no_cue, "rho")
    pulse_rac = -_orientation(pulse, "rac")
    pulse_rho = _orientation(pulse, "rho")
    no_at_10k = max(
        abs(_at(no_cue, no_rac, 10_000.0)),
        abs(_at(no_cue, no_rho, 10_000.0)),
    )
    no_at_20k = max(abs(no_rac[-1]), abs(no_rho[-1]))
    no_cue_stable = bool(
        no_at_20k < 0.02
        and no_at_20k <= 0.9 * max(no_at_10k, 1e-12)
    )
    pulse_persistent = bool(
        min(
            np.min(pulse_rac[pulse.time >= 15_000.0]),
            np.min(pulse_rho[pulse.time >= 15_000.0]),
        ) >= ORIENTATION_THRESHOLD
    )
    return {
        "alpha0": alpha0,
        "beta0": beta0,
        "long_no_cue_rac_20000": float(no_rac[-1]),
        "long_no_cue_rho_20000": float(no_rho[-1]),
        "long_no_cue_maxabs_5000": max(
            abs(_at(no_cue, no_rac, 5_000.0)),
            abs(_at(no_cue, no_rho, 5_000.0)),
        ),
        "long_no_cue_maxabs_10000": no_at_10k,
        "long_no_cue_maxabs_15000": max(
            abs(_at(no_cue, no_rac, 15_000.0)),
            abs(_at(no_cue, no_rho, 15_000.0)),
        ),
        "long_no_cue_maxabs_20000": no_at_20k,
        "long_pulse_rac_20000": float(pulse_rac[-1]),
        "long_pulse_rho_20000": float(pulse_rho[-1]),
        "long_no_cue_stable": no_cue_stable,
        "long_pulse_persistent": pulse_persistent,
        "long_memory_pass": no_cue_stable and pulse_persistent,
    }


def _parallel(function, items, workers: int) -> list[dict]:
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(function, item): item for item in items}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(
                f"alpha0={row['alpha0']:.1f} beta0={row['beta0']:.1f}",
                flush=True,
            )
    return sorted(rows, key=lambda row: (row["alpha0"], row["beta0"]))


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _draw(rows: list[dict], destination: Path) -> None:
    fields = (
        ("no_cue_pass", "No cue: unpolarized"),
        ("persistent_pass", "Brief cue: persistent"),
        ("mechanics_required_pass", "Mechanics-dependent"),
        ("reversal_pass", "Opposing cue: reversible"),
        ("all_pass", "All independent gates"),
    )
    figure, axes = plt.subplots(1, 5, figsize=(13.2, 2.9), constrained_layout=True)
    extent = (VALUES[0] - 0.25, VALUES[-1] + 0.25, VALUES[0] - 0.25, VALUES[-1] + 0.25)
    for axis, (field, title) in zip(axes, fields):
        values = np.zeros((len(VALUES), len(VALUES)))
        for row in rows:
            i = VALUES.index(row["alpha0"])
            j = VALUES.index(row["beta0"])
            values[i, j] = bool(row.get(field, False))
        axis.imshow(
            values, origin="lower", extent=extent, interpolation="nearest",
            cmap="Greys", vmin=0, vmax=1, aspect="equal",
        )
        axis.set_title(title, fontsize=9)
        axis.set_xticks(VALUES[::2])
        axis.set_yticks(VALUES[::2])
        axis.set_xlabel(r"$\beta_0$")
        axis.set_ylabel(r"$\alpha_0$")
    figure.suptitle(
        r"Responsive/polarizable screen at $L=15\,\mu$m; "
        r"$\alpha=0.3\alpha_0$, $\beta=0.3\beta_0$",
        fontsize=10,
    )
    figure.savefig(destination, dpi=300)
    figure.savefig(destination.with_suffix(".pdf"))
    plt.close(figure)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "output_dir", nargs="?", type=Path,
        default=ROOT / "results" / "tmp" / "trial" / "debelly_responsive_polarizable_screen",
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs = [(alpha0, beta0) for alpha0 in VALUES for beta0 in VALUES]
    stage1 = _parallel(_stage1, pairs, args.workers)
    _write_csv(args.output_dir / "stage1.csv", stage1)

    stage1_pairs = [
        (row["alpha0"], row["beta0"]) for row in stage1 if row["stage1_pass"]
    ]
    reversal = _parallel(_reversal_metrics, stage1_pairs, args.workers) if stage1_pairs else []
    reversal_by_pair = {(row["alpha0"], row["beta0"]): row for row in reversal}
    combined = []
    for row in stage1:
        reverse = reversal_by_pair.get((row["alpha0"], row["beta0"]), {})
        merged = {**row, **reverse}
        merged["reversal_pass"] = bool(reverse.get("reversal_pass", False))
        merged["all_pass"] = bool(row["stage1_pass"] and merged["reversal_pass"])
        combined.append(merged)
    all_fields = []
    for row in combined:
        for field in row:
            if field not in all_fields:
                all_fields.append(field)
    with (args.output_dir / "screen.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(combined)
    _draw(combined, args.output_dir / "responsive_polarizable_screen.png")

    passing = [row for row in combined if row["all_pass"]]
    summary = {
        "criteria_are_independent_of_length_dependence": True,
        "length_um": LENGTH,
        "values": list(VALUES),
        "candidate_count": len(passing),
        "candidates": [
            {key: row.get(key) for key in (
                "alpha0", "beta0", "full_rac_5000", "full_rho_5000",
                "local_rac_5000", "local_rho_5000", "reversal_latency_s",
            )}
            for row in passing
        ],
        "selection": None,
        "thresholds": {
            "orientation": ORIENTATION_THRESHOLD,
            "no_cue_tail_max_orientation": 0.02,
            "no_cue_final_orientation": 0.01,
            "no_cue_final_contrast": 0.05,
            "persistence_times_s": [2_000.0, 5_000.0],
            "reversal_max_latency_s": 2_000.0,
        },
    }
    (args.output_dir / "selection.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"passing candidates: {len(passing)}")


if __name__ == "__main__":
    main()
