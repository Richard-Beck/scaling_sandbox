#!/usr/bin/env python3
"""Generate causal and weak-cue diagnostics for the symmetric de Belly extension."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import DeBelly, DeBellyParameters, LocalizedPulse, simulate  # noqa: E402


OUTPUT = ROOT / "results" / "validation" / "debelly_2026_mechanism_activity"
LENGTH = 15.0
CELLS = 75
ALPHA0 = 3.5
BETA0 = 4.0


def _model(*, target: str = "rac", **overrides) -> DeBelly:
    values = {
        "assay_mode": "reflection_symmetric_extension",
        "initialization_mode": "algebraic_fixed_point",
        "input_sampling_mode": "cell_centred",
        "opto_target": target,
        "rac_basal_activation": ALPHA0,
        "rho_basal_activation": BETA0,
        "rac_mechanical_feedback": 0.3 * ALPHA0,
        "rho_mechanical_feedback": 0.3 * BETA0,
    }
    values.update(overrides)
    return DeBelly(DeBellyParameters(**values))


def _pulse(amplitude: float = ALPHA0) -> LocalizedPulse:
    return LocalizedPulse(
        amplitude=amplitude, start=0.0, stop=40.0, taper=20.0,
        width_fraction=0.1, sides=("left",),
    )


def _orientation(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    left = np.mean(values[:, x <= 0.2 * LENGTH], axis=1)
    right = np.mean(values[:, x >= 0.8 * LENGTH], axis=1)
    return (right - left) / np.maximum(np.abs(right) + np.abs(left), 1e-12)


def _trajectory_rows(case: str, condition: str, trajectory) -> list[dict]:
    proximal = trajectory.x <= 0.2 * LENGTH
    distal = trajectory.x >= 0.8 * LENGTH
    parameters = trajectory.metadata["parameters"]
    rows = []
    for index, time in enumerate(trajectory.time):
        row = {"case": case, "condition": condition, "time_s": float(time)}
        for field in ("rac", "rho", "tension", "mca_bound", "mca_unbound"):
            values = trajectory.fields[field][index]
            row[f"proximal_{field}"] = float(np.mean(values[proximal]))
            row[f"distal_{field}"] = float(np.mean(values[distal]))
        row.update(
            {
                "maximum_delta_beta_tension": float(
                    np.max(trajectory.fields["delta_beta_tension"][index])
                ),
                "maximum_delta_alpha_mca": float(
                    np.max(trajectory.fields["delta_alpha_mca"][index])
                ),
                "fraction_tension_above_threshold": float(
                    np.mean(
                        trajectory.fields["tension"][index]
                        > parameters["tension_rho_threshold"]
                    )
                ),
                "fraction_mca_below_threshold": float(
                    np.mean(
                        trajectory.fields["mca_bound"][index]
                        < parameters["mca_rac_threshold"]
                    )
                ),
                "boundary_separation_um": float(
                    trajectory.fields["boundary_separation"][index, 0]
                ),
                "mean_strain": float(trajectory.fields["mean_strain"][index, 0]),
                "tension_range": float(trajectory.fields["tension_range"][index, 0]),
            }
        )
        rows.append(row)
    return rows


def _write_csv(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    summaries: dict[str, dict] = {}

    causal_specs = {
        "rac_pulse_full": _model(target="rac"),
        "rac_pulse_beta_zero": _model(target="rac", rho_mechanical_feedback=0.0),
        "rho_pulse_full": _model(target="rho"),
        "rho_pulse_alpha_zero": _model(target="rho", rac_mechanical_feedback=0.0),
    }
    trajectories = {}
    for name, model in causal_specs.items():
        trajectory = simulate(
            model, protocol=_pulse(), length=LENGTH, cells=CELLS,
            end_time=800.0, output_interval=1.0, step=0.5,
        )
        trajectory.save(OUTPUT / f"{name}.npz")
        trajectories[name] = trajectory
        rows.extend(_trajectory_rows("causal_pulse", name, trajectory))

    full = trajectories["rac_pulse_full"]
    beta_zero = trajectories["rac_pulse_beta_zero"]
    distal = full.x >= 0.8 * LENGTH
    full_rho = np.mean(full.fields["rho"][:, distal], axis=1)
    beta_zero_rho = np.mean(beta_zero.fields["rho"][:, distal], axis=1)
    summaries["left_rac_to_distal_rho"] = {
        "full_peak_change": float(np.max(full_rho - full_rho[0])),
        "beta_zero_peak_change": float(np.max(beta_zero_rho - beta_zero_rho[0])),
        "maximum_full_minus_beta_zero": float(np.max(full_rho - beta_zero_rho)),
    }

    full = trajectories["rho_pulse_full"]
    alpha_zero = trajectories["rho_pulse_alpha_zero"]
    full_rac = np.mean(full.fields["rac"][:, distal], axis=1)
    alpha_zero_rac = np.mean(alpha_zero.fields["rac"][:, distal], axis=1)
    summaries["left_rho_to_distal_rac"] = {
        "full_peak_change": float(np.max(full_rac - full_rac[0])),
        "alpha_zero_peak_change": float(np.max(alpha_zero_rac - alpha_zero_rac[0])),
        "maximum_full_minus_alpha_zero": float(np.max(full_rac - alpha_zero_rac)),
    }

    weak_summary = []
    for amplitude_ratio in (0.1, 0.25, 0.5, 1.0):
        for condition, model in (
            ("full", _model(target="rac")),
            (
                "local_inhibition",
                _model(
                    target="rac", protrusion_speed=0.0,
                    rac_mechanical_feedback=0.0, rho_mechanical_feedback=0.0,
                    active_stress=0.0,
                ),
            ),
        ):
            trajectory = simulate(
                model, protocol=_pulse(amplitude_ratio * ALPHA0), length=LENGTH, cells=CELLS,
                end_time=5_000.0, output_interval=5.0, step=0.5,
            )
            name = f"weak_rac_{condition}_a{amplitude_ratio:g}"
            trajectory.save(OUTPUT / f"{name}.npz")
            rows.extend(_trajectory_rows("weak_cue_withdrawal", name, trajectory))
            orientation = _orientation(trajectory.fields["rac"], trajectory.x)
            item = {
                "amplitude": amplitude_ratio, "condition": condition,
                "maximum_absolute_orientation": float(np.max(np.abs(orientation))),
            }
            for time in (40.0, 200.0, 500.0, 2_000.0, 5_000.0):
                index = int(np.argmin(np.abs(trajectory.time - time)))
                item[f"orientation_at_{time:g}_s"] = float(orientation[index])
            weak_summary.append(item)
    summaries["weak_cue_withdrawal"] = weak_summary

    _write_csv(OUTPUT / "time_courses.csv", rows)
    (OUTPUT / "summary.json").write_text(
        json.dumps(summaries, indent=2, sort_keys=True) + "\n"
    )
    print(OUTPUT / "summary.json")


if __name__ == "__main__":
    main()
