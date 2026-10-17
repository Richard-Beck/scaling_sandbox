#!/usr/bin/env python3
"""Cross-validate the Python de Belly port against native Julia/Gridap trajectories."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
sys.path.insert(0, str(ROOT / "src"))

from polarity1d import (  # noqa: E402
    AlternatingGradient,
    DeBelly,
    DeBellyParameters,
    LocalizedPulse,
    ReversingGradient,
    simulate,
)


CASES = {
    "competition_l5": {
        "length": 5.0,
        "end_time": 200.0,
        "output_interval": 0.5,
        "protocol": LocalizedPulse(
            amplitude=1.0, start=0.0, stop=40.0, taper=20.0,
            width_fraction=0.1, sides=("left", "right"), asymmetry=0.02,
        ),
        "julia_protocol": "localized_pulse",
    },
    "reversal_l15": {
        "length": 15.0,
        "end_time": 2300.0,
        "output_interval": 1.0,
        "protocol": ReversingGradient(slope=1.0 / 15.0, reversal_time=2000.0, transition=20.0),
        "julia_protocol": "reversing_gradient",
        "julia_extra": {"DEBELLY_COMMON_REVERSAL_TIME": "2000"},
    },
    "tracking_l39": {
        "length": 39.0,
        "end_time": 1200.0,
        "output_interval": 2.5,
        "protocol": AlternatingGradient(slope=1.0 / 39.0, period=200.0, start=0.0),
        "julia_protocol": "alternating_gradient",
        "julia_extra": {"DEBELLY_COMMON_PERIOD": "200"},
    },
}


def _read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def _nrmse(reference: np.ndarray, candidate: np.ndarray) -> float:
    span = float(np.max(reference) - np.min(reference))
    return float(np.sqrt(np.mean((candidate - reference) ** 2)) / max(span, 1e-12))


def _relative_l2(reference: np.ndarray, candidate: np.ndarray) -> float:
    return float(np.linalg.norm(candidate - reference) / max(np.linalg.norm(reference), 1e-12))


def _orientation(values: np.ndarray, x: np.ndarray, length: float) -> np.ndarray:
    left = np.mean(values[:, x <= 0.2 * length], axis=1)
    right = np.mean(values[:, x >= 0.8 * length], axis=1)
    return (right - left) / np.maximum(np.abs(left) + np.abs(right), 1e-12)


def run_case(case_name: str, output_root: Path, *, reuse_julia: bool = False) -> Path:
    case = CASES[case_name]
    length = float(case["length"])
    cells = round(length / 0.2)
    destination = output_root / case_name
    julia_output = destination / "julia"
    destination.mkdir(parents=True, exist_ok=True)

    python_model = DeBelly(
        DeBellyParameters(
            initialization_mode="julia_warmup",
            assay_mode="published_one_sided",
            input_sampling_mode="julia_element_right_edge",
        )
    )
    trajectory = simulate(
        python_model,
        protocol=case["protocol"], length=length, cells=cells,
        end_time=float(case["end_time"]),
        output_interval=float(case["output_interval"]), step=0.5,
    )
    trajectory.save(destination / "python_trajectory.npz")

    environment = os.environ.copy()
    environment.update(
        {
            "JULIA": "/app/eb/software/Julia/1.8.2-linux-x86_64/bin/julia",
            "DEBELLY_WRITE_OUTPUTS": "0",
            "DEBELLY_QUIET": "1",
            "DEBELLY_T": str(case["end_time"]),
            "DEBELLY_DT": "0.5",
            "DEBELLY_PARTITION": str(cells),
            "DEBELLY_LENGTH": str(length),
            "DEBELLY_MCA_DENSITY": "0.1",
            "DEBELLY_COMMON_PROTOCOL": str(case["julia_protocol"]),
            "DEBELLY_COMMON_PEAK": "1",
            "DEBELLY_COMMON_ASYMMETRY": "0.02",
            "DEBELLY_TRAJECTORY_OUTPUT": str(julia_output),
        }
    )
    environment.update(case.get("julia_extra", {}))
    runner = (
        REPOSITORY
        / "literature_models/debelly_2026_mechanochemical/scripts/run_debelly.sh"
    )
    if reuse_julia:
        if not (julia_output / "trajectory.csv").is_file():
            raise FileNotFoundError(
                f"cannot reuse missing Julia reference {julia_output / 'trajectory.csv'}"
            )
        julia_exit_code = 0
    else:
        completed = subprocess.run(
            [str(runner), "mechanochemical", "published"], cwd=REPOSITORY,
            env=environment, text=True, capture_output=True, check=False,
        )
        (destination / "julia.stdout").write_text(completed.stdout)
        (destination / "julia.stderr").write_text(completed.stderr)
        if completed.returncode:
            raise RuntimeError(f"Julia failed with exit code {completed.returncode}")
        julia_exit_code = completed.returncode

    chemistry = _read_csv(julia_output / "trajectory.csv")
    result = {
        "case": case_name,
        "model": "2026 de Belly mechanochemical model",
        "length_um": length,
        "cells": cells,
        "dx_um": length / cells,
        "dt_s": 0.5,
        "python_backend": trajectory.metadata["backend"],
        "initialization_mode": "julia_warmup",
        "initialization_scope": "cross-validation only",
        "python_protocol_sampling": (
            "Julia adapter parity: element value j*dx, implemented as a +dx/2 "
            "offset from finite-volume cell centres"
        ),
        "assay_mode": "published_one_sided",
        "julia_backend": "Gridap P1 finite elements, staggered backward Euler",
        "julia_exit_code": julia_exit_code,
        "julia_reference_reused": reuse_julia,
        "fields": {},
    }
    mechanics = _read_csv(julia_output / "mechanics.csv")
    sources = {
        "rac": (chemistry, "rac", True),
        "rho": (chemistry, "rho", True),
        "mca_bound": (chemistry, "mca_bound", True),
        "mca_unbound": (chemistry, "mca_unbound", True),
        "velocity": (mechanics, "velocity", False),
        "displacement": (mechanics, "displacement", False),
        "tension": (_read_csv(julia_output / "tension.csv"), "tension", True),
    }
    for field, (source, source_field, compare_orientation) in sources.items():
        julia_times = np.unique(source["time_s"])
        julia_x = np.unique(source["position_um"])
        raw = source[source_field].reshape(len(julia_times), len(julia_x))
        selected_times = np.asarray(
            [time for time in trajectory.time if np.any(np.isclose(julia_times, time))]
        )
        python_indices = [
            int(np.argmin(np.abs(trajectory.time - time))) for time in selected_times
        ]
        julia_time_indices = [
            int(np.argmin(np.abs(julia_times - time))) for time in selected_times
        ]
        julia_values = raw[julia_time_indices]
        # Mechanics are exported on different interior Gridap spaces. Compare on
        # each Julia field's native coordinates rather than extrapolating its
        # interior values to the Python finite-volume boundary cells.
        native_python_values = np.asarray(
            [
                np.interp(julia_x, trajectory.x, trajectory.fields[field][index])
                for index in python_indices
            ]
        )
        comparison_representation = "native field interpolated to Julia coordinates"
        native_representation_relative_l2 = None
        if field == "tension":
            # The Julia exporter does not write tension at its natural interval
            # midpoint.  It takes backward differences of the interior P1
            # displacement DOFs and labels each value by the right-hand node.
            # Reconstruct that same diagnostic from the Python displacement;
            # retain the direct centred-gradient discrepancy separately.
            mechanics_x = np.unique(mechanics["position_um"])
            python_displacement = np.asarray(
                [
                    np.interp(
                        mechanics_x, trajectory.x,
                        trajectory.fields["displacement"][index],
                    )
                    for index in python_indices
                ]
            )
            python_values = (
                python_model.parameters.membrane_stiffness
                * np.diff(python_displacement, axis=1) / (length / cells)
            )
            if python_values.shape != julia_values.shape:
                raise RuntimeError("Julia tension diagnostic grid did not match mechanics grid")
            native_representation_relative_l2 = _relative_l2(
                julia_values, native_python_values
            )
            comparison_representation = (
                "backward displacement difference on Julia interior-node grid; "
                "direct centred-FV tension error retained separately"
            )
        else:
            python_values = native_python_values
        if compare_orientation:
            julia_orientation = _orientation(julia_values, julia_x, length)
            python_orientation = _orientation(python_values, julia_x, length)
            orientation_rmse = float(
                np.sqrt(np.mean((julia_orientation - python_orientation) ** 2))
            )
        else:
            orientation_rmse = None
        result["fields"][field] = {
            "trajectory_nrmse_span": _nrmse(julia_values, python_values),
            "terminal_nrmse_span": _nrmse(julia_values[-1], python_values[-1]),
            "trajectory_relative_l2": _relative_l2(julia_values, python_values),
            "orientation_rmse": orientation_rmse,
            "julia_minimum": float(np.min(julia_values)),
            "python_minimum": float(np.min(python_values)),
            "comparison_representation": comparison_representation,
            "native_representation_relative_l2": native_representation_relative_l2,
        }
    result["acceptance"] = {
        "all_concentrations_nonnegative": bool(
            all(
                result["fields"][field]["julia_minimum"] >= -1e-8
                and result["fields"][field]["python_minimum"] >= -1e-8
                for field in ("rac", "rho", "mca_bound", "mca_unbound")
            )
        ),
        "maximum_trajectory_relative_l2": float(
            max(values["trajectory_relative_l2"] for values in result["fields"].values())
        ),
        "maximum_orientation_rmse": float(
            max(
                values["orientation_rmse"] for values in result["fields"].values()
                if values["orientation_rmse"] is not None
            )
        ),
        "passed": bool(
            max(values["trajectory_relative_l2"] for values in result["fields"].values()) <= 0.05
            and max(
                values["orientation_rmse"] for values in result["fields"].values()
                if values["orientation_rmse"] is not None
            ) <= 0.05
        ),
        "criteria": (
            "all seven exported state fields relative L2 <= 0.05; Rac, Rho, "
            "tension and both MCA fields orientation RMSE <= 0.05"
        ),
    }
    result_path = destination / "cross_validation.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=tuple(CASES))
    parser.add_argument(
        "--output-root", type=Path,
        default=ROOT / "results" / "cross_validation" / "debelly_2026_v1",
    )
    parser.add_argument(
        "--reuse-julia",
        action="store_true",
        help="compare with an existing Julia trajectory without executing Julia",
    )
    arguments = parser.parse_args(argv)
    print(run_case(arguments.case, arguments.output_root, reuse_julia=arguments.reuse_julia))


if __name__ == "__main__":
    main()
