"""Verify retained native prerequisites for the 20 October signature brief.

This module reads existing files only. It imports no model/solver, runs no
integration, chooses no missing length mapping, and defines no new assay.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


PACKAGE = Path(__file__).resolve().parents[2]
CHEMICAL_LENGTHS = (19., 23., 27., 31., 35., 39.)
SPRING_LENGTH_MAP = {length: length / 15. for length in CHEMICAL_LENGTHS}
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
DB_NAMES = (
    "rac", "rho", "mca_bound", "mca_unbound", "displacement", "velocity",
    "tension", "raw_left_displacement", "raw_right_displacement",
    "left_boundary_displacement", "right_boundary_displacement",
)
STATE_NAMES = {
    "wave_pinning": ("active", "inactive"),
    "goryachev": ("active", "inactive"),
    "holmes_model3": ("rac_active", "rac_inactive", "rho_active", "rho_inactive"),
}


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(path: Path) -> dict:
    return json.loads(path.read_text())


def _relative(path: Path, repository: Path) -> str:
    return str(path.resolve().relative_to(repository.resolve()))


def _array_description(values: np.ndarray) -> dict:
    values = np.asarray(values)
    return {
        "shape": list(values.shape), "dtype": str(values.dtype),
        "finite": bool(np.isfinite(values).all()),
        "minimum": float(values.min()), "maximum": float(values.max()),
        "sha256_c_order_bytes": hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest(),
    }


def _descriptors(path: Path, case: dict) -> dict:
    """Describe exact native state and existing initial integration measures."""
    model, cells, length = case["model"], int(case["cells"]), float(case["length"])
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    bad = [name for name, value in arrays.items()
           if not np.issubdtype(value.dtype, np.number) or not np.isfinite(value).all()]
    if bad:
        raise ValueError("non-numeric/non-finite checkpoint arrays: " + ", ".join(bad))
    x = arrays["x"]
    if x.shape != (cells,) or not np.all(np.diff(x) > 0):
        raise ValueError("native x must contain exactly 201 ordered coordinates")
    if model == "debelly":
        state = {name: arrays["state_" + name] for name in DB_NAMES}
        if any(value.shape != (cells,) for value in state.values()):
            raise ValueError("deBelly state shape mismatch")
        active = state["rac"]
        state_description = {
            "encoding": "dictionary of state_<name> arrays", "names": list(DB_NAMES),
            "arrays": {name: _array_description(value) for name, value in state.items()},
            "native_step": 0.5, "center_boundary_motion": False,
            "boundary_method": "independent_endpoints",
        }
    else:
        state = arrays["state"]
        if model == "spring":
            if state.shape != (2 * cells + 2,):
                raise ValueError("Spring state shape mismatch")
            active = state[:-2:2]
            state_description = {
                "encoding": "interleaved active/inactive followed by two raw endpoints",
                "array": _array_description(state),
                "active_slice": "state[:-2:2]", "inactive_slice": "state[1:-2:2]",
                "raw_endpoints": state[-2:].tolist(),
            }
        else:
            names = STATE_NAMES[model]
            if state.shape != (len(names), cells):
                raise ValueError("reaction-diffusion state shape mismatch")
            active = state[0]
            state_description = {
                "encoding": "state array with field rows", "names": list(names),
                "array": _array_description(state), "active_row": 0,
            }
    np.testing.assert_array_equal(arrays["active"][-1], active)
    if model == "spring":
        actual_length = float(state[-1] - state[-2])
        if actual_length <= 0:
            raise ValueError("Spring endpoints are not ordered")
        q = np.linspace(0., 1., cells)
        np.testing.assert_allclose(x, q * length, rtol=0, atol=1e-12)
        physical_x = float(state[-2]) + q * actual_length
        weights = np.full(cells, actual_length / (cells - 1))
        weights[[0, -1]] *= .5
        coordinate_system = "material endpoint nodes; checkpoint x uses resting length"
        weight_rule = "initial physical trapezoidal weights; half weight at each endpoint"
    else:
        actual_length = length
        np.testing.assert_allclose(x, (np.arange(cells) + .5) * length / cells,
                                   rtol=0, atol=1e-12)
        physical_x = x
        weights = np.full(cells, length / cells)
        coordinate_system = ("fixed nominal solver cell centres; displacement retained separately"
                             if model == "debelly" else "fixed physical cell centres")
        weight_rule = ("nominal solver-cell widths; not a deformed-membrane measure"
                       if model == "debelly" else "physical finite-volume cell widths")
    active_description = _array_description(active)
    active_description.update({"field": "rac" if model == "debelly" else
                               "rac_active" if model == "holmes_model3" else "active",
                               "checkpoint_field_matches_state": True})
    return {
        "all_checkpoint_arrays_finite": True, "checkpoint_arrays_checked": len(arrays),
        "state": state_description, "active": active_description,
        "coordinates": {
            "system": coordinate_system, "checkpoint_x": _array_description(x),
            "checkpoint_x_first_last": [float(x[0]), float(x[-1])],
            "initial_measure_coordinates_first_last": [float(physical_x[0]), float(physical_x[-1])],
            "nominal_length": length, "initial_measure_length": actual_length,
            "source_convention": "revision_five_model_compensation_20261019.regions",
        },
        "weights": {
            "rule": weight_rule, "shape": [cells], "sum": float(weights.sum()),
            "first": float(weights[0]), "interior": float(weights[cells // 2]),
            "last": float(weights[-1]),
            "role": "existing initial integration measure; no challenger zone selected",
        },
    }


def _missing_record(model: str, requested: float | None, reason: str,
                    reference: float | None = None) -> dict:
    return {
        "model": model, "requested_length": requested, "native_length": None,
        "matched_reference_length": reference, "checkpoint": None,
        "expected_sha256": None, "actual_sha256": None, "exists": False,
        "sha256_verified": False, "cells": None, "status": reason,
    }


def load_verified_basals(
    package_root: str | Path | None = None,
    requested_lengths: Sequence[float] = CHEMICAL_LENGTHS,
    extra_lengths: Sequence[float] = (29.,),
    spring_length_map: Mapping[float, float] | None = SPRING_LENGTH_MAP,
) -> dict:
    """Return JSON-ready exact native prerequisites, without filling missing work.

    ``extra_lengths`` requests Wave-pinning/deBelly noise controls only. An
    explicit ``spring_length_map`` maps the brief's reference lengths to native
    resting lengths. The default implements the authorized comparison convention
    L0=L/15, not a calibrated conversion. Pass None to audit unassigned mapping.
    """
    package = Path(package_root).resolve() if package_root else PACKAGE
    repository = package.parents[1]
    source = package / "data/runs/debelly_corrected_three_test_20261008"
    manifest_path = source / "manifest.json"
    manifest = _read(manifest_path)
    provenance = _read(source / "provenance.json")
    terminal = _read(source / "results/status.json")
    validation = _read(source / "results/validation.json")
    actual_manifest_sha = _sha(manifest_path)
    blockers = []
    if actual_manifest_sha != provenance["manifest_sha256"]:
        blockers.append("Corrected source manifest fails its recorded SHA-256.")
    if not (terminal.get("complete") and validation.get("passed")
            and validation.get("complete_execution") and validation.get("independent_audit_passed")):
        blockers.append("Corrected source campaign lacks complete independent validation.")
    if manifest.get("center_boundary_motion") is not False:
        blockers.append("Corrected manifest does not disable centred deBelly endpoints.")
    audit = source / "results/independent_audit.json"
    audit_verified = audit.is_file() and _sha(audit) == validation.get("independent_audit_sha256")
    if not audit_verified:
        blockers.append("Corrected source independent audit hash is unavailable or mismatched.")
    cases = [case for case in manifest["cases"] if case["condition"] == "fixed_physical"]
    requested_lengths = tuple(float(value) for value in requested_lengths)
    extra_lengths = tuple(float(value) for value in extra_lengths)
    model_basals, missing_requests, available = {}, [], {}
    for model in MODELS:
        candidates = [case for case in cases if case["model"] == model and case["eligible"]]
        available[model] = sorted(float(case["length"]) for case in candidates)
        lengths = list(requested_lengths)
        if model in ("wave_pinning", "debelly"):
            lengths += [value for value in extra_lengths if value not in lengths]
        records = []
        for requested in lengths:
            if model == "spring" and spring_length_map is None:
                record = _missing_record(model, None, "resting_length_mapping_required", requested)
                records.append(record)
                missing_requests.append(record.copy())
                continue
            native_length = (float(spring_length_map[requested])
                             if model == "spring" and requested in spring_length_map else
                             None if model == "spring" else requested)
            matches = [case for case in candidates if case["length"] == native_length]
            if len(matches) != 1:
                reason = "exact_basal_missing" if not matches else "ambiguous_exact_basal"
                record = _missing_record(model, requested, reason)
                record["native_length"] = native_length
                records.append(record)
                missing_requests.append(record.copy())
                continue
            case = matches[0]
            checkpoint = package / case["checkpoint"]
            exists = checkpoint.is_file()
            actual = _sha(checkpoint) if exists else None
            record = {
                "model": model, "requested_length": requested, "native_length": native_length,
                "checkpoint": _relative(checkpoint, repository),
                "expected_sha256": case["checkpoint_sha256"], "actual_sha256": actual,
                "exists": exists, "sha256_verified": actual == case["checkpoint_sha256"],
                "cells": int(case["cells"]), "status": "verified",
                "source_case_index": case["case_index"],
                "native_response_peak": case["response_amplitude"],
                "native_background": case["background_stimulus"],
                "fixed_scales": case["fixed_scales"],
                "existing_gates": {
                    name: case.get(name) for name in (
                        "eligible", "preparation_status", "preparation_condition",
                        "initial_orientation", "initial_contrast", "boundary_method",
                        "center_boundary_motion",
                    )
                },
                "gate_scope": "retained preparation eligibility; new withdrawn field tails not evaluated",
            }
            try:
                if not exists or not record["sha256_verified"]:
                    raise ValueError("checkpoint missing or SHA-256 mismatch")
                if record["cells"] != 201:
                    raise ValueError("requested retained discretization is 201 cells")
                if case["preparation_status"] != "polarized":
                    raise ValueError("retained preparation is not polarized")
                if model == "debelly" and (case.get("center_boundary_motion") is not False
                                          or case.get("boundary_method") != "independent_endpoints"):
                    raise ValueError("deBelly independent-endpoint qualification missing")
                record.update(_descriptors(checkpoint, case))
            except (ValueError, KeyError, AssertionError, OSError) as error:
                record["status"] = "verification_failed"
                record["error"] = str(error)
                blockers.append(f"{model} L{requested:g}: {error}")
            records.append(record)
        model_basals[model] = records
    if any(record["model"] == "holmes_model3" for record in missing_requests):
        blockers.append("Requested exact Holmes lengths have no retained qualifying basals; fresh preparation is required.")
    if any(record["model"] == "spring" for record in missing_requests):
        blockers.append("Spring exact mapped resting-length basals need fresh preparation."
                        if spring_length_map is not None else
                        "Spring matched resting lengths are unassigned; no mapping is inferred.")
    checked_sources = [
        manifest_path, source / "provenance.json", source / "results/status.json",
        source / "results/validation.json", audit,
        package / "code/core/native.py", package / "code/core/polarity_establishment.py",
        package / "code/mechanisms/revision_five_model_compensation_20261019.py",
        Path(__file__),
    ]
    return {
        "schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "ready": not blockers and not missing_requests,
        "integrations_run": 0, "solver_imports": 0,
        "requested_chemical_lengths": list(requested_lengths),
        "extra_noise_lengths_wp_debelly": list(extra_lengths),
        "spring_mapping_convention": ("L0=L/15; authorized comparison convention, not calibrated physical conversion"
                                      if spring_length_map == SPRING_LENGTH_MAP else
                                      "unassigned" if spring_length_map is None else "explicit caller-supplied mapping"),
        "model_basals": model_basals, "available_exact_lengths": available,
        "verified_basal_count": sum(record["status"] == "verified"
                                    for records in model_basals.values() for record in records),
        "missing_requests": missing_requests, "blockers": blockers,
        "source_manifest": {
            "path": _relative(manifest_path, repository),
            "expected_sha256": provenance["manifest_sha256"], "actual_sha256": actual_manifest_sha,
            "sha256_verified": actual_manifest_sha == provenance["manifest_sha256"],
            "terminal_complete": terminal.get("complete"),
            "independent_validation_passed": validation.get("independent_audit_passed"),
            "independent_audit_sha256_verified": audit_verified,
        },
        "source_hashes": {_relative(path, repository): _sha(path)
                          for path in checked_sources if path.is_file()},
        "retained_tail_gate_context": {
            "recovery_horizon": manifest["recovery_horizon"],
            "sample_interval": manifest["sample_interval"],
            "gate_window": manifest["gate_window"],
            "new_assay_fate_definition": None, "new_assay_field_tail_evaluated": False,
            "note": "Preparation eligibility, scalar fate and native full-field tail gates remain separate.",
        },
        "existing_cue_geometry": {
            "source": "polarity1d/manuscript/code/core/polarity_establishment.py:41",
            "shape": "right-boundary half-cosine on distance d in [0,w]",
            "signal_above_background": "a * (1 + cos(pi*d/w)) / 2",
            "continuum_integral_if_0_lt_w_le_L": "a*w/2",
            "sampled_integral": "sum(signal_above_background * model_native_weights)",
            "new_brief_area_convention": None,
            "warning": "Brief's native peak times3.8 area may use another cue convention. No shape or factor-of-two conversion has been assumed; domain clipping and discrete integration also matter.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=PACKAGE.parents[1] / "docs/tasks/briefs/20261020_figure2_native_preflight.json")
    args = parser.parse_args()
    report = load_verified_basals()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "ready": report["ready"],
                      "verified_basal_count": report["verified_basal_count"],
                      "missing_requests": len(report["missing_requests"]),
                      "blockers": report["blockers"]}, indent=2))


if __name__ == "__main__":
    main()
