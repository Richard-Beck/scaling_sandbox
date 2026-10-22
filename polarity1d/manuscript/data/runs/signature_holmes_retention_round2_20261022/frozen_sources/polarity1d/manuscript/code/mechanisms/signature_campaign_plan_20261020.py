"""Compile the requested Figure 2 campaign without running native integrations.

This is an intake/planning tool. It has no submission or simulation command.
Unresolved numerical inputs stay null and launch readiness stays false. The
separate native preflight verifies exact retained checkpoints; no nearby length
or assumed Spring length conversion is accepted here.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / "polarity1d/manuscript").is_dir())
PACKAGE = ROOT / "polarity1d/manuscript"
BRIEFS = ROOT / "docs/tasks/briefs"
BRIEF = BRIEFS / "20261020_figure2_signature_assays.txt"
DECISIONS = BRIEFS / "20261020_figure2_decisions.json"
DEFAULT_OUTPUT = BRIEFS / "20261020_figure2_campaign_plan.json"
PREFLIGHT = Path(__file__).with_name("signature_native_preflight_20261020.py")
HARNESS = Path(__file__).with_name("signature_assays_20261020.py")
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
LENGTHS = (19, 23, 27, 31, 35, 39)
WIDTHS = (.95, 1.9, 3.8, 7.6, 15.2)
AREA_FACTORS = (.5, 1., 2.)
PULSE_COUNTS = (2, 4, 8)
GAP_MULTIPLES = (0., .1, .5, 1., 4., 16.)
ON_FACTORS = (1., 1.5, 2.)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def identity(path: Path) -> dict:
    exists = path.is_file()
    return dict(path=relative(path), exists=exists,
                sha256=digest(path) if exists else None)


def identifier(model: str, length: int) -> str:
    return f"{model}:L{length}"


def native_preflight() -> dict:
    if not PREFLIGHT.is_file():
        return dict(ready=False, model_basals={}, missing_requests=[],
                    blockers=["Native preflight module has not been written."],
                    dependency=identity(PREFLIGHT))
    spec = importlib.util.spec_from_file_location("signature_native_preflight", PREFLIGHT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load the native preflight module.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.load_verified_basals(package_root=PACKAGE,
        requested_lengths=LENGTHS, extra_lengths=(29,),
        spring_length_map={float(length): length / 15. for length in LENGTHS})
    result["dependency"] = identity(PREFLIGHT)
    return result


def basal_index(preflight: dict) -> dict:
    result = {}
    for model, records in preflight.get("model_basals", {}).items():
        for record in records:
            length = record.get("requested_length")
            if length is None:
                length = record.get("matched_reference_length")
            if length is not None:
                result[identifier(model, int(length))] = record
    return result


def base_record(model: str, length: int, basals: dict) -> dict:
    basal_id = identifier(model, length)
    basal = basals.get(basal_id, {})
    return dict(model=model, requested_reference_length=length,
                native_length=basal.get("native_length"), basal_id=basal_id,
                basal_status=basal.get("status", "preflight_not_available"),
                baseline_cells=201, status="planned_unresolved_not_executed")


def compile_plan() -> dict:
    if not BRIEF.is_file() or not DECISIONS.is_file():
        raise FileNotFoundError("The saved brief and user decisions are required.")
    decisions = json.loads(DECISIONS.read_text())
    calibration_decision = next((d for d in decisions.get("decisions", [])
        if d.get("topic") == "stage_0_commitment_calibration"), None)
    if not calibration_decision or calibration_decision.get("status") != "user_confirmed":
        raise ValueError("Post-withdrawal calibration requires its recorded user decision.")

    preflight = native_preflight()
    basals = basal_index(preflight)
    calibrations = []
    for model, length in itertools.product(MODELS, LENGTHS):
        calibrations.append(dict(**base_record(model, length, basals),
            calibration_id=identifier(model, length), role="stage_0_core",
            commitment_threshold_D=.5, criterion="post_withdrawal_replacement",
            recovery_duration=None, initial_duration_bracket=None,
            duration_resolution=None, measured_committing_bracket=None,
            field_tail_gate_independent=True))
    for model in ("wave_pinning", "debelly"):
        calibrations.append(dict(**base_record(model, 29, basals),
            calibration_id=identifier(model, 29), role="stage_3_extra_L29",
            commitment_threshold_D=.5, criterion="post_withdrawal_replacement",
            recovery_duration=None, initial_duration_bracket=None,
            duration_resolution=None, measured_committing_bracket=None,
            field_tail_gate_independent=True))

    width_rows = []
    for model, length, width, area in itertools.product(MODELS, LENGTHS, WIDTHS, AREA_FACTORS):
        width_rows.append(dict(**base_record(model, length, basals),
            condition_id=f"S1_{len(width_rows):04d}", requested_width_label=width,
            integrated_area_factor=area,
            brief_reference_area="native baseline peak times 3.8",
            native_width=width / 38. if model == "spring" else width,
            target_integrated_area=None, native_peak=None,
            measured_discrete_integrated_area=None,
            cue_geometry="new_local_native_raised_cosine_exact_discrete_commanded_area",
            challenger_zones=dict(fixed_widest_zone=None, each_cue_own_support=None),
            continuous_cue_horizon=None, replacement_threshold_D=.5,
            full_native_fields_required=True))

    train_rows = []
    for model, length, count, gap, on in itertools.product(
            MODELS, LENGTHS, PULSE_COUNTS, GAP_MULTIPLES, ON_FACTORS):
        train_rows.append(dict(**base_record(model, length, basals),
            condition_id=f"S2_{len(train_rows):04d}", pulse_count=count,
            gap_multiple_of_pulse=gap, total_on_factor_of_calibration=on,
            calibration_id=identifier(model, length),
            pulse_native_duration=None, total_native_on=None,
            actual_switch_times=None, actual_native_on=None,
            cue_peak=None, native_width=None, recovery_duration=None,
            record_D_at_every_pulse_and_gap_end=True,
            record_full_native_field_tail=True, field_tail_gate_independent=True,
            zero_gap_equivalence_candidate=(
                f"{identifier(model, length)}:continuous:{on:g}x" if gap == 0 else None),
            reuse_status="not_authorized_without_exact_input_and_output_checks"))

    noise_rows = []
    for model, length, level, gap, replicate in itertools.product(
            ("wave_pinning", "debelly"), (19, 29, 39),
            ("weak_1", "weak_2"), (1., 4., 16.), range(12)):
        noise_rows.append(dict(**base_record(model, length, basals),
            condition_id=f"S3_{len(noise_rows):04d}", noise_level_label=level,
            noise_level_value=None, replicate_index=replicate, random_seed=None,
            gap_multiple_of_pulse=gap, pulse_count=None, total_on_factor=None,
            calibration_id=identifier(model, length), forcing_law=None,
            forcing_channel=None, forcing_mean=None, forcing_correlation_time=None,
            forcing_correlation_length=None, noise_active_during=None,
            retention_definition=None, recovery_duration=None,
            record_full_native_field_tail=True, field_tail_gate_independent=True))

    unresolved = [
        dict(id="dominance_readout", detail="Verify the newly specified local D formula, incumbent/challenger regions, baseline/excess convention, integration weights and weak-signal handling."),
        dict(id="cue_shape_and_support", detail="Verify the newly specified native raised-cosine cue, support clipping and two challenger-zone conventions. External pilot geometry and numerical results are unavailable and are not reproduced claims."),
        dict(id="equal_area_normalization", detail="Verify exact sampled commanded area, including clipping and dynamic Spring geometry. Use literal native peak times 3.8; Spring native peak times .1."),
        dict(id="missing_holmes_basals", detail="Qualify new exact-L19/23/27/31/35/39 Holmes basals; the retained series contains different lengths."),
        dict(id="spring_length_mapping", detail="Qualify exact Spring resting lengths Lchemical/15, with cue widths scaled by 1/38. This is a declared comparison convention, not optical/biological calibration."),
        dict(id="calibration_bounds", detail="Declare native duration brackets, model-specific recovery/horizons, resolution and local success/failure scans; bisection alone does not establish global monotonicity."),
        dict(id="pulse_train_schedule", detail="Specify cue width/peak and event semantics. Keep exact restart states, realized total ON and full off histories; zero-gap physical duplicates require exact validation before reuse."),
        dict(id="debelly_time_lattice", detail="Declare realization of pulse/gap schedules on the native half-unit lattice; no unreported rounding or off-lattice output is allowed."),
        dict(id="noise_protocol", detail="Specify both weak levels, law/channel, mean, spatial/temporal correlations, activation windows, seeds, pulse count/total ON and retention criterion for linear WP/deBelly."),
        dict(id="tail_and_sampling", detail="Declare full native-field archive cadence, per-model recovery periods, finite-tail windows/scales/tolerances and event/sample output checks."),
        dict(id="retained_assays", detail="Run new local separation, share/dwell withdrawal and initiation overlays. The reported external 102 assays are unavailable and provide no verified evidence."),
        dict(id="numerical_design", detail="Select exact length pair and assay contrasts per model; qualify 401-cell initial states, output sampling, solver sensitivity and deBelly lattice checks."),
        dict(id="resource_pilot", detail="Measure a bounded new-protocol pilot before full native submission; historical resource figures are not an upper bound."),
    ]
    counts = dict(stage_0_core_calibrations=30, stage_0_extra_L29_calibrations=2,
                  stage_1_logical_conditions=len(width_rows),
                  stage_2_logical_conditions=len(train_rows),
                  stage_3_logical_conditions=len(noise_rows),
                  stages_1_to_3_logical_conditions=len(width_rows)+len(train_rows)+len(noise_rows),
                  stage_4_claimed_retained_assays=102,
                  stage_5_brief_approximate_runs=20)
    assert len(calibrations) == 32
    assert len(width_rows) == 450 and len(train_rows) == 1620 and len(noise_rows) == 432
    assert counts["stages_1_to_3_logical_conditions"] == 2502
    assert len({r["condition_id"] for r in width_rows+train_rows+noise_rows}) == 2502
    assert all(r["calibration_id"] in {c["calibration_id"] for c in calibrations}
               for r in train_rows+noise_rows)

    return dict(schema_version=1, generated_utc=datetime.now(timezone.utc).isoformat(),
        status="dry_run_plan_only", readiness=False, launch_enabled=False,
        simulations_run=0, jobs_submitted=0,
        sources=dict(brief=identity(BRIEF), decisions=identity(DECISIONS),
                     compiler=identity(Path(__file__)), named_pilot_harness=identity(HARNESS)),
        reconstruction_authorization="User authorized new local assays because the external harness/pilots are unavailable. Pilot numerical values remain provisional external descriptions, not verified local results.",
        spring_comparison_convention=dict(resting_length_rule="Lchemical/15",
            width_rule="requested width/38", reference_area="native peak times 0.1",
            biological_or_optical_calibration=False),
        confirmed_calibration_decision=calibration_decision,
        scientific_scope="Five native linear models. Ring noise and ring basals do not substitute for these assays.",
        protocol_counts=counts, native_preflight=preflight, basals_by_id=basals,
        calibrations=calibrations, stage_1_equal_area_width=width_rows,
        stage_2_cue_train=train_rows, stage_3_noisy_gap_hold=noise_rows,
        stage_4_retained_assays=dict(status="fresh_local_assays_required", external_claimed_count=102,
            exact_inputs_verified=0, exact_results_verified=0,
            required=["separation versus length", "share/dwell withdrawal map",
                      "initiation half/double strength and fixed-width overlays"],
            missing_data_is_censoring=False),
        stage_5_numerics=dict(status="design_and_execution_incomplete",
            requested_grids=[201, 401], requested_length_pairs_per_model=1,
            approximate_brief_runs=20, conditions=[],
            required=["Exact independently qualified grid-specific basal inputs",
                      "Same native reaction and mechanics parameters",
                      "Native event-time and output-sampling sensitivity",
                      "deBelly half-unit lattice checks",
                      "Readout and full-field-tail comparisons reported separately"]),
        unresolved_requirements=unresolved,
        execution_dependencies=[
            "Resolve pilot identities/readouts/geometry and exact basals",
            "Freeze executable sources, inputs, protocol and numerical acceptance rules",
            "Bounded fidelity/resource pilot and independent source/state/schedule audit",
            "Stage 0 core calibrations plus WP/deBelly L29 calibrations",
            "Stages 1/2 and fully specified Stage 3, with Stage 5 controls",
            "Independent complete native/readout/tail audit before figure claims"],
        native_wrapper_recommendations=dict(
            existing_five_model_platform="polarity1d/manuscript/code/mechanisms/revision_five_model_compensation_20261019.py",
            reusable_functions=["model_for", "sampled_run", "state_at", "save_phase", "save_state"],
            gates="gates() is prior orientation-based evidence only; new D fate gate requires its exact pilot definition.",
            chemical="Native BDF for Wave-pinning/Goryachev/Holmes; preserve native parameters and declare switching max_step.",
            spring="Native LSODA, rtol=2e-7, atol=2e-9, max_step=2.5; preserve raw u/v/boundaries and stretching. Existing output capture assumes one segment.",
            debelly="Corrected independent endpoints, center_boundary_motion=False, original staggered update step=0.5; preserve all coupled fields.",
            noise="Existing ring noise wrappers supply provenance patterns only; they do not implement linear deBelly noise."),
        resources=dict(
            existing_python="/home/4473331/.conda/envs/scaling_sandbox/bin/python",
            recommendation=dict(partition="red", qos="small", ntasks=1,
                cpus_per_task=1, memory="2G", initial_pilot_walltime="00:30:00",
                full_campaign_walltime="set_from_measured_new_pilot", array_throttle=None,
                thread_environment=dict(OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, MKL_NUM_THREADS=1)),
            historical_observations=dict(
                source="746 retained result.json records from 19 October baseline/strength workers",
                grouped_worker_max_wall_seconds=78.39,
                grouped_worker_max_peak_resident_MiB=329.0,
                caution="Those workers share a cue prefix and five withdrawals; new train/noise horizons can be much longer."),
            storage=dict(status="new_protocol_forecast_pending_pilot",
                historical_baseline_and_strength_allocated_bytes=70679977984,
                historical_measurement_utc="2026-10-19T22:09:00Z",
                historical_logical_protocols=3492,
                simplistic_2070_deterministic_conditions_analogy_GiB=(70679977984/3492*2070/1073741824),
                analogy_is_not_forecast_or_upper_bound=True,
                calibration_and_noise_storage_unknown=True,
                record_full_native_field_tails=True,
                manifest_requirement="Update root DATA_MANIFEST.md with measured allocated bytes/date/purpose for new or substantially enlarged ~2 GiB groups before handoff.",
                native_campaigns_stay_local=True)),
        plan_validation=dict(passed=True, unique_conditions=2502,
            exact_stage_counts_checked=True, calibration_dependencies_checked=True,
            numerical_inputs_inferred=False, scientific_results_claimed=False,
            validation_scope="Logical requested grid and dependency inventory only; numerical/scientific readiness is false."))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="Write the planning JSON; no integration or submission is supported.")
    args = parser.parse_args()
    plan = compile_plan()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError("Planning output must be inside this repository.")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(plan, indent=2, allow_nan=False) + "\n")
    temporary.replace(output)
    print(json.dumps(dict(output=relative(output), readiness=plan["readiness"],
        launch_enabled=plan["launch_enabled"], counts=plan["protocol_counts"],
        native_preflight_ready=plan["native_preflight"]["ready"],
        unresolved_requirements=len(plan["unresolved_requirements"])), indent=2))


if __name__ == "__main__":
    main()
