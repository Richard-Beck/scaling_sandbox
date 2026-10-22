"""Corrected complete-active-region Figure 3 A2 and B experiments.

All challenges restart exact independently qualified basals. The beam support
sets the native input only; it never truncates the challenger amount. A2 and B
freeze separately so dependency collection cannot mutate a running producer.
"""
from pathlib import Path
import argparse
import copy
import json
import os
import shlex
import shutil
import time

import numpy as np

import signature_figure3_native_20261020 as native
import signature_figure3_campaign_20261020 as original_campaign
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_regions_20261020 as regions
import signature_figure3_events_20261020 as events

ROOT = native.ROOT
DEFAULT_B = ROOT/"data/runs/signature_figure3_regional_B_20261020"
DEFAULT_A2 = ROOT/"data/runs/signature_figure3_regional_A2_20261020"
REFINEMENT = ROOT/"data/runs/signature_figure3_handover_refinement_20261020"
read, write, sha, now = native.read, native.write, native.sha, native.now


def spatial(adapter, state, placement, detailed=True):
    old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
    return regions.readout(adapter.primary(state), adapter.circumference(state), old,
                           placement["cue_center_fraction"], detailed=detailed)


def measured_gap(adapter, state, placement, fraction=.1):
    C = adapter.circumference(state)
    old = placement["initial_10_percent_contour"]["origin"]*C/placement["basal_circumference"]
    result = geometry.measured_regions(adapter.q*C, C, adapter.primary(state), old,
        placement["cue_center_fraction"]*C, fraction)
    result.update(contour_definition="Global instantaneous floor + fraction × whole-field contrast",
                  threshold=float(adapter.primary(state).min()+fraction*np.ptp(adapter.primary(state))))
    return result


def fate(adapter, state, placement):
    value = spatial(adapter, state, placement, detailed=False)
    labels = {"single_region_cue_associated": "challenger", "single_region_incumbent_associated": "incumbent",
              "single_region_between_sites": "location_ambiguous"}
    return dict(winner=labels.get(value["state_class"], "unresolved"), state_class=value["state_class"],
        resolved_peaks=value["resolved_peaks"], contrast=value["contrast"],
        vector=geometry.activity_vector(adapter.primary(state), adapter.q))


def freeze(run, arm, original, refinement):
    if run.exists():
        raise FileExistsError("Fresh immutable regional campaign directory required")
    primary = original_campaign.verify(original)
    for folder in ("inputs", "frozen_sources", "workers", "results", "logs"):
        (run/folder).mkdir(parents=True)
    sources, inputs = {}, {}
    for relative, digest in primary["source_sha256"].items():
        source, target = original/"frozen_sources"/relative, run/"frozen_sources"/relative
        assert sha(source) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        sources[relative] = digest
    for source in (Path(__file__), Path(regions.__file__), Path(events.__file__)):
        relative = str(source.relative_to(ROOT))
        shutil.copy2(source, run/"frozen_sources"/relative)
        sources[relative] = sha(source)
    for relative, digest in primary["input_sha256"].items():
        source, target = original/relative, run/relative
        assert sha(source) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        inputs[relative] = digest
    jobs, unavailable = [], []
    if arm == "B":
        for placement in primary["placements"]:
            if placement["placement_status"] != "feasible" or placement["classification_10"] != "external":
                continue
            for ratio in (.5, 1., 2.):
                value = dict(slot=len(jobs), arm="B", placement_id=placement["placement_id"], target_ratio=ratio)
                if placement.get("control_variant"):
                    reference = next(j for j in primary["stages"]["controls"] if j["placement_id"] == placement["placement_id"])
                    value.update(control=placement["control_variant"], matched_reference_gap=placement["matched_gap"],
                                 solver_options=reference.get("solver_options"))
                jobs.append(value)
        for old_job in primary["stages"].get("primary", []):
            if old_job.get("control") == "mirror_placement":
                jobs.append(dict(slot=len(jobs), arm="B", placement_id=old_job["placement_id"],
                    target_ratio=1., mirrored=True, control="mirror_placement"))
    else:
        refined = original_campaign.verify(refinement)
        auditpath = refinement/"results/independent_refinement_audit.json"
        audit = read(auditpath)
        assert audit["complete"] and audit["source_manifest_sha256"] == sha(refinement/"manifest.json")
        refined_jobs = {j["slot"]: j for j in refined["jobs"]}
        accepted = {}
        for record in audit["records"]:
            job = refined_jobs[record["slot"]]
            if job.get("control"):
                continue
            assert sha(Path(record["status_path"])) == record["status_sha256"]
            accepted[job["placement"]["placement_id"]] = record
        for placement in primary["placements"]:
            if placement["placement_status"] != "feasible" or placement["classification_10"] not in ("external", "borderline"):
                continue
            record = accepted.get(placement["placement_id"])
            if record is None:
                unavailable.append(dict(placement_id=placement["placement_id"], model=placement["model"],
                    size_parameter=placement["size_parameter"], classification=placement["classification_10"],
                    status="no_continuous_two_region_handover_available", cap_substitution=False))
                continue
            for multiplier in (.5, .75, 1., 1.5, 2.):
                jobs.append(dict(slot=len(jobs), arm="A2", placement_id=placement["placement_id"], multiplier=multiplier,
                    handover_time=record["handover_time"], cue_duration=multiplier*record["handover_time"],
                    source_refinement_status=record["status_path"], source_refinement_status_sha256=record["status_sha256"]))
    plan = copy.deepcopy(primary)
    plan.update(generated_at=now(), experiment="Complete active-region Figure 3 "+arm,
        source_run=str(original), source_manifest_sha256=sha(original/"manifest.json"),
        source_sha256=sources, input_sha256=inputs, jobs=jobs, unavailable=unavailable, arm=arm,
        worker_soft_seconds=3480., root_definition="Complete valley-bounded active excess above common current minimum",
        region_definition_authority="User: yes obviously you should not restrict challenger activity counting to the area under the cue.",
        primary_gap_contour_fraction=.1, other_gap_contours=[.25, .5, .75],
        B_event_honesty="Continuous ratio equality only with two resolved regions throughout a valid dense bracket; topology/birth boundaries explicit",
        B_claim_gate="Continuous ratio event, exactly two resolved associated peaks, positive complete-region amounts, and two separated global10%-contrast contours",
        B_resolution_definition="First single resolved native peak after withdrawal, prominence0.05 and contrast>=0.05; global10/25/50/75 contour coalescence recorded separately",
        geometry="Actual measured region centers and widths; no extrapolation or assumed1/C scaling",
        native_kinetics_changed=False)
    if arm == "A2":
        plan.update(refinement_manifest=str(refinement/"manifest.json"), refinement_manifest_sha256=sha(refinement/"manifest.json"),
                    refinement_audit=str(auditpath), refinement_audit_sha256=sha(auditpath))
    write(run/"manifest.json", plan)
    write(run/"freeze.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    runner = run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command = "exec "+shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", native.PYTHON, str(runner), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    write(run/"slurm_payload.json", dict(generated_at=now(), submitted=False, jobs=len(jobs), no_concurrency_throttle=True,
        sbatch_arguments=["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=4G",
            "--time=01:00:00", "--signal=B:USR1@60", f"--array=0-{len(jobs)-1}", "--job-name=figure3-region-"+arm,
            "--output="+str(run/"logs/%A_%a.out"), "--error="+str(run/"logs/%A_%a.err"), "--wrap", command]))
    return dict(jobs=len(jobs), unavailable=len(unavailable), arm=arm, native_integrations=0)


def trace_arrays(adapter, times, states, placement):
    records = [spatial(adapter, y, placement, detailed=False) for y in states]
    vectors = [geometry.activity_vector(adapter.primary(y), adapter.q) for y in states]
    arrays = dict(time=np.asarray(times), state_class=np.asarray([r["state_class"] for r in records]),
        resolved_peaks=np.asarray([r["resolved_peaks"] for r in records]),
        vector_real=np.asarray([v["real"] for v in vectors]), vector_imag=np.asarray([v["imag"] for v in vectors]),
        vector_magnitude=np.asarray([v["magnitude"] for v in vectors]),
        vector_heading=np.asarray([v["heading_radians"] if v["direction_available"] else np.nan for v in vectors]))
    for key in ("incumbent_excess", "challenger_excess", "challenger_to_incumbent_ratio", "regional_dominance"):
        arrays[key] = np.asarray([r.get(key, np.nan) if r.get(key) is not None else np.nan for r in records])
    return arrays


def resolution_summary(adapter, times, states, placement):
    fates = [fate(adapter, y, placement) for y in states]
    single = np.asarray([v["winner"] in ("incumbent", "challenger", "location_ambiguous") for v in fates])
    primary = original_campaign.first_bracket(times, single)
    found = np.flatnonzero(single)
    first = int(found[0]) if len(found) else None
    contours = {}
    for fraction in (.1, .25, .5, .75):
        groups = [len(geometry.components(adapter.primary(y) >= adapter.primary(y).min()+fraction*np.ptp(adapter.primary(y))))
                  if np.ptp(adapter.primary(y)) >= .05 else 0 for y in states]
        contours[str(fraction)] = dict(initial_components=groups[0], final_components=groups[-1],
            first_one_component=original_campaign.first_bracket(times, np.asarray(groups) == 1)
                if groups[0] == 2 else dict(status="unavailable_not_two_contour_regions_at_withdrawal"),
            definition="Current global contrast contour; coalescence is distinct from disappearance of a resolved activity peak")
    return dict(primary=primary, final_winner=fates[-1], contour_coalescence=contours,
                first_single_region_fate=None if first is None else fates[first],
                one_region_persists_through_cap=False if first is None else bool(single[first:].all()),
                later_multiple_or_ambiguous_state_after_first_single=False if first is None else bool((~single[first:]).any()),
                initial_state_class=fates[0]["state_class"], final_state_class=fates[-1]["state_class"])


def worker(run, slot):
    native.require_slurm()
    plan = original_campaign.verify(run)
    job = plan["jobs"][slot]
    placement = copy.deepcopy(plan["placements"][job["placement_id"]])
    if job.get("mirrored"):
        old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
        placement["cue_center_fraction"] = float((2*old-placement["cue_center_fraction"]) % 1)
    adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
    with np.load(run/placement["input"]) as z:
        initial = z["state"].copy()
    directory = run/"workers"/f"{slot:04d}"
    directory.mkdir(parents=True)
    statuspath = directory/"status.json"
    record = dict(generated_at=now(), status="running", complete=False, job=job, placement=placement,
        slurm_job_id=os.environ["SLURM_JOB_ID"], phases=[])
    write(statuspath, record)
    elapsed, state, history_t, history_y = 0., initial.copy(), [], []
    deadline = time.monotonic()+plan["worker_soft_seconds"]
    def phase(kind, duration, peak, detector=None):
        nonlocal elapsed, state
        record.update(status=kind, generated_at=now())
        write(statuspath, record)
        result = (events.integrate_seed(adapter, state, duration, peak, placement["cue_center_fraction"], detector, deadline=deadline,
                                       solver_options=job.get("solver_options"))
                  if detector is not None else native.integrate(adapter, state, duration, peak, placement["cue_center_fraction"], deadline=deadline,
                                                               solver_options=job.get("solver_options")))
        archive = directory/f'{len(record["phases"]):03d}_{kind}.npz'
        metadata = native.save_phase(archive, result, adapter, peak, placement["cue_center_fraction"])
        metadata.update(path=str(archive.relative_to(run)), absolute_start=elapsed, kind=kind)
        record["phases"].append(metadata)
        use = slice(None) if not history_t else slice(1, None)
        history_t.extend((elapsed+result["time"][use]).tolist()); history_y.extend(result["states"][use])
        elapsed += float(result["time"][-1]); state = result["states"][-1].copy()
        write(statuspath, record)
        if result["status"] not in ("complete", "event"):
            raise RuntimeError("Incomplete native phase: "+result["status"])
        return result
    try:
        basal = spatial(adapter, initial, placement, detailed=False)
        record["basal_complete_excess"] = basal["full_excess_amount"]
        cap = 2000. if adapter.spring else 20000.
        if job["arm"] == "A2":
            assert sha(Path(job["source_refinement_status"])) == job["source_refinement_status_sha256"]
            phase("pulse", job["cue_duration"], placement["cue_peak"])
            record.update(cue_duration=job["cue_duration"], withdrawal_time=elapsed,
                withdrawal_spatial=spatial(adapter, state, placement))
            tail = phase("cue_free_tail", 1300. if adapter.spring else 13000., 0.)
            record.update(final_winner=fate(adapter, state, placement), finite_tail=native.finite_tail(adapter, tail["time"], tail["states"]))
        else:
            detector = events.RegionSeedEvent(adapter, placement, job["target_ratio"], basal["full_excess_amount"])
            seed = phase("seed", cap, placement["cue_peak"], detector)
            record.update(target_ratio=job["target_ratio"], seed_event=seed["event_record"])
            if seed["event_record"] is None:
                record.update(seed_status="target_unreached_at_completed_seed_cap", seed_lower_bound=elapsed,
                    B_distance_claim_eligible=False, B_resolution=dict(status="unavailable_no_withdrawal"))
            else:
                withdrawal, withdrawal_time = state.copy(), elapsed
                spatial_value = spatial(adapter, withdrawal, placement)
                gaps = {str(f): measured_gap(adapter, withdrawal, placement, f) for f in (.1, .25, .5, .75)}
                kind = seed["event_record"]["kind"]
                positive = spatial_value.get("incumbent_excess", 0.) > 1e-10*basal["full_excess_amount"]
                continuous = kind == "continuous_two_region_ratio" and seed["event_record"]["continuous_root"]
                coincident_collapse = "continuous_incumbent_90_percent_excess_loss" in seed["event_record"].get("coincident_event_kinds", []) and continuous
                eligible = bool(continuous and not coincident_collapse and spatial_value["two_region_available"] and positive and gaps["0.1"]["available"])
                if continuous:
                    difference = spatial_value["challenger_excess"]-job["target_ratio"]*spatial_value["incumbent_excess"]
                    assert abs(difference)/basal["full_excess_amount"] <= 1e-7
                    seed_status = "ratio_and_incumbent_collapse_coincident" if coincident_collapse else "ratio_reached_with_two_regions" if eligible else "ratio_reached_but_primary_contours_not_separate_or_amount_ineligible"
                else:
                    seed_status = "takeover_during_cue_incumbent_90_percent_excess_loss" if kind == "continuous_incumbent_90_percent_excess_loss" else kind
                rhs_archive = directory/"withdrawal_exact_rhs.npz"
                np.savez_compressed(rhs_archive, state=withdrawal, cue_free_rhs=adapter.rhs(withdrawal),
                    cue_on_rhs=adapter.rhs(withdrawal, placement["cue_peak"], placement["cue_center_fraction"]),
                    material_positions=adapter.q)
                record.update(seed_status=seed_status, B_distance_claim_eligible=eligible,
                    actual_withdrawal_ratio=spatial_value.get("challenger_to_incumbent_ratio"),
                    withdrawal_time=withdrawal_time, withdrawal_spatial=spatial_value,
                    measured_withdrawal_gap=gaps["0.1"], gap_threshold_sensitivity=gaps,
                    instantaneous_cue_free_exchange=regions.frozen_basin_exchange_budget(adapter, withdrawal, spatial_value),
                    exact_rhs_archive=str(rhs_archive.relative_to(run)), exact_rhs_sha256=sha(rhs_archive))
                tail = phase("cue_free_resolution_tail", cap, 0.)
                resolution = resolution_summary(adapter, tail["time"], tail["states"], placement)
                record.update(B_resolution=resolution if eligible else dict(status="unavailable_for_primary_two_region_distance_claim", reason=seed_status, raw_resolution=resolution),
                    final_winner=resolution["final_winner"], finite_tail=native.finite_tail(adapter, tail["time"], tail["states"]))
        arrays = trace_arrays(adapter, np.asarray(history_t), np.asarray(history_y), placement)
        archive = directory/"complete_region_readouts.npz"
        np.savez_compressed(archive, **arrays)
        record.update(status="complete", complete=True, readouts=str(archive.relative_to(run)), readouts_sha256=sha(archive),
            completed_native_horizon=elapsed, native_field_fate_separate_from_finite_tail=True)
    except Exception as exc:
        record.update(status="failed", complete=True, error_type=type(exc).__name__, error=str(exc))
    record["generated_at"] = now()
    write(statuspath, record)
    print(json.dumps(dict(slot=slot, status=record["status"], arm=job["arm"], seed_status=record.get("seed_status"))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "worker"))
    parser.add_argument("--arm", choices=("A2", "B"), default="B")
    parser.add_argument("--run", type=Path)
    parser.add_argument("--original", type=Path, default=original_campaign.DEFAULT)
    parser.add_argument("--refinement", type=Path, default=REFINEMENT)
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    run = (args.run or (DEFAULT_B if args.arm == "B" else DEFAULT_A2)).resolve()
    result = freeze(run, args.arm, args.original.resolve(), args.refinement.resolve()) if args.action == "freeze" else worker(run, args.slot)
    if result is not None:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
