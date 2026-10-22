"""Refine complete-region A1 handovers from exact retained native brackets.

The short native replay starts at the saved lower state, resolves EC=EI with
the native dense BDF interpolant, and continues to the original upper sample
to verify full-field agreement. Original A1 evidence remains immutable.
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
import signature_figure3_campaign_20261020 as campaign
import signature_figure3_regions_20261020 as regions

ROOT = native.ROOT
DEFAULT = ROOT/"data/runs/signature_figure3_handover_refinement_20261020"
PROPOSALS = ROOT/"data/runs/signature_figure3_region_proposal_20261020"
read, write, sha, now = native.read, native.write, native.sha, native.now


def freeze(run, original, proposal, stage="primary"):
    if run.exists():
        raise FileExistsError("Fresh immutable refinement run required")
    primary = campaign.verify(original)
    proposal_folder = proposal/"results" if stage == "primary" else proposal
    summary = read(proposal_folder/"last_reanalysis.json")
    intended = primary["stages"][stage]
    assert summary["cases"] == len(intended) and summary["source_manifest_sha256"] == sha(original/"manifest.json")
    for folder in ("inputs", "frozen_sources", "workers", "results", "logs"):
        (run/folder).mkdir(parents=True)
    sources, inputs = {}, {}
    for relative, digest in primary["source_sha256"].items():
        source, target = original/"frozen_sources"/relative, run/"frozen_sources"/relative
        assert sha(source) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        sources[relative] = digest
    for source in (Path(__file__), Path(regions.__file__)):
        relative = str(source.relative_to(ROOT))
        target = run/"frozen_sources"/relative
        shutil.copy2(source, target)
        sources[relative] = sha(source)
    jobs, unavailable = [], []
    for job in intended:
        statuspath = original/"workers"/stage/f'{job["slot"]:04d}'/"status.json"
        status = read(statuspath)
        placement = status["placement"]
        path = proposal_folder/f'{stage}_{job["slot"]:04d}_placement_{placement["placement_id"]:03d}_region_proposal.json'
        record = read(path)
        bracket = record["complete_region_handover_first_bracket"]
        if bracket is None:
            unavailable.append(dict(original_slot=job["slot"], placement_id=placement["placement_id"],
                model=placement["model"], size_parameter=placement["size_parameter"],
                classification=placement["classification_10"], control=job.get("control"),
                status=record["complete_region_handover_status"], source_proposal=str(path),
                source_proposal_sha256=sha(path), final_state_class=record["final_state_class"]))
            continue
        saved = {}
        for phase in status["phases"]:
            archive = original/phase["path"]
            assert sha(archive) == phase["sha256"]
            with np.load(archive) as z:
                times = phase["absolute_start"]+z["time"]
                for name in ("lower", "upper"):
                    indices = np.flatnonzero(abs(times-bracket[name]) < 1e-8)
                    if len(indices):
                        saved[name] = z["states"][int(indices[-1])].copy()
        assert set(saved) == {"lower", "upper"}
        target = run/"inputs"/f'original_slot_{job["slot"]:04d}_bracket.npz'
        np.savez_compressed(target, state=saved["lower"], original_upper_state=saved["upper"],
                            original_lower_time=bracket["lower"], original_upper_time=bracket["upper"])
        inputs[str(target.relative_to(run))] = sha(target)
        jobs.append(dict(slot=len(jobs), original_slot=job["slot"], placement=placement,
            control=job.get("control"), bracket=bracket, input=str(target.relative_to(run)),
            input_sha256=sha(target), source_status=str(statuspath), source_status_sha256=sha(statuspath),
            source_proposal=str(path), source_proposal_sha256=sha(path),
            basal_excess=record["initial_complete_basal_excess"]))
    manifest = dict(generated_at=now(), experiment="Native complete-region A1 handover bracket refinement",
        source_run=str(original), source_manifest_sha256=sha(original/"manifest.json"),
        proposal_run=str(proposal), proposal_manifest_sha256=sha(proposal/"manifest.json") if (proposal/"manifest.json").exists() else None,
        proposal_reanalysis_record_sha256=sha(proposal_folder/"last_reanalysis.json"),
        source_sha256=sources, input_sha256=inputs, jobs=jobs, unavailable=unavailable,
        original_cases=len(intended), original_stage=stage, worker_soft_seconds=1680., region_definition_authority="User: yes obviously you should not restrict challenger activity counting to the area under the cue.",
        event="First continuous two-region EC=EI crossing in original saved bracket; no cap or topology substitution",
        replay_upper_field_tolerance=1e-6, native_kinetics_changed=False)
    write(run/"manifest.json", manifest)
    write(run/"freeze.json", dict(manifest_sha256=sha(run/"manifest.json"), generated_at=now()))
    runner = run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command = "exec "+shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", native.PYTHON, str(runner), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    payload = dict(generated_at=now(), submitted=False, jobs=len(jobs), sbatch_arguments=["sbatch", "--parsable",
        "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=2G", "--time=00:30:00", "--signal=B:USR1@60",
        f"--array=0-{len(jobs)-1}", "--job-name=figure3-handover-refine", "--output="+str(run/"logs/%A_%a.out"),
        "--error="+str(run/"logs/%A_%a.err"), "--wrap", command], no_concurrency_throttle=True)
    write(run/"slurm_payload.json", payload)
    return dict(jobs=len(jobs), unavailable=len(unavailable), primary_jobs=sum(not j.get("control") for j in jobs),
        mirror_jobs=sum(j.get("control") == "mirror_placement" for j in jobs),
        control_jobs=sum(bool(j.get("control")) for j in jobs), native_integrations=0)


def worker(run, slot):
    native.require_slurm()
    manifest = campaign.verify(run)
    job = manifest["jobs"][slot]
    placement = job["placement"]
    adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
    directory = run/"workers"/f"{slot:04d}"
    directory.mkdir(parents=True)
    with np.load(run/job["input"]) as z:
        initial, reference_upper = z["state"].copy(), z["original_upper_state"].copy()
    old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
    cue = placement["cue_center_fraction"]
    def observable(state):
        result = regions.readout(adapter.primary(state), adapter.circumference(state), old, cue, detailed=False)
        if not result["two_region_available"]:
            raise ValueError("Topology invalid inside original two-region handover bracket")
        return (result["challenger_excess"]-result["incumbent_excess"])/job["basal_excess"]
    statuspath = directory/"status.json"
    record = dict(generated_at=now(), status="running", complete=False, job=job,
        slurm_job_id=os.environ["SLURM_JOB_ID"], phases=[])
    write(statuspath, record)
    try:
        lower, upper = job["bracket"]["lower"], job["bracket"]["upper"]
        assert observable(initial) < 0 and observable(reference_upper) >= 0
        deadline = time.monotonic()+manifest["worker_soft_seconds"]
        result = native.integrate(adapter, initial, upper-lower, placement["cue_peak"], cue,
                                  event=observable, deadline=deadline)
        archive = directory/"000_refine_handover.npz"
        phase = native.save_phase(archive, result, adapter, placement["cue_peak"], cue)
        phase.update(path=str(archive.relative_to(run)), absolute_start=lower, kind="refine_handover")
        record["phases"].append(phase)
        assert result["status"] == "event" and result["event_time"] is not None
        event_time = lower+result["event_time"]
        state = result["states"][-1]
        assert lower <= event_time <= upper
        root_error = abs(observable(state))
        assert root_error <= 1e-7
        remaining = upper-event_time
        if remaining > 1e-8:
            verify = native.integrate(adapter, state, remaining, placement["cue_peak"], cue, deadline=deadline)
            archive = directory/"001_verify_original_upper.npz"
            phase = native.save_phase(archive, verify, adapter, placement["cue_peak"], cue)
            phase.update(path=str(archive.relative_to(run)), absolute_start=event_time, kind="verify_original_upper")
            record["phases"].append(phase)
            assert verify["status"] == "complete"
            replay_upper = verify["states"][-1]
        else:
            replay_upper = state
        error = float(np.max(abs(replay_upper-reference_upper)/adapter.scales))
        assert error <= manifest["replay_upper_field_tolerance"], error
        event_readout = regions.readout(adapter.primary(state), adapter.circumference(state), old, cue)
        record.update(status="complete", complete=True, handover_time=event_time,
            continuous_two_region_root=True, root_excess_difference_normalized=root_error,
            original_saved_bracket=[lower, upper], upper_replay_maximum_scaled_field_difference=error,
            event_readout=event_readout)
    except Exception as exc:
        record.update(status="failed", complete=True, error_type=type(exc).__name__, error=str(exc))
    record["generated_at"] = now()
    write(statuspath, record)
    print(json.dumps(dict(slot=slot, status=record["status"], handover_time=record.get("handover_time"))))


def audit(run):
    manifest = campaign.verify(run)
    records = []
    for job in manifest["jobs"]:
        statuspath = run/"workers"/f'{job["slot"]:04d}'/"status.json"
        if not statuspath.exists():
            records.append(dict(slot=job["slot"], status="missing")); continue
        record = read(statuspath)
        if record["status"] != "complete":
            records.append(dict(slot=job["slot"], status=record["status"], error=record.get("error"))); continue
        adapter = native.Adapter(job["placement"]["model"], job["placement"]["size_parameter"], job["placement"]["cells"])
        with np.load(run/job["input"]) as z:
            previous, reference_upper = z["state"].copy(), z["original_upper_state"].copy()
        reference_amounts = adapter.amounts(previous)
        root_state = None
        for phase in record["phases"]:
            archive = run/phase["path"]
            assert sha(archive) == phase["sha256"]
            with np.load(archive) as z:
                times, states = z["time"].copy(), z["states"].copy()
            np.testing.assert_array_equal(states[0], previous)
            amounts = np.asarray([adapter.amounts(y) for y in states])
            assert np.max(abs(amounts/reference_amounts-1)) <= 1e-7
            assert min(adapter.concentration_minimum(y) for y in states) >= -5e-10
            previous = states[-1]
            if phase["kind"] == "refine_handover":
                root_state = previous.copy()
                np.testing.assert_allclose(phase["absolute_start"]+times[-1], record["handover_time"], atol=1e-8)
        lower, upper = record["original_saved_bracket"]
        assert lower <= record["handover_time"] <= upper
        error = float(np.max(abs(previous-reference_upper)/adapter.scales))
        assert error <= manifest["replay_upper_field_tolerance"]
        placement = job["placement"]
        spatial = regions.readout(adapter.primary(root_state), adapter.circumference(root_state),
            placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"], placement["cue_center_fraction"])
        assert spatial["two_region_available"]
        root_error = abs(spatial["challenger_excess"]-spatial["incumbent_excess"])/job["basal_excess"]
        assert root_error <= 1e-7
        records.append(dict(slot=job["slot"], original_slot=job["original_slot"], status="verified", control=job.get("control"),
            handover_time=record["handover_time"], status_path=str(statuspath), status_sha256=sha(statuspath),
            root_error=root_error, replay_upper_field_error=error))
    by_slot = {r["slot"]: r for r in records}
    primary = {j["placement"]["placement_id"]: by_slot[j["slot"]] for j in manifest["jobs"] if not j.get("control")}
    mirrors = []
    for job in manifest["jobs"]:
        if job.get("control") != "mirror_placement":
            continue
        reference, actual = primary.get(job["placement"]["placement_id"]), by_slot[job["slot"]]
        valid = reference is not None and reference["status"] == actual["status"] == "verified"
        error = abs(reference["handover_time"]-actual["handover_time"])/max(1., reference["handover_time"], actual["handover_time"]) if valid else None
        if valid:
            assert error <= 1e-5
        mirrors.append(dict(slot=job["slot"], placement_id=job["placement"]["placement_id"], verified=valid, relative_handover_difference=error))
    result = dict(generated_at=now(), requested=len(manifest["jobs"]), verified=sum(r["status"] == "verified" for r in records),
        complete=all(r["status"] == "verified" for r in records), records=records,
        mirror_checks=mirrors, all_mirrors_verified=all(v["verified"] for v in mirrors),
        source_manifest_sha256=sha(run/"manifest.json"), original_cases=manifest["original_cases"],
        unavailable_original_cases=len(manifest["unavailable"]),
        scope="Independent exact native saved fields, original bracket, root equality, full upper-endpoint replay agreement, conservation and positivity; no second trajectory replay")
    write(run/"results/independent_refinement_audit.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "worker", "audit"))
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--original", type=Path, default=campaign.DEFAULT)
    parser.add_argument("--proposal", type=Path, default=PROPOSALS)
    parser.add_argument("--slot", type=int)
    parser.add_argument("--stage", default="primary")
    args = parser.parse_args()
    run = args.run.resolve()
    result = freeze(run, args.original.resolve(), args.proposal.resolve(), args.stage) if args.action == "freeze" else worker(run, args.slot) if args.action == "worker" else audit(run)
    if result is not None:
        print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
