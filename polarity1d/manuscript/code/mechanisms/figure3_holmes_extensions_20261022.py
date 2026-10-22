"""Approved nine-case Holmes cue-free continuation, without changing producers.

Each worker starts at the exact archived 20,000 endpoint and integrates 60,000
more native units. Original withdrawals, complete-region definitions and solver
settings remain frozen. Source campaigns are read only.
"""
from pathlib import Path
import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import time

import numpy as np
import signature_figure3_native_20261020 as native
import signature_figure3_regional_campaign_20261020 as regional

ROOT = native.ROOT
DEFAULT = ROOT / "data/runs/signature_figure3_holmes_extensions_round2_20261022"
CANDIDATES = ROOT / "figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure3/optional_extension_candidates.csv"
read, write, sha, now = native.read, native.write, native.sha, native.now


def verify(run):
    plan = read(run / "manifest.json")
    assert sha(run / "manifest.json") == read(run / "freeze.json")["manifest_sha256"]
    for relative, digest in plan["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    for relative, digest in plan["input_sha256"].items():
        assert sha(run / relative) == digest
    return plan


def freeze(run, candidates):
    assert not run.exists(), "Fresh immutable campaign required; inspect existing campaign before restarting"
    with candidates.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 9
    original = Path(rows[0]["native_status"]).parents[2]
    source = read(original / "manifest.json")
    assert sha(original / "manifest.json") == read(original / "freeze.json")["manifest_sha256"]
    for folder in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run / folder).mkdir(parents=True)
    sources, inputs, jobs = {}, {}, []
    for relative, digest in source["source_sha256"].items():
        old, target = original / "frozen_sources" / relative, run / "frozen_sources" / relative
        assert sha(old) == digest
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(old, target)
        sources[relative] = digest
    relative = "code/mechanisms/" + Path(__file__).name
    shutil.copy2(__file__, run / "frozen_sources" / relative)
    sources[relative] = sha(run / "frozen_sources" / relative)
    shutil.copy2(candidates, run / "inputs/approved_candidates.csv")
    inputs["inputs/approved_candidates.csv"] = sha(run / "inputs/approved_candidates.csv")
    for row in rows:
        oldstatus = Path(row["native_status"])
        status = read(oldstatus)
        assert status["complete"] and status["status"] == "complete"
        assert status["B_distance_claim_eligible"] and status["target_ratio"] == 1
        assert status["placement"]["model"] == "holmes_model3"
        assert status["B_resolution"]["primary"]["upper"] is None
        tail = next(p for p in status["phases"] if p["kind"] == "cue_free_resolution_tail")
        oldtail = original / tail["path"]
        assert sha(oldtail) == tail["sha256"]
        with np.load(oldtail) as z:
            assert z["time"][-1] == 20000.
            endpoint = z["states"][-1].copy()
        slot = len(jobs)
        prefix = f"inputs/{slot:03d}"
        targetstatus = run / (prefix + "_original_status.json")
        shutil.copy2(oldstatus, targetstatus)
        endpointpath = run / (prefix + "_endpoint.npz")
        np.savez_compressed(endpointpath, state=endpoint, elapsed_after_withdrawal=20000.)
        for p in (targetstatus, endpointpath):
            inputs[str(p.relative_to(run))] = sha(p)
        jobs.append(dict(slot=slot, original_slot=int(row["slot"]), size=float(row["size"]),
            measured_gap=float(row["measured_gap"]), original_status=str(oldstatus),
            original_status_sha256=sha(oldstatus), status_input=str(targetstatus.relative_to(run)),
            endpoint_input=str(endpointpath.relative_to(run)), original_tail=str(oldtail),
            original_tail_sha256=sha(oldtail), original_horizon=20000., final_horizon=80000.,
            additional_duration=60000.))
    plan = dict(generated_utc=now(), experiment="Author-approved Figure 3 nine Holmes cue-free extensions",
        approval="manuscript_outline/revision_round2_20261022/urgent_user_steering_20261022T110301Z.txt item 2",
        original_campaign=str(original), original_manifest_sha256=sha(original / "manifest.json"),
        source_sha256=sources, input_sha256=inputs, jobs=jobs, worker_soft_seconds=3480.,
        unchanged_native_kinetics=True, continuation="Exact native endpoint at 20,000; cue remains zero; additional 60,000",
        solver=dict(method="BDF", rtol=1e-9, atol=1e-11, max_step=10.),
        resolution="Original complete-region single resolved peak; 0.05 prominence and contrast gates; sampled interval",
        no_cap_substitution=True, source_campaign_immutable=True)
    write(run / "manifest.json", plan)
    write(run / "freeze.json", dict(generated_utc=now(), manifest_sha256=sha(run / "manifest.json")))
    runner = run / "frozen_sources" / relative
    command = "exec " + shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1", native.PYTHON, str(runner), "worker", "--run", str(run)]) + ' --slot "$SLURM_ARRAY_TASK_ID"'
    arguments = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=4G",
        "--time=01:00:00", "--signal=B:USR1@60", "--array=0-8", "--job-name=figure3-Holmes-80k",
        "--output=" + str(run / "logs/%A_%a.out"), "--error=" + str(run / "logs/%A_%a.err"), "--wrap", command]
    write(run / "slurm_payload.json", dict(generated_utc=now(), sbatch_arguments=arguments, no_concurrency_throttle=True,
        resource_reason="Nine serial CPU-only sparse BDF integrations; existing small QoS suits one CPU/4GB/one hour"))
    (run / "submission_command.txt").write_text(shlex.join(arguments) + "\n")
    verify(run)
    return dict(frozen=True, jobs=len(jobs), run=str(run))


def worker(run, slot):
    native.require_slurm()
    plan = verify(run)
    job = plan["jobs"][slot]
    original = read(run / job["status_input"])
    placement = original["placement"]
    adapter = native.Adapter("holmes_model3", job["size"], placement["cells"])
    directory = run / "workers" / f"{slot:03d}"
    directory.mkdir()
    statuspath = directory / "status.json"
    record = dict(generated_utc=now(), status="running", complete=False, job=job,
        slurm_job_id=os.environ["SLURM_JOB_ID"], placement=placement)
    write(statuspath, record)
    try:
        with np.load(run / job["endpoint_input"]) as z:
            initial = z["state"].copy()
        result = native.integrate(adapter, initial, job["additional_duration"], peak=0.,
            center_fraction=placement["cue_center_fraction"], deadline=time.monotonic()+plan["worker_soft_seconds"])
        archive = directory / "continuation_full_native_fields.npz"
        phase = native.save_phase(archive, result, adapter, 0., placement["cue_center_fraction"])
        phase.update(path=str(archive.relative_to(run)), absolute_start_after_withdrawal=20000., kind="cue_free_resolution_extension")
        record["phase"] = phase
        write(statuspath, record)
        assert result["status"] == "complete", result["status"]
        assert result["time"][-1] == 60000.
        assert sha(Path(job["original_tail"])) == job["original_tail_sha256"]
        with np.load(job["original_tail"]) as z:
            oldt, oldy = z["time"].copy(), z["states"].copy()
        np.testing.assert_array_equal(oldy[-1], result["states"][0])
        times = np.r_[oldt, 20000. + result["time"][1:]]
        states = np.concatenate((oldy, result["states"][1:]))
        resolution = regional.resolution_summary(adapter, times, states, placement)
        trace = directory / "complete_region_readouts.npz"
        np.savez_compressed(trace, **regional.trace_arrays(adapter, times, states, placement))
        record.update(status="complete", complete=True, completed_horizon_after_withdrawal=80000.,
            exact_endpoint_continuity=True, B_resolution=resolution,
            final_winner=resolution["final_winner"], finite_tail=native.finite_tail(adapter, times, states),
            readouts=str(trace.relative_to(run)), readouts_sha256=sha(trace))
    except Exception as exc:
        record.update(status="failed", complete=True, error_type=type(exc).__name__, error=str(exc))
    record["generated_utc"] = now()
    write(statuspath, record)
    return dict(slot=slot, status=record["status"], resolution=record.get("B_resolution", {}).get("primary"))


def submit(run):
    verify(run)
    assert not (run / "submission.json").exists(), "Already submitted; inspect authoritative job state"
    payload = read(run / "slurm_payload.json")
    done = subprocess.run(payload["sbatch_arguments"], capture_output=True, text=True, check=True)
    receipt = dict(submitted_utc=now(), job_id=done.stdout.strip().split(";")[0], stdout=done.stdout, stderr=done.stderr)
    write(run / "submission.json", receipt)
    return receipt


def collect(run):
    plan = verify(run)
    records, rows = [], []
    for job in plan["jobs"]:
        path = run / "workers" / f'{job["slot"]:03d}' / "status.json"
        assert path.exists(), str(path)
        record = read(path)
        assert record["complete"] and record["status"] == "complete", str(path)
        assert record["completed_horizon_after_withdrawal"] == 80000.
        for k in ("original_status", "original_tail"):
            assert sha(Path(job[k])) == job[k+"_sha256"]
        for name in ("phase",):
            assert sha(run / record[name]["path"]) == record[name]["sha256"]
        assert sha(run / record["readouts"]) == record["readouts_sha256"]
        primary = record["B_resolution"]["primary"]
        rows.append(dict(slot=job["original_slot"], extension_slot=job["slot"], model="holmes_model3", size=job["size"],
            measured_gap=job["measured_gap"], resolution_status=primary["status"], resolution_lower=primary["lower"],
            resolution_upper=primary["upper"], extension_status=str(path), extension_status_sha256=sha(path),
            completed_horizon=80000., final_state_class=record["B_resolution"]["final_state_class"],
            finite_tail_passed=record["finite_tail"]["both_late_windows_steady"]))
        records.append(dict(slot=job["slot"], original_slot=job["original_slot"], status_path=str(path), status_sha256=sha(path)))
    with (run / "results/extension_evidence.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    report = dict(generated_utc=now(), complete=True, verified=9, source_manifest_sha256=sha(run / "manifest.json"),
        records=records, resolved=sum(r["resolution_upper"] is not None for r in rows),
        censored=sum(r["resolution_upper"] is None for r in rows), original_source_campaign_unchanged=True)
    write(run / "results/collection.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "verify", "submit", "worker", "collect"))
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--candidates", type=Path, default=CANDIDATES)
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    run = args.run.resolve()
    result = freeze(run, args.candidates) if args.action == "freeze" else worker(run, args.slot) if args.action == "worker" else globals()[args.action](run)
    print(json.dumps(result, indent=2))
