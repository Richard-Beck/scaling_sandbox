"""SLURM-dependent compact collection, adaptive precision and final rebuild.

The controller never claims that a queued campaign has completed. An afterany
dependency lets it record failed workers as well as successful campaigns.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

import figure4_round3_campaign_20261022 as C

ESTIMATOR = C.CODE / "figure4_round3_estimators_20261022.py"
NULL_RECEIPT = C.MS / "figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4/supporting/metadata/null_analysis_receipt.json"


def run_command(command, receipt):
    try:
        p = subprocess.run(command, text=True, capture_output=True)
        record = dict(command=shlex.join(command), returncode=p.returncode,
                      stdout=p.stdout, stderr=p.stderr, generated_utc=C.now())
    except OSError as error:
        record = dict(command=shlex.join(command), returncode=None, error=str(error), generated_utc=C.now())
    C.write(receipt, record)
    return record


def nulls_ready():
    if not NULL_RECEIPT.exists():
        return False
    receipt = C.read(NULL_RECEIPT)
    return "A" in receipt and receipt.get("B", {}).get("conditions") == 615


def finalize(out, rebuild_command):
    """A lightweight scheduled recheck avoids depending on a24h allocation."""
    out = Path(out).resolve()
    if not nulls_ready():
        stamp = str(time.time_ns())
        script = out / "submissions" / f"wait_nulls_{stamp}.sh"
        command = [C.PYTHON, "-B", str(Path(__file__).resolve()), "finalize", "--out", str(out), "--rebuild-command", rebuild_command]
        script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " + shlex.join(command) + "\n")
        script.chmod(0o750)
        submit = ["sbatch", "--parsable", "--partition=red", "--qos=small", "--ntasks=1", "--cpus-per-task=4", "--mem=16G", "--time=04:00:00", "--begin=now+5minutes", "--job-name=F4r3-wait-nullB",
            "--output=" + str(out / "logs/finalize_%j.out"), "--error=" + str(out / "logs/finalize_%j.err")]
        job_id = os.environ.get("SLURM_JOB_ID")
        if job_id and job_id.isdigit():
            submit += ["--dependency=afterany:" + job_id]
        submit += [str(script)]
        receipt = run_command(submit, out / "submissions" / f"wait_nulls_{stamp}.json")
        if receipt["returncode"] != 0:
            raise RuntimeError("Null-B recheck failed to submit; exact receipt retained")
        result = dict(generated_utc=C.now(), status="WAITING_NULL_B", new_compact_complete=True,
            next_job_id=receipt["stdout"].strip(), null_receipt=str(NULL_RECEIPT),
            reason="All615 simulated null-B conditions must be recorded before final figure rebuild; recheck in5minutes.")
        C.write(out / "results/pipeline_status.json", result)
        return result
    receipt = run_command(shlex.split(rebuild_command), out / "results/final_rebuild_execution.json")
    if receipt["returncode"] != 0:
        raise RuntimeError("Final figure rebuild failed; exact execution receipt retained")
    result = C.read(out / "results/pipeline_status.json") if (out / "results/pipeline_status.json").exists() else {}
    result.update(generated_utc=C.now(), status="COMPLETE", final_rebuild=receipt,
        null_receipt_sha256=C.sha(NULL_RECEIPT), old_native_deletion_performed=False)
    C.write(out / "results/pipeline_status.json", result)
    C.update_storage(out)
    return result


def dependency_job(out, jobs, batch, rebuild_command=None, phase=""):
    out = Path(out).resolve()
    jobs = list(dict.fromkeys(str(j).strip().split(";")[0] for j in jobs))
    if not jobs or any(not j.isdigit() for j in jobs):
        raise ValueError("Dependent stage requires actual numeric SLURM job IDs")
    suffix = "_" + phase if phase else ""
    script = out / "submissions" / f"pipeline_stage_batch{batch}{suffix}.sh"
    args = [C.PYTHON, "-B", str(Path(__file__).resolve()), "stage", "--out", str(out), "--batch", str(batch)]
    if rebuild_command:
        args += ["--rebuild-command", rebuild_command]
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " + shlex.join(args) + "\n")
    script.chmod(0o750)
    command = ["sbatch", "--parsable", "--partition=red", "--qos=small", "--ntasks=1", "--cpus-per-task=4", "--mem=16G", "--time=04:00:00",
        "--dependency=afterany:" + ":".join(jobs), "--job-name=F4r3-analysis-b" + str(batch) + suffix,
        "--output=" + str(out / "logs/pipeline_%j.out"), "--error=" + str(out / "logs/pipeline_%j.err"), str(script)]
    receipt = run_command(command, out / "submissions" / f"pipeline_batch{batch}{suffix}.json")
    if receipt["returncode"] != 0:
        raise RuntimeError("Dependent SLURM stage failed to submit: " + receipt.get("stderr", receipt.get("error", "")))
    return dict(job_id=receipt["stdout"].strip().split(";")[0], dependencies=jobs, batch=batch)


def launch(out, rebuild_command=None, additional_dependencies=()):
    out = Path(out).resolve()
    if not out.exists():
        C.prepare(out)
    prior = out / "results/submission_summary.json"
    if prior.exists():
        records = C.read(prior)["records"]
        if any(r.get("submitted") for r in records):
            if not all(r.get("submitted") for r in records):
                raise RuntimeError("Partial initial submission needs targeted review; refusing duplicate seeds")
        else:
            records = C.submit(out)
    else:
        records = C.submit(out)
    if not all(r.get("submitted") for r in records):
        raise RuntimeError("Native submission failed; exact receipts and Richard command are retained")
    existing = out / "submissions/pipeline_batch0.json"
    if existing.exists() and C.read(existing).get("returncode") == 0:
        return dict(status="already_scheduled", receipt=str(existing))
    dependencies = [r["stdout"] for r in records] + list(additional_dependencies)
    control_submission = out / "matched_mean_controls/results/submission.json"
    if control_submission.exists() and C.read(control_submission).get("returncode") == 0:
        dependencies.append(C.read(control_submission)["stdout"])
    return dependency_job(out, dependencies, 0, rebuild_command)


def stage(out, batch, rebuild_command=None):
    out = Path(out).resolve()
    current = C.collect(out)
    if current["missing_seeds"] or current["failed_seeds"] or current["running_seeds"]:
        result = dict(generated_utc=C.now(), status="NATIVE_INCOMPLETE", batch=batch,
            missing=current["missing_seeds"], running=current["running_seeds"], failures=current["failed_seeds"],
            reason="Analysis cannot replace missing/failed seeds. Old native arrays remain preserved.")
        C.write(out / "results/pipeline_status.json", result)
        return result
    precision = out / "results" / f"precision_batch{batch}.json"
    receipt = run_command([C.PYTHON, "-B", str(ESTIMATOR), "precision", "--input", str(out), "--output", str(precision)],
        out / "results" / f"precision_batch{batch}_execution.json")
    if receipt["returncode"] != 0:
        raise RuntimeError("Whole-seed precision analysis failed; see retained execution receipt")
    report = C.read(precision)
    rows = report["conditions"]
    if len(rows) != 260 or len({r["condition_id"] for r in rows}) != 260:
        raise ValueError("Precision analysis must account for every designed condition")
    need = [r for r in rows if r["n_seeds"] < 384 and (
        (r.get("se_log_k") is not None and math.isfinite(r["se_log_k"]) and r["se_log_k"] > .07)
        or r.get("reason") == "bootstrap_valid_fraction_below_0.975")]
    if need:
        if batch >= 3:
            raise RuntimeError("Adaptive batch state inconsistent with384seed maximum")
        native = C.submit(out, batch + 1, precision)
        if not native or not all(r.get("submitted") for r in native):
            raise RuntimeError("Adaptive native submission failed; receipts retained")
        next_job = dependency_job(out, [r["stdout"] for r in native], batch + 1, rebuild_command)
        result = dict(generated_utc=C.now(), status="ADAPTIVE_SUBMITTED", completed_batch=batch,
            selected_conditions=len(need), added_seeds=sum(r["seeds"] for r in native), next_stage=next_job,
            precision=str(precision), undefined_conditions=sum(r.get("se_log_k") is None for r in rows),
            adaptation="Finite whole-seed SE(logk)>.07 or an eligible point lag with<97.5% finite bootstrap logk triggers+96. No eligible point lag remains explicitly undefined without a topup.")
    else:
        final = C.collect(out, precision)
        if not final["complete"]:
            raise RuntimeError("Final cohort failed compact campaign completion audit")
        controls = None
        if (out / "matched_mean_controls/manifest.json").exists():
            import figure4_round3_controls_20261022 as controls_module
            controls = controls_module.collect(out)
            if not controls["complete"]:
                raise RuntimeError("Matched-mean control cohort is incomplete; refusing final figure rebuild")
        candidates_path = out / "results/lower_sigma_candidates.csv"
        candidate_receipt = run_command([C.PYTHON, "-B", str(C.CODE / "figure4_round3_analysis_20261022.py"),
            "--campaign", str(out), "--forcing-candidates", str(candidates_path)], out / "results/forcing_candidate_execution.json")
        if candidate_receipt["returncode"] != 0:
            raise RuntimeError("Forcing candidate classification failed; refusing final figure rebuild")
        with candidates_path.open(newline="") as stream:
            candidates = list(csv.DictReader(stream))
        lower_sigma = None
        if candidates:
            import figure4_round3_sigma_controls_20261022 as sigma_module
            if not (out / "lower_sigma_controls/manifest.json").exists():
                sigma_module.prepare(out, candidates_path)
                native = sigma_module.submit(out)
                if not native.get("submitted"):
                    raise RuntimeError("Targeted lower-sigma diagnostic submission failed")
                next_job = dependency_job(out, [native["stdout"]], batch, rebuild_command, phase="post_forcing")
                result = dict(generated_utc=C.now(), status="TARGETED_E1_DIAGNOSTICS_SUBMITTED", batch=batch,
                    compact_seeds=final["native_seeds_complete"], ambiguous_conditions=len(candidates),
                    paired_trajectories=24 * len(candidates), next_stage=next_job,
                    qualification="Only same-abnormal matched-mean/sigma1 conditions with a false mean baseline trigger12pairedseeds each at sigma0.25/0.5. Original sigma1 seed counts remain unchanged.")
                C.write(out / "results/pipeline_status.json", result)
                C.update_storage(out)
                return result
            lower_sigma = sigma_module.collect(out)
            if not lower_sigma["complete"]:
                raise RuntimeError("Targeted lower-sigma controls remain incomplete; refusing final figure rebuild")
        rebuild = None
        if rebuild_command:
            if not nulls_ready():
                return finalize(out, rebuild_command)
            rebuild = run_command(shlex.split(rebuild_command), out / "results/final_rebuild_execution.json")
            if rebuild["returncode"] != 0:
                raise RuntimeError("Compact campaign COMPLETE but final figure rebuild failed; receipt retained")
        result = dict(generated_utc=C.now(), status="COMPLETE", batch=batch,
            compact_seeds=final["native_seeds_complete"], precision=str(precision),
            undefined_conditions=final["undefined_rate_conditions"], final_rebuild=rebuild,
            matched_mean_controls=controls,
            lower_sigma_controls=lower_sigma,
            old_native_deletion_performed=False,
            optional_extensions="Not automatically run; largest weak-input undefined conditions remain listed for a bounded extension if needed.")
    C.write(out / "results/pipeline_status.json", result)
    C.update_storage(out)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    for name in ("launch", "stage", "finalize"):
        p = sub.add_parser(name)
        p.add_argument("--out", type=Path, default=C.DEFAULT)
        p.add_argument("--rebuild-command")
        if name == "launch":
            p.add_argument("--additional-dependency", action="append", default=[])
        if name == "stage":
            p.add_argument("--batch", type=int, required=True)
    args = parser.parse_args()
    if args.action == "launch":
        result = launch(args.out, args.rebuild_command, args.additional_dependency)
    elif args.action == "stage":
        result = stage(args.out, args.batch, args.rebuild_command)
    else:
        if not args.rebuild_command:
            parser.error("finalize requires --rebuild-command")
        result = finalize(args.out, args.rebuild_command)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
