"""Bounded paired lower-sigma forcing diagnostics for unresolved E1 flags.

Only explicitly listed ambiguous conditions are prepared. Twelve original seed
indices are paired at sigma0.25/0.5; native parameters, means and clocks remain
unchanged. These trajectories never inflate the sigma1 scientific seed count.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time
import traceback
import figure4_round3_campaign_20261022 as C


def prepare(campaign, conditions):
    campaign = Path(campaign).resolve()
    parent = C.verify(campaign)
    with Path(conditions).open(newline="") as stream:
        requested = {r["condition_id"] for r in csv.DictReader(stream)}
    chosen = [c for c in parent["conditions"] if c["condition_id"] in requested]
    if len(chosen) != len(requested) or not chosen:
        raise ValueError("Supply a nonempty list of exact ambiguous designed conditions")
    out = campaign / "lower_sigma_controls"
    if out.exists():
        raise ValueError("Lower-sigma cohort is already frozen; preserve its exact requested set")
    for name in ("workers", "logs", "results", "frozen_control_sources"):
        (out / name).mkdir(parents=True)
    (out / "inputs").symlink_to(campaign / "inputs", target_is_directory=True)
    (out / "frozen_sources").symlink_to(campaign / "frozen_sources", target_is_directory=True)
    sources = copy.deepcopy(parent["source_sha256"])
    for name in (Path(__file__).name, Path(C.__file__).name, "signature_figure4_regions_20261020.py"):
        target = out / "frozen_control_sources" / name
        shutil.copy2(C.CODE / name, target)
        sources[str(target.relative_to(out))] = C.sha(target)
    jobs = [dict(condition=c, sigma=sigma, replicate=replicate,
        seed=202610220000000 + c["condition_index"] * 1000 + replicate)
        for c in chosen for sigma in (.25, .5) for replicate in range(12)]
    m = dict(created_utc=C.now(), campaign="Figure4_round3_paired_E1_diagnostics", jobs=jobs,
        conditions=chosen, native_protocols=parent["native_protocols"], source_sha256=sources,
        input_sha256=parent["input_sha256"], parent_manifest_sha256=C.sha(campaign / "manifest.json"),
        requested_conditions=str(Path(conditions).resolve()), requested_conditions_sha256=C.sha(conditions),
        sigma1_independent_seed_count_added=0, pairing="Same first12 Gaussian paths and exact qualified basal as the fresh sigma1 cohort", native_parameters_retuned=False)
    C.write(out / "manifest.json", m)
    C.write(out / "freeze.json", dict(manifest_sha256=C.sha(out / "manifest.json")))
    script = out / "worker.sh"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " + shlex.join([C.PYTHON, "-B", str(out / "frozen_control_sources" / Path(__file__).name), "worker", "--campaign", str(campaign)]) + ' --task "$SLURM_ARRAY_TASK_ID"\n')
    script.chmod(0o750)
    return dict(conditions=len(chosen), diagnostic_trajectories=len(jobs), out=str(out))


def worker(campaign, task):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Native sigma diagnostics require SLURM")
    import numpy as np
    campaign = Path(campaign).resolve()
    out = campaign / "lower_sigma_controls"
    m = C.verify(out)
    job = m["jobs"][task]
    c, sigma, replicate = job["condition"], job["sigma"], job["replicate"]
    folder = out / "workers" / f'c{c["condition_index"]:03d}' / f"sigma{sigma:g}" / f"seed{replicate:03d}"
    folder.mkdir(parents=True, exist_ok=True)
    status = dict(generated_utc=C.now(), complete=False, status="RUNNING", job=job)
    C.write(folder / "status.json", status)
    start = time.monotonic()
    try:
        t, p, circumference, basal, checks, audit = C.solve(out, c, job["seed"], sigma=sigma)
        original = campaign / "workers" / f'c{c["condition_index"]:03d}' / f"seed{replicate:03d}" / "noise_audit.json"
        if not original.exists() or C.read(original)["gaussian_knots_sha256"] != audit["gaussian_knots_sha256"]:
            raise ValueError("Lower-sigma Gaussian path does not exactly match its sigma1 companion")
        fields, protocol, summary = C.extract(t, p, circumference, basal)
        theta = 2 * np.pi * (np.arange(len(basal)) + .5) / len(basal)
        initial_M1 = c["basal_circumference"] * abs(np.mean(basal * np.exp(1j * theta)))
        fields["R"] = (np.hypot(fields["M1_real"], fields["M1_imag"]) / initial_M1).astype("float32")
        np.savez_compressed(folder / "compact.npz", **fields)
        C.write(folder / "noise_audit.json", audit)
        C.write(folder / "region_audit.json", dict(protocol=protocol, summary=summary))
        status.update(complete=True, status="COMPLETE", compact_sha256=C.sha(folder / "compact.npz"), validation=checks, exact_sigma1_gaussian_path=True)
    except Exception as error:
        status.update(status="FAILED", error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
        raise
    finally:
        status.update(generated_utc=C.now(), wall_seconds=time.monotonic() - start)
        C.write(folder / "status.json", status)
    return status


def submit(campaign):
    out = Path(campaign).resolve() / "lower_sigma_controls"
    m = C.verify(out)
    path = out / "results/submission.json"
    if path.exists() and C.read(path).get("submitted"):
        return C.read(path)
    command = ["sbatch", "--parsable", "--partition=red", "--qos=xxlarge", "--ntasks=1", "--cpus-per-task=1", "--mem=4G", "--time=02:00:00",
        f'--array=0-{len(m["jobs"])-1}', "--job-name=F4r3-E1", "--output=" + str(out / "logs/%A_%a.out"), "--error=" + str(out / "logs/%A_%a.err"), str(out / "worker.sh")]
    p = subprocess.run(command, text=True, capture_output=True, timeout=45)
    result = dict(generated_utc=C.now(), command=shlex.join(command), returncode=p.returncode, stdout=p.stdout, stderr=p.stderr, submitted=p.returncode == 0)
    C.write(path, result)
    return result


def collect(campaign):
    out = Path(campaign).resolve() / "lower_sigma_controls"
    m = C.verify(out)
    rows, incomplete = [], []
    for job in m["jobs"]:
        c, sigma, replicate = job["condition"], job["sigma"], job["replicate"]
        folder = out / "workers" / f'c{c["condition_index"]:03d}' / f"sigma{sigma:g}" / f"seed{replicate:03d}"
        path = folder / "status.json"
        if not path.exists() or not C.read(path).get("complete"):
            incomplete.append(job)
            continue
        s = C.read(path)
        if C.sha(folder / "compact.npz") != s["compact_sha256"]:
            raise ValueError("Diagnostic compact archive changed")
        rows.append(dict(condition_id=c["condition_id"], model=c["model"], size=c["size"], fraction=c["fraction"], sigma=sigma,
            replicate=replicate, seed=job["seed"], compact=str(folder / "compact.npz"), compact_sha256=s["compact_sha256"]))
    fields = ["condition_id", "model", "size", "fraction", "sigma", "replicate", "seed", "compact", "compact_sha256"]
    with (out / "results/compact_sigma_index.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = dict(generated_utc=C.now(), complete=not incomplete, complete_trajectories=len(rows), incomplete_jobs=incomplete,
        index=str(out / "results/compact_sigma_index.csv"))
    C.write(out / "status.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker", "submit", "collect"))
    parser.add_argument("--campaign", type=Path, default=C.DEFAULT)
    parser.add_argument("--conditions", type=Path)
    parser.add_argument("--task", type=int)
    args = parser.parse_args()
    result = prepare(args.campaign, args.conditions) if args.action == "prepare" else worker(args.campaign, args.task) if args.action == "worker" else submit(args.campaign) if args.action == "submit" else collect(args.campaign)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
