"""Exact matched-mean controls for newly designed Figure4 conditions.

One deterministic compact trajectory per missing condition; this is a forcing
diagnostic, not an additional independent stochastic seed or model retuning.
"""
from __future__ import annotations
import argparse
import copy
import csv
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import subprocess
import sys
import time
import traceback

import figure4_round3_campaign_20261022 as C


def prepare(campaign, existing=None, missing=None):
    campaign = Path(campaign).resolve()
    parent = C.verify(campaign)
    out = campaign / "matched_mean_controls"
    if out.exists():
        return dict(out=str(out), already_prepared=True)
    seen = set()
    if existing:
        with Path(existing).open(newline="") as stream:
            for row in csv.DictReader(stream):
                if row.get("exact_native_mean_control", "true").lower() not in ("true", "1"):
                    continue
                seen.add(row.get("condition_id", f'{row["model"]}:C{float(row.get("size", row.get("C_label"))):g}:f{float(row.get("fraction", row.get("f"))):g}'))
    if missing:
        with Path(missing).open(newline="") as stream:
            requested = {r["condition_id"] for r in csv.DictReader(stream)}
        conditions = [c for c in parent["conditions"] if c["condition_id"] in requested]
        if len(conditions) != len(requested):
            raise ValueError("Missing-control list contains an unknown designed condition")
        seen = {c["condition_id"] for c in parent["conditions"] if c["condition_id"] not in requested}
    else:
        conditions = [c for c in parent["conditions"] if c["condition_id"] not in seen]
    for name in ("workers", "logs", "results", "frozen_control_sources"):
        (out / name).mkdir(parents=True)
    (out / "frozen_sources").symlink_to(campaign / "frozen_sources", target_is_directory=True)
    (out / "inputs").symlink_to(campaign / "inputs", target_is_directory=True)
    sources = copy.deepcopy(parent["source_sha256"])
    for name in (Path(__file__).name, Path(C.__file__).name, "signature_figure4_regions_20261020.py"):
        source = C.CODE / name
        target = out / "frozen_control_sources" / name
        shutil.copy2(source, target)
        sources[str(target.relative_to(out))] = C.sha(target)
    manifest = dict(created_utc=C.now(), campaign="Figure4_round3_missing_matched_mean_controls", conditions=conditions,
        parent_campaign=str(campaign), parent_manifest_sha256=C.sha(campaign / "manifest.json"),
        native_protocols=parent["native_protocols"], source_sha256=sources, input_sha256=parent["input_sha256"],
        kind="mean", independent_stochastic_seeds=0, reused_exact_conditions=sorted(seen),
        existing_control_index=str(Path(existing).resolve()) if existing else None,
        existing_control_index_sha256=C.sha(existing) if existing else None,
        missing_control_list=str(Path(missing).resolve()) if missing else None,
        missing_control_list_sha256=C.sha(missing) if missing else None,
        compact_schema=parent["compact_schema"], forcing="Identical constant ensemble mean; no noise field; native kinetics unchanged.",
        forced_diagnostic="Apply the existing round2 matched-mean rule to postburn compact moments/contrast; E1 lower-sigma flags are reused only where exact companions exist.")
    C.write(out / "manifest.json", manifest)
    C.write(out / "freeze.json", dict(manifest_sha256=C.sha(out / "manifest.json")))
    script = out / "worker.sh"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " + shlex.join([C.PYTHON, "-B", str(out / "frozen_control_sources" / Path(__file__).name), "worker", "--campaign", str(campaign)]) + ' --task "$SLURM_ARRAY_TASK_ID"\n')
    script.chmod(0o750)
    return dict(out=str(out), missing_controls=len(conditions), exact_reused_controls=len(seen))


def worker(campaign, task):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Native control integration requires SLURM")
    import numpy as np
    out = Path(campaign).resolve() / "matched_mean_controls"
    m = C.verify(out)
    c = m["conditions"][task]
    folder = out / "workers" / f'c{c["condition_index"]:03d}'
    folder.mkdir(parents=True, exist_ok=True)
    status = dict(generated_utc=C.now(), status="RUNNING", complete=False, condition=c,
        kind="mean", manifest_sha256=C.sha(out / "manifest.json"), slurm_job_id=os.environ["SLURM_JOB_ID"])
    C.write(folder / "status.json", status)
    start = time.monotonic()
    try:
        t, p, circumference, basal, checks, audit = C.solve(out, c, 0, kind="mean")
        fields, protocol, summary = C.extract(t, p, circumference, basal)
        theta = 2 * np.pi * (np.arange(len(basal)) + .5) / len(basal)
        initial_M1 = c["basal_circumference"] * abs(np.mean(basal * np.exp(1j * theta)))
        fields["R"] = (np.hypot(fields["M1_real"], fields["M1_imag"]) / initial_M1).astype("float32")
        np.savez_compressed(folder / "compact.npz", **fields)
        C.write(folder / "noise_audit.json", audit)
        C.write(folder / "region_audit.json", dict(protocol=protocol, summary=summary))
        status.update(status="COMPLETE", complete=True, compact_sha256=C.sha(folder / "compact.npz"),
            native_horizon=float(t[-1]), validation=checks, mean_R=float(fields["R"].mean()),
            low_contrast_fraction=float(fields["low_contrast"].mean()),
            mean_P_abs=float(fields["P_abs"].mean()), region_count_mean=float(fields["region_count"].mean()))
        union = float(np.mean((fields["R"] < .1) | (fields["contrast_ratio"] < .1)))
        status.update(mean_forced_union_fraction=union, forced_mean=union > .5,
            forced_E1=None, forced_combined=True if union > .5 else None,
            forced_qualification="Exact matched-mean baseline; lower-sigma E1 flag unavailable for newly designed conditions.")
    except Exception as error:
        status.update(status="FAILED", error=f"{type(error).__name__}: {error}", traceback=traceback.format_exc())
        raise
    finally:
        status.update(generated_utc=C.now(), wall_seconds=time.monotonic() - start,
            peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        C.write(folder / "status.json", status)
    return status


def submit(campaign):
    out = Path(campaign).resolve() / "matched_mean_controls"
    m = C.verify(out)
    receipt = out / "results/submission.json"
    if receipt.exists() and C.read(receipt).get("returncode") == 0:
        return C.read(receipt)
    if not m["conditions"]:
        return dict(submitted=False, reason="All exact controls already available")
    command = ["sbatch", "--parsable", "--partition=red", "--qos=xxlarge", "--ntasks=1", "--cpus-per-task=1", "--mem=4G", "--time=02:00:00",
        f'--array=0-{len(m["conditions"])-1}', "--job-name=F4r3-mean", "--output=" + str(out / "logs/%A_%a.out"), "--error=" + str(out / "logs/%A_%a.err"), str(out / "worker.sh")]
    p = subprocess.run(command, text=True, capture_output=True, timeout=45)
    result = dict(generated_utc=C.now(), command=shlex.join(command), returncode=p.returncode,
        stdout=p.stdout, stderr=p.stderr, conditions=len(m["conditions"]), submitted=p.returncode == 0)
    C.write(receipt, result)
    return result


def collect(campaign):
    out = Path(campaign).resolve() / "matched_mean_controls"
    m = C.verify(out)
    rows, missing, failures = [], [], []
    for c in m["conditions"]:
        folder = out / "workers" / f'c{c["condition_index"]:03d}'
        path = folder / "status.json"
        if not path.exists():
            missing.append(c["condition_id"])
            continue
        s = C.read(path)
        if not s.get("complete"):
            failures.append(dict(condition_id=c["condition_id"], status=s["status"], error=s.get("error")))
            continue
        if C.sha(folder / "compact.npz") != s["compact_sha256"]:
            raise ValueError("Matched-mean compact archive changed")
        rows.append(dict(condition_id=c["condition_id"], model=c["model"], size=c["size"], fraction=c["fraction"],
            kind="mean", exact_native_mean_control=True, compact=str(folder / "compact.npz"),
            horizon=s["native_horizon"], mean_R=s["mean_R"], low_contrast_fraction=s["low_contrast_fraction"],
            mean_P_abs=s["mean_P_abs"], region_count_mean=s["region_count_mean"], source="round3_new_exact_matched_mean"))
        rows[-1].update(mean_forced_union_fraction=s["mean_forced_union_fraction"], forced_mean=s["forced_mean"], forced_E1=s["forced_E1"], forced_combined=s["forced_combined"])
    fields = ["condition_id", "model", "size", "fraction", "kind", "exact_native_mean_control", "compact", "horizon", "mean_R", "low_contrast_fraction", "mean_P_abs", "region_count_mean", "source", "mean_forced_union_fraction", "forced_mean", "forced_E1", "forced_combined"]
    with (out / "results/compact_control_index.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    result = dict(generated_utc=C.now(), complete=not missing and not failures, complete_controls=len(rows),
        missing=missing, failures=failures, control_index=str(out / "results/compact_control_index.csv"))
    C.write(out / "status.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "submit", "worker", "collect"))
    parser.add_argument("--campaign", type=Path, default=C.DEFAULT)
    parser.add_argument("--existing", type=Path)
    parser.add_argument("--missing", type=Path)
    parser.add_argument("--task", type=int)
    args = parser.parse_args()
    result = prepare(args.campaign, args.existing, args.missing) if args.action == "prepare" else submit(args.campaign) if args.action == "submit" else worker(args.campaign, args.task) if args.action == "worker" else collect(args.campaign)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
