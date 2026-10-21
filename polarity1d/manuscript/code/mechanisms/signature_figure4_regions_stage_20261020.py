"""Prepare read-only A5 observations from complete, audited E1--E4 analyses.

The common arm combines the E2-augmented original sizes with the new E3 sizes.
The longest arm replaces every preselected E4 condition as a whole; it never
pools short and long tracks or silently substitutes a short failed long arm.
This preparer runs no model integration and submits no scheduler jobs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import shlex
import shutil

from signature_figure4_regions_campaign_20261020 import save, sha

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "polarity1d/manuscript").is_dir())
RUNS = ROOT / "polarity1d/manuscript/data/runs"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"


def verified(path, expected):
    path = Path(path).resolve()
    if not expected or sha(path) != expected:
        raise ValueError("A5 stage input identity failed: " + str(path))
    return path


def stage_rows(root):
    root = Path(root).resolve()
    freeze = json.loads((root / "freeze.json").read_text())
    inventory_path = verified(root / "inventory.json", freeze["inventory_sha256"])
    inventory = json.loads(inventory_path.read_text())
    identity_path = root / "results/native_identity_verification.json"
    identity = json.loads(identity_path.read_text())
    status_path = root / "results/analysis/status.json"
    status = json.loads(status_path.read_text())
    if (identity.get("passed") is not True or identity.get("complete") is not True
            or identity.get("inventory_sha256") != sha(inventory_path)
            or status.get("passed") is not True or status.get("complete") is not True
            or status.get("stage") != freeze["stage"]):
        raise ValueError("Complete independently audited native stage analysis is required")
    for name, digest in status["outputs"].items():
        verified(root / "results/analysis" / name, digest)
    arms = json.loads((root / "results/arms.json").read_text())["arms"]
    observed_arms = {float(row["sigma"]): row for row in status["arms"]}
    if set(observed_arms) != {float(row["sigma"]) for row in arms}:
        raise ValueError("Every declared sigma arm needs its own independent audit")
    for arm in arms:
        child = Path(arm["path"])
        receipt = observed_arms[float(arm["sigma"])]
        verified(child / "results/analysis/status.json", receipt["status_sha256"])
        audit_path = verified(child / "results/analysis/independent_estimand_audit.json", receipt["audit_sha256"])
        audit = json.loads(audit_path.read_text())
        if audit.get("passed") is not True or audit.get("condition_groups") != arm["condition_groups"]:
            raise ValueError("An independent stage arm audit is incomplete")
    return freeze["stage"], inventory["cases"], dict(
        root=str(root), stage=freeze["stage"], inventory_sha256=sha(inventory_path),
        identity_sha256=sha(identity_path), analysis_status_sha256=sha(status_path))


def choose_arm(by_stage, arm):
    from signature_figure4_stage_analysis_20261020 import deduplicate
    required = {"E1": {"E1"}, "E2": {"E2"}, "common": {"E2", "E3"},
                "longest": {"E2", "E3", "E4"}}[arm]
    if set(by_stage) != required:
        raise ValueError(f"{arm} requires exactly the complete stages {sorted(required)}")
    if arm in ("E1", "E2"):
        return deduplicate(by_stage[arm])
    common, replacements = deduplicate(by_stage["E2"] + by_stage["E3"])
    if arm == "common":
        if any(row["sigma"] != 1. for row in common if row["kind"] == "noise"):
            raise ValueError("The common arm uses the original sigma-one forcing family")
        return common, replacements
    long_rows, extra = deduplicate(by_stage["E4"])
    selected = {(r["model"], r["C_label"], r["f"])
                for r in long_rows if r["kind"] == "noise"}
    if len(selected) != 52:
        raise ValueError("Every one of the52 preselected E4 conditions is required")
    if any(r["sigma"] != 1. for r in long_rows if r["kind"] == "noise"):
        raise ValueError("E4 cannot replace the original forcing family with another sigma")
    rows = [r for r in common if r["kind"] == "zero"
            or (r["model"], r["C_label"], r["f"]) not in selected] + long_rows
    needed = {(r["model"], r["C_label"], r["native_horizon"])
              for r in rows if r["kind"] == "noise"}
    rows = [r for r in rows if r["kind"] != "zero"
            or (r["model"], r["C_label"], r["native_horizon"]) in needed]
    rows, final = deduplicate(rows)
    horizons = {}
    for row in rows:
        if row["kind"] != "noise":
            continue
        key = (row["model"], row["C_label"], row["f"])
        horizons.setdefault(key, set()).add(row["native_horizon"])
    if any(len(v) != 1 for v in horizons.values()):
        raise ValueError("A longest-arm condition cannot mix short and long realizations")
    return rows, replacements + extra + final


def prepare(out, arm, stages):
    out = Path(out).resolve()
    if not out.is_relative_to(RUNS.resolve()) or out.exists():
        raise ValueError("Use a fresh manuscript/data/runs A5 analysis directory")
    by_stage, authorities = {}, []
    for root in stages:
        stage, rows, authority = stage_rows(root)
        if stage in by_stage:
            raise ValueError("One canonical audited input per native stage is required")
        by_stage[stage] = rows
        authorities.append(authority)
    cases, replacements = choose_arm(by_stage, arm)
    if len({row["id"] for row in cases}) != len(cases):
        raise ValueError("Canonical A5 case identities must be unique")
    for row in cases:
        if not all(row.get(key) is not None for key in ("native_horizon", "profile_cadence", "sigma", "origin")):
            raise ValueError("Every future A5 case must retain its forcing, origin and observation schedule")
    out.mkdir(parents=True)
    for name in ("inputs", "workers", "results", "logs"):
        (out / name).mkdir()
    inventory = dict(created_utc=datetime.now(timezone.utc).isoformat(), cases=cases,
        analysis_arm=arm, source_stages=authorities, counts=dict(Counter(row["kind"] for row in cases)),
        native_integrations=0, canonical_observation_replacements=replacements,
        forcing_family="fixed-mean paired sigma grid" if arm == "E1" else "sigma-one lognormal input",
        horizon_policy="Every selected E4 condition uses only its long arm" if arm == "longest" else "One declared native horizon per condition")
    save(out / "inventory.json", inventory)
    frozen = out / "inputs/frozen_sources"
    frozen.mkdir()
    names = ["signature_figure4_regions_campaign_20261020.py", "signature_figure4_regions_20261020.py",
             "signature_figure4_regions_collect_20261020.py", "signature_figure4_regions_audit_20261020.py",
             "signature_figure4_regions_stage_20261020.py"]
    for name in names:
        shutil.copy2(HERE / name, frozen / name)
    freeze = dict(created_utc=datetime.now(timezone.utc).isoformat(), analysis_arm=arm,
        inventory_sha256=sha(out / "inventory.json"), source_sha256={name: sha(frozen / name) for name in names},
        preparer_sha256=sha(__file__), source_stages=authorities)
    save(out / "inputs/worker_freeze.json", freeze)
    q = shlex.quote
    common = ["#!/bin/bash", "#SBATCH --partition=red", "#SBATCH --qos=small", "#SBATCH --cpus-per-task=1"]
    commands = {
        "measure": (4, "01:00:00", f"#SBATCH --array=0-{len(cases)-1}",
            f"{q(str(frozen / names[0]))} --inventory {q(str(out / 'inventory.json'))} --inventory-sha256 {freeze['inventory_sha256']} --task \"$SLURM_ARRAY_TASK_ID\" --output {q(str(out / 'workers'))}"),
        "audit": (4, "02:00:00", "#SBATCH --array=0-23",
            f"{q(str(frozen / names[3]))} --inventory {q(str(out / 'inventory.json'))} --inventory-sha256 {freeze['inventory_sha256']} --workers {q(str(out / 'workers'))} --output {q(str(out / 'results/independent_all_views'))} --shards 24 --shard \"$SLURM_ARRAY_TASK_ID\""),
        "collect": (8, "04:00:00", None,
            f"{q(str(frozen / names[2]))} --inventory {q(str(out / 'inventory.json'))} --inventory-sha256 {freeze['inventory_sha256']} --workers {q(str(out / 'workers'))} --output {q(str(out / 'results/analysis'))} --worker-freeze {q(str(out / 'inputs/worker_freeze.json'))} --worker-freeze-sha256 {sha(out / 'inputs/worker_freeze.json')}"),
        "audit_collect": (2, "00:30:00", None,
            f"{q(str(frozen / names[4]))} --out {q(str(out))} --collect-audit"),
    }
    for role, (memory, time, array, command) in commands.items():
        lines = common + [f"#SBATCH --job-name=F4-A5-{arm}-{role}", f"#SBATCH --mem={memory}G",
            f"#SBATCH --time={time}", f"#SBATCH --output={out}/logs/{role}_%A_%a.log"]
        if array:
            lines.append(array)
        lines += ["set -euo pipefail", "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1",
                  f"exec {PYTHON} -B {command}"]
        (out / "inputs" / (role + ".slurm")).write_text("\n".join(lines) + "\n")
    return dict(prepared=True, analysis_arm=arm, cases=len(cases), inventory_sha256=freeze["inventory_sha256"],
                output=str(out), submitted=False, native_integrations=0)


def collect_audit(out):
    """Bind every independent shard to the complete canonical measurement set."""
    out = Path(out).resolve()
    freeze = json.loads((out / "inputs/worker_freeze.json").read_text())
    inventory_path = verified(out / "inventory.json", freeze["inventory_sha256"])
    cases = json.loads(inventory_path.read_text())["cases"]
    verified(__file__, freeze["source_sha256"][Path(__file__).name])
    records, shards = [], []
    for shard in range(24):
        path = out / "results/independent_all_views" / f"shard_{shard:03d}.json"
        value = json.loads(path.read_text())
        expected = [case["id"] for i, case in enumerate(cases) if i % 24 == shard]
        actual = [row["case_id"] for row in value["records"]]
        if (value.get("passed") is not True or value.get("inventory_sha256") != freeze["inventory_sha256"]
                or value.get("source_sha256") != freeze["source_sha256"]["signature_figure4_regions_audit_20261020.py"]
                or value.get("shard") != shard or value.get("shards") != 24
                or value.get("cases") != len(expected) or actual != expected):
            raise ValueError("Independent A5 shard identity or canonical case coverage failed")
        for row in value["records"]:
            if row.get("passed") is not True:
                raise ValueError("An independent native-frame audit failed")
            measured = verified(out / "workers" / row["case_id"] / "result.json", row["measurement_sha256"])
            measurement = json.loads(measured.read_text())
            primary = measurement["primary"]
            verified(primary["path"], primary["sha256"])
            expected_views = 24 + int(primary["downsample_factor"] not in (1, 2, 4, 8))
            if row["primary_sha256"] != primary["sha256"] or row["views_checked"] != expected_views:
                raise ValueError("Independent A5 audit does not cover every primary/native view")
        records.extend(value["records"])
        shards.append(dict(path=str(path), sha256=sha(path)))
    if len(records) != len(cases) or len({r["case_id"] for r in records}) != len(cases):
        raise ValueError("Every canonical A5 case requires exactly one independent audit")
    target = out / "results/independent_all_views_audit.json"
    primary_target = out / "results/independent_primary_native_frames_audit.json"
    if target.exists() or primary_target.exists():
        raise FileExistsError("A completed independent A5 audit master is immutable")
    common = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True, complete=True,
        inventory_sha256=freeze["inventory_sha256"], cases=len(cases), native_cases=len(cases),
        source_shards=shards, records=records, source_sha256=sha(__file__), native_integrations=0)
    save(target, dict(common, native_observation_frames_checked=sum(r["all_view_native_frames"] for r in records),
        declared_view_audits=sum(r["views_checked"] for r in records),
        scope="Independent native frame, exact optical/raw baseline, qualification/count/width/budget/event-gate and complete declared-view audit; every measured receipt bound to the canonical inventory."))
    save(primary_target, dict(common, primary_native_frames=sum(r["frames"] for r in records),
        subset_of_all_views=True, all_views_audit_sha256=sha(target),
        scope="Primary-frame projection of the same independent all-view auditor, not a second independent algorithm."))
    return dict(passed=True, complete=True, cases=len(cases), all_views_audit_sha256=sha(target),
                primary_audit_sha256=sha(primary_target))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--arm", choices=("E1", "E2", "common", "longest"))
    parser.add_argument("--stage-analyses", nargs="+", type=Path)
    parser.add_argument("--collect-audit", action="store_true")
    args = parser.parse_args()
    if args.collect_audit:
        if args.arm or args.stage_analyses:
            parser.error("Audit collection uses the already frozen canonical input")
        result = collect_audit(args.out)
    else:
        if not args.arm or not args.stage_analyses:
            parser.error("Preparation requires --arm and --stage-analyses")
        result = prepare(args.out, args.arm, args.stage_analyses)
    print(json.dumps(result, indent=2))
