#!/usr/bin/env python3
"""User-authorized Figure 2 retirement, bounded by current round-3 evidence.

Default: write a reviewable, file-level plan without removing anything.
--execute: verify the sealed plan and unlink only its unchanged regular files.
Retain sources, all inputs/basals, all nonempty metadata, calibration/numerical
checks, all 450 width workers, all 180 plotted train workers, and 14 additional
Holmes noise/control workers supporting the retained supplement. Do not rewrite
native archives or edit the scientific figure builders.
"""
from __future__ import annotations

import argparse
import collections
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "polarity1d/manuscript/data/runs"
DEV = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev"
FIG = DEV / "revision_round3_20261022/figure2"
OUT = ROOT / "reports/figure2_cleanup_round3_20261022"
LEGACY_PREFIXES = (
    "revision_five_model_", "revision_compensation_", "revision_width_",
)


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2) + "\n")


def relative(path):
    return str(path.relative_to(ROOT))


def scopes():
    return sorted(p for p in RUNS.iterdir() if p.is_dir() and not p.is_symlink()
        and (p.name.startswith("signature_") and "20261020" in p.name
             and not p.name.startswith(("signature_figure3_", "signature_figure4_"))
             or p.name.startswith(LEGACY_PREFIXES) and "20261019" in p.name
             or p.name == "signature_holmes_retention_round2_20261022"))


def full_retention(name):
    return name.startswith(("signature_calibration_", "signature_numerics_",
                            "signature_basals_", "signature_width_")) or name in (
        "signature_support_focal_20261020",
        "signature_holmes_retention_round2_20261022",
    )


def selection():
    selected = collections.defaultdict(set)
    for filename in ("C_replacement_times_all450.csv", "D_gap1_vs_pulse2_all180.csv"):
        with (FIG / "supporting/data" / filename).open() as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == (450 if filename.startswith("C") else 180)
        for row in rows:
            selected[Path(row["run"]).name].add(int(row["task"]))
    assert sum(len(v) for k, v in selected.items() if k.startswith("signature_train_")
               or k == "signature_holmes_retention_round2_20261022") == 180
    # Preserve the original main-resolution counterparts of every retained
    # numerical companion, including pulse/width cases outside the main plate.
    pattern = re.compile(r"polarity1d/manuscript/data/runs/([^/\s\"']+)/workers/(\d+)")
    for root in RUNS.glob("signature_numerics_*20261020"):
        for path in (root / "manifest.json", *(root / "results").rglob("*.json")):
            if path.is_file():
                for name, task in pattern.findall(path.read_text()):
                    selected[name].add(int(task))
    return selected


def protected_assets():
    """Pin redraw inputs, current delivery, and protected other main figures."""
    paths = set()
    for s in (ROOT / "manuscript_outline/revision_round3_20261022", FIG,
              DEV / "revision_round2_20261022/figure2",
              DEV / "signature_assays_20261020"):
        paths.update(p for p in s.rglob("*") if p.is_file() and not p.is_symlink())
    for name in ("outline_v5.json", "outline_v6.json", "outline_v6.pdf",
                 "outline_v6.build_manifest.json", "outline_v6.validation.json"):
        paths.add(ROOT / "manuscript_outline" / name)
    for path in read(FIG / "supporting/metadata/source_sha256.json"):
        paths.add(ROOT / path)
    # Immutable publication input manifests may mention further source files.
    for name in ("outline_v6.build_manifest.json",):
        manifest = read(ROOT / "manuscript_outline" / name)
        for row in manifest["inputs"]:
            paths.add(ROOT / row["path"])
    receipt = read(ROOT / "manuscript_outline/revision_round3_20261022/integration.json")
    for section in receipt["protected_sections"].values():
        paths.update(ROOT / p for p in section["asset_sha256"])
    # Every native checkpoint behind panel A is protected and hash-verified.
    with (FIG / "supporting/data/A_incumbent_organization.csv").open() as f:
        for row in csv.DictReader(f):
            path = ROOT / "polarity1d/manuscript" / row["checkpoint"]
            assert sha(path) == row["checkpoint_sha256"], str(path)
            paths.add(path)
    paths.add(Path(__file__).resolve())
    return paths


def plan():
    OUT.mkdir(exist_ok=True)
    selected = selection()
    pins = protected_assets()
    pins_record = {relative(p): sha(p) for p in sorted(pins)}
    roots = scopes()
    rows = []
    kept = []
    by_root = {}
    for root in roots:
        stats = collections.Counter()
        protected = full_retention(root.name)
        for parent, dirs, names in os.walk(root, followlinks=False):
            dirs[:] = [d for d in dirs if not (Path(parent) / d).is_symlink()]
            for name in names:
                path = Path(parent) / name
                s = path.lstat()
                if not stat.S_ISREG(s.st_mode):
                    continue
                parts = path.relative_to(root).parts
                worker = int(parts[1]) if len(parts) > 2 and parts[0] == "workers" and parts[1].isdigit() else None
                is_native = path.suffix in (".npz", ".npy")
                removable = (not protected and is_native and parts[0] in ("workers", "work")
                    and worker not in selected.get(root.name, set())
                    and not any("basal" in part or "preparation" in part for part in parts))
                empty_log = (parts[0] == "logs" and s.st_size == 0
                             and path.suffix in (".out", ".err", ".log"))
                if path in pins:
                    removable = empty_log = False
                record = dict(path=relative(path), allocated_bytes=s.st_blocks * 512,
                    logical_bytes=s.st_size, device=s.st_dev, inode=s.st_ino,
                    nlink=s.st_nlink, mtime_ns=s.st_mtime_ns)
                if removable or empty_log:
                    record["reason"] = "unused Figure 2 native array" if removable else "empty terminal scheduler log"
                    rows.append(record)
                    stats["candidate_files"] += 1
                    stats["candidate_name_allocated_bytes"] += record["allocated_bytes"]
                else:
                    stats["kept_files"] += 1
                    if is_native:
                        record["reason"] = "calibration/numerics/basals/support/width/Holmes supplement" if protected else "plotted worker or retained input/basal"
                        kept.append(record)
                        stats["kept_native_files"] += 1
                        stats["kept_native_name_allocated_bytes"] += record["allocated_bytes"]
        stats["full_native_retention"] = protected
        by_root[relative(root)] = dict(stats)
    baseline = {str(Path(p)): int(n) for n, p in (line.split("\t", 1)
        for line in (OUT / "before.du.tsv").read_text().splitlines())}
    for path, stats in by_root.items():
        stats["before_joint_allocated_bytes"] = baseline[path]
    assert not (set(r["path"] for r in rows) & set(pins_record))
    for root_name, tasks in selected.items():
        for task in tasks:
            folder = RUNS / root_name / "workers" / f"{task:04d}" / "trial"
            if full_retention(root_name) and not folder.exists():
                continue
            assert (folder / "trace.npz").is_file(), str(folder)
            assert list(folder.glob("phase_*.npz")), str(folder)
    # Count a hard-linked inode only if this plan removes every one of its names.
    inode_groups = collections.defaultdict(list)
    for row in rows:
        inode_groups[(row["device"], row["inode"])].append(row)
    reclaimable = sum(v[0]["allocated_bytes"] for v in inode_groups.values()
                      if len(v) == v[0]["nlink"])
    write(OUT / "protected_sha256.json", pins_record)
    write(OUT / "retained_native_inventory.json", kept)
    ledger = OUT / "delete_plan.jsonl"
    with ledger.open("w") as f:
        for row in rows:
            f.write(json.dumps(row) + "\n")
    summary = dict(generated_utc=now(), authorization="User explicitly requested assessment and physical Figure 2 cleanup, retaining the current manuscript_outline/revision_round3_20261022 figure and supporting evidence.",
        method="Unlink only unused native arrays in bounded Figure 2 worker/work trees and empty terminal scheduler logs. Retain all nonempty JSON/scalar results, all frozen sources, inputs, basals, calibration and numerical checks, focal support, current width and temporal workers, and current Holmes supplement workers. No native archive is rewritten.",
        selected_workers={k: sorted(v) for k, v in selected.items()},
        protected_files=len(pins_record), retained_native_files=len(kept),
        candidate_files=len(rows), expected_reclaimable_allocated_bytes=reclaimable,
        expected_reclaimable_GiB=reclaimable / 2**30,
        candidate_inode_names_with_external_links=sum(len(v) for v in inode_groups.values() if len(v) != v[0]["nlink"]),
        roots=by_root, plan_sha256=sha(ledger), protected_inventory_sha256=sha(OUT / "protected_sha256.json"),
        retained_native_inventory_sha256=sha(OUT / "retained_native_inventory.json"),
        scientific_consequence="Unselected historical sweeps lose full native-field re-audit without rerunning; their exact scalar records, original artifact identities, frozen inputs/sources and complete delivered supplements remain. Current plotted native histories remain byte-for-byte unchanged.")
    write(OUT / "assessment.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("selected_workers", "roots")}, indent=2))


def execute():
    summary = read(OUT / "assessment.json")
    assert sha(OUT / "delete_plan.jsonl") == summary["plan_sha256"]
    assert sha(OUT / "protected_sha256.json") == summary["protected_inventory_sha256"]
    assert sha(OUT / "retained_native_inventory.json") == summary["retained_native_inventory_sha256"]
    queue = subprocess.run(["squeue", "-u", "4473331", "-h", "-o", "%i %j %T"],
                           check=True, text=True, capture_output=True).stdout
    assert not any(any(w in line.lower() for w in ("signature", "compensation", "f2-", "round3")) for line in queue.splitlines()), queue
    pins = read(OUT / "protected_sha256.json")
    for path, digest in pins.items():
        assert sha(ROOT / path) == digest, path
    keep = read(OUT / "retained_native_inventory.json")
    for row in keep:
        s = (ROOT / row["path"]).lstat()
        assert stat.S_ISREG(s.st_mode)
        assert (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns) == (row["device"], row["inode"], row["logical_bytes"], row["mtime_ns"]), row["path"]
    rows = [json.loads(line) for line in (OUT / "delete_plan.jsonl").read_text().splitlines()]
    allowed = {ROOT / p for p in summary["roots"]}
    for row in rows:
        p = ROOT / row["path"]
        assert p.is_relative_to(RUNS) and p.parents[len(p.relative_to(RUNS).parts) - 2] in allowed
        assert p.resolve() == p
        assert relative(p) not in pins
        s = p.lstat()
        assert stat.S_ISREG(s.st_mode)
        assert (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_nlink) == (row["device"], row["inode"], row["logical_bytes"], row["mtime_ns"], row["nlink"]), row["path"]
    # All predicates and identities pass before the first unlink.
    started = now()
    with (OUT / "execution.jsonl").open("a", buffering=1) as receipt:
        for index, row in enumerate(rows):
            p = ROOT / row["path"]
            s = p.lstat()
            assert (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns) == (row["device"], row["inode"], row["logical_bytes"], row["mtime_ns"])
            p.unlink()
            receipt.write(json.dumps(dict(path=row["path"], deleted_utc=now())) + "\n")
            if (index + 1) % 2000 == 0:
                print(f"Deleted {index + 1}/{len(rows)} inventoried files", flush=True)
    for path, digest in pins.items():
        assert sha(ROOT / path) == digest, path
    for row in keep:
        s = (ROOT / row["path"]).lstat()
        assert (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns) == (row["device"], row["inode"], row["logical_bytes"], row["mtime_ns"]), row["path"]
    assert all(not (ROOT / r["path"]).exists() for r in rows)
    write(OUT / "execution.json", dict(started_utc=started, finished_utc=now(),
        state="completed", deleted_files=len(rows), all_planned_deletions_confirmed=True,
        protected_hashes_verified=len(pins), retained_native_file_identities_verified=len(keep),
        expected_reclaimable_allocated_bytes=summary["expected_reclaimable_allocated_bytes"],
        queue_before_execution=queue, no_scientific_file_rewrites=True))
    print("Execution complete; retained sources and native files verified.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    execute() if args.execute else plan()
