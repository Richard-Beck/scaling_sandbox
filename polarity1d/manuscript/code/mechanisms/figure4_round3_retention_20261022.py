"""Reviewable, narrowly scoped retirement of old Figure4 native profiles.

Planning never deletes. Execution is gated on the complete new compact cohort,
final whole-seed precision and a successful final figure rebuild. Preserved
scalar/descriptor evidence is explicitly measured rather than called negligible.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import json
import os
from pathlib import Path
import stat

import figure4_round3_campaign_20261022 as C

PATTERNS = ("vision_noise_*20261017*", "revision_noise_*20261019*",
            "revision_spring_noise_*20261019*", "signature_figure4_*20261020*")
NATIVE_NAMES = {"primary.npz", "dense_primary.npz", "full_states.npz", "final_state.npz"}
REPORT = C.ROOT / "reports/figure4_round3_retirement_20261022"


def roots():
    base = C.MS / "data/runs"
    return sorted({p for pattern in PATTERNS for p in base.glob(pattern) if p.is_dir()})


def selected(path, root):
    relative = path.relative_to(root)
    if any(part in ("inputs", "frozen_sources") for part in relative.parts):
        return False
    return path.name in NATIVE_NAMES or (
        "native_field_CSV_archives" in relative.parts and path.name.endswith(".csv.gz"))


def inventory(campaign_roots):
    records, inode_paths = [], defaultdict(list)
    for root in campaign_roots:
        for base, dirs, files in os.walk(root, followlinks=False):
            dirs[:] = [d for d in dirs if not (Path(base) / d).is_symlink()]
            for name in files:
                path = Path(base) / name
                s = path.lstat()
                if not stat.S_ISREG(s.st_mode):
                    continue
                row = dict(path=str(path.relative_to(C.ROOT)), campaign=root.name,
                    allocated_bytes=s.st_blocks * 512, logical_bytes=s.st_size,
                    device=s.st_dev, inode=s.st_ino, nlink=s.st_nlink, mtime_ns=s.st_mtime_ns,
                    selected=selected(path, root))
                records.append(row)
                inode_paths[s.st_dev, s.st_ino].append(row)
    total, remove, retained, reclaim = 0, 0, 0, 0
    for paths in inode_paths.values():
        allocated = paths[0]["allocated_bytes"]
        total += allocated
        targets = [r for r in paths if r["selected"]]
        if targets:
            remove += allocated
        if len(targets) != len(paths):
            retained += allocated
        if len(targets) == paths[0]["nlink"]:
            reclaim += allocated
    return records, dict(total_unique_allocated_bytes=total,
        selected_unique_allocated_bytes=remove, retained_unique_allocated_bytes=retained,
        estimated_reclaimable_bytes=reclaim, files=len(records),
        selected_paths=sum(r["selected"] for r in records),
        hard_link_rule="Each inode counted once across72campaign roots; reclaimed allocation requires every hard link selected. Links outside this inventory reduce the reclaim estimate.")


def document(summary, state="PLANNED_NOT_EXECUTED"):
    begin = "<!-- figure4-round3-old-native-retirement-20261022 begin -->"
    end = "<!-- figure4-round3-old-native-retirement-20261022 end -->"
    execution = summary.get("execution")
    execution_note = (f"Executed after the COMPLETE gate: before {execution['before']['total_unique_allocated_bytes']:,} allocated bytes; after {execution['after']['total_unique_allocated_bytes']:,}; measured reclaim {execution['measured_reclaimed_allocated_bytes']:,} bytes. The execution ledger is authoritative for actual removals."
        if execution else "No native files have been removed by this plan.")
    text = f"""{begin}

## Figure 4 old native-profile retirement plan, 22 October 2026

Measured **{summary['generated_utc']}**. Status **{state}**.
The exact [reviewable plan](reports/figure4_round3_retirement_20261022/summary.json)
and [path/stat ledger](reports/figure4_round3_retirement_20261022/delete_plan.jsonl)
cover {len(summary['roots'])} old Figure 4 noise campaign roots:
**{summary['total_unique_allocated_bytes']:,} allocated bytes**
({summary['total_unique_allocated_bytes']/2**30:.3f} GiB), counting hard links once.
The {summary['selected_paths']:,} selected native profile/coupled-state paths
occupy **{summary['selected_unique_allocated_bytes']:,} unique allocated bytes**;
estimated reclaim is **{summary['estimated_reclaimable_bytes']:,} bytes**
({summary['estimated_reclaimable_bytes']/2**30:.3f} GiB), allowing for preserved hard links.

Selection is limited to `primary.npz`, `dense_primary.npz`, `full_states.npz`,
`final_state.npz`, and full native-field CSV shards beneath
`results/native_field_CSV_archives/` in those named roots. Exact qualified basal
inputs, frozen sources/protocols and all repository round-2 CSVs/figures remain.
The **{summary['retained_unique_allocated_bytes']:,} allocated bytes** of remaining
scalar and descriptor evidence are measured scientific records: legacy
`trace.npz`/scalar histories support historical M1/M1initial and q1=P diagnostics,
matched-mean forced flags and 12-seed comparisons; compressed contour/track JSON
preserves previously measured region identities and sensitivity audits. These
are not retained full native profile histories and are not assumed negligible.

Execution remains gated on the new compact campaign's COMPLETE audit, final
precision receipt and successful final rebuild. No deletion is implied by this
inventory. A separate execution receipt records actual before/after storage
when the authorized operation is performed.
{execution_note}

{end}"""
    for name in ("DATA_MANIFEST.md", "RETENTION.md"):
        path = C.ROOT / name
        content = path.read_text()
        if begin in content:
            a, b = content.index(begin), content.index(end) + len(end)
            content = content[:a] + text + content[b:]
        else:
            content = content.rstrip() + "\n\n" + text + "\n"
        path.write_text(content)


def plan():
    REPORT.mkdir(parents=True, exist_ok=True)
    campaign_roots = roots()
    records, summary = inventory(campaign_roots)
    summary.update(generated_utc=C.now(), status="PLANNED_NOT_EXECUTED",
        roots=[str(p.relative_to(C.ROOT)) for p in campaign_roots],
        selection_names=sorted(NATIVE_NAMES), native_field_archive_subfolder="results/native_field_CSV_archives",
        preserved=["inputs/ and frozen_sources/ everywhere", "trace.npz and scalar CSV histories", "region/contour descriptor JSON", "repository round2 compact CSVs and figures"],
        authorization="User: after COMPLETE new compact campaign, delete OLD Figure4 noise native HPC outputs; no deletion before that gate.")
    with (REPORT / "delete_plan.jsonl").open("w") as stream:
        for row in records:
            if row["selected"]:
                stream.write(json.dumps(row, separators=(",", ":")) + "\n")
    groups = Counter()
    seen = set()
    for row in records:
        key = row["device"], row["inode"]
        if key not in seen and not row["selected"]:
            seen.add(key)
            name = Path(row["path"]).name
            category = "trace.npz scalar histories" if name == "trace.npz" else "compressed contour/tracking descriptors" if name.endswith(".json.gz") else "other preserved inputs/sources/scalars/analysis"
            groups[category] += row["allocated_bytes"]
    summary["retained_measured_categories"] = dict(groups)
    C.write(REPORT / "summary.json", summary)
    document(summary)
    return summary


def execute(campaign):
    campaign = Path(campaign).resolve()
    status = C.read(campaign / "status.json")
    if not status.get("complete") or status.get("status") != "COMPLETE":
        raise RuntimeError("New compact campaign is not COMPLETE; old native deletion is prohibited")
    rebuild = C.read(campaign / "results/final_rebuild_execution.json")
    if rebuild.get("returncode") != 0:
        raise RuntimeError("Final compact figure rebuild has not succeeded")
    precision = C.read(campaign / "results/allocation.json").get("precision_source")
    if not precision or not C.collect(campaign, precision)["complete"]:
        raise RuntimeError("Fresh compact integrity/precision audit did not confirm COMPLETE")
    summary = C.read(REPORT / "summary.json")
    records = [json.loads(line) for line in (REPORT / "delete_plan.jsonl").read_text().splitlines()]
    before = inventory([C.ROOT / p for p in summary["roots"]])[1]
    for row in records:
        path = C.ROOT / row["path"]
        s = path.lstat()
        # st_dev is a host-local mount identifier and can differ on another
        # SLURM node. NFS inode, exact path, size and nanosecond mtime are stable.
        if (s.st_ino, s.st_size, s.st_mtime_ns) != (row["inode"], row["logical_bytes"], row["mtime_ns"]):
            raise RuntimeError("A planned native file changed: " + row["path"])
        if not selected(path, C.MS / "data/runs" / row["campaign"]):
            raise RuntimeError("A ledger path is outside the allowed native selector")
    removed = []
    with (REPORT / "execution_paths.jsonl").open("w") as stream:
        for row in records:
            (C.ROOT / row["path"]).unlink()
            removed.append(row["path"])
            stream.write(json.dumps(dict(path=row["path"], deleted_utc=C.now())) + "\n")
            stream.flush()
    after = inventory([C.ROOT / p for p in summary["roots"]])[1]
    result = dict(generated_utc=C.now(), complete=True, removed_paths=len(removed),
        before=before, after=after, measured_reclaimed_allocated_bytes=before["total_unique_allocated_bytes"] - after["total_unique_allocated_bytes"],
        compact_campaign_manifest_sha256=C.sha(campaign / "manifest.json"),
        compact_complete_seeds=status["native_seeds_complete"], original_plan_sha256=C.sha(REPORT / "delete_plan.jsonl"))
    C.write(REPORT / "execution.json", result)
    summary.update(generated_utc=result["generated_utc"], execution=result)
    C.write(REPORT / "summary.json", summary)
    document(summary, "EXECUTED_AFTER_COMPACT_COMPLETE")
    C.update_storage(campaign)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "execute"))
    parser.add_argument("--campaign", type=Path, default=C.DEFAULT)
    args = parser.parse_args()
    result = plan() if args.action == "plan" else execute(args.campaign)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
