#!/usr/bin/env python3
"""Measure only the October20 Figure2--4 roots and update their manifest section."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "polarity1d/manuscript/data/runs"
DEV = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev"
START = "### Figures 2–4 signature, distance and noise revision additions, 2026-10-20"
END = "### Literature-model implementations and reproduction data"


def git_manifest_row():
    snapshot = ROOT / "docs/tasks/briefs/20261020_revision_git_storage.json"
    if not snapshot.is_file():
        return None
    record = json.loads(snapshot.read_text())
    size = record["allocated_bytes"]
    measured = record["measurement_finished_utc"]
    return (f"| `.git/` after approved revision-payload staging; measured {measured} | "
            f"**{size / 2**30:.3f} GiB** ({size:,} allocated bytes) | "
            "Versioned scientific source and portable figure/evidence exports plus preceding repository history. "
            "These duplicate the exported evidence rather than adding independent experiments. Native archives stay outside Git. "
            "The [measurement record](docs/tasks/briefs/20261020_revision_git_storage.json) identifies this staging snapshot; "
            "small documentation/index/commit objects can follow. No deletion authorized. |")


def main(final=False):
    started = datetime.now(timezone.utc).isoformat()
    roots = []
    for path in sorted(RUNS.glob("signature_*20261020*")):
        if not path.is_dir():
            continue
        group = ("Figure3" if path.name.startswith("signature_figure3_") else
                 "Figure4" if path.name.startswith("signature_figure4_") else "Figure2")
        roots.append((path, group, "native_run_or_analysis"))
    for name, group in (("signature_assays_20261020", "Figure2"),
                        ("revision_program_20261020/figure3_distance", "Figure3"),
                        ("signature_figure4_20261020", "Figure4")):
        path = DEV / name
        if path.is_dir():
            roots.append((path, group, "figure_exports"))
    integrated = ROOT / "manuscript_outline/revision_20261020"
    if integrated.is_dir():
        roots.append((integrated, "Integration", "report_exports"))
    for suffix in (".json", ".pdf", ".build_manifest.json", ".validation.json"):
        path = ROOT / "manuscript_outline" / ("outline_v4" + suffix)
        if path.is_file():
            roots.append((path, "Integration", "report_exports"))
    # One traversal counts each inode once across these disjoint roots.
    result = subprocess.run(["du", "-x", "-s", "-B1", "--", *(str(p) for p, _, _ in roots)],
                            check=True, text=True, capture_output=True)
    measured = {}
    for line in result.stdout.splitlines():
        size, path = line.split("\t", 1)
        measured[Path(path)] = int(size)
    rows, totals, native, exports = [], defaultdict(int), defaultdict(int), defaultdict(int)
    for path, group, category in roots:
        size = measured[path]
        rows.append(dict(path=str(path.relative_to(ROOT)), group=group, category=category,
                         allocated_bytes=size, allocated_GiB=size / 2**30))
        totals[group] += size
        (native if category == "native_run_or_analysis" else exports)[group] += size
    finished = datetime.now(timezone.utc).isoformat()
    record = dict(measurement_started_utc=started, measurement_finished_utc=finished,
        method="One GNU du -x -s -B1 traversal of disjoint native/analysis and figure/report roots; each inode counted once within the joint traversal. Active producers can grow during measurement. No deletion.",
        final_snapshot=final, paths=rows, group_allocated_bytes=dict(totals),
        group_allocated_GiB={k: v / 2**30 for k, v in totals.items()},
        native_allocated_bytes=dict(native), export_allocated_bytes=dict(exports))
    output = ROOT / "docs/tasks/briefs/20261020_revision_storage_current.json"
    output.write_text(json.dumps(record, indent=2) + "\n")
    descriptions = {
        "Figure2": "Five-model signature assays: exact qualified basals, frozen sources and full coupled fields; 36 calibration-worker audits covering 528 native evaluations; 450 equal-commanded-area width runs; 1,458 defined cue trains; 486 bounded-noise trajectories; supporting assays and numerical companions. The final main and 82-page supplement retain 162 calibration-unavailable Holmes cells, Wave-pinning L19 noise losses and deBelly timing/field sensitivities. Superseded pilots and failures remain preserved; native arrays stay local.",
        "Figure3": "Four-model fixed-size distance assays: 12 qualified basals, 176 A1/mirror histories, corrected full-region handover roots, 545 A2 pulse trials, 381 B trials including mirrors, and mesh/time-step controls. Final geometry/direct-heading displays preserve all bounds and numerical exceptions. The 1,150-history atlas has pixel-identical lossless WebP exports; original PNGs and native arrays remain local. Rejected fixed-footprint and discontinuous-boundary analyses remain diagnostic provenance.",
        "Figure4": "Four-model noisy-ring revision: corrected E0 inventory and 110 condition analyses; delivered 34-page interim; strict optical-region measurements and independent native-frame checks. Expanded stages contain 864 E1, 468 E2 and 528 E3 native runs. The literal E4 rule selects 52 noisy conditions, giving 624 noisy realizations plus 64 matched mean/zero controls. Original attempts, failed strict prefix comparisons, exact-cap restart pilots and actual-failure-only replacements are retained; selected replacements do not add independent seeds. S11 full native-field CSV.gz archives preserve exact plotted values on HPC; portable atlas figures and compact inventories retain their paths and hashes. The completed final main, 1,291-page statistical core/index and all 13 native-atlas volumes (578 pages) passed independent source/display review. H1 and H3 remain inconclusive and H2 fails its declared overall criterion. Whole-seed uncertainty, finite-window/coverage qualifications and the original Spring sampling failure remain explicit. Superseded export attempts are retained locally and the exact portable publication allowlist is preserved. Superseded compact-basal, pairwise-only and antipodal-route analyses remain diagnostic provenance. Original 936+26 native sources belong to the separately inventoried 19 October group.",
        "Integration": "Static integrated report exports and exact source/assembly identities; no native integrations are launched by report assembly.",
    }
    patterns = {
        "Figure2": "`polarity1d/manuscript/data/runs/signature_*20261020*/` excluding Figure3/4 roots, and `dev/signature_assays_20261020/` exports",
        "Figure3": "`polarity1d/manuscript/data/runs/signature_figure3_*20261020*/` and `dev/revision_program_20261020/figure3_distance/` exports",
        "Figure4": "`polarity1d/manuscript/data/runs/signature_figure4_*20261020*/` and `dev/signature_figure4_20261020/` exports",
        "Integration": "`manuscript_outline/revision_20261020/` and the `outline_v4` JSON, PDF, build manifest and validation record",
    }
    state = ("All required producers have finished; small documentation/receipt changes can follow this measurement."
             if final else "Figure4 producers and report preparation remain active; this is a growing snapshot, not a final footprint.")
    lines = [START, "", f"Allocated storage was measured from **{started}** to **{finished}** over **{len(roots)} disjoint roots**.",
             state, "The [per-path record](docs/tasks/briefs/20261020_revision_storage_current.json) records allocated bytes, native versus export components and exact paths. Each inode is counted once in this joint traversal. Original failed or superseded work is retained. Inventory does not authorize deletion. Historical parent totals retain their dates.", "",
             "| Path or related group | Allocated GiB | Purpose and importance/status |", "| --- | ---: | --- |"]
    for group, size in totals.items():
        lines.append(f"| {patterns[group]} | **{size / 2**30:.3f} combined** ({size:,} allocated bytes; native/analysis {native[group]:,}, exports {exports[group]:,}) | {descriptions[group]} |")
    git_row = git_manifest_row()
    if git_row:
        lines.append(git_row)
    lines.extend(["", "These are component measurements of the revision and repository metadata. Do not add them to historical manuscript/workspace parent totals as if those totals were current.", "", ""])
    manifest = ROOT / "DATA_MANIFEST.md"
    text = manifest.read_text()
    if text.count(START) != 1 or text.count(END) != 1:
        raise ValueError("Expected one precisely scoped October20 manifest section")
    first, last = text.index(START), text.index(END)
    manifest.write_text(text[:first] + "\n".join(lines) + text[last:])
    print(json.dumps(dict(paths=len(rows), native_allocated_bytes=dict(native),
                         export_allocated_bytes=dict(exports), group_allocated_GiB=record["group_allocated_GiB"]), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final", action="store_true", help="Use only after all required producers are terminal")
    main(parser.parse_args().final)
