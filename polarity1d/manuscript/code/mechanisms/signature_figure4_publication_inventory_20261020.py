"""Build an exact Figure 4 Git allowlist, preserving every local work product.

The allowlist comes from final package receipts and accepted immutable E0
receipts, never a recursive proposal to publish the entire topic directory.
Unlisted files are measured and inventoried for local retention; none is removed.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from html.parser import HTMLParser
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import unquote, urlsplit

import signature_figure4_package_20261020 as Q


def measured(path, digest=None):
    path = Path(path).resolve()
    if digest is not None:
        Q.verified_file(dict(path=str(path), sha256=digest), path.parent)
    return dict(path=str(path), sha256=Q.sha(path), bytes=path.stat().st_size,
                allocated_bytes=path.stat().st_blocks * 512)


def accepted_records(receipt_path):
    """Verify archived and original bytes; only mutable navigation can differ."""
    receipt_path = Path(receipt_path).resolve(); payload = json.loads(receipt_path.read_text())
    artifacts = payload["artifacts"]
    if isinstance(artifacts, dict):
        rows = [dict(original=str(Path(payload["source_topic"]) / name),
                     path=str(receipt_path.parent / name), sha256=digest) for name, digest in artifacts.items()]
    else:
        rows = artifacts
    portable, archive_only = [], []
    for row in rows:
        archived = measured(row["path"], row["sha256"])
        original = Path(row["original"]).resolve()
        if original.name == "README.md" and (not original.exists() or Q.sha(original) != row["sha256"]):
            archive_only.append(dict(archived, reason="Immutable delivered navigation snapshot; current organizer guide is separate"))
        else:
            portable.append(dict(measured(original, row["sha256"]), role="accepted_E0_evidence", archived_path=str(Path(row["path"]).resolve())))
    return portable, archive_only


def local_reason(path):
    name = Path(path).name
    if name.startswith("figure4_S7_E1_v1_"):
        return "Superseded S7 layout v1; exact scientific source unchanged; retained on HPC"
    if name.startswith(("figure4_S7_E1_v2_", "figure4_S8_E2_v1_")):
        return "Completed intermediate stage display/source receipt; final package contains accepted pages and exact plotted CSVs"
    if "layout_fixture" in name:
        return "Synthetic layout fixture; not scientific evidence"
    return "Local work/provenance product outside the receipt-defined publication scope; preserved on HPC"


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore"); writer.writeheader(); writer.writerows(rows)


def verify_public_links(records, topic, repo):
    """A local work file must not become a dead public link via a live guide."""
    class Links(HTMLParser):
        def __init__(self):
            super().__init__(); self.values = []
        def handle_starttag(self, tag, attrs):
            self.values.extend(value for key, value in attrs if key in ("href", "src") and value)
    result = subprocess.run(["git", "-C", str(repo), "ls-files", "-z"], capture_output=True)
    tracked = {str((repo / name.decode()).resolve()) for name in result.stdout.split(b"\0") if name}
    allowed = {str(Path(row["path"]).resolve()) for row in records.values()}
    count = 0
    for row in records.values():
        path = Path(row["path"])
        if path.suffix.lower() not in (".html", ".md"):
            continue
        text = path.read_text(); links = Links(); links.feed(text)
        values = links.values + re.findall(r"\[[^\]]*\]\(([^\n)]+)\)", text)
        for value in values:
            parsed = urlsplit(value.strip().strip("<>"))
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            target = (path.parent / unquote(parsed.path)).resolve()
            in_topic = target == topic or topic in target.parents
            if target.is_dir():
                # Organizer navigation may point at a directory, but a local
                # directory containing only superseded work is not public.
                candidates=allowed|tracked
                public=any(target in Path(candidate).parents and Path(candidate).is_file() for candidate in candidates)
                if not public:raise ValueError('A public guide links an empty or local-only directory: '+str(target))
            elif (not target.is_file() or str(target) not in allowed and (in_topic or str(target) not in tracked)):
                raise ValueError("A public reader/guide links an unapproved local-only or missing artifact: " + str(target))
            count += 1
    return count


def inventory(final_audit_path, accepted_receipts, topic, repo, output_prefix, extras=()):
    topic, repo = Path(topic).resolve(), Path(repo).resolve()
    final_audit_path = Path(final_audit_path).resolve()
    final = json.loads(final_audit_path.read_text())
    if final.get("passed") is not True or final.get("final_export_eligible") is not True or final.get("status") != "final":
        raise ValueError("Publication allowlist requires the actual final eligible export, not a candidate")
    for key in ("candidate_export", "completion_review", "content"):
        Q.verified_file(final[key], final_audit_path.parent)
    content = json.loads(Path(final["content"]["path"]).read_text())
    if content.get("status") != "final":
        raise ValueError("Publication content is not final")
    records = {}

    def add(record, role):
        actual = measured(record["path"], record.get("sha256"))
        if actual["bytes"] >= 100 * 1024 ** 2:
            raise ValueError("No portable Git artifact may be 100MiB or larger")
        actual["repo_path"] = Q.repo_path(actual["path"], repo)
        if "\n" in actual["repo_path"] or "\0" in actual["repo_path"]:
            raise ValueError("Invalid repository path in publication receipt")
        previous = records.get(actual["repo_path"])
        if previous and previous["sha256"] != actual["sha256"]:
            raise ValueError("Two publication receipts disagree on exact artifact bytes")
        actual["roles"] = ";".join(sorted((set((previous or {}).get("roles", "").split(";")) | {role}) - {""}))
        records[actual["repo_path"]] = actual

    for record in final["artifacts"]:
        add(record, "accepted_final_package")
    add(dict(path=str(final_audit_path), sha256=Q.sha(final_audit_path)), "final_export_receipt")
    for key in ("candidate_export", "completion_review", "content"):
        add(final[key], key)
    provenance = topic / "supporting/metadata/accepted_provenance"; provenance.mkdir(parents=True, exist_ok=True)
    archive_only, receipt_records = [], []
    for receipt_path in accepted_receipts:
        receipt_path = Path(receipt_path).resolve()
        portable, local = accepted_records(receipt_path)
        archive_only.extend(local)
        for record in portable:
            add(record, "accepted_E0_evidence")
        target = provenance / (receipt_path.parent.name + "__" + receipt_path.name)
        if target.exists() and Q.sha(target) != Q.sha(receipt_path):
            raise ValueError("Accepted provenance snapshot cannot be replaced")
        if not target.exists():
            shutil.copyfile(receipt_path, target)
        add(dict(path=str(target), sha256=Q.sha(target)), "accepted_E0_receipt")
        receipt_records.append(dict(source=str(receipt_path), source_sha256=Q.sha(receipt_path),
                                    portable_path=Q.repo_path(target, repo)))
    for path in extras:
        add(dict(path=str(Path(path).resolve())), "explicit_navigation_or_reproducibility_metadata")
    public_links_checked = verify_public_links(records, topic, repo)
    output_prefix = Path(output_prefix).resolve(); output_prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {kind: Path(str(output_prefix) + suffix) for kind, suffix in (
        ("manifest", ".json"), ("portable", "_portable.csv"), ("local", "_local_only.csv"), ("git_paths", "_git_paths.nul"))}
    self_paths = {Q.repo_path(path, repo) for path in paths.values()}
    if self_paths & set(records):raise ValueError('Publication ledgers cannot be supplied as already hashed payload artifacts')
    local = []
    for path in sorted(topic.rglob("*")):
        if not path.is_file():
            continue
        relative = Q.repo_path(path, repo)
        if relative not in records and relative not in self_paths:
            local.append(dict(measured(path), repo_path=relative, reason=local_reason(path)))
    portable = sorted(records.values(), key=lambda row: row["repo_path"])
    native_volumes={str((repo/row['path']).resolve()) for row in content.get('supplement_volumes',[])}
    core={str((repo/content[key]).resolve()) for key in ('main_pdf','supplement_pdf') if content.get(key)}
    core.update(str((repo/row['path']).resolve()) for row in content.get('assets',[]))
    groups=('accepted_historical','final_main_and_statistical_core','final_S11_PDF_volumes',
            'compact_companions_navigation_and_provenance')
    for row in portable:
        roles=set(row['roles'].split(';'))
        if row['path'] in native_volumes:group='final_S11_PDF_volumes'
        elif row['path'] in core:group='final_main_and_statistical_core'
        elif roles & {'accepted_E0_evidence','accepted_E0_receipt'} and 'accepted_final_package' not in roles:
            group='accepted_historical'
        else:group='compact_companions_navigation_and_provenance'
        row['storage_group']=group
    storage_groups={group:dict(files=sum(row['storage_group']==group for row in portable),
        logical_bytes=sum(row['bytes'] for row in portable if row['storage_group']==group),
        allocated_bytes=sum(row['allocated_bytes'] for row in portable if row['storage_group']==group)) for group in groups}
    write_csv(paths["portable"], portable, ("repo_path", "sha256", "bytes", "allocated_bytes", "storage_group", "roles", "path"))
    write_csv(paths["local"], local, ("repo_path", "sha256", "bytes", "allocated_bytes", "reason", "path"))
    git_paths = sorted(set(records) | self_paths)
    paths["git_paths"].write_bytes(b"".join(path.encode() + b"\0" for path in git_paths))
    report = dict(passed=True, final_only=True, generated_utc=datetime.now(timezone.utc).isoformat(),
        final_export_audit=dict(path=str(final_audit_path), sha256=Q.sha(final_audit_path)),
        receipt_defined_publication=True, recursive_topic_publication=False, deleted_files=0,
        public_links_checked=public_links_checked, no_local_only_or_missing_public_links=True,
        portable_records=portable, portable_files=len(portable),
        portable_bytes=sum(row["bytes"] for row in portable),
        portable_allocated_bytes=sum(row["allocated_bytes"] for row in portable),
        portable_storage_groups=storage_groups,
        storage_group_scope='Disjoint payload groups; shared-role files counted once. Publication ledger/pathspec metadata are measured separately below.',
        every_portable_file_below_100MiB=True,
        accepted_receipts=receipt_records, immutable_archived_navigation_retained_local=archive_only,
        local_only_files=len(local), local_only_bytes=sum(row["bytes"] for row in local),
        local_only_allocated_bytes=sum(row["allocated_bytes"] for row in local),
        local_only_inventory=measured(paths["local"]),
        native_field_archives="Remain on HPC; exact numeric inventories are portable final package artifacts, never absent public CSV links",
        self_publication_paths=sorted(self_paths),
        portable_inventory=measured(paths["portable"]), git_pathspec=measured(paths["git_paths"]),
        source_sha256=Q.sha(__file__))
    # The Git pathspec includes these four audit files too. Measure their
    # actual bytes, including this JSON's own size, without embedding a
    # self-referential digest. A small fixed-point write stabilizes integer
    # byte-count field widths; every count then equals the written file.
    stable=False
    for _ in range(12):
        before={kind:(path.stat().st_size,path.stat().st_blocks*512) if path.exists() else (0,0)
                for kind,path in paths.items()}
        report['publication_inventory_metadata']=[dict(kind=kind,repo_path=Q.repo_path(paths[kind],repo),
            logical_bytes=size,allocated_bytes=allocated) for kind,(size,allocated) in before.items()]
        report['publication_inventory_metadata_logical_bytes']=sum(value[0] for value in before.values())
        report['publication_inventory_metadata_allocated_bytes']=sum(value[1] for value in before.values())
        report['complete_Git_allowlist_files']=len(portable)+len(paths)
        report['complete_Git_allowlist_logical_bytes']=report['portable_bytes']+report['publication_inventory_metadata_logical_bytes']
        report['complete_Git_allowlist_allocated_bytes']=report['portable_allocated_bytes']+report['publication_inventory_metadata_allocated_bytes']
        report['self_manifest_hash']='Not embedded recursively; obtain the final manifest SHA256 from the external handoff receipt.'
        paths["manifest"].write_text(json.dumps(report,indent=2)+'\n')
        after={kind:(path.stat().st_size,path.stat().st_blocks*512) for kind,path in paths.items()}
        if before==after:stable=True;break
    if not stable:raise ValueError('Publication metadata storage counts did not stabilize')
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-audit", required=True, type=Path)
    parser.add_argument("--accepted-receipt", action="append", required=True, type=Path)
    parser.add_argument("--topic", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--extra-portable", action="append", default=[], type=Path)
    args = parser.parse_args()
    result = inventory(args.final_audit, args.accepted_receipt, args.topic, args.repo, args.output_prefix, args.extra_portable)
    print(json.dumps({key: result[key] for key in ("passed", "portable_files", "portable_bytes", "local_only_files", "deleted_files")}, indent=2))
