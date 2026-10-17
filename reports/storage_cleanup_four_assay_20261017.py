#!/usr/bin/env python3
"""Authorized retirement of deprecated four_assay payload; no model integrations.

All deletions are bounded to explicit output roots below polarity1d/four_assay.
Nonempty worker metadata/logs are losslessly archived and independently checked
before their unpacked copies are removed. Scientific NPZ arrays are retired.
"""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import argparse
import gzip
import hashlib
import io
import json
import os
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'polarity1d/four_assay'
RECORD = ROOT / 'reports/storage_cleanup_four_assay_20261017.json'
ARCHIVES = BASE / 'retained_metadata_20261017'
KEPT_WORK = {
    'work/checkpoints/task_0024.npz',
    'work/checkpoints/task_0025.npz',
    'work/results/trace_0024.npz',
    'work/histories/conditioning_0024.npz',
    'work/histories/release_0025.npz',
}
PRECHECK = {
    'canonical_selected_establishment_npz_sha256_verified': 4952,
    'canonical_selected_diagnostic_npz_sha256_verified': 54,
    'verified_logical_bytes': 564753583,
    'hash_errors': [],
    'data_only_regenerated_corrected_endpoint_rows': 660,
    'new_integrations': 0,
    'main_polarity_png_sha256': 'b9546b972968fc759d00fdea48b1a52bc53f595a9adc80667f97ee40e47684db',
    'supplement_width_response_png_sha256': 'ee4802daaa8f257db99e8aceb2fbed359f7b4ecd68420d9dcdb6fb2120819aeb',
    'both_regenerated_pngs_match_current_canonical_bytes': True,
    'recorded_slurm_job_ids_checked': 37,
    'recorded_slurm_jobs_in_queue': 0,
}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def allocated(path):
    result = subprocess.run(['du', '-s', '-B1', str(path)], text=True,
                            capture_output=True, check=True)
    return int(result.stdout.split()[0])


def write_record(record):
    temporary = RECORD.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(RECORD)


def assert_scope(path):
    # Do not resolve the final path: output trees contain historical symlinks,
    # which must be unlinked rather than followed into their destinations.
    path.absolute().relative_to(BASE.absolute())
    if path == BASE or 'baseline' in path.relative_to(BASE).parts:
        raise RuntimeError(f'Protected or unbounded path: {path}')


def entries(folder):
    """Files and symlinks only; never descend through a symlink."""
    with os.scandir(folder) as stream:
        for ent in stream:
            p = Path(ent.path)
            if ent.is_symlink():
                yield p
            elif ent.is_dir(follow_symlinks=False):
                yield from entries(p)
            elif ent.is_file(follow_symlinks=False):
                yield p
            else:
                raise RuntimeError(f'Unexpected special file: {p}')


def archive_metadata(folder, name):
    destination = ARCHIVES / (name + '.tar.gz')
    if destination.exists():
        raise RuntimeError(f'Refusing to overwrite prior metadata archive: {destination}')
    expected = {}
    omitted_empty = 0
    total_bytes = 0
    links = 0
    with tarfile.open(destination, 'w:gz', compresslevel=6, dereference=False) as archive:
        for p in entries(folder):
            relative = str(p.relative_to(BASE))
            if p.is_symlink():
                info = archive.gettarinfo(str(p), arcname=relative)
                archive.addfile(info)
                expected[relative] = ('symlink', os.readlink(p))
                links += 1
                continue
            if p.suffix == '.npz':
                continue
            if p.stat().st_size == 0:
                omitted_empty += 1
                continue
            # Metadata archives must not silently retain a different huge raw
            # scientific payload in an unexpected format.
            if p.stat().st_size > 128 * 1024 * 1024:
                raise RuntimeError(f'Unexpected large non-NPZ payload: {p}')
            data = p.read_bytes()
            info = archive.gettarinfo(str(p), arcname=relative)
            archive.addfile(info, io.BytesIO(data))
            expected[relative] = ('sha256', hashlib.sha256(data).hexdigest())
            total_bytes += len(data)
    checked = 0
    with tarfile.open(destination, 'r:gz') as archive:
        for member in archive:
            kind, value = expected.pop(member.name)
            if kind == 'symlink':
                if not member.issym() or member.linkname != value:
                    raise RuntimeError(f'Symlink archive mismatch: {member.name}')
            else:
                stream = archive.extractfile(member)
                if stream is None or hashlib.sha256(stream.read()).hexdigest() != value:
                    raise RuntimeError(f'Metadata archive mismatch: {member.name}')
            checked += 1
    if expected:
        raise RuntimeError('Metadata archive has missing entries')
    result = dict(path=str(destination.relative_to(ROOT)), sha256=sha(destination),
                  verified_entries=checked, symlink_entries=links,
                  uncompressed_regular_bytes=total_bytes,
                  empty_files_omitted=omitted_empty,
                  compressed_bytes=destination.stat().st_size)
    print(json.dumps({'metadata_archived': name, **result}), flush=True)
    return result


def protected_hashes():
    result = {}
    for folder in [BASE / 'baseline', BASE / 'reference/chained_v11_20261001',
                   BASE / 'results/assay_explanation_fields']:
        for p in entries(folder):
            if p.is_file() and not p.is_symlink():
                result[str(p.relative_to(ROOT))] = sha(p)
    for rel in KEPT_WORK:
        result[str((BASE / rel).relative_to(ROOT))] = sha(BASE / rel)
    keeper = ROOT / 'polarity1d/archive/four_test_suite_legacy/publication_figures'
    for p in keeper.glob('assay_explanation_v2.*'):
        if p.is_file():
            result[str(p.relative_to(ROOT))] = sha(p)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', required=True)
    parser.parse_args()
    if RECORD.exists():
        raise RuntimeError('Cleanup record already exists; do not execute twice.')
    targets = [BASE / 'work']
    targets += sorted(p for p in (BASE / 'pilots').glob('*/work') if p.is_dir() and not p.is_symlink())
    old = BASE / 'reference/continuation_20261001'
    targets += [old / name for name in ['histories', 'validation', 'inputs', 'profiles', 'logs']]
    for path in targets:
        assert_scope(path)
    record = dict(started_at=now(), status='archiving_metadata_before_deletion',
                  authorization='User explicitly requested removal of large deprecated four-assay outputs.',
                  scope='polarity1d/four_assay only; no changes to manuscript or archives.',
                  preflight=PRECHECK,
                  allocated_bytes_before=allocated(BASE),
                  output_roots=[str(p.relative_to(ROOT)) for p in targets],
                  retained_root_work_npz=sorted(KEPT_WORK),
                  metadata_archives=[])
    preserved = protected_hashes()
    record['protected_files_sha256'] = preserved
    write_record(record)
    ARCHIVES.mkdir(exist_ok=False)
    # All archives must verify before any raw data is deleted.
    for target in targets:
        if not target.exists():
            raise RuntimeError(f'Missing expected output root: {target}')
        name = '__'.join(target.relative_to(BASE).parts)
        record['metadata_archives'].append(archive_metadata(target, name))
        write_record(record)
    record['status'] = 'metadata_verified_deleting_retired_payload'
    write_record(record)
    deleted = Counter()
    inodes = {}

    def remove_file(p, group):
        assert_scope(p)
        if str(p.relative_to(BASE)) in KEPT_WORK:
            return
        st = p.lstat()
        if not p.is_symlink():
            key = (st.st_dev, st.st_ino)
            item = inodes.setdefault(key, [0, st.st_nlink, st.st_blocks * 512])
            item[0] += 1
        p.unlink()
        deleted[group] += 1

    for target in targets:
        group = str(target.relative_to(BASE))
        for p in entries(target):
            remove_file(p, group)
        # Only empty directories are removed; retained selected inputs survive.
        for folder, dirs, files in os.walk(target, topdown=False, followlinks=False):
            p = Path(folder)
            try:
                p.rmdir()
            except OSError:
                pass
        print(json.dumps({'retired_output_root': group, 'unlinked_files': deleted[group]}), flush=True)
    # Duplicated old native-viewer arrays are independent of the canonical
    # manuscript personal diagnostic. Preserve scripts, scalar metadata, previews.
    viewer_dirs = [BASE / 'results/reversal_field_report',
                   BASE / 'results/reversal_field_report_linear_wrong_campaign_archive_20261003']
    record['retired_generated_viewers'] = []
    for folder in viewer_dirs:
        for p in folder.glob('*.npz'):
            remove_file(p, str(folder.relative_to(BASE)))
        index = folder / 'index.html'
        if index.exists():
            old_hash = sha(index)
            old_bytes = index.stat().st_size
            index.write_text('<!doctype html><meta charset="utf-8"><title>Retired historical viewer</title>'
                             '<h1>Historical raw viewer retired, 17 October 2026</h1>'
                             '<p>The user authorized removal of deprecated four-assay large outputs. '
                             'Compact provenance, scripts and previews remain. Native arrays for this '
                             'historical viewer were retired.</p>'
                             '<p><a href="../../../manuscript/figures/diagnostic_reversal_fields/index.html">'
                             'Retained selected personal diagnostic</a> '
                             '(historical experiments, excluded from current manuscript inference).</p>'
                             '<p>See ../../OUTPUT_RETIREMENT_20261017.md for the cleanup record.</p>')
            record['retired_generated_viewers'].append(dict(path=str(index.relative_to(ROOT)),
                                                             old_sha256=old_hash, old_bytes=old_bytes))
    changed = [name for name, digest in preserved.items()
               if not (ROOT / name).is_file() or sha(ROOT / name) != digest]
    if changed:
        raise RuntimeError(f'Protected files changed: {changed}')
    record.update(finished_at=now(), status='complete', unlinked_files_by_output_root=dict(deleted),
                  retired_unique_file_payload_bytes=sum(v[2] for v in inodes.values()),
                  guaranteed_reclaimed_file_payload_bytes=sum(v[2] for v in inodes.values() if v[0] == v[1]),
                  retired_payload_still_linked_elsewhere_bytes=sum(v[2] for v in inodes.values() if v[0] < v[1]),
                  protected_files_reverified=len(preserved),
                  allocated_bytes_after=allocated(BASE))
    record['net_allocated_bytes_reclaimed'] = record['allocated_bytes_before'] - record['allocated_bytes_after']
    write_record(record)
    print(json.dumps({k: record[k] for k in ['status', 'allocated_bytes_before', 'allocated_bytes_after',
                                           'net_allocated_bytes_reclaimed', 'protected_files_reverified']}), flush=True)


if __name__ == '__main__':
    main()
