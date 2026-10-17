#!/usr/bin/env python3
"""Finish verified output retirement with bounded filesystem I/O workers."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from collections import Counter
import json
import os
import shutil
import storage_cleanup_four_assay_20261017 as c


def entries(folder):
    if not folder.exists():
        return
    with os.scandir(folder) as stream:
        items = list(stream)
    for ent in items:
        p = Path(ent.path)
        if ent.is_symlink():
            yield p
        elif ent.is_dir(follow_symlinks=False):
            yield from entries(p)
        else:
            yield p


def units(folder):
    """Separate case directories so workers do not serialize one directory lock."""
    if not folder.exists():
        return
    with os.scandir(folder) as stream:
        items = list(stream)
    directories = [e for e in items if e.is_dir(follow_symlinks=False)]
    contains_keeper = any(folder == c.BASE / rel or folder in (c.BASE / rel).parents
                          for rel in c.KEPT_WORK)
    if not directories and not contains_keeper:
        yield folder, True, len(items)
        return
    for ent in items:
        p = Path(ent.path)
        if ent.is_dir(follow_symlinks=False):
            yield from units(p)
        elif str(p.relative_to(c.BASE)) not in c.KEPT_WORK:
            yield p, False, 1


def main():
    record = json.loads(c.RECORD.read_text())
    assert record['status'] == 'metadata_verified_deleting_retired_payload'
    assert len(record['metadata_archives']) >= len(record['output_roots']) == 19
    record['execution_note'] = ('Serial and file-parallel unlink phases were stopped after all metadata archives verified; '
                                'remaining unlink operations use 16 I/O workers on separate case directories to avoid directory-lock serialization. '
                                'No science integrations or additional deletion scope outside four_assay.')
    extra = [c.BASE / 'pilots/linear_shared_20261001' / name
             for name in ['checkpoints', 'data', 'logs', 'preparations', 'trials']]
    record['additional_output_roots'] = [str(p.relative_to(c.ROOT)) for p in extra]
    for p in extra:
        c.assert_scope(p)
        archive = c.ARCHIVES / ('__'.join(p.relative_to(c.BASE).parts) + '.tar.gz')
        existing = next((r for r in record['metadata_archives'] if r['path'] == str(archive.relative_to(c.ROOT))), None)
        if existing:
            assert c.sha(archive) == existing['sha256']
        else:
            record['metadata_archives'].append(c.archive_metadata(p, '__'.join(p.relative_to(c.BASE).parts)))
        c.write_record(record)
    targets = [c.ROOT / name for name in record['output_roots']] + extra
    counts = Counter()

    def remove(p):
        c.assert_scope(p)
        if str(p.relative_to(c.BASE)) in c.KEPT_WORK:
            return 0
        p.unlink()
        return 1

    def remove_unit(unit):
        p, directory, count = unit
        c.assert_scope(p)
        if directory:
            assert not any(p == c.BASE / rel or p in (c.BASE / rel).parents for rel in c.KEPT_WORK)
            shutil.rmtree(p)
            return count
        return remove(p)

    with ThreadPoolExecutor(max_workers=16) as pool:
        for target in targets:
            counts[str(target.relative_to(c.BASE))] += sum(pool.map(remove_unit, units(target)))
            for folder, dirs, files in os.walk(target, topdown=False, followlinks=False):
                try:
                    Path(folder).rmdir()
                except OSError:
                    pass
            print(json.dumps({'retired_output_root':str(target.relative_to(c.BASE)),
                              'additional_unlinked_files':counts[str(target.relative_to(c.BASE))]}), flush=True)
    viewers = [c.BASE / 'results/reversal_field_report',
               c.BASE / 'results/reversal_field_report_linear_wrong_campaign_archive_20261003']
    record['retired_generated_viewers'] = []
    for folder in viewers:
        raw = list(folder.glob('*.npz'))
        for p in raw:
            remove(p)
        counts[str(folder.relative_to(c.BASE))] += len(raw)
        index = folder / 'index.html'
        old = dict(path=str(index.relative_to(c.ROOT)), old_sha256=c.sha(index), old_bytes=index.stat().st_size)
        temporary = index.with_suffix('.html.retirement.tmp')
        temporary.write_text('<!doctype html><meta charset="utf-8"><title>Retired historical viewer</title>'
                             '<h1>Historical raw viewer retired, 17 October 2026</h1>'
                             '<p>Large deprecated four-assay outputs were removed at the user\'s request. '
                             'Compact records, sources and previews remain.</p>'
                             '<p><a href="../../../manuscript/figures/diagnostic_reversal_fields/index.html">'
                             'Retained selected personal diagnostic</a> (historical, excluded from current manuscript inference).</p>'
                             '<p><a href="../../OUTPUT_RETIREMENT_20261017.md">Retirement details</a></p>')
        temporary.replace(index)
        record['retired_generated_viewers'].append(old)
    leftovers = [str(p.relative_to(c.BASE)) for target in targets for p in entries(target)
                 if str(p.relative_to(c.BASE)) not in c.KEPT_WORK]
    assert not leftovers, leftovers[:20]
    changed = [name for name,digest in record['protected_files_sha256'].items()
               if not (c.ROOT/name).is_file() or c.sha(c.ROOT/name) != digest]
    assert not changed, changed
    record.update(status='complete', finished_at=c.now(), io_workers=16,
                  additional_unlinked_files_by_root=dict(counts),
                  leftover_unselected_output_files=0,
                  protected_files_reverified=len(record['protected_files_sha256']),
                  allocated_bytes_after=c.allocated(c.BASE),
                  scope='Only deprecated four_assay payload retired; manuscript/archive payload unchanged. Retention documentation updated.')
    record['net_allocated_bytes_reclaimed'] = record['allocated_bytes_before'] - record['allocated_bytes_after']
    record['initial_six_root_npz_guaranteed_payload_reclaim_bytes'] = 120004370432
    c.write_record(record)
    print(json.dumps({k:record[k] for k in ['status','allocated_bytes_before','allocated_bytes_after',
                                           'net_allocated_bytes_reclaimed','protected_files_reverified']}), flush=True)


if __name__ == '__main__':
    main()
