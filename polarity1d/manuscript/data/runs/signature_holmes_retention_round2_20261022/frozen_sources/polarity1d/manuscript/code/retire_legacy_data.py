#!/usr/bin/env python3
"""Plan or remove duplicate obsolete manuscript data after extracting diagnostics.

The corrected establishment data remain active. Old linear/width datasets are
removed from the core; the personal native-field report retains only its own
selected original reconstruction inputs and frozen source/record provenance.
No original campaign or canonical corrected scientific data are deleted.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import gzip
import hashlib
import json
import shutil

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
MANIFEST = HERE / 'data/provenance/data_sources.json'

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    doc = json.loads(MANIFEST.read_text())
    if not (HERE / 'data/width_reversal').exists():
        print(json.dumps({'already_completed': True, 'original_campaigns_untouched': True}))
        return
    sources = {r['source']: r for r in doc['files']}
    selected = {}
    case_keys = {}
    for p in (HERE / 'data/reversal_fields').glob('*.json'):
        d = json.loads(p.read_text())
        if {'model','length','condition','duration'} <= d.keys():
            case_keys[(d['model'],d['length'],d['condition'])] = d['duration']
        for source in d.get('sources', {}):
            if source.endswith('.npz'):
                r = dict(sources[source])
                r['old_destination'] = r['destination']
                r['destination'] = 'data/reversal_fields/reconstruction_inputs/' + source
                r['role'] = 'personal diagnostic selected original reconstruction input; excluded manuscript evidence'
                selected[r['destination']] = r
    for r in doc['files']:
        dst = r['destination']
        if (dst.startswith('data/width_reversal/frozen_sources/') and not dst.endswith('.gz')) or dst in {
                'data/width_reversal/manifest.json', 'data/width_reversal/provenance.json',
                'data/width_reversal/status.json'}:
            q = dict(r)
            q['old_destination'] = dst
            q['destination'] = dst.replace('data/width_reversal/', 'data/reversal_fields/original_campaign/', 1)
            q['role'] = 'immutable old width campaign provenance used only by personal diagnostic'
            selected[q['destination']] = q
    moves = ['data/linear_v11', 'data/width_reversal',
             'data/inputs/baseline_records.json.gz', 'data/inputs/normalized_v10_records.json.gz',
             'code/consolidate_data.py']
    plan = {'generated_at': datetime.now(timezone.utc).isoformat(),
            'active_manuscript_evidence': 'corrected basal establishment; corrected reversal pending',
            'delete_duplicate_copies': moves, 'diagnostic_input_copies': len(selected),
            'diagnostic_input_bytes': sum(r['bytes'] for r in selected.values()),
            'selected_native_reconstruction_files': sum('reconstruction_inputs/' in k for k in selected),
            'original_campaigns_untouched': True, 'copy_plan': list(selected.values())}
    (HERE / 'data/provenance/legacy_slimming_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    print(json.dumps({k:v for k,v in plan.items() if k != 'copy_plan'}), flush=True)
    if not args.apply:
        return
    for r in selected.values():
        src = HERE / r['old_destination']
        dst = HERE / r['destination']
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        if sha(dst) != r['sha256']:
            raise RuntimeError(f'Checksum mismatch: {dst}')
    original_records = HERE / 'data/width_reversal/records.json.gz'
    with gzip.open(original_records, 'rt') as f:
        old_records = json.load(f)
    subset = []
    for c in old_records['cases']:
        key = (c['model'],c['length'],c['condition'])
        if key in case_keys:
            c['tests'] = [t for t in c['tests'] if t['duration'] == case_keys[key]]
            if len(c['tests']) != 1:
                raise ValueError(f'Selected diagnostic trial ambiguous: {key}')
            subset.append(c)
    assert len(subset) == len(case_keys) == 42
    subset_path = HERE / 'data/reversal_fields/original_campaign/records.json.gz'
    subset_doc = {'generated_at':old_records.get('generated_at'), 'cases':subset,
                  'subset_provenance': {'source':'polarity1d/four_assay/pilots/width_full_lengths_20261002/results/records.json.gz',
                      'source_sha256':sha(original_records), 'original_case_count':len(old_records['cases']),
                      'retained_case_count':42, 'retained_trial_count':42,
                      'purpose':'personal native-field diagnostic only; excluded from current manuscript inference'}}
    with subset_path.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as z:
            z.write(json.dumps(subset_doc,separators=(',',':')).encode())
    subset_entry = {'source':subset_doc['subset_provenance']['source'],
                    'source_sha256':subset_doc['subset_provenance']['source_sha256'],
                    'destination':str(subset_path.relative_to(HERE)), 'sha256':sha(subset_path),
                    'bytes':subset_path.stat().st_size, 'role':'derived 42-case/42-trial personal diagnostic subset'}
    for relative in moves:
        src = HERE / relative
        if src.exists():
            if src.is_dir():
                shutil.rmtree(src)
            else:
                src.unlink()
    retired = ('data/linear_v11/', 'data/width_reversal/')
    kept = [r for r in doc['files'] if not r['destination'].startswith(retired)
            and r['destination'] not in moves]
    kept.extend({k:v for k,v in r.items() if k != 'old_destination'} for r in selected.values())
    kept.append(subset_entry)
    doc.update({'generated_at': datetime.now(timezone.utc).isoformat(),
                'scope': 'Corrected basal establishment plus personal native-field diagnostic only. Corrected reversal results pending.',
                'files': kept, 'total_bytes': sum(r['bytes'] for r in kept),
                'removed_duplicate_legacy_copies': moves})
    doc['path_map'] = {k:v for k,v in doc.get('path_map',{}).items() if 'establishment_basal' in v}
    doc['diagnostic_reconstruction_map'] = {r['source']:r['destination'] for r in selected.values()}
    MANIFEST.write_text(json.dumps(doc, indent=2) + '\n')
    print(json.dumps({'applied': True, 'active_manifest_files':len(kept), 'active_manifest_bytes':doc['total_bytes']}))

if __name__ == '__main__':
    main()
