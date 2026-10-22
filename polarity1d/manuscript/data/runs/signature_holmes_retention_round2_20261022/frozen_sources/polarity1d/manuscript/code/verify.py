"""Verify immutable retained evidence and selected full-state dependencies."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import sys

BASE = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    errors = []
    manifest = json.loads((BASE / 'data/provenance/data_sources.json').read_text())
    checked = 0
    for entry in manifest['files']:
        if not entry['destination'].startswith('data/'):
            continue
        path = BASE / entry['destination']
        if not path.is_file() or sha(path) != entry['sha256']:
            errors.append(f'Immutable copy changed: {entry["destination"]}')
        checked += 1
    sources = json.loads((BASE / 'data/provenance/core_sources.json').read_text())
    for entry in sources['files']:
        path = BASE / entry['retained']
        if not path.is_file() or sha(path) != entry['retained_sha256']:
            errors.append(f'Core source changed since source manifest: {entry["retained"]}')
    establishment = BASE / 'data/establishment_basal'
    with (establishment / 'results/endpoints.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 3234:
        errors.append('Corrected establishment must contain 3,234 cases')
    legi = [row for row in rows if row['model'] == 'legi']
    if len(legi) != 462 or any(row['status'] != 'unpolarized' for row in legi):
        errors.append('Corrected LEGI basal-withdrawal evidence changed')
    checkpoints = 0
    for row in rows:
        if not row['selected_checkpoint']:
            continue
        path = establishment / row['selected_checkpoint']
        if not path.is_file() or sha(path) != row['selected_checkpoint_sha256']:
            errors.append(f'Establishment checkpoint changed: case {row["case_index"]}')
        checkpoints += 1
    diagnostic = len(list((BASE / 'data/reversal_fields').glob('*_full.npz')))
    if diagnostic != 42:
        errors.append('Personal diagnostic must retain 42 native field trajectories')
    result = {'generated_at': datetime.now(timezone.utc).isoformat(),
              'passed': not errors, 'immutable_files_checked': checked,
              'selected_checkpoints_checked': checkpoints,
              'diagnostic_trajectories': diagnostic, 'errors': errors}
    (BASE / 'data/provenance/verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))
    return int(bool(errors))


if __name__ == '__main__':
    sys.exit(main())
