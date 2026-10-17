"""Publish the final file inventory and current producer/input source snapshot."""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent


def main():
    excluded = {'ARTIFACTS.csv', 'results/source_manifest.json'}
    rows = []
    sources = {}
    for path in sorted(RUN.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts:
            continue
        relative = path.relative_to(RUN).as_posix()
        if relative in excluded:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append(dict(path=relative, bytes=path.stat().st_size, sha256=digest))
        if relative.startswith(('scripts/', 'data/')) or relative == 'README.md':
            sources[relative] = digest
    with (RUN / 'ARTIFACTS.csv').open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=['path', 'bytes', 'sha256'])
        writer.writeheader()
        writer.writerows(rows)
    (RUN / 'results/source_manifest.json').write_text(json.dumps(dict(
        generated=datetime.now(timezone.utc).isoformat(), files=sources,
        interpretation='Final producer/input snapshot. Individual fit records retain their original objective dependency hashes.'
    ), indent=2) + '\n')
    print(f'Inventoried {len(rows)} files; {len(sources)} producer/input hashes.')


if __name__ == '__main__':
    main()
