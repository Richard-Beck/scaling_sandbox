#!/usr/bin/env python3
"""Publish the paired-arm result pages at the familiar Figure 3 review path."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
GROUP = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019'
FIGURE = GROUP / 'figure3'
PAIRED = FIGURE / 'paired_stationary_review'
PDF = FIGURE / 'supporting/vectors/figure3_stationary_review.pdf'
METADATA = FIGURE / 'supporting/metadata'
ARCHIVE = FIGURE / 'preceding_angle_only_review_20261019'
MODELS = ('wave_pinning', 'goryachev', 'holmes', 'spring')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    return str(Path(path).resolve().relative_to(ROOT))


def main():
    generated = datetime.now(timezone.utc).isoformat()
    inputs = [PAIRED / 'supporting/vectors' / f'{model}_paired_vectors.pdf' for model in MODELS]
    assert all(p.is_file() for p in inputs), 'All four model pages must exist before publication'
    assert all(len(PdfReader(p).pages) == 1 for p in inputs)
    angle = FIGURE / 'standard_angles_sizes/supporting/data/stationary_sizes_summary.csv'
    arc = FIGURE / 'standard_arclength_sizes/supporting/data/stationary_arclength_summary.csv'
    for path in (angle, arc):
        table = pd.read_csv(path)
        assert len(table) == 240 and table.status.eq('complete').all()

    preservation = ARCHIVE / 'preservation.json'
    if not preservation.exists():
        assert sha(PDF) == '295b010570175d6c33f1680df0adb7a34f95123d166d0ba1c212282bf8b286d3'
        ARCHIVE.mkdir(parents=True, exist_ok=True)
        originals = [PDF, *METADATA.glob('stationary_review*')]
        records = []
        for old in originals:
            if not old.is_file():
                continue
            new = ARCHIVE / old.name
            shutil.copy2(old, new)
            assert sha(old) == sha(new)
            records.append({'original': rel(old), 'snapshot': rel(new), 'sha256': sha(new)})
        preservation.write_text(json.dumps({'preserved_utc': generated, 'files': records}, indent=2) + '\n')
        (ARCHIVE / 'README.md').write_text(
            '# Preceding angle-only Figure 3 review\n\n'
            'Exact seven-page packet and metadata preserved before publishing the paired-arm results.\n\n'
            '[Original PDF](figure3_stationary_review.pdf) · [Preservation hashes](preservation.json)\n'
        )
    for row in json.loads(preservation.read_text())['files']:
        assert sha(ROOT / row['snapshot']) == row['sha256']

    writer = PdfWriter()
    for page in inputs:
        writer.append(page)
    writer.add_metadata({'/Title': 'Figure 3: fixed-angle and fixed-initial-arclength stationary cues',
                         '/Subject': 'Four models; 480 arm rows, 355 unique protocols; withdrawal and retained activity vectors',
                         '/Author': 'Manuscript figure revision', '/Producer': Path(__file__).name})
    temporary = PDF.with_name(PDF.stem + '.paired.tmp.pdf')
    with temporary.open('wb') as stream:
        writer.write(stream)
    assert len(PdfReader(temporary).pages) == 4
    temporary.replace(PDF)
    source_files = [*inputs, angle, arc, Path(__file__)]
    for name in ('supporting/data/paired_vector_glyphs.csv', 'supporting/metadata/renderer_manifest.json'):
        path = PAIRED / name
        if path.is_file():
            source_files.append(path)
    record = {
        'schema_version': 1, 'generated_utc': generated,
        'document_stage': 'Paired-arm results for author discussion; native data independently validated; figure check in progress',
        'title': 'Figure 3: fixed-angle and fixed-initial-arclength stationary cues',
        'pdf': {'path': rel(PDF), 'sha256': sha(PDF), 'bytes': PDF.stat().st_size, 'pages': 4},
        'inputs': [{'path': rel(p), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in source_files],
        'displayed_arm_rows': 480, 'unique_protocols': 355, 'native_integrations_added_by_assembly': 0,
        'preceding_packet_preservation': rel(preservation),
        'old_stationary_review_json_unchanged': True, 'integrated_outline_unchanged': True,
    }
    (METADATA / 'stationary_review.build_manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    (METADATA / 'stationary_paired_review_delivery.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
