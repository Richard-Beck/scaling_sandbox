"""Verify an outline's six main figures, exact build inputs and readable PDF."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('document', type=Path)
    parser.add_argument('--require-tracked', action='store_true')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    document = json.loads(args.document.read_text())
    manifest_path = args.document.with_suffix('.build_manifest.json')
    manifest = json.loads(manifest_path.read_text())
    sections = document['sections']
    assert [s['id'] for s in sections] == [f'figure{i}' for i in range(1, 7)]
    for section in sections:
        for key in ['title', 'role', 'assets', 'legend', 'narrative', 'measurements', 'limitations', 'sources']:
            assert section[key], (section['id'], key)
    assert 'RESERVED' in sections[-1]['title']
    assert 'no observations' in ' '.join([sections[-1]['role']] + sections[-1]['legend'] + sections[-1]['limitations']).lower()
    visual_sections = sections + [s for s in document.get('supplementary', [])
                                  if isinstance(s, dict) and s.get('assets')]
    expected_assets = Counter(a['path'] for s in visual_sections for a in s['assets'])
    embedded_assets = Counter(a for p in manifest['page_inventory'] for a in p['assets'])
    assert expected_assets == embedded_assets, (expected_assets, embedded_assets)
    for record in manifest['inputs']:
        path = ROOT / record['path']
        assert path.is_file() and sha(path) == record['sha256'], record['path']
    pdf = ROOT / manifest['pdf']['path']
    assert sha(pdf) == manifest['pdf']['sha256']
    reader = PdfReader(pdf)
    assert len(reader.pages) == manifest['pdf']['pages']
    texts = [p.extract_text() or '' for p in reader.pages]
    assert all(t.strip() for t in texts)
    for section in sections:
        notes = [p for p in manifest['page_inventory']
                 if p['section_id'] == section['id'] and p['kind'] == 'figure_notes']
        assert notes, section['id']
        extracted = ' '.join(texts[p['number'] - 1] for p in notes)
        assert 'Detailed legend' in extracted and 'Narrative points' in extracted, section['id']
    assert all(p['bounds_valid'] and p['text_overlap_check'] == 'passed'
               for p in manifest['page_inventory'])
    checked_paths = {str(args.document.resolve().relative_to(ROOT)),
                     str(manifest_path.resolve().relative_to(ROOT)), manifest['pdf']['path']}
    checked_paths.update(r['path'] for r in manifest['inputs'])
    if args.require_tracked:
        tracked = set(subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0'))
        assert checked_paths <= tracked, sorted(checked_paths - tracked)
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
                  stage=document['stage'], main_figures=len(sections),
                  embedded_figure_sheets=sum(expected_assets.values()), pages=len(reader.pages),
                  exact_source_hashes_verified=len(manifest['inputs']),
                  all_figure_legends_and_narratives_verified=True,
                  pending_experiment_placeholder_verified=True,
                  all_build_inputs_tracked=args.require_tracked,
                  scope='Exact retained figure assets, document content and PDF packaging; scientific trial audits are separate.',
                  pdf_sha256=sha(pdf))
    if args.report:
        args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
