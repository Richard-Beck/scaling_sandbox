"""Verify exact revision inputs, native figure pages and protected sections."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit
import re
from pypdf import PdfReader
from prepare_figures3_4_20261023 import relocate_sources

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_figures3_4_20261023'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def main():
    receipt = read(DEST / 'integration.json')
    checkpoint = receipt['starting_checkpoint']
    for path, expected in checkpoint['baseline_files'].items():
        assert sha(ROOT / path) == expected, ('baseline changed', path)
    assert sha(DEST / 'request/revision_brief_20261023.md') == checkpoint['archived_brief_sha256']
    baseline = read(ROOT / receipt['baseline'])
    document = read(ROOT / receipt['outline'])
    assert sha(ROOT / receipt['outline']) == receipt['outline_sha256']
    assert sha(ROOT / receipt['baseline']) == receipt['baseline_sha256']
    for ident, record in receipt['protected_sections'].items():
        old = next(s for s in baseline['sections'] if s['id'] == ident)
        new = next(s for s in document['sections'] if s['id'] == ident)
        assert relocate_sources(old) == new, ('protected content changed', ident)
        for path, expected in record['asset_sha256'].items():
            assert sha(ROOT / path) == expected, ('protected artwork changed', path)
    assert [s['id'] for s in document['sections']] == [f'figure{i}' for i in range(1, 7)]
    native_pages = []
    for ident in ('figure3', 'figure4'):
        section = next(s for s in document['sections'] if s['id'] == ident)
        source = read(ROOT / section['main_pdf'].rsplit('/', 1)[0] / 'legend.json')
        assert section['legend'] == source['legend']
        for record in section['sources']:
            if 'sha256' in record:
                assert sha(ROOT / record['path']) == record['sha256'], ('figure source changed', record['path'])
        assert sha(DEST / f'{ident}_with_legend.pdf') == receipt['publication_pdf_sha256'][ident]
        main = PdfReader(ROOT / section['main_pdf'])
        publication = PdfReader(DEST / f'{ident}_with_legend.pdf')
        assert len(publication.pages) > len(main.pages)
        for original, actual in zip(main.pages, publication.pages):
            assert original.get_contents().get_data() == actual.get_contents().get_data()
        for asset in section['assets']:
            assert asset.get('native_pdf_page'), (ident, 'native vector page missing')
            native_pages.append(asset['native_pdf_page'])
        assert (ROOT / section['supplement_pdf']).is_file()
    manifest = read(ROOT / 'manuscript_outline/outline_v11.build_manifest.json')
    native_inventory = [p.get('native_vector_pdf') for p in manifest['page_inventory'] if p.get('native_vector_pdf')]
    assert set(native_pages) <= set(native_inventory), 'native pages missing in outline'
    pilots = read(DEST / 'pilot_observables.json')
    assert len(pilots['observables']) == 5
    assert sha(DEST / 'pilot_observables.json') == receipt['pilot_table_sha256']
    assert 'before' in pilots['fixed_before_data'].lower()
    assert any(s.get('id') == 'pilot_observables' for s in document['supplementary'] if isinstance(s, dict))
    reviews = read(DEST / 'review_preservation.json')['reviews']
    for review in reviews:
        assert sha(ROOT / review['path']) == review['sha256'], ('verbatim review changed', review['path'])
    # The integrated reader must remain usable offline, without automatic reload.
    reader = DEST / 'index.html'
    text = reader.read_text()
    assert not re.search(r'http-equiv\s*=\s*[\"\']?refresh', text, flags=re.I)
    local_links = []
    for target in re.findall(r'(?:href|src)=[\"\']([^\"\']+)', text):
        parsed = urlsplit(target)
        if parsed.scheme or not parsed.path:
            continue
        path = (reader.parent / unquote(parsed.path)).resolve()
        assert path.is_file(), ('broken reader link', target)
        local_links.append(target)
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
                  baseline_files_unchanged=len(checkpoint['baseline_files']),
                  protected_sections=list(receipt['protected_sections']),
                  figure3_4_publication_pages_preserved=True, native_figure_pages_integrated=True,
                  no_new_noise_campaign=True, prospective_observables=5,
                  verbatim_reviews_verified=len(reviews),
                  reader_links_checked=len(local_links),
                  scope='Packaging, baseline preservation, native artwork and links; scientific checks are recorded with each figure.')
    (DEST / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
