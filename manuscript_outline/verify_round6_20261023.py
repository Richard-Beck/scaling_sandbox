"""Check preservation, source coordinates, native PDF integration and local links."""
from pathlib import Path
from datetime import datetime, timezone
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
import hashlib
import json
import subprocess

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_round6_20261023'
DEV = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
FIGURES = DEV / 'revision_round6_20261023'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src') and value:
                self.links.append(value)


def main():
    starting = read(DEST / 'request/starting_checkpoint.json')
    for name, expected in starting['input_sha256'].items():
        actual = DEST / 'request/revision_brief_round6_20261023.md' if name == 'RevisionBrief.md' else ROOT / name
        assert sha(actual) == expected, name
    original_brief = subprocess.check_output(['git', 'show', starting['baseline_commit'] + ':RevisionBrief.md'], cwd=ROOT)
    assert original_brief == (DEST / 'request/revision_brief_round6_20261023.md').read_bytes()
    assert (ROOT / 'CurrentNarrative.md').read_bytes() == (DEST / 'request/current_narrative_20261023.md').read_bytes()
    baseline = read(ROOT / 'manuscript_outline/outline_v11.json')
    current = read(ROOT / 'manuscript_outline/outline_v12.json')
    for index in (4, 5):
        assert current['sections'][index] == baseline['sections'][index]
    assert next(s for s in current['supplementary'] if s.get('id') == 'pilot_observables') == next(s for s in baseline['supplementary'] if s.get('id') == 'pilot_observables')
    source_pairs = [(ROOT / 'polarity1d/manuscript/data/derived/length_correlations.csv',
                     FIGURES / 'figure1/supporting/data/E_length_correlations.csv')]
    for name in ('A_incumbent_organization.csv', 'B_initiation.csv', 'C_replacement_times_all450.csv',
                 'D_gap1_vs_pulse2_all180.csv', 'E_L39_challenger_share_traces.csv'):
        source_pairs.append((DEV / 'revision_round3_20261022/figure2/supporting/data' / name,
                             FIGURES / 'figure2/supporting/data' / name))
    for old, new in source_pairs:
        assert old.read_bytes() == new.read_bytes(), str(new)
    unchanged = subprocess.run(['git', 'diff', '--exit-code', starting['baseline_commit'], '--',
        'author', 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round5_20261023',
        'manuscript_outline/outline_v11.pdf', 'manuscript_outline/outline_v11.json'], cwd=ROOT, capture_output=True, text=True)
    assert unchanged.returncode == 0, unchanged.stdout
    manifest = read(ROOT / 'manuscript_outline/outline_v12.build_manifest.json')
    pdf = PdfReader(ROOT / 'manuscript_outline/outline_v12.pdf')
    native_count = 0
    for page in manifest['page_inventory']:
        if page.get('native_vector_pdf'):
            original = PdfReader(ROOT / page['native_vector_pdf']).pages[0]
            integrated = pdf.pages[page['number']-1]
            assert original.get_contents().get_data() == integrated.get_contents().get_data()
            native_count += 1
    sources = 0
    for section in current['sections'][:4]:
        folder = FIGURES / section['id']
        assert (folder / 'scientific_readme.md').read_text().startswith('# Unexpected')
        for source in section.get('main_source_data', []):
            assert sha(ROOT / source['path']) == source['sha256']
            sources += 1
        assert Path(ROOT / section['main_pdf']).is_file()
        assert Path(ROOT / section['supplement_pdf']).is_file()
    links_checked = 0
    for path in DEST.glob('*.html'):
        parser = Links()
        parser.feed(path.read_text())
        for link in parser.links:
            parsed = urlsplit(link)
            if not parsed.scheme and parsed.path:
                target = (path.parent / unquote(parsed.path)).resolve()
                assert target.exists(), (path, link)
                links_checked += 1
    preservation = DEST / 'review_preservation.json'
    if preservation.exists():
        for item in read(preservation)['reviews']:
            assert sha(ROOT / item['path']) == item['sha256'], item['path']
    report = {'generated_utc': datetime.now(timezone.utc).isoformat(), 'passed': True,
        'baseline_commit': starting['baseline_commit'], 'baseline_pdf_preserved': True,
        'round5_files_and_author_documents_preserved': True, 'protected_sections': ['figure5', 'figure6', 'pilot_observables'],
        'byte_identical_figure1_2_data_tables': len(source_pairs),
        'main_source_hashes_checked': sources, 'native_pdf_pages_preserved': native_count,
        'local_html_links_checked': links_checked,
        'review_preservation_checked': preservation.exists(),
        'scope': 'Source and packaging audit; model-specific numerical and scientific checks are recorded with Figures 3–4.'}
    (DEST / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
