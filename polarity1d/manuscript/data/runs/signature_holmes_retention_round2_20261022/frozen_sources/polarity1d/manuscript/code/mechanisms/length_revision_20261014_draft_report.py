"""Align the mechanochemical report with the author's 14 October draft.

Rendering only. Read retained native evidence and new endpoint figures; preserve
the preceding report and the first seven chapters. Equations and narrative are
stored together in a reviewable chapter override file.
"""
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import copy
import hashlib
import html
import importlib.util
import json
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'figures/length_mechanisms_20261008'
REV = FIG / 'revision_20261014_draft_alignment'
DATA = ROOT / 'data/runs/length_mechanisms_20261008/revision_20261014_draft_alignment'
PRIOR = FIG / 'archive_20261014_before_draft_alignment'
TITLE = 'Mechanisms of cell-length-dependent polarity reversal'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def layout_renderer():
    path = Path(__file__).with_name('length_revision_20261010_report.py')
    spec = importlib.util.spec_from_file_location('retained_length_report_layout', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.REV, module.DATA = REV, DATA
    equations = equation_renderer()
    module.pdf_report = lambda content, stamp: equations.pdf_report(content, stamp, module)
    return module


@lru_cache(maxsize=1)
def equation_renderer():
    path = Path(__file__).with_name('length_revision_20261014_equations.py')
    spec = importlib.util.spec_from_file_location('publication_equations', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def chapters():
    result = copy.deepcopy(read(PRIOR / 'chapters.json'))
    overrides = read(REV / 'chapter_overrides.json')
    assert set(overrides) == {'8', '9', '10', '11', '12'}
    for number, override in overrides.items():
        n = int(number)
        result[n - 1].update(override)
    for n, chapter in enumerate(result, 1):
        chapter['figure_label'] = str(n)
        chapter['caption'] = re.sub(r'^\*?Figure\s+[^ .\n*]+\*?[. ]*',
                                   f'Figure {n} ', chapter['caption'], count=1)
        chapter.setdefault('sources', [])
    assert len(result) == 12 and len({c['id'] for c in result}) == 12
    assert result[:7] == read(PRIOR / 'chapters.json')[:7]
    return result


def verify_archive():
    records = read(PRIOR / 'preservation.json')
    for name, record in records.items():
        assert sha(PRIOR / name) == record['sha256'], name
    return records


def html_report(content, stamp, evidence):
    esc = html.escape
    prior_document = (PRIOR / 'index.html').read_text()
    style = re.search(r'<style>(.*?)</style>', prior_document, re.S).group(1)
    style += '''
.equation-display{position:relative;margin:22px 0 15px;break-inside:avoid}
img.equation-math{width:100%;height:auto}
.equation-number{position:absolute;right:0;top:50%;transform:translateY(-50%);font-size:16px}
math{font-family:"STIX Two Math","Cambria Math",serif}
'''
    nav = ''.join(f'<a href="#{esc(c["id"], quote=True)}">{n}. {esc(c["title"])}</a>'
                  for n, c in enumerate(content, 1))
    sections = []
    for n, c in enumerate(content, 1):
        paragraphs = ''.join(equation_renderer().html_equation(p) if p.startswith('Eq. ')
                             else '<p>' + esc(p) + '</p>' for p in c['paragraphs'])
        prefix = f'Figure {n}'
        legend = '<em>' + prefix + '</em> ' + esc(c['caption'][len(prefix):].lstrip(' .'))
        links = [(str(Path(c['image']).with_suffix('.pdf')), 'Vector figure'),
                 (c['manifest'], 'Data and provenance')]
        methods = (FIG / c['image']).with_name(Path(c['image']).stem + '_methods.txt')
        if methods.is_file():
            links.append((str(methods.relative_to(FIG)), 'Technical notes'))
        links.extend((url, label) for url, label in c['sources'] if url not in {p[0] for p in links})
        linked = ' · '.join(f'<a href="{esc(url, quote=True)}">{esc(label)}</a>' for url, label in links)
        sections.append(f'<section id="{esc(c["id"], quote=True)}"><div class="number">Figure {n}</div>'
                        f'<h2>{esc(c["title"])}</h2><div class="prose">{paragraphs}</div>'
                        f'<figure><a href="{esc(c["image"], quote=True)}"><img '
                        f'src="{esc(c["image"], quote=True)}" alt="{esc(c["caption"], quote=True)}" '
                        f'loading="{"eager" if n in [1, 8, 9] else "lazy"}"></a>'
                        f'<figcaption>{legend}</figcaption></figure><p class="sources">{linked}</p></section>')
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{TITLE}</title><style>{style}</style></head><body><main><header>
<div class="eyebrow">Polarity manuscript · author-draft revision 14 October 2026</div>
<h1>{TITLE}</h1><p class="intro">Cell length changes the chemical distributions that must be reorganized and the processes that reorganize them. The mechanochemical section follows Rac and Rho activation, cortical motion and local membrane–cortex attachment redistribution, then examines how the initial attachment state and its remodelling depend on cell length.</p>
<a class="download" href="length_mechanisms.pdf">Report PDF · {len(content)} vector figure pages</a><a href="report_text.txt">Plain-text reader</a>
<p class="metadata">Report generated {esc(stamp)}. Core comparison generated {esc(evidence['comparison_generated_at'])}. Native-data producers finished. Figure manifests identify their source and data-generation times. Opening or manually refreshing this page loads the latest generated report.</p></header>
<nav aria-label="Figures">{nav}</nav>{''.join(sections)}
<details><summary>Measurements and reproduction</summary>
<p>Illustrative deBelly reversals use the rightward relative-20% raised-cosine Rac patch, with peak 3.6, applied to each length’s exact corrected basal checkpoint. Other geometries are labelled separately. The selected duration is the shortest tested pulse passing the canonical 1% Rac scalar-strength gate after 10000 s of basal recovery. This duration excludes recovery. The separate native-field tail flag is retained in the figures and exported data.</p>
<p>Bound MCA follows the cortex and unbound MCA follows the membrane. The bound-pool balance distinguishes profile advection, dilution or concentration by the flow gradient, bound diffusion, net binding and uniform replacement. Exact native budgets retain the staggered updates, finite-volume upwind donors and pre-correction concentrations. In the native sequential implementation, the correction replaces the actual mass deficit, including a numerical exchange residual. Its rate therefore differs from boundary outflux divided by length. The residual is documented in the linked global-budget companion.</p>
<p>Endpoint MCA denotes the first or last of the 201 native cells. The near-front regional budgets in Figure 12 average the rightmost 40 cells; they are a different observable. Comparisons at 30, 60, 120 and 180 s use common elapsed times, while withdrawal values use each length’s own selected duration. The coupled native trajectories support an endpoint-remodelling interpretation, without defining a universal MCA shutdown concentration or causally partitioning the cue-duration threshold.</p>
<p>Mechanical source components are superposed at the same coupled state. They do not independently evolve a removed pathway. The corrected model prescribes independent cortical endpoint commands; membrane endpoint velocity is a distinct measured quantity. Primary-paper context and the full equations are in <a href="revision_20261014_draft_alignment/equations_and_paper.md">equations and sources</a>. Earlier centroid and short-cell stability diagnostics remain in the <a href="archive_20261014_before_draft_alignment/reader.html">preserved previous reader</a>.</p>
<p><a href="../../THREE_TESTS.md">Protocol</a> · <a href="../../ASSUMPTIONS.md">Model assumptions</a> · <a href="revision_20261014_draft_alignment/report_provenance.json">Report provenance and audit</a> · <a href="revision_20261014_draft_alignment/revision_notes.md">Revision notes</a> · <a href="revision_20261014_draft_alignment/mechanism_evidence.md">Endpoint mechanism evidence</a></p>
<p>Rebuild with <code>/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_revision_20261014_draft_report.py</code>. Figure renderers and source hashes are recorded in the manifests.</p>
</details><footer>Retained native simulations; local MCA balances and endpoint remodelling.</footer></main></body></html>'''
    (FIG / 'index.html').write_text(document)


def audit(content):
    reader = PdfReader(FIG / 'length_mechanisms.pdf')
    assert len(reader.pages) == len(content) == 12
    searchable = []
    for n, (page, chapter) in enumerate(zip(reader.pages, content), 1):
        text = page.extract_text()
        assert f'FIGURE {n}' in text and f'Figure {n}' in text
        assert len(text) > 450 and '/Font' in page['/Resources']
        assert re.match(r'^Figure ' + str(n) + r' ', chapter['caption'])
        searchable.append(dict(figure=n, characters=len(text)))
    prior = read(PRIOR / 'chapters.json')
    assert content[:7] == prior[:7]
    assert not re.search(r'centroid|cell midpoint|cortical midpoint',
                         ' '.join(p for c in content[7:] for p in c['paragraphs']), re.I)
    document = (FIG / 'index.html').read_text()
    assert 'http-equiv="refresh"' not in document.lower()
    destinations = re.findall(r'(?:href|src)="([^"#]+)"', document)
    for destination in destinations:
        if not destination.startswith(('https:', 'http:', 'mailto:')):
            assert (FIG / html.unescape(destination).split('#', 1)[0]).exists(), destination
    with tempfile.TemporaryDirectory(prefix='length-report-bounds-') as directory:
        path = Path(directory) / 'bounds.html'
        subprocess.run(['pdftotext', '-bbox', str(FIG / 'length_mechanisms.pdf'), str(path)], check=True)
        xml = ET.parse(path)
        bounds, violations = [], []
        for n, page in enumerate(xml.findall('.//{*}page'), 1):
            width, height = float(page.attrib['width']), float(page.attrib['height'])
            words = page.findall('.//{*}word')
            for word in words:
                a = word.attrib
                if float(a['xMin']) < -0.5 or float(a['yMin']) < -0.5 or float(a['xMax']) > width + 0.5 or float(a['yMax']) > height + 0.5:
                    violations.append(dict(figure=n, text=word.text, bounds=a))
            bounds.append(dict(figure=n, words=len(words), width=width, height=height))
        assert not violations, violations
        write(REV / 'page_bounds_audit.json', dict(passed=True, pages=bounds, violations=violations))
    equations = [p for c in content[7:] for p in c['paragraphs'] if p.startswith('Eq. ')]
    assert equations
    assert all('(' + p.split(' ', 2)[1].rstrip('.') + ')' in reader.pages[7 if p in content[7]['paragraphs'] else 8].extract_text()
               for p in equations)
    layout = read(REV / 'equation_typesetting/equation_layout.json')
    assert [r['number'] for r in layout['equations']] == list(range(1, 9))
    assert len({r['number_right_points'] for r in layout['equations']}) == 1
    assert all(row['center_x_points'] == 421.2 for r in layout['equations'] for row in r['rows'])
    previous = PdfReader(FIG / 'archive_20261014_before_equation_typesetting/length_mechanisms.pdf')
    preserved_pages = [n for n in range(1, 13) if n not in {8, 9}]
    assert all(reader.pages[n-1].get_contents().get_data() == previous.pages[n-1].get_contents().get_data()
               for n in preserved_pages)
    return dict(passed=True, pages=len(content), vector_figures=True, searchable_text=searchable,
                local_links_checked=len(destinations), no_timed_refresh=True,
                unchanged_chemical_chapters=True, endpoint_narrative=True,
                equation_paragraphs=len(equations), page_bounds_passed=True,
                equations_centred=True, equation_numbers_right_aligned=True,
                publication_math_font=layout['font'], preserved_pdf_pages=preserved_pages)


def main():
    REV.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    preserved = verify_archive()
    archived_reader = (PRIOR / 'index.html').read_text().replace('<head>', '<head><base href="../">', 1)
    for name in ['length_mechanisms.pdf', 'report_text.txt']:
        archived_reader = archived_reader.replace('href="' + name + '"',
            'href="archive_20261014_before_draft_alignment/' + name + '"')
    (PRIOR / 'reader.html').write_text(archived_reader)
    stamp = datetime.now(timezone.utc).isoformat()
    content = chapters()
    base = layout_renderer()
    evidence = base.verify(content)
    write(REV / 'chapters.json', content)
    write(REV / 'report_provenance.json', {})
    base.pdf_report(content, stamp)
    html_report(content, stamp, evidence)
    plain = [TITLE, 'Generated ' + stamp + ' · native-data producers finished', '']
    for c in content:
        plain.extend(['Figure ' + c['figure_label'] + ' — ' + c['title'], '',
                      '\n\n'.join(c['paragraphs']), '', c['caption'], ''])
    (FIG / 'report_text.txt').write_text('\n'.join(plain))
    checks = audit(content)
    provenance = dict(generated_at=stamp, producer='finished', renderer=str(Path(__file__).relative_to(ROOT)),
                      renderer_sha256=sha(Path(__file__)), **evidence, checks=checks, chapters=content,
                      preserved_prior_report=preserved,
                      pdf_sha256=sha(FIG / 'length_mechanisms.pdf'), html_sha256=sha(FIG / 'index.html'),
                      text_sha256=sha(FIG / 'report_text.txt'),
                      revision_scope=dict(updated_figures=[9, 10], revised_prose=[8, 9, 10, 11, 12],
                                          unchanged_figures=list(range(1, 8)), native_integrations=0))
    provenance['equation_typesetting'] = dict(
        renderer='code/mechanisms/length_revision_20261014_equations.py',
        renderer_sha256=sha(Path(__file__).with_name('length_revision_20261014_equations.py')),
        original_equations_and_prose_unchanged=True, reformatted_pages=[8, 9],
        layout_sha256=sha(REV / 'equation_typesetting/equation_layout.json'),
        svg_sha256={p.name: sha(p) for p in sorted((REV / 'equation_typesetting').glob('equation_*.svg'))},
        preserved_pre_typesetting_report=read(FIG / 'archive_20261014_before_equation_typesetting/preservation.json'))
    for path in [REV / 'report_provenance.json', DATA / 'report_provenance.json',
                 ROOT / 'data/runs/length_mechanisms_20261008/report_provenance.json']:
        write(path, provenance)
    write(DATA / 'report_audit.json', checks)
    old = ROOT / 'figures/three_test_mechanisms/three_test_mechanisms.pdf'
    if old.is_file():
        shutil.copy2(FIG / 'length_mechanisms.pdf', old)
    print(json.dumps(dict(pdf=str(FIG / 'length_mechanisms.pdf'), html=str(FIG / 'index.html'),
                          generated_at=stamp, pages=checks['pages'], passed=True), indent=2))


if __name__ == '__main__':
    main()
