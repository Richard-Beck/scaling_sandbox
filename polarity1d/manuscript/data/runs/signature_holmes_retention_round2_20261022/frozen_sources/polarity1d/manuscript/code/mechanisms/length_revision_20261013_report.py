"""Assemble the annotated deBelly revision from checked, retained evidence.

The chemical chapters and their assets stay at the 10 October revision.
New deBelly chapters supply explicit bound-pool budgets and one illustrative cue.
No simulation inputs or native outcomes are changed by this renderer.
"""
from pathlib import Path
from datetime import datetime, timezone
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
REV = FIG / 'revision_20261013'
DATA = ROOT / 'data/runs/length_mechanisms_20261008/revision_20261013'
PRIOR = FIG / 'archive_20261013_before_debelly_revision'
TITLE = 'Mechanisms of cell-length-dependent polarity reversal'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def previous_renderer():
    path = Path(__file__).with_name('length_revision_20261010_report.py')
    spec = importlib.util.spec_from_file_location('retained_length_report_layout', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.REV, module.DATA = REV, DATA
    return module


def chapters():
    old = read(PRIOR / 'chapters.json')
    chemical = copy.deepcopy(old[:7])
    sequence = read(REV / 'debelly_sequence/chapters.json')
    by_number = {int(c['figure_label']): c for c in sequence}
    by_number[9] = read(REV / 'debelly_figure9/chapter.json')
    by_number[12] = read(REV / 'conclusion/chapter.json')
    assert set(by_number) == {8, 9, 10, 11, 12}, sorted(by_number)
    result = chemical + [by_number[n] for n in sorted(by_number)]
    for n, chapter in enumerate(result, 1):
        chapter['figure_label'] = str(n)
        chapter['id'] = str(chapter['id'])
        chapter.setdefault('manifest', str(Path(chapter['image']).with_suffix('.manifest.json')))
        chapter.setdefault('sources', [])
        chapter['sources'] = [[item['href'], item['label']] if isinstance(item, dict) else item
                              for item in chapter['sources']]
        if not chapter.get('caption'):
            path = FIG / chapter['image']
            chapter['caption'] = path.with_name(path.stem + '_caption.txt').read_text().strip()
        chapter['caption'] = re.sub(r'^\*?Figure\s+[^ .\n*]+\*?[. ]*',
                                    f'Figure {n} ', chapter['caption'], count=1)
    assert result[:7] == chemical
    assert len({c['id'] for c in result}) == len(result)
    return result


def verify_archive():
    records = read(PRIOR / 'preservation.json')
    for name, record in records.items():
        assert sha(PRIOR / name) == record['sha256'], name
    return {name: record['sha256'] for name, record in records.items()}


def html_report(content, stamp, evidence):
    esc = html.escape
    nav = ''.join(f'<a href="#{esc(c["id"], quote=True)}">{n}. {esc(c["title"])}</a>'
                  for n, c in enumerate(content, 1))
    sections = []
    for n, c in enumerate(content, 1):
        paragraphs = ''.join('<p>' + esc(p) + '</p>' for p in c['paragraphs'])
        prefix = 'Figure ' + str(n)
        legend = '<em>' + prefix + '</em> ' + esc(c['caption'][len(prefix):].lstrip(' .'))
        links = [(str(Path(c['image']).with_suffix('.pdf')), 'Vector figure'),
                 (c['manifest'], 'Data and provenance')]
        method_path = (FIG / c['image']).with_name(Path(c['image']).stem + '_methods.txt')
        if method_path.is_file():
            links.append((str(method_path.relative_to(FIG)), 'Technical notes'))
        links.extend((url, label) for url, label in c['sources'] if url not in {x[0] for x in links})
        linked = ' · '.join(f'<a href="{esc(url, quote=True)}">{esc(label)}</a>' for url, label in links)
        sections.append(f'<section id="{esc(c["id"], quote=True)}"><div class="number">Figure {n}</div>'
                        f'<h2>{esc(c["title"])}</h2><div class="prose">{paragraphs}</div>'
                        f'<figure><a href="{esc(c["image"], quote=True)}"><img '
                        f'src="{esc(c["image"], quote=True)}" alt="{esc(c["caption"], quote=True)}" '
                        f'loading="{"eager" if n in [1, 8, 9] else "lazy"}"></a>'
                        f'<figcaption>{legend}</figcaption></figure><p class="sources">{linked}</p></section>')
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{TITLE}</title><style>
:root{{--ink:#173c49;--muted:#586e78;--accent:#197987;--rule:#dce7e9;--wash:#f3f7f7}}
*{{box-sizing:border-box}}body{{margin:0;color:var(--ink);background:white;font:17px/1.65 system-ui,sans-serif}}
main{{max-width:1270px;padding:46px 38px 64px;margin:auto}}header{{border-bottom:2px solid var(--ink);padding-bottom:22px}}
h1{{font-size:clamp(29px,4vw,41px);line-height:1.15;max-width:1000px}}h2{{font-size:clamp(24px,3vw,32px);line-height:1.23;margin:9px 0 20px;max-width:1050px}}
.eyebrow,.number{{font-size:12px;letter-spacing:.06em;font-weight:650;text-transform:uppercase;color:var(--muted)}}
p{{margin:0 0 16px}}.intro{{max-width:950px;font-size:19px}}.metadata,.sources{{font-size:13px;color:var(--muted)}}
a{{color:var(--accent);text-underline-offset:3px}}.download{{display:inline-block;padding:9px 17px;background:var(--ink);color:white;border-radius:4px;text-decoration:none;margin:6px 14px 16px 0}}
nav{{display:grid;grid-template-columns:1fr 1fr;gap:7px 28px;padding:21px 25px;margin:26px 0;background:var(--wash);font-size:14px}}nav a{{text-decoration:none}}
section{{padding:36px 0 30px;border-bottom:1px solid var(--rule)}}.prose{{max-width:980px}}figure{{margin:25px 0 15px}}img{{display:block;width:100%;height:auto}}
figcaption{{font-size:14px;line-height:1.55;color:#425963;max-width:1110px;margin-top:15px}}details{{padding:22px 25px;margin-top:30px;background:var(--wash);font-size:14px}}
summary{{cursor:pointer;font-weight:650}}details p{{margin-top:13px}}footer{{font-size:13px;color:var(--muted);margin-top:22px}}
@media(max-width:700px){{main{{padding:25px 18px}}nav{{grid-template-columns:1fr}}body{{font-size:16px}}}}
@media print{{nav,.download,details{{display:none}}main{{padding:0}}section{{break-before:page}}}}
</style></head><body><main><header><div class="eyebrow">Polarity manuscript · revision 13 October 2026</div>
<h1>{TITLE}</h1><p class="intro">Cell length changes the spatial chemical distributions that must be reorganized and the reactions and mechanical processes that reorganize them. The deBelly section follows the local bound-MCA balance, the redistribution attained during a reversing cue, and the persistent old-end loss and later MCA rebound in the longest cell.</p>
<a class="download" href="length_mechanisms.pdf">Report PDF · {len(content)} vector pages</a><a href="report_text.txt">Plain-text reader</a>
<p class="metadata">Report generated {esc(stamp)}. Core comparison generated {esc(evidence['comparison_generated_at'])}. Native-data producers finished. Figure manifests record source and diagnostic-generation times; opening or manually refreshing this page loads the latest generated report.</p></header>
<nav aria-label="Figures">{nav}</nav>{''.join(sections)}
<details><summary>Measurements and reproduction</summary>
<p>The illustrative deBelly reversals use a raised-cosine patch supported on the rightmost 20% of nominal length, with peak 3.6, applied to each length’s exact basal checkpoint. Figures comparing cue geometries label the fixed-width and constant-slope cases separately. The selected pulse is the shortest tested duration passing the canonical 1% Rac scalar-strength gate after 10000 s of basal recovery. Recovery time is excluded from that duration. Scalar success and the separate native-field tail flag retain their distinct definitions.</p>
<p>Bound MCA follows the cortex; unbound MCA follows the membrane. The bound-pool budget separates transport of an existing concentration gradient, dilution or concentration by a flow gradient, bound diffusion, binding minus unbinding, and the uniform mass correction. Exact native-step budgets retain the solver’s staggering, upwind donor selection and pre-correction concentrations. The total-pool exchange residual is a consequence of sequential updates and is identified separately from physical binding exchange.</p>
<p>Centroid displacement measures redistribution of a changing pool. It includes removal at boundaries, binding and uniform replacement, and does not measure the path length of molecules. Reflected starting profiles supply a geometric reference, not a measured universal withdrawal threshold. Mechanical source allocations use the actual coupled state and do not represent independently evolved pathway removals. The report uses corrected independent deBelly endpoint commands.</p>
<p><a href="../../THREE_TESTS.md">Protocol</a> · <a href="../../ASSUMPTIONS.md">Model assumptions</a> · <a href="revision_20261013/report_provenance.json">Report provenance</a> · <a href="revision_20261013/editorial_review/annotation_response.json">Annotation response</a> · <a href="archive_20261013_before_debelly_revision/reader.html">Preserved 10 October reader</a></p>
<p>Rebuild with <code>/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_revision_20261013_report.py</code>. Dedicated figure generators and full source hashes are recorded in the manifests.</p>
</details><footer>Bound-MCA dynamics and total-pool conservation are reported separately.</footer></main></body></html>'''
    (FIG / 'index.html').write_text(document)


def audit(content):
    reader = PdfReader(FIG / 'length_mechanisms.pdf')
    assert len(reader.pages) == len(content) == 12
    searchable = []
    for n, (page, c) in enumerate(zip(reader.pages, content), 1):
        text = page.extract_text()
        assert f'FIGURE {n}' in text and f'Figure {n}' in text
        assert len(text) > 450 and '/Font' in page['/Resources']
        assert re.match(r'^Figure ' + str(n) + r' ', c['caption'])
        searchable.append(dict(figure=n, characters=len(text)))
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
    old = read(PRIOR / 'chapters.json')
    assert content[:7] == old[:7], 'Chemical chapters changed'
    return dict(passed=True, pages=len(content), vector_figures=True, searchable_text=searchable,
                local_links_checked=len(destinations), no_timed_refresh=True,
                unchanged_chemical_chapters=True, page_bounds_passed=True)


def main():
    REV.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    preserved = verify_archive()
    stamp = datetime.now(timezone.utc).isoformat()
    content = chapters()
    base = previous_renderer()
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
                      renderer_sha256=sha(Path(__file__)), layout_renderer='code/mechanisms/length_revision_20261010_report.py',
                      layout_renderer_sha256=sha(Path(__file__).with_name('length_revision_20261010_report.py')),
                      **evidence, checks=checks, chapters=content, preserved_prior_report=preserved,
                      pdf_sha256=sha(FIG / 'length_mechanisms.pdf'), html_sha256=sha(FIG / 'index.html'),
                      text_sha256=sha(FIG / 'report_text.txt'))
    provenance['review_inputs'] = {
        str(path.relative_to(ROOT)): sha(path) for path in [
            ROOT / 'AGENTS.md',
            ROOT / 'figures/user_annotated_reports/length_mechanisms (2c).pdf',
            REV / 'editorial_review/annotated_comments.json',
            REV / 'editorial_review/annotation_response.json',
        ]
    }
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
