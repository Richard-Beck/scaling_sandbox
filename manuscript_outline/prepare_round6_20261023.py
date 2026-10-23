"""Integrate round6 into the current outline_v11, preserving its later sections."""
from copy import deepcopy
import argparse
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import html
import json
import os
import re
import sys

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import math_to_image
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_round6_20261023'
FIGURES = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023'
sys.path.insert(0, str(ROOT / 'scripts'))
import build_manuscript_outline as prose


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def relative(path):
    return str(path.relative_to(ROOT))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def legend_copy(section):
    ident = section['id']
    path = DEST / f'{ident}_legend.pdf'
    builder = prose.OutlineBuilder({'title': section['title'], 'sections': []}, ROOT, path, FIGURES / ident / 'legend.json')
    with PdfPages(path) as book:
        builder.pdf = book
        builder.section_id = ident
        builder.new_page(section['title'], kind='figure_notes')
        for paragraph in section['legend']:
            builder.paragraph(paragraph, size=10.2)
        for equation in section.get('equations', []):
            builder.equation(equation)
        builder.finish_page()
    writer = PdfWriter()
    writer.append(ROOT / section['main_pdf'])
    writer.append(path)
    writer.write(DEST / f'{ident}_with_legend.pdf')


def html_equations(ident, equations):
    folder = DEST / 'equations'
    folder.mkdir(exist_ok=True)
    body = ''
    for index, equation in enumerate(equations):
        path = folder / f'{ident}_{index+1}.png'
        math_to_image('$' + equation['math'].strip('$') + '$', path,
                      prop=FontProperties(size=16), dpi=180, format='png', color='black')
        number = html.escape(str(equation.get('number', index+1)))
        body += f'<div class="equation"><img src="equations/{path.name}" alt="Equation {number}"><span>({number})</span></div>'
        body += '<p>' + html.escape(equation.get('description', '')) + '</p>'
    return body


def update_navigation(sections, review_candidate=False):
    marker = '<!-- round6-current-edition -->'
    end = '<!-- round6-current-edition end -->'
    for path, prefix in ((ROOT / 'manuscript_outline/README.md', ''),
                         (ROOT / 'current manuscript figs/README.md', '../manuscript_outline/')):
        status = ' Review candidate; the targeted mode campaign is active.' if review_candidate else ''
        intro = f'''{marker}
# Spatial susceptibility and heading persistence

[**Current integrated outline_v12**]({prefix}outline_v12.pdf) · [**Figures 1–4 and legends**]({prefix}revision_round6_20261023/figures1_4_with_legends.pdf) · [**Reader and supporting evidence**]({prefix}revision_round6_20261023/index.html)

Figures 1–2 retain their numerical evidence with separate letters for distinct
quantities. Figure 3 tests the native spatial susceptibility map away from the
edge and exposes its signed contributions to domain drift. Figure 4 shows
all-input closure bounds, size-specific remainder curves and per-model takeover
comparisons; Fourier tests and exact-region-moment checks are recorded alongside.
Figure 5, reserved Figure 6 and the prospective pilot appendix are preserved.{status}

{end}
'''
        old = path.read_text()
        if marker in old:
            before, rest = old.split(marker, 1)
            _, after = rest.split(end, 1)
            path.write_text(before + intro + after)
        else:
            path.write_text(intro + '\n' + old)
    path = ROOT / 'README.md'
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith('[**Current manuscript outline'):
            lines[index] = '[**Current manuscript outline — 23 October 2026**](manuscript_outline/outline_v12.pdf) · [**Revised Figures 1–4 and legends**](manuscript_outline/revision_round6_20261023/figures1_4_with_legends.pdf) · [**Revision reader**](manuscript_outline/revision_round6_20261023/index.html).'
            break
    path.write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--figure4-snapshot', type=Path,
        help='Use an immutable, explicitly interim Figure 4 review snapshot.')
    args = parser.parse_args()
    DEST.mkdir(parents=True, exist_ok=True)
    baseline_path = ROOT / 'manuscript_outline/outline_v11.json'
    baseline = read(baseline_path)
    document = deepcopy(baseline)
    sections = {}
    for n in range(1, 5):
        ident = f'figure{n}'
        source_folder = args.figure4_snapshot.resolve() if n == 4 and args.figure4_snapshot else FIGURES / ident
        section = read(source_folder / 'legend.json')
        if n == 4 and args.figure4_snapshot:
            for asset in section['assets']:
                for key in ('path', 'native_pdf_page'):
                    if key in asset:
                        asset[key] = relative(source_folder / Path(asset[key]).name)
            for key in ('main_pdf', 'supplement_pdf'):
                section[key] = relative(source_folder / Path(section[key]).name)
            section['reader'] = relative(source_folder / 'scientific_readme.md')
            section['sources'] = [
                {'label': 'Snapshot methods', 'path': relative(source_folder / 'methods_provenance.md')},
                {'label': 'Snapshot results and active-campaign status', 'path': section['reader']}]
            section['status'] = 'review_candidate_mode_campaign_active'
            section['narrative'] = ['Larger rings preserve heading longer while their local size exponent weakens as input strengthens. '
                'Completed region diagnostics show all-input closure bounds and remainder curves. '
                'The independent modes-only campaign is active; this review candidate claims no completed mode-collapse result.']
            section['limitations'].append('The targeted modes-only campaign is active in this immutable review snapshot; final precision results will be integrated after it finishes.')
        # Full LaTeX permits a one-token math-font argument; Matplotlib's
        # publication renderer requires braces around that same token.
        for equation in section.get('equations', []):
            equation['math'] = re.sub(r'\\mathrm\s+([A-Za-z])', r'\\mathrm{\1}', equation['math'])
            equation['math'] = re.sub(r'\\overline\s+([A-Za-z])', r'\\overline{\1}', equation['math'])
        section['sources'] = [
            {'label': Path(source).name.replace('_', ' '), 'path': source}
            if isinstance(source, str) and (ROOT / source).is_file() else source
            for source in section.get('sources', [])]
        section['id'] = ident
        section['title'] = f'Figure {n}. ' + section['title'].removeprefix(f'Figure {n}. ')
        document['sections'][n-1] = section
        sections[ident] = section
        legend_copy(section)
    generated = datetime.now(timezone.utc).isoformat()
    document.update(date='23 October 2026', stage='Size, polarity susceptibility and heading persistence',
        subtitle='Polarity replacement, spatial response maps and noise-driven heading memory')
    if args.figure4_snapshot:
        document['stage'] = 'Review candidate; targeted mode campaign active'
    document['overview'] = [baseline['overview'][0], baseline['overview'][1],
        ' '.join(sections['figure3']['narrative']), ' '.join(sections['figure4']['narrative']),
        baseline['overview'][4],
        'This is a working manuscript outline. Figures 1–2 retain their numerical evidence with separate panels for distinct quantities. '
        'The published migration evidence in Figure 5, reserved Figure 6 and prospective pilot appendix retain their preceding content. '
        'Activity polarity, modelled replacement, observed migration and conditional OU velocity memory remain distinct readouts.']
    added = []
    for n in range(1, 5):
        section = sections[f'figure{n}']
        added.append({'title': f'Figure {n}: round6 complete supporting evidence',
            'text': ['Current source data, numerical qualification and preceding supplementary evidence; earlier main artwork is preserved.'],
            'path': section['supplement_pdf'],
            'sources': [{'label': 'Scientific results and rebuild instructions', 'path': section['reader']}]})
    document['supplementary'] = added + document['supplementary']
    path = ROOT / 'manuscript_outline/outline_v12.json'
    write(path, document)
    writer = PdfWriter()
    for ident in sections:
        writer.append(DEST / f'{ident}_with_legend.pdf')
    writer.write(DEST / 'figures1_4_with_legends.pdf')
    protected = {}
    for ident in ('figure5', 'figure6'):
        old = next(s for s in baseline['sections'] if s['id'] == ident)
        new = next(s for s in document['sections'] if s['id'] == ident)
        assert old == new
        protected[ident] = {'section_sha256': hashlib.sha256(json.dumps(old, sort_keys=True).encode()).hexdigest(),
            'asset_sha256': {a['path']: sha(ROOT / a['path']) for a in old['assets']}}
    oldpilot = next(s for s in baseline['supplementary'] if s.get('id') == 'pilot_observables')
    assert oldpilot == next(s for s in document['supplementary'] if s.get('id') == 'pilot_observables')
    write(DEST / 'integration.json', {'generated_utc': generated,
        'baseline_pdf': 'manuscript_outline/outline_v11.pdf',
        'baseline_pdf_sha256': sha(ROOT / 'manuscript_outline/outline_v11.pdf'),
        'baseline_json_sha256': sha(baseline_path), 'outline_json_sha256': sha(path),
        'protected_sections': protected, 'pilot_appendix_unchanged': True,
        'figures1_2_data_unchanged': True, 'round5_files_unchanged': True,
        'review_candidate': bool(args.figure4_snapshot),
        'immutable_figure4_snapshot': relative(args.figure4_snapshot.resolve()) if args.figure4_snapshot else None,
        'archived_request': relative(DEST / 'request/starting_checkpoint.json')})
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Polarity manuscript: round6</title>
<style>body{{max-width:1200px;margin:2rem auto;padding:0 1rem;font:17px/1.55 system-ui;color:#17212b}}img{{max-width:100%;height:auto}}a{{color:#175d98}}.status{{padding:1rem;background:#eef3f7}}.equation{{display:flex;align-items:center;justify-content:center;position:relative;margin:1.5rem 0}}.equation img{{max-width:85%;width:auto}}.equation span{{position:absolute;right:0}}</style>
<h1>Polarity susceptibility and heading persistence</h1><p class="status">Generated {generated}. Document producer finished. Simulation and analysis status are recorded with each figure; this page updates when rebuilt and has no timed refresh.</p>
<p><a href="../outline_v12.pdf">Integrated manuscript outline</a> · <a href="figures1_4_with_legends.pdf">Figures 1–4 with legends</a> · <a href="request/revision_brief_round6_20261023.md">Archived brief</a> · <a href="request/user_prompt_20261023.md">Task prompt</a></p>'''
    if args.figure4_snapshot:
        live = os.path.relpath(FIGURES / 'figure4/index.html', DEST)
        body += f'<p>This is an immutable review candidate with the targeted mode campaign still active. <a href="{live}">Latest Figure 4 producer output and status</a>.</p>'
    for ident, section in sections.items():
        body += '<h2>' + html.escape(section['title']) + '</h2>'
        body += f'<p><a href="{ident}_with_legend.pdf">Figure and legend</a> · <a href="{os.path.relpath(ROOT / section["supplement_pdf"], DEST)}">Complete supplement</a> · <a href="{os.path.relpath(ROOT / section["reader"], DEST)}">Results and methods</a></p>'
        for asset in section['assets']:
            body += f'<img src="{os.path.relpath(ROOT / asset["path"], DEST)}" alt="{html.escape(section["title"], quote=True)}">'
        body += ''.join('<p>' + html.escape(p) + '</p>' for p in section['legend'])
        body += html_equations(ident, section.get('equations', []))
    body += '<p><a href="integration.json">Integration receipt</a></p></html>'
    (DEST / 'index.html').write_text(body)
    (DEST / 'README.md').write_text('''# Round6 manuscript revision

[Integrated outline_v12](../outline_v12.pdf) · [Figures 1–4 and legends](figures1_4_with_legends.pdf) · [Reader](index.html).

The starting point is the current outline_v11 PDF. Its exact hash, the editable
outline, the verbatim brief, the supplied narrative and the task prompt are
archived in [request/](request/). Round5 remains intact. Figures 1–2 change layout
and lettering only; Figures 3–4 report the requested new analyses. Figure 5,
reserved Figure 6 and the prospective pilot appendix are preserved.

Each figure folder exposes one rebuild entry point and a `scientific_readme.md`
beginning with Unexpected. Main figures, compact plotted values and complete
supplements live beside round5 in `revision_round6_20261023/figure1/` through
`figure4/`. Scientific qualifications and campaign completion are recorded there.

After those figure rebuilds, integrate using the existing environment:

```sh
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B manuscript_outline/prepare_round6_20261023.py
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v12.json --output manuscript_outline/outline_v12.pdf
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v12.json --report manuscript_outline/outline_v12.validation.json
```
''')
    update_navigation(sections, bool(args.figure4_snapshot))
    print(path)


if __name__ == '__main__':
    main()
