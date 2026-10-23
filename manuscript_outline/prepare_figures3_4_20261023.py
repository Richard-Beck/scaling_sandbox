"""Integrate the edge-response / domain-birth revision from outline_v10.

Only Figures 3 and 4 change. Prospective pilot observables are a separate
appendix, so the reserved experimental section remains an exact baseline copy.
"""
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
import sys
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import math_to_image
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_figures3_4_20261023'
FIGURES = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round5_20261023'
SOURCE_RELOCATIONS = {'MANUSCRIPT_VISION.md': 'author/vision.md'}
sys.path.insert(0, str(ROOT / 'scripts'))
import build_manuscript_outline as prose


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def relocate_sources(value):
    """Repair a citation to the author document moved before this revision."""
    if isinstance(value, dict):
        return {key: SOURCE_RELOCATIONS.get(item, item) if key == 'path' and isinstance(item, str)
                else relocate_sources(item) for key, item in value.items()}
    if isinstance(value, list):
        return [relocate_sources(item) for item in value]
    return value


def html_equations(ident, equations):
    """Render display mathematics offline; no raw TeX or remote renderer."""
    if not equations:
        return ''
    folder = DEST / 'equations'
    folder.mkdir(exist_ok=True)
    result = ''
    for index, equation in enumerate(equations):
        path = folder / f'{ident}_{index + 1}.svg'
        math_to_image('$' + equation['math'] + '$', path, dpi=180, format='svg',
                      prop=FontProperties(size=14))
        number = html.escape(str(equation.get('number', '')))
        result += (f'<div class="equation"><img src="equations/{path.name}" alt="Model equation {number}">'
                   f'<span>({number})</span></div><p>' + html.escape(equation.get('description', '')) + '</p>')
    return result


def pilot_sheet(pilots):
    """A readable native PDF table, included verbatim in the manuscript PDF."""
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(11.69, 8.27))
    ax = fig.add_axes([.035, .30, .93, .54])
    ax.axis('off')
    widths = [.20, .19, .22, .25, .14]
    headings = ['Cellular observable', 'Model counterpart', 'Drift regime:\nWGD+ / WGD−',
                'Secondary-domain regime:\nWGD+ / WGD−', 'Mechanistic\ninterpretation']
    keys = ['observable', 'counterpart', 'drift_direction', 'secondary_direction', 'status']
    chars = [29, 27, 32, 36, 20]
    rows = [['\n'.join(textwrap.wrap(str(row[key]), width=width,
                                  break_long_words=False, break_on_hyphens=False))
             for key, width in zip(keys, chars)] for row in pilots['observables']]
    table = ax.table(cellText=rows,
                     colLabels=headings,
                     colWidths=widths, cellLoc='left', loc='center', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor('#bcc6cd')
        cell.set_linewidth(.45)
        if row == 0:
            cell.set_text_props(weight='bold')
            cell.set_facecolor('#eaf0f3')
    fig.text(.035, .945, 'Prospective SUM-159 WGD pilot observables', fontsize=16, weight='bold')
    fig.text(.035, .90, pilots['fixed_before_data'], fontsize=10)
    fig.text(.035, .86, '\n'.join(textwrap.wrap(pilots['regime_definition'], 155)), fontsize=9)
    y = .255
    for paragraph in pilots['interpretation_notes']:
        lines = textwrap.wrap(paragraph, 175, break_long_words=False)
        fig.text(.035, y, '\n'.join(lines), fontsize=8.3, va='top', linespacing=1.3)
        y -= .019 * len(lines) + .010
    assert y > .025, 'Pilot notes exceed page'
    fig.savefig(DEST / 'pilot_observables.pdf')
    fig.savefig(DEST / 'pilot_observables.png', dpi=180)
    plt.close(fig)


def publication_copy(section):
    """Use the repository's measured text layout for a flowing legend PDF."""
    ident = section['id']
    legend_path = DEST / f'{ident}_legend.pdf'
    document = {'title': section['title'], 'sections': []}
    source_path = FIGURES / ident / 'legend.json'
    builder = prose.OutlineBuilder(document, ROOT, legend_path, source_path)
    with PdfPages(legend_path) as pdf:
        builder.pdf = pdf
        builder.section_id = ident
        builder.new_page(section['title'], kind='figure_notes')
        for paragraph in section['legend']:
            builder.paragraph(paragraph, size=10.2)
        for equation in section.get('equations', []):
            builder.equation(equation)
        builder.finish_page()
    writer = PdfWriter()
    writer.append(ROOT / section['main_pdf'])
    writer.append(legend_path)
    writer.write(DEST / f'{ident}_with_legend.pdf')


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    baseline_path = ROOT / 'manuscript_outline/outline_v10.json'
    baseline = read(baseline_path)
    document = relocate_sources(deepcopy(baseline))
    sections = {}
    for ident in ('figure3', 'figure4'):
        source = read(FIGURES / ident / 'legend.json')
        assert source['id'] == ident
        section = next(s for s in document['sections'] if s['id'] == ident)
        section.clear()
        section.update(source)
        section['title'] = f'Figure {ident[-1]}. ' + source['title'].removeprefix(f'Figure {ident[-1]}. ')
        sections[ident] = section
        publication_copy(section)
    pilots = read(DEST / 'pilot_observables.json')
    pilot_sheet(pilots)
    generated = datetime.now(timezone.utc).isoformat()
    document.update(date=datetime.now(timezone.utc).strftime('%d %B %Y'),
                    stage='Size, polarity response and heading persistence',
                    subtitle='Size-dependent polarity responses and prospective cellular measurements')
    # The new descriptions come from the delivered sections, not the old null narrative.
    figure2_summary = baseline['overview'][1].split('Cue placement relative')[0].strip()
    document['overview'] = [baseline['overview'][0], figure2_summary,
                           ' '.join(sections['figure3']['narrative']),
                           ' '.join(sections['figure4']['narrative']),
                           'The prospective pilot table fixes model-linked cellular observables before the SUM-159 data arrive. '
                           'Directions are conditional hypotheses; unavailable model measurements remain unavailable. '
                           'Activity heading, observed migration summaries and conditional OU velocity-memory estimates remain separate evidence.',
                           baseline['overview'][-1]]
    pilot_section = dict(
        id='pilot_observables', title='Prospective pilot observables before SUM-159 data',
        role=pilots['fixed_before_data'],
        assets=[dict(path=relative(DEST / 'pilot_observables.png'),
                     native_pdf_page=relative(DEST / 'pilot_observables.pdf'),
                     page_size=[11.69, 8.27], minimum_printed_width_inches=11.69, page_title='')],
        legend=pilots['interpretation_notes'],
        narrative=[pilots['regime_definition']],
        measurements=[r['operational_definition'] for r in pilots['observables']],
        limitations=pilots['limitations'],
        sources=[dict(label='Prospective directions and operational definitions',
                      path=relative(DEST / 'pilot_observables.json'))])
    additions = [pilot_section]
    for ident in ('figure3', 'figure4'):
        section = sections[ident]
        additions.append(dict(title=f'Figure {ident[-1]}: edge-response / secondary-domain revision supplement',
                              text=['New measurements, identification limits, numerical checks and preceding supplementary evidence.'],
                              path=section['supplement_pdf'], sources=section['sources']))
    document['supplementary'] = additions + document['supplementary']
    outline = ROOT / 'manuscript_outline/outline_v11.json'
    write(outline, document)
    protected = {}
    for ident in ('figure1', 'figure2', 'figure5', 'figure6'):
        old = next(s for s in baseline['sections'] if s['id'] == ident)
        new = next(s for s in document['sections'] if s['id'] == ident)
        assert relocate_sources(old) == new
        protected[ident] = dict(section_sha256=hashlib.sha256(json.dumps(old, sort_keys=True).encode()).hexdigest(),
                               asset_sha256={a['path']: sha(ROOT / a['path']) for a in old['assets']})
    writer = PdfWriter()
    for ident in ('figure2', 'figure3', 'figure4'):
        writer.append(ROOT / next(s for s in document['sections'] if s['id'] == ident)['main_pdf'])
    writer.write(DEST / 'figures_2_4.pdf')
    integrated = PdfWriter()
    for ident in ('figure3', 'figure4'):
        integrated.append(DEST / f'{ident}_with_legend.pdf')
    integrated.append(DEST / 'pilot_observables.pdf')
    integrated.write(DEST / 'figures3_4_with_legends_and_pilots.pdf')
    receipt = dict(generated_utc=generated, baseline=relative(baseline_path),
                   baseline_sha256=sha(baseline_path), outline=relative(outline), outline_sha256=sha(outline),
                   starting_checkpoint=read(DEST / 'request/starting_checkpoint.json'),
                   protected_sections=protected, all_other_scientific_content_and_artwork_identical=True,
                   source_relocations=SOURCE_RELOCATIONS,
                   reserved_experimental_content_unchanged=True, no_new_noise_campaign=True,
                   pilot_table_sha256=sha(DEST / 'pilot_observables.json'),
                   publication_pdf_sha256={ident: sha(DEST / f'{ident}_with_legend.pdf') for ident in sections})
    write(DEST / 'integration.json', receipt)
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Figures 3–4 revision</title>
<style>body{{max-width:1150px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}}img{{width:100%}}a{{color:#166aa2}}.equation{{display:flex;align-items:center;justify-content:center;gap:2rem;margin:1.5rem 0}}.equation img{{width:auto;max-width:85%;height:auto}}.equation span{{margin-left:auto}}</style>
<h1>Edge responses, drift and secondary-domain dynamics</h1><p>Generated {generated}; document producer finished. Simulation and analysis completion are recorded with each figure.</p>
<p><a href="../outline_v11.pdf">Integrated manuscript outline</a> · <a href="figures3_4_with_legends_and_pilots.pdf">Figures, legends and pilots</a> · <a href="pilot_observables.pdf">Prospective pilot table</a></p>
'''
    for ident, section in sections.items():
        body += '<h2>' + html.escape(section['title']) + '</h2>'
        body += f'<p><a href="{ident}_with_legend.pdf">Figure and legend</a> · <a href="{os.path.relpath(ROOT / section["supplement_pdf"], DEST)}">Full supplement</a></p>'
        for asset in section['assets']:
            body += f'<img src="{os.path.relpath(ROOT / asset["path"], DEST)}" alt="{html.escape(section["title"], quote=True)}">'
        body += ''.join('<p>' + html.escape(p) + '</p>' for p in section['legend'])
        body += html_equations(ident, section.get('equations', []))
    body += '<p><a href="request/revision_brief_20261023.md">Archived revision brief</a> · <a href="request/user_prompt_20261023.md">Original task prompt</a> · <a href="revision_acceptance.md">Delivery and identification limits</a> · <a href="claude_response.md">Claude reviews and point-by-point responses</a> · <a href="integration.json">Source and protected-section receipt</a></p></html>'
    (DEST / 'index.html').write_text(body)
    readme = '''# Figures 3–4: edge response, drift and secondary-domain dynamics

[Integrated outline_v11](../outline_v11.pdf) · [Figures, legends and pilot table](figures3_4_with_legends_and_pilots.pdf) · [Reader](index.html).

Starting point: the committed outline_v10 PDF and corrected standalone Figure 4, recorded by exact hashes in [the starting checkpoint](request/starting_checkpoint.json). The [revision brief](request/revision_brief_20261023.md) and [task prompt](request/user_prompt_20261023.md) are preserved verbatim.

[Delivery and identification limits](revision_acceptance.md) · [Point-by-point Claude responses](claude_response.md) · [Verbatim review hashes](review_preservation.json). Claude reviewed the starting checkpoint and revision checkpoint `49754ab`; the final deliverables incorporate responses to both reviews. The post-review changes have not received another Claude review.

Figures 1, 2, 5 and the reserved Figure 6 preserve their outline_v10 scientific content and artwork. The Figure 6 and bibliography vision citations are repaired from the deleted root MANUSCRIPT_VISION.md to the existing author/vision.md; no author text changes. No new noisy histories are generated. New deterministic edge measurements and derived compact-record analyses are documented with each figure; missing measures and failed assumptions remain explicit. Prospective observables are in [the pilot table](pilot_observables.json).

Rebuild using the existing scaling_sandbox environment after both figure producers finish:

```sh
python -B manuscript_outline/prepare_figures3_4_20261023.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v11.json --output manuscript_outline/outline_v11.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v11.json --report manuscript_outline/outline_v11.validation.json
python -B manuscript_outline/verify_figures3_4_20261023.py
```
'''
    (DEST / 'README.md').write_text(readme)
    marker = '<!-- Figures 3-4 edge-birth revision 20261023 -->'
    for path, prefix in [(ROOT / 'manuscript_outline/README.md', ''),
                         (ROOT / 'current manuscript figs/README.md', '../manuscript_outline/')]:
        old = path.read_text()
        if old.startswith(marker):
            old = old.split('<!-- Earlier editions preserved below -->\n', 1)[1]
        summary = ' '.join(sections['figure3']['narrative'] + sections['figure4']['narrative'])
        links = (f'[**Integrated outline_v11**]({prefix}outline_v11.pdf) · '
                 f'[**Figures, legends and pilots**]({prefix}revision_figures3_4_20261023/figures3_4_with_legends_and_pilots.pdf) · '
                 f'[**Reader**]({prefix}revision_figures3_4_20261023/index.html)')
        path.write_text(marker + '\n# Figures 3–4: edge responses and secondary-domain dynamics\n\n' +
                        links + '\n\n' + summary + '\n\nFigures 1, 2, 5 and the reserved Figure 6 retain their preceding sections and artwork. '
                        'Prospective pilot directions remain conditional on the available model measurements.\n\n'
                        '<!-- Earlier editions preserved below -->\n' + old)
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
