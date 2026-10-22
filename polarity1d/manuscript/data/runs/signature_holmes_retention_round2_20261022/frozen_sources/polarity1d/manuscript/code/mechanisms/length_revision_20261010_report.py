"""Build the 10 October length-mechanism reader from audited vector figures.

Rendering only: native checkpoints, trajectories and frozen sources are read-only.
Run from the repository with the existing scaling_sandbox Python environment.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import html
import io
import json
import re
import shutil
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath, TextToPath
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'figures/length_mechanisms_20261008'
REV = FIG / 'revision_20261010'
DATA = ROOT / 'data/runs/length_mechanisms_20261008/revision_20261010'
TITLE = 'Mechanisms of cell-length-dependent polarity reversal'
TEXT_METRICS = TextToPath()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def relative(path):
    return str(Path(path).relative_to(FIG))


def asset(stem, group):
    return relative(REV / group / (stem + '.png'))


def caption(image):
    return (FIG / image).with_name(Path(image).stem + '_caption.txt').read_text().strip()


def chapters():
    baseline = read(REV / 'chemical_chapters_baseline.json')
    editorial = read(REV / 'editorial/chapter_overrides.json')
    for c in baseline[1:]:
        c.update(editorial[c['id']])
    c = baseline[0]
    c.update(title='Larger cells carry a greater chemical reorganization burden',
             image=asset('length_responses', 'figure1'),
             paragraphs=[
        'Length changes the shortest cue that leaves a strong reversed pole after withdrawal (Fig. 1A). Each model column puts the three cue geometries on the same axes. A fixed-width patch becomes slower in Wave-pinning and Holmes, and Goryachev slows under all three cues. Relative patches and constant-slope ramps expose different length responses, including decreasing intervals. The chemical and mechanical chapters explain why these trends differ.',
        'The starting-state burden has a common direction: larger cells require more spatial chemical reorganization (Fig. 1B). For each chemical pool, we combine active and inactive forms before comparing its initial distribution with its reflected state. We then add the pool-specific Wasserstein distances. This total increases at every adjacent eligible length in all five models. Holmes includes both total Rac and total Rho; deBelly includes Rac, Rho and bound-plus-unbound MCA. Conversion between forms can change activity locally, but it does not erase this spatial asymmetry of the total pools.',
        'The input itself also matters. A ramp at constant physical slope becomes stronger as the cell lengthens. Matched-peak controls remove reversal at several doses that pass with the growing ramp, and disabling Spring extension feedback removes all four tested reversals. The report therefore follows two connected questions: what must be reorganized, and what sets the rate at which the new pole can replace the old one.'
    ])
    c['caption'] = caption(c['image'])
    c['sources'] = [('revision_20261010/figure1/chemical_reorganization_length_summary.csv', 'Chemical-distance measurements'),
                    ('../../data/runs/three_test_mechanisms_20261006/results/cell_evidence.json', 'Matched-peak and Spring feedback controls')]
    redraft = read(REV / 'debelly/chapter_assets.json')
    narrative = [
        ('overview', '8', 'Rac and Rho turn before MCA memory is rewritten', [
            'The reversal kymographs establish the behavior to explain. At the top of Figure 8, the original Rac-rich front is on the left and the challenged end is on the right. Rac gains at the new end within seconds (a), and Rho shifts toward the future rear (b). At the common 60-s line, these chemical polarities have turned in all three cells, while bound MCA still favors the original rear (c). MCA is the slow memory variable in this reversal sequence.',
            'The mechanical state responds quickly at 23 and 31 µm, but the 39-µm cell retains its original mean cortical-flow direction until about 485 s (d) and its tension bias until about 530 s (e). This long-cell exception is visible in the fields rather than hidden in an average. Bound MCA takes hundreds of seconds to redistribute and can reverse direction again during the pulse and recovery.',
            'Cue withdrawal preserves the inherited chemical and mechanical state. Rac and tension relax, and MCA continues moving. All three displayed pulses ultimately pass the 1% scalar strength gate after basal recovery. The immediate chemical turn is therefore the opening response; the slower task is rebuilding the attachment distribution that supports the replacement pole.'
        ]),
        ('mechanics', '9', 'Cortical motion sweeps MCA out and redistributes it', [
            'Figure 9 follows one cell through cue onset and withdrawal. At rest, membrane springs, cortical viscosity and membrane–cortex attachment (MCA) linkers maintain a polarized mechanical state (A). Rac sets the endpoint protrusion commands, Rho supplies distributed cortical stress, and bound MCA transmits drag between cortex and membrane. The flow components (B) and MCA fluxes (C) show the resulting motion and transport.',
            'Thirty seconds after cue onset, the right endpoint begins driving flow in the opposite direction while the old left command remains active (Fig. 9D–F). Bound MCA is carried by the cortex; unbound MCA follows the membrane and diffuses much faster. Outward advection sweeps MCA through the domain boundaries. The model replaces the lost material uniformly across the cell, as the small distributed arrows indicate. MCA is continually removed and redistributed, rather than accumulating at the downstream boundary.',
            'At withdrawal, the fields remain continuous and retain strong leftward interior flow (Fig. 9G–I). Once old-front Rac falls through its protrusion switch, the left command weakens and the flow turns rightward (J–L). The local mass balances show transport opposing uniform refill throughout. Mechanics sets the transport speed, while sweep-out, replacement and binding set the net redistribution. None of these processes automatically accelerates in proportion to cell length.'
        ]),
        ('distance', '10', 'MCA must move farther, while average speed falls behind', [
            'The distance panels compare the net MCA centroid shift during each cue with two starting-state benchmarks: the distance to the midpoint and the distance to the reflected polarized distribution (Fig. 10A). All three geometries begin from the same checkpoint at each length. From 23 to 39 µm, the reflected-centroid separation grows 1.83-fold. The initial attachment memory therefore spans a greater absolute distance in the larger cell.',
            'The speed panels separate the first 60 s, the final 60 s and the average over the selected cue (Fig. 10B). Over the same length range, average speed changes by factors of 0.90, 0.55 and 1.09 for fixed, relative and ramp cues. None keeps pace with the growing distance. The gray proportional-to-length references make that deficit visible.',
            'The exceptions sharpen the next question. Intermediate lengths are nonmonotonic, and the ramp’s early speed rises almost in proportion to length. Early progress, late progress and whole-cue averages tell different stories. We first explain why the cells start with different MCA distributions, then return to the balance that controls their reversal.'
        ]),
        ('starting_mca', '11', 'The MCA starting state emerges across two polarity transitions', [
            'The steady MCA profiles explain the growing initial distance in Figure 10 (Fig. 11A). Short cells remain nearly uniform. The sampled asymmetric continuation loses polarity between 17.9 and 17.925 µm, while the homogeneous state becomes unstable near 19.53 µm (Fig. 11B). The interval between these boundaries is hard to polarize: a cue-established pole can persist, but a small perturbation of homogeneity decays. Above the second boundary, small asymmetry grows without a cue.',
            'The flux comparison links this change to MCA transport (Fig. 11C). Advection strengthens relative to the opposing combined bound-plus-unbound diffusion as the asymmetric state develops: cortical flow carries bound MCA while diffusion returns material down its concentration gradient. The branch and stability boundaries come from their separate continuation and growth tests; equality of flux magnitudes is a transport diagnostic.',
            'The remaining question is what produces the stronger cortical drive. Figure 12 connects it to the Rac protrusion switch and the way membrane tension feeds back onto chemistry.'
        ]),
        ('steady_mechanism', '12', 'The protrusion switch converts small Rac differences into strong cortical drive', [
            'The MCA profiles (Fig. 12A) and their spatial fluxes (B) establish the state maintained by cortical motion. At 19 µm, opposing advection and diffusion are comparable. At larger lengths, bound-MCA advection dominates the flux magnitude. The flow decomposition identifies endpoint forcing as the growing drive, with a smaller opposing Rho contribution (C).',
            'The Rac and Rho regions occupy similar fractions of each cell, so their absolute widths expand with length (Fig. 12D). Their endpoint concentrations change much less. Those small Rac differences matter because the old front approaches the steep protrusion-switch midpoint (E): the native outward command rises strongly as length increases (F). Larger cells therefore sustain much stronger cortical motion at the front.',
            'That drive loads the membrane. Tension stimulates Rho, which opposes Rac, and tension also directly attenuates protrusion. These feedbacks limit the drive. Longer cells distribute deformation over more membrane: front tension turns down even while the command rises, whereas rear tension grows sharply (Fig. 12G,H). An equal-command control confirms the geometric effect: at 23 and 39 µm, the same endpoint command and displacement give front tensions of 2.762 and 1.214 pN.',
            'The steady-state pattern is clear. A larger MCA distribution is maintained by stronger cortical drive. The longer membrane accommodates that drive with less front tension, but carries much greater tension at the unpolarized end. Reversal inherits both this attachment memory and the mechanics that maintain it.'
        ]),
    ]
    result = baseline
    for key, label, title, paragraphs in narrative:
        a = redraft[key]
        image = a['image']
        result.append(dict(id=key, figure_label=label, title=title,
                           image=image, paragraphs=paragraphs,
                           caption=a.get('caption') or caption(image),
                           sources=a.get('sources', [])))
        if key == 'steady_mechanism':
            result[-1]['sources'].extend([
                ('narrative_20261009/mechanics/debelly_resting_equal_commands.pdf', 'Equal-command mechanical control'),
                ('../../data/runs/length_mechanisms_20261008/narrative_20261009/mechanics/tension_command_measurements.csv', 'Control tension measurements')])
    ending = read(REV / 'conclusion/chapter.json')
    result.append(dict(ending, id='clearance', figure_label='13'))
    for n, c in enumerate(result, 1):
        c['figure_label'] = str(n)
        c['manifest'] = c.get('manifest') or str(Path(c['image']).with_suffix('.manifest.json'))
        if not c.get('caption'):
            c['caption'] = caption(c['image'])
        c['caption'] = re.sub(r'^\*?Figure\s+[^ .\n*]+\*?[. ]*', f'Figure {n} ', c['caption'], count=1)
        c.setdefault('sources', [])
    return result


def resolve_artifact(manifest, record, name):
    p = Path(name)
    if p.is_absolute():
        return p
    choices = [manifest.parent / record.get('artifact_root', '.') / p,
               FIG / p, ROOT / p]
    return next((q for q in choices if q.is_file()), choices[0])


def verify(chapters):
    status = read(ROOT / 'data/runs/debelly_corrected_three_test_20261008/results/status.json')
    assert status['complete'] and status['producer'] == 'finished' and status['execution_errors'] == 0
    assert read(ROOT / 'data/runs/debelly_corrected_three_test_20261008/results/independent_audit.json')['passed']
    sources, artifacts, manifests = {}, {}, []
    for c in chapters:
        image = FIG / c['image']
        manifest = FIG / c['manifest']
        record = read(manifest)
        assert record.get('producer', 'finished') == 'finished'
        assert record.get('passed', record.get('audit_passed', True))
        hashes = record.get('artifacts', record.get('sha256', record.get('figure_sha256', {})))
        assert hashes, manifest
        for name, digest in hashes.items():
            p = resolve_artifact(manifest, record, name)
            assert p.is_file() and sha(p) == digest, (manifest, p)
            artifacts[str(p.relative_to(ROOT))] = digest
        source_hashes = record.get('sources', record.get('source_sha256', record.get('sources_sha256', {})))
        for name, digest in source_hashes.items():
            p = Path(name)
            if not p.is_absolute():
                p = ROOT / p
            assert p.is_file() and sha(p) == digest, (manifest, p)
            sources[str(p.relative_to(ROOT))] = digest
        if record.get('renderer'):
            renderer = Path(record['renderer'])
            if not renderer.is_absolute():
                renderer = next((p for p in [ROOT / renderer, ROOT.parents[1] / renderer]
                                 if p.is_file()), ROOT / renderer)
            assert sha(renderer) == record['renderer_sha256'], renderer
            sources[str(renderer.relative_to(ROOT))] = sha(renderer)
        assert str(image.relative_to(ROOT)) in artifacts, ('unhashed image', image)
        vector = image.with_suffix('.pdf')
        assert vector.is_file() and len(PdfReader(vector).pages) == 1
        assert str(vector.relative_to(ROOT)) in artifacts, ('unhashed vector', vector)
        manifests.append(dict(path=str(manifest.relative_to(ROOT)), sha256=sha(manifest),
                              generated_at=record.get('generated_at')))
        for destination, _ in c['sources']:
            if not destination.startswith(('https:', 'http:')):
                path = (FIG / destination.split('#', 1)[0]).resolve()
                assert path.exists(), destination
                if path.is_file():
                    sources[str(path.relative_to(ROOT))] = sha(path)
    first_draft = ROOT / 'figures/archive/length_mechanisms_first_draft_20261008/length_mechanisms.pdf'
    prior = FIG / 'archive_20261010_before_redraft'
    preserved = read(prior / 'preservation.json')
    for name, digest in preserved.items():
        assert sha(prior / Path(name).name) == digest
    return dict(comparison_generated_at=status['generated_at'], sources=sources,
                artifacts=artifacts, manifests=manifests,
                archived_first_draft_sha256=sha(first_draft))


def wrap(text, width, size, weight='normal'):
    """Wrap using measured glyph widths so scientific units cannot overrun a page."""
    font = FontProperties(family='DejaVu Sans', weight=weight, size=size)
    lines, current = [], ''
    for word in text.split():
        candidate = (current + ' ' + word).strip()
        extent = TEXT_METRICS.get_text_width_height_descent(candidate, font, ismath=False)[0]
        if current and extent > width:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def pdf_report(chapters, stamp):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
    writer, layouts = PdfWriter(), []
    width, margin = 842.4, 38.
    usable = width - 2 * margin
    for n, c in enumerate(chapters, 1):
        vector = (FIG / c['image']).with_suffix('.pdf')
        source = PdfReader(vector).pages[0]
        fw, fh = float(source.mediabox.width), float(source.mediabox.height)
        titles = wrap(c['title'], usable, 17., 'bold')
        bodies = [wrap(p, usable, 11.) for p in c['paragraphs']]
        prefix = 'Figure ' + str(n)
        legend = c['caption'][len(prefix):].lstrip(' .')
        # A small first-line indentation reserves room for the italic prefix.
        legends = wrap(prefix + '  ' + legend, usable, 9.5)
        figure_height = usable * fh / fw
        y = 43 + len(titles) * 23 + 14 + sum(len(lines) * 15.6 + 11 for lines in bodies)
        plot_top = y + 9
        plot_bottom = plot_top + figure_height
        height = plot_bottom + 22 + len(legends) * 13.5 + 48
        fig = plt.figure(figsize=(width / 72, height / 72))
        def put(text, top, size=11., color='#233d48', weight='normal', style='normal', x=margin):
            fig.text(x / width, 1 - top / height, text, va='top', fontsize=size,
                     color=color, weight=weight, style=style, linespacing=1.42)
        put(f'POLARITY REVERSAL  ·  FIGURE {n}', 22, 8.5, '#63747c', 'bold')
        put('\n'.join(titles), 43, 17., '#173c49', 'bold')
        pos = 43 + len(titles) * 23 + 14
        for lines in bodies:
            put('\n'.join(lines), pos)
            pos += len(lines) * 15.6 + 11
        legend_y = plot_bottom + 22
        first = legends[0]
        put(prefix, legend_y, 9.5, '#42565f', style='italic')
        prefix_w = TextPath((0, 0), prefix + '  ', size=9.5,
                            prop=FontProperties(family='DejaVu Sans', style='italic')).get_extents().width + 3
        put(first[len(prefix):].lstrip(), legend_y, 9.5, '#42565f', x=margin + prefix_w)
        if len(legends) > 1:
            put('\n'.join(legends[1:]), legend_y + 13.5, 9.5, '#42565f')
        fig.text(margin / width, 17 / height, 'Generated ' + stamp[:19] + ' UTC · native-data producers finished',
                 fontsize=7.5, color='#63747c')
        fig.text((width - margin) / width, 17 / height, str(n), ha='right', fontsize=8, color='#63747c')
        buffer = io.BytesIO()
        fig.savefig(buffer, format='pdf')
        plt.close(fig)
        buffer.seek(0)
        page = PdfReader(buffer).pages[0]
        scale = usable / fw
        page.merge_transformed_page(source, Transformation().scale(scale).translate(
            margin - float(source.mediabox.left) * scale,
            height - plot_bottom - float(source.mediabox.bottom) * scale))
        writer.add_page(page)
        writer.add_outline_item(f'Figure {n}: {c["title"]}', n - 1)
        layouts.append(dict(figure=n, page_width_points=width, page_height_points=height,
                            vector=str(vector.relative_to(ROOT)),
                            plot_bounds=[margin, height - plot_bottom, width - margin, height - plot_top],
                            narrative_lines=sum(map(len, bodies)), legend_lines=len(legends)))
    writer.add_metadata({'/Title': TITLE, '/Author': 'Polarity manuscript working report',
                         '/Subject': 'Chemical reorganization and membrane–cortex attachment memory'})
    with (FIG / 'length_mechanisms.pdf').open('wb') as stream:
        writer.write(stream)
    (DATA / 'report_layout.json').write_text(json.dumps(layouts, indent=2) + '\n')


def html_report(chapters, stamp, evidence):
    esc = html.escape
    nav = ''.join(f'<a href="#{c["id"]}">{n}. {esc(c["title"])}</a>' for n, c in enumerate(chapters, 1))
    sections = []
    for n, c in enumerate(chapters, 1):
        paragraphs = ''.join('<p>' + esc(p) + '</p>' for p in c['paragraphs'])
        prefix = 'Figure ' + str(n)
        legend = '<em>' + prefix + '</em> ' + esc(c['caption'][len(prefix):].lstrip(' .'))
        links = [(str(Path(c['image']).with_suffix('.pdf')), 'Vector figure'),
                 (c['manifest'], 'Data, methods and hashes')]
        links.extend((url, label) for url, label in c['sources'] if url not in {p[0] for p in links})
        method_path = (FIG / c['image']).with_name(Path(c['image']).stem + '_methods.txt')
        if method_path.exists():
            links.append((relative(method_path), 'Technical notes'))
        linked = ' · '.join(f'<a href="{esc(url, quote=True)}">{esc(label)}</a>' for url, label in links)
        sections.append(f'<section id="{c["id"]}"><div class="number">Figure {n}</div><h2>{esc(c["title"])}</h2>'
                        f'<div class="prose">{paragraphs}</div><figure><a href="{esc(c["image"], quote=True)}">'
                        f'<img src="{esc(c["image"], quote=True)}" alt="{esc(c["caption"], quote=True)}" loading="{"eager" if n == 1 else "lazy"}"></a>'
                        f'<figcaption>{legend}</figcaption></figure><p class="sources">{linked}</p></section>')
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{TITLE}</title><style>
:root{{--ink:#173c49;--muted:#586e78;--accent:#197987;--rule:#dce7e9;--wash:#f3f7f7}}
*{{box-sizing:border-box}}body{{margin:0;color:var(--ink);background:white;font:17px/1.65 system-ui,sans-serif}}
main{{max-width:1270px;padding:46px 38px 64px;margin:auto}}header{{border-bottom:2px solid var(--ink);padding-bottom:22px}}
h1{{font-size:clamp(29px,4vw,41px);line-height:1.15;max-width:1000px}}h2{{font-size:clamp(24px,3vw,32px);line-height:1.23;margin:9px 0 20px;max-width:1050px}}
.eyebrow,.number{{font-size:12px;letter-spacing:.06em;font-weight:650;text-transform:uppercase;color:var(--muted)}}
p{{margin:0 0 16px}}.intro{{max-width:910px;font-size:19px}}.metadata,.sources{{font-size:13px;color:var(--muted)}}
a{{color:var(--accent);text-underline-offset:3px}}.download{{display:inline-block;padding:9px 17px;background:var(--ink);color:white;border-radius:4px;text-decoration:none;margin:6px 14px 16px 0}}
nav{{display:grid;grid-template-columns:1fr 1fr;gap:7px 28px;padding:21px 25px;margin:26px 0;background:var(--wash);font-size:14px}}nav a{{text-decoration:none}}
section{{padding:36px 0 30px;border-bottom:1px solid var(--rule)}}.prose{{max-width:980px}}figure{{margin:25px 0 15px}}img{{display:block;width:100%;height:auto}}
figcaption{{font-size:14px;line-height:1.55;color:#425963;max-width:1110px;margin-top:15px}}details{{padding:22px 25px;margin-top:30px;background:var(--wash);font-size:14px}}
summary{{cursor:pointer;font-weight:650}}details p{{margin-top:13px}}footer{{font-size:13px;color:var(--muted);margin-top:22px}}
@media(max-width:700px){{main{{padding:25px 18px}}nav{{grid-template-columns:1fr}}body{{font-size:16px}}}}
@media print{{nav,.download,details{{display:none}}main{{padding:0}}section{{break-before:page}}}}
</style></head><body><main><header><div class="eyebrow">Polarity manuscript · revision 10 October 2026</div>
<h1>{TITLE}</h1><p class="intro">Larger cells carry more spatial chemical asymmetry to reorganize. The rate of replacement depends on local reaction kinetics, domain growth and, in deBelly, the mechanical transport of membrane–cortex attachments.</p>
<a class="download" href="length_mechanisms.pdf">Report PDF · {len(chapters)} vector pages</a><a href="report_text.txt">Plain-text reader</a>
<p class="metadata">Report generated {esc(stamp)}. Core comparison generated {esc(evidence['comparison_generated_at'])}. Native-data producers finished; each figure links its source generation time and audit. Refreshing this page loads the latest generated report.</p></header>
<nav aria-label="Figures">{nav}</nav>{''.join(sections)}
<details><summary>Measurements and reproduction</summary>
<p>Figure 1 and the deBelly reversals use the shortest tested cue that passes the canonical 1% scalar strength gate, followed by 10000 native time units of basal recovery. The gate requires persistent reversed orientation and contrast in the final 500 units, with first-moment magnitude within 1% of the original. The separate native-field tail flag remains in the data and is shown by open markers in Figure 1. Figure 6 compares an unsuccessful and a successful example pulse; Figure 7 retains uninterrupted forcing.</p>
<p>The chemical distance in Figure 1 is the sum of independently mass-normalized Wasserstein distances for the total chemical pools and their reflected distributions. It measures spatial reorganization distance; amount-weighted costs are exported separately. Spring uses native resting length on the length axis and its actual basal physical coordinates for distance. deBelly centroid distances use nominal reference coordinates. Centroid speed describes a changing distribution, including binding and replacement, rather than the travel of tracked molecules.</p>
<p>The deBelly evidence uses corrected independent endpoint commands. The mechanical and transport decompositions reconstruct the actual coupled trajectories. They are source allocations, while the retained feedback-ablation experiments independently evolve changed dynamics. Early and late Figure 10 windows span 60 s; late and whole-cue windows vary with the selected duration. Figure 13 also uses a common 60–180-s comparison. These distinct windows are identified in the legends and data.</p>
<p><a href="../../THREE_TESTS.md">Protocol</a> · <a href="../../ASSUMPTIONS.md">Model assumptions</a> · <a href="revision_20261010/report_provenance.json">Report provenance and audit</a> · <a href="archive_20261010_before_redraft/reader.html">Preserved previous reader</a></p>
<p>Rebuild with <code>/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_revision_20261010_report.py</code>. Figure renderers and source hashes are recorded in their attached manifests. The earlier full-field, recovery and feedback controls remain accessible from the preserved reader.</p>
</details><footer>Chemical starting state → local resistance and domain growth → MCA transport and mechanical memory.</footer></main></body></html>'''
    (FIG / 'index.html').write_text(document)


def audit(chapters):
    reader = PdfReader(FIG / 'length_mechanisms.pdf')
    assert len(reader.pages) == len(chapters) == 13
    for n, (page, c) in enumerate(zip(reader.pages, chapters), 1):
        text = page.extract_text()
        assert f'FIGURE {n}' in text and f'Figure {n}' in text
        assert len(text) > 600 and '/Font' in page['/Resources']
        assert re.match(r'^Figure ' + str(n) + r' ', c['caption'])
        assert not re.search(r'Fig\. 1[A-Z][–-]', ' '.join(c['paragraphs']))
    document = (FIG / 'index.html').read_text()
    assert 'http-equiv="refresh"' not in document.lower()
    destinations = re.findall(r'(?:href|src)="([^"#]+)"', document)
    for destination in destinations:
        if not destination.startswith(('https:', 'http:', 'mailto:')):
            assert (FIG / destination.split('#', 1)[0]).exists(), destination
    return dict(passed=True, pages=13, vector_figures=True, searchable_text=True,
                consistent_italic_legend_prefix=True, local_links_checked=len(destinations),
                no_timed_refresh=True)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat()
    content = chapters()
    evidence = verify(content)
    prior = FIG / 'archive_20261010_before_redraft'
    archived_html = (prior / 'index.html').read_text()
    archived_html = archived_html.replace('<head>', '<head><base href="../">', 1)
    for name in ['length_mechanisms.pdf', 'report_text.txt']:
        archived_html = archived_html.replace('href="' + name + '"',
            'href="archive_20261010_before_redraft/' + name + '"')
    archived_html = archived_html.replace(
        '../../data/runs/length_mechanisms_20261008/report_provenance.json',
        'archive_20261010_before_redraft/report_provenance.json')
    (prior / 'reader.html').write_text(archived_html)
    (REV / 'chapters.json').write_text(json.dumps(content, indent=2, ensure_ascii=False) + '\n')
    # Create the provenance destination before HTML link verification.
    (REV / 'report_provenance.json').write_text('{}\n')
    pdf_report(content, stamp)
    html_report(content, stamp, evidence)
    plain = [TITLE, 'Generated ' + stamp + ' · native-data producers finished', '']
    for c in content:
        plain.extend(['Figure ' + c['figure_label'] + ' — ' + c['title'], '',
                      '\n\n'.join(c['paragraphs']), '', c['caption'], ''])
    (FIG / 'report_text.txt').write_text('\n'.join(plain))
    checks = audit(content)
    provenance = dict(generated_at=stamp, producer='finished', renderer=str(Path(__file__).relative_to(ROOT)),
                      renderer_sha256=sha(Path(__file__)), **evidence, checks=checks,
                      chapters=content, pdf_sha256=sha(FIG / 'length_mechanisms.pdf'),
                      html_sha256=sha(FIG / 'index.html'), text_sha256=sha(FIG / 'report_text.txt'))
    for p in [REV / 'report_provenance.json', DATA / 'report_provenance.json',
              ROOT / 'data/runs/length_mechanisms_20261008/report_provenance.json']:
        p.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + '\n')
    (DATA / 'report_audit.json').write_text(json.dumps(checks, indent=2) + '\n')
    # The compatibility PDF is already a public entrypoint for this same report.
    old = ROOT / 'figures/three_test_mechanisms/three_test_mechanisms.pdf'
    if old.exists():
        shutil.copy2(FIG / 'length_mechanisms.pdf', old)
    print(json.dumps(dict(pdf=str(FIG / 'length_mechanisms.pdf'), html=str(FIG / 'index.html'),
                         generated_at=stamp, **checks), indent=2))


if __name__ == '__main__':
    main()
