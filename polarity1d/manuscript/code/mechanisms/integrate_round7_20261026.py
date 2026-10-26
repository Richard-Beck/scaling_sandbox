"""Build outline_v13 and preserve the preceding migration figure as Figure 6.

Run ``migration`` once, then ``integrate`` after Figures 3–5 are complete.
All scientific sources and the first two outline sections are immutable.
"""
from copy import deepcopy
from collections import Counter
from datetime import datetime, timezone
import argparse
import hashlib
import html
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
sys.dont_write_bytecode = True
ROOT = next(p for p in Path(__file__).resolve().parents if (p / 'manuscript_outline').is_dir())
DEST = ROOT / 'manuscript_outline/revision_round7_20261026'
FIGURES = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026'
BASELINE = ROOT / 'manuscript_outline/outline_v12.json'
PYTHON = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'manuscript_outline')]
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.mathtext import math_to_image
from pypdf import PdfReader, PdfWriter
import build_manuscript_outline as prose


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def rel(path):
    return str(Path(path).relative_to(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def migration():
    folder = FIGURES / 'figure6'
    for category in ('data', 'metadata', 'vectors', 'viewers'):
        (folder / 'supporting' / category).mkdir(parents=True, exist_ok=True)
    original = ROOT / 'manuscript_outline/assets/migration_v2'
    # The original plotting functions and numerical inputs are reused. Only the
    # figure number in the two suptitles and exported names changes.
    module = load_module('migration_round7_source', ROOT / 'scripts/render_manuscript_migration_v2.py')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.titlesize': 10, 'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.unicode_minus': True})
    source_hashes = {rel(p): sha(p) for p in original.iterdir() if p.is_file()}
    data = module.prepare_data(read(module.CANONICAL), folder / 'supporting/data')
    mapping = {'figure5_migration_v2': 'figure6_round7', 'figure5b_msd_detail': 'figure6B_detail_round7'}

    def save(fig, out, stem):
        fig._suptitle.set_text(fig._suptitle.get_text().replace('Figure 5', 'Figure 6'))
        for suffix in ('png', 'pdf', 'svg'):
            fig.savefig(out / f'{mapping[stem]}.{suffix}', dpi=230, facecolor='white', bbox_inches='tight')
        plt.close(fig)

    module.save = save
    module.integrated(data, folder)
    module.detail(data, folder)
    exported = {}
    for path in sorted((folder / 'supporting/data').glob('*.csv')):
        expected = original / path.name
        assert expected.is_file() and expected.read_bytes() == path.read_bytes(), path.name
        exported[path.name] = sha(path)
    section = deepcopy(next(s for s in read(BASELINE)['sections'] if s['id'] == 'figure5'))
    section.update(id='figure6', title=section['title'].replace('Figure 5', 'Figure 6'))
    for asset, stem in zip(section['assets'], mapping.values()):
        asset.update(path=rel(folder / f'{stem}.png'), native_pdf_page=rel(folder / f'{stem}.pdf'),
            page_size=[float(v) / 72 for v in PdfReader(folder / f'{stem}.pdf').pages[0].mediabox[2:]],
            page_title='', minimum_printed_width_inches=0)
        asset['caption'] = asset['caption'].replace('Figure 5', 'Figure 6')
    section['reader'] = rel(folder / 'scientific_readme.md')
    section['main_pdf'] = rel(folder / 'figure6_round7.pdf')
    write(folder / 'legend.json', section)
    (folder / 'main_legend.md').write_text('# Published enlarged-cell migration supports greater persistence, conditionally on the OU model\n\n' + '\n\n'.join(section['legend']) + '\n')
    (folder / 'scientific_readme.md').write_text('''# Unexpected

No new result, fit or simulation is introduced in Figure 6. The observed Xuan
population MSD crossover (2.32 h) differs from the unchanged OU prediction
(1.11 h). The pooled persistence fold of about 15.8 is conditional on the
hierarchy and priors; it does not establish enlargement as the cause.

## Preservation

This is round-6 Figure 5, renumbered. The existing renderer, all panel numerical
inputs, panel arrangement and caveats are retained. Its enlarged B detail stays
with the main figure. The inherited visual styling is retained under the explicit
instruction to leave this figure otherwise unchanged. Figure-number text and
output filenames are the only presentation changes. All compact plotted CSVs
are byte-identical to the preceding edition. The local provenance records source
and output hashes and the original model validation.

Observed migration summaries, conditional OU velocity-memory estimates, and
the simulated 1D polarity-heading readouts in Figures 3–5 are distinct evidence.

Rebuild using `integrate_round7_20261026.py migration` in the existing
`scaling_sandbox` environment. No external dependency was installed.
''')
    assert all(sha(ROOT / p) == h for p, h in source_hashes.items())
    write(folder / 'supporting/metadata/provenance.json', dict(generated_utc=now(),
        source_figure=5, destination_figure=6, new_fit=False, new_simulation=False,
        change='Figure-number text and filenames only; original plotting functions and data.',
        original_source_sha256=source_hashes, plotted_csv_sha256=exported,
        original_validation=data['validation'],
        output_sha256={p.name: sha(p) for p in folder.iterdir() if p.suffix in ('.png', '.svg', '.pdf')}))
    print(folder)


def legend_pdf(section):
    ident = section['id']
    path = DEST / f'{ident}_legend.pdf'
    builder = prose.OutlineBuilder({'title': section['title'], 'sections': []}, ROOT,
        path, ROOT / 'manuscript_outline/outline_v13.json')
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
    for asset in section['assets']:
        writer.append(ROOT / asset['native_pdf_page'])
    writer.append(path)
    writer.write(DEST / f'{ident}_with_legend.pdf')


def equations_html(section):
    folder = DEST / 'equations'
    folder.mkdir(exist_ok=True)
    body = ''
    for index, equation in enumerate(section.get('equations', []), 1):
        equation = {'math': equation} if isinstance(equation, str) else equation
        path = folder / f'{section["id"]}_{index}.png'
        math_to_image('$' + equation['math'].strip('$') + '$', path,
            prop=FontProperties(size=16), dpi=180, format='png', color='black')
        number = html.escape(str(equation.get('number', index)))
        body += f'<div class="equation"><img src="equations/{path.name}" alt="Equation {number}"><span>({number})</span></div>'
        body += '<p>' + html.escape(equation.get('description', '')) + '</p>'
    return body


def integrate():
    baseline = read(BASELINE)
    document = deepcopy(baseline)
    document.update(date='26 October 2026', stage='Round 7: cytosolic conduction and heading relocations',
        subtitle='Replacement resistance, conduction-limited response and stochastic polarity persistence')
    sections = deepcopy(baseline['sections'][:2])
    for n in range(3, 7):
        section = read(FIGURES / f'figure{n}/legend.json')
        section['id'] = f'figure{n}'
        for equation in section.get('equations', []):
            if isinstance(equation, dict):
                equation['math'] = re.sub(r'\\mathrm\s+([A-Za-z])', r'\\mathrm{\1}', equation['math'])
        sections.append(section)
    reserved = deepcopy(baseline['sections'][5])
    reserved.update(id='figure7', title=reserved['title'].replace('Figure 6', 'Figure 7'))
    for asset in reserved['assets']:
        asset['page_title'] = asset.get('page_title', '').replace('Figure 6', 'Figure 7')
    sections.append(reserved)
    document['sections'] = sections
    document['supplementary'] = []
    document['overview'] = [
        'Working hypothesis: WGD-associated enlargement changes established polarity with local kinetics held fixed. Length represents enlargement. Published migration motivates this hypothesis; model heading and migration persistence remain distinct readouts.',
        'Figures 1–2 retain the preceding evidence unchanged. Larger cells require longer fixed-width and proportional-width cues for replacement. Initiation, domain expansion and pulse-train retention distinguish the response mechanisms; censored results and nonmonotone deBelly responses remain qualified.',
        'Figure 3 separates strong-cue turning from weak boundary response. The conduction test compares exact ring gains with an isolated-front prediction for Wave-pinning, inactive-diffusivity sweeps in Holmes and Spring, and the exceptional steep gain decline in Goryachev. Strong-cue turning times are not inferred from weak-input gain.',
        'Figure 4 opens with drift and jump baselines, then shows heading correlation and the local size exponent across input strengths. Greater persistence of larger rings is tested against geometry rather than identified with cell movement.',
        'Figure 5 compares measured drift with deterministic response predictions, tests drift-plus-takeover closure, and displays unfiltered heading-change distributions. Strong-input patch fields and native/fast-cytosol mode rates test whether long-range transfer contributes to persistence; seed uncertainty and discrepancies with the pilots remain explicit.',
        'Figure 6 is the preceding migration Figure 5, renumbered with its scientific content unchanged. Published migration summaries support greater persistence in enlarged groups, with conditional OU estimates and visible model inadequacy.',
        'Figure 7 reserves the SUM-159 paired WGD motility assays. No new experimental observations are available. Boundary photoactivation, full-perimeter activity and cytosolic mobility are prospective discriminating measurements, not completed assays.'
    ]
    document['references'] = [r for r in baseline['references'] if 'Claude' not in r.get('label', '') and 'response' not in r.get('label', '').lower()]
    document['references'] += [
        {'label': 'Round-7 brief', 'path': rel(DEST / 'request/brief_round7.md')},
        {'label': 'Target narrative; quantitative statements checked against the rebuilt evidence', 'path': rel(DEST / 'request/target_narrative_round7.md')},
        {'label': 'Pilot methods and comparison values; independent reimplementation', 'path': rel(DEST / 'request/pilot_analyses_round6.md')},
    ]
    response = DEST / 'claude_response.md'
    if response.exists():
        document['references'].append({'label': 'Response to round-7 Claude review', 'path': rel(response)})
    path = ROOT / 'manuscript_outline/outline_v13.json'
    write(path, document)
    # Frozen Figures 1–2 retain their exact section records, artwork and data.
    assert document['sections'][:2] == baseline['sections'][:2]
    preserved = {}
    for section in document['sections'][:2]:
        preserved[section['id']] = {'section_sha256': hashlib.sha256(json.dumps(section, sort_keys=True).encode()).hexdigest(),
            'assets': {a['path']: sha(ROOT / a['path']) for a in section['assets']}}
    for section in sections[:6]:
        legend_pdf(section)
    writer = PdfWriter()
    for n in range(1, 7):
        writer.append(DEST / f'figure{n}_with_legend.pdf')
    writer.write(DEST / 'figures1_6_with_legends.pdf')
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Round 7 manuscript</title>
<style>body{{max-width:1400px;margin:2rem auto;padding:0 1rem;font:17px/1.55 system-ui;color:#17212b}}img{{max-width:100%;height:auto}}a{{color:#175d98}}.status{{padding:1rem;background:#eef3f7}}.equation{{display:flex;align-items:center;justify-content:center;position:relative;margin:1.5rem 0}}.equation img{{max-width:85%;width:auto}}.equation span{{position:absolute;right:0}}</style>
<h1>Size, cytosolic conduction and polarity persistence</h1><p class="status">Generated {now()}. Document producer finished; rebuild to load revised data. No timed refresh.</p>
<p><a href="../outline_v13.pdf">Integrated outline_v13</a> · <a href="figures1_6_with_legends.pdf">Figures 1–6 with legends</a> · <a href="request/brief_round7.md">Brief</a></p>'''
    body += '<p>Pilot checks: ' + ' · '.join(
        f'<a href="{os.path.relpath(FIGURES / name, DEST)}">{label}</a>'
        for name, label in (
            ('figure3/pilot_conduction_penalty.png', 'Conduction penalty against size/diffusivity'),
            ('pilots/holmes_mode_rates_native_fast_seed_CI.png', 'Holmes native/fast mode rates with seed intervals'))) + '</p>'
    for section in sections[:6]:
        ident = section['id']
        body += '<h2>' + html.escape(section['title']) + '</h2>'
        body += f'<p><a href="{ident}_with_legend.pdf">Figure and legend</a>'
        if section.get('reader'):
            body += f' · <a href="{os.path.relpath(ROOT / section["reader"], DEST)}">Scientific methods and results</a>'
        body += '</p>'
        for asset in section['assets']:
            body += f'<img src="{os.path.relpath(ROOT / asset["path"], DEST)}" alt="{html.escape(section["title"], quote=True)}">'
        body += ''.join('<p>' + html.escape(p) + '</p>' for p in section['legend'])
        body += equations_html(section)
    body += '<h2>Figure 7 — reserved</h2><p>SUM-159 WGD motility assays; no observations available.</p>'
    body += '<p><a href="integration.json">Integration receipt</a> · <a href="acceptance.json">Validation</a></p>'
    if response.exists():
        body += '<p><a href="claude_response.md">Response to Claude review</a></p>'
    body += '</html>'
    (DEST / 'index.html').write_text(body)
    (DEST / 'README.md').write_text('''# Round 7 manuscript revision

[outline_v13](../outline_v13.pdf) · [Figures 1–6 with legends](figures1_6_with_legends.pdf) · [Reader](index.html).

Figures 3–5 follow the new conduction and relocation narrative. Figures 1–2
retain their exact outline sections, artwork and numerical evidence. The previous
migration Figure 5 is Figure 6; Figure 7 is reserved. Dropped panels and earlier
supplements are absent from this edition. Historical files are preserved.

Figure methods open with Unexpected. Compact CSVs, provenance and execution
records accompany the figures under `revision_round7_20261026`. The pilot PNGs
are checks, not manuscript panels. No output was tuned to pilot values.

Rebuild the individual figures first, then use the existing scaling_sandbox
environment to run `integrate_round7_20261026.py integrate`,
`build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v13.json --output manuscript_outline/outline_v13.pdf`,
and `integrate_round7_20261026.py validate` (paths are repository relative;
the integration script is in `polarity1d/manuscript/code/mechanisms/`).
''')
    write(DEST / 'integration.json', dict(generated_utc=now(), baseline_commit='122d306b',
        baseline_json_sha256=sha(BASELINE), baseline_pdf_sha256=sha(BASELINE.with_suffix('.pdf')),
        outline_json_sha256=sha(path), protected_sections=preserved,
        scientific_sections=7, supplement_included=False, original_data_removed=False))
    update_navigation()
    print(path)


def update_navigation():
    begin, end = '<!-- round7-current-edition -->', '<!-- round7-current-edition end -->'
    for path, prefix in ((ROOT / 'manuscript_outline/README.md', ''),
                         (ROOT / 'current manuscript figs/README.md', '../manuscript_outline/')):
        intro = f'''{begin}
# Size, conduction and polarity persistence

[**Current outline_v13**]({prefix}outline_v13.pdf) · [**Figures 1–6 with legends**]({prefix}revision_round7_20261026/figures1_6_with_legends.pdf) · [**Reader and methods**]({prefix}revision_round7_20261026/index.html).

Figures 1–2 are unchanged. Figures 3–5 test cytosolic conduction, boundary
response, drift and heading relocations. Published migration is renumbered
Figure 6; paired WGD experiments remain reserved as Figure 7.

{end}
'''
        old = path.read_text()
        if begin in old:
            before, rest = old.split(begin, 1)
            _, after = rest.split(end, 1)
            path.write_text(before + intro + after)
        else:
            path.write_text(intro + '\n' + old)
    path = ROOT / 'README.md'
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith('[**Current manuscript outline'):
            lines[index] = '[**Current manuscript outline — 26 October 2026**](manuscript_outline/outline_v13.pdf) · [**Figures 1–6 with legends**](manuscript_outline/revision_round7_20261026/figures1_6_with_legends.pdf) · [**Reader**](manuscript_outline/revision_round7_20261026/index.html).'
            break
    path.write_text('\n'.join(lines) + '\n')


def validate():
    document_path = ROOT / 'manuscript_outline/outline_v13.json'
    document = read(document_path)
    baseline = read(BASELINE)
    assert document['sections'][:2] == baseline['sections'][:2]
    assert [s['id'] for s in document['sections']] == [f'figure{n}' for n in range(1, 8)]
    assert not document['supplementary']
    assert 'RESERVED' in document['sections'][-1]['title']
    manifest = read(document_path.with_suffix('.build_manifest.json'))
    expected_assets = Counter(a['path'] for s in document['sections'] for a in s['assets'])
    embedded_assets = Counter(a for p in manifest['page_inventory'] for a in p['assets'])
    assert expected_assets == embedded_assets
    for row in manifest['inputs']:
        assert sha(ROOT / row['path']) == row['sha256'], row['path']
    pdf = ROOT / manifest['pdf']['path']
    assert sha(pdf) == manifest['pdf']['sha256']
    reader = PdfReader(pdf)
    assert len(reader.pages) == manifest['pdf']['pages']
    assert all((page.extract_text() or '').strip() for page in reader.pages)
    assert all(p['bounds_valid'] and p['text_overlap_check'] == 'passed' for p in manifest['page_inventory'])
    for section in document['sections']:
        assert all(section.get(key) for key in ('title', 'role', 'assets', 'legend', 'narrative', 'measurements', 'limitations', 'sources'))
        assert any(p['section_id'] == section['id'] and p['kind'] == 'figure_notes' for p in manifest['page_inventory'])
    for n in range(3, 7):
        folder = FIGURES / f'figure{n}'
        for suffix in ('png', 'pdf', 'svg'):
            assert list(folder.glob(f'*.{suffix}')), (n, suffix)
        assert (folder / 'main_legend.md').is_file()
        assert (folder / 'scientific_readme.md').read_text().lstrip().startswith('# Unexpected')
        assert list((folder / 'supporting/data').glob('*.csv')), n
    noise = read(FIGURES / 'noise/status.json')
    assert noise['all_jobs_complete'] and noise['production_jobs'] == 1536 and noise['control_jobs'] == 32
    assert noise['seeds_per_condition'] == 64 and noise['paired_noise_hashes_identical']
    assert noise['seed_mode_censoring_counts'] == {'resolved': 24576}
    for n in (4, 5):
        checks = read(FIGURES / f'figure{n}/supporting/metadata/validation.json')
        assert not checks['display_text_outside_page'] and not checks['new_supplement_created']
        if n == 5:
            assert checks['all_finite_lag_pairs_retained'] and checks['raw_heading_histograms'] == 48
    protected = [rel(BASELINE), rel(BASELINE.with_suffix('.pdf'))]
    protected += [a['path'] for s in baseline['sections'][:2] for a in s['assets']]
    protected += [a['native_pdf_page'] for s in baseline['sections'][:2] for a in s['assets']]
    for name in protected:
        original = subprocess.check_output(['git', 'show', f'122d306b:{name}'])
        assert hashlib.sha256(original).hexdigest() == sha(ROOT / name), name
    write(DEST / 'acceptance.json', dict(generated_utc=now(), passed=True,
        figures1_2_exact_sections=True, sections=7, figure7_reserved=True,
        no_supplement=True, build_input_hashes_verified=True, pdf_pages=len(reader.pages),
        pdf_bounds_and_text_overlap_passed=True, figure3_6_exports_and_readmes_present=True,
        production_noise_seeds=1536, numerical_controls=32, paired_forcing_verified=True,
        raw_heading_histograms=48, all_finite_heading_pairs_retained=True,
        protected_files_verified_against_base_commit=len(protected)))
    print(DEST / 'acceptance.json')


def storage():
    started = now()
    groups = [
        (FIGURES, 'Round-7 Figures 3–6, complete conduction and paired strong-input noise experiments, retained restart/frozen sources, compact histories, source CSVs and provenance. Original scientific data preserved.'),
        (DEST, 'Round-7 request, integrated figure-and-legend PDF, reader, acceptance and review response.'),
        (ROOT / 'manuscript_outline/outline_v13.pdf', 'Integrated seven-section outline; Figures 1–2 protected, migration renumbered to Figure 6, Figure 7 reserved; no supplement.'),
        (ROOT / '.git', 'Versioned manuscript source, deliverables, compact campaign evidence and repository history. This is not an independently verified external backup.'),
    ]
    records = []
    for path, purpose in groups:
        if path.exists():
            count = int(subprocess.check_output(['du', '-sx', '-B1', str(path)], text=True).split()[0])
            records.append(dict(path=rel(path), allocated_bytes=count, allocated_GiB=count/2**30, purpose_and_status=purpose))
    result = dict(measurement_started_utc=started, measurement_finished_utc=now(),
        method='du -sx -B1; allocated blocks, targeted non-atomic traversal; parent/child groups overlap.',
        groups=records, original_data_removed=False)
    path = ROOT / 'reports/round7_storage_20261026.json'
    write(path, result)
    begin, end = '<!-- round7-storage-20261026 begin -->', '<!-- round7-storage-20261026 end -->'
    content = f'''{begin}

## Round 7 conduction and heading mechanisms, 26 October 2026

Measured **{result['measurement_finished_utc']}** using allocated storage
(`du -sx -B1`). This targeted non-atomic traversal updates the named groups;
earlier whole-workspace totals remain historical. Inventory does not authorize
deletion. No original scientific outputs were removed.

| Path / related group | Measured allocated storage | Purpose and scientific status |
| --- | ---: | --- |
'''
    for row in records:
        content += f"| `{row['path']}` | {row['allocated_GiB']:.6f} GiB ({row['allocated_bytes']:,} bytes) | {row['purpose_and_status']} |\n"
    content += f'\n[Machine-readable measurement]({rel(path)}). The campaign manifests and seed-level status records distinguish execution completion from statistical precision.\n\n{end}\n'
    path = ROOT / 'DATA_MANIFEST.md'
    old = path.read_text()
    if begin in old:
        before, rest = old.split(begin, 1)
        _, after = rest.split(end, 1)
        path.write_text(before + content + after)
    else:
        path.write_text(old.rstrip() + '\n\n' + content)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('command', choices=('migration', 'integrate', 'validate', 'storage'))
    args = parser.parse_args()
    globals()[args.command]()
