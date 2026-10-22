"""Integrate Figure 4 round 3 into outline_v9, preserving the other five figures."""
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path

from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure4'
DEST = ROOT / 'manuscript_outline/revision_round3_figure4_20261022'


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.relative_to(ROOT))


def prepend(path, content):
    marker = '<!-- Figure 4 round 3 current edition -->'
    old = path.read_text()
    if old.startswith(marker):
        old = old.split('<!-- Earlier editions preserved below -->\n', 1)[1]
    path.write_text(marker + '\n' + content + '\n\n<!-- Earlier editions preserved below -->\n' + old)


def publication_copy(source):
    """Preserve the exact vector main and give its legend a readable separate page."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
    fig = plt.figure(figsize=(8.27, 11.69))
    renderer = fig.canvas.get_renderer()

    def wrap(paragraph, size, weight='normal'):
        font = FontProperties(family='DejaVu Sans', size=size, weight=weight)
        lines, line = [], ''
        for word in paragraph.split():
            trial = (line + ' ' + word).strip()
            if line and renderer.get_text_width_height_descent(trial, font, False)[0] > fig.get_figwidth()*fig.dpi*.85:
                lines.append(line)
                line = word
            else:
                line = trial
        return lines + ([line] if line else [])

    y = .94
    for line in wrap(source['title'], 12, 'bold'):
        fig.text(.075, y, line, fontsize=12, fontweight='bold', va='top', color='black')
        y -= .022
    y -= .022
    for paragraph in source['legend']:
        for line in wrap(paragraph, 9.7):
            fig.text(.075, y, line, fontsize=9.7, va='top', color='black')
            y -= .017
        y -= .014
    assert y > .06, 'Legend exceeds one publication page'
    fig.savefig(DEST / 'figure4_legend.pdf')
    plt.close(fig)
    combined = PdfWriter()
    combined.append(ROOT / source['main_pdf'])
    combined.append(DEST / 'figure4_legend.pdf')
    combined.write(DEST / 'figure4_with_legend.pdf')


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    baseline_path = ROOT / 'manuscript_outline/outline_v8.json'
    baseline = read(baseline_path)
    document = deepcopy(baseline)
    source = read(FIG / 'legend.json')
    section = next(s for s in document['sections'] if s['id'] == 'figure4')
    section.clear()
    section.update(source)
    section['title'] = 'Figure 4. ' + source['title'].removeprefix('Figure 4. ')
    section['role'] = 'Compare physical peak drift and persistent relocations as explanations for circumference-dependent activity-heading decay.'
    section['equations'] = [
        dict(math=r'k(\ell)=-\frac{\ln C(\ell)}{\ell},\qquad k\propto\bar C^{-\alpha}',
             number='4.1', description='C(ℓ) is activity-heading autocorrelation; C̄ is actual mean circumference. Endpoint eligibility uses |P| ≥ 0.1 with P = M₁/M₀.'),
        dict(math=r'\alpha=-\frac{\operatorname{Cov}(\ln\bar C,\ln k)}{\operatorname{Var}(\ln\bar C)}',
             number='4.2', description='The local slope uses the two- or three-size neighborhood and one shared eligible lag; the lag stays fixed in whole-seed bootstrap draws.'),
        dict(math=r'k_{\mathrm{pred}}=\left(\frac{2\pi}{\bar C}\right)^2D_s+\lambda\left(1-\langle\cos\Delta\theta_{\mathrm{jump}}\rangle\right)',
             number='4.3', description='The measured drift and persistent-switch terms form a closure prediction; disagreement with observed decay is retained as a result.')]
    generated = datetime.now(timezone.utc).isoformat()
    status = source.get('status', 'round3_campaign_pending')
    pending = 'pending' in status or 'provisional' in status or 'fallback' in status
    state = ('Figure 4 round 3: provisional round-2 model reanalysis; high-replication campaign pending'
             if pending else 'Figure 4 round 3: size-dependent heading persistence and drift/jump nulls')
    document.update(date='22 October 2026', stage=state,
        subtitle='Cell enlargement, polarity responses and heading persistence with drift/jump nulls')
    old = 'Enlarged cells hold heading longer under fluctuating input — at or above the geometric expectation at weak input — and the advantage declines as input strengthens. (Figure 4). Geometry and multiple-front fractions are descriptive: erosion below the single-front geometric expectation can also occur while a single front persists.'
    replacement = ('Figure 4 compares local size exponents and noise tolerance with drifting-peak and random-relocation nulls. '
        'It asks how much directly measured centroid drift and dominant-region switching explain heading decay. '
        + ('The displayed model estimates currently reuse round-2 seeds; the new 96–384-seed campaign is pending and its unavailable measurements are explicit. '
           if pending else 'Model estimates use the new compact high-replication campaign. ')
        + 'The revision uses no hypothesis gates; the preceding S0 record remains historical.')
    document['overview'] = [p.replace(old, replacement) for p in document['overview']]
    for supporting in document['supplementary']:
        if supporting.get('title') == 'Figure 4: complete revised evidence':
            supporting['title'] = 'Historical Figure 4: previous supplement and decision record'
            supporting['text'] = ['Historical evidence from the preceding Figure 4 edition, preserved unchanged. '
                                  + p for p in supporting['text']]
        if isinstance(supporting.get('text'), str):
            supporting['text'] = supporting['text'].replace(
                'Revised Figure 4 also displays', 'The historical Figure 4 edition also displays')
    protected = {}
    for ident in ('figure1', 'figure2', 'figure3', 'figure5', 'figure6'):
        previous = next(s for s in baseline['sections'] if s['id'] == ident)
        current = next(s for s in document['sections'] if s['id'] == ident)
        assert current == previous
        protected[ident] = dict(section_sha256=hashlib.sha256(json.dumps(previous, sort_keys=True).encode()).hexdigest(),
            asset_sha256={a['path']: sha(ROOT/a['path']) for a in previous['assets']})
    document['supplementary'].insert(0, dict(title='Figure 4: round-3 analysis and null comparisons',
        text=['Pointwise and smooth size exponents, drift/jump cross-checks, sensitivity, valid fractions, replication and analytic null checks. Missing or pending results are identified in the scientific readme.'],
        path=source['supplement_pdf'], sources=[dict(label='Unexpected results, data status and methods',
            path=relative(FIG/'supporting/metadata/scientific_readme.md'))]))
    write(ROOT/'manuscript_outline/outline_v9.json', document)
    publication_copy(section)
    combined = PdfWriter()
    for ident in ('figure2', 'figure3', 'figure4'):
        combined.append(ROOT/next(s for s in document['sections'] if s['id'] == ident)['main_pdf'])
    combined.write(DEST/'figures_2_4.pdf')
    receipt = dict(generated_utc=generated, baseline=relative(baseline_path), baseline_sha256=sha(baseline_path),
        outline='manuscript_outline/outline_v9.json', outline_sha256=sha(ROOT/'manuscript_outline/outline_v9.json'),
        figure4_status=status, campaign_complete=not pending, protected_sections=protected,
        all_other_main_sections_identical=True, previous_editions_preserved=True,
        main_page_dimensions_inches=[float(v)/72 for v in PdfReader(ROOT/section['main_pdf']).pages[0].mediabox[2:]],
        publication_pdf_sha256=sha(DEST/'figure4_with_legend.pdf'))
    write(DEST/'integration.json', receipt)
    relfig = os.path.relpath(FIG, DEST)
    status_text = ('The new model campaign is pending. This is a provisional round-2 model reanalysis and a tested rebuild pipeline; no new high-replication model results are claimed.'
        if pending else 'The compact high-replication campaign is complete.')
    body = f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 4 round 3</title>
<style>body{{max-width:1100px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}}img{{width:100%}}a{{color:#166aa2}}</style>
<h1>Figure 4 round 3 — manuscript integration</h1><p>Generated {generated}. {html.escape(status_text)}</p>
<p><a href="../outline_v9.pdf">Revised manuscript outline</a> · <a href="figure4_with_legend.pdf">Figure 4 and legend</a> · <a href="{relfig}/figure4_supplement.pdf">Supplement</a> · <a href="figures_2_4.pdf">Current main plates 2–4</a></p>
<img src="{relfig}/figure4_main.png" alt="Figure 4 drift and jump nulls and model measurements">
'''
    body += ''.join('<p>'+html.escape(p)+'</p>' for p in section['legend'])
    body += f'<p><a href="{relfig}/supporting/metadata/scientific_readme.md">Scientific readme and unexpected results</a> · <a href="integration.json">Preservation receipt</a></p></html>'
    (DEST/'index.html').write_text(body)
    content = f'''# Figure 4 round 3 — current manuscript outline

[**Revised outline_v9**](../outline_v9.pdf) · [**Figure 4 and legend**](figure4_with_legend.pdf) · [**Reader**](index.html) · [**Supplement**]({relfig}/figure4_supplement.pdf).

{status_text}

Figure 4 compares circumference-dependent heading decay with peak drift and random-relocation nulls. Model hues encode size, and the direct drift/jump decomposition distinguishes measured ingredients from null expectations. Figures 1, 2, 3, 5 and 6 retain their exact outline_v8 sections and artwork. Length represents WGD-associated enlargement; simulated polarity heading, observed migration summaries and conditional OU velocity-memory estimates remain distinct readouts. No hypothesis gates are regenerated.

[Scientific methods and data status]({relfig}/supporting/metadata/scientific_readme.md) · [Integration receipt](integration.json).

Data generation: {generated}; manuscript producer finished. No timed refresh.

Rebuild the figure with its single entry point, then this integration:

```sh
python -B polarity1d/manuscript/code/mechanisms/rebuild_figure4_round3_20261022.py
python -B manuscript_outline/prepare_figure4_round3_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v9.json --output manuscript_outline/outline_v9.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v9.json --report manuscript_outline/outline_v9.validation.json
```
'''
    (DEST/'README.md').write_text(content)
    for path, prefix in ((ROOT/'manuscript_outline/README.md', ''), (ROOT/'current manuscript figs/README.md', '../manuscript_outline/')):
        links = ('[**Revised outline_v9**](outline_v9.pdf) · [**Figure 4 and legend**](revision_round3_figure4_20261022/figure4_with_legend.pdf) · [**Reader and supplement**](revision_round3_figure4_20261022/index.html).'
            if not prefix else '[**Revised outline_v9**](../manuscript_outline/outline_v9.pdf) · [**Figure 4 and legend**](../manuscript_outline/revision_round3_figure4_20261022/figure4_with_legend.pdf) · [**Reader and supplement**](../manuscript_outline/revision_round3_figure4_20261022/index.html).')
        prepend(path, '# Figure 4 round 3 — nulls for size-dependent heading persistence\n\n'+links+'\n\n'+status_text+'\n\nFigures 1, 2, 3, 5 and 6 retain their exact preceding sections and artwork. Length represents WGD-associated enlargement. Observed migration summaries, conditional OU persistence estimates and simulated polarity-heading decay remain distinct readouts.')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
