"""Integrate Figure 2 round 3 as outline_v6, preserving every other main section.

outline_v5 and its completed round-2 package remain an immutable prior edition.
Rendering is a separate read-only assembly command documented in the README.
"""
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'manuscript_outline/revision_round3_20261022'
FIG = ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def relative(path):
    return str(path.relative_to(ROOT))


def prepend(path, text):
    marker = '<!-- Figure 2 round 3 current edition -->'
    previous = path.read_text()
    if previous.startswith(marker):
        previous = previous.split('<!-- Earlier editions preserved below -->\n', 1)[1]
    path.write_text(marker + '\n' + text + '\n\n<!-- Earlier editions preserved below -->\n' + previous)


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    baseline_path = ROOT / 'manuscript_outline/outline_v5.json'
    baseline = read(baseline_path)
    source = read(FIG / 'legend.json')
    assert source['status'] == 'round3_complete_qualified'
    assert read(FIG / 'supporting/metadata/build_receipt.json')['completed']
    assert read(FIG / 'supporting/metadata/supplement_preservation.json')['passed']
    document = deepcopy(baseline)
    document.update(date='22 October 2026',
        stage='Figure 2 round 3 integrated; other main figures retain completed round-2 content',
        subtitle='Working manuscript outline with revised size and pulse-memory signatures and complete retained supplements')
    two = next(s for s in document['sections'] if s['id'] == 'figure2')
    two.clear()
    two.update(source)
    two['title'] = 'Figure 2. ' + source['title']
    narrative = document['overview'][1]
    sentence = 'Enlarged cells resist imposed replacement through distinct initiation, spatial and temporal response signatures (Figure 2).'
    replacement = ('Figure 2 separates size-dependent initiation and cue-width responses from temporal class signatures: '
        'initiation-limited resetting, saturating activation with visible gap retention, and pulse-2 progress after loss of visible activity. '
        'The full Holmes width series shows finite replacement at all tested lengths for sufficiently broad cues; '
        'narrow-cue censoring cannot establish a growing penalty.')
    assert sentence in narrative
    document['overview'][1] = narrative.replace(sentence, replacement)
    for item in document['supplementary']:
        if isinstance(item, dict) and item.get('title') == 'Figure 2: complete revised evidence':
            item.update(title='Figure 2: round-3 evidence and preserved prior panels',
                text=['The complete supplement adds the per-width replacement curves and the deBelly native endpoint audit, then preserves the complete round-2 main and every earlier supplement page. Superseded fate maps, ratios and class tables are historical evidence.'],
                path=relative(FIG / 'figure2_round3_supplement.pdf'),
                sources=[dict(label='Round-3 reader and source data', path=relative(FIG / 'index.html')),
                    dict(label='Exact preservation receipt', path=relative(FIG / 'supporting/metadata/supplement_preservation.json'))])
    protected = {}
    for ident in ('figure1', 'figure3', 'figure4', 'figure5', 'figure6'):
        original = next(s for s in baseline['sections'] if s['id'] == ident)
        current = next(s for s in document['sections'] if s['id'] == ident)
        assert current == original
        protected[ident] = dict(section_sha256=hashlib.sha256(json.dumps(original, sort_keys=True).encode()).hexdigest(),
            asset_sha256={asset['path']: sha(ROOT / asset['path']) for asset in original['assets']})
    write(ROOT / 'manuscript_outline/outline_v6.json', document)
    # A compact figure-only report, with the same editable legend and equations.
    compact = dict(title='Figure 2 — round 3', date=document['date'], stage=document['stage'],
        subtitle='Size responses and pulse-memory signatures', overview=[], sections=[deepcopy(two)],
        supplementary=[dict(title='Complete supplementary evidence', text=['Superseded main panels and all previous supplements are preserved.'],
            sources=[dict(label='Full Figure 2 supplement', path=relative(FIG / 'figure2_round3_supplement.pdf'))])])
    write(DEST / 'figure2_with_legend.json', compact)
    combined = PdfWriter()
    for path in (FIG / 'figure2_round3.pdf',
        ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure3/figure3_main.pdf',
        ROOT / 'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4/figure4_main.pdf'):
        combined.append(path)
    combined.write(DEST / 'figures_2_4.pdf')
    generated = datetime.now(timezone.utc).isoformat()
    receipt = dict(generated_utc=generated, baseline=relative(baseline_path), baseline_sha256=sha(baseline_path),
        outline='manuscript_outline/outline_v6.json', outline_sha256=sha(ROOT / 'manuscript_outline/outline_v6.json'),
        new_figure2=source['main_pdf'], protected_sections=protected, native_runs_submitted=0,
        all_other_main_sections_identical=True, previous_edition_preserved=True,
        combined_main_pages=len(PdfReader(DEST / 'figures_2_4.pdf').pages))
    write(DEST / 'integration.json', receipt)
    rel_fig = os.path.relpath(FIG, DEST)
    body = '<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 2 round 3 in manuscript outline</title><style>body{max-width:1100px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}img{max-width:100%}a{color:#166aa2}</style>'
    body += f'<h1>Figure 2 round 3 integrated into the manuscript outline</h1><p>Generated {generated}. Figure and report producers complete.</p><p><a href="../outline_v6.pdf">Revised manuscript outline</a> · <a href="figure2_with_legend.pdf">Figure 2 and legend</a> · <a href="{rel_fig}/figure2_round3_supplement.pdf">Complete Figure 2 supplement</a> · <a href="figures_2_4.pdf">Main plates 2–4</a></p>'
    body += f'<p>{html.escape(document["overview"][1])}</p><p>Figures 1, 3, 4, 5 and 6 retain their preceding content and artwork.</p><a href="{rel_fig}/index.html"><img src="{rel_fig}/figure2_round3.png" alt="Figure 2 round 3"></a><p><a href="{rel_fig}/README.md">Scientific notes and all qualifications</a> · <a href="integration.json">Exact-source integration receipt</a></p></html>'
    (DEST / 'index.html').write_text(body)
    (DEST / 'README.md').write_text(f'''# Figure 2 — round 3 integration

[Revised outline_v6 PDF](../outline_v6.pdf) · [Editable outline](../outline_v6.json) · [Figure 2 and legend](figure2_with_legend.pdf) · [Reader](index.html) · [Main plates 2–4](figures_2_4.pdf)

The full revised Figure 2 replaces its preceding outline section and legend. A–C show size responses; D–E distinguish temporal classes; F gives per-class prospective biosensor predictions. The plate is 7.2 × 10 inches. All 450 width observations and 180 temporal coordinates are included. Broad Holmes cues already resolve every length, so no new simulations were needed. The deBelly L19/gap16 point was checked against native fields before plotting.

[Complete supplement]({rel_fig}/figure2_round3_supplement.pdf) preserves the superseded main, fate maps, ratios, tables and every earlier supplement page. [Source data and audit]({rel_fig}/README.md) include raw shares, visibility flags, exact endpoint checks and remaining numerical qualifications.

Figures 1, 3, 4, 5 and 6 are exact copies of their outline_v5 sections and use identical image files. The previous outline_v5 edition remains at its original path. [Integration receipt](integration.json).

Rebuild with the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure2_20261022.py
python -B polarity1d/manuscript/code/mechanisms/organize_length_dev.py --group revision_round3_20261022/figure2
python -B manuscript_outline/prepare_figure2_round3_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v6.json --output manuscript_outline/outline_v6.pdf
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/revision_round3_20261022/figure2_with_legend.json --output manuscript_outline/revision_round3_20261022/figure2_with_legend.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v6.json --report manuscript_outline/outline_v6.validation.json
python -B manuscript_outline/verify_figure2_round3_20261022.py
```

Generated UTC: {generated}. No timed browser refresh. No simulations or retuning are part of report assembly.
''')
    prepend(ROOT / 'manuscript_outline/README.md', '''# Figure 2 round 3 — current manuscript outline

[**Revised outline_v6**](outline_v6.pdf) · [**Figure 2 and legend**](revision_round3_20261022/figure2_with_legend.pdf) · [**Reader and complete supplement**](revision_round3_20261022/index.html) · [**Editable outline**](outline_v6.json).

Figure 2 now separates size responses (A–C) from temporal class signatures (D–E), with per-class biosensor predictions (F). It includes every width/area condition, all six deBelly initiation bounds and the independently checked deBelly L19/gap16 first gap. Broad Holmes cues resolve all lengths, so no additional runs or kinetic retuning were needed. The main is 7.2 × 10 inches. Superseded panels and every preceding supplement page are preserved. Figures 1, 3, 4, 5 and 6 retain their exact preceding sections and artwork.

[Integration and rebuild instructions](revision_round3_20261022/README.md).''')
    prepend(ROOT / 'current manuscript figs/README.md', '''# Figure 2 round 3 — current figure and outline

[**Revised manuscript outline**](../manuscript_outline/outline_v6.pdf) · [**Figure 2 and legend**](../manuscript_outline/revision_round3_20261022/figure2_with_legend.pdf) · [**Reader, source data and complete supplement**](../manuscript_outline/revision_round3_20261022/index.html).

Figure 2 now uses per-width replacement times at all three areas, first-gap retention versus pulse-2 headroom gain, challenger-share histories and per-class biosensor predictions. The deBelly L19/gap16 first-gap endpoint was checked against native fields. All 18 Holmes width-15.2 conditions are finite, so no new runs were required. The 7.2 × 10 inch main has uniform black labels and retains all superseded panels in the supplement. Figures 1, 3, 4, 5 and 6 retain their exact preceding outline content and artwork.

Length represents WGD-associated enlargement. Biosensor predictions are prospective; observed migration, conditional OU persistence and simulated polarity readouts remain distinct evidence.''')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
