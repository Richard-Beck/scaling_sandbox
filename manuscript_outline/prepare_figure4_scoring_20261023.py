"""Integrate revised Figure 4 scoring into outline_v10; preserve all other figures."""
from copy import deepcopy
from datetime import datetime,timezone
import hashlib
import html
import json
import os
from pathlib import Path
from pypdf import PdfWriter
import prepare_figure4_round3_20261022 as previous

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4'
DEST=ROOT/'manuscript_outline/revision_figure4_scoring_20261023'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def relative(p):return str(p.relative_to(ROOT))


def main():
    DEST.mkdir(parents=True,exist_ok=True);baseline_path=ROOT/'manuscript_outline/outline_v9.json'
    baseline=read(baseline_path);document=deepcopy(baseline);source=read(FIG/'legend.json')
    receipt=read(FIG/'supporting/metadata/scoring_revision_receipt.json')
    events=next(r for r in receipt['closure_by_regime'] if r['method']=='scorable_events' and r['regime']=='all')
    multiple=next(r for r in receipt['closure_by_regime'] if r['method']=='scorable_events' and r['regime']=='multiple_front_majority')
    mixture=next(r for r in receipt['closure_by_regime'] if r['method']=='detector_free' and r['regime']=='multiple_front_majority')
    status=(f'Figure 4 now normalizes physical drift by the smallest size at the same input and estimates jump angles from scorable events. '
            f'The existing {receipt["counted_replacements"]:,} counted replacements include {receipt["scored_replacements"]:,} scored angles. '
            f'Event-based closure is estimable in {events["estimable"]} conditions; {events["within_15_percent"]} are within 15% of observed decay. '
            f'Among conditions with multiple fronts for most frames, {multiple["within_15_percent"]}/{multiple["estimable"]} event-based and '
            f'{mixture["within_15_percent"]}/{mixture["estimable"]} detector-free comparisons are within 15%. '
            'The decomposition therefore retains its failures in the regime it seeks to explain. No new simulations were run.')
    section=next(s for s in document['sections'] if s['id']=='figure4');old_section=deepcopy(section)
    section.clear();section.update(source);section['title']='Figure 4. '+source['title']
    section['equations']=deepcopy(old_section['equations'])
    section['equations'].append(dict(math=r'\widehat{\lambda}_{\mathrm{eff}}=\frac{N_{\mathrm{all}}}{T_{\mathrm{dom}}}\,\frac{1}{N_{\mathrm{scored}}}\sum_{e\in\mathrm{scored}}(1-\cos\Delta\theta_e)',
        number='4.4',description='All counted replacements supply the rate; scorable events supply the conditional angle distribution. The supplement reports scoring fractions and missing-angle bounds.'))
    document.update(date=datetime.now(timezone.utc).strftime('%d %B %Y'),stage='Figure 4: corrected event scoring and normalized drift',
                    subtitle='Cell enlargement, polarity responses and heading persistence with corrected drift/jump comparisons')
    start='Figure 4 compares local size exponents and noise tolerance with drifting-peak and random-relocation nulls.'
    end='The revision uses no hypothesis gates; the preceding S0 record remains historical.'
    for i,paragraph in enumerate(document['overview']):
        if start in paragraph:
            a=paragraph.index(start);b=paragraph.index(end,a)+len(end)
            document['overview'][i]=paragraph[:a]+status+' '+paragraph[b:]
    assert any(status in p for p in document['overview'])
    document['supplementary'].insert(0,dict(title='Figure 4: corrected event scoring and full closure comparisons',
        text=['Scorable-event coverage, missing-angle bounds, both detector-free cadences and all closure failures. Preceding editions remain linked below.'],
        path=source['supplement_pdf'],sources=[dict(label='Scientific methods and regime-specific counts',path=relative(FIG/'supporting/metadata/scientific_readme.md'))]))
    protected={}
    for ident in ('figure1','figure2','figure3','figure5','figure6'):
        old=next(s for s in baseline['sections'] if s['id']==ident);new=next(s for s in document['sections'] if s['id']==ident)
        assert old==new
        protected[ident]=dict(section_sha256=hashlib.sha256(json.dumps(old,sort_keys=True).encode()).hexdigest(),
                              asset_sha256={a['path']:sha(ROOT/a['path']) for a in old['assets']})
    write(ROOT/'manuscript_outline/outline_v10.json',document)
    previous.DEST=DEST;previous.publication_copy(section)
    combined=PdfWriter()
    for ident in ('figure2','figure3','figure4'):combined.append(ROOT/next(s for s in document['sections'] if s['id']==ident)['main_pdf'])
    combined.write(DEST/'figures_2_4.pdf')
    generated=datetime.now(timezone.utc).isoformat();relfig=os.path.relpath(FIG,DEST)
    integration=dict(generated_utc=generated,baseline=relative(baseline_path),baseline_sha256=sha(baseline_path),
                     outline='manuscript_outline/outline_v10.json',outline_sha256=sha(ROOT/'manuscript_outline/outline_v10.json'),
                     protected_sections=protected,all_other_main_sections_identical=True,no_new_simulations=True,
                     publication_pdf_sha256=sha(DEST/'figure4_with_legend.pdf'),scoring_receipt_sha256=sha(FIG/'supporting/metadata/scoring_revision_receipt.json'))
    write(DEST/'integration.json',integration)
    links=f'[**Revised outline_v10**](../outline_v10.pdf) · [**Figure 4 and legend**](figure4_with_legend.pdf) · [**Reader**](index.html) · [**Supplement**]({relfig}/figure4_supplement.pdf)'
    rebuild='''```sh
python -B polarity1d/manuscript/code/mechanisms/figure4_round4_reanalysis_20261023.py
python -B manuscript_outline/prepare_figure4_scoring_20261023.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v10.json --output manuscript_outline/outline_v10.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v10.json --report manuscript_outline/outline_v10.validation.json
python -B manuscript_outline/verify_figure4_scoring_20261023.py
```'''
    (DEST/'README.md').write_text('# Figure 4: corrected event scoring and normalized drift\n\n'+links+'\n\n'+status+
        '\n\nFigures 1, 2, 3, 5 and 6 preserve their exact outline_v9 sections and artwork. Length represents WGD-associated enlargement; activity heading, migration summaries and conditional OU persistence estimates remain distinct evidence.\n\n'+
        f'[Scientific methods and results]({relfig}/supporting/metadata/scientific_readme.md). Generated {generated}; producer finished.\n\n'+rebuild+'\n')
    (FIG/'README.md').write_text('# Figure 4: corrected event scoring and normalized drift\n\n'+status+
        '\n\n[Main](figure4_main.pdf) · [Supplement](figure4_supplement.pdf) · [Reader](index.html) · [Methods and results](supporting/metadata/scientific_readme.md).\n')
    body=f'''<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 4: revised scoring</title>
<style>body{{max-width:1100px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}}img{{width:100%}}a{{color:#166aa2}}</style>
<h1>Figure 4: corrected event scoring and normalized drift</h1><p>Generated {generated}; producer finished.</p><p>{html.escape(status)}</p>
<p><a href="../outline_v10.pdf">Revised manuscript outline</a> · <a href="figure4_with_legend.pdf">Figure 4 and legend</a> · <a href="{relfig}/figure4_supplement.pdf">Supplement</a></p>
<img src="{relfig}/figure4_main.png" alt="Normalized physical drift, event-based and detector-free jump and closure estimates">
'''
    body+=''.join('<p>'+html.escape(p)+'</p>' for p in source['legend'])
    body+=f'<p><a href="{relfig}/supporting/metadata/scientific_readme.md">Methods and all regime-specific results</a> · <a href="integration.json">Preservation receipt</a></p></html>'
    (DEST/'index.html').write_text(body)
    marker='<!-- Figure 4 corrected scoring current edition -->'
    for path,prefix in [(ROOT/'manuscript_outline/README.md',''),(ROOT/'current manuscript figs/README.md','../manuscript_outline/')]:
        old=path.read_text()
        if old.startswith(marker):old=old.split('<!-- Earlier editions preserved below -->\n',1)[1]
        local_links=f'[**Revised outline_v10**]({prefix}outline_v10.pdf) · [**Figure 4 and legend**]({prefix}revision_figure4_scoring_20261023/figure4_with_legend.pdf) · [**Reader and supplement**]({prefix}revision_figure4_scoring_20261023/index.html)'
        path.write_text(marker+'\n# Figure 4: corrected event scoring and normalized drift\n\n'+local_links+'\n\n'+status+
                        '\n\nFigures 1, 2, 3, 5 and 6 preserve their preceding sections and artwork.\n\n<!-- Earlier editions preserved below -->\n'+old)
    print(json.dumps(integration,indent=2))


if __name__=='__main__':main()
