"""Integrate the single-C Figure 3 layout into outline_v8; preserve every other section."""
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import html
import json
import os
from pathlib import Path
from pypdf import PdfReader, PdfWriter

ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'manuscript_outline/revision_round3_figure3_single_C_20261022'
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3_single_C'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def write(path,value):path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def relative(path):return str(path.relative_to(ROOT))


def prepend(path,text):
    marker='<!-- Figure 3 single-C layout current edition -->'
    previous=path.read_text()
    if previous.startswith(marker):previous=previous.split('<!-- Earlier editions preserved below -->\n',1)[1]
    path.write_text(marker+'\n'+text+'\n\n<!-- Earlier editions preserved below -->\n'+previous)


def compact_legend(source):
    """Two-page publication copy: exact main plate followed by its legend."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
    fig=plt.figure(figsize=(8.27,11.69))
    renderer=fig.canvas.get_renderer()
    def wrap(text,size,bold=False):
        font=FontProperties(family='DejaVu Sans',size=size,weight='bold' if bold else 'normal')
        limit=fig.get_figwidth()*fig.dpi*.85
        lines=[];current=''
        for word in text.split():
            candidate=(current+' '+word).strip()
            width=renderer.get_text_width_height_descent(candidate,font,False)[0]
            if width>limit and current:lines.append(current);current=word
            else:current=candidate
        if current:lines.append(current)
        return lines
    title='Figure 3. '+source['title']
    y=.94
    for line in wrap(title,12,True):
        fig.text(.075,y,line,fontsize=12,fontweight='bold',va='top',color='black')
        y-=.022
    y-=.024
    for paragraph in source['legend']:
        for line in wrap(paragraph,10.2):
            fig.text(.075,y,line,fontsize=10.2,va='top',color='black')
            y-=.018
        y-=.016
    assert y>.07,'Legend exceeds one publication page'
    legend=DEST/'figure3_legend.pdf'
    fig.savefig(legend);plt.close(fig)
    writer=PdfWriter();writer.append(FIG/'figure3_round3.pdf');writer.append(legend)
    writer.write(DEST/'figure3_with_legend.pdf')
    write(DEST/'figure3_with_legend.build_manifest.json',dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        generator=relative(Path(__file__)),
        pdf=dict(path=relative(DEST/'figure3_with_legend.pdf'),pages=2,
            sha256=sha(DEST/'figure3_with_legend.pdf')),
        inputs=[dict(path=relative(p),sha256=sha(p)) for p in
            (Path(__file__),FIG/'legend.json',FIG/'figure3_round3.pdf')],
        layout='Exact native figure page followed by a single publication legend page; font-measured wrapping.'))


def main():
    DEST.mkdir(parents=True,exist_ok=True)
    baseline_path=ROOT/'manuscript_outline/outline_v7.json'
    baseline=read(baseline_path); document=deepcopy(baseline)
    source=read(FIG/'legend.json')
    assert read(FIG/'supporting/metadata/build_receipt.json')['completed']
    document.update(date='22 October 2026',
        stage='Figure 3 single-C landscape layout; Figure 2 round 3 and all other sections preserved',
        subtitle='Size-dependent polarity responses and cue-placement regimes with complete retained evidence')
    section=next(s for s in document['sections'] if s['id']=='figure3')
    section.clear();section.update(source);section['title']='Figure 3. '+source['title']
    old='Distance between competing regions matters only at short range; enlarged cells resist replacement and resolve coexistence slowly mainly because of size itself (Figure 3).'
    replacement='Cue placement relative to the incumbent’s edge switches reorientation between sliding the incumbent and takeover by a new region; enlargement slows both. Figure 3 measures time to 50% of cue-direction heading under a continuous cue and the corresponding mean angular speed; its native profiles distinguish the two regimes.'
    for i,text in enumerate(document['overview']):
        if old in text:document['overview'][i]=text.replace(old,replacement)
        elif 'Spring’s middle size' in text or "Spring's middle size" in text:
            document['overview'][i]=text.replace("Spring's middle size is the supported exception in which extra accessible separation adds a substantial resolution-time contribution. ",'').replace('Spring’s middle size is the supported exception in which extra accessible separation adds a substantial resolution-time contribution. ','')
    # The old Spring sentence can occur in the same overview paragraph.
    document['overview']=[t.replace("Spring's middle size is the supported exception in which extra accessible separation adds a substantial resolution-time contribution. ",'') for t in document['overview']]
    for item in document['supplementary']:
        if isinstance(item,dict) and item.get('title','').startswith('Figure 3:'):
            item.update(title='Figure 3: round-3 evidence and every preserved round-2 panel',
                text=['The supplement retains 90/99%-heading time curves, the preceding four-placement round-3 main and legend including C = 47.6 profiles, and the complete round-2 main and supplement. Original source files remain unchanged.'],
                path=relative(FIG/'figure3_round3_supplement.pdf'),
                sources=[dict(label='Round-3 reader and native source coordinates',path=relative(FIG/'index.html')),
                    dict(label='Exact supplement preservation',path=relative(FIG/'supporting/metadata/supplement_preservation.json'))])
    protected={}
    for ident in ('figure1','figure2','figure4','figure5','figure6'):
        original=next(s for s in baseline['sections'] if s['id']==ident)
        current=next(s for s in document['sections'] if s['id']==ident)
        assert current==original
        protected[ident]=dict(section_sha256=hashlib.sha256(json.dumps(original,sort_keys=True).encode()).hexdigest(),
            asset_sha256={asset['path']:sha(ROOT/asset['path']) for asset in original['assets']})
    write(ROOT/'manuscript_outline/outline_v8.json',document)
    compact=dict(title='Figure 3 — single-C layout',date=document['date'],stage=document['stage'],
        subtitle='Cue-placement regimes and half-turn readouts',overview=[],sections=[deepcopy(section)],supplementary=[])
    write(DEST/'figure3_with_legend.json',compact)
    compact_legend(source)
    combined=PdfWriter()
    for ident in ('figure2','figure3','figure4'):
        sec=next(s for s in document['sections'] if s['id']==ident)
        combined.append(ROOT/sec['main_pdf'])
    combined.write(DEST/'figures_2_4.pdf')
    generated=datetime.now(timezone.utc).isoformat()
    receipt=dict(generated_utc=generated,baseline=relative(baseline_path),baseline_sha256=sha(baseline_path),
        outline='manuscript_outline/outline_v8.json',outline_sha256=sha(ROOT/'manuscript_outline/outline_v8.json'),
        new_figure3=source['main_pdf'],protected_sections=protected,
        all_other_main_sections_identical=True,previous_edition_preserved=True,
        combined_main_pages=len(PdfReader(DEST/'figures_2_4.pdf').pages))
    write(DEST/'integration.json',receipt)
    relfig=os.path.relpath(FIG,DEST)
    body=f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 3 round 3 manuscript integration</title><style>body{{max-width:1100px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}}img{{width:100%}}a{{color:#166aa2}}</style><h1>Figure 3 round 3 integrated</h1><p>Generated {generated}. All source producers finished; layout revision complete without new simulations.</p>'
    body+=f'<p><a href="../outline_v8.pdf">Revised manuscript outline</a> · <a href="figure3_with_legend.pdf">Figure 3 and legend</a> · <a href="{relfig}/figure3_round3_supplement.pdf">Complete supplement</a> · <a href="figures_2_4.pdf">Current main plates 2–4</a></p><h2>{html.escape(source["title"])}</h2><a href="{relfig}/index.html"><img src="{relfig}/figure3_round3.png" alt="Figure 3 round 3"></a>'
    body+=''.join('<p>'+html.escape(p)+'</p>' for p in source['legend'])
    body+=f'<p><a href="{relfig}/README.md">Native data, methods and checks</a> · <a href="integration.json">Preservation receipt</a></p></html>'
    (DEST/'index.html').write_text(body)
    (DEST/'README.md').write_text(f'''# Figure 3 round 3 manuscript integration

[Revised outline_v8](../outline_v8.pdf) · [Figure 3 and legend](figure3_with_legend.pdf) · [reader](index.html) · [editable outline](../outline_v8.json) · [main plates 2–4](figures_2_4.pdf).

Figure 3 places eight native Wave-pinning profiles at C = 30 on the left, with time to 50% heading (B) above mean angular speed (C) on the right. B and C retain all four models and three sizes; i–ii identify the displayed profile placements. This revision reuses the four previously completed relaxation protocols and all existing readouts; no new simulations were run. The C = 47.6 profiles and their 20,000-horizon qualification remain in the previous main and legend appended to the supplement. [Complete supplement]({relfig}/figure3_round3_supplement.pdf) retains every round-2 main and supplement page, plus 90/99% heading-time curves. [Native source data and sampling qualifications]({relfig}/README.md).

Figures 1, 2, 4, 5 and 6 preserve their exact outline_v7 sections and artwork; Figure 2 retains its completed round-3 edition. outline_v7 and all preceding files remain intact. No timed refresh. Data generation: {generated}; producer complete.

Rebuild with the existing scaling_sandbox environment:

```sh
python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_single_C_20261022.py
python -B manuscript_outline/prepare_figure3_round3_single_C_20261022.py
python -B manuscript_outline/build_figure2_round3_outline_20261022.py --input manuscript_outline/outline_v8.json --output manuscript_outline/outline_v8.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v8.json --report manuscript_outline/outline_v8.validation.json
python -B manuscript_outline/verify_figure3_round3_single_C_20261022.py
```
''')
    prepend(ROOT/'manuscript_outline/README.md','''# Figure 3 single-C layout — current manuscript outline

[**Revised outline_v8**](outline_v8.pdf) · [**Figure 3 and legend**](revision_round3_figure3_single_C_20261022/figure3_with_legend.pdf) · [**Reader and complete supplement**](revision_round3_figure3_single_C_20261022/index.html).

Figure 3 now places the eight C = 30 activity profiles on the left and B above C on the right. Only i–ii are highlighted in the main; all three sizes remain in the response curves. The previous four-placement main and legend retain C = 47.6 in the supplement. No new simulations were run. Every round-2 main and supplement page is retained in the supplement. Figure 2 keeps its completed round-3 edition; all other main sections and artwork are identical to outline_v7.

[Methods, data and rebuild instructions](revision_round3_figure3_single_C_20261022/README.md).''')
    prepend(ROOT/'current manuscript figs/README.md','''# Figure 3 single-C layout — current figure and outline

[**Revised manuscript outline**](../manuscript_outline/outline_v8.pdf) · [**Figure 3 and legend**](../manuscript_outline/revision_round3_figure3_single_C_20261022/figure3_with_legend.pdf) · [**Reader, source data and full supplement**](../manuscript_outline/revision_round3_figure3_single_C_20261022/index.html).

Cue placement relative to the incumbent’s edge switches reorientation between sliding and takeover by a new region; enlargement slows both. Figure 3 now shows only C = 30 activity profiles on the left, with B above C on the right; both curve rows retain all models and sizes. The C = 47.6 profiles and their cue-off qualification remain in the supplement. All native results are reused; no new simulations were run for this layout revision. Round 2 and 90/99%-heading times are supplementary. Figure 2 retains round 3; Figures 1, 4, 5 and 6 retain their exact preceding sections and artwork.

Length represents WGD-associated enlargement. Observed migration summaries, conditional OU persistence and simulated polarity reorientation remain distinct readouts.''')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
