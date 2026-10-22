"""Build an explicitly in-progress review; final outline_v5 awaits Holmes runs."""
from copy import deepcopy
from datetime import datetime, timezone
import json, os
from pathlib import Path
from pypdf import PdfWriter
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
OUTLINE=ROOT/'manuscript_outline'; DEST=OUTLINE/'revision_round2_20261022'
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022'
def read(p):return json.loads(p.read_text())
def rel(p):return str(p.relative_to(ROOT))
def write(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n')

def main():
    baseline=read(OUTLINE/'outline_v4.json');doc=deepcopy(baseline)
    fig2=read(FIG/'figure2/legend.json');fig3=read(FIG/'figure3/content.json');four=read(FIG/'figure4/supporting/metadata/main_summary.json')
    narrative=('Enlarged cells resist imposed replacement through distinct initiation, spatial and temporal response signatures (Figure 2); the proposed class interpretation still awaits the specified Holmes large-size retention measurements. '
        'At matched distance, size dominates driven takeover, while enlarged cells reach larger separations where cue-free resolution in mass-conserved models becomes very slow or unresolved (Figure 3). '
        +four['message']+' (Figure 4). Weak-input geometric comparisons describe point estimates; Holmes and Spring intervals also admit sub-geometric gains. '
        'Observed migration summaries, conditional OU velocity-memory estimates and simulated polarity-heading readouts remain separate evidence.')
    doc.update(date='22 October 2026',stage='Round-2 follow-up IN PROGRESS: Figure 2 awaits Holmes SLURM runs',subtitle='Review candidate; final outline_v5 integration is deferred until Figures 2–4 are complete')
    doc['overview'][1]=narrative
    doc['overview'].append('This in-progress review is separate from outline_v5. Figures 3 and 4 are recomposed; Figure 2 has explicit missing Holmes measurements. No class membership changes or optional Figure 3 horizon extensions were made. Scheduler socket permissions prevented the prepared Holmes campaign from being submitted.')
    for s in doc['sections']:
        if s['id']=='figure2':
            s.update(fig2);s['title']='Figure 2. '+fig2['title']+' — IN PROGRESS'
            s['main_pdf']=rel(FIG/'figure2/figure2_round2_pending.pdf');s['supplement_pdf']=rel(FIG/'figure2/figure2_round2_supplement.pdf')
            s['role']='Recomposed mechanism-class candidate; specified Holmes large-size retention measurements are pending.'
            s['narrative']=['Existing response signatures motivate three proposed classes, with unresolved cells clearly marked.','No continuous commitment by 100,000 at Holmes L31, L35 and L39 demonstrates extreme finite-horizon replacement resistance.','Large-size Holmes pulse and noise-retention results are unavailable: the frozen campaign could not contact SLURM.']
        elif s['id']=='figure3':
            s.update(fig3);s['narrative']=[fig3['legend'][0],fig3['legend'][5]]
            s['sources'].append(dict(label='Round-2 decomposition and native-source audit',path=rel(FIG/'figure3/audit.json')))
        elif s['id']=='figure4':
            s.update(title='Figure 4. '+four['message'],status='round2_Case_A_posthoc',role='Decay-rate size advantage and descriptive spatial organization under fluctuating input.',
                assets=[dict(path=rel(FIG/'figure4/figure4_main.png'),caption=four['message'],page_title='Figure 4. Enlarged cells hold heading longer under fluctuating input',page_size=[16.5,21],minimum_printed_width_inches=14)],
                main_pdf=rel(FIG/'figure4/figure4_main.pdf'),supplement_pdf=rel(FIG/'figure4/figure4_supplement.pdf'),legend=[four['legend']],
                narrative=[four['message'],'Post hoc H1 passes in both main-view arms; the frozen rate H1 remains inconclusive. The H1 domain uses inputs below the point estimate of the largest size’s own noise tolerance.','Geometry and multiple-front fractions are descriptive; erosion below geometry also occurs in single-front conditions.'],
                limitations=['The H1 domain amendment is post hoc.','Weak-input geometric comparisons describe point estimates; not all model CIs exclude sub-geometric values.','Non-exponential lag dependence and forcing qualifications remain explicit. Goryachev is retained only in the supplement.','H2, H3 and A6 tau gates are retired; frozen historical calls remain for provenance.','Simulated polarity heading is distinct from migration and OU persistence.'])
            s['sources'].extend([dict(label='Amended and frozen H1 calls',path=rel(FIG/'figure4/supporting/data/S0_posthoc_H1.csv')),dict(label='Main source and render receipt',path=rel(FIG/'figure4/supporting/metadata/main_render_receipt.json'))])
    for ident in ('figure1','figure5','figure6'):
        assert next(s for s in doc['sections'] if s['id']==ident)==next(s for s in baseline['sections'] if s['id']==ident)
    for s in doc['sections']:
        if s['id'] not in ('figure2','figure3','figure4'):continue
        for a in s['assets']:
            w,h=Image.open(ROOT/a['path']).size
            width=max(a['page_size'][0],a['minimum_printed_width_inches']+1.5)
            a['page_size']=[width,max(a['page_size'][1],3.5+(width-1.4)*h/w)]
    write(DEST/'outline_v5_in_progress.json',doc)
    writer=PdfWriter()
    for p in (FIG/'figure2/figure2_round2_pending.pdf',FIG/'figure3/figure3_main.pdf',FIG/'figure4/figure4_main.pdf'):writer.append(p)
    with (DEST/'figures_2_4_in_progress.pdf').open('wb') as f:writer.write(f)
    write(DEST/'followup_integration.json',dict(generated_utc=datetime.now(timezone.utc).isoformat(),status='in_progress_pending_Holmes',final_outline_v5_integrated=False,protected_sections=['figure1','figure5','figure6'],cross_figure_narrative=narrative,figure4_case='A',figure4_H1_posthoc=True,optional_figure3_extensions=9,optional_extensions_submitted=False))
    dest=DEST/'followup.html'
    link=lambda p:os.path.relpath(p,dest.parent)
    body=['<!doctype html><meta charset="utf-8"><title>Round-2 follow-up: in progress</title><style>body{max-width:1450px;margin:2rem auto;font:17px/1.5 system-ui}img{width:100%}</style>',
        '<h1>Figures 2–4 round-2 follow-up — in progress</h1>',
        '<p>Generated '+datetime.now(timezone.utc).isoformat()+'. Producers finished; Holmes simulations are unsubmitted.</p>',
        '<p><a href="outline_v5_in_progress.pdf">In-progress outline review</a> · <a href="figures_2_4_in_progress.pdf">Three main plates</a></p>',
        '<p>'+narrative+'</p>',
        '<p>Final outline_v5 integration awaits the specified Holmes measurements. Nine optional Figure 3 Holmes extensions await author review; none were submitted. Class membership is unchanged.</p>']
    for number,stem in [(2,'figure2_round2_pending'),(3,'figure3_main'),(4,'figure4_main')]:
        folder=FIG/f'figure{number}';body+=['<h2>Figure '+str(number)+'</h2><a href="'+link(folder/'README.md')+'">README and qualifications</a><img src="'+link(folder/(stem+'.png'))+'">']
    dest.write_text('\n'.join(body))
    (DEST/'FOLLOWUP_README.md').write_text('# Round-2 follow-up — in progress\n\n[Review reader](followup.html) · [In-progress outline](outline_v5_in_progress.pdf) · [Main plates](figures_2_4_in_progress.pdf)\n\nFigure 4 Case A and Figure 3 recomposition are complete. Figure 2 includes a frozen and verified 32-worker Holmes campaign, but scheduler socket permissions prevented submission; no new trajectories exist. Class membership remains unchanged. Nine Figure 3 optional extensions are enumerated and unsubmitted. Final outline_v5 integration remains deferred until Figure 2 is complete; the explicitly in-progress review preserves Figures 1, 5 and 6 exactly.\n\nFigure 4 H1 is explicitly post hoc and paired with frozen calls. H2/H3/A6 gates are retired; geometry and multiple-front fractions remain descriptive. All original supplements are preserved.\n')
    print(DEST/'outline_v5_in_progress.json')

if __name__=='__main__':main()
