"""Print-width Figure2 candidate, reconstructed from compact CSV inputs only."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
F2=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2'
OLD=F2/'compensation_candidate_20261019'
CARD=F2/'debelly_width_timecourse_card_20261019'
OUT=F2/'compensation_print_candidate_20261019'
LABELS={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','debelly':'deBelly','holmes_model3':'Holmes 3','spring':'Spring'}
COLORS=dict(zip(LABELS,['#0072b2','#d55e00','#009e73','#7c5ba7','#a6783d']))
SOURCES={}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return list(csv.DictReader(p.open()))
def csvwrite(p,rr):
    keys=list(dict.fromkeys(k for r in rr for k in r))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rr)


def draw(org,width,native):
    plt.rcParams.update({'font.size':8.5,'axes.titlesize':9,'axes.labelsize':8.5,'xtick.labelsize':8,'ytick.labelsize':8,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(7.2,8.6));xx=[.105,.425,.745];ww=.225
    for i,(key,title) in enumerate([('excess_amount_over_reference','A  Localized active signal'),('width_over_reference','B  Incumbent extent'),('peak_over_reference','C  Peak concentration')]):
        ax=fig.add_axes([xx[i],.695,ww,.18])
        for model in LABELS:
            rr=sorted([r for r in org if r['model']==model],key=lambda r:float(r['size_native']))
            ax.plot([float(r['size_over_reference']) for r in rr],[float(r[key]) for r in rr],'-o',ms=2,lw=1.35,color=COLORS[model],label=LABELS[model])
        ax.set(title=title,xlabel='Length / reference',ylabel='Value / reference',xlim=(.1,2.7));ax.grid(alpha=.13)
    fig.legend(*ax.get_legend_handles_labels(),loc='upper center',ncol=5,bbox_to_anchor=(.52,.946),frameon=False,handlelength=1.5,columnspacing=1.,fontsize=8)
    d=fig.add_axes([.105,.355,.865,.235]);d.axhline(1,color='#888',ls='--',lw=.8)
    for model in LABELS:
        rr=sorted([r for r in width if r['model']==model and r['roi']=='organization_core' and r['background_rule']=='current_minimum' and float(r['loss_fraction'])==.5 and r['initial_10pc_core_overlap']=='False'],key=lambda r:float(r['length']))
        finite=[r for r in rr if r['ratio_midpoint']]
        d.plot([float(r['length_over_reference']) for r in finite],[float(r['ratio_midpoint']) for r in finite],'-o',color=COLORS[model],ms=3,lw=1.5)
        for r in finite:d.vlines(float(r['length_over_reference']),float(r['ratio_lower']),float(r['ratio_upper']),color=COLORS[model],lw=.8)
        for r in rr:
            if r['ratio_status']=='narrow_no_sampled_entry_through_horizon':
                x,y=float(r['length_over_reference']),float(r['ratio_upper'])
                d.plot([x-.015,x+.015],[y,y],color=COLORS[model],lw=1)
                d.annotate('',xy=(x,y*.72),xytext=(x,y),arrowprops=dict(arrowstyle='-|>',color=COLORS[model],lw=1.,mutation_scale=7))
            elif r['ratio_status']=='broad_no_sampled_entry_through_horizon':
                x,y=float(r['length_over_reference']),float(r['ratio_lower'])
                d.annotate('',xy=(x,y*1.38),xytext=(x,y),arrowprops=dict(arrowstyle='-|>',color=COLORS[model],lw=1.,mutation_scale=7))
    d.set(title='D  Broadening changes weakening time differently across size',xlabel='Length / reference',ylabel='Broader / narrower\nweakening time',yscale='log',xlim=(.2,2.7),ylim=(.007,2.1))
    d.grid(alpha=.13);d.text(.98,.985,'Ratio 1: equal weakening times',ha='right',va='top',transform=d.transAxes,fontsize=8)
    fig.text(.53,.284,'Corrected deBelly L39: same requested peak 3.6; cue remains ON',ha='center',fontsize=9)
    features=[('challenged_edge_rac','E  Challenged-edge Rac','Native concentration'),('right_endpoint_displacement','F  Right endpoint','Displacement (µm)'),('organization_core_current_minimum_excess_fraction','G  Localized incumbent','Excess signal / initial')]
    for i,(key,title,ylabel) in enumerate(features):
        ax=fig.add_axes([xx[i],.075,ww,.145])
        for w,style in [(1.5,'-'),(5.,'--')]:
            rr=sorted([r for r in native if float(r['width'])==w],key=lambda r:float(r['native_time']))
            ax.plot([float(r['native_time']) for r in rr],[float(r[key]) for r in rr],style,color='#333',lw=1.5,label=f'Nominal width {w:g}')
        ax.set(title=title,xlabel='Cue-on time (s)',ylabel=ylabel,xlim=(0,800));ax.set_xticks([0,400,800]);ax.grid(alpha=.13)
        if i==0:ax.set_ylim(0,4.8)
        if i==2:ax.axhline(.5,color='#888',ls=':',lw=.8);ax.set_ylim(-.03,1.03)
    fig.legend(*ax.get_legend_handles_labels(),loc='upper center',ncol=2,bbox_to_anchor=(.53,.267),frameon=False,handlelength=2.5,fontsize=8)
    fig.text(.53,.985,'Incumbent organization and cue compensation',ha='center',va='top',fontsize=14,weight='bold')
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--refresh-inputs',action='store_true');args=parser.parse_args()
    for n in ['supporting/data','supporting/vectors','supporting/metadata']:(OUT/n).mkdir(parents=True,exist_ok=True)
    orgfile=OUT/'supporting/data/organization_all97.csv';widthfile=OUT/'supporting/data/width_response_all_conventions.csv';nativefile=OUT/'supporting/data/debelly_L39_common_ON_window.csv'
    if args.refresh_inputs or not orgfile.exists():
        org=rows(OLD/'supporting/data/initial_organization_all97.csv');width=rows(OLD/'supporting/data/width_response_all_conventions.csv')
        native=[r for r in rows(CARD/'supporting/data/native_timecourses.csv') if float(r['length'])==39 and r['phase']=='cue_on' and float(r['native_time'])<=800]
        for p,rr in [(orgfile,org),(widthfile,width),(nativefile,native)]:csvwrite(p,rr)
        (OUT/'supporting/metadata/evidence_input_lineage.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    else:org=rows(orgfile);width=rows(widthfile);native=rows(nativefile)
    assert len(org)==97 and len(width)==1746 and len(native)==962
    fig=draw(org,width,native);fig.savefig(OUT/'figure2_print_candidate.png',dpi=240)
    for ext in ['pdf','svg']:fig.savefig(OUT/'supporting/vectors'/('figure2_print_candidate.'+ext))
    plt.close(fig)
    caption='''# Figure2 candidate: incumbent organization and cue compensation

A–C retain all97 exact qualified basals from prior work. The localized amount is active signal above the basal profile minimum integrated to the interpolated half-height boundary, not calibrated molecule number. Four models enlarge mainly through extent, while Goryachev grows taller. Amount, width and peak are normalized to each model's own reference (WP/Goryachev/deBelly L19,Holmes L20,Spring restingL0=1). Spring widths/amounts use its physical domain; deBelly retains nominal solver coordinates. The13 ineligible preparations remain separate.

D is a new within-model native experiment. Fixed requested native peaks are0.1WP,6.4Goryachev,1Holmes,4Spring,3.6deBelly. Broader/narrower supports are5/1.5, exceptSpring0.4/0.1; sampled peaks and integrated inputs are exported. The plotted point is the ratio of saved-entry bracket midpoints for first50% loss of localized excess in the fixed initial organization core, using the current whole-profile minimum as background. Thin bounds use the bracket endpoints. Downward/upward arrows give one-sided bounds on this sampled-entry ratio when one input has no saved entry through its finite continuous-ON horizon (WP/Holmes10000,Goryachev5000,Spring50,deBelly2000). They do not denote infinite latency or an optimized pulse threshold. Of91 nonoverlap size pairs,75 have two observed primary entries,12 have one-sided bounds and4 have neither entry/no ratio;6 initial10%-core-overlap pairs are omitted from the main plot and remain in source tables.

The full sensitivity table preserves two initial-core conventions (fractional organization versus whole-cell worker), three background rules (current minimum,current native-sample5thpercentile,frozen basal minimum) and25/50/75% loss. Missing alternatives are explicitly unavailable, not assigned finite values. The supplement shows completeness and sensitivity. A current-profile floor is a descriptive concentration reference, not a molecular pool. Spring current-floor excess with evolving physical weights is invariant to pure uniform dilution, while frozen basal concentration subtraction is not.

E–G show actual native corrected independent-endpoint deBelly L39 fields for the same requested peak3.6 and nominal support1.5 or5, during a common0–800second cue-ON window. More challenged-edge Rac accompanies an early delay in right endpoint displacement; the wider case later catches up and exceeds the narrow displacement, while50% incumbent-excess loss is later: (370,372.5] versus(447.5,450]. At100s the wider case has edgeRac4.0137 versus2.7684, tension7.2824 versus3.0185 and endpoint displacement0.5087 versus1.1319µm. The native command contains the tension factor1/(1+tension²/10); at this state both Rac switches are saturated, so the conditional broad/narrow command factor is≈0.303. This conditional command differs from actual coupled membrane velocity and accumulated displacement. Full2000ON+10000OFF histories, native L18 contrast and actual physical cue mapping remain in support. Inputs stay fixed in nominal coordinates; this is not a calibrated beam-width prediction.

The targeted500second withdrawal makes both widths return, with native field-tailPASS. Both800/2000second baseline pulses pass retained scalar-strength and native field-tail criteria. Thus earlier50% regional loss does not itself establish a different retention fate at the tested500second pulse. No optimized minimum is inferred or pursued from this pair. Scalar gates remain separate from spatial organization and do not alone prove whole-field replacement.

The association between basal organization and compensation is not an isolated causal mediation estimate or unique molecular identification. Broader inputs change integrated dose; native channels and operating points differ across models. The expanded all-five/all-length strength and duration results, including corrected deBelly, are retained as comprehensive support. Native grids and solvers are unchanged; no new spatial-convergence, reporter-kinetics or common optical-actuator calibration is claimed.\n'''
    detailed=re.sub(r'(?<=[A-Za-z])(?=\d)|(?<=\d)(?=[A-Za-z])',' ',caption)
    (OUT/'DETAILED_METHODS.md').write_text(detailed)
    caption='''# Figure 2. Incumbent organization and cue compensation

**A–C, prior organization evidence.** All 97 qualified basal states contain more localized active signal at greater length, through different combinations of extent and concentration. Amount is signal above the basal minimum integrated to the interpolated half-height boundary. Amount, width and peak are divided by each model's own reference value: WP/Goryachev/deBelly L19, Holmes L20, Spring resting L₀ = 1. Spring uses physical length; deBelly retains nominal solver coordinates. Thirteen ineligible preparations remain separately documented.

**D, new broadening experiments.** Within each model, the same requested native peak is applied over narrower/broader supports 1.5/5, or 0.1/0.4 in Spring. Broader inputs also increase integrated input. “Weakening time” is first sampled 50% loss of positive active signal above the current whole-profile minimum, integrated in the fixed initial organization core. Points are ratios of saved-entry bracket midpoints; thin bounds use bracket endpoints. Below one, broadening accelerates this readout; above one, it delays it. Arrows denote one-sided sampled-entry ratio bounds when one width has no entry by its finite horizon; they do not denote infinite latency or an optimized pulse threshold. Among 91 nonoverlap size pairs, 75 have two observed entries, 12 have one-sided bounds and four have neither entry. Six initial cue/core-overlap pairs remain in the source table. All 97 no-cue controls retain localized excess over their matched horizons, including 25/75% sensitivity.

**E–G, a native coupled-response example.** Corrected independent-endpoint deBelly L39 receives the same requested peak 3.6 over nominal widths 1.5 or 5. The common 0–800 s window remains cue-ON. Wider stimulation produces more challenged-edge Rac, an early delay followed by catch-up in endpoint displacement, and later 50% incumbent-excess loss: (370,372.5] versus (447.5,450] s. These readouts are distinct. Both widths return after a tested 500 s pulse; both 800/2000 s pulses pass retained scalar-strength and native field-tail criteria. Earlier regional weakening therefore does not establish a different retention fate at 500 s. Scalar gates alone do not prove whole-field replacement.

The comprehensive supplement retains full strength/duration outcomes, complete histories, background/core sensitivity and native physical cue mappings. Current-profile floors describe localized signal, not a molecular pool. The example concerns native coupled activity and mechanics, not isolated causal mediation or a calibrated beam-width prediction. Original discretizations and finite observation windows remain explicit in the detailed methods.
'''
    (OUT/'CAPTION.md').write_text(caption)
    (OUT/'README.md').write_text('# Figure2 print candidate\n\n![7.2-inch Figure2 candidate](figure2_print_candidate.png)\n\n[PDF](supporting/vectors/figure2_print_candidate.pdf) · [Caption](CAPTION.md) · [Redraw closure](supporting/metadata/redraw_closure.json).\n\nA–C retain the prior starting-state organization; D tests length-dependent compensation during stimulation; E–G show a coupled native activity/motion example. Nominal input geometry and finite sampled windows remain explicit. Main composition remains provisional for author review.\n')
    compact={str(p.relative_to(ROOT)):sha(p) for p in [orgfile,widthfile,nativefile]};source={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}
    closure=dict(command='python polarity1d/manuscript/code/mechanisms/revision_compensation_print_candidate_20261019.py',working_directory='repository root',dependencies=['numpy','matplotlib'],native_archives_required_for_redraw=False,source_hashes=source,compact_inputs=compact)
    (OUT/'supporting/metadata/redraw_closure.json').write_text(json.dumps(closure,indent=2)+'\n')
    repo=ROOT.parents[1];rel=lambda p:str(p.relative_to(repo))
    content=dict(title='Figure2. Incumbent organization and cue compensation',role='Larger incumbents arise through different spatial organizations and show different native compensation responses.',status='Provisional author-review candidate; numerical/display independent review pending.',
        assets=[dict(path=rel(OUT/'figure2_print_candidate.png'),caption='Incumbent organization, length-dependent broadening response and native deBelly example',page_size=[7.2,8.6],page_title='Figure2. Organization and compensation',minimum_printed_width_inches=7.2)],
        legend=[caption],narrative=['Localized incumbent growth combines extent and concentration differently.','Broadening changes weakening time differently across length, with a corrected deBelly exception.','The native deBelly example separates challenger activity, coupled boundary response, incumbent organization loss and retained fate.'],
        measurements=['Fixed-initial-core localized excess signal and native activity, separately from raw regional activity.','Ratios of saved-entry bracket midpoints, with finite sampled-entry bounds.','Native challenged-edge Rac and actual endpoint displacement during a common cue-ON window.'],
        limitations=['Native/nominal input geometry and integrated dose change are explicit; no calibrated common optical actuation.','Descriptive profile floors and regional redistribution are not molecular depletion.','No isolated causal mediation or unique mechanism classification.','Full finite withdrawal outcomes and field-tail flags remain separate from transient weakening.'],
        sources=[dict(label='Compact reconstruction closure',path=rel(OUT/'supporting/metadata/redraw_closure.json')),dict(label='Detailed methods and native conditional law',path=rel(OUT/'DETAILED_METHODS.md')),dict(label='No-cue regional controls',path=rel(F2/'five_model_compensation_20261019/selection/supporting/data/no_cue_localized_excess_loss_check.csv')),dict(label='All width sensitivities',path=rel(widthfile)),dict(label='Full native deBelly timecourse/support',path=rel(CARD/'README.md')),dict(label='All baseline native outcomes',path=rel(F2/'five_model_compensation_20261019/baseline/README.md')),dict(label='All strength native outcomes',path=rel(F2/'five_model_compensation_20261019/strength/README.md'))])
    (OUT/'supporting/metadata/content.json').write_text(json.dumps(content,indent=2)+'\n')
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),source_hashes={**source,**compact},
        artifact_hashes={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='validation.json'},
        native_replay='Independent full-stage proofs retained separately; this reconstruction uses compact numerical inputs only.'),indent=2)+'\n')
    print(str(OUT/'figure2_print_candidate.png'))


if __name__=='__main__':main()
