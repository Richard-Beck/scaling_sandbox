"""Comprehensive tested-grid supplement; compact CSV-only reconstruction."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parents[2]
F2=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2'
OUT=F2/'five_model_compensation_20261019/supplement'
MODELS=['wave_pinning','goryachev','holmes_model3','spring','debelly']
LABELS=dict(zip(MODELS,['Wave-pinning','Goryachev','Holmes 3','Spring','Corrected linear deBelly']))
SOURCES={}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def rows(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return list(csv.DictReader(p.open()))
def csvwrite(p,rr):
    keys=list(dict.fromkeys(k for r in rr for k in r))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rr)


def extract():
    outcomes=[];events=[]
    for stage,native_name in [('baseline','revision_five_model_compensation_baseline_v2_20261019'),('strength','revision_five_model_compensation_strength_v2_20261019')]:
        native=ROOT/'data/runs'/native_name;nm=read(native/'manifest.json');cm={c['slot']:c for c in nm['cases']}
        deriv=ROOT/'data/runs'/('revision_compensation_excess_'+stage+'_20261019');dm=read(deriv/'manifest.json')
        assert dm['native_manifest_sha256']==sha(native/'manifest.json')
        lookup={}
        for model in MODELS:
            folder=deriv/'workers'/model;result=read(folder/'result.json')
            for rr in rows(folder/'pulse_endpoints.csv'):
                if rr['phase']=='after_withdrawal':lookup[(int(rr['case_slot']),float(rr['pulse_duration']))]=rr
            path=folder/'sampled_loss_events.csv';assert sha(path)==result['derived_artifact_hashes'][str(path.relative_to(deriv))]
            SOURCES[str(path.relative_to(ROOT))]=sha(path)
            for e in csv.DictReader(path.open()):
                c=cm[int(e['case_slot'])]
                if c['case_kind']=='constant_on' and e['phase']=='cue_on' and float(e['pulse_duration'])==max(c['pulse_durations']) and e['observable']=='excess_fraction' and float(e['loss_fraction']) in [.25,.5,.75]:
                    events.append(dict(stage=stage,case_slot=c['slot'],model=model,length=c['length'],width=c['width'],native_peak_multiplier=c['native_peak_multiplier'],
                        requested_native_peak=c['response_amplitude'],sampled_native_peak=c['initial_geometry']['sampled_peak'],
                        initial_10pc_core_overlap_fraction=c['initial_geometry']['initial_core_10pc_footprint_overlap_fraction'],
                        roi=e['roi'],background_rule=e['background_rule'],loss_fraction=float(e['loss_fraction']),
                        observation_horizon=float(e['observation_end']),entry_status=e['status'],entry_lower=float(e['lower']),entry_upper=float(e['upper']) if e['upper'] else None))
        pulsefile=F2/'five_model_compensation_20261019'/stage/'current/data/pulse_outcomes.csv'
        source=rows(pulsefile);assert len(source)==nm['logical_protocol_rows']
        for r in source:
            assert r['state']=='complete';c=cm[int(r['case_slot'])];d=lookup[(int(r['case_slot']),float(r['pulse_duration']))]
            outcomes.append(dict(stage=stage,case_slot=c['slot'],model=c['model'],length=c['length'],width=c['width'],native_peak_multiplier=c['native_peak_multiplier'],
                pulse_duration=float(r['pulse_duration']),requested_native_peak=c['response_amplitude'],sampled_native_peak=c['initial_geometry']['sampled_peak'],
                initial_10pc_core_overlap_fraction=c['initial_geometry']['initial_core_10pc_footprint_overlap_fraction'],
                scalar_orientation_success=r['original_scalar_orientation_success'],scalar_full_strength_success=r['original_scalar_full_polarity_success'],
                native_field_tail_pass=r['native_field_tail_pass'],native_field_tail_max_scaled=r['native_field_tail_max_scaled'],
                no_cue_preparation_valid=r['no_cue_preparation_valid'],raw_initial_worker_core_fraction_final=r['original_core_active_fraction_final'],
                localized_organization_excess_fraction_final=d['organization_core_current_minimum_excess_fraction_last'],
                localized_quantile_excess_fraction_final=d['organization_core_current_sample_quantile_0p05_excess_fraction_last'],
                frozen_floor_excess_fraction_final=d['organization_core_frozen_initial_minimum_excess_fraction_last'],
                actual_off_archive=r['actual_off_archive'],actual_off_archive_sha256=r['actual_off_archive_sha256']))
    assert len(outcomes)==3492 and len(events)==12222
    return outcomes,events


def good(value):return value is True or value=='True'


def page(model,outcomes,events):
    rr=[r for r in outcomes if r['model']==model];lengths=sorted({float(r['length']) for r in rr});mid=.2 if model=='spring' else 3.
    ee=[e for e in events if e['model']==model and e['roi']=='organization_core' and e['background_rule']=='current_minimum' and float(e['loss_fraction'])==.5]
    plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(11.5,10.3));gs=fig.add_gridspec(2,2,height_ratios=[1,1.7],left=.12,right=.975,bottom=.20,top=.87,hspace=.40,wspace=.27)
    palettes=[['#505050','#0072b2','#d55e00'],['#bbbbbb','#8cb7d2','#0072b2','#d55e00','#7b438e']]
    for col,kind in enumerate(['width','strength']):
        ax=fig.add_subplot(gs[0,col]);levels=sorted({float(e['width']) for e in ee if float(e['native_peak_multiplier'])==1}) if kind=='width' else [.03,.1,.3,1.,3.]
        for i,level in enumerate(levels):
            group=[e for e in ee if (float(e['native_peak_multiplier'])==1 and float(e['width'])==level) if kind=='width'] if kind=='width' else [e for e in ee if float(e['width'])==mid and float(e['native_peak_multiplier'])==level]
            group.sort(key=lambda e:float(e['length']));color=palettes[col][i];xx=[float(e['length']) for e in group]
            yy=[(float(e['entry_lower'])+float(e['entry_upper']))/2 if e['entry_upper'] not in ['',None] else np.nan for e in group]
            missing=sum(e['entry_upper'] in ['',None] for e in group)
            label=(f'w={level:g}' if kind=='width' else f'{level:g}×')+f'; no entry {missing}/{len(group)}'
            ax.plot(xx,yy,'-o',color=color,lw=1.4,ms=3,label=label)
            for e in group:
                x=float(e['length']);upper=e['entry_upper']
                if upper in ['',None]:ax.scatter(x,float(e['observation_horizon']),marker='^',s=25,facecolors='none',edgecolors=color)
                elif float(e['initial_10pc_core_overlap_fraction'])>0:ax.scatter(x,(float(e['entry_lower'])+float(upper))/2,marker='o',s=35,facecolors='white',edgecolors=color,zorder=4)
        ax.set(yscale='log',xlabel='Native resting length L₀' if model=='spring' else 'Native nominal length L',ylabel='Time to 50% localized excess loss',title='A  Width at baseline peak' if col==0 else f'B  Peak factor at width {mid:g}')
        ax.grid(alpha=.16)
        legend_options={'loc':'upper left','bbox_to_anchor':(0,.88)} if model=='spring' and col==1 else {}
        ax.legend(frameon=False,fontsize=8,ncol=2,**legend_options)
    colors=['#d8d8d8','#298d78','#eaa756'];cmap=ListedColormap(colors);norm=BoundaryNorm([-.5,.5,1.5,2.5],3)
    for col,kind in enumerate(['width','strength']):
        ax=fig.add_subplot(gs[1,col]);pulses=sorted({float(r['pulse_duration']) for r in rr if float(r['pulse_duration'])>0})
        levels=sorted({float(r['width']) for r in rr if float(r['native_peak_multiplier'])==1 and float(r['pulse_duration'])>0}) if kind=='width' else [.03,.1,.3,1.,3.]
        rowkeys=[(level,t) for level in levels for t in pulses];grid=np.zeros((len(rowkeys),len(lengths)));overlap=[]
        for i,(level,t) in enumerate(rowkeys):
            for j,L in enumerate(lengths):
                r=next(r for r in rr if float(r['length'])==L and float(r['pulse_duration'])==t and ((float(r['width'])==level and float(r['native_peak_multiplier'])==1) if kind=='width' else (float(r['width'])==mid and float(r['native_peak_multiplier'])==level)))
                grid[i,j]=2 if not good(r['native_field_tail_pass']) else 1 if good(r['scalar_full_strength_success']) else 0
                if float(r['initial_10pc_core_overlap_fraction'])>0:overlap.append((j,i))
        ax.imshow(grid,aspect='auto',interpolation='nearest',cmap=cmap,norm=norm,origin='upper')
        if overlap:ax.scatter(*zip(*overlap),marker='.',s=8,color='black')
        ax.set_xticks(range(len(lengths)));ax.set_xticklabels([f'{L:g}' for L in lengths],rotation=90,fontsize=8)
        ax.set_yticks(range(len(rowkeys)));ax.set_yticklabels([f'{"w" if kind=="width" else "×"}={level:g}; T={t:g}' for level,t in rowkeys],fontsize=8)
        ax.set(xlabel='Native qualified lengths; discrete tested conditions',title='C  Width × pulse outcomes' if col==0 else 'D  Peak factor × pulse outcomes')
        for k in range(len(levels)-1):ax.axhline((k+1)*len(pulses)-.5,color='white',lw=2)
    fig.text(.5,.975,LABELS[model]+': complete tested width, strength and pulse series',ha='center',fontsize=15,weight='bold')
    controls=[r for r in rr if float(r['pulse_duration'])==0]
    fig.text(.5,.938,f'{len(lengths)} exact qualified basals; {len(rr)} protocol rows including {len(controls)} valid no-cue controls; 10000 native units after withdrawal.',ha='center',fontsize=10)
    fig.text(.5,.907,'Observed times are saved-entry bracket midpoints; exact brackets and finite no-entry horizons remain in the compact tables.',ha='center',fontsize=9)
    fig.legend(handles=[Patch(color=colors[0],label='Scalar full strength: fail; field tail: pass'),Patch(color=colors[1],label='Scalar full strength + field tail: pass'),Patch(color=colors[2],label='Field tail: fail (scalar flag separately retained)')],loc='upper center',ncol=1,bbox_to_anchor=(.54,.13),frameon=False,fontsize=9)
    fig.text(.09,.047,'Horizon triangles may coincide: curve legends give no-entry counts. Open circles / map dots mark initial cue/core overlap.',fontsize=9)
    fig.text(.09,.027,'Baseline 1× / middle width is displayed in both contexts, but is one protocol identity. Scalar gates do not alone establish spatial replacement.',fontsize=9)
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--extract',action='store_true');args=parser.parse_args()
    for n in ['supporting/data','supporting/vectors','supporting/metadata']:(OUT/n).mkdir(parents=True,exist_ok=True)
    pp=OUT/'supporting/data/all3492_protocol_outcomes.csv';ep=OUT/'supporting/data/all12222_ON_loss_conventions.csv'
    if args.extract:
        outcomes,events=extract();csvwrite(pp,outcomes);csvwrite(ep,events)
        (OUT/'supporting/metadata/evidence_source_hashes.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    else:outcomes=rows(pp);events=rows(ep)
    with PdfPages(OUT/'supporting/vectors/complete_five_model_compensation.pdf') as pdf:
        for model in MODELS:
            fig=page(model,outcomes,events);fig.savefig(OUT/(model+'_response_summary.png'),dpi=170);fig.savefig(OUT/'supporting/vectors'/(model+'_response_summary.svg'));pdf.savefig(fig);plt.close(fig)
    (OUT/'README.md').write_text('# Complete native cue experiment: comprehensive support\n\n[Five-page all-model PDF](supporting/vectors/complete_five_model_compensation.pdf).\n\nEvery one of the 3492 requested baseline/strength protocols is retained in the [compact outcome table](supporting/data/all3492_protocol_outcomes.csv), including 97 no-cue controls. Each model page shows every qualified length, all three baseline widths, all five middle-width peak factors and five coarse pulse durations. Upper observed time points are saved-entry bracket midpoints for conditional first sampled 50% loss of localized incumbent excess during continuous stimulation; exact brackets remain in CSV. No-entry triangles may coincide at the same finite horizon, so every upper-curve legend explicitly reports its no-entry count. Bottom rows retain the distinct scalar-strength and native field-tail flags after 10000 cue-free native units. A lack of sampled success is not an optimized threshold bound, and no duration monotonicity or absence of untested success islands is assumed.\n\n[All 12222 ON loss/core/background records](supporting/data/all12222_ON_loss_conventions.csv) retain 25/50/75% loss, three floor rules and both core conventions, including missing entries and overlap. Inputs and clocks are native/model-specific; larger width changes integrated input. This supplement is a complete experimental record, not a molecular classification or common optical-actuator calibration.\n\n[Two additional targeted deBelly T500 outcomes](../selection/supporting/data/debelly_targeted_T500_outcomes.csv) are retained separately from the core 3492 protocols, for 3494 logical protocols in the expanded study. [Selection/control/source notes](../selection/README.md) document both returns and passing native field tails without inferring a width-specific retained fate or optimized threshold.\n')
    closure=dict(command='python polarity1d/manuscript/code/mechanisms/revision_compensation_supplement_20261019.py',working_directory='repository root',native_archives_required_for_redraw=False,
        source_hashes={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},compact_inputs={str(p.relative_to(ROOT)):sha(p) for p in [pp,ep]})
    (OUT/'supporting/metadata/redraw_closure.json').write_text(json.dumps(closure,indent=2)+'\n')
    repo=ROOT.parents[1];rel=lambda p:str(p.relative_to(repo))
    content=dict(title='Comprehensive linear compensation experiments',role='Full five-model/97-length performed width, peak and pulse campaign.',
        status='Native and derived arrays independently validated; supplementary display review pending.',
        assets=[dict(path=rel(OUT/(model+'_response_summary.png')),caption=LABELS[model]+' complete tested compensation outcomes',page_size=[11.5,10.3],page_title=LABELS[model]+': compensation series',minimum_printed_width_inches=10.) for model in MODELS],
        legend=['Each page shows baseline-width and middle-width peak-factor loss clocks, followed by every discrete tested pulse outcome. Upper observed time points are saved-entry bracket midpoints; exact brackets remain in CSV. Horizon triangles may coincide, with no-entry counts explicitly reported per curve. Scalar full-strength and native field-tail flags remain separate. Open event markers and map dots identify initial cue/core overlap. No sampled entry at a finite horizon is not an infinite latency or optimized duration threshold.'],
        narrative=['All 3492 requested protocols and 12222 loss/core/background convention records remain accessible; the baseline 1× middle-width condition is shown in both contexts with one source identity.'],
        measurements=['First sampled 50% localized-excess loss during continuous cue, with original brackets retained in CSV.','Retained scalar orientation/full-strength, field-tail stability and final native regional/excess quantities after 10000 cue-free native units.'],
        limitations=['Discrete coarse duration outcomes do not establish monotonicity or exclude untested success islands.','Native actuators, integrated inputs and clocks differ across models.','Regional excess is a descriptive organization measure; scalar reversal alone is not spatial proof.'],
        sources=[dict(label='Full all-model PDF',path=rel(OUT/'supporting/vectors/complete_five_model_compensation.pdf')),dict(label='3492 actual protocol outcomes',path=rel(pp)),dict(label='12222 loss conventions',path=rel(ep)),dict(label='CSV-only reconstruction closure',path=rel(OUT/'supporting/metadata/redraw_closure.json'))])
    (OUT/'supporting/metadata/content.json').write_text(json.dumps(content,indent=2)+'\n')
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),protocol_rows=len(outcomes),loss_convention_rows=len(events),
        source_hashes=SOURCES,artifact_hashes={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='validation.json'}),indent=2)+'\n')
    print(str(OUT/'supporting/vectors/complete_five_model_compensation.pdf'))


if __name__=='__main__':main()
