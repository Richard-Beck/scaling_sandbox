"""Compact provisional organization + finite-response candidate, no integration."""
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

ROOT=Path(__file__).resolve().parents[2]
F2=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2'
OUT=F2/'compensation_candidate_20261019'
RUN=ROOT/'data/runs/revision_compensation_excess_baseline_20261019'
ORG=F2/'finite_response_candidate_20261019/supporting/data/organization_all97.csv'
LABELS={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','debelly':'deBelly','holmes_model3':'Holmes 3','spring':'Spring'}
COLORS=dict(zip(LABELS,['#0072b2','#d55e00','#009e73','#7c5ba7','#a6783d']))
REFS=dict(zip(LABELS,[19.,19.,19.,20.,1.]))
FLOORS=['current_minimum','current_sample_quantile_0p05','frozen_initial_minimum']
SOURCES={}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def rows(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return list(csv.DictReader(p.open()))
def csvwrite(p,rr):
    keys=list(dict.fromkeys(k for r in rr for k in r))
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rr)


def comparison(a,b):
    al,bl=float(a['lower']),float(b['lower']);au=float(a['upper']) if a['upper'] else None;bu=float(b['upper']) if b['upper'] else None
    if au is not None and bu is not None and al>0:
        return dict(ratio_status='two_observed_sampled_entries',ratio_lower=bl/au,ratio_upper=bu/al,ratio_midpoint=(bl+bu)/(al+au),
            sign='broader_faster' if bu<al else 'broader_slower' if au<bl else 'overlapping_sample_brackets')
    if bu is not None and au is None:
        return dict(ratio_status='narrow_no_sampled_entry_through_horizon',ratio_lower=None,ratio_upper=bu/al,ratio_midpoint=None,sign='broader_observed_narrow_not_observed')
    if au is not None and bu is None:
        return dict(ratio_status='broad_no_sampled_entry_through_horizon',ratio_lower=bl/au,ratio_upper=None,ratio_midpoint=None,sign='narrow_observed_broader_not_observed')
    return dict(ratio_status='neither_sampled_entry_observed',ratio_lower=None,ratio_upper=None,ratio_midpoint=None,sign='unresolved_finite_window')


def collect():
    m=read(RUN/'manifest.json');lookup={};cases={r['case']['slot']:r['case'] for r in m['cases']};data=[]
    for model in m['models']:
        path=RUN/'workers'/model/'sampled_loss_events.csv';result=read(path.with_name('result.json'))
        assert result['manifest_sha256']==sha(RUN/'manifest.json')
        assert result['derived_artifact_hashes'][str(path.relative_to(RUN))]==sha(path)
        SOURCES[str(path.relative_to(ROOT))]=sha(path)
        for e in csv.DictReader(path.open()):
            c=cases[int(e['case_slot'])]
            if c['case_kind']=='constant_on' and e['phase']=='cue_on' and float(e['pulse_duration'])==max(c['pulse_durations']) and e['observable']=='excess_fraction' and float(e['loss_fraction']) in [.25,.5,.75]:
                lookup[(model,c['length'],c['width'],e['roi'],e['background_rule'],float(e['loss_fraction']))]=e
        widths=[.1,.4] if model=='spring' else [1.5,5.]
        for length in sorted({c['length'] for c in cases.values() if c['model']==model}):
            narrow=next(c for c in cases.values() if c['model']==model and c['length']==length and c['width']==widths[0] and c['case_kind']=='constant_on')
            broad=next(c for c in cases.values() if c['model']==model and c['length']==length and c['width']==widths[1] and c['case_kind']=='constant_on')
            overlap=any(c['initial_geometry']['initial_core_10pc_footprint_overlap_fraction']>0 for c in [narrow,broad])
            for roi in ['organization_core','worker_core']:
                for floor in FLOORS:
                    for loss in [.25,.5,.75]:
                        a=lookup[(model,length,widths[0],roi,floor,loss)];b=lookup[(model,length,widths[1],roi,floor,loss)]
                        data.append(dict(model=model,length=length,length_over_reference=length/REFS[model],narrow_width=widths[0],broad_width=widths[1],
                            requested_native_peak=narrow['response_amplitude'],narrow_sampled_peak=narrow['initial_geometry']['sampled_peak'],
                            broad_sampled_peak=broad['initial_geometry']['sampled_peak'],initial_10pc_core_overlap=overlap,roi=roi,background_rule=floor,loss_fraction=loss,
                            narrow_entry_lower=float(a['lower']),narrow_entry_upper=float(a['upper']) if a['upper'] else None,
                            broad_entry_lower=float(b['lower']),broad_entry_upper=float(b['upper']) if b['upper'] else None,
                            native_ON_horizon=max(narrow['pulse_durations']),**comparison(a,b)))
    return data


def figure(org,data):
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(11.5,9.7));grid=fig.add_gridspec(2,3,height_ratios=[1,1.65],left=.095,right=.975,bottom=.17,top=.88,wspace=.30,hspace=.62)
    for i,(field,title) in enumerate([('excess_amount_over_reference','A  More localized active signal'),('width_over_reference','B  Wider incumbents'),('peak_over_reference','C  Taller in Goryachev')]):
        ax=fig.add_subplot(grid[0,i])
        for model in LABELS:
            rr=sorted([r for r in org if r['model']==model],key=lambda r:float(r['size_native']))
            ax.plot([float(r['size_over_reference']) for r in rr],[float(r[field]) for r in rr],'-o',ms=2.5,lw=1.6,color=COLORS[model],label=LABELS[model])
        ax.set(title=title,xlabel='Native length / reference length',ylabel='Value / reference value');ax.grid(alpha=.15);ax.set_xlim(.2,2.7)
    ax=fig.add_subplot(grid[1,:]);ax.axhline(1,color='#666',ls='--',lw=1)
    for model in LABELS:
        primary=sorted([r for r in data if r['model']==model and r['roi']=='organization_core' and r['background_rule']=='current_minimum' and r['loss_fraction']==.5 and not r['initial_10pc_core_overlap']],key=lambda r:r['length'])
        finite=[r for r in primary if r['ratio_midpoint'] is not None]
        ax.plot([r['length_over_reference'] for r in finite],[r['ratio_midpoint'] for r in finite],'-o',color=COLORS[model],ms=4,lw=1.9,label=LABELS[model])
        for r in primary:
            same=[s for s in data if s['model']==model and s['length']==r['length'] and s['ratio_midpoint'] is not None]
            if r['ratio_midpoint'] is not None:
                lo=min(s['ratio_lower'] for s in same);hi=max(s['ratio_upper'] for s in same)
                ax.vlines(r['length_over_reference'],lo,hi,color=COLORS[model],alpha=.22,lw=5,zorder=1)
                ax.vlines(r['length_over_reference'],r['ratio_lower'],r['ratio_upper'],color=COLORS[model],lw=1.2,zorder=3)
                if len(same)<18:
                    ax.scatter(r['length_over_reference'],r['ratio_midpoint'],marker='s',s=66,facecolors='none',edgecolors=COLORS[model],lw=.8,zorder=4)
            elif r['ratio_status']=='narrow_no_sampled_entry_through_horizon':
                ax.scatter(r['length_over_reference'],r['ratio_upper'],marker='v',s=40,facecolors='white',edgecolors=COLORS[model],zorder=4)
            elif r['ratio_status']=='broad_no_sampled_entry_through_horizon':
                ax.scatter(r['length_over_reference'],r['ratio_lower'],marker='^',s=40,facecolors='white',edgecolors=COLORS[model],zorder=4)
            elif r['ratio_status']=='neither_sampled_entry_observed':
                ax.scatter(r['length_over_reference'],.02,transform=ax.get_xaxis_transform(),marker='x',s=35,color=COLORS[model],zorder=4)
    ax.set(xlabel='Native length / reference length',ylabel='Broader / narrower cue\nTime to 50% localized excess loss',
        title='D  Cue broadening changes weakening time differently across size',yscale='log',xlim=(.2,2.7),ylim=(.006,3.2))
    ax.grid(alpha=.15);ax.text(.98,.96,'Above 1: broader is slower\nBelow 1: broader is faster',transform=ax.transAxes,ha='right',va='top',fontsize=10)
    handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',ncol=5,bbox_to_anchor=(.5,.944),frameon=False)
    fig.text(.5,.985,'Larger incumbents; length-dependent cue compensation',ha='center',va='top',fontsize=18,weight='bold')
    fig.text(.5,.963,'Provisional main-panel candidate · retained organization plus new full-field experiments',ha='center',va='top',fontsize=10,color='#555')
    fig.text(.075,.110,'A–C: prior 97-basal organization evidence. References: WP/Goryachev/deBelly 19; Holmes 20; Spring resting L₀ = 1.',fontsize=9)
    fig.text(.075,.087,'D: fixed native peak; widths 1.5 → 5 (Spring 0.1 → 0.4). Broader input also increases integrated dose; models are not matched optical actuators.',fontsize=9)
    fig.text(.075,.064,'Pale ranges: observed ratios across 25/50/75% loss, three floors and two cores; sensitivity ranges, not confidence intervals.',fontsize=9)
    fig.text(.075,.041,'Open triangles: finite-window bounds. Hollow squares: incomplete sensitivity. × at axis floor: neither cue has a sampled primary loss entry.',fontsize=9)
    return fig


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--redraw',action='store_true');args=parser.parse_args()
    for name in ['supporting/data','supporting/vectors','supporting/metadata']:(OUT/name).mkdir(parents=True,exist_ok=True)
    orgfile=OUT/'supporting/data/initial_organization_all97.csv';widthfile=OUT/'supporting/data/width_response_all_conventions.csv'
    if args.redraw:
        org=rows(orgfile);data=rows(widthfile)
        numeric=['length','length_over_reference','narrow_width','broad_width','requested_native_peak','narrow_sampled_peak','broad_sampled_peak','loss_fraction','narrow_entry_lower','narrow_entry_upper','broad_entry_lower','broad_entry_upper','native_ON_horizon','ratio_lower','ratio_upper','ratio_midpoint']
        for r in data:
            for key in numeric:r[key]=float(r[key]) if r[key] else None
            r['initial_10pc_core_overlap']=r['initial_10pc_core_overlap']=='True'
    else:
        org=rows(ORG);assert len(org)==97;data=collect();csvwrite(orgfile,org);csvwrite(widthfile,data)
        (OUT/'supporting/metadata/native_evidence_sources.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    fig=figure(org,data)
    fig.savefig(OUT/'figure2_compensation_candidate.png',dpi=200)
    for ext in ['pdf','svg']:fig.savefig(OUT/'supporting/vectors'/('figure2_compensation_candidate.'+ext))
    plt.close(fig)
    caption='''# Provisional candidate: organization and length-dependent compensation

The three organization panels retain the prior 97 exact qualified basals. Localized active amount is excess above the basal minimum within its connected half-height core; width uses its interpolated edge. Concentration and integrated signal remain separate. Spring uses physical length; deBelly retains nominal chemical-solver coordinates. Normalization references are model-specific and unchanged from the retained organization panel.

The new response panel compares broader and narrower fixed-native-peak inputs during continuous stimulation. Its primary quantity is time to first saved 50% loss of active excess in the fixed initial organization core, referenced to the current whole-profile minimum. The ratios compare widths 1.5 and 5 in the three chemical models and corrected linear deBelly, or 0.1 and 0.4 in Spring. These widths change integrated input as well as extent. Equal native peaks are not calibrated common optical doses.

Pale ranges span observed ratios under both fractional organization and whole-cell worker cores, three floor rules (current minimum, current native-sample fifth percentile, frozen basal minimum) and 25/50/75% loss. They are sensitivity ranges, not confidence intervals. A hollow square around a primary point means some sensitivity alternatives have no jointly observed ratio; those alternatives are not replaced by finite values in the range. Crosses placed at the axis floor mean neither width has a sampled primary loss entry, so no ratio is plotted. The fifth percentile is a native-grid concentration quantile. Current-floor Spring excess with evolving physical weights is invariant to pure uniform dilution; frozen concentration-floor subtraction is not. Cases with initial cue overlap of the 10% incumbent core are omitted from the plot and remain explicitly labelled in the full tables.

Open downward/upward triangles report a ratio bound for the first *sampled* loss entry when one width has no entry through its finite observation horizon. This is not an optimized pulse threshold or a claim of infinite response time. Native continuous-ON horizons are WP/Holmes 10000, Goryachev 5000, Spring 50 and deBelly 2000. Saved-entry brackets and every source identity remain in the table.

Regional loss can reflect disappearance of localized organization, redistribution or global concentration change. It is not identified molecular damage or proof of lasting replacement. Original scalar-strength and native field-tail outcomes are retained separately in the comprehensive campaign report. This panel tests compensation during stimulation; native paired fields and withdrawal outcomes are being selected as complementary mechanistic illustrations. The accepted native discretizations are retained; no new spatial-convergence or optical-reporter calibration is claimed.
'''
    (OUT/'CAPTION.md').write_text(caption)
    (OUT/'README.md').write_text('# Provisional Figure 2 compensation candidate\n\n![Retained organization and new response comparison](figure2_compensation_candidate.png)\n\n[PDF](supporting/vectors/figure2_compensation_candidate.pdf) · [Caption and limits](CAPTION.md) · [All width comparisons and sensitivity](supporting/data/width_response_all_conventions.csv).\n\nThe organization panels are prior evidence. The response comparison is new full-field cue experimentation, not a repeat of the optimized reversal-duration plot. Main presentation and explanatory native examples remain under review; this candidate preserves all earlier artwork.\n')
    SOURCES[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    redraw=dict(command='python polarity1d/manuscript/code/mechanisms/revision_compensation_candidate_20261019.py --redraw',
        working_directory='repository root',dependencies=['numpy','matplotlib'],native_archives_required_for_redraw=False,
        source_hashes={str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},
        compact_inputs={str(p.relative_to(ROOT)):sha(p) for p in [orgfile,widthfile]},
        evidence_replay='Native/derived archive identity and full-field numerical replay are documented separately in native_evidence_sources.json and independent audit receipts.')
    (OUT/'supporting/metadata/redraw_closure.json').write_text(json.dumps(redraw,indent=2)+'\n')
    repo=ROOT.parents[1];rel=lambda p:str(p.relative_to(repo))
    content=dict(title='Figure 2. Larger incumbents and length-dependent cue compensation',
        role='Prior incumbent organization motivates a new native test of how broadening compensates for enlargement.',
        status='Provisional candidate: strength operating-point and explanatory-field selection in progress; author acceptance open.',
        assets=[dict(path=rel(OUT/'figure2_compensation_candidate.png'),caption='Retained organization and new broadening-response comparison',page_size=[11.5,9.7],page_title='Figure 2. Incumbent organization and cue compensation',minimum_printed_width_inches=10)],
        legend=[caption],narrative=['Localized active amount grows through different combinations of extent and concentration.','A broader native cue accelerates localized incumbent loss across eligible chemical/Spring size ranges, but its relative advantage changes with length.','Corrected deBelly provides an explicit size-dependent exception; descriptive associations do not isolate causal mediation.'],
        measurements=['Fixed-initial-core positive excess activity above named spatial floors, preserving native amplitude and units.','First sampled loss-entry times and bracket-derived broader/narrower ratios.','Separate scalar-strength and full-field-tail outcomes in the comprehensive campaign supplement.'],
        limitations=['Broader input changes integrated dose; native input channel and operating point differ across models.','Current-profile floors are descriptive concentration references, not molecular pools.','Sensitivity and finite sampled observation bounds are not confidence intervals or optimized-duration thresholds.','Native spatial discretizations retained; no new optical/reporter calibration or isolated causal mediation.'],
        sources=[dict(label='Compact redraw closure',path=rel(OUT/'supporting/metadata/redraw_closure.json')),dict(label='All width/floor/ROI/loss comparisons',path=rel(widthfile)),dict(label='Original unchanged organization artwork',path='current manuscript figs/initial_polarity_length.png')])
    (OUT/'supporting/metadata/content.json').write_text(json.dumps(content,indent=2)+'\n')
    artifacts={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='validation.json'}
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),
        source_hashes=SOURCES,artifact_hashes=artifacts,comparison_rows=len(data),retained_organization_rows=97,
        primary='organization fractional core/current profile minimum/50% excess loss',status='provisional_author_review_independent_display_check_pending'),indent=2)+'\n')
    print(str(OUT/'figure2_compensation_candidate.png'))


if __name__=='__main__':main()
