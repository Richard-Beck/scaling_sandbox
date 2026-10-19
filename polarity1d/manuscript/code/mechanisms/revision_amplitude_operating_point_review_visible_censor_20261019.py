"""Compact actual-result display for the bounded native dose extension."""
from pathlib import Path
import csv
import hashlib
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/amplitude_operating_point_20261019'


def rows(name):
    return list(csv.DictReader((OUT/'supporting/data'/name).open()))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    summary=rows('finite_assay_readouts.csv')
    assert len(summary)==32 and all(r['state']=='complete' for r in summary),'Wait for all 32 scientific protocols.'
    entries=rows('regional_sampled_entries.csv')
    channel=rows('realized_channel_quadrature.csv')
    selected={int(r['slot']):r for r in entries if r['phase']=='cue_on' and r['observable']=='fixed_original_core_active_fraction' and float(r['regional_depletion_fraction'])==.25}
    availability={(int(r['slot']),float(r['window_end'])):r for r in channel}
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    fig.subplots_adjust(left=.09,right=.96,bottom=.17,top=.82,hspace=.55,wspace=.32)
    fig.text(.07,.95,'Separating requested dose from available substrate',fontsize=18,weight='bold')
    fig.text(.07,.90,'Performed native extension: fixed physical width 3; unchanged basals, kinetics and input channel.',fontsize=11)
    fig.text(.07,.865,'Eight peak levels per size; WP pulse 240, Goryachev pulse 1600; then 10000 cue-free native units.',fontsize=10)
    displayed=[]
    for col,(model,title,baseline) in enumerate([('wave_pinning','Wave-pinning',.1),('goryachev','Goryachev',6.4)]):
        ax=axes[0,col];av=axes[1,col]
        subset=[r for r in summary if r['model']==model]
        for length,color in zip(sorted(set(float(r['length']) for r in subset)),['#2166ac','#b2182b']):
            group=sorted([r for r in subset if float(r['length'])==length],key=lambda r:float(r['native_peak_multiplier']))
            for r in group:
                slot=int(r['slot']);event=selected[slot];x=float(r['native_peak_multiplier'])
                observed=event['status'] in ('sampled_transition','already_present')
                lo=float(event['lower']);hi=float(event['upper']) if observed else float(r['pulse_duration'])
                if observed:
                    y=(lo+hi)/2
                    ax.errorbar(x,y,yerr=np.array([[y-lo],[hi-y]]),color=color,marker='o',ms=5,lw=1,capsize=2)
                else:
                    peers=[q for q in subset if float(q['native_peak_multiplier'])==x and selected[int(q['slot'])]['status']=='right_censored_at_observed_horizon']
                    plotted_x=x*(.93 if length==min(float(q['length']) for q in peers) else 1.07) if len(peers)>1 else x
                    ax.scatter(plotted_x,hi,marker='^',s=45,facecolors='none',edgecolors=color)
                record=dict(slot=slot,model=model,length=length,peak_multiplier=x,requested_peak=float(r['requested_peak']),
                    plotted_event_peak_multiplier=x if observed else plotted_x,
                    event_status=event['status'],event_lower=lo,event_upper=event['upper'],finite_observation_horizon=float(r['pulse_duration']))
                for window,style in ((1.,'-'),(10.,'--')):
                    a=availability[slot,window]
                    assert a['cue_roi_inactive_end_over_initial']!=''
                    record[f'inactive_fraction_at_{window:g}']=float(a['cue_roi_inactive_end_over_initial'])
                displayed.append(record)
            for window,style in ((1.,'-'),(10.,'--')):
                av.plot([float(r['native_peak_multiplier']) for r in group],
                    [float(availability[int(r['slot']),window]['cue_roi_inactive_end_over_initial']) for r in group],
                    color=color,ls=style,marker='o',ms=4,lw=1.5)
        for a in (ax,av):
            a.set_xscale('log');a.set_xlim(.007,4.2);a.set_xticks([.01,.03,.1,.3,1,3],['.01','.03','.1','.3','1','3'])
            a.axvspan(.5,2,color='.9',zorder=-2);a.grid(axis='y',alpha=.15)
            a.set_xlabel(f'Requested peak / native baseline ({baseline:g})',fontsize=9)
            a.tick_params(labelsize=9)
        ax.set_yscale('log');ax.set_ylim(.5 if col==0 else 1.,350 if col==0 else 2400)
        ax.set_ylabel('Time to 25% original-region activity loss\n(native cue-on units)',fontsize=9)
        ax.set_title(f'{"A" if col==0 else "B"}  {title}: finite weakening time',loc='left',fontsize=11,weight='bold')
        av.set_ylim(0,1.05);av.set_ylabel('Cue-ROI inactive concentration / initial',fontsize=9)
        av.set_title(f'{"C" if col==0 else "D"}  Available inactive at exact saved times',loc='left',fontsize=11,weight='bold')
        handles=[Line2D([],[],color=color,marker='o',label=f'L{length:g}') for length,color in zip(sorted(set(float(r['length']) for r in subset)),['#2166ac','#b2182b'])]
        ax.legend(handles=handles,frameon=False,fontsize=8)
        av.legend(handles=[Line2D([],[],color='black',ls='-',label='t = 1'),Line2D([],[],color='black',ls='--',label='t = 10')],frameon=False,fontsize=8)
    fig.text(.075,.105,'Bars: saved-time brackets, not confidence intervals. Open triangles: no event by pulse end; coincident size markers offset in x.',fontsize=9)
    fig.text(.075,.072,'Grey band: prior 0.5–2 peak range. Exact inactive endpoint samples are observations; legacy channel integrals remain cadence-unqualified.',fontsize=9)
    fig.text(.075,.039,'Native stimulus × inactive; fractional peak matching does not match effective actuation. N = 201; no optical-clock, actuator or spatial-convergence claim.',fontsize=8.5)
    png=OUT/'figure2_amplitude_operating_point.png';pdf=OUT/'supporting/vectors/figure2_amplitude_operating_point.pdf';svg=OUT/'supporting/vectors/figure2_amplitude_operating_point.svg'
    for p in (png,pdf,svg):fig.savefig(p,dpi=170)
    plt.close(fig)
    table=OUT/'supporting/data/operating_point_display.csv'
    with table.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(displayed[0]));w.writeheader();w.writerows(displayed)
    sources={str((OUT/'supporting/data'/name).relative_to(ROOT)):sha(OUT/'supporting/data'/name) for name in ['finite_assay_readouts.csv','regional_sampled_entries.csv','realized_channel_quadrature.csv']}
    sources[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    artifacts={str(p.relative_to(OUT)):sha(p) for p in [png,pdf,svg,table]}
    (OUT/'supporting/metadata/review_validation.json').write_text(json.dumps(dict(source_hashes=sources,artifact_hashes=artifacts,displayed_protocols=32,loss_threshold=.25),indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),displayed_protocols=len(displayed)),indent=2))


if __name__=='__main__':main()
