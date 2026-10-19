"""Paired native Goryachev width histories; no integrations or changed inputs."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import revision_width_fullfield_review_20261019 as native

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure2/paired_width_timecourses_20261019'
COLORS={'narrow':'#2171b5','broad':'#983ca0'}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def selected(time,phase):
    if phase=='cue_on':step=np.full(len(time),.5)
    else:step=np.where(time<=100,.25,np.where(time<=1000,2.5,10.))
    keep=np.isclose(time/step,np.rint(time/step),atol=1e-9,rtol=0)
    keep[[0,-1]]=True
    return np.flatnonzero(keep)


def main():
    for name in ('supporting/data','supporting/vectors','supporting/metadata','supporting/viewers'):(OUT/name).mkdir(parents=True,exist_ok=True)
    cases,sources=native.read_cases()
    cases=sorted([d for d in cases if d['case']['model']=='goryachev'],key=lambda d:(d['case']['length'],d['case']['width']))
    assert len(cases)==4
    fig,axes=plt.subplots(2,3,figsize=(14,9.4))
    fig.subplots_adjust(left=.085,right=.97,bottom=.25,top=.80,wspace=.32,hspace=.57)
    fig.text(.065,.947,'Broader challengers can lose or win after cue withdrawal',fontsize=19,weight='bold')
    fig.text(.065,.90,'Goryachev: exact shared basal and requested peak 6.4 within each pair; the larger cohort receives a longer pulse.',fontsize=11)
    fig.text(.065,.862,'Read the amount, original-region weakening and local concentration together; no single peak or depletion threshold predicts both pairs.',fontsize=10)
    curves=[];endpoints=[];artists=[]
    for row,length in enumerate([9.,25.]):
        pair=[d for d in cases if d['case']['length']==length]
        initial_core=None
        for d in pair:
            c,z=d['case'],d['masks'];dx=c['length']/c['cells'];core=z['incumbent_core_mask'];roi=z['cue_footprint_mask'];distal=z['x']>=c['length']/2
            core0=float(z['initial_state'][0,core].sum()*dx)
            if initial_core is None:initial_core=core0
            assert abs(initial_core-core0)<1e-13
            time_all=[];values=[[],[],[]]
            for phase,p in d['phases'].items():
                ids=selected(p['time'],phase);t=p['time'][ids];u=p['field_active'][ids]
                relative=t-c['cue_duration'] if phase=='cue_on' else t
                yy=[u[:,distal].sum(axis=1)*dx,u[:,core].sum(axis=1)*dx/core0,u[:,roi].max(axis=1)]
                time_all.extend(relative.tolist())
                for col,y in enumerate(yy):values[col].extend(y.tolist())
                for i in range(len(t)):
                    curves.append(dict(length=c['length'],physical_width=c['width'],pulse_duration=c['cue_duration'],phase=phase,
                        phase_local_time=float(t[i]),time_relative_to_withdrawal=float(relative[i]),
                        fixed_distal_active_amount=float(yy[0][i]),original_core_active_fraction=float(yy[1][i]),
                        actual_cue_ROI_active_peak=float(yy[2][i])))
            for col,y in enumerate(values):
                axes[row,col].plot(time_all,y,color=COLORS[c['width_label']],lw=1.9,label=f'w = {c["width"]:g}')
                artists.append(dict(length=length,width=c['width'],observable=col,
                    time_sha256=hashlib.sha256(np.asarray(time_all).tobytes()).hexdigest(),
                    values_sha256=hashlib.sha256(np.asarray(y).tobytes()).hexdigest(),points=len(y)))
            on,off=d['phases']['cue_on'],d['phases']['off']
            endpoints.append(dict(length=c['length'],physical_width=c['width'],pulse_duration=c['cue_duration'],
                requested_peak=c['response_amplitude'],sampled_peak=c['sampled_peak'],
                sampled_spatial_stimulus_integral=c['sampled_spatial_integral'],
                initial_core_active_amount=core0,initial_core_physical_width=float(core.sum()*dx),initial_native_peak=float(z['initial_state'][0].max()),
                distal_region_first_cell_index=int(np.flatnonzero(distal)[0]),distal_region_cells=int(distal.sum()),dx=dx,
                distal_region_left_face=float(np.flatnonzero(distal)[0]*dx),distal_region_right_face=c['length'],
                withdrawal_distal_active_amount=float(on['field_active'][-1,distal].sum()*dx),
                withdrawal_core_active_fraction=float(on['field_active'][-1,core].sum()*dx/core0),
                withdrawal_challenger_peak=float(on['field_active'][-1,roi].max()),
                final_distal_active_amount=float(off['field_active'][-1,distal].sum()*dx),
                final_core_active_fraction=float(off['field_active'][-1,core].sum()*dx/core0),
                final_challenger_peak=float(off['field_active'][-1,roi].max()),
                **{k:d['result']['legacy_diagnostics'][k] for k in ['orientation_success','full_polarity_success','field_tail_pass']}))
        for col,title in enumerate(['Fixed distal active amount','Original-core activity','Challenger peak in actual cue ROI']):
            ax=axes[row,col]
            ax.set_xscale('symlog',linthresh=100);ax.set_xlim(-pair[0]['case']['cue_duration'],10000)
            ax.axvspan(-pair[0]['case']['cue_duration'],0,color='.92',zorder=-2)
            ax.axvline(0,color='.35',ls='--',lw=1)
            ticks=[-pair[0]['case']['cue_duration'],0,100,1000,10000]
            ax.set_xticks(ticks,[f'{v:g}' for v in ticks]);ax.tick_params(labelsize=8)
            ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.15)
            ax.set_title(f'{"ABCDEF"[row*3+col]}  {title}',loc='left',fontsize=10.5,weight='bold')
            ax.set_xlabel('Time relative to withdrawal (native)',fontsize=9)
            ax.set_ylabel(['Native activity amount','Active amount / initial core amount','Native active concentration'][col],fontsize=9)
        axes[row,0].legend(title=f'L{length:g}; pulse {pair[0]["case"]["cue_duration"]:g}',frameon=False,fontsize=8,title_fontsize=9,loc='upper left')
    fig.text(.075,.18,'L9: the broader cue leaves more distal material but a lower peak; that response disappears while the original region recovers.',fontsize=10)
    fig.text(.075,.144,'L25: the broader cue also lowers the challenger peak, but now the original region is depleted more and the challenger is retained.',fontsize=10)
    fig.text(.075,.105,'Within-pair basals and peak are identical. The two rows change length, pulse duration and the width interval; they do not isolate a causal length effect.',fontsize=9)
    fig.text(.075,.074,'Grey = cue on. Time axis is linear for |t| ≤ 100 and logarithmic outside. Curves use selected actual saved frames; all native field tails pass.',fontsize=9)
    fig.text(.075,.043,'Distal ROI = whole cells with centres ≥ L/2; each pair uses identical weights. Native N = 201; input × inactive, dose and clock are not optically calibrated.',fontsize=8.5)
    outputs=[]
    for path in [OUT/'figure2_paired_width_timecourses.png',OUT/'supporting/vectors/figure2_paired_width_timecourses.pdf',OUT/'supporting/vectors/figure2_paired_width_timecourses.svg']:
        fig.savefig(path,dpi=170);outputs.append(path)
    plt.close(fig)
    for name,table in [('paired_native_curves.csv',curves),('paired_native_endpoints.csv',endpoints)]:
        path=OUT/'supporting/data'/name
        with path.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
        outputs.append(path)
    (OUT/'README.md').write_text('''# Paired native width histories: Figure 2 option

[One-sheet PDF](supporting/vectors/figure2_paired_width_timecourses.pdf) · [Native plotted samples](supporting/data/paired_native_curves.csv) · [Withdrawal/final values and source flags](supporting/data/paired_native_endpoints.csv)

![Paired native time courses](figure2_paired_width_timecourses.png)

The same increased width can accompany more distal active material and a lower local peak in both cohorts, but the retention direction differs. At L9 the broader response disappears; at L25 it persists. This illustrates why incumbent weakening and challenger maintenance need to be examined together. The larger native incumbent is part of the context, not an experimentally isolated explanation of these two width contrasts: the pulse and width interval also differ.

The exact basal and requested peak 6.4 are shared within each pair. The distal region uses whole native cells whose centres satisfy $x_i \\ge L/2$; the same cell weights are applied to both widths at that length. Original-core activity is normalized to its exact common initial amount, without independent curve scaling. Peaks are native active concentration in the actual physical cue ROI. Native stimulus times inactive recruits activity locally while nonlinear maintenance and transport continue. Broadening also changes integrated input.

Source fields were saved every 0.1 cue-on unit and every 0.25 cue-free unit through 10000 after withdrawal. Plotted samples use 0.5 cue-on units, 0.25 off units through 100, 2.5 through 1000 and 10 thereafter, including exact endpoints. These are finite histories, not minimum-pulse measurements. All four original scalar-strength and native field-tail flags are retained separately in the endpoint table. The original all-five amount/width/peak organization panel and validated coupled vignette remain unchanged. The [native state-law note](../coupled_state_law_20261019/README.md) provides bounded mechanistic context without assigning the final outcome a unique cause.
''')
    outputs.append(OUT/'README.md');sources[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    (OUT/'supporting/metadata/validation.json').write_text(json.dumps(dict(source_hashes=sources,
        artifact_hashes={str(p.relative_to(OUT)):sha(p) for p in outputs},line_artists=artists,
        plotted_rows=len(curves),native_cases=4,no_integration=True),indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),plotted_rows=len(curves)),indent=2))


if __name__=='__main__':main()
