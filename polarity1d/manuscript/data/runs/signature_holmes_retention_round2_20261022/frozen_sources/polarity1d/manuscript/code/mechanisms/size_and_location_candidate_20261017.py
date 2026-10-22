"""One successor candidate: conditional responses, actual size fate and location.

Only render retained native data. The preceding compact candidate and all
scientific campaigns are immutable inputs, not write targets.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import compact_integrated_candidate_20261017 as base

OUT=base.ROOT/'figures/length_mechanisms_20261008/dev/size_and_location_candidate_20261017'
PREVIOUS=base.ROOT/'figures/length_mechanisms_20261008/dev/compact_integrated_candidate_20261017'
OLD_COLOR='#3976a8'
CUE_COLOR='#c96926'


def goryachev_histories(matrix,ring_manifest):
    cases=[]
    for C in (40,80):
        row=base.pick(matrix,'goryachev',C,.5)[0]
        phases=[]
        for phase,offset in (('cue',0),('recovery',1000)):
            path=Path(row[phase+'_path'])
            assert base.sha(path)==ring_manifest['source_sha256'][str(path)]
            with np.load(base.source(path)) as z:
                time=z['time'].copy()+offset
                primary=(z['states'][:,0] if 'states' in z else z['primary']).copy()
                x=z['x'].copy()
                stimulus=z['stimulus'].copy() if phase=='cue' else phases[0]['stimulus']
            phases.append(dict(time=time,primary=primary,x=x,stimulus=stimulus))
        cue,recovery=phases
        np.testing.assert_array_equal(cue['x'],recovery['x'])
        np.testing.assert_array_equal(cue['primary'][-1],recovery['primary'][0])
        assert np.argmax(cue['primary'][0])==0
        cell=int(float(row['cue_center_cell']))
        time=np.r_[cue['time'],recovery['time'][1:]]
        primary=np.r_[cue['primary'],recovery['primary'][1:]]
        assert time[0]==0 and time[-1]==14000
        assert np.max(stimulus)==6.4
        np.testing.assert_allclose(np.sum(stimulus)*.05,9.6,rtol=0,atol=2e-12)
        trace=dict(time=time,old_coordinate_concentration=primary[:,0],cue_center_concentration=primary[:,cell],
            old_native_coordinate=cue['x'][0],cue_native_coordinate=cue['x'][cell],
            circumference=C,basal_primary=primary[0],final_primary=primary[-1],stimulus=stimulus)
        path=base.layout.dev_output_path(OUT,f'Goryachev_C{C}_exact_fixed_coordinate_histories.npz')
        np.savez_compressed(path,**trace);base.ARTIFACTS.append(path)
        cases.append((row,trace))
    a,b=(r[1]['stimulus'] for r in cases)
    np.testing.assert_allclose(a[a>0],b[b>0],rtol=0,atol=2e-13)
    return cases


def figure1(folds,curves,cases):
    fig=plt.figure(figsize=(8.5,6.5))
    axes=[fig.add_axes(bounds) for bounds in ((.135,.56,.225,.31),(.475,.56,.195,.31),(.785,.56,.195,.31))]
    ax=axes[0]
    for key,marker,label,offset in [('excess_active_fold','o','Excess active amount',-.19),
            ('width_fold','s','Core width',0),('peak_concentration_fold','^','Peak concentration',.19)]:
        ax.scatter([float(r[key]) for r in folds],np.arange(5)+offset,marker=marker,s=25,
            facecolor='#333333' if marker=='o' else 'white',edgecolor='#333333',lw=.8,label=label,zorder=4)
    names=[r['model_label']+'\n'+('L₀ 1→2' if r['model']=='spring' else f"{r['nominal_size_small']}→{r['nominal_size_large']}") for r in folds]
    ax.set(yticks=np.arange(5),yticklabels=names,xlim=(.8,3.2),xticks=[1,2,3],ylim=(4.6,-.6),
        xlabel='Larger / smaller value\n(original interval preparations)')
    ax.tick_params(axis='y',length=0,labelsize=7.2)
    ax.axvline(1,color='#aaaaaa',ls='--',lw=.7)
    ax.set_title('A  Incumbent organization',loc='left',fontsize=9,pad=10)
    ax.legend(loc='upper left',bbox_to_anchor=(-.14,1.47),frameon=False,fontsize=7,handlelength=1.4,labelspacing=.25)
    for ax,key,title in zip(axes[1:],('initial_gain_relative_to_C19','basal_conditional_repair_relative_to_C19'),
            ('B  Initial input gain','C  Basal conditional repair')):
        for model in base.MODELS:
            rr=[r for r in curves if r['model']==model]
            ax.plot([r['circumference'] for r in rr],[r[key] for r in rr],color=base.COLORS[model],
                marker=base.MARKERS[model],ms=3.1,lw=1.2,mec='white',mew=.25)
            for C in (19,39):
                r=next(r for r in rr if r['circumference']==C)
                ax.scatter([C],[r[key]],s=33,facecolors='none',edgecolors=base.COLORS[model],lw=.8,zorder=5)
        ax.axhline(1,color='#aaaaaa',lw=.75,ls='--')
        ax.set(xlim=(7,83),xticks=[10,40,80],xlabel='Ring circumference C',ylabel='Relative to own C = 19 value')
        ax.set_title(title,loc='left',fontsize=9,pad=10)
        ax.grid(axis='y',color='#eeeeee',lw=.6,zorder=0)
    axes[1].set(ylim=(0,2.08),yticks=[0,.5,1,1.5,2])
    axes[2].set_yscale('log',base=2)
    axes[2].set(ylim=(.2,17),yticks=[.25,.5,1,2,4,8,16],yticklabels=['¼','½','1','2','4','8','16'])
    axes[2].set_ylabel('Relative to own C = 19 value\n(log₂ scale)')
    handles=[Line2D([],[],color=base.COLORS[m],marker=base.MARKERS[m],lw=1.2,ms=4,label=base.NAMES[m]) for m in base.MODELS]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.70,.995),ncol=3,frameon=False,fontsize=7.4,
        handlelength=1.4,columnspacing=1.2)
    fig.text(.70,.945,'B–C: native rings; own C19 references and native units',ha='center',fontsize=7.2,color='#555555')
    fig.text(.55,.465,'B: antipodal width-3 input. C: 0.5% inactive-pool hold/return, hold 0.1.',ha='center',fontsize=7.1)
    fig.text(.55,.443,'Positive conditional repair does not establish net core gain or actual retention.',ha='center',fontsize=7.1)
    lower=[fig.add_axes(bounds) for bounds in ((.135,.115,.33,.235),(.59,.115,.33,.235))]
    for ax,(row,z),letter in zip(lower,cases,('D','E')):
        keep=z['time']<=4000
        ax.axvspan(0,1000,color='#eeeeee',zorder=0)
        ax.plot(z['time'][keep],z['old_coordinate_concentration'][keep],color=OLD_COLOR,lw=1.5)
        ax.plot(z['time'][keep],z['cue_center_concentration'][keep],color=CUE_COLOR,lw=1.5)
        ax.axvline(1000,color='#888888',ls='--',lw=.8)
        ax.text(1060,96,'Cue off',color='#555555',fontsize=7,ha='left',va='center')
        title='Replacement' if letter=='D' else 'Incumbent retained'
        ax.set_title(f'{letter}  Goryachev C{int(z["circumference"])}: {title}',loc='left',fontsize=9,pad=8)
        ax.set(xlim=(0,4000),ylim=(-3,125),xticks=[0,1000,2000,4000],yticks=[0,40,80,120],
            xlabel='Time from cue onset (native units)',ylabel='Active concentration\nat fixed coordinate')
        ax.grid(axis='y',color='#eeeeee',lw=.6)
    fig.legend(handles=[Line2D([],[],color=OLD_COLOR,lw=1.5,label='Original peak coordinate'),
            Line2D([],[],color=CUE_COLOR,lw=1.5,label='Cue center')],loc='upper center',bbox_to_anchor=(.55,.409),
        ncol=2,frameon=False,fontsize=7.5,handlelength=2,columnspacing=2.3)
    fig.text(.55,.035,'D–E: same peak-6.4, width-3 cue for 1000 units; size, material and both routes change.\nShared raw amplitude and linear 0–4000 time axes; endpoint outcomes verified through 14,000.',
        ha='center',va='center',fontsize=7,linespacing=1.4)
    base.save(fig,'size_effects_and_actual_competition')


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'xtick.labelsize':7,
        'ytick.labelsize':7,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    OUT.mkdir(parents=True,exist_ok=True)
    base.OUT=OUT
    base.source(__file__);base.source(base.__file__)
    previous=json.loads(base.source(PREVIOUS/'supporting/metadata/manifest.json').read_text())
    generator_rel=str(Path(base.__file__).resolve().relative_to(base.REPO))
    assert base.sha(base.__file__)==previous['source_sha256'][generator_rel]
    ring_manifest=json.loads(base.source(base.RING/'supporting/metadata/ring_functional_report.manifest.json').read_text())
    matrix=base.read_csv(base.RUN/'results/matrix.csv');basal=base.read_csv(base.RUN/'results/basal.csv')
    folds=base.read_csv(PREVIOUS/'supporting/data/all_five_original_interval_folds.csv')
    for name in ('matrix','basal'):
        path=base.RUN/f'results/{name}.csv'
        assert base.sha(path)==ring_manifest['source_sha256'][str(path)]
    assert len(matrix)==867 and sum(r['status']=='complete' for r in matrix)==652
    curves=base.size_curves(matrix,basal)
    base.export_csv('ring_native_size_curves.csv',curves)
    base.export_csv('all_five_original_interval_folds.csv',folds)
    base.export_csv('complete_endpoint_map.csv',[{k:r[k] for k in ('model','circumference','fraction','realized_fraction','short_route','long_route','status','final_location_association','checkpoint_sha256')} for r in matrix])
    cases=goryachev_histories(matrix,ring_manifest)
    figure1(folds,curves,cases)
    selected=base.figure2(matrix)
    for row,z in selected:
        for key in ('cue_path','recovery_path'):
            path=Path(row[key]);assert base.sha(path)==ring_manifest['source_sha256'][str(path)]
    for rel,digest in previous['artifacts'].items():
        assert base.sha(PREVIOUS/rel)==digest,rel
    checks=dict(planned_map_rows=len(matrix),qualified_size_curves=len(curves),native_integrations_performed=0,
        previous_compact_artifacts_unchanged=True,previous_generator_unchanged=True,
        Goryachev_shared_native_cue=dict(peak=6.4,width=3,integral=9.6,duration=1000),
        history_display=dict(readout='native active concentration at original basal-peak coordinate and fixed cue center',
            time_window=[0,4000],complete_history_endpoint=14000,shared_linear_time=True,shared_raw_amplitude_limits=[-3,125]),
        Goryachev_C80_at_withdrawal=dict(old=float(cases[1][1]['old_coordinate_concentration'][np.flatnonzero(cases[1][1]['time']==1000)[0]]),
            challenger=float(cases[1][1]['cue_center_concentration'][np.flatnonzero(cases[1][1]['time']==1000)[0]])),
        Goryachev_endpoint_values={str(int(z['circumference'])):dict(old=float(z['old_coordinate_concentration'][-1]),cue=float(z['cue_center_concentration'][-1])) for row,z in cases},
        Holmes_exact_complete_basal_equality=True,debelly_in_main_ring_panels=False)
    path=base.layout.dev_output_path(OUT,'verification.json');path.write_text(json.dumps(checks,indent=2)+'\n');base.ARTIFACTS.append(path)
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root='../..',source_sha256=base.SOURCES,
        artifacts={str(p.relative_to(OUT)):base.sha(p) for p in base.ARTIFACTS},checks=checks)
    path=base.layout.dev_output_path(OUT,'manifest.json');path.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),artifacts=len(base.ARTIFACTS),checks=checks),indent=2))


if __name__=='__main__': main()
