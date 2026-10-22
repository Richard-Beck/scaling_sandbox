"""Old-pole Rac concentration in canonical native reversals; no integration."""
from pathlib import Path
from datetime import datetime, timezone
import csv, gzip, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'figures/length_mechanisms_20261008/dev/native_reversal'
RUN=ROOT/'data/runs/length_narrative_memory_20261009'
BASE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
STEM='old_pole_Rac_trajectories'
LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
GEOS=['fixed_physical','fixed_relative','fixed_slope']
TITLES=['Fixed patch: 3 µm','Relative patch: 20%','Constant slope: 0.24']
LIMIT=1200.
SOURCES={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def register(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return p
def read(p):return json.loads(register(p).read_text())
def output(name):return layout.dev_output_path(OUT,name)
def write(name,v):output(name).write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')
def csvwrite(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    register(Path(__file__));register(Path(layout.__file__))
    manifest=read(RUN/'manifest.json')
    assert manifest['center_boundary_motion'] is False
    assert read(RUN/'status.json')['audit_passed']
    assert read(BASE/'results/independent_audit.json')['passed']
    prior_audit=read(ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/memory/independent_audit.json')
    assert prior_audit['passed']
    records_path=register(BASE/'results/records.json.gz')
    with gzip.open(records_path,'rt') as f:
        records={r['case_index']:r for r in json.load(f)['cases']}
    traces={};measurements=[];audits=[];arrays={};csvrows=[]
    for job in manifest['jobs']:
        c=job['case'];L=int(c['length']);geo=c['condition'];D=float(job['duration'])
        assert L in LENGTHS and geo in GEOS and D<LIMIT
        assert (geo,L) not in traces
        record=records[c['case_index']]
        assert D==record['result']['full_polarity']['endpoint']
        assert D==min(r['duration'] for r in record['tests'] if r['full_polarity_success'])
        path=register(ROOT/job['source'] if job['source'] else RUN/f"workers/{manifest['missing_jobs'].index(job):03d}.npz")
        meta=read(path.with_suffix('.json'));assert sha(path)==meta['archive_sha256']
        cp=register(ROOT/c['checkpoint']);assert sha(cp)==c['checkpoint_sha256']
        with np.load(cp) as z:initial=z['state_rac'].copy()
        with np.load(path) as z:
            x=z['x'].copy();tp=z['pulse_time'].copy();tr=z['recovery_time'].copy()
            pulse=z['pulse_field_rac'];recovery=z['recovery_field_rac']
            assert np.isfinite(pulse).all() and np.isfinite(recovery).all()
            np.testing.assert_allclose(x,(np.arange(c['cells'])+.5)*L/c['cells'],rtol=0,atol=1e-12)
            np.testing.assert_array_equal(pulse[0],initial)
            np.testing.assert_array_equal(recovery[0],pulse[-1])
            assert tp[0]==0 and tr[0]==0 and abs(tp[-1]-D)<1e-9 and tr[-1]==10000
            assert np.all(np.diff(tp)>0) and np.all(np.diff(tr)>0)
            tail=tr>=9500;rac=recovery[tail]
            polarity=rac@(2*x/L-1)/rac.sum(axis=1);contrast=np.ptp(rac,axis=1)
            strength_error=float(np.max(abs(polarity/abs(c['initial_orientation'])-1)))
            assert polarity.min()>.05 and contrast.min()>=.05 and np.ptp(polarity)<.005 and strength_error<=.01
            assert meta['metrics']['full_polarity_success'] and meta['metrics']['field_tail_pass']==job['expected']['field_tail_pass']
            yp=pulse[:,0].copy();yr=recovery[:,0].copy()
        t=np.r_[tp,D+tr[1:]];y=np.r_[yp,yr[1:]]
        assert len(t)==len(y) and np.all(np.diff(t)>0) and y[0]==initial[0]
        flag=meta['metrics']['field_tail_pass']
        traces[geo,L]=(tp,yp,D+tr,yr,D,flag)
        prefix=f'{geo}_L{L}'
        arrays[prefix+'_time_s']=t;arrays[prefix+'_old_pole_rac']=y
        arrays[prefix+'_cue_duration_s']=np.array(D)
        arrays[prefix+'_old_pole_nominal_x_um']=np.array(x[0])
        arrays[prefix+'_field_tail_pass']=np.array(flag)
        for ti,yi in zip(t,y):
            csvrows.append(dict(condition=geo,length_um=L,time_from_onset_s=float(ti),old_pole_rac=float(yi),phase='cue' if ti<=D else 'basal_recovery'))
        visible=t<=LIMIT;index=np.flatnonzero(visible)[np.argmin(y[visible])]
        measurements.append(dict(condition=geo,length_um=L,old_pole_nominal_x_um=float(x[0]),cue_duration_s=D,
            initial_old_pole_rac=float(yp[0]),old_pole_rac_at_withdrawal=float(yp[-1]),
            old_pole_rac_at_1200_s=float(np.interp(LIMIT,t,y)),final_old_pole_rac=float(yr[-1]),
            minimum_old_pole_rac_in_display=float(y[index]),time_of_display_minimum_s=float(t[index]),
            native_field_tail_pass=flag,maximum_scalar_strength_relative_error=strength_error,
            archive=str(path.relative_to(ROOT)),checkpoint_sha256=c['checkpoint_sha256']))
        audits.append(dict(length_um=L,condition=geo,passed=True,archive_sha256=sha(path),checkpoint_sha256=c['checkpoint_sha256'],
            exact_initial_Rac=True,exact_phase_continuity=True,original_leftmost_cell_used=True,
            shortest_tested_strength_pulse_verified=True,scalar_gate_independently_recomputed=True,
            field_tail_pass=flag,maximum_scalar_strength_relative_error=strength_error))
    assert set(traces)=={(geo,L) for geo in GEOS for L in LENGTHS}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    colors={L:plt.get_cmap('turbo')(.06+.88*i/11) for i,L in enumerate(LENGTHS)}
    fig,axs=plt.subplots(1,3,figsize=(14.2,5.0),layout='constrained',sharex=True,sharey=True)
    for ax,geo,title,letter in zip(axs,GEOS,TITLES,'ABC'):
        ax.set_title(f'{letter}   {title}',loc='left',fontweight='bold',pad=10)
        for L in LENGTHS:
            tp,yp,trec,yr,D,flag=traces[geo,L]
            ax.plot(tp,yp,color=colors[L],lw=1.5)
            ax.plot(trec,yr,color=colors[L],lw=1.5,ls='--')
            ax.plot(D,yp[-1],'o',ms=4.7,mec=colors[L],mfc=colors[L] if flag else 'white',mew=1.1,zorder=4)
        ax.set(xlim=(0,LIMIT),xlabel='Time from reversal-cue onset (s)')
        ax.set_xticks([0,300,600,900,1200]);ax.grid(alpha=.13)
    visible=[np.r_[yp,yr[trec<=LIMIT]] for tp,yp,trec,yr,D,flag in traces.values()]
    lo=min(v.min() for v in visible);hi=max(v.max() for v in visible);pad=.06*(hi-lo)
    axs[0].set_ylim(lo-pad,hi+pad)
    axs[0].set_ylabel('Old-pole Rac (native concentration)')
    fig.suptitle('Rac at the old pole during native deBelly reversal',fontsize=15,fontweight='bold')
    handles=[Line2D([],[],color=colors[L],lw=2,label=str(L)) for L in LENGTHS]
    handles += [Line2D([],[],color='#444444',lw=1.5,label='Cue on'),Line2D([],[],color='#444444',lw=1.5,ls='--',label='Basal recovery'),
        Line2D([],[],color='#444444',marker='o',ms=5,ls='',mfc='#444444',label='Final fields stable'),
        Line2D([],[],color='#444444',marker='o',ms=5,ls='',mfc='white',label='Final fields changing')]
    fig.legend(handles=handles,loc='outside lower center',ncol=8,frameon=False,fontsize=8.5,
               title='Cell length (µm) · circles mark withdrawal of each shortest tested strength-passing pulse')
    hashes={}
    for ext in ['png','pdf','svg']:
        path=output(STEM+'.'+ext);fig.savefig(path,dpi=220,facecolor='white');hashes[str(path.relative_to(OUT))]=sha(path)
    plt.close(fig)
    path=output(STEM+'_summary.csv');csvwrite(path,measurements);hashes[str(path.relative_to(OUT))]=sha(path)
    path=output(STEM+'_samples.csv');csvwrite(path,csvrows);hashes[str(path.relative_to(OUT))]=sha(path)
    path=output(STEM+'_traces.npz');np.savez_compressed(path,**arrays);hashes[str(path.relative_to(OUT))]=sha(path)
    caption=('Absolute Rac concentration at the old left pole during the original native-model deBelly reversals, with corrected independent endpoint commands and native protrusion retained. Three panels compare fixed support 3 µm, relative support 20%, and constant physical slope 0.24. Colors identify all twelve eligible lengths: 18,19,21,23,25,27,29,31,33,35,37,39 µm. The old pole is the first native finite-volume cell, at nominal x=L/(2×201), followed throughout the trial; it becomes the rear after reversal. This is a local concentration, not a whole-cell polarity moment, a moving peak, or a normalized signal. Every trajectory begins at the exact shared corrected basal checkpoint and uses that length and cue geometry’s shortest tested canonical 1%-strength-passing finite pulse. Solid lines indicate cue-on; dashed lines show cue-free basal recovery; circles mark each exact withdrawal. Filled circles denote a passing final native-field tail flag and open circles a changing final field tail; marker fill does not report stability at the moment of withdrawal. Time runs from cue onset on a shared linear 0–1200-s axis, retaining all 36 withdrawal times. Full recorded recovery through 10000 s after withdrawal is exported in NPZ and CSV. All 36 pulse selections, source hashes, exact initial Rac profiles, exact pulse/recovery continuity, native coordinate indexing and final scalar reversal/contrast/strength gates were checked independently. Previously audited native-field tail flags remain separate from scalar success. No new simulations or changes to the model, pulse durations, or main report.')
    path=output(STEM+'_caption.txt');path.write_text(caption+'\n');hashes[str(path.relative_to(OUT))]=sha(path)
    write(STEM+'_audit.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,checks=audits,trajectories=36,new_integrations=0))
    path=output(STEM+'_audit.json');hashes[str(path.relative_to(OUT))]=sha(path)
    write(STEM+'.manifest.json',dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',legend=caption,
        sha256=hashes,source_sha256=SOURCES,independent_audit_passed=True,trajectories=36,new_integrations=0,model='native deBelly with independent endpoints',
        displayed_time_s=[0,LIMIT],selection='shortest tested canonical 1%-strength-passing finite cue',
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))]))
    print(output(STEM+'.png'))

if __name__=='__main__':main()
