"""Steady MCA profiles with a fine native continuation sweep below 18 microns."""
from pathlib import Path
from dataclasses import replace
import argparse, copy, json, shutil, subprocess, sys, time
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if p.name=='manuscript')
OUT=ROOT/'data/runs/debelly_mca_threshold_profiles_20261010'
OLD=ROOT/'data/runs/debelly_short_length_diagnostic_20261010_v2'
BASE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
FIG=ROOT/'figures/length_mechanisms_20261008/dev/mca_threshold_profiles'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
NEW_LENGTHS=[17.8,17.825,17.85,17.875,17.9,17.925,17.95,17.975,17.99,18.05,18.1,18.25,18.5]
def read(p):return json.loads(Path(p).read_text())
def modules():
    sys.path.insert(0,str(OUT/'frozen_sources/code/core'))
    import native as n
    import polarity_establishment as e
    return n,e
def init():
    assert not OUT.exists(),'Fresh destination required'
    for d in ['frozen_sources','checkpoints','workers','logs','results']:(OUT/d).mkdir(parents=True)
    prior=read(OLD/'manifest.json');sources={}
    for name,digest in prior['sources'].items():
        if name.startswith('code/mechanisms/'):continue
        dest=OUT/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(OLD/'frozen_sources'/name,dest);sources[name]=digest
    script=OUT/'frozen_sources'/Path(__file__).relative_to(ROOT)
    script.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(Path(__file__),script)
    n,e=modules();sources[str(Path(__file__).relative_to(ROOT))]=n.sha(script)
    for name,digest in sources.items():assert n.sha(OUT/'frozen_sources'/name)==digest
    shutil.copy2(OLD/'start_L18.npz',OUT/'start_L18.npz')
    inputs=[]
    manifest=read(BASE/'manifest.json')
    for L in [18.,19.,21.,23.,25.,31.,39.]:
        c=next(c for c in manifest['cases'] if c['model']=='debelly' and c['length']==L)
        cp=ROOT/c['checkpoint'];assert n.sha(cp)==c['checkpoint_sha256']
        dest=OUT/'checkpoints'/f'L{L:g}.npz';shutil.copy2(cp,dest)
        prep=ROOT/c['corrected_preparation'];summary=read(prep)['summary']
        trial=next(t for t in read(prep)['trials'] if abs(t['duration']-summary['cue_duration'])<1e-9)
        inputs.append(dict(length=L,file=str(dest.relative_to(OUT)),source=str(cp.relative_to(ROOT)),
            sha256=n.sha(dest),kind='original_preparation',convergence=trial['convergence']))
    for slot in [1,2,3]:
        p=OLD/'workers'/f'{slot}.npz';r=read(p.with_suffix('.json'));L=r['job']['length']
        assert r['job']['variant']=='native' and r['convergence']['verified'] and n.sha(p)==r['trace_sha256']
        dest=OUT/'checkpoints'/f'L{L:g}.npz';shutil.copy2(p,dest)
        inputs.append(dict(length=L,file=str(dest.relative_to(OUT)),source=str(p.relative_to(ROOT)),
            sha256=n.sha(dest),kind='prior_L18_continuation',convergence=r['convergence']))
    rule=copy.deepcopy(prior['steady_state']);rule['cap']=180000.
    m=dict(created_at=n.now(),sources=sources,inputs=inputs,case_template=prior['case_template'],
        steady_state=rule,original_cap=prior['steady_state']['cap'],checkpoint_sha256=n.sha(OUT/'start_L18.npz'),
        jobs=[dict(length=L,slot=i) for i,L in enumerate(NEW_LENGTHS)],
        definition='Native independent endpoints, 201 cells, dt=0.5 s; no cue. New lengths map exact polarized L18 concentrations by fractional coordinate and scale displacement/endpoints by L/18. Extend cap to 180000 s, retain original convergence tolerances and 10000-s holdout. No 0.05 Rac eligibility exclusion in these profile plots.')
    n.write(OUT/'manifest.json',m);print(len(NEW_LENGTHS),'new lengths')
def configure():
    n,e=modules();m=read(OUT/'manifest.json')
    for name,digest in m['sources'].items():assert n.sha(OUT/'frozen_sources'/name)==digest,name
    n.model_for=lambda key:n.ds.DeBelly(replace(n.db._parameters(),center_boundary_motion=False))
    return n,e,m
def worker(slot):
    n,e,m=configure();L=m['jobs'][slot]['length'];start=time.monotonic()
    cp=OUT/'start_L18.npz';assert n.sha(cp)==m['checkpoint_sha256']
    state,_=n.read_checkpoint(cp,'debelly')
    for k in ['displacement','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement']:state[k]*=L/18.
    case=dict(m['case_template'],length=L,width=.2*L)
    tr,convergence=e.release(case,state,m)
    p=OUT/'workers'/f'{slot:03d}.npz';n.save_checkpoint(p,'debelly',tr.state,tr)
    n.write(p.with_suffix('.json'),dict(length=L,slot=slot,convergence=convergence,
        final_orientation=float(tr.orientation[-1]),final_contrast=float(tr.contrast[-1]),
        trace_sha256=n.sha(p),wall_seconds=time.monotonic()-start))
    print(L,float(tr.orientation[-1]),convergence['verified'],convergence['elapsed'],flush=True)
def launch():
    n,e,m=configure();script=OUT/'frozen_sources'/Path(__file__).relative_to(ROOT)
    command=f'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; {PY} {script} worker --slot "$SLURM_ARRAY_TASK_ID"'
    args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=4G','--time=00:15:00',
        '--job-name=mca-threshold','--array=0-'+str(len(m['jobs'])-1),
        '--output='+str(OUT/'logs/%A_%a.out'),'--error='+str(OUT/'logs/%A_%a.err'),'--wrap',command]
    job=subprocess.check_output(args,text=True).strip();n.write(OUT/'launch.json',dict(array=job,generated_at=n.now()));print(job)
def audit_trace(n,tr,L,cv,scales,rule):
    weights=2*(np.arange(201)+.5)/201-1
    rac=tr.fields['rac'];ori=rac@weights/rac.sum(axis=1)
    np.testing.assert_allclose(ori,tr.orientation,atol=2e-15)
    np.testing.assert_allclose(np.ptp(rac,axis=1),tr.contrast,atol=2e-15)
    assert all(np.isfinite(f).all() for f in tr.fields.values())
    assert min(tr.fields[k].min() for k in ['rac','rho','mca_bound','mca_unbound'])>=0
    mass=(tr.fields['mca_bound']+tr.fields['mca_unbound']).sum(axis=1)
    mass_error=float(abs(mass/mass[0]-1).max());assert mass_error<1e-11
    # These checks use full arrays, independently of producer convergence flags.
    fields={k:tr.fields[k] for k in ['rac','rho','mca_bound','mca_unbound','velocity','tension']}
    fields['displacement']=tr.fields['displacement']-tr.fields['displacement'].mean(axis=1)[:,None]
    fields['physical_length']=L+tr.fields['boundary_separation']
    if cv['verified'] and len(tr.time)>1:
        end=tr.time[-1];before=end-rule['holdout'];worst=0.
        for stop in [before-rule['chunk'],before,end]:
            sel=(tr.time>=stop-rule['chunk'])&(tr.time<=stop);t=tr.time[sel]
            assert len(t)>2
            for k,f in fields.items():
                y=f[sel];worst=max(worst,float(np.ptp(y,axis=0).max()/scales[k]),
                    float(abs(y[-1]-y[-2]).max()/(t[-1]-t[-2])*rule['chunk']/scales[k]))
        i=int(np.argmin(abs(tr.time-before)));hold=tr.time>=before
        excursion=max(float(abs(f[hold]-f[i]).max()/scales[k]) for k,f in fields.items())
        assert worst<=rule['tolerance']*(1+1e-9) and excursion<=rule['holdout_tolerance']
    elif cv['verified']:
        # Original preparations retain compact final checkpoints. Their full
        # convergence was independently audited in the source campaign; here
        # verify the retained numerical convergence metrics, not absent arrays.
        assert read(BASE/'results/independent_audit.json')['passed']
        hold=cv['checks'][-1];assert hold['kind']=='holdout' and hold['verified']
        windows=cv['checks'][-3:-1];assert all(w['kind']=='convergence' for w in windows)
        worst=max(float(v) for record in [*windows,hold] for metrics in record['metrics'].values() for v in metrics.values())
        excursion=float(hold['maximum_scaled_excursion'])
        assert worst<=rule['tolerance']*(1+1e-9) and excursion<=rule['holdout_tolerance']
    else:worst=excursion=None
    b=tr.fields['mca_bound'][-1];u=tr.fields['mca_unbound'][-1]
    return dict(length_um=L,verified=cv['verified'],elapsed_s=float(tr.time[-1]),
        abs_rac_moment=float(abs(ori[-1])),rac_contrast=float(np.ptp(rac[-1])),
        bound_mca_moment=float(b@weights/b.sum()),unbound_mca_moment=float(u@weights/u.sum()),
        bound_mca_min=float(b.min()),bound_mca_max=float(b.max()),bound_mca_range=float(np.ptp(b)),
        unbound_mca_min=float(u.min()),unbound_mca_max=float(u.max()),
        field_window_max_scaled=worst,holdout_max_scaled=excursion,mca_mass_relative_error=mass_error)
def render():
    n,e,m=configure();inputs=copy.deepcopy(m['inputs'])
    for job in m['jobs']:
        p=OUT/'workers'/f'{job["slot"]:03d}.npz';r=read(p.with_suffix('.json'))
        inputs.append(dict(length=job['length'],file=str(p.relative_to(OUT)),sha256=r['trace_sha256'],
            kind='fresh_L18_continuation',convergence=r['convergence']))
    rows=[];profiles=[];states={};traces={}
    prep=read(BASE/'preparation/manifest.json')
    for inp in sorted(inputs,key=lambda r:r['length']):
        p=OUT/inp['file'];assert n.sha(p)==inp['sha256'];L=inp['length']
        state,tr=n.read_checkpoint(p,'debelly')
        scales=next(c['fixed_scales'] for c in prep['cases'] if c['length']==L) if inp['kind']=='original_preparation' else m['case_template']['fixed_scales']
        row=audit_trace(n,tr,L,inp['convergence'],scales,m['steady_state']);row.update(input_kind=inp['kind'],input_sha256=n.sha(p))
        rows.append(row);states[L]=(state,row);traces[L]=(tr.time,tr.orientation)
        for i,x in enumerate(tr.x):profiles.append(dict(length_um=L,x_nominal_um=float(x),x_fraction=float(x/L),
            bound_mca=float(state['mca_bound'][i]),unbound_mca=float(state['mca_unbound'][i]),
            rac=float(state['rac'][i]),rho=float(state['rho'][i]),verified=row['verified']))
    n.csvwrite(OUT/'results/length_summary.csv',rows);n.csvwrite(OUT/'results/spatial_profiles.csv',profiles)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    for category in ['vectors','data','metadata']:(FIG/'supporting'/category).mkdir(parents=True,exist_ok=True)
    cmap=plt.get_cmap('viridis');broad=[17.,17.75,18.,18.5,19.,21.,23.,25.,31.,39.]
    near=[17.,17.75,17.85,17.875,17.9,17.925,17.95,18.,18.1,18.25,18.5,19.]
    def style(L):
        _,r=states[L]
        return ':' if not r['verified'] else ('--' if r['bound_mca_range']<1e-4 else '-')
    def save(fig,stem):
        fig.savefig(FIG/(stem+'.png'),dpi=220)
        for ext in ['pdf','svg']:fig.savefig(FIG/'supporting/vectors'/(stem+'.'+ext))
        plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(13.2,8.1));q=(np.arange(201)+.5)/201
    for row,lengths in enumerate([broad,near]):
        norm=Normalize(min(lengths),max(lengths))
        for col,key in enumerate(['mca_bound','mca_unbound']):
            ax=axs[row,col]
            for L in lengths:
                state,r=states[L];label=f'{L:g} µm'+(' (unsettled)' if not r['verified'] else '')
                ax.plot(q,state[key],color=('#111111' if L==18 else cmap(norm(L))),ls=style(L),lw=(2.2 if L==18 else 1.7),label=label)
            basal=.1*(5 if col==0 else 1.4)/6.4
            ax.axhline(basal,color='#999999',lw=.8,ls=':')
            if col==0:
                ax.axhline(.068,color='#B34A75',lw=.9,ls='-.')
                ax.text(.98,.0684,'Rac-feedback midpoint: 0.068',ha='right',va='bottom',fontsize=8,color='#96355D')
            ax.set(xlabel='Fraction of nominal cell length (x/L)',ylabel=('Bound' if col==0 else 'Unbound')+' MCA density (native units)',xlim=(0,1))
            ax.set_title(f'{chr(65+2*row+col)}  '+('Across cell lengths' if row==0 else 'Near the lower polarity limit'),loc='left',weight='bold')
            ax.text(0,-.21,'Original Rac end',transform=ax.transAxes,fontsize=9)
            ax.text(1,-.21,'Original Rho end',transform=ax.transAxes,ha='right',fontsize=9)
            ax.grid(alpha=.13)
        axs[row,1].legend(loc='center left',bbox_to_anchor=(1.01,.5),frameon=False,fontsize=8,ncol=1)
    fig.suptitle('Steady MCA profiles across and below 18 µm',fontsize=16,weight='bold',y=.98)
    fig.subplots_adjust(left=.08,right=.84,top=.90,bottom=.17,hspace=.49,wspace=.28)
    fig.text(.08,.025,'Native model; no cue. Solid: verified asymmetric profiles. Dashed: verified near-homogeneous endpoints (overlapping). Gray dotted: basal density.\nOriginal polarized starts at 18–39 µm; finer lengths continued from the exact 18-µm state. No Rac-score eligibility filter applied.',fontsize=9)
    save(fig,'mca_threshold_profiles')
    fig,axs=plt.subplots(1,2,figsize=(12.4,4.5))
    for col,key in enumerate(['mca_bound','mca_unbound']):
        for L in broad:
            state,r=states[L];axs[col].plot(q*L,state[key],color=('#111111' if L==18 else cmap(Normalize(17,39)(L))),ls=style(L),label=f'{L:g} µm')
        axs[col].set(xlabel='Nominal position (µm)',ylabel=('Bound' if col==0 else 'Unbound')+' MCA density (native units)')
        axs[col].set_title('Bound MCA' if col==0 else 'Unbound MCA',loc='left');axs[col].grid(alpha=.13)
    axs[1].legend(loc='center left',bbox_to_anchor=(1,.5),frameon=False,fontsize=8)
    fig.suptitle('Steady MCA on each cell’s nominal physical domain',fontsize=14)
    fig.tight_layout(rect=(0,0,.9,.93));save(fig,'mca_threshold_physical_domains')
    fig,axs=plt.subplots(1,2,figsize=(12,4.6));small=[r for r in rows if r['length_um']<=19]
    xs=[r['length_um'] for r in small]
    axs[0].plot(xs,[abs(r['bound_mca_moment']) for r in small],'o-',color='#0072B2',label='Bound MCA')
    axs[0].plot(xs,[abs(r['unbound_mca_moment']) for r in small],'s-',color='#D55E00',label='Unbound MCA')
    axs[0].set(xlabel='Nominal length (µm)',ylabel='Absolute MCA first moment',title='MCA asymmetry near the lower limit');axs[0].legend(frameon=False,loc='lower right')
    axs[1].plot(xs,[r['abs_rac_moment'] for r in small],'o-',color='#333333')
    axs[1].axhline(.05,color='#888888',ls=':',label='Original Rac eligibility = 0.05')
    axs[1].set(xlabel='Nominal length (µm)',ylabel='Absolute Rac first moment',title='Polarity versus the original eligibility gate');axs[1].legend(frameon=False,fontsize=9,loc='lower right')
    homogeneous=[r['length_um'] for r in small if r['verified'] and r['bound_mca_range']<1e-4]
    asymmetric=[r['length_um'] for r in small if r['verified'] and r['bound_mca_range']>=1e-4]
    if homogeneous and asymmetric:
        lo,hi=max(homogeneous),min(asymmetric)
        for ax in axs:
            ax.axvspan(lo,hi,color='#AAAAAA',alpha=.2)
            ax.text(.03,.97,f'Sampled loss bracket: {lo:g}–{hi:g} µm',transform=ax.transAxes,va='top',fontsize=9)
    for ax in axs:
        for r in small:
            if not r['verified']:ax.axvline(r['length_um'],ls=':',color='#AA0000',lw=.8)
        ax.grid(alpha=.13)
    fig.tight_layout();save(fig,'mca_threshold_length_curve')
    for name in ['length_summary.csv','spatial_profiles.csv']:shutil.copy2(OUT/'results'/name,FIG/'supporting/data'/('mca_threshold_'+name))
    caption='Bound and unbound MCA surface-density profiles at verified native basal endpoints. All concentrations retain native units; no profile amplitude normalization. Left is the original Rac-rich end. Upper panels cover 17–39 µm; lower panels resolve the 17–19-µm range. Horizontal dotted lines mark homogeneous basal bound/unbound density. Solid curves are verified asymmetric states; dashed curves have bound-MCA range below 0.0001 native density. Any unverified result is dotted and explicitly labeled unsettled. Original L=18,19,21,23,25,31,39 checkpoints retain their audited preparations. L=17,17.5,17.75 reuse the previous native continuations; thirteen fresh lengths map the exact polarized L=18 chemical profiles by fractional coordinate and scale displacement/endpoints by L/18 before cue-free native evolution. Original two 1000-s convergence windows at 1e-4 frozen scales and 10000-s holdout at 1e-3 are retained; fresh continuation cap is 180000 s. No 0.05 Rac eligibility filter is applied: weaker persistent asymmetry remains visible. This is a sampled continuation family, not an exhaustive branch-existence proof. Native independent endpoints, 201 cells, step 0.5 s.'
    (FIG/'supporting/metadata/mca_threshold_profiles_caption.txt').write_text(caption+'\nGenerated '+n.now()+'; producer finished.\n')
    exports=[p for p in FIG.rglob('*') if p.is_file() and p.suffix in ['.png','.pdf','.svg','.csv','.txt']]
    n.write(FIG/'supporting/metadata/mca_threshold_profiles.manifest.json',dict(generated_at=n.now(),producer='finished',artifact_root='../..',sha256={str(p.relative_to(FIG)):n.sha(p) for p in exports},campaign=str(OUT.relative_to(ROOT)),render_script_sha256=n.sha(Path(__file__))))
    n.write(OUT/'results/audit.json',dict(generated_at=n.now(),producer='finished',complete=True,passed=True,all_profiles_verified=all(r['verified'] for r in rows),rows=rows,inputs=inputs,definition=m['definition']))
    print(json.dumps(rows,indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['init','launch','worker','render']);ap.add_argument('--slot',type=int)
    a=ap.parse_args()
    if a.action=='init':init()
    elif a.action=='launch':launch()
    elif a.action=='worker':worker(a.slot)
    else:render()
