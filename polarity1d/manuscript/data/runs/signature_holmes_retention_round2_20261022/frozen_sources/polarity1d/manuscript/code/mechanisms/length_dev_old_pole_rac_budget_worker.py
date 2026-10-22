"""Observe exact staggered Rac budgets in three unchanged native reversals."""
from pathlib import Path
from dataclasses import replace, asdict
import argparse, importlib.util, json, shlex, shutil, subprocess, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'data/runs/debelly_old_pole_rac_budget_20261010'
PARENT=ROOT/'data/runs/length_narrative_memory_20261009'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
def read(p):return json.loads(Path(p).read_text())
def initialize(run):
    sys.path.insert(0,str(ROOT/'code/core'));import native as n
    assert not run.exists()
    for folder in ['workers','logs','frozen_sources']:(run/folder).mkdir(parents=True)
    m=read(PARENT/'manifest.json');assert read(PARENT/'status.json')['audit_passed']
    jobs=[]
    for L in [23,31,39]:
        j=next(j for j in m['jobs'] if j['case']['length']==L and j['case']['condition']=='fixed_relative')
        source=ROOT/j['source'] if j['source'] else PARENT/f"workers/{m['missing_jobs'].index(j):03d}.npz"
        assert n.sha(source)==read(source.with_suffix('.json'))['archive_sha256']
        cp=ROOT/j['case']['checkpoint'];assert n.sha(cp)==j['case']['checkpoint_sha256']
        jobs.append(dict(slot=len(jobs),length_um=L,case=j['case'],duration=j['duration'],limit_s=1200.,
                         source=str(source),source_sha256=n.sha(source),checkpoint=str(cp),checkpoint_sha256=n.sha(cp)))
    sources=n.sources();sources[str(Path(__file__).relative_to(ROOT))]=n.sha(Path(__file__))
    for name,digest in sources.items():
        dest=run/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest);assert n.sha(dest)==digest
    n.write(run/'manifest.json',dict(created_at=n.now(),jobs=jobs,center_boundary_motion=False,
        method='Unmodified corrected native dynamics; observational wrappers capture each staggered half-second Rac solve, including fractional final pulse steps. Absolute and additive baseline-subtracted budgets; no local cue term plotted because exactly zero.'))
    n.write(run/'provenance.json',dict(sources=sources,manifest_sha256=n.sha(run/'manifest.json')))
    script=run/'frozen_sources'/Path(__file__).relative_to(ROOT)
    cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([PYTHON,str(script),'worker','--run',str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job=subprocess.check_output(['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G','--time=00:10:00',
        '--array=0-2','--job-name=old-pole-budget','--chdir='+str(ROOT),'--output='+str(run/'logs/%A_%a.out'),
        '--error='+str(run/'logs/%A_%a.err'),'--wrap',cmd],text=True).strip().split(';')[0]
    n.write(run/'launch.json',dict(array=job,generated_at=n.now()));print(job,flush=True)

def worker(run,slot):
    sys.path.insert(0,str(run/'frozen_sources/code/core'))
    import native as n
    import polarity_establishment as e
    provenance=read(run/'provenance.json');assert n.sha(run/'manifest.json')==provenance['manifest_sha256']
    for name,digest in provenance['sources'].items():assert n.sha(run/'frozen_sources'/name)==digest
    job=read(run/'manifest.json')['jobs'][slot];L=job['length_um'];D=job['duration'];case=job['case']
    assert n.sha(job['source'])==job['source_sha256'] and n.sha(job['checkpoint'])==job['checkpoint_sha256']
    ds=n.db;p=replace(ds._parameters(),center_boundary_motion=False)
    n.model_for=lambda key:ds.DeBelly(p)
    initial,_=n.read_checkpoint(job['checkpoint'],'debelly')
    implicit=ds._implicit_diffusion_decay;derived=ds._derived_rates
    records=[];context={};clock=0.;calls=0;phase='cue';oldrho=initial['rho'].copy()
    maximum_residual=0.;maximum_source_error=0.;maximum_cue=0.
    def rates(parameters,mca,tension,external):
        a,b=derived(parameters,mca,tension,external)
        context.update(mca=mca.copy(),external=np.asarray(external).copy())
        return a,b
    def solve(values,source,diffusion,decay,step,dx):
        nonlocal clock,calls,oldrho,maximum_residual,maximum_source_error,maximum_cue
        new=implicit(values,source,diffusion,decay,step,dx)
        if calls%2==0:
            assert diffusion==p.rac_deactivation*p.rac_diffusion and decay==p.rac_deactivation
            factor=1/(1+oldrho**2)
            switch=.5*(1-np.tanh((context['mca']-p.mca_rac_threshold)/p.mca_switch_scale))
            basal=p.rac_deactivation*p.rac_basal_activation*factor
            mca=p.rac_deactivation*p.rac_mechanical_feedback*p.rho_to_rac_coupling*switch*factor
            cue=p.rac_deactivation*context['external']*factor
            loss=-decay*new
            transport=diffusion*np.diff(np.r_[0.,np.diff(new)/dx,0.])/dx
            actual=(new-values)/step;total=basal+mca+cue+loss+transport
            maximum_residual=max(maximum_residual,float(abs(actual-total).max()))
            maximum_source_error=max(maximum_source_error,float(abs(source-basal-mca-cue).max()))
            maximum_cue=max(maximum_cue,float(abs(cue[0])))
            clock+=step
            records.append(dict(time_s=clock,step_s=step,phase=phase,basal_activation=float(basal[0]),
                mca_activation=float(mca[0]),deactivation=float(loss[0]),diffusion=float(transport[0]),
                net_rate=float(actual[0]),component_sum=float(total[0]),rho_suppression_factor=float(factor[0]),
                old_rho=float(oldrho[0]),updated_bound_mca=float(context['mca'][0]),mca_switch=float(switch[0]),
                local_cue_activation=float(cue[0]),rac_before=float(values[0]),rac_after=float(new[0])))
        else:oldrho=new.copy()
        calls+=1
        return new
    ds._derived_rates=rates;ds._implicit_diffusion_decay=solve
    cue=e.Patch(case['response_amplitude'],case['width'],side='right',width_mode=case['width_mode'],fraction=.2)
    pulse=n.run('debelly',L,initial,cue,D,interval=.5,n=case['cells'])
    assert abs(clock-D)<1e-8
    phase='basal_recovery';clock=D
    recovery=n.run('debelly',L,pulse.state,n.ConstantInput(0),job['limit_s']-D,interval=.5,n=case['cells'])
    ds._derived_rates=derived;ds._implicit_diffusion_decay=implicit
    assert maximum_residual<1e-12 and maximum_source_error<1e-14 and maximum_cue==0
    # The recorder must preserve every retained state at all common sample times.
    maximum_replay_error=0.;compared=0
    with np.load(job['source']) as z:
        for label,tr in [('pulse',pulse),('recovery',recovery)]:
            t=z[label+'_time'];mask=t<=tr.time[-1]+1e-9;times=t[mask]
            indices=np.searchsorted(tr.time,times)
            indices=np.minimum(indices,len(tr.time)-1)
            assert np.max(abs(tr.time[indices]-times))<1e-8
            for key in n.DB_NAMES:
                expected=z[label+'_field_'+key][mask]
                delta=float(abs(tr.fields[key][indices]-expected).max())
                maximum_replay_error=max(maximum_replay_error,delta)
                np.testing.assert_allclose(tr.fields[key][indices],expected,rtol=0,atol=1e-10)
            compared+=len(times)
    prefix=run/'workers'/f'{slot:03d}'
    arrays={k:np.array([r[k] for r in records]) for k in records[0]}
    arrays.update(x=pulse.x,time_rac=np.r_[pulse.time,D+recovery.time[1:]],
                  old_pole_rac=np.r_[pulse.fields['rac'][:,0],recovery.fields['rac'][1:,0]],
                  initial_old_pole_rac=np.array(initial['rac'][0]))
    np.savez_compressed(prefix.with_suffix('.npz'),**arrays)
    n.write(prefix.with_suffix('.json'),dict(generated_at=n.now(),job=job,passed=True,parameters=asdict(p),
        source_sha256=job['source_sha256'],archive_sha256=n.sha(prefix.with_suffix('.npz')),
        maximum_all_cell_budget_residual=maximum_residual,maximum_all_cell_activation_source_residual=maximum_source_error,
        maximum_old_pole_cue_activation=maximum_cue,maximum_native_replay_field_error=maximum_replay_error,
        original_sample_times_compared=compared,recorded_native_Rac_steps=len(records)))
    print(L,'passed',len(records),maximum_residual,maximum_replay_error,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['launch','worker']);ap.add_argument('--run',type=Path,default=RUN);ap.add_argument('--slot',type=int)
    a=ap.parse_args();initialize(a.run.resolve()) if a.action=='launch' else worker(a.run.resolve(),a.slot)
