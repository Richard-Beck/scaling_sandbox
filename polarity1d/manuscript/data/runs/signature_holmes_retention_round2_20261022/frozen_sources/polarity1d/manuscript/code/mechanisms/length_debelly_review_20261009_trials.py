"""Corrected native deBelly replays for the second review; immutable fresh inputs."""
from pathlib import Path
import argparse, copy, csv, gzip, json, shlex, shutil, subprocess, sys
import numpy as np

CODE=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(CODE),str(CODE/'mechanisms')]
import debelly_corrected_reversal as corrected
n,r=corrected.n,corrected.r
import three_test_mechanisms as mechanism
DEFAULT=n.DATA/'runs/length_debelly_review_20261009'


def init(out):
    assert not out.exists(),'Use a fresh destination'
    campaign=corrected.DEFAULT
    cases=r.read(campaign/'manifest.json')['cases']
    with gzip.open(campaign/'results/records.json.gz','rt') as f: records=json.load(f)['cases']
    jobs=[]
    for geo,lengths in [('fixed_slope',[18.,19.,21.,23.,25.,27.,29.,31.,33.,35.,37.,39.]),
                        ('fixed_physical',[25.,27.,29.,31.,33.,35.]),
                        ('fixed_relative',[25.,27.,29.,31.,33.,35.])]:
        for L in lengths:
            c=next(c for c in cases if c['model']=='debelly' and c['condition']==geo and c['length']==L)
            old=next(c for c in records if c['case_index']==c0) if False else next(a for a in records if a['model']=='debelly' and a['condition']==geo and a['length']==L)
            result=old['result']['orientation']
            doses=[('pass',result['endpoint'])]
            if geo!='fixed_slope' or L in [18,23,39]:doses.append(('fail',result['lower_failing_duration']))
            for label,duration in doses:
                expected=next(a for a in old['tests'] if a['duration']==duration)
                jobs.append(dict(case=copy.deepcopy(c),duration=duration,purpose=label,expected=expected))
    for L in [18.,19.,21.]:
        c=next(c for c in cases if c['model']=='debelly' and c['condition']=='fixed_slope' and c['length']==L)
        jobs.append(dict(case=copy.deepcopy(c),duration=1200.,purpose='continuous_and_zero_withdrawal'))
    for f in ['workers','logs','frozen_sources']:(out/f).mkdir(parents=True,exist_ok=True)
    sources=corrected.sources()
    for p in [Path(__file__),Path(mechanism.__file__),campaign/'manifest.json',campaign/'results/independent_audit.json']:
        sources[str(p.relative_to(n.ROOT))]=n.sha(p)
    for path in sources:
        dst=out/'frozen_sources'/path;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/path,dst)
    r.write(out/'manifest.json',dict(jobs=jobs,sources=sources,generated_at=n.now(),center_boundary_motion=False,
        protocol='Exact corrected basal checkpoints. Native selected/adjacent failing orientation pulses, full 10000-unit basal recovery. Extra continuous gradients at L=18,19,21 with recovery after the first sampled MCA zero. Native step0.5; enriched pulse sampling only. No cue/physics fitting.'))
    r.write(out/'status.json',dict(generated_at=n.now(),producer='running',complete=False,jobs=len(jobs)))
    cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)),'worker','--out',str(out)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00',
          '--array=0-'+str(len(jobs)-1),'--job-name=debelly-review','--output='+str(out/'logs/%A_%a.out'),
          '--error='+str(out/'logs/%A_%a.err'),'--wrap',cmd]
    array=subprocess.check_output(opts,text=True).strip()
    r.write(out/'launch.json',dict(array=array,jobs=len(jobs),generated_at=n.now()));print(array,len(jobs))


def worker(out,slot):
    m=r.read(out/'manifest.json')
    for path,digest in m['sources'].items():assert n.sha(n.ROOT/path)==digest,path
    job=m['jobs'][slot];c=job['case'];L=c['length']
    assert n.model_for('debelly').parameters.center_boundary_motion is False
    assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
    state,_=n.read_checkpoint(n.ROOT/c['checkpoint'],'debelly');original=copy.deepcopy(state)
    pulse=n.run('debelly',L,state,corrected.response_cue(c),job['duration'],interval=.5)
    end=pulse.state;D=job['duration'];arrays={'x':pulse.x}
    if job['purpose']=='continuous_and_zero_withdrawal':
        pm=pulse.fields['mca_bound']@(2*pulse.x/L-1)/pulse.fields['mca_bound'].sum(axis=1)
        ix=int(np.flatnonzero(pm<=0)[0]);D=float(pulse.time[ix])
        end={k:pulse.fields[k][ix].copy() for k in n.DB_NAMES}
        arrays['zero_time']=np.array(D)
    recovery=n.run('debelly',L,end,n.ConstantInput(0),10000.,interval=2.5)
    metrics=mechanism.gates(recovery,c)
    if 'expected' in job:
        for key in ['orientation_success','full_polarity_success','field_tail_pass']:
            assert metrics[key]==job['expected'][key],(slot,key,metrics,job['expected'])
        assert abs(metrics['final_orientation']-job['expected']['final_orientation'])<1e-9
    for phase,tr in [('pulse',pulse),('recovery',recovery)]:
        for key in ['time','orientation','contrast']:arrays[phase+'_'+key]=getattr(tr,key)
        arrays.update({phase+'_field_'+k:v for k,v in tr.fields.items() if k in n.DB_NAMES or k=='boundary_separation'})
    arrays.update({'originalstate_'+k:v for k,v in original.items()})
    path=out/f'workers/{slot:03d}.npz';np.savez_compressed(path,**arrays)
    r.write(path.with_suffix('.json'),dict(job=job,archive_sha256=n.sha(path),generated_at=n.now(),metrics=metrics,withdrawal_time=D))
    print(slot,c['condition'],L,job['purpose'],metrics,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker']);p.add_argument('--out',type=Path,default=DEFAULT);p.add_argument('--slot',type=int)
    a=p.parse_args();out=a.out.resolve()
    if a.action=='init':init(out)
    else:worker(out,a.slot)
