"""Frozen single-threshold protrusion ablation: establishment and branch survival."""
from pathlib import Path
import argparse, importlib.util, json, os, shlex, shutil, subprocess, sys, time
from dataclasses import replace, asdict
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'data/runs/debelly_no_protrusion_20261010'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
PYTHON = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'
SCRIPT = Path(__file__).resolve()

def read(p): return json.loads(Path(p).read_text())
def modules(run=None):
    sys.path.insert(0,str((run/'frozen_sources' if run else ROOT)/'code/core'))
    import native as n
    import polarity_establishment as e
    import settling
    return n,e,settling

def initialize(run):
    n,e,_=modules()
    assert not run.exists(), 'Fresh destination required'
    for folder in ['frozen_sources','checkpoints','workers','logs','results']:
        (run/folder).mkdir(parents=True,exist_ok=True)
    sources=n.sources()
    sources[str(SCRIPT.relative_to(ROOT))]=n.sha(SCRIPT)
    for name,digest in sources.items():
        dest=run/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,dest);assert n.sha(dest)==digest
    original=(ROOT/'code/core/models/debelly.py').read_text()
    old='            rac_threshold = 1.3 * rac0'
    new='            rac_threshold = float("inf")'
    assert original.count(old)==1
    variant=run/'frozen_sources/code/mechanisms/debelly_no_protrusion_native.py'
    variant.write_text(original.replace(old,new))
    manifest=read(BASE/'preparation/manifest.json')
    native=read(BASE/'manifest.json')
    assert read(BASE/'results/independent_audit.json')['passed']
    baseline=[];jobs=[]
    for case in manifest['cases']:
        L=case['length'];i=case['case_index']
        result=BASE/'preparation/work/cases'/f'case_{i:04d}/result.json'
        s=read(result)['summary']
        row=dict(length=L,summary=s,result=str(result.relative_to(ROOT)),result_sha256=n.sha(result))
        c=next(c for c in native['cases'] if c['model']=='debelly' and c['length']==L)
        if c['eligible']:
            cp=ROOT/c['checkpoint'];assert n.sha(cp)==c['checkpoint_sha256']
            dest=run/'checkpoints'/f'L{L:g}.npz';shutil.copy2(cp,dest)
            row.update(checkpoint=str(cp.relative_to(ROOT)),checkpoint_sha256=n.sha(cp))
        baseline.append(row)
        for duration in case['base_durations']:
            jobs.append(dict(slot=len(jobs),length=L,kind='from_rest',duration=duration,case=case))
        jobs.append(dict(slot=len(jobs),length=L,kind='small_seed',duration=0.,case=case))
        if c['eligible']:
            jobs.append(dict(slot=len(jobs),length=L,kind='native_branch',duration=0.,case=case,
                             checkpoint_sha256=c['checkpoint_sha256']))
    manifest.update(experiment='rac_protrusion_threshold_infinity',jobs=jobs,native_baseline=baseline,
                    intervention='Only rac_threshold = 1.3 * rac0 becomes rac_threshold = float("inf"); both protrusion switches are exactly zero. All parameters and other equations unchanged.',
                    center_boundary_motion=False,small_seed='Subtract 1e-5*cos(pi*x/L) from bound MCA and add it to unbound MCA; zero pointwise total-MCA change.',
                    interpretation='Existence of a verified asymmetric basal endpoint versus standard eligibility (left Rac moment magnitude >=0.05 and contrast >=0.05). Bounded preparation grid does not exclude all possible states.',
                    created_at=n.now())
    n.write(run/'manifest.json',manifest)
    n.write(run/'provenance.json',dict(sources=sources,variant_sha256=n.sha(variant),manifest_sha256=n.sha(run/'manifest.json'),
                                    source_change=dict(before=old,after=new)))
    n.write(run/'results/status.json',dict(generated_at=n.now(),producer='initialized',complete=False,total=len(jobs)))
    print('initialized',len(jobs),flush=True)

def configure(run):
    n,e,s=modules(run)
    provenance=read(run/'provenance.json')
    assert n.sha(run/'manifest.json')==provenance['manifest_sha256']
    for name,digest in provenance['sources'].items():assert n.sha(run/'frozen_sources'/name)==digest,name
    path=run/'frozen_sources/code/mechanisms/debelly_no_protrusion_native.py'
    assert n.sha(path)==provenance['variant_sha256']
    spec=importlib.util.spec_from_file_location('no_protrusion_native',path)
    ds=importlib.util.module_from_spec(spec);sys.modules[spec.name]=ds;spec.loader.exec_module(ds)
    n.ds=ds
    n.model_for=lambda key:ds.DeBelly(ds.DeBellyParameters(**asdict(replace(n.db._parameters(),center_boundary_motion=False))))
    return n,e,s

def worker(run,slot):
    n,e,_=configure(run);manifest=read(run/'manifest.json');job=manifest['jobs'][slot]
    start=time.monotonic();L=job['length'];case=job['case'];p=n.model_for('debelly').parameters
    assert p.center_boundary_motion is False and p.assay_mode=='reflection_symmetric_extension'
    assert asdict(p)==asdict(replace(n.db._parameters(),center_boundary_motion=False))
    state=None;pulse=None
    if job['kind']=='native_branch':
        cp=run/'checkpoints'/f'L{L:g}.npz';assert n.sha(cp)==job['checkpoint_sha256']
        state,_=n.read_checkpoint(cp,'debelly')
    elif job['kind']=='small_seed':
        initial=n.run('debelly',L,None,n.ConstantInput(0),.5,interval=.5,n=case['cells'])
        state={k:initial.fields[k][0].copy() for k in n.DB_NAMES}
        q=1e-5*np.cos(np.pi*n.rd.Grid1D(L,case['cells']).x/L)
        state['mca_bound']-=q;state['mca_unbound']+=q
    elif job['duration']>0:
        pulse=n.run('debelly',L,None,e.Patch(case['amplitude'],case['width'],width_mode='relative',fraction=.2),
                    job['duration'],interval=min(5.,job['duration']/20),n=case['cells'])
        state=pulse.state
    trace,convergence=e.release(case,state,manifest)
    assert all(np.isfinite(v).all() for v in trace.fields.values())
    prefix=run/'workers'/f'{slot:04d}'
    n.save_checkpoint(prefix.with_suffix('.npz'),'debelly',trace.state,trace)
    if pulse is not None:n.save_checkpoint(prefix.with_suffix('.pulse.npz'),'debelly',pulse.state,pulse)
    result=dict(generated_at=n.now(),job=job,convergence=convergence,final_orientation=float(trace.orientation[-1]),
                final_contrast=float(trace.contrast[-1]),parameters=asdict(p),wall_seconds=time.monotonic()-start,
                trace_sha256=n.sha(prefix.with_suffix('.npz')),pulse_sha256=n.sha(prefix.with_suffix('.pulse.npz')) if pulse else None)
    n.write(prefix.with_suffix('.json'),result)
    print(slot,L,job['kind'],job['duration'],result['final_orientation'],result['final_contrast'],convergence['verified'],result['wall_seconds'],flush=True)

def launch(run):
    manifest=read(run/'manifest.json')
    script=run/'frozen_sources'/SCRIPT.relative_to(ROOT)
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([PYTHON,str(script),'worker','--run',str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G','--time=00:20:00',
          '--job-name=no-protrusion','--array=0-'+str(len(manifest['jobs'])-1),'--chdir='+str(ROOT),
          '--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err')]
    jobid=subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]
    (run/'launch.json').write_text(json.dumps(dict(array=jobid,total=len(manifest['jobs'])),indent=2)+'\n')
    print(jobid)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['init','launch','worker']);ap.add_argument('--run',type=Path,default=DEFAULT);ap.add_argument('--slot',type=int)
    args=ap.parse_args();run=args.run.resolve()
    if args.action=='init':initialize(run)
    elif args.action=='launch':launch(run)
    else:worker(run,args.slot)
