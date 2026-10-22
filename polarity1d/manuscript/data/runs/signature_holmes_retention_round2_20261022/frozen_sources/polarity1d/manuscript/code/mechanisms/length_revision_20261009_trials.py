"""Bounded native replays for the requested length-mechanism figure revision."""
from pathlib import Path
import argparse, copy, json, shlex, shutil, subprocess, sys
import numpy as np

CODE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE))
import three_test_reversal as t
n, r, e = t.n, t.r, t.e
OUT = n.DATA/'runs/length_revision_20261009'


def init():
    assert not OUT.exists(), 'Use a fresh run destination'
    cases = r.read(t.DEFAULT/'manifest.json')['cases']
    jobs = []
    for model in ['goryachev', 'wave_pinning']:
        for length in [9., 21., 39.]:
            c = next(c for c in cases if c['model']==model and c['length']==length and c['condition']=='fixed_relative')
            jobs.append(dict(case=copy.deepcopy(c), duration=5120., recovery=0., purpose='half_mass'))
    c = next(c for c in cases if c['model']=='wave_pinning' and c['length']==15. and c['condition']=='fixed_physical')
    for duration in [20., 60.]:
        jobs.append(dict(case=copy.deepcopy(c), duration=duration, recovery=10000., purpose='withdrawal'))
    for folder in ['logs', 'workers', 'frozen_sources']:
        (OUT/folder).mkdir(parents=True, exist_ok=True)
    sources=n.sources()
    sources.update({str(p.relative_to(n.ROOT)):n.sha(p) for p in [Path(__file__), CODE/'three_test_reversal.py']})
    for path in sources:
        dst=OUT/'frozen_sources'/path
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(n.ROOT/path, dst)
    manifest=dict(jobs=jobs, sources=sources, interval=.5, generated_at=n.now(),
        protocol='Exact retained basal checkpoints; unchanged native kinetics and cues. Six relative-20% continuous-cue histories to 5120; two width-3 Wave-pinning pulses at L=15 followed by 10000 basal units. No fitted inputs.')
    r.write(OUT/'manifest.json',manifest)
    r.write(OUT/'status.json',dict(generated_at=n.now(), producer='running', complete=False))
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)), 'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00',
          '--array=0-7','--job-name=length-revision','--output='+str(OUT/'logs/%A_%a.out'),
          '--error='+str(OUT/'logs/%A_%a.err'),'--wrap',command]
    job=subprocess.check_output(args,text=True).strip()
    r.write(OUT/'launch.json',dict(array=job,generated_at=n.now()))
    print(job)


def worker(slot):
    m=r.read(OUT/'manifest.json')
    for path,digest in m['sources'].items(): assert n.sha(n.ROOT/path)==digest,path
    job=m['jobs'][slot];c=job['case']
    assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
    state,prep=n.read_checkpoint(n.ROOT/c['checkpoint'],c['model'])
    arrays={'originalstate':state.copy()}
    pulse=n.run(c['model'],c['length'],state,t.response_cue(c),job['duration'],interval=m['interval'],n=c['cells'])
    trajectories=[('pulse',pulse)]
    if job['recovery']:
        trajectories.append(('recovery',n.run(c['model'],c['length'],pulse.state,e.withdrawal(c),job['recovery'],interval=2.5,n=c['cells'])))
    for phase,tr in trajectories:
        for key in ['time','x','active','orientation','contrast']:
            arrays[phase+'_'+key]=getattr(tr,key)
        arrays.update({phase+'_field_'+key:value for key,value in tr.fields.items()})
    path=OUT/f'workers/{slot:02d}.npz'
    np.savez_compressed(path,**arrays)
    r.write(path.with_suffix('.json'),dict(job=job,archive_sha256=n.sha(path),generated_at=n.now()))
    print(slot,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker']);p.add_argument('--slot',type=int)
    args=p.parse_args()
    if args.action=='init':init()
    else:worker(args.slot)
