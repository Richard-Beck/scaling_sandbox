"""Prepare immutable targeted Holmes campaign; integrate only in SLURM workers."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, time
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
OUT=ROOT/'polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022'
CODE=ROOT/'polarity1d/manuscript/code'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
def read(p):return json.loads(Path(p).read_text())
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def prepare():
 if OUT.exists():raise FileExistsError(OUT)
 old=ROOT/'polarity1d/manuscript/data/runs/signature_calibration_holmes_extended_20261020'
 calibrations=read(old/'results/calibrations.json')['calibrations']
 cal=next(r for r in calibrations if r['calibration_id']=='holmes_model3:L27')
 duration=cal['calibration']['upper_committing_duration'];assert duration==18125
 basals={r['basal_id']:r for r in read(old/'results/basals.json')['basals']}
 cases=[]
 for length in (31,35,39):
  for gap in (0,.1,.5,1,4,16):cases.append(dict(kind='train',length=length,gap=gap,seed=0,mean=0))
 # One S3-matched condition: 12 independent bounded-input histories plus mean/zero controls.
 for i in range(12):cases.append(dict(kind='noise',length=39,gap=16,seed=202610223900+i,mean=.02))
 cases.extend([dict(kind='mean',length=39,gap=16,seed=0,mean=.02),dict(kind='zero',length=39,gap=16,seed=0,mean=0)])
 for d in ('inputs','frozen_sources','workers','logs','results'):(OUT/d).mkdir(parents=True)
 for length in (31,35,39):
  b=basals[f'holmes_model3:L{length}'];source=Path(b['checkpoint']);digest=b['checkpoint_sha256'];assert sha(source)==digest
  shutil.copy2(source,OUT/f'inputs/L{length}.npz')
  for c in cases:
   if c['length']==length:c.update(checkpoint=f'inputs/L{length}.npz',checkpoint_sha256=digest,fixed_scales=b['fixed_scales'])
 for i,c in enumerate(cases):c.update(task=i,model='holmes_model3',native_length=float(c['length']),requested_reference_length=c['length'],cells=201,baseline_native_peak=1.,cue_width=3.8,cue_area=3.8,recovery_duration=10000.,calibration_duration=duration,pulse_count=4,total_on_factor_of_calibration=1.5,gap_multiple_of_pulse=c['gap'])
 shutil.copytree(CODE,OUT/'frozen_sources/polarity1d/manuscript/code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 for rel in ('data/inputs/fixed_scales.json','data/runs/debelly_corrected_three_test_20261008/manifest.json'):
  dest=OUT/'frozen_sources/polarity1d/manuscript'/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/'polarity1d/manuscript'/rel,dest)
 shutil.copy2(__file__,OUT/'worker.py')
 inputs={str(p.relative_to(OUT)):sha(p) for p in (OUT/'inputs').glob('*')}
 sources={str(p.relative_to(OUT)):sha(p) for p in (OUT/'frozen_sources').rglob('*') if p.is_file()}
 write(OUT/'manifest.json',dict(cases=cases,input_sha256=inputs,source_sha256=sources,worker_sha256=sha(OUT/'worker.py'),reference_calibration=dict(path=str(old/'results/calibrations.json'),sha256=sha(old/'results/calibrations.json'),duration=duration,pulse_duration=duration*1.5/4),protocol=dict(retention_primary='D and raw challenger share at every pulse/gap endpoint; no commitment calibration at target lengths',noise_condition='L39 gap16, mean .02, 12 seeds + mean and zero controls',noise_ell=1.9,noise_tau=10,noise_knot=2.5,noise_sigma=1,noise_transform='bounded_tanh',noise_hard_max=.04),status='prepared_not_submitted'))
 write(OUT/'freeze.json',dict(manifest_sha256=sha(OUT/'manifest.json')))
 worker=OUT/'worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+PYTHON+' -B '+str(OUT/'worker.py')+' worker --out '+str(OUT)+' --task "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
 command=['sbatch','--parsable','--partition=red','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=8G','--time=01:00:00','--array=0-31','--job-name=F2-Holmes-retention','--output='+str(OUT/'logs/%A_%a.out'),'--error='+str(OUT/'logs/%A_%a.err'),str(worker)]
 import shlex
 (OUT/'submission_command.txt').write_text(shlex.join(command)+'\n')
 # Actual submission is authorized. Record failures, never infer scheduler acceptance.
 result=subprocess.run(command,text=True,capture_output=True)
 write(OUT/'results/submission_attempt.json',dict(command=command,returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,submitted=result.returncode==0))
 print(json.dumps(dict(campaign=str(OUT),workers=len(cases),pulse_duration=duration*1.5/4,submission_returncode=result.returncode,stderr=result.stderr)))
def worker(out,task):
 out=Path(out);m=read(out/'manifest.json');assert sha(out/'manifest.json')==read(out/'freeze.json')['manifest_sha256']
 assert sha(out/'worker.py')==m['worker_sha256']
 for rel,digest in {**m['input_sha256'],**m['source_sha256']}.items():assert sha(out/rel)==digest,rel
 if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Native integrations require SLURM')
 sys.path.insert(0,str(out/'frozen_sources/polarity1d/manuscript/code/mechanisms'))
 import signature_campaign_20261020 as c
 import signature_noise_campaign_20261020 as nc
 from signature_analysis_20261020 import pulse_gap_endpoints
 from signature_noise_20261020 import LinearGapNoise
 native,inputs,readouts,_=c.modules();case=m['cases'][task];folder=out/'workers'/f'{task:04d}'
 if (folder/'result.json').exists():raise FileExistsError(folder)
 folder.mkdir(parents=True,exist_ok=True);started=time.monotonic();state,basal=native.read_checkpoint(out/case['checkpoint'],case['model'])
 schedule=inputs.make_train(case['calibration_duration'],1.5,4,case['gap'])
 forcing=None;forcing_audit=None
 if case['kind']!='train':
  latent=LinearGapNoise(case['native_length'],201,schedule.end_time,case['seed'],case['mean'],ell=1.9,tau=10,knot=2.5,sigma=1,kind=case['kind'],transform='bounded_tanh')
  forcing=nc.AppliedInputAudit(latent,1.);forcing_audit=latent.audit()
 trial=c._run_trial(native,inputs,readouts,case,state,basal,folder/'trial',schedule=schedule,campaign_root=out,gap_input=forcing,max_step=1.25 if forcing else None)
 measurement=read(folder/'trial/readouts.json');endpoints=pulse_gap_endpoints(measurement['series'],trial['events'])
 if forcing:
  forcing_audit['actual_queries']=forcing.audit();forcing_audit['hard_bound_pass']=bool(forcing.maximum<=2*case['mean']+1e-14);assert forcing_audit['hard_bound_pass']
 write(folder/'result.json',dict(case=case,trial=trial,endpoints=endpoints,forcing=forcing_audit,wall_seconds=time.monotonic()-started,manifest_sha256=sha(out/'manifest.json'),artifacts={str(p.relative_to(out)):sha(p) for p in (folder/'trial').rglob('*') if p.is_file()}))
 print(json.dumps(dict(task=task,status='complete')))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','worker']);p.add_argument('--out',type=Path,default=OUT);p.add_argument('--task',type=int);a=p.parse_args()
 if a.command=='prepare':prepare()
 else:worker(a.out,a.task)
