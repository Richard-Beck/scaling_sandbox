"""Separate frozen dt-halving controls at strongest input, after production."""
from pathlib import Path
import argparse,importlib.util,json,shutil,subprocess,sys
from datetime import datetime,timezone
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
PARENT=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026/noise'
OUT=PARENT/'strong_controls'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'

def prepare():
    if (OUT/'manifest.json').exists():raise ValueError('Already frozen')
    m=json.loads((PARENT/'manifest.json').read_text());jobs=[]
    for j in m['jobs']:
        if j['f']==.3 and j['C']==40 and j['replicate']<4 and not j.get('control'):
            jobs.append(dict(j,slot=len(jobs),dt=.025,control=True,parent_slot=j['slot']))
    if len(jobs)!=16:raise ValueError('Expected16controls')
    for part in ('frozen_sources','logs','results'):(OUT/part).mkdir(parents=True,exist_ok=True)
    shutil.copytree(PARENT/'inputs',OUT/'inputs')
    shutil.copy2(PARENT/'frozen_sources/round7_noise_20261026.py',OUT/'frozen_sources/round7_noise_20261026.py')
    shutil.copy2(__file__,OUT/'frozen_sources'/Path(__file__).name)
    m.update(jobs=jobs,created_utc=datetime.now(timezone.utc).isoformat(),purpose='Independent supplemental numerical controls at f0.3; not production replicates')
    (OUT/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
    script=OUT/'worker.sh';script.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+PY+' -B '+str(OUT/'frozen_sources'/Path(__file__).name)+' worker --slot "$SLURM_ARRAY_TASK_ID"\n');script.chmod(0o750)
    command=['sbatch','--parsable','--job-name=r7_dtstrong','--qos=small','--cpus-per-task=1','--mem=2G','--time=00:15:00','--array=0-15','--output='+str(OUT/'logs/%A_%a.out'),'--error='+str(OUT/'logs/%A_%a.err'),str(script)]
    r=subprocess.run(command,text=True,capture_output=True,check=True)
    (OUT/'submission.json').write_text(json.dumps(dict(job_id=r.stdout.strip(),command=command),indent=2)+'\n');print(r.stdout)

def worker(slot):
    p=OUT/'frozen_sources/round7_noise_20261026.py';spec=importlib.util.spec_from_file_location('strong_control_worker',p);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.RUN=OUT;mod.worker(slot)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','worker']);ap.add_argument('--slot',type=int);a=ap.parse_args()
    if a.command=='prepare':prepare()
    else:worker(a.slot)
