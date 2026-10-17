#!/usr/bin/env python3
"""Submit all fields and dependent audit/report in one invocation; record job IDs."""
import argparse,json,subprocess,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from common import BASE,REPO,profile_fields,write_json

def main():
    p=argparse.ArgumentParser();p.add_argument('--dry-run',action='store_true');a=p.parse_args()
    (BASE/'results/logs').mkdir(parents=True,exist_ok=True)
    command=['sbatch','--parsable',f'--array=0-{len(profile_fields("full"))-1}',str(BASE/'workflow/run_field.slurm'),'--profile','full']
    if a.dry_run:
        print(json.dumps({'array':command,'finalize':['sbatch','--parsable','--dependency=afterany:ARRAY_JOB_ID',str(BASE/'workflow/finalize.slurm')]},indent=2));return
    record=BASE/'results/full/submission.json'
    if record.exists():raise SystemExit('Submission already recorded; inspect job states before deliberately resubmitting fields and finalizer.')
    job=subprocess.check_output(command,cwd=REPO,text=True).strip().split(';')[0]
    info=dict(array_job=job,fields=profile_fields('full'),resources_per_field={'gpu_a30':1,'cpus':4,'memory_GiB':16,'time_limit_minutes':60},dependency='afterany',finalizer_job=None)
    write_json(record,info)
    final=subprocess.check_output(['sbatch','--parsable',f'--dependency=afterany:{job}',str(BASE/'workflow/finalize.slurm')],cwd=REPO,text=True).strip().split(';')[0]
    info['finalizer_job']=final;write_json(record,info);print(json.dumps(info,indent=2))
if __name__=='__main__':main()
