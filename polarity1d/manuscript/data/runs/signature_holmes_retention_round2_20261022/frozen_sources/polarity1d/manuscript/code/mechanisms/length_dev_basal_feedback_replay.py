"""Reproduce the exact frozen basal campaign in a fresh directory under SLURM."""
from pathlib import Path
import argparse,hashlib,json,shlex,shutil,subprocess
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def prepare(source,destination,prepare_only=False):
 if destination.exists():raise ValueError('A fresh destination is required')
 provenance=json.loads((source/'provenance.json').read_text());manifest=json.loads((source/'manifest.json').read_text())
 assert sha(source/'manifest.json')==provenance['manifest_sha256']
 for name,digest in provenance['frozen_sources'].items():assert sha(source/'frozen_sources'/name)==digest,name
 for job in manifest['jobs']:assert sha(source/'checkpoints'/f"L{job['length']}.npz")==job['checkpoint_sha256']
 destination.mkdir(parents=True)
 for name in ['frozen_sources','checkpoints','audits']:shutil.copytree(source/name,destination/name)
 for name in ['manifest.json','provenance.json']:shutil.copy2(source/name,destination/name)
 for name in ['workers','logs','results']:(destination/name).mkdir()
 command='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n'+shlex.join([PYTHON,str(destination/'frozen_sources/code/mechanisms/length_dev_basal_feedback_worker.py'),'worker','--run',str(destination),'--slot'])+' "$SLURM_ARRAY_TASK_ID"'
 args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=12:00:00','--array=0-11','--job-name=basal-feedback-replay','--output='+str(destination/'logs/%A_%a.out'),'--error='+str(destination/'logs/%A_%a.err'),'--wrap',command]
 record=dict(source=str(source),destination=str(destination),source_manifest_sha256=sha(source/'manifest.json'),frozen_sources_verified=True,checkpoint_hashes_verified=True,slurm_command=args,submitted=False)
 if not prepare_only:record.update(submitted=True,job_id=subprocess.check_output(args,text=True).strip())
 (destination/'replay.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--destination',type=Path,required=True);ap.add_argument('--prepare-only',action='store_true');a=ap.parse_args();prepare(a.source.resolve(),a.destination.resolve(),a.prepare_only)
