"""One global tension→Rho ablation at three exact corrected basal checkpoints."""
from pathlib import Path
import argparse,json,os,shlex,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'code/core'),str(Path(__file__).parent)]
import native as n
import length_dev_global_tension_rho_ablation_worker as worker
PARENT=ROOT/'data/runs/debelly_basal_feedback_20261009_v2'
DEFAULT=ROOT/'data/runs/debelly_global_tension_rho_ablation_20261009'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
def initialize(run):
 if run.exists():raise ValueError('Fresh destination required')
 for name in ['checkpoints','frozen_sources','audits','workers','logs','results']:(run/name).mkdir(parents=True,exist_ok=True)
 parent=json.loads((PARENT/'manifest.json').read_text());provenance=json.loads((PARENT/'provenance.json').read_text())
 assert n.sha(PARENT/'manifest.json')==provenance['manifest_sha256']
 assert json.loads((PARENT/'results/independent_audit.json').read_text())['passed']
 jobs=[];controls=[]
 for L,controlslot in [(23,0),(31,4),(39,8)]:
  old=parent['jobs'][controlslot];assert old['condition']=='native' and old['length']==L
  cp=PARENT/'checkpoints'/f'L{L}.npz';assert n.sha(cp)==old['checkpoint_sha256'];shutil.copy2(cp,run/'checkpoints'/cp.name)
  resultpath=PARENT/'workers'/f'{controlslot:02d}.json';tracepath=resultpath.with_suffix('.npz');result=json.loads(resultpath.read_text())
  assert result['verified_steady_state'] and result['initial_state_exact'] and n.sha(tracepath)==result['trace_sha256']
  controls.append(dict(length=L,trace=str(tracepath.relative_to(ROOT)),trace_sha256=n.sha(tracepath),result=str(resultpath.relative_to(ROOT)),result_sha256=n.sha(resultpath),commands=str((PARENT/'workers'/f'{controlslot:02d}.commands.npz').relative_to(ROOT)),commands_sha256=n.sha(PARENT/'workers'/f'{controlslot:02d}.commands.npz'),checkpoint_sha256=old['checkpoint_sha256']))
  jobs.append(dict(old,slot=len(jobs),condition='global_tension_ablation'))
 sources=n.sources()
 # All scientific dependencies are identical to the already-verified controls.
 for name,digest in sources.items():assert provenance['frozen_sources'][name]==digest,name
 for name in ['code/mechanisms/length_dev_global_tension_rho_ablation_worker.py','code/mechanisms/length_dev_global_tension_rho_ablation_campaign.py']:sources[name]=n.sha(ROOT/name)
 for name,digest in sources.items():
  dest=run/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest);assert n.sha(dest)==digest
 native=ROOT/'code/core/models/debelly.py'
 diagnostic=run/'frozen_sources/code/mechanisms/debelly_basal_feedback_native.py'
 # Inactive marker only: the integration source itself is unchanged native code.
 diagnostic.write_text(native.read_text()+'\nBASAL_DIAGNOSTIC_LEFT_COMMAND = None\n')
 sources[str(diagnostic.relative_to(run/'frozen_sources'))]=n.sha(diagnostic)
 manifest=dict(created_at=n.now(),experiment='single_global_tension_to_rho_ablation',jobs=jobs,conditions=['global_tension_ablation'],native_controls_reused=controls,
  checkpoints=parent['checkpoints'],settling_rule=parent['settling_rule'],center_boundary_motion=False,countercue=False,
  intervention='Return original native delta_alpha_mca and exactly zero delta_beta_tension at every cell; original parameters remain unchanged, including rho_mechanical_feedback=0.2 used for the original homogeneous reference.',
  retained='Native basal Rho activation, mutual Rac/Rho inhibition, MCA→Rac feedback, Rho active stress, direct tension inhibition of both protrusion commands, original endpoint relaxation, all MCA transport/exchange and membrane/cortex mechanics.',
  original_rac_threshold=.899061648345644,original_homogeneous_rac=.899061648345644/1.3,reference_left_command=parent['reference_left_command'],source_parent=str(PARENT.relative_to(ROOT)))
 n.write(run/'manifest.json',manifest);n.write(run/'provenance.json',dict(created_at=n.now(),frozen_sources=sources,manifest_sha256=n.sha(run/'manifest.json'),parent_manifest_sha256=n.sha(PARENT/'manifest.json'),native_control_source_identity=True))
 print(run)
def audit(run):
 m=json.loads((run/'manifest.json').read_text());script=run/'frozen_sources/code/mechanisms/length_dev_global_tension_rho_ablation_worker.py'
 env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
 for j in m['jobs']:
  for action in ['audit','short']:subprocess.run([PYTHON,str(script),action,'--run',str(run),'--slot',str(j['slot'])],check=True,env=env,stdout=subprocess.DEVNULL)
 results=[json.loads((run/'audits'/f"job_{j['slot']:02d}.json").read_text()) for j in m['jobs']]
 assert all(r['single_global_tension_to_rho_feedback_removal'] and r['direct_tension_inhibition_retained'] and r['original_rac_threshold']==m['original_rac_threshold'] for r in results)
 n.write(run/'audits/summary.json',dict(passed=True,tests=results,maximum_rac_budget_residual=max(r['rac_budget_max_residual'] for r in results),original_homogeneous_reference_unchanged=True))
 print(json.dumps(json.loads((run/'audits/summary.json').read_text()),indent=2))
def launch(run):
 assert json.loads((run/'audits/summary.json').read_text())['passed']
 script=run/'frozen_sources/code/mechanisms/length_dev_global_tension_rho_ablation_worker.py'
 command='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n'+shlex.join([PYTHON,str(script),'worker','--run',str(run),'--slot'])+' "$SLURM_ARRAY_TASK_ID"'
 args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=12:00:00','--array=0-2','--job-name=tension-rho-ablation','--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err'),'--wrap',command]
 jobid=subprocess.check_output(args,text=True).strip();n.write(run/'launch.json',dict(submitted_at=n.now(),job_id=jobid,command=args));print(jobid)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['init','audit','launch']);ap.add_argument('--run',type=Path,default=DEFAULT);a=ap.parse_args();{'init':initialize,'audit':audit,'launch':launch}[a.action](a.run.resolve())
