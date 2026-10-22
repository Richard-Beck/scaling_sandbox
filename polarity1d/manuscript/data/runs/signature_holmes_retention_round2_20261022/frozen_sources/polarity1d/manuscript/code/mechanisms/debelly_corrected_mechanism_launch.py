"""Launch fresh advection and audit/render both corrected mechanism figures."""
from pathlib import Path
import subprocess,sys,shlex,json
CODE=Path(__file__).resolve().parent
OUT=CODE.parents[1]/'data/runs/debelly_corrected_advection_20261008'
FIG=CODE.parents[1]/'figures/deBelly_mechanism'
PYTHON=sys.executable

def run(name):subprocess.run([PYTHON,str(CODE/name)],check=True)
def submit(name,command,dependency=None,array=None):
 opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00','--job-name='+name,'--output='+str(OUT/(name+'_%A_%a.out')),'--error='+str(OUT/(name+'_%A_%a.err'))]
 if dependency:opts+=['--dependency=afterok:'+dependency]
 if array:opts+=['--array='+array]
 return subprocess.check_output(opts+['--wrap','export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; '+command],text=True).strip().split(';')[0]
if sys.argv[1]=='launch':
 subprocess.run([PYTHON,str(CODE/'debelly_corrected_advection.py'),'setup'],check=True)
 array=submit('db-corrected-advection',shlex.join([PYTHON,str(CODE/'debelly_corrected_advection.py'),'worker'])+' "$SLURM_ARRAY_TASK_ID"',array='0-8')
 final=submit('db-corrected-mechanism-figures',shlex.join([PYTHON,str(Path(__file__)),'finish']),dependency=array+':20757917')
 (OUT/'launch.json').write_text(json.dumps(dict(array=array,finalizer=final),indent=2))
 print(array,final)
else:
 run('audit_debelly_corrected_advection.py');run('audit_debelly_corrected_zero_crossing.py')
 run('debelly_corrected_advection_figure.py');run('debelly_corrected_zero_crossing_figure.py')
 FIG.mkdir(parents=True,exist_ok=True)
 (FIG/'README.md').write_text('# Corrected deBelly mechanism figures\n\nBoth figures use the reflection-symmetric extension with independent local endpoint commands (`center_boundary_motion=False`). The old in-loop centered extension is flawed and deprecated. Inactive ends remain prescribed constraints; this is not a free-migration model.\n\n- `debelly_zero_crossing.png/pdf/svg`: MCA withdrawal landmark, Rac/Rho relaxation and spontaneous polarizability. An ineligible short-cell preparation is shown as establishment/withdrawal loss, not a reversal trial.\n- `debelly_advection.png/pdf/svg`: MCA polarity and nominal centroid speed, exact instantaneous prescribed-endpoint/direct-Rho/membrane-drag decomposition, physical Rho profiles and a fixed-profile length control.\n\nFresh corrected relative-20% patch/basal preparations start the qualifying reversal traces. Direct force superposition is not a full causal pathway partition. Native source budgets, complete arrays and independent audits are in ../../data/runs/debelly_corrected_{zero_crossing,advection,advection_revision}_20261008/.\n')
 subprocess.run([PYTHON,str(CODE.parent/'build.py'),'index'],check=True)
