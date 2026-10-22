"""Prepare/audit/launch a development-only corrected deBelly basal campaign."""
from pathlib import Path
from dataclasses import asdict
import argparse, json, os, shlex, shutil, subprocess, sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code/core'),str(Path(__file__).parent)]
import native as n
import length_dev_basal_feedback_worker as w
BASE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
DEFAULT=ROOT/'data/runs/debelly_basal_feedback_20261009'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
CONDITIONS=['native','equal_command','mca_feedback','tension_feedback']

def initialize(run):
    if run.exists():raise ValueError('Fresh campaign directory required')
    for name in ['checkpoints','frozen_sources','audits','workers','logs','results']:(run/name).mkdir(parents=True,exist_ok=True)
    source_manifest=json.loads((BASE/'manifest.json').read_text());provenance=json.loads((BASE/'provenance.json').read_text())
    frozen_verified={}
    for name in ['code/core/native.py','code/core/models/debelly.py','code/core/src/polarity1d/debelly.py','code/debelly_corrected_reversal.py']:
        digest=provenance['sources'][name]
        assert n.sha(ROOT/name)==digest and n.sha(BASE/'frozen_sources'/name)==digest
        frozen_verified[name]=digest
    jobs=[];cp_info={}
    for L in [23,31,39]:
        cases=[c for c in source_manifest['cases'] if c['model']=='debelly' and c['length']==L]
        assert len(cases)==3 and all(c['eligible'] and c['center_boundary_motion'] is False for c in cases)
        assert len({c['checkpoint_sha256'] for c in cases})==1 and len({c['checkpoint'] for c in cases})==1
        c=cases[0];cp=ROOT/c['checkpoint'];assert n.sha(cp)==c['checkpoint_sha256']
        shutil.copy2(cp,run/'checkpoints'/f'L{L}.npz');state,tr=n.read_checkpoint(cp,'debelly')
        cp_info[str(L)]=dict(source=c['checkpoint'],sha256=n.sha(cp),preparation_clock=float(tr.time[-1]),all_three_assays_same_checkpoint=True,eligible=True)
        for condition in CONDITIONS:
            jobs.append(dict(slot=len(jobs),length=L,cells=201,condition=condition,checkpoint_sha256=n.sha(cp),fixed_scales=c['fixed_scales']))
    native=ROOT/'code/core/models/debelly.py';text=native.read_text()
    marker='            old_left_boundary_displacement = left_boundary_displacement\n'
    assert text.count(marker)==1
    insertion='            # Development-only diagnostic command override; default preserves native code.\n            if BASAL_DIAGNOSTIC_LEFT_COMMAND is not None:\n                left_protrusion_velocity = float(BASAL_DIAGNOSTIC_LEFT_COMMAND)\n\n'
    diagnostic=text.replace(marker,insertion+marker)
    diagnostic+='\n# Sole solver alteration: an optional prescribed left protrusion command.\nBASAL_DIAGNOSTIC_LEFT_COMMAND = None\n'
    target=run/'frozen_sources/code/mechanisms/debelly_basal_feedback_native.py';target.parent.mkdir(parents=True,exist_ok=True);target.write_text(diagnostic)
    sources=n.sources()
    for name in ['code/mechanisms/length_dev_basal_feedback_worker.py','code/mechanisms/length_dev_basal_feedback_campaign.py']:
        sources[name]=n.sha(ROOT/name)
    for name,digest in sources.items():
        dest=run/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest);assert n.sha(dest)==digest
    sources['code/mechanisms/debelly_basal_feedback_native.py']=n.sha(target)
    # The scientific derivative is frozen only inside this development campaign.
    reference_state,_=n.read_checkpoint(run/'checkpoints/L31.npz','debelly')
    left,_=w.commands(reference_state,n.db._parameters(),n.db)
    rule=json.loads((BASE/'preparation/manifest.json').read_text())['steady_state']
    manifest=dict(created_at=n.now(),experiment='development_only_basal_feedback_controls',countercue=False,center_boundary_motion=False,
        length_units='µm',time_units='native seconds after intervention onset',conditions=CONDITIONS,jobs=jobs,checkpoints=cp_info,
        reference_length=31,reference_left_command=float(left),reference_feedback='L31 pre-cue native nonlinear increments interpolated in physical distance from the left boundary; constant extrapolation below its first cell centre',
        local_window=dict(plateau_um=[0,3],taper_um=[3,5],formula='w=1 for x<=3; w=0.5*(1+cos(pi*(x-3)/2)) for 3<x<5; w=0 for x>=5',application='blend nonlinear native and reference feedback increments, never transported fields'),
        endpoint_control='Left command fixed to initial L31 instantaneous signed command; native independent right command and 0.98 per 0.5 s relaxation retained',
        settling_rule=rule,settling_fields='native settling.monitored: Rac,Rho,bound/unbound MCA,mean-subtracted displacement,velocity,tension,physical length; fixed scales inherited per qualifying preparation; raw endpoints additionally retained and audited',
        diagnostic_native_change=dict(source='code/core/models/debelly.py',source_sha256=n.sha(native),command_override_insertion=insertion),
        canonical_sources_match_corrected_campaign=frozen_verified)
    n.write(run/'manifest.json',manifest)
    n.write(run/'provenance.json',dict(created_at=n.now(),frozen_sources=sources,manifest_sha256=n.sha(run/'manifest.json'),source_corrected_campaign=str(BASE.relative_to(ROOT)),source_manifest_sha256=n.sha(BASE/'manifest.json'),source_preparation_manifest_sha256=n.sha(BASE/'preparation/manifest.json')))
    print(run,flush=True)

def audit(run):
    manifest=w.read(run/'manifest.json');script=run/'frozen_sources/code/mechanisms/length_dev_basal_feedback_worker.py'
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for job in manifest['jobs']:
        subprocess.run([PYTHON,str(script),'audit','--run',str(run),'--slot',str(job['slot'])],env=env,check=True,stdout=subprocess.DEVNULL)
        subprocess.run([PYTHON,str(script),'short','--run',str(run),'--slot',str(job['slot'])],env=env,check=True,stdout=subprocess.DEVNULL)
    results=[w.read(run/'audits'/f"job_{j['slot']:02d}.json") for j in manifest['jobs']]
    assert all(r['initial_state_exact'] and r['left_endpoint_recurrence_passed'] and r['right_endpoint_recurrence_passed'] for r in results)
    assert all(r['native_copy_bit_identical'] for r in results if r['job']['condition']=='native')
    # On the reference length, all three perturbations coincide with native at onset.
    with np.load(run/'short_validation/04.npz') as z:native0={k:z['field_'+k][0] for k in ['rac_activation','rho_activation','displacement','tension']}
    for slot in [5,6,7]:
        with np.load(run/'short_validation'/f'{slot:02d}.npz') as z:
            for k in native0:np.testing.assert_array_equal(z['field_'+k][0],native0[k])
    n.write(run/'audits/summary.json',dict(passed=True,jobs=12,initial_reference_identity=True,
        maximum_rac_budget_residual=max(r['rac_budget_max_residual'] for r in results),
        maximum_short_mca_mass_error=max(w.read(run/'short_validation'/f"{j['slot']:02d}.json")['mca_mass_error'] for j in manifest['jobs']),
        native_unmodified_copy_bit_identical=True,tests=results))
    print(json.dumps(w.read(run/'audits/summary.json'),indent=2),flush=True)

def launch(run):
    assert w.read(run/'audits/summary.json')['passed']
    script=run/'frozen_sources/code/mechanisms/length_dev_basal_feedback_worker.py'
    command='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n'+shlex.join([PYTHON,str(script),'worker','--run',str(run),'--slot'])+' "$SLURM_ARRAY_TASK_ID"'
    args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=12:00:00','--array=0-11','--job-name=basal-feedback','--output='+str(run/'logs/%A_%a.out'),'--error='+str(run/'logs/%A_%a.err'),'--wrap',command]
    jobid=subprocess.check_output(args,text=True).strip()
    n.write(run/'launch.json',dict(submitted_at=n.now(),job_id=jobid,command=args,producer='running'))
    print(jobid,flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['init','audit','launch']);ap.add_argument('--run',type=Path,default=DEFAULT);a=ap.parse_args()
    {'init':initialize,'audit':audit,'launch':launch}[a.action](a.run.resolve())
if __name__=='__main__':main()
