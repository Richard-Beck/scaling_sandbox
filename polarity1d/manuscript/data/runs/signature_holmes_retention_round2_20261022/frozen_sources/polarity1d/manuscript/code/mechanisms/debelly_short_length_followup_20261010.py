"""Attachment-diffusion control for the frozen short-length diagnostic."""
from pathlib import Path
from dataclasses import replace, asdict
import argparse, json, subprocess, sys, time
ROOT=next(p for p in Path(__file__).resolve().parents if p.name=='manuscript')
OUT=ROOT/'data/runs/debelly_short_length_diagnostic_20261010_v2'
sys.path.insert(0,str(OUT/'frozen_sources/code/mechanisms'))
import debelly_short_length_diagnostic_20261010 as s

def linear():
    n,e,m=s.configure();rows=[]
    variants=[('half_unbound_mca_diffusion',dict(mca_unbound_diffusion=.15)),
              ('half_both_mca_diffusion',dict(mca_unbound_diffusion=.15,mca_bound_diffusion=.0005)),
              ('no_gtpase_diffusion',dict(rac_diffusion=0.,rho_diffusion=0.))]
    for label,kw in variants:
        p=replace(s.parameters(n,'native'),**kw)
        for L in [12.,15.,17.,18.,19.5,23.]:
            row=s.tangent(n,L,p);row['variant']=label;rows.append(row)
    from scipy.optimize import brentq
    roots=[]
    for label,p in [('native',s.parameters(n,'native')),
                    ('half_gtpase_diffusion',s.parameters(n,'half_gtpase_diffusion')),
                    ('half_cortical_viscosity',s.parameters(n,'half_cortical_viscosity')),
                    ('half_unbound_mca_diffusion',replace(s.parameters(n,'native'),mca_unbound_diffusion=.15))]:
        root=brentq(lambda L:s.tangent(n,L,p)['rate'],10,25,xtol=1e-5)
        roots.append(dict(variant=label,linear_onset_um=root))
    n.write(OUT/'results/linear_additional.json',dict(generated_at=n.now(),script_sha256=n.sha(Path(__file__)),
            parameter_variants='Only stated diffusion coefficients change',rows=rows,roots=roots))
    print(roots)

def worker(slot):
    n,e,m=s.configure();L=[17.,15.][slot];start=time.monotonic()
    p=replace(s.parameters(n,'native'),mca_unbound_diffusion=.15)
    n.model_for=lambda key:n.ds.DeBelly(p)
    cp=OUT/'start_L18.npz';assert n.sha(cp)==m['checkpoint_sha256']
    state,_=n.read_checkpoint(cp,'debelly')
    for key in ['displacement','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement']:
        state[key]*=L/18.
    case=dict(m['case_template'],length=L,width=.2*L)
    tr,cv=e.release(case,state,m)
    prefix=OUT/'workers'/f'mca_{slot}'
    n.save_checkpoint(prefix.with_suffix('.npz'),'debelly',tr.state,tr)
    n.write(prefix.with_suffix('.json'),dict(job=dict(length=L,variant='half_unbound_mca_diffusion'),
            parameters=asdict(p),convergence=cv,final_orientation=float(tr.orientation[-1]),
            final_contrast=float(tr.contrast[-1]),initial_orientation=float(tr.orientation[0]),
            wall_seconds=time.monotonic()-start,trace_sha256=n.sha(prefix.with_suffix('.npz')),
            script_sha256=n.sha(Path(__file__))))
    print(L,float(tr.orientation[-1]),float(tr.contrast[-1]),cv['verified'])

def launch():
    import shutil
    dest=OUT/'frozen_sources/code/mechanisms'/Path(__file__).name
    assert not dest.exists()
    shutil.copy2(Path(__file__),dest)
    args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G',
          '--time=00:10:00','--job-name=debelly-mca-short','--array=0-1',
          '--output='+str(OUT/'logs/mca_%A_%a.out'),'--error='+str(OUT/'logs/mca_%A_%a.err'),
          '--wrap',f'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; {s.PY} {dest} worker --slot "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(args,text=True).strip();n,e,m=s.configure()
    n.write(OUT/'mca_launch.json',dict(array=job,script_sha256=n.sha(dest),generated_at=n.now()));print(job)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['linear','worker','launch']);ap.add_argument('--slot',type=int)
    a=ap.parse_args()
    if a.action=='linear':linear()
    elif a.action=='launch':launch()
    else:worker(a.slot)
