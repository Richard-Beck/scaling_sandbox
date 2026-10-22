"""Fresh short-domain continuation and native modal counterfactual diagnostics."""
from pathlib import Path
from dataclasses import replace, asdict
import argparse, json, os, shutil, subprocess, sys, time
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if p.name=='manuscript')
OUT = ROOT/'data/runs/debelly_short_length_diagnostic_20261010_v2'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
PY = '/home/4473331/.conda/envs/scaling_sandbox/bin/python'

def read(p): return json.loads(Path(p).read_text())

def modules(frozen=True):
    sys.path.insert(0, str((OUT/'frozen_sources' if frozen else ROOT)/'code/core'))
    import native as n
    import polarity_establishment as e
    return n, e

def init():
    n, e = modules(False)
    assert not OUT.exists(), 'Fresh destination required'
    for d in ['frozen_sources', 'workers', 'logs', 'results']:
        (OUT/d).mkdir(parents=True)
    hashes = n.sources()
    hashes[str(Path(__file__).relative_to(ROOT))] = n.sha(Path(__file__))
    for name, digest in hashes.items():
        dest = OUT/'frozen_sources'/name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/name, dest)
        assert n.sha(dest) == digest
    native = read(BASE/'manifest.json')
    c = next(c for c in native['cases'] if c['model']=='debelly' and c['length']==18)
    cp = ROOT/c['checkpoint']
    assert n.sha(cp)==c['checkpoint_sha256']
    shutil.copy2(cp, OUT/'start_L18.npz')
    prep = read(BASE/'preparation/manifest.json')
    case = next(c for c in prep['cases'] if c['length']==18)
    jobs = [dict(length=L, variant='native') for L in [16.,17.,17.5,17.75,18.]]
    jobs += [dict(length=17., variant=v) for v in ['half_gtpase_diffusion','half_cortical_viscosity']]
    n.write(OUT/'manifest.json', dict(created_at=n.now(), jobs=jobs, case_template=case,
            steady_state=prep['steady_state'], checkpoint_sha256=n.sha(cp), sources=hashes,
            start_definition='Exact native L18 polarized checkpoint; chemical concentrations mapped at equal fractional coordinates; displacement and endpoint displacements scaled by L/18. No cue; unchanged independent-endpoint rules.',
            interpretation='Continuation witnesses, not exhaustive basin search. Variant runs change only Rac/Rho diffusion or cortical viscosity. Original fixed convergence scales and holdout apply.'))
    print(OUT)

def configure():
    n, e = modules()
    m=read(OUT/'manifest.json')
    for name, digest in m['sources'].items():
        assert n.sha(OUT/'frozen_sources'/name)==digest, name
    return n, e, m

def parameters(n, variant):
    p=replace(n.db._parameters(), center_boundary_motion=False)
    if variant=='half_gtpase_diffusion': p=replace(p, rac_diffusion=.25, rho_diffusion=.25)
    elif variant=='half_cortical_viscosity': p=replace(p, cortical_viscosity=5000.)
    elif variant=='no_mca_to_rac': p=replace(p, rac_mechanical_feedback=0.)
    return p

def worker(slot):
    n,e,m=configure(); job=m['jobs'][slot]; L=job['length']; start=time.monotonic()
    p=parameters(n,job['variant'])
    n.model_for=lambda key:n.ds.DeBelly(p)
    cp=OUT/'start_L18.npz'; assert n.sha(cp)==m['checkpoint_sha256']
    state,_=n.read_checkpoint(cp,'debelly')
    for key in ['displacement','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement']:
        state[key]*=L/18.
    case=dict(m['case_template'],length=L,width=.2*L)
    tr,cv=e.release(case,state,m)
    n.save_checkpoint(OUT/'workers'/f'{slot}.npz','debelly',tr.state,tr)
    n.write(OUT/'workers'/f'{slot}.json',dict(job=job,parameters=asdict(p),convergence=cv,
            final_orientation=float(tr.orientation[-1]), final_contrast=float(tr.contrast[-1]),
            initial_orientation=float(tr.orientation[0]), wall_seconds=time.monotonic()-start,
            trace_sha256=n.sha(OUT/'workers'/f'{slot}.npz')))
    print(slot,job,float(tr.orientation[-1]),float(tr.contrast[-1]),cv['verified'],time.monotonic()-start)

def tangent(n,L,p,cells=201,step=.5,mode=1):
    x=n.rd.Grid1D(L,cells).x; co=np.cos(mode*np.pi*x/L); si=np.sin(mode*np.pi*x/L)
    R,H=n.ds._homogeneous_rac_rho(p,step); z=np.zeros(cells);o=np.ones(cells)
    b=dict(rac=o*R,rho=o*H,mca_bound=o*p.mca_density*p.mca_on/(p.mca_on+p.mca_off),
           mca_unbound=o*p.mca_density*p.mca_off/(p.mca_on+p.mca_off),
           **{k:z.copy() for k in ['displacement','velocity','tension','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement']})
    keys=['rac','rho','mca_bound','mca_unbound','displacement','velocity'];J=np.zeros((6,6))
    derivatives=[]
    for j,k in enumerate(keys):
        eps=1e-7 if k.startswith('mca') else 1e-6;ends=[]
        for sign in [-1,1]:
            s={k:v.copy() for k,v in b.items()}
            s[k]+=sign*eps*(si if k in ['displacement','velocity'] else co)
            tr=n.ds.simulate_debelly(n.ds.DeBelly(p),protocol=n.ConstantInput(0),length=L,cells=cells,
                    end_time=step,output_interval=step,step=step,initial_checkpoint=s)
            ends.append(tr.fields)
        d=[]
        for i,target in enumerate(keys):
            diff=(ends[1][target][-1]-ends[0][target][-1])/(2*eps)
            q=si if target in ['displacement','velocity'] else co
            J[i,j]=diff@q/(q@q);d.append(diff)
        derivatives.append(d)
    eig,V=np.linalg.eig(J);j=np.argmax(abs(eig));v=V[:,j];lam=eig[j]
    actual=np.einsum('j,jik->ik',v,np.asarray(derivatives))
    expected=np.array([lam*a*(si if k in ['displacement','velocity'] else co) for k,a in zip(keys,v)])
    return dict(length=L,cells=cells,step=step,mode=mode,rate=float(np.log(abs(lam))/step),
            modal_closure_residual=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected)),
            matrix=J.tolist())

def linear():
    n,e,m=configure();rows=[]
    for variant in ['native','half_gtpase_diffusion','half_cortical_viscosity','no_mca_to_rac']:
        for L in [12.,15.,16.,17.,18.,19.,19.5,19.55,20.,23.]:
            row=tangent(n,L,parameters(n,variant));row['variant']=variant;rows.append(row)
    for cells,step in [(101,.5),(401,.5),(201,.25)]:
        for L in [17.,19.5,19.55,23.]:
            row=tangent(n,L,parameters(n,'native'),cells,step);row['variant']='native';rows.append(row)
    for mode in [2,3,4]:
        row=tangent(n,17.,parameters(n,'native'),mode=mode);row['variant']='native';rows.append(row)
    n.write(OUT/'results/linear.json',dict(generated_at=n.now(),rows=rows))
    for row in rows: print({k:v for k,v in row.items() if k!='matrix'})

def launch():
    n,e,m=configure()
    script=OUT/'frozen_sources'/Path(__file__).relative_to(ROOT)
    args=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=2G',
          '--time=00:10:00','--job-name=debelly-short','--array=0-6',
          '--output='+str(OUT/'logs/%A_%a.out'),'--error='+str(OUT/'logs/%A_%a.err'),
          '--wrap',f'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; {PY} {script} worker --slot "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(args,text=True).strip()
    n.write(OUT/'launch.json',dict(array=job,generated_at=n.now()));print(job)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['init','launch','worker','linear']);ap.add_argument('--slot',type=int)
    a=ap.parse_args()
    if a.action=='init':init()
    elif a.action=='launch':launch()
    elif a.action=='linear':linear()
    else:worker(a.slot)
