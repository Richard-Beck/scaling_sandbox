"""Remove an inactive-profile fraction, evolve briefly, return the exact material."""
from pathlib import Path
import argparse, copy, json, shlex, shutil, subprocess, sys, traceback
import numpy as np

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE))
import three_test_reversal as t
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning
from polarity1d.solvers import simulate_reaction_diffusion, laplacian_no_flux
n=t.n
OUT=n.DATA/'runs/inactive_pool_restoring_response_20261008'
MODELS={'goryachev':Goryachev,'wave_pinning':WavePinning}
AMPLITUDES=[0.,.001,.005,.01]

def read(p):return json.loads(p.read_text())
def write(p,v):n.write(p,v)
def weights(width,L,N):return np.clip(width-np.arange(N)*L/N,0.,L/N)
def core_width(state,L):
    u=state[0];x=(np.arange(len(u))+.5)*L/len(u);h=u.min()+.5*np.ptp(u)
    j=int(np.flatnonzero(u<h)[0]);assert j>0
    return float(x[j-1]+(h-u[j-1])/(u[j]-u[j-1])*(x[j]-x[j-1]))

def submit(out,action,array=None,dependency=None):
    cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)),action,'--out',str(out)])
    if array:cmd+=' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00',
          '--job-name=inactive-return-'+action,'--chdir='+str(n.ROOT),'--output='+str(out/'logs'/f'{action}_%A_%a.out'),
          '--error='+str(out/'logs'/f'{action}_%A_%a.err')]
    if array:opts+=['--array='+array]
    if dependency:opts+=['--dependency=afterany:'+dependency]
    return subprocess.check_output(opts+['--wrap',cmd],text=True).strip().split(';')[0]

def initialize(out,durations):
    assert durations and min(durations)>0
    if out.exists():raise ValueError('Use a fresh destination')
    manifest=read(t.DEFAULT/'manifest.json');assert read(t.DEFAULT/'results/status.json')['complete']
    cases=sorted([c for c in manifest['cases'] if c['model'] in MODELS and c['condition']=='fixed_relative' and c['eligible'] and c['length']>=9],key=lambda c:(c['model'],c['length']))
    jobs=[]
    for c in cases:
        assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
        with np.load(n.ROOT/c['checkpoint']) as z:width=core_width(z['state'],c['length'])
        for D in durations:jobs.append(dict(slot=len(jobs),case=copy.deepcopy(c),width=width,duration=D,amplitudes=AMPLITUDES))
    for folder in ['logs','workers','results','frozen_sources']:(out/folder).mkdir(parents=True,exist_ok=True)
    m=dict(created_at=n.now(),jobs=jobs,durations=durations,amplitudes=AMPLITUDES,recovery=5000.,
        definition='At pulse onset remove epsilon times the original inactive profile pointwise throughout the whole cell into a nonreacting external reservoir. Active is unchanged. Evolve both fields under unchanged native reactions and diffusion for D. At release add back exactly the removed original profile, pointwise. Do not reset the evolved inactive field. Native recovery follows; no cue, clamp or kinetic modification.',
        measurement='Fixed original half-height incumbent core with fractional edge-cell integration. Subtract matched sham. Damage is active core material lost at return. Net reaction excludes diffusion; excess net reaction and total material rate are normalized by that damage.',
        accounting='Active+inactive+reservoir material conserved. Reservoir holds the removed amount during the pulse and is zero after return. Source spatial profile remains fixed while held.',
        solver=dict(backend='bdf',rtol=1e-10,atol=1e-12,max_step=10.),
        sampling='Pulse: 0 plus 100 logarithmic samples 1e-6 to D. Recovery: 0 plus 500 logarithmic samples 1e-5 to 5000.',
        zero_duration='An exact zero-duration removal/return is the identity: no active damage and undefined reaction-per-damage.')
    write(out/'manifest.json',m)
    hashes=n.sources()
    for path in [Path(__file__),CODE/'three_test_reversal.py',t.DEFAULT/'manifest.json']:
        hashes[str(path.relative_to(n.ROOT))]=n.sha(path)
    for name in hashes:
        dest=out/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/name,dest)
    write(out/'provenance.json',dict(sources=hashes,manifest_sha256=n.sha(out/'manifest.json')))
    write(out/'results/status.json',dict(generated_at=n.now(),producer='running',complete=False,jobs=len(jobs)))
    array=submit(out,'worker',f'0-{len(jobs)-1}');final=submit(out,'finish',dependency=array)
    write(out/'launch.json',dict(generated_at=n.now(),array=array,finalizer=final))
    print(json.dumps(dict(jobs=len(jobs),trials=len(jobs)*len(AMPLITUDES),array=array,finalizer=final)))

def verify(out):
    p=read(out/'provenance.json');assert n.sha(out/'manifest.json')==p['manifest_sha256']
    for name,digest in p['sources'].items():
        assert n.sha(n.ROOT/name)==digest,'Source changed: '+name
        assert n.sha(out/'frozen_sources'/name)==digest
    return read(out/'manifest.json')

def run(model,state,L,times):
    return simulate_reaction_diffusion(model,initial_state=state,length=L,cells=state.shape[1],protocol=n.ConstantInput(0.),
        times=times,backend='bdf',rtol=1e-10,atol=1e-12,max_step=10.)

def observables(model,u,v,L,width):
    if isinstance(model,Goryachev):gross=(model.a*model.ec*u*u+model.b*model.ec*u)*v;death=model.c*u
    else:gross=(model.k0+model.feedback*u*u/(model.half_saturation**2+u*u))*v;death=model.decay*u
    reaction=gross-death;diff=model.du*laplacian_no_flux(u,L/u.shape[1]);result={}
    for region,ww in [('core',weights(width,L,u.shape[1])),('half',weights(L/2,L,u.shape[1]))]:
        for name,value in [('active',u),('activation',gross),('deactivation',death),('reaction',reaction),('diffusion',diff),('total_rate',reaction+diff)]:result[region+'_'+name]=value@ww
    return result

def worker(out,slot):
    job=verify(out)['jobs'][slot];c=job['case'];L=c['length'];D=job['duration'];model=MODELS[c['model']]()
    assert n.sha(n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
    with np.load(n.ROOT/c['checkpoint']) as z:original=z['state'].copy()
    directory=out/'workers'/f'{slot:03d}';directory.mkdir(exist_ok=True)
    for epsilon in job['amplitudes']:
        removed=epsilon*original[1];reduced=original.copy();reduced[1]-=removed
        pulse=run(model,reduced,L,np.r_[0.,np.geomspace(1e-6,D,100)])
        returned=np.vstack([pulse.fields['active'][-1],pulse.fields['inactive'][-1]+removed])
        recovery=run(model,returned,L,np.r_[0.,np.geomspace(1e-5,5000.,500)])
        arrays=dict(originalstate=original,removed_inactive=removed,returnedstate=returned,x=pulse.x)
        for phase,tr in [('pulse',pulse),('recovery',recovery)]:
            u=tr.fields['active'];v=tr.fields['inactive'];mass=(u+v).sum(axis=1)+(removed.sum() if phase=='pulse' else 0.)
            assert np.max(abs(mass/original.sum()-1.))<1e-7 and min(u.min(),v.min())>=-1e-9
            arrays.update({phase+'_time':tr.time,phase+'_active':u,phase+'_inactive':v})
            arrays.update({phase+'_'+key:value for key,value in observables(model,u,v,L,job['width']).items()})
        path=directory/f'epsilon_{epsilon:g}.npz';np.savez_compressed(path,**arrays)
        write(path.with_suffix('.json'),dict(generated_at=n.now(),job=job,epsilon=epsilon,archive_sha256=n.sha(path),pulse_solver=pulse.metadata,recovery_solver=recovery.metadata))
    write(directory/'complete.json',dict(generated_at=n.now(),job=job,trials=len(job['amplitudes'])))

def sustained(time,relative,level):
    bad=np.flatnonzero(abs(relative)>level);i=int(bad[-1])+1 if len(bad) else 0
    return float(time[i]) if i<len(time) else None

def finish(out):
    m=verify(out);rows=[];audits=[]
    for job in m['jobs']:
        c=job['case'];L=c['length'];directory=out/'workers'/f'{job["slot"]:03d}'
        assert read(directory/'complete.json')['job']==job
        with np.load(n.ROOT/c['checkpoint']) as cp:original=cp['state'].copy()
        with np.load(directory/'epsilon_0.npz') as s:sham={key:s[key].copy() for key in s.files}
        for eps in job['amplitudes']:
            path=directory/f'epsilon_{eps:g}.npz';metadata=read(path.with_suffix('.json'))
            assert metadata['job']==job and metadata['epsilon']==eps and n.sha(path)==metadata['archive_sha256']
            with np.load(path) as z:
                np.testing.assert_array_equal(z['originalstate'],original)
                np.testing.assert_allclose(z['removed_inactive'],eps*original[1],rtol=0,atol=1e-13)
                np.testing.assert_allclose(z['pulse_active'][0],original[0],rtol=0,atol=1e-12)
                np.testing.assert_allclose(z['pulse_inactive'][0],(1.-eps)*original[1],rtol=0,atol=1e-12)
                np.testing.assert_allclose(z['recovery_active'][0],z['pulse_active'][-1],rtol=0,atol=1e-10)
                np.testing.assert_allclose(z['recovery_inactive'][0],z['pulse_inactive'][-1]+eps*original[1],rtol=0,atol=1e-10)
                balance=0.;dx=L/len(z['x']);w=np.clip(job['width']-np.arange(len(z['x']))*dx,0.,dx)
                for phase in ['pulse','recovery']:
                    u=z[phase+'_active'];v=z[phase+'_inactive'];assert np.isfinite(u).all() and np.isfinite(v).all()
                    assert min(u.min(),v.min())>=-1e-9
                    material=(u+v).sum(axis=1)+(eps*original[1].sum() if phase=='pulse' else 0.)
                    balance=max(balance,float(np.max(abs(material/original.sum()-1.))))
                    if c['model']=='goryachev':activation=(.33/101*u*u+.67/101*u)*v;death=.01733*u;du=.0025
                    else:activation=(.067+u*u/(1.+u*u))*v;death=u;du=.1
                    lap=np.empty_like(u);lap[:,0]=(u[:,1]-u[:,0])/dx**2;lap[:,-1]=(u[:,-2]-u[:,-1])/dx**2
                    lap[:,1:-1]=(u[:,:-2]-2*u[:,1:-1]+u[:,2:])/dx**2
                    for name,expected in [('active',u@w),('activation',activation@w),('deactivation',death@w),('reaction',(activation-death)@w),('diffusion',du*lap@w),('total_rate',(activation-death+du*lap)@w)]:
                        np.testing.assert_allclose(z[phase+'_core_'+name],expected,rtol=1e-9,atol=1e-11)
                assert balance<1e-7
                audits.append(dict(archive=str(path.relative_to(out)),material_plus_reservoir_error=balance,return_identity=True,observables_pass=True))
                if eps==0.:continue
                damage=float(sham['recovery_core_active'][0]-z['recovery_core_active'][0]);assert damage>0
                reaction=z['recovery_core_reaction']-sham['recovery_core_reaction'];total=z['recovery_core_total_rate']-sham['recovery_core_total_rate']
                deficit=(sham['recovery_core_active']-z['recovery_core_active'])/damage
                rows.append(dict(model=c['model'],length=L,duration=job['duration'],epsilon=eps,pole_width=job['width'],
                    core_active_material=float(sham['recovery_core_active'][0]),removed_inactive_amount=float(eps*original[1].sum()*dx),
                    active_damage_at_release=damage,fractional_damage=damage/float(sham['recovery_core_active'][0]),
                    excess_net_activation_at_release=float(reaction[0]),reaction_gain_per_damage=float(reaction[0]/damage),
                    excess_total_material_rate_at_release=float(total[0]),total_gain_per_damage=float(total[0]/damage),
                    half_damage_recovery_time=sustained(z['recovery_time'],deficit,.5),ninety_percent_damage_recovery_time=sustained(z['recovery_time'],deficit,.1),
                    final_deficit_over_initial_damage=float(deficit[-1]),archive=str(path.relative_to(out))))
    n.csvwrite(out/'results/response_vs_length.csv',rows)
    write(out/'results/independent_audit.json',dict(generated_at=n.now(),passed=True,archives=len(audits),checks=audits))
    write(out/'results/status.json',dict(generated_at=n.now(),producer='finished',complete=True,workers=len(m['jobs']),trials=len(audits),audit_passed=True))
    (out/'results/README.md').write_text('# Inactive-profile removal and return\n\n'+m['definition']+'\n\n'+m['measurement']+'\n\n'+m['accounting']+'\n')
    print(json.dumps(dict(complete=True,trials=len(audits))))

def self_check():
    for cls in MODELS.values():
        model=cls();state=np.vstack([np.linspace(1.,2.,201),np.linspace(.1,.08,201)])
        for eps in AMPLITUDES:
            removed=eps*state[1];reduced=state.copy();reduced[1]-=removed;restored=reduced.copy();restored[1]+=removed
            np.testing.assert_allclose(restored,state,rtol=0,atol=1e-15)
            assert abs(reduced.sum()+removed.sum()-state.sum())<1e-10
            x=(np.arange(201)+.5)*39/201
            r0=model.reaction(0,x,state,np.zeros(201));rr=model.reaction(0,x,reduced,np.zeros(201))
            # Both models' activation is linear in inactive concentration.
            death=model.c*state[0] if isinstance(model,Goryachev) else model.decay*state[0]
            np.testing.assert_allclose(rr[0]-r0[0],-eps*(r0[0]+death),rtol=1e-10,atol=1e-14)
    print('Exact zero-duration identity, reservoir accounting and inactive-linear reaction checks passed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker','finish','self-check']);p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--slot',type=int);p.add_argument('--durations',type=float,nargs='+')
    args=p.parse_args();out=args.out.resolve()
    try:
        if args.action=='init':initialize(out,args.durations)
        elif args.action=='worker':worker(out,args.slot)
        elif args.action=='finish':finish(out)
        else:self_check()
    except Exception:
        if out.exists():write(out/'results'/f'{args.action}_{args.slot}.error.json',dict(generated_at=n.now(),error=traceback.format_exc()))
        raise
