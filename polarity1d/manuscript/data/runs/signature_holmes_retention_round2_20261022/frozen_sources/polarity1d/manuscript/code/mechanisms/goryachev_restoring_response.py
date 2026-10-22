"""Small temporary incumbent activation suppression; native restoring response.

Fixed pre-cue half-height region, two pulse durations, three small amplitudes,
and matched sham controls. Native kinetics restored at release. Reaction and
diffusion contributions are recorded separately, with no opposite-end cue.
"""
from pathlib import Path
from dataclasses import dataclass
import argparse, copy, json, os, shlex, shutil, subprocess, sys, traceback
from types import SimpleNamespace
import numpy as np
CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE))
import three_test_reversal as t
g=SimpleNamespace(n=t.n,r=t.r,CODE=CODE,t=t)
from models.goryachev import Goryachev
from polarity1d.solvers import simulate_reaction_diffusion, laplacian_no_flux

OUT=g.n.DATA/'runs/goryachev_restoring_response_20261008'
FIG=g.n.ROOT/'figures/Goryachev_length_mechanism/diagnostics'
AMPLITUDES=[0.,.001,.005,.01]
DURATIONS=[1.,10.]


def core_width(state,L):
    u=state[0];q=(np.arange(len(u))+.5)/len(u)
    threshold=u.min()+.5*np.ptp(u);j=int(np.flatnonzero(u<threshold)[0]);assert j>0
    return float((q[j-1]+(threshold-u[j-1])/(u[j]-u[j-1])*(q[j]-q[j-1]))*L)


def weights(width,L,cells):
    dx=L/cells
    return np.clip(width-np.arange(cells)*dx,0.,dx)


@dataclass(frozen=True)
class Suppressed(Goryachev):
    suppression: float=0.
    width: float=1.
    length: float=9.

    def reaction(self,time,x,state,stimulus):
        u,v=state
        # Native finite-volume averaging of suppression in a partial edge cell.
        mask=weights(self.width,self.length,len(x))/(self.length/len(x))
        transfer=(self.a*self.ec*u*u*v+self.b*self.ec*u*v)*(1.-self.suppression*mask)-self.c*u
        transfer+=stimulus*v
        return np.vstack((transfer,-transfer))


def submit(out,action,array=None,dependency=None):
    command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)),action,'--out',str(out)])
    if array:command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    opts=['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00',
          '--job-name=gory-restoring-'+action,'--chdir='+str(g.n.ROOT),
          '--output='+str(out/'logs'/f'{action}_%A_%a.out'),'--error='+str(out/'logs'/f'{action}_%A_%a.err')]
    if array:opts+=['--array='+array]
    if dependency:opts+=['--dependency=afterany:'+dependency]
    return subprocess.check_output(opts+['--wrap',command],text=True).strip().split(';')[0]


def initialize(out):
    if out.exists():raise ValueError('Use a fresh destination')
    manifest=g.r.read(g.t.DEFAULT/'manifest.json');status=g.r.read(g.t.DEFAULT/'results/status.json')
    assert status['complete'] and status['execution_errors']==0
    cases=sorted([c for c in manifest['cases'] if c['model']=='goryachev' and c['condition']=='fixed_relative' and c['eligible'] and c['length']>=9],key=lambda c:c['length'])
    jobs=[]
    for c in cases:
        with np.load(g.n.ROOT/c['checkpoint']) as z:width=core_width(z['state'],c['length'])
        for D in DURATIONS:jobs.append(dict(slot=len(jobs),case=copy.deepcopy(c),width=width,duration=D,amplitudes=AMPLITUDES))
    for folder in ['logs','workers','results','frozen_sources']:(out/folder).mkdir(parents=True,exist_ok=True)
    m=dict(created_at=g.n.now(),jobs=jobs,lengths=[c['length'] for c in cases],durations=DURATIONS,amplitudes=AMPLITUDES,
       definition='During pulse only, multiply both native positive activation terms by 1-epsilon inside the fixed pre-cue half-height pole core. Partial edge-cell averaging. No added cue. Equal/opposite active/inactive transfers conserve material. Native kinetics and diffusion restored afterward.',
       measurement='Integrate gross activation, deactivation, net reaction, active diffusion and total active rate over the SAME fixed pre-cue core. Additional left-half measurements. Subtract matched sham at the same time; normalize post-release excess net reaction and total active rate by active core material lost at release.',
       recovery=5000.,solver=dict(backend='bdf',rtol=1e-10,atol=1e-12,max_step=10.),
       sampling='Pulse: 0 and 100 logarithmic samples 1e-5 to D. Recovery: 0 and 500 logarithmic samples 1e-5 to 5000. No temporal interpolation for plotted samples.',
       qualification='Tests the chosen local activation-suppression perturbation. It does not measure a global stability eigenvalue or resistance to finite reversal.')
    g.r.write(out/'manifest.json',m)
    hashes=g.n.sources()
    for p in [Path(__file__),g.CODE/'three_test_reversal.py',g.t.DEFAULT/'manifest.json']:
        hashes[str(p.relative_to(g.n.ROOT))]=g.n.sha(p)
    for name in hashes:
        dest=out/'frozen_sources'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(g.n.ROOT/name,dest)
    g.r.write(out/'provenance.json',dict(sources=hashes,manifest_sha256=g.n.sha(out/'manifest.json')))
    g.r.write(out/'results/status.json',dict(generated_at=g.n.now(),complete=False,producer='running',jobs=len(jobs)))
    array=submit(out,'worker',f'0-{len(jobs)-1}');final=submit(out,'finish',dependency=array)
    g.r.write(out/'launch.json',dict(generated_at=g.n.now(),array=array,finalizer=final))
    print(json.dumps(dict(jobs=len(jobs),trials=len(jobs)*len(AMPLITUDES),array=array,finalizer=final)))


def verify(out):
    p=g.r.read(out/'provenance.json');assert g.n.sha(out/'manifest.json')==p['manifest_sha256']
    for name,digest in p['sources'].items():
        assert g.n.sha(g.n.ROOT/name)==digest,'Source changed: '+name
        assert g.n.sha(out/'frozen_sources'/name)==digest
    return g.r.read(out/'manifest.json')


def run(model,state,L,times,max_step=10.,rtol=1e-10,atol=1e-12):
    return simulate_reaction_diffusion(model,initial_state=state,length=L,cells=state.shape[1],
          protocol=g.n.ConstantInput(0.),times=times,backend='bdf',rtol=rtol,atol=atol,max_step=max_step)


def observables(u,v,L,width,epsilon=0.):
    model=Goryachev();dx=L/u.shape[1];w=weights(width,L,u.shape[1]);mask=w/dx
    gross=(model.a*model.ec*u*u*v+model.b*model.ec*u*v)*(1.-epsilon*mask)
    death=model.c*u;net=gross-death;diff=model.du*laplacian_no_flux(u,dx)
    halfw=weights(L/2,L,u.shape[1]);arrays={}
    for region,ww in [('core',w),('half',halfw)]:
        for name,values in [('active',u),('activation',gross),('deactivation',death),('reaction',net),('diffusion',diff),('total_rate',net+diff)]:arrays[region+'_'+name]=values@ww
    return arrays


def worker(out,slot):
    job=verify(out)['jobs'][slot];c=job['case'];L=c['length'];D=job['duration'];width=job['width']
    assert g.n.sha(g.n.ROOT/c['checkpoint'])==c['checkpoint_sha256']
    with np.load(g.n.ROOT/c['checkpoint']) as z:original=z['state'].copy()
    directory=out/'workers'/f'{slot:03d}';directory.mkdir(exist_ok=True)
    for epsilon in job['amplitudes']:
        pulse=run(Suppressed(suppression=epsilon,width=width,length=L),original,L,np.r_[0.,np.geomspace(1e-5,D,100)])
        state=np.vstack([pulse.fields['active'][-1],pulse.fields['inactive'][-1]])
        recovery=run(Goryachev(),state,L,np.r_[0.,np.geomspace(1e-5,5000.,500)])
        arrays=dict(originalstate=original,x=pulse.x)
        for phase,tr in [('pulse',pulse),('recovery',recovery)]:
            u=tr.fields['active'];v=tr.fields['inactive'];mass=(u+v).sum(axis=1)
            assert np.max(abs(mass/original.sum()-1.))<1e-7 and min(u.min(),v.min())>=-1e-9
            arrays.update({phase+'_time':tr.time,phase+'_active':u,phase+'_inactive':v})
            arrays.update({phase+'_'+key:value for key,value in observables(u,v,L,width,epsilon if phase=='pulse' else 0.).items()})
        path=directory/f'epsilon_{epsilon:g}.npz';np.savez_compressed(path,**arrays)
        g.r.write(path.with_suffix('.json'),dict(generated_at=g.n.now(),job=job,epsilon=epsilon,archive_sha256=g.n.sha(path),
                    pulse_solver=pulse.metadata,recovery_solver=recovery.metadata))
    g.r.write(directory/'complete.json',dict(generated_at=g.n.now(),job=job,trials=len(job['amplitudes'])))


def sustained_return(time,relative,level):
    # Earliest sampled entry after which all later samples meet the tolerance.
    bad=np.flatnonzero(abs(relative)>level)
    i=(int(bad[-1])+1) if len(bad) else 0
    return float(time[i]) if i<len(time) else None


def finish(out):
    m=verify(out);rows=[];audits=[]
    for job in m['jobs']:
        c=job['case'];L=c['length'];directory=out/'workers'/f'{job["slot"]:03d}'
        assert g.r.read(directory/'complete.json')['job']==job
        with np.load(g.n.ROOT/c['checkpoint']) as cp:original=cp['state'].copy()
        with np.load(directory/'epsilon_0.npz') as z:sham={k:z[k].copy() for k in z.files}
        for eps in job['amplitudes']:
            path=directory/f'epsilon_{eps:g}.npz';metadata=g.r.read(path.with_suffix('.json'))
            assert metadata['job']==job and metadata['epsilon']==eps and g.n.sha(path)==metadata['archive_sha256']
            with np.load(path) as z:
                np.testing.assert_array_equal(z['originalstate'],original)
                np.testing.assert_allclose(z['pulse_active'][0],original[0],rtol=0,atol=1e-12)
                np.testing.assert_allclose(z['pulse_inactive'][0],original[1],rtol=0,atol=1e-12)
                for field in ['active','inactive']:np.testing.assert_allclose(z['recovery_'+field][0],z['pulse_'+field][-1],rtol=0,atol=1e-10)
                balance=0.
                # Independent integration of the instantaneous full RHS.
                dx=L/len(z['x']);w=np.clip(job['width']-np.arange(len(z['x']))*dx,0.,dx)
                for phase in ['pulse','recovery']:
                    u=z[phase+'_active'];v=z[phase+'_inactive'];assert np.isfinite(u).all() and np.isfinite(v).all()
                    assert min(u.min(),v.min())>=-1e-9
                    mass=(u+v).sum(axis=1);balance=max(balance,float(np.max(abs(mass/original.sum()-1.))))
                    activation=(.33/101*u*u*v+.67/101*u*v)*(1.-(eps if phase=='pulse' else 0.)*w/dx)
                    reaction=activation-.01733*u
                    lap=np.empty_like(u);lap[:,0]=(u[:,1]-u[:,0])/dx**2;lap[:,-1]=(u[:,-2]-u[:,-1])/dx**2
                    lap[:,1:-1]=(u[:,:-2]-2*u[:,1:-1]+u[:,2:])/dx**2
                    for name,expected in [('activation',activation@w),('deactivation',.01733*u@w),('reaction',reaction@w),('diffusion',.0025*lap@w),('total_rate',(reaction+.0025*lap)@w),('active',u@w)]:
                        np.testing.assert_allclose(z[phase+'_core_'+name],expected,rtol=1e-9,atol=1e-11)
                assert balance<1e-7
                audits.append(dict(archive=str(path.relative_to(out)),mass_error=balance,passed=True))
                if eps==0.:continue
                damage=float(sham['recovery_core_active'][0]-z['recovery_core_active'][0]);assert damage>0
                deficit=(sham['recovery_core_active']-z['recovery_core_active'])/damage
                reaction=z['recovery_core_reaction']-sham['recovery_core_reaction']
                total=z['recovery_core_total_rate']-sham['recovery_core_total_rate']
                initial_amount=float(sham['recovery_core_active'][0])
                rows.append(dict(length=L,duration=job['duration'],epsilon=eps,pole_width=job['width'],
                    core_active_material=initial_amount,active_damage_at_release=damage,fractional_damage=damage/initial_amount,
                    excess_net_activation_at_release=float(reaction[0]),reaction_gain_per_damage=float(reaction[0]/damage),
                    excess_total_material_rate_at_release=float(total[0]),total_gain_per_damage=float(total[0]/damage),
                    half_damage_recovery_time=sustained_return(z['recovery_time'],deficit,.5),
                    ninety_percent_damage_recovery_time=sustained_return(z['recovery_time'],deficit,.1),
                    final_deficit_over_initial_damage=float(deficit[-1]),
                    maximum_absolute_late_deficit_over_initial_damage=float(np.max(abs(deficit[z['recovery_time']>=1000.]))),
                    archive=str(path.relative_to(out))))
    g.n.csvwrite(out/'results/response_vs_length.csv',rows)
    g.r.write(out/'results/independent_audit.json',dict(generated_at=g.n.now(),passed=True,archives=len(audits),checks=audits))
    render(out,m,rows)
    g.r.write(out/'results/status.json',dict(generated_at=g.n.now(),producer='finished',complete=True,workers=len(m['jobs']),trials=len(audits),audit_passed=True))
    print(json.dumps(dict(trials=len(audits),complete=True)))


def render(out,m,rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    FIG.mkdir(parents=True,exist_ok=True)
    colors=['#0072B2','#D55E00','#009E73'];pdf=PdfPages(FIG/'restoring_response.pdf')
    for D in DURATIONS:
        fig,axes=plt.subplots(2,2,figsize=(12,9))
        specs=[('fractional_damage','Active core material lost (%)',100.),
               ('excess_net_activation_at_release','Excess net activation immediately after release\n(concentration × length / time)',1.),
               ('reaction_gain_per_damage','Excess net activation / active material lost\n(1 / native time)',1.),
               ('half_damage_recovery_time','Time to sustained recovery of half the damage',1.)]
        for ax,(key,label,factor) in zip(axes.flat,specs):
            for eps,color in zip(AMPLITUDES[1:],colors):
                rs=sorted([r for r in rows if r['duration']==D and r['epsilon']==eps],key=lambda r:r['length'])
                ax.plot([r['length'] for r in rs],[r[key]*factor if r[key] is not None else np.nan for r in rs],'o-',color=color,label=f'{100*eps:g}% suppression',ms=4)
            ax.set(xlabel='Cell length L',ylabel=label);ax.grid(alpha=.2)
        axes[0,0].legend(frameon=False,fontsize=10)
        fig.suptitle(f'Incumbent activation suppression for {D:g} native time unit(s); then native recovery',fontsize=14)
        fig.tight_layout(rect=[0,.03,1,.96]);fig.text(.07,.015,'Fixed half-height pole core. Responses subtract matched sham controls. Net activation excludes diffusion.',fontsize=10)
        fig.savefig(FIG/f'length_response_D{int(D)}.png',dpi=170);fig.savefig(FIG/f'length_response_D{int(D)}.svg');pdf.savefig(fig);plt.close(fig)
    for D in DURATIONS:
        fig,axes=plt.subplots(1,3,figsize=(15,5))
        for L,color in zip([9.,21.,39.],colors):
            row=next(r for r in rows if r['length']==L and r['duration']==D and r['epsilon']==.005)
            job=next(j for j in m['jobs'] if j['case']['length']==L and j['duration']==D)
            directory=out/'workers'/f'{job["slot"]:03d}'
            with np.load(out/row['archive']) as z,np.load(directory/'epsilon_0.npz') as sham:
                time=z['recovery_time'];damage=row['active_damage_at_release']
                values=[(sham['recovery_core_active']-z['recovery_core_active'])/damage,
                        (z['recovery_core_reaction']-sham['recovery_core_reaction'])/damage,
                        (z['recovery_core_total_rate']-sham['recovery_core_total_rate'])/damage]
                for ax,v in zip(axes,values):ax.plot(time,v,color=color,label=f'L={L:g}',lw=1.8)
        for ax,label in zip(axes,['Remaining active deficit / deficit at release','Excess net activation / initial deficit\n(1 / native time)','Excess material recovery rate / initial deficit\n(1 / native time)']):
            ax.set(xlabel='Time since native activation restored',ylabel=label,xscale='symlog',xlim=(0,1000.));ax.set_xscale('symlog',linthresh=.001);ax.grid(alpha=.2);ax.axhline(0,color='#777777',lw=.8)
        axes[0].axhline(.5,color='#777777',ls=':',lw=.8);axes[0].legend(frameon=False)
        fig.suptitle(f'0.5% activation suppression for {D:g} unit(s): native restoring trajectories')
        fig.tight_layout();fig.savefig(FIG/f'recovery_trajectories_D{int(D)}.png',dpi=170);pdf.savefig(fig);plt.close(fig)
    pdf.close();stamp=g.n.now()
    csvlink=os.path.relpath(out/'results/response_vs_length.csv',FIG)
    (FIG/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Goryachev incumbent restoring response</title></head><body style="font:17px system-ui;max-width:1400px;margin:30px auto"><h1>Incumbent restoring response</h1><p>Generated '+stamp+'; producer finished; independent audits passed. Temporary 0.1%, 0.5%, or 1% activation suppression in the fixed incumbent half-height core. Exact saved native starts; no countercue. Two pulse durations; native recovery afterward.</p><p><a href="restoring_response.pdf">Four-page PDF</a> · <a href="'+csvlink+'">Measurements</a></p>'+''.join('<img style="width:100%" src="'+s+'">' for s in ['length_response_D1.png','recovery_trajectories_D1.png','length_response_D10.png','recovery_trajectories_D10.png'])+'</body></html>')
    (out/'results/README.md').write_text('# Goryachev incumbent restoring response\n\nGenerated '+stamp+'; producer finished.\n\n'+m['definition']+'\n\n'+m['measurement']+'\n\nThe activation-suppression percentages are rate perturbations, not prescribed fractions of material removed. The damage and initial restoring response are measured at release. Sustained half-recovery is the first sampled time after which the absolute remaining core deficit stays below half its release value through the 5000-unit observation. Net reaction is activation minus deactivation; total material recovery additionally includes active diffusion across the fixed boundary. Matched sham subtraction removes the nonzero steady reaction/diffusion balance of the core. Tight BDF tolerances resolve the small response. All original inputs and checkpoints remain unchanged.\n\n'+m['qualification']+'\n')


def self_check():
    x=(np.arange(201)+.5)*39/201;state=np.vstack([np.linspace(100,.001,201),np.linspace(.04,.1,201)])
    m=Suppressed(length=39.,width=.8)
    np.testing.assert_array_equal(m.reaction(0,x,state,np.zeros(201)),Goryachev().reaction(0,x,state,np.zeros(201)))
    from dataclasses import replace
    for eps in AMPLITUDES:
        f=replace(m,suppression=eps).reaction(0,x,state,np.zeros(201));np.testing.assert_array_equal(f.sum(axis=0),np.zeros(201))
        mask=weights(.8,39.,201)/(39/201)
        np.testing.assert_array_equal(f[:,mask==0],Goryachev().reaction(0,x,state,np.zeros(201))[:,mask==0])
    print('Native identity, reaction conservation, fixed-core mask and untouched exterior checks passed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['init','worker','finish','self-check']);p.add_argument('--out',type=Path,default=OUT);p.add_argument('--slot',type=int)
    args=p.parse_args();out=args.out.resolve()
    try:
        if args.action=='self-check':self_check()
        elif args.action=='init':initialize(out)
        elif args.action=='worker':worker(out,args.slot)
        else:finish(out)
    except Exception:
        if out.exists():g.r.write(out/'results'/('worker_'+str(args.slot)+'.error.json' if args.action=='worker' else 'error.json'),dict(generated_at=g.n.now(),error=traceback.format_exc()))
        raise
