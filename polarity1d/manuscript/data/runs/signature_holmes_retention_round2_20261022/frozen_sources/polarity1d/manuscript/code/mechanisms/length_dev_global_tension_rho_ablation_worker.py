"""Global tension-to-Rho ablation only; original parameters and threshold retained.

Derived from length_dev_basal_feedback_worker.py. The frozen native solver is
unchanged; a feedback wrapper returns native delta-alpha and zero delta-beta.
"""
from pathlib import Path
from dataclasses import asdict
import argparse, importlib.util, json, sys, time
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'code/core'), str(Path(__file__).parent)]
import native as n
import settling

def read(p): return json.loads(Path(p).read_text())
def write(p,v): n.write(Path(p),v)
def sha(p): return n.sha(Path(p))

def load_solver(run):
    path=run/'frozen_sources/code/mechanisms/debelly_basal_feedback_native.py'
    spec=importlib.util.spec_from_file_location('basal_diagnostic_native',path)
    ds=importlib.util.module_from_spec(spec);sys.modules[spec.name]=ds;spec.loader.exec_module(ds)
    n.ds=ds
    def model_for(key):
        assert key=='debelly'
        values=asdict(n.db._parameters());values['center_boundary_motion']=False
        return ds.DeBelly(ds.DeBellyParameters(**values))
    n.model_for=model_for
    return ds

def window(x):
    x=np.asarray(x)
    return np.where(x<=3.,1.,np.where(x<5.,.5*(1+np.cos(np.pi*(x-3.)/2.)),0.))

def commands(fields,p,ds,override=None):
    R,_=ds._homogeneous_rac_rho(p,.5)
    left=p.protrusion_speed*ds.smooth_threshold(fields['rac'][...,0],p.protrusion_switch_scale,1.3*R)/(1+fields['tension'][...,0]**2/p.protrusion_tension_scale)
    right=-p.protrusion_speed*ds.smooth_threshold(fields['rac'][...,-1],p.protrusion_switch_scale,1.3*R)/(1+fields['tension'][...,-1]**2/p.protrusion_tension_scale)
    if override is not None: left=np.full_like(np.asarray(left),override)
    return left,right

def configure(ds,run,job,manifest):
    p=n.model_for('debelly').parameters
    assert job['condition']=='global_tension_ablation'
    assert p.rho_mechanical_feedback==.2 and p.center_boundary_motion is False
    original=ds._feedback_terms
    def feedback(parameters,mca_bound,tension):
        alpha,beta=original(parameters,mca_bound,tension)
        return alpha,np.zeros_like(beta)
    ds._feedback_terms=feedback
    ds.BASAL_DIAGNOSTIC_LEFT_COMMAND=None
    R,H=ds._homogeneous_rac_rho(p,.5)
    assert abs(1.3*R-.899061648345644)<1e-13
    return p,original,np.ones(job['cells']),np.zeros(job['cells']),np.zeros(job['cells'])

def laplacian(values,dx):
    flux=np.diff(values)/dx
    return np.diff(np.r_[0.,flux,0.])/dx

def audit_step(run,job,manifest):
    ds=load_solver(run);p,original,w,reference_da,reference_db=configure(ds,run,job,manifest)
    state,cp=n.read_checkpoint(run/'checkpoints'/f"L{int(job['length'])}.npz",'debelly')
    for k in n.DB_NAMES:np.testing.assert_array_equal(state[k],cp.fields[k][-1])
    native_da,native_db=original(p,state['mca_bound'],state['tension'])
    controlled_da,controlled_db=ds._feedback_terms(p,state['mca_bound'],state['tension'])
    frozen_rac=p.rac_deactivation*(controlled_da-native_da)/(1+state['rho']**2)
    frozen_rho=p.rho_deactivation*(controlled_db-native_db)/(1+state['rac']**2)
    np.testing.assert_array_equal(controlled_da,native_da)
    np.testing.assert_array_equal(controlled_db,np.zeros_like(native_db))
    R0,H0=ds._homogeneous_rac_rho(p,.5)
    Rnative,Hnative=n.db._homogeneous_rac_rho(n.db._parameters(),.5)
    assert R0==Rnative and H0==Hnative
    capture=[];implicit=ds._implicit_diffusion_decay
    def recording(values,source,diffusion,decay,step,dx):
        result=implicit(values,source,diffusion,decay,step,dx)
        capture.append(dict(old=values.copy(),source=source.copy(),diffusion=diffusion,decay=decay,step=step,dx=dx,new=result.copy()))
        return result
    ds._implicit_diffusion_decay=recording
    tr=n.run('debelly',job['length'],state,n.ConstantInput(0),.5,interval=.5,n=job['cells'])
    ds._implicit_diffusion_decay=implicit
    assert len(capture)==2
    R=capture[0];alpha=tr.fields['rac_activation'][1]
    activation=p.rac_deactivation*alpha
    inhibition=-activation*state['rho']**2/(1+state['rho']**2)
    decay=-R['decay']*R['new'];diffusion=R['diffusion']*laplacian(R['new'],R['dx'])
    rate=(R['new']-R['old'])/R['step'];residual=rate-activation-inhibition-decay-diffusion
    source_error=float(np.max(np.abs(R['source']-activation-inhibition)))
    assert source_error<1e-14 and np.max(np.abs(residual))<1e-13
    left,right=commands(state,p,ds,ds.BASAL_DIAGNOSTIC_LEFT_COMMAND)
    expected_left=p.boundary_relaxation*state['raw_left_displacement'][0]+.5*left
    expected_right=p.boundary_relaxation*state['raw_right_displacement'][0]+.5*right
    np.testing.assert_allclose(tr.state['raw_left_displacement'][0],expected_left,rtol=0,atol=1e-15)
    np.testing.assert_allclose(tr.state['raw_right_displacement'][0],expected_right,rtol=0,atol=1e-15)
    initial_mass=np.sum(state['mca_bound']+state['mca_unbound'])
    mass_error=float(np.max(np.abs(np.sum(tr.fields['mca_bound']+tr.fields['mca_unbound'],axis=1)/initial_mass-1)))
    assert mass_error<1e-12
    # A native copied solve must be bit-identical to the unmodified native solve.
    baseline_identical=None
    if job['condition']=='native':
        modified=n.ds;n.ds=n.db
        from dataclasses import replace
        n.model_for=lambda key:n.db.DeBelly(replace(n.db._parameters(),center_boundary_motion=False))
        native=n.run('debelly',job['length'],state,n.ConstantInput(0),.5,interval=.5,n=job['cells'])
        for k in tr.fields:np.testing.assert_array_equal(tr.fields[k],native.fields[k])
        baseline_identical=True;n.ds=modified
    arrays=dict(x=tr.x,global_ablation_mask=w,prescribed_delta_beta=reference_db,
        native_delta_alpha=native_da,native_delta_beta=native_db,controlled_delta_alpha=controlled_da,
        controlled_delta_beta=controlled_db,frozen_rac_derivative_change=frozen_rac,frozen_rho_derivative_change=frozen_rho,
        rac_activation=activation,rho_inhibition=inhibition,rac_deactivation=decay,rac_diffusion=diffusion,
        actual_rac_rate=rate,rac_budget_residual=residual)
    np.savez_compressed(run/'audits'/f"job_{job['slot']:02d}.npz",**arrays)
    result=dict(job=job,initial_state_exact=True,native_copy_bit_identical=baseline_identical,
        single_global_tension_to_rho_feedback_removal=True,original_homogeneous_rac=R0,original_homogeneous_rho=H0,original_rac_threshold=1.3*R0,direct_tension_inhibition_retained=True,native_rho_basal_activation_retained=p.rho_basal_activation==1.,rac_budget_max_residual=float(abs(residual).max()),
        rac_source_max_residual=source_error,mca_mass_error=mass_error,
        left_command=float(left),right_command=float(right),left_endpoint_recurrence_passed=True,right_endpoint_recurrence_passed=True,
        frozen_left_rac_derivative_change=float(frozen_rac[0]),frozen_left_rho_derivative_change=float(frozen_rho[0]))
    write(run/'audits'/f"job_{job['slot']:02d}.json",result)
    return result

def integrate(run,slot,short=None):
    manifest=read(run/'manifest.json');job=manifest['jobs'][slot]
    for name,digest in read(run/'provenance.json')['frozen_sources'].items():assert sha(run/'frozen_sources'/name)==digest,name
    assert sha(run/'checkpoints'/f"L{int(job['length'])}.npz")==job['checkpoint_sha256']
    ds=load_solver(run);p,_,_,_,_=configure(ds,run,job,manifest)
    state,_=n.read_checkpoint(run/'checkpoints'/f"L{int(job['length'])}.npz",'debelly')
    start=time.monotonic();rule=manifest['settling_rule'];scales=job['fixed_scales'];checks=[];trace=None
    elapsed=0.;consecutive=0;verified=False
    # Fine early sampling resolves the endpoint relaxation and the staggered response.
    early=n.run('debelly',job['length'],state,n.ConstantInput(0),10.,interval=.5,n=job['cells'])
    for k in n.DB_NAMES:np.testing.assert_array_equal(early.fields[k][0],state[k])
    trace=early;state=early.state;elapsed=10.
    if short is None:
        while elapsed+rule['chunk']<=rule['cap']:
            tr=n.run('debelly',job['length'],state,n.ConstantInput(0),rule['chunk'],interval=rule['sample_interval'],n=job['cells'])
            state=tr.state
            assert all(np.isfinite(v).all() for v in tr.fields.values())
            metric=settling.diagnostics(tr,'debelly',job['length'],scales,rule)
            trace=n.concatenate(trace,tr,elapsed);elapsed+=rule['chunk']
            checks.append(dict(time=elapsed,kind='convergence',**metric));consecutive=consecutive+1 if metric['steady'] else 0
            write(run/'workers'/f'{slot:02d}.status.json',dict(producer='running',elapsed=elapsed,consecutive=consecutive,latest=checks[-1]))
            if consecutive<rule['consecutive_windows']:continue
            if elapsed+rule['holdout']>rule['cap']:break
            target={k:v[-1].copy() for k,v in settling.monitored(trace,'debelly',job['length']).items()}
            hold=n.run('debelly',job['length'],state,n.ConstantInput(0),rule['holdout'],interval=rule['sample_interval'],n=job['cells'])
            state=hold.state
            excursion=float(settling.deviations(settling.monitored(hold,'debelly',job['length']),target,scales).max())
            final=settling.diagnostics(hold,'debelly',job['length'],scales,rule)
            verified=excursion<=rule['holdout_tolerance'] and final['steady']
            trace=n.concatenate(trace,hold,elapsed);elapsed+=rule['holdout']
            checks.append(dict(time=elapsed,kind='holdout',maximum_scaled_excursion=excursion,verified=bool(verified),**final))
            if verified:break
            consecutive=0
    mass=np.sum(trace.fields['mca_bound']+trace.fields['mca_unbound'],axis=1)
    mass_error=float(np.max(abs(mass/mass[0]-1)));assert mass_error<1e-10
    left,right=commands(trace.fields,p,ds,ds.BASAL_DIAGNOSTIC_LEFT_COMMAND)
    folder=run/('short_validation' if short is not None else 'workers');folder.mkdir(exist_ok=True)
    n.save_checkpoint(folder/f'{slot:02d}.npz','debelly',state,trace)
    np.savez_compressed(folder/f'{slot:02d}.commands.npz',time=trace.time,left_command=left,right_command=right)
    result=dict(job=job,producer='finished',verified_steady_state=bool(verified),censored_at_cap=bool(short is None and not verified),elapsed=elapsed,wall_seconds=time.monotonic()-start,settling_checks=checks,mca_mass_error=mass_error,initial_state_exact=True,
        final_orientation=float(trace.orientation[-1]),final_contrast=float(trace.contrast[-1]),retained_polarity=bool(trace.orientation[-1]<-.05 and trace.contrast[-1]>=.05),parameters=asdict(p),trace_sha256=sha(folder/f'{slot:02d}.npz'))
    write(folder/f'{slot:02d}.json',result)
    print(json.dumps({k:result[k] for k in ['job','verified_steady_state','censored_at_cap','elapsed','wall_seconds','mca_mass_error']}),flush=True)
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['audit','worker','short']);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--slot',type=int,required=True);a=ap.parse_args()
    if a.action=='audit':print(json.dumps(audit_step(a.run,read(a.run/'manifest.json')['jobs'][a.slot],read(a.run/'manifest.json')),indent=2))
    else:integrate(a.run,a.slot,10 if a.action=='short' else None)
if __name__=='__main__':main()
