"""Frozen deterministic edge pulses and untuned linear response, Figure 3 round 5.
All native integrations require SLURM. No stochastic realization is generated.
"""
from pathlib import Path
import os,sys,json,hashlib,time,shutil,subprocess,argparse
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np
from scipy.integrate import solve_ivp
from scipy import sparse
from scipy.sparse.linalg import eigs
from scipy.optimize._numdiff import approx_derivative
import figure4_round3_campaign_20261022 as source
ROOT=source.ROOT; MS=source.MS; CODE=source.CODE
ORIGINAL=source.DEFAULT
RUN=MS/'data/runs/figure3_round5_edge_pulses_v2_20261023'
OUT=MS/'figures/length_mechanisms_20261008/dev/revision_round5_20261023/figure3'
FRACTIONS=(.0001,.0003,.001,.003,.005,.01)
REF=source.REFERENCES

def write(p,v):source.write(p,v)
def read(p):return source.read(p)
def sha(p):return source.sha(p)
def now():return source.now()

def edges(u,C):
    """Linearly interpolated rising/falling current half-contrast crossings."""
    u=np.asarray(u,dtype=float); n=len(u); level=(u.min()+u.max())/2
    if np.ptp(u)<1e-6:return None
    high=u>=level
    rising=np.flatnonzero(high & ~np.roll(high,1)); falling=np.flatnonzero(~high & np.roll(high,1))
    if len(rising)!=1 or len(falling)!=1:return None
    answer=[]
    for j in (rising[0],falling[0]):
        i=(j-1)%n; a=(level-u[i])/(u[j]-u[i]); answer.append(((i+.5+a)*C/n)%C)
    return np.array(answer)

def differences(e,reference,C):return (np.asarray(e)-reference+C/2)%C-C/2

def prepare():
    if RUN.exists():raise ValueError('Fresh run directory required')
    m=source.verify(ORIGINAL)
    for name in ('inputs','workers','logs','frozen_sources','results'): (RUN/name).mkdir(parents=True)
    for base in ('frozen_sources',):shutil.copytree(ORIGINAL/base,RUN/base,dirs_exist_ok=True)
    jobs=[]
    for b in m['basals']:
        shutil.copy2(ORIGINAL/b['input'],RUN/b['input'])
        for f in FRACTIONS:
            jobs.append(dict(slot=len(jobs),model=b['model'],size=b['size'],cells=b['cells'],platform=('spring' if b['model']=='spring_ring' else 'chemical'),fraction=f,mean=REF[b['model']]*f,input=b['input'],input_sha256=b['input_sha256'],basal_circumference=b['basal_circumference']))
    shutil.copy2(Path(__file__),RUN/'frozen_sources'/Path(__file__).name)
    # Dispatch helper expects this source in its original nested frozen layout.
    for name in ('figure4_round3_campaign_20261022.py','signature_figure4_campaign_20261020.py'):
        shutil.copy2(CODE/name,RUN/'frozen_sources'/name)
    manifest=dict(generated_utc=now(),jobs=jobs,native_protocols=m['native_protocols'],source_campaign=str(ORIGINAL),source_manifest_sha256=sha(ORIGINAL/'manifest.json'),
        amplitude_fraction_of_reference=.0001,amplitude_multipliers=[.5,1.,2.],pulse_shape='periodic nearest-arc Gaussian exp(-d^2/(2 ell^2)); peak at interpolated half-contrast edge',pulse_duration='tau',observation_horizon='20 tau',mean_preparation='uniform declared mean from original basal; 16000 chemical /1600 Spring; no noise',
        numerical_solver=dict(rtol=2e-10,atol=2e-12),resources=dict(qos='small',cpus=1,memory='4G',walltime='00:45:00'),
        input_sha256={j['input']:j['input_sha256'] for j in jobs},source_sha256={str(p.relative_to(RUN)):sha(p) for p in (RUN/'frozen_sources').rglob('*.py')})
    write(RUN/'manifest.json',manifest);write(RUN/'freeze.json',dict(manifest_sha256=sha(RUN/'manifest.json')))
    worker=RUN/'worker.sh';worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+source.PYTHON+' -B '+str(RUN/'frozen_sources'/Path(__file__).name)+' worker --slot "$SLURM_ARRAY_TASK_ID"\n');worker.chmod(0o750)
    return dict(run=str(RUN),jobs=len(jobs))

def covariance_integral(n,C,ell,tau,knot,sigma=1.,phases=24,steps_per_knot=12):
    """Integral_0^infty of phase-averaged native input covariance / mean^2.
    Exact lognormal image; discrete circular Gaussian covariance; analytic
    interpolation covariance. Integration tail at 20 tau and deterministic
    quadrature (no sampled forcing). Returns one periodic covariance row.
    """
    offsets=np.arange(n)*C/n; images=int(np.ceil(8*ell/C))+1
    k=sum(np.exp(-((offsets+j*C)/ell)**2) for j in range(-images,images+1));k/=np.linalg.norm(k)
    space=np.fft.irfft(abs(np.fft.rfft(k))**2,n=n)
    dt=knot/steps_per_knot; lags=np.arange(int(round(20*tau/dt))+1)*dt
    rho=np.exp(-knot/tau); values=np.zeros((len(lags),n))
    for phase in (np.arange(phases)+.5)/phases:
        first=phase; second=phase+lags/knot; index=np.floor(second).astype(int);b=second-index
        scale1=np.sqrt((1-first)**2+first**2+2*first*(1-first)*rho)
        scale2=np.sqrt((1-b)**2+b*b+2*b*(1-b)*rho)
        corr=((1-first)*(1-b)*rho**abs(index)+(1-first)*b*rho**abs(index+1)+first*(1-b)*rho**abs(index-1)+first*b*rho**abs(index))/(scale1*scale2)
        values+=np.expm1(sigma*sigma*corr[:,None]*space[None,:])/phases
    return np.trapezoid(values,lags,axis=0),space

def phase_kernel(rhs,jac,y,n,C,spring=False):
    """Native translation-phase susceptibility to each native input cell.
    Adjoint eigenvector of the translation mode, selected using the spatial
    derivative. Conserved-pool adjoints are orthogonal to this derivative.
    """
    species=(len(y)-1)//n if spring else len(y)//n
    fields=y[:species*n].reshape(species,n)
    # derivative of y(x-s) with respect to s; native periodic centered FD.
    q=(-(np.roll(fields,-1,axis=1)-np.roll(fields,1,axis=1))/(2*C/n)).ravel()
    if spring:q=np.r_[q,0.]
    J=jac(y)
    # Paired left/right eigenvectors. Small positive shift avoids exact mass
    # conservation singularity without regularizing the native dynamics.
    ev,vr=eigs(J,k=min(8,len(y)-2),sigma=1e-7,which='LM',tol=1e-10)
    order=np.argsort(-np.abs(vr.T.conj()@q)/np.linalg.norm(vr,axis=0))
    idx=int(order[0]); lam=ev[idx]
    el,vl=eigs(J.T,k=min(8,len(y)-2),sigma=1e-7,which='LM',tol=1e-10)
    near=np.argsort(abs(el-lam))[:min(4,len(el))]
    # Translation has large overlap with q; conservation left modes do not.
    ix=max(near,key=lambda i:abs(np.vdot(vl[:,i],q)))
    v=np.real(vl[:,ix]);v/=float(v@q)
    eps=1e-5
    zero=np.full(n,0.)
    # Native input is pointwise and acts as active/inactive exchange. Direct
    # finite difference of a uniform input recovers every local coefficient.
    B=(rhs(y,zero+eps)-rhs(y,zero))/eps
    br=B.reshape(-1)[:species*n].reshape(species,n)
    g=np.sum(v[:species*n].reshape(species,n)*br,axis=0)
    adjres=np.linalg.norm(J.T@v-float(lam.real)*v)/np.linalg.norm(v)
    return g,dict(translation_eigenvalue=float(lam.real),translation_eigenvalue_imag=float(lam.imag),adjoint_residual=float(adjres),right_translation_overlap=float(abs(np.vdot(vr[:,idx],q))/np.linalg.norm(vr[:,idx])/np.linalg.norm(q)),phase_normalization=float(v@q))

def worker(slot):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('SLURM required')
    m=read(RUN/'manifest.json');assert sha(RUN/'manifest.json')==read(RUN/'freeze.json')['manifest_sha256']
    c=m['jobs'][slot];folder=RUN/'workers'/f'{slot:03d}';folder.mkdir(exist_ok=True)
    if (folder/'status.json').exists() and read(folder/'status.json').get('complete'):return
    start=time.monotonic(); status=dict(generated_utc=now(),condition=c,complete=False,slurm_job_id=os.environ['SLURM_JOB_ID'])
    write(folder/'status.json',status)
    try:
        native,helper=source.native_module(RUN,c['platform'])
        with np.load(RUN/c['input']) as z:y0=z['state'].copy().ravel()
        assert sha(RUN/c['input'])==c['input_sha256']
        n=c['cells'];spring=c['platform']=='spring';spec=m['native_protocols'][c['platform']];ell=spec['noise']['ell'];tau=spec['noise']['tau'];duration=tau;horizon=20*tau
        if spring:
            closure=native.Closure(c['size'],n)
            rhs=lambda y,a:closure.rhs(y,a)
            jac=lambda y:sparse.csc_matrix(closure.jacobian(y,c['mean']))
            sparsity=None
        else:
            grid=native.ring.RingGrid(c['size'],n);model=native.ring.native_model(c['model']);shape=(len(model.state_names),n)
            rhs=lambda y,a:native.ring.native_rhs(model,grid,y.reshape(shape),np.broadcast_to(a,(n,))).ravel()
            sparsity=native.ring.periodic_sparsity(shape[0],n)
            jac=lambda y:approx_derivative(lambda v:rhs(v,c['mean']),y,method='3-point',sparsity=sparsity).tocsc()
        kwargs=dict(method='BDF',rtol=2e-10,atol=2e-12)
        if spring:kwargs['jac']=lambda t,y:jac(y)
        else:kwargs['jac_sparsity']=sparsity
        def solve(initial,a,end,times,maxstep):
            options=dict(kwargs)
            if spring:options['jac']=lambda t,y:sparse.csc_matrix(closure.jacobian(y,a))
            sol=solve_ivp(lambda t,y:rhs(y,a),(0.,end),initial,t_eval=times,max_step=maxstep,**options)
            if not sol.success:raise RuntimeError(sol.message)
            return sol.y.T
        warm_end=1600. if spring else 16000.
        warm=solve(y0,c['mean'],warm_end,np.array([0.,warm_end*.9,warm_end]),50. if spring else 500.)
        y=warm[-1];C=float(y[-1]) if spring else c['size'];initial_edges=edges(y[:n],C)
        residual=float(np.max(abs(rhs(y,c['mean']))));drift=float(np.max(abs(warm[-1]-warm[-2]))/max(np.max(abs(y)),1.))
        status.update(mean_stationary_rhs=residual,mean_tail_scaled_change=drift,mean_circumference=C,mean_contrast=float(np.ptp(y[:n])),mean_negative_min=float(y[:-1].min() if spring else y.min()))
        np.savez_compressed(folder/'mean_state.npz',state=y)
        if initial_edges is None or np.ptp(y[:n])<.05:
            status.update(complete=True,defined=False,reason='No unique polarized half-contrast boundary at the declared mean input');return
        if residual>1e-6 or drift>1e-4:
            status.update(complete=True,defined=False,reason='Mean-input state did not become stationary within the bounded preparation');return
        times=np.unique(np.r_[np.linspace(0,duration,51),np.linspace(duration,horizon,381)])
        control=solve(y,c['mean'],horizon,times,tau/2)
        Cs=control[:,-1] if spring else np.full(len(times),C)
        controls=np.array([edges(v[:n],cc) for v,cc in zip(control,Cs)])
        # Displacements are material fractions times the frozen mean-input C
        # for Spring; uniform expansion cannot masquerade as edge motion.
        controls=controls/Cs[:,None]*C
        curves=[];a0=REF[c['model']]*m['amplitude_fraction_of_reference']
        x=(np.arange(n)+.5)*C/n
        for edge_index in (0,1):
            d=(x-initial_edges[edge_index]+C/2)%C-C/2
            profile=np.exp(-.5*(d/ell)**2)
            for mult in m['amplitude_multipliers']:
                amp=a0*mult
                on_times=times[times<=duration];on=solve(y,c['mean']+amp*profile,duration,on_times,tau/20)
                off_times=times[times>duration]-duration
                off=solve(on[-1],c['mean'],horizon-duration,np.r_[0.,off_times],tau/2)[1:]
                states=np.vstack([on,off]);circ=states[:,-1] if spring else np.full(len(times),C)
                all_edges=[edges(v[:n],cc) for v,cc in zip(states,circ)]
                if any(e is None for e in all_edges):raise ValueError('Weak pulse lost unique edge topology')
                e=np.array(all_edges)/circ[:,None]*C
                delta=differences(e,controls,C)
                area=np.trapezoid(delta,times,axis=0)/(amp*duration)
                gain=delta[-1]/(amp*duration)
                curves.append(dict(pulse_edge=edge_index,multiplier=mult,amplitude=amp,delta=delta,integrated=area,gain=gain,relaxation_tail_change=np.max(abs(delta[-1]-delta[-11]))))
        g,linear=phase_kernel(rhs,jac,y,n,C,spring)
        ci,sp=covariance_integral(n,c['basal_circumference'] if spring else C,ell,tau,spec['noise']['knot'])
        # g maps input cell values to phase velocity. Circular convolution
        # includes the exact cell discretization; no dx factor is missing.
        I=float(np.dot(g,np.fft.irfft(np.fft.rfft(ci)*np.fft.rfft(g),n=n)))
        De=c['mean']**2*I
        linear.update(D_left=De,D_right=De,Ds_independent=De/2,Ds_full_crosscovariance=De,
            spatial_kernel='Adjoint of unchanged native linearization; both edges share the translation mode at zero frequency',covariance_time_integral_unit_mean=I,
            edge_gain_validation=[])
        for e in (0,1):
            d=(x-initial_edges[e]+C/2)%C-C/2;profile=np.exp(-.5*(d/ell)**2)
            predicted=float(g@profile)
            observed=next(v['gain'][e] for v in curves if v['pulse_edge']==e and v['multiplier']==1)
            linear['edge_gain_validation'].append(dict(edge=e,predicted_phase_gain=predicted,measured_edge_gain=float(observed),relative_error=float(abs(observed-predicted)/max(abs(predicted),1e-20))))
        linearity=[]
        for e in (0,1):
            base=next(v for v in curves if v['pulse_edge']==e and v['multiplier']==1)['delta']
            for mult in (.5,2.):
                test=next(v for v in curves if v['pulse_edge']==e and v['multiplier']==mult)['delta']/mult
                err=float(np.max(abs(test-base))/max(np.max(abs(base)),1e-16))
                linearity.append(dict(pulse_edge=e,multiplier=mult,relative_max_curve_error=err))
        np.savez_compressed(folder/'pulse_responses.npz',time=times,control_edges=controls,initial_edges=initial_edges,
            pulse_edges=np.array([r['pulse_edge'] for r in curves]),multipliers=np.array([r['multiplier'] for r in curves]),amplitudes=np.array([r['amplitude'] for r in curves]),
            displacement=np.array([r['delta'] for r in curves]),integrated=np.array([r['integrated'] for r in curves]),terminal_gain=np.array([r['gain'] for r in curves]),phase_sensitivity=g,input_covariance_integral=ci,gaussian_spatial_covariance=sp)
        status.update(complete=True,defined=True,linear_response=linear,amplitude=a0,linearity=linearity,maximum_linearity_error=max(r['relative_max_curve_error'] for r in linearity),pulse_relaxation_tail_change=max(float(r['relaxation_tail_change']) for r in curves),pulse_duration=duration,observation_horizon=horizon)
    except Exception as exc:
        import traceback
        status.update(error=str(exc),traceback=traceback.format_exc());raise
    finally:
        status.update(generated_utc=now(),wall_seconds=time.monotonic()-start);write(folder/'status.json',status)

def submit():
    cmd=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time=00:45:00','--job-name=f3edge5','--array=0-119','--output='+str(RUN/'logs/%A_%a.out'),str(RUN/'worker.sh')]
    p=subprocess.run(cmd,text=True,capture_output=True,check=True)
    write(RUN/'submission.json',dict(generated_utc=now(),command=cmd,stdout=p.stdout,stderr=p.stderr));return p.stdout

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','submit','worker'));parser.add_argument('--slot',type=int);a=parser.parse_args()
    print(prepare() if a.action=='prepare' else submit() if a.action=='submit' else worker(a.slot))
