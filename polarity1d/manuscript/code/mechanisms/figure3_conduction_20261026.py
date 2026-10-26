"""Round7 deterministic conduction perturbations and isolated-front prediction.
Heavy integration is admitted only in SLURM. Native kinetics remain unchanged.
"""
from pathlib import Path
import os,sys,json,hashlib,argparse,time
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
from scipy.integrate import solve_ivp,quad,cumulative_trapezoid
from scipy.optimize import brentq
from scipy.optimize._numdiff import approx_derivative
from scipy import sparse
from dataclasses import replace
import figure3_round5_edge_20261023 as edge
source=edge.source; ROOT=edge.ROOT; MS=edge.MS; CODE=edge.CODE
OUT=MS/'figures/length_mechanisms_20261008/dev/revision_round7_20261026/figure3'
DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
DV={'wave_pinning':[2.5,5,10,20,40,100,1000], 'holmes_model3':[5,12.5,25,50,200,5000], 'spring_ring':[.1,.3,1,3,10,100,1000]}
NATIVE={'wave_pinning':10.,'holmes_model3':50.,'spring_ring':10.}
def write(p,v):source.write(p,v)
def front():
    """Maxwell front, quadratures of its translational adjoint; no ring fit."""
    a=.067+.0001;du=.1;ell=1.5
    A=lambda u:a+u*u/(1+u*u)
    def roots(v):return sorted(x.real for x in np.roots([1,-v*(1+a),1,-v*a]) if abs(x.imag)<1e-8)
    def area(v):
        r=roots(v)
        return quad(lambda u:v*A(u)-u,r[0],r[-1],epsabs=1e-13)[0]
    v=brentq(area,1.78,1.83,xtol=1e-14);lo,mid,hi=roots(v)
    def integrate(n):
        u=np.linspace(lo+1e-7,hi-1e-7,n)
        F=v*((1+a)*(u-lo)-(np.arctan(u)-np.arctan(lo)))-(u*u-lo*lo)/2
        p=np.sqrt(np.maximum(-2/du*F,1e-30));x=cumulative_trapezoid(1/p,u,initial=0);x-=np.interp((lo+hi)/2,u,x)
        I=np.trapezoid(p,u);kap=quad(A,lo,hi,epsabs=1e-13)[0]/I;lam=v*np.trapezoid(np.exp(-x*x/(2*ell*ell)),u)/I
        return kap,lam,x,u
    k,l,x,u=integrate(40001);kc,lc,_,_=integrate(20001)
    record=dict(v_star=v,low=lo,unstable=mid,high=hi,delta_u=hi-lo,kappa=k,lambda_push=l,mean_input=.0001,Du=du,pulse_sd=ell,maxwell_residual=area(v),quadrature_relative_change=max(abs(k/kc-1),abs(l/lc-1)),method='Maxwell equal-area root; u_x=sqrt(-2 integral f/Du); kappa=integral f_v du/integral u_x du; lambda=v*integral Gaussian(x(u)) du/integral u_x du; no ring data fitted')
    META.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    write(META/'isolated_front.json',record)
    take=np.unique(np.linspace(0,len(u)-1,1001,dtype=int));pd.DataFrame(dict(x=x[take],u=u[take])).to_csv(DATA/'isolated_front.csv',index=False)
    print(json.dumps(record,indent=2));return record

def worker(model):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Native integrations require SLURM')
    DATA.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    manifest=source.read(edge.RUN/'manifest.json');jobs=[c for c in manifest['jobs'] if c['model']==model and c['fraction']==.001]
    native,helper=source.native_module(edge.RUN,'spring' if model=='spring_ring' else 'chemical')
    rows=[];snapshots={};started=time.monotonic();spring=model=='spring_ring'
    for c in jobs:
        n=c['cells'];size=c['size'];path=edge.RUN/'workers'/f"{c['slot']:03d}"/'mean_state.npz'
        with np.load(path) as z:ybase=z['state'].copy().ravel()
        mean=.0004 if spring else c['mean'];ell=.15 if spring else 1.5
        for D in sorted(DV[model],key=lambda z:abs(np.log(z/NATIVE[model]))):
            t0=time.monotonic()
            if spring:
                closure=native.Closure(size,n)
                def correction(y):
                    C=y[-1];v=y[n:-1];return (D-10)*(np.roll(v,1)-2*v+np.roll(v,-1))*(n/C)**2
                def rhs(y,a):
                    v=closure.rhs(y,a);v[n:-1]+=correction(y);return v
                def jac(y):
                    J=closure.jacobian(y,mean);C=y[-1];ii=np.arange(n);lap=sparse.diags([-2*np.ones(n),np.ones(n-1),np.ones(n-1)],[0,-1,1],shape=(n,n)).toarray();lap[0,-1]=lap[-1,0]=1;J[n:-1,n:-1]+=(D-10)*(n/C)**2*lap;J[n:-1,-1]+=-2*correction(y)/C;return sparse.csc_matrix(J)
                shape=None;sparsity=None
            else:
                grid=native.ring.RingGrid(size,n);original=native.ring.native_model(model);mod=replace(original,**({'dv':D} if model=='wave_pinning' else {'inactive_diffusion':D}));shape=(len(mod.state_names),n)
                def rhs(y,a):return native.ring.native_rhs(mod,grid,y.reshape(shape),np.broadcast_to(a,(n,))).ravel()
                sparsity=native.ring.periodic_sparsity(shape[0],n)
                def jac(y):return approx_derivative(lambda v:rhs(v,mean),y,method='3-point',sparsity=sparsity).tocsc()
            end=1600 if spring else 16000
            sol=solve_ivp(lambda t,y:rhs(y,mean),(0,end),ybase,method='BDF',t_eval=[end*.9,end],rtol=2e-10,atol=2e-12,jac=lambda t,y:jac(y),max_step=50 if spring else 500)
            if not sol.success:raise RuntimeError(sol.message)
            y=sol.y[:,-1];C=float(y[-1]) if spring else size;e=edge.edges(y[:n],C);res=float(np.max(abs(rhs(y,mean))));row=dict(model=model,size=size,C=C,D=D,mean=mean,cells=n,contrast=float(np.ptp(y[:n])),rhs_residual=res,negative_min=float(y.min()),input_state=str(path.relative_to(ROOT)),input_state_sha256=source.sha(path),slurm_job_id=os.environ['SLURM_JOB_ID'],defined=e is not None and np.ptp(y[:n])>.05 and res<1e-6)
            if row['defined']:
                g,quality=edge.phase_kernel(rhs,jac,y,n,C,spring)
                x=(np.arange(n)+.5)*C/n;d=(x-e[1]+C/2)%C-C/2;h=np.exp(-d*d/(2*ell*ell));gain=float(g@h);row.update(gain=abs(gain),theta=float(((e[1]-e[0])%C)/C),**quality)
                snapshots[f'D{D:g}_C{size:g}']=y
            else:row.update(gain=None,theta=None,reason='polarized branch absent at tested D' if e is None or np.ptp(y[:n])<=.05 else 'bounded steady preparation unresolved')
            row['elapsed_seconds']=time.monotonic()-t0;rows.append(row)
            pd.DataFrame(rows).to_csv(DATA/f'conduction_{model}.csv',index=False)
            write(META/f'{model}_progress.json',dict(complete=False,finished=len(rows),total=len(jobs)*len(DV[model]),elapsed_seconds=time.monotonic()-started))
            print(model,size,D,row.get('gain'),res,flush=True)
    np.savez_compressed(DATA/f'conduction_states_{model}.npz',**snapshots)
    write(META/f'{model}_progress.json',dict(complete=True,finished=len(rows),total=len(rows),elapsed_seconds=time.monotonic()-started))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['front','worker']);p.add_argument('--model',choices=list(DV));a=p.parse_args()
    front() if a.command=='front' else worker(a.model)
