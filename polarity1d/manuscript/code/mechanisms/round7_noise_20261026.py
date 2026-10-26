"""Round-7 strong-input native/fast-cytosol experiment and Figure 5 E/F.

No kinetic parameters are fitted. The source and exact restart fields are frozen
before SLURM execution. Native periodic finite differences use an IMEX Euler
step: explicit native reaction, implicit diffusion solved by a cyclic Thomas
algorithm. The forcing is the unchanged repository NoiseInput implementation.
"""
from pathlib import Path
import argparse, csv, hashlib, importlib.util, json, math, os, shutil, subprocess, sys, time
from datetime import datetime, timezone
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
from scipy.signal import fftconvolve
from numba import njit

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript'
BASE=MS/'data/runs/signature_figure4_round3_compact_20261022'
OUT=MS/'figures/length_mechanisms_20261008/dev/revision_round7_20261026'
RUN=OUT/'noise'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
MODELS=('wave_pinning','holmes_model3')
NSEED=64

def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def read(p):return json.loads(Path(p).read_text())
def table(p,rows):
    rows=list(rows);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as s:
        w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m

def prepare():
    if (RUN/'manifest.json').exists():raise ValueError('Already frozen')
    for part in ('supporting/data','supporting/metadata','supporting/vectors','supporting/viewers','results','logs','frozen_sources','inputs'):
        (RUN/part).mkdir(parents=True,exist_ok=True)
    frozen=RUN/'frozen_sources'/Path(__file__).name
    shutil.copy2(__file__,frozen)
    # Native package remains immutable and is addressed by exact input hashes.
    sources={str(p.relative_to(ROOT)):sha(p) for p in sorted((BASE/'frozen_sources').rglob('*.py'))}
    inputs={}
    for m in MODELS:
        for C in (20,40,60):
            p=BASE/f'inputs/{m}_C{C}.npz';d=RUN/'inputs'/p.name;shutil.copy2(p,d);inputs[str(d.relative_to(RUN))]=sha(d)
    jobs=[]
    for m in MODELS:
        for f in (.1,.3):
            for C in (20,40,60):
                for D in ((10.,1000.) if m=='wave_pinning' else (50.,1000.)):
                    for seed in range(NSEED):
                        jobs.append(dict(slot=len(jobs),model=m,f=f,C=C,D=D,replicate=seed,seed=202610260000+seed,
                            dt=.05,end=32000. if m=='wave_pinning' else 20000.,burn=16000. if m=='wave_pinning' else 4000.))
    # Paired temporal controls do not enter the production seed bootstrap.
    for m in MODELS:
        for D in ((10.,1000.) if m=='wave_pinning' else (50.,1000.)):
            for seed in range(4):
                jobs.append(dict(slot=len(jobs),model=m,f=.1,C=40,D=D,replicate=seed,seed=202610260000+seed,
                    dt=.025,end=32000. if m=='wave_pinning' else 20000.,burn=16000. if m=='wave_pinning' else 4000.,control=True))
    manifest=dict(created_utc=now(),jobs=jobs,seeds_per_condition=NSEED,source_sha256=sources,input_sha256=inputs,
        worker_source_sha256=sha(frozen),solver='IMEX Euler, periodic finite differences dx=.05, cyclic Thomas implicit diffusion',
        forcing='Unchanged native NoiseInput, Gaussian spatial covariance SD 1.5, OU tau50, knots5 coupled at2.5, lognormal sigma1',
        design='Same seed and noise realization paired across native/fast D; independently seeded replicates; all 3 sizes equally included',
        patch_thresholds=[.55,.65,.75],patch_minimum_grid_cells=10,frame_interval=5.)
    write(RUN/'manifest.json',manifest)
    script=RUN/'worker.sh';script.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+PY+' -B '+str(frozen)+' worker --slot "$SLURM_ARRAY_TASK_ID"\n');script.chmod(0o750)
    print(json.dumps(dict(jobs=len(jobs),manifest=str(RUN/'manifest.json'))))

@njit
def factors(n,r):
    b=1+2*r; gamma=-b
    d=np.empty(n);cp=np.empty(n)
    d[0]=2*b;cp[0]=-r/d[0]
    for i in range(1,n):
        bi=b+r*r/b if i==n-1 else b
        d[i]=bi+r*cp[i-1];cp[i]=-r/d[i]
    z=np.zeros(n);z[0]=gamma/d[0]
    for i in range(1,n):z[i]=((-r if i==n-1 else 0.)+r*z[i-1])/d[i]
    for i in range(n-2,-1,-1):z[i]-=cp[i]*z[i+1]
    return d,cp,z

@njit
def integrate(state,knots,mean,model,diff,dt,end,burn,dx):
    nf,n=state.shape;steps=int(round(end/dt));stride=int(round(5/dt));start=int(round(burn/dt))
    frames=(steps-start)//stride+1;profiles=np.empty((frames,n),np.float32)
    ds=np.empty((nf,n));cps=np.empty((nf,n));zs=np.empty((nf,n))
    for s in range(nf):ds[s],cps[s],zs[s]=factors(n,diff[s]*dt/dx**2)
    buf=np.empty((nf,n));out=np.empty((nf,n));stim=np.empty(n)
    sum0=np.empty(2);sum0[0]=np.sum(state[0]+state[1]);sum0[1]=np.sum(state[2]+state[3]) if nf==4 else 0.
    mass_error=0.;minimum=1e100;k=0;rho=np.exp(-5/50)
    for step in range(steps+1):
        if step>=start and (step-start)%stride==0:
            profiles[k,:]=state[0,:];k+=1
        if step==steps:break
        t=step*dt;j=min(int(t/5),len(knots)-2);a=t/5-j
        scale=np.sqrt((1-a)**2+a*a+2*a*(1-a)*rho)
        for i in range(n):stim[i]=mean*np.exp(((1-a)*knots[j,i]+a*knots[j+1,i])/scale-.5)
        for i in range(n):
            u=state[0,i];v=state[1,i]
            if model==0:
                transfer=v*(.067+u*u/(1+u*u)+stim[i])-u
            else:
                h=state[2,i];hi=state[3,i]
                transfer=(3.5/(2*(1+(h/1.25)**3))+stim[i])*v/7.5-u
                transfer2=(6.6/(2*(1+u**3)))*hi/3.1-h
                buf[2,i]=h+dt*transfer2;buf[3,i]=hi-dt*transfer2
            buf[0,i]=u+dt*transfer;buf[1,i]=v-dt*transfer
        for s in range(nf):
            r=diff[s]*dt/dx**2
            out[s,0]=buf[s,0]/ds[s,0]
            for i in range(1,n):out[s,i]=(buf[s,i]+r*out[s,i-1])/ds[s,i]
            for i in range(n-2,-1,-1):out[s,i]-=cps[s,i]*out[s,i+1]
            gamma=-(1+2*r)
            fact=(out[s,0]-r*out[s,n-1]/gamma)/(1+zs[s,0]-r*zs[s,n-1]/gamma)
            wanted=0.;found=0.
            for i in range(n):
                state[s,i]=out[s,i]-fact*zs[s,i]
                wanted+=buf[s,i];found+=state[s,i]
            # Diffusion preserves each species' zero mode exactly; restore
            # roundoff lost in the ill-conditioned fast-diffusion solve.
            correction=(wanted-found)/n
            for i in range(n):
                state[s,i]+=correction
                minimum=min(minimum,state[s,i])
        if step%1000==0:
            err=abs(np.sum(state[0]+state[1])-sum0[0])/sum0[0]
            if nf==4:err=max(err,abs(np.sum(state[2]+state[3])-sum0[1])/sum0[1])
            mass_error=max(mass_error,err)
    return profiles,mass_error,minimum

def components(mask):
    n=len(mask)
    if mask.all():return [np.arange(n)]
    starts=np.flatnonzero(mask & ~np.roll(mask,1));out=[]
    for start in starts:
        stop=start
        while mask[stop%n]:stop+=1
        if stop-start>=10:out.append(np.arange(start,stop)%n)
    return out

def patch_stats(p,C,threshold):
    n=p.shape[1];dx=C/n;angle=2*np.pi*(np.arange(n)+.5)/n
    count=[];width=[];vel=[];previous=[]
    for field in p:
        regions=components(field>threshold);count.append(len(regions));width.extend(len(c)*dx for c in regions)
        headings=[np.angle(np.exp(1j*angle[c]).sum()) for c in regions]
        # One-to-one mutual-best overlap excludes ambiguous splits and mergers.
        if regions and previous:
            overlap=np.array([[len(np.intersect1d(a[0],b,assume_unique=True)) for b in regions] for a in previous])
            for i,(oldcells,oldheading) in enumerate(previous):
                j=int(overlap[i].argmax())
                if overlap[i,j]>0 and int(overlap[:,j].argmax())==i:
                    delta=np.angle(np.exp(1j*(headings[j]-oldheading)))
                    vel.append(C*delta/(2*np.pi*5))
        previous=list(zip(regions,headings))
    return dict(threshold=threshold,patch_density=float(np.mean(count)/C),mean_patch_count=float(np.mean(count)),
        width=float(np.mean(width)) if width else None,speed_rms=float(np.sqrt(np.mean(np.square(vel)))) if vel else None,
        multipatch_fraction=float(np.mean(np.array(count)>1)),matched_transitions=len(vel))

def acf(x,maxlag=800):
    # Uncentered stationary spatial modes; normalization uses all lag-zero frames.
    x=np.asarray(x);n=len(x)
    raw=fftconvolve(x,np.conj(x[::-1]),mode='full')[n-1:n+maxlag]
    return (raw/np.arange(n,n-len(raw),-1)).real/float(np.mean(abs(x)**2))

def rate(c):
    c=np.asarray(c);valid=np.flatnonzero((c[1:]<=.5)&np.isfinite(c[1:]))+1
    if not len(valid):return None
    first=int(valid[0]);indices=[first-1,first] if first>1 else [first]
    j=min(indices,key=lambda i:abs(c[i]-.5))
    return float(-np.log(max(c[j],1e-12))/(5*j))

def worker(slot,benchmark=False):
    manifest=read(RUN/'manifest.json');job=manifest['jobs'][slot];target=RUN/('benchmark' if benchmark else 'results')/f'{slot:04d}'
    if (target/'status.json').exists() and read(target/'status.json').get('complete'):return
    target.mkdir(parents=True,exist_ok=True)
    frozen=RUN/'frozen_sources'/Path(__file__).name
    if sha(frozen)!=manifest['worker_source_sha256']:raise ValueError('Frozen worker changed')
    p=RUN/f"inputs/{job['model']}_C{job['C']}.npz"
    if sha(p)!=manifest['input_sha256'][str(p.relative_to(RUN))]:raise ValueError('Restart changed')
    helper=load_module('round7_native_helper',BASE/'frozen_sources/signature_figure4_campaign_20261020.py')
    native=helper.native_module(BASE,'chemical')
    with np.load(p) as z:state=z['state'].copy()
    model=0 if job['model']=='wave_pinning' else 1
    diff=np.array([.1,job['D']] if model==0 else [.1,job['D'],.1,job['D']])
    end=100. if benchmark else job['end'];burn=0. if benchmark else job['burn']
    reference=.1 if model==0 else 1.
    generator=native.NoiseInput(native.ring.RingGrid(job['C'],state.shape[1]),end,job['seed'],job['f']*reference)
    started=time.time()
    profiles,merr,minimum=integrate(state,np.asarray(generator.knots),job['f']*reference,model,diff,job['dt'],end,burn,.05)
    elapsed=time.time()-started
    if not np.isfinite(profiles).all() or minimum < -1e-8:raise ValueError('Numerically invalid trajectory')
    if merr>1e-7:raise ValueError('Mass conservation failed')
    coefficients=np.fft.rfft(profiles,axis=1)[:,:17]/profiles.shape[1]
    heading=coefficients[:,1]/abs(coefficients[:,1])
    coherence=abs(coefficients[:,1])/coefficients[:,0].real
    curves=np.stack([acf(coefficients[:,n],min(800,len(profiles)//4)) for n in range(1,17)])
    hcurve=acf(heading,min(800,len(profiles)//4))
    stats=[patch_stats(profiles,job['C'],t) for t in manifest['patch_thresholds']] if model==0 else []
    np.savez_compressed(target/'compact.npz',mode_acf=curves.astype('float32'),heading_acf=hcurve.astype('float32'),
        modes=coefficients.astype('complex64'),P_abs=coherence.astype('float32'))
    if job['replicate']<2 or job.get('control'):
        # All post-burn-in frames, physical dx=.05; float16 is illustration only.
        np.savez_compressed(target/'kymograph.npz',time=burn+np.arange(len(profiles))*5,field=profiles.astype('float16'),C=job['C'])
    write(target/'status.json',dict(job=job,complete=True,generated_utc=now(),elapsed_seconds=elapsed,
        mass_relative_error=merr,minimum_concentration=minimum,patch_statistics=stats,
        mode1_rate=rate(curves[0]),heading_rate=rate(hcurve),P_abs_quantiles=np.quantile(coherence,[0,.05,.25,.5,.75,.95,1]).tolist(),
        P_abs_mean=float(coherence.mean()),noise_knots_sha256=hashlib.sha256(np.asarray(generator.knots,dtype='<f8').tobytes()).hexdigest(),
        noise_reference='Frozen repository NoiseInput, unchanged',benchmark=benchmark))
    print(json.dumps(dict(slot=slot,elapsed=elapsed,mass_error=merr,minimum=minimum)))

def submit():
    manifest=read(RUN/'manifest.json')
    command=['sbatch','--parsable','--job-name=r7_noise','--qos=xxlarge','--cpus-per-task=1','--mem=2G','--time=02:00:00',
        '--array=0-'+str(len(manifest['jobs'])-1),'--output='+str(RUN/'logs/%A_%a.out'),'--error='+str(RUN/'logs/%A_%a.err'),str(RUN/'worker.sh')]
    result=subprocess.run(command,text=True,capture_output=True,check=True)
    write(RUN/'submission.json',dict(created_utc=now(),command=command,job_id=result.stdout.strip()));print(result.stdout)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','worker','submit']);ap.add_argument('--slot',type=int);ap.add_argument('--benchmark',action='store_true');a=ap.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='submit':submit()
    else:worker(a.slot,a.benchmark)
