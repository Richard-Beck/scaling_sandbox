from pathlib import Path
import sys,os,csv,json,hashlib,time,platform,argparse
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor
from model import simulate,metrics,msd,transition,rayleigh_p
ROOT=Path(__file__).resolve().parent;sys.path.append(str(ROOT.parent))
from data import DATA,BYID,Q,METRICS
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
# Shared controls; experimental summaries remain separate constraints.
GROUPS={'O':['O'],'N':['N','HN0','RN0'],'P':['P','HP0','RP0'],
        'HN':['HN1'],'HP':['HP1'],'RN':['RN1'],'RP':['RP1']}
OBS_TO_GROUP={id:g for g,ids in GROUPS.items() for id in ids}
PAIRS=[('N','P','Baseline PACC/non-PACC'),('N','HN','HIF · non-PACC'),('P','HP','HIF · PACC'),('N','RN','RhoA · non-PACC'),('P','RP','RhoA · PACC')]
TAU=np.geomspace(.02,32,49);BETA=np.linspace(0,2,41)
BANKMET=('path','directness','parallel','endpoint','perpendicular','angle_parallel')

def write_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def bank_task(args):
    i,n=args;p=simulate(1,TAU[i],0,n=n,seed=280926);times=np.arange(33)*.5
    out=np.empty((len(BETA),len(BANKMET),3))
    for j,beta in enumerate(BETA):
        pp=p.copy();pp[:,:,1]+=beta*times[None,:];m=metrics(pp)
        out[j]=[np.quantile(m[k],Q) for k in BANKMET]
    return i,out

def bank(n,workers):
    dest=OUT/'bank.npz'
    if dest.exists():
        b=np.load(dest)
        if int(b['cells'])!=n:raise ValueError('Bank cell-count mismatch')
        return b['q']
    qq=np.empty((len(TAU),len(BETA),len(BANKMET),3))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,res in pool.map(bank_task,[(i,n) for i in range(len(TAU))]):
            qq[i]=res
            if i%8==0:print(f'OU bank {i+1}/{len(TAU)} persistence slices',flush=True)
    np.savez_compressed(dest,q=qq,tau=TAU,beta=BETA,cells=n)
    return qq

class Loss:
    def __init__(self,q):self.interp=RegularGridInterpolator((np.log(TAU),BETA),q,bounds_error=False,fill_value=np.nan)
    def quantities(self,lt,beta):return self.interp([lt,beta])[0]
    def score(self,g,lt,beta,sigma=None,mode='full',weights=(1,1,1)):
        if not np.log(TAU[0])<=lt<=np.log(TAU[-1]) or not 0<=beta<=BETA[-1]:return 1e6,1.
        qq=self.quantities(lt,beta);w=np.array(weights);rows=[BYID[id] for id in GROUPS[g]]
        if sigma is None:
            numer=denom=0.
            for r in rows:
                target=np.array(r['path']);scale=max(target[2]-target[0],1.)
                numer+=np.sum(w*qq[0]*target)/scale**2;denom+=np.sum(w*qq[0]**2)/scale**2
            sigma=float(numer/denom)
        score=0.
        for r in rows:
            for metric in ('path','directness','parallel'):
                if metric=='parallel' and mode=='omit_direction':continue
                idx=BANKMET.index('angle_parallel' if metric=='parallel' and mode=='angle_normalization' else metric)
                predicted=qq[idx]*(sigma if metric=='path' else 1)
                target=np.array(r[metric]);scale=max(target[2]-target[0],1 if metric=='path' else .1)
                score+=np.sum(w*((predicted-target)/scale)**2)/w.sum()
        return float(score/len(rows)),float(sigma)
    def fit(self,g,mode='full',weights=(1,1,1)):
        bounds=[(float(np.log(TAU[0])),float(np.log(TAU[-1]))),(0.,2.)]
        if g=='O':
            rr=differential_evolution(lambda z:self.score(g,z[0],0,mode=mode,weights=weights)[0],[bounds[0]],seed=21,tol=1e-9,popsize=15,maxiter=180,polish=True)
            lt,beta=rr.x[0],0.
        else:
            rr=differential_evolution(lambda z:self.score(g,*z,mode=mode,weights=weights)[0],bounds,seed=21,tol=1e-8,popsize=15,maxiter=200,polish=True)
            lt,beta=rr.x
        sc,sigma=self.score(g,lt,beta,mode=mode,weights=weights)
        return dict(group=g,tau=float(np.exp(lt)),sigma=sigma,random_rms=float(np.sqrt(2)*sigma),u=float(beta*sigma),beta=float(beta),score=sc,
                    source_ids=';'.join(GROUPS[g]),tau_at_bound=bool(np.isclose(np.exp(lt),TAU[0],rtol=.01) or np.isclose(np.exp(lt),TAU[-1],rtol=.01)))

def reduced(loss,ga,gb,model,fits):
    ls=(np.log(.2),np.log(150));lt=(np.log(TAU[0]),np.log(TAU[-1]));ub=(0,40)
    if model=='amplitude_only':
        bounds=[ls,ls,lt,ub]
        def decode(z):return [(np.exp(z[0]),z[2],z[3]),(np.exp(z[1]),z[2],z[3])]
    elif model=='amplitude_and_bias':
        bounds=[ls,ls,lt,ub,ub]
        def decode(z):return [(np.exp(z[0]),z[2],z[3]),(np.exp(z[1]),z[2],z[4])]
    elif model=='amplitude_and_memory':
        bounds=[ls,ls,lt,lt,ub]
        def decode(z):return [(np.exp(z[0]),z[2],z[4]),(np.exp(z[1]),z[3],z[4])]
    else:raise ValueError(model)
    counts=[len(GROUPS[ga]),len(GROUPS[gb])]
    def f(z):
        total=0.
        for g,n,(s,t,u) in zip((ga,gb),counts,decode(z)):total+=n*loss.score(g,t,u/s,sigma=s)[0]
        return total/sum(counts)/3
    rr=differential_evolution(f,bounds,seed=42,tol=1e-7,popsize=20,maxiter=300,polish=True)
    pp=[dict(group=g,sigma=float(s),tau=float(np.exp(t)),u=float(u)) for g,(s,t,u) in zip((ga,gb),decode(rr.x))]
    return float(rr.fun),pp

def tests():
    checks=[]
    p,v=simulate(4,.8,2,n=32768,duration=4,seed=97,return_velocity=True)
    d=p[:,-1];sq=(d*d).sum(axis=1);expected=float(msd(4,4,.8,2));se=sq.std(ddof=1)/np.sqrt(len(sq))
    checks.append(dict(test='Integrated OU MSD',actual=float(sq.mean()),expected=expected,passed=bool(abs(sq.mean()-expected)<5*se)))
    mean=d[:,1].mean();se=d[:,1].std(ddof=1)/np.sqrt(len(d))
    checks.append(dict(test='Directed mean displacement',actual=float(mean),expected=8.,passed=bool(abs(mean-8)<5*se)))
    cov=v[:,0,0]*v[:,1,0];expected=16*np.exp(-.5/.8);se=cov.std(ddof=1)/np.sqrt(len(cov))
    checks.append(dict(test='Velocity memory',actual=float(cov.mean()),expected=float(expected),passed=bool(abs(cov.mean()-expected)<5*se)))
    m=metrics(p);err=np.max(abs(m['directness']**2-m['parallel']**2-m['perpendicular']**2))
    checks.append(dict(test='Directionality identity',actual=float(err),expected=0,passed=bool(err<1e-12)))
    p=simulate(0,.8,2,n=16);m=metrics(p)
    checks.append(dict(test='Deterministic drift limit',actual=float(np.max(abs(m['path']-32))),expected=0,passed=bool(np.allclose(m['path'],32))))
    # Integrating at half the sampling interval then downsampling gives same distribution.
    p1=simulate(4,.8,2,n=32768,seed=101);p2=simulate(4,.8,2,n=32768,cadence=.25,seed=102)[:,::2]
    x=metrics(p1)['path'];y=metrics(p2)['path'];se=np.sqrt(x.var()/len(x)+y.var()/len(y))
    checks.append(dict(test='Exact integration cadence consistency',actual=float(x.mean()-y.mean()),expected=0,passed=bool(abs(x.mean()-y.mean())<5*se)))
    write_csv(OUT/'tests.csv',checks)
    if not all(c['passed'] for c in checks):raise AssertionError(checks)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cells',type=int,default=4096);ap.add_argument('--workers',type=int,default=4);args=ap.parse_args()
    start=time.time();tests();q=bank(args.cells,args.workers);loss=Loss(q)
    fits={};sensitivity=[]
    for g in GROUPS:
        fits[g]=loss.fit(g);print('Shared OU fit',fits[g],flush=True)
        for mode,w,label in [('full',(0,1,0),'medians_only'),('omit_direction',(1,1,1),'without_directional_index'),('angle_normalization',(1,1,1),'index_divided_by_net_distance')]:
            sensitivity.append(dict(sensitivity=label,**loss.fit(g,mode,w)))
    comp=[]
    for ga,gb,label in PAIRS:
        full=sum(len(GROUPS[g])*fits[g]['score'] for g in (ga,gb))/sum(len(GROUPS[g]) for g in (ga,gb))/3
        for model in ('amplitude_only','amplitude_and_bias','amplitude_and_memory'):
            score,pp=reduced(loss,ga,gb,model,fits)
            comp.append(dict(pair=label,model=model,score=score,excess=score-full,parameters=json.dumps(pp)))
        comp.append(dict(pair=label,model='all_three',score=full,excess=0,parameters=json.dumps([fits[ga],fits[gb]])))
    write_csv(OUT/'fits.csv',list(fits.values()));write_csv(OUT/'sensitivity.csv',sensitivity);write_csv(OUT/'model_comparison.csv',comp)
    predictions=[];ray=[];timepred=[];validation=[];tracks={};counter=[];samples={}
    for g,f in fits.items():
        p,v=simulate(f['sigma'],f['tau'],f['u'],n=16384,seed=300926,return_velocity=True)
        m=metrics(p);samples[g]=m;tracks[g]=p[:24];qq={k:np.quantile(m[k],Q) for k in BANKMET}
        for k in METRICS:
            for j,quant in enumerate(Q):predictions.append(dict(group=g,metric=k,quantile=quant,predicted=float(qq[k][j]),used_in_fit=k in ('path','directness','parallel')))
        qb=loss.quantities(np.log(f['tau']),f['beta']).copy();qb[[0,3]]*=f['sigma']
        for ik,k in enumerate(BANKMET):
            for j,quant in enumerate(Q):validation.append(dict(group=g,metric=k,quantile=quant,bank=float(qb[ik,j]),independent=float(qq[k][j])))
        # Directional test checks: retained numerically; no standalone p-value panel.
        rng=np.random.default_rng(131);ind=rng.integers(0,len(p),(3000,60));ps=rayleigh_p(m['angle'][ind])
        for id in GROUPS[g]:ray.append(dict(source=id,group=g,reported_p=BYID[id]['rayleigh'],p05=float(np.quantile(ps,.05)),p50=float(np.median(ps)),p95=float(np.quantile(ps,.95))))
        for t in np.linspace(.5,32,64):timepred.append(dict(group=g,time_h=float(t),msd=float(msd(t,f['sigma'],f['tau'],f['u'])),mean_forward_displacement=float(f['u']*t),velocity_correlation=float(np.exp(-t/f['tau']))))
        if g in ('O','N','P'):
            # Per-cell time-varying speeds for model specification illustration.
            tracks[g+'_velocity']=v[:6]
            p0=simulate(f['sigma'],f['tau'],0,n=16384,seed=531);m0=metrics(p0)
            counter.append(dict(group=g,drift_removed_directness=float(np.median(m0['directness'])),drift_removed_endpoint=float(np.median(m0['endpoint']))))
    # Replace redundant control predictions with a source -> shared model mapping.
    mapped=[]
    for r in DATA:
        for k in METRICS:
            for j,quant in enumerate(Q):mapped.append(dict(source=r['id'],group=OBS_TO_GROUP[r['id']],metric=k,quantile=quant,observed=float(r[k][j])))
    # Qualitative endpoint-class check, retained in CSV only, not an unjustified numerical validation panel.
    sub=[]
    for group,kd,block in [('P','HP','HIF'),('P','RP','RhoA')]:
        threshold=np.quantile(samples[group]['endpoint'],.75)
        for g in (group,kd):sub.append(dict(block=block,group=g,threshold=float(threshold),positive_high_fraction=float(np.mean((samples[g]['endpoint']>threshold)&(samples[g]['parallel']>0)))))
    for name,rows in [('predictions',predictions),('source_observations',mapped),('rayleigh_checks',ray),('time_predictions',timepred),('numerical_checks',validation),('drift_removed',counter),('qualitative_subpopulation_check',sub)]:write_csv(OUT/(name+'.csv'),rows)
    np.savez_compressed(OUT/'tracks.npz',**tracks)
    # Profile amplitude/bias while varying memory, to show identifiability without claiming confidence intervals.
    profiles=[]
    for g in ('O','N','P'):
        for lt in np.log(TAU):
            if g=='O':beta=0.
            else:
                rr=differential_evolution(lambda z:loss.score(g,lt,z[0])[0],[(0,2)],seed=12,popsize=8,maxiter=60,tol=1e-6,polish=True);beta=float(rr.x[0])
            sc,sig=loss.score(g,lt,beta)
            profiles.append(dict(group=g,tau=float(np.exp(lt)),sigma=sig,u=float(sig*beta),score=sc/3))
    write_csv(OUT/'memory_profiles.csv',profiles)
    source=[ROOT/'model.py',ROOT/'run.py',ROOT.parent/'data.py',ROOT.parents[2]/'literature_review/hosny_2026.md']
    meta=dict(arguments=vars(args),groups=GROUPS,OU='dV=-(V-u*e_y)dt/tau + sqrt(2*sigma^2/tau)dW; dr=Vdt',
              sampling='exact joint Gaussian position/velocity updates; 0.5h cadence; 16h observation; stationary velocities',
              objective='Mean squared quantile residuals scaled by each reported IQR, averaged over datasets and metrics; not a likelihood',
              fit_quantiles=list(Q),fit_metrics=['path','directness','parallel'],bias='jointly fitted mean velocity u, not angular steering kappa',
              cells_independent_validation=16384,wall_seconds=time.time()-start,python=platform.python_version(),numpy=np.__version__,
              source_sha256={str(p.relative_to(ROOT.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in source})
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('OU revision finished',flush=True)
if __name__=='__main__':main()
