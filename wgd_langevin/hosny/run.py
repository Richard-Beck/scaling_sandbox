"""Summary-quantile optimization and prediction; supplied Hosny notes only."""
from pathlib import Path
import os,csv,json,hashlib,platform,time,argparse
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import differential_evolution, minimize
from scipy.special import i0,i1
from concurrent.futures import ProcessPoolExecutor
from model import simulate,metrics,rayleigh_p
from data import DATA,BYID,PAIRS,METRICS,Q

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
TAUS=np.geomspace(.04,12,33)
AS=np.array([0,.025,.05,.08,.12,.18,.25,.35,.45,.6,.8,1.,1.3,1.7,2.2,3.])
BANKMET=('path','directness','parallel','endpoint','perpendicular','angle_parallel','angle_perpendicular')


def write_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def task(args):
    i,j,n=args;p=simulate(1,TAUS[i],AS[j],n=n,seed=190926)
    m=metrics(p)
    quant=np.array([np.quantile(m[k],Q) for k in BANKMET])
    return i,j,quant

def bank(workers,n):
    path=OUT/'simulation_bank.npz'
    if path.exists():
        b=np.load(path)
        if int(b['n'])!=n or not np.array_equal(b['tau'],TAUS) or not np.array_equal(b['a'],AS):
            raise RuntimeError('Existing bank settings mismatch; choose a new output directory')
        return b['quantiles']
    simulate(1,1,.3,n=2,duration=.5)
    quant=np.empty((len(TAUS),len(AS),len(BANKMET),3))
    tasks=[(i,j,n) for i in range(len(TAUS)) for j in range(len(AS))]
    start=time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for k,(i,j,q) in enumerate(pool.map(task,tasks,chunksize=1),1):
            quant[i,j]=q
            if k%32==0: print(f'Bank {k}/{len(tasks)}; {time.time()-start:.0f}s',flush=True)
    np.savez_compressed(path,quantiles=quant,tau=TAUS,a=AS,n=n)
    return quant

class Objective:
    def __init__(self,quant):
        self.quant=quant
        self.interp=RegularGridInterpolator((np.log(TAUS),AS),quant,bounds_error=False,fill_value=np.nan)
    def score(self,id,x,mode='full',weights=(1,1,1)):
        logtau,a=x
        if logtau<np.log(TAUS[0]) or logtau>np.log(TAUS[-1]) or a<0 or a>AS[-1]:return 1e6,1.
        q=self.interp([logtau,a])[0]
        r=BYID[id];w=np.array(weights,float)
        target=np.array(r['path']);u=q[0]
        v=float(np.sum(w*u*target)/np.sum(w*u*u))
        score=np.sum(w*((v*u-target)/max(target[2]-target[0],1))**2)/w.sum()
        for k in ('directness','parallel'):
            if mode=='omit_fmi' and k=='parallel':continue
            target=np.array(r[k]);idx=BANKMET.index('angle_parallel' if mode=='angle_fmi' and k=='parallel' else k)
            score+=np.sum(w*((q[idx]-target)/max(target[2]-target[0],.1))**2)/w.sum()
        return float(score),v
    def optimize(self,id,mode='full',weights=(1,1,1),fixed_tau=None,fixed_kappa=None):
        bounds=[(float(np.log(TAUS[0])),float(np.log(TAUS[-1]))),(0.,float(AS[-1]))]
        if id=='O':
            f=lambda z:self.score(id,[float(z[0]),0],mode,weights)[0]
            rr=differential_evolution(f,[bounds[0]],seed=12,tol=1e-8,maxiter=100,polish=True)
            x=np.array([rr.x[0],0.])
        else:
            f=lambda x:self.score(id,x,mode,weights)[0]
            rr=differential_evolution(f,bounds,seed=17,tol=1e-7,maxiter=150,popsize=12,polish=True)
            x=rr.x
        sc,v=self.score(id,x,mode,weights)
        return dict(id=id,tau=float(np.exp(x[0])),a=float(x[1]),kappa=float(x[1]/np.exp(x[0])),v=v,score=sc,
                    tau_at_bound=bool(np.isclose(np.exp(x[0]),TAUS[0],rtol=.02) or np.isclose(np.exp(x[0]),TAUS[-1],rtol=.02)),
                    a_at_upper_bound=bool(x[1]>AS[-1]-.02))


def joint(obj,ids,model):
    # Sharing steering means kappa, not the dimensionless a=kappa*tau.
    lt=(float(np.log(TAUS[0])),float(np.log(TAUS[-1])))
    if model=='speed_only':
        bounds=[lt,(0,3)]
        def xs(z):return [z,z]
    elif model=='speed_and_steering':
        bounds=[lt,(0,3),(0,3)]
        def xs(z):return [[z[0],z[1]],[z[0],z[2]]]
    elif model=='speed_and_persistence':
        bounds=[lt,lt,(0,3/TAUS[0])]
        def xs(z):return [[z[0],z[2]*np.exp(z[0])],[z[1],z[2]*np.exp(z[1])]]
    else:raise ValueError(model)
    def objective(z):return sum(obj.score(id,x)[0] for id,x in zip(ids,xs(z)))
    rr=differential_evolution(objective,bounds,seed=32,tol=1e-7,popsize=15,maxiter=200,polish=True)
    xx=xs(rr.x);params=[]
    for id,x in zip(ids,xx):
        sc,v=obj.score(id,x);params.append(dict(id=id,tau=float(np.exp(x[0])),kappa=float(x[1]/np.exp(x[0])),v=v,a=float(x[1])))
    return float(rr.fun),params


def validation_tests():
    checks=[]
    p=simulate(4,.8,0,n=8192,duration=4,seed=101)
    end=np.sum(p[:,-1]**2,axis=1)
    expected=2*4**2*.8**2*(4/.8-1+np.exp(-4/.8))
    se=end.std(ddof=1)/np.sqrt(len(end))
    checks.append(dict(check='Unbiased MSD',observed=float(end.mean()),expected=expected,se=float(se),passed=bool(abs(end.mean()-expected)<5*se)))
    p=simulate(4,.8,.7,n=8192,duration=4,seed=102)
    expected=4*4*i1(.7)/i0(.7);vals=p[:,-1,1];se=vals.std(ddof=1)/np.sqrt(len(vals))
    checks.append(dict(check='Stationary directed displacement',observed=float(vals.mean()),expected=expected,se=float(se),passed=bool(abs(vals.mean()-expected)<5*se)))
    m=metrics(p);err=np.max(np.abs(m['parallel']**2+m['perpendicular']**2-m['directness']**2))
    checks.append(dict(check='FMI/directness identity',observed=float(err),expected=0,se=0,passed=bool(err<1e-12)))
    passed=all(r['passed'] for r in checks)
    write_csv(OUT/'validation_tests.csv',checks)
    if not passed:raise AssertionError(checks)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--bank-cells',type=int,default=1536)
    args=ap.parse_args();start=time.time()
    validation_tests();q=bank(args.workers,args.bank_cells);obj=Objective(q)
    fits={};sens=[]
    for r in DATA:
        fit=obj.optimize(r['id']);fits[r['id']]=fit
        print('Fit',r['id'],fit,flush=True)
        for mode,weights,label in [('full',(0,1,0),'medians_only'),('omit_fmi',(1,1,1),'omit_FMI'),('angle_fmi',(1,1,1),'FMI_normalized_by_net_distance')]:
            sens.append(dict(sensitivity=label,**obj.optimize(r['id'],mode,weights)))
    comparisons=[]
    for a,b,label in PAIRS:
        for model in ('speed_only','speed_and_steering','speed_and_persistence'):
            sc,pp=joint(obj,(a,b),model)
            comparisons.append(dict(pair=label,model=model,score=sc,score_per_metric=sc/6,parameters=json.dumps(pp)))
        sc=fits[a]['score']+fits[b]['score']
        comparisons.append(dict(pair=label,model='all_three',score=sc,score_per_metric=sc/6,parameters=json.dumps([fits[a],fits[b]])))
    write_csv(OUT/'fits.csv',list(fits.values()));write_csv(OUT/'sensitivity.csv',sens);write_csv(OUT/'model_comparison.csv',comparisons)
    # Independent larger ensembles: all predictions come from these unrefitted fits.
    prediction=[];mc=[];ray=[];subclasses=[];trajectories={};timecurves=[]
    for r in DATA:
        id=r['id'];f=fits[id]
        p=simulate(f['v'],f['tau'],f['a'],n=8192,seed=270926,resolution=1)
        trajectories[id]=p;m=metrics(p)
        qsim=np.array([np.quantile(m[k],Q) for k in BANKMET])
        qb=obj.interp([np.log(f['tau']),f['a']])[0].copy();qb[0]*=f['v'];qb[3]*=f['v']
        for ik,k in enumerate(METRICS):
            ii=BANKMET.index(k)
            for j,level in enumerate(Q):
                prediction.append(dict(id=id,metric=k,quantile=level,observed=float(r[k][j]),predicted=float(qsim[ii,j]),fitted=k in ('path','directness','parallel')))
                mc.append(dict(id=id,metric=k,quantile=level,bank=float(qb[ii,j]),independent=float(qsim[ii,j])))
        # Independent 60-track batches, conditional on independent cells (not replicate-level CI).
        rng=np.random.default_rng(388);idx=rng.integers(0,len(p),(3000,60));pp=rayleigh_p(m['angle'][idx])
        ray.append(dict(id=id,reported_p=r['rayleigh'],predicted_p05=float(np.quantile(pp,.05)),predicted_p50=float(np.median(pp)),predicted_p95=float(np.quantile(pp,.95)),fraction_p_lt_05=float(np.mean(pp<.05))))
        # Full duration predictions and intervention-free comparison.
        for k in range(1,p.shape[1]):
            mm=metrics(p[:,:k+1]);d=p[:,k]-p[:,0]
            timecurves.append(dict(id=id,time_h=k*.5,msd=float(np.mean(np.sum(d*d,axis=1))),median_endpoint=float(np.median(mm['endpoint'])),mean_parallel=float(np.mean(mm['parallel'])),median_directness=float(np.median(mm['directness']))))
        # Resolution and initial-heading sensitivity for baseline; no refitting.
        if r['block']=='baseline':
            for kind,reso,initial in [('double_resolution',2,'stationary'),('isotropic_start',1,'isotropic')]:
                pp=simulate(f['v'],f['tau'],f['a'],n=8192,seed=270927,resolution=reso,initial=initial);mm=metrics(pp)
                for k in METRICS:
                    vv=np.quantile(mm[k],Q)
                    for j,lev in enumerate(Q):mc.append(dict(id=id,metric=k,quantile=lev,bank=float(qsim[BANKMET.index(k),j]),independent=float(vv[j]),check=kind))
    for block in ('HIF','RhoA'):
        control=next(r for r in DATA if r['block']==block and r['group']=='PACC' and r['treatment']=='control')['id']
        cm=metrics(trajectories[control]);thresholds=np.quantile(cm['endpoint'],[.25,.75])
        for r in DATA:
            if r['block']!=block:continue
            mm=metrics(trajectories[r['id']]);states=np.digitize(mm['endpoint'],thresholds)
            for state in range(3):
                for sign in (-1,1):
                    frac=np.mean((states==state)&((mm['parallel']>=0) if sign==1 else (mm['parallel']<0)))
                    subclasses.append(dict(id=r['id'],class_name=('low','intermediate','high')[state],direction='positive' if sign==1 else 'negative',fraction=float(frac),low_threshold=float(thresholds[0]),high_threshold=float(thresholds[1])))
    # Quantitative counterfactual: steering removed, propulsion/tau unchanged.
    counter=[]
    for id in ('N','P'):
        f=fits[id];p=simulate(f['v'],f['tau'],0,n=8192,seed=711)
        m=metrics(p)
        counter.append(dict(id=id,median_directness=float(np.median(m['directness'])),median_endpoint=float(np.median(m['endpoint'])),median_parallel=float(np.median(m['parallel']))))
    write_csv(OUT/'predictions.csv',prediction);write_csv(OUT/'numerical_checks.csv',mc);write_csv(OUT/'rayleigh_checks.csv',ray)
    write_csv(OUT/'time_predictions.csv',timecurves);write_csv(OUT/'subpopulation_predictions.csv',subclasses);write_csv(OUT/'steering_removed.csv',counter)
    np.savez_compressed(OUT/'illustrative_tracks.npz',**{id:p[:32] for id,p in trajectories.items()})
    meta=dict(arguments=vars(args),wall_seconds=time.time()-start,seed_bank=190926,seed_validation=270926,validation_cells=8192,
              python=platform.python_version(),numpy=np.__version__,fit_quantiles=list(Q),fit_metrics=['path','directness','parallel'],
              population='homogeneous',FMI_assumption='net displacement projection / accumulated distance',initial_headings='stationary von Mises',
              objective='Sum of per-metric mean squared quantile errors scaled by observed IQR; not a likelihood or confidence interval',
              inputs={str(p.relative_to(ROOT.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'run.py',ROOT/'model.py',ROOT/'data.py',ROOT.parent.parent/'literature_review/hosny_2026.md']})
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('Finished optimization and predictions',flush=True)
if __name__=='__main__':main()
