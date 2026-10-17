from pathlib import Path
import sys,json,time,hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import differential_evolution,minimize
ROOT=Path(__file__).resolve().parent;PARENT=ROOT.parent;EX=PARENT/'exhaustive'
sys.path.insert(0,str(EX));import fit_frontier as ff
from scale_stats import score_and_scale
from run import write_csv
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
PARAMS=('tau','rate','amplitude','kappa')
RATIOS=np.unique(np.r_[np.geomspace(.03,30,25),1.])

def init():
    global REPS,REF,BOUNDS,MAPPING,AMPIDS,ZREF
    ff.init_worker();allfits=json.load((EX/'results/fits.json').open())
    REPS={b:next(f for f in allfits if f['base']==b and f['hmask']==8 and f['rmask']==8) for b in ('N','P')}
    REF=np.mean([f['bank_score'] for f in REPS.values()]);BOUNDS,MAPPING,AMPIDS=ff.schema(8,8)
    ZREF=np.r_[REPS['N']['z'],REPS['P']['z']]

def stats(z,base):return score_and_scale(z,MAPPING,AMPIDS,ff.BANK,ff.XGRID,ff.YGRID,ff.KGRID,ff.TARGETS[base],ff.SCALES[base],ff.OWNERS)

def joint(z,parameter,ratio,return_parameters=False):
    sn,nn,dn=stats(z[:5],'N');sp,np_,dp=stats(z[5:],'P')
    an=nn/dn;ap=np_/dp;score=(sn+sp)/2
    if parameter=='amplitude':
        a=(nn+ratio*np_)/(dn+ratio*ratio*dp)
        score+=(dn*(a-an)**2+dp*(ratio*a-ap)**2)/90
        an=a;ap=a*ratio
    if not return_parameters:return score
    pars=[]
    for base,zz,amplitude in [('N',z[:5],an),('P',z[5:],ap)]:
        _,fit=ff.LOSS.evaluate([base,'H'+base,'R'+base],zz[MAPPING],AMPIDS,return_fits=True)
        for f in fit:f['amplitude']=float(amplitude)
        pars += [{key:f[key] for key in ('group','tau','rate','amplitude','kappa')} for f in fit]
    return score,pars

def task(args):
    parameter,ratio=args;bounds=list(BOUNDS)*2;fixed_index={'tau':0,'rate':1,'kappa':2}.get(parameter)
    if fixed_index is None:
        def decode(w):return w
        start=ZREF.copy()
    else:
        j=fixed_index;lo,hi=BOUNDS[j]
        if j in (0,1):bounds[j]=(max(lo,lo-np.log(ratio)),min(hi,hi-np.log(ratio)))
        else:bounds[j]=(0.,float(np.log1p(8/max(1,ratio))))
        del bounds[5+j]
        def decode(w):
            z=np.insert(w,5+j,0.)
            z[5+j]=z[j]+np.log(ratio) if j in (0,1) else np.log1p(np.expm1(z[j])*ratio)
            return z
        start=np.delete(ZREF,5+j)
    start=np.array([np.clip(x,a,b) for x,(a,b) in zip(start,bounds)])
    def obj(w):return joint(decode(w),parameter,ratio)
    candidates=[(obj(start),start,False)]
    for seed in (31,53):
        rr=differential_evolution(obj,bounds,seed=seed,popsize=12,maxiter=320,tol=1e-7,polish=True,x0=start)
        candidates.append((float(rr.fun),rr.x,bool(rr.success)))
    rr=minimize(obj,start,method='L-BFGS-B',bounds=bounds,options={'maxiter':300,'ftol':1e-12});candidates.append((float(rr.fun),rr.x,bool(rr.success)))
    score,w,success=min(candidates,key=lambda v:v[0]);score,parameters=joint(decode(w),parameter,ratio,True)
    actual=parameters[3][parameter]/parameters[0][parameter]
    assert abs(np.log(actual/ratio))<1e-8,(parameter,ratio,actual)
    return dict(parameter=parameter,ratio=float(ratio),score=float(score),excess=float(score-REF),converged=success,parameters=parameters)

def main():
    start=time.time();init();fits=json.load((EX/'results/fits.json').open());ns={(f['hmask'],f['rmask']):f for f in fits if f['base']=='N'};ps={(f['hmask'],f['rmask']):f for f in fits if f['base']=='P'}
    rows=[]
    for key,n in ns.items():
        p=ps[key];rows.append(dict(hmask=key[0],rmask=key[1],HIF=n['HIF'],RhoA=n['RhoA'],complexity_per_population=n['complexity'],total_KD_parameters=2*n['complexity'],total_parameters=8+2*n['complexity'],score=(n['independent_score']+p['independent_score'])/2,bank_score=(n['bank_score']+p['bank_score'])/2,nonPACC_score=n['independent_score'],PACC_score=p['independent_score']))
    current=np.inf
    for k in (0,2,3,4,5,6,7,8):
        candidates=[r for r in rows if r['complexity_per_population']==k];best=min(candidates,key=lambda r:r['score'])
        for r in candidates:r['best_at_complexity']=r is best;r['pareto']=False
        if best['score']<current:best['pareto']=True;current=best['score']
    write_csv(OUT/'all_shared_models.csv',rows);write_csv(OUT/'shared_frontier.csv',[r for r in rows if r['pareto']])
    ratios=[];allpars=[]
    for r in rows:
        if not r['pareto']:continue
        key=(r['hmask'],r['rmask']);a=ns[key]['parameters'][0];b=ps[key]['parameters'][0]
        for name in PARAMS:ratios.append(dict(complexity=r['complexity_per_population'],HIF=r['HIF'],RhoA=r['RhoA'],parameter=name,nonPACC=a[name],PACC=b[name],ratio=b[name]/a[name]))
        for group in (ns[key],ps[key]):allpars += [dict(complexity=r['complexity_per_population'],HIF=r['HIF'],RhoA=r['RhoA'],**f) for f in group['parameters']]
    write_csv(OUT/'frontier_ratios.csv',ratios);write_csv(OUT/'frontier_parameters.csv',allpars)
    # Verify exposed quadratic reproduces the unchanged loss and direct shared-scale loss.
    checks=[]
    for base in ('N','P'):
        z=np.array(REPS[base]['z']);s,num,den=stats(z,base)
        assert abs(s-REPS[base]['bank_score'])<1e-10
        checks.append(dict(test='Representative score equivalence '+base,error=abs(s-REPS[base]['bank_score']),passed=True))
    test_ratio=1.7;score,pars=joint(ZREF,'amplitude',test_ratio,True);direct=0.
    for base,z,fit in [('N',ZREF[:5],pars[:3]),('P',ZREF[5:],pars[3:])]:
        q=ff.LOSS.quantities(z[MAPPING])[:,:3].copy();q[:,0]*=fit[0]['amplitude']
        for i,owner in enumerate(ff.OWNERS):direct+=np.mean(((q[owner]-ff.TARGETS[base][i])/ff.SCALES[base][i,:,None])**2)/10
    assert abs(score-direct)<1e-10
    checks.append(dict(test='Constrained amplitude ratio loss equivalence',error=abs(score-direct),passed=True));write_csv(OUT/'tests.csv',checks)
    profiles=[];dest=OUT/'profiles.json'
    if dest.exists():profiles=json.load(dest.open())
    done={(r['parameter'],r['ratio']) for r in profiles}
    with ProcessPoolExecutor(max_workers=4,initializer=init) as pool:
        for i,r in enumerate(pool.map(task,[(p,float(v)) for p in PARAMS for v in RATIOS if (p,float(v)) not in done])):
            profiles.append(r)
            if i%8==0:dest.write_text(json.dumps(profiles,indent=2));print('Profiles completed',len(profiles),'/104',flush=True)
    dest.write_text(json.dumps(profiles,indent=2));write_csv(OUT/'ratio_profiles.csv',[{k:v for k,v in r.items() if k!='parameters'} for r in profiles])
    # Predictions for the representative shared mechanism, not population-specific winners.
    preds=[];pars=[]
    for base in ('N','P'):
        _,p,_,_=ff.evaluate_parameters(base,REPS[base]['parameters'],65536,99223)
        preds += [dict(base=base,**r) for r in p];pars+=REPS[base]['parameters']
    write_csv(OUT/'representative_predictions.csv',preds);write_csv(OUT/'representative_parameters.csv',pars)
    sources=[ROOT/'run_shared.py',ROOT/'scale_stats.py',ROOT/'plot.py',ROOT/'validate_equal.py',EX/'fit_frontier.py',PARENT/'run.py',PARENT/'model.py',PARENT.parent/'data.py']
    meta=dict(constraint='Identical HIF/RhoA parameter subsets across populations; population-specific values and change magnitudes retained.',
        combination='Arithmetic mean of matching-subset population scores: each population has five equally weighted source datasets. Continuous fitting separates exactly under this subset-only constraint; no new unrestricted fits required.',
        representative='HIF kappa and RhoA kappa: best shared mechanism at minimum nonzero tested complexity. Presentation choice by parsimony, not statistical model selection.',
        complexity='k denotes KD-specific parameters per population; across both populations there are 2k KD parameters plus eight baseline parameters.',
        profiles='Representative model only. Fix a baseline PACC/nonPACC parameter ratio, jointly refit all other baseline/KD parameters in both populations. Shared-amplitude ratio handled analytically. No CI thresholds or model averaging.',
        profile_range=[.03,30],profile_seeds=[31,53],profile_bank_reference=REF,profile_numerics='Two DE searches plus original-fit local start, existing simulation bank; small negative excesses can reflect numerical optimization variation.',
        limits='Parameter contrasts are conditional on stationary homogeneous model, source-summary weighting and unresolved directional-index normalization. Same subset is not same effect size.',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},input_sha256=hashlib.sha256((EX/'results/fits.json').read_bytes()).hexdigest(),wall_seconds=time.time()-start)
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2));print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
