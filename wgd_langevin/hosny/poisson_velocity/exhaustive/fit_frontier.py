"""Exhaustive joint HIF/RhoA parameter-sharing models; no change-sign restrictions."""
from pathlib import Path
import sys,json,csv,hashlib,time,argparse
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from numba import njit
from scipy.optimize import differential_evolution,minimize
ROOT=Path(__file__).resolve().parent;PARENT=ROOT.parent
sys.path.insert(0,str(PARENT))
from run import Loss,LT,LR,KC,GROUPS,KM,write_csv
from model import simulate,metrics,moments,msd,rayleigh_p
from tests import run_tests
from data import BYID,Q
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
PARAMS=('tau','rate','amplitude','kappa');SYMBOLS=('tau','lambda','A','kappa')
SHAPE_PARAMS=(0,1,3);BOUNDS=(LT,LR,KC)
PAIRS=[(0,0)]+[(h,r) for h in range(1,16) for r in range(1,16)]

@njit(cache=True)
def fast_objective(z,mapping,ampids,bank,xgrid,ygrid,kgrid,targets,scales,owners):
    qq=np.zeros((3,3,3))
    for g in range(3):
        x=z[mapping[g,0]];y=z[mapping[g,1]];k=z[mapping[g,2]]
        ix=min(max(np.searchsorted(xgrid,x)-1,0),len(xgrid)-2)
        iy=min(max(np.searchsorted(ygrid,y)-1,0),len(ygrid)-2)
        ik=min(max(np.searchsorted(kgrid,k)-1,0),len(kgrid)-2)
        ax=(x-xgrid[ix])/(xgrid[ix+1]-xgrid[ix]);ay=(y-ygrid[iy])/(ygrid[iy+1]-ygrid[iy]);ak=(k-kgrid[ik])/(kgrid[ik+1]-kgrid[ik])
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    weight=(ax if a else 1-ax)*(ay if b else 1-ay)*(ak if c else 1-ak)
                    for m in range(3):
                        for q in range(3):qq[g,m,q]+=weight*bank[ix+a,iy+b,ik+c,m,q]
    num=np.zeros(3);den=np.zeros(3)
    for d in range(len(owners)):
        g=owners[d];aid=ampids[g]
        for q in range(3):
            num[aid]+=qq[g,0,q]*targets[d,0,q]/scales[d,0]**2
            den[aid]+=qq[g,0,q]**2/scales[d,0]**2
    amps=np.ones(3)
    for a in range(3):
        if den[a]>0:amps[a]=num[a]/den[a]
    total=0.
    for d in range(len(owners)):
        g=owners[d]
        for m in range(3):
            for q in range(3):
                pred=qq[g,m,q]*(amps[ampids[g]] if m==0 else 1.)
                total+=((pred-targets[d,m,q])/scales[d,m])**2
    return total/(len(owners)*9)

def schema(h,r):
    bounds=list(BOUNDS);mapping=np.tile(np.arange(3),(3,1))
    for g,mask in [(1,h),(2,r)]:
        for j,bit in enumerate(SHAPE_PARAMS):
            if mask&(1<<bit):mapping[g,j]=len(bounds);bounds.append(BOUNDS[j])
    ampids=np.array([0,1 if h&4 else 0,2 if r&4 else 0],dtype=np.int64)
    return bounds,mapping,ampids

def encode(parameters,h,r):
    p=parameters
    shapes=np.array([[np.log(f['tau']),np.log(f['rate']),np.log1p(f['kappa'])] for f in p])
    bounds,mapping,ampids=schema(h,r);z=np.empty(len(bounds));z[:3]=shapes[0]
    for g in (1,2):
        for j in range(3):
            if mapping[g,j]>=3:z[mapping[g,j]]=shapes[g,j]
    return z

def subset(mask):return '+'.join(SYMBOLS[i] for i in range(4) if mask&(1<<i)) or 'none'
def key(base,h,r):return f'{base}_H{h:02d}_R{r:02d}'

def init_worker():
    global BANK,XGRID,YGRID,KGRID,LOSS,TARGETS,SCALES,OWNERS
    b=np.load(PARENT/'results/bank.npz');BANK=b['q'];XGRID=np.log(b['tau']);YGRID=np.log(b['rate']);KGRID=np.log1p(b['kappa'])
    LOSS=Loss(BANK);TARGETS={};SCALES={};OWNERS=np.array([0,0,0,1,2],dtype=np.int64)
    for base in ('N','P'):
        ids=GROUPS[base]+GROUPS['H'+base]+GROUPS['R'+base]
        TARGETS[base]=np.array([[BYID[id][m] for m in KM[:3]] for id in ids])
        SCALES[base]=np.maximum(TARGETS[base][:,:,2]-TARGETS[base][:,:,0],[1,.1,.1])
    bounds,mapping,aa=schema(8,1)
    fast_objective(np.array([(a+b)/2 for a,b in bounds]),mapping,aa,BANK,XGRID,YGRID,KGRID,TARGETS['N'],SCALES['N'],OWNERS)

def objective_for(base,h,r):
    bounds,mapping,aa=schema(h,r)
    def obj(z):return fast_objective(z,mapping,aa,BANK,XGRID,YGRID,KGRID,TARGETS[base],SCALES[base],OWNERS)
    return bounds,mapping,aa,obj

def fit_task(task):
    base,h,r,parents,seeds=task;bounds,mapping,aa,objective=objective_for(base,h,r)
    candidates=[];runs=[]
    for p in parents:
        z=encode(p['parameters'],h,r)
        # A valid smaller-model fit remains an upper bound on the released model loss.
        candidates.append((objective(z),z,'inherited',False))
        rr=minimize(objective,z,bounds=bounds,method='L-BFGS-B',options={'maxiter':180,'ftol':1e-12})
        candidates.append((float(rr.fun),rr.x,'local',bool(rr.success)))
    for seed in seeds:
        rr=differential_evolution(objective,bounds,seed=seed,popsize=12,maxiter=260,tol=2e-7,polish=True)
        candidates.append((float(rr.fun),rr.x,'DE',bool(rr.success)))
        runs.append(dict(seed=seed,score=float(rr.fun),converged=bool(rr.success),iterations=int(rr.nit)))
    score,z,method,converged=min(candidates,key=lambda v:v[0])
    actual,parameters=LOSS.evaluate([base,'H'+base,'R'+base],z[mapping],aa,return_fits=True)
    assert abs(score-actual)<1e-10,(score,actual)
    return dict(id=key(base,h,r),base=base,hmask=h,rmask=r,HIF=subset(h),RhoA=subset(r),complexity=h.bit_count()+r.bit_count(),total_parameters=4+h.bit_count()+r.bit_count(),bank_score=actual,parameters=parameters,z=z.tolist(),selected_method=method,selected_run_converged=converged,runs=runs)

def evaluate_parameters(base,parameters,n,seed,save_samples=False):
    total=0.;pred=[];samples={};tracks={}
    for f in parameters:
        g=f['group'];p,v=simulate(f['tau'],f['rate'],f['amplitude'],f['kappa'],n=n,seed=seed,return_velocity=True)
        m=metrics(p);q=np.array([np.quantile(m[k],Q) for k in KM])
        for id in GROUPS[g]:
            target=np.array([BYID[id][k] for k in KM[:3]]);scale=np.maximum(target[:,2]-target[:,0],[1,.1,.1]);total+=float(np.mean(((q[:3]-target)/scale[:,None])**2))
        for im,k in enumerate(KM):
            for iq,quant in enumerate(Q):pred.append(dict(group=g,metric=k,quantile=quant,predicted=float(q[im,iq])))
        if save_samples:samples[g]=m;tracks[g]=p[:16];tracks[g+'_velocity']=v[:4]
    return total/5,pred,samples,tracks

def validate_task(args):
    fit,n,seed=args;score,pred,_,_=evaluate_parameters(fit['base'],fit['parameters'],n,seed)
    return dict(id=fit['id'],independent_score=score,predictions=pred)

def checks():
    init_worker();rng=np.random.default_rng(414);err=0.
    for h,r in PAIRS:
        bounds,mapping,aa,obj=objective_for('P',h,r)
        z=np.array([rng.uniform(a,b) for a,b in bounds]);expected,_=LOSS.evaluate(['P','HP','RP'],z[mapping],aa,return_fits=True)
        err=max(err,abs(obj(z)-expected))
    assert err<1e-9,err
    rows=run_tests();rows.append(dict(test='Accelerated loss agrees with reference for all 226 sharing structures',actual=err,expected=0,tolerance=1e-9,passed=True))
    write_csv(OUT/'tests.csv',rows)

def summarize(fits,refs):
    for base in ('N','P'):
        best_so_far=np.inf
        for k in (0,2,3,4,5,6,7,8):
            rows=[f for f in fits if f['base']==base and f['complexity']==k]
            win=min(rows,key=lambda r:r['independent_score'])
            for r in rows:r['best_at_complexity']=r['id']==win['id'];r['pareto']=False
            if win['independent_score']<best_so_far:win['pareto']=True;best_so_far=win['independent_score']
    keys=['id','base','hmask','rmask','HIF','RhoA','complexity','total_parameters','bank_score','independent_score','best_at_complexity','pareto','selected_method','selected_run_converged']
    write_csv(OUT/'all_models.csv',[{k:r[k] for k in keys} for r in fits])
    write_csv(OUT/'frontier.csv',[{k:r[k] for k in keys} for r in fits if r['pareto']])
    write_csv(OUT/'best_at_complexity.csv',[{k:r[k] for k in keys} for r in fits if r['best_at_complexity']])
    write_csv(OUT/'original_reference.csv',[dict(base=r['base'],bank_score=r['score'],independent_score=r['independent_score']) for r in refs])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--validation-cells',type=int,default=32768);args=ap.parse_args();start=time.time()
    # Validate provenance before reusing the established simulation bank.
    oldmeta=json.load((PARENT/'results/metadata.json').open())
    repo=PARENT.parents[2]
    for relative,h in oldmeta['source_sha256'].items():assert hashlib.sha256((repo/relative).read_bytes()).hexdigest()==h,relative
    checks();old=json.load((PARENT/'results/optimization.json').open());refs=[r for r in old if r['model']=='primary' and r['base'] in ('N','P')]
    fits=[];cache=OUT/'fits.json'
    if cache.exists():fits=json.load(cache.open())
    done={(f['base'],f['hmask'],f['rmask']):f for f in fits}
    with ProcessPoolExecutor(max_workers=args.workers,initializer=init_worker) as pool:
        for complexity in (0,2,3,4,5,6,7,8):
            tasks=[]
            for base in ('N','P'):
                for h,r in PAIRS:
                    if h.bit_count()+r.bit_count()!=complexity or (base,h,r) in done:continue
                    parents=[]
                    # Immediate tested subsets embed exactly; old fits provide additional starts.
                    candidates=[f for f in fits if f['base']==base and f['hmask']&h==f['hmask'] and f['rmask']&r==f['rmask']]
                    immediate=[f for f in candidates if f['complexity']==complexity-1]
                    if not immediate and candidates:immediate=[min(candidates,key=lambda f:f['bank_score'])]
                    parents+=immediate
                    parents += [f for f in old if f['base']==base and f['model'] in ('primary','full')]
                    tasks.append((base,h,r,parents,(17,29,43)))
            for f in pool.map(fit_task,tasks):
                fits.append(f);done[(f['base'],f['hmask'],f['rmask'])]=f
                cache.write_text(json.dumps(fits,indent=2))
                if len(fits)%20==0:print(f'Fitted {len(fits)}/452 joint models',flush=True)
            print('Completed complexity',complexity,'models',len(fits),flush=True)
        # Test all subset/superset relations, not just adjacent chosen winners.
        violations=[]
        for base in ('N','P'):
            rr=[r for r in fits if r['base']==base]
            for a in rr:
                for b in rr:
                    if a['hmask']&b['hmask']==a['hmask'] and a['rmask']&b['rmask']==a['rmask'] and b['bank_score']>a['bank_score']+1e-9:violations.append((a['id'],b['id']))
        assert not violations,violations[:10]
        pending=[(f,args.validation_cells,90221) for f in fits if 'independent_score' not in f]
        for i,res in enumerate(pool.map(validate_task,pending)):
            f=next(f for f in fits if f['id']==res['id']);f.update({k:v for k,v in res.items() if k!='id'})
            if i%20==0:cache.write_text(json.dumps(fits,indent=2));print(f'Independent simulation {i+1}/{len(pending)}',flush=True)
    cache.write_text(json.dumps(fits,indent=2))
    for r in refs:r['independent_score'],r['predictions'],_,_=evaluate_parameters(r['base'],r['parameters'],args.validation_cells,90221)
    (OUT/'original_reference.json').write_text(json.dumps(refs,indent=2))
    summarize(fits,refs)
    # Additional independent Monte Carlo replicate for every best-at-complexity model.
    repeat=[]
    for f in fits:
        if f['best_at_complexity']:
            score,_,_,_=evaluate_parameters(f['base'],f['parameters'],args.validation_cells,90222)
            repeat.append(dict(id=f['id'],base=f['base'],complexity=f['complexity'],first_score=f['independent_score'],repeat_score=score))
    write_csv(OUT/'frontier_repeat.csv',repeat)
    # Detailed panels use the best one-parameter-per-KD model, a fair comparison to the original hypothesis.
    detailed=[];tracks={};sub=[];ray=[];times=[];parameters=[];counter=[]
    for base in ('N','P'):
        selected=next(f for f in fits if f['base']==base and f['complexity']==2 and f['best_at_complexity'])
        for label,f in [('best_two',selected),('original',next(r for r in refs if r['base']==base))]:
            score,pred,samples,tr=evaluate_parameters(base,f['parameters'],65536,99223,True)
            detailed += [dict(model=label,base=base,**r) for r in pred]
            parameters += [dict(model=label,base=base,**p) for p in f['parameters']]
            if label=='best_two':
                tracks.update(tr)
                for p in f['parameters']:
                    g=p['group'];m=samples[g];rng=np.random.default_rng(131);ps=rayleigh_p(m['angle'][rng.integers(0,len(m['angle']),(3000,60))])
                    for id in GROUPS[g]:ray.append(dict(group=g,source=id,reported_p=BYID[id]['rayleigh'],predicted_significance_fraction=float(np.mean(ps<.05))))
                    for t in np.geomspace(.1,32,100):times.append(dict(group=g,time_h=t,msd=float(msd(t,p['tau'],p['rate'],p['amplitude'],p['kappa']))))
                if base=='P':
                    threshold=float(np.quantile(samples['P']['endpoint'],.75))
                    for g in ('P','HP','RP'):sub.append(dict(group=g,threshold=threshold,positive_high_fraction=float(np.mean((samples[g]['endpoint']>threshold)&(samples[g]['parallel']>0)))))
                p=f['parameters'][0];m0=metrics(simulate(p['tau'],p['rate'],p['amplitude'],0,n=32768,seed=984))
                counter.append(dict(group=base,directness=float(np.median(m0['directness']))))
    write_csv(OUT/'detailed_predictions.csv',detailed);write_csv(OUT/'detailed_parameters.csv',parameters);np.savez_compressed(OUT/'tracks.npz',**tracks)
    write_csv(OUT/'time_predictions.csv',times);write_csv(OUT/'rayleigh_checks.csv',ray);write_csv(OUT/'qualitative_subpopulation_check.csv',sub);write_csv(OUT/'cue_removed.csv',counter)
    (OUT/'fits.json').write_text(json.dumps(fits,indent=2))
    sources=[ROOT/'fit_frontier.py',PARENT/'model.py',PARENT/'run.py',PARENT.parent/'data.py',repo/'literature_review/hosny_2026.md']
    meta=dict(arguments=vars(args),model_count_per_class=226,subset_search='All 15 nonempty HIF masks x 15 nonempty RhoA masks, plus joint zero-change null. One-KD-only null models not included.',
        fitting='Joint baseline + HIF + RhoA; all selected parameters can increase or decrease; nonselected parameters equal baseline; each source dataset weighted equally.',
        objective='Unchanged mean squared quartile residual / reported IQR^2 across path, directness and forward/path; 16h sampled at 0.5h; same supplied observations.',
        optimization='Same established 2048-cell simulation bank; accelerated numerically equivalent loss; 3 differential evolution seeds plus polishing and embedded subset starts.',
        complexity='Number of KD-specific scalar parameters in addition to four baseline parameters; independently profiled amplitude(s) counted.',
        frontier='Nondominated by independent 32768-cell error and complexity; best per complexity also retained; descriptive, no likelihood, experimental holdout or CI.',
        detailed_panels='Best model at complexity 2 (one parameter per KD) for each class, displayed against original constrained hypothesis with 65536-cell predictions. Not an overall model selection claim.',
        numerical_checks='All subset/superset bank loss ordering checked; separate repeat simulations for all best-at-complexity models.',
        model_limits='Homogeneous stationary populations; FMI normalization ambiguous; sub-frame tau/amplitude tradeoff persists; no spatial oxygen response or area selection model.',
        source_sha256={str(p.relative_to(repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},bank_sha256=hashlib.sha256((PARENT/'results/bank.npz').read_bytes()).hexdigest(),wall_seconds=time.time()-start)
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2));print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
