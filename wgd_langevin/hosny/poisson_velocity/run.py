"""Simulation-bank optimization of jointly constrained Poisson-kick populations."""
from pathlib import Path
import sys,csv,json,time,hashlib,argparse
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import differential_evolution
from concurrent.futures import ProcessPoolExecutor
from model import simulate,metrics,moments,msd,rayleigh_p
from tests import run_tests
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT.parent))
from data import DATA,BYID,Q,METRICS
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
GROUPS={'O':['O'],'N':['N','HN0','RN0'],'P':['P','HP0','RP0'],'HN':['HN1'],'HP':['HP1'],'RN':['RN1'],'RP':['RP1']}
TAU=np.r_[np.geomspace(.005,.04,7)[:-1],np.geomspace(.04,16,21)];RATE=np.geomspace(.12,24,21);KAPPA=np.r_[0,np.geomspace(.01,8,24)]
KM=('path','directness','parallel','endpoint','perpendicular','angle_parallel')
LT=(float(np.log(TAU[0])),float(np.log(TAU[-1])));LR=(float(np.log(RATE[0])),float(np.log(RATE[-1])));KC=(0.,float(np.log1p(KAPPA[-1])))
MODELS=['primary','HIF_memory','RhoA_direction','both','full']

def write_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def bank_task(args):
    it,ir,n=args;out=np.empty((len(KAPPA),len(KM),3))
    for ik,k in enumerate(KAPPA):
        m=metrics(simulate(TAU[it],RATE[ir],1,k,n=n,seed=221001))
        out[ik]=[np.quantile(m[key],Q) for key in KM]
    return it,ir,out

def bank(n,workers):
    dest=OUT/'bank.npz'
    qq=np.empty((len(TAU),len(RATE),len(KAPPA),len(KM),3))
    tasks=[(i,j,n) for i in range(len(TAU)) for j in range(len(RATE))]
    if dest.exists():
        b=np.load(dest)
        assert int(b['cells'])==n and np.array_equal(b['rate'],RATE) and np.array_equal(b['kappa'],KAPPA)
        existing={float(t):i for i,t in enumerate(b['tau'])}
        for i,t in enumerate(TAU):
            if float(t) in existing:qq[i]=b['q'][existing[float(t)]]
        tasks=[task for task in tasks if float(TAU[task[0]]) not in existing]
        if not tasks:return qq
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for count,(it,ir,q) in enumerate(pool.map(bank_task,tasks)):
            qq[it,ir]=q
            if count%21==0:print(f'Bank: {count+1}/{len(tasks)} new tau/rate slices',flush=True)
    np.savez_compressed(dest,q=qq,tau=TAU,rate=RATE,kappa=KAPPA,cells=n)
    return qq

class Loss:
    def __init__(self,qq):
        self.interp=RegularGridInterpolator((np.log(TAU),np.log(RATE),np.log1p(KAPPA)),qq,bounds_error=False,fill_value=np.nan)
        self.targets={g:np.array([[BYID[id][k] for k in KM[:3]] for id in ids]) for g,ids in GROUPS.items()}
        self.scales={g:np.maximum(t[:,:,2]-t[:,:,0],[1,.1,.1]) for g,t in self.targets.items()}
    def quantities(self,shapes):return self.interp(np.asarray(shapes))
    def evaluate(self,groups,shapes,amplitude_ids,mode='standard',return_fits=False):
        q=self.quantities(shapes)
        if not np.isfinite(q).all():return (1e6,[]) if return_fits else 1e6
        weights=np.array([0,1,0]) if mode=='medians' else np.ones(3)
        amps={}
        for a in set(amplitude_ids):
            num=den=0.
            for i,g in enumerate(groups):
                if amplitude_ids[i]!=a:continue
                t=self.targets[g][:,0];scale=self.scales[g][:,0]
                num+=np.sum(weights*q[i,0]*t/scale[:,None]**2)
                den+=np.sum(weights*q[i,0]**2/scale[:,None]**2)
            amps[a]=float(num/den)
        total=0.;fits=[]
        for i,g in enumerate(groups):
            amp=amps[amplitude_ids[i]];pred=q[i,:3].copy();pred[0]*=amp
            if mode=='angle':pred[2]=q[i,5]
            res=((pred-self.targets[g])/self.scales[g][:,:,None])**2
            if mode=='no_direction':res=res[:,:2]
            score=float(np.sum(res*weights)/(weights.sum()*res.shape[1]))
            total+=score
            if return_fits:
                tau,lam,k=np.exp(shapes[i][0]),np.exp(shapes[i][1]),np.expm1(shapes[i][2]);u,cv=moments(tau,lam,amp,k)
                fits.append(dict(group=g,tau=float(tau),gamma=float(1/tau),rate=float(lam),amplitude=amp,kappa=float(k),mean_forward_velocity=float(u),velocity_sd_total=float(np.sqrt(cv)),score=score/len(GROUPS[g]),source_ids=';'.join(GROUPS[g])))
        score=total/sum(len(GROUPS[g]) for g in groups)
        return (score,fits) if return_fits else score

def specification(base,model):
    groups=[base,'H'+base,'R'+base]
    if model=='full':
        bounds=[LT,LR,KC]*3
        def decode(z):return np.array(z).reshape(3,3),[0,1,2]
    else:
        # HIF concentration and RhoA tau are monotonic fractions of baseline.
        # Tau fraction is in log coordinates between the numerical lower bound and baseline.
        bounds=[LT,LR,KC,(0.,1.),(0.,1.)]
        if model in ('HIF_memory','both'):bounds.append(LT)
        if model in ('RhoA_direction','both'):bounds.append(KC)
        def decode(z):
            lt,lr,k,hfrac,rfrac=z[:5]
            sh=np.array([[lt,lr,k],[lt,lr,np.log1p(np.expm1(k)*hfrac)],[LT[0]+rfrac*(lt-LT[0]),lr,k]])
            j=5
            if model in ('HIF_memory','both'):sh[1,0]=z[j];j+=1
            if model in ('RhoA_direction','both'):sh[2,2]=z[j]
            return sh,[0,0,0]
    return groups,bounds,decode

def fit(loss,base,model,mode='standard',fixed_rate=None,fixed_tau=None,seeds=(17,29,43),starts=()):
    if base=='O':
        groups=['O'];bounds=[LT,LR]
        def decode(z):return [[z[0],z[1],0]],[0]
    else:groups,bounds,decode=specification(base,model)
    bounds=list(bounds)
    if fixed_rate is not None:bounds[1]=(np.log(fixed_rate),np.log(fixed_rate))
    if fixed_tau is not None:bounds[0]=(np.log(fixed_tau),np.log(fixed_tau))
    def objective(z):
        sh,aa=decode(z);return loss.evaluate(groups,sh,aa,mode)
    candidates=[]
    for z in starts:
        if len(z)==len(bounds):candidates.append((objective(z),np.array(z),True))
    for seed in seeds:
        rr=differential_evolution(objective,bounds,seed=seed,popsize=12,maxiter=240 if model=='full' else 160,tol=2e-6,polish=True)
        candidates.append((float(rr.fun),rr.x,bool(rr.success)))
    score,z,success=min(candidates,key=lambda x:x[0]);sh,aa=decode(z)
    score,fits=loss.evaluate(groups,sh,aa,mode,True)
    return dict(base=base,model=model,score=score,parameters=fits,z=z.tolist(),optimizer_converged=success,mode=mode)

def embedded_start(prev,newmodel):
    p=prev['parameters'];z=prev['z'][:5]
    if newmodel in ('HIF_memory','both'):z.append(np.log(p[1]['tau']))
    if newmodel in ('RhoA_direction','both'):z.append(np.log1p(p[2]['kappa']))
    if newmodel=='full':z=[val for f in p for val in (np.log(f['tau']),np.log(f['rate']),np.log1p(f['kappa']))]
    return z

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cells',type=int,default=2048);ap.add_argument('--workers',type=int,default=4);args=ap.parse_args()
    start=time.time();write_csv(OUT/'tests.csv',run_tests());loss=Loss(bank(args.cells,args.workers))
    primary={};comparisons=[];allfits={};optim=[]
    for base in ('O','N','P'):
        r=fit(loss,base,'primary');primary[base]=r;comparisons.append(r);allfits.update({f['group']:f for f in r['parameters']});optim.append(r)
        print('PRIMARY',json.dumps(r),flush=True)
        (OUT/'optimization_partial.json').write_text(json.dumps(optim,indent=2))
        if base!='O':
            previous=r
            for model in MODELS[1:]:
                starts=[embedded_start(previous,model)]
                # Both expansions contain the primary model and either one-parameter extension.
                if model=='both':starts=[embedded_start(v,model) for v in optim if v['base']==base and v['model'] in ('primary','HIF_memory','RhoA_direction')]
                if model=='full':starts=[embedded_start(v,model) for v in optim if v['base']==base]
                rr=fit(loss,base,model,starts=starts);optim.append(rr);comparisons.append(rr)
                if model in ('HIF_memory','RhoA_direction'):previous=r
                else:previous=rr
                print('EXPANSION',base,model,rr['score'],flush=True)
                (OUT/'optimization_partial.json').write_text(json.dumps(optim,indent=2))
    (OUT/'optimization.json').write_text(json.dumps(optim,indent=2))
    write_csv(OUT/'fits.csv',list(allfits.values()))
    write_csv(OUT/'model_comparison.csv',[dict(base=r['base'],model=r['model'],score=r['score'],n_parameters=(3 if r['base']=='O' else {'primary':6,'HIF_memory':7,'RhoA_direction':7,'both':8,'full':12}[r['model']]),optimizer_converged=r['optimizer_converged']) for r in comparisons])
    # Independent simulations validate every compared fit, not just the selected primary.
    predictions=[];validation=[];tracks={};samples={};ind_scores=[]
    for r in optim:
        score=0.;count=0
        for f in r['parameters']:
            g=f['group'];p,v=simulate(f['tau'],f['rate'],f['amplitude'],f['kappa'],n=32768,seed=77019,return_velocity=True)
            m=metrics(p);qq=np.array([np.quantile(m[k],Q) for k in KM]);expected=loss.quantities([[np.log(f['tau']),np.log(f['rate']),np.log1p(f['kappa'])]])[0];expected[[0,3]]*=f['amplitude']
            for id in GROUPS[g]:
                target=np.array([BYID[id][k] for k in KM[:3]]);scale=np.maximum(target[:,2]-target[:,0],[1,.1,.1]);score+=float(np.mean(((qq[:3]-target)/scale[:,None])**2));count+=1
            for ik,k in enumerate(KM):
                for j,quant in enumerate(Q):validation.append(dict(base=r['base'],model=r['model'],group=g,metric=k,quantile=quant,bank=float(expected[ik,j]),independent=float(qq[ik,j])))
            if r['model']=='primary':
                tracks[g]=p[:24];tracks[g+'_velocity']=v[:6];samples[g]=m
                for ik,k in enumerate(KM):
                    for j,quant in enumerate(Q):predictions.append(dict(group=g,metric=k,quantile=quant,predicted=float(qq[ik,j]),used_in_fit=ik<3))
        ind_scores.append(dict(base=r['base'],model=r['model'],bank_score=r['score'],independent_score=score/count))
    write_csv(OUT/'predictions.csv',predictions);write_csv(OUT/'numerical_checks.csv',validation);write_csv(OUT/'independent_scores.csv',ind_scores)
    np.savez_compressed(OUT/'tracks.npz',**tracks)
    print('Independent prediction validation complete',flush=True)
    # Profiles jointly reoptimize all other primary parameters at each fixed baseline rate.
    profiles=[]
    for base in ('O','N','P'):
        for rate in np.geomspace(RATE[0],RATE[-1],13):
            z=np.array(primary[base]['z']);z[1]=np.log(rate)
            r=fit(loss,base,'primary',fixed_rate=rate,seeds=(61,),starts=[z])
            f=r['parameters'][0];profiles.append(dict(base=base,rate=rate,score=r['score'],amplitude=f['amplitude'],tau=f['tau'],kappa=f['kappa']))
        print('Rate profile complete',base,flush=True)
    write_csv(OUT/'rate_profiles.csv',profiles)
    memory=[]
    for base in ('O','N','P'):
        for tau in np.geomspace(TAU[0],4,13):
            z=np.array(primary[base]['z']);z[0]=np.log(tau)
            r=fit(loss,base,'primary',fixed_tau=tau,seeds=(71,),starts=[z])
            f=r['parameters'][0];memory.append(dict(base=base,tau=tau,rate=f['rate'],amplitude=f['amplitude'],score=r['score']))
        print('Memory profile complete',base,flush=True)
    write_csv(OUT/'memory_profiles.csv',memory)
    sensitivity=[]
    for base in ('N','P'):
        for mode in ('angle','medians','no_direction'):
            r=fit(loss,base,'primary',mode=mode,seeds=(41,),starts=[primary[base]['z']]);sensitivity.append(r)
    (OUT/'sensitivity.json').write_text(json.dumps(sensitivity,indent=2))
    ray=[];sub=[];timepred=[];counter=[]
    for g,f in allfits.items():
        m=samples[g];rng=np.random.default_rng(131);ps=rayleigh_p(m['angle'][rng.integers(0,len(m['angle']),(3000,60))])
        for id in GROUPS[g]:ray.append(dict(source=id,group=g,reported_p=BYID[id]['rayleigh'],predicted_p05=float(np.quantile(ps,.05)),predicted_p50=float(np.median(ps)),predicted_p95=float(np.quantile(ps,.95)),fraction_p_below_05=float(np.mean(ps<.05))))
        for t in np.geomspace(.1,32,100):timepred.append(dict(group=g,time_h=t,msd=float(msd(t,f['tau'],f['rate'],f['amplitude'],f['kappa'])),mean_forward_displacement=f['mean_forward_velocity']*t,velocity_correlation=float(np.exp(-t/f['tau']))))
        if g in ('N','P'):
            m0=metrics(simulate(f['tau'],f['rate'],f['amplitude'],0,n=32768,seed=810));counter.append(dict(group=g,directness=float(np.median(m0['directness'])),path=float(np.median(m0['path'])),endpoint=float(np.median(m0['endpoint']))))
    threshold=float(np.quantile(samples['P']['endpoint'],.75))
    for g in ('P','HP','RP'):sub.append(dict(group=g,threshold=threshold,positive_high_fraction=float(np.mean((samples[g]['endpoint']>threshold)&(samples[g]['parallel']>0)))))
    write_csv(OUT/'rayleigh_checks.csv',ray);write_csv(OUT/'qualitative_subpopulation_check.csv',sub);write_csv(OUT/'time_predictions.csv',timepred);write_csv(OUT/'cue_removed.csv',counter)
    sources=[ROOT/'model.py',ROOT/'run.py',ROOT/'tests.py',ROOT.parent/'data.py',ROOT.parents[2]/'literature_review/hosny_2026.md']
    meta=dict(arguments=vars(args),groups=GROUPS,parameters='tau=1/gamma (h), rate (/h), Rayleigh kick scale (um/h), von Mises concentration (dimensionless)',
        model='dv=-v dt/tau plus independent Poisson velocity jumps; dr=v dt; Rayleigh magnitudes and von Mises angles from +y',
        initialization='Near stationary: integrate past kicks to 14 tau plus the first older event; omitted mean tail bounded by exp(-14), variance by exp(-28)',
        sampling='Exact event integration; 16 h tracks sampled each 0.5 h',
        primary='Joint baseline + HIF + RhoA. HIF concentration <= baseline; RhoA tau <= baseline. Rate and amplitude shared across all three. Normoxia kappa=0.',
        expansions='HIF_memory adds an unconstrained HIF tau; RhoA_direction adds unconstrained RhoA kappa; both includes both; full releases all KD sharing and monotonic restrictions. All independently jointly refitted.',
        loss='Per-dataset mean over 3 metrics x 3 quartiles of squared residual / observed IQR^2; floors 1 um and 0.1 for fractions; equal dataset weights; shared amplitude analytically profiled.',
        normalization='Parallel displacement / path assumed; angle-normalization sensitivity retained.',
        bounds=dict(tau=TAU[[0,-1]].tolist(),rate=RATE[[0,-1]].tolist(),kappa=KAPPA[[0,-1]].tolist()),
        optimization='Trilinear interpolation in log tau, log rate, log(1+kappa); three DE seeds plus polish; nested starting candidates retained; independent 32768-cell validation.',
        caveats=['Descriptive loss, no likelihood/CI or complexity penalty.','Approximate literature quantiles, no raw tracks or replicate covariance.','Constant local oxygen direction; no spatial/temporal response or cell area model.','Rayleigh checks assume 60 independent cells; biological clustering unknown.'],
        wall_seconds=time.time()-start,source_sha256={str(p.relative_to(ROOT.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2))
    print('DONE',time.time()-start,flush=True)
if __name__=='__main__':main()
