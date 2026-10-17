"""Independent-seed checks, resolution checks and explicit alternative explanation."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from model import simulate, metrics, msd
from run_analysis import read_data, surface, interval, write_csv
ROOT=Path(__file__).resolve().parent
out=ROOT/'results'
summary=json.loads((out/'summary.json').read_text())
data=read_data()
validation=[]
for study in ('Xuan2018','Mallin2023'):
    r=next(r for r in data if r['study']==study)
    cadence,duration=r['cadence_h'],r['duration_h']
    bank=list(csv.DictReader((out/f'{study}_simulation_bank.csv').open()))
    grid=np.array([float(r['tau_h']) for r in bank])
    for group in ('reference','polyploid'):
        fit=summary[study][group]['best']; t=fit['tau_h']; v=fit['propulsion_um_h']
        base_substeps=max(32,int(np.ceil(8*cadence/t)))
        for substeps in (base_substeps,2*base_substeps):
            p=simulate(v,t,cadence,duration,n=8192,substeps=substeps,seed=91923)
            m=metrics(p,cadence)
            for metric in ('speed','directness','endpoint'):
                original=np.interp(np.log(t),np.log(grid),[float(r[metric]) for r in bank])
                if metric != 'directness': original*=v
                validation.append(dict(study=study,group=group,substeps=substeps,metric=metric,
                                       fitted_bank_prediction=original,independent_mean=float(m[metric].mean()),
                                       independent_MC_SE=float(m[metric].std(ddof=1)/np.sqrt(len(p))),
                                       relative_difference=float(m[metric].mean()/original-1)))
write_csv(out/'numerical_validation.csv',validation)

# Fit shared tau plus different mean speed and a polyploid CV.
# Ref CV is fixed at zero; this is a competing phenomenological explanation.
alternatives=[]
for study in ('Jemaa2023','Jemaa2017','Xuan2018','Mallin2023'):
    bank=list(csv.DictReader((out/f'{study}_simulation_bank.csv').open()))
    knots=np.array([float(r['tau_h']) for r in bank]);tau=np.geomspace(knots[0],knots[-1],321)
    unit={k:np.interp(np.log(tau),np.log(knots),[float(r[k]) for r in bank]) for k in ('speed','directness','endpoint')}
    speeds=np.geomspace(.1,500,1401)
    targets=[[r for r in data if r['study']==study and r['group']==g and r['role']=='fit'] for g in ('reference','polyploid')]
    a,ma=surface(targets[0],tau,unit,speeds)
    winner=None; any_compatible=False
    for cv in np.linspace(0,2.5,51):
        b,mb=surface(targets[1],tau,unit,speeds,cv=cv)
        any_compatible |= bool((ma.any(axis=1)&mb.any(axis=1)).any())
        costs=a.min(axis=1)+b.min(axis=1);i=int(costs.argmin());j=int(a[i].argmin());k=int(b[i].argmin())
        if winner is None or costs[i]<winner['distance_score']:
            winner=dict(study=study,model='speed_and_polyploid_heterogeneity_shared_tau',distance_score=float(costs[i]),
                        reference_speed=float(speeds[j]),polyploid_speed=float(speeds[k]),shared_tau_h=float(tau[i]),
                        reference_cv=0,polyploid_cv=float(cv))
    winner['compatible_with_all_fit_intervals']=any_compatible
    alternatives.append(winner)
write_csv(out/'heterogeneity_model_comparison.csv',alternatives)

# Does any baseline fitting-compatible pair allow the reported crossover?
study='Xuan2018'
bank=list(csv.DictReader((out/f'{study}_simulation_bank.csv').open()))
knots=np.array([float(r['tau_h']) for r in bank]);tau=np.geomspace(knots[0],knots[-1],321)
unit={k:np.interp(np.log(tau),np.log(knots),[float(r[k]) for r in bank]) for k in ('speed','directness','endpoint')}
speeds=np.geomspace(.1,500,1401)
allowed=[]
for group in ('reference','polyploid'):
    targets=[r for r in data if r['study']==study and r['group']==group and r['role']=='fit']
    _,mask=surface(targets,tau,unit,speeds)
    ii,jj=np.where(mask);allowed.append((tau[ii],speeds[jj]))
# Check for crossing within [2,3] by signs at its endpoints for every pair.
ta,va=allowed[0];tb,vb=allowed[1]
d2=msd(2,vb,tb)[None,:]-msd(2,va,ta)[:,None]
d3=msd(3,vb,tb)[None,:]-msd(3,va,ta)[:,None]
pairs=np.argwhere((d2<0)&(d3>0))
result={'fitting_compatible_pairs':int(len(ta)*len(tb)), 'pairs_crossing_between_2_and_3h':int(len(pairs)),
        'interpretation':'Existence/sensitivity check only; no probability assigned to grid points.'}
if len(pairs):
    i,j=pairs[0];result['example']={'reference_tau_h':float(ta[i]),'reference_speed':float(va[i]),'polyploid_tau_h':float(tb[j]),'polyploid_speed':float(vb[j])}
(out/'crossover_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
print('Independent validation and alternative-model analysis complete.')

# Small trajectory illustration using primary fits and matched axes within study.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,2,figsize=(9,9))
for row,study in enumerate(('Xuan2018','Mallin2023')):
    obs=next(r for r in data if r['study']==study)
    limit=0
    for col,group in enumerate(('reference','polyploid')):
        fit=summary[study][group]['best']
        nsub=max(32,int(np.ceil(8*obs['cadence_h']/fit['tau_h'])))
        p=simulate(fit['propulsion_um_h'],fit['tau_h'],obs['cadence_h'],obs['duration_h'],n=24,substeps=nsub,seed=212)
        for track in p: axes[row,col].plot(track[:,0],track[:,1],lw=.8,alpha=.75)
        limit=max(limit,float(np.abs(p).max())*1.05)
        axes[row,col].set(title=f'{study}: {group}',xlabel='x (µm)',ylabel='y (µm)',aspect='equal')
    for ax in axes[row]: ax.set(xlim=(-limit,limit),ylim=(-limit,limit))
fig.suptitle('Simulated cell-centre tracks; same scale within each study\nMallin reference parameters are weakly constrained')
fig.tight_layout();fig.savefig(out/'example_trajectories.png',dpi=160);plt.close(fig)
(out/'validation_metadata.json').write_text(json.dumps({'seed':91923,'cells':8192,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'notes':'Independent-seed simulation; base and twice-base internal resolution; alternative model uses original fitting bank.'},indent=2)+'\n')
