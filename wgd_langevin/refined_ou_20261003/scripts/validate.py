"""Frozen analytical-fit checks; no validation-based parameter selection."""
import json
from datetime import datetime,timezone
import numpy as np
from fit import *
def sampled_predictions(z,n,seed):
    out=predict(z,n,seed)
    pars,_,_=decode(z)
    for i,study in enumerate(STUDIES):
        if study=='Hosny_baseline':continue
        cad,T=SCHEDULE[study]
        for g in (0,1):
            tau,sigma=pars[i,g]
            inc=unit_increments(float(tau),cad,T,n,seed+100*i)*sigma
            path=np.linalg.norm(inc,axis=2).sum(axis=1)
            end=np.linalg.norm(inc.sum(axis=1),axis=1)
            p=out[study][g]
            p.update(speed=float(path.mean()/T),path_length=float(path.mean()),endpoint=float(end.mean()),directness=float((end/path).mean()),early_msd=float((inc[:,0,:]**2).sum(axis=1).mean()))
            # Exact marginal displacement sample at the specified MSD lag.
            v=float(displacement_variance(11.1111111111,np.array([1/tau]),np.array([sigma*sigma*tau/2])))
            normals=normal_bank(1,n,seed+100*i+51).reshape(n,2)
            p['msd']=float(v*(normals**2).sum(axis=1).mean())
            if study=='Xuan2018':
                lags=np.sort(CURVE.lag_h.unique())
                variances=displacement_variance(lags,np.array([1/tau]),np.array([sigma*sigma*tau/2]))
                p['curve_lags_h']=lags.tolist();p['msd_curve']=(variances*(normals**2).sum(axis=1).mean()).tolist()
    return out

def validate():
    fits=[json.loads(p.read_text()) for p in sorted(OUT.glob('fit_*.json'))]
    assert len(fits)==8*len(SCENARIOS),'Incomplete campaign'
    for scenario in SCENARIOS:
        best=min((f for f in fits if f['scenario']==scenario),key=lambda f:f['total_score'])
        z=best['z'];a=analytic(z)
        sampled=[sampled_predictions(z,16384,seed) for seed in (74101,92317)]
        comparisons=[]
        for checks in (False,True):
            _,ra=data_residuals(a,checks)
            rr=[data_residuals(p,checks)[1] for p in sampled]
            for j,row in enumerate(ra):
                held=checks or (row['study']=='Hosny_baseline' and row.get('quantile')!=1)
                # Constant-path approximation provides no path quartile prediction.
                ap=None if row['study']=='Hosny_baseline' and row['metric']=='path_length' and row.get('quantile')!=1 else row['predicted']
                comparisons.append(dict(study=row['study'],group=row['group'],metric=row['metric'],quantile=row.get('quantile',''),obs=row['obs'],analytic=ap,tracks_seed1=rr[0][j]['predicted'],tracks_seed2=rr[1][j]['predicted'],role='correlated curve check' if checks and row['study']=='Xuan2018' and row['metric']=='msd' else 'held out' if held else 'fit',lag_h=row.get('lag_h'),lag_bin=row.get('lag_bin'),objective_weight=row.get('objective_weight'),statistic=row.get('statistic'),id=row.get('id')))
        atomic(OUT/f'validation_{scenario}.json',dict(selected_start=best['start'],z=z,comparison=comparisons,n=16384,seeds=[74101,92317],generated=datetime.now(timezone.utc).isoformat()))
        # Save example tracks with unaltered fitted parameters, including Hosny drift.
        pars,drifts,_=decode(z);tracks={}
        for i,s in enumerate(STUDIES):
            cad,T=SCHEDULE[s]
            for g in (0,1):
                inc=unit_increments(float(pars[i,g,0]),cad,T,32,17001+i)*pars[i,g,1]
                if s=='Hosny_baseline':inc[:,:,1]+=drifts[g]*cad
                tracks[f'{s}_{g}']=np.concatenate([np.zeros((32,1,2)),np.cumsum(inc,axis=1)],axis=1)
        np.savez_compressed(OUT/f'tracks_{scenario}.npz',**tracks)
        unit_increments.cache_clear()
        normal_bank.cache_clear()
    summarize(fits)

def summarize(fits):
    selected={s:min((f for f in fits if f['scenario']==s),key=lambda f:f['total_score']) for s in SCENARIOS}
    atomic(OUT/'selected_fits.json',selected)
    lines=['# Refined OU fit results','', 'Generated '+datetime.now(timezone.utc).isoformat()+'. All 40 bounded starts and frozen trajectory validation finished.','', '|Scenario|Start|Data score|Pooling|Prior|Total|Converged starts|','|---|---:|---:|---:|---:|---:|---:|']
    for s,b in selected.items():
        n=sum(f['success'] for f in fits if f['scenario']==s)
        lines.append(f"|{s}|{b['start']}|{b['data_score']:.5g}|{b['pooling_penalty']:.5g}|{b['persistence_prior_penalty']:.5g}|{b['total_score']:.5g}|{n}/8|")
    b=selected['primary'];pars,drift,_=decode(b['z'])
    lines+=['','Primary estimates (hours; noise in micrometers/hour^(3/2)):','', '|Study|tau reference|tau WGD|tau fold|noise reference|noise WGD|noise fold|','|---|---:|---:|---:|---:|---:|---:|']
    for i,study in enumerate(STUDIES):
        a,c=pars[i];lines.append(f'|{study}|{a[0]:.5g}|{c[0]:.5g}|{c[0]/a[0]:.5g}|{a[1]:.5g}|{c[1]:.5g}|{c[1]/a[1]:.5g}|')
    lines+=['',f'Hosny signed drifts: {drift[0]:.5g}, {drift[1]:.5g} micrometers/hour.','', 'Xuan curve deviations (log ratio RMS; each regime contributes equally):','', '|Group|Regime|Points|RMS log(predicted/observed)|Mean log ratio|Inside raster envelope|','|---|---|---:|---:|---:|---:|']
    for g in (0,1):
        for regime in (0,1,2):
            r=[x for x in b['rows'] if x['metric']=='msd_curve' and x['group']==g and x['lag_bin']==regime]
            v=np.log([x['predicted']/x['obs'] for x in r]);lines.append(f"|{g}|{['early','intermediate','late'][regime]}|{len(r)}|{np.sqrt(np.mean(v*v)):.5g}|{v.mean():.5g}|{sum(x['within'] for x in r)}/{len(r)}|")
    d=b['diagnostics'];lines+=['',f"Primary regularized Jacobian rank: {d['local_rank']}/22. Active numerical bounds: {d['active_bounds']}. These are numerical sensitivity diagnostics, not uncertainty intervals.",'','Validation files retain two frozen 16,384-track Sobol ensembles per scenario; all fitting curves are analytical. Sampled Xuan MSD validation uses exact Gaussian marginal displacements at every curve lag with common standardized draws, so agreement tests the marginal formula, not independent full-path correlations. Movie trajectory diagnostics are separate. Hosny nonmedian path predictions are undefined under the analytical approximation and only sampled trajectories are used for those checks. Mallin central statistic remains unspecified; the historical mean approximation is retained. Early 600-second Xuan checks overlap the fitted curve and are correlated checks.','', 'Objective tolerances are unchanged historical descriptive scales recentered on reviewed values; raster bounds and SEM are retained separately. Xuan native digitized points are normalized within three equally spaced log-lag regimes; each regime contributes one third of a group curve block, and each group curve block has total weight one, matching one speed target. Raster envelopes are not population confidence intervals.']
    (OUT/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    atomic(OUT/'status.json',dict(generated=datetime.now(timezone.utc).isoformat(),state='complete',fits=40,validations=5))
if __name__=='__main__':validate()
