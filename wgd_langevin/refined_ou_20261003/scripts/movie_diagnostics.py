"""Representative Xuan movie checks; never contributes fitting residuals.

Exact integrated stationary Gaussian OU draws at the actual 73 frame intervals.
Centroid allowance is an assumed per-coordinate Gaussian SD for sensitivity only.
"""
import argparse, csv, hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import toeplitz, cholesky, solve_triangular
from scipy.stats import norm, rayleigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
SOURCE = ROOT / 'literature_review/migration_data/xuan_2018/extraction'
OUT = HERE / 'diagnostics'
IDS = ['N010','N015','N016','N017','N018','N021','N024','N025']
SCENARIOS = {'methods_10min':1/6, '12h_endpoint':12/73, '16h_endpoint':16/73}
COLORS = ['#0072b2','#d55e00']
N = 1024

def covariance(tau, sigma, cad, allowance):
    t = np.arange(74)*cad
    z = t/tau
    shape = z+np.expm1(-z)
    small = z < 1e-3
    shape[small] = z[small]**2/2-z[small]**3/6+z[small]**4/24-z[small]**5/120
    v = sigma*sigma*tau**3*shape
    c = toeplitz(np.r_[v[1], .5*(v[2:]+v[:-2]-2*v[1:-1])])
    c += allowance**2*toeplitz(np.r_[2., -1., np.zeros(71)])
    d=cad/tau
    closed=np.r_[sigma*sigma*tau**3*(d+np.expm1(-d)),
                 sigma*sigma*tau**3/2*np.expm1(-d)**2*np.exp(-np.arange(72)*d)]
    closed+=allowance**2*np.r_[2.,-1.,np.zeros(71)]
    assert np.max(np.abs(c[0]-closed)) < max(c[0,0],1e-30)*1e-9
    return c

def acf(inc, maxlag=12):
    # Raw (zero-model-mean) dot-products, avoiding short-track demean bias.
    var = np.mean(np.sum(inc*inc,axis=-1),axis=-1)
    return np.stack([np.mean(np.sum(inc[:,:-k]*inc[:,k:],axis=-1),axis=-1)/var for k in range(1,maxlag+1)],axis=-1)

def longest_run(mask):
    best = np.zeros(mask.shape[0],dtype=int); now = best.copy()
    for col in mask.T:
        now = (now+1)*col; best = np.maximum(best,now)
    return best

def statistics(inc, root):
    speed = np.linalg.norm(inc,axis=-1)
    standardized = speed/root[0,0]
    whitened = solve_triangular(root,inc.transpose(1,0,2).reshape(73,-1),lower=True).reshape(73,-1,2).transpose(1,0,2)
    return dict(step=speed, standardized_step=standardized, acf=acf(inc), residual=whitened,
                slow=longest_run(standardized < rayleigh.ppf(.25)),
                fast=longest_run(standardized > rayleigh.ppf(.75)))

def write_csv(name, rows):
    pd.DataFrame(rows).to_csv(OUT/name,index=False)

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);a=p.parse_args()
    a.fit=a.fit.resolve()
    fit=json.loads(a.fit.read_text())
    if 'primary' in fit: fit=fit['primary']
    pars=np.exp(np.asarray(fit['z'][:20]).reshape(5,2,2))[0]
    tracks=pd.read_csv(SOURCE/'tracks.csv');qc=pd.read_csv(SOURCE/'track_qc.csv')
    observed=[]
    for cell in IDS:
        t=tracks[tracks.cell_id==cell].sort_values('frame_index')
        assert len(t)==74 and np.array_equal(t.frame_index.values,np.arange(74))
        assert set(t.qc_status)=={'reviewed_full_window'}
        observed.append(t[['x_um','y_um']].values)
    observed=np.array(observed); increments=np.diff(observed,axis=1)
    OUT.mkdir(exist_ok=True)
    summary=[]; corr=[]; cdf=[]; examples=[]; cache={}; validation=[]
    for scenario,cad in SCENARIOS.items():
        for px in (0,3):
            allowance=px*25/78
            for group in (0,1):
                tau,sigma=pars[group];cov=covariance(tau,sigma,cad,allowance)
                root=cholesky(cov,lower=True)
                rng=np.random.default_rng(202610030+group)
                sim=np.einsum('ij,njk->nik',root,rng.normal(size=(N,73,2)))
                ids=[0] if group==1 else list(range(1,8))
                obs=increments[ids]
                sm=statistics(sim,root);om=statistics(obs,root)
                cache[scenario,px,group]=(sim,sm,om,ids,cov)
                # Finite ensemble verification, independent of biological comparison.
                expected=2*cov[0,0]; realized=np.mean(sm['step']**2)
                se=np.std(np.mean(sm['step']**2,axis=1),ddof=1)/np.sqrt(N)
                assert abs(realized-expected)<6*se
                validation.append(dict(scenario=scenario,allowance_px=px,group=group,expected_step2=expected,simulated_step2=float(realized),mc_se=float(se)))
                for j,idx in enumerate(ids):
                    residual=om['residual'][j]
                    summary.append(dict(cell_id=IDS[idx],comparison_group='WGD+ model (appearance proxy)' if group else 'WGD- model (appearance proxy)',timing_scenario=scenario,cadence_h=cad,centroid_allowance_px=px,tau_h=tau,sigma=sigma,
                        mean_step_um=float(om['step'][j].mean()),scale_sensitivity_low_um=float(om['step'][j].mean()*.982),scale_sensitivity_high_um=float(om['step'][j].mean()*1.018),conditional_speed_um_h=float(om['step'][j].mean()/cad),model_mean_step_um=float(np.sqrt(np.pi*cov[0,0]/2)),
                        residual_rms=float(np.sqrt(np.mean(residual**2))),residual_mean_x=float(residual[:,0].mean()),residual_mean_y=float(residual[:,1].mean()),
                        residual_lag1=float(acf(residual[None],1)[0,0]),longest_slow_frames=int(om['slow'][j]),longest_fast_frames=int(om['fast'][j]),
                        simulated_slow_q025=float(np.quantile(sm['slow'],.025)),simulated_slow_q975=float(np.quantile(sm['slow'],.975)),simulated_fast_q025=float(np.quantile(sm['fast'],.025)),simulated_fast_q975=float(np.quantile(sm['fast'],.975))))
                # Match seven-track mean envelope to smaller-nucleus diagnostic.
                simacf=sm['acf'] if group else sm['acf'][:1022].reshape(146,7,12).mean(axis=1)
                for k in range(12):
                    corr.append(dict(scenario=scenario,allowance_px=px,group=group,lag_frames=k+1,observed=float(om['acf'][:,k].mean()),model=float(cov[0,k+1]/cov[0,0]),simulated_q025=float(np.quantile(simacf[:,k],.025)),simulated_q975=float(np.quantile(simacf[:,k],.975))))
                for q in np.linspace(.01,.99,99):
                    cdf.append(dict(scenario=scenario,allowance_px=px,group=group,quantile=q,observed_step_um=float(np.quantile(om['step'],q)),model_step_um=float(rayleigh.ppf(q,scale=np.sqrt(cov[0,0]))),observed_whitened=float(np.quantile(om['residual'],q)),model_whitened=float(norm.ppf(q))))
                if scenario=='methods_10min' and px==3:
                    pos=np.concatenate([np.zeros((8,1,2)),np.cumsum(sim[:8],axis=1)],axis=1)
                    for track,xy in enumerate(pos):
                        for frame,point in enumerate(xy):examples.append(dict(group=group,simulation_index=track,frame=frame,x_um=point[0],y_um=point[1]))
    write_csv('movie_track_checks.csv',summary);write_csv('movie_correlations.csv',corr);write_csv('movie_quantiles.csv',cdf);write_csv('movie_simulated_examples.csv',examples)
    # Measured trajectories and the first unselected model realizations on shared axes.
    fig,axes=plt.subplots(2,2,figsize=(10,9),layout='constrained')
    limit=0
    for group in (0,1):
        sim,sm,om,ids,cov=cache['methods_10min',3,group]
        real=observed[ids]-observed[ids,:1]
        fake=np.concatenate([np.zeros((len(ids),1,2)),np.cumsum(sim[:len(ids)],axis=1)],axis=1)
        for row,paths in enumerate((real,fake)):
            for j,path in enumerate(paths):
                axes[row,group].plot(path[:,0],path[:,1],lw=1.2,label=IDS[ids[j]] if row==0 else f'Draw {j+1}')
                axes[row,group].scatter(*path[-1],s=18)
            limit=max(limit,np.abs(paths).max());axes[row,group].scatter(0,0,c='black',marker='x')
            axes[row,group].set_title(('Measured: ' if row==0 else 'OU simulated: ')+('giant nucleus' if group else '7 smaller nuclei'))
            axes[row,group].set(xlabel='x displacement (µm)',ylabel='y displacement (µm)',aspect='equal');axes[row,group].legend(fontsize=8)
    for ax in axes.flat:ax.set(xlim=(-limit*1.05,limit*1.05),ylim=(-limit*1.05,limit*1.05));ax.grid(alpha=.2)
    fig.suptitle('Representative movie tracks: 74 frames; common scale\nOU comparison conditional on 10-minute cadence; 3-pixel centroid sensitivity')
    for ext in ('png','pdf'):fig.savefig(OUT/f'xuan_movie_trackplots.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    labels=['7 smaller nuclei / WGD− model','1 giant nucleus / WGD+ model']
    for group in (0,1):
        sim,sm,om,ids,cov=cache['methods_10min',3,group];col=COLORS[group]
        vals=np.sort(om['step'].ravel());axes[0,0].plot(vals,np.arange(1,len(vals)+1)/len(vals),c=col,label=labels[group])
        x=np.linspace(0,max(vals.max(),rayleigh.ppf(.995,scale=np.sqrt(cov[0,0]))),200);axes[0,0].plot(x,rayleigh.cdf(x,scale=np.sqrt(cov[0,0])),c=col,ls='--')
        rows=[r for r in corr if r['scenario']=='methods_10min' and r['allowance_px']==3 and r['group']==group]
        ks=np.arange(1,13);axes[0,1].plot(ks,[r['observed'] for r in rows],c=col);axes[0,1].plot(ks,[r['model'] for r in rows],c=col,ls='--');axes[0,1].fill_between(ks,[r['simulated_q025'] for r in rows],[r['simulated_q975'] for r in rows],color=col,alpha=.15)
        q=np.linspace(.01,.99,99);axes[0,2].plot(norm.ppf(q),np.quantile(om['residual'],q),c=col)
        for j,idx in enumerate(ids):axes[1,0].plot(np.arange(73),om['step'][j],lw=.8,c=col,alpha=.7)
        for k,key in enumerate(('slow','fast')):
            ax=axes[1,1+k];positions=np.arange(len(ids))+(0 if group else 1)
            ax.scatter(positions,om[key],c=col)
            lo,hi=np.quantile(sm[key],[.025,.975]);ax.axhspan(lo,hi,color=col,alpha=.12);ax.axhline(np.median(sm[key]),c=col,ls='--')
    axes[0,0].set(title='One-frame displacement CDF',xlabel='Step (µm)',ylabel='Cumulative fraction');axes[0,0].legend(fontsize=8)
    axes[0,1].set(title='Step-vector temporal correlation',xlabel='Lag (frames)',ylabel='Raw dot-product correlation')
    axes[0,2].plot([-3,3],[-3,3],c='gray',ls='--');axes[0,2].set(title='Forward OU prediction residuals',xlabel='Standard-normal quantile',ylabel='Whitened residual quantile')
    axes[1,0].set(title='Slow and fast fluctuations',xlabel='Frame interval',ylabel='Step (µm)')
    for k,label in enumerate(('Slow (<model Q25)','Fast (>model Q75)')):
        axes[1,1+k].set(title=f'Longest {label} episode',ylabel='Consecutive frame intervals',xticks=range(8),xticklabels=IDS);axes[1,1+k].tick_params(axis='x',rotation=45)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Movie diagnostics, conditional 10-minute cadence / 3-pixel sensitivity\nSolid: measured; dashed: OU expectation; shading: simulation 2.5–97.5% range (not biological confidence intervals)',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(OUT/f'xuan_movie_diagnostics.{ext}',dpi=180)
    plt.close(fig)
    paths=[a.fit,Path(__file__),SOURCE/'tracks.csv',SOURCE/'track_qc.csv',SOURCE/'METHODS.md',SOURCE/'protocol.json']
    meta=dict(generated=datetime.now(timezone.utc).isoformat(),fit_file=str(a.fit.relative_to(ROOT)),provenance={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},
        frame_count=74,accepted_track_ids=IDS,accepted_tracks=8,total_segments=len(qc),excluded_status_counts=qc[~qc.cell_id.isin(IDS)].qc_status.value_counts().to_dict(),
        timing_scenarios=SCENARIOS,simulations_per_group=N,seed_by_group=[202610030,202610031],allowance_px=[0,3],validation=validation,
        limitations=['Ploidy is not measured: groups are appearance proxies, not experimental WGD labels.','78 segments are not independent cells; only eight visually reviewed equal-window tracks enter checks.','Partial, merged, boundary-censored, ambiguous and unreviewed segments excluded, without interpolation.','No drift subtraction; no fixed specimen landmarks support it.','3 px treated as per-coordinate Gaussian SD only for sensitivity, not measured uncertainty or confidence interval.','1.8% scale reading sensitivity adds proportional position uncertainty; it is not biological variation.','Three physical timing scenarios unresolved; no acquisition-time claim.','Movie diagnostics never used as population fitting targets; finite trajectory envelopes are conditional model ranges.'])
    (OUT/'movie_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(f'Wrote movie diagnostics to {OUT}')

if __name__=='__main__':main()
