"""Evaluate the established WGD−-fixed persistence/noise surface; no refitting."""
from pathlib import Path
import json, csv, hashlib
from datetime import datetime, timezone
import numpy as np
import fit
from figure_panels import load_fit_scenarios
RUN=Path(__file__).resolve().parents[1]
def main():
    selected=load_fit_scenarios()['primary'];base=np.asarray(selected['z']);pars=base[:20].reshape(5,2,2)
    fold=np.exp((pars[:,1,:]-pars[:,0,:]).mean(axis=0));plus=pars[:,1,:]
    limits=[]
    for k,(lo,hi) in enumerate(((.5,30.),(.05,2.))):
        model_lo=fit.LOW[:20].reshape(5,2,2)[:,1,k];model_hi=fit.HIGH[:20].reshape(5,2,2)[:,1,k]
        limits.append((max(np.log(lo/fold[k]),float(np.max(model_lo-plus[:,k]))),min(np.log(hi/fold[k]),float(np.min(model_hi-plus[:,k])))))
    axes=[np.unique(np.r_[np.geomspace(fold[k]*np.exp(limits[k][0]),fold[k]*np.exp(limits[k][1]),61),1.,fold[k]]) for k in range(2)]
    rows=[]
    for tau in axes[0]:
        for sigma in axes[1]:
            z=base.copy();p=pars.copy();p[:,1,:]+=np.log(np.array([tau,sigma])/fold);z[:20]=p.ravel()
            r,rr=fit.data_residuals(fit.analytic(z));mask=np.asarray([q['study']!='Hosny_baseline' or q.get('quantile')==1 for q in rr]);data=float(r[mask]@r[mask]);pool=fit.pooling_residuals(z,fit.SCENARIOS['primary']);prior=fit.persistence_residuals(z,'primary');total=data+pool@pool+prior@prior
            rows.append(dict(tau_ratio=float(tau),sigma_ratio=float(sigma),fixed_drift_delta_data=data-selected['record']['data_score'],fixed_drift_delta_total=total-selected['record']['total_score'],status='valid'))
    path=RUN/'results/main_relative_surface.csv'
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    meta=dict(created_utc=datetime.now(timezone.utc).isoformat(),selected_fit=selected['path'].name,selected_fit_sha256=hashlib.sha256(selected['path'].read_bytes()).hexdigest(),folds=fold.tolist(),rows=len(rows),definition='All WGD− values, Hosny drifts and study deviations remain fixed; only a common WGD+ log shift is swept. This is a heuristic loss sensitivity, not a likelihood.')
    path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':main()
