"""Bounded background Xuan memory fit with both WGD+ times greater than WGD−."""
import os,time,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit
from campaign import OUT,observations,modes,displacement_variance,score,save_json
HERE=Path(__file__).resolve().parent
DEST=OUT/'xuan_memory_ordered'
CAP=12.;FLOOR=1e-6;STARTS=32;MAX_NFEV=500
NAMES=['tp_reference','tm_reference','tp_gap_logit','tm_gap_logit','log_r']

def decode(z):
    a=np.exp(z[:2]);b=a+(CAP-a)*expit(z[2:4]);r=float(np.exp(z[4]))
    return a,b,r

def predict(z,early=False):
    a,b,r=decode(z);pred=[]
    for tp,tm in (a,b):
        theta=np.array([tp,tm,r]);rates,w=modes(theta,'memory')
        cad=1/6;lag=cad if early else 11.1111111111
        pred.append(dict(speed=float(np.sqrt(np.pi/2*displacement_variance(cad,rates,w))/cad),
                         msd=float(2*displacement_variance(lag,rates,w))))
    return pred

def result(z,fit,starts,tick,complete=False):
    obs=observations('Xuan2018');res,amp,rows,feasible=score(predict(z),obs,'Xuan2018')
    a,b,r=decode(z)
    checks=score(predict(z,True),observations('Xuan2018',True),'Xuan2018',amp)[2]
    assert np.all(a>0) and np.all(b>a) and np.all(b<CAP)
    return dict(study='Xuan2018',family='memory',complete=complete,generated=datetime.now(timezone.utc).isoformat(),
        parameters=dict(tp=float(a[0]),tm=float(a[1]),tp_poly=float(b[0]),tm_poly=float(b[1]),r=r),
        ratios=dict(tp=float(b[0]/a[0]),tm=float(b[1]/a[1])),shared_amplitude=amp,score=float(res@res),rows=rows,checks=checks,
        all_reading_intervals=all(o['within'] for o in rows),amplitude_interval_feasible=feasible,
        constraint_audit=dict(passes=True,time_lower_open=0,time_upper_open=12,both_poly_times_strictly_greater=True),
        numerical_bounds=dict(reference_time_floor_h=FLOOR,reference_time_ceiling_h=CAP-1e-8,gap_logit=[-16.,16.],feedback_r=[.005,.99]),
        numerical_floor_note='Positive open interval implemented with a 1e-6 h numerical floor, not the previous 0.05 h floor. Strict ordering uses a logistic fraction of the remaining interval up to 12 h.',
        starts_completed=len(starts),starts_planned=STARTS,max_nfev_per_start=MAX_NFEV,starts=starts,
        optimizer_success=bool(fit.success),z=z.tolist(),local_jacobian_singular_values=np.linalg.svd(fit.jac,compute_uv=False).tolist(),
        wall_seconds=time.monotonic()-tick,slurm_job=os.environ.get('SLURM_JOB_ID'),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        note='Analytical speed and MSD objective; shared r and amplitude. No Monte Carlo approximation or polarity-model inputs. Descriptive bounded search, not a proof of global impossibility.')

def publish(record):
    previous=json.loads((OUT/'memory_persistence_Xuan2018.json').read_text())
    parts=['<!doctype html><meta charset="utf-8"><title>Xuan memory: ordered persistence times</title><style>body{font:16px system-ui;max-width:1100px;margin:30px auto}td,th{padding:9px;border:1px solid #ccd5df}table{border-collapse:collapse}</style><h1>Xuan memory: both WGD+ times larger</h1>',
        f"<p>Generated {record['generated']}. Completed {record['starts_completed']}/{STARTS} starts. Campaign complete: {record['complete']}. Refresh manually for new checkpoints.</p>",
        '<p>All times are hours. Constraints: 0 &lt; tp− &lt; tp+ &lt; 12 and 0 &lt; tm− &lt; tm+ &lt; 12. Feedback and physical readout/noise scale are shared. This is a separate fit; the original fit and example tracks are unchanged.</p>',
        f"<p>Current best score: {record['score']:.6g}. Original unconstrained-order score: {previous['score']:.6g}. All fitted reading intervals met: {record['all_reading_intervals']}. These are descriptive losses, not likelihoods.</p>",
        '<table><tr><th>Parameter</th><th>Original</th><th>Ordered fit</th></tr>']
    for k,v in record['parameters'].items():parts.append(f"<tr><td>{k}</td><td>{previous['parameters'][k]:.7g}</td><td>{v:.7g}</td></tr>")
    parts.append(f"<tr><td>shared amplitude</td><td>{previous['shared_amplitude']:.7g}</td><td>{record['shared_amplitude']:.7g}</td></tr></table>")
    parts.append('<h2>Fitted measurements</h2><table><tr><th>Group</th><th>Metric</th><th>Observed</th><th>Interval</th><th>Predicted</th></tr>')
    for o in record['rows']:parts.append(f"<tr><td>{'WGD+' if o['group'] else 'WGD−'}</td><td>{o['metric']}</td><td>{o['obs']:g}</td><td>{o['low']:g}–{o['high']:g}</td><td>{o['predicted']:.6g}</td></tr>")
    parts.append('</table><p><a href="best.json">Current best parameters, starts and excluded checks</a> · <a href="../xuan_memory_tracks.html">Original fitted example tracks</a> · <a href="../report.html">Main report</a></p>')
    parts.append('<p>The open interval uses a 1e-6 h numerical floor. Gaps use a bounded logistic parameterization, so almost-equal times remain possible while strict inequalities hold. The feedback bounds remain 0.005–0.99. Small gaps or boundary optima should not be interpreted as identified effects.</p>')
    temp=DEST/'report.html.tmp';temp.write_text('\n'.join(parts));temp.replace(DEST/'report.html')

def main():
    DEST.mkdir(exist_ok=True);tick=time.monotonic();rng=np.random.default_rng(202609272)
    lo=np.array([np.log(FLOOR),np.log(FLOOR),-16.,-16.,np.log(.005)])
    hi=np.array([np.log(CAP-1e-8),np.log(CAP-1e-8),16.,16.,np.log(.99)])
    starts=[];best=None;best_fit=None;obs=observations('Xuan2018')
    def objective(z):return score(predict(z),obs,'Xuan2018')[0]
    for i in range(STARTS):
        if i==0:z=np.array([np.log(.05),np.log(.4),-5.,-1.,np.log(.897)])
        elif i<16:z=np.array([rng.uniform(np.log(.002),np.log(8.)),rng.uniform(np.log(.002),np.log(8.)),rng.uniform(-8,4),rng.uniform(-8,4),rng.uniform(np.log(.02),np.log(.98))])
        else:z=rng.uniform(lo,hi)
        fit=least_squares(objective,z,bounds=(lo,hi),max_nfev=MAX_NFEV,ftol=1e-10,xtol=1e-10,gtol=1e-10,diff_step=1e-4)
        starts.append(dict(index=i,z=fit.x.tolist(),score=float(2*fit.cost),success=bool(fit.success),nfev=fit.nfev))
        if best_fit is None or fit.cost<best_fit.cost:best_fit=fit;best=fit.x.copy()
        record=result(best,best_fit,starts,tick,complete=i==STARTS-1)
        save_json(DEST/'best.json',record);publish(record)
        print(i,2*fit.cost,fit.success,'best',record['score'],flush=True)
    print('Complete',flush=True)
if __name__=='__main__':main()
