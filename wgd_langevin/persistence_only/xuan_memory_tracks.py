"""Stationary example trajectories from the saved Xuan linear-memory fit."""
import csv,json,hashlib
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky,toeplitz
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from campaign import OUT,decode,modes,displacement_variance,save_json

FIT=OUT/'memory_persistence_Xuan2018.json'
COLORS=['#0072b2','#d55e00'];LABELS=['WGD−','WGD+']
DURATION=12.;FINE=1/60;CADENCE=1/6;N=32;SEEDS=[202609270,202609271]

def covariance(cadence,duration,rates,weights):
    variance=displacement_variance(np.arange(round(duration/cadence)+1)*cadence,rates,weights)
    first=np.r_[variance[1],.5*(variance[2:]+variance[:-2]-2*variance[1:-1])]
    return toeplitz(first)

def sample(cov,n,seed,amplitude):
    root=cholesky(cov,lower=True)
    rng=np.random.default_rng(seed)
    increments=(root@rng.normal(size=(len(root),n*2))).reshape(len(root),n,2).transpose(1,0,2)*amplitude
    return np.concatenate([np.zeros((n,1,2)),np.cumsum(increments,axis=1)],axis=1)

def main():
    fit=json.loads(FIT.read_text());amp=fit['shared_amplitude'];records=[];tracks=[]
    for group in (0,1):
        theta,_,_=decode(fit['z'],fit['names'],group);rates,w=modes(theta,'memory')
        finecov=covariance(FINE,DURATION,rates,w);coarsecov=covariance(CADENCE,DURATION,rates,w)
        # Exact summation of fine increments must recover the fitting cadence covariance.
        aggregated=finecov.reshape(72,10,72,10).sum(axis=(1,3))
        err=float(np.max(np.abs(aggregated-coarsecov))/np.max(np.abs(coarsecov)))
        assert err<1e-8
        positions=sample(finecov,N,SEEDS[group],amp);tracks.append(positions)
        validation=sample(coarsecov,8192,SEEDS[group]+1000,amp)
        speed=np.linalg.norm(np.diff(validation,axis=1),axis=2).sum(axis=1)/DURATION
        endpoint2=np.sum(validation[:,-1]**2,axis=1)
        speed_exact=float(amp*np.sqrt(np.pi/2*displacement_variance(CADENCE,rates,w))/CADENCE)
        msd_exact=float(2*amp**2*displacement_variance(DURATION,rates,w))
        checks={}
        for name,values,target in [('speed',speed,speed_exact),('msd_12h',endpoint2,msd_exact)]:
            se=float(values.std(ddof=1)/np.sqrt(len(values)));difference=float(values.mean()-target)
            assert abs(difference)<5*se
            checks[name]=dict(simulated=float(values.mean()),analytic=target,monte_carlo_se=se,z=difference/se)
        record=dict(group=LABELS[group],tp_h=float(theta[0]),tm_h=float(theta[1]),feedback_r=float(theta[2]),
                    rates_per_h=rates.tolist(),decay_times_h=(1/rates).tolist(),normalized_weights=(w/w.sum()).tolist(),
                    correlation_1h=float(np.sum(w*np.exp(-rates))/w.sum()),integrated_correlation_h=float(np.sum(w/rates)/w.sum()),
                    speed_um_h=speed_exact,msd_12h_um2=msd_exact,validation=checks,covariance_aggregation_relative_error=err)
        records.append(record)
    times=np.arange(721)*FINE
    with (OUT/'xuan_memory_example_tracks.csv').open('w') as handle:
        writer=csv.writer(handle);writer.writerow(['condition','track','time_h','x_um','y_um','camera_frame'])
        for group,positions in enumerate(tracks):
            for i,position in enumerate(positions):
                for j,(t,xy) in enumerate(zip(times,position)):
                    writer.writerow([LABELS[group],i+1,float(t),float(xy[0]),float(xy[1]),j%10==0])
    # All example tracks, at the original imaging cadence, with identical axes.
    fig,axes=plt.subplots(1,2,figsize=(11,5.8),layout='constrained')
    limit=np.ceil(max(np.abs(p[:,::10]).max() for p in tracks)/10)*10
    for group,ax in enumerate(axes):
        for p in tracks[group][:,::10]:
            ax.plot(p[:,0],p[:,1],color=COLORS[group],alpha=.42,lw=.9)
            ax.scatter(*p[-1],color=COLORS[group],s=14,alpha=.75)
        ax.scatter(0,0,c='black',marker='x',s=55,zorder=5)
        ax.set(xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)',aspect='equal')
        r=records[group]
        ax.set_title(f"{LABELS[group]}: tp = {60*r['tp_h']:.2f} min; tm = {r['tm_h']:.3g} h")
        ax.grid(alpha=.15)
    fig.suptitle('Xuan fitted memory model: 32 independent example tracks per condition\n12 hours, sampled every 10 minutes; common origin and identical axes',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(OUT/f'xuan_memory_tracks.{ext}',dpi=180)
    plt.close(fig)
    # First six draws, without selection by displacement or visual appearance.
    fig,axes=plt.subplots(2,6,figsize=(15,6),layout='constrained')
    limit=np.ceil(max(np.abs(p[:6]).max() for p in tracks)/10)*10
    norm=Normalize(0,DURATION)
    for group in (0,1):
        for i in range(6):
            ax=axes[group,i];p=tracks[group][i]
            segments=np.stack([p[:-1],p[1:]],axis=1)
            line=LineCollection(segments,cmap='viridis',norm=norm,linewidths=1.2);line.set_array(times[:-1]);ax.add_collection(line)
            ax.scatter(p[::10,0],p[::10,1],c=times[::10],cmap='viridis',norm=norm,s=5)
            ax.scatter(0,0,c='black',marker='x',s=24);ax.scatter(*p[-1],c='red',s=15)
            ax.set(xlim=(-limit,limit),ylim=(-limit,limit),aspect='equal',title=f'{LABELS[group]} · track {i+1}')
            if group==1:ax.set_xlabel('x (µm)')
            if i==0:ax.set_ylabel('y (µm)')
            ax.grid(alpha=.15)
    fig.colorbar(line,ax=axes,shrink=.75,label='Time (hours)')
    fig.suptitle('First six simulated tracks per condition, without selection\nLines: 1-minute resolution; dots: 10-minute camera frames; black ×: start; red: endpoint',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(OUT/f'xuan_memory_individual_tracks.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
    for group,r in enumerate(records):
        rates=np.array(r['rates_per_h']);weights=np.array(r['normalized_weights']);w=weights*r['tp_h']/2
        for ax,ts in zip(axes[:2],[np.linspace(0,.5,500),np.linspace(0,12,500)]):
            acf=np.sum(weights[:,None]*np.exp(-rates[:,None]*ts),axis=0)
            ax.plot(ts,acf,color=COLORS[group],label=LABELS[group],lw=2)
        ts=np.linspace(0,12,500);axes[2].plot(ts,2*amp**2*displacement_variance(ts,rates,w),color=COLORS[group],label=LABELS[group],lw=2)
        target=next(o for o in fit['rows'] if o['group']==group and o['metric']=='msd')
        axes[2].errorbar(target['lag_h'],target['obs'],yerr=[[target['obs']-target['low']],[target['high']-target['obs']]],fmt='s',color=COLORS[group],capsize=5)
    axes[0].set(title='Rapid velocity decorrelation',xlabel='Lag (hours)',ylabel='Normalized velocity autocorrelation',ylim=(0,1))
    axes[1].set(title='Hour-scale correlation tail (zoom)',xlabel='Lag (hours)',ylabel='Normalized velocity autocorrelation',ylim=(0,.12))
    axes[2].set(title='Expected displacement grows differently',xlabel='Lag (hours)',ylabel='MSD (µm²)')
    for ax in axes:ax.grid(alpha=.2);ax.legend(frameon=False)
    fig.suptitle('Same fitted model: shorter tp and tm, but a stronger slow-component weight in WGD+',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(OUT/f'xuan_memory_dynamics.{ext}',dpi=180)
    plt.close(fig)
    def msd_difference(t):
        vals=[]
        for r in records:
            rates=np.array(r['rates_per_h']);w=np.array(r['normalized_weights'])*r['tp_h']/2
            vals.append(float(displacement_variance(t,rates,w)))
        return vals[1]-vals[0]
    crossing=brentq(msd_difference,.2,12)
    generated=datetime.now(timezone.utc).isoformat()
    save_json(OUT/'xuan_memory_track_summary.json',dict(generated=generated,fit_file=FIT.name,fit_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        examples_per_condition=N,duration_h=DURATION,camera_cadence_h=CADENCE,display_resolution_h=FINE,seeds=SEEDS,
        sampling='Exact stationary Gaussian integrated-velocity covariance; independent seeds for the two conditions; no selected tracks or endpoint alignment.',
        shared_amplitude=amp,groups=records,msd_crossing_h=crossing,validation_passes=True))
    headers=['Condition','tp (h)','tm (h)','Fast decay (h)','Slow decay (h)','Slow weight','Correlation at 1 h','Mean sampled speed (µm/h)']
    rows=[]
    for r in records:
        rows.append([r['group'],f"{r['tp_h']:.5g}",f"{r['tm_h']:.5g}",f"{r['decay_times_h'][1]:.4g}",f"{r['decay_times_h'][0]:.4g}",f"{100*r['normalized_weights'][0]:.3g}%",f"{r['correlation_1h']:.4g}",f"{r['speed_um_h']:.3g}"])
    table='<table><tr>'+''.join('<th>'+h+'</th>' for h in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+v+'</td>' for v in row)+'</tr>' for row in rows)+'</table>'
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Xuan memory model: example tracks</title><style>body{font:16px system-ui;max-width:1450px;margin:30px auto;padding:0 20px;color:#172735}img{width:100%}table{border-collapse:collapse}td,th{padding:10px;border:1px solid #ccd5df}th{background:#eef3f7}a{color:#075ca4}p{max-width:1050px}</style><h1>Xuan: tracks from the fitted memory model</h1>'''
    document+=f'<p>Generated {generated}. These are simulated examples, not experimental tracks. No refitting was performed.</p>'
    document+='<p><a href="report.html">Main fit report</a> · <a href="xuan_memory_example_tracks.csv">Track coordinates CSV</a> · <a href="xuan_memory_track_summary.json">Parameters and numerical checks</a> · <a href="memory_persistence_Xuan2018.json">Original fitted parameters</a></p>'+table
    document+=f'''<p>Both time constants decrease: tp ratio = {fit['ratios']['tp']:.3g}, tm ratio = {fit['ratios']['tm']:.3g}. tp describes the instantaneous relaxation of p; tm describes the relaxation of the latent memory state. Neither alone is the decay time of the coupled velocity process. Its autocorrelation is a weighted sum of two exponentials. WGD+ has a shorter slow decay time, but approximately 16 times as much weight in that component. This gives stronger hour-scale correlations while lowering sampled speed. The predicted MSD curves cross at approximately {crossing:.2f} hours.</p>
<p>The shared feedback is r = {records[0]['feedback_r']:.5g}; the shared readout scale is {amp:.5g}. The stationary covariance closure is exactly the one used for fitting. Positions start at a common origin with stationary velocities; there is no imposed directional cue. Seeds are independent between conditions, so equal track numbers do not represent paired cells.</p>'''
    for filename,title in [('xuan_memory_tracks','All 32 example tracks'),('xuan_memory_individual_tracks','Six unselected individual examples per condition'),('xuan_memory_dynamics','Correlation and displacement explain the tracks')]:
        document+=f'<h2>{title}</h2><a href="{filename}.png"><img src="{filename}.png" alt="{title}"></a><p><a href="{filename}.pdf">PDF</a> · <a href="{filename}.png">PNG</a></p>'
    document+='''<p>The longest decay component is approximately 117 h in WGD− but has only 0.62% weight; a 12 h experiment does not establish that extrapolated tail. The integral of the normalized correlation over infinite time is about 0.804 h for WGD− and 0.484 h for WGD+. Thus WGD+ is not more persistent under every possible definition: its stronger correlation over the measured hour-scale window is what matters here.</p>
<p>These are one best-fit parameter set, not uniquely identified estimates. WGD+ tp touches the 0.05 h search floor and WGD− tm touches the 12 h ceiling. Those limits and the limited summary data qualify literal interpretation of the time constants. Track sampling was checked against analytical mean speed and 12 h MSD using 8,192 independent tracks per condition; the fine-to-camera covariance aggregation was also checked.</p></html>'''
    (OUT/'xuan_memory_tracks.html').write_text(document)
    print(json.dumps(dict(groups=records,msd_crossing_h=crossing,validation_passes=True),indent=2))

if __name__=='__main__':main()
