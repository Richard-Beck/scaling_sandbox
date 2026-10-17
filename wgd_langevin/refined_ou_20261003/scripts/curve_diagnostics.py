"""Plot population MSD shape and residual checks using the selected frozen fit."""
from pathlib import Path
import csv, hashlib, json
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from figure_panels import load_fit_scenarios

RUN=Path(__file__).resolve().parents[1]
COLORS=['#168b87','#d87529']
def main():
    selected=load_fit_scenarios()['primary']; record=selected['record']
    fig,axs=plt.subplots(2,2,figsize=(11,7),sharex='col',height_ratios=[1.5,1],layout='constrained')
    exported=[]
    for g,color in enumerate(COLORS):
        rows=sorted([r for r in record['rows'] if r['metric']=='msd_curve' and r['group']==g],key=lambda r:r['lag_h'])
        t=np.array([r['lag_h'] for r in rows]); obs=np.array([r['obs'] for r in rows]); pred=np.array([r['predicted'] for r in rows])
        weight=np.array([r['objective_weight'] for r in rows]); residual=np.array([r['unweighted_residual'] for r in rows])
        axs[0,g].fill_between(t,[r['low'] for r in rows],[r['high'] for r in rows],color=color,alpha=.16,label='Raster extraction bounds')
        axs[0,g].plot(t,obs,color=color,lw=1.2,label='Published population MSD')
        axs[0,g].plot(t,pred,color='#19394d',lw=1.5,ls='--',label='Selected OU fit')
        axs[0,g].set(xscale='log',yscale='log',ylabel='MSD (µm²)',title=f'WGD{"−" if g==0 else "+"}: full population curve')
        axs[0,g].legend(frameon=False,fontsize=8)
        axs[1,g].plot(t,np.log(pred/obs),color=color,lw=1.5)
        axs[1,g].axhline(0,color='#19394d',lw=.8)
        axs[1,g].set(xscale='log',xlabel='Lag (h; published population clock)',ylabel='Log(predicted / observed)')
        edges=np.exp(np.linspace(np.log(t.min()),np.log(t.max()),4))
        for ax in axs[:,g]:
            for x in edges[1:-1]:ax.axvline(x,color='#6d7c86',ls=':',lw=.8)
            ax.grid(alpha=.15)
        contrib=[float(np.sum(weight[np.array([r['lag_bin'] for r in rows])==b]*residual[np.array([r['lag_bin'] for r in rows])==b]**2)) for b in range(3)]
        axs[1,g].text(.02,.04,'Early / middle / late loss: '+ ' / '.join(f'{c:.3f}' for c in contrib),transform=axs[1,g].transAxes,fontsize=8)
        exported.extend([dict(group=g,lag_h=r['lag_h'],lag_bin=r['lag_bin'],observed=r['obs'],predicted=r['predicted'],log_prediction_over_observed=float(np.log(r['predicted']/r['obs'])),objective_weight=r['objective_weight'],weighted_squared_residual=r['residual']**2) for r in rows])
    fig.suptitle('Xuan population MSD: shape residuals under the balanced OU fit',fontsize=14,color='#19394d')
    for suffix in ('png','pdf','svg'):fig.savefig(RUN/'figures'/f'xuan_population_MSD_diagnostics.{suffix}',dpi=220)
    plt.close(fig)
    target=RUN/'figures/xuan_population_MSD_residuals.csv'
    with target.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(exported[0]));w.writeheader();w.writerows(exported)
    caption='The full published population MSD curves and their fitted OU predictions are shown on logarithmic lag axes. Shading shows raster extraction bounds, not population uncertainty. Residuals use the log predicted/observed ratio; vertical guides split equal log-lag regimes. Native digitized points share one third of each group curve block within each regime. Each group curve has total objective weight one, matching one sampled-speed target. Early/middle/late displayed losses are weighted squared standardized residual sums. These systematic shape residuals test fit adequacy; they are not independent per-pixel sampling errors.'
    (RUN/'figures/xuan_population_MSD_diagnostics_caption.md').write_text(caption+'\n')
    inputs=[selected['path'],RUN/'data/xuan_msd_curve.csv',Path(__file__)]
    meta=dict(created_utc=datetime.now(timezone.utc).isoformat(),selected_fit=selected['path'].name,source_sha256={str(p.relative_to(RUN)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},rows=len(exported),caption=caption)
    (RUN/'figures/xuan_population_MSD_diagnostics_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':main()
