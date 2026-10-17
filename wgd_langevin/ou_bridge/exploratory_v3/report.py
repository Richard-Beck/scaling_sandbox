"""Small visual reports, updated from completed results without browser timers."""
import sys,json,csv,time,argparse
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from free_fit import OUT,BASE,LENGTHS,FITS
sys.path.insert(0,str(BASE/'refit_v2'))
import campaign as c
from reduced import NAMES,active
LABELS={'wave_pinning':'Wave pinning','otsuji':'Otsuji','goryachev':'Goryachev','legi':'LEGI','debelly':'DeBelly'}
FLABELS=['Gaussian','Memory','Coupled','Self-activation'];COLORS=['#0072B2','#E69F00','#009E73','#CC79A7']
STYLE='body{font:16px/1.5 system-ui;margin:28px auto;padding:0 24px;max-width:1500px;color:#172735}img{width:100%;height:auto}figure{margin:30px 0}a{color:#075ca4}summary{cursor:pointer;font-weight:bold}table{border-collapse:collapse}td,th{padding:6px 10px;border:1px solid #ddd}.stamp{color:#657381;font-size:13px}nav{display:flex;gap:20px;flex-wrap:wrap}'

def shell(title,body):return '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+title+'</title><style>'+STYLE+'</style><h1>'+title+'</h1><p class="stamp">Generated '+datetime.now(timezone.utc).isoformat()+'</p>'+body+'</html>'
def write(name,text):
    p=OUT/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(text);t.replace(p)

def spatial():
    fits=[json.loads(p.read_text())for p in (OUT/'fits').glob('*.json')]
    body=f'<p><b>{len(fits)}/180 independent fits complete.</b> Nine lengths per spatial source; all seven experiments used at every length. No held-out conditions and no imposed length dependence.</p><nav><a href="msd.html">Reference-to-WGD MSD predictions</a>'+''.join(f'<a href="#{m}">{LABELS[m]}</a>'for m in c.MODELS)+'</nav>'
    body+='<p>Each image shows the freely fitted parameter trends for one spatial source. Colors distinguish the four OU candidates. Large points and connecting lines show the lowest-error fit at each length; small translucent points show alternative starts within 5% of that objective. These alternatives are not confidence intervals. Lines connect fits, not fitted power laws. Crosses mark a best start that stopped without optimizer convergence. Gray dashed lines show parameter bounds.</p><p>tp and tm use each source’s numerical clock (100 s for the four Jilkine models; 1000 s for DeBelly). Compare trends within a source. The latent-state conventions and bounds are unchanged; “free” means independently fitted across lengths, not unbounded.</p>'
    body+='<details><summary>Fit quality: all experiments, no held-outs</summary><p>Post-event normalized RMSE, pooled by RMS across all experiments and fitted lengths. Smaller is better. Parameter trends from poor fits should not be interpreted as measured mechanisms.</p><table><tr><th>Spatial source</th>'+''.join('<th>'+f+'</th>'for f in FLABELS)+'</tr>'
    for m in c.MODELS:
        body+='<tr><td>'+LABELS[m]+'</td>'
        for family in c.FAMILIES:
            values=[metric['post_nrmse']for r in fits if r['model']==m and r['family']==family for metric in r['metrics']]
            body+='<td>'+(f'{np.sqrt(np.mean(np.square(values))):.3f}'if values else 'pending')+'</td>'
        body+='</tr>'
    body+='</table></details>'
    allrows=[]
    for m in c.MODELS:
        rs=[r for r in fits if r['model']==m]
        if not rs:continue
        fig,axes=plt.subplots(2,4,figsize=(17,8),layout='constrained')
        for i,ax in enumerate(axes.flat):
            for f,color,label in zip(c.FAMILIES,COLORS,FLABELS):
                if i not in active(f):continue
                subset=sorted([r for r in rs if r['family']==f],key=lambda r:r['length'])
                if not subset:continue
                ax.plot([r['length']for r in subset],[r['theta'][i]for r in subset],marker='o',color=color,label=label,ms=4,lw=1.2,ls=['-','--','-.',':'][c.FAMILIES.index(f)])
                for r in subset:
                    near=[s['theta'][i]for s in r['starts']if s['cost']<=1.05*r['cost']+1e-8]
                    ax.scatter([r['length']]*len(near),near,color=color,alpha=.25,s=12)
                    best=min(r['starts'],key=lambda s:s['cost'])
                    if not best['success']:ax.scatter(r['length'],r['theta'][i],marker='x',color='#b00020',s=45)
                    allrows.append(dict(model=m,family=f,length=r['length'],parameter=NAMES[i],value=r['theta'][i],alternative_min=min(near),alternative_max=max(near),cost=r['cost'],best_converged=best['success']))
            from reduced import LOW,HIGH
            ax.axhline(LOW[i],color='#aaa',ls='--',lw=.7);ax.axhline(HIGH[i],color='#aaa',ls='--',lw=.7)
            if i==2:ax.axhline(.995,color=COLORS[1],ls=':',lw=.8)
            if i in (2,6):ax.set_yscale('symlog',linthresh=.01);ax.set_ylim(bottom=0)
            else:ax.set_yscale('log')
            ax.set_title(NAMES[i]);ax.set_xlabel('Length (µm)');ax.grid(alpha=.15)
        fig.suptitle(LABELS[m]+f' · free parameter trends ({len(rs)}/36 fits)',fontsize=18)
        fig.legend(handles=[Line2D([0],[0],color=color,marker='o',ls=ls,label=label)for color,label,ls in zip(COLORS,FLABELS,['-','--','-.',':'])],loc='outside lower center',ncol=4,frameon=False)
        name=f'parameters_{m}.png';fig.savefig(OUT/name,dpi=150);plt.close(fig)
        body+=f'<section id="{m}"><h2>{LABELS[m]}</h2><a href="{name}"><img src="{name}" alt="Free fitted parameter trends across nine lengths for {LABELS[m]}"></a></section>'
    if allrows:
        with (OUT/'parameter_trends.csv').open('w',newline='')as h:
            w=csv.DictWriter(h,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    body+='<p><a href="parameter_trends.csv">All fitted parameters and alternative ranges (CSV)</a></p><details><summary>Fitting details</summary><p>Eight starts per length and candidate, including nearby old fits as initial guesses. The eight equally weighted observables are seven phase-balanced trajectories and the weak-plus minus weak-minus response. No parameters are shared between lengths. All sources retain their prior tested range, with four intermediate lengths added.</p></details>'
    write('spatial.html',shell('Free parameter trends across length',body));return len(fits)

def msd():
    results=[json.loads(p.read_text())for p in (OUT/'reference').glob('*.json')]
    body=f'<p><b>{len(results)}/20 reference calibrations complete.</b> For each study, calibrate the four OU families on non-WGD observations only; then apply the five spatial-source power laws to each family. This produces 20 candidate predictions per study. The baseline calibration is shared where the reference inference problem is identical.</p><nav><a href="spatial.html">Free spatial parameter trends</a>'+''.join(f'<a href="#{s}">{s}</a>'for s in c.STUDIES)+'</nav>'
    body+='<p><b>Black: fitted reference. Orange: predicted WGD at sqrt(2) times the reference length.</b> The selected power-law exponents remain fixed; each affected parameter is multiplied by sqrt(2) raised to its exponent. Every other parameter, the noise scale and the physical readout scale stay fixed. No WGD observation enters calibration or curve selection.</p><p>All 20 facets within a study share the same axes. Shaded ranges show alternative near-optimal reference fits found by multistart optimization, not confidence intervals. Speed-only studies cannot determine an MSD curve uniquely. Weakly supported spatial laws remain exploratory. This report uses the existing 20 laws; the new free-length analysis has not been converted into replacement laws.</p>'
    csvrows=[]
    for study in c.STUDIES:
        found={r['family']:r for r in results if r['study']==study}
        if not found:continue
        fig,axes=plt.subplots(5,4,figsize=(18,16),layout='constrained',sharex=True,sharey=True)
        for i,m in enumerate(c.MODELS):
            for j,f in enumerate(c.FAMILIES):
                ax=axes[i,j];r=found.get(f)
                if i==0:ax.set_title(FLABELS[j],fontsize=14)
                if j==0:ax.set_ylabel(LABELS[m]+'\nMSD (µm²)')
                ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.18);ax.set_xlabel('Lag (hours)')
                if not r:ax.text(.5,.5,'Reference fit pending',transform=ax.transAxes,ha='center');continue
                p=next(p for p in r['predictions']if p['model']==m);lags=np.array(r['lags'])
                for best,alts,color in [(r['reference_msd'],r['alternative_reference_msd'],'#222'),(p['wgd_msd'],p['alternative_wgd_msd'],'#d76b00')]:
                    if alts:
                        values=np.array([best]+alts);ax.fill_between(lags,values.min(0),values.max(0),color=color,alpha=.13,lw=0)
                    ax.plot(lags,best,color=color,lw=2)
                if study=='Xuan2018':
                    obs=c.OBS[(c.OBS.study==study)&(c.OBS.metric=='msd')&(c.OBS.role.isin(['fit','check']))]
                    for o in obs.itertuples():
                        color='#222'if o.group=='reference'else '#d76b00'
                        ax.errorbar(o.lag_h,o.value,yerr=[[o.value-o.low],[o.high-o.value]],fmt='o'if o.role=='fit'else 's',ms=4,color=color,capsize=2,markerfacecolor=color if o.role=='fit'else 'white')
                text='; '.join(k+' ×'+f'{v:.2f}'for k,v in p['multipliers'].items())
                if p['inactive']:text+='\nGain law has no effect without cue'
                if p['spatial_calibration_improvement']<.1:text+='\nSpatial law weakly supported'
                ax.text(.03,.97,text,transform=ax.transAxes,va='top',fontsize=8,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
                for lag,ref,wgd in zip(lags,r['reference_msd'],p['wgd_msd']):csvrows.append(dict(study=study,model=m,family=f,length_ratio=r['length_ratio'],lag_h=float(lag),reference_msd=ref,wgd_msd=wgd))
        fig.suptitle(study+' · reference fit → fixed length-change prediction',fontsize=18)
        fig.legend(handles=[Line2D([0],[0],color='#222',label='Reference fit'),Line2D([0],[0],color='#d76b00',label='WGD prediction: length ×sqrt(2)')],loc='outside lower center',ncol=2,frameon=False)
        name=f'msd_{study}.png';fig.savefig(OUT/name,dpi=150);plt.close(fig)
        body+=f'<section id="{study}"><h2>{study}</h2>'
        body+='<p>'+('Observed MSD points: filled circles are the long-lag summaries; open squares are early-lag checks. Only reference fitting observations were used. WGD points are comparisons, never fitting targets.'if study=='Xuan2018'else 'No numerical observed MSD curve is available in the transcribed data. These curves are predictions conditioned on the available reference summaries.')+'</p>'
        if study.startswith('Jemaa'):body+='<p><b>Speed-only calibration: persistence and MSD shape are underdetermined.</b></p>'
        unfinished=[f for f,r in found.items() if not r['polish_success']]
        if unfinished:body+='<p><b>Reference optimization reached its iteration limit:</b> '+', '.join(unfinished)+'. These curves are provisional.</p>'
        missed=[f for f,r in found.items() if any(row.get('within') is False for row in r['reference_rows']+r['reference_validation'])]
        if missed:body+='<p><b>Reference fit or its independent simulation misses a recorded reading interval:</b> '+', '.join(missed)+'. Interpret their WGD curves conditional on that imperfect baseline fit.</p>'
        body+=f'<a href="{name}"><img src="{name}" alt="Twenty reference and WGD MSD comparisons for {study}"></a><details><summary>What was fitted to the reference?</summary><table><tr><th>Candidate</th><th>Observable</th><th>Observed</th><th>Fitted</th><th>Independent simulation</th></tr>'
        for f in c.FAMILIES:
            if f not in found:continue
            r=found[f]
            for idx,row in enumerate(r['reference_rows']):
                observable=row['metric']+(' '+row['block']+' Q'+str(row['quantile'])if 'block'in row else '')
                body+='<tr>'+''.join('<td>'+str(x)+'</td>'for x in [f,observable,f"{row['obs']:.3g}",f"{row['predicted']:.3g}",f"{r['reference_validation'][idx]['predicted']:.3g}"])+'</tr>'
        body+='</table></details></section>'
    if csvrows:
        with (OUT/'msd_predictions.csv').open('w',newline='')as h:
            w=csv.DictWriter(h,fieldnames=list(csvrows[0]));w.writeheader();w.writerows(csvrows)
    body+='<p><a href="msd_predictions.csv">All MSD curves (CSV)</a></p><details><summary>Model and numerical details</summary><p>Gaussian and memory curves use exact stationary MSD formulas. Nonlinear curves use independent 2D trajectories and time-averaged squared displacement, with at least 96 hours of burn-in and step at most 0.01 hours. Nonlinear q/p closure and linear-memory noise covariance match the previous campaign. Six bounded reference starts plus polishing; physical amplitude profiled using reference observations only. Hosny uses non-PACC control quartiles and a constant directional cue, so its MSD includes drift. Other studies are unforced, making cue gain unobservable. The gray/orange spread reflects up to three alternative good reference fits; absence of spread does not establish identifiability.</p></details>'
    write('msd.html',shell('Reference-calibrated WGD predictions: MSD versus lag',body));return len(results)

def run(watch=False):
    previous=(-1,-1)
    while True:
        counts=(len(list((OUT/'fits').glob('*.json'))),len(list((OUT/'reference').glob('*.json'))))
        if counts[0]!=previous[0]:spatial()
        if counts[1]!=previous[1]:msd()
        previous=counts
        print(datetime.now(timezone.utc).isoformat(),counts,flush=True)
        if not watch or counts==(180,20):break
        time.sleep(30)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--watch',action='store_true');run(p.parse_args().watch)
