"""Per-experiment errors; independent fits are unavailable at untrained lengths."""
import sys,json,csv
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from design import OUT as OLD,MODELS,FAMILIES,ASSAYS,CALIBRATION,lengths,Protocol
from fit import load_cases,prepare,evaluate
from reduced import integrate
from scaling import parameters
OUT=Path(__file__).resolve().parent/'results'

def build():
    rows=[];viewer={};objectives=[]
    for m in MODELS:
        for f in FAMILIES:
            s=json.loads((OLD/'scaling'/f'{m}_{f}.json').read_text())
            laws={'constant':next(r for r in s['laws']if not r['indices']), 'power':min((r for r in s['laws']if r['indices']),key=lambda r:r['cost'])}
            for L in lengths(m):
                cases=load_cases(m,L);training=L in lengths(m,True)
                prepared,scale=prepare([c for c in cases if c['assay']in CALIBRATION])
                for regime in ['independent','constant','power']:
                    if regime=='independent':
                        if not training:continue
                        record=json.loads((OLD/'fits'/f'{m}_{L:g}_{f}.json').read_text());theta=record['theta']
                    else:theta=parameters(laws[regime],L)
                    metrics,pred=evaluate(theta,f,cases,scale)
                    for r in metrics:
                        rows.append(dict(model=m,family=f,length=L,length_role='training' if training else 'held_out',regime=regime,**r))
                    if training:
                        for assay,(t,y,spec,w) in zip(CALIBRATION+['weak_difference'],prepared):
                            val=integrate(theta,t,spec,FAMILIES.index(f))[:,0] if spec.ndim==1 else integrate(theta,t,spec[0],FAMILIES.index(f))[:,0]-integrate(theta,t,spec[1],FAMILIES.index(f))[:,0]
                            objectives.append(dict(model=m,family=f,length=L,regime=regime,assay=assay,error=float(np.sqrt(5*np.sum(((val-y)*w)**2)))))
                    for c in cases:
                        key='|'.join([m,str(L),c['assay'],regime]);protocol=Protocol(m,c['assay'])
                        focus=0 if c['assay']in ('pulse','competition')else protocol.condition-.1*(protocol.end-protocol.condition)
                        onset=np.searchsorted(c['time'],focus)
                        idx=np.unique(np.r_[np.linspace(0,onset,min(50,onset+1)).astype(int),np.linspace(onset,len(c['time'])-1,min(600,len(c['time'])-onset)).astype(int)])
                        if key not in viewer:viewer[key]=dict(time=c['time'][idx].tolist(),target=c['p'][idx].tolist(),families={},role='training length' if training else 'held-out length')
                        viewer[key]['families'][f]=pred[c['assay']][idx,0].tolist()
        print(m,flush=True)
    for name,table in [('spatial_errors.csv',rows),('training_objective_errors.csv',objectives)]:
        with (OUT/name).open('w',newline='') as h:
            w=csv.DictWriter(h,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    render(rows,viewer)

LABELS={'wave_pinning':'Wave pinning','otsuji':'Otsuji','goryachev':'Goryachev','legi':'LEGI','debelly':'DeBelly'}
FLABELS=['Gaussian','Memory','Coupled','Self-act.']
ALABELS={'pulse':'Pulse / relaxation','weak_plus':'Weak + perturbation','weak_minus':'Weak − perturbation','tracking_train':'Calibration tracking','tracking_test':'Unseen tracking period','reversal':'Cue reversal','competition':'Competing cues'}
REGIMES=['independent','constant','power']
COLORS=['#0072B2','#E69F00','#009E73','#CC79A7']
STYLES=['-','--','-.',':']

def heat(ax, matrix, labels, title, columns=FLABELS):
    cmap=plt.get_cmap('YlOrRd').copy();cmap.set_bad('#eeeeee')
    ax.imshow(matrix,cmap=cmap,vmin=0,vmax=1,aspect='auto')
    ax.set_xticks(range(len(columns)),columns,rotation=25,ha='right',fontsize=9)
    ax.set_yticks(range(len(labels)),labels,fontsize=9);ax.set_title(title,fontsize=12,pad=10)
    for i in range(len(matrix)):
        for j in range(len(columns)):
            v=matrix[i,j]
            ax.text(j,i,f'{v:.2f}' if np.isfinite(v) else '—',ha='center',va='center',fontsize=9,color='white' if v>.6 else '#172735')

def render(rows,viewer):
    def error(m,f,r,role,assays):
        vals=[x['post_nrmse'] for x in rows if x['model']==m and x['family']==f and x['regime']==r and x['length_role']==role and x['assay']in assays]
        return np.sqrt(np.mean(np.square(vals))) if vals else np.nan
    groups=[('Fitted experiments · fitted lengths','training',CALIBRATION),
            ('New experiments · fitted lengths','training',ASSAYS[4:]),
            ('All experiments · unseen lengths','held_out',ASSAYS)]
    fig,axes=plt.subplots(3,3,figsize=(15,10),layout='constrained')
    for i,(title,role,assays) in enumerate(groups):
        for j,regime in enumerate(REGIMES):
            ax=axes[i,j]
            if i==2 and j==0:
                ax.axis('off');ax.text(.5,.5,'Independent fits cannot predict\nunseen lengths.',ha='center',va='center',transform=ax.transAxes,color='#66717c');continue
            mat=np.array([[error(m,f,regime,role,assays)for f in FAMILIES]for m in MODELS])
            heat(ax,mat,[LABELS[m]for m in MODELS],regime.capitalize()+'\n'+title)
    fig.suptitle('Response error overview · lower is better',fontsize=18)
    fig.savefig(OUT/'error_overview.png',dpi=160);plt.close(fig)
    # Two optional detailed figures; the optimization objective stays in CSV only.
    for name,role,assays,regimes,title in [
        ('fit_detail','training',CALIBRATION,REGIMES,'Which fitted experiments are difficult?'),
        ('heldout_detail','held_out',ASSAYS,['constant','power'],'Which experiments fail at unseen lengths?')]:
        fig,axes=plt.subplots(1,len(regimes),figsize=(17,10 if name=='fit_detail' else 14),layout='constrained')
        for ax,regime in zip(axes,regimes):
            mat=np.array([[error(m,f,regime,role,[a])for f in FAMILIES]for m in MODELS for a in assays])
            heat(ax,mat,[LABELS[m]+' / '+ALABELS[a]for m in MODELS for a in assays],regime.capitalize())
            for y in range(len(assays),len(mat),len(assays)):ax.axhline(y-.5,color='#172735',lw=1)
        fig.suptitle(title,fontsize=18);fig.savefig(OUT/(name+'.png'),dpi=150);plt.close(fig)
    from matplotlib.lines import Line2D
    legend=[Line2D([0],[0],color='#111',lw=2.5,label='Spatial data')]+[Line2D([0],[0],color=c,ls=ls,lw=1.7,label=lab)for c,ls,lab in zip(COLORS,STYLES,FLABELS)]
    for m in MODELS:
        for group,assays in [('calibration',CALIBRATION),('heldout',ASSAYS[4:])]:
            for regime in REGIMES:
                fig,axes=plt.subplots(len(assays),5,figsize=(18,2.5*len(assays)+1.1),layout='constrained',squeeze=False)
                for i,a in enumerate(assays):
                    protocol=Protocol(m,a)
                    origin=0 if a in ('pulse','competition') else protocol.condition
                    # Pulse/competition retain the complete initial transient. Other
                    # assays focus on the perturbation after long conditioning.
                    start=0 if origin==0 else origin-.1*(protocol.end-origin)
                    rowvalues=[]
                    for L in lengths(m):
                        d=viewer.get('|'.join([m,str(L),a,regime]))
                        if d:
                            mask=np.array(d['time'])>=start
                            rowvalues.extend(np.array(v)[mask]for v in [d['target']]+list(d['families'].values()))
                    low=min(-.01,min(float(v.min())for v in rowvalues));high=max(.01,max(float(v.max())for v in rowvalues));pad=.06*(high-low)
                    for j,L in enumerate(lengths(m)):
                        ax=axes[i,j];key='|'.join([m,str(L),a,regime]);d=viewer.get(key)
                        known=L in lengths(m,True)
                        if i==0:ax.set_title(f'{L:g} µm · '+('fitted length'if known else 'unseen length'),fontsize=11,color='#172735'if known else '#884600',fontweight='bold')
                        if j==0:ax.set_ylabel(ALABELS[a]+'\nPolarity p',fontsize=10)
                        ax.set_ylim(low-pad,high+pad);ax.set_xlim(start-origin,protocol.end-origin)
                        ax.grid(alpha=.18);ax.tick_params(labelsize=8)
                        if not d:
                            ax.set_facecolor('#eeeeee');ax.text(.5,.5,'No independent fit\nat unseen length',ha='center',va='center',transform=ax.transAxes,fontsize=9,color='#66717c');continue
                        t=np.array(d['time']);mask=t>=start
                        ax.plot(t[mask]-origin,np.array(d['target'])[mask],color='#111',lw=2.3,alpha=.85)
                        for f,c,ls in zip(FAMILIES,COLORS,STYLES):ax.plot(t[mask]-origin,np.array(d['families'][f])[mask],color=c,ls=ls,lw=1.5)
                        if origin:ax.axvline(0,color='#aaa',lw=.8)
                        else:ax.axvspan(0,protocol.pulse,color='#999',alpha=.12)
                        ax.set_xlabel('Time (s)' if origin==0 else 'Seconds from cue change',fontsize=8)
                fig.suptitle(LABELS[m]+' · '+('fitted experiments'if group=='calibration'else 'held-out experiments')+' · '+regime+' fits',fontsize=17)
                fig.legend(handles=legend,loc='outside lower center',ncol=5,frameon=False,fontsize=11)
                fig.savefig(OUT/f'trajectories_{m}_{group}_{regime}.png',dpi=150);plt.close(fig)
        print('Facets:',m,flush=True)
    now=datetime.now(timezone.utc).isoformat()
    doc='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Spatial fits: overview and trajectory sheets</title>
<style>body{font:16px/1.5 system-ui;margin:28px auto;padding:0 22px;max-width:1500px;color:#172735}img{width:100%;height:auto}a{color:#075ca4}h1,h2{line-height:1.2}section{margin:42px 0}figure{margin:22px 0}figcaption{color:#475663}summary{cursor:pointer;font-weight:600;padding:12px;background:#eef3f6}details{margin:18px 0}button{font:inherit;padding:8px 18px;margin:4px;border:1px solid #8193a3;border-radius:4px;background:white;cursor:pointer}button[aria-pressed=true]{background:#172735;color:white}button:focus-visible,a:focus-visible{outline:3px solid #dc8600}.controls{position:sticky;top:0;background:#fffef6;padding:9px;z-index:2;border-bottom:1px solid #ccd4dc}.links{display:flex;flex-wrap:wrap;gap:18px}.note{max-width:1050px}.stamp{font-size:13px;color:#657381}</style>
<h1>Spatial fits: overview and trajectory sheets</h1>'''
    doc+=f'<p class="stamp">Generated {now}. Existing fitted parameters and simulations; presentation updated only.</p>'
    doc+='<nav class="links"><a href="report.html">WGD report</a><a href="#errors">Error overview</a><a href="#trajectories">Trajectory sheets</a></nav>'
    doc+='''<section id="errors"><h2>Where do the fits fail?</h2><p class="note">Read down to distinguish fitting errors from failures on new experiments or new lengths. Each cell condenses a spatial source and OU candidate across the relevant experiments and lengths. Read across to compare independent, constant and power-law fits.</p>
<p class="note"><b>One error metric throughout:</b> post-event RMSE normalized by the calibration-amplitude scale, then RMS across included experiments and lengths. Lower is better; colors cap at 1, printed values do not. The bottom row pools all seven experiments at the two unseen lengths.</p>
<a href="error_overview.png"><img src="error_overview.png" alt="Condensed response-error overview: fitting, new experiments and unseen lengths, across independent, constant and power-law fits"></a>
<details><summary>Experiment breakdown: fitted experiments</summary><p>Four calibration experiments per spatial source; RMS across the three fitted lengths.</p><a href="fit_detail.png"><img loading="lazy" src="fit_detail.png" alt="Training errors broken down by source and calibration experiment"></a></details>
<details><summary>Experiment breakdown: unseen lengths</summary><p>All seven experiments per spatial source; RMS across the two unseen lengths.</p><a href="heldout_detail.png"><img loading="lazy" src="heldout_detail.png" alt="Held-out-length errors broken down by source and experiment"></a></details>
<p class="note">Independent fits use separate parameters at each fitted length. Constant fits share parameters across lengths. Power fits use the best tested nonconstant law selected on calibration error; selection does not establish evidence for scaling. The optimizer-objective heatmap has been removed: its different weighting obscured comparison and duplicated the training view.</p>
<p><a href="spatial_errors.csv">All per-experiment errors (CSV)</a> · <a href="training_objective_errors.csv">Optimizer-objective components (CSV)</a></p></section>
<section id="trajectories"><h2>Ten trajectory sheets</h2><p class="note">Each spatial source has two images: <b>20 calibration facets</b> and <b>15 held-out-experiment facets</b>. Rows are experiments; columns are lengths. Each facet overlays the spatial trajectory and all four OU candidates. Change the fit once to update all ten sheets. Click any image to inspect it at full resolution.</p>
<p class="note">Vertical scales match across lengths within each experiment and fit. Pulse and competition show the complete trajectory; other experiments show the response around the cue change after conditioning. Unseen lengths have amber column titles. Independent fits leave those columns empty.</p>
<div class="controls" role="group" aria-label="Fit regime for all trajectory sheets">Show all sheets with: <button data-regime="independent" aria-pressed="false">Independent</button><button data-regime="constant" aria-pressed="false">Constant</button><button data-regime="power" aria-pressed="true">Power law</button><span id="fit-status" aria-live="polite">Power-law fits shown</span></div><nav class="links">'''
    doc+=''.join(f'<a href="#sheets-{m}">{LABELS[m]}</a>'for m in MODELS)+'</nav>'
    for m in MODELS:
        doc+=f'<section id="sheets-{m}"><h3>{LABELS[m]}</h3>'
        for group,label in [('calibration','Fitted experiments · 20 facets'),('heldout','Held-out experiments · 15 facets')]:
            stem=f'trajectories_{m}_{group}'
            doc+=f'<figure><figcaption><b>{label}</b></figcaption><a class="sheet-link" href="{stem}_power.png"><img class="sheet" data-stem="{stem}" loading="lazy" src="{stem}_power.png" alt="{LABELS[m]}: {label}; five lengths with spatial and four OU trajectories"></a></figure>'
        doc+='</section>'
    doc+='''</section><script>document.querySelectorAll('button[data-regime]').forEach(button=>button.addEventListener('click',()=>{const r=button.dataset.regime;document.querySelectorAll('button[data-regime]').forEach(b=>b.setAttribute('aria-pressed',String(b===button)));document.querySelectorAll('img.sheet').forEach(img=>{const path=img.dataset.stem+'_'+r+'.png';img.src=path;img.parentElement.href=path;});document.getElementById('fit-status').textContent=button.textContent+' fits shown';}));</script></html>'''
    temp=OUT/'spatial.html.tmp';temp.write_text(doc);temp.replace(OUT/'spatial.html')

if __name__=='__main__':build()
