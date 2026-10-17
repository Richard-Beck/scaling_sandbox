"""Post-processing for normalized gradient reversal and finite countercue."""
from pathlib import Path
from datetime import datetime,timezone
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
from matplotlib.colors import TwoSlopeNorm
HERE=Path('/home/4473331/projects/scaling_sandbox/polarity1d/four_assay')
import sys
sys.path.insert(0,str(HERE))
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,default=HERE/'results');args=ap.parse_args()
OUT=args.results.resolve()
MODELS=['wave_pinning','otsuji','goryachev','legi','debelly','holmes_model3','spring']
LABELS=['Wave-pinning','Otsuji','Goryachev','LEGI','deBelly','Holmes 3','Spring']
TIME_COLORS=['#176d9c','#db7b24','#27804b']
MODES={'ineligible':('#d9d9d9','Ineligible state'), 'unsettled':('#c8b6e2','Preparation unsettled'),
       'censored':('#f0c6ab','Endpoint not observed'), 'not_applicable':('#eeeeee','Not applicable')}
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
METRICS=[('settle_conditioning','Conditioning arrival'),('settle_finite_zero','Release arrival'),
         ('gradient_reversal','Gradient response'),('finite_countercue','Countercue duration')]

def readcsv(p):return list(csv.DictReader(p.open()))
def writecsv(p,rows):
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def number(x):return float(x) if x not in [None,''] else np.nan

def load():
    import sys
    sys.path.insert(0,str(HERE))
    import runner
    return runner.export_tables(OUT)


def intervals(lengths):
    x=np.asarray(sorted(set(lengths)),dtype=float)
    step=(x[1]-x[0])/2 if len(x)>1 else .5
    return x,np.r_[x[0]-step,(x[:-1]+x[1:])/2,x[-1]+step]

def shade(ax,x,modes,ymin=0,ymax=1):
    xx,b=intervals(x)
    for L,lo,hi in zip(xx,b[:-1],b[1:]):
        mode=modes.get(L,'not_applicable')
        if mode!='success':
            color=MODES.get(mode,MODES['not_applicable'])[0]
            ax.axvspan(lo,hi,ymin=ymin,ymax=ymax,color=color,lw=0,zorder=-3,
                       hatch=None)
    ax.set_xlim(b[0],b[-1])

def regression(model,metric,x,y):
    x=np.asarray(x,float);y=np.asarray(y,float);good=np.isfinite(x)&np.isfinite(y);x=x[good];y=y[good]
    rec=dict(model=model,metric=metric,n_valid=len(x),fit_eligible=len(x)>3,slope=None,intercept=None,r_squared=None,
             min_valid_length=float(x.min()) if len(x) else None,max_valid_length=float(x.max()) if len(x) else None)
    if len(x)>3:
        xc=x-x.mean();slope=np.dot(xc,y-y.mean())/np.dot(xc,xc);intercept=y.mean()-slope*x.mean()
        fitted=intercept+slope*x;sse=np.sum((y-fitted)**2);sst=np.sum((y-y.mean())**2)
        rec.update(slope=float(slope),intercept=float(intercept),r_squared=float(1-sse/sst) if sst>1e-20 else None)
    return rec

def plot_series(ax,model,metric,x,y,color,records,label=None):
    ax.plot(x,y,'o-',color=color,lw=1.05,ms=2.5,label=label)
    fit=regression(model,metric,x,y);records[(model,metric)]=fit
    if fit['fit_eligible']:
        xx=np.linspace(fit['min_valid_length'],fit['max_valid_length'],100)
        ax.plot(xx,fit['intercept']+fit['slope']*xx,'--',color=color,lw=1.5,alpha=.95)

def time_panel(ax,model,metric,lengths,byL,records,preparation=False):
    if preparation:
        modes={L:('success' if byL[L]['settled']=='True' else 'unsettled') if L in byL else 'not_applicable' for L in lengths}
        y=[number(byL[L]['time_to_settle']) if modes[L]=='success' else np.nan for L in lengths]
    else:
        modes={L:byL[L]['mode'] for L in lengths}
        y=[number(byL[L]['endpoint']) if modes[L]=='success' else np.nan for L in lengths]
    shade(ax,lengths,modes)
    plot_series(ax,model,metric,lengths,y,TIME_COLORS[0],records)
    ax.set_ylim(bottom=0)
    if not np.isfinite(y).any():ax.set_ylim(0,1)
    units='native units' if model=='spring' else 'model units' if model in MODELS[:4] else 's'
    ax.set_ylabel('Time ('+units+')');ax.tick_params(labelsize=7)


def decorate(fig,axes,titles,stamp,heading):
    for ci,title in enumerate(titles):axes[0,ci].set_title(title,fontsize=10,pad=8)
    for ri,label in enumerate(LABELS):axes[ri,0].set_ylabel(label+'\n'+axes[ri,0].get_ylabel(),fontsize=9)
    for ax in axes[-1,:]:ax.set_xlabel('Resting L₀ (spring); length µm otherwise',fontsize=7)
    handles=[Patch(facecolor=color,label=label,hatch=None) for mode,(color,label) in MODES.items() if mode!='not_applicable']
    fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,fontsize=9,bbox_to_anchor=(.5,.008))
    fig.suptitle(f'{heading} | generated {stamp}\nSettling: first persistent 1% arrival. Responses start from verified settled profiles. Dashed: OLS, >3 valid points.',fontsize=11)
    fig.subplots_adjust(left=.07,right=.99,bottom=.055,top=.94,wspace=.34,hspace=.43)


def export(fig,name,dpi=180):
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=dpi)
    plt.close(fig)


def figures(rows,settling,stamp):
    records={}
    titles=['Conditioning: 1% arrival (shared)','Release: 1% arrival (countercue)','Gradient reversal response','Countercue pulse duration']
    fig,axes=plt.subplots(7,4,figsize=(20,19))
    for i,m in enumerate(MODELS):
        lengths=sorted({number(r['length']) for r in rows if r['model']==m})
        ss=[r for r in settling if r['model']==m]
        gradient={number(r['length']):r for r in ss if r['assay']=='gradient_reversal' and r['phase_index']=='0'}
        finite={number(r['length']):r for r in ss if r['assay']=='finite_countercue' and r['phase_index']=='0'}
        for L in lengths:
            if L not in gradient or L not in finite:raise ValueError(f'Missing conditioning stage: {m} L={L}')
            for field in ['settled','time_to_settle']:
                if gradient[L][field]!=finite[L][field]:raise ValueError(f'Conditioning is not shared: {m} L={L} {field}')
        time_panel(axes[i,0],m,'settle_conditioning',lengths,gradient,records,preparation=True)
        zero={number(r['length']):r for r in ss if r['assay']=='finite_countercue' and r['cue']=='zero'}
        time_panel(axes[i,1],m,'settle_finite_zero',lengths,zero,records,preparation=True)
        for ci,assay in enumerate(['gradient_reversal','finite_countercue'],start=2):
            ds={number(r['length']):r for r in rows if r['model']==m and r['assay']==assay}
            time_panel(axes[i,ci],m,assay,lengths,ds,records)
    decorate(fig,axes,titles,stamp,'Response cues A*x; original A*x/L preparations')
    export(fig,'length_response_and_settling',200)
    return [records[(m,k)] for m in MODELS for k,_ in METRICS]


def normalize(records):
    for r in records:
        r.update(characteristic_length=None,characteristic_response=None,normalized_slope=None)
        if r['fit_eligible']:
            length=(r['min_valid_length']+r['max_valid_length'])/2
            response=r['intercept']+r['slope']*length
            r.update(characteristic_length=length,characteristic_response=response)
            if np.isfinite(response) and abs(response)>1e-12:
                r['normalized_slope']=r['slope']*length/response


def heatmap(records,stamp):
    lookup={(r['model'],r['metric']):r for r in records}
    values=[r['normalized_slope'] for r in records if r['normalized_slope'] is not None]
    lim=max(max(abs(v) for v in values) if values else 0,1e-12)
    norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim);cmap=plt.get_cmap('RdBu_r')
    fig=plt.figure(figsize=(11,7))
    grid=fig.add_gridspec(1,len(METRICS),wspace=.12)
    for j,(metric,title) in enumerate(METRICS):
        ax=fig.add_subplot(grid[0,j])
        for i,m in enumerate(MODELS):
            value=lookup[(m,metric)]['normalized_slope'];ok=value is not None and np.isfinite(value)
            ax.add_patch(Rectangle((0,i),1,1,facecolor=cmap(norm(value)) if ok else '#d9d9d9',edgecolor='white',lw=1))
            ax.text(.5,i+.5,f'{value:.3g}' if ok else '—',ha='center',va='center',fontsize=10,
                    color='white' if ok and abs(value)/lim>.65 else '#333333')
        ax.set(xlim=(0,1),ylim=(7,0),xticks=[],yticks=np.arange(7)+.5)
        ax.set_yticklabels(LABELS if j==0 else []);ax.tick_params(axis='y',length=0)
        ax.set_title(title.replace(' ','\n',1),fontsize=10)
        for spine in ax.spines.values():spine.set_visible(False)
    cax=fig.add_axes([.25,.09,.55,.035])
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal')
    cax.set_xlabel('Normalized slope = slope × characteristic length / fitted response (dimensionless)',fontsize=9)
    fig.suptitle(f'Normalized length slopes | generated {stamp}\n'
                  'Characteristic length: midpoint of each fit’s valid range. One color scale.\n'
                  'Gray: ≤3 valid points, unavailable fit, or zero reference response. Raw slopes retained in CSV.',fontsize=11)
    fig.subplots_adjust(left=.14,right=.99,bottom=.23,top=.77)
    export(fig,'length_response_gradients',220)


def main():
    stamp=datetime.now(timezone.utc).isoformat();rows,settling=load()
    records=figures(rows,settling,stamp)
    normalize(records);heatmap(records,stamp)
    writecsv(OUT/'endpoints.csv',rows);writecsv(OUT/'length_regressions.csv',records)
    counts={}
    for r in rows:counts[r['status']]=counts.get(r['status'],0)+1
    provenance=dict(generated_at=stamp,plot_sha256=sha(Path(__file__)),records_sha256=sha(OUT/'records.json.gz'),
                    counts=counts,eligible_regressions=sum(r['fit_eligible'] for r in records),
                    regression='Unweighted OLS with intercept, >3 valid observations; shared conditioning fitted once; release settling and assay responses fitted separately; no cumulative fits',
                    normalization='slope * characteristic_length / fitted_response; characteristic_length is midpoint of each fit valid range; zero reference response unavailable',
                    response_time='Gradient endpoint excludes conditioning; finite endpoint is passing pulse duration, excluding conditioning, release settling and recovery observation',
                    settling_time='First persistent 1% neighborhood entry on fixed field scales; shared conditioning from native input; release from settled conditioning',
                    response_inputs='Gradient: settled conditioning; finite countercue: settled release; verified continuation plus holdout',
                    heatmap='Four dimensionless component slopes; shared conditioning once; one symmetric color scale; unavailable fit or zero reference gray',source_revision=11)
    (OUT/'visualization_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    companion = ''
    if (OUT/'length_response_profiles_provenance.json').exists():
        profile_info = json.loads((OUT/'length_response_profiles_provenance.json').read_text())
        companion = (f'<h2>Companion spatial profiles</h2><p>Generated {profile_info["generated_at"]}; '
                     f'results generation {profile_info["results_generated_at"]}; producer finished.</p>'
                     '<p><a href="length_response_profiles.pdf">PDF</a> · '
                     '<a href="LENGTH_RESPONSE_PROFILES_NOTES.md">Notes</a></p>'
                     '<img style="width:100%" src="length_response_profiles.png">')
    (OUT/'index.html').write_text(f'<html><body><h1>Normalized gradient reversal and finite countercue</h1><p>Generated {stamp}; producer finished.</p><img style="width:100%" src="length_response_and_settling.png"><img style="width:100%" src="length_response_gradients.png">{companion}</body></html>')
    exports=[f for f in OUT.iterdir() if f.suffix in ['.png','.pdf','.svg','.csv']]
    (OUT/'artifact_hashes.json').write_text(json.dumps({f.name:sha(f) for f in exports},indent=2)+'\n')
    print(json.dumps(provenance))
if __name__=='__main__':main()
