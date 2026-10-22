"""Review-only Figure4 round-2 plots; this does not compose the main figure."""
from __future__ import annotations
import json, math, shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import revision_round2_figure4_rates_20261022 as R

OUT=R.OUT;DATA=OUT/'supporting/data';VEC=OUT/'supporting/vectors';META=OUT/'supporting/metadata'
def n(r,k):return R.num(r.get(k))
def save(fig,stem,png=True):
    fig.savefig(VEC/(stem+'.pdf'),bbox_inches='tight');fig.savefig(VEC/(stem+'.svg'),bbox_inches='tight')
    if png:fig.savefig(OUT/(stem+'.png'),dpi=160,bbox_inches='tight')
def finite(x):return np.isfinite(x)
def bound_marker(ax,x,lo,hi,color,offset=1):
    x=x*offset
    if lo>0 and finite(lo) and finite(hi) and hi>lo:
        ax.plot([x,x],[lo,hi],color=color,lw=1)
        ax.plot(x,lo,'^',mfc='none',color=color,ms=5);ax.plot(x,hi,'v',mfc='none',color=color,ms=5)
        return
    if lo>0 and finite(lo):
        ax.plot(x,lo,'^',mfc='none',color=color,ms=5);ax.annotate('',xy=(x,lo*1.7),xytext=(x,lo),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
    elif finite(hi) and hi>0:
        ax.plot(x,hi,'v',mfc='none',color=color,ms=5);ax.annotate('',xy=(x,hi*.6),xytext=(x,hi),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
    else:ax.plot(x,.025,'x',color=color,ms=4)

def ratio_plot(rows,metric,view='main'):
    fig,axes=plt.subplots(1,4,figsize=(15,4),sharey=True)
    for ax,model in zip(axes,R.MODELS):
        for arm,style,offset in [('common','--',.985),('longest','-',1.015)]:
            rr=sorted([r for r in rows if r['model']==model and r['arm']==arm and r['view']==view],key=lambda r:n(r,'f'))
            f=np.array([n(r,'f') for r in rr]);y=np.array([n(r,metric) for r in rr]);ax.plot(f,y,style,color=R.COLORS[model],lw=1.2,label=arm)
            for r,x,val in zip(rr,f,y):
                color='#999999' if R.yes(r.get('forced_combined')) else R.COLORS[model]
                lo=n(r,metric+'_ci_lower');hi=n(r,metric+'_ci_upper')
                if finite(val) and val>0:
                    if finite(lo) and finite(hi) and lo>0:ax.plot([x*offset]*2,[lo,hi],color=color,alpha=.45,lw=1)
                    marker={'single_front':'o','multiple_fronts':'s','low_contrast':'D'}.get(r.get('large_regime'),'v')
                    ax.plot(x*offset,val,marker,ms=4,color=color,mfc='white' if R.yes(r.get('non_exponential')) else color)
                    if math.isinf(hi):ax.annotate('',xy=(x*offset,min(val*2,1e7)),xytext=(x*offset,val),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
                elif r.get('status')=='one_sided_lower_bound':
                    b=n(r,'rho_k_lower')/(n(r,'geometry_expectation') if metric=='rho_geo' else 1)
                    if finite(b) and b>0:ax.plot(x,b,'^',mfc='none',color=color);ax.annotate('',xy=(x,b*2),xytext=(x,b),arrowprops=dict(arrowstyle='->',color=color))
                if not R.yes(r.get('evaluable')):ax.plot(x,.025,'x',color='#aaaaaa',ms=3)
        ax.axhline(1,color='black',lw=.8)
        if metric=='rho_k':
            expectation=[n(r,'geometry_expectation') for r in rows if r['model']==model and r['view']==view]
            expectation=[x for x in expectation if finite(x)]
            ax.axhline(np.median(expectation) if expectation else 9,color='#888888',lw=.7,ls=':')
        ax.set(xscale='log',yscale='log',xlabel='Input fraction f',title=R.NAMES[model]);ax.grid(alpha=.15)
    axes[0].set_ylabel(r'$\rho_k = k_{small}/k_{large}$' if metric=='rho_k' else r'$\rho_{geo}$ relative to single-front geometry')
    axes[0].legend(fontsize=8);fig.suptitle(('Heading decay-rate size ratio' if metric=='rho_k' else 'Size advantage relative to single-front geometric expectation')+' — '+view)
    fig.text(.5,-.01,'Circles: single front; squares: multiple fronts; diamonds: low contrast; grey: forced. Open: >1.5× lag drift. Cross: unavailable. Upward arrows: unbounded CI. Dotted: actual geometry (nominal 9).',ha='center',fontsize=7)
    fig.tight_layout();return fig

def tolerance_plot(rows,view='main',inclusive=False):
    fig,axes=plt.subplots(1,4,figsize=(14,4),sharey=True)
    for ax,model in zip(axes,R.MODELS):
        for arm,style in [('common','--'),('longest','-')]:
            rr=sorted([r for r in rows if r['model']==model and r['arm']==arm and r['view']==view and R.yes(r['forced_inclusive'])==inclusive],key=lambda r:n(r,'C_label'))
            xs=[n(r,'C_label') for r in rr];ys=[n(r,'fhalf') for r in rr];ax.plot(xs,ys,style,color=R.COLORS[model],label=arm)
            for r,x,y in zip(rr,xs,ys):
                lo=n(r,'fhalf_ci_lower');hi=n(r,'fhalf_ci_upper')
                if finite(y):
                    ax.plot(x,y,'o',color=R.COLORS[model]);
                    if lo>0 and finite(hi):ax.plot([x,x],[lo,hi],color=R.COLORS[model],alpha=.5)
                    elif math.isinf(hi):ax.annotate('',xy=(x,y*2),xytext=(x,y),arrowprops=dict(arrowstyle='->',color=R.COLORS[model]))
                else:
                    lower=n(r,'fhalf_lower');upper=n(r,'fhalf_upper');bound_marker(ax,x,lower,upper,R.COLORS[model])
        ax.set(title=R.NAMES[model],xlabel='Designed circumference label',yscale='log');ax.grid(alpha=.15)
    axes[0].set_ylabel(r'$f_{1/2}$: input at $C(\ell_{ref})=0.5$');axes[0].legend(fontsize=8)
    fig.suptitle(f'Noise tolerance — {view}; forced conditions '+('included' if inclusive else 'excluded'))
    fig.text(.5,-.01,'Spring labels denote resting/design circumference. Open triangles+arrows: bounds; paired triangles: finite interval. Crosses: unavailable.',ha='center',fontsize=8)
    fig.tight_layout();return fig

def correlations(rows,pdf,page_selection=None):
    for arm in ('common','longest'):
        for view in R.VIEWS:
            for model in R.MODELS:
                if page_selection is not None and (arm,view,model) not in page_selection:continue
                selected=[r for r in rows if r['arm']==arm and r['view']==view and r['model']==model]
                sizes=sorted({n(r,'C_label') for r in selected});fs=sorted({n(r,'f') for r in selected})
                fig,axes=plt.subplots(3,len(sizes),figsize=(11,9),squeeze=False)
                for col,C in enumerate(sizes):
                    for i,f in enumerate(fs):
                        rr=sorted([r for r in selected if n(r,'C_label')==C and n(r,'f')==f],key=lambda r:n(r,'lag'))
                        lag=[n(r,'lag') for r in rr];color=plt.get_cmap('viridis')(i/max(1,len(fs)-1))
                        for row,key in enumerate(('C','k','coverage')):
                            ax=axes[row,col];ax.plot(lag,[n(r,key) for r in rr],'.-',color=color,label=f'{f:g}',lw=.8,ms=3)
                            if key in ('C','k'):ax.fill_between(lag,[n(r,key+'_ci_lower') for r in rr],[n(r,key+'_ci_upper') for r in rr],color=color,alpha=.08)
                            ax.set_xscale('log');ax.grid(alpha=.12)
                    axes[0,col].set_title(f'Design circumference {C:g}')
                    axes[0,col].axhline(.98,c='grey',lw=.6,ls=':');axes[0,col].set_ylim(-1.05,1.05)
                    axes[1,col].set_yscale('log');axes[2,col].axhline(.5,c='grey',lw=.6,ls=':');axes[2,col].set_ylim(-.03,1.03);axes[2,col].set_xlabel('Lag (native units)')
                for row,label in enumerate(('Heading autocorrelation C','Rotational decay rate k','Eligible endpoint coverage')):axes[row,0].set_ylabel(label)
                handles,labels=axes[0,-1].get_legend_handles_labels();fig.legend(handles,labels,title='Input f',loc='center right',fontsize=7)
                fig.suptitle(f'{R.NAMES[model]} — {arm}, {view}\nWhole-seed 95% intervals; R≥0.1 at both endpoints. Blank k: C≤0; CI upper unbounded if C lower≤0.',fontsize=11)
                fig.tight_layout(rect=[0,0,.91,.94]);pdf.savefig(fig);plt.close(fig);print('C/k page',arm,view,model,flush=True)

def main():
    ratios=R.rows(DATA/'size_ratios.csv');tols=R.rows(DATA/'noise_tolerance.csv');lags=R.rows(DATA/'condition_lag_metrics.csv')
    for metric in ('rho_k','rho_geo'):
        fig=ratio_plot(ratios,metric);save(fig,'figure4_'+metric);plt.close(fig)
    fig=tolerance_plot(tols);save(fig,'figure4_noise_tolerance');plt.close(fig)
    pending_supplement=VEC/'figure4_reanalysis_supplement.rendering.pdf'
    with PdfPages(pending_supplement) as pdf:
        for view in R.VIEWS:
            for metric in ('rho_k','rho_geo'):
                fig=ratio_plot(ratios,metric,view);pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
            for inc in (False,True):
                fig=tolerance_plot(tols,view,inc);pdf.savefig(fig);plt.close(fig)
        correlations(lags,pdf)
    pending_supplement.replace(VEC/'figure4_reanalysis_supplement.pdf')
    old=R.rows(R.RUN/'results/final_analysis/size_ratios.csv');om={(r['arm'],r['model'],r['view'],r['f']):r for r in old};comparison=[]
    for r in ratios:
        o=om[(r['arm'],r['model'],r['view'],r['f'])];comparison.append(dict(r,old_tau_ratio=o.get('tau_ratio'),old_tau_ratio_lower=o.get('tau_ratio_lower'),old_tau_ratio_upper=o.get('tau_ratio_upper'),old_tau_ratio_ci_lower=o.get('tau_ratio_ci_lower'),old_tau_ratio_ci_upper=o.get('tau_ratio_ci_upper'),old_tau_evaluable=o.get('evaluable')))
    R.write(DATA/'old_tau_vs_new_rate.csv',comparison)
    fig,axes=plt.subplots(1,4,figsize=(14,4),sharey=True)
    for ax,model in zip(axes,R.MODELS):
        rr=[r for r in comparison if r['arm']=='longest' and r['view']=='main' and r['model']==model]
        ax.plot([n(r,'f') for r in rr],[n(r,'old_tau_ratio') for r in rr],'o--',label='Old tau point',color='#888888')
        ax.plot([n(r,'f') for r in rr],[n(r,'rho_k') for r in rr],'s-',label='New rate ratio',color=R.COLORS[model]);ax.axhline(1,color='black',lw=.6)
        for r in rr:
            if not finite(n(r,'old_tau_ratio')):bound_marker(ax,n(r,'f'),n(r,'old_tau_ratio_lower'),n(r,'old_tau_ratio_upper'),'#888888',.98)
            if not finite(n(r,'rho_k')):bound_marker(ax,n(r,'f'),n(r,'rho_k_lower'),n(r,'rho_k_upper'),R.COLORS[model],1.02)
        ax.set(xscale='log',yscale='log',title=R.NAMES[model],xlabel='Input f');ax.grid(alpha=.15)
    axes[0].legend(fontsize=8);axes[0].set_ylabel('Size ratio');fig.suptitle('Old persistence and new rate estimands — longest, main\nCensored values: open triangles and arrows. Crosses at 0.025: unavailable/uninformative interval.')
    fig.tight_layout();save(fig,'figure4_old_new');plt.close(fig)
    for name in ('criteria.csv','size_ratios.csv','noise_tolerance.csv','window_qualification.csv','horizon_qualification.csv'):
        p=R.RUN/'results/final_analysis'/name
        if p.exists():shutil.copy2(p,DATA/('legacy_'+name))
    source=R.ROOT/'figures/length_mechanisms_20261008/dev/signature_figure4_20261020'
    shutil.copy2(source/'supporting/vectors/figure4_main_selected_components_20261020_panelC.pdf',VEC/'legacy_seven_view_robustness.pdf')
    shutil.copy2(source/'supporting/vectors/figure4_S0_S12_final_decisions_20261020.pdf',VEC/'legacy_S0_S12_decisions.pdf')
    sensitivity=[]
    for p in source.rglob('*seven*'):
        if p.is_file():sensitivity.append(str(p.relative_to(OUT.parent)))
    (META/'render_receipt.json').write_text(json.dumps(dict(review_only=True,no_main_figure_built=True,rendered_utc=R.datetime.now(R.timezone.utc).isoformat(),source_analysis_sha256=R.sha(META/'analysis_receipt.json'),outputs={str(p.relative_to(OUT)):R.sha(p) for p in [*VEC.glob('figure4*'),*OUT.glob('figure4*.png')]}),indent=2)+'\n')

if __name__=='__main__':main()
