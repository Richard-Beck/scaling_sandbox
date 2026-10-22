"""Author-authorized Case A composition and explicitly post hoc H1 amendment.

Read-only native inputs; the frozen interim decisions are never overwritten.
"""
from __future__ import annotations
import json, math, csv, gzip
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfWriter
import revision_round2_figure4_rates_20261022 as R
import figure4_round2_decisions_20261022 as D

OUT=R.OUT; DATA=OUT/'supporting/data'; META=OUT/'supporting/metadata'; VEC=OUT/'supporting/vectors'
MODELS=('wave_pinning','holmes_model3','spring_ring')
MESSAGE='Enlarged cells hold heading longer under fluctuating input — at or above the geometric expectation at weak input — and the advantage declines as input strengthens.'
def n(r,k): return R.num(r.get(k))
def dump(path,value): path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=True)+'\n')
def edges(t):
    mid=(t[:-1]+t[1:])/2
    return np.r_[t[0]-(mid[0]-t[0]),mid,t[-1]+(t[-1]-mid[-1])]

def amendment(ratios,tols):
    frozen=D.read(DATA/'S0_old_new_criteria.csv'); lookup={D.key(r):r for r in frozen}; result=[]; eligibility=[]
    for arm in ('common','longest'):
        for view in R.VIEWS:
            for inc in (False,True):
                calls=[]
                for model in R.MODELS:
                    ts=[r for r in tols if (r['arm'],r['view'],r['model'],R.yes(r['forced_inclusive']))==(arm,view,model,inc)]
                    ts.sort(key=lambda r:n(r,'C_label')); cutoff=n(ts[-1],'fhalf')
                    allr=[r for r in ratios if (r['arm'],r['view'],r['model'])==(arm,view,model)]
                    rr=[D.mapped_ratio(r) for r in allr if n(r,'f')<cutoff]
                    if math.isfinite(cutoff): call=D.F.h1_model(rr,ts[0],ts[-1],inclusive=inc)
                    else: call=D.F.decision('H1_model','inconclusive','Large-cell fhalf is unresolved; amended domain cannot be selected')
                    calls.append(call)
                    previous=lookup[arm,view,inc,model,'H1_model']
                    result.append(dict(arm=arm,view=view,forced_inclusive=inc,model=model,**call,
                        frozen_rate_status=previous['status'],old_tau_status=previous['old_tau_status'],post_hoc=True,
                        domain='Gains and no-reversal checks restricted to f < point estimate of largest size fhalf; CI and tolerance rules unchanged',large_fhalf=cutoff,
                        retained_inputs=[n(r,'f') for r in rr]))
                    for r in allr: eligibility.append(dict(arm=arm,view=view,forced_inclusive=inc,model=model,f=r['f'],large_fhalf=cutoff,
                        below_large_fhalf=math.isfinite(cutoff) and n(r,'f')<cutoff,evaluable=r['evaluable'],forced_combined=r['forced_combined']))
                overall=D.F.at_least('H1_overall',calls,3)
                previous=lookup[arm,view,inc,'all','H1_overall']
                result.append(dict(arm=arm,view=view,forced_inclusive=inc,model='all',**overall,
                    frozen_rate_status=previous['status'],old_tau_status=previous['old_tau_status'],post_hoc=True))
    D.write(DATA/'S0_posthoc_H1.csv',result);R.write(DATA/'S0_posthoc_input_domain.csv',eligibility)
    retired=[dict(r,status='retired',frozen_rate_status=r['status'],reason='Primary estimand changed; tau-based H2/H3/A6 gates no longer apply')
        for r in frozen if r['criterion'].startswith(('H2','H3'))]
    D.write(DATA/'S0_retired_gates.csv',retired)
    # A6 was nested in the frozen H2 record. Preserve it explicitly as provenance.
    a6=[]
    for r in frozen:
        if r['criterion']=='H2_overall':
            for c in json.loads(r['components']):
                if c['criterion']=='H2_A6_predictor_gate':a6.append(dict(arm=r['arm'],view=r['view'],forced_inclusive=r['forced_inclusive'],**dict(c,frozen_status=c['status'],status='retired')))
    D.write(DATA/'S0_retired_A6.csv',a6)
    with PdfPages(OUT/'figure4_S0_posthoc.pdf') as pdf:
        for arm in ('common','longest'):
            for inc in (False,True):
                rr=[r for r in result if r['arm']==arm and r['forced_inclusive']==inc]
                fig,ax=plt.subplots(figsize=(13,6));ax.axis('off')
                cells=[]
                for m in (*R.MODELS,'all'):
                    cells.append([R.NAMES.get(m,'Overall')]+[next(r for r in rr if r['model']==m and r['view']==v)['frozen_rate_status']+' → '+next(r for r in rr if r['model']==m and r['view']==v)['status'] for v in R.VIEWS])
                table=ax.table(cellText=cells,colLabels=['H1']+[v.replace('analysis_','').replace('_',' ') for v in R.VIEWS],loc='center',cellLoc='center');table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,2.8)
                ax.set_title(f'S0: frozen rate H1 → post hoc amended H1\n{arm}; forced conditions '+('included' if inc else 'excluded'),fontsize=14)
                fig.text(.07,.13,'Post hoc amendment: evaluate H1 gains only at f strictly below the largest size’s estimated f½.\nNoise-tolerance CI separation and all other H1 rules remain unchanged. Original τ calls remain in the frozen supplement.\nH2, H3 and A6 gates are retired. ρgeo and multiple-front fractions are descriptive outcomes.',fontsize=10)
                pdf.savefig(fig);plt.close(fig)
    return result

def native_examples():
    manifest=R.read(R.ROOT/'data/runs/signature_figure4_E2_analysis_20261020/arms/sigma_1p0/analysis_manifest.json')
    examples=[]
    for model,C,f in [('wave_pinning',20,.003),('wave_pinning',60,.003),('holmes_model3',20,.003),('holmes_model3',60,.003),('wave_pinning',60,.1)]:
        group=next(g for g in manifest['groups'] if g['model']==model and g['C_label']==C and g['f']==f)
        case=sorted(group['noise'],key=lambda r:int(r['seed']))[0]
        p=Path(case['profile_path']);assert R.sha(p)==case['profile_sha256']
        with np.load(p,allow_pickle=False) as z: t=z['time'];a=z['primary'];circ=float(np.broadcast_to(z['circumference'],t.shape).mean())
        # Preserve all saved spatial cells and native time samples in the viewport.
        keep=t>=t[-1]*.5;t=t[keep];a=a[keep]
        examples.append(dict(model=model,C=C,f=f,time=t,field=a,circ=circ,source=case))
    dump(META/'main_native_examples.json',[dict(model=e['model'],C=e['C'],f=e['f'],source=e['source'],selection='Lowest numeric seed in pre-existing E2 group; second half of native history',shape=list(e['field'].shape)) for e in examples])
    with gzip.open(DATA/'main_native_kymographs.csv.gz','wt',newline='') as stream:
        w=csv.writer(stream);w.writerow(['model','C_label','f','time','native_cell','arc_position','activity'])
        for e in examples:
            for t,a in zip(e['time'],e['field']):
                w.writerows((e['model'],e['C'],e['f'],t,i,(i+.5)*e['circ']/len(a),v) for i,v in enumerate(a))
    return examples

def marker(ax,x,y,lo,hi,color,shape='o',fill=None):
    if math.isfinite(y) and y>0:
        ax.plot([x,x],[max(lo,.3) if math.isfinite(lo) else y,hi if math.isfinite(hi) else y],color=color,alpha=.4,lw=.9)
        ax.plot(x,y,shape,color=color,mfc=fill or color,ms=4)
        if math.isinf(hi):ax.annotate('',(x,y*1.8),(x,y),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
        if lo==0:ax.annotate('',(x,.31),(x,max(.5,y*.7)),arrowprops=dict(arrowstyle='->',color=color,lw=.8))

def lag_consistency_supplement():
    """Display every usable lag, retaining undefined logarithms explicitly."""
    rows=R.rows(DATA/'all_lag_size_ratios.csv'); pages=[]
    with PdfPages(OUT/'figure4_all_lag_ratios.pdf') as pdf:
        for arm in ('common','longest'):
            for view in R.VIEWS:
                for model in R.MODELS:
                    rr=[r for r in rows if (r['arm'],r['view'],r['model'])==(arm,view,model)]
                    inputs=sorted({n(r,'f') for r in rr})
                    fig,axes=plt.subplots(1,2,figsize=(14,8))
                    display_limits={}
                    for ax,metric in zip(axes,('rho_k','rho_geo')):
                        points=[n(r,metric) for r in rr if R.yes(r['usable']) and np.isfinite(n(r,metric)) and n(r,metric)>0]
                        lower=min(.02,min(points)/2) if points else .02
                        upper=max(1e4,max(points)*2) if points else 1e4
                        ax.set_ylim(lower,upper);display_limits[metric]=[lower,upper]
                    audit=[]
                    for i,f in enumerate(inputs):
                        condition=sorted([r for r in rr if n(r,'f')==f],key=lambda r:n(r,'lag'))
                        usable=[r for r in condition if R.yes(r['usable'])]
                        color=plt.get_cmap('viridis')(i/max(1,len(inputs)-1))
                        forced=condition and R.yes(condition[0]['forced_combined'])
                        drift=condition and R.yes(condition[0]['non_exponential'])
                        label=f'{f:g}'+(' *' if drift else '')+(' [forced]' if forced else '')
                        for ax,metric in zip(axes,('rho_k','rho_geo')):
                            x=[n(r,'lag') for r in usable]
                            y=[n(r,metric) for r in usable]
                            ax.plot(x,y,'s--' if drift else 'o-',color='#999999' if forced else color,
                                    label=label,lw=1,ms=3)
                            for r in usable:
                                point=n(r,metric);lo=n(r,metric+'_ci_lower');hi=n(r,metric+'_ci_upper')
                                if np.isfinite(point) and point>0:
                                    bottom,top=display_limits[metric];lag=n(r,'lag')
                                    ax.plot([lag]*2,[max(lo,bottom) if np.isfinite(lo) else bottom,
                                            min(hi,top) if np.isfinite(hi) else top],color=color,alpha=.2,lw=.7)
                                    if not np.isfinite(hi) or hi>top:
                                        ax.plot(lag,top,'^',mfc='none',color=color,ms=4,clip_on=False)
                                        ax.annotate('',(lag,top),(lag,top/2),arrowprops=dict(arrowstyle='->',color=color,lw=.7))
                                    if not np.isfinite(lo) or lo<bottom:
                                        ax.plot(lag,bottom,'v',mfc='none',color=color,ms=4,clip_on=False)
                                        ax.annotate('',(lag,bottom),(lag,bottom*2),arrowprops=dict(arrowstyle='->',color=color,lw=.7))
                        audit.append(dict(f=f,usable_lags=[n(r,'lag') for r in usable],
                                          finite_ratios=sum(bool(np.isfinite(n(r,'rho_k'))) for r in usable),
                                          undefined_log_ratios=sum(not np.isfinite(n(r,'rho_k')) for r in usable),
                                          non_exponential=bool(drift),forced=bool(forced)))
                    for ax in axes:
                        ax.set(xscale='log',yscale='log',xlabel='Lag');ax.grid(alpha=.15)
                        ticks=sorted({n(r,'lag') for r in rr})
                        ax.set_xticks(ticks,[f'{x:g}' for x in ticks]);ax.xaxis.set_minor_formatter(NullFormatter())
                        ax.axhline(1,color='black',lw=.7)
                    axes[0].set_ylabel(r'Decay-rate ratio $\rho_k$')
                    axes[0].axhline(9,color='grey',ls=':',lw=.7)
                    axes[1].set_ylabel(r'$\rho_{\mathrm{geo}}$: relative to single-front geometry')
                    axes[1].legend(title='Input f (*: >1.5× lag drift)',fontsize=7,loc='upper left',bbox_to_anchor=(1.01,1))
                    fig.suptitle(f'{R.NAMES[model]}: all usable lags — {arm}, {view}',fontsize=13)
                    counts='; '.join(f"f={a['f']:g}: {len(a['usable_lags'])} usable, {a['undefined_log_ratios']} undefined" for a in audit)
                    import textwrap
                    fig.text(.07,.06,'\n'.join(textwrap.wrap(counts,145)),fontsize=8)
                    fig.text(.07,.015,'Pointwise whole-seed 95% intervals; every finite positive point is within the axes. Open edge triangles and arrows denote CI tails beyond the axes; full intervals are in CSV.\nA usable lag can have C ≤ 0, where −ln C and its ratio are undefined; all such counts are retained above. Asterisks/squares flag >1.5× drift. Geometry is not site diffusivity.',fontsize=8)
                    fig.subplots_adjust(left=.07,right=.81,bottom=.24,top=.90,wspace=.25)
                    pdf.savefig(fig);plt.close(fig)
                    pages.append(dict(arm=arm,view=view,model=model,conditions=audit,display_limits=display_limits,finite_points_clipped=0))
    dump(META/'all_lag_render_receipt.json',dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=R.sha(DATA/'all_lag_size_ratios.csv'),pages=pages,all_usable_lags_accounted=True,
        usable_lags=sum(R.yes(r['usable']) for r in rows),
        undefined_at_usable_lags=sum(R.yes(r['usable']) and not np.isfinite(n(r,'rho_k')) for r in rows)))

def main():
    ratios=R.rows(DATA/'size_ratios.csv');tols=R.rows(DATA/'noise_tolerance.csv');lags=R.rows(DATA/'condition_lag_metrics.csv')
    calls=amendment(ratios,tols);examples=native_examples()
    mainrat=[r for r in ratios if r['arm']=='common' and r['view']=='main' and r['model'] in MODELS]
    maintol=[r for r in tols if r['arm']=='common' and r['view']=='main' and r['model'] in MODELS and not R.yes(r['forced_inclusive'])]
    mainlag=[r for r in lags if r['arm']=='common' and r['view']=='main' and r['model'] in MODELS and n(r,'f')==.003 and n(r,'C_label') in (20,60,1,3)]
    metrics=R.rows(R.RUN/'results/final_analysis/condition_metrics.csv')
    mult=[r for r in metrics if r['arm']=='common' and r['view']=='main' and r['model'] in MODELS and n(r,'C_label') in ((1,3) if r['model']=='spring_ring' else (20,60))]
    for stem,rr in [('main_size_ratios',mainrat),('main_noise_tolerance',maintol),('main_autocorrelation',mainlag),('main_multiple_fronts',mult)]:R.write(DATA/(stem+'.csv'),rr)
    plt.rcParams.update({'font.size':9,'axes.titlesize':10,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(15,17));gs=fig.add_gridspec(4,1,height_ratios=[1.5,1,1.15,.9],hspace=.44)
    fig.suptitle('Enlarged cells hold heading longer under fluctuating input',fontsize=16,y=.99)
    fig.text(.5,.966,'Gains reach or exceed geometry at each model’s weakest measurable input; the advantage declines as input strengthens',ha='center',fontsize=13)
    # A: small panels have exactly one third the large ring's physical width.
    grid=gs[0].subgridspec(2,7,width_ratios=[1,3,.4,1,3,.4,3],height_ratios=[1,.055],hspace=.55,wspace=.35)
    for j,(model,i) in enumerate([('wave_pinning',0),('holmes_model3',2)]):
        pair=examples[i:i+2];vmax=max(e['field'].max() for e in pair)
        for col,e in enumerate(pair):
            ax=fig.add_subplot(grid[0,j*3+col]);ax.pcolormesh(np.linspace(0,e['circ'],e['field'].shape[1]+1),edges(e['time']),e['field'],shading='flat',cmap='viridis',vmin=0,vmax=vmax,rasterized=True)
            ax.set_title(f"{R.NAMES[model]} · {'small' if col==0 else 'large'}\nf = 0.003",fontsize=9);ax.set_xlabel('Arc position');ax.invert_yaxis()
            if col==0:ax.set_ylabel('Time')
            else:ax.tick_params(labelleft=False)
    e=examples[-1];ax=fig.add_subplot(grid[0,6]);ax.pcolormesh(np.linspace(0,e['circ'],e['field'].shape[1]+1),edges(e['time']),e['field'],shading='flat',cmap='viridis',vmin=0,vmax=max(examples[0]['field'].max(),examples[1]['field'].max(),e['field'].max()),rasterized=True);ax.set_title('Wave-pinning · large\nf = 0.1, multiple fronts');ax.invert_yaxis();ax.set_xlabel('Arc position')
    # Bottom row of A remains available for the explicitly shared colour scales.
    for j,i in enumerate((0,2)):
        cax=fig.add_subplot(grid[1,j*3:j*3+2]);pair=examples[i:i+2];vmax=max(e['field'].max() for e in pair+([examples[-1]] if i==0 else []))
        # Ensure all WP images, including strong input, share a model scale.
        for imageax in fig.axes:
            if imageax.collections and imageax.get_title().startswith(R.NAMES[pair[0]['model']]): imageax.collections[0].set_clim(0,vmax)
        gradient=np.linspace(0,vmax,256)[None,:];cax.imshow(gradient,extent=[0,vmax,0,1],aspect='auto',cmap='viridis');cax.set_yticks([]);cax.set_xlabel(R.NAMES[pair[0]['model']]+' activity (shared scale)');cax.set_box_aspect(.025)
    fig.text(.055,.943,'A  Native activity histories',weight='bold',fontsize=12)
    bg=gs[1].subgridspec(1,3,wspace=.3)
    for j,m in enumerate(MODELS):
        ax=fig.add_subplot(bg[j]);sizes=(1,3) if m=='spring_ring' else (20,60)
        for C,style,alpha in zip(sizes,('-',':'),(.5,1)):
            rr=sorted([r for r in mainlag if r['model']==m and n(r,'C_label')==C],key=lambda r:n(r,'lag'))
            x=np.array([n(r,'lag') for r in rr]);y=np.array([n(r,'C') for r in rr]);lo=np.array([n(r,'C_ci_lower') for r in rr]);hi=np.array([n(r,'C_ci_upper') for r in rr]);keep=y>0
            ax.plot(x[keep],y[keep],style,color=R.COLORS[m],alpha=alpha,label='small' if C==sizes[0] else 'large');ax.fill_between(x,np.maximum(lo,1e-4),hi,color=R.COLORS[m],alpha=.1)
        # A common ordinate would conceal resolved weak-input slopes.
        ax.set(yscale='log',xlabel='Lag',title=R.NAMES[m],ylim=({'wave_pinning':(.9,1.005),'holmes_model3':(.75,1.015),'spring_ring':(.01,1.05)}[m]));ax.grid(alpha=.15)
        if j==0:ax.set_ylabel('B  Heading autocorrelation');ax.legend(fontsize=8)
    cd=gs[2].subgridspec(2,1,hspace=.12,height_ratios=[1.15,1]);axc=fig.add_subplot(cd[0]);axd=fig.add_subplot(cd[1],sharex=axc)
    for m in MODELS:
        rr=sorted([r for r in mainrat if r['model']==m],key=lambda r:n(r,'f'));color=R.COLORS[m]
        axc.plot([n(r,'f') for r in rr],[n(r,'rho_k') for r in rr],color=color,lw=1,label=R.NAMES[m])
        for r in rr:
            col='#999999' if R.yes(r['forced_combined']) else color
            fill={'single_front':'white','multiple_fronts':col,'low_contrast':'#aaaaaa'}.get(r['large_regime'],'white')
            marker(axc,n(r,'f'),n(r,'rho_k'),n(r,'rho_k_ci_lower'),n(r,'rho_k_ci_upper'),col,'s' if R.yes(r['non_exponential']) else 'o',fill)
            if r['status']=='one_sided_lower_bound':
                y=n(r,'rho_k_lower');axc.plot(n(r,'f'),y,'^',mfc='none',color=col);axc.annotate('',(n(r,'f'),y*1.8),(n(r,'f'),y),arrowprops=dict(arrowstyle='->',color=col))
        for C,style,alpha in zip((1,3) if m=='spring_ring' else (20,60),('-',':'),(.5,1)):
            ss=sorted([r for r in mult if r['model']==m and n(r,'C_label')==C],key=lambda r:n(r,'f'));x=[n(r,'f') for r in ss]
            axd.plot(x,[n(r,'multiple_fronts_fraction') for r in ss],style,color=color,alpha=alpha,lw=1.4)
            axd.fill_between(x,[n(r,'multiple_fronts_ci_lower') for r in ss],[n(r,'multiple_fronts_ci_upper') for r in ss],color=color,alpha=.06)
    axc.axhline(1,color='black',lw=.8);axc.axhline(9,color='grey',ls='--',lw=.8);axc.set(xscale='log',yscale='log',ylabel=r'C  Decay-rate ratio $\rho_k$',ylim=(.3,100));axc.legend(ncol=3,fontsize=8,loc='upper right');axc.grid(alpha=.15);axc.tick_params(labelbottom=False)
    axd.set(xscale='log',xlabel='Input fraction f',ylabel='D  Multiple-front fraction',ylim=(-.02,1.02));axd.grid(alpha=.15)
    eg=gs[3].subgridspec(1,3,wspace=.3)
    for j,m in enumerate(MODELS):
        ax=fig.add_subplot(eg[j]);rr=sorted([r for r in maintol if r['model']==m],key=lambda r:n(r,'C_label'))
        x=[n(r,'C_label') for r in rr];y=[n(r,'fhalf') for r in rr];ax.plot(x,y,'o-',color=R.COLORS[m],ms=4)
        ax.fill_between(x,[n(r,'fhalf_ci_lower') for r in rr],[n(r,'fhalf_ci_upper') for r in rr],color=R.COLORS[m],alpha=.15)
        ax.set(title=R.NAMES[m],xlabel='Designed circumference',yscale='log');ax.grid(alpha=.15)
        if j==0:ax.set_ylabel(r'E  Noise tolerance $f_{1/2}$')
    fig.subplots_adjust(left=.09,right=.97,bottom=.075,top=.9)
    fig.text(.09,.018,'Small: lighter solid; large: darker dotted. C: open fill single front; filled multiple fronts; grey fill low contrast; square >1.5× lag drift.\nPointwise whole-seed 95% CIs. Horizontal references: ratio 1 and nominal geometry 9. D is descriptive; it does not establish a cause of erosion.\nH1 uses a post hoc f < large-cell f½ domain; frozen H1 remains inconclusive. Goryachev retained in supplement only.',fontsize=8)
    for ext in ('png','pdf','svg'):fig.savefig(OUT/('figure4_main.'+ext),dpi=180)
    plt.close(fig)
    legend=(MESSAGE+' The geometric statement describes model-specific weakest measurable point estimates: f = 0.003 for Wave-pinning and Holmes, and f = 0.0003–0.001 for Spring. Confidence intervals also admit sub-geometric gains in Holmes and Spring. At the shared displayed input f = 0.003, Spring is already below its actual geometric expectation (rate ratio 7.55 versus geometry 8.55), while its large-cell regime is still classified as single front. (A) Native activity kymographs for small and large Wave-pinning and Holmes at f = 0.003, plus a large Wave-pinning example at f = 0.1; width is proportional to circumference, with shared activity scales within each model. '
        '(B) Heading autocorrelation at f = 0.003 on model-specific logarithmic ordinates with whole-seed 95% bootstrap bands; intervals extending below zero continue beyond the positive plotting range. (C) Small/large rotational decay-rate ratio at the shortest usable lag satisfying the frozen lag rule, with pointwise 95% CIs. Reference lines show unity and nominal single-front geometry (9); actual-circumference normalization is reported in the supplement. Marker fill describes the large-cell regime; squares flag non-exponential lag dependence. '
        '(D) Multiple-front fractions with bootstrap intervals, reported descriptively: Holmes and Spring can fall below geometry while still single-front, so these data do not establish multiple fronts as the cause of erosion. (E) Fixed-reference-lag noise tolerance increases with size; smallest/largest CIs separate in all three shown models. Small cells are lighter solid, large cells darker dotted. Native units are used. '
        'Post hoc H1 amendment: evaluate gains and supported reversals only at inputs strictly below the point estimate of the largest cell’s own f½; strict CI and size noise-tolerance rules are unchanged. Both main-view observation arms pass the amended H1 overall, while frozen rate H1 is inconclusive in both; all calls appear alongside in S0. H2/H3/A6 τ gates are retired. Goryachev is excluded from the main because its ratios derive from forced or non-exponential conditions, and remains complete in the supplement. Simulated polarity heading is distinct from migration persistence and conditional OU estimates.')
    (META/'main_legend.md').write_text(legend+'\n')
    maincalls=[r for r in calls if r['view']=='main' and not r['forced_inclusive'] and r['model']=='all']
    dump(META/'main_summary.json',dict(message=MESSAGE,case='A',post_hoc=True,H1_main=maincalls,retired_gates=['H2','H3','A6'],legend=legend,generated_utc=datetime.now(timezone.utc).isoformat()))
    assert len(maincalls)==2 and all(r['status']=='pass' for r in maincalls),maincalls
    lag_consistency_supplement()
    # All existing supplement material remains, including superseded main artwork.
    old=R.ROOT/'figures/length_mechanisms_20261008/dev/signature_figure4_20261020'
    parts=[OUT/'figure4_S0_posthoc.pdf',OUT/'figure4_S0_round2_old_new.pdf',OUT/'figure4_all_lag_ratios.pdf',VEC/'figure4_reanalysis_supplement.pdf',VEC/'figure4_old_new.pdf',VEC/'legacy_seven_view_robustness.pdf']
    parts += [old/'supporting/vectors/figure4_revision_main.pdf',old/'supporting/vectors/figure4_revision_complete_supplement.pdf']
    writer=PdfWriter()
    for p in parts:writer.append(p)
    with (OUT/'figure4_supplement.pdf').open('wb') as f:writer.write(f)
    dump(META/'main_render_receipt.json',dict(generated_utc=datetime.now(timezone.utc).isoformat(),no_simulations=True,post_hoc=True,main_panels=5,main_models=list(MODELS),palette=R.COLORS,inputs={str(p):R.sha(p) for p in [DATA/'size_ratios.csv',DATA/'noise_tolerance.csv',DATA/'condition_lag_metrics.csv']},outputs={str(p.relative_to(OUT)):R.sha(p) for p in [*OUT.glob('figure4_main.*'),OUT/'figure4_supplement.pdf']},supplement_parts=[str(p) for p in parts]))
    (OUT/'README.md').write_text('# Figure 4 — round 2, author-authorized Case A\n\n'+MESSAGE+'\n\n[Main PNG](figure4_main.png) · [PDF](figure4_main.pdf) · [SVG](figure4_main.svg) · [Supplement](figure4_supplement.pdf) · [Plotted CSVs](supporting/data/)\n\nSince 20 October: replaces τ with decay-rate ratios, adds absolute fixed-lag noise tolerance and descriptive geometry/multiple-front outcomes, composes Case A, and removes Goryachev from the main. The H1 domain amendment is **post hoc**; frozen and amended calls appear alongside in S0, for both arms, all views and both forcing choices. H2/H3/A6 gates are retired. Frozen interim tables remain unchanged. Non-exponential conditions, unavailable lags and censoring remain explicit; open triangles and arrows denote bounds.\n\n'+legend+'\n')
    (OUT/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Figure 4 round 2</title><style>body{font:17px/1.5 system-ui;max-width:1400px;margin:auto}img{width:100%}</style><h1>Figure 4: heading persistence under fluctuating input</h1><p>Generated '+datetime.now(timezone.utc).isoformat()+'. Analysis and rendering complete.</p><p>'+MESSAGE+'</p><p><a href="figure4_main.pdf">Main PDF</a> · <a href="figure4_main.svg">SVG</a> · <a href="figure4_supplement.pdf">Complete round-2 supplement</a> · <a href="supporting/data/">CSVs</a></p><img src="figure4_main.png"><p>'+legend+'</p>')
    print(json.dumps(dict(case='A',H1_main=[(r['arm'],r['frozen_rate_status'],r['status']) for r in maincalls])))

if __name__=='__main__':main()
