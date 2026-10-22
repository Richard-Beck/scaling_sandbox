"""Requested Figure 4 diagnostics, using actual compact estimator contracts."""
from __future__ import annotations
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from figure4_round3_render_20261022 import (MODELS,NAMES,SIZES,primary,read_csv,
    write_csv,write_json,n,finite,equivalent_size,size_color,style_axis,log_f,
    facet,alpha_panel,NULLS,NULL_NAMES,first)


def main_rows(rows):
    return [r for r in primary(rows) if n(r,'flicker_frames',default=5)==5 and math.isclose(n(r,'MSD_range_fraction',default=.25),.25)]


def smooth_fit(data):
    """Fixed degree-two tensor surface: a descriptive fit, never a main gate."""
    rows=primary(read_csv(data/'analysis_lag_metrics.csv'));curves=[];receipt=[]
    for model in MODELS:
        rr=[r for r in rows if r.get('model')==model and n(r,'k')>0 and n(r,'C_mean')>0]
        if len(rr)<9:
            receipt.append(dict(model=model,status='unavailable_fewer_than_nine_positive_conditions'));continue
        x=np.log([n(r,'C_mean') for r in rr]);y=np.log([n(r,'f') for r in rr]);z=np.log([n(r,'k') for r in rr])
        xm,ym=x.mean(),y.mean();xs,ys=x.std(),y.std()
        x=(x-xm)/xs;y=(y-ym)/ys
        powers=[(i,j) for i in range(3) for j in range(3)]
        design=np.column_stack([x**i*y**j for i,j in powers]);beta,residuals,rank,_=np.linalg.lstsq(design,z,rcond=None)
        if rank<len(powers):
            receipt.append(dict(model=model,status='unavailable_rank_deficient',rank=int(rank)));continue
        for size in SIZES:
            actual=[n(r,'C_mean') for r in rr if math.isclose(equivalent_size(r),size)]
            if not actual:continue
            C=float(np.mean(actual));xx=(math.log(C)-xm)/xs
            for f in np.geomspace(min(n(r,'f') for r in rr),max(n(r,'f') for r in rr),81):
                yy=(math.log(f)-ym)/ys
                derivative=sum(b*i*xx**(i-1)*yy**j/xs for b,(i,j) in zip(beta,powers) if i)
                curves.append(dict(model=model,C_label=size/20 if model=='spring_ring' else size,C_mean=C,f=f,alpha=-float(derivative),status='diagnostic_fixed_degree2_tensor_fit',source_conditions=len(rr)))
        receipt.append(dict(model=model,status='resolved',fit='unweighted least squares tensor powers0..2 in ln(Cmean),ln(f); nine fixed coefficients; no tuning',source_conditions=len(rr),coefficients=beta.tolist(),lnC_center=xm,lnf_center=ym,lnC_scale=xs,lnf_scale=ys,residual_sum_squares=float(np.sum((design@beta-z)**2))))
    write_csv(data/'smooth_alpha.csv',curves)
    write_json(data.parent/'metadata/smooth_fit_receipt.json',receipt)
    return curves


def supplement(data, output):
    plots=[];main=main_rows(read_csv(data/'drift_decomposition.csv'))
    all_tracks=primary(read_csv(data/'drift_decomposition.csv'))
    detector=primary(read_csv(data/'detector_free_mle.csv'))
    with PdfPages(output) as pdf:
        def finish(fig,title):
            fig.suptitle(title,fontsize=10,y=.98)
            fig.subplots_adjust(left=.1,right=.97,bottom=.11,top=.90,wspace=.42,hspace=.47)
            pdf.savefig(fig);plt.close(fig);plots.append(title)
        def notice(title,body,height=6.):
            fig,ax=plt.subplots(figsize=(7.0866,height))
            ax.axis('off');ax.text(.02,.9,body,transform=ax.transAxes,va='top',fontsize=8,wrap=True,linespacing=1.6)
            finish(fig,title)
        def grid():return plt.subplots(2,2,figsize=(7.0866,6.8))
        def trend(rows,metric,title,*,logy=False,limits=None):
            fig,axes=grid()
            for ax,model in zip(axes.flat,MODELS):
                for size in SIZES:
                    rr=sorted([r for r in rows if r.get('model')==model and math.isclose(equivalent_size(r),size)],key=lambda r:n(r,'f'))
                    yy=[n(r,metric) if (n(r,metric)>0 if logy else finite(n(r,metric))) else math.nan for r in rr]
                    ax.plot([n(r,'f') for r in rr],yy,'o-',color=size_color(model,size),ms=2,lw=.65)
                log_f(ax);facet(ax,NAMES[model]);ax.set_ylabel(metric.replace('_',' '));style_axis(ax)
                if logy:ax.set_yscale('log')
                if limits:ax.set_ylim(limits)
            finish(fig,title)

        smooth=smooth_fit(data)
        if smooth:
            fig,axes=grid()
            for ax,model in zip(axes.flat,MODELS):alpha_panel(ax,smooth,model,ylabel=True)
            fig.text(.1,.025,'Fixed degree-two tensor polynomial in ln k versus ln circumference and ln input; diagnostic only. No missing local estimate is filled in the main figure.',fontsize=6,wrap=True)
            finish(fig,'Smooth diagnostic size–input surface')
        else:notice('Smooth diagnostic size–input surface','Insufficient positive rate conditions for the fixed nine-coefficient tensor fit; no curve is imputed.')

        key=lambda r:(r.get('model'),n(r,'C_label'),n(r,'f'))
        tracked={key(r):r for r in main}
        for metric in ['Ds','lambda_eff']:
            matched=[]
            for d in detector:
                if n(d,'detector_free_step',default=1)!=1:continue
                t=tracked.get(key(d))
                if t and n(t,metric)>0 and n(d,metric)>0:matched.append((t,d))
            if not matched:
                notice('Detector-free versus tracked '+metric,'No positive matched estimates are available. Zero switch contributions are empirical zeros and cannot appear on logarithmic axes. Undefined conditions remain in the source records.');continue
            fig,axes=grid()
            for ax,model in zip(axes.flat,MODELS):
                points=[(t,d) for t,d in matched if t['model']==model]
                for t,d in points:ax.plot(n(t,metric),n(d,metric),'o',color=size_color(model,equivalent_size(t)),ms=3)
                vals=[n(r,metric) for pair in points for r in pair]
                limits=(min(vals)/2,max(vals)*2) if vals else (1e-8,1.)
                ax.plot(limits,limits,color='.65',lw=.6);ax.set(xscale='log',yscale='log',xlim=limits,ylim=limits,xlabel='Tracked '+metric,ylabel='Detector-free '+metric)
                facet(ax,NAMES[model]);style_axis(ax)
            finish(fig,'Detector-free versus tracked '+metric+'; five circumference shades')
            write_csv(data/('supplement_scatter_'+metric+'.csv'),[dict(model=t['model'],C_label=t['C_label'],f=t['f'],tracked=n(t,metric),detector_free=n(d,metric)) for t,d in matched])
        for metric in ['Ds','lambda_eff']:
            for step in [1,2]:
                rr=[r for r in detector if n(r,'detector_free_step',default=1)==step]
                if rr:trend(rr,metric,f'Detector-free {metric}; cadence multiplier {step}',logy=True)

        if main and all_tracks:
            for metric in ['Ds','lambda_eff']:
                fig,axes=grid()
                ratios=[]
                for ax,model in zip(axes.flat,MODELS):
                    for size in SIZES:
                        for hold,fraction,ls in [(2,.25,'--'),(10,.25,':'),(5,.125,'-.'),(5,.5,'-')]:
                            rr=sorted([r for r in all_tracks if r.get('model')==model and math.isclose(equivalent_size(r),size) and n(r,'flicker_frames')==hold and math.isclose(n(r,'MSD_range_fraction'),fraction)],key=lambda r:n(r,'f'))
                            xx=[];yy=[]
                            for r in rr:
                                baseline=tracked.get(key(r),{});value=n(r,metric)/n(baseline,metric) if n(baseline,metric)>0 else math.nan
                                xx.append(n(r,'f'));yy.append(value if value>0 else math.nan)
                                ratios.append(dict(model=model,C_label=r['C_label'],f=r['f'],metric=metric,hold=hold,MSD_fraction=fraction,ratio=value))
                            ax.plot(xx,yy,ls,color=size_color(model,size),lw=.6)
                    ax.axhline(1,color='.65',lw=.5);log_f(ax);ax.set_yscale('log');ax.set_ylabel('Sensitivity / primary');facet(ax,NAMES[model]);style_axis(ax)
                fig.text(.1,.025,'Dashed: hold 2; dotted: hold 10; dash-dot: MSD range 1/8; solid: MSD range 1/2. Primary hold 5, range 1/4.',fontsize=6)
                finish(fig,'Tracked sensitivity: '+metric)
                write_csv(data/('supplement_sensitivity_'+metric+'.csv'),ratios)
        else:notice('Tracked sensitivity','The declared 2/10-frame persistence and 1/8,1/2-tenure MSD sensitivities await tracked compact records.')

        valid=read_csv(data/'valid_fraction.csv')
        selected_rates={key(r):r for r in primary(read_csv(data/'analysis_lag_metrics.csv'))}
        if not valid:
            condition=primary(read_csv(data/'analysis_lag_metrics.csv'))
            valid=[dict(model=r['model'],C_label=r['C_label'],f=r['f'],valid_fraction=n(r,'valid_fraction'),bootstrap_valid_fraction=n(r,'bootstrap_valid_fraction'),status=r.get('status')) for r in condition]
            write_csv(data/'valid_fraction.csv',valid)
        for row in valid:
            selected=selected_rates.get(key(row),{})
            row['bootstrap_valid_fraction']=selected.get('bootstrap_valid_fraction')
            row['analysis_status']=selected.get('status')
            row['analysis_lag']=selected.get('analysis_lag')
        write_csv(data/'supplement_valid_fraction.csv',valid)
        for metric in ['valid_fraction','bootstrap_valid_fraction']:
            trend(valid,metric,'Every-condition '+metric.replace('_',' '),limits=(-.035,1.035))
        if main:
            trend(main,'MSD_exponent','Within-tenure free MSD power exponent')
            msds=read_csv(data/'within_tenure_msd.csv')
            for model in MODELS:
                rr=[r for r in msds if r.get('model')==model]
                if not rr:continue
                fig,axes=plt.subplots(1,5,figsize=(10.5,3.8))
                fs=sorted({n(r,'f') for r in rr})
                for ax,size in zip(axes,SIZES):
                    for j,f in enumerate(fs):
                        q=sorted([r for r in rr if math.isclose(equivalent_size(r),size) and n(r,'f')==f],key=lambda r:n(r,'lag'))
                        ax.plot([n(r,'lag') for r in q],[n(r,'MSD') for r in q],color=plt.get_cmap('viridis')(j/max(1,len(fs)-1)),lw=.55)
                    ax.set(xscale='log',yscale='log',xlabel='Lag',ylabel='Physical MSD',title=f'C={size/20 if model=="spring_ring" else size:g}');style_axis(ax)
                finish(fig,NAMES[model]+': within-tenure MSD; increasing input from purple to yellow')
        else:notice('Within-tenure MSD exponents','No tracked physical-centroid tenures support MSD fits in these records.')

        comparison=read_csv(data/'round2_comparison.csv')
        if comparison:
            if any('measure' in r for r in comparison):
                for measure,label in [('k_at_fixed_reference','Fixed-lag decay rate'),('fhalf_first_downcrossing','Noise tolerance')]:
                    selected=[r for r in comparison if r.get('measure')==measure]
                    fig,axes=grid()
                    for ax,model in zip(axes.flat,MODELS):
                        q=[r for r in selected if r.get('model')==model and n(r,'matched_P_round2_value')>0 and n(r,'round3_first12_value')>0]
                        for r in q:ax.plot(n(r,'matched_P_round2_value'),n(r,'round3_first12_value'),'o',ms=3,color=size_color(model,equivalent_size(r)))
                        vals=[n(r,k) for r in q for k in ('matched_P_round2_value','round3_first12_value')];limits=(min(vals)/2,max(vals)*2) if vals else (1e-5,1.)
                        ax.plot(limits,limits,color='.6',lw=.6);ax.set(xscale='log',yscale='log',xlim=limits,ylim=limits,xlabel='Round 2, matched P definition',ylabel='Round 3, first 12 histories');facet(ax,NAMES[model]);style_axis(ax)
                    fig.text(.1,.025,'Rate comparisons use the actual old 12-history subset. Tolerance crosses an input grid whose old replication can vary; all seed counts and original/matched-P baselines remain in CSV.',fontsize=6,wrap=True)
                    finish(fig,'Round-2 versus round-3: '+label)
            # Find explicit old/new columns while preserving the full CSV.
            for old,new,label in [('round2_k','round3_k','Round-2 versus round-3 k'),('round2_fhalf','round3_fhalf','Round-2 versus round-3 noise tolerance')]:
                if not any(n(r,old)>0 and n(r,new)>0 for r in comparison):continue
                fig,axes=grid()
                for ax,model in zip(axes.flat,MODELS):
                    q=[r for r in comparison if r.get('model')==model and n(r,old)>0 and n(r,new)>0]
                    for r in q:ax.plot(n(r,old),n(r,new),'o',ms=3,color=size_color(model,equivalent_size(r)))
                    vals=[n(r,k) for r in q for k in (old,new)];limits=(min(vals)/2,max(vals)*2) if vals else (1e-5,1.)
                    ax.plot(limits,limits,color='.6',lw=.6);ax.set(xscale='log',yscale='log',xlim=limits,ylim=limits,xlabel='Round 2',ylabel='Round 3');facet(ax,NAMES[model]);style_axis(ax)
                fig.text(.1,.025,'Actual old/new independent-history counts are retained in the compact comparison CSV; sparse old ensembles are not relabelled as 12-seed ensembles.',fontsize=6,wrap=True)
                finish(fig,label+' comparison')
        else:notice('Round-2 versus round-3 12-seed comparison','The new from-scratch compact campaign is pending. A reanalysis of the same retained records is not an independent campaign replication; the historical round-2 evidence is preserved.')
        extended=[r for r in read_csv(data/'extended_conditions.csv') if r.get('status')!='none']
        receipt_path=data.parent/'metadata/analysis_receipt.json'
        import json
        completed=json.loads(receipt_path.read_text()).get('completed_campaign',False) if receipt_path.exists() else False
        absent='The completed compact campaign has no extended horizon or condition above 96 histories.' if completed else 'No completed new compact condition, horizon extension or condition above 96 seeds is represented by this provisional figure.'
        if extended:
            lines=[f"{r.get('model')}, C={r.get('C_label')}, f={r.get('f')}: n={r.get('n_seeds')}, horizon={first(r,'max_horizon','horizon')}" for r in extended]
            for start in range(0,len(lines),35):notice(f'Extended horizons and conditions above 96 seeds ({start+1}–{min(start+35,len(lines))})','\n'.join(lines[start:start+35]),height=8.)
        else:notice('Extended horizons and conditions above 96 seeds',absent)

        checks=read_csv(data/'null_B_analytic_agreement.csv')
        if checks:
            fig,axes=plt.subplots(1,3,figsize=(8.,3.8))
            for ax,kind in zip(axes,NULLS):
                rr=[r for r in checks if r.get('model')==kind]
                for size in SIZES:
                    q=sorted([r for r in rr if math.isclose(equivalent_size(r),size)],key=lambda r:n(r,'f'))
                    ax.plot([n(r,'f') for r in q],[n(r,'relative_deviation') for r in q],'o-',color=size_color('null',size),ms=1.8,lw=.65)
                ax.axhline(0,color='.7',lw=.5);log_f(ax);facet(ax,kind.replace('_',' '));ax.set_ylabel('(Simulated / analytic k) − 1');style_axis(ax)
            finish(fig,'Null B analytic agreement; analytic rates are never substituted into main curves')
        acf=read_csv(data/'null_A_acf.csv')
        if acf:
            fig,axes=plt.subplots(1,3,figsize=(8.,3.8))
            for ax,kind in zip(axes,NULLS):
                rr=[r for r in acf if r.get('model')==kind]
                for size in (20,60):
                    q=sorted([r for r in rr if n(r,'C_label')==size],key=lambda r:n(r,'lag'))
                    ax.plot([n(r,'lag')/1000 for r in q],[n(r,'C')-n(r,'analytic_C') for r in q],color=size_color('null',size),lw=.8)
                ax.axhline(0,color='.7',lw=.5);ax.set(xlabel='Lag (×10³)',ylabel='Simulated − analytic C',title=kind.replace('_',' '));style_axis(ax)
            finish(fig,'Null A analytic autocorrelation agreement')
        if not checks:notice('Null B analytic agreement','Null B simulations have not yet completed. No analytic replacement is plotted.')

        alphas=primary(read_csv(data/'local_alpha.csv'))
        # A graphical census preserves all 185/260 designed-condition statuses.
        fig,axes=grid()
        for ax,model in zip(axes.flat,MODELS):
            q=[r for r in alphas if r.get('model')==model]
            for r in q:ax.plot(n(r,'f'),equivalent_size(r),'o' if finite(n(r,'alpha')) else 'x',ms=3,color=size_color(model,equivalent_size(r)))
            log_f(ax);ax.set_ylabel('C (Spring ×20)');facet(ax,NAMES[model]);style_axis(ax)
        finish(fig,'Every local neighborhood: circle resolved, cross unavailable; reasons in compact CSV')
    write_json(data.parent/'metadata/supplement_receipt.json',dict(pages=len(plots),diagnostics=plots,no_tables=True,undefined_not_imputed=True))
