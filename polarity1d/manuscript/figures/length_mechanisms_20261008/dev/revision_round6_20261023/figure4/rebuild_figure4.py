#!/home/4473331/.conda/envs/scaling_sandbox/bin/python
"""Single rebuild entry point. Read source evidence, collect finished jobs, render.

No model integrations or scheduler submissions occur here. Missing rerun evidence
is explicitly recorded; re-running this entry point incorporates newly finished
condition analyses and mode-only precision outputs.
"""
import os,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'scripts'))
from common import *
import math,gzip,shutil,subprocess
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import figure4_round3_render_20261022 as R
import figure4_round3_estimators_20261022 as E
import exemplars
from mode_analysis import model_test

GENERATED=datetime.now(timezone.utc).isoformat()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7.3,'axes.titlesize':8,
   'xtick.labelsize':6.3,'ytick.labelsize':6.3,'legend.fontsize':6.2,'axes.linewidth':.5,'pdf.fonttype':42})

def primary(rows):return R.primary(rows)
def sizes(model):return (1.,1.5,2.,2.5,3.) if model=='spring_ring' else (20.,30.,40.,50.,60.)
def color(model,size):return R.size_color(model,size*20 if model=='spring_ring' else size)
def readold(name):return csvread(OLD/'supporting/data'/f'{name}.csv')
def forcing_flag(model,size,f):
    return any(r.get('forced_combined')=='True' for r in readold('forcing_conditions') if r['model']==model and num(r,'C_label')==size and num(r,'f')==f)
def style(ax,model=None,logf=False):
    R.style_axis(ax)
    ax.tick_params(labelsize=6.3)
    if model:ax.set_title(NAMES[model],color='black',fontsize=8,pad=3)
    if logf:R.log_f(ax)
def label(fig,spec,letter):
    p=spec.get_position(fig);fig.text(max(.014,p.x0-.055),p.y1+.005,letter,fontsize=11,fontweight='bold',color='black',ha='left',va='bottom')
def panel_axes(fig,spec,models=MODELS):
    g=spec.subgridspec(1,len(models),wspace=.36)
    return [fig.add_subplot(g[j]) for j in range(len(models))]
def lines(ax,rows,model,key,ci=False,logf=True):
    for size in sizes(model):
        rr=sorted([r for r in rows if r['model']==model and num(r,'C_label')==size],key=lambda r:num(r,'f'))
        x=[num(r,'f') for r in rr];y=[num(r,key) for r in rr];co=color(model,size)
        ax.plot(x,y,color=co,lw=.85,marker='o',ms=1.7)
        if ci:
            for r,xx,yy in zip(rr,x,y):
                lo=num(r,key+'_ci_lower');hi=num(r,key+'_ci_upper')
                if np.isfinite(yy) and np.isfinite(lo+hi):ax.plot([xx,xx],[lo,hi],color=co,lw=.4)
    style(ax,model,logf)

def stage_sources():
    DATA.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    names=['local_alpha','noise_tolerance_folds','noise_tolerance','heading_amplitude_condition_metrics','birth_heading_condition_metrics',
      'takeover_by_condition_distance','takeover_by_size_distance','null_A_acf','model_D_selection','condition_lag_metrics','analysis_lag_metrics','jump_scoring_coverage']
    for name in names:
        shutil.copyfile(OLD/'supporting/data'/f'{name}.csv',DATA/f'{name}.csv')
    rows=[r for r in primary(readold('drift_decomposition')) if num(r,'flicker_frames')==5 and num(r,'MSD_range_fraction')==.25]
    coverage={(r['model'],num(r,'C_label'),num(r,'f')):r for r in readold('jump_scoring_coverage')}
    closure=[]
    for r in rows:
        cov=coverage[r['model'],num(r,'C_label'),num(r,'f')];k=num(r,'kobs');q=2*np.pi/num(r,'C_mean');drift=q*q*num(r,'Ds')
        finite_drift=np.isfinite(drift)
        lo=(drift if finite_drift else 0.)+num(cov,'lambda_eff_missing_angle_lower')
        hi=drift+num(cov,'lambda_eff_missing_angle_upper') if finite_drift else np.inf
        availability='resolved_interval' if finite_drift else 'unbounded_missing_Ds'
        if not np.isfinite(k) or k<=0:availability+=';unavailable_observed_rate'
        closure.append(dict(r,kpred_lower=lo,kpred_upper=hi,ratio_lower=lo/k if k>0 else np.nan,ratio_upper=hi/k if k>0 else np.nan,
            n_unscored=int(num(cov,'n_switches')-num(cov,'n_scored_switches')),closure_availability=availability,
            interval='all unscored angular losses independently bounded in [0,2]; missing Ds has nonnegative lower bound and unbounded upper limit'))
    csvwrite(DATA/'G_closure_intervals.csv',closure);csvwrite(DATA/'F_drift_ratio.csv',rows)
    # I selection is per model. Gate on all certified births at every size,
    # exactly as requested; remove d<=0 only when computing plotted probabilities.
    take=readold('takeover_by_condition_distance');tol=readold('noise_tolerance');selected=[];out=[]
    allseed=readold('takeover_seed_counts')
    for model in MODELS:
        ref=200. if model=='spring_ring' else 2000.
        tolerance=next(r for r in tol if r['model']==model and num(r,'C_label')==sizes(model)[0] and num(r,'reference_lag')==ref)
        half=num(tolerance,'fhalf');fs=sorted({num(r,'f') for r in take if r['model']==model});eligible=[]
        for f in fs:
            counts=[num(next(r for r in readold('birth_heading_condition_metrics') if r['model']==model and num(r,'C_label')==s and num(r,'f')==f),'certified_births') for s in sizes(model)]
            if min(counts)>=200:eligible.append((abs(np.log(f/half)),f,min(counts)))
        if not eligible:
            selected.append(dict(model=model,status='unavailable_no_input_with_200_births_at_all_sizes',smallest_fhalf=half));continue
        _,f,minbirths=min(eligible)
        usable=[sum(num(r,'births') for r in take if r['model']==model and num(r,'C_label')==s and num(r,'f')==f and num(r,'distance_bin')>0) for s in sizes(model)]
        selected.append(dict(model=model,f=f,smallest_fhalf=half,minimum_births=minbirths,minimum_usable_positive_distance_births=min(usable),reference_lag=ref,
              rule='nearest designed f in log distance to smallest fhalf; >=200 total certified births at all five sizes; d<=0 excluded only from plotted probabilities'))
        for size in sizes(model):
            for binid in (1,2):
                r=dict(next(r for r in take if r['model']==model and num(r,'C_label')==size and num(r,'f')==f and num(r,'distance_bin')==binid))
                seeds=[s for s in allseed if s['model']==model and num(s,'C_label')==size and num(s,'f')==f and num(s,'distance_bin')==binid]
                births=np.array([num(s,'births') for s in seeds]);wins=np.array([num(s,'takeovers') for s in seeds]);n=len(seeds)
                w=np.random.default_rng(E.stable_seed(('round6_takeover',model,size,f,binid))).multinomial(n,np.full(n,1/n),2000)
                with np.errstate(divide='ignore',invalid='ignore'):boot=(w@wins)/(w@births)
                valid=np.isfinite(boot);r['bootstrap_valid_fraction']=float(valid.mean())
                if valid.mean()<.975:
                    r['ci_lower'],r['ci_upper']=0.,1.;r['interval']='uninformative probability bounds [0,1]; bootstrap defined in fewer than97.5% draws'
                else:
                    r['ci_lower'],r['ci_upper']=E.percentile(boot,True);r['interval']='whole independent histories, 2000 bootstrap draws'
                out.append(r)
    csvwrite(DATA/'I_takeover_selection.csv',selected);csvwrite(DATA/'I_takeover_points.csv',out)
    return closure

def kymographs(fig,spec,null=False,export=True):
    prefix='null_A' if null else 'model_D';records=R.load_kymographs(OLD/'supporting/data',prefix)
    kinds=R.NULLS if null else MODELS;grid=spec.subgridspec(1,len(kinds),wspace=.38)
    index=[]
    for j,model in enumerate(kinds):
        pair=grid[j].subgridspec(1,2,width_ratios=[1,3],wspace=.18)
        selected=[r for r in records if (R.null_kind(r)==model if null else r['model']==model)]
        vmax=1 if null else max(r['field'].max() for r in selected)
        for k,size in enumerate((20,60)):
            r=next(r for r in selected if R.equivalent_size(r)==size);field=r['field'].copy();t=r['time']
            if not null:
                nb=field.shape[1];theta=2*np.pi*(np.arange(nb)+.5)/nb
                phase=float(np.angle(field[0]@np.exp(1j*theta))%(2*np.pi));shift=int(np.rint(nb*(.5-phase/(2*np.pi))))
                field=np.roll(field,shift,axis=1)
            # Rendering uses all saved bins and frames, exactly as the source.
            row=dict(r,field=field);ax=fig.add_subplot(pair[k]);native_size=size/20 if model=='spring_ring' else size
            vv=(9.2 if k==0 else 83.) if model=='goryachev' else vmax
            R.kymo(ax,row,vv,ylabel=k==0,label=f'C = {native_size:g}',elapsed=not null)
            ax.tick_params(labelsize=5.8)
            if not null and model=='goryachev':ax.text(.5,.96,f'0–{vv:g}',transform=ax.transAxes,ha='center',va='top',fontsize=5.8,bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.6))
            filename=f'{"A" if null else "C"}_kymograph_{model}_C{native_size:g}.csv.gz'
            index.append(dict(panel='A' if null else 'C',model=model,C_label=native_size,f=r.get('f',''),file=filename,
                field_shape=str(list(field.shape)),vmin=0,vmax=vv,source=str(r['npz_path']),
                plotted_time_origin=0 if null else float(t[0]),plotted_time_divisor=1000.,
                plotted_time_formula='(CSV time - plotted_time_origin)/plotted_time_divisor'))
            if export and not (DATA/filename).exists():
                with gzip.open(DATA/filename,'wt',newline='') as handle:
                    writer=csv.writer(handle);writer.writerow(['time','arc_position','activity'])
                    x=.5*(r['space_edges'][:-1]+r['space_edges'][1:])
                    for ti,values in zip(t,field):writer.writerows((format(ti,'.9g'),format(xi,'.7g'),format(v,'.6g')) for xi,v in zip(x,values))
        pos=grid[j].get_position(fig)
        title=R.NULL_NAMES[model].replace('\n',' ') if null else NAMES[model]+f', f = {selected[0]["f"]}'
        fig.text((pos.x0+pos.x1)/2,pos.y1+.005,title,ha='center',fontsize=7.6,color='black')
    if export:csvwrite(DATA/('A_kymograph_index.csv' if null else 'C_kymograph_index.csv'),index)

def closure_plot(axes,rows):
    for ax,model in zip(axes,MODELS):
        unbounded=[];unavailable=[]
        for size in sizes(model):
            rr=sorted([r for r in rows if r['model']==model and num(r,'C_label')==size],key=lambda r:num(r,'f'))
            x=np.array([num(r,'f') for r in rr]);lo=np.array([num(r,'ratio_lower') for r in rr]);hi=np.array([num(r,'ratio_upper') for r in rr]);co=color(model,size)
            ax.fill_between(x,lo,hi,color=co,alpha=.1,lw=0)
            ax.plot(x,lo,color=co,lw=.7);ax.plot(x,hi,color=co,lw=.6,ls='--')
            ax.vlines(x,lo,hi,color=co,lw=.45)
            unbounded.extend((num(r,'f'),co) for r in rr if np.isinf(num(r,'ratio_upper')))
            unavailable.extend((num(r,'f'),co) for r in rr if not np.isfinite(num(r,'kobs')))
        ax.axhline(1,color='.55',lw=.65);ax.set_yscale('log');style(ax,model,True)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ymin,ymax=ax.get_ylim()
        for x,co in unbounded:
            ax.plot(x,ymax/1.2,'^',mfc='white',mec=co,ms=3,mew=.65)
            ax.annotate('',xy=(x,ymax),xytext=(x,ymax/1.2),arrowprops=dict(arrowstyle='->',lw=.5,color=co))
        for x,co in unavailable:ax.plot(x,.965,'x',transform=ax.get_xaxis_transform(),color=co,ms=2.4,mew=.6)
        if unavailable:ax.text(.02,.97,'× observed rate unavailable',transform=ax.transAxes,va='top',fontsize=5.8,color='.35')
    axes[0].set_ylabel(r'$k_{\rm pred}/k_{\rm obs}$')

def takeover_plot(axes,rows):
    selection=csvread(DATA/'I_takeover_selection.csv')
    for ax,model in zip(axes,MODELS):
        for binid,ls,mk in ((1,'-','o'),(2,'--','s')):
            rr=sorted([r for r in rows if r['model']==model and num(r,'distance_bin')==binid],key=lambda r:num(r,'C_label'))
            x=[num(r,'C_label') for r in rr];y=[num(r,'probability') for r in rr]
            ax.plot(x,y,color=R.COLORS[model],ls=ls,lw=.9,marker=mk,ms=3,mfc='white' if binid==2 else R.COLORS[model],mew=.7)
            for r,xx,yy in zip(rr,x,y):
                lo,hi=num(r,'ci_lower'),num(r,'ci_upper')
                if np.isfinite(lo+hi+yy):ax.plot([xx,xx],[lo,hi],color=R.COLORS[model],lw=.6)
        chosen=next((r for r in selection if r['model']==model),{})
        style(ax);ax.set_title(NAMES[model]+f', f = {num(chosen,"f"):g}',fontsize=7.8,pad=3)
        ax.set_ylim(-.02,1.02);ax.set_xticks(sizes(model)[::2]);ax.set_xlabel('Circumference C')
    axes[0].set_ylabel('Takeover probability (20τ)')
    axes[0].legend(handles=[Line2D([0],[0],color='.25',marker='o',ms=3,lw=.8,label='0 < d ≤ C/6'),Line2D([0],[0],color='.25',marker='s',mfc='white',ls='--',ms=3,lw=.8,label='d > C/6')],frameon=False,fontsize=6,loc='upper right')

def main_figure(closure):
    fig=plt.figure(figsize=(12,16))
    grid=fig.add_gridspec(8,1,left=.075,right=.975,bottom=.045,top=.975,height_ratios=[1.1,1.05,1.4,1.2,1.2,1.2,1.2,1.2],hspace=.62)
    label(fig,grid[0],'A');kymographs(fig,grid[0],True)
    label(fig,grid[1],'B');axes=panel_axes(fig,grid[1],R.NULLS)
    acf=readold('null_A_acf')
    for j,(ax,kind) in enumerate(zip(axes,R.NULLS)):
        R.correlation(ax,acf,kind,null=True,legend=j==0);ax.text(.97,.88,r'$\alpha='+('2','0','-1')[j]+'$',transform=ax.transAxes,ha='right',fontsize=9)
        ax.set_title(R.NULL_NAMES[kind].replace('\n',' '),fontsize=8,pad=4)
    axes[0].set_ylabel(r'Heading $C(\ell)$');csvwrite(DATA/'B_null_heading_acf.csv',acf)
    label(fig,grid[2],'C');kymographs(fig,grid[2])
    label(fig,grid[3],'D');axes=panel_axes(fig,grid[3]);alpha=primary(readold('local_alpha'))
    for j,(ax,m) in enumerate(zip(axes,MODELS)):R.alpha_panel(ax,alpha,m,ylabel=j==0);ax.text(.02,.95,'C: light → dark',transform=ax.transAxes,fontsize=6,va='top')
    csvwrite(DATA/'D_local_alpha_points.csv',alpha)
    pair=grid[4].subgridspec(1,2,width_ratios=[1.05,3.4],wspace=.5)
    label(fig,pair[0],'E');ax=fig.add_subplot(pair[0]);fold=primary(readold('noise_tolerance_folds'));R.fold_panel(ax,fold);ax.set_ylim(.95,5.4);ax.set_xlabel('C (Spring ×20)');csvwrite(DATA/'E_tolerance_folds.csv',fold)
    label(fig,pair[1],'F');axes=panel_axes(fig,pair[1]);drift=csvread(DATA/'F_drift_ratio.csv')
    for ax,m in zip(axes,MODELS):
        lines(ax,drift,m,'Ds_over_smallest',True);ax.set_yscale('log');ax.axhline(1,color='.65',lw=.5)
        if m!='goryachev':ax.set_ylim(.2,8)
        else:ax.set_ylim(1e-7,2)
        ax.set_yticks([1e-6,1e-3,1] if m=='goryachev' else [.25,1,4],['$10^{-6}$','$10^{-3}$','1'] if m=='goryachev' else ['0.25','1','4'])
        ax.yaxis.set_minor_formatter(NullFormatter())
    axes[0].set_ylabel(r'$D_s/D_s(C_{\min})$')
    label(fig,grid[5],'G');closure_plot(panel_axes(fig,grid[5]),closure)
    label(fig,grid[6],'H');axes=panel_axes(fig,grid[6]);heading=readold('heading_amplitude_condition_metrics')
    for ax,m in zip(axes,MODELS):
        lines(ax,heading,m,'moment_residual_fraction',True);ax.axhline(0,color='.7',lw=.5)
    axes[0].set_ylabel('Remainder share R');csvwrite(DATA/'H_remainder_points.csv',heading)
    label(fig,grid[7],'I');takeover_plot(panel_axes(fig,grid[7]),csvread(DATA/'I_takeover_points.csv'))
    fig.savefig(OUT/'figure4_main.tmp.pdf');fig.savefig(OUT/'figure4_main.tmp.png',dpi=180);plt.close(fig)
    (OUT/'figure4_main.tmp.pdf').replace(OUT/'figure4_main.pdf');(OUT/'figure4_main.tmp.png').replace(OUT/'figure4_main.png')

def collect_modes():
    exemplars.collect()
    rates=[];curves=[]
    for p in sorted((RUN/'precision').glob('c*_rates.csv')):rates.extend(csvread(p))
    for p in sorted((RUN/'precision').glob('c*_curves.csv')):curves.extend(csvread(p))
    csvwrite(DATA/'targeted_mode_rates.csv',rates);csvwrite(DATA/'targeted_mode_curves.csv',curves)
    exemplar_rates=csvread(DATA/'exemplar_rates.csv')
    combined=[dict(r,mode_source='targeted independent histories') for r in rates]
    targeted_ids={r['condition_id'] for r in rates}
    combined.extend(dict(r,mode_source='existing exemplars, two seeds') for r in exemplar_rates if r['condition_id'] not in targeted_ids)
    forcing={(r['model'],num(r,'C_label'),num(r,'f')):r.get('forced_combined','unknown') for r in readold('forcing_conditions')}
    for r in combined:r['forced_combined']=forcing.get((r['model'],num(r,'C_label'),num(r,'f')),'unknown')
    csvwrite(DATA/'all_input_mode_rates.csv',combined)
    # Directly test the extra assumption linking complex mode turnover to
    # heading turnover. Their weighting and |P| eligibility are different.
    phase_curves=curves+[r for r in csvread(DATA/'exemplar_curves.csv') if r['condition_id'] not in targeted_ids]
    phase_by={}
    for r in phase_curves:
        if num(r,'n')==1:phase_by.setdefault(r['condition_id'],[]).append(r)
    inherited={(r['model'],num(r,'C_label'),num(r,'f')):num(r,'k') for r in readold('analysis_lag_metrics')}
    phase_rates=[]
    for r in combined:
        if num(r,'n')!=1:continue
        rr=phase_by.get(r['condition_id'],[])
        admissible=[s for s in rr if .03<num(s,'phase_only_a1')<.995 and num(s,'phase_eligible_pairs')>0]
        gamma_phase=phase_lag=np.nan;crossed=False
        if admissible:
            nearest=min(admissible,key=lambda s:abs(num(s,'phase_only_a1')-.5));phase_lag=num(nearest,'lag')
            gamma_phase=-np.log(num(nearest,'phase_only_a1'))/phase_lag
            crossed=any(num(s,'phase_only_a1')>=.5 for s in rr) and any(num(s,'phase_only_a1')<=.5 for s in rr)
        k=inherited.get((r['model'],num(r,'C_label'),num(r,'f')),np.nan)
        ratio=num(r,'gamma')/gamma_phase if gamma_phase>0 else np.nan
        phase_rates.append(dict(model=r['model'],C_label=r['C_label'],f=r['f'],condition_id=r['condition_id'],mode_source=r['mode_source'],
            gamma_complex_a1=num(r,'gamma'),complex_status=r.get('status'),gamma_phase_a1=gamma_phase,phase_analysis_lag=phase_lag,
            complex_half_crossing=r.get('half_crossing_observed'),phase_half_crossing=crossed,
            complex_over_phase=ratio if crossed and str(r.get('half_crossing_observed'))=='True' and r.get('status')=='resolved' else np.nan,
            inherited_full_campaign_heading_k=k,complex_over_inherited_k=num(r,'gamma')/k if k>0 else np.nan,
            comparison='complex: all endpoints, amplitude weighted; phase: both |P|>=0.1, unit weight; inherited k also has a different independent-history cohort and original lag grid'))
    csvwrite(DATA/'mode1_complex_phase_rate_comparison.csv',phase_rates)
    comparisons,collapse=model_test(combined,primary(readold('local_alpha')))
    for r in collapse:
        r['passes_rate_screen']=num(r,'n_comparisons')>=10 and num(r,'log_gamma_RMSE')<=np.log(1.5)
        r['passes_slope_screen']=num(r,'n_alpha_comparisons')>=10 and num(r,'alpha_slope_RMSE')<=.5
    sources={(r['model'],num(r,'C_label'),num(r,'f')):r for r in combined}
    for r in comparisons:
        source=sources[r['model'],num(r,'C_label'),num(r,'f')]
        r['mode_source']=source['mode_source'];r['forced_combined']=source['forced_combined']
    by_input=[]
    for model in MODELS:
        for f in sorted({num(r,'f') for r in comparisons if r['model']==model}):
            rr=[r for r in comparisons if r['model']==model and num(r,'f')==f];aa=[r for r in rr if np.isfinite(num(r,'alpha'))]
            by_input.append(dict(model=model,f=f,n_size_predictions=len(rr),n_alpha_comparisons=len(aa),
                log_gamma_RMSE=np.sqrt(np.mean([num(r,'log_prediction_error')**2 for r in rr])) if rr else np.nan,
                alpha_slope_RMSE=np.sqrt(np.mean([(num(r,'alpha')-num(r,'local_mode_slope'))**2 for r in aa])) if aa else np.nan,
                stratum=rr[0]['input_stratum']))
    csvwrite(DATA/'K_alpha_mode_slope_points.csv',comparisons);csvwrite(DATA/'mode_collapse_summary.csv',collapse)
    csvwrite(DATA/'mode_collapse_by_input.csv',by_input)
    # Only a completed independent-seed analysis could support a main-figure
    # explanation; no promotion while native jobs or precision topups are pending.
    targets=read(RUN/'manifest.json')['conditions'] if (RUN/'manifest.json').exists() else []
    complete=(len(rates)==len(targets)*33 and bool(targets))
    precision=csvread(META/'modes_precision.csv')
    precision_ok=complete and len(precision)==len(targets) and all(num(r,'n_seeds')>=24 and (num(r,'SE_log_gamma')<=.1 or num(r,'n_seeds')==96) for r in precision)
    supported=precision_ok and len(collapse)==4 and all(num(r,'n_comparisons')>=10 and num(r,'n_alpha_comparisons')>=10 and num(r,'log_gamma_RMSE')<=np.log(1.5) and num(r,'alpha_slope_RMSE')<=.5 for r in collapse)
    sparse=any(num(r,'n_comparisons')<10 or num(r,'n_alpha_comparisons')<10 for r in collapse)
    final_statement=('Independent mode data remain inconclusive because too few supported cross-size rate/slope comparisons are available; the strong regime remains unexplained.' if sparse else
       'Independent mode data fail the declared transition-and-stronger collapse/slope checks; the strong-input regime remains unexplained.')
    decision=dict(generated_utc=GENERATED,targeted_rate_conditions=len(rates)//33,designed_targeted_conditions=len(targets),precision_complete=precision_ok,
       supports_mode_explanation=bool(supported),main_panels='A–K' if supported else 'A–I',J_K_location='main' if supported else 'supplement',
       decision=('Independent mode-collapse test not yet assessed: targeted precision evidence is incomplete; the strong regime is currently unexplained.' if not precision_ok else
          final_statement if not supported else 'Independent mode data support the declared transition-and-stronger collapse and slope checks; J–K are included in the main figure.'),
       transparent_checks='all targeted conditions evaluated; SE ln gamma1 <=0.1 or96 cap; fixed round 5 transition and stronger input stratum; per-model leave-one-size-out ln-gamma RMSE<=ln1.5; alpha-slope RMSE<=0.5; >=10 supported rate and alpha comparisons/model',
       caveat='Numerical tolerances are a declared descriptive screen, not a significance test; main J/K require scientific support, not a fitted curve.')
    write(META/'mode_explanation_decision.json',decision)
    return rates,exemplar_rates,decision

def supplement(rates,exemplar_rates):
    promoted=read(META/'mode_explanation_decision.json')['supports_mode_explanation']
    temporary=OUT/'figure4_supplement.tmp.pdf'
    with PdfPages(temporary) as pdf:
        def page(title,caption='',rows=1,height=4.0):
            fig,axes=plt.subplots(rows,4,figsize=(11.7,height),squeeze=False)
            fig.subplots_adjust(left=.07,right=.975,bottom=.18,top=.83,wspace=.38,hspace=.65)
            fig.text(.055,.965,title,ha='left',va='top',fontsize=12,color='black')
            if caption:fig.text(.06,.035,caption,ha='left',va='bottom',fontsize=7,wrap=True)
            return fig,axes
        def save(fig,letter=None):
            pdf.savefig(fig)
            if promoted and letter:
                for item in list(fig.texts):item.remove()
                fig.text(.015,.955,letter,color='black',fontweight='bold',fontsize=12)
                vector=OUT/'supporting/vectors'/f'main_{letter}.pdf';vector.parent.mkdir(parents=True,exist_ok=True)
                fig.savefig(vector);fig.savefig(OUT/f'figure4_main_{letter}.png',dpi=180)
            plt.close(fig)
        birth=readold('birth_heading_condition_metrics');fig,axes=page('S4.1  Certified persistent-secondary birth rates',
           'Rates are certified lower bounds per actual perimeter–time. Uncertified onsets can contain split descendants; native local parameters are unchanged.')
        for ax,m in zip(axes[0],MODELS):lines(ax,birth,m,'certified_birth_rate_lower_bound',True);ax.set_yscale('symlog',linthresh=1e-8)
        axes[0,0].set_ylabel('Certified births / perimeter / time');save(fig)
        fig,axes=page('S4.1b  Native heading correlations at the transition input',
            'Superseded round 5 model heading-correlation row, at the unchanged selected input. These are the complete independent-history estimates; source CSV condition_lag_metrics.csv.')
        for j,(ax,m) in enumerate(zip(axes[0],MODELS)):
            chosen=num(next(r for r in readold('model_D_selection') if r['model']==m),'f')
            rr=[r for r in primary(readold('condition_lag_metrics')) if r['model']==m and num(r,'f')==chosen]
            R.correlation(ax,rr,m,legend=j==0);ax.set_title(NAMES[m]+f', f={chosen:g}',fontsize=8)
        axes[0,0].set_ylabel('Heading correlation');save(fig)
        heading=readold('heading_amplitude_condition_metrics');fig,axes=page('S4.2  Motion, amount and membership shares',
            'Round5 native-amplitude symmetric Shapley convention. Signed fractions can be negative or exceed one; no change to source profiles, eligibility or events.',rows=3,height=8.5)
        for i,key in enumerate(['motion_fraction','coexisting_amount_fraction','birth_loss_fraction']):
            for ax,m in zip(axes[i],MODELS):lines(ax,heading,m,key,True);ax.axhline(0,color='.7',lw=.5)
            axes[i,0].set_ylabel(['Motion share','Amount share','Membership share'][i])
        save(fig)
        fig,axes=page('S4.3  Pooled finite-horizon takeovers',
            'Historical pooled f >0.01 diagnostic, shown only for positive distance. Main I uses a separately selected input per model and >=200 births at every size.')
        takeover_plot(axes[0],readold('takeover_by_size_distance'))
        for ax,m in zip(axes[0],MODELS):ax.set_title(NAMES[m]+', pooled f >0.01',fontsize=8)
        save(fig)
        # All inputs: same plot, four facets; zero/short follow-up exclusions are
        # retained in source CSV and are never treated as unsuccessful births.
        take=readold('takeover_by_condition_distance')
        for binid in (1,2):
            fig,axes=page('S4.4  Takeover probability at every input: '+('0 < d ≤ C/6' if binid==1 else 'd > C/6'),
                'Common follow-up20τ; each curve is one size. Certified no-overlap subset only. Inputs with no certified births remain missing.')
            for ax,m in zip(axes[0],MODELS):lines(ax,[r for r in take if num(r,'distance_bin')==binid],m,'probability');ax.set_ylim(-.02,1.02)
            axes[0,0].set_ylabel('Takeover probability');save(fig)
        exact=csvread(DATA/'exemplar_exact_masks.csv')
        for conv in ['point_matched_mass','exact_full_mask','exact_contour_mask','round5_point_native_amplitude','round5_exact_native_amplitude']:
            rr=[]
            for model in MODELS:
                for size in sizes(model):
                    for f in sorted({num(r,'f') for r in exact if r['model']==model}):
                        ss=[r for r in exact if r['model']==model and num(r,'C_label')==size and num(r,'f')==f and r['convention']==conv]
                        den=sum(num(s,'total_heading_loss') for s in ss)
                        if ss:rr.append(dict(model=model,C_label=size,f=f,remainder_fraction=sum(num(s,'remainder_fraction')*num(s,'total_heading_loss') for s in ss if num(s,'total_heading_loss')>0)/den if den else np.nan))
            csvwrite(DATA/('exact_mask_curve_'+conv+'.csv'),rr)
            fig,axes=page('S4.5  Own-mask moment check: '+conv.replace('_',' '),
                'Existing seeds 0/1; saved binned float16 profiles. Native-mass denominator in first 3 variants; last 2 retain round 5 region normalization and native amplitude scaling. Shape enters amount/shape in first 3 and motion/shape in last 2. Algebraic reconstruction is not causal validation.')
            for ax,m in zip(axes[0],MODELS):lines(ax,rr,m,'remainder_fraction');ax.axhline(0,color='.7',lw=.5)
            axes[0,0].set_ylabel('Remainder share');save(fig)
        # Main-candidate J: both weak and round 5 transition input, one rate plot
        # with line style distinguishing the inputs; faceting only by model.
        jr=rates or exemplar_rates;fig,axes=page('J  Fourier-mode decorrelation versus physical wavenumber',
             'Solid weak f=0.001 (existing exemplars); dashed round 5 transition input. Open triangles with arrows bound rates when half crossing is outside the saved lag range. Targeted records replace strong curves after precision analysis.',height=6)
        selection=readold('model_D_selection');jp=[]
        for ax,m in zip(axes[0],MODELS):
            transition=num(next(r for r in selection if r['model']==m),'f')
            for f,ls in ((.001,'-'),(transition,'--')):
                available=exemplar_rates if f==.001 else jr
                for size in sizes(m):
                    rr=sorted([r for r in available if r['model']==m and num(r,'f')==f and num(r,'C_label')==size and num(r,'n')>0],key=lambda r:num(r,'q'))
                    resolved=[num(r,'gamma') if r.get('status')=='resolved' and str(r.get('half_crossing_observed'))=='True' else np.nan for r in rr]
                    forced=forcing_flag(m,size,f)
                    ax.plot([num(r,'q') for r in rr],resolved,ls=ls,color=color(m,size),lw=.8,marker='o',ms=1.5,mfc='white' if forced else color(m,size),mew=.5)
                    for r in rr:
                        bound=num(r,'crossing_rate_bound');direction=r.get('crossing_bound_direction')
                        if np.isfinite(bound) and bound>0 and direction in ('lower','upper'):
                            ax.plot(num(r,'q'),bound,marker='^' if direction=='lower' else 'v',mfc='white',mec=color(m,size),ms=2.2,mew=.5)
                            ax.annotate('',xy=(num(r,'q'),bound*(1.35 if direction=='lower' else .74)),xytext=(num(r,'q'),bound),arrowprops=dict(arrowstyle='->',lw=.35,color=color(m,size)))
                    jp.extend(dict(r,input_role='weak' if f==.001 else 'transition',forced_combined=forced) for r in rr)
            ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Physical wavenumber q');style(ax,m)
            rr=[r for r in jp if r['model']==m and num(r,'gamma')>0]
            if rr:
                q0=np.median([num(r,'q') for r in rr]);g0=np.median([num(r,'gamma') for r in rr]);xx=np.array([q0/2,q0*2])
                ax.plot(xx,g0*(xx/q0)**2,color='.5',lw=.7,ls=':');ax.plot(xx,[g0]*2,color='.5',lw=.7,ls='-.')
                ax.text(xx[-1],g0*4,'2',fontsize=7,color='.4');ax.text(xx[-1],g0,'0',fontsize=7,color='.4')
            else:ax.set_xlim(.05,10);ax.set_ylim(1e-6,.1)
        axes[0,0].set_ylabel('Decorrelation rate γ');csvwrite(DATA/'J_mode_rate_points.csv',jp);save(fig,'J')
        comp=csvread(DATA/'K_alpha_mode_slope_points.csv');fig,axes=page('K  Size exponent versus independently predicted local mode slope',
             'Mode slope uses the nearest 8 mode points from other sizes at the same model and input; the target size is excluded. One-sided predictions are open squares. Forced conditions have open circles. This tests collapse rather than replotting the same size slope.',height=6)
        for ax,m in zip(axes[0],MODELS):
            rr=[r for r in comp if r['model']==m];pts=[]
            for r in rr:
                if not np.isfinite(num(r,'alpha')+num(r,'local_mode_slope')):continue
                ax.plot(num(r,'local_mode_slope'),num(r,'alpha'),marker='s' if r.get('extrapolated')=='True' else 'o',mfc='white' if r.get('extrapolated')=='True' or r.get('forced_combined')=='True' else color(m,num(r,'C_label')),mec=color(m,num(r,'C_label')),ms=2.8,mew=.6)
                pts.extend([num(r,'alpha'),num(r,'local_mode_slope')])
            limit=max(3.,max(abs(np.array(pts))) if pts else 3.);limit=min(limit,20)
            ax.plot([-limit,limit],[-limit,limit],color='.6',lw=.7);ax.set_xlim(-limit,limit);ax.set_ylim(-limit,limit);ax.set_xlabel('Local slope of γ(q)');style(ax,m)
        axes[0,0].set_ylabel('Figure4 local α');save(fig,'K')
        # Equal-time spectra every designed input; missing conditions appear as
        # empty facets instead of a surrogate spectrum.
        fs=sorted({num(r,'f') for r in exemplar_rates})
        for f in fs:
            fig,axes=page(f'S4.6  Equal-time spectrum, f = {f:g}',
                r'$S(q)=\langle|a_n|^2\rangle/\langle a_0^2\rangle$; n=1–32, two existing exemplars/condition. Spring q uses mean C. n* and C/n* are in exemplar_rates.csv.')
            for ax,m in zip(axes[0],MODELS):
                for size in sizes(m):
                    rr=sorted([r for r in exemplar_rates if r['model']==m and num(r,'C_label')==size and num(r,'f')==f and num(r,'n')>0],key=lambda r:num(r,'q'))
                    ax.plot([num(r,'q') for r in rr],[num(r,'S') for r in rr],color=color(m,size),lw=.8)
                ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Physical wavenumber q');style(ax,m)
                available=[r for r in exemplar_rates if r['model']==m and num(r,'f')==f and num(r,'S')>0 and num(r,'n')>0]
                if not available:ax.set_xlim(.05,10);ax.set_ylim(1e-12,1)
            axes[0,0].set_ylabel('Equal-time spectrum S(q)');save(fig)
        # Phase-only source equivalence is shown as a difference curve, directly
        # identifying finite binning/precision effects.
        phase=csvread(DATA/'exemplar_phase_check.csv');fig,axes=page('S4.7  First-mode phase reproduces native heading correlation',
            'Difference between Fourier a1 phase correlation and native compact P heading, at identical exact saved endpoints. Eligibility uses each stored |P|≥0.1; nonzero differences include threshold crossing and bin/float16 precision.')
        for ax,m in zip(axes[0],MODELS):
            for size in sizes(m):
                rr=[r for r in phase if r['model']==m and num(r,'C_label')==size and np.isfinite(num(r,'difference'))]
                ax.scatter([num(r,'lag') for r in rr],[num(r,'difference') for r in rr],s=2,color=color(m,size),alpha=.5)
            ax.axhline(0,color='.6',lw=.6);ax.set_xscale('log');ax.set_xlim(20 if m=='spring_ring' else 200,2000 if m=='spring_ring' else 20000);ax.set_xlabel('Lag');style(ax,m)
        axes[0,0].set_ylabel('Phase C − native heading C');save(fig)
        fig,axes=page('S4.7b  Complex-mode and phase-only first-mode rates',
             'Same mode-history cohort and nearest-half lag rule. Complex correlation is amplitude weighted over all endpoints; phase-only correlation uses |P|≥0.1. Only observed half crossings enter the ratio; all statuses and comparisons with inherited full-campaign k are in CSV.')
        phase_rates=csvread(DATA/'mode1_complex_phase_rate_comparison.csv')
        for ax,m in zip(axes[0],MODELS):
            lines(ax,phase_rates,m,'complex_over_phase');ax.set_yscale('log');ax.axhline(1,color='.55',lw=.6)
        axes[0,0].set_ylabel('γ₁ complex / γ₁ phase');save(fig)
        # Normalized complex correlations themselves, all modes and conditions.
        curves=csvread(DATA/'exemplar_curves.csv')
        for m in MODELS:
            inputs=sorted({num(r,'f') for r in curves if r['model']==m});nrows=max(3,int(np.ceil(len(inputs)/5)))
            fig,axes=plt.subplots(nrows,5,figsize=(11.7,8.3 if nrows==3 else 10.8));fig.subplots_adjust(left=.065,right=.98,bottom=.09,top=.9,hspace=.65,wspace=.42)
            fig.text(.055,.97,'S4.8  Normalized complex autocorrelations: '+NAMES[m],fontsize=12,va='top')
            for ax,f in zip(axes.flat,inputs):
                for size in sizes(m):
                    # First, second, fourth and eighth modes display the same
                    # measure across scales; complete0–32 curves are in CSV.
                    for mode,ls in ((1,'-'),(2,'--'),(4,':'),(8,'-.')):
                        rr=sorted([r for r in curves if r['model']==m and num(r,'f')==f and num(r,'C_label')==size and num(r,'n')==mode],key=lambda r:num(r,'lag'))
                        ax.plot([num(r,'lag') for r in rr],[num(r,'correlation') for r in rr],ls=ls,color=color(m,size),lw=.55)
                ax.axhline(.5,color='.7',lw=.4);ax.set_xscale('log');ax.set_ylim(-.25,1.1);ax.set_title(f'f = {f:g}',fontsize=7);style(ax)
            for ax in list(axes.flat)[len(inputs):]:ax.axis('off')
            fig.text(.06,.025,'n=1 solid,2 dashed,4 dotted,8 dash-dot; size light→dark. All33 modes including n=0 are retained in compact CSV.',fontsize=7)
            save(fig)
    temporary.replace(OUT/'figure4_supplement.pdf')

def documents(decision):
    selection=csvread(DATA/'I_takeover_selection.csv');phase=csvread(DATA/'exemplar_phase_check.csv');ex=csvread(DATA/'exemplar_exact_masks.csv')
    diffs=[abs(num(r,'difference')) for r in phase if np.isfinite(num(r,'difference'))]
    phase_max=max(diffs) if diffs else np.nan
    comparisons=[]
    for r in ex:
        if r['convention']!='exact_full_mask':continue
        p=next((s for s in ex if s['condition_id']==r['condition_id'] and s['replicate']==r['replicate'] and s['convention']=='point_matched_mass'),{})
        if np.isfinite(num(r,'remainder_fraction')+num(p,'remainder_fraction')):comparisons.append(abs(num(r,'remainder_fraction'))-abs(num(p,'remainder_fraction')))
    negative=sum(v<0 for v in comparisons);original_comparisons=[];paired=[]
    for r in ex:
        if r['convention'] not in ('exact_full_mask','round5_exact_native_amplitude'):continue
        baseline='point_matched_mass' if r['convention']=='exact_full_mask' else 'round5_point_native_amplitude'
        p=next((s for s in ex if s['condition_id']==r['condition_id'] and s['replicate']==r['replicate'] and s['convention']==baseline),{})
        a,b=num(r,'remainder_fraction'),num(p,'remainder_fraction')
        paired.append(dict(condition_id=r['condition_id'],model=r['model'],C_label=r['C_label'],f=r['f'],replicate=r['replicate'],
            exact_convention=r['convention'],point_convention=baseline,point_R=b,exact_R=a,absolute_R_change=abs(a)-abs(b)))
        if r['convention']=='round5_exact_native_amplitude' and np.isfinite(a+b):original_comparisons.append(abs(a)-abs(b))
    csvwrite(DATA/'exact_mask_paired_comparison.csv',paired)
    original_shrink=sum(v<0 for v in original_comparisons)
    spectral_summary=[]
    er=csvread(DATA/'exemplar_rates.csv')
    for model in MODELS:
        rr=[r for r in er if r['model']==model and num(r,'n')==1]
        spectral_summary.append(dict(model=model,conditions=len(rr),nonzero_nstar_above1=sum(num(r,'n_star')>1 for r in rr),
            nstar_min=min([num(r,'n_star') for r in rr],default=np.nan),nstar_max=max([num(r,'n_star') for r in rr],default=np.nan),
            lambda_star_min=min([num(r,'lambda_star') for r in rr],default=np.nan),lambda_star_max=max([num(r,'lambda_star') for r in rr],default=np.nan)))
    csvwrite(DATA/'spectrum_dominant_mode_summary.csv',spectral_summary)
    collapse_rows=csvread(DATA/'mode_collapse_summary.csv')
    screen_table='| Model | Supported size predictions | ln-rate RMSE | α–slope RMSE | Rate screen | Slope screen |\n|---|---:|---:|---:|---|---|\n'+'\n'.join(
        '| '+NAMES[r['model']]+f' | {num(r,"n_comparisons"):g} | {num(r,"log_gamma_RMSE"):.3g} | {num(r,"alpha_slope_RMSE"):.3g} | '+str(r.get('passes_rate_screen','pending'))+' | '+str(r.get('passes_slope_screen','pending'))+' |' for r in collapse_rows)
    text=[
      'Figure 4. Enlargement preserves heading as noisy polarity changes regime.',
      'A, Detector-free null kymographs at C20 and C60: identical physical diffusion, size-independent jumping, or perimeter-proportional jumping. Each facet shows the same activity measure.',
      'B, Null heading autocorrelation; the corresponding persistence size exponents are 2, 0 and −1. C20 is light and C60 dark.',
      'C, Native-model kymographs at the unchanged round 5 transition inputs: Wave-pinning 0.03, Holmes 0.02, Spring 0.005 and Goryachev 0.01. Positions are centred at the initial heading; time is elapsed after burn-in. Goryachev uses the explicitly labelled per-size activity ranges.',
      'D, Local persistence exponent α versus input fraction f, with the inherited fixed three-size neighbourhood (two size values at endpoints), whole-history 95% intervals and light-to-dark size ramps. In D and F–H, sizes are C = 20, 30, 40, 50 and 60 (Spring 1, 1.5, 2, 2.5 and 3), from light to dark. Open circles retain forced-condition flags; missing rates are unavailable, not zero.',
      'E, Noise-tolerance fold at fixed heading-correlation lag 2000 (Spring 200), relative to the smallest ring. Spring circumference is multiplied by 20 for plotting only.',
      'F, Physical within-tenure drift coefficient normalized by the smallest size at the same input. Native equations and the five-frame dominant-tenure rule are unchanged.',
      'G, Drift plus counted five-frame-dominant takeover closure over every designed input. Solid lower and dashed upper limits include every unscored takeover angle with decorrelation loss between 0 and 2; shaded intervals are partial-identification bounds, not confidence intervals. Grey line: unity. Open upward triangles with arrows mark an unbounded upper limit where drift is unavailable; drift is bounded below by zero there. Crosses at the top mark unavailable observed heading rates. No mean angle is imputed to unscored events.',
      'H, Signed remainder share versus input at each size, using the round 5 native-amplitude point-region convention. It contains amplitude, within-region shape, unresolved activity and background; it is descriptive. Exact finite-width moments can be attenuated again by the original native-amplitude scaling, so their nonshrinking remainder does not establish failure of a domain mechanism. Exact own-mask checks appear in the supplement.',
      'I, Certified takeover probability within 20 input correlation times, separately selected per model. Positive birth-to-incumbent distances only; filled circles 0 < d ≤ C/6, open squares d > C/6. '+ '; '.join(NAMES[r['model']]+f' f={num(r,"f"):g}, minimum total births={num(r,"minimum_births"):g}' for r in selection)+'. The selected designed input is nearest in log distance to the smallest ring’s f½ among those with at least 200 total certified births at every size. Intervals resample independent histories 2000 times; [0,1] bounds replace intervals when fewer than 97.5% of resamples have a defined probability. The two distance bins can have very different denominators. Certification and finite follow-up limit inference about all births.',
      ('J, Fourier-mode rates against physical wavenumber at weak and transition inputs. K, Figure4 size exponent against local mode slope estimated from other sizes; the target size is excluded. '+decision['decision']) if decision['supports_mode_explanation'] else
      'The strong-input regime is currently unexplained. Candidate J–K, spectra and exact-mask diagnostics remain supplementary; two exemplar seeds alone cannot validate collapse. '+decision['decision']]
    if decision['supports_mode_explanation']:
        text[-1]='J, Fourier-mode decorrelation rates against physical wavenumber, with all five sizes overlaid at weak and transition inputs. Solid weak curves use existing exemplars; dashed transition curves use independent targeted histories. Grey guides have slopes 2 and 0. Half-crossing bounds are open triangles with arrows; forced conditions have open circles.'
        text.extend(['K, Figure4 local persistence exponent against the local logarithmic mode-rate slope predicted from other sizes at the same input; the target size is excluded. Grey line: equality. One-sided predictions are open squares; forced conditions have open circles. Weak-input points use existing exemplars, targeted inputs use independent mode-only histories.',decision['decision']])
    text.insert(1,'Model panels use unchanged positive lognormal input with mean f times the native reference: Wave-pinning 0.1, Holmes 1, Spring 4 and Goryachev 6.4. The Gaussian driving field has unit variance, lognormal noise amplitude σ = 1, spatial parameter 1.5 and OU correlation time 50 (Spring basal-physical 0.15 and 5). Spring input is stationary in material coordinate, so its physical correlation length scales with circumference. Native clocks are not calibrated between models. Existing statistics comprise 260 conditions and 45,696 histories; targeted independent mode-only histories are additional and do not alter A–I.')
    (OUT/'legend.md').write_text('\n\n'.join(text)+'\n')
    readme=f'''# Unexpected

The main figure now exposes the strong-input closure intervals and signed remainder at every input. {decision['decision']}

Preserving the original round 5 normalization and native-amplitude convention, exact own-mask moments reduce the magnitude of R in only {original_shrink}/{len(original_comparisons)} eligible seed-condition comparisons. They do not generally remove the main-panel remainder. Inserting finite-width moments into a normalization built for unit-modulus point regions can attenuate them again through native-amplitude scaling; this nonshrinking remainder does not show that a domain mechanism fails. With a native-total denominator, the exact-mask moment reconstructs all represented excess by algebra; its smaller remainder cannot validate a causal domain explanation.

The birth-count requirement shifts the takeover comparison well above the smallest ring’s half-memory input: {', '.join(NAMES[r['model']]+' f='+str(num(r,'f')) for r in selection)}. These are per-model common inputs; pooled strong-input probabilities no longer occupy the main figure.

## Existing exemplar evidence

Completed conditions: {read(META/'exemplar_status.json')['complete_conditions']}/260, two saved profiles each. Phase-only $a_1$ and same-endpoint native heading correlation differ by at most {phase_max:.6g}; finite float16/binning and $|P|$ eligibility differences remain recorded. Full 0–32 complex correlations and equal-time spectra are exported. Dominant nonzero mode $n^*$ and wavelength $C/n^*$ are reported in spectrum_dominant_mode_summary.csv and every condition's mode-rate rows.

Across all designed inputs, { '; '.join(NAMES[r['model']]+f': n*={r["nstar_min"]:g}–{r["nstar_max"]:g}, C/n*={r["lambda_star_min"]:.3g}–{r["lambda_star_max"]:.3g} native length units; n*>1 in {r["nonzero_nstar_above1"]}/{r["conditions"]} conditions' for r in spectral_summary)}. A dominant first harmonic does not establish a single detected domain. The separate complex/phase rate check measures whether the amplitude-weighted complex first mode has the same rate as unit-weight heading correlation; their equality is tested, not assumed (mode1_complex_phase_rate_comparison.csv and supplement S4.7b).

Exact full-basin own-mask moments reduced the magnitude of the matched-denominator signed remainder in {negative}/{len(comparisons)} seed-condition comparisons. This is a bounded-exemplar sensitivity, not a recovery of discarded full-campaign masks. If the full masks cover all activity excess, their exact moments reconstruct P by algebra and R is zero by construction; this does not validate a domain mechanism. A change of moment denominator and movement of within-region shape into the amount/shape factor are declared explicitly; it must not be confused with the round 5 amplitude-scaled point convention in main H.

## Current targeted campaign

Generated {GENERATED}. Initial design 200 conditions × 24 seeds; up to 96 per condition where SE ln γ1 >0.1 (unresolved estimates also top up). Same generator, qualified basal, BDF numerics, half-horizon burn-in and full horizons as round 3. Only modes 0–32, |P|, time and circumference are durable. No new native profile history is written. Native integration runs through SLURM, using existing scaling_sandbox Conda.

Targeted precision-analyzed conditions: {decision['targeted_rate_conditions']}/200. See supporting/metadata/modes_submission.json and modes_precision.csv. The independent precision controller updates mode evidence and rebuilds artifacts after the authorized additions. The source round 5 deliverables are unchanged.

The following point-estimate diagnostics use each model's fixed round 5 transition input and all stronger inputs; all-input results are retained separately. These descriptive screens are not significance tests, and the table remains provisional while allocated cohorts are incomplete.

{screen_table}

## Rebuild

`{PYTHON} -B {OUT.relative_to(ROOT)}/rebuild_figure4.py`

This one entry point reads source data, collects completed analyses, exports every plotted point and rebuilds main PDF/PNG, supplement PDF and prose. Scheduler submission is separate and explicit: `{PYTHON} -B {OUT.relative_to(ROOT)}/scripts/modes_campaign.py submit` (do not repeat while arrays are active).

Biological scope: length represents WGD-associated enlargement. Simulated one-dimensional activity-heading persistence is distinct from observed migration summaries and conditional OU velocity-persistence estimates.
'''
    (OUT/'scientific_readme.md').write_text(readme)
    methods=r'''# Methods and provenance

The scientific baseline is manuscript_outline/outline_v11.pdf and JSON, with the completed round 5 Figure4 and unchanged compact round 3 native campaign. Rebuild uses the existing scaling_sandbox Conda environment. Source paths and hashes are in supporting/metadata/provenance.json. Existing evidence comprises 260 conditions and 45,696 native histories; no source data are overwritten.

## Fourier representation

For saved cell-centre phases θj, use the positive complex convention that matches the pre-existing first activity moment:

$$a_n(t)=N^{-1}\sum_j u_j(t)e^{\mathrm{i}n\theta_j},\quad P(t)=a_1(t)/a_0(t),\quad q_n=2\pi n/\overline{C}. \tag{S4.1}$$

Spring uses material phases and mean physical circumference. Modes0–32 are stored, bounded by Nyquist. Equal-time spectrum is the mean squared mode amplitude divided by mean squared zero-mode amplitude. n* maximizes the nonzero spectrum and λ*=C/n*. There is no spatial or temporal mean subtraction.

$$C_n(\ell)=\frac{\sum_{t,s}\operatorname{Re}[a_{n,s}(t+\ell)\overline{a_{n,s}(t)}]}{\frac12\sum_{t,s}(|a_{n,s}(t+\ell)|^2+|a_{n,s}(t)|^2)},\qquad \gamma_n=-\frac{\ln C_n(\ell_*)}{\ell_*}. \tag{S4.2}$$

Use the inherited sqrt2-spaced lag rule, preserving all original heading lags 250–16000 (Spring 25–1600) and extending the same geometric grid down to the saved cadence for faster higher modes. All lags are snapped to exact saved endpoints; no interpolation or refitting of native trajectories occurs. The admissible correlation lies strictly between0.03 and0.995; choose the lag nearest0.5, with the original finite-lag estimator. `half_crossing_observed` flags whether both sides of0.5 are observed. A rate cannot be resolved if no lag is admissible. Exemplar frame spacing is nonuniform: exact matching of saved times replaces an incorrect median-cadence assumption. Phase-only a1 uses |P|≥0.1 at both endpoints and is compared with compact native heading at identical endpoints. Complex-first-mode and phase-only-first-mode rates are also compared at their respective nearest-half lags within the same history cohort; their amplitude weighting and endpoint eligibility differ. Equality is part of the tested hypothesis, not a guaranteed identity. Complex mode correlations do not require domain detection or a |P| threshold.

Whole histories are resampled 2000 times. At the selected lag, SE lnγ is the sample standard deviation of valid log-rate draws. Fewer than 97.5% valid draws gives an unresolved precision status. Initial 24 histories top up to 96 when SE >0.1 or precision is unresolved; no more than 96 are authorized. Time, circumference,33 complex64 modes and float32 |P| cost approximately 4.32 GB initial or 17.29 GB maximum before compression and metadata, with a 34.57 GB reserve. These are payload estimates, not strict allocated-storage upper bounds; measured allocated blocks are authoritative. Every seed is a fresh unchanged native realization. No full profiles are durable.

## Missing-angle closure

$$k_{\mathrm{pred}}^{\mathrm{lo}}=q_1^2D_s+\frac{\sum_{\mathrm{scored}}(1-\cos\Delta\theta)}{T},\quad k_{\mathrm{pred}}^{\mathrm{hi}}=k_{\mathrm{pred}}^{\mathrm{lo}}+2\frac{N_{\mathrm{unscored}}}{T}. \tag{S4.3}$$

T is inherited eligible dominant-tenure exposure; use five-frame persistent global dominant changes and pre/post whole-P angle windows. Bounds contain every unscored event, including low-coherence intervals. They are identification intervals; drift/rate sampling confidence intervals are separate. Only unavailable observed rates prevent division. No weak-input or rare-multiple-domain filter is applied in main G.

## Exact region-mask sensitivity

Existing float16 exemplars retain at most400 spatial bins and2000 frames. Re-run the frozen half-contrast detector and maximum-overlap tracker at that saved resolution, with unchanged basal contrast. For every region integrate its own complete-basin mask, and separately its contour mask; never replace the region by a point at its centroid:

$$z_i(t)=\frac{\Delta x_t\sum_{j\in\mathcal M_i(t)}[u_j(t)-\min u(t)]e^{\mathrm i\theta_j}}{M_0(t)}. \tag{S4.4}$$

Native cell quadrature is exact for the saved piecewise-constant sampled profile and own mask, not a claim to recover the discarded native grid. Region shared-identity phase changes form motion; complex magnitudes include amount and within-region shape; membership inserts/removes endpoint identities. Symmetric three-factor Shapley changes add to the reconstructed moment change. Remainder closes the native P change. Project each change onto the unit-heading chord using the same signed loss identity as round 5. For a fair point-reduction comparison, the point variant uses the same masks and native-total denominator. Main H keeps the original round 5 native-amplitude convention; these changes of convention are a sensitivity, not silent redefinition. The two extra point/exact comparisons preserve region-total normalization and native-amplitude scaling. Finite-width exact moments then attenuate a representation built for unit-modulus point regions; scaling again can reintroduce an amplitude artefact. Their nonshrinking remainder does not establish failure of a domain mechanism. In this additional comparison, within-region complex shape changes are assigned with motion, and amount retains the complete-mask activity weight. Every eligible exact endpoint pair at the inherited campaign lag is used, with native |P|≥0.1; unobserved endpoint pairs are unavailable.

## Takeover selection and scope

The inherited no-overlap enclosure certifies a subset of persistent secondary births; exact parents and split classification are absent from full-campaign compact records. Use common20τ remaining record and five-frame global dominance for success. Select separately for each model the designed input nearest the smallest size’s f½ in log distance, requiring ≥200 total certified births at every size, then removing d≤0 from plotted probabilities. Report the two positive distance bins separately; distance retains the inherited complete-amount/global-contrast plateau approximation. Intervals resample whole independent histories including zero-event histories. These probabilities concern the certified subset over finite follow-up.

## Mode-collapse test

The tested claim is a size-independent γ(q) at fixed model and input, and α matching its local log slope at q1. Predictions use the nearest 8 mode points from other sizes; the target size is excluded. One-sided extrapolation remains flagged. No rates or native parameters are fitted to promote agreement. Completeness and a declared descriptive screen (per-model leave-one-size-out log-rate RMSE≤ln1.5, α-slope RMSE≤0.5,≥10 comparisons/model) are recorded. This screen is not a formal significance test. Until independent histories support the claim, J–K remain in the supplement and the strong regime is explicitly unexplained.
'''
    (OUT/'methods_provenance.md').write_text(methods)
    equations=read(OLD/'legend.json')['equations']
    equations[-1]['description']='C(ℓ) is activity-heading autocorrelation; mean C is actual circumference. Endpoint eligibility uses |P|≥0.1. Main G evaluates every designed input; unavailable observed rates remain explicit.'
    equations.extend([
      dict(math=r'a_n(t)=\frac{1}{N}\sum_{j=1}^{N}u_j(t)e^{\mathrm{i}n\theta_j},\qquad P(t)=\frac{a_1(t)}{a_0(t)},\qquad q_n=\frac{2\pi n}{\overline{C}}',number='4.5',description='Fourier modes use cell-centre material phases and the positive convention matching the native first moment. Heading is the phase of P; Spring physical wavenumber uses mean circumference.'),
      dict(math=r'\gamma_n=-\frac{\ln C_n(\ell_*)}{\ell_*},\qquad k(C)=\gamma(2\pi/C),\qquad \alpha=\frac{\mathrm{d}\ln\gamma}{\mathrm{d}\ln q}',number='4.6',description='The last two equalities are the detector-free hypothesis under a common size-independent mode curve, not an identity assumed in analysis. Complex correlations use the nearest-half lag rule; omitted half crossings are explicit bounds. '+decision['decision']),
      dict(math=r'k_{\mathrm{pred}}^{\mathrm{lo}}=q_1^2D_s+\frac{\sum_{\mathrm{scored}}(1-\cos\Delta\theta)}{T},\qquad k_{\mathrm{pred}}^{\mathrm{hi}}=k_{\mathrm{pred}}^{\mathrm{lo}}+\frac{2N_{\mathrm{unscored}}}{T}',number='4.7',description='All counted takeovers enter the interval: each unscored angular loss lies between 0 and 2. T is inherited eligible exposure. Missing drift has a nonnegative lower bound and unbounded upper limit; an unavailable observed heading rate prevents a numerical ratio.')])
    rel=str(OUT.relative_to(ROOT));section=dict(id='figure4',title=text[0],role='Size-dependent heading persistence, region-based diagnostic limits and a detector-free mode test.',
        assets=[dict(path=rel+'/figure4_main.png',native_pdf_page=rel+'/figure4_main.pdf',caption='Round6 A–I; one plot per letter; all-input closure bounds and remainder curves.',page_size=[12,16],page_title='',minimum_printed_width_inches=12)],
        legend=text[1:],narrative=['Larger rings preserve heading longer, while the size exponent weakens as input strengthens. Region-based closure and remainder now expose their limits across the full input range. The strong-input regime remains unexplained by the tested mode-collapse description.'],
        measurements=['Heading persistence and local size exponent, not migration persistence.','Takeover probabilities are per-model, finite-horizon certified-subset quantities.','Existing exemplar Fourier spectra and exact own-mask diagnostic; targeted mode-only independent histories retain unchanged native equations.'],
        limitations=['Missing takeover angles are bounds, not imputed means.','Exact birth/split population classification remains partly unidentifiable.','Two exemplar seeds cannot validate universal mode collapse.','Strong-input mechanism remains unexplained.'],
        sources=[rel+'/methods_provenance.md',rel+'/scientific_readme.md'],equations=equations,main_pdf=rel+'/figure4_main.pdf',supplement_pdf=rel+'/figure4_supplement.pdf',reader=rel+'/index.html')
    if decision['supports_mode_explanation']:
        section['assets'][0]['native_pdf_page']=rel+'/figure4_A_I.pdf'
        for letter in ('J','K'):
            section['assets'].append(dict(path=rel+f'/figure4_main_{letter}.png',native_pdf_page=rel+f'/supporting/vectors/main_{letter}.pdf',
              caption='Detector-free mode test '+letter,page_size=[11.7,6],page_title='',minimum_printed_width_inches=11.7))
        section['narrative']=['Larger rings preserve heading longer while their local size exponent weakens as input strengthens. Independent mode-only histories support the declared collapse and local-slope tests in the predeclared transition-and-stronger input stratum. Region-based closure and remainder remain visible.']
        section['limitations']=[x for x in section['limitations'] if x!='Strong-input mechanism remains unexplained.']
    write(OUT/'legend.json',section)
    (OUT/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Figure4 round6</title></head><body><h1>Enlargement and polarity heading</h1><p>Data generation: '+GENERATED+'</p><p>'+decision['decision']+'</p><p><a href="figure4_main.pdf">Main PDF</a> · <a href="figure4_supplement.pdf">Supplement PDF</a> · <a href="scientific_readme.md">Scientific readme</a> · <a href="methods_provenance.md">Methods</a> · <a href="legend.md">Legend</a></p><img src="figure4_main.png" style="max-width:100%" alt="Figure4 A through I"></body></html>')
    write(META/'provenance.json',dict(generated_utc=GENERATED,baseline_outline='manuscript_outline/outline_v11.pdf',
       source_manifest=str(SOURCE/'manifest.json'),source_manifest_sha256=digest(SOURCE/'manifest.json'),round5=str(OLD),
       scripts={str(p.relative_to(ROOT)):digest(p) for p in sorted((OUT/'scripts').glob('*.py'))},
       rebuild_sha256=digest(Path(__file__)),rebuild_command=shlex_command(),phase_C_max_abs_difference=phase_max,
       exact_mask_lower_signed_remainder_count=negative,exact_mask_compared_seed_conditions=len(comparisons)))
    write(META/'artifact_hashes.json',{p.name:digest(p) for p in [OUT/'figure4_main.pdf',OUT/'figure4_main.png',OUT/'figure4_supplement.pdf',OUT/'legend.md']})
    snapshot=META/'targeted_campaign';snapshot.mkdir(exist_ok=True)
    for source in [RUN/'manifest.json',RUN/'precision/status.json']+list((RUN/'submissions').glob('*.json')):
        if source.exists():
            name=source.name if source.parent==RUN else source.parent.name+'_'+source.name
            shutil.copyfile(source,snapshot/name)

def shlex_command():return PYTHON+' -B '+str(Path(__file__).resolve())

if __name__=='__main__':
    closure=stage_sources();rates,exemplar_rates,decision=collect_modes();main_figure(closure);supplement(rates,exemplar_rates);documents(decision)
    if decision['supports_mode_explanation']:
        from pypdf import PdfWriter
        shutil.copyfile(OUT/'figure4_main.pdf',OUT/'figure4_A_I.pdf')
        writer=PdfWriter();writer.append(OUT/'figure4_main.pdf')
        for letter in ('J','K'):writer.append(OUT/'supporting/vectors'/f'main_{letter}.pdf')
        temp=OUT/'figure4_main.tmp.pdf'
        with temp.open('wb') as f:writer.write(f)
        temp.replace(OUT/'figure4_main.pdf')
        write(META/'artifact_hashes.json',{p.name:digest(p) for p in [OUT/'figure4_main.pdf',OUT/'figure4_main.png',OUT/'figure4_supplement.pdf',OUT/'legend.md']})
    subprocess.run([PYTHON,'-B',str(OUT/'scripts/validate.py')],check=True)
    print(json.dumps(dict(main=str(OUT/'figure4_main.pdf'),supplement=str(OUT/'figure4_supplement.pdf'),decision=decision),indent=2))
