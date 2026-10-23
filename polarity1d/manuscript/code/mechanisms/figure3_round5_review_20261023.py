"""Reviewer-facing Figure 3 derivative: shared translation and explicit size test.
Reuses all completed deterministic/noise-derived scalar data. No integration.
The reviewed 49754ab main/supplement remain immutable historical artifacts.
"""
from pathlib import Path
import os,sys,shutil,json,html,textwrap
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader,PdfWriter
import figure3_round5_render_20261023 as base

OUT=base.OUT;DATA=base.DATA;META=base.META;ROOT=base.ROOT;RUN=base.RUN
HISTORY=OUT/'supporting/historical/reviewed_49754ab'
STEM='figure3_round5_review'
MODELS=base.MODELS;NAMES=base.NAMES;COLORS=base.COLORS;MARKERS=['o','s','^','v','D','P']


def archive_reviewed():
    """Capture before replacing prose; all reviewed main/supp files stay put."""
    HISTORY.mkdir(parents=True,exist_ok=True)
    names=['legend.json','legend.md','README.md','index.html','METHODS_RESULTS.md']
    for name in names:
        target=HISTORY/name
        if not target.exists():shutil.copy2(OUT/name,target)
    keep=[OUT/'figure3_round5.pdf',OUT/'figure3_round5.png',OUT/'figure3_round5.svg',OUT/'figure3_round5_supplement.pdf',OUT/'figure3_round5_new_supplement.pdf']
    records=[dict(path=base.rel(p),sha256=base.edge.sha(p),bytes=p.stat().st_size) for p in keep]
    scientific=[DATA/name for name in ('edge_prediction_comparison.csv','pulse_response_summary.csv','edge_response_curves.csv.gz','kernel_pulse_validation.csv','edge_response_size_folds.csv','edge_independence_summary.csv','edge_MSD_curves.csv')]
    receipt=META/'reviewed_49754ab_preservation.json'
    record=dict(generated_utc=base.edge.now(),reviewed_commit='49754ab',reviewed_main_and_supplement=records,unchanged_scientific_data=[dict(path=base.rel(p),sha256=base.edge.sha(p)) for p in scientific])
    if receipt.exists():
        original=base.read(receipt)
        assert original['reviewed_main_and_supplement']==records
        assert original['unchanged_scientific_data']==record['unchanged_scientific_data']
    else:base.write(receipt,record)
    return record


def normalized_data(frame,pulses):
    drift=pd.read_csv(base.DRIFT);drift=drift[(drift.flicker_frames==5)&(drift.MSD_range_fraction==.25)]
    records=[]
    for model in MODELS:
        group=frame[frame.model==model];minimum=group['size'].min()
        for strength,q in group.groupby('fraction'):
            reference=q[q['size']==minimum].iloc[0]
            for r in q.itertuples():
                source=drift[(drift.model==model)&(drift.C_label==r.size)&(drift.f==strength)].iloc[0]
                pred_ref=reference.Ds_predicted_full_covariance
                defined=bool(r.defined and reference.defined and np.isfinite(pred_ref) and pred_ref>0)
                records.append(dict(model=model,size=r.size,fraction=strength,reference_size=minimum,observed_ratio=r.observed_Ds/reference.observed_Ds,observed_ratio_ci_lower=source.Ds_over_smallest_ci_lower,observed_ratio_ci_upper=source.Ds_over_smallest_ci_upper,shared_translation_predicted_ratio=r.Ds_predicted_full_covariance/pred_ref if defined else np.nan,independent_predicted_ratio=r.Ds_predicted_independent/reference.Ds_predicted_independent if defined else np.nan,defined=defined,undefined_reason='' if defined else ('Smallest-size mean state lacks a polarized boundary' if not reference.defined else r.reason),observed_Ds=r.observed_Ds,observed_reference_Ds=reference.observed_Ds,shared_translation_predicted_Ds=r.Ds_predicted_full_covariance,shared_translation_predicted_reference_Ds=pred_ref,valid_fraction=r.valid_fraction,multiple_front_fraction=r.multiple_front_fraction,n_seeds=r.n_seeds))
    normal=pd.DataFrame(records);normal.to_csv(DATA/'size_normalized_prediction_comparison.csv',index=False)
    q=pulses[(pulses.fraction==.001)&(pulses.pulse_edge==1)&(pulses.tracked_edge==1)&(pulses.multiplier==1)].copy()
    q['integrated_response_over_smallest']=np.nan
    for model in MODELS:
        sel=q.model==model;reference=q[sel].sort_values('size').iloc[0].integrated_displacement_per_pulse_dose
        q.loc[sel,'integrated_response_over_smallest']=q.loc[sel,'integrated_displacement_per_pulse_dose']/reference
    q.to_csv(DATA/'normalized_integrated_edge_response.csv',index=False)
    return normal,q


def main_figure(frame,curves,normal,areas):
    base.setup();fig=plt.figure(figsize=(16.2,15.5))
    ga=fig.add_gridspec(4,2,left=.055,right=.345,bottom=.645,top=.936,wspace=.28,hspace=.20)
    gb=fig.add_gridspec(1,4,left=.425,right=.955,bottom=.727,top=.936,wspace=.57)
    gc=fig.add_gridspec(1,4,left=.07,right=.955,bottom=.433,top=.575,wspace=.40)
    gd=fig.add_gridspec(1,4,left=.07,right=.955,bottom=.263,top=.378,wspace=.40)
    gn=fig.add_gridspec(1,4,left=.07,right=.955,bottom=.077,top=.203,wspace=.40)
    base.panel_A(fig,ga);base.panel_B(fig,gb)
    for label,x,y in [('A',.017,.967),('B',.385,.967),('C',.025,.607),('D',.025,.408)]:fig.text(x,y,label,fontsize=13,fontweight='bold',color='black')
    fig.legend(handles=[Patch(facecolor='.87',label='cue support'),Line2D([],[],marker='D',mfc='white',mec='black',ls='none',label='saddle'),Line2D([],[],marker='|',color='black',ls='none',label='heading')],loc='upper left',bbox_to_anchor=(.052,.986),frameon=False,ncol=3,fontsize=7,handlelength=1)
    fig.legend(handles=[Patch(facecolor='.92',label='heading undefined at ≥1 size'),Line2D([],[],ls=':',color='.25',label='first saddle at 50% (C style)'),Line2D([],[],marker='^',mfc='white',mec='black',ls='none',label='time lower bound'),Line2D([],[],marker='o',mfc=COLORS['wave_pinning'],mec='black',ls='none',label='i–ii: A placements')],loc='upper left',bbox_to_anchor=(.421,.701),frameon=False,ncol=2,fontsize=7.4,columnspacing=1.3,handlelength=1.6)
    fig.text(.425,.651,'B: C = 30, 37.8, 47.6; Spring 2.5, 3.15, 3.97.\nC–D: chemical C = 20–60; Spring 1–3. These are not matched size grids.',fontsize=8,color='black',va='top')
    for j,model in enumerate(MODELS):
        ax=fig.add_subplot(gc[0,j]);ax.set_title(NAMES[model],pad=22);tau=5 if model=='spring_ring' else 50;ax.axvspan(0,tau,color='.94',zorder=0)
        q=curves[(curves.model==model)&(curves.pulse_edge==1)&(curves.tracked_edge==1)&(curves.multiplier==1)]
        sizes=sorted(q['size'].unique());handles=[]
        for i,size in enumerate(sizes):
            z=q[q['size']==size];col=base.color(model,i);ax.plot(z.time,z.displacement_per_amplitude,color=col,lw=1.35);handles.append(Line2D([],[],color=col,label=f'{size:g}'))
        ax.set_ylim(0,q.displacement_per_amplitude.max()*2.05);ax.set_xlabel('Time from pulse onset');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3));ax.tick_params(length=3,pad=2)
        if j==0:ax.set_ylabel(r'Pulsed-edge displacement / $a$')
        ax.text(-.01,1.063,'C',transform=ax.transAxes,ha='right',fontsize=6.5)
        ax.legend(handles=handles,frameon=False,fontsize=6.4,ncol=5,loc='upper left',bbox_to_anchor=(0,1.105),handlelength=.85,columnspacing=.7,handletextpad=.3,labelspacing=.2,borderpad=.2)
        ins=ax.inset_axes([.50,.53,.47,.43]);g=areas[areas.model==model].sort_values('size')
        ins.axhline(1,color='.8',lw=.55,zorder=0);ins.plot(g['size'],g.integrated_response_over_smallest,color=COLORS[model],lw=1,marker='o',ms=2.5);ins.set_ylim(0,1.1);ins.set_yticks([0,.5,1]);ins.tick_params(labelsize=5.7,length=2,pad=1);ins.set_xlabel('C',fontsize=6,labelpad=1);ins.set_ylabel(r'$\mathcal{A}/\mathcal{A}_{\min}$',fontsize=7,labelpad=1)
        # D(i): absolute physical diffusivity, retaining the rejected candidate.
        ax=fig.add_subplot(gd[0,j]);ax.set_title(NAMES[model]);ax.set_xscale('log');ax.set_yscale('log');g=frame[(frame.model==model)&frame.defined&(frame.observed_Ds>0)]
        lo=min(g.Ds_predicted_independent.min(),g.observed_Ds.min())*.4;hi=max(g.Ds_predicted_full_covariance.max(),g.observed_Ds.max())*2
        ax.plot([lo,hi],[lo,hi],color='.35',lw=.9,zorder=0)
        for i,size in enumerate(sorted(frame[frame.model==model]['size'].unique())):
            for k,strength in enumerate(base.edge.FRACTIONS):
                z=g[(g['size']==size)&(g.fraction==strength)]
                if not len(z):continue
                r=z.iloc[0];col=base.color(model,i)
                for name,fill in [('Ds_predicted_independent','white'),('Ds_predicted_full_covariance',col)]:
                    y=r[name]
                    if pd.notna(r.observed_Ds_ci_lower) and r.observed_Ds_ci_lower>0:ax.plot([r.observed_Ds_ci_lower,r.observed_Ds_ci_upper],[y,y],color=col,lw=.5,zorder=1)
                    ax.plot(r.observed_Ds,y,marker=MARKERS[k],ms=4.25,mfc=fill,mec=col if fill=='white' else COLORS[model],mew=.65 if fill=='white' else .3,ls='none',zorder=2 if fill=='white' else 3)
        ax.set_xlim(lo,hi);ax.set_ylim(lo,hi);ax.set_xlabel(r'Measured $D_s$ (Fig. 4)');ax.tick_params(length=3,pad=2)
        if j==0:ax.set_ylabel(r'Predicted $D_s$')
        if model=='goryachev':ax.text(.02,.98,'2 mean states undefined',transform=ax.transAxes,va='top',fontsize=6.5)
        # D(ii): the size comparison, normalized within each f and model.
        ax=fig.add_subplot(gn[0,j]);ax.set_title(NAMES[model]);g=normal[(normal.model==model)&normal.defined]
        if model=='goryachev':
            ax.set_xscale('log');ax.set_yscale('log');low=min(g.observed_ratio.min(),g.shared_translation_predicted_ratio.min())*.5;high=1.5;ax.set_xlim(low,high);ax.set_ylim(low,high)
        else:low=0.;high=1.55 if model=='spring_ring' else 1.2;ax.set_xlim(low,high);ax.set_ylim(low,high);ax.set_xticks([0,.5,1] if high<1.3 else [0,.5,1,1.5]);ax.set_yticks([0,.5,1] if high<1.3 else [0,.5,1,1.5])
        ax.plot([low if low else 0,high],[low if low else 0,high],color='.35',lw=.9,zorder=0)
        for k,strength in enumerate(base.edge.FRACTIONS):
            z=g[g.fraction==strength].sort_values('size');ax.plot(z.observed_ratio,z.shared_translation_predicted_ratio,color='.83',lw=.45,zorder=0)
        for i,size in enumerate(sorted(frame[frame.model==model]['size'].unique())):
            for k,strength in enumerate(base.edge.FRACTIONS):
                z=g[(g['size']==size)&(g.fraction==strength)]
                if not len(z):continue
                r=z.iloc[0];col=base.color(model,i)
                if pd.notna(r.observed_ratio_ci_lower) and r.observed_ratio_ci_lower>0:ax.plot([r.observed_ratio_ci_lower,r.observed_ratio_ci_upper],[r.shared_translation_predicted_ratio]*2,color=col,lw=.5,zorder=1)
                ax.plot(r.observed_ratio,r.shared_translation_predicted_ratio,marker=MARKERS[k],ms=4.5,mfc=col,mec=COLORS[model],mew=.35,ls='none',zorder=3)
        ax.set_xlabel(r'Measured $D_s(C)/D_s(C_{\min})$');ax.tick_params(length=3,pad=2)
        if j==0:
            ax.set_ylabel(r'Predicted $D_s(C)/D_s(C_{\min})$')
            r=g[(g['size']==60)&(g.fraction==.001)].iloc[0];ax.plot(r.observed_ratio,r.shared_translation_predicted_ratio,'o',ms=8,mfc='none',mec='black',mew=.8,zorder=4)
            ax.text(.035,.97,'f = 0.001, C60/C20\n0.49 predicted; 0.47 measured',transform=ax.transAxes,va='top',fontsize=7.2)
        if model=='goryachev':ax.text(.025,.975,'f = 0.01 reference undefined',transform=ax.transAxes,va='top',fontsize=6.5)
    fig.text(.07,.405,r'C: $f=0.001$; $a=10^{-4}$ × native reference. Insets: finite integrated displacement, normalized to the smallest size (shared 0–1.1 scale).',fontsize=8)
    fig.text(.07,.226,'D: above, absolute diffusivity; below, size ratios at the same input. Lines: 1:1. Horizontal bars: whole-seed-bootstrap 95% CI.',fontsize=8)
    fig.legend(handles=[Line2D([],[],marker='o',mfc='.5',mec='.25',ls='none',label='shared translation, no fit'),Line2D([],[],marker='o',mfc='white',mec='.5',ls='none',label='independent edges (rejected)')],loc='lower left',bbox_to_anchor=(.067,.034),ncol=2,frameon=False,fontsize=8,handlelength=1.1,columnspacing=2)
    fig.legend(handles=[Line2D([],[],marker=m,mfc='.65',mec='.25',ls='none',label=f'{f:g}') for f,m in zip(base.edge.FRACTIONS,MARKERS)],title='f',title_fontsize=7.5,loc='lower left',bbox_to_anchor=(.585,.022),ncol=3,frameon=False,fontsize=7.1,handlelength=1,columnspacing=1.1)
    for ext in ('pdf','png','svg'):fig.savefig(OUT/f'{STEM}.{ext}',dpi=220)
    plt.close(fig)
    recipe=dict(generated_utc=base.edge.now(),figure_size=[16.2,15.5],D_absolute=dict(filled='Ds_predicted_full_covariance',open='Ds_predicted_independent',candidate_status='tested and rejected',parameter_fits=0),D_size_normalized=dict(source='size_normalized_prediction_comparison.csv',reference='same model and f, fixed smallest grid size',filled='shared_translation_predicted_ratio',rows=int(normal.defined.sum()),undefined_reference_rows=int((~normal.defined).sum()),WP_f001_C60_over_C20=dict(predicted=.4939004889085995,observed=.47324320218752636)),C_inset=dict(source='normalized_integrated_edge_response.csv',y_limits=[0,1.1],normalization='same model, f=.001, smallest size, finite area20tau'))
    base.write(META/'review_render_recipe.json',recipe)


def supplement(frame,normal):
    base.setup();new=OUT/f'{STEM}_new_supplement.pdf'
    with PdfPages(new) as pdf:
        base.page_text(pdf,'Shared boundary translation and the explicit size comparison',[
            'The deterministic data, measured drift estimates, pulse-gain calibration and exact forcing covariance are unchanged. Filled main-panel D markers now retain the shared neutral translation; open markers preserve the tested independent-edge candidate. No correction is fitted to measured drift. Near-unity late-lag boundary co-motion refutes independence in Wave-pinning, Holmes and Spring.',
            'The second D row compares predicted and measured size ratios at each fixed input. The same model-specific smallest size is used in both denominators. This removes the dominant f-squared scaling from the display. In Wave-pinning at f=0.001, the C60/C20 ratio is 0.4939 predicted and 0.4732 measured; the measured terminal susceptibility ratio is 0.7007. This is a concrete size test, not an inference from the absolute five-decade span.',
            'C insets show finite integrated displacement divided by the smallest-size value, on identical zero-based 0–1.1 axes. The terminal gain falls strongly in Goryachev and moderately in Wave-pinning, while Holmes and Spring change much less. The two perturbation assays are distinct: varying-offset continuous photoactivation tests sliding versus takeover and half-turn timing; a weak brief boundary pulse tests small-signal boundary displacement and common translation. The edge assay does not uniquely identify the Figure 2 mechanism classes.',
            'Goryachev C20 and C30 at f=0.01 lack a stationary polarized boundary at the mean input. Their absolute predictions are undefined. All five Goryachev size-normalized predictions at that f are undefined because the fixed C20 reference is undefined; a larger size is not substituted. Goryachev disagreement remains a result, and physical tracked domain drift does not explain its entire heading decay.',
            'The historical reviewed main page and complete 25-page reviewed supplement follow these current supplementary pages. They remain dated records of the independent-edge presentation at commit 49754ab. The current main and legend supersede their interpretation, while all numerical source values and scientific limitations are retained.'])
        fig,axes=plt.subplots(2,4,figsize=(14,7.8))
        for j,model in enumerate(MODELS):
            g=frame[(frame.model==model)&frame.defined&(frame.observed_Ds>0)]
            for i,size in enumerate(sorted(g['size'].unique())):
                q=g[g['size']==size].sort_values('fraction');col=base.color(model,i)
                axes[0,j].plot(q.fraction,q.full_covariance_over_observation,color=col,marker='o',ms=3,lw=1)
            axes[0,j].axhline(1,color='.5',ls=':',lw=.8);axes[0,j].set_xscale('log');axes[0,j].set_yscale('log');axes[0,j].set_title(NAMES[model]);axes[0,j].set_xlabel('f')
            q=normal[(normal.model==model)&normal.defined]
            for strength in base.edge.FRACTIONS:
                z=q[q.fraction==strength].sort_values('size');axes[1,j].plot(z['size'],z.shared_translation_predicted_ratio,color=COLORS[model],alpha=.5,lw=.9);axes[1,j].plot(z['size'],z.observed_ratio,color=COLORS[model],alpha=.5,lw=.9,ls='--')
            axes[1,j].set_xlabel('C');axes[1,j].set_ylim(bottom=0)
        axes[0,0].set_ylabel('Full prediction / observed drift');axes[1,0].set_ylabel('Size ratio at fixed input');fig.subplots_adjust(left=.07,right=.98,bottom=.10,top=.90,wspace=.36,hspace=.42);fig.text(.05,.96,'Absolute agreement and size dependence; no fitted correction',fontsize=12,fontweight='bold');fig.text(.07,.025,'Bottom: solid prediction, dashed observation, one curve per input. Source tables retain condition-specific values, missing references and exposure.',fontsize=9);pdf.savefig(fig);plt.close(fig)
    reviewed=[OUT/'figure3_round5.pdf',OUT/'figure3_round5_supplement.pdf']
    w=PdfWriter();w.append(new)
    for p in reviewed:w.append(p)
    target=OUT/f'{STEM}_supplement.pdf';w.write(target)
    base.write(META/'review_supplement_preservation.json',dict(generated_utc=base.edge.now(),new_pages=len(PdfReader(new).pages),preserved=[dict(path=base.rel(p),sha256=base.edge.sha(p),pages=len(PdfReader(p).pages)) for p in reviewed],total_pages=len(PdfReader(target).pages)))


def documents(frame,normal,areas):
    stats=base.read(META/'results.json')['stats'];prior=base.read(HISTORY/'legend.json')
    title='Enlargement slows both turning routes and selectively suppresses boundary susceptibility.'
    legend=[
        'A. Wave-pinning activity at C = 30 for a cue straddling the incumbent boundary (i, gap −1) or separated from it (ii, gap 4). Rows show 0%, 50% and 90% of cue-direction heading, then the relaxed cue-free profile. Grey: cue support; diamonds: inter-region saddles; ticks: heading.',
        'B. Half-turn time versus cue-edge gap from the incumbent’s 10%-contrast boundary. Increasing sizes use lighter solid, intermediate dashed and darker dotted curves; chemical C = 30, 37.8, 47.6, Spring 2.5, 3.15, 3.97. Vertical lines mark the first resolved two-region saddle at half-turn; grey bands mark undefined heading at ≥1 size; triangles are time lower bounds. Enlargement slows both turning routes. C–D use different size grids, so these panels do not establish a matched-size relation between strong-cue turning and weak-pulse response.',
        'C. Control-subtracted displacement of the pulsed half-contrast boundary divided by peak input a, at f = 0.001. A Gaussian pulse has SD ℓ = 1.5 and duration τ = 50 (Spring 0.15 and 5); peak a = 10⁻⁴ of the native input reference. Shades denote chemical C = 20–60 or Spring resting C = 1–3. Insets show displacement integrated to 20τ, divided by pulse dose aτ, then normalized to the smallest size on a common 0–1.1 axis. Largest/smallest terminal gains are 0.701 Wave-pinning, 0.914 Holmes, 0.989 Spring and 0.0878 Goryachev. Both edges and half/double amplitudes are checked in the supplement.',
        'D. Above: predicted versus Figure 4 physical tracked drift for f ≤ 0.01. Filled markers retain shared translation; open markers retain the rejected independent-edge candidate. Measured pulse gain sets response amplitude; the unchanged native linearization supplies spatial shape, and the exact input covariance supplies forcing statistics. Nothing is fitted to drift. Below: predicted versus measured Dₛ(C)/Dₛ(Cmin) at the same f, removing the dominant input-strength scaling. Wave-pinning at f = 0.001 gives C60/C20 = 0.49 predicted and 0.47 measured. Lines: 1:1; bars: whole-seed-bootstrap 95% CI. Goryachev C20/C30 at f = 0.01 have no polarized mean-state boundary; the normalized prediction at that f also lacks its C20 reference.',
        'Wave-pinning, Holmes and Spring have shared-translation prediction/observation medians 1.021, 1.080 and 1.069. Resolved quarter-window edge correlations are 0.980, 0.983 and 0.997 in the two exemplar seeds, refuting edge independence. Goryachev retains large quantitative disagreement. Tracked drift is conditional on trackable domains and is distinct from total heading decay.',
        'Two cellular assays follow: vary photoactivation offset from straddling to separated placement and measure half-turn time; separately apply a weak brief boundary pulse and track both boundaries through withdrawal. The second assay measures susceptibility and common translation; it does not uniquely assign Figure 2 mechanism classes. Native units are not calibrated across models.']
    narrative=[
        'Cue placement separates sliding of the incumbent from takeover by a new region; enlargement slows both routes. Weak boundary susceptibility changes differently across models: terminal response at the largest size is 70% of the smallest-size response in Wave-pinning, 9% in Goryachev, 91% in Holmes and 99% in Spring.',
        'The deterministic boundary response predicts physical domain drift without fitting when shared translation is retained. Wave-pinning’s C60/C20 drift ratio at f = 0.001 is 0.494 predicted and 0.473 measured. Holmes and Spring show much weaker susceptibility changes. The independent-edge approximation is refuted by boundary co-motion; Goryachev’s drift prediction remains quantitatively inconsistent and does not explain its full heading persistence.'
    ]
    section=dict(id='figure3',title=title,status='round5_review_complete',legend=legend,role='Size-dependent susceptibility links deterministic boundary perturbation to weak-input domain drift; cue position separates two slower turning routes in enlarged rings.',narrative=narrative,measurements=['Vary photoactivation offset from straddling to separated placement; measure time to half-turn.','Apply a weak brief pulse across the leading-edge boundary; measure both boundary displacements, common translation and covariance. This readout complements rather than uniquely identifies the Figure 2 mechanism classes.'],limitations=['Panel B and C–D use different size grids and perturbation strengths; a matched-size comparison is not established. One pulse does not identify spatial response shape; native linearization supplies that explicit assumption, while measured pulse gain supplies amplitude. Two float16 exemplar seeds support the covariance check, with short-lag precision limits. Goryachev’s unresolved mean-state boundaries and quantitative drift disagreement remain visible. Physical domain drift is not total heading decay or empirical cell migration persistence.'],equations=[dict(math=r'D_s^{\mathrm{ind}}=\frac{D_L+D_R}{4},\qquad D_s=\frac{D_L+D_R+2D_{LR}}{4}',number='3.1',description='Open markers test independent edges; filled markers retain the untuned shared-translation covariance.'),dict(math=r'\widehat g_{e,j}=\frac{M_e}{Q_e}g_j,\qquad M_e=\frac{\Delta x_e(20\tau)}{a\tau},\quad Q_e=\sum_jg_jh_{e,j}',number='3.2',description='Measured terminal displacement per pulse dose supplies amplitude; the unchanged native response supplies spatial shape. The covariance derivation is in the methods and supplement.')],assets=[dict(path=base.rel(OUT/f'{STEM}.png'),native_pdf_page=base.rel(OUT/f'{STEM}.pdf'),page_size=[16.2,15.5],minimum_printed_width_inches=16.2,page_title='')],main_pdf=base.rel(OUT/f'{STEM}.pdf'),supplement_pdf=base.rel(OUT/f'{STEM}_supplement.pdf'),reader=base.rel(OUT/'index.html'),sources=[dict(label='Physical diffusion prediction and observed drift',path=base.rel(DATA/'edge_prediction_comparison.csv')),dict(label='Size-normalized prediction and observed drift',path=base.rel(DATA/'size_normalized_prediction_comparison.csv')),dict(label='Normalized integrated boundary response',path=base.rel(DATA/'normalized_integrated_edge_response.csv')),dict(label='Methods and covariance derivation',path=base.rel(OUT/'METHODS_RESULTS_review.md')),dict(label='Both boundary assays and numerical checks',path=base.rel(META/'numerical_checks.json')),dict(label='Reviewed main/supplement preservation',path=base.rel(META/'reviewed_49754ab_preservation.json'))])
    base.write(OUT/'legend.json',section);(OUT/'legend.md').write_text('**Figure 3. '+title+'**\n\n'+'\n\n'.join(legend)+'\n')
    # Retain the full tested derivation, updating only which untuned prediction
    # is prominent. Relaxation details and exact covariance stay here.
    methods=(HISTORY/'METHODS_RESULTS.md').read_text()
    first=methods.index('## Derivation')
    methods='# Figure 3: size-dependent boundary response and shared translation\n\nGenerated '+base.edge.now()+'. All native producers finished; this derivative uses existing data only.\n\n'+'\n\n'.join(narrative)+'\n\n## Assays and size matching\n\nThe nominal weak-pulse peaks are 10⁻⁵ Wave-pinning, 10⁻⁴ Holmes, 4 × 10⁻⁴ Spring and 6.4 × 10⁻⁴ Goryachev, each 10⁻⁴ of its model’s fixed native input reference. Gaussian SD and pulse duration are 1.5 and 50 chemical units, or 0.15 and 5 Spring units. Each starts from the stationary polarized state under the same condition’s mean input, without noise.\n\nA–B use finite-strength continuous cues placed across or away from the incumbent boundary, measuring time to half-turn. C uses a separate weak brief boundary pulse. B uses chemical sizes 30/37.8/47.6 (Spring 2.5/3.15/3.97), whereas C–D use 20/30/40/50/60 (Spring 1/1.5/2/2.5/3). A matched-size link between strong-cue reorientation and small-signal susceptibility is therefore not established; the two assays are complementary and do not uniquely sort the Figure 2 mechanism classes.\n\n'+methods[first:]
    methods=methods.replace('Keeping covariance is a separate, untuned sensitivity; we retain the failing independence prediction in the main panel.','The main panel shows both the shared-translation prediction and the rejected independent-edge candidate. Retaining the cross term uses the same native model and measured pulse gain; it is not a fitted factor-of-two correction.')
    methods+='''\n## Fixed-input size test\n\nFor each model and input, both observed and predicted diffusivity are divided by their own value at that model’s fixed smallest grid size. The same reference size is retained at every input. The resulting 115 defined prediction pairs are exported with observed bootstrap ratio intervals in `supporting/data/size_normalized_prediction_comparison.csv`. Five Goryachev ratios at f=0.01 are unavailable because the C20 mean-state response is undefined; C40 is not substituted. At f=0.001, Wave-pinning’s C60/C20 ratio is 0.49390049 predicted and 0.47324320 measured. The size-normalized candidate ratios coincide to numerical precision because the rejected independence assumption changes the absolute diffusivity by essentially one common factor; therefore the size subrow shows the covariance-complete prediction, while the absolute subrow retains both candidates.\n\nC insets use the finite integrated response through 20τ divided by the smallest-size response, with the same zero-based 0–1.1 range for all models. Their source is `supporting/data/normalized_integrated_edge_response.csv`. These insets compare integrated displacement; the reported 0.701/0.914/0.989/0.0878 folds compare terminal displacement per pulse dose. They are distinct summaries.\n\n## Relaxation row and historical evidence\n\nFor A’s cue-free row, the cue is switched off at 95% heading or at 20,000 if 95% is not reached. Cue-free relaxation ends at the native stationarity tolerance or at four times the cue-on duration. Example i withdraws at 790.5 and relaxes 258.6; example ii withdraws at 227.6 and relaxes 291.3. Both reach the stationarity criterion 2.2683×10⁻⁷ in active-species RHS. Profiles and all A/B numerical coordinates remain unchanged.\n\nThe reviewed main and complete supplement at commit 49754ab are preserved unchanged and appended after the new supplementary pages; their earlier single-candidate interpretation is historical. The current main exposes both candidates and the fixed-input size test. All covariance integration, conservation and amplitude checks remain valid because no scientific numerical source changed.\n'''
    (OUT/'METHODS_RESULTS_review.md').write_text(methods)
    (OUT/'README.md').write_text('# Figure3: size-dependent boundary susceptibility and drift\n\n[Main PDF]('+STEM+'.pdf) · [supplement]('+STEM+'_supplement.pdf) · [legend](legend.md) · [methods/results](METHODS_RESULTS_review.md) · [reader](index.html).\n\n'+'\n\n'.join(narrative)+'\n\nNo new native integrations, stochastic histories or parameter changes. The main displays the untuned shared-translation prediction and keeps the rejected independent-edge candidate. C insets use a shared normalized zero-based scale; D includes a fixed-input size comparison. The reviewed edition remains intact and is appended in the supplement.\n\nRebuild:\n\n```sh\npython -B polarity1d/manuscript/code/mechanisms/figure3_round5_review_20261023.py\npython -B polarity1d/manuscript/code/mechanisms/figure3_round5_review_checks_20261023.py\n```\n')
    generated=base.edge.now();body=f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure3: boundary susceptibility and size-dependent drift</title><style>body{{max-width:1300px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}}img{{width:100%}}a{{color:#166aa2}}</style><h1>{html.escape(title)}</h1><p>Data generated {generated}; all producers complete. No new native runs.</p><p><a href="{STEM}.pdf">Main PDF</a> · <a href="{STEM}_supplement.pdf">Supplement</a> · <a href="METHODS_RESULTS_review.md">Methods/results</a> · <a href="supporting/data/size_normalized_prediction_comparison.csv">Size comparison</a></p><img src="{STEM}.png" alt="Two turning routes, normalized edge response, shared-translation drift and fixed-input size comparison">'
    body+=''.join('<p>'+html.escape(p)+'</p>' for p in legend);body+='</html>';(OUT/'index.html').write_text(body)
    base.write(META/'review_build_receipt.json',dict(generated_utc=generated,completed=True,source_scientific_data_unchanged=True,native_reruns=0,retuned_parameters=0,main_sha256=base.edge.sha(OUT/f'{STEM}.pdf'),legend_sha256=base.edge.sha(OUT/'legend.json'),main_pages=len(PdfReader(OUT/f'{STEM}.pdf').pages),main_path=base.rel(OUT/f'{STEM}.pdf'),preserved_reviewed_commit='49754ab'))


def main():
    archive_reviewed();frame=pd.read_csv(DATA/'edge_prediction_comparison.csv');pulses=pd.read_csv(DATA/'pulse_response_summary.csv');curves=pd.read_csv(DATA/'edge_response_curves.csv.gz');normal,areas=normalized_data(frame,pulses);main_figure(frame,curves,normal,areas);supplement(frame,normal);documents(frame,normal,areas);print(OUT/f'{STEM}.pdf')
if __name__=='__main__':main()
