#!/home/4473331/.conda/envs/scaling_sandbox/bin/python
"""Presentation-only supplement renderer. Reads the 62 frozen numeric exports.

No estimation, pooling, bootstrap, simulations, or numeric file writes occur.
Plot helpers are imported from the immutable executed-analysis source snapshot.
Historical K_/L_ CSV filenames remain unchanged; rendered labels are S4.7/S4.8.
"""
import importlib.util
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import common
spec=importlib.util.spec_from_file_location('_executed_analysis_plot_helpers',common.META/'final_analysis_source_snapshot/rebuild_figure4.py')
helpers=importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)
for key,value in vars(helpers).items():
    if not key.startswith('__'):globals()[key]=value

def render_cached_supplement():
    exemplar_rates=csvread(DATA/'exemplar_rates.csv')
    temporary=OUT/'figure4_supplement.tmp.pdf'
    with PdfPages(temporary) as pdf:
        def page(title,caption='',rows=1,height=4.0):
            fig,axes=plt.subplots(rows,4,figsize=(11.7,height),squeeze=False)
            fig.subplots_adjust(left=.07,right=.975,bottom=.18,top=.83,wspace=.38,hspace=.65)
            fig.text(.055,.965,title,ha='left',va='top',fontsize=12,color='black')
            if caption:fig.text(.06,.035,caption,ha='left',va='bottom',fontsize=7,wrap=True)
            return fig,axes
        def save(fig):
            pdf.savefig(fig)
            plt.close(fig)
        birth=readold('birth_heading_condition_metrics');fig,axes=page('S4.1  Certified persistent-secondary birth rates',
           'Rates are certified lower bounds per actual perimeter–time. Uncertified onsets can contain split descendants; native local parameters are unchanged.')
        for ax,m in zip(axes[0],MODELS):lines(ax,birth,m,'certified_birth_rate_lower_bound',True);ax.set_yscale('symlog',linthresh=1e-8);ax.set_ylim(bottom=0)
        axes[0,0].set_ylabel('Certified births / perimeter / time');save(fig)
        fig,axes=page('S4.1b  Native heading correlations across all five sizes',
            'Transition inputs match main C–D. Whole-history 95% intervals; light to dark denotes increasing circumference. Endpoint coherence |P|≥0.1. All points are in condition_lag_metrics.csv.')
        for ax,m in zip(axes[0],MODELS):
            chosen=num(next(r for r in readold('model_D_selection') if r['model']==m),'f')
            for size in sizes(m):
                rr=sorted([r for r in primary(readold('condition_lag_metrics')) if r['model']==m and num(r,'f')==chosen and num(r,'C_label')==size],key=lambda r:num(r,'lag'))
                x=[num(r,'lag')/1000 for r in rr];co=color(m,size)
                ax.plot(x,[num(r,'C') for r in rr],color=co,lw=.8)
                ax.fill_between(x,[num(r,'C_ci_lower') for r in rr],[num(r,'C_ci_upper') for r in rr],color=co,alpha=.1,lw=0)
            ax.set_xlabel(r'Lag ($\times10^3$)');ax.set_ylim(-.08,1.05);style(ax,m)
        axes[0,0].set_ylabel('Heading correlation');save(fig)
        heading=readold('heading_amplitude_condition_metrics');fig,axes=page('S4.2  Motion, amount and membership shares',
            'Native-amplitude symmetric Shapley convention. Signed fractions can be negative or exceed one; Eligibility requires coherent native heading at both endpoints.',rows=3,height=8.5)
        for i,key in enumerate(['motion_fraction','coexisting_amount_fraction','birth_loss_fraction']):
            for ax,m in zip(axes[i],MODELS):lines(ax,heading,m,key,True);ax.axhline(0,color='.7',lw=.5)
            axes[i,0].set_ylabel(['Motion share','Amount share','Membership share'][i])
        save(fig)
        fig,axes=page('S4.3  Pooled finite-horizon takeovers',
            'Pooled f > 0.01 diagnostic, shown only for positive distance. Main J uses a separately selected input per model and >=200 births at every size.')
        takeover_plot(axes[0],readold('takeover_by_size_distance'))
        for ax,m in zip(axes[0],MODELS):ax.set_title(NAMES[m]+', pooled f >0.01',fontsize=8)
        save(fig)
        # All inputs: same plot, four facets; zero/short follow-up exclusions are
        # retained in source CSV and are never treated as unsuccessful births.
        take=readold('takeover_by_condition_distance')
        for binid in (1,2):
            fig,axes=page('S4.4  Takeover probability at every input: '+('0 < d ≤ C/6' if binid==1 else 'd > C/6'),
                'Common follow-up 20τ; each curve is one size. Certified no-overlap subset only. Inputs with no certified births remain missing.')
            for ax,m in zip(axes[0],MODELS):lines(ax,[r for r in take if num(r,'distance_bin')==binid],m,'probability');ax.set_ylim(-.02,1.02)
            axes[0,0].set_ylabel('Takeover probability');save(fig)
        for conv in ['point_matched_mass','exact_full_mask','exact_contour_mask','round5_point_native_amplitude','round5_exact_native_amplitude']:
            rr=csvread(DATA/('exact_mask_curve_'+conv+'.csv'))
            fig,axes=page('S4.5  Own-mask moment check: '+conv.replace('_',' '),
                'Existing seeds 0/1; saved binned float16 profiles. Native-mass denominator in first 3 variants; last 2 use region-total normalization and native amplitude scaling. Shape enters amount/shape in first 3 and motion/shape in last 2. Algebraic reconstruction is not causal validation.')
            for ax,m in zip(axes[0],MODELS):lines(ax,rr,m,'remainder_fraction');ax.axhline(0,color='.7',lw=.5)
            axes[0,0].set_ylabel('Remainder share');save(fig)
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
        # Cached mode points and guides retain the executed numerical analysis.
        points=csvread(DATA/'K_mode_rate_points.csv');guides=csvread(DATA/'K_mode_guides.csv')
        fig,axes=page('S4.7  Fourier-mode decorrelation versus physical wavenumber',
             'Solid weak f = 0.001 (existing exemplars); dashed transition input. Open triangles show inverse-half-time window bounds, conditional on a monotone first crossing. Transition curves use the independent mode-only cohort; source and resolution flags are exported.',height=6)
        for ax,m in zip(axes[0],MODELS):
            for role,ls in (('weak','-'),('transition','--')):
                for size in sizes(m):
                    rr=sorted([r for r in points if r['model']==m and r['input_role']==role and num(r,'C_label')==size],key=lambda r:num(r,'q'))
                    resolved=[num(r,'gamma') if r.get('status')=='resolved' and r.get('half_crossing_observed')=='True' else np.nan for r in rr]
                    forced=bool(rr and rr[0].get('forced_combined')=='True')
                    ax.plot([num(r,'q') for r in rr],resolved,ls=ls,color=color(m,size),lw=.8,marker='o',ms=1.5,mfc='white' if forced else color(m,size),mew=.5)
                    for r in rr:
                        bound=num(r,'crossing_rate_bound');direction=r.get('crossing_bound_direction')
                        if np.isfinite(bound) and bound>0 and direction in ('lower','upper'):
                            ax.plot(num(r,'q'),bound,marker='^' if direction=='lower' else 'v',mfc='white',mec=color(m,size),ms=2.2,mew=.5)
                            ax.annotate('',xy=(num(r,'q'),bound*(1.35 if direction=='lower' else .74)),xytext=(num(r,'q'),bound),arrowprops=dict(arrowstyle='->',lw=.35,color=color(m,size)))
            ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Physical wavenumber q');style(ax,m)
            for slope,ls in ((0,'-.'),(2,':')):
                gg=sorted([r for r in guides if r['model']==m and num(r,'guide_slope')==slope],key=lambda r:num(r,'q'))
                if gg:
                    ax.plot([num(r,'q') for r in gg],[num(r,'gamma') for r in gg],color='.5',lw=.7,ls=ls)
                    ax.text(num(gg[-1],'q'),num(gg[-1],'gamma'),str(slope),fontsize=7,color='.4')
            if not any(r['model']==m for r in guides):ax.set_xlim(.05,10);ax.set_ylim(1e-6,.1)
        axes[0,0].set_ylabel('Decorrelation rate γ');save(fig)
        comp=csvread(DATA/'L_alpha_mode_slope_points.csv');equality=csvread(DATA/'L_equality_guide.csv')
        fig,axes=page('S4.8  Size exponent versus independently predicted local mode slope',
             'Mode slope uses the nearest 8 mode points from other sizes at the same model and input; the target size is excluded. One-sided predictions are open squares. Forced conditions have open circles. This tests collapse rather than replotting the same size slope.',height=6)
        for ax,m in zip(axes[0],MODELS):
            for r in [r for r in comp if r['model']==m]:
                if not np.isfinite(num(r,'alpha')+num(r,'local_mode_slope')):continue
                ax.plot(num(r,'local_mode_slope'),num(r,'alpha'),marker='s' if r.get('extrapolated')=='True' else 'o',mfc='white' if r.get('extrapolated')=='True' or r.get('forced_combined')=='True' else color(m,num(r,'C_label')),mec=color(m,num(r,'C_label')),ms=2.8,mew=.6)
            gg=sorted([r for r in equality if r['model']==m],key=lambda r:num(r,'local_mode_slope'))
            xx=[num(r,'local_mode_slope') for r in gg];yy=[num(r,'alpha') for r in gg]
            ax.plot(xx,yy,color='.6',lw=.7);ax.set_xlim(min(xx),max(xx));ax.set_ylim(min(yy),max(yy));ax.set_xlabel('Local slope of γ(q)');style(ax,m)
        axes[0,0].set_ylabel('Figure 4 local α');save(fig)
        # Phase-only source equivalence is shown as a difference curve, directly
        # identifying finite binning/precision effects.
        phase=csvread(DATA/'exemplar_phase_check.csv');fig,axes=page('S4.9  First-mode phase reproduces native heading correlation',
            'Difference between Fourier a1 phase correlation and native compact P heading, at identical exact saved endpoints. Eligibility uses each stored |P|≥0.1; nonzero differences include threshold crossing and bin/float16 precision.')
        for ax,m in zip(axes[0],MODELS):
            for size in sizes(m):
                rr=[r for r in phase if r['model']==m and num(r,'C_label')==size and np.isfinite(num(r,'difference'))]
                ax.scatter([num(r,'lag') for r in rr],[num(r,'difference') for r in rr],s=2,color=color(m,size),alpha=.5)
            ax.axhline(0,color='.6',lw=.6);ax.set_xscale('log');ax.set_xlim(20 if m=='spring_ring' else 200,2000 if m=='spring_ring' else 20000);ax.set_xlabel('Lag');style(ax,m)
        axes[0,0].set_ylabel('Phase C − native heading C');save(fig)
        fig,axes=page('S4.10  Complex-mode and phase-only first-mode rates',
             'Same cohort, nearest-half lag rule. Complex: amplitude-weighted all endpoints; phase: |P|≥0.1, point estimates without bootstrap precision. Crosses flag phase coverage <0.5. The extreme Spring ratio uses only 4 eligible terminal-lag pairs; all point estimates remain exported.')
        phase_rates=csvread(DATA/'mode1_complex_phase_rate_comparison.csv')
        for ax,m in zip(axes[0],MODELS):
            lines(ax,phase_rates,m,'complex_over_phase');ax.set_yscale('log');ax.axhline(1,color='.55',lw=.6)
            sparse=[r for r in phase_rates if r['model']==m and r.get('phase_coverage_below_original_alpha_gate')=='True' and np.isfinite(num(r,'complex_over_phase'))]
            ax.scatter([num(r,'f') for r in sparse],[num(r,'complex_over_phase') for r in sparse],s=12,marker='x',color='.25',linewidths=.6,zorder=6)
        axes[0,0].set_ylabel('γ₁ complex / γ₁ phase');save(fig)
        # Normalized complex correlations themselves, all modes and conditions.
        curves=csvread(DATA/'exemplar_curves.csv')
        for m in MODELS:
            inputs=sorted({num(r,'f') for r in curves if r['model']==m});nrows=max(3,int(np.ceil(len(inputs)/4)))
            fig,axes=plt.subplots(nrows,4,figsize=(11.7,8.3 if nrows==3 else 10.8));fig.subplots_adjust(left=.065,right=.98,bottom=.09,top=.9,hspace=.65,wspace=.42)
            fig.text(.055,.97,'S4.11  Normalized complex autocorrelations: '+NAMES[m],fontsize=12,va='top')
            for ax,f in zip(axes.flat,inputs):
                for size in sizes(m):
                    # First, second, fourth and eighth modes display the same
                    # measure across scales; complete 0–32 curves are in CSV.
                    for mode,ls in ((1,'-'),(2,'--'),(4,':'),(8,'-.')):
                        rr=sorted([r for r in curves if r['model']==m and num(r,'f')==f and num(r,'C_label')==size and num(r,'n')==mode],key=lambda r:num(r,'lag'))
                        ax.plot([num(r,'lag') for r in rr],[num(r,'correlation') for r in rr],ls=ls,color=color(m,size),lw=.55)
                ax.axhline(.5,color='.7',lw=.4);ax.set_xscale('log');ax.set_ylim(-.25,1.1);ax.set_title(f'f = {f:g}',fontsize=7);style(ax)
            for ax in list(axes.flat)[len(inputs):]:ax.axis('off')
            fig.supxlabel('Lag (native time)',y=.052,fontsize=8);fig.supylabel('Normalized complex correlation',x=.015,fontsize=8)
            fig.text(.06,.02,'n = 1 solid, 2 dashed, 4 dotted, 8 dash-dot; size light→dark. All 33 modes including n = 0 are retained in compact CSV.',fontsize=7)
            save(fig)
    temporary.replace(OUT/'figure4_supplement.pdf')
