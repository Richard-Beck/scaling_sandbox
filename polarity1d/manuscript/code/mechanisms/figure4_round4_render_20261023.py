"""Revised Figure 4E: normalized drift and two explicit decompositions."""
import math
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import figure4_round3_render_20261022 as R


def decomposition(fig,grid,data,exported,omissions):
    rows=[r for r in R.primary(R.read_csv(data/'drift_decomposition.csv'))
          if R.n(r,'flicker_frames')==5 and R.n(r,'MSD_range_fraction')==.25]
    detector=[r for r in R.read_csv(data/'detector_free_decomposition.csv') if R.n(r,'detector_free_step')==1]
    source=R.primary(R.read_csv(data/'analysis_lag_metrics.csv'))
    R.write_csv(data/'multiple_front_fraction.csv',source)
    keys=['Ds_over_smallest','lambda_eff','kpred_over_kobs','multiple_front_fraction']
    labels=[r'$D_s/D_s(C_{\min})$',r'$\lambda_{\mathrm{eff}}$',r'$k_{\mathrm{pred}}/k_{\mathrm{obs}}$','Multiple-front\nfraction']
    plotted=[]
    for i,(quantity,label) in enumerate(zip(keys,labels)):
        for j,model in enumerate(R.MODELS):
            ax=fig.add_subplot(grid[i,j]);rr=[r for r in rows if r['model']==model]
            if i==0:
                ax.set_yscale('log',base=2 if model!='goryachev' else 10)
                ax.set_ylim((.25,8) if model!='goryachev' else (1e-5,2))
                if model!='goryachev':ax.set_yticks([.25,.5,1,2,4,8],['¼','½','1','2','4','8'])
                else:ax.set_yticks([1e-4,1e-2,1],[r'$10^{-4}$',r'$10^{-2}$','1'])
                ax.yaxis.set_minor_formatter(NullFormatter());ax.axhline(1,color='.6',lw=.5)
            elif i==1:
                values=[R.n(r,quantity) for r in rr+detector if r['model']==model and R.finite(R.n(r,quantity))]
                top=max(values)*1.8 if values and max(values)>0 else 1.
                ax.set_yscale('symlog',linthresh=1e-7,linscale=.6);ax.set_ylim(0,top)
                ax.set_yticks([v for v in [0,1e-6,1e-4,1e-2,1] if v<top])
            elif i==2:
                ax.set_yscale('log');ax.set_ylim(.003,1e4)
                ax.set_yticks([.01,1,100,1e4],[r'$10^{-2}$','1',r'$10^2$',r'$10^4$'])
                ax.axhline(1,color='.6',lw=.5)
            else:ax.set_ylim(0,1);ax.set_yticks([0,.5,1])
            methods=[('scorable_events',rr,'-')]
            if i in (1,2):methods.append(('detector_free',[r for r in detector if r['model']==model],'--'))
            for method,records,ls in methods:
                for size in R.SIZES:
                    relevant=sorted([r for r in records if math.isclose(R.equivalent_size(r),size)],key=lambda r:R.n(r,'f'))
                    x=[];y=[];color=R.size_color(model,size)
                    for row in relevant:
                        value=R.n(row,quantity);x.append(R.n(row,'f'))
                        defined=R.finite(value) and (value>=0 if i in (1,3) else value>0)
                        y.append(value if defined else math.nan)
                        if not defined:
                            reason=('Observed decay rate unavailable' if i==2 and not R.finite(R.n(row,'kobs')) else
                                    'Drift or smallest-size reference unavailable' if i==0 else row.get('jump_status',row.get('status','undefined')))
                            omissions.append(dict(panel='E',model=model,C_label=row['C_label'],f=row['f'],quantity=quantity,method=method,reason=reason))
                            continue
                        lower=R.n(row,quantity+'_ci_lower');upper=R.n(row,quantity+'_ci_upper')
                        forced=R.yes(row.get('forced_combined'))
                        plotted.append(dict(panel='E',model=model,C_label=row['C_label'],f=row['f'],quantity=quantity,method=method,
                                            value=value,ci_lower=lower,ci_upper=upper,forced_combined=forced,
                                            n_switches=row.get('n_switches'),n_scored_switches=row.get('n_scored_switches'),
                                            jump_scored_fraction=row.get('jump_scored_fraction'),analysis_lag=row.get('analysis_lag'),
                                            n_seeds=row.get('n_seeds'),source='drift_decomposition.csv' if method=='scorable_events' else 'detector_free_decomposition.csv'))
                        if method=='scorable_events':R.bounded_point(ax,R.n(row,'f'),value,lower,upper,color,open_marker=forced,ms=1.6,lw=.35)
                        elif i==2 and not .003<=value<=1e4:
                            ax.plot(R.n(row,'f'),.003 if value<.003 else 1e4,'v' if value<.003 else '^',ms=2.5,mfc='white',mec=color,mew=.5,clip_on=False)
                    ax.plot(x,y,ls,color=color,lw=.7,zorder=2 if method=='detector_free' else 1)
            R.log_f(ax)
            if i<3:ax.tick_params(labelbottom=False);ax.set_xlabel('')
            if j==0:ax.set_ylabel(label,labelpad=2)
            if i==0:R.facet(ax,R.NAMES[model])
            if i==1 and j==0:
                ax.legend([Line2D([],[],color='black',lw=.7),Line2D([],[],color='black',lw=.7,ls='--')],
                          ['Scored events','Mixture fit'],frameon=False,fontsize=5,loc='center left',bbox_to_anchor=(.02,.4),handlelength=1.3,borderpad=.1)
            ax.tick_params(axis='y',labelsize=5.3);R.style_axis(ax)
    exported['E']=plotted


def supplement_pages(data,pdf,finish,grid,trend):
    """Expose scoring coverage, missing-angle bounds, and full closure failures."""
    coverage=R.read_csv(data/'jump_scoring_coverage.csv')
    if not coverage:return
    trend(coverage,'jump_scored_fraction','Fraction of counted replacements with scorable jump angles',limits=(-.035,1.035))
    import matplotlib.pyplot as plt
    fig,axes=grid()
    for ax,model in zip(axes.flat,R.MODELS):
        for size in R.SIZES:
            rr=sorted([r for r in coverage if r['model']==model and math.isclose(R.equivalent_size(r),size)],key=lambda r:R.n(r,'f'))
            x=np.array([R.n(r,'f') for r in rr]);lo=np.array([R.n(r,'lambda_eff_missing_angle_lower') for r in rr]);hi=np.array([R.n(r,'lambda_eff_missing_angle_upper') for r in rr])
            color=R.size_color(model,size);ax.fill_between(x,lo,hi,color=color,alpha=.13)
            ax.plot(x,[R.n(r,'lambda_eff') for r in rr],color=color,lw=.7)
        ax.set_yscale('symlog',linthresh=1e-7);ax.set_ylim(bottom=0);R.log_f(ax);R.facet(ax,R.NAMES[model]);R.style_axis(ax)
        ax.set_ylabel(r'$\lambda_{\mathrm{eff}}$')
    fig.text(.1,.025,'Shading: every unscored event assigned loss 0 or 2. Lines: all-event rate × mean loss among scorable events. These identification bounds are distinct from bootstrap intervals.',fontsize=6,wrap=True)
    finish(fig,'Missing jump angles: assumption-free loss bounds and conditional estimate')
    rows=[r for r in R.read_csv(data/'drift_decomposition.csv') if R.n(r,'flicker_frames')==5 and R.n(r,'MSD_range_fraction')==.25]
    detector=R.read_csv(data/'detector_free_decomposition.csv')
    for method,rr in [('Scorable events',rows),('Mixture at shortest cadence',[r for r in detector if R.n(r,'detector_free_step')==1]),
                      ('Mixture at twice cadence',[r for r in detector if R.n(r,'detector_free_step')==2])]:
        trend(rr,'kpred_over_kobs',method+': full predicted / observed rate range',logy=True)
