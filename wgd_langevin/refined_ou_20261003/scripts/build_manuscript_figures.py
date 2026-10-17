"""Regenerate the established integrated OU figure with refined-data fits."""
from pathlib import Path
import sys, json, csv, hashlib
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator, NullFormatter
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
CAMPAIGN=ROOT/"wgd_langevin/refined_ou_20261003"
HERE = Path(__file__).resolve().parent
OUT = CAMPAIGN / 'figures'
REVIEW = OUT
sys.path.insert(0, str(HERE))
import figure_panels as panels
import figure_tracks as tracks

MINUS, PLUS, NAVY, GRAY = tracks.MINUS, tracks.PLUS, tracks.NAVY, tracks.GRAY
FIT_DIR = CAMPAIGN / 'results'
SURFACE = FIT_DIR / 'main_relative_surface.csv'

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,
                     'axes.spines.top':False,'axes.spines.right':False})

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def add_letter(ax, letter):
    ax.text(-.055, 1.04, letter, transform=ax.transAxes, color=NAVY,
            fontsize=14, fontweight='bold', va='bottom')

def build_figure3():
    REVIEW.mkdir(parents=True, exist_ok=True)
    scenarios = panels.load_fit_scenarios()
    primary = scenarios['primary']
    observation_csv, observation_meta = panels.write_observation_table()
    rows, row_prov = panels.build_observation_rows(primary_fit=primary)
    if len(rows) != 9 or not {'Jemaa2017','Jemaa2023'}.issubset({r['study'] for r in rows}):
        raise ValueError('Observed/predicted panel does not contain the nine required rows')

    fit = primary['record']
    z = np.asarray(primary['z'], dtype=float)
    pars = np.exp(z[:20].reshape(5,2,2))
    names = ['Xuan2018','Mallin2023','Jemaa2017','Jemaa2023','Hosny_baseline']
    scenarios_order = ['primary','prior5','prior20','strong','weak']
    group_label = ['WGD−','WGD+']

    # Selected model-generated tracks and the exact fitted OU MSD expectation.
    track_archive = np.load(FIT_DIR/'tracks_primary.npz')
    chosen, track_metrics, selected_indices = [], [], {}
    for group, color in ((0,MINUS),(1,PLUS)):
        path, metric, selection = tracks.choose_representative(track_archive[f'Xuan2018_{group}'])
        chosen.append(path); track_metrics.append(metric)
        selected_indices[group_label[group]] = selection['index']
    tau_eta = pars[0]
    lag = max(float(r['lag_h']) for r in fit['rows'] if r['metric']=='msd_curve')
    expected_msd = [float(tracks.msd_expected(lag, *tau_eta[g])) for g in (0,1)]
    fig = plt.figure(figsize=(16,12))
    gs = fig.add_gridspec(3, 2, left=.075,right=.965,bottom=.12,top=.92,
                          height_ratios=[1.0,1.15,.78], width_ratios=[1.08,1.0],
                          hspace=.38,wspace=.24)
    # A: all observed/predicted folds in one panel.
    axA = fig.add_subplot(gs[0,0])
    panels.draw_observed_predicted(axA, rows)
    add_letter(axA,'A')
    axA.get_legend().remove()
    axA.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor='#5b7380',markeredgecolor='white',label='Observed fold',markersize=6),
                        Line2D([],[],marker='D',color='none',markerfacecolor='#d36c32',markeredgecolor='white',label='Predicted fold',markersize=5.5),
                        Line2D([],[],color='#33454e',linestyle='--',label='Equal groups')],
               frameon=False,ncol=3,fontsize=7,loc='lower left',bbox_to_anchor=(0,1.02))
    axA.xaxis.set_major_locator(FixedLocator([.5,1,2,3]))
    axA.xaxis.set_major_formatter(FixedFormatter(['0.5','1','2','3']))
    axA.xaxis.set_minor_locator(NullLocator());axA.xaxis.set_minor_formatter(NullFormatter())
    axA.tick_params(labelsize=8)

    # B: logarithmic axes expose the crossover without compressing short lags.
    axmsd = fig.add_subplot(gs[0,1])
    add_letter(axmsd,'B')
    curves = [sorted([r for r in fit['rows'] if r['study']=='Xuan2018' and r['metric']=='msd_curve' and int(r['group'])==g],key=lambda r:r['lag_h']) for g in (0,1)]
    lower=max(curve[0]['lag_h'] for curve in curves)
    upper=min(curve[-1]['lag_h'] for curve in curves)
    tt = np.geomspace(lower,upper,500)
    observed = []
    for g,(label,color) in enumerate(((group_label[0],MINUS),(group_label[1],PLUS))):
        curve=curves[g]
        lag_values=np.array([r['lag_h'] for r in curve])
        values=np.array([r['obs'] for r in curve])
        observed.append((lag_values,values))
        axmsd.plot(lag_values,values,color=color,lw=1.7,label=f'{label} population curve')
        axmsd.plot(tt,tracks.msd_expected(tt,*tau_eta[g]),color=color,lw=1.8,ls='--',label=f'{label} OU fit')
    native=np.unique(np.concatenate([a[(a>=lower)&(a<=upper)] for a,_ in observed]))
    difference=lambda t:np.interp(t,*observed[0])-np.interp(t,*observed[1])
    difference_values=difference(native)
    crossings=[float(brentq(difference,native[i],native[i+1])) for i in np.flatnonzero(difference_values[:-1]*difference_values[1:]<0)]
    fit_difference=lambda t:float(tracks.msd_expected(t,*tau_eta[0])-tracks.msd_expected(t,*tau_eta[1]))
    fit_values=np.array([fit_difference(t) for t in tt])
    fitted_crossings=[float(brentq(fit_difference,tt[i],tt[i+1])) for i in np.flatnonzero(fit_values[:-1]*fit_values[1:]<0)]
    for crossing in crossings:
        value=float(np.interp(crossing,*observed[0]))
        axmsd.scatter(crossing,value,s=38,facecolors='white',edgecolors=NAVY,zorder=5)
        axmsd.annotate(f'Population crossover\n{crossing:.2f} h',xy=(crossing,value),xytext=(.50,.30),textcoords='axes fraction',fontsize=8,color=NAVY,
                       arrowprops=dict(arrowstyle='->',color=NAVY,lw=.9),ha='center')
    axmsd.set(xscale='log',yscale='log',xlim=(lower*.92,upper*1.07),xlabel='Lag time (h; log scale)',ylabel='MSD (µm²; log scale)')
    axmsd.set_xticks([.2,.5,1,2,5,10],['0.2','0.5','1','2','5','10'])
    axmsd.xaxis.set_minor_formatter(NullFormatter())
    axmsd.grid(alpha=.15)
    axmsd.legend(frameon=False,ncol=2,fontsize=7,loc='upper left')
    axmsd.tick_params(labelsize=8)

    # C: enlarged chosen WGD−/WGD+ simulated paths, with no MSD subplot.
    sub = gs[1,0].subgridspec(1,2,wspace=.24)
    lim = float(max(np.abs(t).max() for t in chosen)*1.13)
    for g,(label,color) in enumerate(((group_label[0],MINUS),(group_label[1],PLUS))):
        ax = fig.add_subplot(sub[0,g])
        if g == 0: add_letter(ax,'C')
        tracks.axes_track(ax,chosen[g],color,linewidth=1.5)
        tracks.style_axis(ax,lim)
        tracks.arrow(ax,chosen[g][0],chosen[g][-1],NAVY,1.0)
        for j in (18,36,54): tracks.arrow(ax,chosen[g][j],chosen[g][j+1],color,.65)
        for a,b in ((0,12),(60,72)):
            ax.plot(chosen[g][[a,b],0],chosen[g][[a,b],1],ls='--',color=GRAY,lw=.85)
        ax.set_title(f"{label} · {track_metrics[g]['mean_sampled_speed_um_h']:.1f} µm/h\n"
                     f"{track_metrics[g]['path_length_um']:.0f} µm path · directness {track_metrics[g]['directness']:.2f}",
                     loc='left',color=color,fontweight='bold',fontsize=9,pad=6)
        ax.set_xlabel('x (µm)',fontsize=8); ax.set_ylabel('y (µm)',fontsize=8)
        ax.tick_params(labelsize=7)
    # D: paired data and regularized loss surfaces, with a shared raw-loss scale.
    surface_grid = gs[1,1].subgridspec(1,3,width_ratios=[1,1,.055],wspace=.38)
    axC1,axC2 = fig.add_subplot(surface_grid[0,0]),fig.add_subplot(surface_grid[0,1])
    cax=fig.add_subplot(surface_grid[0,2])
    with SURFACE.open() as surface_handle:
        surface_rows=list(csv.DictReader(surface_handle))
    scores=[float(r[k]) for r in surface_rows for k in ('fixed_drift_delta_data','fixed_drift_delta_total')]
    norm = matplotlib.colors.AsinhNorm(linear_width=1., vmin=min(-.01,min(scores)), vmax=max(scores))
    m1 = panels.draw_relative_surface(axC1,csv_path=SURFACE,score_column='fixed_drift_delta_data',norm=norm,
                                      title='Data misfit',cmap='magma')
    panels.draw_relative_surface(axC2,csv_path=SURFACE,score_column='fixed_drift_delta_total',norm=norm,
                                 title='Total objective',cmap='magma')
    axC2.set_ylabel('')
    for ax in (axC1,axC2):
        ax.set_xlabel('WGD+ / WGD− persistence fold',fontsize=8)
        ax.set_ylabel('WGD+ / WGD− noise fold',fontsize=8)
        ax.tick_params(labelsize=7)
    cb=fig.colorbar(m1,cax=cax)
    cb.set_label('Change in score from fitted point',fontsize=8)
    cb.set_ticks([0,1,10,100,1000,10000])
    cb.ax.tick_params(labelsize=6)
    add_letter(axC1,'D')
    axC1.text(0,1.13,'WGD− fixed; WGD+ persistence and noise shifted',transform=axC1.transAxes,
              ha='left',va='bottom',color=NAVY,fontsize=9,fontweight='bold')
    axC2.legend(handles=[Line2D([],[],marker='o',color='none',markerfacecolor='white',markeredgecolor='#242f35',label='Equal group means',markersize=5),
                         Line2D([],[],marker='*',color='none',markerfacecolor='#111111',markeredgecolor='white',label='Selected fit',markersize=8)],
                frameon=False,fontsize=6.5,loc='upper right')

    # E: primary fit paired persistence and velocity-noise parameters.
    dgs=gs[2,0].subgridspec(1,2,wspace=.35)
    d_tau_ax=fig.add_subplot(dgs[0,0]); d_sigma_ax=fig.add_subplot(dgs[0,1])
    for k,(ax,label,unit) in enumerate(((d_tau_ax,'Persistence','τ (min)'),
                                         (d_sigma_ax,'Velocity noise','σ (µm h⁻³ᐟ²)'))):
        offsets=[-.12,.12]
        for j,name in enumerate(names):
            ax.plot([pars[j,0,k]* (60 if k==0 else 1),pars[j,1,k]*(60 if k==0 else 1)],
                    [j+offsets[0],j+offsets[1]],color='#b8c2c9',lw=1,zorder=1)
            for g,(c,off) in enumerate(((MINUS,offsets[0]),(PLUS,offsets[1]))):
                ax.scatter(pars[j,g,k]*(60 if k==0 else 1),j+off,s=25,color=c,
                           edgecolor='white',linewidth=.45,zorder=3)
        ax.set_yticks(np.arange(5),names);ax.invert_yaxis()
        ax.set_xscale('log')
        ax.set_xlabel(unit,fontsize=8);ax.set_title(label,loc='left',color=NAVY,fontsize=9)
        ax.grid(axis='x',alpha=.15);ax.tick_params(labelsize=7)
    dlegend=[Line2D([],[],marker='o',color='none',markerfacecolor=c,label=l,markersize=5)
             for l,c in ((group_label[0],MINUS),(group_label[1],PLUS))]
    d_tau_ax.legend(handles=dlegend,frameon=False,fontsize=7,loc='lower right')
    add_letter(d_tau_ax,'E')

    # F: original prior/pooling matrix of actual WGD+ / WGD− persistence folds.
    axE=fig.add_subplot(gs[2,1])
    scenario_folds=[]
    for scenario in scenarios_order:
        sp=np.exp(np.asarray(scenarios[scenario]['z'][:20]).reshape(5,2,2))
        scenario_folds.append(sp[:,:,0][:,1]/sp[:,:,0][:,0])
    folds=np.asarray(scenario_folds).T
    pooled=np.exp(np.mean(np.log(folds),axis=0))
    display=np.vstack([folds,pooled])
    logfold=np.log2(display)
    im=axE.imshow(logfold,aspect='auto',cmap='YlOrRd',vmin=0.,vmax=max(1.,float(np.max(logfold))))
    axE.set_xticks(range(5),['Primary','5-min prior','20-min prior','Strong pooling','Weak pooling'],rotation=22,ha='right')
    axE.set_yticks(range(6),names+['Pooled geometric mean'])
    axE.set_title('WGD+ / WGD− persistence fold',loc='left',color=NAVY,fontsize=9)
    for i in range(6):
        for j in range(5): axE.text(j,i,f'{display[i,j]:.1f}',ha='center',va='center',fontsize=6.5,
                                     color='white' if logfold[i,j]>4.6 else '#40230f')
    cb2=fig.colorbar(im,ax=axE,fraction=.046,pad=.025,ticks=None);cb2.set_label('log₂ persistence fold',fontsize=7);cb2.ax.tick_params(labelsize=7)
    add_letter(axE,'F')

    fig.suptitle('Figure 3 | Migration summaries, simulated tracks and fitted OU dynamics',
                 fontsize=15,color=NAVY,fontweight='bold')
    caption = ("A, observed and primary-fit OU-predicted WGD+ / WGD− folds for nine migration summaries (the last shared population lag summarizes the Xuan MSD fold, with native observed/predicted traces interpolated to that shared lag); "
               "Hosny uses group medians. B, full digitized population MSD curves (solid) and fitted OU ensemble MSD (dashed), on logarithmic lag and MSD axes to expose the crossover. "
               "The open circle marks the population-curve crossover interpolated between native trace samples, approximately "+f"{crossings[0]:.2f}"+" hours; this is a descriptive digitized crossing, not a precise population estimate. "
               "C, enlarged representative primary-fit model-generated Xuan2018 tracks. D, with WGD− fixed, relative WGD+ persistence/noise data-misfit and regularized "
               "objective surfaces. E, paired primary-fit persistence and velocity-noise estimates. "
               "F, WGD+ / WGD− persistence folds across prior/pooling scenarios, including their pooled geometric mean. Track paths are simulations, "
               "not observed individual cells. Surface scores use the model's heuristic weighted objective; "
               "their scale is not a calibrated likelihood or biological unit. Prior/pooling contrasts are "
               "conditional sensitivity analyses, not uncertainty intervals.")
    fig.text(.075,.035,'A observed vs primary-fit folds · B population MSD + OU curves · C simulated tracks · D WGD−-fixed relative surfaces · E primary parameters · F prior/pooling sensitivity',
             ha='left',va='bottom',fontsize=7.6,color='#394b56')
    fig.text(.075,.015,'C speed = path / 12 h; directness = net / path. Dashed 2 h chords contribute squared displacements to track MSD. Full population MSD curves are fitting targets; simulated tracks are illustrations.',
             ha='left',va='bottom',fontsize=7.2,color='#394b56')
    figure_path=REVIEW/'main_03_migration_OU.png'
    for suffix in ('png','pdf','svg'):
        fig.savefig(figure_path.with_suffix('.'+suffix),dpi=300,bbox_inches='tight')
    plt.close(fig)

    canonical=figure_path
    # Persistence export retained under its established name, now with all five scenarios.
    with (OUT/'figure3_persistence.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['study','reference_minutes','polyploid_minutes',*scenarios_order])
        for j,name in enumerate(names):
            w.writerow([name,pars[j,0,0]*60,pars[j,1,0]*60,*[np.exp(np.asarray(scenarios[s]['z'][:20]).reshape(5,2,2))[j,1,0]/np.exp(np.asarray(scenarios[s]['z'][:20]).reshape(5,2,2))[j,0,0] for s in scenarios_order]])
        scenario_fold_arrays=[np.exp(np.asarray(scenarios[s]['z'][:20]).reshape(5,2,2))[:,:,0][:,1]/
                              np.exp(np.asarray(scenarios[s]['z'][:20]).reshape(5,2,2))[:,:,0][:,0]
                              for s in scenarios_order]
        pooled_folds=[float(np.exp(np.mean(np.log(fold)))) for fold in scenario_fold_arrays]
        pooled_reference=float(np.exp(np.mean(np.log(pars[:,0,0]*60))))
        pooled_polyploid=float(np.exp(np.mean(np.log(pars[:,1,0]*60))))
        w.writerow(['Pooled geometric mean',pooled_reference,pooled_polyploid,*pooled_folds])
    with (OUT/'figure3_primary_noise.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['study','reference_sigma_um_h_1p5','polyploid_sigma_um_h_1p5','fold'])
        for j,name in enumerate(names): w.writerow([name,pars[j,0,1],pars[j,1,1],pars[j,1,1]/pars[j,0,1]])

    files=[ROOT/'wgd_langevin/persistence_only/campaign.py', HERE/'main_surface.py', SURFACE.with_suffix('.json'), CAMPAIGN/'data/adapter.py', CAMPAIGN/'data/observations.csv', CAMPAIGN/'data/xuan_msd_curve.csv', CAMPAIGN/'data/hosny_observations.csv', primary['path'],HERE/'fit.py',HERE/'model.py',FIT_DIR/'tracks_primary.npz',SURFACE,
           observation_csv,observation_meta,Path(__file__),HERE/'figure_panels.py',HERE/'figure_tracks.py']
    provenance={
        'created_utc':datetime.now(timezone.utc).isoformat(),
        'builder':str(Path(__file__).relative_to(ROOT)),
        'figure_png':str(figure_path.relative_to(ROOT)),
        'canonical_png':str(canonical.relative_to(ROOT)),
        'caption':caption,
        'selected_primary_fit':primary['path'].name,
        'selected_fit_source_hashes_recorded':True,
        'source_sha256':{str(p.relative_to(ROOT)):digest(p) for p in files},
        'observation_rows':9,
        'track_selection_indices':selected_indices,
        'fitted_msd_at_final_curve_lag':dict(zip(group_label,expected_msd)),
        'final_curve_lag_h':lag,
        'full_population_curves_plotted':True,
        'panel_layout':'A folds / B MSD; C tracks / D surfaces; E parameters / F sensitivity',
        'msd_axes':'logarithmic lag and logarithmic MSD',
        'population_curve_crossovers_h':crossings,
        'fitted_ou_crossovers_h':fitted_crossings,
        'surface_source':str(SURFACE.relative_to(ROOT)),
        'surface_score_scale':'shared asinh color normalization; underlying score units are the heuristic weighted objective',
        'scenario_order':scenarios_order,
        'scenario_fit_files':{k:{'path':str(v['path'].relative_to(ROOT)),'sha256':digest(v['path'])} for k,v in scenarios.items()},
        'observation_provenance':row_prov,
        'notes':'PNG, PDF and SVG. Track paths are model simulations. Surface scores are not likelihoods; prior/pooling results are conditional sensitivity analyses.'}
    (REVIEW/'main_03_migration_OU_provenance.json').write_text(json.dumps(provenance,indent=2))
    (REVIEW/'main_03_migration_OU_caption.md').write_text('# Figure 3 caption\n\n'+caption+'\n')
    return provenance

if __name__=='__main__':
    print(json.dumps(build_figure3(),indent=2))
