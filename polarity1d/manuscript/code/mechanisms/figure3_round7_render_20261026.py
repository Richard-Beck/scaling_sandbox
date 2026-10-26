"""Figure3 A–F and compact scientific provenance for round7; no native runs."""
from pathlib import Path
import os,sys,json,shutil,hashlib
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
import figure3_conduction_20261026 as campaign
ROOT=campaign.ROOT;MS=campaign.MS;OUT=campaign.OUT;DATA=campaign.DATA;META=campaign.META
OLD=MS/'figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure3'
MODELS=['wave_pinning','holmes_model3','spring_ring','goryachev']
NAMES=dict(wave_pinning='Wave-pinning',holmes_model3='Holmes',spring_ring='Spring',goryachev='Goryachev')
COLORS=dict(wave_pinning='#236894',holmes_model3='#387a50',spring_ring='#9a741d',goryachev='#af3e45')
SOURCES={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    SOURCES[str(Path(p).relative_to(ROOT))]=sha(p)
    return json.loads(Path(p).read_text())
def csv(p):
    SOURCES[str(Path(p).relative_to(ROOT))]=sha(p)
    return pd.read_csv(p)
def rel(p):return str(Path(p).relative_to(ROOT))
def write(p,v):campaign.write(p,v)
def shade(model,i,n):return np.array(to_rgb(COLORS[model]))*(.45+.55*i/max(n-1,1))+(1-(.45+.55*i/max(n-1,1)))
def prepare():
    for folder in ['data','metadata','vectors','viewers']:(OUT/'supporting'/folder).mkdir(parents=True,exist_ok=True)
    for name in ['half_heading_readouts.csv','regime_lines.csv','undefined_heading_bands.csv','susceptibility_and_contributions.csv','map_amplitudes.csv','deterministic_experiments.csv','linear_range_limits.csv','matched_B_susceptibility.csv']:
        csv(OLD/'supporting/data'/name);shutil.copy2(OLD/'supporting/data'/name,DATA/name)
    fronts=read(META/'isolated_front.json');frozen_manifest=read(campaign.edge.RUN/'manifest.json');
    for name,digest in frozen_manifest['source_sha256'].items():
        path=campaign.edge.RUN/name;assert sha(path)==digest;SOURCES[rel(path)]=digest
    tables=[]
    for model in MODELS[:3]:
        assert read(META/f'{model}_progress.json')['complete'],model
        tables.append(csv(DATA/f'conduction_{model}.csv'))
    gain=pd.concat(tables,ignore_index=True);gain.to_csv(DATA/'conduction_sweep.csv',index=False)
    # Native theta from independent stationary rings; lambda/kappa are isolated-front values.
    wp=gain[(gain.model=='wave_pinning')&(gain.D==10)&gain.defined].copy();theta=float(wp.theta.median())
    cwp=20/(fronts['kappa']*fronts['delta_u']*theta*(1-theta))
    fits=[];cross={'wave_pinning':cwp};slope_rows=[]
    for model in ['holmes_model3','spring_ring']:
        q=gain[(gain.model==model)&gain.defined].copy();Dmax=q.D.max();ref=q[q.D==Dmax].set_index('size').gain
        q['free_gain']=q['size'].map(ref);q['penalty']=q.free_gain/q.gain-1
        use=q[(q.D>=1)&(q.D<Dmax)].copy()
        x=use.C/use.D;beta=float(np.sum(x*use.penalty)/np.sum(x*x));Cstar=campaign.NATIVE[model]/beta
        native=q[q.D==campaign.NATIVE[model]].sort_values('C'); nativecoef=np.polyfit(native.C,1/native.gain,1); nativeCstar=float(nativecoef[1]/nativecoef[0]); cross[model]=Cstar if model=='spring_ring' else nativeCstar
        fits.append(dict(model=model,C_star=cross[model],apparent_native_C_star=nativeCstar,penalty_C_star=Cstar,native_inverse_gain_intercept=float(nativecoef[1]),native_inverse_gain_slope=float(nativecoef[0]),beta=beta,D_native=campaign.NATIVE[model],fit=('Spring plotted transport C*: through-origin G(Dmax,C)/G(D,C)-1 against C/D over D>=1; apparent native inverse-gain crossover retained separately' if model=='spring_ring' else 'Holmes plotted C*: native inverse-gain intercept/slope over five sizes; diffusivity penalty fit retained separately') ,D_reference=Dmax,n=len(use),rms_penalty_error=float(np.sqrt(np.mean((use.penalty-beta*x)**2)))))
        for D,r in q.groupby('D'):
            coef=np.polyfit(r.C,1/r.gain,1);slope_rows.append(dict(model=model,D=D,intercept=coef[1],slope=coef[0],C_star=coef[1]/coef[0],Cstar_over_D=coef[1]/coef[0]/D))
    fits.insert(0,dict(model='wave_pinning',C_star=cwp,theta=theta,kappa=fronts['kappa'],delta_u=fronts['delta_u'],lambda_push=fronts['lambda_push'],D_native=10,fit='isolated-front prediction 2 D/(kappa delta_u theta(1-theta)); theta from stationary native ring; no ring gain fit'))
    write(META/'conduction_law_parameters.json',dict(models=fits));pd.DataFrame(slope_rows).to_csv(DATA/'inverse_gain_fits.csv',index=False)
    penalties=[]
    for model in MODELS[:3]:
        q=gain[(gain.model==model)&gain.defined];ref=q[q.D==q.D.max()].set_index('size').gain
        for r in q.itertuples():penalties.append(dict(model=model,size=r.size,C=r.C,D=r.D,C_over_D=r.C/r.D,penalty=float(ref.loc[r.size]/r.gain-1),D_reference=q.D.max(),law_C_star=cross[model]))
    penalty=pd.DataFrame(penalties);penalty.to_csv(DATA/'conduction_penalty.csv',index=False)
    # Deterministic kymographs: exact retained cue-on + withdrawal histories.
    kymos=[]
    for label in ['i','ii']:
        folder=MS/'data/runs/figure3_round3_relaxation_v2_20261022/workers'/label;status=read(folder/'status.json');p=status['job']['placement'];C=p['size_parameter'];n=p['cells'];origin=p['initial_10_percent_contour']['origin'];pieces=[];times=[]
        for part in ['cue_on','cue_off']:
            path=folder/f'{part}.npz';SOURCES[rel(path)]=sha(path)
            with np.load(path) as z:
                tt=z['time'].astype(float)+(status['cue_on_duration'] if part=='cue_off' else 0);uu=z['states'][:,:n].astype(float)
            if part=='cue_off':tt,uu=tt[1:],uu[1:]
            times.extend(tt);pieces.extend(uu)
        x=(np.arange(n)+.5)*C/n;x=(x-origin+C/2)%C-C/2;order=np.argsort(x);t=np.asarray(times);u=np.asarray(pieces)[:,order];x=x[order]
        # Compact sampled coordinates support exact displayed color mesh.
        ti=np.unique(np.linspace(0,len(t)-1,min(350,len(t)),dtype=int));xi=np.arange(0,n,2);t=t[ti];u=u[ti][:,xi];x=x[xi]
        np.savez_compressed(DATA/f'A_kymograph_{label}.npz',time=t,x=x,activity=u,withdrawal=status['cue_on_duration'],cue_left=p['cue_near_edge_relative'],cue_right=p['cue_far_edge_relative'])
        kymos.append((label,t,x,u,status,p))
    return gain,penalty,cross,fronts,kymos

def render(gain,penalty,cross,fronts,kymos):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','text.color':'black','axes.titlecolor':'black'})
    fig=plt.figure(figsize=(18,16));positions={}
    # A: one activity field in two placement facets.
    ga=fig.add_gridspec(1,2,left=.065,right=.45,bottom=.73,top=.95,wspace=.25)
    for j,(label,t,x,u,status,p) in enumerate(kymos):
        ax=fig.add_subplot(ga[0,j]);mesh=ax.pcolormesh(x,t,u,cmap='magma',shading='auto',rasterized=True,vmin=.17,vmax=1.15);ax.axhline(status['cue_on_duration'],ls='--',color='white',lw=1);ax.plot([p['cue_near_edge_relative'],p['cue_far_edge_relative']],[0,0],color='cyan',lw=4,clip_on=False);ax.set(xlabel='Arc position',ylabel='Time' if j==0 else '',title='Overlapping cue' if j==0 else 'Separated cue')
    cb=fig.colorbar(mesh,ax=fig.axes[:2],fraction=.035,pad=.02);cb.set_label('Active concentration u')
    gb=fig.add_gridspec(1,4,left=.52,right=.97,bottom=.73,top=.95,wspace=.55)
    f=csv(DATA/'half_heading_readouts.csv');bands=csv(DATA/'undefined_heading_bands.csv')
    for j,model in enumerate(MODELS):
        key='spring' if model=='spring_ring' else model;q=f[(f.model==key)&(~f.no_turn)];ax=fig.add_subplot(gb[0,j]);ax.set_title(NAMES[model]);sizes=sorted(q.size_parameter.unique())
        for b in bands[bands.model==key].itertuples():ax.axvspan(b.left,b.right,color='.92',zorder=-2)
        for i,size in enumerate(sizes):
            z=q[q.size_parameter==size].sort_values('gap');col=shade(model,i,len(sizes));valid=z[z.t50.notna()];ax.plot(valid.gap,valid.t50,'o-',ms=2.5,color=col,label=f'{size:g}')
            for r in z[z.status=='right_censored_saved_heading_at_completed_horizon'].itertuples():
                ax.plot(r.gap,r.lower,'^',mfc='white',mec=col,ms=5);ax.annotate('',(r.gap,r.lower*1.6),(r.gap,r.lower),arrowprops=dict(arrowstyle='->',color=col,lw=.8))
        ax.set(xlabel='Cue edge gap',yscale='log');ax.legend(title='Size',fontsize=6.5,title_fontsize=7,frameon=False,loc='upper left');ax.set_ylabel('Half-turn time' if j==0 else '')
    gc=fig.add_gridspec(1,4,left=.065,right=.97,bottom=.48,top=.64,wspace=.5);maps=csv(DATA/'susceptibility_and_contributions.csv')
    for j,model in enumerate(MODELS):
        size=2 if model=='spring_ring' else 40;q=maps[(maps.model==model)&(maps['size']==size)].sort_values('x');ax=fig.add_subplot(gc[0,j]);ax.set_title(NAMES[model]);col=COLORS[model];ax.plot(q.x,q.standard_pulse_gain,color=col,lw=1.7,label='Gaussian-pulse susceptibility');ax.axhline(0,color='.7',lw=.6);ax.set(xlabel='Arc position',ylabel='Displacement / pulse dose' if j==0 else '')
        twin=ax.twinx();twin.fill_between(q.x,q.activity,color='.6',alpha=.18);twin.plot(q.x,q.activity,color='.55',lw=.7);twin.set_ylabel('Mean active concentration' if j==3 else '',color='.35');twin.tick_params(labelsize=7,colors='.35');ax.set_zorder(twin.get_zorder()+1);ax.patch.set_visible(False)
    gd=fig.add_gridspec(1,1,left=.065,right=.32,bottom=.075,top=.35);ax=fig.add_subplot(gd[0,0]);exp=csv(DATA/'deterministic_experiments.csv');q=exp[(exp.model=='wave_pinning')&(exp.kind=='amplitude_sweep')]
    # Retained input files spell kind 'amplitude' in some revisions.
    if q.empty:q=exp[(exp.model=='wave_pinning')&exp.kind.str.contains('amplitude')]
    curves=[]
    for i,size in enumerate([20,60]):
        z=q[q['size']==size].sort_values('amplitude');amp=z.amplitude.to_numpy();disp=z.shared_gain.to_numpy()*z.dose.to_numpy();col=shade('wave_pinning',i,2);ax.plot(amp/.1,disp,'o-',ms=3,color=col,label=f'C{size}');g=csv(DATA/'map_amplitudes.csv');g=float(g[(g.model=='wave_pinning')&(g['size']==size)].edge_gain.iloc[0]);xx=np.linspace(0,max(amp),100);duration=float(z.duration.iloc[0]);ax.plot(xx/.1,g*xx*duration,'--',color=col,lw=1);curves.extend(dict(size=size,amplitude_over_reference=float(a/.1),permanent_displacement=float(d)) for a,d in zip(amp,disp))
    pd.DataFrame(curves).to_csv(DATA/'D_linearity_plot.csv',index=False);ax.set(xlabel='Pulse peak / reference input',ylabel='Permanent shared displacement');ax.legend(frameon=False);ax.text(.03,.97,'Solid: measured\nDashed: weak-input slope',va='top',transform=ax.transAxes,fontsize=8)
    ge=fig.add_gridspec(1,4,left=.39,right=.97,bottom=.225,top=.35,wspace=.6);amp=csv(DATA/'map_amplitudes.csv');E=[]
    for j,model in enumerate(MODELS):
        ax=fig.add_subplot(ge[0,j]);ax.set_title(NAMES[model]);q=amp[amp.model==model].sort_values('size');small=q.iloc[0]['size'];gg=q.edge_gain.to_numpy()/q.edge_gain.iloc[0];x=q['size']/small;ax.plot(x,gg,'o',ms=3,color=COLORS[model]);ax.set(xlabel='Size / smallest',ylim=(0,1.08),ylabel='Gain / smallest' if j==0 else '');E.extend(dict(model=model,size=s,relative_size=xx,relative_gain=y) for s,xx,y in zip(q['size'],x,gg))
        if model in cross:
            # Spring physical C includes extension; rescale nominal resting grid using corresponding steady circumference.
            physical=gain[(gain.model==model)&(gain.D==campaign.NATIVE[model])].sort_values('size');C0=float(physical.iloc[0].C);cc=np.linspace(float(physical.C.min()),float(physical.C.max()),100);xx=np.interp(cc,physical.C,physical['size'])/small;curve=(1+C0/cross[model])/(1+cc/cross[model]);
            if model=='spring_ring':
                free=gain[(gain.model==model)&(gain.D==1000)].sort_values('C');ratio=np.interp(cc,free.C,free.gain)/free.iloc[0].gain;ax.plot(xx,curve,'--',color=COLORS[model],lw=1,label='fitted conduction');curve=curve*ratio;ax.plot(xx,curve,'-.',color=COLORS[model],lw=1.3,label='fit × measured rise');ax.legend(fontsize=5.7,frameon=False,loc='lower left')
            else:
                ax.plot(xx,curve,'-' if model=='wave_pinning' else '--',color=COLORS[model],lw=1.3,label='predicted (isolated front)' if model=='wave_pinning' else 'fitted (native gains)');ax.legend(fontsize=5.7,frameon=False,loc='lower left')
    pd.DataFrame(E).to_csv(DATA/'E_native_gain_size.csv',index=False)
    gf=fig.add_gridspec(1,1,left=.39,right=.97,bottom=.075,top=.16);ax=fig.add_subplot(gf[0,0]);q=gain[gain.model=='wave_pinning'];pal=plt.cm.viridis(np.linspace(.05,.95,7))
    for col,D in zip(pal,sorted(q.D.unique())):
        z=q[(q.D==D)&q.defined].sort_values('C');ax.plot(z.C,z.gain,'o-',color=col,ms=3,lw=1.3,label=f'{D:g}')
        lost=q[(q.D==D)&(~q.defined)]
        for r in lost.itertuples():ax.plot(r.C,.06,'x',color=col,ms=5)
    ax.set(xlabel='Circumference C',ylabel='Edge gain G');ax.legend(title='Inactive Dᵥ',title_fontsize=7,fontsize=7,ncol=7,frameon=False,loc='lower right',bbox_to_anchor=(1,1.02),borderaxespad=0)
    for letter,x,y in [('A',.025,.97),('B',.485,.97),('C',.025,.665),('D',.025,.385),('E',.35,.385),('F',.35,.185)]:fig.text(x,y,letter,color='black',weight='bold',fontsize=15)
    fig.text(.065,.688,'A: Wave-pinning C30, native continuous cue; dashed white: withdrawal. B: shades denote size; grey: heading undefined at ≥1 size.',fontsize=9)
    fig.text(.065,.425,'C: f = 0.001, C40 (Spring 2). Colored line: permanent displacement per unit Gaussian-pulse dose; grey: mean profile.',fontsize=9)
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure3_round7.{ext}',dpi=200,facecolor='white')
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    for ax,model in zip(axes,MODELS[:3]):
        q=penalty[(penalty.model==model)&(penalty.D>=1)&(penalty.D<penalty.D_reference)];ax.set_title(NAMES[model]);
        for D,z in q.groupby('D'):ax.plot(z.C_over_D,z.penalty,'o-',ms=3,label=f'{D:g}')
        ax.set(xlabel='C / inactive diffusivity',ylabel='Free-conduction gain / gain − 1',xscale='log',yscale='log');ax.legend(title='D',fontsize=6,frameon=False,ncol=2)
    fig.tight_layout();fig.savefig(OUT/'pilot_conduction_penalty.png',dpi=180);plt.close(fig)
    return gain,cross

def documents(gain,cross,fronts):
    params=read(META/'conduction_law_parameters.json');wp=gain[(gain.model=='wave_pinning')&(gain.D==10)&gain.defined].sort_values('size');pred=fronts['lambda_push']/2/(1+wp.C/cross['wave_pinning']);errors=pred/wp.gain-1
    report=dict(isolated_front=fronts,C_star=cross,WP_front_prediction_relative_errors=[float(x) for x in errors],WP_gain_C60_C20=float(wp.iloc[-1].gain/wp.iloc[0].gain),native_sweep_counts={m:int(sum(gain.model==m)) for m in MODELS[:3]},undefined=gain[~gain.defined].reindex(columns=['model','size','D','reason']).to_dict('records'),source_hashes=SOURCES,producer_state='finished',generated_utc=campaign.source.now())
    write(META/'scientific_results.json',report)
    legend=[f'Fast cytosol removes Wave-pinning’s size-dependent boundary gain; an isolated-front prediction gives its absolute native-ring gains within 6.9–8.0% without fitting them. Native gain declines less in Holmes and Spring. Goryachev’s gain falls strongly, but no conduction-law comparison is made for that model.',
    'A, Native Wave-pinning C30 active-concentration kymographs for a cue overlapping the incumbent’s 10%-contrast contour (edge gap −1) or separated from it (gap 4). Time starts at cue onset; cyan bars show spatial cue support. White dashed lines mark withdrawal after 95% of the target heading change; subsequent single-domain relaxation is included. Native cue peak 0.1 and support 3 are unchanged.',
    'B, Time to turn halfway toward a continuous cue against its edge gap, at C30, 37.8 and 47.6 (Spring resting sizes 2.5, 3.15 and 3.97). Shades denote size within each model. Open triangles with upward arrows are right-censored lower bounds; grey spans mark gaps at which heading was undefined at one or more sizes. These strong-cue turning times are distinct from weak-pulse susceptibility.',
    'C, Susceptibility to a weak Gaussian pulse, permanent displacement per peak×duration dose, against position at C40 (Spring resting size 2). Colored lines use the native translation-mode adjoint at f=0.001 (Spring mean input 0.004); grey curves and right axes show the same stationary mean active profile. The Gaussian SD is 1.5 (Spring 0.15). Signed gain is localized to the two lateral boundaries.',
    'D, Wave-pinning shared permanent boundary displacement against edge-pulse peak divided by reference input 0.1, at C20 and C60; duration 50. Dashed lines use the weak-input gain. Ten-percent departures are bracketed between 5–20% of reference input in Wave-pinning, 2.6–5% in Holmes, 1.3–2.6% in Spring and 2.6–20% in Goryachev (neighboring tested amplitudes, pooled size ranges). Dose-matched duration sweeps agree within 0.2%.',
    f'E, Native edge gain normalized to each model’s smallest size, against relative resting size. Wave-pinning’s solid curve is predicted from the isolated front: G=(λ/2)/(1+C/C*), C*={cross["wave_pinning"]:.3g}, with stationary domain fraction θ and no fitted ring gains. Holmes’s dashed curve is fitted to native inverse gain, 1/G=a+bC, C*=a/b≈{cross["holmes_model3"]:.3g}. Spring’s dashed conduction curve uses C*≈{cross["spring_ring"]:.2g}, fitted to the free-conduction penalty G(1000,C)/G(D,C)−1 against C/D for D≥1; the dash-dot curve multiplies that fit by the measured free-conduction gain rise. This empirical correction and the apparent native crossover (about 140) are distinct from an isolated-front prediction. Goryachev has points only. Spring’s actual circumference includes extension; its native points use mean input 0.004, the sweep 0.0004. Fit-definition sensitivity is reported in the readme.',
    'F, Wave-pinning adjoint edge gain against circumference at inactive diffusivities 2.5, 5, 10 (native), 20, 40, 100 and 1000. Every point uses the same weak mean input and unchanged reaction kinetics. All 35 Wave-pinning preparations retain a polarized branch, including Dᵥ=2.5 at C60 where the pilot reported domain loss. Fast cytosol largely removes the native size dependence. All quantities use model-native units; none is a calibrated migration velocity or biochemical dose.']
    (OUT/'main_legend.md').write_text('\n\n'.join(legend)+'\n')
    readme=f'''# Unexpected

The native Dᵥ=2.5 C60 preparation retains a polarized Wave-pinning domain (G=0.777482, contrast=0.915221, maximum RHS=7.46×10⁻¹³), although the pilot reported domain loss. All tested preparations retain polarized branches. The mechanism remains supported: fast cytosol nearly removes Wave-pinning’s size dependence. Spring also has an intrinsic approximately 5% gain rise with size at fast cytosol; this partly cancels its conduction penalty, so its apparent native inverse-gain crossover cannot be interpreted as a transport length. Both effects remain visible rather than being absorbed into a fit.

## Results and definitions

Boundary gain G is terminal displacement per unit peak-amplitude×duration dose for a periodic Gaussian pulse of SD 1.5 (Spring 0.15). C* is the crossover at which front and cytosolic conduction resistances are equal. Native local reaction kinetics stay fixed. Chemical mean f=0.001. The Spring diffusivity sweep uses mean input 0.0004, as in the pilot; C’s retained spatial maps and E’s retained native points use round6 f=0.001 (mean 0.004). This weak-input difference is disclosed; the kinetics are identical. Goryachev’s native susceptibility curve is retained without fitting a conduction law.

The isolated front solves the Maxwell equal-area condition at weak uniform mean input, then computes front-speed and Gaussian-push coefficients from scalar translation-mode quadratures. Grid refinement from 20001 to 40001 concentration nodes changes either coefficient by at most {fronts['quadrature_relative_change']:.3g}. This is an isolated clamped-cytosol front, not a ring-gain fit. It predicts absolute native-ring gains within {100*abs(errors).min():.2f}–{100*abs(errors).max():.2f}% across C20–60. Its κ={fronts['kappa']:.6f} and λ={fronts['lambda_push']:.6f} differ by {100*(fronts['kappa']/1.2816-1):+.3f}% and {100*(fronts['lambda_push']/5.2958-1):+.3f}% from the reported pilot coefficients. These small differences do not affect the mechanistic interpretation; rounded pilot gain comparisons are retained compactly in pilot_gain_comparison.csv and were never optimization targets.

$$G=\\frac{{\\lambda/2}}{{1+C/C^*}},\\qquad C^*=\\frac{{2D_v}}{{\\kappa\\,\\Delta u\\,\\theta(1-\\theta)}}.\\tag{{3.1}}$$

Wave-pinning C*={cross['wave_pinning']:.6g} (pilot about 76); Holmes C*={cross['holmes_model3']:.6g} (pilot about 400); Spring C*={cross['spring_ring']:.6g} (pilot about 43, range 30–80). Plotted Holmes C* uses the native-D inverse-gain fit (intercept/slope of 1/G versus physical C), and its alternative free-conduction penalty fit gives 469.055. Spring’s mechanistic transport C*=49.0025 uses the through-origin penalty fit over D≥1. Its apparent native inverse-gain crossover is 142.299: the intrinsic small-ring free-conduction gain increase offsets the native resistance decline. This apparent crossover is not a transport length. In E, Spring’s dashed curve shows the fitted conduction factor; its dash-dot curve multiplies this factor by the measured free-conduction baseline G(1000,C)/G(1000,Cmin), which reveals the intrinsic rise rather than absorbing it into C*. These alternatives are fit-definition sensitivity, not uncertainty intervals, and no noisy persistence result is fitted. Per-diffusivity inverse-gain fits, including nonlinearity or negative slopes at nearly free conduction, remain in inverse_gain_fits.csv. Spring’s small-ring free-cytosol increase is retained rather than normalized away in the native plotted points. D=0.1 and 0.3 are reported but excluded from the Spring fit. Even over D≥1, Spring does not exhibit the clean C/D collapse seen in Holmes: its pure conduction factor misses the intrinsic gain rise, and the measured-rise correction is empirical. The fitted Spring description is not an isolated-front prediction or evidence for a universal cross-model law.

## Preparation and numerical verification

All sweeps restart the exact f=0.001 round5 deterministic mean states, except Spring’s mean is changed from 0.004 to 0.0004 as specified in the pilot. BDF tolerances are 2×10⁻¹⁰ relative and 2×10⁻¹² absolute. Chemical relaxations last 16000 native units; Spring 1600. Each condition must retain a unique half-contrast boundary, contrast >0.05 and max RHS <10⁻⁶ before adjoint evaluation. Derivatives include the complete Spring circumference contribution to the altered inactive diffusion. The full Jacobian is differentiated using the same native finite-volume operator. Translation eigenvalues, adjoint residuals, pool-preserving native sources, actual C and stationary residuals accompany every gain.

A reuses exact native cue-on and withdrawal trajectories from the qualified round3 relaxation assay and displays sampled native activity fields, with no new deterministic turning integration. B reuses immutable half-turn brackets. C–D and E’s native points reuse round6 susceptibility, pulse and amplitude-sweep tables. D’s accurate linearity brackets replace the target narrative’s provisional 3–25%. These strong-cue and weak-cue assays use different size grids and do not establish a single timescale bridge.

## Cellular measurements and WGD predictions

Three prospective assays connect these model mechanisms to a matched WGD+/− pair. They are predictions, not completed cellular experiments.

1. Apply photoactivation at offsets from straddling to separated relative to the incumbent leading-edge boundary, matching cue dose and boundary offset between WGD− and WGD+. The models predict longer half-turn times in WGD+; straddling versus separated cues distinguish sliding of the incumbent from challenger replacement. Goryachev does not exhibit the same turning-regime switch. These strong-cue measurements are separate from the weak linear-response assay.
2. Apply a weak, brief photoactivation pulse to one lateral boundary and measure the permanent common displacement of both lateral edges. A matched dose is predicted to move WGD+ boundaries less, markedly in Wave-pinning- and Goryachev-like states and little in Holmes- and Spring-like states. This response contrast identifies a susceptibility class; by itself it does not distinguish conduction from inactive-pool depletion. Pair it with cytosolic-pool mobility and the Figure 2 first-rise assay to test those mechanisms.
3. Measure spontaneous input/activity fluctuations at both lateral edges and in the active plateau separately. In the drift regime, edge fluctuations should account for heading drift, with little direct contribution from plateau input. Provided boundary-response shape, edge co-motion and input statistics are shared, the squared WGD+/− gain ratio predicts the ratio of physical boundary drift. Angular decorrelation additionally includes the perimeter factor:

$$\\frac{{D_{{s,+}}}}{{D_{{s,-}}}}=\\left(\\frac{{G_+}}{{G_-}}\\right)^2,\\qquad\\frac{{k_{{\\mathrm{{drift}},+}}}}{{k_{{\\mathrm{{drift}},-}}}}=\\left(\\frac{{G_+}}{{G_-}}\\right)^2\\left(\\frac{{C_-}}{{C_+}}\\right)^2.\\tag{{3.2}}$$

This conditional weak-drift prediction is supported for Wave-pinning and the eligible Holmes/Spring regimes; it does not explain the full Goryachev heading readout. A jump-dominated Goryachev-like state requires separate first-rise and relocation measurements. Physical boundary drift, model angular heading decorrelation, empirical migration summaries and conditional OU velocity persistence remain distinct observables.

## Data and provenance

`supporting/data/` contains every plotted coordinate, compact native states, the full diffusivity sweeps, isolated front and preserved linearity brackets. `supporting/metadata/scientific_results.json` records source SHA256 hashes, numerical outcomes, pilot discrepancies and completed producer status. The separate pilot_conduction_penalty.png is a check and is not a manuscript panel. No dropped panel was moved to a supplement.

Rebuild with the existing scaling_sandbox environment: run figure3_conduction_20261026.py front; submit its worker --model choices through SLURM, then run figure3_round7_render_20261026.py. Each worker uses one CPU, 4 GB, small QoS, one-hour cap. Sources are fixed during production; no new environment is required.
'''
    (OUT/'scientific_readme.md').write_text(readme)
    section=dict(id='figure3',title='Figure 3. Enlarged rings turn more slowly and respond less to weak boundary input',role='Separates strong-cue turning from linear boundary susceptibility and tests the cytosolic conduction mechanism.',status='round7_complete',legend=legend,narrative=['Larger rings take longer to turn toward a continuous cue across the tested positions. Overlap distinguishes incumbent sliding from nucleation and replacement by a separate domain; further separation has relatively little effect outside Goryachev. These strong-cue responses are distinct from weak-pulse gain.', 'Weak-input translation susceptibility is localized to the two lateral boundaries. Native gain falls substantially with size in Wave-pinning and Goryachev, and changes little in Holmes and Spring.', 'Fast cytosol removes Wave-pinning’s size dependence: the C60/C20 gain ratio rises from about 0.70 at native diffusivity to 1.00 at fast diffusivity, with unchanged reaction kinetics. This perturbation identifies a cytosolic-transport contribution to boundary gain.', 'An isolated-front conduction law predicts absolute native Wave-pinning gains within 6.9–8.0% without fitting ring gains. Holmes is represented by a native inverse-gain fit. Spring’s fitted transport factor is offset by an independently measured free-conduction gain rise, so its corrected curve is empirical and the simple C/D law does not collapse all Spring data.'],measurements=['At matched cue dose and boundary offset in a WGD+/− pair, vary photoactivation from straddling to separated relative to the incumbent boundary and measure half-turn time, incumbent sliding and challenger replacement. The models predict longer WGD+ half-turn times; Goryachev lacks the same overlap-to-separation regime switch. This strong-cue assay is distinct from the linear response.', 'Apply a weak brief pulse at one lateral boundary and measure permanent common displacement of both edges. At matched dose, WGD+ displacement is predicted to decrease strongly for Wave-pinning- and Goryachev-like susceptibility classes, but little for Holmes- and Spring-like classes. This contrast alone cannot separate conduction from pool depletion; combine it with inactive-pool mobility and Figure 2 first-rise measurements.', 'Compare spontaneous lateral-edge and plateau input/activity fluctuations in WGD− and WGD+. Edge fluctuations should dominate weak heading drift. Under shared input statistics, response shape and correlated edge motion, the squared gain ratio predicts physical drift ratio; angular drift additionally includes the squared inverse perimeter ratio (Eq. 3.2). This conditional drift prediction does not establish the Goryachev jump mechanism or a migration-velocity law.'],limitations=['Strong turning and weak response use different amplitudes and size grids; susceptibility does not explain the strong turning times.','Spring’s low-diffusivity branch changes lie outside the reported conduction law. Goryachev has no fitted conduction curve.','Model time, input and size have no cross-model biological calibration.'],equations=[dict(math=r'G=\frac{\lambda/2}{1+C/C^*},\qquad C^*=\frac{2D_v}{\kappa\Delta u\theta(1-\theta)}',number='3.1',description='Wave-pinning crossover uses the isolated front and stationary domain fraction; Holmes uses the native inverse-gain regression and Spring uses the diffusivity-penalty regression, with intrinsic free-conduction rise shown explicitly.'),dict(math=r'\frac{D_{s,+}}{D_{s,-}}=\left(\frac{G_+}{G_-}\right)^2,\qquad\frac{k_{\mathrm{drift},+}}{k_{\mathrm{drift},-}}=\left(\frac{G_+}{G_-}\right)^2\left(\frac{C_-}{C_+}\right)^2',number='3.2',description='Conditional prediction for shared input statistics, response shape and correlated edge motion. The first ratio is physical boundary drift; angular decorrelation additionally includes geometry.')],assets=[dict(path=rel(OUT/'figure3_round7.png'),native_pdf_page=rel(OUT/'figure3_round7.pdf'),page_size=[18,16],minimum_printed_width_inches=18,page_title='')],main_pdf=rel(OUT/'figure3_round7.pdf'),reader=rel(OUT/'scientific_readme.md'),sources=[dict(label='Scientific methods, departures from pilots and numerical scope',path=rel(OUT/'scientific_readme.md')),dict(label='Diffusivity sweep and gains',path=rel(DATA/'conduction_sweep.csv')),dict(label='Front quadrature and fitted crossover parameters',path=rel(META/'conduction_law_parameters.json'))])
    write(OUT/'legend.json',section)
    artifacts={rel(p):sha(p) for p in OUT.rglob('*') if p.is_file() and p.suffix not in ['.log'] and p.name!='provenance.json'}
    write(META/'provenance.json',dict(generated_utc=campaign.source.now(),producer_state='finished',baseline_commit='122d306b',source_hashes=SOURCES,output_hashes=artifacts,code_sha256={rel(Path(__file__)):sha(__file__),rel(Path(campaign.__file__)):sha(campaign.__file__)}))

if __name__=='__main__':
    gain,penalty,cross,fronts,kymos=prepare();render(gain,penalty,cross,fronts,kymos);documents(gain,cross,fronts)
