"""Figure 3 round5: immutable A/B data, deterministic edge response, drift test."""
from pathlib import Path
import sys,os,copy,json,html,textwrap,hashlib
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgb
from pypdf import PdfReader,PdfWriter
import figure3_round5_edge_20261023 as edge
ROOT=edge.ROOT;MS=edge.MS;RUN=edge.RUN;OUT=edge.OUT;DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
OLD=MS/'figures/length_mechanisms_20261008/dev/revision_round3_20261022/figure3_single_C'
DRIFT=MS/'figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4/supporting/data/drift_decomposition.csv'
MODELS=('wave_pinning','holmes_model3','spring_ring','goryachev')
NAMES=dict(wave_pinning='Wave-pinning',holmes_model3='Holmes',spring_ring='Spring',goryachev='Goryachev')
COLORS=dict(wave_pinning='#17619c',holmes_model3='#28783b',spring_ring='#aa7a00',goryachev='#b93434')

def color(model,i,count=5):
    frac=.38+.62*i/max(count-1,1);return tuple(np.array(to_rgb(COLORS[model]))*frac+1-frac)
def read(p):return edge.read(p)
def write(p,v):
    def clean(x):
        if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
        if isinstance(x,np.generic):return x.item()
        return x
    edge.write(p,clean(v))
def rel(p):return str(p.relative_to(ROOT))
def setup():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.4,'axes.labelsize':8.5,'axes.titlesize':9,'xtick.labelsize':7.2,'ytick.labelsize':7.2,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def collect():
    DATA.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    import shutil
    for name in ('activity_profiles.csv.gz','activity_profile_metadata.json','half_heading_readouts.csv','regime_lines.csv','undefined_heading_bands.csv'):
        shutil.copy2(OLD/'supporting/data'/name,DATA/name)
    m=read(RUN/'manifest.json');summaries=[];pulses=[];curves=[];gains=[]
    drift=pd.read_csv(DRIFT);drift=drift[(drift.flicker_frames==5)&(drift.MSD_range_fraction==.25)]
    for c in m['jobs']:
        folder=RUN/'workers'/f"{c['slot']:03d}"
        retry=RUN.with_name('figure3_round5_edge_numeric_retry_20261023')/'workers'/f"{c['slot']:03d}"
        if not read(folder/'status.json').get('complete') and (retry/'status.json').exists() and read(retry/'status.json').get('complete'):folder=retry
        s=read(folder/'status.json')
        if not s.get('complete'):raise ValueError(f"Incomplete required deterministic condition {c['slot']}")
        d=drift[(drift.model==c['model'])&(drift.C_label==c['size'])&(drift.f==c['fraction'])]
        assert len(d)==1
        d=d.iloc[0]
        r=dict(model=c['model'],size=c['size'],fraction=c['fraction'],mean=c['mean'],defined=s.get('defined'),reason=s.get('reason',''),mean_circumference=s.get('mean_circumference'),mean_contrast=s.get('mean_contrast'),mean_stationary_rhs=s.get('mean_stationary_rhs'),mean_tail_scaled_change=s.get('mean_tail_scaled_change'),observed_Ds=d.Ds,observed_Ds_ci_lower=d.Ds_ci_lower,observed_Ds_ci_upper=d.Ds_ci_upper,n_seeds=d.n_seeds,valid_fraction=d.valid_fraction,multiple_front_fraction=d.multiple_front_fraction,observed_MSD_exponent=d.MSD_exponent,observed_MSD_max_lag=d.MSD_max_lag,observed_status=d.status,analysis_lag=d.analysis_lag,source_status=rel(folder/'status.json'))
        if s.get('defined'):
            linear=s['linear_response']
            calibrations=[v['measured_edge_gain']/v['predicted_phase_gain'] for v in linear['edge_gain_validation']]
            DL=linear['D_left']*calibrations[0]**2;DR=linear['D_right']*calibrations[1]**2
            cross=linear['D_left']*calibrations[0]*calibrations[1]
            measured_pred=(DL+DR)/4;full_pred=(DL+DR+2*cross)/4
            r.update(left_measured_gain_scale=calibrations[0],right_measured_gain_scale=calibrations[1],raw_native_Ds_independent=linear['Ds_independent'],amplitude=s['amplitude'],pulse_duration=s['pulse_duration'],observation_horizon=s['observation_horizon'],maximum_linearity_error=s['maximum_linearity_error'],maximum_pulse_gain_error=max(x['relative_error'] for x in linear['edge_gain_validation']),translation_eigenvalue=linear['translation_eigenvalue'],adjoint_residual=linear['adjoint_residual'],Ds_predicted_independent=measured_pred,Ds_predicted_full_covariance=full_pred,D_left=DL,D_right=DR,prediction_over_observation=measured_pred/d.Ds if d.Ds>0 else np.nan,full_covariance_over_observation=full_pred/d.Ds if d.Ds>0 else np.nan)
            with np.load(folder/'pulse_responses.npz') as z:
                t=z['time'];dis=z['displacement'];amps=z['amplitudes'];mults=z['multipliers'];es=z['pulse_edges'];area=z['integrated'];gain=z['terminal_gain'];initial=z['initial_edges'];sens=z['phase_sensitivity']
                for k in range(len(amps)):
                    for tracked in (0,1):
                        pulses.append(dict(model=c['model'],size=c['size'],fraction=c['fraction'],pulse_edge=int(es[k]),tracked_edge=tracked,multiplier=mults[k],amplitude=amps[k],terminal_displacement=dis[k,-1,tracked],terminal_gain=gain[k,tracked],integrated_displacement_per_pulse_dose=area[k,tracked],pulse_duration=s['pulse_duration'],observation_horizon=s['observation_horizon']))
                        if c['fraction']==.001:
                            for i,tt in enumerate(t):curves.append(dict(model=c['model'],size=c['size'],fraction=c['fraction'],pulse_edge=int(es[k]),tracked_edge=tracked,multiplier=mults[k],amplitude=amps[k],time=tt,displacement=dis[k,i,tracked],displacement_per_amplitude=dis[k,i,tracked]/amps[k],displacement_per_pulse_dose=dis[k,i,tracked]/amps[k]/s['pulse_duration']))
                for q in linear['edge_gain_validation']:gains.append(dict(model=c['model'],size=c['size'],fraction=c['fraction'],**q))
        summaries.append(r)
    frame=pd.DataFrame(summaries);frame.to_csv(DATA/'edge_prediction_comparison.csv',index=False)
    pulse_frame=pd.DataFrame(pulses);pulse_frame.to_csv(DATA/'pulse_response_summary.csv',index=False)
    folds=[]
    for model in MODELS:
        q=pulse_frame[(pulse_frame.model==model)&(pulse_frame.fraction==.001)&(pulse_frame.pulse_edge==1)&(pulse_frame.tracked_edge==1)&(pulse_frame.multiplier==1)].sort_values('size')
        first,last=q.iloc[0],q.iloc[-1]
        folds.append(dict(model=model,fraction=.001,small_size=first['size'],large_size=last['size'],amplitude=first.amplitude,pulse_duration=first.pulse_duration,terminal_gain_small=first.terminal_gain,terminal_gain_large=last.terminal_gain,terminal_gain_large_over_small=last.terminal_gain/first.terminal_gain,integrated_response_large_over_small=last.integrated_displacement_per_pulse_dose/first.integrated_displacement_per_pulse_dose))
    pd.DataFrame(folds).to_csv(DATA/'edge_response_size_folds.csv',index=False)
    pd.DataFrame(curves).to_csv(DATA/'edge_response_curves.csv.gz',index=False)
    pd.DataFrame(gains).to_csv(DATA/'kernel_pulse_validation.csv',index=False)
    # Exact source hashes establish preservation of A/B quantitative data.
    preservation={name:dict(previous=rel(OLD/'supporting/data'/name),sha256=edge.sha(DATA/name),identical=edge.sha(DATA/name)==edge.sha(OLD/'supporting/data'/name)) for name in ('activity_profiles.csv.gz','activity_profile_metadata.json','half_heading_readouts.csv','regime_lines.csv','undefined_heading_bands.csv')}
    write(META/'preserved_AB.json',dict(generated_utc=edge.now(),files=preservation,all_identical=all(x['identical'] for x in preservation.values())))
    return frame,pd.DataFrame(pulses),pd.DataFrame(curves),pd.DataFrame(gains)

def panel_A(fig,grid):
    p=pd.read_csv(DATA/'activity_profiles.csv.gz');metadata=read(DATA/'activity_profile_metadata.json')
    axes=[]
    for j,(label,gap) in enumerate((('i',-1),('ii',4))):
        ymax=p[p.label==label].active.max()*1.1
        for stage in range(4):
            ax=fig.add_subplot(grid[stage,j]);axes.append(ax)
            m=next(x for x in metadata if x['label']==label and x['stage']==stage);q=p[(p.label==label)&(p.stage==stage)].sort_values('x')
            ax.axvspan(m['cue_left'],m['cue_right'],color='.87');ax.plot(q.x,q.active,color=COLORS['wave_pinning'],lw=1.3)
            for saddle in m['saddles']:ax.plot(saddle,np.interp(saddle,q.x,q.active),'D',mfc='white',mec='black',ms=3.5)
            ax.plot(m['heading'],0,'|',color='black',ms=8,mew=1.4,clip_on=False)
            ax.set_ylim(0,ymax);ax.set_xlim(-15,15);ax.tick_params(length=2.5,pad=2)
            if stage==0:ax.set_title(f'({label})  gap {gap:g}\n'+('straddling' if gap<0 else 'separated'),pad=7)
            if stage<3:ax.tick_params(labelbottom=False)
            else:ax.set_xlabel('Arc position')
            if j==0:ax.set_ylabel('Activity, u')
            ax.text(.96,.88,'cue off' if stage==3 else f'{int(m["fraction"]*100)}%',transform=ax.transAxes,ha='right',va='top',fontsize=7.5)
    return axes

def panel_B(fig,grid):
    f=pd.read_csv(DATA/'half_heading_readouts.csv');regimes=pd.read_csv(DATA/'regime_lines.csv');bands=pd.read_csv(DATA/'undefined_heading_bands.csv');axes=[]
    mapping=dict(wave_pinning='wave_pinning',holmes_model3='holmes_model3',spring_ring='spring',goryachev='goryachev')
    for j,model in enumerate(MODELS):
        key=mapping[model];ax=fig.add_subplot(grid[0,j]);axes.append(ax);ax.set_title(NAMES[model]);ax.set_yscale('log')
        for r in bands[bands.model==key].itertuples():ax.axvspan(r.left,r.right,color='.92',zorder=-3,lw=0)
        sizes=sorted(f[f.model==key].size_parameter.unique())
        for i,size in enumerate(sizes):
            g=f[(f.model==key)&(f.size_parameter==size)&(~f.no_turn)].sort_values('gap');rg=regimes[(regimes.model==key)&(regimes.size_parameter==size)].iloc[0];col=color(model,i,3);style=('-','--',':')[i]
            if pd.notna(rg.first_saddle_gap):ax.axvline(rg.first_saddle_gap,color=col,ls=style,lw=1.8-.4*i,zorder=-1)
            ax.plot(g.gap,g.t50,color=col,ls=style,lw=1.3,marker='o',ms=2.4)
            for r in g[g.status=='right_censored_saved_heading_at_completed_horizon'].itertuples():ax.plot(r.gap,r.lower,'^',mfc='white',mec=col,ms=5)
            valid=g[g.t50.notna()]
            if len(valid):
                r=valid.iloc[-1];ax.annotate(f'C={size:g}',(r.gap,r.t50),xytext=(4,(-8,0,9)[i]),textcoords='offset points',color=col,fontsize=6.6,clip_on=False)
            if model=='wave_pinning' and size==30:
                for label,gap in (('i',-1),('ii',4)):
                    r=g[np.isclose(g.gap,gap)].iloc[0];ax.plot(r.gap,r.t50,'o',ms=5,mfc=col,mec='black',mew=.5);ax.annotate(label,(r.gap,r.t50),xytext=(0,-11),textcoords='offset points',ha='center',fontsize=7)
        g=f[(f.model==key)&(~f.no_turn)];span=g.gap.max()-g.gap.min();ax.set_xlim(g.gap.min()-.025*span,g.gap.max()+.14*span);lo,hi=ax.get_ylim();ax.set_ylim(lo/1.2,hi*1.4);ax.set_xlabel('Cue edge gap')
        if j==0:ax.set_ylabel(r'Half-turn time, $t_{50}$')
    return axes

def main_figure(frame,pulses,curves):
    setup();fig=plt.figure(figsize=(16.2,12.9))
    ga=fig.add_gridspec(4,2,left=.055,right=.345,bottom=.57,top=.92,wspace=.28,hspace=.20)
    gb=fig.add_gridspec(1,4,left=.425,right=.955,bottom=.67,top=.92,wspace=.57)
    gc=fig.add_gridspec(1,4,left=.07,right=.955,bottom=.33,top=.50,wspace=.40)
    gd=fig.add_gridspec(1,4,left=.07,right=.955,bottom=.075,top=.245,wspace=.40)
    panel_A(fig,ga);panel_B(fig,gb)
    fig.text(.017,.962,'A',fontsize=13,fontweight='bold');fig.text(.385,.962,'B',fontsize=13,fontweight='bold');fig.text(.025,.53,'C',fontsize=13,fontweight='bold');fig.text(.025,.275,'D',fontsize=13,fontweight='bold')
    fig.legend(handles=[Patch(facecolor='.87',label='cue support'),Line2D([],[],marker='D',mfc='white',mec='black',ls='none',label='saddle'),Line2D([],[],marker='|',color='black',ls='none',label='heading')],loc='upper left',bbox_to_anchor=(.052,.982),frameon=False,ncol=3,fontsize=7,handlelength=1)
    fig.legend(handles=[Patch(facecolor='.92',label='heading undefined at ≥1 size'),Line2D([],[],ls=':',color='.25',label='first saddle at 50% (C style)'),Line2D([],[],marker='^',mfc='white',mec='black',ls='none',label='time lower bound'),Line2D([],[],marker='o',mfc=COLORS['wave_pinning'],mec='black',ls='none',label='i–ii: A placements')],loc='upper left',bbox_to_anchor=(.421,.635),frameon=False,ncol=2,fontsize=7.4,columnspacing=1.3,handlelength=1.6)
    fig.text(.425,.574,'B: C = 30, 37.8, 47.6; Spring 2.5, 3.15, 3.97.\nC–D: Figure 4 sizes; native units within each model.',fontsize=8,color='black',va='top')
    markers=['o','s','^','v','D','P']
    for j,model in enumerate(MODELS):
        ax=fig.add_subplot(gc[0,j]);ax.set_title(NAMES[model],pad=22);ax.axvspan(0,5 if model=='spring_ring' else 50,color='.94',zorder=0)
        q=curves[(curves.model==model)&(curves.pulse_edge==1)&(curves.tracked_edge==1)&(curves.multiplier==1)]
        sizes=sorted(q['size'].unique()); lines=[]
        for i,size in enumerate(sizes):
            z=q[q['size']==size];ax.plot(z.time,z.displacement_per_amplitude,color=color(model,i),lw=1.35,label=f'{size:g}');lines.append(Line2D([],[],color=color(model,i),label=f'{size:g}'))
        ax.set_xlabel('Time from pulse onset');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3));ax.tick_params(length=3,pad=2)
        if j==0:ax.set_ylabel(r'Pulsed-edge displacement / $a$')
        ax.set_ylim(0,q.displacement_per_amplitude.max()*2.05)
        ax.text(-.01,1.063,'C',transform=ax.transAxes,ha='right',fontsize=6.5)
        ax.legend(handles=lines,frameon=False,fontsize=6.4,ncol=5,loc='upper left',bbox_to_anchor=(0,1.105),handlelength=.85,columnspacing=.7,handletextpad=.3,labelspacing=.2,borderpad=.2)
        inset=ax.inset_axes([.50,.53,.47,.43]);g=pulses[(pulses.model==model)&(pulses.fraction==.001)&(pulses.pulse_edge==1)&(pulses.tracked_edge==1)&(pulses.multiplier==1)].sort_values('size')
        inset.plot(g['size'],g.integrated_displacement_per_pulse_dose,color=COLORS[model],lw=1,marker='o',ms=2.5);inset.tick_params(labelsize=5.7,length=2,pad=1);inset.set_xlabel('C',fontsize=6,labelpad=1);inset.set_ylabel(r'$\mathcal{A}_{20\tau}$',fontsize=7,labelpad=1);inset.ticklabel_format(axis='y',style='sci',scilimits=(-2,3));inset.yaxis.offsetText.set_fontsize(5)
        dx=fig.add_subplot(gd[0,j]);dx.set_title(NAMES[model]);dx.set_xscale('log');dx.set_yscale('log');g=frame[(frame.model==model)&(frame.defined)&(frame.observed_Ds>0)]
        lo=min(g.Ds_predicted_independent.min(),g.observed_Ds.min())*.4;hi=max(g.Ds_predicted_independent.max(),g.observed_Ds.max())*2
        dx.plot([lo,hi],[lo,hi],color='.35',lw=.9,zorder=0)
        for i,size in enumerate(sorted(frame[frame.model==model]['size'].unique())):
            for k,f in enumerate(edge.FRACTIONS):
                z=g[(g['size']==size)&(g.fraction==f)]
                if not len(z):continue
                r=z.iloc[0];col=color(model,i)
                if pd.notna(r.observed_Ds_ci_lower) and r.observed_Ds_ci_lower>0:dx.plot([r.observed_Ds_ci_lower,r.observed_Ds_ci_upper],[r.Ds_predicted_independent]*2,color=col,lw=.55,zorder=1)
                dx.plot(r.observed_Ds,r.Ds_predicted_independent,marker=markers[k],ms=4.3,mfc=col,mec=COLORS[model],mew=.35,ls='none')
        dx.set_xlim(lo,hi);dx.set_ylim(lo,hi);dx.set_xlabel(r'Measured $D_s$ (Fig. 4)');dx.tick_params(length=3,pad=2)
        if j==0:dx.set_ylabel(r'Predicted $D_s$, independent edges')
        if model=='goryachev':dx.text(.02,.98,'2 undefined mean states',transform=dx.transAxes,va='top',fontsize=6.8)
    fig.text(.07,.297,r'C: $f=0.001$; $a=10^{-4}$ × native reference; shaded interval: pulse on. Inset: displacement integrated to $20\tau$, divided by pulse dose $a\tau$.',fontsize=8)
    fig.legend(handles=[Line2D([],[],marker=m,mfc='.65',mec='.25',ls='none',label=f'{f:g}') for f,m in zip(edge.FRACTIONS,markers)],title='D: f',title_fontsize=7.5,loc='lower left',bbox_to_anchor=(.065,.013),ncol=6,frameon=False,fontsize=7.4,handlelength=1.3,columnspacing=1.5)
    fig.text(.70,.025,'Line: 1:1; horizontal bars: seed-bootstrap 95% CI.',fontsize=7.4)
    for ext in ('pdf','png','svg'):fig.savefig(OUT/f'figure3_round5.{ext}',dpi=220)
    plt.close(fig)

def page_text(pdf,title,paragraphs):
    fig=plt.figure(figsize=(11.7,8.3));fig.text(.06,.94,title,fontsize=13,fontweight='bold',va='top');y=.87
    for text in paragraphs:
        for line in textwrap.wrap(text,135,break_long_words=False,break_on_hyphens=False):fig.text(.06,y,line,fontsize=9,va='top');y-=.026
        y-=.02
    pdf.savefig(fig);plt.close(fig)

def supplement(frame,pulses,curves,gains):
    setup();path=OUT/'figure3_round5_new_supplement.pdf'
    indep=pd.read_csv(DATA/'edge_independence_summary.csv');msds=pd.read_csv(DATA/'edge_MSD_curves.csv')
    with PdfPages(path) as pdf:
        page_text(pdf,'Figure 3 supplementary methods: linear response and the neutral mode',[
            'Every deterministic condition starts from its exact Figure 4 basal checkpoint. Uniform mean input is applied without noise and allowed to relax for 16,000 chemical time units or 1,600 Spring time units. A half-contrast edge must be unique, contrast at least 0.05, maximum native RHS no greater than 1e-6, and final 10% preparation change no greater than 1e-4 of the maximum state scale. These gates define a stationary polarized response; failure is recorded rather than replaced by another preparation.',
            'The Gaussian pulse is centered at the interpolated rising or falling half-contrast crossing; its spatial standard deviation equals the declared noise covariance length ell (1.5 chemical, 0.15 Spring), and its duration is tau (50 chemical, 5 Spring). Peak a is 1e-4 of the unchanged native reference: 1e-5 Wave-pinning, 1e-4 Holmes, 4e-4 Spring, 6.4e-4 Goryachev. Each edge receives independent restarts at a/2, a and 2a. The other edge is tracked as well. A no-pulse continuation is subtracted at matching times. Spring edge locations use material fractions times mean-input circumference to remove uniform expansion.',
            'A translated polarized profile is a neutral state. A brief pulse can therefore leave a persistent positional offset after the transient has relaxed: displacement does not have to return to zero. The plotted response is displacement divided by pulse amplitude. The inset is the time integral of that displacement through the explicitly finite 20-tau window, divided by the pulse dose a*tau. It is not an infinite-time integral. Terminal displacement divided by a*tau is the zero-frequency velocity gain.',
            'A single spatial pulse does not identify susceptibility to every location on a ring. The additional spatial response shape is calculated from the unchanged native Jacobian at that same stationary mean-input state. Each edge kernel is normalized to unit gain under that exact pulse shape and multiplied by its measured terminal displacement per pulse dose. The unscaled native gain is separately validated against the measurement. No measured diffusion is used in this calculation. Both edges have the same zero-frequency translation susceptibility; transient shape/width modes need not be the same. The main prediction deliberately applies the requested independent-edge formula to those marginal edge diffusivities.'])
        fig=plt.figure(figsize=(11.7,8.3));axes=fig.subplots(2,2);y=.91
        equations=[
            (r'$\delta\dot{\mathbf{y}}=J\,\delta\mathbf{y}+B\,\delta\mathbf{u},\quad J\mathbf{q}=0,\quad J^{\mathsf{T}}\mathbf{p}=0,\quad\mathbf{p}^{\mathsf{T}}\mathbf{q}=1.$','(1)'),
            (r'$\dot s=\sum_j g_j\,\delta u_j,\qquad g_j=(\mathbf{p}^{\mathsf{T}}B)_j.$','(2)'),
            (r'$K_u(j,k;t,t+v)=\mu^2\{\exp[\sigma^2 C_x(j-k)C_t(t,t+v)]-1\}.$','(3)'),
            (r'$D_e=\sum_{j,k}g_j g_k\int_0^\infty\overline{K_u(j,k;t,t+v)}^{\,t}\,\mathrm{d}v.$','(4)'),
            (r'$D_s^{\mathrm{ind}}=\frac{D_L+D_R}{4},\qquad D_s=\frac{D_L+D_R+2D_{LR}}{4}.$','(5)'),
            (r'$\mathcal{A}_{20\tau}=\frac{1}{a\tau}\int_0^{20\tau}[x_e^{\mathrm{pulse}}(t)-x_e^{\mathrm{control}}(t)]\,\mathrm{d}t.$','(6)'),(r'$\widehat g_{e,j}=\frac{M_e}{Q_e}g_j,\quad M_e=\frac{\Delta x_e(20\tau)}{a\tau},\quad Q_e=\sum_jg_jh_{e,j}.$','(7)')]
        plt.close(fig);fig=plt.figure(figsize=(11.7,8.3));fig.text(.06,.95,'Linear-response definitions and exact imposed covariance',fontsize=13,fontweight='bold')
        yy=.86
        for equation,number in equations:fig.text(.5,yy,equation,ha='center',fontsize=14);fig.text(.94,yy,number,ha='right',fontsize=10);yy-=.095
        fig.text(.06,.07,'q: derivative of the translated native profile; p: its adjoint, excluding conserved-pool modes.\nCx: discrete circular filter autocorrelation. Ct: analytically variance-normalized OU-knot interpolation.\nThe overbar averages the knot phase. Native sigma = 1; lognormal input covariance is not Gaussian OU covariance.',fontsize=9,va='bottom')
        pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,4,figsize=(14,7.5));
        for j,model in enumerate(MODELS):
            g=frame[(frame.model==model)&(frame.defined)]
            for i,size in enumerate(sorted(g['size'].unique())):
                z=g[g['size']==size].sort_values('fraction');axes[0,j].plot(z.fraction,z.maximum_linearity_error,color=color(model,i),marker='o',ms=3)
            axes[0,j].set_xscale('log');axes[0,j].set_yscale('log');axes[0,j].set_title(NAMES[model]);axes[0,j].axhline(.05,color='.5',ls=':');axes[0,j].set_xlabel('f')
            q=gains[gains.model==model];axes[1,j].scatter(q.predicted_phase_gain,q.measured_edge_gain,s=12,color=COLORS[model]);v=np.r_[q.predicted_phase_gain,q.measured_edge_gain];lo,hi=v.min(),v.max();axes[1,j].plot([lo,hi],[lo,hi],color='.5',lw=.8);axes[1,j].set_xlabel('Native spatial-kernel gain')
        axes[0,0].set_ylabel('Max response scaling error');axes[1,0].set_ylabel('Measured terminal gain');fig.subplots_adjust(left=.07,right=.98,bottom=.10,top=.90,wspace=.36,hspace=.4);fig.text(.05,.96,'Amplitude checks and pulse validation; neither uses observed drift',fontsize=12,fontweight='bold');pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,4,figsize=(14,7.5))
        for j,model in enumerate(MODELS):
            q=curves[(curves.model==model)&(curves.multiplier==1)]
            for pulse_edge in (0,1):
                ax=axes[pulse_edge,j];ax.set_title(NAMES[model]+(' · rising-edge pulse' if pulse_edge==0 else ' · falling-edge pulse'))
                for i,size in enumerate(sorted(q['size'].unique())):
                    for tracked_edge in (0,1):
                        z=q[(q['size']==size)&(q.pulse_edge==pulse_edge)&(q.tracked_edge==tracked_edge)];ax.plot(z.time,z.displacement_per_amplitude,color=color(model,i),ls='-' if tracked_edge==pulse_edge else '--',lw=1)
                ax.set_xlabel('Native time');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3))
        for ax in axes[:,0]:ax.set_ylabel('Edge displacement / a')
        fig.subplots_adjust(left=.07,right=.98,bottom=.12,top=.90,wspace=.35,hspace=.4);fig.text(.05,.96,'Both half-contrast edges move after a local pulse; f = 0.001',fontsize=12,fontweight='bold');fig.text(.07,.03,'Solid: pulsed edge; dashed: opposite edge. Each shade is a Figure 4 size. A persistent common displacement is the neutral translation mode.',fontsize=9);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,4,figsize=(14,7.5))
        for j,model in enumerate(MODELS):
            for k,kind in enumerate(('nearest_tau','quarter_horizon')):
                ax=axes[k,j];q=indep[(indep.model==model)&(indep.lag_kind==kind)];
                for i,size in enumerate(sorted(q['size'].unique())):
                    z=q[q['size']==size];resolved=z.resolved_above_quantization
                    ax.scatter(z.fraction[resolved],z.increment_correlation[resolved],s=18,color=color(model,i),edgecolors=COLORS[model],linewidths=.25)
                    ax.scatter(z.fraction[~resolved],z.increment_correlation[~resolved],s=18,facecolors='none',edgecolors=color(model,i),linewidths=.7)
                ax.axhline(0,color='.5',lw=.7);ax.set_xscale('log');ax.set_ylim(-1.05,1.05);ax.set_xlabel('f');ax.set_title(NAMES[model]+(' · lag ≈ τ' if k==0 else ' · lag ¼ window'))
        for ax in axes[:,0]:ax.set_ylabel('Edge increment correlation')
        fig.subplots_adjust(left=.07,right=.98,bottom=.13,top=.90,wspace=.35,hspace=.4);fig.text(.05,.96,'Independence check: Figure 4 exemplar seeds 0 and 1',fontsize=12,fontweight='bold');fig.text(.07,.03,'Filled: both edge MSDs >10× conservative float16 crossing-error bound²; open: saved-profile precision may dominate.\nEach point is a seed/size/input, only uninterrupted single-region pairs. These two exemplars do not give campaign confidence intervals.',fontsize=9);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(1,4,figsize=(14,4.3))
        for j,model in enumerate(MODELS):
            ax=axes[j];q=msds[(msds.model==model)&np.isclose(msds.fraction,.001)]
            for i,size in enumerate(sorted(q['size'].unique())):
                z=q[(q['size']==size)&(q.replicate==0)].sort_values('lag');ax.loglog(z.lag,z.left_MSD,color=color(model,i),lw=1.2);ax.loglog(z.lag,z.right_MSD,color=color(model,i),ls='--',lw=1.2)
            ax.set_title(NAMES[model]);ax.set_xlabel('Native lag')
        axes[0].set_ylabel('Half-contrast edge MSD');fig.subplots_adjust(left=.07,right=.98,bottom=.19,top=.83,wspace=.35);fig.text(.05,.95,'Edge MSDs at f = 0.001; solid rising, dashed falling, seed 0',fontsize=12,fontweight='bold');fig.text(.07,.04,'All seeds, fractions, pair counts, quantization bounds and covariance factors are provided in the plotted-data tables.',fontsize=9);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(1,4,figsize=(14,4.5))
        for j,model in enumerate(MODELS):
            ax=axes[j];g=frame[(frame.model==model)&(frame.defined)&(frame.observed_Ds>0)];ax.set_xscale('log');ax.set_yscale('log')
            lo=min(g.Ds_predicted_full_covariance.min(),g.observed_Ds.min())*.4;hi=max(g.Ds_predicted_full_covariance.max(),g.observed_Ds.max())*2;ax.plot([lo,hi],[lo,hi],color='.5',lw=.9)
            for i,size in enumerate(sorted(g['size'].unique())):
                z=g[g['size']==size];ax.scatter(z.observed_Ds,z.Ds_predicted_full_covariance,s=15,color=color(model,i),edgecolors=COLORS[model],linewidths=.25)
            ax.set_xlim(lo,hi);ax.set_ylim(lo,hi);ax.set_title(NAMES[model]);ax.set_xlabel('Measured tracked drift Ds')
        axes[0].set_ylabel('Full-covariance linear response');fig.subplots_adjust(left=.07,right=.98,bottom=.2,top=.83,wspace=.35);fig.text(.05,.95,'Full covariance sensitivity: retain the shared neutral translation',fontsize=12,fontweight='bold');fig.text(.07,.04,'No parameter fit: the full correlated prediction is approximately twice the prescribed independent-edge prediction because the two measured-gain kernels nearly coincide.\nThis is a prediction of physical tracked drift, not an explanation of all heading decay, especially for Goryachev.',fontsize=9);pdf.savefig(fig);plt.close(fig)
    # Preserve the entire preceding Figure3 main and supplement, after new pages.
    w=PdfWriter();w.append(path);w.append(OLD/'figure3_round3.pdf');w.append(OLD/'figure3_round3_supplement.pdf');w.write(OUT/'figure3_round5_supplement.pdf')
    write(META/'supplement_preservation.json',dict(generated_utc=edge.now(),new_pages=len(PdfReader(path).pages),preserved=[dict(path=rel(p),sha256=edge.sha(p),pages=len(PdfReader(p).pages)) for p in (OLD/'figure3_round3.pdf',OLD/'figure3_round3_supplement.pdf')],total_pages=len(PdfReader(OUT/'figure3_round5_supplement.pdf').pages)))

def documents(frame):
    previous=read(OLD/'legend.json');indep=pd.read_csv(DATA/'edge_independence_summary.csv');stats={}
    for model in MODELS:
        g=frame[(frame.model==model)&(frame.defined)&(frame.observed_Ds>0)]
        q=indep[(indep.model==model)&(indep.lag_kind=='quarter_horizon')&indep.resolved_above_quantization]
        stats[model]=dict(n_prediction_pairs=len(g),independent_ratio_min=float(g.prediction_over_observation.min()),independent_ratio_median=float(g.prediction_over_observation.median()),independent_ratio_max=float(g.prediction_over_observation.max()),full_ratio_median=float(g.full_covariance_over_observation.median()),max_linearity_error=float(g.maximum_linearity_error.max()),max_kernel_pulse_error=float(g.maximum_pulse_gain_error.max()),late_edge_correlation_median=float(q.increment_correlation.median()),late_edge_correlation_n=len(q))
    missing=frame[~frame.defined].to_dict('records');write(META/'results.json',dict(generated_utc=edge.now(),stats=stats,undefined_conditions=missing))
    title='Local edge pulses expose coupled boundary motion and test an untuned prediction of polarity-domain drift.'
    b=previous['legend'][1].split(' C. Mean angular speed')[0]+'. Time and length retain each model’s native units. Curves, gaps, brackets and three Figure 3 sizes are unchanged; lighter solid, intermediate dashed and darker dotted lines denote increasing C. Filled i–ii markers identify A.'
    legend=[previous['legend'][0],b,previous['legend'][2].replace('and matching open down-triangles mark speed upper bounds; ','').replace('Open up-triangles mark censored time lower bounds, ','Open up-triangles mark censored time lower bounds; '),previous['legend'][3],
        'C. Control-subtracted displacement of the pulsed falling half-contrast edge, divided by peak input a, at f = 0.001. Shading denotes Figure 4 circumference: chemical C = 20, 30, 40, 50, 60; Spring resting C = 1, 1.5, 2, 2.5, 3. The stationary polarized state at each declared constant mean is used, without noise or retuning. A Gaussian pulse centered on the boundary has spatial SD ℓ = 1.5 and duration τ = 50 (Spring 0.15 and 5). Peak a is 10⁻⁴ of the native reference: 10⁻⁵ Wave-pinning, 10⁻⁴ Holmes, 4 × 10⁻⁴ Spring and 6.4 × 10⁻⁴ Goryachev. The grey interval is pulse-on. Inset A₂₀τ is displacement integrated through 20τ divided by pulse dose aτ. Both edges, a/2 and 2a are checked in the supplement. A residual positional offset is allowed because translation is neutral; this is not an infinite displacement integral. Largest/smallest terminal displacement per pulse dose is 0.701 Wave-pinning, 0.914 Holmes, 0.989 Spring and 0.0878 Goryachev.',
        'D. The prescribed independent-edge prediction, Dₛ = (Dleft + Dright)/4, against physical tracked domain drift measured in Figure 4 at f ≤ 0.01 (markers). Colors and shading match C. Each measured terminal edge gain per pulse dose scales the native spatial response shape, normalized to give unit response to that exact localized pulse; this extra spatial-shape assumption is explicit. The native linearization supplies its shape only, while measured C amplitudes supply its gain; the imposed discrete spatial covariance and phase-averaged OU-knot interpolation enter the exact lognormal covariance μ²[exp(σ²CₓCₜ) − 1], with σ = 1. There is no diffusion fitting or parameter adjustment. Line: 1:1; horizontal bars: Figure 4 whole-seed-bootstrap 95% CI. This predicts Figure 4’s tracked drift, not total heading decay. The independent prediction/observation medians are '+', '.join(NAMES[k]+' '+f"{stats[k]['independent_ratio_median']:.3f}" for k in MODELS)+'. Goryachev spans 0.00359–4.01 across conditions, so a single bias does not summarize its failure. Goryachev C20/C30 at f = 0.01 has no polarized boundary at the declared mean, so its prediction is undefined.',
        'The edge-independence supplement measures simultaneous half-contrast displacements from existing Figure 4 exemplar profiles, seeds 0 and 1, using uninterrupted single-region intervals and exposing float16 precision limits. At quarter-window lag, resolved median correlations are '+', '.join(NAMES[k]+' '+f"{stats[k]['late_edge_correlation_median']:.3f}" for k in MODELS)+'. Positive common motion violates edge independence. Retaining the native shared translation covariance doubles the prediction; this sensitivity is supplementary and does not adjust any parameter. Goryachev’s disagreement and short valid-region exposure remain explicit.',
        'Cellular counterpart of C: photoactivation straddling a leading-edge boundary, tracking both flanking boundaries and their displacement after withdrawal. Native units are not calibrated between models. deBelly is excluded from rings.']
    section=dict(id='figure3',title=title,status='round5_complete',legend=legend,role='Cue placement separates sliding from takeover; local boundary response tests a no-retuning prediction of physical domain drift.',narrative=['Cue placement relative to the incumbent boundary separates sliding from takeover, and enlargement slows half-turn reorientation. In the weak-input deterministic response the two boundaries share a long-lived translation: the prescribed independent-edge approximation gives median prediction/observation ratios 0.511, 0.540 and 0.534 in Wave-pinning, Holmes and Spring. Existing noise exemplars confirm correlated boundary motion. Goryachev does not give a general drift-based explanation of heading persistence.'],measurements=['Photoactivate across a leading-edge half-contrast boundary; measure both boundary responses, fluctuations and covariance.'],limitations=['The local pulse alone does not identify spatial susceptibility; the unchanged native Jacobian supplies it and is separately checked against both terminal pulse gains. Independence is a tested approximation and fails at long lag. Exemplar edge correlations use two seeds and bounded float16 profiles; unresolved short-lag displacements are flagged. Goryachev C20/C30 f=.01 has no stationary polarized boundary at constant mean. Physical tracked drift is conditional on trackable single-region intervals and is distinct from global heading decay.'],assets=[dict(path=rel(OUT/'figure3_round5.png'),native_pdf_page=rel(OUT/'figure3_round5.pdf'),page_size=[16.2,12.9],minimum_printed_width_inches=16.2,page_title='')],main_pdf=rel(OUT/'figure3_round5.pdf'),supplement_pdf=rel(OUT/'figure3_round5_supplement.pdf'),reader=rel(OUT/'index.html'),sources=[dict(label='Preserved A/B source tables',path=rel(META/'preserved_AB.json')),dict(label='Deterministic frozen protocol',path=rel(RUN/'manifest.json')),dict(label='Prediction and measured drift',path=rel(DATA/'edge_prediction_comparison.csv')),dict(label='Both edges and amplitude checks',path=rel(DATA/'pulse_response_summary.csv')),dict(label='Exemplar edge independence',path=rel(DATA/'edge_independence_summary.csv')),dict(label='Methods and results',path=rel(OUT/'METHODS_RESULTS.md'))])
    section['equations']=[
        dict(math=r'D_s^{\mathrm{ind}}=\frac{D_L+D_R}{4},\qquad D_s=\frac{D_L+D_R+2D_{LR}}{4}',number='3.1',description='The main panel applies edge independence; the supplement retains the measured positive boundary coupling.'),
        dict(math=r'\widehat g_{e,j}=\frac{M_e}{Q_e}g_j,\qquad M_e=\frac{\Delta x_e(20\tau)}{a\tau},\quad Q_e=\sum_j g_jh_{e,j}',number='3.2',description='Measured terminal displacement per pulse dose supplies gain. The unchanged native linearization supplies only the spatial shape g; h is the actual localized pulse shape. No measured diffusion enters.'),
        dict(math=r'K_u(j,k;t,t+v)=\mu^2\{\exp[\sigma^2C_x(j-k)C_t(t,t+v)]-1\}',number='3.3',description='This is the declared lognormal covariance, with the exact discrete circular spatial kernel and variance-normalized OU-knot interpolation; the knot phase is averaged.'),
        dict(math=r'D_e=\sum_{j,k}\widehat g_{e,j}\widehat g_{e,k}\int_0^\infty\overline{K_u(j,k;t,t+v)}^{\,t}\,\mathrm{d}v',number='3.4',description='Native linear response and the convention MSD = 2Dt give each marginal edge diffusivity. The time integral is evaluated without generating new stochastic inputs.')]
    write(OUT/'legend.json',section);(OUT/'legend.md').write_text('**Figure 3. '+title+'**\n\n'+'\n\n'.join(legend)+'\n')
    methods='''# Figure 3: deterministic boundary response and no-retuning drift\n\nGenerated '''+edge.now()+'''; deterministic and profile-analysis producers complete.\n\nFigures A/B retain exactly the preceding quantitative source files, with source SHA-256 checks in `supporting/metadata/preserved_AB.json`. C uses Figure 4’s size grid, not B’s three sizes. No native parameters or noise realization changed.\n\n'''
    methods+='\n\n'.join(legend[4:])+'\n\n'
    methods+='''## Derivation\n\nThe supplementary PDF renders the numbered equations. For a stationary native state, linear perturbations satisfy\n\n$$\n\\delta\\dot{\\mathbf{y}}=J\\delta\\mathbf{y}+B\\delta\\mathbf{u},\\qquad J\\mathbf{q}=0,\\quad J^{\\mathsf{T}}\\mathbf{p}=0,\\quad\\mathbf{p}^{\\mathsf{T}}\\mathbf{q}=1. \\tag{1}\n$$\n\nHere **q** is the derivative of the translated complete profile. Conserved-pool adjoints have zero overlap with it. The native spatial shape for the velocity response to input cell j is\n\n$$\n\\dot{s}=\\sum_j g_j\\delta u_j,\\qquad g_j=(\\mathbf{p}^{\\mathsf{T}}B)_j. \\tag{2}\n$$\n\nThe native gain under the actual edge-pulse shape $h_e$ is $Q_e=\\sum_j g_jh_{e,j}$. The measured finite-window terminal displacement per dose is $M_e$. The kernel used for D is therefore $\\widehat g_{e,j}=(M_e/Q_e)g_j$: its amplitude comes from the pulse measurement, and only its unfitted spatial shape comes from the native Jacobian. The raw native gain $Q_e$ is an independent numerical check. A localized pulse alone cannot determine every spatial component, so the native linearization is an explicit additional assumption, not hidden identification from one pulse. The covariance of the actual lognormal forcing is\n\n$$\nK_u(j,k;t,t+v)=\\mu^2\\left\\{\\exp\\left[\\sigma^2 C_x(j-k)C_t(t,t+v)\\right]-1\\right\\}. \\tag{3}\n$$\n\nThe periodic Gaussian spatial filter is discretely L2 normalized exactly as in the frozen campaign. Cx is its circular autocorrelation, including the zero spatial mode. Ct sums the four OU-knot covariances weighted by the two variance-normalized interpolation coefficients at each time. We integrate v from zero to 20τ at knot/12 spacing and average 24 midpoint knot phases. The omitted tail is bounded by the exponentially small latent covariance; numerical refinement is checked independently. No Gaussian small-sigma replacement is used at sigma=1. With the convention MSD=2Dt, the marginal edge diffusion is\n\n$$\nD_e=\\sum_{j,k}g_jg_k\\int_0^{\\infty}\\overline{K_u(j,k;t,t+v)}^{\\,t}\\,\\mathrm{d}v. \\tag{4}\n$$\n\nBoth zero-frequency edge kernels coincide because both boundaries share translation. Transient width and shape changes decay and are not assigned a nonzero asymptotic diffusion. The requested approximation and its covariance-complete counterpart are\n\n$$\nD_s^{\\mathrm{ind}}=\\frac{D_L+D_R}{4},\\qquad D_s=\\frac{D_L+D_R+2D_{LR}}{4}. \\tag{5}\n$$\n\nThe prescribed independence assumption therefore undercounts the shared neutral translation. Keeping covariance is a separate, untuned sensitivity; we retain the failing independence prediction in the main panel. The integrated-displacement summary is explicitly finite:\n\n$$\n\\mathcal A_{20\\tau}=\\frac{1}{a\\tau}\\int_0^{20\\tau}\\left[x_e^{\\mathrm{pulse}}(t)-x_e^{\\mathrm{control}}(t)\\right]\\,\\mathrm{d}t. \\tag{6}\n$$\n\nIt cannot be replaced by an integral to infinity when the pulse leaves a permanent position shift.

The kernel used for each marginal diffusivity is explicitly measurement-scaled:

$$
\\widehat g_{e,j}=\\frac{M_e}{Q_e}g_j,\\qquad M_e=\\frac{\\Delta x_e(20\\tau)}{a\\tau},\\quad Q_e=\\sum_jg_jh_{e,j}. \\tag{7}
$$

Equation (4) uses this measured-gain kernel in place of the raw native g. Neither $M_e$ nor the normalized spatial shape uses observed drift.\n\n## Numerical preparation and checks\n\nOriginal Figure 4 basal states are copied with their hashes. Each is allowed to settle at its exact condition mean for 16000 chemical units or 1600 Spring units. Eligibility requires a unique half-contrast region, contrast ≥0.05, max absolute native RHS ≤1e-6, and tail change ≤1e-4 of max(1, maximum state). Pulses use BDF rtol 2e-10, atol 2e-12, exact restart at withdrawal. Both edges receive independent half/nominal/double pulses; a matched unpulsed continuation removes any preparation drift. Current half-contrast crossings use periodic linear interpolation. Spring positions are material fractions times mean-input circumference; pure extension is excluded.\n\nNative sparse chemical Jacobians use centered numerical differentiation; Spring uses its complete analytic native Jacobian. Translation eigenvectors are selected by overlap with the native profile derivative; adjoint conservation modes are excluded. The spatial input derivative is measured directly from the unchanged native RHS. Every comparison table retains linearity, pulse validation, eigenvalue and stationary residuals. Amplitudes and horizons were fixed before diffusion comparisons.\n\n## Independence diagnostics and precision\n\nThe pre-existing float16 exemplars are the only saved profiles; at most 400 bins and 2000 frames each. Current half-contrast edges are identified only when exactly one region exists and current contrast is at least 0.1 of the original native basal contrast. Displacements are unwrapped within uninterrupted unique-region segments; no interpolation joins segments. For every lag we report left/right MSD, cross-product, centered increment correlation, and actual pair count. We retain the empirical cross-product in the variance identity, including systematic increments, rather than replacing it by the centered covariance. Filled supplement points require both MSDs greater than ten times a conservative squared float16 positional error bound, derived from half ULP of profile values and local crossing slope. Low-input short-lag results often fail this precision gate and remain open symbols. Two exemplar seeds describe covariance but do not establish a campaign confidence interval.\n\n## Undefined cases and scope\n\n'''
    methods+='\n'.join('- '+r['model']+f" C={r['size']:g}, f={r['fraction']:g}: "+r['reason'] for r in missing)+'\n\n'
    methods+='No missing execution is represented as censoring. Goryachev tracked drift can be very small and retains substantial disagreement. None of this equates tracked domain motion with whole-heading decay or with empirical migration persistence. Length represents WGD-associated enlargement; observed migration summaries and conditional OU fits remain distinct from these 1D polarity readouts.\n'
    methods+='\nThe maximum change in any terminal response over its final 0.5τ window is below 0.026% of the condition’s maximum terminal displacement; the transient is therefore settled within the observation window while a neutral positional offset remains. The covariance quadrature refinement changes the integrated kernel by 0.0026% (chemical) and 0.0007% (Spring). The declared interpolation differs from a continuous OU covariance and is retained exactly. Numerical checks are saved under `supporting/metadata/numerical_checks.json`.\n'
    (OUT/'METHODS_RESULTS.md').write_text(methods)
    (OUT/'README.md').write_text('''# Figure 3 round 5\n\n[Main PDF](figure3_round5.pdf) · [complete supplement](figure3_round5_supplement.pdf) · [legend](legend.md) · [methods/results](METHODS_RESULTS.md) · [reader](index.html).\n\n'''+section['narrative'][0]+'\n\n'+f"Deterministic campaign: `{rel(RUN)}`. SLURM array 21332598, 120 independent one-core tasks, 4 GB each, small QoS, 45-minute limit. Existing-profile edge analysis: SLURM 21332996 (cohort-corrected; first pass 21332981), one core / 4 GB. No new noise campaign. First failed API attempt is preserved under `figure3_round5_edge_pulses_20261023`, plus the incomplete preparation under `_prepare_failed`; scientific v2 inputs are frozen and valid. An unused solver-only retry in `figure3_round5_edge_numeric_retry_20261023` was canceled after the original strict run completed; no retry output supplies a figure datum.\n\n"+'''Rebuild using the existing scaling_sandbox Python:\n\n```sh\npython -B polarity1d/manuscript/code/mechanisms/figure3_round5_render_20261023.py\npython -B polarity1d/manuscript/code/mechanisms/figure3_round5_checks_20261023.py\n```\n\nEvery preceding supplement and main page is appended unchanged after the new supplement pages. Source hashes establish A/B preservation. No timed refresh. All data tables retain undefined statuses.\n''')
    generated=edge.now();body=f'<!doctype html><html lang="en"><meta charset="utf-8"><title>Figure 3 boundary response</title><style>body{{max-width:1300px;margin:2rem auto;padding:0 1rem;font:17px/1.5 system-ui;color:#182330}}img{{width:100%}}a{{color:#166aa2}}table{{border-collapse:collapse}}td,th{{padding:.4em;border:1px solid #ccc}}</style><h1>Figure 3: boundary response and predicted drift</h1><p>Data generated {generated}. Producers complete.</p><p><a href="figure3_round5.pdf">Main PDF</a> · <a href="figure3_round5_supplement.pdf">Complete supplement</a> · <a href="METHODS_RESULTS.md">Methods/results</a> · <a href="supporting/data/edge_prediction_comparison.csv">Prediction data</a> · <a href="supporting/data/edge_independence_summary.csv">Edge correlations</a></p><img src="figure3_round5.png" alt="Preserved turning assays, deterministic edge pulse and no-retuning drift comparison">'
    body+=''.join('<p>'+html.escape(t)+'</p>' for t in legend);body+='<h2>Numerical and covariance summary</h2>'+pd.DataFrame(stats).T.to_html(float_format=lambda x:f'{x:.4g}')+'</html>';(OUT/'index.html').write_text(body)
    write(META/'build_receipt.json',dict(generated_utc=generated,completed=True,native_conditions=len(frame),defined=int(frame.defined.sum()),undefined=len(missing),untuned=True,main_pages=len(PdfReader(OUT/'figure3_round5.pdf').pages),main_sha256=edge.sha(OUT/'figure3_round5.pdf'),legend_sha256=edge.sha(OUT/'legend.json'),source_manifest_sha256=edge.sha(RUN/'manifest.json')))

if __name__=='__main__':
    frame,pulses,curves,gains=collect();main_figure(frame,pulses,curves);supplement(frame,pulses,curves,gains);documents(frame);print(OUT)
