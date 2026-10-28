"""Rebuild R1 from portable compact CSVs; --extract refreshes from HPC records.

All scientific estimators are imported unchanged from round 3. No simulation,
fit or retuning is performed by the default portable rebuild.
"""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys
from datetime import datetime, timezone
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR',f'/tmp/ring-r1r2-mpl-{os.getuid()}')
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figure4_round3_estimators_20261022 as E
import ring_r1_campaign_20261027 as campaign
from ring_figure_revision1_20261028 import SIZE_COLORS, pair_rotation, rotate_display

ROOT=campaign.ROOT;CODE=campaign.CODE;BASE=campaign.BASE;RUN=campaign.RUN
DEV=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
R7=DEV/'revision_round7_20261026'
OUT=DEV/'ring_r1r2_20261027/figureR1'
DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
BLUE='#2878b5';SOURCES={}
def now():return datetime.now(timezone.utc).isoformat()
def rel(p):return str(Path(p).resolve().relative_to(ROOT))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(p):SOURCES[rel(p)]=sha(p);return p
def read(p):return pd.read_csv(source(p))
def write(p,v):campaign.write(p,v)
def setup():
    for name in ('data','metadata','vectors','viewers'):(OUT/'supporting'/name).mkdir(parents=True,exist_ok=True)
def markdown_table(frame):
    def cell(v):
        if isinstance(v,(float,np.floating)):
            return f'{v:.7g}' if np.isfinite(v) else str(v)
        return str(v)
    rows=['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    rows.extend('| '+' | '.join(cell(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None))
    return '\n'.join(rows)
def schematic_paths():
    """One unbranched cytosolic loop, with a resistor in each successive arc."""
    paths=[]
    def arc(name,start,stop,radius=1.,resistor=False):
        theta=np.linspace(start,stop,13 if resistor else 81)
        r=np.full(len(theta),radius)
        if resistor:
            r[2:-2]+=.065*np.where(np.arange(len(theta)-4)%2,1,-1)
        paths.extend(dict(path=name,order=i,x=rr*np.cos(t),y=rr*np.sin(t))
                     for i,(t,rr) in enumerate(zip(theta,r)))
    # Break the circle beneath each zigzag: no chord, rail or bypass branch.
    arc('upper_right',.20*np.pi,.37*np.pi)
    arc('domain_resistor',.37*np.pi,.63*np.pi,resistor=True)
    arc('upper_left',.63*np.pi,.80*np.pi)
    arc('back_left',.80*np.pi,1.35*np.pi)
    arc('back_resistor',1.35*np.pi,1.65*np.pi,resistor=True)
    arc('back_right',1.65*np.pi,2.20*np.pi)
    arc('domain',.20*np.pi,.80*np.pi,radius=1.17)
    return pd.DataFrame(paths)
def seed_summary(paths,identity,dense=False):
    moments=[];circ=[]
    for p in paths:
        with np.load(p) as z:
            t=z['time'].astype(float);P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float);cad=float(np.median(np.diff(t)))
            if dense:
                steps=np.unique(np.r_[1,np.arange(2,101,2),[120,160,200,250,320,400,500,640,800,1000,1280,1600,2000,2500,3000,3200]])
                lags=steps*cad;target=lags
            else:target,lags=E.lag_grid('wave_pinning',t[-1]-t[0],cad)
            moments.append(E.seed_moments(P,abs(P)>=.1,cad,lags));circ.append(float(z['circumference'].mean()))
    return E.summarize(moments,circ,identity,target,lags)
def summary_rows(s,**base):
    rows=[]
    for i,lag in enumerate(s['lags']):
        lo,hi=E.percentile(s['Cboot'][:,i]);klo,khi=E.percentile(s['kboot'][:,i],True)
        rows.append(dict(**base,lag=lag,C=s['C'][i],C_ci_lower=lo,C_ci_upper=hi,k=s['k'][i],k_ci_lower=klo,k_ci_upper=khi,
                         C_mean=s['circumference'],n_seeds=s['n_seeds'],coverage=s['coverage'][i],
                         eligible_pairs=s['eligible_pairs'][i],possible_pairs=s['possible_pairs'][i]))
    return rows
def extract():
    setup();source(RUN/'manifest.json');sources=[]
    for C in (20,60):
        p=source(R7/f'figure4/supporting/data/A_drift_C{C}.csv.gz')
        shutil.copy2(p,DATA/f'A_drift_C{C}.csv.gz')
    # C/D f=.01 are refreshed directly from the native and Campaign 2
    # compact seed archives by refresh_cd(), after this all-panel extraction.
    lags=np.linspace(0,24000,481)
    pd.DataFrame([dict(C_label=C,lag=l,Ds=.002,C=np.exp(-(2*np.pi/C)**2*.002*l)) for C in (20,60) for l in lags]).to_csv(DATA/'B_analytic_drift.csv',index=False)
    native=read(R7/'figure4/supporting/data/condition_lag_metrics.csv')
    native=native[(native.model=='wave_pinning')&np.isclose(native.f,.01)&native.C_label.isin([20,60])].copy();native['D']=10.
    acfs=native.to_dict('records')
    native_alpha=read(R7/'figure4/supporting/data/local_alpha.csv')
    native_alpha=native_alpha[(native_alpha.model=='wave_pinning')&native_alpha.f.isin([.003,.005,.01])].copy();native_alpha['D']=10.
    alphas=native_alpha.to_dict('records')
    fast_index=RUN/'results/compact_index.csv'
    fast=pd.read_csv(fast_index) if fast_index.exists() and fast_index.stat().st_size>1 else pd.DataFrame()
    fast_lags=[];checks=[]
    if len(fast):
        source(fast_index)
        for f in (.003,.005,.01):
            groups={}
            for C in (20,30,40,50,60):
                q=fast[(fast.D==1000)&np.isclose(fast.f,f)&np.isclose(fast.C,C)].sort_values('replicate')
                if len(q):
                    groups[float(C)]=seed_summary(q.compact,('R1fast',f,C));fast_lags.extend(summary_rows(groups[float(C)],D=1000.,f=f,C_label=C))
                    sources.extend(q[['slot','D','C','f','replicate','seed','compact','compact_sha256','native_compact','native_compact_sha256','paired_noise_verified']].to_dict('records'))
                if f==.01 and C in (20,60) and len(q):
                    dense=seed_summary(q.compact,('R1fast_dense',f,C),True);acfs.extend(summary_rows(dense,D=1000.,f=f,C_label=C))
                    q0=q[q.replicate==0]
                    if len(q0) and q0.iloc[0].exemplar:
                        p=Path(q0.iloc[0].exemplar);source(p)
                        with np.load(p) as z:
                            keep=np.linspace(0,len(z['time'])-1,min(800,len(z['time'])),dtype=int);t=z['time'][keep];u=z['primary'][keep];x=z['material_positions']
                        pd.DataFrame(dict(time=np.repeat(t,len(x)),relative_arc_position=np.tile(x,len(t)),activity=u.ravel())).to_csv(DATA/f'C_fast_C{C}.csv.gz',index=False)
            for C in (20,30,40,50,60):alphas.append(dict(model='wave_pinning',D=1000.,f=f,**E.local_alpha(groups,float(C),'wave_pinning')))
        for C in (20,60):
            q=fast[(fast.D==10)&np.isclose(fast.C,C)].sort_values('replicate')
            if len(q):
                a=seed_summary(q.native_compact,('R1solver',C));b=seed_summary(q.compact,('R1solver',C));j=E.analysis_lag(a)
                if j is not None:
                    ratio=b['k'][j]/a['k'][j];lo,hi=E.percentile(b['kboot'][:,j]/a['kboot'][:,j],True)
                    checks.append(dict(C=C,n_seeds=len(q),lag=a['lags'][j],native_k=a['k'][j],round7_k=b['k'][j],ratio=ratio,ratio_ci_lower=lo,ratio_ci_upper=hi))
                    sources.extend(q[['slot','D','C','f','replicate','seed','compact','compact_sha256','native_compact','native_compact_sha256','paired_noise_verified']].to_dict('records'))
    if not len(fast):
        for f in (.003,.005,.01):
            for C in (20,30,40,50,60):alphas.append(dict(model='wave_pinning',D=1000.,f=f,C_label=C,status='unavailable_campaign_pending'))
    pd.DataFrame(acfs).to_csv(DATA/'D_heading_acf.csv',index=False)
    guides=[]
    for D in (10.,1000.):
        q=pd.DataFrame(acfs);q=q[(q.D==D)&(q.C_label==20)]
        for r in q.itertuples():guides.append(dict(D=D,lag=r.lag,C=np.power(r.C,1/9) if r.C>=0 else np.nan,source_C20=r.C,status='defined' if r.C>=0 else 'undefined_negative_C20'))
    pd.DataFrame(guides).to_csv(DATA/'D_drift_guides.csv',index=False)
    pd.DataFrame(alphas).to_csv(DATA/'H_local_alpha.csv',index=False)
    pd.DataFrame(fast_lags,columns=list(fast_lags[0]) if fast_lags else ['D','f','C_label','lag','C','C_ci_lower','C_ci_upper','k','k_ci_lower','k_ci_upper','C_mean','n_seeds','coverage','eligible_pairs','possible_pairs']).to_csv(DATA/'H_fast_condition_lag_metrics.csv',index=False)
    pd.DataFrame(checks,columns=['C','n_seeds','lag','native_k','round7_k','ratio','ratio_ci_lower','ratio_ci_upper']).to_csv(DATA/'solver_check.csv',index=False)
    pulse=read(R7/'figure3/supporting/data/D_linearity_plot.csv');pulse=pulse[pulse['size'].isin([20,60])].copy()
    pulse['dose']=pulse.amplitude_over_reference*.1*50; pulse.to_csv(DATA/'E_edge_pulse.csv',index=False)
    fit=[];curves=[]
    for C,q in pulse.groupby('size'):
        # Fit only the declared weak, linear dose range (f<=.001); no model parameter is changed.
        w=q[q.amplitude_over_reference<=.001];gain=float(w.dose@w.permanent_displacement/(w.dose@w.dose));fit.append(dict(C=C,gain=gain,max_fit_fraction=.001,n_fit=len(w)))
        for dose in np.linspace(0,q.dose.max(),101):curves.append(dict(C=C,dose=dose,displacement=gain*dose))
    pd.DataFrame(fit).to_csv(DATA/'E_fitted_gain.csv',index=False);pd.DataFrame(curves).to_csv(DATA/'E_fitted_lines.csv',index=False)
    front=json.loads(source(R7/'figure3/supporting/metadata/isolated_front.json').read_text());shutil.copy2(R7/'figure3/supporting/metadata/isolated_front.json',META/'isolated_front.json')
    shutil.copy2(source(R7/'figure3/supporting/data/isolated_front.csv'),DATA/'G_isolated_front.csv')
    conduction=read(R7/'figure3/supporting/data/conduction_wave_pinning.csv');law=[]
    for r in conduction.itertuples():
        Cstar=2*r.D/(front['kappa']*front['delta_u']*r.theta*(1-r.theta));prediction=(front['lambda_push']/2)/(1+r.C/Cstar)
        law.append(dict(C=r.C,D=r.D,gain=r.gain,theta=r.theta,Cstar=Cstar,predicted_gain=prediction,relative_error=prediction/r.gain-1,defined=r.defined))
    law=pd.DataFrame(law);law.to_csv(DATA/'G_conduction_points.csv',index=False)
    line=[]
    for D,q in law.groupby('D'):
        q=q.sort_values('C')
        for C in np.linspace(20,60,161):
            theta=np.interp(C,q.C,q.theta);Cstar=2*D/(front['kappa']*front['delta_u']*theta*(1-theta))
            line.append(dict(C=C,D=D,theta=theta,Cstar=Cstar,gain=(front['lambda_push']/2)/(1+C/Cstar)))
    line=pd.DataFrame(line);line.to_csv(DATA/'G_conduction_law.csv',index=False)
    prediction=line[line.D==10].copy();prediction['alpha']=2+2*(prediction.C/prediction.Cstar)/(1+prediction.C/prediction.Cstar)
    prediction[['C','Cstar','alpha']].to_csv(DATA/'H_native_prediction.csv',index=False)
    schematic_paths().to_csv(DATA/'F_schematic_paths.csv',index=False)
    pd.DataFrame(sources).to_csv(META/'campaign_sources.csv',index=False)
    manifest=json.loads((RUN/'manifest.json').read_text());status=json.loads((RUN/'status.json').read_text()) if (RUN/'status.json').exists() else dict(complete=False,completed=0,expected=len(manifest['jobs']))
    write(META/'campaign_manifest.json',manifest);shutil.copy2(RUN/'submission.json',META/'campaign_submission.json')
    write(META/'extraction_provenance.json',dict(generated_utc=now(),producer_status=status,input_sha256=SOURCES,
        code_sha256=sha(Path(__file__)),estimators=str(E.__file__),estimators_sha256=sha(E.__file__),bootstrap_resamples=2000,
        prediction='Closed-form isolated-front law, no gain fitting; row theta from exact ring states; theta linearly interpolated between known C',
        pulse_slope='Through-origin least squares over f<=.001, edge dose=50*0.1*f; nonlinear higher doses retained',
        missing='Uncompleted campaign conditions marked unavailable, never imputed',new_profiles_on_disk=False))

def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,'axes.titlesize':7,'xtick.labelsize':6,'ytick.labelsize':6,
        'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,'lines.linewidth':1,'lines.markersize':3,
        'xtick.major.width':.55,'ytick.major.width':.55,'xtick.major.size':2,'ytick.major.size':2,
        'pdf.fonttype':42,'svg.fonttype':'none','text.color':'black','axes.titlecolor':'black'})
def letter(fig,spec,label):
    b=spec.get_position(fig);fig.text(b.x0-.040,b.y1+.012,label,weight='bold',size=10,color='black')
def facet(ax,label):ax.set_title(label,pad=4,color='black')
def missing(ax,text='Campaign pending'):
    ax.text(.5,.5,text,transform=ax.transAxes,ha='center',va='center',color='.4',size=7)
def refresh_cd():
    """Read existing f=.01 records only; preserve H and the solver check."""
    setup();index=pd.read_csv(RUN/'results/compact_index.csv')
    native_index=pd.read_csv(BASE/'results/compact_index.csv');acfs=[];inputs={}
    for C in (20,60):
        q=index[(index.D==1000)&np.isclose(index.f,.01)&(index.C==C)].sort_values('replicate')
        assert len(q)==384
        acfs.extend(summary_rows(seed_summary(q.compact,('R1fast_dense',.01,C),True),D=1000.,f=.01,C_label=C))
        nq=native_index[(native_index.model=='wave_pinning')&np.isclose(native_index.fraction,.01)&(native_index['size']==C)].sort_values('replicate')
        assert len(nq)==384
        acfs.extend(summary_rows(seed_summary(nq.compact,('R1native_dense',.01,C),True),D=10.,f=.01,C_label=C))
        for kind,path in (('native',Path(nq.iloc[0].exemplar)),('fast',Path(q.iloc[0].exemplar))):
            inputs[rel(path)]=sha(path)
            with np.load(path) as z:
                t=z['time'];u=z['primary'];x=z['material_positions']
            assert np.isclose(t[0],16000) and np.isclose(t[-1],32000)
            pd.DataFrame(dict(time=np.repeat(t,len(x)),relative_arc_position=np.tile(x,len(t)),activity=u.ravel())).to_csv(DATA/f'C_{kind}_C{C}.csv.gz',index=False)
    d=pd.DataFrame(acfs);d.to_csv(DATA/'D_heading_acf.csv',index=False)
    guides=[]
    for r in d[d.C_label==20].itertuples():
        guides.append(dict(D=r.D,lag=r.lag,C=r.C**(1/9) if r.C>=0 else np.nan,source_C20=r.C,status='defined' if r.C>=0 else 'undefined_negative_C20'))
    pd.DataFrame(guides).to_csv(DATA/'D_drift_guides.csv',index=False)
    provenance=json.loads((META/'extraction_provenance.json').read_text())
    provenance.update(generated_utc=now(),revision1_cd=dict(f=.01,post_burn_window=[16000,32000],n_seeds_per_condition=384,input_sha256=inputs,code_sha256=sha(__file__),scope='C/D only; H and solver-check records preserved'))
    write(META/'extraction_provenance.json',provenance)

def matrix(path):
    d=pd.read_csv(path);t=np.sort(d.time.unique());return t,d.activity.to_numpy().reshape(len(t),-1)

def kymo(ax,p,C,kind,vmax=1.,ylabel=False,rotation=0.):
    if p.exists():
        t,u=matrix(p);u=rotate_display(u,rotation)
        im=ax.imshow(u,origin='lower',extent=[0,C,float(t[0]),float(t[-1])],aspect='auto',cmap='viridis',vmin=0,vmax=vmax,interpolation='nearest',rasterized=True)
        # Highlight the periodic display seam when the dominant path crosses it.
        if np.any(abs(np.diff(np.argmax(u,axis=1)))>u.shape[1]/2):
            for seam in (0,C):ax.axvline(seam,color='white',ls=':',lw=.5,zorder=4,clip_on=False)
    else:missing(ax);im=None;ax.set_ylim(16000,32000)
    ax.set_xlim(0,C);ax.set_xticks([0,C]);ax.set_xlabel('Arc position');facet(ax,f'C{C}')
    if ylabel:ax.set_ylabel('Time')
    else:ax.tick_params(labelleft=False)
    return im
def render(export=True):
    setup();style();fig=plt.figure(figsize=(180/25.4,262/25.4))
    rows=fig.add_gridspec(5,1,left=.09,right=.96,bottom=.055,top=.97,height_ratios=[1.05,1.75,1.05,1.18,1.20],hspace=.66)
    r0=rows[0].subgridspec(1,2,width_ratios=[1,1],wspace=.40)
    ga=r0[0].subgridspec(1,2,width_ratios=[1,3],wspace=.22);letter(fig,r0[0],'A')
    rotation=pair_rotation([matrix(DATA/f'A_drift_C{C}.csv.gz')[1] for C in (20,60)])
    for i,C in enumerate((20,60)):
        ax=fig.add_subplot(ga[i]);im=kymo(ax,DATA/f'A_drift_C{C}.csv.gz',C,'drift',ylabel=i==0,rotation=rotation);ax.set_yticks([0,12000,24000])
    # Shared activity scale is also the null's exact peak amplitude.
    cb=fig.colorbar(im,ax=fig.axes[:2],fraction=.025,pad=.035,ticks=[0,1]);cb.set_label('Activity')
    ax=fig.add_subplot(r0[1]);letter(fig,r0[1],'B');d=pd.read_csv(DATA/'B_analytic_drift.csv')
    for C,col in ((20,SIZE_COLORS[20]),(60,SIZE_COLORS[60])):
        q=d[d.C_label==C];ax.plot(q.lag,q.C,color=col,label=f'C{C}')
    ax.set(xlabel='Lag',ylabel='Heading correlation',ylim=(-.02,1.03),xlim=(0,24000));ax.set_xticks([0,12000,24000]);ax.legend(frameon=False,fontsize=6,title='Size',title_fontsize=6)
    r1=rows[1].subgridspec(1,2,width_ratios=[1.14,.86],wspace=.40);letter(fig,r1[0],'C');letter(fig,r1[1],'D')
    gc=r1[0].subgridspec(2,1,hspace=.55);kaxes=[]
    wave_vmax=max(matrix(DATA/f'C_{kind}_C{C}.csv.gz')[1].max() for kind in ('native','fast') for C in (20,60))
    for j,D in enumerate((10,1000)):
        kind='native' if D==10 else 'fast'
        rotation=pair_rotation([matrix(DATA/f'C_{kind}_C{C}.csv.gz')[1] for C in (20,60)])
        sub=gc[j].subgridspec(1,2,width_ratios=[1,3],wspace=.22)
        for i,C in enumerate((20,60)):
            ax=fig.add_subplot(sub[i]);kaxes.append(ax);im2=kymo(ax,DATA/f"C_{'native' if D==10 else 'fast'}_C{C}.csv.gz",C,'wave',vmax=wave_vmax,ylabel=i==0,rotation=rotation)
            ax.set_yticks([16000,24000,32000])
            if j==0:ax.set_xlabel('');ax.tick_params(labelbottom=False)
        pos=gc[j].get_position(fig);fig.text(pos.x0,pos.y1+.022,f'$D_v$ = {D}'+(' (native)' if D==10 else ''),size=7)
    if any(a.images for a in kaxes):
        im2=next(a.images[0] for a in kaxes if a.images);cb=fig.colorbar(im2,ax=kaxes,fraction=.025,pad=.035,ticks=[0,wave_vmax]);cb.set_label('Activity')
    gd=r1[1].subgridspec(2,1,hspace=.55);d=pd.read_csv(DATA/'D_heading_acf.csv');g=pd.read_csv(DATA/'D_drift_guides.csv')
    for j,D in enumerate((10,1000)):
        ax=fig.add_subplot(gd[j]);facet(ax,f'$D_v$ = {D}')
        for C,col in ((20,SIZE_COLORS[20]),(60,SIZE_COLORS[60])):
            q=d[(d.D==D)&(d.C_label==C)].sort_values('lag');ax.plot(q.lag,q.C,color=col,label=f'C{C}')
            if len(q):ax.fill_between(q.lag,q.C_ci_lower,q.C_ci_upper,color=col,alpha=.13,linewidth=0)
        q=g[g.D==D];ax.plot(q.lag,q.C,'--',color='.5',label='Fixed physical drift')
        if not len(d[d.D==D]):missing(ax)
        ax.set(ylim=(-.05,1.01),xlim=(0,16000),ylabel='Heading correlation');ax.set_xticks([0,8000,16000]);ax.set_xlabel('Lag' if j else '')
        if j==0:
            handles,labels=ax.get_legend_handles_labels();handles.append(Line2D([0],[0],lw=4,color='.4',alpha=.2));labels.append('95% seed CI')
            ax.legend(handles,labels,frameon=False,fontsize=5.8,loc='lower left',labelspacing=.15)
    r2=rows[2].subgridspec(1,2,wspace=.40);letter(fig,r2[0],'E');letter(fig,r2[1],'F')
    ax=fig.add_subplot(r2[0]);d=pd.read_csv(DATA/'E_edge_pulse.csv');gains=pd.read_csv(DATA/'E_fitted_gain.csv');fit_lines=[]
    for C,col in ((20,SIZE_COLORS[20]),(60,SIZE_COLORS[60])):
        q=d[d['size']==C];ax.plot(q.dose,q.permanent_displacement,'o',color=col,ms=2.4,label=f'C{C}')
        doses=np.geomspace(d.dose.min(),d.dose.max(),200);gain=float(gains[gains.C==C].gain.iloc[0]);ax.plot(doses,gain*doses,color=col,lw=.9);fit_lines.extend(dict(C=C,dose=x,displacement=gain*x) for x in doses)
    if export:pd.DataFrame(fit_lines).to_csv(DATA/'E_fitted_lines.csv',index=False)
    weak=d[d.amplitude_over_reference<=.001];ax.axvspan(d.dose.min(),weak.dose.max(),color='.5',alpha=.10,lw=0,zorder=0)
    ax.set_xscale('log');ax.set_yscale('log')
    ax.set(xlabel='Edge-pulse dose',ylabel='Permanent displacement',xlim=(d.dose.min()*.8,d.dose.max()*1.2));ax.legend(frameon=False,fontsize=6,title='Lines: fitted gain (slope 1)',title_fontsize=6)
    # Schematic geometry is generated on portable rebuild as well as extraction.
    if export:schematic_paths().to_csv(DATA/'F_schematic_paths.csv',index=False)
    ax=fig.add_subplot(r2[1]);d=pd.read_csv(DATA/'F_schematic_paths.csv')
    for p,q in d.groupby('path',sort=False):
        ax.plot(q.x,q.y,color=BLUE if p=='domain' else '.45',lw=4 if p=='domain' else .9,zorder=2 if p=='domain' else 1)
    edge=(np.cos(.20*np.pi),np.sin(.20*np.pi))
    ax.annotate('Pulse',xy=edge,xytext=(1.40,1.03),arrowprops=dict(arrowstyle='->',color='black',lw=.8),ha='center',fontsize=7)
    ax.plot([edge[0],-edge[0]],[edge[1],edge[1]],'o',color=BLUE,ms=3)
    ax.text(0,1.31,'Domain',ha='center',size=7,color='black')
    ax.text(0,.30,'Across domain',ha='center',size=6.3)
    ax.text(0,-.50,'Around back',ha='center',size=6.3)
    for start,stop in ((.60*np.pi,.44*np.pi),(1.70*np.pi,1.86*np.pi)):
        ax.annotate('',xy=(np.cos(stop),np.sin(stop)),xytext=(np.cos(start),np.sin(start)),arrowprops=dict(arrowstyle='->',color='black',lw=.8,connectionstyle='arc3,rad=-.15'),zorder=5)
    ax.set(xlim=(-1.22,1.84),ylim=(-1.13,1.52),aspect='equal');ax.axis('off')
    ax=fig.add_subplot(rows[3]);letter(fig,rows[3],'G');d=pd.read_csv(DATA/'G_conduction_points.csv');l=pd.read_csv(DATA/'G_conduction_law.csv')
    Ds=sorted(d.D.unique())
    for i,D in enumerate(Ds):
        col='.25';q=d[d.D==D].sort_values('C');ax.plot(q.C,q.gain,'o',color=col,ms=3,label=f'{D:g}')
        q=l[l.D==D].sort_values('C');ax.plot(q.C,q.gain,color=col,lw=.9)
        ax.text(60.7,q.gain.iloc[-1],f'{D:g}',color=col,fontsize=6,va='center')
    ax.set(xlabel='Circumference, C',ylabel='Edge gain, G',xlim=(19,65),ylim=(.45,2.9));ax.set_xticks([20,30,40,50,60])
    ax.legend(handles=[Line2D([],[],color='.25',marker='o',lw=0,label='Exact ring'),Line2D([],[],color='.25',label='Law')],frameon=False,fontsize=6,ncol=2,loc='lower left');ax.text(1.,1.02,'$D_v$',transform=ax.transAxes,ha='right',fontsize=7)
    gh=rows[4].subgridspec(1,3,wspace=.28);letter(fig,rows[4],'H');d=pd.read_csv(DATA/'H_local_alpha.csv');pred=pd.read_csv(DATA/'H_native_prediction.csv')
    for j,f in enumerate((.003,.005,.01)):
        ax=fig.add_subplot(gh[j]);facet(ax,f'f = {f:g}');ax.axhline(2,color='.70',ls=':',lw=.8,label='Geometric drift: 2')
        ax.plot(pred.C,pred.alpha,color='.5',lw=.9,label='Native prediction')
        for D,marker in ((10,'o'),(1000,'s')):
            q=d[(d.D==D)&np.isclose(d.f,f)&(d.status=='resolved')&d.C_label.isin([30,40,50])].sort_values('C_label')
            if len(q):
                for r in q.itertuples():
                    ax.errorbar(r.C_label+(-.65 if D==10 else .65),r.alpha,yerr=[[max(0,r.alpha-r.alpha_ci_lower)],[max(0,r.alpha_ci_upper-r.alpha)]],fmt=marker,ms=3,color='black',mfc='white' if D==1000 else 'black',elinewidth=.65,capsize=1.5,label=f'$D_v$ = {D} (95% CI)' if r.Index==q.index[0] else None)
            if not len(q) and D==1000:ax.text(.5,.08,'Fast cytosol pending',transform=ax.transAxes,ha='center',color='.4',fontsize=5.5)
        ax.set(xlabel='Circumference, C',xlim=(17,63),ylim=(0,5.5),xticks=[20,40,60]);ax.set_ylabel('Local exponent, α' if j==0 else '')
        if j:ax.tick_params(labelleft=False)
        if j==0:ax.legend(frameon=False,fontsize=5.8,loc='upper left',labelspacing=.1)
    if not export:return fig
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'figureR1_main.{ext}',dpi=300,facecolor='white')
    plt.close(fig);documents()

def documents():
    provenance=json.loads((META/'extraction_provenance.json').read_text());state=provenance['producer_status'];complete=state.get('complete',False)
    law=pd.read_csv(DATA/'G_conduction_points.csv');alpha=pd.read_csv(DATA/'H_local_alpha.csv');solver=pd.read_csv(DATA/'solver_check.csv');gain=pd.read_csv(DATA/'E_fitted_gain.csv')
    native=alpha[(alpha.D==10)&(alpha.status=='resolved')&alpha.C_label.isin([30,40,50])];fast=alpha[(alpha.D==1000)&(alpha.status=='resolved')&alpha.C_label.isin([30,40,50])]
    def rng(x):return f'{x.min():.3g}–{x.max():.3g}' if len(x) else 'unavailable'
    g20=float(law[(law.D==1000)&(law.C==20)].gain.iloc[0]);g60=float(law[(law.D==1000)&(law.C==60)].gain.iloc[0])
    afcf=pd.read_csv(DATA/'D_heading_acf.csv');guides=pd.read_csv(DATA/'D_drift_guides.csv');q=afcf[(afcf.D==1000)&(afcf.C_label==60)]
    discrepancy=None
    if len(q):
        m=q.merge(guides[guides.D==1000],on='lag',suffixes=('_observed','_guide'));discrepancy=float(np.max(abs(m.C_observed-m.C_guide)))
    result_C='Fast-cytosol sample paths are pending.' if not (DATA/'C_fast_C60.csv.gz').exists() else 'Both sizes retain one moving activity domain; the paired fast-cytosol paths are shown without seed selection.'
    result_D='The fast-cytosol condition is pending.' if discrepancy is None else f'The largest absolute C60 departure from the fast-cytosol fixed-drift guide is {discrepancy:.3g} over the observed lag grid.'
    result_H=f'Native point estimates span {rng(native.alpha)}.'+(f' Fast-cytosol estimates span {rng(fast.alpha)}.' if len(fast) else ' Fast-cytosol estimates are pending.')
    legend=[
        '# Figure R1. Cytosolic transport suppresses Wave-pinning domain drift beyond the geometric size effect',
        r'**A.** Drift-null activity kymographs at circumference $C=20$ and $60$, with panel width proportional to circumference, time upward and a shared linear viridis activity scale with its colourbar. The fixed physical diffusivity is $D_s=0.002$. A larger ring converts the same physical wandering into smaller heading changes. For each size pair, one common angular rotation places the initial dominant domain near mid-panel; the identical angular shift is used for both sizes at all times. Thin dotted seam lines mark any remaining wrap. All quantities in this figure use native units.',
        r'**B.** Closed-form heading autocorrelation $C(\ell)=\exp[-(2\pi/C)^2D_s\ell]$ for A; orange (#E69F00) is C20 and blue (#0072B2) C60. The heading-decorrelation rate is $k=-\log C(\ell)/\ell$. The geometric drift null therefore gives a ninefold smaller rate at C60.',
        r'**C.** Wave-pinning activity at input fraction $f=0.01$ (mean forcing divided by native reference input 0.1), with inactive diffusivity $D_v=10$ (native) or 1000. Each size uses native replicate 0 and the identical frozen Gaussian noise realization across diffusivities; the full post-burn window 16000–32000 is shown, time and linear viridis activity scales are shared across both diffusivities and sizes, and width is proportional to circumference. '+result_C,
        r'**D.** Whole-activity heading autocorrelation for C over the full lag horizon 0–16000, orange C20 and blue C60, with whole-seed bootstrap 95% bands. Heading is $\arg P$, with $P=M_1/M_0$; eligible endpoints have $|P|\geq0.1$. Dashed grey curves are the size-independent physical-drift guide $C_{20}(\ell)^{(20/60)^2}$. '+result_D,
        r'**E.** Permanent shared domain displacement against edge-pulse dose (peak amplitude times duration 50), with orange C20 and blue C60. Log–log axes retain the full dose grid. Points are retained deterministic pulses; each slope-1 line is its through-origin fitted weak-dose gain (input fractions at most 0.001); light shading marks that fit range, defining edge gain $G$. The C60 slope is smaller, and higher doses depart below the linear extrapolation.',
        r'**F.** Ring schematic with one active domain, an edge pulse and the two cytosolic paths between its edges, across the domain and around the back, shown as resistor segments in series around the closed loop. Arrows on both arcs show material moving from the far edge to the pushed edge. The schematic identifies the two transport distances that enter the conduction law.',
        r'**G.** Edge gain against circumference for seven inactive diffusivities. Every series is dark grey and directly labelled at its right end by $D_v$; the small key defines points as exact-ring adjoint calculations and lines as $G=(\lambda/2)/(1+C/C^*)$, where $C^*=2D_v/[\kappa\Delta u\,\theta(1-\theta)]$. Isolated-front values of $\kappa$, $\lambda$ and $\Delta u$ are fixed and $\theta$ comes from ring states; nothing is fitted. Native gain falls with size while the $D_v=1000$ C60/C20 gain ratio is '+f'{g60/g20:.4f}.',
        r'**H.** Local size exponent $\alpha=-d\log k/d\log C$ from the unchanged three-size round-3 estimator, at the three weak inputs available at all five sizes. Estimates are centred at C30, C40 and C50; each uses its three adjacent sizes, so all five sizes contribute. Black filled circles denote native cytosol ($D_v=10$) and black open squares fast cytosol ($D_v=1000$), with 95% whole-seed bootstrap intervals; the solid grey line is $\alpha=2+2(C/C^*)/(1+C/C^*)$ at native diffusivity and the dotted guide is 2. Inputs are retained in separate shared-scale facets. '+result_H]
    (META/'main_legend.md').write_text('\n\n'.join(legend)+'\n')
    unexpected=[f'The full seven-diffusivity law does not satisfy a universal 10% agreement check: at Dᵥ=2.5 the relative error is {rng(law[law.D==2.5].relative_error*100)}%, and at Dᵥ=5 it is {rng(law[law.D==5].relative_error*100)}%. The native Dᵥ=10 error is {rng(law[law.D==10].relative_error*100)}%; the Dᵥ=1000 error is {rng(law[law.D==1000].relative_error*100)}%. The per-row stationary θ and isolated-front constants were checked; no parameters were changed.']
    for r in solver.itertuples():
        if r.ratio_ci_upper<1 or r.ratio_ci_lower>1:
            unexpected.append(f'The C{r.C:g} native solver check has k(round 7)/k(round 3)={r.ratio:.6f} (95% paired seed CI {r.ratio_ci_lower:.6f}–{r.ratio_ci_upper:.6f}), a {(r.ratio-1)*100:.3g}% difference whose interval excludes exact equality. The seed, Gaussian field, observation window and estimator are identical; the fixed IMEX time discretization differs from native BDF. This small solver difference is reported without adjustment.')
    if not complete:unexpected.append(f"Campaign 2 is incomplete at extraction: {state.get('completed',0)}/{state.get('expected',5776)} records completed. Missing curves and estimates are labelled pending; no missing seed is imputed. SLURM submission and frozen manifest are beside this readme.")
    if discrepancy is not None and discrepancy>.1:unexpected.append(f'The fast C60 curve differs from the fixed-physical-drift guide by as much as {discrepancy:.3g}, exceeding a 0.1 correlation-unit descriptive check; the full measured curve is retained.')
    if len(fast) and (fast.alpha.min()<1.5 or fast.alpha.max()>2.5):unexpected.append('At least one fast-cytosol local point estimate lies outside 1.5–2.5; see the complete estimates and intervals below. No value was adjusted toward 2.')
    unavailable=alpha[alpha.status!='resolved']
    if len(unavailable):unexpected.append(f'{len(unavailable)} local-α estimates are unavailable. Their condition and estimator status are retained in H_local_alpha.csv; missing estimates are not censoring bounds and receive no numeric marker.')
    text=['# Unexpected','',*['- '+s for s in unexpected],'','## Key numbers per panel','',
        'A–B: Dₛ = 0.002, C20/C60, T = 24000, sample-path cadence 20. Analytic k(C20) = '+f'{(2*np.pi/20)**2*.002:.8g}, k(C60) = {(2*np.pi/60)**2*.002:.8g}; ratio = 9. No null-ensemble simulation was added.',
        '',f'C–D: f = 0.01, post-burn window 16000–32000, observation cadence 5, 384 native seeds per condition. {result_D}',
        '', 'E: weak-dose slopes '+', '.join(f'C{r.C:g}: G={r.gain:.8g}' for r in gain.itertuples())+'. Higher-dose observations remain plotted.',
        '', 'F: qualitative transport-path schematic; no inferred measurement.',
        '',f'G: exact-ring gain C60/C20 at Dᵥ=1000 = {g60/g20:.8g}; '+', '.join(f'Dᵥ={D:g}: law relative error [{q.relative_error.min():.5g}, {q.relative_error.max():.5g}]' for D,q in law.groupby('D'))+'.',
        '', 'H: '+result_H,'',markdown_table(alpha[alpha.C_label.isin([30,40,50])][['D','f','C_label','status','alpha','alpha_ci_lower','alpha_ci_upper','lag','n_seeds']]),
        '', '## Methods','',
        'No model retuning. Campaign 2 uses the unchanged round7_noise_20261026.integrate kernel (IMEX Euler dt 0.05, periodic native dx 0.05). Fast runs start from each exact native basal state and use the native per-condition seed, horizon, burn-in, NoiseInput parameters and Gaussian knots. Each completed worker verifies the exact latent SHA256 against its native noise audit before integration. Native histories exist only in RAM and are discarded after the unchanged round-3 region tracker, compact descriptor extraction and exemplar reduction. Replicates 0 and 1 alone retain ≤2000 frames and ≤400 bins, as float16 illustration.',
        '', 'All autocorrelations use the unchanged round-3 endpoint cosine estimator and |P|≥0.1 endpoint rule. H uses the unchanged lag grid, fixed adjacent-size neighborhoods and shared-lag selector nearest correlation 0.5, OLS of log rate against log actual circumference, and 2000 whole-seed bootstrap resamples. H plots only three-size interior neighborhoods (C20/30/40, C30/40/50 and C40/50/60); the two-size endpoint estimates remain in the extraction table but are excluded from the requested three-size plot. H keeps f separate, avoiding a pooling assumption. Bands in D are pointwise 95% whole-seed intervals. Dashed D guides and B/H predictions are closed form, with no simulated null ensembles. A reuses the original fixed drift sample paths.',
        '', 'The deterministic conduction calculation uses isolated-front κ=1.281489794, λ=5.304509755, Δu=0.960849846. Each plotted exact ring supplies its own θ; the prediction line interpolates θ between retained circumference states, without fitting gains. E fits only the declared weak-dose range through the origin, keeping every larger-dose point.',
        '', '### Solver check','', 'Native Dᵥ=10 is rerun with the same round-7 worker for original replicates 0–7 at C20/C60, f=0.005. k is compared at the round-3 cohort analysis lag, using paired whole-seed bootstrap weights for the k ratio.',
        '', markdown_table(solver) if len(solver) else 'Pending native solver-check records; no solver-equivalence claim is made.',
        '', '### Rebuild and provenance','',
        f'Data generation time: **{provenance["generated_utc"]}**. Producer complete: **{complete}**. Portable rebuild: `'+rel(CODE/'rebuild_figureR1_20261027.py')+'` with no flags. To refresh compact data after campaign completion, run `ring_r1_campaign_20261027.py collect`, then `rebuild_figureR1_20261027.py --extract` under SLURM. Rebuild writes PNG, PDF and SVG at 180 mm width. Extraction hashes, exact source tables, frozen campaign design and submission receipts are beside this file.',
        '', 'Biological scope: size denotes WGD-associated enlargement. These simulated activity-heading and edge-motion readouts remain distinct from measured migration summaries and conditional OU velocity-persistence estimates. No migration-velocity law is assumed.']
    (META/'scientific_readme.md').write_text('\n'.join(text)+'\n')
    code_commit=subprocess.run(['git','log','-1','--format=%H','--',rel(Path(__file__))],cwd=ROOT,text=True,capture_output=True).stdout.strip() or 'not yet committed; integration_commit.json records final code commit'
    outputs={rel(p):sha(p) for p in sorted(DATA.glob('*')) if p.is_file()}
    write(META/'provenance.json',dict(generated_utc=now(),data_generated_utc=provenance['generated_utc'],code_commit=code_commit,
        code_path=rel(Path(__file__)),code_sha256=sha(Path(__file__)),display_conventions_sha256=sha(CODE/'ring_figure_revision1_20261028.py'),campaign_manifest=rel(META/'campaign_manifest.json'),
        campaign_manifest_sha256=sha(META/'campaign_manifest.json'),output_sha256=outputs,width_mm=180,height_mm=262,
        producer_complete=complete,panels='A B C D E F G H',scope='New R1 only; native Figures 1–6 unchanged'))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--extract',action='store_true');p.add_argument('--refresh-cd',action='store_true');a=p.parse_args()
    if a.extract:
        extract();refresh_cd()
    if a.refresh_cd:refresh_cd()
    render()
