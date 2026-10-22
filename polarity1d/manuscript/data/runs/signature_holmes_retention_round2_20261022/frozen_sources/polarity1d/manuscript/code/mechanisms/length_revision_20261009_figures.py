"""Review revision: figure-only graphics, attached manifests and native-array checks.

Run after length_revision_20261009_trials.py; no integrations in this renderer.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import LogLocator, NullFormatter, MaxNLocator, FuncFormatter

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core'),str(ROOT/'code/core/src')]
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning
from polarity1d.solvers import laplacian_no_flux
import figure_builder as fb
fb.CONDITIONS=['fixed_physical','fixed_relative','fixed_slope']

RUN=ROOT/'data/runs/length_revision_20261009'
FIG=ROOT/'figures/length_mechanisms_20261008/review_20261009_quality'
DATA=ROOT/'data/runs/length_mechanisms_20261008/review_20261009_quality'
POOL=ROOT/'data/runs/inactive_pool_figure_exports_20261008'
SINK=ROOT/'data/runs/inactive_sink_loss_20261008'
MODELS={'goryachev':('Goryachev',Goryachev(),'#8B3F82'), 'wave_pinning':('Wave-pinning',WavePinning(),'#007D73')}
LC={9:'#0072B2',21:'#8B3F82',39:'#D55E00'}
sources={}; measurements={}; evidence={}; checks={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':8,
    'pdf.fonttype':42,'svg.fonttype':'none'})


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def register(p): sources[str(p.relative_to(ROOT))]=sha(p)
def read(p): register(p);return json.loads(p.read_text())
def load(p):
    register(p)
    with np.load(p) as z:return {k:z[k].copy() for k in z.files}
def rows(p): register(p);return list(csv.DictReader(p.open()))
def csvwrite(p, data):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def panel(ax, letter, title=None):
    ax.text(-.13,1.04,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans')
    if title:ax.set_title(title,pad=10)
    ax.grid(alpha=.12,lw=.6)
def save(fig,stem,caption,methods):
    outputs={}
    for ext in ['png','pdf','svg']:
        p=FIG/f'{stem}.{ext}';fig.savefig(p,dpi=220,facecolor='white',bbox_inches='tight',pad_inches=.16);outputs[p.name]=sha(p)
    plt.close(fig)
    (FIG/f'{stem}_caption.txt').write_text(caption+'\n')
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',
        description=caption,legend=caption,methods=methods,sha256=outputs,
        source_sha256=dict(sources),renderer_sha256=sha(Path(__file__)),
        reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__)),
                      '/home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(ROOT/'code/mechanisms/length_mechanisms_report_20261008.py')])
    (FIG/f'{stem}.manifest.json').write_text(json.dumps(record,indent=2)+'\n')


def native_f(model,u,v):
    u,v=np.broadcast_arrays(u,v)
    return model.reaction(0.,np.zeros(u.size),np.stack([u.ravel(),v.ravel()]),np.zeros(u.size))[0].reshape(u.shape)


def audit_trials():
    m=read(RUN/'manifest.json')
    for path,digest in m['sources'].items():assert sha(ROOT/path)==digest,path
    audits=[];trials={}
    for slot,j in enumerate(m['jobs']):
        c=j['case'];p=RUN/f'workers/{slot:02d}.npz';meta=read(p.with_suffix('.json'));z=load(p)
        assert meta['job']==j and sha(p)==meta['archive_sha256']
        cp=load(ROOT/c['checkpoint']);assert sha(ROOT/c['checkpoint'])==c['checkpoint_sha256']
        np.testing.assert_array_equal(z['originalstate'],cp['state'])
        # BDF dense-output evaluation at t=0 can round the supplied state.
        np.testing.assert_allclose(z['pulse_active'][0],cp['state'][0],rtol=2e-13,atol=1e-15)
        np.testing.assert_allclose(z['pulse_field_inactive'][0],cp['state'][1],rtol=2e-13,atol=1e-15)
        errors=[]
        for phase in ['pulse','recovery']:
            if phase+'_active' not in z:continue
            u,v=z[phase+'_active'],z[phase+'_field_inactive'];x=z[phase+'_x'];L=c['length']
            assert u.shape==v.shape and np.isfinite(u).all() and min(u.min(),v.min())>-1e-8
            error=float(np.max(abs((u+v).sum(axis=1)/cp['state'].sum()-1.)))
            assert error<1e-6,error
            errors.append(error)
            np.testing.assert_allclose(u@(2*x/L-1)/u.sum(axis=1),z[phase+'_orientation'],atol=2e-13)
        result=dict(slot=slot,model=c['model'],length=c['length'],checkpoint=True,
                    mass_error=max(errors),native_polarity_recomputed=True)
        if j['recovery']:
            np.testing.assert_allclose(z['recovery_active'][0],z['pulse_active'][-1],rtol=2e-13,atol=1e-15)
            np.testing.assert_allclose(z['recovery_field_inactive'][0],z['pulse_field_inactive'][-1],rtol=2e-13,atol=1e-15)
            tail=z['recovery_time']>=9500
            target=cp['state'] if j['duration']==20 else cp['state'][:,::-1]
            final=np.stack([z['recovery_active'][-1],z['recovery_field_inactive'][-1]])
            scales=np.array([c['fixed_scales'][k] for k in ['active','inactive']])[:,None]
            match=float(np.max(abs(final-target)/scales));assert match<1e-3,match
            drift=max(float(np.ptp(z['recovery_'+key][tail],axis=0).max())/c['fixed_scales'][field]
                      for key,field in [('active','active'),('field_inactive','inactive')])
            assert drift<1e-3 and abs(z['recovery_orientation'][-1])>.4
            result.update(final_full_field_match=match,field_tail_range_scaled=drift,
                          outcome='original pole' if j['duration']==20 else 'reversed pole')
        audits.append(result);trials[slot]=z
    record=dict(passed=True,generated_at=datetime.now(timezone.utc).isoformat(),cases=audits,
        checks='Checkpoint identity, archive hashes, finite nonnegative fields, conserved total, independent native polarity, phase continuity and fully relaxed polarized withdrawal endpoints.')
    (RUN/'independent_audit.json').write_text(json.dumps(record,indent=2)+'\n')
    (RUN/'status.json').write_text(json.dumps(dict(generated_at=record['generated_at'],complete=True,producer='finished',jobs=8,audit_passed=True),indent=2)+'\n')
    checks['new_trials']=record
    return trials


def duration_figure():
    campaign=ROOT/'data/runs/debelly_corrected_three_test_20261008'
    status=read(campaign/'results/status.json');assert status['complete']
    assert read(campaign/'results/independent_audit.json')['passed']
    rs=fb.read_endpoint_rows(campaign/'results/endpoints.csv');register(campaign/'results/endpoints.csv')
    fig=plt.figure(figsize=(11.5,14.7),layout='constrained')
    grid=fig.add_gridspec(6,15,height_ratios=[1,1,1,1,1,1.05])
    axs=np.array([[fig.add_subplot(grid[i,5*j:5*j+5]) for j in range(3)] for i in range(5)])
    decreases=[]
    tests=['fixed_physical','fixed_relative','fixed_slope']
    for i,m in enumerate(fb.MODELS):
        for j,g in enumerate(tests):
            ax=axs[i,j];panel(ax,chr(65+3*i+j))
            chosen=sorted([r for r in rs if r['model']==m and r['condition']==g],key=lambda r:r['length'])
            negative=[];valid=[];caps=[]
            for gate,color,marker,_ in fb.GATES:
                rr=[r for r in chosen if r['definition']==gate]
                x=np.array([r['length'] for r in rr]);y=np.array([r['endpoint'] if r['status']=='success' else np.nan for r in rr])
                ax.plot(x,y,color=color,ls='-' if gate=='orientation' else '--',lw=1.4)
                for a,b in zip(rr[:-1],rr[1:]):
                    if a['status']==b['status']=='success' and b['endpoint']<a['endpoint']:
                        negative.append((a['length'],b['length']))
                        decreases.append(dict(model=m,condition=g,gate=gate,length_start=a['length'],length_end=b['length'],
                            duration_start=a['endpoint'],duration_end=b['endpoint'],secant=(b['endpoint']-a['endpoint'])/(b['length']-a['length'])))
                        ax.plot([a['length'],b['length']],[a['endpoint'],b['endpoint']],color=color,lw=3,alpha=.75)
                for r in rr:
                    if r['status']=='success':
                        valid.append(r['endpoint']);ax.plot(r['length'],r['endpoint'],marker,color=color,ms=3.5,mfc=color if r['selected_field_tail_pass'] else 'white')
                    elif r['status']=='censored':
                        caps.append(r['maximum_tested_duration']);ax.plot(r['length'],r['maximum_tested_duration'],'^',ms=4,color=color,mfc='white')
                    elif gate=='orientation':ax.plot(r['length'],-.16,'x',color='#888888',transform=ax.get_xaxis_transform(),clip_on=False,ms=4)
            spans=[]
            for a,b in sorted(set(negative)):
                if spans and a<=spans[-1][1]:spans[-1]=(spans[-1][0],max(b,spans[-1][1]))
                else:spans.append((a,b))
            for a,b in spans:ax.axvspan(a,b,color='#909090',alpha=.12,zorder=-1,lw=0)
            ax.set_yscale('log');ax.set_ylim(min(valid)/1.5,max(valid+caps)*1.6)
            lo,hi=ax.get_ylim()
            if hi/lo>=20:
                ax.yaxis.set_major_locator(LogLocator(base=10,numticks=4))
            else:
                ax.yaxis.set_major_locator(MaxNLocator(nbins=3,min_n_ticks=2))
                ax.yaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.tick_params(labelsize=8)
            if i==0:ax.set_title(['Fixed width 3 (Spring 0.2)','Relative width 20%','Constant physical slope'][j],pad=12)
            if j==0:ax.set_ylabel(fb.LABELS[m]+'\nCue duration')
            ax.set_xlabel('Resting length L0' if m=='spring' else 'Cell length L',labelpad=14)
    handles=[Line2D([],[],color=c,ls='-' if g=='orientation' else '--',marker=m,label=label) for g,c,m,label in fb.GATES]
    handles += [Patch(facecolor='#eeeeee',label='Sampled decreasing interval'),Line2D([],[],color='#888888',marker='o',mfc='white',ls='none',label='Field flag fails'),Line2D([],[],color='#888888',marker='^',mfc='white',ls='none',label='Tested limit'),Line2D([],[],color='#888888',marker='x',ls='none',label='Ineligible')]
    transport=transport_measurements()
    for col,model in enumerate(fb.MODELS):
        ax=fig.add_subplot(grid[5,3*col:3*col+3]);panel(ax,'PQRST'[col],fb.LABELS[model])
        rr=[r for r in transport if r['model']==model]
        ax.plot([r['length'] for r in rr],[r['full_active_W1'] for r in rr],color='#0072B2',lw=1.7)
        ax.plot([r['length'] for r in rr],[r['excess_active_W1'] for r in rr],color='#D55E00',ls='--',lw=1.7)
        ax.set(xlabel='Resting L0' if model=='spring' else 'Cell length L',ylim=(0,None))
        ax.tick_params(labelsize=8)
        if col==0:ax.set_ylabel('Mean reflection\ntransport distance')
    handles += [Line2D([],[],color='#0072B2',label='P–T: full active field'),Line2D([],[],color='#D55E00',ls='--',label='P–T: above-background field')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    csvwrite(DATA/'negative_intervals.csv',decreases);measurements['negative_intervals']=decreases
    save(fig,'length_responses','Figure 1. A–O: five models (rows), three cue geometries (columns). Blue/orange denote orientation/strength gates. Duration axes are logarithmic. Gray bands and heavier segments mark sampled decreases; open markers retain field-tail failures, triangles tested limits, and crosses ineligible preparations. P–T: mean earth mover distance W1 from each initial active profile to its reflected target; solid blue includes basal background and dashed orange subtracts the field minimum. Distance is per unit transported active material, in native length units; Spring uses its physical domain and deBelly nominal coordinates.',
         'Unchanged 660 endpoint rows and 97 exact initial checkpoints. W1 is the integral of the absolute difference between normalized cumulative mass distributions, evaluated with native quadrature (finite-volume weights, or trapezoidal Spring node weights) and independently checked with scipy.stats.wasserstein_distance. The reflected target has identical active amount. Amount×W1 is also exported. This describes a target transport burden, not the realized route of active/inactive chemical interconversion. Gray bands mark decreases in either gate.')


def transport_measurements():
    from scipy.stats import wasserstein_distance
    rr=rows(ROOT/'data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv');result=[]
    for row in rr:
        model=row['model'];L=float(row['domain_length']);physical=float(row['physical_length']);cp=load(ROOT/row['checkpoint'])
        assert sha(ROOT/row['checkpoint'])==row['checkpoint_sha256']
        if model=='spring':
            u=cp['state'][:-2:2];x=np.linspace(0,physical,len(u));w=np.full(len(u),physical/(len(u)-1));w[[0,-1]]*=.5
        else:
            u=cp['state_rac'] if model=='debelly' else cp['state'][0]
            x=(np.arange(len(u))+.5)*L/len(u);w=np.full(len(u),L/len(u))
        measured={}
        for name,y in [('full_active',u),('excess_active',u-u.min())]:
            mass=y*w;target=mass[::-1];amount=float(mass.sum());assert amount>0
            cdf=np.cumsum(mass/amount-target/amount)
            W1=float(np.sum(abs(cdf[:-1])*np.diff(x)))
            np.testing.assert_allclose(W1,wasserstein_distance(x,x,u_weights=mass,v_weights=target),atol=2e-12,rtol=2e-12)
            measured[name+'_W1']=W1;measured[name+'_amount']=amount;measured[name+'_amount_times_W1']=W1*amount
        result.append(dict(model=model,length=L,physical_length=physical,**measured))
    csvwrite(DATA/'transport_distance.csv',result);measurements['transport']=result
    checks['transport_W1_independently_verified']=True
    return result


def initial_figure():
    rs=rows(ROOT/'data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv')
    p=read(ROOT/'data/runs/initial_polarity_length_corrected_20261008/report_provenance.json')
    fig,axs=plt.subplots(4,5,figsize=(12,8.5),sharex='col',layout='constrained')
    metrics=[('old_pole_total_active','old_pole_excess_active','Active material'),('old_pole_width',None,'Pole width'),
        ('old_pole_mean_active','old_pole_peak_active','Active concentration'),('opposite_inactive',None,'Inactive concentration')]
    for col,m in enumerate(fb.MODELS):
        rr=sorted([r for r in rs if r['model']==m],key=lambda r:float(r['domain_length']));x=[float(r['domain_length']) for r in rr]
        all_lengths=x+[float(v) for v in p['excluded_lengths'][m]]
        padding=.04*(max(all_lengths)-min(all_lengths))
        for row,(key,extra,label) in enumerate(metrics):
            ax=axs[row,col]
            if col==0:panel(ax,'ABCD'[row]);ax.set_ylabel(label)
            else:ax.grid(alpha=.16)
            if row==0:ax.set_title(fb.LABELS[m],pad=12)
            if key=='opposite_inactive' and m=='debelly':ax.set_yticks([])
            else:
                primary=[float(r[key]) for r in rr]
                secondary=[float(r[extra]) for r in rr] if extra else []
                ax.plot(x,primary,color='#0072B2',marker='o',ms=2.5,lw=1.4)
                if extra:ax.plot(x,secondary,color='#D55E00',ls='--',lw=1.4)
                ax.set_ylim(0,1.08*max(primary+secondary))
            ex=p['excluded_lengths'][m]
            ax.scatter(ex,[-.05]*len(ex),transform=ax.get_xaxis_transform(),marker='x',color='#888888',clip_on=False,s=15)
            ax.set_xlim(min(all_lengths)-padding,max(all_lengths)+padding)
            if row==3:ax.set_xlabel('Resting length L0' if m=='spring' else 'Cell length L',labelpad=14)
    fig.legend(handles=[Line2D([],[],color='#0072B2',label='A: total; C: mean'),Line2D([],[],color='#D55E00',ls='--',label='A: above background; C: peak'),Line2D([],[],color='#888888',marker='x',ls='none',label='Ineligible')],loc='outside lower center',ncol=3,frameon=False)
    save(fig,'initial_polarity_length','Figure 2. Exact initial basal states: active pole material (A), connected half-contrast width (B), pole mean and peak concentration (C), and opposite-end inactive concentration (D). Rows share panel letters across the five model columns. deBelly has no explicit inactive Rac field; its D axis is empty. Units are model-native; amount is a one-dimensional concentration integral.',
         'All 97 current checkpoint measurements, including 12 corrected deBelly states; values unchanged. Width and amount use the boundary-connected half-contrast core. Spring uses its physical domain and resting length on the horizontal axis.')


def chemical_figure():
    z=load(POOL/'evidence_D0p1.npz');pp=read(POOL/'provenance.json')
    for p,h in pp['inputs'].items():assert sha(POOL/p)==h
    fig,axs=plt.subplots(4,2,figsize=(11.8,14.1),layout='constrained');quant=[]
    for col,(m,(label,model,mc)) in enumerate(MODELS.items()):
        register(ROOT/f'code/core/models/{m}.py')
        ax=axs[0,col];right=ax.twinx();right.spines['right'].set_visible(True)
        for L,c in LC.items():
            u,v=z[f'{m}_L{L}_steady'];x=(np.arange(len(u))+.5)*L/len(u)
            ax.plot(x,u,color=c,lw=1.8);right.plot(x,v,color=c,ls='--',lw=1.5)
        panel(ax,'AB'[col],label);ax.set(xlabel='Position x',ylabel='Active u (solid)',xlim=(0,40),ylim=(0,None));right.set(ylabel='Inactive v (dashed)',ylim=(0,None))
        umax=max(z[f'{m}_L{L}_steady'][0].max() for L in LC)*1.08
        vmax=max(z[f'{m}_L{L}_steady'][1].max() for L in LC)*1.15
        ug=np.linspace(0,umax,501);vg=np.linspace(0,vmax,401);U,V=np.meshgrid(ug,vg);f=native_f(model,U,V)
        # Verify native calls with a separately written closed-form expression.
        direct=(model.ec*V*(model.a*U**2+model.b*U)-model.c*U if m=='goryachev' else V*(model.k0+model.feedback*U**2/(model.half_saturation**2+U**2))-model.decay*U)
        np.testing.assert_allclose(f,direct,atol=2e-14,rtol=2e-13)
        ax=axs[1,col]
        # The same quantitative map applies to all three length trajectories.
        # A signed logarithmic map reveals native rates near the incumbent without
        # clipping the larger rates elsewhere in the concentration plane.
        limit=10**np.ceil(np.log10(abs(f).max()))
        linear=.01 if m=='goryachev' else .001
        norm=SymLogNorm(linthresh=linear,linscale=1,vmin=-limit,vmax=limit,base=10)
        mesh=ax.pcolormesh(U,V,f,norm=norm,cmap='RdBu_r',shading='auto',rasterized=True)
        ticks=np.r_[-10.**np.arange(np.log10(limit),np.log10(linear)-.1,-1),0,
                     10.**np.arange(np.log10(linear),np.log10(limit)+.1,1)]
        cb=fig.colorbar(mesh,ax=ax,orientation='horizontal',pad=.13,fraction=.06,aspect=35,ticks=ticks)
        cb.ax.tick_params(labelsize=8)
        cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
        cb.set_label('Basal net activation rate f(u, v)',fontsize=9)
        ax.contour(U,V,f,levels=[0],colors=['#333333'],linestyles=['--'],linewidths=1.2)
        for L,c in LC.items():
            u,v=z[f'{m}_L{L}_steady'];line,=ax.plot(u,v,color=c,lw=2.1)
            line.set_path_effects([pe.Stroke(linewidth=3.6,foreground='white'),pe.Normal()])
            ax.plot(u[0],v[0],'^',color=c,ms=6,mec='white');ax.plot(u[-1],v[-1],'o',color=c,ms=5,mec='white')
            k=len(u)//5
            ax.annotate('',xy=(u[k+5],v[k+5]),xytext=(u[k],v[k]),arrowprops=dict(arrowstyle='->',color=c,lw=1.3))
        panel(ax,'CD'[col]);ax.set(xlabel='Active u(x)',ylabel='Inactive v(x)',xlim=(0,umax),ylim=(0,vmax))
        for L,c in LC.items():
            u,v=z[f'{m}_L{L}_steady'];x=(np.arange(len(u))+.5)*L/len(u);fr=native_f(model,u,v)
            axs[2,col].plot(x,fr,color=c,lw=1.8)
            frozen=native_f(model,ug,np.full_like(ug,v[-1]));axs[3,col].plot(ug,frozen,color=c,lw=2.1 if L==9 else 1.8,ls={9:'-',21:'--',39:':'}[L])
            if m=='goryachev':threshold=(model.c/(model.ec*v[-1])-model.b)/model.a
            else:
                roots=np.roots([model.decay,-v[-1]*(model.k0+model.feedback),model.decay*model.half_saturation**2,-v[-1]*model.k0*model.half_saturation**2])
                threshold=float(sorted(r.real for r in roots if abs(r.imag)<1e-8)[1])
            axs[3,col].plot(threshold,0,'o',color=c,ms=5)
            quant.append(dict(model=m,length=L,inactive_at_new_end=float(v[-1]),unstable_local_threshold=threshold,
                              incumbent_peak=float(u.max()),max_incumbent_net_activation=float(fr.max()),min_incumbent_net_activation=float(fr.min()),
                              cue_free_rate_at_u1=float(native_f(model,np.array(1.),np.array(v[-1])))))
            evidence[f'{m}_L{L}_frozen_active']=ug;evidence[f'{m}_L{L}_frozen_reaction']=frozen
            evidence[f'{m}_L{L}_spatial_reaction']=fr
        ax=axs[2,col];panel(ax,'EF'[col]);ax.axhline(0,color='#888888',lw=.7);ax.set(xlabel='Position across the initial pole, x',ylabel='Basal net activation rate f',xlim=(0,5 if m=='goryachev' else 16))
        ax=axs[3,col];panel(ax,'GH'[col]);ax.axhline(0,color='#888888',lw=.7);ax.set(xlabel='Competing peak active concentration u',ylabel='Cue-free net activation rate f',xlim=(0,65 if m=='goryachev' else umax))
        if m=='goryachev':
            # A linear rate axis makes the increasingly deep/extended loss
            # region directly comparable, rather than visually compressing it.
            ax.set_ylim(-.32,2.6)
        evidence[f'{m}_rate_map_active']=ug;evidence[f'{m}_rate_map_inactive']=vg
        evidence[f'{m}_rate_map_f']=f
    handles=[Line2D([],[],color=c,lw=2,label=f'L = {L}') for L,c in LC.items()]
    handles += [Line2D([],[],color='#333333',ls='--',label='C–D: zero net reaction'),Line2D([],[],color='#333333',marker='^',ls='none',label='C–D: original pole'),Line2D([],[],color='#333333',marker='o',ls='none',label='C–D: opposite end; G–H: unstable threshold')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    csvwrite(DATA/'chemical_measurements.csv',quant);measurements['chemical']=quant
    save(fig,'chemical_hurdles','Figure 3. Goryachev (left) and Wave-pinning (right). A–B: active u (solid, left axis) and inactive v (dashed, right axis) profiles on physical position axes; concentration axes are separately scaled. C–D: colored paths trace the paired concentrations (u(x),v(x)) across the same initial poles; triangles mark the incumbent end, circles the opposite end and arrows point along increasing x. The quantitative background gives native basal net activation f(u,v): red is activation, blue deactivation, dashed black zero. Numeric colorbars retain absolute rate amplitude with a signed logarithmic normalization, linear within ±0.01 (C) or ±0.001 (D). The maps have model-specific native scales and are shared across lengths. E–F: basal reaction along each spatial profile, on linear rate axes, excluding diffusion. G–H: cue-free reaction of a competing peak at the original opposite-end inactive concentration; dots mark unstable local thresholds. G focuses on u = 0–65 on linear axes. Length colors are shared throughout; nearly coincident Wave-pinning curves overlap.',
        'Exact native reaction calls independently checked against closed-form kinetics. All C–D pixels contain f(u,v), not its sign or a per-length normalized rate. G–H hold the inactive pool and incumbent fixed and ignore diffusion while the new peak forms; this is a frozen-substrate local diagnostic, not an evolved competing peak or a spatial survival criterion. E–F use the spatially varying native v(x). All rates are model-native concentration per native time.')
    checks['chemical_closed_form']=True


def restoration_figure():
    z=load(ROOT/'data/runs/length_mechanisms_20261008/comparative_revision/compact_evidence.npz')
    rs=rows(POOL/'response_vs_length.csv')
    assert read(ROOT/'data/runs/length_mechanisms_20261008/comparative_revision/independent_observable_audit.json')['passed']
    fig,axs=plt.subplots(3,2,figsize=(11.5,9.5),layout='constrained')
    for col,(m,(label,model,mc)) in enumerate(MODELS.items()):
        for L,c in LC.items():
            p=f'{m}_L{L}_';axs[0,col].plot(z[p+'time'],100*z[p+'fractional_core_damage'],color=c,lw=1.8)
            axs[1,col].plot(z[p+'recovery_time'],z[p+'reaction_per_initial_damage'],color=c,lw=1.8)
        ax=axs[0,col];panel(ax,'AB'[col],label);ax.axvspan(-.1,0,color='#dddddd');ax.axvline(0,color='#888888',lw=.7)
        ax.set_xscale('symlog',linthresh=.1);ax.set(xlim=(-.1,100),xlabel='Time from inactive return',ylabel='Active core material lost (%)')
        ax.set_xticks([-.1,0,.1,1,10,100],['−0.1','0','0.1','1','10','100'])
        ax=axs[1,col];panel(ax,'CD'[col]);ax.set(xlim=(0,2),xlabel='Time after inactive return',ylabel='Excess reaction\n/ damage at return')
        rr=sorted([r for r in rs if r['model']==m and float(r['duration'])==.1 and float(r['epsilon'])==.005],key=lambda r:float(r['length']))
        axs[2,0].plot([float(r['length']) for r in rr],[100*float(r['fractional_damage']) for r in rr],'o-',color=mc,lw=1.8,ms=4,label=label)
        axs[2,1].plot([float(r['length']) for r in rr],[float(r['reaction_gain_per_damage']) for r in rr],'o-',color=mc,lw=1.8,ms=4,label=label)
    panel(axs[2,0],'E');axs[2,0].set(xlabel='Cell length L',ylabel='Active core material lost (%)',yscale='log')
    panel(axs[2,1],'F');axs[2,1].set(xlabel='Cell length L',ylabel='Immediate excess reaction\n/ damage at return')
    handles=[Line2D([],[],color=c,label=f'A–D: L = {L}') for L,c in LC.items()]
    handles += [Line2D([],[],color=mc,label=f'E–F: {label}') for label,model,mc in MODELS.values()]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    save(fig,'incumbent_restoration','Figure 4. Active-core loss through depletion and return (A–B), post-return restoring reaction per unit initial damage (C–D), and length summaries of damage and immediate restoring response (E–F). Gray shading marks 0.5% whole-inactive-profile removal for 0.1 native units, returned exactly at time zero. Colors in A–D encode length; colors in E–F encode model. Reaction excludes diffusion.',
         'Unchanged independently audited native perturbation histories and 18-length summaries. Fixed original half-height core; matched sham subtraction. Immediate reaction gain does not measure late recovery time.')


def sink_figure():
    rs=rows(SINK/'results/peak_and_amount_measurements.csv');m=read(SINK/'manifest.json')
    assert read(SINK/'results/audit.json')['passed']
    selected=[r for r in rs if float(r['strength'])==1000 and float(r['length'])>=9]
    fig,axs=plt.subplots(4,2,figsize=(11.5,13),layout='constrained');penalties=[];plateau=[];snapshots=[]
    for name,(label,model,color) in MODELS.items():
        for g,ls in [('fixed_gap','-'),('opposite','--')]:
            rr=sorted([r for r in selected if r['model']==name and r['geometry']==g],key=lambda r:float(r['length']))
            x=np.array([float(r['length']) for r in rr]);amount=np.array([float(r['fractional_loss_rate_10_to_30']) for r in rr])
            peak=np.array([float(r['peak_fractional_loss_rate_10_to_30']) for r in rr])
            axs[0,0].plot(x,peak,ls,color=color,lw=1.8);axs[0,1].plot(x,amount,ls,color=color,lw=1.8)
            for r0 in rr:
                low=next(r for r in rs if r['model']==name and r['geometry']==g and r['length']==r0['length'] and float(r['strength'])==100.)
                plateau.append(dict(model=name,length=float(r0['length']),geometry=g,amount_increase_100_to_1000=float(r0['fractional_loss_rate_10_to_30'])/float(low['fractional_loss_rate_10_to_30'])-1,peak_increase_100_to_1000=float(r0['peak_fractional_loss_rate_10_to_30'])/float(low['peak_fractional_loss_rate_10_to_30'])-1))
        for L in x:
            a=next(r for r in selected if r['model']==name and float(r['length'])==L and r['geometry']=='fixed_gap')
            b=next(r for r in selected if r['model']==name and float(r['length'])==L and r['geometry']=='opposite')
            tf=.2/float(a['fractional_loss_rate_10_to_30']);to=.2/float(b['fractional_loss_rate_10_to_30'])
            penalties.append(dict(model=name,length=L,fixed_gap_10_to_30_interval=tf,opposite_10_to_30_interval=to,extra_interval=to-tf,extra_over_fixed=(to-tf)/tf))
        pp=[r for r in penalties if r['model']==name]
        axs[1,0].plot(x,[r['fixed_gap_10_to_30_interval'] for r in pp],'-',color=color,lw=2)
        axs[1,0].plot(x,[r['opposite_10_to_30_interval'] for r in pp],'--',color=color,lw=2)
        axs[1,1].plot(x,[r['extra_interval'] for r in pp],color=color,lw=2)
    for col,key in enumerate(['Peak-height','Active-amount']):
        panel(axs[0,col],'AB'[col]);axs[0,col].set(yscale='log',ylim=(1e-4,.5),xlabel='Cell length L',ylabel=key+' loss rate\n(initial fraction / time)')
    panel(axs[1,0],'C');axs[1,0].set(yscale='log',xlabel='Cell length L',ylabel='10–30% active-loss interval\n(native time)')
    panel(axs[1,1],'D');axs[1,1].set(yscale='log',xlabel='Cell length L',ylabel='Added interval: opposite − fixed gap\n(native time)')
    for col,(name,(label,model,color)) in enumerate(MODELS.items()):
        for row,L in enumerate([9,39]):
            job=next(j for j in m['jobs'] if j['case']['model']==name and j['case']['length']==L and j['geometry']=='fixed_gap')
            p=SINK/f'workers/{job["slot"]:03d}/k_12.npz';z=load(p);meta=read(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256']
            ax=axs[row+2,col];right=ax.twinx();right.spines['right'].set_visible(True)
            r0=next(r for r in selected if r['model']==name and float(r['length'])==L and r['geometry']=='fixed_gap')
            for target,c in [(0.,'#777777'),(float(r0['t10_loss']),'#0072B2'),(float(r0['t30_loss']),'#D55E00')]:
                # Interpolate native profiles to the same interpolated loss crossing.
                k=int(np.searchsorted(z['time'],target));i=max(0,k-1);k=min(k,len(z['time'])-1)
                weight=0 if k==i else (target-z['time'][i])/(z['time'][k]-z['time'][i])
                u=z['active'][i]*(1-weight)+z['active'][k]*weight;v=z['inactive'][i]*(1-weight)+z['inactive'][k]*weight
                ax.plot(z['x'],u,color=c,lw=1.7);right.plot(z['x'],v,color=c,ls='--',lw=1.5)
                snapshots.append(dict(model=name,length=L,time=target,core_loss_fraction=0 if target==0 else .1 if target==float(r0['t10_loss']) else .3))
            ax.axvspan(job['sink_left'],job['sink_right'],color='#dddddd',zorder=-1)
            panel(ax,'EFGH'[2*row+col],f'{label}, L = {L}');ax.set(xlim=(0,L),ylim=(0,125 if name=='goryachev' else 1.25),xlabel='Physical position x',ylabel='Active u (solid)')
            ax.set_xticks([0,3,6,9] if L==9 else [0,10,20,30,39])
            right.set(ylabel='Inactive v (dashed)',ylim=(0,.3 if name=='goryachev' else 2.1))
    handles=[Line2D([],[],color=c,label=label) for label,_,c in MODELS.values()]
    handles += [Line2D([],[],color='#333333',ls='-',label='A–C: fixed gap 3'),Line2D([],[],color='#333333',ls='--',label='A–C: opposite end'),Line2D([],[],color='#777777',label='E–H: initial'),Line2D([],[],color='#0072B2',label='E–H: 10% core loss'),Line2D([],[],color='#D55E00',label='E–H: 30% core loss'),Patch(facecolor='#dddddd',label='E–H: sink support')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    csvwrite(DATA/'sink_distance_contrast.csv',penalties);csvwrite(DATA/'sink_saturation.csv',plateau);csvwrite(DATA/'sink_snapshots.csv',snapshots)
    measurements['sink']=dict(penalties=penalties,saturation=plateau,snapshots=snapshots)
    save(fig,'sink_distance_and_resistance','Figure 5. Peak-height (A) and initial-core active-amount (B) loss rates, with common quantitative axes for both models. C: absolute 10–30% active-loss intervals under fixed-gap and opposite-end sinks. D: their paired difference, opposite minus fixed gap, on a common time axis. A–D use logarithmic vertical axes and model-native time. E–H: fixed-gap active/inactive profiles at 0%, 10% and 30% core loss; L = 9 in E–F and L = 39 in G–H. Each row spans the actual cell, with identical physical x axes across models. Concentration axes are fixed across lengths within each model; solid active and dashed inactive use separate scales. Gray shading is the width-1 sink, positioned 3 beyond the original pole edge. Snapshot times (initial, 10%, 30%) are Goryachev L = 9: 0, 7.69, 23.51; L = 39: 0, 40.53, 118.08; Wave-pinning L = 9: 0, 0.82, 1.46; L = 39: 0, 1.75, 3.67.',
         'Width-1 persistent sink from t=0 to 20000 native units; no reversal cue. Inactive equation: dv/dt = Dv d²v/dx² − f(u,v) − k χsink(x) v, with finite-volume overlap fractions at sink boundaries. k=1000, compared with k=100 as a sampled near-plateau check. Loss rate=0.2/(t30−t10). Extra interval=(t30−t10)opposite−(t30−t10)fixed; an exact paired contrast, not an additive causal partition of coupled dynamics. Core and sink position remain fixed.')


def withdrawal_figure(trials):
    fig,axs=plt.subplots(2,2,figsize=(11.5,7),layout='constrained');quant=[]
    for col,slot in enumerate([6,7]):
        z=trials[slot];D=float(z['pulse_time'][-1]);x=z['pulse_x']
        for target,c in [(0,'#777777'),(D/2,'#0072B2'),(D,'#D55E00')]:
            k=np.argmin(abs(z['pulse_time']-target));axs[0,col].plot(x,z['pulse_active'][k],color=c,lw=1.8,label=f'Cue time {z["pulse_time"][k]:g}')
        axs[0,col].axvspan(12,15,color='#eeeeee')
        for target,c in [(0,'#D55E00'),(25,'#CC79A7'),(250,'#0072B2'),(10000,'#222222')]:
            k=np.argmin(abs(z['recovery_time']-target));axs[1,col].plot(x,z['recovery_active'][k],color=c,lw=1.8,label=f'After withdrawal {z["recovery_time"][k]:g}')
            quant.append(dict(duration=D,recovery_time=float(z['recovery_time'][k]),orientation=float(z['recovery_orientation'][k]),peak=float(z['recovery_active'][k].max())))
        for row in range(2):
            ax=axs[row,col];panel(ax,'ABCD'[2*row+col],f'Pulse duration {D:g}' if row==0 else None)
            ax.set(xlim=(0,15),ylim=(0,1.6),xlabel='Physical position x',ylabel='Active concentration u')
    handles=[Line2D([],[],color='#777777',label='A–B: initial'),Line2D([],[],color='#0072B2',label='A–B: halfway through cue'),
             Line2D([],[],color='#D55E00',label='A–D: withdrawal'),Patch(facecolor='#eeeeee',label='A–B: fixed cue support'),
             Line2D([],[],color='#CC79A7',label='C–D: 25 after withdrawal'),Line2D([],[],color='#0072B2',label='C–D: 250 after withdrawal'),
             Line2D([],[],color='#222222',label='C–D: 10000; relaxed endpoint')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,frameon=False)
    csvwrite(DATA/'withdrawal_measurements.csv',quant);measurements['withdrawal']=quant
    save(fig,'wave_pinning_withdrawal','Figure 6. Wave-pinning at L = 15 with fixed physical cue width 3 and peak 0.1. A–B: starting profile, halfway through the pulse and at withdrawal for durations 20 and 60. Gray shading is cue support. C–D: the corresponding post-withdrawal profiles at 0, 25, 250 and 10000. The short pulse restores the original pole; the long pulse reaches a fully relaxed reversed pole, with both fields matching the original or mirrored basal state within 0.1% of their fixed scales.',
         'Two fresh native replays from the exact same basal checkpoint; standard 10000-unit basal recovery. Independent full-field endpoint and final-500 stability checks. The L=15 strength threshold in Figure 1 is 35.5, between the two illustrative doses. All four panels use identical axes.')


def half_mass(u,L):
    dx=L/u.shape[1];cum=np.cumsum(u,axis=1)*dx;half=cum[:,-1]/2
    ix=(cum>=half[:,None]).argmax(axis=1);a=np.arange(len(u));prev=np.where(ix>0,cum[a,np.maximum(ix-1,0)],0)
    result=ix*dx+(half-prev)/u[a,ix]
    # Independent direct finite-volume integral to the returned coordinate.
    w=np.clip(result[:,None]-np.arange(u.shape[1])[None,:]*dx,0,dx)
    np.testing.assert_allclose((u*w).sum(axis=1),half,atol=3e-12,rtol=3e-13)
    return result


def half_mass_figure(trials):
    fig,axs=plt.subplots(2,2,figsize=(11.5,7.8),layout='constrained');quant=[]
    endpoints=rows(ROOT/'data/runs/debelly_corrected_three_test_20261008/results/endpoints.csv')
    for col,(name,(label,model,mc)) in enumerate(MODELS.items()):
        for j,(L,c) in enumerate(LC.items()):
            slot=col*3+j;z=trials[slot];t=z['pulse_time'];u=z['pulse_active'];x=half_mass(u,L)
            displacement=x-x[0];speed=np.diff(x)/np.diff(t);velocity_time=(t[:-1]+t[1:])/2
            # Native PDE rate gives a second, independent local quantile velocity.
            dx=L/u.shape[1];v=z['pulse_field_inactive'];xx=z['pulse_x']
            case=json.loads((RUN/f'workers/{slot:02d}.json').read_text())['job']['case']
            from three_test_reversal import response_cue
            cue=response_cue(case)(0,xx,L)
            reaction=model.reaction(0.,np.tile(xx,len(t)),np.stack([u.ravel(),v.ravel()]),np.tile(cue,len(t)))[0].reshape(u.shape)
            udot=reaction+model.du*laplacian_no_flux(u,dx)
            ix=np.minimum((x/dx).astype(int),len(xx)-1);a=np.arange(len(t))
            cc=np.cumsum(udot,axis=1)*dx;left=np.where(ix>0,cc[a,np.maximum(ix-1,0)],0)
            exactspeed=(.5*cc[:,-1]-left-(x-ix*dx)*udot[a,ix])/u[a,ix]
            axs[0,col].plot(t,displacement,color=c,lw=1.8,label=f'L = {L}')
            axs[1,col].plot(velocity_time,speed,color=c,lw=1.5)
            threshold=next(r for r in endpoints if r['model']==name and float(r['length'])==L and r['condition']=='fixed_relative' and r['definition']=='orientation')
            duration=float(threshold['endpoint']);field_pass=threshold['selected_field_tail_pass']=='True'
            axs[0,col].plot(duration,np.interp(duration,t,displacement),'o',color=c,ms=7,mfc=c if field_pass else 'white',mew=1.5)
            mirrored=L-2*x[0];ts={}
            for q in [.1,.5,.9]:
                hits=np.flatnonzero(displacement>=q*mirrored)
                ts[str(q)]=float(t[hits[0]]) if len(hits) else None
            quant.append(dict(model=name,length=L,initial_half_mass=float(x[0]),final_half_mass=float(x[-1]),
                t10=ts['0.1'],t50=ts['0.5'],t90=ts['0.9'],max_sampled_velocity=float(speed.max()),
                max_native_instantaneous_velocity=float(exactspeed.max()),time_max_sampled_velocity=float(velocity_time[speed.argmax()]),
                initial_velocity=float(exactspeed[0]),median_velocity_first10=float(np.median(speed[velocity_time<=10])),
                minimum_orientation_cue_duration=duration,selected_field_tail_pass=field_pass))
            for key,val in [('time',t),('half_mass',x),('velocity_time',velocity_time),('velocity',speed),('native_velocity',exactspeed)]:evidence[f'{name}_L{L}_{key}']=val
        panel(axs[0,col],'AB'[col],label);axs[0,col].set(xlabel='Time with cue on',ylabel='Active half-mass displacement Δx50',ylim=(0,40))
        panel(axs[1,col],'CD'[col]);axs[1,col].set(xlabel='Time with cue on',ylabel='Mean half-mass velocity\n(0.5-unit window)',ylim=(-.5,80))
        axs[1,col].set_yscale('symlog',linthresh=.001)
        for ax in axs[:,col]:
            ax.set_xscale('symlog',linthresh=10);ax.set_xlim(0,5120);ax.set_xticks([0,10,100,1000,5120],['0','10','100','1000','5120'])
    fig.legend(handles=[Line2D([],[],color=c,label=f'L = {L}') for L,c in LC.items()]+[Line2D([],[],color='#333333',marker='o',ls='none',label='Minimum orientation pulse'),
                        Line2D([],[],color='#333333',marker='o',mfc='white',ls='none',label='Minimum pulse; field flag fails')],
               loc='outside lower center',ncol=2,frameon=False)
    csvwrite(DATA/'half_mass_measurements.csv',quant);measurements['half_mass']=quant;checks['half_mass_integral']=True
    save(fig,'half_mass_dynamics','Figure 7. Active half-mass displacement from its initial position (A–B), and 0.5-unit mean velocity (C–D), at L = 9, 21 and 39. Circles locate the shortest orientation pulses from Figure 1 on the continuously forced histories; open circles retain failed field-tail flags. Goryachev is left and Wave-pinning right; paired axes have identical scales. The relative-20% cue stays on throughout 5120 native units, with length-independent model-specific peak (6.4 and 0.1). Time is linear through 10 and logarithmic thereafter; velocity is symmetric-logarithmic near zero. x50 divides the full active integral in half and includes basal background.',
         'Six fresh native histories, common 0.5-unit sampling. x50 is the exact finite-volume median, with within-cell linear interpolation; plotted velocity is the displacement over each 0.5-unit window divided by 0.5, located at the window midpoint. Peak window means depend on that declared window; they are not resolved instantaneous maxima. A separate native PDE-derived instantaneous quantile rate at recorded samples is exported. Rapid median shifts between separated poles can occur when their mass shares cross one half; this is not interface speed or direct proof of commitment. Different model-native time and cue scales prevent interpreting absolute cross-model rates as a calibrated biological comparison.')


def main():
    FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
    trials=audit_trials()
    duration_figure();initial_figure();chemical_figure();restoration_figure();sink_figure();withdrawal_figure(trials);half_mass_figure(trials)
    np.savez_compressed(DATA/'plot_evidence.npz',**evidence)
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',passed=True,
        sources=sources,renderer_sha256=sha(Path(__file__)),checks=checks,measurements=measurements,
        figures={str(p.relative_to(FIG)):sha(p) for p in FIG.iterdir() if p.suffix in ['.png','.pdf','.svg']})
    (DATA/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print(FIG)


if __name__=='__main__':main()
