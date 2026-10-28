"""Rebuild R2 from compact CSV records; `extract` reads immutable HPC modes.

No evolution in this entry point. Rate and heading estimators are imported
unchanged. Fixed first replicate supplies every reconstructed exemplar.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,subprocess,sys
import shutil
from pathlib import Path
from datetime import datetime,timezone
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR',f'/tmp/ring-r1r2-mpl-{os.getuid()}')
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import ring_r2_campaign_20261027 as C
import figure4_round3_estimators_20261022 as H
from ring_figure_revision1_20261028 import SIZE_COLORS, pair_rotation, rotate_display

ROOT=C.ROOT;MS=C.MS;RUN=C.RUN;SOURCE=C.SOURCE;OUT=C.OUT
DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
LEVELS=('low','medium','high')
MODELS=('goryachev','wave_pinning','holmes_model3','spring_ring')
NAMES=dict(goryachev='Goryachev',wave_pinning='Wave-pinning',holmes_model3='Holmes',spring_ring='Spring')
R7=MS/'figures/length_mechanisms_20261008/dev/revision_round7_20261026/figure4/supporting'
def setup():
    for p in (DATA,META):p.mkdir(parents=True,exist_ok=True)
def save(name,rows):
    if rows:pd.DataFrame(rows).to_csv(DATA/name,index=False)
def estimate_text(r):
    if r.status=='resolved':return f'α = {r.alpha:.4g}'
    if r.status=='lower_bound':return f'α ≥ {r.alpha_lower:.4g}'
    if r.status=='upper_bound':return f'α ≤ {r.alpha_upper:.4g}'
    if r.status=='execution_pending':return 'α pending'
    return 'α unbounded'
def heading_lag_grid(maxlag,cadence=.5):
    """Keep the initial grid; extend its logarithmic sampling beyond 64,000."""
    transition=64000.
    lags=np.unique(np.r_[cadence,np.geomspace(1,min(transition,maxlag),100)]/cadence).astype(int)*cadence
    if maxlag>transition:
        ratio=transition**(1/99)
        extra=transition*ratio**np.arange(1,int(math.ceil(math.log(maxlag/transition)/math.log(ratio)))+1)
        extra=np.r_[extra[extra<maxlag],maxlag]
        lags=np.r_[lags,np.unique((extra/cadence).astype(int)*cadence)]
    return lags
def latest_level(run,cohort):
    target=int(C.read(run/'expected_levels.json')[str(cohort['slots'][0])])
    for level in range(target,-1,-1):
        if all((C.segment(run,s,level)/'status.json').exists() and C.read(C.segment(run,s,level)/'status.json').get('complete') and (run!=RUN or C.read(C.segment(run,s,level)/'status.json').get('reduced_record')) for s in cohort['slots']):return level
    return None
def record(run,job,level):
    slot=job['slot'];parts=[C.segment(run,slot,k) for k in range(level+1)]
    if run==RUN:
        P=np.concatenate([np.load(p/'P.npy') for p in parts]);t=np.concatenate([np.load(p/'heading_time.npy') for p in parts])
    else:
        a=[np.load(p/'u_modes.npy',mmap_mode='r') for p in parts]
        P=np.concatenate([(z[:,1].conj()/z[:,0]).astype('complex64') for z in a]);t=np.concatenate([np.load(p/'time.npy') for p in parts])
    circles=np.concatenate([np.load(p/'circumference.npy') for p in parts]) if job['model']=='spring_ring' else np.full(len(P),job['size'])
    return t,P,circles
def event_times(job,end):
    """Exact state-independent event sequence from the three frozen streams."""
    streams=np.random.SeedSequence(job['seed']).spawn(3)
    wait,marks,pos=[np.random.default_rng(s) for s in streams]
    rate=job.get('cell_rate',job['rho']*job['size']);t=0.;events=[]
    while True:
        t+=float(wait.exponential())/rate
        if t>end+1e-9:break
        events.append((t,float(marks.exponential(job['amplitude_mean'])),float(pos.uniform())))
    return np.array(events).reshape(-1,3),dict(wait=wait.bit_generator.state,marks=marks.bit_generator.state,positions=pos.bit_generator.state)
class WindowScorer:
    """Original circular-window observable, evaluated with prefix sums."""
    def __init__(self,t,P):
        self.t=t
        unit=P/np.maximum(abs(P),1e-30)
        valid=np.isfinite(unit)
        self.prefix=np.r_[0j,np.cumsum(np.where(valid,unit,0),dtype=np.complex128)]
        self.invalid=np.r_[0,np.cumsum(~valid)]
    def scores(self,centres,w):
        left=np.searchsorted(self.t,centres,side='left');right=np.searchsorted(self.t,centres,side='right')
        valid=(left>=w)&(right+w<=len(self.t));left=left[valid];right=right[valid]
        a=(self.prefix[left]-self.prefix[left-w])/w;b=(self.prefix[right+w]-self.prefix[right])/w
        valid=(self.invalid[left]==self.invalid[left-w])&(self.invalid[right+w]==self.invalid[right])&(abs(a)>1e-12)&(abs(b)>1e-12)
        return 1-np.real(b[valid]*a[valid].conjugate())/(abs(a[valid])*abs(b[valid]))

def event_scores(t,P,events,w):
    return WindowScorer(t,P).scores(events[:,0],w)

def nulls():
    import figure4_round3_analysis_20261022 as A
    rows=[];receipt=[]
    for beta,kind in ((1,'perimeter'),(0,'same_rate'),(-1,'inverse')):
        for size in (20.,60.):
            name=f'A_beta{beta:+d}_C{size:g}.csv.gz';destination=DATA/name
            if beta in (0,1):
                q=pd.read_csv(R7/f'data/A_{kind}_C{size:g}.csv.gz');q['x']=q.pop('relative_arc_position')*size;q.to_csv(destination,index=False)
                receipt.append(dict(beta=beta,size=size,source=str((R7/f'data/A_{kind}_C{size:g}.csv.gz').relative_to(ROOT)),sha256=C.sha(R7/f'data/A_{kind}_C{size:g}.csv.gz')))
            else:
                seed=H.stable_seed('A:inverse');theta=A.simulate_null(1,size,0.,(1/7000)*(size/40)**beta,24000.,20.,seed,bernoulli=True)[0]
                x=(np.arange(240)+.5)*size/240;position=np.mod(theta,2*np.pi)*size/(2*np.pi)
                distance=(x[None,:]-position[:,None]+size/2)%size-size/2
                field=np.exp(-.5*(distance/1.5)**2).astype('float32')
                # Preserve the established compact kymograph column layout.
                q=pd.DataFrame(dict(time=np.repeat(np.arange(len(theta))*20,240),x=np.tile(x,len(theta)),activity=field.ravel()))
                q.to_csv(destination,index=False)
                receipt.append(dict(beta=beta,size=size,seed=seed,generator='figure4_round3_analysis_20261022.simulate_null',transition='original literal Bernoulli lambda*dt, unchanged',lambda0=1/7000))
            for lag in np.linspace(0,24000,241):
                lam=(1/7000)*(size/40)**beta;rows.append(dict(beta=beta,size=size,lag=lag,lambda0=1/7000,lambda_C=lam,C=math.exp(-lam*lag)))
    save('B_null_heading_acf.csv',rows);C.write(META/'null_provenance.json',receipt)

def extract():
    setup();nulls();o=C.legacy();e=o.engine();rates=[];coverage=[];raw_cache={};headings=[];eventrows=[];eventchecks=[];regenerated_logs=[]
    shutil.copy2(SOURCE/'manifest.json',META/'existing_campaign_manifest.json')
    if (RUN/'manifest.json').exists():shutil.copy2(RUN/'manifest.json',META/'per_cell_campaign_manifest.json')
    for src,dest in ((RUN/'validation/cadence.csv',DATA/'record_cadence_validation.csv'),(RUN/'validation/selection.json',META/'record_cadence_selection.json'),(RUN/'results/storage_estimate.json',META/'per_cell_storage_estimate.json')):
        if src.exists():shutil.copy2(src,dest)
    source_state=C.read(SOURCE/'results/status.json') if (SOURCE/'results/status.json').exists() else {}
    for run,arm in ((SOURCE,'per_perimeter'),(RUN,'per_cell')):
        if not (run/'manifest.json').exists():continue
        m=C.read(run/'manifest.json')
        for cohort in m['conditions']:
            job=m['jobs'][cohort['slots'][0]];level=latest_level(run,cohort);key=(arm,cohort['model'],cohort['rate_level'],cohort['amplitude_level'],cohort['size'])
            coverage.append(dict(arm=arm,model=job['model'],size=job['size'],rate_level=job['rate_level'],amplitude_level=job['amplitude_level'],complete_seeds=4 if level is not None else sum((C.segment(run,s,0)/'status.json').exists() and C.read(C.segment(run,s,0)/'status.json').get('complete',False) for s in cohort['slots']),extension_level=level,status='complete_cohort' if level is not None else 'execution_pending'))
            if level is None:continue
            cadence=job['frame_interval']*job.get('mode_stride',1)
            raw=np.array([e.acf_raw(C.series(run,s,level)) for s in cohort['slots']]);curve=raw.mean(0)/raw[:,0].mean()
            g,kind,b=e.half_rate(curve,cadence);lo,hi,los,his,p=o.bootstrap_rates(raw,cadence)
            rec=[record(run,m['jobs'][s],level) for s in cohort['slots']];circles=np.array([r[2].mean() for r in rec])
            raw_cache[key]=(raw,cadence,circles)
            rates.append(dict(arm=arm,model=job['model'],size=job['size'],mean_physical_C=float(circles.mean()),rate_level=job['rate_level'],amplitude_level=job['amplitude_level'],amplitude_mean=job['amplitude_mean'],gamma=g,status=kind,bound=b,gamma_ci_low=lo,gamma_ci_high=hi,ci_low_status=los,ci_high_status=his,analysis_lag=(len(curve)-1)*cadence,crossing_half_lag=math.log(2)/g if g else None,record_duration=job['sample_duration']*2**level,extension_level=level,seeds=4,source_campaign=str(run.relative_to(ROOT))))
            # D uses the existing heading estimator with the original coherence
            # endpoint mask. Its logarithmic display grid spans the record.
            requested=(job['model']=='goryachev' and job['amplitude_level']=='low') or (job['model']=='wave_pinning' and job['amplitude_level']=='low')
            if arm=='per_perimeter' and requested and job['rate_level']=='low' and job['size'] in (20.,60.):
                maxlag=(min(len(t) for t,P,c in rec)//2)*job['frame_interval'];lags=heading_lag_grid(maxlag,job['frame_interval'])
                moments=[H.seed_moments(P,abs(P)>=.1,job['frame_interval'],lags) for t,P,c in rec]
                summary=H.summarize(moments,circles,('R2D',key),lags,lags)
                for ix,lag in enumerate(lags):
                    l,h=H.percentile(summary['Cboot'][:,ix]);headings.append(dict(model=job['model'],size=job['size'],rate_level='low',amplitude_level=job['amplitude_level'],lag=lag,C=summary['C'][ix],C_ci_low=l,C_ci_high=h,coverage=summary['coverage'][ix],eligible_pairs=int(summary['eligible_pairs'][ix]),seeds=4))
                # C: reconstruct first fixed replicate directly from retained
                # Fourier coefficients; no model rerun or seed search.
                slot=cohort['slots'][0];a=np.concatenate([np.load(C.segment(run,slot,k)/'u_modes.npy',mmap_mode='r') for k in range(level+1)])
                tt=rec[0][0];ix=np.unique(np.linspace(0,len(a)-1,321).astype(int));q=(np.arange(240)+.5)/240
                phase=np.exp(2j*np.pi*np.arange(a.shape[1])[:,None]*q)
                field=a[ix,0].real[:,None]+2*np.real(a[ix,1:]@phase[1:]);rr=[]
                for ii,jj in enumerate(ix):
                    for jjx,x in enumerate(q*job['size']):rr.append(dict(time=float(tt[jj]-tt[0]),x=float(x),activity=float(field[ii,jjx])))
                save(f'C_{job["model"]}_C{job["size"]:g}.csv.gz',rr)
            # G uses the existing low-rate, high-amplitude outer-size runs.
            # Chemical times are regenerated exactly; Spring uses logged
            # per-cell events because its original perimeter hazard is dynamic.
            outer=(1.,3.) if job['model']=='spring_ring' else (20.,60.)
            selected=(job['rate_level']=='low' and job['amplitude_level']=='high' and job['size'] in outer and ((arm=='per_perimeter' and job['model']!='spring_ring') or (arm=='per_cell' and job['model']=='spring_ring')))
            if selected:
                windows=(2,3,4,6,8,12,16,24,32,48,64)
                scores={w:[] for w in windows};controls={w:[] for w in windows};counts=[];control_seeds=[]
                for si,s in enumerate(cohort['slots']):
                    jj=m['jobs'][s];t,P,circ=rec[si]
                    if arm=='per_cell':
                        events=np.concatenate([np.loadtxt(C.segment(run,s,k)/'events.csv',delimiter=',',ndmin=2) for k in range(level+1)])
                    else:
                        events,rng_state=event_times(jj,float(t[-1]))
                        regenerated_logs.extend(dict(model=job['model'],size=job['size'],amplitude_level=job['amplitude_level'],replicate=jj['replicate'],seed=jj['seed'],time=float(v[0]),amplitude=float(v[1]),material_position=float(v[2])) for v in events)
                    checkpoint=C.read(C.segment(run,s,level)/'checkpoint.json');counts.append(len(events)==checkpoint['event_count'])
                    if arm=='per_perimeter':
                        retained=np.array(checkpoint['events']);tail=events[-len(retained):] if len(retained) else events[:0]
                        if len(retained) and not np.allclose(tail,retained,rtol=0,atol=1e-6):raise ValueError('Regenerated event tail differs from native checkpoint')
                        if rng_state!=checkpoint['rng']:raise ValueError('Regenerated RNG states differ from native checkpoint')
                        eventchecks.append(dict(arm=arm,model=job['model'],size=job['size'],amplitude=job['amplitude_level'],slot=s,event_count=len(events),checkpoint_event_count=checkpoint['event_count'],first_generated_event_time=float(events[0,0]) if len(events) else None,last_generated_event_time=float(events[-1,0]) if len(events) else None,checkpoint_last_event_time=float(retained[-1,0]) if len(retained) else None,matching_rng_all_three_streams=True,matching_retained_event_tail=True))
                    events=events[(events[:,0]>t[0])&(events[:,0]<t[-1])]
                    scorer=WindowScorer(t,P)
                    control_seed=H.stable_seed(('R2G-random-windows',key,jj['seed']));control_seeds.append(control_seed)
                    for w in windows:
                        observed=scorer.scores(events[:,0],w);scores[w].append(observed)
                        rng=np.random.default_rng(np.random.SeedSequence([control_seed,w]))
                        # Uniform timing over the same run and fully observed
                        # window support. Incidental real events are retained.
                        random_times=rng.uniform(t[w],t[-w-1],100*len(observed))
                        controls[w].append(scorer.scores(random_times,w))
                eventchecks.append(dict(arm=arm,model=job['model'],size=job['size'],amplitude=job['amplitude_level'],all_event_counts_match_checkpoint=all(counts),random_control_seeds=control_seeds,controls_per_event=100,window_frames=list(windows)))
                weights=np.random.default_rng(H.stable_seed(('R2G-excess',key))).multinomial(4,np.full(4,.25),size=2000)
                for w in windows:
                    sums=np.array([v.sum() for v in scores[w]]);ns=np.array([len(v) for v in scores[w]])
                    rs=np.array([v.sum() for v in controls[w]]);rn=np.array([len(v) for v in controls[w]])
                    den=weights@ns;rden=weights@rn
                    boot=np.divide(weights@sums,den,out=np.full(2000,np.nan),where=den>0)-np.divide(weights@rs,rden,out=np.full(2000,np.nan),where=rden>0)
                    l,h=H.percentile(boot);observed=float(sums.sum()/ns.sum()) if ns.sum() else np.nan;baseline=float(rs.sum()/rn.sum()) if rn.sum() else np.nan
                    eventrows.append(dict(arm=arm,model=job['model'],size=job['size'],mean_physical_C=float(circles.mean()),rate_level='low',amplitude_level='high',amplitude_mean=job['amplitude_mean'],window_frames=w,window_native_time=w*job['frame_interval'],mean_heading_change=observed,mean_random_heading_change=baseline,excess=observed-baseline,ci_low=l,ci_high=h,scored_events=int(ns.sum()),random_windows=int(rn.sum()),seeds=4,seed_event_sums=json.dumps(sums.tolist()),seed_event_counts=json.dumps(ns.tolist()),seed_random_sums=json.dumps(rs.tolist()),seed_random_counts=json.dumps(rn.tolist()),control_seeds=json.dumps(control_seeds)))
    save('mode1_rates.csv',rates);save('coverage.csv',coverage);save('D_heading_acf.csv',headings);save('G_event_heading_change.csv',eventrows);save('G_regenerated_event_logs.csv.gz',regenerated_logs);C.write(META/'event_regeneration_checks.json',eventchecks)
    alpha=[];bootstrap_rows=[]
    for arm in ('per_perimeter','per_cell'):
        for model in MODELS:
            if arm=='per_cell' and model not in ('wave_pinning','spring_ring'):continue
            for rate in LEVELS:
                for amplitude in LEVELS:
                    sizes=(1.,2.,3.) if model=='spring_ring' else (20.,40.,60.)
                    keys=[('per_perimeter' if size in (40.,2.) else arm,model,rate,amplitude,size) for size in sizes]
                    base=dict(arm=arm,model=model,rate_level=rate,amplitude_level=amplitude,amplitude_fraction={'low':.3,'medium':1.,'high':3.}[amplitude],seeds_per_size=4)
                    if any(k not in raw_cache for k in keys):alpha.append(dict(base,status='execution_pending'));continue
                    ss=[raw_cache[k] for k in keys];x=np.log([s[2].mean() for s in ss]);a=-(x-x.mean())/np.sum((x-x.mean())**2)
                    point=[e.half_rate(s[0].mean(0)/s[0][:,0].mean(),s[1]) for s in ss]
                    low=[];high=[]
                    for g,k,b in point:
                        low.append(math.log(g if g else b) if k!='beyond_analysis_lag' else -math.inf)
                        high.append(math.log(g if g else b) if k!='below_frame_resolution' else math.inf)
                    lower=sum(ai*(l if ai>0 else h) for ai,l,h in zip(a,low,high) if abs(ai)>1e-14)
                    upper=sum(ai*(h if ai>0 else l) for ai,l,h in zip(a,low,high) if abs(ai)>1e-14)
                    resolved=all(k=='resolved' for _,k,_ in point);value=float(sum(ai*math.log(g) for ai,(g,k,b) in zip(a,point))) if resolved else None
                    rng=np.random.default_rng(H.stable_seed(('R2alpha',arm,model,rate,amplitude)));draws=[];bounded=0;bl=[];bh=[]
                    count_vectors,_=o.bootstrap_counts();count_vectors=count_vectors.astype(int)
                    lookup={tuple(c):ix for ix,c in enumerate(count_vectors)};cached=[]
                    for raw,cad,circle in ss:
                        result=[]
                        for counts in count_vectors:
                            rr=counts@raw;result.append(e.half_rate(rr/rr[0],cad))
                        cached.append(result)
                    # 2000 whole-seed bootstrap resamples, with original gamma
                    # estimator called unchanged on each pooled raw numerator.
                    for draw in range(2000):
                        ys=[];cs=[];okay=True;draw_status=[];rl=[];rh=[]
                        for si,(raw,cad,circle) in enumerate(ss):
                            counts=rng.multinomial(4,np.full(4,.25));g,k,b=cached[si][lookup[tuple(counts)]]
                            if k!='resolved':okay=False
                            ys.append(math.log(g if g else b));cs.append(math.log(counts@circle/4));draw_status.append(k)
                            rl.append(math.log(g if g else b) if k!='beyond_analysis_lag' else -math.inf);rh.append(math.log(g if g else b) if k!='below_frame_resolution' else math.inf)
                        xx=np.array(cs);xx-=xx.mean();coef=-xx/(xx@xx)
                        dl=sum(ai*(ll if ai>0 else hh) for ai,ll,hh in zip(coef,rl,rh) if abs(ai)>1e-14);dh=sum(ai*(hh if ai>0 else ll) for ai,ll,hh in zip(coef,rl,rh) if abs(ai)>1e-14)
                        bl.append(dl);bh.append(dh)
                        if okay:
                            draws.append(-float(xx@ys/(xx@xx)))
                        else:bounded+=1
                        bootstrap_rows.append(dict(arm=arm,model=model,rate_level=rate,amplitude_level=amplitude,draw=draw,alpha_lower=dl,alpha_upper=dh,status='resolved' if okay else 'bounded',rate_statuses=';'.join(draw_status)))
                    l=float(np.quantile(bl,.025,method='inverted_cdf'));h=float(np.quantile(bh,.975,method='inverted_cdf'))
                    status='resolved' if resolved else 'lower_bound' if math.isfinite(lower) and not math.isfinite(upper) else 'upper_bound' if math.isfinite(upper) and not math.isfinite(lower) else 'unbounded' if not math.isfinite(lower) and not math.isfinite(upper) else 'interval_bound'
                    alpha.append(dict(base,status=status,alpha=value,alpha_lower=lower,alpha_upper=upper,alpha_ci_low=l,alpha_ci_high=h,bootstrap_resolved_fraction=len(draws)/2000,rate_statuses=';'.join(k for g,k,b in point),mean_circumferences=';'.join(f'{math.exp(xx):.12g}' for xx in x),analysis_lags=';'.join(str((s[0].shape[1]-1)*s[1]) for s in ss),crossing_lags=';'.join(str(math.log(2)/g if g else None) for g,k,b in point)))
    save('EF_local_alpha.csv',alpha)
    save('EF_whole_seed_bootstrap.csv.gz',bootstrap_rows)
    provenance=dict(generated_utc=C.now(),code_commit=subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True).stdout.strip(),code_sha256={Path(__file__).name:C.sha(__file__),Path(C.__file__).name:C.sha(C.__file__),Path(H.__file__).name:C.sha(H.__file__)},source_manifest=str((SOURCE/'manifest.json').relative_to(ROOT)),source_manifest_sha256=C.sha(SOURCE/'manifest.json'),new_manifest=str((RUN/'manifest.json').relative_to(ROOT)),new_manifest_sha256=C.sha(RUN/'manifest.json') if (RUN/'manifest.json').exists() else None,existing_producer_status=source_state,new_producer_status=C.read(RUN/'results/status.json') if (RUN/'results/status.json').exists() else {'state':'RUNNING'},estimators='Imported unchanged: native acf_raw/half_rate, production.bootstrap_rates; round3 seed_moments/summarize. Alpha is negative three-size OLS slope and 2000 whole-seed resamples.',exemplar_selection='replicate0, fixed before inspection; no new simulation; full retained mode decomposition reconstructed at240 spatial positions',source_records_read_only=True)
    C.write(META/'provenance.json',provenance);documents()

def frame_matrix(path):
    q=pd.read_csv(path)
    if {'time','x'}.issubset(q.columns):
        val=next(c for c in ('activity','u','field','value') if c in q.columns);p=q.pivot(index='time',columns='x',values=val);return p.index.to_numpy(),p.columns.to_numpy(),p.to_numpy()
    return q.iloc[:,0].to_numpy(),np.arange(len(q.columns)-1),q.iloc[:,1:].to_numpy()

def render(export=True):
    setup();import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullLocator, ScalarFormatter
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.5,'axes.linewidth':.6,'lines.linewidth':1,'xtick.major.size':2,'ytick.major.size':2,'svg.fonttype':'none','pdf.fonttype':42})
    fig=plt.figure(figsize=(180/25.4,11.9))
    outer=fig.add_gridspec(5,1,height_ratios=[2.7,2.8,2.0,3.0,2.1],hspace=.53,left=.085,right=.975,bottom=.05,top=.955)
    def label(gs,text):
        b=gs.get_position(fig);fig.text(b.x0-.045,b.y1+.008,text,fontweight='bold',fontsize=10,ha='left',va='bottom')
    def axis(gs,title=None):
        ax=fig.add_subplot(gs);ax.spines[['top','right']].set_visible(False)
        if title:ax.set_title(title,fontsize=7,pad=5)
        return ax
    def key(gs,handles,ncol):
        b=gs.get_position(fig)
        fig.legend(handles=handles,frameon=False,fontsize=5.2,ncol=ncol,loc='lower right',bbox_to_anchor=(b.x1,b.y1+.025),borderaxespad=0,handlelength=1.9,handletextpad=.4,columnspacing=.9)
    def ky(gs,files,sizes,model=None,showx=True,showtitle=True,colorbar=True):
        sub=gs.subgridspec(1,3,width_ratios=[*sizes,1.35],wspace=.12)
        mats=[frame_matrix(p) for p in files]
        rotation=pair_rotation([m[2] for m in mats]);vmax=1. if model is None else max(a[2].max() for a in mats)
        for i,((t,x,v),size) in enumerate(zip(mats,sizes)):
            ax=axis(sub[i]);v=rotate_display(v,rotation)
            # Display clipping only: compact reconstructions and modes stay signed.
            im=ax.imshow(np.maximum(v,0),origin='lower',aspect='auto',extent=[0,size,t[0],t[-1]],cmap='viridis',vmin=0,vmax=vmax,rasterized=True)
            ax.set_ylim(0,min(m[0][-1] for m in mats));ax.set_xticks([0,size]);ax.set_xlabel((f'C{size:g} position' if model is None else 'Position') if showx else '')
            if showtitle:ax.set_title(f'C{size:g}',fontsize=6.5,pad=3)
            if not showx:ax.tick_params(labelbottom=False)
            if i==0:ax.set_ylabel('Time')
            else:ax.tick_params(labelleft=False)
        if colorbar:
            cb=fig.colorbar(im,cax=fig.add_subplot(sub[2]),ticks=[0,vmax]);cb.set_label('Activity',fontsize=5.3,labelpad=1);cb.ax.tick_params(labelsize=5.1,length=1,pad=2);cb.outline.set_linewidth(.35)
            cb.minorticks_off();cb.set_ticks([0,vmax],labels=['0',f'{vmax:.3g}'])
    top=outer[0].subgridspec(1,2,width_ratios=[1.3,1],wspace=.3);label(top[0],'A');label(top[1],'B')
    aa=top[0].subgridspec(3,1,hspace=.65);bb=top[1].subgridspec(3,1,hspace=.65);null=pd.read_csv(DATA/'B_null_heading_acf.csv')
    for row,beta in enumerate((1,0,-1)):
        ky(aa[row],[DATA/f'A_beta{beta:+d}_C{s}.csv.gz' for s in (20,60)],(20,60),showx=row==2,showtitle=False,colorbar=row==2)
        box=aa[row].get_position(fig);fig.text((box.x0+box.x1)/2,box.y1+.006,f'β = {beta}',fontsize=7,ha='center')
        ax=axis(bb[row],f'β = {beta}'+(' (curves coincide)' if beta==0 else ''))
        for size in (20,60):
            q=null[(null.beta==beta)&(null['size']==size)];ax.plot(q.lag,q.C,color=SIZE_COLORS[size],label=f'C{size}',zorder=3 if size==20 else 2)
        ax.set(ylim=(-.03,1.04),xlim=(0,24000),ylabel='C(ℓ)');ax.set_xticks([0,12000,24000]);ax.set_xlabel('Lag ℓ' if row==2 else '');ax.legend(frameon=False,fontsize=5.5,ncol=2,loc='upper right')
        if row<2:ax.tick_params(labelbottom=False)
    second=outer[1].subgridspec(1,2,width_ratios=[1.3,1],wspace=.3);label(second[0],'C');label(second[1],'D');cc=second[0].subgridspec(2,1,hspace=.8);dd=second[1].subgridspec(2,1,hspace=.8)
    heading=pd.read_csv(DATA/'D_heading_acf.csv');guides=[]
    for i,model in enumerate(('goryachev','wave_pinning')):
        files=[DATA/f'C_{model}_C{s}.csv.gz' for s in (20,60)]
        if all(p.exists() for p in files):ky(cc[i],files,(20,60),model=model,showx=i==1)
        else:
            ax=axis(cc[i]);ax.text(.5,.5,'Execution pending',ha='center',transform=ax.transAxes)
        box=cc[i].get_position(fig);fig.text(box.x0,box.y1+.025,NAMES[model],fontsize=7)
        ax=axis(dd[i],NAMES[model])
        for size in (20,60):
            color=SIZE_COLORS[size];q=heading[(heading.model==model)&(heading['size']==size)].sort_values('lag');ax.plot(q.lag,q.C,color=color,label=f'C{size}');ax.fill_between(q.lag,q.C_ci_low,q.C_ci_high,color=color,alpha=.12,linewidth=0)
        q=heading[(heading.model==model)&(heading['size']==20)].sort_values('lag');values=q.C.where(q.C>0)**(1/9)
        ax.plot(q.lag,values,'--',color='.45',linewidth=.8,label='Drift guide')
        guides.extend(dict(model=model,guide='Drift null',lag=l,C=v) for l,v in zip(q.lag,values))
        ax.set(xlim=(0,12000),ylim=(-.1,1.05),ylabel='C(ℓ)',xlabel='Lag ℓ');ax.set_xticks([0,4000,8000,12000])
        handles,labels=ax.get_legend_handles_labels();handles.append(Line2D([],[],color='.5',lw=4,alpha=.2,label='95% seed CI'))
        if i==0:key(second[1],handles,2)
    if export:save('D_analytic_guides.csv',guides)
    alpha=pd.read_csv(DATA/'EF_local_alpha.csv');eg=outer[2].subgridspec(1,4,wspace=.35);label(outer[2],'E')
    rate_styles={'low':('-', 'o',True),'medium':('--','s',True),'high':(':','^',True)}
    arm_styles={'per_perimeter':('-','o',True),'per_cell':('--','o',False)}
    def alpha_plot(gs,q,styles,group,title,ylabel=False,showxlabel=True):
        sub=gs.subgridspec(2,1,height_ratios=[1,.11],hspace=.07);ax=axis(sub[0],title);band=axis(sub[1])
        band.set_facecolor('.96');band.spines[['left','right','top','bottom']].set_visible(False);band.set_yticks([])
        for k,(name,(style,marker,filled)) in enumerate(styles.items()):
            qq=q[q[group]==name].sort_values('amplitude_fraction');resolved=qq[qq.status=='resolved']
            ax.plot(resolved.amplitude_fraction,resolved.alpha,linestyle=style,marker=marker,ms=2.5,color='black',mfc='black' if filled else 'white')
            if len(resolved):
                low=np.maximum(0,resolved.alpha-resolved.alpha_ci_low);high=np.maximum(0,resolved.alpha_ci_high-resolved.alpha)
                ax.errorbar(resolved.amplitude_fraction,resolved.alpha,yerr=[low,high],fmt='none',color='black',linewidth=.6,capsize=1.4)
            for r in qq[qq.status!='resolved'].itertuples(index=False):
                if r.status=='execution_pending':band.plot(r.amplitude_fraction,(k+.5)/len(styles),'x',color='black',ms=3);continue
                low=getattr(r,'alpha_lower',np.nan);high=getattr(r,'alpha_upper',np.nan)
                if np.isfinite(low) and not np.isfinite(high):endpoints=((low,1,'^'),)
                elif np.isfinite(high) and not np.isfinite(low):endpoints=((high,-1,'v'),)
                elif np.isfinite(low) and np.isfinite(high):
                    ax.vlines(r.amplitude_fraction,low,high,color='black',lw=.6);endpoints=((low,1,'^'),(high,-1,'v'))
                else:endpoints=((-1.,-1,'v'),(1.,1,'^'))
                for value,sgn,symbol in endpoints:
                    ax.plot(r.amplitude_fraction,value,symbol,color='black',mfc='white',mec='black',ms=4)
                    ax.annotate('',xy=(r.amplitude_fraction,value+sgn*.8),xytext=(r.amplitude_fraction,value),arrowprops=dict(arrowstyle='->',color='black',lw=.7,shrinkA=0,shrinkB=0,mutation_scale=5.5))
        for y in (2,0,-1):ax.axhline(y,color='.7',lw=.6,ls='--',zorder=0)
        for a in (ax,band):
            a.set_xscale('log');a.set_xlim(.26,3.5);a.set_xticks([.3,1,3],labels=['0.3','1','3']);a.xaxis.set_minor_locator(NullLocator())
        ax.set(ylim=(-4,7));ax.tick_params(labelbottom=False);band.set_ylim(0,1)
        if showxlabel:band.set_xlabel('Amplitude / reference',labelpad=1)
        else:band.tick_params(labelbottom=False)
        if ylabel:ax.set_ylabel('Local α')
        else:ax.tick_params(labelleft=False)
    for i,model in enumerate(MODELS):
        alpha_plot(eg[i],alpha[(alpha.model==model)&(alpha.arm=='per_perimeter')],rate_styles,'rate_level',NAMES[model],ylabel=i==0)
    key(outer[2],[Line2D([],[],color='black',ls=style,marker=marker,ms=3,label=name+' rate') for name,(style,marker,fill) in rate_styles.items()]+[Line2D([],[],color='black',lw=0,marker='^',mfc='white',label='Bound + arrow'),Line2D([],[],color='black',lw=0,marker='x',label='Pending (band)')],5)
    fg=outer[3].subgridspec(2,3,wspace=.28,hspace=.70);label(outer[3],'F')
    for i,model in enumerate(('wave_pinning','spring_ring')):
        for j,rate in enumerate(LEVELS):
            alpha_plot(fg[i,j],alpha[(alpha.model==model)&(alpha.rate_level==rate)],arm_styles,'arm',f'{NAMES[model]}; {rate} rate',ylabel=j==0,showxlabel=i==1)
    key(outer[3],[Line2D([],[],color='black',ls=style,marker=marker,mfc='black' if filled else 'white',ms=3,label=name.replace('_',' ')) for name,(style,marker,filled) in arm_styles.items()],2)
    gg=outer[4].subgridspec(1,4,wspace=.42);label(outer[4],'G');events=pd.read_csv(DATA/'G_event_heading_change.csv');display_segments=[]
    for i,model in enumerate(MODELS):
        ax=axis(gg[i],NAMES[model]);q=events[(events.model==model)&(events.amplitude_level=='high')]
        for size in ((1.,3.) if model=='spring_ring' else (20.,60.)):
            z=q[q['size']==size].sort_values('window_frames');color=SIZE_COLORS[size]
            ax.plot(z.window_frames,z.excess,'o-',ms=2.4,color=color,label=f'Spring {size:g}' if model=='spring_ring' else f'C{size:g}')
            if len(z):
                ax.fill_between(z.window_frames,z.ci_low,z.ci_high,color=color,alpha=.14,lw=0)
                for r in z.itertuples():display_segments.append(dict(model=model,size=size,arm=r.arm,amplitude_level='high',window_frames=r.window_frames,excess=r.excess,source_csv='G_event_heading_change.csv'))
        ax.axhline(0,color='.65',lw=.6);ax.set_xscale('log',base=2);ax.set(xlim=(2,64),xlabel='Window w (frames)');ax.set_xticks([2,4,8,16,32,64],labels=['2','4','8','16','32','64']);ax.xaxis.set_minor_locator(NullLocator());ax.ticklabel_format(axis='y',useOffset=False,style='sci',scilimits=(-2,2))
        ax.legend(frameon=False,fontsize=5.2,loc='upper left')
        if i==0:ax.set_ylabel('Event-triggered excess')
    if export:pd.DataFrame(display_segments).to_csv(DATA/'G_display_line_segments.csv',index=False)
    if export:save('EF_analytic_guides.csv',[dict(amplitude_fraction=x,alpha=y,label=f'α = {y}') for y in (2,0,-1) for x in (.3,3.)])
    if not export:return fig
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'figureR2_main.{ext}',dpi=300)
    plt.close(fig)
    provenance=C.read(META/'provenance.json');provenance.update(rendered_utc=C.now(),renderer_sha256=C.sha(__file__),display_conventions_sha256=C.sha(MS/'code/mechanisms/ring_figure_revision1_20261028.py'),output_sha256={p.name:C.sha(p) for p in sorted(OUT.glob('figureR2_main.*'))})
    C.write(META/'provenance.json',provenance);documents()

def documents():
    setup();rates=pd.read_csv(DATA/'mode1_rates.csv');alpha=pd.read_csv(DATA/'EF_local_alpha.csv');cov=pd.read_csv(DATA/'coverage.csv');events=pd.read_csv(DATA/'G_event_heading_change.csv')
    provenance=C.read(META/'provenance.json')
    both_complete=all(provenance.get(key,{}).get('state')=='COMPLETE' for key in ('existing_producer_status','new_producer_status'))
    unexpected=[];comparisons=[]
    pending=cov[(cov.complete_seeds<4)|(cov.status!='complete_cohort')]
    if len(pending):unexpected.append(f'{len(pending)} conditions have incomplete execution; pending estimates are × on a separate thin band below the α axes, with no assigned α value.')
    if not both_complete:unexpected.append('Campaign execution or extension checks remain unfinished in the extraction provenance; this snapshot is not final inference. The existing dependency finalizer rebuilds R2 again after both producers finish.')
    elif not len(pending):unexpected.append('Both campaign producers report COMPLETE with no incomplete conditions; extension checks are complete. Remaining bounds represent finite-record censoring.')
    bounded=alpha[~alpha.status.isin(['resolved','execution_pending'])]
    unexpected.append(f'{len(bounded)} α estimates are censored. Open triangles and arrows preserve one-sided bounds; paired opposing triangles/arrows show two-sided or unbounded slopes.')
    q=alpha[(alpha.arm=='per_perimeter')&(alpha.rate_level=='low')&(alpha.amplitude_level=='low')]
    expected=dict(goryachev=4.,wave_pinning=1.6,holmes_model3=.9,spring_ring=0.)
    for r in q.itertuples(index=False):
        unexpected.append(f'Pilot check {NAMES[r.model]}: {estimate_text(r)} ({r.status}); earlier check '+('above 4' if r.model=='goryachev' else f'about {expected[r.model]:g}')+f'. Bootstrap resolved fraction {getattr(r,"bootstrap_resolved_fraction",float("nan")):.4g}. No parameter or seed was changed.')
    for model in ('wave_pinning','spring_ring'):
        pp=alpha[(alpha.model==model)&(alpha.arm=='per_perimeter')&(alpha.status=='resolved')];pc=alpha[(alpha.model==model)&(alpha.arm=='per_cell')&(alpha.status=='resolved')]
        matched=pp.merge(pc,on=['model','rate_level','amplitude_level'],suffixes=('_perimeter','_cell'))
        if len(matched):
            delta=matched.alpha_cell-matched.alpha_perimeter
            comparisons.append(f'{NAMES[model]} per-cell minus per-perimeter α ranges from {delta.min():.4g} to {delta.max():.4g} across {len(matched)} matched resolved conditions (earlier check about +1).')
    unexpected.extend(comparisons)
    unexpected.extend([
        'Existing Spring event times depend on its changing physical circumference and cannot be regenerated exactly from seeds alone. G therefore retains the logged per-cell outer-size histories for Spring; chemical models use exact regenerated per-perimeter events. Supply arms differ across these model facets.',
        'The shared Spring C40 (size 2) centre retains its original state-dependent per-perimeter hazard. The per-cell outer sizes instead have constant intensity equal to ρ times the exact size-2 basal physical circumference; the centre is not an exactly constant-intensity trajectory.',
        'The β = −1 sample paths retain the original literal Bernoulli arrival transition at interval 20, while B shows the requested continuous-time analytic exponential. This finite-step mismatch is retained.',
        'The prelaunch cadence selection used 46 completed target-model cohorts, including the shared centre and every rate/amplitude combination at C20. Eight high-rate target cohorts were still executing and were excluded from that selection.'
    ])
    for model in ('goryachev','wave_pinning'):
        low=min(pd.read_csv(DATA/f'C_{model}_C{s}.csv.gz').activity.min() for s in (20,60))
        if low<0:unexpected.append(f'{NAMES[model]} Fourier reconstruction reaches {low:.6g} activity because of its finite mode cutoff. Display values are clipped at zero on the shared linear viridis scale only; signed compact reconstruction data, native Fourier coefficients and all heading calculations remain unchanged.')
    for model in MODELS:
        q=events[events.model==model];sizes=(1.,3.) if model=='spring_ring' else (20.,60.)
        lo=q[q['size']==sizes[0]].set_index('window_frames');hi=q[q['size']==sizes[1]].set_index('window_frames')
        windows=lo.index.intersection(hi.index)
        if len(windows):
            delta=hi.loc[windows,'excess']-lo.loc[windows,'excess']
            unexpected.append(f'G {NAMES[model]}: large-minus-small event-triggered excess ranges from {delta.min():.6g} to {delta.max():.6g} across w=2–64 original frames'+('; its sign changes with window width.' if delta.min()<0<delta.max() else '.'))
        if len(q) and (q.excess<0).any():unexpected.append(f'G {NAMES[model]} has negative excess at {(q.excess<0).sum()} size/window combinations: random windows have larger mean heading change there. All signed estimates and intervals are retained, without clipping or retuning.')
        if set(q['size'])!=set(sizes):unexpected.append(f'G {NAMES[model]} outer-size coverage is incomplete; missing histories are not imputed.')
    if (META/'scheduler_dependency_repair.json').exists():unexpected.append('Two original allocation records remain stuck in failed prolog on hpctpa3pc0067 after successful replacement runs completed. The existing initial controller dependency was repaired to retain every remaining live worker and exclude only those stale records; no simulation or scientific record was cancelled or changed. The exact receipt is scheduler_dependency_repair.json.')
    cadence_postlaunch=META/'record_cadence_postlaunch_summary.json'
    if cadence_postlaunch.exists():
        audit=C.read(cadence_postlaunch)
        unexpected.append(f'The additive cadence audit checked {audit["cohorts_checked"]} complete target cohorts: {audit["passes"]} passed, {len(audit["failures"])} failed, and {len(audit["pending_cohorts"])} were pending at that audit. Frozen prelaunch selection unchanged: {audit.get("frozen_selection_unchanged")}. Full comparisons remain in record_cadence_postlaunch_validation.csv.')
    cadence_final=META/'record_cadence_finalization.json'
    if cadence_final.exists() and not C.read(cadence_final)['passed']:unexpected.append('The final complete-grid cadence validation failed or is incomplete; its exact status remains in record_cadence_finalization.json. No cadence or model parameter was changed to close the gap.')
    numbers=[f'Data generation time: **{provenance["generated_utc"]}**. Producers complete: **{both_complete}**. Exact execution/extension status is in provenance.json and coverage.csv.',
        'A–B: λ₀=1/7000, duration 24000, original sample interval 20, physical peak SD 1.5. Widths follow C20:C60=1:3. β=0 curves coincide.',
        'C: Existing-arm replicate 0 only, all retained active Fourier modes, 240 spatial positions, 321 displayed times. Both models use low rate and low amplitude; no example search or rerun.',
        'D: Whole-profile heading arg P, P=M₁/M₀, endpoint |P|≥0.1, 2000 whole-seed resamples. Dashed guide is C20 correlation to the power 1/9; shared display lag axis is 0–12000, and full retained-record curves remain in CSV.'
    ]
    for panel,arm in (('E','per_perimeter'),('F','per_cell')):
        q=alpha[alpha.arm==arm];numbers.append(f'{panel}: {len(q)} model/rate/amplitude combinations; {sum(q.status=="resolved")} resolved, {sum(q.status=="execution_pending")} pending, {sum(~q.status.isin(["resolved","execution_pending"]))} bounded.')
    for model in MODELS:
        q=events[(events.model==model)&(events.window_frames==2)]
        numbers.append(f'G {NAMES[model]} at w=2: {int(q.scored_events.sum())} event windows, {int(q.random_windows.sum())} random control windows over both outer sizes; four seeds per size.')
    methods=[
        'All quantities use native model units. Size represents WGD-associated enlargement. Observed migration summaries, conditional OU velocity persistence and simulated activity-heading persistence remain distinct readouts.',
        'The frozen Fourier estimator is unchanged: actual-pair-count uncentered mode autocorrelation with raw numerators and powers pooled over four seeds; γ₁ is ln 2 divided by the first linearly interpolated half-height lag. All four seeds extend together when the half-height lag exceeds record/8, through the existing fourfold cap. Original finite-record statuses are preserved.',
        'E/F use the negative OLS slope of log γ₁ against log actual mean physical circumference over three sizes. The 2000 bootstrap draws independently resample whole seeds at each size and apply the unchanged rate estimator to pooled numerators. Signed OLS coefficients propagate log-rate intervals into α bounds; bootstrap censoring is retained. Pending symbols have their own band and no numerical α coordinate.',
        'G retains the original event-window observable. Circular mean headings are computed from the w original stored frames immediately before and after each event, excluding a sample exactly at the centre. Every finite nonzero circular mean enters without a coherence gate. The widths are 2,3,4,6,8,12,16,24,32,48,64 frames. This corresponds to model-specific native times (chemical cadence 0.5; Spring cadence 0.05). Overlapping windows remain eligible.',
        'For each seed and each width, 100 times as many control centres as scorable events are sampled uniformly from the same run’s fully observed window support. These random windows retain incidental real events and therefore estimate the same-run background, rather than an event-free trajectory. Control timing uses independent deterministic seeds recorded with seed-level event/control sufficient statistics in G_event_heading_change.csv. Neither event timing nor the model trajectory is changed.',
        r'The plotted event-triggered excess is the difference between the pooled event mean and pooled random-window mean:',
        r'$$\mathcal{E}(w)=\left\langle1-\cos\Delta\theta_w\right\rangle_{\mathrm{event}}-\left\langle1-\cos\Delta\theta_w\right\rangle_{\mathrm{random}}.\tag{1}$$',
        'Confidence bands are pointwise 95% intervals from 2000 whole-seed bootstrap draws. The same seed multiplicities weight event and random sums/counts, preserving their covariance. They are conditional on the sampled random centres; 100 controls per event reduce timing Monte Carlo variability. Excess remains signed and is an event association, not a causal pathway attribution.',
        'Chemical events are regenerated from the three frozen SeedSequence streams, with exact event counts, retained event tails and all RNG states checked against continuation checkpoints. Spring uses its saved per-cell event logs. No new simulation campaign is launched.',
        'One common angular rotation per kymograph pair is determined from the circular mean of the two initial dominant-domain positions. The same angle is applied at every time to both sizes, rounded only to the display grid. Every kymograph uses linear viridis with a colourbar; the same model within the figure shares one activity scale. Negative Fourier ringing is clipped only for display, with signed plotted-source CSVs unchanged. D uses the same low-rate, low-amplitude conditions as C; the fractional-power guide is missing where the measured C20 correlation is nonpositive.',
        'Rebuild portable data with `python polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render`. Read immutable HPC evidence under SLURM using command `all`. The existing dependency finalizer loads the revised renderer after both campaigns and their extensions finish. Frozen simulation sources, seeds, dynamics, amplitudes and rate estimators are unchanged.'
    ]
    text='# Figure R2 scientific readme\n\n## Unexpected\n\n'+'\n\n'.join(unexpected)+'\n\n## Key numbers by panel\n\n'+'\n\n'.join(numbers)+'\n\n## Methods\n\n'+'\n\n'.join(methods)+'\n\n## Event-triggered excess endpoints\n\n'
    text+='| Model | Size | w (frames) | Event mean | Random mean | Excess | 95% seed CI |\n|---|---:|---:|---:|---:|---:|---|\n'
    for r in events[events.window_frames.isin([2,64])].itertuples(index=False):
        text+=f'| {NAMES[r.model]} | {r.size:g} | {r.window_frames} | {r.mean_heading_change:.6g} | {r.mean_random_heading_change:.6g} | {r.excess:.6g} | [{r.ci_low:.6g}, {r.ci_high:.6g}] |\n'
    (META/'scientific_readme.md').write_text(text)
    result='; '.join(f'{NAMES[r.model]} {estimate_text(r)} ({r.status})' for r in alpha[(alpha.arm=='per_perimeter')&(alpha.rate_level=='low')&(alpha.amplitude_level=='low')].itertuples(index=False))+'.'
    legend=[
        '# Jump supply, event response and size-dependent heading memory',
        '**Encoding.** C20, C40 and C60 use orange #E69F00, green #009E73 and blue #0072B2 respectively; Spring sizes 1,2,3 use the same hues. Model identity is given by facets. Each kymograph uses linear viridis activity with its colourbar; kymographs of the same model share their scale. Every C20/C60 pair receives one common angular rotation to place its initial dominant domain near mid-panel; the same angle is used for both sizes throughout the window.',
        r'**A.** Jump-null kymographs at C20 and C60 with widths proportional to circumference and time upward. β appears above each facet. Uniform relocations arrive at $\lambda(C)=\lambda_0(C/40)^\beta$, $\lambda_0=1/7000$. β=1 represents perimeter supply and constant success; β=0 perimeter supply with success proportional to 1/C (or per-cell supply with constant success); β=−1 perimeter supply with success proportional to 1/C². All quantities use native units.',
        r'**B.** Analytic heading correlation $C(\ell)=\exp[-\lambda(C)\ell]$ for A, in the fixed size colours. The β=0 curves coincide; enlargement speeds loss at β=1 and slows it at β=−1.',
        '**C.** Existing per-perimeter campaign activity reconstructed from stored active Fourier modes at fixed replicate 0. Both Goryachev and Wave-pinning use low rate and low amplitude. Both sizes share one model-specific linear activity scale and sampled time axis. Fourier ringing is clipped at zero for display only; signed reconstruction CSVs and the original Fourier data are unchanged.',
        r'**D.** Heading correlation for the exact conditions in C, pooled over four seeds, with endpoint $|P|\geq0.1$, $P=M_1/M_0$, and pointwise 95% whole-seed bootstrap bands. Orange is C20 and blue C60. Dashed grey is the fixed-physical-drift guide $C_{20}(\ell)^{1/9}$, undefined at nonpositive C20 correlation. Both facets retain the common lag axis 0–12000; full curves extend to the retained-record horizon in the source CSV.',
        r'**E.** Local exponent $\alpha$ against event amplitude/reference for per-perimeter supply, faceted by model. It is minus the three-size OLS slope of log mode-1 decorrelation rate against log actual mean circumference. Rates use the first half-height lag. All series are black: low rate solid/circle, medium dashed/square, high dotted/triangle. Filled markers and intervals denote resolved rates with 95% whole-seed uncertainty; open triangles with arrows denote bounds, including opposing arrows for two-sided or unbounded slopes. Pending estimates are × on the thin band below each axis, never at an assigned α. Grey guides are α=2,0,−1. '+result,
        '**F.** The same exponent compares supply arms, in black: per-perimeter solid with filled circles, per-cell dashed with open circles. Facets identify model and event rate; bounds, pending band, intervals and guides match E. Per-cell outer sizes use the fixed C40 arrival intensity. The shared Spring C40 centre (size 2) keeps its original state-dependent hazard; it is not rerun. '+' '.join(comparisons),
        r'**G.** Event-triggered excess $\mathcal E(w)$, the event-centred mean of $1-\cos\Delta\theta$ minus its random-window mean in the same runs, versus window width w=2–64 original frames on a log axis. Each Δθ compares the circular mean headings of the w frames immediately before and after a centre; a sample exactly at the centre is excluded. Conditions are low rate/high amplitude. Orange and blue denote C20/C60 (Spring 1/3). Chemical models use per-perimeter events regenerated and verified against checkpoints; Spring uses logged per-cell outer-size events. Each seed contributes 100 uniformly timed controls per scorable event at each w, over the fully observed run support, retaining incidental events. Lines join measured excesses; bands are pointwise 95% whole-seed bootstrap intervals, resampling event and control contributions jointly. Zero means equal event and background averages; signed excess and window-dependent size differences remain visible. No shared-centre value is required.'
    ]
    (META/'main_legend.md').write_text('\n\n'.join(legend)+'\n')
    (OUT/'README.md').write_text('# Figure R2\n\n![Figure R2](figureR2_main.png)\n\n[PDF](figureR2_main.pdf) · [SVG](figureR2_main.svg) · [Legend](supporting/metadata/main_legend.md) · [Scientific readme](supporting/metadata/scientific_readme.md) · [Provenance](supporting/metadata/provenance.json)\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=('extract','render','all'));a=ap.parse_args()
    if a.command in ('extract','all'):extract()
    if a.command in ('render','all'):render()
