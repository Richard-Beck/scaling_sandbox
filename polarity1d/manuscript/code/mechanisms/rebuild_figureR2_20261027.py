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
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import ring_r2_campaign_20261027 as C
import figure4_round3_estimators_20261022 as H

ROOT=C.ROOT;MS=C.MS;RUN=C.RUN;SOURCE=C.SOURCE;OUT=C.OUT
DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
LEVELS=('low','medium','high')
MODELS=('goryachev','wave_pinning','holmes_model3','spring_ring')
NAMES=dict(goryachev='Goryachev',wave_pinning='Wave-pinning',holmes_model3='Holmes',spring_ring='Spring')
HUES=dict(wave_pinning='#2878b5',holmes_model3='#36944c',spring_ring='#c59a25',goryachev='#c83b3b')
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
def event_scores(t,P,events,w):
    unit=P/np.maximum(abs(P),1e-30);vals=[]
    for event,amplitude,position in events:
        # Window width is w original saved frames, with an event excluded from
        # either window if it lies exactly on a saved frame.
        left=np.searchsorted(t,event,side='left');right=np.searchsorted(t,event,side='right')
        if left<w or right+w>len(t):continue
        before=unit[left-w:left];after=unit[right:right+w]
        if not (np.isfinite(before).all() and np.isfinite(after).all()):continue
        a=before.mean();b=after.mean()
        if abs(a)<=1e-12 or abs(b)<=1e-12:continue
        vals.append(float(1-np.real(b*a.conjugate())/(abs(a)*abs(b))))
    return np.array(vals)

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
            requested=(job['model']=='goryachev' and job['amplitude_level']=='low') or (job['model']=='wave_pinning' and job['amplitude_level']=='medium')
            if arm=='per_perimeter' and requested and job['rate_level']=='low' and job['size'] in (20.,60.):
                maxlag=min(64000.,rec[0][0][-1]-rec[0][0][0]);lags=np.unique(np.r_[.5,np.geomspace(1,maxlag,100)]/.5).astype(int)*.5
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
            # Chemical old-arm event times do not depend on state. Spring old
            # times do, therefore Spring G uses only new per-cell outer sizes.
            if job['rate_level']=='low' and ((arm=='per_perimeter' and job['model']!='spring_ring') or (arm=='per_cell' and job['model']=='spring_ring')):
                scores={w:[] for w in (2,5,20)};counts=[]
                for si,s in enumerate(cohort['slots']):
                    jj=m['jobs'][s];t,P,circ=rec[si]
                    if arm=='per_cell':
                        events=np.concatenate([np.loadtxt(C.segment(run,s,k)/'events.csv',delimiter=',',ndmin=2) for k in range(level+1)])
                    else:
                        events,rng_state=event_times(jj,float(t[-1]))
                        regenerated_logs.extend(dict(model=job['model'],size=job['size'],amplitude_level=job['amplitude_level'],replicate=jj['replicate'],seed=jj['seed'],time=float(v[0]),amplitude=float(v[1]),material_position=float(v[2])) for v in events)
                    checkpoint=C.read(C.segment(run,s,level)/'checkpoint.json');counts.append(len(events)==checkpoint['event_count'])
                    if arm=='per_perimeter':
                        retained=np.array(checkpoint['events'])
                        tail=events[-len(retained):] if len(retained) else events[:0]
                        if len(retained) and not np.allclose(tail,retained,rtol=0,atol=1e-6):raise ValueError('Regenerated event tail differs from native checkpoint')
                        if rng_state!=checkpoint['rng']:raise ValueError('Regenerated RNG states differ from native checkpoint')
                        eventchecks.append(dict(arm=arm,model=job['model'],size=job['size'],amplitude=job['amplitude_level'],slot=s,event_count=len(events),checkpoint_event_count=checkpoint['event_count'],first_generated_event_time=float(events[0,0]) if len(events) else None,last_generated_event_time=float(events[-1,0]) if len(events) else None,checkpoint_last_event_time=float(retained[-1,0]) if len(retained) else None,matching_rng_all_three_streams=True,matching_retained_event_tail=True))
                    events=events[(events[:,0]>t[0])&(events[:,0]<t[-1])]
                    for w in (2,5,20):scores[w].append(event_scores(t,P,events,w))
                eventchecks.append(dict(arm=arm,model=job['model'],size=job['size'],amplitude=job['amplitude_level'],all_event_counts_match_checkpoint=all(counts)))
                weights=np.random.default_rng(H.stable_seed(('R2G',key))).multinomial(4,np.full(4,.25),size=2000)
                for w in (2,5,20):
                    sums=np.array([v.sum() for v in scores[w]]);ns=np.array([len(v) for v in scores[w]])
                    den=weights@ns;boot=np.divide(weights@sums,den,out=np.full(2000,np.nan),where=den>0);l,h=H.percentile(boot)
                    eventrows.append(dict(arm=arm,model=job['model'],size=job['size'],mean_physical_C=float(circles.mean()),rate_level='low',amplitude_level=job['amplitude_level'],amplitude_mean=job['amplitude_mean'],window_frames=w,window_native_time=w*job['frame_interval'],mean_heading_change=float(sums.sum()/ns.sum()) if ns.sum() else np.nan,ci_low=l,ci_high=h,scored_events=int(ns.sum()),seeds=4))
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

def render():
    setup();import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.colors import to_rgb, AsinhNorm
    from matplotlib.ticker import NullLocator
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.5,'axes.linewidth':.6,'lines.linewidth':1,'xtick.major.size':2,'ytick.major.size':2,'svg.fonttype':'none','pdf.fonttype':42})
    fig=plt.figure(figsize=(180/25.4,11.9));outer=fig.add_gridspec(5,1,height_ratios=[2.7,2.8,1.8,2.8,2.0],hspace=.48,left=.085,right=.975,bottom=.04,top=.965)
    def label(gs,text):
        b=gs.get_position(fig);fig.text(b.x0-.045,b.y1+.008,text,fontweight='bold',fontsize=10,ha='left',va='bottom')
    def axis(gs,title=None):
        ax=fig.add_subplot(gs);ax.spines[['top','right']].set_visible(False)
        if title:ax.set_title(title,fontsize=7,pad=5)
        return ax
    def ramp(hue,n):
        rgb=np.array(to_rgb(hue));return [tuple(1-(1-rgb)*v) for v in np.linspace(.5,1,n)]
    def shared_key(gs,text,names):
        b=gs.get_position(fig)
        fig.text(b.x0,b.y1+.030,text,fontsize=5.7)
        colors=ramp('#333333',len(names))
        handles=[Line2D([],[],color=color,marker='o',ms=2.2,lw=1,label=name) for name,color in zip(names,colors)]
        fig.legend(handles=handles,frameon=False,fontsize=5.0,ncol=len(names),loc='lower right',bbox_to_anchor=(b.x1,b.y1+.025),borderaxespad=0,handlelength=1.5,handletextpad=.4,columnspacing=.9)
    def ky(gs,files,sizes,hue,model=None,showx=True,showtitle=True,colorbar=True):
        # Reserve the identical colourbar column in every facet, including
        # facets without a displayed bar, so physical widths stay 1:3.
        sub=gs.subgridspec(1,3,width_ratios=[*sizes,1.35],wspace=.12);mats=[frame_matrix(p) for p in files];vmax=max(a[2].max() for a in mats);vmin=min(a[2].min() for a in mats) if model else 0.
        axes=[]
        for i,((t,x,v),s) in enumerate(zip(mats,sizes)):
            ax=axis(sub[i]);axes.append(ax)
            norm=AsinhNorm(linear_width=1.,vmin=vmin,vmax=vmax) if model=='goryachev' else None
            im=ax.imshow(v,origin='lower',aspect='auto',extent=[0,s,t[0],t[-1]],cmap='Greys',norm=norm,vmin=None if norm else vmin,vmax=None if norm else vmax,rasterized=True)
            ax.set_ylim(0,min(m[0][-1] for m in mats));ax.set_xticks([0,s]);ax.set_xlabel('Position' if showx else '')
            if showtitle:ax.set_title(f'C{s:g}',fontsize=6.5,pad=3)
            if not showx:ax.tick_params(labelbottom=False)
            if i==0:ax.set_ylabel('Time')
            else:ax.tick_params(labelleft=False)
        if colorbar:
            cb=fig.colorbar(im,cax=fig.add_subplot(sub[2]),ticks=[vmin,vmax]);cb.set_label('Activity (asinh)' if model=='goryachev' else 'Activity',fontsize=5.3,labelpad=1);cb.ax.tick_params(labelsize=5.1,length=1,pad=2);cb.outline.set_linewidth(.35)
            cb.minorticks_off();cb.set_ticks([vmin,vmax],labels=[f'{vmin:.3g}',f'{vmax:.3g}'])
    top=outer[0].subgridspec(1,2,width_ratios=[1.3,1],wspace=.3);label(top[0],'A');label(top[1],'B')
    aa=top[0].subgridspec(3,1,hspace=.55);bb=top[1].subgridspec(3,1,hspace=.55);null=pd.read_csv(DATA/'B_null_heading_acf.csv')
    for row,beta in enumerate((1,0,-1)):
        ky(aa[row],[DATA/f'A_beta{beta:+d}_C{s}.csv.gz' for s in (20,60)],(20,60),'#555555',showx=row==2,showtitle=row==0,colorbar=row==2)
        box=aa[row].get_position(fig);supply={1:'Perimeter × constant',0:'Perimeter × 1/C (or cell × constant)',-1:'Perimeter × 1/C²'}[beta]
        fig.text(box.x0+.005,box.y0+.007,f'β = {beta}: {supply}',fontsize=5.3,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.5))
        ax=axis(bb[row],f'β = {beta}');cols=['#aaaaaa','#333333']
        for color,size in zip(cols,(20,60)):
            q=null[(null.beta==beta)&(null['size']==size)];ax.plot(q.lag,q.C,color=color,label=f'C{size}')
        ax.set(ylim=(-.03,1.04),xlim=(0,24000),ylabel='C(ℓ)');ax.set_xticks([0,12000,24000]);ax.set_xlabel('Lag ℓ' if row==2 else '');ax.legend(frameon=False,fontsize=5.5,ncol=2,loc='upper right')
        if row<2:ax.tick_params(labelbottom=False)
    second=outer[1].subgridspec(1,2,width_ratios=[1.3,1],wspace=.3);label(second[0],'C');label(second[1],'D');cc=second[0].subgridspec(2,1,hspace=.8);dd=second[1].subgridspec(2,1,hspace=.8)
    heading=pd.read_csv(DATA/'D_heading_acf.csv');guides=[]
    for i,model in enumerate(('goryachev','wave_pinning')):
        files=[DATA/f'C_{model}_C{s}.csv.gz' for s in (20,60)]
        if all(p.exists() for p in files):ky(cc[i],files,(20,60),HUES[model],model=model,showx=i==1)
        else:axis(cc[i]).text(.5,.5,'Execution pending',ha='center',transform=fig.axes[-1].transAxes)
        box=cc[i].get_position(fig);fig.text(box.x0,box.y1+.025,NAMES[model]+('; low rate, low amplitude' if i==0 else '; low rate, medium amplitude'),fontsize=6.3)
        ax=axis(dd[i],NAMES[model]+'; low rate');cols=ramp(HUES[model],2)
        for size,color in zip((20,60),cols):
            q=heading[(heading.model==model)&(heading['size']==size)].sort_values('lag');ax.plot(q.lag,q.C,color=color,label=f'C{size}');ax.fill_between(q.lag,q.C_ci_low,q.C_ci_high,color=color,alpha=.12,linewidth=0)
        q=heading[(heading.model==model)&(heading['size']==20)].sort_values('lag');positive=q.C.where(q.C>0)
        for style,power,name in ((':',1.,'No size effect'),('--',(20/60)**2,'Drift null')):
            values=q.C if power==1 else positive**power
            ax.plot(q.lag,values,style,color='#777777',linewidth=.8,label=name)
            guides.extend(dict(model=model,guide=name,lag=l,C=v) for l,v in zip(q.lag,values))
        ax.set(xlim=(0,12000),ylim=(-.1,1.05),ylabel='C(ℓ)',xlabel='Lag ℓ');ax.set_xticks([0,4000,8000,12000]);ax.legend(frameon=False,fontsize=5.1,ncol=2,loc='upper right');ax.text(.02,.03,'Bands: 95% seed CI',fontsize=5,transform=ax.transAxes)
    save('D_analytic_guides.csv',guides)
    alpha=pd.read_csv(DATA/'EF_local_alpha.csv');eg=outer[2].subgridspec(1,4,wspace=.35);label(outer[2],'E')
    def alpha_plot(ax,q,colors,group,legend=False):
        for k,(name,col) in enumerate(colors.items()):
            qq=q[q[group]==name].sort_values('amplitude_fraction');resolved=qq[qq.status=='resolved'];ax.plot(resolved.amplitude_fraction,resolved.alpha,'o-',ms=2.3,color=col,label=name.replace('_',' '))
            if len(resolved):ax.errorbar(resolved.amplitude_fraction,resolved.alpha,yerr=[resolved.alpha-resolved.alpha_ci_low,resolved.alpha_ci_high-resolved.alpha],fmt='none',color=col,linewidth=.6,capsize=1.4)
            for r in qq[qq.status!='resolved'].itertuples(index=False):
                if r.status=='execution_pending':ax.plot(r.amplitude_fraction,-3.5+k*.35,'x',color=col,ms=3);continue
                low=getattr(r,'alpha_lower',np.nan);high=getattr(r,'alpha_upper',np.nan)
                if np.isfinite(low) and not np.isfinite(high):v=low;direction=1;mark='^'
                elif np.isfinite(high) and not np.isfinite(low):v=high;direction=-1;mark='v'
                else:
                    # Two-sided or unbounded slopes are never invented points.
                    # Paired open triangles show the two interval directions.
                    if np.isfinite(low) and np.isfinite(high):
                        ax.vlines(r.amplitude_fraction,low,high,color=col,lw=.6)
                        endpoints=((low,1,'^'),(high,-1,'v'))
                    else:endpoints=((-1.,-1,'v'),(1.,1,'^'))
                    for value,sgn,symbol in endpoints:
                        ax.plot(r.amplitude_fraction,value,symbol,mfc='white',mec=col,ms=4)
                        ax.annotate('',xy=(r.amplitude_fraction,value+sgn*.8),xytext=(r.amplitude_fraction,value),arrowprops=dict(arrowstyle='->',color=col,lw=.7,shrinkA=0,shrinkB=0,mutation_scale=5.5))
                    continue
                ax.plot(r.amplitude_fraction,v,mark,mfc='white',mec=col,ms=4)
                if direction:ax.annotate('',xy=(r.amplitude_fraction,v+direction*.8),xytext=(r.amplitude_fraction,v),arrowprops=dict(arrowstyle='->',color=col,lw=.7,shrinkA=0,shrinkB=0,mutation_scale=5.5))
        for y in (2,0,-1):ax.axhline(y,color='#bbbbbb',lw=.6,ls='--',zorder=0)
        ax.set_xscale('log');ax.set_xticks([.3,1,3],labels=['0.3','1','3']);ax.xaxis.set_minor_locator(NullLocator());ax.set(xlabel='Amplitude / reference',ylim=(-4,7),xlim=(.26,3.5))
        if legend:ax.legend(frameon=False,fontsize=5.3,loc='upper left')
    for i,model in enumerate(MODELS):
        ax=axis(eg[i],NAMES[model]);alpha_plot(ax,alpha[(alpha.model==model)&(alpha.arm=='per_perimeter')],dict(zip(LEVELS,ramp(HUES[model],3))),'rate_level')
        if i==0:ax.set_ylabel('Local α')
        else:ax.tick_params(labelleft=False)
    shared_key(outer[2],'Colour: event rate; guides: α = 2, 0, −1; △ + arrow: bound; × pending',LEVELS)
    fg=outer[3].subgridspec(2,3,wspace=.28,hspace=.65);label(outer[3],'F')
    for i,model in enumerate(('wave_pinning','spring_ring')):
        for j,rate in enumerate(LEVELS):
            ax=axis(fg[i,j],f'{NAMES[model]}; {rate} rate');alpha_plot(ax,alpha[(alpha.model==model)&(alpha.rate_level==rate)],dict(zip(('per_perimeter','per_cell'),ramp(HUES[model],2))),'arm')
            if i==0:ax.set_xlabel('')
            if j==0:ax.set_ylabel('Local α')
            else:ax.tick_params(labelleft=False)
    shared_key(outer[3],'Colour: supply arm; guides: α = 2, 0, −1; △ + arrow: bound; × pending',('per perimeter','per cell'))
    gg=outer[4].subgridspec(1,4,wspace=.42);label(outer[4],'G');events=pd.read_csv(DATA/'G_event_heading_change.csv');rate_table=pd.read_csv(DATA/'mode1_rates.csv')
    for i,model in enumerate(MODELS):
        ax=axis(gg[i],NAMES[model]);q=events[(events.model==model)&(events.window_frames==5)]
        for level,color in zip(LEVELS,ramp(HUES[model],3)):
            z=q[q.amplitude_level==level].sort_values('mean_physical_C')
            # The unavailable shared Spring centre remains a visible gap.
            line=z.set_index('size').reindex([1.,2.,3.]) if model=='spring_ring' else z
            ax.plot(line.mean_physical_C,line.mean_heading_change,'o-',ms=2.5,color=color,label=level)
            if len(z):ax.errorbar(z.mean_physical_C,z.mean_heading_change,yerr=[z.mean_heading_change-z.ci_low,z.ci_high-z.mean_heading_change],fmt='none',color=color,lw=.6,capsize=1.4)
        ax.set_xlabel('Mean circumference C');ax.set_ylim(bottom=0)
        known=rate_table[rate_table.model==model].mean_physical_C
        if len(known):
            low,high=float(known.min()),float(known.max());pad=.08*(high-low);ax.set_xlim(low-pad,high+pad)
        ax.ticklabel_format(axis='x',useOffset=False,style='plain');ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xticks([1,2,3] if model=='spring_ring' else [20,40,60])
        if i==0:ax.set_ylabel('⟨1 − cos Δθ⟩')
        if model=='spring_ring':
            missing=[str(int(s)) for s in (1.,3.) if s not in set(q['size'])]
            note='Per cell\ncentre unavailable'+('\nsize '+', '.join(missing)+' pending' if missing else '')
            ax.text(.55,.95,note,fontsize=5.1,ha='center',va='top',transform=ax.transAxes,bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=1))
    shared_key(outer[4],'Colour: event amplitude; low rate; w = 5 original frames',LEVELS)
    save('EF_analytic_guides.csv',[dict(amplitude_fraction=x,alpha=y,label=f'α = {y}') for y in (2,0,-1) for x in (.3,3.)])
    for ax in fig.axes:
        ticks=ax.get_xticklabels()
        if ticks:ticks[0].set_ha('left');ticks[-1].set_ha('right')
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'figureR2_main.{ext}',dpi=300)
    plt.close(fig);documents()

def documents():
    setup();rates=pd.read_csv(DATA/'mode1_rates.csv') if (DATA/'mode1_rates.csv').exists() else pd.DataFrame();alpha=pd.read_csv(DATA/'EF_local_alpha.csv') if (DATA/'EF_local_alpha.csv').exists() else pd.DataFrame();cov=pd.read_csv(DATA/'coverage.csv') if (DATA/'coverage.csv').exists() else pd.DataFrame();events=pd.read_csv(DATA/'G_event_heading_change.csv') if (DATA/'G_event_heading_change.csv').exists() else pd.DataFrame()
    unexpected=[];comparisons=[]
    if len(cov):
        pending=cov[cov.complete_seeds<4]
        unexpected.append(f'{len(pending)} conditions have incomplete execution; these are marked × pending, not treated as rate censoring. Existing-arm extension checks must finish before this is final inference.')
    if len(alpha):
        bounded=alpha[~alpha.status.isin(['resolved','execution_pending'])];unexpected.append(f'{len(bounded)} α estimates are censored; open triangles/arrows preserve one-sided bounds, and paired opposing open triangles/arrows identify two-sided or unbounded slopes. No horizon caps or omitted bounds.')
        q=alpha[(alpha.arm=='per_perimeter')&(alpha.rate_level=='low')&(alpha.amplitude_level=='low')]
        expected=dict(goryachev=4.,wave_pinning=1.6,holmes_model3=.9,spring_ring=0.)
        for r in q.itertuples(index=False):unexpected.append(f'Pilot check {NAMES[r.model]}: {estimate_text(r)}, status {r.status}; earlier check '+('above 4' if r.model=='goryachev' else f'about {expected[r.model]:g}')+f'. Bootstrap resolved fraction {getattr(r,"bootstrap_resolved_fraction",float("nan")):.4g}. No parameter or seed was changed.')
        for model in ('wave_pinning','spring_ring'):
            pp=alpha[(alpha.model==model)&(alpha.arm=='per_perimeter')&(alpha.status=='resolved')]
            pc=alpha[(alpha.model==model)&(alpha.arm=='per_cell')&(alpha.status=='resolved')]
            matched=pp.merge(pc,on=['model','rate_level','amplitude_level'],suffixes=('_perimeter','_cell'))
            if len(matched):
                delta=matched.alpha_cell-matched.alpha_perimeter
                comparisons.append(f'{NAMES[model]} per-cell minus per-perimeter α ranges from {delta.min():.4g} to {delta.max():.4g} across {len(matched)} matched resolved conditions (earlier check about +1)')
            else:comparisons.append(f'{NAMES[model]} per-cell versus per-perimeter α comparison remains pending')
        unexpected.extend(comparisons)
    unexpected+=['Existing Spring event times depend on model state through its changing physical circumference and cannot be regenerated exactly from seeds alone. Panel G therefore uses the new per-cell Spring outer sizes only, with the shared centre unavailable; chemical panels use regenerated existing-arm events.',
      'Spring fixed per-cell intensity is ρ times the exact basal physical circumference at size 2. The shared centre keeps its original state-dependent per-perimeter hazard, as required by the no-rerun instruction; it is not an exactly constant-intensity centre trajectory.',
      'The β = −1 sample paths retain the original round-7 literal Bernoulli arrival transition at interval 20. Panel B is the requested continuous-time analytic exponential, so this tiny finite-step discrepancy is explicit.',
      'Cadence validation used the 46 completed target-model cohorts available before launch, including the shared centre and every rate/amplitude combination at C20. Eight high-rate target-model cohorts were still executing and were not used to choose cadence.',
      'The finite Fourier cutoff produces small negative ringing in reconstructed Goryachev activity (minimum about −0.423). These band-limited reconstruction values are not negative native solver concentrations. Saved Fourier coefficients, including mode zero M₀, remain unchanged, and no reconstruction values are clipped. Its shared C20/C60 grayscale uses an asinh transform with linear width 1 native activity unit, explicitly keyed; neither ring is normalized separately.']
    cadence_postlaunch=META/'record_cadence_postlaunch_summary.json'
    if cadence_postlaunch.exists():
        audit=C.read(cadence_postlaunch)
        preservation='The frozen prelaunch selection and validation records remain unchanged.' if audit.get('frozen_selection_unchanged') else 'The frozen-selection preservation check failed and requires review.'
        unexpected.append(f'An additive post-launch cadence audit checks {audit["cohorts_checked"]} currently complete existing-arm cohorts at the selected cadences: {audit["passes"]} pass, {len(audit["failures"])} fail, and {len(audit["pending_cohorts"])} cohorts await completion. {preservation} Full comparisons are retained in record_cadence_postlaunch_validation.csv; snapshot identities are beside this readme.')
    cadence_final=META/'record_cadence_finalization.json'
    if cadence_final.exists() and not C.read(cadence_final)['passed']:
        unexpected.append('The final complete-grid cadence validation failed or is incomplete. The full available figure is retained; the exact audit, exit status and remaining discrepancy are recorded in record_cadence_finalization.json. No cadence or model parameter was changed to close this gap.')
    numbers=['A: λ₀ = 1/7000; duration 24,000; original sample interval 20; physical peak SD 1.5. C20/C60 widths are 1:3.','B: Exact continuous-time λ(C) = λ₀(C/40)^β and C(ℓ) = exp[−λ(C)ℓ]. No null ensembles were run.','C: Existing-arm replicate 0 only, reconstructed from all stored active Fourier modes at 240 spatial positions; 321 displayed times. Goryachev low rate/low amplitude; Wave-pinning low rate/medium amplitude.','D: Whole-profile heading θ = arg P, P = M₁/M₀; endpoint |P| ≥ 0.1; 2000 whole-seed bootstrap draws. Guides use the C20 curve directly (dotted) and its 1/9 power (dashed).']
    if len(alpha):
        for panel,arm in (('E','per_perimeter'),('F','per_cell')):
            q=alpha[alpha.arm==arm];numbers.append(f'{panel}: {len(q)} model/rate/amplitude combinations; {sum(q.status=="resolved")} resolved, {sum(q.status=="execution_pending")} pending, {sum(~q.status.isin(["resolved","execution_pending"]))} bounded.')
    if len(events):
        for model in MODELS:
            q=events[(events.model==model)&(events.window_frames==5)];numbers.append(f'G {NAMES[model]}: {int(q.scored_events.sum()) if len(q) else 0} scored event-window pairs across displayed amplitude/size conditions.')
    methods=['All times, lengths, rates and activities use native units. The models retain their distinct native spatial and temporal scales. Biological migration persistence, conditional OU estimates and these 1D activity-heading readouts are distinct.',
      'The existing Fourier rate estimator is imported unchanged from the frozen Poisson campaign: uncentered actual-pair-count mode autocorrelation, raw numerators and powers pooled over four seeds, first linearly interpolated half-height crossing. The record extends together for all four seeds when its half-height lag exceeds record/8, through the original fourfold maximum. Original upper/lower rate statuses are retained.',
      'Mode 1 is sufficient for γ₁ by definition; active mode 0 is retained to compute P. The validation compares cadence strides 1,2,4,8,16,32,64,128 against unchanged full-cadence rates and exact four-seed bootstrap CIs. Selected strides are 64 (Wave-pinning, cadence32) and32 (Spring, cadence1.6). P remains at original cadence .5/.05 to score w=2,5,20 original-frame windows. No inactive modes or full profiles are added to the existing campaign.',
      'For E/F, α is minus the OLS slope of log γ₁ against log actual mean physical circumference over the three sizes. 2000 whole-seed bootstrap draws independently resample the four seeds at each size. The unchanged mode estimator is applied to each resampled pooled numerator. A censored γ gives its mathematical log-rate interval; signed OLS coefficients propagate it to a one-sided or unbounded α interval. Bootstrap censoring is recorded rather than replaced by a finite rate.',
      'G compares circular mean headings in [t−w,t) and (t,t+w] around each forcing event using w original stored frames, excluding any sample exactly at the event. Every finite nonzero circular mean enters, without a coherence gate. The reported value is pooled over scorable events; confidence intervals resample whole seeds. Overlapping event windows are retained. Chemical event regeneration uses the three frozen independent SeedSequence streams and checks total event counts against continuation checkpoints.',
      'R2 kymographs use the first predeclared replicate with no seed search. Fourier reconstruction uses the retained half-cell phase convention directly: the nonnegative-index active coefficients and their conjugate partners reconstruct the spatial signal. Finite retained mode cutoff makes this a band-limited reconstruction, never a rerun. The stored coefficients and their mode-zero contribution are unchanged; small negative reconstructed values are Fourier ringing rather than native concentration negativity. Goryachev uses one shared asinh activity scale with linear width 1; all negative ringing is retained. Nulls and Wave-pinning use linear grayscale. D displays the common unrescaled lag range 0–12,000, with full curves retained in CSV. Its fractional-power guide is undefined at nonpositive measured C20 and those guide values remain missing, with no fitting or smoothing.',
      'Rebuild from compact records: python polarity1d/manuscript/code/mechanisms/rebuild_figureR2_20261027.py render. Re-extract immutable HPC evidence under SLURM with the same entry point and command all. Campaign manifests, current producer state and code SHA are in provenance.json.']
    generated=C.read(META/'provenance.json').get('generated_utc','unavailable') if (META/'provenance.json').exists() else 'unavailable'
    numbers.insert(0,f'Data generation time: {generated}. Current campaign and extension states are retained in provenance.json; each pending condition is explicitly listed in coverage.csv.')
    text='# Figure R2 scientific readme\n\n## Unexpected\n\n'+'\n\n'.join(unexpected)+'\n\n## Key numbers by panel\n\n'+'\n\n'.join(numbers)+'\n\n## Methods\n\n'+'\n\n'.join(methods)+'\n\n## Event-window sensitivity (readme only)\n\n'
    if len(events):
        text+='| Model | Arm | Size | Amplitude | w=2 | w=5 | w=20 |\n|---|---|---:|---|---:|---:|---:|\n'
        for key,g in events.groupby(['model','arm','size','amplitude_level'],sort=False):
            values={int(r.window_frames):r.mean_heading_change for r in g.itertuples(index=False)}
            text+='| '+' | '.join([NAMES[key[0]],key[1],str(key[2]),key[3]]+[f'{values.get(w,float("nan")):.6g}' for w in (2,5,20)])+' |\n'
    (META/'scientific_readme.md').write_text(text)
    result='Execution and extension checks are still running; pending values are marked and recorded.'
    if len(alpha):
        q=alpha[(alpha.arm=='per_perimeter')&(alpha.rate_level=='low')&(alpha.amplitude_level=='low')]
        result='; '.join(f'{NAMES[r.model]} {estimate_text(r)} ({r.status})' for r in q.itertuples(index=False))+'.'
    legend=['# Jump supply, event response and size-dependent heading memory',
      '**A.** Jump-null kymographs at C20 and C60, with width proportional to circumference, grayscale activity and time upward, faceted by β. Uniform peak relocations arrive at λ(C)=λ₀(C/40)^β, λ₀=1/7000. β=1 denotes supply per perimeter and constant success; β=0 denotes perimeter supply with success proportional to 1/C, or per-cell supply with constant success; β=−1 denotes perimeter supply with success proportional to 1/C². Quantities use native units. The fixed paths show different size dependence of relocation frequency without changing peak width.',
      '**B.** Closed-form heading autocorrelation C(ℓ)=exp[−λ(C)ℓ] for the same nulls. Light and dark gray indicate C20 and C60. Larger rings lose heading faster at β=1, equally at β=0 and slower at β=−1.',
      '**C.** Event-noise kymographs reconstructed from the existing per-perimeter campaign’s stored Fourier modes, using fixed replicate 0: Goryachev at low event rate and low amplitude, and Wave-pinning at low rate and medium amplitude. Each C20/C60 pair shares its activity scale and common sampled time axis; widths follow circumference. Goryachev uses a shared asinh grayscale with linear width 1 activity unit, preserving small negative Fourier ringing; Wave-pinning uses a linear scale. Goryachev retains a stationary C60 peak while C20 relocates; Wave-pinning exhibits moving and deforming activity at both sizes.',
      '**D.** Heading autocorrelation of the same conditions, pooled over four seeds with |P|≥0.1 at both endpoints and 95% whole-seed bootstrap bands. P=M₁/M₀ defines the whole-profile heading vector; light/dark model-hue shades indicate C20/C60. Dotted and dashed gray guides are the C20 curve and its 1/9 power, respectively. Goryachev C60 retains near-unit correlation while C20 decorrelates; both Wave-pinning sizes approach zero across the displayed early-lag window, below the physical-drift guide.',
      '**E.** Local size exponent α against mean event amplitude/reference, faceted by model for the per-perimeter arm. α is minus the three-size OLS slope of log mode-1 decorrelation rate γ₁ against log actual mean circumference; γ₁ is ln2 divided by the first mode-1 autocorrelation half-height lag. Colour intensity encodes event rate (light low, dark high), and guides are α=2,0,−1. Filled points have resolved rates and 95% whole-seed bootstrap intervals; open triangles and arrows show bounds; paired opposing open triangles/arrows denote two-sided or unbounded slopes; × means execution pending. '+result,
      '**F.** The same exponent against event amplitude compares per-perimeter (light) and per-cell (dark) supply, faceted by Wave-pinning/Spring and event rate; guides and censor symbols match E. Per-cell arrivals use the fixed C40 rate, with shared existing C40 trajectories; for Spring this uses the exact size-2 basal physical circumference and retains the original dynamic-hazard shared centre. '+('; '.join(comparisons)+'.' if comparisons else 'The supply-arm comparison remains pending.'),
      '**G.** Mean event-associated heading change ⟨1−cosΔθ⟩ against mean physical circumference at low rate, with amplitude encoded light to dark. Δθ compares circular mean headings in the five-original-frame windows immediately before and after each forcing event; bars are 95% whole-seed bootstrap intervals. Chemical models use exact seed-regenerated existing-arm event logs. Spring uses new per-cell outer sizes only because existing Spring event times depend on its evolving circumference; its shared centre is unavailable. Window sensitivities at two and twenty frames appear in the readme. At high amplitude the event-associated change decreases from C20 to C60 in all three chemical models.']
    (META/'main_legend.md').write_text('\n\n'.join(legend)+'\n')
    (OUT/'README.md').write_text('# Figure R2\n\n![Figure R2](figureR2_main.png)\n\n[PDF](figureR2_main.pdf) · [SVG](figureR2_main.svg) · [Legend](supporting/metadata/main_legend.md) · [Scientific readme](supporting/metadata/scientific_readme.md) · [Provenance](supporting/metadata/provenance.json)\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=('extract','render','all'));a=ap.parse_args()
    if a.command in ('extract','all'):extract()
    if a.command in ('render','all'):render()
