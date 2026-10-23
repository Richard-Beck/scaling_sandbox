"""Bounded existing-profile Fourier and exact own-mask region checks."""
from common import *
import argparse,subprocess,shlex,itertools,time
import numpy as np
from mode_analysis import fourier,sufficient,aggregate_modes
import signature_figure4_regions_20261020 as regions
import figure4_round5_analysis_20261023 as A

def exact_decomposition(t,p,C,basal,lag):
    """Matched native-mass denominator; exact moments include each region shape.

    Integrate excess activity on the region's full basin mask using saved cell
    quadrature. Shared-region motion changes phase; the coexisting factor
    changes complex magnitude (amount and shape), membership changes IDs.
    Therefore within-region shape is allocated explicitly instead of discarded.
    An additional contour-mask sensitivity is also evaluated.
    """
    original=regions.frame_regions;captured=[]
    def capture(profile,circumference,*args,**kw):
        result=original(profile,circumference,*args,**kw);captured.append(result);return result
    regions.frame_regions=capture
    try:tracked=regions.track_profiles(t,p,C,basal_contrast=float(np.ptp(basal)),turning_intervals=())
    finally:regions.frame_regions=original
    n=p.shape[1];phase=np.exp(2j*np.pi*(np.arange(n)+.5)/n);mean=p.mean(1)
    P=p@phase/n/mean
    maps={k:[] for k in ('exact_full_mask','exact_contour_mask','point_matched_mass')};weights=[]
    for frame,(profile,circ,desc) in enumerate(zip(p,C,captured)):
        by={k:{} for k in maps};ww={};floor=profile.min();dx=circ/n;M0=mean[frame]*circ
        for reg in desc['regions']:
            identity=reg['track_id'];excess=profile-floor
            for key,mask in [('exact_full_mask',reg['full_mask']),('exact_contour_mask',reg['mask'])]:
                by[key][identity]=complex(np.sum(excess[mask]*phase[mask])*dx/M0)
            by['point_matched_mass'][identity]=complex(np.sum(excess[reg['full_mask']])*dx/M0*np.exp(1j*reg['heading'])) if reg['heading'] is not None else complex(np.nan)
            ww[identity]=float(np.sum(excess[reg['full_mask']])*dx/M0)
        for key in maps:maps[key].append(by[key])
        weights.append(ww)
    maps['round5_point_native_amplitude']=maps['point_matched_mass']
    maps['round5_exact_native_amplitude']=maps['exact_full_mask']
    step0=np.arange(len(t));step1=np.searchsorted(t,t+lag)
    ok=step1<len(t);step0=step0[ok];step1=step1[ok];ok=np.isclose(t[step1]-t[step0],lag,atol=1e-5)
    step0=step0[ok];step1=step1[ok]
    outputs=[]
    for convention,frames in maps.items():
        sums=np.zeros(4);total=0.;eligible=0;maxerr=0.;native_amplitude=convention.startswith('round5_')
        for i,j in zip(step0,step1):
            if min(abs(P[i]),abs(P[j]))<.1:continue
            a,b=frames[i],frames[j]
            if not a or not b or any(not np.isfinite(z) for z in list(a.values())+list(b.values())):continue
            common=set(a)&set(b);lost=sum(a[k] for k in set(a)-common);born=sum(b[k] for k in set(b)-common)
            states={};native_amplitude=convention.startswith('round5_')
            for m,w,membership in itertools.product((0,1),repeat=3):
                if native_amplitude:
                    w0,w1=weights[i],weights[j]
                    z=sum((w1[k] if w else w0[k])*((b[k]/w1[k]) if m else (a[k]/w0[k])) for k in common)
                    totalW=sum(w1[k] if w else w0[k] for k in common)+sum((w1 if membership else w0)[k] for k in (set(b) if membership else set(a))-common)
                    states[m,w,membership]=(z+(born if membership else lost))/totalW if totalW>1e-12 else complex(np.nan)
                else:
                    z=sum(abs(b[k] if w else a[k])*np.exp(1j*np.angle(b[k] if m else a[k])) for k in common)
                    states[m,w,membership]=z+(born if membership else lost)
            if any(not np.isfinite(v) for v in states.values()):continue
            components=[]
            for factor in range(3):
                z=0j;others=[k for k in range(3) if k!=factor]
                for x,y in itertools.product((0,1),repeat=2):
                    s=[0,0,0];s[others[0]]=x;s[others[1]]=y;before=tuple(s);s[factor]=1
                    z+=(1/3 if x==y else 1/6)*(states[tuple(s)]-states[before])
                components.append(z)
            if native_amplitude:components=[z*.5*(abs(P[i])+abs(P[j])) for z in components]
            components.append(P[j]-P[i]-sum(components))
            h=P[j]/abs(P[j])-P[i]/abs(P[i]);loss=np.real(np.array(components)*np.conj(h))/(abs(P[j])+abs(P[i]))
            actual=.5*abs(h)**2;maxerr=max(maxerr,abs(sum(loss)-actual));sums+=loss;total+=actual;eligible+=1
        outputs.append(dict(convention=convention,eligible_pairs=eligible,total_heading_loss=total,
            motion_fraction=sums[0]/total if total else np.nan,amount_shape_fraction=sums[1]/total if total else np.nan,
            membership_fraction=sums[2]/total if total else np.nan,remainder_fraction=sums[3]/total if total else np.nan,
            max_closure_error=maxerr,reconstructed_P_error_rms=float(np.sqrt(np.mean([abs(P[f]-sum(z.values()))**2 for f,z in enumerate(frames)]))),
            motion_definition='shared phase and within-region shape; native-amplitude convention' if native_amplitude else 'shared phase; within-region magnitude shape allocated to amount/shape'))
    # Reproduce the round5 convention on precisely these bounded profiles and
    # endpoints, allowing comparison with a changed denominator to remain clear.
    frames=tracked['frames'];rr=[r for f in frames for r in f['regions']]
    z=dict(time=t,M0=mean*C,M1_real=(P*mean*C).real,M1_imag=(P*mean*C).imag,
      region_offset=np.r_[0,np.cumsum([len(f['regions']) for f in frames])],region_id=np.array([r['track_id'] for r in rr]),
      region_amount=np.array([r['excess'] for r in rr]),region_heading=np.array([r['heading'] if r['heading'] is not None else np.nan for r in rr]),
      region_amount_available=np.array([r['amount_available'] for r in rr]))
    # Native step grid is nonuniform in saved exemplars; evaluate the inherited
    # rule only for regular endpoints below, do not silently assume a median dt.
    return outputs

def worker(ci):
    c=read(SOURCE/'manifest.json')['conditions'][ci];seeds=[];mask=[];phasecheck=[];Cmeans=[]
    oldlag=next(r for r in csvread(OLD/'supporting/data/analysis_lag_metrics.csv') if r['model']==c['model'] and num(r,'C_label')==c['size'] and num(r,'f')==c['fraction'])
    lag=num(oldlag,'analysis_lag');started=time.monotonic()
    with np.load(SOURCE/c['input']) as z:basal=z['state'][:c['cells']] if c['platform']=='spring' else z['state'][0]
    for rep in (0,1):
        folder=SOURCE/'workers'/f'c{ci:03d}'/f'seed{rep:03d}'
        with np.load(folder/'exemplar.npz') as z:
            t=z['time'].astype(float);p=z['primary'].astype(float);C=z['circumference'].astype(float);indices=z['native_frame_indices']
        a=fourier(p);s=sufficient(t,a,c['model']);seeds.append(s);Cmeans.append(float(C.mean()))
        with np.load(folder/'compact.npz') as z:
            nativeP=z['P_real'][indices]+1j*z['P_imag'][indices]
        phase=np.angle(a[:,1]);nativephase=np.angle(nativeP)
        angleerr=np.angle(np.exp(1j*(phase-nativephase)))
        phasecheck.append(dict(condition_id=c['condition_id'],model=c['model'],C_label=c['size'],f=c['fraction'],replicate=rep,
            heading_phase_RMS_radians=np.sqrt(np.mean(angleerr**2)),heading_phase_max_radians=np.max(abs(angleerr)),
            source=str(folder/'exemplar.npz'),source_sha256=digest(folder/'exemplar.npz'),native_cells=c['cells'],saved_bins=p.shape[1],frames=len(t)))
        # Exact same endpoint eligibility compares binned-profile and native
        # phase correlations; bin/float16 error is reported instead of hidden.
        native_s=sufficient(t,np.column_stack([np.ones(len(t)),nativeP]),c['model'],abs(nativeP))
        for j,l in enumerate(s['lags']):
            x=s['phase_num'][j]/s['phase_den'][j] if s['phase_den'][j] else np.nan
            y=native_s['phase_num'][j]/native_s['phase_den'][j] if native_s['phase_den'][j] else np.nan
            phasecheck.append(dict(condition_id=c['condition_id'],model=c['model'],C_label=c['size'],f=c['fraction'],replicate=rep,
                 lag=l,phase_correlation=x,native_heading_correlation=y,difference=x-y,source=str(folder/'exemplar.npz')))
        if np.isfinite(lag):
            for r in exact_decomposition(t,p,C,basal,lag):mask.append(dict(condition_id=c['condition_id'],model=c['model'],C_label=c['size'],f=c['fraction'],replicate=rep,analysis_lag=lag,**r))
    c=dict(c,C_mean=float(np.mean(Cmeans)));rates,curves=aggregate_modes(seeds,c,identity=('exemplar',ci))
    dest=OUT/'exemplar_analysis'/f'c{ci:03d}';dest.mkdir(parents=True,exist_ok=True)
    for name,rows in [('rates',rates),('curves',curves),('exact_masks',mask),('phase_check',phasecheck)]:csvwrite(dest/(name+'.csv'),rows)
    write(dest/'complete.json',dict(condition=c,wall_seconds=time.monotonic()-started,
         algorithm_sha256=digest(Path(__file__)),mode_algorithm_sha256=digest(Path(__file__).with_name('mode_analysis.py'))))
    print(c['condition_id'],time.monotonic()-started,flush=True)

def collect():
    for name in ('rates','curves','exact_masks','phase_check'):
        rows=[]
        for p in sorted((OUT/'exemplar_analysis').glob('c*/'+name+'.csv')):rows.extend(csvread(p))
        csvwrite(DATA/('exemplar_'+name+'.csv'),rows)
    write(META/'exemplar_status.json',dict(complete_conditions=len(list((OUT/'exemplar_analysis').glob('c*/complete.json'))),designed_conditions=260))

def submit():
    logs=OUT/'logs';logs.mkdir(exist_ok=True)
    cmd=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=2G','--time=00:30:00','--array=0-259',
        '--job-name=F4r6-existing',f'--output={logs}/%A_%a.out',f'--error={logs}/%A_%a.err',
        '--wrap',shlex.join([PYTHON,'-B',str(Path(__file__).resolve()),'worker','--condition'])+' "$SLURM_ARRAY_TASK_ID"']
    p=subprocess.run(cmd,text=True,capture_output=True);write(META/'exemplar_submission.json',dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,
         qos_reason='260 independent bounded-profile analyses, burst throughput; one CPU per condition, no throttle'))
    return p.stdout,p.stderr

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['worker','collect','submit']);p.add_argument('--condition',type=int);a=p.parse_args()
    if a.action=='worker':worker(a.condition)
    else:print(globals()[a.action]())
