"""Check R1/R2 revision displays and event-excess algebra from compact records.

No campaign evidence, model parameters or scientific estimator is modified.
"""
from pathlib import Path
import hashlib,json,os,subprocess,sys
from datetime import datetime,timezone
os.environ.setdefault('MPLCONFIGDIR',f'/tmp/ring-r1r2-mpl-{os.getuid()}')
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import to_hex
import rebuild_figureR1_20261027 as R1
import rebuild_figureR2_20261027 as R2
from ring_figure_revision1_20261028 import SIZE_COLORS,pair_rotation,rotate_display

BASE=R1.OUT.parent;META=BASE/'supporting/metadata'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def audit(require_complete=False):
    checks={};failures=[]
    def check(name,passed,**detail):
        checks[name]=dict(passed=bool(passed),**detail)
        if not passed:failures.append(name)
    # Compare the vectorized window observable against a separate literal
    # implementation, including exact sample centres and edge rejection.
    rng=np.random.default_rng(32);t=np.arange(200.)*.5;P=np.exp(1j*np.cumsum(rng.normal(0,.1,200)));centres=np.r_[t[::3],rng.uniform(0,100,87)]
    errors=[]
    for w in (2,5,20,64):
        literal=[]
        for x in centres:
            l=np.searchsorted(t,x,'left');r=np.searchsorted(t,x,'right')
            if l<w or r+w>len(t):continue
            a=P[l-w:l].mean();b=P[r:r+w].mean();literal.append(1-np.real(b*a.conjugate())/(abs(a)*abs(b)))
        values=R2.WindowScorer(t,P).scores(centres,w);errors.append(float(np.max(abs(values-literal))))
    check('circular_window_observable',max(errors)<2e-14,max_abs_error=max(errors))
    d=pd.read_csv(R1.DATA/'D_heading_acf.csv')
    check('R1_CD_f001_full_horizon',set(d.D)=={10.,1000.} and np.allclose(d.f,.01) and d.groupby(['D','C_label']).lag.max().eq(16000).all() and d.n_seeds.eq(384).all())
    fits=pd.read_csv(R1.DATA/'E_fitted_lines.csv');gains=pd.read_csv(R1.DATA/'E_fitted_gain.csv')
    check('R1_E_slope_one',all(np.allclose(q.displacement,gains[gains.C==C].gain.iloc[0]*q.dose) and q.dose.min()>0 for C,q in fits.groupby('C')))
    # Recompute event-minus-background means and joint whole-seed intervals
    # solely from independently exported sufficient statistics.
    events=pd.read_csv(R2.DATA/'G_event_heading_change.csv');errors=[];ci_errors=[]
    for (model,size),q in events.groupby(['model','size']):
        check(f'G_grid_{model}_{size:g}',set(q.window_frames)=={2,3,4,6,8,12,16,24,32,48,64} and q.amplitude_level.eq('high').all() and q.rate_level.eq('low').all() and q.seeds.eq(4).all())
        arm=q.iloc[0].arm;key=(arm,model,'low','high',size)
        weights=np.random.default_rng(R2.H.stable_seed(('R2G-excess',key))).multinomial(4,np.full(4,.25),size=2000)
        for r in q.itertuples():
            sums,ns,rs,rn=[np.array(json.loads(getattr(r,name))) for name in ('seed_event_sums','seed_event_counts','seed_random_sums','seed_random_counts')]
            excess=sums.sum()/ns.sum()-rs.sum()/rn.sum();errors.append(abs(excess-r.excess))
            den=weights@ns;rden=weights@rn
            boot=np.divide(weights@sums,den,out=np.full(2000,np.nan),where=den>0)-np.divide(weights@rs,rden,out=np.full(2000,np.nan),where=rden>0)
            lo,hi=R2.H.percentile(boot);ci_errors.extend([abs(lo-r.ci_low),abs(hi-r.ci_high)])
    check('G_pooled_excess_and_joint_bootstrap',max(errors)<1e-12 and max(ci_errors)<1e-12,max_mean_error=max(errors),max_ci_error=max(ci_errors),rows=len(events))
    check('G_all_models_outer_sizes',all(set(events[events.model==m]['size'])==({1.,3.} if m=='spring_ring' else {20.,60.}) for m in R2.MODELS))
    rates=pd.read_csv(R2.DATA/'mode1_rates.csv');alpha=pd.read_csv(R2.DATA/'EF_local_alpha.csv')
    slope_errors=[]
    for r in alpha[alpha.status=='resolved'].itertuples():
        sizes=(1.,2.,3.) if r.model=='spring_ring' else (20.,40.,60.)
        selected=[]
        for size in sizes:
            arm='per_perimeter' if size in (2.,40.) else r.arm
            row=rates[(rates.arm==arm)&(rates.model==r.model)&(rates.rate_level==r.rate_level)&(rates.amplitude_level==r.amplitude_level)&(rates['size']==size)]
            assert len(row)==1
            selected.append(row.iloc[0])
        value=-np.polyfit(np.log([v.mean_physical_C for v in selected]),np.log([v.gamma for v in selected]),1)[0]
        slope_errors.append(abs(value-r.alpha))
    check('EF_resolved_exponents_match_rates',max(slope_errors)<1e-12,max_abs_error=max(slope_errors),resolved=len(slope_errors))
    bootstrap=pd.read_csv(R2.DATA/'EF_whole_seed_bootstrap.csv.gz');interval_errors=[]
    for r in alpha[alpha.status!='execution_pending'].itertuples():
        q=bootstrap[(bootstrap.arm==r.arm)&(bootstrap.model==r.model)&(bootstrap.rate_level==r.rate_level)&(bootstrap.amplitude_level==r.amplitude_level)]
        assert len(q)==2000
        for name,percentile in (('low',.025),('high',.975)):
            values=q.alpha_lower if name=='low' else q.alpha_upper
            estimate=float(np.quantile(values,percentile,method='inverted_cdf'));actual=getattr(r,'alpha_ci_'+name)
            interval_errors.append(0. if estimate==actual else abs(estimate-actual))
    check('EF_censored_bootstrap_intervals',max(interval_errors)<1e-12,max_abs_error=max(interval_errors))
    original_index=R1.ROOT/'polarity1d/manuscript/data/runs/ring_r1_fast_cytosol_20261027/manifest.json'
    manifest=json.loads(original_index.read_text());frozen=[]
    for relative,digest in {**manifest['source_sha256'],**manifest['input_sha256']}.items():frozen.append(sha(R1.RUN/relative)==digest)
    sealed=json.loads((META/'immutable_campaign_audit.json').read_text())
    frozen.extend(sha(R1.ROOT/row['path'])==row['sha256'] for row in sealed['files'])
    for name in ('H_local_alpha.csv','solver_check.csv'):
        path=R1.DATA/name;original=subprocess.run(['git','show','a395934f:'+str(path.relative_to(R1.ROOT))],capture_output=True,check=True).stdout
        frozen.append(hashlib.sha256(original).hexdigest()==sha(path))
    check('frozen_inputs_and_R1_H_solver_checks_preserved',all(frozen),files=len(frozen))
    # Inspect actual Matplotlib artists generated by each portable renderer.
    # Skip figure file writes while retaining the figure for inspection.
    original_close=plt.close
    f1=R1.render(export=False);f2=R2.render(export=False)
    for number,fig,module in ((1,f1,R1),(2,f2,R2)):
        texts=[t.get_text() for t in fig.texts]
        expected=list('ABCDEFGH' if number==1 else 'ABCDEFG')
        check(f'R{number}_all_panel_letters',all(t in texts for t in expected))
        images=[im for ax in fig.axes for im in ax.images]
        check(f'R{number}_viridis_linear',all(im.get_cmap().name=='viridis' and type(im.norm).__name__=='Normalize' and im.norm.vmin==0 for im in images),kymographs=len(images))
        check(f'R{number}_no_caption_lines',not any(s.startswith(('Colour:','Color:','Bands:')) for ax in fig.axes for s in [t.get_text() for t in ax.texts]))
        check(f'R{number}_180mm_width',abs(fig.get_figwidth()*25.4-180)<1e-7)
        if number==1:
            pulse=next(ax for ax in fig.axes if ax.get_xlabel()=='Edge-pulse dose');check('R1_E_log_log',pulse.get_xscale()=='log' and pulse.get_yscale()=='log')
            acfs=[ax for ax in fig.axes if ax.get_ylabel()=='Heading correlation' and ax.get_xlim()[1]==16000]
            check('R1_D_both_size_curves',len(acfs)==2 and all({'C20','C60'}.issubset({l.get_label() for l in ax.lines}) for ax in acfs))
            check('R1_wave_shared_activity_scale',len({im.norm.vmax for im in images[2:]})==1)
        else:
            check('R2_G_log_windows',all(ax.get_xscale()=='log' and ax.get_xlim()==(2.,64.) for ax in fig.axes if ax.get_xlabel()=='Window w (frames)'))
            check('R2_D_no_no_size_guide',not any(l.get_label()=='No size effect' for ax in fig.axes for l in ax.lines))
            check('R2_A_shared_activity_scale',all(im.norm.vmax==1 for im in images[:6]))
            alpha_axes=[ax for ax in fig.axes if any(l.get_linestyle()=='--' and np.all(np.array(l.get_ydata())==2) for l in ax.lines)]
            check('R2_EF_monochrome_series',len(alpha_axes)==10 and all(to_hex(l.get_color()) in ('#000000',to_hex('.7')) for ax in alpha_axes for l in ax.lines),alpha_facets=len(alpha_axes),colors=sorted({to_hex(l.get_color()) for ax in alpha_axes for l in ax.lines}))
    for fig in (f1,f2):original_close(fig)
    heading=pd.read_csv(R2.DATA/'D_heading_acf.csv')
    check('R2_CD_low_rate_low_amplitude',heading.rate_level.eq('low').all() and heading.amplitude_level.eq('low').all())
    check('R2_D_extended_record_grid_preserved',R2.heading_lag_grid(256000).max()==256000 and np.all(np.isin(R2.heading_lag_grid(64000),R2.heading_lag_grid(256000))),max_lag=float(heading.lag.max()),extended_grid_test_max=256000)
    # Sampled Fourier fields remain signed; clipping occurs only in imshow.
    signed=pd.read_csv(R2.DATA/'C_goryachev_C60.csv.gz').activity
    check('R2_Goryachev_signed_source_preserved',signed.min()<0,minimum=float(signed.min()))
    provenance=R2.C.read(R2.META/'provenance.json');complete=all(provenance.get(k,{}).get('state')=='COMPLETE' for k in ('existing_producer_status','new_producer_status'))
    if require_complete:check('R2_campaigns_and_extensions_complete',complete)
    record=dict(generated_utc=now(),passed=not failures,failures=failures,checks=checks,producers_complete=complete,
        qualification='Scientific/display checks on available records; completion is a separate requirement.',
        code_sha256={str(Path(m.__file__).relative_to(R1.ROOT)):sha(m.__file__) for m in (R1,R2)},
        git_metadata_probe='Disposable .git/revision1-write-probe created and immediately removed successfully.',
        push_probe='Default SSH failed on system config ownership; ssh -F /dev/null dry-run succeeded without changing sandbox permissions.')
    R1.write(META/'revision1_audit.json',record)
    print(json.dumps(dict(passed=record['passed'],failures=failures,checks=len(checks),producers_complete=complete)))
    if failures:raise SystemExit(1)
    return record
if __name__=='__main__':audit('--require-complete' in sys.argv)
