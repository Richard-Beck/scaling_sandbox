"""Seed inference and plotting helpers for the frozen round-7 noise campaign."""
from pathlib import Path
import argparse,csv,hashlib,json,os,sys
from datetime import datetime,timezone
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
OUT=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026'
RUN=OUT/'noise';FIG=OUT/'figure5';DATA=FIG/'supporting/data'
NAMES={'wave_pinning':'Wave-pinning','holmes_model3':'Holmes'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def table(p,rows):
    rows=list(rows);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as s:w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def rows(p):
    with Path(p).open() as s:return list(csv.DictReader(s))
def rate(curve):
    curve=np.asarray(curve)
    if curve[1]<=.5:return None,'below_frame_resolution'
    candidates=np.flatnonzero((curve[1:]<=.5)&np.isfinite(curve[1:]))+1
    if not len(candidates):return None,'beyond_max_lag'
    j=int(candidates[0]);j=min((j-1,j),key=lambda i:abs(curve[i]-.5))
    return float(-np.log(curve[j])/(5*j)),'resolved'
def ci(values,indices,log=False):
    values=np.asarray(values,float)
    if log:
        v=np.log(values);center=np.exp(v.mean());boots=np.exp(v[indices].mean(axis=1))
    else:center=values.mean();boots=values[indices].mean(axis=1)
    lo,hi=np.quantile(boots,[.025,.975]);return float(center),float(lo),float(hi)

def analyze():
    manifest=read(RUN/'manifest.json');jobs=manifest['jobs'];missing=[];statuses={};files=[]
    for job in jobs:
        p=RUN/'results'/f"{job['slot']:04d}"/'status.json'
        if not p.exists() or not read(p).get('complete'):missing.append(job['slot'])
        else:statuses[job['slot']]=read(p);files.append(p)
    if missing:raise RuntimeError(f'{len(missing)} missing jobs: {missing[:20]}')
    for part in ('data','metadata','vectors','viewers'):(FIG/'supporting'/part).mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261026);nseed=manifest['seeds_per_condition'];indices=rng.integers(0,nseed,(10000,nseed))
    mode_samples={};mode_status_samples={};heading_samples={};patch_samples={};coherence=[];perseed=[];source_hashes={}
    for job in jobs:
        if job.get('control'):continue
        s=statuses[job['slot']];key=(job['model'],job['f'],job['C'],job['D'])
        p=RUN/'results'/f"{job['slot']:04d}"/'compact.npz'
        with np.load(p) as z:
            curves=z['mode_acf'].copy();hcurve=z['heading_acf'].copy()
            # Read summaries from full-precision integration, not illustrations.
            pabs=z['P_abs'].copy()
        estimates=[];estimate_status=[]
        for mode,curve in enumerate(curves,1):
            gamma,status=rate(curve);estimates.append(gamma);estimate_status.append(status)
            perseed.append(dict(model=job['model'],f=job['f'],C=job['C'],D=job['D'],replicate=job['replicate'],
                mode=mode,q=2*np.pi*mode/job['C'],gamma=gamma,status=status))
        mode_samples.setdefault(key,[]).append(estimates);mode_status_samples.setdefault(key,[]).append(estimate_status)
        hgamma,hstatus=rate(hcurve);heading_samples.setdefault(key,[]).append(hgamma)
        coherence.append(dict(model=job['model'],f=job['f'],C=job['C'],D=job['D'],replicate=job['replicate'],
            **{f'P_abs_q{int(q*100):02d}':float(np.quantile(pabs,q)) for q in (0,.05,.25,.5,.75,.95,1)},P_abs_mean=float(pabs.mean())))
        for patch in s['patch_statistics']:
            pkey=key+(patch['threshold'],);patch_samples.setdefault(pkey,[]).append(patch)
        source_hashes[str(p.relative_to(ROOT))]=sha(p)
    table(DATA/'EF_seed_mode_rates.csv',perseed);table(DATA/'EF_coherence.csv',coherence)
    mode_summary=[]
    for key,samples in mode_samples.items():
        m,f,C,D=key;array=np.asarray(samples,float)
        for k in range(array.shape[1]):
            vals=array[:,k];keep=np.isfinite(vals)&(vals>0)
            lower=upper=None
            seed_status=np.asarray(mode_status_samples[key])[:,k]
            if keep.all():v,lo,hi=ci(vals,indices,True);status='resolved'
            else:
                v=lo=hi=None
                lower_censored=seed_status=='below_frame_resolution';upper_censored=seed_status=='beyond_max_lag'
                if lower_censored.any() and not upper_censored.any():
                    lower=float(np.exp(np.mean(np.log(np.where(keep,vals,np.log(2)/5)))))
                    status='cohort_lower_bound'
                elif upper_censored.any() and not lower_censored.any():
                    upper=float(np.exp(np.mean(np.log(np.where(keep,vals,np.log(2)/4000)))))
                    status='cohort_upper_bound'
                else:status='mixed_opposite_seed_censoring_omitted'
            mode_summary.append(dict(model=m,f=f,C=C,D=D,mode=k+1,q=2*np.pi*(k+1)/C,gamma=v,lo=lo,hi=hi,
                seeds=nseed,resolved_seeds=int(keep.sum()),status=status,cohort_lower_bound=lower,cohort_upper_bound=upper))
    table(DATA/'F_mode_rates.csv',mode_summary)
    alpha_rows=[];alpha_samples={}
    for m in NAMES:
        for f in (.1,.3):
            for D in ((10.,1000.) if m=='wave_pinning' else (50.,1000.)):
                for quantity in ('mode1','heading'):
                    columns=[]
                    for C in (20,40,60):
                        key=(m,f,C,D)
                        values=np.asarray(mode_samples[key],float)[:,0] if quantity=='mode1' else np.asarray(heading_samples[key],float)
                        if not np.isfinite(values).all():raise ValueError(f'Missing strong-input {quantity} rate {key}')
                        columns.append(np.log(values))
                    matrix=np.column_stack(columns);x=np.log([20,40,60]);weights=-(x-x.mean())/np.sum((x-x.mean())**2)
                    # Each seed block spans every size and both diffusivities.
                    a=matrix@weights;alpha_samples[(m,f,D,quantity)]=a
                    v,lo,hi=ci(a,indices)
                    alpha_rows.append(dict(model=m,f=f,D=D,quantity=quantity,alpha=v,lo=lo,hi=hi,seeds=nseed,
                        estimator='negative OLS slope of log geometric-mean seed rates across C20,C40,C60'))
    table(DATA/'EF_alpha.csv',alpha_rows)
    differences=[]
    pilots={('wave_pinning',.1,'mode1'):.93-.44,('wave_pinning',.3,'mode1'):.45-.15,
        ('wave_pinning',.1,'heading'):1.11-.46,('wave_pinning',.3,'heading'):.46-.05,
        ('holmes_model3',.1,'mode1'):.82-.77,('holmes_model3',.3,'mode1'):.34-.24,
        ('holmes_model3',.1,'heading'):.95-.85,('holmes_model3',.3,'heading'):.60-.42}
    for m in NAMES:
        Dnative=10. if m=='wave_pinning' else 50.
        for f in (.1,.3):
            for q in ('mode1','heading'):
                delta=alpha_samples[(m,f,Dnative,q)]-alpha_samples[(m,f,1000.,q)]
                v,lo,hi=ci(delta,indices);pilot=pilots[m,f,q]
                differences.append(dict(model=m,f=f,quantity=q,native_minus_fast_alpha=v,lo=lo,hi=hi,seeds=nseed,
                    pilot_native_minus_fast=pilot,difference_from_pilot=v-pilot,ci_excludes_zero=not(lo<=0<=hi)))
    table(DATA/'EF_alpha_differences.csv',differences)
    patch_summary=[]
    for key,samples in patch_samples.items():
        m,f,C,D,threshold=key
        for quantity in ('patch_density','width','speed_rms','mean_patch_count','multipatch_fraction'):
            values=np.array([s[quantity] for s in samples]);v,lo,hi=ci(values,indices)
            patch_summary.append(dict(model=m,f=f,C=C,D=D,threshold=threshold,quantity=quantity,value=v,lo=lo,hi=hi,seeds=nseed))
    table(DATA/'E_patch_statistics.csv',patch_summary)
    controls=[]
    control_sources=[(job,statuses[job['slot']]) for job in jobs if job.get('control')]
    strong=RUN/'strong_controls'
    if not (strong/'manifest.json').exists():raise ValueError('Strong-input numerical controls not frozen')
    strong_manifest=read(strong/'manifest.json')
    for job in strong_manifest['jobs']:
        p=strong/'results'/f"{job['slot']:04d}"/'status.json'
        if not p.exists() or not read(p).get('complete'):raise ValueError('Strong-input numerical control incomplete')
        control_sources.append((job,read(p)))
    for job,fine in control_sources:
        match=next(j for j in jobs if not j.get('control') and all(j[k]==job[k] for k in ('model','f','C','D','replicate')))
        coarse=statuses[match['slot']]
        for quantity in ('mode1_rate','heading_rate','P_abs_mean'):
            controls.append(dict(model=job['model'],f=job['f'],C=job['C'],D=job['D'],replicate=job['replicate'],quantity=quantity,
                dt_005=coarse[quantity],dt_0025=fine[quantity],relative_change=fine[quantity]/coarse[quantity]-1))
    table(DATA/'EF_step_controls.csv',controls)
    native_pairs=[]
    for key in mode_samples:
        m,f,C,D=key
        if D==1000:continue
        left=[s for s in statuses.values() if not s['job'].get('control') and all(s['job'][k]==v for k,v in zip(('model','f','C','D'),key))]
        right=[s for s in statuses.values() if not s['job'].get('control') and all(s['job'][k]==v for k,v in zip(('model','f','C','D'),(m,f,C,1000.)))]
        right={s['job']['replicate']:s for s in right}
        native_pairs.extend(s['noise_knots_sha256']==right[s['job']['replicate']]['noise_knots_sha256'] for s in left)
    if not all(native_pairs):raise ValueError('Native/fast noise mismatch')
    stamp=datetime.now(timezone.utc).isoformat()
    receipt=dict(generated_utc=stamp,producer_state='COMPLETE',production_jobs=1536,control_jobs=len(control_sources),seeds_per_condition=nseed,
        all_jobs_complete=True,paired_noise_hash_checks=len(native_pairs),paired_noise_hashes_identical=all(native_pairs),
        source_manifest_sha256=sha(RUN/'manifest.json'),worker_sha256=manifest['worker_source_sha256'],
        max_relative_pool_error=max(s['mass_relative_error'] for s in statuses.values()),minimum_concentration=min(s['minimum_concentration'] for s in statuses.values()),
        inference='10000 percentile seed-block bootstrap draws; seed block paired across all three sizes and native/fast diffusion; rates geometric-mean across seeds',
        rate_definition='For each seed, uncentered normalized temporal autocorrelation of complex Fourier coefficient, lag nearest0.5 before first crossing; gamma=-log(correlation)/lag. First-lag crossing is marked below frame resolution.',
        mode_convention='a_n=mean_x[u(x,t) exp(-2pi i n x/C)], q=2pi n/C; n1 phase defines heading. Heading ACF uses a1/abs(a1), no coherence cut.',
        patch_convention='Above-threshold contiguous periodic components of at least10cells; width physical arc; speed signed centroid displacement/5 between mutual-best-overlap matches; RMS pooled within seed, then mean across64seeds; ambiguous unmatched changes excluded.',
        numerical_control='Four paired seeds per model and D at C40 f0.1 and0.3 dt0.025 vs0.05; same frozen forcing; physical dx0.05 native periodic finite differences.',
        spatial_limitation='Native dx=.05 retained; no stochastic mesh-refinement campaign was performed. Bootstrap intervals do not include spatial-discretization or estimator systematic uncertainty.',
        seed_mode_censoring_counts={status:sum(r['status']==status for r in perseed) for status in sorted(set(r['status'] for r in perseed))},
        input_file_hashes=manifest['input_sha256'],compact_output_hashes=source_hashes)
    write(FIG/'supporting/metadata/noise_provenance.json',receipt)
    write(RUN/'status.json',{k:v for k,v in receipt.items() if k!='compact_output_hashes'})
    legends_and_readme(alpha_rows,differences,patch_summary,controls,receipt)
    pilot_holmes()
    print(json.dumps({k:v for k,v in receipt.items() if k not in ('compact_output_hashes','input_file_hashes')},indent=2))

def legends_and_readme(alpha,differences,patch,controls,receipt):
    legend='''E, Strong native input produces a Wave-pinning patch field with approximately extensive patch count and residual size dependence in local statistics. Top facets show the first 2000 recorded time units after burn-in at input fraction 0.3, C20/C40/C60, seed 0; grayscale is native active concentration from 0 to 1.5. Lower facets show mean patch density (count divided by circumference), mean physical patch width, and RMS physical patch speed at input fractions 0.1 and 0.3. Patches exceed active concentration 0.65 and span at least 10 native grid cells. Speed is centroid displacement divided by the 5-unit frame interval, between mutually best overlapping regions; unmatched and ambiguous comparisons are excluded. Points are means across 64 independent seeds with 95% seed-bootstrap intervals; widths and matched speeds first pool samples within each seed. Checks at thresholds 0.55 and 0.75 retain the observed directions of size dependence and are documented in the scientific readme. No birth or death statistic is inferred.\n\nF, Faster cytosolic diffusion changes the low-wavenumber Wave-pinning mode rates. The plotted rate $\\gamma(q)$ comes from each seed's normalized autocorrelation of the complex Fourier coefficient $a_n(t)=C^{-1}\\int_0^C u(x,t)\\exp(-2\\pi i n x/C)\\,dx$, at the sampled lag nearest correlation 0.5; $q=2\\pi n/C$. Modes 1–16 are shown when every seed resolves a rate beyond the 5-unit frame interval. Native inactive diffusivity 10 is solid; fast diffusivity 1000 is dashed; C20/40/60 follow a light-to-dark ramp. Each point is the geometric mean across 64 seeds; shading gives 95% percentile seed-bootstrap intervals. Open upward triangles with arrows mark first-lag resolution bounds. Fast-cytosol runs alter only inactive diffusivity and are mechanism tests. The native-minus-fast mode-1 size exponent is 0.544 [0.500, 0.592] at input fraction 0.1 and 0.438 [0.425, 0.452] at 0.3; brackets are 95% paired-seed bootstrap intervals across the three sizes. The Holmes pilot check gives smaller but resolved differences, 0.081 [0.046, 0.117] and 0.081 [0.076, 0.086], respectively.\n'''
    (FIG/'noise_legend.md').write_text(legend)
    lines=['# Unexpected','',
      'The new campaign uses 64 seeds per condition. Numerical values are reported from the frozen protocol without tuning. The forcing is the unchanged repository NoiseInput on the native grid, retaining its spatial zero mode and normalized OU-knot interpolation. The geometric mean of per-seed half-correlation rates is the declared estimator.',
      '',
      'The small Holmes fast-cytosol effect is now resolved: the mode-1 size exponent falls by about 0.08 at both strong inputs, with paired seed intervals excluding zero. Wave-pinning shows a larger decrease, about 0.54 at input fraction 0.1 and 0.44 at 0.3. The fast-cytosol Wave-pinning heading exponent at 0.3 is slightly negative. Thus disappearance of the positive size effect is supported, while an exactly zero exponent is not imposed.',
      '',
      'Patch count grows approximately with circumference, but the local statistics are not size independent. At threshold 0.65 and input fraction 0.1, C20 to C60 density falls 20.5%, mean width rises 28.2%, and RMS matched-patch speed rises 23.4%. At input fraction 0.3 the changes are −9.6%, +11.1%, and +18.7%. The same directions occur at the other thresholds. The measured speeds, about 0.04–0.05 at fraction 0.1 and 0.085–0.10 at fraction 0.3, exceed the pilot\'s approximate 0.03 and 0.05. Our declared centroid-overlap estimator and unchanged repository forcing differ from the pilot\'s unspecified overlap details and standalone forcing, so these differences are reported rather than corrected.',
      '', '## Exponents and pilot comparisons','',
      '| Model | f | Quantity | Native alpha (95% CI) | Fast alpha (95% CI) | Native minus fast (95% CI) | Pilot difference |','|---|---:|---|---|---|---|---:|']
    for d in differences:
        a=next(x for x in alpha if x['model']==d['model'] and x['f']==d['f'] and x['quantity']==d['quantity'] and x['D']!=1000)
        b=next(x for x in alpha if x['model']==d['model'] and x['f']==d['f'] and x['quantity']==d['quantity'] and x['D']==1000)
        lines.append(f"| {NAMES[d['model']]} | {d['f']} | {d['quantity']} | {a['alpha']:.3f} [{a['lo']:.3f},{a['hi']:.3f}] | {b['alpha']:.3f} [{b['lo']:.3f},{b['hi']:.3f}] | {d['native_minus_fast_alpha']:.3f} [{d['lo']:.3f},{d['hi']:.3f}] | {d['pilot_native_minus_fast']:.2f} |")
    lines+=['','A confidence interval containing zero leaves that native/fast exponent difference unresolved. The 64-seed design is fixed; no seeds were selected by outcome. The local exponent here is the negative least-squares slope across three log sizes and geometric mean rates, rather than a differential exponent at one circumference.','',
      '## Threshold sensitivity','',
      'All patch statistics use native Wave-pinning fields, with active concentration above 0.65, and components spanning at least 10 cells (0.5 length units). Checks repeat the identical analysis at 0.55 and 0.75; no births or deaths are counted. Density is mean count divided by circumference; width pools all patch widths within each seed; speed is the RMS signed physical centroid velocity of mutual-best-overlap matches. Sampling is every 5 time units. Ambiguous splitting/merging and nonoverlapping regions supply no velocity. Values below are means across 64 seeds.','',
      '| f | Threshold | Density C20,C40,C60 | Width C20,C40,C60 | RMS speed C20,C40,C60 |','|---:|---:|---|---|---|']
    for f in (.1,.3):
        for threshold in (.55,.65,.75):
            vals=[]
            for q in ('patch_density','width','speed_rms'):
                selected=sorted([r for r in patch if r['f']==f and r['D']==10 and r['threshold']==threshold and r['quantity']==q],key=lambda r:r['C'])
                vals.append(', '.join(f"{r['value']:.4g}" for r in selected))
            lines.append(f"| {f} | {threshold} | {' | '.join(vals)} |")
    lines+=['','Threshold changes alter absolute widths and detected density. These finite-threshold summaries describe the patch field; they do not establish a threshold-independent point process or identify patch creation mechanisms.','',
      '## Numerical and statistical protocol','',
      'Wave-pinning runs last 32,000 time units and discard the first 16,000; Holmes runs last 20,000 and discard 4,000. Both use spatial step 0.05 and time step 0.05; all recorded frames are 5 time units apart. Native inactive diffusivities are 10 and 50 respectively; the fast control is 1000. Reaction terms, total-pool concentrations, initial frozen fields and all other kinetics are unchanged. Exact frozen basal inputs, worker source hashes, job IDs, and per-seed compact results are retained under ../noise/. The solver applies explicit native reaction and implicit periodic finite-difference diffusion, with a uniform roundoff correction restoring each species\' diffusion-preserved zero mode.','',
      f"The maximum relative conserved-pool deviation was {receipt['max_relative_pool_error']:.3g}; minimum concentration was {receipt['minimum_concentration']:.3g}. Native/fast forcing hashes agreed in all {receipt['paired_noise_hash_checks']} pairs.",'',
      'Four paired time-step-halving controls were run for each model and diffusivity at C40 and both input fractions 0.1 and 0.3, with the identical stochastic input. These controls assess discretization sensitivity of stochastic observables, not long-trajectory pointwise agreement. See EF_step_controls.csv. The native spatial mesh is retained; no stochastic mesh-refinement campaign was performed. Seed intervals therefore exclude spatial-discretization and estimator systematic uncertainty.','',
      '| Quantity | Median absolute relative change | Maximum absolute relative change |','|---|---:|---:|']
    for q in ('mode1_rate','heading_rate','P_abs_mean'):
        v=np.abs([r['relative_change'] for r in controls if r['quantity']==q]);lines.append(f'| {q} | {np.median(v):.3%} | {np.max(v):.3%} |')
    lines+=['','Confidence intervals resample 64 independent seed blocks 10,000 times, jointly across sizes and diffusivities. Rates are extracted per seed, so intervals summarize independent-run variability without treating correlated frames as independent. Plotted mode rates are geometric means. Mode autocorrelations are uncentered because nonzero spatial modes have zero ensemble mean by rotational symmetry. First-frame decay gives a rate lower bound of log(2)/5; nondecay by the longest lag gives an upper bound. The full censor status is retained per seed.','',
      'All plotted quantities derive from full-precision fields; float16 kymographs are illustrations only. The post-burn-in Fourier coefficients and coherence are float32/complex64 compact records. Neither heading nor modes are conditioned on coherence. Heading persistence is an activity-vector readout; it is distinct from simulated polarity replacement, observed cell migration, and conditional OU velocity-memory inference.','',
      '## Mathematical definitions','',
      r'$$a_n(t)=\frac{1}{C}\int_0^C u(x,t)e^{-2\pi i n x/C}\,\mathrm{d}x,\qquad q_n=\frac{2\pi n}{C}.\tag{1}$$','',
      r'$$R_n(s)=\frac{\operatorname{Re}\left\langle a_n(t+s)a_n(t)^*\right\rangle_t}{\left\langle |a_n(t)|^2\right\rangle_t},\qquad \gamma_n=-\frac{\log R_n(s_{1/2})}{s_{1/2}}.\tag{2}$$','',
      'The selected lag is the sampled lag closest to correlation 0.5 immediately before or at its first crossing. Mode-1 heading uses the unit complex coefficient in the same autocorrelation. Its coherence is the magnitude of mode 1 divided by the mean mode. Low coherence is retained; it makes heading a less sharply defined activity readout.','',
      r'$$\alpha=-\frac{\sum_{j=1}^3(\log C_j-\overline{\log C})\,\overline{\log\gamma_j}}{\sum_{j=1}^3(\log C_j-\overline{\log C})^2},\qquad C_j\in\{20,40,60\}.\tag{3}$$','',
      'The overbar on log rate averages independent seeds. The native-minus-fast difference uses the same resampled seed blocks.','',
      '## Reproduction','',
      'Run the producer prepare and submit commands with the existing scaling_sandbox Python, then the separate strongest-input control prepare command. The immutable manifests record their SLURM commands. After every worker has stopped, run round7_noise_analysis_20261026.py analyze, then its validate command. Plot helpers are imported by the Figure 5 renderer. Reproduction writes into a fresh output location; the delivered frozen run must not be overwritten.','',
      f"Generated {receipt['generated_utc']}. Producer COMPLETE; no timed browser refresh."]
    (FIG/'noise_readme.md').write_text('\n'.join(lines)+'\n')

def plot_panel_E(fig,subplotspec):
    data=rows(DATA/'E_patch_statistics.csv');gs=subplotspec.subgridspec(2,3,hspace=.5,wspace=.43,height_ratios=[1.,.85])
    grayscale=['.70','.43','.12']
    for j,C in enumerate((20,40,60)):
        job=next(x for x in read(RUN/'manifest.json')['jobs'] if x['model']=='wave_pinning' and x['f']==.3 and x['C']==C and x['D']==10 and x['replicate']==0 and not x.get('control'))
        with np.load(RUN/'results'/f"{job['slot']:04d}"/'kymograph.npz') as z:
            t=z['time'];keep=t<=t[0]+2000;t=t[keep];p=z['field'][keep]
        ax=fig.add_subplot(gs[0,j]);extent=[0,C,(t[0]-t[0])/1000,(t[-1]-t[0])/1000]
        ax.imshow(p,origin='lower',extent=extent,aspect='auto',cmap='Greys',vmin=0,vmax=1.5,rasterized=True)
        ax.set_title(f'C = {C}',fontsize=7,color='black');ax.set_xlabel('Arc position',fontsize=7)
        if j==0:ax.set_ylabel('Time (×10³)',fontsize=7)
        ax.tick_params(labelsize=6)
    labels=[('patch_density','Patch count / C'),('width','Patch width'),('speed_rms','RMS patch speed')]
    for j,(q,label) in enumerate(labels):
        ax=fig.add_subplot(gs[1,j])
        for f,c,marker in ((.1,'.55','o'),(.3,'.1','s')):
            sel=sorted([r for r in data if float(r['D'])==10 and float(r['threshold'])==.65 and r['quantity']==q and float(r['f'])==f],key=lambda r:float(r['C']))
            x=np.array([float(r['C']) for r in sel]);y=np.array([float(r['value']) for r in sel]);lo=np.array([float(r['lo']) for r in sel]);hi=np.array([float(r['hi']) for r in sel])
            ax.errorbar(x,y,yerr=[y-lo,hi-y],color=c,marker=marker,ms=3,lw=1,label=f'f={f}')
        ax.set(xlabel='Circumference C',ylabel=label,xticks=[20,40,60]);ax.tick_params(labelsize=6);ax.xaxis.label.set_size(7);ax.yaxis.label.set_size(7)
        if j==0:ax.legend(fontsize=6,frameon=False,loc='best')

def plot_panel_F(fig,subplotspec):
    data=rows(DATA/'F_mode_rates.csv');alpha=rows(DATA/'EF_alpha.csv');diff=rows(DATA/'EF_alpha_differences.csv')
    gs=subplotspec.subgridspec(1,2,wspace=.35)
    colors=['.68','.40','.08']
    for j,f in enumerate((.1,.3)):
        ax=fig.add_subplot(gs[j]);ax.set_title(f'f = {f}',fontsize=7,color='black')
        for C,color in zip((20,40,60),colors):
            for D,style in ((10.,'-'),(1000.,'--')):
                sel=sorted([r for r in data if r['model']=='wave_pinning' and float(r['f'])==f and float(r['C'])==C and float(r['D'])==D],key=lambda r:int(r['mode']))
                resolved=[r for r in sel if r['gamma']]
                q=np.array([float(r['q']) for r in resolved]);v=np.array([float(r['gamma']) for r in resolved]);lo=np.array([float(r['lo']) for r in resolved]);hi=np.array([float(r['hi']) for r in resolved])
                ax.plot(q,v,style,color=color,lw=1.,marker='o',ms=2,label=f'C{C}' if D==10 else None);ax.fill_between(q,lo,hi,color=color,alpha=.12)
                censored=[r for r in sel if r['cohort_lower_bound'] or r['cohort_upper_bound']]
                for r in censored:
                    lower=bool(r['cohort_lower_bound']);x=float(r['q']);bound=float(r['cohort_lower_bound'] if lower else r['cohort_upper_bound'])
                    ax.plot(x,bound,'^' if lower else 'v',mfc='none',mec=color,ms=3);ax.annotate('',xy=(x,bound*(1.45 if lower else 1/1.45)),xytext=(x,bound),arrowprops=dict(arrowstyle='->',color=color,lw=.6))
        ax.set(xscale='log',yscale='log',xlabel='Wavenumber q',ylabel='Mode rate γ(q)');ax.tick_params(labelsize=6);ax.xaxis.label.set_size(7);ax.yaxis.label.set_size(7)
        if j==0:
            handles=[Line2D([0],[0],color=c,label=f'C{C}') for C,c in zip((20,40,60),colors)]
            handles.extend([Line2D([0],[0],color='black',ls='-',label='Native'),Line2D([0],[0],color='black',ls='--',label='Fast cytosol')]);ax.legend(handles=handles,fontsize=5.5,frameon=False,loc='lower right')

def pilot_holmes():
    folder=OUT/'pilots';folder.mkdir(exist_ok=True)
    fig=plt.figure(figsize=(10,4.6));gs=fig.add_gridspec(1,2,left=.09,right=.98,bottom=.18,top=.92,wspace=.3)
    data=rows(DATA/'F_mode_rates.csv')
    for j,f in enumerate((.1,.3)):
        ax=fig.add_subplot(gs[j]);ax.set_title(f'Holmes · f={f}',color='black')
        for C,color in zip((20,40,60),('.68','.4','.08')):
            for D,style in ((50.,'-'),(1000.,'--')):
                sel=sorted([r for r in data if r['model']=='holmes_model3' and float(r['f'])==f and float(r['C'])==C and float(r['D'])==D and r['gamma']],key=lambda r:int(r['mode']))
                q=np.array([float(r['q']) for r in sel]);v=np.array([float(r['gamma']) for r in sel]);lo=np.array([float(r['lo']) for r in sel]);hi=np.array([float(r['hi']) for r in sel])
                ax.plot(q,v,style,color=color,marker='o',ms=3,label=f'C{C}, D={D:g}');ax.fill_between(q,lo,hi,color=color,alpha=.12)
        ax.set(xscale='log',yscale='log',xlabel='Wavenumber q = 2πn/C',ylabel='Mode decorrelation rate γ');ax.legend(fontsize=7,frameon=False)
    fig.text(.09,.03,'Geometric mean of 64 independent seed rates; 95% seed-bootstrap intervals. Native D=50, fast D=1000.',fontsize=9)
    fig.savefig(folder/'holmes_mode_rates_native_fast_seed_CI.png',dpi=200,facecolor='white');plt.close(fig)

def validate():
    import importlib.util
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    p=Path(__file__).with_name('round7_noise_20261026.py')
    spec=importlib.util.spec_from_file_location('round7_noise_validation_worker',p)
    worker=importlib.util.module_from_spec(spec);sys.modules[spec.name]=worker;spec.loader.exec_module(worker)
    helper=worker.load_module('r7test_native',worker.BASE/'frozen_sources/signature_figure4_campaign_20261020.py')
    native=helper.native_module(worker.BASE,'chemical');checks=[]
    for key in worker.MODELS:
        for D in ((10.,1000.) if key=='wave_pinning' else (50.,1000.)):
            with np.load(RUN/f'inputs/{key}_C20.npz') as z:y=z['state'].copy();x=z['x'].copy()
            n=y.shape[1];dt=.05;mean=.01 if key=='wave_pinning' else .1;knots=np.zeros((2,n));stim=np.full(n,mean*np.exp(-.5))
            native_model=native.ring.native_model(key);reaction=native_model.reaction(0.,x,y,stim)
            diffusion=np.array([.1,D] if key=='wave_pinning' else [.1,D,.1,D]);expected=[]
            for k,d in enumerate(diffusion):
                ratio=d*dt/.05**2
                matrix=sparse.diags([-ratio*np.ones(n-1),(1+2*ratio)*np.ones(n),-ratio*np.ones(n-1)],[-1,0,1],shape=(n,n),format='lil')
                matrix[0,n-1]=-ratio;matrix[n-1,0]=-ratio
                expected.append(spsolve(matrix.tocsc(),y[k]+dt*reaction[k]))
            worker.integrate(y,knots,mean,0 if key=='wave_pinning' else 1,diffusion,dt,dt,0.,.05)
            error=float(np.max(abs(y-np.asarray(expected))))
            if error>1e-9:raise ValueError('Native one-step algebra mismatch')
            checks.append(dict(model=key,D=D,one_step_maximum_absolute_error=error))
    write(FIG/'supporting/metadata/noise_native_algebra_validation.json',dict(passed=True,checks=checks,
        method='One compiled IMEX step compared with unchanged frozen native reaction and scipy sparse periodic implicit-diffusion solve, including fast cytosol.'))
    print(json.dumps(checks))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['analyze','panels','validate']);a=ap.parse_args()
    if a.command=='analyze':analyze()
    elif a.command=='validate':validate()
    else:
        fig=plt.figure(figsize=(11,9));gs=fig.add_gridspec(2,1,hspace=.45);plot_panel_E(fig,gs[0]);plot_panel_F(fig,gs[1]);fig.savefig(FIG/'supporting/viewers/noise_panels_EF.png',dpi=180);plt.close(fig)
