"""Replicated lambda screen and geometry-adjusted membrane angular diagnostics."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validation_campaign import frame
from lambda_sweep import LAMBDAS, NAMES

HERE=Path(__file__).resolve().parent
ROOT=HERE/'results'/'lambda_bias_01'
SEEDS=(2101,2102,2103,2104,2105)
THETA=(np.arange(36)+.5)*2*np.pi/36
HARMONICS=(1,2,6)
WIDTHS=(1,3,5)
LABELS=NAMES+['F total','B total','P total']+[f'{g}{m} covariant residual' for g in ('F','B','P') for m in range(6)]

def geometry(ids):
    q=ids%500;r=ids//500
    cq=np.angle(np.exp(2j*np.pi*q/500).mean())*500/(2*np.pi)
    cr=np.angle(np.exp(2j*np.pi*r/500).mean())*500/(2*np.pi)
    dq=(q-cq+250)%500-250;dr=(r-cr+250)%500-250
    x=.1*(dq+.5*dr);y=.1*np.sqrt(3)/2*dr
    neighbors=np.array([((r+b)%500)*500+(q+a)%500 for a,b in ((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))])
    occupied=np.zeros(250000,bool);occupied[ids]=True
    inner=occupied.copy();depth=np.zeros(len(ids),int)
    for layer in range(1,6):
        retained=inner[ids]&np.all(inner[neighbors],axis=0)
        depth[inner[ids]&~retained]=layer
        inner[:]=False;inner[ids[retained]]=True
    sectors=(np.mod(np.arctan2(y,x),2*np.pi)*36/(2*np.pi)).astype(int)%36
    return x,y,depth,sectors

def angular(ids,v):
    x,y,depth,sectors=geometry(ids)
    values=np.concatenate((v,np.stack([v[a:a+6].sum(0) for a in (7,13,19)])))
    profiles=[];coefficients=[];orientation=[];counts=[]
    for width in WIDTHS:
        layer_profiles=[];layer_counts=[]
        for layer in range(1,width+1):
            mask=depth==layer;n=np.bincount(sectors[mask],minlength=36)
            if np.any(n==0):raise ValueError('Empty angular sector: adjust binning')
            layer_counts.append(n)
            layer_profiles.append(np.array([np.bincount(sectors[mask],weights=z[mask],minlength=36)/n for z in values]))
        # Equal angular-sector and depth-layer weighting: constant fields give zero bias,
        # regardless of sector pixel counts or a purely depth-dependent profile.
        p=np.mean(layer_profiles,axis=0);norm=p.mean(axis=1)
        normalized=np.divide(p,norm[:,None],out=np.zeros_like(p),where=norm[:,None]>0)
        residuals=[]
        for start in (7,13,19):
            # Rotate orientation m back by m*60 degrees, preserving relative abundance.
            group=p[start:start+6];group_mean=group.mean()
            aligned=np.array([np.roll(row,-6*m) for m,row in enumerate(group)])
            common=aligned.mean(axis=0)
            residuals.extend([np.roll(aligned[m]-common,6*m)/group_mean if group_mean>0 else np.zeros(36) for m in range(6)])
        all_profiles=np.concatenate((normalized,np.array(residuals)))
        z=np.array([(all_profiles*np.exp(1j*k*THETA)).mean(axis=1) for k in HARMONICS]).T
        # Global orientation imbalance, using geometry-adjusted sector/layer means.
        orient=[]
        for start in (7,13,19):
            shares=p[start:start+6].mean(axis=1);total=shares.sum()
            orient.append([np.sum(shares*np.exp(1j*k*np.arange(6)*np.pi/3))/total if total>0 else 0j for k in (1,2)])
        profiles.append(all_profiles);coefficients.append(z);orientation.append(orient);counts.append(np.array(layer_counts).sum(axis=0))
    return np.array(profiles),np.array(coefficients),np.array(orientation),np.array(counts)

def render(folder,ids,v,lam,seed):
    caps=np.load(ROOT/'color_limits.npy');x,y,_,_=geometry(ids)
    fig,axes=plt.subplots(5,5,figsize=(16,15),layout='constrained')
    for j,ax in enumerate(axes.flat):
        sc=ax.scatter(x,y,c=v[j],s=1,vmin=0,vmax=caps[j],cmap='viridis',rasterized=True)
        ax.set(xlim=(-14,14),ylim=(-14,14),aspect='equal',title=NAMES[j],xlabel='µm',ylabel='µm')
        fig.colorbar(sc,ax=ax,shrink=.7)
    fig.suptitle(f'λ={lam:.4g} · seed {seed} · 120 s · {len(ids)} sites')
    fig.savefig(folder/'final_all_fields.png',dpi=100);plt.close(fig)

def run(index):
    li,si=divmod(index,len(SEEDS));lam=float(LAMBDAS[li]);seed=SEEDS[si]
    folder=ROOT/f'lambda_{li:02d}_seed{seed}';folder.mkdir(exist_ok=False)
    cmd=[str(ROOT/'maree2006'),'--steps','1200','--seed',str(seed),'--no-pulse','--xi','.06',
        '--area-mode','cpm_lattice_area','--lambda-cpm',str(lam),'--binary-snapshots','--snapshot-seconds','120','--output',str(folder)]
    start=time.monotonic()
    with (folder/'simulation.log').open('w') as log:code=subprocess.call(cmd,stdout=log,stderr=subprocess.STDOUT)
    (folder/'execution.json').write_text(json.dumps(dict(command=cmd,exit_code=code,wall_seconds=time.monotonic()-start),indent=2))
    if code:raise RuntimeError(f'Simulation failed {code}')
    for p in folder.glob('frame_*.bin'):
        t,ids,v=frame(p);np.savez_compressed(p.with_suffix('.npz'),time_s=t,ids=ids,fields=v)
        with np.load(p.with_suffix('.npz')) as check:assert np.array_equal(check['fields'],v) and np.array_equal(check['ids'],ids)
        p.unlink()
    with np.load(folder/'frame_1200.npz') as z:ids=z['ids'].copy();v=z['fields'].copy()
    profiles,coefficients,orientation,counts=angular(ids,v)
    np.savez_compressed(folder/'angular.npz',profiles=profiles,coefficients=coefficients,orientation=orientation,
                        sector_counts=counts,widths=WIDTHS,harmonics=HARMONICS,labels=LABELS,theta=THETA)
    d=np.genfromtxt(folder/'metrics.csv',delimiter=',',names=True);late=d['time_s']>=90
    assert d['time_s'][-1]==120
    pools=[]
    for a,b in [('C_mass','Ci_mass'),('R_mass','Ri_mass'),('rho_mass','rhoi_mass')]:
        pool=d[a]+d[b];pools.append(float(np.max(np.abs(pool/pool[0]-1))))
    s=dict(lambda_index=li,lambda_cpm=lam,seed=seed,folder=folder.name,final_area_percent=float(d['sites'][-1]/300),
        late_drift_percent_per_min=float(np.polyfit(d['time_s'][late],d['sites'][late]/300,1)[0]*60),
        late_copies_per_second=float((d['extensions'][late]+d['retractions'][late]).sum()/(late.sum()*.1)),
        pool_drift=max(pools),nonfinite=int((~np.isfinite(v)).sum()),minimum=float(v.min()),
        color_clipped_values=int((v>np.load(ROOT/'color_limits.npy')[:,None]).sum()))
    render(folder,ids,v,lam,seed)
    (folder/'summary.json').write_text(json.dumps(s,indent=2))

def summarize(z,rng):
    mean=z.mean();amps=np.abs(z);n=len(z)
    # Seed bootstrap, never pixel bootstrap; no significance calls or multiplicity claims.
    means=z[rng.integers(0,n,size=(4000,n))].mean(axis=1)
    return dict(coherent_amplitude=float(abs(mean)),mean_seed_amplitude=float(amps.mean()),
        coherence=float(abs(mean)/amps.mean()) if amps.mean()>1e-12 else None,
        phase_degrees=float(np.angle(mean,deg=True)) if abs(mean)>1e-12 else None,
        real=float(mean.real),imag=float(mean.imag),
        real_ci95=np.quantile(means.real,[.025,.975]).tolist(),imag_ci95=np.quantile(means.imag,[.025,.975]).tolist())

def report():
    runs=[];arrays={};stats=[];orientation_stats=[];rng=np.random.default_rng(92817)
    for p in sorted(ROOT.glob('lambda_*/summary.json')):
        s=json.loads(p.read_text());runs.append(s)
        with np.load(p.parent/'angular.npz') as z:arrays[s['folder']]={k:z[k].copy() for k in ('coefficients','profiles','orientation')}
    if not runs:raise RuntimeError('No complete runs')
    for li,lam in enumerate(LAMBDAS):
        group=[s for s in runs if s['lambda_index']==li]
        if not group:continue
        coeff=np.array([arrays[s['folder']]['coefficients'] for s in group])
        ori=np.array([arrays[s['folder']]['orientation'] for s in group])
        for wi,width in enumerate(WIDTHS):
            for fi,label in enumerate(LABELS):
                for ki,k in enumerate(HARMONICS):
                    stats.append(dict(lambda_index=li,lambda_cpm=float(lam),n=len(group),width=width,field=label,harmonic=k,
                                      **summarize(coeff[:,wi,fi,ki],rng)))
            for gi,g in enumerate(('F','B','P')):
                for ki,k in enumerate((1,2)):
                    orientation_stats.append(dict(lambda_index=li,lambda_cpm=float(lam),n=len(group),width=width,field=g,harmonic=k,
                                                   **summarize(ori[:,wi,gi,ki],rng)))
    (ROOT/'summaries.json').write_text(json.dumps(runs,indent=2))
    (ROOT/'angular_statistics.json').write_text(json.dumps(stats,indent=2))
    (ROOT/'orientation_statistics.json').write_text(json.dumps(orientation_stats,indent=2))
    with (ROOT/'angular_statistics.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(stats[0]));writer.writeheader();writer.writerows(stats)
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for ax,key,title in zip(axes,('final_area_percent','late_drift_percent_per_min','late_copies_per_second'),
                           ('Final area (% target)','90–120 s area drift (% target/min)','90–120 s accepted copies/s')):
        for li,lam in enumerate(LAMBDAS):
            vals=[s[key] for s in runs if s['lambda_index']==li]
            if vals:ax.errorbar(lam,np.mean(vals),yerr=np.std(vals,ddof=1) if len(vals)>1 else 0,fmt='o',color='tab:blue',capsize=3)
        ax.set(xscale='log',xlabel='λ',ylabel=title)
    axes[0].axhspan(95,105,color='green',alpha=.12);axes[1].axhline(0,color='gray',lw=.6)
    fig.savefig(ROOT/'area.png',dpi=140);plt.close(fig)
    displayed=NAMES[:7]+['F total','B total','P total']
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    for ki,k in enumerate(HARMONICS):
        for row,key in enumerate(('coherent_amplitude','mean_seed_amplitude')):
            matrix=np.full((len(displayed),21),np.nan)
            for st in stats:
                if st['width']==3 and st['harmonic']==k and st['field'] in displayed:
                    matrix[displayed.index(st['field']),st['lambda_index']]=st[key]*100
            ax=axes[row,ki];im=ax.imshow(matrix,aspect='auto',vmin=0,cmap='magma')
            ax.set_yticks(range(len(displayed)),displayed);ax.set_xticks(range(0,21,4),[f'{LAMBDAS[i]:.0e}' for i in range(0,21,4)],rotation=45)
            ax.set(title=f"k={k}: {'consistent across seeds' if row==0 else 'mean within-seed magnitude'}",xlabel='λ')
            fig.colorbar(im,ax=ax,label='Normalized harmonic magnitude (%)')
    fig.savefig(ROOT/'bias.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for gi,g in enumerate(('F','B','P')):
        for ki,k in enumerate((1,2)):
            st=[s for s in orientation_stats if s['width']==3 and s['field']==g and s['harmonic']==k]
            axes[gi].semilogx([s['lambda_cpm'] for s in st],[100*s['coherent_amplitude'] for s in st],'o-',label=f'k={k} coherent')
            axes[gi].semilogx([s['lambda_cpm'] for s in st],[100*s['mean_seed_amplitude'] for s in st],'--',label=f'k={k} seed mean')
        axes[gi].set(title=f'{g}: orientation-class imbalance',xlabel='λ',ylabel='Magnitude (%)');axes[gi].legend(fontsize=8)
    fig.savefig(ROOT/'orientations.png',dpi=140);plt.close(fig)
    # Complete field-by-field tables, including directional fields and covariance residuals.
    tables=[]
    for li,lam in enumerate(LAMBDAS):
        group=[s for s in runs if s['lambda_index']==li]
        if not group:continue
        rows=[]
        for label in LABELS:
            sel=[s for s in stats if s['lambda_index']==li and s['width']==3 and s['field']==label]
            cells=''.join(f"<td>{100*s['coherent_amplitude']:.3g} / {100*s['mean_seed_amplitude']:.3g}</td>" for s in sel)
            one=sel[0];direction='—' if one['phase_degrees'] is None else f"{one['phase_degrees']:.1f}°"
            rows.append(f'<tr><td>{label}</td>{cells}<td>{direction}</td></tr>')
        links=''.join(f'<p>Seed {s["seed"]}: area {s["final_area_percent"]:.2f}%; <a href="{s["folder"]}/final_all_fields.png">all-field snapshot</a> · <a href="{s["folder"]}/frame_1200.npz">raw final fields</a> · <a href="{s["folder"]}/angular.npz">angular profiles</a> · <a href="{s["folder"]}/metrics.csv">metrics</a></p>' for s in group)
        tables.append(f'<details><summary>λ={lam:.4g} — {len(group)} seeds</summary><table><tr><th>Field</th><th>k=1 coherent / seed mean %</th><th>k=2 coherent / seed mean %</th><th>k=6 coherent / seed mean %</th><th>k=1 direction</th></tr>{"".join(rows)}</table>{links}</details>')
    primary=[s for s in stats if s['width']==3 and s['harmonic']==1]
    actin_max=max(s['coherent_amplitude'] for s in primary if s['field'] in ('F total','B total','P total'))
    gtp_max=max(s['coherent_amplitude'] for s in primary if s['field'] in NAMES[:6])
    candidate=[s for s in runs if s['lambda_index']==16]
    findings=f'Across λ, the largest coherent first-harmonic bias of the actin totals is {100*actin_max:.3f}%; the GTPase fields are below {100*gtp_max:.2g}% (numerical roundoff). These measurements show no large common-direction effect, but do not establish absence of small systematic bias.'
    if candidate:findings+=f" At λ=10⁻⁵, final areas span {min(s['final_area_percent'] for s in candidate):.2f}–{max(s['final_area_percent'] for s in candidate):.2f}% of target."
    page=f'''<!doctype html><meta charset="utf-8"><title>Replicated λ screen: membrane angular bias</title>
<style>body{{font:16px system-ui;max-width:1400px;margin:25px auto;padding:20px}}img{{width:100%}}td,th{{padding:7px;border-bottom:1px solid #ddd}}details{{margin:12px 0;padding:12px;border:1px solid #ccc}}</style>
<h1>Resting cells: replicated λ screen and membrane angular bias</h1><p>{len(runs)}/105 runs complete: 21 λ values × five seeds (2101–2105), 120 s each, no cue, ξ=0.06, same model as the <a href="../lambda_sweep_01/index.html">original screen</a>. All simulations submitted in one unthrottled parallel array. “Fixed area” means a fixed target and finite area penalty, not a clamped boundary.</p>
<p><strong>{findings}</strong></p>
<h2>Area and boundary activity</h2><p>Points are seed means; error bars are seed SD, not confidence intervals.</p><img src="area.png">
<h2>Membrane angular bias</h2><p>Primary band: three inner hex-lattice layers (roughly 0.3 µm); one- and five-layer sensitivity results are included in the downloadable statistics. Angles are measured about the final cell centroid, +x = orientation 0, counterclockwise positive. Within each depth layer, field means are computed in 36 equal angular sectors, then layers receive equal weight. This avoids interpreting unequal pixel counts or purely depth-dependent concentration as angular bias.</p>
<p>For each field, zₖ = meanθ[(sector concentration / angular mean) exp(ikθ)]. k=1 measures a preferred side, k=2 an axis, and k=6 a hexagonal imprint. The coherent statistic is |mean across seeds zₖ|; the within-seed statistic is mean across seeds |zₖ|. A common direction reinforces the first; unrelated directions cancel. For a sinusoid 1+a cos(kθ), the reported magnitude is |a|/2. Numerical values below 10⁻¹² should be treated as zero; phase/coherence are omitted there.</p>
<img src="bias.png"><p>With only five seeds, a nonzero coherent magnitude alone is not evidence of systematic error. Component-wise 95% seed-bootstrap intervals (4,000 draws) and phase/coherence are supplied in JSON/CSV. They are exploratory, unadjusted for multiple comparisons, and not a definitive significance test. Seeds are independent within each λ but paired across λ. Only final-time bias is measured; persistence in time is not established.</p>
<h2>Orientation-resolved actin</h2><p>A single Fm/Bm/Pm field can have a preferred membrane side in a rotationally symmetric resting cell. Raw harmonics of those fields are therefore descriptive, not a failure criterion. We additionally rotate each orientation profile back by m×60°, subtract their common aligned profile, then rotate the residual back. These covariance residuals retain orientation-class abundance differences and test whether equivalent orientations behave equivalently. The orientation-class plot below tests unequal mean abundance among the six classes, independent of their expected spatial patterns. Geometry and lattice effects can still contribute.</p>
<img src="orientations.png"><h2>Full statistics and final snapshots</h2><p><a href="angular_statistics.csv">All fields, λ values, widths, harmonics: CSV</a> · <a href="angular_statistics.json">JSON with bootstrap intervals</a> · <a href="orientation_statistics.json">Orientation-class statistics</a> · <a href="summaries.json">Run summaries</a>. Snapshot color limits are shared across all new runs, set from the prior sweep with 20% headroom; {sum(s['color_clipped_values'] for s in runs)} plotted values exceed those limits. Raw float64 snapshots are unclipped.</p>{''.join(tables)}
<p>Maximum normalized GTPase pool drift: {max(s['pool_drift'] for s in runs):.3g}; nonfinite final values: {sum(s['nonfinite'] for s in runs)}. No model parameters were changed outside the prescribed λ/seed sweep.</p>'''
    (ROOT/'index.html').write_text(page)

def self_test():
    with np.load(HERE/'results/lambda_sweep_01/lambda_00/frame_1200.npz') as z:ids=z['ids']
    x,y,depth,sectors=geometry(ids)
    assert np.all(np.bincount(depth)[1:6]>0)
    v=np.ones((25,len(ids)));_,c,_,_=angular(ids,v)
    assert np.max(np.abs(c))<1e-12
    # Arbitrary depth-only radial gradient remains angularly unbiased.
    v[:]=1+depth;_,c,_,_=angular(ids,v);assert np.max(np.abs(c))<1e-12
    v[:]=1+.2*np.cos(THETA[sectors]);_,c,_,_=angular(ids,v)
    assert np.allclose(c[:,0,0],.1,atol=1e-12)
    v[:]=1
    for start in (7,13,19):
        for m in range(6):v[start+m]=1+.4*np.cos(THETA[sectors]-m*np.pi/3)
    _,c,o,_=angular(ids,v)
    assert np.max(np.abs(c[:,28:]))<1e-12
    assert np.max(np.abs(o))<1e-12
    assert np.allclose(np.abs(c[:,7:13,0]),.2,atol=1e-12)
    print('PASS: constant/radial fields, known dipole, rotational covariance, orientation balance')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['run','report','self-test']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.action=='run':run(a.index)
    elif a.action=='report':report()
    else:self_test()
