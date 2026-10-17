"""Resting-only 120-s area-penalty screen; report never starts simulations."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validation_campaign import frame

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'results' / 'lambda_sweep_01'
LAMBDAS = np.logspace(-9, -4, 21)
NAMES = ['Cdc42', 'Rac', 'Rho', 'Cdc42 inactive', 'Rac inactive', 'Rho inactive', 'Arp2/3'] + [f'{k}{i}' for k in ('F', 'B', 'P') for i in range(6)]

def run(index):
    lam = float(LAMBDAS[index]); folder = ROOT / f'lambda_{index:02d}'
    folder.mkdir(parents=True, exist_ok=False)
    cmd = [str(ROOT/'maree2006'), '--steps', '1200', '--seed', '1101', '--no-pulse',
           '--xi', '0.06', '--area-mode', 'cpm_lattice_area', '--lambda-cpm', str(lam),
           '--binary-snapshots', '--snapshot-seconds', '120', '--output', str(folder)]
    start = time.monotonic()
    with (folder/'simulation.log').open('w') as log:
        code = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)
    (folder/'execution.json').write_text(json.dumps(dict(lambda_cpm=lam, index=index,
        command=cmd, exit_code=code, wall_seconds=time.monotonic()-start), indent=2))
    if code: raise RuntimeError(f'Simulation exited {code}')
    for p in folder.glob('frame_*.bin'):
        t, ids, v = frame(p)
        np.savez_compressed(p.with_suffix('.npz'), time_s=t, ids=ids, fields=v)
        with np.load(p.with_suffix('.npz')) as z:
            assert np.array_equal(z['fields'], v) and np.array_equal(z['ids'], ids)
        p.unlink()

def coords(ids):
    q=ids%500; r=ids//500
    # Initial cells are central and resting: no periodic seam is approached.
    return .1*(q+.5*r-375), .1*np.sqrt(3)/2*(r-250)

def report():
    records=[]; fields=[]; missing=[]
    for index,lam in enumerate(LAMBDAS):
        folder=ROOT/f'lambda_{index:02d}'
        if not (folder/'frame_1200.npz').exists():
            missing.append(index); continue
        d=np.genfromtxt(folder/'metrics.csv', delimiter=',', names=True)
        assert d['time_s'][-1]==120
        with np.load(folder/'frame_1200.npz') as z:
            ids=z['ids'].copy(); v=z['fields'].copy()
        late=d['time_s']>=90
        masses=[]
        for a,b in [('C_mass','Ci_mass'),('R_mass','Ri_mass'),('rho_mass','rhoi_mass')]:
            pool=d[a]+d[b]; masses.append(float(np.max(np.abs(pool/pool[0]-1))))
        s=dict(index=index, lambda_cpm=float(lam), final_area_percent=float(d['sites'][-1]/300),
            min_area_percent=float(d['sites'].min()/300), max_area_percent=float(d['sites'].max()/300),
            late_area_drift_percent_per_min=float(np.polyfit(d['time_s'][late],d['sites'][late]/300,1)[0]*60),
            extensions=int(d['extensions'].sum()), retractions=int(d['retractions'].sum()),
            late_copies_per_second=float((d['extensions'][late]+d['retractions'][late]).sum()/(late.sum()*.1)),
            max_pool_relative_drift=max(masses), min_field=float(v.min()), nonfinite=int((~np.isfinite(v)).sum()),
            mean_late_extension_area_dH=float(d['mean_extension_area_dH'][late].mean()),
            mean_late_retraction_area_dH=float(d['mean_retraction_area_dH'][late].mean()))
        records.append(s); fields.append((folder,ids,v,d))
    if not records: raise RuntimeError('No completed runs')
    # Per-field limits shared across every lambda. All 25 fields remain separate.
    caps=np.max(np.stack([v.max(axis=1) for _,_,v,_ in fields]),axis=0)
    caps=np.maximum(caps,1e-30)
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for s,(_,_,_,d) in zip(records,fields):
        axes[0,0].plot(d['time_s'],d['sites']/300,label=f"{s['lambda_cpm']:.2g}")
    axes[0,0].set(xlabel='Time (s)',ylabel='Area (% target)',title='Resting-cell area')
    axes[0,0].axhspan(95,105,color='green',alpha=.12)
    axes[0,0].legend(ncol=3,fontsize=6)
    xs=[s['lambda_cpm'] for s in records]
    for ax,key,label in [(axes[0,1],'final_area_percent','Final area (% target)'),
                          (axes[1,0],'late_area_drift_percent_per_min','Last 30 s drift (% target / min)'),
                          (axes[1,1],'late_copies_per_second','Last 30 s accepted copies / s')]:
        ax.semilogx(xs,[s[key] for s in records],'o-');ax.set(xlabel='λ (site-count Hamiltonian)',ylabel=label)
    axes[0,1].axhspan(95,105,color='green',alpha=.12)
    axes[1,0].axhline(0,color='gray',lw=.7)
    fig.savefig(ROOT/'overview.png',dpi=140);plt.close(fig)
    cards=[]
    # All-species snapshots, with fixed coordinates and fixed per-field color limits.
    for s,(folder,ids,v,d) in zip(records,fields):
        x,y=coords(ids)
        fig,axes=plt.subplots(5,5,figsize=(16,15),layout='constrained')
        for j,ax in enumerate(axes.flat):
            sc=ax.scatter(x,y,c=v[j],s=1,vmin=0,vmax=caps[j],cmap='viridis',rasterized=True)
            ax.set(xlim=(-14,14),ylim=(-14,14),aspect='equal',title=NAMES[j],xlabel='µm',ylabel='µm')
            fig.colorbar(sc,ax=ax,shrink=.7)
        fig.suptitle(f"λ={s['lambda_cpm']:.4g} · 120 s · {len(ids)} sites · all 25 fields")
        fig.savefig(folder/'final_all_fields.png',dpi=120);plt.close(fig)
        cards.append(f'<details><summary>λ={s["lambda_cpm"]:.4g}: final area {s["final_area_percent"]:.2f}%</summary><a href="{folder.name}/final_all_fields.png">Full-resolution snapshot</a> · <a href="{folder.name}/frame_1200.npz">All final fields (NPZ)</a> · <a href="{folder.name}/metrics.csv">Metrics</a> · <a href="{folder.name}/manifest.txt">Parameters</a><img loading="lazy" src="{folder.name}/final_all_fields.png"></details>')
    # Compact final-cell comparison; full field panels are linked below.
    fig,axes=plt.subplots(3,7,figsize=(18,8),layout='constrained')
    for ax in axes.flat: ax.set_visible(False)
    for s,(_,ids,v,_) in zip(records,fields):
        ax=axes.flat[s['index']];ax.set_visible(True);x,y=coords(ids)
        ax.scatter(x,y,c=v[0],s=.4,vmin=0,vmax=caps[0],cmap='viridis',rasterized=True)
        ax.set(xlim=(-14,14),ylim=(-14,14),aspect='equal',title=f"λ={s['lambda_cpm']:.2g}\n{s['final_area_percent']:.1f}% area")
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Final Cdc42 fields at 120 s — identical position, spatial scale and color scale')
    fig.savefig(ROOT/'final_comparison.png',dpi=150);plt.close(fig)
    rows=''.join(f"<tr><td>{s['lambda_cpm']:.4g}</td><td>{s['final_area_percent']:.2f}</td><td>{s['late_area_drift_percent_per_min']:.3f}</td><td>{s['extensions']}</td><td>{s['retractions']}</td><td>{s['late_copies_per_second']:.2f}</td></tr>" for s in records)
    page=f'''<!doctype html><meta charset="utf-8"><title>Marée area-penalty screen</title>
<style>body{{font:16px system-ui;max-width:1250px;margin:30px auto;padding:15px}}img{{width:100%}}td,th{{padding:8px;border-bottom:1px solid #ddd}}details{{margin:15px 0;padding:10px;border:1px solid #ccc}}</style>
<h1>Resting-cell λ screen</h1><p>{len(records)}/21 complete. Missing indices: {missing}.</p>
<p><strong>Update: <a href="../lambda_bias_01/index.html">five-seed replicated sweep with membrane angular-bias statistics and all final-field snapshots</a>.</strong></p>
<p>120 s, no polarization cue, ξ=0.06, seed 1101, 30,000-site target, 500×500 lattice. Only λ varies; 21 simultaneous array tasks, four CPU threads each. No moving-cell runs, longer confirmations, or GIFs.</p>
<p>This is a short empirical screen, not identification of the published coefficient. A near-target area with continuing drift is not equilibrium; suppressed copies can indicate an overconstrained boundary. Resting controls alone cannot establish that migration remains possible.</p>
<img src="overview.png"><table><tr><th>λ</th><th>Final area %</th><th>Late drift % target/min</th><th>Extensions</th><th>Retractions</th><th>Late copies/s</th></tr>{rows}</table>
<p>Max GTPase pool drift: {max(s['max_pool_relative_drift'] for s in records):.3g}; nonfinite final values: {sum(s['nonfinite'] for s in records)}. Conservation is post-normalization, not independent validation of remapping. <a href="summaries.json">Numerical summaries</a>.</p>
<h2>Final-cell comparison</h2><img src="final_comparison.png">
<h2>All final fields, each λ</h2><p>Seven soluble fields in µM; F/P orientations in nm⁻¹; B orientations in nm⁻². Each field has the same color limits across λ. Coordinates in µm, fixed about the initial cell center. All 25 float64 fields and site IDs are retained at 0 and 120 s.</p>{''.join(cards)}'''
    (ROOT/'index.html').write_text(page)
    (ROOT/'summaries.json').write_text(json.dumps(records,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['run','report']);p.add_argument('--index',type=int)
    a=p.parse_args()
    if a.action=='run': run(a.index)
    else: report()
