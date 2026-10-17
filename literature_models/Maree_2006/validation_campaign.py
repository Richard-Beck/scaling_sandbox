"""Parallel campaign runner and offline report; no simulations during reporting."""
import argparse, csv, html, json, os, subprocess, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE/'results'/'validation_01'
CONDITIONS=[('rest',.06,False,300)]+[(f'xi_{x:.2f}',x,True,600) for x in (0,.01,.03,.06,.09,.12)]
TASKS=[dict(name=f'{name}_seed{seed}',condition=name,xi=x,pulse=pulse,seconds=seconds,seed=seed)
       for name,x,pulse,seconds in CONDITIONS for seed in (1101,1102,1103)]

def frame(path):
    with open(path,'rb') as f:
        step,n,k=np.fromfile(f,dtype='<i4',count=3)
        ids=np.fromfile(f,dtype='<i4',count=n)
        values=np.fromfile(f,dtype='<f8',count=n*k).reshape(k,n)
    assert len(ids)==n and k==25
    return step*.1,ids,values

def run(index):
    task=TASKS[index]; folder=ROOT/task['name'];folder.mkdir(parents=True,exist_ok=True)
    (folder/'task.json').write_text(json.dumps(task,indent=2))
    cmd=[str(ROOT/'maree2006'), '--steps',str(task['seconds']*10),'--seed',str(task['seed']),
         '--xi',str(task['xi']),'--area-mode','cpm_lattice_area','--lambda-cpm','1e-9',
         '--binary-snapshots','--snapshot-seconds','5','--output',str(folder)]
    if not task['pulse']:cmd+=['--no-pulse']
    start=time.monotonic()
    with (folder/'simulation.log').open('w') as log:
        code=subprocess.call(cmd,stdout=log,stderr=subprocess.STDOUT)
    (folder/'execution.json').write_text(json.dumps(dict(exit_code=code,wall_seconds=time.monotonic()-start,command=cmd)))
    if code:raise RuntimeError(f'simulation failed: {code}')
    analyze(task)

def analyze(task):
    folder=ROOT/task['name']
    data=np.genfromtxt(folder/'metrics.csv',delimiter=',',names=True)
    paths=sorted(folder.glob('frame_*.bin'),key=lambda p:int(p.stem.split('_')[1]))
    images=[];trajectory=[];polarity=[];minimum=1e100;nonfinite=0;frame_areas=[]
    dims=['Cdc42 (µM)','Rac (µM)','Rho (µM)','F total (nm⁻¹)','B total (nm⁻²)','P total (nm⁻¹)']
    caps=[2.4,7.5,3.1,.5,5e-5,.15]
    for path in paths:
        t,ids,v=frame(path);q=ids%500;r=ids//500
        minimum=min(minimum,float(np.nanmin(v)));nonfinite+=int((~np.isfinite(v)).sum())
        cq=np.angle(np.exp(2j*np.pi*q/500).mean())*500/(2*np.pi)%500
        cr=np.angle(np.exp(2j*np.pi*r/500).mean())*500/(2*np.pi)%500
        dq=(q-cq+250)%500-250;dr=(r-cr+250)%500-250
        x=.1*(dq+.5*dr);y=.1*np.sqrt(3)/2*dr
        trajectory.append((t,cq,cr))
        frame_areas.append(len(ids))
        cent=[]
        for f in (0,1,2):cent.append(np.array([np.sum(x*v[f]),np.sum(y*v[f])])/v[f].sum())
        polarity.append(float(np.linalg.norm(cent[0]-cent[2])))
        fig,axes=plt.subplots(2,3,figsize=(10,6),layout='constrained')
        fields=[v[0],v[1],v[2],v[7:13].sum(0),v[13:19].sum(0),v[19:25].sum(0)]
        for ax,z,title,cap in zip(axes.flat,fields,dims,caps):
            sc=ax.scatter(x,y,c=z,s=1,vmin=0,vmax=cap,cmap='viridis',rasterized=True)
            ax.set(xlim=(-22,22),ylim=(-19,19),aspect='equal',title=title,xlabel='µm',ylabel='µm')
            fig.colorbar(sc,ax=ax,shrink=.7)
        fig.suptitle(f"{task['name']} · {t:.0f} s · {len(ids)} sites · cell-centered")
        fig.canvas.draw();images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert('RGB'))
        if path==paths[-1]:fig.savefig(folder/'final.png',dpi=110)
        plt.close(fig)
    images[0].save(folder/'fields.gif',save_all=True,append_images=images[1:],duration=120,loop=0)
    tr=np.array(trajectory);uq=np.unwrap(tr[:,1]*2*np.pi/500)*500/(2*np.pi);ur=np.unwrap(tr[:,2]*2*np.pi/500)*500/(2*np.pi)
    xy=np.column_stack((.1*(uq+.5*ur),.1*np.sqrt(3)/2*ur));xy-=xy[0]
    speeds=np.linalg.norm(np.diff(xy,axis=0),axis=1)/np.diff(tr[:,0])*1000
    late=tr[1:,0]>=task['seconds']/2
    area=data['sites'];mass_error={}
    for name,a,b in [('C','C_mass','Ci_mass'),('Rac','R_mass','Ri_mass'),('Rho','rho_mass','rhoi_mass')]:
        pool=data[a]+data[b];mass_error[name]=float(np.max(np.abs(pool/pool[0]-1)))
    copies=int((data['extensions']+data['retractions']).sum())
    record=dict(**task,completed_seconds=float(data['time_s'][-1]),copies=copies,copies_per_second=copies/task['seconds'],
        area_min=int(area.min()),area_max=int(area.max()),area_cv=float(area.std()/area.mean()),
        area_mean_bias=float(area.mean()/30000-1),late_speed_nm_s=float(speeds[late].mean()),
        displacement_um=float(np.linalg.norm(xy[-1])),final_polarity_um=polarity[-1],
        max_pool_relative_error=mass_error,min_field=minimum,nonfinite=nonfinite,
        trajectory_time_s=tr[:,0].tolist(),trajectory_um=xy.tolist(),speeds_nm_s=speeds.tolist(),polarity_um=polarity,
        frame_sites=frame_areas,area_within_five_percent=bool(np.all(np.abs(area/30000-1)<=.05)))
    (folder/'summary.json').write_text(json.dumps(record,indent=2))
    # Compress losslessly after analysis, retaining all 25 orientation-resolved fields.
    for path in paths:
        t,ids,v=frame(path);np.savez_compressed(path.with_suffix('.npz'),time_s=t,ids=ids,fields=v)
        with np.load(path.with_suffix('.npz')) as check:assert np.array_equal(check['fields'],v,equal_nan=True)
        path.unlink()

def report():
    ROOT.mkdir(parents=True,exist_ok=True);records=[];cards=[]
    for task in TASKS:
        folder=ROOT/task['name'];sp=folder/'summary.json'
        if not sp.exists():
            cards.append(f'<p>{task["name"]}: incomplete; see task logs.</p>');continue
        s=json.loads(sp.read_text())
        metrics=np.genfromtxt(folder/'metrics.csv',delimiter=',',names=True)
        s['frame_sites']=np.interp(s['trajectory_time_s'],metrics['time_s'],metrics['sites']).tolist()
        s['area_within_five_percent']=bool(np.all(np.abs(metrics['sites']/30000-1)<=.05))
        xy=np.array(s['trajectory_um']);times=np.array(s['trajectory_time_s'])
        late_xy=xy[times>=s['seconds']/2];late_times=times[times>=s['seconds']/2]
        path_length=np.linalg.norm(np.diff(late_xy,axis=0),axis=1).sum()
        net_length=np.linalg.norm(late_xy[-1]-late_xy[0])
        s['late_net_speed_nm_s']=float(1000*net_length/(late_times[-1]-late_times[0]))
        s['late_persistence']=float(net_length/path_length) if path_length else 0.
        s['final_area_fraction']=float(metrics['sites'][-1]/30000)
        s['simulation_wall_seconds']=json.loads((folder/'execution.json').read_text())['wall_seconds']
        manifest={}
        if (folder/'manifest.txt').exists():
            manifest=dict(line.split('=',1) for line in (folder/'manifest.txt').read_text().splitlines() if '=' in line)
        s['max_geometry_correction']=max([float(v) for k,v in manifest.items() if k.startswith('max_pre_correction_')],default=0.)
        records.append(s)
        expanded=' open' if s['name'] in ('rest_seed1101','xi_0.06_seed1101') else ''
        cards.append(f'''<details{expanded}><summary>{s['name']} · late speed {s['late_speed_nm_s']:.1f} nm/s · area {s['area_min']}–{s['area_max']} · copies {s['copies']}</summary>
<p>Five-second sampling; color scales shared across every run. Coordinates follow the cell; trajectories below retain translation.</p>
<img loading="lazy" src="{s['name']}/final.png"><p><a href="{s['name']}/fields.gif">Open animated fields</a> · <a href="{s['name']}/metrics.csv">Every-MCS metrics</a> · <a href="{s['name']}/summary.json">Summary</a> · <a href="{s['name']}/manifest.txt">Parameters and geometry audit</a></p>
<img loading="lazy" src="{s['name']}/fields.gif"></details>''')
    fig,axes=plt.subplots(3,2,figsize=(12,13),layout='constrained')
    colors={name:plt.cm.tab10(i) for i,(name,*_) in enumerate(CONDITIONS)}
    for s in records:
        xy=np.array(s['trajectory_um']);t=np.array(s['trajectory_time_s']);label=s['condition'] if s['seed']==1101 else '_nolegend_';color=colors[s['condition']]
        axes[0,0].plot(xy[:,0],xy[:,1],label=label,color=color,alpha=.75)
        axes[0,1].plot(t,s['polarity_um'],label=label,color=color,alpha=.75)
        axes[1,0].scatter(s['xi'],s['late_speed_nm_s'],marker='x' if not s['pulse'] else 'o',c='gray' if not s['pulse'] else 'tab:blue')
        axes[1,1].scatter(s['xi'],100*s['area_mean_bias'],marker='x' if not s['pulse'] else 'o',c='gray' if not s['pulse'] else 'tab:blue')
        axes[2,0].plot(t,np.array(s['frame_sites'])/30000,label=label,color=color,alpha=.75)
        axes[2,1].plot(t[1:],s['speeds_nm_s'],label=label,color=color,alpha=.75)
    axes[0,0].set(title='Periodic-unwrapped centroid trajectories',xlabel='µm',ylabel='µm',aspect='equal')
    axes[0,1].set(title='Cdc42–Rho concentration-centroid separation',xlabel='s',ylabel='µm')
    axes[1,0].set(title='Late-half mean 5-s centroid speed (gray × = resting)',xlabel='ξ',ylabel='nm/s');axes[1,0].axhspan(70,110,color='green',alpha=.15)
    axes[1,1].set(title='Mean area bias',xlabel='ξ',ylabel='% of target')
    axes[2,0].set(title='Area over time',xlabel='s',ylabel='Fraction of target');axes[2,0].axhspan(.95,1.05,color='green',alpha=.15)
    axes[2,1].set(title='5-second centroid speed',xlabel='s',ylabel='nm/s')
    axes[0,1].legend(fontsize=8)
    fig.savefig(ROOT/'overview.png',dpi=150);plt.close(fig)
    rows=''.join(f"<tr><td>{s['name']}</td><td>{s['late_speed_nm_s']:.1f}</td><td>{s['late_net_speed_nm_s']:.1f}</td><td>{s['late_persistence']:.2f}</td><td>{100*s['area_cv']:.2f}</td><td>{100*s['area_mean_bias']:.2f}</td><td>{s['copies_per_second']:.1f}</td><td>{s['final_polarity_um']:.2f}</td><td>{max(s['max_pool_relative_error'].values()):.2g}</td></tr>" for s in records)
    groups=[]
    for name,*_ in CONDITIONS:
        group=[s for s in records if s['condition']==name]
        if group:
            speeds=[s['late_speed_nm_s'] for s in group]
            groups.append(f"<tr><td>{name}</td><td>{len(group)}</td><td>{np.mean(speeds):.1f} ± {np.std(speeds,ddof=1) if len(group)>1 else 0:.1f}</td><td>{min(s['area_min'] for s in group)/300:.1f}–{max(s['area_max'] for s in group)/300:.1f}%</td><td>{np.mean([s['final_polarity_um'] for s in group]):.2f}</td></tr>")
    bad=[s for s in records if not s['area_within_five_percent']]
    findings=f"{len(bad)} of {len(records)} analyzed runs leave the ±5% area window. "
    if bad:findings+='Area instability limits interpretation of centroid speed as sustained cell migration.'
    audits=''
    if records:
        audits=f"""<h2>Piggyback numerical audits</h2><p>Maximum relative GTPase-pool drift: {max(max(s['max_pool_relative_error'].values()) for s in records):.3g}.
        Nonfinite values in saved fields: {sum(s['nonfinite'] for s in records)}; minimum saved concentration/density: {min(s['min_field'] for s in records):.3g}.
        Largest pre-normalization geometry discrepancy: {max(s['max_geometry_correction'] for s in records):.3g} (relative).
        Numerical completion and conserved pools do not establish correct mechanics or local moving-boundary transfers.</p>
        <p>Simulation wall time per run: {min(s['simulation_wall_seconds'] for s in records)/60:.1f}–{max(s['simulation_wall_seconds'] for s in records)/60:.1f} min, four CPU threads each; excludes plotting/compression. <a href="summaries.json">Combined machine-readable results</a> · <a href="README.md">Reproduction and snapshot format</a></p>"""
    outcomes=[]
    rest=[s for s in records if s['condition']=='rest']
    standard=[s for s in records if s['condition']=='xi_0.06']
    if rest:
        outcomes.append(f"Resting controls: late centroid speeds {min(s['late_speed_nm_s'] for s in rest):.2f}–{max(s['late_speed_nm_s'] for s in rest):.2f} nm/s; final area {min(s['final_area_fraction'] for s in rest)*100:.1f}–{max(s['final_area_fraction'] for s in rest)*100:.1f}% of target. Maximum final biochemical centroid separation {max(s['final_polarity_um'] for s in rest):.2g} µm.")
    if standard:
        outcomes.append(f"Standard polarized runs (ξ=0.06): late speeds {min(s['late_speed_nm_s'] for s in standard):.1f}–{max(s['late_speed_nm_s'] for s in standard):.1f} nm/s; late net speeds {min(s['late_net_speed_nm_s'] for s in standard):.1f}–{max(s['late_net_speed_nm_s'] for s in standard):.1f} nm/s; final area {min(s['final_area_fraction'] for s in standard)*100:.1f}–{max(s['final_area_fraction'] for s in standard)*100:.1f}% of target.")
    body=f'''<!doctype html><meta charset="utf-8"><title>Marée 2006 validation</title><style>body{{font:16px system-ui;max-width:1200px;margin:35px auto;padding:20px;color:#203040}}img{{width:100%}}table{{border-collapse:collapse}}td,th{{padding:8px;border-bottom:1px solid #ddd}}details{{padding:16px;border:1px solid #ddd;margin:12px 0}}summary{{cursor:pointer}}</style>
<h1>Marée 2006: resting, polarization, and contractility</h1><p>{len(records)} / {len(TASKS)} runs analyzed.</p>
<p><strong>{findings}</strong></p>
{''.join('<p>'+item+'</p>' for item in outcomes)}
<p>The <a href="../../spec_addendum_01.md">area addendum</a> explicitly leaves the numerical rescaling of λ unresolved. The lattice-site interpretation restores boundary copies, but these longer tests assess whether λ=10⁻⁹ actually regulates size. Area drift is consistent with insufficient restoring strength; this campaign alone does not isolate that from other mechanics or coupling errors. No parameters were retuned during these tests.</p>
<p>For scale, 10,000 excess sites produce an extension area penalty λ[2(N−N₀)+1] ≈ 2×10⁻⁵, compared with CPM temperature 0.008. This is a weak restoring signal. The recorded per-copy energy diagnostics in metrics.csv permit further investigation without rerunning this campaign.</p>
<p>CPU implicit multirate model; 500×500 hex lattice; 30,000-site target; lattice-area λ=10⁻⁹. Three seeds per condition. Resting controls run 300 s with no cue. Polarized runs and the ξ sweep run 600 s, with a 15% Cdc42 activation gradient for the first 10 s. Standard polarized case: ξ=0.06. All simulations submitted together in an unthrottled SLURM array.</p>
<p>The 70–110 nm/s band is the specification's approximate migration target. Speed is the mean distance per 5-s interval over the final half of each run, not instantaneous CPM jitter. Area CV includes the complete run and can reflect drift. Polarity is Cdc42–Rho spatial centroid separation. Short-run agreement does not establish the published phenotype; inspect the movies and resting replicates.</p>
<p>Expected contractility fingerprint: ξ≤0.01, little net translation and wobble; ξ≈0.03, irregular slow movement; ξ≈0.06, persistent keratocyte-like migration; ξ≥0.08, stronger contraction and triangular/unstable shapes. Resting cells should remain circular, unpolarized and stationary. These are qualitative specification targets, not fitted classifications. Three seeds are exploratory rather than a statistical validation.</p>
<table><tr><th>Condition</th><th>Replicates</th><th>Late speed mean ± SD, nm/s</th><th>Area range, % target</th><th>Final polarity, µm</th></tr>{''.join(groups)}</table>
<img src="overview.png"><p>Net speed uses late-half endpoint displacement; persistence is endpoint displacement / path length over that same interval. Both can still be affected by asymmetric deformation. The mean 5-s speed includes the interval ending at the halfway point.</p><table><tr><th>Run</th><th>Late speed nm/s</th><th>Net speed nm/s</th><th>Persistence</th><th>Area CV %</th><th>Area bias %</th><th>Copies/s</th><th>Polarity µm</th><th>Pool drift</th></tr>{rows}</table>
{audits}
<h2>Detailed fields and movies</h2>{''.join(cards)}
<p>Raw compressed NPZ snapshots in each run directory contain site IDs and all 25 fields: C, Rac, Rho, Ci, Ri, rhoi, Arp2/3, F0–5, B0–5, P0–5. Snapshots every 5 s; scalar diagnostics every 0.1 s. Geometry normalization magnitudes are in manifests. GTPase conservation is audited on these same simulations; orientation-specific geometry transfers are not independently validated by a post-normalization pool check.</p>'''
    (ROOT/'index.html').write_text(body);(ROOT/'summaries.json').write_text(json.dumps(records,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['run','report']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.action=='run':run(a.index)
    else:report()
