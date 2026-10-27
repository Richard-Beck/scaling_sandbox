"""Immutable preliminary snapshots of completed four-seed Poisson cohorts."""
from __future__ import annotations
import argparse,csv,hashlib,importlib.util,json,math,os,shlex,shutil,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript';RUN=MS/'data/runs/poisson_ring_modes_20261027'
OUT=MS/'figures/length_mechanisms_20261008/dev/poisson_ring_modes_20261027'
PY='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
MODELS=('wave_pinning','goryachev','holmes_model3','spring_ring')
NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes',spring_ring='Spring')


def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp')
    q.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n');q.replace(p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def table(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def helper():
    p=RUN/'frozen_sources/production.py';spec=importlib.util.spec_from_file_location('preliminary_legacy_analysis',p)
    O=importlib.util.module_from_spec(spec);sys.modules[spec.name]=O;spec.loader.exec_module(O);O.RUN=RUN
    O.verify();return O


def prepare():
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');snap=RUN/'preliminary'/stamp
    for part in ('results','logs'):(snap/part).mkdir(parents=True)
    m=read(RUN/'manifest.json');coverage=[];ready=[];count=0
    for c in m['conditions']:
        completed=[]
        for slot in c['slots']:
            p=RUN/'workers'/f'{slot:03d}'/'segments/00/status.json'
            if p.exists() and read(p).get('complete'):completed.append(slot)
        count+=len(completed)
        coverage.append(dict(condition=c['index'],model=c['model'],size=c['size'],rate_level=c['rate_level'],
                             amplitude_level=c['amplitude_level'],finished=len(completed),required=4,
                             plotted=len(completed)==4))
        if len(completed)==4:ready.append(c['index'])
    if not ready:raise RuntimeError('No complete four-seed conditions yet')
    shutil.copy2(__file__,snap/'producer.py')
    write(snap/'manifest.json',dict(generated_utc=now(),stamp=stamp,producer_sha256=sha(snap/'producer.py'),
          source_manifest_sha256=sha(RUN/'manifest.json'),coverage=coverage,ready=ready,
          finished_runs=count,expected_runs=432,extension_level=0,
          scope='Completed four-seed initial-window cohorts only; no extensions or incomplete cohorts; preliminary snapshot, not final inference'))
    cmd=[PY,'-B',str(snap/'producer.py'),'worker','--snapshot',str(snap),'--slot']
    script=snap/'worker.sh';script.write_text('#!/bin/bash\nset -euo pipefail\nexec '+shlex.join(cmd)+' "$SLURM_ARRAY_TASK_ID"\n');script.chmod(0o750)
    base=['sbatch','--parsable','--partition=red','--qos=small','--cpus-per-task=1','--mem=4G','--time=01:00:00']
    command=base+['--job-name=poisson_prelim','--array='+','.join(map(str,ready)),
                  '--output='+str(snap/'logs/%A_%a.out'),'--error='+str(snap/'logs/%A_%a.err'),str(script)]
    r=subprocess.run(command,capture_output=True,text=True,check=True);jid=r.stdout.strip().split(';')[0]
    render=base+['--job-name=poisson_prelim_plot','--dependency=afterany:'+jid,'--output='+str(snap/'logs/%j.out'),
                 '--error='+str(snap/'logs/%j.err'),'--wrap',shlex.join([PY,'-B',str(snap/'producer.py'),'render','--snapshot',str(snap)])]
    r=subprocess.run(render,capture_output=True,text=True,check=True);plotjob=r.stdout.strip().split(';')[0]
    receipt=dict(snapshot=str(snap),analysis_job_id=jid,plot_job_id=plotjob,complete_cohorts=len(ready),finished_runs=count,commands=[command,render])
    write(snap/'submission.json',receipt);write(RUN/'preliminary/latest.json',receipt);print(json.dumps(receipt))


def worker(snap,ci):
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Use SLURM for Fourier analysis')
    s=read(snap/'manifest.json');assert ci in s['ready'] and sha(snap/'producer.py')==s['producer_sha256']
    O=helper();E=O.engine();m=read(RUN/'manifest.json');c=m['conditions'][ci];j=m['jobs'][c['slots'][0]]
    rows=[];powers=[];data=[];C=c['size'];sources={}
    if c['model']=='spring_ring':C=float(np.mean([np.load(O.segment(slot,0)/'circumference.npy',mmap_mode='r').mean() for slot in c['slots']]))
    for slot in c['slots']:
        folder=O.segment(slot,0);assert read(folder/'status.json')['complete']
        sources[str((folder/'status.json').relative_to(ROOT))]=sha(folder/'status.json')
        data.append(np.load(folder/'u_modes.npy',mmap_mode='r'))
    for n in range(j['maximum_mode']+1):
        raw=np.array([E.acf_raw(z[:,n]) for z in data]);curve=raw.mean(0)/raw[:,0].mean()
        gamma,kind,bound=E.half_rate(curve,j['frame_interval'])
        lo,hi,los,his,p=O.bootstrap_rates(raw,j['frame_interval'])
        rows.append(dict(condition=ci,model=c['model'],size=c['size'],mean_physical_C=C,rate_level=c['rate_level'],
                         amplitude_level=c['amplitude_level'],seeds=4,mode=n,wavelength=C/n if n else None,
                         gamma=gamma,status=kind,bound=bound,gamma_ci_low=lo,gamma_ci_high=hi,
                         ci_low_status=los,ci_high_status=his,bootstrap_resolved_probability=p,
                         record_duration=j['sample_duration'],analysis_lag=(len(curve)-1)*j['frame_interval'],
                         mode_power=float(raw[:,0].mean()),source='initial_window_preliminary'))
    table(snap/'results'/f'{ci:03d}.csv',rows)
    write(snap/'results'/f'{ci:03d}.json',dict(complete=True,generated_utc=now(),condition=ci,
          source_status_sha256=sources,frames_per_seed=[len(z) for z in data],modes=len(rows)))


def render(snap):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullFormatter
    s=read(snap/'manifest.json');stamp=s['stamp'];rows=[];missing=[]
    for ci in s['ready']:
        p=snap/'results'/f'{ci:03d}.json'
        if not p.exists() or not read(p).get('complete'):missing.append(ci);continue
        with (snap/'results'/f'{ci:03d}.csv').open() as f:rows+=list(csv.DictReader(f))
    if missing:
        write(snap/'status.json',dict(state='INCOMPLETE',missing=missing,generated_utc=now()))
        raise RuntimeError('Missing preliminary analysis, not censoring: '+str(missing))
    for part in ('data','metadata','vectors','viewers'):(OUT/'supporting'/part).mkdir(parents=True,exist_ok=True)
    prefix='poisson_ring_modes_preliminary_'+stamp
    datafile=OUT/'supporting/data'/f'{prefix}_rates.csv';coveragefile=OUT/'supporting/data'/f'{prefix}_coverage.csv'
    table(datafile,rows);table(coveragefile,s['coverage']);artifacts=[datafile,coveragefile];plots=[]
    levels=('low','medium','high');cmap=plt.get_cmap('Blues',5)
    fig,axes=plt.subplots(2,2,figsize=(12,5.6),layout='constrained')
    for ax,model in zip(axes.flat,MODELS):
        sizes=sorted({r['size'] for r in s['coverage'] if r['model']==model});arr=np.zeros((3,9),int)
        for r in s['coverage']:
            if r['model']==model:arr[sizes.index(r['size']),3*levels.index(r['rate_level'])+levels.index(r['amplitude_level'])]=r['finished']
        im=ax.imshow(arr,vmin=0,vmax=4,cmap=cmap,aspect='auto')
        for yy in range(3):
            for xx in range(9):ax.text(xx,yy,str(arr[yy,xx]),ha='center',va='center',color='white' if arr[yy,xx]>=3 else '#183848',fontsize=10)
        ax.set_xticks(range(9),[f'{r[0].upper()}/{a[0].upper()}' for r in levels for a in levels],fontsize=8)
        ax.set_yticks(range(3),[f'C={v:g}' for v in sizes]);ax.set_title(NAMES[model]);ax.set_xlabel('Arrival rate / mean amplitude')
    fig.suptitle(f"PRELIMINARY snapshot: {s['finished_runs']}/432 finished; entries are completed seeds out of 4",fontsize=12)
    coveragepng=OUT/f'{prefix}_coverage.png';fig.savefig(coveragepng,dpi=170);plt.close(fig);artifacts.append(coveragepng)
    summary=[]
    for model in MODELS:
        selected=[r for r in rows if r['model']==model and int(r['mode'])>0]
        sizes=sorted({float(r['size']) for r in s['coverage'] if r['model']==model})
        colors=dict(zip(sizes,('#6ca6c9','#287bb0','#133955')))
        fig,axes=plt.subplots(3,3,figsize=(11,9.5),sharex=True,sharey=True,layout='constrained')
        for ir,rate in enumerate(levels):
            for ia,amp in enumerate(levels):
                ax=axes[ir,ia];q=[r for r in selected if r['rate_level']==rate and r['amplitude_level']==amp]
                if not q:ax.set_facecolor('#f4f4f4');ax.text(.5,.45,'Waiting for complete\nfour-seed cohorts',ha='center',va='center',transform=ax.transAxes,fontsize=9,color='#666')
                for r in q:
                    size=float(r['size']);resolved=r['status']=='resolved';x=float(r['wavelength']);y=float(r['gamma'] if resolved else r['bound'])
                    marker='o' if resolved else 'v' if r['status']=='beyond_analysis_lag' else '^'
                    if r['ci_low_status']==r['ci_high_status']=='resolved':ax.vlines(x,float(r['gamma_ci_low']),float(r['gamma_ci_high']),color=colors[size],alpha=.45,lw=.75)
                    ax.plot(x,y,marker,ms=3.5,color=colors[size],mfc=colors[size] if resolved else 'white',mew=.75)
                cov=[r for r in s['coverage'] if r['model']==model and r['rate_level']==rate and r['amplitude_level']==amp]
                label='  '.join(f"C{r['size']:g}: {r['finished']}/4" for r in sorted(cov,key=lambda r:r['size']))
                ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.15)
                ticks=[.2,.5,1,2,3] if model=='spring_ring' else [2,5,10,20,50]
                ax.set_xticks(ticks,[f'{v:g}' for v in ticks],fontsize=9)
                ax.xaxis.set_minor_formatter(NullFormatter());ax.tick_params(axis='y',labelsize=9)
                ax.set_title(f'{rate} rate; {amp} amplitude\n{label}',fontsize=8.5)
                if ir==2:ax.set_xlabel('Wavelength λ (native length)')
                if ia==0:ax.set_ylabel('Decorrelation rate γ (native time⁻¹)')
        handles=[Line2D([],[],marker='o',ls='',color=colors[v],label=f'Nominal C={v:g}',ms=4) for v in sizes]
        handles+=[Line2D([],[],marker='v',ls='',color='#555',mfc='white',label='Upper rate bound',ms=5),Line2D([],[],marker='^',ls='',color='#555',mfc='white',label='Lower rate bound',ms=5)]
        fig.legend(handles=handles,loc='outside lower center',ncol=5,fontsize=8,frameon=False)
        ncohort=sum(r['plotted'] for r in s['coverage'] if r['model']==model)
        fig.suptitle(f"PRELIMINARY — {NAMES[model]} — {ncohort}/27 complete four-seed cohorts\nInitial sampling windows; extension checks pending",fontsize=12)
        png=OUT/f'{prefix}_{model}.png';pdf=OUT/'supporting/vectors'/f'{prefix}_{model}.pdf'
        fig.savefig(png,dpi=170);fig.savefig(pdf);plt.close(fig);plots.append((model,png,pdf));artifacts+=[png,pdf]
        summary.append(dict(model=model,complete_cohorts=ncohort,resolved_points=sum(r['status']=='resolved' for r in selected),bound_points=sum(r['status']!='resolved' for r in selected)))
    reader=OUT/'supporting/viewers'/f'{prefix}.html'
    html=['<!doctype html><html><head><meta charset="utf-8"><title>Preliminary Poisson ring spectra</title><style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;line-height:1.5}img{width:100%;height:auto}.status{background:#f0f4f7;padding:12px}math{font-size:1.15em}.eq{display:flex;justify-content:center;position:relative;margin:20px 0}.eq>span{position:absolute;right:0;top:50%;transform:translateY(-50%)}</style></head><body>',
          '<h1>Preliminary Poisson ring decorrelation spectra</h1>',
          f"<p class='status'>Snapshot: {s['generated_utc']}. Plot generation: {now()}. Preliminary producer finished; simulations continue. No timed refresh. {s['finished_runs']}/432 histories were complete at snapshot; {len(s['ready'])}/108 complete four-seed conditions enter these plots.</p>",
          '<p>Only conditions with all four initial-window repeats complete are plotted. Counts show missing cohorts explicitly. Completion speed can depend on input and model dynamics, so this is an incomplete experimental grid. All selected seeds enter the pooled raw numerator and zero-lag power; per-seed rates are not averaged. Whole-seed bootstrap intervals are pointwise. Error bars are drawn when both interval endpoints are resolved; endpoint censor statuses remain in the CSV. Open downward/upward triangles show upper/lower rate bounds, not measured rates. No asymptotic or multi-exponential decay fit is imposed.</p>',
          '<div class="eq"><math display="block"><mrow><msub><mi>γ</mi><mi>n</mi></msub><mo>=</mo><mfrac><mrow><mi>ln</mi><mn>2</mn></mrow><msub><mi>ℓ</mi><mfrac><mn>1</mn><mn>2</mn></mfrac></msub></mfrac><mo>,</mo><mspace width="2em"/><msub><mi>λ</mi><mi>n</mi></msub><mo>=</mo><mfrac><mi>C</mi><mi>n</mi></mfrac></mrow></math><span>(1)</span></div>',
          '<p>Chemical rings discard 16,000 and sample 128,000 native time units at interval 0.5. Spring discards 1,600 and samples 12,800 at interval 0.05. The maximum analysis lag is half the sampled record. The half-height lag is linearly interpolated between its adjacent saved-frame lags. Modes 1 through floor(C/ell) are plotted; ell=1.5 chemical and .15 Spring. Spring wavelength uses mean sampled physical circumference; legend C is its nominal resting circumference. Native times are not cross-model biological calibrations.</p>',
          f'<p><a href="../data/{datafile.name}">Plotted rates and bounds CSV</a> · <a href="../data/{coveragefile.name}">Completion coverage CSV</a></p>',
          f'<img src="../../{coveragepng.name}" alt="Number of finished repeats at each condition">']
    for model,png,pdf in plots:
        html+=[f'<h2>{NAMES[model]}</h2><p><a href="../vectors/{pdf.name}">Vector PDF</a></p><img src="../../{png.name}" alt="Preliminary {NAMES[model]} wavelength decorrelation plots">']
    html+=['</body></html>'];reader.write_text('\n'.join(html));artifacts.append(reader)
    receipt=dict(generated_utc=now(),scope='preliminary initial-window snapshot',finished_runs=s['finished_runs'],complete_cohorts=len(s['ready']),by_model=summary,
                 reader=str(reader),plots=[str(png) for _,png,_ in plots],coverage_plot=str(coveragepng))
    write(snap/'status.json',dict(state='COMPLETE',**receipt))
    write(OUT/'supporting/metadata'/f'{prefix}.manifest.json',dict(artifact_root='../..',generated_utc=now(),sha256={str(p.relative_to(OUT)):sha(p) for p in artifacts}))
    write(OUT/'supporting/metadata/preliminary_latest.json',receipt)
    with (OUT/'README.md').open('a') as f:f.write(f'\n[Preliminary initial-window plots: {stamp}](supporting/viewers/{reader.name})\n')
    subprocess.run([PY,'-B',str(MS/'code/mechanisms/organize_length_dev.py'),'--group','poisson_ring_modes_20261027'],check=True,capture_output=True,text=True)
    print(json.dumps(receipt))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=('prepare','worker','render'))
    ap.add_argument('--snapshot',type=Path);ap.add_argument('--slot',type=int);a=ap.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='worker':worker(a.snapshot,a.slot)
    else:render(a.snapshot)
