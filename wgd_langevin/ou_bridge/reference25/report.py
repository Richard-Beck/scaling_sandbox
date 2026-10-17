"""Publish the current fixed-length fit review, without browser timers."""
import csv
import fcntl
import html
import json
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from campaign import OUT, MODELS, FAMILIES, suite

COLORS=['#277da8','#f8961e','#43aa8b','#9b5de5']
LABELS=['Gaussian OU','Linear memory','Nonlinear coupled','Self-activation']


def main():
    OUT.mkdir(exist_ok=True)
    with (OUT/'.report.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        sources=sorted((OUT/'fits').glob('*.json'))+sorted((OUT/'spatial').glob('*/*.npz'))
        signature=str([(str(p),p.stat().st_mtime_ns) for p in sources])
        signature+=str(__import__('pathlib').Path(__file__).stat().st_mtime_ns)
        marker=OUT/'.report_signature'
        if (marker.exists() and marker.read_text()==signature
                and (OUT/'report.html').exists()
                and marker.stat().st_mtime_ns >= (OUT/'report.html').stat().st_mtime_ns):
            print('Report already current.',flush=True)
            return
        build()
        marker.write_text(signature)


def build():
    generated=datetime.now(timezone.utc).isoformat()
    records={}
    for m in MODELS:
        for f in FAMILIES:
            path=OUT/'fits'/f'{m}_{f}.json'
            if path.exists(): records[m,f]=json.loads(path.read_text())
    images=OUT/'figures'
    images.mkdir(exist_ok=True)
    parts=[f'''<!doctype html><html><head><meta charset="utf-8"><title>25 µm polarity review</title>
    <style>body{{font:16px system-ui;max-width:1550px;margin:auto;padding:24px;color:#24313c}}
    img{{width:100%}}table{{border-collapse:collapse}}td,th{{padding:10px;border:1px solid #ccd}}
    summary{{cursor:pointer;font-weight:bold}}.note{{background:#eef5f8;padding:12px}}</style></head><body>
    <h1>Polarity representation at 25 µm</h1><p>Report generated: {generated}</p>
    <p>{len(records)} / 20 joint fits complete. All 18 stimuli enter each fit. Manual refresh loads the latest completed results.</p>
    <p class="note">p = 2 × x50/L − 1, with x50 the active-field half-mass position; dX/dt = v0 × p.
    LEGI uses its response field. Zero-total-mass frames are assigned p = 0 and flagged in the raw data.
    Fits use one parameter set per spatial model and reduced family. Five stimulus groups have equal weight,
    then equal weight per trajectory and per stimulus segment. No trajectory amplitude normalization,
    WGD fitting, length scaling, or stochastic-noise inference.</p>
    <p>Black: spatial model. Blue: Gaussian OU. Orange: linear memory. Green: nonlinear coupled.
    Purple: self-activation. Upper panels show p; lower panels show the signed cue projection.
    Conditioned protocols show time relative to the perturbation; their complete conditioning histories are fitted.
    Polarity axes share limits within each stimulus group and are not amplitude-normalized.</p>
    <h2>Weighted polarity RMSE</h2><table><tr><th>Spatial model</th>''']
    parts += [f'<th>{label}</th>' for label in LABELS]
    parts.append('</tr>')
    for m in MODELS:
        parts.append(f'<tr><th>{m}</th>')
        for f in FAMILIES:
            r=records.get((m,f))
            if r:
                best=min(r['starts'],key=lambda s:s['cost'])
                flags=('; optimizer limit' if not best['success'] else '')
                if any(x['out_of_bounds'] for x in r['metrics']): flags+='; |p| > 1'
                parts.append(f'<td>{r["weighted_rmse"]:.4f}{flags}</td>')
            else: parts.append('<td>pending</td>')
        parts.append('</tr>')
    parts.append('</table><p>Descriptive fit errors on the absolute polarity scale; no statistical pass/fail threshold.</p>')
    rows=[]
    for m in MODELS:
        predictions={}
        for f in FAMILIES:
            if (m,f) in records:
                with np.load(OUT/'fits'/f'{m}_{f}.npz') as d:
                    predictions[f]={k:d[k] for k in d.files}
                rows += [dict(model=m,family=f,**v) for v in records[m,f]['metrics']]
        cases=suite(m)
        present=[c for c in cases if (OUT/'spatial'/m/(c['name']+'.npz')).exists()]
        if not present: continue
        parts.append(f'<h2>{m}</h2>')
        for group in ['pulse','weak','counter','tracking','competition']:
            selected=[c for c in present if c['group']==group]
            if not selected: continue
            fig,axes=plt.subplots(2,len(selected),figsize=(4.4*len(selected),4.4),
                squeeze=False,gridspec_kw={'height_ratios':[3,1]},layout='constrained')
            limits=[]
            for j,c in enumerate(selected):
                with np.load(OUT/'spatial'/m/(c['name']+'.npz')) as d:
                    t=d['time']; p=d['p']; cue=d['cue']
                mask=t>=c['focus']
                x=t[mask]-c['focus']
                axes[0,j].plot(x,p[mask],color='black',lw=2,label='Spatial')
                limits.extend([p[mask].min(),p[mask].max()])
                for f,color,label in zip(FAMILIES,COLORS,LABELS):
                    if f in predictions:
                        axes[0,j].plot(x,predictions[f][c['name']][mask],color=color,lw=1.3,label=label)
                        v=predictions[f][c['name']][mask]
                        limits.extend([v.min(),v.max()])
                axes[0,j].axhline(0,color='#bbb',lw=.6)
                axes[0,j].set_title(c['name'])
                axes[0,j].set_ylabel('Polarity p')
                axes[1,j].plot(x,cue[mask],color='#777',lw=1)
                axes[1,j].set_ylabel('Cue u')
                axes[1,j].set_xlabel('Time from perturbation (s)' if c['focus'] else 'Time (s)')
                # A logarithmic tail still displays the short pulse and its onset.
                if group in ('pulse','counter','competition'):
                    for ax in axes[:,j]: ax.set_xscale('symlog',linthresh=c['unit']/10)
                for ax in axes[:,j]:
                    ax.grid(alpha=.15)
                    ax.set_xlim(0,x[-1])
            low,high=min(limits),max(limits)
            pad=max(.005,.08*(high-low))
            for ax in axes[0]: ax.set_ylim(low-pad,high+pad)
            axes[0,0].legend(fontsize=7)
            name=f'{m}_{group}.png'
            fig.savefig(images/name,dpi=140)
            plt.close(fig)
            parts.append(f'<h3>{group}</h3><a href="figures/{name}"><img src="figures/{name}"></a>')
        # Competition spatial context: same scalar p need not mean same field.
        competition=[c for c in present if c['group']=='competition']
        if competition:
            fig,axes=plt.subplots(1,len(competition),figsize=(10,3),squeeze=False,layout='constrained')
            for ax,c in zip(axes[0],competition):
                with np.load(OUT/'spatial'/m/(c['name']+'.npz')) as d:
                    pulse=c['segments'][0][1]
                    for t in [pulse, min(10*pulse,c['end']),c['end']]:
                        k=np.argmin(abs(d['time']-t))
                        ax.plot(d['x'],d['profiles'][k],label=f'{d["time"][k]:g} s')
                ax.set(title=c['name'],xlabel='Position (µm)',ylabel='Active field')
                ax.legend(fontsize=8)
            name=f'{m}_competition_profiles.png';fig.savefig(images/name,dpi=140);plt.close(fig)
            parts.append(f'<details><summary>Competition spatial profiles</summary><img src="figures/{name}"></details>')
        parts.append('<details><summary>Parameters and optimizer status</summary><pre>')
        params={f:{'parameters':dict(zip(records[m,f]['parameter_names'],records[m,f]['theta'])),
                   'starts':[dict(seed=s['seed'],cost=s['cost'],success=s['success'],nfev=s['nfev'])
                             for s in records[m,f]['starts']]}
                for f in FAMILIES if (m,f) in records}
        parts.append(html.escape(json.dumps(params,indent=2)))
        parts.append('</pre></details>')
    if rows:
        with (OUT/'metrics.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    parts.append('<p><a href="manifest.json">Campaign specification</a> · <a href="metrics.csv">Per-stimulus fit errors</a></p></body></html>')
    temp=OUT/'report.html.tmp';temp.write_text('\n'.join(parts));temp.replace(OUT/'report.html')
    print(f'Report: {len(records)}/20 fits; {generated}',flush=True)


if __name__=='__main__': main()
