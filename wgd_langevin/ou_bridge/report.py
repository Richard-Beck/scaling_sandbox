"""Publish current completed data without browser refresh timers."""
from datetime import datetime,timezone
import csv
import html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from design import *
from reduced import *
from fit import load_cases,prepare,residual,evaluate
from scaling import parameters

def csvwrite(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def build():
    from characterize import build as characterize
    characterize()
    now=datetime.now(timezone.utc).isoformat()
    fits=[json.loads(p.read_text()) for p in sorted((OUT/'fits').glob('*.json'))]
    scales=[json.loads(p.read_text()) for p in sorted((OUT/'scaling').glob('*.json'))]
    literature=[json.loads(p.read_text()) for p in sorted((OUT/'literature').glob('*.json'))]
    validation=[json.loads(p.read_text()) for p in sorted((OUT/'literature_validation').glob('*.json'))]
    validation_rows=[dict(model=r['model'],family=r['family'],**c)for r in validation for c in r['checks']]
    pilot=[];diagnostics=[];audits=[];lengthrows=[];litrows=[];viewer={}
    for f in fits:
        common=dict(model=f['model'],family=f['family'],length=f['length'])
        pilot.extend([dict(**common,**r)for r in f['metrics']])
        prepared,_=prepare(load_cases(f['model'],f['length'],True))
        rr=residual(encode(f['theta'],f['family']),f['family'],prepared)
        recost=float(rr@rr/2)
        sv=min(f['starts'],key=lambda r:r['cost'])['jacobian_singular_values']
        diagnostics.append(dict(**common,cost=f['cost'],recomputed_cost=recost,
             objective_agrees=abs(recost-f['cost'])<=1e-6*max(1.,f['cost']),
             converged_starts=sum(s['success']for s in f['starts']),
             near_best_starts=sum(s['cost']<=1.01*f['cost']+1e-10 for s in f['starts']),
             jacobian_condition=sv[0]/max(sv[-1],1e-30),max_step_halving_error=f['max_step_halving_error']))
    for model in MODELS:
        middle=lengths(model,True)[1]
        for assay in ASSAYS:
            fine=OUT/'refinement'/f'{model}_{middle:g}'/f'{assay}.npz'
            if not fine.exists():continue
            a=np.load(OUT/'spatial'/f'{model}_{middle:g}'/f'{assay}.npz');b=np.load(fine)
            delta=np.interp(a['time'],b['time'],b['p'])-a['p']
            audits.append(dict(model=model,length=middle,assay=assay,
                               max_coordinate_error=float(np.max(np.abs(delta))),
                               rms_coordinate_error=float(np.sqrt(np.mean(delta**2))),
                               final_error=float(abs(delta[-1]))))
    for s in scales:
        selected=next(x for x in s['laws'] if x['law']==s['selected_law'])
        for law in s['laws']:
            for row in law['validation']:
                lengthrows.append(dict(model=s['model'],family=s['family'],law=law['law'],
                                  selected=law['law']==s['selected_law'],cost=law['cost'],**row))
        for L in lengths(s['model']):
            cases=load_cases(s['model'],L)
            _,pred=evaluate(parameters(selected,L),s['family'],cases,s['scale'])
            for c in cases:
                key=f"{s['model']} | {L:g} um | {c['assay']}"
                idx=np.unique(np.linspace(0,len(c['time'])-1,min(500,len(c['time']))).astype(int))
                if key not in viewer:
                    viewer[key]=dict(time=c['time'][idx].tolist(),target=c['p'][idx].tolist(),
                                     cue=[Protocol(s['model'],c['assay']).projected(t)for t in c['time'][idx]],
                                     families={},role='calibration length' if L in lengths(s['model'],True) else 'held-out length')
                viewer[key]['families'][s['family']]=pred[c['assay']][idx,0].tolist()
    for item in literature:
        for r in item['results']:
            litrows.append(dict(model=item['model'],family=item['family'],study=r['study'],
                 hypothesis=r['hypothesis'],score=r['score'],status=r['status'],
                 reference_length=r['reference_length'],polyploid_length=r['polyploid_length'],
                 k_hours_per_unit=r['k_hours_per_unit'],noise=r['noise'],
                 compatible_grid_points=r.get('compatible_grid_points',''),
                 simulation_valid=r['simulation_valid']))
    csvwrite('pilot_metrics.csv',pilot);csvwrite('fit_diagnostics.csv',diagnostics)
    csvwrite('spatial_refinement.csv',audits);csvwrite('length_validation.csv',lengthrows)
    csvwrite('literature_summary.csv',litrows)
    csvwrite('literature_independent_checks.csv',validation_rows)
    matrix=np.full((5,4),np.nan)
    for i,m in enumerate(MODELS):
        for j,f in enumerate(FAMILIES):
            values=[r['post_nrmse']for r in lengthrows if r['model']==m and r['family']==f and r['selected']
                    and r['length_role']=='held_out' and r['assay']in ['tracking_test','reversal']]
            if values:matrix[i,j]=np.mean(values)
    fig,ax=plt.subplots(figsize=(9,5));im=ax.imshow(matrix,cmap='magma_r',vmin=0,vmax=1)
    ax.set_xticks(range(4),['Gaussian\nOU','Linear\nmemory','Coupled','Coupled +\nactivation'])
    ax.set_yticks(range(5),MODELS)
    for i in range(5):
        for j in range(4):
            if np.isfinite(matrix[i,j]):ax.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',color='white' if matrix[i,j]>.5 else 'black')
    ax.set_title('Untuned response error at held-out lengths\nMean normalized RMSE: reversal and unused tracking period')
    fig.colorbar(im,ax=ax,label='NRMSE (lower is better; display capped at 1)')
    fig.tight_layout();fig.savefig(OUT/'performance_matrix.png',dpi=170);plt.close(fig)
    # All candidate laws are available in CSV; this figure shows selected-law transfer.
    fig,axes=plt.subplots(5,4,figsize=(15,15),sharey=True)
    colors={'pulse':'#1464a0','reversal':'#e4572e','tracking_test':'#559e54'}
    for i,m in enumerate(MODELS):
        for j,f in enumerate(FAMILIES):
            ax=axes[i,j]
            for assay,color in colors.items():
                rows=sorted([r for r in lengthrows if r['model']==m and r['family']==f and r['selected'] and r['assay']==assay],key=lambda r:r['length'])
                ax.plot([r['length']for r in rows],[r['post_nrmse']for r in rows],color=color,label=assay)
                for r in rows:ax.scatter(r['length'],r['post_nrmse'],marker='o' if r['length_role']=='calibration' else 's',color=color,s=25)
            ax.set_title(f'{m}: {f}',fontsize=9);ax.set_yscale('symlog',linthresh=.01);ax.grid(alpha=.15)
            ax.set_xlabel('Length (um)')
            if j==0:ax.set_ylabel('Normalized response error')
    axes[0,0].legend(fontsize=7);fig.suptitle('Joint length laws: circles calibration lengths, squares held-out lengths')
    fig.tight_layout();fig.savefig(OUT/'length_transfer.png',dpi=150);plt.close(fig)
    retained=json.loads((OUT/'retained_transfer.json').read_text()) if (OUT/'retained_transfer.json').exists() else {'rows':[]}
    csvwrite('retained_debelly_transfer.csv',retained['rows'])
    text=[f'# OU bridge campaign\n\nGenerated {now}. No prospective experiment selection performed.\n',
          f'Completed: {len(fits)}/60 pilot fits; {len(scales)}/20 length-law comparisons; {len(literature)}/20 literature screens.\n',
          '[Interpretation and conclusions](CONCLUSIONS.md) summarize what this campaign does and does not establish.\n',
          '## Reading the results\n',
          'The three model families are Gaussian OU, a specified linear exponential-memory realization, and nonlinear coupled feedback with self-activation off/on. These are effective polarity reductions, not replacements for accepted WGD fits.\n',
          'The literature screen fixes the selected length law and scans a common time conversion. Effective length pairs are scenarios, not independently measured sizes. No-compatible-grid-point is not a proof of incompatibility. Hosny is a descriptive baseline quantile comparison, with unchanged knockdowns held out.\n',
          '## Held-out length performance\n',
          '| Spatial model | Gaussian | Memory | Coupled | Self-activation |\n|---|---:|---:|---:|---:|']
    for i,m in enumerate(MODELS):text.append('| '+m+' | '+' | '.join(f'{v:.3f}'for v in matrix[i])+' |')
    text+=['\nValues are mean normalized response RMSE for reversal and the unused tracking period. They do not score full spatial profiles or establish statistical significance. See `pilot_metrics.csv` and `length_validation.csv` for retention, reversal classes, physical-range violations, and competition diagnostics.\n',
           '## Numerical and identification audit\n',
           f'- Pilot objective rechecks: {sum(r["objective_agrees"]for r in diagnostics)}/{len(diagnostics)} agree.',
           f'- Largest reduced-model step-halving coordinate difference: {max((r["max_step_halving_error"]for r in diagnostics),default=0):.6g}.',
           f'- Largest spatial-refinement coordinate difference: {max((r["max_coordinate_error"]for r in audits),default=0):.6g}.',
           f'- Fits with at least one unconverged start: {sum(r["converged_starts"]<8 for r in diagnostics)}/{len(diagnostics)}.',
           '- Near-best start counts and Jacobian singular values diagnose ambiguity; they are not confidence intervals. Model rankings are conditional on this bounded search.',
           '\n## Literature screen\n',
           '| Model source | Family | Study | Hypothesis | Status | Score |\n|---|---|---|---|---|---:|']
    for r in litrows:text.append('| '+' | '.join(str(r[k])for k in ['model','family','study','hypothesis','status'])+f' | {r["score"]:.3g} |')
    text+=[f'\nIndependent literature checks completed for {len(validation)}/20 model/source combinations. See `literature_independent_checks.csv`; initial grid hits are not robust if independent predictions leave the supplied intervals or stationarity fails.\n',
           '\nExact parameter bounds, stochastic closures and scope limits are in `../README.md`. Nonlinear stationarity audits and complete predictions are saved per model in `literature/`. Grid compatibility remains provisional until independent simulation checks pass.\n',
           '## Files\n\n- `report.html`: offline interactive trajectory inspection.\n- `performance_matrix.png`, `length_transfer.png`: shareable figures.\n- `effective_dynamics.csv`: zero-cue fixed points and scalar local relaxation times across length; these are not calibrated 2D turning-persistence times.\n- `fits/`, `scaling/`, `literature/`: full parameters, starts, validation and provenance.\n- `retained_debelly_transfer.csv`: untuned transfer to the original retained Figure 2 protocols.\n']
    (OUT/'REPORT.md').write_text('\n'.join(text))
    payload=json.dumps(viewer,allow_nan=False)
    page='''<!doctype html><meta charset="utf-8"><title>OU bridge campaign</title>
<style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;color:#202b35}select{padding:8px;width:100%}svg{width:100%;background:#fafafa}small{color:#555}.legend{display:flex;gap:20px;flex-wrap:wrap}img{max-width:100%}</style>
<h1>OU bridge: polarity responses and literature consistency</h1><p>GENERATED</p>
<p>Three families; nonlinear feedback is shown with and without self-activation. No prospective experiment selection. This is an effective-model campaign, not a replacement for the accepted WGD analyses.</p>
<p><a href="CONCLUSIONS.md">Interpretation</a> · <a href="REPORT.md">Full report</a> · <a href="length_validation.csv">Length validation CSV</a> · <a href="literature_summary.csv">Literature screen CSV</a> · <a href="../README.md">Methods and limitations</a></p>
<img src="performance_matrix.png" alt="Held-out response error matrix">
<h2>Inspect trajectories</h2><select id="choice"></select><p id="role"></p><div class="legend" id="legend"></div>
<svg id="plot" viewBox="0 0 1100 430"></svg><small>Time in seconds; signed half-mass polarity. No clipping. Cue acts directly on the polarity coordinate. Latent q is not equated to measured MCA.</small>
<h2>Length transfer</h2><img src="length_transfer.png" alt="Length transfer">
<p>Literature compatibility is conditional on the stated noise and size scenarios. Failed stationarity checks mean unresolved simulation; a finite-grid miss does not prove biological incompatibility. Refresh manually to load newly generated results.</p>
<script>const D=PAYLOAD;const colors={target:'#111',gaussian:'#2274a5',memory:'#dd7826',coupled:'#3a923a',self_activation:'#9b52ad'};
let select=document.getElementById('choice');Object.keys(D).sort().forEach(k=>{let o=document.createElement('option');o.textContent=k;select.append(o)});
function draw(){let d=D[select.value];if(!d)return;document.getElementById('role').textContent=d.role;
let series={target:d.target,...d.families};let vals=Object.values(series).flat();let lo=Math.min(-.02,...vals),hi=Math.max(.02,...vals);let pad=.05*(hi-lo);lo-=pad;hi+=pad;
let X=t=>65+990*(t-d.time[0])/(d.time.at(-1)-d.time[0]);let Y=v=>370-330*(v-lo)/(hi-lo);
let s='<line x1="65" y1="370" x2="1055" y2="370" stroke="#aaa"/>';
for(let i=0;i<=4;i++){let t=d.time[0]+i*(d.time.at(-1)-d.time[0])/4,v=lo+i*(hi-lo)/4;s+=`<text x="${X(t)}" y="398" font-size="12">${t.toFixed(0)}</text><text x="5" y="${Y(v)}" font-size="12">${v.toFixed(3)}</text><line x1="65" x2="1055" y1="${Y(v)}" y2="${Y(v)}" stroke="#eee"/>`}
for(let [name,values]of Object.entries(series)){let points=values.map((v,i)=>X(d.time[i])+','+Y(v)).join(' ');s+=`<polyline points="${points}" fill="none" stroke="${colors[name]}" stroke-width="${name==='target'?3:1.8}"/>`}
document.getElementById('plot').innerHTML=s;document.getElementById('legend').innerHTML=Object.keys(series).map(k=>`<span style="color:${colors[k]}">${k==='target'?'Spatial target':k}</span>`).join('');}
select.onchange=draw;draw();</script>'''
    (OUT/'report.html').write_text(page.replace('GENERATED',html.escape('Data generated '+now+'. No timed refresh.')).replace('PAYLOAD',payload))
    save_json(OUT/'report_status.json',dict(generated=now,pilot_fits=len(fits),length_fits=len(scales),literature_screens=len(literature),literature_checks=len(validation),
              all_objectives_agree=all(r['objective_agrees']for r in diagnostics)))
    print('Report generated',now,flush=True)

if __name__=='__main__':build()
