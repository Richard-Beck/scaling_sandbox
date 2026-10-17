"""Atomic, manually refreshed report; safe for concurrent completing workers."""
import json, html, os, fcntl
from pathlib import Path
from datetime import datetime, timezone
HERE=Path(__file__).resolve().parent;OUT=HERE/'results'
def publish():
    with (OUT/'.publish.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        manifest=json.loads((OUT/'manifest.json').read_text());data={}
        for j in manifest['jobs']:
            p=OUT/(j['key']+'.json')
            if p.exists():data[j['key']]=json.loads(p.read_text())
        now=datetime.now(timezone.utc).isoformat();rows=[];details=[]
        for source in manifest['selections']:
            for study in ['Xuan2018','Mallin2023','Jemaa2017','Jemaa2023','Hosny2026']:
                mask=[i for i in source['indices']if study=='Hosny2026' or i!=5]
                key=source['family']+'_'+('-'.join(map(str,mask))if mask else 'shared')+'_'+study
                r=data.get(key)
                if not r:continue
                comparisons=[]
                for i,e in zip(source['indices'],source['exponents']):
                    name={0:'tp',1:'tm',5:'force'}[i]
                    if name not in r['ratios']:comparisons.append(name+': unobserved without cue');continue
                    ratio=r['ratios'][name];weak=abs(e*__import__('math').log(source['length_span']))<.05
                    agree=(ratio-1)*e>0
                    comparisons.append(f'{name}: {ratio:.3g}×; '+('negligible polarity trend' if weak else 'same direction' if agree else 'opposite direction'))
                validation_ok=all(all(x['within'] for x in v['rows'])for v in r['validations']) if study!='Hosny2026' else None
                status='speed only: underdetermined' if r['limited_information'] else ('all fit + validation intervals' if r['all_reading_intervals'] and validation_ok else 'misses ≥1 reading interval') if study!='Hosny2026' else 'quartile fit: descriptive'
                null=data.get(source['family']+'_shared_'+study)
                import math
                near=[k for k,v in r['parameters'].items() if min(math.log(v/r['bounds'][k][0]),math.log(r['bounds'][k][1]/v)) < .03*math.log(r['bounds'][k][1]/r['bounds'][k][0])]
                if near: status+='; near bounds: '+', '.join(near)
                rows.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>'for v in [source['model'],source['family'],source['law'],f"{100*source['improvement']:.1f}%",study,status,f"{r['score']:.2f}",f"{null['score']:.2f}" if null else 'pending','; '.join(comparisons)])+'</tr>')
        for key,r in data.items():
            table='<tr><th>Group/block</th><th>Metric</th><th>Observed</th><th>Fitted</th><th>Independent</th><th>Finer step / longer burn</th></tr>'
            for j,o in enumerate(r['rows']):
                label=str(o.get('block',o['group']))+(' Q'+str(o['quantile'])if 'quantile'in o else '')
                vals=[label,o['metric'],f"{o['obs']:.4g}",f"{o['predicted']:.4g}"]+[f"{v['rows'][j]['predicted']:.4g}"for v in r['validations']]
                table+='<tr>'+''.join('<td>'+html.escape(v)+'</td>'for v in vals)+'</tr>'
            details.append(f'<details><summary>{key} — score {r["score"]:.2f}</summary><p>Parameters: {html.escape(str(r["parameters"]))}. Shared amplitude {r["shared_amplitude"]:.4g}. Optimizer: {sum(s["success"]for s in r["starts"])}/4 search starts converged; polish converged: {r["polish_success"]}.</p><table>{table}</table><p>Correlated checks: {html.escape(str(r["checks"]))}</p><a href="{key}.json">Full result, bounds, starts and checks</a></details>')
        document='''<!doctype html><meta charset="utf-8"><title>OU bridge: revised joint WGD fits</title><style>body{font:16px system-ui;margin:30px;max-width:1500px;color:#172735}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border:1px solid #ccd4dc;text-align:left}th{background:#e9f1f6}details{margin:15px 0}summary{cursor:pointer}img{max-width:100%}a{color:#075ca4}</style><h1>OU bridge: joint WGD refits</h1>'''
        document+=f'<p>Generated {now}. Completed {len(data)}/{len(manifest["jobs"])} distinct fits. Refresh manually for newly completed results.</p>'
        document+='''<p><a href="spatial.html">Experiment-specific training/held-out heatmaps and grouped trajectory inspection</a> · <a href="FINDINGS.md">First-pass findings</a> · <a href="comparison.csv">Comparison CSV</a> · <a href="audit.json">Numerical audit</a> · <a href="manifest.json">Frozen selection manifest</a> · <a href="../../results/report.html">Original campaign (superseded literature screen)</a></p>
<p>Each study is fitted independently. Both groups share every identifiable parameter except those named by the selected nonconstant polarity law. All values are reoptimized; no numerical prior and no directional constraint is imposed. Constant laws are shared-population nulls. Fits with identical permitted differences are reused across spatial sources. A larger effective length is assumed for the polyploid/PACC group when comparing signs.</p>
<p>We select the lowest calibration-error nonconstant law per source/candidate, even when its advantage over constant is weak or absent. Less than 10% improvement is exploratory. A polarity change below 5% across tested lengths is labelled negligible. No length-law magnitude is transferred. Parameter values and bounds are available below. Latent ps=1 fixes a scale gauge; the shared readout amplitude is fitted analytically. Gain is unidentifiable without imposed cue; Hosny fits a constant force (gain × cue), not the two factors separately.</p>
<p>Xuan/Mallin success requires all fitted summaries inside the recorded figure-reading intervals, with the same fitted amplitude in independent simulations and the finer-step check. These intervals are not confidence intervals. A matching best-fit parameter direction alone is provisional, not identified mechanistic agreement. Jemaà speed-only fits cannot establish dynamics. Hosny fits control-group quartiles with residuals divided by the observed IQR; its score is descriptive, not a statistical acceptance test. Knockdown responses are not fitted here.</p>
<p>Near-bound parameters are flagged when within 3% of the allowed log-range from either limit. A strong boundary dependence is a reason to revisit bounds and the model before drawing a mechanistic conclusion.</p><p>Numerics: four bounded starts, common random numbers; search 64 independent cells for nonlinear models, polish 256, validation 1024. Linear models use exact Gaussian displacement covariance (256/1024/4096 cells; speed and MSD analytic). Nonlinear validation changes seed, halves the time step and doubles the burn-in. Polarity-derived starts are supplemented by matched and perturbed starts. Finite bounds and local multistart optimization do not establish impossibility when a fit fails.</p>
<table><tr><th>Spatial source</th><th>OU candidate</th><th>Nonconstant law</th><th>Calibration improvement</th><th>WGD study</th><th>Adequacy / boundary diagnostic</th><th>Score</th><th>Shared-null score</th><th>Poly/reference parameter ratio and direction</th></tr>'''+''.join(rows)+'</table><h2>Distinct fits and observable predictions</h2>'+''.join(details)
        temp=OUT/'report.html.tmp';temp.write_text(document);temp.replace(OUT/'report.html')
if __name__=='__main__':publish()
