"""Publish each completed fit atomically, with manual refresh and honest diagnostics."""
import csv,json,html,fcntl,io
from pathlib import Path
from datetime import datetime,timezone
OUT=Path(__file__).resolve().parent/'results'
FAMILIES=['gaussian','memory','coupled','self_activation']
def atomic(name,value):
    p=OUT/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(value);t.replace(p)
def esc(x):return html.escape(str(x))
def table(headers,rows):
    return '<table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table>'
def status(r):
    if r['limited_information']:return 'Speed only: underdetermined'
    if r['study'].startswith('Hosny'):return 'Descriptive quartile fit'
    return 'All fitting + validation intervals met' if r['all_reading_intervals'] and all(all(o['within'] for o in v['rows']) for v in r['validations']) else 'At least one fitting/validation interval missed'
def publish():
    with (OUT/'.publish.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        manifest=json.loads((OUT/'manifest.json').read_text());data={}
        for j in manifest['jobs']:
            p=OUT/(j['key']+'.json')
            if p.exists():data[j['key']]=json.loads(p.read_text())
        figmeta=json.loads((OUT/'figure_manifest.json').read_text()) if (OUT/'figure_manifest.json').exists() else {}
        now=datetime.now(timezone.utc).isoformat();studies=list(dict.fromkeys(j['study'] for j in manifest['jobs']))
        summaries=[];detail=[];csvrows=[];obsrows=[];parameter_rows=[]
        for key,r in data.items():
            a,b=r['constraint_audit']['group_parameters']
            for name in a:
                if name in ('gain','ps') or (name in ('tm','r') and r['family']=='gaussian'):continue
                if name in ('vs','self_a','self_K') and name not in r['parameters']:continue
                if name=='force' and not r['study'].startswith('Hosny'):continue
                parameter_rows.append(dict(study=r['study'],family=r['family'],hypothesis='persistence_only' if r['indices'] else 'shared_null',
                    state='complete',parameter=name,reference=a[name],wgd_plus=b[name],ratio=b[name]/a[name] if a[name] else None,
                    shared=a[name]==b[name],fit_score=r['score'],result_file=key+'.json',generated=r['generated']))
        for study in studies:
            sr=[];blocks=[];times=[];shared=[]
            for family in FAMILIES:
                r=data.get(family+'_persistence_'+study);null=data.get(family+'_shared_'+study)
                if r is None:
                    sr.append([family,'Pending','—', '—',f"{null['score']:.3g}" if null else 'pending','—','—','—']);continue
                a,b=r['constraint_audit']['group_parameters']
                fm=lambda x:f'{x:.5g}'
                times.append([family,fm(a['tp']),fm(b['tp']),fm(b['tp']/a['tp']),
                    fm(a['tm']) if family!='gaussian' else '—',fm(b['tm']) if family!='gaussian' else '—',fm(b['tm']/a['tm']) if family!='gaussian' else '—'])
                shared.append([family]+[fm(r['parameters'][k]) if k in r['parameters'] else ('1 (fixed)' if k=='eta' else '—') for k in ['r','eta','vs','self_a','self_K','force']]+[fm(r['shared_amplitude'])])
                diag=r['identifiability'];ratios='; '.join(f'{k}: {v:.3g}' for k,v in r['ratios'].items())
                improvement=100*(1-r['score']/null['score']) if null and null['score']>1e-10 else None
                check_status=(f"{sum(o['within'] for o in r['checks'])}/{len(r['checks'])} intervals" if r['checks'] and not study.startswith('Hosny') else 'see details' if r['checks'] else 'none available')
                sr.append([family,status(r),check_status,f"{r['score']:.3g}",f"{null['score']:.3g}" if null else '—',ratios,
                           f"{diag['local_rank']}/{diag['shape_parameters']}",', '.join(diag['near_bounds']) or 'none'])
                csvrows.append(dict(study=study,family=family,score=r['score'],null_score=null['score'] if null else None,
                    improvement_percent=improvement,status=status(r),tp_ratio=r['ratios'].get('tp'),tm_ratio=r['ratios'].get('tm'),
                    local_rank=diag['local_rank'],shape_parameters=diag['shape_parameters'],near_bounds=';'.join(diag['near_bounds']),
                    validation_min=min(v['score'] for v in r['validations']),validation_max=max(v['score'] for v in r['validations'])))
                rr=[]
                for i,o in enumerate(r['rows']):
                    interval=f"{o['low']:g}–{o['high']:g}" if 'low' in o else f"IQR scale {o['width']:g}"
                    vals=[f"{v['rows'][i]['predicted']:.4g}" for v in r['validations']]
                    rr.append(['WGD+' if o['group'] else 'WGD−',o['metric']+(f" Q{[25,50,75][o['quantile']]}" if 'quantile' in o else ''),f"{o['obs']:.4g}",interval,f"{o['predicted']:.4g}",'; '.join(vals)])
                    obsrows.append(dict(study=study,family=family,role='fit',group=o['group'],metric=o['metric'],quantile=o.get('quantile'),observed=o['obs'],low=o.get('low'),high=o.get('high'),predicted=o['predicted']))
                cc=[]
                for o in r['checks']:
                    cc.append(['WGD+' if o['group'] else 'WGD−',o['metric']+(f" Q{[25,50,75][o['quantile']]}" if 'quantile' in o else ''),f"{o['obs']:.4g}",f"{o['predicted']:.4g}",o.get('within')])
                    obsrows.append(dict(study=study,family=family,role='check',group=o['group'],metric=o['metric'],quantile=o.get('quantile'),observed=o['obs'],low=o.get('low'),high=o.get('high'),predicted=o['predicted']))
                par=[]
                for k,v in r['parameters'].items():par.append([k,f'{v:.5g}',' – '.join(f'{x:g}' for x in r['bounds'][k])])
                blocks.append('<details><summary>'+esc(family)+f" — score {r['score']:.3g}</summary><p>Shared readout amplitude: {r['shared_amplitude']:.5g}. "
                    +f"Converged search starts: {sum(s['success'] for s in r['starts'])}/8; polish converged: {r['polish_success']}. "
                    +f"Constraint audit: {r['constraint_audit']['passes']}. Only changed parameters: {esc(r['constraint_audit']['changed_parameters'])}.</p>"
                    +table(['Group','Measurement','Observed','Reading interval / scale','Fitted','Independent; half-step; longer burn'],rr)
                    +'<h3>Excluded-measurement checks</h3>'+table(['Group','Metric','Observed','Predicted','Within reading interval'],cc)
                    +'<h3>Parameters</h3>'+table(['Parameter','Value','Search bounds'],par)
                    +'<p>Near-best search-start ratio ranges (not confidence intervals): '+esc(diag['nearbest_start_ratio_ranges'])+'.</p>'
                    +f'<p><a href="{r["key"]}.json">Full result and numerical diagnostics</a></p></details>')
            figure=f'{study}.png'
            summaries.append('<section class="study" data-study="'+esc(study)+'"><h2>'+esc(study)+'</h2>'
                +table(['Variant','Fitted-target assessment','Excluded checks','Fit score','Shared-null score','WGD+/WGD− time ratios','Local rank / parameters','Near bounds'],sr)
                +'<h3>Current best persistence times (hours)</h3>'+table(['Variant','tp WGD−','tp WGD+','tp ratio','tm WGD−','tm WGD+','tm ratio'],times)
                +'<h3>Parameters shared within the pair</h3>'+table(['Variant','Feedback r','Noise eta','Saturation vs','Self-activation a','Self-activation K','Directional force','Readout scale'],shared)
                +(f'<p><a href="{figure}"><img src="{figure}" alt="Observed versus predicted {esc(study)}"></a></p>' if (OUT/figure).exists() and figmeta.get(study)==sorted(k for k,r in data.items() if r['study']==study and r['indices']) else '')+''.join(blocks)+'</section>')
        for name,rows in [('comparison.csv',csvrows),('predictions.csv',obsrows),('current_best_parameters.csv',parameter_rows)]:
            if rows:
                s=io.StringIO();w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows);atomic(name,s.getvalue())
        document='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>WGD: persistence-only OU fits</title><style>body{font:16px system-ui;max-width:1450px;margin:30px auto;padding:0 20px;color:#172735}table{border-collapse:collapse;font-size:13px;margin:12px 0}th,td{border:1px solid #ccd5df;padding:8px;text-align:left}th{background:#eef3f7}details{margin:15px 0}summary{cursor:pointer;font-weight:600}img{max-width:100%}select{font:inherit;padding:8px}a{color:#075ca4}</style><h1>Direct WGD fitting: only persistence times differ</h1>'''
        document+=f'<p>Report generated {now}. Fits complete: <b>{len(data)}/{len(manifest["jobs"])}</b>. Latest fit: {max((r["generated"] for r in data.values()),default="none")}. Open or refresh manually to load newly completed results.</p>'
        document+='''<p>Each experiment is fitted independently. Within each WGD−/WGD+ pair, every parameter is shared except tp (all variants) and tm (memory, coupled and self-activation). The null also shares those times. Free noise parameters, feedback, saturation, self-activation, directional forcing and readout amplitude cannot differ between groups. Linear memory retains its existing covariance closure: the derived noise covariance depends on tp/tm, without an extra group-specific noise parameter. No polarity-model fits, lengths, exponents or transferred starting values enter this campaign.</p>
<p><a href="../README.md">Methods and limits</a> · <a href="comparison.csv">Model comparison CSV</a> · <a href="predictions.csv">Predictions CSV</a> · <a href="current_best_parameters.csv">Current best parameters CSV</a> · <a href="manifest.json">Provenance and job manifest</a> · <a href="FINDINGS.md">Findings</a></p>
<p>Scores are descriptive residual sums, comparable between variants within an experiment. Figure-reading intervals are not sampling confidence intervals; Hosny residuals use observed IQRs. No p-values or AIC are inferred. Xuan early MSD, Mallin endpoint, and Hosny endpoint/perpendicular summaries are excluded from fitting. They are correlated checks from the same experiments, not independent biological validation. Hosny control and knockdown pairs are fitted independently; no treatment mechanism is imposed. Unpaired normoxia is excluded.</p>
<p>Eight fresh starts, two best starts polished plus the shared null as a nested candidate. Independent simulations fix the fitted amplitude. Nonlinear checks also halve the step and double burn-in. Local sensitivity rank, parameter bounds and near-best starts diagnose weak identification; they do not establish confidence intervals. Speed-only Jemaà data cannot identify persistence. tp/tm are model times, not necessarily the measured autocorrelation time of a feedback system.</p>'''
        if (OUT/'xuan_memory_ordered/report.html').exists():document+='<p><a href="xuan_memory_ordered/report.html">Background Xuan memory fit: both WGD+ times greater</a></p>'
        if (OUT/'xuan_memory_tracks.html').exists():document+='<p><a href="xuan_memory_tracks.html">Xuan memory model: example tracks, correlation curves and MSD</a></p>'
        for figure in ['comparison.png','persistence_ratios.png']:
            if (OUT/figure).exists() and figmeta.get('all_keys')==sorted(data):document+=f'<p><a href="{figure}"><img src="{figure}" alt="{figure}"></a></p>'
        document+='<p>The parameter tables below show the best completed persistence-only fit for each variant. Pending fits have no parameter estimate displayed. The CSV also includes shared-null parameter sets. Values update automatically on the server when a fit completes; refresh this page manually. These are point estimates, not confidence intervals.</p><label>Experiment <select id="study"><option value="all">All experiments</option>'+''.join('<option>'+esc(s)+'</option>' for s in studies)+'</select></label>'
        document+=''.join(summaries)+'''<script>document.getElementById('study').addEventListener('change',e=>document.querySelectorAll('.study').forEach(s=>s.hidden=e.target.value!=='all'&&s.dataset.study!==e.target.value));</script></html>'''
        atomic('report.html',document)
if __name__=='__main__':publish()
