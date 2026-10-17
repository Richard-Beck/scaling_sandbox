"""Read-only publisher for atomic ordered-campaign checkpoints; no browser timer."""
import csv,json,html,fcntl,io,os,getpass,subprocess,time,signal,argparse
from pathlib import Path
from datetime import datetime,timezone
OUT=Path(__file__).resolve().parent/'results'
FAMILIES=['gaussian','memory','coupled','self_activation']
STOP=False

def now():return datetime.now(timezone.utc).isoformat()
def atomic(name,value):
    p=OUT/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(value);t.replace(p)
def esc(value):return html.escape(str(value))
def fmt(x):return '—' if x is None else f'{x:.5g}'
def table(head,rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'
def write_csv(name,fields,rows):
    s=io.StringIO();w=csv.DictWriter(s,fieldnames=fields);w.writeheader();w.writerows(rows);atomic(name,s.getvalue())

def snapshot():
    manifest=json.loads((OUT/'manifest.json').read_text());data={}
    for job in manifest['jobs']:
        p=OUT/(job['key']+'.json')
        if p.exists():
            r=json.loads(p.read_text());r.setdefault('checks',[]);r.setdefault('validations',[])
            data[job['key']]=r
    return manifest,data

def scheduler(execution):
    arrays={execution['null_array'],execution['ordered_array']}
    try:
        result=subprocess.run(['squeue','-h','-r','-u',getpass.getuser(),'-o','%F|%K|%T|%r'],capture_output=True,text=True,timeout=15,check=True)
        tasks={}
        for line in result.stdout.splitlines():
            fields=line.split('|',3)
            if len(fields)==4 and fields[0] in arrays:tasks[fields[0]+'_'+fields[1]]=dict(state=fields[2],reason=fields[3])
        return dict(checked=now(),available=True,tasks=tasks,active=len(tasks),error=None)
    except (OSError,subprocess.SubprocessError) as exc:return dict(checked=now(),available=False,tasks={},active=None,error=str(exc))

def state(job,r,monitor,execution):
    if r and r.get('complete'):return 'Complete'
    array=execution['ordered_array'] if job['ordered'] else execution['null_array']
    task=monitor.get('tasks',{}).get(array+'_'+str(job['index']%36))
    stage=(r or {}).get('stage','no checkpoint')
    if task:return task['state'].capitalize()+' · '+stage
    if monitor.get('available'):return 'Stopped / incomplete · '+stage
    return 'Scheduler unknown · '+stage

def assessment(r):
    if not r.get('complete'):return 'Provisional; validation unfinished'
    expected=1 if r['family'] in ('gaussian','memory') else 3
    validations=r.get('validations',[])
    if len(validations)!=expected or any(not v.get('numerically_resolved',False) for v in validations):return 'Numerical validation unresolved'
    if r['study'].startswith('Jemaa'):return 'Speed only: underdetermined'
    if r['study'].startswith('Hosny'):return 'Descriptive quartile fit'
    return 'All fitted + validation intervals met' if r.get('all_reading_intervals') and all(all(o.get('within') for o in v['rows']) for v in validations) else 'At least one fitted/validation interval missed'

def audit_sharing(r):
    a,b=r['constraint_audit']['group_parameters'];times=['tp'] if r['family']=='gaussian' else ['tp','tm']
    shared=all(a[k]==b[k] for k in a if k not in times)
    bounded=all(0<a[k]<12 and 0<b[k]<12 for k in times)
    ordered=all(a[k]<b[k] for k in times) if r['hypothesis']=='ordered_times' else a==b
    return shared and bounded and ordered

def publish(monitor=None,watcher_state='One-time snapshot'):
    with (OUT/'.publish.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        manifest,data=snapshot();execution=json.loads((OUT/'execution.json').read_text())
        if monitor is None:monitor=scheduler(execution)
        generated=now();studies=list(dict.fromkeys(j['study'] for j in manifest['jobs']));jobmap={j['key']:j for j in manifest['jobs']}
        from figures import render
        plotdata={k:dict(r,validations=[v for v in r['validations'] if v.get('numerically_resolved') and 'rows' in v]) for k,r in data.items()}
        render(plotdata,studies)
        counts=dict(complete=sum(r['complete'] for r in data.values()),provisional=sum(not r['complete'] for r in data.values()),awaiting_checkpoint=len(manifest['jobs'])-len(data))
        parameters=[];predictions=[];comparison=[];progress=[];sections=[]
        for job in manifest['jobs']:
            r=data.get(job['key']);progress.append(dict(key=job['key'],study=job['study'],family=job['family'],ordered=job['ordered'],state=state(job,r,monitor,execution),stage=r['stage'] if r else '',complete=r['complete'] if r else False,starts_completed=r.get('starts_completed',0) if r else 0,checkpoint_time=r['generated'] if r else '',score=r['score'] if r else None))
        for key,r in data.items():
            a,b=r['constraint_audit']['group_parameters'];st=state(jobmap[key],r,monitor,execution)
            for name in a:
                if name=='ps' or (name in ('tm','r') and r['family']=='gaussian'):continue
                if name in ('vs','self_a','self_K') and name not in r['parameters']:continue
                if name=='force' and not r['study'].startswith('Hosny'):continue
                parameters.append(dict(study=r['study'],family=r['family'],hypothesis=r['hypothesis'],state=st,stage=r['stage'],complete=r['complete'],parameter=name,reference=a[name],wgd_plus=b[name],ratio=b[name]/a[name] if a[name] else None,shared=a[name]==b[name],fit_score=r['score'],checkpoint_time=r['generated'],result_file=key+'.json'))
            for role in ('rows','checks'):
                for o in r[role]:predictions.append(dict(study=r['study'],family=r['family'],hypothesis=r['hypothesis'],stage=r['stage'],complete=r['complete'],role='fit' if role=='rows' else 'excluded_check',group=o['group'],metric=o['metric'],quantile=o.get('quantile'),observed=o['obs'],low=o.get('low'),high=o.get('high'),predicted=o['predicted'],checkpoint_time=r['generated']))
        for study in studies:
            overview=[];times=[];shared=[];details=[];nulls=[]
            for family in FAMILIES:
                key=family+'_ordered_'+study;nk=family+'_shared_'+study;r=data.get(key);null=data.get(nk)
                nulls.append([family,state(jobmap[nk],null,monitor,execution),fmt(null['score']) if null else '—',null.get('starts_completed',0) if null else 0,null['generated'] if null else '—'])
                if not r:
                    overview.append([family,state(jobmap[key],None,monitor,execution),'Awaiting checkpoint','—',fmt(null['score']) if null else '—','—','—']);continue
                validations=[v for v in r['validations'] if v.get('numerically_resolved') and 'rows' in v]
                st=state(jobmap[key],r,monitor,execution);desc=assessment(r)
                checks=(f"{sum(bool(o['within']) for o in r['checks'])}/{len(r['checks'])} reading intervals" if r['checks'] and not study.startswith('Hosny') else 'see details' if r['checks'] else 'pending' if not r['complete'] else 'unavailable' if 'checks_unresolved' in r else 'none available')
                overview.append([family,st,desc,fmt(r['score']),fmt(null['score']) if null else '—',checks,r['numerically_unresolved_evaluations']])
                a,b=r['constraint_audit']['group_parameters'];times.append([family,r['stage'],fmt(a['tp']),fmt(b['tp']),fmt(b['tp']/a['tp']),fmt(a['tm']) if family!='gaussian' else '—',fmt(b['tm']) if family!='gaussian' else '—',fmt(b['tm']/a['tm']) if family!='gaussian' else '—'])
                shared.append([family]+[fmt(r['parameters'][k]) if k in r['parameters'] else '1 (fixed)' if k=='eta' else '—' for k in ['r','eta','vs','self_a','self_K','force']]+[fmt(r['shared_amplitude'])])
                comparison.append(dict(study=study,family=family,state=st,stage=r['stage'],complete=r['complete'],score=r['score'],null_score=null['score'] if null else None,null_complete=null['complete'] if null else None,assessment=desc,tp_ratio=r['ratios']['tp'],tm_ratio=r['ratios'].get('tm'),unresolved_evaluations=r['numerically_unresolved_evaluations'],validation_min=min((v['score'] for v in validations),default=None),validation_max=max((v['score'] for v in validations),default=None),constraint_check=audit_sharing(r),checkpoint_time=r['generated']))
                observed=[]
                for i,o in enumerate(r['rows']):
                    label=o['metric']+(f" Q{[25,50,75][o['quantile']]}" if 'quantile' in o else '')
                    intervals=f"{o['low']:g}–{o['high']:g}" if 'low' in o else f"IQR scale {o['width']:g}"
                    values='; '.join(fmt(v['rows'][i]['predicted']) for v in validations) or 'pending / unavailable'
                    observed.append(['WGD+' if o['group'] else 'WGD−',label,fmt(o['obs']),intervals,fmt(o['predicted']),values])
                excluded=[['WGD+' if o['group'] else 'WGD−',o['metric']+(f" Q{[25,50,75][o['quantile']]}" if 'quantile' in o else ''),fmt(o['obs']),fmt(o['predicted']),o.get('within')] for o in r['checks']]
                validation_status=[['seed '+str(v['seed']),v['max_dt'],v['burn'],v.get('numerically_resolved'),fmt(v.get('score')),v.get('reason','')] for v in r['validations']]
                details.append(f'<details><summary>{esc(family)} — {esc(st)} — score {fmt(r["score"])}</summary><p>Checkpoint: {esc(r["generated"])}. Search starts completed: {r["starts_completed"]}/8. Search starts reporting convergence: {sum(s["success"] for s in r["starts"])}. Simulation ensemble: {r["sample_n"]}. Independently checked sharing/order constraints: {audit_sharing(r)}.</p>'
                    +table(['Group','Measurement','Observed','Interval / scale','Current prediction','Validation predictions'],observed)
                    +'<h3>Numerical validation</h3>'+table(['Seed','Maximum dt (h)','Burn-in (h)','Resolved','Score','Note'],validation_status)
                    +'<h3>Excluded measurements</h3>'+table(['Group','Metric','Observed','Predicted','Within reading interval'],excluded)
                    +('<p>Excluded checks unresolved: '+esc(r['checks_unresolved'])+'</p>' if 'checks_unresolved' in r else '')
                    +f'<p><a href="{key}.json">Full checkpoint, parameters and starts</a></p></details>')
            sections.append('<section class="study" data-study="'+esc(study)+'"><h2>'+esc(study)+'</h2>'
                +table(['Variant','Progress','Assessment','Current score','Shared-null score','Excluded checks','Unresolved evaluations'],overview)
                +'<h3>Current persistence times (hours)</h3>'+table(['Variant','Stage','tp WGD−','tp WGD+','tp ratio','tm WGD−','tm WGD+','tm ratio'],times)
                +'<h3>Parameters shared within each pair</h3>'+table(['Variant','Feedback r','Noise eta','Saturation vs','Self-activation a','Self-activation K','Force','Readout scale'],shared)
                +(f'<p><a href="{study}.png"><img src="{study}.png" alt="Observed and predicted {esc(study)}"></a> <a href="{study}.pdf">PDF</a></p>' if (OUT/(study+'.png')).exists() else '')
                +''.join(details)+'<details><summary>Shared-null benchmarks</summary>'+table(['Variant','Progress','Score','Starts completed','Checkpoint time'],nulls)+'</details></section>')
        exports=[('comparison.csv',comparison,['study','family','state','stage','complete','score','null_score','null_complete','assessment','tp_ratio','tm_ratio','unresolved_evaluations','validation_min','validation_max','constraint_check','checkpoint_time']),
                 ('predictions.csv',predictions,['study','family','hypothesis','stage','complete','role','group','metric','quantile','observed','low','high','predicted','checkpoint_time']),
                 ('current_best_parameters.csv',parameters,['study','family','hypothesis','state','stage','complete','parameter','reference','wgd_plus','ratio','shared','fit_score','checkpoint_time','result_file']),
                 ('progress.csv',progress,['key','study','family','ordered','state','stage','complete','starts_completed','checkpoint_time','score'])]
        for filename,rows,fields in exports:write_csv(filename,fields,rows)
        latest=max((r['generated'] for r in data.values()),default='none')
        monitor_record=dict(generated=generated,latest_checkpoint=latest,watcher_state=watcher_state,publisher_job=os.environ.get('SLURM_JOB_ID'),counts=counts,scheduler=monitor,all_complete=counts['complete']==len(manifest['jobs']))
        atomic('report_status.json',json.dumps(monitor_record,indent=2)+'\n')
        document='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>WGD: ordered persistence-time fits</title><style>body{font:16px system-ui;max-width:1450px;margin:30px auto;padding:0 20px;color:#172735}table{border-collapse:collapse;font-size:13px;margin:12px 0}th,td{border:1px solid #ccd5df;padding:8px;text-align:left}th{background:#eef3f7}details{margin:15px 0}summary{cursor:pointer;font-weight:600}img{max-width:100%}select{font:inherit;padding:8px}a{color:#075ca4}.scroll{overflow-x:auto}.status{padding:14px;background:#eef3f7;border-left:4px solid #2171a5}</style><h1>WGD fits: both persistence times increase</h1>'''
        producer=('Running or queued: '+str(monitor['active'])+' tasks') if monitor['available'] and monitor['active'] else 'No fitting tasks active' if monitor['available'] else 'Scheduler status unavailable'
        document+=f'<div class="status"><b>{counts["complete"]}/{len(manifest["jobs"])} fits complete</b> · {counts["provisional"]} provisional checkpoints · {counts["awaiting_checkpoint"]} awaiting a checkpoint.<br>Report generated: {generated}. Latest fit checkpoint: {latest}.<br>Fit producer: {esc(producer)}. Publisher: {esc(watcher_state)}. Scheduler checked: {monitor["checked"]}.</div>'
        document+='''<p>The server checks for new checkpoints every 30 seconds and updates this report, figures and CSVs. Open or manually refresh to load the latest snapshot; the browser does not refresh itself. Initial and search checkpoints are provisional. Scores from different optimization stages can use different simulation ensemble sizes.</p>
<p>All times are hours. Ordered alternatives enforce 0 &lt; tp− &lt; tp+ &lt; 12 and, where present, 0 &lt; tm− &lt; tm+ &lt; 12. All non-time parameters are shared within each experiment. Equal-time nulls are separate benchmarks. The numerical floor is 1e-6 h. Near-equality is allowed, so ratios close to one need not imply a detectable effect.</p>
<p><a href="../README.md">Methods</a> · <a href="comparison.csv">Comparison CSV</a> · <a href="current_best_parameters.csv">Current parameters CSV</a> · <a href="predictions.csv">Predictions CSV</a> · <a href="progress.csv">All task progress</a> · <a href="manifest.json">Manifest</a> · <a href="report_status.json">Publisher status</a> · <a href="../../persistence_only/results/report.html">Previous unordered-fit report</a></p>
<p>Losses are descriptive, not likelihoods. Compare variants within an experiment. Reading intervals are not confidence intervals. Numerical resource limits are reported as unresolved evaluations, not biological infeasibility. A completed search is assessed only after its expected validation checks are present and numerically resolved. Jemaà speed-only fits remain underdetermined; Hosny has no statistical adequacy threshold. Excluded measurements are correlated checks from the same experiments.</p>'''
        for figure in ('comparison','persistence_ratios'):
            document+=f'<p><a href="{figure}.png"><img src="{figure}.png" alt="{figure.replace("_"," ")}"></a> <a href="{figure}.pdf">PDF</a></p>'
        document+='<label>Experiment <select id="study"><option value="all">All experiments</option>'+''.join('<option>'+esc(s)+'</option>' for s in studies)+'</select></label>'+''.join(sections)
        document+='''<script>document.getElementById('study').addEventListener('change',e=>document.querySelectorAll('.study').forEach(s=>s.hidden=e.target.value!=='all'&&s.dataset.study!==e.target.value));</script></html>'''
        atomic('report.html',document)
        print(json.dumps(dict(generated=generated,counts=counts,watcher=watcher_state,active=monitor['active'])),flush=True)
        return monitor_record

def stop(signum,frame):
    global STOP
    STOP=True

def watch(interval=30):
    signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
    while True:
        execution=json.loads((OUT/'execution.json').read_text());monitor=scheduler(execution)
        _,data=snapshot();complete=len(data)==72 and all(r['complete'] for r in data.values())
        terminal=complete or (monitor['available'] and monitor['active']==0)
        label='Stopped: all fits completed' if complete else 'Stopped: fitting jobs no longer active' if terminal else 'Stopped by signal' if STOP else 'Running; 30-second checkpoint polling'
        publish(monitor,label)
        if terminal or STOP:return
        time.sleep(interval)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--watch',action='store_true');args=parser.parse_args()
    if args.watch:watch()
    else:publish()
