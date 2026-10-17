"""Supplement stage-1 audit by separating approximately monotone end-polarized fields."""
from pathlib import Path
import csv,json,gzip
from collections import defaultdict
import numpy as np
import stage1_endpoint_audit as a

OUT=a.OUT
scales={(r['model'],r['length']):r['fixed_scales'][a.FIELD.get(r['model'],'active')] for r in json.loads((a.CODE/'reference/continuation_20261001/fixed_scales.json').read_text()) if r['event']=='release'}
groups=defaultdict(list)
for campaign,(work,results) in a.CAMPAIGNS.items():
 with gzip.open(results/'records.json.gz','rt') as f:records=json.load(f)
 tasks={c['task_index']:c for c in records['cases']}
 for job in json.loads((work/'response_manifest.json').read_text()):
  task=tasks[job['task_index']]
  if task['assay']!='finite_countercue':continue
  p=work/'trials'/f"chunk_{job['slot']:05d}.json"
  row=json.loads(p.read_text())
  with np.load(p.with_suffix('.npz')) as z:
   ys=z['recovery_final_active'];os=z['recovery_tail_orientation']
   for i,t in enumerate(row['tests']):
    if not t['settled']:continue
    o=os[i]
    if np.min(o)>.05:y=ys[i];side='right'
    elif np.max(o)<-.05:y=ys[i][::-1];side='left'
    else:continue
    key=(task['model'],task['length']);scale=scales[key]
    drawdown=float(np.max(np.maximum.accumulate(y)-y)/scale)
    groups[key].append(dict(y=y.copy(),drawdown=drawdown,campaign=campaign,duration=t['duration'],side=side,orientation=t['final_orientation'],source=str(p.with_suffix('.npz').relative_to(a.ROOT)),row=i,orientation_range=t['orientation_range_last_500']))
rows=[];witnesses={}
for key in sorted(groups,key=lambda k:(a.MODELS.index(k[0]),k[1])):
 model,L=key;items=groups[key];monotone=[t for t in items if t['drawdown']<=.01]
 for comparison,sel in [('all',monotone),('right_only',[t for t in monotone if t['side']=='right'])]:
  diameter=None;pair=None
  if len(sel)>=2:
   m=np.stack([t['y'] for t in sel]);env=np.ptp(m,axis=0);col=int(np.argmax(env));diameter=float(env[col]/scales[key]);pair=(int(np.argmin(m[:,col])),int(np.argmax(m[:,col])))
  rows.append(dict(model=model,length=L,comparison=comparison,total_gate_passing=len(items),approximately_monotone=len(sel),nonmonotone=sum(t['drawdown']>.01 for t in items),monotone_scaled_diameter=diameter,monotone_above_1_percent=diameter is not None and diameter>.01))
  if pair and comparison=='all':
   witnesses[f'{model}_L{L:g}']=dict(diameter_scaled=diameter,first={k:v for k,v in sel[pair[0]].items() if k!='y'},second={k:v for k,v in sel[pair[1]].items() if k!='y'})
a.csvwrite(OUT/'monotone_endpoint_diameters.csv',rows)
a.write(OUT/'monotone_witnesses.json',witnesses)
summary=[]
for model in a.MODELS:
 r=[t for t in rows if t['model']==model and t['comparison']=='all']
 summary.append(dict(model=model,lengths_with_monotone_comparison=sum(t['approximately_monotone']>=2 for t in r),lengths_monotone_above_1_percent=sum(t['monotone_above_1_percent'] for t in r),maximum_monotone_diameter=max((t['monotone_scaled_diameter'] for t in r if t['monotone_scaled_diameter'] is not None),default=None),lengths_with_nonmonotone_endpoints=sum(t['nonmonotone']>0 for t in r),gate_passing_nonmonotone_endpoints=sum(t['nonmonotone'] for t in r)))
a.csvwrite(OUT/'monotone_model_summary.csv',summary)
a.write(OUT/'monotone_summary.json',dict(generated_at=a.now(),source_sha256=a.sha(Path(__file__).resolve()),classification='After orientation alignment, maximum backward excursion from running maximum <=1% of existing fixed field scale. A morphological screen, not a count of fronts or stability test.',models=summary))
print(json.dumps(summary),flush=True)
