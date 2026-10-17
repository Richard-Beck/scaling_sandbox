"""Audit saved results without refitting and write the campaign findings."""
import json,math,argparse
from datetime import datetime,timezone
import numpy as np
from campaign import OUT,jobs,jobkey,provenance,save_json
from publish import publish,status,atomic

def main(require_complete=False):
    data={};issues=[];diagnostics=[]
    for f,mask,s in jobs():
        key=jobkey(f,mask,s);p=OUT/(key+'.json')
        if not p.exists():continue
        r=json.loads(p.read_text());data[key]=r
        if r['source_hashes']!=provenance():issues.append(key+': changed sources')
        a,b=r['constraint_audit']['group_parameters'];changed={k for k in a if a[k]!=b[k]}
        if not changed<=({'tp','tm'} if mask else set()):issues.append(key+': group sharing violated')
        if not mask and a!=b:issues.append(key+': null not shared')
        if any(k.endswith('_poly') and k not in ('tp_poly','tm_poly') for k in r['names']):issues.append(key+': ineligible split')
        if s.startswith('Hosny'):
            score=sum(((o['predicted']-o['obs'])/o['width'])**2 for o in r['rows'])
        else:score=sum((math.log(o['predicted']/o['obs'])/(math.log(o['high']/o['low'])/2))**2 for o in r['rows'])
        if not np.isclose(score,r['score'],rtol=1e-9,atol=1e-9):issues.append(key+': score recomputation mismatch')
        for k,v in r['parameters'].items():
            lo,hi=r['bounds'][k]
            if not lo*(1-1e-8)<=v<=hi*(1+1e-8):issues.append(key+': parameter outside bounds')
        if mask:
            nullpath=OUT/(jobkey(f,(),s)+'.json')
            null=json.loads(nullpath.read_text())
            if r['score']>null['score']+1e-4:issues.append(key+': alternative worse than its shared null')
        # Changes are normalized by the same descriptive scales used in fitting.
        deltas=[]
        for validation in r['validations']:
            delta=[]
            for o,v in zip(r['rows'],validation['rows']):
                delta.append(abs((v['predicted']-o['predicted'])/o['width']) if s.startswith('Hosny') else abs(math.log(v['predicted']/o['predicted'])/(math.log(o['high']/o['low'])/2)))
            deltas.append(max(delta))
        diagnostics.append(dict(key=key,max_validation_shift_scaled=max(deltas),near_bounds=r['identifiability']['near_bounds'],status=status(r)))
    complete=len(data)==len(jobs())
    audit=dict(generated=datetime.now(timezone.utc).isoformat(),completed=len(data),expected=len(jobs()),complete=complete,issues=issues,passes=not issues and complete,
               constraint_checks_pass=not any('sharing' in x or 'split' in x for x in issues),numerical_diagnostics=diagnostics,
               note='Scaled numerical shifts include finite simulation noise; no probabilistic threshold is imposed.')
    save_json(OUT/'audit.json',audit)
    lines=['# Direct WGD fits: persistence-only differences','',f'Generated {audit["generated"]}.','',f'Completed {len(data)}/{len(jobs())} fits. '+('All result/provenance/constraint audits pass.' if complete and not issues else 'Campaign is still running.' if not complete else 'Audit issues: '+str(issues)), '',
        'Every experiment is fitted independently. Within a pair, only tp (Gaussian) or tp and tm (other variants) can change. No polarity-derived values, length laws or eligibility rules were used.', '',
        '[Interactive report](report.html) · [Model comparison](comparison.csv) · [Predictions](predictions.csv) · [Numerical audit](audit.json)', '',
        '## Xuan and Mallin', '', '| Experiment | Variant | Fit score | Null score | Assessment |', '|---|---|---:|---:|---|']
    for s in ('Xuan2018','Mallin2023'):
        for f in ('gaussian','memory','coupled','self_activation'):
            r=data.get(jobkey(f,(0,),s));n=data.get(jobkey(f,(),s))
            if r:lines.append(f"| {s} | {f} | {r['score']:.3g} | {n['score']:.3g} | {status(r)} |")
    lines+=['', 'Passing the fitted intervals does not imply agreement with excluded measurements. Early Xuan MSD and Mallin endpoint are displayed as separate checks in the report. These checks use a different simulation seed and fixed fitted amplitude.', '', '## Hosny pairs', '', '| Experiment | Lowest validation score | Fit score | Independent validation score range |', '|---|---|---:|---:|']
    studies=[s for s in dict.fromkeys(s for _,_,s in jobs()) if s.startswith('Hosny')]
    for s in studies:
        rr=[data[jobkey(f,(0,),s)] for f in ('gaussian','memory','coupled','self_activation') if jobkey(f,(0,),s) in data]
        if rr:
            r=min(rr,key=lambda r:np.mean([v['score'] for v in r['validations']]))
            values=[v['score'] for v in r['validations']]
            lines.append(f"| {s} | {r['family']} | {r['score']:.3g} | {min(values):.3g}–{max(values):.3g} |")
    lines+=['', '## Excluded-measurement checks', '', '| Experiment | Variant | Reading intervals met | Missed checks |', '|---|---|---:|---|']
    for s in ('Xuan2018','Mallin2023'):
        for f in ('gaussian','memory','coupled','self_activation'):
            r=data.get(jobkey(f,(0,),s))
            if not r:continue
            cc=r['checks'];miss=[]
            for o in cc:
                if not o['within']:
                    group='WGD+' if o['group'] else 'WGD−'
                    miss.append(f"{group} {o['metric']}: {o['predicted']:.3g} versus {o['low']:g}–{o['high']:g}")
            lines.append(f"| {s} | {f} | {sum(o['within'] for o in cc)}/{len(cc)} | {'; '.join(miss) or 'none'} |")
    lines+=['', 'Hosny comparisons are descriptive and close scores do not establish a preferred mechanism. Models are optimized on path, directness and parallel quartiles; endpoints and perpendicular motion remain excluded checks. Separate treatment fits do not identify knockdown mechanisms.', '', '## What is identified', '',
        'Jemaà supplies only two speed measurements per pair. All models have more free parameters than these measurements can identify. The report therefore labels every Jemaà result underdetermined even when speeds fit well. Other nonlinear fits also have more free parameters than the available Xuan/Mallin summaries.', '',
        'Time ratios are best-fit descriptions, not uniquely identified WGD effects. The report gives local numerical sensitivity rank, near-bound parameters, and variation across near-best multistart solutions. These are not confidence intervals. Finite bounds, local optimization and simulation sensitivity preclude claims of global exclusion.', '',
        'Gaussian OU with a shared noise scale links increased persistence to increased sampled speed and MSD. Xuan asks for lower speed but larger late MSD, and Mallin for shorter paths but larger directness, making those contrasts restrictive for the single-time Gaussian variant.', '',
        '## Validation', '',f'{len(data)} saved results checked for exact parameter sharing, bounds, data/code hashes, and recomputed objectives. Nonlinear fits use independent-seed simulations, half-step integration and doubled burn-in. These checks test numerical sensitivity; no statistical confidence or general stationarity guarantee is inferred.', '']
    atomic('FINDINGS.md','\n'.join(lines));publish()
    print(json.dumps({k:audit[k] for k in ('completed','expected','complete','passes','issues')},indent=2))
    if issues or (require_complete and not complete):raise SystemExit(1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--require-complete',action='store_true');main(p.parse_args().require_complete)
