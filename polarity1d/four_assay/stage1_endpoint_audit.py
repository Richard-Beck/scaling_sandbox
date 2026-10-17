"""Read-only comparison of existing cue-free finite-pulse recovery endpoints."""
from pathlib import Path
import csv, gzip, hashlib, json
from collections import Counter, defaultdict
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[1]
OUT = CODE / 'pilots' / 'stage1_endpoint_audit_20261002'
CAMPAIGNS = {
    'Ax_over_L': (CODE / 'work', CODE / 'results'),
    'Ax': (CODE / 'pilots/physical_cue_20261002/work', CODE / 'pilots/physical_cue_20261002/results'),
}
MODELS = ['wave_pinning', 'otsuji', 'goryachev', 'legi', 'debelly', 'holmes_model3', 'spring']
LABELS = ['Wave-pinning', 'Otsuji', 'Goryachev', 'LEGI', 'deBelly', 'Holmes 3', 'Spring']
FIELD = {'legi':'response', 'debelly':'rac', 'holmes_model3':'rac_active'}
THRESHOLDS = [.001, .01, .05]

def now(): return datetime.now(timezone.utc).isoformat()
def write(path, obj): path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576), b''): digest.update(block)
    return digest.hexdigest()
def csvwrite(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def main():
    OUT.mkdir(exist_ok=False)
    write(OUT/'status.json',dict(generated_at=now(),complete=False,stage='reading_existing_diagnostics'))
    inputs = {}
    def track(path):
        inputs[str(path.relative_to(ROOT))] = sha(path)
        return path
    track(Path(__file__).resolve())
    groups = defaultdict(list)
    counters = {}
    source_records = {}
    plans = {}
    for campaign, (work, results) in CAMPAIGNS.items():
        status=json.loads(track(results/'status.json').read_text())
        assert status['complete']
        with gzip.open(track(results/'records.json.gz'),'rt') as f: records=json.load(f)
        source_records[campaign]=records
        plan=json.loads(track(work/'response_manifest.json').read_text()); plans[campaign]=plan
        tasks={c['task_index']:c for c in records['cases']}
        expected={c['task_index']:c['tests'] for c in records['cases'] if c['assay']=='finite_countercue' and c['eligible']}
        seen=defaultdict(list)
        count=Counter()
        for job in plan:
            task=tasks[job['task_index']]
            if task['assay']!='finite_countercue': continue
            path=work/'trials'/f"chunk_{job['slot']:05d}.json"
            row=json.loads(track(path).read_text())
            with np.load(track(path.with_suffix('.npz'))) as z:
                assert set(z.files)=={'durations','pulse_final_active','recovery_final_active','recovery_tail_orientation'}, z.files
                profiles=z['recovery_final_active']; tails=z['recovery_tail_orientation']; durations=z['durations']
                assert len(row['tests'])==len(durations)==len(profiles)==len(tails)
                for i,t in enumerate(row['tests']):
                    assert float(durations[i])==t['duration']
                    seen[task['task_index']].append(t)
                    count['trials']+=1
                    ori=tails[i]
                    np.testing.assert_allclose(ori[-1],t['final_orientation'],rtol=1e-10,atol=1e-12)
                    np.testing.assert_allclose(np.ptp(ori),t['orientation_range_last_500'],rtol=1e-10,atol=1e-12)
                    stable=bool(np.ptp(ori)<.005 and t['minimum_contrast_last_500']>=.05)
                    assert stable==t['settled']
                    right=bool(stable and np.min(ori)>.05)
                    left=bool(stable and np.max(ori)<-.05)
                    assert right==t['success']
                    if not (left or right):
                        count['excluded_nonpolar_or_failed_gate']+=1
                        continue
                    side='right' if right else 'left'
                    count['accepted_'+side]+=1
                    y=np.asarray(profiles[i],float)
                    assert np.all(np.isfinite(y))
                    groups[(task['model'],task['length'])].append(dict(campaign=campaign,duration=t['duration'],side=side,orientation=t['final_orientation'],orientation_range=t['orientation_range_last_500'],minimum_contrast=t['minimum_contrast_last_500'],profile=y,aligned=y if right else y[::-1],source=str(path.with_suffix('.npz').relative_to(ROOT)),row=i))
        assert set(seen)==set(expected)
        for idx,tests in seen.items(): assert sorted(tests,key=lambda t:t['duration'])==expected[idx]
        counters[campaign]=dict(count)
        print(campaign,counters[campaign],flush=True)
    assert source_records['Ax']['preparations']==source_records['Ax_over_L']['preparations']
    scales_path=track(CODE/'reference/continuation_20261001/fixed_scales.json')
    scales={(r['model'],r['length']):r['fixed_scales'][FIELD.get(r['model'],'active')] for r in json.loads(scales_path.read_text()) if r['event']=='release'}
    preps={ (p['model'],p['length']):p for p in source_records['Ax_over_L']['preparations'] if p['assay']=='finite_countercue'}
    references={}
    for key,p in preps.items():
        cp=track(ROOT/p['checkpoint'])
        assert inputs[str(cp.relative_to(ROOT))]==p['checkpoint_sha256']
        if not p['eligible']:continue
        with np.load(cp) as z:
            y=z['active'][-1].copy(); orientation=float(z['orientation'][-1])
        assert orientation<-.05
        references[key]=dict(campaign='release_reference',duration=None,side='left',orientation=orientation,orientation_range=None,minimum_contrast=None,profile=y,aligned=y[::-1],source=str(cp.relative_to(ROOT)),row=None)
    rows=[]; witnesses={}; examples={}
    def diameter(items, field, scale):
        if len(items)<2: return 0.,None
        matrix=np.stack([r[field] for r in items])
        envelope=np.ptp(matrix,axis=0); col=int(np.argmax(envelope))
        i=int(np.argmin(matrix[:,col])); j=int(np.argmax(matrix[:,col]))
        return float(envelope[col]/scale),(items[i],items[j],col)
    def metadata(item): return {k:v for k,v in item.items() if k not in ['profile','aligned']}
    for key in sorted(preps,key=lambda k:(MODELS.index(k[0]),k[1])):
        model,L=key; items=groups[key]; scale=scales[key]
        for campaign in ['Ax_over_L','Ax','pooled']:
            selected=[x for x in items if campaign=='pooled' or x['campaign']==campaign]
            for side in ['left','right','aligned']:
                subset=[x for x in selected if side=='aligned' or x['side']==side]
                field='aligned' if side=='aligned' else 'profile'
                diam,pair=diameter(subset,field,scale)
                ref=references.get(key)
                vs_ref=max((float(np.max(abs(x['aligned']-ref['aligned']))/scale) for x in subset),default=0.) if ref else None
                row=dict(model=model,length=L,campaign=campaign,comparison=side,n_endpoints=len(subset),fixed_field_scale=scale,diameter_scaled=diam,maximum_difference_from_aligned_reference=vs_ref,above_0_1_percent=diam>.001,above_1_percent=diam>.01,above_5_percent=diam>.05,maximum_orientation_range=max((x['orientation_range'] for x in subset),default=None))
                rows.append(row)
                if pair and campaign=='pooled' and side=='aligned':
                    a,b,col=pair
                    wid=f'{model}_L{L:g}'
                    witnesses[wid]=dict(diameter_scaled=diam,coordinate_index=col,first=metadata(a),second=metadata(b),scale=scale)
                    examples[key]=(a,b,ref,diam)
    csvwrite(OUT/'endpoint_diameters.csv',rows)
    write(OUT/'witnesses.json',witnesses)
    summary=[]
    for model in MODELS:
        selected=[r for r in rows if r['model']==model and r['campaign']=='pooled' and r['comparison']=='aligned' and r['n_endpoints']>=2]
        right=[r for r in rows if r['model']==model and r['campaign']=='pooled' and r['comparison']=='right' and r['n_endpoints']>=2]
        rec=dict(model=model,lengths_with_comparable_endpoints=len(selected),lengths_above_0_1_percent=sum(r['above_0_1_percent'] for r in selected),lengths_above_1_percent=sum(r['above_1_percent'] for r in selected),lengths_above_5_percent=sum(r['above_5_percent'] for r in selected),maximum_scaled_diameter=max((r['diameter_scaled'] for r in selected),default=None),right_only_lengths_above_1_percent=sum(r['above_1_percent'] for r in right),right_only_maximum_scaled_diameter=max((r['diameter_scaled'] for r in right),default=None),maximum_difference_from_reference=max((r['maximum_difference_from_aligned_reference'] for r in selected if r['maximum_difference_from_aligned_reference'] is not None),default=None))
        summary.append(rec)
    csvwrite(OUT/'model_summary.csv',summary)
    # Plot the largest between-endpoint difference at any length in each model.
    fig,axes=plt.subplots(7,2,figsize=(12,19),gridspec_kw={'width_ratios':[1,1.4]})
    saved={}
    for i,(model,label) in enumerate(zip(MODELS,LABELS)):
        ax,profileax=axes[i]
        for campaign,color in [('Ax_over_L','#176d9c'),('Ax','#db7b24'),('pooled','#27804b')]:
            sel=[r for r in rows if r['model']==model and r['campaign']==campaign and r['comparison']=='aligned' and r['n_endpoints']>=2]
            if sel: ax.plot([r['length'] for r in sel],[100*r['diameter_scaled'] for r in sel],'.-',label=campaign,color=color)
        ax.axhline(1,color='gray',ls='--',lw=.8);ax.set_ylabel(label+'\nDiameter (% fixed scale)');ax.set_xlabel('Nominal length (Spring: resting length)')
        choices=[(k,v) for k,v in examples.items() if k[0]==model]
        if not choices:
            profileax.text(.5,.5,'No comparable gate-passing polarized endpoints',ha='center',va='center',transform=profileax.transAxes,wrap=True)
        else:
            key,(a,b,ref,diam)=max(choices,key=lambda kv:kv[1][3])
            for item,color in [(a,'#176d9c'),(b,'#db7b24')]:
                profileax.plot(np.linspace(0,1,len(item['aligned'])),item['aligned'],color=color,label=f"{item['campaign']}, pulse {item['duration']:g}, {item['side']}")
            if ref: profileax.plot(np.linspace(0,1,len(ref['aligned'])),ref['aligned'],'k--',lw=.9,label='Released preparation, reflected')
            profileax.set_title(f'{label}: L={key[1]:g}, largest diameter {100*diam:.4g}%')
            profileax.legend(fontsize=7)
            saved[f'{model}_first']=a['aligned']; saved[f'{model}_second']=b['aligned']
            if ref:saved[f'{model}_reference']=ref['aligned']
        profileax.set_xlabel('Relative coordinate (aligned to right)');profileax.set_ylabel('Native active/response field')
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Stage 1: saved cue-free recovery endpoints across pulse durations\nOrientation/contrast gate passed; active-field stationarity and full-state attraction are untested',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.97))
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'endpoint_profile_audit.{ext}',dpi=160)
    plt.close(fig)
    np.savez_compressed(OUT/'witness_profiles.npz',**saved)
    write(OUT/'input_hashes.json',inputs)
    result=dict(generated_at=now(),complete=True,passed=True,campaign_counts=counters,models=summary,input_file_count=len(inputs),methods=dict(comparison='Exact maximum pairwise L-infinity difference, computed from coordinate-wise envelope; no amplitude normalization',normalization='Existing model/length release fixed field scales',alignment='Reverse reference/material grid for left endpoints; also report same-orientation comparisons',thresholds=THRESHOLDS,reference='Original verified release checkpoint; preparation records identical across campaigns'),limitations=['Saved trial diagnostics contain only the final active/response field and final 500-unit orientation history, not the full final state or active-field time history.','Gate-passing differences identify candidate retained profiles, not proofs of multiple attracting equilibria.','Finite recovery horizons may leave slow transients.','For mechanical models reflection compares chemical fields on reference/material coordinates, not complete mechanically equivalent states.','Gradient endpoints have different maintained inputs and are excluded from the cue-free uniqueness comparison.','Trials sample only states reachable from the existing preparation; agreement cannot prove uniqueness.'])
    write(OUT/'summary.json',result);write(OUT/'status.json',dict(generated_at=now(),complete=True,passed=True,producer='finished'))
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
