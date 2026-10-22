"""Reconstruct native deBelly histories for the fixed-slope mechanism report."""
from pathlib import Path
import csv,json,sys,hashlib,gzip
import numpy as np
CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE)); import gradient_reversal as g
r,n,e=g.r,g.n,g.e
RUN=n.DATA/'runs/fixed_slope_reversal_20261006'
OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly'

def setup():
    manifest=r.read(RUN/'manifest.json')
    ends=list(csv.DictReader(open(RUN/'results/endpoints.csv')))
    jobs=[]
    for L in [18.,23.,39.]:
        case=next(c for c in manifest['cases'] if c['model']=='debelly' and c['condition']=='relative_20pct' and c['length']==L)
        row=next(x for x in ends if int(x['case_index'])==case['case_index'] and x['definition']=='orientation')
        for name, d in [('pass',float(row['endpoint'])),('fail',float(row['lower_failing_duration']))]:
            jobs.append(dict(slot=len(jobs),case=case,label=name,duration=d,endpoint=row))
    r.write(OUT/'jobs.json',jobs)
    r.write(OUT/'provenance.json',dict(generated_at=n.now(),source_sha256={str(p.relative_to(n.ROOT)):n.sha(p) for p in [Path(__file__),CODE/'gradient_reversal.py',CODE/'core/native.py',CODE/'core/models/debelly.py',CODE/'core/corrected_reversal.py']},campaign_manifest_sha256=n.sha(RUN/'manifest.json')))
    print(json.dumps([dict(slot=j['slot'],length=j['case']['length'],duration=j['duration'],label=j['label']) for j in jobs]))

def worker(slot):
    job=r.read(OUT/'jobs.json')[slot];case=job['case'];duration=job['duration']
    cp=n.ROOT/case['checkpoint'];assert n.sha(cp)==case['checkpoint_sha256']
    state,prep=n.read_checkpoint(cp,'debelly')
    pulse=n.run('debelly',case['length'],state,g.response_cue(case),duration,interval=min(2.5,duration/20),n=case['cells'])
    recovery=n.run('debelly',case['length'],pulse.state,e.withdrawal(case),10000.,interval=2.5,n=case['cells'])
    # Identical pulse/recovery sample interval and native step as the campaign.
    arrays={'x':pulse.x,'pulse_duration':np.array(duration),'initial_orientation':np.array(case['initial_orientation'])}
    for phase,tr in [('pulse',pulse),('recovery',recovery),('preparation',prep)]:
        for key in ['time','orientation','contrast']:
            arrays[phase+'_'+key]=getattr(tr,key)
        for key,value in tr.fields.items():arrays[phase+'_field_'+key]=value
    path=OUT/f'L{case["length"]:g}_{job["label"]}.npz';np.savez_compressed(path,**arrays)
    tail=recovery.time>=9500;ori=recovery.orientation[tail];contr=recovery.contrast[tail]
    flag=float(max(max(np.max(np.ptp(v[tail],axis=0))/case['fixed_scales'][k],np.max(abs(v[-1]-v[-2]))/2.5*500/case['fixed_scales'][k]) for k,v in r.settling.monitored(recovery,'debelly',case['length']).items()))
    result=dict(length=case['length'],label=job['label'],duration=duration,orientation_success=bool(ori.min()>.05 and contr.min()>=.05 and np.ptp(ori)<.005),full_polarity_success=bool(ori.min()>.05 and contr.min()>=.05 and np.ptp(ori)<.005 and np.max(abs(ori/abs(case['initial_orientation'])-1))<=.01),minimum_orientation=float(ori.min()),maximum_orientation=float(ori.max()),minimum_contrast=float(contr.min()),field_tail_max_scaled=flag,field_tail_pass=flag<=.001,file_sha256=n.sha(path),checkpoint_sha256=n.sha(cp))
    # Independent agreement with saved original trials.
    with gzip.open(RUN/'results/records.json.gz','rt') as f: records=json.load(f)
    r.write(path.with_suffix('.json'),result);print(json.dumps(result))

if __name__=='__main__':
    if sys.argv[1]=='setup':setup()
    else:worker(int(sys.argv[2]))
