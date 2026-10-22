"""Exact Holmes finite-gradient restarts and native-field diagnostic figures.

Use the retained manuscript solver, parameter regime, slope and checkpoints.
Only output sampling is enriched. This script never changes original evidence.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, sys
from dataclasses import replace
from pathlib import Path
import numpy as np

CODE=Path(__file__).resolve().parents[1]
ROOT=CODE.parent
sys.path.insert(0,str(CODE))
sys.path.insert(0,str(CODE/'core'))
import native as n
import gradient_reversal as g
from polarity1d.simulation import simulate
from polarity1d.solvers import laplacian_no_flux
from polarity1d.contract import Grid1D

CAMPAIGN=ROOT/'data/runs/fixed_slope_reversal_20261006'
OUT=ROOT/'data/runs/fixed_slope_mechanisms_20261006/holmes_model3'
FIG=ROOT/'figures/fixed_slope_mechanisms/holmes_model3'
LENGTHS=[10.,28.,50.]

def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,default=n.convert,allow_nan=False)+'\n')

def trials():
    manifest=json.loads((CAMPAIGN/'manifest.json').read_text())
    cases={c['length']:c for c in manifest['cases'] if c['model']=='holmes_model3' and c['condition']=='relative_20pct'}
    rows=list(csv.DictReader((CAMPAIGN/'results/endpoints.csv').open()))
    tasks=[]
    for length in LENGTHS:
        for gate in ['orientation','full_polarity']:
            row=next(r for r in rows if r['model']=='holmes_model3' and r['condition']=='relative_20pct' and float(r['length'])==length and r['definition']==gate)
            for result,col in [('pass','endpoint'),('fail','lower_failing_duration')]:
                duration=float(row[col])
                if any(t['length']==length and t['duration']==duration for t in tasks):continue
                tasks.append(dict(length=length,gate=gate,result=result,duration=duration,case=cases[length],row=row,
                                  tag=f'L{int(length)}_{gate}_{result}'))
    return tasks

def history(tr):
    states=np.stack([tr.fields[k] for k in n.model_for('holmes_model3').state_names],axis=1)
    coord=2*tr.x/(tr.x[0]+tr.x[-1])-1
    ori=states[:,0]@coord/states[:,0].sum(axis=1)
    return states,ori

def worker(slot):
    task=trials()[slot];case=task['case'];length=task['length'];duration=task['duration']
    checkpoint=ROOT/case['checkpoint']
    if n.sha(checkpoint)!=case['checkpoint_sha256']:raise RuntimeError('checkpoint changed')
    state,prep=n.read_checkpoint(checkpoint,'holmes_model3')
    model=n.model_for('holmes_model3')
    kwargs=dict(length=length,cells=case['cells'],backend='bdf',max_step=10.,rtol=2e-6,atol=1e-9)
    pulse=simulate(model,initial_state=state,protocol=g.response_cue(case),end_time=duration,
                   output_interval=min(.025,duration/400),**kwargs)
    p,p_o=history(pulse)
    times=np.unique(np.r_[np.arange(0,20,.025),np.arange(20,100,.25),np.arange(100,10000,2.5),10000.])
    recovery=simulate(model,initial_state=p[-1],protocol=n.ConstantInput(0),times=times,**kwargs)
    r,r_o=history(recovery)
    tail=times>=9500;contrast=np.ptp(r[:,0],axis=1)
    og=bool(np.min(r_o[tail])>.05 and np.ptp(r_o[tail])<.005 and np.min(contrast[tail])>=.05)
    fg=bool(og and np.max(abs(r_o[tail]/abs(case['initial_orientation'])-1))<=.01)
    field_metrics={k:dict(range_scaled=float(np.max(np.ptp(r[tail,i],axis=0)))/case['fixed_scales'][k],
                         projected_drift_500=float(np.max(abs(r[-1,i]-r[-2,i])))/2.5*500/case['fixed_scales'][k])
                   for i,k in enumerate(model.state_names)}
    maximum=max(v for item in field_metrics.values() for v in item.values())
    meta={k:v for k,v in task.items() if k!='case'}
    meta.update(generated_at=n.now(),checkpoint=case['checkpoint'],checkpoint_sha256=n.sha(checkpoint),
                slope=case['gradient_slope'],initial_orientation=case['initial_orientation'],
                orientation_success=og,full_polarity_success=fg,field_tail_pass=maximum<=.001,
                field_tail_metrics=field_metrics,pulse_solver=pulse.metadata,recovery_solver=recovery.metadata,
                maximum_conservation_error=max(np.max(abs(r[:,[0,1]].sum(axis=(1,2))/state[[0,1]].sum()-1)),
                                               np.max(abs(r[:,[2,3]].sum(axis=(1,2))/state[[2,3]].sum()-1))),
                source_hashes=n.sources(),script_sha256=n.sha(__file__))
    OUT.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(OUT/(task['tag']+'.npz'),x=pulse.x,pulse_time=pulse.time,pulse_state=p,
                        pulse_orientation=p_o,recovery_time=times,recovery_state=r,recovery_orientation=r_o,
                        preparation_state=state)
    write(OUT/(task['tag']+'.json'),meta)
    print(task['tag'],duration,'orientation',og,'strength',fg,'field',maximum<=.001,flush=True)

def mixing_trial(length,duration):
    """Labeled transport perturbation: accelerate both inactive pools 100-fold.

    Exact original native state, original cue and reactions. This is not a new
    baseline curve and zero controls expose any preparation mismatch.
    """
    case=next(t['case'] for t in trials() if t['length']==length)
    state,prep=n.read_checkpoint(ROOT/case['checkpoint'],'holmes_model3')
    model=replace(n.model_for('holmes_model3'),inactive_diffusion=5000.)
    kwargs=dict(length=length,cells=201,backend='bdf',max_step=10.,rtol=2e-6,atol=1e-9)
    start=state
    if duration:
        p=simulate(model,initial_state=state,protocol=g.response_cue(case),end_time=duration,output_interval=duration/100,**kwargs)
        start=np.vstack([p.fields[k][-1] for k in model.state_names])
    r=simulate(model,initial_state=start,protocol=n.ConstantInput(0),end_time=10000.,output_interval=2.5,**kwargs)
    states,o=history(r);tail=r.time>=9500;contrast=np.ptp(states[:,0],axis=1)
    og=bool(np.min(o[tail])>.05 and np.ptp(o[tail])<.005 and np.min(contrast[tail])>=.05)
    fg=bool(og and np.max(abs(o[tail]/abs(case['initial_orientation'])-1))<=.01)
    final=states[-1]
    write(OUT/'mixing_control'/f'L{int(length)}_T{duration:g}.json',dict(length=length,duration=duration,
          inactive_diffusion=5000.,native_inactive_diffusion=50.,orientation_success=og,full_polarity_success=fg,
          final_orientation=float(o[-1]),initial_orientation=case['initial_orientation'],
          zero_scaled_excursion=float(np.max([np.max(abs(states[:,i]-state[i]))/case['fixed_scales'][k] for i,k in enumerate(model.state_names)])) if duration==0 else None,
          checkpoint_sha256=n.sha(ROOT/case['checkpoint']),generated_at=n.now()))
    path=OUT/'mixing_control'/f'L{int(length)}_T{duration:g}.npz'
    np.savez_compressed(path,x=r.x,final_state=final,time=r.time,orientation=o)
    print('mixing_control',length,duration,og,fg,float(o[-1]),flush=True)
    return og,fg

def mixing_control(slot):
    durations=[0.,1.,2.,3.,4.,6.,8.,10.,12.]
    return mixing_trial(LENGTHS[slot//len(durations)],durations[slot%len(durations)])

def mixing_refine(slot):
    length=LENGTHS[slot]
    rows=sorted([json.loads(p.read_text()) for p in (OUT/'mixing_control').glob(f'L{int(length)}_T*.json')],key=lambda r:r['duration'])
    hi=min(r['duration'] for r in rows if r['full_polarity_success'])
    lo=max(r['duration'] for r in rows if r['duration']<hi and not r['full_polarity_success'])
    for _ in range(8):
        if (hi-lo)/hi<=.01:break
        middle=(hi+lo)/2
        _,passed=mixing_trial(length,middle)
        if passed:hi=middle
        else:lo=middle
    write(OUT/'mixing_control'/f'bracket_L{int(length)}.json',dict(length=length,lower_failed_duration=lo,shortest_tested_full_strength_duration=hi,
          relative_bracket_width=(hi-lo)/hi,caution='Local diagnostic refinement after coarse successes; this control is not a globally exhaustive native-protocol duration search.'))

def continuation(slot):
    tags=['L28_orientation_pass','L28_orientation_fail','L50_orientation_pass','L50_orientation_fail','L28_full_polarity_fail','L50_full_polarity_fail']
    tag=tags[slot];task=next(t for t in trials() if t['tag']==tag);case=task['case']
    z=np.load(OUT/(tag+'.npz'));state=z['recovery_state'][-1]
    times=np.unique(np.r_[np.arange(0,89500,100.),np.arange(89500,90000,2.5),90000.])
    tr=simulate(n.model_for('holmes_model3'),initial_state=state,protocol=n.ConstantInput(0),length=task['length'],
                cells=201,times=times,backend='bdf',max_step=10.,rtol=2e-6,atol=1e-9)
    states,o=history(tr);tail=times>=89500;contrast=np.ptp(states[:,0],axis=1)
    og=bool(np.min(o[tail])>.05 and np.ptp(o[tail])<.005 and np.min(contrast[tail])>=.05)
    fg=bool(og and np.max(abs(o[tail]/abs(case['initial_orientation'])-1))<=.01)
    path=OUT/'continuation';path.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path/(tag+'.npz'),x=tr.x,time=10000+times,state=states,orientation=o)
    write(path/(tag+'.json'),dict(tag=tag,length=task['length'],duration=task['duration'],final_recovery_time=100000.,
          initial_orientation=case['initial_orientation'],final_orientation=float(o[-1]),orientation_success=og,
          full_polarity_success=fg,generated_at=n.now(),checkpoint_sha256=case['checkpoint_sha256']))
    print('continuation',tag,og,fg,float(o[-1]),flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['list','worker','mixing-control','mixing-refine','continuation','analyze']);p.add_argument('--slot',type=int)
    args=p.parse_args()
    if args.action=='list':
        print(json.dumps([{k:v for k,v in t.items() if k not in ['case','row']} for t in trials()],indent=2))
    elif args.action=='worker':worker(args.slot)
    elif args.action=='mixing-control':mixing_control(args.slot)
    elif args.action=='mixing-refine':mixing_refine(args.slot)
    elif args.action=='continuation':continuation(args.slot)
    elif args.action=='analyze':
        from holmes_figures import main as figures_main
        figures_main()

if __name__=='__main__':main()
