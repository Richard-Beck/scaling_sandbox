"""Native fixed-slope Wave-pinning reconstruction and bounded mechanism controls."""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, sys
import numpy as np

MANUSCRIPT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MANUSCRIPT/'code'))
import gradient_reversal as g
n, r, e = g.n, g.r, g.e
from polarity1d.solvers import laplacian_no_flux

OUT = MANUSCRIPT/'data/runs/fixed_slope_mechanisms_20261006/wave_pinning'
CAMPAIGN = MANUSCRIPT/'data/runs/fixed_slope_reversal_20261006'

def load_records():
    return json.load(gzip.open(CAMPAIGN/'results/records.json.gz', 'rt'))['cases']

def initialize():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'logs').mkdir(exist_ok=True)
    cases=[c for c in load_records() if c['model']=='wave_pinning' and c['condition']=='relative_20pct' and c['length'] in [5.,13.,39.]]
    jobs=[]
    for case in cases:
        for gate in r.DEFINITIONS:
            selected=min((t for t in case['tests'] if t['duration']>0 and t[gate+'_success']), key=lambda t:t['duration'])
            failed=max((t for t in case['tests'] if t['duration']<selected['duration'] and not t[gate+'_success']),key=lambda t:t['duration'])
            for category,t in [('pass',selected),('fail',failed)]:
                if any(j['case']['case_index']==case['case_index'] and j['duration']==t['duration'] for j in jobs):
                    next(j for j in jobs if j['case']['case_index']==case['case_index'] and j['duration']==t['duration'])['gates'].append(gate+'_'+category)
                    continue
                jobs.append(dict(slot=len(jobs),case={k:v for k,v in case.items() if k not in ['tests','result']},
                                 duration=t['duration'],gates=[gate+'_'+category],original_row=t,extra_recovery=40000. if case['length']==39. else 0.))
    n.write(OUT/'jobs.json',jobs)
    sources={str(p.relative_to(MANUSCRIPT)):n.sha(p) for p in [Path(__file__),MANUSCRIPT/'code/gradient_reversal.py',MANUSCRIPT/'code/core/native.py',MANUSCRIPT/'code/core/models/wave_pinning.py',MANUSCRIPT/'code/core/src/polarity1d/solvers.py',CAMPAIGN/'manifest.json',CAMPAIGN/'results/records.json.gz']}
    n.write(OUT/'provenance.json',dict(generated_at=n.now(),sources=sources,jobs_sha256=n.sha(OUT/'jobs.json')))
    print(json.dumps([dict(slot=j['slot'],length=j['case']['length'],duration=j['duration'],gates=j['gates']) for j in jobs],indent=2))

def gate_metrics(tr,case):
    tail=tr.time>=r.RECOVERY-r.TAIL
    o=tr.orientation[tail]; contrast=tr.contrast[tail]
    orientation=bool(np.min(o)>.05 and np.min(contrast)>=.05 and np.ptp(o)<.005)
    error=float(np.max(abs(o/abs(case['initial_orientation'])-1)))
    fields=r.settling.monitored(tr,case['model'],case['length'])
    fieldmax=0.
    for name,y in fields.items():
        a=y[tail];scale=case['fixed_scales'][name]
        fieldmax=max(fieldmax,float(np.max(np.ptp(a,axis=0)))/scale,float(np.max(abs(a[-1]-a[-2])))/r.STEP*r.TAIL/scale)
    return dict(orientation_success=orientation,full_polarity_success=bool(orientation and error<=.01),final_orientation=float(o[-1]),minimum_contrast_last_500=float(contrast.min()),orientation_range_last_500=float(np.ptp(o)),maximum_relative_polarity_error_last_500=error,field_tail_pass=bool(fieldmax<=.001),field_tail_max_scaled=fieldmax)

def worker(slot):
    job=json.loads((OUT/'jobs.json').read_text())[slot];c=job['case'];D=job['duration'];L=c['length']
    checkpoint=n.ROOT/c['checkpoint']
    if n.sha(checkpoint)!=c['checkpoint_sha256']:raise RuntimeError('checkpoint hash mismatch')
    state,prep=n.read_checkpoint(checkpoint,'wave_pinning')
    # Output sampling changes, while native BDF tolerances and max step are preserved.
    pulse=n.run('wave_pinning',L,state,g.response_cue(c),D,interval=.1,n=c['cells'])
    recovery=n.run('wave_pinning',L,pulse.state,e.withdrawal(c),r.RECOVERY,interval=r.STEP,n=c['cells'])
    metrics=gate_metrics(recovery,c)
    original=job['original_row']
    for key in ['orientation_success','full_polarity_success','field_tail_pass']:
        if metrics[key]!=original[key]:raise RuntimeError(f'reconstruction {key} differs')
    if abs(metrics['final_orientation']-original['final_orientation'])>1e-9:raise RuntimeError('reconstruction scalar differs')
    times=np.r_[pulse.time,D+recovery.time[1:]]
    u=np.r_[pulse.fields['active'],recovery.fields['active'][1:]]
    v=np.r_[pulse.fields['inactive'],recovery.fields['inactive'][1:]]
    o=np.r_[pulse.orientation,recovery.orientation[1:]]
    if job['extra_recovery']:
        extra=n.run('wave_pinning',L,recovery.state,e.withdrawal(c),job['extra_recovery'],interval=10.,n=c['cells'])
        times=np.r_[times,D+r.RECOVERY+extra.time[1:]]
        u=np.r_[u,extra.fields['active'][1:]];v=np.r_[v,extra.fields['inactive'][1:]];o=np.r_[o,extra.orientation[1:]]
        metrics['orientation_after_50000_recovery']=float(extra.orientation[-1])
        metrics['contrast_after_50000_recovery']=float(extra.contrast[-1])
    x=pulse.x;dx=L/c['cells'];m=n.model_for('wave_pinning')
    cue=np.where(times[:,None]<=D,c['gradient_slope']*x[None,:],0.)
    basal=v*(m.k0+m.feedback*u*u/(m.half_saturation**2+u*u))-m.decay*u
    stimulus_transfer=v*cue
    net=basal+stimulus_transfer
    flux_u=-m.du*np.diff(u,axis=1)/dx;flux_v=-m.dv*np.diff(v,axis=1)/dx
    diffusion_u=m.du*laplacian_no_flux(u,dx);diffusion_v=m.dv*laplacian_no_flux(v,dx)
    mass=(u+v).sum(axis=1)*dx
    metrics['mass_relative_error']=float(np.max(abs(mass/mass[0]-1)))
    metrics['checkpoint_sha256']=c['checkpoint_sha256']
    archive=OUT/f'trajectory_{slot:02d}.npz'
    np.savez_compressed(archive,time=times,x=x,active=u,inactive=v,orientation=o,contrast=np.ptp(u,axis=1),
                        pulse_duration=D,length=L,slope=c['gradient_slope'],cue=cue,basal_transfer=basal,
                        stimulus_transfer=stimulus_transfer,net_transfer=net,active_flux=flux_u,inactive_flux=flux_v,
                        active_diffusion=diffusion_u,inactive_diffusion=diffusion_v,total_mass=mass,
                        mean_inactive=v.mean(axis=1),mean_active=u.mean(axis=1),prep_active=prep.fields['active'][-1],prep_inactive=prep.fields['inactive'][-1])
    n.write(archive.with_suffix('.json'),dict(job=job,metrics=metrics,generated_at=n.now(),archive_sha256=n.sha(archive)))
    print(json.dumps(dict(slot=slot,length=L,duration=D,metrics=metrics),indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['init','worker']);parser.add_argument('--slot',type=int)
    args=parser.parse_args()
    if args.action=='init':initialize()
    else:worker(args.slot)
