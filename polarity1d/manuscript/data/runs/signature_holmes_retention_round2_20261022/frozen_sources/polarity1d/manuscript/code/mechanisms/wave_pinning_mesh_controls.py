"""Bounded native-parameter401-cell check of long-cell interior competition."""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import wave_pinning_mechanism as w
n,g,e=w.n,w.g,w.e

def run(slot):
    source,duration=[(4,39.),(7,49.25)][slot]
    c=dict(json.loads((w.OUT/'jobs.json').read_text())[source]['case'])
    initial,prep=n.read_checkpoint(n.ROOT/c['checkpoint'],'wave_pinning')
    oldmean=float(initial.sum(axis=0).mean());newx=n.Grid1D(c['length'],401).x
    state=np.array([np.interp(newx,prep.x,a) for a in initial]);correction=oldmean-float(state.sum(axis=0).mean());state[1]+=correction
    initialori=float((state[0]@(2*newx/(newx[0]+newx[-1])-1))/state[0].sum())
    c['cells']=401;c['initial_orientation']=initialori
    pulse=n.run('wave_pinning',c['length'],state,g.response_cue(c),duration,interval=.1,n=401)
    recovery=n.run('wave_pinning',c['length'],pulse.state,e.withdrawal(c),10000.,interval=2.5,n=401)
    metrics=w.gate_metrics(recovery,c)
    extra=n.run('wave_pinning',c['length'],recovery.state,e.withdrawal(c),40000.,interval=10.,n=401)
    times=np.r_[pulse.time,duration+recovery.time[1:],duration+10000+extra.time[1:]]
    active=np.r_[pulse.fields['active'],recovery.fields['active'][1:],extra.fields['active'][1:]]
    inactive=np.r_[pulse.fields['inactive'],recovery.fields['inactive'][1:],extra.fields['inactive'][1:]]
    orientation=np.r_[pulse.orientation,recovery.orientation[1:],extra.orientation[1:]]
    path=w.OUT/f'mesh_control_{slot:02d}.npz'
    np.savez_compressed(path,time=times,x=newx,active=active,inactive=inactive,orientation=orientation,length=c['length'],duration=duration,cells=401)
    n.write(path.with_suffix('.json'),dict(generated_at=n.now(),source_slot=source,case=c,duration=duration,metrics=metrics,orientation_after50000=float(extra.orientation[-1]),inactive_mass_correction=correction,initial_total_mean=oldmean,
        archive_sha256=n.sha(path),code_sha256=n.sha(Path(__file__)),note='Diagnostic401-cell mesh control. Native preparation interpolated, with constant inactive correction to retain exacttotalmass. Native physics and same slope/pulse retained. This is not an exact campaign checkpoint or a threshold refinement.'))

if __name__=='__main__':run(int(sys.argv[1]))
