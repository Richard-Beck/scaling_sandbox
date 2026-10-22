"""Bounded ramp-strength controls with native physics and exact preparations."""
from pathlib import Path
import json,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import wave_pinning_mechanism as w
n,g,e=w.n,w.g,w.e

def run(slot):
    source,factor,duration=[(0,2.6,13.875),(4,1./3.,39.),(4,.5,39.),(4,2.,39.),(0,1.,13.875)][slot]
    jobs=json.loads((w.OUT/'jobs.json').read_text());c=dict(jobs[source]['case']);checkpoint=n.ROOT/c['checkpoint']
    if n.sha(checkpoint)!=c['checkpoint_sha256']:raise RuntimeError('checkpoint changed')
    initial,prep=n.read_checkpoint(checkpoint,'wave_pinning');c['gradient_slope']*=factor
    pulse=n.run('wave_pinning',c['length'],initial,g.response_cue(c),duration,interval=.1,n=c['cells'])
    recovery=n.run('wave_pinning',c['length'],pulse.state,e.withdrawal(c),10000.,interval=2.5,n=c['cells'])
    metrics=w.gate_metrics(recovery,c)
    p=w.OUT/f'cue_control_{slot:02d}.npz'
    np.savez_compressed(p,time=np.r_[pulse.time,duration+recovery.time[1:]],x=pulse.x,
                        active=np.r_[pulse.fields['active'],recovery.fields['active'][1:]],
                        inactive=np.r_[pulse.fields['inactive'],recovery.fields['inactive'][1:]],
                        orientation=np.r_[pulse.orientation,recovery.orientation[1:]],length=c['length'],duration=duration,slope=c['gradient_slope'])
    n.write(p.with_suffix('.json'),dict(generated_at=n.now(),source_slot=source,cue_slope_factor=factor,case=c,duration=duration,metrics=metrics,archive_sha256=n.sha(p),code_sha256=n.sha(Path(__file__)),
        note='Diagnostic single-duration ramp-strength intervention. Native reaction/diffusion, exact native preparation and standard10000 recovery unchanged. Not a new threshold search or revised campaign.'))

if __name__=='__main__':run(int(sys.argv[1]))
