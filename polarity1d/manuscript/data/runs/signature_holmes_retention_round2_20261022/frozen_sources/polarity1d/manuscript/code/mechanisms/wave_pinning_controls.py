"""Post-cutoff mobility controls from identical native fixed-slope pulse states."""
from pathlib import Path
from dataclasses import replace
import json,sys
import numpy as np
MANUSCRIPT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(MANUSCRIPT/'code/core'))
import native as n
OUT=MANUSCRIPT/'data/runs/fixed_slope_mechanisms_20261006/wave_pinning'

def run(slot):
    kind,source=[('inactive_fast',4),('active_fast',4),('inactive_fast',7),('active_fast',7),('inactive_10x',4),('inactive_10x',7),('inactive_10x_prep',4),('inactive_fast_prep',4)][slot]
    z=np.load(OUT/f'trajectory_{source:02d}.npz');D=float(z['pulse_duration']);L=float(z['length'])
    i=int(np.argmin(abs(z['time']-D)));state=np.vstack([z['prep_active'],z['prep_inactive']]) if kind.endswith('_prep') else np.vstack([z['active'][i],z['inactive'][i]])
    m=n.model_for('wave_pinning');m=replace(m,dv=1000.) if kind.startswith('inactive_fast') else replace(m,dv=100.) if kind.startswith('inactive_10x') else replace(m,du=1.)
    tr=n._sim(m,state,L,201,n.ConstantInput(0.),10000.,interval=2.5)
    ori=n._orientation(tr);u=tr.fields['active'];v=tr.fields['inactive'];x=tr.x
    # Counterfactual changes only one mobility after stimulus withdrawal.
    path=OUT/f'control_{slot:02d}.npz'
    np.savez_compressed(path,recovery_time=tr.time,active=u,inactive=v,orientation=ori,x=x,
                        native_pulse_duration=D,length=L,du=m.du,dv=m.dv,source_slot=source)
    n.write(path.with_suffix('.json'),dict(generated_at=n.now(),kind=kind,source_slot=source,source_sha256=n.sha(OUT/f'trajectory_{source:02d}.npz'),archive_sha256=n.sha(path),code_sha256=n.sha(Path(__file__)),du=m.du,dv=m.dv,
        source_state='native preparation' if kind.endswith('_prep') else 'native pulse cutoff',
        note='Exact native state, no cue. Only declared diffusion coefficient changes during diagnostic continuation; this is not a revised native result.'))

if __name__=='__main__':run(int(sys.argv[1]))
