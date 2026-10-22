"""A labeled mesh check; baseline evidence still uses its exact 201-cell states.

The finer-grid diagnostic remaps original fields and corrects the spatial mean
of each inactive pool to preserve its original conserved-group mass. No-cue
controls measure drift introduced by remapping; this does not replace any prep.
"""
import argparse,json
import numpy as np
import holmes_mechanism as h
from polarity1d.contract import Grid1D
from polarity1d.simulation import simulate

def tasks():
    out=[]
    for L in [10,28,50]:
        for gate in ['zero','orientation','full_polarity']:
            t=next(t for t in h.trials() if t['length']==L and t['gate']==('orientation' if gate=='zero' or L==10 else gate) and t['result']=='pass')
            out.append(dict(length=L,gate=gate,duration=0. if gate=='zero' else t['duration'],case=t['case']))
    return out

def worker(slot):
    task=tasks()[slot];case=task['case'];L=task['length'];T=task['duration'];N=401
    state,prep=h.n.read_checkpoint(h.ROOT/case['checkpoint'],'holmes_model3');newx=Grid1D(L,N).x
    remap=np.vstack([np.interp(newx,prep.x,y) for y in state])
    for a,b in [(0,1),(2,3)]:remap[b]+=float((state[a]+state[b]).mean()-(remap[a]+remap[b]).mean())
    model=h.n.model_for('holmes_model3');kwargs=dict(length=L,cells=N,backend='bdf',max_step=10.,rtol=2e-6,atol=1e-9)
    final=remap
    if T:
        tr=simulate(model,initial_state=remap,protocol=h.g.response_cue(case),end_time=T,output_interval=T/40,**kwargs)
        final=np.vstack([tr.fields[k][-1] for k in model.state_names])
    tr=simulate(model,initial_state=final,protocol=h.n.ConstantInput(0),end_time=10000.,output_interval=2.5,**kwargs)
    s,o=h.history(tr);tail=tr.time>=9500;contrast=np.ptp(s[:,0],axis=1)
    og=bool(np.min(o[tail])>.05 and np.ptp(o[tail])<.005 and np.min(contrast[tail])>=.05)
    fg=bool(og and np.max(abs(o[tail]/abs(case['initial_orientation'])-1))<=.01)
    path=h.OUT/'mesh_check';path.mkdir(parents=True,exist_ok=True);tag=f'L{L}_{task["gate"]}_N{N}'
    np.savez_compressed(path/(tag+'.npz'),x=tr.x,initial=remap,final_state=s[-1],time=tr.time,orientation=o)
    h.write(path/(tag+'.json'),dict(length=L,gate=task['gate'],duration=T,cells=N,baseline_cells=201,
           orientation_success=og,full_polarity_success=fg,initial_original_orientation=case['initial_orientation'],final_orientation=float(o[-1]),
           zero_scaled_excursion=float(max(np.max(abs(s[:,i]-remap[i]))/case['fixed_scales'][k] for i,k in enumerate(model.state_names))) if T==0 else None,
           checkpoint_sha256=case['checkpoint_sha256'],script_sha256=h.n.sha(__file__),generated_at=h.n.now(),
           diagnostic='Conservative-group-preserving field interpolation to401cells; original201-cellstatesremaincanonical.'))
    print(tag,T,og,fg,float(o[-1]),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--slot',type=int,required=True);worker(p.parse_args().slot)
