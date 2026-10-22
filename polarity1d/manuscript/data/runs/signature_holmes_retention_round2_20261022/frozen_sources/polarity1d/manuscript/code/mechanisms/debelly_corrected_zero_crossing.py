"""Dedicated native deBelly MCA-withdrawal and basal perturbation experiment."""
from pathlib import Path
import sys,json,numpy as np
from debelly_controls import n,r,base,g
import debelly_corrected_reversal as corrected

OUT=n.DATA/'runs/debelly_corrected_zero_crossing_20261008'
LENGTHS=[12.,16.,18.,19.,19.5,20.,21.,23.,27.,33.,39.]
KEYS=['rac','rho','mca_bound','mca_unbound','displacement','velocity']
def pol(u,x,L):return u@(2*x/L-1)/u.sum(axis=1)
def save(path,tr):
 np.savez_compressed(path,time=tr.time,x=tr.x,orientation=tr.orientation,**{'field_'+k:v for k,v in tr.fields.items()})
def tangent(L):
 x=n.rd.Grid1D(L,201).x;c=np.cos(np.pi*x/L);s=np.sin(np.pi*x/L);b=base(L);J=np.zeros((6,6))
 for j,k in enumerate(KEYS):
  eps=1e-6 if k in ['rac','rho','displacement','velocity'] else 1e-7;ends=[]
  for sign in [-1,1]:
   state={k:v.copy() for k,v in b.items()};state[k]+=sign*eps*(s if k in ['displacement','velocity'] else c)
   ends.append(n.run('debelly',L,state,n.ConstantInput(0),.5,interval=.5).state)
  for i,target in enumerate(KEYS):
   q=s if target in ['displacement','velocity'] else c
   J[i,j]=((ends[1][target]-ends[0][target])/(2*eps))@q/(q@q)
 return float(np.log(max(abs(np.linalg.eigvals(J))))/.5)
def setup():
 OUT.mkdir(parents=True,exist_ok=True)
 jobs=[dict(kind='seed',length=L,epsilon=eps) for L in LENGTHS for eps in [1e-5,1e-6]]
 jobs += [dict(kind='withdrawal',length=L) for L in [18.,23.,39.]]
 jobs += [dict(kind='linear')]
 r.write(OUT/'jobs.json',jobs)
 r.write(OUT/'provenance.json',dict(generated_at=n.now(),script_sha256=n.sha(Path(__file__)),model_sha256=n.sha(Path(n.ds.__file__)),adapter_sha256=n.sha(Path(corrected.__file__)),center_boundary_motion=False,definition='Homogeneous algebraic basal rest; zero-mean cosine bound-MCA seeds 1e-5 and 1e-6; no external cue; 40000-unit horizon. Withdrawal uses exact current relative20 basal checkpoints and corrected independent-endpoint slope0.24 gradient; fresh corrected relative20 patch/basal preparations.'))
 print(len(jobs))
def worker(slot):
 assert n.sha(Path(__file__))==r.read(OUT/'provenance.json')['script_sha256']
 assert n.sha(Path(corrected.__file__))==r.read(OUT/'provenance.json')['adapter_sha256']
 job=r.read(OUT/'jobs.json')[slot]
 if job['kind']=='linear':
  grid=np.r_[np.arange(12,39.01,.5),19.55,19.6]
  r.write(OUT/'linear.json',[dict(length=float(L),rate=tangent(float(L))) for L in sorted(set(grid))]);return
 L=job['length']
 if job['kind']=='seed':
  state=base(L);state['mca_bound']+=job['epsilon']*np.cos(np.pi*n.rd.Grid1D(L,201).x/L)
  # Equal total MCA, positive fields; relaxed homogeneous start before perturbation.
  tr=n.run('debelly',L,state,n.ConstantInput(0),40000.,interval=20)
  save(OUT/f'seed_L{L:g}_e{job["epsilon"]:g}.npz',tr)
 else:
  cases=r.read(corrected.DEFAULT/'manifest.json')['cases']
  case=next(c for c in cases if c['model']=='debelly' and c['condition']=='fixed_slope' and c['length']==L)
  assert case['eligible'], ('No qualified corrected basal state',L)
  cp=n.ROOT/case['checkpoint'];assert n.sha(cp)==case['checkpoint_sha256']
  state,_=n.read_checkpoint(cp,'debelly')
  cue=corrected.response_cue(case);full=n.run('debelly',L,state,cue,1200.,interval=.5)
  pp=pol(full.fields['mca_bound'],full.x,L);ix=np.flatnonzero(pp<=0)[0];cross=float(full.time[ix])
  records=[]
  for delta in [-100.,-50.,0.]:
   D=max(.5,cross+delta);i=int(np.searchsorted(full.time,D));end={k:full.fields[k][i].copy() for k in n.DB_NAMES}
   rec=n.run('debelly',L,end,n.ConstantInput(0),10000.,interval=2.5)
   arrays={'x':full.x,'duration':np.array(D),'uninterrupted_crossing':np.array(cross)}
   for phase,tt,fields in [('pulse',full.time[:i+1],{k:v[:i+1] for k,v in full.fields.items()}),('recovery',rec.time,rec.fields)]:
    arrays[phase+'_time']=tt
    arrays.update({phase+'_field_'+k:v for k,v in fields.items()})
   path=OUT/f'withdrawal_L{L:g}_offset{delta:g}.npz';np.savez_compressed(path,**arrays)
   records.append(dict(length=L,offset=delta,duration=D,uninterrupted_crossing=cross,checkpoint_sha256=n.sha(cp),file_sha256=n.sha(path)))
  r.write(OUT/f'withdrawal_L{L:g}.json',records)
 (OUT/f'done_{slot}.json').write_text(json.dumps(dict(job=job,finished_at=n.now())))
if __name__=='__main__':
 if sys.argv[1]=='setup':setup()
 else:worker(int(sys.argv[2]))
