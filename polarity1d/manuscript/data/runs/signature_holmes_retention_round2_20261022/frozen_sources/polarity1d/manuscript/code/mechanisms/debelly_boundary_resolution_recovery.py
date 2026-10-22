"""Withdraw each audited variant at its own recorded MCA zero crossing."""
from pathlib import Path
import sys,numpy as np
from debelly_boundary_resolution import OUT,n,r,run,save,polarity

def setup():
 jobs=r.read(OUT/'jobs.json');assert all((OUT/f'{j["variant"]}_L{j["length"]:g}_summary.json').exists() for j in jobs)
 r.write(OUT/'recovery_provenance.json',dict(generated_at=n.now(),script_sha256=n.sha(Path(__file__)),original_sources=r.read(OUT/'provenance.json')['sources'],cue_sources={str(OUT/f'{j["variant"]}_L{j["length"]:g}_cue.npz'):n.sha(OUT/f'{j["variant"]}_L{j["length"]:g}_cue.npz') for j in jobs}))
def worker(slot):
 prov=r.read(OUT/'recovery_provenance.json');assert n.sha(Path(__file__))==prov['script_sha256']
 for path,d in prov['original_sources'].items():assert n.sha(Path(path))==d
 job=r.read(OUT/'jobs.json')[slot];v=job['variant'];L=job['length'];summary=r.read(OUT/f'{v}_L{L:g}_summary.json');path=OUT/f'{v}_L{L:g}_cue.npz';assert n.sha(path)==prov['cue_sources'][str(path)]
 if summary['crossing_time'] is None:r.write(OUT/f'{v}_L{L:g}_recovery.json',dict(available=False,reason='No MCA zero crossing'));return
 with np.load(path) as z:
  i=np.searchsorted(z['time'],summary['crossing_time']);state={k:z['field_'+k][i].copy() for k in n.DB_NAMES}
 tr=run(v,L,state,n.ConstantInput(0),10000.,interval=5.);save(OUT/f'{v}_L{L:g}_recovery.npz',tr);P=polarity(tr.fields['rac'],tr.x,L);M=polarity(tr.fields['mca_bound'],tr.x,L);tail=tr.time>=9500
 r.write(OUT/f'{v}_L{L:g}_recovery.json',dict(generated_at=n.now(),available=True,length=L,variant=v,pulse_duration=summary['crossing_time'],final_rac_polarity=float(P[-1]),final_mca_polarity=float(M[-1]),retained_right_orientation=bool(P[tail].min()>.05 and np.ptp(P[tail])<.005 and np.ptp(tr.fields['rac'][tail],axis=1).min()>=.05),rac_tail_min=float(P[tail].min()),rac_tail_range=float(np.ptp(P[tail])),field_tail_range={k:float(np.max(np.ptp(tr.fields[k][tail],axis=0))) for k in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension']},sha256=n.sha(OUT/f'{v}_L{L:g}_recovery.npz')))
if __name__=='__main__':
 if sys.argv[1]=='setup':setup()
 else:worker(int(sys.argv[2]))
