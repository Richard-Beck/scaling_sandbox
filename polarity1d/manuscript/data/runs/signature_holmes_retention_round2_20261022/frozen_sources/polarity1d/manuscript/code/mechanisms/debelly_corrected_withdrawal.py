"""Corrected withdrawal diagnostics, including explicit ineligible short-cell branch."""
from pathlib import Path
import sys,numpy as np
import debelly_corrected_zero_crossing as z
n,r,OUT=z.n,z.r,z.OUT

def worker(slot):
 L=[18.,23.,39.][slot-22]
 cases=r.read(z.corrected.DEFAULT/'manifest.json')['cases']
 case=next(c for c in cases if c['model']=='debelly' and c['condition']=='fixed_slope' and c['length']==L)
 if case['eligible']:return z.worker(slot)
 assert L==18.,('Unexpected missing polarized long-cell preparation',L)
 records=[]
 for delta,D in zip([-100.,-50.,0.],[100.,500.,2000.]):
  cue=z.corrected.e.Patch(case['amplitude'],.2*L,side='left',width_mode='relative',fraction=.2)
  full=n.run('debelly',L,None,cue,D,interval=.5)
  rec=n.run('debelly',L,full.state,n.ConstantInput(0),10000.,interval=2.5)
  arrays={'x':full.x,'duration':np.array(D),'uninterrupted_crossing':np.array(np.nan),'preparation_ineligible':np.array(True)}
  for phase,tr in [('pulse',full),('recovery',rec)]:
   arrays[phase+'_time']=tr.time
   arrays.update({phase+'_field_'+k:v for k,v in tr.fields.items()})
  path=OUT/f'withdrawal_L{L:g}_offset{delta:g}.npz';np.savez_compressed(path,**arrays)
  records.append(dict(length=L,offset=delta,duration=D,kind='left_establishment_withdrawal_diagnostic',qualified_reversal=False,file_sha256=n.sha(path)))
 r.write(OUT/f'withdrawal_L{L:g}.json',records)
 r.write(OUT/f'done_{slot}.json',dict(generated_at=n.now(),kind='ineligible_short_cell_establishment_diagnostic',qualified_reversal=False))

if __name__=='__main__':worker(int(sys.argv[1]))
