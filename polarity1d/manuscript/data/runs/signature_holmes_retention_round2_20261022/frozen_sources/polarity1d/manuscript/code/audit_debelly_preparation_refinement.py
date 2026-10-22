"""Read-only provenance and native scalar selection audit for extended preparations."""
from pathlib import Path
import json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/runs/debelly_corrected_three_test_20261008/preparation'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=read(OUT/'refinement_provenance.json');assert sha(OUT/'refinement_plan.json')==p['plan_sha256']
 for rel,digest in p['sources'].items():assert sha(ROOT/rel)==digest==sha(OUT/'refinement_sources'/rel)
 rows=[]
 for f in sorted((OUT/'work/cases').glob('case_*/result.json')):
  d=read(f);s=d['summary'];old=f.with_name('initial_eight_step_result.json')
  if old.exists():
   audited=read(f.with_name('refinement_audit.json'));assert sha(old)==audited['original_result_sha256'] and sha(f)==audited['refined_result_sha256']
  assert not s['duration_boundary_unresolved']
  if s['status']=='polarized':
   valid=[t for t in d['trials'] if t['verified'] and t['final_contrast']>=.05]
   best=max(t['cue_aligned_strength'] for t in valid)
   passing=[t for t in valid if t['cue_aligned_strength']>=.99*best and t['polarity_arrival'] is not None]
   selected=min(passing,key=lambda t:t['duration']);assert selected['duration']==s['cue_duration']
   lower=max(t['duration'] for t in d['trials'] if t['duration']<selected['duration']) if selected['duration'] else 0
   assert (selected['duration']-lower)/selected['duration']<=.010000001 if selected['duration'] else True
   for t in passing:
    with np.load(OUT/t['trace']) as z:
     assert abs(z['release_orientation'][-1]-t['final_orientation'])<1e-12
   assert sha(OUT/selected['checkpoint'])==s['selected_checkpoint_sha256']
  rows.append(dict(length=s['length'],status=s['status'],selected_duration=s['cue_duration'],lower_bracket_resolved=not s['duration_boundary_unresolved'],extended=old.exists(),trials=s['tested_durations']))
 assert len(rows)==22
 (OUT/'refinement_independent_audit.json').write_text(json.dumps(dict(passed=True,source_provenance_verified=True,all_22_lengths=True,rows=rows),indent=2)+'\n')
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
