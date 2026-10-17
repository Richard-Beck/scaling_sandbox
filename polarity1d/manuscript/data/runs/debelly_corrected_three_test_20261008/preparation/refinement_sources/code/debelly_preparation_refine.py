"""Continue unresolved preparation-duration brackets without changing dynamics.

The inherited eight threshold steps may leave a zero-to-positive bracket open
when an unstable homogeneous state responds to extremely small cues. Continue
bounded local refinement; do not call such a bracket a reversal failure.
"""
from pathlib import Path
import sys,copy,shutil,json
import debelly_corrected_reversal as d
n,r,e=d.n,d.r,d.e
OUT=d.DEFAULT

def setup():
 manifest=d.prep_verify(OUT);manifest=copy.deepcopy(manifest)
 manifest['refinement']['threshold_steps']=56
 r.write(OUT/'preparation/refinement_plan.json',manifest)
 sources=d.sources();sources[str(Path(__file__).relative_to(n.ROOT))]=n.sha(Path(__file__))
 for rel in sources:
  p=OUT/'preparation/refinement_sources'/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(n.ROOT/rel,p)
 r.write(OUT/'preparation/refinement_provenance.json',dict(generated_at=n.now(),sources=sources,plan_sha256=n.sha(OUT/'preparation/refinement_plan.json'),reason='Continue inherited unresolved preparation lower brackets; all scientific dynamics and verified-release gates unchanged. Preparation pulse duration is not the reversal readout.'))

def worker(slot):
 v=r.read(OUT/'preparation/refinement_provenance.json')
 for rel,digest in v['sources'].items():assert n.sha(n.ROOT/rel)==digest==n.sha(OUT/'preparation/refinement_sources'/rel)
 assert n.sha(OUT/'preparation/refinement_plan.json')==v['plan_sha256']
 manifest=r.read(OUT/'preparation/refinement_plan.json');p=OUT/'preparation/work/cases'/f'case_{slot:04d}/result.json'
 original=r.read(p)
 if not original['summary']['duration_boundary_unresolved']:return
 shutil.copy2(p,p.with_name('initial_eight_step_result.json'))
 e.verify=lambda _:manifest
 e.worker(OUT/'preparation',slot)
 result=r.read(p)
 if result['summary']['duration_boundary_unresolved']:raise RuntimeError('Preparation lower boundary remains unresolved after bounded extension')
 r.write(p.with_name('refinement_audit.json'),dict(generated_at=n.now(),passed=True,original_result_sha256=n.sha(p.with_name('initial_eight_step_result.json')),refined_result_sha256=n.sha(p),provenance_sha256=n.sha(OUT/'preparation/refinement_provenance.json')))
if __name__=='__main__':
 if sys.argv[1]=='setup':setup()
 else:worker(int(sys.argv[2]))
