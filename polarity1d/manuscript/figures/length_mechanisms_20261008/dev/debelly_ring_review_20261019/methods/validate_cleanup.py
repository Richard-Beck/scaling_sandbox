"""Verify the reviewed retirement without deleting or regenerating anything."""
from pathlib import Path
from datetime import datetime,timezone
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import hashlib,json,os,subprocess

PACKAGE=Path(__file__).resolve().parents[1]
ROOT=next(p for p in PACKAGE.parents if (p/'polarity1d/manuscript/AGENTS.md').exists())
PLAN_SHA='73d408ba1109b9ebdcbc6c5bfa6f880a94379f1bb7a66d8f56db1c91d1c82bfb'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def resolve(label):
 p=Path(label);return p if p.is_absolute() else ROOT/p
def main():
 hashes={}
 def pin(p):hashes[str(p.relative_to(ROOT))]=sha(p)
 def verify(p,digest):assert p.is_file() and sha(p)==digest,str(p)
 plan=read(PACKAGE/'cleanup_plan.json');assert sha(PACKAGE/'cleanup_plan.json')==PLAN_SHA
 receipt=read(PACKAGE/'cleanup_receipt.json');assert receipt['passed'] and receipt['plan_sha256']==PLAN_SHA
 removed=[]
 for row in plan['candidate_delete']:
  p=resolve(row['path']);assert not p.exists() and not p.is_symlink(),str(p);removed.append(row['path'])
 assert {r['path'] for r in receipt['removed']}==set(removed)
 protected={}
 for row in plan['protected']:
  p=resolve(row['path']);assert p.exists(),str(p)
  if 'sha256' in row:protected[row['path']]=row['sha256']
  protected.update(row.get('protected_hash_samples',{}))
 for name,digest in protected.items():verify(resolve(name),digest)
 packaging=read(PACKAGE/'methods/viewer_packaging_manifest.json')
 for name,row in packaging['files'].items():verify(PACKAGE/name,row['sha256'])
 evidence=read(PACKAGE/'evidence_manifest.json')
 for row in evidence['assets']:verify(PACKAGE/row['retained_path'],row['retained_sha256'])
 prior=read(PACKAGE/'package_validation.json');assert prior['passed']
 for name,digest in prior['source_sha256'].items():
  if name=='index.html':continue # Root adds final retirement/storage text only.
  verify(PACKAGE/name,digest)
 assert prior['checks']['all18_topology_maps_preserve_all7_piecewise_constant_profiles_exactly']
 assert prior['checks']['all9_late_windows_continue_prior_phase_direction_for_every60_native_halfsteps']
 browser=read(PACKAGE/'methods/post_cleanup_viewer_result.json')
 assert browser['passed'] and len(browser['checks'])==21 and not browser['runtime_errors']
 assert all(r['all_fields']==7 for r in browser['checks'])
 class Links(HTMLParser):
  def __init__(self):super().__init__();self.links=[]
  def handle_starttag(self,tag,attrs):
   for k,v in attrs:
    if k in ('href','src') and v:self.links.append(v)
 links=[]
 for p in PACKAGE.rglob('*.html'):
  parser=Links();text=p.read_text();parser.feed(text)
  assert 'http-equiv="refresh"' not in text.lower() and 'location.reload' not in text
  for url in parser.links:
   u=urlsplit(url)
   if u.scheme or u.netloc or not u.path:continue
   dest=(p.parent/unquote(u.path)).resolve();assert dest.is_relative_to(PACKAGE.resolve()),str(dest)
   if dest!=PACKAGE/'cleanup_validation.json':assert dest.exists(),str(dest)
   links.append(dict(source=str(p.relative_to(PACKAGE)),target=str(dest.relative_to(PACKAGE))))
 # Shared entry points must route to the retained review, not the retired topic.
 shared=[ROOT/'DATA_MANIFEST.md',ROOT/'polarity1d/manuscript/.gitignore',ROOT/'polarity1d/manuscript/code/mechanisms/organize_length_dev.py',PACKAGE.parent/'README.md']
 obsolete='debelly_ring_development_20261019'
 for p in shared:
  assert p.is_file(),str(p)
  text=p.read_text();assert obsolete not in text,str(p)
  if p.name!='.gitignore':assert 'debelly_ring_review_20261019' in text,str(p)
  pin(p)
 assert not (ROOT/'polarity1d/manuscript/DEBELLY_RING_STATUS_20261019.md').exists()
 dated_runs=list((ROOT/'polarity1d/manuscript/data/runs').glob('*debelly*20261019*'))
 dated_sources=[p for p in (ROOT/'polarity1d/manuscript/code/mechanisms').iterdir() if p.is_file() and 'debelly' in p.name.lower() and '20261019' in p.name]
 assert not dated_runs and not dated_sources
 queue=subprocess.run(['squeue','-u',os.environ.get('USER','4473331'),'-h','-o','%i|%j|%T|%R'],text=True,capture_output=True)
 assert queue.returncode==0,queue.stderr
 jobs=queue.stdout.splitlines();assert not any('debelly' in j.lower() for j in jobs)
 producers=[]
 for line in subprocess.check_output(['ps','-eo','uid,pid,args','--no-headers'],text=True).splitlines():
  x=line.split(None,2)
  if len(x)!=3 or x[0]!=str(os.getuid()) or int(x[1])==os.getpid():continue
  command=x[2]
  if 'debelly' not in command.lower() or str(PACKAGE.relative_to(ROOT)) in command:continue
  if any(row['path'] in command for row in plan['candidate_delete']):producers.append(dict(pid=int(x[1]),command=command))
 assert not producers,producers
 for name in ['cleanup_plan.json','cleanup_receipt.json','package_validation.json','executive_assessment.json','evidence_manifest.json','index.html','supporting/metadata/storage_footprint.json','methods/viewer_packaging_manifest.json','methods/post_cleanup_viewer_result.json','methods/report_browser_checks_result.json']:
  pin(PACKAGE/name)
 pin(Path(__file__))
 result=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,mode='INDEPENDENT_POST_RETIREMENT_VERIFICATION',scope='Consolidation/cleanup complete. The scientific attempt remains unsuccessful as a validated periodic deBelly replacement.',checks=dict(all304_reviewed_candidate_paths_absent=True,all177_protected_paths_present=True,all391_original_protected_hash_samples_unchanged=True,all70_frozen_viewer_method_files_unchanged=True,all23_curated_evidence_assets_unchanged=True,retained_package_validation_and_event_claims_still_hash_valid=True,all_local_HTML_links_resolve_within_package=True,shared_entry_points_and_storage_manifest_updated=True,no_dated_debelly_runs_or_standalone_sources_remaining=True,no_named_debelly_SLURM_jobs_or_candidate_producers=True),removed_paths=len(removed),workspace_allocated_bytes_retired=receipt['old_workspace_allocated_bytes'],system_temporary_allocated_bytes_retired=receipt['old_system_temporary_allocated_bytes'],protected_paths=len(plan['protected']),protected_hash_samples=len(protected),retained_viewer_method_files=packaging['file_count'],curated_evidence_assets=len(evidence['assets']),local_links=links,live_state=dict(squeue_jobs=jobs,matching_candidate_producers=producers),source_sha256=hashes,limits=['Full retired scientific histories are intentionally unavailable. Historical audit paths remain provenance identifiers.','The retained event viewer and its exporter use local raw inputs; package verification does not establish physical closure, convergence or a validated ring assay.','Independent Spring/revision work and established corrected linear evidence are protected and unchanged.'])
 target=PACKAGE/'cleanup_validation.json';target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 assert all((PACKAGE/r['target']).exists() for r in links)
 print(json.dumps(dict(passed=True,path=str(target.relative_to(ROOT)),sha256=sha(target),removed=len(removed),protected_hashes=len(protected))))
if __name__=='__main__':main()
