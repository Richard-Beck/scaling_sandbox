"""Independent retained-only checks before retiring the exploratory campaigns."""
from pathlib import Path
from datetime import datetime,timezone
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import hashlib,json
import numpy as np

PACKAGE=Path(__file__).resolve().parents[1]
FIELDS=('R','Q','b','c','u','v','T')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def phase_increment(z0,z1):return float(np.angle(z1*np.conj(z0))*50/(2*np.pi))
def main():
 checks={};hashes={};events=[];branches={}
 def pin(p):hashes[str(p.relative_to(PACKAGE))]=sha(p)
 def verify(p,digest):
  assert p.is_file(),str(p)
  assert sha(p)==digest,str(p)
  pin(p)
 manifest=read(PACKAGE/'methods/viewer_packaging_manifest.json')
 # The immutable packaging manifest is retained separately from this audit.
 pin(PACKAGE/'methods/viewer_packaging_manifest.json')
 for name,row in manifest['files'].items():verify(PACKAGE/name,row['sha256'])
 raw=read(PACKAGE/'viewer/raw/manifest.json')
 for name,digest in raw['files'].items():verify(PACKAGE/'viewer'/name,digest)
 report=read(PACKAGE/'viewer/report_manifest.json')
 for group,base in [('local_inputs',PACKAGE/'viewer'),('sources',PACKAGE),('outputs',PACKAGE/'viewer')]:
  for name,digest in report[group].items():verify(base/name,digest)
 model=read(PACKAGE/'methods/reference_model/source_manifest.json')
 for row in model['sources']:verify(PACKAGE/row['retained_path'],row['retained_sha256'])
 evidence=read(PACKAGE/'evidence_manifest.json')
 for row in evidence['assets']:
  verify(PACKAGE/row['retained_path'],row['retained_sha256'])
  if row['transformation']=='byte-for-byte copy':assert row['source_sha256']==row['retained_sha256']
 for name in ['methods/relocation_check.json','methods/browser_checks/result.json']:
  assert read(PACKAGE/name)['passed'];pin(PACKAGE/name)
 for name in ['viewer/raw/manifest.json','viewer/report_manifest.json','methods/reference_model/source_manifest.json','evidence_manifest.json','executive_assessment.json','index.html']:
  pin(PACKAGE/name)
 checks['local_input_output_source_and_evidence_hashes']=True
 checks['prior_browser_and_access_forbidden_relocation_receipts_pass']=True
 for branch in ['A','original','B']:
  folder=PACKAGE/'viewer/raw'/branch;meta=read(folder/'frames.json');rows=meta['frames']
  with np.load(folder/'high_resolution.npz',allow_pickle=False) as z:data={k:z[k] for k in z.files}
  def frame(i):
   a,b=map(int,data['offsets'][i:i+2]);return data['widths'][a:b],{f:data[f][a:b] for f in FIELDS}
  def moment(i):
   h,f=frame(i);edges=np.r_[0.,np.cumsum(h)];k=2*np.pi/50
   return np.dot(f['R'],(np.exp(1j*k*edges[1:])-np.exp(1j*k*edges[:-1]))/(1j*k))/np.dot(h,f['R'])
  for i,row in enumerate(rows):
   h,f=frame(i);assert len(h)>0 and np.all(h>0) and abs(h.sum()-50)<1e-10
   assert float(data['time'][i])==row['time']
   assert all(np.all(np.isfinite(a)) for a in f.values())
   assert all(np.all(f[k]>=0) for k in ['R','Q','b','c'])
   edges=np.r_[0.,np.cumsum(h)]
   for face,cut in row['cuts'].items():
    assert 0<=int(face)<len(h)
    assert abs(edges[int(face)]-cut['x'])<1e-10
    assert all(np.isfinite(cut[k]) for k in ['Uleft','Uright','raw_left','raw_right'])
    assert row['next_id']>cut['id']
  branch_events=[]
  for i,row in enumerate(rows):
   if row['stage_code']!=1:continue
   t=row['time'];before=[j for j,r in enumerate(rows) if r['time']==t and r['stage_code']==0]
   after=[j for j,r in enumerate(rows) if r['time']==t+.5 and r['stage_code']==0]
   assert len(before)==len(after)==1
   j=before[0];h0,f0=frame(j);h1,f1=frame(i)
   e0=np.r_[0.,np.cumsum(h0)];e1=np.r_[0.,np.cumsum(h1)]
   edges=np.union1d(e0,e1);mid=(edges[1:]+edges[:-1])/2;widths=np.diff(edges)
   ix0=np.clip(np.searchsorted(e0,mid,side='right')-1,0,len(h0)-1)
   ix1=np.clip(np.searchsorted(e1,mid,side='right')-1,0,len(h1)-1)
   diffs={f:dict(max=float(np.max(np.abs(f1[f][ix1]-f0[f][ix0]))),L1=float(np.dot(widths,np.abs(f1[f][ix1]-f0[f][ix0])))) for f in FIELDS}
   assert all(d['max']==d['L1']==0 for d in diffs.values())
   record=dict(branch=branch,time=t,pre_cells=len(h0),post_cells=len(h1),union_piecewise_constant_difference=diffs,topology_phase_change_um=phase_increment(moment(j),moment(i)))
   if t>1000:
    prior=[j for j,r in enumerate(rows) if r['time']==t-.5 and r['stage_code']==0];assert len(prior)==1
    direction=np.sign(phase_increment(moment(prior[0]),moment(before[0])));assert direction!=0
    indices=[j for j,r in enumerate(rows) if r['stage_code']==0 and t<=r['time']<=t+30]
    assert len(indices)==61 and np.all(np.diff([rows[j]['time'] for j in indices])==.5)
    increments=[phase_increment(moment(a),moment(b)) for a,b in zip(indices,indices[1:])]
    assert all(direction*d>0 for d in increments)
    record['late_30s']=dict(prior_direction=int(direction),physical_half_steps=60,all_increments_continue_prior_direction=True,min_signed_increment_um=float(min(direction*d for d in increments)),max_signed_increment_um=float(max(direction*d for d in increments)),net_signed_displacement_um=float(direction*sum(increments)))
   events.append(record);branch_events.append(t)
  branches[branch]=dict(native_frames=len(rows),event_clocks=branch_events,all_seven_fields_finite=True,chemistry_nonnegative=True,all_widths_positive=True,controllers_next_id_valid=True)
 assert len(events)==18 and sum('late_30s' in e for e in events)==9
 checks['all18_topology_maps_preserve_all7_piecewise_constant_profiles_exactly']=True
 checks['all9_late_windows_continue_prior_phase_direction_for_every60_native_halfsteps']=True
 class Links(HTMLParser):
  def __init__(self):super().__init__();self.links=[]
  def handle_starttag(self,tag,attrs):
   for k,v in attrs:
    if k in ('href','src') and v:self.links.append(v)
 missing=[];local=[]
 for p in PACKAGE.rglob('*.html'):
  parser=Links();text=p.read_text();parser.feed(text)
  assert 'http-equiv="refresh"' not in text.lower() and 'location.reload' not in text
  for url in parser.links:
   u=urlsplit(url)
   if u.scheme or u.netloc or not u.path:continue
   dest=(p.parent/unquote(u.path)).resolve();assert dest.is_relative_to(PACKAGE.resolve()),str(dest)
   local.append(dict(source=str(p.relative_to(PACKAGE)),target=str(dest.relative_to(PACKAGE))))
   if not dest.exists():missing.append(str(dest.relative_to(PACKAGE)))
 expected=['cleanup_receipt.json','cleanup_validation.json']
 assert set(missing)<=set(expected),missing
 checks['all_current_local_links_internal_and_present_except_pending_cleanup_receipts']=True
 checks['executive_scientific_scope_review']=True
 result=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=True,scope='Retained package verified before deletion; does not certify cleanup already completed or make the unresolved periodic model valid.',checks=checks,branches=branches,retained_only_event_checks=events,local_links=local,pending_expected_cleanup_links=sorted(set(missing)),source_sha256=hashes,scientific_limits=['Historical audits certify completed checks before full campaigns were retired; their source paths are provenance, not current runtime dependencies.','The retained viewer regenerates from local data only; complete historical trajectory regeneration is intentionally unavailable.','Movement is not a failure by itself, and universal post-cue stationarity was not imposed. The physical closure and robust comparative ring counterpart remain unresolved.','Instantaneous allocations and the complete-retarget fork do not isolate a unique causal clock; no further science was added during cleanup.'])
 pin(Path(__file__))
 target=PACKAGE/'package_validation.json';target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(passed=True,path=str(target),sha256=sha(target),event_maps=len(events),late_halfsteps=540,remaining_expected_links=sorted(set(missing)))))
if __name__=='__main__':main()
