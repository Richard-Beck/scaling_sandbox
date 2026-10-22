"""Paired-length native field viewer for the frozen localized-patch campaign."""
from pathlib import Path
from datetime import datetime,timezone
import sys,gzip,json,hashlib,argparse
import numpy as np
SCRIPT=Path(__file__).resolve()
ROOT=SCRIPT.parents[2]
CODE=ROOT/'code/core'
HERE=ROOT/'data/reversal_fields'
OUTPUT=ROOT/'figures/diagnostic_reversal_fields'
TEMPLATE=SCRIPT.parent/'template.html'
CAMPAIGN=ROOT/'data/reversal_fields/original_campaign'
FROZEN=CODE
sys.path.insert(0,str(FROZEN))
import native as r
import settling
# Canonical copied adapters preserve the equations; immutable original frozen sources remain in data/width_reversal.
r.ROOT=ROOT
PAIRS={'wave_pinning':[9.,25.],'goryachev':[9.,21.],'debelly':[19.,35.],'holmes_model3':[14.,24.],'spring':[.8,2.4],'otsuji':[9.,25.],'legi':[9.,23.]}
LABELS={'wave_pinning':'Wave-pinning','goryachev':'Goryachev','debelly':'deBelly','holmes_model3':'Holmes 3','spring':'Spring','otsuji':'Otsuji','legi':'LEGI'}
CONDITIONS=['fractional_peak_matched','physical_peak_matched','physical_integral_matched']
JOBS=[(m,cond,L) for m,Ls in PAIRS.items() for cond in CONDITIONS for L in Ls]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def spring_observables():
 return {'tension_feedback':{'definition':'80 * (physical_length - resting_length_L0); delta = 3 + tension_feedback','label':'Global spring/tension feedback','units':'uncalibrated author deactivation-rate units','nature':'spatially uniform tension-coupled biochemical rate increment, not a calibrated mechanical force','frozen_source':str((FROZEN/'models/spring.py').relative_to(ROOT)),'source_sha256':sha(FROZEN/'models/spring.py'),'source_lines':[24,25,28,117,118]},'spring_restoring_term':{'definition':'0.01 * (physical_length - resting_length_L0)','label':'Drag-scaled spring restoring term','units':'uncalibrated author endpoint-velocity units','nature':'positive elastic restoring contribution at left endpoint, negative at right; no separate stiffness/drag force calibration','frozen_source':str((FROZEN/'models/spring.py').relative_to(ROOT)),'source_sha256':sha(FROZEN/'models/spring.py'),'source_lines':[24,25]}}
from functools import cache
@cache
def records():return json.load(gzip.open(CAMPAIGN/'records.json.gz','rt'))
def filename(m,cond,L):return f'{m}_{cond}_L{L:g}'
class Cue:
 """Frozen historical right patch; retained solely for the personal diagnostic."""
 def __init__(self,case,duration):self.case=case;self.stop=float(duration);self.event_times=(self.stop,)
 def __call__(self,t,x,L):
  x=np.asarray(x,float);c=self.case
  if t>=self.stop:return np.zeros_like(x)
  w=.2*(L if c['width_mode']=='fractional' else c['reference_domain_length'])
  if c.get('physical_patch_policy')=='restrict' and w>L/2+1e-8:raise ValueError('Physical patch exceeds half the instantaneous domain')
  peak=c['A']
  if c['dose_mode']=='integral_matched':
   if L>=w:peak*=L/c['reference_domain_length']
   else:peak*=.1*L/(.5*(L+w/np.pi*np.sin(np.pi*L/w)))
  distance=(L-x)/w
  return np.where((distance>=0)&(distance<=1),peak*.5*(1+np.cos(np.pi*np.clip(distance,0,1))),0.)

def reconstruct(slot):
 m,cond,L=JOBS[slot];c=next(c for c in records()['cases'] if c['model']==m and c['length']==L and c['condition']==cond)
 assert c['eligible'] and not c['restricted']
 selection=json.loads((ROOT/'data/reversal_fields'/(filename(m,cond,L)+'.json')).read_text())
 definition=selection['selected_definition'];D=float(selection['duration']);test=next(t for t in c['tests'] if t['duration']==D)
 checkpoint=ROOT/'data/reversal_fields/reconstruction_inputs'/c['checkpoint'];assert sha(checkpoint)==c['checkpoint_sha256'];state,prep=r.read_checkpoint(checkpoint,m)
 cue=Cue(c,D);pulse=r.run(m,L,state,cue,D,interval=min(2.5,D/20),n=c['cells'])
 recovery=r.run(m,L,pulse.state,r.ConstantInput(0.),10000,interval=2.5,n=c['cells'])
 keep=recovery.time>=9500;o=recovery.orientation[keep];contrast=recovery.contrast[keep]
 orientation=bool(np.ptp(o)<.005 and np.min(contrast)>=.05 and np.min(o)>.05)
 full_error=float(np.max(abs(o/abs(c['initial_orientation'])-1)))
 full=bool(orientation and full_error<=.01)
 assert orientation==test['success'] and full==test['full_success']
 np.testing.assert_allclose(o[-1],test['final_orientation'],rtol=1e-7,atol=1e-8)
 monitored=settling.monitored(recovery,m,L);tailmetrics={}
 for name,y in monitored.items():
  scale=c['fixed_scales'][name];tail=y[keep]
  tailmetrics[name]={'range_scaled':float(np.max(np.ptp(tail,axis=0))/scale),'projected_drift_500':float(np.max(abs(tail[-1]-tail[-2]))/2.5*500/scale)}
 fieldtail=max(max(v.values()) for v in tailmetrics.values())
 assert (fieldtail<=.001)==test['tail_field_pass_0_1_percent']
 np.testing.assert_allclose(fieldtail,test['tail_field_max_scaled'],rtol=1e-6,atol=1e-10)
 # Check every originally retained spatial field, and original full states, not just the observable.
 trace=ROOT/'data/reversal_fields/reconstruction_inputs'/test['trace'];prefix=f'trial{test["trace_row"]}_';errors={}
 with np.load(trace) as z:
  for phase,tr in [('pulse',pulse),('recovery',recovery)]:
   expected_times=z[prefix+phase+'_time'];ids=np.searchsorted(tr.time,expected_times)
   np.testing.assert_allclose(tr.time[ids],expected_times,rtol=0,atol=1e-9)
   for name,values in tr.fields.items():
    expected=z[prefix+phase+'_field_'+name]
    np.testing.assert_allclose(values[ids],expected,rtol=1e-7,atol=1e-8)
    errors[phase+'_'+name]=float(np.max(abs(values[ids]-expected)))
   expected=z[prefix+phase+'_active'];np.testing.assert_allclose(tr.active[ids],expected,rtol=1e-7,atol=1e-8)
  np.testing.assert_allclose(o,z[prefix+'recovery_tail_orientation'],rtol=1e-7,atol=1e-8)
 np.testing.assert_allclose(pulse.active[0],prep.active[-1],rtol=1e-9,atol=1e-10)
 tr=r.concatenate(pulse,recovery,D)
 fields={k:np.asarray(v) for k,v in tr.fields.items()}
 applied=[]
 for i,t in enumerate(tr.time):
  domain=float(fields['length'][i]) if m=='spring' else L
  x=tr.x/L*domain if m=='spring' else tr.x
  if m=='debelly' and r.model_for(m).parameters.input_sampling_mode=='julia_element_right_edge':x=np.minimum(x+.5*L/c['cells'],L)
  applied.append(cue(float(t),x,domain))
 fields['applied_patch']=np.asarray(applied)
 if m=='spring':
  fields['tension_feedback']=80.0*(fields['length']-L)
  fields['spring_restoring_term']=0.01*(fields['length']-L)
 assert all(np.isfinite(v).all() for v in fields.values())
 file=filename(m,cond,L)
 np.savez_compressed(HERE/(file+'_full.npz'),time=tr.time,x=tr.x,orientation=tr.orientation,contrast=tr.contrast,**{'field_'+k:v for k,v in fields.items()})
 meta=dict(model=m,label=LABELS[m],condition=cond,length=L,duration=D,amplitude=c['A'],reference_length=c['reference_length'],reference_domain_length=c['reference_domain_length'],initial_domain_length=c['initial_domain_length'],physical_patch_policy=c['physical_patch_policy'],selected_definition=definition,selected_orientation_pass=orientation,selected_full_polarity_pass=full,selected_full_polarity_relative_error=full_error,selected_tail_field_pass=fieldtail<=.001,selected_tail_field_max_scaled=fieldtail,tail_field_metrics=tailmetrics,morphology=test['morphology'],initial_orientation=c['initial_orientation'],final_orientation=test['final_orientation'],results=c['result'],cells=c['cells'],case_index=c['case_index'],canonical_index=c['canonical_index'],full_frames=len(tr.time),checks={'checkpoint_sha256':sha(checkpoint),'all_original_native_field_max_error':max(errors.values()),'original_native_field_errors':errors,'orientation_gate_matches':True,'full_polarity_gate_matches':True,'field_tail_flag_matches':True},sources={str(path.relative_to(ROOT)):sha(path) for path in [checkpoint,trace,CAMPAIGN/'records.json.gz',CAMPAIGN/'manifest.json']})
 if m=='spring':meta['derived_observables']=spring_observables()
 (HERE/(file+'.json')).write_text(json.dumps(meta,separators=(',',':')))
 print(json.dumps({'model':m,'condition':cond,'length':L,'duration':D,'orientation':orientation,'full':full,'fieldtail':fieldtail<=.001,'all_native_field_error':max(errors.values())}),flush=True)
def interpolate(times,values,requested):
 a=np.asarray(values);flat=a.reshape(len(a),-1)
 y=np.stack([np.interp(requested,times,v) for v in flat.T],axis=1).reshape((len(requested),)+a.shape[1:])
 return y

def bundle():
 experiments=[]
 for m,Ls in PAIRS.items():
  for cond in CONDITIONS:
   trajectories=[]
   for L in Ls:
    name=filename(m,cond,L);meta=json.loads((HERE/(name+'.json')).read_text())
    with np.load(HERE/(name+'_full.npz')) as z:raw={k:z[k].copy() for k in z.files}
    trajectories.append((meta,raw))
   cutoffs=[meta['duration'] for meta,raw in trajectories];ends=[raw['time'][-1] for meta,raw in trajectories]
   # Common actual time; dense pulse coverage for both separate cutoffs, logarithmic recovery and explicit gates.
   wanted=np.unique(np.r_[*[np.linspace(0,D,41) for D in cutoffs],*[D+np.r_[0,np.geomspace(2.5,10000,95),9500,10000] for D in cutoffs]])
   time=wanted[(wanted>=0)&(wanted<=max(ends))]
   facets=[];scalars=[];prepared=[]
   names=[k[6:] for k in trajectories[0][1] if k.startswith('field_')]
   for meta,raw in trajectories:
    valid=time<=raw['time'][-1]+1e-9
    prepared.append({'valid':valid,'arrays':{k:interpolate(raw['time'],v,time[valid]) for k,v in raw.items() if k!='x'},'x':raw['x']/meta['length']})
   for name in names:
    vs=[raw['field_'+name] for meta,raw in trajectories]
    scalar=all(v.ndim==1 or (name!='applied_patch' and np.all(np.ptp(v,axis=1)==0)) for v in vs)
    lo=min(float(v.min()) for v in vs);hi=max(float(v.max()) for v in vs);margin=max((hi-lo)*.06,max(abs(lo),abs(hi))*1e-5,1e-8)
    entry={'name':name,'range':[lo-margin,hi+margin],'series':[]}
    for j,(meta,raw) in enumerate(trajectories):
     item=prepared[j];vals=item['arrays']['field_'+name]
     if scalar and vals.ndim>1:vals=vals[:,0]
     # The patch is evaluated exactly at common time (not interpolated across cue withdrawal), and follows the current Spring domain.
     if name=='applied_patch':
      vals=[];case=next(c for c in records()['cases'] if c['case_index']==meta['case_index']);cue=Cue(case,meta['duration'])
      for k,t in enumerate(time[item['valid']]):
       domain=float(item['arrays']['field_length'][k]) if m=='spring' else meta['length'];x=item['x']*domain if m=='spring' else raw['x']
       if m=='debelly' and r.model_for(m).parameters.input_sampling_mode=='julia_element_right_edge':x=np.minimum(x+.5*meta['length']/meta['cells'],meta['length'])
       vals.append(cue(float(t),x,domain))
      vals=np.asarray(vals)
     rounded=np.vectorize(lambda q:float(f'{q:.9g}'))(vals).tolist()
     aligned=[];cursor=0
     for valid in item['valid']:
      if valid:aligned.append(rounded[cursor]);cursor+=1
      else:aligned.append(None)
     initial=raw['field_'+name][0]
     if scalar and np.ndim(initial)>0:initial=float(initial[0])
     entry['series'].append({'position':item['x'].tolist(),'values':aligned,'initial':np.asarray(initial).tolist()})
    (scalars if scalar else facets).append(entry)
   ors=[];cs=[]
   for j,(meta,raw) in enumerate(trajectories):
    valid=prepared[j]['valid'];ors.append([float(np.interp(t,raw['time'],raw['orientation'])) if flag else None for t,flag in zip(time,valid)]);cs.append([float(np.interp(t,raw['time'],raw['contrast'])) if flag else None for t,flag in zip(time,valid)])
   experiments.append({'model':m,'label':LABELS[m],'condition':cond,'id':m+'_'+cond,'time':time.tolist(),'lengths':[meta for meta,raw in trajectories],'orientation':ors,'contrast':cs,'facets':facets,'scalars':scalars})
 # Keep facet axes consistent across all three cue conditions, both lengths and all native times.
 for model in PAIRS:
  es=[ex for ex in experiments if ex['model']==model]
  for category in ['facets','scalars']:
   for name in [f['name'] for f in es[0][category]]:
    fs=[next(f for f in ex[category] if f['name']==name) for ex in es]
    common=[min(f['range'][0] for f in fs),max(f['range'][1] for f in fs)]
    for f in fs:f['range']=common
 stamp=datetime.now(timezone.utc).isoformat();payload={'generated_at':stamp,'campaign':'Historical localized-patch reversal, linear-v11 prepared; personal diagnostic; data/reversal_fields/original_campaign','experiments':experiments}
 compressed=gzip.compress(json.dumps(payload,separators=(',',':')).encode())
 import base64
 OUTPUT.mkdir(parents=True,exist_ok=True)
 (OUTPUT/'index.html').write_text(TEMPLATE.read_text().replace('__PAYLOAD__',base64.b64encode(compressed).decode()))
 sources={str(path.relative_to(ROOT)):sha(path) for path in [SCRIPT,TEMPLATE,FROZEN/'protocol.json',*FROZEN.glob('*.py'),*FROZEN.glob('models/*.py'),*(FROZEN/'src/polarity1d').glob('*.py')]}
 provenance={'generated_at':stamp,'campaign':payload['campaign'],'pairs':PAIRS,'conditions':CONDITIONS,'spring_derived_observables':spring_observables(),'spatial_mass_median':{'definition':'Coordinate bisecting the trapezoidal integral of the current nonnegative field on the actual plotted native grid; exact integral inversion inside the piecewise-linear segment.','spatial_axis':'x/nominal L; Spring material fraction. Uniform per-frame coordinate Jacobian cancels. No extension beyond first/last plotted spatial sample.','negative_policy':'Any value below -(1e-12 + 1e-10 * max_abs_value) makes x50 undefined. Only negatives within that roundoff tolerance are set to zero for this diagnostic; native displayed and saved profiles remain unchanged. No absolute-value conversion.','negligible_integral':'x50 undefined unless total > 1e-12 * plotted_coordinate_span * max(1,max_abs_value). Zero/off-cue and unavailable profiles have no marker.','marker':'Per-length blue/orange dotted verticals with current numeric facet readouts; starting profiles remain dashed, and global scalar histories get no spatial median.'},'experiments':[{k:v for k,v in ex.items() if k not in ['facets','scalars','orientation','contrast','time']} for ex in experiments],'frozen_sources':sources,'sampling':'Common REAL time from simultaneous pulse onset. Native state fields pointwise linearly interpolated between full native solver output samples; exact cutoff and observation-end times included. Applied patch evaluated directly at common time, with interpolated actual Spring length. Profiles become unavailable beyond each recorded observation end; no extrapolation or held future states. Native grid points retained; x/nominal L, Spring material coordinate. No concentration normalization. Fixed axes span both complete trajectories.','archive':'../reversal_field_report_linear_wrong_campaign_archive_20261003 is the numerically faithful but wrong-campaign first report.','outputs':{path.name:sha(path) for path in [OUTPUT/'index.html',*HERE.glob('*_full.npz')]}}
 provenance['retention']={'input_manifest':'data/provenance/data_sources.json','numerical_inputs':'Retained full native fields and original records, unchanged','working_sources':'Relocated canonical core; original frozen campaign sources preserved separately'}
 (OUTPUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 print('Bundled',len(experiments),'paired experiments; HTML bytes',(OUTPUT/'index.html').stat().st_size,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--slot',type=int);ap.add_argument('--bundle',action='store_true');args=ap.parse_args()
 if args.bundle:bundle()
 else:
  if args.slot is None:ap.error('Choose --bundle or --slot 0..41; reconstruction requires SLURM')
  if not 0 <= args.slot < len(JOBS):ap.error('Slot must be in 0..41')
  HERE=ROOT/'data/reconstructed_reversal_fields';HERE.mkdir(parents=True,exist_ok=True)
  reconstruct(args.slot)
