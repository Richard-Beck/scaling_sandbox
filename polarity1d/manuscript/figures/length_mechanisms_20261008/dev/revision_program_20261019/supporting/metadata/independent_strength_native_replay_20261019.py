"""Independent read-only native/regional replay; single CPU, immutable hash cache."""
from pathlib import Path
import argparse,collections,datetime,hashlib,importlib.util,json,sys
import numpy as np
P=argparse.ArgumentParser();P.add_argument('--run',required=True,type=Path);P.add_argument('--limit-new',type=int);a=P.parse_args();D=a.run.resolve();R=D.parents[2];B=R/'figures/length_mechanisms_20261008/dev/revision_program_20261019';M=B/'supporting/metadata';stage='baseline' if 'baseline' in D.name else 'strength';expected={'baseline':'f44a5e28ad4e6e1d3680d66806f33de530e9dea369fe21f99790cd7c6e23dc9b','strength':'62cda77004a549a9b3bdb8310ce733c0cb47321ed3b093e615a12e69ba524ff8'}[stage]
read=lambda p:json.loads(Path(p).read_text())
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for part in iter(lambda:f.read(8*1024*1024),b''):h.update(part)
 return h.hexdigest()
def write(p,z):
 q=p.with_suffix('.tmp.json');q.write_text(json.dumps(z,indent=2,allow_nan=False)+'\n');q.replace(p)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
m=read(D/'manifest.json');mh=sha(D/'manifest.json');assert mh==expected;pr=read(D/'provenance.json');assert pr['manifest_sha256']==mh
for p,h in pr['frozen_sources'].items():assert sha(D/'frozen_sources'/p)==h
for p,h in pr['frozen_inputs'].items():assert sha(D/p)==h
sys.dont_write_bytecode=True;sys.path.insert(0,str(D/'frozen_sources/code/mechanisms'));import revision_five_model_compensation_20261019 as engine
NAMES={'wave_pinning':['active','inactive'],'goryachev':['active','inactive'],'holmes_model3':['rac_active','rac_inactive','rho_active','rho_inactive']};DB=['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement']
cachepath=M/f'FIGURE2_COMPENSATION_{stage.upper()}_NATIVE_CASE_CACHE.json';cache=read(cachepath) if cachepath.exists() else dict(manifest_sha256=mh,cases={});assert cache['manifest_sha256']==mh;numeric=0;maximum=0.
def ck(a,b,k,tol=4e-10):
 global numeric,maximum
 a=np.asarray(a,float);b=np.asarray(b,float);assert a.shape==b.shape,(k,a.shape,b.shape);assert np.array_equal(np.isnan(a),np.isnan(b)),k;e=float(np.nanmax(abs(a-b))) if a.size and not np.isnan(a).all() else 0.;numeric+=a.size;maximum=max(maximum,e);assert e<tol,(k,e)
def load(p,h=None):
 if h:assert sha(p)==h,(p,'hash')
 with np.load(p) as z:arrays={k:z[k].copy() for k in z.files}
 for key,value in arrays.items():
  if key.startswith('field_'):assert np.isfinite(value).all(),(p,key)
 return arrays
def state(z,key,i):
 if key=='debelly':return {k:z['field_'+k][i] for k in DB}
 if key=='spring':
  y=np.empty(404);y[:-2:2]=z['field_active'][i];y[1:-2:2]=z['field_inactive'][i];y[-2:]=z['raw_boundaries'][i];return y
 return np.stack([z['field_'+k][i] for k in NAMES[key]])
def statefile(p,key,h):
 z=load(p,h);return {k:z['state_'+k] for k in DB} if key=='debelly' else z['state']
def statecmp(aa,bb):
 xs=list(aa.values()) if isinstance(aa,dict) else [aa];ys=list(bb.values()) if isinstance(bb,dict) else [bb];e=max(float(np.max(abs(x-y))) for x,y in zip(xs,ys));bound=32*np.finfo(float).eps*max(1,max(float(np.max(abs(y))) for y in ys));assert e<=bound,(e,bound);return e,bound
def monitor(z,key,L):
 out={k[6:]:v for k,v in z.items() if k.startswith('field_')}
 if key!='debelly':return out
 out={k:out[k] for k in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','boundary_separation']};out['displacement']=out['displacement']-out['displacement'].mean(axis=1)[:,None];out['physical_length']=L+out.pop('boundary_separation');return out
def event(t,good):
 hit=np.flatnonzero(good)
 if not len(hit):return dict(status='no_sampled_entry_through_finite_horizon',lower=float(t[-1]),upper=None)
 i=int(hit[0]);bad=np.flatnonzero(~good);last=int(bad[-1]+1) if len(bad) else 0
 return dict(status='already_present' if i==0 else 'sampled_entry',lower=float(t[max(0,i-1)]),upper=float(t[i]),retained_status='not_retained_through_last_sample' if last==len(t) else 'present_all_samples' if last==0 else 'sampled_retained_entry',retained_lower=None if last==len(t) else float(t[max(0,last-1)]),retained_upper=None if last==len(t) else float(t[last]))
def regional(z,c,g):
 key=c['model'];u=z['field_'+('rac' if key=='debelly' else 'rac_active' if key=='holmes_model3' else 'active')];t=z['time'];core=g['initial_core_mask'];q=g['material_fraction'];L=np.diff(z['raw_boundaries'],axis=1)[:,0] if key=='spring' else np.full(len(t),c['length'])
 if key=='spring':weights=np.broadcast_to(g['initial_weights']/float(g['initial_physical_length']),u.shape)*L[:,None];fp=(1-q[None,:])*L[:,None]<=c['width'];ori=np.trapezoid(u*(2*q-1),q,axis=1)/np.trapezoid(u,q,axis=1)
 else:
  weights=np.broadcast_to(g['initial_weights'],u.shape);fp=np.broadcast_to(g['initial_cue_footprint'],u.shape)
  if key=='debelly':ori=np.sum(u*(2*q-1),axis=1)/u.sum(axis=1);L=c['length']+z['field_boundary_separation'][:,0]
  else:up=np.maximum(u,0);ori=up@(2*z['x']/(z['x'][0]+z['x'][-1])-1)/up.sum(axis=1)
 ck(z['orientation'],ori,'orientation',3e-12);contrast=np.ptp(np.maximum(u,0),axis=1);ck(z['contrast'],contrast,'contrast',3e-12);initial=np.sum(g['initial_active']*core*g['initial_weights']);old=np.sum(u*weights*core,axis=1);overlap=np.sum(weights*core*fp,axis=1)/np.sum(weights*core,axis=1);peak=np.max(np.where(fp,u,-np.inf),axis=1)
 out=dict(time=t,original_core_active_amount=old,original_core_active_fraction=old/initial,actual_cue_region_active_peak=peak,actual_cue_region_active_amount=np.sum(u*weights*fp,axis=1),fixed_distal_material_region_active_amount=np.sum(u*weights*(q>=.5),axis=1),actual_cue_original_core_overlap_fraction=overlap,actual_domain_length=L,original_core_physical_or_nominal_extent=np.sum(weights*core,axis=1),global_orientation=ori,active_contrast=contrast)
 inactive='rac_inactive' if key=='holmes_model3' else 'inactive'
 if 'field_'+inactive in z:out['actual_cue_region_inactive_mean']=np.sum(z['field_'+inactive]*weights*fp,axis=1)/np.sum(weights*fp,axis=1)
 return out
new_this_execution=0
for c in m['cases']:
 cid=str(c['slot']);rr=c['exact_reuse'];F=D/'workers'/f'{c["slot"]:04d}'
 if rr:
  assert sha(R/rr['source_result'])==rr['source_result_sha256']
  if 'numerical_screen' in rr:assert sha(R/rr['numerical_screen'])==rr['numerical_screen_sha256']
  if cid not in cache['cases']:
   source=read(R/rr['source_result']);assert source['case']['model']==c['model'] and source['case']['length']==c['length'] and source['case']['checkpoint_sha256']==c['checkpoint_sha256'];outcomes=[]
   if c['case_kind']=='constant_on':assert sha(R/rr['cue_archive'])==rr['cue_archive_sha256']
   for w in rr['withdrawals']:
    for field,hkey in [('archive','archive_sha256'),('restart','restart_sha256')]:assert sha(R/w[field])==w[hkey]
    outcomes.append(dict(duration=w['duration'],orientation_success=w['orientation_success'],full_polarity_success=w['full_polarity_success'],field_tail_pass=w['field_tail_pass'],no_cue_preparation_valid=w['no_cue_preparation_valid']))
   prior=M/('FIGURE2_FIVE_MODEL_PILOT_V2_TERMINAL_INDEPENDENT_CHECK.json' if stage=='baseline' else 'FIGURE2_FIVE_MODEL_HIGH3_SCREEN_TERMINAL_INDEPENDENT_CHECK.json');pv=read(prior);assert pv['passed']
   for x in outcomes:
    br=next(y for y in pv['branch_outcomes'] if y['model']==c['model'] and y['length']==c['length'] and y['pulse_duration']==x['duration']);assert x['full_polarity_success']==br['full_polarity_success'] and x['field_tail_pass']==br['field_tail_pass']
   cache['cases'][cid]=dict(case=c,result_sha256=rr['source_result_sha256'],execution='exact_previously_audited_reuse',outcomes=outcomes,proof_sha256=sha(prior),native_frames_replayed=0)
  continue
 result=F/'result.json'
 if not result.exists():continue
 res=read(result);rh=sha(result);assert res['case']==c and res['manifest_sha256']==mh and read(F/'status.json')['state']=='complete'
 if cid in cache['cases']:assert cache['cases'][cid]['result_sha256']==rh;continue
 initial,prep=engine.n.read_checkpoint(D/c['frozen_checkpoint'],c['model']);g=load(D/c['regional_input'],c['regional_input_sha256']);on=None;caseframes=0;casepool=0.;casearchives={};events=[];outcomes=[];startnumeric=numeric
 if res['on']:
  w=res['on'];on=load(D/w['archive'],w['archive_sha256']);casearchives[w['archive']]=w['archive_sha256'];np.testing.assert_array_equal(on['time'],engine.sample_times(c['model'],max(c['pulse_durations'])));statecmp(state(on,c['model'],0),initial);reg=regional(on,c,g);saved=load(D/w['regional'],w['regional_sha256'])
  for key,v in reg.items():ck(saved[key],v,'ONregional_'+key)
  for threshold in [.9,.75,.5]:events.append(dict(phase='cue_on',observable='fixed_original_core_active_fraction',threshold=threshold,**event(on['time'],reg['original_core_active_fraction']<=threshold)))
  if np.any(g['initial_core_mask']&g['initial_cue_footprint']):events.append(dict(phase='cue_on',observable='separated_challenger_peak_above_initial_halfcontrast',threshold=float(g['initial_halfcontrast_threshold']),status='unavailable_initial_cue_core_overlap',lower=None,upper=None))
  else:events.append(dict(phase='cue_on',observable='separated_challenger_peak_above_initial_halfcontrast',threshold=float(g['initial_halfcontrast_threshold']),**event(on['time'],(reg['actual_cue_region_active_peak']>=g['initial_halfcontrast_threshold'])&(reg['actual_cue_original_core_overlap_fraction']==0))))
  caseframes+=len(on['time'])
 for j,w in enumerate(res['withdrawals']):
  key=c['model'];duration=w['duration'];restart=statefile(D/w['restart'],key,w['restart_sha256']);wanted=initial if duration==0 else state(on,key,int(np.flatnonzero(on['time']==duration)[0]));e,bound=statecmp(restart,wanted);assert e==0;off=load(D/w['archive'],w['archive_sha256']);casearchives[w['archive']]=w['archive_sha256'];np.testing.assert_array_equal(off['time'],engine.sample_times(key,10000.,off=True));e,bound=statecmp(state(off,key,0),restart);ck(e,w['saved_initial_max_absolute_roundoff'],'restartoutput');ck(bound,w['saved_initial_roundoff_bound'],'restartbound');caseframes+=len(off['time']);reg=regional(off,c,g);saved=load(D/w['regional'],w['regional_sha256'])
  for field,v in reg.items():ck(saved[field],v,'OFFregional_'+field);ck(w['regional_endpoint'][field],v[-1],'regionalendpoint') if field!='time' else None
  for threshold in [.9,.75,.5]:events.append(dict(duration=duration,phase='off',observable='fixed_original_core_active_fraction',threshold=threshold,**event(off['time'],reg['original_core_active_fraction']<=threshold)))
  if np.any(g['initial_core_mask']&g['initial_cue_footprint']):events.append(dict(duration=duration,phase='off',observable='separated_challenger_peak_above_initial_halfcontrast',threshold=float(g['initial_halfcontrast_threshold']),status='unavailable_initial_cue_core_overlap',lower=None,upper=None))
  else:events.append(dict(duration=duration,phase='off',observable='separated_challenger_peak_above_initial_halfcontrast',threshold=float(g['initial_halfcontrast_threshold']),**event(off['time'],(reg['actual_cue_region_active_peak']>=g['initial_halfcontrast_threshold'])&(reg['actual_cue_original_core_overlap_fraction']==0))))
  target=np.arange(0,10001.25,2.5);ix=np.searchsorted(off['time'],target);np.testing.assert_array_equal(off['time'][ix],target);ix=ix[target>=9500];o=reg['global_orientation'][ix];co=reg['active_contrast'][ix];os=bool(o.min()>.05 and co.min()>=.05 and np.ptp(o)<.005);er=float(np.max(abs(o/abs(c['initial_orientation'])-1)));fs=bool(os and er<=.01);fm=monitor(off,key,c['length']);maximumtail=0.;exc=0.;prefields={n:np.asarray(v) for n,v in prep.fields.items()}
  if key=='debelly':prefields['displacement']=prefields['displacement']-prefields['displacement'].mean(axis=1)[:,None];prefields['physical_length']=c['length']+prefields['boundary_separation']
  for name,v in fm.items():
   assert np.isfinite(v).all(),(cid,name);vv=v[ix];scale=c['fixed_scales'][name];variation=float(np.max(np.ptp(vv,axis=0))/scale);drift=float(np.max(abs(vv[-1]-vv[-2]))/2.5*500/scale);maximumtail=max(maximumtail,variation,drift);ck(variation,w['field_tail_metrics'][name]['range_scaled'],'tailrange');ck(drift,w['field_tail_metrics'][name]['projected_drift_500'],'taildrift')
   if duration==0:exc=max(exc,float(np.max(abs(v-prefields[name][-1]))/scale))
  for name,value in [('orientation_success',os),('full_polarity_success',fs),('field_tail_pass',maximumtail<=.001)]:assert value==w[name],(cid,duration,name)
  ck(o[-1],w['final_orientation'],'finalorientation');ck(er,w['maximum_relative_polarity_error_last_500'],'strengtherror');ck(maximumtail,w['field_tail_max_scaled'],'tailmax')
  valid=None
  if duration==0:valid=bool(exc<=.01 and not os and o.max()<-.05 and np.ptp(o)<.005 and co.min()>=.05);assert valid==w['no_cue_preparation_valid'];ck(exc,w['no_cue_maximum_scaled_excursion'],'zerocontrol')
  if key=='debelly':assert np.all(off['time']/.5==np.round(off['time']/.5))
  else:
   if key=='spring':pools=[off['field_length']*np.trapezoid(off['field_active']+off['field_inactive'],g['material_fraction'],axis=1)]
   else:pools=[sum(off['field_'+v] for v in pair).sum(axis=1) for pair in ([NAMES[key][:2],NAMES[key][2:]] if key=='holmes_model3' else [NAMES[key]])]
   pe=max(float(np.max(abs(v/v[0]-1))) for v in pools);assert pe<2e-5,(cid,duration,'pool',pe);casepool=max(casepool,pe)
  outcomes.append(dict(duration=duration,orientation_success=os,full_polarity_success=fs,field_tail_pass=maximumtail<=.001,field_tail_max_scaled=maximumtail,no_cue_preparation_valid=valid,regional_endpoint={key:float(v[-1]) for key,v in reg.items() if key!='time'}));del off,reg,saved,fm
 assert len(events)==len(res['events'])
 for aa,bb in zip(events,res['events']):
  for key,value in aa.items():
   if value is None or isinstance(value,str):assert bb[key]==value,(cid,key)
   else:ck(bb[key],value,'event_'+key)
 cache['cases'][cid]=dict(case=c,result_sha256=rh,execution='new_native',outcomes=outcomes,native_frames_replayed=caseframes,numeric_values_checked=numeric-startnumeric,maximum_pool_relative_error=casepool,native_archive_sha256=casearchives)
 write(cachepath,cache);print(json.dumps(dict(generated_at=now(),stage=stage,verified_groups=len(cache['cases']),requested_groups=len(m['cases']),last_slot=c['slot'],model=c['model'])),flush=True)
 if on is not None:del on
 new_this_execution+=1
 if a.limit_new and new_this_execution>=a.limit_new:break
write(cachepath,cache);missing=[c['slot'] for c in m['cases'] if str(c['slot']) not in cache['cases']];items=list(cache['cases'].values());outcomes=[dict(model=x['case']['model'],length=x['case']['length'],width=x['case']['width'],peak_multiplier=x['case']['native_peak_multiplier'],**w) for x in items for w in x['outcomes']];passednotmissing=not missing
report=dict(generated_at=now(),passed=True,native_scope_complete=passednotmissing,manifest_sha256=mh,stage=stage,requested_groups=len(m['cases']),groups_independently_verified=len(items),new_native_groups_verified=sum(x['execution']=='new_native' for x in items),reused_groups_preserved=sum(x['execution']!='new_native' for x in items),logical_branch_outcomes_verified=len(outcomes),pending_case_slots=missing,raw_frames_replayed=sum(x['native_frames_replayed'] for x in items),new_numeric_values_checked_this_execution=numeric,maximum_numeric_replay_error_this_execution=maximum,all_no_cue_controls_valid=(all(x['no_cue_preparation_valid'] for x in outcomes if x['peak_multiplier']==0) if any(x['peak_multiplier']==0 for x in outcomes) else None),no_cue_controls_evaluated=sum(x['peak_multiplier']==0 for x in outcomes),full_field_tail_failures=[{k:v for k,v in x.items() if k!='regional_endpoint'} for x in outcomes if not x['field_tail_pass']],cache_path=str(cachepath.relative_to(R)),cache_sha256=sha(cachepath),audit_script_sha256=sha(Path(__file__)),no_native_integrations=True,scope='All new native canonical arrays/restarts/regional sampled-readout arrays and endpoint/scalar/full-field gates are replayed. Reused canonical groups preserve exact hashes and reconcile to the earlier independent pilot/highscreen branch proof, without repeating those full native audits. Actual full table/art mapping and refinement are separate.',limitations=['Fixed original core is a material or nominal region, not site identity. Overlap disables a separated-challenger event.','Original deBelly eleven-state half-unit dynamics and nominal actuator coordinates remain distinct from physical endpoint/membrane motion.','Field-tail failures remain measured finite-horizon states. No optimized minimum or absent unobserved success island is inferred from five pulses.'])
out=M/f'FIGURE2_COMPENSATION_{stage.upper()}_NATIVE_PROGRESS_INDEPENDENT_CHECK.json';write(out,report);print(json.dumps(report,indent=2),flush=True)
