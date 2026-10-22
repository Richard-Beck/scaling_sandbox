"""Round-2 Figure 2, including independently audited native Holmes retention."""
from pathlib import Path
import csv, json, hashlib, sys
from datetime import datetime,timezone
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
CORE_PREVIEW='--core-preview' in sys.argv
OUT=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure2'
OLD=OUT.parents[1]/'signature_assays_20261020'
RUN=ROOT/'polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022'
sys.path.insert(0,str(ROOT/'polarity1d/manuscript/code/mechanisms'))
from signature_analysis_20261020 import replacement_time_ratio
from signature_figure2_supplement_20261020 import display_indices, exact_binomial_interval
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader,PdfWriter
ORDER=('goryachev','wave_pinning','holmes_model3','debelly','spring')
NAMES=dict(goryachev='Goryachev',wave_pinning='Wave-pinning',holmes_model3='Holmes',debelly='deBelly',spring='Spring')
COLORS=dict(wave_pinning='#2878b5',goryachev='#c83b3b',debelly='#82599c',holmes_model3='#36944c',spring='#c59a25')
LENGTHS=(19,23,27,31,35,39);GAPS=(0,.1,.5,1,4,16)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def csvout(name,rows):
 destination=OUT/'supporting/core_preview' if CORE_PREVIEW else OUT
 destination.mkdir(parents=True,exist_ok=True)
 with (destination/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def clean(ax):ax.spines[['top','right']].set_visible(False)
for d in ('supporting','supporting/data','supporting/metadata'):(OUT/d).mkdir(exist_ok=True)
coverage=read(OLD/'supporting/metadata/coverage.json')
assert coverage['final_export_eligible'] and not coverage['audit_errors']
identities=read(OLD/'supporting/metadata/artifact_sha256.json')
for rel in ('supporting/data/width_ratio_sampled_bounds.json','supporting/data/selected_train_endpoints.json','supporting/vectors/figure2_signature_assays.pdf'):
 assert sha(OLD/rel)==identities[rel]
old_display=read(OLD/'supporting/metadata/S1_S2_display_audit.json')
rel='supporting/data/figure2_S2_all_native_pulse_gap_endpoints.json'
assert sha(OLD/rel)==old_display['artifact_sha256'][rel]
supplement_proof=read(OLD/'supporting/metadata/combined_supplement_audit.json')
assert supplement_proof['passed'] and sha(supplement_proof['pdf']['path'])==supplement_proof['pdf']['sha256']
definition=read(OLD/'supporting/metadata/reporter_initiation_definition.json')
reporter_dir=Path(definition['path']);audit=read(reporter_dir/'reporter_sensitivity_audit.json')
assert sha(reporter_dir/'reporter_sensitivity_audit.json')==definition['audit_sha256']
for name,digest in audit['output_sha256'].items():assert sha(reporter_dir/name)==digest
reporters=read(reporter_dir/'all450_fixed_reporter_events.json')
orgfile=Path(coverage['retained_organization']['path']);assert sha(orgfile)==coverage['retained_organization']['sha256']
org=list(csv.DictReader(orgfile.open()))
ratios=read(OLD/'supporting/data/width_ratio_sampled_bounds.json')
trains=read(OLD/'supporting/data/figure2_S2_all_native_pulse_gap_endpoints.json')
examples=read(OLD/'supporting/data/selected_train_endpoints.json')
assert len(trains)==1458 and len(reporters)==450
debelly_proof=read(OLD/'supporting/metadata/S5_debelly_temporal_display_audit.json')
assert debelly_proof['passed'] and sha(debelly_proof['source_summary'])==debelly_proof['source_summary_sha256']
debelly_numerics=read(debelly_proof['source_summary'])
debelly_comparison=next(r for r in debelly_numerics['records'] if r['assay']=='train_4_gap4_on1p5' and r['requested_reference_length']==39)
assert debelly_comparison['companion_no_cue_control_pass']
assert all(r['classification_agrees'] and r['main_challenger_bias_lost'] and r['companion_challenger_bias_lost'] for r in debelly_comparison['first_two_gap_bias_loss'])
assert debelly_comparison['main201_final']['final_d']>=.5 and debelly_comparison['companion_final_D']>=.5
csvout('debelly_numerical_qualification.csv',[dict(endpoint=e['name'],time=e['time'],main_step_D=e['main201']['D_half'],quarter_step_D=e['companion']['D_half'],bias_sign_agrees=e['challenger_bias_sign_agrees'],within_D_tolerance=e['within_D_tolerance']) for e in debelly_comparison['phase_comparison']['endpoints']])
if CORE_PREVIEW:
 rows=[read(RUN/'results/independent_audit/workers'/f'{task:04d}.json') for task in range(18)]
 assert all(r['passed'] for r in rows)
 proof=dict(passed=True,workers=18,all_native_phase_join_audits_pass=True,manifest_sha256=sha(RUN/'manifest.json'),trajectories=rows)
else:proof=read(RUN/'results/independent_native_readout_audit.json')
assert proof['passed'] and proof['workers']==(18 if CORE_PREVIEW else 32) and proof['all_native_phase_join_audits_pass']
assert proof['manifest_sha256']==sha(RUN/'manifest.json')
new=proof['trajectories']
manifest=read(RUN/'manifest.json')
assert sha(manifest['reference_calibration']['path'])==manifest['reference_calibration']['sha256']
calibrations=read(manifest['reference_calibration']['path'])['calibrations']
for length in (31,35,39):
 calibration=next(r['calibration'] for r in calibrations if r['calibration_id']==f'holmes_model3:L{length}')
 assert not calibration['train_ready'] and calibration['upper_committing_duration'] is None
 assert calibration['finite_duration_horizon']==100000
assert {(r['case']['length'],r['case']['gap']) for r in new if r['case']['kind']=='train'}=={(length,gap) for length in (31,35,39) for gap in GAPS}
if not CORE_PREVIEW:assert len({r['case']['seed'] for r in new if r['case']['kind']=='noise'})==12
for row in new:
 assert sha(row['result_path'])==row['result_sha256']
holmes=next(r for r in new if r['case']['kind']=='train' and r['case']['length']==39 and r['case']['gap']==4)
examples=[r for r in examples if r['model']!='holmes_model3']+[dict(model='holmes_model3',displayed_length=39,source=dict(run_root=str(RUN),task=holmes['task']),schedule=holmes['trial']['schedule'],final_d=holmes['trial']['final_d'],field_tail_pass=holmes['trial']['field_tail_pass'],endpoints=holmes['endpoints'])]
# Ratios retain both measurement intervals; a zero lower denominator gives no finite upper bound.
init=[]
for model in ORDER:
 selected=sorted([r for r in reporters if r['model']==model and r['area']==1 and r['width']==3.8],key=lambda r:r['length'])
 ref=selected[0]['fixed_reporter_events']['0.1']
 for r in selected:
  event=r['fixed_reporter_events']['0.1'];ratio=replacement_time_ratio(event,ref)
  if r['length']==19:ratio.update(lower=1.,upper=1.,estimate_from_interpolation=1.,status='reference_identity')
  unresolved=ref['lower_elapsed']==0 and r['length']!=19
  init.append(dict(model=model,length=r['length'],time_lower=event['lower_elapsed'],time_upper=event['upper_elapsed'],ratio_lower=ratio['lower'],ratio_upper=ratio['upper'],ratio_estimate=None if unresolved else ratio.get('estimate_from_interpolation'),status='normalization_unresolved' if unresolved else ratio['status']))
csvout('initiation_normalized.csv',init)
csvout('equal_area_width_ratios.csv',[dict(model=r['model'],length=r['length'],area=r['area'],lower=r['ratio'].get('lower'),upper=r['ratio'].get('upper'),estimate=r['ratio'].get('estimate_from_interpolation'),informative=r['ratio'].get('informative'),status=r['ratio']['status']) for r in ratios])
fates=[];endrows=[]
for r in trains:
 c=r['case']
 if c['pulse_count']!=4 or c['total_on_factor_of_calibration']!=1.5:continue
 fates.append(dict(model=c['model'],length=c['requested_reference_length'],gap=c['gap_multiple_of_pulse'],value=r['final_d'],eligible=r['readout_eligible'],measurement='post-recovery D',tail_pass=r['field_tail_pass'],status='observed'))
 for e in r['endpoints']:
  if e['kind'] in ('on','gap'):
   endrows.append(dict(model=c['model'],length=c['requested_reference_length'],gap=c['gap_multiple_of_pulse'],kind=e['kind'],time=e['end']['time'],D=e['end']['D_half'],share=e['end']['share_half'],retention_ratio=e.get('primary_retention_ratio')))
holmes_endrows=[]
for row in new:
 c=row['case'];trial=row['trial'];gaps=[r for r in row['endpoints'] if r['kind']=='gap']
 assert len(gaps)==3
 if c['kind']=='train':
  fates.append(dict(model='holmes_model3',length=c['length'],gap=c['gap'],value=gaps[-1]['end']['D_half'],eligible=gaps[-1]['end']['readout_eligible'],measurement='last inter-pulse gap-end D',tail_pass=trial['field_tail_pass'],status='observed_retention_without_target_calibration'))
 for e in row['endpoints']:
  if e['kind'] not in ('on','gap'):continue
  endpoint=dict(model='holmes_model3',length=c['length'],gap=c['gap'],kind=e['kind'],time=e['end']['time'],D=e['end']['D_half'],share=e['end']['share_half'],retention_ratio=e.get('primary_retention_ratio'))
  if c['kind']=='train':endrows.append(endpoint)
  holmes_endrows.append(dict(task=c['task'],forcing_kind=c['kind'],seed=c['seed'],**endpoint,pulse_D=e['start']['D_half'],pulse_share=e['start']['share_half'],gained_share_retention=e.get('retention_ratio'),S3_pulse_eligible=e['primary_retention_pulse_eligible'],S3_denominator_eligible=e['primary_retention_denominator_eligible'],S3_success=e['primary_retention_success'],final_D=trial['final_d'],tail_pass=trial['field_tail_pass']))
assert len(fates)==180
csvout('train_fates_and_retention.csv',fates);csvout('all_selected_train_endpoints.csv',endrows);csvout('holmes_large_L_endpoints.csv',holmes_endrows)
# Numeric evidence for every model x four signature cells; unresolved claims stay explicit.
evidence=[]
for m in ORDER:
 i=[r for r in init if r['model']==m and r['length']==39][0]
 w=[r for r in ratios if r['model']==m and r['length'] in (19,39)]
 e=next(r for r in examples if r['model']==m);g=[r for r in e['endpoints'] if r['kind']=='gap'];first=g[0]
 raw=first['end']['share_half']/first['start']['share_half']
 evidence.extend([
 dict(model=m,signature='initiation slows with L',lower=i['ratio_lower'],upper=i['ratio_upper'],value=i['ratio_estimate'],qualification='deBelly first-update intervals do not resolve a normalized size effect' if m=='debelly' else 'L39/L19, central area, fixed 10% L19 reporter'),
 dict(model=m,signature='narrow-cue penalty and growth',lower=w[-1]['ratio'].get('lower'),upper=w[-1]['ratio'].get('upper'),value=w[-1]['ratio'].get('estimate_from_interpolation'),qualification='L19 '+str(w[0]['ratio'])+'; L39 '+str(w[-1]['ratio'])+('; Holmes growth verified L19 ~29.03 to L23 >133.16; growth through L39 unresolved, not an assigned infinite ratio' if m=='holmes_model3' else '; modest non-growing penalty is distinct from the large saturating-model penalty' if m in ('spring','goryachev') else '')),
 dict(model=m,signature='share held through gaps',lower=None,upper=None,value=raw,qualification=f"L39 gap4 first-gap D {first['start']['D_half']:.8g} to {first['end']['D_half']:.8g}; share {first['start']['share_half']:.8g} to {first['end']['share_half']:.8g}; all raw retention ratios "+str([r['end']['share_half']/r['start']['share_half'] for r in g])+('; first two gap D<0 fail the secondary S3 D>=0 criterion but retain >84% raw share' if m=='holmes_model3' else '')),
 dict(model=m,signature='commits after in-gap collapse',lower=None,upper=None,value=e['final_d'],qualification=f"L{e['displayed_length']} gap4 gap-end D "+str([r['end']['D_half'] for r in g])+('; deBelly first two collapse and final replacement agree with quarter-step final D '+str(debelly_comparison['companion_final_D'])+'; later-gap timing, third-gap bias sign and endpoint fields remain sensitive' if m=='debelly' else ''))])
classes=dict(goryachev='Initiation-limited',wave_pinning='Mass-conserved, saturating',holmes_model3='Mass-conserved, saturating',debelly='Slow non-activity state',spring='Slow non-activity state')
for row in evidence:
 row['class']=classes[row['model']]
 m=row['model'];example=next(r for r in examples if r['model']==m);gaps=[e for e in example['endpoints'] if e['kind']=='gap']
 row.update(evidence_length=39,reference_length=19,reference_value=None,native_time_upper=None,native_reference_time_upper=None,pulse_share=None,gap_end_share=None,gap1_D=None,gap2_D=None,gap3_D=None,post_recovery_D=None)
 if row['signature']=='initiation slows with L':
  endpoint=next(r for r in init if r['model']==m and r['length']==39);reference=next(r for r in init if r['model']==m and r['length']==19)
  row.update(reference_value=1.,native_time_upper=endpoint['time_upper'],native_reference_time_upper=reference['time_upper'])
  row['decision']='yes, size-dependent initiation' if row['model']=='goryachev' else 'unresolved within native first-update intervals' if row['model']=='debelly' else 'near-unity L39/L19, no substantial initiation slowing'
 elif row['signature']=='narrow-cue penalty and growth':
  reference=next(r for r in ratios if r['model']==m and r['length']==19)
  row['reference_value']=reference['ratio']['estimate_from_interpolation']
  if m=='holmes_model3':
   moderate=next(r for r in ratios if r['model']==m and r['length']==23)
   row.update(evidence_length=23,lower=moderate['ratio']['lower'],upper=moderate['ratio']['upper'],value=moderate['ratio']['estimate_from_interpolation'])
   row['qualification']='Deciding comparison is L23 versus L19; extreme L39 remains uninformative. '+row['qualification']
  row['decision']='large growing penalty over the resolved length range; extreme sizes retain bounds' if row['model'] in ('wave_pinning','holmes_model3') else 'modest non-growing penalty' if row['model'] in ('goryachev','spring') else 'no narrow-cue penalty'
 elif row['signature']=='share held through gaps':
  row.update(reference_length=39,reference_value=gaps[0]['start']['share_half'],pulse_share=gaps[0]['start']['share_half'],gap_end_share=gaps[0]['end']['share_half'],gap1_D=gaps[0]['end']['D_half'],gap2_D=gaps[1]['end']['D_half'],gap3_D=gaps[2]['end']['D_half'])
  row['decision']='retained raw share with incremental gains across pulse-gap cycles' if row['model'] in ('wave_pinning','holmes_model3') else 'resets' if row['model']=='goryachev' else 'early share collapses'
 else:
  row.update(reference_length=39,reference_value=gaps[0]['end']['D_half'],gap1_D=gaps[0]['end']['D_half'],gap2_D=gaps[1]['end']['D_half'],gap3_D=gaps[2]['end']['D_half'],post_recovery_D=example['final_d'])
  row['decision']='yes, replacement despite early collapse' if row['model'] in ('debelly','spring') else 'no, final D below replacement' if row['model']=='goryachev' else 'not applicable: share does not collapse'
 row['membership_changed']=False
csvout('class_membership_check.csv',evidence)
# End-gap noise counts use independent seeds, separately for each gap; no gap pseudoreplication.
noise_counts=[]
for kind in (() if CORE_PREVIEW else ('noise','mean','zero')):
 selected=[r for r in new if r['case']['kind']==kind]
 assert len(selected)==(12 if kind=='noise' else 1)
 for index in range(3):
  gaps=[[e for e in r['endpoints'] if e['kind']=='gap'][index] for r in selected]
  eligible=[e for e in gaps if e['primary_retention_denominator_eligible']]
  success=sum(bool(e['primary_retention_success']) for e in eligible)
  low,high=exact_binomial_interval(success,len(eligible)) if kind=='noise' else (None,None)
  noise_counts.append(dict(kind=kind,gap_index=index+1,seeds=len(selected),eligible=len(eligible),successes=success,fraction=success/len(eligible) if eligible else None,CI_lower=low,CI_upper=high,minimum_raw_share_retention=min(e['end']['share_half']/e['start']['share_half'] for e in gaps),maximum_raw_share_retention=max(e['end']['share_half']/e['start']['share_half'] for e in gaps)))
if noise_counts:csvout('holmes_L39_bounded_noise_retention.csv',noise_counts)
if not CORE_PREVIEW:
 noise_fates=[]
 for r in new:
  if r['case']['kind']=='train':continue
  recovery=next(e for e in r['endpoints'] if e['kind']=='recovery')['end']
  noise_fates.append(dict(task=r['task'],kind=r['case']['kind'],seed=r['case']['seed'],final_D=r['trial']['final_d'],final_replacement=recovery['replacement'],tail_pass=r['trial']['field_tail_pass']))
 csvout('holmes_L39_noise_fates.csv',noise_fates)
 noisy=[r for r in noise_counts if r['kind']=='noise'];nf=[r for r in noise_fates if r['kind']=='noise']
 noise_summary=dict(S3_successes=[r['successes'] for r in noisy],eligible=[r['eligible'] for r in noisy],minimum_raw_share_retention=min(r['minimum_raw_share_retention'] for r in noisy),maximum_raw_share_retention=max(r['maximum_raw_share_retention'] for r in noisy),minimum_final_D=min(r['final_D'] for r in nf),maximum_final_D=max(r['final_D'] for r in nf),final_replacements=sum(r['final_replacement'] for r in nf),independent_seeds=12)
plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'svg.fonttype':'none','pdf.fonttype':42})
fig=plt.figure(figsize=(12.5,14.5),layout='constrained');grid=fig.add_gridspec(6,1,height_ratios=(1.2,1.35,1.6,1.65,1.6,1.45))
fig.suptitle('Replacement resistance in enlarged cells arises through three mechanism classes;\ncue shape and timing distinguish their signatures',fontsize=13,fontweight='bold')
# A
sg=grid[0].subgridspec(1,3)
for j,(metric,title) in enumerate([('excess_amount','Incumbent amount'),('width','Incumbent extent'),('peak','Incumbent peak')]):
 ax=fig.add_subplot(sg[j])
 for m in ORDER:
  rows=sorted([r for r in org if r['model']==m],key=lambda r:float(r['size_native']))
  ax.plot([float(r['size_over_reference']) for r in rows],[float(r[metric+'_over_reference']) for r in rows],color=COLORS[m],label=NAMES[m])
 ax.set(title=('A  ' if j==0 else '')+title,xlabel='Length / own reference',ylabel='Fold of own reference');clean(ax)
 if j==2:ax.legend(frameon=False,fontsize=7,ncol=2)
# B
sg=grid[1].subgridspec(1,5)
for j,m in enumerate(ORDER):
 ax=fig.add_subplot(sg[j]);rows=[r for r in init if r['model']==m]
 if m=='debelly':
  ax.set_axis_off()
  inset=ax.inset_axes([.22,.37,.72,.55])
  inset.plot([19,39],[.5,.5],marker='v',mfc='white',mec=COLORS[m],ls='none',ms=5)
  for length in (19,39):inset.annotate('',xy=(length,.10),xytext=(length,.48),arrowprops={'arrowstyle':'->','color':COLORS[m]})
  inset.set(xlim=(15,43),ylim=(0,.65),xticks=[19,39],yticks=[0,.5],xlabel='Length L',ylabel='Initiation time')
  inset.set_title('Native-time upper bounds',fontsize=7.5);inset.tick_params(labelsize=7.5);inset.xaxis.label.set_size(7.5);inset.yaxis.label.set_size(7.5);clean(inset)
  ax.text(.5,.02,'Normalized ratio unavailable',ha='center',transform=ax.transAxes,fontsize=7.5)
 else:
  x=[];y=[]
  for r in rows:
   lo,hi=r['ratio_lower'],r['ratio_upper'];v=r['ratio_estimate'];assert lo<=v<=hi;ax.errorbar(r['length'],v,yerr=[[v-lo],[hi-v]],fmt='o',color=COLORS[m],ms=3,capsize=2);x.append(r['length']);y.append(v)
  ax.plot(x,y,color=COLORS[m])
 ax.set(title=('B  ' if j==0 else '')+NAMES[m],xlabel='Length L',yscale='log',ylim=(.75,4),xticks=[19,29,39]);clean(ax)
 if m!='debelly':ax.axhline(1,color='.6',ls=':',lw=.8)
 if j==0:ax.set_ylabel('Initiation time / own L19')
# C
ax=fig.add_subplot(grid[2])
for m in ORDER:
 x=[];y=[]
 for r in [r for r in ratios if r['model']==m]:
  lo,hi=r['ratio'].get('lower'),r['ratio'].get('upper')
  if not r['ratio'].get('informative',True) or lo is None:
   ax.text(r['length'],.68,'?',color=COLORS[m],ha='center');continue
  if hi is None:
   ax.plot(r['length'],lo,'^',mfc='white',mec=COLORS[m]);ax.annotate('',xy=(r['length'],lo*1.45),xytext=(r['length'],lo),arrowprops={'arrowstyle':'->','color':COLORS[m]})
  elif lo>0:
   v=r['ratio']['estimate_from_interpolation'];assert lo<=v<=hi;ax.errorbar(r['length'],v,yerr=[[v-lo],[hi-v]],fmt='o',ms=4,capsize=2,color=COLORS[m]);x.append(r['length']);y.append(v)
 if x:ax.plot(x,y,color=COLORS[m],label=NAMES[m])
ax.axhline(1,color='.6',ls=':');ax.set(title='C  Narrow cues impose the strongest equal-area penalty in mass-conserved models',xlabel='Length L',ylabel='Narrow / broad replacement time',yscale='log',ylim=(.55,200));clean(ax);ax.legend(frameon=False,ncol=5,fontsize=7.5)
# D
sg=grid[3].subgridspec(1,6,width_ratios=[1,1,1,1,1,.08])
for j,m in enumerate(ORDER):
 ax=fig.add_subplot(sg[j]);values=np.full((6,6),np.nan)
 for r in fates:
  if r['model']!=m:continue
  xi=LENGTHS.index(r['length']);yi=GAPS.index(r['gap'])
  if r['eligible']:values[yi,xi]=r['value']
  else:ax.text(xi,yi,'×',ha='center',va='center')
  if r['status']=='observed_retention_without_target_calibration':ax.plot(xi,yi,'D',mfc='none',mec='black',ms=6,mew=.8)
  if r['tail_pass'] is False:ax.plot(xi,yi,'x',color='black',ms=3,mew=.65)
 im=ax.imshow(values,origin='lower',vmin=-1,vmax=1,cmap='coolwarm',aspect='auto');ax.set(xticks=[0,2,5],xticklabels=[19,27,39],yticks=range(6),yticklabels=GAPS,title=('D  ' if j==0 else '')+NAMES[m],xlabel='Length L')
 ax.title.set_color(COLORS[m])
 if j==0:ax.set_ylabel('Gap / pulse')
fig.colorbar(im,cax=fig.add_subplot(sg[5]),label='D (◇ gap end; otherwise fate)')
# E
sg=grid[4].subgridspec(1,5);tracecsv=[]
for j,m in enumerate(ORDER):
 ax=fig.add_subplot(sg[j]);e=next(r for r in examples if r['model']==m)
 source=e['source'];folder=Path(source['run_root'])/'workers'/f"{source['task']:04d}"/'trial';measurement=read(folder/'readouts.json');pulse=e['schedule']['pulse_duration'];stop=e['schedule']['segments'][-1][2]
 time=np.array(measurement['series']['time']);d=np.array([np.nan if v is None else v for v in measurement['series']['D_half']]);keep=np.ones(len(time),bool) if m=='debelly' else time<=stop+pulse
 select=display_indices(time[keep],d[keep],endpoints=[ep['end']['time'] for ep in e['endpoints'] if ep['kind'] in ('on','gap')],bins=350)
 ax.plot(time[keep][select]/pulse,d[keep][select],color=COLORS[m],lw=1.2,ls=':')
 tracecsv.extend(dict(model=m,time=float(t),pulse_time=float(t/pulse),D=None if not np.isfinite(v) else float(v),share=None if not np.isfinite(v) else float((v+1)/2)) for t,v in zip(time[keep][select],d[keep][select]))
 for kind,start,end in e['schedule']['segments']:
  if kind=='on':ax.axvspan(start/pulse,end/pulse,color=COLORS[m],alpha=.15)
 for ep in e['endpoints']:
  if ep['kind'] in ('on','gap') and ep['end']['D_half'] is not None:ax.plot(ep['end']['time']/pulse,ep['end']['D_half'],'o' if ep['kind']=='on' else 's',color=COLORS[m],ms=3)
 if m=='debelly':ax.plot(time[-1]/pulse,d[-1],'o',color=COLORS[m],ms=4)
 gap=next(ep for ep in e['endpoints'] if ep['kind']=='gap')
 retention=gap['end']['share_half']/gap['start']['share_half']
 ax.text(.98,.05,f'First-gap share retained: {retention:.1%}',ha='right',transform=ax.transAxes,fontsize=7.5)
 ax.text(.02,.94,f'Post-recovery D = {e["final_d"]:.3f}',ha='left',transform=ax.transAxes,fontsize=7.5)
 ax.axhline(.5,color='.6',ls=':',lw=.7);ax.set(title=('E  ' if j==0 else '')+NAMES[m]+' L39'+('*' if m=='debelly' else ''),xlabel='Time / pulse',ylim=(-1.05,1.05));clean(ax)
 if j==0:ax.set_ylabel('Activity dominance D')
csvout('gap_retention_traces.csv',tracecsv)
# F: every signature has a numeric row in class_membership_check.csv.
ax=fig.add_subplot(grid[5]);ax.axis('off');ax.set_title('F  Shape and timing signatures distinguish mechanism classes',loc='left')
measurement='Beam width at fixed photon dose;\npulse trains; gap-end biosensor share'
classrows=[['Initiation-limited\nGoryachev','Yes (2.90×)','Modest (1.12×)','Resets','No',measurement],['Mass-conserved, saturating\nWave-pinning, Holmes','No resolved slowing','Large; grows over\nresolved L range','Share retained\n(ratchet)','—',measurement],['Slow non-activity state\ndeBelly*, Spring','Unresolved* / no','Modest (~1×)','Collapses early','Yes',measurement]]
table=ax.table(cellText=classrows,colLabels=['Class / models','Initiation slows','Equal-area penalty','Share through gaps','Commit after collapse','Experimental measurement'],cellLoc='center',bbox=[0,.23,1,.70],colWidths=[.22,.12,.17,.16,.12,.21]);table.auto_set_font_size(False);table.set_fontsize(7.5)
ax.text(0,.12,'◇ Holmes L31–L39: last inter-pulse gap-end D, with no target-length calibration. × failed field tail; ? uninformative ratio.',transform=ax.transAxes,fontsize=7.2)
ax.text(0,.03,'* deBelly initiation is interval-censored; later-gap timing and endpoint fields carry numerical qualifications (complete supplement).',transform=ax.transAxes,fontsize=7.3)
stem='figure2_round2_preview' if CORE_PREVIEW else 'figure2_round2'
for ext in ('png','pdf','svg'):fig.savefig(OUT/f'{stem}.{ext}',dpi=180,facecolor='white')
plt.close(fig)
if CORE_PREVIEW:
 print(json.dumps(dict(status='core_only_preview',audited_Holmes_workers=18,main=str(OUT/f'{stem}.pdf'),noise_ensemble_complete=False)))
 sys.exit(0)
# Supplement additions, followed by all preceding supplement pages and main plate.
with PdfPages(OUT/'supporting/figure2_round2_additions.pdf') as pdf:
 fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title('Holmes large-L retention: completed frozen design and qualifications',loc='left',fontsize=15)
 text='No continuous commitment by 100,000 is an observed finite-horizon finding.\n\nL27 upper committing duration = 18,125; four pulses each 6,796.875.\nTarget lengths: 31, 35, 39; gaps: 0, 0.1, 0.5, 1, 4, 16 × pulse.\nPrimary: D and challenger share at every pulse and gap endpoint. Fate is secondary.\n\nL39 bounded-noise check: gap16, mean 2% native peak, maximum 4%;\n12 independent seeds plus matched-mean and zero controls.\nS3-matched nonperiodic spatial OU input: ell 1.9, tau 10, knot 2.5.\n\nAll 32 native SLURM trajectories and independent numeric/phase/restart/tail\naudits completed. No kinetics or calibration thresholds were retuned.\nLarge-L Holmes D can fall below zero while >84% of raw share is retained.\nThe secondary S3 criterion requires both D ≥ 0 and raw retention ≥80%.\n\nClass membership is unchanged. deBelly initiation remains interval-unresolved;\nHolmes width-penalty growth is established only over measured moderate sizes.'
 text=text.replace('Large-L Holmes D can fall below zero while >84% of raw share is retained.','Clean L39 gap4 retains >84% raw share even with slightly negative gap-end D.')
 text+='\n\nThe theoretical mean gap-integrated commanded input is 2.46 × ON exposure;\nthis long-gap assay is weak instantaneously, not a matched-total-dose test.\nHeterogeneous input yields variable gap retention and 10/12 final replacements;\nthe matched-mean and zero controls remain below replacement.'
 ax.text(.02,.92,text,va='top',linespacing=1.5,fontsize=10.5);pdf.savefig(fig);plt.close(fig)
 # All 18 targeted deterministic histories, including every exact endpoint.
 for length in (31,35,39):
  fig,axes=plt.subplots(2,3,figsize=(11.7,8.3),layout='constrained',sharey=True)
  fig.suptitle(f'Holmes L{length}: retained challenger share across all six gap conditions',fontsize=15)
  for ax,gap in zip(axes.flat,GAPS):
   row=next(r for r in new if r['case']['kind']=='train' and r['case']['length']==length and r['case']['gap']==gap)
   pulse=row['trial']['schedule']['pulse_duration'];series=read(Path(row['result_path']).parent/'trial/readouts.json')['series'];time=np.asarray(series['time']);D=np.asarray(series['D_half'],float)
   index=display_indices(time,D,[e['end']['time'] for e in row['endpoints']],bins=300)
   ax.plot(time[index]/pulse,D[index],color=COLORS['holmes_model3'],lw=.9)
   for kind,a,b in row['trial']['schedule']['segments']:
    if kind=='on':ax.axvspan(a/pulse,b/pulse,color=COLORS['holmes_model3'],alpha=.12)
   for e in row['endpoints']:
    if e['kind'] in ('on','gap'):ax.plot(e['end']['time']/pulse,e['end']['D_half'],'o' if e['kind']=='on' else 's',color=COLORS['holmes_model3'],ms=3)
   ax.axhline(0,color='.65',ls=':',lw=.7);ax.axhline(.5,color='.65',ls='--',lw=.7)
   ax.set(title=f"Gap / pulse = {gap:g}; final D = {row['trial']['final_d']:.4f}",xlabel='Time / pulse',ylabel='Activity dominance D',ylim=(-1.05,1.05));clean(ax)
  pdf.savefig(fig);plt.close(fig)
 # Independent noise histories; each gap is evaluated separately, with exact binomial CIs.
 fig,axes=plt.subplots(2,2,figsize=(11.7,8.3),layout='constrained');axes=axes.ravel()
 fig.suptitle('Holmes L39 gap16: bounded-input retention and separate dominance criterion',fontsize=14)
 for row in new:
  kind=row['case']['kind']
  if kind=='train':continue
  gaprows=[e for e in row['endpoints'] if e['kind']=='gap'];style=dict(color=COLORS['holmes_model3'] if kind=='noise' else '.3' if kind=='mean' else '.65',alpha=.45 if kind=='noise' else 1,lw=.7 if kind=='noise' else 1.5,ls='-' if kind=='noise' else '--' if kind=='mean' else ':')
  label='12 independent bounded inputs' if row['case']['task']==18 else 'Matched-mean control' if kind=='mean' else 'Zero control' if kind=='zero' else None
  axes[0].plot([1,2,3],[e['end']['D_half'] for e in gaprows],marker='o',ms=3,label=label,**style)
  axes[1].plot([1,2,3],[e['end']['share_half']/e['start']['share_half'] for e in gaprows],marker='o',ms=3,**style)
 for j,kind in enumerate(('noise','mean','zero')):
  group=[r for r in noise_counts if r['kind']==kind]
  for r in group:
   x=r['gap_index']+.08*(j-1)
   if r['fraction'] is not None:
    axes[2].plot(x,r['fraction'],'o' if kind=='noise' else 's' if kind=='mean' else 'x',color=COLORS['holmes_model3'] if kind=='noise' else '.3' if kind=='mean' else '.65',ms=4)
    if r['CI_lower'] is not None:axes[2].vlines(x,r['CI_lower'],r['CI_upper'],color=COLORS['holmes_model3'])
 axes[0].axhline(0,color='.5',ls=':');axes[1].axhline(.8,color='.5',ls=':')
 for ax,title,ylabel in zip(axes,('End-gap dominance','Raw challenger-share retention','Secondary S3 criterion'),('D at gap end','Share at gap end / pulse end','Conditional success fraction')):
  ax.set(title=title,xlabel='Gap index',ylabel=ylabel,xticks=[1,2,3]);clean(ax)
 axes[0].set_ylim(-1.05,1.05);axes[2].set_ylim(-.05,1.05);axes[0].legend(frameon=False,fontsize=7,loc='lower right')
 for index,row in enumerate(noise_fates):
  axes[3].plot(index,row['final_D'],'o' if row['kind']=='noise' else 's',color=COLORS['holmes_model3'] if row['kind']=='noise' else '.3' if row['kind']=='mean' else '.65',mfc=COLORS['holmes_model3'] if row['kind']=='noise' else 'white')
 axes[3].axhline(.5,color='.5',ls=':');axes[3].set(title='Secondary final fate: 10/12 noisy histories replace',xlabel='Independent history / control',ylabel='Post-recovery D',ylim=(-1.05,1.05),xticks=[0,3,6,9,12,13],xticklabels=['1','4','7','10','Mean','Zero']);clean(axes[3])
 pdf.savefig(fig);plt.close(fig)
 # Compact endpoint values for the complete deterministic grid, four pulses and three gaps.
 for length in (31,35,39):
  fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title(f'Holmes L{length}: all pulse and gap endpoint D values',loc='left',fontsize=15)
  table_rows=[]
  for gap in GAPS:
   row=next(r for r in new if r['case']['kind']=='train' and r['case']['length']==length and r['case']['gap']==gap)
   endpoint=[e for e in row['endpoints'] if e['kind'] in ('on','gap')]
   table_rows.append([f'{gap:g}']+[f"{e['end']['D_half']:.7f}" for e in endpoint]+[f"{row['trial']['final_d']:.7f}"])
  table=ax.table(cellText=table_rows,colLabels=['Gap / pulse','Pulse 1','Gap 1','Pulse 2','Gap 2','Pulse 3','Gap 3','Pulse 4','Recovery'],cellLoc='center',bbox=[0,.28,1,.55]);table.auto_set_font_size(False);table.set_fontsize(8)
  ax.text(0,.13,'Gap0 has three exact zero-duration inter-pulse boundaries. Full raw-share values, eligibility,\nS3 success, field tails and final fate are in holmes_large_L_endpoints.csv. Values are observed\nendpoints; the target-length commitment calibration is unavailable.',fontsize=10)
  pdf.savefig(fig);plt.close(fig)
 for m in ORDER:
  fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title('Class-membership evidence: '+NAMES[m],loc='left',fontsize=15)
  import textwrap
  text='\n\n'.join(r['signature']+f"\nValue: {r['value']}; bounds: [{r['lower']}, {r['upper']}]\n"+textwrap.fill(r['qualification'],100) for r in evidence if r['model']==m)
  ax.text(.02,.92,text,va='top',fontsize=9,linespacing=1.5);pdf.savefig(fig);plt.close(fig)
writer=PdfWriter()
for p in (OUT/'supporting/figure2_round2_additions.pdf',OLD/'supporting/vectors/figure2_signature_assays.pdf',OLD/'supporting/vectors/figure2_signature_assays_complete_supplement.pdf'):
 writer.append(str(p))
writer.write(str(OUT/'figure2_round2_supplement.pdf'))
combined_pages=PdfReader(OUT/'figure2_round2_supplement.pdf').pages
previous_pages=list(PdfReader(OLD/'supporting/vectors/figure2_signature_assays.pdf').pages)+list(PdfReader(OLD/'supporting/vectors/figure2_signature_assays_complete_supplement.pdf').pages)
assert all(a.get_contents().get_data()==b.get_contents().get_data() for a,b in zip(combined_pages[-len(previous_pages):],previous_pages))
legend=[
'Enlarged cells resist replacement through distinct initiation, cue-width and gap-retention signatures, supporting three mechanism classes in the deterministic cue assays.',
'A, Incumbent amount, extent and peak across the retained 97 states. B, Initiation normalized to each model\u2019s L19 value at central area and the fixed 10% L19 reporter threshold. Brackets propagate saved-sample intervals. deBelly L19 and L39 initiation both lie in (0, 0.5] native units; its normalized ratio is unavailable. Its separate native-time inset places open downward triangles at the measured upper bounds 0.5 and at the actual lengths 19 and 39, with arrows toward smaller times.',
'C, Narrow/broad replacement-time ratio at fixed commanded spatial area. Open triangles and arrows are one-sided bounds; question marks denote uninformative ratios. Wave-pinning has a growing large penalty. Holmes grows from about 29 at L19 to above 133 at L23; growth through L39 is unresolved because both width arms are censored. Goryachev and Spring have modest penalties rather than exact width invariance. Commanded area is not calibrated photon dose.',
'D, Four equal pulses with total ON 1.5 times the model-length commitment calibration; final D after recovery is the colour. For Holmes L31, L35 and L39, diamonds identify the last inter-pulse gap-end D from the targeted L27-referenced protocol; this endpoint is distinct from final fate. Gap0 uses the zero-duration boundary after the third pulse. No continuous commitment by 100,000 at these three lengths is extreme finite-horizon replacement resistance. Crosses preserve failed field tails, including '+str(sum(r['model']=='debelly' and r['tail_pass'] is False for r in fates))+' of the 36 deBelly fate cells.',
'E, L39 gap4 histories. Shading marks cue ON; circles and squares mark pulse and gap endpoints. The deBelly trace includes full recovery to its measured final fate; other traces show through one pulse duration of recovery, with full post-recovery D stated in each plot. Holmes retains more than 84% of raw challenger share through all three gaps, although D is slightly negative in the first two and final D remains below 0.5. Retention of share differs from the secondary S3 criterion requiring D at least zero. F, The class table has a numerical evidence row for every model and signature; membership is unchanged. deBelly early collapse and final replacement survive its separate quarter-step check, while later-gap timing and endpoint fields remain sensitive. All native kinetics are unchanged. These polarity readouts are distinct from observed migration and conditional OU persistence.'
]
legend.append('The separate Holmes L39 gap16 transfer has heterogeneous bounded input during gaps: S3 successes are '+', '.join(f"{r['successes']}/{r['eligible']}" for r in noisy)+' after gaps 1–3. Raw share ratios span '+f"{noise_summary['minimum_raw_share_retention']:.5f}–{noise_summary['maximum_raw_share_retention']:.5f}"+'; final D spans '+f"{noise_summary['minimum_final_D']:.6f}–{noise_summary['maximum_final_D']:.6f}"+'. Ten of twelve noisy histories replace, while matched-mean and zero controls do not. The class table describes gaps without added input; this separate transfer demonstrates variable retention and fate under spatially heterogeneous input.')
limitations=['deBelly normalized initiation unresolved within first-update intervals','Holmes width-penalty growth through L39 uninformative; moderate-length growth verified','deBelly later-gap timing and endpoint fields remain numerically sensitive','Holmes share retention is distinct from S3 D>=0 retention success and from commitment','Modest Goryachev/Spring width penalties qualify binary no-penalty wording; no class change','Selected mechanistic signatures do not causally allocate replacement resistance']
limitations.append('Bounded heterogeneous Holmes gap input gives S3 successes 8/12, 9/12, 10/12; some histories lose share or replace, unlike deterministic/mean controls. Integrated gap exposure is 2.46 times ON exposure; deterministic class membership does not imply universal noisy-gap retention.')
write(OUT/'legend.json',dict(id='figure2',title='Replacement resistance in enlarged cells arises through three mechanism classes; cue shape and timing distinguish their signatures',status='qualified_complete',assets=[dict(path=str((OUT/'figure2_round2.png').relative_to(ROOT)),caption='Figure 2 round2 with completed Holmes retention',page_size=[12.5,14.5],minimum_printed_width_inches=10,page_title='Figure 2. Shape and timing distinguish mechanism classes')],legend=legend,limitations=limitations,sources=[dict(label='Class evidence CSV',path=str((OUT/'class_membership_check.csv').relative_to(ROOT))),dict(label='Holmes endpoints',path=str((OUT/'holmes_large_L_endpoints.csv').relative_to(ROOT)))]))
csvout('design_parameters.csv',[dict(incumbent_states=len(org),pulse_count=4,total_on_factor=1.5,Holmes_reference_calibration=18125,Holmes_pulse_duration=6796.875,Holmes_continuous_horizon=100000,Holmes_deterministic_workers=18,Holmes_noise_seeds=12,Holmes_controls=2,Holmes_noise_mean=.02,Holmes_noise_max=.04,Holmes_noise_ell=1.9,Holmes_noise_tau=10,Holmes_noise_knot=2.5,Holmes_noise_gap_total_duration=326250.,Holmes_noise_expected_gap_integrated_input=254475.,Holmes_on_integrated_input=103312.5,Holmes_noise_expected_gap_to_on_input_ratio=254475/103312.5,debelly_failed_fate_tails=sum(r['model']=='debelly' and r['tail_pass'] is False for r in fates),debelly_fate_cells=36)])
generated=datetime.now(timezone.utc).isoformat()
pages=len(PdfReader(OUT/'figure2_round2_supplement.pdf').pages)
readme=f'''# Figure 2 round 2 — replacement-resistance signatures

Enlarged cells resist replacement through distinct initiation, cue-width and gap-retention signatures, supporting three mechanism classes in deterministic cue assays. Class membership is unchanged.

- Main: [PNG](figure2_round2.png), [PDF](figure2_round2.pdf), [SVG](figure2_round2.svg); [legend and integration assets](legend.json).
- [{pages}-page complete supplement](figure2_round2_supplement.pdf): new Holmes train/noise and numerical class-check pages, followed by all 83 preceding main/supplement pages. Their decoded page contents are verified unchanged, preserving the per-model table, profile inset, reporter sensitivity and numerical qualifications.
- Plotted values: [organization](incumbent_organization.csv), [initiation](initiation_normalized.csv), [width ratios](equal_area_width_ratios.csv), [fate/retention maps](train_fates_and_retention.csv), [endpoints](all_selected_train_endpoints.csv), [traces](gap_retention_traces.csv).
- [Every class cell](class_membership_check.csv); [all Holmes endpoints](holmes_large_L_endpoints.csv); [noise counts/CIs](holmes_L39_bounded_noise_retention.csv); [individual noisy/control fates](holmes_L39_noise_fates.csv); [design](design_parameters.csv); [deBelly numerical comparison](debelly_numerical_qualification.csv); [verification](supporting/metadata/verification.json); [brief acceptance](supporting/metadata/brief_acceptance.json).

Compared with 20 October, the main has six message-led panels: compressed organization, common log initiation axes normalized to own L19, one width-ratio plot, class-ordered fate maps and L39 traces, and a class table. All targeted Holmes measurements are complete. Diamonds distinguish their gap-end D from post-recovery fate. Crosses restore failed field tails, including 31/36 deBelly cells. The deBelly native-time inset contains actual upper bounds, and its trace now includes full recovery and late replacement. Pending/preview plates remain superseded history.

The numerical class check supports Goryachev initiation slowing (~2.90-fold L39/L19), large width penalties and share retention in Wave-pinning/Holmes, and early share collapse followed by replacement in deBelly/Spring. Three qualifications remain essential. deBelly initiation at both endpoint sizes lies in (0, 0.5], leaving its normalized ratio unresolved. Holmes width growth is established from L19 (~29.03) to L23 (>133.16), while both L39 arms are censored. Goryachev/Spring have modest, non-growing width penalties (~1.12/~1.16 at L39), rather than exact invariance. deBelly's first two gap losses and final replacement survive its quarter-step companion (final D 0.916459 versus 0.920853), but third-gap bias sign, later timing and endpoint fields remain sensitive. These signatures do not causally allocate resistance among molecular terms.

Holmes has no continuous commitment by 100,000 at L31/L35/L39. Its 18 targeted trains use the fixed L27 calibration: four pulses of 6,796.875 and the six specified gap multiples. At L39 gap4, all three gaps retain >84% of raw challenger share, while the first two gap-end D values are slightly negative and final D remains below replacement. Raw-share retention is therefore distinct from S3 success (D≥0 plus ≥80% retention), and from commitment. At gap0 the displayed diamond is the exact boundary after pulse3.

The S3-matched L39 gap16 check adds 12 independent bounded-input histories and matched-mean/zero controls; native kinetics are unchanged. All actual queried inputs satisfy the 4% native-peak bound. Each gap has a separate exact seed-binomial CI. The theoretical mean gap-integrated commanded input is 2.46 times the ON exposure: this is an instantaneously weak-input check, not a matched-total-dose test or empirical noise calibration. S3 successes are 8/12, 9/12 and 10/12 after gaps 1–3; individual raw-share retention spans about 0.02152–1.78070, and final D spans −0.063410–0.999944. Ten of twelve noisy histories replace, while both controls remain below replacement. The deterministic class table does not claim universal noisy-gap retention. Full settings and counts are linked above.

All 32 outputs pass independent frozen-source readout, phase, exact-basal restart, artifact-hash and field-tail reconstruction. Native array: 21279666; frozen final audit: 21279960. The existing scaling_sandbox environment and red/QoS small were used without a concurrency throttle. Palette: WP #2878b5, Goryachev #c83b3b, deBelly #82599c, Holmes #36944c, Spring #c59a25. L39 traces are dotted; bounds are open triangles with arrows, and uninformative ratios remain explicit. Commanded spatial area is not calibrated photon dose. WGD-associated enlargement motivates length changes; polarity readouts remain distinct from observed migration and conditional OU persistence.

Generated UTC: {generated}. Redraw: `python -B polarity1d/manuscript/code/mechanisms/revision_round2_figure2_20261022.py`. Native/audit source: `polarity1d/manuscript/data/runs/signature_holmes_retention_round2_20261022/`. Redraw launches no integration or timed refresh.
'''
(OUT/'README.md').write_text(readme)
write(OUT/'supporting/metadata/verification.json',dict(generated_utc=generated,status='qualified_complete',producer_state='stopped',main_panels=list('ABCDEF'),large_size_line_style='dotted',debelly_native_time_bounds_separate=True,readme_sha256=sha(OUT/'README.md'),reporter_source_audit_verified=True,organization_source_digest_verified=True,source_train_count=len(trains),source_reporter_count=len(reporters),new_Holmes_workers=32,noise_transfer_summary=noise_summary,independent_native_audit_sha256=sha(RUN/'results/independent_native_readout_audit.json'),class_evidence_cells=len(evidence),fate_or_retention_cells=len(fates),debelly_failed_field_tails_visible=31,complete_supplement_pages=pages,bound_convention='open triangle with arrow; uninformative values explicitly marked',palette=COLORS,class_membership_changed=False,complete=True,qualifications=limitations,sha256={str(p.relative_to(OUT)):sha(p) for p in OUT.glob('*') if p.suffix in ('.csv','.pdf','.svg','.png','.json')}))
write(OUT/'supporting/metadata/brief_acceptance.json',dict(generated_utc=generated,complete=True,requirements={
 '3.1_Holmes_18_targeted_trains':{'passed':True,'evidence':'holmes_large_L_endpoints.csv'},
 '3.1_Holmes_bounded_noise_12_seeds_and_controls':{'passed':True,'evidence':'holmes_L39_bounded_noise_retention.csv'},
 '3.1_variable_noise_retention_and_fates':{'passed':True,'evidence':'holmes_L39_bounded_noise_retention.csv; holmes_L39_noise_fates.csv','summary':noise_summary},
 '3.1_primary_exact_D_and_share_endpoints':{'passed':True,'evidence':'holmes_large_L_endpoints.csv'},
 '3.1_continuous_commitment_horizon_finding':{'passed':True,'evidence':'design_parameters.csv'},
 '3.2_six_main_panels_A_to_F':{'passed':True,'evidence':'figure2_round2.svg'},
 '3.2_native_censoring_without_ratio_point_estimates':{'passed':True,'evidence':'initiation_normalized.csv; equal_area_width_ratios.csv'},
 '3.2_Holmes_large_L_distinct_gap_endpoint_diamonds':{'passed':True,'evidence':'train_fates_and_retention.csv'},
 '3.2_L39_examples_for_every_model':{'passed':True,'evidence':'gap_retention_traces.csv'},
 '3.2_all_20_class_signature_cells':{'passed':True,'evidence':'class_membership_check.csv','qualification':'Unresolved deBelly initiation and extreme-Holmes width growth stay explicit; no membership change'},
 '3.3_new_Holmes_train_noise_pages_and_class_checks':{'passed':True,'evidence':'supporting/figure2_round2_additions.pdf'},
 '3.3_prior_per_model_table_and_complete_supplement_preserved':{'passed':True,'evidence':'figure2_round2_supplement.pdf','original_pages_preserved':len(previous_pages),'decoded_page_content_identity_verified':True},
 'author_steering_deBelly_failed_field_tails':{'passed':True,'evidence':'train_fates_and_retention.csv; legend.json','count':31},
 '4_PNG_PDF_SVG_supplement_CSV_README_legend':{'passed':True,'evidence':'supporting/metadata/verification.json'}},
 model_palette=COLORS,class_membership_changed=False,no_retuning=True))
print(json.dumps(dict(main=str(OUT/'figure2_round2.pdf'),supplement_pages=pages,status='qualified_complete')))
