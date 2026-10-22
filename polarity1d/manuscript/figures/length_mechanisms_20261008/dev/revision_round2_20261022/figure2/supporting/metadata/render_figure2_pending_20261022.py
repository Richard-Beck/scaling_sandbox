"""Round2 Figure2 recomposition from verified existing values; Holmes gaps remain pending."""
from pathlib import Path
import csv, json, hashlib, sys
from datetime import datetime,timezone
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
OUT=Path(__file__).resolve().parent
OLD=OUT.parents[1]/'signature_assays_20261020'
sys.path.insert(0,str(ROOT/'polarity1d/manuscript/code/mechanisms'))
from signature_analysis_20261020 import replacement_time_ratio
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
 with (OUT/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def clean(ax):ax.spines[['top','right']].set_visible(False)
for d in ('supporting','supporting/data','supporting/metadata'):(OUT/d).mkdir(exist_ok=True)
coverage=read(OLD/'supporting/metadata/coverage.json')
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
# Ratios retain both measurement intervals; a zero lower denominator gives no finite upper bound.
init=[]
for model in ORDER:
 selected=sorted([r for r in reporters if r['model']==model and r['area']==1 and r['width']==3.8],key=lambda r:r['length'])
 ref=selected[0]['fixed_reporter_events']['0.1']
 for r in selected:
  event=r['fixed_reporter_events']['0.1'];ratio=replacement_time_ratio(event,ref)
  if r['length']==19:ratio.update(lower=1.,upper=1.,estimate_from_interpolation=1.,status='reference_identity')
  init.append(dict(model=model,length=r['length'],time_lower=event['lower_elapsed'],time_upper=event['upper_elapsed'],ratio_lower=ratio['lower'],ratio_upper=ratio['upper'],ratio_estimate=ratio.get('estimate_from_interpolation'),status='normalization_unresolved' if ref['lower_elapsed']==0 and r['length']!=19 else ratio['status']))
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
for l in (31,35,39):
 for g in GAPS:fates.append(dict(model='holmes_model3',length=l,gap=g,value=None,eligible=None,measurement='end-gap D requested',tail_pass=None,status='pending_SLURM_submission'))
csvout('train_fates_and_pending_retention.csv',fates);csvout('all_selected_train_endpoints.csv',endrows)
# Numeric evidence for every proposed table cell. No class assignment is changed.
evidence=[]
for m in ORDER:
 i=[r for r in init if r['model']==m and r['length']==39][0]
 w=[r for r in ratios if r['model']==m and r['length'] in (19,39)]
 e=next(r for r in examples if r['model']==m);g=[r for r in e['endpoints'] if r['kind']=='gap'];first=g[0]
 evidence.extend([
 dict(model=m,signature='initiation slows with L',lower=i['ratio_lower'],upper=i['ratio_upper'],value=i['ratio_estimate'],qualification='deBelly first-update intervals do not resolve a normalized size effect' if m=='debelly' else 'L39/L19, central area, fixed 10% L19 reporter'),
 dict(model=m,signature='narrow-cue penalty and growth',lower=w[-1]['ratio'].get('lower'),upper=w[-1]['ratio'].get('upper'),value=w[-1]['ratio'].get('estimate_from_interpolation'),qualification='L19 '+str(w[0]['ratio'])+'; L39 '+str(w[-1]['ratio'])),
 dict(model=m,signature='share held through gaps',lower=None,upper=None,value=first.get('primary_retention_ratio'),qualification=f"L{e['displayed_length']} gap4 first-gap D {first['start']['D_half']:.8g} to {first['end']['D_half']:.8g}; final D {e['final_d']:.8g}; Holmes L39 pending" if m=='holmes_model3' else f"L39 gap4 first-gap D {first['start']['D_half']:.8g} to {first['end']['D_half']:.8g}"),
 dict(model=m,signature='commits after in-gap collapse',lower=None,upper=None,value=e['final_d'],qualification=f"L{e['displayed_length']} gap4 gap-end D "+str([r['end']['D_half'] for r in g])+('; deBelly later-gap timing sensitive; first two collapse and final fate agree' if m=='debelly' else ''))])
csvout('class_membership_check.csv',evidence)
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
  ax.plot(19,1,'o',color=COLORS[m],ms=3)
  ax.text(.5,.65,'Both endpoint sizes:\n0 < initiation ≤ 0.5\n(native units)',ha='center',transform=ax.transAxes,fontsize=8)
  ax.plot([.3,.7],[.38,.38],marker='v',mfc='white',mec=COLORS[m],ls='none',transform=ax.transAxes)
  for xx in (.3,.7):ax.annotate('',xy=(xx,.25),xytext=(xx,.37),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':COLORS[m]})
  ax.text(.5,.05,'L19 normalization unresolved',ha='center',transform=ax.transAxes,fontsize=7)
 else:
  x=[];y=[]
  for r in rows:
   lo,hi=r['ratio_lower'],r['ratio_upper'];v=np.sqrt(lo*hi);ax.errorbar(r['length'],v,yerr=[[v-lo],[hi-v]],fmt='o',color=COLORS[m],ms=3,capsize=2);x.append(r['length']);y.append(v)
  ax.plot(x,y,color=COLORS[m])
 ax.set(title=('B  ' if j==0 else '')+NAMES[m],xlabel='Length L',yscale='log',ylim=(.75,4),xticks=[19,29,39]);ax.axhline(1,color='.6',ls=':',lw=.8);clean(ax)
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
   v=np.sqrt(lo*hi);ax.errorbar(r['length'],v,yerr=[[v-lo],[hi-v]],fmt='o',ms=4,capsize=2,color=COLORS[m]);x.append(r['length']);y.append(v)
 if x:ax.plot(x,y,color=COLORS[m],label=NAMES[m])
ax.axhline(1,color='.6',ls=':');ax.set(title='C  Narrow cues impose the strongest equal-area penalty in mass-conserved models',xlabel='Length L',ylabel='Narrow / broad replacement time',yscale='log',ylim=(.55,200));clean(ax);ax.legend(frameon=False,ncol=5,fontsize=7)
# D
sg=grid[3].subgridspec(1,6,width_ratios=[1,1,1,1,1,.08])
for j,m in enumerate(ORDER):
 ax=fig.add_subplot(sg[j]);values=np.full((6,6),np.nan)
 for r in fates:
  if r['model']!=m:continue
  xi=LENGTHS.index(r['length']);yi=GAPS.index(r['gap'])
  if r['status']=='observed' and r['eligible']:values[yi,xi]=r['value']
  elif r['status']!='observed':ax.plot(xi,yi,'D',mfc='white',mec=COLORS[m],ms=5)
  else:ax.text(xi,yi,'×',ha='center',va='center')
 im=ax.imshow(values,origin='lower',vmin=-1,vmax=1,cmap='coolwarm',aspect='auto');ax.set(xticks=[0,2,5],xticklabels=[19,27,39],yticks=range(6),yticklabels=GAPS,title=('D  ' if j==0 else '')+NAMES[m],xlabel='Length L')
 if j==0:ax.set_ylabel('Gap / pulse')
fig.colorbar(im,cax=fig.add_subplot(sg[5]),label='Post-recovery D')
# E
sg=grid[4].subgridspec(1,5);tracecsv=[]
for j,m in enumerate(ORDER):
 ax=fig.add_subplot(sg[j]);e=next(r for r in examples if r['model']==m)
 if m=='holmes_model3':ax.text(.5,.5,'L39 retention train\nprepared; SLURM blocked\n◇ endpoint D pending',ha='center',va='center',transform=ax.transAxes,color=COLORS[m])
 else:
  source=e['source'];folder=Path(source['run_root'])/'workers'/f"{source['task']:04d}"/'trial';measurement=read(folder/'readouts.json');pulse=e['schedule']['pulse_duration'];stop=e['schedule']['segments'][-1][2]
  time=np.array(measurement['series']['time']);d=np.array([np.nan if v is None else v for v in measurement['series']['D_half']]);keep=time<=stop+pulse
  ax.plot(time[keep]/pulse,d[keep],color=COLORS[m],lw=1)
  tracecsv.extend(dict(model=m,time=float(t),pulse_time=float(t/pulse),D=None if not np.isfinite(v) else float(v)) for t,v in zip(time[keep],d[keep]))
  for kind,start,end in e['schedule']['segments']:
   if kind=='on':ax.axvspan(start/pulse,end/pulse,color=COLORS[m],alpha=.15)
  for ep in e['endpoints']:
   if ep['kind'] in ('on','gap') and ep['end']['D_half'] is not None:ax.plot(ep['end']['time']/pulse,ep['end']['D_half'],'o' if ep['kind']=='on' else 's',color=COLORS[m],ms=3)
 ax.axhline(.5,color='.6',ls=':',lw=.7);ax.set(title=('E  ' if j==0 else '')+NAMES[m]+' L39'+('*' if m=='debelly' else ''),xlabel='Time / pulse',ylim=(-1.05,1.05));clean(ax)
 if j==0:ax.set_ylabel('Activity dominance D')
csvout('gap_retention_traces.csv',tracecsv)
# F: classes proposed by brief; explicit unverified cells remain qualified.
ax=fig.add_subplot(grid[5]);ax.axis('off');ax.set_title('F  Shape and timing signatures distinguish proposed mechanism classes',loc='left')
classrows=[['Initiation-limited\nGoryachev','Yes','Modest (~1.1×)','Resets','No','Beam width at fixed dose;\npulse trains; gap-end share'],['Mass-conserved, saturating\nWave-pinning, Holmes','No','Yes; Holmes growth\nunresolved at L39','Ratchet; Holmes\nlarge L pending','—','Beam width at fixed dose;\npulse trains; gap-end share'],['Slow non-activity state\ndeBelly*, Spring','No resolved slowing*','Modest (~1×)','Collapses early','Yes','Beam width at fixed dose;\npulse trains; gap-end share']]
table=ax.table(cellText=classrows,colLabels=['Class / models','Initiation slows','Equal-area penalty','Share through gaps','Commit after collapse','Measurement'],cellLoc='center',bbox=[0,.23,1,.70],colWidths=[.22,.12,.17,.16,.12,.21]);table.auto_set_font_size(False);table.set_fontsize(7)
ax.text(0,.12,'◇ Holmes large-L retention is unmeasured; ? width ratio uninformative. Open triangles + arrows denote bounds.',transform=ax.transAxes,fontsize=7.5)
ax.text(0,.03,'* deBelly initiation is interval-censored; later-gap timing and endpoint fields carry numerical qualifications (complete supplement).',transform=ax.transAxes,fontsize=7.3)
for ext in ('png','pdf','svg'):fig.savefig(OUT/f'figure2_round2_pending.{ext}',dpi=180,facecolor='white')
plt.close(fig)
# Supplement additions, followed by all preceding supplement pages and main plate.
with PdfPages(OUT/'supporting/figure2_round2_additions.pdf') as pdf:
 fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title('Holmes large-L retention: frozen design and execution status',loc='left',fontsize=15)
 text='No continuous commitment by 100,000 is an observed finite-horizon finding.\n\nL27 upper committing duration = 18,125; four pulses each 6,796.875.\nTarget lengths: 31, 35, 39; gaps: 0, 0.1, 0.5, 1, 4, 16 × pulse.\nPrimary: D and challenger share at each pulse and gap endpoint. Fate is secondary.\n\nL39 bounded-noise check: gap16, mean 2% native peak, maximum 4%;\n12 independent seeds plus matched-mean and zero controls.\nS3-matched nonperiodic spatial OU input: ell 1.9, tau 10, knot 2.5.\n\n32 SLURM tasks prepared, 1 CPU each, QoS small, 8 GB, 1 hour; no throttle.\nSubmission blocked by sandbox socket permissions. No native trajectories\nwere produced and no simulated retention measurements are assigned.\n\nClass membership is unchanged. The main plate is a pending-data review draft.'
 ax.text(.02,.92,text,va='top',linespacing=1.7,fontsize=11);pdf.savefig(fig);plt.close(fig)
 for m in ORDER:
  fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off');ax.set_title('Class-membership evidence: '+NAMES[m],loc='left',fontsize=15)
  import textwrap
  text='\n\n'.join(r['signature']+f"\nValue: {r['value']}; bounds: [{r['lower']}, {r['upper']}]\n"+textwrap.fill(r['qualification'],100) for r in evidence if r['model']==m)
  ax.text(.02,.92,text,va='top',fontsize=9,linespacing=1.5);pdf.savefig(fig);plt.close(fig)
writer=PdfWriter()
for p in (OUT/'supporting/figure2_round2_additions.pdf',OLD/'supporting/vectors/figure2_signature_assays.pdf',OLD/'supporting/vectors/figure2_signature_assays_complete_supplement.pdf'):
 writer.append(str(p))
writer.write(str(OUT/'figure2_round2_supplement.pdf'))
legend=[
'Enlarged cells resist replacement through distinct initiation, spatial and temporal signatures; the proposed three-class interpretation remains conditional on the pending Holmes large-length retention runs.',
'A, Incumbent amount, extent and peak across the retained 97 states. B, Initiation normalized to each model\u2019s L19 value at central area and the fixed 10% L19 reporter threshold. Brackets propagate saved-sample intervals. deBelly L19 and L39 initiation both lie in (0, 0.5] native units; its ratio is unresolved and is not assigned a point estimate. The open downward symbols show those upper bounds schematically, outside the ratio scale.',
'C, Narrow/broad replacement-time ratio at fixed commanded spatial area. Open triangles and arrows are one-sided bounds; question marks denote uninformative ratios. Holmes large-length width curves do not establish a growing penalty when both arms are censored. Commanded area is not calibrated photon dose.',
'D, Four equal pulses with total ON 1.5 times the model-length commitment calibration; final D after recovery is the colour. Open diamonds identify unexecuted Holmes large-length retention conditions, whose future endpoint D will be distinct from post-recovery fate. No continuous commitment by 100,000 at Holmes L31, L35 and L39 is a finding of extreme finite-horizon replacement resistance.',
'E, L39 gap4 histories. Holmes L39 is pending, with no L27 substitution. Shading marks cue ON; circles and squares mark pulse and gap endpoints. F, Proposed classes with explicit unresolved cells. No membership changes have been made. deBelly later-gap timing and endpoint-field sensitivity remain qualified; its early collapse and final replacement survive the separate quarter-step check. All native kinetics are unchanged. Polarity readouts are distinct from observed migration and conditional OU persistence.'
]
write(OUT/'legend.json',dict(id='figure2',title='Replacement resistance in enlarged cells arises through three mechanism classes; cue shape and timing distinguish their signatures',status='in_progress_pending_Holmes_SLURM',assets=[dict(path=str((OUT/'figure2_round2_pending.png').relative_to(ROOT)),caption='Figure 2 round2 pending Holmes retention',page_size=[12.5,14.5],minimum_printed_width_inches=10,page_title='Figure 2. Shape and timing distinguish mechanism classes')],legend=legend,limitations=['32 targeted Holmes workers prepared but submission blocked; no large-L retention data','deBelly normalized initiation unresolved within first-update intervals','Holmes large-L equal-area growth uninformative; no membership change'],sources=[dict(label='Class evidence CSV',path=str((OUT/'class_membership_check.csv').relative_to(ROOT)))]))
write(OUT/'supporting/metadata/verification.json',dict(generated_utc=datetime.now(timezone.utc).isoformat(),status='qualified_in_progress',reporter_source_audit_verified=True,organization_source_digest_verified=True,source_train_count=len(trains),source_reporter_count=len(reporters),bound_convention='open triangle with arrow; uninformative and unavailable values explicitly marked',palette=COLORS,class_membership_changed=False,complete=False,pending='32 Holmes native workers; SLURM blocked',sha256={str(p.relative_to(OUT)):sha(p) for p in OUT.glob('*') if p.suffix in ('.csv','.pdf','.svg','.png','.json')}))
print(json.dumps(dict(main=str(OUT/'figure2_round2_pending.pdf'),supplement_pages=len(PdfReader(OUT/'figure2_round2_supplement.pdf').pages),status='in_progress')))
