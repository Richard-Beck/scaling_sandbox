"""Native corrected MCA/Rac review: immutable SLURM replays and figure manifests.

launch: submit the 36 uninterrupted-cue records and missing selected pulse histories.
render: independent native checks and publication figures; no integration.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse,copy,gzip,hashlib,json,shlex,subprocess,sys,csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core'),str(ROOT/'code/mechanisms')]
import debelly_corrected_reversal as corrected
import three_test_mechanisms as mech
n=corrected.n
RUN=ROOT/'data/runs/length_quality_debelly_traces_20261009'
OLD=ROOT/'data/runs/length_debelly_review_20261009'
BASE=corrected.DEFAULT
FIG=ROOT/'figures/length_mechanisms_20261008/quality_20261009_debelly'
DATA=ROOT/'data/runs/length_mechanisms_20261008/quality_20261009_debelly'
GEOS=['fixed_physical','fixed_relative','fixed_slope']
TITLES=['Fixed width: 3 µm','Relative width: 20%','Constant slope: 0.24']
LENGTHS=[18,19,21,23,25,27,29,31,33,35,37,39]
COLORS={L:plt.get_cmap('turbo')(.06+.88*i/11) for i,L in enumerate(LENGTHS)}
SOURCES={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':8,'pdf.fonttype':42,'svg.fonttype':'none'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def write(p,a):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(a,indent=2)+'\n')
def load(p,keys=None):
 SOURCES[str(p.relative_to(ROOT))]=sha(p)
 with np.load(p) as z:return {k:z[k].copy() for k in (keys or z.files)}
def csvwrite(p,rs):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def polarity(y,x,L):return y@(2*x/L-1)/y.sum(axis=-1)
def median(y,x,L):
 y=np.atleast_2d(y);s=np.cumsum(y,axis=1);s/=s[:,-1,None];edges=np.r_[0,x+L/len(x)/2]
 return np.array([np.interp(.5,np.r_[0,c],edges) for c in s])
def firstzero(t,y,positive=True):
 inds=np.flatnonzero(y>=0 if positive else y<=0)
 if not len(inds):return None
 i=inds[0]
 return float(t[i]) if i==0 else float(np.interp(0,[y[i-1],y[i]] if positive else [-y[i-1],-y[i]],[t[i-1],t[i]]))
def panel(ax,label,title=None):
 ax.text(-.13,1.04,label,transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans')
 if title:ax.set_title(title,pad=10)
 ax.grid(alpha=.13);ax.axhline(0,color='#8c8c8c',lw=.7)
def save(fig,stem,caption,methods):
 out={}
 for ext in ['png','pdf','svg']:
  p=FIG/f'{stem}.{ext}';fig.savefig(p,dpi=220,facecolor='white');out[p.name]=sha(p)
 plt.close(fig);(FIG/f'{stem}_caption.txt').write_text(caption+'\n')
 write(FIG/f'{stem}.manifest.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,methods=methods,sha256=out,source_sha256=dict(SOURCES),renderer_sha256=sha(Path(__file__)),reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))+' render']))

def launch():
 assert not RUN.exists(),'Fresh destination required'
 old=read(OLD/'manifest.json');m=read(BASE/'manifest.json')
 with gzip.open(BASE/'results/records.json.gz','rt') as f:records=json.load(f)['cases']
 jobs=[]
 for geo in GEOS:
  for L in LENGTHS:
   c=next(c for c in m['cases'] if c['model']=='debelly' and c['condition']==geo and c['length']==L)
   rec=next(a for a in records if a['case_index']==c['case_index'])
   result=rec['result']['orientation'];D=result['endpoint'];expected=next(a for a in rec['tests'] if a['duration']==D)
   previous=next((i for i,j in enumerate(old['jobs']) if j['case']['condition']==geo and j['case']['length']==L and j['purpose']=='pass'),None)
   jobs.append(dict(case=c,duration=D,expected=expected,existing_pulse_slot=previous))
 sources=corrected.sources();sources[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 for name,digest in sources.items():
  src=ROOT/name;dst=RUN/'frozen_sources'/name;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(src.read_bytes())
 for d in ['workers','logs']:(RUN/d).mkdir(parents=True,exist_ok=True)
 write(RUN/'manifest.json',dict(jobs=jobs,sources=sources,center_boundary_motion=False,continuous_horizon=2000.,protocol='Uninterrupted native cue at all 12 eligible corrected lengths in all three geometries. Replay missing selected pulses with 10000 s basal recovery; reuse 24 immutable exactly replayed selected histories. Native step0.5; output interval1.'))
 cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(Path(__file__)),'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"'
 array=subprocess.check_output(['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00','--array=0-35','--job-name=debelly-quality','--output='+str(RUN/'logs/%A_%a.out'),'--error='+str(RUN/'logs/%A_%a.err'),'--wrap',cmd],text=True).strip()
 write(RUN/'launch.json',dict(array=array,jobs=36,generated_at=n.now()));print(array,flush=True)

def worker(slot):
 m=read(RUN/'manifest.json');j=m['jobs'][slot];c=j['case'];L=c['length']
 for name,digest in m['sources'].items():assert sha(ROOT/name)==digest,name
 assert n.model_for('debelly').parameters.center_boundary_motion is False
 assert sha(ROOT/c['checkpoint'])==c['checkpoint_sha256']
 state,_=n.read_checkpoint(ROOT/c['checkpoint'],'debelly');tr=n.run('debelly',L,copy.deepcopy(state),corrected.response_cue(c),2000.,interval=1.)
 arrays={'x':tr.x,'continuous_time':tr.time}
 for k in n.DB_NAMES:arrays['continuous_field_'+k]=tr.fields[k]
 arrays.update({'originalstate_'+k:v for k,v in state.items()})
 metrics=None
 if j['existing_pulse_slot'] is None:
  tr=n.run('debelly',L,copy.deepcopy(state),corrected.response_cue(c),j['duration'],interval=1.)
  rec=n.run('debelly',L,tr.state,n.ConstantInput(0),10000.,interval=2.5)
  metrics=mech.gates(rec,c)
  for key in ['orientation_success','full_polarity_success','field_tail_pass']:assert metrics[key]==j['expected'][key],(slot,key)
  assert abs(metrics['final_orientation']-j['expected']['final_orientation'])<1e-9
  for phase,obj in [('pulse',tr),('recovery',rec)]:
   arrays[phase+'_time']=obj.time
   for k in n.DB_NAMES:arrays[phase+'_field_'+k]=obj.fields[k]
   arrays[phase+'_field_boundary_separation']=obj.fields['boundary_separation']
 p=RUN/f'workers/{slot:03d}.npz';np.savez_compressed(p,**arrays)
 write(p.with_suffix('.json'),dict(job=j,archive_sha256=sha(p),metrics=metrics,generated_at=n.now()));print(slot,c['condition'],L,'finished',flush=True)

def pulse(j,slot):
 if j['existing_pulse_slot'] is None:p=RUN/f'workers/{slot:03d}.npz';meta=read(p.with_suffix('.json'));metrics=meta['metrics']
 else:p=OLD/f'workers/{j["existing_pulse_slot"]:03d}.npz';meta=read(p.with_suffix('.json'));metrics=meta['metrics'];assert abs(meta['withdrawal_time']-j['duration'])<1e-10
 assert sha(p)==meta['archive_sha256'];keys=['x']+[s+'_'+k for s in ['pulse','recovery'] for k in ['time','field_rac','field_mca_bound']]
 return load(p,keys),metrics

def render():
 FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
 for p in [Path(__file__),ROOT/'code/debelly_corrected_reversal.py',ROOT/'code/core/native.py']:SOURCES[str(p.relative_to(ROOT))]=sha(p)
 assert read(BASE/'results/independent_audit.json')['passed'];assert read(OLD/'independent_audit.json')['passed']
 m=read(RUN/'manifest.json');assert m['center_boundary_motion'] is False
 for name,digest in m['sources'].items():assert sha(RUN/'frozen_sources'/name)==digest,name
 execution=dict(frozen_sources_verified=True,execution_source_sha256=m['sources']['code/mechanisms/'+Path(__file__).name],current_renderer_sha256=sha(Path(__file__)),note='Native execution used the retained frozen source snapshot. Subsequent renderer changes affect plots and audits only; reproduce native histories using prepare-native-reproduction.')
 write(DATA/'execution_source_audit.json',execution)
 results=[];traces={};audit=[]
 for slot,j in enumerate(m['jobs']):
  c=j['case'];L=int(c['length']);geo=c['condition'];p=RUN/f'workers/{slot:03d}.npz';meta=read(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256'];assert meta['job']==j
  z=load(p,['x','continuous_time','continuous_field_rac','continuous_field_mca_bound']+['originalstate_'+k for k in n.DB_NAMES]);cp=load(ROOT/c['checkpoint'])
  for k in n.DB_NAMES:np.testing.assert_array_equal(z['originalstate_'+k],cp['state_'+k])
  x=z['x'];t=z['continuous_time'];r=z['continuous_field_rac'];b=z['continuous_field_mca_bound'];assert np.isfinite(r).all() and np.isfinite(b).all() and min(r.min(),b.min())>=-1e-10
  np.testing.assert_allclose(r[0],cp['state_rac'],atol=1e-12);np.testing.assert_allclose(b[0],cp['state_mca_bound'],atol=1e-12)
  a,metrics=pulse(j,slot)
  for field in ['rac','mca_bound']:
   np.testing.assert_allclose(a['pulse_field_'+field][0],cp['state_'+field],atol=1e-12)
   np.testing.assert_allclose(a['recovery_field_'+field][0],a['pulse_field_'+field][-1],atol=1e-12)
  pr=polarity(r,x,L);pm=polarity(b,x,L);bm=median(b,x,L);D=j['duration'];t_m=firstzero(t,pm,False);t_r=firstzero(t,pr);t_half=firstzero(t,bm-L/2,False)
  # Independent orientation/strength gate from archived full recovery.
  tail=a['recovery_time']>=9500;o=polarity(a['recovery_field_rac'][tail],x,L);q=np.ptp(a['recovery_field_rac'][tail],axis=1)
  success=bool(o.min()>.05 and q.min()>=.05 and np.ptp(o)<.005);assert success==metrics['orientation_success'] and success
  # Recompute the separate full-field tail diagnostic on the same frozen scales.
  qp=RUN/f'workers/{slot:03d}.npz' if j['existing_pulse_slot'] is None else OLD/f"workers/{j['existing_pulse_slot']:03d}.npz"
  maximum=0.
  with np.load(qp) as full:
   for key in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','physical_length']:
    values=(full['recovery_field_boundary_separation']+L if key=='physical_length' else full['recovery_field_'+key])[tail]
    assert np.isfinite(values).all()
    if key=='displacement':values=values-values.mean(axis=1)[:,None]
    scale=c['fixed_scales'][key]
    maximum=max(maximum,float(np.ptp(values,axis=0).max()/scale),float(abs(values[-1]-values[-2]).max()/2.5*500/scale))
  assert bool(maximum<=.001)==metrics['field_tail_pass']
  assert abs(maximum-metrics['field_tail_max_scaled'])<1e-10
  cut=polarity(a['pulse_field_mca_bound'][-1],x,L)
  rr=dict(condition=geo,length=L,minimum_reversing_cue=D,rac_zero=t_r,mca_zero=t_m,mca_halfmass_zero=t_half,initial_mca_polarity=pm[0],mca_polarity_at_cut=float(cut),initial_mca_halfmass=float(bm[0]),initial_halfmass_distance=abs(float(bm[0])-L/2),centroid_distance_to_zero=L/2*pm[0],centroid_displacement_under_minimum_cue=L/2*(pm[0]-cut),mean_centroid_speed_to_zero=L/2*pm[0]/t_m if t_m else None,field_tail_pass=metrics['field_tail_pass'],final_rac=float(o[-1]))
  results.append(rr);traces[geo,L]=(z,a,metrics);audit.append(dict(slot=slot,length=L,condition=geo,exact_checkpoint=True,phase_continuity=True,orientation_recomputed=True,field_tail_recomputed=True,field_tail_max_scaled=maximum,archive_sha256=sha(p)))
 csvwrite(DATA/'memory_measurements.csv',results);write(DATA/'independent_audit.json',dict(passed=True,generated_at=n.now(),checks=audit))
 memory_figure(results,traces);continuous_figure(results,traces);recovery_figure(results,traces);flow_figure(results);patch_figure(m)
 write(DATA/'provenance.json',dict(generated_at=n.now(),producer='finished',passed=True,new_continuous_trials=36,new_selected_pulse_replays=12,reused_exact_pulse_histories=24,center_boundary_motion=False,source_sha256=SOURCES,figures={p.name:sha(p) for p in FIG.glob('*.png')}))
 write(RUN/'status.json',dict(generated_at=n.now(),producer='finished',complete=True,jobs=36,independent_audit_passed=True))
 print(FIG,flush=True)

def length_legend(fig,style=False):
 handles=[Line2D([],[],color=COLORS[L],lw=2,label=str(L)) for L in LENGTHS]
 if style:handles +=[Line2D([],[],color='#333333',lw=1.6,label='Cue'),Line2D([],[],color='#333333',lw=1.6,ls='--',label='Recovery'),Line2D([],[],color='#333333',marker='o',ls='',mfc='#333333',label='Stable fields'),Line2D([],[],color='#333333',marker='o',ls='',mfc='white',label='Changing fields')]
 fig.legend(handles=handles,loc='outside lower center',ncol=6 if not style else 8,frameon=False,title='Cell length (µm)'+(' · line style gives protocol phase' if style else ''))
def memory_figure(rr,traces):
 fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained',sharex=True,sharey='row')
 for col,geo in enumerate(GEOS):
  for row,field in enumerate(['rac','mca_bound']):
   ax=axs[row,col];panel(ax,chr(65+row*3+col),TITLES[col] if row==0 else None)
   for L in LENGTHS:
    z,a,metrics=traces[geo,L];rec=next(r for r in rr if r['condition']==geo and r['length']==L);D=rec['minimum_reversing_cue'];sign=1 if field=='rac' else -1
    y=sign*polarity(a['pulse_field_'+field],a['x'],L);yp=sign*polarity(a['recovery_field_'+field],a['x'],L)
    ax.plot(a['pulse_time'],y,color=COLORS[L],lw=1.35);ax.plot(D+a['recovery_time'],yp,color=COLORS[L],lw=1.25,ls='--');ax.plot([D],[y[-1]],'o',ms=3.8,mec=COLORS[L],mfc=COLORS[L] if metrics['field_tail_pass'] else 'white',mew=1,zorder=3)
   ax.set_xlim(0,1800);ax.set_xticks([0,500,1000,1500])
   if row==0:ax.set_ylim(-.16,.46)
   else:ax.set_ylim(-.17,.17)
   if col==0:ax.set_ylabel('Rac polarity toward new front' if row==0 else '−MCA polarity toward new rear')
   if row==1:ax.set_xlabel('Time from cue onset (s)')
 length_legend(fig,True)
 save(fig,'debelly_memory','Figure 8. Rac and attachment redistribution during reversal at all 12 eligible lengths. A–C: Rac polarity; D–F: negative bound-MCA polarity, so a positive value indicates the attachment bias appropriate to the new rear. Columns compare a 3 µm patch, a 20% patch and a gradient of slope 0.24. Each color is one cell length. Solid segments are the shortest tested successful orientation pulse; circles mark withdrawal (filled: stable native-field tail; open: changing fields); dashed segments follow basal recovery. Both rows use the same absolute linear time axis, with the first 1800 s displayed. Reversal is assessed over the final 500 s of a 10000 s recovery, exported in the attached table; the full recovery is separately displayed. MCA evolves on the slower commitment scale after Rac has responded to the cue. Polarities are normalized first moments in nominal reference coordinates. MCA initially marks the old rear, opposite the Rac front. All preparations and endpoint commands use the corrected independent-boundary model.','36 selected pulse histories; 24 existing exact native replays and 12 fresh exact native replays; no MCA-zero withdrawal experiment is included. Every selected pulse passes the independently recomputed orientation, contrast and tail-range gates. Full-field stability remains a separate table column.')
def continuous_figure(rr,traces):
 fig,axs=plt.subplots(2,3,figsize=(12,7.7),layout='constrained',sharex=True,sharey='row')
 for col,geo in enumerate(GEOS):
  for row,field in enumerate(['rac','mca_bound']):
   ax=axs[row,col];panel(ax,chr(65+row*3+col),TITLES[col] if row==0 else None)
   for L in LENGTHS:
    z,a,metrics=traces[geo,L];sign=1 if row==0 else -1;ax.plot(z['continuous_time'],sign*polarity(z['continuous_field_'+field],z['x'],L),color=COLORS[L],lw=1.4)
   ax.set_xlim(0,1200);ax.set_xticks([0,300,600,900,1200]);
   if col==0:ax.set_ylabel('Rac polarity toward new front' if row==0 else '−MCA polarity toward new rear')
   if row==1:ax.set_xlabel('Time under uninterrupted cue (s)')
 length_legend(fig)
 save(fig,'debelly_continuous','Supplement. Rac (A–C) and negative MCA (D–F) polarity under uninterrupted reversal cues at every eligible length. Colors and cue columns match Figure 8; all time axes are linear. Native cues remain on for 2000 s; the first 1200 s is displayed. These trajectories define the cue-only zero-crossing clocks in Figure 9 and are distinct from the finite pulses followed by recovery in Figure 8.','Fresh SLURM integration at 36 corrected checkpoint–cue combinations; native step0.5 and output interval1 s. First crossings are linearly interpolated between native outputs; absence through 2000 s is retained as missing rather than extrapolated.')
def recovery_figure(rr,traces):
 fig,axs=plt.subplots(2,3,figsize=(12,7.7),layout='constrained',sharex=True,sharey='row')
 for col,geo in enumerate(GEOS):
  for row,field in enumerate(['rac','mca_bound']):
   ax=axs[row,col];panel(ax,chr(65+row*3+col),TITLES[col] if row==0 else None)
   for L in LENGTHS:
    z,a,metrics=traces[geo,L];sign=1 if row==0 else -1;ax.plot(a['recovery_time'],sign*polarity(a['recovery_field_'+field],a['x'],L),color=COLORS[L],lw=1.3)
   ax.set_xlim(0,10000);ax.set_xticks([0,2500,5000,7500,10000]);
   if col==0:ax.set_ylabel('Rac polarity toward new front' if row==0 else '−MCA polarity toward new rear')
   if row==1:ax.set_xlabel('Time since withdrawal (s)')
 length_legend(fig)
 save(fig,'debelly_recovery','Supplement. Complete basal recovery after the shortest tested reversing pulse at every eligible length. A–C show Rac polarity, D–F sign-reversed bound-MCA polarity; colors and cue columns match Figure 8. All selected trials pass the orientation gate over the final 500 s. This scalar persistence does not require all native chemical and mechanical fields to have equilibrated: their independent tail flags remain in the measurement table.','All 10000 s of archived native recovery, sampled every2.5 s. No time-axis compression or extrapolation.')

def flow_figure(rr):
 from matplotlib.ticker import ScalarFormatter
 FLOW=ROOT/'data/runs/debelly_corrected_advection_20261008';REV=ROOT/'data/runs/debelly_corrected_advection_revision_20261008';OLDMEAS=ROOT/'data/runs/length_mechanisms_20261008/review_20261009_debelly'
 assert read(FLOW/'independent_audit.json')['passed'];assert read(REV/'controlled_calculation_audit.json')['passed']
 fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained');a=[r for r in rr if r['condition']=='fixed_slope'];xx=[r['length'] for r in a]
 for ax,letter in zip(axs.flat,'ABCDEF'):panel(ax,letter)
 for key,color,ls,marker,label in [('rac_zero','#777777','-',None,'Rac zero'),('mca_zero','#762a83','-','o','MCA zero'),('minimum_reversing_cue','#0072B2','-','s','Minimum reversing cue')]:
  axs[0,0].plot(xx,[r[key] for r in a],color=color,ls=ls,marker=marker,ms=4,label=label)
 axs[0,0].set_ylabel('Time (s)');axs[0,0].legend(frameon=False,fontsize=8)
 for key,color,label in [('initial_halfmass_distance','#009E73','Initial half-mass distance to midpoint'),('centroid_distance_to_zero','#762a83','Initial centroid distance to midpoint'),('centroid_displacement_under_minimum_cue','#0072B2','Centroid shift during reversing cue')]:axs[0,1].plot(xx,[r[key] for r in a],color=color,marker='o',ms=4,label=label)
 axs[0,1].set_ylabel('MCA distance / centroid shift (µm)');axs[0,1].legend(loc='lower left',bbox_to_anchor=(0,.05),frameon=True,framealpha=1,facecolor='white',edgecolor='none',fontsize=7.5)
 axs[0,2].plot(xx,[r['mean_centroid_speed_to_zero']*1000 for r in a],color='#009E73',marker='o',ms=4);axs[0,2].set_ylabel('Mean speed to zero (10⁻³ µm/s)')
 for ax in axs[0]:ax.set_xlabel('Cell length (µm)');ax.set_xlim(17,40)
 p=n.model_for('debelly').parameters;b0=.078125;xi=p.linker_friction*b0;z=load(FLOW/'L23_profiles.npz');refx=z['zero_x_fraction']*23;refrho=z['zero_rho_concentration'];checks=[]
 for L,color in zip([23,39,60,100],['#0072B2','#009E73','#D55E00','#762a83']):
  speeds=[]
  for N in [2001,4001]:
   dx=L/N;x=(np.arange(N)+.5)*dx;rho=np.interp(x,refx,refrho);force=p.active_stress*p.rho_to_rac_coupling*n.ds._gradient(rho,dx);v=-n.ds._dirichlet_elliptic(np.full(N,xi),p.cortical_viscosity,force,dx,0.,0.);speeds.append(np.interp(5,x,v))
   if N==2001:axs[1,0].plot(np.r_[0,x,L],np.r_[0,v,0],color=color,label=f'{L} µm')
  assert abs(speeds[0]/speeds[1]-1)<.001;checks.append(dict(length=L,relative_mesh_change=abs(speeds[0]/speeds[1]-1)))
 axs[1,0].axvline(5,color='#999999',ls=':',lw=.8);axs[1,0].set_xlabel('Distance from original left end (µm)');axs[1,0].set_ylabel('Rho-only leftward speed (µm/s)');axs[1,0].legend(frameon=False,title='Extended domain',fontsize=7)
 def rows(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return list(csv.DictReader(p.open()))
 controls=rows(REV/'fixed_reference_rho_extension.csv')
 for friction,color,ls,label in [('True','#0072B2','-','MCA friction'),('False','#777777','--','Zero friction')]:
  b=[r for r in controls if r['friction']==friction];axs[1,1].plot([float(r['length']) for r in b],[float(r['speed_at_x5']) for r in b],color=color,ls=ls,label=label)
 axs[1,1].set_xlabel('Extended domain length (µm)');axs[1,1].set_ylabel('Rho-only speed at x = 5 µm (µm/s)');axs[1,1].legend(frameon=False,fontsize=8)
 budget=rows(OLDMEAS/'zero_rate_budget.csv');assert max(abs(float(r['budget_residual'])) for r in budget)<1e-11
 for key,color,ls,label in [('rate_boundary','#D55E00','-','Endpoint advection'),('rate_rho','#0072B2','-','Rho advection'),('other','#999999','-','Other terms'),('reorientation_rate','#222222','--','Net rate')]:
  yy=[sum(float(r['rate_'+k]) for k in ['membrane','diffusion','exchange','replenishment']) if key=='other' else float(r[key]) for r in budget];axs[1,2].plot([float(r['length']) for r in budget],np.array(yy)*1e4,color=color,ls=ls,label=label)
 axs[1,2].set_xlabel('Cell length (µm)');axs[1,2].set_ylabel('−dP(MCA)/dt (10⁻⁴ s⁻¹)');axs[1,2].legend(frameon=False,fontsize=7.5)
 cap='Figure 9. Attachment travel, speed and mechanical screening. A–C use the constant-slope cue and all 12 eligible lengths. A compares the first uninterrupted-cue Rac and MCA zero crossings with the shortest tested reversing pulse. B compares initial MCA half-mass and centroid distances to the midpoint with the centroid displacement during that actual pulse. C divides distance to MCA zero by its uninterrupted-cue crossing time; it is a mean distribution speed, not a material-particle velocity. D–E are frozen mechanical controls: the L = 23 Rho profile at MCA zero is held fixed in physical reference coordinates and its low-Rho tail is extended. Membrane and endpoint velocities are zero; uniform basal bound MCA fixes friction. D shows the resulting cortical speed profiles. E samples the same static solution at x = 5 µm with MCA friction (blue) or zero friction (gray). Friction screens sensitivity to the remote end. F is the independently audited instantaneous native MCA polarity-rate decomposition at the first MCA zero. Endpoint- and Rho-driven advection, plus membrane motion, diffusion, binding/unbinding and restoration (grouped as other), sum to the net rate with changing MCA normalization included. Positive values erase the old attachment bias. F is an additive source identity evaluated on each native state; D–E are separately solved frozen controls. Initial preparations below 18 µm are ineligible; every eligible length has a reversing tested pulse. Full-field tail flags are in the attached table.'
 save(fig,'debelly_flow',cap,dict(static_control_mesh_checks=checks,native_source_budget_max_residual=max(abs(float(r['budget_residual'])) for r in budget),native_budget_source='Previously independently audited corrected native zero-rate budget, unchanged.'))
 # Compare the scalar clocks without selecting a favorable subset of lengths.
 errors=[];fig,axs=plt.subplots(1,3,figsize=(12,3.9),layout='constrained',sharey=True)
 for col,geo in enumerate(GEOS):
  ax=axs[col];panel(ax,chr(65+col),TITLES[col]);b=[r for r in rr if r['condition']==geo];x=[r['length'] for r in b]
  for key,color,ls,label in [('rac_zero','#777777','-','Rac zero'),('mca_zero','#762a83','-','MCA zero'),('mca_halfmass_zero','#009E73',':','MCA half-mass at midpoint'),('minimum_reversing_cue','#0072B2','-','Minimum reversing cue')]:
   ax.plot(x,[r[key] if r[key] is not None else np.nan for r in b],color=color,ls=ls,marker='s' if key=='minimum_reversing_cue' else None,label=label)
   if key!='minimum_reversing_cue':
    valid=[r for r in b if r[key] is not None and r[key]>0];d=np.array([abs(np.log(r[key]/r['minimum_reversing_cue'])) for r in valid]);errors.append(dict(condition=geo,readout=key,available_lengths=len(valid),median_absolute_log_ratio=float(np.median(d)),median_factor_difference=float(np.exp(np.median(d))),root_mean_square_time_error=float(np.sqrt(np.mean([(r[key]-r['minimum_reversing_cue'])**2 for r in valid])))))
  ax.set_xlabel('Cell length (µm)');ax.set_xlim(17,40);ax.set_ylim(bottom=0)
 axs[0].set_ylabel('Time (s)');axs[1].legend(frameon=False,fontsize=7.5)
 csvwrite(DATA/'scalar_clock_comparison.csv',errors)
 save(fig,'debelly_clock_comparison','Supplement. Approximate redistribution clocks under all three cues. A–C compare uninterrupted-cue Rac zero, bound-MCA zero, MCA half-mass arrival at the midpoint, and the minimum reversing pulse duration. MCA clocks are much closer to the duration scale than the rapid Rac sign change, while their quantitative agreement depends on cue and length. No clock is fitted to outcomes. The half-mass landmark is an additional candidate scalar, evaluated by the same absolute log-ratio timing error at all available eligible lengths; the table reports missing crossings and both relative and absolute errors. A reversal pulse is assessed after recovery, whereas each crossing here is measured during an uninterrupted cue.','Independent native normalized first moments and piecewise-uniform finite-volume median; native outputs every1 s, crossings linearly interpolated. Timing-error summaries are descriptive, not calibrated predictors.')

def patch_figure(m):
 old=read(OLD/'manifest.json');rows=[];fig,axs=plt.subplots(3,2,figsize=(11,10),layout='constrained',sharex=True)
 for col,geo in enumerate(GEOS[:2]):
  data=[]
  for L in [25,27,29,31,33,35]:
   slot=next(i for i,j in enumerate(m['jobs']) if j['case']['condition']==geo and j['case']['length']==L);j=m['jobs'][slot];z,metrics=pulse(j,slot);t=z['pulse_time'];pm=polarity(z['pulse_field_mca_bound'],z['x'],L);post=polarity(z['recovery_field_mca_bound'],z['x'],L);v=L/2*(np.interp(199,t,pm)-np.interp(201,t,pm))/2
   failure=next(i for i,a in enumerate(old['jobs']) if a['case']['condition']==geo and a['case']['length']==L and a['purpose']=='fail');failmeta=read(OLD/f'workers/{failure:03d}.json')
   row=dict(condition=geo,length=L,duration=j['duration'],failing_duration=failmeta['withdrawal_time'],mca_at_cut=float(pm[-1]),mca_at_200=float(np.interp(200,t,pm)),centroid_speed_at_200=float(v),centroid_shift_first_100_postcue=float(L/2*(post[0]-np.interp(100,z['recovery_time'],post))),field_tail_pass=metrics['field_tail_pass']);data.append(row);rows.append(row)
  x=[r['length'] for r in data]
  for row in range(3):panel(axs[row,col],chr(65+row*2+col),TITLES[col] if row==0 else None)
  ax=axs[0,col];ax.plot(x,[r['duration'] for r in data],color='#0072B2',label='Minimum reversing cue');
  for r in data:ax.plot(r['length'],r['duration'],'o',mec='#0072B2',mfc='#0072B2' if r['field_tail_pass'] else 'white',ms=5)
  ax.plot(x,[r['failing_duration'] for r in data],color='#D55E00',marker='x',ls='',label='Adjacent failing cue');ax.set_ylabel('Pulse duration (s)');ax.set_ylim((290,375) if col==0 else (235,300));ax.set_yticks([300,325,350,375] if col==0 else [240,260,280,300]);ax.legend(frameon=False,fontsize=8)
  ax=axs[1,col];ax.plot(x,[r['mca_at_200'] for r in data],color='#777777',marker='s',label='At common cue time: 200 s');ax.plot(x,[r['mca_at_cut'] for r in data],color='#762a83',marker='o',label='At reversing-cue withdrawal');ax.set_ylabel('MCA polarity (old rear bias)');ax.legend(frameon=False,fontsize=7.5)
  ax=axs[2,col];ax.plot(x,[r['centroid_shift_first_100_postcue'] for r in data],color='#D55E00',marker='o',label='MCA centroid shift after withdrawal');ax.set_ylabel('MCA shift in next 100 s (µm)');ax.legend(frameon=False,fontsize=7.5)
  for ax in axs[:,col]:ax.set_xticks(x);ax.set_xlabel('Cell length (µm)');ax.set_xlim(24,36)
 csvwrite(DATA/'patch_nonmonotonicity.csv',rows)
 save(fig,'debelly_patch_turnovers','Figure 11. Attachment diagnostics around the patch-duration decreases at L = 25–35 µm. Columns show fixed 3 µm and relative 20% patches. A–B compare the shortest tested reversing cue with its adjacent failing duration; open successful markers retain a failed native-field tail flag. C–D compare MCA polarity at the common cue time of 200 s and at the successful cue withdrawal. Both the common-time residual attachment bias and its value at successful withdrawal vary with length; the cutoff is not a fixed MCA polarity. E–F show MCA centroid displacement during the first 100 s of basal recovery after that successful cutoff. The fixed-width dip coincides with increased post-cue progress; the relative-width dip has different post-cue behavior. These separate trends locate the turnover within coupled cue and recovery dynamics without assigning it to one cause. Full-field tail flags and the additional common-time centroid speeds are exported in the attached table.','Same exact corrected native selected and adjacent-failure replays as Figure 8; no new intervention, extrapolation or fitted turnover. Piecewise linear interpolation of native measurements at 199,200,201 and100 s.')

def prepare_native_reproduction(destination):
 """Reconstitute the original executable directory layout without changing evidence."""
 assert not destination.exists(),'Fresh reproduction root required'
 m=read(RUN/'manifest.json')
 for name,digest in m['sources'].items():assert sha(RUN/'frozen_sources'/name)==digest
 shutil.copytree(RUN/'frozen_sources',destination)
 # The imported gate helper was not in the launch snapshot; independently match
 # its prior immutable source archive before including the dependency closure.
 refs=read(OLD/'manifest.json')['sources']
 support={}
 for name in ['code/mechanisms/three_test_mechanisms.py','code/three_test_reversal.py']:
  src=ROOT/name
  if name in refs:assert sha(src)==refs[name]
  else:
   prior=ROOT/'data/runs/three_test_mechanisms_20261006/frozen_sources'/name
   assert prior.exists() and sha(prior)==sha(src),name
  target=destination/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,target);support[name]=sha(src)
 for j in m['jobs']:
  name=j['case']['checkpoint'];assert sha(ROOT/name)==j['case']['checkpoint_sha256']
  dst=destination/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dst)
 out=destination/'data/runs/length_quality_debelly_traces_20261009'
 for d in ['workers','logs']:(out/d).mkdir(parents=True,exist_ok=True)
 shutil.copy2(RUN/'manifest.json',out/'manifest.json')
 runner=destination/'code/mechanisms'/Path(__file__).name
 command='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(runner),'worker'])+' --slot "$SLURM_ARRAY_TASK_ID"'
 submit=shlex.join(['sbatch','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00','--array=0-35','--output='+str(out/'logs/%A_%a.out'),'--error='+str(out/'logs/%A_%a.err'),'--wrap',command])
 write(destination/'REPRODUCTION.json',dict(original_manifest_sha256=sha(RUN/'manifest.json'),execution_sources_verified=True,supplementary_imported_sources_matched_to_prior_immutable_archives=support,submission_command=submit,note='Original frozen executable layout and exact corrected checkpoints. Submission is explicit and writes only this fresh destination.'))
 print(submit)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['launch','worker','render','prepare-native-reproduction']);p.add_argument('--slot',type=int);p.add_argument('--reproduction-root',type=Path);a=p.parse_args()
 if a.action=='launch':launch()
 elif a.action=='worker':worker(a.slot)
 elif a.action=='prepare-native-reproduction':prepare_native_reproduction(a.reproduction_root or RUN/'native_reproduction')
 else:render()
