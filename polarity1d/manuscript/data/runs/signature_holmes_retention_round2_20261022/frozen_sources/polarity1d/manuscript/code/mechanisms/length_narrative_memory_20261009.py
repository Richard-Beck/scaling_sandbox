"""Strength-selected native Rac/MCA trajectories and redistribution clocks.

`launch` runs only missing exact pulse histories, with immutable native sources.
`render` audits and plots the retained/native histories without long integration.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, copy, csv, gzip, hashlib, json, shlex, shutil, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core'),str(ROOT/'code/mechanisms')]
import length_quality_debelly_traces_20261009 as previous
import debelly_corrected_advection as ad
n=previous.n
RUN=ROOT/'data/runs/length_narrative_memory_20261009'
FIG=ROOT/'figures/length_mechanisms_20261008/narrative_20261009/memory'
DATA=ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/memory'
LENGTHS=previous.LENGTHS; GEOS=previous.GEOS; TITLES=previous.TITLES
COLORS={L:plt.get_cmap('turbo')(.06+.88*i/11) for i,L in enumerate(LENGTHS)}
SOURCES={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','legend.fontsize':8})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def register(p):SOURCES[str(Path(p).relative_to(ROOT))]=sha(p)
def read(p):register(p);return json.loads(p.read_text())
def write(p,a):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(a,indent=2)+'\n')
def load(p,keys=None):
 register(p)
 with np.load(p) as z:return {k:z[k].copy() for k in (keys or z.files)}
def csvwrite(p,rows):
 with p.open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
def select():
 old=read(previous.RUN/'manifest.json');register(previous.BASE/'results/records.json.gz')
 with gzip.open(previous.BASE/'results/records.json.gz','rt') as stream: records=json.load(stream)['cases']
 jobs=[]
 for slot,j in enumerate(old['jobs']):
  c=j['case'];r=next(a for a in records if a['case_index']==c['case_index']);D=r['result']['full_polarity']['endpoint']
  expected=next(a for a in r['tests'] if a['duration']==D)
  assert D==min(a['duration'] for a in r['tests'] if a['full_polarity_success'])
  if D==j['duration']:
   source=previous.RUN/f'workers/{slot:03d}.npz' if j['existing_pulse_slot'] is None else previous.OLD/f'workers/{j["existing_pulse_slot"]:03d}.npz'
  else:source=None
  jobs.append(dict(case=c,duration=D,expected=expected,quality_slot=slot,source=str(source.relative_to(ROOT)) if source else None))
 return jobs
def launch():
 assert not RUN.exists(),'Use fresh run destination'
 jobs=select();missing=[j for j in jobs if j['source'] is None]
 for folder in ['workers','logs']:(RUN/folder).mkdir(parents=True)
 # A complete code snapshot avoids imports from mutable live sources.
 shutil.copytree(ROOT/'code',RUN/'frozen_sources/code',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 snapshots={str(p.relative_to(RUN/'frozen_sources')):sha(p) for p in (RUN/'frozen_sources/code').rglob('*.py')}
 write(RUN/'manifest.json',dict(jobs=jobs,missing_jobs=missing,source_sha256=snapshots,center_boundary_motion=False,native_step=.5,selection='Shortest tested 1%-strength-passing cue; exact corrected basal checkpoint; 10000 s recovery.'))
 worker='''from pathlib import Path
import sys,json,copy,hashlib
import numpy as np
run=Path(sys.argv[1]);root=Path(sys.argv[2]);slot=int(sys.argv[3])
sys.path[:0]=[str(run/'frozen_sources/code'),str(run/'frozen_sources/code/core'),str(run/'frozen_sources/code/mechanisms')]
import debelly_corrected_reversal as c
import three_test_mechanisms as mech
m=json.loads((run/'manifest.json').read_text());j=m['missing_jobs'][slot];case=j['case'];L=case['length']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name,digest in m['source_sha256'].items():assert sha(run/'frozen_sources'/name)==digest,name
cp=root/case['checkpoint'];assert sha(cp)==case['checkpoint_sha256'];assert not c.n.model_for('debelly').parameters.center_boundary_motion
state,_=c.n.read_checkpoint(cp,'debelly');a=c.n.run('debelly',L,copy.deepcopy(state),c.response_cue(case),j['duration'],interval=.5)
b=c.n.run('debelly',L,a.state,c.n.ConstantInput(0),10000.,interval=2.5);metrics=mech.gates(b,case)
for key in ['orientation_success','full_polarity_success','field_tail_pass']:assert metrics[key]==j['expected'][key],(slot,key)
assert abs(metrics['final_orientation']-j['expected']['final_orientation'])<1e-9
arrays={'x':a.x};arrays.update({'originalstate_'+k:v for k,v in state.items()})
for phase,tr in [('pulse',a),('recovery',b)]:
 arrays[phase+'_time']=tr.time
 for key in c.n.DB_NAMES+['boundary_separation']:arrays[phase+'_field_'+key]=tr.fields[key]
p=run/f'workers/{slot:03d}.npz';np.savez_compressed(p,**arrays)
p.with_suffix('.json').write_text(json.dumps(dict(job=j,metrics=metrics,archive_sha256=sha(p),generated_at=c.n.now()),indent=2)+'\\n')
print(slot,case['condition'],L,'finished',flush=True)
'''
 workerpath=RUN/'worker.py';workerpath.write_text(worker)
 cmd='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1; '+shlex.join([sys.executable,str(workerpath),str(RUN),str(ROOT)])+' "$SLURM_ARRAY_TASK_ID"'
 array=subprocess.check_output(['sbatch','--parsable','--qos=small','--ntasks=1','--cpus-per-task=1','--mem=3G','--time=01:00:00',f'--array=0-{len(missing)-1}','--job-name=debelly-strength','--output='+str(RUN/'logs/%A_%a.out'),'--error='+str(RUN/'logs/%A_%a.err'),'--wrap',cmd],text=True).strip()
 write(RUN/'launch.json',dict(array=array,jobs=len(missing),reused=36-len(missing),worker_sha256=sha(workerpath),generated_at=n.now()));print(array,flush=True)
def pol(y,x,L):return y@(2*x/L-1)/y.sum(axis=-1)
def endpoint(y):return y[:,0]-y[:,-1]
def annotate(ax,label,title=None):
 ax.text(-.14,1.025,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
 if title:ax.set_title(title,pad=12)
 ax.axhline(0,lw=.65,color='#8a8a8a');ax.grid(alpha=.12)
def save(fig,stem,caption,methods):
 hashes={}
 for ext in ['pdf','svg','png']:
  p=FIG/f'{stem}.{ext}';fig.savefig(p,dpi=190,facecolor='white');hashes[p.name]=sha(p)
 plt.close(fig);(FIG/f'{stem}_caption.txt').write_text(caption+'\n')
 write(FIG/f'{stem}.manifest.json',dict(producer='finished',passed=True,generated_at=n.now(),legend=caption,description=caption,sha256=hashes,source_sha256=dict(SOURCES),renderer_sha256=sha(Path(__file__)),methods=methods,reproduction=['OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))+' render']))
def color_legend(fig):
 fig.legend(handles=[Line2D([],[],color=COLORS[L],lw=2,label=str(L)) for L in LENGTHS],loc='outside lower center',ncol=12,frameon=False,title='Cell length (µm)',handlelength=1.35,columnspacing=.9)
def trajectories(traces,complete=False,continuous=False):
 fig,axs=plt.subplots(3,3,figsize=(12,9.5),sharex=True,sharey='row',layout='constrained')
 for col,geo in enumerate(GEOS):
  for row in range(3):
   ax=axs[row,col];annotate(ax,chr(65+3*row+col),TITLES[col] if row==0 else None)
   for L in LENGTHS:
    a,z,j,metrics=traces[geo,L];D=j['duration']
    calc=lambda y:pol(y,a['x'],L)*(1 if row==0 else -1) if row<2 else endpoint(y)
    key='rac' if row==0 else 'mca_bound'
    if continuous:ax.plot(z['continuous_time'],calc(z['continuous_field_'+key]),color=COLORS[L],lw=1.2)
    elif complete:
     ax.plot(a['recovery_time'],calc(a['recovery_field_'+key]),color=COLORS[L],lw=1.2)
    else:
     y=calc(a['pulse_field_'+key]);ax.plot(a['pulse_time'],y,color=COLORS[L],lw=1.3)
     ax.plot(D+a['recovery_time'],calc(a['recovery_field_'+key]),color=COLORS[L],ls=(0,(3,2)),lw=1.15)
     ax.plot(D,y[-1],'o',ms=3.4,mfc='white',mec=COLORS[L],mew=1,zorder=4)
   if complete:
    ax.set_xscale('symlog',linthresh=2.5,linscale=.6);ax.set_xlim(0,10000);ax.set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000'])
   else:ax.set_xlim(0,1200);ax.set_xticks([0,300,600,900,1200])
   if col==0:ax.set_ylabel(['Rac polarity\n(new front →)','−MCA spatial polarity\n(new rear →)','Endpoint MCA difference\nleft − right (native conc.)'][row])
   if row==2:ax.set_xlabel('Time since withdrawal (s)' if complete else 'Time under uninterrupted cue (s)' if continuous else 'Time from cue onset (s)')
 color_legend(fig)
 if not complete and not continuous:
  fig.legend(handles=[Line2D([],[],color='#333333',lw=1.4,label='Cue on'),Line2D([],[],color='#333333',ls=(0,(3,2)),lw=1.4,label='Basal recovery'),Line2D([],[],color='#333333',marker='o',ls='',mfc='white',ms=4,label='Withdrawal')],loc='outside upper center',ncol=3,frameon=False)
 common='Columns compare fixed support 3 µm, relative support 20%, and constant slope 0.24. Colors identify all twelve eligible lengths. A–C: normalized Rac first moment toward the new right front. D–F: negative normalized bound-MCA first moment, positive toward the new left rear. G–I: left-minus-right bound MCA at the two endpoint-adjacent finite-volume samples, in native concentration units. The spatial moment measures whole-cell redistribution; the endpoint difference measures local end-state exchange. All runs start at their exact corrected relative-20% basal checkpoint with independent endpoint commands.'
 if complete:
  stem='debelly_recovery';cap='Complete basal recovery after each shortest tested 1%-strength-passing cue. '+common+' Time after withdrawal is linear from 0 to 2.5 s, then logarithmic; all 10000 s are retained.'
 elif continuous:
  stem='debelly_continuous';cap='Uninterrupted-cue histories defining the redistribution clocks. '+common+' Cues stay on for 2000 s; the first 1200 s is displayed on linear time. First crossings are interpolated between retained 1-s samples.'
 else:
  stem='debelly_memory';cap='Local endpoint exchange accompanies slower spatial attachment redistribution. '+common+' Solid curves show each length’s shortest tested cue satisfying the canonical 1% scalar strength gate; circles mark withdrawal, and dashed curves show basal recovery. The first 1200 s use linear time, retaining every selected withdrawal. The complete 10000-s recovery is supplied separately on logarithmic time. Selection requires Rac orientation, contrast, final-window persistence and reflected polarity strength within 1%; native-field tail flags are separately tabulated. No field is reset at withdrawal.'
 save(fig,stem,cap,dict(selection='Shortest tested full_polarity_success (canonical scalar strength within1%).',replay_jobs=21,reused_histories=15,full_recovery_s=10000,coordinates='nominal reference cell-centers',endpoint_measure='leftmost bound-MCA concentration minus rightmost concentration'))
 return stem,cap
def render():
 FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
 for p in [Path(__file__),Path(previous.__file__),Path(ad.__file__),Path(n.__file__),Path(n.ds.__file__)]:register(p)
 m=read(RUN/'manifest.json');jobs=select();assert jobs==m['jobs'];assert not n.model_for('debelly').parameters.center_boundary_motion
 assert read(previous.BASE/'results/independent_audit.json')['passed']
 assert read(previous.RUN/'status.json')['independent_audit_passed']
 for key,digest in m['source_sha256'].items():assert sha(RUN/'frozen_sources'/key)==digest
 results=[];traces={};audits=[]
 keys=['x']+[p+'_'+k for p in ['pulse','recovery'] for k in ['time','field_rac','field_mca_bound']]
 for j in jobs:
  c=j['case'];L=int(c['length']);geo=c['condition'];cp=load(ROOT/c['checkpoint']);assert sha(ROOT/c['checkpoint'])==c['checkpoint_sha256']
  if j['source']:p=ROOT/j['source']
  else:p=RUN/f'workers/{m["missing_jobs"].index(j):03d}.npz'
  meta=read(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256'];a=load(p,keys);metrics=meta['metrics']
  q=previous.RUN/f'workers/{j["quality_slot"]:03d}.npz';qm=read(q.with_suffix('.json'));assert sha(q)==qm['archive_sha256']
  z=load(q,['continuous_time','continuous_field_rac','continuous_field_mca_bound'])
  for key in ['rac','mca_bound']:
   np.testing.assert_allclose(a['pulse_field_'+key][0],cp['state_'+key],atol=1e-12)
   np.testing.assert_allclose(a['recovery_field_'+key][0],a['pulse_field_'+key][-1],atol=1e-12)
   np.testing.assert_allclose(z['continuous_field_'+key][0],cp['state_'+key],atol=1e-12)
  tail=a['recovery_time']>=9500;r=pol(a['recovery_field_rac'][tail],a['x'],L);contrast=np.ptp(a['recovery_field_rac'][tail],axis=1)
  maximum=0.
  with np.load(p) as full:
   for key in n.DB_NAMES:
    np.testing.assert_allclose(full['pulse_field_'+key][0],cp['state_'+key],atol=1e-12)
    np.testing.assert_allclose(full['recovery_field_'+key][0],full['pulse_field_'+key][-1],atol=1e-12)
   for key in ['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','physical_length']:
    values=(full['recovery_field_boundary_separation']+L if key=='physical_length' else full['recovery_field_'+key])[tail]
    assert np.isfinite(values).all()
    if key=='displacement':values=values-values.mean(axis=1)[:,None]
    scale=c['fixed_scales'][key]
    maximum=max(maximum,float(np.ptp(values,axis=0).max()/scale),float(abs(values[-1]-values[-2]).max()/2.5*500/scale))
  assert bool(maximum<=.001)==metrics['field_tail_pass']
  assert abs(maximum-metrics['field_tail_max_scaled'])<1e-10
  orientation=bool(r.min()>.05 and contrast.min()>=.05 and np.ptp(r)<.005);error=float(abs(r/abs(c['initial_orientation'])-1).max());assert orientation and error<=.01
  assert abs(r[-1]-j['expected']['final_orientation'])<1e-9
  assert metrics['full_polarity_success'] and metrics['field_tail_pass']==j['expected']['field_tail_pass']
  t=z['continuous_time'];rac=pol(z['continuous_field_rac'],a['x'],L);mca=pol(z['continuous_field_mca_bound'],a['x'],L);end=endpoint(z['continuous_field_mca_bound']);mzero=previous.firstzero(t,mca,False)
  cut=float(pol(a['pulse_field_mca_bound'][-1],a['x'],L));bm=previous.median(z['continuous_field_mca_bound'][0],a['x'],L)[0]
  row=dict(condition=geo,length=L,minimum_strength_cue=j['duration'],rac_zero=previous.firstzero(t,rac),mca_zero=mzero,endpoint_mca_zero=previous.firstzero(t,end),initial_mca_polarity=float(mca[0]),mca_polarity_at_cut=cut,initial_halfmass_distance=float(abs(bm-L/2)),centroid_distance_to_zero=float(L/2*mca[0]),centroid_shift_during_strength_cue=float(L/2*(mca[0]-cut)),mean_centroid_speed_to_zero=float(L/2*mca[0]/mzero),field_tail_pass=metrics['field_tail_pass'],maximum_relative_strength_error=error,final_rac=float(r[-1]),source_archive=str(p.relative_to(ROOT)),checkpoint_sha256=c['checkpoint_sha256'])
  results.append(row);traces[geo,L]=(a,z,j,metrics)
  audits.append(dict(length=L,condition=geo,exact_checkpoint=True,phase_continuity=True,strength_gate_independently_recomputed=True,maximum_relative_strength_error=error,minimum_tested_strength_duration_verified=True,field_tail_independently_recomputed=True,field_tail_max_scaled=maximum,archive_sha256=sha(p)))
 csvwrite(DATA/'memory_measurements.csv',results);write(DATA/'independent_audit.json',dict(passed=True,generated_at=n.now(),checks=audits,frozen_execution_sources_verified=True))
 for complete,continuous in [(False,False),(True,False),(False,True)]:trajectories(traces,complete,continuous)
 flow_cap=flow(results,traces);budget_cap=budgets(traces)
 write(DATA/'chapter_overrides.json',dict(mca=dict(image='narrative_20261009/memory/debelly_memory.png',renderer=Path(__file__).name,title='Rac switches before attachment redistribution is complete',paragraphs=['The cue first recruits Rac at the opposite end. Bound MCA responds more slowly: local attachment exchange at the two ends and redistribution across the whole cell progress on different trajectories. Across all three cue geometries, the length-colored histories make this separation visible while retaining the same strength-based definition of successful reversal.', 'Cue withdrawal removes the external input without resetting the inherited chemical or mechanical fields. The complete recovery shows how the reversed Rac pole and its opposing attachment bias subsequently rebuild. Longer cells therefore enter and leave the cue with distinct attachment distributions as well as distinct resting mechanical states.'],caption=(FIG/'debelly_memory_caption.txt').read_text().strip(),sources=links('debelly_memory')+links('debelly_recovery')+links('debelly_continuous')+[['../../data/runs/length_mechanisms_20261008/narrative_20261009/memory/memory_measurements.csv','Strength selection, clocks and full-field flags']]),flow=dict(image='narrative_20261009/memory/debelly_flow.png',renderer=Path(__file__).name,title='Attachment redistribution carries a length-dependent distance and time',paragraphs=['Under the uninterrupted constant-slope cue, the spatial MCA crossing gives a useful clock for the length-dependent redistribution burden. Its initial distance to the midpoint grows with length, while its mean distribution speed also changes. The endpoint-MCA crossing adds a local clock to this whole-cell comparison.', 'The shortest strength-passing finite pulse is compared with these uninterrupted-cue clocks. The distance and speed panels describe the evolving MCA distribution, and the paired source budgets show how endpoint forcing and Rho-driven cortical flow contribute at different stages of the cue.'],caption=flow_cap,sources=links('debelly_flow')+[['../../data/runs/length_mechanisms_20261008/narrative_20261009/memory/memory_measurements.csv','All twelve lengths and measured clocks']]),flow_budget=dict(image='narrative_20261009/memory/debelly_budget.png',renderer=Path(__file__).name,title='Endpoint and Rho contributions evolve during stimulation',paragraphs=['The source budget resolves the changing attachment distribution into endpoint-driven advection, Rho-driven advection and the remaining native terms. Comparing an early cue state with the state just before withdrawal exposes the changing balance during the same strength-selected pulses.', 'These contributions describe how the saved state changes over one native step. Their sum recovers the total MCA polarity rate; the comparison complements the independently evolved tension-feedback ablation by identifying the immediate transport terms in the intact model.'],caption=budget_cap,sources=links('debelly_budget')+[['../../data/runs/length_mechanisms_20261008/narrative_20261009/memory/cue_timing_budget.csv','Budget values and actual sampling times']])))
 overrides=json.loads((DATA/'chapter_overrides.json').read_text())
 overrides['mca_recovery']=dict(image='narrative_20261009/memory/debelly_recovery.png',renderer=Path(__file__).name,title='The reversed pole rebuilds through the complete basal recovery',paragraphs=['Recovery extends well beyond the rapid cue response. The whole-cell MCA moment and the local endpoint difference follow the inherited state while Rac approaches its reversed strength. The logarithmic time axis resolves the early response and the full 10000-s recovery within the same comparison.', 'Each trajectory belongs to the shortest tested strength-passing pulse shown in the preceding figure. Its final 500 s satisfy the same scalar reversal criterion; the separate native-field stability measurements remain available with the data.'],caption=(FIG/'debelly_recovery_caption.txt').read_text().strip(),sources=links('debelly_recovery')+[['../../data/runs/length_mechanisms_20261008/narrative_20261009/memory/independent_audit.json','Checkpoint, continuity and independently recomputed gates']])
 write(DATA/'chapter_overrides.json',overrides)
 write(RUN/'status.json',dict(producer='finished',complete=True,audit_passed=True,jobs=21,reused=15,generated_at=n.now()));print(FIG,flush=True)
def links(stem):return [[f'narrative_20261009/memory/{stem}.pdf','Vector figure'],[f'narrative_20261009/memory/{stem}.manifest.json','Legend, hashes and reproduction']]
def flow(rows,traces):
 rr=[r for r in rows if r['condition']=='fixed_slope'];xx=[r['length'] for r in rr]
 fig,axs=plt.subplots(1,3,figsize=(12,4.3),layout='constrained')
 for ax,label,title in zip(axs,'ABC',['Cue and attachment clocks','Attachment redistribution distance','Mean distribution speed']):annotate(ax,label,title);ax.set_xlabel('Cell length (µm)');ax.set_xticks([18,23,27,31,35,39]);ax.set_xlim(17,40)
 for key,color,ls,mark,label in [('minimum_strength_cue','#0072B2','-','s','Strength-passing cue'),('mca_zero','#762a83','-','o','Spatial MCA zero'),('endpoint_mca_zero','#009E73','-','^','Endpoint MCA zero'),('rac_zero','#777777','--',None,'Rac zero')]:axs[0].plot(xx,[r[key] for r in rr],color=color,ls=ls,marker=mark,ms=4,label=label)
 axs[0].set_ylabel('Time (s)');axs[0].legend(frameon=False,fontsize=8,loc='upper center')
 for key,color,ls,mark,label in [('initial_halfmass_distance','#009E73','--','^','Initial half-mass distance'),('centroid_distance_to_zero','#762a83','-','o','Initial centroid distance'),('centroid_shift_during_strength_cue','#0072B2','-','s','Shift during selected pulse')]:axs[1].plot(xx,[r[key] for r in rr],color=color,ls=ls,marker=mark,ms=4,label=label)
 axs[1].set_ylabel('Nominal distance / shift (µm)');axs[1].legend(frameon=False,fontsize=8,loc='upper left')
 axs[2].plot(xx,[r['mean_centroid_speed_to_zero']*1000 for r in rr],color='#762a83',marker='o',ms=4);axs[2].set_ylabel('Centroid shift / time (10⁻³ µm/s)')
 cap='Attachment redistribution burden under the constant-slope cue (0.24), at all twelve eligible corrected lengths. A: shortest tested 1%-strength-passing finite cue (blue), compared with first Rac-polarity zero (gray), spatial bound-MCA polarity zero (purple), and left-minus-right endpoint bound-MCA zero (green). All three crossing clocks use the same uninterrupted-cue histories and linear interpolation between recorded 1-s samples. B: initial half-mass and centroid distances from the nominal midpoint, and centroid displacement during the actual selected finite pulse. C: initial MCA centroid distance divided by its uninterrupted spatial-zero time. C quantifies movement of a distribution; it is neither cortical velocity nor tracked material speed. Spatial coordinates are nominal reference coordinates. Every eligible length has a strength-passing tested pulse; native-field tail flags are retained in the table.'
 save(fig,'debelly_flow',cap,dict(clock_protocol='Uninterrupted native cues at all12 lengths; no mixed withdrawal histories.',cue_selection='Canonical scalar1% strength gate',speed='Distribution centroid displacement/time, nominal coordinates'));return cap
def budgets(traces):
 # Keep the unchanged native half-second source identity, sampled from dense
 # uninterrupted histories so duration−60 is consistently cue-on.
 old=read(previous.OLD/'manifest.json');flowroot=ROOT/'data/runs/debelly_corrected_advection_20261008';rows=[];checks=[]
 from types import SimpleNamespace
 for L in LENGTHS:
  a,z,j,metrics=traces['fixed_slope',L];c=j['case'];D=j['duration']
  if L<=21:
   slot=next(i for i,v in enumerate(old['jobs']) if v['case']['length']==L and v['case']['condition']=='fixed_slope' and v['purpose']=='continuous_and_zero_withdrawal');p=previous.OLD/f'workers/{slot:03d}.npz';meta=read(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256'];prefix='pulse_'
  else:p=flowroot/f'L{L}_native.npz';meta=read(flowroot/f'L{L}_summary.json');assert sha(p)==meta['native_sha256'];prefix=''
  register(p)
  with np.load(p) as native:
   t=native[prefix+'time'];x=native['x'];assert np.allclose(np.diff(t),.5,atol=1e-10)
   for label,target in [('A',60.),('B',D-60.)]:
    i=int(np.argmin(abs(t-target)));assert 0<i<len(t) and abs(t[i]-target)<=.25+1e-9 and 0<t[i]<D
    fields={key:native[prefix+'field_'+key][i-1:i+1].copy() for key in n.DB_NAMES}
    tr=SimpleNamespace(x=x,time=t[i-1:i+1],fields=fields);r,_=ad.decompose(tr,1,L,previous.corrected.response_cue(c));other=sum(r['rate_'+key] for key in ['membrane','diffusion','exchange','replenishment'])
    assert abs(r['rate_boundary']+r['rate_rho']+other-r['reorientation_rate'])<1e-11
    r.update(panel=label,length=L,duration=D,requested_time=target,actual_time=float(t[i]),timing_offset=float(t[i]-target),rate_other=other,source_archive=str(p.relative_to(ROOT)));rows.append(r);checks.append(dict(length=L,panel=label,budget_residual=r['budget_residual'],mechanical_residual=r['mechanical_residual'],timing_offset=r['timing_offset']))
 csvwrite(DATA/'cue_timing_budget.csv',rows)
 fig,axs=plt.subplots(1,2,figsize=(11.8,4.4),layout='constrained',sharex=True,sharey=True)
 curves=[('rate_boundary','#D55E00','-','Endpoint advection'),('rate_rho','#0072B2','-','Rho advection'),('rate_other','#999999','-','Other terms'),('reorientation_rate','#222222','--','Net rate')]
 for ax,label,title in zip(axs,'AB',['60 s after cue onset','60 s before cue withdrawal']):
  annotate(ax,label,title);rr=[r for r in rows if r['panel']==label]
  for key,color,ls,name in curves:ax.plot([r['length'] for r in rr],[r[key]*1e4 for r in rr],color=color,ls=ls,lw=1.8,label=name)
  ax.set_xlabel('Cell length (µm)');ax.set_xlim(17,40);ax.set_xticks([18,23,27,31,35,39])
 axs[0].set_ylabel('−dP(MCA)/dt (10⁻⁴ s⁻¹)');fig.legend(handles=[Line2D([],[],color=c,ls=ls,lw=1.8,label=name) for _,c,ls,name in curves],loc='outside lower center',ncol=4,frameon=False)
 cap='Native bound-MCA polarity-rate budget during each shortest tested 1%-strength-passing constant-slope pulse (0.24). A: 60 s after onset; B: 60 s before withdrawal. Positive rates erase the original attachment bias. Orange: endpoint-driven cortical advection; blue: Rho-gradient-driven advection; gray: membrane motion, diffusion, MCA exchange and uniform restoration; black dashed: their net sum. The exact changing-normalization identity uses one native 0.5-s step with frozen full-flow upwind donors. This is an instantaneous source decomposition at saved states, not independently evolved pathway removal. Dense uninterrupted histories match the finite pulse before withdrawal. Actual nearest 0.5-s times and offsets are tabulated (maximum offset '+f'{max(abs(r["timing_offset"]) for r in rows):.6f}'+' s). All 24 sampled states are cue-on.'
 save(fig,'debelly_budget',cap,dict(audit=checks,maximum_budget_residual=max(abs(r['budget_residual']) for r in rows),selection='Shortest strength-passing pulse; distinct from prior orientation-selected budget'));return cap
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('action',choices=['launch','render']);args=parser.parse_args()
 launch() if args.action=='launch' else render()
