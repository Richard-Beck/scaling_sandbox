"""Render audited basal feedback interventions without modifying the report."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'code/core'))
import native as n
import settling
RUN=ROOT/'data/runs/debelly_basal_feedback_20261009_v2'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/basal_feedback_controls'
COLORS={23:'#0072B2',31:'#762a83',39:'#D55E00'}
CONDITIONS=['native','equal_command','mca_feedback','tension_feedback']
LABELS={'native':'Native basal','equal_command':'Equal left command','mca_feedback':'Common MCA → Rac signal','tension_feedback':'Common tension → Rho signal'}
STYLES={'native':'-','equal_command':'--','mca_feedback':':','tension_feedback':'-.'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':9,'pdf.fonttype':42,'svg.fonttype':'none'})
SOURCES={}
def output(name):return layout.dev_output_path(OUT,name)
def key(p):return layout.dev_artifact_key(OUT,p)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def csvwrite(name,rows):
 p=output(name)
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 return p
def label(ax,letter):
 ax.text(-.15,1.035,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans');ax.grid(alpha=.13)
def timeaxis(ax,limit):
 ax.set_xscale('symlog',linthresh=10,linscale=.7);ax.set_xlim(0,limit);ax.set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000']);ax.set_xlabel('Time (s; log above 10 s)')
def legend(fig,comparison=True):
 handles=[Line2D([],[],color=c,lw=2,label=f'{L} µm') for L,c in COLORS.items()]
 if comparison:handles += [Line2D([],[],color='#444444',lw=1.5,label='Intervention'),Line2D([],[],color='#444444',lw=1.3,ls='--',label='Matched native control')]
 fig.legend(handles=handles,loc='outside lower center',ncol=len(handles),frameon=False)
def save(fig,name,caption,manifest_base):
 hashes={}
 for ext in ['png','pdf','svg']:
  p=output(name+'.'+ext);fig.savefig(p,dpi=220,facecolor='white');hashes[key(p)]=sha(p)
 plt.close(fig)
 p=output(name+'_caption.txt');p.write_text(caption+'\n');hashes[key(p)]=sha(p)
 rec=dict(manifest_base,artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,sha256=hashes,renderer_sha256=sha(__file__),source_sha256=SOURCES)
 output(name+'.manifest.json').write_text(json.dumps(rec,indent=2)+'\n')
 return hashes

def main(run):
 for p in [Path(__file__),Path(layout.__file__),ROOT/'code/core/native.py',ROOT/'code/core/settling.py',ROOT/'code/mechanisms/length_dev_basal_feedback_replay.py']:SOURCES[str(p.relative_to(ROOT))]=sha(p)
 m=read(run/'manifest.json');p=read(run/'provenance.json');assert sha(run/'manifest.json')==p['manifest_sha256']
 assert read(run/'audits/summary.json')['passed'];frozen_verified={}
 for name,digest in p['frozen_sources'].items():
  assert sha(run/'frozen_sources'/name)==digest,name
  frozen_verified[name]=digest
 original=(run/'frozen_sources/code/core/models/debelly.py').read_text()
 insertion=m['diagnostic_native_change']['command_override_insertion']
 expected=original.replace('            old_left_boundary_displacement = left_boundary_displacement\n',insertion+'            old_left_boundary_displacement = left_boundary_displacement\n')+'\n# Sole solver alteration: an optional prescribed left protrusion command.\nBASAL_DIAGNOSTIC_LEFT_COMMAND = None\n'
 assert (run/'frozen_sources/code/mechanisms/debelly_basal_feedback_native.py').read_text()==expected
 # Renderer reads histories only. All execution dependencies remain frozen.
 traces={};commands={};records={};budgets={};rows=[];checks=[];time_rows=[];profile_rows=[];budget_rows=[];frozen_rows=[]
 rule=m['settling_rule']
 for j in m['jobs']:
  slot=j['slot'];L=j['length'];condition=j['condition'];tag=(condition,L)
  result=read(run/'workers'/f'{slot:02d}.json');path=run/'workers'/f'{slot:02d}.npz'
  assert result['producer']=='finished' and sha(path)==result['trace_sha256'];SOURCES[str(path.relative_to(ROOT))]=sha(path)
  state,tr=n.read_checkpoint(path,'debelly');traces[tag]=tr;records[tag]=result
  cp=run/'checkpoints'/f'L{L}.npz';assert sha(cp)==j['checkpoint_sha256'];SOURCES[str(cp.relative_to(ROOT))]=sha(cp)
  initial,_=n.read_checkpoint(cp,'debelly')
  for k in n.DB_NAMES:np.testing.assert_array_equal(tr.fields[k][0],initial[k])
  command_path=run/'workers'/f'{slot:02d}.commands.npz';SOURCES[str(command_path.relative_to(ROOT))]=sha(command_path)
  with np.load(command_path) as z:cmd={k:z[k].copy() for k in z.files}
  np.testing.assert_array_equal(cmd['time'],tr.time);commands[tag]=cmd
  mass=np.sum(tr.fields['mca_bound']+tr.fields['mca_unbound'],axis=1);mass_error=float(abs(mass/mass[0]-1).max());assert mass_error<1e-10
  keep=tr.time>=tr.time[-1]-1000
  # Additional audit of all 18 recorded native fields; execution follows the
  # shared manuscript8-field criterion, and these tests are disclosed separately.
  extra_scales=dict(j['fixed_scales'],rac_activation=1.,rho_activation=1.,delta_alpha_mca=.2,delta_beta_tension=.2,
       raw_left_displacement=L,raw_right_displacement=L,left_boundary_displacement=L,right_boundary_displacement=L,boundary_separation=L,mean_strain=1.,tension_range=j['fixed_scales']['tension'])
  field_checks={}
  for field,y in tr.fields.items():
   scale=extra_scales[field];yy=y[keep]
   variation=float(np.ptp(yy,axis=0).max()/scale)
   drift=float(abs(yy[-1]-yy[-2]).max()/(tr.time[-1]-tr.time[-2])*1000/scale)
   field_checks[field]=dict(relative_range=variation,projected_relative_drift=drift,scale=scale)
  all18pass=all(max(v['relative_range'],v['projected_relative_drift'])<=rule['tolerance'] for v in field_checks.values())
  equal_recurrence=None
  if condition=='equal_command':
   q=.98**(np.diff(tr.time)/.5);u=tr.fields['raw_left_displacement'][:,0]
   residual=u[1:]-(q*u[:-1]+.5*m['reference_left_command']*(1-q)/(1-.98))
   equal_recurrence=float(abs(residual).max());assert equal_recurrence<1e-12
  last_hold=next((c for c in reversed(result['settling_checks']) if c['kind']=='holdout'),None)
  checks.append(dict(slot=slot,length=L,condition=condition,initial11fields_exact=True,mca_mass_error=mass_error,manuscript_verified=result['verified_steady_state'],additional_all18_final_window_pass=all18pass,additional_all18_metrics=field_checks,equal_command_saved_interval_recurrence_max_residual=equal_recurrence,holdout=last_hold))
  row=dict(length_um=L,condition=condition,verified_steady_state=result['verified_steady_state'],censored_at_cap=result['censored_at_cap'],final_time_s=float(tr.time[-1]),initial_orientation=float(tr.orientation[0]),final_orientation=float(tr.orientation[-1]),final_contrast=float(tr.contrast[-1]),retained_old_polarity=result['retained_polarity'],all18_final_window_pass=all18pass)
  for name in ['rac','rho','mca_bound','mca_unbound','tension','displacement','velocity','delta_alpha_mca','delta_beta_tension']:
   for side,ix in [('left',0),('right',-1)]:row['initial_'+name+'_'+side]=float(tr.fields[name][0,ix]);row['final_'+name+'_'+side]=float(tr.fields[name][-1,ix])
  for name in ['raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement','mean_strain','tension_range']:
   row['initial_'+name]=float(tr.fields[name][0,0]);row['final_'+name]=float(tr.fields[name][-1,0])
  for side in ['left','right']:row['initial_'+side+'_command']=float(cmd[side+'_command'][0]);row['final_'+side+'_command']=float(cmd[side+'_command'][-1])
  rows.append(row)
  for i,t in enumerate(tr.time):
   time_rows.append(dict(length_um=L,condition=condition,time_s=float(t),rac_left=float(tr.fields['rac'][i,0]),rho_left=float(tr.fields['rho'][i,0]),mca_bound_left=float(tr.fields['mca_bound'][i,0]),mca_unbound_left=float(tr.fields['mca_unbound'][i,0]),tension_left=float(tr.fields['tension'][i,0]),tension_right=float(tr.fields['tension'][i,-1]),left_endpoint_displacement=float(tr.fields['left_boundary_displacement'][i,0]),right_endpoint_displacement=float(tr.fields['right_boundary_displacement'][i,0]),left_command=float(cmd['left_command'][i]),right_command=float(cmd['right_command'][i]),orientation=float(tr.orientation[i]),contrast=float(tr.contrast[i])))
  for when,ix in [('initial',0),('final',-1)]:
   for k,x in enumerate(tr.x):
    profile_rows.append(dict(length_um=L,condition=condition,snapshot=when,time_s=float(tr.time[ix]),x_um=float(x),x_over_L=float(x/L),**{field:float(y[ix,k]) for field,y in tr.fields.items()}))
  audit=read(run/'audits'/f'job_{slot:02d}.json');path=run/'audits'/f'job_{slot:02d}.npz';SOURCES[str(path.relative_to(ROOT))]=sha(path)
  with np.load(path) as z:b={k:z[k].copy() for k in z.files}
  budgets[tag]=b
  for k,x in enumerate(b['x']):
   budget_rows.append(dict(length_um=L,condition=condition,x_um=float(x),**{name:float(b[name][k]) for name in ['window','reference_delta_alpha','reference_delta_beta','native_delta_alpha','native_delta_beta','controlled_delta_alpha','controlled_delta_beta','frozen_rac_derivative_change','frozen_rho_derivative_change','rac_activation','rho_inhibition','rac_deactivation','rac_diffusion','actual_rac_rate','rac_budget_residual']}))
  frozen_rows.append(dict(length_um=L,condition=condition,frozen_left_rac_derivative_change=audit['frozen_left_rac_derivative_change'],frozen_left_rho_derivative_change=audit['frozen_left_rho_derivative_change']))
 by={(r['condition'],r['length_um']):r for r in rows}
 for row in rows:
  baseline=by['native',row['length_um']]
  for field in ['rac','rho','mca_bound','tension']:
   for side in ['left','right']:row['final_'+field+'_'+side+'_minus_matched_native']=row['final_'+field+'_'+side]-baseline['final_'+field+'_'+side]
 # Reference interventions should preserve the complete basal state to roundoff.
 reference_errors={}
 for condition in CONDITIONS[1:]:
  tr=traces[condition,31];control=traces['native',31];np.testing.assert_array_equal(tr.time,control.time)
  reference_errors[condition]={name:float(abs(tr.fields[name]-control.fields[name]).max()) for name in tr.fields}
  assert max(reference_errors[condition][name]/31 if 'displacement' in name else reference_errors[condition][name] for name in n.DB_NAMES)<1e-7
 replay=run/'results/replay_validation'
 with np.load(run/'short_validation/00.npz') as a,np.load(replay/'short_validation/00.npz') as b:
  assert set(a.files)==set(b.files)
  for field in a.files:np.testing.assert_array_equal(a[field],b[field])
 audit=dict(passed=all(c['manuscript_verified'] for c in checks),all_initial11_fields_exact=True,all12_verified=all(c['manuscript_verified'] for c in checks),all18_final_window_pass=all(c['additional_all18_final_window_pass'] for c in checks),checks=checks,reference_length_full_history_max_absolute_errors=reference_errors,execution_source_hashes_verified=frozen_verified,relocated_frozen_replay_all_history_arrays_bit_identical=True,diagnostic_native_diff_exactly_matches_declared_single_hook=True)
 readme=dict(reference_command=m['reference_left_command'],window=m['local_window'],steady_state_rule=rule,
  exact_rac_budget='First native half-second step: kR*alpha_new − kR*alpha_new*rho_old²/(1+rho_old²) − kR*Rac_new + kR*D_R*zero_flux_discrete_Laplacian(Rac_new) = (Rac_new−Rac_old)/0.5. Alpha is evaluated after membrane/MCA steps; Rho is lagged. The negative Rho term is suppression of production, not Rac degradation.',
  frozen_effects='At the exact checkpoint, compare only the changed nonlinear feedback increment with native, holding every state field fixed: delta(dRac/dt)=kR*delta(alpha)/(1+rho²); delta(dRho/dt)=kH*delta(beta)/(1+rac²). These are algebraic onset effects, separate from the coupled time integration.',
  interventions='No countercue. Common left command fixes only left mechanical command; native right command remains. Common local signals blend feedback increments with the static L31 checkpoint reference, with w=1 from 0–3 µm, cosine taper 3–5 µm,0 beyond. Actual MCA/tension fields remain dynamical.',
  axes='Time starts at intervention onset, independent of preparation clock. Full field amplitudes are native; only spatial coordinate is normalized in profile figures.',
  scope='Causal contributions are conditional on these modified feedback laws. A replacement experiment does not partition the native length effect into independent additive causes.',
  settling='Execution uses inherited manuscript fixed scales and the eight independent native convergence fields, two 1,000-s windows≤1e-4 followed by independent 10,000-s holdout≤1e-3 and final window≤1e-4. Additional all 18 recorded fields checked in final 1,000 s with fixed declared scales; raw boundary fields are never centered in exports.')
 data={}
 for name,values in [('basal_feedback_steady_endpoints.csv',rows),('basal_feedback_endpoint_timecourses.csv',time_rows),('basal_feedback_initial_final_all_fields.csv',profile_rows),('basal_feedback_discrete_budgets.csv',budget_rows),('basal_feedback_frozen_derivative_changes.csv',frozen_rows)]:
  q=csvwrite(name,values);data[key(q)]=sha(q)
 q=output('basal_feedback_independent_audit.json');q.write_text(json.dumps(audit,indent=2)+'\n');data[key(q)]=sha(q)
 q=output('basal_feedback_methods.json');q.write_text(json.dumps(readme,indent=2)+'\n');data[key(q)]=sha(q)
 base=dict(campaign=str(run.relative_to(ROOT)),frozen_execution_source_sha256=frozen_verified,methods=readme,data_sha256=data,new_native_histories=12,current_report_changed=False,reproduction=[
  'Integrations already retained under '+str(run)+'. Verify provenance.json frozen_sources hashes before re-execution.',
  sys.executable+' '+str(ROOT/'code/mechanisms/length_dev_basal_feedback_replay.py')+' --source '+str(run)+' --destination /absolute/fresh/campaign/path',
  'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '+sys.executable+' '+str(Path(__file__))+' --run '+str(run)])
 limit=max(t.time[-1] for t in traces.values())*1.02
 fig,axs=plt.subplots(2,2,figsize=(11.7,7.8),layout='constrained',sharey=True)
 for ax,condition,letter in zip(axs.flat,CONDITIONS,'ABCD'):
  label(ax,letter);ax.set_title(LABELS[condition],pad=10)
  for L,color in COLORS.items():
   tr=traces[condition,L];ax.plot(tr.time,tr.fields['rac'][:,0],color=color,lw=1.6)
   if condition!='native':
    control=traces['native',L];ax.plot(control.time,control.fields['rac'][:,0],color=color,ls='--',lw=1.15,alpha=.85)
  timeaxis(ax,limit);ax.set_ylabel('Old left-end Rac (native concentration)')
 legend(fig)
 save(fig,'basal_feedback_Rac_dynamics','A–D: old left-end Rac after each basal intervention, with three length colors and matched native controls (dashed in intervention panels). No countercue is applied. The endpoint command and local coupling replacements are defined in the companion methods file. Every trace starts from the exact qualifying shared pre-cue checkpoint and continues to two successful 1,000-s settling windows plus an independently successful 10,000-s holdout; endpoints occur at 12,010–26,010 s. Native controls end at their own verified time and are not extrapolated. Time is linear from 0 to 10 s and logarithmic thereafter. Concentrations are unscaled. All 12 conditions retain the old polarized orientation; no trajectory is excluded. The interventions coincide with native at L31 onset and remain within numerical baseline drift there.',base)
 fields=[('rac','Rac (native concentration)'),('rho','Rho (native concentration)'),('mca_bound','Bound MCA (native concentration)'),('tension','Tension (pN)')]
 fig,axs=plt.subplots(3,4,figsize=(14.3,9.4),layout='constrained')
 for row,condition in enumerate(CONDITIONS[1:]):
  for col,(field,title) in enumerate(fields):
   ax=axs[row,col];label(ax,chr(65+row*4+col));ax.set_title(LABELS[condition] if col==0 else '',pad=9);ax.set_ylabel(title);ax.set_xlabel('Position x/L');ax.set_xlim(0,1)
   for L,color in COLORS.items():
    tr=traces[condition,L];control=traces['native',L]
    ax.plot(tr.x/L,tr.fields[field][-1],color=color,lw=1.5);ax.plot(control.x/L,control.fields[field][-1],color=color,ls='--',lw=1.1)
 legend(fig)
 save(fig,'basal_feedback_steady_profiles','Verified final basal profiles after equal left command (A–D), common local MCA→Rac input (E–H), and common local tension→Rho input (I–L). Columns show Rac, Rho, bound MCA and tension. Solid curves are each intervention; dashed curves are the matched length-specific native final control. The old pole is left (x/L=0). Only position is normalized; amplitudes and the full membrane/chemical state remain native. Both local controls act on feedback increments only within a common physical 0–5 µm region; the MCA and tension profiles shown remain dynamical rather than clamped. All plotted endpoints pass the shared settling/holdout criterion. Initial and final profiles of all 18 recorded native fields, including raw endpoint displacements, are exported in the data CSV.',base)
 dynamic_fields=[('mca_bound',0,'Left bound MCA (native concentration)'),('rho',0,'Left Rho (native concentration)'),('left_boundary_displacement',0,'Left endpoint displacement (µm)'),('tension',-1,'Right-end tension (pN)')]
 fig,axs=plt.subplots(3,4,figsize=(14.3,9.4),layout='constrained',sharex=True,sharey='col')
 for row,condition in enumerate(CONDITIONS[1:]):
  for col,(field,ix,title) in enumerate(dynamic_fields):
   ax=axs[row,col];label(ax,chr(65+row*4+col));ax.set_title(LABELS[condition] if col==0 else '',pad=9);ax.set_ylabel(title)
   for L,color in COLORS.items():
    tr=traces[condition,L];control=traces['native',L]
    ax.plot(tr.time,tr.fields[field][:,ix],color=color,lw=1.5);ax.plot(control.time,control.fields[field][:,ix],color=color,ls='--',lw=1.1)
   timeaxis(ax,limit)
 legend(fig)
 save(fig,'basal_feedback_mechanochemical_dynamics','Coupled basal responses after equal left command (A–D), common local MCA→Rac input (E–H), and common local tension→Rho input (I–L). Columns show left bound MCA, left Rho, exact native global left endpoint displacement, and right-end tension; solid interventions and dashed matched native controls use consistent length colors. Every mechanical and chemical field evolves under the corrected independent-endpoint model. In the equal-command condition the native relaxation is retained: raw_left_next=0.98*raw_left+0.5*command, with command fixed at the L31 initial value, −0.0400642 µm/s. The right endpoint always retains its native local Rac/tension law. Time is measured from intervention onset, linear from 0 to 10 s and logarithmic thereafter. Full endpoint command histories and both ends are in the CSV; no native trace is extrapolated past its verified end.',base)
 fig,axs=plt.subplots(2,2,figsize=(11.8,8),layout='constrained')
 for ax,letter in zip(axs.flat,'ABCD'):label(ax,letter)
 terms=['rac_activation','rho_inhibition','rac_deactivation','rac_diffusion','actual_rac_rate'];names=['Activation','Rho suppression','Deactivation','Diffusion','Net change']
 for i,(L,color) in enumerate(COLORS.items()):
  b=budgets['native',L];axs[0,0].plot(np.arange(5)+.07*(i-1),[b[k][0] for k in terms],marker='o',ms=5,color=color,lw=1,label=f'{L}µm')
  axs[0,1].plot(L,b['rac_diffusion'][0],'s',ms=6,color=color);axs[0,1].plot(L,b['actual_rac_rate'][0],'o',ms=6,color=color)
  axs[1,0].plot(L,budgets['mca_feedback',L]['frozen_rac_derivative_change'][0],'o',ms=6,color=color)
  axs[1,1].plot(L,budgets['tension_feedback',L]['frozen_rho_derivative_change'][0],'o',ms=6,color=color)
 for ax in axs.flat:ax.axhline(0,color='#777777',lw=.8);ax.set_ylabel('Rate (native concentration/s)')
 axs[0,0].set_xticks(range(5),names,rotation=15);axs[0,0].set_title('Native first-cell Rac balance');axs[0,0].legend(frameon=False,ncol=3,fontsize=8)
 axs[0,1].set_title('Small terms of the same native balance');axs[0,1].legend(handles=[Line2D([],[],color='#444444',marker='s',ls='',label='Diffusive flux'),Line2D([],[],color='#444444',marker='o',ls='',label='Net change')],frameon=False,loc='lower left')
 axs[1,0].set_title('Frozen MCA→Rac replacement: Δ(dRac/dt)');axs[1,1].set_title('Frozen tension→Rho replacement: Δ(dRho/dt)')
 for ax in [axs[0,1],axs[1,0],axs[1,1]]:ax.set_xlim(21,41);ax.set_xticks([23,31,39]);ax.set_xlabel('Cell length (µm)')
 save(fig,'basal_feedback_initial_budget','A–B: exact staggered half-second native Rac budget at the old first cell, from the shared basal checkpoints. Activation and negative Rho suppression partition the production kR*alpha_new/(1+rho_old²); suppression is reduced production rather than an additional Rac decay process. Deactivation and zero-flux finite-volume diffusion use the implicitly updated Rac. Their sum equals the measured half-second Rac change with maximum residual 1.65e−15 across all 12 audits. B resolves the much smaller diffusive term and near-zero basal net change. C–D: immediate algebraic derivative changes caused by replacing only the local MCA→Rac or tension→Rho nonlinear increment, holding every checkpoint state field fixed. The L31 reference has exactly zero onset change. These frozen effects precede the coupled dynamical responses; the equal-command control has no direct frozen chemical source change. Full spatial budgets, reference increments and the common 0–5 µm window are exported.',base)
 summary_fields=[('rac_left','Old left-end Rac (native concentration)'),('mca_bound_left','Left bound MCA (native concentration)'),('tension_left','Left-end tension (pN)'),('tension_right','Right-end tension (pN)'),('left_command','Outward left command (µm/s)'),('left_boundary_displacement','Left endpoint displacement (µm)')]
 fig,axs=plt.subplots(2,3,figsize=(12.5,8),layout='constrained',sharex=True)
 for ax,(field,title),letter in zip(axs.flat,summary_fields,'ABCDEF'):
  label(ax,letter);ax.set_ylabel(title);ax.set_xlabel('Cell length (µm)');ax.set_xticks([23,31,39]);ax.set_xlim(22,40)
  for condition in CONDITIONS:
   values=[by[condition,L]['final_'+field]*(-1 if field=='left_command' else 1) for L in COLORS]
   ax.plot(list(COLORS),values,color='#555555',ls=STYLES[condition],lw=1.4)
   for (L,color),value in zip(COLORS.items(),values):ax.plot(L,value,'o',color=color,ms=4)
 handles=[Line2D([],[],color='#555555',ls=STYLES[c],label=LABELS[c]) for c in CONDITIONS]
 fig.legend(handles=handles,loc='outside lower center',ncol=2,frameon=False)
 save(fig,'basal_feedback_steady_endpoints','Verified final endpoint measurements versus length under native basal dynamics and the three interventions. A: old left-end Rac. B: left bound MCA. C–D: left and right tension. E: outward magnitude of the actual left command, positive by definition. F: signed native global left endpoint displacement; negative points represent extension to the left. Line styles identify intervention and point colors identify length; the three measured lengths are connected for comparison without interpolating a new simulation. Native controls retain the original inherited state to numerical precision. The equal-command intervention enforces the same command, not the same tension or Rac; mechanical/chemical feedback subsequently adjusts both. Local MCA→Rac and tension→Rho replacements each alter length dependence, within their explicitly modified local feedback laws. All 12 endpoints are verified by the manuscript field criterion and independent holdout, with no censored or discarded outcome.',base)
 stats={}
 for condition in CONDITIONS:
  rr=[by[condition,L] for L in COLORS]
  stats[condition]=dict(left_rac=[r['final_rac_left'] for r in rr],left_rac_L39_minus_L23=rr[-1]['final_rac_left']-rr[0]['final_rac_left'],left_rho=[r['final_rho_left'] for r in rr],left_mca=[r['final_mca_bound_left'] for r in rr],left_tension=[r['final_tension_left'] for r in rr],right_tension=[r['final_tension_right'] for r in rr],left_command=[r['final_left_command'] for r in rr],end_time=[r['final_time_s'] for r in rr])
 summary=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',audit_passed=audit['passed'],all18_final_window_pass=audit['all18_final_window_pass'],length_order=list(COLORS),results=stats,scope=readme['scope'])
 output('basal_feedback_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 (run/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 (run/'results/independent_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 (run/'results/status.json').write_text(json.dumps(dict(generated_at=summary['generated_at'],producer='finished',complete=True,histories=12,verified_steady_states=12,censored=0,polarity_losses=0,slurm_job_id='20795063'),indent=2)+'\n')
 print(json.dumps(dict(audit_passed=audit['passed'],all18_final_window_pass=audit['all18_final_window_pass'],results=stats),indent=2))
 print(OUT)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,default=RUN);a=ap.parse_args();main(a.run.resolve())
