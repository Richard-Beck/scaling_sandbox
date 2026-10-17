"""Audited single global tension→Rho ablation; no report or solver changes."""
from pathlib import Path
from dataclasses import replace,asdict
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
RUN=ROOT/'data/runs/debelly_global_tension_rho_ablation_20261009'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/tension_rho_ablation'
COLORS={23:'#0072B2',31:'#762a83',39:'#D55E00'}
SOURCES={};DATA={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':9,'pdf.fonttype':42,'svg.fonttype':'none'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def output(name):return layout.dev_output_path(OUT,name)
def key(p):return layout.dev_artifact_key(OUT,p)
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def csvwrite(name,rows):
 p=output(name)
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 DATA[key(p)]=sha(p)
def label(ax,letter):ax.text(-.16,1.03,letter,transform=ax.transAxes,fontfamily='DejaVu Sans',fontsize=12,fontweight='bold');ax.grid(alpha=.13)
def legend(fig):
 handles=[Line2D([],[],color=color,lw=2,label=f'{L} µm') for L,color in COLORS.items()]+[Line2D([],[],color='#444444',lw=1.5,label='Tension→Rho removed'),Line2D([],[],color='#444444',ls='--',lw=1.3,label='Native basal')]
 fig.legend(handles=handles,loc='outside lower center',ncol=5,frameon=False)
def save(fig,name,caption,base):
 hashes={}
 for ext in ['png','pdf','svg']:
  p=output(name+'.'+ext);fig.savefig(p,dpi=220,facecolor='white');hashes[key(p)]=sha(p)
 plt.close(fig);p=output(name+'_caption.txt');p.write_text(caption+'\n');hashes[key(p)]=sha(p)
 output(name+'.manifest.json').write_text(json.dumps(dict(base,artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,sha256=hashes,data_sha256=DATA,source_sha256=SOURCES,renderer_sha256=sha(__file__)),indent=2)+'\n')

def main(run):
 for p in [Path(__file__),Path(layout.__file__),ROOT/'code/core/native.py',ROOT/'code/core/settling.py',ROOT/'code/mechanisms/length_dev_global_tension_rho_ablation_replay.py']:SOURCES[str(p.relative_to(ROOT))]=sha(p)
 m=read(run/'manifest.json');provenance=read(run/'provenance.json');assert sha(run/'manifest.json')==provenance['manifest_sha256'];assert read(run/'audits/summary.json')['passed']
 for name,digest in provenance['frozen_sources'].items():assert sha(run/'frozen_sources'/name)==digest,name
 native_source=(run/'frozen_sources/code/core/models/debelly.py').read_text()
 assert (run/'frozen_sources/code/mechanisms/debelly_basal_feedback_native.py').read_text()==native_source+'\nBASAL_DIAGNOSTIC_LEFT_COMMAND = None\n'
 p=replace(n.db._parameters(),center_boundary_motion=False);R,H=n.db._homogeneous_rac_rho(p,.5);threshold=1.3*R
 assert threshold==m['original_rac_threshold'];traces={};commands={};results={};rows=[];time_rows=[];profiles=[];budgets=[];rho_budgets=[];audits=[]
 for job,control in zip(m['jobs'],m['native_controls_reused']):
  L=job['length'];slot=job['slot'];cp=run/'checkpoints'/f'L{L}.npz';assert sha(cp)==job['checkpoint_sha256'];SOURCES[str(cp.relative_to(ROOT))]=sha(cp);state,initial=n.read_checkpoint(cp,'debelly')
  for condition,path,result_path,command_path,digest in [
   ('ablation',run/'workers'/f'{slot:02d}.npz',run/'workers'/f'{slot:02d}.json',run/'workers'/f'{slot:02d}.commands.npz',None),
   ('native',ROOT/control['trace'],ROOT/control['result'],ROOT/control['commands'],control['trace_sha256'])]:
   record=read(result_path);assert record['producer']=='finished' and record['initial_state_exact']
   if condition=='native':assert sha(result_path)==control['result_sha256'] and sha(command_path)==control['commands_sha256']
   assert sha(path)==record['trace_sha256'] and (digest is None or sha(path)==digest)
   SOURCES[str(path.relative_to(ROOT))]=sha(path);SOURCES[str(command_path.relative_to(ROOT))]=sha(command_path)
   final,tr=n.read_checkpoint(path,'debelly');traces[condition,L]=tr;results[condition,L]=record
   with np.load(command_path) as z:cmd={k:z[k].copy() for k in z.files}
   np.testing.assert_array_equal(cmd['time'],tr.time);commands[condition,L]=cmd
   assert record['parameters']==asdict(p)
   for field in n.DB_NAMES:np.testing.assert_array_equal(tr.fields[field][0],state[field])
   # Sixteen recorded initial fields are identical; the two derived Rho-rate
   # fields change immediately because the intervention changes their algebra.
   for field in tr.fields:
    if condition=='ablation' and field in ['rho_activation','delta_beta_tension']:continue
    np.testing.assert_array_equal(tr.fields[field][0],initial.fields[field][-1])
   native_da,native_db=n.db._feedback_terms(p,tr.fields['mca_bound'],tr.fields['tension'])
   np.testing.assert_array_equal(tr.fields['delta_alpha_mca'],native_da)
   if condition=='ablation':
    np.testing.assert_array_equal(tr.fields['delta_beta_tension'],np.zeros_like(native_db))
    np.testing.assert_array_equal(tr.fields['rho_activation'],np.ones_like(native_db))
   else:np.testing.assert_array_equal(tr.fields['delta_beta_tension'],native_db)
   expected_left=p.protrusion_speed*n.db.smooth_threshold(tr.fields['rac'][:,0],p.protrusion_switch_scale,threshold)/(1+tr.fields['tension'][:,0]**2/p.protrusion_tension_scale)
   expected_right=-p.protrusion_speed*n.db.smooth_threshold(tr.fields['rac'][:,-1],p.protrusion_switch_scale,threshold)/(1+tr.fields['tension'][:,-1]**2/p.protrusion_tension_scale)
   np.testing.assert_array_equal(cmd['left_command'],expected_left);np.testing.assert_array_equal(cmd['right_command'],expected_right)
   mass=(tr.fields['mca_bound']+tr.fields['mca_unbound']).sum(axis=1);mass_error=float(abs(mass/mass[0]-1).max());assert mass_error<1e-10
   assert tr.time[1]==.5
   oldrho=tr.fields['rho'][0];newrho=tr.fields['rho'][1];newrac=tr.fields['rac'][1];beta=tr.fields['rho_activation'][1];dx=L/job['cells']
   activation=p.rho_deactivation*beta;suppression=-activation*newrac**2/(1+newrac**2);decay=-p.rho_deactivation*newrho
   diffusion=p.rho_deactivation*p.rho_diffusion*np.diff(np.r_[0.,np.diff(newrho)/dx,0.])/dx
   rate=(newrho-oldrho)/.5;residual=rate-activation-suppression-decay-diffusion;rho_error=float(abs(residual).max());assert rho_error<1e-13
   for i,x in enumerate(tr.x):rho_budgets.append(dict(length_um=L,condition=condition,x_um=float(x),activation=float(activation[i]),rac_suppression=float(suppression[i]),rho_deactivation=float(decay[i]),rho_diffusion=float(diffusion[i]),actual_rho_rate=float(rate[i]),budget_residual=float(residual[i])))
   scales=dict(job['fixed_scales'],rac_activation=1.,rho_activation=1.,delta_alpha_mca=.2,delta_beta_tension=.2,raw_left_displacement=L,raw_right_displacement=L,left_boundary_displacement=L,right_boundary_displacement=L,boundary_separation=L,mean_strain=1.,tension_range=job['fixed_scales']['tension'])
   keep=tr.time>=tr.time[-1]-1000;field_metrics={}
   for field,y in tr.fields.items():
    yy=y[keep];field_metrics[field]=dict(scale=scales[field],relative_range=float(np.ptp(yy,axis=0).max()/scales[field]),projected_relative_drift=float(abs(yy[-1]-yy[-2]).max()/(tr.time[-1]-tr.time[-2])*1000/scales[field]))
   all18pass=all(max(v['relative_range'],v['projected_relative_drift'])<=1e-4 for v in field_metrics.values())
   hold=next(c for c in reversed(record['settling_checks']) if c['kind']=='holdout')
   audits.append(dict(length_um=L,condition=condition,initial11_restart_fields_exact=True,other16_recorded_initial_fields_exact=True,only_derived_rho_rate_fields_changed_at_onset=condition=='ablation',parameters_original=True,threshold_original=True,direct_tension_inhibition_retained=True,mca_feedback_native=True,mca_mass_error=mass_error,rho_discrete_budget_max_residual=rho_error,manuscript_verified=record['verified_steady_state'],additional_all18_final_window_pass=all18pass,additional_all18_final_window_metrics=field_metrics,holdout=hold))
   row=dict(length_um=L,condition=condition,verified_steady_state=record['verified_steady_state'],censored_at_cap=record['censored_at_cap'],final_time_s=float(tr.time[-1]),rac_polarity=float(tr.orientation[-1]),rac_contrast=float(tr.contrast[-1]),old_orientation_gate_pass=bool(tr.orientation[-1]<-.05 and tr.contrast[-1]>=.05),left_biased_nonuniform_state=bool(tr.orientation[-1]<0 and tr.contrast[-1]>=.05),rac_threshold=threshold,left_command=float(cmd['left_command'][-1]),outward_left_command=-float(cmd['left_command'][-1]),right_command=float(cmd['right_command'][-1]))
   for field in ['rac','rho','mca_bound','mca_unbound','tension','velocity','displacement']:
    for side,ix in [('left',0),('right',-1)]:row[field+'_'+side]=float(tr.fields[field][-1,ix])
   for field in ['left_boundary_displacement','right_boundary_displacement','raw_left_displacement','raw_right_displacement','mean_strain','tension_range']:row[field]=float(tr.fields[field][-1,0])
   rows.append(row)
   for i,t in enumerate(tr.time):
    time_rows.append(dict(length_um=L,condition=condition,time_s=float(t),rac_polarity=float(tr.orientation[i]),rac_contrast=float(tr.contrast[i]),rac_left=float(tr.fields['rac'][i,0]),rho_left=float(tr.fields['rho'][i,0]),mca_bound_left=float(tr.fields['mca_bound'][i,0]),tension_left=float(tr.fields['tension'][i,0]),tension_right=float(tr.fields['tension'][i,-1]),left_endpoint_displacement=float(tr.fields['left_boundary_displacement'][i,0]),right_endpoint_displacement=float(tr.fields['right_boundary_displacement'][i,0]),left_command=float(cmd['left_command'][i]),right_command=float(cmd['right_command'][i])))
   for snapshot,ix in [('initial',0),('final',-1)]:
    for i,x in enumerate(tr.x):profiles.append(dict(length_um=L,condition=condition,snapshot=snapshot,time_s=float(tr.time[ix]),x_um=float(x),x_over_L=float(x/L),**{field:float(y[ix,i]) for field,y in tr.fields.items()}))
  auditpath=run/'audits'/f'job_{slot:02d}.npz';SOURCES[str(auditpath.relative_to(ROOT))]=sha(auditpath)
  with np.load(auditpath) as z:
   for i,x in enumerate(z['x']):budgets.append(dict(length_um=L,x_um=float(x),**{field:float(z[field][i]) for field in z.files if field!='x'}))
 by={(r['condition'],r['length_um']):r for r in rows}
 native_u=np.array([by['native',L]['left_boundary_displacement'] for L in COLORS]);ablated_u=np.array([by['ablation',L]['left_boundary_displacement'] for L in COLORS])
 summary=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',length_order=list(COLORS),native_signed_left_displacements_um=native_u.tolist(),ablated_signed_left_displacements_um=ablated_u.tolist(),native_displacement_spread_um=float(np.ptp(native_u)),ablated_displacement_spread_um=float(np.ptp(ablated_u)),displacement_spread_reduction_fraction=float(1-np.ptp(ablated_u)/np.ptp(native_u)),all_three_ablations_verified=all(results['ablation',L]['verified_steady_state'] for L in COLORS),all_three_left_biased_nonuniform=all(by['ablation',L]['left_biased_nonuniform_state'] for L in COLORS),ablated_orientation_gate_passes=[L for L in COLORS if by['ablation',L]['old_orientation_gate_pass']],results=rows)
 with np.load(run/'short_validation/00.npz') as a,np.load(run/'results/replay_validation/short_validation/00.npz') as b:
  assert set(a.files)==set(b.files)
  for field in a.files:np.testing.assert_array_equal(a[field],b[field])
 audit=dict(passed=all(a['manuscript_verified'] and a['additional_all18_final_window_pass'] for a in audits),native_solver_source_unchanged=True,global_zero_delta_beta_exact=True,original_parameters_and_homogeneous_reference_unchanged=True,original_homogeneous_rac=R,original_homogeneous_rho=H,original_rac_threshold=threshold,all_fields_initial_identity_except_required_derived_rho_rates=True,checks=audits,frozen_execution_source_sha256=provenance['frozen_sources'],relocated_frozen_replay_all_history_arrays_bit_identical=True)
 csvwrite('global_tension_rho_ablation_endpoints.csv',rows);csvwrite('global_tension_rho_ablation_endpoint_timecourses.csv',time_rows);csvwrite('global_tension_rho_ablation_initial_final_all_fields.csv',profiles);csvwrite('global_tension_rho_ablation_discrete_budgets.csv',budgets);csvwrite('global_tension_rho_ablation_rho_budgets.csv',rho_budgets)
 for name,value in [('summary',summary),('independent_audit',audit)]:
  q=output('global_tension_rho_ablation_'+name+'.json');q.write_text(json.dumps(value,indent=2)+'\n');DATA[key(q)]=sha(q);(run/'results'/f'{name}.json').write_text(json.dumps(value,indent=2)+'\n')
 methods=dict(intervention=m['intervention'],retained=m['retained'],parameters=asdict(p),original_homogeneous_rac=R,original_homogeneous_rho=H,original_protrusion_threshold=threshold,settling=m['settling_rule'],additional_audit='All 18 recorded fields also tested in the final 1,000 s, using fixed scales explicitly stored in independent_audit.json. Original 11 restart fields exact; only the two derived Rho-rate fields necessarily change at onset.',budget='Exact staggered Rac budget: activation kR*alpha_new, Rho suppression −kR*alpha_new*rho_old²/(1+rho_old²), decay −kR*Rac_new, zero-flux implicit finite-volume diffusion kR*D_R*laplacian(Rac_new). Their sum equals (Rac_new−Rac_old)/0.5. Rho uses beta_new and Rac_new: kH*beta_new − kH*beta_new*Rac_new²/(1+Rac_new²) − kH*Rho_new + kH*D_H*laplacian(Rho_new) = (Rho_new−Rho_old)/0.5; the Rho budget CSV retains every term for ablation and native controls. The frozen-state intervention effect is delta(dRho/dt)=−kH*delta_beta_native/(1+Rac_old²).',polarity='Rac polarity is its normalized first spatial moment; negative means old left bias. Spatial contrast is max(Rac)−min(Rac). The separate numerical orientation gate requires polarity<−0.05 and contrast>=0.05; falling slightly below that moment threshold does not imply a spatially uniform Rac profile.',commands='Instantaneous native endpoint law evaluated at each recorded state; signed left command=−0.15*S(Rac_left)/(1+T_left²/10). Direct tension inhibition remains. Native half-second relaxation u_next=.98*u+.5*command is unchanged; actual global boundary displacements are used, never first interior-cell displacement.',scope='One global feedback ablation at three lengths. This tests the modified model; it does not uniquely partition native causal contributions.')
 q=output('global_tension_rho_ablation_methods.json');q.write_text(json.dumps(methods,indent=2)+'\n');DATA[key(q)]=sha(q)
 base=dict(campaign=str(run.relative_to(ROOT)),slurm_job='20799436',new_native_histories=3,native_controls_reused=3,current_report_changed=False,methods=methods,frozen_execution_source_sha256=provenance['frozen_sources'],reproduction=[
  sys.executable+' '+str(ROOT/'code/mechanisms/length_dev_global_tension_rho_ablation_replay.py')+' --source '+str(run)+' --destination /absolute/fresh/campaign/path',
  'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 '+sys.executable+' '+str(Path(__file__))+' --run '+str(run)])
 fields=[('left_boundary_displacement','Signed left endpoint displacement (µm)'),('outward_left_command','Outward left command (µm/s)'),('rac_polarity','Rac polarity (negative: left)'),('rac_contrast','Rac spatial contrast (native conc.)'),('rac_left','Left Rac (native concentration)'),('rho_left','Left Rho (native concentration)'),('mca_bound_left','Left bound MCA (native concentration)'),('tension_left','Left-end tension (pN)'),('tension_right','Right-end tension (pN)')]
 fig,axs=plt.subplots(3,3,figsize=(12.5,10.5),layout='constrained',sharex=True)
 for ax,(field,title),letter in zip(axs.flat,fields,'ABCDEFGHI'):
  label(ax,letter);ax.set_ylabel(title);ax.set_xlabel('Cell length (µm)');ax.set_xlim(22,40);ax.set_xticks([23,31,39])
  for condition,style in [('native','--'),('ablation','-')]:
   yy=[by[condition,L][field] for L in COLORS];ax.plot(list(COLORS),yy,color='#555555',ls=style,lw=1.5)
   for (L,color),y in zip(COLORS.items(),yy):ax.plot(L,y,'o',color=color,ms=5)
  if field=='tension_left':ax.axhline(0,color='#888888',lw=.7)
 legend(fig)
 caption='Verified basal endpoints after removing the tension-dependent increment to Rho activation globally, at the same corrected L = 23, 31, 39 µm pre-cue checkpoints. Solid: ablation; dashed: exact reused native controls; colors identify length. A–B show signed global left endpoint displacement and instantaneous outward command, with native direct tension inhibition and the 0.98-per-half-second relaxation retained. C–D show normalized Rac polarity and spatial contrast. E–I show left Rac/Rho/bound MCA and both end tensions. All three ablations pass two 1,000-s field-settling windows and an independent 10,000-s holdout; endpoints occur at 15,010, 16,010 and 18,010 s. The left displacement spread decreases from 1.968526 to 0.976774 µm but remains length dependent. All three final Rac profiles remain left-biased with contrast 0.205–0.210; the absolute polarity moment is 0.05316, 0.04969 and 0.04743, so only L23 passes the separate 0.05 orientation gate. Negative left tension at L39 is retained. Original model parameters, homogeneous references and protrusion threshold 0.8990616483 are unchanged; only delta_beta_tension is zeroed, while MCA→Rac feedback, mutual inhibition and Rho active stress remain native.'
 save(fig,'global_tension_rho_ablation_endpoints',caption,base)
 fig,axs=plt.subplots(2,3,figsize=(12.5,7.8),layout='constrained')
 timefields=[('left_boundary_displacement','Signed left endpoint displacement (µm)'),('orientation','Rac polarity (negative: left)'),('contrast','Rac spatial contrast (native conc.)')]
 limit=max(t.time[-1] for t in traces.values())*1.02
 for col,(field,title) in enumerate(timefields):
  ax=axs[0,col];label(ax,chr(65+col));ax.set_ylabel(title)
  for L,color in COLORS.items():
   for condition,style in [('native','--'),('ablation','-')]:
    tr=traces[condition,L];y=getattr(tr,field) if field in ['orientation','contrast'] else tr.fields[field][:,0];ax.plot(tr.time,y,color=color,ls=style,lw=1.5 if condition=='ablation' else 1.15)
  ax.set_xscale('symlog',linthresh=10,linscale=.7);ax.set_xlim(0,limit);ax.set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000']);ax.set_xlabel('Time after ablation (s; log above 10 s)')
 for col,(field,title) in enumerate([('rac','Rac (native concentration)'),('rho','Rho (native concentration)'),('mca_bound','Bound MCA (native concentration)')]):
  ax=axs[1,col];label(ax,chr(68+col));ax.set_ylabel(title);ax.set_xlabel('Position x/L');ax.set_xlim(0,1)
  for L,color in COLORS.items():
   for condition,style in [('native','--'),('ablation','-')]:
    tr=traces[condition,L];ax.plot(tr.x/L,tr.fields[field][-1],color=color,ls=style,lw=1.5 if condition=='ablation' else 1.15)
 legend(fig)
 save(fig,'global_tension_rho_ablation_dynamics_profiles','Global tension→Rho ablation from the exact native basal checkpoint, without a countercue. A–C: signed left endpoint displacement, Rac polarity and Rac contrast from ablation onset through verified steady state; time is linear from 0 to 10 s and logarithmic thereafter. Native controls are shown only over their own recorded verified interval. D–F: final full Rac, Rho and bound-MCA profiles; the old pole remains on the left. Solid ablations and dashed native controls share length colors. Concentrations and displacements are absolute native values; only spatial position is divided by L. The ablated states retain spatial polarity while its normalized first-moment magnitude decreases. The intervention removes only tension-dependent Rho activation globally; it does not remove direct tension inhibition of protrusion, Rho-driven cortical stress, or MCA feedback. All 18 initial/final fields and endpoint timecourses are exported; all restart state fields are unchanged at onset.',base)
 (run/'results/status.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',complete=True,ablations=3,native_controls_reused=3,verified_ablations=3,censored=0,left_biased_nonuniform_states=3,old_orientation_gate_passes=1,audit_passed=audit['passed']),indent=2)+'\n')
 print(json.dumps(summary,indent=2));print(OUT)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,default=RUN);a=ap.parse_args();main(a.run.resolve())
