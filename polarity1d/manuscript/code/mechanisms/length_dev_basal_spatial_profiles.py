"""Development-only complete basal spatial-field comparison at L23,31,39.

Uses verified shared corrected checkpoints. Only x is normalized by length;
all saved field amplitudes and raw displacement retain their native values.
"""
from pathlib import Path
from datetime import datetime,timezone
from dataclasses import asdict
import csv,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import length_dev_output_layout as layout
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core')]
import debelly_corrected_reversal as corrected
n=corrected.n
BASE=corrected.DEFAULT
OUT=ROOT/'figures/length_mechanisms_20261008/dev/resting_state'
STEM='basal_spatial_profiles_L23_L31_L39'
GLOBSTEM='basal_global_fields_L23_L31_L39'
LENGTHS=[23,31,39]
COLORS={23:'#0072B2',31:'#762a83',39:'#D55E00'}
GEOS=['fixed_physical','fixed_relative','fixed_slope']
SPATIAL=[('rac','Rac','Native concentration'),('rho','Rho','Native concentration'),('mca_bound','Bound MCA','Native MCA density'),('mca_unbound','Unbound MCA','Native MCA density'),('displacement','Membrane displacement u','Displacement (µm)'),('velocity','Cortical velocity v','Velocity (µm/s)'),('tension','Membrane tension T','Tension (pN)'),('rac_activation','Rac activation α','Native activation units'),('rho_activation','Rho activation β','Native activation units'),('delta_alpha_mca','MCA contribution Δα','Native activation units'),('delta_beta_tension','Tension contribution Δβ','Native activation units')]
GLOBALS=[('raw_left_displacement','Raw left displacement','Displacement (µm)'),('raw_right_displacement','Raw right displacement','Displacement (µm)'),('left_boundary_displacement','Left endpoint displacement','Displacement (µm)'),('right_boundary_displacement','Right endpoint displacement','Displacement (µm)'),('boundary_separation','Endpoint separation change','Length change (µm)'),('mean_strain','Mean strain','Strain'),('tension_range','Tension range','Tension (pN)')]
SOURCES={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':10,'pdf.fonttype':42,'svg.fonttype':'none'})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):SOURCES[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
def output(name):return layout.dev_output_path(OUT,name)
def key(p):return layout.dev_artifact_key(OUT,p)
def csvwrite(path,rows):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def panel(ax,letter,title,ylabel):
 ax.text(-.14,1.04,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',fontfamily='DejaVu Sans');ax.set_title(title,pad=10);ax.set_ylabel(ylabel);ax.grid(alpha=.13)
def main():
 p=n.model_for('debelly').parameters;assert p.center_boundary_motion is False and p.assay_mode=='reflection_symmetric_extension'
 frozen=read(BASE/'provenance.json');verified={}
 for name in ['code/core/native.py','code/core/models/debelly.py','code/core/src/polarity1d/debelly.py','code/debelly_corrected_reversal.py']:
  digest=frozen['sources'][name];assert sha(BASE/'frozen_sources'/name)==digest and sha(ROOT/name)==digest,name;SOURCES[name]=digest;verified[name]=digest
 for path in [Path(__file__),Path(layout.__file__)]:SOURCES[str(path.relative_to(ROOT))]=sha(path)
 manifest=read(BASE/'manifest.json');assert manifest['center_boundary_motion'] is False
 assert read(BASE/'results/independent_audit.json')['passed']
 profiles={};globalrows=[];spatialrows=[];summary=[];audits=[];npz={}
 for L in LENGTHS:
  cases=[next(c for c in manifest['cases'] if c['model']=='debelly' and c['condition']==geo and c['length']==L) for geo in GEOS]
  assert all(c['eligible'] and c['preparation_status']=='polarized' and c['center_boundary_motion'] is False for c in cases)
  assert len({c['checkpoint'] for c in cases})==1 and len({c['checkpoint_sha256'] for c in cases})==1
  c=cases[0];cp=ROOT/c['checkpoint'];assert sha(cp)==c['checkpoint_sha256'];SOURCES[str(cp.relative_to(ROOT))]=sha(cp)
  prep=read(ROOT/c['corrected_preparation']);assert prep['summary']['status']=='polarized'
  with np.load(cp) as z:
   x=z['x'].copy();time=float(z['time'][-1]);fields={k[6:]:z[k][-1].copy() for k in z.files if k.startswith('field_')}
   assert set(fields)=={a[0] for a in SPATIAL+GLOBALS}
   for name in n.DB_NAMES:np.testing.assert_array_equal(z['state_'+name],fields[name])
   assert np.allclose(x,(np.arange(201)+.5)*L/201,atol=1e-12)
  for name,value in fields.items():assert np.isfinite(value).all() and value.shape==(201,)
  ori=float(fields['rac']@(2*x/L-1)/fields['rac'].sum());contrast=float(np.ptp(fields['rac']))
  assert ori<-.05 and contrast>=.05;assert abs(ori-c['initial_orientation'])<1e-12 and abs(contrast-c['initial_contrast'])<1e-12
  alpha,beta=n.ds._derived_rates(p,fields['mca_bound'],fields['tension'],np.zeros(201));da,db=n.ds._feedback_terms(p,fields['mca_bound'],fields['tension'])
  derived={'rac_activation':alpha,'rho_activation':beta,'delta_alpha_mca':da,'delta_beta_tension':db}
  residuals={}
  for name,value in derived.items():residuals[name]=float(np.max(abs(value-fields[name])));np.testing.assert_allclose(value,fields[name],rtol=0,atol=1e-14)
  np.testing.assert_allclose(fields['rac_activation'],p.rac_basal_activation+fields['delta_alpha_mca'],rtol=0,atol=1e-14);np.testing.assert_allclose(fields['rho_activation'],p.rho_basal_activation+fields['delta_beta_tension'],rtol=0,atol=1e-14)
  gl={}
  for name,title,units in GLOBALS:assert np.ptp(fields[name])==0;gl[name]=float(fields[name][0])
  np.testing.assert_allclose(gl['boundary_separation'],gl['right_boundary_displacement']-gl['left_boundary_displacement'],rtol=0,atol=1e-14)
  np.testing.assert_allclose(gl['mean_strain'],gl['boundary_separation']/L,rtol=0,atol=1e-14)
  np.testing.assert_allclose(gl['tension_range'],np.ptp(fields['tension']),rtol=0,atol=1e-14)
  globalrows.append(dict(length_um=L,preparation_checkpoint_time_s=time,countercue_time_s=0.,initial_rac_polarity=ori,initial_rac_contrast=contrast,**gl,physical_endpoint_separation_um=L+gl['boundary_separation']))
  for i in range(201):spatialrows.append(dict(length_um=L,preparation_checkpoint_time_s=time,countercue_time_s=0.,x_um=float(x[i]),x_over_L=float(x[i]/L),**{name:float(fields[name][i]) for name,title,units in SPATIAL}))
  for name,title,units in SPATIAL:summary.append(dict(length_um=L,field=name,units=units,left_cell=float(fields[name][0]),right_cell=float(fields[name][-1]),minimum=float(fields[name].min()),maximum=float(fields[name].max())))
  profiles[L]=(x,fields,time);npz[f'L{L}_x_um']=x;npz[f'L{L}_x_over_L']=x/L;npz[f'L{L}_preparation_checkpoint_time_s']=np.array(time)
  for name,value in fields.items():npz[f'L{L}_'+name]=value
  audits.append(dict(length_um=L,checkpoint=str(cp.relative_to(ROOT)),checkpoint_sha256=sha(cp),all_three_assays_share_exact_checkpoint=True,preparation_qualified=True,native_checkpoint_time_s=time,spatial_fields=11,uniform_native_globals=7,derived_rate_max_absolute_residual=residuals,global_identities_pass=True,raw_displacement_preserved=True))
 csvwrite(output(STEM+'.csv'),spatialrows);csvwrite(output(GLOBSTEM+'.csv'),globalrows);csvwrite(output(STEM+'_summary.csv'),summary);np.savez_compressed(output(STEM+'.npz'),**npz)
 handles=[Line2D([],[],color=COLORS[L],lw=2,label=f'L = {L} µm') for L in LENGTHS]
 fig,axs=plt.subplots(4,3,figsize=(12,12.3),layout='constrained',sharex=True)
 for index,(name,title,units) in enumerate(SPATIAL):
  ax=axs.flat[index];panel(ax,chr(65+index),title,units)
  for L in LENGTHS:x,f,t=profiles[L];ax.plot(x/L,f[name],color=COLORS[L],lw=1.8)
  ax.set_xlim(0,1);ax.set_xticks([0,.25,.5,.75,1]);ax.set_xlabel('Position x/L')
  if name in ['displacement','velocity','tension']:ax.axhline(0,color='#999999',lw=.7)
 axs.flat[11].set_axis_off();axs.flat[11].legend(handles=handles,loc='center',frameon=False,title='Shared pre-countercue basal checkpoints',fontsize=11,title_fontsize=10)
 outputs={}
 for ext in ['png','pdf','svg']:
  path=output(STEM+'.'+ext);fig.savefig(path,dpi=220,facecolor='white');outputs[key(path)]=sha(path)
 plt.close(fig)
 fig,axs=plt.subplots(2,4,figsize=(12,6.8),layout='constrained')
 for index,(name,title,units) in enumerate(GLOBALS):
  ax=axs.flat[index];panel(ax,chr(65+index),title,units);values=[next(r[name] for r in globalrows if r['length_um']==L) for L in LENGTHS]
  ax.plot(LENGTHS,values,color='#bbbbbb',lw=1.2,zorder=1)
  for L,v in zip(LENGTHS,values):ax.plot(L,v,'o',color=COLORS[L],ms=6,zorder=2)
  ax.set_xlim(21,41);ax.set_xticks(LENGTHS);ax.set_xlabel('Cell length (µm)');ax.ticklabel_format(axis='y',style='plain',useOffset=False)
  if name in ['raw_right_displacement','right_boundary_displacement']:ax.set_ylim(-.01,.01)
 axs.flat[7].set_axis_off();axs.flat[7].legend(handles=handles,loc='center',frameon=False,title='Same checkpoints as spatial profiles')
 for ext in ['png','pdf','svg']:
  path=output(GLOBSTEM+'.'+ext);fig.savefig(path,dpi=220,facecolor='white');outputs[key(path)]=sha(path)
 plt.close(fig)
 caption='Complete pre-countercue spatial profiles at the exact qualifying corrected basal checkpoints for L = 23, 31 and 39 µm. Each facet overlays all three lengths with the same colors: blue, purple and orange. A–G show native Rac, Rho, bound/unbound MCA, raw membrane displacement u, cortical velocity v and membrane tension T. H–K show saved Rac/Rho activation α/β and their MCA/tension increments Δα/Δβ. All 201 native cell-center values are used; only x is divided by nominal L. Field amplitudes retain native absolute values and displacement is not centered. x/L = 0 is the old left Rac pole; x/L = 1 is the end targeted by the subsequent countercue. Preparation checkpoint clocks are 49010, 42010 and 38010 s respectively; each exact state becomes countercue time zero. These are the previously verified basal endpoints, with no additional settling or integration. All three response assays share each checkpoint. Saved derived rates are independently verified at zero extra cue: α = basal Rac activation + Δα and β = basal Rho activation + Δβ. The seven uniformly stored native globals appear in the companion figure.'
 globalcaption='All seven native global fields stored with the same L = 23, 31 and 39 µm pre-countercue checkpoints. A–D show raw left/right displacement and the actual left/right endpoint displacement; the pairs agree because the corrected model does not center endpoint motion. E is the displacement-induced change in endpoint separation, not the full cell length. F is that change divided by nominal L; G is maximum minus minimum native tension. These fields are constant across the spatial grid and are therefore compared at three cell lengths rather than repeated as horizontal spatial lines. Colors match the companion full spatial profiles. Values are copied from the exact basal checkpoints; no amplitude normalization, integration or new settling.'
 output(STEM+'_caption.txt').write_text(caption+'\n');output(GLOBSTEM+'_caption.txt').write_text(globalcaption+'\n')
 for name in [STEM+'.csv',GLOBSTEM+'.csv',STEM+'_summary.csv',STEM+'.npz',STEM+'_caption.txt',GLOBSTEM+'_caption.txt']:
  path=output(name);outputs[key(path)]=sha(path)
 methods=dict(spatial_coverage=[dict(native_name=name,title=title,units=units) for name,title,units in SPATIAL],global_coverage=[dict(native_name=name,title=title,units=units) for name,title,units in GLOBALS],saved_field_coverage='Every native field_ array in each checkpoint is included once: 11 spatial fields and 7 uniform globals.',normalization='Only spatial coordinate x/L. Raw membrane displacement and every field amplitude remain uncentered native values.',time='Native preparation checkpoint time is retained separately from countercue time zero.',derived_identities='Native _derived_rates and _feedback_terms evaluated algebraically from each saved bound-MCA and tension profile at zero extra cue. No solver step.',native_globals='Saved arrays verified uniform; separation=right endpoint displacement−left endpoint displacement; strain=separation/L; tension range=max(T)−min(T).')
 record=dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',description=caption,legend=caption,companion_description=globalcaption,methods=methods,sha256=outputs,source_sha256=SOURCES,renderer_sha256=sha(Path(__file__)),parameters=asdict(p),frozen_parameter_source_identity=verified,independent_audit=audits,reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python '+str(Path(__file__))],new_integrations=0,current_report_changed=False)
 output(STEM+'.manifest.json').write_text(json.dumps(record,indent=2)+'\n')
 companion=dict(record);companion.update(description=globalcaption,legend=globalcaption,spatial_figure_description=caption);output(GLOBSTEM+'.manifest.json').write_text(json.dumps(companion,indent=2)+'\n')
 print(json.dumps(globalrows,indent=2));print(output(STEM+'.png'));print(output(GLOBSTEM+'.png'))
if __name__=='__main__':main()
