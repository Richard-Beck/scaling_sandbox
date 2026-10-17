"""Complete corrected deBelly L29 relative20 native history, without integrations.
'Full' describes field coverage. Selection retains the canonical scalar 1% gate.
"""
from pathlib import Path
import length_dev_output_layout as devlayout
from datetime import datetime,timezone
import csv,gzip,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'code/core')]
import debelly_corrected_reversal as corrected
n=corrected.n
BASE=ROOT/'data/runs/debelly_corrected_three_test_20261008'
RUN=ROOT/'data/runs/length_debelly_review_20261009'
OUT=ROOT/'figures/length_mechanisms_20261008/dev/L29_fullstate'
def output(name):return devlayout.dev_output_path(OUT,name)
def artifact_key(path):return devlayout.dev_artifact_key(OUT,path)
OUT.mkdir(parents=True,exist_ok=True)
PREFIX='L29_relative20_fullstate'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':200})
SOURCES={str(Path(devlayout.__file__).relative_to(ROOT)):hashlib.sha256(Path(devlayout.__file__).read_bytes()).hexdigest()}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def register(p):SOURCES[str(p.relative_to(ROOT))]=sha(p)
def read(p):register(p);return json.loads(p.read_text())
def write(p,a):p.write_text(json.dumps(a,indent=2)+'\n')
def csvwrite(p,rs):
 with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def panel(ax,label,title):ax.text(-.12,1.04,label,transform=ax.transAxes,fontfamily='DejaVu Sans',fontsize=12,fontweight='bold');ax.set_title(title,pad=10)
def timeaxis(ax,axis='y'):
 getattr(ax,'set_'+axis+'scale')('symlog',linthresh=.5,linscale=.5)
 getattr(ax,'set_'+axis+'ticks')([0,1,10,100,1000,10000],['0','1','10','100','1000','10000'])
 getattr(ax,'set_'+axis+'label')('Time from cue onset (s; symlog, linear 0–0.5 s)')
def save(fig,stem,pdf):
 exports={}
 for ext in ['png','pdf','svg']:
  p=output(f'{PREFIX}_{stem}.{ext}');fig.savefig(p,bbox_inches='tight');exports[artifact_key(p)]=sha(p)
 pdf.savefig(fig,bbox_inches='tight');plt.close(fig);return exports
# Exact canonical selection, independently recomputed from retained scalar histories.
register(BASE/'results/records.json.gz')
with gzip.open(BASE/'results/records.json.gz','rt') as f:case=next(c for c in json.load(f)['cases'] if c['model']=='debelly' and c['condition']=='fixed_relative' and c['length']==29)
assert case['center_boundary_motion'] is False and case['boundary_method']=='independent_endpoints'
assert read(BASE/'results/independent_audit.json')['passed']
checks=[];selected=None
for trial in sorted(case['tests'],key=lambda t:t['duration']):
 path=BASE/trial['trace'];register(path);meta=read(path.with_suffix('.json'));assert sha(path)==meta['trace_sha256'];p=trial['trace_prefix']
 with np.load(path) as z:
  time=z[p+'recovery_time'];tail=time>=9500;o=z[p+'recovery_orientation'][tail];contrast=z[p+'recovery_contrast'][tail]
  assert len(o)==201
  orientation=bool(o.min()>.05 and contrast.min()>=.05 and np.ptp(o)<.005)
  error=float(np.max(abs(o/abs(case['initial_orientation'])-1)));passed=orientation and error<=.01
  assert passed==trial['full_polarity_success']
  checks.append(dict(duration=trial['duration'],orientation_success=orientation,scalar_relative_error_final500=error,canonical_full_polarity_success=passed))
  if passed and selected is None:selected=trial
assert selected['duration']==case['result']['full_polarity']['endpoint']
D=float(selected['duration']);lower=max(t['duration'] for t in checks if t['duration']<D);assert (D-lower)/D<=.01
m=read(RUN/'manifest.json');slot,job=next((i,j) for i,j in enumerate(m['jobs']) if j['case']['length']==29 and j['case']['condition']=='fixed_relative' and j['purpose']=='pass')
assert job['duration']==D and job['case']['checkpoint_sha256']==case['checkpoint_sha256']
path=RUN/f'workers/{slot:03d}.npz';meta=read(path.with_suffix('.json'));assert sha(path)==meta['archive_sha256'];register(path)
assert read(RUN/'independent_audit.json')['passed']
with np.load(path) as z:archive={k:z[k].copy() for k in z.files}
cp_path=ROOT/case['checkpoint'];register(cp_path);assert sha(cp_path)==case['checkpoint_sha256']
with np.load(cp_path) as cp:
 initial={k[6:]:cp[k].copy() for k in cp.files if k.startswith('state_')}
 prefields={k[6:]:cp[k][-1].copy() for k in cp.files if k.startswith('field_')}
x=archive['x'];L=29.;p=n.model_for('debelly').parameters;assert not p.center_boundary_motion
initial_da,initial_db=n.ds._feedback_terms(p,initial['mca_bound'],initial['tension'])
initial_alpha,initial_beta=n.ds._derived_rates(p,initial['mca_bound'],initial['tension'],np.zeros_like(x))
for name,value in [('rac_activation',initial_alpha),('rho_activation',initial_beta),('delta_alpha_mca',initial_da),('delta_beta_tension',initial_db)]:
 np.testing.assert_allclose(value,prefields[name],rtol=1e-13,atol=1e-13)
fields={};phase_fields={}
for phase in ['pulse','recovery']:
 t=archive[phase+'_time'];ff={k:archive[phase+'_field_'+k].copy() for k in n.DB_NAMES}
 for k in n.DB_NAMES:
  assert np.isfinite(ff[k]).all()
  if phase=='pulse':np.testing.assert_allclose(ff[k][0],initial[k],rtol=1e-12,atol=1e-12)
  else:np.testing.assert_allclose(ff[k][0],archive['pulse_field_'+k][-1],rtol=1e-12,atol=1e-12)
 cue=corrected.response_cue(case) if phase=='pulse' else n.ConstantInput(0)
 external=np.array([cue(tt,x,L) for tt in t]);da,db=n.ds._feedback_terms(p,ff['mca_bound'],ff['tension']);alpha,beta=n.ds._derived_rates(p,ff['mca_bound'],ff['tension'],external)
 ff.update(rac_activation=alpha,rho_activation=beta,delta_alpha_mca=da,delta_beta_tension=db)
 sep=ff['right_boundary_displacement']-ff['left_boundary_displacement']
 ff.update(boundary_separation=sep,mean_strain=sep/L,tension_range=np.repeat(np.ptp(ff['tension'],axis=1)[:,None],len(x),axis=1))
 np.testing.assert_allclose(sep,archive[phase+'_field_boundary_separation'],atol=1e-13)
 phase_fields[phase]=ff
assert set(phase_fields['pulse'])==set(prefields),set(phase_fields['pulse'])^set(prefields)
t=np.r_[archive['pulse_time'][:-1],D+archive['recovery_time']]
assert t[0]==0 and np.all(np.diff(t)>0) and t[-1]==D+10000
fields={k:np.concatenate([phase_fields['pulse'][k][:-1],phase_fields['recovery'][k]],axis=0) for k in prefields}
# Reflection respects tensor parity and independently commanded endpoints.
def reflected(name,a):
 other=name.replace('left','RIGHT').replace('right','left').replace('RIGHT','right') if ('left' in name or 'right' in name) else name
 return prefields[other][::-1]*(-1 if name in ['displacement','velocity'] or name in ['raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement'] else 1)
target={k:reflected(k,v) for k,v in prefields.items()}
SPATIAL=['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','rac_activation','rho_activation','delta_alpha_mca','delta_beta_tension']
GLOBALS=['raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement','boundary_separation','mean_strain','tension_range']
assert set(SPATIAL+GLOBALS)==set(fields)
for k in GLOBALS:assert np.max(np.ptp(fields[k],axis=1))<1e-13
GLOBALS_EXTRA=['physical_length','mean_displacement','total_bound_mca','total_unbound_mca','total_mca']
globals_={k:fields[k][:,0] for k in GLOBALS}
globals_.update(physical_length=L+globals_['boundary_separation'],mean_displacement=fields['displacement'].mean(axis=1),total_bound_mca=fields['mca_bound'].sum(axis=1)*L/len(x),total_unbound_mca=fields['mca_unbound'].sum(axis=1)*L/len(x))
globals_['total_mca']=globals_['total_bound_mca']+globals_['total_unbound_mca']
np.testing.assert_allclose(globals_['total_mca'],globals_['total_mca'][0],atol=5e-13)
gtarget={k:target[k][0] for k in GLOBALS}
gtarget.update(physical_length=L+gtarget['boundary_separation'],mean_displacement=target['displacement'].mean(),total_bound_mca=target['mca_bound'].sum()*L/len(x),total_unbound_mca=target['mca_unbound'].sum()*L/len(x));gtarget['total_mca']=gtarget['total_bound_mca']+gtarget['total_unbound_mca']
# Secondary full-field diagnostics do not change the canonical scalar selection.
boundnorm=max(abs(initial['left_boundary_displacement']).max(),abs(initial['right_boundary_displacement']).max());tail=t>=D+9500
errors=[]
for k,a in fields.items():
 norm=boundnorm if k in GLOBALS[:4] else max(float(abs(target[k]).max()),1e-8)
 errors.append(dict(field=k,target_Linf_normalization=norm,normalization_floor=1e-8,maximum_relative_target_error_final500=float(abs(a[tail]-target[k]).max()/norm),relative_target_error_endpoint=float(abs(a[-1]-target[k]).max()/norm),secondary_diagnostic_only=True))
# Saved sampled arrays and both limits of the withdrawal source discontinuity.
arrays=dict(time=t,x=x,cue_duration=np.array(D),recovery_horizon=np.array(10000.),withdrawal_leftlimit_rac_activation=phase_fields['pulse']['rac_activation'][-1],withdrawal_rightlimit_rac_activation=phase_fields['recovery']['rac_activation'][0])
arrays.update({'field_'+k:a for k,a in fields.items()});arrays.update({'initial_'+k:a for k,a in prefields.items()});arrays.update({'reflected_target_'+k:a for k,a in target.items()});arrays.update({'global_'+k:a for k,a in globals_.items()})
np.savez_compressed(output(f'{PREFIX}_fields.npz'),**arrays)
csvwrite(output(f'{PREFIX}_duration_selection.csv'),checks);csvwrite(output(f'{PREFIX}_target_errors.csv'),errors)
labels={'rac':'Rac','rho':'Rho','mca_bound':'Bound MCA','mca_unbound':'Unbound MCA','displacement':'Membrane displacement','velocity':'Cortical velocity','tension':'Membrane tension','rac_activation':'Rac activation rate α','rho_activation':'Rho activation rate β','delta_alpha_mca':'MCA contribution to α','delta_beta_tension':'Tension contribution to β','raw_left_displacement':'Raw left endpoint displacement','raw_right_displacement':'Raw right endpoint displacement','left_boundary_displacement':'Left endpoint displacement','right_boundary_displacement':'Right endpoint displacement','boundary_separation':'Boundary separation change','mean_strain':'Mean membrane strain','tension_range':'Tension range','physical_length':'Physical cell length','mean_displacement':'Mean membrane displacement','total_bound_mca':'Total bound MCA','total_unbound_mca':'Total unbound MCA','total_mca':'Total MCA'}
units={k:'native concentration units' for k in ['rac','rho','mca_bound','mca_unbound']};units.update(displacement='µm',velocity='µm/s',tension='pN',rac_activation='native activation units',rho_activation='native activation units',delta_alpha_mca='native activation units',delta_beta_tension='native activation units',raw_left_displacement='µm',raw_right_displacement='µm',left_boundary_displacement='µm',right_boundary_displacement='µm',boundary_separation='µm',mean_strain='dimensionless',tension_range='pN',physical_length='µm',mean_displacement='µm',total_bound_mca='concentration × µm',total_unbound_mca='concentration × µm',total_mca='concentration × µm')
exports={};combined=output(f'{PREFIX}.pdf')
with PdfPages(combined) as pdf:
 fig,axs=plt.subplots(3,4,figsize=(16,12.5),layout='constrained')
 for i,k in enumerate(SPATIAL):
  ax=axs.flat[i];a=fields[k];signed=k in ['displacement','velocity','tension'];vmin,vmax=float(a.min()),float(a.max())
  if signed:
   lim=max(abs(vmin),abs(vmax),1e-12);norm=Normalize(-lim,lim);cmap='RdBu_r'
  else:norm=Normalize(vmin,vmax);cmap='viridis'
  mesh=ax.pcolormesh(x,t,a,shading='nearest',cmap=cmap,norm=norm,rasterized=True)
  panel(ax,chr(65+i),labels[k]);ax.axhline(D,color='white',lw=2,ls='--');ax.axhline(D,color='#333333',lw=.8,ls='--');timeaxis(ax);ax.set(xlabel='Reference position x (µm)',xlim=(0,L),ylim=(0,t[-1]));fig.colorbar(mesh,ax=ax,shrink=.88,pad=.025,label=units[k])
 axs.flat[-1].axis('off');fig.legend(handles=[Line2D([],[],color='#333333',ls='--',label=f'Cue withdrawal: {D:.6f} s')],loc='outside lower center',frameon=False)
 exports.update(save(fig,'kymographs',pdf))
 fig,axs=plt.subplots(3,4,figsize=(16,11),layout='constrained')
 for i,k in enumerate(GLOBALS+GLOBALS_EXTRA):
  ax=axs.flat[i];a=globals_[k];panel(ax,chr(65+i),labels[k]);ax.plot(t,a,color='#0072B2',lw=1.3);ax.axhline(gtarget[k],color='#222222',ls=':',lw=1.1);ax.axvline(D,color='#D55E00',ls='--',lw=1);timeaxis(ax,'x');ax.set(ylabel=units[k],xlim=(0,t[-1]));ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.grid(alpha=.14)
 fig.legend(handles=[Line2D([],[],color='#0072B2',label='Native trajectory'),Line2D([],[],color='#222222',ls=':',label='Reflected pre-cue target'),Line2D([],[],color='#D55E00',ls='--',label='Cue withdrawal')],loc='outside lower center',ncol=3,frameon=False)
 exports.update(save(fig,'globals',pdf))
 fig,axs=plt.subplots(3,4,figsize=(16,11),layout='constrained')
 for i,k in enumerate(SPATIAL):
  ax=axs.flat[i];panel(ax,chr(65+i),labels[k]);ax.plot(x,prefields[k],color='#999999',lw=1.2);ax.plot(x,target[k],color='#222222',ls='--',lw=1.6);ax.plot(x,fields[k][np.searchsorted(t,D+9500)],color='#D55E00',lw=1.1);ax.plot(x,fields[k][-1],color='#0072B2',lw=1.1);ax.set(xlabel='Reference position x (µm)',ylabel=units[k],xlim=(0,L));ax.grid(alpha=.14)
 axs.flat[-1].axis('off');fig.legend(handles=[Line2D([],[],color='#999999',label='Pre-cue basal state'),Line2D([],[],color='#222222',ls='--',label='Reflected pre-cue target'),Line2D([],[],color='#D55E00',label='9500 s after withdrawal'),Line2D([],[],color='#0072B2',label='10000 s after withdrawal')],loc='outside lower center',ncol=4,frameon=False)
 exports.update(save(fig,'profiles',pdf))
exports[artifact_key(combined)]=sha(combined)
caption=(f'Corrected deBelly, L = 29 µm, relative-width20% right patch (width5.8 µm, peak3.6 native units). Cue duration {D:.12f} s is the shortest of79 tested durations passing the canonical scalar strength gate: final500 s of10000 s basal recovery have right Rac orientation>0.05, Rac contrast≥0.05, orientation range<0.005, and Rac first-moment magnitude within1% of its pre-cue magnitude. The neighboring tested failure is {lower:.12f} s; bracket width {(D-lower)/D*100:.4f}%. Full describes field coverage, not an additional success definition. Kymographs show all11 spatial native recorded fields; seven spatially constant native fields are line graphs, alongside five explicitly derived global readouts. Raw membrane displacement is retained. Time uses a true symlog coordinate:0–0.5 s is linear, later time is logarithmic; time0 is included. Dashed lines mark withdrawal. The native0.5-s step, pulse0.5-s output and recovery2.5-s output are unchanged. Times are measured from cue onset; final time is {D+10000:.12f} s. Derived activation fields use the actual cue before withdrawal and basal input at/after withdrawal; at the duplicated withdrawal sample the right limit is displayed and both rate limits are exported. Rac/MCA/Rho and scalar derived fields reflect evenly; displacement/velocity reflect oddly; left/right endpoint displacements swap and change sign. Every field is in native reference coordinates. Secondary pointwise target errors use field-specific reflected-target maximum norms, with a common endpoint-pair amplitude for zero endpoint targets and1e−8 floors; these diagnostics do not select the cue. The separate canonical0.1% field-tail stability flag is false for this selected scalar-success trial.')
# Exact native record identities for derived fields and full field coverage.
audit=dict(passed=True,selection='canonical scalar full_polarity_success, independently recomputed for all79 retained trials',duration=D,lower_failing_duration=lower,relative_bracket_width=(D-lower)/D,recorded_native_fields=sorted(fields),spatial_native_fields=SPATIAL,global_native_fields=GLOBALS,auxiliary_global_readouts=GLOBALS_EXTRA,native_fields_omitted=[],checkpoint_identity=True,phase_continuity=True,mca_mass_max_error=float(abs(globals_['total_mca']-globals_['total_mca'][0]).max()),time_monotonic=True,pulse_output_interval=.5,recovery_output_interval=2.5,number_of_time_samples=len(t),spatial_cells=len(x),canonical_scalar_gate_pass=True,separate_field_tail_pass=meta['metrics']['field_tail_pass'],secondary_allfield_max_relative_target_error_final500=max(a['maximum_relative_target_error_final500'] for a in errors),reconstructed_fields=['rac_activation','rho_activation','delta_alpha_mca','delta_beta_tension','boundary_separation','mean_strain','tension_range'],reconstruction='Exact current native record definitions, including cue-on/withdrawn input and bound MCA/tension feedback.',field_units=units)
register(Path(__file__));register(Path(n.ds.__file__));register(Path(corrected.__file__))
write(output(f'{PREFIX}_audit.json'),audit)
(output(f'{PREFIX}_caption.txt')).write_text(caption+'\n')
manifest=dict(artifact_root='../..',generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',passed=True,description='Complete corrected deBelly L29 relative20 shortest canonical scalar1% reversal',legend=caption,figure_sha256=exports,data_sha256={artifact_key(p):sha(p) for p in [output(f'{PREFIX}_fields.npz'),output(f'{PREFIX}_duration_selection.csv'),output(f'{PREFIX}_target_errors.csv'),output(f'{PREFIX}_audit.json')]},source_sha256=SOURCES,reproduction=['OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_dev_fullstate_kymographs.py'],parameters=vars(p),activation_definitions={'delta_alpha_mca':'rac_mechanical_feedback * rho_to_rac_coupling * 0.5 * [1 − tanh((bound_MCA − mca_rac_threshold)/mca_switch_scale)]','delta_beta_tension':'rho_mechanical_feedback * rac_to_rho_coupling * 0.5 * [1 + tanh((tension − tension_rho_threshold)/tension_switch_scale)]','rac_activation':'rac_basal_activation + cue (opto_target=rac) + delta_alpha_mca','rho_activation':'rho_basal_activation + delta_beta_tension','cue':'Right cosine patch, peak3.6, support5.8µm; identically0 after withdrawal.'},audit=audit)
write(output(f'{PREFIX}.manifest.json'),manifest)
for stem in ['kymographs','globals','profiles']:
 attached=dict(manifest);attached['figure_sha256']={name:digest for name,digest in exports.items() if Path(name).name.startswith(PREFIX+'_'+stem+'.')};attached['view']=stem
 write(output(f'{PREFIX}_{stem}.manifest.json'),attached)
 (output(f'{PREFIX}_{stem}_caption.txt')).write_text(caption+'\n')
html='<html lang="en"><meta charset="utf-8"><title>deBelly L29 full native history</title><style>body{max-width:1500px;margin:30px auto;padding:0 24px;font:16px/1.5 system-ui}img{width:100%;height:auto}p{max-width:1100px}</style><h1>deBelly: L = 29 µm, relative20% cue</h1><p>'+caption+'</p><p><a href="'+devlayout.dev_viewer_link(OUT,PREFIX+'.pdf')+'">Three-page PDF</a> · <a href="'+devlayout.dev_viewer_link(OUT,PREFIX+'.manifest.json')+'">Manifest</a> · <a href="'+devlayout.dev_viewer_link(OUT,PREFIX+'_fields.npz')+'">All native fields</a> · <a href="'+devlayout.dev_viewer_link(OUT,PREFIX+'_target_errors.csv')+'">Secondary target errors</a> · <a href="'+devlayout.dev_viewer_link(OUT,PREFIX+'_duration_selection.csv')+'">Duration selection</a></p>'+''.join('<h2>'+title+'</h2><a href="'+devlayout.dev_viewer_link(OUT,PREFIX+'_'+stem+'.pdf')+'"><img src="'+devlayout.dev_viewer_link(OUT,PREFIX+'_'+stem+'.png')+'" alt="'+title+'"></a>' for stem,title in [('kymographs','Every spatial native field'),('globals','Every native global scalar and derived totals'),('profiles','Pre-cue, reflected target and late recovery profiles')])+'</html>'
(output(f'{PREFIX}.html')).write_text(html)
print(json.dumps(dict(output=str(output(f'{PREFIX}.html')),duration=D,native_field_count=len(fields),worst_secondary_relative_target_error=audit['secondary_allfield_max_relative_target_error_final500'],passed=True),indent=2))
