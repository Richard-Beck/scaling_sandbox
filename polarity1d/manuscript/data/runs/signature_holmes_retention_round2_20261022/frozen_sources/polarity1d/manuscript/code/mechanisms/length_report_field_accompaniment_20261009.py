"""Report-only integration assets for L29 fields and the paired Fig9F budget.
Reuses immutable native histories; no long integrations or native budget changes.
"""
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'code/mechanisms'),str(ROOT/'code'),str(ROOT/'code/core')]
import length_quality_debelly_traces_20261009 as original
from length_dev_output_layout import dev_output_path
n=original.n
FIG=ROOT/'figures/length_mechanisms_20261008'
OUT=FIG/'quality_20261009_field_accompaniment'
DATA=ROOT/'data/runs/length_mechanisms_20261008/quality_20261009_field_accompaniment'
OUT.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
write=lambda p,z:p.write_text(json.dumps(z,indent=2)+'\n')
FLOW_CAPTION='Figure 9A–E. Attachment travel and mechanical screening under the constant-slope cue (0.24) at all twelve eligible lengths. A compares uninterrupted-cue Rac/MCA zero with the shortest tested orientation-reversing pulse. B shows initial MCA half-mass/centroid distance to the midpoint and centroid displacement during that pulse. C is distance to MCA zero divided by its crossing time, a distribution speed. D–E freeze the L = 23 Rho profile at MCA zero in reference coordinates and extend its low tail, with zero membrane and endpoint velocities. D shows the uniform-basal-MCA cortical solutions; E compares speed at x = 5 µm with and without MCA friction. Cue-relative native rate budgets are shown in Fig. 9F-A,B.'
# Preserve the original five plots and their data exactly; discard the obsolete F axis.
def save_flow(fig,stem,caption,methods):
 if stem!='debelly_flow':plt.close(fig);return
 assert len(fig.axes)==6;fig.delaxes(fig.axes[5]);fig.set_layout_engine(None)
 positions=[[.072,.59,.235,.32],[.392,.59,.235,.32],[.712,.59,.235,.32],[.17,.12,.29,.32],[.61,.12,.29,.32]]
 for ax,bounds in zip(fig.axes,positions):ax.set_position(bounds)
 exports={}
 for ext in ['png','pdf','svg']:
  p=OUT/('debelly_flow_AE.'+ext);fig.savefig(p,dpi=220,facecolor='white');exports[p.name]=sha(p)
 plt.close(fig)
 sources=dict(original.SOURCES)
 for p in [Path(__file__),Path(original.__file__),ROOT/'code/mechanisms/length_dev_output_layout.py']:sources[str(p.relative_to(ROOT))]=sha(p)
 write(OUT/'debelly_flow_AE.manifest.json',dict(producer='finished',generated_at=datetime.now(timezone.utc).isoformat(),passed=True,description='Unchanged Figure9A–E with old zero-crossing budget panel removed.',legend=FLOW_CAPTION,sha256=exports,source_sha256=sources,renderer_sha256=sha(Path(__file__)),reproduction=['OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_report_field_accompaniment_20261009.py'],methods=methods))
 (OUT/'debelly_flow_AE_caption.txt').write_text(FLOW_CAPTION+'\n')
original.save=save_flow
measure=ROOT/'data/runs/length_mechanisms_20261008/quality_20261009_debelly/memory_measurements.csv'
original.SOURCES[str(measure.relative_to(ROOT))]=sha(measure)
rows=[]
for row in csv.DictReader(measure.open()):
 for k,v in row.items():
  if k not in ['condition','field_tail_pass']:
   try:row[k]=float(v)
   except (ValueError,TypeError):pass
 rows.append(row)
original.flow_figure(rows)
# Audit the narrated times directly on the archived field export.
source=dev_output_path(FIG/'dev','L29_relative20_fullstate_fields.npz')
with np.load(source) as z:a={k:z[k].copy() for k in z.files}
t=a['time'];D=float(a['cue_duration']);params=n.model_for('debelly').parameters
rac0,rho0=n.ds._homogeneous_rac_rho(params,.5);threshold=1.3*rac0
switch=n.ds.smooth_threshold(a['field_rac'],params.protrusion_switch_scale,threshold)
readouts=[]
for requested in [70.5,D-60,D,430.2551020408163,435.2551020408163,442.7551020408163,500.,750.,8817.755102040816]:
 i=int(np.argmin(abs(t-requested)))
 readouts.append(dict(requested_time=requested,actual_time=float(t[i]),time_after_withdrawal=float(t[i]-D),left_rac=float(a['field_rac'][i,0]),right_rac=float(a['field_rac'][i,-1]),left_switch=float(switch[i,0]),right_switch=float(switch[i,-1]),left_boundary_displacement=float(a['field_left_boundary_displacement'][i,0]),right_boundary_displacement=float(a['field_right_boundary_displacement'][i,0]),left_tension=float(a['field_tension'][i,0]),maximum_tension=float(a['field_tension'][i].max())))
with (DATA/'L29_endpoint_story.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(readouts[0]));w.writeheader();w.writerows(readouts)
# Continuous exact switch crossing from retained samples, without altering states.
idx=np.flatnonzero((t>D)&(a['field_rac'][:,-1]<threshold))[0]
crossing=float(np.interp(threshold,a['field_rac'][idx-1:idx+1,-1][::-1],t[idx-1:idx+1][::-1]))
write(DATA/'endpoint_story_audit.json',dict(passed=True,source_sha256=sha(source),duration=D,rac_rest=rac0,rac_switch_threshold=threshold,right_switch_threshold_crossing_from_onset=crossing,endpoint_relaxation_timescale=-.5/np.log(.98),readouts=readouts,units='seconds,micrometres,pN; chemical fields in native concentration units',no_new_integrations=True))
KCAP='Figure 8S1. Complete spatial history accompanying Figure 8: L = 29 µm, relative-width 20% right cosine patch (support 5.8 µm, peak 3.6), withdrawn at 282.755102 s, then 10000 s basal recovery. A–D: Rac, Rho, bound and unbound MCA. E–G: absolute membrane displacement, cortical velocity and membrane tension. H–K: total Rac/Rho activation rates and their MCA/tension contributions. Every spatial native field is shown; spatially uniform native fields are in Figure 8S2. Colorbars retain absolute native ranges. Time coordinates are symlog (linear from 0 to 0.5 s, logarithmic thereafter); dashed lines mark withdrawal. Pulse and recovery outputs are 0.5 and 2.5 s apart. Activation rates use basal input at the withdrawal sample; both rate limits are exported. Cue selection uses the canonical scalar 1% polarity gate.'
GCAP='Figure 8S2. Global variables from the exact Figure 8S1 history, with the same symlog time and withdrawal marker. A–D: raw and physical endpoint displacements; they coincide under independent endpoint commands. E–G: boundary separation change, mean strain and tension range. H–I: physical length and mean membrane displacement. J–L: bound, unbound and total MCA amounts. A–G are all seven spatially uniform native fields; H–L are derived readouts. Dotted lines show the reflected pre-cue target. Displacement retains its absolute offset; no centering or state reset is applied. Total MCA is conserved. The first persistent scalar 1% arrival is 8817.755 s from cue onset; it is distinct from cue duration.'
FCAP='Figure 9F. Native bound-MCA polarity-rate budget under the constant-slope cue (0.24) at all twelve eligible lengths. A: 60 s after cue onset. B: 60 s before each length’s shortest tested orientation-reversing pulse ends. Both states are cue-on. Orange endpoint advection, blue Rho advection and gray remaining terms sum to the black dashed net rate; positive rates erase the old attachment bias. Remaining terms include membrane motion, diffusion, exchange and uniform restoration. The changing polarity denominator is included. A uses 60 s exactly; B uses the nearest retained 0.5-s sample (maximum offset 0.229592 s; all actual times are tabulated). These are fixed-state native source identities with frozen upwind donors, not independently evolved interventions. This gradient comparison uses different cues and trials from the L = 29 relative-patch example in Figure 8S1.'
# Report-width redraws keep labels readable without changing a sampled value.
plt.rcParams.update({'font.size':9.5,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8})
def save_fields(fig,stem,caption,panels):
 exports={}
 for ext in ['png','pdf','svg']:
  path=OUT/(stem+'.'+ext);fig.savefig(path,dpi=220,facecolor='white');exports[path.name]=sha(path)
 plt.close(fig)
 sources={str(path.relative_to(ROOT)):sha(path) for path in [Path(__file__),source,dev_output_path(FIG/'dev','L29_relative20_fullstate.manifest.json'),Path(n.ds.__file__),ROOT/'code/mechanisms/length_dev_output_layout.py']}
 write(OUT/(stem+'.manifest.json'),dict(passed=True,producer='finished',generated_at=datetime.now(timezone.utc).isoformat(),legend=caption,sha256=exports,source_sha256=sources,renderer_sha256=sha(Path(__file__)),panels=panels,reproduction=['OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_report_field_accompaniment_20261009.py'],methods='Exact existing native samples; absolute values; true symlog coordinate with0–0.5s linear; no native integrations.'))
 (OUT/(stem+'_caption.txt')).write_text(caption+'\n')
def time_axis(ax,direction='y'):
 getattr(ax,'set_'+direction+'scale')('symlog',linthresh=.5,linscale=.5)
 getattr(ax,'set_'+direction+'ticks')([0,1,10,100,1000,10000],['0','1','10','100','1000','10000'])
 getattr(ax,'set_'+direction+'label')('Time from cue onset (s)')
def panel_label(ax,i,title):
 ax.text(-.13,1.04,chr(65+i),transform=ax.transAxes,fontfamily='DejaVu Sans',fontsize=12,fontweight='bold');ax.set_title(title,loc='left',pad=9)
spatial=['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','rac_activation','rho_activation','delta_alpha_mca','delta_beta_tension']
titles=['Rac','Rho','Bound MCA','Unbound MCA','Membrane displacement','Cortical velocity','Membrane tension','Rac activation α','Rho activation β','MCA contribution to α','Tension contribution to β']
fig,axs=plt.subplots(3,4,figsize=(12,10.8),layout='constrained')
for i,(key,title) in enumerate(zip(spatial,titles)):
 ax=axs.flat[i];values=a['field_'+key];signed=key in ['displacement','velocity','tension'];lim=float(abs(values).max())
 options=dict(vmin=-lim,vmax=lim,cmap='RdBu_r') if signed else dict(vmin=float(values.min()),vmax=float(values.max()),cmap='viridis')
 mesh=ax.pcolormesh(a['x'],t,values,shading='nearest',rasterized=True,**options);panel_label(ax,i,title);time_axis(ax)
 ax.axhline(D,color='white',lw=1.8,ls='--');ax.axhline(D,color='#333333',lw=.7,ls='--');ax.set(xlabel='Position (µm)',xlim=(0,29),ylim=(0,t[-1]),xticks=[0,10,20,29]);ax.tick_params(labelsize=8)
 unit={'displacement':'µm','velocity':'µm/s','tension':'pN'}.get(key,'native conc.' if i<4 else 'activation units')
 fig.colorbar(mesh,ax=ax,pad=.025,shrink=.88,label=unit).ax.tick_params(labelsize=8)
axs.flat[-1].axis('off');save_fields(fig,'L29_fields',KCAP,spatial)
global_names=['raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement','boundary_separation','mean_strain','tension_range','physical_length','mean_displacement','total_bound_mca','total_unbound_mca','total_mca']
global_titles=['Raw left endpoint','Raw right endpoint','Left endpoint','Right endpoint','Separation change','Mean strain','Tension range','Physical cell length','Mean displacement','Total bound MCA','Total unbound MCA','Total MCA']
fig,axs=plt.subplots(3,4,figsize=(12,9.5),layout='constrained')
for i,(key,title) in enumerate(zip(global_names,global_titles)):
 ax=axs.flat[i];values=a['global_'+key];panel_label(ax,i,title);ax.plot(t,values,color='#0072B2',lw=1.2)
 if key in ['physical_length','mean_displacement','total_bound_mca','total_unbound_mca','total_mca']:
  target_value={'physical_length':29+a['reflected_target_boundary_separation'][0],'mean_displacement':a['reflected_target_displacement'].mean(),'total_bound_mca':a['reflected_target_mca_bound'].sum()*29/len(a['x']),'total_unbound_mca':a['reflected_target_mca_unbound'].sum()*29/len(a['x']),'total_mca':(a['reflected_target_mca_bound']+a['reflected_target_mca_unbound']).sum()*29/len(a['x'])}[key]
 else:target_value=a['reflected_target_'+key][0]
 ax.axhline(target_value,color='#222222',ls=':',lw=1);ax.axvline(D,color='#D55E00',ls='--',lw=1);time_axis(ax,'x');ax.set_xlim(0,t[-1]);ax.tick_params(labelsize=8);ax.grid(alpha=.14);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
 ax.set_ylabel('dimensionless' if key=='mean_strain' else 'pN' if key=='tension_range' else 'conc. × µm' if key.startswith('total_') else 'µm')
fig.legend(handles=[plt.Line2D([],[],color='#0072B2',label='Native trajectory'),plt.Line2D([],[],color='#222222',ls=':',label='Reflected pre-cue target'),plt.Line2D([],[],color='#D55E00',ls='--',label='Cue withdrawal')],loc='outside lower center',ncol=3,frameon=False)
save_fields(fig,'L29_globals',GCAP,global_names)
base_dev=[['dev/L29_fullstate/supporting/viewers/L29_relative20_fullstate.html','Complete field-history explorer'],['dev/L29_fullstate/supporting/data/L29_relative20_fullstate_fields.npz','Every native field and reflected target'],['dev/L29_fullstate/supporting/metadata/L29_relative20_fullstate_audit.json','Field coverage and independent scalar-gate audit'],['dev/L29_fullstate/supporting/vectors/L29_relative20_fullstate_profiles.pdf','Initial, reflected and recovered profiles'],['dev/L29_fullstate/supporting/data/L29_relative20_fullstate_target_errors.csv','Secondary pointwise target comparisons'],['../../data/runs/length_mechanisms_20261008/quality_20261009_field_accompaniment/L29_endpoint_story.csv','Endpoint switch and mechanics measurements']]
overrides={
'mca_fields':dict(id='mca_fields',figure_label='8S1',title='A successful cue leaves a long mechanical and chemical recovery',image='quality_20261009_field_accompaniment/L29_fields.png',renderer='length_report_field_accompaniment_20261009.py',caption_override=True,caption=KCAP,paragraphs=[
 'The L = 29 µm relative-patch example supplies the fields behind the polarity traces. Its 282.755102-s cue is the shortest of 79 tested durations passing the canonical scalar 1% gate after 10000 s recovery; 280.229592 s fails.',
 'During cue application, Rac rises at the right while Rho falls (Fig. 8S1A,B). Bound and unbound MCA redistribute more slowly (Fig. 8S1C,D). The old left front’s protrusion switch falls below 0.001 by 70.5 s. At withdrawal, right Rac is 4.45; left boundary displacement is nearly zero, but left tension remains 8.86 pN. A stationary endpoint therefore does not imply an unstrained membrane (Fig. 8S1E–G).',
 'Withdrawal removes external activation without resetting any state (Fig. 8S1H–K). Right Rac subsequently falls through the 0.89906 endpoint switch threshold between 430 and 435 s from onset. Prescribed protrusion then weakens and endpoint relaxation releases much of the strain: maximum tension falls from 11.4 pN near 430 s to about 0.74 pN at 500 s and 0.071 pN at 750 s (Fig. 8S1E–G). Residual Rho gradients support small interior flow. Both ends remain weakly protrusive for a long interval before the new chemical and mechanical pole develops; scalar 1% arrival occurs at 8817.755 s. These observations trace the coupled model response, rather than isolate a single cause.'],sources=[['quality_20261009_field_accompaniment/L29_fields.pdf','Full-resolution field kymographs'],['quality_20261009_field_accompaniment/L29_fields.manifest.json','Native field legend, hashes and reproduction']]+base_dev),
'mca_globals':dict(id='mca_globals',figure_label='8S2',title='Endpoint relaxation is a prescribed mechanical rule',image='quality_20261009_field_accompaniment/L29_globals.png',renderer='length_report_field_accompaniment_20261009.py',caption_override=True,caption=GCAP,paragraphs=[
 'The endpoint law explains the rapid mechanical decline after the right front weakens. Its smooth Rac switch has midpoint 0.89906 native units and width 0.01. Each 0.5-s native step updates endpoint displacement as 0.98 times its previous value plus 0.5 times its local, tension-inhibited protrusion velocity. With protrusion suppressed, displacement relaxes on a 24.75-s timescale. Endpoint cortical velocity also follows the local protrusion command. This is prescribed relaxation, not a passive stress-free boundary.',
 'Raw and physical endpoint histories coincide because the two commands are independent (Fig. 8S2A–D). Separation, mean strain, tension range and physical length show the same contraction and late recovery (Fig. 8S2E–H). The mean displacement is retained, so neither translation removal nor recentering hides the endpoint response (Fig. 8S2I). Bound and unbound MCA amounts change slightly through their exchange, while their sum remains conserved (Fig. 8S2J–L).'],sources=[['quality_20261009_field_accompaniment/L29_globals.pdf','Full-resolution global trajectories'],['quality_20261009_field_accompaniment/L29_globals.manifest.json','Global definitions, hashes and reproduction']]+base_dev),
'flow':dict(image='quality_20261009_field_accompaniment/debelly_flow_AE.png',renderer='length_report_field_accompaniment_20261009.py',caption_override=True,caption=FLOW_CAPTION,sources=[['quality_20261009_field_accompaniment/debelly_flow_AE.pdf','Figure 9A–E vector plots'],['quality_20261009_field_accompaniment/debelly_flow_AE.manifest.json','Legend, hashes and reproduction']]),
'flow_budget':dict(id='flow_budget',figure_label='9F',title='Endpoint and Rho contributions depend on when the cue is sampled',image='dev/native_reversal/Fig9F_cue_timing.png',renderer='length_dev_fig9f_cue_timing.py',caption_override=True,caption=FCAP,paragraphs=[
 'At 60 s after onset, Rho-driven advection and endpoint advection both erase the old MCA bias across the sampled lengths; the grouped remaining terms oppose them (Fig. 9F-A). At 60 s before withdrawal, Rho remains positive, while the endpoint term varies in sign with length (Fig. 9F-B). The net rate remains positive in both panels. Thus endpoint forcing cannot be assigned one sign throughout the cue from its value at MCA zero.',
 'These are the shortest orientation-reversing constant-slope pulses used for the length comparison, not the relative-patch example above. Every selected pulse exceeds 60 s. The exact half-second source decomposition includes the changing MCA normalization; its contributions sum to the native rate. They are instantaneous source bookkeeping, rather than causal fractions of reversal time.'],sources=[['dev/native_reversal/supporting/vectors/Fig9F_cue_timing.pdf','Paired Figure 9F vector comparison'],['dev/native_reversal/supporting/metadata/Fig9F_cue_timing.manifest.json','Source identity, hashes and timing audit'],['dev/native_reversal/supporting/data/Fig9F_cue_timing.csv','Requested and actual native sampling times']])}
write(DATA/'chapter_overrides.json',overrides)
print('Generated report accompaniment overrides and Figure9A–E; no native integrations.')
