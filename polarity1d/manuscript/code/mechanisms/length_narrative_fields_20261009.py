"""Narrative-ready complete fields and compact endpoint diagnostics; no integration."""
from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'code'), str(ROOT/'code/core')]
import debelly_corrected_reversal as corrected
n = corrected.n
FIG = ROOT/'figures/length_mechanisms_20261008'
OUT = FIG/'narrative_20261009/fields'
DATA = ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/fields'
OUT.mkdir(parents=True, exist_ok=True); DATA.mkdir(parents=True, exist_ok=True)
SOURCE = FIG/'dev/L29_fullstate/supporting/data/L29_relative20_fullstate_fields.npz'
AUDIT = FIG/'dev/L29_fullstate/supporting/metadata/L29_relative20_fullstate_audit.json'
SOURCE_MANIFEST = FIG/'dev/L29_fullstate/supporting/metadata/L29_relative20_fullstate.manifest.json'
SELECT = FIG/'dev/L29_fullstate/supporting/data/L29_relative20_fullstate_duration_selection.csv'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p, obj): p.write_text(json.dumps(obj, indent=2)+'\n')
def csvwrite(p, rows):
    with p.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with np.load(SOURCE) as z: a={k:z[k].copy() for k in z.files}
aud = json.loads(AUDIT.read_text()); assert aud['passed'] and aud['canonical_scalar_gate_pass']
selection=list(csv.DictReader(SELECT.open()))
passing=[float(row['duration']) for row in selection if row['canonical_full_polarity_success']=='True']
t=a['time'];x=a['x'];D=float(a['cue_duration']);assert D==min(passing)
p=n.model_for('debelly').parameters; assert not p.center_boundary_motion
rac0,_=n.ds._homogeneous_rac_rho(p,.5); threshold=1.3*rac0
rac=a['field_rac'][:,[0,-1]]; tension=a['field_tension'][:,[0,-1]]
switch=n.ds.smooth_threshold(rac,p.protrusion_switch_scale,threshold)
command=p.rac_to_rho_coupling*p.protrusion_speed*switch/(1+tension**2/p.protrusion_tension_scale)
command[:,1]*=-1
u=np.stack([a['global_left_boundary_displacement'],a['global_right_boundary_displacement']],axis=1)
# Stored cortical endpoint velocities are exactly the independently prescribed commands.
# The cell-centred velocity field is not substituted for an endpoint command.
idx=np.flatnonzero((t>D)&(rac[:,1]<threshold))[0]
Tswitch=float(np.interp(threshold,rac[idx-1:idx+1,1][::-1],t[idx-1:idx+1][::-1]))
orientation=n.db._orientation(a['field_rac']); strength=abs(orientation[0])
within=(orientation>.05)&(abs(orientation/strength-1)<=.01)
last_bad=int(np.flatnonzero(~within)[-1]); arrival=float(t[last_bad+1])
assert arrival>=D and within[last_bad+1:].all()
assert np.allclose(a['global_raw_left_displacement'],u[:,0]) and np.allclose(a['global_raw_right_displacement'],u[:,1])
# Check the exact 0.5-s endpoint update on contiguous cue samples, independently.
regular=np.flatnonzero(np.isclose(np.diff(t),.5))
errors=u[regular+1]-(p.boundary_relaxation*u[regular]+.5*command[regular])
np.testing.assert_allclose(errors,0,atol=2e-13,rtol=1e-12)
relaxation=-.5/np.log(p.boundary_relaxation)
MARKS=[(0.,'Pre-cue','#666666','-'),(D,'Withdrawal','#D55E00','--'),(Tswitch,'Right switch decline','#0072B2',':'),(arrival,'Recovered polarity','#762A83','-.')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':220})
SOURCES={str(z.relative_to(ROOT)):sha(z) for z in [Path(__file__),SOURCE,AUDIT,SOURCE_MANIFEST,SELECT,Path(n.ds.__file__),Path(corrected.__file__)]}
COMMON='Exact corrected independent-endpoint history at L = 29 µm: right cosine patch, support 5.8 µm (20% of length), peak 3.6, withdrawn at 282.755102 s, followed by 10000 s at basal input. This is the shortest of 79 tested pulses passing the canonical 1%-strength gate; 280.229592 s fails. Native step 0.5 s; stored cue/recovery intervals 0.5/2.5 s. No fields are reset or centered at withdrawal.'
def save(fig,stem,caption,panels,methods):
    exports={}
    for ext in ['pdf','svg','png']:
        path=OUT/(stem+'.'+ext);fig.savefig(path,facecolor='white');exports[path.name]=sha(path)
    plt.close(fig)
    manifest=dict(passed=True,producer='finished',generated_at=datetime.now(timezone.utc).isoformat(),legend=caption,sha256=exports,figure_sha256=exports,source_sha256=SOURCES,renderer_sha256=sha(Path(__file__)),panels=panels,methods=methods,reproduction=['OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_narrative_fields_20261009.py'],data_sha256={str(z.relative_to(ROOT)):sha(z) for z in sorted(DATA.glob('*')) if z.suffix in ['.csv','.npz']})
    write(OUT/(stem+'.manifest.json'),manifest);(OUT/(stem+'_caption.txt')).write_text(caption+'\n')
def panel(ax,label,title):
    ax.set_title(title,loc='left',pad=8)
    ax.text(-.11,1.035,label,transform=ax.transAxes,fontweight='bold',fontsize=12)
def hybrid(ax,direction):
    getattr(ax,'set_'+direction+'scale')('symlog',linthresh=500,linscale=.5)
    getattr(ax,'set_'+direction+'ticks')([0,500,1000,3000,10000],['0','500','1000','3000','10000'])
def marks(ax,direction,limit):
    for when,_,color,style in MARKS[1:]:
        if when<=limit:
            getattr(ax,'ax'+('hline' if direction=='y' else 'vline'))(when,color=color,ls=style,lw=1.1,alpha=.9)
SPATIAL=['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','rac_activation','rho_activation','delta_alpha_mca','delta_beta_tension']
TITLES=['Rac','Rho','Bound MCA','Unbound MCA','Membrane displacement','Cortical flow','Membrane tension','Rac activation','Rho activation','MCA → Rac activation','Tension → Rho activation']
UNITS=['native conc.']*4+['µm','µm/s','pN']+['activation units']*4
snap_indices=[0,int(np.argmin(abs(t-D))),int(np.argmin(abs(t-Tswitch))),int(np.argmin(abs(t-arrival)))]
rows=[]
for requested,label,_,_ in MARKS:
    i=int(np.argmin(abs(t-requested)))
    row=dict(event=label,requested_time=requested,actual_time=float(t[i]),time_after_withdrawal=float(t[i]-D),Rac_orientation=float(orientation[i]))
    for side,j in [('left',0),('right',1)]:
        row.update({side+'_Rac':float(rac[i,j]),side+'_Rac_over_threshold':float(rac[i,j]/threshold),side+'_switch':float(switch[i,j]),side+'_protrusion_command':float(command[i,j]),side+'_displacement':float(u[i,j]),side+'_tension':float(tension[i,j])})
    row['maximum_tension']=float(a['field_tension'][i].max());rows.append(row)
csvwrite(DATA/'event_measurements.csv',rows)
csvwrite(DATA/'endpoint_diagnostics.csv',[dict(time=float(tt),**{side+'_'+key:float(values[i,j]) for side,j in [('left',0),('right',1)] for key,values in [('Rac',rac),('Rac_over_threshold',rac/threshold),('switch',switch),('protrusion_command',command),('displacement',u),('tension',tension)]}) for i,tt in enumerate(t)])
np.savez_compressed(DATA/'selected_spatial_snapshots.npz',x=x,event_time=t[snap_indices],**{key:a['field_'+key][snap_indices] for key in SPATIAL})
fig,axs=plt.subplots(3,4,figsize=(12.8,9.0),layout='constrained')
for i,(key,title,unit) in enumerate(zip(SPATIAL,TITLES,UNITS)):
    ax=axs.flat[i];values=a['field_'+key]
    if key in ['displacement','velocity','tension']:
        vlim=float(abs(values).max());norm=Normalize(-vlim,vlim);cmap='RdBu_r'
    else:norm=Normalize(float(values.min()),float(values.max()));cmap='viridis'
    mesh=ax.pcolormesh(x/29,t,values,shading='nearest',rasterized=True,norm=norm,cmap=cmap)
    panel(ax,chr(65+i),title);hybrid(ax,'y');marks(ax,'y',t[-1]);ax.set(xlim=(0,1),ylim=(0,t[-1]),xticks=[0,.5,1],xlabel='Position x/L')
    if i%4==0:ax.set_ylabel('Time from cue onset (s)')
    ax.tick_params(labelsize=9)
    cb=fig.colorbar(mesh,ax=ax,pad=.025,shrink=.9,aspect=24);cb.set_label(unit,fontsize=8);cb.ax.tick_params(labelsize=8)
    cb.ax.yaxis.get_offset_text().set_fontsize(8)
ax=axs.flat[-1];ax.axis('off')
ax.text(.04,.94,'One continuous history',transform=ax.transAxes,fontweight='bold',fontsize=12,va='top')
for j,(when,label,color,style) in enumerate(MARKS):
    yy=.77-j*.145;ax.plot([.05,.22],[yy,yy],transform=ax.transAxes,color=color,ls=style,lw=1.8)
    ax.text(.27,yy,label+'\n'+('0 s' if when==0 else f'{when:.1f} s'),transform=ax.transAxes,va='center',fontsize=10)
ax.text(.04,.10,'Linear time to 500 s;\nlogarithmic thereafter.\nAbsolute native field ranges.',transform=ax.transAxes,fontsize=9,va='top',color='#444444')
KCAP='Complete spatial history of a successful pulse. '+COMMON+' A–D show Rac, Rho and bound/unbound MCA; E–G show absolute membrane displacement, cortical velocity and tension; H–K show total activation rates and their MCA/tension feedback contributions. Every one of the eleven recorded spatial fields is retained. Coordinates are the native reference position x/L. Time is linear through 500 s and logarithmic thereafter. Common lines mark withdrawal, the right endpoint Rac switch midpoint (433.016 s, interpolated between retained samples), and the first persistent scalar 1%-strength arrival (8817.755 s); the pre-cue chemical and mechanical state lies at the lower edge. Activation maps begin with the cue-on rate at onset. The switch crossing is displayed at its continuous interpolated time; all field values remain native samples. Activation at the withdrawal sample uses basal input; both cue-on and cue-off rate limits are retained in the native archive. The full recovery is displayed; scalar arrival does not impose an extra field-convergence gate.'
save(fig,'L29_fields',KCAP,SPATIAL,'All retained spatial samples; no interpolation of field values; common hybrid symlog time (linear to 500 s).')
# Compact endpoint mechanism with early linear and complete recovery views.
fig,axs=plt.subplots(5,2,figsize=(12.2,9.0),sharex='col',layout='constrained')
series=[rac/threshold,switch,command,u,tension]
ylabels=['Rac / threshold','Switch activation','Command (µm/s)','Displacement (µm)','Local tension (pN)']
colors=['#D55E00','#0072B2'];side_labels=['Left (old front)','Right (challenged end)']
for r,(values,ylabel) in enumerate(zip(series,ylabels)):
    for c in range(2):
        ax=axs[r,c];limit=800 if c==0 else t[-1]
        ax.axvspan(0,D,color='#E8E8E8',alpha=.8,zorder=0)
        for j in range(2):
            ax.plot(t,values[:,j],color=colors[j],lw=1.65,label=side_labels[j])
            ax.plot(0,values[0,j],marker='o',ms=3,color=colors[j],clip_on=False,zorder=5)
        marks(ax,'x',limit);ax.set_xlim(0,limit);ax.grid(axis='y',alpha=.16);ax.set_ylabel(ylabel)
        if r in [2,3,4]:ax.axhline(0,color='#999999',lw=.7,zorder=0)
        if r==0:ax.axhline(1,color='#555555',ls='--',lw=.9)
        if r==1:ax.set_ylim(-.04,1.04)
        if c==1:hybrid(ax,'x')
        else:ax.set_xticks([0,200,400,600,800])
        ax.text(.985,.94,chr(65+2*r+c),transform=ax.transAxes,fontweight='bold',fontsize=11,ha='right',va='top',bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.6))
axs[0,0].set_title('Cue and mechanical relaxation',loc='left',fontweight='bold',pad=13)
axs[0,1].set_title('Complete basal recovery',loc='left',fontweight='bold',pad=13)
for ax in axs[-1]:ax.set_xlabel('Time from cue onset (s)')
handles=[Line2D([],[],color=color,lw=2,label=label) for color,label in zip(colors,side_labels)]
handles.extend([Line2D([],[],color='#666666',lw=6,alpha=.2,label='Cue on'),Line2D([],[],color=MARKS[2][2],ls=MARKS[2][3],label='Right switch decline'),Line2D([],[],color=MARKS[3][2],ls=MARKS[3][3],label='Recovered polarity')])
fig.legend(handles=handles,loc='outside lower center',ncol=5,frameon=False,fontsize=9)
GCAP='The endpoint switch links chemical decline to fast mechanical relaxation. '+COMMON+' Rows compare both ends: Rac divided by its protrusion threshold (0.899062 native units), smooth switch activation, the actual signed tension-inhibited protrusion command, absolute endpoint displacement, and local membrane tension. Left-column time is linear; right-column time is linear to 500 s and logarithmic thereafter. Gray shading marks cue application; the vertical orange line is withdrawal, blue dotted line the right switch midpoint, and purple line recovered scalar polarity. Circles mark the exact pre-cue state. Commands are recomputed from native Rac and tension: left command = −0.15 × switch / (1 + tension²/10), with the opposite sign at right; they are the prescribed cortical boundary velocities, not the cell-centred flow or the net displacement velocity. Each native 0.5-s step updates each endpoint independently as displacement(new) = 0.98 × displacement(old) + 0.5 × command. When the switch closes, the inherited displacement therefore relaxes with a 24.749-s exponential timescale (approximately 25 s). The old endpoint can become almost stationary while local tension persists. Raw and physical endpoints coincide in this corrected history. All seven uniform native fields and all global derived readouts remain in the linked complete native archive.'
save(fig,'L29_endpoints',GCAP,['Rac/threshold','switch','actual protrusion command','endpoint displacement','local tension'],'Exact native endpoint states, native smooth threshold and signed tension-inhibited commands. Endpoint updates independently checked on all contiguous half-second output samples.')
# Linked profile companion makes the four common events spatially explicit.
fig,axs=plt.subplots(2,3,figsize=(12.2,5.5),layout='constrained')
for i,key in enumerate(['rac','rho','mca_bound','displacement','velocity','tension']):
    ax=axs.flat[i];panel(ax,chr(65+i),TITLES[SPATIAL.index(key)])
    for j,idx in enumerate(snap_indices):
        _,label,color,style=MARKS[j];ax.plot(x/29,a['field_'+key][idx],color=color,ls=style,lw=1.65,label=label)
    if key=='rac':ax.axhline(threshold,color='#999999',lw=.8,ls='--')
    ax.set(xlim=(0,1),xlabel='Position x/L',ylabel=UNITS[SPATIAL.index(key)],xticks=[0,.5,1]);ax.grid(alpha=.13)
fig.legend(handles=[Line2D([],[],color=color,ls=style,label=label) for _,label,color,style in MARKS],loc='outside lower center',ncol=4,frameon=False)
SCAP='Spatial states accompanying the common history markers. '+COMMON+' The four states are exact retained samples at 0, 282.755102, 432.755102 and 8817.755102 s (the nearest sample to the interpolated right-switch midpoint is 0.261288 s earlier). A–F: Rac, Rho, bound MCA, absolute membrane displacement, cortical velocity and membrane tension. The Rac threshold is the gray dashed line. Coordinates and concentrations retain their native values; no state is reflected or recentered.'
save(fig,'L29_snapshots',SCAP,['rac','rho','mca_bound','displacement','velocity','tension'],'Exact nearest samples for four marker times; actual times in event_measurements.csv.')
write(DATA/'audit.json',dict(passed=True,source_sha256=sha(SOURCE),selected_duration=D,shortest_strength_passing_duration=min(passing),tested_durations=len(selection),full_recovery=True,eleven_spatial_fields=SPATIAL,complete_global_native_fields_retained=True,endpoint_switch_threshold=threshold,right_switch_midpoint_interpolated=Tswitch,first_persistent_scalar_1pct_arrival=arrival,endpoint_relaxation_timescale=relaxation,half_second_update_checked_samples=len(regular),max_endpoint_update_error=float(abs(errors).max()),no_new_integrations=True,exact_marker_samples=rows,source_manifest=str(SOURCE_MANIFEST.relative_to(ROOT))))
base_sources=[['dev/L29_fullstate/supporting/data/L29_relative20_fullstate_fields.npz','Complete native spatial and global history'],['dev/L29_fullstate/supporting/metadata/L29_relative20_fullstate.manifest.json','Native history selection, field audit and provenance'],['../../data/runs/length_mechanisms_20261008/narrative_20261009/fields/event_measurements.csv','Marker states and actual sampling times'],['../../data/runs/length_mechanisms_20261008/narrative_20261009/fields/endpoint_diagnostics.csv','Both endpoint diagnostics throughout the complete history'],['narrative_20261009/fields/L29_snapshots.pdf','Spatial snapshots at the four history markers']]
overrides={
'mca_fields':dict(title='A successful pulse carries mechanical memory into recovery',image='narrative_20261009/fields/L29_fields.png',renderer=Path(__file__).name,caption_override=True,caption=KCAP,paragraphs=['A successful pulse first raises Rac at the challenged right end and redistributes its MCA support. The old left front weakens early: its protrusion switch is below 0.001 by 70.5 s, and its displacement is nearly zero at withdrawal, even though local tension remains 8.86 pN. Membrane tension is spatially distributed, so a nearly stationary endpoint can retain substantial mechanical loading.','Withdrawal removes activation while preserving the inherited fields. Right Rac falls through its protrusion-switch midpoint at 433.0 s from onset; the endpoint command then declines and membrane strain relaxes rapidly. Much later, the reversed chemical pole and its mechanical state rebuild together. The complete history separates fast mechanical relaxation from the slower recovery of a strong reversed pole.'],sources=[['narrative_20261009/fields/L29_fields.pdf','Complete spatial history as vector PDF'],['narrative_20261009/fields/L29_fields.manifest.json','Legend, hashes and reproduction']]+base_sources),
'mca_globals':dict(title='The local endpoint switch releases inherited displacement',image='narrative_20261009/fields/L29_endpoints.png',renderer=Path(__file__).name,caption_override=True,caption=GCAP,paragraphs=['The endpoint diagnostic connects the Rac switch to the actual command that moves each end. Tension attenuates that command even when the switch is fully open. Once Rac falls below the switch threshold, the command collapses and the prescribed endpoint rule relaxes displacement on an approximately 25-s timescale.','This mechanical response precedes the late rebuilding of the reversed pole. The old end and challenged end follow their own local commands throughout. The two local histories show how mechanical relaxation can be fast while attachment and chemical polarity continue to recover.'],sources=[['narrative_20261009/fields/L29_endpoints.pdf','Compact endpoint diagnostic as vector PDF'],['narrative_20261009/fields/L29_endpoints.manifest.json','Command definitions, endpoint update audit and reproduction']]+base_sources)}
write(DATA/'chapter_overrides.json',overrides)
print(json.dumps(dict(passed=True,output=str(OUT),duration=D,switch_crossing=Tswitch,arrival=arrival,endpoint_update_max_error=float(abs(errors).max())),indent=2))
