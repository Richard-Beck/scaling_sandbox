"""Six-panel length/advection figure from completed native trajectories."""
from pathlib import Path
import json,textwrap,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from debelly_corrected_advection import OUT,LENGTHS,STAGES,n,r,moment
FIG=n.ROOT/'figures/deBelly_mechanism'
REV=n.DATA/'runs/debelly_corrected_advection_revision_20261008'
SC={'early':'#293A49','mid':'#BC7922','zero':'#8954A3'}
SL={'early':'50 s after cue onset','mid':'MCA polarity ≈ +0.05','zero':'MCA zero crossing'}
CC={'boundary':'#D55E00','rho':'#0072B2','membrane':'#009E73','diffusion':'#8F969C','exchange':'#CC79A7','replenishment':'#E6AB02','total':'#202020'}
CN={'boundary':'Prescribed-endpoint advection','rho':'Rho advection','membrane':'Membrane drag → cortex','diffusion':'Bound-MCA diffusion','exchange':'Bound/unbound exchange','replenishment':'Uniform replenishment','total':'Total'}

def readall():
 summaries=[];rows=[]
 for L in LENGTHS:
  a=r.read(OUT/f'L{L:g}_summary.json');assert n.sha(OUT/f'L{L:g}_native.npz')==a['native_sha256'];assert n.sha(OUT/f'L{L:g}_profiles.npz')==a['profiles_sha256'];summaries.append(a);rows+=a['rows']
 assert max(abs(a['budget_residual']) for a in rows)<1e-11 and max(a['mechanical_residual'] for a in rows)<1e-11
 return summaries,rows

def main():
 FIG.mkdir(parents=True,exist_ok=True);REV.mkdir(parents=True,exist_ok=True);assert r.read(OUT/'independent_audit.json')['passed'];
 for path,digest in r.read(OUT/'provenance.json')['sources'].items():assert n.sha(Path(path))==digest
 summaries,rows=readall();lookup={(a['length'],a['stage']):a for a in rows};p=n.model_for('debelly').parameters
 plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig=plt.figure(figsize=(18,16.8));outer=fig.add_gridspec(3,2,left=.065,right=.98,bottom=.19,top=.91,hspace=.42,wspace=.27)
 lengthcolors=plt.cm.viridis(np.linspace(.12,.9,len(LENGTHS)))
 def title(spec,text):
  pos=spec.get_position(fig);fig.text(pos.x0,pos.y1+.018,text,fontsize=13,weight='bold')
 def style(ax):ax.grid(alpha=.15);ax.tick_params(labelsize=8)
 def stageplot(ax,key,scale=1):
  for stage in STAGES:
   yy=[lookup[L,stage][key]*scale for L in LENGTHS];ax.plot(LENGTHS,yy,color=SC[stage],marker={'early':'o','mid':'s','zero':'^'}[stage],ms=4,lw=1.7,label=SL[stage])
  ax.set_xlabel('Cell length (µm)');style(ax)
 # A: native time traces and signed rate, sampled every .5 s.
 title(outer[0,0],'A   MCA memory reorients while the cue remains on')
 sub=outer[0,0].subgridspec(2,1,hspace=.23,height_ratios=[1,.8]);aa=[fig.add_subplot(sub[i]) for i in range(2)]
 for L,color,a in zip(LENGTHS,lengthcolors,summaries):
  with np.load(OUT/f'L{L:g}_native.npz') as z:
   t=z['time'];P=moment(z['field_mca_bound'],L,z['x']);rate=-np.gradient(P,t);mask=t<=800
   aa[0].plot(t[mask],P[mask],color=color,lw=1.6,label=f'{L:g}');aa[0].scatter(a['crossing_time'],0,color=color,s=16,zorder=3)
   aa[1].plot(t[mask],rate[mask]*1e4,color=color,lw=1.5)
 for ax in aa:ax.axhline(0,color='#777',lw=.8);ax.set_xlim(0,800);style(ax)
 aa[0].set_ylabel('Bound-MCA polarity');aa[0].legend(title='L (µm)',frameon=False,ncol=5,fontsize=8,loc='upper right');aa[0].tick_params(labelbottom=False)
 aa[1].set_ylabel('Reorientation rate\n(10⁻⁴ s⁻¹)');aa[1].set_xlabel('Time since cue onset (s)');aa[0].text(.03,.07,'Dots: MCA zero crossing',transform=aa[0].transAxes,fontsize=8)
 # B: finite native-step centroid speed, stage-matched length dependence.
 title(outer[0,1],'B   Length dependence of nominal MCA centroid speed')
 b=fig.add_subplot(outer[0,1]);stageplot(b,'centroid_speed');b.set_ylabel('Nominal bound-MCA centroid speed (µm/s)');b.legend(frameon=False,fontsize=9,loc='upper left');b.set_ylim(bottom=0)
 rr=lookup[23.,'early'];b.plot(LENGTHS,[rr['centroid_speed']*L/23 for L in LENGTHS],color='#AAA',ls=':',lw=1.5,label='Reference ∝ L')
 b.text(.96,.14,'Reference ∝ L\nGray dotted line anchored at L = 23',transform=b.transAxes,ha='right',fontsize=9,color='#59656D')
 b.text(.03,.74,'Nominal speed = (L / 2) × reorientation rate',transform=b.transAxes,ha='left',fontsize=9)
 # C: exact native mechanical superposition, short/long and early/zero.
 title(outer[1,0],'C   Decompose cortical flow at the same native state')
 sub=outer[1,0].subgridspec(2,2,wspace=.2,hspace=.70)
 for row,stage in enumerate(['early','zero']):
  for col,L in enumerate([23.,39.]):
   ax=fig.add_subplot(sub[row,col]);z=np.load(OUT/f'L{L:g}_profiles.npz');x=z[stage+'_x_fraction']
   for component in ['boundary','rho','membrane','total']:ax.plot(x,z[stage+'_'+component],color=CC[component],lw=1.5,ls='--' if component=='total' else '-')
   ax.axhline(0,color='#AAA',lw=.6);ax.set_ylim(-.17,.17);ax.set_xlim(0,1);ax.set_title(f'L = {L:g} µm · '+('50 s' if stage=='early' else 'MCA zero'),fontsize=10);style(ax)
   if col==0:ax.set_ylabel('Cortical velocity (µm/s)')
   else:ax.tick_params(labelleft=False)
   ax.set_xlabel('Position / L (left → right)')
 pos=outer[1,0].get_position(fig);fig.legend([Line2D([0],[0],color=CC[k],ls='--' if k=='total' else '-') for k in ['boundary','rho','membrane','total']],['Prescribed endpoints','Rho forcing','Membrane drag','Total'],loc='upper left',bbox_to_anchor=(pos.x0,pos.y0-.034),ncol=4,frameon=False,fontsize=8,borderaxespad=0)
 # D: retain boundary mechanism and show the native Rho transition in nominal model coordinates in micrometres.
 title(outer[1,1],'D   Prescribed endpoints and the localized Rho transition')
 sub=outer[1,1].subgridspec(2,2,wspace=.3,hspace=.72)
 d1=fig.add_subplot(sub[0,0]);d2=fig.add_subplot(sub[0,1]);
 for stage in STAGES:
  d1.plot(LENGTHS,[lookup[L,stage]['left_boundary_speed'] for L in LENGTHS],color=SC[stage],lw=1.5)
  d1.plot(LENGTHS,[lookup[L,stage]['right_boundary_speed'] for L in LENGTHS],color=SC[stage],lw=1.5,ls='--')
 d1.set_xlabel('Cell length (µm)');d1.set_ylabel('|Local endpoint command| (µm/s)');d1.set_title('Left: solid · Right: dashed',fontsize=8);style(d1)
 for stage in STAGES:
  d2.plot(LENGTHS,[lookup[L,stage]['left_tension'] for L in LENGTHS],color=SC[stage],lw=1.5)
  d2.plot(LENGTHS,[lookup[L,stage]['right_tension'] for L in LENGTHS],color=SC[stage],lw=1.5,ls='--')
 d2.set_ylabel('Endpoint tension');d2.set_xlabel('Cell length (µm)');d2.set_title('Left: solid · Right: dashed',fontsize=8);style(d2)
 rawprofiles={}
 for L in LENGTHS:
  with np.load(OUT/f'L{L:g}_native.npz') as z:
   for stage in ['early','zero']:
    i=int(np.searchsorted(z['time'],lookup[L,stage]['time']));rawprofiles[L,stage]=(z['x'].copy(),z['field_rho'][i-1].copy())
 for col,stage in enumerate(['early','zero']):
  ax=fig.add_subplot(sub[1,col]);widths=[]
  for L in [23.,31.,39.]:
   x,rho=rawprofiles[L,stage];color=lengthcolors[LENGTHS.index(L)];ax.plot(x,rho,color=color,lw=1.8,label=f'{L:g} µm')
   level=(rho[0]+rho[-1])/2;j=np.flatnonzero(rho<=level)[0];w=float(np.interp(level,rho[j-1:j+1][::-1],x[j-1:j+1][::-1]));widths.append(w);ax.scatter(w,level,color=color,s=18,zorder=3)
  ax.set_xlim(0,39);ax.set_ylim(0,.62);ax.set_title('Rho: '+('50 s' if stage=='early' else 'at MCA zero'),fontsize=9);ax.set_xlabel('Reference distance from left end (µm)');style(ax)
  if col==0:ax.set_ylabel('Rho concentration');ax.legend(frameon=False,fontsize=7,ncol=3,loc='upper right')
  ax.text(.96,.06,f'Half-drop: {widths[0]:.2f} → {widths[-1]:.2f} µm',transform=ax.transAxes,ha='right',fontsize=7.5)
 pos=outer[1,1].get_position(fig);fig.text(pos.x0,pos.y0-.035,'Independent local endpoint commands; dots mark the reference-grid Rho half-drop.',fontsize=8.5,color='#43525B')
 # E: hold the physical Rho profile fixed, extend only its low-Rho region.
 title(outer[2,0],'E   Same local Rho forcing; a longer cortex moves faster')
 sub=outer[2,0].subgridspec(1,2,wspace=.32);e1=fig.add_subplot(sub[0]);e2=fig.add_subplot(sub[1]);x0,r0=rawprofiles[23.,'zero'];b0=.1*5/6.4;xi=p.linker_friction*b0;diagnostics=[]
 def controlled(L,N=2001,friction=True):
  dx=L/N;x=(np.arange(N)+.5)*dx;rho=np.interp(x,x0,r0,left=r0[0],right=r0[-1]);force=p.active_stress*n.ds._gradient(rho,dx)
  v=n.ds._dirichlet_elliptic(np.full(N,xi if friction else 0.),p.cortical_viscosity,force,dx,0.,0.)
  return x,rho,-v
 for L,color in zip([23.,39.,60.,100.],['#482878','#35B779','#E69F00','#B55D60']):
  x,rho,v=controlled(L);e1.plot(np.r_[0,x,L],np.r_[0,v,0],color=color,lw=1.8,label=f'{L:g} µm');e1.scatter(L,0,color=color,s=18,zorder=3)
 e1.axvspan(0,23,color='#B4BEC5',alpha=.17);e1.axvline(5,color='#777',ls=':',lw=1);e1.set_xlim(0,100);e1.set_ylim(0,.019);e1.set_xlabel('Reference distance from left end (µm)');e1.set_ylabel('Rho-only speed toward left (µm/s)');e1.set_title('Fixed reference-grid Rho profile; basal MCA friction',fontsize=8.5);e1.legend(frameon=False,fontsize=8,ncol=2,loc='upper left');style(e1)
 e1.text(.04,.04,'Same forcing region ≤23 µm\nVelocity = 0 at both endpoints',transform=e1.transAxes,fontsize=8)
 inset=e1.inset_axes([.62,.60,.35,.23]);inset.plot(np.r_[0,x0,100],np.r_[r0[0],r0,r0[-1]],color='#555',lw=1.2);inset.set_xlim(0,100);inset.set_ylim(0,.6);inset.set_xticks([0,23,100]);inset.set_yticks([0,.5]);inset.tick_params(labelsize=6);inset.set_title('Rho profile held fixed',fontsize=7)
 grid=np.linspace(23,150,90);curves={}
 for friction in [True,False]:
  yy=[]
  for L in grid:
   x,rho,v=controlled(float(L),friction=friction);speed=float(np.interp(5,x,v));yy.append(speed);diagnostics.append(dict(length=float(L),friction=friction,speed_at_x5=speed,peak_speed=float(v.max()),rms_speed=float(np.sqrt(np.mean(v*v)))))
  curves[friction]=np.array(yy);e2.plot(grid,yy,color='#0072B2' if friction else '#777',ls='-' if friction else '--',lw=2,label='With MCA friction' if friction else 'Without friction')
 native=[]
 for L in LENGTHS:
  z=np.load(OUT/f'L{L:g}_profiles.npz');v=-z['zero_rho'];x=z['zero_x_fraction']*L;native.append(float(np.interp(5,x,v)))
 e2.scatter(LENGTHS,native,color='#202020',s=23,zorder=4,label='Native, at MCA zero');e2.axvspan(23,39,color='#E69F00',alpha=.1);e2.set_xlim(23,150);e2.set_ylim(bottom=0);e2.set_xlabel('Cell length (µm)');e2.set_ylabel('Rho-only speed at x = 5 µm (µm/s)');e2.set_title('Longer return distance reduces viscous resistance',fontsize=8.5);e2.legend(frameon=False,fontsize=8,loc='lower right');style(e2)
 e2.text(.04,.50,'Same Rho force throughout\nZero protrusion / membrane forcing',transform=e2.transAxes,va='top',fontsize=8)
 e2.text(.04,.34,f'MCA screening length ≈ {np.sqrt(p.cortical_viscosity/xi):.1f} µm',transform=e2.transAxes,fontsize=8)
 checks=[]
 for L in [23.,39.,100.,150.]:
  for friction in [True,False]:
   x,rho,v=controlled(L,friction=friction);xf,rf,vf=controlled(L,N=4001,friction=friction);rel=abs(np.interp(5,x,v)/np.interp(5,xf,vf)-1);assert rel<.001;checks.append(dict(length=L,friction=friction,relative_mesh_change=float(rel)))
 n.csvwrite(REV/'fixed_reference_rho_extension.csv',diagnostics)
 n.csvwrite(REV/'native_fixed_position_speeds.csv',[dict(length=L,speed_at_x5=v) for L,v in zip(LENGTHS,native)])
 r.write(REV/'controlled_calculation_audit.json',dict(generated_at=n.now(),passed=True,mesh_checks=checks,source_native_sha256=n.sha(OUT/'L23_native.npz'),rho_template='L23 at MCA zero on native nominal coordinates; extend terminal low-Rho concentration without stretching profile',coordinate_convention='reference grid in micrometres; not displacement-mapped laboratory coordinates',bound_mca=b0,fixed_position=5.,cells=2001,zero_boundary_and_membrane_forcing=True))
 # F: exact bound-MCA native-step polarity budget; same donor selection throughout.
 title(outer[2,1],'F   Flow components contribute differently as MCA moves')
 sub=outer[2,1].subgridspec(1,2,wspace=.3);ff=[]
 for col,stage in enumerate(['early','zero']):
  ax=fig.add_subplot(sub[col]);ff.append(ax)
  for component in ['boundary','rho']:
   ax.plot(LENGTHS,[lookup[L,stage]['rate_'+component]*1e4 for L in LENGTHS],color=CC[component],lw=1.6,ls=':' if component in ['diffusion','membrane'] else '-',label=CN[component])
  remainder=[sum(lookup[L,stage]['rate_'+k] for k in ['membrane','diffusion','exchange','replenishment'])*1e4 for L in LENGTHS]
  ax.plot(LENGTHS,remainder,color='#90969C',lw=1.6,label='Other terms combined')
  ax.plot(LENGTHS,[lookup[L,stage]['reorientation_rate']*1e4 for L in LENGTHS],color=CC['total'],ls='--',lw=2,label='Total')
  ax.axhline(0,color='#999',lw=.8);ax.set_title('50 s after cue onset' if stage=='early' else 'At MCA zero crossing',fontsize=10);ax.set_xlabel('Cell length (µm)');style(ax)
  if col==0:ax.set_ylabel('Contribution to MCA reorientation\n(10⁻⁴ s⁻¹; positive helps reversal)')
  else:ax.tick_params(labelleft=False)
 shared_limits=(min(ax.get_ylim()[0] for ax in ff),max(ax.get_ylim()[1] for ax in ff))
 for ax in ff:ax.set_ylim(*shared_limits)
 pos=outer[2,1].get_position(fig);fig.legend(*ff[0].get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(pos.x0,pos.y0-.035),ncol=2,frameon=False,fontsize=8,borderaxespad=0)
 fig.suptitle('deBelly: how cell length changes mechanical flow and MCA reorientation',fontsize=20,weight='bold',y=.975)
 fig.text(.065,.944,'Independent endpoint commands and direct Rho stress drive flow; a fixed-profile control isolates mechanical domain-length effects.',fontsize=12,color='#43525B')
 caption=('A–B retain native MCA polarity, its negative derivative and nominal-coordinate centroid speed under the continuous gradient (slope 0.24; peak increases with length). Corrected reflection-symmetric extension: no in-loop centering; inactive endpoints remain constrained. All nine cells start from fresh corrected relative-20% patch/basal checkpoints; native regime, 201 cells, step 0.5 s. C: boundary, Rho and membrane-drag solves reconstruct the same native cortical velocity; positive velocity is rightward. D: each endpoint has its own Rac/tension command (left solid, right dashed); Rho half-drop markers show the measured transition width in nominal micrometres. E: preserve the L=23 Rho profile at MCA zero in nominal model coordinates in micrometres and extend only its terminal low-Rho region. Uniform basal MCA friction, zero boundary/membrane forcing; the Rho-only component vanishes at both endpoints. A farther endpoint allows a more gradual velocity return and lower viscous resistance; MCA friction limits the range of this coupling (screening length shown). Curves use 2001 cells; speed at the fixed reference position x=5 avoids changing-domain RMS effects. Native markers are instantaneous Rho components at MCA zero, not separately evolved interventions. F: exact MCA-polarity budget; other terms combine membrane-drag advection, bound diffusion, exchange and uniform MCA replenishment. The complete term budget remains available with the original figure and stage CSV. This is instantaneous source superposition at fixed state, not a full causal pathway partition. Decomposition uses native staggered forcing and total-flow upwind donors; contribution signs depend on the MCA distribution as well as flow speed.')
 fig.text(.065,.105,'\n'.join(textwrap.wrap(caption,190)),fontsize=9,va='top',linespacing=1.4)
 for ext in ['png','pdf','svg']:fig.savefig(FIG/f'debelly_advection.{ext}',dpi=200)
 plt.close(fig);n.csvwrite(REV/'stage_measurements.csv',rows)
 (FIG/'debelly_advection_caption.txt').write_text(caption+'\n')
 audit=dict(generated_at=n.now(),producer='finished',complete=True,lengths=LENGTHS,states=len(rows),max_velocity_reconstruction_residual=max(a['mechanical_residual'] for a in rows),max_polarity_budget_residual=max(abs(a['budget_residual']) for a in rows),figure_sha256={ext:n.sha(FIG/f'debelly_advection.{ext}') for ext in ['png','pdf','svg']},renderer_sha256=n.sha(Path(__file__)))
 audit['native_independent_audit_sha256']=n.sha(OUT/'independent_audit.json');audit['controlled_audit_sha256']=n.sha(REV/'controlled_calculation_audit.json');r.write(REV/'figure_audit.json',audit);print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
