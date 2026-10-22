"""Six-panel length/advection figure from completed native trajectories."""
from pathlib import Path
import json,textwrap,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from debelly_advection import OUT,LENGTHS,STAGES,n,r,moment
FIG=n.ROOT/'figures/debelly_advection'
SC={'early':'#293A49','mid':'#BC7922','zero':'#8954A3'}
SL={'early':'50 s after cue onset','mid':'MCA polarity ≈ +0.05','zero':'MCA zero crossing'}
CC={'boundary':'#D55E00','rho':'#0072B2','membrane':'#009E73','diffusion':'#8F969C','exchange':'#CC79A7','replenishment':'#E6AB02','total':'#202020'}
CN={'boundary':'Boundary advection','rho':'Rho advection','membrane':'Membrane drag → cortex','diffusion':'Bound-MCA diffusion','exchange':'Bound/unbound exchange','replenishment':'Uniform replenishment','total':'Total'}

def readall():
 summaries=[];rows=[]
 for L in LENGTHS:
  a=r.read(OUT/f'L{L:g}_summary.json');assert n.sha(OUT/f'L{L:g}_native.npz')==a['native_sha256'];assert n.sha(OUT/f'L{L:g}_profiles.npz')==a['profiles_sha256'];summaries.append(a);rows+=a['rows']
 assert max(abs(a['budget_residual']) for a in rows)<1e-11 and max(a['mechanical_residual'] for a in rows)<1e-11
 return summaries,rows

def main():
 FIG.mkdir(parents=True,exist_ok=True);summaries,rows=readall();lookup={(a['length'],a['stage']):a for a in rows};p=n.model_for('debelly').parameters
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
 title(outer[0,1],'B   Early speed grows with length; near-zero speed converges')
 b=fig.add_subplot(outer[0,1]);stageplot(b,'centroid_speed');b.set_ylabel('Bound-MCA centroid speed toward new orientation (µm/s)');b.legend(frameon=False,fontsize=9,loc='upper left');b.set_ylim(bottom=0)
 rr=lookup[23.,'early'];b.plot(LENGTHS,[rr['centroid_speed']*L/23 for L in LENGTHS],color='#AAA',ls=':',lw=1.5,label='Reference ∝ L')
 b.text(.96,.14,'Early reference: speed ∝ L\nGray dotted line anchored at L = 23',transform=b.transAxes,ha='right',fontsize=9,color='#59656D')
 b.text(.96,.04,'Centroid speed = (L / 2) × reorientation rate',transform=b.transAxes,ha='right',fontsize=9)
 # C: exact native mechanical superposition, short/long and early/zero.
 title(outer[1,0],'C   Decompose cortical flow at the same native state')
 sub=outer[1,0].subgridspec(2,2,wspace=.2,hspace=.70)
 for row,stage in enumerate(['early','zero']):
  for col,L in enumerate([23.,39.]):
   ax=fig.add_subplot(sub[row,col]);z=np.load(OUT/f'L{L:g}_profiles.npz');x=z[stage+'_x_fraction']
   for component in ['boundary','rho','membrane','total']:ax.plot(x,z[stage+'_'+component],color=CC[component],lw=1.5,ls='--' if component=='total' else '-')
   ax.axhline(0,color='#AAA',lw=.6);ax.set_ylim(-.065,.065);ax.set_xlim(0,1);ax.set_title(f'L = {L:g} µm · '+('50 s' if stage=='early' else 'MCA zero'),fontsize=10);style(ax)
   if col==0:ax.set_ylabel('Cortical velocity (µm/s)')
   else:ax.tick_params(labelleft=False)
   ax.set_xlabel('Position / L (left → right)')
 pos=outer[1,0].get_position(fig);fig.legend([Line2D([0],[0],color=CC[k],ls='--' if k=='total' else '-') for k in ['boundary','rho','membrane','total']],['Boundary forcing','Rho forcing','Membrane drag','Total'],loc='upper left',bbox_to_anchor=(pos.x0,pos.y0-.034),ncol=4,frameon=False,fontsize=8,borderaxespad=0)
 # D: actual boundary speed and local tensions, distinguish saturated Rac vs tension.
 title(outer[1,1],'D   Longer cells sustain faster boundary motion')
 sub=outer[1,1].subgridspec(1,2,wspace=.32);d1=fig.add_subplot(sub[0]);d2=fig.add_subplot(sub[1]);stageplot(d1,'boundary_speed');d1.set_ylabel('Outward cortical boundary speed (µm/s)')
 for stage in STAGES:
  d2.plot(LENGTHS,[lookup[L,stage]['left_tension'] for L in LENGTHS],color=SC[stage],lw=1.6,ls='-',marker='o',ms=3)
  d2.plot(LENGTHS,[lookup[L,stage]['right_tension'] for L in LENGTHS],color=SC[stage],lw=1.6,ls='--',marker='s',ms=3)
 d2.set_ylabel('Endpoint tension (native units)');d2.set_xlabel('Cell length (µm)');style(d2);d2.set_title('Left: solid · Right: dashed',fontsize=9)
 d1.text(.96,.05,'Boundary speed stays larger\nat MCA zero crossing',transform=d1.transAxes,ha='right',fontsize=8)
 pos=outer[1,1].get_position(fig);fig.text(pos.x0,pos.y0-.038,'Rac protrusion switches are saturated; tension limits the imposed speed.',fontsize=9,color='#43525B')
 # E: actual gradient/flow trends and controlled mechanics-only template.
 title(outer[2,0],'E   A weaker Rho gradient can still drive faster flow')
 sub=outer[2,0].subgridspec(2,2,wspace=.28,hspace=.52,height_ratios=[1,1]);e1=fig.add_subplot(sub[0,0]);e2=fig.add_subplot(sub[0,1]);e3=fig.add_subplot(sub[1,0]);e4=fig.add_subplot(sub[1,1]);stageplot(e1,'rho_difference');e1.set_ylabel('Left–right Rho difference');e1.set_ylim(.3,.56);stageplot(e2,'stress_gradient_rms');e2.set_ylabel('Rho stress-gradient RMS\n(native force-density units)');stageplot(e3,'rho_velocity_rms');e3.set_ylabel('Rho-driven velocity RMS (µm/s)')
 grid=np.geomspace(3,1000,160);b0=.1*5/6.4;xi=p.linker_friction*b0;screen=np.sqrt(p.cortical_viscosity/xi);yy=[]
 for L in grid:
  N=201;force=np.ones(N)*(-p.active_stress*.5/L);v=n.ds._dirichlet_elliptic(np.ones(N)*xi,p.cortical_viscosity,force,L/N,0,0);yy.append(np.sqrt(np.mean(v*v)))
 yy=np.array(yy);e4.loglog(grid,yy,color='#42566A',lw=2);e4.loglog(grid,yy[0]*grid/grid[0],color='#AAA',ls=':',lw=1);e4.loglog(grid,yy[-1]*grid[-1]/grid,color='#AAA',ls='--',lw=1);e4.axvline(screen,color='#8954A3',ls=':',lw=1);e4.axvspan(23,39,color='#E69F00',alpha=.2);e4.set_ylim(1e-3,.04);e4.set_yticks([.001,.003,.01,.03]);e4.set_yticklabels(['0.001','0.003','0.01','0.03']);e4.minorticks_off();e4.set_xlim(3,1000);e4.set_xticks([10,100,1000]);e4.set_xticklabels(['10','100','1,000']);e4.set_xlabel('Length (µm), log scale');e4.set_ylabel('Controlled Rho-flow RMS (µm/s)');style(e4)
 e4.text(.04,.95,'Fixed ΔRho; uniform MCA\nNo boundary / membrane forcing',transform=e4.transAxes,va='top',fontsize=7.5);e4.text(.96,.06,f'Screening length ≈ {screen:.1f} µm\nOrange: sampled native lengths',transform=e4.transAxes,ha='right',fontsize=7)
 # F: exact bound-MCA native-step polarity budget; same donor selection throughout.
 title(outer[2,1],'F   Flow components contribute differently as MCA moves')
 sub=outer[2,1].subgridspec(1,2,wspace=.3);ff=[]
 for col,stage in enumerate(['early','zero']):
  ax=fig.add_subplot(sub[col]);ff.append(ax)
  for component in ['boundary','rho','membrane','diffusion','exchange','replenishment']:
   ax.plot(LENGTHS,[lookup[L,stage]['rate_'+component]*1e4 for L in LENGTHS],color=CC[component],lw=1.6,ls=':' if component in ['diffusion','membrane'] else '-',label=CN[component])
  ax.plot(LENGTHS,[lookup[L,stage]['reorientation_rate']*1e4 for L in LENGTHS],color=CC['total'],ls='--',lw=2,label='Total')
  ax.axhline(0,color='#999',lw=.8);ax.set_ylim(-3.6,7.3);ax.set_title('50 s after cue onset' if stage=='early' else 'At MCA zero crossing',fontsize=10);ax.set_xlabel('Cell length (µm)');style(ax)
  if col==0:ax.set_ylabel('Contribution to MCA reorientation\n(10⁻⁴ s⁻¹; positive helps reversal)')
  else:ax.tick_params(labelleft=False)
 ff[1].text(.04,.96,'Rho helps reorientation;\nboundary term opposes it here',transform=ff[1].transAxes,va='top',fontsize=8)
 pos=outer[2,1].get_position(fig);fig.legend(*ff[0].get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(pos.x0,pos.y0-.035),ncol=3,frameon=False,fontsize=8,borderaxespad=0)
 fig.suptitle('deBelly: mechanical advection and the length dependence of MCA reorientation',fontsize=20,weight='bold',y=.975)
 fig.text(.065,.944,'Early physical redistribution speeds rise with length; later convergence is not caused by boundary motion switching off.',fontsize=12,color='#43525B')
 caption=('Continuous rightward gradient, slope 0.24, from exact current basal-polarized checkpoints; native parameter regime, 201 cells, step 0.5 s. A: MCA polarity is the normalized bound-MCA first moment; reorientation rate is its negative derivative. B: centroid speed describes distribution movement, not individual-linker speed. Stages compare 50 s, first moment ≈ +0.05, and first recorded zero crossing. C–E: the linear cortical force-balance solve is decomposed with the same instantaneous MCA friction; boundary velocities, Rho stress gradient and membrane drag reconstruct the native velocity. Staggered forcing uses the preceding half-step. These are instantaneous contributions, not separately evolved causal interventions. E inset: a controlled linear Rho profile with fixed contrast 0.5, uniform basal MCA, and zero boundary/membrane forcing illustrates gradual viscous-to-friction crossover; gray guides show speed ∝ L and ∝ 1/L, not fits to native data. F: exact native-step bound-MCA polarity budget, including exchange and the solver’s uniform MCA mass replenishment. Component advection uses the total-flow upwind donors. Positive terms help reorientation; redistribution can change a term’s sign despite sustained flow. Cue peak increases with L under the fixed-slope protocol; this is not a fixed-peak scaling experiment.')
 fig.text(.065,.105,'\n'.join(textwrap.wrap(caption,190)),fontsize=9,va='top',linespacing=1.4)
 for ext in ['png','pdf','svg']:fig.savefig(FIG/f'debelly_advection.{ext}',dpi=200)
 plt.close(fig);n.csvwrite(OUT/'stage_measurements.csv',rows)
 (FIG/'caption.txt').write_text(caption+'\n')
 audit=dict(generated_at=n.now(),producer='finished',complete=True,lengths=LENGTHS,states=len(rows),max_velocity_reconstruction_residual=max(a['mechanical_residual'] for a in rows),max_polarity_budget_residual=max(abs(a['budget_residual']) for a in rows),figure_sha256={ext:n.sha(FIG/f'debelly_advection.{ext}') for ext in ['png','pdf','svg']},renderer_sha256=n.sha(Path(__file__)))
 r.write(OUT/'figure_audit.json',audit);print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
