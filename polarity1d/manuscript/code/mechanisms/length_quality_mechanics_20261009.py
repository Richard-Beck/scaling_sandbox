"""Native, deformed spring/dashpot snapshots and declared screening controls.
No integrations of the mechanochemical model; all snapshots are corrected native outputs.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,csv,sys,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.linalg import solve_banded
ROOT=Path(__file__).resolve().parents[2]
FIG=ROOT/'figures/length_mechanisms_20261008/quality_20261009_mechanics'
DATA=ROOT/'data/runs/length_mechanisms_20261008/quality_20261009_mechanics'
RUN=ROOT/'data/runs/length_debelly_review_20261009'
FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
BLUE='#0072B2';RED='#D55E00';PURPLE='#762a83';GREY='#646464'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
sources={}
def register(p):sources[os.path.relpath(p,ROOT)]=sha(p)
def panel(ax,label,title):
 ax.text(-.085,1.035,label,transform=ax.transAxes,fontfamily='DejaVu Sans',fontsize=12,fontweight='bold',va='bottom');ax.set_title(title,loc='left',pad=9)
def spring(ax,a,b,y,color):
 x=np.r_[a,np.linspace(a+.15*(b-a),b-.15*(b-a),11),b];yy=np.full_like(x,y);yy[2:-2]+=.075*np.resize([1,-1],len(x)-4);ax.plot(x,yy,color=color,lw=1.6)
def dashpot(ax,p0,p1,color,lw=1.1,size=.07):
 p0,p1=np.asarray(p0),np.asarray(p1);d=p1-p0;unit=d/np.linalg.norm(d);n=np.array([-unit[1],unit[0]]);b=p0+.36*d;c=p0+.68*d
 # Open cylinder and a separate piston; tilt and piston offset follow the slip geometry.
 ax.plot(*np.array([p0,b]).T,color=color,lw=lw)
 ax.plot(*np.array([b+size*n,c+size*n,c-size*n,b-size*n]).T,color=color,lw=lw)
 piston=p0+.51*d;ax.plot(*np.array([piston-size*n,piston+size*n]).T,color=color,lw=lw)
 ax.plot(*np.array([piston,p1]).T,color=color,lw=lw)
def interp(x,z,xx):return np.interp(xx,x,z)
def load_case(L):
 m=json.loads((RUN/'manifest.json').read_text());assert m['center_boundary_motion'] is False
 i=next(i for i,j in enumerate(m['jobs']) if j['case']['length']==L and j['case']['condition']=='fixed_slope' and j['purpose']=='pass')
 p=RUN/f'workers/{i:03d}.npz';j=json.loads(p.with_suffix('.json').read_text());assert sha(p)==j['archive_sha256'];register(p);register(p.with_suffix('.json'))
 case=j['job']['case'];assert case['center_boundary_motion'] is False and case['boundary_method']=='independent_endpoints'
 checkpoint=ROOT/case['checkpoint'];assert sha(checkpoint)==case['checkpoint_sha256'];register(checkpoint)
 with np.load(p) as a:z={k:a[k].copy() for k in a.files if k=='x' or k.startswith('pulse_')}
 with np.load(checkpoint) as cp:
  for name in ['rac','rho','mca_bound','displacement','velocity','tension']:
   np.testing.assert_allclose(z['pulse_field_'+name][0],cp['state_'+name],rtol=1e-12,atol=1e-12)
 return z,j
CASES={L:load_case(L) for L in [23,39]}
T_NORM=Normalize(-14,14);B_NORM=Normalize(0,11);TCMAP=plt.get_cmap('coolwarm');BCMAP=plt.get_cmap('YlOrBr')
records=[];arrays={};audits=[]
def draw_native(ax,L,t):
 z,meta=CASES[L];time=z['pulse_time'];i=int(np.argmin(abs(time-t)));assert abs(time[i]-t)<1e-12
 j=int(np.argmin(abs(time-(t-10))));assert abs(time[j]-(t-10))<1e-12
 x=z['x'];xx=np.linspace(x[0],x[-1],17);u=z['pulse_field_displacement'][i];u0=z['pulse_field_displacement'][j]
 # Eulerian local dashpot slip is the time integral of its relative velocity.
 # It is a connector kinematic diagnostic, not a permanent material/bond identity.
 intv=np.trapezoid(z['pulse_field_velocity'][j:i+1],time[j:i+1],axis=0)
 slip=intv-(u-u0);xm=xx+interp(x,u,xx);xc=xm+interp(x,slip,xx)
 assert np.all(np.diff(xm)>0) and np.all(np.diff(xc)>0)
 b=interp(x,z['pulse_field_mca_bound'][i],xx);T=interp(x,z['pulse_field_tension'][i],(xx[:-1]+xx[1:])/2)
 vel=interp(x,z['pulse_field_velocity'][i],xx)
 # Stored tension must follow the same native displacement derivative.
 dx=L/len(x);left=z['pulse_field_left_boundary_displacement'][i,0];right=z['pulse_field_right_boundary_displacement'][i,0]
 derived=np.gradient(u,dx);derived[0]=(u[1]-left)/(1.5*dx);derived[-1]=(right-u[-2])/(1.5*dx)
 # Interior source gradient is central; endpoint stencil is checked separately in source manifest.
 error=float(np.max(abs(100*derived[1:-1]-z['pulse_field_tension'][i,1:-1])));assert error<1e-9
 for a,bb,tt in zip(xm[:-1],xm[1:],T):spring(ax,a,bb,1.25,TCMAP(T_NORM(tt)))
 for a,bb in zip(xc[:-1],xc[1:]):dashpot(ax,[a,0],[bb,0],PURPLE,lw=1.1,size=.09)
 for a,bb,bi in zip(xm,xc,b):dashpot(ax,[bb,.04],[a,1.21],BCMAP(B_NORM(100*bi)),lw=1.25,size=.07)
 ax.scatter(xm,np.full(17,1.25),s=14,facecolor='white',edgecolor=GREY,zorder=5,lw=.7)
 ax.scatter(xc,np.zeros(17),s=10,marker='s',color=PURPLE,zorder=5)
 ax.scatter(xx,np.full(17,1.6),s=5,color='#bdbdbd',zorder=1)
 # Purple velocity arrows: 0.1 micrometres/second is 1.8 micrometres on the plot.
 ax.quiver(xc,np.full(17,-.27),vel*18,np.zeros(17),angles='xy',scale_units='xy',scale=1,color=PURPLE,width=.003,headwidth=3.5)
 xp=x+u;rac=z['pulse_field_rac'][i];rho=z['pulse_field_rho'][i]
 ax.fill_between(xp,-1.1,-1.1+.58*rac/10,color=RED,alpha=.18);ax.plot(xp,-1.1+.58*rac/10,color=RED,lw=1.4)
 ax.fill_between(xp,-2,-2+.58*rho,color=BLUE,alpha=.16);ax.plot(xp,-2+.58*rho,color=BLUE,lw=1.4)
 for y in [-1.1,-2]:ax.axhline(y,color='#dddddd',lw=.5,zorder=-1)
 ax.set(xlim=(-4,42),ylim=(-2.2,1.85),xticks=[0,10,20,30,40],yticks=[1.6,1.25,0,-.81,-1.71],yticklabels=['Reference','Membrane','Cortex','Rac / 10','Rho / 1'])
 ax.spines[['left','top','right']].set_visible(False);ax.tick_params(axis='y',length=0,labelsize=8);ax.set_xlabel('Horizontal position (µm)')
 p=f'L{L}_t{t}';arrays.update({p+'_reference_x':xx,p+'_membrane_x':xm,p+'_cortex_registered_x':xc,p+'_local_slip':interp(x,slip,xx),p+'_tension':T,p+'_mca_friction':100*b,p+'_cortex_velocity':vel,p+'_native_x':x,p+'_native_displacement':u,p+'_native_rac':rac,p+'_native_rho':rho})
 records.append(dict(length=L,time=t,nodes=17,slip_window_seconds=10,minimum_slip=float(slip.min()),maximum_slip=float(slip.max()),minimum_membrane_displacement=float(u.min()),maximum_membrane_displacement=float(u.max()),minimum_tension=float(T.min()),maximum_tension=float(T.max()),max_interior_tension_identity_error=error,withdrawal_time=meta['withdrawal_time']))

def elliptic(L,N,rho,xi=7.8125):
 dx=L/N;x=(np.arange(N)+.5)*dx;force=100*np.gradient(rho,dx);a=np.zeros((3,N));r=10000/dx**2;a[1]=xi+2*r;a[1,[0,-1]]=xi+3*r;a[0,1:]=-r;a[2,:-1]=-r
 return x,solve_banded((1,1),a,force)
fig=plt.figure(figsize=(12,10.6));g=fig.add_gridspec(3,2,left=.092,right=.865,bottom=.17,top=.955,hspace=.48,wspace=.38,height_ratios=[1.1,1.1,.95])
for row,t in enumerate([30,100]):
 for col,L in enumerate([23,39]):
  ax=fig.add_subplot(g[row,col]);panel(ax,'ABCD'[2*row+col],f'L = {L} µm; {t} s after cue onset');draw_native(ax,L,t)
cax=fig.add_axes([.905,.715,.014,.20]);fig.colorbar(plt.cm.ScalarMappable(norm=T_NORM,cmap=TCMAP),cax=cax,label='Membrane tension k ∂u/∂x (pN)',ticks=[-14,0,14])
cax=fig.add_axes([.905,.425,.014,.20]);fig.colorbar(plt.cm.ScalarMappable(norm=B_NORM,cmap=BCMAP),cax=cax,label='MCA friction μ b (pN s / µm²)',ticks=[0,5,10])
ell=np.sqrt(10000/7.8125);lengths=[23,39,80];colors={23:'#009E73',39:'#0072B2',80:'#333333'}
ax=fig.add_subplot(g[2,0]);panel(ax,'E','Rac-driven endpoint: isolated screening')
for L in lengths:
 d=np.linspace(0,L,501);y=np.sinh((L-d)/ell)/np.sinh(L/ell);ax.plot(d,y,color=colors[L],label=f'L = {L} µm');arrays[f'boundary_L{L}_distance']=d;arrays[f'boundary_L{L}_fraction']=y
ax.plot(np.linspace(0,80,501),np.exp(-np.linspace(0,80,501)/ell),'--',color='#aaaaaa',lw=1,label='Long-domain envelope')
ax.set(xlim=(0,80),ylim=(0,1.04),xlabel='Distance from driven endpoint (µm)',ylabel='Cortical speed / endpoint speed');ax.legend(frameon=False,fontsize=7,loc='upper right')
ax=fig.add_subplot(g[2,1]);panel(ax,'F','Rho gradient: identical local forcing')
z,_=CASES[23];ix=np.argmin(abs(z['pulse_time']-100));refx=z['x'];ref_rho=z['pulse_field_rho'][ix]
for L in lengths:
 vs=[]
 for N in [1001,2001]:
  xx=(np.arange(N)+.5)*L/N;rho=np.interp(xx,refx,ref_rho);x,v=elliptic(L,N,rho);vs.append(float(np.interp(5,x,v)))
  if N==2001:ax.plot(np.r_[0,x,L],np.r_[0,v,0],color=colors[L]);arrays[f'rho_L{L}_x']=x;arrays[f'rho_L{L}_velocity']=v;arrays[f'rho_L{L}_rho']=rho
 error=abs(vs[0]-vs[1])/max(abs(vs[1]),1e-8);assert error<.005;audits.append(dict(length=L,relative_mesh_change_at_x5=error))
ax.axhline(0,color='#aaaaaa',lw=.5);ax.set(xlim=(0,80),xlabel='Distance from original left endpoint (µm)',ylabel='Rho-driven cortical velocity (µm/s)')
handles=[Line2D([],[],color=GREY,label='Membrane springs'),Line2D([],[],color=PURPLE,label='Cortical axial dashpots'),Line2D([],[],color='#ad6418',label='MCA slip dashpots'),Line2D([],[],color=PURPLE,marker='>',ls='',label='Velocity: 1.8 µm arrow = 0.1 µm/s'),Line2D([],[],color=RED,label='Rac concentration / 10'),Line2D([],[],color=BLUE,label='Rho concentration / 1')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.072),ncol=3,frameon=False,fontsize=8,columnspacing=1.5)
outputs={}
for ext in ['png','pdf','svg']:
 p=FIG/('debelly_native_mechanics.'+ext);fig.savefig(p);outputs[p.name]=sha(p)
plt.close(fig)
np.savez_compressed(DATA/'mechanics_fields.npz',**arrays)
with (DATA/'native_snapshot_measurements.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
register(Path(__file__));register(ROOT/'code/core/src/polarity1d/debelly.py');register(ROOT.parents[1]/'literature_models/debelly_2026_mechanochemical/scripts/generate_debelly_dashpot_visuals.py')
caption=('Figure 10. The corrected deBelly mechanical fields as a coupled spring/dashpot system. A–D: the same rightward constant-slope cue (0.24 concentration units per µm) at 30 and 100 s, before withdrawal, in 23- and 39-µm cells. Each chain samples 17 positions. Membrane node positions are native x + u; spring color gives signed native tension, k ∂u/∂x. Purple axial dashpots represent cortical viscosity; arrows give native cortical velocity. Gold MCA dashpots carry local friction μ b, colored on one shared absolute scale. Their inclination shows accumulated local relative slip, integral of v dt minus change in u over the preceding 10 s, with no geometric magnification. The cortex row starts registered to the membrane at the beginning of that window. This fixed-coordinate kinematic rendering does not identify permanent cortex particles or bonds: MCA linkages exchange. Gray dots mark undeformed reference positions. Native Rac and Rho profiles lie below each chain, scaled by 10 and 1 concentration units respectively. E: uniform basal MCA, a stationary membrane, one imposed cortical endpoint speed and a fixed opposite endpoint give v(d)/v_boundary = sinh((L − d)/ell)/sinh(L/ell), where d is distance from the driven endpoint; ell = 35.8 µm. F: the L = 23, t = 100 s Rho profile is fixed in physical reference coordinates and its end value extended for longer control domains; both endpoint velocities and membrane motion are zero. Colors in E–F identify domain length. E–F are isolated quasistatic controls, not snapshots of an additional reversal experiment.')
manifest=dict(passed=True,data_sha256={'mechanics_fields.npz':sha(DATA/'mechanics_fields.npz'),'native_snapshot_measurements.csv':sha(DATA/'native_snapshot_measurements.csv')},generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',caption=caption,figure_sha256=outputs,sources_sha256=sources,reproduction=['OPENBLAS_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_quality_mechanics_20261009.py'],native_equations={'membrane':'μ b (u_t − v) = k u_xx','cortex':'μ b (v − u_t) − η v_xx = ζ Rho_x','tension':'T = k u_x','slip':'s(x,t) = integral from t−10 to t of v(x,q) dq − [u(x,t) − u(x,t−10)]'},parameters=dict(membrane_stiffness=100,cortical_viscosity=10000,linker_friction=100,active_stress=100,basal_bound_mca=.078125,screening_length=ell,center_boundary_motion=False),snapshot_checks=records,mesh_checks=audits,limitations=['Dashpots resist relative velocity, not relative displacement; the displayed connector angle visualizes recent integrated slip.','Cortex points are registered local kinematic markers, not permanently paired material particles.','Prescribed independent endpoint commands and boundary relaxation are retained; this is not free migration.','The 17 drawing nodes resample native 201-cell fields and do not define a new coarse mechanical simulation.'])
(FIG/'debelly_native_mechanics_caption.txt').write_text(caption+'\n')
(FIG/'debelly_native_mechanics.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(DATA/'audit.json').write_text(json.dumps(dict(passed=True,native_snapshots=records,control_mesh=audits,data_sha256=sha(DATA/'mechanics_fields.npz')),indent=2)+'\n')
chapter={'screening':dict(title='Membrane strain and MCA slip make mechanical memory visible',image='quality_20261009_mechanics/debelly_native_mechanics.png',shape='wide',paragraphs=['For node spacing h, the membrane springs have stiffness k/h, cortical axial dashpots have drag eta/h, and MCA linkers have drag mu*b*h. Their continuum balances are mu*b*(u_t − v) = k*u_xx and mu*b*(v − u_t) − eta*v_xx = zeta*Rho_x. Thus Rho supplies a distributed stress gradient, Rac sets endpoint commands, and MCA resists relative velocity; membrane tension is k*u_x.', 'Figure 10A–D shows corrected native fields at matched times under the same gradient cue. The long cell retains a strongly asymmetric strain and flow pattern at 100 s. Springs deform with the measured membrane displacement and carry their measured signed tension. MCA dashpot color shows the local friction coefficient. Their inclination records the preceding 10 s of integrated local slip; it is a kinematic illustration, not elastic bond extension or permanent particle pairing. All four panels share physical position, slip, arrow and color scales.', 'To isolate propagation, Figure 10E freezes the membrane and sets MCA to its uniform basal value. The endpoint response decays over ell = sqrt(eta / (mu b)) = 35.8 µm. Figure 10F holds the local Rho profile fixed and extends its tail: a longer domain changes the response to the same force. These controls explain the mechanical screening operator; they do not assign a single screening length to the evolving, nonuniform cell.'],sources=[['quality_20261009_mechanics/debelly_native_mechanics.pdf','Vector figure'],['quality_20261009_mechanics/debelly_native_mechanics.manifest.json','Legend, data hashes and reproduction'],['../../data/runs/length_mechanisms_20261008/quality_20261009_mechanics/mechanics_fields.npz','Rendered native fields and screening arrays'],['../../data/runs/length_mechanisms_20261008/quality_20261009_mechanics/native_snapshot_measurements.csv','Native snapshot measurements'],['../../data/runs/length_mechanisms_20261008/quality_20261009_mechanics/audit.json','Checkpoint, force identity and mesh audit']],caption=caption)}
(DATA/'chapter_overrides.json').write_text(json.dumps(chapter,indent=2)+'\n')
print(FIG/'debelly_native_mechanics.png');print('Native snapshots and mesh checks passed')
