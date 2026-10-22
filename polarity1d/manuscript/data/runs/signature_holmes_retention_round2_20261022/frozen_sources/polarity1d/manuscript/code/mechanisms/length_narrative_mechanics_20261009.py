"""Resting native and equal-command spring/dashpot comparison; no integration."""
from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.linalg import solve_banded

ROOT = Path(__file__).resolve().parents[2]
FIGROOT = ROOT/'figures/length_mechanisms_20261008'
FIG = FIGROOT/'narrative_20261009/mechanics'
DATA = ROOT/'data/runs/length_mechanisms_20261008/narrative_20261009/mechanics'
BASE = ROOT/'data/runs/debelly_corrected_three_test_20261008'
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(ROOT/'code'), str(ROOT/'code/core')]
import debelly_corrected_reversal as corrected
ds = corrected.n.ds
p = corrected.n.model_for('debelly').parameters
assert not p.center_boundary_motion
BLUE, RED, PURPLE, GRAY = '#0072B2', '#D55E00', '#762a83', '#666666'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
 'axes.titlesize':10, 'axes.labelsize':9, 'axes.spines.top':False,
 'axes.spines.right':False, 'pdf.fonttype':42, 'svg.fonttype':'none'})
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
sources = {}
def register(path):
 sources[str(Path(path).relative_to(ROOT))] = sha(path)

def independent_solve(L, b, rho, displacement_left, command_left):
 """Direct stationary finite-volume solves, separate assembly from native helper."""
 N = len(b); dx = L/N
 def solve(diag, coefficient, rhs, left, right):
  r = coefficient/dx**2
  ab = np.zeros((3,N)); ab[1] = diag+2*r
  ab[1,0] += r; ab[1,-1] += r
  ab[0,1:] = -r; ab[2,:-1] = -r
  rr = rhs.copy(); rr[0] += 2*r*left; rr[-1] += 2*r*right
  result = solve_banded((1,1), ab, rr)
  residual = (diag+2*r)*result
  residual[0] += r*result[0]; residual[-1] += r*result[-1]
  residual[:-1] -= r*result[1:]; residual[1:] -= r*result[:-1]
  err = float(np.max(abs(residual-rr))/max(1,float(np.max(abs(rr)))))
  assert err < 1e-10
  return result, err
 friction = p.linker_friction*np.maximum(b,1e-10)
 grad_rho = np.gradient(rho,dx)
 v, ev = solve(friction,p.cortical_viscosity,p.active_stress*grad_rho,command_left,0.)
 u, eu = solve(np.zeros(N),p.membrane_stiffness,friction*v,displacement_left,0.)
 tension = p.membrane_stiffness*ds._gradient(u,dx,displacement_left,0.)
 # Independent agreement with exact native operators at the stationary limit.
 np.testing.assert_allclose(v,ds._dirichlet_elliptic(friction,p.cortical_viscosity,p.active_stress*ds._gradient(rho,dx),dx,command_left,0.),rtol=2e-11,atol=1e-12)
 np.testing.assert_allclose(u,ds._dirichlet_elliptic(np.zeros(N),p.membrane_stiffness,friction*v,dx,displacement_left,0.),rtol=2e-11,atol=1e-12)
 return u,v,tension,dict(velocity_residual=ev,displacement_residual=eu)

manifest = json.loads((BASE/'manifest.json').read_text()); register(BASE/'manifest.json')
register(BASE/'results/independent_audit.json')
assert json.loads((BASE/'results/independent_audit.json').read_text())['passed']
LENGTHS = sorted({int(c['length']) for c in manifest['cases'] if c['eligible'] and c['model']=='debelly'})
CASES = {}; rows = []; arrays = {}; checks=[]
COMMAND = -.04
DISPLACEMENT = .5*COMMAND/(1-p.boundary_relaxation)
assert abs(DISPLACEMENT+1)<1e-12
for L in LENGTHS:
 cs = [c for c in manifest['cases'] if c['length']==L and c['eligible'] and c['model']=='debelly']
 assert len({c['checkpoint_sha256'] for c in cs})==1
 c = cs[0]; cp=ROOT/c['checkpoint']; register(cp)
 assert sha(cp)==c['checkpoint_sha256'] and not c['center_boundary_motion']
 with np.load(cp) as z:
  native={name:z['state_'+name].copy() for name in ['rac','rho','mca_bound','displacement','velocity','tension','left_boundary_displacement','right_boundary_displacement']}
 x = (np.arange(len(native['rac']))+.5)*L/len(native['rac'])
 rac0,_=ds._homogeneous_rac_rho(p,.5)
 q = p.protrusion_speed*float(ds.smooth_threshold(native['rac'][0],p.protrusion_switch_scale,1.3*rac0))/(1+native['tension'][0]**2/p.protrusion_tension_scale)
 u,v,t,audit=independent_solve(L,native['mca_bound'],native['rho'],DISPLACEMENT,COMMAND)
 eq = dict(native,displacement=u,velocity=v,tension=t,left_boundary_displacement=np.array([DISPLACEMENT]),right_boundary_displacement=np.array([0.]))
 # Frozen chemistry predicts native resting mechanics with its own native command.
 un,vn,tn,_ = independent_solve(L,native['mca_bound'],native['rho'],float(native['left_boundary_displacement'][0]),q)
 native_error = max(float(np.max(abs(un-native['displacement']))),float(np.max(abs(vn-native['velocity']))))
 assert native_error<2e-5,(L,native_error)
 # Uniform attachment/flat Rho is the concise geometry-only reference.
 uc,vc,tc,_=independent_solve(L,np.full_like(x,.078125),np.ones_like(x),DISPLACEMENT,COMMAND)
 CASES[L] = dict(x=x,native=native,equal=eq,command=q)
 rows.append(dict(length_um=L,native_command_um_per_s=q,native_endpoint_displacement_um=float(native['left_boundary_displacement'][0]),native_front_tension_pN=float(native['tension'][0]),equal_command_um_per_s=COMMAND,equal_endpoint_displacement_um=DISPLACEMENT,equal_front_tension_pN=float(t[0]),uniform_control_front_tension_pN=float(tc[0]),native_static_reconstruction_max_error=native_error))
 checks.append(dict(length=L,checkpoint_sha256=sha(cp),native_static_reconstruction_max_error=native_error,**audit))
 for kind,fields in [('native',native),('equal',eq)]:
  for name,value in fields.items(): arrays[f'L{L}_{kind}_{name}']=value
 arrays[f'L{L}_x']=x

def spring(ax,a,b,y,color):
 xx=np.r_[a,np.linspace(a+.12*(b-a),b-.12*(b-a),13),b]
 yy=np.full_like(xx,y); yy[2:-2]+=.070*np.resize([1,-1],len(xx)-4)
 ax.plot(xx,yy,color=color,lw=1.6,solid_capstyle='round')
def dashpot(ax,p0,p1,color,lw=1.15,size=.07):
 p0,p1=np.asarray(p0),np.asarray(p1); d=p1-p0
 n=np.array([-d[1],d[0]])/np.linalg.norm(d)
 b=p0+.32*d; c=p0+.70*d; piston=p0+.50*d
 ax.plot(*np.array([p0,b]).T,color=color,lw=lw)
 ax.plot(*np.array([b+size*n,c+size*n,c-size*n,b-size*n]).T,color=color,lw=lw)
 ax.plot(*np.array([piston-size*n,piston+size*n]).T,color=color,lw=lw)
 ax.plot(*np.array([piston,p1]).T,color=color,lw=lw)
TNORM=TwoSlopeNorm(vmin=-3.5,vcenter=0,vmax=9); BNORM=Normalize(0,11)
TCMAP=plt.get_cmap('coolwarm'); BCMAP=plt.get_cmap('YlOrBr')
def panel(ax,letter,title):
 ax.text(-.10,1.15 if letter in 'EF' else 1.04,letter,transform=ax.transAxes,fontsize=12,fontweight='bold')
 ax.set_title(title,loc='left',pad=10)
def chain(ax,L,kind):
 case=CASES[L]; z=case[kind]; x=case['x']
 nodes=np.linspace(0,L,17)
 u=np.interp(nodes,np.r_[0,x,L],np.r_[z['left_boundary_displacement'][0],z['displacement'],z['right_boundary_displacement'][0]])
 pos=nodes+u; assert np.all(np.diff(pos)>0)
 b=np.interp(nodes,x,z['mca_bound']); tension=np.interp(.5*(nodes[:-1]+nodes[1:]),x,z['tension'])
 vv=np.interp(nodes,np.r_[0,x,L],np.r_[case['command'] if kind=='native' else COMMAND,z['velocity'],0.])
 ax.scatter(nodes,np.full(17,1.60),s=7,color='#bbbbbb',zorder=1)
 for a,bb,tt in zip(pos[:-1],pos[1:],tension):spring(ax,a,bb,1.20,TCMAP(TNORM(tt)))
 for a,bb in zip(pos[:-1],pos[1:]):dashpot(ax,[a,0],[bb,0],PURPLE,size=.085)
 for a,bi in zip(pos,b):dashpot(ax,[a,.05],[a,1.15],BCMAP(BNORM(p.linker_friction*bi)),lw=1.35,size=.07)
 ax.scatter(pos,np.full(17,1.20),s=16,facecolor='white',edgecolor=GRAY,lw=.7,zorder=5)
 ax.scatter(pos,np.zeros(17),s=9,color=PURPLE,marker='s',zorder=5)
 ax.quiver(pos,np.full(17,-.25),vv*30,np.zeros(17),angles='xy',scale_units='xy',scale=1,color=PURPLE,width=.004,headwidth=3.2)
 xp=x+z['displacement']
 for name,base,color in [('rac',-1.12,RED),('rho',-1.92,BLUE)]:
  ax.axhline(base,color='#dddddd',lw=.6)
  ax.fill_between(xp,base,base+.60*z[name],color=color,alpha=.15)
  ax.plot(xp,base+.60*z[name],color=color,lw=1.4)
  # One native concentration unit displayed with the same height in both rows.
  ax.plot([41,41],[base,base+.60],color=color,lw=1)
  ax.text(41.4,base+.30,'1',color=color,fontsize=7,va='center')
 q=case['command'] if kind=='native' else COMMAND
 ax.text(.015,.995,f'Left endpoint: u = {z["left_boundary_displacement"][0]:.2f} µm; command = {q:.3f} µm/s',transform=ax.transAxes,fontsize=8,color=GRAY,va='top')
 ax.set(xlim=(-3.5,43),ylim=(-2.04,1.93),xticks=[0,10,20,30,40],yticks=[1.60,1.20,0,-.82,-1.62],yticklabels=['Reference','Membrane','Cortex','Rac','Rho'])
 ax.spines[['left','top','right']].set_visible(False)
 ax.tick_params(axis='y',length=0,labelsize=8)
 ax.set_xlabel('Physical position (µm)')

fig=plt.figure(figsize=(12,10.2))
g=fig.add_gridspec(3,2,left=.092,right=.867,bottom=.137,top=.95,wspace=.33,hspace=.47,height_ratios=[1.05,1.05,.80])
for row,kind in enumerate(['native','equal']):
 for col,L in enumerate([23,39]):
  ax=fig.add_subplot(g[row,col]);panel(ax,'ABCD'[row*2+col],f'{"Native resting state" if kind=="native" else "Equal-command mechanics"} · L = {L} µm');chain(ax,L,kind)
cax=fig.add_axes([.91,.704,.014,.198]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=TNORM,cmap=TCMAP),cax=cax,label='Membrane tension (pN)',ticks=[-3,0,3,6,9]);cb.solids.set_rasterized(False)
cax=fig.add_axes([.91,.407,.014,.198]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=BNORM,cmap=BCMAP),cax=cax,label='MCA drag density (pN s / µm²)',ticks=[0,5,10]);cb.solids.set_rasterized(False)
ax=fig.add_subplot(g[2,0]);panel(ax,'E','Equal commands change the tension trend')
ls=np.array(LENGTHS)
ax.plot(ls,[r['native_front_tension_pN'] for r in rows],color=GRAY,marker='o',ms=3.8,label='Native resting state')
ax.plot(ls,[r['equal_front_tension_pN'] for r in rows],color=BLUE,marker='o',ms=3.8,label='Equal command; native Rho + MCA frozen')
ax.plot(ls,[r['uniform_control_front_tension_pN'] for r in rows],color='#009E73',ls='--',lw=1.5,label='Equal command; uniform MCA, flat Rho')
ax.axhline(0,color='#aaaaaa',lw=.6)
ax.set(xlim=(17,40),xlabel='Reference length (µm)',ylabel='Front (left-end) tension (pN)',xticks=[18,23,29,35,39])
ax.legend(frameon=False,fontsize=7.3,loc='upper right')
ax=fig.add_subplot(g[2,1]);panel(ax,'F','Geometry alone spreads endpoint-driven flow')
ell=np.sqrt(p.cortical_viscosity/(p.linker_friction*.078125))
for L,color in [(23,'#009E73'),(39,BLUE)]:
 d=np.linspace(0,L,501); y=np.sinh((L-d)/ell)/np.sinh(L/ell)
 ax.plot(d,y,color=color,lw=1.7,label=f'L = {L} µm')
 arrays[f'screening_L{L}_distance']=d;arrays[f'screening_L{L}_relative_flow']=y
 # Independent finite-volume verification of the analytical shape.
 xx=(np.arange(2001)+.5)*L/2001
 _,vv,_,_=independent_solve(L,np.full_like(xx,.078125),np.ones_like(xx),DISPLACEMENT,COMMAND)
 err=float(np.max(abs(vv/COMMAND-np.sinh((L-xx)/ell)/np.sinh(L/ell))))
 assert err<1e-6
 checks.append(dict(control_length=L,uniform_velocity_analytical_max_error=err))
ax.set(xlim=(0,40),ylim=(0,1.03),xlabel='Distance from front (µm)',ylabel='Cortical speed / command')
ax.legend(frameon=False,fontsize=8,loc='upper right')
ax.text(.97,.51,f'Uniform MCA; flat Rho\nScreening length {ell:.1f} µm',transform=ax.transAxes,fontsize=8,ha='right',color=GRAY)
handles=[Line2D([],[],color=GRAY,label='Membrane springs'),Line2D([],[],color=PURPLE,label='Cortical dashpots'),Line2D([],[],color='#ad6418',label='MCA dashpots'),Line2D([],[],color=PURPLE,marker='<',ls='',label='Flow: 3 µm arrow = 0.1 µm/s'),Line2D([],[],color=RED,label='Rac'),Line2D([],[],color=BLUE,label='Rho')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.50,.034),ncol=3,frameon=False,fontsize=8,columnspacing=2)
outputs={}
for ext in ['png','pdf','svg']:
 path=FIG/('debelly_resting_equal_commands.'+ext);fig.savefig(path,dpi=240,facecolor='white');outputs[path.name]=sha(path)
plt.close(fig)
np.savez_compressed(DATA/'mechanics_fields.npz',**arrays)
with (DATA/'tension_command_measurements.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
caption=('A–B: exact corrected pre-cue resting states at L = 23 and 39 µm. C–D: stationary mechanical responses with each length’s native Rho and bound MCA profiles frozen, but identical left cortical command −0.040 µm/s and left endpoint displacement −1.00 µm; both right values are zero. The native endpoint update is u_next = 0.98 u + 0.5 command, so these equal boundary values obey its stationary relation u = 25 command. No chemistry is evolved in C–D. Each chain samples 17 nodes from the 201-cell native grid. Membrane springs are placed at measured x + u and colored by signed tension T = k du/dx; gray dots show undeformed reference nodes. MCA dashpot color gives local drag density mu b. Purple axial dashpots depict cortical viscosity and arrows show cortical velocity on one shared physical scale. The cortex drawing is registered locally beneath the membrane; vertical MCA connectors do not depict elastic extension or permanent molecular pairing. Rac and Rho retain their native concentration profiles and share the same vertical concentration scale (unit bars at right). E: first-cell front tension at all 12 eligible lengths; the blue control retains each native Rho/MCA distribution, whereas the green reference has uniform basal MCA b = 0.078125 and flat Rho. F: the same uniform reference gives normalized flow sinh((L − d)/ell)/sinh(L/ell), ell = sqrt(eta/(mu b)) = 35.8 µm. Stationary membrane balance is −k u_xx = mu b v; cortical balance is mu b v − eta v_xx = zeta Rho_x. Spring stiffness, cortical drag and MCA drag for reference node spacing h are k/h, eta/h and mu b h. All mechanical controls preserve independent prescribed endpoints. Front tension refers to the old left pole, not the more highly tensed right end of long native resting cells.')
(FIG/'debelly_resting_equal_commands_caption.txt').write_text(caption+'\n')
for path in [Path(__file__),ROOT/'code/core/src/polarity1d/debelly.py',BASE/'provenance.json']:register(path)
datahash={name:sha(DATA/name) for name in ['mechanics_fields.npz','tension_command_measurements.csv']}
audit=dict(passed=True,checks=checks,no_new_integrations=True,data_sha256=datahash)
(DATA/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
record=dict(passed=True,generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',caption=caption,legend=caption,sha256=outputs,figure_sha256=outputs,source_sha256=sources,renderer_sha256=sha(Path(__file__)),data_sha256=datahash,reproduction=['OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/length_narrative_mechanics_20261009.py'],control=dict(left_command_um_per_s=COMMAND,left_endpoint_displacement_um=DISPLACEMENT,right_command_um_per_s=0,right_endpoint_displacement_um=0,native_Rho_and_MCA_frozen=True,stationary_membrane=True),checks=checks)
(FIG/'debelly_resting_equal_commands.manifest.json').write_text(json.dumps(record,indent=2)+'\n')
chapter={'screening':dict(title='Cell length changes the mechanical response to an equal command',image='narrative_20261009/mechanics/debelly_resting_equal_commands.png',shape='wide',renderer='length_narrative_mechanics_20261009.py',paragraphs=['The resting spring/dashpot chains make the starting mechanical difference concrete. The long cell combines greater left-end displacement and cortical flow with a strongly uneven membrane tension. Rac sets the endpoint command, Rho supplies distributed cortical stress, and MCA transmits drag between the moving cortex and membrane. These measured fields therefore describe a coupled mechanical state that a reversal cue must reshape.','Equalizing the endpoint commands separates that native starting-state difference from the mechanical response. With Rho and MCA frozen at their native profiles, imposing the same stationary endpoint displacement and cortical command produces different tension distributions at different lengths, and changes the front-tension trend. A uniform-MCA, flat-Rho reference shows how geometry alone spreads endpoint-driven flow. Together these comparisons explain why similar local activity need not produce similar whole-cell mechanics.'],caption=caption,sources=[['narrative_20261009/mechanics/debelly_resting_equal_commands.pdf','Vector figure'],['narrative_20261009/mechanics/debelly_resting_equal_commands.manifest.json','Legend, hashes and reproduction'],['../../data/runs/length_mechanisms_20261008/narrative_20261009/mechanics/tension_command_measurements.csv','Native/equal-command tension measurements'],['../../data/runs/length_mechanisms_20261008/narrative_20261009/mechanics/mechanics_fields.npz','Full native and control fields'],['../../data/runs/length_mechanisms_20261008/narrative_20261009/mechanics/audit.json','Checkpoint, stationary reconstruction and analytical checks']])}
(DATA/'chapter_overrides.json').write_text(json.dumps(chapter,indent=2)+'\n')
print(FIG/'debelly_resting_equal_commands.png')
print('All checkpoint, stationary reconstruction, native-operator and analytical checks passed.')
