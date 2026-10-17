"""Scientific profile/kinetics and removal/return mechanism figures; native data only."""
from pathlib import Path
from datetime import datetime, timezone
import csv,hashlib,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import NullFormatter

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code/core'),str(ROOT/'code/core/src')]
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning
OUT=ROOT/'data/runs/length_mechanisms_20261008/comparative_revision'
FIG=ROOT/'figures/length_mechanisms_20261008/comparative_revision'
POOL=ROOT/'data/runs/inactive_pool_figure_exports_20261008'
RUN=ROOT/'data/runs/inactive_pool_restoring_response_20261008'
OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
MODELS={'goryachev':('Goryachev',Goryachev(),'#8c398f'),'wave_pinning':('Wave-pinning',WavePinning(),'#168577')}
LC={9:'#376caf',21:'#bc8127',39:'#a53f47'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
 'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
sources={};audit=[];compact={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def register(p):sources[str(p.relative_to(ROOT))]=sha(p)
def load(p):
 register(p)
 with np.load(p) as z:return {k:z[k].copy() for k in z.files}
def record(p):register(p);return json.loads(p.read_text())
def local_reaction(model,u,v,cue=0.):
 # The native contract accepts two one-dimensional fields; evaluate pointwise.
 u=np.asarray(u);v=np.asarray(v);flat=np.stack([u.ravel(),v.ravel()])
 return model.reaction(0.,np.zeros(u.size),flat,np.full(u.size,cue))[0].reshape(u.shape)
def save(fig,stem,caption):
 fig.subplots_adjust(left=.085,right=.97,bottom=.115,top=.91,hspace=.48,wspace=.28)
 fig.text(.085,.075,caption,fontsize=8.5,va='top',color='#42525d')
 for suffix in ['png','pdf','svg']:fig.savefig(FIG/f'{stem}.{suffix}',dpi=195)
 (FIG/f'{stem}_caption.txt').write_text(caption+'\n');plt.close(fig)

provenance=record(POOL/'provenance.json')
for p,h in provenance['inputs'].items():assert sha(POOL/p)==h
z=load(POOL/'evidence_D0p1.npz')
path=POOL/'response_vs_length.csv';register(path);rows=list(csv.DictReader(path.open()))
lengths=sorted({int(float(r['length'])) for r in rows})
manifest=record(RUN/'manifest.json');status=record(RUN/'results/status.json');assert status['complete'] and status['audit_passed']
assert record(RUN/'results/independent_audit.json')['passed']

fig,axes=plt.subplots(3,2,figsize=(12.3,11.5))
fig.suptitle('Length-dependent concentration profiles and local reaction kinetics',fontsize=16,y=.985)
fig.text(.5,.945,'An added local cue contributes Δf = s v: its initial effect depends on inactive concentration.',ha='center',fontsize=11,color='#42525d')
init_rows=[]
for col,(name,(label,model,color)) in enumerate(MODELS.items()):
 register(ROOT/f'code/core/models/{name}.py')
 assert sha(ROOT/f'code/core/models/{name}.py')==provenance['models'][name]
 for L,c in LC.items():
  u,v=z[f'{name}_L{L}_steady'];x=(np.arange(len(u))+.5)*L/len(u)
  axes[0,col].plot(x,u,color=c,lw=2,label=f'L = {L}')
  axes[1,col].plot(x,v,color=c,lw=2,label=f'L = {L}')
  np.testing.assert_allclose(model.reaction(0.,x,np.stack([u,v]),np.ones_like(x))[0]-model.reaction(0.,x,np.stack([u,v]),np.zeros_like(x))[0],v,rtol=1e-12,atol=1e-12)
 axes[0,col].set(title=f'{"AB"[col]}  {label}: active concentration',xlabel='Position from initial pole',ylabel='Active concentration, u',xlim=(0,40),ylim=(0,None))
 axes[1,col].set(title=f'{"CD"[col]}  {label}: inactive concentration',xlabel='Position from initial pole',ylabel='Inactive concentration, v',xlim=(0,40),ylim=(0,None))
 axes[0,col].legend(frameon=False,ncol=3);axes[1,col].legend(frameon=False,ncol=3)
 umax=max(float(z[f'{name}_L{L}_steady'][0].max()) for L in LC)*1.25
 vmax=max(float(z[f'{name}_L{L}_steady'][1].max()) for L in LC)*1.15
 ug=np.linspace(0,umax,301);vg=np.linspace(0,vmax,261);U,V=np.meshgrid(ug,vg)
 f=local_reaction(model,U,V)
 extent=float(np.max(abs(f)));linthresh=extent*.003
 ax=axes[2,col];heat=ax.pcolormesh(U,V,f,cmap='RdBu_r',norm=SymLogNorm(linthresh=linthresh,vmin=-extent,vmax=extent),shading='auto',rasterized=True)
 if name=='goryachev':null=model.c/(model.ec*(model.a*ug+model.b))
 else:null=model.decay*ug/(model.k0+model.feedback*ug**2/(model.half_saturation**2+ug**2))
 keep=null<=vmax;ax.plot(ug[keep],null[keep],color='#111111',ls='--',lw=1.5,label='Zero net reaction, f = 0')
 np.testing.assert_allclose(model.reaction(0.,ug[keep],np.stack([ug[keep],null[keep]]),np.zeros_like(ug[keep]))[0],0.,atol=2e-14)
 for L,c in LC.items():
  u,v=z[f'{name}_L{L}_steady'];ax.plot(u,v,color=c,lw=1.7)
  ax.plot(u[-1],v[-1],'o',color=c,ms=6,markeredgecolor='white',markeredgewidth=.8)
  ax.plot(u[0],v[0],'^',color=c,ms=5,markeredgecolor='white',markeredgewidth=.6)
 ax.set(title=f'{"EF"[col]}  {label}: basal net reaction, f(u, v)',xlabel='Active concentration, u',ylabel='Inactive concentration, v',xlim=(0,umax),ylim=(0,vmax))
 bar=fig.colorbar(heat,ax=ax,pad=.015,fraction=.045);bar.set_label('Net reaction: production (+), loss (−)',fontsize=8);bar.ax.tick_params(labelsize=7)
 ticks=[-extent,-extent*.01,0.,extent*.01,extent];bar.set_ticks(ticks);bar.set_ticklabels([f'{v:.2g}' for v in ticks])
 bar.ax.set_title('Native rate',fontsize=7,pad=7)
 inset=ax.inset_axes([.58,.54,.37,.30] if name=='goryachev' else [.40,.19,.40,.30])
 vv=np.array([z[f'{name}_L{L}_steady'][1,-1] for L in lengths]);ratios=vv/vv[0]
 inset.plot(lengths,ratios,color=color,lw=1.8);inset.scatter([9,21,39],[ratios[lengths.index(L)] for L in [9,21,39]],c=list(LC.values()),s=18,zorder=3)
 inset.set(xlim=(9,39),ylim=(0,1.06),xticks=[9,21,39],yticks=[0,.5,1],xlabel='Cell length',ylabel='s v / (s v at L = 9)')
 inset.set_title('Equal-cue initiation',fontsize=8);inset.tick_params(labelsize=7);inset.xaxis.label.set_size(7);inset.yaxis.label.set_size(7);inset.grid(alpha=.15)
 for L,vv0,ratio in zip(lengths,vv,ratios):init_rows.append(dict(model=name,length=L,inactive_opposite_end=float(vv0),equal_local_cue_rate_relative_to_L9=float(ratio)))
 compact[f'{name}_reaction_u']=ug;compact[f'{name}_reaction_v']=vg;compact[f'{name}_basal_net_reaction']=f
for ax in axes[:2,:].flat:ax.grid(alpha=.15)
save(fig,'equal_cue_initiation','A–D: exact basal spatial profiles at three cell lengths. E–F: native basal net reaction; added cue is zero and diffusion is excluded.\nSolid curves show spatial profile loci, circles the opposite end, triangles the initial pole; dashed curves show f = 0. Spatial steady states include diffusion.\nInsets evaluate the added source s v at the opposite end for identical local cue s. Concentration and time units are native to each model.')

traces={};restore_rows=[]
for name,(label,model,color) in MODELS.items():
 for L in LC:
  job=next(j for j in manifest['jobs'] if j['case']['model']==name and j['case']['length']==L and j['duration']==.1)
  directory=RUN/'workers'/f'{job["slot"]:03d}'
  a=load(directory/'epsilon_0.005.npz');b=load(directory/'epsilon_0.npz')
  for p in [directory/'epsilon_0.005.npz',directory/'epsilon_0.npz']:
   meta=record(p.with_suffix('.json'));assert sha(p)==meta['archive_sha256']
  checkpoint=ROOT/job['case']['checkpoint'];register(checkpoint);assert sha(checkpoint)==job['case']['checkpoint_sha256']
  x=a['x'];dx=L/len(x);w=np.clip(job['width']-np.arange(len(x))*dx,0.,dx)
  loss=b['recovery_core_active'][0]-a['recovery_core_active'][0]
  assert loss>0
  for phase in ['pulse','recovery']:
   for archive in [a,b]:
    reaction=local_reaction(model,archive[phase+'_active'],archive[phase+'_inactive'])
    np.testing.assert_allclose(reaction@w,archive[phase+'_core_reaction'],rtol=2e-10,atol=2e-12)
  np.testing.assert_allclose(a['pulse_active'][0],a['originalstate'][0],rtol=0,atol=0)
  np.testing.assert_allclose(a['recovery_inactive'][0]-a['pulse_inactive'][-1],a['removed_inactive'],atol=5e-16,rtol=2e-13)
  t=np.r_[a['pulse_time']-.1,a['recovery_time'][1:]]
  fraction=np.r_[(b['pulse_core_active']-a['pulse_core_active'])/b['pulse_core_active'],(b['recovery_core_active'][1:]-a['recovery_core_active'][1:])/b['recovery_core_active'][1:]]
  gain=(a['recovery_core_reaction']-b['recovery_core_reaction'])/loss
  actualgain=(a['recovery_core_total_rate']-b['recovery_core_total_rate'])/loss
  rr=next(r for r in rows if r['model']==name and float(r['length'])==L and float(r['duration'])==.1 and float(r['epsilon'])==.005)
  np.testing.assert_allclose(gain[0],float(rr['reaction_gain_per_damage']),rtol=1e-10)
  np.testing.assert_allclose(loss/b['recovery_core_active'][0],float(rr['fractional_damage']),rtol=1e-10)
  traces[name,L]=(t,fraction,a['recovery_time'],gain,actualgain)
  for key,value in [('time',t),('fractional_core_damage',fraction),('recovery_time',a['recovery_time']),('reaction_per_initial_damage',gain),('actual_repair_per_initial_damage',actualgain)]:compact[f'{name}_L{L}_{key}']=value
  audit.append(dict(model=name,length=L,source_slot=job['slot'],core_width=job['width'],native_reaction_integrals_passed=True,pointwise_return_passed=True,checkpoint_hash_passed=True))

fig,axes=plt.subplots(3,2,figsize=(12.3,11.4))
fig.suptitle('Incumbent-pole damage and restoration after transient inactive-pool reduction',fontsize=15,y=.985)
fig.text(.5,.946,'A 0.5% inactive-profile reduction lasts 0.1 native time units; the exact removed profile is returned pointwise at t = 0.',ha='center',fontsize=10.5,color='#42525d')
for col,(name,(label,model,color)) in enumerate(MODELS.items()):
 for L,c in LC.items():
  t,fraction,rt,gain,actual=traces[name,L]
  axes[0,col].plot(t,100*fraction,color=c,lw=2,label=f'L = {L}')
  axes[1,col].plot(rt,gain,color=c,lw=2,label=f'L = {L}')
  restore_rows.append(dict(model=name,length=L,core_damage_percent_at_return=float(100*fraction[np.argmin(abs(t))]),reaction_per_damage_at_return=float(gain[0]),actual_material_repair_per_damage_at_return=float(actual[0])))
 ax=axes[0,col];ax.axvspan(-.1,0,color='#9ca8b0',alpha=.23);ax.axvline(0,color='#555555',lw=.8)
 ax.set(title=f'{"AB"[col]}  {label}: fractional active-core loss',xlabel='Time relative to inactive-profile return',ylabel='Active core material lost (%)',xscale='symlog',xlim=(-.1,100),ylim=(-.0001,.0011) if name=='goryachev' else (-.005,.055))
 ax.set_xscale('symlog',linthresh=.1);ax.set_xticks([-.1,0,.1,1,10,100],['−0.1','0','0.1','1','10','100']);ax.legend(frameon=False,loc='upper right')
 ax.text(.03,.96,'Inactive held externally',transform=ax.transAxes,fontsize=8,va='top',color='#53616b')
 ax=axes[1,col];ax.axhline(0,color='#888888',lw=.6)
 ax.set(title=f'{"CD"[col]}  {label}: immediate restoring reaction',xlabel='Time after inactive-profile return',ylabel='Excess net reaction / material lost at return\n(1 / native time)',xlim=(0,2),ylim=(-.1,16) if name=='goryachev' else (-.02,.85));ax.legend(frameon=False)
 rr=sorted([r for r in rows if r['model']==name and float(r['duration'])==.1 and float(r['epsilon'])==.005],key=lambda r:float(r['length']))
 axes[2,0].plot([float(r['length']) for r in rr],[100*float(r['fractional_damage']) for r in rr],'o-',color=color,ms=4,lw=2,label=label)
 axes[2,1].plot([float(r['length']) for r in rr],[float(r['reaction_gain_per_damage']) for r in rr],'o-',color=color,ms=4,lw=2,label=label)
axes[2,0].set(title='E  Fractional damage at profile return',xlabel='Cell length',ylabel='Active core material lost (%)',yscale='log');axes[2,0].legend(frameon=False)
axes[2,0].set_yticks([.001,.01,.05],['0.001','0.01','0.05']);axes[2,0].yaxis.set_minor_formatter(NullFormatter())
axes[2,1].set(title='F  Restoring response per unit damage',xlabel='Cell length',ylabel='Immediate excess net reaction / material lost\n(1 / native time)');axes[2,1].legend(frameon=False)
for ax in axes.flat:ax.grid(alpha=.15)
save(fig,'incumbent_restoration','The original half-height incumbent core is fixed. Matched untreated trajectories are subtracted; native kinetics and diffusion are unchanged.\nA–B: shaded removal interval; time is linear near zero and logarithmic beyond 0.1. C–D: post-return net reaction excludes diffusion.\nE–F: 18 lengths, 0.5% removal. Damage fractions compare length dependence within each model; native units are not calibrated between models.')

for stem,data in [('equal_cue_initiation',init_rows),('incumbent_restoration',restore_rows)]:
 with (OUT/f'{stem}.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
np.savez_compressed(OUT/'compact_evidence.npz',**compact)
(OUT/'independent_observable_audit.json').write_text(json.dumps(dict(passed=True,native_reaction_map_calls=True,cue_source_identity_passed=True,nullclines_independently_evaluated=True,cases=audit),indent=2)+'\n')
(OUT/'provenance.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',integration='none',renderer_sha256=sha(Path(__file__)),sources=sources,figures={p.name:sha(p) for p in FIG.iterdir() if p.suffix in ['.png','.pdf','.svg']}),indent=2)+'\n')
print(FIG)
