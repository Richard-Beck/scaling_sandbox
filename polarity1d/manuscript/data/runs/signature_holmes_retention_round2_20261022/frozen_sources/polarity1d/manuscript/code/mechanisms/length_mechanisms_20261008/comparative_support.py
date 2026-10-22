"""Compact native evidence for the length-mechanism reader report; no integration."""
from pathlib import Path
import csv, hashlib, json
from datetime import datetime, timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MANUSCRIPT = Path(__file__).resolve().parents[3]
DATA = MANUSCRIPT/'data/runs/length_mechanisms_20261008/comparative_support'
FIG = MANUSCRIPT/'figures/length_mechanisms_20261008/support'
POOL = MANUSCRIPT/'data/runs/inactive_pool_figure_exports_20261008'
NATIVE = MANUSCRIPT/'data/runs/three_test_mechanisms_20261006/trajectories'
DATA.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
COLORS = {'goryachev':'#8c398f', 'wave_pinning':'#168577'}
LABELS = {'goryachev':'Goryachev', 'wave_pinning':'Wave-pinning'}
LC = {9:'#376caf',21:'#bc8127',39:'#a53f47'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
 'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,
 'legend.fontsize':8,'pdf.fonttype':42})
sources={}
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read_npz(path):
 sources[str(path.relative_to(MANUSCRIPT))]=sha(path)
 with np.load(path) as z: return {k:z[k].copy() for k in z.files}
def save(fig,stem,caption):
 fig.tight_layout(rect=(0,.11,1,.94))
 fig.text(.055,.075,caption,va='top',fontsize=8.5,color='#42525d')
 for suffix in ['png','pdf','svg']:fig.savefig(FIG/f'{stem}.{suffix}',dpi=190)
 plt.close(fig);(FIG/f'{stem}_caption.txt').write_text(caption+'\n')
def width(u,x,side='left',threshold=.6):
 if side=='right': return width(u[::-1],x[-1]+x[0]-x[::-1], 'left', threshold)
 if u[0]<=threshold:return 0.
 stop=np.flatnonzero(u<=threshold)
 if not len(stop):return x[-1]+x[0]
 j=stop[0]
 return float(x[j-1]+(threshold-u[j-1])/(u[j]-u[j-1])*(x[j]-x[j-1]))

z=read_npz(POOL/'evidence_D0p1.npz')
pool_provenance=json.loads((POOL/'provenance.json').read_text())
for filename,digest in pool_provenance['inputs'].items():assert sha(POOL/filename)==digest
path=POOL/'response_vs_length.csv';sources[str(path.relative_to(MANUSCRIPT))]=sha(path)
rows=list(csv.DictReader(path.open()))
lengths=sorted({int(float(r['length'])) for r in rows})
quant=[]

fig,axes=plt.subplots(1,3,figsize=(12,4.4))
fig.suptitle('An equal local cue activates less material when inactive concentration is lower',fontsize=14)
for col,m in enumerate(COLORS):
 for L in [9,21,39]:
  v=z[f'{m}_L{L}_steady'][1];x=(np.arange(len(v))+.5)*L/len(v)
  axes[col].plot(x,v,color=LC[L],lw=2,label=f'L = {L}')
 axes[col].set(title=f'{"AB"[col]}  {LABELS[m]} inactive pool',xlabel='Position from incumbent end',ylabel='Inactive concentration',xlim=(0,40),ylim=(0,None))
 axes[col].legend(frameon=False)
 values=[float(z[f'{m}_L{L}_steady'][1,-1]) for L in lengths]
 ratios=np.array(values)/values[0]
 axes[2].plot(lengths,ratios,'o-',color=COLORS[m],lw=2,ms=4,label=LABELS[m])
 for L,v,ratio in zip(lengths,values,ratios):quant.append(dict(model=m,length=L,inactive_opposite_end=v,equal_local_cue_rate_relative_to_L9=float(ratio)))
axes[2].set(title='C  Initial activation by the same local cue',xlabel='Cell length',ylabel='Added local activation rate / rate at L = 9',ylim=(0,1.07));axes[2].legend(frameon=False)
for ax in axes:ax.grid(alpha=.15)
save(fig,'equal_cue_initiation','Native pre-cue checkpoints. Added local activation = cue strength × inactive concentration.\nEach model is normalized to its own L = 9 rate in C; this compares length effects under equal local input.')

fig,axes=plt.subplots(2,2,figsize=(11,7.5))
fig.suptitle('Taller Goryachev peaks resist a brief loss of inactive substrate',fontsize=14)
for col,m in enumerate(COLORS):
 p=f'{m}_L21_';x=z[p+'x'];ax=axes[0,col]
 active=z[f'{m}_L21_steady'][0];shade=ax.twinx();shade.fill_between(x,0,active/active.max(),color='#777777',alpha=.11)
 shade.set(ylim=(0,1.1),yticks=[]);shade.spines['right'].set_visible(False)
 ax.plot(x,z[p+'deficit_release'],color='#3869a3',lw=2,label='At inactive return')
 ax.plot(x,z[p+'deficit_later'],color='#c57d26',lw=2,label=f'{float(z[p+"later_time"]):.2f} after return')
 ax.axvline(float(z[p+'core_width']),color='#888888',ls=':',lw=1)
 ax.set(title=f'{"AB"[col]}  {LABELS[m]}: location of residual damage',xlabel='Position from incumbent end',ylabel='Active deficit / initial core material lost',xlim=(0,14))
 ax.axhline(0,color='#aaa',lw=.6);ax.legend(frameon=False)
 rr=sorted([r for r in rows if r['model']==m and float(r['duration'])==.1 and float(r['epsilon'])==.005],key=lambda r:float(r['length']))
 axes[1,0].plot([float(r['length']) for r in rr],[100*float(r['fractional_damage']) for r in rr],'o-',color=COLORS[m],label=LABELS[m],ms=4,lw=2)
 axes[1,1].plot([float(r['length']) for r in rr],[float(r['reaction_gain_per_damage']) for r in rr],'o-',color=COLORS[m],label=LABELS[m],ms=4,lw=2)
axes[1,0].set(title='C  Incumbent damage',xlabel='Cell length',ylabel='Active core material lost (%)',yscale='log');axes[1,0].legend(frameon=False)
axes[1,1].set(title='D  Immediate restoring response',xlabel='Cell length',ylabel='Net reaction / active core material lost\n(1 / native time)');axes[1,1].legend(frameon=False)
for ax in axes.flat:ax.grid(alpha=.15)
save(fig,'incumbent_restoration','Remove 0.5% of the entire inactive profile for 0.1 native time units; return that exact material pointwise.\nNative reactions and diffusion remain unchanged. Matched sham subtraction; fixed initial half-height core.\nA–B: L = 21, gray = original active shape, dotted = core edge. D measures reaction, excluding diffusion.')

tr={s:read_npz(NATIVE/f'trajectory_{s:04d}.npz') for s in [18,21,22,23]}
for s in tr:
 p=NATIVE/f'trajectory_{s:04d}.json';sources[str(p.relative_to(MANUSCRIPT))]=sha(p)
 metadata=json.loads(p.read_text());assert sha(p.with_suffix('.npz'))==metadata['archive_sha256']
 checkpoint=MANUSCRIPT/metadata['job']['case']['checkpoint']
 assert sha(checkpoint)==metadata['checkpoint_sha256']
q=tr[23];x=q['pulse_x'];D=float(q['pulse_time'][-1]);L=x[-1]+x[0]
fig,axes=plt.subplots(2,2,figsize=(11.8,8.2));fig.suptitle('Wave-pinning: lasting redistribution comes from moving the high–low interfaces',fontsize=14)
for target,c in zip([0,25,200,D],['#aaa','#bb9b55','#498ca5','#8f4148']):
 i=np.argmin(abs(q['pulse_time']-target));u=q['pulse_active'][i]
 axes[0,0].plot(x,u,color=c,lw=2,label=f'Cue time {q["pulse_time"][i]:g}')
axes[0,0].axvspan(L*.8,L,color='#b1c4cb',alpha=.15)
axes[0,0].set(title='A  Old domain narrows as the new domain expands',xlabel='Position',ylabel='Active concentration',xlim=(0,L),ylim=(0,1.68));axes[0,0].legend(frameon=False)
for target,c in zip([0,25,1000,10000],['#8f4148','#bc9946','#478984','#283e61']):
 i=np.argmin(abs(q['recovery_time']-target));axes[0,1].plot(x,q['recovery_active'][i],color=c,lw=2,label=f'After withdrawal {q["recovery_time"][i]:,.0f}')
axes[0,1].set(title='B  Cue-induced height relaxes; retained growth is wider',xlabel='Position',ylabel='Active concentration',xlim=(0,L),ylim=(0,1.68));axes[0,1].legend(frameon=False)
for side,color,label in [('left','#786182','Old high-active domain'),('right','#168577','New high-active domain')]:
 ww=np.array([width(u,x,side) for u in q['pulse_active']]);axes[1,0].plot(q['pulse_time'],ww,color=color,lw=2,label=label)
axes[1,0].set(title='C  Both poles change mainly by edge displacement',xlabel='Time while cue is on',ylabel='Contiguous domain width (active > 0.6)',xlim=(0,D));axes[1,0].legend(frameon=False)
for s,c,ls,label in [(21,'#3869a3','-','Pulse 204: just passes orientation'),(22,'#3869a3','--','Pulse 202: fails orientation'),(23,'#8f4148','-','Pulse 615: recovers initial strength')]:
 a=tr[s];t=a['recovery_time'];axes[1,1].plot(t,a['recovery_orientation'],color=c,ls=ls,lw=2,label=label)
axes[1,1].axhline(.05,color='#888',ls=':',lw=1)
axes[1,1].set(title='D  Lasting polarity requires sufficient redistribution',xlabel='Time after withdrawal',ylabel='Active polarity (first moment)',xscale='symlog',xlim=(0,10000),ylim=(0,.46));axes[1,1].set_xticks([0,10,100,1000,10000],['0','10','100','1,000','10,000']);axes[1,1].legend(frameon=False,loc='upper left')
for ax in axes.flat:ax.grid(alpha=.15)
save(fig,'wave_pinning_interface','Native L = 39, relative-20% patch, exact shared basal checkpoint. A–C use the 615-unit strength pulse.\nShading in A marks cue support. Domain width uses the same active = 0.6 threshold at every time.\nD retains adjacent orientation pass/fail doses; dotted line = orientation gate. Native time units.')

fig,axes=plt.subplots(1,2,figsize=(11,4.8));fig.suptitle('Interface growth increases less than the length-dependent material burden',fontsize=14)
length_measures=[]
for s,c in [(18,'#bc8127'),(23,'#a53f47')]:
 a=tr[s];xx=a['pulse_x'];LL=float(xx[-1]+xx[0]);tt=a['pulse_time'];u=a['pulse_active'];initial_width=width(u[0],xx,'left')
 for target,ls in [(25,'--'),(60,'-')]:
  k=int(np.argmin(abs(tt-target)));axes[0].plot(LL-xx,u[k],ls=ls,color=c,lw=2,label=f'L = {LL:g}, cue time {tt[k]:g}')
 keep=(tt>=25)&(tt<=60);ww=np.array([width(uu,xx,'right') for uu in u[keep]])
 axes[1].plot(tt[keep],ww/initial_width,color=c,lw=2,label=f'L = {LL:g}')
 i25=int(np.argmin(abs(tt-25)));i60=int(np.argmin(abs(tt-60)));w25=width(u[i25],xx,'right');w60=width(u[i60],xx,'right');rate=(w60-w25)/(tt[i60]-tt[i25])
 length_measures.append(dict(length=LL,initial_incumbent_width=initial_width,new_width_at_25=w25,new_width_at_60=w60,average_interface_advance_25_to_60=float(rate),normalized_average_interface_advance=float(rate/initial_width)))
axes[0].set(title='A  New-domain growth during the same cue interval',xlabel='Distance from the new pole',ylabel='Active concentration',xlim=(0,10),ylim=(0,1.5));axes[0].legend(frameon=False)
axes[1].set(title='B  Width gained relative to the starting incumbent',xlabel='Time while cue is on',ylabel='New domain width / initial incumbent width',xlim=(25,60),ylim=(.32,.56));axes[1].legend(frameon=False)
for ax in axes:ax.grid(alpha=.15)
rate_ratio=length_measures[1]['average_interface_advance_25_to_60']/length_measures[0]['average_interface_advance_25_to_60']
burden_ratio=length_measures[1]['initial_incumbent_width']/length_measures[0]['initial_incumbent_width']
save(fig,'wave_pinning_length_growth',f'Native relative-20% patch, equal local amplitude. Domain width = active > 0.6. Common cue times 25–60.\nInterface advance is {rate_ratio:.2f}× faster at L = 39, while the starting incumbent is {burden_ratio:.2f}× wider: the fractional gain is slower.')
(DATA/'wave_pinning_length_growth.json').write_text(json.dumps(length_measures,indent=2)+'\n')

with (DATA/'equal_cue_initiation.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(quant[0]));w.writeheader();w.writerows(quant)
measures=[]
for s,a in tr.items():
 for target in [0.,25.,1000.,10000.]:
  i=int(np.argmin(abs(a['recovery_time']-target)));u=a['recovery_active'][i];t=float(a['recovery_time'][i]);xx=a['recovery_x'];LL=xx[-1]+xx[0];dx=LL/len(xx)
  measures.append(dict(slot=s,length=float(LL),pulse_duration=float(a['pulse_time'][-1]),recovery_time=t,new_pole_peak=float(u[len(xx)//2:].max()),new_domain_width=width(u,xx,'right'),old_domain_width=width(u,xx,'left'),new_half_active_amount=float(u[len(xx)//2:].sum()*dx),orientation=float(a['recovery_orientation'][i])))
with (DATA/'wave_pinning_interface_measurements.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(measures[0]));w.writeheader();w.writerows(measures)
(DATA/'provenance.json').write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',integration='none; native audited histories only',sources=sources,renderer_sha256=sha(Path(__file__)),figures={p.name:sha(p) for p in FIG.iterdir() if p.suffix in ['.png','.pdf','.svg']},measurements=measures),indent=2)+'\n')
print(FIG)
