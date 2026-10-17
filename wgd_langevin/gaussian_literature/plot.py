from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
def read(n):return list(csv.DictReader((OUT/n).open()))
fits=read('fits.csv');surfaces=read('parameter_surfaces.csv');pred=read('predictions.csv');times=read('time_predictions.csv');slow=read('slow_fraction_sensitivity.csv');comp=read('model_comparison.csv');tracks=np.load(OUT/'tracks.npz');checks=json.load((OUT/'checks.json').open())
F={(f['study'],f['group']):f for f in fits};COL={'reference':'#286FA7','polyploid':'#C5522E'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':12,'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(24,25),facecolor='white')
fig.text(.05,.978,'WGD migration with Gaussian velocity noise: two parameters per population',fontsize=23,weight='bold',va='top')
fig.text(.05,.958,'Velocity memory tau and Wiener noise coefficient η; RMS speed = η √tau. Exact stationary Gaussian OU sampling; homogeneous, unbiased populations.',fontsize=11.5,color='#555555')
fig.text(.05,.942,'Original fitted summaries and reading intervals retained. Reference = parental/diploid comparator; polyploid = PGCC/PACC or selected tetraploid clone, depending on the study.',fontsize=10.5,color='#555555')
fig.legend(handles=[Line2D([],[],color=COL[g],lw=3,label=l) for g,l in [('reference','Reference'),('polyploid','Polyploid')]],ncol=2,bbox_to_anchor=(.96,.935),loc='upper right',frameon=False)

# A: two-parameter compatibility regions in physically readable coordinates.
for j,study in enumerate(('Xuan2018','Mallin2023')):
 ax=fig.add_axes([.075+j*.31,.762,.24,.14]);ax.set_title(('A   ' if j==0 else '')+study+' · fitted parameters',loc='left',fontsize=12)
 missing=[]
 for g in COL:
  rr=[r for r in surfaces if r['study']==study and r['group']==g];x=np.array([float(r['tau']) for r in rr]);lo=np.array([float(r['rms_low']) for r in rr]);hi=np.array([float(r['rms_high']) for r in rr]);ok=np.array([r['compatible']=='True' for r in rr])
  if ok.any():ax.fill_between(x,lo,hi,where=ok,color=COL[g],alpha=.23)
  else:missing.append(g)
  f=F[study,g];ax.plot(float(f['tau']),float(f['rms_velocity']),'*',color=COL[g],ms=12)
 ax.set(xscale='log',yscale='log',xlabel='Velocity memory tau (h)',ylabel='Stationary RMS speed (µm/h)',xlim=(.0004,48))
 ax.text(0,-.29,'Stars: optimal fits. Shading: both fitting intervals satisfied.\nReading tolerances, not confidence regions.'+('\nNo compatible region: '+', '.join(missing) if missing else ''),transform=ax.transAxes,va='top',fontsize=9,color='#555555')

ax=fig.add_axes([.745,.762,.21,.14]);ax.set_title('B   Which parameters must differ?',loc='left',fontsize=12)
mat=np.array([[float(next(r['score'] for r in comp if r['study']==s and int(r['mask'])==m)) for m in range(4)] for s in ('Xuan2018','Mallin2023')]);ax.imshow(np.log1p(mat),aspect='auto',cmap='YlOrRd')
for i in range(2):
 for j in range(4):ax.text(j,i,f'{mat[i,j]:.1f}',ha='center',va='center',fontsize=7.5)
ax.set(xticks=range(4),xticklabels=['None','tau','η','Both'],yticks=[0,1],yticklabels=['Xuan','Mallin']);ax.tick_params(length=0,labelsize=8);ax.set_frame_on(False)
ax.text(0,-.29,'Joint reference/polyploid refits. Columns name free differences.\nUnlisted parameters shared and jointly optimized.\nScores are reading-scaled distances, not likelihoods.',transform=ax.transAxes,fontsize=9,color='#555555',va='top')

fig.text(.05,.690,'C   Fitted parameters and equivalent RMS speed · two free quantities per population',fontsize=14,weight='bold')
for j,(key,title,unit) in enumerate([('tau','Velocity memory','h'),('noise_coefficient','Wiener noise coefficient η','µm / h^(3/2)'),('rms_velocity','Equivalent stationary RMS speed','µm/h')]):
 ax=fig.add_axes([.075+j*.31,.578,.24,.090])
 for gi,g in enumerate(COL):
  vals=[float(F[s,g][key]) for s in ('Xuan2018','Mallin2023')];bars=ax.bar(np.arange(2)+(gi-.5)*.30,vals,width=.30,color=COL[g])
  for b,v in zip(bars,vals):ax.annotate(f'{v:.3g}'+(' †' if j<3 and g=='reference' and b.get_x()>.5 else ''),(b.get_x()+b.get_width()/2,v),xytext=(0,3),textcoords='offset points',ha='center',fontsize=9)
 ax.set(xticks=[0,1],xticklabels=['Xuan','Mallin'],ylabel=unit,ylim=(0,None));ax.margins(y=.25);ax.set_title(title,loc='left')
fig.text(.075,.552,'RMS speed is derived from tau and η, not a third fitted parameter. † Mallin reference is at the tau lower bound; its tau, η and RMS speed are unresolved.',fontsize=10,color='#555555')

fig.text(.05,.523,'D   Original fitting targets and the Mallin endpoint check',fontsize=14,weight='bold')
for j,(study,metric,title) in enumerate([('Xuan2018','speed','Xuan sampled speed'),('Xuan2018','msd','Xuan late MSD'),('Mallin2023','path_length','Mallin path length'),('Mallin2023','directness','Mallin directness'),('Mallin2023','endpoint','Mallin endpoint · check')]):
 ax=fig.add_axes([.075+j*.184,.422,.14,.077]);rr=[r for r in pred if r['study']==study and r['metric']==metric and (r['role']=='fit' or metric=='endpoint')]
 for i,g in enumerate(COL):
  r=next(r for r in rr if r['group']==g);v=float(r['observed']);lo=float(r['low']);hi=float(r['high']);ax.errorbar(i-.1,v,yerr=[[v-lo],[hi-v]],fmt='D',mfc='white',color='black',ms=5,capsize=3);ax.plot(i+.1,float(r['predicted']),'o',color=COL[g],ms=6)
 ax.set(xticks=[0,1],xticklabels=['Ref','Poly'],xlim=(-.5,1.5),ylim=(0,None));ax.set_title(title,fontsize=11,loc='left');ax.grid(axis='y',alpha=.1)
fig.text(.075,.395,'Open diamonds/bars: literature means and reading intervals (not CIs). Coloured dots: exact Gaussian means; directness uses 65,536 independent simulated cells. No whole MSD curve was fitted.',fontsize=10,color='#555555')

ax=fig.add_axes([.075,.234,.25,.112]);ax.set_title('E   Xuan MSD, early-lag check and crossover',loc='left')
for g in COL:
 rr=[r for r in times if r['study']=='Xuan2018' and r['group']==g];t=np.array([float(r['time']) for r in rr]);y=np.array([float(r['msd']) for r in rr]);keep=t<=12;ax.plot(t[keep],y[keep],color=COL[g],lw=2)
 for r in [r for r in pred if r['study']=='Xuan2018' and r['group']==g and r['metric']=='msd']:
  lag=11.1111111111 if r['role']=='fit' else 1/6;v=float(r['observed']);ax.errorbar(lag,v,yerr=[[v-float(r['low'])],[float(r['high'])-v]],fmt='D',mfc='white',color=COL[g],ms=4)
ax.axvspan(2,3,color='#bbbbbb',alpha=.2);ax.set(xscale='log',yscale='log',xlabel='Lag (h)',ylabel='MSD (µm²)')
ax.text(0,-.31,'Grey: reported crossover interval. Early abscissa is only\napproximate (600–1000 s). Shared tau gives proportional\nMSD curves, so changing noise strength alone cannot cross.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')

ax=fig.add_axes([.40,.234,.25,.112]);ax.set_title('F   Xuan slow fraction: threshold sensitivity',loc='left')
for g in COL:
 for mode,style in [('net_displacement','-'),('sampled_path','--')]:
  rr=[r for r in slow if r['group']==g and r['interpretation']==mode];ax.plot([float(r['threshold_um']) for r in rr],[float(r['slow_fraction']) for r in rr],color=COL[g],ls=style,lw=1.4)
 lo,hi=(.4,.5) if g=='reference' else (.65,.74);ax.axhspan(lo,hi,color=COL[g],alpha=.12)
ax.set(xlim=(0,60),ylim=(0,1),xlabel='Common one-hour displacement threshold (µm)',ylabel='Fraction of time classified slow')
ax.text(0,-.31,'Shading: reported slow fractions. Solid: net displacement;\ndashed: sampled path. Cell length and exact implementation\nare unknown. No group-specific length threshold is assumed.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')

ax=fig.add_axes([.735,.234,.22,.112]);ax.set_title('G   Jemaà speed-only families',loc='left')
for study,style in [('Jemaa2023','-'),('Jemaa2017','--')]:
 for g in COL:
  rr=[r for r in surfaces if r['study']==study and r['group']==g];x=[float(r['tau']) for r in rr];ax.plot(x,[float(r['rms']) for r in rr],color=COL[g],ls=style,lw=1.4);ax.fill_between(x,[float(r['rms_low']) for r in rr],[float(r['rms_high']) for r in rr],color=COL[g],alpha=.1)
ax.set(xscale='log',yscale='log',xlabel='Velocity memory tau (h)',ylabel='Stationary RMS speed (µm/h)')
ax.text(0,-.31,'Solid: 2023 RKO; dashed: 2017 sarcoma.\nSpeed alone constrains a curve, not a unique point.\nCamera cadence is included in the mean speed.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')

fig.text(.05,.143,'H   Homogeneous simulated trajectories · short tracks are prefixes of long tracks',fontsize=14,weight='bold')
for col,study in enumerate(('Xuan2018','Mallin2023')):
 for sub,duration in enumerate((1,12 if study=='Xuan2018' else 24)):
  ax=fig.add_axes([.075+(col*2+sub)*.23,.039,.155,.081]);cad=1/6 if study=='Xuan2018' else .5
  limit=max(10,float(np.ceil(max(np.abs(tracks[study+'_'+g][:12,:round(duration/cad)+1]).max() for g in COL)/10)*10))
  for g in COL:
   for p in tracks[study+'_'+g][:12,:round(duration/cad)+1]:ax.plot(p[:,0],p[:,1],color=COL[g],alpha=.65,lw=.7)
  ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)');ax.set_title(f'{study} · {duration} h',loc='left',fontsize=11);ax.tick_params(labelsize=8)
fig.text(.05,.010,'Gaussian OU for these studies only; Hosny remains unchanged. Previous models are archived. Ambiguous thresholds, figure-read means and finite-track sampling limit inference.',fontsize=9.5,color='#555555')
fig.savefig(OUT/'literature_gaussian_overview.png',dpi=180);print(OUT/'literature_gaussian_overview.png')
