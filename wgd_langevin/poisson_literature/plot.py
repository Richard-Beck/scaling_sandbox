from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
def read(n):return list(csv.DictReader((OUT/n).open()))
fits=read('fits.csv');surfaces=read('parameter_surfaces.csv');pred=read('predictions.csv');times=read('time_predictions.csv');slow=read('slow_fraction_sensitivity.csv');profiles=read('rate_profiles.csv');comp=read('model_comparison.csv');tracks=np.load(OUT/'tracks.npz');checks=json.load((OUT/'checks.json').open())
F={(f['study'],f['group']):f for f in fits};COL={'reference':'#286FA7','polyploid':'#C5522E'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':12,'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(24,25),facecolor='white')
fig.text(.05,.978,'The earlier WGD literature, revisited with the same Poisson-kick velocity model',fontsize=23,weight='bold',va='top')
fig.text(.05,.958,'Same simulator as Hosny: velocity decays over tau; kicks arrive at rate λ with Rayleigh magnitude scale A. Uniform-medium assays use κ = 0; all cells share parameters.',fontsize=11.5,color='#555555')
fig.text(.05,.942,'Original fitted summaries and reading intervals retained. Reference = parental/diploid comparator; polyploid = PGCC/PACC or selected tetraploid clone, depending on the study.',fontsize=10.5,color='#555555')
fig.legend(handles=[Line2D([],[],color=COL[g],lw=3,label=l) for g,l in [('reference','Reference'),('polyploid','Polyploid')]],ncol=2,bbox_to_anchor=(.96,.935),loc='upper right',frameon=False)

# A: full three-parameter ambiguity, rather than a unique propulsion/persistence point.
for j,study in enumerate(('Xuan2018','Mallin2023')):
 ax=fig.add_axes([.075+j*.31,.762,.24,.14]);ax.set_title(('A   ' if j==0 else '')+study+' · compatible parameter families',loc='left',fontsize=12)
 missing=[]
 for g in COL:
  rr=[r for r in surfaces if r['study']==study and r['group']==g];tx=np.unique([float(r['tau']) for r in rr]);ry=np.unique([float(r['rate']) for r in rr]);mask=np.array([[next(r['compatible']=='True' for r in rr if float(r['tau'])==t and float(r['rate'])==v) for t in tx] for v in ry])
  if mask.any():ax.contourf(tx,ry,mask,levels=[.5,1.5],colors=[COL[g]],alpha=.17);ax.contour(tx,ry,mask,levels=[.5],colors=[COL[g]],linewidths=1)
  else:missing.append(g)
  f=F[study,g];ax.plot(float(f['tau']),float(f['rate']),'*',color=COL[g],ms=12)
 ax.set(xscale='log',yscale='log',xlabel='Velocity memory tau (h)',ylabel='Kick rate λ (/h)',xlim=(tx[0],tx[-1]),ylim=(ry[0],ry[-1]))
 ax.text(0,-.29,'Shading: some A meets every fitting interval on the grid.\nStars: illustrative numerical minima, not unique parameter estimates.'+ ('\nNo compatible grid region: '+', '.join(missing) if missing else ''),transform=ax.transAxes,va='top',fontsize=9,color='#555555')

ax=fig.add_axes([.745,.762,.21,.14]);ax.set_title('B   Which parameters must differ?',loc='left',fontsize=12)
mat=np.array([[float(next(r['score'] for r in comp if r['study']==s and int(r['mask'])==m)) for m in range(8)] for s in ('Xuan2018','Mallin2023')]);ax.imshow(np.log1p(mat),aspect='auto',cmap='YlOrRd')
for i in range(2):
 for j in range(8):ax.text(j,i,f'{mat[i,j]:.1f}',ha='center',va='center',fontsize=7.5)
ax.set(xticks=range(8),xticklabels=['None','tau','λ','tau\nλ','A','tau\nA','λ\nA','All'],yticks=[0,1],yticklabels=['Xuan','Mallin']);ax.tick_params(length=0,labelsize=8);ax.set_frame_on(False)
ax.text(0,-.29,'Joint reference/polyploid refits. Columns name free differences.\nOther parameters shared; A is reoptimized jointly.\nScores are reading-scaled distances, not likelihoods.',transform=ax.transAxes,fontsize=9,color='#555555',va='top')

fig.text(.05,.690,'C   Illustrative parameter values · three free quantities per population, often only two fitted summaries',fontsize=14,weight='bold')
for j,(key,title,unit) in enumerate([('tau','Velocity memory','h'),('rate','Kick rate','events/h'),('amplitude','Kick magnitude scale','µm/h')]):
 ax=fig.add_axes([.075+j*.31,.578,.24,.090])
 for gi,g in enumerate(COL):
  vals=[float(F[s,g][key]) for s in ('Xuan2018','Mallin2023')];bars=ax.bar(np.arange(2)+(gi-.5)*.30,vals,width=.30,color=COL[g])
  for b,v in zip(bars,vals):ax.annotate(f'{v:.3g}',(b.get_x()+b.get_width()/2,v),xytext=(0,3),textcoords='offset points',ha='center',fontsize=9)
 ax.set(xticks=[0,1],xticklabels=['Xuan','Mallin'],ylabel=unit,ylim=(0,None));ax.margins(y=.25);ax.set_title(title,loc='left')
fig.text(.075,.552,'Directional concentration κ = 0 by assay symmetry, not inferred from directional observations. Large A at short tau is not a measured instantaneous cell speed.',fontsize=10,color='#555555')

fig.text(.05,.523,'D   Original fitting targets and the Mallin endpoint check',fontsize=14,weight='bold')
for j,(study,metric,title) in enumerate([('Xuan2018','speed','Xuan sampled speed'),('Xuan2018','msd','Xuan late MSD'),('Mallin2023','path_length','Mallin path length'),('Mallin2023','directness','Mallin directness'),('Mallin2023','endpoint','Mallin endpoint · check')]):
 ax=fig.add_axes([.075+j*.184,.422,.14,.077]);rr=[r for r in pred if r['study']==study and r['metric']==metric and (r['role']=='fit' or metric=='endpoint')]
 for i,g in enumerate(COL):
  r=next(r for r in rr if r['group']==g);v=float(r['observed']);lo=float(r['low']);hi=float(r['high']);ax.errorbar(i-.1,v,yerr=[[v-lo],[hi-v]],fmt='D',mfc='white',color='black',ms=5,capsize=3);ax.plot(i+.1,float(r['predicted']),'o',color=COL[g],ms=6)
 ax.set(xticks=[0,1],xticklabels=['Ref','Poly'],xlim=(-.5,1.5),ylim=(0,None));ax.set_title(title,fontsize=11,loc='left');ax.grid(axis='y',alpha=.1)
fig.text(.075,.395,'Open diamonds/bars: literature means and reading intervals (not CIs). Coloured dots: independent 65,536-cell means; MSD is analytic. No whole MSD curve was fitted.',fontsize=10,color='#555555')

ax=fig.add_axes([.075,.234,.25,.112]);ax.set_title('E   Xuan MSD, early-lag check and crossover',loc='left')
for g in COL:
 rr=[r for r in times if r['study']=='Xuan2018' and r['group']==g];t=np.array([float(r['time']) for r in rr]);y=np.array([float(r['msd']) for r in rr]);keep=t<=12;ax.plot(t[keep],y[keep],color=COL[g],lw=2)
 for r in [r for r in pred if r['study']=='Xuan2018' and r['group']==g and r['metric']=='msd']:
  lag=11.1111111111 if r['role']=='fit' else 1/6;v=float(r['observed']);ax.errorbar(lag,v,yerr=[[v-float(r['low'])],[float(r['high'])-v]],fmt='D',mfc='white',color=COL[g],ms=4)
ax.axvspan(2,3,color='#bbbbbb',alpha=.2);ax.set(xscale='log',yscale='log',xlabel='Lag (h)',ylabel='MSD (µm²)')
ax.text(0,-.31,'Grey: reported crossover interval. Early abscissa is only\napproximate (600–1000 s). Shared tau gives proportional\nMSD curves, so changing rate/amplitude alone cannot cross.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')

ax=fig.add_axes([.40,.234,.25,.112]);ax.set_title('F   Xuan slow fraction: threshold sensitivity',loc='left')
for g in COL:
 for mode,style in [('net_displacement','-'),('sampled_path','--')]:
  rr=[r for r in slow if r['group']==g and r['interpretation']==mode];ax.plot([float(r['threshold_um']) for r in rr],[float(r['slow_fraction']) for r in rr],color=COL[g],ls=style,lw=1.4)
 lo,hi=(.4,.5) if g=='reference' else (.65,.74);ax.axhspan(lo,hi,color=COL[g],alpha=.12)
ax.set(xlim=(0,60),ylim=(0,1),xlabel='Common one-hour displacement threshold (µm)',ylabel='Fraction of time classified slow')
ax.text(0,-.31,'Shading: reported slow fractions. Solid: net displacement;\ndashed: sampled path. Cell length and exact implementation\nare unknown. No group-specific length threshold is assumed.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')

ax=fig.add_axes([.735,.234,.22,.112]);ax.set_title('G   Jemaà speed-only families',loc='left')
for study,style in [('Jemaa2023','-'),('Jemaa2017','--')]:
 rates=np.unique([float(r['rate']) for r in surfaces if r['study']==study]);rate=rates[np.argmin(abs(np.log(rates)))];
 for g in COL:
  rr=sorted([r for r in surfaces if r['study']==study and r['group']==g and float(r['rate'])==rate],key=lambda r:float(r['tau']));x=[float(r['tau']) for r in rr];ax.plot(x,[float(r['amplitude']) for r in rr],color=COL[g],ls=style,lw=1.4);ax.fill_between(x,[float(r['amplitude_low']) for r in rr],[float(r['amplitude_high']) for r in rr],color=COL[g],alpha=.1)
ax.set(xscale='log',yscale='log',xlabel='Velocity memory tau (h)',ylabel='Kick magnitude scale A (µm/h)')
ax.text(0,-.31,f'Slice at λ = {rate:.2g}/h, not an estimated rate.\nSolid: 2023 RKO; dashed: 2017 sarcoma.\nEvery rate generates another viable family.',transform=ax.transAxes,va='top',fontsize=9,color='#555555')

fig.text(.05,.143,'H   Homogeneous simulated trajectories · short tracks are prefixes of long tracks',fontsize=14,weight='bold')
for col,study in enumerate(('Xuan2018','Mallin2023')):
 for sub,duration in enumerate((1,12 if study=='Xuan2018' else 24)):
  ax=fig.add_axes([.075+(col*2+sub)*.23,.039,.155,.081]);cad=1/6 if study=='Xuan2018' else .5
  limit=max(10,float(np.ceil(max(np.abs(tracks[study+'_'+g][:12,:round(duration/cad)+1]).max() for g in COL)/10)*10))
  for g in COL:
   for p in tracks[study+'_'+g][:12,:round(duration/cad)+1]:ax.plot(p[:,0],p[:,1],color=COL[g],alpha=.65,lw=.7)
  ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)');ax.set_title(f'{study} · {duration} h',loc='left',fontsize=11);ax.tick_params(labelsize=8)
fig.text(.05,.010,'Same Hosny Poisson simulator; κ = 0; no permanent heterogeneity. Previous constant-speed results are archived. Ambiguous thresholds, figure-read means and finite-track sampling limit inference.',fontsize=9.5,color='#555555')
fig.savefig(OUT/'literature_poisson_overview.png',dpi=180);print(OUT/'literature_poisson_overview.png')
