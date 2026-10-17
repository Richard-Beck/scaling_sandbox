from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
def read(name):return list(csv.DictReader((OUT/name).open()))
allm=read('all_shared_models.csv');front=read('shared_frontier.csv');ratios=read('frontier_ratios.csv');profiles=read('ratio_profiles.csv');pars=read('representative_parameters.csv')
base={p['group']:p for p in pars if p['group'] in ('N','P')};NC='#286FA7';PC='#C5522E';RC='#6b4ca5';GOLD='#d49b19'
PARAMS=[('tau','Velocity memory, tau','h'),('rate','Kick rate, λ','events/h'),('amplitude','Kick magnitude scale, A','µm/h'),('kappa','Kick directionality, κ','dimensionless')]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'axes.titlesize':13,'xtick.labelsize':10,'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(24,21),facecolor='white')
fig.text(.05,.98,'PACC versus non-PACC: which fitted parameter differences survive model choice?',fontsize=23,weight='bold',va='top')
fig.text(.05,.956,'HIF and RhoA must change the same parameter subsets in both populations. Parameter values and change magnitudes remain population-specific.',fontsize=12,color='#444444')
fig.text(.05,.939,'Representative model: HIF changes directionality; RhoA changes directionality. This is the best simplest nonzero shared mechanism, chosen for presentation—not as a unique explanation.',fontsize=11.5,color='#444444')

ax=fig.add_axes([.075,.77,.39,.135]);ax.set_title('A   Shared-mechanism Pareto frontier',loc='left',weight='bold',fontsize=15)
x=np.array([int(r['complexity_per_population']) for r in allm]);y=np.array([float(r['score']) for r in allm]);rng=np.random.default_rng(613);dx=rng.uniform(-.085,.085,len(x));dx[x==0]=0
ax.scatter(x+dx,y,s=15,color='#777777',alpha=.25,label='All 225 shared subset pairs + null')
ff=sorted(front,key=lambda r:int(r['complexity_per_population']));ax.plot([int(r['complexity_per_population']) for r in ff],[float(r['score']) for r in ff],'o-',color=RC,lw=2,ms=5,label='Pareto frontier')
rep=next(r for r in front if int(r['complexity_per_population'])==2);ax.plot(2,float(rep['score']),'o',ms=12,mfc='#fff5d8',mec=GOLD,mew=2,label='Representative: HIF κ, RhoA κ')
ax.set(xticks=[0,2,3,4,5,6,7,8],xlabel='KD-specific parameters per population, k',ylabel='Combined independent-simulation error',ylim=(0,.067),xlim=(-.3,8.3));ax.grid(axis='y',alpha=.12);ax.legend(frameon=False,fontsize=9,loc='upper right')
ax.text(0,-.29,'Total parameters across both populations = 8 baseline + 2k KD parameters.\nThe score weights the ten source datasets equally; it is not a likelihood.',transform=ax.transAxes,fontsize=9.5,color='#555555',va='top')

ax=fig.add_axes([.55,.761,.41,.145]);ax.axis('off');ax.text(0,1.04,'B   Mechanisms on the shared frontier',fontsize=15,weight='bold')
def sym(x):return x.replace('tau','τ').replace('lambda','λ').replace('kappa','κ').replace('+',' + ')
rows=[[r['complexity_per_population'],sym(r['HIF']),sym(r['RhoA']),f"{float(r['score']):.4f}"] for r in ff]
t=ax.table(cellText=rows,colLabels=['k','HIF changes','RhoA changes','Error'],cellLoc='center',colWidths=[.08,.36,.38,.18],bbox=[0,.12,1,.85]);t.auto_set_font_size(False);t.set_fontsize(10)
for (i,j),cell in t.get_celld().items():
    cell.set_edgecolor('#dddddd');cell.set_linewidth(.4)
    if i==0:cell.set_facecolor('#eeeeee');cell.set_text_props(weight='bold')
    elif rows[i-1][0]=='2':cell.set_facecolor('#fff5dc')
ax.text(0,.025,'Same mechanism means the same parameters can vary, not equal KD effect sizes.\nNo parameter averaging across models; no formal preferred complexity inferred.',fontsize=9.5,color='#555555',va='top')

fig.text(.05,.696,'C   Representative model: the four fitted baseline parameters',fontsize=15,weight='bold')
fig.text(.05,.678,'Blue: non-PACC; orange: PACC. These are conditional point estimates. Panel E checks how much they can change while maintaining a similar fit.',fontsize=11,color='#555555')
for j,(key,title,unit) in enumerate(PARAMS):
    ax=fig.add_axes([.075+j*.232,.553,.18,.105]);values=[float(base[g][key]) for g in ('N','P')]
    ax.plot([0,1],values,color='#bbbbbb',lw=1.1,zorder=1)
    for i,(v,c) in enumerate(zip(values,(NC,PC))):ax.plot(i,v,'o',ms=8,color=c);ax.annotate(f'{v:.3g}',(i,v),xytext=(8,7),textcoords='offset points',fontsize=11,color=c)
    lo,hi=min(values),max(values);span=max(hi-lo,hi*.3);ax.set(xlim=(-.4,1.7),ylim=(max(0,lo-span*.8),hi+span),xticks=[0,1],xticklabels=['non-PACC','PACC'],ylabel=unit)
    ax.set_title(title,loc='left',fontsize=12);ax.text(.02,.94,f'PACC/non-PACC = {values[1]/values[0]:.2f}',transform=ax.transAxes,va='top',fontsize=10,weight='bold',color='#555555')
fig.text(.075,.523,'A is a velocity-jump scale, not average cell speed. Memory is well below the 30-minute observation interval in this representative fit.',fontsize=10,color='#555555')

fig.text(.05,.490,'D   Between-model sensitivity: baseline PACC/non-PACC ratios across the nonzero shared frontier',fontsize=15,weight='bold')
for j,(key,title,unit) in enumerate(PARAMS):
    ax=fig.add_axes([.075+j*.232,.353,.18,.107]);rr=sorted([r for r in ratios if r['parameter']==key and int(r['complexity'])>0],key=lambda r:int(r['complexity']));x=np.array([int(r['complexity']) for r in rr]);y=np.array([float(r['ratio']) for r in rr])
    ax.axhline(1,color='#888888',ls='--',lw=1);ax.plot(x,y,'o-',color=RC,lw=1.6,ms=5);ax.plot(x[0],y[0],'o',mfc='#fff5d8',mec=GOLD,mew=2,ms=9)
    low=min(.8,y.min()*.65);high=max(1.3,y.max()*1.5)
    ax.set(xticks=[2,3,4,5,6,7,8],xlim=(1.7,8.3),yscale='log',ylim=(low,high),xlabel='KD parameters per population, k',ylabel='Baseline PACC / non-PACC');ax.set_title(title,loc='left',fontsize=12)
    ticks={'tau':[1,2,5,10,20],'rate':[.5,.7,1],'amplitude':[.1,.2,.5,1],'kappa':[1,1.5,2,3]}[key]
    ax.set_yticks(ticks);ax.set_yticklabels([f'{v:g}' for v in ticks]);ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.text(.02,.95,f'Range: {y.min():.2f}–{y.max():.2f}',transform=ax.transAxes,va='top',fontsize=10,color='#555555')
fig.text(.075,.319,'Dashed line: equal values. Gold marker: the representative model. A consistent sign across these point estimates is not, by itself, proof of an identified difference.',fontsize=10,color='#555555')

fig.text(.05,.284,'E   Within-model sensitivity: fix each baseline ratio, then jointly refit everything else',fontsize=15,weight='bold')
for j,(key,title,unit) in enumerate(PARAMS):
    ax=fig.add_axes([.075+j*.232,.132,.18,.116]);rr=sorted([r for r in profiles if r['parameter']==key],key=lambda r:float(r['ratio']));x=np.array([float(r['ratio']) for r in rr]);y=np.array([float(r['excess']) for r in rr]);rep_ratio=float(base['P'][key])/float(base['N'][key])
    ax.axvline(1,color='#888888',ls='--',lw=1);ax.axhline(0,color='#bbbbbb',lw=.8);ax.plot(x,y,'o-',color=RC,lw=1.3,ms=3);ax.plot(rep_ratio,0,'*',color=GOLD,ms=12,zorder=5)
    ax.set(xscale='log',yscale='symlog',xlim=(.025,36),ylim=(-.0002,.045),xlabel='Fixed baseline PACC / non-PACC',ylabel='Increase in combined fitting error');ax.set_yscale('symlog',linthresh=.002);ax.set_yticks([0,.001,.01,.04]);ax.set_yticklabels(['0','.001','.01','.04'])
    ax.set_xticks([.03,.1,1,10,30]);ax.set_xticklabels(['.03','.1','1','10','30']);ax.set_title(title,loc='left',fontsize=12)
    eq=float(next(r['excess'] for r in rr if np.isclose(float(r['ratio']),1)));ax.text(0,-.36,f'Equal-value penalty: {eq:.4f}',transform=ax.transAxes,fontsize=10,color='#555555')
fig.text(.075,.055,'Profiles use the representative HIF κ / RhoA κ mechanism only. All other baseline and KD parameters are refitted jointly; tau and A may compensate for one another.',fontsize=10,color='#555555')
fig.text(.075,.038,'All four profiles use the same vertical scale and the optimization simulation bank. Flat curves indicate weak identification; these are not confidence intervals.',fontsize=10,color='#555555')
fig.text(.05,.015,'Source: supplied hosny_2026.md only. Homogeneous stationary Poisson-kick populations; exact event integration. Directional-index normalization remains ambiguous; no spatial oxygen response or cell-area model.',fontsize=9.5,color='#555555')
fig.savefig(OUT/'hosny_shared_parameter_comparison.png',dpi=180);print(OUT/'hosny_shared_parameter_comparison.png')
