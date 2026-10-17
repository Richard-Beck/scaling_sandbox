from pathlib import Path
import sys, csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, MaxNLocator
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from model import simulate, metrics
fits=json.loads((ROOT/'results/summary.json').read_text())
obs=list(csv.DictReader((ROOT/'data/observations.csv').open()))
for r in obs:
    for k in ('value','low','high','lag_h'):
        r[k]=float(r[k]) if r[k] else None
C={'reference':'#286DA8','polyploid':'#C85834'}
G=('reference','polyploid')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,
                     'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'axes.linewidth':.8,'savefig.facecolor':'white'})
paths={};met={};msds={}
for study,cad,duration in [('Xuan2018',1/6,12),('Mallin2023',.5,24)]:
    for group in G:
        f=fits[study][group]['best']
        sub=max(48,int(np.ceil(10*cad/f['tau_h'])))
        p=simulate(f['propulsion_um_h'],f['tau_h'],cad,duration,n=4096,substeps=sub,seed=491031,cv=0)
        paths[study,group]=p;met[study,group]=metrics(p,cad)
        if study=='Xuan2018':
            msds[group]=np.array([np.mean(np.sum((p[:,k:]-p[:,:-k])**2,axis=2)) for k in range(1,p.shape[1])])
print('Simulation complete',flush=True)

fig=plt.figure(figsize=(22,18),facecolor='white')
fig.text(.055,.967,'A common propulsion–persistence map of polyploid cell migration',fontsize=23,weight='bold',va='top')
fig.text(.055,.938,'Existing parameter fits • homogeneous populations • experimental sampling retained • no refitting',fontsize=12,color='#555555')
fig.legend(handles=[Line2D([],[],color=C[g],lw=3,label=l) for g,l in zip(G,('Reference cells','Polyploid cells'))]+[
    Line2D([],[],color='black',marker='D',mfc='white',ls='none',label='Literature estimate'),
    Line2D([],[],color='black',marker='o',ls='none',label='Simulated mean')],
    loc='upper right',bbox_to_anchor=(.975,.941),frameon=False,ncol=4,fontsize=10,handlelength=1.5,columnspacing=1.4)

def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=12)
    ax.text(-.12,1.055,letter,transform=ax.transAxes,fontsize=17,weight='bold',va='bottom')
    ax.tick_params(direction='out',length=4,color='#555555')

def record(study,group,metric):
    return next(r for r in obs if r['study']==study and r['group']==group and r['metric']==metric)

# A: reuse original fitted banks and existing observation intervals for Jemaa.
ax=fig.add_axes([.065,.633,.405,.258])
panel(ax,'A','Shared parameter space')
for study,ls in [('Jemaa2023','-'),('Jemaa2017','--')]:
    bank=list(csv.DictReader((ROOT/f'results/{study}_simulation_bank.csv').open()))
    t=np.array([float(r['tau_h']) for r in bank]);atten=np.array([float(r['speed']) for r in bank])
    for group in G:
        r=next(r for r in obs if r['study']==study and r['group']==group and r['role']=='fit')
        ax.fill_between(t,r['low']/atten,r['high']/atten,color=C[group],alpha=.13)
        ax.plot(t,r['value']/atten,ls=ls,color=C[group],lw=1.5,alpha=.8)
# Indicate Xuan compatibility projections lightly; do not create new fit contours.
for study,marker,sz in [('Xuan2018','X',120),('Mallin2023','+',210)]:
    for group in G:
        f=fits[study][group]['best'];t=f['tau_h'];v=f['propulsion_um_h']
        ax.scatter(t,v,s=sz+100,marker=marker,color='white',linewidths=5,zorder=5)
        ax.scatter(t,v,s=sz,marker=marker,color=C[group],linewidths=2.8,zorder=6)
        offsets={('Xuan2018','reference'):(18,-15),('Xuan2018','polyploid'):(14,-4),
                 ('Mallin2023','reference'):(10,22),('Mallin2023','polyploid'):(-65,-39)}
        name='Xuan' if study=='Xuan2018' else 'Mallin'
        ax.annotate(f'{name}\n{v:.1f} µm/h, {t:.3g} h',(t,v),xytext=offsets[study,group],textcoords='offset points',
                    fontsize=10,color=C[group],weight='normal',zorder=7,
                    bbox=dict(boxstyle='round,pad=.2',fc='white',ec='none',alpha=.85))
for study in ('Xuan2018','Mallin2023'):
    a=fits[study]['reference']['best'];b=fits[study]['polyploid']['best']
    ax.annotate('',xy=(b['tau_h'],b['propulsion_um_h']),xytext=(a['tau_h'],a['propulsion_um_h']),
                arrowprops=dict(arrowstyle='->',color='#888888',lw=1.15,alpha=.7),zorder=3)
ax.set(xscale='log',yscale='log',xlim=(.007,55),ylim=(2,205),xlabel='Persistence time (h)',ylabel='Intrinsic propulsion speed (µm/h)')
ax.set_xticks([.01,.1,1,10]);ax.set_xticklabels(['0.01','0.1','1','10'])
ax.set_yticks([2,5,10,20,50,100,200]);ax.yaxis.set_major_formatter(ScalarFormatter())
ax.legend(handles=[Line2D([],[],color='#555555',ls='-',label='Jemaà 2023 viable band'),
                   Line2D([],[],color='#555555',ls='--',label='Jemaà 2017 viable band'),
                   Line2D([],[],color='#555555',marker='X',ls='none',label='Xuan fitted point'),
                   Line2D([],[],color='#555555',marker='+',ls='none',label='Mallin fitted point')],
          loc='upper right',frameon=True,facecolor='white',edgecolor='none',fontsize=9)
ax.text(0,-.21,'Bands: original mean-speed constraints. Crosses: original homogeneous optima.\nMallin reference is weakly constrained; crosses are not unique parameter estimates.',
        transform=ax.transAxes,fontsize=9.5,color='#555555',va='top')

# B/C: same cells in short and long views. No selection for appearance.
track_rects=[(.555,.778,.172,.115),(.79,.778,.172,.115),(.555,.615,.172,.115),(.79,.615,.172,.115)]
for row,(study,cad,duration,letter,name) in enumerate([('Xuan2018',1/6,12,'B','Xuan'),('Mallin2023',.5,24,'C','Mallin')]):
    for col,time in enumerate((1,duration)):
        ax=fig.add_axes(track_rects[row*2+col]);lim=32 if col==0 else 25*np.ceil(max(np.abs(paths[st,g][:16]).max() for st in ('Xuan2018','Mallin2023') for g in G)/25)
        for group in G:
            p=paths[study,group][:16,:int(round(time/cad))+1]
            for track in p:
                ax.plot(track[:,0],track[:,1],color=C[group],lw=.8,alpha=.55)
                ax.plot(track[-1,0],track[-1,1],'.',color=C[group],ms=2)
        ax.plot(0,0,'k+',ms=5,mew=.9)
        ax.set(xlim=(-lim,lim),ylim=(-lim,lim),aspect='equal')
        ax.set_xticks([-lim+2,0,lim-2] if col==0 else [-lim,0,lim]);ax.set_yticks([-30,0,30] if col==0 else [-lim,0,lim])
        ax.set_title(f'{name} · {time:g} h',loc='left',fontsize=12,pad=5)
        ax.text(-.24,1.07,f'{letter}{col+1}',transform=ax.transAxes,fontsize=15,weight='bold')
        ax.set_xlabel('x (µm)',labelpad=1,fontsize=9);ax.set_ylabel('y (µm)',labelpad=0,fontsize=9)
        ax.tick_params(labelsize=8,length=3)
fig.text(.55,.58,'16 cells per group; short views are the first hour of the same long tracks.\nPosition sampling: Xuan 10 min; Mallin 30 min. Identical spatial scales within each duration column.',fontsize=9.5,color='#555555')

# Xuan: every reported motility summary, independent of fitting status.
fig.text(.055,.546,'XUAN 2018  /  Compare the same fitted realizations with every reported random-migration summary',fontsize=13,weight='bold')
xs=[.065,.302,.539,.776];w=.182;yb=.338;h=.17
ax=fig.add_axes([xs[0],yb,w,h]);panel(ax,'D','Mean speed  ·  used in fit')
for j,g in enumerate(G):
    values=met['Xuan2018',g]['speed'];r=record('Xuan2018',g,'speed')
    vp=ax.violinplot(values,[j-.1],widths=.28,showextrema=False)
    vp['bodies'][0].set_facecolor(C[g]);vp['bodies'][0].set_edgecolor(C[g]);vp['bodies'][0].set_alpha(.23)
    ax.plot(j-.1,values.mean(),'o',color=C[g],ms=7)
    ax.errorbar(j+.13,r['value'],yerr=[[r['value']-r['low']],[r['high']-r['value']]],fmt='D',mfc='white',mec=C[g],ecolor=C[g],capsize=4,ms=7)
ax.set(xticks=[0,1],xticklabels=['Reference','PGCC'],ylabel='Sampled speed (µm/h)',ylim=(0,31),xlim=(-.5,1.5))
ax.text(0,-.22,'Violin: simulated per-cell mean speeds.\nAll cells share one intrinsic speed within a group.',transform=ax.transAxes,fontsize=9,color='#555555',va='top')

ax=fig.add_axes([xs[1],yb,w,h]);panel(ax,'E','MSD across time')
lags=np.arange(1,73)/6
for group in G:
    ax.plot(lags,msds[group],color=C[group],lw=2)
    rr=[r for r in obs if r['study']=='Xuan2018' and r['group']==group and r['metric']=='msd']
    for r in rr:
        if r['role']=='fit':
            ax.errorbar(r['lag_h'],r['value'],yerr=[[r['value']-r['low']],[r['high']-r['value']]],fmt='D',mfc='white',mec=C[group],ecolor=C[group],ms=6,capsize=3)
        else:
            # Local note supplies an early-lag rectangle, not an exact coordinate.
            ax.add_patch(Rectangle((1/6,r['low']),1000/3600-1/6,r['high']-r['low'],facecolor=C[group],edgecolor=C[group],alpha=.19,lw=1))
ax.axvspan(2,3,color='#888888',alpha=.13)
ax.axvline(1.332386584290902,color='#777777',ls=':',lw=1)
ax.set(xscale='log',yscale='log',xlim=(.14,14),ylim=(4,9000),xlabel='Lag (h)',ylabel='MSD (µm²)')
ax.set_xticks([1/6,1,3,12]);ax.set_xticklabels(['1/6','1','3','12'])
ax.text(.04,.98,'Grey: reported crossover ≈2–3 h\nDotted: fitted crossover ≈1.33 h',transform=ax.transAxes,fontsize=8.5,va='top')
ax.text(0,-.22,'Late diamond: fitted. Early boxes: unfitted.\nEarly lag is only localized to 600–1000 s.',transform=ax.transAxes,fontsize=9,color='#555555',va='top')

ax=fig.add_axes([xs[2],yb,w,h]);panel(ax,'F','MSD slope  ·  not fitted')
starts=np.linspace(1/6,2,50)
for group in G:
    r=record('Xuan2018',group,'msd_slope')
    slopes=[]
    for start in starts:
        mask=(lags>=start)&(lags<=6)
        slopes.append(np.polyfit(np.log(lags[mask]),np.log(msds[group][mask]),1)[0])
    ax.plot(starts,slopes,color=C[group],lw=2)
    ax.axhspan(r['low'],r['high'],color=C[group],alpha=.13)
    ax.axhline(r['value'],color=C[group],ls='--',lw=1)
ax.set(xlim=(1/6,2),ylim=(.98,1.8),xlabel='Lower lag in slope fit (h)',ylabel='Log–log MSD slope')
ax.text(.04,.95,'Upper lag fixed at 6 h',transform=ax.transAxes,fontsize=9,va='top')
ax.text(0,-.22,'Bands/dashes: reported slopes. Curves: simulations.\nThe original regression window is unspecified.',transform=ax.transAxes,fontsize=9,color='#555555',va='top')

ax=fig.add_axes([xs[3],yb,w,h]);panel(ax,'G','Slow / fast fraction  ·  not fitted')
L=np.linspace(0,30,601);displacements={};intervals={}
for group in G:
    p=paths['Xuan2018',group]
    d=np.linalg.norm(p[:,6:]-p[:,:-6],axis=2).ravel();displacements[group]=d
    r=record('Xuan2018',group,'slow_fraction')
    ordered=np.sort(d);cdf=np.searchsorted(ordered,L,side='right')/len(d)
    ax.plot(L,cdf,color=C[group],lw=2)
    ax.axhline(r['value'],color=C[group],ls='--',lw=1)
    ax.axhspan(r['low'],r['high'],color=C[group],alpha=.07)
    lo,hi=np.quantile(d,[r['low'],r['high']]);intervals[group]=(lo,hi)
    ax.axvspan(lo,hi,color=C[group],alpha=.17)
    ax.plot(np.quantile(d,r['value']),r['value'],'D',mfc='white',mec=C[group],ms=6,zorder=5)
    print(group,'threshold range for slow fraction',lo,hi,flush=True)
Lref=float(np.quantile(displacements['reference'],.45));polyfrac=float(np.mean(displacements['polyploid']<Lref))
ax.set(xlim=(0,30),ylim=(0,1),xlabel='Assumed cutoff length L (µm)',ylabel='Fraction classified slow')
ax.set_yticks([0,.25,.5,.75,1]);sec=ax.secondary_yaxis('right',functions=(lambda x:1-x,lambda x:1-x));sec.set_ylabel('Fast fraction',fontsize=9);sec.set_yticks([0,.5,1]);sec.tick_params(labelsize=9)
ax.text(.04,.95,'Slow: 1 h net displacement < L',transform=ax.transAxes,fontsize=8.7,va='top')
ax.text(0,-.22,'Horizontal bands: observed slow fractions.\nVertical bands: L needed to match each group.\nDefinition assumed; the numerical cutoff is missing.',transform=ax.transAxes,fontsize=9,color='#555555',va='top')
print('Reference-matched L',Lref,'PGCC predicted slow',polyfrac,flush=True)

# Mallin full distributions, central values, and qualitative tail summaries.
fig.text(.055,.264,'MALLIN 2023  /  Compare the same fitted realizations with path length, displacement and directness',fontsize=13,weight='bold')
for index,(metric,title,ylabel,ylim) in enumerate([
    ('path_length','Accumulated distance  ·  used in fit','Total sampled path length (µm)',(0,1150)),
    ('endpoint','Endpoint displacement  ·  not fitted','Start-to-end distance (µm)',(0,320)),
    ('directness','Directness  ·  used in fit','Endpoint distance / path length',(0,1.02))]):
    ax=fig.add_axes([.065+index*.314,.087,.253,.135]);panel(ax,chr(ord('H')+index),title)
    for j,group in enumerate(G):
        vals=met['Mallin2023',group][metric];r=record('Mallin2023',group,metric)
        vp=ax.violinplot(vals,[j-.15],widths=.32,showextrema=False)
        for b in vp['bodies']:
            b.set_facecolor(C[group]);b.set_edgecolor(C[group]);b.set_alpha(.3)
        ax.plot(j-.15,vals.mean(),'o',color=C[group],ms=7,zorder=5)
        ax.errorbar(j+.15,r['value'],yerr=[[r['value']-r['low']],[r['high']-r['value']]],fmt='D',mfc='white',mec=C[group],ecolor=C[group],capsize=4,ms=7,zorder=6)
        if metric=='path_length':
            lo,hi=(700,750) if group=='reference' else (700,1100)
            ax.plot([j+.15]*2,[lo,hi],color=C[group],lw=5,alpha=.3,solid_capstyle='butt')
            ax.annotate('reported\nhigh tail',(j+.15,(lo+hi)/2),xytext=(14,0),textcoords='offset points',fontsize=8.5,color=C[group],va='center')
        elif metric=='endpoint':
            if group=='polyploid':
                ax.plot([j+.15]*2,[200,280],color=C[group],lw=5,alpha=.3,solid_capstyle='butt')
                ax.annotate('reported\nupper tail',(j+.15,250),xytext=(15,0),textcoords='offset points',fontsize=8.5,color=C[group],va='center')
    ax.set(xticks=[0,1],xticklabels=['Reference','PACC'],xlim=(-.6,1.8),ylim=ylim,ylabel=ylabel)
    if metric=='directness':
        ax.text(.97,.96,'Computed per cell, then averaged',transform=ax.transAxes,fontsize=9,ha='right',va='top')
    notes={'path_length':'Violins: simulated distributions. Diamonds: extracted central values.\nHigh-tail bars are qualitative literature ranges, not error bars.',
           'endpoint':'The homogeneous model predicts the central displacement shift,\nbut does not reproduce all reported extremes.',
           'directness':'Parental directness is slightly overpredicted.\nThe same 30 min / 24 h sampling is applied to both groups.'}
    ax.text(0,-.27,notes[metric],transform=ax.transAxes,fontsize=9,color='#555555',va='top')

fig.text(.055,.026,'Provided-data comparison only: Xuan Fig. 4 and Mallin Fig. 2A–D. Numerical literature values are approximate readings; displayed ranges are not confidence intervals.',fontsize=9.5,color='#555555')
fig.text(.055,.012,'Homogeneous = one propulsion speed and one persistence time per group; trajectories remain stochastic. No translational noise, speed heterogeneity, switching or interaction terms.',fontsize=9.5,color='#555555')
out=ROOT/'results/unified_homogeneous_mockup.png'
fig.savefig(out,dpi=200)
print(out,flush=True)
