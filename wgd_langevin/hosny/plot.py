"""Single standalone PNG: Hosny optimization, validation and predictions."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.ticker import NullFormatter
from data import DATA,BYID,PAIRS,RADIUS,RADIAL_FMI,RADIAL_SE
from model import metrics
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results'
def load(name):return list(csv.DictReader((OUT/name).open()))
fits={r['id']:{k:(float(v) if k in ('tau','a','kappa','v','score') else v) for k,v in r.items()} for r in load('fits.csv')}
pred=load('predictions.csv');comp=load('model_comparison.csv');sens=load('sensitivity.csv')
ray=load('rayleigh_checks.csv');times=load('time_predictions.csv');classes=load('subpopulation_predictions.csv');cf=load('steering_removed.csv')
tracks=np.load(OUT/'illustrative_tracks.npz')
C={'PACC':'#C85632','non-PACC':'#276EA5','normoxia':'#777777'}
colors={r['id']:C[r['group']] for r in DATA}
ids=[r['id'] for r in DATA]
labels={'O':'Normoxia','N':'Hypoxic non-PACC','P':'Hypoxic PACC','HN0':'HIF · non-PACC · ctrl','HN1':'HIF · non-PACC · KD','HP0':'HIF · PACC · ctrl','HP1':'HIF · PACC · KD','RN0':'RhoA · non-PACC · ctrl','RN1':'RhoA · non-PACC · KD','RP0':'RhoA · PACC · ctrl','RP1':'RhoA · PACC · KD'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12.5,'axes.labelsize':10.5,
                     'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(24,25),facecolor='white')
fig.text(.057,.976,'Hosny 2026: separating propulsion, persistence and oxygen-directed steering',fontsize=23,weight='bold',va='top')
fig.text(.057,.953,'11 conditions • 16 h trajectories sampled every 30 min • homogeneous active Langevin cells • supplied literature notes only',fontsize=12,color='#555555')
fig.legend(handles=[Line2D([],[],color=c,lw=3,label=g) for g,c in C.items()]+[
    Line2D([],[],color='black',marker='D',mfc='white',ls='-',label='Reported median + IQR'),
    Line2D([],[],color='black',marker='o',ls='-',label='Simulated median + IQR')],
    ncol=5,loc='upper left',bbox_to_anchor=(.055,.944),frameon=False,fontsize=10)

def panel(rect,letter,title):
    ax=fig.add_axes(rect);ax.set_title(title,loc='left',pad=11)
    ax.text(-.13,1.055,letter,transform=ax.transAxes,fontsize=17,weight='bold',va='bottom')
    return ax

def note(ax,text,y=-.2):ax.text(0,y,text,transform=ax.transAxes,fontsize=9,color='#555555',va='top')

def simq(id,metric):
    rows=sorted([r for r in pred if r['id']==id and r['metric']==metric],key=lambda r:float(r['quantile']))
    return np.array([float(r['predicted']) for r in rows])

# A/B parameter estimates with paired perturbation arrows.
for rect,letter,metric,title,ylabel in [([.075,.806,.255,.107],'A','v','Propulsion–persistence map','Intrinsic propulsion (µm/h)'),
                                       ([.39,.806,.22,.107],'B','kappa','Steering is fitted separately','Steering strength κ (per h)')]:
    ax=panel(rect,letter,title)
    for a,b,name in PAIRS[1:]:
        ax.annotate('',xy=(fits[b]['tau'],fits[b][metric]),xytext=(fits[a]['tau'],fits[a][metric]),
                    arrowprops=dict(arrowstyle='->',lw=1.3,color=colors[a],linestyle='-' if a.startswith('H') else '--',alpha=.7))
    for id in ids:
        f=fits[id];base=id in ('O','N','P');mk='*' if base else ('v' if id.endswith('1') else 'o')
        ax.scatter(f['tau'],f[metric],s=130 if base else 40,marker=mk,color=colors[id],facecolors=colors[id] if base or id.endswith('1') else 'white',zorder=4)
        if base:ax.annotate(id,(f['tau'],f[metric]),xytext=(6,5),textcoords='offset points',color=colors[id],weight='bold')
    ax.set(xscale='log',xlabel='Intrinsic persistence τ (h)',ylabel=ylabel,xlim=(.12,3.8))
    ax.set_xticks([.2,.5,1,2,3]);ax.set_xticklabels(['0.2','0.5','1','2','3']);ax.xaxis.set_minor_formatter(NullFormatter())
    if metric=='v':ax.set(ylim=(4,24))
    else:ax.set(ylim=(-.01,.235))
    note(ax,'Stars: baseline O / N / P. Open circles → triangles: control → KD.\nSolid arrows: HIF-1α; dashed: RhoA. Normoxic steering fixed to zero.',y=-.27)

ax=panel([.727,.806,.233,.107],'C','Which parameter changes improve the fit?')
models=['speed_only','speed_and_steering','speed_and_persistence','all_three']
mat=np.array([[float(next(r['score_per_metric'] for r in comp if r['pair']==name and r['model']==m)) for m in models] for a,b,name in PAIRS])
# Relative score increases keep the baseline irreducible error visible in text.
delta=mat-mat[:,-1,None]
im=ax.imshow(delta,aspect='auto',cmap='YlOrRd',vmin=0,vmax=.05)
for i in range(5):
    for j in range(4):ax.text(j,i,f'{mat[i,j]:.3f}',ha='center',va='center',fontsize=10,color='black')
ax.set_yticks(range(5));ax.set_yticklabels(['Baseline P/N','HIF · non-PACC','HIF · PACC','RhoA · non-PACC','RhoA · PACC'],fontsize=9)
ax.set_xticks(range(4));ax.set_xticklabels(['v only','v + κ','v + τ','v + τ + κ'],fontsize=10)
ax.tick_params(length=0);ax.set_frame_on(False)
note(ax,'Numbers: mean scaled quantile error; lower is better.\nColour: excess over all-three fit. No likelihood or significance test.',y=-.27)

# D five metrics including reserved observables. All groups, no pooling.
fig.text(.057,.747,'D   Observations versus independent realizations of the fitted model',fontsize=16,weight='bold')
fig.text(.057,.731,'Bars are population quartiles, not confidence intervals. Endpoints and perpendicular FMI were excluded from optimization.',fontsize=10.5,color='#555555')
metric_info=[('path','Accumulated distance','µm',(0,315)),('directness','Directness','endpoint / path',(0,.85)),('parallel','Parallel FMI','projection / path',(-.4,.65)),('endpoint','Net displacement','µm',(0,145)),('perpendicular','Perpendicular FMI','projection / path',(-.45,.45))]
for j,(metric,title,xlabel,xlim) in enumerate(metric_info):
    ax=fig.add_axes([.155+j*.163,.516,.137,.193])
    ax.set_title(title+'\n'+('FITTED' if j<3 else 'PREDICTED CHECK'),fontsize=11.5,loc='left',pad=12)
    for i,id in enumerate(ids):
        q=np.array(BYID[id][metric]);s=simq(id,metric);c=colors[id]
        ax.plot([s[0],s[2]],[i+.13,i+.13],color=c,lw=2.3)
        ax.plot(s[1],i+.13,'o',color=c,ms=4.5)
        if np.isfinite(q[[0,2]]).all():ax.plot([q[0],q[2]],[i-.13,i-.13],color='#444444',lw=1.1)
        ax.plot(q[1],i-.13,'D',mfc='white',mec='#333333',ms=4.5)
    ax.set(ylim=(10.7,-.8),xlim=xlim,xlabel=xlabel)
    ax.set_yticks(range(11));ax.set_yticklabels([labels[id] for id in ids] if j==0 else ['']*11,fontsize=9.5)
    for sep in (2.5,6.5):ax.axhline(sep,color='#d6d6d6',lw=.8)
    if metric in ('parallel','perpendicular'):ax.axvline(0,color='#cccccc',ls=':',lw=.8)
    ax.grid(axis='x',alpha=.12)
fig.text(.155,.496,'The homogeneous model misses the broad path-length IQRs. Some net-displacement checks also miss; fitting the medians is not sufficient to reproduce the distributions.',fontsize=11,weight='normal',color='#8b3c2a')

# E width adequacy: main structural limitation.
ax=panel([.075,.336,.19,.117],'E','Population spread is not explained')
for i,id in enumerate(ids):
    oq=np.array(BYID[id]['path']);sq=simq(id,'path');ratio=(sq[2]-sq[0])/(oq[2]-oq[0])
    ax.bar(i,ratio,color=colors[id],width=.7)
ax.axhline(1,color='#555555',ls='--',lw=1)
ax.set(ylim=(0,1.08),ylabel='Simulated / reported path IQR',xticks=range(11),xticklabels=ids)
ax.tick_params(axis='x',labelsize=8)
note(ax,'One speed per condition produces too little cell-to-cell spread.\nHeterogeneity or temporal speed changes were not fitted.',y=-.24)

# F sensitivity: baseline persistence ratio depends on FMI and summary choice.
ax=panel([.335,.336,.18,.117],'F','How conditional is the persistence result?')
rows=[('Primary fit',fits['P']['tau']/fits['N']['tau'])]
for label,mode in [('Medians only','medians_only'),('FMI omitted','omit_FMI'),('FMI / net distance','FMI_normalized_by_net_distance')]:
    p=next(r for r in sens if r['sensitivity']==mode and r['id']=='P');n=next(r for r in sens if r['sensitivity']==mode and r['id']=='N')
    rows.append((label,float(p['tau'])/float(n['tau'])))
for i,(label,value) in enumerate(rows):
    ax.plot(value,i,'o',color='#333333' if i else C['PACC'],ms=7)
    ax.annotate(f'{value:.2g}×',(value,i),xytext=(7,0),textcoords='offset points',va='center',fontsize=10)
ax.axvline(1,color='#aaaaaa',ls='--');ax.set(xscale='log',xlim=(.8,150),ylim=(3.5,-.5),yticks=range(4),yticklabels=[r[0] for r in rows],xlabel='PACC / non-PACC intrinsic persistence')
ax.set_xticks([1,3,10,30,100]);ax.set_xticklabels(['1','3','10','30','100'])
note(ax,'Omitting FMI drives the non-PACC fit to the τ search floor.\nThe alternative FMI normalization is a sensitivity, not a definition.',y=-.24)

# G directional-statistic check at exact sample size.
ax=panel([.605,.336,.15,.117],'G','Directional-test check')
for i,id in enumerate(ids):
    r=next(r for r in ray if r['id']==id)
    lo,hi=float(r['predicted_p05']),float(r['predicted_p95'])
    ax.plot([lo,hi],[i,i],color=colors[id],lw=2)
    ax.plot(float(r['predicted_p50']),i,'o',color=colors[id],ms=4)
    ax.plot(float(r['reported_p']),i,'D',mfc='white',mec='black',ms=4)
ax.axvline(.05,color='#999999',ls='--',lw=1)
ax.set(xscale='log',xlim=(1e-9,1.5),ylim=(10.7,-.7),yticks=range(11),yticklabels=ids,xlabel='Rayleigh p, 60 tracks')
ax.set_xticks([1e-8,1e-4,.05,1]);ax.set_xticklabels(['10⁻⁸','10⁻⁴','.05','1'])
note(ax,'Bars: simulated 5–95% range; 3,000 sets of 60 independent cells.\nNo replicate correlation modeled. PACC directionality is too strong.',y=-.24)

# H specific perturbation prediction following the supplied threshold rule.
ax=panel([.825,.336,.135,.117],'H','Positive-FMI, high-motility PACCs')
for k,block in enumerate(('HIF','RhoA')):
    pre='HP' if block=='HIF' else 'RP';vals=[]
    for suffix in ('0','1'):
        r=next(r for r in classes if r['id']==pre+suffix and r['class_name']=='high' and r['direction']=='positive');vals.append(float(r['fraction']))
    ax.plot([k-.17,k+.17],vals,'o-',color=C['PACC'],lw=2,ms=6)
    for x,val in zip([k-.17,k+.17],vals):ax.annotate(f'{100*val:.1f}%',(x,val),xytext=(0,8),textcoords='offset points',ha='center',fontsize=9)
ax.set(xticks=[0,1],xticklabels=['HIF-1α','RhoA'],ylabel='Fraction of simulated PACCs',ylim=(0,.34),xlim=(-.5,1.5))
ax.text(.98,.97,'Control → KD',transform=ax.transAxes,ha='right',va='top',fontsize=9)
note(ax,'High motility: above simulated control-PACC endpoint Q75.\nObserved direction of change: decrease; no numerical target.',y=-.24)

# I two trajectory views with shared tracks and physical scales across groups.
fig.text(.057,.278,'I   Baseline realizations',fontsize=15,weight='bold')
for col,(duration,rect,limit) in enumerate([(1,[.075,.167,.135,.085],40),(16,[.24,.167,.135,.085],210)]):
    ax=fig.add_axes(rect)
    for id in ('O','N','P'):
        p=tracks[id][:12,:int(duration/.5)+1]
        for track in p:ax.plot(track[:,0],track[:,1],color=colors[id],alpha=.65,lw=.8)
    ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)')
    ax.set_title(f'{duration} h',loc='left',fontsize=11)
    ax.set_xticks([-limit,0,limit]);ax.set_yticks([-limit,0,limit]);ax.tick_params(labelsize=8)
    ax.annotate('O₂',xy=(.94,.92),xytext=(.94,.55),xycoords='axes fraction',textcoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#333333'),ha='center',fontsize=9)
fig.text(.075,.141,'12 cells per group; 1 h views are prefixes of the same 16 h tracks.\nLocal oxygen-rich direction is +y. Normoxic cells have no steering.',fontsize=9,color='#555555')

# J model-generated curves: explicitly no supplied MSD target.
ax=panel([.435,.167,.205,.085],'J','Predicted displacement across time')
for id in ('O','N','P'):
    rr=[r for r in times if r['id']==id]
    t=[float(r['time_h']) for r in rr];y=[float(r['msd']) for r in rr]
    ax.plot(t,y,color=colors[id],lw=2,label=labels[id])
ax.set(xscale='log',yscale='log',xlabel='Time (h)',ylabel='MSD from starting point (µm²)')
ax.set_xticks([.5,1,4,16]);ax.set_xticklabels(['.5','1','4','16'])
note(ax,'Prediction only: the notes contain no MSD or turning-angle curves.',y=-.33)

# K counterfactual remove directional guidance.
ax=panel([.72,.167,.24,.085],'K','Prediction: remove oxygen-directed steering')
for j,id in enumerate(('N','P')):
    before=simq(id,'directness')[1];after=float(next(r['median_directness'] for r in cf if r['id']==id))
    ax.plot([j-.13,j+.13],[before,after],'o-',color=colors[id],lw=2,ms=6)
    ax.annotate(f'{before:.2f} → {after:.2f}',(j,(before+after)/2),xytext=(0,15),textcoords='offset points',ha='center',fontsize=10)
ax.set(xticks=[0,1],xticklabels=['non-PACC','PACC'],xlim=(-.6,1.6),ylim=(0,.7),ylabel='Median directness')
note(ax,'Set κ = 0; retain fitted propulsion and intrinsic persistence.\nThe model predicts that much of the directness difference remains.',y=-.33)

# L model cannot describe radial and temporal heterogeneity with constant coefficients.
ax=panel([.075,.036,.285,.062],'L','Spatial observation outside this optimization')
ax.errorbar(RADIUS/1000,RADIAL_FMI,yerr=RADIAL_SE,fmt='o-',color='#555555',ms=4,lw=1,capsize=3)
ax.set(xlabel='Radius (mm)',ylabel='Mean parallel FMI',ylim=(-.03,.48))
ax.text(.98,.97,'Provided radial-bin means ± SEM',transform=ax.transAxes,ha='right',va='top',fontsize=9)

# Compact self-contained scientific interpretation and assumptions.
ax=fig.add_axes([.435,.026,.525,.089]);ax.axis('off')
ratio=fits['P']['tau']/fits['N']['tau'];vr=fits['P']['v']/fits['N']['v']
text=(f'M   Interpretation and limits\n'
      f'Baseline PACC/non-PACC: propulsion {vr:.2f}×; intrinsic persistence {ratio:.2f}×. Shared steering fits almost as well.\n'
      'Knockdown optima reduce propulsion and persistence; precise steering changes depend on which summaries are fitted.\n'
      'The model does not reproduce population spread or every prediction check. These are effective parameters, not unique mechanisms.\n\n'
      'Model: propulsion v; random-turning diffusion 1/τ; restoring angular drift −κ sin(heading) toward +y.\n'
      'Assumptions: constant local gradient; stationary initial headings; no speed heterogeneity or temporal switching.\n'
      'FMI is assumed to be projected net displacement / accumulated distance; the notes leave its normalization unresolved.\n'
      'No radial O₂ response was fitted; the supplied spatial profile and qualitative rise-and-fall in per-step FMI need a richer model.\n'
      'Optimization: nine quantiles per condition, equal metric weight after IQR scaling; errors are descriptive, not likelihoods.\n'
      'Validation: 8,192 new cells/condition; physical identities passed; baseline timestep and initial-heading sensitivities retained.\n'
      'Source: literature_review/hosny_2026.md. PACC labels are retrospective growth/size classes, not a clean WGD manipulation.')
ax.text(0,1,text,va='top',fontsize=9.5,linespacing=1.35,color='#333333')
fig.savefig(OUT/'hosny_optimization.png',dpi=190)
print(OUT/'hosny_optimization.png')
