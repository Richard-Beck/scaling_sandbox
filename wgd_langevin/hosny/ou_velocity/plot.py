from pathlib import Path
import sys,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,Patch
from matplotlib.ticker import NullFormatter
ROOT=Path(__file__).resolve().parent;sys.path.append(str(ROOT.parent))
from data import DATA,BYID
OUT=ROOT/'results'
def load(name):return list(csv.DictReader((OUT/name).open()))
fits={r['group']:{k:(float(v) if k in ('tau','sigma','random_rms','u','score') else v) for k,v in r.items()} for r in load('fits.csv')}
pred=load('predictions.csv');comp=load('model_comparison.csv');obs=load('source_observations.csv');counter=load('drift_removed.csv');times=load('time_predictions.csv');tracks=np.load(OUT/'tracks.npz')
G={'O':['O'],'N':['N','HN0','RN0'],'P':['P','HP0','RP0'],'HN':['HN1'],'HP':['HP1'],'RN':['RN1'],'RP':['RP1']}
C={'O':'#797979','N':'#286FA7','P':'#C5522E','HN':'#286FA7','HP':'#C5522E','RN':'#286FA7','RP':'#C5522E'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':12,'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(23,21),facecolor='white')
fig.text(.055,.975,'Hosny 2026 revisited: fluctuating OU velocities with shared baseline controls',fontsize=22,weight='bold',va='top')
fig.text(.055,.951,'Speed and direction fluctuate within each cell. Seven parameter sets are fitted jointly against 11 experimental summaries; all cells within a condition follow the same dynamics.',fontsize=11,color='#555555')
fig.legend(handles=[Line2D([],[],color=C[g],lw=3,label=l) for g,l in [('O','Normoxia'),('N','Hypoxic non-PACC'),('P','Hypoxic PACC')]],ncol=3,loc='upper left',bbox_to_anchor=(.055,.94),frameon=False,fontsize=10)

def panel(rect,letter,title):
    ax=fig.add_axes(rect);ax.set_title(title,loc='left',pad=11);ax.text(-.12,1.055,letter,transform=ax.transAxes,weight='bold',fontsize=16,va='bottom');return ax

def note(ax,text,y=-.24):ax.text(0,y,text,transform=ax.transAxes,fontsize=9,color='#555555',va='top')

def q(g,k):return np.array([float(r['predicted']) for r in sorted([r for r in pred if r['group']==g and r['metric']==k],key=lambda r:float(r['quantile']))])

def treatment_marker(id):return 'D' if len(id)==1 else ('o' if id.startswith('H') else 's')

def source_marker(ax,x,y,id,color='black',size=4):
    filled=id.endswith('1')
    ax.plot(x,y,marker=treatment_marker(id),mfc=color if filled else 'white',mec=color,ms=size,ls='none',zorder=6)

ax=panel([.075,.778,.245,.129],'A','Three baselines, with knockdown changes')
for g in ('O','N','P'):
    f=fits[g];ax.scatter(f['tau'],f['random_rms'],s=90,color=C[g],zorder=4)
    ax.annotate({'O':'Normoxia','N':'non-PACC','P':'PACC'}[g],(f['tau'],f['random_rms']),xytext=(7,7),textcoords='offset points',color=C[g],fontsize=10)
for base,kd,label,style in [('N','HN','HIF','-'),('P','HP','HIF','-'),('N','RN','RhoA','--'),('P','RP','RhoA','--')]:
    a,b=fits[base],fits[kd]
    ax.annotate('',xy=(b['tau'],b['random_rms']),xytext=(a['tau'],a['random_rms']),arrowprops=dict(arrowstyle='->',color=C[base],ls=style,lw=1.5))
    ax.scatter(b['tau'],b['random_rms'],marker='v' if label=='HIF' else 's',s=48,color=C[base])
    offsets={'HN':(-20,-17),'HP':(-12,12),'RN':(8,-14),'RP':(-35,14)}
    ax.annotate(label,(b['tau'],b['random_rms']),xytext=offsets[kd],textcoords='offset points',color=C[base],fontsize=9)
ax.set(xscale='log',xlim=(.11,3),ylim=(5.5,28),xlabel='Velocity-memory time (h)',ylabel='Random-motion RMS speed (µm/h)')
ax.set_xticks([.15,.3,.5,1,2]);ax.set_xticklabels(['.15','.3','.5','1','2']);ax.xaxis.set_minor_formatter(NullFormatter())
note(ax,'Hypoxic baselines combine untreated + both siNTC datasets.\nKnockdown arrows start from these shared baseline fits.')

ax=panel([.39,.778,.20,.129],'B','Oxygen-directed velocity is fitted jointly')
for j,g in enumerate(('O','N','P')):
    ax.plot(j,fits[g]['u'],'o',color=C[g],ms=8)
    if g!='O':
        for offset,kd,mk in [(-.12,'H'+g,'v'),(.12,'R'+g,'s')]:
            ax.plot([j,j+offset],[fits[g]['u'],fits[kd]['u']],color=C[g],lw=1,alpha=.6)
            ax.plot(j+offset,fits[kd]['u'],mk,color=C[g],ms=6)
ax.set(xticks=[0,1,2],xticklabels=['Normoxia','non-PACC','PACC'],xlim=(-.45,2.45),ylim=(-.15,3.2),ylabel='Mean velocity toward oxygen (µm/h)')
ax.legend(handles=[Line2D([],[],color='#444444',marker=m,ls='none',label=l) for m,l in [('o','Baseline'),('v','HIF KD'),('s','RhoA KD')]],loc='upper left',frameon=False,fontsize=9)
note(ax,'OU bias is a mean velocity, not the old angular steering rate.\nAmplitude, memory and bias were optimized simultaneously.')

ax=panel([.713,.778,.247,.129],'C','Independent refits of constrained models')
names=['Baseline PACC/non-PACC','HIF · non-PACC','HIF · PACC','RhoA · non-PACC','RhoA · PACC'];models=['amplitude_only','amplitude_and_bias','amplitude_and_memory','all_three']
mat=np.array([[float(next(r['score'] for r in comp if r['pair']==p and r['model']==m)) for m in models] for p in names])
ax.imshow(mat-mat[:,-1,None],cmap='YlOrRd',aspect='auto',vmin=0,vmax=max(.025,np.max(mat-mat[:,-1,None])))
for i in range(5):
    for j in range(4):ax.text(j,i,f'{mat[i,j]:.3f}',ha='center',va='center',fontsize=10)
ax.set(xticks=range(4),yticks=range(5));ax.set_xticklabels(['Amplitude\nonly','Amplitude\n+ bias','Amplitude\n+ memory','All three'],fontsize=9);ax.set_yticklabels(['Baseline P/N','HIF · non-PACC','HIF · PACC','RhoA · non-PACC','RhoA · PACC'],fontsize=9)
ax.tick_params(length=0);ax.set_frame_on(False)
note(ax,'Columns name parameters allowed to differ. Others are shared, not zeroed.\nNumbers: quantile error (lower is better); colour: excess over full fit.\nDescriptive scores, not likelihoods or significance tests.',y=-.30)

# D: one predicted distribution per biological/treatment state, separate source symbols.
fig.text(.055,.7,'D   One simulated baseline per population; each experiment remains a separate observation',fontsize=15,weight='bold')
fig.text(.055,.681,'Coloured bands: simulated Q25–Q75, with a median tick. Black symbols and bars: experimental medians and IQRs. IQRs describe population spread, not uncertainty.',fontsize=10.5,color='#555555')
fig.legend(handles=[Line2D([],[],color='black',marker='D',mfc='white',ls='none',label='Untreated'),
                    Line2D([],[],color='black',marker='o',mfc='white',ls='none',label='HIF siNTC'),
                    Line2D([],[],color='black',marker='s',mfc='white',ls='none',label='RhoA siNTC'),
                    Line2D([],[],color='black',marker='o',mfc='black',ls='none',label='HIF knockdown'),
                    Line2D([],[],color='black',marker='s',mfc='black',ls='none',label='RhoA knockdown')],
           ncol=5,loc='upper left',bbox_to_anchor=(.145,.675),frameon=False,fontsize=9,handletextpad=.45,columnspacing=1.9)
ypos={'O':0,'N':1.3,'HN':2.15,'RN':2.95,'P':4.35,'HP':5.2,'RP':6.0}
rowlabels={'O':'Normoxia','N':'non-PACC · baseline','HN':'non-PACC · HIF KD','RN':'non-PACC · RhoA KD','P':'PACC · baseline','HP':'PACC · HIF KD','RP':'PACC · RhoA KD'}
metric_info=[('path','Accumulated distance','µm',(0,315)),('directness','Track directness','net displacement / path',(0,.85)),('parallel','Along oxygen gradient','forward displacement / path',(-.4,.62)),('endpoint','Net displacement','µm',(0,145)),('perpendicular','Across oxygen gradient','sideways displacement / path',(-.45,.45))]
for j,(metric,title,xlabel,limits) in enumerate(metric_info):
    ax=fig.add_axes([.155+j*.162,.449,.135,.185])
    ax.set_title(title+'\n'+('USED IN FIT' if j<3 else 'PREDICTION CHECK'),loc='left',fontsize=10.8,pad=10)
    for g,y in ypos.items():
        sq=q(g,metric);height=.50 if g in ('N','P') else .34
        ax.add_patch(Rectangle((sq[0],y-height/2),sq[2]-sq[0],height,fc=C[g],alpha=.20,ec='none'))
        ax.plot([sq[1],sq[1]],[y-height/2,y+height/2],color=C[g],lw=2)
        offsets=[-.17,0,.17] if len(G[g])==3 else [0]
        for id,dy in zip(G[g],offsets):
            oq=np.array(BYID[id][metric]);yy=y+dy
            if np.isfinite(oq[[0,2]]).all():ax.plot([oq[0],oq[2]],[yy,yy],color='#555555',lw=.8,alpha=.8)
            source_marker(ax,oq[1],yy,id,size=4)
    ax.set(ylim=(6.6,-.65),xlim=limits,xlabel=xlabel,yticks=list(ypos.values()))
    ax.set_yticklabels(list(rowlabels.values()) if j==0 else ['']*len(ypos),fontsize=9)
    for sep in (.62,3.65):ax.axhline(sep,color='#dddddd',lw=.7)
    if metric in ('parallel','perpendicular'):ax.axvline(0,color='#cccccc',ls=':',lw=.8)
    ax.grid(axis='x',alpha=.10)
fig.text(.155,.408,'Directional fractions use the assumed forward-migration-index normalization: displacement projection / total path. The notes do not resolve that denominator.',fontsize=10,color='#555555')

# E demonstrate literal variable-speed Langevin dynamics.
ax=panel([.075,.252,.245,.115],'E','Speed now fluctuates within a cell')
for g in ('O','N','P'):
    vv=tracks[g+'_velocity'];t=np.arange(vv.shape[1])*.5
    for i in range(2):ax.plot(t,np.linalg.norm(vv[i],axis=1),color=C[g],lw=1,alpha=.75)
ax.set(xlabel='Time (h)',ylabel='Instantaneous speed (µm/h)',xlim=(0,16),ylim=(0,None))
note(ax,'Two random cell realizations per population.\nParameters are shared; the realized velocities are not constant.',y=-.29)

# F quantified distributional adequacy, no conflation with SEM.
ax=panel([.395,.252,.235,.115],'F','Does OU reproduce the path-length spread?')
for g in G:
    sq=q(g,'path');swidth=sq[2]-sq[0]
    for id in G[g]:
        oq=np.array(BYID[id]['path']);width=oq[2]-oq[0]
        source_marker(ax,width,swidth,id,color=C[g],size=6)
ax.plot([0,125],[0,125],'--',color='#999999',lw=1)
ax.set(xlim=(0,125),ylim=(0,125),xlabel='Reported accumulated-distance IQR (µm)',ylabel='Predicted IQR (µm)')
ax.text(.04,.96,'Dashed line: exact agreement',transform=ax.transAxes,va='top',fontsize=9)
note(ax,'Speed variability is real, but the OU model still underpredicts\nmost reported spreads. Quartiles were included in the loss.',y=-.29)

# G trajectories, matched prefixes, bias in +y.
fig.text(.697,.391,'G   Realizations from the shared baseline fits',fontsize=13,weight='bold')
for duration,rect,limit in [(1,[.70,.254,.105,.113],35),(16,[.855,.254,.105,.113],210)]:
    ax=fig.add_axes(rect)
    for g in ('O','N','P'):
        pp=tracks[g][:12,:int(duration/.5)+1]
        for p in pp:ax.plot(p[:,0],p[:,1],color=C[g],lw=.7,alpha=.65)
    ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)')
    ax.set_xticks([-limit,0,limit]);ax.set_yticks([-limit,0,limit]);ax.set_title(f'{duration} h',loc='left',fontsize=11);ax.tick_params(labelsize=8)
    ax.annotate('O₂',xy=(.92,.93),xytext=(.92,.56),xycoords='axes fraction',textcoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#333333'),ha='center',fontsize=9)
fig.text(.70,.215,'12 cells/group. Short tracks are the first hour of the long tracks.\nOxygen-rich direction is +y; normoxia has zero directional bias.',fontsize=9,color='#555555')

ax=panel([.075,.082,.245,.095],'H','Prediction: displacement across time')
for g in ('O','N','P'):
    rr=[r for r in times if r['group']==g];tt=np.array([float(r['time_h']) for r in rr]);yy=np.array([float(r['msd']) for r in rr]);mask=tt<=16
    ax.plot(tt[mask],yy[mask],color=C[g],lw=2);ax.plot(tt[~mask],yy[~mask],color=C[g],lw=1.5,ls='--')
ax.axvline(16,color='#cccccc',ls=':',lw=1)
ax.set(xscale='log',yscale='log',xlabel='Time (h)',ylabel='Mean squared displacement (µm²)')
ax.set_xticks([.5,1,4,16,32]);ax.set_xticklabels(['.5','1','4','16','32'])
note(ax,'No MSD data were supplied. Dashed curves extend beyond 16 h\nand are extrapolations of the stationary fitted model.',y=-.36)

ax=panel([.395,.082,.235,.095],'I','Prediction: remove the directional cue')
for j,g in enumerate(('N','P')):
    before=q(g,'directness')[1];after=float(next(r['drift_removed_directness'] for r in counter if r['group']==g))
    ax.plot([j-.12,j+.12],[before,after],'o-',color=C[g],lw=2,ms=6)
    ax.annotate(f'{before:.2f} → {after:.2f}',(j,(before+after)/2),xytext=(0,13),textcoords='offset points',ha='center',fontsize=10)
ax.set(xticks=[0,1],xticklabels=['non-PACC','PACC'],xlim=(-.5,1.5),ylim=(0,.7),ylabel='Median track directness')
note(ax,'Set directed velocity to zero, retain amplitude and memory.\nA residual directness difference is predicted.',y=-.36)

ax=fig.add_axes([.697,.041,.27,.153]);ax.axis('off')
ax.text(0,1,'J   Fitted parameters · current Gaussian OU model',va='top',fontsize=12,weight='bold')
order=['O','N','P','HN','RN','HP','RP']
labels={'O':'Normoxia','N':'non-PACC baseline','P':'PACC baseline','HN':'non-PACC · HIF KD','RN':'non-PACC · RhoA KD','HP':'PACC · HIF KD','RP':'PACC · RhoA KD'}
rows=[]
for g in order:
    f=fits[g]
    bias='0 fixed' if g=='O' else ('<0.001' if f['u']<.001 else f"{f['u']:.3f}")
    rows.append([labels[g],f"{f['sigma']:.2f}",f"{f['tau']:.3f}",bias])
tab=ax.table(cellText=rows,colLabels=['Condition','sigma\nµm/h','tau\nh','u\nµm/h'],cellLoc='center',colLoc='center',colWidths=[.49,.17,.17,.17],bbox=[0,.26,1,.64])
tab.auto_set_font_size(False);tab.set_fontsize(9)
for (i,j),cell in tab.get_celld().items():
    cell.set_edgecolor('#dddddd');cell.set_linewidth(.4)
    if i==0:cell.set_facecolor('#eeeeee');cell.set_text_props(weight='bold')
    elif j==0:cell.set_text_props(ha='left',color=C[order[i-1]])
ax.text(0,.20,'sigma: stationary velocity SD per axis; tau: velocity memory.\nu: mean oxygen-directed velocity (not a noise-angle parameter).\nPoint estimates; controls share the indicated baselines.\nParameter uncertainty is not represented by these rounded values.',va='top',fontsize=9,color='#444444',linespacing=1.3)
fig.text(.055,.015,'Exact integrated OU; 16,384 independent cells per fitted state; six validation checks passed. Source: literature_review/hosny_2026.md. Spatial and temporal changes in oxygen response are not modeled.',fontsize=9,color='#555555')
fig.savefig(OUT/'hosny_ou_shared_baselines.png',dpi=190)
print(OUT/'hosny_ou_shared_baselines.png')
