from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import sys
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';SHARED=ROOT.parent
sys.path.insert(0,str(SHARED.parent.parent))
from data import BYID

def read(path):return list(csv.DictReader(path.open()))
models=read(OUT/'selected_models.csv');pred=read(OUT/'predictions.csv');pars=read(OUT/'parameters.csv');times=read(OUT/'time_predictions.csv');tracks=np.load(OUT/'tracks.npz')
allm=read(SHARED/'results/all_shared_models.csv');front=read(SHARED/'results/shared_frontier.csv');ratios=read(SHARED/'results/frontier_ratios.csv');profiles=read(SHARED/'results/ratio_profiles.csv')
ORDER=['null','direction','memory','full'];COLOR={'null':'#858585','direction':'#d2a021','memory':'#008e98','full':'#8558ac'};MARK={'null':'s','direction':'o','memory':'^','full':'D'}
SHORT={'null':'K=0','direction':'K=2: κ / κ','memory':'K=2: κ / tau','full':'K=8'}
LABEL={'null':'K=0 · no KD changes','direction':'K=2 · HIF directionality / RhoA directionality','memory':'K=2 · HIF directionality / RhoA memory','full':'K=8 · all four parameters change for each KD'}
NC='#286FA7';PC='#C5522E';PURPLE='#6b4ca5';PARAMS=[('tau','Velocity memory, tau','h'),('rate','Kick rate, λ','events/h'),('amplitude','Kick magnitude scale, A','µm/h'),('kappa','Kick directionality, κ','dimensionless')]
P={(r['model'],r['group']):r for r in pars}
GROUPS={'O':['O'],'N':['N','HN0','RN0'],'HN':['HN1'],'RN':['RN1'],'P':['P','HP0','RP0'],'HP':['HP1'],'RP':['RP1']}
STATE=['O','N','HN','RN','P','HP','RP']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':10,'axes.titlesize':12,'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(26,37),facecolor='white')
fig.text(.05,.982,'Shared knockdown mechanisms: four representative explanations of PACC migration',fontsize=25,weight='bold',va='top')
fig.text(.05,.967,'Identical HIF/RhoA parameter subsets across PACC and non-PACC; parameter values remain population-specific. Every fit jointly uses baseline, control and KD observations.',fontsize=12,color='#444444')

# A: same frontier, four highlighted models; the legend replaces the old mechanism table.
ax=fig.add_axes([.075,.843,.40,.102]);ax.set_title('A   Shared-mechanism Pareto frontier · four choices',loc='left',weight='bold',fontsize=15)
x=np.array([int(r['complexity_per_population']) for r in allm]);y=np.array([float(r['score']) for r in allm]);rng=np.random.default_rng(13);dx=rng.uniform(-.08,.08,len(x));dx[x==0]=0
ax.scatter(x+dx,y,s=13,color='#aaaaaa',alpha=.35)
ff=sorted(front,key=lambda r:int(r['complexity_per_population']));ax.plot([int(r['complexity_per_population']) for r in ff],[float(r['score']) for r in ff],color='#777777',lw=1.3)
for r in models:
    m=r['model'];ax.plot(int(r['K']),float(r['score']),marker=MARK[m],color=COLOR[m],ms=9,ls='none',zorder=6)
ax.set(xticks=[0,2,3,4,5,6,7,8],xlabel='KD-specific parameters per population, K',ylabel='Combined fitting error',xlim=(-.3,8.3),ylim=(0,.067));ax.grid(axis='y',alpha=.12)
fig.legend(handles=[Line2D([],[],marker=MARK[m],color=COLOR[m],lw=3,ms=8,label=LABEL[m]) for m in ORDER],loc='upper left',bbox_to_anchor=(.535,.95),frameon=False,fontsize=12,labelspacing=1.1)
fig.text(.55,.860,'Model colours identify parameter-facet rows and data comparisons.\nBoth K=2 models allow the selected parameters to increase or decrease.\nK counts extra parameters per population: total fitted scalars = 8 + 2K.\nThe grey line is the frontier; no likelihood or statistical model ranking is implied.',fontsize=10.5,color='#555555',va='top',linespacing=1.5)

# B: four representative rows, four parameter columns. Each facet has its own scale.
# Add five inches here; preserve the physical size of all other panels below.
B_AXES=[];B_TEXT=[];RATIO=32/37;INSERT=5/32

def btext(x,y,text,**kwargs):
    artist=fig.text(x,y*RATIO,text,**kwargs);B_TEXT.append(artist);return artist

btext(.05,.807+INSERT,'B   Parameter values: adjacent non-PACC and PACC bars within each representative',fontsize=16,weight='bold')
btext(.05,.793+INSERT,'Blue: non-PACC; orange: PACC. Each facet has its own linear scale starting at zero; compare bar heights within a facet and read values across facets.',fontsize=11,color='#555555')
norm=P['reference','O']
btext(.05,.781+INSERT,f"Common normoxia reference: tau = {float(norm['tau']):.3g} h; λ = {float(norm['rate']):.3g} /h; A = {float(norm['amplitude']):.3g} µm/h; κ = {float(norm['kappa']):.3g}.",fontsize=10.5,color='#777777')
for row,m in enumerate(ORDER):
    bottom=.711+INSERT-row*.064
    btext(.035,bottom+.023,SHORT[m],fontsize=12,color=COLOR[m],weight='bold',rotation=90,va='center',ha='center')
    for j,(key,title,unit) in enumerate(PARAMS):
        ax=fig.add_axes([.075+j*.231,bottom*RATIO,.184,.047*RATIO]);B_AXES.append(ax)
        width=.32;maximum=0.
        for pop,groups,color,offset in [('N',['N','HN','RN'],NC,-width/2),('P',['P','HP','RP'],PC,width/2)]:
            values=[float(P[m,g][key]) for g in groups];maximum=max(maximum,max(values))
            bars=ax.bar(np.arange(3)+offset,values,width=width,color=color,edgecolor='white',linewidth=.35)
            for bar,value in zip(bars,values):
                label='0' if value==0 else (f'{value:.1e}' if value<.001 else f'{value:.3g}')
                ax.annotate(label,(bar.get_x()+bar.get_width()/2,value),xytext=(0,3),textcoords='offset points',ha='center',va='bottom',fontsize=8.5,color=color)
        ax.set(xticks=[0,1,2],xticklabels=['Baseline','HIF KD','RhoA KD'],ylabel=unit,xlim=(-.6,2.6),ylim=(0,max(maximum*1.27,1e-8)))
        if row==0:ax.set_title(title,loc='left',fontsize=12,pad=8)
        ax.grid(axis='y',alpha=.10);ax.tick_params(axis='x',labelsize=10)
btext(.075,.653,'Bars are fitted point estimates without confidence intervals. A is a velocity-jump scale; short memory and large kicks can compensate below the sampling interval.',fontsize=10,color='#555555')

# C: every trajectory summary, all four models, actual source controls remain separate.
fig.text(.05,.626,'C   Fits and held-out checks: model quartiles versus every reported trajectory summary',fontsize=16,weight='bold')
fig.text(.05,.612,'Coloured bands/ticks: model Q25–Q75/median, top to bottom in the same four-model order. Black observations sit just below each model stack; their bars are population IQRs.',fontsize=10.5,color='#555555')
fig.legend(handles=[Line2D([],[],color='black',marker=m,mfc=fc,ls='none',label=l) for m,fc,l in [('D','white','Untreated'),('o','white','HIF siNTC'),('s','white','RhoA siNTC'),('o','black','HIF KD'),('s','black','RhoA KD')]],ncol=5,loc='upper left',bbox_to_anchor=(.17,.607),frameon=False,fontsize=10)
def q(m,g,k):
    rr=[r for r in pred if r['model']==m and r['group']==g and r['metric']==k]
    return np.array([float(r['predicted']) for r in sorted(rr,key=lambda r:float(r['quantile']))])
ypos={'O':0,'N':1.5,'HN':2.8,'RN':4.1,'P':5.7,'HP':7.,'RP':8.3}
labels={'O':'Normoxia','N':'non-PACC · baseline','HN':'non-PACC · HIF KD','RN':'non-PACC · RhoA KD','P':'PACC · baseline','HP':'PACC · HIF KD','RP':'PACC · RhoA KD'}
for j,(metric,title,xlabel,limits) in enumerate([('path','Accumulated distance','µm',(0,315)),('directness','Track directness','net displacement / path',(0,.9)),('parallel','Along oxygen gradient','forward displacement / path',(-.42,.65)),('endpoint','Net displacement','µm',(0,155)),('perpendicular','Across oxygen gradient','sideways displacement / path',(-.48,.48))]):
    ax=fig.add_axes([.155+j*.162,.444,.135,.138]);ax.set_title(title+'\n'+('USED IN FIT' if j<3 else 'HELD-OUT CHECK'),loc='left',fontsize=11,pad=8)
    for g,y in ypos.items():
        for i,m in enumerate(['reference'] if g=='O' else ORDER):
            yy=y if g=='O' else y+(i-1.5)*.18;sq=q(m,g,metric);color='#999999' if g=='O' else COLOR[m]
            ax.plot(sq[[0,2]],[yy]*2,color=color,lw=5,alpha=.35,solid_capstyle='butt');ax.plot([sq[1]]*2,[yy-.07,yy+.07],color=color,lw=1.6)
        offsets=[-.09,0,.09] if len(GROUPS[g])==3 else [0]
        for id,dy in zip(GROUPS[g],offsets):
            oq=np.array(BYID[id][metric]);yy=y+(.22 if g=='O' else .48)+dy
            if np.isfinite(oq[[0,2]]).all():ax.plot(oq[[0,2]],[yy]*2,color='#444444',lw=.8)
            ax.plot(oq[1],yy,marker='D' if len(id)==1 else ('o' if id[0]=='H' else 's'),mfc='black' if id.endswith('1') else 'white',mec='black',ls='none',ms=3.6,zorder=7)
    ax.set(ylim=(9.15,-.4),xlim=limits,xlabel=xlabel,yticks=list(ypos.values()));ax.set_yticklabels([labels[g] for g in ypos] if j==0 else ['']*7,fontsize=9.5)
    for sep in (.85,4.95):ax.axhline(sep,color='#dddddd',lw=.7)
    if metric in ('parallel','perpendicular'):ax.axvline(0,color='#dddddd',ls=':',lw=.8)
    ax.grid(axis='x',alpha=.1)
fig.text(.155,.425,'65,536 independent homogeneous cells/state; 16 h sampled every 30 min. Directional-index normalization (projection/path) remains an assumption unresolved in the notes.',fontsize=10,color='#555555')

# D: 4 models x 2 durations, with consistent spatial scales across choices.
fig.text(.05,.403,'D   Baseline realizations: 1 h and 16 h',fontsize=15,weight='bold')
fig.text(.05,.392,'Blue: non-PACC; orange: PACC. Twelve cells per population.',fontsize=10,color='#555555')
for row,duration in enumerate((1,16)):
    limit=max(10,float(np.ceil(max(np.abs(tracks[m+'_'+g][:12,:int(duration/.5)+1]).max() for m in ORDER for g in ('N','P'))/10)*10))
    for col,m in enumerate(ORDER):
        ax=fig.add_axes([.07+col*.111,.315-row*.079,.088,.061])
        for g,c in [('N',NC),('P',PC)]:
            for p in tracks[m+'_'+g][:12,:int(duration/.5)+1]:ax.plot(p[:,0],p[:,1],color=c,lw=.7,alpha=.65)
        ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit));ax.set_xticks([-limit,0,limit]);ax.set_yticks([-limit,0,limit]);ax.tick_params(labelsize=7)
        if row==0:ax.set_title(SHORT[m],fontsize=9.5,color=COLOR[m],loc='left',pad=5)
        if col==0:ax.set_ylabel(f'{duration} h · y (µm)',fontsize=9)
        if row==1:ax.set_xlabel('x (µm)',fontsize=8)
        ax.annotate('O₂',xy=(.90,.91),xytext=(.90,.60),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.7),ha='center',fontsize=7)
fig.text(.07,.221,'Short tracks are exact prefixes; spatial limits match across models at each duration.',fontsize=9.5,color='#555555')

# E: all fitted conditions' analytic MSD, faceted by selected model.
fig.text(.555,.403,'E   Predicted MSD across all fitted hypoxic conditions',fontsize=15,weight='bold')
fig.text(.555,.392,'Blue: non-PACC; orange: PACC. Baseline —; HIF – –; RhoA ···',fontsize=10,color='#555555')
msdvals=[float(r['msd']) for r in times];ylim=(10**np.floor(np.log10(min(msdvals))),10**np.ceil(np.log10(max(msdvals))))
for i,m in enumerate(ORDER):
    row,col=divmod(i,2);ax=fig.add_axes([.585+col*.214,.319-row*.079,.17,.056])
    for g in ('N','HN','RN','P','HP','RP'):
        rr=[r for r in times if r['model']==m and r['group']==g];t=np.array([float(r['time']) for r in rr]);y=np.array([float(r['msd']) for r in rr]);style='-' if len(g)==1 else ('--' if g[0]=='H' else ':')
        ax.plot(t,y,color=NC if g.endswith('N') else PC,ls=style,lw=1.2,alpha=.85)
    ax.axvspan(16,32,color='#dddddd',alpha=.35);ax.set(xscale='log',yscale='log',xlim=(.1,32),ylim=ylim);ax.set_title(SHORT[m],loc='left',fontsize=10,color=COLOR[m],pad=4)
    if col==0:ax.set_ylabel('MSD (µm²)',fontsize=9)
    if row==1:ax.set_xlabel('Lag / elapsed time (h)',fontsize=9)
    ax.tick_params(labelsize=8)
fig.text(.585,.221,'No supplied MSD targets. Shading marks extrapolation beyond 16 h.',fontsize=9.5,color='#555555')

# F/G: compact paired 2x2 facets, preserving both kinds of sensitivity.
fig.text(.05,.201,'F   Baseline PACC/non-PACC ratios across model choices',fontsize=14,weight='bold')
fig.text(.05,.190,'Line: shared frontier; symbols: the four highlighted choices.',fontsize=9.5,color='#555555')
fig.text(.555,.201,'G   Within-model ratio profiles · κ / κ model only',fontsize=14,weight='bold')
fig.text(.555,.190,'Fix a baseline ratio; jointly refit the remaining parameters.',fontsize=9.5,color='#555555')
for i,(key,title,unit) in enumerate(PARAMS):
    row,col=divmod(i,2)
    ax=fig.add_axes([.075+col*.222,.125-row*.084,.169,.050])
    rr=sorted([r for r in ratios if r['parameter']==key and int(r['complexity'])>0],key=lambda r:int(r['complexity']));xx=np.array([int(r['complexity']) for r in rr]);yy=np.array([float(r['ratio']) for r in rr]);ax.plot(xx,yy,'o-',color='#aaaaaa',lw=1,ms=3)
    vals=list(yy)
    for m in ORDER:
        r=float(P[m,'P'][key])/float(P[m,'N'][key]);vals.append(r);k=int(next(v['K'] for v in models if v['model']==m));ax.plot(k,r,marker=MARK[m],ls='none',color=COLOR[m],ms=6,zorder=5)
    ax.axhline(1,color='#888888',ls='--',lw=.8);ax.set(yscale='log',ylim=(min(.7,min(vals)*.65),max(1.4,max(vals)*1.5)),xlim=(-.3,8.3),xticks=[0,2,4,6,8]);ax.set_title(title,loc='left',fontsize=10,pad=4)
    ticks={'tau':[1,2,5,10,20],'rate':[.5,.7,1],'amplitude':[.1,.2,.5,1],'kappa':[1,1.5,2,3]}[key];ax.set_yticks(ticks);ax.set_yticklabels([f'{t:g}' for t in ticks]);ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());ax.tick_params(labelsize=8)
    if col==0:ax.set_ylabel('PACC / non-PACC',fontsize=9)
    if row==1:ax.set_xlabel('KD parameters per population, K',fontsize=9)
    ax=fig.add_axes([.585+col*.214,.125-row*.084,.17,.050]);rr=sorted([r for r in profiles if r['parameter']==key],key=lambda r:float(r['ratio']));x=np.array([float(r['ratio']) for r in rr]);y=np.array([float(r['excess']) for r in rr]);r=float(P['direction','P'][key])/float(P['direction','N'][key])
    ax.axvline(1,color='#888888',ls='--',lw=.8);ax.axhline(0,color='#bbbbbb',lw=.7);ax.plot(x,y,'o-',color=COLOR['direction'],lw=1,ms=2.5);ax.plot(r,0,'*',color=COLOR['direction'],ms=9)
    ax.set(xscale='log',xlim=(.025,36),ylim=(-.0002,.045));ax.set_yscale('symlog',linthresh=.002);ax.set_yticks([0,.001,.01,.04]);ax.set_yticklabels(['0','.001','.01','.04']);ax.set_xticks([.03,.1,1,10,30]);ax.set_xticklabels(['.03','.1','1','10','30']);ax.set_title(title,loc='left',fontsize=10,pad=4);ax.tick_params(labelsize=8)
    if col==0:ax.set_ylabel('Extra fitting error',fontsize=9)
    if row==1:ax.set_xlabel('Fixed baseline PACC / non-PACC',fontsize=9)
fig.text(.05,.015,'Supplied hosny_2026.md data only; no new fits in this presentation. Profiles are descriptive, not confidence intervals; all G facets share the same error scale. Spatial oxygen response and cell area are outside this stationary model.',fontsize=9.5,color='#555555')
# Insert canvas space for the expanded B grid without squeezing C–G.
for ax in fig.axes:
    if ax in B_AXES:continue
    pos=ax.get_position();shift=INSERT if pos.y0>=.81 else 0
    ax.set_position([pos.x0,(pos.y0+shift)*RATIO,pos.width,pos.height*RATIO])
for artist in fig.texts:
    if artist in B_TEXT:continue
    x,y=artist.get_position();artist.set_position((x,(y+(INSERT if y>=.81 else 0))*RATIO))
for legend in fig.legends:
    box=legend.get_bbox_to_anchor().transformed(fig.transFigure.inverted())
    legend.set_bbox_to_anchor((box.x0,(box.y0+(INSERT if box.y0>=.81 else 0))*RATIO),transform=fig.transFigure)
fig.savefig(OUT/'hosny_four_representatives.png',dpi=180);print(OUT/'hosny_four_representatives.png')
