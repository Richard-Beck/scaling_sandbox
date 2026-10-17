from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';PARENT=ROOT.parent
sys.path.insert(0,str(PARENT.parent))
from data import BYID

def read(name):return list(csv.DictReader((OUT/name).open()))
models=json.load((OUT/'fits.json').open());refs=json.load((OUT/'original_reference.json').open());pred=read('detailed_predictions.csv');pp=read('detailed_parameters.csv');tracks=np.load(OUT/'tracks.npz');times=read('time_predictions.csv');repeat=read('frontier_repeat.csv')
oldpred=list(csv.DictReader((PARENT/'results/predictions.csv').open()));oldpars=list(csv.DictReader((PARENT/'results/fits.csv').open()))
best={g:next(f for f in models if f['base']==g and f['complexity']==2 and f['best_at_complexity']) for g in ('N','P')}
C={'O':'#777777','N':'#286FA7','P':'#C5522E','HN':'#286FA7','HP':'#C5522E','RN':'#286FA7','RP':'#C5522E'}
G={'O':['O'],'N':['N','HN0','RN0'],'P':['P','HP0','RP0'],'HN':['HN1'],'HP':['HP1'],'RN':['RN1'],'RP':['RP1']}
labels={'O':'Normoxia','N':'non-PACC baseline','P':'PACC baseline','HN':'non-PACC · HIF KD','RN':'non-PACC · RhoA KD','HP':'PACC · HIF KD','RP':'PACC · RhoA KD'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,'axes.titlesize':13,'xtick.labelsize':10,'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(24,24),facecolor='white')
fig.text(.05,.979,'Which parameter changes explain HIF and RhoA knockdown most economically?',fontsize=24,weight='bold',va='top')
fig.text(.05,.956,'Poisson-kick velocity model · 225 joint subset combinations + a no-change reference per PACC class · every model refits its baseline and both KDs jointly',fontsize=12,color='#444444')
fig.text(.05,.939,'tau: velocity decay time; λ: kick rate; A: Rayleigh kick magnitude scale; κ: von Mises directionality. Selected KD parameters can increase or decrease; all others equal baseline.',fontsize=11.5,color='#444444')

# A: all models, independent-simulation frontier, and original signed hypothesis.
for i,base in enumerate(('N','P')):
    ax=fig.add_axes([.075+i*.465,.782,.385,.131]);ax.set_title(('A   ' if i==0 else '')+('non-PACC' if base=='N' else 'PACC'),loc='left',weight='bold',fontsize=15)
    rr=[r for r in models if r['base']==base];rng=np.random.default_rng(430)
    x=np.array([r['complexity'] for r in rr]);y=np.array([r['independent_score'] for r in rr]);jitter=rng.uniform(-.085,.085,len(rr));jitter[x==0]=0
    ax.scatter(x+jitter,y,s=13,color=C[base],alpha=.23,edgecolors='none',label='All tested models')
    front=sorted([r for r in rr if r['pareto']],key=lambda r:r['complexity']);ax.plot([r['complexity'] for r in front],[r['independent_score'] for r in front],'o-',color=C[base],lw=2,ms=6,label='Pareto frontier',zorder=5)
    ref=next(r for r in refs if r['base']==base);ax.plot(2,ref['independent_score'],'X',ms=10,color='#222222',label='Original HIF κ↓ / RhoA tau↓',zorder=7)
    win=best[base];ax.plot(2,win['independent_score'],'o',ms=11,mfc='white',mec=C[base],mew=2,zorder=6,label='Best one parameter per KD')
    for r in front:ax.annotate(str(r['complexity']),(r['complexity'],r['independent_score']),xytext=(4,6),textcoords='offset points',fontsize=9,color=C[base])
    ax.set(xlim=(-.3,8.35),xticks=[0,2,3,4,5,6,7,8],xlabel='KD-specific parameters, in addition to four baseline parameters',ylabel='Independent simulated quantile error',ylim=(0,max(y.max(),ref['independent_score'])*1.12))
    ax.grid(axis='y',alpha=.12)
    if i==0:ax.legend(frameon=False,fontsize=9,loc='upper right',ncol=1)
    ax.text(0,-.28,'A point is dominated if an equally or less complex model fits better. Lower error is better.',transform=ax.transAxes,fontsize=9.5,color='#555555')

# B: winning mechanisms at each budget, with true frontier membership distinguished.
fig.text(.05,.720,'B   Best parameter subset at each complexity · filled dots identify the Pareto frontier',fontsize=15,weight='bold')
ax=fig.add_axes([.05,.613,.91,.096]);ax.axis('off')
def sym(text):return text.replace('lambda','λ').replace('kappa','κ').replace('tau','τ').replace('+',' + ')
rows=[];winners=[]
for k in (0,2,3,4,5,6,7,8):
    rr=[next(r for r in models if r['base']==b and r['complexity']==k and r['best_at_complexity']) for b in ('N','P')];winners.append(rr)
    row=[str(k)]
    for r in rr:row += [sym(r['HIF']),sym(r['RhoA']),('● ' if r['pareto'] else '')+f"{r['independent_score']:.4f}"]
    rows.append(row)
table=ax.table(cellText=rows,colLabels=['Extra pars','non-PACC · HIF','non-PACC · RhoA','Error','PACC · HIF','PACC · RhoA','Error'],cellLoc='center',colWidths=[.08,.17,.17,.10,.19,.19,.10],bbox=[0,0,1,1]);table.auto_set_font_size(False);table.set_fontsize(10)
for (i,j),cell in table.get_celld().items():
    cell.set_edgecolor('#dddddd');cell.set_linewidth(.4)
    if i==0:cell.set_facecolor('#eeeeee');cell.set_text_props(weight='bold')
    elif i==2:cell.set_facecolor('#fff5dc')
    if i>0 and j in (1,2,3):cell.set_text_props(color=C['N'])
    if i>0 and j in (4,5,6):cell.set_text_props(color=C['P'])
fig.text(.05,.599,'Subsets need not be nested. The four best non-PACC two-parameter tau/A combinations are numerically near-tied. These are descriptive fits, not significance tests.',fontsize=10,color='#555555')

# C: equal-complexity detailed comparison; fresh 65536-cell simulations for both.
fig.text(.05,.573,'C   At the original two-parameter budget: compare the best subset directly with the original hypothesis',fontsize=15,weight='bold')
selection=';  '.join(('non-PACC' if b=='N' else 'PACC')+': HIF '+sym(best[b]['HIF'])+', RhoA '+sym(best[b]['RhoA']) for b in ('N','P'))
fig.text(.05,.556,selection+'. These detailed panels illustrate an equal-complexity comparison, not an overall selected optimum.',fontsize=11,color='#555555')
fig.legend(handles=[Line2D([],[],color='#888888',lw=7,alpha=.4,label='Best subset: IQR + median'),Line2D([],[],color='#444444',ls='--',lw=1.2,label='Original hypothesis: IQR + median'),*[Line2D([],[],color='black',marker=m,mfc=fc,ls='none',label=l) for m,fc,l in [('D','white','Untreated'),('o','white','HIF siNTC'),('s','white','RhoA siNTC'),('o','black','HIF KD'),('s','black','RhoA KD')]]],ncol=7,loc='upper left',bbox_to_anchor=(.10,.548),frameon=False,fontsize=9,handletextpad=.4,columnspacing=1)

def q(g,k,model='best_two'):
    rr=[r for r in (oldpred if g=='O' else pred) if r['group']==g and r['metric']==k and (g=='O' or r['model']==model)]
    return np.array([float(r['predicted']) for r in sorted(rr,key=lambda r:float(r['quantile']))])
def point(ax,x,y,id,color='black',size=4):ax.plot(x,y,marker='D' if len(id)==1 else ('o' if id[0]=='H' else 's'),mfc=color if id.endswith('1') else 'white',mec=color,ms=size,ls='none',zorder=6)
ypos={'O':0,'N':1.3,'HN':2.15,'RN':2.95,'P':4.35,'HP':5.2,'RP':6.0}
for j,(metric,title,xlabel,limits) in enumerate([('path','Accumulated distance','µm',(0,315)),('directness','Track directness','net displacement / path',(0,.9)),('parallel','Along oxygen gradient','forward displacement / path',(-.42,.65)),('endpoint','Net displacement','µm',(0,155)),('perpendicular','Across oxygen gradient','sideways displacement / path',(-.48,.48))]):
    ax=fig.add_axes([.155+j*.162,.330,.135,.185]);ax.set_title(title+'\n'+('USED IN FIT' if j<3 else 'HELD-OUT CHECK'),loc='left',fontsize=11,pad=9)
    for g,y in ypos.items():
        sq=q(g,metric);height=.5 if g in ('N','P') else .34
        ax.add_patch(Rectangle((sq[0],y-height/2),sq[2]-sq[0],height,fc=C[g],alpha=.22,ec='none'));ax.plot([sq[1]]*2,[y-height/2,y+height/2],color=C[g],lw=2)
        if g!='O':
            oldq=q(g,metric,'original');yy=y+.30 if g in ('N','P') else y+.23
            ax.plot(oldq[[0,2]],[yy]*2,ls='--',color=C[g],lw=1.1,alpha=.9);ax.plot(oldq[1],yy,'|',color=C[g],ms=6)
        for id,dy in zip(G[g],[-.17,0,.17] if len(G[g])==3 else [0]):
            oq=np.array(BYID[id][metric]);yy=y+dy
            if np.isfinite(oq[[0,2]]).all():ax.plot(oq[[0,2]],[yy]*2,color='#555555',lw=.8)
            point(ax,oq[1],yy,id)
    ax.set(ylim=(6.65,-.65),xlim=limits,xlabel=xlabel,yticks=list(ypos.values()));ax.set_yticklabels([labels[g] for g in ypos] if j==0 else ['']*7,fontsize=9.5)
    for sep in (.62,3.65):ax.axhline(sep,color='#dddddd',lw=.7)
    if metric in ('parallel','perpendicular'):ax.axvline(0,color='#cccccc',ls=':',lw=.8)
    ax.grid(axis='x',alpha=.1)
fig.text(.155,.303,'16 h tracks, sampled every 30 min. Loss fits Q25/median/Q75, scaled by reported IQR; each source dataset has equal weight. IQRs describe population spread.',fontsize=10,color='#555555')
fig.text(.155,.289,'The directional fraction is assumed to mean displacement projection / total path; the supplied notes leave its denominator ambiguous. Normoxia is the unchanged reference fit.',fontsize=10,color='#555555')

# D: explicit fitted parameters of the highlighted comparison, not hidden separate controls.
ax=fig.add_axes([.05,.117,.35,.125]);ax.axis('off');ax.text(0,1.07,'D   Fitted parameters · best two-parameter subsets',fontsize=13,weight='bold')
pars={p['group']:p for p in pp if p['model']=='best_two'};pars['O']=next(p for p in oldpars if p['group']=='O');order=['O','N','HN','RN','P','HP','RP'];cells=[]
for g in order:cells.append([labels[g]]+[f'{float(pars[g][p]):.3g}' for p in ('tau','rate','amplitude','kappa')])
tab=ax.table(cellText=cells,colLabels=['Condition','tau (h)','λ (/h)','A (µm/h)','κ'],cellLoc='center',colWidths=[.46,.135,.135,.135,.135],bbox=[0,.12,1,.87]);tab.auto_set_font_size(False);tab.set_fontsize(9.5)
for (i,j),cell in tab.get_celld().items():
    cell.set_edgecolor('#dddddd');cell.set_linewidth(.4)
    if i==0:cell.set_facecolor('#eeeeee');cell.set_text_props(weight='bold')
    elif j==0:cell.set_text_props(ha='left',color=C[order[i-1]])
    elif len(order[i-1])==2:
        g=order[i-1];mask=best[g[-1]]['hmask' if g[0]=='H' else 'rmask']
        if mask&(1<<(j-1)):cell.set_facecolor('#fff1c4')
ax.text(0,.02,'Yellow: KD-specific parameters. All other entries share baseline values.\nPoint estimates only: sub-frame memory and kick size trade off.',fontsize=9.5,color='#555555',va='top')

# E: tracks from the same-budget winners, not the original or full model.
fig.text(.445,.255,'E   Predictions: baseline tracks from the highlighted fits',fontsize=13,weight='bold')
for duration,rect in [(1,[.46,.125,.105,.109]),(16,[.605,.125,.105,.109])]:
    ax=fig.add_axes(rect);allp=[]
    for g in ('N','P'):
        pp0=tracks[g][:12,:int(duration/.5)+1];allp.append(pp0)
        for p in pp0:ax.plot(p[:,0],p[:,1],color=C[g],lw=.8,alpha=.7)
    limit=max(10,float(np.ceil(max(abs(np.concatenate(allp).ravel()))/10)*10));ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)');ax.set_title(f'{duration} h',loc='left',fontsize=11);ax.set_xticks([-limit,0,limit]);ax.set_yticks([-limit,0,limit]);ax.tick_params(labelsize=8)
    ax.annotate('O₂',xy=(.90,.91),xytext=(.90,.6),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'),ha='center',fontsize=9)
fig.text(.46,.085,'12 homogeneous cells per hypoxic baseline. The 1 h tracks\nare prefixes of the same 16 h realizations.',fontsize=9.5,color='#555555')

ax=fig.add_axes([.795,.125,.165,.109]);ax.set_title('F   Predicted MSD',loc='left',fontsize=13,weight='bold')
for g in ('N','P'):
    rr=[r for r in times if r['group']==g];t=np.array([float(r['time_h']) for r in rr]);y=np.array([float(r['msd']) for r in rr]);mask=t<=16
    ax.plot(t[mask],y[mask],color=C[g],lw=2);ax.plot(t[~mask],y[~mask],color=C[g],ls='--',lw=1.5)
ax.axvline(16,color='#cccccc',ls=':');ax.set(xscale='log',yscale='log',xlabel='Lag / elapsed time (h)',ylabel='MSD (µm²)')
fig.text(.795,.085,'No supplied MSD targets. Dashed curves\nextend beyond the 16 h observation.',fontsize=9.5,color='#555555')

# Compact explicit limits and additional checks.
maxrepeat=max(abs(float(r['first_score'])-float(r['repeat_score'])) for r in repeat)
sub=read('qualitative_subpopulation_check.csv');sf={r['group']:float(r['positive_high_fraction']) for r in sub}
fig.text(.05,.055,f'Validation: 32,768 fresh cells/state for all 452 models; repeated frontier checks differed by at most {maxrepeat:.4f} in error. Detailed comparisons use 65,536 cells/state.',fontsize=10,color='#555555')
fig.text(.05,.039,f"Qualitative PACC check: positive/high-endpoint fraction {sf['P']:.0%} baseline → {sf['HP']:.0%} HIF, {sf['RP']:.0%} RhoA. Notes report reductions after both KDs; no numerical fractions were supplied.",fontsize=10,color='#555555')
fig.text(.05,.023,'Only supplied hosny_2026.md data. No likelihood or parameter CIs; a numerical frontier does not establish unique mechanisms. Spatial oxygen response, temporal adaptation and cell area are outside the model.',fontsize=9.5,color='#555555')
fig.savefig(OUT/'hosny_poisson_pareto.png',dpi=180)
print(OUT/'hosny_poisson_pareto.png')
