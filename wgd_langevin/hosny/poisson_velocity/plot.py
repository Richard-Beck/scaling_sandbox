from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results';sys.path.insert(0,str(ROOT.parent))
from data import BYID

def read(name):return list(csv.DictReader((OUT/name).open()))
fits={r['group']:{k:(v if k in ('group','source_ids') else float(v)) for k,v in r.items()} for r in read('fits.csv')}
pred=read('predictions.csv');comp=read('independent_scores.csv');profiles=read('rate_profiles.csv');memory_profiles=read('memory_profiles.csv');times=read('time_predictions.csv');ray=read('rayleigh_checks.csv');sub=read('qualitative_subpopulation_check.csv');cue=read('cue_removed.csv');tracks=np.load(OUT/'tracks.npz')
G={'O':['O'],'N':['N','HN0','RN0'],'P':['P','HP0','RP0'],'HN':['HN1'],'HP':['HP1'],'RN':['RN1'],'RP':['RP1']}
C={'O':'#777777','N':'#286FA7','P':'#C5522E','HN':'#286FA7','HP':'#C5522E','RN':'#286FA7','RP':'#C5522E'}
labels={'O':'Normoxia','N':'non-PACC baseline','P':'PACC baseline','HN':'non-PACC · HIF KD','RN':'non-PACC · RhoA KD','HP':'PACC · HIF KD','RP':'PACC · RhoA KD'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':12,'xtick.labelsize':9,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig=plt.figure(figsize=(24,23),facecolor='white')
fig.text(.05,.98,'Hosny 2026: Poisson kicks, jointly fitted baseline and knockdown mechanisms',fontsize=23,weight='bold',va='top')
fig.text(.05,.956,'Between kicks: dv/dt = −v/tau. Poisson rate λ; Rayleigh kick magnitude scale A; von Mises angle concentration κ toward +y. Position follows velocity.',fontsize=12,color='#444444')
fig.text(.05,.938,'Primary hypothesis: HIF reduces κ only; RhoA reduces tau only. Baseline and KD observations jointly determine shared parameters. Homogeneous populations throughout.',fontsize=12,color='#444444')
fig.legend(handles=[Line2D([],[],color=C[g],lw=3,label=l) for g,l in [('O','Normoxia'),('N','Hypoxic non-PACC'),('P','Hypoxic PACC')]],ncol=3,loc='upper left',bbox_to_anchor=(.045,.925),frameon=False,fontsize=10)

def panel(rect,letter,title):
    ax=fig.add_axes(rect);ax.set_title(title,loc='left',pad=12);ax.text(-.10,1.045,letter,transform=ax.transAxes,weight='bold',fontsize=16,va='bottom');return ax

def note(ax,text,y=-.26):ax.text(0,y,text,transform=ax.transAxes,fontsize=9.5,color='#555555',va='top',linespacing=1.35)

def q(g,k):return np.array([float(r['predicted']) for r in sorted([r for r in pred if r['group']==g and r['metric']==k],key=lambda r:float(r['quantile']))])

def marker(id):return 'D' if len(id)==1 else ('o' if id.startswith('H') else 's')

def point(ax,x,y,id,color='black',size=4):ax.plot(x,y,marker=marker(id),mfc=color if id.endswith('1') else 'white',mec=color,ms=size,ls='none',zorder=6)

ax=fig.add_axes([.05,.768,.34,.13]);ax.axis('off');ax.text(0,1.05,'A   Fitted parameters · primary joint model',fontsize=14,weight='bold')
order=['O','N','HN','RN','P','HP','RP']
rows=[]
for g in order:
    f=fits[g];rows.append([labels[g],f"{f['tau']:.3g}",f"{f['rate']:.3g}",f"{f['amplitude']:.3g}",'0 fixed' if g=='O' else f"{f['kappa']:.3g}"])
tab=ax.table(cellText=rows,colLabels=['Condition','tau (h)','λ (/h)','A (µm/h)','κ'],cellLoc='center',colWidths=[.46,.135,.135,.135,.135],bbox=[0,.10,1,.88]);tab.auto_set_font_size(False);tab.set_fontsize(10)
for (i,j),cell in tab.get_celld().items():
    cell.set_edgecolor('#dddddd');cell.set_linewidth(.4)
    if i==0:cell.set_facecolor('#eeeeee');cell.set_text_props(weight='bold')
    elif j==0:cell.set_text_props(ha='left',color=C[order[i-1]])
    elif (order[i-1].startswith('H') and j==4) or (order[i-1].startswith('R') and j==1):cell.set_facecolor('#fff1c4')
ax.text(0,.01,'Highlighted: the only KD-specific parameters. γ = 1/tau.\nA is a velocity jump scale. Sub-frame tau estimates are weakly resolved.',va='top',fontsize=9.5,color='#555555')

ax=panel([.455,.788,.18,.105],'B','The jointly fitted interventions')
for base in ('N','P'):
    f=fits[base];h=fits['H'+base];r=fits['R'+base]
    for target,mk in [(h,'v'),(r,'s')]:
        ax.annotate('',xy=(target['tau'],target['kappa']),xytext=(f['tau'],f['kappa']),arrowprops=dict(arrowstyle='->',color=C[base],lw=1.7))
        ax.plot(target['tau'],target['kappa'],mk,color=C[base],ms=7)
    ax.plot(f['tau'],f['kappa'],'o',color=C[base],ms=8)
ax.set(xlabel='Velocity memory, tau (h)',ylabel='Kick directionality, κ',xlim=(0,None),ylim=(-.02,None))
ax.legend(handles=[Line2D([],[],color='#444444',marker=m,ls='none',label=l) for m,l in [('o','Baseline'),('v','HIF KD'),('s','RhoA KD')]],frameon=False,fontsize=9)
note(ax,'Rate and kick amplitude are shared with both KDs.\nArrows depict jointly fitted parameter differences.',y=-.29)

ax=panel([.75,.788,.205,.105],'C','Do extra KD parameters improve the fit?')
models=['primary','HIF_memory','RhoA_direction','both','full'];mat=np.array([[float(next(r['independent_score'] for r in comp if r['base']==b and r['model']==m)) for m in models] for b in ('N','P')])
ax.imshow(mat,cmap='YlOrRd',aspect='auto',vmin=0,vmax=mat.max()*1.25)
for i in range(2):
    for j in range(5):ax.text(j,i,f'{mat[i,j]:.3f}',ha='center',va='center',fontsize=10)
ax.set(xticks=range(5),yticks=[0,1],yticklabels=['non-PACC','PACC']);ax.set_xticklabels(['Primary\n6 pars','+ HIF\ntau · 7','+ RhoA\nκ · 7','Both\n8 pars','All free\n12 pars'],fontsize=9);ax.set_frame_on(False);ax.tick_params(length=0)
note(ax,'All models refitted jointly, including their baselines.\nIndependent-simulation quantile error; lower is better.\nNo complexity penalty or statistical significance implied.',y=-.34)

fig.text(.05,.730,'Main limitation: RhoA memory reduction retains oxygen bias. Allowing RhoA directionality to change cuts PACC error from 0.051 to 0.025.',fontsize=11.5,weight='bold',color='#555555')
fig.text(.05,.706,'D   Every reported trajectory summary, compared with one model distribution per condition',fontsize=16,weight='bold')
fig.text(.05,.687,'Coloured bands: independent simulated Q25–Q75 and median. Black symbols/bars: experimental medians/IQRs. IQRs show population spread, not uncertainty.',fontsize=11,color='#555555')
fig.legend(handles=[Line2D([],[],color='black',marker=m,mfc=fc,ls='none',label=l) for m,fc,l in [('D','white','Untreated'),('o','white','HIF siNTC'),('s','white','RhoA siNTC'),('o','black','HIF KD'),('s','black','RhoA KD')]],ncol=5,loc='upper left',bbox_to_anchor=(.145,.680),frameon=False,fontsize=10)
ypos={'O':0,'N':1.3,'HN':2.15,'RN':2.95,'P':4.35,'HP':5.2,'RP':6.0}
for j,(metric,title,xlabel,limits) in enumerate([('path','Accumulated distance','µm',(0,315)),('directness','Track directness','net displacement / path',(0,.9)),('parallel','Along oxygen gradient','forward displacement / path',(-.42,.65)),('endpoint','Net displacement','µm',(0,155)),('perpendicular','Across oxygen gradient','sideways displacement / path',(-.48,.48))]):
    ax=fig.add_axes([.155+j*.162,.455,.135,.192]);ax.set_title(title+'\n'+('USED IN FIT' if j<3 else 'HELD-OUT CHECK'),loc='left',fontsize=11,pad=9)
    for g,y in ypos.items():
        sq=q(g,metric);h=.5 if g in ('N','P') else .34
        ax.add_patch(Rectangle((sq[0],y-h/2),sq[2]-sq[0],h,fc=C[g],alpha=.22,ec='none'));ax.plot([sq[1]]*2,[y-h/2,y+h/2],color=C[g],lw=2)
        for id,dy in zip(G[g],[-.17,0,.17] if len(G[g])==3 else [0]):
            oq=np.array(BYID[id][metric]);yy=y+dy
            if np.isfinite(oq[[0,2]]).all():ax.plot(oq[[0,2]],[yy]*2,color='#555555',lw=.8)
            point(ax,oq[1],yy,id)
    ax.set(ylim=(6.6,-.65),xlim=limits,xlabel=xlabel,yticks=list(ypos.values()));ax.set_yticklabels([labels[g] for g in ypos] if j==0 else ['']*7,fontsize=10)
    for sep in (.62,3.65):ax.axhline(sep,color='#dddddd',lw=.7)
    if metric in ('parallel','perpendicular'):ax.axvline(0,color='#cccccc',ls=':',lw=.8)
    ax.grid(axis='x',alpha=.1)
fig.text(.155,.420,'Fit: 16 h tracks sampled every 30 min; Q25, median and Q75; equal weights per source dataset and metric, with residuals divided by the reported IQR.',fontsize=10,color='#555555')
fig.text(.155,.405,'Directional-fraction normalization is assumed: projection / total path. The supplied notes leave its denominator ambiguous; alternative-normalization fits are retained.',fontsize=10,color='#555555')

ax=panel([.07,.262,.25,.102],'E','Speeds fluctuate within homogeneous populations')
for g in ('O','N','P'):
    v=tracks[g+'_velocity'];t=np.arange(v.shape[1])*.5
    for i in range(2):ax.plot(t,np.linalg.norm(v[i],axis=1),color=C[g],alpha=.75,lw=1)
ax.set(xlabel='Time (h)',ylabel='Instantaneous speed (µm/h)',xlim=(0,16),ylim=(0,None))
note(ax,'Two cells per baseline, shown at the experimental frame times.\nKicks occur between frames too; straight lines join sampled speeds.',y=-.31)

ax=panel([.40,.262,.23,.102],'F','Can the model reproduce path-length variability?')
for g in G:
    sq=q(g,'path');width=sq[2]-sq[0]
    for id in G[g]:
        oq=np.array(BYID[id]['path']);point(ax,oq[2]-oq[0],width,id,C[g],6)
ax.plot([0,125],[0,125],'--',color='#999999',lw=1);ax.set(xlim=(0,125),ylim=(0,125),xlabel='Reported path-length IQR (µm)',ylabel='Simulated path-length IQR (µm)')
note(ax,'Dashed line: equal observed and predicted spread.\nQuartiles constrain this variability during optimization.',y=-.31)

fig.text(.70,.384,'G   Parameter profiles · jointly refit the rest',fontsize=13,weight='bold')
for idx,(rrall,param,xlabel) in enumerate([(profiles,'rate','Kick rate λ (/h)'),(memory_profiles,'tau','Memory tau (h)')]):
    ax=fig.add_axes([.725+idx*.127,.262,.100,.102])
    for g in ('O','N','P'):
        rr=[r for r in rrall if r['base']==g];x=np.array([float(r[param]) for r in rr]);y=np.array([float(r['score']) for r in rr]);ref=float(next(r['bank_score'] for r in comp if r['base']==g and r['model']=='primary'))
        ax.plot(x,y-ref,'o-',ms=2.5,lw=1.2,color=C[g]);ax.plot(fits[g][param],0,'*',ms=9,color=C[g])
    ax.set(xscale='log',yscale='symlog',xlabel=xlabel,ylim=(-.002,10));ax.set_yscale('symlog',linthresh=.005)
    ax.set_yticks([0,.01,.1,1,10]);ax.set_yticklabels(['0','.01','.1','1','10'] if idx==0 else ['']*5);ax.axhline(0,color='#bbbbbb',lw=.7)
    if idx==0:ax.set_ylabel('Excess joint quantile error')
fig.text(.725,.220,'Stars: fitted values. Profiles reoptimize shared and KD parameters.\nFlat regions indicate weak identification; these are not CIs.\nVertical scale expands errors near zero.',fontsize=9.5,color='#555555',va='top')

fig.text(.05,.184,'H   Predictions: the same tracks over 1 h and 16 h',fontsize=13,weight='bold')
for duration,rect in [(1,[.07,.071,.108,.095]),(16,[.225,.071,.108,.095])]:
    ax=fig.add_axes(rect);allp=[]
    for g in ('N','P'):
        pp=tracks[g][:12,:int(duration/.5)+1];allp.append(pp)
        for p in pp:ax.plot(p[:,0],p[:,1],color=C[g],lw=.8,alpha=.7)
    limit=float(np.ceil(max(abs(np.concatenate(allp).ravel()))/10)*10);limit=max(limit,10)
    ax.set(aspect='equal',xlim=(-limit,limit),ylim=(-limit,limit),xlabel='x (µm)',ylabel='y (µm)');ax.set_title(f'{duration} h',loc='left',fontsize=11);ax.set_xticks([-limit,0,limit]);ax.set_yticks([-limit,0,limit]);ax.tick_params(labelsize=8)
    ax.annotate('O₂',xy=(.91,.91),xytext=(.91,.6),xycoords='axes fraction',arrowprops=dict(arrowstyle='->'),ha='center',fontsize=9)
fig.text(.07,.027,'12 cells per hypoxic baseline; short tracks are prefixes.\nAll plots use the primary constrained fit.',fontsize=9.5,color='#555555')

ax=panel([.40,.074,.23,.094],'I','Prediction: mean squared displacement')
for g in ('O','N','P'):
    rr=[r for r in times if r['group']==g];t=np.array([float(r['time_h']) for r in rr]);y=np.array([float(r['msd']) for r in rr]);mask=t<=16
    ax.plot(t[mask],y[mask],color=C[g],lw=2);ax.plot(t[~mask],y[~mask],color=C[g],ls='--',lw=1.5)
ax.axvline(16,color='#cccccc',ls=':');ax.set(xscale='log',yscale='log',xlabel='Lag / elapsed time (h)',ylabel='MSD (µm²)')
note(ax,'Analytic stationary-model prediction; no supplied MSD targets.\nDashed curves extrapolate beyond the 16 h observation.',y=-.32)

ax=fig.add_axes([.71,.039,.265,.146]);ax.axis('off');ax.text(0,1,'J   Additional reported observations and predictions',fontsize=12.5,weight='bold',va='top')
fraction={r['group']:float(r['positive_high_fraction']) for r in sub}
lines=['Reported: both KDs reduce positive/high-motility PACCs.',f"Model: baseline {fraction['P']:.1%} → HIF {fraction['HP']:.1%}; RhoA {fraction['RP']:.1%}.",'High = above baseline PACC endpoint Q75; positive = forward > 0.','', 'Rayleigh check: probability of p < 0.05 in 60 simulated cells:']
for base in ('N','P'):
    probs=[float(next(r['fraction_p_below_05'] for r in ray if r['group']==g)) for g in (base,'H'+base,'R'+base)]
    lines.append(('non-PACC' if base=='N' else 'PACC')+f': {probs[0]:.0%} baseline; {probs[1]:.0%} HIF; {probs[2]:.0%} RhoA.')
lines+=['Reported: baseline/control p < 0.05; both KDs p > 0.05.','Simulation assumes independent cells; replicate effects unavailable.','','Cue-removal prediction (κ = 0; other parameters retained):']
for g in ('N','P'):
    after=float(next(r['directness'] for r in cue if r['group']==g));lines.append(('non-PACC' if g=='N' else 'PACC')+f" median directness {q(g,'directness')[1]:.2f} → {after:.2f}.")
ax.text(0,.86,'\n'.join(lines),va='top',fontsize=9.5,linespacing=1.35,color='#444444')
fig.text(.05,.008,'Source: supplied hosny_2026.md only. 32,768 independent prediction cells/state; exact event integration; nine simulator checks passed. Spatial oxygen responses, temporal adaptation and cell area are outside this model.',fontsize=9,color='#555555')
fig.savefig(OUT/'hosny_poisson_joint.png',dpi=180)
print(OUT/'hosny_poisson_joint.png')
