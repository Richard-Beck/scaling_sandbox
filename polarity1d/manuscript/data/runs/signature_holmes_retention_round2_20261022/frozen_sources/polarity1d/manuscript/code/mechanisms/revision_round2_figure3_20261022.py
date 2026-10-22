"""Figure 3 round 2: existing native evidence, interval-safe decomposition.

No integrations or scheduler submissions. Rebuild in scaling_sandbox environment.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, shutil, textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import to_rgb
import signature_figure3_native_20261020 as native
import signature_figure3_synthesis_20261020 as synthesis
import signature_figure3_render_20261020 as previous

ROOT=native.ROOT
OUT=ROOT/'figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure3'
OLD=previous.OUT
COLORS=dict(wave_pinning='#2878b5',goryachev='#c83b3b',holmes_model3='#36944c',spring='#c59a25')
NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes',spring='Spring')
TITLE='Separation slows the resolution of coexisting regions more than it slows driven takeover'
STYLES=['-','--',':']

def color(model,i):
    c=np.array(to_rgb(COLORS[model]));return tuple(c+(1-c)*[.55,.25,0][i])

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(ax):previous.clean(ax)
def finite(v):return np.isfinite(v)
def bound(ax,x,y,c,horizontal=False):
    ax.plot(x,y,marker='>' if horizontal else '^',ms=5,mfc='white',mec=c,zorder=6)
    end=(x+.16,y) if horizontal else (x,y*1.5)
    ax.annotate('',xy=end,xytext=(x,y),arrowprops=dict(arrowstyle='->',color=c,lw=.9))

def query(frame,gap):
    """Interpolate endpoint intervals; never replace an infinite endpoint by cap."""
    d=frame.sort_values('x');x=d.x.to_numpy()
    assert x.min()-1e-8<=gap<=x.max()+1e-8
    hit=np.flatnonzero(abs(x-gap)<1e-8)
    if len(hit):
        r=d.iloc[int(hit[0])];return dict(lower=float(r.lo),upper=float(r.hi),slots=[int(r.slot)],fraction=0.)
    j=int(np.searchsorted(x,gap));a,b=d.iloc[j-1],d.iloc[j]
    f=float((gap-a.x)/(b.x-a.x))
    lo=(1-f)*a.lo+f*b.lo
    hi=np.inf if not finite(a.hi) or not finite(b.hi) else (1-f)*a.hi+f*b.hi
    return dict(lower=float(lo),upper=float(hi),slots=[int(a.slot),int(b.slot)],fraction=f)

def ratio(a,b):
    lo=a['lower']/b['upper'] if finite(b['upper']) else 0.
    hi=a['upper']/b['lower'] if b['lower']>0 else np.inf
    status=('indeterminate_both_censored' if lo==0 and not finite(hi) else 'lower_bound' if not finite(hi) else 'upper_bound' if lo==0 else 'exact' if abs(hi-lo)<1e-10 else 'finite_interval')
    point=np.sqrt(lo*hi) if lo>0 and finite(hi) else np.nan
    return dict(lower=lo,upper=hi,point=point,status=status)

def evidence():
    A=pd.read_csv(synthesis.DESTINATION/'A1_complete_region_evidence.csv')
    B=pd.read_csv(synthesis.DESTINATION/'B_complete_region_evidence_and_geometry.csv')
    A.to_csv(OUT/'A1_plotted_values.csv',index=False)
    B[B.target_ratio==1].to_csv(OUT/'B_ratio1_plotted_values.csv',index=False)
    a=A[A.control.isna() & (A.classification=='external')].copy()
    a['x']=a.gap;a['lo']=a.handover_time.fillna(a.handover_lower_bound);a['hi']=a.handover_time
    # Missing topology is unavailable; finite lower bounds alone are censoring.
    a.loc[a.lo.notna() & a.hi.isna(),'hi']=np.inf
    b=B[B.control.isna() & (B.target_ratio==1) & B.primary_claim_eligible].copy()
    b['x']=b.measured_gap;b['lo']=b.resolution_lower;b['hi']=b.resolution_upper.fillna(np.inf)
    records=[];queries=[];summary=[];candidates=[]
    for arm,d in [('A1',a),('B',b)]:
      for model in NAMES:
        sizes=native.SIZES[model];small=d[d.model.eq(model)&d.size_parameter.eq(sizes[0])]
        # The measured endpoint is the actual gap at the initial antipodal placement.
        gs=float(small.x.max());s=query(small,gs)
        for index,size in enumerate(sizes[1:],1):
            large=d[d.model.eq(model)&d.size_parameter.eq(size)];gl=float(large.x.max());gm=min(gs,gl)
            sm=query(small,gm);lm=query(large,gm);la=query(large,gl)
            assert abs(gm-gs)<1e-8, 'Generalized total requires small endpoint equals common gap'
            for name,v in [('small_matched',sm),('large_matched',lm),('large_antipode',la)]:
                queries.append(dict(arm=arm,model=model,size=size,query=name,gap=gm if name!='large_antipode' else gl,
                    lower=v['lower'],upper=v['upper'],source_slots=json.dumps(v['slots']),interpolation_fraction=v['fraction']))
            terms={key:ratio(n,dn) for key,n,dn in [('size_at_matched_gap',lm,sm),('extra_separation',la,lm),('total',la,s)]}
            row=dict(arm=arm,model=model,reference_size=sizes[0],size=size,size_index=index,matched_gap=gm,small_antipode_gap=gs,large_antipode_gap=gl)
            for term,result in terms.items():
                records.append(dict(**row,term=term,**result))
                for k,v in result.items():row[term+'_'+k]=v
            if all(finite(v['point']) for v in terms.values()):
                np.testing.assert_allclose(terms['size_at_matched_gap']['point']*terms['extra_separation']['point'],terms['total']['point'],rtol=1e-12)
            summary.append(row)
        if arm=='B':
            # Optional 2.4: only smallest/middle sizes inside the common matched range.
            for size in sizes[:2]:
                for r in d[d.model.eq(model)&d.size_parameter.eq(size)&d.hi.eq(np.inf)&d.x.le(gs+1e-8)].to_dict('records'):
                    candidates.append(dict(model=model,size=size,slot=r['slot'],measured_gap=r['x'],matched_gap_max=gs,
                        completed_horizon=r['lo'],proposed_horizon=4*r['lo'],native_status=r['native_status'],submission_status='NOT_SUBMITTED_REQUIRES_USER_REVIEW'))
    D=pd.DataFrame(records);S=pd.DataFrame(summary);Q=pd.DataFrame(queries);X=pd.DataFrame(candidates)
    D.to_csv(OUT/'decomposition_terms.csv',index=False);S.to_csv(OUT/'decomposition_summary.csv',index=False);Q.to_csv(OUT/'interpolation_queries.csv',index=False);X.to_csv(OUT/'optional_extension_candidates.csv',index=False)
    return A,B,D,S,Q,X

def curves(ax,frame,model,arm):
    for i,size in enumerate(native.SIZES[model]):
        c=color(model,i)
        d=frame[frame.model.eq(model)&frame.size_parameter.eq(size)&frame.control.isna()].copy()
        if arm=='A1':
            d=d[d.classification.ne('internal')].sort_values('gap');x='gap';lo='handover_time';hi='handover_time'
            valid=d[d.classification.eq('external')]
            ax.plot(valid.gap,valid.handover_time,STYLES[i],color=c,lw=1.3)
            for r in d.to_dict('records'):
                if finite(r['handover_time']):ax.plot(r[x],r[lo],'o',ms=3,mfc='white' if r['classification']=='borderline' else c,mec=c)
                elif finite(r['handover_lower_bound']):bound(ax,r[x],r['handover_lower_bound'],c)
                else:ax.plot(r[x],.025,'x',transform=ax.get_xaxis_transform(),color=c,ms=3)
            end=float(valid.gap.max())
        else:
            d=d[d.target_ratio.eq(1)].sort_values('initial_gap');eligible=d[d.primary_claim_eligible]
            for r in eligible.to_dict('records'):
                x0=r['measured_gap'];l=r['resolution_lower'];h=r['resolution_upper']
                if not finite(h):bound(ax,x0,l,c)
                elif l>0:
                    mid=np.sqrt(l*h);ax.errorbar(x0,mid,yerr=[[mid-l],[h-mid]],fmt='o',ms=3,color=c,capsize=1,lw=.7)
            for a,b in zip(d.to_dict('records')[:-1],d.to_dict('records')[1:]):
                if all(r['primary_claim_eligible'] and finite(r['resolution_upper']) and r['resolution_lower']>0 for r in (a,b)):
                    ax.plot([a['measured_gap'],b['measured_gap']],[np.sqrt(r['resolution_lower']*r['resolution_upper']) for r in (a,b)],STYLES[i],color=c,lw=1.1)
            end=float(eligible.measured_gap.max())
        ax.plot(end,.985,marker='|',ms=8,mew=1.4,color=c,transform=ax.get_xaxis_transform(),clip_on=False)
    ax.set_yscale('log');ax.set_title(NAMES[model],color=COLORS[model],fontsize=10,weight='bold');clean(ax)
    ax.set_xlabel('Initial edge gap' if arm=='A1' else 'Withdrawal edge gap')
    ax.set_ylabel('Handover time' if arm=='A1' else 'Resolution time')

def b_kymographs(axes,B):
    rows=[];matrices=[]
    for slot in [67,94]:
        row=B[B.slot.eq(slot)&B.control.isna()].iloc[0];status=native.read(Path(row.native_status));tt=[];yy=[]
        for phase in status['phases']:
            if phase['kind']!='cue_free_resolution_tail':continue
            archive=synthesis.B/phase['path'];assert sha(archive)==phase['sha256']
            with np.load(archive) as z:t=z['time'].copy();y=z['states'][:,:status['placement']['cells']].copy()
            tt.extend((t+phase['absolute_start']-status['withdrawal_time']).tolist());yy.extend(y)
        matrices.append((np.array(tt),np.array(yy)));rows.append(row)
    vmax=max(y.max() for t,y in matrices)
    for ax,row,(t,y) in zip(axes,rows,matrices):
        n=y.shape[1];x=(np.arange(n)+.5)*47.6/n
        im=ax.pcolormesh(x,t,y,cmap='viridis',vmin=0,vmax=vmax,shading='nearest',rasterized=True)
        ax.set_ylim(20000,0);ax.set_title(f'Wave-pinning · g={row.measured_gap:.2f}',fontsize=9)
        if finite(row.resolution_upper):ax.axhline(row.resolution_upper,color='white',ls='--',lw=.8)
        else:ax.text(.5,.95,'unresolved ≥20,000',ha='center',transform=ax.transAxes,color='white',fontsize=8)
        ax.set_xlabel('Native ring arc');ax.set_ylabel('After withdrawal');clean(ax)
    return im

def decomposition(ax,D,model):
    d=D[D.model.eq(model)]
    for j,(arm,index) in enumerate([('A1',1),('A1',2),('B',1),('B',2)]):
        c=color(model,index);start=0
        for term,hatch in [('size_at_matched_gap',''),('extra_separation','////')]:
            r=d[d.arm.eq(arm)&d.size_index.eq(index)&d.term.eq(term)].iloc[0]
            if r.status=='indeterminate_both_censored':
                ax.text(max(start,.02),j+(.2 if term=='extra_separation' else 0),'? '+('extra' if term=='extra_separation' else 'size'),fontsize=7,va='center',color='#555555');continue
            if r.lower>0 and not finite(r.upper):
                width=np.log10(r.lower);ax.barh(j,width,left=start,height=.56,color=c,edgecolor=c,hatch=hatch,alpha=.65 if hatch else 1)
                bound(ax,start+width,j,c,True);start+=width
            elif finite(r.point):
                width=np.log10(r.point);ax.barh(j,width,left=start,height=.56,color=c,edgecolor='white',hatch=hatch,alpha=.65 if hatch else 1)
                if r.status=='finite_interval':ax.plot([start+np.log10(r.lower),start+np.log10(r.upper)],[j,j],color='#333333',lw=1.1)
                start+=width
            else:ax.text(.01,j,'upper bound',fontsize=6)
    ax.axvline(0,color='#555555',lw=.6);ax.set_yticks(range(4),['A1 middle','A1 large','B middle','B large'],fontsize=7)
    ax.invert_yaxis();ax.set_xlim(-.15,1.65);ax.set_xticks([0,.5,1,1.5]);ax.set_xlabel('log₁₀ time ratio');ax.set_title(NAMES[model],color=COLORS[model],fontsize=10,weight='bold');clean(ax)

def main_figure(A,B,D):
    plt.rcParams.update({'font.size':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(12,16));gs=fig.add_gridspec(5,12,height_ratios=[.8,1.65,1.3,1.3,1.35],hspace=.80,wspace=1.8)
    fig.subplots_adjust(top=.95,bottom=.075,left=.075,right=.95)
    fig.suptitle(TITLE,fontsize=13,weight='bold',y=.988)
    fig.text(.075,.964,'A  Rotating a cue separates distance from native cell size',fontsize=11,weight='bold')
    ref=A[A.model.eq('goryachev')&A.size_parameter.eq(37.8)].iloc[0]
    for sl,g,title in [(slice(0,2),2.,'Near'),(slice(2,4),8.,'Far')]:
        previous.ring(fig.add_subplot(gs[0,sl]),1.,ref.initial_10_percent_width,3.,37.8,g,title)
    for sl,C in [(slice(4,6),30.),(slice(6,8),47.6)]:
        r=A[A.model.eq('wave_pinning')&A.size_parameter.eq(C)].iloc[0]
        ax=fig.add_subplot(gs[0,sl]);previous.ring(ax,C/30,r.initial_10_percent_width,3,C,4,f'C={C:g}');ax.set_xlim(-2.7,2.7);ax.set_ylim(-2.7,2.7)
    ax=fig.add_subplot(gs[0,8:]);ax.axis('off');ax.text(0,.98,'Same physical gap across sizes',weight='bold',va='top',fontsize=9)
    ax.text(0,.72,'A1: continuous cue → handover\nB: equal complete-region amounts,\ncue withdrawal → single-peak resolution\nBlue incumbent; orange cue; g is edge gap.',va='top',fontsize=8,linespacing=1.6)
    kg=gs[1,:].subgridspec(1,6,width_ratios=[1,1,.055,1,1,.055],wspace=.55)
    axes=[fig.add_subplot(kg[0,i]) for i in [0,1,3,4]]
    im=previous.kymographs(axes[:2],A,'goryachev',[2,8]);fig.colorbar(im,cax=fig.add_subplot(kg[0,2]),label='Active signal')
    im=b_kymographs(axes[2:],B);fig.colorbar(im,cax=fig.add_subplot(kg[0,5]),label='Active signal')
    axes[1].set_ylabel('');axes[3].set_ylabel('')
    for i,m in enumerate(NAMES):curves(fig.add_subplot(gs[2,i*3:(i+1)*3]),A,m,'A1')
    for i,m in enumerate(NAMES):curves(fig.add_subplot(gs[3,i*3:(i+1)*3]),B,m,'B')
    for i,m in enumerate(NAMES):decomposition(fig.add_subplot(gs[4,i*3:(i+1)*3]),D,m)
    # Position labels in row margins using actual gridspec bounds.
    for row,label in [(1,'B  Native histories: driven handover (left); cue-free near/far resolution (right)'),(2,'C  Driven handover against initial edge gap'),(3,'D  Cue-free resolution against measured withdrawal gap · equal amounts'),(4,'E  Size at matched gap (solid) + extra separation (hatched)')]:
        fig.text(.075,gs[row,0].get_position(fig).y1+.029,label,fontsize=11,weight='bold')
    handles=[Line2D([],[],color=color('wave_pinning',i),ls=STYLES[i],label=s) for i,s in enumerate(['small','middle','large'])]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.27,.025),ncol=3,frameon=False)
    fig.text(.62,.035,'Open triangle + arrow: lower bound.  ? Both times censored.\nTop ticks: each size’s antipodal-placement gap. ×: unavailable handover.',fontsize=8,ha='center')
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure3_main.{ext}',dpi=200)
    plt.close(fig)
