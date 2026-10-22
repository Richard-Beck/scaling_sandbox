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
    ax.set_yscale('log')
    if arm=='B':
        eligible=frame[frame.model.eq(model)&frame.control.isna()&frame.target_ratio.eq(1)&frame.primary_claim_eligible]
        ax.set_ylim(top=float(eligible.resolution_lower.max())*2.)
    ax.set_title(NAMES[model],color=COLORS[model],fontsize=10,weight='bold');clean(ax)
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
    ax.set_ylim(3.65,-.65);ax.set_xlim(-.15,1.65);ax.set_xticks([0,.5,1,1.5]);ax.set_xlabel('log₁₀ time ratio');ax.set_title(NAMES[model],color=COLORS[model],fontsize=10,weight='bold');clean(ax)

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
    im=previous.kymographs(axes[:2],A,'goryachev',[2,8]);cb=fig.colorbar(im,cax=fig.add_subplot(kg[0,2]));cb.ax.set_title('Signal',fontsize=7,pad=7)
    im=b_kymographs(axes[2:],B);fig.colorbar(im,cax=fig.add_subplot(kg[0,5]),label='Active signal')
    axes[1].set_ylabel('');axes[2].set_ylabel('');axes[3].set_ylabel('')
    for i,m in enumerate(NAMES):curves(fig.add_subplot(gs[2,i*3:(i+1)*3]),A,m,'A1')
    for i,m in enumerate(NAMES):curves(fig.add_subplot(gs[3,i*3:(i+1)*3]),B,m,'B')
    for i,m in enumerate(NAMES):decomposition(fig.add_subplot(gs[4,i*3:(i+1)*3]),D,m)
    # Position labels in row margins using actual gridspec bounds.
    for row,label in [(1,'B  Native histories: driven handover (left); time after withdrawal (right)'),(2,'C  Driven handover against initial edge gap'),(3,'D  Cue-free resolution against measured withdrawal gap · equal amounts'),(4,'E  Size at matched gap (solid) + extra separation (hatched)')]:
        fig.text(.075,gs[row,0].get_position(fig).y1+.029,label,fontsize=11,weight='bold')
    handles=[Line2D([],[],color=color('wave_pinning',i),ls=STYLES[i],label=s) for i,s in enumerate(['small','middle','large'])]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.27,.025),ncol=3,frameon=False)
    fig.text(.62,.035,'Open triangle + arrow: lower bound.  ? Both times censored.\nTop ticks: each size’s antipodal-placement gap. ×: unavailable handover.',fontsize=8,ha='center')
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'figure3_main.{ext}',dpi=200)
    plt.close(fig)

def describe(r):
    if r.status=='indeterminate_both_censored':return 'indeterminate (both censored)'
    if r.status=='lower_bound':return f'≥{r.lower:.3g}'
    if r.status=='upper_bound':return f'≤{r.upper:.3g}'
    if r.status=='exact':return f'{r.point:.3g}'
    return f'{r.lower:.3g}–{r.upper:.3g}'

def supplement(D,S,X):
    from pypdf import PdfReader,PdfWriter
    new=OUT/'supplement_additions.pdf'
    with PdfPages(new) as pdf:
        fig=plt.figure(figsize=(11,12));fig.text(.07,.95,'Figure 3 round 2 · decomposition and censoring',fontsize=15,weight='bold')
        paragraphs=[
            'At matched edge gap, larger native size accounts for most of the increase in driven handover time. Extra accessible separation contributes relatively little. Cue-free competition has a stronger separation response and remains unresolved in several mass-conserved-model conditions.',
            'For each model and assay, the matched gap is the largest gap within both measured size curves. In this dataset it is the small-size antipodal-placement endpoint. In B, the coordinate is the measured gap at withdrawal, not the initial cue gap. The endpoint called antipode is the trajectory seeded by the initial antipodal cue: its withdrawal region centers can shift.',
            'A1 uses native continuous complete-region handover roots. B uses saved-sample lower and upper time brackets. Linear interpolation is restricted to adjacent measured eligible placements. Interval endpoints are interpolated separately; any censored upper endpoint stays infinite. Missing or topology-ineligible values are never converted into bounds.',
            'A ratio with a censored numerator and resolved denominator has a lower bound. A resolved numerator and censored denominator gives an upper bound. Two censored times give no finite one-sided ratio bound; the result is marked indeterminate. In particular, dividing two horizon caps is invalid. Finite bars use the geometric center of the resulting interval, with the full interval recorded in CSV.',
            'The exact product identity uses the same interpolated matched-gap time in both terms. For finite B intervals the representative values obey the identity; multiplying independently propagated marginal interval limits would ignore this shared denominator and numerator. Total intervals are computed directly from endpoint times.',
            f'Optional extensions awaiting user review: {len(X)} Holmes ratio-1 histories (5 small, 4 middle), each completed to 20,000 and eligible for a proposed 80,000 horizon. They lie within the common withdrawal-gap range. No extension has been submitted. These extensions may still leave large-size matched-gap or antipodal times censored.',
            'The following pages retain BOTH superseded main pages, including the complete translation page, then the entire seven-page existing supplement. All four native-profile appendix volumes remain linked in the existing reader; no model inputs, histories, or prior supplements are replaced.'
        ]
        y=.90
        for p in paragraphs:
            text='\n'.join(textwrap.wrap(p,105));fig.text(.07,y,text,va='top',fontsize=9,linespacing=1.5);y-=.024*(text.count('\n')+1)+.024
        pdf.savefig(fig);plt.close(fig)
        fig,ax=plt.subplots(figsize=(11,10));ax.axis('off');rows=[]
        for r in S.itertuples():
            part=D[D.arm.eq(r.arm)&D.model.eq(r.model)&D['size'].eq(r.size)]
            vals=[describe(part[part.term.eq(t)].iloc[0]) for t in ['size_at_matched_gap','extra_separation','total']]
            rows.append([NAMES[r.model],r.arm,f'{r.size:g} / {r.reference_size:g}',f'{r.matched_gap:.4g}',*vals])
        table=ax.table(cellText=rows,colLabels=['Model','Assay','Size / small','Matched gap','Size term','Extra separation','Total'],loc='center',cellLoc='left',colWidths=[.14,.07,.11,.11,.17,.22,.17])
        table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.9)
        fig.suptitle('All decomposition values · time ratios with interval bounds',fontsize=13)
        fig.text(.05,.07,'Intervals reflect saved-time brackets, not biological uncertainty. ≥ is a one-sided lower bound.\nAll interval endpoints, interpolation weights and source slots appear in the companion CSVs.',fontsize=9)
        pdf.savefig(fig);plt.close(fig)
    writer=PdfWriter()
    for p in [new,OLD/'supporting/vectors/figure3_distance_main.pdf',OLD/'supporting/vectors/figure3_distance_supplement.pdf']:
        reader=PdfReader(p)
        for page in reader.pages:writer.add_page(page)
    with (OUT/'figure3_supplement.pdf').open('wb') as f:writer.write(f)
    return len(writer.pages)

def documents(A,B,D,S,X,pages):
    h=B[B.control.isna()&B.target_ratio.eq(1)&B.primary_claim_eligible&B.model.eq('holmes_model3')]
    n=int(h.resolution_upper.isna().sum());den=len(h)
    legend=f'''Separation slows cue-free resolution more strongly than driven handover; enlarged cells can sustain larger separations where coexisting regions persist. Times and gaps use native units.

A. A native cue is rotated at fixed cell size. The initial gap runs from the incumbent 10%-contrast edge to the cue support near edge. The inset holds the physical gap at 4 across small and large Wave-pinning rings; changing native widths also changes center separation. A1 measures continuous-cue complete-region handover. B withdraws the cue at complete-region challenger/incumbent excess ratio 1 and measures the first single resolved peak.

B. Left: Goryachev continuous-cue histories at circumference 37.8 and gaps 2 and 8; white dashed lines mark native continuous handover roots. Right: Wave-pinning at circumference 47.6, measured withdrawal gaps 2.78 and 14.84; the near condition resolves within 1140–1160, while the far condition remains unresolved through 20,000. Colour scales are shared within each model; the different assay time windows are shown explicitly. Native species are active chemical species. These maps display fields, not cell trajectories.

C. Goryachev handover rises with gap, whereas Holmes is nearly flat. Wave-pinning and Spring rise near the edge and then plateau; Wave-pinning does have a distance effect. Open circles are borderline placements. An × at the axis floor denotes an unavailable genuine two-region handover, not a time estimate. Size is consistently encoded by lightness and line style. The coloured top ticks mark initial antipodal-placement gaps.

D. Equal-amount cue-free resolution rises steeply with gap in Wave-pinning, Holmes and Spring, with unresolved large-size conditions retained as open triangles and upward arrows. Holmes has {n}/{den} eligible ratio-1 main-curve times censored. Finite brackets are drawn as intervals, and lines do not bridge censored or ineligible placements. Top ticks mark each size's measured withdrawal gap in the initial-antipodal-cue trajectory; they do not assert exactly antipodal withdrawal centers. Goryachev includes a nonmonotone near-edge response. All alternate, merged and unseeded cases remain in the CSV and original native reader.

E. Each larger size is compared with the smallest size at the largest shared measured gap, using piecewise-linear interpolation inside measured support. Solid bars represent size at matched gap; hatched bars represent extra separation to the larger-size antipodal-placement endpoint. Bars give log₁₀ time ratios; finite B interval limits are marked. Open-ended bars retain one-sided bounds. A question mark means both relevant times are censored and the ratio is indeterminate. Such bars cannot be stacked into a measured total. Total ratios and all interpolation sources are tabulated separately. Size at matched gap dominates the driven-handover time increase in all four models.

The complete translation page, earlier main plate, all three B seeding ratios, A2 pulse brackets, 50%-edge reclassification, direct heading checks and numerical controls remain in the supplement and original reader. The finer-mesh Goryachev B ratio-1 gap-4 comparison has about 4% resolution-time sensitivity. Finite horizons establish lower bounds rather than permanent arrest; scalar/spatial fate and native-field stationarity are separate. No model retuning or new simulation is included. Simulated activity heading, observed migration and conditional OU persistence remain distinct.
'''
    (OUT/'legend.md').write_text(legend)
    (OUT/'README.md').write_text(f'''# Figure 3 · round 2 (22 October 2026)

**Message:** separation slows cue-free resolution more than driven takeover. At matched gap, size dominates the driven-handover increase. Enlarged cells access farther-separated coexisting regions.

[Main PNG](figure3_main.png) · [PDF](figure3_main.pdf) · [SVG](figure3_main.svg) · [complete supplement ({pages} pages)](figure3_supplement.pdf) · [legend](legend.md).

Since 20 October, the main figure has five panels, a same-gap size inset, the requested Goryachev/Wave-pinning native contrasts, separate driven and cue-free curves, and the new matched-gap decomposition. Model colours now match Figure 4 round 2. The translation page and both superseded main pages are retained in the supplement, followed by the complete existing supplement. The [existing reader](../../revision_program_20261020/figure3_distance/index.html) retains all four native-profile appendix volumes and 1,150 native-history images.

**Qualifications:** cue-free B times are saved-sample intervals, not exact roots. Withdrawal gaps differ from initial cue gaps; an antipodal-placement endpoint need not have exactly antipodal withdrawal centers. Interpolation is confined to adjacent eligible measured placements and assumes piecewise-linear curves. Numerator censoring yields a lower bound only with a resolved denominator; two censored times leave the ratio indeterminate. Holmes B has {n}/{den} eligible main-curve bounds; its decomposition is indeterminate. Wave-pinning large B has a lower-bounded size/total effect and indeterminate extra-separation ratio. Spring large B has lower-bounded size/total terms and an indeterminate extra-separation ratio because matched-gap interpolation touches a censored neighbor. Main bars never substitute horizon caps for unresolved times. Goryachev B includes a nonmonotone response and an approximately 4% fine-mesh timing qualification. Topology-ineligible/unseeded cases remain explicit in the evidence CSVs; handover, native-field stationarity, activity heading, migration summaries and conditional OU estimates remain separate.

**Review stop:** {len(X)} optional Holmes extensions qualify under section 2.4 (5 smallest, 4 middle), from 20,000 to proposed 80,000. [Candidate CSV](optional_extension_candidates.csv). **No extension submitted.** Large-size terms may remain censored even after these extensions.

**CSVs:** [decomposition terms](decomposition_terms.csv), [summary](decomposition_summary.csv), [interpolation provenance](interpolation_queries.csv), [A1 plotted values](A1_plotted_values.csv), [B ratio-1 evidence](B_ratio1_plotted_values.csv). Infinite interval upper endpoints mean censoring; zero-to-infinity ratio intervals mean indeterminate. B finite representative values are geometric centers of brackets; full endpoints are retained. A1 unavailable topology is never censoring.

Rebuild: `/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/revision_round2_figure3_20261022.py` from repository root. Existing completed outputs only; no SLURM work. Generated at {datetime.now(timezone.utc).isoformat()}.
''')
    report=dict(generated_utc=datetime.now(timezone.utc).isoformat(),source_sha256={n:sha(synthesis.DESTINATION/n) for n in ['A1_complete_region_evidence.csv','B_complete_region_evidence_and_geometry.csv']},
        model_colors=COLORS,main_panels=5,supplement_pages=pages,existing_supplement_pages_preserved=True,optional_extension_count=len(X),extensions_submitted=0,
        holmes_B_censored=n,holmes_B_eligible=den,decomposition_rows=len(D),finite_product_identity='passed',bounds='interval arithmetic; no horizon point estimates',generator_sha256=sha(__file__))
    (OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    sources=[]
    for arm,frame,slots in [('A1',A,[48,54]),('B',B,[67,94])]:
        if arm=='A1':
            part=A[A.model.eq('goryachev')&A.size_parameter.eq(37.8)&A.control.isna()&A.gap.isin([2.,8.])]
        else:part=B[B.slot.isin(slots)&B.control.isna()]
        for row in part.itertuples():
            status=native.read(Path(row.native_status))
            for phase in status['phases']:
                if arm=='B' and phase['kind']!='cue_free_resolution_tail':continue
                root=synthesis.PRIMARY if arm=='A1' else synthesis.B
                sources.append(dict(assay=arm,model=row.model,slot=row.slot,archive=str(root/phase['path']),sha256=phase['sha256'],absolute_start=phase['absolute_start'],plotted_cap=500 if arm=='A1' else 20000))
    pd.DataFrame(sources).to_csv(OUT/'kymograph_sources.csv',index=False)
    topic=OUT.relative_to(ROOT.parents[1])
    content=dict(id='figure3',title='Figure 3. '+TITLE,status='round2_review',main_pages=1,supplement_pages=pages,
        assets=[dict(path=str(topic/'figure3_main.png'),caption=TITLE,page_size=[12,16],minimum_printed_width_inches=10,page_title='Figure 3. '+TITLE)],
        main_pdf=str(topic/'figure3_main.pdf'),supplement_pdf=str(topic/'figure3_supplement.pdf'),reader=str(topic/'index.html'),
        legend=[p.strip() for p in legend.split('\n\n') if p.strip()],
        role='Separate matched-gap size effects from extra accessible separation for driven handover and cue-free resolution.',
        purpose='Existing-output recomposition and interval-safe decomposition without model retuning.',
        limitations=['Censored ratios remain bounds or indeterminate; no horizon point substitution.','Nine optional Holmes horizon extensions await user review and have not been submitted.'],
        provenance=report)
    (OUT/'content.json').write_text(json.dumps(content,indent=2)+'\n')
    import html
    body=''.join('<p>'+html.escape(p)+'</p>' for p in content['legend'])
    (OUT/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Figure 3 round 2</title><style>body{max-width:1100px;margin:2rem auto;font:17px/1.5 system-ui;padding:1rem}img{width:100%}a{color:#2878b5}</style></head><body><h1>'+html.escape(TITLE)+'</h1><p>Generated '+report['generated_utc']+'; producer complete.</p><p><a href="figure3_main.pdf">Main PDF</a> · <a href="figure3_main.svg">SVG</a> · <a href="figure3_supplement.pdf">Complete supplement</a> · <a href="decomposition_terms.csv">Decomposition CSV</a> · <a href="optional_extension_candidates.csv">9 optional extensions: not submitted</a></p><img src="figure3_main.png" alt="Figure 3 main">'+body+'<p><a href="../../revision_program_20261020/figure3_distance/index.html">Preserved original reader and full native appendices</a></p></body></html>')


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    A,B,D,S,Q,X=evidence()
    main_figure(A,B,D)
    pages=supplement(D,S,X)
    documents(A,B,D,S,X,pages)
    print(json.dumps(dict(output=str(OUT),optional_extensions=len(X),supplement_pages=pages)))
