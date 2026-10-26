"""Recompose round7 Figures4 and5 from frozen data and new threshold-free headings."""
from __future__ import annotations
import argparse, json, os, shutil, re
from pathlib import Path
from datetime import datetime, timezone
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[_k]='1'
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
import round7_persistence_20261026 as A
ROOT=A.ROOT;DEV=A.DEV;OUT=A.OUT;BASE=A.BASE;RUN=A.RUN
MODELS=A.MODELS;NAMES=A.NAMES;COLORS=A.COLORS;WEAK=A.WEAK
SOURCES={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'text.color':'black','axes.titlecolor':'black','axes.labelcolor':'black'})
def read(p):SOURCES[A.rel(p)]=A.sha(p);return pd.read_csv(p)
def num(row,key,default=np.nan):
 try:return float(row[key])
 except (TypeError,ValueError,KeyError):return default
def gray(size,model):
 eq=size*20 if model=='spring_ring' else size
 return str(np.interp(eq,[20,60],[.65,.06]))
def primary(d):return d[(d.arm.fillna('common')=='common')&(d['view'].fillna('main')=='main')] if 'arm' in d else d

def letter(fig,spec,label):
 p=spec.get_position(fig);fig.text(max(.012,p.x0-.045),p.y1+.014,label,fontsize=14,fontweight='bold',color='black',va='bottom')
def facet(ax,text):ax.set_title(text,fontweight='normal',color='black',pad=4)
def size_legend(ax,sizes,model,loc='upper right'):
 ax.legend(handles=[Line2D([],[],color=gray(c,model),lw=1.5,label=f'C = {c:g}') for c in sizes],frameon=False,fontsize=6.7,ncol=len(sizes),loc=loc,handlelength=1,columnspacing=.7)
def export_kymo(folder,name,t,field,C):
 t=np.asarray(t);field=np.asarray(field);x=(np.arange(field.shape[1])+.5)/field.shape[1]
 tt=np.repeat(t,len(x));xx=np.tile(x,len(t))
 dest=folder/'supporting/data'/f'{name}.csv.gz'
 if not dest.exists():pd.DataFrame(dict(time=tt,relative_arc_position=xx,activity=field.ravel())).to_csv(dest,index=False)
 return dict(time=t,field=field,C=float(C))
def draw_kymo(ax,r,label,*,time_scale=1000,vmax=None,ylabel=False):
 t=r['time']-r['time'][0];field=r['field'];C=r['C'];vmax=float(field.max()) if vmax is None else vmax
 im=ax.imshow(field,origin='lower',extent=[0,1,t[0]/time_scale,t[-1]/time_scale],aspect='auto',cmap='Greys',vmin=0,vmax=vmax,rasterized=True,interpolation='nearest')
 ax.set(xticks=[0,.5,1],xlabel='Relative arc position')
 facet(ax,label)
 if ylabel:ax.set_ylabel('Elapsed time / '+str(time_scale))
 return im

def prepare4():
 folder=A.setup(4);data=folder/'supporting/data'
 for name in ('B_null_heading_acf.csv','condition_lag_metrics.csv','local_alpha.csv'):
  d=read(BASE/'figure4/supporting/data'/name)
  if name=='condition_lag_metrics.csv':d=primary(d);d=d[np.isclose(d.f,.005)|np.isclose(d.f,.3)];d=d[d.C_label.isin([1,2,3,20,40,60])]
  if name=='condition_lag_metrics.csv':
   workers=sorted((data/'dense_acf_workers').glob('[0-9][0-9].csv'))
   if len(workers)!=24:raise RuntimeError('Dense early-lag ACF campaign incomplete')
   d=pd.concat([read(p) for p in workers],ignore_index=True)
  d.to_csv(data/name,index=False)
 idx=read(BASE/'figure4/supporting/data/A_kymograph_index.csv');ky=[]
 for r in idx.itertuples(index=False):
  path=Path(r.source);SOURCES[A.rel(path)]=A.sha(path)
  with np.load(path) as z:ky.append(dict(kind=r.model,size=float(r.C_label),record=export_kymo(folder,f'A_{r.model}_C{r.C_label:g}',z['time'],z['field'],z['circumference'])))
 compact=read(RUN/'results/compact_index.csv');examples=[]
 for f in (.005,.3):
  for size in (20,40,60):
   q=compact[(compact.model=='wave_pinning')&np.isclose(compact['size'],size)&np.isclose(compact.fraction,f)&(compact.replicate==0)].iloc[0]
   p=Path(q.exemplar);SOURCES[A.rel(p)]=A.sha(p)
   with np.load(p) as z:
    keep=np.linspace(0,len(z['time'])-1,min(800,len(z['time'])),dtype=int)
    examples.append(dict(f=f,size=size,record=export_kymo(folder,f'D_wave_pinning_f{f:g}_C{size}',z['time'][keep],z['primary'][keep],np.mean(z['circumference']))))
 return folder,ky,examples

def figure4():
 folder,ky,examples=prepare4();data=folder/'supporting/data'
 fig=plt.figure(figsize=(16,21.5));outer=fig.add_gridspec(4,1,left=.065,right=.975,top=.965,bottom=.055,height_ratios=[1.03,1.65,1.22,.95],hspace=.36)
 top=outer[0].subgridspec(1,2,width_ratios=[1.45,1],wspace=.24)
 ga=top[0].subgridspec(1,3,wspace=.35);letter(fig,top[0],'A')
 for j,kind in enumerate(('drift','same_rate','perimeter')):
  pair=ga[j].subgridspec(1,2,wspace=.15)
  for k,c in enumerate((20,60)):
   ax=fig.add_subplot(pair[k]);r=next(x['record'] for x in ky if x['kind']==kind and x['size']==c)
   draw_kymo(ax,r,f'C = {c}',vmax=1,ylabel=j==0 and k==0)
   ax.set_xlabel('')
   if k:ax.tick_params(labelleft=False)
  pos=ga[j].get_position(fig);fig.text((pos.x0+pos.x1)/2,pos.y1+.025,{'drift':'Drift','same_rate':'Same-rate jumps','perimeter':'Perimeter-rate jumps'}[kind],ha='center',fontsize=8)
  fig.text((pos.x0+pos.x1)/2,pos.y0-.025,'Relative arc position',ha='center',fontsize=8)
 gb=top[1].subgridspec(3,1,hspace=.8);letter(fig,top[1],'B');null=read(data/'B_null_heading_acf.csv')
 for j,kind in enumerate(('drift','same_rate','perimeter')):
  ax=fig.add_subplot(gb[j]);facet(ax,{'drift':'Drift, α = 2','same_rate':'Same-rate jumps, α = 0','perimeter':'Perimeter-rate jumps, α = −1'}[kind])
  for c,ls in ((20,'-'),(60,'--')):
   q=null[(null.model==kind)&np.isclose(null.C_label,c)].sort_values('lag');col=gray(c,'wave_pinning')
   ax.plot(q.lag/1000,q.C,color=col,ls=ls,lw=1.2,label=f'C = {c}')
   ax.plot(q.lag/1000,q.analytic_C,color=col,ls=':',lw=.7)
  ax.set(ylim=(-.06,1.04),xlim=(0,24),xticks=[0,10,20],xlabel='Lag / 1000',ylabel='Heading correlation')
  if j<2:ax.set_xlabel('')
  if j==0:ax.legend(frameon=False,fontsize=7,ncol=2,loc='upper right')
 gc=outer[1].subgridspec(2,4,wspace=.32,hspace=.42);letter(fig,outer[1],'C');acf=read(data/'condition_lag_metrics.csv')
 for j,model in enumerate(MODELS):
  for i,f in enumerate((.005,.3)):
   ax=fig.add_subplot(gc[i,j]);facet(ax,NAMES[model]+f', f = {f:g}')
   sizes=(1,2,3) if model=='spring_ring' else (20,40,60)
   for c in sizes:
    q=acf[(acf.model==model)&np.isclose(acf.f,f)&np.isclose(acf.C_label,c)].sort_values('lag');col=gray(c,model)
    ax.plot(q.lag/(5 if model=='spring_ring' else 50),q.C,color=col,lw=1.2)
    ax.fill_between(q.lag.to_numpy()/(5 if model=='spring_ring' else 50),q.C_ci_lower.to_numpy(),q.C_ci_upper.to_numpy(),color=col,alpha=.12,lw=0)
   ax.set(ylim=(-.08,1.04),xlabel='Lag / τ')
   if i==0:ax.set_xlim(0,320)
   else:ax.set_xlim(0,80 if model=='goryachev' else 12)
   if j==0:ax.set_ylabel('Heading correlation')
   if i==0:size_legend(ax,sizes,model,loc='lower left')
 gd=outer[2].subgridspec(2,3,wspace=.22,hspace=.43);letter(fig,outer[2],'D')
 for i,f in enumerate((.005,.3)):
  for j,c in enumerate((20,40,60)):
   ax=fig.add_subplot(gd[i,j]);r=next(r['record'] for r in examples if r['f']==f and r['size']==c)
   vmax=max(z['record']['field'].max() for z in examples if z['f']==f)
   im=draw_kymo(ax,r,f'Wave-pinning, C = {c}, f = {f:g}',vmax=vmax,ylabel=j==0)
   if j:ax.tick_params(labelleft=False)
   if j==2:fig.colorbar(im,ax=ax,pad=.015,fraction=.02,label='Activity u')
 ge=outer[3].subgridspec(1,4,wspace=.32);letter(fig,outer[3],'E');alpha=read(data/'local_alpha.csv')
 for j,model in enumerate(MODELS):
  ax=fig.add_subplot(ge[j]);facet(ax,NAMES[model])
  for guide in (-1,0,2):ax.axhline(guide,color='.78',lw=.7,ls='--')
  for c in sorted(alpha[alpha.model==model].C_label.unique()):
   q=alpha[(alpha.model==model)&np.isclose(alpha.C_label,c)].sort_values('f');col=gray(c,model)
   ax.plot(q.f,q.alpha,color=col,lw=.9,marker='o',ms=2.5)
   ok=q.alpha.notna();v=q[ok];ax.vlines(v.f,v.alpha_ci_lower,v.alpha_ci_upper,color=col,lw=.5,alpha=.7)
  ax.set(xscale='log',xlabel='Mean input fraction f',ylim=(-1.5,4.3));ax.set_xticks([.0001,.01,1],['10⁻⁴','10⁻²','1'])
  if j==0:ax.set_ylabel('Local size exponent α')
  size_legend(ax,(1,2,3) if model=='spring_ring' else (20,40,60),model,loc='lower left')
 for ext in ('png','pdf','svg'):fig.savefig(folder/f'figure4_main.{ext}',dpi=200,facecolor='white')
 display_audit(fig,folder);plt.close(fig)
 docs4(folder)
 return folder

def display_audit(fig,folder):
 fig.canvas.draw();renderer=fig.canvas.get_renderer();bad=[];ignored=set()
 for ax in fig.axes:
  for axis,lim in ((ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())):
   for tick in list(axis.get_major_ticks())+list(axis.get_minor_ticks()):
    if not min(lim)<=tick.get_loc()<=max(lim):ignored.update((id(tick.label1),id(tick.label2)))
 for text in fig.findobj(matplotlib.text.Text):
  if not text.get_visible() or not text.get_text() or id(text) in ignored:continue
  b=text.get_window_extent(renderer).transformed(fig.transFigure.inverted())
  if b.x0<-.005 or b.x1>1.005 or b.y0<-.005 or b.y1>1.005:bad.append(dict(text=text.get_text(),bounds=list(b.bounds)))
 A.write(folder/'supporting/metadata/display_audit.json',dict(generated_utc=A.now(),page_inches=list(fig.get_size_inches()),text_outside_page=bad))
 if bad:raise ValueError(bad)

def docs4(folder):
 legend=[
 'Larger native rings retain heading longer under fluctuating input, with the largest size exponents in the weak-input regime. The baseline examples distinguish geometric persistence from a size-dependent change in physical dynamics.',
 'A, Synthetic activity kymographs for a single drifting peak, same-rate random relocations and relocations whose rate scales with circumference, at C20 and C60. Position is relative arc length; grayscale denotes activity from 0 (white) to 1 (black). Time increases upward. These pedagogical baselines have no biochemical kinetics.',
 'B, Heading correlation in the same three baselines, from 1000 seeds; solid and dashed lines indicate C20 and C60, dotted lines the analytic curves. A fixed physical drift diffusivity gives k proportional to C⁻² (α = 2); size-independent jumps give α = 0, and jumps per unit perimeter give α = −1. The same-rate curves coincide.',
 'C, Native-model heading correlation at common weak f = 0.005 and strong f = 0.3, for C20, C40 and C60 (Spring 1, 2 and 3). Size is a light-to-dark grayscale ramp within each facet. Shading is the 95% confidence interval from 2000 seed-bootstrap resamples of the complete retained native histories, recomputed on a dense early-lag grid. Weak-input facets show 0–320τ; strong-input facets show 0–12τ, Goryachev 0–80τ, to resolve their different correlation scales. τ = 50 native time units, Spring 5. Heading correlation is the mean cosine of the change in arg(P), where P is the first Fourier activity moment divided by total activity. These retained round6 estimates require |P| ≥ 0.1 at both lag endpoints; their coverage is in condition_lag_metrics.csv. They differ from the ungated distributions in Figure5D.',
 'D, Native Wave-pinning examples at f = 0.005 and 0.3 and C20, C40 and C60, selected as replicate 0 without a seed search. Every panel shows the entire retained post-burn-in interval 16000–32000, on relative arc position. Activity scales are common across sizes at each input; the color bars state the separate weak- and strong-input ranges. Weak input mainly shifts a coherent domain; strong input produces a changing patch field. The kymographs are illustrative, not an event classification.',
 'E, Local α, the negative fitted slope of log heading-decorrelation rate against log circumference in adjacent-size neighborhoods, against mean input f. Lines are size neighborhoods centered on C20, C30, C40, C50 and C60 (Spring 1, 1.5, 2, 2.5 and 3); light to dark follows size. Vertical bars are seed-bootstrap 95% CIs. Missing estimates interrupt curves. Horizontal guides are α = 2, 0 and −1. Native f multiplies reference input 0.1, 1, 4 and 6.4 in Wave-pinning, Holmes, Spring and Goryachev respectively; doses and clocks are not calibrated across models.',
 'All simulations reuse the complete, frozen round6 compact campaign (96–384 independent seeds per condition). Local kinetics, basal preparations and solver outputs are unchanged. α and activity-heading persistence are model observables, distinct from migration persistence or fitted OU velocity-memory parameters.'
 ]
 narrative=['Geometry alone increases heading persistence when a domain has unchanged physical drift. Size-independent random jumps do not provide that increase.','At common weak and strong inputs, the displayed native heading correlations decay more slowly for larger rings, conditional on the stated coherence eligibility.','Weak-input size exponents are near the drift baseline in Holmes and Spring and exceed it in Wave-pinning where resolved. The advantage narrows as input strengthens; unavailable estimates and their reasons remain in the compact source tables.']
 limitations=['The common f values compare model regimes, not matched biochemical receptor stimulation; native clocks and Spring length units remain model specific.','Autocorrelation in C uses the original |P| ≥ 0.1 endpoint eligibility mask and is conditional on retained coherence; Figure5D exposes raw headings without that mask.','Kymographs do not establish that every heading change is a discrete relocation. Strong-input patch fields can have small first-mode coherence.','Seed-bootstrap intervals reflect the simulated seed ensemble and do not calibrate biological variation.']
 readme='''# Unexpected

The common weak input f = 0.005 lies near the C20 relocation transition in Goryachev but well below it in Wave-pinning and Holmes. Weak input therefore has very different correlation scales across facets. Strong-input activity is not always a single relocating domain; the Wave-pinning examples show a patch field. This figure describes heading persistence and leaves the drift/jump interpretation to Figure5.

## Evidence and definitions

No new simulations are used. Baseline curves and α confidence intervals preserve the round6 estimates. Model ACFs are recomputed from the same complete native compact histories at dense early lags, with the original |P| ≥ 0.1 endpoint rule and 2000 seed-bootstrap resamples; their original 96–384 seeds are preserved. Baseline fields and seed-0 Wave-pinning examples are re-exported as compact CSVs. C uses the original |P| ≥ 0.1 endpoint eligibility mask, while Figure5D intentionally has none. Native f is a fraction of the model reference input, not a common cellular dose.

## Baselines and model readouts

The three synthetic baselines have α = 2, 0 and −1 by construction. Physical drift is fixed across size in the first case. Jump rate is respectively fixed or proportional to C in the other two. P = M1/M0 is the normalized complex first moment of activity; heading is arg(P). The native local exponent is estimated in the original adjacent-size neighborhoods at a shared lag with eligible endpoint correlation and coverage. Missing α estimates are not censored values and receive no marker.

## Reproduction and provenance

Run `round7_persistence_render_20261026.py --figure4` in the existing scaling_sandbox Python environment. `supporting/metadata/provenance.json` hashes the unchanged source tables and exemplar archives. Data-generation timestamps are stored there. No dropped round6 panels or supplement are generated.
'''
 section(folder,4,'Larger cells maintain heading longer under fluctuating input',legend,narrative,limitations,readme,[16,21.5])

def section(folder,n,title,legend,narrative,limitations,readme,size):
 equations=[]
 if n==5:
  equations=[dict(number='5.1',math=r'a_n(t)=\frac{1}{C}\int_0^C u(s,t)e^{-i q_n s}\,\mathrm{d}s,\qquad q_n=\frac{2\pi n}{C}',description='Complex spatial Fourier coefficient of active concentration u on a ring of physical circumference C; n is the integer mode index and qₙ its physical wavenumber.')]
 math_md=''.join('\n\n$$\n'+e['math']+r'\tag{'+e['number']+'}\n$$\n\n'+e['description'] for e in equations)
 (folder/'main_legend.md').write_text('\n\n'.join(legend)+math_md+'\n');(folder/'scientific_readme.md').write_text(readme)
 stem=folder/f'figure{n}_main'
 obj=dict(id=f'figure{n}',title=f'Figure {n}. '+title,role='Native polarity-heading persistence and mechanistic comparisons',assets=[dict(path=A.rel(stem.with_suffix('.png')),native_pdf_page=A.rel(stem.with_suffix('.pdf')),page_size=size,minimum_printed_width_inches=size[0],page_title='')],main_pdf=A.rel(stem.with_suffix('.pdf')),legend=legend,narrative=narrative,measurements=['Proposed cellular assays measure full-perimeter polarity activity, coherence, boundary displacement and cytosolic-pool mobility in matched WGD− and WGD+ cells. These assays are prospective; no experiment is inferred from the model trajectories.'],limitations=limitations,sources=[dict(label='Main source data and definitions',path=A.rel(folder/'scientific_readme.md')),dict(label='Exact source hashes and generation time',path=A.rel(folder/'supporting/metadata/provenance.json'))],equations=equations)
 A.write(folder/'legend.json',obj)
 A.write(folder/'supporting/metadata/provenance.json',dict(generated_utc=A.now(),source_sha256=SOURCES,code=A.rel(Path(__file__)),code_sha256=A.sha(Path(__file__)),outputs={p.name:A.sha(p) for p in folder.glob('figure*_main.*')},new_native_integration=False))

def prepare5():
 folder=A.setup(5);data=folder/'supporting/data'
 mapping={'A_drift_ratio.csv':BASE/'figure3/supporting/data/raw_native_size_normalized_prediction.csv','B_alpha_bridge.csv':BASE/'figure3/supporting/data/alpha_bridge.csv','C_closure_intervals.csv':BASE/'figure4/supporting/data/H_closure_intervals.csv'}
 for name,src in mapping.items():read(src).to_csv(data/name,index=False)
 # Use the original flicker-suppressed relocation rate only as a lag diagnostic.
 dr=read(BASE/'figure4/supporting/data/G_drift_ratio.csv');raw=read(data/'D_relocation_lag_diagnostics.csv');rows=[]
 for r in raw.to_dict('records'):
  q=dr[(dr.model==r['model'])&np.isclose(dr.C_label,r['size'])&np.isclose(dr.f,r['f'])].iloc[0]
  rate=float(q.lambda_s);r.update(counted_rate=rate,counted_mean_wait=1/rate if rate>0 else np.nan,counted_median_tenure=float(q.median_tenure),lag16_over_counted_mean_wait=r['lag16']*rate,lag16_not_well_below_counted_wait=bool(r['lag16']*rate>=.1),counted_event_definition='round6 dominant-track switch; 5-frame flicker removal, contrast eligibility; not a threshold-free relocation detector')
  rows.append(r)
 pd.DataFrame(rows).to_csv(data/'D_relocation_lag_diagnostics_with_counted.csv',index=False)
 return folder

def panel5A(fig,spec,folder):
 data=folder/'supporting/data';d=read(data/'A_drift_ratio.csv');grid=spec.subgridspec(1,4,wspace=.35)
 for j,model in enumerate(MODELS):
  ax=fig.add_subplot(grid[j]);facet(ax,NAMES[model]);q=d[(d.model==model)&(d.fraction<=.01)&(d.multiple_front_fraction<=.05)].copy()
  ax.axhline(1,color='.65',ls='--',lw=.8)
  if model!='goryachev':
   for k,f in enumerate(sorted(q.fraction.unique())):
    r=q[np.isclose(q.fraction,f)].sort_values('size');color=plt.cm.Blues(.3+.65*k/5)
    ax.plot(r['size']/r.reference_size,r.observed_ratio,color=color,lw=.8,marker='o',ms=2.5)
    ok=r.observed_ratio.notna();ax.vlines((r['size']/r.reference_size)[ok],r.observed_ratio_ci_lower[ok],r.observed_ratio_ci_upper[ok],color=color,lw=.4)
    ax.plot(r['size']/r.reference_size,r.shared_translation_predicted_ratio,color=color,ls=':',lw=1)
   ax.set(ylim=(.25,1.65))
  else:
   ax.text(.5,.52,'No qualified\ndrift bridge',ha='center',va='center',transform=ax.transAxes,fontsize=8)
   ax.set(ylim=(.25,1.65))
  ax.set(xlabel='C / Cmin',xlim=(.9,3.1))
  if j==0:ax.set_ylabel('Physical drift ratio Ds / Ds(Cmin)')
  if j==1:ax.legend(handles=[Line2D([],[],color='.25',marker='o',lw=.8,label='Measured'),Line2D([],[],color='.25',ls=':',lw=1,label='Predicted')],frameon=False,fontsize=6.6,loc='lower left')


def panel5B(fig,spec,folder):
 d=read(folder/'supporting/data/B_alpha_bridge.csv');ax=fig.add_subplot(spec)
 ax.plot([1.1,3.4],[1.1,3.4],color='.5',lw=.9)
 for model in MODELS:
  q=d[(d.model==model)&d.bridge_eligible.fillna(False)]
  ax.vlines(q.predicted_alpha,q.observed_alpha_ci_lower,q.observed_alpha_ci_upper,color=COLORS[model],lw=.55,alpha=.6)
  ax.scatter(q.predicted_alpha,q.observed_alpha,s=17,c=COLORS[model],marker='o',edgecolors='none',label=f'{NAMES[model]} (n = {len(q)})')
 ax.set(xlabel='Predicted local α',ylabel='Observed local α',xlim=(1.1,3.4),ylim=(.5,4.2));ax.legend(frameon=False,fontsize=7,loc='upper left',ncol=2)
 ax.text(.98,.02,'f pooled; one symbol per eligible condition',ha='right',transform=ax.transAxes,fontsize=7)


def interval(ax,x,lo,hi,color):
 ymin,ymax=ax.get_ylim()
 if not np.isfinite(lo) and not np.isfinite(hi):return
 if np.isfinite(lo) and np.isfinite(hi) and lo>0 and hi>0:
  a,b=max(lo,ymin),min(hi,ymax)
  if b>=a:ax.vlines(x,a,b,color=color,lw=.85)
  mid=np.sqrt(lo*hi)
  if ymin<=mid<=ymax:ax.plot(x,mid,'o',ms=2.5,color=color)
  for bound,direction in ((lo,-1),(hi,1)):
   if bound<ymin or bound>ymax:
    y=ymin*1.08 if bound<ymin else ymax/1.08;sign=-1 if bound<ymin else 1
    ax.plot(x,y,'v' if sign<0 else '^',mfc='white',mec=color,ms=4)
    ax.annotate('',xy=(x,y/1.6 if sign<0 else y*1.6),xytext=(x,y),arrowprops=dict(arrowstyle='->',color=color,lw=.65),annotation_clip=False)
 elif np.isfinite(hi) and hi>0:
  y=np.clip(hi,ymin*1.2,ymax/1.2);ax.plot(x,y,'v',mfc='white',mec=color,ms=4)
  ax.annotate('',xy=(x,max(y/2,ymin*.8)),xytext=(x,y),arrowprops=dict(arrowstyle='->',color=color,lw=.65),annotation_clip=False)
 elif np.isfinite(lo) and lo>0:
  y=np.clip(lo,ymin*1.2,ymax/1.2);ax.plot(x,y,'^',mfc='white',mec=color,ms=4)
  ax.annotate('',xy=(x,min(y*2,ymax*1.15)),xytext=(x,y),arrowprops=dict(arrowstyle='->',color=color,lw=.65),annotation_clip=False)


def panel5C(fig,spec,folder):
 d=read(folder/'supporting/data/C_closure_intervals.csv');grid=spec.subgridspec(1,4,wspace=.35)
 for j,model in enumerate(MODELS):
  ax=fig.add_subplot(grid[j]);facet(ax,NAMES[model]);ax.set(xscale='log',yscale='log',ylim=(.008,180),xlim=(7e-5,1.5),xlabel='Mean input fraction f')
  ax.axhspan(.85,1.15,color='.92');ax.axhline(1,color='.6',ls='--',lw=.7)
  q=d[d.model==model]
  sizes=sorted(q.C_label.unique())
  for k,c in enumerate(sizes):
   z=q[np.isclose(q.C_label,c)];col=gray(c,model);off=np.exp(np.linspace(-.07,.07,len(sizes))[k])
   for r in z.to_dict('records'):interval(ax,r['f']*off,num(r,'ratio_lower'),num(r,'ratio_upper'),col)
  ax.set_xticks([.0001,.01,1],['10⁻⁴','10⁻²','1'])
  if j==0:ax.set_ylabel('(Drift + counted losses) / observed rate')


def panel5D(fig,spec,folder):
 d=read(folder/'supporting/data/D_raw_heading_distributions.csv.gz');grid=spec.subgridspec(4,4,wspace=.35,hspace=.68)
 lagcols=['#486aa0','#a76329','#2f3339'];styles=['-','--',':']
 for i,model in enumerate(MODELS):
  sizes=(1,3) if model=='spring_ring' else (20,60)
  for j,(f,c) in enumerate([(WEAK[model],sizes[0]),(WEAK[model],sizes[1]),(.3,sizes[0]),(.3,sizes[1])]):
   ax=fig.add_subplot(grid[i,j]);q=d[(d.model==model)&np.isclose(d.f,f)&np.isclose(d['size'],c)]
   for k,mult in enumerate((1,4,16)):
    z=q[q.lag_over_tau==mult];valid=z.rescaled_density>0
    # Nonuniform physical bins resolve the narrow core without discarding the broad tail.
    ax.plot(z.rescaled_displacement,z.rescaled_density.where(valid),color=lagcols[k],ls=styles[k],lw=.9,label=f'{mult}τ')
   limit=float(np.max(abs(q[q['count']>0].rescaled_displacement)))*1.05;threshold=c*1e-4/np.sqrt(5 if model=='spring_ring' else 50)
   ax.set(xscale='symlog',yscale='log',xlim=(-limit,limit));ax.set_xscale('symlog',linthresh=threshold)
   ax.set_xticks([-limit/2,0,limit/2],[f'{-limit/2:.2g}','0',f'{limit/2:.2g}'])
   facet(ax,f'{NAMES[model]}, f = {f:g}, C = {c:g}')
   if i==3:ax.set_xlabel('Signed physical displacement / √lag')
   if j==0:ax.set_ylabel('Probability density')
   if i==0 and j==0:ax.legend(frameon=False,fontsize=7,ncol=3,loc='upper center',handlelength=1.3,columnspacing=.9)


def figure5(*,allow_pending=False):
 folder=prepare5();fig=plt.figure(figsize=(16,23.5));outer=fig.add_gridspec(5,1,left=.065,right=.975,top=.97,bottom=.045,height_ratios=[1,1.04,2.75,1.55,1.05],hspace=.40)
 first=outer[0].subgridspec(1,2,width_ratios=[2.7,1],wspace=.23)
 letter(fig,first[0],'A');panel5A(fig,first[0],folder)
 letter(fig,first[1],'B');panel5B(fig,first[1],folder)
 letter(fig,outer[1],'C');panel5C(fig,outer[1],folder)
 letter(fig,outer[2],'D');panel5D(fig,outer[2],folder)
 pos=outer[2].get_position(fig);fig.legend(handles=[Line2D([],[],color=c,ls=st,lw=1.2,label=label) for c,st,label in zip(['#486aa0','#a76329','#2f3339'],['-','--',':'],['τ','4τ','16τ'])],loc='lower center',bbox_to_anchor=(.52,pos.y1+.01),ncol=3,frameon=False,fontsize=8)
 from round7_noise_analysis_20261026 import plot_panel_E,plot_panel_F
 letter(fig,outer[3],'E');plot_panel_E(fig,outer[3])
 letter(fig,outer[4],'F');plot_panel_F(fig,outer[4])
 for ext in ('png','pdf','svg'):fig.savefig(folder/f'figure5_main.{ext}',dpi=200,facecolor='white')
 display_audit(fig,folder);plt.close(fig);docs5(folder)
 return folder

def docs5(folder):
 noiselegend=(folder/'noise_legend.md').read_text() if (folder/'noise_legend.md').exists() else ''
 noisereadme=(folder/'noise_readme.md').read_text() if (folder/'noise_readme.md').exists() else ''
 noiselegend=re.sub(r'\$a_n\(t\)=.*?\$', 'aₙ(t), defined below', noiselegend)
 noiselegend=noiselegend.replace(r'$\gamma(q)$','γ(q)').replace(r'$q=2\pi n/C$','q = 2πn/C')
 noisereadme=noisereadme.replace('# Unexpected','### Strong-input findings').replace('\n## ','\n### ')
 legend=[
 'Weak-input physical drift is suppressed beyond the geometric baseline in Wave-pinning, while drift predictions fail quantitatively for Goryachev. Raw heading distributions and strong-input patch statistics expose the limits of treating every change as translation of one coherent domain.',
 'A, Physical drift diffusivity Ds relative to the smallest ring at the same f, against relative circumference. Solid curves and circles are retained native trajectory estimates; dotted curves are the deterministic susceptibility prediction with shared translation and full input covariance. All designed f ≤ 0.01 with multiple-domain fraction ≤ 0.05 enter; color intensity indexes f (10⁻⁴–10⁻²). Vertical bars are retained seed-bootstrap 95% CIs. The geometric null is 1. Goryachev has no qualified drift bridge and is explicitly unavailable here, rather than treated as a drifting domain.',
 'B, Predicted versus observed local α, pooled over eligible f ≤ 0.01, with a single circle symbol throughout and model color. Prediction comes from deterministic native susceptibility and the same adjacent-size neighborhoods as the observed α; no noisy-data fitting or pooling average is applied. Vertical bars are observed seed-bootstrap 95% CIs; the diagonal is equality. Eligibility requires a defined native mean state, linear-response validation, resolved observed α and neighborhood multiple-domain fraction ≤ 0.05. Goryachev contributes no eligible point, a quantitative failure retained in B_alpha_bridge.csv.',
 'C, Drift plus counted takeover angular loss divided by observed heading decorrelation, at every designed input and size. For each counted takeover whose angle is unavailable, 1 − cos(Δθ) is bounded independently in [0,2]. Vertical segments span the resulting partial-identification interval, not a confidence interval; circles mark its geometric midpoint when positive finite endpoints exist. Gray band is 0.85–1.15. Open triangles with arrows show one-sided intervals or portions extending past the displayed range. Light to dark denotes C20–C60 (Spring 1–3). Missing observed rates cannot supply a ratio. If drift is missing its nonnegative lower bound can be used, but the upper limit is unbounded.',
 'D, Raw heading-change probability density at lags τ (blue solid), 4τ (brown dashed) and 16τ (black dotted), with no coherence, contrast or track-identity cut. Weak input near the C20 transition is f = 0.03, 0.02, 0.005 and 0.005 in Wave-pinning, Holmes, Spring and Goryachev; strong input is f = 0.3 for all. Facets separate model, input and C20/C60 (Spring 1/3). Signed Δθ is the principal difference in arg(P), and physical displacement is mean endpoint circumference × Δθ/(2π). Binning is in physical displacement before x and probability density are transformed by 1/√lag and √lag respectively. The signed x axis is symmetric-logarithmic to reveal the narrow core and broad tails; y is logarithmic. τ = 50, Spring 5, in native time units. Zero-count bins are omitted without a pseudocount. All retained 96–384 independent seeds per condition enter. Collapse is descriptive, not a test passed by the data; broad tails can also reflect weak first-mode coherence or overlapping relocations.',
 ]
 if noiselegend:legend.extend([p.strip() for p in noiselegend.split('\n\n') if p.strip()])
 legend.append('Physical polarity drift, counted track switching, raw activity heading, complex Fourier-mode relaxation and migration velocity memory are distinct quantities. Fast inactive diffusion is a mechanistic perturbation; all native comparisons retain native local kinetics.')
 coh=read(folder/'supporting/data/D_coherence_distributions.csv');wait=read(folder/'supporting/data/D_relocation_lag_diagnostics_with_counted.csv')
 lines=['# Unexpected','','At strong input in Goryachev C60, the median |P| is '+f'{coh[(coh.model=="goryachev")&(coh.f==.3)&(coh["size"]==60)].iloc[0]["median"]:.3f}'+'. The raw phase is therefore poorly localized even though it is mathematically defined; its broad displacement density cannot be identified uniquely as relocation of a coherent dominant domain. The distributions are shown without a cut or a pass criterion. The new panel does not by itself prove that Goryachev heading loss is exclusively jump-mediated.','','## Raw heading-change method','','Every retained compact seed for each displayed condition enters (96–384 seeds). Heading is arg(P). Changes are principal signed differences, so the readout is the shortest wrapped displacement over a lag, not cumulative path length. Physical circumference is averaged between endpoints. There is no |P|, contrast, region-count or track-identity mask. Only nonfinite numerical differences would be excluded, and their counts are exported. Variable physical-width bins are symmetric asinh spaced with scale C×10⁻⁴; centers and bounds are in D_raw_heading_distributions.csv.gz. Density normalization uses all finite pairs and the exact bin widths. Rescaling by √lag changes both coordinates and density. The x-axis symmetric-logarithmic scale preserves signed values and the core.','','## Coherence distribution per displayed condition','','Quantiles are over all retained frames and seeds; these are descriptive frame distributions, not independent-sample confidence intervals.','','|Model|f|C label|Seeds|Mean coherence|5%|25%|Median|75%|95%|Fraction coherence<0.1|','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
 for r in coh.to_dict('records'):lines.append(f'|{NAMES[r["model"]]}|{r["f"]:g}|{r["size"]:g}|{r["seeds"]}|{r["mean"]:.3f}|{r["q05"]:.3f}|{r["q25"]:.3f}|{r["median"]:.3f}|{r["q75"]:.3f}|{r["q95"]:.3f}|{r["fraction_below_0_1"]:.3f}|')
 lines.extend(['','## Lag relative to counted relocation interval','','The following diagnostics use the frozen round6 dominant-track event detector after removal of 5-frame flicker. They are separate from the ungated heading distributions. The mean waiting interval is reciprocal pooled counted-switch rate; it is a detector-dependent proxy, not a threshold-free biological relocation time. We call 16τ not well below this interval when 16τ/mean wait ≥ 0.1, purely as a caution about overlapping events. Unscored switches and low coherence limit this comparison, especially strong-input Goryachev C60. No counted switch gives no finite waiting estimate, not proof of no relocations. The companion CSV also reports raw consecutive-frame track changes to expose detector sensitivity.','','|Model|f|C label|16τ|Counted mean wait|16τ / mean wait|Not well below|','|---|---:|---:|---:|---:|---:|---|'])
 for r in wait.to_dict('records'):
  w=r['counted_mean_wait'];ws=f'{w:.1f}' if np.isfinite(w) else 'unresolved (zero counted)'
  lines.append(f'|{NAMES[r["model"]]}|{r["f"]:g}|{r["size"]:g}|{r["lag16"]:g}|{ws}|{r["lag16_over_counted_mean_wait"]:.3f}|'+('yes' if r['lag16_not_well_below_counted_wait'] else 'no / unresolved')+'|')
 lines.extend(['','## Drift, α and partial closure','','A uses native deterministic response with shared translation, no gain recalibration to noisy data. B pools the eligible point set across f but does not average away input dependence. The exclusion inventory is retained in B_alpha_bridge.csv; Goryachev has no eligible bridge estimate. C uses the original complete condition inventory and identifies missing drift as a nonnegative, unbounded quantity. Its intervals are angular-loss bounds, not sampling uncertainty. Domain-translation closure fails in many strong-input conditions; a numerical interval containing one does not validate the event mechanism.','','## New strong-input campaign','',noisereadme,'','## Reproduction and provenance','','`round7_persistence_20261026.extract()` writes D from frozen compact archives; run it with SLURM for bulk file reads. `round7_persistence_render_20261026.py --figure5` assembles all six panels after the new strong-input campaign analysis completes. Exact compact source paths and hashes are in supporting/metadata/D_compact_sources.csv, indexed by the immutable complete-campaign CSV hash. No dropped round6 panel is reproduced as a supplement.'])
 narrative=['Wave-pinning drift falls below unchanged-physical-drift geometry and broadly follows deterministic susceptibility. Holmes and Spring sit closer to size-independent physical drift; Goryachev has no qualified drift bridge.','Raw heading distributions are presented at three lags without a coherence cut. Their interpretation is constrained by broad low-coherence phase changes and by lags that overlap counted relocations.','Strong-input patch-field and mode-rate experiments test inactive-pool transport without changing the native main-comparison kinetics.']
 limitations=['Collapsing lag-rescaled distributions is a descriptive drift signature; a noncollapsing tail is not uniquely a track relocation, especially when |P| is small.','The largest lag can be comparable with or exceed the typical counted-switch interval. Zero detected switches do not prove a long biological waiting time.','C contains detector-conditioned counted events and partial-identification intervals; it is not an exhaustive causal decomposition or a confidence bound.','Native model reference doses, length units and time units are not experimentally calibrated across models.']
 section(folder,5,'Physical drift, raw heading changes and cytosolic transport constrain the persistence mechanism',legend,narrative,limitations,'\n'.join(lines)+'\n',[16,23.5])


def main():
 p=argparse.ArgumentParser();p.add_argument('--figure4',action='store_true');p.add_argument('--figure5',action='store_true');args=p.parse_args()
 if args.figure4:figure4()
 if args.figure5:figure5()
 if not (args.figure4 or args.figure5):figure4();figure5()
if __name__=='__main__':main()
