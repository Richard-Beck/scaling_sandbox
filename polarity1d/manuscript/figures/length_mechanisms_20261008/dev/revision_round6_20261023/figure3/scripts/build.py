"""Collect completed deterministic experiments and render Figure 3 round6."""
from pathlib import Path
import os,sys,json,shutil,hashlib,html,textwrap
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import native_analysis as new
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader,PdfWriter
import figure3_round5_render_20261023 as base
import figure4_round3_estimators_20261022 as estimators
OUT=new.OUT;DATA=OUT/'supporting/data';META=OUT/'supporting/metadata';OLD=base.OUT;OLD_DATA=base.DATA;MS=base.MS
MODELS=base.MODELS;NAMES=base.NAMES;COLORS=base.COLORS
base.DATA=DATA
F4=MS/'figures/length_mechanisms_20261008/dev/revision_round5_20261023/figure4/supporting/data'
STEM='figure3_round6'

def write(p,v):new.write(p,v)
def rel(p):return str(Path(p).relative_to(new.ROOT))
def shade(model,size,sizes=None):
 if sizes is None:sizes=[1,1.5,2,2.5,3] if model=='spring_ring' else [20,30,40,50,60]
 i=float(np.interp(size,sizes,np.arange(len(sizes))));return base.color(model,i,len(sizes))

def collect():
 sources=[]
 for name in ('activity_profiles.csv.gz','activity_profile_metadata.json','half_heading_readouts.csv','regime_lines.csv','undefined_heading_bands.csv','edge_prediction_comparison.csv','size_normalized_prediction_comparison.csv','edge_response_size_folds.csv','pulse_response_summary.csv','edge_response_curves.csv.gz','kernel_pulse_validation.csv','edge_independence_summary.csv','edge_MSD_curves.csv'):
  p=OLD_DATA/name;shutil.copy2(p,DATA/name);sources.append(dict(path=rel(p),sha256=new.old.sha(p),copied_to=rel(DATA/name)))
 manifest=new.read(META/'native_manifest.json');maps=[];experiments=[];curves=[];fractions=[];limits=[];states=[];validations=[];folds=[]
 for c in manifest['jobs']:
  folder=new.NATIVE/f"{c['slot']:03d}";s=new.read(folder/'status.json');assert s.get('complete') and s.get('defined'),s
  C=s['mean_circumference'];ell=.15 if c['platform']=='spring' else 1.5
  z=np.load(folder/'map.npz');x=z['x'];center=float(z['center']);xx=new.wrap(x-center,C);order=np.argsort(xx);e=new.wrap(z['edges']-center,C);g=z['g'];con=z['contribution'];total=float(z['Ds']);noise_cov=z['input_covariance_integral']
  # Raw native and round5 measured-gain versions are distinguished explicitly.
  raw_con=g*np.fft.irfft(np.fft.rfft(noise_cov)*np.fft.rfft(g),n=len(x))*c['mean']**2
  historical_total=total;con=raw_con;total=float(raw_con.sum())
  delta_left=abs(new.wrap(x-z['edges'][0],C));delta_right=abs(new.wrap(x-z['edges'][1],C));left=(delta_left<=2*ell)&(delta_left<=delta_right);right=(delta_right<=2*ell)&(delta_right<delta_left);in_active=(x-z['edges'][0])%C<((z['edges'][1]-z['edges'][0])%C);plateau=in_active&~left&~right;basal=~in_active&~left&~right
  fractions.append(dict(model=c['model'],size=c['size'],mean_circumference=C,Ds=total,raw_native_Ds=raw_con.sum(),left_edge_fraction=con[left].sum()/total,right_edge_fraction=con[right].sum()/total,edge_union_fraction=con[left|right].sum()/total,plateau_interior_fraction=con[plateau].sum()/total,basal_arc_fraction=con[basal].sum()/total,partition_sum=con[left|right|plateau|basal].sum()/total,native_gain_scale=float(z['scale']),historical_calibrated_Ds=historical_total))
  for i in order:maps.append(dict(model=c['model'],size=c['size'],fraction=.001,cells=c['cells'],x=float(xx[i]),native_x=float(x[i]),mean_circumference=C,activity=float(z['activity'][i]),raw_g=float(g[i]),standard_pulse_gain=float(z['pulse_map'][i]),contribution=float(con[i]),normalized_contribution=float(con[i]/total),raw_native_contribution=float(raw_con[i]),raw_native_Ds=float(raw_con.sum()),Ds=total,native_gain_scale=float(z['scale']),region='left_edge' if left[i] else 'right_edge' if right[i] else 'plateau_interior' if plateau[i] else 'basal_arc',left_edge=float(e[0]),right_edge=float(e[1]),matched_B_size=bool(c['size'] in ([2.5,3.15,3.97] if c['platform']=='spring' else [30,37.8,47.6]))))
  ex=pd.read_csv(folder/'experiments.csv');ex.insert(0,'size',c['size']);ex.insert(0,'model',c['model']);ex['fraction']=.001;ex['normalized_map_error']=abs(ex.shared_gain-ex.predicted_gain)/max(abs(z['pulse_map']));experiments.append(ex)
  cv=pd.read_csv(folder/'curves.csv.gz');cv=cv.merge(ex[['experiment','kind','label','amplitude','duration','dose']],on='experiment',validate='many_to_one');cv.insert(0,'size',c['size']);cv.insert(0,'model',c['model']);curves.append(cv)
  for k in ('mean_stationary_rhs','mean_tail_scaled_change','mean_contrast'):validations.append(dict(model=c['model'],size=c['size'],check=k,value=s[k]))
  if 'amplitude_limit' in s:limits.append(dict(model=c['model'],size=c['size'],**s['amplitude_limit']))
  folds.append(dict(model=c['model'],size=c['size'],map_amplitude=max(abs(z['pulse_map'])),edge_gain=float(ex[(ex.kind=='map_validation')&(ex.label=='right_edge')].shared_gain.iloc[0])))
  states.append(dict(model=c['model'],size=c['size'],state_sha256=new.old.sha(folder/'mean_state.npz'),map_sha256=new.old.sha(folder/'map.npz'),status_sha256=new.old.sha(folder/'status.json'),source_slot=c['source_slot'],source_size=c['source_size'],exact_round5_state=c['exact'],mean_circumference=C,native_gain_scale=float(z['scale'])))
 mp=pd.DataFrame(maps);ex=pd.concat(experiments,ignore_index=True);cv=pd.concat(curves,ignore_index=True);fr=pd.DataFrame(fractions);li=pd.DataFrame(limits);fo=pd.DataFrame(folds)
 existing=pd.read_csv(DATA/'pulse_response_summary.csv');existing=existing[(existing.fraction==.001)&(existing.pulse_edge==1)&(existing.tracked_edge==1)&(existing.multiplier==1)].copy()
 existing['centered_position']=[mp[(mp.model==r.model)&(mp['size']==r.size)].right_edge.iloc[0] for r in existing.itertuples()];existing.to_csv(DATA/'existing_edge_pulses_on_maps.csv',index=False)
 mp.to_csv(DATA/'susceptibility_and_contributions.csv',index=False,float_format='%.10g');ex.to_csv(DATA/'deterministic_experiments.csv',index=False,float_format='%.10g');cv.to_csv(DATA/'deterministic_time_courses.csv.gz',index=False,float_format='%.10g');fr.to_csv(DATA/'drift_contribution_fractions.csv',index=False);li.to_csv(DATA/'linear_range_limits.csv',index=False);fo.to_csv(DATA/'map_amplitudes.csv',index=False);pd.DataFrame(states).to_csv(DATA/'mean_state_provenance.csv',index=False)
 # Weak-input alpha bridge follows the exact fixed-neighbor design of Fig4.
 pred=pd.read_csv(DATA/'edge_prediction_comparison.csv');pred['raw_native_Ds']=2*pred.raw_native_Ds_independent
 # Exact-native normalized prediction replaces empirical pulse normalization on main.
 normal=pd.read_csv(DATA/'size_normalized_prediction_comparison.csv')
 for i,r in normal.iterrows():
  group=pred[(pred.model==r['model'])&(pred.fraction==r['fraction'])];ref=group[group['size']==r['reference_size']].iloc[0];own=group[group['size']==r['size']].iloc[0]
  normal.loc[i,'historical_calibrated_predicted_ratio']=r['shared_translation_predicted_ratio'];normal.loc[i,'raw_native_predicted_Ds']=own.raw_native_Ds
  normal.loc[i,'shared_translation_predicted_ratio']=own.raw_native_Ds/ref.raw_native_Ds if r['defined'] else np.nan
 normal.to_csv(DATA/'raw_native_size_normalized_prediction.csv',index=False)
 obs=pd.read_csv(F4/'local_alpha.csv');obs=obs[(obs.arm=='common')&(obs.view=='main')&(obs.f<=.001)]
 bridge=[]
 for r in obs.itertuples():
  sizes=estimators.fixed_neighborhood(r.C_label,r.model);q=pred[(pred.model==r.model)&(pred.fraction==r.f)&pred['size'].isin(sizes)].sort_values('size')
  good=len(q)==len(sizes) and q.defined.all() and (q.Ds_predicted_full_covariance>0).all();alpha=2-np.polyfit(np.log(q.mean_circumference),np.log(q.raw_native_Ds),1)[0] if good else np.nan
  bridge.append(dict(model=r.model,size=r.C_label,fraction=r.f,neighborhood=str(sizes),n_neighbors=len(sizes),predicted_alpha=alpha,observed_alpha=r.alpha,observed_alpha_ci_lower=r.alpha_ci_lower,observed_alpha_ci_upper=r.alpha_ci_upper,observed_status=r.status,prediction_defined=good))
 br=pd.DataFrame(bridge);br.to_csv(DATA/'alpha_bridge.csv',index=False)
 # Prespecified operational separated-cue plateau: furthest three observed
 # certified two-region half-turns, without smoothing or extrapolation.
 half=pd.read_csv(DATA/'half_heading_readouts.csv');matched=[]
 for model in MODELS:
  key='spring' if model=='spring_ring' else model;sizes=[2.5,3.15,3.97] if model=='spring_ring' else [30.,37.8,47.6]
  q=fo[(fo.model==model)&fo['size'].isin(sizes)].sort_values('size');ref=q.iloc[0]
  records=[]
  for size in sizes:
   b=half[(half.model==key)&np.isclose(half.size_parameter,size)&(half.t50.notna())&(half.topology=='two_resolved_spatially_associated_regions')&(~half.status.str.contains('censored'))].sort_values('gap').tail(3)
   assert len(b),('No separated plateau',model,size)
   row=q[np.isclose(q['size'],size)].iloc[0];records.append(dict(model=model,size=size,plateau_t50=float(b.t50.median()),plateau_min=float(b.t50.min()),plateau_max=float(b.t50.max()),n_plateau_points=len(b),plateau_gaps=','.join(f'{x:g}' for x in b.gap),plateau_selection='Median furthest three finite observed two-region half-turns',map_amplitude=row.map_amplitude,edge_gain=row.edge_gain,map_amplitude_fold=row.map_amplitude/ref.map_amplitude,edge_gain_fold=row.edge_gain/ref.edge_gain))
  for row in records:row['half_turn_time_fold']=row['plateau_t50']/records[0]['plateau_t50']
  matched.extend(records)
 ma=pd.DataFrame(matched);ma.to_csv(DATA/'matched_B_susceptibility.csv',index=False)
 write(META/'provenance.json',dict(generated_utc=new.old.now(),sources=sources,starting_outline=dict(path='manuscript_outline/outline_v11.json',sha256=new.old.sha(new.ROOT/'manuscript_outline/outline_v11.json'),pdf_sha256=new.old.sha(new.ROOT/'manuscript_outline/outline_v11.pdf')),native_manifest_sha256=new.old.sha(META/'native_manifest.json'),worker_source=dict(recorded_sha256=manifest['source_code_sha256'][rel(HERE/'native_analysis.py')],current_sha256=new.old.sha(HERE/'native_analysis.py'),matches_run_record=new.old.sha(HERE/'native_analysis.py')==manifest['source_code_sha256'][rel(HERE/'native_analysis.py')],frozen_copy=rel(OUT/'supporting/native/frozen_sources/native_analysis.py')),downstream_assembler=dict(path=rel(Path(__file__)),sha256=new.old.sha(Path(__file__)),raw_contribution_recomputation='Native worker saved pulse-normalized historical contribution. Assembler recomputes raw native g times its covariance convolution, without that normalization; saved native worker outputs are unchanged.'),new_stochastic_runs=0,retuned_parameters=0,immutable_round5=True))
 return mp,ex,cv,fr,li,br,ma

def setup():
 base.setup();plt.rcParams.update({'font.size':9,'axes.labelsize':9,'axes.titlesize':10,'xtick.labelsize':8,'ytick.labelsize':8})
def label(fig,text,x,y):fig.text(x,y,text,fontsize=15,fontweight='bold',color='black',va='top')
def facets(fig,left=.075,right=.965,bottom=.1,top=.3):return fig.add_gridspec(1,4,left=left,right=right,bottom=bottom,top=top,wspace=.42)
def endpoints(model):return [1,3] if model=='spring_ring' else [20,60]
def map_row(fig,grid,mp,ex,quantity='standard_pulse_gain',all_sizes=False,points=True):
 axes=[]
 for j,model in enumerate(MODELS):
  ax=fig.add_subplot(grid[0,j]);axes.append(ax);ax.set_title(NAMES[model],color='black')
  sizes=sorted(mp[mp.model==model]['size'].unique()) if all_sizes else endpoints(model)
  for size in sizes:
   q=mp[(mp.model==model)&(mp['size']==size)].sort_values('x');col=shade(model,size);ax.plot(q.x,q[quantity],color=col,lw=.75 if quantity=='raw_g' else 1.6,label=f'C={size:g}')
   if points and quantity=='standard_pulse_gain':
    v=ex[(ex.model==model)&(ex['size']==size)&(ex.kind=='map_validation')];ax.scatter(v.centered_position,v.shared_gain,s=16,c=[col],edgecolors='black',linewidths=.35,zorder=4)
    oldp=pd.read_csv(DATA/'pulse_response_summary.csv');oldp=oldp[(oldp.model==model)&(oldp['size']==size)&(oldp.fraction==.001)&(oldp.pulse_edge==1)&(oldp.tracked_edge==1)&(oldp.multiplier==1)]
    if len(oldp):ax.scatter([q.right_edge.iloc[0]],oldp.terminal_gain,s=34,marker='s',facecolors='white',edgecolors=[col],linewidths=1,zorder=5)
   if quantity=='activity':
    # Activity gets a separate letter and axis directly beneath the map.
    ax.lines[-1].set_color(str(.72-.32*sizes.index(size)/max(1,len(sizes)-1)))
   for e in (q.left_edge.iloc[0],q.right_edge.iloc[0]):ax.axvline(e,color=col,ls=':',lw=.6,alpha=.6,zorder=-1)
  ax.axhline(0,color='.7',lw=.5);ax.set_xlabel('Arc position from domain centre');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3));ax.yaxis.offsetText.set_fontsize(7)
  if quantity!='activity':ax.legend(frameon=False,fontsize=8,loc='best',ncol=2 if len(sizes)<4 else 3)
  if j==0:ax.set_ylabel({'standard_pulse_gain':'Terminal shift / dose','activity':'Mean-state activity','normalized_contribution':r'Signed contribution, $c_j/D_s$','raw_g':r'Raw native $g_j$'}[quantity])
 return axes

def time_panel(fig,spec,cv,ex,all_models=False):
 models=MODELS if all_models else ['wave_pinning'];grid=spec.subgridspec(1,len(models)*(1 if all_models else 2),wspace=.35)
 for j,model in enumerate(models):
  sizes=endpoints(model)
  for k,size in enumerate(sizes):
   if all_models and k:continue
   ax=fig.add_subplot(grid[0,j if all_models else k]);chosen=sizes if all_models else [size]
   for size in chosen:
    q=cv[(cv.model==model)&(cv['size']==size)&(cv.kind=='map_validation')&(cv.label=='right_edge')].sort_values('time');dose=q.dose.iloc[0];col=shade(model,size)
    ax.plot(q.time,q.right_displacement/dose,color=col,lw=1.7,label=f'C={size:g}, pulsed edge');ax.plot(q.time,q.left_displacement/dose,color=col,lw=1.2,ls='--',label=f'C={size:g}, opposite edge')
   ax.axvspan(0,5 if model=='spring_ring' else 50,color='.93',zorder=-2);ax.set_xlabel('Time from pulse onset');ax.set_title(NAMES[model]+('' if all_models else f', C={size:g}'));ax.set_ylabel('Edge displacement / dose');ax.legend(frameon=False,fontsize=7,loc='lower right');ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,3))

def drift_panel(fig,grid,normal):
 markers=['o','s','^','v','D','P'];strengths=sorted(normal.fraction.unique())
 for j,model in enumerate(MODELS):
  ax=fig.add_subplot(grid[0,j]);ax.set_title(NAMES[model]);q=normal[(normal.model==model)&normal.defined]
  lo=min(q.observed_ratio.min(),q.shared_translation_predicted_ratio.min())*.65;hi=max(q.observed_ratio.max(),q.shared_translation_predicted_ratio.max())*1.5;ax.plot([lo,hi],[lo,hi],color='.5',lw=.8)
  for r in q.itertuples():
   col=shade(model,r.size);ax.errorbar(r.observed_ratio,r.shared_translation_predicted_ratio,xerr=[[r.observed_ratio-r.observed_ratio_ci_lower],[r.observed_ratio_ci_upper-r.observed_ratio]],fmt=markers[strengths.index(r.fraction)],color=col,ms=4.8,mec=COLORS[model],mew=.4,elinewidth=.55)
  ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(lo,hi);ax.set_ylim(lo,hi);ax.set_xlabel(r'Measured $D_s/D_s(C_{\min})$');
  if j==0:ax.set_ylabel(r'Predicted $D_s/D_s(C_{\min})$')
  if model=='goryachev':ax.text(.03,.98,'own scale; model fails',ha='left',va='top',transform=ax.transAxes,fontsize=8)

def alpha_panel(fig,grid,br):
 for j,model in enumerate(MODELS):
  ax=fig.add_subplot(grid[0,j]);ax.set_title(NAMES[model]);q=br[(br.model==model)&br.observed_alpha.notna()];p=br[(br.model==model)&br.prediction_defined]
  lim=[min(1.2,p.predicted_alpha.min()-.3),max(3.2,p.predicted_alpha.max()+.3)];ax.plot(lim,lim,color='.5',lw=.8)
  for r in q.itertuples():ax.errorbar(r.predicted_alpha,r.observed_alpha,yerr=[[r.observed_alpha-r.observed_alpha_ci_lower],[r.observed_alpha_ci_upper-r.observed_alpha]],fmt={.0001:'o',.0003:'s',.001:'^'}[r.fraction],color=shade(model,r.size),mec=COLORS[model],mew=.4,ms=5,elinewidth=.7)
  if not len(q):ax.text(.5,.48,'No qualified weak-input\nFigure 4 local α',ha='center',va='center',transform=ax.transAxes,fontsize=9)
  ax.set_xlim(lim);ax.set_ylim(min(0.5,lim[0]),max(3.6,lim[1]));ax.set_xlabel('α implied by predicted drift');
  if j==0:ax.set_ylabel('Figure 4 weak-input α')

def main_figure(mp,ex,cv,fr,li,br,ma):
 setup();fig=plt.figure(figsize=(16.2,19.5))
 # A/B reuse exact source coordinate files and the native existing renderers.
 ga=fig.add_gridspec(4,2,left=.065,right=.335,bottom=.735,top=.965,wspace=.28,hspace=.25);base.panel_A(fig,ga)
 gb=facets(fig,left=.425,right=.96,bottom=.795,top=.965);baxes=base.panel_B(fig,gb)
 half=pd.read_csv(DATA/'half_heading_readouts.csv')
 for ax,model in zip(baxes,MODELS):
  key='spring' if model=='spring_ring' else model
  for r in half[(half.model==key)&half.status.str.contains('right_censored')&(~half.no_turn)].itertuples():ax.annotate('',xy=(r.gap,r.lower*1.2),xytext=(r.gap,r.lower),arrowprops=dict(arrowstyle='->',color=shade(model,r.size_parameter,sorted(half[half.model==key].size_parameter.unique())),lw=.7))
 label(fig,'A',.018,.989);label(fig,'B',.38,.989)
 fig.text(.067,.995,'Grey: cue support; diamonds: saddle; | heading',fontsize=8,color='black')
 fig.text(.425,.775,'B: light → dark = C 30, 37.8, 47.6; Spring 2.5, 3.15, 3.97.\nOpen triangles: time lower bounds. Dotted lines: first saddle. Grey: heading undefined.',fontsize=8,va='top')
 for letter,quantity,bottom,top in [('C','standard_pulse_gain',.585,.695),('D','activity',.502,.548),('E','normalized_contribution',.355,.455)]:
  grid=facets(fig,bottom=bottom,top=top);map_row(fig,grid,mp,ex,quantity,points=True);label(fig,letter,.02,top+.025)
 fig.text(.075,.708,'Lines: native map; circles: new pulse tests; open squares: existing edge pulse. Dotted verticals: half-contrast edges.',fontsize=8)
 label(fig,'F',.02,.312);gs=fig.add_gridspec(1,1,left=.075,right=.49,bottom=.211,top=.29);time_panel(fig,gs[0,0],cv,ex)
 fig.text(.555,.285,'Permanent shift, both edges together',fontsize=11,color='black')
 folds=pd.read_csv(DATA/'edge_response_size_folds.csv');gain=dict(zip(folds.model,folds.terminal_gain_large_over_small))
 t=f"Terminal gain Cmax/Cmin:\nWave-pinning {gain['wave_pinning']:.3f}; Holmes {gain['holmes_model3']:.3f}\nSpring {gain['spring_ring']:.3f}; Goryachev {gain['goryachev']:.3f}"
 fig.text(.555,.261,t,fontsize=9,va='top')
 fig.text(.555,.224,'Weak input anywhere shifts the domain by g(x) × dose.\nNoise is many such inputs; its covariance and g predict drift.\nNo drift data or model parameters are fitted.',fontsize=9,va='top')
 label(fig,'G',.02,.183);normal=pd.read_csv(DATA/'raw_native_size_normalized_prediction.csv');drift_panel(fig,facets(fig,bottom=.077,top=.16),normal)
 # A compact second page would separate quantities; keep one complete main plate
 # and give alpha enough vertical space by increasing the canvas below zero.
 # Figure native page is intentionally tall; alpha is its final own row.
 fig.set_size_inches(16.2,23.0)
 # Move every existing axes upward, scaling heights, to reserve the last row.
 for ax in fig.axes:
  pos=ax.get_position();ax.set_position([pos.x0,.16+.84*pos.y0,pos.width,.84*pos.height])
 for txt in fig.texts:
  x,y=txt.get_position();txt.set_position((x,.16+.84*y))
 label(fig,'H',.02,.185);alpha_panel(fig,facets(fig,bottom=.045,top=.16),br)
 fig.savefig(OUT/f'{STEM}.pdf');fig.savefig(OUT/f'{STEM}.png',dpi=155);fig.savefig(OUT/f'{STEM}.svg');plt.close(fig)

def supplement(mp,ex,cv,fr,li,br,ma):
 setup();new_pdf=OUT/f'{STEM}_new_supplement.pdf'
 with PdfPages(new_pdf) as pdf:
  def page(title,func,size=(14,7.5)):
   fig=plt.figure(figsize=size);fig.text(.06,.96,title,fontsize=13,color='black');func(fig);pdf.savefig(fig);plt.close(fig)
  for quantity,title in [('standard_pulse_gain','S1. Standard Gaussian-pulse susceptibility: all five original sizes and B-matched additions'),('raw_g','S2. Raw native cell-input susceptibility (thin lines): all sizes'),('normalized_contribution','S3. Signed contributions to shared-translation drift: all sizes')]:
   page(title,lambda fig,q=quantity:map_row(fig,facets(fig,bottom=.13,top=.83),mp,ex,q,all_sizes=True,points=q=='standard_pulse_gain'))
  page('S4. Native mean-state activity profiles: all original and matched sizes',lambda fig:map_row(fig,facets(fig,bottom=.13,top=.83),mp,ex,'activity',all_sizes=True,points=False))
  page('S5. Edge co-motion and permanent translation: endpoint sizes',lambda fig:time_panel(fig,fig.add_gridspec(1,1,left=.07,right=.97,bottom=.14,top=.83)[0,0],cv,ex,True))
  for kind,title,xcol in [('amplitude','S6. Right-edge terminal gain against amplitude: first >10% departure included','amplitude'),('duration','S7. Fixed-dose right-edge gain against duration','duration')]:
   def draw(fig):
    grid=facets(fig,bottom=.15,top=.82)
    for j,model in enumerate(MODELS):
     ax=fig.add_subplot(grid[0,j]);ax.set_title(NAMES[model]);q=ex[(ex.model==model)&(ex.kind==kind)]
     for size in sorted(q['size'].unique()):
      v=q[q['size']==size].sort_values(xcol);nom=ex[(ex.model==model)&(ex['size']==size)&(ex.kind=='map_validation')&(ex.label=='right_edge')].shared_gain.iloc[0];ax.plot(v[xcol],v.shared_gain/nom,'o-',ms=3,color=shade(model,size),label=f'C={size:g}')
     ax.axhspan(.9,1.1,color='.93',zorder=-2);ax.axhline(1,color='.5',lw=.7);ax.set_xscale('log');ax.set_xlabel('Pulse peak amplitude' if kind=='amplitude' else 'Pulse duration (fixed dose)');ax.legend(frameon=False,fontsize=7)
     if j==0:ax.set_ylabel('Terminal gain / nominal gain')
   page(title,draw)
  def matched(fig):
   grid=facets(fig,bottom=.17,top=.82)
   for j,model in enumerate(MODELS):
    ax=fig.add_subplot(grid[0,j]);ax.set_title(NAMES[model]);q=ma[ma.model==model].sort_values('size');ax.plot(q.map_amplitude_fold,q.half_turn_time_fold,color=COLORS[model],lw=1)
    for r in q.itertuples():ax.scatter(r.map_amplitude_fold,r.half_turn_time_fold,c=[shade(model,r.size,q['size'].tolist())],s=30,edgecolors='black',linewidths=.4);ax.annotate(f'C={r.size:g}',(r.map_amplitude_fold,r.half_turn_time_fold),xytext=(3,3),textcoords='offset points',fontsize=7)
    ax.set_xlabel('Map-amplitude fold');ax.set_yscale('log')
    if j==0:ax.set_ylabel('Separated-cue half-turn-time fold')
   fig.text(.07,.05,'Map amplitude = max |standard pulse response|. Plateau = median of furthest three observed two-region half-turns.\nExact size/source gaps and plateau ranges remain in matched_B_susceptibility.csv; this definition is explicit and does not claim a fitted asymptote.',fontsize=9)
  page('S8. Matched B sizes link separated-cue turning to weak susceptibility',matched)
  def absolute(fig):
   grid=facets(fig,bottom=.17,top=.82);frame=pd.read_csv(DATA/'edge_prediction_comparison.csv')
   for j,model in enumerate(MODELS):
    ax=fig.add_subplot(grid[0,j]);ax.set_title(NAMES[model]);q=frame[(frame.model==model)&frame.defined&(frame.observed_Ds>0)]
    lo=min(q.observed_Ds.min(),q.Ds_predicted_independent.min())*.5;hi=max(q.observed_Ds.max(),q.Ds_predicted_full_covariance.max())*2;ax.plot([lo,hi],[lo,hi],color='.5',lw=.8)
    for r in q.itertuples():
     col=shade(model,r.size);ax.plot(r.observed_Ds,2*r.raw_native_Ds_independent,'o',ms=4,color=col);ax.plot(r.observed_Ds,r.raw_native_Ds_independent,'o',ms=4,mfc='white',mec=col);ax.plot(r.observed_Ds,r.Ds_predicted_full_covariance,'+',ms=4,color=col,mew=.5)
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(lo,hi);ax.set_ylim(lo,hi);ax.set_xlabel(r'Measured absolute $D_s$')
    if j==0:ax.set_ylabel(r'Predicted absolute $D_s$')
   fig.text(.07,.05,'Filled: raw native shared translation. Open: raw native independent-edge candidate, rejected by co-motion. Plus: historical pulse-normalized shared prediction. Each model has its own scale.',fontsize=9)
  page('S9. Absolute native diffusion, historical pulse-normalization and rejected independent edges',absolute)
  def integration(fig):
   grid=facets(fig,bottom=.15,top=.83)
   for j,model in enumerate(MODELS):
    ax=fig.add_subplot(grid[0,j]);ax.set_title(NAMES[model]);q=ex[(ex.model==model)&(ex.kind=='map_validation')&(ex.label=='right_edge')].sort_values('size');ax.plot(q['size'],q.integrated_shared_gain,'o-',color=COLORS[model]);ax.set_xlabel('Circumference parameter')
    if j==0:ax.set_ylabel(r'Integrated displacement / dose, through $20\tau$')
  page('S10. Finite integrated response retained separately from terminal susceptibility',integration)
  page('S11. Weak-input alpha bridge: qualified observations and explicit missingness',lambda fig:alpha_panel(fig,facets(fig,bottom=.17,top=.83),br))
 # The entire preceding supplement is kept after new analysis pages.
 writer=PdfWriter();writer.append(str(new_pdf));historical=OLD/'figure3_round5_review_supplement.pdf';writer.append(str(historical));writer.write(str(OUT/f'{STEM}_supplement.pdf'))
 write(META/'supplement_preservation.json',dict(generated_utc=new.old.now(),new_pages=len(PdfReader(new_pdf).pages),historical_path=rel(historical),historical_sha256=new.old.sha(historical),historical_pages=len(PdfReader(historical).pages),total_pages=len(PdfReader(OUT/f'{STEM}_supplement.pdf').pages)))

def documents(mp,ex,cv,fr,li,br,ma):
 folds=pd.read_csv(DATA/'edge_response_size_folds.csv');normal=pd.read_csv(DATA/'raw_native_size_normalized_prediction.csv');wp=normal[(normal.model=='wave_pinning')&(normal['size']==60)&(normal.fraction==.001)].iloc[0]
 maxerr=float(ex[ex.kind=='map_validation'].normalized_map_error.max());limits=[]
 for model in MODELS:
  q=li[li.model==model].sort_values('size');a=q.iloc[0];b=q.iloc[-1];limits.append(f'{NAMES[model]}: first >10% departures {a.first_stopping_amplitude:g} at C={a["size"]:g} and {b.first_stopping_amplitude:g} at C={b["size"]:g}, size fold {b.first_stopping_amplitude/a.first_stopping_amplitude:g}; last tested linear amplitudes {a.last_linear_amplitude:g} and {b.last_linear_amplitude:g}')
 fraction_text=[]
 for model in MODELS:
  for size in endpoints(model):
   r=fr[(fr.model==model)&(fr['size']==size)].iloc[0];fraction_text.append(f'{NAMES[model]} C{size:g}: left edge {r.left_edge_fraction:.3f}, right edge {r.right_edge_fraction:.3f}, plateau interior {r.plateau_interior_fraction:.3f}, basal arc {r.basal_arc_fraction:.3f}')
 terminal=', '.join(f'{NAMES[r.model]} {r.terminal_gain_large_over_small:.3f}' for r in folds.itertuples())
 title='Size-dependent susceptibility links deterministic perturbations to weak-input domain drift'
 legend=[
 'A. Exact preceding Wave-pinning activity profiles, C=30, for a cue straddling the incumbent edge (i) or separated from it (ii). Rows are initial, half-turn, 90% turn and cue-free relaxation; the grey band is cue support, diamonds are saddles and the black tick is heading. All A data and timings are unchanged.',
 'B. Exact preceding half-turn time against cue-edge gap for the three matched sizes, light to dark: C=30,37.8,47.6; Spring 2.5, 3.15, 3.97. Lines mark first saddle appearance; grey marks undefined heading at one or more sizes. Open upward triangles and arrows denote lower bounds. Continuous-cue data and native parameters are unchanged.',
 'C. A weak input anywhere shifts the domain by g(x) times dose; noise is many such inputs; drift follows from g and the input covariance; nothing is fitted. Lines are the terminal displacement per dose predicted by the linearized native model for the standard Gaussian pulse centred at each arc position: SD 1.5, duration 50 (Spring 0.15, 5), at f=0.001. Light and dark are C20/C60 (Spring 1/3). New nominal-amplitude tests from settled states are filled circles; the existing edge pulse is an open square. Both subtract a matched unpulsed control. The pulse confirms g at the edge; the map elsewhere comes from the linearized native model and is tested at domain centre, halfway to the edge, edge±1/3 pulse SD, quarter-basal arc and antipode. Maximum absolute gain error normalized to map amplitude is '+f'{100*maxerr:.2f}%. Terminal edge gain large/small: {terminal}.',
 'D. Mean-state activity in light grey on its own axis, sharing the centred physical arc coordinate with C and E. Dotted verticals locate half-contrast edges. Spring positions are material fractions times the mean-state circumference.',
 'E. Signed per-position contribution c_j/D_s at the same sizes and input, using raw native g and the unchanged noise covariance, with no measured-gain rescaling. Contributions sum to one; negative values are retained. Edge zones are within two input correlation lengths of the nearest edge; plateau interior and basal arc exclude these zones. Fractions of D_s, left/right edge, plateau interior, basal arc respectively: '+'; '.join(fraction_text)+'. Raw-native and measurement-normalized totals are both exported.',
 'F. Wave-pinning C20 and C60: pulsed falling edge (solid) and opposite rising edge (dashed) move together and retain a permanent offset after the grey pulse-on interval. Displacement is divided by pulse dose. Width and centroid histories, other models, every amplitude and duration are supplementary.',
 'G. Predicted against measured physical domain diffusion after each is divided by its own fixed-smallest-size value at the same input. Only the shared-translation prediction appears here; the grey line is equality. The main prediction uses raw native g without measured-gain rescaling; the historical pulse-normalized variant is retained in the supplement. Shape encodes f=0.0001,0.0003,0.001,0.003,0.005,0.01 (circle,square,up/down triangle,diamond,plus); size retains the model colour ramp. Horizontal intervals are measured bootstrap 95% intervals. The ratio is approximately the squared edge-gain fold: Wave-pinning 0.70²≈0.49; the full spatial covariance gives '+f'{wp.shared_translation_predicted_ratio:.3f}, versus {wp.observed_ratio:.3f} measured at f=0.001. Shape changes make the square relation approximate, not an exact identity. Absolute predictions and independent-edge candidate are supplementary. Goryachev stays on its own scale and fails quantitatively.',
 'H. Local α=2−dlnD_s/dlnC implied by the predicted drift against Figure 4 weak-input local α (f≤0.001), with identical fixed-neighbour size fits, two neighbours at endpoints and three inside. Grey is equality. Vertical intervals are Figure 4 bootstrap 95% intervals. '+f'{br.observed_alpha.notna().sum()} of {len(br)} requested pairs qualify; Wave-pinning and Goryachev have no qualified weak-input observed α under the fixed shared-lag rule. No missing estimate is replaced. Marker shapes distinguish the three weak inputs as in G. '
 ]
 readme='''# Unexpected

'''+f'''The requested map tests support the native linear response away from the edge: the maximum error is {100*maxerr:.2f}% of the condition's map amplitude. Symmetric centre/antipode pulses give essentially zero shared translation, while off-centre pulses shift the complete domain. This does not repair Goryachev's absolute drift disagreement.

The covariance contribution maps localize 99.59–100.00% of predicted drift within two input correlation lengths of the two edges at the endpoint sizes. Small signed plateau/basal contributions and all negative bins are retained.

The requested weak-input alpha bridge is only partially observable: {br.observed_alpha.notna().sum()}/{len(br)} rows have a qualified Figure 4 exponent. Wave-pinning and Goryachev have none under the existing shared-lag rule. Their deterministic alpha predictions remain in the CSV; they are not plotted against invented observations.

The matched-size comparison is now explicit. The operational separated-cue plateau is the median of the furthest three finite observed two-region half-turns, with every selected gap and range exported. It is not a fitted asymptote.

# Scientific results

Terminal largest/smallest edge-gain folds are {terminal}. Wave-pinning's f=0.001 size-normalized drift ratio remains {wp.shared_translation_predicted_ratio:.6f} predicted and {wp.observed_ratio:.6f} measured. The full spatial prediction is approximately the squared terminal edge-gain fold; it is not mathematically identical when the susceptibility shape and covariance change with size.

The linear range is bracketed by the last tested amplitude within 10% of nominal gain and the first stopping amplitude. All 20 native-size sweeps stop on gain departure, before the continuous-cue peak and without losing the unique half-contrast domain:

'''+ '\n'.join('- '+x for x in limits)+r'''

Thresholds are geometric-grid brackets, not precisely located bifurcations. Fixed-dose durations test 12.5, 25, 50, 100, 200 chemical units and Spring/10 at the unchanged nominal dose.

# Methods and provenance

Starting reference: `manuscript_outline/outline_v11.pdf` and `.json`. A–B source coordinate files are byte-identical to round 5 and all previous editions remain intact. C–H use deterministic assays and existing observed drift/heading tables. No stochastic realization, kinetic parameter fit or retuning was introduced.

The 28 conditions include all 20 original f=0.001 sizes and 8 added B-matched sizes. Original states, adjoints and covariance kernels are reused from the frozen round 5 campaign. Added states periodically resample the nearest original basal onto unchanged native cell density and settle at the exact mean input for 16000 chemical or 1600 Spring units. Eligibility uses the unchanged unique-region, contrast≥0.05, native RHS≤10⁻⁶ and tail-change≤10⁻⁴ gates. All 28 pass. BDF uses rtol 2×10⁻¹⁰, atol 2×10⁻¹², exact pulse withdrawal restart, and the unchanged native sparse Jacobian. Every integration ran under SLURM 21334289 on one CPU, small QoS, 4 GB, 30-minute ceiling; actual workers took under 36 seconds.

For the full native state y, translation mode q, and normalized adjoint p:

$$
Jq=0,\qquad J^{\mathsf T}p=0,\qquad p^{\mathsf T}q=1,\qquad g_j=(p^{\mathsf T}B)_j. \tag{1}
$$

The plotted standard-pulse map is the convolution of raw native g with the exact periodic Gaussian pulse shape h, per peak-amplitude×duration dose:

$$
G(x_0)=\sum_j g_j h(x_j-x_0),\qquad h(d)=\exp[-d^2/(2\ell^2)]. \tag{2}
$$

Raw cell-input g appears as thin lines in the supplement; it is not substituted for G on the main map. Points measure the terminal common shift of the two half-contrast edges divided by pulse dose, subtracting matched unpulsed controls. Every test retains left/right edges, activity centroid, width and region count through 20 nominal pulse durations. Spring uses material fractions multiplied by its mean circumference, excluding uniform extension.

Main C/E/G/H use raw native g without any empirical amplitude normalization. The historical round 5 prediction used the mean of the two independently measured edge/native gain factors s, close to one; that variant is exported separately and retained in the historical supplement. Observed stochastic drift is never used. The raw-native contribution partition uses

$$
c_j=g_j\sum_k g_k\int_0^\infty\overline{K_u(j,k;t,t+v)}^{\,t}\,\mathrm{d}v,\qquad D_s=\sum_j c_j. \tag{3}
$$

The input covariance is the exact native lognormal image of the discrete spatial Gaussian filter with variance-normalized OU-knot interpolation, zero spatial mode retained. Round5's deterministic covariance quadrature is reused. Edge zones use two input correlation lengths; where zones overlap the nearest edge assigns the cell. Plateau and basal remainders exclude edge zones, so the four signed fractions sum to one. No negative contribution is clipped. Every map bin and contribution is in `supporting/data/susceptibility_and_contributions.csv`.

At the right edge, amplitudes double from nominal/2 until gain departs by >10%, topology changes, or the native continuous-cue peak is reached. The first stopping point is retained. Duration tests hold nominal aτ constant and retain the same final 20τ nominal observation horizon. No response threshold or amplitude was tuned to improve agreement.

$$
\alpha_{\mathrm{pred}}=2-\frac{\mathrm{d}\ln D_s}{\mathrm{d}\ln \bar C}. \tag{4}
$$

The local derivative uses the exact fixed adjacent-size neighborhoods of Figure 4 (two members at grid endpoints, three at interior sizes) and mean-state physical circumference. The observed exponent retains its existing coverage and shared-lag eligibility rules; unavailable estimates remain rows with their reason.

Length denotes WGD-associated enlargement. Observed migration summaries, conditional OU persistence estimates and simulated 1D polarity-heading readouts remain distinct. Domain drift is not identified with total heading decay or empirical migration persistence.

# Outputs and rebuild

Main: `figure3_round6.pdf`, `.png`, `.svg`. Complete supplement appends the unchanged round 5 reviewed supplement after 11 new pages. Compact CSVs include every plotted source point, raw/smoothed maps, signed contribution bins, pulse tests/time courses, amplitude brackets, matched-size comparison and alpha bridge. `supporting/metadata/provenance.json` records hashes and exact sources.

Run `./rebuild.sh` to verify/collect/render from completed workers. Run `./rebuild.sh --native` only to prepare/submit any missing deterministic workers and then render; it requires SLURM and does not repeat completed conditions. The output tree contains only small derivative/native scalar archives; no large full-field campaign was created.
'''
 # Backslashes in Python text are kept literal; no raw Markdown equations in chat.
 (OUT/'scientific_readme.md').write_text(readme);(OUT/'METHODS_RESULTS.md').write_text(readme)
 (OUT/'legend.md').write_text('# '+title+'\n\n'+'\n\n'.join(legend)+'\n\nLinear range: '+'; '.join(limits)+'. Full response curves, duration tests, matched-size comparison, raw g and absolute drift are supplementary.\n')
 narrative=[f'Weak standard-pulse susceptibility changes with size in a model-dependent manner; its spatial map is validated at edge, active-interior and basal-arc positions (maximum error {100*maxerr:.2f}% of map amplitude). The covariance-weighted map predicts domain drift without fitting stochastic data.',f'Wave-pinning f=0.001 C60/C20 drift remains {wp.shared_translation_predicted_ratio:.3f} predicted and {wp.observed_ratio:.3f} measured; Goryachev remains quantitatively inconsistent. The weak-input alpha bridge retains unavailable observed exponents explicitly.','Linear range at the right edge (first tested >10% gain departure, with last linear amplitude one doubling below): '+'; '.join(limits)+'.']
 section=dict(id='figure3',title=title,status='round6_complete',legend=legend,role='Spatial susceptibility and signed covariance contributions expose the deterministic input-to-drift bridge.',narrative=narrative,measurements=['Map terminal domain shift per dose with local weak pulses at edge, interior and basal positions.','Measure both boundaries, centroid and width; compare deterministic susceptibility and covariance-based drift at matched sizes.'],limitations=['Goryachev absolute drift prediction fails despite small deterministic-map validation error.','Existing weak-input alpha estimates are unavailable for Wave-pinning and Goryachev under the fixed shared-lag rule.','Linear-limit amplitudes are geometric-grid brackets; separated-cue plateau is an operational finite-data median.','Deterministic drift, total heading decay and empirical migration persistence remain distinct.'],equations=[dict(math=r'G(x_0)=\sum_j g_j h(x_j-x_0)',number='3.1',description='The main map is the standard Gaussian-pulse response per dose, not raw cell-input g.'),dict(math=r'c_j=g_j\sum_k g_k\int_0^\infty\overline{K_u(j,k;t,t+v)}^{\,t}\,\mathrm{d}v,\qquad D_s=\sum_j c_j',number='3.2',description='Signed per-position contributions retain the exact native input covariance.'),dict(math=r'\alpha_{\mathrm{pred}}=2-\frac{\mathrm{d}\ln D_s}{\mathrm{d}\ln \bar C}',number='3.3',description='Same fixed adjacent-size neighborhoods as Figure 4; missing observed alpha remains unavailable.')],assets=[dict(path=rel(OUT/f'{STEM}.png'),native_pdf_page=rel(OUT/f'{STEM}.pdf'),page_size=[16.2,23],minimum_printed_width_inches=16.2,page_title='')],main_pdf=rel(OUT/f'{STEM}.pdf'),supplement_pdf=rel(OUT/f'{STEM}_supplement.pdf'),reader=rel(OUT/'index.html'),sources=[dict(label='Scientific methods and unexpected results',path=rel(OUT/'scientific_readme.md')),dict(label='Susceptibility maps and signed contributions',path=rel(DATA/'susceptibility_and_contributions.csv')),dict(label='Every deterministic assay',path=rel(DATA/'deterministic_experiments.csv')),dict(label='Alpha bridge and unavailable rows',path=rel(DATA/'alpha_bridge.csv'))])
 write(OUT/'legend.json',section)
 timestamp=new.old.now();(OUT/'README.md').write_text('# Figure 3 round6\n\n[Main PDF](figure3_round6.pdf) · [complete supplement](figure3_round6_supplement.pdf) · [legend](legend.md) · [scientific readme](scientific_readme.md).\n\nAll 28 deterministic workers complete. One rebuild entry point: `./rebuild.sh`.\n')
 (OUT/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Figure 3 round6</title><style>body{max-width:1300px;margin:2rem auto;font:17px/1.5 system-ui;color:#17212b}img{max-width:100%}</style><h1>'+html.escape(title)+'</h1><p>Data generated '+timestamp+'; all native producers complete.</p><p><a href="figure3_round6.pdf">Main PDF</a> · <a href="figure3_round6_supplement.pdf">Complete supplement</a> · <a href="scientific_readme.md">Scientific methods and results</a></p><img src="figure3_round6.png" alt="Figure 3 susceptibility, signed contributions and drift bridge">'+''.join('<p>'+html.escape(p)+'</p>' for p in legend))
 checks=dict(generated_utc=timestamp,complete=True,workers=28,mean_states_eligible=28,maximum_map_amplitude_normalized_error=maxerr,contribution_partition_max_abs_error=float(abs(fr.partition_sum-1).max()),weak_alpha_pairs=int(br.observed_alpha.notna().sum()),weak_alpha_requested=len(br),AB_data_identical=all(new.old.sha(DATA/name)==new.old.sha(OLD_DATA/name) for name in ('activity_profiles.csv.gz','half_heading_readouts.csv','regime_lines.csv','undefined_heading_bands.csv')),main_pages=len(PdfReader(OUT/f'{STEM}.pdf').pages),legend_sha256=new.old.sha(OUT/'legend.json'),assembler_sha256=new.old.sha(Path(__file__)),worker_source_matches_manifest=new.old.sha(HERE/'native_analysis.py')==new.read(META/'native_manifest.json')['source_code_sha256'][rel(HERE/'native_analysis.py')],main_pdf_sha256=new.old.sha(OUT/f'{STEM}.pdf'),main_png_sha256=new.old.sha(OUT/f'{STEM}.png'),allocated_bytes=new.old.source.allocated(OUT),native_stochastic_runs=0,retuned_parameters=0)
 assert checks['AB_data_identical'] and checks['contribution_partition_max_abs_error']<1e-12;write(META/'verification.json',checks)

def main():
 values=collect();main_figure(*values);supplement(*values);documents(*values);print(OUT/f'{STEM}.pdf')
if __name__=='__main__':main()
