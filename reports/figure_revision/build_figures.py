"""Build revised standalone figures only. Never import a report renderer."""
from pathlib import Path
import sys,json,csv,gzip,importlib.util,hashlib
from collections import defaultdict
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;OUT=HERE/'figures';OUT.mkdir(exist_ok=True)
DATA=ROOT/'polarity1d/results/data';OLD=ROOT/'reports/trial_report_assets'
NAVY='#19394d';TEAL='#078b8a';ORANGE='#d87529';COLORS=['#169b91','#526f93','#83aa45','#d9ad3a'];MODELS=['wave_pinning','otsuji','goryachev'];NAMES=['Wave-pinning','Otsuji','Goryachev'];LENGTHS=[9,12,15,18,21,27,33,39]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
EXPORT={}
def rows(p):return list(csv.DictReader(open(p)))
def module(path):
 spec=importlib.util.spec_from_file_location(path.stem,path);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m

def save(fig,name):
 fig.savefig(OUT/(name+'.png'),dpi=240,bbox_inches='tight');fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def title(ax,letter,text):
 ax.set_title('',loc='center');ax.set_title(text,loc='left',pad=10,color=NAVY);ax.text(-.14,1.07,letter,transform=ax.transAxes,fontsize=15,fontweight='bold',color=NAVY)

def poolrows():
 return [r for p in (HERE/'data').glob('pool_*.json') for r in json.loads(p.read_text())]

def trajectories2():
 out=defaultdict(list)
 with gzip.open(DATA/'figure2/figure2_v2_polarity_trajectories.csv.gz','rt') as f:
  for r in csv.DictReader(f):
   if (r['panel']=='A' and r['variant']=='full') or (r['panel']=='C' and r['field']=='rac') or r['panel']=='E':
    out[r['panel'],r['variant'],r['task_id'],float(r['length_um']),r['field']].append((float(r['display_time_s']),float(r['x50_over_L'])))
 return {k:np.array(v) for k,v in out.items()}
def field_kymo(ax,t,x,z,label,cmap='viridis',relative=True):
 stride=max(1,len(t)//1400);idx=np.unique(np.r_[np.arange(0,len(t),stride),len(t)-1]);v=z[idx]
 im=ax.pcolormesh(t[idx],x,v.T,shading='nearest',cmap=cmap,rasterized=True)
 ax.set_ylim(x[0],x[-1]);ax.set_xlabel('Time from cue onset (s)');ax.set_ylabel('Position / L' if relative else 'Position');ax.set_title(label,loc='left',fontsize=9)
 return im

def figure2():
 tr=trajectories2();steady=pd.read_csv(DATA/'figure2/figure2_v2_steady_state_metrics.csv');threshold=pd.read_csv(DATA/'figure2/figure2_v2_minimum_countercue_duration.csv')
 response=json.loads((ROOT/'polarity1d/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_linear50_response.json').read_text());taskmap={t['task_id']:t for t in response['tasks']}
 fig=plt.figure(figsize=(16,15));gs=fig.add_gridspec(3,3,left=.06,right=.93,bottom=.065,top=.95,wspace=.35,hspace=.38,height_ratios=[1.8,1,1.4])
 # A retains original full-model Rac/MCA trajectories and adds full-horizon fields at L=21.
 ga=gs[0,:2].subgridspec(3,2,hspace=.70,wspace=.30)
 camp=ROOT/'polarity1d/results/campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000';r=next(r for r in rows(camp/'task_summary.csv') if r['coupling_variant']=='full' and float(r['length_um'])==21)
 z=np.load(camp/'tasks'/r['task_id']/'trajectory.npz');t=z['time']-200;x=z['x']/21
 norm=matplotlib.colors.Normalize(1,39);cmap=plt.get_cmap('viridis')
 for j,(field,label,cm) in enumerate([('rac','Rac','viridis'),('mca_bound','Bound MCA','cividis')]):
  ax=fig.add_subplot(ga[0,j]);im=field_kymo(ax,t,x,z['field__'+field],f'L = 21 µm: {label}',cm)
  ax.set_xscale('symlog',linthresh=200);ax.set_xlim(-200,9800)
  for boundary in [0,160]:ax.axvline(boundary,c='white',ls='--',lw=.8)
  fig.colorbar(im,ax=ax,pad=.02,fraction=.025).set_label(label)
  a=next(v for k,v in tr.items() if k[0]=='A' and k[3]==21 and k[4]==field)
  ax.plot(a[:,0],a[:,1],c='white',lw=1.4);ax.plot(a[:,0],a[:,1],c='black',lw=.6)
  if j==0:title(ax,'A','Full-model fields and polarity')
  ax=fig.add_subplot(ga[2,j])
  for k,a in tr.items():
   if k[0]=='A' and k[4]==field:
    stride=max(1,len(a)//1000);ax.plot(a[::stride,0],a[::stride,1],c=cmap(norm(k[3])),lw=.8)
    st=steady[(steady.panel=='A')&(steady.variant=='full')&(steady.field==field)&(steady.length_um==k[3])]
    if len(st) and bool(st.iloc[0].polarized_steady_state):
     tt=st.iloc[0].steady_state_time_s-200;ax.scatter(tt,np.interp(tt,a[:,0],a[:,1]),s=12,c=[cmap(norm(k[3]))])
  ax.set_xscale('symlog',linthresh=100);ax.axvspan(0,160,color='#edb63e',alpha=.2);ax.axhline(.5,c='.6',ls=':',lw=.8);ax.set(xlim=(-200,10000),xlabel='Time from cue onset (s)',ylabel=f'{label} half-mass position / L');ax.grid(alpha=.1)
 for j,(field,label,cm) in enumerate([('rho','Rho','magma'),('tension','Tension','coolwarm')]):
  ax=fig.add_subplot(ga[1,j]);im=field_kymo(ax,t,x,z['field__'+field],f'L = 21 µm: {label}',cm)
  ax.set_xscale('symlog',linthresh=200);ax.set_xlim(-200,9800)
  for boundary in [0,160]:ax.axvline(boundary,c='white',ls='--',lw=.8)
  fig.colorbar(im,ax=ax,pad=.02,fraction=.025).set_label(label)
 ax=fig.add_subplot(gs[0,2]);a=steady[(steady.panel=='A')&(steady.variant=='full')&(steady.field=='rac')&(steady.polarized_steady_state==True)&(steady.steady_state_reached==True)].sort_values('length_um')
 ax.plot(a.length_um,a.steady_state_time_s-200,'o-',c=NAVY,ms=4);ax.set(xlabel='Domain length L (µm)',ylabel='Rac settling time from cue onset (s)');title(ax,'B','Full-model polarization time');ax.grid(alpha=.15)
 gc=gs[1,:2].subgridspec(1,1,wspace=.24)
 feedback_colors=['#d87529','#654195','#169b69']
 for i in [1]:
  ax=fig.add_subplot(gc[0,0])
  for k,a in tr.items():
   if k[0]=='C' and k[1]==f'beta_{i}x':
    stride=max(1,len(a)//1000);ax.plot(a[::stride,0],a[::stride,1],c=cmap(norm(k[3])),lw=.8)
    st=steady[(steady.panel=='C')&(steady.variant==k[1])&(steady.field=='rac')&(steady.length_um==k[3])].iloc[0]
    if st.steady_state_reached:
     tt=st.steady_state_time_s;ax.scatter(tt,np.interp(tt,a[:,0],a[:,1]),s=13,marker='o' if st.polarized_steady_state else 'x',c=[cmap(norm(k[3]))])
  ax.axvspan(0,80,color='#edb63e',alpha=.2);ax.axhline(.5,c='.6',ls=':',lw=.8);ax.set_xscale('symlog',linthresh=100);ax.set(xlim=(0,10000),xlabel='Time from initial cue (s)',title=f'Rho feedback {i}×');ax.grid(alpha=.1)
  ax.set_ylabel('Rac half-mass position / L');title(ax,'C','Full model: default Rho feedback, 80-s initial cue')
  ax.legend(handles=[Line2D([0],[0],marker='o',ls='',c='.3',label='Polarized steady state'),Line2D([0],[0],marker='x',ls='',c='.3',label='Unpolarized steady state')],frameon=False,fontsize=8,loc='lower left')
 ax=fig.add_subplot(gs[1,2])
 for i,c in enumerate(feedback_colors):
  a=steady[(steady.panel=='C')&(steady.variant==f'beta_{i}x')&(steady.field=='rac')&(steady.polarized_steady_state==True)&(steady.steady_state_reached==True)].sort_values('length_um')
  ax.plot(a.length_um,a.steady_state_time_s,'o-',c=c,ms=4,label=f'Rho feedback {i}×')
 ax.set(xlabel='Domain length L (µm)',ylabel='Rac settling time (s)');title(ax,'D','Settling of polarized states');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
 ge=gs[2,:2].subgridspec(3,2,wspace=.25,hspace=.38)
 for i,duration in enumerate([150,300,600]):
  for j,(field,label) in enumerate([('rac','Rac'),('mca_bound','MCA')]):
   ax=fig.add_subplot(ge[i,j])
   for k,a in tr.items():
    if k[0]=='E' and k[4]==field and taskmap[k[2]]['countercue_duration_s']==duration:
     stride=max(1,len(a)//1000);ax.plot(a[::stride,0],a[::stride,1],c=cmap(norm(k[3])),lw=.8)
   ax.set_xscale('symlog',linthresh=300);ax.axvspan(0,duration,color='#edb63e',alpha=.23);ax.axhline(.5,c='.6',ls=':',lw=.8);ax.set(xlim=(-1000,10000),ylim=(.15,.85),ylabel=f'{label} position / L');ax.grid(alpha=.1);ax.text(.02,.08,f'{duration}-s opposing cue',transform=ax.transAxes,fontsize=8)
   if i==2:ax.set_xlabel('Time from opposing cue onset (s)')
   else:ax.set_xticklabels([])
   if i==0:title(ax,'E' if j==0 else '',f'Established polarity: {label}')
 ax=fig.add_subplot(gs[2,2])
 for i,c in enumerate(feedback_colors):
  a=threshold[(threshold.field=='rac')&(threshold.variant==f'beta_{i}x')].set_index('length_um');xx=np.arange(1,40,2);yy=[a.loc[L,'minimum_successful_countercue_duration_s'] if L in a.index else np.nan for L in xx]
  ax.plot(xx,yy,'o-',c=c,ms=4,label=f'Rho feedback {i}×')
 ax.axvspan(1,17,color='.9');ax.text(8.5,530,'No eligible\ninitial polarity',ha='center',fontsize=9,color='.4');ax.set(xlim=(1,40),ylim=(170,850),xlabel='Domain length L (µm)',ylabel='Shortest successful opposing cue (s)');title(ax,'F','Lasting Rac reversal');ax.legend(frameon=False,fontsize=8,loc='upper left');ax.grid(alpha=.15)
 cax=fig.add_axes([.96,.24,.012,.50]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,label='Domain length L (µm)')
 fig.suptitle('Mechanochemical polarity: establishment, coordination and lasting reversal',weight='bold',fontsize=16,color=NAVY,y=.995);save(fig,'figure2')
 EXPORT['figure2']={'representative_length_um':21,'kymograph_source':str(camp/'tasks'/r['task_id']/'trajectory.npz'),'trajectory_source':str(DATA/'figure2/figure2_v2_polarity_trajectories.csv.gz')}

def figure3():
 from build_figure3_updated import build_figure3
 EXPORT['figure3']=build_figure3()

def figure1():
 m=module(ROOT/'polarity1d/scripts/build_wp1_manuscript_figures.py')
 fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(3,3,left=.065,right=.975,bottom=.07,top=.93,wspace=.35,hspace=.45)
 pool=poolrows();fixed=[r for p in (HERE/'data').glob('fixed_*.json') for r in json.loads(p.read_text())];primary=[r for r in pool if r.get('n')==201 and r.get('status')!='no_matching_root']
 for row,assay,origin,times,letters in [(0,'competition','absolute',(20.,200.,20000.),'ABC'),(1,'reversal','reversal',(0.,20.,100.),'DEF')]:
  ax=fig.add_subplot(gs[row,0]);m._cue_kymograph(ax,assay);m._panel(ax,letters[0]);
  sub=gs[row,1].subgridspec(2,2,hspace=.2,wspace=.15);axes=[fig.add_subplot(sub[i,j]) for i in range(2) for j in range(2)]
  m._focused_profile_snapshots(axes,assay,times,OUT/f'figure1_{assay}_profiles.csv',origin=origin,zero_from_left=(assay=='reversal'));title(axes[0],letters[1],'');axes[1].get_legend().remove();axes[1].legend(frameon=False,fontsize=6,loc='upper left',bbox_to_anchor=(1.02,1.05),title='Time')
  ax=fig.add_subplot(gs[row,2]);key='competition_resolution_latency_s' if assay=='competition' else 'reversal_latency_s';prefix='competing' if row==0 else 'gradient'
  for model,c,label in zip(MODELS,COLORS,NAMES):
   for arm,ls in [('concentration','-'),('fixed_peak','--')]:
    xx=[];yy=[]
    for L in LENGTHS:
     rr=[r for r in (fixed if arm=='fixed_peak' else primary) if r.get('n')==201 and r['model']==model and r['L']==L and r['arm']==arm and r['assay'].startswith(prefix)]
     if arm=='concentration' and L>=15:rr=[r for r in fixed if r['n']==201 and r['model']==model and r['L']==L and r['arm']=='native_restart' and r['assay'].startswith(prefix)]
     key_used='competition_latency_from_original_withdrawal_s' if row==0 and L>=15 else key
     vals=[r[key_used] for r in rr if r.get(key_used) is not None];xx.append(L);yy.append(np.median(vals) if vals else np.nan)
     if rr and not vals:
      eligible=any(r.get('competition_established' if row==0 else 'reversal_conditioned',False) for r in rr)
      ax.scatter([L],[24000 if eligible else 7],marker='^' if eligible else 'x',c=c,s=20)
    ax.plot(xx,yy,ls=ls,marker='o' if arm=='concentration' else 's',ms=3,c=c,label=label if arm=='concentration' else None)
  # LEGI remains a cue-following comparator in persistent gradient reversal.
  if row==1:
   rr=rows(DATA/'manuscript/figure1_long_reversal.csv');a=[r for r in rr if r['model']=='jilkine_legi' and r['reversal_status']=='success'];ax.plot([float(r['length_um']) for r in a],[float(r['reversal_latency_s']) for r in a],'o-',c=COLORS[3],ms=3,label='LEGI')
  ax.set(yscale='log',ylim=(5,33000),xlabel='Domain length L (µm)',ylabel='Competition resolution latency (s)' if row==0 else 'Reversal latency (s)');ax.axvline(15,c='.5',ls=':',lw=.8);ax.grid(alpha=.12);title(ax,letters[2],'Competition' if row==0 else 'Persistent-gradient reversal');ax.legend(frameon=False,fontsize=7,loc='upper left')
  ax.text(.02,.02,'▲ horizon-censored; × not established',transform=ax.transAxes,fontsize=6.5)
 # G finite opposing cue schematic, H example traces, I duration brackets.
 ax=fig.add_subplot(gs[2,0]);t=np.r_[np.linspace(-200,0,80,endpoint=False),np.linspace(0,320,160),np.linspace(321,1600,160)];x=np.linspace(0,1,100);cue=np.where((t>=0)&(t<=320),1.,0.)[:,None]*x[None,:];im=ax.pcolormesh(t,x,cue.T,cmap='YlOrBr',shading='nearest',rasterized=True);ax.set(xlabel='Time from opposing cue (s)',ylabel='Position / L');ax.axvline(320,c='black',ls='--',lw=.8);title(ax,'G','Finite cue, then stimulus-free recovery');ax.text(.04,.93,'Conditioned, autonomous polarity',transform=ax.transAxes,fontsize=8)
 ax=fig.add_subplot(gs[2,1]);example=[]
 for model,c,label in zip(MODELS,COLORS,NAMES):
  p=HERE/f'data/finite_{model}_L15_n201.json';r=json.loads(p.read_text());assays=r['assays']
  if not assays:continue
  yes=[a for a in assays if a['success']];d=min(yes,key=lambda a:a['duration'])['duration'] if yes else assays[-1]['duration'];z=np.load(HERE/f'data/finite_{model}_L15_n201_d{d:g}.npz')
  ax.plot(z['time'],z['polarity'],c=c,label=f'{label}: {d:g} s');ax.axvline(d,c=c,ls=':',lw=.7);example.append(dict(model=model,L=15,duration=d))
 ax.set_xscale('symlog',linthresh=100);ax.axhline(0,c='.5',ls='--',lw=.8);ax.set(xlim=(-500,10000),ylim=(-1,1),xlabel='Time from opposing cue (s)',ylabel='Active-GTPase first moment / total');title(ax,'H','Replacement must persist after withdrawal');ax.legend(frameon=False,fontsize=7,loc='lower right');ax.grid(alpha=.12)
 ax=fig.add_subplot(gs[2,2]);summ=[]
 for model,c,label in zip(MODELS,COLORS,NAMES):
  xx=[];yy=[]
  for L in LENGTHS:
   r=json.loads((HERE/f'data/finite_{model}_L{L}_n201.json').read_text());yes=[a['duration'] for a in r['assays'] if a['success']]
   if not r['eligible']:ax.scatter([L],[6],marker='x',c=c,s=20);continue
   if not yes:ax.scatter([L],[6500],marker='^',c=c,s=20);continue
   hi=min(yes);no=[a['duration'] for a in r['assays'] if not a['success'] and a['duration']<hi];lo=max(no) if no else 0;xx.append(L);yy.append(hi);ax.plot([L,L],[max(lo,5),hi],c=c,lw=2);summ.append(dict(model=model,L=L,lower=lo,upper=hi))
  ax.plot(xx,yy,'o-',c=c,ms=3,label=label)
 ax.set(yscale='log',ylim=(5,9000),xlabel='Domain length L (µm)',ylabel='Shortest tested lasting-reversal cue (s)');title(ax,'I','Resistance to finite opposing cues');ax.text(.02,.02,'▲ no lasting reversal by 5,120-s cue; × no eligible state',transform=ax.transAxes,fontsize=6.3);ax.grid(alpha=.12);ax.legend(frameon=False,fontsize=7)
 fig.legend(handles=[Line2D([0],[0],c='.3',ls='-',label='Fixed total concentration'),Line2D([0],[0],c='.3',ls='--',marker='s',ms=3,label='Same physical peaks + background pool in added domain')],loc='upper center',bbox_to_anchor=(.55,.985),ncol=2,frameon=False)
 save(fig,'figure1');EXPORT['figure1']={'finite_examples':example,'finite_thresholds':summ,'control':'Physical peak preservation at the restart: competition 200 s; reversal 2000 s. Gradient cross-fade begins 5 s after restart.'}

def supplement1():
 d=pd.read_csv(OLD/'chiou_diagnostics.csv');d=d[(d.cells==201)&(d.time==200)&d.two_peaks&((d.latency+40>d.time)|d.censored)];pool=[r for r in poolrows() if r.get('n')==201 and r['seed']==1 and r.get('status')!='no_matching_root']
 fixed=[r for p in (HERE/'data').glob('fixed_*.json') for r in json.loads(p.read_text()) if r['n']==201 and r['seed']==1]
 fig,axs=plt.subplots(3,3,figsize=(14,12));fig.subplots_adjust(left=.075,right=.97,bottom=.13,top=.93,wspace=.35,hspace=.55)
 for name,label,c in zip(MODELS,NAMES,COLORS):
  a=d[d.model==name].copy();a['pool_fraction']=a.inactive_at_peak/a.initial_inactive_mean
  for ax,key in zip(axs[0],['smaller_peak_ratio','width_mean','pool_fraction']):
   q=a.groupby('L')[key].agg(['mean','min','max']);ax.plot(q.index,q['mean'],'o-',c=c,ms=4,label=label);ax.fill_between(q.index,q['min'],q['max'],color=c,alpha=.15)
 for i,(ax,name,y) in enumerate(zip(axs[0],['Peak approach to saturation','Peak broadening','Inactive pool under the peak'],['Smaller peak / saturation height','Mean half-height width (µm)','Inactive / initial concentration'])):
  title(ax,'ABC'[i],name);ax.set(xlabel='Domain length L (µm)',ylabel=y);ax.grid(alpha=.12)
 axs[0,0].axhline(1,c='.5',ls=':');axs[0,0].legend(frameon=False,fontsize=8)
 for i,(name,label,c) in enumerate(zip(MODELS,NAMES,COLORS)):
  ax=axs[1,i]
  for assay,marker in [('competing','o'),('gradient','s')]:
   for arm,ls in [('concentration','-'),('matched_pool',':'),('fixed_peak','--')]:
    rr=sorted([r for r in (fixed if arm=='fixed_peak' else pool) if r['model']==name and r['arm']==arm and r['assay'].startswith(assay)],key=lambda r:r['L']);key='competition_resolution_latency_s' if assay=='competing' else 'reversal_latency_s'
    ax.plot([r['L'] for r in rr],[r.get('competition_latency_from_original_withdrawal_s' if arm=='fixed_peak' and assay=='competing' else key) if r.get('competition_latency_from_original_withdrawal_s' if arm=='fixed_peak' and assay=='competing' else key) is not None else np.nan for r in rr],ls=ls,marker=marker,ms=4,c=c,alpha=1 if assay=='competing' else .55)
  ax.set(yscale='log',xlabel='Domain length L (µm)',ylabel='Latency (s)');title(ax,'DEF'[i],label+': pool controls');ax.grid(alpha=.12)
 fig.legend(handles=[Line2D([0],[0],ls=s,c='.3',label=l) for s,l in [('-','Fixed concentration'),(':','Matched inactive mean'),('--','Same physical peaks')]]+[Line2D([0],[0],marker=m,c='.3',ls='',label=l) for m,l in [('o','Competition'),('s','Gradient reversal')]],loc='lower center',ncol=5,frameon=False,fontsize=8)
 for i,(model,label) in enumerate(zip(MODELS,NAMES)):
  ax=axs[2,i]
  for arm,ls in [('concentration','-'),('matched_pool',':'),('fixed_peak','--')]:
   rr=sorted([r for r in (fixed if arm=='fixed_peak' else pool) if r['model']==model and r['arm']==arm and r['assay'].startswith('competing')],key=lambda r:r['L'])
   peak=[r['initial_active_peak_mass']/r['reference_active_peak_mass'] if arm=='fixed_peak' else r['anchor_peak_excess_mass']/r['reference_peak_excess_mass'] for r in rr]
   inactive=[r['initial_inactive_mean']/r['reference_inactive_mean'] if arm=='fixed_peak' else r['anchor_inactive_mean']/r['target_inactive_mean'] for r in rr]
   ax.plot([r['L'] for r in rr],peak,ls=ls,c=TEAL,marker='o',ms=3)
   ax.plot([r['L'] for r in rr],inactive,ls=ls,c=ORANGE,marker='s',ms=3)
  ax.axhline(1,c='.6',ls=':',lw=.7);ax.set(xlabel='Domain length L (µm)',ylabel='Assay-start quantity / L=15 value');title(ax,'GHI'[i],label+': what the control preserves');ax.grid(alpha=.12)
  if i==0:ax.legend(handles=[Line2D([0],[0],c=TEAL,marker='o',label='Active peak excess mass'),Line2D([0],[0],c=ORANGE,marker='s',label='Inactive mean')],frameon=False,fontsize=7)
 save(fig,'supplement1')
 with (OUT/'fixed_peak_control_diagnostics.csv').open('w') as f:
  fields=['model','L','n','actual_cells','seed','assay','arm','anchor_s','initial_active_peak_mass','reference_active_peak_mass','initial_peak_height','reference_peak_height','initial_inactive_mean','reference_inactive_mean','mass_error','competition_resolution_latency_s','competition_latency_from_original_withdrawal_s','competition_established','reversal_latency_s','reversal_conditioned'];w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows([r for p in (HERE/'data').glob('fixed_*.json') for r in json.loads(p.read_text())])
 # Explicit validation of the control, rather than claiming constant peaks from the mass law.
 with (OUT/'pool_control_diagnostics.csv').open('w') as f:
  fields=['model','L','n','seed','assay','arm','anchor_s','anchor_inactive_mean','target_inactive_mean','inactive_mean_relative_error','anchor_peak_excess_mass','reference_peak_excess_mass','anchor_peak_height','mass_error','competition_resolution_latency_s','reversal_latency_s'];w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(poolrows())


def supplement_tracking():
 m=module(ROOT/'polarity1d/scripts/build_wp1_manuscript_figures.py');d=rows(DATA/'manuscript/figure1_dynamic_tracking.csv');fig=plt.figure(figsize=(14,11));gs=fig.add_gridspec(3,2,left=.07,right=.96,bottom=.07,top=.95,hspace=.5,wspace=.32)
 wave=pd.read_csv(DATA/'manuscript/figureS2_tracking_waveforms.csv')
 for i,T in enumerate([100,400]):
  ax=fig.add_subplot(gs[0,i]);a=wave[wave.period_s==T];ax.plot(a.time_over_period,a.cue_orientation,c='.6',ls='--',label='Stimulus');ax.plot(a.time_over_period,a.response_orientation,c=TEAL,label='Polarity');ax.plot(a.time_over_period,a.first_harmonic,c=ORANGE,ls=':',label='First harmonic');title(ax,'AB'[i],f'Wave-pinning, L = 15 µm, T = {T} s');ax.set(xlabel='Time / period',ylabel='Signed orientation');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.12)
 sub=gs[1:,:].subgridspec(4,2,hspace=.10,wspace=.16);axes=[[fig.add_subplot(sub[i,j]) for j in range(2)] for i in range(4)];m._tracking_facets(axes,d);m._panel(axes[0][0],'C');save(fig,'supplement3_tracking')


def supplement2():
 spring=OLD/'zmurchok2020';threshold=pd.read_csv(spring/'results/thresholds.csv');allr=[json.loads(p.read_text()) for p in (HERE/'data').glob('spring_*.json')]
 fig=plt.figure(figsize=(15,11));gs=fig.add_gridspec(3,3,left=.07,right=.96,bottom=.08,top=.94,wspace=.4,hspace=.6)
 ax=fig.add_subplot(gs[0,0])
 for fb,c,label in [(0,ORANGE,'No tension feedback'),(80,TEAL,'Tension feedback')]:
  z=np.load(spring/f'results/validation_201_{fb}.npz');u=z['y'][-1,:-2:2];ax.plot(np.linspace(0,1,len(u)),u,c=c,lw=2,label=label)
 title(ax,'A','Author-code validation');ax.set(xlabel='Position / current length',ylabel='Active Rac at t = 500');ax.legend(frameon=False,fontsize=8)
 ax=fig.add_subplot(gs[0,1])
 for n,c,marker in [(101,'.6','o'),(201,NAVY,'s')]:
  rr=threshold[threshold.n==n];mid=(rr.lower+rr.upper)/2;ax.errorbar(rr.L0,mid,yerr=np.vstack([mid-rr.lower,rr.upper-mid]),fmt=marker+'-',c=c,ms=4,capsize=3,label=f'{n} points')
 title(ax,'B','Established-polarity reversal');ax.set(xlabel='Resting length / reference',ylabel='Opposing-cue threshold (model time)');ax.legend(frameon=False,fontsize=8)
 ax=fig.add_subplot(gs[0,2])
 for n,c,marker in [(101,'.6','o'),(201,NAVY,'s')]:
  rr=sorted([r for r in allr if r['n']==n],key=lambda r:r['L']);ax.plot([r['L'] for r in rr],[r['resolution_latency'] if r['resolution_latency'] is not None else np.nan for r in rr],marker+'-',c=c,ms=4,label=f'{n} points')
  for r in rr:
   if r['resolution_latency'] is None:ax.scatter(r['L'],3200 if r['established'] else 3,c=c,marker='^' if r['established'] else 'x')
 title(ax,'C','Peak competition after dual cues');ax.set(yscale='log',ylim=(1,5000),xlabel='Resting length / reference',ylabel='Resolution latency (model time)');ax.legend(frameon=False,fontsize=8);ax.text(.02,.04,'▲ unresolved; × no two peaks',transform=ax.transAxes,fontsize=7)
 z=np.load(HERE/'data/spring_L1_n201.npz');t=z['time'];y=z['y'];x=np.linspace(0,1,201)
 for j,(field,label,cmap) in enumerate([('stimulus','Applied dual-cue stimulus','YlOrBr'),('active','Active Rac','viridis'),('tension','Global spring/tension feedback','magma')]):
  ax=fig.add_subplot(gs[1,j]);values=y[:,:-2:2] if field=='active' else (np.tile(z[field][:,None],(1,201)) if field=='tension' else z[field]);im=field_kymo(ax,t,x,values,label,cmap);ax.set_xlabel('Model time from dual cue');ax.set_xscale('symlog',linthresh=40);ax.axvline(40,c='white',ls='--',lw=.8);fig.colorbar(im,ax=ax,pad=.02,fraction=.035);title(ax,'DEF'[j],label)
 ax=fig.add_subplot(gs[2,0]);ax.plot(t,z['protrusion_left'],c=TEAL,label='Left protrusion');ax.plot(t,z['protrusion_right'],c=ORANGE,label='Right protrusion');ax.plot(t,z['retraction'],c=NAVY,ls='--',label='Spring retraction');ax.set_xscale('symlog',linthresh=40);ax.set_xlim(0,3040);ax.set(xlabel='Model time from dual cue',ylabel='Boundary speed contributions');ax.legend(frameon=False,fontsize=8);title(ax,'G','Protrusion and retraction')
 ax=fig.add_subplot(gs[2,1]);ax.plot(t,z['polarity'],c=TEAL,label='Rac polarity');ax.axhline(0,c='.6',ls=':');ax.set_xscale('symlog',linthresh=40);ax.set_xlim(0,3040);ax.set(xlabel='Model time from dual cue',ylabel='Rac first moment / active mass');title(ax,'H','Polarity emerging from competition')
 ax=fig.add_subplot(gs[2,2]);ax.plot(t,z['loser_ratio'],c=NAVY);ax.axhline(.1,c=ORANGE,ls='--',label='Resolution cutoff');ax.axvline(40,c='.6',ls=':');ax.set_xscale('symlog',linthresh=40);ax.set_xlim(0,3040);ax.set(xlabel='Model time from dual cue',ylabel='Smaller / larger peak excess');ax.legend(frameon=False,fontsize=8);title(ax,'I','Two peaks to a sustained winner')
 for ax in fig.axes:
  if ax.get_ylabel() not in ['']:ax.grid(alpha=.1)
 fig.suptitle('Independent spring/tension model: reversal, competition and the mechanical response',fontsize=15,fontweight='bold',color=NAVY,y=.995);save(fig,'supplement2')
 EXPORT['supplement2']={'competition_conditions':allr,'example':'spring_L1_n201.npz','mechanics':'protrusion = author F(active Rac at boundary); retraction = 0.01*(length-L0); tension feedback = 80*(length-L0)'}

def main():
 selected=sys.argv[1:] or ['1','2','3','S1','S2','S3']
 existing=OUT/'provenance.json'
 if existing.exists():
  try: EXPORT.update(json.loads(existing.read_text()))
  except (json.JSONDecodeError,OSError): pass
 for name in selected:{'1':figure1,'2':figure2,'3':figure3,'S1':supplement1,'S2':supplement2,'S3':supplement_tracking}[name]();print('Built',name,flush=True)
 EXPORT['source_code_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.py')}
 EXPORT['generated_utc']=datetime.now(timezone.utc).isoformat();(OUT/'provenance.json').write_text(json.dumps(EXPORT,indent=2))
if __name__=='__main__':main()
