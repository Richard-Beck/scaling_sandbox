"""Plots for finite-countercue width and independent mean/slope pilot."""
from pathlib import Path
import gzip,json,html
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import countercue_pilot as p
import runner as r

OUT=p.RESULTS
MODELS=r.MODELS
LABELS=['Wave-pinning','Otsuji','Goryachev','LEGI','deBelly','Holmes 3','Spring']
plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def save(fig,name):
 fig.savefig(OUT/(name+'.png'),dpi=160)
 fig.savefig(OUT/(name+'.pdf'))
 fig.savefig(OUT/(name+'.svg'))
 plt.close(fig)

def mark_series(ax,rows,label,color,linestyle='-'):
 rows=sorted(rows,key=lambda c:c['length']);x=[c['length'] for c in rows];y=[c['result']['endpoint'] if c['result']['status']=='success' else np.nan for c in rows]
 ax.plot(x,y,linestyle,color=color,label=label,lw=1.1)
 for c in rows:
  q=c['result'];L=c['length']
  if q['status']=='censored':ax.plot(L,q['maximum_tested_duration'],'^',mfc='none',mec=color,ms=6)
  else:
   marker='v' if q['lower_grid_limited'] else ('s' if q['selected_morphology']=='nonmonotone' else 'o')
   ax.plot(L,q['endpoint'],marker,color=color,ms=4,mfc='none' if q['lower_grid_limited'] else color)
 ax.set_yscale('log');ax.set_xlabel('Resting L0 (Spring); nominal length otherwise');ax.set_ylabel('Passing pulse duration')


def width_figure(cases,stamp):
 fig,axes=plt.subplots(7,2,figsize=(13,18))
 for i,(model,label) in enumerate(zip(MODELS,LABELS)):
  for j,dose in enumerate(['peak_matched','integral_matched']):
   ax=axes[i,j]
   for fraction,color in [(.1,'#176d9c'),(.2,'#db7b24')]:
    for mode,ls in [('fractional','--'),('physical','-')]:
     chosen=[c for c in cases if c['model']==model and c['family']=='width' and c['fraction']==fraction and c['width_mode']==mode and c['dose_mode']==('peak_matched' if mode=='fractional' else dose)]
     mark_series(ax,chosen,f'{mode}, {100*fraction:g}%',color,ls)
   original=[]
   for L in p.LENGTHS[model]:
    c=next(c for c in cases if c['model']==model and c['length']==L)
    original.append(dict(length=L,result=dict(status=c['result']['original_status'],endpoint=c['result']['original_endpoint'],maximum_tested_duration=max(r.PROTOCOL['duration_grids'][model]['refined_durations']),lower_grid_limited=False,selected_morphology='original')))
   mark_series(ax,original,'Original A*x/L','#777777',':')
   ax.set_title(label+(' — fixed peak A' if j==0 else ' — equal spatial integral'))
   if i==0:ax.legend(fontsize=7,ncol=2)
 fig.suptitle('Finite countercue: physical versus fractional patch width\nHollow upward triangle: no passing pulse. Hollow downward triangle: lower-grid limit. Square: nonmonotone final field.\nGenerated '+stamp+'; exact original cue-free release states; no length fits.',fontsize=10)
 fig.tight_layout(rect=(0,0,1,.94));save(fig,'width_thresholds')


def mean_figure(cases,stamp):
 fig,axes=plt.subplots(7,3,figsize=(14,18))
 for i,(model,label) in enumerate(zip(MODELS,LABELS)):
  selected=[c for c in cases if c['model']==model and c['family']=='mean_slope']
  values=[c['result']['endpoint'] for c in selected if c['result']['status']=='success']
  low=min(values) if values else .1;high=max(values) if values else 1.
  if high<=low:high=low*1.01
  norm=LogNorm(low,high);cmap=plt.get_cmap('viridis')
  for j,L in enumerate(p.LENGTHS[model]):
   ax=axes[i,j]
   for mi,mean in enumerate([.5,1.,2.]):
    for si,slope in enumerate([.25,.5,1.]):
     c=next(c for c in selected if c['length']==L and c['mean_factor']==mean and c['slope_factor']==slope);q=c['result']
     color=cmap(norm(q['endpoint'])) if q['status']=='success' else '#dddddd'
     ax.add_patch(Rectangle((si,mi),1,1,facecolor=color,edgecolor='white',hatch='//' if q['status']!='success' else None))
     text=f"{q['endpoint']:.3g}" if q['status']=='success' else 'none'
     if q['lower_grid_limited']:text='≤'+text
     if q['status']=='success' and q['selected_morphology']=='nonmonotone':text+=' *'
     ax.text(si+.5,mi+.5,text,ha='center',va='center',fontsize=9,color='white' if q['status']=='success' and norm(q['endpoint'])<.65 else '#222222')
   ax.set(xlim=(0,3),ylim=(3,0),xticks=np.arange(3)+.5,yticks=np.arange(3)+.5)
   ax.set_xticklabels(['0.25','0.5','1']);ax.set_yticklabels(['0.5','1','2']);ax.set_xlabel('Slope / (A / length bound)');ax.set_ylabel('Mean / A');ax.set_title(f'{label}: L={L:g}')
 fig.suptitle('Finite countercue: independently varied mean and physical slope\nCells show shortest passing tested duration; none = no passing pulse; * = nonmonotone final field.\nColor scale shared across lengths within each model. Generated '+stamp,fontsize=10)
 fig.tight_layout(rect=(0,0,1,.94));save(fig,'mean_slope_thresholds')


def selected_arrays(c):
 q=c['result'];path=r.ROOT/q['selected_trace'];prefix=f"trial{q['selected_trace_row']}_"
 with np.load(path) as z:return {name[len(prefix):]:z[name].copy() for name in z.files if name.startswith(prefix)}


def profiles_figure(cases,stamp):
 fig,axes=plt.subplots(7,3,figsize=(15,18))
 for i,(model,label) in enumerate(zip(MODELS,LABELS)):
  L=p.LENGTHS[model][-1]
  chosen=[next(c for c in cases if c['model']==model and c['length']==L and c['family']=='width' and c['fraction']==.2 and c['width_mode']==mode and c['dose_mode']=='peak_matched') for mode in ['fractional','physical']]
  chosen.append(next(c for c in cases if c['model']==model and c['length']==L and c['family']=='mean_slope' and c['mean_factor']==1 and c['slope_factor']==1))
  for j,c in enumerate(chosen):
   ax=axes[i,j];q=c['result'];z=selected_arrays(c);x=z['x']/L
   ax.plot(x,z['pulse_active'][0],':',color='#777777',label='Initial release state')
   ax.plot(x,z['pulse_active'][-1],'--',color='#db7b24',label='At cue withdrawal')
   ax.plot(x,z['recovery_active'][-1],'-',color='#176d9c',label='After 10000 recovery')
   title=['Fractional width 20%, peak A','Physical width 20% of reference, peak A','Mean A, slope A / bound'][j]
   ax.set_title(f'{label}, L={L:g}: {title}\nD={q["selected_duration"]:.4g}; {q["status"]}');ax.set_xlabel('Relative native/material position');ax.set_ylabel('Native active/response field')
   if i==0 and j==0:ax.legend(fontsize=7)
 fig.suptitle('Spatial outcomes at the longest pilot length\nShortest passing tested pulse; maximum tested pulse where no success. Generated '+stamp,fontsize=10)
 fig.tight_layout(rect=(0,0,1,.95));save(fig,'spatial_outcomes')


def recovery_figure(cases,stamp):
 fig,axes=plt.subplots(7,2,figsize=(13,18))
 for i,(model,label) in enumerate(zip(MODELS,LABELS)):
  L=p.LENGTHS[model][-1]
  chosen=[next(c for c in cases if c['model']==model and c['length']==L and c['family']=='width' and c['fraction']==.2 and c['width_mode']==mode and c['dose_mode']=='peak_matched') for mode in ['fractional','physical']]
  chosen.append(next(c for c in cases if c['model']==model and c['length']==L and c['family']=='mean_slope' and c['mean_factor']==1 and c['slope_factor']==1))
  for c,color in zip(chosen,['#176d9c','#db7b24','#27804b']):
   z=selected_arrays(c);q=c['result'];title=c['condition']
   axes[i,0].plot(z['recovery_time'],z['recovery_orientation'],'.-',color=color,label=title,ms=2)
   y=z['recovery_active'];scale=c['fixed_scales'][{'legi':'response','debelly':'rac','holmes_model3':'rac_active'}.get(model,'active')]
   distance=np.max(abs(y-y[-1]),axis=1)/scale
   axes[i,1].plot(z['recovery_time'],np.maximum(distance,1e-14),'.-',color=color,ms=2)
  axes[i,0].axhline(.05,color='gray',ls='--',lw=.7);axes[i,0].set_xscale('symlog',linthresh=2.5);axes[i,0].set_title(f'{label}, L={L:g}: recovery polarity');axes[i,0].set_xlabel('Time after cue withdrawal');axes[i,0].set_ylabel('Orientation')
  axes[i,1].set_xscale('symlog',linthresh=2.5);axes[i,1].set_yscale('log');axes[i,1].set_title('Active field distance from final profile');axes[i,1].set_xlabel('Time after withdrawal');axes[i,1].set_ylabel('Max difference / fixed scale')
  if i==0:axes[i,0].legend(fontsize=6)
 fig.suptitle('Recovery trajectories: stable orientation versus spatial relaxation\nSampled histories; same three cues as spatial-outcome figure. Generated '+stamp,fontsize=10)
 fig.tight_layout(rect=(0,0,1,.95));save(fig,'recovery_behaviour')


def main():
 records=json.load(gzip.open(OUT/'records.json.gz','rt'));cases=records['cases'];stamp=r.now()
 width_figure(cases,stamp);mean_figure(cases,stamp);profiles_figure(cases,stamp);recovery_figure(cases,stamp)
 groups=[]
 for model in MODELS:
  for family in ['width','mean_slope']:
   cs=[c for c in cases if c['model']==model and c['family']==family]
   groups.append(dict(model=model,family=family,conditions=len(cs),success=sum(c['result']['status']=='success' for c in cs),censored=sum(c['result']['status']=='censored' for c in cs),lower_grid_limited=sum(c['result']['lower_grid_limited'] for c in cs),nonmonotonic=sum(c['result']['nonmonotonic_success_pattern'] for c in cs),successful_nonmonotone_profile=sum(c['result']['status']=='success' and c['result']['selected_morphology']=='nonmonotone' for c in cs)))
 r.csvwrite(OUT/'model_summary.csv',groups)
 write=p.write
 write(OUT/'plot_provenance.json',dict(generated_at=stamp,source_sha256=r.sha(Path(__file__).resolve()),records_sha256=r.sha(OUT/'records.json.gz'),no_length_regressions=True,figures=['width_thresholds','mean_slope_thresholds','spatial_outcomes','recovery_behaviour']))
 page=f'<html><head><meta charset="utf-8"><title>Finite-countercue pilot</title><style>body{{font-family:sans-serif;max-width:1500px;margin:30px auto;line-height:1.5}}img{{width:100%}}</style></head><body><h1>Finite-countercue width and mean/slope pilot</h1><p>Generated {stamp}; producer finished. Exact original cue-free release checkpoints; no gradient reversal.</p><p>Thresholds retain the original orientation/contrast gate. Full-field tail diagnostics and morphology are reported separately. No length regressions on these three-point pilot series.</p><p><a href="endpoints.csv">Endpoints and diagnostics</a> · <a href="validation.json">Validation</a> · <a href="../manifest.json">Design and sampled-history definitions</a></p>'
 for name,title in [('width_thresholds','Physical versus fractional width'),('mean_slope_thresholds','Independent mean and slope'),('spatial_outcomes','Spatial outcomes'),('recovery_behaviour','Recovery behaviour')]:page+=f'<h2>{title}</h2><a href="{name}.pdf">PDF</a><img src="{name}.png">'
 page+='</body></html>'
 (OUT/'index.html').write_text(page)
 print(json.dumps(dict(generated_at=stamp,figures=4,conditions=len(cases))),flush=True)

if __name__=='__main__':main()
