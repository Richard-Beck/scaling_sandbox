"""All-length threshold figures, refreshed explicitly while the producer runs."""
from pathlib import Path
import html,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import runner as r
OUT=Path(__file__).resolve().parent/'pilots/width_full_lengths_20261002/results'
LABELS=['Wave-pinning','Otsuji','Goryachev','LEGI','deBelly','Holmes 3','Spring']
EXPERIMENTS=[('fractional','peak_matched','Fractional width 20%, peak A'),('physical','peak_matched','Fixed physical width, peak A'),('physical','integral_matched','Fixed physical width, integral matched')]
DEFINITIONS=[('orientation','#0072B2','o','Original orientation gate'),('full_polarity','#D55E00','s','Full polarity within 1%')]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def plot(cases,status):
 fig,axes=plt.subplots(7,3,figsize=(17,23));complete=status['complete']
 for i,(model,label) in enumerate(zip(r.MODELS,LABELS)):
  for j,(mode,dose,title) in enumerate(EXPERIMENTS):
   ax=axes[i,j];chosen=sorted([c for c in cases if c['model']==model and c['width_mode']==mode and c['dose_mode']==dose],key=lambda c:c['length'])
   maximum=max(r.PROTOCOL['duration_grids'][model]['refined_durations']);values=[]
   for name,color,marker,text in DEFINITIONS:
    lengths=[c['length'] for c in chosen];ys=[c['result'][name]['endpoint'] if c['result'][name]['status']=='success' else np.nan for c in chosen]
    ax.plot(lengths,ys,color=color,lw=1.4,ls='-' if name=='orientation' else '--',zorder=2)
    for c in chosen:
     q=c['result'][name]
     if q['status']=='success':
      values.append(q['endpoint']);ax.plot(c['length'],q['endpoint'],marker,ms=4,color=color,mfc=color if q['selected_tail_field_pass_0_1_percent'] else 'white',mew=1,zorder=3)
     elif q['status']=='censored':ax.plot(c['length'],maximum,'^',ms=5,mfc='white',mec=color,mew=1.0,zorder=3)
   ax.set_yscale('log');ax.set_xlim(min(lengths)-.03*(max(lengths)-min(lengths)),max(lengths)+.03*(max(lengths)-min(lengths)))
   if not values:ax.set_ylim(maximum/100,maximum*1.6)
   else:ax.set_ylim(min(values)/1.5,max(max(values),maximum if any(c['result'][name]['status']=='censored' for c in chosen for name,_,_,_ in DEFINITIONS) else max(values))*1.6)
   # A status strip in axes-relative height avoids implying a duration.
   for c in chosen:
    q=c['result']['orientation'];kind=q['status']
    if kind in ['ineligible','restricted']:
     ax.plot(c['length'],-.07,marker='x' if kind=='ineligible' else '|',color='#999999' if kind=='ineligible' else '#555555',ms=5 if kind=='ineligible' else 9,mew=1.2,transform=ax.get_xaxis_transform(),clip_on=False)
    elif kind=='pending':ax.plot(c['length'],-.07,'.',color='#bbbbbb',ms=3,transform=ax.get_xaxis_transform(),clip_on=False)
   ax.set_xlabel('Resting L0' if model=='spring' else 'Nominal cell length',labelpad=14,fontsize=8)
   if j==0:ax.set_ylabel(label+'\nShortest passing tested pulse',fontsize=10)
   ax.grid(axis='y',color='#eeeeee',lw=.6)
   ax.tick_params(labelsize=8)
   if i==0:ax.set_title(title,fontsize=11,pad=12)
   ineligible=sum(not c['eligible'] for c in chosen);restricted=sum(c['eligible'] and c['restricted'] for c in chosen)
   if ineligible or restricted:ax.text(.98,.04,f'{ineligible} ineligible; {restricted} restricted',transform=ax.transAxes,ha='right',fontsize=7,color='#666666')
 fig.suptitle('Finite countercue across all lengths: two reorientation definitions',fontsize=16,y=.985)
 fig.text(.5,.962,'Exact settled release states; three cues; recovery 10000 units. Full polarity must stay within 0.99–1.01 of reversed initial strength in the final 500 units.',ha='center',fontsize=10)
 handles=[Line2D([0],[0],color=color,marker=marker,ls='-' if name=='orientation' else '--',label=text) for name,color,marker,text in DEFINITIONS]
 handles += [Line2D([0],[0],color='#555555',marker='^',mfc='white',ls='none',label='No qualifying tested pulse'),Line2D([0],[0],color='#999999',marker='x',ls='none',label='Ineligible initial state (below axis)'),Line2D([0],[0],color='#555555',marker='|',ms=9,ls='none',label='Fixed patch exceeds half-cell (below axis)')]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.951),ncol=3,frameon=False,fontsize=9)
 stamp=status['generated_at'];state='finished' if complete else 'running — provisional endpoints'
 fig.text(.5,.025,'Physical width = 20% of the same pilot reference domain. Both gates share exactly the same tested durations.\nHollow circles/squares flag native-field tail variation >0.1% of fixed scales; polarity strength does not measure complete pattern replacement.\nMissing endpoints break lines. Native time units and y ranges differ between models. Generated '+stamp+'; producer '+state+'.',ha='center',fontsize=9)
 fig.subplots_adjust(left=.09,right=.985,top=.91,bottom=.08,hspace=.5,wspace=.20)
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'length_response_two_definitions.{ext}',dpi=170)
 plt.close(fig)
 # Paired figures retain the two-gate comparison while making each definition easy to inspect.
 counts=status.get('status_counts_by_definition',{})
 text=f'''<!doctype html><html><head><meta charset="utf-8"><title>All-length countercue comparison</title><style>body{{font-family:sans-serif;max-width:1700px;margin:30px auto;line-height:1.5}}img{{width:100%}}</style></head><body><h1>Three patch cues across all seven models and all lengths</h1><p>Generated {stamp}; producer {state}. Stage: {status['stage']}.</p><p>All 154 canonical model–length pairs: 462 cue conditions and 924 readouts. Ineligible settled states are retained. Fixed-width tests are restricted to patches fitting within one half of the cell; for Spring the nominal resting length supplies a safe minimum domain.</p><p>The original orientation/contrast/stability gate and the 1% polarity gate are scored from the same trials. The latter requires every original sample in the final 500 recovery units to have polarity between 0.99 and 1.01 times the magnitude of initial settled polarity, plus the original gate. Recovery lasts 10000 units. These are shortest passing tested durations; no global monotonicity or universal scaling law is assumed.</p><p>Reference lengths stay fixed at the pilot values: wave-pinning/Otsuji/Goryachev 21, LEGI 17, deBelly 27, Holmes 30, Spring L0=1.6 (physical release length defines its patch width). Fractional width is 20% of the instantaneous domain. Equal-integral physical cue matches that fractional cue at each length.</p><p>{html.escape(json.dumps(counts))}</p><p><a href="endpoints.csv">Endpoints</a> · <a href="status.json">Current status</a> · <a href="../manifest.json">Design and provenance</a> · <a href="validation.json">Final validation</a></p><a href="length_response_two_definitions.pdf">PDF</a> · <a href="length_response_two_definitions.svg">SVG</a><img src="length_response_two_definitions.png"></body></html>'''
 (OUT/'index.html').write_text(text)
