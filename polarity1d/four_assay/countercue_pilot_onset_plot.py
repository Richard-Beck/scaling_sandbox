"""Cue-onset polarity traces from completed, immutable finite-countercue data."""
from pathlib import Path
import gzip,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import countercue_pilot as p
import runner as r
from countercue_pilot_recovery_plot import selected_arrays

OUT=p.RESULTS
LABELS=['Wave-pinning','Otsuji','Goryachev','LEGI','deBelly','Holmes 3','Spring']
COLORS=['#0072B2','#D55E00','#009E73']
EXPERIMENTS=[('fractional','peak_matched','Fractional width 20%, peak A'),('physical','peak_matched','Fixed physical width, peak A'),('physical','integral_matched','Fixed physical width, integral matched')]
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})

def main():
 p.verify()
 before=r.sha(OUT/'records.json.gz');cases=json.load(gzip.open(OUT/'records.json.gz','rt'))['cases'];stamp=r.now()
 fig,axes=plt.subplots(7,3,figsize=(16,22),sharey='row')
 checks=[];table=[]
 for i,(model,label) in enumerate(zip(r.MODELS,LABELS)):
  lengths=p.LENGTHS[model];chosen=[c for c in cases if c['model']==model and c['family']=='width' and c['fraction']==.2]
  minimum_duration=min(c['result']['selected_duration'] for c in chosen)
  maximum_duration=max(c['result']['selected_duration'] for c in chosen)
  row_low,row_high=-1.,1.
  for j,(mode,dose,title) in enumerate(EXPERIMENTS):
   ax=axes[i,j]
   ax.axhline(-1,color='#bbbbbb',ls=':',lw=.7);ax.axhline(1,color='#bbbbbb',ls=':',lw=.7)
   ax.axhline(0,color='#dddddd',lw=.7)
   for L,color in zip(lengths,COLORS):
    c=next(c for c in chosen if c['length']==L and c['width_mode']==mode and c['dose_mode']==dose)
    q=c['result'];z=selected_arrays(c);D=float(q['selected_duration'])
    initial=float(z['pulse_orientation'][0]);scale=abs(initial)
    assert initial<-.05
    np.testing.assert_allclose(z['pulse_time'][[0,-1]],[0,D],rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(z['recovery_time'][[0,-1]],[0,10000],rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(z['pulse_orientation'][-1],z['recovery_orientation'][0],rtol=1e-8,atol=1e-10)
    time=np.r_[z['pulse_time'],D+z['recovery_time'][1:]]
    polarity=np.r_[z['pulse_orientation'],z['recovery_orientation'][1:]]/scale
    assert np.all(np.diff(time)>0) and np.all(np.isfinite(polarity))
    row_low=min(row_low,float(np.min(polarity)));row_high=max(row_high,float(np.max(polarity)))
    success=q['status']=='success';ls='-' if success else '--'
    ax.plot(time,polarity,color=color,ls=ls,lw=1.6,label=f"{'L0' if model=='spring' else 'L'}={L:g}; D={D:.4g}"+(' (no pass)' if not success else ''))
    ax.plot(D,z['pulse_orientation'][-1]/scale,marker='|',color=color,ms=13,mew=2,ls='none',zorder=5)
    checks.append(dict(case_index=c['case_index'],model=model,length=L,experiment=dose+'_'+mode,selected_duration=D,status=q['status'],initial_orientation=initial,normalized_initial_polarity=float(polarity[0]),withdrawal_time=D,withdrawal_normalized_polarity=float(z['pulse_orientation'][-1]/scale),last_time=float(time[-1]),samples=len(time)))
    for t,y in zip(time,polarity):table.append(dict(case_index=c['case_index'],model=model,length=L,width_mode=mode,dose_mode=dose,status=q['status'],duration=D,time_from_cue_onset=float(t),normalized_polarity=float(y),initial_orientation=initial))
   ax.set_xscale('symlog',linthresh=minimum_duration/5)
   ax.set_xlim(0,maximum_duration+10000)
   lo,hi=ax.get_ylim();ax.set_ylim(min(lo,-1.12),max(hi,1.12))
   ax.legend(loc='upper left',fontsize=8,framealpha=.88,handlelength=1.8)
   ax.set_xlabel('Time from countercue onset (native units)',fontsize=8)
   ax.tick_params(labelsize=8)
   if j==0:ax.set_ylabel(label+'\nPolarity / |initial polarity|',fontsize=10)
   if i==0:ax.set_title(title,fontsize=11,pad=12)
  row_span=row_high-row_low
  axes[i,0].set_ylim(row_low-.06*row_span,row_high+.08*row_span)
 fig.suptitle('Countercue and recovery: polarity relative to the initial settled state',fontsize=16,y=.987)
 fig.text(.5,.966,'Three lengths per facet. Each trace starts at −1; +1 is equally strong reversed polarity. Coloured ticks mark cue withdrawal.',ha='center',fontsize=10)
 handles=[Line2D([0],[0],color=color,lw=2,label=label) for color,label in zip(COLORS,['Shortest length','Middle length','Longest length'])]
 handles += [Line2D([0],[0],color='#444444',lw=1.6,label='Shortest passing tested pulse'),Line2D([0],[0],color='#444444',lw=1.6,ls='--',label='Maximum tested pulse; none passed'),Line2D([0],[0],color='#444444',marker='|',ms=12,mew=2,ls='none',label='Cue withdrawal')]
 fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.953),ncol=3,frameon=False,fontsize=9)
 fig.text(.5,.025,'Fixed physical width = 20% of the reference-length domain; integral matching is to the fractional cue at the same length.\nRecovery continues for 10000 units after each pulse. Time axes compress late times; ranges and units differ between models.\nPolarity is an orientation measure, not a measure of complete spatial pattern replacement. Generated '+stamp,ha='center',fontsize=9)
 fig.subplots_adjust(left=.09,right=.985,top=.905,bottom=.075,wspace=.15,hspace=.40)
 for ext in ['png','pdf','svg']:fig.savefig(OUT/f'recovery_behaviour.{ext}',dpi=180)
 plt.close(fig)
 assert len(checks)==63
 assert before==r.sha(OUT/'records.json.gz')
 r.csvwrite(OUT/'recovery_behaviour_data.csv',table)
 p.write(OUT/'recovery_plot_provenance.json',dict(generated_at=stamp,source_sha256=r.sha(Path(__file__).resolve()),records_sha256=before,simulation_data_unchanged=True,normalization='orientation / abs(orientation at cue onset)',time_origin='countercue onset',withdrawal_mark='vertical tick at selected pulse duration',facets=21,traces=63,checks_passed=True,traces_checked=checks))
 from countercue_pilot_plot_revision import dashboard
 dashboard(stamp)
 p.write(OUT/'artifact_hashes.json',{f.name:r.sha(f) for f in sorted(OUT.iterdir()) if f.suffix in ['.png','.pdf','.svg','.csv']})
 print(json.dumps(dict(generated_at=stamp,facets=21,traces=63,checks_passed=True,simulation_data_unchanged=True)))

if __name__=='__main__':main()
