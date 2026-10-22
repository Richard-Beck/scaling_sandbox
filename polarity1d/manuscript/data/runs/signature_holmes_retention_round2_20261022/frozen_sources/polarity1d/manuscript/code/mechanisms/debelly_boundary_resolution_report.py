"""Boundary audit comparison figure; all original figures remain intact."""
from pathlib import Path
import json,csv,textwrap,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from debelly_boundary_resolution import OUT,LENGTHS,n,r,run,base,polarity
FIG=n.ROOT/'figures/debelly_boundary_resolution'

def main():
 assert r.read(OUT/'independent_audit.json')['passed'] and r.read(OUT/'budget_audit.json')['passed'];FIG.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 colors={'legacy_centered':'#777777','independent_endpoints':'#0072B2'};labels={'legacy_centered':'Legacy centered','independent_endpoints':'Independent ends'}
 fig,axes=plt.subplots(3,2,figsize=(15,13.8));fig.subplots_adjust(left=.09,right=.97,top=.9,bottom=.16,hspace=.52,wspace=.3)
 state=base(39.);state['rac'][-1]*=2
 for v in colors:
  tr=run(v,39.,state,n.ConstantInput(0),.5,interval=.5);x=tr.x
  axes[0,0].plot(x,tr.fields['velocity'][-1],color=colors[v],lw=2,label=labels[v]);axes[0,1].plot(x,tr.fields['tension'][-1],color=colors[v],lw=2)
 axes[0,0].axhline(0,color='#AAA',lw=.8);axes[0,0].set_title('A   Right Rac active; left Rac inactive',loc='left',weight='bold');axes[0,0].set_ylabel('Cortical velocity (µm/s)');axes[0,0].set_xlabel('Nominal position (µm)');axes[0,0].legend(frameon=False,fontsize=9)
 axes[0,1].set_title('B   Centering changes local tension',loc='left',weight='bold');axes[0,1].set_ylabel('Tension (native units)');axes[0,1].set_xlabel('Nominal position (µm)');axes[0,1].text(.04,.9,'Same end separation; different internal mechanics',transform=axes[0,1].transAxes,fontsize=9)
 summaries={(v,L):r.read(OUT/f'{v}_L{L:g}_summary.json') for v in colors for L in LENGTHS}
 for v in colors:
  axes[1,0].plot(LENGTHS,[summaries[v,L]['crossing_time'] for L in LENGTHS],color=colors[v],marker='o',lw=2,label=labels[v])
 axes[1,0].set_title('C   MCA zero-crossing time changes',loc='left',weight='bold');axes[1,0].set_xlabel('Length (µm)');axes[1,0].set_ylabel('Continuous cue time to MCA zero (s)');axes[1,0].legend(frameon=False,fontsize=9)
 for v,ls in [('legacy_centered','--'),('independent_endpoints','-')]:
  for stage,color in [('early','#BC7922'),('zero','#8954A3')]:
   vals=[next(a for a in summaries[v,L]['rows'] if a['stage']==stage)['centroid_speed'] for L in LENGTHS];axes[1,1].plot(LENGTHS,vals,color=color,ls=ls,marker='o' if ls=='-' else None,lw=1.8,label=('50 s' if stage=='early' else 'At MCA zero')+' · '+labels[v])
 axes[1,1].set_title('D   Nominal centroid speeds depend on the boundary rule',loc='left',weight='bold',fontsize=11);axes[1,1].set_xlabel('Length (µm)');axes[1,1].set_ylabel('Nominal MCA-centroid speed (µm/s)');axes[1,1].legend(frameon=False,fontsize=8)
 for v,ls in [('legacy_centered','--'),('independent_endpoints','-')]:
  for key,color,label in [('direct_rho_rms','#0072B2','Direct Rho stress'),('endpoint_rms','#D55E00','Prescribed endpoints')]:
   axes[2,0].plot(LENGTHS,[next(a for a in summaries[v,L]['rows'] if a['stage']=='zero')[key] for L in LENGTHS],color=color,ls=ls,lw=2,label=label+' · '+labels[v])
 axes[2,0].set_title('E   Direct flow contributions at MCA zero',loc='left',weight='bold',fontsize=11);axes[2,0].set_xlabel('Length (µm)');axes[2,0].set_ylabel('Component cortical-velocity RMS (µm/s)');axes[2,0].legend(frameon=False,fontsize=8)
 for L,color in zip(LENGTHS,plt.cm.viridis(np.linspace(.12,.9,len(LENGTHS)))):
  with np.load(OUT/f'independent_endpoints_L{L:g}_recovery.npz') as z:
   P=polarity(z['field_rac'],z['x'],L);reference=abs(summaries['independent_endpoints',L]['preparation_left_polarity']);axes[2,1].plot(z['time'],P/reference,color=color,lw=1.8,label=f'{L:g} µm')
 axes[2,1].axhline(1,color='#777',ls=':',lw=1);axes[2,1].set_xscale('symlog',linthresh=100);axes[2,1].set_xticks([0,100,1000,10000]);axes[2,1].set_xticklabels(['0','100','1,000','10,000']);axes[2,1].set_title('F   Withdrawal at MCA zero retains reversal',loc='left',weight='bold');axes[2,1].set_xlabel('Time since withdrawal (s)');axes[2,1].set_ylabel('Right Rac polarity / initial left magnitude');axes[2,1].legend(frameon=False,fontsize=8,ncol=3)
 for ax in axes.flat:ax.grid(alpha=.15)
 fig.suptitle('deBelly boundary audit: remove the dynamical centering projection',fontsize=18,weight='bold',y=.97)
 fig.text(.09,.933,'Independent local end commands recover the published one-sided limit; inactive ends remain constrained.',fontsize=11,color='#43525B')
 caption=('A–B: one 0.5 s step from homogeneous basal rest with right-end Rac activated; no external cue. In-loop centering imposes left motion and redistributes strain, so it is not a display-frame change. C–F: each variant independently polarizes from identical small MCA seeds over 60000 basal units, including a stable final 10000-unit holdout at all five lengths. The same slope-0.24 continuous gradient starts each variant’s own state. C–D distinguish legacy centered and independent prescribed endpoints; speeds derived from MCA first moments use nominal coordinates. E is frozen-state source superposition, not causal removal of Rho/protrusion feedback. F withdraws each corrected trial at its own first recorded MCA zero and recovers 10000 units; all five retain right orientation. The corrected branch preserves the original prescribed-boundary model’s inactive-end constraint, boundary relaxation, quasi-static mechanics and uniform MCA mass restoration. It is not a stress-free migrating-cell model. Original scientific sources and retained comparisons are unchanged.')
 fig.text(.09,.10,'\n'.join(textwrap.wrap(caption,165)),fontsize=9,va='top',linespacing=1.4)
 for ext in ['png','pdf','svg']:fig.savefig(FIG/f'boundary_resolution.{ext}',dpi=200)
 plt.close(fig);r.write(OUT/'report_provenance.json',dict(generated_at=n.now(),producer='finished',renderer_sha256=n.sha(Path(__file__)),native_audit_sha256=n.sha(OUT/'independent_audit.json'),budget_audit_sha256=n.sha(OUT/'budget_audit.json'),figure_sha256={ext:n.sha(FIG/f'boundary_resolution.{ext}') for ext in ['png','pdf','svg']}));(FIG/'caption.txt').write_text(caption+'\n')
 print(FIG)
if __name__=='__main__':main()
