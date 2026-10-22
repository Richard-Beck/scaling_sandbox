"""Descriptive checks of existing reconstructed patch trajectories; no new simulations."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'data/reversal_fields'
OUTPUT=ROOT/'figures/diagnostic_reversal_fields'
PAIRS={'wave_pinning':[9,25],'holmes_model3':[14,24],'spring':[.8,2.4],'debelly':[19,35]}
LABELS={'wave_pinning':'Wave-pinning','holmes_model3':'Holmes 3','spring':'Spring','debelly':'deBelly'}
CONDS=['physical_peak_matched','fractional_peak_matched','physical_integral_matched']
rows=[]
fig,axes=plt.subplots(4,3,figsize=(15,13),layout='constrained')
for i,(model,lengths) in enumerate(PAIRS.items()):
 for condition in CONDS:
  for j,L in enumerate(lengths):
   stem=f'{model}_{condition}_L{L:g}'
   meta=json.loads((HERE/(stem+'.json')).read_text())
   with np.load(HERE/(stem+'_full.npz')) as z:
    t=z['time'];q=z['x']/L
    u=z['field_'+{'holmes_model3':'rac_active','debelly':'rac'}.get(model,'active')]
    initial=u[0];height=initial.max();base=initial.min();contrast=height-base
    left=u[:,q<.05].mean(1);right=u[:,q>.95].mean(1)
    valley=u[:,(q>=.05)&(q<=.95)].min(1)
    row=dict(model=model,condition=condition,length=L,pulse_duration=meta['duration'],initial_peak=height,initial_baseline=base,initial_half_height_width=np.trapezoid((initial>=base+.5*contrast).astype(float),q)*L,field_tail_pass=meta['selected_tail_field_pass'])
    for threshold in [.1,.25,.5]:
     both=(left-valley>threshold*contrast)&(right-valley>threshold*contrast)
     # Conservative sampled duration: count interval only if both endpoints pass.
     row[f'two_pole_duration_threshold_{threshold}']=np.sum(np.diff(t)*(both[1:]&both[:-1]))
    rows.append(row)
    if condition!='physical_peak_matched':continue
    color=['#2474b5','#dd7c25'][j];label=f'L={L:g}'
    ax=axes[i,0];ax.plot(t,left,color=color,ls='--',label=label+' old end');ax.plot(t,right,color=color,label=label+' new end');ax.axvline(meta['duration'],color=color,alpha=.4)
    ax.set_xscale('symlog',linthresh=max(.01,min(meta['duration'],1)))
    ax.set_ylabel(LABELS[model]+'\nActive/Rac, native units');ax.set_xlabel('Time; vertical lines = cue withdrawal');ax.legend(fontsize=7)
    ax=axes[i,1];ax.plot(q,initial,color=color,label=label);ax.set_xlabel('x / nominal length (Spring: material fraction)');ax.set_ylabel('Initial active/Rac, native units');ax.legend(fontsize=8)
    ax=axes[i,2]
    if model=='spring':
     feedback=z['field_tension_feedback'];ax.plot(t,feedback,color=color,label=label);ax.set_ylabel('Global tension-coupled deactivation increment')
    elif model=='holmes_model3':
     feedback=z['field_rho_active'];ax.plot(t,feedback[:,q<.05].mean(1),color=color,ls='--',label=label+' old end');ax.plot(t,feedback[:,q>.95].mean(1),color=color,label=label+' new end');ax.set_ylabel('Active Rho, native units')
    elif model=='debelly':
     feedback=z['field_delta_alpha_mca'];ax.plot(t,feedback[:,q<.05].mean(1),color=color,ls='--',label=label+' old end');ax.plot(t,feedback[:,q>.95].mean(1),color=color,label=label+' new end');ax.set_ylabel('MCA increment to Rac activation')
    else:
     feedback=z['field_inactive'];ax.plot(t,feedback.mean(1),color=color,label=label);ax.set_ylabel('Mean inactive substrate, native units')
    ax.axvline(meta['duration'],color=color,alpha=.4);ax.set_xscale('symlog',linthresh=max(.01,min(meta['duration'],1)));ax.set_xlabel('Native time');ax.legend(fontsize=7)
for ax in axes.flat:ax.grid(alpha=.15)
for ax,title in zip(axes[0],['Old/new pole amplitudes','Starting profiles','Substrate / inhibitory feedback']):ax.set_title(title)
fig.suptitle('Quick mechanism check: fixed physical width and peak; existing shortest full-polarity passing trials\nDifferent selected pulse durations; dashed = old end, solid = challenged end. Correlations do not establish causation.',fontsize=12)
fig.savefig(OUTPUT/'quick_mechanism_check.png',dpi=150)
plt.close(fig)
with (OUTPUT/'quick_mechanism_check.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
(OUTPUT/'quick_mechanism_check_definitions.json').write_text(json.dumps(dict(source='Existing 42-trial reconstructed report; analysis uses 24 successful examples from four models and three conditions',two_poles='Mean active concentration in first and last 5% of domain each exceeds the minimum interior concentration by threshold times INITIAL full-profile contrast; thresholds 0.1, 0.25, 0.5; count intervals only when both sampled endpoints pass',half_height_width='Trapezoidal integral of initial concentration >= minimum + half initial contrast; nominal physical length, Spring native nominal units',limits=['Chosen shortest passing pulses differ between lengths; trajectories are not matched-duration controls','Two-pole indicator is descriptive, not an eigenvalue or attractor analysis','Spring recovery samples every 2.5 native units; fast early changes have limited resolution','Scalar full-polarity success does not guarantee final complete-field stability','No model-specific saturation point is estimated','MCA and inhibitor traces establish association only'],generated_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()),indent=2))
print('Wrote quick_mechanism_check.png, .csv, and _definitions.json')
