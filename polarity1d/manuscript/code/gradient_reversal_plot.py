"""Standalone manuscript supplement for the finite fixed-slope countergradient."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import figure_builder as f
from matplotlib.lines import Line2D


def render(campaign):
    result = campaign/'results'
    status = json.loads((result/'status.json').read_text())
    validation = json.loads((result/'validation.json').read_text())
    assert status['complete'] and validation['passed'] and validation['independent_audit_passed']
    manifest = json.loads((campaign/'manifest.json').read_text())
    assert manifest['experiment']=='fixed_physical_slope_countergradient'
    rows = f.read_endpoint_rows(result/'endpoints.csv')
    fits = f.fit_summary(rows)
    for row in fits:row['gradient_slope']=manifest['slope_by_model'][row['model']]
    with (result/'length_trends.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(fits[0]));writer.writeheader();writer.writerows(fits)
    output=f.ROOT/'figures/supplement_gradient_response';output.mkdir(parents=True,exist_ok=True)
    fig,axes=f.plt.subplots(5,3,figsize=(14.5,15.5))
    for i,model in enumerate(f.MODELS):
        for j,condition in enumerate(f.CONDITIONS):
            ax=axes[i,j];chosen=sorted([r for r in rows if r['model']==model and r['condition']==condition],key=lambda r:r['length'])
            for gate,color,marker,label in f.GATES:
                series=[r for r in chosen if r['definition']==gate]
                ax.plot([r['length'] for r in series],[r['endpoint'] if r['status']=='success' else np.nan for r in series],color=color,ls='-' if gate=='orientation' else '--',lw=1.3)
                for r in series:
                    if r['status']=='success':ax.plot(r['length'],r['endpoint'],marker,ms=4,color=color,mfc=color if r['selected_field_tail_pass'] else 'white')
                    elif r['status']=='censored':ax.plot(r['length'],r['maximum_tested_duration'],'^',ms=4,mfc='white',mec=color)
                    elif gate=='orientation':ax.plot(r['length'],-.1,'x' if r['status']=='ineligible' else '+',color='#999999',transform=ax.get_xaxis_transform(),clip_on=False)
            ax.set_yscale('log');ax.grid(axis='y',color='#eeeeee');ax.tick_params(labelsize=7)
            if i==0:ax.set_title('Preparation: '+f.CONDITION_LABELS[j],fontsize=10)
            if j==0:ax.set_ylabel(f.LABELS[model]+f'\ng = {manifest["slope_by_model"][model]:.5g}\nShortest passing tested duration',fontsize=9)
            ax.set_xlabel('Native resting length L0' if model=='spring' else 'Nominal cell length',fontsize=8)
    fig.suptitle('Fixed physical slope · finite countergradient from settled basal preparations',fontsize=14,y=.98)
    fig.text(.5,.955,'Same corrected checkpoints, recovery, scalar gates, duration grids and native-field diagnostics as the patch comparison',ha='center',fontsize=9)
    handles=[Line2D([0],[0],color=color,marker=marker,ls='-' if gate=='orientation' else '--',label=label) for gate,color,marker,label in f.GATES]
    handles += [Line2D([0],[0],color=f.GRAY,marker='^',mfc='white',ls='none',label='Censored within tested durations'),Line2D([0],[0],color=f.GRAY,marker='o',mfc='white',ls='none',label='Field-tail flag fails'),Line2D([0],[0],color=f.GRAY,marker='x',ls='none',label='Ineligible preparation')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.935),ncol=3,fontsize=8,frameon=False)
    fig.text(.5,.035,'Added stimulus g*x during the pulse, then basal recovery for 10,000 native units. No gradient reconditioning or crossfade.\nOne slope per model across every length and preparation; peak and mean increase with length. Spring uses instantaneous physical position.\nPreparation columns identify the original left patch, not the gradient geometry. Full strength compares scalar polarity, not reflected native profiles.\nSource generated '+status['generated_at']+'; producer finished. Open markers are scalar successes with failing field diagnostics.',ha='center',fontsize=7.5,linespacing=1.5)
    fig.subplots_adjust(left=.09,right=.985,top=.86,bottom=.135,hspace=.5,wspace=.25)
    f.save(fig,output,'fixed_slope_reversal',dpi=180)
    stamp=status['generated_at']
    (output/'index.html').write_text(f'<h1>Fixed-slope countergradient supplement</h1><p>Data generated {stamp}; producer finished. Independent trial audit passed.</p><a href="fixed_slope_reversal.pdf">PDF</a><img width="100%" src="fixed_slope_reversal.png"><p><a href="../../data/runs/{campaign.name}/results/length_trends.csv">Fits and membership</a> · <a href="../../data/runs/{campaign.name}/results/endpoints.csv">Readouts</a> · <a href="../../GRADIENT_REVERSAL.md">Protocol</a></p>')
    (result/'index.html').write_text(f'<h1>Fixed-slope countergradient</h1><p>Data generated {stamp}; producer finished. All trials passed independent audit.</p><a href="../../../../figures/supplement_gradient_response/index.html">Final supplement</a> · <a href="endpoints.csv">Readouts</a> · <a href="length_trends.csv">Fits</a> · <a href="independent_audit.json">Audit</a>')
    import build
    build.index()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign',type=Path,required=True)
    render(parser.parse_args().campaign.resolve())
