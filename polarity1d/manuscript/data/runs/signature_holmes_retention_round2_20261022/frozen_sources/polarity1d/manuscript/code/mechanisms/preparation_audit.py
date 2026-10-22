"""Compare the three saved basal preparations without resimulating them."""
from pathlib import Path
import sys
import json
import csv
import itertools
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE/'core'))
import native as n
import settling

OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/preparation_audit'
FIG=n.ROOT/'figures/fixed_slope_mechanisms'
MODELS=['wave_pinning','debelly','holmes_model3']
CONDITIONS=['fixed_10pct','fixed_20pct','relative_20pct']


def main():
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    campaign=n.DATA/'runs/fixed_slope_reversal_20261006'
    manifest=json.loads((campaign/'manifest.json').read_text())
    endpoint=list(csv.DictReader((campaign/'results/endpoints.csv').open()))
    groups={}
    for c in manifest['cases']:
        if c['model'] in MODELS:groups.setdefault((c['model'],c['length']),[]).append(c)
    pairs=[];aggregate=[]
    for (model,length),cases in groups.items():
        if not all(c['eligible'] for c in cases):continue
        fields={};states={}
        for c in cases:
            path=n.ROOT/c['checkpoint']
            if n.sha(path)!=c['checkpoint_sha256']:raise RuntimeError('Input checkpoint hash changed')
            state,tr=n.read_checkpoint(path,model)
            fields[c['condition']]={k:v[-1] for k,v in settling.monitored(tr,model,length).items()}
            states[c['condition']]=state
        maximum=0.;worst_field=None
        for a,b in itertools.combinations(cases,2):
            for field in fields[a['condition']]:
                absolute=float(np.max(np.abs(fields[a['condition']][field]-fields[b['condition']][field])))
                scaled=absolute/a['fixed_scales'][field]
                pairs.append(dict(model=model,length=length,condition_a=a['condition'],condition_b=b['condition'],field=field,maximum_absolute_difference=absolute,maximum_scaled_difference=scaled))
                if scaled>=maximum:maximum=scaled;worst_field=field
        duration_differences={}
        for gate in ['orientation','full_polarity']:
            values=[float(r['endpoint']) for r in endpoint if r['model']==model and float(r['length'])==length and r['definition']==gate]
            duration_differences[gate+'_duration_range']=max(values)-min(values)
        aggregate.append(dict(model=model,length=length,maximum_scaled_preparation_difference=maximum,worst_field=worst_field,
                              initial_orientation_range=float(np.ptp([c['initial_orientation'] for c in cases])),
                              **duration_differences))
    for name,rows in [('pairwise_fields.csv',pairs),('summary.csv',aggregate)]:n.csvwrite(OUT/name,rows)
    summary={}
    for model in MODELS:
        rows=[r for r in aggregate if r['model']==model]
        summary[model]=dict(eligible_lengths=len(rows),maximum_scaled_preparation_difference=max(r['maximum_scaled_preparation_difference'] for r in rows),
                            maximum_initial_orientation_range=max(r['initial_orientation_range'] for r in rows),
                            maximum_orientation_duration_range=max(r['orientation_duration_range'] for r in rows),
                            maximum_full_polarity_duration_range=max(r['full_polarity_duration_range'] for r in rows))
    n.write(OUT/'audit.json',dict(generated_at=n.now(),complete=True,producer='finished',models=summary,
              method='Pointwise pairwise comparison of every monitored native field on frozen scales; deBelly displacement mean-centered and physical length included using the existing field convention.',
              source_manifest_sha256=n.sha(campaign/'manifest.json'),source_endpoints_sha256=n.sha(campaign/'results/endpoints.csv'),script_sha256=n.sha(Path(__file__))))
    fig,axes=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
    labels=['Wave-pinning','deBelly','Holmes 3']
    for ax,model,label in zip(axes,MODELS,labels):
        rows=[r for r in aggregate if r['model']==model]
        ax.semilogy([r['length'] for r in rows],[max(r['maximum_scaled_preparation_difference'],1e-16) for r in rows],'o-',color='#0072B2',ms=4)
        ax.axhline(1e-4,color='#D55E00',ls='--',lw=1,label='Preparation window tolerance')
        ax.set_title(label);ax.set_xlabel('Cell length');ax.set_ylabel('Largest difference / frozen field scale');ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Three preparations converge to nearly the same native polarized state',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(FIG/f'preparation_equivalence.{ext}',dpi=220)
    plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
