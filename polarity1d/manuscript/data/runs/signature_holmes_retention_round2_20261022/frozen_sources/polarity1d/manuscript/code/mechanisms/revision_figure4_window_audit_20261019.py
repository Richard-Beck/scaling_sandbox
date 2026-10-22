"""Finite-window F4 sensitivity, without a stationarity certification test."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os

for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS=ROOT/'polarity1d/manuscript'
OUT=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/postburn_signature_candidate_20261019'
MODELS=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes 3',spring_ring='Spring ring')
METRICS=('qualified_multiregion_fraction','low_contrast_fraction','mean_R','mean_A','mean_Q','mean_contrast_ratio')
BAGS=2000


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def corr(r):
    eligible=r.pairs_eligible.sum();available=r.pairs_available.sum()
    return float(r.direction_cosine_sum.sum()/eligible) if eligible>=100 and eligible/available>=.01 else np.nan


def interval(a):
    take=np.isfinite(a);fraction=float(take.mean())
    if fraction<.95:return np.nan,np.nan,fraction
    lo,hi=np.quantile(a[take],[.025,.975]);return float(lo),float(hi),fraction


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=OUT);args=parser.parse_args();out=args.out
    data=out/'supporting/data';s=pd.read_csv(data/'postburn_scalar_case_windows.csv');m=pd.read_csv(data/'postburn_memory_case_windows.csv')
    records=[]
    for ordinal,(key,g) in enumerate(s.groupby(['model','size','fraction','kind'],sort=True)):
        ids=g.sort_values('seed').case_id.drop_duplicates().tolist();n=len(ids)
        rng=np.random.default_rng(202610195499+ordinal);bag=rng.integers(0,n,size=(BAGS,n)) if key[3]=='noise' and n>=3 else None
        a=g[g.window=='quarter_3'].set_index('case_id').loc[ids];b=g[g.window=='quarter_4'].set_index('case_id').loc[ids]
        base=dict(zip(['model','size','fraction','kind'],key),seeds=n,contrast='quarter4_minus_quarter3')
        for metric in METRICS:
            delta=b[metric].to_numpy()-a[metric].to_numpy()
            lo,hi,defined=interval(delta[bag].mean(axis=1)) if bag is not None else (np.nan,np.nan,0.)
            records.append(dict(**base,metric=metric,lag=np.nan,estimate=float(delta.mean()),
                ci_low=lo,ci_high=hi,bootstrap_defined_fraction=defined,paired_seed_min=float(delta.min()),paired_seed_max=float(delta.max())))
        lag=200. if key[0]=='spring_ring' else 2000.
        mm=m[(m.model==key[0])&(m['size']==key[1])&(m.fraction==key[2])&(m.kind==key[3])&(m.lag==lag)]
        a=mm[mm.window=='quarter_3'].set_index('case_id').loc[ids];b=mm[mm.window=='quarter_4'].set_index('case_id').loc[ids]
        point=corr(b)-corr(a)
        pairs=0.;sample=np.full(BAGS,np.nan)
        if bag is not None:
            values=[]
            for rows in (a,b):
                count=rows.pairs_eligible.to_numpy()[bag].sum(axis=1);total=rows.pairs_available.to_numpy()[bag].sum(axis=1)
                sums=rows.direction_cosine_sum.to_numpy()[bag].sum(axis=1)
                valid=(count>=100)&(count/total>=.01);x=np.full(BAGS,np.nan);x[valid]=sums[valid]/count[valid];values.append(x)
            sample=values[1]-values[0]
        lo,hi,defined=interval(sample) if bag is not None else (np.nan,np.nan,0.)
        seeded=b.direction_correlation.to_numpy()-a.direction_correlation.to_numpy();supported=a.support.to_numpy()&b.support.to_numpy();seeded=seeded[supported]
        records.append(dict(**base,metric='direction_cosine',lag=lag,estimate=point,ci_low=lo,ci_high=hi,
            bootstrap_defined_fraction=defined,paired_seed_min=float(np.nanmin(seeded)) if len(seeded) else np.nan,
            paired_seed_max=float(np.nanmax(seeded)) if len(seeded) else np.nan,
            quarter3_eligible_fraction=float(a.pairs_eligible.sum()/a.pairs_available.sum()),
            quarter4_eligible_fraction=float(b.pairs_eligible.sum()/b.pairs_available.sum()),
            quarter3_direction=corr(a),quarter4_direction=corr(b)))
    d=pd.DataFrame(records);d.to_csv(data/'late_quarter_sensitivity.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
    with PdfPages(out/'supporting/vectors/late_quarter_sensitivity.pdf') as pdf:
        for metrics,title,ylim in ((('qualified_multiregion_fraction','low_contrast_fraction'),'Late-quarter organization sensitivity',(-1.05,1.05)),
                                  (('mean_R','direction_cosine'),'Late-quarter magnitude and direction sensitivity',(-2.05,2.05))):
            fig,axes=plt.subplots(4,3,figsize=(12.8,12.7),sharey=True)
            fig.subplots_adjust(left=.075,right=.98,top=.84,bottom=.12,hspace=.62,wspace=.2)
            fig.suptitle(title,x=.075,y=.974,ha='left',fontsize=18)
            fig.text(.075,.939,'Quarter4 minus quarter3: chemicals24000–32000 minus16000–24000; Spring2400–3200 minus1600–2400.',fontsize=10.1)
            colors=('#8b4f9f','#d08029') if metrics[0]=='qualified_multiregion_fraction' else ('#3d6b9e','#16887d')
            names={'qualified_multiregion_fraction':'Qualified multiple-region occupancy','low_contrast_fraction':'Low-contrast occupancy','mean_R':'Mean R','direction_cosine':'Direction cosine at native lag2000/200'}
            fig.legend(handles=[Line2D([],[],color=c,lw=1.6,label=names[k]) for k,c in zip(metrics,colors)],loc='upper left',bbox_to_anchor=(.068,.922),ncol=2,frameon=False,fontsize=10)
            fig.text(.075,.862,'Filled squares:n=12; open circles:n=2. Bars:95% whole-seed intervals where n≥3 and direction support≥95%; unavailable directions stay absent.',fontsize=9.2)
            for row,(model,label) in enumerate(MODELS.items()):
                for column,size in enumerate((1.,2.,3.) if model=='spring_ring' else (20.,40.,60.)):
                    ax=axes[row,column]
                    for metric,color in zip(metrics,colors):
                        q=d[(d.model==model)&(d['size']==size)&(d.kind=='noise')&(d.metric==metric)].sort_values('fraction')
                        ax.plot(q.fraction,q.estimate,color=color,lw=1.35)
                        ax.vlines(q.fraction,q.ci_low,q.ci_high,color=color,lw=.8,alpha=.7)
                        for _,r in q.iterrows():
                            if np.isfinite(r.estimate):ax.plot(r.fraction,r.estimate,marker='s' if r.seeds==12 else 'o',mfc=color if r.seeds==12 else 'white',mec=color,ls='none',ms=3.8)
                    ax.axhline(0,color='#8b969b',lw=.8);ax.set_xscale('log');ax.set_ylim(*ylim)
                    ax.set_xticks([.0001,.001,.01,.1,1.],['10⁻⁴','10⁻³','10⁻²','10⁻¹','1'])
                    ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=8)
                    ax.set_title(f'{label} · {"C₀" if model=="spring_ring" else "C"}={size:g}',loc='left',fontsize=10)
                    ax.set_xlabel('Whole-signal strength f',fontsize=9)
                    if column==0:ax.set_ylabel('Quarter4 − quarter3',fontsize=9)
            fig.text(.075,.072,'These are finite-window differences, not a stationarity pass/fail test. Intervals resample independent whole seeds with quarter pairing preserved.',fontsize=9.3)
            fig.text(.075,.046,'Quarter direction uses shorter native lag2000/200; main lag8000/800 leaves only one pair per seed in a quarter and cannot meet the reporting support gate.',fontsize=9.3)
            fig.text(.075,.02,'Use the exported quarter-specific pair coverage and seed counts. No continuous region lifetime or migration-velocity inference is made.',fontsize=9.3)
            stem='late_quarter_'+metrics[0];fig.savefig(out/'supporting/viewers'/(stem+'.png'),dpi=140);pdf.savefig(fig);plt.close(fig)
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),source=str(Path(__file__).relative_to(ROOT)),source_sha256=sha(__file__),
        input_scalar_sha256=sha(data/'postburn_scalar_case_windows.csv'),input_memory_sha256=sha(data/'postburn_memory_case_windows.csv'),
        rows=len(d),resamples=BAGS,base_rng_seed=202610195499,paired_quarters=True,stationarity_certification=False,native_integrations=0)
    (out/'supporting/metadata/LATE_QUARTER_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
