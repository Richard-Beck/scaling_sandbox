from pathlib import Path
import json,sys,numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT.parent/'exhaustive'))
import fit_frontier as ff
from run import write_csv
OUT=ROOT/'results';profiles=json.load((OUT/'profiles.json').open());models=json.load((ROOT.parent/'exhaustive/results/fits.json').open());rows=[]
for seed in (44112,44113):
    ref=0.
    for b in ('N','P'):
        f=next(f for f in models if f['base']==b and f['hmask']==8 and f['rmask']==8)
        ref+=ff.evaluate_parameters(b,f['parameters'],32768,seed)[0]/2
    for parameter in ('tau','rate','amplitude','kappa'):
        f=next(f for f in profiles if f['parameter']==parameter and f['ratio']==1);score=0.
        for b in ('N','P'):
            pars=[p for p in f['parameters'] if p['group'] in (b,'H'+b,'R'+b)]
            score+=ff.evaluate_parameters(b,pars,32768,seed)[0]/2
        rows.append(dict(parameter=parameter,seed=seed,reference=ref,equality_score=score,penalty=score-ref))
write_csv(OUT/'equal_ratio_validation.csv',rows)
print(rows)
