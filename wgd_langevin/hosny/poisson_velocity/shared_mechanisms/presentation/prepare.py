from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent;SHARED=ROOT.parent;POISSON=SHARED.parent
sys.path.insert(0,str(POISSON/'exhaustive'))
import fit_frontier as ff
from run import write_csv
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
CHOICES=[('null',0,0,'K=0 · no KD changes'),('direction',8,8,'K=2 · HIF κ / RhoA κ'),('memory',8,1,'K=2 · HIF κ / RhoA tau'),('full',15,15,'K=8 · all parameters')]
def main():
    fits=json.load((POISSON/'exhaustive/results/fits.json').open());selected=[];pred=[];tracks={};pars=[];time=[]
    for model,h,r,label in CHOICES:
        score=0.
        for base in ('N','P'):
            f=next(f for f in fits if f['base']==base and f['hmask']==h and f['rmask']==r)
            score+=f['independent_score']/2
            _,qq,_,tt=ff.evaluate_parameters(base,f['parameters'],65536,77124,True)
            pred += [dict(model=model,**q) for q in qq];pars += [dict(model=model,**p) for p in f['parameters']]
            tracks.update({model+'_'+key:v for key,v in tt.items() if not key.endswith('_velocity')})
            for p in f['parameters']:
                for t in np.geomspace(.1,32,90):time.append(dict(model=model,group=p['group'],time=t,msd=float(ff.msd(t,p['tau'],p['rate'],p['amplitude'],p['kappa']))))
        selected.append(dict(model=model,label=label,hmask=h,rmask=r,K=h.bit_count()+r.bit_count(),score=score))
    # A common normoxic reference is independent of mechanism selection.
    old=json.load((POISSON/'results/optimization.json').open());p=next(r for r in old if r['base']=='O')['parameters'][0]
    m=ff.metrics(ff.simulate(p['tau'],p['rate'],p['amplitude'],p['kappa'],n=65536,seed=77124))
    pars.append(dict(model='reference',**p))
    for key in ff.KM:
        for q,v in zip(ff.Q,np.quantile(m[key],ff.Q)):pred.append(dict(model='reference',group='O',metric=key,quantile=q,predicted=float(v)))
    write_csv(OUT/'selected_models.csv',selected);write_csv(OUT/'parameters.csv',pars);write_csv(OUT/'predictions.csv',pred);write_csv(OUT/'time_predictions.csv',time);np.savez_compressed(OUT/'tracks.npz',**tracks)
    meta=dict(selection=CHOICES,fit_changes='None. All four choices reuse the current exhaustive fits with identical subsets across populations; selected KD values may increase or decrease.',
              sampling='65536 homogeneous cells/state; 16h sampled at 0.5h; seed 77124; stationary Poisson-kick model. Normoxia is one common reference.',
              trajectories='12 baseline tracks per population/model; short curves are exact first-hour prefixes; common axis limits across models at each duration.',
              MSD='Analytic predictions for baseline and both KDs in both populations. No supplied MSD data.',
              sensitivity='Reused shared-frontier ratios and within-model ratio profiles. Profiles belong only to HIF kappa/RhoA kappa representative.',
              input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (POISSON/'exhaustive/results/fits.json',SHARED/'results/ratio_profiles.csv',SHARED/'results/frontier_ratios.csv')})
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2));print('Prepared four representatives; no refitting.')
if __name__=='__main__':main()
