"""Dense early-lag heading ACF from unchanged complete native compact histories."""
from pathlib import Path
import os,sys,hashlib,json
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np,pandas as pd
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[4];RUN=ROOT/'polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022';OUT=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026/figure4/supporting'
MODELS=('wave_pinning','holmes_model3','spring_ring','goryachev')
conditions=[(m,f,c) for m in MODELS for f in (.005,.3) for c in ((1.,2.,3.) if m=='spring_ring' else (20.,40.,60.))]
def main(slot):
 model,f,size=conditions[slot];tau=5 if model=='spring_ring' else 50
 path=RUN/'results/compact_index.csv';idx=pd.read_csv(path);q=idx[(idx.model==model)&np.isclose(idx['size'],size)&np.isclose(idx.fraction,f)].sort_values('replicate')
 steps=np.unique(np.r_[np.arange(0,101,2),[1,120,160,200,250,320,400,500,640,800,1000,1280,1600,2000,2500,3000]])
 sums=[];counts=[];possible=[];circ=[]
 for r in q.itertuples(index=False):
  with np.load(r.compact) as z:
   P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float);t=z['time'];cad=float(np.median(np.diff(t)));good=np.isfinite(P)&(np.abs(P)>=.1);unit=np.divide(P,np.abs(P),where=np.abs(P)>0,out=np.zeros_like(P));C=z['circumference'];s=[];n=[];pos=[]
   for step in steps:
    a=unit[step:];b=unit[:-step] if step else unit;ok=good[step:]&(good[:-step] if step else good);s.append(float(np.real(a[ok]*np.conjugate(b[ok])).sum()));n.append(int(ok.sum()));pos.append(len(a))
   sums.append(s);counts.append(n);possible.append(pos);circ.append(float(C.mean()))
 sums=np.array(sums);counts=np.array(counts);possible=np.array(possible);estimate=sums.sum(0)/counts.sum(0)
 rng=np.random.default_rng(20261026+slot);weights=rng.multinomial(len(q),np.ones(len(q))/len(q),size=2000);den=weights@counts;boot=(weights@sums)/den;lo,hi=np.quantile(boot,[.025,.975],axis=0)
 out=pd.DataFrame(dict(model=model,f=f,C_label=size,lag=steps*cad,lag_over_tau=steps*cad/tau,C=estimate,C_ci_lower=lo,C_ci_upper=hi,n_seeds=len(q),coverage=counts.sum(0)/possible.sum(0),eligible_pairs=counts.sum(0),possible_pairs=possible.sum(0),C_mean=np.mean(circ),tau=tau))
 folder=OUT/'data/dense_acf_workers';folder.mkdir(parents=True,exist_ok=True);out.to_csv(folder/f'{slot:02d}.csv',index=False)
 src=q[['model','size','fraction','replicate','compact','compact_sha256']].copy();src['compact']=src.compact.map(lambda p:str(Path(p).relative_to(ROOT)));src.to_csv(folder/f'{slot:02d}_sources.csv',index=False)
 if slot==0:
  meta=dict(generated_utc=datetime.now(timezone.utc).isoformat(),compact_index=str(path.relative_to(ROOT)),compact_index_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),endpoint_eligibility='finite P and abs(P)>=0.1 at both endpoints',estimand='pooled eligible-pair mean cosine of principal heading difference; no contrast mask',bootstrap='2000 multinomial resamples of independent whole seeds, deterministic seed 20261026+slot',lags='0,1,even2–100,120,160,200,250,320,400,500,640,800,1000,1280,1600,2000,2500,3000 frame intervals; every frame interval=tau/10',new_simulations=False)
  (OUT/'metadata').mkdir(parents=True,exist_ok=True);(OUT/'metadata/dense_acf_provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':main(int(os.environ.get('SLURM_ARRAY_TASK_ID',sys.argv[1] if len(sys.argv)>1 else '0')))
