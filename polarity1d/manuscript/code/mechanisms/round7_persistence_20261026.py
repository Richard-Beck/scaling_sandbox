"""Round7 Figures4–5 persistence analysis and display, retaining native evidence."""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys
from datetime import datetime, timezone
from pathlib import Path
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[_k]='1'
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[4]
DEV=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev'
OUT=DEV/'revision_round7_20261026'
BASE=DEV/'revision_round6_20261023'
RUN=ROOT/'polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022'
MODELS=('wave_pinning','holmes_model3','spring_ring','goryachev')
NAMES=dict(wave_pinning='Wave-pinning',holmes_model3='Holmes',spring_ring='Spring',goryachev='Goryachev')
WEAK=dict(wave_pinning=.03,holmes_model3=.02,spring_ring=.005,goryachev=.005)
COLORS=dict(wave_pinning='#17619c',holmes_model3='#28783b',spring_ring='#aa7a00',goryachev='#b93434')
SOURCES={}
def now(): return datetime.now(timezone.utc).isoformat()
def rel(p): return str(Path(p).resolve().relative_to(ROOT))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    SOURCES[rel(p)]=sha(p)
    return pd.read_csv(p)
def write(p,v):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def setup(n):
    folder=OUT/f'figure{n}'
    for name in ('data','metadata','vectors','viewers'): (folder/'supporting'/name).mkdir(parents=True,exist_ok=True)
    return folder

def extract():
    folder=setup(5); data=folder/'supporting/data';meta=folder/'supporting/metadata'
    idx=read(RUN/'results/compact_index.csv')
    bins_rows=[];coh_rows=[];sources=[];waiting=[]
    for model in MODELS:
        tau=5 if model=='spring_ring' else 50
        for f in (WEAK[model],.3):
            for size in ((1.,3.) if model=='spring_ring' else (20.,60.)):
                q=idx[(idx.model==model)&np.isclose(idx['size'],size)&np.isclose(idx.fraction,f)].sort_values('replicate')
                # Every complete retained seed enters. No coherence, contrast or region-identity gate.
                hist=np.zeros((3,600),dtype=np.int64);bounds=size*1e-4*np.sinh(np.linspace(-np.arcsinh(20000),np.arcsinh(20000),601))
                coherences=[];circumferences=[];counts=np.zeros(3,dtype=np.int64);nonfinite=0;switches=0;exposure=0.
                for r in q.itertuples(index=False):
                    with np.load(r.compact,allow_pickle=False) as z:
                        t=z['time'].astype(float);circ=z['circumference'].astype(float)
                        P=z['P_real'].astype(float)+1j*z['P_imag'].astype(float)
                        coherence=np.abs(P);theta=np.angle(P);ids=z['dominant_track_id'];cad=float(np.median(np.diff(t)))
                        coherences.append(coherence);circumferences.append(circ)
                        # Raw detector-based switches are reported only as a lag diagnostic; they do not gate the histograms.
                        good=(ids[:-1]>=0)&(ids[1:]>=0)
                        switches+=int(np.count_nonzero(good&(ids[1:]!=ids[:-1])))
                        exposure+=float(t[-1]-t[0])
                        for k,mult in enumerate((1,4,16)):
                            step=int(round(mult*tau/cad)); angle=np.angle(np.exp(1j*(theta[step:]-theta[:-step])))
                            displacement=.5*(circ[step:]+circ[:-step])*angle/(2*np.pi)
                            valid=np.isfinite(displacement);nonfinite+=int(np.count_nonzero(~valid))
                            hist[k]+=np.histogram(displacement[valid],bounds)[0];counts[k]+=int(valid.sum())
                    sources.append(dict(model=model,f=f,size=size,replicate=int(r.replicate),compact=rel(r.compact),compact_sha256=r.compact_sha256))
                coh=np.concatenate(coherences);circ=np.concatenate(circumferences);quant=np.quantile(coh,[0,.01,.05,.25,.5,.75,.95,.99,1])
                coh_rows.append(dict(model=model,f=f,size=size,seeds=len(q),frames=len(coh),mean=float(np.mean(coh)),sd=float(np.std(coh)),minimum=float(quant[0]),q01=float(quant[1]),q05=float(quant[2]),q25=float(quant[3]),median=float(quant[4]),q75=float(quant[5]),q95=float(quant[6]),q99=float(quant[7]),maximum=float(quant[8]),fraction_below_0_05=float(np.mean(coh<.05)),fraction_below_0_1=float(np.mean(coh<.1)),nonfinite_lag_pairs=nonfinite,mean_circumference=float(np.mean(circ))))
                mean_wait=exposure/switches if switches else None
                # Counted switches are detector-dependent and can include patch/track fragmentation.
                waiting.append(dict(model=model,f=f,size=size,tau=tau,lag16=16*tau,raw_track_switches=switches,exposure=exposure,raw_track_switch_rate=switches/exposure,raw_track_mean_wait=mean_wait,lag16_over_mean_wait=16*tau/mean_wait if mean_wait else 0.,lag16_not_well_below_wait=bool(mean_wait and 16*tau>=.1*mean_wait)))
                for k,mult in enumerate((1,4,16)):
                    lag=mult*tau;width=np.diff(bounds);centres=.5*(bounds[1:]+bounds[:-1]);density=hist[k]/counts[k]/width
                    for i in range(len(centres)):
                        bins_rows.append(dict(model=model,f=f,size=size,tau=tau,lag=lag,lag_over_tau=mult,physical_bin_left=float(bounds[i]),physical_bin_right=float(bounds[i+1]),physical_bin_center=float(centres[i]),rescaled_displacement=float(centres[i]/np.sqrt(lag)),rescaled_density=float(density[i]*np.sqrt(lag)),count=int(hist[k,i]),pairs=int(counts[k]),seeds=len(q)))
                print(model,f,size,len(q),'complete',flush=True)
    pd.DataFrame(bins_rows).to_csv(data/'D_raw_heading_distributions.csv.gz',index=False)
    pd.DataFrame(coh_rows).to_csv(data/'D_coherence_distributions.csv',index=False)
    pd.DataFrame(waiting).to_csv(data/'D_relocation_lag_diagnostics.csv',index=False)
    pd.DataFrame(sources).to_csv(meta/'D_compact_sources.csv',index=False)
    write(meta/'D_analysis_provenance.json',dict(generated_utc=now(),input_sha256=SOURCES,source_count=len(sources),coherence_gate='none',contrast_gate='none',heading='arg(P), signed principal difference in (-pi,pi]',physical_displacement='mean endpoint circumference × principal heading difference / (2pi)',histogram='600 symmetric asinh-spaced physical displacement bins from -2 nominal C to +2 nominal C (covers instantaneous Spring circumference), asinh scale C*1e-4; transform x by 1/sqrt(lag), density by sqrt(lag)',lags='tau,4tau,16tau',weak_inputs=WEAK,strong_input=.3,selection='all retained complete seeds; no seed selection',nonfinite='only nonfinite numeric differences excluded; count retained per condition'))

