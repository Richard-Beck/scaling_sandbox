"""Describe existing raw Figure5D histograms without classification or gating.

Quantiles and fixed-distance tail weights integrate the displayed physical-bin
histograms using a uniform density inside each original bin. This is a compact
descriptive reanalysis, not a new simulation or a drift/relocation pass test.
"""
from pathlib import Path
import csv,gzip,hashlib,json
from datetime import datetime,timezone
import numpy as np
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round7_20261026/figure5'
SOURCE=FIG/'supporting/data/D_raw_heading_distributions.csv.gz'
DEST=FIG/'supporting/data/D_descriptive_lag_quantiles.csv'

def main():
    with gzip.open(SOURCE,'rt') as s:raw=list(csv.DictReader(s))
    groups={}
    for row in raw:
        key=tuple(row[k] for k in ('model','f','size','tau','lag','lag_over_tau','pairs','seeds'))
        groups.setdefault(key,[]).append(row)
    results=[]
    for key,rows in groups.items():
        model,f,size,tau,lag,mult,pairs,seeds=key
        left=np.array([float(r['physical_bin_left']) for r in rows]);right=np.array([float(r['physical_bin_right']) for r in rows]);count=np.array([int(r['count']) for r in rows]);width=right-left
        if int(count.sum())!=int(pairs):raise ValueError('Incomplete raw histogram')
        def absolute_cdf(radius):
            overlap=np.maximum(0.,np.minimum(right,radius)-np.maximum(left,-radius))
            return float(np.sum(count*overlap/width)/count.sum())
        def quantile(p):
            lo=0.;hi=float(np.max(np.maximum(abs(left),abs(right))))
            for _ in range(64):
                mid=(lo+hi)/2
                if absolute_cdf(mid)>=p:hi=mid
                else:lo=mid
            return (lo+hi)/2
        median=quantile(.5);q95=quantile(.95);threshold=.15 if model=='spring_ring' else 1.5
        results.append(dict(model=model,f=float(f),size=float(size),tau=float(tau),lag=float(lag),lag_over_tau=int(mult),
            pairs=int(pairs),seeds=int(seeds),physical_abs_median=median,physical_abs_q95=q95,
            rescaled_abs_median=median/np.sqrt(float(lag)),rescaled_abs_q95=q95/np.sqrt(float(lag)),
            fixed_physical_tail_threshold=threshold,tail_probability_beyond_fixed_physical_displacement=1-absolute_cdf(threshold),
            method='original-bin uniform-density integration; no coherence, contrast, or track-identity cut'))
    results.sort(key=lambda r:(r['model'],r['f'],r['size'],r['lag_over_tau']))
    with DEST.open('w',newline='') as s:w=csv.DictWriter(s,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    receipt=dict(generated_utc=datetime.now(timezone.utc).isoformat(),input_path=str(SOURCE.relative_to(ROOT)),
        input_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),output_path=str(DEST.relative_to(ROOT)),
        output_sha256=hashlib.sha256(DEST.read_bytes()).hexdigest(),histograms=len(results),all_original_counts_retained=True,
        estimator='Integrate piecewise-uniform density inside original physical bins; absolute displacement median and95th percentile, then divide by square root of physical lag. Tail mass uses fixedphysical displacement1.5(Spring0.15).',
        interpretation='Descriptive summaries only; no passcriterion, no coherence cut, no claim tail implies a tracked relocation; overlapping lags are not independent observations.')
    (FIG/'supporting/metadata/D_descriptive_lag_quantiles_provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
    for r in results:
        print(f"{r['model']:14s} f{r['f']:g} C{r['size']:g} {r['lag_over_tau']:2d}τ median{r['rescaled_abs_median']:.5g} q95{r['rescaled_abs_q95']:.5g} tail{r['tail_probability_beyond_fixed_physical_displacement']:.5g}")

if __name__=='__main__':main()
