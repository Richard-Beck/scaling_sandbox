"""Edges from the existing bounded float16 Figure 4 exemplars; no new noise."""
from pathlib import Path
import os,sys
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
sys.dont_write_bytecode=True
import numpy as np,pandas as pd
import figure3_round5_edge_20261023 as edge
OUT=edge.OUT; DATA=OUT/'supporting/data'; META=OUT/'supporting/metadata'
DRIFT=edge.MS/'figures/length_mechanisms_20261008/dev/revision_round4_20261023/figure4/supporting/data/drift_decomposition.csv'

def profile_edges(p,C):
    ans=edge.edges(p,C)
    if ans is None:return None,np.nan
    # Bound from half-ULP errors in two crossing samples and profile extrema;
    # conservative positional uncertainty, not a fitted measurement variance.
    error=.5*np.max(abs(np.spacing(np.asarray(p,dtype=np.float16)).astype(float)))
    level=(p.min()+p.max())/2;n=len(p);high=p>=level
    j=np.flatnonzero(high!=np.roll(high,1));slopes=abs(p[j]-p[(j-1)%n])/(C/n)
    return ans,float(np.max(2*error/np.maximum(slopes,1e-20)))

def main():
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('SLURM required')
    DATA.mkdir(parents=True,exist_ok=True);META.mkdir(parents=True,exist_ok=True)
    index=pd.read_csv(edge.ORIGINAL/'results/compact_index.csv'); index=index[(index.replicate.isin([0,1]))&(index.fraction<=.01)].copy()
    drift=pd.read_csv(DRIFT);drift=drift[(drift.flicker_frames==5)&(drift.MSD_range_fraction==.25)]
    summaries=[];curves=[];coverage=[]
    manifest=edge.read(edge.ORIGINAL/'manifest.json');basals={(b['model'],b['size']):b for b in manifest['basals']}
    for row in index.itertuples():
        path=Path(row.exemplar)
        if not path.is_file():coverage.append(dict(condition_id=row.condition_id,replicate=row.replicate,status='missing_exemplar'));continue
        with np.load(path) as z:
            t=z['time'].astype(float); C=z['circumference'].astype(float);p=z['primary'].astype(float)
        basal=basals[row.model,row.size]
        with np.load(edge.ORIGINAL/basal['input']) as z:
            state=z['state'];bp=state[:basal['cells']] if row.model=='spring_ring' else state[0]
        contrast_cutoff=.1*float(np.ptp(bp))
        Cref=float(C.mean()); positions=np.full((len(t),2),np.nan);error=np.full(len(t),np.nan)
        for i in range(len(t)):
            if np.ptp(p[i])<contrast_cutoff:continue
            e,q=profile_edges(p[i],C[i]);error[i]=q
            if e is not None:positions[i]=e/C[i]*Cref
        valid=np.isfinite(positions).all(axis=1); segments=np.cumsum(~valid)
        # Unwrap each uninterrupted unique-domain segment separately.
        for s in np.unique(segments[valid]):
            ids=np.flatnonzero(valid&(segments==s));positions[ids]=np.unwrap(positions[ids]*2*np.pi/Cref,axis=0)*Cref/(2*np.pi)
        cadence=float(np.median(np.diff(t)));maxlag=max(1,int(len(t)/4));lags=np.unique(np.r_[1,np.round(np.geomspace(1,maxlag,32)).astype(int)])
        model_drift=drift[(drift.model==row.model)&(drift.C_label==row.size)&(drift.f==row.fraction)]
        target=float(model_drift.iloc[0].analysis_lag) if len(model_drift) else np.nan
        # analysis_lag is sometimes unavailable because the decay rate is not
        # resolved; edge diagnostics remain defined at tau and 1/4 horizon.
        tau=5. if row.model=='spring_ring' else 50.
        lags=np.unique(np.r_[lags,max(1,int(round(tau/cadence)))])
        if np.isfinite(target):lags=np.unique(np.r_[lags,min(maxlag,max(1,int(round(target/cadence))))])
        records=[]
        for lag in lags:
            ok=valid[:-lag]&valid[lag:]&(segments[:-lag]==segments[lag:]); count=int(ok.sum())
            if count<5:continue
            dis=positions[lag:][ok]-positions[:-lag][ok];left,right=dis.T
            cov=float(np.mean((left-left.mean())*(right-right.mean())));vari=np.var(dis,axis=0)
            rho=cov/np.sqrt(vari.prod()) if min(vari)>0 else np.nan
            dl,dr=np.mean(dis*dis,axis=0);cross=float(np.mean(left*right));dt=float(np.mean((t[lag:]-t[:-lag])[ok]))
            bound=float(np.nanmean((error[:-lag][ok]+error[lag:][ok])**2))
            rec=dict(condition_id=row.condition_id,model=row.model,size=row.size,fraction=row.fraction,replicate=row.replicate,lag=dt,pairs=count,left_MSD=float(dl),right_MSD=float(dr),edge_cross_product=cross,increment_correlation=rho,domain_MSD=(dl+dr+2*cross)/4,independent_domain_MSD=(dl+dr)/4,correlation_factor=(dl+dr+2*cross)/(dl+dr) if dl+dr>0 else np.nan,quantization_displacement_bound_squared=bound,resolved_above_quantization=bool(min(dl,dr)>10*bound),unique_domain_frame_fraction=float(valid.mean()),exemplar_bins=p.shape[1],exemplar_cadence=cadence)
            records.append(rec);curves.append(rec)
        if records:
            selected=min(records,key=lambda r:abs(r['lag']-tau))
            late=max(records,key=lambda r:r['lag'])
            summaries.append(dict(selected,lag_kind='nearest_tau',analysis_lag_requested=target))
            summaries.append(dict(late,lag_kind='quarter_horizon',analysis_lag_requested=target))
        coverage.append(dict(condition_id=row.condition_id,replicate=row.replicate,status='measured' if records else 'no_unique_domain_pairs',unique_domain_frame_fraction=float(valid.mean()),exemplar=str(path),exemplar_sha256=edge.sha(path),rows=len(records)))
    pd.DataFrame(summaries).to_csv(DATA/'edge_independence_summary.csv',index=False)
    pd.DataFrame(curves).to_csv(DATA/'edge_MSD_curves.csv',index=False)
    pd.DataFrame(coverage).to_csv(DATA/'edge_independence_coverage.csv',index=False)
    edge.write(META/'edge_independence_receipt.json',dict(generated_utc=edge.now(),producer='COMPLETE',n_exemplars=len(index),measured=sum(c['status']=='measured' for c in coverage),producer_sha256=edge.sha(Path(__file__)),input_index_sha256=edge.sha(edge.ORIGINAL/'results/compact_index.csv'),quantization='Float16 profile values; conservative crossing bound from half ULP and local crossing slope. Low-amplitude short-lag correlations may be below saved-profile precision.',selection='Both edges only in uninterrupted frames with exactly one current half-contrast region and contrast >= 0.1 native basal contrast; no interpolation across multiple-domain or homogeneous frames.',spring_coordinate='Fractional material-edge displacement times mean postburn physical circumference; expansion excluded.'))
    print('complete',len(summaries),len(curves))
if __name__=='__main__':main()
