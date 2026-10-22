"""Weak-global-vector context for the existing endpoint-concurrency exposure.

Adds no exclusion, threshold change or event classifier. The original endpoint
source, protocol and counts remain unchanged.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np
import pandas as pd
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
FIG=ROOT/'polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4'
OUT=FIG/'endpoint_concurrency_20261019';DATA=OUT/'supporting/data'
RUN=ROOT/'polarity1d/manuscript/data/runs/revision_noise_strength_replication_scoped_20261019'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    coverage=pd.read_csv(DATA/'case_exposure_coverage.csv').set_index('case_id')
    assert coverage.undefined_initial_centres.sum()==0
    geometry_path=FIG/'signature_assembly_20261019/supporting/data/postburn_region_samples.npz'
    cached=json.loads((RUN/'results/partial/case_cache.json').read_text());rows=[]
    with np.load(geometry_path) as geometry:
        for case_id,c in cached.items():
            d=c['descriptors'][0];t=geometry[case_id+'__time'];qual=geometry[case_id+'__qualified'];regions=geometry[case_id+'__regions']
            keep=(t[:-2]>=t[-1]/2)&qual[:-2]&(regions[:-2]==1)
            assert int(keep.sum())==int(coverage.loc[case_id,'eligible_initial_frames'])
            assert sha(c['trace'])==c['trace_sha256']
            with np.load(c['trace']) as z:
                ix=np.searchsorted(z['time'],t);np.testing.assert_array_equal(z['time'][ix],t)
                R=np.hypot(z['P_real'][ix],z['P_imag'][ix])
            n=int(keep.sum());first=R[:-2][keep];last=R[2:][keep]
            rows.append(dict(case_id=case_id,**{k:d[k] for k in ('model','size','fraction','kind','seed')},eligible_initial_frames=n,
                initial_weak_R_frames=int((first<.1).sum()),endpoint_weak_R_frames=int((last<.1).sum()),
                initial_weak_R_fraction=float(np.mean(first<.1)) if n else np.nan,
                endpoint_weak_R_fraction=float(np.mean(last<.1)) if n else np.nan,
                trace=c['trace'],trace_sha256=c['trace_sha256']))
    pd.DataFrame(rows).to_csv(DATA/'case_global_vector_context.csv',index=False)
    receipt=dict(generated_at=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),cases=len(rows),
        native_integrations=0,global_R_used_as_exclusion=False,
        original_endpoint_counts_and_protocol_unchanged=True,
        eligibility_reconciliation='All936 exact; no undefined initial centres, so saved qualified single-region masks exactly reproduce existing exposure',
        source_geometry_sha256=sha(geometry_path),output_sha256=sha(DATA/'case_global_vector_context.csv'))
    (OUT/'supporting/metadata/GLOBAL_VECTOR_CONTEXT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
