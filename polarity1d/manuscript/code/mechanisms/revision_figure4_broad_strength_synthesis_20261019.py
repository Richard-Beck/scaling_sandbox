"""Read-only memory/organization synthesis and an unsubmitted replication proposal."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[4]
MS=ROOT/'polarity1d/manuscript'
RUN=MS/'data/runs/revision_noise_strength_broad_20261019'
OUT=MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/broad_strength_campaign_20261019/review_complete_scientific'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    partial=RUN/'results/partial';data=OUT/'supporting/data';meta=OUT/'supporting/metadata'
    cond=pd.read_csv(data/'condition_window_snapshot.csv')
    sources={}
    for name in ('case_window_memory.csv','case_window_descriptors.csv','worker_inventory.csv'):
        raw=(partial/name).read_bytes();(data/name).write_bytes(raw);sources[name]=hashlib.sha256(raw).hexdigest()
    memory=pd.read_csv(data/'case_window_memory.csv');cases=pd.read_csv(data/'case_window_descriptors.csv')
    assert memory.case_id.nunique()==cases.case_id.nunique()==576
    rr=[]
    for key,g in memory.groupby(['model','size','fraction','kind','window','lag'],sort=True):
        count=int(g.pairs_eligible.sum());total=int(g.pairs_available.sum())
        corr=float(g.direction_cosine_sum.sum()/count) if count else np.nan
        supported=(count>=100 and count/total>=.01)
        energ=float(np.sqrt(g.energy_a.sum()*g.energy_b.sum()))
        seedvals=g[g.support].direction_correlation
        rr.append(dict(zip(['model','size','fraction','kind','window','lag'],key),cases=len(g),
            direction_pairs=count,all_pairs=total,direction_coverage=count/total,seeds_with_pairs=int(np.sum(g.pairs_eligible>0)),
            seeds_with_declared_support=int(g.support.sum()),group_support=supported,
            direction_cosine_pooled=corr,direction_cosine_reported=corr if supported else np.nan,
            supported_seed_cosine_min=float(seedvals.min()) if len(seedvals) else np.nan,
            supported_seed_cosine_max=float(seedvals.max()) if len(seedvals) else np.nan,
            raw_vector_correlation=float(g.raw_product_sum.sum()/energ) if energ else np.nan))
    aggregate=pd.DataFrame(rr)
    aggregate.to_csv(data/'memory_all_conditions_windows_lags.csv',index=False)
    terminal=aggregate[(aggregate.kind=='noise')&(((aggregate.model=='spring_ring')&(aggregate.lag==800))|((aggregate.model!='spring_ring')&(aggregate.lag==8000)))]
    joint=cond[cond.kind=='noise'].merge(terminal,on=['model','size','fraction','kind','window'],validate='one_to_one')
    joint.to_csv(data/'organization_memory_joint_windows.csv',index=False)
    full=joint[joint.window=='full'].copy()
    for win in ('first_half','last_half'):
        w=joint[joint.window==win]
        full=full.merge(w[['model','size','fraction','qualified_multiregion_fraction','direction_cosine_reported','direction_coverage','low_contrast_fraction']],on=['model','size','fraction'],suffixes=('',f'_{win}'),validate='one_to_one')
    full['halfwindow_occupancy_difference']=abs(full.qualified_multiregion_fraction_first_half-full.qualified_multiregion_fraction_last_half)
    full['halfwindow_direction_difference']=abs(full.direction_cosine_reported_first_half-full.direction_cosine_reported_last_half)
    full['halfwindow_coverage_difference']=abs(full.direction_coverage_first_half-full.direction_coverage_last_half)
    full.to_csv(data/'organization_memory_seed_window_sensitivity.csv',index=False)
    longest=full[((full.model=='spring_ring')&(full['size']==3))|((full.model!='spring_ring')&(full['size']==60))]
    longest.to_csv(data/'largest_size_memory_context.csv',index=False)
    episode=[]
    for key,g in cases[cases.kind=='noise'].groupby(['model','size','fraction','window']):
        q=g.qualified_multiregion_run_q90.dropna();span=g.qualified_multiregion_longest_observed_span.dropna()
        episode.append(dict(zip(['model','size','fraction','window'],key),cases=len(g),cases_with_sampled_positive_runs=len(q),
            q90_seed_min=float(q.min()) if len(q) else np.nan,q90_seed_max=float(q.max()) if len(q) else np.nan,
            longest_span_seed_min=float(span.min()) if len(span) else np.nan,longest_span_seed_max=float(span.max()) if len(span) else np.nan,
            cases_left_window_censored=int(g.qualified_multiregion_left_window_censored.sum()),cases_right_window_censored=int(g.qualified_multiregion_right_window_censored.sum()),
            scalar_cadence=float(g.scalar_cadence.iloc[0])))
    pd.DataFrame(episode).to_csv(data/'sampled_qualified_run_spans.csv',index=False)
    manifest=json.loads((RUN/'manifest.json').read_text());proposal=[];planned=[]
    for model in ('wave_pinning','goryachev','holmes_model3','spring_ring'):
        for size in ((1.,2.,3.) if model=='spring_ring' else (20.,40.,60.)):
            for f in (.003,.03,.1):
                jobs=[j for j in manifest['jobs'] if j['model']==model and j['size']==size and j['fraction']==f and j['kind']=='noise' and j['variant']=='standard']
                assert len(jobs)==2 and all(not j['dense_observation'] for j in jobs)
                elapsed=[];allocated=[]
                for j in jobs:
                    folder=RUN/'workers'/f"{j['slot']:04d}"
                    status=json.loads((folder/'status.json').read_text());elapsed.append(status['elapsed_seconds'])
                    allocated.append(sum(p.stat().st_blocks*512 for p in folder.rglob('*') if p.is_file()))
                row=dict(model=model,size=size,fraction=f,existing_seeds=2,target_seeds=12,additional_workers=10,
                    existing_elapsed_min=min(elapsed),existing_elapsed_max=max(elapsed),existing_elapsed_mean=float(np.mean(elapsed)),
                    estimated_worker_hours=10*float(np.mean(elapsed))/3600,existing_allocated_mean=float(np.mean(allocated)),
                    estimated_additional_allocated_bytes=10*float(np.mean(allocated)),input=jobs[0]['input'],input_sha256=jobs[0]['input_sha256'])
                proposal.append(row)
                first=min(j['seed'] for j in jobs)
                assert sorted(j['seed'] for j in jobs)==[first,first+1]
                for rep in range(2,12):planned.append(dict(model=model,size=size,fraction=f,replicate=rep,seed=first+rep,input=jobs[0]['input'],input_sha256=jobs[0]['input_sha256'],submitted=False))
    pr=pd.DataFrame(proposal);pr.to_csv(data/'replication_proposal_conditions.csv',index=False)
    pd.DataFrame(planned).to_csv(data/'replication_proposal_unsubmitted_cases.csv',index=False)
    summary=dict(generated_at=datetime.now(timezone.utc).isoformat(),submitted=False,source_sha256=sha(__file__),inputs=sources,
        condition_snapshot_sha256=sha(data/'condition_window_snapshot.csv'),models=4,sizes_per_model=3,common_strengths=[.003,.03,.1],
        condition_groups=36,additional_noisy_workers=360,target_seeds=12,additional_mean_controls=0,additional_basals=0,
        estimated_worker_hours=float(pr.estimated_worker_hours.sum()),measured_input_worker_max_seconds=float(pr.existing_elapsed_max.max()),
        estimated_additional_allocated_bytes=int(pr.estimated_additional_allocated_bytes.sum()),storage_reserve_30pct_bytes=int(pr.estimated_additional_allocated_bytes.sum()*1.3),
        proposed_resources=dict(qos='small',cpus=1,memory_GiB=3,walltime='02:00:00',array_throttle=None),
        qualification='Proposal only; requires Spring numerical/observation acceptance and root/independent preflight; author steering is welcome without an additional authorization gate; no job submitted',
        limits='Timing/output estimates use two actual existing cases per group, not upper prediction bounds; retain half-window sensitivity and n2 endpoint uncertainty')
    (meta/'REPLICATION_PROPOSAL.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
