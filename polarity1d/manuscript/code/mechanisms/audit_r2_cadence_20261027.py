"""Additive post-launch audit of frozen R2 mode cadences; no simulations.

The prelaunch choice and manifests stay immutable. This independent audit uses
the original full-record mode-rate estimator and exact four-seed bootstrap.
"""
from pathlib import Path
import argparse, json, os, sys
from datetime import datetime, timezone
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import ring_r2_campaign_20261027 as C

def now():return datetime.now(timezone.utc).isoformat()
def audit(out):
    out=Path(out).resolve()
    if out.exists():raise ValueError('Choose a fresh additive audit directory')
    out.mkdir(parents=True)
    original=C.legacy();engine=original.engine();manifest=C.read(C.SOURCE/'manifest.json')
    selection=C.read(C.RUN/'validation/selection.json');pre=pd.read_csv(C.RUN/'validation/cadence.csv')
    preconditions=set(pre.condition.astype(int))
    levels=C.read(C.SOURCE/'expected_levels.json');cohorts=[];pending=[]
    for condition in manifest['conditions']:
        if condition['model'] not in ('wave_pinning','spring_ring'):continue
        expected=[int(levels[str(s)]) for s in condition['slots']]
        if len(set(expected))!=1:raise ValueError('Mixed expected levels in cohort')
        level=expected[0]
        statuses=[C.segment(C.SOURCE,s,level)/'status.json' for s in condition['slots']]
        if all(p.exists() and C.read(p).get('complete') for p in statuses):
            cohorts.append(dict(condition=condition,level=level,statuses=[str(p.relative_to(C.ROOT)) for p in statuses],
                                status_sha256=[C.sha(p) for p in statuses]))
        else:pending.append(dict(condition=condition['index'],model=condition['model'],size=condition['size'],rate=condition['rate_level'],amplitude=condition['amplitude_level'],expected_level=level))
    C.write(out/'snapshot.json',dict(generated_utc=now(),expected_levels=levels,cohorts=cohorts,pending=pending,
        source_manifest_sha256=C.sha(C.SOURCE/'manifest.json'),selection_sha256=C.sha(C.RUN/'validation/selection.json'),
        prelaunch_cadence_sha256=C.sha(C.RUN/'validation/cadence.csv'),code_sha256=C.sha(__file__),
        engine_sha256=C.sha(C.SOURCE/'frozen_sources/engine.py'),production_sha256=C.sha(C.SOURCE/'frozen_sources/production.py')))
    rows=[]
    for entry in cohorts:
        condition=entry['condition'];level=entry['level'];job=manifest['jobs'][condition['slots'][0]]
        factor=int(selection['factors'][condition['model']]);cadence=job['frame_interval']
        histories=[C.series(C.SOURCE,s,level) for s in condition['slots']]
        raw=np.asarray([engine.acf_raw(x) for x in histories]);full=raw.mean(0)/raw[:,0].mean()
        gamma,status,bound=engine.half_rate(full,cadence)
        low,high,low_status,high_status,p_resolved=original.bootstrap_rates(raw,cadence)
        reduced_raw=np.asarray([engine.acf_raw(x[::factor]) for x in histories]);reduced=reduced_raw.mean(0)/reduced_raw[:,0].mean()
        reduced_gamma,reduced_status,reduced_bound=engine.half_rate(reduced,factor*cadence)
        value=reduced_gamma if reduced_gamma is not None else reduced_bound
        same_status=reduced_status==status
        inside_ci=bool(low<=value<=high)
        bound_preserved=bool(np.isclose(reduced_bound,bound,rtol=.002)) if status!='resolved' and same_status else status=='resolved'
        passed=bool(same_status and (inside_ci if status=='resolved' else bound_preserved))
        rows.append(dict(condition=condition['index'],model=condition['model'],size=condition['size'],rate=condition['rate_level'],
            amplitude=condition['amplitude_level'],extension_level=level,new_since_prelaunch=condition['index'] not in preconditions,
            seeds=len(condition['slots']),factor=factor,mode_cutoff=1,full_cadence=cadence,reduced_cadence=factor*cadence,
            full_gamma=gamma,full_status=status,full_bound=bound,full_ci_low=low,full_ci_high=high,
            full_ci_low_status=low_status,full_ci_high_status=high_status,full_bootstrap_resolved_fraction=p_resolved,
            reduced_gamma=reduced_gamma,reduced_status=reduced_status,reduced_bound=reduced_bound,
            full_analysis_lag=(len(full)-1)*cadence,reduced_analysis_lag=(len(reduced)-1)*factor*cadence,
            same_status=same_status,reduced_inside_full_ci=inside_ci,bound_preserved=bound_preserved,passes=passed))
        C.table(out/'cadence.csv',rows)
        print(json.dumps(rows[-1]),flush=True)
    C.write(out/'summary.json',dict(generated_utc=now(),slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        complete=True,cohorts_checked=len(rows),prelaunch_cohorts=len(preconditions),new_cohorts=sum(r['new_since_prelaunch'] for r in rows),
        prelaunch_conditions_not_currently_complete=sorted(preconditions-{r['condition'] for r in rows}),
        pending_cohorts=pending,passes=sum(r['passes'] for r in rows),failures=[r for r in rows if not r['passes']],
        selected_factors=selection['factors'],criterion='Same original half-rate censor status; resolved reduced gamma inside full-record exact four-seed 95% bootstrap CI; noncrossing rate bound preserved to prelaunch rtol .002.',
        frozen_selection_unchanged=True,simulations_performed=False,source_records_read_only=True,
        selection_sha256=C.sha(C.RUN/'validation/selection.json'),prelaunch_cadence_sha256=C.sha(C.RUN/'validation/cadence.csv')))
    print(json.dumps(C.read(out/'summary.json')),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();audit(a.output)
