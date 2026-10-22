"""Exploratory fixed-lag geographic concurrency from existing F4 profiles.

No integration, lineage, birth, establishment, lifetime or event-rate classifier.
The main signatures do not depend on this optional derived supplement.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS = ROOT/'polarity1d/manuscript'
RUN = MS/'data/runs/revision_noise_strength_replication_scoped_20261019'
PILOT = MS/'data/runs/revision_noise_strength_broad_20261019'
OUT = MS/'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/endpoint_concurrency_20261019'
PROTOCOL = dict(
    window='Second half; both endpoints within [T/2,T]',
    lag=dict(chemical=80., spring=8.),
    native_profile_cadence=dict(chemical=40., spring=4.),
    physical_distance_bin_width=dict(chemical=2., spring=.1),
    region='Periodic connected cells at or above current min+0.5*(max-min); flat contrast<=1e-12 gives zero',
    amplitude_qualification='Current contrast >=0.1 of initial basal contrast, at each endpoint',
    initial_exposure='Amplitude-qualified initial frame with exactly one half-contrast region and defined local circular centre',
    centre='Circular centre weighted by primary-current minimum inside the connected region; valid if local phasor magnitude/weight sum >1e-6',
    outcome='At endpoint, amplitude qualification, at least one region overlapping the initial region mask, and at least one separate nonoverlapping region with defined centre in the distance bin',
    probability='Per initial frame and distance bin: at least one such endpoint region; several bins may co-occur. Never component count divided by frames.',
    bin_exposure='Initial eligible frame with a positive length of initial exterior cells in that distance bin; all such frames remain in denominator, irrespective of endpoint outcome',
    distance='Shortest periodic centre separation; actual initial-frame C times wrapped material-coordinate difference; initial Spring C used, endpoint stretch exported separately',
    undefined='Undefined initial centre excludes exposure but remains in coverage; undefined secondary centres remain explicit unresolved outcomes with lower/upper probability bounds',
    weighting='Per-case frame probabilities; condition means give each contributing independent case equal weight, with observed seed ranges and explicit contributor counts',
    cadence_sensitivity='Four largest-size f=1 standard paths: same physical lag on each of ten decimations of saved dense profiles. These are descriptive phase-sampling checks, not an extra numerical acceptance gate.',
    scope='Exploratory geographic endpoint concurrency only. Geographic mask overlap is not incumbent identity, birth, survival, continuous coexistence or establishment.',
    original_Spring_run_q90_gate='Failed and unchanged; no continuous lifetime or event rate is inferred')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
    if isinstance(x,np.generic):return x.item()
    return x


def write(path,x):Path(path).write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n')


def components(mask):
    if not mask.any():return []
    if mask.all():return [np.arange(len(mask))]
    starts=np.flatnonzero(mask & ~np.roll(mask,1))
    result=[]
    for start in starts:
        n=0
        while mask[(start+n)%len(mask)]:n+=1
        result.append((start+np.arange(n))%len(mask))
    return result


def local_centre(cells,u,q):
    weights=u[cells]-u.min()
    total=float(weights.sum())
    z=np.sum(weights*np.exp(2j*np.pi*q[cells]))
    coherence=float(abs(z)/total) if total>0 else 0.
    return ((float(np.angle(z)/(2*np.pi))%1.) if coherence>1e-6 else np.nan),coherence


def read_profile(path):
    with np.load(path) as z:
        t=z['time'].copy();u=z['primary'].copy()
        C=np.broadcast_to(z['circumference'],len(t)).copy()
        q=z['material_positions'].copy() if 'material_positions' in z else z['x']/C[0]
    return t,u,C,q


def assess(path,platform,edges,phase=0,stride=1):
    t,u,C,q=read_profile(path)
    dt=float(t[1]-t[0]);lag=PROTOCOL['lag'][platform]
    offset=int(round(lag/dt));assert abs(offset*dt-lag)<1e-8
    contrast=np.ptp(u,axis=1);qualified=(contrast>1e-12)&(contrast>=.1*contrast[0])
    masks=u>=u.min(axis=1)[:,None]+.5*contrast[:,None]
    groups=[components(mask) if contrast[i]>1e-12 else [] for i,mask in enumerate(masks)]
    start_indices=np.arange(phase,len(t)-offset,stride)
    start_indices=start_indices[(t[start_indices]>=t[-1]/2)&(t[start_indices+offset]<=t[-1])]
    bins=len(edges)-1
    exposure=np.zeros(bins,dtype=int);events=exposure.copy();unresolved=exposure.copy()
    opportunity=np.zeros(bins);follow_qualified=exposure.copy();follow_overlap=exposure.copy()
    all_pairs=len(start_indices);single=eligible=undefined=geographic=0;maxstretch=0.
    initial_local=[];stretches=[]
    for i in start_indices:
        if not qualified[i] or len(groups[i])!=1:continue
        single+=1
        initial=masks[i];centre,coherence=local_centre(groups[i][0],u[i],q)
        if not np.isfinite(centre):undefined+=1;continue
        eligible+=1;initial_local.append(coherence)
        cell_distance=C[i]*abs((q-centre+.5)%1-.5)
        outside=~initial
        arc=np.histogram(cell_distance[outside],bins=edges)[0]*C[i]/len(q)
        risk=arc>0
        exposure+=risk;opportunity+=arc
        j=i+offset;stretch=(C[j]-C[i])/C[i];stretches.append(stretch)
        maxstretch=max(maxstretch,abs(stretch))
        if not qualified[j]:continue
        follow_qualified+=risk
        overlap=[np.any(initial[cells]) for cells in groups[j]]
        if not any(overlap):continue
        follow_overlap+=risk;geographic+=1
        secondary=[cells for cells,hit in zip(groups[j],overlap) if not hit]
        distances=[];bad=False
        for cells in secondary:
            end,quality=local_centre(cells,u[j],q)
            if not np.isfinite(end):bad=True;continue
            distances.append(C[i]*abs((end-centre+.5)%1-.5))
        hits=np.histogram(distances,bins=edges)[0]>0
        events+=hits & risk
        if bad:unresolved+=risk & ~hits
    rows=[]
    for b in range(bins):
        n=int(exposure[b]);k=int(events[b]);missing=int(unresolved[b])
        rows.append(dict(distance_lo=edges[b],distance_hi=edges[b+1],bin_exposure=n,
            event_frames=k,unresolved_secondary_frames=missing,
            endpoint_qualified_frames=int(follow_qualified[b]),endpoint_overlap_frames=int(follow_overlap[b]),
            endpoint_qualified_coverage=float(follow_qualified[b]/n) if n else np.nan,
            endpoint_overlap_coverage=float(follow_overlap[b]/n) if n else np.nan,
            probability_lower=float(k/n) if n else np.nan,probability_upper=float((k+missing)/n) if n else np.nan,
            available_exterior_arc_sum=float(opportunity[b]),
            mean_available_exterior_arc=float(opportunity[b]/n) if n else np.nan))
    coverage=dict(profile_cadence=dt,sampling_stride=stride,sampling_phase=phase,lag=lag,
        potential_initial_frames=all_pairs,qualified_single_region_frames=single,
        undefined_initial_centres=undefined,eligible_initial_frames=eligible,
        eligible_fraction=float(eligible/all_pairs),
        endpoint_qualified_geographic_overlap_frames=geographic,
        minimum_initial_local_coherence=min(initial_local) if initial_local else np.nan,
        mean_relative_C_change=float(np.mean(stretches)) if stretches else np.nan,
        maximum_absolute_relative_C_change=maxstretch if stretches else np.nan)
    return rows,coverage


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,default=OUT)
    p.add_argument('--freeze-only',action='store_true');a=p.parse_args();out=a.out
    data=out/'supporting/data';metadata=out/'supporting/metadata'
    for d in (data,metadata):d.mkdir(parents=True,exist_ok=True)
    protocol=metadata/'PROTOCOL.json'
    if protocol.exists():assert json.loads(protocol.read_text())['definitions']==PROTOCOL
    else:write(protocol,dict(frozen_at=datetime.now(timezone.utc).isoformat(),definitions=PROTOCOL,source_sha256=sha(__file__)))
    if a.freeze_only:print(protocol);return
    assert os.environ.get('SLURM_JOB_ID'),'Whole-experiment derived analysis runs on SLURM'
    start=time.monotonic();status=json.loads((RUN/'results/partial/status.json').read_text())
    assert status['complete'] and status['scientific_cases_available']==936
    cachepath=RUN/'results/partial/case_cache.json';cache=json.loads(cachepath.read_text());assert len(cache)==936
    rows=[];coverage=[];identities=[]
    for ordinal,(case_id,case) in enumerate(cache.items()):
        d=case['descriptors'][0];base={k:d[k] for k in ('model','platform','size','fraction','kind','seed')};base['case_id']=case_id
        path=Path(case['profile']);assert sha(path)==case['profile_sha256']
        step=PROTOCOL['physical_distance_bin_width'][d['platform']]
        maximum=30. if d['platform']=='chemical' else 1.6
        edges=np.arange(0,round(maximum/step)+1)*step
        rr,cc=assess(path,d['platform'],edges)
        rows.extend(dict(**base,**r) for r in rr);coverage.append(dict(**base,**cc))
        identities.append(dict(**base,path=str(path),sha256=sha(path)))
        if ordinal%96==0:print(json.dumps(dict(cases=ordinal+1,elapsed_seconds=time.monotonic()-start)),flush=True)
    cases=pd.DataFrame(rows);cases.to_csv(data/'case_distance_probabilities.csv',index=False)
    pd.DataFrame(coverage).to_csv(data/'case_exposure_coverage.csv',index=False)
    pd.DataFrame(identities).to_csv(data/'source_profiles.csv',index=False)
    metrics=['probability_lower','probability_upper','endpoint_qualified_coverage','endpoint_overlap_coverage','mean_available_exterior_arc']
    grouped=[]
    keys=['model','size','fraction','kind','distance_lo','distance_hi']
    for key,g in cases.groupby(keys,sort=True):
        row=dict(zip(keys,key),cases=len(g),exposure_frames=int(g.bin_exposure.sum()),
            contributing_cases=int((g.bin_exposure>0).sum()),event_frames=int(g.event_frames.sum()),
            unresolved_frames=int(g.unresolved_secondary_frames.sum()))
        for metric in metrics:
            v=g[metric].dropna();row[metric]=v.mean() if len(v) else np.nan
            row[metric+'_seed_min']=v.min() if len(v) else np.nan;row[metric+'_seed_max']=v.max() if len(v) else np.nan
        grouped.append(row)
    pd.DataFrame(grouped).to_csv(data/'condition_distance_probabilities.csv',index=False)
    manifest=json.loads((PILOT/'manifest.json').read_text());dense=[];densecoverage=[];denseidentities=[]
    for job in manifest['jobs']:
        if job.get('variant')!='standard' or not job.get('dense_observation'):continue
        path=PILOT/'workers'/f'{job["slot"]:04d}'/'dense_primary.npz';assert path.exists()
        step=PROTOCOL['physical_distance_bin_width'][job['platform']]
        edges=np.arange(0,round((30. if job['platform']=='chemical' else 1.6)/step)+1)*step
        base={k:job[k] for k in ('model','size','fraction','seed')}
        for phase in range(10):
            rr,cc=assess(path,job['platform'],edges,phase=phase,stride=10)
            dense.extend(dict(**base,phase=phase,**r) for r in rr)
            densecoverage.append(dict(**base,**cc))
        denseidentities.append(dict(**base,path=str(path),sha256=sha(path)))
    assert len(denseidentities)==4
    pd.DataFrame(dense).to_csv(data/'dense_phase_distance_probabilities.csv',index=False)
    pd.DataFrame(densecoverage).to_csv(data/'dense_phase_exposure_coverage.csv',index=False)
    pd.DataFrame(denseidentities).to_csv(data/'dense_profile_sources.csv',index=False)
    receipt=dict(generated_at=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,
        native_integrations=0,cases=len(cache),source_sha256=sha(__file__),protocol_sha256=sha(protocol),
        cache_sha256=sha(cachepath),native_generation=status['generated_at'],independent_review='pending',
        outputs={str(p.relative_to(out)):sha(p) for p in sorted(data.iterdir())})
    write(metadata/'ANALYSIS_RECEIPT.json',receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
