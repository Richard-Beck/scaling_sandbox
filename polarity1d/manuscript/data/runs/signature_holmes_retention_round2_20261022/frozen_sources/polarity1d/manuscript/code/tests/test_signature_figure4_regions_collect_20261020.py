from pathlib import Path
import sys
import csv
import gzip
import json

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure4_regions_collect_20261020 import (
    aggregate, average, collect, estimands, rmst_bounds, sha, survival_arrays, weighted_survival)
from signature_figure4_regions_20261020 import track_profiles


def movie():
    q=(np.arange(100)+.5)/100
    old=np.exp(-.5*(((q-.2+.5) % 1-.5)/.06)**2)
    new=np.exp(-.5*(((q-.7+.5) % 1-.5)/.06)**2)
    profiles=np.array([.1+old,.1+old+new,.1+.1*old+new,.1+.1*old+new])
    return track_profiles(np.arange(4.)*80,profiles,20.,basal_contrast=1.)


def test_bootstrap_unit_is_seed_and_missing_seed_estimands_are_reportable():
    values=np.array([1.,np.nan,3.])
    weights=np.array([[1.,1.,1.],[0.,3.,0.],[0.,0.,3.]])
    result=average(values,weights)
    assert result[0]==2 and np.isnan(result[1]) and result[2]==3


def test_positive_km_tail_after_last_censor_remains_unidentified():
    result=dict(tracks=[dict(left_censored=False,right_censored=True,duration=2.)])
    grid,risk,event,censor=survival_arrays([result],clock="observed",horizon=10.)
    curve=weighted_survival(risk,event,np.ones((1,1)))
    assert np.isnan(curve[0,-1])
    lower,upper=rmst_bounds(grid,curve,10.)
    assert lower[0]==2 and upper[0]==10


def test_extinct_km_tail_is_identified_to_fixed_horizon():
    result=dict(tracks=[dict(left_censored=False,right_censored=False,duration=2.)])
    grid,risk,event,censor=survival_arrays([result],clock="observed",horizon=10.)
    curve=weighted_survival(risk,event,np.ones((1,1)))
    lower,upper=rmst_bounds(grid,curve,10.)
    assert lower[0]==upper[0]==2 and curve[0,-1]==0


def test_missing_previous_dominant_is_unscorable_even_with_no_recorded_takeovers():
    result=movie()
    result["tracks"][1]["previous_dominant_id"]=None
    assert estimands(result)["ever_unscorable_takeover_fraction"]==1


def test_condition_distributions_and_survival_keep_censoring_and_whole_seed_counts():
    result=movie()
    metadata=dict(model="wave_pinning",C_label=20.)
    metrics,distributions,survival=aggregate([result,result],metadata=metadata,draws=20)
    count=next(x for x in metrics if x["statistic"]=="region_count_mean")
    assert count["independent_seeds"]==2 and count["limited_two_seed_ensemble"]
    assert count["ci_lower"]==pytest.approx(count["ci_upper"])
    assert sum(x["probability"] for x in distributions if x["statistic"]=="simultaneous_region_count")==pytest.approx(1.)
    assert all(x["at_risk"]<=2 for x in survival)
    assert {x["lifetime_clock"] for x in survival}=={"observed","lower","upper"}


def test_deterministic_controls_are_not_bootstrapped_as_independent_seeds():
    result=movie()
    metadata=dict(model="wave_pinning",C_label=20.,kind="mean")
    metrics,_,_=aggregate([result],metadata=metadata,draws=20)
    assert all(row["independent_seeds"]==0 and row["deterministic_controls"]==1
               and row["bootstrap_resamples"]==0 for row in metrics)
    with pytest.raises(ValueError,match="deterministic"):
        aggregate([result,result],metadata=metadata,draws=20)


def test_long_survival_keeps_full_observation_range_and_fixed_comparison_rmst():
    result=movie()
    for track in result["tracks"]:
        if not track["left_censored"]:
            track.update(right_censored=True,duration=20000.)
    metadata=dict(model="wave_pinning",C_label=20.,kind="noise",native_horizon=128000.)
    metrics,_,survival=aggregate([result],metadata=metadata,draws=20)
    assert max(row["time"] for row in survival)==64000.
    rmst=[row for row in metrics if row["statistic"].startswith("restricted_mean_lifetime_")]
    assert all(row["integration_horizon"]==16000. and row["estimate"]==16000. for row in rmst)


def test_full_collector_preserves_condition_kind_across_every_variant_and_table(tmp_path,monkeypatch):
    # This boundary regression covers the sink-loop shadowing error: every
    # variant after the first used to receive kind='survival'.
    import signature_figure4_regions_collect_20261020 as module
    cases=[dict(id=kind,model="wave_pinning",C_label=20.,f=.01,kind=kind,sigma=sigma)
           for kind,sigma in (("noise",1.),("mean",0.))]
    inventory=tmp_path/"inventory.json"
    inventory.write_text(json.dumps(dict(cases=cases)))
    workers=tmp_path/"workers"
    for case in cases:
        root=workers/case["id"];root.mkdir(parents=True)
        variants=[]
        for i,(threshold,factor,smoothed) in enumerate(((.3,1,False),(.5,2,True))):
            path=root/f"v{i}.json.gz"
            with gzip.open(path,"wt") as stream:json.dump(dict(result={}),stream)
            variants.append(dict(threshold=threshold,downsample_factor=factor,smoothed=smoothed,
                                 path=str(path),sha256=sha(path)))
        (root/"result.json").write_text(json.dumps(dict(passed=True,case=case,
            inventory_sha256=sha(inventory),source_sha256="adapter",tracking_source_sha256="tracking",
            sensitivities=variants,primary=variants[1],raw_profile_cadence=40.,common_imaging_cadence=80.)))
    def fake_aggregate(results,*,metadata,seed):
        return ([dict(metadata,statistic="region_count_mean",estimate=1.)],
                [dict(metadata,statistic="count",probability=1.)],
                [dict(metadata,time=0.,survival=1.)])
    monkeypatch.setattr(module,"aggregate",fake_aggregate)
    output=tmp_path/"export"
    collect(inventory,workers,output,expected_inventory_sha256=sha(inventory),worker_sources={
        "signature_figure4_regions_campaign_20261020.py":"adapter",
        "signature_figure4_regions_20261020.py":"tracking"})
    for table in ("metrics","distributions","survival"):
        rows=list(csv.DictReader((output/(table+".csv")).open()))
        assert len(rows)==4
        assert {r["kind"] for r in rows}=={"noise","mean"}
        primary=[r for r in rows if r["primary"]=="True"]
        assert len(primary)==2 and {r["kind"] for r in primary}=={"noise","mean"}


def test_same_seed_at_distinct_horizons_is_separate_but_same_horizon_companion_is_rejected(tmp_path,monkeypatch):
    import signature_figure4_regions_collect_20261020 as module
    cases=[dict(id=str(i),model="wave_pinning",C_label=20.,f=.01,kind="noise",sigma=1.,
                seed=123,native_horizon=horizon) for i,horizon in enumerate((32000.,128000.))]
    inventory=tmp_path/"inventory.json"
    workers=tmp_path/"workers"
    def inputs():
        inventory.write_text(json.dumps(dict(cases=cases,analysis_arm="longest")))
        for case in cases:
            root=workers/case["id"];root.mkdir(parents=True,exist_ok=True)
            path=root/"v.json.gz"
            with gzip.open(path,"wt") as stream:json.dump(dict(result={}),stream)
            variant=dict(threshold=.5,downsample_factor=2,smoothed=True,path=str(path),sha256=sha(path))
            (root/"result.json").write_text(json.dumps(dict(passed=True,case=case,
                inventory_sha256=sha(inventory),source_sha256="adapter",tracking_source_sha256="tracking",
                sensitivities=[variant],primary=variant,raw_profile_cadence=40.,common_imaging_cadence=80.)))
    def fake_aggregate(results,*,metadata,seed):
        return ([dict(metadata,statistic="region_count_mean",estimate=1.)],[],[])
    monkeypatch.setattr(module,"aggregate",fake_aggregate)
    sources={"signature_figure4_regions_campaign_20261020.py":"adapter",
             "signature_figure4_regions_20261020.py":"tracking"}
    inputs()
    collect(inventory,workers,tmp_path/"good",expected_inventory_sha256=sha(inventory),worker_sources=sources)
    rows=list(csv.DictReader((tmp_path/"good/metrics.csv").open()))
    assert {float(r["native_horizon"]) for r in rows}=={32000.,128000.}
    assert {r["experiment"] for r in rows}=={"longest"}
    cases[1]["native_horizon"]=32000.;inputs()
    with pytest.raises(ValueError,match="companions"):
        collect(inventory,workers,tmp_path/"bad",expected_inventory_sha256=sha(inventory),worker_sources=sources)
