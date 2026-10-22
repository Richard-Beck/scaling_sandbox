"""Independent native-to-CSV and decision audit for the completed Figure 4.

Recomputes heading cosines via phase differences, rather than calling the
production extraction or rate helpers. Original native inputs are read only.
"""
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from pypdf import PdfReader

csv.field_size_limit(50_000_000)

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures/length_mechanisms_20261008/dev/revision_round2_20261022/figure4"
DATA = FIG / "supporting/data"
META = FIG / "supporting/metadata"
VIEWS = {"main":(.5,1,1),"analysis_first_half":(.5,.75,1),
    "analysis_second_half":(.75,1,1),"late_quarter":(.875,1,1),
    "cadence_x2":(.5,1,2),"cadence_x4":(.5,1,4),"cadence_x8":(.5,1,8)}


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(4 << 20),b""):
            value.update(block)
    return value.hexdigest()


def number(value):
    return float(value) if value not in (None,"") else math.nan


def equal(actual, expected):
    np.testing.assert_allclose(actual,expected,rtol=2e-10,atol=2e-10,equal_nan=True)


def key(row):
    return row["arm"],row["model"],float(row["C_label"]),float(row["f"]),row["view"]


def main():
    receipt=read(META/"analysis_receipt.json")
    assert receipt["existing_outputs_only"] and receipt["no_simulations"]
    for item in [receipt["input_condition_metrics"],*receipt["manifests"],*receipt["scalar_sources"]]:
        assert sha(item["path"])==item["sha256"],item["path"]
    for name,digest in receipt["outputs"].items():
        assert sha(DATA/name)==digest,name
    condition=rows(DATA/"condition_lag_metrics.csv")
    seeds=rows(DATA/"seed_lag_summaries.csv")
    summaries={(key(r),r["case_id"],float(r["lag"])):r for r in seeds}
    assert len(summaries)==len(seeds)
    metrics=rows(DATA/"old_tau_condition_metrics.csv")
    manifests={}
    for item in receipt["manifests"]:
        for group in read(item["path"])["groups"]:
            manifests[group["stage"],group["model"],float(group["C_label"]),float(group["f"])]=group
    used={}
    for r in metrics:
        origin=r.get("common_observation_arm") or r["selected_observation_arm"]
        stage="E4" if r["selected_observation_arm"]=="preselected_E4" else "E3" if origin=="common_new_E3_size" else "E2"
        groupkey=stage,r["model"],float(r["C_label"]),float(r["f"])
        used[groupkey]=manifests[groupkey]
    # One raw extraction supports every common/longest row selecting that cohort.
    raw={};raw_circumference={};raw_count=0
    for groupkey,group in used.items():
        lagset=np.array([250,500,1000,2000,4000,8000],float)/(10 if group["model"]=="spring_ring" else 1)
        for case in group["noise"]:
            with np.load(case["scalar_path"],allow_pickle=False) as z:
                time=z["time"];heading=z["P_real"]+1j*z["P_imag"]
                circumference=np.broadcast_to(z["circumference"],time.shape)
                for view,(start,stop,stride) in VIEWS.items():
                    indices=np.flatnonzero((time>=case["native_horizon"]*start-1e-8)&(time<=case["native_horizon"]*stop+1e-8))[::stride]
                    p=heading[indices];dt=float(time[indices[1]]-time[indices[0]])
                    raw_circumference[groupkey,view,case["id"]]=float(circumference[indices].mean())
                    valid=np.abs(p)>=.1
                    for lag in lagset:
                        step=round(lag/dt);value=(math.nan,math.nan,math.nan)
                        if step>=1 and step<len(p) and math.isclose(step*dt,lag,abs_tol=1e-7):
                            mask=valid[:-step]&valid[step:]
                            phases=np.angle(p[step:][mask]*np.conj(p[:-step][mask]))
                            value=(float(np.cos(phases).sum()),int(mask.sum()),len(p)-step)
                        raw[groupkey,view,case["id"],float(lag)]=value;raw_count+=1
    selected=defaultdict(list)
    for r in metrics:
        origin=r.get("common_observation_arm") or r["selected_observation_arm"]
        stage="E4" if r["selected_observation_arm"]=="preselected_E4" else "E3" if origin=="common_new_E3_size" else "E2"
        groupkey=stage,r["model"],float(r["C_label"]),float(r["f"])
        group=used[groupkey]
        for case in group["noise"]:
            for lag in np.array([250,500,1000,2000,4000,8000],float)/(10 if r["model"]=="spring_ring" else 1):
                source=summaries[key(r),case["id"],float(lag)]
                values=raw[groupkey,r["view"],case["id"],float(lag)]
                equal([number(source[k]) for k in ("cosine_sum","eligible_pairs","possible_pairs")],values)
                selected[key(r),float(lag)].append(values)
    boot=np.load(DATA/"whole_seed_bootstrap.npz",allow_pickle=False)
    bootstrap_conditions=0
    for r in metrics:
        origin=r.get("common_observation_arm") or r["selected_observation_arm"]
        stage="E4" if r["selected_observation_arm"]=="preselected_E4" else "E3" if origin=="common_new_E3_size" else "E2"
        groupkey=stage,r["model"],float(r["C_label"]),float(r["f"])
        group=used[groupkey];cases=group["noise"]
        lags=np.array([250,500,1000,2000,4000,8000],float)/(10 if r["model"]=="spring_ring" else 1)
        values=np.array([[raw[groupkey,r["view"],case["id"],float(lag)] for lag in lags] for case in cases])
        pairing=group.get("bootstrap_pairing_key",f"{group['model']}:{group['C_label']}:{group['f']}")
        seed=2026102004+int(hashlib.sha256(pairing.encode()).hexdigest()[:8],16)
        weights=np.random.default_rng(seed).multinomial(len(cases),np.full(len(cases),1/len(cases)),size=2000)
        with np.errstate(divide="ignore",invalid="ignore"):
            expected=(weights@values[:,:,0])/(weights@values[:,:,1])
        bkey="|".join([r["arm"],r["model"],str(float(r["C_label"])),str(float(r["f"])),r["view"]])
        equal(boot[bkey+"|C"],expected)
        mean_circumference=np.array([raw_circumference[groupkey,r["view"],case["id"]] for case in cases])
        equal(boot[bkey+"|circumference"],weights@mean_circumference/len(cases))
        bootstrap_conditions+=1
    for r in condition:
        values=np.asarray(selected[key(r),float(r["lag"])])
        total=values.sum(axis=0)
        with np.errstate(invalid="ignore",divide="ignore"):
            corr=total[0]/total[1];coverage=total[1]/total[2]
        equal([number(r[k]) for k in ("C","coverage","eligible_pairs","possible_pairs")],
              [corr,coverage,total[1],total[2]])
        rate=-math.log(min(corr,1))/float(r["lag"]) if corr>0 else math.nan
        equal(number(r["k"]),rate)
        bkey="|".join([r["arm"],r["model"],str(float(r["C_label"])),str(float(r["f"])),r["view"]])+"|C"
        lagindex=list(np.array([250,500,1000,2000,4000,8000],float)/(10 if r["model"]=="spring_ring" else 1)).index(float(r["lag"]))
        draws=boot[bkey][:,lagindex];valid=np.isfinite(draws)
        if valid.mean()>=.975 or not math.isfinite(corr):
            limits=np.quantile(draws[valid],[.025,.975],method="inverted_cdf") if valid.any() else [math.nan,math.nan]
        else:
            limits=[-1,1]
        equal([number(r["C_ci_lower"]),number(r["C_ci_upper"])],limits)
        if math.isfinite(corr):
            low=-math.log(min(limits[1],1))/float(r["lag"]) if limits[1]>0 else math.nan
            high=-math.log(min(limits[0],1))/float(r["lag"]) if limits[0]>0 else math.inf
        else:
            low=high=math.nan
        equal([number(r["k_ci_lower"]),number(r["k_ci_upper"])],[low,high])
    boot.close()
    ratios=rows(DATA/"size_ratios.csv");tolerances=rows(DATA/"noise_tolerance.csv")
    all_lags=rows(DATA/"all_lag_size_ratios.csv")
    for r in ratios:
        if r["status"]=="unavailable_unmatched_designed_input":continue
        rr=[x for x in all_lags if (x["arm"],x["model"],x["view"],x["f"])==(r["arm"],r["model"],r["view"],r["f"])]
        eligible=[x for x in rr if x["usable"]=="True" and 0<number(x["C_small"])<=.98 and number(x["C_large"])>0]
        if not eligible:
            assert r["evaluable"]=="False";continue
        assert float(r["primary_lag"])==min(float(x["lag"]) for x in eligible)
        if math.isfinite(number(r["rho_k"])):
            equal(number(r["rho_k"]),number(r["k_small"])/number(r["k_large"]))
            equal(number(r["rho_geo"]),number(r["rho_k"])/number(r["geometry_expectation"]))
        drift=[number(x["rho_k"]) for x in rr if x["usable"]=="True" and math.isfinite(number(x["rho_k"])) and number(x["rho_k"])>0]
        assert (r["non_exponential"]=="True")==bool(len(drift)>1 and max(drift)/min(drift)>1.5)
    calls=rows(DATA/"S0_posthoc_H1.csv")
    assert len(calls)==140
    callkey=lambda r:(r["arm"],r["view"],r["forced_inclusive"],r["model"])
    for call in calls:
        if call["model"]=="all":
            continue
        ts=sorted([r for r in tolerances if callkey(r)==callkey(call)],key=lambda r:float(r["C_label"]))
        cutoff=number(ts[-1]["fhalf"])
        rr=[r for r in ratios if (r["arm"],r["view"],r["model"])==(call["arm"],call["view"],call["model"]) and float(r["f"])<cutoff]
        assert json.loads(call["retained_inputs"])==[float(r["f"]) for r in rr]
        if not math.isfinite(cutoff):
            assert call["status"]=="inconclusive";continue
        eligible=[r for r in rr if r["evaluable"]=="True" and (call["forced_inclusive"]=="True" or r["forced_combined"]=="False")]
        unknown=call["forced_inclusive"]=="False" and any(r["evaluable"]=="True" and r["forced_combined"] not in ("True","False") for r in rr)
        reverse=any(number(r["rho_k_ci_upper"])<1 for r in eligible)
        gain=bool(eligible) and not unknown and all(number(r["rho_k_ci_lower"])>1 for r in eligible)
        tolpass=number(ts[-1]["fhalf_ci_lower"])>number(ts[0]["fhalf_ci_upper"])
        tolfail=number(ts[-1]["fhalf_ci_upper"])<=number(ts[0]["fhalf_ci_lower"])
        assert call["status"]==("fail" if reverse or tolfail else "pass" if gain and tolpass else "inconclusive"),callkey(call)
    for call in calls:
        if call["model"]!="all":continue
        component=[r for r in calls if callkey(r)[:3]==callkey(call)[:3] and r["model"]!="all"]
        yes=sum(r["status"]=="pass" for r in component);unknown=sum(r["status"]=="inconclusive" for r in component)
        assert call["status"]==("pass" if yes>=3 else "fail" if yes+unknown<3 else "inconclusive")
    retired=rows(DATA/"S0_retired_gates.csv")+rows(DATA/"S0_retired_A6.csv")
    assert len(retired)==196 and all(r["status"]=="retired" for r in retired)
    lag=read(META/"all_lag_render_receipt.json")
    assert len(lag["pages"])==56 and lag["input_sha256"]==sha(DATA/"all_lag_size_ratios.csv")
    assert len(PdfReader(FIG/"figure4_all_lag_ratios.pdf").pages)==56
    render=read(META/"main_render_receipt.json")
    for name,digest in render["outputs"].items():assert sha(FIG/name)==digest,name
    result=dict(passed=True,generated_utc=datetime.now(timezone.utc).isoformat(),
        scope="Every accepted native scalar hash and every seed/lag cosine/count; every whole-seed correlation/circumference bootstrap draw; every condition correlation/rate/CI; independent amended H1 arithmetic; all-lag presentation coverage.",
        native_scalar_hashes=len(receipt["scalar_sources"]),native_cohorts=len(used),raw_seed_view_lags=raw_count,
        exported_seed_lags=len(seeds),condition_lag_rows=len(condition),bootstrap_conditions=bootstrap_conditions,bootstrap_resamples=2000,H1_calls=len(calls),retired_gate_rows=len(retired),
        all_lag_pages=56,main_pages=len(PdfReader(FIG/"figure4_main.pdf").pages),
        supplement_pages=len(PdfReader(FIG/"figure4_supplement.pdf").pages))
    (META/"completion_scientific_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))


if __name__=="__main__":
    main()
