"""Frozen gap-noise campaign for linear Wave-pinning and corrected deBelly.

Preparation never submits. Native integrations run only through explicit SLURM
workers. A six-case bounded-input pilot precedes the complete 432 stochastic,
36 matched-mean and 18 zero-input cases. The original 12-case lognormal pilot
is a separate intermittent-input diagnostic. Recovery is always noise free.
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import platform
import resource
import shlex
import shutil
import sys
import time
import traceback

import signature_campaign_20261020 as campaign

ROOT=campaign.ROOT
PACKAGE=campaign.PACKAGE
MODELS=("wave_pinning","debelly")
LENGTHS=(19,29,39)
LEVELS=(.005,.02)
GAPS=(1.,4.,16.)


def condition_grid(mode="full"):
    if mode not in ("pilot","bounded_pilot","full"):
        raise ValueError("Expected pilot or full noise mode")
    rows=[]
    for mi,model in enumerate(MODELS):
        for li,length in enumerate(LENGTHS):
            if mode=="pilot" and length==29:continue
            if mode=="bounded_pilot" and length!=39:continue
            for gi,gap in enumerate(GAPS):
                if mode in ("pilot","bounded_pilot") and gap!=16.:continue
                for fi,level in enumerate(LEVELS):
                    if mode in ("pilot","bounded_pilot") and level!=.02:continue
                    for replicate in range(1 if mode in ("pilot","bounded_pilot") else 12):
                        seed=202610203000000+mi*1000000+li*100000+fi*10000+gi*100+replicate
                        rows.append(dict(model=model,requested_reference_length=length,
                            kind="noise",mean_fraction=level,replicate_index=replicate,seed=seed,
                            gap_multiple_of_pulse=gap))
                    rows.append(dict(model=model,requested_reference_length=length,
                        kind="mean",mean_fraction=level,replicate_index=0,seed=0,
                        gap_multiple_of_pulse=gap))
                rows.append(dict(model=model,requested_reference_length=length,
                    kind="zero",mean_fraction=0.,replicate_index=0,seed=0,
                    gap_multiple_of_pulse=gap))
    expected=12 if mode=="pilot" else 6 if mode=="bounded_pilot" else 486
    if len(rows)!=expected:
        raise ValueError("Noise condition grid has incorrect coverage")
    for index,row in enumerate(rows):row["task"]=index
    return rows


def protocol(transform="bounded_tanh"):
    if transform not in ("bounded_tanh","lognormal"):
        raise ValueError("Unknown gap-input transform")
    return dict(models=list(MODELS),reference_lengths=list(LENGTHS),cells=201,
        pulse_count=4,total_on_factor=1.5,gap_multiples=list(GAPS),
        mean_fractions=list(LEVELS),independent_seeds=12,
        ell=1.9,tau=dict(wave_pinning=10.,debelly=5.),
        knot=dict(wave_pinning=2.5,debelly=.5),sigma=1.,
        maximum_step=dict(wave_pinning=1.25,debelly=.5),
        input_transform=transform,
        input_formula="mean*(1+tanh(sigma*Z))" if transform=="bounded_tanh" else "mean*exp(sigma*Z-sigma^2/2)",
        hard_maximum_fraction_of_native_peak=.04 if transform=="bounded_tanh" else None,
        forcing="Stationary nonperiodic spatial-exponential/temporal-OU Gaussian knots; variance-normalized interpolation; explicitly declared nonlinear activation-input map.",
        covariance_scope="Gaussian ell and tau describe latent noise. Transformed input correlations are measured separately and are not assumed exponential.",
        application="Inter-pulse gaps only; latent input uses whole-train absolute time. Cue ON and final recovery have no gap input.",
        recovery=campaign.RECOVERY,primary_retention=dict(
            denominator="Pulse endpoint eligible and D_pulse>=0",
            success="Gap eligible, D_gap>=0 and share_gap/share_pulse>=threshold",
            threshold=.8,sensitivity_thresholds=[.7,.8,.9],
            observed_ineligible_gap="Known loss, retained in eligible-pulse denominator",
            missing_native_endpoint="Error; never a retention loss or finite censoring"),
        supplemental_retention="Basal-gain-normalized share ratio from the independent analysis helper",
        field_tail="Full native fields in noise-free final recovery; tail diagnostic does not define replacement or gap retention.",
        physical_calibration=False)


class AppliedInputAudit:
    """Track input actually queried during native gap solves, without clipping."""
    event_times=()

    def __init__(self,forcing,native_peak):
        self.forcing,self.native_peak=forcing,float(native_peak)
        self.calls,self.values,self.above_reference=0,0,0
        self.minimum,self.maximum,self.sum= float("inf"),float("-inf"),0.
        self.maximum_time=None

    def __call__(self,time,x,length):
        import numpy as np
        value=self.forcing(time,x,length)
        minimum,maximum=float(value.min()),float(value.max())
        self.minimum=min(self.minimum,minimum)
        if maximum>self.maximum:self.maximum,self.maximum_time=maximum,float(time)
        self.calls+=1;self.values+=value.size;self.sum+=float(np.sum(value))
        self.above_reference+=int(np.count_nonzero(value>self.native_peak))
        return value

    def audit(self):
        return dict(native_rhs_calls=self.calls,queried_spatial_values=self.values,
            minimum=self.minimum if self.values else None,maximum=self.maximum if self.values else None,
            maximum_global_time=self.maximum_time,
            maximum_over_native_reference_peak=self.maximum/self.native_peak if self.values else None,
            fraction_queried_values_above_native_reference_peak=self.above_reference/self.values if self.values else None,
            evaluation_weighted_mean=self.sum/self.values if self.values else None,
            scope="Actual solver queries during gaps, including adaptive/repeated evaluations; this mean is not a time average.",
            clipping_or_realized_normalization=False)


def sampled_gap_distribution(forcing,schedule,native_peak):
    """Audit the applied-gap marginal on half-knot time samples and native cells."""
    import numpy as np
    gaps=[(start,stop) for kind,start,stop in schedule.segments if kind=="gap"]
    if not gaps:raise ValueError("Noise protocol requires positive inter-pulse gaps")
    total_duration=sum(stop-start for start,stop in gaps)
    x=(np.arange(forcing.cells)+.5)*forcing.length/forcing.cells
    sample_interval=forcing.knot/2.
    count=sum(len(np.arange(start,stop,sample_interval))+1 for start,stop in gaps)
    if forcing.kind!="noise":
        value=forcing.mean if forcing.kind=="mean" else 0.
        return dict(sample_interval=sample_interval,temporal_samples=count,spatial_samples=forcing.cells,
            gap_duration=total_duration,time_averaged_spatial_mean=value,minimum=value,maximum=value,
            quantiles={key:value for key in ("0.5","0.9","0.99","0.999")},
            maximum_over_native_reference_peak=value/native_peak,
            fraction_samples_above_native_reference_peak=float(value>native_peak),
            scope="Exact deterministic applied input; no random generator consumed.")
    parts=[];integrated_mean=0.
    for start,stop in gaps:
        times=np.r_[np.arange(start,stop,sample_interval),stop]
        values=np.asarray([forcing(float(t),x,forcing.length) for t in times])
        integrated_mean+=float(np.trapezoid(values.mean(axis=1),times))
        parts.append(values.ravel())
    values=np.concatenate(parts)
    quantiles=np.quantile(values,[.5,.9,.99,.999])
    return dict(sample_interval=sample_interval,temporal_samples=count,spatial_samples=forcing.cells,
        gap_duration=total_duration,time_averaged_spatial_mean=integrated_mean/total_duration,
        minimum=float(values.min()),maximum=float(values.max()),
        quantiles={str(q):float(value) for q,value in zip((.5,.9,.99,.999),quantiles)},
        maximum_over_native_reference_peak=float(values.max()/native_peak),
        fraction_samples_above_native_reference_peak=float(np.mean(values>native_peak)),
        scope="Applied gaps only, half-knot temporal samples and native spatial cells; quantiles are sample-weighted and maximum is sampled, not a continuous-time supremum.",
        clipping_or_realized_normalization=False)


def prepare(out,mode,calibrations_path,pilot_audit=None,transform="bounded_tanh"):
    out=Path(out).resolve();calibrations_path=Path(calibrations_path).resolve()
    if out.exists() or not out.is_relative_to((PACKAGE/"data/runs").resolve()):
        raise ValueError("Choose a fresh campaign directory under manuscript/data/runs")
    collection=campaign.verified_calibration_collection(calibrations_path)
    calibrations={row["calibration_id"]:row for row in collection["calibrations"]}
    source_root=calibrations_path.parent.parent
    basals={row["basal_id"]:row for row in campaign.read(source_root/"results/basals.json")["basals"]}
    if mode=="full":
        if pilot_audit is None or campaign.read(Path(pilot_audit)).get("native_execution_ready") is not True:
            raise ValueError("Full noise preparation requires a completed, valid native pilot audit")
        if campaign.read(Path(pilot_audit)).get("input_transform","lognormal")!=transform:
            raise ValueError("Full noise forcing must match its independently completed native pilot")
    config=protocol(transform);cases=[];checkpoints={}
    for row in condition_grid(mode):
        model=row["model"];length=row["requested_reference_length"];basal_id=f"{model}:L{length}"
        calibration=calibrations.get(basal_id,{})
        basal=basals.get(basal_id,{})
        if not calibration.get("calibration",{}).get("train_ready") or not basal.get("qualified"):
            raise ValueError("Missing qualified post-withdrawal calibration/basal: "+basal_id)
        source=Path(basal["checkpoint"])
        digest=calibration["basal_sha256"]
        if not source.is_file() or campaign.sha(source)!=digest or digest!=basal["checkpoint_sha256"]:
            raise ValueError("Noise and calibration must use the exact same basal")
        relative=f"inputs/basals/{model}_L{length}.npz"
        checkpoints[relative]=dict(source=source,sha256=digest)
        peak=campaign.NATIVE_PEAK[model]
        cases.append(dict(row,stage="gap_noise",mode=mode,basal_id=basal_id,
            native_length=float(length),cells=201,checkpoint=relative,checkpoint_sha256=digest,
            fixed_scales=basal["fixed_scales"],baseline_native_peak=peak,
            pulse_count=4,total_on_factor_of_calibration=1.5,
            calibration_duration=calibration["calibration"]["upper_committing_duration"],
            cue_width=3.8,cue_area=peak*3.8,recovery_duration=campaign.RECOVERY[model],
            noise_mean=peak*row["mean_fraction"],noise_ell=1.9,noise_tau=config["tau"][model],
            noise_knot=config["knot"][model],noise_sigma=1.,
            noise_transform=transform,
            max_step=1.25 if model=="wave_pinning" else None))
    for directory in ("inputs/basals","frozen_sources","workers","logs","results"):
        (out/directory).mkdir(parents=True,exist_ok=True)
    inputs={};sources={}
    for relative,record in checkpoints.items():
        shutil.copy2(record["source"],out/relative)
        if campaign.sha(out/relative)!=record["sha256"]:raise ValueError("Basal changed while freezing")
        inputs[relative]=record["sha256"]
    runtime_data=(PACKAGE/"data/inputs/fixed_scales.json",
        PACKAGE/"data/runs/debelly_corrected_three_test_20261008/manifest.json")
    for source in [path for path in sorted(campaign.CODE.rglob("*"))
            if path.is_file() and "__pycache__" not in path.parts and path.suffix in (".py",".json",".txt",".md")]+list(runtime_data):
        relative=str(source.relative_to(ROOT));target=out/"frozen_sources"/relative
        target.parent.mkdir(parents=True,exist_ok=True);digest=campaign.sha(source)
        shutil.copy2(source,target)
        if campaign.sha(source)!=digest or campaign.sha(target)!=digest:
            raise ValueError("Runtime source changed during freeze: "+relative)
        sources[relative]=digest
    for source in (calibrations_path,ROOT/"docs/tasks/briefs/20261020_figure2_decisions.json"):
        target=out/"inputs"/source.name;digest=campaign.sha(source);shutil.copy2(source,target)
        if campaign.sha(target)!=digest or campaign.sha(source)!=digest:raise ValueError("Input collection changed")
        inputs[str(target.relative_to(out))]=digest
    if pilot_audit:
        source=Path(pilot_audit).resolve();target=out/"inputs/pilot_audit.json"
        shutil.copy2(source,target);inputs[str(target.relative_to(out))]=campaign.sha(target)
    campaign.write(out/"config.json",config);inputs["config.json"]=campaign.sha(out/"config.json")
    manifest=dict(created_utc=campaign.now(),stage="gap_noise",mode=mode,status="frozen_not_submitted",
        result_schema=2,cases=cases,workers=len(cases),logical_protocols=len(cases),
        source_sha256=sources,input_sha256=inputs,
        full_brief_noise_counts=dict(noise=432,mean=36,zero=18),
        original_repository_root=str(ROOT),root_review_required=True,automatic_submission=False,
        calibration_source=dict(path=str(calibrations_path),sha256=campaign.sha(calibrations_path)),
        complete_brief_coverage=False,
        scope="Fresh gap-noise assay; pilot numerical execution is not full-ensemble or all-numerics qualification.")
    campaign.write(out/"manifest.json",manifest)
    campaign.write(out/"freeze.json",dict(manifest_sha256=campaign.sha(out/"manifest.json"),generated_utc=campaign.now()))
    worker=out/"worker.sh";frozen=out/"frozen_sources"/Path(__file__).relative_to(ROOT)
    worker.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec "+shlex.quote(campaign.PYTHON)+" -B "+shlex.quote(str(frozen))+" worker --out "+shlex.quote(str(out))+" --task \"$SLURM_ARRAY_TASK_ID\"\n")
    worker.chmod(0o750)
    command=["sbatch","--parsable","--partition=red","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=8G","--time=01:00:00",f"--array=0-{len(cases)-1}","--job-name=F2-gap-noise-"+mode,
        "--output="+str(out/"logs/%A_%a.out"),"--error="+str(out/"logs/%A_%a.err"),str(worker)]
    (out/"submission_command.txt").write_text(shlex.join(command)+"\n")
    campaign.verify(out)
    return dict(output=str(out),mode=mode,workers=len(cases),manifest_sha256=campaign.sha(out/"manifest.json"),submitted=False)


def physical_diagnostics(path,case):
    import numpy as np
    model=case["model"];fields=("active","inactive") if model=="wave_pinning" else ("rac","rho","mca_bound","mca_unbound")
    minimum={};amount={}
    with np.load(path,allow_pickle=False) as archive:
        for name in fields:
            values=archive["field_"+name]
            if not np.isfinite(values).all():raise ValueError("Nonfinite native density field")
            minimum[name]=float(values.min())
            integrated=values.sum(axis=1)*(case["native_length"]/case["cells"])
            amount[name]=dict(initial=float(integrated[0]),final=float(integrated[-1]),
                minimum=float(integrated.min()),maximum=float(integrated.max()))
        conservation=None
        if model=="wave_pinning":
            total=(archive["field_active"]+archive["field_inactive"]).sum(axis=1)*(case["native_length"]/case["cells"])
            drift=float(np.max(abs(total-total[0]))/total[0])
            conservation=dict(initial=float(total[0]),maximum_relative_drift=drift,tolerance=1e-6,passed=bool(drift<=1e-6),pool="active plus inactive")
    positive=all(value>=-1e-7 for value in minimum.values())
    return dict(density_minimum=minimum,positivity_tolerance=1e-7,positivity_pass=positive,
        nominal_activity_amounts=amount,conservation=conservation,
        valid=bool(positive and (conservation is None or conservation["passed"])),
        debelly_scope="Rac/Rho/MCA amount summaries are not claimed conserved; mechanics fields are signed and excluded from density positivity.")


def worker(out,task):
    started=time.monotonic();manifest=campaign.verify(out)
    if not 0<=task<len(manifest["cases"]):raise ValueError("Task is outside frozen noise array")
    case=manifest["cases"][task];directory=out/"workers"/f"{task:04d}"
    if (directory/"result.json").exists():raise FileExistsError("Refusing completed worker overwrite")
    directory.mkdir(parents=True,exist_ok=True);context=campaign.scheduler_context()
    campaign.write(directory/"status.json",dict(generated_utc=campaign.now(),state="running",scheduler=context))
    native,inputs,readouts,_=campaign.modules()
    from signature_noise_20261020 import LinearGapNoise
    from signature_analysis_20261020 import pulse_gap_endpoints
    import numpy as np
    import scipy
    state,basal=native.read_checkpoint(out/case["checkpoint"],case["model"])
    schedule=inputs.make_train(case["calibration_duration"],1.5,4,case["gap_multiple_of_pulse"],
        native_step=.5 if case["model"]=="debelly" else None)
    forcing=LinearGapNoise(case["native_length"],201,schedule.end_time,case["seed"],case["noise_mean"],
        ell=case["noise_ell"],tau=case["noise_tau"],knot=case["noise_knot"],sigma=1.,kind=case["kind"],
        transform=case.get("noise_transform","lognormal"))
    forcing_audit=forcing.audit()
    forcing_audit["sampled_applied_gap_distribution"]=sampled_gap_distribution(forcing,schedule,case["baseline_native_peak"])
    forcing_audit["environment"]=dict(python=sys.version,executable=sys.executable,
        numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),
        threads={key:os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")})
    campaign.write(directory/"forcing_audit.json",forcing_audit)
    audited_input=AppliedInputAudit(forcing,case["baseline_native_peak"])
    trial=campaign._run_trial(native,inputs,readouts,case,state,basal,directory/"trial",
        schedule=schedule,campaign_root=out,gap_input=audited_input,max_step=case["max_step"])
    forcing_audit["actual_native_gap_input_evaluations"]=audited_input.audit()
    if forcing.transform=="bounded_tanh":
        upper=2*case["noise_mean"]
        observed=max(forcing_audit["sampled_applied_gap_distribution"]["maximum"],audited_input.maximum)
        forcing_audit["hard_bound_check"]=dict(passed=bool(observed<=upper+1e-14),
            analytic_upper_bound=upper,observed_sampled_or_queried_maximum=observed,
            upper_bound_fraction_of_native_peak=upper/case["baseline_native_peak"])
        if not forcing_audit["hard_bound_check"]["passed"]:
            raise ValueError("Bounded gap input exceeded its declared hard limit")
    campaign.write(directory/"forcing_audit.json",forcing_audit)
    measurement=campaign.read(directory/"trial/readouts.json")
    endpoints=pulse_gap_endpoints(measurement["series"],trial["events"])
    campaign.write(directory/"gap_retention.json",dict(endpoints=endpoints,
        ineligible_sample_fraction=1.-sum(measurement["series"]["readout_eligible"])/len(measurement["series"]["time"]),
        primary="Raw gap/pulse challenger share ratio, conditional eligible challenger-biased pulse; known gap inactivity is loss.",
        supplemental="Basal-gain-normalized ratio retained separately."))
    diagnostics=physical_diagnostics(directory/"trial/trace.npz",case)
    campaign.write(directory/"physical_diagnostics.json",diagnostics)
    if not diagnostics["valid"]:raise ValueError("Native noise trajectory failed density/conservation validity")
    result=dict(generated_utc=campaign.now(),result_schema=2,case=case,scheduler=context,
        manifest_sha256=campaign.sha(out/"manifest.json"),
        basal=dict(qualified=True,basal_id=case["basal_id"],model=case["model"],
            requested_reference_length=case["requested_reference_length"],native_length=case["native_length"],
            checkpoint=str((out/case["checkpoint"]).resolve()),checkpoint_sha256=case["checkpoint_sha256"],
            fixed_scales=case["fixed_scales"]),
        trial=trial,forcing_audit=forcing_audit,physical_diagnostics=diagnostics,
        gap_retention_file=campaign.artifact_record(directory/"gap_retention.json",out),
        wall_seconds=time.monotonic()-started,
        peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        allocated_bytes_before_final_records=campaign.allocated_bytes(directory))
    result["artifact_sha256"]={str(path.relative_to(out)):campaign.sha(path) for path in campaign.scientific_files(directory)}
    campaign.write(directory/"result.json",result)
    campaign.write(directory/"status.json",dict(generated_utc=campaign.now(),state="complete",scheduler=context))
    return dict(task=task,state="complete",kind=case["kind"],wall_seconds=result["wall_seconds"])


def collect(out,job_id=None,verified_shards=None):
    status=campaign.collect(out,job_id=job_id,verified_shards=verified_shards)
    if not status["complete"]:return status
    manifest=campaign.verify(out);kinds={};valid=True;rows=[];bounds_valid=True;trajectories=[]
    for case in manifest["cases"]:
        directory=out/"workers"/f"{case['task']:04d}";result=campaign.read(directory/"result.json")
        valid=valid and result["physical_diagnostics"]["valid"]
        if case.get("noise_transform")=="bounded_tanh":
            bounds_valid=bounds_valid and result["forcing_audit"]["hard_bound_check"]["passed"]
        kinds[case["kind"]]=kinds.get(case["kind"],0)+1
        sampled=result["forcing_audit"]["sampled_applied_gap_distribution"]
        gap_exposure=(sampled["time_averaged_spatial_mean"]*sampled["gap_duration"]*case["native_length"])
        cue_exposure=case["cue_area"]*result["trial"]["schedule"]["actual_total_on"]
        trajectories.append(dict(task=case["task"],model=case["model"],reference_length=case["requested_reference_length"],
            kind=case["kind"],mean_fraction=case["mean_fraction"],seed=case["seed"],gap_multiple=case["gap_multiple_of_pulse"],
            result_path=str(directory/"result.json"),result_sha256=campaign.sha(directory/"result.json"),
            sampled_applied_gap_distribution=sampled,
            actual_native_gap_input_evaluations=result["forcing_audit"]["actual_native_gap_input_evaluations"],
            hard_bound_check=result["forcing_audit"].get("hard_bound_check"),
            gap_integrated_commanded_exposure=gap_exposure,
            theoretical_mean_gap_integrated_exposure=(case["noise_mean"]*sampled["gap_duration"]*case["native_length"] if case["kind"]!="zero" else 0.),
            on_integrated_commanded_exposure=cue_exposure,
            gap_exposure_over_on_exposure=gap_exposure/cue_exposure,
            exposure_scope="Native nominal-coordinate commanded coefficient integrated in space and time; sampled trapezoidal gap integration. This is not integrated biochemical activation or calibrated photon dose.",
            physical_diagnostics=result["physical_diagnostics"],final_d=result["trial"]["final_d"],
            final_field_tail_pass=result["trial"]["field_tail_pass"]))
        retention=campaign.read(directory/"gap_retention.json")
        for gap_index,gap in enumerate(row for row in retention["endpoints"] if row["kind"]=="gap"):
            rows.append(dict(task=case["task"],model=case["model"],reference_length=case["requested_reference_length"],
                kind=case["kind"],mean_fraction=case["mean_fraction"],seed=case["seed"],
                gap_multiple=case["gap_multiple_of_pulse"],gap_index=gap_index+1,
                pulse_eligible=gap["primary_retention_pulse_eligible"],
                denominator_eligible=gap["primary_retention_denominator_eligible"],
                retention_ratio=gap["primary_retention_ratio"],success=gap["primary_retention_success"],
                sensitivity=gap["primary_retention_success_by_threshold"],
                supplemental_gain_ratio=gap["retention_ratio"]))
    summary=dict(generated_utc=campaign.now(),mode=manifest["mode"],workers=len(manifest["cases"]),
        kinds=kinds,all_density_conservation_guards_pass=valid,
        input_transform=campaign.read(out/"config.json").get("input_transform","lognormal"),
        all_hard_input_bounds_pass=bounds_valid,
        native_execution_ready=bool(valid and bounds_valid and status["all_scientific_artifact_identities_verified"]),
        manifest_sha256=campaign.sha(out/"manifest.json"),rows=rows,trajectories=trajectories,
        weak_scope="Weak means instantaneous bounded input relative to the native reference peak; long gaps can accumulate substantial commanded exposure, which is reported separately.",
        statistical_scope="Independent-seed ensemble is the unit of replication; each trajectory's gaps are not independent replicates.",
        field_tail_filters_retention=False,full_ensemble_complete=manifest["mode"]=="full",
        qualification_scope="Native execution/artifact/validity coverage; numerical sensitivity and scientific interpretation remain separate.")
    campaign.write(out/"results/noise_summary.json",summary)
    return {key:value for key,value in summary.items() if key not in ("rows","trajectories")}


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest="command",required=True)
    prep=sub.add_parser("prepare");prep.add_argument("--out",type=Path,required=True)
    prep.add_argument("--mode",choices=("pilot","bounded_pilot","full"),default="bounded_pilot")
    prep.add_argument("--transform",choices=("bounded_tanh","lognormal"),default="bounded_tanh")
    prep.add_argument("--calibrations",type=Path,required=True);prep.add_argument("--pilot-audit",type=Path)
    work=sub.add_parser("worker");work.add_argument("--out",type=Path,required=True);work.add_argument("--task",type=int,required=True)
    for name in ("verify","collect","submission-command"):
        command=sub.add_parser(name);command.add_argument("--out",type=Path,required=True)
        if name=="collect":
            command.add_argument("--job-id")
            command.add_argument("--verified-shards",type=Path)
    args=parser.parse_args();out=args.out.resolve()
    if args.command=="prepare":value=prepare(out,args.mode,args.calibrations,args.pilot_audit,args.transform)
    elif args.command=="verify":
        manifest=campaign.verify(out);value=dict(verified=True,workers=manifest["workers"],submitted=False)
    elif args.command=="collect":value=collect(out,args.job_id,args.verified_shards)
    elif args.command=="submission-command":
        campaign.verify(out);print((out/"submission_command.txt").read_text().strip());return
    else:
        try:value=worker(out,args.task)
        except Exception:
            failure=dict(generated_utc=campaign.now(),state="failed",task=args.task,
                scheduler=campaign.scheduler_context(),traceback=traceback.format_exc())
            campaign.write(out/"workers"/f"{args.task:04d}"/"failed.json",failure)
            campaign.write(out/"workers"/f"{args.task:04d}"/"status.json",failure)
            raise
    print(json.dumps(value,indent=2,allow_nan=False))


if __name__=="__main__":main()
