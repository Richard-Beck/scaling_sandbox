"""Immutable staged native Figure 2 campaign; submission is a separate action.

prepare writes frozen source/input/config closure and a SLURM command, never
submits. worker runs exactly one frozen case. collect requires every native
worker to finish and does not convert missing work into right censoring.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import shlex
import shutil
import socket
import subprocess
import sys
import time
import traceback

sys.dont_write_bytecode = True
for thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[thread_variable] = "1"

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / "polarity1d/manuscript").is_dir())
PACKAGE = ROOT / "polarity1d/manuscript"
CODE = PACKAGE / "code"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
LENGTHS = (19, 23, 27, 31, 35, 39)
NATIVE_PEAK = dict(wave_pinning=.1, goryachev=6.4, debelly=3.6,
                   holmes_model3=1., spring=4.)
RECOVERY = {model: 4000. if model == "debelly" else 10000. for model in MODELS}
CALIBRATION_HORIZON = dict(wave_pinning=10000., goryachev=10000.,
                          holmes_model3=10000., spring=100., debelly=2000.)
WIDTH_HORIZON = dict(wave_pinning=6000., goryachev=10000.,
                    holmes_model3=10000., spring=100., debelly=2000.)
PILOT_ON = dict(wave_pinning=100., goryachev=100., holmes_model3=100.,
                spring=1., debelly=50.)
STAGES = ("resource_pilot", "basals", "calibration", "width", "train_resource_pilot", "train")
TRAIN_PILOT_SPECS=((8,16.,2.),(4,4.,1.5))


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def allocated_bytes(path: Path) -> int:
    seen = set()
    total = 0
    for file in path.rglob("*"):
        if file.is_file() and not file.is_symlink():
            stat = file.stat()
            key = stat.st_dev, stat.st_ino
            if key not in seen:
                seen.add(key)
                total += stat.st_blocks * 512
    return total


def scientific_files(directory: Path):
    """Stable worker evidence, excluding producer/exception bookkeeping."""
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError("Scientific worker artifacts must not be symlinks")
        if not path.is_file() or path.name in ("status.json","result.json","failed.json"):
            continue
        if ".tmp" in path.name:
            raise ValueError("Completed worker retains an unfinished temporary artifact")
        yield path


def artifact_record(path: Path, root: Path) -> dict:
    resolved=path.resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError("Artifact escapes its frozen campaign")
    return dict(path=str(resolved.relative_to(root.resolve())),sha256=sha(resolved))


def stat_fingerprint(path: Path) -> dict:
    stat=path.stat()
    return dict(size=stat.st_size,mtime_ns=stat.st_mtime_ns,ctime_ns=stat.st_ctime_ns,
                device=stat.st_dev,inode=stat.st_ino)


def scheduler_context() -> dict:
    return dict(host=socket.gethostname(),pid=os.getpid(),
        slurm_job_id=os.environ.get("SLURM_JOB_ID"),
        slurm_array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID"),
        slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"))


def scheduler_states(job_id: str | None) -> dict:
    """Observe a supplied job handle; state files alone never imply liveness."""
    if job_id is None:
        return dict(available=False,reason="No authoritative job handle supplied",states={})
    if not all(character.isdigit() or character in "_," for character in job_id):
        raise ValueError("Invalid SLURM job handle")
    records={}
    errors=[]
    for command in (["sacct","-X","-n","-P","-j",job_id,"--format=JobID%64,State,ExitCode"],
                    ["squeue","-h","-j",job_id,"-o","%i|%T"]):
        try:
            completed=subprocess.run(command,capture_output=True,text=True,timeout=15,check=True)
            for line in completed.stdout.splitlines():
                columns=line.strip().split("|")
                if len(columns)>=2 and columns[0]:
                    records[columns[0]]=dict(state=columns[1].split()[0].rstrip("+"),
                        exit_code=columns[2] if len(columns)>2 else None,source=command[0])
        except (OSError,subprocess.SubprocessError) as error:
            errors.append(str(error))
    return dict(available=bool(records),job_id=job_id,observed_utc=now(),states=records,errors=errors)


def modules():
    directory = str(Path(__file__).resolve().parent)
    if directory not in sys.path:
        sys.path.insert(0, directory)
    return tuple(importlib.import_module(name) for name in (
        "signature_native_20261020", "signature_inputs_20261020",
        "signature_readouts_20261020", "signature_calibration_20261020"))


def config(*, calibration_horizon_overrides=None) -> dict:
    value = dict(recovery=dict(RECOVERY), calibration_horizon=dict(CALIBRATION_HORIZON),
        width_horizon=WIDTH_HORIZON, pilot_on=PILOT_ON,
        native_peak=NATIVE_PEAK, baseline_cells=201,
        recovery_clock_scope="Declared finite comparison horizons in each model's native units; not matched kinetics.",
        primary_dominance="D_half from disjoint reference half domains and current-minimum positive excess",
        post_withdrawal_replacement_threshold=.5,
        field_tail=dict(window=500., scaled_tolerance=.001),
        calibration=dict(relative_tolerance=.01, max_bisections=32,
            absolute_tolerance=dict(wave_pinning=.25,goryachev=.25,holmes_model3=.25,spring=.001,debelly=.5),
            scans=dict(wave_pinning=[1.,10.,100.,1000.,3000.,6000.,10000.],
                goryachev=[1.,10.,100.,1000.,3000.,6000.,10000.],
                holmes_model3=[1.,10.,100.,1000.,3000.,6000.,10000.],
                spring=[.01,.1,1.,10.,30.,60.,100.],
                debelly=[.5,2.,10.,50.,200.,800.,2000.])),
        spring_comparison=dict(resting_length="Lchemical/15", width="chemical_width/38",
            reference_area="native_peak*0.1", biological_calibration=False),
        cue=dict(shape="native endpoint raised cosine",reference_area="native_peak*3.8; Spring native_peak*0.1",
            normalization="exact discrete weighted commanded spatial area at every RHS evaluation"),
        noise_stage="separate later frozen payload; not implemented by this deterministic worker",
        retained_assays="Fresh local Stage 4 assays remain required; external 102 are unavailable.")
    import math
    for model, horizon in (calibration_horizon_overrides or {}).items():
        if model not in MODELS or not math.isfinite(horizon) or horizon < CALIBRATION_HORIZON[model]:
            raise ValueError("Extended calibration horizons must be finite and at least the original horizon")
        value["calibration_horizon"][model] = float(horizon)
        scans = value["calibration"]["scans"][model]
        candidate = 3. * CALIBRATION_HORIZON[model]
        while candidate < horizon:
            scans.append(candidate)
            candidate *= 2.
        value["calibration"]["scans"][model] = sorted(set(scans + [float(horizon)]))
    value["calibration_horizon_extension"] = dict(
        overrides=dict(calibration_horizon_overrides or {}),
        scope="Explicit finite-duration extension; native kinetics, cue and recovery remain unchanged.")
    return value


def _basal_records(preflight: dict) -> dict:
    result = {}
    for model, rows in preflight.get("model_basals", {}).items():
        for row in rows:
            length = row.get("requested_length")
            if length is None:
                length = row.get("matched_reference_length")
            if length is not None:
                result[f"{model}:L{int(length)}"] = row
    return result


def verified_calibration_collection(path: Path) -> dict:
    """Check an already-collected calibration against its exact worker receipts."""
    collection=read(path)
    source=path.resolve().parent.parent
    manifest=verify(source)
    if collection.get("manifest_sha256")!=sha(source/"manifest.json"):
        raise ValueError("Calibration collection and source manifest differ")
    status=read(source/"results/status.json")
    controls=read(source/"results/control_flags.json")
    audit=read(source/"results/artifact_verification.json")
    if not (status.get("complete") is True
            and status.get("all_scientific_artifact_identities_verified") is True
            and controls.get("all_observed_controls_pass") is True
            and len(controls["controls"])==len(manifest["cases"])):
        raise ValueError("Calibration source lacks complete artifact and basal-control verification")
    receipts={row["task"]:row for row in audit["workers"]}
    rows={row["calibration_id"]:row for row in collection["calibrations"]}
    if len(rows)!=len(manifest["cases"]):
        raise ValueError("Calibration collection coverage differs from frozen workers")
    for case in manifest["cases"]:
        result_path=source/"workers"/f"{case['task']:04d}"/"result.json"
        receipt=receipts.get(case["task"],{})
        if (not receipt.get("passed") or not receipt.get("all_scientific_artifacts_hashed")
                or sha(result_path)!=receipt.get("worker_result_sha256")):
            raise ValueError("Calibration worker no longer matches its artifact receipt")
        result=read(result_path)
        row=rows[case["basal_id"]]
        if (row["calibration"]!=result["calibration"]
                or row["basal_sha256"]!=result["basal"]["checkpoint_sha256"]
                or row["native_length"]!=case["native_length"]):
            raise ValueError("Calibration collection changed a native worker measurement")
    return collection


def prepare(out: Path, stage: str, *, basals_path: Path | None = None,
            calibrations_path: Path | None = None, models=None,
            reference_lengths=None, calibration_horizon=None) -> dict:
    if out.exists():
        raise FileExistsError("Choose a fresh campaign directory; frozen payloads are immutable.")
    if not out.resolve().is_relative_to((PACKAGE / "data/runs").resolve()):
        raise ValueError("Native campaigns belong under manuscript/data/runs.")
    # Importing the engine is a dependency/interface preflight, not integration.
    native, inputs, readouts, calibration = modules()
    for function in ("run_native", "run_phases", "read_checkpoint", "save_phase", "prepare_basal", "native_field_tail_gate"):
        if not callable(getattr(native, function, None)):
            raise RuntimeError("Native adapter is not complete: " + function)
    for module, function in ((inputs,"AreaCue"),(inputs,"make_train"),
                             (readouts,"analyze_trace"),(calibration,"calibrate")):
        if not callable(getattr(module,function,None)):
            raise RuntimeError("Incomplete source interface: " + function)
    from signature_campaign_plan_20261020 import compile_plan
    plan = compile_plan()
    source_proof=plan["native_preflight"].get("source_manifest",{})
    if not all(source_proof.get(key) is True for key in (
            "sha256_verified","terminal_complete","independent_validation_passed",
            "independent_audit_sha256_verified")):
        raise ValueError("Retained source campaign lacks complete verified qualification provenance")
    basals = _basal_records(plan["native_preflight"])
    external_inputs = {}
    if basals_path:
        external_inputs["basal_collection"] = dict(path=str(basals_path.resolve()),sha256=sha(basals_path))
        for record in read(basals_path).get("basals", []):
            if not record.get("qualified"):
                raise ValueError("Unqualified collected basal")
            basals[record["basal_id"]] = record
    calibration_records = {}
    if calibrations_path:
        external_inputs["calibration_collection"] = dict(path=str(calibrations_path.resolve()),sha256=sha(calibrations_path))
        calibration_records = {r["calibration_id"]:r for r in verified_calibration_collection(calibrations_path).get("calibrations",[])}
    if stage in ("train","train_resource_pilot") and not calibration_records:
        raise ValueError("Train preparation requires the complete qualified calibration collection.")
    selected_models=tuple(models or MODELS)
    selected_lengths=tuple(reference_lengths or LENGTHS)
    if calibration_horizon is not None and stage != "calibration":
        raise ValueError("A calibration horizon override is only valid for calibration workers")
    configuration = config(calibration_horizon_overrides=(
        {model: calibration_horizon for model in selected_models}
        if calibration_horizon is not None else None))
    if (len(set(selected_models))!=len(selected_models)
            or not set(selected_models)<=set(MODELS)
            or len(set(selected_lengths))!=len(selected_lengths)
            or not set(selected_lengths)<=set(LENGTHS)):
        raise ValueError("Selection must contain unique requested models/reference lengths")

    if stage == "resource_pilot":
        selected = [r for r in plan["calibrations"] if r["requested_reference_length"] in (19,39)]
    elif stage in ("basals", "calibration"):
        selected = plan["calibrations"]
    elif stage == "width":
        selected = plan["stage_1_equal_area_width"]
    elif stage in ("train","train_resource_pilot"):
        selected = plan["stage_2_cue_train"]
        if stage=="train_resource_pilot":
            selected=[row for row in selected if (row["pulse_count"],row["gap_multiple_of_pulse"],
                row["total_on_factor_of_calibration"]) in TRAIN_PILOT_SPECS]
    else:
        raise ValueError("Unsupported stage")
    selected=[r for r in selected if r["model"] in selected_models and (
        r["requested_reference_length"] in selected_lengths or
        (stage in ("basals","calibration") and r.get("role")=="stage_3_extra_L29"))]
    if not selected:
        raise ValueError("Selection contains no requested conditions")
    expected_main=(len(selected_models)*len(set(selected_lengths)&{19,39}) if stage=="resource_pilot"
        else len(selected_models)*len(selected_lengths)*(15 if stage=="width" else
            54 if stage=="train" else len(TRAIN_PILOT_SPECS) if stage=="train_resource_pilot" else 1))
    expected_extra=sum(model in ("wave_pinning","debelly") for model in selected_models) if stage in ("basals","calibration") else 0
    if len(selected)!=expected_main+expected_extra:
        raise ValueError("Selected grid does not match its explicit condition count")

    cases = []
    checkpoints = {}
    for row in selected:
        model = row["model"]
        reference_length = int(row["requested_reference_length"])
        basal_id = f"{model}:L{reference_length}"
        basal = basals.get(basal_id,{})
        case = dict(row,task=len(cases),stage=stage,basal_id=basal_id,
            native_length=reference_length/15. if model=="spring" else float(reference_length),
            cells=201,baseline_native_peak=NATIVE_PEAK[model],recovery_duration=RECOVERY[model],
            cue_width=.1 if model=="spring" else 3.8,
            cue_area=NATIVE_PEAK[model]*(.1 if model=="spring" else 3.8),
            checkpoint=None,checkpoint_sha256=None,allow_fresh_basal=stage in ("basals","resource_pilot"))
        checkpoint = basal.get("checkpoint")
        if checkpoint and (basal.get("status")=="verified" or basal.get("qualified") is True) and basal.get("sha256_verified",basal.get("qualified",False)):
            source = Path(checkpoint)
            if not source.is_absolute():
                source = ROOT / source
            expected = basal.get("actual_sha256") or basal.get("checkpoint_sha256") or basal.get("expected_sha256")
            if not source.is_file() or sha(source) != expected:
                raise ValueError("Exact basal checksum failed: " + basal_id)
            relative = f"inputs/basals/{model}_L{reference_length}.npz"
            checkpoints[relative] = dict(source=str(source.resolve()),sha256=expected)
            case.update(checkpoint=relative,checkpoint_sha256=expected,
                        fixed_scales=basal.get("fixed_scales"))
        elif not case["allow_fresh_basal"]:
            raise ValueError("Prepare and qualify missing exact basal first: " + basal_id)
        if stage == "width":
            case["cue_width"] = float(row["requested_width_label"])/(38. if model=="spring" else 1.)
            case["cue_area"] *= float(row["integrated_area_factor"])
            case["on_duration"] = WIDTH_HORIZON[model]
        elif stage == "resource_pilot":
            case["on_duration"] = PILOT_ON[model]
        elif stage in ("train","train_resource_pilot"):
            record = calibration_records.get(basal_id)
            if not record or not record["calibration"].get("train_ready"):
                raise ValueError("Missing/failed/nonmonotone calibration: " + basal_id)
            if record.get("basal_sha256") != case["checkpoint_sha256"]:
                raise ValueError("Train and calibration use different exact basal: " + basal_id)
            case["calibration_duration"] = record["calibration"]["upper_committing_duration"]
        cases.append(case)

    for directory in ("inputs/basals","frozen_sources","workers","logs","results"):
        (out/directory).mkdir(parents=True,exist_ok=True)
    input_hashes = {}
    for relative, record in checkpoints.items():
        target = out/relative
        shutil.copy2(record["source"],target)
        if sha(target) != record["sha256"]:
            raise RuntimeError("Input changed during freeze")
        input_hashes[relative] = record["sha256"]
    source_hashes = {}
    # Keep package-relative runtime layout. No trajectory/output payload enters
    # the frozen executable closure.
    for source in sorted(CODE.rglob("*")):
        if not source.is_file() or "__pycache__" in source.parts or source.suffix not in (".py",".json",".txt",".md"):
            continue
        relative = str(source.relative_to(ROOT))
        target = out/"frozen_sources"/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        before = sha(source)
        shutil.copy2(source,target)
        if sha(source) != before or sha(target) != before:
            raise RuntimeError("Source changed during freeze: " + relative)
        source_hashes[relative] = before
    for source in (PACKAGE/"data/inputs/fixed_scales.json",
                   PACKAGE/"data/runs/debelly_corrected_three_test_20261008/manifest.json"):
        relative=str(source.relative_to(ROOT))
        target=out/"frozen_sources"/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        before=sha(source)
        shutil.copy2(source,target)
        if sha(source)!=before or sha(target)!=before:
            raise RuntimeError("Runtime data changed during freeze: "+relative)
        source_hashes[relative]=before
    write(out/"config.json",configuration)
    input_hashes["config.json"] = sha(out/"config.json")
    for name, record in external_inputs.items():
        target=out/"inputs"/(name+".json")
        shutil.copy2(record["path"],target)
        if sha(target)!=record["sha256"]:
            raise RuntimeError("Collection changed during freeze")
        input_hashes[str(target.relative_to(out))]=record["sha256"]
    for source in (ROOT/"docs/tasks/briefs/20261020_figure2_signature_assays.txt",
                   ROOT/"docs/tasks/briefs/20261020_figure2_decisions.json"):
        target=out/"inputs"/source.name
        before=sha(source)
        shutil.copy2(source,target)
        if sha(source)!=before or sha(target)!=before:
            raise RuntimeError("Brief/decision changed during freeze")
        input_hashes[str(target.relative_to(out))]=before
    manifest = dict(created_utc=now(),stage=stage,status="frozen_not_submitted",result_schema=2,
        cases=cases,workers=len(cases),
        logical_protocols=2*len(cases) if stage=="resource_pilot" else len(cases),
        source_sha256=source_hashes,input_sha256=input_hashes,
        external_inputs=external_inputs,config_sha256=sha(out/"config.json"),
        selection=dict(models=list(selected_models),reference_lengths=list(selected_lengths),
            extra_noise_calibrations=expected_extra,expected_workers=expected_main+expected_extra,
            train_pilot_specs=list(TRAIN_PILOT_SPECS) if stage=="train_resource_pilot" else None,
            complete_model_length_coverage=(set(selected_models)==set(MODELS)
                and set(selected_lengths)==set(LENGTHS) and stage!="resource_pilot"),
            complete_brief_coverage=False,
            scope="This staged subset never replaces the complete five-model/six-length brief."),
        original_repository_root=str(ROOT),source_closure="package code sources, retained JSON/text/runtime metadata",
        native_outputs_local=True,automatic_submission=False,
        scope="New local assays; no external pilot results are claimed reproduced.")
    write(out/"manifest.json",manifest)
    write(out/"freeze.json",dict(manifest_sha256=sha(out/"manifest.json"),generated_utc=now()))
    worker = out/"worker.sh"
    frozen_script=out/"frozen_sources"/Path(__file__).relative_to(ROOT)
    worker.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec "+shlex.quote(PYTHON)+" -B "+shlex.quote(str(frozen_script))+" worker --out "+shlex.quote(str(out))+" --task \"$SLURM_ARRAY_TASK_ID\"\n")
    worker.chmod(0o750)
    memory="8G" if stage=="train_resource_pilot" else "4G" if stage=="train" else "2G"
    command=["sbatch","--parsable","--partition=red","--qos=small","--ntasks=1","--cpus-per-task=1","--mem="+memory,"--time=00:30:00",f"--array=0-{len(cases)-1}",
        "--job-name=F2-signature-"+stage,"--output="+str(out/"logs/%A_%a.out"),
        "--error="+str(out/"logs/%A_%a.err"),str(worker)]
    (out/"submission_command.txt").write_text(shlex.join(command)+"\n")
    write(out/"results/status.json",dict(generated_utc=now(),state="frozen_not_submitted",completed=0,total=len(cases)))
    verify(out)
    return dict(output=str(out),stage=stage,workers=len(cases),manifest_sha256=sha(out/"manifest.json"),submitted=False)


def verify(out: Path) -> dict:
    frozen=read(out/"freeze.json")
    if sha(out/"manifest.json")!=frozen["manifest_sha256"]:
        raise ValueError("Frozen manifest changed")
    manifest=read(out/"manifest.json")
    for relative,digest in manifest["source_sha256"].items():
        if sha(out/"frozen_sources"/relative)!=digest:
            raise ValueError("Frozen source changed: "+relative)
    for relative,digest in manifest["input_sha256"].items():
        if sha(out/relative)!=digest:
            raise ValueError("Frozen input changed: "+relative)
    if not manifest["cases"] or [c["task"] for c in manifest["cases"]]!=list(range(len(manifest["cases"]))):
        raise ValueError("Incomplete or unordered worker cases")
    return manifest


def _trace_mapping(trace,model: str) -> dict:
    import numpy as np
    get=lambda key:trace[key] if isinstance(trace,dict) else getattr(trace,key)
    result=dict(time=get("time"),x=get("x"),active=get("active"),sample_kind="nodes" if model=="spring" else "cells")
    if model=="spring":
        length=get("fields")["length"]
        result["physical_x"]=np.asarray(length)[:,None]*np.linspace(0.,1.,len(get("x")))[None,:]
    else:
        x=np.asarray(get("x"))
        dx=x[1]-x[0]
        result["reference_edges"]=np.r_[x[0]-.5*dx,x+.5*dx]
    return result


def _analyze(trace,basal,case,readouts) -> dict:
    reference_length=case["native_length"]
    fixed_width=15.2/38. if case["model"]=="spring" else 15.2
    basal_physical_length=float(basal.fields["length"][-1]) if case["model"]=="spring" else reference_length
    own_reference_width=case["cue_width"]*reference_length/basal_physical_length
    fixed_reference_width=fixed_width*reference_length/basal_physical_length
    parameters=dict(challenger_own_intervals=[(max(0.,reference_length-own_reference_width),reference_length)],
        challenger_fixed_intervals=[(max(0.,reference_length-fixed_reference_width),reference_length)],
        challenger_own_width_physical=case["cue_width"],
        challenger_side="right",denominator_atol=1e-12,contrast_atol=.05)
    return readouts.analyze_trace(_trace_mapping(trace,case["model"]),parameters,
                                 _trace_mapping(basal,case["model"]))


def _run_trial(native,inputs,readouts,case,state,basal,folder,on_duration=None,schedule=None,
               *,campaign_root,gap_input=None,max_step=None):
    model=case["model"]
    cue=inputs.AreaCue(case["cue_area"],case["cue_width"],grid_kind="node" if model=="spring" else "cell")
    zero=inputs.AreaCue(0.,case["cue_width"],grid_kind="node" if model=="spring" else "cell")
    phases=[]
    if schedule is not None:
        for index,(kind,start,end) in enumerate(schedule.segments):
            if end>start:
                stimulus=cue if kind=="on" else zero
                if kind=="gap" and gap_input is not None:
                    from signature_noise_20261020 import ShiftedInput
                    stimulus=ShiftedInput(gap_input,start)
                phases.append(dict(name=f"{index:02d}_{kind}",duration=end-start,cue=stimulus))
        schedule_metadata=schedule.as_dict()
    else:
        duration=float(on_duration or 0.)
        if duration>0:
            phases.append(dict(name="cue_on",duration=duration,cue=cue))
        schedule_metadata=dict(actual_total_on=duration,requested_total_on=duration)
    phases.append(dict(name="recovery",duration=case["recovery_duration"],cue=zero))
    trace,events=native.run_phases(model,case["native_length"],state,phases,
        output_dir=folder,cells=case["cells"],max_step=max_step)
    measurement=_analyze(trace,basal,case,readouts)
    tail=native.native_field_tail_gate(trace,model,case["native_length"],
        case.get("fixed_scales"),window=500.,tolerance=.001)
    write(folder/"readouts.json",measurement)
    write(folder/"native_tail.json",tail)
    import numpy as np
    times=np.asarray(measurement["series"]["time"])
    endpoints=[]
    for event in events:
        for label in ("start","stop"):
            target=float(event[label])
            hits=np.flatnonzero(np.isclose(times,target,rtol=0,atol=1e-8))
            if len(hits)!=1:
                raise ValueError("A native phase endpoint lacks one exact retained scalar sample")
            index=int(hits[0])
            endpoints.append(dict(phase=event["name"],endpoint=label,time=target,
                scalar_sample_index=index,
                readouts={key:values[index] for key,values in measurement["series"].items()}))
    on_intervals=[[float(event["start"]),float(event["stop"])] for event in events
        if event["name"]=="cue_on" or event["name"].endswith("_on")]
    recovery_events=[event for event in events if event["name"]=="recovery"]
    if len(recovery_events)!=1:
        raise ValueError("Exactly one declared final recovery phase is required")
    recovery_interval=[float(recovery_events[0]["start"]),float(recovery_events[0]["stop"])]
    write(folder/"phase_evidence.json",dict(on_intervals=on_intervals,
        recovery_interval=recovery_interval,endpoints=endpoints))
    result=dict(final_d=measurement["final_d"],readout_eligible=measurement["readout_eligible"],
        ineligibility_reason=measurement["ineligibility_reason"],
        field_tail_pass=tail["field_tail_pass"],
        field_tail=tail,schedule=schedule_metadata,events=events,
        trial_directory=str(folder.relative_to(campaign_root)),
        phase_on_intervals=on_intervals,recovery_interval=recovery_interval,
        phase_endpoint_readouts=endpoints,
        readouts=str((folder/"readouts.json").name),
        full_trace_sha256=sha(folder/"trace.npz"),
        artifacts={name:artifact_record(folder/filename,campaign_root) for name,filename in (
            ("readouts","readouts.json"),("trace","trace.npz"),
            ("native_tail","native_tail.json"),("phase_evidence","phase_evidence.json"))})
    if schedule is None and not on_duration:
        control=native.no_cue_control_gate(trace,basal,case["fixed_scales"],tolerance=.01)
        dominant=[v for v in measurement["series"]["D_half"] if v is not None]
        flipped=any(v>=.5 for v in dominant)
        result["no_cue_control"]=dict(control,spontaneous_D_replacement=flipped,
            passed=bool(not flipped and control["no_cue_control_pass"] and tail["field_tail_pass"]))
    return result


def worker(out: Path,task: int) -> dict:
    started=time.monotonic()
    manifest=verify(out)
    if not 0<=task<len(manifest["cases"]):
        raise ValueError("Worker task outside frozen array")
    case=manifest["cases"][task]
    directory=out/"workers"/f"{task:04d}"
    if (directory/"result.json").exists():
        raise FileExistsError("Refusing to overwrite a completed worker")
    directory.mkdir(parents=True,exist_ok=True)
    native,inputs,readouts,calibration=modules()
    context=scheduler_context()
    write(directory/"status.json",dict(generated_utc=now(),state="running",task=task,scheduler=context))
    if case["checkpoint"]:
        checkpoint=out/case["checkpoint"]
        basal_record=dict(qualified=True,qualification="exact previously qualified retained checkpoint",
                          fixed_scales=case["fixed_scales"])
    else:
        if not case["allow_fresh_basal"]:
            raise ValueError("Frozen case has no qualified exact basal")
        basal_record=native.prepare_basal(case["model"],case["native_length"],
            cells=case["cells"],output_dir=directory/"basal")
        if not basal_record["qualified"]:
            raise ValueError("Fresh exact basal did not qualify")
        checkpoint=Path(basal_record["checkpoint"])
    state,basal=native.read_checkpoint(checkpoint,case["model"])
    basal_record.update(basal_id=case["basal_id"],model=case["model"],
        requested_reference_length=case["requested_reference_length"],
        native_length=case["native_length"],checkpoint=str(checkpoint.resolve()),
        checkpoint_sha256=sha(checkpoint),actual_sha256=sha(checkpoint),sha256_verified=True,cells=case["cells"])
    case=dict(case,fixed_scales=basal_record.get("fixed_scales") or case.get("fixed_scales"))
    result=dict(generated_utc=now(),case=case,basal=basal_record,result_schema=2,scheduler=context,
                manifest_sha256=sha(out/"manifest.json"))
    if case["stage"]=="basals":
        result["zero_control"]=_run_trial(native,inputs,readouts,case,state,basal,directory/"zero",on_duration=0.,campaign_root=out)
        if not result["zero_control"]["no_cue_control"]["passed"]:
            write(directory/"zero_control_failed.json",result["zero_control"])
            raise ValueError("Basal failed its independent no-cue full-field/D control")
    elif case["stage"]=="resource_pilot":
        result["zero_control"]=_run_trial(native,inputs,readouts,case,state,basal,directory/"zero",on_duration=0.,campaign_root=out)
        if not result["zero_control"]["no_cue_control"]["passed"]:
            write(directory/"zero_control_failed.json",result["zero_control"])
            raise ValueError("Resource-pilot basal failed no-cue field/D control")
        result["challenge"]=_run_trial(native,inputs,readouts,case,state,basal,directory/"challenge",on_duration=case["on_duration"],campaign_root=out)
    elif case["stage"]=="calibration":
        configuration=read(out/"config.json")
        settings=configuration["calibration"]
        search=calibration.DurationSearch(horizon=configuration["calibration_horizon"][case["model"]],
            scan_durations=tuple(settings["scans"][case["model"]]),
            absolute_tolerance=settings["absolute_tolerance"][case["model"]],
            relative_tolerance=settings["relative_tolerance"],
            native_step=.5 if case["model"]=="debelly" else None,
            max_bisections=settings["max_bisections"])
        def evaluate(duration):
            folder=directory/"evaluations"/format(duration,".12g").replace(".","p")
            value=_run_trial(native,inputs,readouts,case,state,basal,folder,on_duration=duration,campaign_root=out)
            value["actual_duration"]=duration
            if duration==0. and not value["no_cue_control"]["passed"]:
                write(folder/"invalid_basal_control.json",value)
                raise ValueError("Calibration basal failed its no-cue full-field/D qualification")
            return value
        result["calibration"]=calibration.calibrate(search,evaluate)
    elif case["stage"]=="width":
        result["trial"]=_run_trial(native,inputs,readouts,case,state,basal,directory/"trial",on_duration=case["on_duration"],campaign_root=out)
    elif case["stage"] in ("train","train_resource_pilot"):
        schedule=inputs.make_train(case["calibration_duration"],case["total_on_factor_of_calibration"],
            case["pulse_count"],case["gap_multiple_of_pulse"],
            native_step=.5 if case["model"]=="debelly" else None)
        result["trial"]=_run_trial(native,inputs,readouts,case,state,basal,directory/"trial",schedule=schedule,campaign_root=out)
    else:
        raise ValueError("Unsupported worker stage")
    if "zero_control" in result:
        basal_record["no_cue_control"]=result["zero_control"]["no_cue_control"]
    elif "calibration" in result:
        basal_record["no_cue_control"]=next(e["no_cue_control"] for e in
            result["calibration"]["evaluations"] if e["duration"]==0.)
    result.update(wall_seconds=time.monotonic()-started,
        peak_resident_memory_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        allocated_bytes_before_final_records=allocated_bytes(directory))
    result["artifact_sha256"]={str(path.relative_to(out)):sha(path) for path in scientific_files(directory)}
    write(directory/"result.json",result)
    write(directory/"status.json",dict(generated_utc=now(),state="complete",task=task,scheduler=context))
    return dict(task=task,status="complete",wall_seconds=result["wall_seconds"])


def verify_worker_artifacts(out: Path,case: dict,record: dict,*,require_complete=True) -> dict:
    directory=out/"workers"/f"{case['task']:04d}"
    if record["manifest_sha256"]!=sha(out/"manifest.json"):
        raise ValueError("Worker manifest identity differs from frozen campaign")
    for key,value in case.items():
        if key=="fixed_scales" and value is None:
            continue
        if record["case"].get(key)!=value:
            raise ValueError("Worker changed frozen case field: "+key)
    checkpoint=Path(record["basal"]["checkpoint"])
    if not checkpoint.resolve().is_relative_to(out.resolve()):
        raise ValueError("Worker basal is outside frozen input/worker evidence")
    if sha(checkpoint)!=record["basal"]["checkpoint_sha256"]:
        raise ValueError("Worker exact basal hash mismatch")
    if case["checkpoint_sha256"] and record["basal"]["checkpoint_sha256"]!=case["checkpoint_sha256"]:
        raise ValueError("Worker substituted its retained basal")
    checks=0
    if "artifact_sha256" in record:
        expected=record["artifact_sha256"]
        actual_paths={str(path.relative_to(out)) for path in scientific_files(directory)}
        if actual_paths!=set(expected):
            raise ValueError("Worker scientific artifact inventory is incomplete or changed")
        for relative,digest in expected.items():
            path=out/relative
            if not path.resolve().is_relative_to(directory.resolve()) or sha(path)!=digest:
                raise ValueError("Worker scientific artifact changed: "+relative)
            checks+=1
        return dict(passed=True,all_scientific_artifacts_hashed=True,artifact_hashes_verified=checks,
                    worker_result_sha256=sha(directory/"result.json"))
    if require_complete:
        raise ValueError("Worker lacks the required complete scientific artifact hash inventory")
    # Legacy pilot histories have phase and combined native hashes, but no
    # original readout JSON identity. Preserve that limitation explicitly.
    def visit(value):
        nonlocal checks
        if isinstance(value,dict):
            if "archive" in value and isinstance(value["archive"],dict):
                archive=value["archive"]
                if sha(Path(archive["path"]))!=archive["sha256"]:
                    raise ValueError("Legacy native phase hash mismatch")
                checks+=1
            for nested in value.values():visit(nested)
        elif isinstance(value,list):
            for nested in value:visit(nested)
    visit(record)
    legacy_trials=[]
    for label,folder_name in (("zero_control","zero"),("challenge","challenge"),("trial","trial")):
        if label in record:legacy_trials.append((record[label],directory/folder_name))
    for evaluation in record.get("calibration",{}).get("evaluations",[]):
        folder=directory/"evaluations"/format(evaluation["duration"],".12g").replace(".","p")
        legacy_trials.append((evaluation,folder))
    for trial,folder in legacy_trials:
        if "full_trace_sha256" in trial:
            if sha(folder/"trace.npz")!=trial["full_trace_sha256"]:
                raise ValueError("Legacy combined native trace hash mismatch")
            checks+=1
    return dict(passed=True,all_scientific_artifacts_hashed=False,artifact_hashes_verified=checks,
        unhashed_legacy_evidence="Readout JSON and some compact records require independent reconstruction from hashed native fields.",
        worker_result_sha256=sha(directory/"result.json"))


def audit_shard(out: Path,index: int,shards: int) -> dict:
    """Hash a disjoint worker subset once and retain exact stat snapshots."""
    if shards<1 or not 0<=index<shards:raise ValueError("Invalid audit shard")
    manifest=verify(out);rows=[];missing=[]
    for case in manifest["cases"]:
        if case["task"]%shards!=index:continue
        directory=out/"workers"/f"{case['task']:04d}";result_path=directory/"result.json"
        if not result_path.is_file():missing.append(case["task"]);continue
        result=read(result_path)
        paths=list(scientific_files(directory))+[result_path]
        before={str(path.relative_to(out)):stat_fingerprint(path) for path in paths}
        receipt=verify_worker_artifacts(out,case,result,require_complete=True)
        after={str(path.relative_to(out)):stat_fingerprint(path) for path in paths}
        if before!=after:raise ValueError("Scientific artifact changed during shard hashing")
        snapshots={relative:dict(fingerprint=metadata,
            sha256=(receipt["worker_result_sha256"] if relative.endswith("/result.json") else result["artifact_sha256"][relative]))
            for relative,metadata in after.items()}
        rows.append(dict(task=case["task"],**receipt,artifact_snapshots=snapshots))
    report=dict(generated_utc=now(),manifest_sha256=sha(out/"manifest.json"),
        driver_sha256=sha(Path(__file__)),index=index,shards=shards,workers=rows,
        missing_tasks=missing,complete=not missing,
        scope="Disjoint native artifact hash pass with pre/post stat identity; no native integration or readout recomputation.")
    target=out/"results/artifact_shards"/f"shard_{index:03d}.json"
    write(target,report)
    return dict(index=index,shards=shards,complete=report["complete"],workers=len(rows),output=str(target))


def verified_shard_receipts(out: Path,directory: Path) -> dict:
    reports=[read(path) for path in sorted(directory.glob("shard_*.json"))]
    if not reports:raise ValueError("No scientific artifact shards found")
    count=reports[0]["shards"]
    if len(reports)!=count or {r["index"] for r in reports}!=set(range(count)):
        raise ValueError("Scientific artifact shard coverage is incomplete")
    rows={}
    for report in reports:
        if not report["complete"] or report["manifest_sha256"]!=sha(out/"manifest.json"):
            raise ValueError("Incomplete/mismatched scientific artifact shard")
        for row in report["workers"]:
            if row["task"] in rows:raise ValueError("Overlapping scientific artifact shards")
            directory_path=out/"workers"/f"{row['task']:04d}"
            actual={str(path.relative_to(out)) for path in scientific_files(directory_path)}|{str((directory_path/"result.json").relative_to(out))}
            if actual!=set(row["artifact_snapshots"]):raise ValueError("Artifact inventory changed after shard audit")
            for relative,snapshot in row["artifact_snapshots"].items():
                # NFS st_dev identifies a client's local mount, so the same
                # exported file has different device numbers on compute and
                # login nodes. The canonical campaign path plus inode, size
                # and both timestamps identify the unchanged hashed artifact.
                current=stat_fingerprint(out/relative)
                expected=snapshot["fingerprint"]
                if any(current[key]!=expected[key] for key in
                       ("size","mtime_ns","ctime_ns","inode")):
                    raise ValueError("Scientific artifact stat identity changed after hashing: "+relative)
            rows[row["task"]]=dict(row,hash_proof_reused_from_shard=True,
                stat_identity_scope="Canonical campaign path, inode, size, mtime and ctime; NFS client-local device number excluded.")
    return rows


def collect(out: Path,job_id=None,verified_shards=None) -> dict:
    manifest=verify(out)
    records=[]
    missing=[]
    failed=[]
    verification=[]
    controls=[]
    scheduler=scheduler_states(job_id)
    shard_receipts=verified_shard_receipts(out,Path(verified_shards)) if verified_shards else None
    for case in manifest["cases"]:
        directory=out/"workers"/f"{case['task']:04d}"
        path=directory/"result.json"
        if path.is_file():
            record=read(path)
            if shard_receipts is None:
                verification.append(dict(task=case["task"],**verify_worker_artifacts(out,case,record,
                    require_complete=manifest.get("result_schema",1)>=2)))
            else:
                if case["task"] not in shard_receipts:raise ValueError("Missing requested worker artifact shard receipt")
                verification.append(shard_receipts[case["task"]])
            if "zero_control" in record:
                controls.append(dict(task=case["task"],basal_id=case["basal_id"],
                    role=case["stage"]+"_basal_qualification",**record["zero_control"]["no_cue_control"]))
            if "calibration" in record:
                zero=[r for r in record["calibration"]["evaluations"] if r["duration"]==0.]
                if len(zero)!=1 or "no_cue_control" not in zero[0]:
                    raise ValueError("Calibration lacks one explicit no-cue basal qualification")
                controls.append(dict(task=case["task"],basal_id=case["basal_id"],
                    role="calibration_basal_qualification",**zero[0]["no_cue_control"]))
            records.append(record)
        elif (directory/"failed.json").is_file():
            failed.append(dict(task=case["task"],failure=read(directory/"failed.json")))
        else:
            missing.append(case["task"])
    terminal_states={"FAILED","CANCELLED","TIMEOUT","OUT_OF_MEMORY","NODE_FAIL","PREEMPTED","BOOT_FAIL","DEADLINE","REVOKED"}
    terminal_failures=[]
    live_tasks=[]
    terminal_missing=[]
    for task in missing+[entry["task"] for entry in failed]:
        handle=f"{job_id}_{task}" if job_id and "_" not in job_id else job_id
        observed=scheduler["states"].get(handle,{})
        if observed.get("state") in terminal_states:terminal_failures.append(task)
        elif observed.get("state")=="COMPLETED":terminal_missing.append(task)
        elif observed.get("state") in {"RUNNING","PENDING","COMPLETING","CONFIGURING","SUSPENDED"}:live_tasks.append(task)
    status=dict(generated_utc=now(),stage=manifest["stage"],completed=len(records),
        total=len(manifest["cases"]),missing_tasks=missing,recorded_failed_tasks=[r["task"] for r in failed],
        complete=len(records)==len(manifest["cases"]),
        state="complete" if len(records)==len(manifest["cases"]) else "recorded_failures" if failed else "incomplete",
        scheduler=scheduler,
        terminal_failed_tasks=terminal_failures,live_tasks=live_tasks,
        scheduler_complete_but_result_missing_tasks=terminal_missing,
        artifact_hashes_verified=sum(r["artifact_hashes_verified"] for r in verification),
        all_scientific_artifact_identities_verified=bool(records) and all(r["all_scientific_artifacts_hashed"] for r in verification),
        allocated_bytes=allocated_bytes(out),measurement_scope="This campaign root, local allocated file blocks; no directory blocks.")
    write(out/"results/artifact_verification.json",dict(generated_utc=now(),workers=verification,
        recorded_failures=failed,scope="Scientific file identity; scientific interpretation requires independent native/readout checks."))
    write(out/"results/control_flags.json",dict(generated_utc=now(),controls=controls,
        all_observed_controls_pass=all(c["passed"] for c in controls),
        response_field_tail_filters_replacement=False,
        scope="No-cue controls qualify the basal; response field stability remains independent of D replacement."))
    write(out/"results/status.json",status)
    if not status["complete"]:
        return status
    basals={r["basal"]["basal_id"]:r["basal"] for r in records}
    write(out/"results/basals.json",dict(generated_utc=now(),manifest_sha256=sha(out/"manifest.json"),basals=list(basals.values())))
    if manifest["stage"]=="calibration":
        calibrations=[dict(calibration_id=r["case"]["basal_id"],model=r["case"]["model"],
            requested_reference_length=r["case"]["requested_reference_length"],
            native_length=r["case"]["native_length"],basal_sha256=r["basal"]["checkpoint_sha256"],
            calibration=r["calibration"]) for r in records]
        write(out/"results/calibrations.json",dict(generated_utc=now(),
            manifest_sha256=sha(out/"manifest.json"),calibrations=calibrations,
            all_train_ready=all(r["calibration"]["train_ready"] for r in calibrations)))
        unavailable=[]
        for record in records:
            case=record["case"];calibration=record["calibration"]
            if case["requested_reference_length"] not in LENGTHS or calibration["train_ready"]:
                continue
            result_path=out/"workers"/f"{case['task']:04d}"/"result.json"
            receipt=next(row for row in verification if row["task"]==case["task"])
            for pulses in (2,4,8):
                for gap in (0.,.1,.5,1.,4.,16.):
                    for total_on_factor in (1.,1.5,2.):
                        unavailable.append(dict(model=case["model"],
                            requested_reference_length=case["requested_reference_length"],
                            native_length=case["native_length"],pulse_count=pulses,
                            gap_multiple_of_pulse=gap,total_on_factor_of_calibration=total_on_factor,
                            status="calibration_unavailable",reason=calibration["status"],
                            searched_on_cap=calibration["finite_duration_horizon"],
                            recovery_duration=case["recovery_duration"],
                            exact_basal_sha256=record["basal"]["checkpoint_sha256"],
                            calibration_result_path=str(result_path.resolve()),
                            calibration_result_sha256=receipt["worker_result_sha256"],
                            calibration_artifact_hashes_verified=receipt["all_scientific_artifacts_hashed"],
                            calibration_manifest_sha256=sha(out/"manifest.json"),
                            calibration_duration=None,train_executed=False))
        write(out/"results/unavailable_trains.json",dict(generated_utc=now(),
            manifest_sha256=sha(out/"manifest.json"),conditions=unavailable,
            unavailable_conditions=len(unavailable),native_missing_or_failed=False,
            scope="Dependent standardized trains are undefined without a qualified calibration; this is not an integration failure or an infinite-time impossibility claim."))
    write(out/"results/resources.json",dict(generated_utc=now(),workers=len(records),
        summed_wall_seconds=sum(r["wall_seconds"] for r in records),
        maximum_wall_seconds=max(r["wall_seconds"] for r in records),
        maximum_peak_resident_memory_KiB=max(r["peak_resident_memory_KiB"] for r in records),
        allocated_bytes=allocated_bytes(out),measurement_utc=now(),
        scope="Measured completed worker resources; not a guarantee for longer/new assay arms."))
    return status


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="command",required=True)
    prep=sub.add_parser("prepare")
    prep.add_argument("--out",type=Path,required=True)
    prep.add_argument("--stage",choices=STAGES,required=True)
    prep.add_argument("--basals",type=Path)
    prep.add_argument("--calibrations",type=Path)
    prep.add_argument("--models",nargs="+",choices=MODELS,
                      help="Freeze a staged model subset; full brief coverage remains required.")
    prep.add_argument("--reference-lengths",nargs="+",type=int,choices=LENGTHS,
                      help="Requested six-length subset; calibration adds selected WP/deBelly L29.")
    prep.add_argument("--calibration-horizon",type=float,
                      help="Explicit finite ON-horizon extension for selected calibration models.")
    for name in ("verify","collect","submission-command"):
        command=sub.add_parser(name)
        command.add_argument("--out",type=Path,required=True)
        if name=="collect":
            command.add_argument("--job-id",help="Live/terminal SLURM handle for authoritative producer-state observation.")
            command.add_argument("--verified-shards",type=Path,help="Reuse complete hash shards after exact current file-stat checks.")
    audit=sub.add_parser("audit-shard")
    audit.add_argument("--out",type=Path,required=True)
    audit.add_argument("--index",type=int,required=True)
    audit.add_argument("--shards",type=int,required=True)
    work=sub.add_parser("worker")
    work.add_argument("--out",type=Path,required=True)
    work.add_argument("--task",type=int,required=True)
    args=parser.parse_args()
    out=args.out.resolve()
    if args.command=="prepare":
        result=prepare(out,args.stage,basals_path=args.basals,calibrations_path=args.calibrations,
                       models=args.models,reference_lengths=args.reference_lengths,
                       calibration_horizon=args.calibration_horizon)
    elif args.command=="verify":
        manifest=verify(out)
        result=dict(verified=True,stage=manifest["stage"],workers=manifest["workers"],submitted=False)
    elif args.command=="collect":
        result=collect(out,job_id=args.job_id,verified_shards=args.verified_shards)
        result=dict(result)
        observation=result.pop("scheduler")
        result["scheduler_summary"]=dict(available=observation["available"],
            job_id=observation.get("job_id"),observed_utc=observation.get("observed_utc"),
            observed_records=len(observation["states"]),errors=observation.get("errors",[]))
    elif args.command=="submission-command":
        verify(out)
        print((out/"submission_command.txt").read_text().strip())
        return
    elif args.command=="audit-shard":
        result=audit_shard(out,args.index,args.shards)
    else:
        try:
            result=worker(out,args.task)
        except Exception:
            context=scheduler_context()
            write(out/"workers"/f"{args.task:04d}"/"failed.json",dict(generated_utc=now(),
                state="failed",task=args.task,scheduler=context,traceback=traceback.format_exc()))
            write(out/"workers"/f"{args.task:04d}"/"status.json",dict(generated_utc=now(),
                state="failed",task=args.task,scheduler=context))
            raise
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
