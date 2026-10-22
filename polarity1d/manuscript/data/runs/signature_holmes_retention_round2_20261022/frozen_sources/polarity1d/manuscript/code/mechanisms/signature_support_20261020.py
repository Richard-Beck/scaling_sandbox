"""Fresh local S4 assays: position, sampled-share withdrawal and initiation.

External pilot/102-run results are unavailable and are not reconstructed here.
``prepare`` freezes a reviewable payload and prints commands; it never submits.
Native execution belongs in SLURM after review. Existing scientific inputs and
the main boundary-cue comparison are unchanged.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shlex
import shutil
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import signature_inputs_20261020 as inputs
import signature_readouts_20261020 as readouts
from signature_campaign_20261020 import NATIVE_PEAK, RECOVERY, WIDTH_HORIZON

PACKAGE = Path(__file__).resolve().parents[2]
ROOT = PACKAGE.parent.parent
FOCAL_MODELS = ("wave_pinning", "goryachev", "debelly")
LENGTHS = (19, 23, 27, 31, 35, 39)
SHARES = (.5, .75, .9)
DWELL_FACTORS = (0., .25, 1.)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


class CenteredAreaCue:
    """New full-support centered cosine, with exact native commanded area."""
    event_times = ()

    def __init__(self, area, width, center):
        self.area, self.width, self.center = map(float, (area, width, center))
        if not all(math.isfinite(v) for v in (self.area, self.width, self.center)) or self.area < 0 or self.width <= 0:
            raise ValueError("Finite nonnegative area, positive width and finite center required")

    def profile(self, x, length):
        x = np.asarray(x, dtype=float)
        weights = inputs.integration_weights(x, length, "cell")
        distance = 2 * np.abs(x - self.center) / self.width
        raw = np.where(distance < 1, .5 * (1 + np.cos(np.pi * np.clip(distance, 0, 1))), 0.)
        raw_area = float(raw @ weights)
        if raw_area <= 0 and self.area > 0:
            raise ValueError("Native grid does not resolve the interior cue")
        multiplier = self.area / raw_area if self.area else 0.
        return raw * multiplier, weights, multiplier

    def __call__(self, time, x, length):
        del time
        return self.profile(x, length)[0]

    def profile_diagnostics(self, x, length):
        profile, weights, peak = self.profile(x, length)
        return dict(shape="new full-support centered cosine", center=self.center,
                    requested_width=self.width, requested_area=self.area,
                    sampled_area=float(profile @ weights), sampled_peak=float(profile.max()),
                    reconstructed_peak=peak, support=[self.center - self.width / 2, self.center + self.width / 2],
                    support_truncated=bool(self.center - self.width / 2 < 0 or self.center + self.width / 2 > length),
                    area_definition="Spatial native commanded area per ON time; no cumulative-dose calibration")


class RelativePeakCue:
    """Relative-20% boundary support at a fixed commanded peak; dose varies."""
    event_times = ()

    def __init__(self, peak, fraction=.2):
        self.peak, self.fraction = float(peak), float(fraction)
        if not math.isfinite(self.peak) or self.peak < 0 or not 0 < self.fraction <= 1:
            raise ValueError("Nonnegative finite peak and support fraction in (0,1] required")

    def __call__(self, time, x, length):
        del time
        distance = (length - np.asarray(x)) / (length * self.fraction)
        return np.where((distance >= 0) & (distance <= 1),
                        self.peak * .5 * (1 + np.cos(np.pi * np.clip(distance, 0, 1))), 0.)

    def profile_diagnostics(self, x, length):
        profile = self(0, x, length)
        weights = inputs.integration_weights(x, length, "cell")
        return dict(shape="relative-20% endpoint cosine; fixed peak", requested_width=self.fraction * length,
                    reconstructed_boundary_peak=self.peak, sampled_peak=float(profile.max()),
                    sampled_area=float(profile @ weights), support_truncated=False,
                    area_definition="Dose-changing supporting diagnostic; not an equal-area shape contrast")


def position_centers(length, width=3.8):
    left, right = .5 * length + .5 * width, length - .5 * width
    if left > right:
        raise ValueError("The full centered support cannot fit inside the challenger half")
    return np.linspace(left, right, 3).tolist()


def actual_dwell(model, requested):
    requested = float(requested)
    if not math.isfinite(requested) or requested < 0:
        raise ValueError("Finite nonnegative additional cue dwell required")
    return math.ceil(requested / .5 - 1e-12) * .5 if model == "debelly" else requested


def _trace_mapping(trace):
    x = np.asarray(trace.x)
    dx = x[1] - x[0]
    return dict(time=trace.time, x=x, active=trace.active, sample_kind="cells",
                reference_edges=np.r_[x[0] - .5 * dx, x + .5 * dx])


def _parameters(case, basal):
    length, width = case["native_length"], case["cue_width"]
    support = ([case["cue_center"] - width / 2, case["cue_center"] + width / 2]
               if case["assay"] == "position" else [length - width, length])
    return dict(challenger_own_intervals=[support],
                challenger_fixed_intervals=[[max(0., length - 15.2), length]],
                challenger_side="right", denominator_atol=1e-12, contrast_atol=.05)


def _incumbent_tenth_support(basal):
    mapping = _trace_mapping(basal)
    x, active = mapping["x"], np.asarray(basal.active[-1])
    excess = active - active.min()
    if excess.max() <= 0:
        raise ValueError("Homogeneous basal has no incumbent 10% support")
    selected = excess >= .1 * excess.max()
    peak = int(np.argmax(excess))
    lo = hi = peak
    while lo > 0 and selected[lo - 1]:
        lo -= 1
    while hi + 1 < len(x) and selected[hi + 1]:
        hi += 1
    edges = mapping["reference_edges"]
    return [float(edges[lo]), float(edges[hi + 1])]


def share_event(measurement, target):
    """Locate a first eligible upper sample; never interpolate a native state."""
    from signature_analysis_20261020 import threshold_event
    series = measurement["series"]
    event = threshold_event(series, series["time"][0], series["time"][-1], target,
                            value_key="share_half", eligibility_key="D_half_eligible")
    if event.get("observed"):
        time = np.asarray(series["time"])
        indices = np.flatnonzero(np.abs(time - event["upper_time"]) <= 1e-10)
        if len(indices) != 1:
            raise ValueError("Observed share endpoint has no exact native sample")
        event.update(upper_sample_index=int(indices[0]), actual_upper_sample_share=series["share_half"][int(indices[0])],
                     native_state_interpolated=False)
    return event


def build_plan(calibration_root, width_root):
    calibration_root, width_root = Path(calibration_root).resolve(), Path(width_root).resolve()
    basals_file = calibration_root / "results/basals.json"
    calibration_file = calibration_root / "results/calibrations.json"
    basals = {row["basal_id"]: row for row in read(basals_file)["basals"]}
    calibration = {row["calibration_id"]: row for row in read(calibration_file)["calibrations"]}
    width_manifest = read(width_root / "manifest.json")
    prefixes, cases, reused = [], [], []
    for model in FOCAL_MODELS:
        for length in LENGTHS:
            identifier = f"{model}:L{length}"
            basal, calibrated = basals[identifier], calibration[identifier]
            checkpoint = Path(basal["checkpoint"])
            digest = sha(checkpoint)
            if not basal["qualified"] or digest != basal["checkpoint_sha256"] or digest != calibrated["basal_sha256"]:
                raise ValueError("Unqualified or mismatched exact basal: " + identifier)
            if not calibrated["calibration"]["train_ready"]:
                raise ValueError("Unresolved/nonmonotone commitment calibration: " + identifier)
            common = dict(model=model, native_length=float(length), requested_reference_length=length,
                stage="supporting_assays",
                basal_id=identifier, cells=201, checkpoint=str(checkpoint), checkpoint_sha256=digest,
                fixed_scales=basal["fixed_scales"], cue_width=3.8,
                cue_area=NATIVE_PEAK[model] * 3.8, recovery_duration=RECOVERY[model],
                on_duration=WIDTH_HORIZON[model],
                calibration_duration=calibrated["calibration"]["upper_committing_duration"],
                calibration_collection_sha256=sha(calibration_file))
            prefix_task = len(prefixes)
            prefixes.append(dict(common, task=prefix_task, assay="share_prefix", shares=list(SHARES)))
            for center in position_centers(length):
                cases.append(dict(common, task=len(cases), assay="position", cue_center=center))
            for target in SHARES:
                for factor in DWELL_FACTORS:
                    requested = factor * common["calibration_duration"]
                    cases.append(dict(common, task=len(cases), assay="share_withdrawal",
                        prefix_task=prefix_task, target_share=target, dwell_factor_of_calibration=factor,
                        requested_additional_dwell=requested, additional_dwell=actual_dwell(model, requested)))
            for factor in (.5, 1., 2.):
                cases.append(dict(common, task=len(cases), assay="relative_width_initiation",
                    cue_width=.2 * length, cue_peak=2 * NATIVE_PEAK[model] * factor,
                    strength_factor=factor, diagnostic_dose_changes=True))
                matches = [row for row in width_manifest["cases"] if row["model"] == model
                    and row["requested_reference_length"] == length
                    and row["requested_width_label"] == 3.8
                    and row["integrated_area_factor"] == factor]
                if len(matches) != 1:
                    raise ValueError("Missing/duplicate retained fixed-width strength arm")
                path = width_root / "workers" / f"{matches[0]['task']:04d}" / "result.json"
                result = read(path)
                if result["basal"]["checkpoint_sha256"] != digest or not result["trial"]["readout_eligible"]:
                    raise ValueError("Retained fixed-width result has wrong basal or unqualified readout")
                reused.append(dict(model=model, requested_reference_length=length, strength_factor=factor,
                    result=str(path), result_sha256=sha(path),
                    field_tail_pass=result["trial"]["field_tail_pass"],
                    qualification="Reuse eligible scalar readout regardless of separate finite-tail result"))
    return dict(schema_version=1, scientific_status="Proposed fresh local S4 implementation; no results implied",
        source_collections={str(basals_file): sha(basals_file), str(calibration_file): sha(calibration_file),
                            str(width_root / "manifest.json"): sha(width_root / "manifest.json")},
        counts=dict(shared_prefixes=18, position_trials=54, share_withdrawal_trials=162,
                    relative_width_initiation_trials=54, fresh_trials=270, reused_fixed_width_strength=54),
        prefix_cases=prefixes, cases=cases, reused_stage1_results=reused,
        scope="Three focal models and six core lengths; five-model class scope comes from main stages 1-2.",
        interpretation_limits=["Position uses a new centered full-support cue; not the main boundary cue or a geometry-only causal partition.",
            "Location contrasts hold size/native basal/kinetics fixed; length changes multiple basal properties.",
            "Threshold withdrawal restarts an exact first eligible upper sample, not an interpolated threshold state.",
            "Censored share targets cause no fabricated withdrawal; all native-field-tail flags remain separate.",
            "Relative-width initiation holds peak rather than area; dose changes explicitly."])


def _native():
    import signature_native_20261020
    return signature_native_20261020


def _write_trial(folder, case, state, basal, cue, *, on_duration=None):
    native = _native()
    phases = []
    on_duration = case["on_duration"] if on_duration is None else on_duration
    if on_duration > 0:
        phases.append(dict(name="cue_on", duration=on_duration, cue=cue))
    phases.append(dict(name="recovery", duration=case["recovery_duration"],
                       cue=inputs.AreaCue(0., case["cue_width"])))
    trace, events = native.run_phases(case["model"], case["native_length"], state,
        phases, output_dir=folder, cells=case["cells"])
    parameters = _parameters(case, basal)
    measured = readouts.analyze_trace(_trace_mapping(trace), parameters, _trace_mapping(basal))
    tail = native.native_field_tail_gate(trace, case["model"], case["native_length"],
        case["fixed_scales"], window=500., tolerance=.001)
    write(folder / "readouts.json", measured)
    write(folder / "native_tail.json", tail)
    support = parameters["challenger_own_intervals"][0]
    incumbent = _incumbent_tenth_support(basal)
    gap = support[0] - incumbent[1]
    overlap = max(0., min(support[1], incumbent[1]) - max(support[0], incumbent[0]))
    return dict(final_d=measured["final_d"], readout_eligible=measured["readout_eligible"],
        ineligibility_reason=measured["ineligibility_reason"], final_replacement=measured["final_replacement"],
        field_tail_pass=tail["field_tail_pass"], field_tail=tail, events=events,
        phase_on_intervals=[[e["start"], e["stop"]] for e in events if e["name"] == "cue_on"],
        recovery_interval=[events[-1]["start"], events[-1]["stop"]],
        readout_parameters=parameters, incumbent_10pct_support=incumbent,
        cue_support=support, incumbent_10pct_edge_to_cue_gap=gap,
        incumbent_10pct_cue_overlap_extent=overlap,
        coordinate_scope="Nominal cell-centred solver coordinates; deBelly mechanics retained independently",
        phase_readout_endpoints=[dict(name=e["name"], start=e["start"], stop=e["stop"],
            first=next(i for i, t in enumerate(measured["series"]["time"]) if abs(t - e["start"]) <= 1e-9),
            last=next(i for i, t in enumerate(measured["series"]["time"]) if abs(t - e["stop"]) <= 1e-9)) for e in events])


def run_prefix(case, folder):
    native = _native()
    folder = Path(folder)
    if (folder / "result.json").exists():
        raise FileExistsError("Refusing to overwrite completed shared prefix")
    if sha(case["checkpoint"]) != case["checkpoint_sha256"]:
        raise ValueError("Exact basal changed")
    state, basal = native.read_checkpoint(case["checkpoint"], case["model"])
    cue = inputs.AreaCue(case["cue_area"], case["cue_width"])
    trace, events = native.run_phases(case["model"], case["native_length"], state,
        [dict(name="cue_on", duration=case["on_duration"], cue=cue)],
        output_dir=folder, cells=case["cells"])
    measured = readouts.analyze_trace(_trace_mapping(trace), _parameters(case, basal), _trace_mapping(basal))
    write(folder / "readouts.json", measured)
    result = dict(case=case, events=events, shares={str(target): share_event(measured, target) for target in SHARES},
                  trace_sha256=sha(folder / "trace.npz"), readouts_sha256=sha(folder / "readouts.json"))
    result["artifact_sha256"] = {str(p.relative_to(folder)): sha(p) for p in folder.rglob("*")
                                if p.is_file() and p.name != "result.json"}
    write(folder / "result.json", result)
    return result


def run_case(case, folder, prefix_folder=None):
    native = _native()
    folder = Path(folder)
    if (folder / "result.json").exists():
        raise FileExistsError("Refusing to overwrite completed support case")
    if sha(case["checkpoint"]) != case["checkpoint_sha256"]:
        raise ValueError("Exact basal changed")
    state, basal = native.read_checkpoint(case["checkpoint"], case["model"])
    shared = None
    if case["assay"] == "share_withdrawal":
        if prefix_folder is None:
            raise ValueError("Share withdrawal requires its completed shared prefix")
        prefix_folder = Path(prefix_folder)
        shared = read(prefix_folder / "result.json")
        if (sha(prefix_folder / "trace.npz") != shared["trace_sha256"] or
                sha(prefix_folder / "readouts.json") != shared["readouts_sha256"] or
                shared["case"]["checkpoint_sha256"] != case["checkpoint_sha256"]):
            raise ValueError("Shared exact native prefix identity mismatch")
        event = shared["shares"][str(case["target_share"])]
        recomputed_event = share_event(read(prefix_folder / "readouts.json"), case["target_share"])
        if event != recomputed_event:
            raise ValueError("Shared prefix threshold event differs from its recorded native readout")
        for key in ("model", "native_length", "cue_width", "cue_area", "on_duration"):
            if shared["case"][key] != case[key]:
                raise ValueError("Shared prefix protocol mismatch: " + key)
        if not event.get("observed"):
            result = dict(case=case, status="share_target_not_observed_in_finite_native_prefix",
                          withdrawal_executed=False, prefix_event=event,
                          prefix_result_sha256=sha(prefix_folder / "result.json"))
            write(folder / "result.json", result)
            return result
        _, prefix_trace = native.read_checkpoint(prefix_folder / "trace.npz", case["model"])
        state = native.state_at(prefix_trace, case["model"], event["upper_sample_index"])
        cue = inputs.AreaCue(case["cue_area"], case["cue_width"])
        result = _write_trial(folder, case, state, basal, cue, on_duration=case["additional_dwell"])
        result.update(prefix_event=event, prefix_result_sha256=sha(prefix_folder / "result.json"),
                      actual_share_at_threshold_sample=event["actual_upper_sample_share"],
                      prefix_sample_time=event["upper_time"], additional_dwell=case["additional_dwell"],
                      total_cue_on_to_withdrawal=event["upper_time"] + case["additional_dwell"],
                      native_state_interpolated=False)
    else:
        cue = (CenteredAreaCue(case["cue_area"], case["cue_width"], case["cue_center"])
               if case["assay"] == "position" else RelativePeakCue(case["cue_peak"]))
        result = _write_trial(folder, case, state, basal, cue)
    result.update(case=case, status="complete", withdrawal_executed=True)
    result["artifact_sha256"] = {str(p.relative_to(folder)): sha(p) for p in folder.rglob("*")
                                if p.is_file() and p.name != "result.json"}
    write(folder / "result.json", result)
    return result


def prepare(out, calibration_root, width_root):
    out = Path(out).resolve()
    if out.exists():
        raise FileExistsError("Choose a fresh immutable supporting payload directory")
    if not out.is_relative_to((PACKAGE / "data/runs").resolve()):
        raise ValueError("Supporting native payloads belong under manuscript/data/runs")
    plan = build_plan(calibration_root, width_root)
    out.mkdir(parents=True)
    inputs_hash, sources_hash = {}, {}
    for row in plan["prefix_cases"]:
        target = out / "inputs/basals" / (row["basal_id"].replace(":", "_") + ".npz")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(row["checkpoint"], target)
        if sha(target) != row["checkpoint_sha256"]:
            raise ValueError("Basal changed during freeze")
        inputs_hash[str(target.relative_to(out))] = sha(target)
        for case in plan["prefix_cases"] + plan["cases"]:
            if case["basal_id"] == row["basal_id"]:
                case["checkpoint"] = str(target.relative_to(out))
    for source in sorted((PACKAGE / "code").rglob("*")):
        if source.is_file() and "__pycache__" not in source.parts and source.suffix in (".py", ".json", ".txt", ".md"):
            relative = source.relative_to(ROOT)
            target = out / "frozen_sources" / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            before = sha(source)
            shutil.copy2(source, target)
            if sha(source) != before or sha(target) != before:
                raise ValueError("Source changed during freeze")
            sources_hash[str(relative)] = before
    for source in (PACKAGE / "data/inputs/fixed_scales.json",
                   PACKAGE / "data/runs/debelly_corrected_three_test_20261008/manifest.json"):
        relative = source.relative_to(ROOT)
        target = out / "frozen_sources" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        before = sha(source)
        shutil.copy2(source, target)
        if sha(source) != before or sha(target) != before:
            raise ValueError("Runtime scientific input changed during freeze")
        sources_hash[str(relative)] = before
    for name, digest in plan["source_collections"].items():
        source = Path(name)
        target = out / "inputs" / ("calibration_" + source.name if source.name != "manifest.json" else "width_manifest.json")
        shutil.copy2(source, target)
        if sha(target) != digest:
            raise ValueError("Collection changed during freeze")
        inputs_hash[str(target.relative_to(out))] = digest
    plan.update(source_sha256=sources_hash, input_sha256=inputs_hash, result_schema=2,
                created_utc=datetime.now(timezone.utc).isoformat(), execution_status="frozen_for_review_not_submitted")
    write(out / "manifest.json", plan)
    write(out / "freeze.json", dict(manifest_sha256=sha(out / "manifest.json")))
    source = out / "frozen_sources" / Path(__file__).resolve().relative_to(ROOT)
    python = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
    commands = {kind: shlex.join([python, "-B", str(source), "worker", "--out", str(out),
        "--kind", kind, "--task", "SLURM_ARRAY_TASK_ID"]) for kind in ("prefix", "trial")}
    write(out / "proposed_worker_commands.json", dict(commands=commands,
        dependency="All 18 qualified prefix workers must finish before the 162 share branches; root review required before submission.",
        prefix_tasks=18, trial_tasks=270, submitted=False))
    return dict(output=str(out), counts=plan["counts"], submitted=False,
                manifest_sha256=sha(out / "manifest.json"))


def worker(out, kind, task):
    out = Path(out).resolve()
    manifest = read(out / "manifest.json")
    if sha(out / "manifest.json") != read(out / "freeze.json")["manifest_sha256"]:
        raise ValueError("Frozen support manifest changed")
    for name, digest in manifest["source_sha256"].items():
        if sha(out / "frozen_sources" / name) != digest:
            raise ValueError("Frozen supporting source changed: " + name)
    for name, digest in manifest["input_sha256"].items():
        if sha(out / name) != digest:
            raise ValueError("Frozen supporting input changed: " + name)
    rows = manifest["prefix_cases"] if kind == "prefix" else manifest["cases"]
    if not 0 <= int(task) < len(rows):
        raise ValueError("Supporting task outside declared frozen array")
    case = dict(rows[int(task)])
    case["checkpoint"] = str(out / case["checkpoint"])
    folder = out / ("prefixes" if kind == "prefix" else "workers") / f"{int(task):04d}"
    if kind == "prefix":
        result = run_prefix(case, folder)
        result.update(manifest_sha256=sha(out / "manifest.json"), record_kind="shared_native_prefix")
        write(folder / "result.json", result)
        return result
    prefix = out / "prefixes" / f"{case['prefix_task']:04d}" if case["assay"] == "share_withdrawal" else None
    observed = run_case(case, folder, prefix)
    result = dict(generated_utc=datetime.now(timezone.utc).isoformat(), case=case, result_schema=2,
                  manifest_sha256=sha(out / "manifest.json"), support_status=observed["status"],
                  withdrawal_executed=observed["withdrawal_executed"],
                  basal=dict(qualified=True, basal_id=case["basal_id"], checkpoint=case["checkpoint"],
                      checkpoint_sha256=case["checkpoint_sha256"], fixed_scales=case["fixed_scales"]))
    if observed["withdrawal_executed"]:
        trial = {key: value for key, value in observed.items()
                 if key not in ("case", "artifact_sha256", "status")}
        trial.update(trial_directory=str(folder.relative_to(out)), readouts="readouts.json",
                     full_trace_sha256=sha(folder / "trace.npz"))
        trial["artifacts"] = {name: dict(path=str((folder / filename).relative_to(out)), sha256=sha(folder / filename))
            for name, filename in (("trace", "trace.npz"), ("readouts", "readouts.json"), ("native_tail", "native_tail.json"))}
        result["trial"] = trial
    else:
        result.update(prefix_event=observed["prefix_event"],
                      prefix_result_sha256=observed["prefix_result_sha256"])
    result["artifact_sha256"] = {str(path.relative_to(out)): sha(path)
        for path in folder.rglob("*") if path.is_file() and path.name != "result.json"}
    write(folder / "result.json", result)
    return result


def audit_prefixes(out):
    """Rebuild all threshold readouts from the immutable full native histories."""
    from signature_analysis_20261020 import _same_numeric_tree, audit_native_phase_join
    out=Path(out).resolve();manifest=read(out/'manifest.json');digest=sha(out/'manifest.json')
    if digest!=read(out/'freeze.json')['manifest_sha256']:
        raise ValueError('Frozen support manifest changed')
    for name,expected in manifest['source_sha256'].items():
        if sha(out/'frozen_sources'/name)!=expected:raise ValueError('Frozen support source changed')
    for name,expected in manifest['input_sha256'].items():
        if sha(out/name)!=expected:raise ValueError('Frozen support input changed')
    records=[];missing=[];native=_native()
    for original in manifest['prefix_cases']:
        folder=out/'prefixes'/f"{original['task']:04d}";path=folder/'result.json'
        if not path.exists():missing.append(original['task']);continue
        case=dict(original,checkpoint=str(out/original['checkpoint']));result=read(path)
        if result['manifest_sha256']!=digest or result['case']!=case:
            raise ValueError('Shared native prefix identity changed')
        scientific={str(p.relative_to(folder)) for p in folder.rglob('*')
                    if p.is_file() and p.name!='result.json'}
        if scientific!=set(result['artifact_sha256']):raise ValueError('Incomplete shared-prefix artifact inventory')
        for name,expected in result['artifact_sha256'].items():
            if sha(folder/name)!=expected:raise ValueError('Shared native prefix artifact changed')
        _,trace=native.read_checkpoint(folder/'trace.npz',case['model'])
        _,basal=native.read_checkpoint(case['checkpoint'],case['model'])
        measurement=readouts.analyze_trace(_trace_mapping(trace),_parameters(case,basal),_trace_mapping(basal))
        _same_numeric_tree(read(folder/'readouts.json'),measurement,'prefix_readouts')
        events={str(target):share_event(measurement,target) for target in SHARES}
        _same_numeric_tree(result['shares'],events,'share_threshold_events')
        phases=audit_native_phase_join(folder/'trace.npz',result['events'],_trace_mapping(trace))
        records.append(dict(task=original['task'],basal_id=case['basal_id'],passed=True,
            events=events,native_phase_join=phases,result_sha256=sha(path),
            archived_native_restart_checkpoints=True))
    summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=not missing,
        manifest_sha256=digest,analysis_source_sha256=sha(Path(__file__)),
        expected_prefixes=len(manifest['prefix_cases']),audited_prefixes=len(records),missing_tasks=missing,
        native_history_and_readout_reconstruction=True,records=records)
    write(out/'results/prefix_independent_audit.json',summary)
    return summary


def audit_restarts(out):
    """Check withdrawal states against the exact archived shared-prefix samples."""
    out=Path(out).resolve();manifest=read(out/'manifest.json');native=_native()
    records=[];missing=[];cached=None;prefix_trace=None
    for case in manifest['cases']:
        if case['assay']!='share_withdrawal':continue
        folder=out/'workers'/f"{case['task']:04d}";path=folder/'result.json'
        if not path.exists():missing.append(case['task']);continue
        result=read(path);prefix=out/'prefixes'/f"{case['prefix_task']:04d}"
        parent=read(prefix/'result.json');event=parent['shares'][str(case['target_share'])]
        if result['manifest_sha256']!=sha(out/'manifest.json'):raise ValueError('Withdrawal manifest identity changed')
        if not event['observed']:
            if result['withdrawal_executed'] or (folder/'trace.npz').exists():
                raise ValueError('Unobserved share target generated a fabricated withdrawal trace')
            records.append(dict(task=case['task'],passed=True,withdrawal_executed=False,
                                classification='finite_prefix_nonreach_retained'))
            continue
        if not result['withdrawal_executed']:raise ValueError('Observed threshold lost its withdrawal branch')
        if cached!=case['prefix_task']:
            _,prefix_trace=native.read_checkpoint(prefix/'trace.npz',case['model']);cached=case['prefix_task']
        expected=native.state_at(prefix_trace,case['model'],event['upper_sample_index'])
        _,trace=native.read_checkpoint(folder/'trace.npz',case['model'])
        trace.active=trace.active[:1];trace.fields={key:value[:1] for key,value in trace.fields.items()}
        if hasattr(trace,'raw_boundaries'):trace.raw_boundaries=trace.raw_boundaries[:1]
        gate=native.restart_error(expected,trace)
        if not gate['passed']:raise ValueError('Withdrawal branch did not restart the exact shared native sample')
        trial=result['trial']
        if trial['prefix_result_sha256']!=sha(prefix/'result.json') or trial['prefix_sample_time']!=event['upper_time']:
            raise ValueError('Withdrawal branch shared-prefix reference changed')
        records.append(dict(task=case['task'],passed=True,withdrawal_executed=True,
            prefix_task=case['prefix_task'],prefix_sample_time=event['upper_time'],restart_check=gate))
    summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=not missing,
        manifest_sha256=sha(out/'manifest.json'),analysis_source_sha256=sha(Path(__file__)),
        expected_withdrawals=manifest['counts']['share_withdrawal_trials'],audited_withdrawals=len(records),
        missing_tasks=missing,records=records)
    write(out/'results/shared_restart_independent_audit.json',summary)
    return summary


def combine_audits(out):
    """Index valid WP/G and corrected DB audits without rewriting either history."""
    out=Path(out).resolve();manifest=read(out/'manifest.json');records=[];missing=[];counts={}
    prefix=out/'results/prefix_independent_audit.json';restarts=out/'results/shared_restart_independent_audit.json'
    if not read(prefix)['passed'] or not read(restarts)['passed']:
        raise ValueError('Shared native-prefix and exact-restart audits must pass')
    corrected=out/'results/analysis_debelly_corrected_20261020/status.json'
    if not corrected.exists() or not read(corrected)['declared_subset_complete']:
        raise ValueError('Corrected deBelly native/readout audit is incomplete')
    executed=0
    for case in manifest['cases']:
        name='analysis_debelly_corrected_20261020' if case['model']=='debelly' else 'analysis_20261020'
        path=out/'results'/name/'workers'/f"{case['task']:04d}.json"
        if not path.exists():missing.append(case['task']);continue
        analysis=read(path)
        if analysis['task']!=case['task'] or analysis['manifest_sha256']!=sha(out/'manifest.json'):
            raise ValueError('Independent support audit changed its scientific case identity')
        for trial in analysis['trials']:
            if not trial['native_readout_numeric_audit'] or not trial['native_phase_join_audit']['passed']:
                raise ValueError('Native readout/phase audit did not pass')
        result=read(out/'workers'/f"{case['task']:04d}"/'result.json')
        executed+=result['withdrawal_executed'];label=case['assay']
        counts[label]=counts.get(label,0)+1
        records.append(dict(task=case['task'],model=case['model'],assay=label,
            requested_reference_length=case['requested_reference_length'],
            withdrawal_executed=result['withdrawal_executed'],support_status=result['support_status'],
            audited_worker=str(path.relative_to(out)),audited_worker_sha256=sha(path),
            analysis_source_sha256=analysis['analysis_source_sha256']))
    summary=dict(generated_utc=datetime.now(timezone.utc).isoformat(),passed=not missing,
        declared_subset_complete=not missing,manifest_sha256=sha(out/'manifest.json'),
        audit_index_source_sha256=sha(Path(__file__)),expected_workers=len(manifest['cases']),
        audited_workers=len(records),executed_trials=executed,finite_prefix_nonreaches=len(records)-executed,
        missing_tasks=missing,assay_counts=counts,records=records,
        shared_prefix_audit=dict(path=str(prefix.relative_to(out)),sha256=sha(prefix)),
        shared_restart_audit=dict(path=str(restarts.relative_to(out)),sha256=sha(restarts)),
        earlier_checker_history='Original full audit retained; its cue-dependent DB diagnostic continuity errors are superseded only by corrected DB audits.')
    write(out/'results/combined_independent_audit.json',summary);return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--out", type=Path, required=True)
    prep.add_argument("--calibration-root", type=Path, required=True)
    prep.add_argument("--width-root", type=Path, required=True)
    work = sub.add_parser("worker")
    work.add_argument("--out", type=Path, required=True)
    work.add_argument("--kind", choices=("prefix", "trial"), required=True)
    work.add_argument("--task", type=int, required=True)
    audit=sub.add_parser('audit-prefixes')
    audit.add_argument('--out',type=Path,required=True)
    restarts=sub.add_parser('audit-restarts')
    restarts.add_argument('--out',type=Path,required=True)
    combined=sub.add_parser('combine-audits')
    combined.add_argument('--out',type=Path,required=True)
    arguments = parser.parse_args()
    if arguments.command=='prepare':result=prepare(arguments.out,arguments.calibration_root,arguments.width_root)
    elif arguments.command=='audit-prefixes':result=audit_prefixes(arguments.out)
    elif arguments.command=='audit-restarts':result=audit_restarts(arguments.out)
    elif arguments.command=='combine-audits':result=combine_audits(arguments.out)
    else:result=worker(arguments.out,arguments.kind,arguments.task)
    print(json.dumps(result, indent=2, allow_nan=False))
