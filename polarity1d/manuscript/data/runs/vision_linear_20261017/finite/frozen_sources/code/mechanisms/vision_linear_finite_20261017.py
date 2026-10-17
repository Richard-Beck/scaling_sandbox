"""Frozen all-five finite linear exposures; preflight never integrates or submits.

The primary and anchor protocols were authorized before execution. Native model
sources and retained evidence are immutable. Worker execution uses this file in
the campaign's frozen_sources tree, never a changing development import tree.
"""
from pathlib import Path
from datetime import datetime, timezone
from dataclasses import replace
from types import SimpleNamespace
import argparse
import copy
import csv
import gzip
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/core"))
import native as n
import polarity_establishment as establishment
from polarity1d.simulation import simulate

DEFAULT = ROOT / "data/runs/vision_linear_20261017/finite"
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/vision_linear_20261017/supporting"
CANONICAL = ROOT / "data/runs/debelly_corrected_three_test_20261008"
METRICS = ROOT / "data/runs/initial_polarity_length_corrected_20261008/steady_state_metrics.csv"
PILOT = ROOT / "data/runs/all_five_incumbent_recovery_pilot_20261017/manifest.json"
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
HORIZONS = (3000., 8000., 10000., 13000.)
CASE_KEYS = ("case_index", "model", "length", "cells", "width", "width_mode",
             "width_fraction_at_reference", "background_stimulus", "response_amplitude",
             "response_factor", "checkpoint", "checkpoint_sha256", "eligible",
             "initial_orientation", "initial_contrast", "fixed_scales", "condition")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def clean(value):
    if isinstance(value, dict): return {key: clean(v) for key, v in value.items()}
    if isinstance(value, (tuple, list)): return [clean(v) for v in value]
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def csvwrite(path, rows):
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(clean(rows))


def material_weights(nodes, edge):
    """Exact piecewise-linear nodal integral over material interval [0, edge]."""
    q = np.linspace(0., 1., nodes); weights = np.zeros(nodes)
    for i in range(nodes - 1):
        width = q[i+1] - q[i]; covered = min(max(edge-q[i], 0.), width)
        weights[i] += covered-covered**2/(2*width)
        weights[i+1] += covered**2/(2*width)
    np.testing.assert_allclose(weights.sum(), edge, rtol=1e-13, atol=1e-14)
    return weights


def geometry(case, state, x):
    model = case["model"]; length = case["length"]
    u = state["rac"] if model == "debelly" else state[:-2:2] if model == "spring" else state[0]
    physical = float(state[-1]-state[-2]) if model == "spring" else length
    q = np.linspace(0., 1., len(u)) if model == "spring" else x/length
    b0 = float(u.min()); threshold = b0 + .5*float(np.ptp(u))
    below = np.flatnonzero(u < threshold)
    if not len(below) or below[0] == 0: raise ValueError("No boundary-connected incumbent core")
    i = int(below[0]); edge = float(q[i-1]+(threshold-u[i-1])/(u[i]-u[i-1])*(q[i]-q[i-1]))
    weights = material_weights(len(u), edge) if model == "spring" else np.clip(edge-np.arange(len(u))/len(u), 0., 1/len(u))
    initial_raw = float(u @ weights * physical)
    initial_excess = float((u-b0) @ weights * physical)
    return dict(background=b0, threshold=threshold, edge_material=edge,
                core_width_native=edge*physical, physical_length=physical,
                initial_raw=initial_raw, initial_excess=initial_excess,
                initial_peak=float(u.max()), weights=weights,
                cue_core_overlap=bool(1-case["width"]/physical < edge),
                core_mirror_overlap=bool(2*edge > 1.))


def native_model(model):
    result = n.model_for(model)
    # The corrected canonical entry point applies this explicit replacement;
    # native.py alone deliberately retains the historical default constructor.
    return n.ds.DeBelly(replace(result.parameters, center_boundary_motion=False)) if model == "debelly" else result


def source_cases():
    manifest = read(CANONICAL / "manifest.json")
    with gzip.open(CANONICAL / "results/records.json.gz", "rt") as handle:
        records = json.load(handle)["cases"]
    by_case = {c["case_index"]: c for c in records}
    basal = {(r["model"], float(r["domain_length"])): r for r in csv.DictReader(METRICS.open())}
    anchor_sizes = {r["model"]: r["selected_sizes"] for r in read(PILOT)["selection"]}
    jobs = []; ineligible = []
    for original in manifest["cases"]:
        if original["condition"] != "fixed_physical": continue
        case = {key: original.get(key) for key in CASE_KEYS}
        if not case["eligible"]:
            ineligible.append(dict(model=case["model"], size_native=case["length"],
                                   status="preparation_ineligible", canonical_case_index=case["case_index"]))
            continue
        row = basal[case["model"], case["length"]]
        assert case["checkpoint"] == row["checkpoint"] and case["checkpoint_sha256"] == row["checkpoint_sha256"]
        case["basal_physical_length"] = float(row["physical_length"])
        case["basal_core_width"] = float(row["old_pole_width"])
        case["basal_background"] = float(row["background_active"])
        case["basal_threshold"] = float(row["pole_threshold"])
        record = by_case[case["case_index"]]
        zero = next(t for t in record["tests"] if t["duration"] == 0)
        assert zero["no_cue_preparation_valid"]
        anchor = case["length"] in anchor_sizes[case["model"]]
        durations = [50.] if case["model"] == "spring" else [1000.]
        if anchor: durations += [25., 100.] if case["model"] == "spring" else [500., 2000.]
        for duration in durations:
            expected = next((t for t in record["tests"] if t["duration"] == duration), None)
            jobs.append(dict(slot=len(jobs), case=copy.deepcopy(case), duration=duration,
                             primary=duration == (50. if case["model"] == "spring" else 1000.),
                             anchor=anchor, anchor_label=("short", "middle", "long")[anchor_sizes[case["model"]].index(case["length"])] if anchor else "",
                             expected_10000=copy.deepcopy(expected),
                             old_zero_cue=dict(maximum_scaled_excursion=zero["no_cue_maximum_scaled_excursion"],
                                               passed=zero["no_cue_preparation_valid"],
                                               trace=zero["trace"], trace_prefix=zero["trace_prefix"])))
    assert len(jobs) == 127 and sum(j["primary"] for j in jobs) == 97 and len(ineligible) == 13
    return jobs, ineligible, anchor_sizes


def sample_times(model, end, phase):
    early = .5 if model == "debelly" else .02 if model == "spring" else .1
    pieces = [np.arange(0., min(10., end)+early/2, early), np.arange(10., min(100., end)+1e-8, 2.5),
              np.arange(100., end+1e-8, 10.), np.array([0., end])]
    if phase == "cue":
        pieces += [end*np.array([.01, .1, .25, .5, .75, 1.]), np.arange(max(0., end-100.), end+1e-8, 2.5)]
    else:
        pieces += [np.array([.5, 2.5, 10., 50., 100., 500., 1000., *HORIZONS])]
        pieces += [np.arange(h-500., h+1e-8, 2.5) for h in HORIZONS]
    times = np.unique(np.round(np.concatenate(pieces), 10))
    times = times[(times >= 0.) & (times <= end)]
    assert times[0] == 0. and times[-1] == end and np.all(np.diff(times) > 0)
    if model == "debelly":
        np.testing.assert_allclose(times/.5, np.round(times/.5), atol=1e-11, rtol=0)
    return times


def preflight(destination=None):
    jobs, ineligible, anchors = source_cases(); seen = set(); checks = []
    assert read(CANONICAL / "results/independent_audit.json")["passed"]
    assert native_model("debelly").parameters.center_boundary_motion is False
    for job in jobs:
        c = job["case"]; path = ROOT / c["checkpoint"]
        if c["checkpoint"] not in seen:
            assert sha(path) == c["checkpoint_sha256"]; seen.add(c["checkpoint"])
        state, tr = n.read_checkpoint(path, c["model"]); g = geometry(c, state, tr.x)
        for key, expected in (("background", c["basal_background"]), ("threshold", c["basal_threshold"]),
                              ("core_width_native", c["basal_core_width"]), ("physical_length", c["basal_physical_length"])):
            np.testing.assert_allclose(g[key], expected, rtol=2e-11, atol=2e-12)
        cue = establishment.Patch(c["response_amplitude"], c["width"], side="right", width_mode="physical", background=0.)
        L = g["physical_length"] if c["model"] == "spring" else c["length"]
        assert cue(0., np.array([L]), L)[0] == c["response_amplitude"]
        assert cue(0., np.array([max(0., L-c["width"])]), L)[0] <= 1e-14
        for phase, end in (("cue", job["duration"]), ("post", 13000.)):
            times = sample_times(c["model"], end, phase)
            assert len(times) > 2
        checks.append(dict(slot=job["slot"], model=c["model"], size_native=c["length"],
                           cue_duration=job["duration"], primary=job["primary"],
                           cue_core_overlap=g["cue_core_overlap"], native_checkpoint_hash_verified=True))
    result = dict(generated_at=now(), passed=True, integrations=0, submitted_jobs=0,
                  eligible_basals=len(seen), primary_branches=97, sensitivity_branches=30,
                  ineligible_preparations=len(ineligible), anchor_sizes=anchors,
                  corrected_debelly=True, all_debelly_output_times_native_half_units=True, cases=checks,
                  script_sha256=sha(Path(__file__)))
    if destination: write(destination, result)
    return result


def initialize(run):
    if run.exists(): raise ValueError("Choose a fresh finite campaign directory")
    check = preflight(); jobs, ineligible, anchors = source_cases()
    for folder in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/folder).mkdir(parents=True)
    sources = [Path(__file__)] + sorted((ROOT/"code/core").glob("*.py")) + sorted((ROOT/"code/core").glob("*.json"))
    sources += sorted((ROOT/"code/core/models").glob("*.py")) + sorted((ROOT/"code/core/src/polarity1d").glob("*.py"))
    sources += [METRICS, PILOT, CANONICAL/"manifest.json", CANONICAL/"results/records.json.gz",
                CANONICAL/"results/independent_audit.json", TOPIC/"metadata/FINITE_DESIGN.md", TOPIC/"metadata/LONG_CUE_INVENTORY.md"]
    hashes = {}
    for source in sources:
        relative = source.relative_to(ROOT); target = run/"frozen_sources"/relative
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target); hashes[str(relative)] = sha(target)
    inputs = {}
    for job in jobs:
        c = job["case"]; key = c["checkpoint"]
        if key not in inputs:
            target = run/"inputs"/f'basal_{len(inputs):03d}.npz'; shutil.copy2(ROOT/key, target)
            assert sha(target) == c["checkpoint_sha256"]; inputs[key] = str(target.relative_to(run))
        job["input"] = inputs[key]
    manifest = dict(generated_at=now(), experiment="All-five established linear finite local-cue competition",
                    geometry="linear_fixed_support", primary_duration_native=1000., spring_primary_duration_native=50.,
                    withdrawal_horizon_native=13000., readout_horizons=list(HORIZONS),
                    cue="Right raised cosine; support3(Spring0.2), peaks0.1/6.4/3.6/1/4, zero extra-input background",
                    jobs=jobs, ineligible=ineligible, anchor_sizes=anchors, frozen_source_sha256=hashes,
                    native_numerics=dict(RD_backend="BDF", RD_rtol=2e-6, RD_atol=1e-9, RD_max_step=10.,
                                         Spring_backend="LSODA", Spring_rtol=2e-7, Spring_atol=2e-9, Spring_max_step=2.5,
                                         deBelly_step=.5, deBelly_center_boundary_motion=False),
                    old_expectations="Original10000-withdrawal scalar samples and separate field-tail flags; no new shams",
                    preflight=check)
    write(run/"manifest.json", manifest)
    write(run/"provenance.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    write(run/"results/status.json", dict(generated_at=now(), producer="initialized", complete=False, jobs=127))
    return manifest


def verify(run):
    manifest = read(run/"manifest.json")
    assert sha(run/"manifest.json") == read(run/"provenance.json")["manifest_sha256"]
    for relative, digest in manifest["frozen_source_sha256"].items():
        assert sha(run/"frozen_sources"/relative) == digest, relative
    for job in manifest["jobs"]:
        assert sha(run/job["input"]) == job["case"]["checkpoint_sha256"]
    return manifest


def integrate(case, state, cue, times):
    """Only output times differ from canonical runs; native dynamics are unchanged."""
    model = case["model"]; L = case["length"]; cells = case["cells"]; initial_roundoff = 0.
    if model == "debelly":
        native = native_model(model); assert native.parameters.center_boundary_motion is False
        tr = n.ds.simulate_debelly(native, protocol=cue, length=L, cells=cells,
                                  end_time=times[-1], times=times, step=.5, initial_checkpoint=copy.deepcopy(state))
        fields = {key: value.copy() for key, value in tr.fields.items()}
        final = {key: fields[key][-1].copy() for key in n.DB_NAMES}
    elif model == "spring":
        q = np.linspace(0., 1., cells)
        def rhs(t, y):
            physical = y[-1]-y[-2]; profile = np.asarray(cue(t, q*physical, physical))
            result = n.sp.rhs(t, y, L0=L, d1=80., b=4., cue=0.)
            if np.any(profile):
                u = y[:-2:2]; v = y[1:-2:2]; delta = 3.+80.*(physical-L)
                transfer = n.sp.auth["f"](u, v, 4.+profile, 5, 6, 2, delta)-n.sp.auth["f"](u, v, 4., 5, 6, 2, delta)
                result[:-2:2] += transfer; result[1:-2:2] -= transfer
            return result
        sol = solve_ivp(rhs, (0., times[-1]), state, method="LSODA", t_eval=times,
                        rtol=2e-7, atol=2e-9, max_step=2.5)
        if not sol.success: raise RuntimeError(sol.message)
        y = sol.y.T.copy()
        # LSODA's t_eval interpolation can return t=0 one ulp away from y0.
        # The exact passed packed state is the authoritative initial sample;
        # bound/report interpolation roundoff before recording that known state.
        initial_roundoff = float(np.max(abs(y[0]-state)))
        allowance = 8*np.finfo(float).eps*max(1., float(np.max(abs(state))))
        if initial_roundoff > allowance: raise ArithmeticError("Spring initial output exceeds eight-ulp scale allowance")
        y[0] = state
        fields = dict(active=y[:, :-2:2], inactive=y[:, 1:-2:2], length=y[:, -1]-y[:, -2],
                      left_boundary=y[:, -2], right_boundary=y[:, -1])
        final = y[-1].copy()
    else:
        native = n.model_for(model)
        if np.min(state) < 0: raise ValueError("Negative restart; no replacement/clipping permitted")
        tr = simulate(native, initial_state=state, length=L, cells=cells, protocol=cue,
                      end_time=times[-1], times=times, backend="bdf", max_step=10.)
        fields = {key: value.copy() for key, value in tr.fields.items()}
        allowance = 8*np.finfo(float).eps*max(1., float(np.max(abs(state))))
        for i, key in enumerate(native.state_names):
            error = float(np.max(abs(fields[key][0]-state[i])))
            initial_roundoff = max(initial_roundoff, error)
            if error > allowance: raise ArithmeticError("RD initial output exceeds eight-ulp scale allowance")
            fields[key][0] = state[i]
        final = np.vstack([fields[key][-1] for key in native.state_names])
    if any(not np.isfinite(y).all() for y in fields.values()): raise ValueError("Nonfinite native field")
    return dict(time=np.array(times), fields=fields, state=final,
                exact_input_recorded_at_zero=True, initial_output_roundoff_max=initial_roundoff)


def primary_field(model):
    return "rac" if model == "debelly" else "rac_active" if model == "holmes_model3" else "active"


def coupled_fields(case):
    if case["model"] == "debelly": return list(n.DB_NAMES)
    if case["model"] == "spring": return ["active", "inactive", "length", "left_boundary", "right_boundary"]
    return list(native_model(case["model"]).state_names)


def native_orientation(case, fields, x):
    # Retain each model's actual canonical scalar convention; regional
    # integrals always use the unmodified raw activity stored in fields.
    raw = fields[primary_field(case["model"])]; u = np.maximum(raw, 0.)
    if case["model"] == "spring":
        q = np.linspace(0., 1., case["cells"])
        return np.trapezoid(raw*(2*q-1), q, axis=1)/np.trapezoid(raw, q, axis=1)
    if case["model"] == "debelly": return n.ds._orientation(raw)
    q = x/(x[0]+x[-1])
    return (u @ (2*q-1))/u.sum(axis=1)


def monitored(case, fields):
    model = case["model"]
    if model == "debelly":
        result = {key: fields[key] for key in ("rac", "rho", "mca_bound", "mca_unbound", "velocity", "tension")}
        result["displacement"] = fields["displacement"]-fields["displacement"].mean(axis=1)[:, None]
        result["physical_length"] = case["length"]+fields["boundary_separation"]
        return result
    return {key: fields[key] for key in case["fixed_scales"]}


def gates(case, phase, x, horizon):
    times = phase["time"]; keep = (times >= horizon-500.-1e-9) & (times <= horizon+1e-9)
    t = times[keep]; np.testing.assert_allclose(t, np.arange(horizon-500., horizon+1e-8, 2.5), atol=1e-9)
    ori = native_orientation(case, phase["fields"], x)[keep]
    contrast = np.ptp(np.maximum(phase["fields"][primary_field(case["model"])] , 0.), axis=1)[keep]
    orientation = bool(ori.min() > .05 and contrast.min() >= .05 and np.ptp(ori) < .005)
    error = float(np.max(abs(ori/abs(case["initial_orientation"])-1)))
    maximum = 0.; field_metrics = {}
    for key, field in monitored(case, phase["fields"]).items():
        y = field[keep]; scale = case["fixed_scales"][key]
        variation = float(np.ptp(y, axis=0).max()/scale)
        drift = float(abs(y[-1]-y[-2]).max()/2.5*500/scale)
        maximum = max(maximum, variation, drift); field_metrics[key] = dict(range_scaled=variation, projected_drift_500=drift)
    return dict(orientation_success=orientation, full_polarity_success=bool(orientation and error <= .01),
                final_orientation=float(ori[-1]), minimum_contrast_last_500=float(contrast.min()),
                orientation_range_last_500=float(np.ptp(ori)), maximum_relative_polarity_error_last_500=error,
                field_tail_pass=bool(maximum <= .001), field_tail_max_scaled=maximum,
                field_tail_metrics=field_metrics)


def regional(case, phase, g, x):
    u = phase["fields"][primary_field(case["model"])]; q = np.linspace(0., 1., case["cells"]) if case["model"] == "spring" else x/case["length"]
    physical = phase["fields"]["length"] if case["model"] == "spring" else np.full(len(u), case["length"])
    old_raw = u @ g["weights"] * physical; new_raw = u @ g["weights"][::-1] * physical
    old_excess = (u-g["background"]) @ g["weights"] * physical
    new_excess = (u-g["background"]) @ g["weights"][::-1] * physical
    old_peak = np.max(u[:, q <= g["edge_material"]], axis=1)
    new_peak = np.max(u[:, q >= 1-g["edge_material"]], axis=1)
    high = u >= g["threshold"]
    sites = high[:, 0].astype(int)+np.sum(high[:, 1:] & ~high[:, :-1], axis=1)
    return dict(physical_length=physical, old_raw=old_raw, challenger_raw=new_raw,
                old_excess=old_excess, challenger_excess=new_excess,
                old_raw_ratio=old_raw/g["initial_raw"], challenger_raw_ratio=new_raw/g["initial_raw"],
                old_excess_ratio=old_excess/g["initial_excess"], challenger_excess_ratio=new_excess/g["initial_excess"],
                old_peak=old_peak, challenger_peak=new_peak,
                old_peak_threshold_ratio=(old_peak-g["background"])/(g["threshold"]-g["background"]),
                challenger_peak_threshold_ratio=(new_peak-g["background"])/(g["threshold"]-g["background"]),
                sites_above_initial_halfcontrast=sites,
                native_orientation=native_orientation(case, phase["fields"], x),
                native_contrast=np.ptp(np.maximum(u, 0.), axis=1),
                minimum_raw_active=u.min(axis=1))


def first_crossing(times, values, level, above):
    good = values >= level if above else values <= level
    indices = np.flatnonzero(good)
    if not len(indices): return dict(time=None, left=None, right=None)
    i = int(indices[0]); value = float(times[i])
    if i and values[i] != values[i-1]:
        value = float(times[i-1]+(level-values[i-1])/(values[i]-values[i-1])*(times[i]-times[i-1]))
    return dict(time=value, left=float(times[max(0, i-1)]), right=float(times[i]))


def worker(run, slot):
    started = time.perf_counter(); manifest = verify(run); job = manifest["jobs"][slot]; c = job["case"]
    state, checkpoint = n.read_checkpoint(run/job["input"], c["model"])
    g = geometry(c, state, checkpoint.x)
    cue = establishment.Patch(c["response_amplitude"], c["width"], side="right", width_mode="physical", background=0.)
    pulse = integrate(c, state, cue, sample_times(c["model"], job["duration"], "cue"))
    post = integrate(c, pulse["state"], establishment.withdrawal(c), sample_times(c["model"], 13000., "post"))
    initial_fields = {key: checkpoint.fields[key][-1] for key in coupled_fields(c) if key in checkpoint.fields}
    for key, initial in initial_fields.items(): np.testing.assert_array_equal(pulse["fields"][key][0], initial)
    # Cue-dependent derived activation rates may jump at withdrawal; the
    # complete dynamical state and mechanical controllers remain continuous.
    for key in coupled_fields(c): np.testing.assert_array_equal(post["fields"][key][0], pulse["fields"][key][-1])
    if c["model"] == "debelly":
        for key in n.DB_NAMES: np.testing.assert_array_equal(pulse["fields"][key][0], state[key])
    measurements = {str(int(h)): gates(c, post, checkpoint.x, h) for h in HORIZONS}
    expected = job["expected_10000"]; comparison = None
    if expected:
        observed = measurements["10000"]
        comparisons = {key: observed[key] == expected[key] for key in ("orientation_success", "full_polarity_success", "field_tail_pass")}
        comparison = dict(flags=comparisons, all_flags_equal=all(comparisons.values()),
                          final_orientation_absolute_error=abs(observed["final_orientation"]-expected["final_orientation"]),
                          field_tail_absolute_error=abs(observed["field_tail_max_scaled"]-expected["field_tail_max_scaled"]))
        if not all(comparisons.values()): raise RuntimeError("Archived10000 scalar/field flag differs: "+str(comparison))
    arrays = dict(x=checkpoint.x, material_coordinate=np.linspace(0., 1., c["cells"]) if c["model"] == "spring" else checkpoint.x/c["length"],
                  core_weights_material=g["weights"])
    regions = {}
    for phase, tr in (("cue", pulse), ("post", post)):
        regions[phase] = regional(c, tr, g, checkpoint.x)
        arrays[phase+"_time"] = tr["time"]
        arrays.update({phase+"_field_"+key: value for key, value in tr["fields"].items()})
        arrays.update({phase+"_regional_"+key: value for key, value in regions[phase].items()})
    if isinstance(state, dict): arrays.update({"originalstate_"+key: value for key, value in state.items()})
    else: arrays["originalstate"] = state
    clocks = {}
    full_t = np.r_[pulse["time"], job["duration"]+post["time"][1:]]
    for name, key, level, above in (("first_incumbent_excess_half_loss", "old_excess_ratio", .5, False),
                                   ("first_challenger_peak_initial_halfcontrast", "challenger_peak_threshold_ratio", 1., True),
                                   ("first_challenger_excess_half_initial", "challenger_excess_ratio", .5, True)):
        values = np.r_[regions["cue"][key], regions["post"][key][1:]]
        clocks[name] = first_crossing(full_t, values, level, above)
        clocks[name]["phase"] = "not_reached" if clocks[name]["time"] is None else "cue_on" if clocks[name]["time"] <= job["duration"] else "postwithdrawal"
    destination = run/"workers"/f'{slot:03d}.npz'
    temporary = destination.with_suffix(".tmp.npz"); np.savez_compressed(temporary, **arrays); temporary.replace(destination)
    write(destination.with_suffix(".json"), dict(generated_at=now(), job=job, geometry={key:value for key,value in g.items() if key != "weights"},
          gates=measurements, clocks=clocks, canonical_10000_comparison=comparison,
          exact_checkpoint_start=True, exact_withdrawal_continuity=True,
          initial_output_roundoff_max={"cue": pulse["initial_output_roundoff_max"], "post": post["initial_output_roundoff_max"]},
          archive_sha256=sha(destination), worker_seconds=time.perf_counter()-started,
          manifest_sha256=sha(run/"manifest.json"), execution_script_sha256=sha(Path(__file__))))
    print(json.dumps(dict(slot=slot, model=c["model"], size=c["length"], duration=job["duration"], complete=True)), flush=True)


def submit(run):
    manifest = verify(run)
    if (run/"launch.json").exists(): raise ValueError("Campaign already submitted")
    runner = run/"frozen_sources/code/mechanisms"/Path(__file__).name
    command = shlex.join([sys.executable, str(runner), "worker", "--out", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G",
          "--time=00:30:00", f'--array=0-{len(manifest["jobs"])-1}', "--job-name=vision-linear-finite",
          "--chdir="+str(ROOT), "--output="+str(run/"logs/%A_%a.out"), "--error="+str(run/"logs/%A_%a.err"),
          "--wrap", command], text=True).strip()
    write(run/"launch.json", dict(generated_at=now(), array=job, jobs=len(manifest["jobs"])))
    write(run/"results/status.json", dict(generated_at=now(), producer="running", complete=False, array=job, jobs=len(manifest["jobs"])))
    print(job, flush=True)


def identity(job):
    c = job["case"]
    return dict(slot=job["slot"], model=c["model"], size_native=c["length"], physical_length=c["basal_physical_length"],
                cue_duration_native=job["duration"], cue_peak=c["response_amplitude"], cue_support_native=c["width"],
                geometry="linear_fixed_support", fixed_support=True, primary=job["primary"], anchorflag=job["anchor"],
                anchor_label=job["anchor_label"], checkpoint=c["checkpoint"], checkpoint_sha256=c["checkpoint_sha256"])


def finish(run, topic):
    manifest = verify(run); summaries = []; samples = []; audits = []; allocated = 0
    for job in manifest["jobs"]:
        p = run/"workers"/f'{job["slot"]:03d}.npz'; meta = read(p.with_suffix(".json"))
        assert meta["job"] == job and sha(p) == meta["archive_sha256"]
        assert meta["manifest_sha256"] == sha(run/"manifest.json")
        assert meta["execution_script_sha256"] == manifest["frozen_source_sha256"]["code/mechanisms/"+Path(__file__).name]
        c = job["case"]; state, initial = n.read_checkpoint(run/job["input"], c["model"]); g = geometry(c, state, initial.x)
        summary = identity(job); summary.update(status="complete", cue_core_overlap=g["cue_core_overlap"],
                 core_mirror_overlap=g["core_mirror_overlap"], basal_background=g["background"], initial_halfcontrast_threshold=g["threshold"],
                 core_width_native=g["core_width_native"], initial_incumbent_raw=g["initial_raw"],
                 initial_incumbent_excess=g["initial_excess"], archive=str(p.relative_to(run)), archive_sha256=meta["archive_sha256"],
                 manifest_sha256=sha(run/"manifest.json"), source_sha256=meta["execution_script_sha256"],
                 old_zero_cue_maximum_scaled_excursion=job["old_zero_cue"]["maximum_scaled_excursion"])
        with np.load(p) as archive:
            rebuilt = {}
            for phase in ("cue", "post"):
                prefix = phase+"_field_"
                fields = {key[len(prefix):]: archive[key] for key in archive.files if key.startswith(prefix)}
                tr = dict(time=archive[phase+"_time"], fields=fields)
                values = regional(c, tr, g, initial.x); rebuilt[phase] = (tr, values)
                for key, value in values.items(): np.testing.assert_allclose(value, archive[phase+"_regional_"+key], rtol=1e-12, atol=1e-12)
                wanted = job["duration"]*np.array([0., .01, .1, .25, .5, .75, 1.]) if phase == "cue" else np.array([0., .5, 2.5, 10., 50., 100., 500., 1000., *HORIZONS])
                for target in wanted:
                    i = int(np.argmin(abs(tr["time"]-target))); assert abs(tr["time"][i]-target) < 1e-8
                    sample = dict(summary, phase=phase, phase_time_native=float(tr["time"][i]),
                                  absolute_time_native=float(tr["time"][i]+(job["duration"] if phase == "post" else 0.)))
                    sample.update({key: float(value[i]) for key, value in values.items()})
                    samples.append(sample)
                if phase == "cue":
                    for key, value in values.items(): summary["withdrawal_"+key] = float(value[-1])
            for key in coupled_fields(c):
                np.testing.assert_array_equal(rebuilt["cue"][0]["fields"][key][-1], rebuilt["post"][0]["fields"][key][0])
            for h in HORIZONS:
                result = gates(c, rebuilt["post"][0], initial.x, h)
                for key in result:
                    if key != "field_tail_metrics": summary[f'post{int(h)}_'+key] = result[key]
                i = np.flatnonzero(abs(rebuilt["post"][0]["time"]-h) < 1e-8)[0]
                summary.update({f'post{int(h)}_'+key: float(value[i]) for key,value in rebuilt["post"][1].items()})
                assert clean(result) == meta["gates"][str(int(h))]
        for name, clock in meta["clocks"].items():
            summary[name+"_native"] = clock["time"]; summary[name+"_phase"] = clock["phase"]
            summary[name+"_bracket_left_native"] = clock["left"]; summary[name+"_bracket_right_native"] = clock["right"]
        if job["expected_10000"]:
            assert meta["canonical_10000_comparison"]["all_flags_equal"]
        summaries.append(summary); allocated += p.stat().st_blocks*512
        audits.append(dict(slot=job["slot"], checkpoint_and_withdrawal=True, all_regional_values_reconstructed=True,
                           all_horizon_gates_reconstructed=True, canonical_10000_comparison=meta["canonical_10000_comparison"]))
    for row in manifest["ineligible"]:
        summaries.append(dict(row, primary=True, cue_duration_native=50. if row["model"] == "spring" else 1000.,
                              cue_support_native=.2 if row["model"] == "spring" else 3., geometry="linear_fixed_support", fixed_support=True))
    csvwrite(run/"results/metrics.csv", summaries); csvwrite(run/"results/samples.csv", samples)
    write(run/"results/verification.json", dict(generated_at=now(), passed=True, cases=audits,
          complete_primary=97, complete_sensitivity=30, ineligible_visible=13, allocated_worker_bytes=allocated))
    write(run/"results/status.json", dict(generated_at=now(), producer="finished", complete=True, jobs=127,
          verification_passed=True, allocated_worker_bytes=allocated))
    render(run, topic, summaries, samples)


def render(run, topic, rows, samples):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    topic.mkdir(parents=True, exist_ok=True)
    vectors = topic/"supporting/vectors"; vectors.mkdir(parents=True, exist_ok=True)
    report = topic/"supporting/metadata/FINITE_RESULTS.md"; report.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 5, figsize=(15, 6.2), constrained_layout=True)
    labels = ("Wave-pinning", "Goryachev", "deBelly", "Holmes 3", "Spring")
    for i, model in enumerate(MODELS):
        group = sorted([r for r in rows if r.get("model") == model and r.get("primary") and r["status"] == "complete"], key=lambda r:r["size_native"])
        x = [r["size_native"] for r in group]
        for name, color, label in (("old_excess_ratio", "#325c9b", "Incumbent"), ("challenger_excess_ratio", "#b74b32", "Challenger")):
            axes[0, i].plot(x, [r["withdrawal_"+name] for r in group], color=color, marker="o", ms=3, label=label)
            axes[1, i].plot(x, [r["post13000_"+name] for r in group], color=color, marker="o", ms=3, label=label)
        axes[0, i].set_title(labels[i]); axes[1, i].set_xlabel("Native resting length" if model == "spring" else "Interval length")
        for ax in axes[:, i]: ax.axhline(0, lw=.6, color="gray"); ax.grid(alpha=.15)
    axes[0, 0].set_ylabel("At withdrawal\nCore excess / initial incumbent")
    axes[1, 0].set_ylabel("13,000 cue-free units\nCore excess / initial incumbent")
    axes[0, 0].legend(frameon=False)
    fig.savefig(topic/"vision_linear_finite_competition.png", dpi=180, facecolor="white")
    for ext in ("pdf", "svg"): fig.savefig(vectors/f'vision_linear_finite_competition.{ext}', facecolor="white")
    plt.close(fig)
    counts = []
    for model in MODELS:
        group = [r for r in rows if r.get("model") == model and r.get("primary") and r["status"] == "complete"]
        counts.append((model, len(group), sum(r["post10000_orientation_success"] for r in group),
                       sum(r["post10000_full_polarity_success"] for r in group),
                       sum(r["post13000_orientation_success"] for r in group),
                       sum(r["post13000_full_polarity_success"] for r in group)))
    table = "\n".join(f'|{m}|{count}|{o}|{s}|{o13}|{s13}|' for m,count,o,s,o13,s13 in counts)
    text = f'''# Fixed local-cue competition in all five established linear models

Generated {now()}; producer finished. These development results keep interval and ring evidence separate. No movement law is assumed.

![All-five finite competition](../../vision_linear_finite_competition.png)

The same model-specific peak and physical support are applied across every retained size. The primary pulse lasts1000 native units (Spring50), followed by13,000 cue-free units. The top row shows activity at withdrawal; the bottom row shows what remains after withdrawal. Signed excess in fixed basal cores is normalized by initial incumbent excess; raw integrals, peaks, complete native fields and all unavailable crossing clocks are retained separately. Negative excess is displayed without clipping. These local-region amounts are not unique-site or migration classifications.

|Model|Eligible|10,000 orientation|10,000 strength|13,000 orientation|13,000 strength|
|---|---:|---:|---:|---:|---:|
{table}

The independent old10,000-unit scalar outcomes supply fixed expectations for all97 primary branches. The new histories reproduce those three separate original flags, including the field-tail diagnostic. A scalar strength pass does not establish whole-field convergence. The long finite exposure is matched by input geometry, not by producing a desired fate; Holmes' reference strength clock exceeds1000 and larger-size failures remain visible. Full nonmonotone size responses belong in the interpretation.

Experimental counterparts are local reporter amounts, peak intensities, spatial extent, challenger onset and incumbent loss/rebound, measured during a fixed local stimulus and after its removal. These measurements motivate polarity mechanisms; no cell velocity or migration persistence is inferred from the first moment.

All97 qualifying basals and13 unavailable preparations are retained in [metrics.csv]({os.path.relpath(run/'results/metrics.csv', report.parent)}). [Declared samples]({os.path.relpath(run/'results/samples.csv', report.parent)}) include cue-on fractions, exact withdrawal and postcue3000/8000/10,000/13,000. The30 short/middle/long anchor sensitivities change cue duration only. [Verification]({os.path.relpath(run/'results/verification.json', report.parent)}) reconstructs every regional value, gate and native phase boundary from archived fields. [Protocol and frozen sources]({os.path.relpath(run/'manifest.json', report.parent)}) identify the complete execution closure. Fixed-support/core overlap is flagged, particularly GoryachevL3. Spring core membership uses its initial material coordinate while integrals use its own evolving physical length. deBelly uses the corrected independent endpoints; spatial integrals remain in nominal native coordinates.
'''
    report.write_text(text)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("preflight", "init", "submit", "worker", "finish"))
    parser.add_argument("--out", type=Path, default=DEFAULT); parser.add_argument("--slot", type=int)
    parser.add_argument("--output", type=Path); parser.add_argument("--topic", type=Path, default=TOPIC.parent)
    args = parser.parse_args(); run = args.out.resolve()
    if args.action == "preflight": print(json.dumps(clean(preflight(args.output)), indent=2))
    elif args.action == "init": initialize(run)
    elif args.action == "submit": submit(run)
    elif args.action == "worker": worker(run, args.slot)
    else: finish(run, args.topic.resolve())


if __name__ == "__main__":
    main()
