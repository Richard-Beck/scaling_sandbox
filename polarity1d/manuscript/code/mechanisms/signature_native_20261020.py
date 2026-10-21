"""Native integration/checkpoint adapters for the 20 October signature campaign.

No campaign is submitted by this module. Integration calls belong in SLURM
workers. Native kinetics and the corrected independent deBelly endpoints are
unchanged. Spatial readouts and commanded cues live in separate modules.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import copy
import hashlib
import json
import shutil
import sys

import numpy as np

PACKAGE = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PACKAGE / "code/core"))
import native as n
import polarity_establishment as establishment
import settling

OFF_HORIZONS = {"wave_pinning": 10000., "goryachev": 10000.,
                "holmes_model3": 10000., "spring": 10000., "debelly": 4000.}
DB_NAMES = tuple(n.DB_NAMES)


def model_for(model):
    native_model = n.model_for(model)
    if model == "debelly":
        native_model = n.ds.DeBelly(replace(native_model.parameters, center_boundary_motion=False))
        assert native_model.parameters.center_boundary_motion is False
    return native_model


def active_key(model):
    return "rac" if model == "debelly" else "rac_active" if model == "holmes_model3" else "active"


def state_at(trace, model=None, index=-1):
    model = model or trace.model
    if model == "debelly":
        return {name: trace.fields[name][index].copy() for name in DB_NAMES}
    if model == "spring":
        state = np.empty(2 * trace.active.shape[1] + 2)
        state[:-2:2] = trace.fields["active"][index]
        state[1:-2:2] = trace.fields["inactive"][index]
        state[-2:] = trace.raw_boundaries[index]
        return state
    return np.vstack([trace.fields[name][index] for name in model_for(model).state_names])


def sample_times(model, duration, interval=2.5, *, early=True, debelly_step=.5):
    """Explicit output clock; deBelly outputs never subdivide its native steps."""
    duration, interval = float(duration), float(interval)
    if duration <= 0 or interval <= 0:
        raise ValueError("positive duration and output interval required")
    if model == "debelly":
        lattice = np.array([duration, interval]) / debelly_step
        if not np.allclose(lattice, np.rint(lattice), rtol=0, atol=1e-9):
            raise ValueError("deBelly duration/output interval must lie on the native update lattice")
        early_values = np.arange(0., min(100., duration) + debelly_step / 2, debelly_step)
    elif model == "spring":
        early_values = np.r_[np.arange(0., min(.1, duration), .001),
                             np.arange(0., min(1., duration), .01),
                             np.arange(0., min(10., duration), .1)]
    else:
        early_values = np.r_[np.arange(0., min(1., duration), .001),
                             np.arange(0., min(10., duration), .01),
                             np.arange(0., min(100., duration), .1)]
    values = np.r_[np.arange(0., duration, interval), early_values if early else []]
    values = np.unique(np.round(values, 12))
    values = values[(values > 0) & (values < duration) & (np.abs(values - duration) > 1e-12)]
    # Keep the requested endpoint itself: rounded display/output grids must not
    # alter a bisection duration or a pulse/gap boundary by even a grid rounding.
    return np.r_[0., values, duration]


def _finish(trace, model, length):
    trace.model, trace.length = model, float(length)
    trace.cells = trace.active.shape[1]
    trace.contrast = np.ptp(trace.active, axis=1)
    base = np.min(trace.active, axis=1)
    split = max(2, trace.cells // 3)
    trace.left = np.max(trace.active[:, :split], axis=1) - base
    trace.right = np.max(trace.active[:, -split:], axis=1) - base
    trace.physical_length = (trace.fields["length"] if model == "spring" else
                             np.full(len(trace.time), float(length)))
    trace.sample_kind = "nodes" if model == "spring" else "cells"
    trace.state = state_at(trace)
    if not all(np.isfinite(value).all() for value in trace.fields.values()):
        raise FloatingPointError("nonfinite full native field")
    return trace


def run_native(model, length, state, cue, times, cells=201, max_step=None, debelly_step=.5):
    """Run a constant-stage native adapter, retaining the full coupled state.

    Exact pulse discontinuities are implemented by ``run_phases`` restarts.
    Spring's requested output times change t_eval alone; its native LSODA
    tolerances and max_step=2.5 remain unchanged.
    """
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError("times must start at zero and increase strictly")
    if model == "debelly":
        if not np.allclose(times / debelly_step, np.rint(times / debelly_step), rtol=0, atol=1e-9):
            raise ValueError("deBelly outputs must not change the native integration-step lattice")
        native = n.ds.simulate_debelly(model_for(model), protocol=cue, length=float(length),
                    cells=int(cells), end_time=float(times[-1]), times=times,
                    step=float(debelly_step), initial_checkpoint=state)
        fields = {name: np.asarray(value) for name, value in native.fields.items()}
        trace = SimpleNamespace(time=native.time, x=native.x, fields=fields,
                                active=fields["rac"], orientation=n.ds._orientation(fields["rac"]))
    elif model == "spring":
        import scipy.integrate
        original = scipy.integrate.solve_ivp
        captures = []

        def capture(fun, span, y0, *args, **kwargs):
            assert kwargs["method"] == "LSODA" and kwargs["rtol"] == 2e-7 and kwargs["atol"] == 2e-9
            assert kwargs["max_step"] == 2.5
            if max_step is not None:
                kwargs["max_step"] = float(max_step)
            selected = times[(times >= span[0]) & (times <= span[1])]
            kwargs["t_eval"] = np.unique(np.r_[span[0], selected, span[1]])
            result = original(fun, span, y0, *args, **kwargs)
            captures.append(result)
            return result

        scipy.integrate.solve_ivp = capture
        try:
            native = n.sp.run_spring(length=float(length), stimulus=cue, end_time=float(times[-1]),
                          output_interval=2.5, cells=int(cells), initial_state=state)
        finally:
            scipy.integrate.solve_ivp = original
        raw = np.concatenate([item.y.T if index == 0 else item.y.T[1:]
                              for index, item in enumerate(captures)], axis=0)
        np.testing.assert_allclose(native["time"], times, rtol=0, atol=1e-12)
        trace = SimpleNamespace(time=native["time"], x=native["x_fraction"] * float(length),
                    fields={name: native[name] for name in ("active", "inactive", "length")},
                    active=native["active"], orientation=native["orientation"],
                    raw_boundaries=raw[:, -2:].copy())
    else:
        native_model = model_for(model)
        initial = np.asarray(state, dtype=float)
        if np.min(initial) < -1e-8:
            raise ValueError("restart violates native nonnegative invariant")
        native = n.simulate(native_model, initial_state=np.maximum(initial, 0.), length=float(length),
                    cells=int(cells), protocol=cue, times=times, end_time=float(times[-1]),
                    max_step=10. if max_step is None else float(max_step), backend="bdf")
        fields = {name: native.fields[name] for name in native_model.state_names}
        trace = SimpleNamespace(time=native.time, x=native.x, fields=fields,
                    active=fields[active_key(model)],
                    orientation=n.ho._orientation(native) if model == "holmes_model3" else n.rd._orientation(native))
    trace = _finish(trace, model, length)
    trace.restart_check = restart_error(state, trace)
    if not trace.restart_check["passed"]:
        raise ValueError("native output at time zero changed the restart state beyond roundoff")
    return trace


def restart_error(initial, trace):
    """Check the entire native restart against its first recorded full state."""
    first = state_at(trace, index=0)
    expected = initial if isinstance(initial, dict) else {"state": initial}
    observed = first if isinstance(first, dict) else {"state": first}
    details = {}
    for name, value in expected.items():
        value = np.asarray(value)
        actual = np.asarray(observed[name])
        if actual.shape != value.shape:
            raise ValueError("restart shape changed")
        allowance = 32 * np.finfo(float).eps * np.maximum(1., np.abs(value))
        difference = np.abs(actual - value)
        details[name] = {"maximum_absolute_error": float(difference.max()),
                         "maximum_roundoff_multiple": float(np.max(difference / allowance)),
                         "passed": bool(np.all(difference <= allowance))}
    return {"passed": all(item["passed"] for item in details.values()), "fields": details,
            "rule": "32 binary64 eps times max(1, abs(initial)) per component"}


def native_weights(trace, index=None):
    cells = trace.cells
    if trace.model == "spring":
        weights = np.repeat((trace.physical_length / (cells - 1))[:, None], cells, axis=1)
        weights[:, [0, -1]] *= .5
    else:
        weights = np.full((len(trace.time), cells), trace.length / cells)
    return weights if index is None else weights[index]


def profile_geometry(trace, index=-1):
    if trace.model == "spring":
        q = np.linspace(0., 1., trace.cells)
        return dict(sample_kind="nodes", reference_x=trace.x,
                    physical_x=q * trace.physical_length[index])
    edges = np.linspace(0., trace.length, trace.cells + 1)
    return dict(sample_kind="cells", reference_x=trace.x, reference_edges=edges,
                physical_x=trace.x, physical_edges=edges)


def save_phase(path, trace):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {name: np.asarray(getattr(trace, name)) for name in
            ("time", "x", "active", "orientation", "contrast", "left", "right", "physical_length")}
    data.update({"field_" + name: value for name, value in trace.fields.items()})
    data.update(model=np.array(trace.model), length=np.array(trace.length), cells=np.array(trace.cells))
    if hasattr(trace, "raw_boundaries"):
        data["raw_boundaries"] = trace.raw_boundaries
    if hasattr(trace, "restart_check"):
        data["restart_check_json"] = np.array(json.dumps(trace.restart_check))
    if trace.model == "debelly":
        data.update({"state_" + name: value for name, value in trace.state.items()})
    else:
        data["state"] = trace.state
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **data)
    temporary.replace(path)
    return {"path": str(path), "sha256": n.sha(path), "frames": len(trace.time),
            "fields": list(trace.fields), "full_native_fields": True}


def save_checkpoint(path, model, state, trace):
    """Save core-compatible checkpoint fields, retaining exact raw endpoints."""
    last = copy.copy(trace)
    for name in ("time", "active", "orientation", "contrast", "left", "right", "physical_length"):
        setattr(last, name, np.asarray(getattr(trace, name))[-1:])
    last.fields = {name: value[-1:] for name, value in trace.fields.items()}
    last.state, last.model = state, model
    if hasattr(trace, "raw_boundaries"):
        last.raw_boundaries = trace.raw_boundaries[-1:]
    return save_phase(path, last)


def read_checkpoint(path, model):
    state, old = n.read_checkpoint(Path(path), model)
    fields = old.fields
    with np.load(path, allow_pickle=False) as archive:
        length = float(archive["length"]) if "length" in archive else (
            float(old.x[-1]) if model == "spring" else float(old.x[0] + old.x[-1]))
        if model == "spring":
            if "raw_boundaries" not in archive and len(old.time) != 1:
                # A legacy checkpoint may retain diagnostic history but only
                # its terminal complete raw state. Use that exact terminal
                # state; never invent a moving-boundary history for the past.
                for name in ("time", "active", "orientation", "contrast", "left", "right"):
                    setattr(old, name, np.asarray(getattr(old, name))[-1:])
                old.fields = {name: value[-1:] for name, value in old.fields.items()}
                fields = old.fields
                old.checkpoint_history_terminal_only = True
            boundaries = (archive["raw_boundaries"].copy() if "raw_boundaries" in archive else
                          np.array([state[-2:]], dtype=float))
            if len(boundaries) != len(old.time):
                raise ValueError("full Spring checkpoint history lacks exact raw endpoint history")
            old.raw_boundaries = boundaries
    old.fields, old.state = fields, state
    return state, _finish(old, model, length)


load_checkpoint = read_checkpoint


def concatenate(traces):
    if not traces:
        raise ValueError("at least one phase required")
    result = copy.copy(traces[0])
    offset, time_parts = 0., []
    field_parts = {name: [] for name in result.fields}
    attribute_parts = {name: [] for name in ("active", "orientation")}
    raw = []
    for index, trace in enumerate(traces):
        selection = slice(None) if index == 0 else slice(1, None)
        time_parts.append(offset + trace.time[selection])
        for name in field_parts:
            field_parts[name].append(trace.fields[name][selection])
        for name in attribute_parts:
            attribute_parts[name].append(getattr(trace, name)[selection])
        if hasattr(trace, "raw_boundaries"):
            raw.append(trace.raw_boundaries[selection])
        offset += float(trace.time[-1])
    result.time = np.concatenate(time_parts)
    result.fields = {name: np.concatenate(parts) for name, parts in field_parts.items()}
    for name, parts in attribute_parts.items():
        setattr(result, name, np.concatenate(parts))
    if raw:
        result.raw_boundaries = np.concatenate(raw)
    return _finish(result, result.model, result.length)


def run_phases(model, length, state, phases, output_dir=None, cells=201,
               max_step=None, debelly_step=.5):
    """Restart from every exact phase endpoint; archive every full native field."""
    traces, events, offset = [], [], 0.
    for index, phase in enumerate(phases):
        duration = float(phase["duration"])
        if duration == 0:
            events.append(dict(index=index, name=phase.get("name", str(index)),
                               start=offset, stop=offset, zero_duration=True))
            continue
        times = phase.get("times")
        if times is None:
            times = sample_times(model, duration, phase.get("sample_interval", 2.5),
                                 debelly_step=debelly_step)
        trace = run_native(model, length, state, phase["cue"], times, cells=cells,
                           max_step=max_step, debelly_step=debelly_step)
        event = dict(index=index, name=phase.get("name", str(index)), start=offset,
                     stop=offset + duration, frames=len(trace.time), zero_duration=False,
                     restart_check=trace.restart_check)
        if hasattr(phase["cue"], "profile_diagnostics"):
            diagnostics = []
            for frame in (0, -1):
                domain = float(trace.physical_length[frame])
                coordinates = (np.linspace(0., domain, cells) if model == "spring" else trace.x)
                diagnostics.append(phase["cue"].profile_diagnostics(coordinates, domain))
            event["cue_initial_final"] = diagnostics
        if output_dir is not None:
            event["archive"] = save_phase(Path(output_dir) / f"phase_{index:03d}.npz", trace)
        state, offset = trace.state, offset + duration
        traces.append(trace)
        events.append(event)
    combined = concatenate(traces)
    if output_dir is not None:
        save_phase(Path(output_dir) / "trace.npz", combined)
    return combined, events


def native_field_tail_gate(trace, model=None, length=None, baseline_scales=None,
                           window=500., tolerance=.001):
    """Finite sampled full-field range/drift gate, independent of scalar fate."""
    model, length = model or trace.model, length or trace.length
    if baseline_scales is None:
        raise ValueError("inherited fixed native-field scales are required")
    if trace.time[-1] - trace.time[0] < window:
        return dict(field_tail_pass=None, available=False, reason="history_shorter_than_window")
    keep = trace.time >= trace.time[-1] - window - 1e-10
    times = trace.time[keep]
    if len(times) < 2:
        return dict(field_tail_pass=None, available=False, reason="fewer_than_two_tail_samples")
    fields = settling.monitored(trace, model, length)
    details = {}
    for name, values in fields.items():
        scale = float(baseline_scales[name])
        if scale <= 0 or not np.isfinite(values).all():
            raise ValueError("positive scale and finite full-field samples required")
        values = values[keep]
        details[name] = {
            "range_scaled": float(np.max(np.ptp(values, axis=0)) / scale),
            "projected_drift_500": float(np.max(abs(values[-1] - values[-2])) /
                                          (times[-1] - times[-2]) * window / scale),
        }
    maximum = max(max(value.values()) for value in details.values())
    return dict(field_tail_pass=bool(maximum <= tolerance), available=True,
                field_tail_max_scaled=maximum, field_tail_metrics=details,
                window=float(window), tolerance=float(tolerance),
                tail_samples=len(times), tail_start=float(times[0]), tail_stop=float(times[-1]),
                interpretation="finite sampled full-field diagnostic; independent of D fate")


def tail_gate(trace, model, length=None, scales=None, window=500., tolerance=.001):
    if isinstance(model, dict):
        case = model
        model, length, scales = case["model"], case["length"], case["fixed_scales"]
    return native_field_tail_gate(trace, model, length, scales, window, tolerance)


def no_cue_control_gate(trace, basal_trace, baseline_scales, tolerance=.01):
    """All saved full-field excursions from the exact basal, separately from D."""
    fields = settling.monitored(trace, trace.model, trace.length)
    reference = settling.monitored(basal_trace, basal_trace.model, basal_trace.length)
    detail = {}
    for name, values in fields.items():
        scale = float(baseline_scales[name])
        if scale <= 0 or not np.isfinite(values).all():
            raise ValueError("positive fixed scales and finite control fields required")
        detail[name] = float(np.max(np.abs(values - reference[name][-1])) / scale)
    maximum = max(detail.values())
    return dict(no_cue_control_pass=bool(maximum <= tolerance),
                no_cue_maximum_scaled_excursion=maximum,
                no_cue_scaled_excursions=detail, tolerance=float(tolerance),
                scope="every saved monitored full native field; scalar fate assessed separately")


sampling_times = sample_times


def select_preparation_endpoint(trials, manifest):
    """Select a stable native basal without claiming a minimum seed duration."""
    passed, best = establishment.qualifying(trials, manifest)
    if not passed:
        raise RuntimeError("Fresh preparation has no independently converged polarized candidate")
    selected = passed[0]
    ordered = sorted(trials, key=lambda item: item["duration"])
    lower = max((item["duration"] for item in trials if item["duration"] < selected["duration"]), default=0.)
    bracket = ((selected["duration"] - lower) / selected["duration"] if selected["duration"] else 0.)
    upper = ordered[-1]
    minimum_contrast = manifest["selection"]["minimum_contrast"]
    previous = [item["cue_aligned_strength"] for item in ordered[:-1]
                if item["verified"] and item["final_contrast"] >= minimum_contrast]
    upper_unresolved = bool(not upper["verified"] or (previous and upper["cue_aligned_strength"] >= .05 and
                           upper["cue_aligned_strength"] > max(previous) * 1.01))
    if upper_unresolved:
        raise RuntimeError("Fresh preparation upper observed-strength optimum remains unresolved")
    return dict(selected=selected, best_observed_strength=best,
                selected_relative_duration_bracket=bracket,
                preparation_duration_unresolved=bool(bracket > .010000001),
                upper_optimum_unresolved=False)


def prepare_basal(model, length=None, cells=201, output_dir=None, checkpoint_out=None):
    """Fresh exact Holmes/Spring basal using the existing establishment recipe.

    Every trial starts from native rest, applies the nominal relative-20% left
    cosine, and uses the existing convergence/holdout and endpoint selection.
    This function runs integrations and is for frozen SLURM workers only.
    """
    if isinstance(model, dict):
        specification = model
        # Also accept prepare_basal(case, checkpoint_out) for worker convenience.
        if length is not None and not isinstance(length, (int, float)):
            checkpoint_out = length
        model, length, cells = specification["model"], specification["length"], specification.get("cells", cells)
    if model not in ("holmes_model3", "spring"):
        raise ValueError("fresh preparation is scoped to Holmes/Spring; reuse exact retained other basals")
    if checkpoint_out is not None:
        checkpoint_out = Path(checkpoint_out)
        output_dir = output_dir or checkpoint_out.parent / (checkpoint_out.stem + "_preparation")
    if output_dir is None:
        raise ValueError("fresh campaign preparation output directory required")
    output_dir = Path(output_dir)
    n.protect_retained_output(output_dir)
    retained_corrected = (PACKAGE / "data/runs/debelly_corrected_three_test_20261008").resolve()
    if output_dir.resolve() == retained_corrected or retained_corrected in output_dir.resolve().parents:
        raise ValueError("Fresh preparations cannot overwrite corrected retained evidence")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = establishment.make_manifest()
    source = json.loads((PACKAGE / "data/runs/debelly_corrected_three_test_20261008/manifest.json").read_text())
    source_cases = [case for case in source["cases"] if case["model"] == model and case["eligible"]]
    scales = {name: max(float(case["fixed_scales"][name]) for case in source_cases)
              for name in source_cases[0]["fixed_scales"]}
    if model == "spring":
        scales["length"] = float(length)
    lo, hi = (.002, 50.) if model == "spring" else (.05, 5120.)
    case = dict(case_index=0, model=model, length=float(length), cells=int(cells),
                width=.2 * float(length), condition="relative_20pct", width_mode="relative",
                width_fraction_at_reference=.2, reference_length=1. if model == "spring" else 15.,
                amplitude_factor=1., amplitude=n.PROTOCOL["normalization"]["amplitude_by_model"][model],
                fixed_scales=scales, background_stimulus=0.,
                base_durations=[0.] + np.geomspace(lo, hi, 17).tolist())
    manifest["cases"] = [case]
    n.write(output_dir / "preparation_manifest.json", manifest)
    trials = []

    def evaluate(duration):
        result = establishment.trial(output_dir, case, duration, manifest)
        if not any(item["duration"] == duration for item in trials):
            trials.append(result)
        return result

    for duration in case["base_durations"]:
        evaluate(duration)
    for _ in range(manifest["refinement"]["upper_extensions"]):
        eligible = sorted([item for item in trials if item["verified"] and
                           item["final_contrast"] >= manifest["selection"]["minimum_contrast"]],
                          key=lambda item: item["duration"])
        if len(eligible) < 2 or eligible[-1]["cue_aligned_strength"] < .05:
            break
        last = eligible[-1]
        if last["cue_aligned_strength"] <= max(item["cue_aligned_strength"] for item in eligible[:-1]) * 1.01:
            break
        duration = min(last["duration"] * 4, manifest["refinement"]["duration_cap"])
        if duration <= max(item["duration"] for item in trials):
            break
        evaluate(duration)
    for _ in range(manifest["refinement"]["strength_rounds"]):
        for duration in establishment.strength_candidates(trials, manifest):
            evaluate(duration)
    for _ in range(manifest["refinement"]["threshold_steps"]):
        durations = establishment.threshold_candidates(trials, manifest)
        if not durations:
            break
        for duration in durations:
            evaluate(duration)
    try:
        endpoint_selection = select_preparation_endpoint(trials, manifest)
    except RuntimeError:
        _, best = establishment.qualifying(trials, manifest)
        receipt = dict(status="unqualified", qualified=False, model=model, length=float(length), cells=cells,
                       best_observed_strength=best, trials=trials, fixed_scales=scales)
        n.write(output_dir / "qualification.json", receipt)
        raise
    selected = endpoint_selection["selected"]
    bracket = endpoint_selection["selected_relative_duration_bracket"]
    source_checkpoint = output_dir / selected["checkpoint"]
    checkpoint_out = checkpoint_out or output_dir / "qualified_basal.npz"
    checkpoint_out.parent.mkdir(parents=True, exist_ok=True)
    if checkpoint_out.resolve() != source_checkpoint.resolve():
        shutil.copy2(source_checkpoint, checkpoint_out)
    receipt = dict(status="qualified", qualified=True, model=model, length=float(length), cells=int(cells),
        checkpoint=str(checkpoint_out), checkpoint_sha256=n.sha(checkpoint_out), fixed_scales=scales,
        preparation_status="polarized", selected=selected, trials=trials,
        selected_relative_duration_bracket=bracket, upper_optimum_unresolved=False,
        preparation_duration_unresolved=bool(bracket > .010000001),
        preparation_duration_interpretation="Declared shortest tested qualifying seed pulse; no minimum nonzero seed duration is claimed.",
        preparation_recipe="existing nominal relative20 establishment from independent native rest; unchanged kinetics",
        scale_rule="maximum retained same-model scales across sizes; Spring length scale equals requested L0",
        independent_review_required=True)
    n.write(output_dir / "qualification.json", receipt)
    return receipt


def qualify_existing_preparation(preparation_dir, output_dir, *, terminal_failure_path):
    """Qualify archived native preparations without changing/rerunning their source.

    This bounded recovery is for terminal workers rejected solely because a
    positive symmetry-breaking seed has an unresolved lower duration boundary.
    It independently checks archived convergence arithmetic, scalar histories,
    complete checkpoint hashes and selected native endpoint fields. Sparse
    preparation profile snapshots cannot replay every historical field gate;
    their source metrics remain explicit, and the new no-cue control is required.
    """
    preparation_dir, output_dir = Path(preparation_dir), Path(output_dir)
    failure_path = Path(terminal_failure_path)
    failure = json.loads(failure_path.read_text())
    if "Fresh preparation duration boundary or upper optimum remains unresolved" not in failure.get("traceback", ""):
        raise ValueError("Requalification is scoped to the recorded terminal preparation-boundary rejection")
    if output_dir.exists():
        raise FileExistsError("Choose a new requalification directory; source evidence stays immutable")
    n.protect_retained_output(output_dir)
    manifest_path = preparation_dir / "preparation_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    case = manifest["cases"][0]
    if len(manifest["cases"]) != 1 or case["model"] not in ("holmes_model3", "spring"):
        raise ValueError("Expected exactly one fresh Holmes/Spring preparation")
    paths = list((preparation_dir / "work/cases/case_0000").glob("duration_*.json"))
    trials = sorted([json.loads(path.read_text()) for path in paths], key=lambda row: row["duration"])
    if not trials or len({row["duration"] for row in trials}) != len(trials):
        raise ValueError("Missing or duplicate trial durations")
    if not set(case["base_durations"]).issubset({row["duration"] for row in trials}):
        raise ValueError("Original declared preparation scan is incomplete")
    selection, steady = manifest["selection"], manifest["steady_state"]
    checks, eligible = [], []
    source_hashes = {str(manifest_path): n.sha(manifest_path), str(failure_path): n.sha(failure_path)}

    def persistent(time, values, target, tolerance):
        bad = np.flatnonzero(~np.isfinite(values) | (abs(values - target) > tolerance))
        index = int(bad[-1] + 1) if len(bad) else 0
        if index >= len(time) or time[-1] - time[index] < 1000.:
            return None
        return float(time[index])

    for path, trial in [(next(path for path in paths if json.loads(path.read_text())["duration"] == item["duration"]), item)
                        for item in trials]:
        source_hashes[str(path)] = n.sha(path)
        verified = bool(trial["verified"])
        detail = dict(duration=trial["duration"], recorded_verified=verified)
        if not verified:
            checks.append(detail)
            continue
        convergence = trial["convergence"]
        history = convergence["checks"]
        if not convergence["verified"] or len(history) < steady["consecutive_windows"] + 1:
            raise ValueError("Verified trial lacks convergence/holdout evidence")
        last = history[-1]
        for item in history:
            computed_steady = all(max(metric.values()) <= steady["tolerance"]
                                  for metric in item["metrics"].values())
            if computed_steady != item["steady"]:
                raise ValueError("Stored convergence boolean disagrees with per-field metrics")
        previous = history[-1 - steady["consecutive_windows"]:-1]
        if (last["kind"] != "holdout" or not last["verified"] or not last["steady"]
                or last["maximum_scaled_excursion"] > steady["holdout_tolerance"]
                or any(item["kind"] != "convergence" or not item["steady"] for item in previous)):
            raise ValueError("Final two-window/holdout criteria are not proved")
        if last["time"] - previous[-1]["time"] != steady["holdout"]:
            raise ValueError("Holdout duration differs from the declared native rule")
        trace_path = path.with_suffix(".npz")
        checkpoint = preparation_dir / trial["checkpoint"]
        if n.sha(checkpoint) != trial["checkpoint_sha256"]:
            raise ValueError("Trial checkpoint digest mismatch")
        source_hashes[str(trace_path)], source_hashes[str(checkpoint)] = n.sha(trace_path), n.sha(checkpoint)
        state, endpoint = read_checkpoint(checkpoint, case["model"])
        monitored = settling.monitored(endpoint, case["model"], case["length"])
        with np.load(trace_path, allow_pickle=False) as archive:
            if any(not np.isfinite(archive[name]).all() for name in archive.files):
                raise ValueError("Preparation archive has nonfinite arrays")
            time = archive["release_time"]
            orientation, contrast = archive["release_orientation"], archive["release_contrast"]
            distances = archive["maximum_scaled_field_distance"]
            if np.any(np.diff(time) <= 0) or time[-1] != convergence["elapsed"]:
                raise ValueError("Release native time/elapsed evidence differs")
            np.testing.assert_allclose(orientation[-1], trial["final_orientation"], rtol=0, atol=1e-14)
            np.testing.assert_allclose(contrast[-1], trial["final_contrast"], rtol=0, atol=1e-14)
            for name, values in monitored.items():
                np.testing.assert_array_equal(archive["endpoint_" + name], values[-1])
                np.testing.assert_array_equal(archive["profile_" + name][-1], values[-1])
            arrival = persistent(time, orientation, trial["final_orientation"], trial["polarity_arrival_tolerance"])
            field_arrival = persistent(time, distances, 0., .01)
            if arrival != trial["polarity_arrival"] or field_arrival != trial["field_arrival"]:
                raise ValueError("Recorded persistent arrival differs from full saved scalar history")
        detail.update(convergence_arithmetic_pass=True, checkpoint_hash_pass=True,
            checkpoint_native_endpoint_fields_pass=True, scalar_arrival_recomputed=True,
            holdout_maximum_scaled_excursion=last["maximum_scaled_excursion"],
            final_field_metrics=last["metrics"], polarity_arrival=arrival)
        checks.append(detail)
        if (trial["final_contrast"] >= selection["minimum_contrast"]
                and trial["polarity_arrival"] is not None):
            eligible.append(trial)
    if not eligible:
        raise ValueError("No stable observed preparation endpoint")
    best = max(row["cue_aligned_strength"] for row in eligible)
    candidates = [row for row in eligible if best >= selection["minimum_strength"]
                  and row["cue_aligned_strength"] >= best * (1 - selection["relative_strength_tolerance"])]
    if not candidates:
        raise ValueError("No qualifying left-polarized endpoint")
    selected = min(candidates, key=lambda row: row["duration"])
    lower = max((row["duration"] for row in trials if row["duration"] < selected["duration"]), default=0.)
    bracket = (selected["duration"] - lower) / selected["duration"] if selected["duration"] else 0.
    upper = trials[-1]
    earlier = [row["cue_aligned_strength"] for row in trials[:-1]
               if row["verified"] and row["final_contrast"] >= selection["minimum_contrast"]]
    upper_unresolved = bool(not upper["verified"] or (earlier and upper["cue_aligned_strength"] >= .05
                            and upper["cue_aligned_strength"] > max(earlier) * 1.01))
    if upper_unresolved:
        raise ValueError("This is not solely a lower seed-duration-boundary rejection")
    checkpoint = preparation_dir / selected["checkpoint"]
    output_dir.mkdir(parents=True)
    destination = output_dir / "qualified_basal.npz"
    shutil.copy2(checkpoint, destination)
    shutil.copy2(Path(__file__), output_dir / "qualification_source.py")
    if n.sha(destination) != selected["checkpoint_sha256"]:
        raise ValueError("Requalified checkpoint changed during copy")
    receipt = dict(generated_utc=n.now(), status="qualified", qualified=True, model=case["model"], length=case["length"],
        cells=case["cells"], checkpoint=str(destination.resolve()), checkpoint_sha256=n.sha(destination),
        fixed_scales=case["fixed_scales"], preparation_status="polarized", selected=selected,
        trials=trials, selected_relative_duration_bracket=bracket, upper_optimum_unresolved=False,
        preparation_duration_unresolved=bool(bracket > .010000001),
        preparation_duration_interpretation="Declared shortest tested qualifying seed pulse; no minimum nonzero seed duration is claimed.",
        requalification_integrations=0, original_terminal_failure=str(failure_path.resolve()),
        original_failure_preserved=True, source_preparation_preserved=True,
        source_hashes=source_hashes, independent_record_checks=checks,
        audit_scope="All verified trial native checkpoint digests/endpoints, full scalar arrival histories, and per-field convergence/holdout metric arithmetic.",
        audit_limitation="Sparse historical preparation profiles cannot independently replay every field sample. The complete source metrics are checked and a separate new full-field no-cue control remains required.",
        no_cue_control_required=True, independent_review_required=True,
        requalification_source_sha256=n.sha(output_dir / "qualification_source.py"))
    n.write(output_dir / "qualification.json", receipt)
    return receipt
