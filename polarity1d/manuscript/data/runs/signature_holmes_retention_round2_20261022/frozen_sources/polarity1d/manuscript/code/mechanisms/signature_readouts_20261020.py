"""Native-coordinate, positive-excess readouts for the Figure 2 signature assays.

This module runs no solver. Primary replacement uses a disjoint reference-half
partition. Regional cue amounts are separate, potentially overlapping diagnostics.
Cell samples are piecewise constant; nodal samples are piecewise linear. Both
receive exact fractional support integration with optional moving physical weights.
The current profile minimum is a descriptive floor, not a conserved-species pool.
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np


def _vector(value, name):
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or not len(result) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a nonempty finite one-dimensional array")
    return result


def _increasing(value, name):
    value = _vector(value, name)
    if len(value) < 2 or not np.all(np.diff(value) > 0):
        raise ValueError(f"{name} must be strictly increasing with at least two points")
    return value


def _geometry(reference_x, active, sample_kind, reference_edges,
              physical_x, physical_edges):
    x = _vector(reference_x, "reference_x")
    y = _vector(active, "active")
    if x.shape != y.shape:
        raise ValueError("reference_x and active must have identical shape")
    if sample_kind == "cells":
        if reference_edges is None:
            raise ValueError("Cell samples require explicit reference_edges")
        ref = _increasing(reference_edges, "reference_edges")
        if len(ref) != len(y) + 1:
            raise ValueError("reference_edges must contain one more value than active")
        if not np.all((x > ref[:-1]) & (x < ref[1:])):
            raise ValueError("Each cell sample must lie inside its reference edges")
        if physical_x is not None and physical_edges is None:
            raise ValueError("Cell physical weights require physical_edges, not only centres")
        phys = ref if physical_edges is None else _increasing(physical_edges, "physical_edges")
    elif sample_kind == "nodes":
        ref = _increasing(x, "reference_x")
        if physical_edges is not None:
            raise ValueError("Nodal physical weights require physical_x")
        phys = ref if physical_x is None else _increasing(physical_x, "physical_x")
    else:
        raise ValueError("sample_kind must be 'nodes' or 'cells'")
    if phys.shape != ref.shape:
        raise ValueError("Physical and reference coordinate arrays must have identical shape")
    return x, y, ref, phys


def _intervals(intervals, lo, hi):
    """Union clipped reference intervals; touching intervals count only once."""
    if intervals is None:
        intervals = [(lo, hi)]
    items = np.asarray(intervals, dtype=float)
    if items.size == 0:
        return []
    if items.ndim == 1 and items.shape == (2,):
        items = items.reshape(1, 2)
    if items.ndim != 2 or items.shape[1] != 2 or not np.all(np.isfinite(items)):
        raise ValueError("Intervals must be finite pairs of reference coordinates")
    if np.any(items[:, 1] < items[:, 0]):
        raise ValueError("Interval endpoints must be ordered")
    clipped = sorted((max(lo, float(a)), min(hi, float(b)))
                     for a, b in items if min(hi, b) > max(lo, a))
    union = []
    for a, b in clipped:
        if union and a <= union[-1][1]:
            union[-1][1] = max(union[-1][1], b)
        else:
            union.append([a, b])
    return union


def _integral(y, ref, phys, intervals, sample_kind):
    total = 0.0
    widths = np.diff(ref)
    physical_widths = np.diff(phys)
    for a, b in _intervals(intervals, ref[0], ref[-1]):
        left = np.maximum(ref[:-1], a)
        right = np.minimum(ref[1:], b)
        span = np.maximum(right - left, 0)
        if sample_kind == "cells":
            total += np.sum(y * physical_widths * span / widths)
        else:
            slope = np.diff(y) / widths
            y_left = y[:-1] + slope * (left - ref[:-1])
            y_right = y[:-1] + slope * (right - ref[:-1])
            total += np.sum((y_left + y_right) * .5 * physical_widths * span / widths)
    return float(total)


def integrate_profile(reference_x, values, intervals=None, *, sample_kind="nodes",
                      reference_edges=None, physical_x=None, physical_edges=None):
    """Integrate a profile on the union of reference intervals, in physical units.

    Reference coordinates define support membership. Physical coordinates supply
    the integration Jacobian; shifting a whole physical domain changes no amount.
    Overlapping input intervals are unioned before integration.
    """
    _, y, ref, phys = _geometry(reference_x, values, sample_kind, reference_edges,
                               physical_x, physical_edges)
    return _integral(y, ref, phys, intervals, sample_kind)


def _overlap(first, second):
    return float(sum(max(0.0, min(b, d) - max(a, c))
                     for a, b in first for c, d in second))


def _share(challenger, incumbent, denominator_atol):
    total = challenger + incumbent
    if total < denominator_atol or total <= 0:
        return None, None
    share = challenger / total
    return float(2 * share - 1), float(share)


def _halfheight_intervals(y, ref, sample_kind):
    threshold = .5 * float(np.max(y))
    if threshold <= 0:
        return []
    spans = []
    if sample_kind == "cells":
        spans = [(ref[i], ref[i + 1]) for i in range(len(y)) if y[i] >= threshold]
    else:
        for i, (first, last) in enumerate(zip(y[:-1], y[1:])):
            left, right = ref[i], ref[i + 1]
            if first >= threshold and last >= threshold:
                spans.append((left, right))
            elif first >= threshold or last >= threshold:
                crossing = left + (threshold - first) * (right - left) / (last - first)
                spans.append((left, crossing) if first >= threshold else (crossing, right))
    return _intervals(spans, ref[0], ref[-1])


def describe_basal(reference_x, active, *, sample_kind="nodes", reference_edges=None,
                   physical_x=None, physical_edges=None, incumbent_side="left",
                   denominator_atol=1e-12, contrast_atol=.05):
    """Freeze the connected basal half-height support containing its main peak.

    Half-height is measured above the current basal minimum. Multiple disconnected
    half-height regions remain in the descriptor, so a multi-site basal is visible.
    """
    if incumbent_side not in ("left", "right"):
        raise ValueError("incumbent_side must be 'left' or 'right'")
    x, y, ref, phys = _geometry(reference_x, active, sample_kind, reference_edges,
                               physical_x, physical_edges)
    floor = float(np.min(y))
    excess = np.maximum(y - floor, 0)
    intervals = _halfheight_intervals(excess, ref, sample_kind)
    max_indices = np.flatnonzero(excess == np.max(excess))
    peak_index = int(max_indices[0] if incumbent_side == "left" else max_indices[-1])
    peak_x = float(x[peak_index])
    incumbent = [span for span in intervals if span[0] <= peak_x <= span[1]]
    total = _integral(excess, ref, phys, None, sample_kind)
    scale = abs(float(np.max(y))) + abs(floor)
    contrast = float(np.max(excess) / scale) if scale > 0 else 0.0
    mid = .5 * (ref[0] + ref[-1])
    expected_half = peak_x <= mid if incumbent_side == "left" else peak_x >= mid
    return {
        "incumbent_intervals": incumbent,
        "all_halfheight_intervals": intervals,
        "incumbent_side": incumbent_side,
        "peak_reference_x": peak_x,
        "current_min_floor": floor,
        "active_peak": float(np.max(y)),
        "excess_peak": float(np.max(excess)),
        "E_total": total,
        "E_incumbent_fixed": _integral(excess, ref, phys, incumbent, sample_kind),
        "incumbent_extent_reference": float(sum(b - a for a, b in incumbent)),
        "incumbent_extent_physical": _integral(np.ones_like(y), ref, phys, incumbent, sample_kind),
        "contrast": contrast,
        "contrast_eligible": bool(contrast >= contrast_atol),
        "amount_eligible": bool(total > denominator_atol),
        "peak_in_expected_half": bool(expected_half),
        "multiple_halfheight_regions": len(intervals) > 1,
        "sample_kind": sample_kind,
        "reference_domain": [float(ref[0]), float(ref[-1])],
        "physical_domain": [float(phys[0]), float(phys[-1])],
        "denominator_atol": float(denominator_atol),
        "contrast_atol": float(contrast_atol),
    }


def measure_profile(reference_x, active, *, incumbent_intervals,
                    challenger_own_intervals, challenger_fixed_intervals,
                    sample_kind="nodes", reference_edges=None, physical_x=None,
                    physical_edges=None, challenger_side="right",
                    denominator_atol=1e-12, contrast_atol=.05,
                    basal_E_total=None, basal_contrast=None,
                    amount_fraction_of_basal=.01, contrast_fraction_of_basal=.1):
    """Return primary half-partition fate and separate regional amount diagnostics."""
    if challenger_side not in ("left", "right"):
        raise ValueError("challenger_side must be 'left' or 'right'")
    if denominator_atol < 0 or not 0 <= contrast_atol <= 1:
        raise ValueError("Require nonnegative amount tolerance and contrast tolerance in [0,1]")
    if amount_fraction_of_basal < 0 or contrast_fraction_of_basal < 0:
        raise ValueError("Basal-scaled eligibility fractions must be nonnegative")
    if basal_E_total is not None and (not np.isfinite(basal_E_total) or basal_E_total < 0):
        raise ValueError("basal_E_total must be finite and nonnegative")
    if basal_contrast is not None and (not np.isfinite(basal_contrast) or not 0 <= basal_contrast <= 1):
        raise ValueError("basal_contrast must be finite and within [0,1]")
    amount_threshold = max(denominator_atol, amount_fraction_of_basal * (basal_E_total or 0))
    contrast_threshold = max(contrast_atol, contrast_fraction_of_basal * (basal_contrast or 0))
    _, y, ref, phys = _geometry(reference_x, active, sample_kind, reference_edges,
                               physical_x, physical_edges)
    floor, peak = float(np.min(y)), float(np.max(y))
    excess = np.maximum(y - floor, 0)
    midpoint = .5 * (ref[0] + ref[-1])
    left = _integral(excess, ref, phys, [(ref[0], midpoint)], sample_kind)
    right = _integral(excess, ref, phys, [(midpoint, ref[-1])], sample_kind)
    total = _integral(excess, ref, phys, None, sample_kind)
    incumbent = _intervals(incumbent_intervals, ref[0], ref[-1])
    own = _intervals(challenger_own_intervals, ref[0], ref[-1])
    fixed = _intervals(challenger_fixed_intervals, ref[0], ref[-1])
    regional_incumbent = _integral(excess, ref, phys, incumbent, sample_kind)
    regional_own = _integral(excess, ref, phys, own, sample_kind)
    regional_fixed = _integral(excess, ref, phys, fixed, sample_kind)
    challenger_half, incumbent_half = (right, left) if challenger_side == "right" else (left, right)
    dominance_raw, share_raw = _share(challenger_half, incumbent_half, denominator_atol)
    dominance, share = _share(challenger_half, incumbent_half, amount_threshold)
    D_own, share_own = _share(regional_own, regional_incumbent, amount_threshold)
    D_fixed, share_fixed = _share(regional_fixed, regional_incumbent, amount_threshold)
    own_overlap, fixed_overlap = _overlap(incumbent, own), _overlap(incumbent, fixed)
    scale = abs(peak) + abs(floor)
    contrast = (peak - floor) / scale if scale > 0 else 0.0
    amount_ok = total >= amount_threshold and total > 0
    contrast_ok = contrast >= contrast_threshold
    eligible = amount_ok and contrast_ok
    reasons = []
    if not amount_ok:
        reasons.append("positive_excess_below_frozen_basal_amount_gate")
    if not contrast_ok:
        reasons.append("contrast_below_frozen_basal_contrast_gate")
    return {
        "D": dominance if eligible else None,
        "D_half": dominance if eligible else None,
        "share_half": share if eligible else None,
        "D_half_raw": dominance_raw, "share_half_raw": share_raw,
        "E_left": left, "E_right": right, "E_total": total,
        "E_challenger_half": challenger_half, "E_incumbent_half": incumbent_half,
        "E_incumbent_fixed": regional_incumbent,
        "E_challenger_own": regional_own, "E_challenger_fixed": regional_fixed,
        "D_roi_own": D_own, "D_roi_fixed": D_fixed,
        "share_roi_own": share_own, "share_roi_fixed": share_fixed,
        "roi_own_overlap": own_overlap > 0, "roi_fixed_overlap": fixed_overlap > 0,
        "roi_own_overlap_reference_extent": own_overlap,
        "roi_fixed_overlap_reference_extent": fixed_overlap,
        "roi_own_eligible": bool(eligible and own_overlap == 0),
        "roi_fixed_eligible": bool(eligible and fixed_overlap == 0),
        "current_min_floor": floor, "active_peak": peak,
        "excess_peak": peak - floor, "raw_active_amount": _integral(y, ref, phys, None, sample_kind),
        "contrast": float(contrast), "contrast_eligible": bool(contrast_ok),
        "amount_eligible": bool(amount_ok), "D_half_eligible": bool(eligible),
        "homogeneous": bool(peak == floor), "negative_active_present": bool(floor < 0),
        "replacement": bool(dominance >= .5) if eligible else False,
        "eligible": bool(eligible),
        "readout_eligible": bool(eligible),
        "ineligibility_reason": ";".join(reasons) if reasons else None,
        "partition_sum_error": float(left + right - total),
        "physical_length": float(phys[-1] - phys[0]),
        "denominator_atol": float(denominator_atol), "contrast_atol": float(contrast_atol),
        "amount_threshold": float(amount_threshold),
        "contrast_threshold": float(contrast_threshold),
        "basal_E_total": None if basal_E_total is None else float(basal_E_total),
        "basal_contrast": None if basal_contrast is None else float(basal_contrast),
        "amount_fraction_of_basal": float(amount_fraction_of_basal),
        "contrast_fraction_of_basal": float(contrast_fraction_of_basal),
    }


def first_crossing(time, values, threshold, *, direction="up"):
    """First sampled threshold hit with a bracket; interpolation is only an estimate.

    Missing values never bridge an event. An initial or post-missing successful
    sample is left-censored, because its preceding threshold state is unknown.
    """
    t = _vector(time, "time")
    v = np.asarray(values, dtype=float)
    if v.shape != t.shape or np.any(np.diff(t) <= 0):
        raise ValueError("time and values require equal shape and increasing time")
    if direction not in ("up", "down"):
        raise ValueError("direction must be 'up' or 'down'")
    hits = v >= threshold if direction == "up" else v <= threshold
    hits &= np.isfinite(v)
    positions = np.flatnonzero(hits)
    if not len(positions):
        return {"observed": False, "threshold": float(threshold), "direction": direction,
                "lower_time": None, "upper_time": None, "estimated_time": None,
                "left_censored": False, "right_censored_at": float(t[-1])}
    i = int(positions[0])
    previous_known = i > 0 and np.isfinite(v[i - 1])
    estimate = None
    if previous_known:
        estimate = float(t[i - 1] + (threshold - v[i - 1]) *
                         (t[i] - t[i - 1]) / (v[i] - v[i - 1]))
    return {"observed": True, "threshold": float(threshold), "direction": direction,
            "lower_time": float(t[i - 1]) if previous_known else None,
            "upper_time": float(t[i]), "estimated_time": estimate,
            "sample_index": i, "left_censored": not previous_known,
            "right_censored_at": None}


def _frame(value, index, count, name):
    array = np.asarray(value, dtype=float)
    if array.ndim == 1:
        return array
    if array.ndim == 2 and array.shape[0] == count:
        return array[index]
    raise ValueError(f"{name} must be static one-dimensional or have one row per time")


def _frame_geometry(trace, index, count):
    result = {"sample_kind": trace.get("sample_kind", "nodes")}
    for key in ("reference_edges", "physical_x", "physical_edges"):
        if trace.get(key) is not None:
            result[key] = _frame(trace[key], index, count, key)
    return result


def _physical_width_interval(x, geometry, width, side):
    """Map an endpoint-anchored physical width into this frame's reference map."""
    if not np.isfinite(width) or width < 0:
        raise ValueError("Physical challenger support width must be finite and nonnegative")
    if geometry.get("sample_kind", "nodes") == "cells":
        ref = np.asarray(geometry["reference_edges"], dtype=float)
        phys = np.asarray(geometry.get("physical_edges", ref), dtype=float)
    else:
        ref = np.asarray(x, dtype=float)
        phys = np.asarray(geometry.get("physical_x", ref), dtype=float)
    if side == "right":
        start = float(np.interp(max(phys[0], phys[-1] - width), phys, ref))
        return [[start, float(ref[-1])]]
    end = float(np.interp(min(phys[-1], phys[0] + width), phys, ref))
    return [[float(ref[0]), end]]


def _case_intervals(value, index, count, name):
    array = np.asarray(value, dtype=float)
    if array.ndim == 3:
        if array.shape[0] != count or array.shape[2] != 2:
            raise ValueError(f"{name} frame-specific intervals require shape (time, intervals, 2)")
        return array[index].tolist()
    return array.tolist()


def analyze_trace(trace: Mapping[str, Any], case: Mapping[str, Any], basal_trace=None):
    """Analyze a native trace into JSON-compatible scalar series and final fate.

    Required trace keys: time, x (frozen reference coordinates), active (T,N).
    Cells additionally require reference_edges. Optional physical_x/physical_edges
    can be static or frame-specific. The case supplies frozen reference intervals
    challenger_own_intervals/challenger_fixed_intervals and optional incumbent
    intervals, tolerances and challenger_side. No field-tail qualification is
    inferred here; native engine diagnostics remain separate.
    """
    t = _vector(trace["time"], "time")
    if np.any(np.diff(t) <= 0):
        raise ValueError("Trace time must be strictly increasing")
    active = np.asarray(trace["active"], dtype=float)
    if active.ndim != 2 or active.shape[0] != len(t):
        raise ValueError("Trace active must have shape (time, spatial samples)")
    x = _vector(trace["x"] if "x" in trace else trace["xref"], "x")
    side = case.get("challenger_side", "right")
    tolerance = {"denominator_atol": case.get("denominator_atol", 1e-12),
                 "contrast_atol": case.get("contrast_atol", .05)}
    baseline = trace if basal_trace is None else basal_trace
    baseline_active = np.asarray(baseline["active"], dtype=float)
    baseline_count = 1 if baseline_active.ndim == 1 else len(baseline_active)
    baseline_index = 0 if baseline is trace else baseline_count - 1
    baseline_profile = baseline_active if baseline_active.ndim == 1 else baseline_active[baseline_index]
    baseline_x = baseline["x"] if "x" in baseline else baseline["xref"]
    baseline_geometry = _frame_geometry(baseline, baseline_index, baseline_count)
    descriptor = describe_basal(
        baseline_x, baseline_profile,
        incumbent_side="left" if side == "right" else "right",
        **baseline_geometry, **tolerance)
    basal_scale = {
        "basal_E_total": descriptor["E_total"],
        "basal_contrast": descriptor["contrast"],
        "amount_fraction_of_basal": case.get("amount_fraction_of_basal", .01),
        "contrast_fraction_of_basal": case.get("contrast_fraction_of_basal", .1),
    }
    incumbent = case.get("incumbent_intervals", descriptor["incumbent_intervals"])
    own, fixed = case["challenger_own_intervals"], case["challenger_fixed_intervals"]
    fixed_material = (_physical_width_interval(
        baseline_x, baseline_geometry, case["challenger_fixed_width_physical"], side)
        if "challenger_fixed_width_physical" in case else None)
    rows, own_by_frame, fixed_by_frame = [], [], []
    for i, profile in enumerate(active):
        geometry = _frame_geometry(trace, i, len(t))
        own_current = (_physical_width_interval(x, geometry, case["challenger_own_width_physical"], side)
                       if "challenger_own_width_physical" in case else
                       _case_intervals(own, i, len(t), "challenger_own_intervals"))
        fixed_current = (fixed_material if fixed_material is not None else
                         _case_intervals(fixed, i, len(t), "challenger_fixed_intervals"))
        own_by_frame.append(own_current)
        fixed_by_frame.append(fixed_current)
        rows.append(measure_profile(
            x, profile, incumbent_intervals=incumbent,
            challenger_own_intervals=own_current, challenger_fixed_intervals=fixed_current,
            challenger_side=side, **geometry, **tolerance, **basal_scale))
    series = {"time": t.tolist()}
    series.update({key: [row[key] for row in rows] for key in rows[0]})
    event_values = [row["D_half"] if row["D_half_eligible"] else None for row in rows]
    result = {"series": series, "basal_descriptor": descriptor,
              "first_continuous_D_half_ge_0_5_event": first_crossing(t, event_values, .5),
              "final": rows[-1], "final_replacement": rows[-1]["replacement"],
              "final_d": rows[-1]["D_half"], "replacement": rows[-1]["replacement"],
              "eligible": rows[-1]["eligible"],
              "readout_eligible": rows[-1]["readout_eligible"],
              "ineligibility_reason": rows[-1]["ineligibility_reason"],
              "primary_D_definition": "Disjoint full-reference-domain halves; current-min positive excess; physical integration weights.",
              "field_tail_gate": "Not evaluated by readouts; use independent native-field diagnostic.",
              "incumbent_intervals": incumbent,
              "challenger_own_intervals": own,
              "challenger_fixed_intervals": fixed,
              "challenger_own_intervals_by_frame": own_by_frame,
              "challenger_fixed_intervals_by_frame": fixed_by_frame}
    return result


# The short alias is convenient for phase engines that need only scalar series.
measure_trace = analyze_trace
