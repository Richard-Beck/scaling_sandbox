"""Bounded duration calibration for the 20 October signature-assay campaign.

This module does not define dominance, run a solver, or choose native units.
The native worker supplies a post-withdrawal readout. Field stability remains
an independent diagnostic and never changes the replacement classification.
The returned threshold is a sampled bracket, not proof of global minimality.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Mapping, Any


@dataclass(frozen=True)
class DurationSearch:
    horizon: float
    scan_durations: tuple[float, ...]
    absolute_tolerance: float
    relative_tolerance: float = 0.01
    native_step: float | None = None
    max_bisections: int = 32
    dominance_threshold: float = 0.5


def calibrate(config: DurationSearch,
              evaluate: Callable[[float], Mapping[str, Any]]) -> dict[str, Any]:
    """Scan, then refine the earliest observed failure-to-success transition.

    ``evaluate`` must return ``final_d`` after withdrawal/recovery and a boolean
    ``field_tail_pass``. Explicit low-activity ineligibility permits null D only
    with ``readout_eligible=False`` and an ``ineligibility_reason``. It may include archive hashes and native gate details;
    these are retained verbatim. Every requested duration must be realizable on
    the declared native lattice. A missing/nonfinite dominance aborts the search
    instead of being silently counted as failure. Evaluation errors propagate.
    """
    if not math.isfinite(config.horizon) or config.horizon <= 0:
        raise ValueError("A finite positive duration horizon is required")
    if (not math.isfinite(config.absolute_tolerance)
            or config.absolute_tolerance <= 0
            or not math.isfinite(config.relative_tolerance)
            or config.relative_tolerance < 0
            or config.max_bisections < 1):
        raise ValueError("Invalid duration-search tolerance or iteration bound")
    step = config.native_step
    if step is not None and (not math.isfinite(step) or step <= 0):
        raise ValueError("Native lattice spacing must be finite and positive")

    def valid_time(value: float) -> float:
        value = float(value)
        if not math.isfinite(value) or not 0 <= value <= config.horizon:
            raise ValueError("Scan duration lies outside the finite horizon")
        if step is not None:
            tick = round(value / step)
            if not math.isclose(value, tick * step, rel_tol=0, abs_tol=1e-9):
                raise ValueError("Duration does not lie on the native lattice")
            value = tick * step
        return value

    durations = sorted({valid_time(t) for t in
                        (0.0, *config.scan_durations, config.horizon)})
    records: dict[float, dict[str, Any]] = {}

    def run(duration: float) -> dict[str, Any]:
        duration = valid_time(duration)
        if duration not in records:
            result = dict(evaluate(duration))
            eligible = result.get("readout_eligible", True)
            if type(eligible) is not bool:
                raise ValueError("Readout eligibility must be boolean")
            if result["final_d"] is None and not eligible and result.get("ineligibility_reason"):
                final_d = None
            else:
                final_d = float(result["final_d"])
                if not math.isfinite(final_d) or not -1.000000001 <= final_d <= 1.000000001:
                    raise ValueError("Post-withdrawal dominance must be finite in [-1, 1]")
            if not eligible and not result.get("ineligibility_reason"):
                raise ValueError("Ineligible readouts require an explicit reason")
            if type(result["field_tail_pass"]) is not bool:
                raise ValueError("Record a separate boolean native-field-tail gate")
            if ("actual_duration" in result and
                    not math.isclose(float(result["actual_duration"]), duration,
                                     rel_tol=0, abs_tol=1e-9)):
                raise ValueError("Native worker changed the requested duration")
            result.update(duration=duration, final_d=final_d, readout_eligible=eligible,
                          replacement=bool(eligible and final_d is not None
                                           and final_d >= config.dominance_threshold))
            records[duration] = result
        return records[duration]

    # Scan the complete declared grid even after a first success, so a later
    # observed loss of commitment remains visible rather than being discarded.
    for duration in durations:
        run(duration)
    successes = [t for t in durations if records[t]["replacement"]]
    lower = upper = None
    refinements = 0
    resolved = False
    if records[0.0]["replacement"]:
        status = "zero_duration_already_replaces_invalid_incumbent_control"
    elif not successes:
        status = "no_sampled_commitment_through_finite_horizon"
    else:
        upper = successes[0]
        lower = durations[durations.index(upper) - 1]
        for _ in range(config.max_bisections):
            tolerance = max(config.absolute_tolerance,
                            config.relative_tolerance * upper)
            if upper - lower <= tolerance:
                resolved = True
                break
            middle = (lower + upper) / 2.0
            if step is not None:
                low_tick, high_tick = round(lower / step), round(upper / step)
                if high_tick - low_tick <= 1:
                    resolved = True
                    break
                middle = ((low_tick + high_tick) // 2) * step
            if not lower < middle < upper:
                break
            if run(middle)["replacement"]:
                upper = middle
            else:
                lower = middle
            refinements += 1
        tolerance = max(config.absolute_tolerance,
                        config.relative_tolerance * upper)
        resolved = resolved or upper - lower <= tolerance
        if step is not None:
            resolved = resolved or round((upper - lower) / step) <= 1
        status = ("sampled_commitment_bracket_resolved" if resolved
                  else "sampled_commitment_bracket_iteration_limit")

    ordered = [records[t] for t in sorted(records)]
    transitions = [
        dict(lower=a["duration"], upper=b["duration"],
             from_replacement=a["replacement"], to_replacement=b["replacement"])
        for a, b in zip(ordered, ordered[1:])
        if a["replacement"] != b["replacement"]
    ]
    reversal = any(t["from_replacement"] and not t["to_replacement"]
                   for t in transitions)
    return dict(
        status=status, lower_failing_duration=lower,
        upper_committing_duration=upper, bracket_resolved=resolved,
        dominance_threshold=config.dominance_threshold,
        fate_definition="Post-withdrawal D at the declared recovery horizon",
        field_tail_used_to_define_fate=False,
        native_step=step, finite_duration_horizon=config.horizon,
        observed_nonmonotone_commitment=reversal,
        train_ready=bool(resolved and not reversal),
        scan_durations=durations, bisections=refinements,
        sampled_transitions=transitions, evaluations=ordered,
        limitation="Earliest observed success bracket only; unsampled success islands and global minimality remain unproven.",
    )
