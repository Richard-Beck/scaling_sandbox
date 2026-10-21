"""Topology-aware whole-region cue-withdrawal events in native BDF evolution.

Native equations, tolerances and accepted-step conservation gates are the
retained Figure 3 implementation. A discontinuity in the *measurement* when
a second resolved peak appears is never presented as a continuous ratio root.
"""
from __future__ import annotations

import time
import numpy as np
from scipy.integrate import BDF
from scipy.optimize import brentq

import signature_figure3_native_20261020 as native
import signature_figure3_regions_20261020 as regions


class RegionSeedEvent:
    def __init__(self, adapter, placement, ratio, basal_excess):
        self.adapter, self.placement = adapter, placement
        self.ratio, self.basal_excess = float(ratio), float(basal_excess)
        self.old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
        self.cue = placement["cue_center_fraction"]
        self.ever_two = False
        self.topology_transitions = []

    def readout(self, state):
        return regions.readout(self.adapter.primary(state), self.adapter.circumference(state),
                               self.old, self.cue, detailed=False)

    def ratio_value(self, spatial):
        if not spatial["two_region_available"]:
            raise ValueError("Ratio root is unavailable without two resolved regions")
        return (spatial["challenger_excess"]-self.ratio*spatial["incumbent_excess"])/self.basal_excess

    def locate(self, first_t, first, final_t, final, dense):
        """First eligible event in one accepted BDF interval, or no event.

        Dense subdivisions are required when topology changes or either ratio
        or collapse can cross. Brent's method only sees a same-topology pair.
        A topology/birth boundary uses a declared binary-predicate bracket,
        not a fabricated root of a discontinuous ratio.
        """
        changed = first["state_class"] != final["state_class"]
        two_first, two_final = first["two_region_available"], final["two_region_available"]
        needs_check = changed
        if two_first and two_final:
            needs_check |= self.ratio_value(first) < 0 <= self.ratio_value(final)
            needs_check |= first["incumbent_excess"] > .1*self.basal_excess >= final["incumbent_excess"]
        same_single_old = first["state_class"] == final["state_class"] == "single_region_incumbent_associated"
        if same_single_old:
            needs_check |= first["incumbent_excess"] > .1*self.basal_excess >= final["incumbent_excess"]
        if not needs_check:
            self.ever_two |= two_final
            return None
        ts = np.linspace(first_t, final_t, 17)
        values = [first]+[self.readout(dense(t)) for t in ts[1:-1]]+[final]
        ever_two = self.ever_two or two_first
        for left, right, before, after in zip(ts[:-1], ts[1:], values[:-1], values[1:]):
            if before["state_class"] != after["state_class"]:
                self.topology_transitions.append(dict(lower=float(left), upper=float(right),
                    previous=before["state_class"], current=after["state_class"]))
            before_two, after_two = before["two_region_available"], after["two_region_available"]
            candidates = []
            same_single_old = before["state_class"] == after["state_class"] == "single_region_incumbent_associated"
            if (before_two and after_two) or same_single_old:
                functions = [("continuous_incumbent_90_percent_excess_loss", lambda v: .1-v["incumbent_excess"]/self.basal_excess)]
                if before_two and after_two:
                    functions.insert(0, ("continuous_two_region_ratio", self.ratio_value))
                for name, function in functions:
                    if function(before) < 0 <= function(after):
                        def root_function(t):
                            spatial = self.readout(dense(t))
                            if not spatial["two_region_available"] and not (same_single_old and spatial["state_class"] == "single_region_incumbent_associated"):
                                raise ValueError("Topology changed inside the proposed continuous root")
                            return function(spatial)
                        try:
                            root = float(brentq(root_function, left, right, xtol=1e-9, rtol=1e-12))
                        except ValueError:
                            candidates.append(dict(time=float(right),
                                kind="target_attainment_across_unresolved_topology" if name == "continuous_two_region_ratio" else "incumbent_collapse_across_unresolved_topology",
                                continuous_root=False, measurement_boundary_bracket=[float(left), float(right)],
                                note="Dense ratio/collapse root lost valid topology; preserve alternate outcome without asserting equality"))
                            continue
                        at_root = self.readout(dense(root))
                        candidates.append(dict(time=root, kind=name, continuous_root=True,
                            root_value=float(function(at_root)), accepted_step_lower=float(first_t),
                            accepted_step_upper=float(final_t), dense_bracket=[float(left), float(right)]))
            if after_two and not before_two and self.ratio_value(after) >= 0:
                def attained(t):
                    spatial = self.readout(dense(t))
                    return spatial["two_region_available"] and self.ratio_value(spatial) >= 0
                low, high = float(left), float(right)
                while high-low > 1e-9:
                    middle = (low+high)/2
                    if attained(middle): high = middle
                    else: low = middle
                candidates.append(dict(time=high, kind="target_attainment_at_resolved_region_birth",
                    continuous_root=False, measurement_boundary_bracket=[low, high],
                    note="First available two-region ratio already exceeds target; no continuous ratio equality asserted"))
            cue_single = after["state_class"] == "single_region_cue_associated"
            low_contrast = after["state_class"] == "homogeneous_or_unresolved_low_contrast"
            if (cue_single or low_contrast) and before["state_class"] != after["state_class"]:
                def alternate(t):
                    value = self.readout(dense(t))
                    return value["state_class"] == after["state_class"]
                low, high = float(left), float(right)
                while high-low > 1e-9:
                    middle = (low+high)/2
                    if alternate(middle): high = middle
                    else: low = middle
                candidates.append(dict(time=high,
                    kind="takeover_during_cue_single_region" if cue_single else "activity_collapse_during_cue",
                    continuous_root=False, measurement_boundary_bracket=[low, high],
                    note="Topology outcome, not a two-region ratio root"))
            if candidates:
                candidates.sort(key=lambda item: item["time"])
                result = candidates[0]
                result["coincident_event_kinds"] = [v["kind"] for v in candidates if abs(v["time"]-result["time"]) <= 1e-7]
                self.ever_two = ever_two or before_two or after_two
                return result
            ever_two |= before_two or after_two
        self.ever_two = ever_two
        return None


def integrate_seed(adapter, initial, duration, peak, center_fraction, detector,
                   deadline=None, solver_options=None):
    if duration <= 0:
        raise ValueError("Need positive seeding duration")
    wanted = native.samples(duration)
    tolerances = dict(rtol=1e-9, atol=1e-11, max_step=10.)
    tolerances.update(solver_options or {})
    solver = BDF(lambda t, y: adapter.rhs(y, peak, center_fraction), 0.,
        np.array(initial, copy=True), float(duration), **tolerances,
        **adapter.jacobian_options(peak, center_fraction))
    tt, yy, index = [0.], [np.array(initial, copy=True)], 1
    reference = adapter.amounts(initial)
    minimum, mass_error = adapter.concentration_minimum(initial), 0.
    status, message, event_record = "complete", "", None
    previous = detector.readout(initial)
    while solver.status == "running":
        if native.STOP or (deadline is not None and time.monotonic() >= deadline):
            status, message = "timeout", "SLURM advance signal or worker soft deadline"
            break
        previous_t = solver.t
        message = solver.step() or ""
        if solver.status == "failed":
            status = "numerical_error"; break
        minimum = min(minimum, adapter.concentration_minimum(solver.y))
        mass_error = max(mass_error, float(np.max(np.abs(adapter.amounts(solver.y)/reference-1))))
        if minimum < -5e-10 or mass_error > 1e-7 or not np.isfinite(solver.y).all():
            status, message = "numerical_error", "Positivity, finiteness or native amount gate failed"; break
        dense = solver.dense_output()
        current = detector.readout(solver.y)
        event_record = detector.locate(previous_t, previous, solver.t, current, dense)
        cutoff = event_record["time"] if event_record is not None else solver.t
        while index < len(wanted) and wanted[index] <= cutoff:
            tt.append(float(wanted[index])); yy.append(dense(wanted[index]).copy()); index += 1
        if event_record is not None:
            status = "event"
            if tt[-1] != cutoff:
                tt.append(float(cutoff)); yy.append(dense(cutoff).copy())
            minimum = min(minimum, adapter.concentration_minimum(yy[-1]))
            mass_error = max(mass_error, float(np.max(np.abs(adapter.amounts(yy[-1])/reference-1))))
            assert minimum >= -5e-10 and mass_error <= 1e-7
            break
        previous = current
    if event_record is None and tt[-1] != solver.t:
        tt.append(float(solver.t)); yy.append(solver.y.copy())
    return dict(time=np.asarray(tt), states=np.asarray(yy), status=status, message=message,
        event_time=None if event_record is None else float(event_record["time"]), event_record=event_record,
        topology_transitions=detector.topology_transitions, max_amount_relative_error=mass_error,
        minimum_concentration=minimum, nfev=int(solver.nfev), njev=int(solver.njev),
        nlu=int(solver.nlu), requested_duration=float(duration))
