"""Frozen native ring adapter and bounded accepted-step archives for Figure 3."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import signal
import time

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"
import numpy as np
from scipy.integrate import BDF
from scipy.optimize import brentq
import ring_distance_prototype_20261017 as chemical
from spring_ring_model_20261019 import Closure, physical_raised_cosine

ROOT = Path(__file__).resolve().parents[2]
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring")
SIZES = {k: (30., 37.8, 47.6) for k in MODELS[:-1]}
SIZES["spring"] = (2.5, 3.15, 3.97)
PREP_PEAKS = dict(wave_pinning=.1, goryachev=.1, holmes_model3=1., spring=4.)
PEAKS = dict(wave_pinning=.1, goryachev=6.4, holmes_model3=1., spring=4.)
SCALES = dict(wave_pinning=[2.2683, 2.2683], goryachev=[2.75033, 2.75033],
              holmes_model3=[7.5, 7.5, 3.1, 3.1])
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def require_slurm():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Full native integration requires SLURM")
    signal.signal(signal.SIGUSR1, lambda *_: globals().__setitem__("STOP", True))
    signal.signal(signal.SIGTERM, lambda *_: globals().__setitem__("STOP", True))


class Adapter:
    def __init__(self, model, size, cells=None):
        self.key, self.size = model, float(size)
        self.spring = model == "spring"
        self.cells = int(round(size*(100 if self.spring else 20))) if cells is None else int(cells)
        self.model = Closure(self.size, self.cells) if self.spring else chemical.native_model(model)
        self.grid = None if self.spring else chemical.RingGrid(self.size, self.cells)
        self.q = (np.arange(self.cells)+.5)/self.cells
        self.support = .2 if self.spring else 3.
        self.scales = np.r_[np.full(2*self.cells, 2.), self.size] if self.spring else np.repeat(np.asarray(SCALES[model]), self.cells)

    def initial(self):
        return self.model.initial() if self.spring else chemical.native_rest(self.key, self.grid).ravel()

    def circumference(self, state):
        return float(state[-1]) if self.spring else self.size

    def primary(self, state):
        return np.asarray(state)[:self.cells]

    def amounts(self, state):
        if self.spring:
            return np.asarray([self.model.amount(state)])
        return chemical.pool_amounts(self.model, self.grid, np.asarray(state).reshape(-1, self.cells))

    def concentration_minimum(self, state):
        return float(np.min(np.asarray(state)[:-1] if self.spring else state))

    def cue(self, state, peak, center_fraction):
        if peak == 0:
            return np.zeros(self.cells), np.zeros(self.cells)
        if self.spring:
            return physical_raised_cosine(self.model, state, center_fraction, peak, self.support)
        distance = np.abs((self.q-center_fraction+.5) % 1-.5)*self.size
        cue = np.where(distance < self.support/2,
                       peak*.5*(1+np.cos(2*np.pi*distance/self.support)), 0.)
        return cue, np.zeros(self.cells)

    def rhs(self, state, peak=0., center_fraction=0.):
        cue, _ = self.cue(state, peak, center_fraction)
        if self.spring:
            return self.model.rhs(state, cue)
        return chemical.native_rhs(self.model, self.grid, np.asarray(state).reshape(-1, self.cells), cue).ravel()

    def jacobian_options(self, peak, center_fraction):
        if self.spring:
            return dict(jac=lambda t, y: self.model.jacobian(y, *self.cue(y, peak, center_fraction)))
        return dict(jac_sparsity=chemical.periodic_sparsity(len(self.model.state_names), self.cells))


def samples(end, dense=False):
    if dense:
        return np.unique(np.r_[np.arange(0, end, min(.1, end/100)), end])
    return np.unique(np.r_[0., np.arange(.1, min(10., end), .1),
        np.arange(10., min(100., end), 1.), np.arange(100., min(1000., end), 5.),
        np.arange(1000., end, 20.), end])


def integrate(adapter, initial, duration, peak=0., center_fraction=0.,
              deadline=None, event=None, interval=None, solver_options=None):
    """Complete fields at saved samples, exact endpoints and event root.

    Events detect a negative-to-nonnegative crossing on accepted BDF steps and
    solve the dense interpolant to the crossing. A supplied event must be
    continuous. All invariant checks include accepted steps, not only output.
    Timeout and numerical failure preserve a partial archive for diagnosis;
    neither is scientific right censoring at the requested horizon.
    """
    if duration <= 0:
        raise ValueError("Need positive integration duration")
    wanted = samples(duration) if interval is None else np.unique(np.r_[np.arange(0, duration, interval), duration])
    tolerances = dict(rtol=1e-9, atol=1e-11, max_step=10.)
    tolerances.update(solver_options or {})
    solver = BDF(lambda t, y: adapter.rhs(y, peak, center_fraction), 0.,
                 np.array(initial, copy=True), float(duration),
                 **tolerances, **adapter.jacobian_options(peak, center_fraction))
    tt, yy, index = [0.], [np.array(initial, copy=True)], 1
    reference = adapter.amounts(initial)
    minimum, mass_error = adapter.concentration_minimum(initial), 0.
    status, message, event_time = "complete", "", None
    previous_event = float(event(initial)) if event else None
    if previous_event is not None and previous_event >= 0:
        return dict(time=np.asarray(tt), states=np.asarray(yy), status="event_at_start", message="",
                    event_time=0., max_amount_relative_error=0., minimum_concentration=minimum,
                    nfev=0, njev=0, nlu=0, requested_duration=float(duration))
    while solver.status == "running":
        if STOP or (deadline is not None and time.monotonic() >= deadline):
            status, message = "timeout", "SLURM advance signal or worker soft deadline"
            break
        previous_t = solver.t
        message = solver.step() or ""
        if solver.status == "failed":
            status = "numerical_error"
            break
        minimum = min(minimum, adapter.concentration_minimum(solver.y))
        mass_error = max(mass_error, float(np.max(np.abs(adapter.amounts(solver.y)/reference-1))))
        if minimum < -5e-10 or mass_error > 1e-7 or not np.isfinite(solver.y).all():
            status, message = "numerical_error", "Positivity, finiteness or native conserved-amount gate failed"
            break
        dense = solver.dense_output()
        cutoff = solver.t
        if event is not None:
            value = float(event(solver.y))
            if previous_event < 0 <= value:
                event_time = float(brentq(lambda t: event(dense(t)), previous_t, solver.t, xtol=1e-9, rtol=1e-12))
                cutoff, status = event_time, "event"
            previous_event = value
        while index < len(wanted) and wanted[index] <= cutoff:
            tt.append(float(wanted[index])); yy.append(dense(wanted[index]).copy()); index += 1
        if event_time is not None:
            if tt[-1] != event_time:
                tt.append(event_time); yy.append(dense(event_time).copy())
            break
    if event_time is None and tt[-1] != solver.t:
        tt.append(float(solver.t)); yy.append(solver.y.copy())
    return dict(time=np.asarray(tt), states=np.asarray(yy), status=status, message=message,
                event_time=event_time, max_amount_relative_error=mass_error,
                minimum_concentration=minimum, nfev=int(solver.nfev), njev=int(solver.njev),
                nlu=int(solver.nlu), requested_duration=float(duration))


def save_phase(path, result, adapter, peak=0., center_fraction=0.):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, time=result["time"], states=result["states"],
        material_positions=adapter.q, model=adapter.key, size_parameter=adapter.size,
        cells=adapter.cells, cue_peak=peak, cue_center_fraction=center_fraction,
        cue_physical_support=adapter.support,
        state_layout="u[0:N], v[N:2N], circumference[-1]" if adapter.spring else "native species-major flat fields",
        state_names=np.asarray(["active", "inactive", "circumference"] if adapter.spring else adapter.model.state_names))
    return dict(path=str(path), sha256=sha(path),
                **{k: v for k, v in result.items() if k not in ("time", "states")},
                stored_full_native_field_samples=len(result["time"]))


def window_metric(times, states, scales, window=1000., tolerance=1e-4):
    keep = times >= max(0., times[-1]-window)
    t, y = times[keep], states[keep]
    if len(t) < 2:
        return dict(available=False, passed=False)
    variation = float(np.max(np.ptp(y, axis=0)/scales))
    projected = float(np.max(np.abs(y[-1]-y[-2])/scales)/(t[-1]-t[-2])*window)
    return dict(available=True, passed=max(variation, projected) <= tolerance,
                scaled_range=variation, projected_last_saved_rate=projected,
                window=window, tolerance=tolerance, samples=len(t))


def finite_tail(adapter, times, states):
    window = 100. if adapter.spring else 1000.
    records = []
    for end in (times[-1]-window, times[-1]):
        keep = (times >= end-window) & (times <= end)
        records.append(window_metric(times[keep], states[keep], adapter.scales, window=window)
                       if keep.sum() >= 2 else dict(available=False, passed=False))
    return dict(windows=records, both_late_windows_steady=all(r["passed"] for r in records),
                interpretation="Finite full native-field gate, independent of winner or scalar takeover")
