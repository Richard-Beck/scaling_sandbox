"""Reusable development adapter for the Spring circular extension analogue.

The dynamics are inherited unchanged from the 17 October development closure.
They retain native local chemistry, diffusivities, Rac force switch, global
extension feedback and dilution. The mean-force circumference law and uniform
material deformation are explicit geometric assumptions, not native 2D
mechanics. This file adds an exact dense Jacobian and periodic assay/readout
helpers; it does not replace the active interval model or retained evidence.
"""
from __future__ import annotations

from dataclasses import dataclass
import os

import numpy as np

from spring_ring_closure_20261017 import Closure as _OriginalClosure


class Closure(_OriginalClosure):
    """Identical Spring ring dynamics with the complete analytic Jacobian."""

    def jacobian(self, state, cue, cue_circumference_derivative=0.0):
        """Return d(rhs)/d(state), including all global dependencies.

        ``cue`` is an additive activation profile. It may depend on the current
        circumference, with derivative supplied in the third argument. The
        default derivative zero is correct only for a state-independent cue
        profile. In particular, a fixed physical width in material coordinates
        requires the derivative returned by :func:`physical_raised_cosine`.
        No dependence of the imposed cue on chemical concentrations is assumed.
        """
        u, v, circumference = self.fields(state)
        n = self.cells
        signal = np.broadcast_to(np.asarray(cue, dtype=float), (n,))
        signal_dC = np.broadcast_to(
            np.asarray(cue_circumference_derivative, dtype=float), (n,)
        )
        if not np.isfinite(signal).all() or not np.isfinite(signal_dC).all():
            raise ValueError("Cue and circumference derivative must be finite")

        extension = circumference - self.resting_circumference
        # These derivative constants match the fixed native laws used by the
        # inherited RHS: F(u)=0.001/(1+exp(-20*(u-1))) and Hill exponent six.
        from spring_ring_closure_20261017 import native
        force = native.auth["F"](u)
        speed = 2.0 * (np.mean(force) - 0.01 * extension)
        speed_du = 40.0 / n * force * (1.0 - force / 0.001)
        speed_dC = -0.02
        activation = 4.0 + signal + 5.0 * u**6 / (1.0 + u**6)
        reaction_du = 30.0 * u**5 / (1.0 + u**6)**2 * v
        reaction_du -= 3.0 + 80.0 * extension
        reaction_dC = signal_dC * v - 80.0 * u

        index = np.arange(n)
        laplacian = np.zeros((n, n))
        laplacian[index, index] = -2.0
        laplacian[index, (index - 1) % n] = 1.0
        laplacian[index, (index + 1) % n] = 1.0
        laplacian *= (n / circumference)**2
        diffusion_u = 0.01 * (np.roll(u, 1) - 2.0*u + np.roll(u, -1))
        diffusion_u *= (n / circumference)**2
        diffusion_v = 10.0 * (np.roll(v, 1) - 2.0*v + np.roll(v, -1))
        diffusion_v *= (n / circumference)**2

        jac = np.zeros((2*n + 1, 2*n + 1))
        # Mean force changes the dilution of every material cell, including
        # inactive cells far from the perturbed active concentration.
        jac[:n, :n] = -(u / circumference)[:, None] * speed_du[None, :]
        jac[:n, :n] += 0.01 * laplacian
        jac[index, index] += reaction_du - speed / circumference
        jac[index, n + index] = activation
        jac[n:-1, :n] = -(v / circumference)[:, None] * speed_du[None, :]
        jac[n + index, index] -= reaction_du
        jac[n:-1, n:-1] = 10.0 * laplacian
        jac[n + index, n + index] -= activation + speed / circumference

        dilution_dC = -speed_dC / circumference + speed / circumference**2
        jac[:n, -1] = reaction_dC - 2.0*diffusion_u/circumference + u*dilution_dC
        jac[n:-1, -1] = -reaction_dC - 2.0*diffusion_v/circumference + v*dilution_dC
        jac[-1, :n] = speed_du
        jac[-1, -1] = speed_dC
        return jac


def physical_raised_cosine(closure, state, center_fraction, peak, total_support):
    """Return a periodic cue and its derivative with respect to circumference.

    The centre is fixed in material coordinates; ``total_support`` is the full
    instantaneous physical arc width. The peak is additive native activation.
    This is a circular cue, not a dose match to the interval endpoint half-patch.
    At the support boundary both the cue and its first derivative are zero.
    """
    _, _, circumference = closure.fields(state)
    if not np.isfinite(center_fraction) or not np.isfinite(peak):
        raise ValueError("Cue centre and peak must be finite")
    if not np.isfinite(total_support) or not 0.0 < total_support <= circumference:
        raise ValueError("Total cue support must be positive and at most the circumference")
    material_distance = np.abs(
        (closure.material_positions - float(center_fraction) + 0.5) % 1.0 - 0.5
    )
    distance = material_distance * circumference
    inside = distance < total_support / 2.0
    phase = 2.0 * np.pi * distance[inside] / total_support
    signal = np.zeros(closure.cells)
    derivative = np.zeros(closure.cells)
    signal[inside] = float(peak) * 0.5 * (1.0 + np.cos(phase))
    derivative[inside] = (
        -float(peak) * np.pi / total_support * np.sin(phase) * material_distance[inside]
    )
    return signal, derivative


def _active_vector_history(closure, active):
    """First circular moment for raw histories, including a zero active pool."""
    active = np.asarray(active, dtype=float)
    total = np.sum(active, axis=-1)
    numerator = active @ np.exp(2j*np.pi*closure.material_positions)
    vector = np.zeros_like(numerator, dtype=complex)
    np.divide(numerator, total, out=vector, where=total != 0.0)
    return vector


def polarity_vector(closure, state, direction_floor=1e-10):
    """Normalized first circular moment of active Rac, without a velocity law.

    Uniform cell weights reflect the uniform material grid. Circumference
    cancels between numerator and denominator. Report direction as unavailable
    when the normalized vector magnitude is at or below ``direction_floor``;
    this avoids assigning a meaningful angle to a uniform or cancelling field.
    Magnitude alone does not count local peaks or establish whole-field polarity.
    """
    u, _, _ = closure.fields(state)
    if direction_floor < 0 or not np.isfinite(direction_floor):
        raise ValueError("Direction floor must be finite and nonnegative")
    if np.min(u) < 0.0:
        raise ValueError("The active-weighted polarity vector requires nonnegative activity")
    total = float(np.sum(u))
    if total == 0.0:
        return dict(x=0.0, y=0.0, magnitude=0.0, direction_radians=None)
    vector = complex(_active_vector_history(closure, u))
    magnitude = float(np.abs(vector))
    return dict(x=float(vector.real), y=float(vector.imag), magnitude=magnitude,
                direction_radians=float(np.angle(vector)) if magnitude > direction_floor else None)


@dataclass(frozen=True)
class CuePhase:
    """One constant-amplitude circular cue phase, including zero-cue recovery.

    ``center_fraction`` stays at a material location. ``total_support`` is an
    instantaneous physical arc width; it is ignored when ``peak`` is zero.
    Every phase boundary is an exact solver restart, including withdrawal.
    Durations and output intervals are uncalibrated native Spring time units.
    """

    duration: float
    peak: float = 0.0
    center_fraction: float = 0.0
    total_support: float = 0.2
    output_interval: float = 5.0

    def __post_init__(self):
        values = (self.duration, self.peak, self.center_fraction,
                  self.total_support, self.output_interval)
        if not np.isfinite(values).all():
            raise ValueError("Phase values must be finite")
        if self.duration <= 0 or self.total_support <= 0 or self.output_interval <= 0:
            raise ValueError("Phase duration, support and output interval must be positive")


def run_protocol(closure, phases, *, initial_state=None, rtol=1e-9,
                 atol=1e-11, max_step=10.0):
    """Run a restartable cue protocol via SLURM and return full ring histories.

    This reusable entry point uses the same inherited equations as the retained
    feasibility pilot, with the complete analytic Jacobian. Native/frozen
    sources are never changed. Submit this call in a SLURM job; freezing sources
    and recording campaign provenance belong to the caller. Pass any returned
    final ``states[-1]`` as ``initial_state`` to continue that exact state.
    """
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Spring ring trajectory integration must run via SLURM")
    from scipy.integrate import solve_ivp

    phases = tuple(phases)
    if not phases or not all(isinstance(phase, CuePhase) for phase in phases):
        raise ValueError("Provide at least one CuePhase")
    if not np.isfinite((rtol, atol, max_step)).all() or min(rtol, atol, max_step) <= 0:
        raise ValueError("Solver tolerances and maximum step must be finite and positive")
    state = closure.initial() if initial_state is None else np.array(initial_state, dtype=float, copy=True)
    closure.fields(state)
    if np.min(state[:-1]) < 0 or closure.amount(state) <= 0:
        raise ValueError("Initial concentrations must be nonnegative with positive physical amount")
    reference_amount = closure.amount(state)
    times, states, history = [0.0], [state.copy()], []
    elapsed = 0.0

    for phase in phases:
        def cue_at(y):
            if phase.peak == 0.0:
                return 0.0, 0.0
            return physical_raised_cosine(closure, y, phase.center_fraction,
                                         phase.peak, phase.total_support)

        def function(t, y):
            return closure.rhs(y, cue_at(y)[0])

        def jacobian(t, y):
            return closure.jacobian(y, *cue_at(y))

        sample = np.arange(0.0, phase.duration, phase.output_interval)
        sample = np.r_[sample, phase.duration]
        solution = solve_ivp(function, (0.0, phase.duration), state,
                             method="BDF", t_eval=sample, jac=jacobian,
                             rtol=rtol, atol=atol, max_step=max_step)
        if not solution.success:
            raise RuntimeError("Spring ring integration failed: " + solution.message)
        if not np.isfinite(solution.y).all():
            raise RuntimeError("Spring ring integration returned a nonfinite state")
        times.extend(elapsed + solution.t[1:])
        states.extend(solution.y.T[1:])
        history.append(dict(start=elapsed, stop=elapsed + phase.duration,
            peak=phase.peak, center_fraction=phase.center_fraction,
            total_physical_support=phase.total_support if phase.peak else None,
            output_interval=phase.output_interval, nfev=solution.nfev,
            njev=solution.njev, nlu=solution.nlu))
        state = solution.y[:, -1].copy()
        elapsed += phase.duration

    states = np.asarray(states)
    circumference = states[:, -1]
    active = states[:, :closure.cells]
    inactive = states[:, closure.cells:-1]
    amount = circumference * np.mean(active + inactive, axis=1)
    if float(np.min(states[:, :-1])) < -5e-10:
        raise RuntimeError("Spring ring saved concentration is negative beyond roundoff")
    amount_error = float(np.max(np.abs(amount/reference_amount - 1.0)))
    if amount_error > 1e-7:
        raise RuntimeError("Spring ring physical-amount error exceeds 1e-7")
    # Tiny negative roundoff must never be hidden by clipping. Compute vector
    # directly here while retaining the raw fields and minimum concentration.
    vector = _active_vector_history(closure, active)
    magnitude = np.abs(vector)
    direction = np.where(magnitude > 1e-10, np.angle(vector), np.nan)
    return dict(model="Spring circular global-extension analogue", time=np.asarray(times),
        states=states, active=active, inactive=inactive, circumference=circumference,
        material_positions=closure.material_positions,
        physical_positions=circumference[:, None]*closure.material_positions[None, :],
        physical_amount=amount, max_saved_amount_relative_error=amount_error,
        minimum_saved_concentration=float(np.min(states[:, :-1])),
        polarity_x=vector.real, polarity_y=vector.imag, polarity_magnitude=magnitude,
        polarity_direction_radians=direction, phases=history,
        native_units="uncalibrated author-model length and time", native_2d_equivalence=False)
