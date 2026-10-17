"""Canonical Holmes parameter mapping and homogeneous stability analysis."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from scipy.integrate import solve_ivp

from migration1d.holmes import HolmesParameters


class GaussianRacPulse:
    """Localized Rac-GEF pulse in relative material coordinates."""

    def __init__(
        self,
        amplitude: float = 1.0,
        center_fraction: float = 0.1,
        width_fraction: float = 0.05,
        start: float = 10.0,
        stop: float = 30.0,
    ) -> None:
        if amplitude < 0.0 or width_fraction <= 0.0 or stop <= start:
            raise ValueError("pulse amplitude/width must be positive and stop must follow start")
        if not 0.0 <= center_fraction <= 1.0:
            raise ValueError("pulse center_fraction must be in [0, 1]")
        self.amplitude = amplitude
        self.center_fraction = center_fraction
        self.width_fraction = width_fraction
        self.start = start
        self.stop = stop

    def __call__(self, time: float, position: np.ndarray) -> np.ndarray:
        if not self.start <= time < self.stop:
            return np.zeros_like(position)
        fraction = (position - position[0]) / (position[-1] - position[0])
        distance = (fraction - self.center_fraction) / self.width_fraction
        return self.amplitude * np.exp(-0.5 * distance * distance)


class RacPulseTrain:
    """Sum of localized Rac-GEF pulses, used by the reversal protocol."""

    def __init__(self, pulses: tuple[GaussianRacPulse, ...]) -> None:
        if not pulses:
            raise ValueError("a pulse train requires at least one pulse")
        self.pulses = pulses

    def __call__(self, time: float, position: np.ndarray) -> np.ndarray:
        result = np.zeros_like(position)
        for pulse in self.pulses:
            result += pulse(time, position)
        return result


def canonical_holmes_parameters(
    *, kappa_r: float = 0.2 / 7.5, hill_exponent: float = 3.0
) -> HolmesParameters:
    """Return the GTPase-only Holmes baseline stated in parameter_mapping.

    ``kappa_r = I_R1 / (d_R R_tot)`` and ``d_R=1 / s`` in this reduction.
    The legacy migration example uses feedback 0.55 and a continuously applied
    signal; the standardized screen instead uses the mapped 0.65 / s feedback
    and supplies signal only through the explicit pulse.
    """
    if kappa_r < 0.0:
        raise ValueError("kappa_r must be non-negative")
    base = HolmesParameters()
    return replace(
        base,
        rac_basal_activation=kappa_r * base.rac_total,
        hill_exponent=hill_exponent,
        cdc42_to_rac_feedback=0.65,
        signal_baseline=0.0,
        active_diffusivity=0.1,
        inactive_diffusivity=50.0,
    )


def active_reaction(active: np.ndarray, inactive: np.ndarray, p: HolmesParameters) -> np.ndarray:
    """Homogeneous active-field reaction rates without an external cue."""
    cdc42, rac, rho = active
    cdc42_i, rac_i, rho_i = inactive
    accessibility = p.inactive_accessibility
    return np.asarray((
        p.cdc42_activation
        / (1.0 + (rho / p.rho_inhibits_cdc42_half_level) ** p.hill_exponent)
        * accessibility * cdc42_i / p.cdc42_total - cdc42,
        (p.rac_basal_activation + p.cdc42_to_rac_feedback * cdc42)
        * accessibility * rac_i / p.rac_total - rac,
        p.rho_activation
        / (1.0 + (rac / p.rac_inhibits_rho_half_level) ** p.hill_exponent)
        * accessibility * rho_i / p.rho_total - rho,
    ))


def homogeneous_equilibrium(
    parameters: HolmesParameters, *, duration: float = 1000.0
) -> np.ndarray:
    """Relax the conserved homogeneous circuit to its low-activity rest state."""
    totals = parameters.totals
    initial = np.minimum(np.asarray((0.05, 0.05, 0.05)), 0.05 * totals)

    def rhs(_time: float, active: np.ndarray) -> np.ndarray:
        return active_reaction(active, totals - active, parameters)

    solution = solve_ivp(
        rhs, (0.0, duration), initial, method="BDF",
        rtol=1.0e-10, atol=1.0e-12,
    )
    if not solution.success:
        raise RuntimeError(f"homogeneous Holmes relaxation failed: {solution.message}")
    active = solution.y[:, -1]
    residual = np.max(np.abs(rhs(duration, active)))
    if residual > 1.0e-7:
        raise RuntimeError(f"homogeneous Holmes relaxation residual is {residual:g}")
    return active


def _jacobian(function, state: np.ndarray) -> np.ndarray:
    state = np.asarray(state, dtype=float)
    result = np.empty((len(state), len(state)))
    for column, value in enumerate(state):
        step = 1.0e-6 * max(1.0, abs(float(value)))
        left, right = state.copy(), state.copy()
        left[column] -= step
        right[column] += step
        result[:, column] = (function(right) - function(left)) / (2.0 * step)
    return result


def homogeneous_stability(
    parameters: HolmesParameters,
    active: np.ndarray,
    domain_length: float,
    *, maximum_mode: int = 64,
) -> dict[str, float | int]:
    """Return constrained homogeneous and nonzero spatial growth rates.

    The zero-mode calculation uses the three conserved constraints
    ``inactive = total - active``. Nonzero Neumann modes use the full six-field
    Jacobian and add ``-D (m*pi/L)^2`` to each diagonal entry.
    """
    totals = parameters.totals
    local_jacobian = _jacobian(
        lambda values: active_reaction(values, totals - values, parameters), active
    )
    local_growth = float(np.max(np.real(np.linalg.eigvals(local_jacobian))))

    state = np.concatenate((active, totals - active))

    def full_reaction(values: np.ndarray) -> np.ndarray:
        rates = active_reaction(values[:3], values[3:], parameters)
        return np.concatenate((rates, -rates))

    reaction_jacobian = _jacobian(full_reaction, state)
    diffusion = np.asarray(
        [parameters.active_diffusivity] * 3 + [parameters.inactive_diffusivity] * 3
    )
    spatial_growth = -np.inf
    fastest_mode = 1
    for mode in range(1, maximum_mode + 1):
        wave_number = mode * np.pi / domain_length
        matrix = reaction_jacobian - np.diag(diffusion * wave_number**2)
        growth = float(np.max(np.real(np.linalg.eigvals(matrix))))
        if growth > spatial_growth:
            spatial_growth, fastest_mode = growth, mode
    return {
        "homogeneous_growth_rate": local_growth,
        "maximum_spatial_growth_rate": spatial_growth,
        "fastest_spatial_mode": fastest_mode,
    }
