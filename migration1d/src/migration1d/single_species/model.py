"""Dimensionless controls for single-species strip wave-pinning."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from scipy.optimize import brentq


class ActivationProtocol(Protocol):
    def values(self, time: float, position: np.ndarray) -> np.ndarray: ...


@dataclass(frozen=True)
class NoActivation:
    def values(self, time: float, position: np.ndarray) -> np.ndarray:
        return np.zeros_like(position)


@dataclass(frozen=True)
class ActivationPulse:
    """Gaussian activation-rate pulse in relative material coordinates."""

    amplitude: float = 1.0
    center_fraction: float = 0.1
    width_fraction: float = 0.05
    start: float = 10.0
    stop: float = 30.0

    def __post_init__(self) -> None:
        if self.amplitude < 0.0:
            raise ValueError("pulse amplitude must be non-negative")
        if not 0.0 <= self.center_fraction <= 1.0:
            raise ValueError("pulse center_fraction must be in [0, 1]")
        if self.width_fraction <= 0.0 or self.stop <= self.start:
            raise ValueError("pulse width must be positive and stop must follow start")

    def values(self, time: float, position: np.ndarray) -> np.ndarray:
        if not self.start <= time < self.stop:
            return np.zeros_like(position)
        fraction = (position - position[0]) / (position[-1] - position[0])
        distance = (fraction - self.center_fraction) / self.width_fraction
        return self.amplitude * np.exp(-0.5 * distance * distance)


@dataclass(frozen=True)
class PulseTrain:
    pulses: tuple[ActivationPulse, ...]

    def values(self, time: float, position: np.ndarray) -> np.ndarray:
        values = np.zeros_like(position)
        for pulse in self.pulses:
            values += pulse.values(time, position)
        return values


@dataclass(frozen=True)
class DualPulseInitialCondition:
    """Two Gaussian active enrichments with an optional antisymmetric bias."""

    amplitude: float = 2.0
    left_center: float = 0.2
    right_center: float = 0.8
    width_fraction: float = 0.04
    asymmetry: float = 0.0

    def __post_init__(self) -> None:
        if self.amplitude < 0.0 or self.width_fraction <= 0.0:
            raise ValueError("initial pulse amplitude and width must be positive")
        if not 0.0 < self.left_center < self.right_center < 1.0:
            raise ValueError("initial pulse centers must be ordered inside the domain")
        if not 0.0 <= self.asymmetry < 1.0:
            raise ValueError("initial pulse asymmetry must be in [0, 1)")

    def profile(self, position: np.ndarray) -> np.ndarray:
        fraction = (position - position[0]) / (position[-1] - position[0])
        left = np.exp(-0.5 * ((fraction - self.left_center) / self.width_fraction) ** 2)
        right = np.exp(-0.5 * ((fraction - self.right_center) / self.width_fraction) ** 2)
        return self.amplitude * (
            (1.0 + self.asymmetry) * left + (1.0 - self.asymmetry) * right
        )


@dataclass(frozen=True)
class SingleSpeciesParameters:
    """Primary controls for the membrane--cytoplasm strip model.

    ``total_abundance`` is ``mu = mean(u + h*v)/K`` and
    ``dimensionless_length`` is ``Lambda = L*sqrt(delta/D_u)``.
    """

    total_abundance: float = 2.2683312
    dimensionless_length: float = 31.622776601683793
    cytoplasm_depth: float = 1.0
    diffusion_ratio: float = 100.0
    basal_ratio: float = 0.067
    feedback_ratio: float = 1.0
    hill_coefficient: float = 2.0
    half_saturation: float = 1.0
    inactivation_rate: float = 1.0
    domain_length: float = 10.0

    def __post_init__(self) -> None:
        positive = (
            self.total_abundance,
            self.dimensionless_length,
            self.cytoplasm_depth,
            self.diffusion_ratio,
            self.feedback_ratio,
            self.hill_coefficient,
            self.half_saturation,
            self.inactivation_rate,
            self.domain_length,
        )
        if any(value <= 0.0 for value in positive) or self.basal_ratio < 0.0:
            raise ValueError("single-species parameters must be positive")

    @property
    def membrane_diffusivity(self) -> float:
        return self.inactivation_rate * (
            self.domain_length / self.dimensionless_length
        ) ** 2

    @property
    def cytoplasm_diffusivity(self) -> float:
        return self.diffusion_ratio * self.membrane_diffusivity

    @property
    def total_density(self) -> float:
        return self.total_abundance * self.half_saturation


def reaction_flux(
    active: np.ndarray | float,
    inactive: np.ndarray | float,
    parameters: SingleSpeciesParameters,
    extra_activation: np.ndarray | float = 0.0,
) -> np.ndarray:
    """Membrane reaction flux ``J(u,v)`` in concentration/time units."""
    u = np.asarray(active, dtype=float)
    v = np.asarray(inactive, dtype=float)
    p = parameters
    scaled = np.maximum(u, 0.0) / p.half_saturation
    hill = scaled**p.hill_coefficient
    feedback = p.feedback_ratio * hill / (1.0 + hill)
    activation = p.inactivation_rate * (
        p.basal_ratio + feedback + np.asarray(extra_activation, dtype=float)
    )
    return v * activation - p.inactivation_rate * u


def _roots_on_interval(function, lower: float, upper: float, samples: int = 2049) -> list[float]:
    grid = np.linspace(lower, upper, samples)
    values = np.asarray([function(value) for value in grid])
    roots: list[float] = []
    for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if f_left == 0.0:
            root = float(left)
        elif f_left * f_right < 0.0:
            root = float(brentq(function, left, right))
        else:
            continue
        if not roots or abs(root - roots[-1]) > 1.0e-7:
            roots.append(root)
    if values[-1] == 0.0 and (not roots or abs(upper - roots[-1]) > 1.0e-7):
        roots.append(float(upper))
    return roots


def local_reaction_roots(
    inactive: float, parameters: SingleSpeciesParameters, maximum_active: float | None = None
) -> list[dict[str, float | bool]]:
    """Find and classify roots of ``J(u,v)=0`` for fixed cytoplasmic ``v``."""
    upper = maximum_active or max(5.0 * parameters.half_saturation, parameters.total_density)
    function = lambda active: float(reaction_flux(active, inactive, parameters))
    roots = _roots_on_interval(function, 0.0, upper)
    records = []
    step = max(1.0e-6, 1.0e-5 * parameters.half_saturation)
    for root in roots:
        derivative = (function(root + step) - function(max(0.0, root - step))) / (
            step if root < step else 2.0 * step
        )
        records.append({"active": root, "stable": bool(derivative < 0.0)})
    return records


def homogeneous_equilibria(parameters: SingleSpeciesParameters) -> list[dict[str, float | bool]]:
    """Equilibria on ``u + h*v = mu*K`` with homogeneous stability labels."""
    p = parameters
    total = p.total_density

    def constrained_flux(active: float) -> float:
        inactive = (total - active) / p.cytoplasm_depth
        return float(reaction_flux(active, inactive, p))

    roots = _roots_on_interval(constrained_flux, 0.0, total)
    records = []
    step = max(1.0e-7, 1.0e-5 * total)
    for active in roots:
        inactive = (total - active) / p.cytoplasm_depth
        left, right = max(0.0, active - step), min(total, active + step)
        derivative = (constrained_flux(right) - constrained_flux(left)) / (right - left)
        records.append({
            "active": active,
            "inactive": inactive,
            "stable": bool(derivative < 0.0),
        })
    return records


class StripWavePinningChemistry:
    """One active membrane species coupled to a volumetric inactive pool."""

    species = ("active_gtpase",)

    def __init__(
        self,
        parameters: SingleSpeciesParameters | None = None,
        activation: ActivationProtocol | None = None,
        initial_condition: DualPulseInitialCondition | None = None,
    ) -> None:
        self.parameters = parameters or SingleSpeciesParameters()
        self.activation = activation or NoActivation()
        self.initial_condition = initial_condition
        equilibria = [item for item in homogeneous_equilibria(self.parameters) if item["stable"]]
        if not equilibria:
            raise ValueError("parameter set has no stable homogeneous equilibrium")
        resting = min(equilibria, key=lambda item: float(item["active"]))
        self.initial_membrane = np.asarray([resting["active"]], dtype=float)
        self.initial_cytoplasm = np.asarray([resting["inactive"]], dtype=float)
        self.membrane_diffusion = np.asarray(
            [self.parameters.membrane_diffusivity], dtype=float
        )
        self.cytoplasm_diffusion = np.asarray(
            [self.parameters.cytoplasm_diffusivity], dtype=float
        )

    def initial_profiles(self, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        active = np.full((1, len(position)), self.initial_membrane[0])
        inactive = np.full((1, len(position)), self.initial_cytoplasm[0])
        if self.initial_condition is not None:
            transfer = self.initial_condition.profile(position)
            maximum = 0.9 * self.parameters.cytoplasm_depth * inactive[0]
            transfer = np.minimum(transfer, maximum)
            active[0] += transfer
            inactive[0] -= transfer / self.parameters.cytoplasm_depth
        return active, inactive

    def reaction(self, time, position, membrane, cytoplasm, strain, stress):
        _ = strain, stress
        flux = reaction_flux(
            membrane[0], cytoplasm[0], self.parameters,
            self.activation.values(time, position),
        )
        membrane_rate = flux[np.newaxis, :]
        cytoplasm_rate = -membrane_rate / self.parameters.cytoplasm_depth
        return membrane_rate, cytoplasm_rate

    def protrusive_activity(self, membrane, cytoplasm):
        """Use the sole active GTPase as the model's Rac-equivalent activity."""

        _ = cytoplasm
        return np.maximum(membrane[0], 0.0)

    def tension_deactivation(self, membrane, cytoplasm, added_rate):
        """Transfer tension-inactivated GTPase without changing ``u + h*v``."""

        _ = cytoplasm
        loss = np.asarray(added_rate, dtype=float) * np.maximum(membrane[0], 0.0)
        membrane_rate = -loss[np.newaxis, :]
        cytoplasm_rate = -membrane_rate / self.parameters.cytoplasm_depth
        return membrane_rate, cytoplasm_rate
