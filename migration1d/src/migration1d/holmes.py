"""Reduced six-field Holmes Cdc42--Rac--Rho polarity circuit.

Each GTPase has a slowly diffusing active membrane field and a rapidly
diffusing combined inactive field.  The reaction expressions and default
parameters match the reduced 1D Holmes reference in ``wgd_v1.polarity_models``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chemistry import Cue, NoCue


@dataclass(frozen=True)
class HolmesParameters:
    cdc42_total: float = 2.4
    rac_total: float = 7.5
    rho_total: float = 3.1
    cdc42_activation: float = 2.95
    rac_basal_activation: float = 0.2
    rho_activation: float = 6.6
    rho_inhibits_cdc42_half_level: float = 1.25
    rac_inhibits_rho_half_level: float = 1.0
    hill_exponent: float = 3.0
    cdc42_to_rac_feedback: float = 0.55
    signal_baseline: float = 0.3
    inactive_accessibility: float = 1.0
    active_diffusivity: float = 0.1
    inactive_diffusivity: float = 50.0
    # Prestimulus homogeneous equilibrium (external Rac signal equal to zero).
    # As in the repository's reduced Holmes reference, signal_baseline is
    # switched on when the spatial simulation begins.
    initial_active: tuple[float, float, float] = (
        0.43933939943501754,
        0.4170771291151181,
        2.0614794296784096,
    )

    @property
    def totals(self) -> np.ndarray:
        return np.asarray(
            (self.cdc42_total, self.rac_total, self.rho_total), dtype=float
        )

    def __post_init__(self) -> None:
        positive = (
            self.cdc42_total, self.rac_total, self.rho_total,
            self.cdc42_activation, self.rho_activation,
            self.rho_inhibits_cdc42_half_level,
            self.rac_inhibits_rho_half_level, self.hill_exponent,
            self.inactive_accessibility,
        )
        if any(value <= 0.0 for value in positive):
            raise ValueError(
                "Holmes totals, rates, half-levels, hill exponent, and "
                "accessibility must be positive"
            )
        if self.rac_basal_activation < 0.0 or self.cdc42_to_rac_feedback < 0.0:
            raise ValueError("Holmes Rac activation terms must be non-negative")
        if self.active_diffusivity < 0.0 or self.inactive_diffusivity < 0.0:
            raise ValueError("Holmes diffusivities must be non-negative")
        active = np.asarray(self.initial_active, dtype=float)
        if (
            active.shape != (3,)
            or np.any(active < 0.0)
            or np.any(active >= self.totals)
        ):
            raise ValueError(
                "initial_active must contain three values between zero and "
                "their totals"
            )


class HolmesChemistry:
    """Canonical active/combined-inactive six-field Holmes reduction.

    Array row order is Cdc42, Rac, Rho.  ``membrane`` stores active fields and
    ``cytoplasm`` stores the combined inactive pools; the latter name is the
    migration framework's historical two-compartment label.
    """

    species = ("cdc42", "rac", "rho")

    def __init__(
        self,
        parameters: HolmesParameters | None = None,
        cue: Cue | None = None,
    ) -> None:
        self.parameters = parameters or HolmesParameters()
        self.cue = cue or NoCue()
        p = self.parameters
        self.membrane_diffusion = np.full(3, p.active_diffusivity, dtype=float)
        self.cytoplasm_diffusion = np.full(3, p.inactive_diffusivity, dtype=float)
        self.initial_membrane = np.asarray(p.initial_active, dtype=float)
        self.initial_cytoplasm = p.totals - self.initial_membrane

    def reaction(self, time, position, membrane, cytoplasm, strain, stress):
        _ = strain, stress
        p = self.parameters
        cdc42, rac, rho = membrane
        cdc42_i, rac_i, rho_i = cytoplasm
        signal = p.signal_baseline + self.cue(time, position)

        cdc42_rate = (
            p.cdc42_activation
            / (1.0 + (rho / p.rho_inhibits_cdc42_half_level) ** p.hill_exponent)
            * p.inactive_accessibility * cdc42_i / p.cdc42_total
            - cdc42
        )
        rac_rate = (
            (p.rac_basal_activation + p.cdc42_to_rac_feedback * cdc42 + signal)
            * p.inactive_accessibility * rac_i / p.rac_total
            - rac
        )
        rho_rate = (
            p.rho_activation
            / (1.0 + (rac / p.rac_inhibits_rho_half_level) ** p.hill_exponent)
            * p.inactive_accessibility * rho_i / p.rho_total
            - rho
        )
        active_rate = np.stack((cdc42_rate, rac_rate, rho_rate))
        return active_rate, -active_rate

    def protrusive_activity(self, membrane, cytoplasm):
        _ = cytoplasm
        return np.maximum(membrane[1], 0.0)

    def tension_deactivation(self, membrane, cytoplasm, added_rate):
        _ = cytoplasm
        loss = np.asarray(added_rate, dtype=float) * np.maximum(membrane[1], 0.0)
        active_rate = np.zeros_like(membrane)
        active_rate[1] = -loss
        return active_rate, -active_rate
