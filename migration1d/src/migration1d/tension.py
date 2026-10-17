"""Optional tension-to-polarity feedback for the moving 1D model."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class GlobalTensionFeedback:
    """Uniform Rac-equivalent deactivation driven by elastic extension.

    The solver uses the positive part of its existing spring force,
    ``max(0, k * (L - L0))``.  Dashpot force and compression are excluded.
    """

    alpha_t: float = 0.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.alpha_t) or self.alpha_t < 0.0:
            raise ValueError("alpha_t must be finite and non-negative")


@dataclass(frozen=True)
class SpatialTensionFeedback:
    """Diffusive excess-tension field coupled locally to polarity activity.

    The integrated field is ``theta = T - baseline`` and obeys
    ``D_T theta_xx - theta/tau_T + beta*P(a)`` with no-flux boundaries.
    """

    diffusivity: float
    relaxation_time: float
    generation_rate: float
    alpha_t: float
    baseline: float = 0.0

    def __post_init__(self) -> None:
        finite = (
            self.diffusivity,
            self.relaxation_time,
            self.generation_rate,
            self.alpha_t,
            self.baseline,
        )
        if not all(math.isfinite(value) for value in finite):
            raise ValueError("spatial tension parameters must be finite")
        if self.diffusivity < 0.0 or self.generation_rate < 0.0 or self.alpha_t < 0.0:
            raise ValueError(
                "tension diffusivity, generation_rate, and alpha_t must be non-negative"
            )
        if self.relaxation_time <= 0.0:
            raise ValueError("relaxation_time must be positive")
        if self.baseline < 0.0:
            raise ValueError("baseline tension must be non-negative")

    @property
    def characteristic_length(self) -> float:
        """Communication length ``sqrt(D_T*tau_T)``."""

        return math.sqrt(self.diffusivity * self.relaxation_time)


TensionFeedback = GlobalTensionFeedback | SpatialTensionFeedback
