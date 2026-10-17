"""One public simulation function for every WP1 model."""

from __future__ import annotations

from typing import Any

from .contract import ReactionDiffusionModel, Trajectory
from .debelly import DeBelly, simulate_debelly
from .solvers import simulate_reaction_diffusion


def simulate(model: ReactionDiffusionModel | DeBelly, **kwargs: Any) -> Trajectory:
    """Dispatch through the shared model contract while retaining best solver per class."""
    if isinstance(model, DeBelly):
        unsupported = {key for key in ("backend", "rtol", "atol", "max_step") if key in kwargs}
        if unsupported:
            raise TypeError(f"De Belly specialized solver does not accept {sorted(unsupported)}")
        return simulate_debelly(model, **kwargs)
    if isinstance(model, ReactionDiffusionModel):
        return simulate_reaction_diffusion(model, **kwargs)
    raise TypeError(f"unsupported model type {type(model).__name__}")

