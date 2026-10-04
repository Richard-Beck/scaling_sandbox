"""Reproducible initial-condition perturbations for noise qualification."""

from __future__ import annotations

import numpy as np

from .contract import Grid1D, ReactionDiffusionModel


def noisy_initial_state(
    model: ReactionDiffusionModel,
    grid: Grid1D,
    *,
    amplitude: float = 0.01,
    seed: int = 1,
) -> np.ndarray:
    """Perturb the observable with zero-mean uniform noise.

    For declared conserved groups the opposite perturbation is applied to the
    second pool, preserving total material pointwise as well as globally.
    This resolves the ambiguity in the printed Jilkine noise formula while
    retaining its intended concentration scale.
    """
    if amplitude < 0:
        raise ValueError("amplitude must be nonnegative")
    state = np.asarray(model.initial_state(grid), dtype=float).copy()
    if amplitude == 0:
        return state
    random = np.random.default_rng(seed).uniform(-0.5, 0.5, grid.cells)
    random -= np.mean(random)
    perturbation = amplitude * random
    groups = tuple(getattr(model, "conserved_groups", ()))
    observable = model.state_names.index(model.observable_name)
    paired = None
    for group in groups:
        if observable in group and len(group) == 2:
            paired = group[1] if group[0] == observable else group[0]
            break
    if paired is not None:
        limit = 0.9 * min(float(np.min(state[observable])), float(np.min(state[paired])))
        maximum = float(np.max(np.abs(perturbation)))
        if maximum > limit:
            perturbation *= limit / maximum
        state[observable] += perturbation
        state[paired] -= perturbation
    else:
        limit = 0.9 * float(np.min(state[observable]))
        if limit <= 0:
            # Add one-sided noise to a zero initial field without inventing negatives.
            perturbation -= np.min(perturbation)
        elif np.max(np.abs(perturbation)) > limit:
            perturbation *= limit / float(np.max(np.abs(perturbation)))
        state[observable] += perturbation
    return state

