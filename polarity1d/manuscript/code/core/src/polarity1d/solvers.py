"""Fast and independently checkable solvers for the canonical 1D models."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Literal

import numpy as np
from scipy import sparse
from scipy.fft import dct, idct
from scipy.integrate import solve_ivp

from .contract import Grid1D, ReactionDiffusionModel, Stimulus, Trajectory, output_times
from .protocols import NoStimulus


Backend = Literal["bdf", "split", "rk45"]


def laplacian_no_flux(values: np.ndarray, dx: float) -> np.ndarray:
    """Conservative cell-centred second-order Laplacian."""
    result = np.empty_like(values)
    result[..., 0] = (values[..., 1] - values[..., 0]) / dx**2
    result[..., -1] = (values[..., -2] - values[..., -1]) / dx**2
    result[..., 1:-1] = (
        values[..., :-2] - 2.0 * values[..., 1:-1] + values[..., 2:]
    ) / dx**2
    return result


def _jacobian_sparsity(species: int, cells: int) -> sparse.csr_matrix:
    pattern = sparse.lil_matrix((species * cells, species * cells), dtype=np.int8)
    for field in range(species):
        rows = field * cells + np.arange(cells)
        pattern[rows, rows] = 1
        pattern[rows[1:], rows[:-1]] = 1
        pattern[rows[:-1], rows[1:]] = 1
        for other in range(species):
            pattern[rows, other * cells + np.arange(cells)] = 1
    return pattern.tocsr()


def _validate_initial(
    model: ReactionDiffusionModel,
    grid: Grid1D,
    provided: np.ndarray | None = None,
) -> np.ndarray:
    state = np.asarray(model.initial_state(grid) if provided is None else provided, dtype=float)
    expected = (len(model.state_names), grid.cells)
    if state.shape != expected:
        raise ValueError(f"initial state has shape {state.shape}, expected {expected}")
    if not np.all(np.isfinite(state)):
        raise ValueError("initial state contains non-finite values")
    if model.nonnegative and np.min(state) < 0:
        raise ValueError("nonnegative model has a negative initial state")
    if len(model.diffusivities) != expected[0] or min(model.diffusivities) < 0:
        raise ValueError("diffusivities do not match the state")
    return state


def _simulate_mol(
    model: ReactionDiffusionModel,
    grid: Grid1D,
    times: np.ndarray,
    stimulus: Stimulus,
    initial: np.ndarray,
    *,
    method: Literal["BDF", "RK45"],
    rtol: float,
    atol: float,
    max_step: float,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    species = initial.shape[0]
    diffusion = np.asarray(model.diffusivities)[:, None]

    def rhs(time: float, flattened: np.ndarray) -> np.ndarray:
        state = flattened.reshape(species, grid.cells)
        signal = np.asarray(stimulus(time, grid.x, grid.length), dtype=float)
        local = model.reaction(time, grid.x, state, signal)
        return (diffusion * laplacian_no_flux(state, grid.dx) + local).reshape(-1)

    kwargs = {}
    if method == "BDF":
        kwargs["jac_sparsity"] = _jacobian_sparsity(species, grid.cells)
    boundaries = np.asarray(
        [times[0]]
        + [event for event in stimulus.event_times if times[0] < event < times[-1]]
        + [times[-1]],
        dtype=float,
    )
    event_gaps = np.diff(np.unique(boundaries))
    effective_max_step = max_step
    if event_gaps.size:
        # Prevent an adaptive solver from stepping over an entire short pulse.
        effective_max_step = min(max_step, float(np.min(event_gaps)) / 4.0)
    start = perf_counter()
    solution = solve_ivp(
        rhs,
        (float(times[0]), float(times[-1])),
        initial.reshape(-1),
        method=method,
        t_eval=times,
        rtol=rtol,
        atol=atol,
        max_step=effective_max_step,
        **kwargs,
    )
    elapsed = perf_counter() - start
    if not solution.success:
        raise RuntimeError(solution.message)
    states = solution.y.T.reshape(len(times), species, grid.cells)
    if model.nonnegative and np.min(states) < -max(50.0 * atol, 1e-7):
        raise FloatingPointError(f"{method} produced a materially negative state")
    return states, {
        "backend": method.lower(),
        "wall_seconds": elapsed,
        "rhs_evaluations": int(solution.nfev),
        "jacobian_evaluations": int(getattr(solution, "njev", 0)),
        "linear_decompositions": int(getattr(solution, "nlu", 0)),
        "maximum_step": effective_max_step,
    }


@dataclass
class _DiffusionCache:
    grid: Grid1D
    diffusivities: np.ndarray

    def __post_init__(self) -> None:
        modes = np.arange(self.grid.cells)
        self.eigenvalues = -4.0 * np.sin(np.pi * modes / (2.0 * self.grid.cells)) ** 2 / self.grid.dx**2
        self._multipliers: dict[float, np.ndarray] = {}

    def half_step(self, state: np.ndarray, step: float) -> np.ndarray:
        key = float(np.round(step, 14))
        multiplier = self._multipliers.get(key)
        if multiplier is None:
            multiplier = np.exp(
                0.5 * step * self.diffusivities[:, None] * self.eigenvalues[None, :]
            )
            self._multipliers[key] = multiplier
        coefficients = dct(state, type=2, axis=1, norm="ortho")
        return idct(coefficients * multiplier, type=2, axis=1, norm="ortho")


def _reaction_rk4(
    model: ReactionDiffusionModel,
    stimulus: Stimulus,
    grid: Grid1D,
    state: np.ndarray,
    time: float,
    step: float,
) -> np.ndarray:
    def evaluate(at: float, values: np.ndarray) -> np.ndarray:
        signal = np.asarray(stimulus(at, grid.x, grid.length), dtype=float)
        return model.reaction(at, grid.x, values, signal)

    k1 = evaluate(time, state)
    k2 = evaluate(time + 0.5 * step, state + 0.5 * step * k1)
    k3 = evaluate(time + 0.5 * step, state + 0.5 * step * k2)
    k4 = evaluate(time + step, state + step * k3)
    return state + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def _simulate_split(
    model: ReactionDiffusionModel,
    grid: Grid1D,
    times: np.ndarray,
    stimulus: Stimulus,
    initial: np.ndarray,
    *,
    step: float,
) -> tuple[np.ndarray, dict[str, float | int | str]]:
    state = initial.copy()
    cache = _DiffusionCache(grid, np.asarray(model.diffusivities, dtype=float))
    states = np.empty((len(times), *state.shape))
    states[0] = state
    time = float(times[0])
    accepted = 0
    rejected = 0
    events = np.asarray(
        sorted(event for event in stimulus.event_times if times[0] < event < times[-1]),
        dtype=float,
    )
    start = perf_counter()
    for frame in range(1, len(times)):
        target = float(times[frame])
        while time < target - 1e-13:
            trial_step = min(step, target - time)
            future = events[events > time + 1e-13]
            if future.size:
                trial_step = min(trial_step, float(future[0] - time))
            while True:
                trial = cache.half_step(state, trial_step)
                trial = _reaction_rk4(model, stimulus, grid, trial, time, trial_step)
                trial = cache.half_step(trial, trial_step)
                invalid = not np.all(np.isfinite(trial))
                negative = model.nonnegative and np.min(trial) < -1e-9
                if not invalid and not negative:
                    break
                trial_step *= 0.5
                rejected += 1
                if trial_step < 1e-8:
                    raise FloatingPointError("split solver could not find a stable positive step")
            state = trial
            # Only eliminate roundoff-scale undershoots; substantive negatives retry above.
            if model.nonnegative:
                state[(state < 0) & (state > -1e-9)] = 0.0
            time += trial_step
            accepted += 1
        states[frame] = state
    elapsed = perf_counter() - start
    return states, {
        "backend": "split",
        "wall_seconds": elapsed,
        "accepted_steps": accepted,
        "rejected_steps": rejected,
        "fixed_step": step,
    }


def simulate_reaction_diffusion(
    model: ReactionDiffusionModel,
    *,
    protocol: Stimulus | None = None,
    length: float = 10.0,
    cells: int = 101,
    end_time: float = 300.0,
    output_interval: float = 1.0,
    times: np.ndarray | None = None,
    backend: Backend = "bdf",
    step: float | None = None,
    rtol: float = 2e-6,
    atol: float = 1e-9,
    max_step: float = 5.0,
    initial_state: np.ndarray | None = None,
) -> Trajectory:
    """Run any canonical model with an adaptive or split-step backend."""
    grid = Grid1D(length=length, cells=cells)
    initial = _validate_initial(model, grid, initial_state)
    protocol = NoStimulus() if protocol is None else protocol
    if times is None:
        times = output_times(end_time, output_interval)
    else:
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
            raise ValueError("times must start at zero and be strictly increasing")
    if backend == "split":
        states, diagnostics = _simulate_split(
            model,
            grid,
            times,
            protocol,
            initial,
            step=model.preferred_step if step is None else step,
        )
    elif backend in {"bdf", "rk45"}:
        states, diagnostics = _simulate_mol(
            model,
            grid,
            times,
            protocol,
            initial,
            method="BDF" if backend == "bdf" else "RK45",
            rtol=rtol,
            atol=atol,
            max_step=max_step,
        )
    else:
        raise ValueError(f"unknown backend {backend!r}")
    fields = {name: states[:, index] for index, name in enumerate(model.state_names)}
    metadata = {
        "framework": "polarity1d",
        "length_um": length,
        "cells": cells,
        "dx_um": grid.dx,
        "parameters": dict(model.metadata()),
        **diagnostics,
    }
    conservation = {}
    for group in getattr(model, "conserved_groups", ()):
        total = np.sum(states[:, group, :], axis=(1, 2))
        conservation["+".join(model.state_names[index] for index in group)] = float(
            np.max(np.abs(total / total[0] - 1.0))
        )
    if conservation:
        metadata["conservation_errors"] = conservation
    return Trajectory(
        model=model.name,
        time=times,
        x=grid.x,
        fields=fields,
        observable_name=model.observable_name,
        metadata=metadata,
    )
