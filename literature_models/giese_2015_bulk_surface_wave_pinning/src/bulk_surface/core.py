"""Common bulk--surface state, reaction, integration, and output interfaces."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Callable, Protocol

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import factorized

from .geometry import FixedGeometry
from .protocols import NoStimulus


BulkKernel = Callable[[float, np.ndarray, np.ndarray, np.ndarray], None]
SurfaceKernel = Callable[
    [float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    None,
]


class Stimulus(Protocol):
    def values(self, time: float, membrane_fraction: np.ndarray,
               out: np.ndarray) -> None: ...


def _zero_bulk_kernel(time: float, bulk: np.ndarray,
                      parameters: np.ndarray, out: np.ndarray) -> None:
    _ = time, bulk, parameters
    out.fill(0.0)


def _zero_surface_kernel(
    time: float,
    surface: np.ndarray,
    bulk_trace: np.ndarray,
    stimulus: np.ndarray,
    parameters: np.ndarray,
    surface_local_out: np.ndarray,
    exchange_out: np.ndarray,
) -> None:
    _ = time, surface, bulk_trace, stimulus, parameters
    surface_local_out.fill(0.0)
    exchange_out.fill(0.0)


@dataclass(frozen=True)
class ReactionSystem:
    """Configurable local reactions and membrane exchange channels.

    An exchange channel is an event rate per membrane length.  Its bulk and
    surface stoichiometry columns specify the change caused by one event.
    """

    bulk_species: tuple[str, ...]
    surface_species: tuple[str, ...]
    parameters: np.ndarray = field(default_factory=lambda: np.empty(0))
    bulk_stoichiometry: np.ndarray = field(
        default_factory=lambda: np.empty((0, 0))
    )
    surface_stoichiometry: np.ndarray = field(
        default_factory=lambda: np.empty((0, 0))
    )
    bulk_kernel: BulkKernel = _zero_bulk_kernel
    surface_kernel: SurfaceKernel = _zero_surface_kernel
    stimulus: Stimulus = field(default_factory=NoStimulus)

    def __post_init__(self) -> None:
        parameters = np.ascontiguousarray(self.parameters, dtype=float)
        bulk_stoich = np.ascontiguousarray(self.bulk_stoichiometry, dtype=float)
        surface_stoich = np.ascontiguousarray(self.surface_stoichiometry, dtype=float)
        if len(set(self.bulk_species + self.surface_species)) != (
            len(self.bulk_species) + len(self.surface_species)
        ):
            raise ValueError("species names must be unique")
        if bulk_stoich.ndim != 2 or bulk_stoich.shape[0] != len(self.bulk_species):
            raise ValueError("bulk_stoichiometry has the wrong shape")
        if surface_stoich.ndim != 2 or surface_stoich.shape[0] != len(
            self.surface_species
        ):
            raise ValueError("surface_stoichiometry has the wrong shape")
        if bulk_stoich.shape[1] != surface_stoich.shape[1]:
            raise ValueError("bulk and surface stoichiometry need equal channel counts")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "bulk_stoichiometry", bulk_stoich)
        object.__setattr__(self, "surface_stoichiometry", surface_stoich)

    @property
    def n_exchange(self) -> int:
        return self.bulk_stoichiometry.shape[1]


@dataclass(frozen=True)
class ConservationPool:
    name: str
    bulk_coefficients: np.ndarray
    surface_coefficients: np.ndarray

    def validate(self, reactions: ReactionSystem) -> None:
        if np.asarray(self.bulk_coefficients).shape != (len(reactions.bulk_species),):
            raise ValueError(f"pool {self.name!r} has wrong bulk coefficients")
        if np.asarray(self.surface_coefficients).shape != (
            len(reactions.surface_species),
        ):
            raise ValueError(f"pool {self.name!r} has wrong surface coefficients")
        exchange_balance = (
            np.asarray(self.bulk_coefficients) @ reactions.bulk_stoichiometry
            + np.asarray(self.surface_coefficients) @ reactions.surface_stoichiometry
        )
        if not np.allclose(exchange_balance, 0.0, atol=1e-14):
            raise ValueError(f"pool {self.name!r} is not exchange-conserved")


@dataclass(frozen=True)
class InitialState:
    bulk: np.ndarray
    surface: np.ndarray


@dataclass(frozen=True)
class SimulationConfig:
    end_time: float
    dt: float
    output_interval: float

    def __post_init__(self) -> None:
        if self.end_time <= 0 or self.dt <= 0 or self.output_interval <= 0:
            raise ValueError("end_time, dt, and output_interval must be positive")
        for value, label in (
            (self.end_time / self.dt, "end_time"),
            (self.output_interval / self.dt, "output_interval"),
        ):
            if not np.isclose(value, round(value), rtol=0.0, atol=1e-10):
                raise ValueError(f"{label} must be divisible by dt")


@dataclass
class SimulationResult:
    times: np.ndarray
    bulk: np.ndarray
    surface: np.ndarray
    bulk_species: tuple[str, ...]
    surface_species: tuple[str, ...]
    pool_totals: dict[str, np.ndarray]
    geometry: FixedGeometry

    def save(self, output: str | Path) -> Path:
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            directory / "trajectory.npz",
            times=self.times,
            bulk=self.bulk,
            surface=self.surface,
            x=self.geometry.x,
            y=self.geometry.y,
            phi=self.geometry.phi,
            membrane_points=self.geometry.membrane_points,
            membrane_fraction=self.geometry.membrane_fraction,
            bulk_species=np.asarray(self.bulk_species),
            surface_species=np.asarray(self.surface_species),
            **{f"pool_{name}": values for name, values in self.pool_totals.items()},
        )
        metadata = {
            "bulk_species": list(self.bulk_species),
            "surface_species": list(self.surface_species),
            "bulk_shape": list(self.geometry.bulk_shape),
            "area": self.geometry.area,
            "perimeter": self.geometry.perimeter,
            "n_surface": self.geometry.n_surface,
            "maximum_relative_pool_errors": {
                name: float(np.max(np.abs(values / values[0] - 1.0)))
                for name, values in self.pool_totals.items()
                if values[0] != 0
            },
        }
        (directory / "metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return directory / "trajectory.npz"


class _DiffusionStepper:
    def __init__(self, geometry: FixedGeometry, dt: float):
        self.geometry = geometry
        self.dt = dt
        self._bulk: dict[float, Callable[[np.ndarray], np.ndarray]] = {}
        self._surface: dict[float, Callable[[np.ndarray], np.ndarray]] = {}

    @staticmethod
    def _make_factor(weights: np.ndarray, stiffness: sparse.csr_matrix,
                     dt: float, diffusion: float):
        operator = sparse.diags(weights) + dt * diffusion * stiffness
        return factorized(operator.tocsc())

    def bulk(self, values: np.ndarray, diffusion: float) -> np.ndarray:
        if diffusion == 0:
            return values.copy()
        if diffusion < 0:
            raise ValueError("diffusion coefficients must be nonnegative")
        solve = self._bulk.get(diffusion)
        if solve is None:
            solve = self._make_factor(
                self.geometry.bulk_weights,
                self.geometry.bulk_stiffness,
                self.dt,
                diffusion,
            )
            self._bulk[diffusion] = solve
        return np.asarray(solve(self.geometry.bulk_weights * values))

    def surface(self, values: np.ndarray, diffusion: float) -> np.ndarray:
        if diffusion == 0:
            return values.copy()
        if diffusion < 0:
            raise ValueError("diffusion coefficients must be nonnegative")
        solve = self._surface.get(diffusion)
        if solve is None:
            solve = self._make_factor(
                self.geometry.surface_weights,
                self.geometry.surface_stiffness,
                self.dt,
                diffusion,
            )
            self._surface[diffusion] = solve
        return np.asarray(solve(self.geometry.surface_weights * values))


class Simulator:
    """Strang-split reaction/exchange and implicit diffusion simulator."""

    def __init__(
        self,
        geometry: FixedGeometry,
        reactions: ReactionSystem,
        *,
        bulk_diffusion: np.ndarray,
        surface_diffusion: np.ndarray,
        conservation_pools: tuple[ConservationPool, ...] = (),
    ):
        self.geometry = geometry
        self.reactions = reactions
        self.bulk_diffusion = np.asarray(bulk_diffusion, dtype=float)
        self.surface_diffusion = np.asarray(surface_diffusion, dtype=float)
        if self.bulk_diffusion.shape != (len(reactions.bulk_species),):
            raise ValueError("bulk_diffusion has the wrong shape")
        if self.surface_diffusion.shape != (len(reactions.surface_species),):
            raise ValueError("surface_diffusion has the wrong shape")
        if np.any(self.bulk_diffusion < 0) or np.any(self.surface_diffusion < 0):
            raise ValueError("diffusion coefficients must be nonnegative")
        self.conservation_pools = conservation_pools
        for pool in conservation_pools:
            pool.validate(reactions)

    def _validate_state(self, initial: InitialState) -> tuple[np.ndarray, np.ndarray]:
        bulk = np.ascontiguousarray(initial.bulk, dtype=float).copy()
        surface = np.ascontiguousarray(initial.surface, dtype=float).copy()
        expected_bulk = (len(self.reactions.bulk_species), self.geometry.n_bulk)
        expected_surface = (
            len(self.reactions.surface_species), self.geometry.n_surface
        )
        if bulk.shape != expected_bulk:
            raise ValueError(f"bulk state must have shape {expected_bulk}")
        if surface.shape != expected_surface:
            raise ValueError(f"surface state must have shape {expected_surface}")
        if not np.all(np.isfinite(bulk)) or not np.all(np.isfinite(surface)):
            raise ValueError("initial state must be finite")
        return bulk, surface

    def _pool_total(self, pool: ConservationPool, bulk: np.ndarray,
                    surface: np.ndarray) -> float:
        bulk_mass = bulk @ self.geometry.bulk_weights
        surface_mass = surface @ self.geometry.surface_weights
        return float(
            np.asarray(pool.bulk_coefficients) @ bulk_mass
            + np.asarray(pool.surface_coefficients) @ surface_mass
        )

    def _rhs(
        self,
        time: float,
        bulk: np.ndarray,
        surface: np.ndarray,
        bulk_out: np.ndarray,
        surface_out: np.ndarray,
        bulk_trace: np.ndarray,
        stimulus: np.ndarray,
        exchange: np.ndarray,
    ) -> None:
        bulk_out.fill(0.0)
        surface_out.fill(0.0)
        exchange.fill(0.0)
        for species in range(bulk.shape[0]):
            bulk_trace[species] = self.geometry.trace @ bulk[species]
        self.reactions.stimulus.values(
            time, self.geometry.membrane_fraction, stimulus
        )
        self.reactions.bulk_kernel(
            time, bulk, self.reactions.parameters, bulk_out
        )
        self.reactions.surface_kernel(
            time,
            surface,
            bulk_trace,
            stimulus,
            self.reactions.parameters,
            surface_out,
            exchange,
        )
        if self.reactions.n_exchange == 0:
            return
        surface_out += self.reactions.surface_stoichiometry @ exchange
        bulk_boundary_rates = self.reactions.bulk_stoichiometry @ exchange
        weighted = bulk_boundary_rates * self.geometry.surface_weights
        for species in range(bulk.shape[0]):
            amount_rate = self.geometry.trace.T @ weighted[species]
            bulk_out[species] += amount_rate / self.geometry.bulk_weights

    def _reaction_step(
        self,
        time: float,
        step: float,
        bulk: np.ndarray,
        surface: np.ndarray,
        workspace: tuple[np.ndarray, ...],
    ) -> None:
        (bulk_k1, surface_k1, bulk_k2, surface_k2, bulk_mid, surface_mid,
         bulk_trace, stimulus, exchange) = workspace
        self._rhs(
            time, bulk, surface, bulk_k1, surface_k1,
            bulk_trace, stimulus, exchange,
        )
        bulk_mid[:] = bulk + 0.5 * step * bulk_k1
        surface_mid[:] = surface + 0.5 * step * surface_k1
        self._rhs(
            time + 0.5 * step, bulk_mid, surface_mid, bulk_k2, surface_k2,
            bulk_trace, stimulus, exchange,
        )
        bulk += step * bulk_k2
        surface += step * surface_k2

    def run(self, initial: InitialState, config: SimulationConfig) -> SimulationResult:
        bulk, surface = self._validate_state(initial)
        n_steps = round(config.end_time / config.dt)
        save_every = round(config.output_interval / config.dt)
        save_steps = np.arange(0, n_steps + 1, save_every, dtype=int)
        if save_steps[-1] != n_steps:
            save_steps = np.append(save_steps, n_steps)
        n_save = save_steps.size
        times = save_steps.astype(float) * config.dt
        bulk_history = np.empty((n_save, *bulk.shape), dtype=float)
        surface_history = np.empty((n_save, *surface.shape), dtype=float)
        pool_totals = {
            pool.name: np.empty(n_save, dtype=float) for pool in self.conservation_pools
        }
        bulk_history[0] = bulk
        surface_history[0] = surface
        for pool in self.conservation_pools:
            pool_totals[pool.name][0] = self._pool_total(pool, bulk, surface)

        workspace = (
            np.empty_like(bulk), np.empty_like(surface),
            np.empty_like(bulk), np.empty_like(surface),
            np.empty_like(bulk), np.empty_like(surface),
            np.empty((bulk.shape[0], self.geometry.n_surface)),
            np.empty(self.geometry.n_surface),
            np.empty((self.reactions.n_exchange, self.geometry.n_surface)),
        )
        diffusion = _DiffusionStepper(self.geometry, config.dt)
        save_index = 1
        next_save = save_steps[save_index] if save_index < n_save else -1
        half = 0.5 * config.dt
        for step_index in range(1, n_steps + 1):
            time = (step_index - 1) * config.dt
            self._reaction_step(time, half, bulk, surface, workspace)
            for species, coefficient in enumerate(self.bulk_diffusion):
                bulk[species] = diffusion.bulk(bulk[species], float(coefficient))
            for species, coefficient in enumerate(self.surface_diffusion):
                surface[species] = diffusion.surface(
                    surface[species], float(coefficient)
                )
            self._reaction_step(time + half, half, bulk, surface, workspace)
            if step_index == next_save:
                bulk_history[save_index] = bulk
                surface_history[save_index] = surface
                for pool in self.conservation_pools:
                    pool_totals[pool.name][save_index] = self._pool_total(
                        pool, bulk, surface
                    )
                save_index += 1
                next_save = save_steps[save_index] if save_index < n_save else -1

        return SimulationResult(
            times=times,
            bulk=bulk_history,
            surface=surface_history,
            bulk_species=self.reactions.bulk_species,
            surface_species=self.reactions.surface_species,
            pool_totals=pool_totals,
            geometry=self.geometry,
        )
