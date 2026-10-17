"""Dynamic phase-field bulk coupled to an explicit moving 1D membrane."""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, pi
from pathlib import Path
import json

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from .core import ConservationPool, ReactionSystem, SimulationConfig
from .diffuse import (
    DiffuseGeometry,
    DiffuseGeometryFields,
    PhaseFieldMechanics,
    _PhaseFieldStepper,
    _advect,
    _bulk_diffusion_rhs,
    _project_material,
    _repair_small_negatives,
)
from .geometry import _bilinear_trace, _surface_operators


@dataclass
class MovingMembrane:
    """Ordered material markers representing a closed one-dimensional curve."""

    points: np.ndarray

    @classmethod
    def circle(cls, radius: float, n_points: int = 256) -> "MovingMembrane":
        if radius <= 0 or n_points < 8:
            raise ValueError("radius must be positive and n_points at least eight")
        theta = 2.0 * pi * np.arange(n_points, dtype=float) / n_points
        return cls(np.column_stack((radius * np.cos(theta), radius * np.sin(theta))))

    def __post_init__(self) -> None:
        self.points = np.ascontiguousarray(self.points, dtype=float)
        if self.points.ndim != 2 or self.points.shape[1] != 2:
            raise ValueError("membrane points must have shape (n_points, 2)")
        if self.points.shape[0] < 8:
            raise ValueError("a moving membrane needs at least eight points")
        self.weights, self.stiffness = _surface_operators(self.points)
        self.fraction = _curve_fraction(self.points)

    @property
    def n_points(self) -> int:
        return self.points.shape[0]

    @property
    def perimeter(self) -> float:
        return float(np.sum(self.weights))

    def update(self, points: np.ndarray) -> None:
        points = np.ascontiguousarray(points, dtype=float)
        if points.shape != self.points.shape:
            raise ValueError("updated membrane points have the wrong shape")
        self.points = points
        self.weights, self.stiffness = _surface_operators(points)
        self.fraction = _curve_fraction(points)


def _curve_fraction(points: np.ndarray) -> np.ndarray:
    edges = np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(edges)))
    return cumulative[:-1] / cumulative[-1]


@dataclass(frozen=True)
class MovingCurveInitialState:
    bulk: np.ndarray
    surface: np.ndarray


@dataclass
class MovingCurveResult:
    times: np.ndarray
    x: np.ndarray
    y: np.ndarray
    phi: np.ndarray
    bulk: np.ndarray
    surface: np.ndarray
    membrane_points: np.ndarray
    bulk_material: np.ndarray
    surface_material: np.ndarray
    bulk_species: tuple[str, ...]
    surface_species: tuple[str, ...]
    pool_totals: dict[str, np.ndarray]
    area: np.ndarray
    perimeter: np.ndarray
    centroid: np.ndarray
    maximum_edge_ratio: np.ndarray

    def save(self, output: str | Path) -> Path:
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "trajectory.npz"
        np.savez_compressed(
            path,
            times=self.times, x=self.x, y=self.y, phi=self.phi,
            bulk=self.bulk, surface=self.surface,
            membrane_points=self.membrane_points,
            bulk_material=self.bulk_material,
            surface_material=self.surface_material,
            bulk_species=np.asarray(self.bulk_species),
            surface_species=np.asarray(self.surface_species),
            area=self.area, perimeter=self.perimeter, centroid=self.centroid,
            maximum_edge_ratio=self.maximum_edge_ratio,
            **{f"pool_{name}": values for name, values in self.pool_totals.items()},
        )
        metadata = {
            "format": "bulk_surface_moving_curve_v1",
            "bulk_species": list(self.bulk_species),
            "surface_species": list(self.surface_species),
            "n_membrane_points": int(self.membrane_points.shape[1]),
            "maximum_relative_area_error": float(
                np.max(np.abs(self.area / self.area[0] - 1.0))
            ),
            "maximum_marker_edge_ratio": float(np.max(self.maximum_edge_ratio)),
            "maximum_relative_pool_errors": {
                name: float(np.max(np.abs(values / values[0] - 1.0)))
                for name, values in self.pool_totals.items() if values[0] != 0
            },
        }
        (directory / "metadata.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return path


def _grid_drive(points: np.ndarray, values: np.ndarray,
                geometry: DiffuseGeometry) -> np.ndarray:
    xx, yy = np.meshgrid(geometry.x, geometry.y)
    query = np.column_stack((xx.ravel(), yy.ravel()))
    nearest = cKDTree(points).query(query, k=1)[1]
    drive = values[nearest].reshape(geometry.shape)
    return gaussian_filter(drive, 0.65, mode="nearest")


def _project_to_interface(points: np.ndarray, geometry: DiffuseGeometry,
                          iterations: int = 4) -> np.ndarray:
    """Newton-project material markers to the new ``phi=0.5`` contour."""
    projected = points.copy()
    gy, gx = np.gradient(geometry.phi, geometry.dx, edge_order=2)
    for _ in range(iterations):
        trace = _bilinear_trace(geometry.x, geometry.y, projected)
        value = np.asarray(trace @ geometry.phi.ravel())
        grad_x = np.asarray(trace @ gx.ravel())
        grad_y = np.asarray(trace @ gy.ravel())
        norm2 = grad_x * grad_x + grad_y * grad_y
        if np.min(norm2) <= 1e-14:
            raise RuntimeError("cannot project membrane marker through a flat phase field")
        correction = (value - 0.5) / norm2
        projected[:, 0] -= correction * grad_x
        projected[:, 1] -= correction * grad_y
        if (np.min(projected[:, 0]) <= geometry.x[0]
                or np.max(projected[:, 0]) >= geometry.x[-1]
                or np.min(projected[:, 1]) <= geometry.y[0]
                or np.max(projected[:, 1]) >= geometry.y[-1]):
            raise RuntimeError("moving membrane reached the computational box boundary")
    return projected


class MovingCurveSimulator:
    """Phase-field bulk with a material, explicitly one-dimensional membrane."""

    def __init__(
        self,
        geometry: DiffuseGeometry,
        membrane: MovingMembrane,
        reactions: ReactionSystem,
        *,
        bulk_diffusion: np.ndarray,
        surface_diffusion: np.ndarray,
        mechanics: PhaseFieldMechanics,
        conservation_pools: tuple[ConservationPool, ...] = (),
        concentration_floor: float = 1e-6,
    ):
        self.geometry = geometry
        self.membrane = membrane
        self.reactions = reactions
        self.bulk_diffusion = np.asarray(bulk_diffusion, dtype=float)
        self.surface_diffusion = np.asarray(surface_diffusion, dtype=float)
        self.mechanics = mechanics
        if self.bulk_diffusion.shape != (len(reactions.bulk_species),):
            raise ValueError("bulk_diffusion has the wrong shape")
        if self.surface_diffusion.shape != (len(reactions.surface_species),):
            raise ValueError("surface_diffusion has the wrong shape")
        if mechanics.surface_force_coefficients.shape != (
            len(reactions.surface_species),
        ):
            raise ValueError("mechanical force coefficients have the wrong shape")
        if np.any(self.bulk_diffusion < 0) or np.any(self.surface_diffusion < 0):
            raise ValueError("diffusion coefficients must be nonnegative")
        self.conservation_pools = conservation_pools
        for pool in conservation_pools:
            pool.validate(reactions)
        self.concentration_floor = concentration_floor

    def _validate_state(self, initial: MovingCurveInitialState):
        bulk = np.ascontiguousarray(initial.bulk, dtype=float)
        surface = np.ascontiguousarray(initial.surface, dtype=float)
        expected_bulk = (len(self.reactions.bulk_species), *self.geometry.shape)
        expected_surface = (
            len(self.reactions.surface_species), self.membrane.n_points
        )
        if bulk.shape != expected_bulk or surface.shape != expected_surface:
            raise ValueError(
                f"state shapes must be {expected_bulk} and {expected_surface}"
            )
        if ((bulk.size and np.min(bulk) < 0)
                or (surface.size and np.min(surface) < 0)):
            raise ValueError("initial concentrations must be nonnegative")
        return (
            self.geometry.phi[None, :, :] * bulk,
            surface * self.membrane.weights[None, :],
        )

    def _concentrations(self, q_bulk: np.ndarray, surface_amount: np.ndarray):
        bulk = q_bulk / np.maximum(
            self.geometry.phi, self.concentration_floor
        )[None, :, :]
        surface = surface_amount / self.membrane.weights[None, :]
        return bulk, surface

    def _chemistry_rhs(self, time: float, q_bulk: np.ndarray,
                       surface_amount: np.ndarray):
        bulk, surface = self._concentrations(q_bulk, surface_amount)
        n_grid = self.geometry.size
        bulk_flat = bulk.reshape(bulk.shape[0], n_grid)
        trace = _bilinear_trace(
            self.geometry.x, self.geometry.y, self.membrane.points
        )
        bulk_at_membrane = np.vstack([
            trace @ bulk_flat[species] for species in range(bulk.shape[0])
        ]) if bulk.shape[0] else np.empty((0, self.membrane.n_points))
        bulk_local = np.zeros_like(bulk_flat)
        surface_local = np.zeros_like(surface)
        exchange = np.zeros(
            (self.reactions.n_exchange, self.membrane.n_points), dtype=float
        )
        stimulus = np.empty(self.membrane.n_points, dtype=float)
        self.reactions.stimulus.values(time, self.membrane.fraction, stimulus)
        self.reactions.bulk_kernel(
            time, bulk_flat, self.reactions.parameters, bulk_local
        )
        self.reactions.surface_kernel(
            time, surface, bulk_at_membrane, stimulus,
            self.reactions.parameters, surface_local, exchange,
        )
        bulk_rate = self.geometry.phi[None, :, :] * bulk_local.reshape(q_bulk.shape)
        surface_rate = self.membrane.weights[None, :] * surface_local
        if self.reactions.n_exchange:
            bulk_events = self.reactions.bulk_stoichiometry @ exchange
            surface_events = self.reactions.surface_stoichiometry @ exchange
            surface_rate += self.membrane.weights[None, :] * surface_events
            weighted = bulk_events * self.membrane.weights[None, :]
            for species in range(q_bulk.shape[0]):
                amount_rate = trace.T @ weighted[species]
                bulk_rate[species] += amount_rate.reshape(
                    self.geometry.shape
                ) / self.geometry.dx**2
        return bulk_rate, surface_rate

    def _chemistry_step(self, time: float, dt: float, q_bulk: np.ndarray,
                        surface_amount: np.ndarray):
        b1, s1 = self._chemistry_rhs(time, q_bulk, surface_amount)
        b2, s2 = self._chemistry_rhs(
            time + 0.5 * dt,
            q_bulk + 0.5 * dt * b1,
            surface_amount + 0.5 * dt * s1,
        )
        return q_bulk + dt * b2, surface_amount + dt * s2

    def _diffusion_step(self, dt: float, q_bulk: np.ndarray,
                        surface_amount: np.ndarray):
        maximum_bulk = (
            float(np.max(self.bulk_diffusion)) if self.bulk_diffusion.size else 0.0
        )
        bulk_steps = max(
            1, int(ceil(maximum_bulk * dt / (0.20 * self.geometry.dx**2)))
        )
        for _ in range(bulk_steps):
            bulk, _ = self._concentrations(q_bulk, surface_amount)
            q_bulk += dt / bulk_steps * _bulk_diffusion_rhs(
                bulk, self.geometry.phi, self.bulk_diffusion, self.geometry.dx
            )
            q_bulk = _repair_small_negatives(q_bulk)

        edges = np.linalg.norm(
            np.roll(self.membrane.points, -1, axis=0) - self.membrane.points,
            axis=1,
        )
        maximum_surface = (
            float(np.max(self.surface_diffusion))
            if self.surface_diffusion.size else 0.0
        )
        surface_steps = max(
            1, int(ceil(maximum_surface * dt / (0.20 * np.min(edges)**2)))
        )
        for _ in range(surface_steps):
            _, surface = self._concentrations(q_bulk, surface_amount)
            for species, coefficient in enumerate(self.surface_diffusion):
                surface_amount[species] -= (
                    dt / surface_steps * coefficient
                    * (self.membrane.stiffness @ surface[species])
                )
            surface_amount = _repair_small_negatives(
                surface_amount[:, None, :]
            )[:, 0, :]
        return q_bulk, surface_amount

    def _move_geometry(self, phase_stepper: _PhaseFieldStepper,
                       q_bulk: np.ndarray, surface_amount: np.ndarray,
                       dt: float):
        if self.mechanics.frozen:
            zeros = np.zeros_like(self.geometry.phi)
            return zeros, zeros, surface_amount
        _, surface = self._concentrations(q_bulk, surface_amount)
        force = self.mechanics.surface_force_coefficients @ surface
        drive = _grid_drive(self.membrane.points, force, self.geometry)
        old_trace = _bilinear_trace(
            self.geometry.x, self.geometry.y, self.membrane.points
        )
        new_phi, vx, vy = phase_stepper.step_from_drive(drive)
        marker_velocity = np.column_stack((
            old_trace @ vx.ravel(), old_trace @ vy.ravel()
        ))
        moved = self.membrane.points + dt * marker_velocity
        self.geometry.phi = new_phi
        moved = _project_to_interface(moved, self.geometry)
        self.membrane.update(moved)
        edges = np.linalg.norm(
            np.roll(moved, -1, axis=0) - moved, axis=1
        )
        if float(np.max(edges) / np.min(edges)) > 1.25:
            old_totals = np.sum(surface_amount, axis=1)
            surface = surface_amount / self.membrane.weights[None, :]
            old_fraction = self.membrane.fraction
            closed_fraction = np.concatenate((old_fraction, [1.0]))
            target_fraction = np.arange(self.membrane.n_points) / self.membrane.n_points
            closed_points = np.vstack((moved, moved[0]))
            redistributed_points = np.column_stack([
                np.interp(target_fraction, closed_fraction, closed_points[:, axis])
                for axis in range(2)
            ])
            redistributed_points = _project_to_interface(
                redistributed_points, self.geometry
            )
            redistributed_surface = np.vstack([
                np.interp(
                    target_fraction,
                    closed_fraction,
                    np.concatenate((values, [values[0]])),
                )
                for values in surface
            ])
            self.membrane.update(redistributed_points)
            surface_amount = redistributed_surface * self.membrane.weights[None, :]
            new_totals = np.sum(surface_amount, axis=1)
            for species, total in enumerate(old_totals):
                surface_amount[species] *= total / new_totals[species]
        return vx, vy, surface_amount

    def _pool_total(self, pool: ConservationPool, q_bulk: np.ndarray,
                    surface_amount: np.ndarray) -> float:
        bulk_amount = np.sum(q_bulk, axis=(1, 2)) * self.geometry.dx**2
        membrane_amount = np.sum(surface_amount, axis=1)
        return float(
            np.asarray(pool.bulk_coefficients) @ bulk_amount
            + np.asarray(pool.surface_coefficients) @ membrane_amount
        )

    def run(self, initial: MovingCurveInitialState,
            config: SimulationConfig) -> MovingCurveResult:
        q_bulk, surface_amount = self._validate_state(initial)
        n_steps = round(config.end_time / config.dt)
        save_every = round(config.output_interval / config.dt)
        save_steps = np.arange(0, n_steps + 1, save_every, dtype=int)
        if save_steps[-1] != n_steps:
            save_steps = np.append(save_steps, n_steps)
        n_save = save_steps.size
        phi_history = np.empty((n_save, *self.geometry.shape), dtype=np.float32)
        bulk_history = np.empty(
            (n_save, q_bulk.shape[0], *self.geometry.shape), dtype=np.float32
        )
        surface_history = np.empty(
            (n_save, surface_amount.shape[0], self.membrane.n_points),
            dtype=np.float32,
        )
        point_history = np.empty(
            (n_save, self.membrane.n_points, 2), dtype=np.float32
        )
        bulk_material = np.empty((n_save, q_bulk.shape[0]), dtype=float)
        membrane_material = np.empty((n_save, surface_amount.shape[0]), dtype=float)
        area = np.empty(n_save)
        perimeter = np.empty(n_save)
        centroid = np.empty((n_save, 2))
        edge_ratio = np.empty(n_save)
        pools = {
            pool.name: np.empty(n_save) for pool in self.conservation_pools
        }

        def save(index: int) -> None:
            fields = self.geometry.fields()
            bulk, surface = self._concentrations(q_bulk, surface_amount)
            edges = np.linalg.norm(
                np.roll(self.membrane.points, -1, axis=0) - self.membrane.points,
                axis=1,
            )
            phi_history[index] = self.geometry.phi
            bulk_history[index] = bulk
            surface_history[index] = surface
            point_history[index] = self.membrane.points
            bulk_material[index] = np.sum(q_bulk, axis=(1, 2)) * self.geometry.dx**2
            membrane_material[index] = np.sum(surface_amount, axis=1)
            area[index] = fields.area
            perimeter[index] = self.membrane.perimeter
            centroid[index] = (fields.centroid_x, fields.centroid_y)
            edge_ratio[index] = float(np.max(edges) / np.min(edges))
            for pool in self.conservation_pools:
                pools[pool.name][index] = self._pool_total(
                    pool, q_bulk, surface_amount
                )

        save(0)
        phase_stepper = _PhaseFieldStepper(self.geometry, self.mechanics, config.dt)
        half = 0.5 * config.dt
        next_save_index = 1
        for step in range(1, n_steps + 1):
            time = (step - 1) * config.dt
            q_bulk, surface_amount = self._chemistry_step(
                time, half, q_bulk, surface_amount
            )
            q_bulk, surface_amount = self._diffusion_step(
                half, q_bulk, surface_amount
            )
            bulk_totals = np.sum(q_bulk, axis=(1, 2))
            vx, vy, surface_amount = self._move_geometry(
                phase_stepper, q_bulk, surface_amount, config.dt
            )
            q_bulk = _advect(q_bulk, vx, vy, self.geometry.dx, config.dt)
            q_bulk = _project_material(
                q_bulk, self.geometry.phi, self.concentration_floor
            )
            current = np.sum(q_bulk, axis=(1, 2))
            for species, total in enumerate(bulk_totals):
                if current[species] != 0:
                    q_bulk[species] *= total / current[species]
            q_bulk, surface_amount = self._diffusion_step(
                half, q_bulk, surface_amount
            )
            q_bulk, surface_amount = self._chemistry_step(
                time + half, half, q_bulk, surface_amount
            )
            if (not np.all(np.isfinite(q_bulk))
                    or not np.all(np.isfinite(surface_amount))):
                raise RuntimeError("moving-curve state became non-finite")
            if ((q_bulk.size and np.min(q_bulk) < -1e-8)
                    or (surface_amount.size and np.min(surface_amount) < -1e-8)):
                raise RuntimeError("moving-curve state became negative")
            if next_save_index < n_save and step == save_steps[next_save_index]:
                save(next_save_index)
                next_save_index += 1

        return MovingCurveResult(
            times=save_steps.astype(float) * config.dt,
            x=self.geometry.x.copy(), y=self.geometry.y.copy(),
            phi=phi_history, bulk=bulk_history, surface=surface_history,
            membrane_points=point_history,
            bulk_material=bulk_material, surface_material=membrane_material,
            bulk_species=self.reactions.bulk_species,
            surface_species=self.reactions.surface_species,
            pool_totals=pools, area=area, perimeter=perimeter,
            centroid=centroid, maximum_edge_ratio=edge_ratio,
        )
