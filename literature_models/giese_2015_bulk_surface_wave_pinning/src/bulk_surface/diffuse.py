"""Dynamic diffuse-domain bulk--surface reaction--diffusion mechanics.

This module complements the fixed polygon/Cartesian backend in :mod:`core`.
The cell is represented on a fixed laboratory grid by a phase field, while
bulk and membrane species are stored as the conservative densities
``q_bulk = phi*c`` and ``q_surface = delta(phi)*u``.  Bulk--surface exchange
is evaluated once and applied with paired stoichiometry on the same diffuse
interface, so declared pools remain conservative while the cell deforms.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
from math import ceil, pi, sqrt

import numpy as np
from scipy import sparse
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.sparse.linalg import factorized

from .core import ConservationPool, InitialState, ReactionSystem, SimulationConfig


@dataclass
class DiffuseGeometry:
    """A mutable phase field on a fixed square laboratory grid."""

    x: np.ndarray
    y: np.ndarray
    phi: np.ndarray
    dx: float
    epsilon: float
    target_area: float

    @classmethod
    def circle(
        cls,
        radius: float = 5.0,
        *,
        dx: float = 0.25,
        interface_width: float | None = None,
        box_half_width: float | None = None,
        padding_widths: float = 8.0,
    ) -> "DiffuseGeometry":
        if radius <= 0 or dx <= 0:
            raise ValueError("radius and dx must be positive")
        epsilon = float(interface_width if interface_width is not None else dx)
        if epsilon <= 0 or padding_widths < 4:
            raise ValueError("interface_width must be positive and padding_widths >= 4")
        half = float(box_half_width or radius + padding_widths * epsilon)
        if half <= radius + 3.0 * epsilon:
            raise ValueError("the computational box needs at least three interface widths of padding")
        n_half = ceil(half / dx)
        x = np.arange(-n_half, n_half + 1, dtype=float) * dx
        y = x.copy()
        xx, yy = np.meshgrid(x, y)
        phi = 0.5 * (
            1.0 - np.tanh((np.hypot(xx, yy) - radius) / (sqrt(2.0) * epsilon))
        )
        target = float(np.sum(phi) * dx * dx)
        return cls(x=x, y=y, phi=phi, dx=dx, epsilon=epsilon,
                   target_area=target)

    @property
    def shape(self) -> tuple[int, int]:
        return self.phi.shape

    @property
    def size(self) -> int:
        return self.phi.size

    def fields(self) -> "DiffuseGeometryFields":
        return geometry_fields(self.phi, self.x, self.y, self.dx)


@dataclass(frozen=True)
class DiffuseGeometryFields:
    delta: np.ndarray
    normal_x: np.ndarray
    normal_y: np.ndarray
    curvature: np.ndarray
    projection_xx: np.ndarray
    projection_xy: np.ndarray
    projection_yy: np.ndarray
    membrane_fraction: np.ndarray
    area: float
    perimeter: float
    centroid_x: float
    centroid_y: float


def _centered_gradient(values: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray]:
    gy, gx = np.gradient(values, dx, edge_order=2)
    return gx, gy


def geometry_fields(
    phi: np.ndarray, x: np.ndarray, y: np.ndarray, dx: float
) -> DiffuseGeometryFields:
    """Derive one mutually consistent set of diffuse-interface quantities."""
    gx, gy = _centered_gradient(phi, dx)
    delta = np.hypot(gx, gy)
    regularizer = max(float(np.max(delta)) * 1e-8, 1e-14)
    magnitude = np.sqrt(gx * gx + gy * gy + regularizer * regularizer)
    # phi is one inside and decreases outwards, hence -grad(phi) is outward.
    nx = -gx / magnitude
    ny = -gy / magnitude
    nx_x, _ = _centered_gradient(nx, dx)
    _, ny_y = _centered_gradient(ny, dx)
    curvature = nx_x + ny_y
    pxx = 1.0 - nx * nx
    pxy = -nx * ny
    pyy = 1.0 - ny * ny
    cell_area = float(np.sum(phi) * dx * dx)
    xx, yy = np.meshgrid(x, y)
    cx = float(np.sum(phi * xx) * dx * dx / cell_area)
    cy = float(np.sum(phi * yy) * dx * dx / cell_area)
    angle = np.mod(np.arctan2(yy - cy, xx - cx), 2.0 * pi)
    return DiffuseGeometryFields(
        delta=delta,
        normal_x=nx,
        normal_y=ny,
        curvature=curvature,
        projection_xx=pxx,
        projection_xy=pxy,
        projection_yy=pyy,
        membrane_fraction=angle.ravel() / (2.0 * pi),
        area=cell_area,
        perimeter=float(np.sum(delta) * dx * dx),
        centroid_x=cx,
        centroid_y=cy,
    )


@dataclass(frozen=True)
class PhaseFieldMechanics:
    """Overdamped interface mechanics with optional Kelvin--Voigt tension.

    ``surface_force_coefficients`` maps surface concentrations to local
    outward force. Positive entries are protrusive and negative entries are
    contractile. The Allen--Cahn energy already supplies baseline curvature
    tension. ``perimeter_elasticity`` and ``perimeter_viscosity`` add a
    separate global cortical tension when a Kelvin--Voigt analogue is wanted.
    """

    surface_force_coefficients: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=float)
    )
    relaxation_time: float = 1.0
    surface_tension: float = 1.0
    area_penalty: float = 0.1
    remove_mean_force: bool = True
    enforce_area: bool = True
    perimeter_elasticity: float = 0.0
    perimeter_viscosity: float = 0.0
    rest_perimeter: float | None = None
    velocity_smoothing: float = 1.0
    frozen: bool = False

    def __post_init__(self) -> None:
        coefficients = np.ascontiguousarray(
            self.surface_force_coefficients, dtype=float
        )
        object.__setattr__(self, "surface_force_coefficients", coefficients)
        if self.relaxation_time <= 0 or self.surface_tension < 0:
            raise ValueError("relaxation_time must be positive and tension nonnegative")
        if self.area_penalty < 0 or self.perimeter_elasticity < 0:
            raise ValueError("area and perimeter elasticities must be nonnegative")
        if self.perimeter_viscosity < 0 or self.velocity_smoothing < 0:
            raise ValueError("viscosity and velocity smoothing must be nonnegative")
        if self.rest_perimeter is not None and self.rest_perimeter <= 0:
            raise ValueError("rest_perimeter must be positive")


@dataclass(frozen=True)
class DiffuseInitialState:
    """Initial bulk and membrane concentrations on the Cartesian grid."""

    bulk: np.ndarray
    surface: np.ndarray


@dataclass
class DiffuseSimulationResult:
    times: np.ndarray
    x: np.ndarray
    y: np.ndarray
    phi: np.ndarray
    bulk: np.ndarray
    surface: np.ndarray
    bulk_material: np.ndarray
    surface_material: np.ndarray
    bulk_species: tuple[str, ...]
    surface_species: tuple[str, ...]
    pool_totals: dict[str, np.ndarray]
    area: np.ndarray
    perimeter: np.ndarray
    centroid: np.ndarray

    def save(self, output: str | Path) -> Path:
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / "trajectory.npz"
        np.savez_compressed(
            path,
            times=self.times,
            x=self.x,
            y=self.y,
            phi=self.phi,
            bulk=self.bulk,
            surface=self.surface,
            bulk_material=self.bulk_material,
            surface_material=self.surface_material,
            bulk_species=np.asarray(self.bulk_species),
            surface_species=np.asarray(self.surface_species),
            area=self.area,
            perimeter=self.perimeter,
            centroid=self.centroid,
            **{f"pool_{name}": values for name, values in self.pool_totals.items()},
        )
        metadata = {
            "format": "bulk_surface_diffuse_v1",
            "bulk_species": list(self.bulk_species),
            "surface_species": list(self.surface_species),
            "maximum_relative_area_error": float(
                np.max(np.abs(self.area / self.area[0] - 1.0))
            ),
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
        return path


def _neumann_stiffness(shape: tuple[int, int], dx: float) -> sparse.csr_matrix:
    """Positive graph stiffness representing ``-laplacian``."""
    ny, nx = shape
    ids = np.arange(ny * nx, dtype=np.int64).reshape(shape)
    left = np.concatenate((ids[:, :-1].ravel(), ids[:-1, :].ravel()))
    right = np.concatenate((ids[:, 1:].ravel(), ids[1:, :].ravel()))
    conductance = np.full(left.size, 1.0 / (dx * dx))
    rows = np.concatenate((left, right, left, right))
    cols = np.concatenate((left, right, right, left))
    values = np.concatenate((conductance, conductance, -conductance, -conductance))
    return sparse.coo_matrix(
        (values, (rows, cols)), shape=(ny * nx, ny * nx)
    ).tocsr()


def _force_area(phi: np.ndarray, delta: np.ndarray, target: float, dx: float) -> np.ndarray:
    """Shift a phase field along its interface to recover the target area."""
    current = float(np.sum(phi) * dx * dx)
    if np.isclose(current, target, rtol=0.0, atol=1e-13 * max(target, 1.0)):
        return np.clip(phi, 0.0, 1.0)
    direction = 1.0 if current < target else -1.0
    lo, hi = (0.0, direction)
    for _ in range(30):
        trial = np.clip(phi + hi * delta, 0.0, 1.0)
        if (float(np.sum(trial) * dx * dx) - target) * direction >= 0:
            break
        hi *= 2.0
    if hi < lo:
        lo, hi = hi, lo
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        trial = np.clip(phi + mid * delta, 0.0, 1.0)
        if float(np.sum(trial) * dx * dx) < target:
            lo = mid
        else:
            hi = mid
    return np.clip(phi + 0.5 * (lo + hi) * delta, 0.0, 1.0)


class _PhaseFieldStepper:
    def __init__(self, geometry: DiffuseGeometry, mechanics: PhaseFieldMechanics,
                 dt: float):
        self.geometry = geometry
        self.mechanics = mechanics
        self.dt = dt
        coefficient = (
            dt * mechanics.surface_tension * geometry.epsilon**2
            / mechanics.relaxation_time
        )
        operator = sparse.eye(geometry.size, format="csc") + coefficient * (
            _neumann_stiffness(geometry.shape, geometry.dx).tocsc()
        )
        self.solve = factorized(operator)
        fields = geometry.fields()
        self.previous_perimeter = fields.perimeter
        self.rest_perimeter = mechanics.rest_perimeter or fields.perimeter

    def step_from_drive(self, drive: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Advance the interface from an outward-force field on the grid."""
        old = self.geometry.phi
        if self.mechanics.frozen:
            zeros = np.zeros_like(old)
            return old.copy(), zeros, zeros
        fields = self.geometry.fields()
        drive = np.asarray(drive, dtype=float).copy()
        if drive.shape != old.shape:
            raise ValueError("interface drive must match the phase-field grid")
        delta_sum = float(np.sum(fields.delta))
        if self.mechanics.remove_mean_force and delta_sum > 0:
            drive -= float(np.sum(fields.delta * drive) / delta_sum)
        area_error = fields.area - self.geometry.target_area
        drive -= self.mechanics.area_penalty * area_error

        strain = (fields.perimeter - self.rest_perimeter) / self.rest_perimeter
        strain_rate = (
            (fields.perimeter - self.previous_perimeter)
            / (self.dt * self.rest_perimeter)
        )
        cortical_tension = (
            self.mechanics.perimeter_elasticity * strain
            + self.mechanics.perimeter_viscosity * strain_rate
        )
        if cortical_tension:
            drive -= cortical_tension * fields.curvature

        # W(phi)=phi^2(1-phi)^2/2. The Laplacian is treated implicitly.
        well_prime = old * (1.0 - old) * (1.0 - 2.0 * old)
        rhs = old + self.dt / self.mechanics.relaxation_time * (
            -self.mechanics.surface_tension * well_prime
            + self.geometry.epsilon * fields.delta * drive
        )
        new = np.asarray(self.solve(rhs.ravel())).reshape(old.shape)
        new = np.clip(new, 0.0, 1.0)
        if self.mechanics.enforce_area:
            new_fields = geometry_fields(new, self.geometry.x, self.geometry.y,
                                         self.geometry.dx)
            new = _force_area(new, new_fields.delta, self.geometry.target_area,
                              self.geometry.dx)
        phi_t = (new - old) / self.dt
        vx, vy = _extend_interface_velocity(
            fields, phi_t, self.mechanics.velocity_smoothing
        )
        self.previous_perimeter = fields.perimeter
        return new, vx, vy

    def step(self, surface: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Advance from diffuse-surface concentrations (legacy backend)."""
        coefficients = self.mechanics.surface_force_coefficients
        if coefficients.shape != (surface.shape[0],):
            raise ValueError(
                "surface_force_coefficients must have one entry per surface species"
            )
        drive = np.tensordot(coefficients, surface, axes=(0, 0))
        return self.step_from_drive(drive)


def _extend_interface_velocity(
    fields: DiffuseGeometryFields, phi_t: np.ndarray, smoothing: float
) -> tuple[np.ndarray, np.ndarray]:
    """Extend the measured normal interface velocity across the fixed grid."""
    delta = fields.delta
    peak = float(np.max(delta))
    if peak <= 0:
        return np.zeros_like(delta), np.zeros_like(delta)
    band = delta >= 0.25 * peak
    vn = np.zeros_like(delta)
    vn[band] = phi_t[band] / np.maximum(delta[band], 1e-14)
    indices = distance_transform_edt(~band, return_distances=False,
                                     return_indices=True)
    vx = (vn * fields.normal_x)[tuple(indices)]
    vy = (vn * fields.normal_y)[tuple(indices)]
    if smoothing > 0:
        vx = gaussian_filter(vx, smoothing, mode="nearest")
        vy = gaussian_filter(vy, smoothing, mode="nearest")
    return vx, vy


def _advect_once(q: np.ndarray, vx: np.ndarray, vy: np.ndarray,
                 dx: float, dt: float) -> np.ndarray:
    """Conservative first-order upwind advection with closed box faces."""
    result = q.copy()
    ux = 0.5 * (vx[:, :-1] + vx[:, 1:])
    uy = 0.5 * (vy[:-1, :] + vy[1:, :])
    for species in range(q.shape[0]):
        values = q[species]
        flux_x = ux * np.where(ux >= 0.0, values[:, :-1], values[:, 1:])
        flux_y = uy * np.where(uy >= 0.0, values[:-1, :], values[1:, :])
        result[species, :, :-1] -= dt / dx * flux_x
        result[species, :, 1:] += dt / dx * flux_x
        result[species, :-1, :] -= dt / dx * flux_y
        result[species, 1:, :] += dt / dx * flux_y
    return result


def _advect(q: np.ndarray, vx: np.ndarray, vy: np.ndarray,
            dx: float, dt: float) -> np.ndarray:
    speed = float(np.max(np.abs(vx)) + np.max(np.abs(vy)))
    substeps = max(1, int(ceil(speed * dt / (0.45 * dx))))
    result = q
    for _ in range(substeps):
        result = _advect_once(result, vx, vy, dx, dt / substeps)
    return result


def _project_material(q: np.ndarray, weight: np.ndarray, floor: float) -> np.ndarray:
    """Confine advected material to a new diffuse support without mass loss."""
    totals = np.sum(q, axis=(1, 2))
    concentration = q / np.maximum(weight, floor)[None, :, :]
    projected = concentration * weight[None, :, :]
    projected_totals = np.sum(projected, axis=(1, 2))
    for species, total in enumerate(totals):
        if projected_totals[species] > 0:
            projected[species] *= total / projected_totals[species]
    return projected


def _repair_small_negatives(q: np.ndarray, *, maximum_fraction: float = 1e-6) -> np.ndarray:
    """Remove small conservative-scheme undershoots without changing amount."""
    for species in range(q.shape[0]):
        values = q[species]
        if np.min(values) >= 0:
            continue
        total = float(np.sum(values))
        negative = float(-np.sum(values[values < 0]))
        if total <= 0 or negative > maximum_fraction * total + 1e-14:
            raise RuntimeError(
                "diffusion produced material negative mass fraction "
                f"{negative / max(abs(total), 1e-300):.6g}; reduce the timestep"
            )
        values[values < 0] = 0.0
        positive_total = float(np.sum(values))
        values *= total / positive_total
    return q


def _bulk_diffusion_rhs(
    concentration: np.ndarray, phi: np.ndarray, diffusion: np.ndarray, dx: float
) -> np.ndarray:
    rhs = np.zeros_like(concentration)
    invdx = 1.0 / dx
    phi_x = 0.5 * (phi[:, :-1] + phi[:, 1:])
    phi_y = 0.5 * (phi[:-1, :] + phi[1:, :])
    for species, coefficient in enumerate(diffusion):
        if coefficient == 0:
            continue
        c = concentration[species]
        flux_x = coefficient * phi_x * (c[:, 1:] - c[:, :-1]) * invdx
        flux_y = coefficient * phi_y * (c[1:, :] - c[:-1, :]) * invdx
        rhs[species, :, :-1] += flux_x * invdx
        rhs[species, :, 1:] -= flux_x * invdx
        rhs[species, :-1, :] += flux_y * invdx
        rhs[species, 1:, :] -= flux_y * invdx
    return rhs


def _surface_diffusion_rhs(
    concentration: np.ndarray,
    fields: DiffuseGeometryFields,
    diffusion: np.ndarray,
    dx: float,
    normal_regularization: float,
) -> np.ndarray:
    rhs = np.zeros_like(concentration)
    delta = fields.delta
    invdx = 1.0 / dx
    for species, coefficient in enumerate(diffusion):
        if coefficient == 0:
            continue
        u = concentration[species]
        gx, gy = _centered_gradient(u, dx)
        tensor_xx = fields.projection_xx + normal_regularization * (
            1.0 - fields.projection_xx
        )
        tensor_xy = (1.0 - normal_regularization) * fields.projection_xy
        tensor_yy = fields.projection_yy + normal_regularization * (
            1.0 - fields.projection_yy
        )
        vector_x = delta * coefficient * (tensor_xx * gx + tensor_xy * gy)
        vector_y = delta * coefficient * (tensor_xy * gx + tensor_yy * gy)
        flux_x = 0.5 * (vector_x[:, :-1] + vector_x[:, 1:])
        flux_y = 0.5 * (vector_y[:-1, :] + vector_y[1:, :])
        rhs[species, :, :-1] += flux_x * invdx
        rhs[species, :, 1:] -= flux_x * invdx
        rhs[species, :-1, :] += flux_y * invdx
        rhs[species, 1:, :] -= flux_y * invdx
    return rhs


class DiffuseDomainSimulator:
    """Coupled phase-field, bulk, and diffuse-membrane simulator."""

    def __init__(
        self,
        geometry: DiffuseGeometry,
        reactions: ReactionSystem,
        *,
        bulk_diffusion: np.ndarray,
        surface_diffusion: np.ndarray,
        mechanics: PhaseFieldMechanics | None = None,
        conservation_pools: tuple[ConservationPool, ...] = (),
        concentration_floor: float = 1e-6,
        surface_floor_fraction: float = 1e-4,
        normal_regularization: float = 1e-3,
    ):
        self.geometry = geometry
        self.reactions = reactions
        self.bulk_diffusion = np.asarray(bulk_diffusion, dtype=float)
        self.surface_diffusion = np.asarray(surface_diffusion, dtype=float)
        self.mechanics = mechanics or PhaseFieldMechanics(
            surface_force_coefficients=np.zeros(len(reactions.surface_species)),
            frozen=True,
        )
        if self.bulk_diffusion.shape != (len(reactions.bulk_species),):
            raise ValueError("bulk_diffusion has the wrong shape")
        if self.surface_diffusion.shape != (len(reactions.surface_species),):
            raise ValueError("surface_diffusion has the wrong shape")
        if np.any(self.bulk_diffusion < 0) or np.any(self.surface_diffusion < 0):
            raise ValueError("diffusion coefficients must be nonnegative")
        if self.mechanics.surface_force_coefficients.shape != (
            len(reactions.surface_species),
        ):
            raise ValueError("mechanical force coefficients have the wrong shape")
        if concentration_floor <= 0 or surface_floor_fraction <= 0:
            raise ValueError("concentration floors must be positive")
        if not 0 <= normal_regularization <= 1:
            raise ValueError("normal_regularization must lie in [0, 1]")
        self.concentration_floor = concentration_floor
        self.surface_floor_fraction = surface_floor_fraction
        self.normal_regularization = normal_regularization
        self.conservation_pools = conservation_pools
        for pool in conservation_pools:
            pool.validate(reactions)

    def _validate_state(self, initial: DiffuseInitialState | InitialState):
        bulk = np.ascontiguousarray(initial.bulk, dtype=float)
        surface = np.ascontiguousarray(initial.surface, dtype=float)
        expected_bulk = (len(self.reactions.bulk_species), *self.geometry.shape)
        expected_surface = (len(self.reactions.surface_species), *self.geometry.shape)
        if bulk.shape != expected_bulk:
            raise ValueError(f"bulk concentrations must have shape {expected_bulk}")
        if surface.shape != expected_surface:
            raise ValueError(f"surface concentrations must have shape {expected_surface}")
        if (bulk.size and np.min(bulk) < 0) or (surface.size and np.min(surface) < 0):
            raise ValueError("initial concentrations must be nonnegative")
        fields = self.geometry.fields()
        return (
            self.geometry.phi[None, :, :] * bulk,
            fields.delta[None, :, :] * surface,
        )

    def _concentrations(self, q_bulk: np.ndarray, q_surface: np.ndarray,
                        fields: DiffuseGeometryFields):
        bulk = q_bulk / np.maximum(
            self.geometry.phi, self.concentration_floor
        )[None, :, :]
        surface_floor = max(
            float(np.max(fields.delta)) * self.surface_floor_fraction,
            self.concentration_floor,
        )
        surface = q_surface / np.maximum(fields.delta, surface_floor)[None, :, :]
        return bulk, surface, surface_floor

    def _chemistry_rhs(self, time: float, q_bulk: np.ndarray,
                       q_surface: np.ndarray, fields: DiffuseGeometryFields):
        bulk, surface, _ = self._concentrations(q_bulk, q_surface, fields)
        n_grid = self.geometry.size
        bulk_flat = bulk.reshape(bulk.shape[0], n_grid)
        surface_flat = surface.reshape(surface.shape[0], n_grid)
        bulk_local = np.zeros_like(bulk_flat)
        surface_local = np.zeros_like(surface_flat)
        exchange = np.zeros((self.reactions.n_exchange, n_grid), dtype=float)
        stimulus = np.empty(n_grid, dtype=float)
        self.reactions.stimulus.values(
            time, fields.membrane_fraction, stimulus
        )
        self.reactions.bulk_kernel(
            time, bulk_flat, self.reactions.parameters, bulk_local
        )
        self.reactions.surface_kernel(
            time, surface_flat, bulk_flat, stimulus,
            self.reactions.parameters, surface_local, exchange,
        )
        bulk_rate = self.geometry.phi[None, :, :] * bulk_local.reshape(q_bulk.shape)
        surface_rate = fields.delta[None, :, :] * surface_local.reshape(
            q_surface.shape
        )
        if self.reactions.n_exchange:
            bulk_exchange = (
                self.reactions.bulk_stoichiometry @ exchange
            ).reshape(q_bulk.shape)
            surface_exchange = (
                self.reactions.surface_stoichiometry @ exchange
            ).reshape(q_surface.shape)
            bulk_rate += fields.delta[None, :, :] * bulk_exchange
            surface_rate += fields.delta[None, :, :] * surface_exchange
        return bulk_rate, surface_rate

    def _chemistry_step(self, time: float, dt: float, q_bulk: np.ndarray,
                        q_surface: np.ndarray, fields: DiffuseGeometryFields):
        b1, s1 = self._chemistry_rhs(time, q_bulk, q_surface, fields)
        bm = q_bulk + 0.5 * dt * b1
        sm = q_surface + 0.5 * dt * s1
        b2, s2 = self._chemistry_rhs(time + 0.5 * dt, bm, sm, fields)
        return q_bulk + dt * b2, q_surface + dt * s2

    def _diffusion_step(self, dt: float, q_bulk: np.ndarray,
                        q_surface: np.ndarray, fields: DiffuseGeometryFields):
        maximum_diffusion = max(
            float(np.max(self.bulk_diffusion)) if self.bulk_diffusion.size else 0.0,
            float(np.max(self.surface_diffusion)) if self.surface_diffusion.size else 0.0,
        )
        # The tensor surface flux is explicit. Subcycling also keeps the bulk
        # path useful for fast cytosolic species without forcing the mechanics
        # and reactions onto the diffusion CFL timestep.
        substeps = max(
            1,
            int(ceil(maximum_diffusion * dt / (0.20 * self.geometry.dx**2))),
        )
        substep = dt / substeps
        for _ in range(substeps):
            bulk, surface, _ = self._concentrations(q_bulk, q_surface, fields)
            q_bulk = q_bulk + substep * _bulk_diffusion_rhs(
                bulk, self.geometry.phi, self.bulk_diffusion, self.geometry.dx
            )
            q_surface = q_surface + substep * _surface_diffusion_rhs(
                surface, fields, self.surface_diffusion, self.geometry.dx,
                self.normal_regularization,
            )
            q_bulk = _repair_small_negatives(q_bulk)
            q_surface = _repair_small_negatives(q_surface)
        return q_bulk, q_surface

    def _pool_total(self, pool: ConservationPool, q_bulk: np.ndarray,
                    q_surface: np.ndarray) -> float:
        scale = self.geometry.dx**2
        bulk_amount = np.sum(q_bulk, axis=(1, 2)) * scale
        surface_amount = np.sum(q_surface, axis=(1, 2)) * scale
        return float(
            np.asarray(pool.bulk_coefficients) @ bulk_amount
            + np.asarray(pool.surface_coefficients) @ surface_amount
        )

    def run(self, initial: DiffuseInitialState | InitialState,
            config: SimulationConfig) -> DiffuseSimulationResult:
        q_bulk, q_surface = self._validate_state(initial)
        n_steps = round(config.end_time / config.dt)
        save_every = round(config.output_interval / config.dt)
        save_steps = np.arange(0, n_steps + 1, save_every, dtype=int)
        if save_steps[-1] != n_steps:
            save_steps = np.append(save_steps, n_steps)
        n_save = save_steps.size
        shape = (n_save, *self.geometry.shape)
        phi_history = np.empty(shape, dtype=np.float32)
        bulk_history = np.empty((n_save, q_bulk.shape[0], *self.geometry.shape),
                                dtype=np.float32)
        surface_history = np.empty(
            (n_save, q_surface.shape[0], *self.geometry.shape), dtype=np.float32
        )
        bulk_material = np.empty((n_save, q_bulk.shape[0]), dtype=float)
        surface_material = np.empty((n_save, q_surface.shape[0]), dtype=float)
        area = np.empty(n_save, dtype=float)
        perimeter = np.empty(n_save, dtype=float)
        centroid = np.empty((n_save, 2), dtype=float)
        pools = {
            pool.name: np.empty(n_save, dtype=float) for pool in self.conservation_pools
        }

        def save(index: int) -> None:
            fields = self.geometry.fields()
            bulk, surface, _ = self._concentrations(q_bulk, q_surface, fields)
            phi_history[index] = self.geometry.phi
            bulk_history[index] = bulk
            surface_history[index] = surface
            bulk_material[index] = np.sum(q_bulk, axis=(1, 2)) * self.geometry.dx**2
            surface_material[index] = (
                np.sum(q_surface, axis=(1, 2)) * self.geometry.dx**2
            )
            area[index] = fields.area
            perimeter[index] = fields.perimeter
            centroid[index] = (fields.centroid_x, fields.centroid_y)
            for pool in self.conservation_pools:
                pools[pool.name][index] = self._pool_total(
                    pool, q_bulk, q_surface
                )

        save(0)
        phase_stepper = _PhaseFieldStepper(
            self.geometry, self.mechanics, config.dt
        )
        next_save_index = 1
        next_save = save_steps[1] if n_save > 1 else -1
        half = 0.5 * config.dt
        for step in range(1, n_steps + 1):
            time = (step - 1) * config.dt
            fields = self.geometry.fields()
            q_bulk, q_surface = self._chemistry_step(
                time, half, q_bulk, q_surface, fields
            )
            q_bulk, q_surface = self._diffusion_step(
                half, q_bulk, q_surface, fields
            )
            _, surface, _ = self._concentrations(q_bulk, q_surface, fields)
            new_phi, vx, vy = phase_stepper.step(surface)
            bulk_totals = np.sum(q_bulk, axis=(1, 2))
            surface_totals = np.sum(q_surface, axis=(1, 2))
            q_bulk = _advect(
                q_bulk, vx, vy, self.geometry.dx, config.dt
            )
            q_surface = _advect(
                q_surface, vx, vy, self.geometry.dx, config.dt
            )
            self.geometry.phi = new_phi
            new_fields = self.geometry.fields()
            q_bulk = _project_material(
                q_bulk, new_phi, self.concentration_floor
            )
            surface_floor = max(
                float(np.max(new_fields.delta)) * self.surface_floor_fraction,
                self.concentration_floor,
            )
            q_surface = _project_material(q_surface, new_fields.delta, surface_floor)
            # Projection is conservative analytically; restore roundoff-level
            # drift separately for each material field.
            for material, totals in ((q_bulk, bulk_totals),
                                     (q_surface, surface_totals)):
                current = np.sum(material, axis=(1, 2))
                for species, total in enumerate(totals):
                    if current[species] != 0:
                        material[species] *= total / current[species]
            q_bulk, q_surface = self._diffusion_step(
                half, q_bulk, q_surface, new_fields
            )
            q_bulk, q_surface = self._chemistry_step(
                time + half, half, q_bulk, q_surface, new_fields
            )
            if not np.all(np.isfinite(q_bulk)) or not np.all(np.isfinite(q_surface)):
                raise RuntimeError("diffuse-domain state became non-finite")
            bulk_negative = q_bulk.size and np.min(q_bulk) < -1e-8
            surface_negative = q_surface.size and np.min(q_surface) < -1e-8
            if bulk_negative or surface_negative:
                raise RuntimeError(
                    "diffuse-domain state became negative at step "
                    f"{step}: bulk minimum="
                    f"{float(np.min(q_bulk)) if q_bulk.size else 0.0:.6g}, "
                    "surface minimum="
                    f"{float(np.min(q_surface)) if q_surface.size else 0.0:.6g}; "
                    "reduce the timestep"
                )
            if step == next_save:
                save(next_save_index)
                next_save_index += 1
                next_save = (
                    save_steps[next_save_index]
                    if next_save_index < n_save else -1
                )

        return DiffuseSimulationResult(
            times=save_steps.astype(float) * config.dt,
            x=self.geometry.x.copy(),
            y=self.geometry.y.copy(),
            phi=phi_history,
            bulk=bulk_history,
            surface=surface_history,
            bulk_material=bulk_material,
            surface_material=surface_material,
            bulk_species=self.reactions.bulk_species,
            surface_species=self.reactions.surface_species,
            pool_totals=pools,
            area=area,
            perimeter=perimeter,
            centroid=centroid,
        )
