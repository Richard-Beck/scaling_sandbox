"""Fixed 2D bulk geometry coupled to a separate closed 1D membrane.

The bulk uses a phase-field-weighted Cartesian finite-volume discretization.
The membrane is an ordered periodic polygon.  A sparse bilinear trace maps
bulk concentrations to membrane nodes; its mass-weighted transpose is used by
the simulator to spread exchange fluxes back to the bulk conservatively.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, pi, sqrt

import numpy as np
from scipy import sparse


def _graph_laplacian(
    n: int, left: np.ndarray, right: np.ndarray, conductance: np.ndarray
) -> sparse.csr_matrix:
    rows = np.concatenate((left, right, left, right))
    cols = np.concatenate((left, right, right, left))
    values = np.concatenate((conductance, conductance, -conductance, -conductance))
    return sparse.coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()


def _bulk_stiffness(phi: np.ndarray, dx: float) -> sparse.csr_matrix:
    """Return the positive stiffness for ``-div(phi grad(c))``."""
    ny, nx = phi.shape
    ids = np.arange(nx * ny, dtype=np.int64).reshape(ny, nx)
    horizontal_left = ids[:, :-1].ravel()
    horizontal_right = ids[:, 1:].ravel()
    vertical_lower = ids[:-1, :].ravel()
    vertical_upper = ids[1:, :].ravel()
    # Square cells make face length / centre distance equal to one.  Keeping
    # dx in the signature documents the finite-volume construction and makes
    # the mass/stiffness units explicit.
    _ = dx
    horizontal_phi = 0.5 * (phi[:, :-1] + phi[:, 1:]).ravel()
    vertical_phi = 0.5 * (phi[:-1, :] + phi[1:, :]).ravel()
    left = np.concatenate((horizontal_left, vertical_lower))
    right = np.concatenate((horizontal_right, vertical_upper))
    conductance = np.concatenate((horizontal_phi, vertical_phi))
    return _graph_laplacian(nx * ny, left, right, conductance)


def _surface_operators(points: np.ndarray) -> tuple[np.ndarray, sparse.csr_matrix]:
    following = np.roll(points, -1, axis=0)
    edge_lengths = np.linalg.norm(following - points, axis=1)
    if np.any(edge_lengths <= 0):
        raise ValueError("membrane points must be distinct and ordered")
    weights = 0.5 * (edge_lengths + np.roll(edge_lengths, 1))
    left = np.arange(points.shape[0], dtype=np.int64)
    right = np.roll(left, -1)
    stiffness = _graph_laplacian(
        points.shape[0], left, right, 1.0 / edge_lengths
    )
    return weights, stiffness


def _bilinear_trace(x: np.ndarray, y: np.ndarray,
                    points: np.ndarray) -> sparse.csr_matrix:
    nx = x.size
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    for row, (px, py) in enumerate(points):
        ix = int(np.searchsorted(x, px) - 1)
        iy = int(np.searchsorted(y, py) - 1)
        ix = min(max(ix, 0), x.size - 2)
        iy = min(max(iy, 0), y.size - 2)
        tx = (px - x[ix]) / (x[ix + 1] - x[ix])
        ty = (py - y[iy]) / (y[iy + 1] - y[iy])
        node_ids = (
            iy * nx + ix,
            iy * nx + ix + 1,
            (iy + 1) * nx + ix,
            (iy + 1) * nx + ix + 1,
        )
        node_weights = (
            (1.0 - tx) * (1.0 - ty),
            tx * (1.0 - ty),
            (1.0 - tx) * ty,
            tx * ty,
        )
        for node, weight in zip(node_ids, node_weights):
            rows.append(row)
            cols.append(node)
            values.append(weight)
    return sparse.coo_matrix(
        (values, (rows, cols)), shape=(points.shape[0], x.size * y.size)
    ).tocsr()


@dataclass(frozen=True)
class FixedGeometry:
    """Operators and measures for one fixed cell cross-section."""

    x: np.ndarray
    y: np.ndarray
    phi: np.ndarray
    membrane_points: np.ndarray
    membrane_fraction: np.ndarray
    bulk_weights: np.ndarray
    surface_weights: np.ndarray
    bulk_stiffness: sparse.csr_matrix
    surface_stiffness: sparse.csr_matrix
    trace: sparse.csr_matrix

    @classmethod
    def circle(
        cls,
        radius: float = 5.0,
        *,
        dx: float = 0.25,
        interface_width: float | None = None,
        n_surface: int | None = None,
        padding_widths: float = 5.0,
    ) -> "FixedGeometry":
        if radius <= 0 or dx <= 0:
            raise ValueError("radius and dx must be positive")
        epsilon = float(interface_width if interface_width is not None else dx)
        if epsilon <= 0 or padding_widths < 2:
            raise ValueError("interface_width must be positive and padding_widths >= 2")
        half_width = radius + padding_widths * epsilon
        n_half = ceil(half_width / dx)
        x = np.arange(-n_half, n_half + 1, dtype=float) * dx
        y = x.copy()
        xx, yy = np.meshgrid(x, y)
        phi = 0.5 * (
            1.0 - np.tanh((np.hypot(xx, yy) - radius) / (sqrt(2.0) * epsilon))
        )
        # A finite-width tanh interface has a small O(epsilon^2) excess area.
        # Apply one constant quadrature correction so the declared radius and
        # the bulk integration measure agree.  Scaling both mass and stiffness
        # by this same factor leaves the diffusion equation unchanged.
        measure_scale = (pi * radius * radius) / (np.sum(phi) * dx * dx)
        if n_surface is None:
            n_surface = max(32, int(ceil(2.0 * pi * radius / dx)))
        if n_surface < 8:
            raise ValueError("n_surface must be at least 8")
        theta = 2.0 * pi * np.arange(n_surface, dtype=float) / n_surface
        points = np.column_stack((radius * np.cos(theta), radius * np.sin(theta)))
        surface_weights, surface_stiffness = _surface_operators(points)
        bulk_weights = (measure_scale * phi * dx * dx).ravel()
        if np.min(bulk_weights) <= 0:
            raise ValueError("phase field underflowed; reduce padding_widths")
        trace = _bilinear_trace(x, y, points)
        geometry = cls(
            x=x,
            y=y,
            phi=phi,
            membrane_points=points,
            membrane_fraction=theta / (2.0 * pi),
            bulk_weights=bulk_weights,
            surface_weights=surface_weights,
            bulk_stiffness=measure_scale * _bulk_stiffness(phi, dx),
            surface_stiffness=surface_stiffness,
            trace=trace,
        )
        geometry._validate()
        return geometry

    @classmethod
    def ellipse(
        cls,
        semi_major: float,
        semi_minor: float,
        *,
        dx: float = 0.25,
        interface_width: float | None = None,
        n_surface: int | None = None,
        padding_widths: float = 5.0,
    ) -> "FixedGeometry":
        """Return a smooth elliptical bulk with an ordered periodic membrane.

        The phase-field coordinate is the radial ellipse coordinate multiplied
        by the minor semi-axis.  It is exact on the boundary and provides a
        smooth, well-resolved transition without requiring a closed-form
        Euclidean signed distance to an ellipse.
        """
        if semi_major <= 0 or semi_minor <= 0 or dx <= 0:
            raise ValueError("ellipse axes and dx must be positive")
        if semi_minor > semi_major:
            raise ValueError("semi_minor must not exceed semi_major")
        epsilon = float(interface_width if interface_width is not None else dx)
        if epsilon <= 0 or padding_widths < 2:
            raise ValueError("interface_width must be positive and padding_widths >= 2")

        nx_half = ceil((semi_major + padding_widths * epsilon) / dx)
        ny_half = ceil((semi_minor + padding_widths * epsilon) / dx)
        x = np.arange(-nx_half, nx_half + 1, dtype=float) * dx
        y = np.arange(-ny_half, ny_half + 1, dtype=float) * dx
        xx, yy = np.meshgrid(x, y)
        radial_distance = semi_minor * (
            np.sqrt((xx / semi_major) ** 2 + (yy / semi_minor) ** 2) - 1.0
        )
        phi = 0.5 * (1.0 - np.tanh(radial_distance / (sqrt(2.0) * epsilon)))
        target_area = pi * semi_major * semi_minor
        measure_scale = target_area / (np.sum(phi) * dx * dx)

        if n_surface is None:
            # Ramanujan's second approximation is easily accurate enough for
            # choosing a boundary resolution; the actual polygonal perimeter
            # is computed below from its edge lengths.
            h = ((semi_major - semi_minor) / (semi_major + semi_minor)) ** 2
            perimeter = pi * (semi_major + semi_minor) * (
                1.0 + 3.0 * h / (10.0 + sqrt(4.0 - 3.0 * h))
            )
            n_surface = max(32, int(ceil(perimeter / dx)))
        if n_surface < 8:
            raise ValueError("n_surface must be at least 8")

        # Equal-angle nodes are ordered and the surface mass/stiffness
        # operators account for their nonuniform physical spacing.
        theta = 2.0 * pi * np.arange(n_surface, dtype=float) / n_surface
        points = np.column_stack(
            (semi_major * np.cos(theta), semi_minor * np.sin(theta))
        )
        surface_weights, surface_stiffness = _surface_operators(points)
        cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(
            np.diff(np.vstack((points, points[0])), axis=0), axis=1
        ))))
        membrane_fraction = cumulative[:-1] / cumulative[-1]
        bulk_weights = (measure_scale * phi * dx * dx).ravel()
        if np.min(bulk_weights) <= 0:
            raise ValueError("phase field underflowed; reduce padding_widths")
        trace = _bilinear_trace(x, y, points)
        geometry = cls(
            x=x,
            y=y,
            phi=phi,
            membrane_points=points,
            membrane_fraction=membrane_fraction,
            bulk_weights=bulk_weights,
            surface_weights=surface_weights,
            bulk_stiffness=measure_scale * _bulk_stiffness(phi, dx),
            surface_stiffness=surface_stiffness,
            trace=trace,
        )
        geometry._validate()
        return geometry

    @property
    def bulk_shape(self) -> tuple[int, int]:
        return self.phi.shape

    @property
    def n_bulk(self) -> int:
        return self.phi.size

    @property
    def n_surface(self) -> int:
        return self.membrane_points.shape[0]

    @property
    def area(self) -> float:
        return float(np.sum(self.bulk_weights))

    @property
    def perimeter(self) -> float:
        return float(np.sum(self.surface_weights))

    def trace_bulk(self, bulk: np.ndarray) -> np.ndarray:
        bulk = np.asarray(bulk, dtype=float)
        if bulk.shape[-1] != self.n_bulk:
            raise ValueError("bulk array has the wrong final dimension")
        if bulk.ndim == 1:
            return np.asarray(self.trace @ bulk)
        return np.vstack([self.trace @ field for field in bulk])

    def reshape_bulk(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        if values.shape[-1] != self.n_bulk:
            raise ValueError("bulk array has the wrong final dimension")
        return values.reshape(*values.shape[:-1], *self.bulk_shape)

    def _validate(self) -> None:
        ones_bulk = np.ones(self.n_bulk)
        ones_surface = np.ones(self.n_surface)
        if not np.allclose(self.trace @ ones_bulk, ones_surface, atol=1e-14):
            raise ValueError("trace must preserve constant fields")
        if not np.allclose(self.bulk_stiffness @ ones_bulk, 0.0, atol=1e-12):
            raise ValueError("bulk diffusion operator must preserve constants")
        if not np.allclose(
            self.surface_stiffness @ ones_surface, 0.0, atol=1e-12
        ):
            raise ValueError("surface diffusion operator must preserve constants")
