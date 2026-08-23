"""P1 finite-element solver for the Giese bulk--surface wave-pinning model."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import pi
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable

import numpy as np
import meshio
from scipy import sparse
from scipy.sparse.linalg import factorized, spsolve
from skfem import Basis, ElementTriP1, MeshTri, asm
from skfem.models.poisson import laplace, mass

from .metrics import membrane_metrics


@dataclass(frozen=True)
class GieseParameters:
    """Published WP parameters; lengths are micrometres and time is seconds."""

    k0: float = 0.067
    gamma: float = 1.0
    K: float = 0.1
    delta: float = 1.0
    Dm: float = 0.015
    Dc: float = 3.0
    u0: float = 0.026
    v0: float = 0.2


@dataclass(frozen=True)
class Stimulus:
    """A constant membrane stimulus on a circular arc."""

    strength: float
    center_angle: float
    membrane_fraction: float = 0.05
    start_time: float = 0.0
    end_time: float = 10.0

    def values(self, angles: np.ndarray, time: float) -> np.ndarray:
        if time < self.start_time or time >= self.end_time:
            return np.zeros_like(angles)
        difference = np.angle(np.exp(1j * (angles - self.center_angle)))
        half_width = pi * self.membrane_fraction
        return self.strength * (np.abs(difference) <= half_width)


@dataclass
class CircularSystem:
    diameter: float
    target_h: float
    refinement: int
    mesh: MeshTri
    bulk_mass: sparse.csc_matrix
    bulk_stiffness: sparse.csc_matrix
    boundary_nodes: np.ndarray
    boundary_angles: np.ndarray
    edge_lengths: np.ndarray
    membrane_mass: sparse.csc_matrix
    membrane_stiffness: sparse.csc_matrix
    quadrature_angles: np.ndarray
    quadrature_shapes: np.ndarray
    quadrature_weights: np.ndarray

    @property
    def perimeter(self) -> float:
        return float(self.edge_lengths.sum())

    @property
    def area(self) -> float:
        ones = np.ones(self.mesh.p.shape[1])
        return float(ones @ self.bulk_mass @ ones)

    @property
    def max_edge_length(self) -> float:
        facets = self.mesh.facets
        return float(np.linalg.norm(
            self.mesh.p[:, facets[1]] - self.mesh.p[:, facets[0]], axis=0
        ).max())


@dataclass
class SimulationResult:
    system: CircularSystem
    parameters: GieseParameters
    times: np.ndarray
    membrane: np.ndarray
    bulk: np.ndarray
    metrics: list[dict[str, float]]
    total_mass: np.ndarray
    selected_bulk: dict[float, np.ndarray] = field(default_factory=dict)
    selected_membrane: dict[float, np.ndarray] = field(default_factory=dict)
    maximum_step_flux_imbalance: float = 0.0


def _boundary_order(mesh: MeshTri) -> np.ndarray:
    facets = mesh.facets[:, mesh.boundary_facets()]
    nodes = np.unique(facets)
    angles = np.mod(np.arctan2(mesh.p[1, nodes], mesh.p[0, nodes]), 2.0 * pi)
    ordered = nodes[np.argsort(angles)]
    actual_edges = {tuple(sorted(edge)) for edge in facets.T.tolist()}
    expected_edges = {
        tuple(sorted((int(ordered[i]), int(ordered[(i + 1) % ordered.size]))))
        for i in range(ordered.size)
    }
    if actual_edges != expected_edges:
        raise ValueError("Boundary facets do not form one angle-ordered circular cycle")
    return ordered


def _membrane_matrices(edge_lengths: np.ndarray) -> tuple[sparse.csc_matrix,
                                                            sparse.csc_matrix]:
    n = edge_lengths.size
    rows: list[int] = []
    cols: list[int] = []
    mass_data: list[float] = []
    stiffness_data: list[float] = []
    for i, length in enumerate(edge_lengths):
        j = (i + 1) % n
        for row, col, mass_value, stiffness_value in (
            (i, i, length / 3.0, 1.0 / length),
            (i, j, length / 6.0, -1.0 / length),
            (j, i, length / 6.0, -1.0 / length),
            (j, j, length / 3.0, 1.0 / length),
        ):
            rows.append(row)
            cols.append(col)
            mass_data.append(mass_value)
            stiffness_data.append(stiffness_value)
    shape = (n, n)
    return (
        sparse.coo_matrix((mass_data, (rows, cols)), shape=shape).tocsc(),
        sparse.coo_matrix((stiffness_data, (rows, cols)), shape=shape).tocsc(),
    )


def make_circular_system(diameter: float, target_h: float = 0.25,
                         max_refinement: int = 8) -> CircularSystem:
    """Construct a circular triangular mesh at approximately fixed physical h."""
    if diameter <= 0.0 or target_h <= 0.0:
        raise ValueError("diameter and target_h must be positive")
    radius = diameter / 2.0
    refinement = 0
    while refinement <= max_refinement:
        unit = MeshTri.init_circle(nrefs=refinement, smoothed=True)
        mesh = MeshTri(unit.p * radius, unit.t)
        facets = mesh.facets
        max_h = np.linalg.norm(
            mesh.p[:, facets[1]] - mesh.p[:, facets[0]], axis=0
        ).max()
        if max_h <= target_h:
            break
        refinement += 1
    else:
        raise ValueError("target_h requires a refinement above max_refinement")

    return _assemble_circular_system(
        mesh, diameter=diameter, target_h=target_h, refinement=refinement
    )


def make_gmsh_circular_system(diameter: float, boundary_h: float = 0.06,
                              bulk_h: float = 0.18,
                              refinement_distance_fraction: float = 0.3
                              ) -> CircularSystem:
    """Construct an irregular Gmsh disk refined toward the outer membrane."""
    if diameter <= 0.0 or boundary_h <= 0.0 or bulk_h <= 0.0:
        raise ValueError("diameter and mesh sizes must be positive")
    if boundary_h > bulk_h:
        raise ValueError("boundary_h must not exceed bulk_h")
    environment_gmsh = Path(sys.executable).with_name("gmsh")
    gmsh = shutil.which("gmsh") or (
        str(environment_gmsh) if environment_gmsh.is_file() else None
    )
    if gmsh is None:
        raise RuntimeError("Gmsh executable is required for boundary-refined meshes")
    radius = diameter / 2.0
    refinement_distance = radius * refinement_distance_fraction
    geometry = f"""SetFactory(\"OpenCASCADE\");
Disk(1) = {{0, 0, 0, {radius}, {radius}}};
Mesh.Algorithm = 6;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.RandomFactor = 1e-9;
Field[1] = Distance;
Field[1].CurvesList = {{1}};
Field[1].Sampling = 200;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = {boundary_h};
Field[2].SizeMax = {bulk_h};
Field[2].DistMin = 0;
Field[2].DistMax = {refinement_distance};
Background Field = 2;
Mesh.MshFileVersion = 2.2;
"""
    with tempfile.TemporaryDirectory(prefix="giese_gmsh_") as directory:
        directory_path = Path(directory)
        geo_path = directory_path / "disk.geo"
        mesh_path = directory_path / "disk.msh"
        geo_path.write_text(geometry)
        completed = subprocess.run(
            [gmsh, str(geo_path), "-2", "-format", "msh2", "-o", str(mesh_path),
             "-v", "0"],
            check=True, capture_output=True, text=True,
        )
        if completed.stderr.strip():
            raise RuntimeError(f"Gmsh reported an error: {completed.stderr.strip()}")
        generated = meshio.read(mesh_path)
    triangles = generated.cells_dict.get("triangle")
    if triangles is None:
        raise RuntimeError("Gmsh output did not contain first-order triangles")
    mesh = MeshTri(generated.points[:, :2].T, triangles.T).remove_unused_nodes()
    return _assemble_circular_system(
        mesh, diameter=diameter, target_h=bulk_h, refinement=-1
    )


def _assemble_circular_system(mesh: MeshTri, diameter: float, target_h: float,
                              refinement: int) -> CircularSystem:
    basis = Basis(mesh, ElementTriP1())
    bulk_mass = asm(mass, basis).tocsc()
    bulk_stiffness = asm(laplace, basis).tocsc()
    boundary_nodes = _boundary_order(mesh)
    xy = mesh.p[:, boundary_nodes]
    next_xy = np.roll(xy, -1, axis=1)
    edge_lengths = np.linalg.norm(next_xy - xy, axis=0)
    membrane_mass, membrane_stiffness = _membrane_matrices(edge_lengths)
    boundary_angles = np.mod(np.arctan2(xy[1], xy[0]), 2.0 * pi)

    # Two-point Gauss rule on every boundary edge.  The same quadrature load
    # vector is used by both compartments with opposite signs.
    xi = np.array([0.5 - 0.5 / np.sqrt(3.0),
                   0.5 + 0.5 / np.sqrt(3.0)])
    shapes = np.column_stack((1.0 - xi, xi))
    qxy = (xy[:, :, None] * shapes[None, None, :, 0]
           + next_xy[:, :, None] * shapes[None, None, :, 1])
    quadrature_angles = np.mod(np.arctan2(qxy[1], qxy[0]), 2.0 * pi)
    quadrature_weights = edge_lengths[:, None] * 0.5

    return CircularSystem(
        diameter=diameter, target_h=target_h, refinement=refinement, mesh=mesh,
        bulk_mass=bulk_mass, bulk_stiffness=bulk_stiffness,
        boundary_nodes=boundary_nodes, boundary_angles=boundary_angles,
        edge_lengths=edge_lengths, membrane_mass=membrane_mass,
        membrane_stiffness=membrane_stiffness,
        quadrature_angles=quadrature_angles, quadrature_shapes=shapes,
        quadrature_weights=quadrature_weights,
    )


class GieseSolver:
    """First-order IMEX solver for coupled bulk and membrane P1 systems."""

    def __init__(self, system: CircularSystem,
                 parameters: GieseParameters | None = None):
        self.system = system
        self.parameters = parameters or GieseParameters()
        self._factorizations: dict[float, tuple] = {}

    def initial_state(self) -> tuple[np.ndarray, np.ndarray]:
        return (
            np.full(self.system.boundary_nodes.size, self.parameters.u0),
            np.full(self.system.mesh.p.shape[1], self.parameters.v0),
        )

    def total_mass(self, u: np.ndarray, v: np.ndarray) -> float:
        membrane = np.ones(u.size) @ self.system.membrane_mass @ u
        bulk = np.ones(v.size) @ self.system.bulk_mass @ v
        return float(membrane + bulk)

    def reaction_load(self, u: np.ndarray, v: np.ndarray, time: float,
                      stimuli: Iterable[Stimulus] = ()) -> np.ndarray:
        system = self.system
        p = self.parameters
        shape = system.quadrature_shapes
        u_next = np.roll(u, -1)
        v_boundary = v[system.boundary_nodes]
        v_next = np.roll(v_boundary, -1)
        uq = u[:, None] * shape[None, :, 0] + u_next[:, None] * shape[None, :, 1]
        vq = (v_boundary[:, None] * shape[None, :, 0]
              + v_next[:, None] * shape[None, :, 1])
        flux = vq * (p.k0 + p.gamma * uq**2 / (p.K**2 + uq**2)) - p.delta * uq
        weighted = system.quadrature_weights * flux
        load = np.zeros_like(u)
        load += np.sum(weighted * shape[None, :, 0], axis=1)
        np.add.at(load, (np.arange(u.size) + 1) % u.size,
                  np.sum(weighted * shape[None, :, 1], axis=1))
        for item in stimuli:
            load += self.stimulus_reaction_load(v, item, time)
        return load

    def _coefficient_matrix(self, coefficient: np.ndarray) -> sparse.csc_matrix:
        """Assemble an edge matrix from a coefficient at standard quadrature."""
        system = self.system
        shape = system.quadrature_shapes
        weighted = system.quadrature_weights * coefficient
        diagonal_left = np.sum(weighted * shape[None, :, 0] ** 2, axis=1)
        off_diagonal = np.sum(
            weighted * shape[None, :, 0] * shape[None, :, 1], axis=1
        )
        diagonal_right = np.sum(weighted * shape[None, :, 1] ** 2, axis=1)
        left = np.arange(system.boundary_nodes.size)
        right = (left + 1) % left.size
        rows = np.concatenate((left, left, right, right))
        cols = np.concatenate((left, right, left, right))
        data = np.concatenate((diagonal_left, off_diagonal,
                               off_diagonal, diagonal_right))
        return sparse.coo_matrix(
            (data, (rows, cols)), shape=(left.size, left.size)
        ).tocsc()

    def _stimulus_coefficient_matrix(self, stimulus: Stimulus,
                                     time: float) -> sparse.csc_matrix:
        """Exactly integrate a discontinuous stimulus times P1 basis products."""
        system = self.system
        n = system.boundary_nodes.size
        if time < stimulus.start_time or time >= stimulus.end_time:
            return sparse.csc_matrix((n, n))
        perimeter = system.perimeter
        center = (stimulus.center_angle % (2.0 * pi)) / (2.0 * pi) * perimeter
        half_width = 0.5 * stimulus.membrane_fraction * perimeter
        edge_starts = np.concatenate(([0.0], np.cumsum(system.edge_lengths[:-1])))
        gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []
        for edge, (start, length) in enumerate(zip(edge_starts, system.edge_lengths)):
            end = start + length
            local_matrix = np.zeros((2, 2))
            for shift in (-perimeter, 0.0, perimeter):
                left = max(start, center - half_width + shift)
                right = min(end, center + half_width + shift)
                if right <= left:
                    continue
                xi_left = (left - start) / length
                xi_right = (right - start) / length
                xi = 0.5 * (xi_left + xi_right) + 0.5 * (xi_right - xi_left) * gauss
                weights = np.full(2, 0.5 * length * (xi_right - xi_left))
                shapes = np.column_stack((1.0 - xi, xi))
                for row in range(2):
                    for col in range(2):
                        local_matrix[row, col] += stimulus.strength * np.dot(
                            weights, shapes[:, row] * shapes[:, col]
                        )
            next_edge = (edge + 1) % n
            for row, global_row in enumerate((edge, next_edge)):
                for col, global_col in enumerate((edge, next_edge)):
                    if local_matrix[row, col] != 0.0:
                        rows.append(global_row)
                        cols.append(global_col)
                        data.append(local_matrix[row, col])
        return sparse.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsc()

    def _activation_matrix(self, u: np.ndarray, time: float,
                           stimuli: Iterable[Stimulus]) -> sparse.csc_matrix:
        system = self.system
        p = self.parameters
        shape = system.quadrature_shapes
        uq = (u[:, None] * shape[None, :, 0]
              + np.roll(u, -1)[:, None] * shape[None, :, 1])
        activation = p.k0 + p.gamma * uq**2 / (p.K**2 + uq**2)
        matrix = self._coefficient_matrix(activation)
        for stimulus in stimuli:
            matrix += self._stimulus_coefficient_matrix(stimulus, time)
        return matrix

    def _reaction_u_jacobian(self, u: np.ndarray, v: np.ndarray) -> sparse.csc_matrix:
        system = self.system
        p = self.parameters
        shape = system.quadrature_shapes
        uq = (u[:, None] * shape[None, :, 0]
              + np.roll(u, -1)[:, None] * shape[None, :, 1])
        v_boundary = v[system.boundary_nodes]
        vq = (v_boundary[:, None] * shape[None, :, 0]
              + np.roll(v_boundary, -1)[:, None] * shape[None, :, 1])
        derivative = (
            vq * p.gamma * (2.0 * uq * p.K**2) / (p.K**2 + uq**2) ** 2
            - p.delta
        )
        return self._coefficient_matrix(derivative)

    def stimulus_reaction_load(self, v: np.ndarray, stimulus: Stimulus,
                               time: float) -> np.ndarray:
        """Integrate ``v * kS`` exactly over the prescribed membrane fraction.

        The stimulus is discontinuous at its patch ends, so applying the
        background Gauss rule without splitting an intersected edge would make
        its physical coverage mesh-dependent.  Linear trace fields times P1
        test functions are quadratic and are exact under the split two-point
        rule used here.
        """
        system = self.system
        load = np.zeros(system.boundary_nodes.size)
        if time < stimulus.start_time or time >= stimulus.end_time:
            return load
        perimeter = system.perimeter
        center = (stimulus.center_angle % (2.0 * pi)) / (2.0 * pi) * perimeter
        half_width = 0.5 * stimulus.membrane_fraction * perimeter
        edge_starts = np.concatenate(([0.0], np.cumsum(system.edge_lengths[:-1])))
        v_boundary = v[system.boundary_nodes]
        gauss = np.array([-1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0)])
        for edge, (start, length) in enumerate(zip(edge_starts, system.edge_lengths)):
            end = start + length
            for shift in (-perimeter, 0.0, perimeter):
                patch_start = center - half_width + shift
                patch_end = center + half_width + shift
                left = max(start, patch_start)
                right = min(end, patch_end)
                if right <= left:
                    continue
                xi_left = (left - start) / length
                xi_right = (right - start) / length
                xi = 0.5 * (xi_left + xi_right) + 0.5 * (xi_right - xi_left) * gauss
                weights = np.full(2, 0.5 * length * (xi_right - xi_left))
                shapes = np.column_stack((1.0 - xi, xi))
                next_edge = (edge + 1) % load.size
                vq = (v_boundary[edge] * shapes[:, 0]
                      + v_boundary[next_edge] * shapes[:, 1])
                local = stimulus.strength * weights * vq
                load[edge] += np.dot(local, shapes[:, 0])
                load[next_edge] += np.dot(local, shapes[:, 1])
        return load

    def _solvers(self, dt: float):
        key = float(dt)
        if key not in self._factorizations:
            p = self.parameters
            membrane_operator = (
                self.system.membrane_mass + dt * p.Dm * self.system.membrane_stiffness
            ).tocsc()
            bulk_operator = (
                self.system.bulk_mass + dt * p.Dc * self.system.bulk_stiffness
            ).tocsc()
            self._factorizations[key] = (
                factorized(membrane_operator), factorized(bulk_operator)
            )
        return self._factorizations[key]

    def step(self, u: np.ndarray, v: np.ndarray, time: float, dt: float,
             stimuli: Iterable[Stimulus] = (),
             scheme: str = "imex") -> tuple[np.ndarray, np.ndarray, float]:
        if scheme == "giese_semi_implicit":
            return self._semi_implicit_step(u, v, time, dt, tuple(stimuli))
        if scheme != "imex":
            raise ValueError(f"unknown time-integration scheme: {scheme}")
        reaction = self.reaction_load(u, v, time, stimuli)
        membrane_solve, bulk_solve = self._solvers(dt)
        membrane_rhs = self.system.membrane_mass @ u + dt * reaction
        bulk_rhs = self.system.bulk_mass @ v
        bulk_rhs = np.asarray(bulk_rhs).copy()
        bulk_rhs[self.system.boundary_nodes] -= dt * reaction
        next_u = membrane_solve(membrane_rhs)
        next_v = bulk_solve(bulk_rhs)
        if not np.all(np.isfinite(next_u)) or not np.all(np.isfinite(next_v)):
            raise FloatingPointError("non-finite concentration encountered")
        if min(float(next_u.min()), float(next_v.min())) < -1e-8:
            raise FloatingPointError(
                "materially negative concentration; reduce the time step"
            )
        return next_u, next_v, float(reaction.sum())

    def _semi_implicit_step(self, u: np.ndarray, v: np.ndarray, time: float,
                            dt: float, stimuli: tuple[Stimulus, ...]
                            ) -> tuple[np.ndarray, np.ndarray, float]:
        """Giese scheme (48): implicit u with lagged v, then implicit v with lagged u."""
        system = self.system
        p = self.parameters
        membrane_linear = (
            system.membrane_mass + dt * p.Dm * system.membrane_stiffness
        ).tocsc()
        membrane_rhs = system.membrane_mass @ u
        next_u = u.copy()
        for _ in range(20):
            reaction = self.reaction_load(next_u, v, time, stimuli)
            residual = membrane_linear @ next_u - membrane_rhs - dt * reaction
            jacobian = membrane_linear - dt * self._reaction_u_jacobian(next_u, v)
            update = spsolve(jacobian, -residual)
            next_u += update
            if np.linalg.norm(update, ord=np.inf) <= 1e-11 * max(
                    1.0, np.linalg.norm(next_u, ord=np.inf)):
                break
        else:
            raise RuntimeError("semi-implicit membrane Newton iteration did not converge")

        activation = self._activation_matrix(u, time, stimuli)
        boundary_rows, boundary_cols = activation.nonzero()
        activation_values = activation[boundary_rows, boundary_cols].A1
        bulk_reaction = sparse.coo_matrix(
            (activation_values,
             (system.boundary_nodes[boundary_rows],
              system.boundary_nodes[boundary_cols])),
            shape=system.bulk_mass.shape,
        ).tocsc()
        bulk_operator = (
            system.bulk_mass + dt * p.Dc * system.bulk_stiffness
            + dt * bulk_reaction
        ).tocsc()
        bulk_rhs = np.asarray(system.bulk_mass @ v).copy()
        bulk_rhs[system.boundary_nodes] += dt * p.delta * (
            system.membrane_mass @ u
        )
        next_v = spsolve(bulk_operator, bulk_rhs)
        if not np.all(np.isfinite(next_u)) or not np.all(np.isfinite(next_v)):
            raise FloatingPointError("non-finite concentration encountered")
        if min(float(next_u.min()), float(next_v.min())) < -1e-8:
            raise FloatingPointError(
                "materially negative concentration; reduce the time step"
            )
        surface_reaction = self.reaction_load(next_u, v, time, stimuli)
        return next_u, next_v, float(surface_reaction.sum())

    def run(self, end_time: float, dt: float = 0.05,
            stimuli: Iterable[Stimulus] = (), output_interval: float = 1.0,
            selected_times: Iterable[float] = (),
            record_bulk: bool = False,
            scheme: str = "imex",
            initial: tuple[np.ndarray, np.ndarray] | None = None) -> SimulationResult:
        if end_time < 0.0 or dt <= 0.0 or output_interval <= 0.0:
            raise ValueError("end_time must be nonnegative; dt and output_interval positive")
        stimuli = tuple(stimuli)
        u, v = self.initial_state() if initial is None else (
            np.asarray(initial[0], dtype=float).copy(),
            np.asarray(initial[1], dtype=float).copy(),
        )
        selected = sorted({float(t) for t in selected_times if 0.0 <= t <= end_time})
        record_targets = np.arange(0.0, end_time + 0.5 * output_interval,
                                   output_interval)
        if record_targets.size == 0 or record_targets[-1] < end_time - 1e-12:
            record_targets = np.append(record_targets, end_time)
        event_times = sorted(set(record_targets.tolist() + selected + [end_time]))

        times: list[float] = []
        membrane: list[np.ndarray] = []
        bulk: list[np.ndarray] = []
        metrics: list[dict[str, float]] = []
        masses: list[float] = []
        selected_bulk: dict[float, np.ndarray] = {}
        selected_membrane: dict[float, np.ndarray] = {}

        def save(time: float) -> None:
            if np.any(np.isclose(time, record_targets, atol=1e-10, rtol=0.0)):
                times.append(time)
                membrane.append(u.copy())
                if record_bulk:
                    bulk.append(v.copy())
                metrics.append(membrane_metrics(u, self.system.edge_lengths))
                masses.append(self.total_mass(u, v))
            for target in selected:
                if abs(time - target) <= 1e-10:
                    selected_bulk[target] = v.copy()
                    selected_membrane[target] = u.copy()

        time = 0.0
        save(time)
        maximum_imbalance = 0.0
        for event in event_times[1:]:
            while time < event - 1e-12:
                step_dt = min(dt, event - time)
                old_mass = self.total_mass(u, v)
                u, v, reaction_integral = self.step(
                    u, v, time, step_dt, stimuli, scheme=scheme
                )
                time += step_dt
                new_mass = self.total_mass(u, v)
                maximum_imbalance = max(maximum_imbalance, abs(new_mass - old_mass))
                # Retain this explicit assertion: the shared reaction integral
                # must be gained by the membrane and lost by the bulk.
                if not np.isfinite(reaction_integral):
                    raise FloatingPointError("non-finite integrated reaction flux")
            time = event
            save(time)

        return SimulationResult(
            system=self.system, parameters=self.parameters,
            times=np.asarray(times), membrane=np.asarray(membrane),
            bulk=(np.asarray(bulk) if bulk else
                  np.empty((0, self.system.mesh.p.shape[1]))),
            metrics=metrics, total_mass=np.asarray(masses),
            selected_bulk=selected_bulk, selected_membrane=selected_membrane,
            maximum_step_flux_imbalance=maximum_imbalance,
        )
