"""Fast conservative 1D port of the 2026 De Belly mechanochemical model.

The released Gridap code uses staggered backward-Euler finite-element solves.
This port retains that update structure and the published variables, but uses
cell-centred finite volumes, M-matrix implicit diffusion/upwind transport, and
banded linear algebra.  The choices make positivity and MCA conservation
explicit while avoiding repeated finite-element assembly during screens.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter
from typing import Any

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import least_squares

from .contract import Grid1D, Stimulus, Trajectory, output_times
from .protocols import NoStimulus


def smooth_threshold(values: np.ndarray | float, scale: float, threshold: float):
    """Released smooth switch ``0.5 * (tanh((x-threshold)/scale) + 1)``."""
    return 0.5 * (np.tanh((np.asarray(values) - threshold) / scale) + 1.0)


@dataclass(frozen=True)
class DeBellyParameters:
    """Parameters of the released 1D model in seconds, microns, and pN."""

    linker_friction: float = 100.0
    cortical_viscosity: float = 10_000.0
    membrane_stiffness: float = 100.0
    active_stress: float = 100.0
    mca_on: float = 5.0
    mca_off: float = 1.4
    mca_bound_diffusion: float = 0.001
    mca_unbound_diffusion: float = 0.3
    mca_density: float = 0.1
    rac_deactivation: float = 0.04
    rho_deactivation: float = 0.04
    rac_basal_activation: float = 1.0
    rho_basal_activation: float = 1.0
    rac_mechanical_feedback: float = 0.2
    rho_mechanical_feedback: float = 0.2
    rac_diffusion: float = 0.5
    rho_diffusion: float = 0.5
    mca_rac_threshold: float = 0.068
    mca_switch_scale: float = 0.01
    tension_rho_threshold: float = 5.0
    tension_switch_scale: float = 5.0
    protrusion_speed: float = -0.15
    protrusion_switch_scale: float = 0.01
    protrusion_tension_scale: float = 10.0
    boundary_relaxation: float = 0.98
    rac_to_rho_coupling: float = 1.0
    rho_to_rac_coupling: float = 1.0
    opto_target: str = "rac"
    initialization_mode: str = "algebraic_fixed_point"
    assay_mode: str = "published_one_sided"
    center_boundary_motion: bool = True
    passive_boundary: str = "none"
    input_sampling_mode: str = "cell_centred"

    def __post_init__(self) -> None:
        if self.passive_boundary not in {"none", "left", "right", "local_rac"}:
            raise ValueError("passive_boundary must be none, left, right, or local_rac")
        if self.passive_boundary != "none" and (self.center_boundary_motion or self.assay_mode != "reflection_symmetric_extension"):
            raise ValueError("A passive edge requires uncentered symmetric endpoint rules")
        positive = (
            self.linker_friction,
            self.cortical_viscosity,
            self.membrane_stiffness,
            self.mca_on,
            self.mca_off,
            self.mca_density,
            self.rac_deactivation,
            self.rho_deactivation,
            self.mca_switch_scale,
            self.tension_switch_scale,
            self.protrusion_switch_scale,
            self.protrusion_tension_scale,
        )
        if min(positive) <= 0:
            raise ValueError("positive De Belly parameters must be greater than zero")
        nonnegative = (
            self.active_stress,
            self.mca_bound_diffusion,
            self.mca_unbound_diffusion,
            self.rac_basal_activation,
            self.rho_basal_activation,
            self.rac_mechanical_feedback,
            self.rho_mechanical_feedback,
            self.rac_diffusion,
            self.rho_diffusion,
            self.rac_to_rho_coupling,
            self.rho_to_rac_coupling,
        )
        if min(nonnegative) < 0 or not 0 <= self.boundary_relaxation <= 1:
            raise ValueError("couplings/diffusion must be nonnegative; relaxation is in [0,1]")
        if self.opto_target not in {"rac", "rho"}:
            raise ValueError("opto_target must be 'rac' or 'rho'")
        if self.initialization_mode not in {
            "algebraic_fixed_point",
            "julia_warmup",
        }:
            raise ValueError(
                "initialization_mode must be 'algebraic_fixed_point' or 'julia_warmup'"
            )
        if self.assay_mode not in {
            "published_one_sided",
            "reflection_symmetric_extension",
        }:
            raise ValueError(
                "assay_mode must be 'published_one_sided' or "
                "'reflection_symmetric_extension'"
            )
        if self.input_sampling_mode not in {
            "cell_centred",
            "julia_element_right_edge",
        }:
            raise ValueError(
                "input_sampling_mode must be 'cell_centred' or "
                "'julia_element_right_edge'"
            )


@dataclass(frozen=True)
class DeBelly:
    parameters: DeBellyParameters = DeBellyParameters()

    observable_name = "rac"
    preferred_step = 0.5

    @property
    def name(self) -> str:
        if self.parameters.passive_boundary != "none":
            return "debelly_2026_mechanochemical_passive_" + self.parameters.passive_boundary
        if self.parameters.assay_mode == "published_one_sided":
            return "debelly_2026_mechanochemical_published"
        if not self.parameters.center_boundary_motion:
            return "debelly_2026_mechanochemical_symmetric_uncentered"
        return "debelly_2026_mechanochemical_reflection_symmetric_extension"

    def metadata(self) -> dict[str, Any]:
        return asdict(self.parameters)


def _neumann_transport_matrix(
    cells: int,
    dx: float,
    step: float,
    diffusion: float,
    velocity: np.ndarray,
    boundary_velocity: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    """Banded matrix for backward-Euler diffusion and conservative upwinding.

    Diffusion has homogeneous Neumann conditions. Advective boundary fluxes
    use the transported concentration at the boundary, matching the strong
    divergence term in the released Gridap weak form. In particular, the
    prescribed negative leading-edge velocity carries MCA out of the domain;
    treating that face as zero flux reverses the MCA redistribution.
    """
    matrix = np.zeros((3, cells), dtype=float)
    ratio = step * diffusion / dx**2
    matrix[1] = 1.0 + 2.0 * ratio
    matrix[1, 0] = matrix[1, -1] = 1.0 + ratio
    matrix[0, 1:] = -ratio
    matrix[2, :-1] = -ratio

    face_velocity = 0.5 * (velocity[:-1] + velocity[1:])
    scale = step / dx
    for face, speed in enumerate(face_velocity, start=1):
        left, right = face - 1, face
        coefficient = scale * speed
        if speed >= 0:
            matrix[1, left] += coefficient
            matrix[2, left] -= coefficient
        else:
            matrix[0, right] += coefficient
            matrix[1, right] -= coefficient
    left_velocity, right_velocity = boundary_velocity
    matrix[1, 0] -= scale * left_velocity
    matrix[1, -1] += scale * right_velocity
    return matrix


def _transport(
    values: np.ndarray,
    velocity: np.ndarray,
    diffusion: float,
    step: float,
    dx: float,
    *,
    decay: float = 0.0,
    source: np.ndarray | float = 0.0,
    boundary_velocity: tuple[float, float] = (0.0, 0.0),
) -> np.ndarray:
    matrix = _neumann_transport_matrix(
        len(values), dx, step, diffusion, velocity, boundary_velocity
    )
    matrix[1] += step * decay
    result = solve_banded(
        (1, 1), matrix, values + step * source,
        overwrite_ab=True, check_finite=False,
    )
    if np.min(result) < -1e-10:
        raise FloatingPointError("MCA transport produced a negative concentration")
    return np.maximum(result, 0.0)


def _dirichlet_elliptic(
    diagonal: np.ndarray,
    coefficient: float,
    rhs: np.ndarray,
    dx: float,
    left_value: float,
    right_value: float,
) -> np.ndarray:
    """Solve ``(diag - coefficient*L_D) u = rhs`` at cell centres."""
    cells = len(diagonal)
    ratio = coefficient / dx**2
    matrix = np.zeros((3, cells), dtype=float)
    matrix[1] = diagonal + 2.0 * ratio
    matrix[1, 0] = diagonal[0] + 3.0 * ratio
    matrix[1, -1] = diagonal[-1] + 3.0 * ratio
    matrix[0, 1:] = -ratio
    matrix[2, :-1] = -ratio
    adjusted = rhs.copy()
    adjusted[0] += 2.0 * ratio * left_value
    adjusted[-1] += 2.0 * ratio * right_value
    return solve_banded((1, 1), matrix, adjusted, overwrite_ab=True, check_finite=False)


def _implicit_diffusion_decay(
    values: np.ndarray,
    source: np.ndarray,
    diffusion: float,
    decay: float,
    step: float,
    dx: float,
) -> np.ndarray:
    cells = len(values)
    ratio = step * diffusion / dx**2
    matrix = np.zeros((3, cells), dtype=float)
    matrix[1] = 1.0 + step * decay + 2.0 * ratio
    matrix[1, 0] = matrix[1, -1] = 1.0 + step * decay + ratio
    matrix[0, 1:] = -ratio
    matrix[2, :-1] = -ratio
    result = solve_banded(
        (1, 1),
        matrix,
        values + step * source,
        overwrite_ab=True,
        check_finite=False,
    )
    if np.min(result) < -1e-10:
        raise FloatingPointError("Rac/Rho solve produced a negative concentration")
    return np.maximum(result, 0.0)


def _gradient(values: np.ndarray, dx: float, left: float | None = None, right: float | None = None) -> np.ndarray:
    """Second-order cell-centred gradient with optional boundary-face values."""
    result = np.empty_like(values)
    result[1:-1] = (values[2:] - values[:-2]) / (2.0 * dx)
    if left is None:
        result[0] = (values[1] - values[0]) / dx
    else:
        result[0] = (values[1] - (2.0 * left - values[0])) / (2.0 * dx)
    if right is None:
        result[-1] = (values[-1] - values[-2]) / dx
    else:
        result[-1] = ((2.0 * right - values[-1]) - values[-2]) / (2.0 * dx)
    return result


def _homogeneous_rac_rho(
    parameters: DeBellyParameters, step: float
) -> tuple[float, float]:
    bound = parameters.mca_density * parameters.mca_on / (parameters.mca_on + parameters.mca_off)
    rac_mechanical = (
        parameters.rac_mechanical_feedback
        * parameters.rho_to_rac_coupling
        * 0.5
        * (1.0 - np.tanh((bound - parameters.mca_rac_threshold) / parameters.mca_switch_scale))
    )
    rho_mechanical = (
        parameters.rho_mechanical_feedback
        * parameters.rac_to_rho_coupling
        * float(smooth_threshold(0.0, parameters.tension_switch_scale, parameters.tension_rho_threshold))
    )
    alpha = parameters.rac_basal_activation + rac_mechanical
    beta = parameters.rho_basal_activation + rho_mechanical
    rac, rho = 0.0, 0.62
    if parameters.initialization_mode == "julia_warmup":
        # Gridap performs one initialization solve followed by 250 staggered
        # backward-Euler warm-up solves. This path exists specifically for
        # trajectory-level validation against the released implementation.
        for _ in range(251):
            rac = (
                rac + step * parameters.rac_deactivation * alpha / (1.0 + rho * rho)
            ) / (1.0 + step * parameters.rac_deactivation)
            rho = (
                rho + step * parameters.rho_deactivation * beta / (1.0 + rac * rac)
            ) / (1.0 + step * parameters.rho_deactivation)
    else:
        # Production assays begin at the homogeneous algebraic steady state,
        # making the initial condition independent of the integration step.
        for _ in range(10_000):
            previous_rac, previous_rho = rac, rho
            rac = alpha / (1.0 + rho * rho)
            rho = beta / (1.0 + rac * rac)
            if max(abs(rac - previous_rac), abs(rho - previous_rho)) < 1e-14:
                break
        else:
            # Plain Picard iteration can enter a two-cycle even when the
            # algebraic homogeneous state is well-defined (for example,
            # alpha=beta=2 has the exact solution Rac=Rho=1).  Resolve that
            # numerical failure without changing the requested equilibrium.
            solution = least_squares(
                lambda state: np.asarray(
                    [
                        state[0] - alpha / (1.0 + state[1] ** 2),
                        state[1] - beta / (1.0 + state[0] ** 2),
                    ]
                ),
                x0=np.asarray([max(rac, 1e-8), max(rho, 1e-8)]),
                bounds=(0.0, np.inf),
                xtol=1e-14,
                ftol=1e-14,
                gtol=1e-14,
                max_nfev=10_000,
            )
            residual = np.max(np.abs(solution.fun))
            if not solution.success or residual > 1e-11:
                raise RuntimeError("homogeneous Rac-Rho fixed point did not converge")
            rac, rho = solution.x
    return float(rac), float(rho)


def _derived_rates(
    parameters: DeBellyParameters,
    mca_bound: np.ndarray,
    tension: np.ndarray,
    external: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    delta_alpha_mca, delta_beta_tension = _feedback_terms(
        parameters, mca_bound, tension
    )
    alpha = (
        parameters.rac_basal_activation
        + (external if parameters.opto_target == "rac" else 0.0)
        + delta_alpha_mca
    )
    beta = (
        parameters.rho_basal_activation
        + (external if parameters.opto_target == "rho" else 0.0)
        + delta_beta_tension
    )
    return alpha, beta


def _feedback_terms(
    parameters: DeBellyParameters,
    mca_bound: np.ndarray,
    tension: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the two mechanical increments to Rac and Rho activation."""
    mca_signal = 0.5 * (
        1.0
        - np.tanh(
            (mca_bound - parameters.mca_rac_threshold) / parameters.mca_switch_scale
        )
    )
    tension_signal = smooth_threshold(
        tension, parameters.tension_switch_scale, parameters.tension_rho_threshold
    )
    delta_alpha_mca = (
        parameters.rac_mechanical_feedback
        * parameters.rho_to_rac_coupling
        * mca_signal
    )
    delta_beta_tension = (
        parameters.rho_mechanical_feedback
        * parameters.rac_to_rho_coupling
        * tension_signal
    )
    return delta_alpha_mca, delta_beta_tension


def simulate_debelly(
    model: DeBelly,
    *,
    protocol: Stimulus | None = None,
    length: float = 20.0,
    cells: int = 100,
    end_time: float = 350.0,
    output_interval: float = 1.0,
    times: np.ndarray | None = None,
    step: float | None = None,
    initial_fields: dict[str, np.ndarray] | None = None,
    substrate_clutch=None,
) -> Trajectory:
    """Run the mechanochemical model with staggered conservative banded solves."""
    step = model.preferred_step if step is None else step
    if step <= 0:
        raise ValueError("step must be positive")
    grid = Grid1D(length=length, cells=cells)
    parameters = model.parameters
    if parameters.passive_boundary != "none" and substrate_clutch is None:
        raise ValueError("The passive-edge branch requires the substrate clutch force balance")
    protocol = NoStimulus() if protocol is None else protocol
    if times is None:
        times = output_times(end_time, output_interval)
    else:
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
            raise ValueError("times must start at zero and be strictly increasing")

    equilibrium_fraction = parameters.mca_on / (parameters.mca_on + parameters.mca_off)
    mca_bound = np.full(cells, parameters.mca_density * equilibrium_fraction)
    mca_unbound = np.full(cells, parameters.mca_density * (1.0 - equilibrium_fraction))
    rac0, rho0 = _homogeneous_rac_rho(parameters, step)
    rac = np.full(cells, rac0)
    rho = np.full(cells, rho0)
    if initial_fields is not None:
        unknown = set(initial_fields) - {"rac", "rho"}
        if unknown:
            raise ValueError(
                "initial_fields only accepts 'rac' and 'rho'; got "
                + ", ".join(sorted(unknown))
            )
        for name, target in (("rac", rac), ("rho", rho)):
            if name not in initial_fields:
                continue
            values = np.asarray(initial_fields[name], dtype=float)
            if values.shape != (cells,):
                raise ValueError(
                    f"initial {name} must have shape ({cells},), got {values.shape}"
                )
            if not np.all(np.isfinite(values)) or np.min(values) < 0:
                raise ValueError(f"initial {name} must be finite and nonnegative")
            target[:] = values
    displacement = np.zeros(cells)
    velocity = np.zeros(cells)
    tension = np.zeros(cells)
    raw_left_displacement = 0.0
    raw_right_displacement = 0.0
    left_boundary_displacement = 0.0
    right_boundary_displacement = 0.0
    boundary_diagnostics = {key: 0.0 for key in (
        "passive_stress_residual", "passive_kinematic_residual", "global_force_residual",
        "left_total_boundary_stress", "right_total_boundary_stress",
        "passive_membrane_stress", "passive_cortical_stress",
        "passive_membrane_velocity", "passive_cortical_velocity")}
    maximum_boundary_residuals = {key: 0.0 for key in (
        "passive_stress_residual", "passive_kinematic_residual", "global_force_residual")}
    total_mca = float(np.sum(mca_bound + mca_unbound))
    if parameters.input_sampling_mode == "julia_element_right_edge":
        # The released common-protocol adapter assigns one forcing value per
        # finite element at j*dx.  This mode exists for native trajectory
        # validation; production retains the FV cell-centred representation.
        input_x = np.minimum(grid.x + 0.5 * grid.dx, grid.length)
    else:
        input_x = grid.x

    names = (
        "rac",
        "rho",
        "mca_bound",
        "mca_unbound",
        "displacement",
        "velocity",
        "tension",
        "rac_activation",
        "rho_activation",
        "delta_alpha_mca",
        "delta_beta_tension",
        "left_boundary_displacement",
        "right_boundary_displacement",
        "boundary_separation",
        "mean_strain",
        "tension_range",
    )
    if substrate_clutch is not None:
        if parameters.assay_mode != "reflection_symmetric_extension":
            raise ValueError("Free substrate-clutch motion requires the symmetric assay")
        if substrate_clutch.ids.size:
            raise ValueError("Use a fresh substrate clutch for each simulation")
        substrate_clutch.initialize(np.linspace(0, length, cells + 1))
        names += ("clutch_bound", "clutch_extension", "clutch_traction",
                  "clutch_traction_lab_density", "cell_translation", "cell_velocity",
                  "cortical_velocity_lab")
    storage = {name: np.empty((len(times), cells)) for name in names}
    if parameters.passive_boundary != "none":
        storage.update({key: np.empty((len(times),cells)) for key in boundary_diagnostics})
        boundary_diagnostics["left_total_boundary_stress"] = parameters.active_stress * parameters.rho_to_rac_coupling * rho0
        boundary_diagnostics["right_total_boundary_stress"] = boundary_diagnostics["left_total_boundary_stress"]
        boundary_diagnostics["passive_cortical_stress"] = boundary_diagnostics["left_total_boundary_stress"]

    def record(frame: int, at_time: float) -> None:
        external = np.asarray(protocol(at_time, input_x, grid.length), dtype=float)
        alpha, beta = _derived_rates(parameters, mca_bound, tension, external)
        delta_alpha_mca, delta_beta_tension = _feedback_terms(
            parameters, mca_bound, tension
        )
        separation = right_boundary_displacement - left_boundary_displacement
        scalar_fields = {
            "left_boundary_displacement": left_boundary_displacement,
            "right_boundary_displacement": right_boundary_displacement,
            "boundary_separation": separation,
            "mean_strain": separation / grid.length,
            "tension_range": float(np.max(tension) - np.min(tension)),
        }
        for name, values in (
            ("rac", rac),
            ("rho", rho),
            ("mca_bound", mca_bound),
            ("mca_unbound", mca_unbound),
            ("displacement", displacement),
            ("velocity", velocity),
            ("tension", tension),
            ("rac_activation", alpha),
            ("rho_activation", beta),
            ("delta_alpha_mca", delta_alpha_mca),
            ("delta_beta_tension", delta_beta_tension),
        ):
            storage[name][frame] = values
        for name, value in scalar_fields.items():
            storage[name][frame] = value
        if parameters.passive_boundary != "none":
            for key, value in boundary_diagnostics.items():
                storage[key][frame] = value
        if substrate_clutch is not None:
            for name, values in substrate_clutch.fields.items():
                storage[name][frame] = values
            storage["cell_translation"][frame] = substrate_clutch.translation
            storage["cell_velocity"][frame] = substrate_clutch.speed
            storage["cortical_velocity_lab"][frame] = velocity + substrate_clutch.speed
            substrate_clutch.record(at_time)

    record(0, float(times[0]))
    time = float(times[0])
    accepted_steps = 0
    events = np.asarray(
        sorted(event for event in protocol.event_times if times[0] < event < times[-1]),
        dtype=float,
    )
    start_clock = perf_counter()
    for frame in range(1, len(times)):
        target = float(times[frame])
        while time < target - 1e-13:
            dt = min(step, target - time)
            future = events[events > time + 1e-13]
            if future.size:
                dt = min(dt, float(future[0] - time))
            external = np.asarray(protocol(time, input_x, grid.length), dtype=float)

            # Rac-dependent protrusion. The published assay has only a moving
            # left endpoint. The reflection-symmetric Jilkine extension uses
            # the same tension-inhibited protrusion law at both ends, with
            # optional mean subtraction controlled by center_boundary_motion.
            rac_threshold = 1.3 * rac0
            left_switch = float(
                smooth_threshold(rac[0], parameters.protrusion_switch_scale, rac_threshold)
            )
            left_protrusion_velocity = (
                parameters.rac_to_rho_coupling
                * parameters.protrusion_speed
                * left_switch
                / (1.0 + tension[0] ** 2 / parameters.protrusion_tension_scale)
            )
            if parameters.assay_mode == "reflection_symmetric_extension":
                right_switch = float(
                    smooth_threshold(
                        rac[-1], parameters.protrusion_switch_scale, rac_threshold
                    )
                )
                right_protrusion_velocity = (
                    -parameters.rac_to_rho_coupling
                    * parameters.protrusion_speed
                    * right_switch
                    / (1.0 + tension[-1] ** 2 / parameters.protrusion_tension_scale)
                )
            else:
                right_protrusion_velocity = 0.0

            old_left_boundary_displacement = left_boundary_displacement
            old_right_boundary_displacement = right_boundary_displacement
            # The retained Figure 2 run uses 0.98 per 0.5 s. Preserve that
            # physical relaxation timescale when refining the clutch timestep.
            relaxation = parameters.boundary_relaxation
            if substrate_clutch is not None:
                relaxation = relaxation ** (dt / 0.5)
            raw_left_displacement = (
                relaxation * raw_left_displacement
                + dt * left_protrusion_velocity
            )
            raw_right_displacement = (
                relaxation * raw_right_displacement
                + dt * right_protrusion_velocity
            )
            if (parameters.assay_mode == "reflection_symmetric_extension"
                    and parameters.center_boundary_motion):
                mean_boundary_displacement = 0.5 * (
                    raw_left_displacement + raw_right_displacement
                )
                left_boundary_displacement = (
                    raw_left_displacement - mean_boundary_displacement
                )
                right_boundary_displacement = (
                    raw_right_displacement - mean_boundary_displacement
                )
                mean_protrusion_velocity = 0.5 * (
                    left_protrusion_velocity + right_protrusion_velocity
                )
                left_cortical_boundary_velocity = (
                    left_protrusion_velocity - mean_protrusion_velocity
                )
                right_cortical_boundary_velocity = (
                    right_protrusion_velocity - mean_protrusion_velocity
                )
            else:
                left_boundary_displacement = raw_left_displacement
                right_boundary_displacement = raw_right_displacement
                left_cortical_boundary_velocity = left_protrusion_velocity
                right_cortical_boundary_velocity = right_protrusion_velocity
            left_membrane_velocity = (
                left_boundary_displacement - old_left_boundary_displacement
            ) / dt
            right_membrane_velocity = (
                right_boundary_displacement - old_right_boundary_displacement
            ) / dt

            # Membrane force balance, then material velocity.
            friction_density = parameters.linker_friction * np.maximum(mca_bound, 1e-10)
            diagonal = friction_density / dt
            rhs_x = diagonal * (displacement + dt * velocity)
            old_displacement = displacement
            if parameters.passive_boundary != "none":
                # Solve both mechanical fields and the passive endpoint together.
                # The substrate contact quadrature uses the previous footprint.
                previous_edges = np.linspace(0,length,cells+1) + np.r_[
                    old_left_boundary_displacement,
                    .5*(old_displacement[:-1]+old_displacement[1:]),
                    old_right_boundary_displacement] + substrate_clutch.translation
                passive_right = parameters.passive_boundary == "right"
                if parameters.passive_boundary == "local_rac":
                    displacement, velocity, boundary_u, boundary_v, boundary_diagnostics = substrate_clutch.advance_local_rac_edges(
                        edges=previous_edges, dx=grid.dx, dt=dt, old_u=old_displacement,
                        old_boundaries=(old_left_boundary_displacement, old_right_boundary_displacement),
                        activated=(rac[0] > rac_threshold, rac[-1] > rac_threshold),
                        commanded_velocity=(left_protrusion_velocity, right_protrusion_velocity),
                        friction=friction_density, stiffness=parameters.membrane_stiffness,
                        viscosity=parameters.cortical_viscosity,
                        active_stress=parameters.active_stress*parameters.rho_to_rac_coupling,
                        rho=rho, resting_rho=rho0)
                else:
                    displacement, velocity, boundary_u, boundary_v, boundary_diagnostics = substrate_clutch.advance_passive_edge(
                        edges=previous_edges, dx=grid.dx, dt=dt, old_u=old_displacement,
                        old_boundary=old_right_boundary_displacement if passive_right else old_left_boundary_displacement,
                        active_displacement=left_boundary_displacement if passive_right else right_boundary_displacement,
                        active_velocity=left_cortical_boundary_velocity if passive_right else right_cortical_boundary_velocity,
                        passive_side=parameters.passive_boundary, friction=friction_density,
                        stiffness=parameters.membrane_stiffness, viscosity=parameters.cortical_viscosity,
                        active_stress=parameters.active_stress*parameters.rho_to_rac_coupling,
                        rho=rho, resting_rho=rho0,
                    )
                left_boundary_displacement, right_boundary_displacement = boundary_u
                left_cortical_boundary_velocity, right_cortical_boundary_velocity = boundary_v
                raw_left_displacement, raw_right_displacement = boundary_u
                left_membrane_velocity = (left_boundary_displacement-old_left_boundary_displacement)/dt
                right_membrane_velocity = (right_boundary_displacement-old_right_boundary_displacement)/dt
                for key in maximum_boundary_residuals:
                    maximum_boundary_residuals[key] = max(maximum_boundary_residuals[key],abs(boundary_diagnostics[key]))
                new_edges = np.linspace(0,length,cells+1) + np.r_[
                    left_boundary_displacement,.5*(displacement[:-1]+displacement[1:]),right_boundary_displacement]
                if np.min(np.diff(new_edges)) <= 0:
                    raise FloatingPointError(f"Passive-edge deformation folded the coordinate map at t={time+dt}")
            else:
                displacement = _dirichlet_elliptic(
                    diagonal,
                    parameters.membrane_stiffness,
                    rhs_x,
                    grid.dx,
                    left_boundary_displacement,
                    right_boundary_displacement,
                )
            membrane_velocity = (displacement - old_displacement) / dt
            tension = parameters.membrane_stiffness * _gradient(
                displacement,
                grid.dx,
                left_boundary_displacement,
                right_boundary_displacement,
            )

            # Published staggered MCA update: bound uses the old unbound field,
            # then unbound uses the newly solved bound field.  A final scalar
            # correction removes only the global mass drift of that split solve.
            old_unbound = mca_unbound
            mca_bound = _transport(
                mca_bound,
                velocity,
                parameters.mca_bound_diffusion,
                dt,
                grid.dx,
                decay=parameters.mca_off,
                source=parameters.mca_on * old_unbound,
                boundary_velocity=(
                    left_cortical_boundary_velocity,
                    right_cortical_boundary_velocity,
                ),
            )
            mca_unbound = _transport(
                old_unbound,
                membrane_velocity,
                parameters.mca_unbound_diffusion,
                dt,
                grid.dx,
                decay=parameters.mca_on,
                source=parameters.mca_off * mca_bound,
                boundary_velocity=(left_membrane_velocity, right_membrane_velocity),
            )
            # The released code compensates boundary-advected MCA with a
            # spatially uniform source, partitioned by the resting exchange
            # fractions. Apply the corresponding exact roundoff/mass closure
            # here; multiplicative renormalization spuriously preserves and
            # amplifies the MCA orientation generated at the moving boundary.
            deficit_per_cell = (
                total_mca - float(np.sum(mca_bound + mca_unbound))
            ) / cells
            mca_bound += equilibrium_fraction * deficit_per_cell
            mca_unbound += (1.0 - equilibrium_fraction) * deficit_per_cell
            if min(float(np.min(mca_bound)), float(np.min(mca_unbound))) < -1e-10:
                raise FloatingPointError("MCA mass closure produced a negative concentration")
            mca_bound = np.maximum(mca_bound, 0.0)
            mca_unbound = np.maximum(mca_unbound, 0.0)

            # Cortical force balance driven by the Rho active-stress gradient.
            active_gradient = (
                parameters.active_stress
                * parameters.rho_to_rac_coupling
                * _gradient(rho, grid.dx)
            )
            updated_friction_density = parameters.linker_friction * np.maximum(
                mca_bound, 1e-10
            )
            rhs_velocity = updated_friction_density * membrane_velocity + active_gradient
            if parameters.passive_boundary != "none":
                pass  # Both mechanical fields were already solved together above.
            elif substrate_clutch is None:
                velocity = _dirichlet_elliptic(
                    updated_friction_density,
                    parameters.cortical_viscosity,
                    rhs_velocity,
                    grid.dx,
                    left_cortical_boundary_velocity,
                    right_cortical_boundary_velocity,
                )
            else:
                edge_displacement = np.r_[left_boundary_displacement,
                    0.5 * (displacement[:-1] + displacement[1:]), right_boundary_displacement]
                physical_edges = (np.linspace(0, length, cells+1) + edge_displacement
                                  + substrate_clutch.translation)
                if np.min(np.diff(physical_edges)) <= 0:
                    raise FloatingPointError("Clutch coordinate map folded")
                velocity = substrate_clutch.advance(
                    edges=physical_edges, dx=grid.dx, dt=dt,
                    diagonal=updated_friction_density, viscosity=parameters.cortical_viscosity,
                    rhs=rhs_velocity,
                    boundary=(left_cortical_boundary_velocity, right_cortical_boundary_velocity),
                    solve=_dirichlet_elliptic,
                )

            # Staggered positive backward-Euler Rac and Rho updates.
            alpha, beta = _derived_rates(parameters, mca_bound, tension, external)
            rac_source = parameters.rac_deactivation * alpha / (1.0 + rho * rho)
            rac = _implicit_diffusion_decay(
                rac,
                rac_source,
                parameters.rac_deactivation * parameters.rac_diffusion,
                parameters.rac_deactivation,
                dt,
                grid.dx,
            )
            rho_source = parameters.rho_deactivation * beta / (1.0 + rac * rac)
            rho = _implicit_diffusion_decay(
                rho,
                rho_source,
                parameters.rho_deactivation * parameters.rho_diffusion,
                parameters.rho_deactivation,
                dt,
                grid.dx,
            )
            if not all(
                np.all(np.isfinite(values))
                for values in (rac, rho, mca_bound, mca_unbound, displacement, velocity)
            ):
                raise FloatingPointError("De Belly simulation became non-finite")
            time += dt
            accepted_steps += 1
        record(frame, target)

    elapsed = perf_counter() - start_clock
    initial_mass = float(np.sum(storage["mca_bound"][0] + storage["mca_unbound"][0]))
    mass_series = np.sum(storage["mca_bound"] + storage["mca_unbound"], axis=1)
    mass_error = float(np.max(np.abs(mass_series / initial_mass - 1.0)))
    metadata = {
        "framework": "polarity1d",
        "backend": "finite_volume_banded_staggered_be",
        "source_lineage": (
            "2026 de Belly Rac-Rho-MCA mechanochemical polarity model; "
            "released 1D Gridap staggered backward-Euler implementation; "
            "membrane-cortex mechanics inherited from 2023"
        ),
        "length_um": length,
        "cells": cells,
        "dx_um": grid.dx,
        "fixed_step": step,
        "accepted_steps": accepted_steps,
        "wall_seconds": elapsed,
        "mca_conservation_error": mass_error,
        "custom_initial_fields": (
            sorted(initial_fields) if initial_fields is not None else []
        ),
        "parameters": model.metadata(),
    }
    if substrate_clutch is not None:
        metadata["substrate_clutch"] = asdict(substrate_clutch.parameters)
        metadata["max_net_substrate_force_pN"] = substrate_clutch.max_force_residual
        metadata["boundary_relaxation_reference_step_s"] = 0.5
        if parameters.passive_boundary != "none":
            metadata["backend"] = "finite_volume_coupled_membrane_cortex_banded_be"
            metadata["maximum_boundary_residuals"] = maximum_boundary_residuals
            metadata["resting_total_stress"] = parameters.active_stress * parameters.rho_to_rac_coupling * rho0
            metadata["passive_edge_scheme"] = (
                "Joint backward-Euler membrane/cortex solve with lagged MCA/Rho; "
                "conservative face-stress discretization; inactive edge follows cortex and carries resting total stress. "
                "Active edge keeps prescribed displacement and cortical velocity. "
                "The passive side is specified for the assay and does not switch with Rac."
            )
        metadata["substrate_clutch_scheme"] = (
            "Fixed laboratory adhesion sites; lagged slip-bond off-rate; implicit n and q=n*e; "
            "implicit cortical flow and zero-net-traction translation; first-order footprint update. "
            "New sites unbound; sites leaving contact detach. Forces deposited per reference length."
        )
    if parameters.passive_boundary == "local_rac":
        metadata["passive_edge_scheme"] = (
            "Both ends at resting total stress. Local Rac above 1.3*resting Rac imposes "
            "outward tension-inhibited cortex velocity; otherwise cortex follows membrane endpoint. "
            "Endpoint displacement is solved, without relaxation or centering. Substrate-fixed reference; U=0.")
        metadata["substrate_clutch_scheme"] = (
            "Fixed laboratory adhesion sites; lagged slip-bond off-rate; implicit n and q=n*e; "
            "joint membrane/cortex mechanics in substrate frame; first-order footprint update.")
    return Trajectory(
        model=model.name,
        time=times,
        x=grid.x,
        fields=storage,
        observable_name=model.observable_name,
        metadata=metadata,
    )
