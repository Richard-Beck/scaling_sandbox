"""Discrete spring--dashpot--viscous-cortex mechanics.

This module evolves the pre-continuum mechanical rules.  ``x[i]`` is the
displacement of membrane node i, ``v[i]`` is the velocity of its paired
cortical element, and ``c[i]`` is the integrated cortical displacement used
only to render the material elements.

At every free membrane node,

    F_spring[i] + mu[i] * (v[i] - xdot[i]) = 0.

At every free cortical node,

    F_viscous[i] + F_active[i] + mu[i] * (xdot[i] - v[i]) = 0.

Adding the two balances gives the tridiagonal problem used to recover cortical
velocity.  No continuum PDE is solved in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded


ScalarProtocol = Callable[[float], tuple[float, float]]
VelocityProtocol = Callable[[float], tuple[float, float]]
StressProtocol = Callable[[float, np.ndarray], np.ndarray]


def zero_protocol(_time: float) -> tuple[float, float]:
    """Return zero value and zero time derivative."""
    return 0.0, 0.0


def smooth_step(time: float, *, amplitude: float, ramp_time: float) -> tuple[float, float]:
    """Half-cosine ramp followed by a hold; returns value and derivative."""
    if time <= 0:
        return 0.0, 0.0
    if time >= ramp_time:
        return amplitude, 0.0
    phase = np.pi * time / ramp_time
    return amplitude * 0.5 * (1.0 - np.cos(phase)), amplitude * 0.5 * np.pi / ramp_time * np.sin(phase)


def smooth_pulse(
    time: float,
    *,
    amplitude: float,
    ramp_time: float,
    hold_until: float,
    stop_time: float,
) -> tuple[float, float]:
    """Half-cosine rise, hold, and half-cosine fall."""
    if time <= 0 or time >= stop_time:
        return 0.0, 0.0
    if time < ramp_time:
        return smooth_step(time, amplitude=amplitude, ramp_time=ramp_time)
    if time <= hold_until:
        return amplitude, 0.0
    duration = stop_time - hold_until
    phase = np.pi * (time - hold_until) / duration
    value = amplitude * 0.5 * (1.0 + np.cos(phase))
    derivative = -amplitude * 0.5 * np.pi / duration * np.sin(phase)
    return value, derivative


@dataclass(frozen=True)
class DiscreteParameters:
    nodes: int = 21
    spacing: float = 1.0
    spring_stiffness: float = 1.0
    cortical_viscosity: float = 5.0
    linker_friction: float | np.ndarray = 1.0

    def friction_profile(self) -> np.ndarray:
        profile = np.broadcast_to(np.asarray(self.linker_friction, dtype=float), (self.nodes,)).copy()
        if np.any(profile <= 0):
            raise ValueError("all linker frictions must be positive")
        return profile


@dataclass(frozen=True)
class DiscreteProtocol:
    """Boundary and active-stress inputs for a discrete experiment."""

    membrane_left: ScalarProtocol | None = None
    membrane_right: ScalarProtocol | None = None
    cortex_left: VelocityProtocol = zero_protocol
    cortex_right: VelocityProtocol = zero_protocol
    active_stress: StressProtocol | None = None


@dataclass
class DiscreteTrajectory:
    time: np.ndarray
    reference: np.ndarray
    membrane_displacement: np.ndarray
    cortex_displacement: np.ndarray
    membrane_velocity: np.ndarray
    cortex_velocity: np.ndarray
    spring_tension: np.ndarray
    spring_force: np.ndarray
    linker_force: np.ndarray
    active_stress: np.ndarray
    friction: np.ndarray


def spring_forces(displacement: np.ndarray, stiffness: float) -> tuple[np.ndarray, np.ndarray]:
    """Return link tensions and their equal-and-opposite nodal forces."""
    tension = stiffness * np.diff(displacement)
    force = np.empty_like(displacement)
    force[0] = tension[0]
    force[1:-1] = tension[1:] - tension[:-1]
    force[-1] = -tension[-1]
    return tension, force


def active_nodal_forces(link_stress: np.ndarray, nodes: int) -> np.ndarray:
    """Discrete divergence of link-centred active cortical stress."""
    if link_stress.shape != (nodes - 1,):
        raise ValueError("active stress must have one value per cortical link")
    force = np.empty(nodes)
    force[0] = link_stress[0]
    force[1:-1] = link_stress[1:] - link_stress[:-1]
    force[-1] = -link_stress[-1]
    return force


def _cortex_velocity(
    displacement: np.ndarray,
    time: float,
    parameters: DiscreteParameters,
    protocol: DiscreteProtocol,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve the discrete summed force balance for cortical velocity."""
    n = parameters.nodes
    _, membrane_force = spring_forces(displacement, parameters.spring_stiffness)
    link_coordinate = (np.arange(n - 1) + 0.5) * parameters.spacing
    link_stress = (
        np.zeros(n - 1)
        if protocol.active_stress is None
        else np.asarray(protocol.active_stress(time, link_coordinate), dtype=float)
    )
    active_force = active_nodal_forces(link_stress, n)
    left_velocity, _ = protocol.cortex_left(time)
    right_velocity, _ = protocol.cortex_right(time)

    # eta * (v[i+1] - 2v[i] + v[i-1]) = -(F_mem + F_active)
    interior = n - 2
    eta = parameters.cortical_viscosity
    banded = np.zeros((3, interior))
    banded[0, 1:] = eta
    banded[1, :] = -2.0 * eta
    banded[2, :-1] = eta
    rhs = -(membrane_force[1:-1] + active_force[1:-1])
    rhs[0] -= eta * left_velocity
    rhs[-1] -= eta * right_velocity
    velocity = np.empty(n)
    velocity[0], velocity[-1] = left_velocity, right_velocity
    velocity[1:-1] = solve_banded((1, 1), banded, rhs)
    return velocity, link_stress, membrane_force


def _impose_membrane_values(
    displacement: np.ndarray,
    time: float,
    protocol: DiscreteProtocol,
) -> tuple[np.ndarray, dict[int, float]]:
    displacement = displacement.copy()
    prescribed_rates: dict[int, float] = {}
    if protocol.membrane_left is not None:
        displacement[0], prescribed_rates[0] = protocol.membrane_left(time)
    if protocol.membrane_right is not None:
        displacement[-1], prescribed_rates[len(displacement) - 1] = protocol.membrane_right(time)
    return displacement, prescribed_rates


def simulate_discrete(
    parameters: DiscreteParameters,
    protocol: DiscreteProtocol,
    times,
    *,
    initial_membrane=None,
) -> DiscreteTrajectory:
    """Integrate a discrete experiment and return every mechanical field."""
    times = np.asarray(times, dtype=float)
    if parameters.nodes < 4 or np.any(np.diff(times) <= 0) or times[0] < 0:
        raise ValueError("use at least four nodes and strictly increasing nonnegative times")
    friction = parameters.friction_profile()
    initial_x = np.zeros(parameters.nodes) if initial_membrane is None else np.asarray(initial_membrane, dtype=float).copy()
    if initial_x.shape != (parameters.nodes,):
        raise ValueError("initial_membrane has the wrong shape")
    initial_x, _ = _impose_membrane_values(initial_x, times[0], protocol)
    initial = np.concatenate([initial_x, np.zeros(parameters.nodes)])

    def rhs(time, state):
        x, prescribed_rates = _impose_membrane_values(state[: parameters.nodes], time, protocol)
        velocity, _, membrane_force = _cortex_velocity(x, time, parameters, protocol)
        xdot = velocity + membrane_force / friction
        for index, rate in prescribed_rates.items():
            xdot[index] = rate
        return np.concatenate([xdot, velocity])

    solution = solve_ivp(
        rhs,
        (times[0], times[-1]),
        initial,
        t_eval=times,
        method="BDF",
        rtol=3e-8,
        atol=1e-10,
        max_step=max((times[-1] - times[0]) / 300.0, 1e-4),
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    x = solution.y[: parameters.nodes].T
    cortex_displacement = solution.y[parameters.nodes :].T
    velocity = np.empty_like(x)
    xdot = np.empty_like(x)
    tensions = np.empty((len(times), parameters.nodes - 1))
    membrane_forces = np.empty_like(x)
    linker_forces = np.empty_like(x)
    stresses = np.empty((len(times), parameters.nodes - 1))
    for frame, time in enumerate(times):
        x[frame], prescribed_rates = _impose_membrane_values(x[frame], time, protocol)
        velocity[frame], stresses[frame], membrane_forces[frame] = _cortex_velocity(
            x[frame], time, parameters, protocol
        )
        tensions[frame], _ = spring_forces(x[frame], parameters.spring_stiffness)
        xdot[frame] = velocity[frame] + membrane_forces[frame] / friction
        for index, rate in prescribed_rates.items():
            xdot[frame, index] = rate
        linker_forces[frame] = friction * (velocity[frame] - xdot[frame])

    return DiscreteTrajectory(
        time=times,
        reference=np.arange(parameters.nodes) * parameters.spacing,
        membrane_displacement=x,
        cortex_displacement=cortex_displacement,
        membrane_velocity=xdot,
        cortex_velocity=velocity,
        spring_tension=tensions,
        spring_force=membrane_forces,
        linker_force=linker_forces,
        active_stress=stresses,
        friction=friction,
    )
