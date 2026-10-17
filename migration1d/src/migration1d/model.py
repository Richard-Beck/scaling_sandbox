"""Flexible reaction--diffusion fields on a Kelvin--Voigt moving interval.

The only mechanical degrees of freedom are the two cell endpoints. A spring
and dashpot connect them, substrate drag acts at each endpoint, and arbitrary
callables map endpoint field samples to applied forces. Interior grid points
follow the affine interpolation of the endpoint motion; this is a consequence
of the endpoint model, not a continuum active-gel approximation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

from .chemistry import Chemistry
from .forces import EndpointForce, EndpointState, ZeroForce
from .single_species.model import StripWavePinningChemistry
from .tension import GlobalTensionFeedback, SpatialTensionFeedback, TensionFeedback

Array = NDArray[np.float64]


@dataclass(frozen=True)
class MechanicsParameters:
    spring_stiffness: float = 2.5
    dashpot_viscosity: float = 0.8
    reference_length: float = 10.0
    left_friction: float = 1.0
    right_friction: float = 1.0
    minimum_length: float = 1.0e-5

    def __post_init__(self) -> None:
        if self.spring_stiffness < 0.0 or self.dashpot_viscosity < 0.0:
            raise ValueError("spring stiffness and dashpot viscosity must be non-negative")
        if self.reference_length <= 0.0 or self.minimum_length <= 0.0:
            raise ValueError("reference and minimum lengths must be positive")
        if self.left_friction <= 0.0 or self.right_friction <= 0.0:
            raise ValueError("endpoint friction coefficients must be positive")


@dataclass(frozen=True)
class SimulationConfig:
    nodes: int = 81
    duration: float = 80.0
    output_interval: float = 0.5
    seed: int = 1
    noise_amplitude: float = 0.03
    rtol: float = 2.0e-5
    atol: float = 1.0e-8
    max_step: float = 0.2
    integrator: str = "BDF"
    mechanics_activation_time: float | None = 0.0
    jacobian_sparsity: bool = False

    def __post_init__(self) -> None:
        if self.nodes < 3:
            raise ValueError("nodes must be at least 3")
        if self.duration < 0.0 or self.output_interval <= 0.0 or self.max_step <= 0.0:
            raise ValueError("duration must be non-negative; time intervals must be positive")
        if self.noise_amplitude < 0.0:
            raise ValueError("noise_amplitude must be non-negative")
        if self.mechanics_activation_time is not None and self.mechanics_activation_time < 0.0:
            raise ValueError("mechanics_activation_time must be non-negative or None")


@dataclass(frozen=True)
class MigrationModel:
    mechanics: MechanicsParameters = field(default_factory=MechanicsParameters)
    chemistry: Chemistry = field(default_factory=StripWavePinningChemistry)
    left_force: EndpointForce = field(default_factory=ZeroForce)
    right_force: EndpointForce = field(default_factory=ZeroForce)
    tension_feedback: TensionFeedback | None = None


@dataclass
class SimulationResult:
    model: MigrationModel
    config: SimulationConfig
    times: Array
    reference: Array
    position: Array
    velocity: Array
    membrane: Array
    cytoplasm: Array
    applied_edge_force: Array
    spring_force: Array
    dashpot_force: Array
    tension_excess: Array | None = None

    @property
    def rear(self) -> Array:
        return self.position[:, 0]

    @property
    def front(self) -> Array:
        return self.position[:, -1]

    @property
    def length(self) -> Array:
        return self.front - self.rear

    @property
    def center(self) -> Array:
        return 0.5 * (self.front + self.rear)

    @property
    def edge_velocity(self) -> Array:
        return self.velocity[:, [0, -1]]

    @property
    def strain(self) -> Array:
        return self.length / self.model.mechanics.reference_length - 1.0

    @property
    def internal_force(self) -> Array:
        return self.spring_force + self.dashpot_force

    @property
    def elastic_tension(self) -> Array:
        """Positive elastic force; compression and dashpot force are excluded."""

        return np.maximum(self.spring_force, 0.0)

    @property
    def tension(self) -> Array | None:
        """Spatial tension ``T0 + theta``, or ``None`` outside spatial mode."""

        if self.tension_excess is None:
            return None
        feedback = self.model.tension_feedback
        if not isinstance(feedback, SpatialTensionFeedback):
            return None
        return feedback.baseline + self.tension_excess

    def compartment_totals(self) -> tuple[Array, Array]:
        weights = np.asarray([control_volumes(x) for x in self.position])
        return (
            np.sum(self.membrane * weights[:, None, :], axis=2),
            np.sum(self.cytoplasm * weights[:, None, :], axis=2),
        )


def control_volumes(position: Array) -> Array:
    intervals = np.diff(position)
    weights = np.empty_like(position)
    weights[0] = 0.5 * intervals[0]
    weights[-1] = 0.5 * intervals[-1]
    weights[1:-1] = 0.5 * (intervals[:-1] + intervals[1:])
    return weights


def _positions(endpoints: Array, nodes: int) -> Array:
    return np.linspace(endpoints[0], endpoints[1], nodes)


def initial_state(model: MigrationModel, config: SimulationConfig) -> Array:
    """Pack endpoint positions and conservative amounts for every field."""
    endpoints = np.asarray([0.0, model.mechanics.reference_length])
    position = _positions(endpoints, config.nodes)
    weights = control_volumes(position)
    species = len(model.chemistry.species)
    profile_factory = getattr(model.chemistry, "initial_profiles", None)
    if profile_factory is None:
        membrane = np.broadcast_to(
            model.chemistry.initial_membrane[:, None], (species, config.nodes)
        ).copy()
        cytoplasm = np.broadcast_to(
            model.chemistry.initial_cytoplasm[:, None], (species, config.nodes)
        ).copy()
    else:
        membrane, cytoplasm = (
            np.asarray(values, dtype=float).copy()
            for values in profile_factory(position)
        )
        expected = (species, config.nodes)
        if membrane.shape != expected or cytoplasm.shape != expected:
            raise ValueError(f"chemistry initial profile arrays must have shape {expected}")
        if not np.isfinite(membrane).all() or not np.isfinite(cytoplasm).all():
            raise ValueError("chemistry initial profiles must be finite")
    if config.noise_amplitude:
        rng = np.random.default_rng(config.seed)
        coordinate = np.linspace(0.0, np.pi, config.nodes)
        long_mode = rng.choice((-1.0, 1.0), size=(species, 1)) * np.cos(coordinate)
        transfer = config.noise_amplitude * membrane * (
            long_mode + 0.15 * rng.standard_normal((species, config.nodes))
        )
        transfer -= np.sum(transfer * weights, axis=1)[:, None] / np.sum(weights)
        transfer = np.maximum(transfer, np.maximum(-0.9 * membrane, -0.9 * cytoplasm))
        membrane += transfer
        cytoplasm -= transfer
    pieces = [endpoints, (membrane * weights).ravel(), (cytoplasm * weights).ravel()]
    if isinstance(model.tension_feedback, SpatialTensionFeedback):
        # Excess tension is a material signal, not a conserved concentration.
        pieces.append(np.zeros(config.nodes))
    return np.concatenate(pieces)


def unpack_state(state: Array, species: int, nodes: int) -> tuple[Array, Array, Array]:
    size = species * nodes
    return (
        state[:2],
        state[2 : 2 + size].reshape(species, nodes),
        state[2 + size : 2 + 2 * size].reshape(species, nodes),
    )


def _chemical_state_size(species: int, nodes: int) -> int:
    return 2 + 2 * species * nodes


def _validate_tension_coupling(model: MigrationModel) -> None:
    if model.tension_feedback is None:
        return
    chemistry = model.chemistry
    if not callable(getattr(chemistry, "tension_deactivation", None)):
        raise TypeError(
            "tension feedback requires chemistry.tension_deactivation(...)"
        )
    if (
        isinstance(model.tension_feedback, SpatialTensionFeedback)
        and not callable(getattr(chemistry, "protrusive_activity", None))
    ):
        raise TypeError(
            "spatial tension feedback requires chemistry.protrusive_activity(...)"
        )


def endpoint_forces(
    model: MigrationModel,
    time: float,
    position: Array,
    membrane: Array,
    cytoplasm: Array,
) -> Array:
    """Evaluate the two user mappings in laboratory-frame sign convention."""
    values = []
    for side, mapping in (("left", model.left_force), ("right", model.right_force)):
        state = EndpointState(
            time=time, side=side, position=position, membrane=membrane,
            cytoplasm=cytoplasm, species=tuple(model.chemistry.species),
        )
        value = float(mapping(state))
        if not np.isfinite(value):
            raise ValueError(f"{side} endpoint force mapping returned {value}")
        values.append(value)
    return np.asarray(values)


def solve_endpoint_mechanics(
    model: MigrationModel,
    time: float,
    position: Array,
    membrane: Array,
    cytoplasm: Array,
) -> tuple[Array, Array, float, float]:
    """Solve endpoint drag balanced by applied and Kelvin--Voigt forces.

    If ``F=(F_left,F_right)`` are lab-frame applied forces and
    ``T=k(L-L0)+eta(V_right-V_left)``, then
    ``gamma_left V_left=F_left+T`` and
    ``gamma_right V_right=F_right-T``.
    """
    p = model.mechanics
    length = position[-1] - position[0]
    if length <= p.minimum_length:
        raise RuntimeError("cell interval collapsed")
    applied = endpoint_forces(model, time, position, membrane, cytoplasm)
    spring = p.spring_stiffness * (length - p.reference_length)
    matrix = np.asarray([
        [p.left_friction + p.dashpot_viscosity, -p.dashpot_viscosity],
        [-p.dashpot_viscosity, p.right_friction + p.dashpot_viscosity],
    ])
    rhs = np.asarray([applied[0] + spring, applied[1] - spring])
    edge_velocity = np.linalg.solve(matrix, rhs)
    velocity = np.linspace(edge_velocity[0], edge_velocity[1], len(position))
    dashpot = p.dashpot_viscosity * (edge_velocity[1] - edge_velocity[0])
    return velocity, applied, float(spring), float(dashpot)


def diffusion_amount_rate(concentration: Array, position: Array, diffusivity: Array) -> Array:
    """Conservative finite-volume diffusion with no flux at both cell edges."""
    flux = -diffusivity[:, None] * np.diff(concentration, axis=1) / np.diff(position)[None, :]
    rate = np.zeros_like(concentration)
    rate[:, 0] = -flux[:, 0]
    rate[:, -1] = flux[:, -1]
    rate[:, 1:-1] = flux[:, :-1] - flux[:, 1:]
    return rate


def mechanics_state(
    model: MigrationModel,
    config: SimulationConfig,
    time: float,
    position: Array,
    membrane: Array,
    cytoplasm: Array,
) -> tuple[Array, Array, float, float]:
    release = config.mechanics_activation_time
    if release is not None and time >= release:
        return solve_endpoint_mechanics(model, time, position, membrane, cytoplasm)
    applied = endpoint_forces(model, time, position, membrane, cytoplasm)
    spring = model.mechanics.spring_stiffness * (
        position[-1] - position[0] - model.mechanics.reference_length
    )
    return np.zeros_like(position), applied, float(spring), 0.0


def _rhs(model: MigrationModel, config: SimulationConfig) -> Callable[[float, Array], Array]:
    species, nodes = len(model.chemistry.species), config.nodes
    spatial_feedback = (
        model.tension_feedback
        if isinstance(model.tension_feedback, SpatialTensionFeedback)
        else None
    )
    chemical_size = _chemical_state_size(species, nodes)

    def evaluate(time: float, state: Array) -> Array:
        endpoints, membrane_amount, cytoplasm_amount = unpack_state(
            state[:chemical_size], species, nodes
        )
        position = _positions(endpoints, nodes)
        weights = control_volumes(position)
        membrane = membrane_amount / weights
        cytoplasm = cytoplasm_amount / weights
        velocity, _, spring, dashpot = mechanics_state(
            model, config, time, position, membrane, cytoplasm
        )
        strain = np.full(
            nodes,
            (endpoints[1] - endpoints[0]) / model.mechanics.reference_length - 1.0,
        )
        stress = np.full(nodes, spring + dashpot)
        membrane_reaction, cytoplasm_reaction = model.chemistry.reaction(
            time, position, membrane, cytoplasm, strain, stress
        )
        expected = (species, nodes)
        if membrane_reaction.shape != expected or cytoplasm_reaction.shape != expected:
            raise ValueError(f"chemistry reaction arrays must have shape {expected}")
        tension_rate = None
        feedback = model.tension_feedback
        if isinstance(feedback, GlobalTensionFeedback):
            added_deactivation = np.full(nodes, feedback.alpha_t * max(0.0, spring))
        elif spatial_feedback is not None:
            tension_excess = state[chemical_size:]
            if tension_excess.shape != (nodes,):
                raise ValueError(f"spatial tension field must have shape {(nodes,)}")
            activity = np.asarray(
                model.chemistry.protrusive_activity(membrane, cytoplasm), dtype=float
            )
            if activity.shape != (nodes,) or not np.isfinite(activity).all():
                raise ValueError(f"protrusive activity must be finite with shape {(nodes,)}")
            tension_diffusion = diffusion_amount_rate(
                tension_excess[None, :], position,
                np.asarray([spatial_feedback.diffusivity]),
            )[0] / weights
            tension_rate = (
                tension_diffusion
                - tension_excess / spatial_feedback.relaxation_time
                + spatial_feedback.generation_rate * np.maximum(activity, 0.0)
            )
            added_deactivation = (
                spatial_feedback.alpha_t * np.maximum(tension_excess, 0.0)
            )
        else:
            added_deactivation = None
        if added_deactivation is not None:
            correction_m, correction_c = model.chemistry.tension_deactivation(
                membrane, cytoplasm, added_deactivation
            )
            correction_m = np.asarray(correction_m, dtype=float)
            correction_c = np.asarray(correction_c, dtype=float)
            if correction_m.shape != expected or correction_c.shape != expected:
                raise ValueError(
                    f"tension deactivation arrays must have shape {expected}"
                )
            membrane_reaction = membrane_reaction + correction_m
            cytoplasm_reaction = cytoplasm_reaction + correction_c
        dm = diffusion_amount_rate(
            membrane, position, model.chemistry.membrane_diffusion
        ) + membrane_reaction * weights
        dc = diffusion_amount_rate(
            cytoplasm, position, model.chemistry.cytoplasm_diffusion
        ) + cytoplasm_reaction * weights
        pieces = [velocity[[0, -1]], dm.ravel(), dc.ravel()]
        if tension_rate is not None:
            pieces.append(tension_rate)
        return np.concatenate(pieces)

    return evaluate


def reaction_diffusion_jacobian_sparsity(
    species: int, nodes: int, spatial_tension: bool = False
):
    """Sparsity pattern for local chemistry and nearest-neighbor diffusion.

    Endpoint rows are dense because force mappings are user-defined. Chemical
    rows also include all edge fields because endpoint-derived stress may enter
    an otherwise local reaction law.
    """
    fields = 2 * species + int(spatial_tension)
    size = 2 + fields * nodes
    pattern = lil_matrix((size, size), dtype=bool)
    pattern[:2, :] = True
    edge_columns = []
    for field in range(fields):
        edge_columns.extend((2 + field * nodes, 2 + field * nodes + nodes - 1))
    for field in range(fields):
        for node in range(nodes):
            row = 2 + field * nodes + node
            pattern[row, :2] = True
            pattern[row, edge_columns] = True
            for coupled_field in range(fields):
                pattern[row, 2 + coupled_field * nodes + node] = True
            for neighbor in (node - 1, node + 1):
                if 0 <= neighbor < nodes:
                    pattern[row, 2 + field * nodes + neighbor] = True
    return pattern.tocsr()


def run_simulation(
    model: MigrationModel | None = None,
    config: SimulationConfig | None = None,
) -> SimulationResult:
    model = model or MigrationModel()
    config = config or SimulationConfig()
    _validate_tension_coupling(model)
    species, nodes = len(model.chemistry.species), config.nodes
    if len(model.chemistry.initial_membrane) != species or len(model.chemistry.initial_cytoplasm) != species:
        raise ValueError("initial concentration arrays must match chemistry.species")
    if len(model.chemistry.membrane_diffusion) != species or len(model.chemistry.cytoplasm_diffusion) != species:
        raise ValueError("diffusion arrays must match chemistry.species")
    times = np.arange(
        0.0, config.duration + 0.5 * config.output_interval, config.output_interval
    )
    if times[-1] < config.duration:
        times = np.append(times, config.duration)
    times[-1] = config.duration
    state0 = initial_state(model, config)
    if config.duration == 0.0:
        state_matrix = state0[:, None]
    else:
        solver_options = {}
        if config.jacobian_sparsity and config.integrator in {"BDF", "Radau"}:
            solver_options["jac_sparsity"] = reaction_diffusion_jacobian_sparsity(
                species, nodes,
                spatial_tension=isinstance(
                    model.tension_feedback, SpatialTensionFeedback
                ),
            )
        solution = solve_ivp(
            _rhs(model, config), (0.0, config.duration), state0,
            method=config.integrator, t_eval=times, rtol=config.rtol,
            atol=config.atol, max_step=config.max_step,
            **solver_options,
        )
        if not solution.success:
            raise RuntimeError(f"integration failed: {solution.message}")
        state_matrix = solution.y
    position = np.empty((len(times), nodes))
    velocity = np.empty_like(position)
    membrane = np.empty((len(times), species, nodes))
    cytoplasm = np.empty_like(membrane)
    applied = np.empty((len(times), 2))
    spring = np.empty(len(times))
    dashpot = np.empty(len(times))
    tension_excess = (
        np.empty((len(times), nodes))
        if isinstance(model.tension_feedback, SpatialTensionFeedback)
        else None
    )
    chemical_size = _chemical_state_size(species, nodes)
    for index, (time, state) in enumerate(zip(times, state_matrix.T)):
        endpoints, ma, ca = unpack_state(state[:chemical_size], species, nodes)
        position[index] = _positions(endpoints, nodes)
        weights = control_volumes(position[index])
        membrane[index], cytoplasm[index] = ma / weights, ca / weights
        velocity[index], applied[index], spring[index], dashpot[index] = mechanics_state(
            model, config, time, position[index], membrane[index], cytoplasm[index]
        )
        if tension_excess is not None:
            tension_excess[index] = state[chemical_size:]
    reference = np.linspace(0.0, model.mechanics.reference_length, nodes)
    return SimulationResult(
        model=model, config=config, times=times, reference=reference,
        position=position, velocity=velocity, membrane=membrane,
        cytoplasm=cytoplasm, applied_edge_force=applied,
        spring_force=spring, dashpot_force=dashpot,
        tension_excess=tension_excess,
    )
