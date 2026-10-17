"""Explicit membrane--cytosol Holmes Rac/Cdc42/Rho circuit."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

from .core import ConservationPool, InitialState, ReactionSystem, Simulator
from .geometry import FixedGeometry


FAMILIES = ("cdc42", "rac", "rho")
ACTIVE_SPECIES = tuple(f"{name}_active" for name in FAMILIES)
MEMBRANE_INACTIVE_SPECIES = tuple(f"{name}_inactive_membrane" for name in FAMILIES)
BULK_SPECIES = tuple(f"{name}_inactive_cytosol" for name in FAMILIES)
SURFACE_SPECIES = ACTIVE_SPECIES + MEMBRANE_INACTIVE_SPECIES


@dataclass(frozen=True)
class HolmesParameters:
    cdc42_total: float = 2.4
    rac_total: float = 7.5
    rho_total: float = 3.1
    cdc42_activation: float = 2.95
    rac_basal_activation: float = 0.2
    rho_activation: float = 6.6
    rho_inhibits_cdc42_half_level: float = 1.25
    rac_inhibits_rho_half_level: float = 1.0
    hill: float = 3.0
    cdc42_to_rac_feedback: float = 0.55
    membrane_diffusion: float = 0.1
    effective_inactive_diffusion: float = 50.0
    baseline_membrane_fraction: float = 0.5
    dissociation_rate: float = 1.0
    signal_mean: float = 0.3
    signal_half_range: float = 0.1

    @property
    def totals(self) -> np.ndarray:
        return np.asarray(
            [self.cdc42_total, self.rac_total, self.rho_total], dtype=float
        )

    @property
    def cytosolic_diffusion(self) -> float:
        c0 = self.baseline_membrane_fraction
        if not 0.0 < c0 < 1.0:
            raise ValueError("baseline_membrane_fraction must be in (0, 1)")
        return (
            self.effective_inactive_diffusion - c0 * self.membrane_diffusion
        ) / (1.0 - c0)

    def as_dict(self) -> dict[str, float]:
        values = asdict(self)
        values["cytosolic_diffusion"] = self.cytosolic_diffusion
        return values


@dataclass(frozen=True)
class CartesianGradientStimulus:
    """A static Rac-GEF gradient evaluated at membrane coordinates."""

    values_at_membrane: np.ndarray

    @classmethod
    def for_geometry(
        cls, geometry: FixedGeometry, *, mean: float, half_range: float
    ) -> "CartesianGradientStimulus":
        x = geometry.membrane_points[:, 0]
        scale = float(np.max(np.abs(x)))
        if scale <= 0.0:
            raise ValueError("geometry has no x extent")
        return cls(np.ascontiguousarray(mean + half_range * x / scale))

    def values(self, time: float, membrane_fraction: np.ndarray,
               out: np.ndarray) -> None:
        _ = time, membrane_fraction
        if out.shape != self.values_at_membrane.shape:
            raise ValueError("stimulus and membrane sizes differ")
        out[:] = self.values_at_membrane


@njit(cache=True)
def holmes_surface_kernel(
    time, surface, bulk_trace, stimulus, parameters,
    surface_local_out, exchange_out,
):
    """Activation/inactivation plus inactive membrane--cytosol cycling."""
    _ = time
    surface_local_out.fill(0.0)
    exchange_out.fill(0.0)
    (ct, rt, rhot, ic, ir, irho, a1, a2, hill, feedback,
     c0, association, dissociation) = parameters
    for node in range(surface.shape[1]):
        cdc42 = surface[0, node]
        rac = surface[1, node]
        rho = surface[2, node]
        cdc42_i_m = surface[3, node]
        rac_i_m = surface[4, node]
        rho_i_m = surface[5, node]

        i_cdc42 = ic / (1.0 + (rho / a1) ** hill)
        i_rac = ir + feedback * cdc42 + stimulus[node]
        i_rho = irho / (1.0 + (rac / a2) ** hill)
        f_cdc42 = i_cdc42 * (cdc42_i_m / c0) / ct - cdc42
        f_rac = i_rac * (rac_i_m / c0) / rt - rac
        f_rho = i_rho * (rho_i_m / c0) / rhot - rho

        surface_local_out[0, node] = f_cdc42
        surface_local_out[1, node] = f_rac
        surface_local_out[2, node] = f_rho
        surface_local_out[3, node] = -f_cdc42
        surface_local_out[4, node] = -f_rac
        surface_local_out[5, node] = -f_rho
        exchange_out[0, node] = association * bulk_trace[0, node] - (
            dissociation * cdc42_i_m
        )
        exchange_out[1, node] = association * bulk_trace[1, node] - (
            dissociation * rac_i_m
        )
        exchange_out[2, node] = association * bulk_trace[2, node] - (
            dissociation * rho_i_m
        )


@dataclass(frozen=True)
class HolmesCase:
    simulator: Simulator
    initial: InitialState
    parameters: HolmesParameters
    association_rate: float
    equilibrium_membrane_fraction: float
    projected_totals: np.ndarray
    abundance_rule: str


def _association_rate(
    baseline_geometry: FixedGeometry, parameters: HolmesParameters
) -> float:
    c0 = parameters.baseline_membrane_fraction
    return (
        parameters.dissociation_rate
        * c0 / (1.0 - c0)
        * baseline_geometry.area / baseline_geometry.perimeter
    )


def _reaction(active: np.ndarray, membrane_inactive: np.ndarray,
              signal: float, parameters: HolmesParameters) -> np.ndarray:
    cdc42, rac, rho = active
    c_i, rac_i, rho_i = membrane_inactive
    c0 = parameters.baseline_membrane_fraction
    return np.asarray([
        parameters.cdc42_activation
        / (1.0 + (rho / parameters.rho_inhibits_cdc42_half_level) ** parameters.hill)
        * (c_i / c0) / parameters.cdc42_total - cdc42,
        (parameters.rac_basal_activation
         + parameters.cdc42_to_rac_feedback * cdc42 + signal)
        * (rac_i / c0) / parameters.rac_total - rac,
        parameters.rho_activation
        / (1.0 + (rac / parameters.rac_inhibits_rho_half_level) ** parameters.hill)
        * (rho_i / c0) / parameters.rho_total - rho,
    ])


def homogeneous_equilibrium(
    geometry: FixedGeometry,
    parameters: HolmesParameters,
    association_rate: float,
    *,
    density_reference: InitialState | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Relax the homogeneous direct model before applying the Rac cue.

    With ``density_reference``, its mean membrane densities and bulk
    concentrations are retained as the pre-relaxation abundance rule. This is
    used for the enlarged condition.
    """
    load = geometry.area / geometry.perimeter
    c0 = parameters.baseline_membrane_fraction
    if density_reference is None:
        active0 = np.full(3, 0.1)
        inactive = parameters.totals - active0
        membrane0 = c0 * inactive
        bulk0 = (1.0 - c0) * inactive / load
    else:
        active0 = np.mean(density_reference.surface[:3], axis=1)
        membrane0 = np.mean(density_reference.surface[3:], axis=1)
        bulk0 = np.mean(density_reference.bulk, axis=1)
    projected_totals = active0 + membrane0 + load * bulk0
    state0 = np.concatenate((active0, membrane0, bulk0))

    def rhs(time: float, state: np.ndarray) -> np.ndarray:
        _ = time
        active = state[:3]
        membrane = state[3:6]
        bulk = state[6:9]
        reaction = _reaction(active, membrane, 0.0, parameters)
        exchange = association_rate * bulk - parameters.dissociation_rate * membrane
        return np.concatenate((reaction, -reaction + exchange, -exchange / load))

    result = solve_ivp(
        rhs, (0.0, 600.0), state0, method="BDF",
        rtol=1e-10, atol=1e-12,
    )
    if not result.success:
        raise RuntimeError(f"homogeneous Holmes relaxation failed: {result.message}")
    relaxed = result.y[:, -1]
    if np.min(relaxed) < -1e-10:
        raise RuntimeError("homogeneous Holmes relaxation produced a negative state")
    return relaxed[:3], relaxed[3:6], relaxed[6:9], projected_totals


def holmes_direct_case(
    geometry: FixedGeometry,
    *,
    baseline_geometry: FixedGeometry | None = None,
    parameters: HolmesParameters | None = None,
    density_reference: InitialState | None = None,
    abundance_rule: str = "projected_total_density",
) -> HolmesCase:
    """Construct the explicit nine-field Holmes bulk--surface system."""
    p = parameters or HolmesParameters()
    reference_geometry = baseline_geometry or geometry
    association = _association_rate(reference_geometry, p)
    fraction = association / (
        association + p.dissociation_rate * geometry.area / geometry.perimeter
    )
    active, membrane, bulk, totals = homogeneous_equilibrium(
        geometry, p, association, density_reference=density_reference
    )
    initial = InitialState(
        bulk=np.repeat(bulk[:, None], geometry.n_bulk, axis=1),
        surface=np.repeat(
            np.concatenate((active, membrane))[:, None], geometry.n_surface, axis=1
        ),
    )

    parameters_array = np.asarray([
        p.cdc42_total, p.rac_total, p.rho_total,
        p.cdc42_activation, p.rac_basal_activation, p.rho_activation,
        p.rho_inhibits_cdc42_half_level, p.rac_inhibits_rho_half_level,
        p.hill, p.cdc42_to_rac_feedback, p.baseline_membrane_fraction,
        association, p.dissociation_rate,
    ])
    bulk_stoichiometry = -np.eye(3)
    surface_stoichiometry = np.zeros((6, 3))
    surface_stoichiometry[3:, :] = np.eye(3)
    stimulus = CartesianGradientStimulus.for_geometry(
        geometry, mean=p.signal_mean, half_range=p.signal_half_range
    )
    reactions = ReactionSystem(
        bulk_species=BULK_SPECIES,
        surface_species=SURFACE_SPECIES,
        parameters=parameters_array,
        bulk_stoichiometry=bulk_stoichiometry,
        surface_stoichiometry=surface_stoichiometry,
        surface_kernel=holmes_surface_kernel,
        stimulus=stimulus,
    )
    pools = []
    for index, family in enumerate(FAMILIES):
        bulk_coefficients = np.zeros(3)
        surface_coefficients = np.zeros(6)
        bulk_coefficients[index] = 1.0
        surface_coefficients[index] = 1.0
        surface_coefficients[index + 3] = 1.0
        pools.append(ConservationPool(
            f"{family}_total", bulk_coefficients, surface_coefficients
        ))
    simulator = Simulator(
        geometry, reactions,
        bulk_diffusion=np.full(3, p.cytosolic_diffusion),
        surface_diffusion=np.full(6, p.membrane_diffusion),
        conservation_pools=tuple(pools),
    )
    return HolmesCase(
        simulator=simulator,
        initial=initial,
        parameters=p,
        association_rate=association,
        equilibrium_membrane_fraction=fraction,
        projected_totals=totals,
        abundance_rule=abundance_rule,
    )


def baseline_holmes_geometry(*, dx: float = 0.2,
                             n_surface: int | None = None) -> FixedGeometry:
    return FixedGeometry.ellipse(
        10.0, 4.0, dx=dx, n_surface=n_surface
    )


def enlarged_holmes_geometry(*, dx: float = 0.2,
                             n_surface: int | None = None) -> FixedGeometry:
    scale = 2.0 ** (1.0 / 3.0)
    return FixedGeometry.ellipse(
        10.0 * scale, 4.0 * scale, dx=dx, n_surface=n_surface
    )
