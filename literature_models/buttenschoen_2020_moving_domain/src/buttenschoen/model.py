"""Full 1D moving-domain model from Buttenschön, Liu & Edelstein-Keshet.

The cell is represented on the fixed material coordinate ``xbar in [-1, 1]``.
Active and inactive GTPase concentrations are coupled to Kelvin--Voigt endpoint
mechanics through Eqs. 12a--e of Buttenschön et al. (2020).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

Coupling = Literal["protrusive", "contractile"]


@dataclass(frozen=True)
class ChemistryParameters:
    """Local wave-pinning chemistry in the paper's native units."""

    basal_activation: float = 0.1       # beta, s^-1
    feedback_activation: float = 25.0   # gamma, s^-1
    inactivation: float = 1.0           # I, s^-1
    hill_exponent: float = 4.0          # n
    active_diffusion: float = 0.1       # D, um^2 s^-1
    inactive_diffusion: float = 100.0   # Di, um^2 s^-1


@dataclass(frozen=True)
class MechanicalParameters:
    """Kelvin--Voigt and GTPase-force parameters.

    ``elastic_modulus`` is the 1D spring coefficient E, ``cell_viscosity`` is
    the internal dashpot eta, and ``endpoint_friction`` is the external drag
    mu. The coupling selects outward Rac-like or inward Rho-like forces.
    """

    coupling: Coupling = "protrusive"
    rest_length: float = 50.0           # L0, um
    elastic_modulus: float = 1250.0     # E, pN um^-1
    cell_viscosity: float = 60.0        # eta, pN s um^-1
    endpoint_friction: float = 100000.0 # mu, pN s um^-1
    maximum_force: float = 50000.0      # Fmax, pN
    force_half_level: float = 1.25      # Gc, concentration
    force_hill_exponent: float = 4.0    # m


@dataclass(frozen=True)
class InitialCondition:
    """Edge-biased active-GTPase step with a prescribed physical total."""

    total_amount: float = 25.0          # integral(G + Gi) dx, concentration*um
    active_low: float = 0.01
    active_high: float = 1.5
    active_edge_fraction: float = 0.02
    active_edge: Literal["left", "right"] = "right"
    centroid: float = 0.0


@dataclass(frozen=True)
class SimulationParameters:
    n_cells: int = 161
    end_time: float = 360.0
    output_interval: float = 1.0
    rtol: float = 1e-5
    atol: float = 1e-7
    method: Literal["BDF", "Radau"] = "BDF"
    minimum_radius: float = 1e-3


@dataclass(frozen=True)
class SimulationResult:
    times: np.ndarray
    xbar: np.ndarray
    active: np.ndarray
    inactive: np.ndarray
    centroid: np.ndarray
    radius: np.ndarray
    chemistry: ChemistryParameters
    mechanics: MechanicalParameters
    initial_condition: InitialCondition
    simulation: SimulationParameters
    solver_statistics: dict[str, int | str]

    @property
    def left_edge(self) -> np.ndarray:
        return self.centroid - self.radius

    @property
    def right_edge(self) -> np.ndarray:
        return self.centroid + self.radius

    @property
    def length(self) -> np.ndarray:
        return 2.0 * self.radius

    @property
    def physical_positions(self) -> np.ndarray:
        return self.centroid[:, None] + self.radius[:, None] * self.xbar

    @property
    def total_amount(self) -> np.ndarray:
        dxbar = 2.0 / self.xbar.size
        return self.radius * dxbar * np.sum(self.active + self.inactive, axis=1)

    @property
    def relative_mass_error(self) -> np.ndarray:
        mass = self.total_amount
        return np.abs(mass / mass[0] - 1.0)

    def save(self, output: str | Path) -> Path:
        """Save arrays and a readable parameter/solver metadata record."""
        directory = Path(output)
        directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            directory / "trajectory.npz",
            times=self.times,
            xbar=self.xbar,
            active=self.active,
            inactive=self.inactive,
            centroid=self.centroid,
            radius=self.radius,
            total_amount=self.total_amount,
        )
        metadata = {
            "chemistry": asdict(self.chemistry),
            "mechanics": asdict(self.mechanics),
            "initial_condition": asdict(self.initial_condition),
            "simulation": asdict(self.simulation),
            "solver_statistics": self.solver_statistics,
            "maximum_relative_mass_error": float(np.max(self.relative_mass_error)),
            "minimum_length": float(np.min(self.length)),
            "maximum_length": float(np.max(self.length)),
            "final_length": float(self.length[-1]),
            "final_centroid": float(self.centroid[-1]),
            "final_active_minimum": float(np.min(self.active[-1])),
            "final_active_maximum": float(np.max(self.active[-1])),
            "final_active_contrast": float(
                (np.max(self.active[-1])-np.min(self.active[-1]))
                / max(np.max(self.active[-1])+np.min(self.active[-1]), 1e-15)
            ),
        }
        path = directory / "metadata.json"
        path.write_text(json.dumps(metadata, indent=2) + "\n")
        return path

    @classmethod
    def load(cls, directory: str | Path) -> "SimulationResult":
        """Load a trajectory and its explicit parameter metadata."""
        directory = Path(directory)
        metadata = json.loads((directory/"metadata.json").read_text())
        with np.load(directory/"trajectory.npz", allow_pickle=False) as raw:
            arrays = {name: raw[name].copy() for name in (
                "times", "xbar", "active", "inactive", "centroid", "radius"
            )}
        return cls(
            **arrays,
            chemistry=ChemistryParameters(**metadata["chemistry"]),
            mechanics=MechanicalParameters(**metadata["mechanics"]),
            initial_condition=InitialCondition(**metadata["initial_condition"]),
            simulation=SimulationParameters(**metadata["simulation"]),
            solver_statistics=metadata["solver_statistics"],
        )


def _validate(chemistry: ChemistryParameters, mechanics: MechanicalParameters,
              initial: InitialCondition, simulation: SimulationParameters) -> None:
    positive = {
        "rest_length": mechanics.rest_length,
        "elastic_modulus": mechanics.elastic_modulus,
        "cell_viscosity": mechanics.cell_viscosity,
        "endpoint_friction": mechanics.endpoint_friction,
        "maximum_force": mechanics.maximum_force,
        "force_half_level": mechanics.force_half_level,
        "force_hill_exponent": mechanics.force_hill_exponent,
        "total_amount": initial.total_amount,
        "n_cells": simulation.n_cells,
        "end_time": simulation.end_time,
        "output_interval": simulation.output_interval,
    }
    invalid = [name for name, value in positive.items() if value <= 0]
    if invalid:
        raise ValueError("parameters must be positive: " + ", ".join(invalid))
    if mechanics.coupling not in {"protrusive", "contractile"}:
        raise ValueError("coupling must be 'protrusive' or 'contractile'")
    if simulation.n_cells < 5:
        raise ValueError("n_cells must be at least 5")
    if not 0 < initial.active_edge_fraction < 0.5:
        raise ValueError("active_edge_fraction must lie between 0 and 0.5")
    if initial.active_low < 0 or initial.active_high <= initial.active_low:
        raise ValueError("initial active concentrations are invalid")
    chemistry_values = asdict(chemistry)
    if any(value <= 0 for value in chemistry_values.values()):
        raise ValueError("all chemistry parameters must be positive")


def initial_state(mechanics: MechanicalParameters, initial: InitialCondition,
                  n_cells: int) -> tuple[np.ndarray, np.ndarray]:
    """Construct a conservative step-profile initial state and material grid."""
    dxbar = 2.0 / n_cells
    xbar = -1.0 + (np.arange(n_cells) + 0.5) * dxbar
    if initial.active_edge == "right":
        patch = xbar >= 1.0 - 2.0 * initial.active_edge_fraction
    else:
        patch = xbar <= -1.0 + 2.0 * initial.active_edge_fraction
    active = np.where(patch, initial.active_high, initial.active_low)
    mean_total = initial.total_amount / mechanics.rest_length
    inactive_level = mean_total - float(np.mean(active))
    if inactive_level <= 0:
        raise ValueError(
            "initial active profile exceeds the requested total amount; "
            "reduce active_high/edge_fraction or increase total_amount"
        )
    inactive = np.full(n_cells, inactive_level)
    state = np.concatenate((active, inactive,
                            [initial.centroid, mechanics.rest_length / 2.0]))
    return xbar, state


def endpoint_forces(active: np.ndarray,
                    mechanics: MechanicalParameters) -> tuple[float, float]:
    """Return signed forces at the left and right endpoints."""
    exponent = mechanics.force_hill_exponent
    edge = np.maximum(active[[0, -1]], 0.0) ** exponent
    response = edge / (mechanics.force_half_level**exponent + edge)
    direction = 1.0 if mechanics.coupling == "protrusive" else -1.0
    left = -direction * mechanics.maximum_force * response[0]
    right = direction * mechanics.maximum_force * response[1]
    return float(left), float(right)


def _laplacian_no_flux(values: np.ndarray, dx: float) -> np.ndarray:
    result = np.empty_like(values)
    result[0] = (values[1] - values[0]) / dx**2
    result[-1] = (values[-2] - values[-1]) / dx**2
    result[1:-1] = (values[:-2] - 2.0*values[1:-1] + values[2:]) / dx**2
    return result


def _jacobian_sparsity(n: int):
    """Sparsity pattern used by SciPy's finite-difference BDF Jacobian."""
    size = 2*n + 2
    pattern = lil_matrix((size, size), dtype=np.int8)
    for field_start in (0, n):
        for index in range(n):
            row = field_start + index
            pattern[row, field_start + max(0, index-1):field_start + min(n, index+2)] = 1
            pattern[row, index] = 1
            pattern[row, n + index] = 1
            pattern[row, 2*n + 1] = 1
    # Mechanics depend on the two active edge values and radius.
    pattern[2*n:, 0] = 1
    pattern[2*n:, n-1] = 1
    pattern[2*n:, 2*n + 1] = 1
    return pattern.tocsr()


def simulate(
    chemistry: ChemistryParameters | None = None,
    mechanics: MechanicalParameters | None = None,
    initial: InitialCondition | None = None,
    simulation: SimulationParameters | None = None,
) -> SimulationResult:
    """Integrate the spatial mechanochemical model on a moving interval."""
    chemistry = chemistry or ChemistryParameters()
    mechanics = mechanics or MechanicalParameters()
    initial = initial or InitialCondition()
    simulation = simulation or SimulationParameters()
    _validate(chemistry, mechanics, initial, simulation)
    n = simulation.n_cells
    xbar, state0 = initial_state(mechanics, initial, n)
    dxbar = 2.0 / n

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        active = state[:n]
        inactive = state[n:2*n]
        radius = state[-1]
        if radius <= simulation.minimum_radius:
            raise RuntimeError("cell radius reached the configured lower bound")
        force_left, force_right = endpoint_forces(active, mechanics)
        dcentroid = (force_left + force_right) / (2.0*mechanics.endpoint_friction)
        dradius = (
            force_right - force_left
            - 2.0*mechanics.elastic_modulus*(2.0*radius-mechanics.rest_length)
        ) / (2.0*(mechanics.endpoint_friction + 2.0*mechanics.cell_viscosity))
        dilution = dradius / radius
        positive_active = np.maximum(active, 0.0)
        feedback = positive_active**chemistry.hill_exponent
        activation = (chemistry.basal_activation
                      + chemistry.feedback_activation*feedback/(1.0+feedback))
        reaction = activation*inactive - chemistry.inactivation*active
        dactive = (chemistry.active_diffusion/radius**2
                   * _laplacian_no_flux(active, dxbar)
                   + reaction - dilution*active)
        dinactive = (chemistry.inactive_diffusion/radius**2
                     * _laplacian_no_flux(inactive, dxbar)
                     - reaction - dilution*inactive)
        return np.concatenate((dactive, dinactive, [dcentroid, dradius]))

    def jacobian(_time: float, state: np.ndarray):
        """Exact sparse Jacobian for the stiff method-of-lines system."""
        active = state[:n]
        inactive = state[n:2*n]
        radius = state[-1]
        force_left, force_right = endpoint_forces(active, mechanics)
        mechanical_drag = mechanics.endpoint_friction + 2.0*mechanics.cell_viscosity
        dradius = (
            force_right-force_left
            - 2.0*mechanics.elastic_modulus*(2.0*radius-mechanics.rest_length)
        ) / (2.0*mechanical_drag)
        dilution = dradius/radius

        exponent = mechanics.force_hill_exponent
        gc_power = mechanics.force_half_level**exponent
        edge_active = np.maximum(active[[0, -1]], 0.0)
        force_slope = (exponent*gc_power*edge_active**(exponent-1.0)
                       / (gc_power+edge_active**exponent)**2)
        direction = 1.0 if mechanics.coupling == "protrusive" else -1.0
        dforce_left = -direction*mechanics.maximum_force*force_slope[0]
        dforce_right = direction*mechanics.maximum_force*force_slope[1]
        dr_da_left = -dforce_left/(2.0*mechanical_drag)
        dr_da_right = dforce_right/(2.0*mechanical_drag)
        dr_dr = -2.0*mechanics.elastic_modulus/mechanical_drag
        ddilution_da = (dr_da_left/radius, dr_da_right/radius)
        ddilution_dr = dr_dr/radius - dradius/radius**2

        positive_active = np.maximum(active, 0.0)
        feedback = positive_active**chemistry.hill_exponent
        activation = (chemistry.basal_activation
                      + chemistry.feedback_activation*feedback/(1.0+feedback))
        activation_slope = np.zeros_like(active)
        positive = active > 0
        activation_slope[positive] = (
            chemistry.feedback_activation*chemistry.hill_exponent
            * positive_active[positive]**(chemistry.hill_exponent-1.0)
            / (1.0+feedback[positive])**2
        )
        reaction_slope = activation_slope*inactive-chemistry.inactivation
        lap_active = _laplacian_no_flux(active, dxbar)
        lap_inactive = _laplacian_no_flux(inactive, dxbar)
        matrix = lil_matrix((2*n+2, 2*n+2), dtype=float)
        active_diffusion = chemistry.active_diffusion/radius**2/dxbar**2
        inactive_diffusion = chemistry.inactive_diffusion/radius**2/dxbar**2
        for index in range(n):
            neighbors = ((index, -1.0), (1, 1.0)) if index == 0 else (
                ((index, -1.0), (n-2, 1.0)) if index == n-1 else
                ((index, -2.0), (index-1, 1.0), (index+1, 1.0))
            )
            for column, coefficient in neighbors:
                matrix[index, column] += active_diffusion*coefficient
                matrix[n+index, n+column] += inactive_diffusion*coefficient
            matrix[index, index] += reaction_slope[index]-dilution
            matrix[index, n+index] += activation[index]
            matrix[n+index, index] += -reaction_slope[index]
            matrix[n+index, n+index] += -activation[index]-dilution
            for column, slope in ((0, ddilution_da[0]),
                                  (n-1, ddilution_da[1])):
                matrix[index, column] += -active[index]*slope
                matrix[n+index, column] += -inactive[index]*slope
            matrix[index, 2*n+1] = (
                -2.0*chemistry.active_diffusion/radius**3*lap_active[index]
                - active[index]*ddilution_dr
            )
            matrix[n+index, 2*n+1] = (
                -2.0*chemistry.inactive_diffusion/radius**3*lap_inactive[index]
                - inactive[index]*ddilution_dr
            )
        matrix[2*n, 0] = dforce_left/(2.0*mechanics.endpoint_friction)
        matrix[2*n, n-1] = dforce_right/(2.0*mechanics.endpoint_friction)
        matrix[2*n+1, 0] = dr_da_left
        matrix[2*n+1, n-1] = dr_da_right
        matrix[2*n+1, 2*n+1] = dr_dr
        return matrix.tocsr()

    times = np.arange(0.0, simulation.end_time + 0.5*simulation.output_interval,
                      simulation.output_interval)
    if times[-1] < simulation.end_time:
        times = np.append(times, simulation.end_time)
    result = solve_ivp(
        rhs, (0.0, simulation.end_time), state0,
        method=simulation.method, t_eval=times,
        rtol=simulation.rtol, atol=simulation.atol,
        jac=jacobian,
    )
    if not result.success:
        raise RuntimeError(result.message)
    states = result.y.T
    return SimulationResult(
        times=result.t,
        xbar=xbar,
        active=states[:, :n],
        inactive=states[:, n:2*n],
        centroid=states[:, 2*n],
        radius=states[:, 2*n+1],
        chemistry=chemistry,
        mechanics=mechanics,
        initial_condition=initial,
        simulation=simulation,
        solver_statistics={
            "status": str(result.message),
            "rhs_evaluations": int(result.nfev),
            "jacobian_evaluations": int(result.njev),
            "lu_decompositions": int(result.nlu),
        },
    )


def figure9_regimes() -> dict[str, tuple[MechanicalParameters, InitialCondition]]:
    """Explicit reconstruction presets for the four regimes in paper Fig. 9.

    Figure 9 gives E and total GTPase for each panel. The remaining defaults
    come from Tables 1--2 and are intentionally exposed in the dataclasses.
    """
    base_mechanics = MechanicalParameters()
    stiff_initial = InitialCondition(total_amount=25.0)
    soft_initial = replace(stiff_initial, total_amount=10.0)
    return {
        "contractile_stiff": (
            replace(base_mechanics, coupling="contractile", elastic_modulus=1500.0),
            stiff_initial,
        ),
        "contractile_soft": (
            replace(base_mechanics, coupling="contractile", elastic_modulus=1250.0),
            soft_initial,
        ),
        "protrusive_stiff": (
            replace(base_mechanics, coupling="protrusive", elastic_modulus=1250.0),
            stiff_initial,
        ),
        "protrusive_soft": (
            replace(base_mechanics, coupling="protrusive", elastic_modulus=250.0),
            soft_initial,
        ),
    }


def scale_wgd_condition(
    mechanics: MechanicalParameters,
    initial: InitialCondition,
    volume_ratio: float = 2.0,
) -> tuple[MechanicalParameters, InitialCondition]:
    """Apply the matched isotropic WGD scaling requested for this model.

    Local biochemical concentrations and kinetics are unchanged. Because both
    chemical fields occupy the model's effective 1D support, preserving their
    concentrations scales the conserved native amount with line length, ``s``.
    The initial physical centroid is held fixed while the rest length is
    enlarged symmetrically about it.
    """
    if volume_ratio <= 0:
        raise ValueError("volume_ratio must be positive")
    scale = volume_ratio**(1.0/3.0)
    scaled_mechanics = replace(
        mechanics,
        rest_length=scale*mechanics.rest_length,
        elastic_modulus=scale*mechanics.elastic_modulus,
        cell_viscosity=scale*mechanics.cell_viscosity,
        endpoint_friction=scale**2*mechanics.endpoint_friction,
        maximum_force=scale**2*mechanics.maximum_force,
    )
    scaled_initial = replace(
        initial,
        total_amount=scale*initial.total_amount,
        centroid=initial.centroid,
    )
    return scaled_mechanics, scaled_initial
