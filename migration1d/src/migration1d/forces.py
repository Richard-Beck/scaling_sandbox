"""Composable mappings from endpoint field samples to endpoint forces."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal, Protocol

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]
Side = Literal["left", "right"]
Compartment = Literal["membrane", "cytoplasm"]
Direction = Literal["lab", "outward", "inward"]


@dataclass(frozen=True)
class EndpointState:
    """Read-only state supplied to an endpoint force mapping.

    Custom mappings may inspect every field. :meth:`sample` provides the
    common operation of selecting one field at the endpoint being evaluated.
    Returned forces use laboratory-frame signs: positive points to the right.
    """

    time: float
    side: Side
    position: Array
    membrane: Array
    cytoplasm: Array
    species: tuple[str, ...]

    @property
    def index(self) -> int:
        return 0 if self.side == "left" else -1

    def sample(self, field: str | int, compartment: Compartment = "membrane") -> float:
        if isinstance(field, str):
            try:
                field_index = self.species.index(field)
            except ValueError as error:
                raise KeyError(f"unknown field {field!r}; choices are {self.species}") from error
        else:
            field_index = field
        values = self.membrane if compartment == "membrane" else self.cytoplasm
        try:
            return float(values[field_index, self.index])
        except IndexError as error:
            raise IndexError(f"field index {field_index} is out of range") from error


class EndpointForce(Protocol):
    """Any callable returning one laboratory-frame endpoint force."""

    def __call__(self, state: EndpointState) -> float: ...


@dataclass(frozen=True)
class ZeroForce:
    def __call__(self, state: EndpointState) -> float:
        return 0.0


def identity(value: float) -> float:
    return value


@dataclass(frozen=True)
class FieldForce:
    """Map one named/indexed field sample into force.

    ``response`` may be any scalar user function. ``direction='outward'``
    automatically uses a negative sign at the left edge and positive at the
    right. With ``direction='lab'``, gain and offset directly set signed force.
    """

    field: str | int
    compartment: Compartment = "membrane"
    response: Callable[[float], float] = identity
    gain: float = 1.0
    offset: float = 0.0
    direction: Direction = "outward"

    def __post_init__(self) -> None:
        if self.compartment not in {"membrane", "cytoplasm"}:
            raise ValueError("compartment must be 'membrane' or 'cytoplasm'")
        if self.direction not in {"lab", "outward", "inward"}:
            raise ValueError("direction must be 'lab', 'outward', or 'inward'")

    def __call__(self, state: EndpointState) -> float:
        value = self.offset + self.gain * float(
            self.response(state.sample(self.field, self.compartment))
        )
        if self.direction == "lab":
            return value
        outward_sign = -1.0 if state.side == "left" else 1.0
        return value * (outward_sign if self.direction == "outward" else -outward_sign)


@dataclass(frozen=True)
class SumForces:
    """Add several endpoint mappings without defining a new callable class."""

    terms: tuple[EndpointForce, ...]

    def __call__(self, state: EndpointState) -> float:
        return float(sum(term(state) for term in self.terms))


@dataclass(frozen=True)
class RacRhoForce:
    """Outward Rac activity minus inward Rho activity at an endpoint.

    The signed outward activity is
    ``rac_gain*h(Rac) - rho_gain*h(Rho)``.  It is converted to laboratory
    coordinates using the endpoint side, so front Rac and rear Rho cooperate
    to move a polarized cell in the same direction.
    """

    rac_field: str | int = "rac"
    rho_field: str | int = "rho"
    rac_gain: float = 0.08
    rho_gain: float = 0.08
    rac_half_level: float = 1.0
    rho_half_level: float = 1.0
    hill_exponent: float = 4.0

    def __post_init__(self) -> None:
        if self.rac_gain < 0.0 or self.rho_gain < 0.0:
            raise ValueError("Rac and Rho force gains must be non-negative")
        if self.rac_half_level <= 0.0 or self.rho_half_level <= 0.0:
            raise ValueError("Rac and Rho force half-levels must be positive")
        if self.hill_exponent <= 0.0:
            raise ValueError("force hill exponent must be positive")

    def _hill(self, value: float, half_level: float) -> float:
        positive = max(value, 0.0)
        power = positive**self.hill_exponent
        return power / (half_level**self.hill_exponent + power)

    def __call__(self, state: EndpointState) -> float:
        rac = self._hill(
            state.sample(self.rac_field, "membrane"), self.rac_half_level
        )
        rho = self._hill(
            state.sample(self.rho_field, "membrane"), self.rho_half_level
        )
        outward_activity = self.rac_gain * rac - self.rho_gain * rho
        outward_sign = -1.0 if state.side == "left" else 1.0
        return outward_sign * outward_activity
