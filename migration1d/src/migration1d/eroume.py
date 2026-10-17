"""Nine-field Eroume--Maree chemistry adapted to a 1D moving interval."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from maree_migration.model import (
    ChemistryParameters,
    SPECIES,
    reaction_rates,
)


@dataclass(frozen=True)
class EroumeGeometry:
    """Coordinate mapping from the successful local 2D WGD reference.

    ``native_axis_length`` is the physical x span of the cell mask in the
    saved 1x checkpoint. ``source_x_*`` are its baseline-equivalent VCML cue
    coordinates. A shorter 1D interval rescales diffusion by length squared,
    preserving the reaction--diffusion timescale in relative coordinates.
    """

    native_axis_length: float = 70.54013
    source_x_min: float = 0.3125
    source_x_max: float = 49.6875


@dataclass(frozen=True)
class EroumeTwoEdgeCue:
    """Competing endpoint Cdc42 inputs on top of uniform basal ``Ic``.

    Endpoint widths are integrated over the nodal control volumes, so the
    represented cue width is independent of grid alignment and cannot wrap
    around the finite interval.
    """

    left_level: float = 5.10
    right_level: float = 4.65
    basal_level: float = 2.95
    width_fraction: float = 0.064
    start: float = 0.0
    stop: float = 10.0

    def __post_init__(self) -> None:
        if min(self.left_level, self.right_level, self.basal_level) < 0.0:
            raise ValueError("Ic levels must be non-negative")
        if not 0.0 < self.width_fraction <= 0.5:
            raise ValueError("width_fraction must be in (0, 0.5]")
        if self.stop <= self.start:
            raise ValueError("cue stop must be after start")

    def values(self, time: float, position: np.ndarray) -> np.ndarray:
        values = np.full_like(position, self.basal_level, dtype=float)
        if not self.start <= time < self.stop:
            return values
        fraction = (position-position[0])/(position[-1]-position[0])
        faces = np.empty(len(fraction)+1)
        faces[0], faces[-1] = 0.0, 1.0
        faces[1:-1] = 0.5*(fraction[:-1]+fraction[1:])
        widths = faces[1:]-faces[:-1]
        left_overlap = np.maximum(
            0.0, np.minimum(faces[1:], self.width_fraction)-faces[:-1]
        )
        right_overlap = np.maximum(
            0.0,
            faces[1:]-np.maximum(faces[:-1], 1.0-self.width_fraction),
        )
        values += (self.left_level-self.basal_level)*left_overlap/widths
        values += (self.right_level-self.basal_level)*right_overlap/widths
        return values


class EroumeChemistry:
    """Exact nine reactions with the source model's transient Cdc42 cue.

    All nine nontrivial fields occupy ``membrane`` in migration1d's paired
    storage API. The companion ``cytoplasm`` rows are identically zero. This
    is a storage adapter only: native initial concentrations, reaction terms,
    and relative diffusion timescales are preserved.
    """

    species = tuple(SPECIES)

    def __init__(
        self,
        parameters: ChemistryParameters | None = None,
        *,
        reference_length: float = 10.0,
        geometry: EroumeGeometry | None = None,
        cue: EroumeTwoEdgeCue | None = None,
    ) -> None:
        if reference_length <= 0.0:
            raise ValueError("reference_length must be positive")
        self.parameters = parameters or ChemistryParameters()
        self.geometry = geometry or EroumeGeometry()
        self.cue = cue
        if self.geometry.native_axis_length <= 0.0:
            raise ValueError("native Eroume axis length must be positive")
        self.reference_length = float(reference_length)
        coordinate_scale = self.reference_length / self.geometry.native_axis_length
        self.diffusion_scale = coordinate_scale**2
        self.membrane_diffusion = (
            np.asarray(self.parameters.diffusion, dtype=float)*self.diffusion_scale
        )
        self.cytoplasm_diffusion = np.zeros(9)
        # Do not pre-relax: P1 has a retained source, so relaxation changes the
        # PI state and no longer reproduces the successful VCML protocol.
        self.initial_membrane = np.asarray(self.parameters.initial, dtype=float)
        self.initial_cytoplasm = np.zeros(9)

    def stimulus(self, time: float, position: np.ndarray) -> np.ndarray:
        if self.cue is not None:
            return self.cue.values(time, position)
        if time <= 10.0:
            relative = (position-position[0])/(position[-1]-position[0])
            source_x = (
                self.geometry.source_x_min
                + relative*(self.geometry.source_x_max-self.geometry.source_x_min)
            )
            return 2.6 + 0.05*source_x
        return np.full_like(position, 2.95)

    def reaction(self, time, position, membrane, cytoplasm, strain, stress):
        _ = cytoplasm, strain, stress
        rate = reaction_rates(
            membrane, self.stimulus(time, position), self.parameters
        )
        return rate, np.zeros_like(rate)

    def protrusive_activity(self, membrane, cytoplasm):
        _ = cytoplasm
        return np.maximum(membrane[3], 0.0)

    def tension_deactivation(self, membrane, cytoplasm, added_rate):
        _ = cytoplasm
        loss = np.asarray(added_rate, dtype=float) * np.maximum(membrane[3], 0.0)
        membrane_rate = np.zeros_like(membrane)
        membrane_rate[2] = loss
        membrane_rate[3] = -loss
        return membrane_rate, np.zeros_like(membrane_rate)
