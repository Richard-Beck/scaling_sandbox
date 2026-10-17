"""Holmes et al. (2012) minimal Rac--Rho mutual-inhibition model (Model 3).

This adapter implements the explicit no-phosphoinositide two-GTPase version,
using the conserved active/composite-inactive reduction in their Eq. (1) and
the activation terms in Eqs. (5)--(6). See ``source`` in :meth:`metadata` for
the primary source and the choices needed to use its defaults in our common
1D assay protocol.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np

from .contract import Grid1D


@dataclass(frozen=True)
class HolmesModel3:
    """Two conserved active/inactive Rac and Rho pools, without PI feedback.

    Concentrations are μM, time is seconds, and length is μm. ``stimulus`` is
    an added Rac GEF activation rate (μM/s), following the external signal
    term in the paper's Model 4 formulation. The paper's tabulated membrane
    cycling constants are unavailable; therefore the tabulated effective
    inactive diffusivity at L0 is held fixed at all lengths, and the
    membrane-accessibility factor is held at its L0 reference value (1).
    """

    rac_total: float = 7.5
    rho_total: float = 3.1
    # Retained as the first-pass low-end reference value. Table 1 lists the
    # Model 4 coefficients I_R1 and I_R2; it does not identify a unique Model 3
    # I_hat_R, so this numeric legacy default is not a published Model 3 value.
    rac_activation: float = 0.2
    rho_activation: float = 6.6
    rho_half_inhibition_of_rac: float = 1.25
    rac_half_inhibition_of_rho: float = 1.0
    hill_n: float = 3.0
    rac_decay: float = 1.0
    rho_decay: float = 1.0
    active_diffusion: float = 0.1
    inactive_diffusion: float = 50.0
    reference_length: float = 20.0

    name = "holmes_model3"
    state_names = ("rac_active", "rac_inactive", "rho_active", "rho_inactive")
    observable_name = "rac_active"
    preferred_step = 0.02
    nonnegative = True
    conserved_groups = ((0, 1), (2, 3))

    @property
    def diffusivities(self) -> tuple[float, float, float, float]:
        return (self.active_diffusion, self.inactive_diffusion,
                self.active_diffusion, self.inactive_diffusion)

    def _homogeneous_active(self) -> tuple[float, float]:
        """Return the low-activity homogeneous fixed point by fixed-point iteration."""
        rac, rho = 0.0, 0.0
        for _ in range(10000):
            rac_alpha = self.rac_activation / (
                2.0 * (1.0 + (rho / self.rho_half_inhibition_of_rac) ** self.hill_n)
            )
            rho_alpha = self.rho_activation / (
                2.0 * (1.0 + (rac / self.rac_half_inhibition_of_rho) ** self.hill_n)
            )
            new_rac = self.rac_total * rac_alpha / (rac_alpha + self.rac_decay * self.rac_total)
            new_rho = self.rho_total * rho_alpha / (rho_alpha + self.rho_decay * self.rho_total)
            if max(abs(new_rac - rac), abs(new_rho - rho)) < 1e-13:
                return new_rac, new_rho
            rac, rho = new_rac, new_rho
        raise RuntimeError("Holmes Model 3 homogeneous fixed point did not converge")

    def initial_state(self, grid: Grid1D) -> np.ndarray:
        rac, rho = self._homogeneous_active()
        return np.vstack((
            np.full(grid.cells, rac), np.full(grid.cells, self.rac_total - rac),
            np.full(grid.cells, rho), np.full(grid.cells, self.rho_total - rho),
        ))

    def reaction(self, time, x, state, stimulus) -> np.ndarray:
        del time, x
        rac, rac_i, rho, rho_i = state
        # Original Eq. (5): Rac is suppressed by active Rho; external cue adds
        # to Rac GEF. Both Eq. (5) and Eq. (6) include the leading factor 2.
        rac_rate = self.rac_activation / (
            2.0 * (1.0 + (rho / self.rho_half_inhibition_of_rac) ** self.hill_n)
        ) + stimulus
        # Original Eq. (6): Rho is suppressed by active Rac.
        rho_rate = self.rho_activation / (
            2.0 * (1.0 + (rac / self.rac_half_inhibition_of_rho) ** self.hill_n)
        )
        rac_transfer = rac_rate * rac_i / self.rac_total - self.rac_decay * rac
        rho_transfer = rho_rate * rho_i / self.rho_total - self.rho_decay * rho
        return np.vstack((rac_transfer, -rac_transfer, rho_transfer, -rho_transfer))

    def metadata(self) -> Mapping[str, Any]:
        return {
            **asdict(self),
            "source": {
                "citation": "Holmes WR, Lin B, Levchenko A, Edelstein-Keshet L (2012), PLoS Comput Biol 8(6):e1002366",
                "doi": "10.1371/journal.pcbi.1002366",
                "url": "https://doi.org/10.1371/journal.pcbi.1002366",
                "version": "Model 3, no-PI Rac-Rho mutual-inhibition subsystem",
                "equations": "Original Eq. (1) with G=Rac,Rho and Eqs. (5)-(6), including factor 2 in both Hill denominators; external Rac cue follows the additive S term in Eq. (8).",
                "parameter_table": "Table 1 lists the Model 4 Rac activation components I_R1 and I_R2, each 0.2; it does not specify a unique basal Rac activation I_hat_R for Model 3.",
            },
            "assumptions": [
                "Use the Model 3 subsystem without PI feedback, as specified for the two-GTPase comparison.",
                "Hold effective inactive diffusivity Dmc at its tabulated L0 value (50 μm²/s) across lengths because membrane cycling constants are not tabulated.",
                "Hold the normalized membrane accessibility factor c-tilde(L) at its reference value 1; length-dependent membrane cycling cannot be evaluated from the published defaults.",
                "Initial state is the low-activity homogeneous fixed point of the no-cue equations.",
                "The default rac_activation=0.2 is retained only as a legacy first-pass reference value; it is not represented as a published Model 3 default.",
            ],
            "units": {"concentration": "μM", "time": "s", "length": "μm", "diffusion": "μm²/s"},
        }
