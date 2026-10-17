"""Conservative diffuse-domain Eroume/Maree migration model.

The nine biochemical reactions and parameters are a direct transcription of
``R/vcell_model.cpp`` and ``scenarios/wgd_spatial_scale.vcml``.  The boundary
coupling is deliberately separate and phenomenological; it is not attributed
to Eroume et al.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Callable

import numpy as np
from numba import njit

from .stimulus import StimulusProtocol

SPECIES = ("cdc42_inact", "cdc42_act", "Rac_inact", "Rac_act",
           "Rho_inact", "Rho_act", "P1", "P2", "P3")


@dataclass(frozen=True)
class ChemistryParameters:
    # Order and values match ReactionParameters in R/vcell_model.cpp.
    n: float = 3.0
    a1: float = 1.25
    a2: float = 1.0
    alpha: float = 4.5
    beta: float = 0.3
    cdc42_total: float = 2.4
    rac_total: float = 7.5
    rho_total: float = 3.1
    d_cdc42: float = 1.0
    d_rac: float = 1.0
    d_rho: float = 1.0
    lipid_feedback_fraction: float = 0.4
    I_p: float = 3.3
    I_p1: float = 10.5
    I_r: float = 0.5
    k21: float = 0.014
    k_pi3k: float = 7.2e-4
    k_pi5k: float = 0.084
    k_pten: float = 0.432
    p3_baseline: float = 0.05
    rac_baseline: float = 3.0
    rho_baseline: float = 1.25
    delta_p1: float = 0.21
    diffusion: tuple[float, ...] = (50.0, 0.1, 50.0, 0.1, 50.0, 0.1,
                                    5.0, 5.0, 5.0)
    initial: tuple[float, ...] = (1.4, 1.0, 4.5, 3.0, 1.85, 1.25,
                                  55.0, 30.0, 0.05)

    def values(self) -> tuple[float, ...]:
        return (self.n, self.a1, self.a2, self.alpha, self.beta,
                self.cdc42_total, self.rac_total, self.rho_total,
                self.d_cdc42, self.d_rac, self.d_rho,
                self.lipid_feedback_fraction, self.I_p, self.I_p1,
                self.I_r, self.k21, self.k_pi3k, self.k_pi5k,
                self.k_pten, self.p3_baseline, self.rac_baseline,
                self.rho_baseline, self.delta_p1)


@dataclass(frozen=True)
class MotilityParameters:
    epsilon: float = 1.5
    interface_mobility: float = 1.0
    surface_tension: float = 1.0
    area_penalty: float = 0.02
    protrusion: float = 0.10
    cdc42_weight: float = 0.35
    retraction: float = 0.10
    concentration_floor: float = 1.0e-3


@dataclass(frozen=True)
class SimulationConfig:
    duration: float = 600.0
    dt: float = 0.004
    dx: float = 1.0
    box_size: float = 96.0
    cell_radius: float = 12.5
    center_x: float | None = None
    center_y: float | None = None
    output_interval: float = 2.0
    seed: int = 1
    noise_amplitude: float = 0.01
    freeze_boundary: bool = False
    chemistry: ChemistryParameters = field(default_factory=ChemistryParameters)
    motility: MotilityParameters = field(default_factory=MotilityParameters)
    stimulus: StimulusProtocol | None = None

    def __post_init__(self) -> None:
        if self.duration < 0 or self.dt <= 0 or self.dx <= 0 or self.box_size <= 0:
            raise ValueError("duration must be non-negative and dt/dx/box_size positive")
        if self.output_interval <= 0 or self.cell_radius <= 0:
            raise ValueError("output_interval and cell_radius must be positive")
        if max(self.chemistry.diffusion) * self.dt / self.dx**2 > 0.24:
            raise ValueError("explicit diffusion stability requires max(D)*dt/dx^2 <= 0.24")


@dataclass
class SimulationResult:
    config: SimulationConfig
    x: np.ndarray
    y: np.ndarray
    times: np.ndarray
    phi: np.ndarray
    concentrations: np.ndarray
    material: np.ndarray
    exterior_material: np.ndarray


def reaction_rates(concentrations: np.ndarray, stimulus: np.ndarray | float,
                   parameters: ChemistryParameters | None = None) -> np.ndarray:
    """Vectorized public reference version of the exact nine reaction rates."""
    p = parameters or ChemistryParameters()
    c = concentrations
    feedback = (1.0 - p.lipid_feedback_fraction
                + p.lipid_feedback_fraction * c[8] / p.p3_baseline)
    j0 = (stimulus / (1.0 + (c[5] / p.a1) ** p.n)
          * (c[0] / p.cdc42_total) * feedback - p.d_cdc42 * c[1])
    j1 = ((p.I_r + p.alpha * c[1]) * (c[2] / p.rac_total) * feedback
          - p.d_rac * c[3])
    j2 = ((p.I_p + p.beta * c[3]) / (1.0 + (c[1] / p.a2) ** p.n)
          * (c[4] / p.rho_total) - p.d_rho * c[5])
    j3 = (-p.k21 * c[7]
          + p.k_pi5k / 2.0 * (1.0 + c[3] / p.rac_baseline) * c[6])
    j4 = (p.k_pi3k / 2.0 * (1.0 + c[3] / p.rac_baseline) * c[7]
          - p.k_pten / 2.0 * (1.0 + c[5] / p.rho_baseline) * c[8])
    # This source is intentionally independent of the represented PI pool.
    j5 = p.I_p1 - p.delta_p1 * c[6]
    return np.stack((-j0, j0, -j1, j1, -j2, j2, -j3 + j5, j3 - j4, j4))


@njit(cache=True, fastmath=True)
def _step(phi, q, stimulus, dx, dt, chem, diffusion, motility, target_area,
          freeze_boundary, phi_next, q_next, concentrations):
    ny, nx = phi.shape
    (n, a1, a2, alpha, beta, ctot, ractot, rhotot, dcdc, drac, drho,
     feedback_fraction, ip, ip1, ir, k21, kpi3k, kpi5k, kpten, p3b,
     rb, rhob, delta_p1) = chem
    (epsilon, mobility, surface_tension, area_penalty, protrusion,
     cdc_weight, retraction, floor) = motility
    invdx2 = 1.0 / (dx * dx)
    area_sum = 0.0
    interface_sum = 0.0
    signal_sum = 0.0
    for j in range(ny):
        for i in range(nx):
            ph = phi[j, i]
            area_sum += ph
            denom = ph if ph > floor else floor
            for k in range(9):
                concentrations[k, j, i] = q[k, j, i] / denom
    # The interface-weighted mean removes uniform biochemical pressure; the
    # area constraint then acts only on residual numerical/curvature changes.
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            gx = 0.5 * (phi[j, i + 1] - phi[j, i - 1]) / dx
            gy = 0.5 * (phi[j + 1, i] - phi[j - 1, i]) / dx
            grad = math.sqrt(gx * gx + gy * gy)
            signal = (protrusion * (concentrations[3, j, i]
                                     + cdc_weight * concentrations[1, j, i])
                      - retraction * concentrations[5, j, i])
            interface_sum += grad
            signal_sum += grad * signal
    mean_signal = signal_sum / max(interface_sum, floor)
    area_error = area_sum * dx * dx - target_area
    for k in range(9):
        for j in range(ny):
            q_next[k, j, 0] = q[k, j, 0]
            q_next[k, j, nx - 1] = q[k, j, nx - 1]
        for i in range(nx):
            q_next[k, 0, i] = q[k, 0, i]
            q_next[k, ny - 1, i] = q[k, ny - 1, i]
    for j in range(ny):
        phi_next[j, 0] = phi[j, 0]
        phi_next[j, nx - 1] = phi[j, nx - 1]
    for i in range(nx):
        phi_next[0, i] = phi[0, i]
        phi_next[ny - 1, i] = phi[ny - 1, i]
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            ph = phi[j, i]
            c0 = concentrations[0, j, i]; c1 = concentrations[1, j, i]
            c2 = concentrations[2, j, i]; c3 = concentrations[3, j, i]
            c4 = concentrations[4, j, i]; c5 = concentrations[5, j, i]
            c6 = concentrations[6, j, i]; c7 = concentrations[7, j, i]
            c8 = concentrations[8, j, i]
            feedback = (1.0 - feedback_fraction) + feedback_fraction * c8 / p3b
            j0 = stimulus[j, i] / (1.0 + (c5 / a1) ** n) * (c0 / ctot) * feedback - dcdc * c1
            j1 = (ir + alpha * c1) * (c2 / ractot) * feedback - drac * c3
            j2 = (ip + beta * c3) / (1.0 + (c1 / a2) ** n) * (c4 / rhotot) - drho * c5
            j3 = -k21 * c7 + kpi5k / 2.0 * (1.0 + c3 / rb) * c6
            j4 = kpi3k / 2.0 * (1.0 + c3 / rb) * c7 - kpten / 2.0 * (1.0 + c5 / rhob) * c8
            j5 = ip1 - delta_p1 * c6
            rates0 = -j0; rates1 = j0; rates2 = -j1; rates3 = j1
            rates4 = -j2; rates5 = j2; rates6 = -j3 + j5
            rates7 = j3 - j4; rates8 = j4
            for k in range(9):
                center = concentrations[k, j, i]
                # Closed laboratory-box faces. Every internal face appears
                # with opposite signs in its two cells, so diffusion is
                # conservative to roundoff.
                right = 0.0 if i == nx - 2 else 0.5 * (ph + phi[j, i + 1]) * (concentrations[k, j, i + 1] - center)
                left = 0.0 if i == 1 else 0.5 * (ph + phi[j, i - 1]) * (center - concentrations[k, j, i - 1])
                top = 0.0 if j == ny - 2 else 0.5 * (ph + phi[j + 1, i]) * (concentrations[k, j + 1, i] - center)
                bottom = 0.0 if j == 1 else 0.5 * (ph + phi[j - 1, i]) * (center - concentrations[k, j - 1, i])
                flux = (right - left + top - bottom) * invdx2
                rate = rates0
                if k == 1: rate = rates1
                elif k == 2: rate = rates2
                elif k == 3: rate = rates3
                elif k == 4: rate = rates4
                elif k == 5: rate = rates5
                elif k == 6: rate = rates6
                elif k == 7: rate = rates7
                elif k == 8: rate = rates8
                value = q[k, j, i] + dt * (diffusion[k] * flux + ph * rate)
                q_next[k, j, i] = max(0.0, value)
            if freeze_boundary:
                phi_next[j, i] = ph
            else:
                gx = 0.5 * (phi[j, i + 1] - phi[j, i - 1]) / dx
                gy = 0.5 * (phi[j + 1, i] - phi[j - 1, i]) / dx
                grad = math.sqrt(gx * gx + gy * gy)
                lap = (phi[j, i + 1] + phi[j, i - 1] + phi[j + 1, i]
                       + phi[j - 1, i] - 4.0 * ph) * invdx2
                # Standard double-well derivative for the initialized
                # tanh(z/(sqrt(2)*epsilon)) interface profile.
                gp = ph * (1.0 - ph) * (1.0 - 2.0 * ph)
                signal = (protrusion * (c3 + cdc_weight * c1) - retraction * c5)
                rhs = (surface_tension * (lap - gp / (epsilon * epsilon))
                       - area_penalty * area_error * grad
                       + (signal - mean_signal) * grad) / mobility
                value = ph + dt * rhs
                phi_next[j, i] = min(1.0, max(0.0, value))


def _smooth_noise(rng: np.random.Generator, shape: tuple[int, int]) -> np.ndarray:
    noise = rng.normal(size=shape)
    for _ in range(5):
        noise = (4 * noise + np.roll(noise, 1, 0) + np.roll(noise, -1, 0)
                 + np.roll(noise, 1, 1) + np.roll(noise, -1, 1)) / 8.0
    noise -= noise.mean()
    scale = noise.std()
    return noise / scale if scale else noise


def run_simulation(config: SimulationConfig = SimulationConfig(), *,
                   progress: Callable[[int, int], None] | None = None) -> SimulationResult:
    """Run a trajectory and retain only requested float32 movie snapshots."""
    center_x = config.box_size / 2 if config.center_x is None else config.center_x
    center_y = config.box_size / 2 if config.center_y is None else config.center_y
    x = np.arange(0.0, config.box_size + 0.5 * config.dx, config.dx)
    y = x.copy()
    X, Y = np.meshgrid(x, y)
    epsilon = config.motility.epsilon
    radius = np.hypot(X - center_x, Y - center_y)
    phi = 0.5 * (1.0 - np.tanh((radius - config.cell_radius) /
                                (math.sqrt(2.0) * epsilon)))
    c = np.empty((9, *phi.shape), dtype=np.float64)
    for k, value in enumerate(config.chemistry.initial):
        c[k] = value
    if config.noise_amplitude:
        perturbation = config.noise_amplitude * _smooth_noise(
            np.random.default_rng(config.seed), phi.shape)
        # Pairwise transfers preserve each GTPase total at every pixel.
        for inactive, active in ((0, 1), (2, 3), (4, 5)):
            delta = perturbation * config.chemistry.initial[active]
            c[active] = np.maximum(0.0, c[active] + delta)
            c[inactive] = np.maximum(0.0, c[inactive] - delta)
    q = phi[None, :, :] * c
    phi_next = np.empty_like(phi)
    q_next = np.empty_like(q)
    concentration_workspace = np.empty_like(q)
    target_area = float(phi.sum() * config.dx**2)
    steps = int(math.ceil(config.duration / config.dt))
    every = max(1, int(round(config.output_interval / config.dt)))
    save_steps = list(range(0, steps + 1, every))
    if save_steps[-1] != steps:
        save_steps.append(steps)
    phis = np.empty((len(save_steps), *phi.shape), dtype=np.float32)
    concentrations = np.empty((len(save_steps), 9, *phi.shape), dtype=np.float32)
    material = np.empty((len(save_steps), 9), dtype=np.float64)
    exterior_material = np.empty((len(save_steps), 9), dtype=np.float64)
    times = np.empty(len(save_steps))

    def save(index: int, step: int) -> None:
        phis[index] = phi
        concentrations[index] = q / np.maximum(phi, config.motility.concentration_floor)
        material[index] = q.sum(axis=(1, 2)) * config.dx**2
        # Deep exterior excludes the diffuse interface itself; material at
        # phi<0.05 is the relevant leakage diagnostic.
        exterior_material[index] = q[:, phi < 0.05].sum(axis=1) * config.dx**2
        times[index] = min(step * config.dt, config.duration)

    save(0, 0)
    next_save = 1
    protocol = config.stimulus or StimulusProtocol(
        mode="free", pivot_x=center_x, pivot_y=center_y,
        cell_center_x=center_x, cell_center_y=center_y)
    chem_values = np.asarray(config.chemistry.values(), dtype=np.float64)
    diffusion = np.asarray(config.chemistry.diffusion, dtype=np.float64)
    motility_values = np.asarray((
        config.motility.epsilon, config.motility.interface_mobility,
        config.motility.surface_tension, config.motility.area_penalty,
        config.motility.protrusion, config.motility.cdc42_weight,
        config.motility.retraction, config.motility.concentration_floor,
    ), dtype=np.float64)
    stimulus = np.asarray(protocol.evaluate(X, Y, 0.0), dtype=np.float64)
    for step in range(1, steps + 1):
        time = (step - 1) * config.dt
        # All built-in protocols are piecewise constant in time. Rebuild the
        # laboratory-coordinate field only at a protocol transition.
        crossed_rotate = (protocol.mode == "rotate" and time >= protocol.rotate_at
                          and time - config.dt < protocol.rotate_at)
        crossed_transient = (protocol.mode == "transient"
                             and time >= protocol.transient_duration
                             and time - config.dt < protocol.transient_duration)
        crossed_legacy = (protocol.mode == "legacy" and time > 10.0
                          and time - config.dt <= 10.0)
        if crossed_rotate or crossed_transient or crossed_legacy:
            stimulus = np.asarray(protocol.evaluate(X, Y, time), dtype=np.float64)
        _step(phi, q, stimulus, config.dx, config.dt, chem_values, diffusion,
              motility_values, target_area, config.freeze_boundary,
              phi_next, q_next, concentration_workspace)
        phi, phi_next = phi_next, phi
        q, q_next = q_next, q
        if next_save < len(save_steps) and step == save_steps[next_save]:
            save(next_save, step)
            if progress is not None:
                progress(step, steps)
            next_save += 1
    return SimulationResult(config, x, y, times, phis, concentrations, material,
                            exterior_material)
