"""Deforming-cell Taniguchi model (SI Eqs. S8--S9 and Table S1)."""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numba import njit


@dataclass(frozen=True)
class PhaseFieldParameters:
    D_u: float = 0.05
    D_v: float = 0.20
    alpha: float = 8.0
    beta: float = 4.8
    K_k: float = 3.5
    K_p: float = 5.2
    S: float = 1.0
    gamma: float = 0.15
    mu: float = 1.0
    chi_u: float = 50.0
    chi_v: float = 50.0
    tau: float = 0.83
    epsilon: float = 1.0
    eta: float = 1.0
    area0: float = 25.0 * np.pi
    area_penalty: float = 0.5
    protrusion_a: float = 10.0
    retraction_b: float = 0.030
    theta: float = 1.0e-3
    noise_d: float = 0.07
    noise_sigma: float = 5.7


PHASE_PARAMETER_SETS = {
    "4F": PhaseFieldParameters(),
    "4G": replace(PhaseFieldParameters(), alpha=12.0, beta=3.0, K_k=5.8,
                  K_p=3.0, gamma=0.28, chi_u=80.0, chi_v=80.0,
                  protrusion_a=14.0, noise_d=0.06, noise_sigma=7.3),
    "4H": replace(PhaseFieldParameters(), alpha=12.0, beta=4.0, K_k=6.0,
                  K_p=3.0, protrusion_a=12.0),
    "4I": replace(PhaseFieldParameters(), alpha=12.0, beta=5.0, K_k=6.0,
                  K_p=2.8, gamma=0.20, protrusion_a=12.0,
                  retraction_b=0.035, theta=0.85e-3, noise_d=0.10,
                  noise_sigma=4.0),
}


@dataclass
class PhaseFieldResult:
    x: np.ndarray
    y: np.ndarray
    parameters: PhaseFieldParameters
    times: np.ndarray
    phi: np.ndarray
    U: np.ndarray
    V: np.ndarray
    event_count: int


@njit(cache=True)
def _phase_rhs(phi, q_u, q_v, dx, values):
    (D_u, D_v, alpha, beta, K_k, K_p, S, gamma, mu, chi_u, chi_v,
     tau, epsilon, eta, area0, area_penalty, protrusion_a,
     retraction_b) = values
    ny, nx = phi.shape
    # Concentrations are undefined in the vanishing phase-field tail.  A
    # small diffuse-domain cutoff prevents division by an exponentially tiny
    # phi from feeding exterior roundoff back through the flux stencil.
    floor = 1e-3
    U = q_u / np.maximum(phi, floor)
    V = q_v / np.maximum(phi, floor)
    area = np.sum(phi) * dx * dx
    mean_u = np.sum(phi * U) * dx * dx / area
    mean_v2 = np.sum(phi * V * V) * dx * dx / area
    grad2 = np.zeros_like(phi)
    lap_phi = np.zeros_like(phi)
    div_u = np.zeros_like(phi)
    div_v = np.zeros_like(phi)
    inv2dx = 0.5 / dx
    invdx2 = 1.0 / (dx * dx)
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            gx = (phi[j, i + 1] - phi[j, i - 1]) * inv2dx
            gy = (phi[j + 1, i] - phi[j - 1, i]) * inv2dx
            grad2[j, i] = gx * gx + gy * gy
            lap_phi[j, i] = (phi[j, i + 1] + phi[j, i - 1]
                             + phi[j + 1, i] + phi[j - 1, i]
                             - 4.0 * phi[j, i]) * invdx2
            # Face-averaged phi gives central, conservative div(phi grad c).
            div_u[j, i] = (
                0.5 * (phi[j, i + 1] + phi[j, i]) * (U[j, i + 1] - U[j, i])
                - 0.5 * (phi[j, i] + phi[j, i - 1]) * (U[j, i] - U[j, i - 1])
                + 0.5 * (phi[j + 1, i] + phi[j, i]) * (U[j + 1, i] - U[j, i])
                - 0.5 * (phi[j, i] + phi[j - 1, i]) * (U[j, i] - U[j - 1, i])
            ) * invdx2
            div_v[j, i] = (
                0.5 * (phi[j, i + 1] + phi[j, i]) * (V[j, i + 1] - V[j, i])
                - 0.5 * (phi[j, i] + phi[j, i - 1]) * (V[j, i] - V[j, i - 1])
                + 0.5 * (phi[j + 1, i] + phi[j, i]) * (V[j + 1, i] - V[j, i])
                - 0.5 * (phi[j, i] + phi[j - 1, i]) * (V[j, i] - V[j - 1, i])
            ) * invdx2
    surface_norm = np.sum(grad2) * dx * dx
    phi_t = np.zeros_like(phi)
    du = np.zeros_like(phi)
    dv = np.zeros_like(phi)
    area_error = area - area0
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            ph = phi[j, i]
            u, v = U[j, i], V[j, i]
            grad = np.sqrt(grad2[j, i])
            # G'(phi), G=18 phi^2 (1-phi)^2.
            gp = 36.0 * ph * (1.0 - ph) * (1.0 - 2.0 * ph)
            phi_t[j, i] = (eta * (lap_phi[j, i] - gp / (epsilon * epsilon))
                           - area_penalty * area_error * grad
                           + (protrusion_a * v - retraction_b * u) * grad) / tau
            uv = u * v
            forward = alpha * uv * v / (K_k + mean_v2)
            reverse = beta * uv / (K_p + mean_u)
            leak_u = chi_u * u * grad2[j, i] / max(surface_norm, floor)
            leak_v = chi_v * v * grad2[j, i] / max(surface_norm, floor)
            du[j, i] = (D_u * div_u[j, i]
                        + ph * (-forward + reverse + S - gamma * u) - leak_u)
            dv[j, i] = (D_v * div_v[j, i]
                        + ph * (forward - reverse - mu * v) - leak_v)
    return phi_t, du, dv


@njit(cache=True)
def _phase_step(phi, q_u, q_v, dx, dt, values):
    p1, u1, v1 = _phase_rhs(phi, q_u, q_v, dx, values)
    pm = phi + 0.5 * dt * p1
    um = q_u + 0.5 * dt * u1
    vm = q_v + 0.5 * dt * v1
    p2, u2, v2 = _phase_rhs(pm, um, vm, dx, values)
    phi += dt * p2
    q_u += dt * u2
    q_v += dt * v2
    # The continuum phase field is bounded. Tiny explicit overshoots otherwise
    # amplify q/phi in the nominally exterior region.
    for j in range(phi.shape[0]):
        for i in range(phi.shape[1]):
            if phi[j, i] < 0.0:
                phi[j, i] = 0.0
            elif phi[j, i] > 1.0:
                phi[j, i] = 1.0
            if phi[j, i] < 1e-3:
                q_u[j, i] = 0.0
                q_v[j, i] = 0.0


def _initial_phi(X: np.ndarray, Y: np.ndarray, p: PhaseFieldParameters) -> np.ndarray:
    radius = np.sqrt(p.area0 / np.pi)
    return 0.5 * (1.0 - np.tanh((np.hypot(X, Y) - radius)
                                / (np.sqrt(2.0) * p.epsilon)))


def simulate_phase_field(
    end_time: float,
    *,
    parameter_set: str = "4F",
    parameters: PhaseFieldParameters | None = None,
    dx: float = 0.1,
    dt: float = 8e-5,
    box_half_width: float = 10.0,
    output_interval: float = 0.05,
    seed: int = 1,
    stochastic: bool = True,
    initial_transfers: tuple[tuple[float, float, float, float], ...] = (),
) -> PhaseFieldResult:
    """Simulate the deformable model on the specified fixed computational box."""
    if parameters is None:
        try:
            p = PHASE_PARAMETER_SETS[parameter_set]
        except KeyError as exc:
            raise ValueError(f"unknown parameter set {parameter_set!r}") from exc
    else:
        p = parameters
    x = np.arange(-box_half_width, box_half_width + 0.5 * dx, dx)
    y = x.copy()
    X, Y = np.meshgrid(x, y)
    phi = _initial_phi(X, Y, p)
    U = np.full_like(phi, 1.0 / p.gamma)
    V = np.zeros_like(phi)
    for xc, yc, amplitude, width in initial_transfers:
        transfer = amplitude * np.exp(-((X - xc) ** 2 + (Y - yc) ** 2)
                                      / (2.0 * width * width))
        U -= transfer
        V += transfer
    q_u, q_v = phi * U, phi * V
    values = (p.D_u, p.D_v, p.alpha, p.beta, p.K_k, p.K_p, p.S, p.gamma,
              p.mu, p.chi_u, p.chi_v, p.tau, p.epsilon, p.eta, p.area0,
              p.area_penalty, p.protrusion_a, p.retraction_b)
    steps = int(np.ceil(end_time / dt))
    every = max(1, int(round(output_interval / dt)))
    save_steps = list(range(0, steps + 1, every))
    if save_steps[-1] != steps:
        save_steps.append(steps)
    shape = (len(save_steps), *phi.shape)
    phis = np.empty(shape, dtype=np.float32)
    Us = np.empty(shape, dtype=np.float32)
    Vs = np.empty(shape, dtype=np.float32)
    times = np.empty(len(save_steps))

    def save(index: int, time: float) -> None:
        phis[index] = phi
        Us[index] = q_u / np.maximum(phi, 1e-3)
        Vs[index] = q_v / np.maximum(phi, 1e-3)
        times[index] = time

    save(0, 0.0)
    next_save = 1
    rng = np.random.default_rng(seed)
    events = 0
    for step in range(1, steps + 1):
        _phase_step(phi, q_u, q_v, dx, dt, values)
        if stochastic:
            area = float(phi.sum() * dx * dx)
            n_events = rng.poisson(p.theta * area * dt)
            for _ in range(n_events):
                # "Intracellular" is reconstructed as the sharp interior of
                # the diffuse domain.  Sampling with phi weights would place
                # a material fraction of events in the lossy interface.
                inside = np.flatnonzero(phi.ravel() >= 0.5)
                index = inside[rng.integers(inside.size)]
                j, i = np.unravel_index(index, phi.shape)
                amplitude = rng.exponential(p.noise_sigma)
                transfer = amplitude * np.exp(-((X - X[j, i]) ** 2
                                                + (Y - Y[j, i]) ** 2)
                                               / (2.0 * p.noise_d ** 2))
                q_u -= phi * transfer
                q_v += phi * transfer
            events += n_events
        if next_save < len(save_steps) and step == save_steps[next_save]:
            save(next_save, min(step * dt, end_time))
            next_save += 1
    return PhaseFieldResult(x, y, p, times, phis, Us, Vs, events)
