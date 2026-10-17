"""Fixed circular-domain Taniguchi PIP2/PIP3 model (SI Eq. S6).

The mesh is an equilateral triangular lattice.  Missing neighbours at the
circular mask are reflected to the centre value, implementing zero normal
flux on the staircase boundary.  Time integration is explicit midpoint RK2.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
from numba import njit


@dataclass(frozen=True)
class FixedParameters:
    D_u: float = 0.05
    D_v: float = 0.20
    alpha: float = 8.0
    beta: float = 2.0
    K_k: float = 5.7
    K_p: float = 5.0
    S: float = 1.0
    gamma: float = 0.15
    mu: float = 1.0
    theta: float = 1.5e-2
    noise_d: float = 0.2
    noise_sigma: float = 0.5


@dataclass(frozen=True)
class TriangularGrid:
    x: np.ndarray
    y: np.ndarray
    mask: np.ndarray
    dx: float
    radius: float
    node_area: float

    @classmethod
    def circle(cls, radius: float = 5.0, dx: float = 5.0 / 31.0) -> "TriangularGrid":
        dy = sqrt(3.0) * dx / 2.0
        ny = int(np.ceil(2.0 * radius / dy)) + 3
        nx = int(np.ceil(2.0 * radius / dx)) + 3
        rows = np.arange(ny)
        cols = np.arange(nx)
        y = (rows - (ny - 1) / 2.0) * dy
        x = (cols[None, :] - (nx - 1) / 2.0) * dx
        x = x + 0.5 * dx * ((rows[:, None] & 1) - 0.5)
        x = np.broadcast_to(x, (ny, nx)).copy()
        yy = np.broadcast_to(y[:, None], (ny, nx)).copy()
        mask = x * x + yy * yy <= radius * radius
        return cls(x=x, y=yy, mask=mask, dx=dx, radius=radius,
                   node_area=sqrt(3.0) * dx * dx / 2.0)

    @property
    def area(self) -> float:
        return float(self.mask.sum() * self.node_area)


@dataclass
class FixedResult:
    grid: TriangularGrid
    parameters: FixedParameters
    times: np.ndarray
    U: np.ndarray
    V: np.ndarray
    event_count: int


@njit(cache=True)
def _rhs(U, V, mask, dx, D_u, D_v, alpha, beta, K_k, K_p, S, gamma, mu):
    ny, nx = U.shape
    sum_u = 0.0
    sum_v2 = 0.0
    count = 0
    for j in range(ny):
        for i in range(nx):
            if mask[j, i]:
                sum_u += U[j, i]
                sum_v2 += V[j, i] * V[j, i]
                count += 1
    mean_u = sum_u / count
    mean_v2 = sum_v2 / count
    du = np.zeros_like(U)
    dv = np.zeros_like(V)
    factor = 2.0 / (3.0 * dx * dx)
    for j in range(ny):
        odd = j & 1
        for i in range(nx):
            if not mask[j, i]:
                continue
            lap_u = 0.0
            lap_v = 0.0
            # Horizontal plus four diagonal neighbours on an offset-row grid.
            ni = (i - 1, i + 1, i - 1 + odd, i + odd, i - 1 + odd, i + odd)
            nj = (j, j, j - 1, j - 1, j + 1, j + 1)
            for k in range(6):
                ii, jj = ni[k], nj[k]
                if 0 <= ii < nx and 0 <= jj < ny and mask[jj, ii]:
                    lap_u += U[jj, ii] - U[j, i]
                    lap_v += V[jj, ii] - V[j, i]
            lap_u *= factor
            lap_v *= factor
            uv = U[j, i] * V[j, i]
            forward = alpha * uv * V[j, i] / (K_k + mean_v2)
            reverse = beta * uv / (K_p + mean_u)
            du[j, i] = D_u * lap_u - forward + reverse + S - gamma * U[j, i]
            dv[j, i] = D_v * lap_v + forward - reverse - mu * V[j, i]
    return du, dv


@njit(cache=True)
def _rk2_step(U, V, mask, dx, dt, values):
    du1, dv1 = _rhs(U, V, mask, dx, *values)
    Um = U + 0.5 * dt * du1
    Vm = V + 0.5 * dt * dv1
    du2, dv2 = _rhs(Um, Vm, mask, dx, *values)
    U += dt * du2
    V += dt * dv2


def gaussian_transfer(U: np.ndarray, V: np.ndarray, grid: TriangularGrid,
                      x_c: float, y_c: float, amplitude: float,
                      width: float) -> None:
    """Apply an instantaneous, exactly balanced U-to-V Gaussian transfer."""
    amount = amplitude * np.exp(-((grid.x - x_c) ** 2 + (grid.y - y_c) ** 2)
                                / (2.0 * width * width))
    amount *= grid.mask
    U -= amount
    V += amount


def _random_inside(rng: np.random.Generator, radius: float) -> tuple[float, float]:
    r = radius * np.sqrt(rng.random())
    angle = 2.0 * np.pi * rng.random()
    return r * np.cos(angle), r * np.sin(angle)


def simulate_fixed(
    end_time: float,
    *,
    parameters: FixedParameters | None = None,
    grid: TriangularGrid | None = None,
    dt: float = 5e-4,
    output_interval: float = 0.1,
    seed: int = 1,
    stochastic: bool = True,
    initial_transfers: tuple[tuple[float, float, float, float], ...] = (),
    initial_U: np.ndarray | None = None,
    initial_V: np.ndarray | None = None,
) -> FixedResult:
    """Simulate SI Eq. S6; transfer tuples are ``(x, y, amplitude, width)``."""
    p = parameters or FixedParameters()
    g = grid or TriangularGrid.circle()
    U = np.zeros(g.mask.shape, dtype=np.float64)
    V = np.zeros_like(U)
    U[g.mask] = 1.0 / p.gamma
    if (initial_U is None) != (initial_V is None):
        raise ValueError("initial_U and initial_V must be supplied together")
    if initial_U is not None:
        if initial_U.shape != g.mask.shape or initial_V.shape != g.mask.shape:
            raise ValueError("initial fields must match the grid shape")
        U[g.mask] = np.asarray(initial_U, dtype=float)[g.mask]
        V[g.mask] = np.asarray(initial_V, dtype=float)[g.mask]
    for transfer in initial_transfers:
        gaussian_transfer(U, V, g, *transfer)

    steps = int(np.ceil(end_time / dt))
    every = max(1, int(round(output_interval / dt)))
    save_steps = list(range(0, steps + 1, every))
    if save_steps[-1] != steps:
        save_steps.append(steps)
    times = np.empty(len(save_steps))
    Us = np.empty((len(save_steps), *U.shape), dtype=np.float32)
    Vs = np.empty_like(Us)
    times[0], Us[0], Vs[0] = 0.0, U, V
    next_save = 1
    rng = np.random.default_rng(seed)
    event_count = 0
    values = (p.D_u, p.D_v, p.alpha, p.beta, p.K_k, p.K_p,
              p.S, p.gamma, p.mu)
    for step in range(1, steps + 1):
        _rk2_step(U, V, g.mask, g.dx, dt, values)
        if stochastic:
            n_events = rng.poisson(p.theta * g.area * dt)
            for _ in range(n_events):
                xc, yc = _random_inside(rng, g.radius)
                gaussian_transfer(U, V, g, xc, yc,
                                  rng.exponential(p.noise_sigma), p.noise_d)
            event_count += n_events
        if next_save < len(save_steps) and step == save_steps[next_save]:
            times[next_save] = min(step * dt, end_time)
            Us[next_save], Vs[next_save] = U, V
            next_save += 1
    return FixedResult(g, p, times, Us, Vs, event_count)
