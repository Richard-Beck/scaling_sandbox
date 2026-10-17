"""Numerical building blocks for the one-dimensional two-layer mechanics.

The functions here intentionally stop before Rac/Rho kinetics.  They solve the
linear mechanical limits used in the explanatory report and expose every
derived field needed to check the force balances.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


def _validate_times(times):
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or times[0] < 0:
        raise ValueError("times must be an increasing one-dimensional array")
    if np.any(np.diff(times) <= 0):
        raise ValueError("times must be strictly increasing")
    return times


def membrane_pull(*, mu, times, length=10.0, stiffness=1.0, pull=1.0, n=241):
    """Step-displace the left end of an elastic membrane over a fixed cortex.

    The PDE is ``x_t = (k / mu) x_ss`` with ``x(0)=pull`` and ``x(L)=0``.
    Returned tension is ``k x_s``.  The initial boundary jump represents a
    sudden experimental pull; output begins at the requested first time.
    """
    times = _validate_times(times)
    if mu <= 0 or length <= 0 or stiffness <= 0 or n < 5:
        raise ValueError("mu, length and stiffness must be positive; n >= 5")
    s = np.linspace(0.0, length, n)
    ds = s[1] - s[0]
    diffusivity = stiffness / mu
    interior_n = n - 2

    def rhs(_time, interior):
        full = np.empty(n)
        full[0], full[-1] = pull, 0.0
        full[1:-1] = interior
        return diffusivity * (full[:-2] - 2.0 * full[1:-1] + full[2:]) / ds**2

    sol = solve_ivp(
        rhs,
        (0.0, times[-1]),
        np.zeros(interior_n),
        t_eval=times,
        method="BDF",
        rtol=2e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    displacement = np.zeros((len(times), n))
    displacement[:, 0] = pull
    displacement[:, 1:-1] = sol.y.T
    tension = stiffness * np.gradient(displacement, ds, axis=1, edge_order=2)
    return {
        "s": s,
        "time": times,
        "displacement": displacement,
        "tension": tension,
        "diffusivity": diffusivity,
    }


def active_stress_response(
    *,
    mu,
    times,
    length=10.0,
    stiffness=1.0,
    viscosity=100.0,
    stress_amplitude=1.0,
    stress_center=1.5,
    stress_width=0.55,
    pulse_stop=np.inf,
    n=256,
):
    """Solve the periodic coupled mechanics driven by prescribed active stress.

    For every nonzero Fourier mode, force balance gives

    ``v_q = -(k/eta)x_q + i sigma_q/(eta q)``

    and membrane balance gives

    ``(x_q)_t = i sigma_q/(eta q) - k(1/eta + q^2/mu)x_q``.

    The zero mode is set to zero because uniform translation carries no strain.
    A Gaussian active stress is mean-subtracted and switched off at pulse_stop.
    """
    times = _validate_times(times)
    if min(mu, length, stiffness, viscosity, stress_width) <= 0 or n < 8:
        raise ValueError("physical parameters must be positive; n >= 8")
    s = np.linspace(0.0, length, n, endpoint=False)
    ds = length / n
    distance = (s - stress_center + length / 2.0) % length - length / 2.0
    sigma_on = stress_amplitude * np.exp(-0.5 * (distance / stress_width) ** 2)
    sigma_on -= sigma_on.mean()
    sigma_hat = np.fft.fft(sigma_on)
    q = 2.0 * np.pi * np.fft.fftfreq(n, d=ds)
    nonzero = q != 0
    forcing = np.zeros(n, dtype=complex)
    forcing[nonzero] = 1j * sigma_hat[nonzero] / (viscosity * q[nonzero])
    decay = stiffness * (1.0 / viscosity + q**2 / mu)
    decay[0] = 0.0

    x_hat = np.zeros((len(times), n), dtype=complex)
    for j, time in enumerate(times):
        on_time = min(time, pulse_stop)
        on_factor = np.zeros(n)
        on_factor[nonzero] = -np.expm1(-decay[nonzero] * on_time) / decay[nonzero]
        state_at_stop = forcing * on_factor
        if time > pulse_stop:
            state_at_stop *= np.exp(-decay * (time - pulse_stop))
        x_hat[j] = state_at_stop

    displacement = np.fft.ifft(x_hat, axis=1).real
    tension = np.fft.ifft(1j * stiffness * q[None, :] * x_hat, axis=1).real
    sigma = np.where((times <= pulse_stop)[:, None], sigma_on[None, :], 0.0)
    sigma_hat_t = np.where(
        (times <= pulse_stop)[:, None], sigma_hat[None, :], 0.0
    )
    velocity_hat = -(stiffness / viscosity) * x_hat
    velocity_hat[:, nonzero] += (
        1j * sigma_hat_t[:, nonzero] / (viscosity * q[nonzero])
    )
    velocity = np.fft.ifft(velocity_hat, axis=1).real
    xdot_hat = velocity_hat - (stiffness / mu) * q[None, :] ** 2 * x_hat
    membrane_velocity = np.fft.ifft(xdot_hat, axis=1).real
    slip = velocity - membrane_velocity
    friction = mu * slip
    return {
        "s": s,
        "time": times,
        "active_stress": sigma,
        "displacement": displacement,
        "tension": tension,
        "cortex_velocity": velocity,
        "membrane_velocity": membrane_velocity,
        "slip": slip,
        "friction": friction,
        "lambda": np.sqrt(viscosity / mu),
        "tau": viscosity / stiffness,
    }


def variable_friction_response(
    *,
    mu_profile,
    times,
    length=10.0,
    stiffness=1.0,
    pull=0.0,
    cortex_velocity=None,
):
    """Solve ``x_t = v + (k/mu(s)) x_ss`` on a fixed interval.

    Displacement is fixed to ``pull`` at the left boundary and zero at the
    right.  This is the direct local-friction form of membrane force balance,
    not a Fickian diffusion equation with ``(D x_s)_s``.
    """
    times = _validate_times(times)
    mu_profile = np.asarray(mu_profile, dtype=float)
    if mu_profile.ndim != 1 or len(mu_profile) < 5 or np.any(mu_profile <= 0):
        raise ValueError("mu_profile must be a positive one-dimensional array")
    n = len(mu_profile)
    s = np.linspace(0.0, length, n)
    ds = s[1] - s[0]
    if cortex_velocity is None:
        velocity = np.zeros(n)
    elif callable(cortex_velocity):
        velocity = np.asarray(cortex_velocity(s), dtype=float)
    else:
        velocity = np.asarray(cortex_velocity, dtype=float)
    if velocity.shape != (n,):
        raise ValueError("cortex_velocity must have the shape of mu_profile")

    def rhs(_time, interior):
        full = np.empty(n)
        full[0], full[-1] = pull, 0.0
        full[1:-1] = interior
        lap = (full[:-2] - 2.0 * full[1:-1] + full[2:]) / ds**2
        return velocity[1:-1] + stiffness * lap / mu_profile[1:-1]

    sol = solve_ivp(
        rhs,
        (0.0, times[-1]),
        np.zeros(n - 2),
        t_eval=times,
        method="BDF",
        rtol=2e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    displacement = np.zeros((len(times), n))
    displacement[:, 0] = pull
    displacement[:, 1:-1] = sol.y.T
    tension = stiffness * np.gradient(displacement, ds, axis=1, edge_order=2)
    return {"s": s, "time": times, "displacement": displacement, "tension": tension}
