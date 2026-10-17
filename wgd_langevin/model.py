"""Overdamped, two-dimensional active Brownian cells. Units: micrometres, hours.

Cell-specific speed is constant along a track; lognormal_cv is between-cell
heterogeneity, not temporally fluctuating propulsion. No forces or interactions.
"""
import numpy as np


def msd(t, speed, tau, cv=0.0, diffusion=0.0):
    t = np.asarray(t, dtype=float)
    q = t / tau
    return 4 * diffusion * t + 2 * speed**2 * (1 + cv**2) * tau**2 * (q + np.expm1(-q))


def simulate(speed, tau, cadence, duration, n=1024, substeps=32, seed=20260921,
             cv=0.0, diffusion=0.0, localization_sd=0.0):
    """Midpoint heading integration; exact Gaussian heading increments.

    Returns observed coordinates including independent localization noise.
    dt must resolve tau; pipeline uses dt/tau <= 1/8.
    """
    if min(speed, diffusion, cv, localization_sd) < 0 or min(tau, cadence, duration) <= 0:
        raise ValueError('Invalid physical parameters')
    if n < 1 or substeps < 1:
        raise ValueError('n and substeps must be positive')
    frames = int(round(duration / cadence))
    if not np.isclose(frames * cadence, duration):
        raise ValueError('Duration must be a multiple of cadence')
    rng = np.random.default_rng(seed)
    theta = rng.uniform(-np.pi, np.pi, n)
    sigma = np.sqrt(np.log1p(cv**2))
    speeds = speed * np.exp(sigma * rng.standard_normal(n) - sigma**2 / 2)
    dt = cadence / substeps
    positions = np.zeros((n, frames + 1, 2))
    x = np.zeros((n, 2))
    for frame in range(frames):
        for _ in range(substeps):
            # Independent half increments give a Brownian heading at the midpoint.
            theta += np.sqrt(dt / tau) * rng.standard_normal(n)
            x[:, 0] += speeds * np.cos(theta) * dt
            x[:, 1] += speeds * np.sin(theta) * dt
            theta += np.sqrt(dt / tau) * rng.standard_normal(n)
            if diffusion:
                x += np.sqrt(2 * diffusion * dt) * rng.standard_normal((n, 2))
        positions[:, frame + 1] = x
    if localization_sd:
        positions += localization_sd * rng.standard_normal(positions.shape)
    return positions


def metrics(paths, cadence):
    increments = np.diff(paths, axis=1)
    step = np.linalg.norm(increments, axis=2)
    length = step.sum(axis=1)
    end = np.linalg.norm(paths[:, -1] - paths[:, 0], axis=1)
    directness = np.divide(end, length, out=np.zeros_like(end), where=length > 0)
    return {'speed': length / (cadence * step.shape[1]), 'path_length': length,
            'endpoint': end, 'directness': directness,
            'step_msd': (step**2).mean(axis=1)}
