"""Standalone 1D polarity runners extracted from the public report models."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp


def laplacian_no_flux(values: np.ndarray, dx: float) -> np.ndarray:
    out = np.empty_like(values)
    out[0] = (values[1] - values[0]) / dx**2
    out[-1] = (values[-2] - values[-1]) / dx**2
    out[1:-1] = (values[:-2] - 2 * values[1:-1] + values[2:]) / dx**2
    return out


def _integrate(rhs: Callable, initial: np.ndarray, end_time: float,
               output_interval: float) -> tuple[np.ndarray, np.ndarray]:
    times = np.arange(0.0, end_time + 0.5 * output_interval, output_interval)
    if times[-1] < end_time:
        times = np.append(times, end_time)
    result = solve_ivp(rhs, (0.0, end_time), initial, method="BDF",
                       t_eval=times, rtol=2e-6, atol=1e-8)
    if not result.success:
        raise RuntimeError(result.message)
    return result.t, result.y.T


def simulate_mori(length: float, concentration_scale: float = 1.0,
                  n: int = 81, end_time: float = 300.0,
                  output_interval: float = 1.0) -> dict[str, np.ndarray]:
    dx = length / n
    x = (np.arange(n) + 0.5) * dx
    a0, b0 = 0.2683 * concentration_scale, 2.0 * concentration_scale
    initial = np.concatenate((np.full(n, a0), np.full(n, b0)))
    stimulus_width = 0.04 * length

    def rhs(time: float, state: np.ndarray) -> np.ndarray:
        a, b = state[:n], state[n:]
        stimulus = np.where(
            (time >= 50.0) & (time < 70.0) & (x <= stimulus_width), 10.0, 0.0
        )
        reaction = b * (0.067 + a * a / (1.0 + a * a) + stimulus) - a
        return np.concatenate((0.1 * laplacian_no_flux(a, dx) + reaction,
                               10.0 * laplacian_no_flux(b, dx) - reaction))

    times, states = _integrate(rhs, initial, end_time, output_interval)
    return {"times": times, "x": x, "active": states[:, :n],
            "inactive": states[:, n:], "support_length": np.asarray(length)}


HOLMES_BASE = {
    "Ct": 2.4, "Rt": 7.5, "rhot": 3.1, "IC": 2.95, "IR": 0.2,
    "Irho": 6.6, "a1": 1.25, "a2": 1.0, "hill": 3.0,
    "feedback": 0.55, "Dm": 0.1, "Dmc0": 50.0,
}


def holmes_accessibility(volume_ratio: float, baseline_fraction: float) -> float:
    current = 1.0 / (1.0 + ((1.0 - baseline_fraction) / baseline_fraction)
                     * volume_ratio ** (2.0 / 3.0))
    return current / baseline_fraction


def holmes_inactive_diffusion(volume_ratio: float,
                               baseline_fraction: float) -> float:
    p = HOLMES_BASE
    cytosolic = (p["Dmc0"] - baseline_fraction * p["Dm"]) / (1-baseline_fraction)
    current = baseline_fraction * holmes_accessibility(volume_ratio, baseline_fraction)
    return current * p["Dm"] + (1-current) * cytosolic


def _holmes_reaction(C, Rac, Rho, Ci, Raci, Rhoi, accessibility, signal,
                     totals):
    p = HOLMES_BASE
    i_c = p["IC"] / (1 + (Rho / p["a1"]) ** p["hill"])
    i_rac = p["IR"] + p["feedback"] * C + signal
    i_rho = p["Irho"] / (1 + (Rac / p["a2"]) ** p["hill"])
    return (
        i_c * accessibility * Ci / totals[0] - C,
        i_rac * accessibility * Raci / totals[1] - Rac,
        i_rho * accessibility * Rhoi / totals[2] - Rho,
    )


def simulate_holmes(length: float, volume_ratio: float,
                    concentration_scale: float = 1.0,
                    baseline_fraction: float = 0.5, n: int = 41,
                    end_time: float = 300.0,
                    output_interval: float = 5.0,
                    accessibility_override: float | None = None,
                    inactive_diffusion_override: float | None = None,
                    signal_mean: float = 0.3,
                    signal_half_range: float = 0.1) -> dict[str, np.ndarray]:
    p = HOLMES_BASE
    totals = concentration_scale * np.asarray([p["Ct"], p["Rt"], p["rhot"]])
    accessibility = holmes_accessibility(volume_ratio, baseline_fraction)
    d_inactive = holmes_inactive_diffusion(volume_ratio, baseline_fraction)
    if accessibility_override is not None:
        accessibility = float(accessibility_override)
    if inactive_diffusion_override is not None:
        d_inactive = float(inactive_diffusion_override)

    homogeneous0 = np.asarray([0.1, 0.1, 0.1,
                               totals[0]-0.1, totals[1]-0.1, totals[2]-0.1])

    def homogeneous_rhs(time, state):
        reactions = _holmes_reaction(*state, accessibility, 0.0, totals)
        return np.asarray([*reactions, *(-np.asarray(reactions))])

    relaxed = solve_ivp(homogeneous_rhs, (0, 300), homogeneous0, method="BDF",
                        rtol=1e-9, atol=1e-11).y[:, -1]
    dx = length / n
    x = (np.arange(n) + 0.5) * dx
    signal = signal_mean + signal_half_range * (2*x/length - 1)
    initial = np.concatenate([np.full(n, value) for value in relaxed])

    def rhs(time, state):
        fields = [state[i*n:(i+1)*n] for i in range(6)]
        reactions = _holmes_reaction(*fields, accessibility, signal, totals)
        active = [p["Dm"] * laplacian_no_flux(fields[i], dx) + reactions[i]
                  for i in range(3)]
        inactive = [d_inactive * laplacian_no_flux(fields[i+3], dx)-reactions[i]
                    for i in range(3)]
        return np.concatenate(active + inactive)

    times, states = _integrate(rhs, initial, end_time, output_interval)
    fields = states.reshape(states.shape[0], 6, n)
    return {"times": times, "x": x,
            # Legacy generic keys are Rac active/inactive and remain stable.
            "active": fields[:, 1, :], "inactive": fields[:, 4, :],
            "cdc42": fields[:, 0, :], "rho": fields[:, 2, :],
            # The inactive Cdc42 and Rho pools complete the six-field output.
            "cdc42_inactive": fields[:, 3, :],
            "rho_inactive": fields[:, 5, :],
            "signal": signal,
            "support_length": np.asarray(length),
            "accessibility": np.asarray(accessibility),
            "inactive_diffusion": np.asarray(d_inactive)}


def simulate_holmes_explicit_exchange(
        length: float, volume_ratio: float, exchange_scale: float,
        concentration_scale: float = 1.0, baseline_fraction: float = 0.5,
        n: int = 41, end_time: float = 300.0,
        output_interval: float = 5.0, signal_mean: float = 0.3,
        signal_half_range: float = 0.1) -> dict[str, np.ndarray]:
    """Simulate the unreduced active/membrane-inactive/cytosolic network.

    ``exchange_scale`` multiplies both association and dissociation rates while
    their ratio remains fixed. The projected exchange flux is constructed so
    its rapid-equilibrium membrane fraction is exactly the fraction used by
    :func:`holmes_accessibility`. This makes the singular limit an explicit
    numerical check of the six-field quasi-steady-state reduction.
    """

    if exchange_scale <= 0:
        raise ValueError("exchange_scale must be positive")
    p = HOLMES_BASE
    totals = concentration_scale * np.asarray([p["Ct"], p["Rt"], p["rhot"]])
    geometry_load = volume_ratio ** (2.0 / 3.0)
    association = exchange_scale * baseline_fraction / (1.0-baseline_fraction)
    dissociation = exchange_scale
    qss_fraction = association / (association + geometry_load*dissociation)
    cytosolic_diffusion = (
        p["Dmc0"] - baseline_fraction*p["Dm"]
    ) / (1.0-baseline_fraction)

    def reactions(active, membrane_inactive, signal):
        # GEFs act on membrane-inactive GTPase. Dividing by the baseline
        # fraction preserves the activation-rate normalization in the reduced
        # equations, where accessibility is expressed relative to c0.
        accessible_equivalent = membrane_inactive / baseline_fraction
        return np.asarray(_holmes_reaction(
            active[0], active[1], active[2],
            accessible_equivalent[0], accessible_equivalent[1],
            accessible_equivalent[2], 1.0, signal, totals,
        ))

    # Obtain the common unstimulated equilibrium from the reduced homogeneous
    # system, then split its inactive pool at the exact exchange equilibrium.
    # This avoids giving slow-exchange cases an incompletely relaxed initial
    # condition merely because lambda is small.
    reduced0 = np.asarray([
        0.1, 0.1, 0.1, totals[0]-0.1, totals[1]-0.1, totals[2]-0.1,
    ])

    def reduced_homogeneous_rhs(time, state):
        reaction = np.asarray(_holmes_reaction(
            *state, qss_fraction/baseline_fraction, 0.0, totals,
        ))
        return np.concatenate((reaction, -reaction))

    reduced_relaxed = solve_ivp(
        reduced_homogeneous_rhs, (0.0, 300.0), reduced0, method="BDF",
        rtol=1e-9, atol=1e-11,
    ).y[:, -1]
    relaxed = np.concatenate((
        reduced_relaxed[:3],
        qss_fraction*reduced_relaxed[3:],
        (1.0-qss_fraction)*reduced_relaxed[3:],
    ))
    dx = length/n
    x = (np.arange(n)+0.5)*dx
    signal = signal_mean + signal_half_range*(2*x/length-1)
    initial = np.concatenate([np.full(n, value) for value in relaxed])

    def rhs(time, state):
        fields = [state[i*n:(i+1)*n] for i in range(9)]
        active = np.asarray(fields[:3])
        membrane_inactive = np.asarray(fields[3:6])
        cytosolic_inactive = np.asarray(fields[6:9])
        reaction = reactions(active, membrane_inactive, signal)
        exchange = (
            association*cytosolic_inactive
            - geometry_load*dissociation*membrane_inactive
        )
        d_active = np.asarray([
            p["Dm"]*laplacian_no_flux(active[i], dx)+reaction[i]
            for i in range(3)
        ])
        d_membrane = np.asarray([
            p["Dm"]*laplacian_no_flux(membrane_inactive[i], dx)
            - reaction[i]+exchange[i]
            for i in range(3)
        ])
        d_cytosolic = np.asarray([
            cytosolic_diffusion*laplacian_no_flux(cytosolic_inactive[i], dx)
            - exchange[i]
            for i in range(3)
        ])
        return np.concatenate((
            d_active.reshape(-1),
            d_membrane.reshape(-1),
            d_cytosolic.reshape(-1),
        ))

    times, states = _integrate(rhs, initial, end_time, output_interval)
    fields = states.reshape(states.shape[0], 9, n)
    return {
        "times": times,
        "x": x,
        "active": fields[:, :3, :],
        "membrane_inactive": fields[:, 3:6, :],
        "cytosolic_inactive": fields[:, 6:9, :],
        "signal": signal,
        "support_length": np.asarray(length),
        "exchange_scale": np.asarray(exchange_scale),
        "association_rate": np.asarray(association),
        "dissociation_rate": np.asarray(dissociation),
        "geometry_load": np.asarray(geometry_load),
        "qss_membrane_fraction": np.asarray(qss_fraction),
        "cytosolic_diffusion": np.asarray(cytosolic_diffusion),
    }


def simulate_buttenschoen(length: float, concentration_scale: float = 1.0,
                          n: int = 81, end_time: float = 1200.0,
                          output_interval: float = 5.0) -> dict[str, np.ndarray]:
    dx_bar = 2.0 / n
    xbar = -1 + (np.arange(n) + 0.5) * dx_bar
    radius = length / 2.0
    mean_total = 2.5 * concentration_scale
    active_profile = np.where(xbar <= -0.75, 1.5, 0.1)
    active0 = concentration_scale * active_profile
    inactive0 = np.full(
        n, concentration_scale * (2.5 - np.mean(active_profile))
    )
    if np.min(inactive0) <= 0:
        raise ValueError("Buttenschoen fixed-total initial inactive field is nonpositive")
    initial = np.concatenate((active0, inactive0))

    def rhs(time, state):
        active, inactive = state[:n], state[n:]
        activation = 0.1 + 2.0 * active**4 / (1 + active**4)
        reaction = activation * inactive - active
        return np.concatenate((0.01/radius**2 * laplacian_no_flux(active, dx_bar)
                               + reaction,
                               1.0/radius**2 * laplacian_no_flux(inactive, dx_bar)
                               - reaction))

    times, states = _integrate(rhs, initial, end_time, output_interval)
    return {"times": times, "x": xbar * radius, "active": states[:, :n],
            "inactive": states[:, n:], "support_length": np.asarray(length)}
