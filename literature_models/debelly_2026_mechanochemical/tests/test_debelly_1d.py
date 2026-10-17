from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from debelly_1d.model import active_stress_response, membrane_pull  # noqa: E402
from debelly_1d.discrete import (  # noqa: E402
    DiscreteParameters,
    DiscreteProtocol,
    active_nodal_forces,
    simulate_discrete,
    smooth_step,
    zero_protocol,
)


def test_membrane_pull_converges_to_uniform_tension():
    result = membrane_pull(mu=1.0, times=[1.0, 400.0], length=10.0, n=121)
    final = result["tension"][-1]
    assert np.max(np.abs(final + 0.1)) < 2e-4


def test_fourier_mode_decay_rate_is_one_plus_lambda_squared_q_squared():
    length = 10.0
    stiffness = 2.0
    viscosity = 80.0
    mu = 5.0
    tau = viscosity / stiffness
    lam = np.sqrt(viscosity / mu)
    q = 2.0 * np.pi / length
    expected = (1.0 + lam**2 * q**2) / tau
    direct = stiffness * (1.0 / viscosity + q**2 / mu)
    assert np.isclose(direct, expected)


def test_internal_friction_matches_membrane_tension_gradient():
    result = active_stress_response(
        mu=4.0, times=[1.0, 20.0], length=10.0, n=512
    )
    tension = result["tension"][-1]
    ds = result["s"][1] - result["s"][0]
    tension_gradient = np.fft.ifft(
        1j
        * 2.0
        * np.pi
        * np.fft.fftfreq(len(tension), d=ds)
        * np.fft.fft(tension)
    ).real
    # Membrane balance: T_s = mu (x_t - v) = -mu(v-x_t).
    assert np.max(np.abs(tension_gradient + result["friction"][-1])) < 2e-9

    q = 2.0 * np.pi * np.fft.fftfreq(len(tension), d=ds)
    cortex_viscous_force = np.fft.ifft(
        -(q**2) * 100.0 * np.fft.fft(result["cortex_velocity"][-1])
    ).real
    active_stress_gradient = np.fft.ifft(
        1j * q * np.fft.fft(result["active_stress"][-1])
    ).real
    cortex_residual = cortex_viscous_force + active_stress_gradient - result["friction"][-1]
    assert np.max(np.abs(cortex_residual)) < 2e-9


def test_discrete_solver_satisfies_both_layer_force_balances():
    params = DiscreteParameters(nodes=17, linker_friction=2.0)

    def active_stress(time, links):
        return (1 - np.exp(-time)) * np.exp(-0.5 * ((links - 8.0) / 1.2) ** 2)

    trajectory = simulate_discrete(
        params,
        DiscreteProtocol(active_stress=active_stress),
        np.linspace(0, 4, 9),
    )
    membrane_residual = trajectory.spring_force + trajectory.linker_force
    assert np.max(np.abs(membrane_residual)) < 2e-10

    v = trajectory.cortex_velocity
    viscous = params.cortical_viscosity * (v[:, 2:] - 2 * v[:, 1:-1] + v[:, :-2])
    active = np.array(
        [active_nodal_forces(stress, params.nodes)[1:-1] for stress in trajectory.active_stress]
    )
    cortex_residual = viscous + active - trajectory.linker_force[:, 1:-1]
    assert np.max(np.abs(cortex_residual)) < 2e-10


def test_discrete_friction_reverses_membrane_pull_and_cortex_loading_response():
    times = np.linspace(0, 10, 41)

    def pulled_end(time):
        return smooth_step(time, amplitude=2.0, ramp_time=1.0)

    def driven_cortex(time):
        value, _ = smooth_step(time, amplitude=0.35, ramp_time=1.0)
        return value, 0.0

    responses = {}
    for friction in [0.05, 5.0]:
        params = DiscreteParameters(nodes=21, linker_friction=friction)
        pulled = simulate_discrete(
            params,
            DiscreteProtocol(membrane_left=zero_protocol, membrane_right=pulled_end),
            times,
        )
        driven = simulate_discrete(
            params,
            DiscreteProtocol(cortex_right=driven_cortex),
            times,
        )
        responses[friction] = {
            "pull_transmission": abs(pulled.spring_tension[-1, 0] / pulled.spring_tension[-1, -1]),
            "drive_tension": np.linalg.norm(driven.spring_tension[-1]),
        }
    assert responses[0.05]["pull_transmission"] > responses[5.0]["pull_transmission"]
    assert responses[0.05]["drive_tension"] < responses[5.0]["drive_tension"]
