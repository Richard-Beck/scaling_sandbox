from __future__ import annotations

import numpy as np

from taniguchi.fixed import (
    FixedParameters,
    TriangularGrid,
    gaussian_transfer,
    simulate_fixed,
)
from taniguchi.phase_field import PHASE_PARAMETER_SETS, simulate_phase_field


def test_table_s1_parameter_sets_are_complete():
    assert set(PHASE_PARAMETER_SETS) == {"4F", "4G", "4H", "4I"}
    assert PHASE_PARAMETER_SETS["4G"].chi_u == 80.0
    assert PHASE_PARAMETER_SETS["4I"].theta == 0.85e-3
    assert np.isclose(PHASE_PARAMETER_SETS["4F"].area0, 25 * np.pi)


def test_gaussian_impulse_is_exactly_balanced():
    grid = TriangularGrid.circle(radius=1.0, dx=0.2)
    U = np.full(grid.mask.shape, 4.0)
    V = np.zeros_like(U)
    before = (U + V).copy()
    gaussian_transfer(U, V, grid, 0.1, -0.1, 0.7, 0.2)
    np.testing.assert_allclose((U + V)[grid.mask], before[grid.mask], atol=1e-14)


def test_homogeneous_rest_is_stable_without_firing():
    p = FixedParameters()
    result = simulate_fixed(0.05, parameters=p,
                            grid=TriangularGrid.circle(radius=1.0, dx=0.2),
                            stochastic=False, output_interval=0.05)
    np.testing.assert_allclose(result.U[-1][result.grid.mask], 1 / p.gamma,
                               atol=1e-12)
    np.testing.assert_allclose(result.V[-1][result.grid.mask], 0.0, atol=1e-12)


def test_subthreshold_and_suprathreshold_responses_differ():
    grid = TriangularGrid.circle(radius=2.0, dx=0.2)
    low = simulate_fixed(1.0, grid=grid, stochastic=False, output_interval=0.05,
                         initial_transfers=((0, 0, 0.2, 0.3),))
    high = simulate_fixed(1.0, grid=grid, stochastic=False, output_interval=0.05,
                          initial_transfers=((0, 0, 1.0, 0.3),))
    assert low.V.max() < 0.25
    assert high.V.max() > 2.0


def test_phase_field_short_run_is_finite_and_bounded():
    result = simulate_phase_field(8e-4, dx=0.2, dt=8e-5,
                                  box_half_width=6.5, stochastic=False,
                                  output_interval=8e-4)
    assert np.isfinite(result.phi).all()
    assert np.isfinite(result.U).all()
    assert np.isfinite(result.V).all()
    assert result.phi.min() >= 0.0
    assert result.phi.max() <= 1.0
