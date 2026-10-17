import numpy as np

from polarity1d.contract import Grid1D
from polarity1d.holmes import HolmesModel3


def test_original_holmes_eq5_eq6_factor_two_and_pool_scaling():
    model = HolmesModel3(rac_activation=2.0)
    rac, rho = 1.0, 1.25
    state = np.array([[rac] * 5,
                      [model.rac_total - rac] * 5,
                      [rho] * 5,
                      [model.rho_total - rho] * 5])
    derivative = model.reaction(0.0, np.arange(5.0), state, np.zeros(5))

    # At rho=a1 and Rac=a2, each source-paper Hill denominator is 2*(1+1)=4.
    rac_rate = (derivative[0, 0] + model.rac_decay * rac) / (
        (model.rac_total - rac) / model.rac_total
    )
    rho_rate = (derivative[2, 0] + model.rho_decay * rho) / (
        (model.rho_total - rho) / model.rho_total
    )
    assert np.isclose(rac_rate, 2.0 / 4.0)
    assert np.isclose(rho_rate, 6.6 / 4.0)
    assert np.allclose(derivative[0] + derivative[1], 0.0)
    assert np.allclose(derivative[2] + derivative[3], 0.0)


def test_holmes_model3_initial_homogeneous_state_is_no_cue_equilibrium():
    model = HolmesModel3()
    grid = Grid1D(20.0, 11)
    state = model.initial_state(grid)
    derivative = model.reaction(0.0, grid.x, state, np.zeros(grid.cells))
    assert np.max(np.abs(derivative)) < 1e-10
