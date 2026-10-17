import unittest
import numpy as np
from simulate import DimensionalParameters, Parameters, drift, simulate


class ModelChecks(unittest.TestCase):
    def test_vector_self_activation_is_radial(self):
        from dataclasses import replace
        params = DimensionalParameters(a=.02, K=.1, n=2)
        p = np.array([[.06, .08], [0., 0.]])
        v = np.zeros_like(p)
        dp, _ = drift(p, v, params)
        base, _ = drift(p, v, replace(params, a=0))
        np.testing.assert_allclose(dp-base, [[.006, .008], [0., 0.]], atol=1e-15)

    def test_dimensional_drift_matches_change_of_units(self):
        d = DimensionalParameters(tau_p=2., tau_v=7., alpha=0.8, beta=1.3, ps=3., vs=4.)
        n = Parameters(gp=d.beta*d.tau_p*d.vs/d.ps,
                       gv=d.alpha*d.tau_v*d.ps/d.vs, r=d.tau_p/d.tau_v)
        p, v = np.array([[0.8, -0.3], [0.2, 1.1]]), np.array([[1.2, 0.1], [-0.5, 0.4]])
        dp, dv = drift(p, v, d)
        dP, dV = drift(p/d.ps, v/d.vs, n)
        np.testing.assert_allclose(dp, dP*d.ps/d.tau_p)
        np.testing.assert_allclose(dv, dV*d.vs/d.tau_p)

    def test_ungated_velocity_timescale(self):
        p, v = np.array([[1., 0.]]), np.array([[0.25, 0.]])
        params = Parameters(r=2.)
        dp, dv = drift(p, v, params)
        np.testing.assert_allclose(dp, [[-0.4, 0.]])
        np.testing.assert_allclose(dv, [[2.5, 0.]])

    def test_disagreement_does_not_suppress_feedback(self):
        dp, dv = drift(np.array([[1., 0.]]), np.array([[-1., 0.]]), DimensionalParameters())
        np.testing.assert_allclose(dp, [[-1.5, 0.]])
        np.testing.assert_allclose(dv, [[0.6, 0.]])

    def test_uncoupled_deterministic_decay_and_position(self):
        p0, v0 = np.array([[1., -2.]]), np.array([[2., 1.]])
        dt, steps, r = 0.01, 100, 2.
        _, p, v, x = simulate(Parameters(gp=0., gv=0., r=r), p0, v0, 0., 0., 1., dt, 1)
        np.testing.assert_allclose(p[-1], p0 * (1-dt)**steps)
        np.testing.assert_allclose(v[-1], v0 * (1-r*dt)**steps)
        np.testing.assert_allclose(x[-1], v0 * (1-(1-r*dt)**steps)/r)

    def test_noise_scaling_and_independence(self):
        initial = np.zeros((30000, 2))
        _, p, v, _ = simulate(Parameters(r=3.), initial, initial, 0.4, 0.7, 0.02, 0.02, 18)
        self.assertAlmostEqual(np.var(p[-1]) / (0.4**2 * 0.02), 1., delta=0.025)
        self.assertAlmostEqual(np.var(v[-1]) / (0.7**2 * 0.02), 1., delta=0.025)
        self.assertLess(abs(np.corrcoef(p[-1].ravel(), v[-1].ravel())[0, 1]), 0.02)


if __name__ == '__main__':
    unittest.main()
