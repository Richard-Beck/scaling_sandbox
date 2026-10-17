import unittest
import numpy as np
from simulate import DimensionalParameters, Parameters, drift, gate, simulate


class ModelChecks(unittest.TestCase):
    def test_dimensional_drift_matches_change_of_units(self):
        d = DimensionalParameters(tau_p=2., tau_v=7., alpha=0.8, beta=1.3, ps=3., vs=4., k=0.6, z0=0.9)
        n = Parameters(gp=d.beta*d.tau_p*d.vs/d.ps,
                       gv=d.alpha*d.tau_v*d.ps/d.vs, r=d.tau_p/d.tau_v,
                       K=d.k*d.ps*d.vs, Z0=d.z0/(d.ps*d.vs))
        p, v = np.array([[0.8, -0.3], [0.2, 1.1]]), np.array([[1.2, 0.1], [-0.5, 0.4]])
        dp, dv = drift(p, v, d)
        dP, dV = drift(p/d.ps, v/d.vs, n)
        np.testing.assert_allclose(dp, dP*d.ps/d.tau_p)
        np.testing.assert_allclose(dv, dV*d.vs/d.tau_p)

    def test_gate_and_velocity_timescale(self):
        p, v = np.array([[1., 0.]]), np.array([[0.25, 0.]])
        params = Parameters(r=2.)
        np.testing.assert_allclose(gate(p, v, params), 0.5)
        dp, dv = drift(p, v, params)
        np.testing.assert_allclose(dp, [[-0.7, 0.]])
        np.testing.assert_allclose(dv, [[1., 0.]])

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
