import unittest
import numpy as np
from campaign import integrate, pulses_from_task, reversal, settling


class ReductionChecks(unittest.TestCase):
    def test_exact_uncoupled_decay(self):
        theta = np.array([3., 7., 0., 0., 1., 1., 1., 1., 0.])
        t = np.linspace(0, 30, 101)
        y = integrate(theta, t, np.empty((0, 4)), np.array([.2, -.3]), .1)
        np.testing.assert_allclose(y[:, 0], .2*np.exp(-t/3), atol=1e-8)
        np.testing.assert_allclose(y[:, 1], -.3*np.exp(-t/7), atol=1e-8)

    def test_event_alignment_and_reflection(self):
        theta = np.array([3., 7., 0., 0., 1., 1., 1., 1., .1])
        t = np.array([0., .9, 1.3, 2.])
        pulse = np.array([[.13, .77, 0., 2.]])
        y = integrate(theta, t, pulse, np.zeros(2), .05)
        expected = .2*3*(1-np.exp(-(.77-.13)/3))*np.exp(-(t[-1]-.77)/3)
        self.assertAlmostEqual(y[-1, 0], expected, places=8)
        pulse[:, 3] *= -1
        np.testing.assert_allclose(integrate(theta, t, pulse, np.zeros(2), .05), -y)

    def test_cue_projection(self):
        task = {'protocol': {'name': 'localized_pulse', 'amplitude': 3.6,
                             'width_fraction': .2, 'sides': ['left','right'],
                             'asymmetry': .15, 'start': 0., 'stop': 80., 'taper': 40.}}
        self.assertAlmostEqual(pulses_from_task(task)[0, 3], -.108)

    def test_reversal_and_confirmation(self):
        t = np.arange(0., 2001., 5.)
        y = np.where(t < 200, -.1, .1)
        self.assertTrue(reversal(t, y, 200., -.1)[0])
        self.assertFalse(reversal(t, y*.01, 200., -.001)[0])
        y[-1] = .15
        self.assertIsNone(settling(t, y, 200.)[0])


if __name__ == '__main__':
    unittest.main()
