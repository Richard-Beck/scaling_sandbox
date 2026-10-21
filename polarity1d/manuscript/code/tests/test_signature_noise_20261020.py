"""Distribution and reproducibility checks for the declared linear noise law."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / "mechanisms/signature_noise_20261020.py"
SPEC = importlib.util.spec_from_file_location("signature_noise", SOURCE)
NOISE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NOISE)


class LinearNoiseTests(unittest.TestCase):
    def test_stationary_marginals_and_nonperiodic_covariance(self):
        field = NOISE.LinearGapNoise(380., 200, 50000., 213, .01,
                                    ell=1.9, tau=5., knot=5.)
        values = field.knots
        self.assertAlmostEqual(float(values.mean()), 0., delta=.01)
        self.assertAlmostEqual(float(values.var()), 1., delta=.02)
        self.assertAlmostEqual(float(np.mean(values[:, :-1] * values[:, 1:])),
                               np.exp(-1), delta=.015)
        self.assertAlmostEqual(float(np.mean(values[:-1] * values[1:])),
                               np.exp(-1), delta=.015)
        # A periodic implementation would correlate these distant endpoints.
        self.assertAlmostEqual(float(np.mean(values[:, 0] * values[:, -1])),
                               0., delta=.05)
        self.assertAlmostEqual(float(values[:, 0].var()), 1., delta=.05)

    def test_interpolation_preserves_gaussian_variance(self):
        field = NOISE.LinearGapNoise(380., 200, 5000., 789, .01,
                                    ell=1.9, tau=5., knot=5.)
        samples = np.vstack([field.gaussian(t) for t in np.arange(2.5, 5000., 5.)])
        self.assertAlmostEqual(float(samples.mean()), 0., delta=.02)
        self.assertAlmostEqual(float(samples.var()), 1., delta=.025)
        lognormal = .01 * np.exp(samples - .5)
        self.assertAlmostEqual(float(lognormal.mean()), .01, delta=.0004)

    def test_seed_and_audit_reproduce_the_same_realization(self):
        first = NOISE.LinearGapNoise(19., 201, 100., 99, .001)
        second = NOISE.LinearGapNoise(19., 201, 100., 99, .001)
        self.assertEqual(first.audit()["gaussian_knots_sha256"],
                         second.audit()["gaussian_knots_sha256"])
        np.testing.assert_array_equal(first.knots, second.knots)

    def test_mean_and_zero_controls_do_not_consume_random_fields(self):
        x = (np.arange(201) + .5) * 19 / 201
        for kind, target in (("mean", .001), ("zero", 0.)):
            field = NOISE.LinearGapNoise(19., 201, 100., 5, .001, kind=kind)
            self.assertIsNone(field.knots)
            np.testing.assert_array_equal(field(50., x, 19.), np.full(201, target))

    def test_bounded_transform_has_stationary_mean_and_hard_limits(self):
        field = NOISE.LinearGapNoise(380., 200, 50000., 213, .02,
                                    ell=1.9, tau=5., knot=5., transform="bounded_tanh")
        for values in (field.knots[:2000],field.knots[4000:6000],field.knots[-2000:]):
            mapped=field.map_gaussian(values)
            self.assertAlmostEqual(float(mapped.mean()),.02,delta=.0003)
            self.assertGreaterEqual(float(mapped.min()),0.)
            self.assertLessEqual(float(mapped.max()),.04)
        extremes=field.map_gaussian(np.array([-1e9,0.,1e9]))
        np.testing.assert_array_equal(extremes,[0.,.02,.04])
        x=(np.arange(field.cells)+.5)*field.length/field.cells
        interpolated=np.vstack([field(t,x,field.length) for t in np.arange(2.5,10000.,5.)])
        self.assertAlmostEqual(float(interpolated.mean()),.02,delta=.0003)
        self.assertGreaterEqual(float(interpolated.min()),0.)
        self.assertLessEqual(float(interpolated.max()),.04)
        audit=field.audit()
        self.assertEqual(audit["theoretical_hard_input_range"],[0.,.04])
        self.assertAlmostEqual(audit["theoretical_input_mean"],.02)
        self.assertAlmostEqual(audit["transformed_knot_statistics"]["variance"],
                               audit["theoretical_input_variance"],delta=4e-6)

    def test_bounded_input_correlation_matches_transformed_gaussian_law(self):
        field = NOISE.LinearGapNoise(380.,200,50000.,213,.02,
                                    ell=1.9,tau=5.,knot=5.,transform="bounded_tanh")
        # Integrate the nonlinear covariance against a correlated unit-Gaussian
        # pair. The transformed input must not inherit the OU correlation label.
        nodes,weights=np.polynomial.hermite.hermgauss(60)
        z=np.sqrt(2.)*nodes;rho=np.exp(-1.)
        first=np.tanh(z[:,None])
        second=np.tanh(rho*z[:,None]+np.sqrt(1-rho*rho)*z[None,:])
        covariance=float(np.sum(weights[:,None]*weights[None,:]*first*second)/np.pi)
        variance=float(np.sum(weights*np.tanh(z)**2)/np.sqrt(np.pi))
        expected=covariance/variance
        stats=field.audit()["transformed_knot_statistics"]
        self.assertAlmostEqual(stats["adjacent_spatial_correlation"],expected,delta=.015)
        self.assertAlmostEqual(stats["adjacent_temporal_correlation"],expected,delta=.015)
        self.assertAlmostEqual(float(field.knots.var()),1.,delta=.02)


if __name__ == "__main__":
    unittest.main()
