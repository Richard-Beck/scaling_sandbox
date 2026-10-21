"""Scientific checks for the supplementary saved-interval heading rates."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure3_reader_20261020 import interval_rates


class RateTests(unittest.TestCase):
    def test_constant_rotation_recovers_native_angular_rate_across_angle_seam(self):
        t=np.arange(10,dtype=float)*.5
        middle,rates=interval_rates(t,.5*np.exp(1j*(2.9+.2*t)))
        np.testing.assert_allclose(middle,.5*(t[:-1]+t[1:]))
        for rate in rates.values():np.testing.assert_allclose(rate,.2,atol=1e-14)

    def test_opposite_visible_endpoints_do_not_invent_rate_across_cancellation(self):
        _,rates=interval_rates([0.,1.,2.],np.array([.5,-.5,-.5j]))
        for rate in rates.values():
            self.assertTrue(np.isnan(rate[0]))
            self.assertAlmostEqual(rate[1],np.pi/2)

    def test_quality_gates_keep_raw_direction_and_exclude_weak_reporter(self):
        _,rates=interval_rates([0.,1.],.03*np.exp(1j*np.array([0.,.1])))
        self.assertAlmostEqual(rates['1e-10'][0],.1)
        for gate in ['0.05','0.1','0.2']:self.assertTrue(np.isnan(rates[gate][0]))


if __name__=='__main__':unittest.main()
