import unittest
import numpy as np
from fit import *
class RefinedTests(unittest.TestCase):
    def test_integrated_ou_units_and_limits(self):
        tau=.3;sigma=7.;t=np.array([1e-8,.02,1.,20.])
        v=displacement_variance(t,np.array([1/tau]),np.array([sigma*sigma*tau/2]))
        self.assertAlmostEqual(v[0]/(sigma*sigma*tau*t[0]**2/2),1.,places=6)
        expected=sigma*sigma*tau*tau*(t[1:]-tau*(-np.expm1(-t[1:]/tau)))
        np.testing.assert_allclose(v[1:],expected,rtol=1e-12)
    def test_curve_contribution_balanced(self):
        self.assertEqual(len(CURVE_OBS),510)
        for g in (0,1):
            for b in range(3):
                self.assertAlmostEqual(sum(o['objective_weight'] for o in CURVE_OBS if o['group']==g and o['lag_bin']==b),1/3)
    def test_curve_prediction_uses_each_exact_lag(self):
        z=initial(0);p=analytic(z);pars,_,_=decode(z)
        _,rows=data_residuals(p)
        for row in rows:
            if row['metric']!='msd_curve':continue
            tau,sigma=pars[0,row['group']]
            expected=2*float(displacement_variance(row['lag_h'],np.array([1/tau]),np.array([sigma*sigma*tau/2])))
            self.assertAlmostEqual(row['predicted']/expected,1.,places=12)
    def test_signed_drift_tails_finite(self):
        z=initial(0);z[20:]=[-200,200]
        r=residual(z,'primary')
        self.assertTrue(np.isfinite(r).all())
if __name__=='__main__':unittest.main()
