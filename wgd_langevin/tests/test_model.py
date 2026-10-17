"""Physical invariants and stochastic moment checks, independent of fitting."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model import msd, simulate, metrics

class ModelTests(unittest.TestCase):
    def test_ballistic_limit(self):
        p=simulate(3,1e18,.25,2,n=64,seed=11)
        m=metrics(p,.25)
        np.testing.assert_allclose(m['speed'],3,rtol=1e-8)
        np.testing.assert_allclose(m['directness'],1,rtol=1e-8)
        self.assertAlmostEqual(float(msd(1,3,1e5)),9,places=4)

    def test_active_msd_moment(self):
        p=simulate(4,.7,.25,4,n=8192,substeps=32,seed=13)
        empirical=np.sum(p[:,-1]**2,axis=1)
        expected=float(msd(4,4,.7))
        self.assertLess(abs(empirical.mean()-expected),5*empirical.std(ddof=1)/np.sqrt(len(p)))

    def test_diffusion_and_localization(self):
        p=simulate(0,1,.25,2,n=8192,substeps=4,seed=17,diffusion=2,localization_sd=.4)
        d=np.sum((p[:,-1]-p[:,0])**2,axis=1)
        expected=4*2*2+4*.4**2
        self.assertLess(abs(d.mean()-expected),5*d.std(ddof=1)/np.sqrt(len(p)))

    def test_heterogeneity_second_moment(self):
        p=simulate(2,.5,.5,4,n=8192,substeps=32,seed=23,cv=.7)
        d=np.sum(p[:,-1]**2,axis=1)
        expected=float(msd(4,2,.5,cv=.7))
        self.assertLess(abs(d.mean()-expected),5*d.std(ddof=1)/np.sqrt(len(p)))

    def test_invalid_duration(self):
        with self.assertRaises(ValueError): simulate(1,1,.3,1)

if __name__=='__main__': unittest.main()
