import unittest
import numpy as np
from campaign import specification,decode,jobs,predictions,observations,score,ensemble,modes,displacement_variance
class TestCampaign(unittest.TestCase):
    def test_only_persistence_can_change(self):
        self.assertEqual(len(jobs()),72)
        for f,mask,s in jobs():
            names,lo,hi=specification(f,mask,s)
            self.assertEqual(set(k for k in names if k.endswith('_poly')),set() if not mask else {'tp_poly'} if f=='gaussian' else {'tp_poly','tm_poly'})
            z=(lo+hi)/2
            for j,k in enumerate(names):
                if k.endswith('_poly'):z[j]+=.2
            a=decode(z,names,0);b=decode(z,names,1)
            np.testing.assert_array_equal(a[0][2:],b[0][2:]);self.assertEqual(a[1:],b[1:])
    def test_shared_groups_exactly_equal(self):
        for f in ['gaussian','memory','coupled','self_activation']:
            names,lo,hi=specification(f,(),'Hosny_baseline');p=predictions((lo+hi)/2,names,f,'Hosny_baseline',n=24)
            self.assertEqual(p[0],p[1])
    def test_gaussian_analytic_speed_and_msd(self):
        z=np.log([.7]);p=predictions(z,['tp'],'gaussian','Xuan2018',n=16)[0]
        cad=1/6;t=11.1111111111;tau=.7
        var=lambda t:tau*tau*(t-tau*(-np.expm1(-t/tau)))
        self.assertAlmostEqual(p['speed'],np.sqrt(np.pi/2*var(cad))/cad)
        self.assertAlmostEqual(p['msd'],2*var(t))
    def test_memory_covariance_closure(self):
        from scipy.linalg import expm
        for tp,tm,r in [(.05,12.,.99),(12.,.05,.005),(.7,1.4,.6)]:
            theta=np.array([tp,tm,r,1.,1.,1.,0.,.1])
            a=np.array([[-1/tp,r/tp],[1/tm,-1/tm]])
            cov=np.diag([tp/2,tp*tp/(2*r*tm)])
            noise=-(a@cov+cov@a.T)
            self.assertGreaterEqual(np.linalg.eigvalsh(noise).min(),-1e-10)
            self.assertAlmostEqual(noise[0,0],1.)
            rates,weights=modes(theta,'memory')
            for t in (0.,.1,2.,12.):
                self.assertAlmostEqual(float((expm(a*t)@cov)[0,0]),float(np.sum(weights*np.exp(-rates*t))),places=9)
    def test_nonlinear_integrator_ou_limit(self):
        theta=np.array([.7,1.,0.,1.,1.,1.,0.,.1])
        inc,lag=ensemble(theta,1.,0.,.5,4.,4.,4096,987,.01,8.)
        expected=2*displacement_variance(4.,*modes(theta,'gaussian'))
        self.assertLess(abs(lag.mean()/expected-1),.065)
        self.assertLess(abs(np.mean(inc)),.03)
    def test_holdouts(self):
        for s in ['Xuan2018','Mallin2023','Hosny_HIF_KD']:
            fit=observations(s);check=observations(s,True)
            self.assertTrue(check)
            a={(r['metric'],r.get('lag_h')) for r in fit};b={(r['metric'],r.get('lag_h')) for r in check}
            self.assertFalse(a&b)
    def test_shared_amplitude(self):
        z=np.log([.7,2.]);p=predictions(z,['tp','tp_poly'],'gaussian','Jemaa2017',n=16)
        _,amp,rows,_=score(p,observations('Jemaa2017'),'Jemaa2017')
        for r in rows:self.assertAlmostEqual(r['predicted']/p[r['group']]['speed'],amp)
if __name__=='__main__':unittest.main()
