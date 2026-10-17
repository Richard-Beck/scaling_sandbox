import unittest
import numpy as np
from campaign import jobs,specification,physical,theta_values,predict,core,CAP
class OrderedTests(unittest.TestCase):
    def test_constraints(self):
        rng=np.random.default_rng(181)
        self.assertEqual(len(jobs()),72)
        for f,mask,study in jobs():
            names,lo,hi=specification(f,mask,study)
            for i in range(30):
                z=rng.uniform(lo,hi);d,pn,pz,t=theta_values(z,names)
                times=['tp'] if f=='gaussian' else ['tp','tm']
                for k in times:
                    self.assertTrue(0<d[k]<CAP)
                    if mask:self.assertTrue(d[k]<d[k+'_poly']<CAP)
                np.testing.assert_array_equal(t[0][0][2:],t[1][0][2:])
                self.assertEqual(t[0][1:],t[1][1:])
    def test_null_exact(self):
        for f in ('gaussian','memory','coupled','self_activation'):
            names,lo,hi=specification(f,(),'Hosny_baseline')
            d=dict(tp=.6,tm=2.,r=.6,vs=.5,eta=.5,self_a=.1,self_K=.2,force=.1)
            z=np.log([d[k] for k in names]);p=predict(z,names,f,'Hosny_baseline',n=12)
            self.assertEqual(p[0],p[1])
    def test_exact_linear_targets(self):
        for f in ('gaussian','memory'):
            names,lo,hi=specification(f,(0,),'Xuan2018');z=(lo+hi)/2
            d,pn,pz,_=theta_values(z,names)
            old=core.predictions(pz,pn,f,'Xuan2018',n=16);new=predict(z,names,f,'Xuan2018',n=16)
            for a,b in zip(old,new):
                self.assertAlmostEqual(a['speed'],b['speed']);self.assertAlmostEqual(a['msd'],b['msd'])
    def test_short_time_ou_numerics(self):
        theta=np.array([.001,.001,0.,1.,1.,1.,0.,.1])
        inc,lag=core.ensemble(theta,1.,0.,.05,.5,.5,2048,763,.000125,.05)
        target=2*core.displacement_variance(.5,*core.modes(theta,'gaussian'))
        self.assertLess(abs(lag.mean()/target-1),.07)
if __name__=='__main__':unittest.main()
