import unittest
import numpy as np
from scipy.linalg import expm
from design import *
from reduced import *
from spatial import halfmass

class NumericalChecks(unittest.TestCase):
    def test_halfmass_reflection_and_flat(self):
        x=(np.arange(100)+.5)/100
        profiles=np.array([np.ones(100),np.exp(-20*x),np.exp(-20*(1-x))])
        p=halfmass(profiles)
        self.assertAlmostEqual(p[0],0.)
        self.assertAlmostEqual(p[1],-p[2])

    def test_continuous_cue_projection(self):
        x=(np.arange(20000)+.5)/20000
        for m in MODELS:
            for assay in ASSAYS:
                p=Protocol(m,assay)
                for t in [p.pulse/3,p.condition+.25*p.unit,p.condition+2.1*p.unit]:
                    z=p(t,x,1.)
                    self.assertAlmostEqual(z[10000:].mean()-z[:10000].mean(),p.projected(t),places=7)

    def test_memory_analytical_step(self):
        theta=np.array([.7,2.,.8,.1,.1,.2,0.,.1])
        spec=np.array([1.,100.,.4,.2,.1,1.])
        times=np.linspace(0,10,101)
        A=jacobian(theta,'memory');u=np.array([theta[5]*cue(0.,spec),0.])
        expected=np.array([np.linalg.solve(A,(expm(A*t)-np.eye(2))@u) for t in times])
        np.testing.assert_allclose(integrate(theta,times,spec,1),expected,atol=2e-7)
        self.assertTrue(np.all(np.linalg.eigvals(A)<0))

    def test_self_zero_nested(self):
        theta=np.array([.7,2.,3.,.1,.2,.2,0.,.1])
        t=np.linspace(0,10,100);s=np.array([0.,10.,.4,.2,.1,1.])
        np.testing.assert_array_equal(integrate(theta,t,s,2),integrate(theta,t,s,3))

    def test_latent_rescaling_preserves_old_model(self):
        # Old dq_old=-q_old/tm+alpha*p/(1+abs(p)/ps).
        tp,tm,alpha,beta,ps,vs=.7,2.,.3,.4,.1,.2
        factor=alpha*tm
        theta=np.array([tp,tm,beta*factor*tp,ps,vs/factor,.5,0.,.1])
        p,q=.13,-.09
        dp,dq=drift(p,q/factor,theta,2,.2)
        self.assertAlmostEqual(dp,-p/tp+beta*q/(1+abs(q)/vs)+.1)
        self.assertAlmostEqual(dq*factor,-q/tm+alpha*p/(1+abs(p)/ps))

    def test_stationary_memory_closure_and_camera_variance(self):
        from literature import modes,displacement_variance,linear_tracks
        theta=np.array([.7,2.,.8,.1,.1,.2,0.,.1])
        tp,tm,r=theta[:3];B=tp/2
        A=jacobian(theta,'memory');S=np.diag([B,B*tp/(r*tm)])
        Q=-(A@S+S@A.T)
        self.assertGreater(np.linalg.eigvalsh(Q).min(),0)
        rates,w=modes(theta,'memory')
        for t in [.01,1.,20.]:
            self.assertAlmostEqual(float(np.sum(w*np.exp(-rates*t))),float((expm(A*t)@S)[0,0]))
        inc=linear_tracks(theta,'memory',.2,.5,2.,n=20000,seed=141)
        expected=2*.2**2*displacement_variance(2/.2,rates,w)
        empirical=np.sum(inc.sum(1)**2,axis=1)
        self.assertLess(abs(empirical.mean()-expected),5*empirical.std()/np.sqrt(len(empirical)))

    def test_nonlinear_simulator_ou_limit(self):
        from literature import nonlinear_bank
        theta=np.array([.7,2.,0.,.1,.1,.2,0.,.1])
        _,p=nonlinear_bank(theta,.02,.3,2**17,631,2**15)
        self.assertLess(abs(np.mean(np.sum(p*p,axis=1))-.3**2*.7),.1*.3**2*.7)
        correlation=np.sum(p[1:]*p[:-1])/np.sum(p[:-1]*p[:-1])
        self.assertLess(abs(correlation-np.exp(-16*.02/.7)),.05)

if __name__=='__main__':unittest.main()
