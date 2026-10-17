"""Active Brownian point cells with local oxygen-directed steering.

Heading measured from +y; dtheta=-kappa*sin(theta)dt+sqrt(2/tau)dW.
A stationary von Mises heading distribution has concentration a=kappa*tau.
All cells within a fitted condition share propulsion, tau and kappa.
"""
import numpy as np
from numba import njit

@njit(cache=True)
def _integrate(v,tau,kappa,cadence,duration,initial,substeps,seed):
    np.random.seed(seed)
    n=len(initial);frames=int(round(duration/cadence));dt=cadence/substeps
    out=np.zeros((n,frames+1,2))
    noise=np.sqrt(dt/tau)
    for cell in range(n):
        theta=initial[cell];x=0.;y=0.
        for frame in range(frames):
            for sub in range(substeps):
                theta += -.5*kappa*np.sin(theta)*dt + noise*np.random.randn()
                x += v*np.sin(theta)*dt
                y += v*np.cos(theta)*dt
                theta += -.5*kappa*np.sin(theta)*dt + noise*np.random.randn()
            out[cell,frame+1,0]=x;out[cell,frame+1,1]=y
    return out

def simulate(v,tau,a,n=2048,cadence=.5,duration=16,seed=220926,resolution=1,initial='stationary'):
    if v<0 or tau<=0 or a<0: raise ValueError('Invalid parameters')
    rng=np.random.default_rng(seed)
    angles=rng.vonmises(0,a,n) if initial=='stationary' else rng.uniform(-np.pi,np.pi,n)
    kappa=a/tau
    sub=int(np.ceil(max(24,10*cadence/tau,10*cadence*kappa)*resolution))
    return _integrate(v,tau,kappa,cadence,duration,angles,sub,seed+1)

def metrics(p):
    length=np.linalg.norm(np.diff(p,axis=1),axis=2).sum(axis=1)
    displacement=p[:,-1]-p[:,0];distance=np.linalg.norm(displacement,axis=1)
    safe=np.maximum(length,1e-15);den=np.maximum(distance,1e-15)
    return {'path':length,'endpoint':distance,'directness':distance/safe,
            'parallel':displacement[:,1]/safe,'perpendicular':displacement[:,0]/safe,
            'angle_parallel':displacement[:,1]/den,'angle_perpendicular':displacement[:,0]/den,
            'angle':np.arctan2(displacement[:,0],displacement[:,1])}

def rayleigh_p(angles):
    """Conventional finite-n approximation; only diagnostic at n=60."""
    n=angles.shape[-1]
    r=np.abs(np.mean(np.exp(1j*angles),axis=-1));z=n*r*r
    p=np.exp(-z)*(1+(2*z-z*z)/(4*n)-(24*z-132*z*z+76*z**3-9*z**4)/(288*n*n))
    return np.clip(p,0,1)
