"""Exactly sampled integrated 2D OU velocity, lengths um and times h.

 dr = V dt
 dV = -(V - u e_y)dt/tau + sqrt(2 sigma**2/tau) dW
 sigma is stationary per-axis velocity SD; u is oxygen-directed mean velocity.
 Initial velocities are stationary. Each cell shares the same parameters.
"""
import numpy as np

def transition(tau,dt,sigma=1.):
    if tau<=0 or dt<=0 or sigma<0:raise ValueError('Invalid OU parameters')
    q=dt/tau;e=-np.expm1(-q);a=1-e;b=tau*e
    vv=sigma**2*e*(2-e)
    xv=sigma**2*tau*e**2
    xx=sigma**2*tau**2*(2*q-2*e-e**2)
    return a,b,vv,xv,max(xx,0.)

def simulate(sigma,tau,u,n=8192,cadence=.5,duration=16,seed=230926,initial='stationary',return_velocity=False):
    if min(sigma,u)<0:raise ValueError('sigma and u must be nonnegative')
    frames=int(round(duration/cadence))
    if not np.isclose(frames*cadence,duration):raise ValueError('Nonintegral frame count')
    a,b,vv,xv,xx=transition(tau,cadence,sigma)
    rng=np.random.default_rng(seed);w=sigma*rng.standard_normal((n,2)) if initial=='stationary' else np.zeros((n,2))
    p=np.zeros((n,frames+1,2));velocity=np.zeros_like(p) if return_velocity else None
    if return_velocity:velocity[:,0]=w;velocity[:,0,1]+=u
    sv=np.sqrt(vv);cross=xv/sv if sv else 0.;sr=np.sqrt(max(xx-cross**2,0.))
    for k in range(frames):
        z=rng.standard_normal((n,2));z2=rng.standard_normal((n,2))
        dr=b*w+cross*z+sr*z2;dr[:,1]+=u*cadence
        p[:,k+1]=p[:,k]+dr
        w=a*w+sv*z
        if return_velocity:velocity[:,k+1]=w;velocity[:,k+1,1]+=u
    return (p,velocity) if return_velocity else p

def msd(t,sigma,tau,u):
    t=np.asarray(t,dtype=float)
    return (u*t)**2+4*sigma**2*tau**2*(t/tau+np.expm1(-t/tau))

def metrics(p):
    step=np.diff(p,axis=1);path=np.linalg.norm(step,axis=2).sum(axis=1)
    end=p[:,-1]-p[:,0];dist=np.linalg.norm(end,axis=1)
    denom=np.maximum(path,1e-15)
    return dict(path=path,endpoint=dist,directness=dist/denom,parallel=end[:,1]/denom,
                perpendicular=end[:,0]/denom,angle_parallel=end[:,1]/np.maximum(dist,1e-15),
                angle=np.arctan2(end[:,0],end[:,1]))

def rayleigh_p(angles):
    n=angles.shape[-1];r=np.abs(np.mean(np.exp(1j*angles),axis=-1));z=n*r*r
    p=np.exp(-z)*(1+(2*z-z*z)/(4*n)-(24*z-132*z*z+76*z**3-9*z**4)/(288*n*n))
    return np.clip(p,0,1)
