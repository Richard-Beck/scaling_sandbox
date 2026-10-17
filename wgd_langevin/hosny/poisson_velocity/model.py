"""Linear velocity relaxation driven by marked Poisson kicks; h and micrometres.

Between events dv/dt=-v/tau. Events have rate lam (/h); a kick's magnitude
is Rayleigh(scale=A) in um/h and its angle from +y is von Mises(kappa).
No mean drift parameter. All cells share parameters. Event integration is exact;
stationary initialization neglects only events older than 14*tau.
"""
import numpy as np
from numba import njit
from scipy.special import i0e,i1e
from scipy.integrate import cumulative_trapezoid

def angle_table(kappa):
    if kappa<0:raise ValueError('Negative concentration')
    x=np.linspace(-np.pi,np.pi,32769)
    c=cumulative_trapezoid(np.exp(kappa*(np.cos(x)-1)),x,initial=0);c/=c[-1]
    return np.interp(np.linspace(0,1,32769),c,x)

@njit(cache=True)
def kick(table):
    z=np.random.random()*(len(table)-1);i=min(int(z),len(table)-2)
    theta=table[i]+(z-i)*(table[i+1]-table[i])
    mag=np.sqrt(-2*np.log(max(np.random.random(),1e-300)))
    return mag*np.sin(theta),mag*np.cos(theta)

@njit(cache=True)
def _simulate(tau,lam,A,table,n,cadence,frames,seed,warmup):
    np.random.seed(seed)
    p=np.zeros((n,frames+1,2));v=np.zeros_like(p)
    for cell in range(n):
        vx=0.;vy=0.
        age=np.random.exponential(1/lam)
        while True:
            jx,jy=kick(table);decay=np.exp(-age/tau);vx+=jx*decay;vy+=jy*decay
            if age>warmup*tau:break
            age+=np.random.exponential(1/lam)
        v[cell,0,0]=vx*A;v[cell,0,1]=vy*A
        t=0.;px=0.;py=0.;event=np.random.exponential(1/lam)
        for frame in range(1,frames+1):
            end=frame*cadence
            while event<end:
                dt=event-t;decay=np.exp(-dt/tau);integ=-tau*np.expm1(-dt/tau)
                px+=vx*integ;py+=vy*integ;vx*=decay;vy*=decay
                jx,jy=kick(table);vx+=jx;vy+=jy;t=event
                event+=np.random.exponential(1/lam)
            dt=end-t;decay=np.exp(-dt/tau);integ=-tau*np.expm1(-dt/tau)
            px+=vx*integ;py+=vy*integ;vx*=decay;vy*=decay;t=end
            p[cell,frame,0]=px*A;p[cell,frame,1]=py*A
            v[cell,frame,0]=vx*A;v[cell,frame,1]=vy*A
    return p,v

def simulate(tau,lam,A,kappa,n=4096,cadence=.5,duration=16,seed=220926,warmup=14,return_velocity=False):
    if tau<=0 or lam<=0 or A<0:raise ValueError('Invalid parameters')
    frames=int(round(duration/cadence))
    if cadence<=0 or duration<=0 or not np.isclose(frames*cadence,duration):raise ValueError('Invalid observation schedule')
    p,v=_simulate(tau,lam,A,angle_table(kappa),n,cadence,frames,seed,warmup)
    return (p,v) if return_velocity else p

def moments(tau,lam,A,kappa):
    ratio=i1e(kappa)/i0e(kappa)
    return lam*tau*A*np.sqrt(np.pi/2)*ratio,lam*tau*A*A

def msd(t,tau,lam,A,kappa):
    t=np.asarray(t);mean,trcov=moments(tau,lam,A,kappa)
    return (mean*t)**2+2*trcov*tau*(t+tau*np.expm1(-t/tau))

def metrics(p):
    path=np.linalg.norm(np.diff(p,axis=1),axis=2).sum(axis=1)
    end=p[:,-1]-p[:,0];dist=np.linalg.norm(end,axis=1);denom=np.maximum(path,1e-300)
    return dict(path=path,endpoint=dist,directness=dist/denom,parallel=end[:,1]/denom,
                perpendicular=end[:,0]/denom,angle_parallel=end[:,1]/np.maximum(dist,1e-300),angle=np.arctan2(end[:,0],end[:,1]))

def rayleigh_p(angles):
    n=angles.shape[-1];r=np.abs(np.mean(np.exp(1j*angles),axis=-1));z=n*r*r
    return np.clip(np.exp(-z)*(1+(2*z-z*z)/(4*n)-(24*z-132*z*z+76*z**3-9*z**4)/(288*n*n)),0,1)
