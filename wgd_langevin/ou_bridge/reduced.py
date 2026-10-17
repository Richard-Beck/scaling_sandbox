"""Reduced deterministic models; time measured in each spatial model's clock.

q's scale is fixed by dq=(p/(1+abs(p)/ps)-q)/tm. This removes a latent-state
scaling degeneracy, not a feedback mechanism. q is not asserted to be MCA.
The linear memory realization is a response hypothesis, distinct from its
specified stationary substitute-kernel stochastic closure in literature.py.
"""
import numpy as np
from numba import njit
from design import FAMILIES

NAMES = ['tp', 'tm', 'r', 'ps', 'vs', 'gain', 'self_a', 'self_K']
LOW = np.array([.01, .01, 0., .003, .003, 1e-5, 0., .003])
HIGH = np.array([1000., 1000., 30., 2., 2., 20., 5., .5])
LOG = np.array([True, True, False, True, True, True, False, True])

def active(family):
    return {'gaussian': [0, 5], 'memory': [0, 1, 2, 5],
            'coupled': [0, 1, 2, 3, 4, 5], 'self_activation': list(range(8))}[family]

def bounds(family):
    lo, hi = LOW.copy(), HIGH.copy()
    if family == 'memory': hi[2] = .995
    lo[LOG], hi[LOG] = np.log(lo[LOG]), np.log(hi[LOG])
    return lo[active(family)], hi[active(family)]

def encode(theta, family):
    z = np.array(theta).copy()
    z[LOG] = np.log(z[LOG])
    return z[active(family)]

def decode(z, family):
    theta = np.array([1., 1., 0., .1, .1, .1, 0., .05])
    for i, value in zip(active(family), z):
        theta[i] = np.exp(value) if LOG[i] else value
    return theta

@njit(cache=True)
def cue(t, spec):
    assay, condition, pulse, peak, width, period = spec
    if assay == 0 or assay == 6:
        if t >= pulse: return 0.
        envelope = 1. if t <= pulse/2 else .5*(1+np.cos(np.pi*(t-pulse/2)/(pulse/2)))
        return -peak*width*envelope*(.02 if assay == 6 else 1.)
    direction = -1.
    if t >= condition:
        if assay == 3 or assay == 4:
            direction = -np.cos(2*np.pi*(t-condition)/period)
        elif assay == 5:
            direction = 1.
        elif t < condition+1.:
            direction += .05 if assay == 1 else -.05
    return .4*peak*direction

@njit(cache=True)
def drift(p, q, theta, family, u):
    tp, tm, r, ps, vs, gain, a, K = theta
    if family == 0: return -p/tp+gain*u, 0.
    if family == 1: return (-p+r*q)/tp+gain*u, (p-q)/tm
    return ((-p+r*q/(1+abs(q)/vs))/tp + a*p*abs(p)/(K*K+p*p)+gain*u,
            (p/(1+abs(p)/ps)-q)/tm)

@njit(cache=True)
def integrate(theta, times, spec, family, refinement=1.):
    result = np.zeros((len(times), 2))
    p = q = 0.
    t = times[0]
    dt = min(.1, .25*theta[0], .25*theta[1], .25/(theta[6]+1e-6))/refinement
    events = np.array([spec[2]/2, spec[2], spec[1], spec[1]+1.])
    for i in range(1, len(times)):
        while t < times[i]-1e-10:
            h = min(dt, times[i]-t)
            for event in events:
                if event > t+1e-10: h = min(h, event-t)
            a,b = drift(p,q,theta,family,cue(t+1e-11,spec))
            c,d = drift(p+h*a/2,q+h*b/2,theta,family,cue(t+h/2,spec))
            e,f = drift(p+h*c/2,q+h*d/2,theta,family,cue(t+h/2,spec))
            g,j = drift(p+h*e,q+h*f,theta,family,cue(t+h-1e-11,spec))
            p += h*(a+2*c+2*e+g)/6
            q += h*(b+2*d+2*f+j)/6
            t += h
            if abs(p)>1e6 or abs(q)>1e6:
                result[i:,:] = 1e6
                return result
        result[i,0],result[i,1] = p,q
    return result

def spec_for(protocol):
    from design import ASSAYS
    return np.array([ASSAYS.index(protocol.assay), protocol.condition/protocol.unit,
                     protocol.pulse/protocol.unit, protocol.peak, protocol.width,
                     protocol.period/protocol.unit])

def jacobian(theta, family, p=0., q=0.):
    fid = FAMILIES.index(family)
    eps = 1e-6
    base = np.array([p, q])
    return np.column_stack([(np.array(drift(*(base+eps*np.eye(2)[i]),theta,fid,0.))-
                             np.array(drift(*(base-eps*np.eye(2)[i]),theta,fid,0.)))/(2*eps)
                            for i in range(2)])
