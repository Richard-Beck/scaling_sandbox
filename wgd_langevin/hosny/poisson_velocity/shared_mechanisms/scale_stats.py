"""Same representative-model loss with amplitude quadratic coefficients exposed.
Only for shared amplitude across baseline and both KDs (mask 8/8).
"""
import numpy as np
from numba import njit

@njit(cache=True)
def score_and_scale(z,mapping,ampids,bank,xgrid,ygrid,kgrid,targets,scales,owners):
    qq=np.zeros((3,3,3))
    for g in range(3):
        x=z[mapping[g,0]];y=z[mapping[g,1]];k=z[mapping[g,2]]
        ix=min(max(np.searchsorted(xgrid,x)-1,0),len(xgrid)-2)
        iy=min(max(np.searchsorted(ygrid,y)-1,0),len(ygrid)-2)
        ik=min(max(np.searchsorted(kgrid,k)-1,0),len(kgrid)-2)
        ax=(x-xgrid[ix])/(xgrid[ix+1]-xgrid[ix]);ay=(y-ygrid[iy])/(ygrid[iy+1]-ygrid[iy]);ak=(k-kgrid[ik])/(kgrid[ik+1]-kgrid[ik])
        for a in range(2):
            for b in range(2):
                for c in range(2):
                    weight=(ax if a else 1-ax)*(ay if b else 1-ay)*(ak if c else 1-ak)
                    for m in range(3):
                        for q in range(3):qq[g,m,q]+=weight*bank[ix+a,iy+b,ik+c,m,q]
    num=np.zeros(3);den=np.zeros(3)
    for d in range(len(owners)):
        g=owners[d];aid=ampids[g]
        for q in range(3):
            num[aid]+=qq[g,0,q]*targets[d,0,q]/scales[d,0]**2
            den[aid]+=qq[g,0,q]**2/scales[d,0]**2
    amps=np.ones(3)
    for a in range(3):
        if den[a]>0:amps[a]=num[a]/den[a]
    total=0.
    for d in range(len(owners)):
        g=owners[d]
        for m in range(3):
            for q in range(3):
                pred=qq[g,m,q]*(amps[ampids[g]] if m==0 else 1.)
                total+=((pred-targets[d,m,q])/scales[d,m])**2
    return total/(len(owners)*9),num[0],den[0]

