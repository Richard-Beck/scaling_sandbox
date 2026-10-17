"""Effective zero-cue fixed points and scalar local response times, not turning times."""
import csv
import numpy as np
from scipy.optimize import brentq
from design import *
from reduced import drift,jacobian
from scaling import parameters

def build():
    rows=[]
    for source in sorted((OUT/'scaling').glob('*.json')):
        d=json.loads(source.read_text());family=d['family'];model=d['model']
        law=next(r for r in d['laws']if r['law']==d['selected_law'])
        for L in np.linspace(min(lengths(model)),max(lengths(model)),41):
            theta=parameters(law,L);roots=[0.]
            if family in ['coupled','self_activation']:
                def f(p):return drift(p,p/(1+p/theta[3]),theta,FAMILIES.index(family),0.)[0]
                grid=np.geomspace(1e-8,100.,300)
                for a,b in zip(grid[:-1],grid[1:]):
                    if f(a)*f(b)<0:roots.append(brentq(f,a,b))
            for p in roots:
                q=p/(1+p/theta[3]) if family in ['coupled','self_activation']else p
                eigen=np.array([-1/theta[0]])if family=='gaussian'else np.linalg.eigvals(jacobian(theta,family,p,q))
                stable=bool(np.max(eigen.real)<-1e-9)
                unit=clocks(model)[0]
                row=dict(model=model,family=family,law=law['law'],length=float(L),p_fixed=p,q_fixed=q,
                         stable=stable,physical_p=abs(p)<=.5,
                         fast_relaxation_s=float(unit/(-eigen.real).max())if stable else None,
                         slow_relaxation_s=float(unit/(-eigen.real).min())if stable else None,
                         oscillatory=bool(np.max(np.abs(eigen.imag))>1e-8),
                         tp_s=float(theta[0]*unit),tm_s=float(theta[1]*unit),feedback_ratio=float(theta[2]))
                rows.append(row)
    if rows:
        with (OUT/'effective_dynamics.csv').open('w',newline='')as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows

if __name__=='__main__':build()
