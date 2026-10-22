"""Audit modal tangent closure on the dominant eigenvector rather than zero entries."""
from pathlib import Path
import json,sys,numpy as np
from debelly_controls import n,r,OUT,base
keys=['rac','rho','mca_bound','mca_unbound','displacement','velocity']
d=r.read(OUT/'linear_mode_growth.json');checks=[]
for L in [18.,23.,39.]:
 row=next(a for a in d['rows'] if a['length']==L)['modes'][0];J=np.array(row['matrix']);eig,V=np.linalg.eig(J);j=np.argmax(abs(eig));lam=eig[j].real;v=V[:,j].real
 x=n.rd.Grid1D(L,201).x;c=np.cos(np.pi*x/L);s=np.sin(np.pi*x/L);eps=1e-6;b=base(L);traces=[]
 for sign in [-1,1]:
  state={k:y.copy() for k,y in b.items()}
  for key,a in zip(keys,v):state[key]+=sign*eps*a*(s if key in ['displacement','velocity'] else c)
  traces.append(n.run('debelly',L,state,n.ConstantInput(0.),.5,interval=.5).state)
 expected=[];actual=[]
 for key,a in zip(keys,v):
  q=s if key in ['displacement','velocity'] else c;expected.append(lam*a*q);actual.append((traces[1][key]-traces[0][key])/(2*eps))
 expected=np.array(expected);actual=np.array(actual)
 z=np.load(OUT/f'L{L:g}_homogeneous_seed.npz');t=z['time'];m=z['field_mca_bound'];a=m@c/(c@c)
 mask=(t>=1000)&(t<=4000);rate=np.polyfit(t[mask],np.log(abs(a[mask])),1)[0]
 checks.append(dict(length=L,dominant_growth_rate=row['dominant_growth_rate'],dominant_mode_full_native_step_relative_residual=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected)),measured_seed_growth_rate_1000_to_4000=float(rate),note='Large per-entry relative residuals in linear_mode_growth.json involve near-zero tangent components; this full eigenvector residual is the meaningful modal closure check.'))
r.write(OUT/'linear_mode_validation.json',checks);print(json.dumps(checks,indent=2))
