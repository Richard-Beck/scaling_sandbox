"""Explicit diagnostic controls; canonical response curves are unchanged."""
from pathlib import Path
import json,sys
from dataclasses import replace
import numpy as np
CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE));import gradient_reversal as g
r,n,e=g.r,g.n,g.e
OUT=n.DATA/'runs/fixed_slope_mechanisms_20261006/debelly'

def base(L):
    p=n.model_for('debelly').parameters;R,H=n.ds._homogeneous_rac_rho(p,.5)
    z=np.zeros(201);o=np.ones(201)
    return dict(rac=o*R,rho=o*H,mca_bound=o*.1*5/6.4,mca_unbound=o*.1*1.4/6.4,displacement=z.copy(),velocity=z.copy(),tension=z.copy(),raw_left_displacement=z.copy(),raw_right_displacement=z.copy(),left_boundary_displacement=z.copy(),right_boundary_displacement=z.copy())

def eigens():
    result=[];keys=['rac','rho','mca_bound','mca_unbound','displacement','velocity']
    for L in np.arange(16.,40.):
        row={'length':L,'modes':[]}
        for mode in range(1,5):
            x=n.rd.Grid1D(L,201).x;c=np.cos(mode*np.pi*x/L);s=np.sin(mode*np.pi*x/L)
            b=base(L);J=np.zeros((6,6));err=0.
            for j,key in enumerate(keys):
                eps=1e-6 if key in ['rac','rho','displacement','velocity'] else 1e-7
                ts=[]
                for sign in [-1,1]:
                    state={k:v.copy() for k,v in b.items()};state[key]+=sign*eps*(s if key in ['displacement','velocity'] else c)
                    tr=n.run('debelly',L,state,n.ConstantInput(0.),.5,interval=.5)
                    ts.append(tr.state)
                for i,target in enumerate(keys):
                    diff=(ts[1][target]-ts[0][target])/(2*eps)
                    basis=s if target in ['displacement','velocity'] else c
                    J[i,j]=np.dot(diff,basis)/np.dot(basis,basis)
                    err=max(err,float(np.linalg.norm(diff-J[i,j]*basis)/(1e-12+np.linalg.norm(diff))))
            eig,vec=np.linalg.eig(J);idx=np.argmax(abs(eig));v=np.real(vec[:,idx]);v*=1 if v[0]>0 else -1
            rates=np.log(eig.astype(complex))/.5
            row['modes'].append(dict(mode=mode,dominant_growth_rate=float(np.max(rates.real)),eigenvalues=eig.real.tolist(),dominant_vector={k:float(t) for k,t in zip(keys,v)},matrix=J.tolist(),max_projection_residual=err))
        result.append(row)
    r.write(OUT/'linear_mode_growth.json',dict(generated_at=n.now(),definition='Native half-unit step tangent projected onto Neumann cosine chemistry and Dirichlet sine mechanics, about homogeneous basal rest. Matrix is a diagnostic local linearization; no manuscript checkpoints replaced.',rows=result))
    for row in result:print(row['length'],[round(q['dominant_growth_rate'],7) for q in row['modes']])

def worker(slot):
    # Slots 0-5 continue saved endpoints an extra 20,000 basal units.
    if slot<6:
        job=r.read(OUT/'jobs.json')[slot];L=job['case']['length'];label=job['label']
        with np.load(OUT/f'L{L:g}_{label}.npz') as z:
            state={k:z['recovery_field_'+k][-1].copy() for k in n.DB_NAMES}
        tr=n.run('debelly',L,state,n.ConstantInput(0.),20000.,interval=10.)
        path=OUT/f'L{L:g}_{label}_continuation.npz'
    else:
        L=[18.,23.,39.][slot-6];state=base(L);q=np.cos(np.pi*n.rd.Grid1D(L,201).x/L)
        state['mca_bound']+=1e-5*q
        # Small conserved zero-mean adhesion perturbation, original basal equations.
        tr=n.run('debelly',L,state,n.ConstantInput(0.),10000.,interval=10.)
        path=OUT/f'L{L:g}_homogeneous_seed.npz'
    arrays={k:getattr(tr,k) for k in ['x','time','orientation','contrast']}
    arrays.update({'field_'+k:v for k,v in tr.fields.items()});np.savez_compressed(path,**arrays)
    r.write(path.with_suffix('.json'),dict(generated_at=n.now(),length=L,diagnostic='Basal continuation beyond fixed assay horizon' if slot<6 else 'Small cosine bound-MCA perturbation about homogeneous basal rest; does not replace any exact manuscript preparation',sha256=n.sha(path),final_orientation=float(tr.orientation[-1])))
    print(path,float(tr.orientation[-1]))

if __name__=='__main__':
    if sys.argv[1]=='eigens':eigens()
    else:worker(int(sys.argv[2]))
