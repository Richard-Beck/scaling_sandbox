"""Native seven-model solver and checkpoint contract for corrected basal assays.

This module contains no experiment launcher or legacy preparation chain.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace
from datetime import datetime,timezone
import csv,hashlib,json,sys
import numpy as np
CODE=Path(__file__).resolve().parent
ROOT=CODE.parents[1]
DATA=ROOT/'data'
sys.path.insert(0,str(CODE/'src'))
from models import wave_pinning,otsuji,goryachev,legi,holmes_model3 as ho,debelly as db,spring as sp
from polarity1d.contract import Grid1D
from polarity1d.simulation import simulate
PROTOCOL=json.loads((CODE/'protocol.json').read_text())
MODELS=['wave_pinning','otsuji','goryachev','legi','debelly','holmes_model3','spring']
RD_MODELS=set(MODELS[:4]+['holmes_model3'])
DB_NAMES=['rac','rho','mca_bound','mca_unbound','displacement','velocity','tension','raw_left_displacement','raw_right_displacement','left_boundary_displacement','right_boundary_displacement']
def resolve_path(name):
    """Resolve immutable legacy record paths into this self-contained package."""
    name = str(name)
    maps = {
        'polarity1d/four_assay/work/': 'data/linear_v11/work/',
        'polarity1d/four_assay/results/': 'data/linear_v11/',
        'polarity1d/four_assay/pilots/width_full_lengths_20261002/': 'data/width_reversal/',
        'polarity1d/four_assay/pilots/width_mean_slope_20261002/': 'data/width_reversal/support/width_pilot/',
        'polarity1d/four_assay/pilots/full_reorientation_width_20261002/': 'data/width_reversal/support/full_reorientation/',
        'polarity1d/four_assay/pilots/polarity_establishment_ref15_relative20_basal_20261003/': 'data/establishment_basal/',
    }
    for before, after in maps.items():
        if name.startswith(before):
            return ROOT / (after + name[len(before):])
    return ROOT / name

def protect_retained_output(path):
    """Never mutate the original retained evidence with adapted campaign code."""
    path = Path(path).resolve()
    for label in ['linear_v11', 'width_reversal', 'establishment_basal', 'reversal_fields']:
        retained = DATA / label
        if path == retained or retained in path.parents:
            raise ValueError('Choose a fresh campaign output under data/runs; retained evidence is immutable.')

def make_model(key):
    return {'wave_pinning':wave_pinning,'otsuji':otsuji,'goryachev':goryachev,'legi':legi}[key].make_model()

def _sim(model, state, length, cells, stimulus, end, interval=2.5, max_step=10.0):
    state = np.asarray(state, dtype=float)
    if model.nonnegative:
        minimum = float(np.min(state))
        if minimum < -1e-8:
            raise ValueError(f"restart state violates nonnegative invariant by {minimum:g}")
        state = np.maximum(state, 0.0)
    return simulate(model, initial_state=state, length=float(length), cells=int(cells),
                    protocol=stimulus, end_time=float(end), output_interval=float(interval),
                    max_step=float(max_step), backend="bdf")

def _active(tr):
    if "active" in tr.fields:
        return tr.fields["active"]
    return tr.observable

def _orientation(tr):
    u = np.maximum(_active(tr), 0.0)
    center = 2.0 * tr.x / (tr.x[0] + tr.x[-1]) - 1.0
    denom = u.sum(axis=1)
    return np.divide(u @ center, denom, out=np.full(len(denom), np.nan), where=denom > 1e-12)

def model_for(key):
    if key=='holmes_model3':return replace(ho.HolmesModel3Eq56(),rac_activation=ho.ACTIVATION)
    if key=='debelly':return ds.DeBelly(ds.DeBellyParameters(**db._parameters().__dict__))
    if key=='spring':return sp
    return rd.make_model(key)

def initial(key,L,n=201):
    m=model_for(key)
    if key=='debelly':return None
    if key=='spring':return m.initial(n,float(L),perturb=1e-3)
    if key=='holmes_model3':return ho._low_state(m,L,n)
    grid=rd.Grid1D(L,n)
    return m.initial_state(grid)

def run(key,L,state,cue,end,interval=5,n=201,max_step=10):
    m=model_for(key)
    if key=='debelly':
        tr=ds.simulate_debelly(m,protocol=cue,length=L,cells=n,end_time=end,output_interval=interval,step=.5,initial_checkpoint=state)
        fields=tr.fields;active=fields['rac'];x=tr.x;t=tr.time
        final={name:fields[name][-1].copy() for name in DB_NAMES}
        ori=db._orientation(active)
    elif key=='spring':
        tr=sp._run(L,state,cue,end,sample=interval,cells=n)
        active=tr['active'];x=tr['x_fraction']*L;t=tr['time'];ori=tr['orientation']
        fields={name:tr[name] for name in ['active','inactive','length']};final=sp._end(tr)
    else:
        sim=ho._sim if key=='holmes_model3' else rd._sim
        tr=sim(m,state,L,n,cue,end,interval=interval,max_step=max_step)
        fields={name:tr.fields[name] for name in m.state_names};active=rd._active(tr)
        x=tr.x;t=tr.time;ori=ho._orientation(tr) if key=='holmes_model3' else rd._orientation(tr)
        final=np.vstack([fields[name][-1] for name in m.state_names])
    active=np.maximum(active,0)
    split=max(2,n//3);base=np.min(active,axis=1)
    return SimpleNamespace(time=np.asarray(t),fields=fields,active=active,x=np.asarray(x),state=final,
        orientation=np.asarray(ori),contrast=np.ptp(active,axis=1),
        left=np.max(active[:,:split],axis=1)-base,right=np.max(active[:,-split:],axis=1)-base)

def concatenate(a,b,offset):
    if a is None:
        b.time=b.time+offset;return b
    a.time=np.concatenate([a.time,offset+b.time[1:]])
    for name in ['active','orientation','contrast','left','right']:
        setattr(a,name,np.concatenate([getattr(a,name),getattr(b,name)[1:]],axis=0))
    a.fields={name:np.concatenate([a.fields[name],b.fields[name][1:]],axis=0) for name in a.fields}
    a.state=b.state;return a

def save_checkpoint(path,key,state,tr):
    contents=dict(time=tr.time,x=tr.x,active=tr.active,orientation=tr.orientation,contrast=tr.contrast,left=tr.left,right=tr.right)
    contents.update({'field_'+name:y for name,y in tr.fields.items()})
    if key=='debelly':contents.update({'state_'+name:value for name,value in state.items()})
    else:contents['state']=state
    np.savez_compressed(path,**contents)

def read_checkpoint(path,key):
    with np.load(path) as z:
        state={name:z['state_'+name].copy() for name in DB_NAMES} if key=='debelly' else z['state'].copy()
        tr=SimpleNamespace(**{name:z[name].copy() for name in ['time','x','active','orientation','contrast','left','right']},
                           fields={name[6:]:z[name].copy() for name in z.files if name.startswith('field_')},state=state)
    return state,tr

def sustained(time,good,dwell,start=0):
    for i in np.flatnonzero((time>=start)&good):
        j=np.searchsorted(time,time[i]+dwell,side='left')
        if j<len(time) and np.all(good[i:j+1]):return float(time[i]-start)
    return None

def save_trace(path,tr):
    fields={name:tr.fields[name] for name in ['length','boundary_separation','mca_bound','rho','tension'] if name in tr.fields}
    np.savez_compressed(path,time=tr.time,x=tr.x,active=tr.active,orientation=tr.orientation,
                        contrast=tr.contrast,left=tr.left,right=tr.right,**fields)

def now():return datetime.now(timezone.utc).isoformat()

def convert(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    raise TypeError(type(v).__name__)

def write(path,value):
    path=Path(path)
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w',prefix=path.name+'.',suffix='.tmp',dir=path.parent,delete=False) as handle:
        json.dump(value,handle,indent=2,default=convert,allow_nan=False)
        handle.write('\n')
        temporary=Path(handle.name)
    temporary.replace(path)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def csvwrite(path,rows):
    with path.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

rd=SimpleNamespace(make_model=make_model,Grid1D=Grid1D,_sim=_sim,_active=_active,_orientation=_orientation,RD_MODELS=RD_MODELS)
ds=db

def tasks():
    return PROTOCOL['tasks']

def sources():
    paths=list(CODE.glob('*.py'))+[CODE/'protocol.json',CODE/'environment.json']
    paths+=list((CODE/'models').glob('*.py'))+list((CODE/'src/polarity1d').glob('*.py'))
    paths.append(DATA/'inputs/fixed_scales.json')
    return {str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}

class ConstantInput:
    event_times=()
    def __init__(self,value=0.):self.value=float(value)
    def __call__(self,time,x,length):return np.full_like(np.asarray(x,dtype=float),self.value)
