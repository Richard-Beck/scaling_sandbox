"""Cheap native solver check; writes no campaign outputs and submits no jobs."""
import json
import numpy as np
import native as r

rows=[]
for model in r.MODELS:
    length=1.0 if model=='spring' else 15.0
    cells=31
    state=r.initial(model,length,cells)
    trace=r.run(model,length,state,r.ConstantInput(.1 if model=='legi' else 0.),1.0,interval=.5,n=cells)
    assert len(trace.time)>=3
    assert np.isfinite(trace.active).all()
    assert np.isfinite(trace.state).all() if not isinstance(trace.state,dict) else all(np.isfinite(v).all() for v in trace.state.values())
    rows.append(dict(model=model,cells=cells,end_time=float(trace.time[-1]),active_shape=list(trace.active.shape),finite=True))
print(json.dumps(dict(passed=True,native_solvers=rows,imports_root=str(r.CODE)),indent=2))
