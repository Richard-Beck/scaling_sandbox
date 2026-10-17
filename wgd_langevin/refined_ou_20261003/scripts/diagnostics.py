"""Independent local data sensitivity; explicitly separate from regularizers."""
import json
from datetime import datetime,timezone
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from fit import *

def diagnostic(z):
    z=np.asarray(z)
    def data(x):
        r,rows=data_residuals(analytic(x));return r[[o['study']!='Hosny_baseline' or o.get('quantile')==1 for o in rows]]
    output=[]
    for step in (1e-4,1e-5,1e-6):
        jac=approx_derivative(data,z,method='3-point',rel_step=step)
        sv=np.linalg.svd(jac,compute_uv=False)
        output.append(dict(relative_step=step,singular_values=sv.tolist(),rank_relative_1e6=int(np.sum(sv>sv[0]*1e-6))))
    return dict(data_jacobian=output,regularized_dimension=22,note='Data-only numerical sensitivity, not confidence limits. Each Jemaà study fits only two mean speeds: one scalar per group cannot identify both persistence and noise without regularization.')

def main():
    selected=json.loads((OUT/'selected_fits.json').read_text())
    obj={s:diagnostic(b['z']) for s,b in selected.items()}
    atomic(OUT/'identifiability.json',dict(generated=datetime.now(timezone.utc).isoformat(),scenarios=obj))
    path=OUT/'RESULTS.md';doc=path.read_text();doc+='\nData-only Jacobian ranks (relative singular-value threshold 1e-6; central finite differences at three relative steps):\n\n|Scenario|step 1e-4|step 1e-5|step 1e-6|\n|---|---:|---:|---:|\n'
    for s,d in obj.items():doc+='|'+s+'|'+ '|'.join(str(x['rank_relative_1e6'])+'/22' for x in d['data_jacobian'])+'|\n'
    doc+='\nThe remaining rank is supplied by regularization. Each Jemaà assay observes only speed, leaving a persistence/noise tradeoff per group. Fitted times are therefore conditional on prior and pooling assumptions.\n\nThe Xuan reference curve shows systematic early overprediction and late underprediction (mean log ratios +0.267 and -0.257): the fitted OU curve shape does not reproduce the full digitized curve within its raster envelope. This discrepancy is larger than digitization precision and survives the completed optimization; completion is not evidence of model adequacy.\n'
    path.write_text(doc)
if __name__=='__main__':main()
