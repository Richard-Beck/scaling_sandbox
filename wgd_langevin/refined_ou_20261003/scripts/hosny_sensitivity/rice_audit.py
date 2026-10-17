"""Check stable Rice means against a fresh unpatched SciPy distribution."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json
import numpy as np
from scipy.stats._continuous_distns import rice_gen
import analyze as analysis

def main():
    x=np.array([0.,.1,1.,10.])
    original=rice_gen(a=0,name='independent_rice')
    expected=original.mean(x);stable=analysis.rice_mean_stable(x)
    error=float(np.max(np.abs(expected-stable)))
    if error>=1e-8:raise AssertionError(f'Stable Rice mean differs from true SciPy: {error}')
    meta=dict(generated_utc=datetime.now(timezone.utc).isoformat(),test_noncentralities=x.tolist(),true_scipy_means=expected.tolist(),scaled_bessel_means=stable.tolist(),maximum_absolute_difference=error,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),note='Fresh scipy rice_gen instance avoids global distribution instance monkey patches in the refined fit module. Inherited analysis validator compares two stable implementations; this independent check supplies genuine agreement with unpatched SciPy.')
    (analysis.OUT/'rice_independent_validation.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':main()
