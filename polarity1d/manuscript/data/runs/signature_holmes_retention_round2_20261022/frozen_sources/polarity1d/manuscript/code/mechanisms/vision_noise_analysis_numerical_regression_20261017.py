"""Direct-product regression for the separately versioned numerical amendment."""
from pathlib import Path
import importlib.util
import json
import hashlib
import shutil
import sys
from datetime import datetime,timezone
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/"data/runs/vision_noise_main_20261017"
SOURCE=Path(__file__).with_name("vision_noise_analysis_stable_20261017.py")
OUT=RUN/"results/analysis_numerical_amendment"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def check(a,vector,maxlag=1600):
    parts=a.correlation_parts(vector,maxlag)
    corr=a.ratio_correlation(*parts[:3])
    error_c=error_num=error_e0=error_e1=0.;undefined=0
    for lag in range(maxlag+1):
        earlier=vector[:len(vector)-lag];later=vector[lag:]
        numerator=float(np.vdot(earlier,later).real/len(earlier))
        e0=float(np.mean(abs(earlier)**2));e1=float(np.mean(abs(later)**2))
        denom=np.sqrt(e0)*np.sqrt(e1)
        if denom==0:
            assert not np.isfinite(corr[lag]);undefined+=1;continue
        direct=numerator/denom
        assert np.isfinite(corr[lag])
        error_c=max(error_c,abs(corr[lag]-direct))
        error_num=max(error_num,abs(parts[0][lag]-numerator)/denom)
        error_e0=max(error_e0,abs(parts[1][lag]-e0)/e0)
        error_e1=max(error_e1,abs(parts[2][lag]-e1)/e1)
    assert error_c<1e-8,(error_c,error_num,error_e0,error_e1)
    assert error_e0<3e-12 and error_e1<3e-12
    return dict(lags=maxlag+1,max_absolute_correlation_error=error_c,
                max_normalized_numerator_error=error_num,
                max_relative_origin_energy_error=error_e0,max_relative_later_energy_error=error_e1,
                direct_fallback_lags=int(parts[4]["direct_fallback"].sum()),undefined_energy_lags=undefined,
                no_correlation_clipping=True)


def main():
    a=module(SOURCE,"amended_noise_analysis");manifest=json.loads((RUN/"manifest.json").read_text())
    rows=[]
    for job in manifest["jobs"]:
        if job["kind"]=="noise":continue
        path=RUN/"workers"/f"{job['slot']:04d}"/"trace.npz";status=json.loads(path.with_name("status.json").read_text())
        assert status["complete"] and status["producer"]=="terminal"
        with np.load(path) as z:t=z["time"].copy();P=z["P_real"]+1j*z["P_imag"]
        for window,(start,end) in {"full":(0,32000),"first_half":(0,16000),"last_half":(16000,32000)}.items():
            result=check(a,P[(t>=start)&(t<=end)])
            rows.append(dict(slot=job["slot"],model=job["model"],circumference=job["circumference"],kind=job["kind"],fraction=job["mean_fraction"],window=window,source=str(path),source_sha256=sha(path),**result))
    t=np.arange(6401)
    fixtures={"decaying_plus_tiny_oscillation":np.exp(-t/80)+1e-16*np.exp(.17j*t),
              "impulse_plus_tiny_tail":np.r_[1.+0j,1e-120*np.exp(.003j*t[1:])],
              "all_tiny_resolved":1e-120*np.exp(.003j*t),
              "all_large_resolved":1e100*np.exp(.003j*t),
              "genuinely_underflowed_energy":np.r_[1.+0j,1e-170*np.exp(.003j*t[1:])]}
    fixture_rows={name:check(a,vector) for name,vector in fixtures.items()}
    OUT.mkdir(parents=True,exist_ok=True)
    copy=OUT/"corrected_analysis_source.py";shutil.copy2(SOURCE,copy)
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,
                original_frozen_source=str(RUN/"protocol/analysis_source.py"),original_frozen_sha256=sha(RUN/"protocol/analysis_source.py"),
                corrected_source=str(copy),corrected_sha256=sha(copy),regression_source=str(Path(__file__).resolve()),regression_sha256=sha(__file__),
                unchanged_producer_manifest_sha256=sha(RUN/"manifest.json"),unchanged_configuration_sha256=sha(RUN/"configuration.json"),
                statistical_definition_changed=False,simulation_or_noise_law_changed=False,correlation_clipping=False,
                correction="Positive suffix sums avoid subtracting large prefix energies; direct products resolve weak FFT endpoint products; sqrt endpoint energies separately avoids product underflow",
                actual_control_cases=27,actual_control_windows=len(rows),all_prescribed_lags_checked=1601,
                maximum_actual_correlation_error=max(r["max_absolute_correlation_error"] for r in rows),
                controls=rows,extreme_fixtures=fixture_rows,original_selfcheck=a.selfcheck())
    path=OUT/"REGRESSION.json";path.write_text(json.dumps(record,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:record[k] for k in ("passed","actual_control_cases","actual_control_windows","maximum_actual_correlation_error","corrected_sha256")},indent=2))

if __name__=="__main__":main()
