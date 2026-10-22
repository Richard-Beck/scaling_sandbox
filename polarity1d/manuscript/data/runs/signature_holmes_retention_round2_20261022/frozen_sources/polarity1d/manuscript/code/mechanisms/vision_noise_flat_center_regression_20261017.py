"""Validate the saved-center observability adapter against every main archive."""
from pathlib import Path
import importlib.util,inspect,json,hashlib,shutil
from datetime import datetime,timezone
import numpy as np

ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/"data/runs/vision_noise_main_20261017"
CODE=Path(__file__).parent;SOURCE=CODE/"vision_noise_analysis_stable_v2_20261017.py"
OUT=RUN/"results/flat_center_validation_amendment"

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    new=module(SOURCE,"v2");old=module(RUN/"results/statistics_workflow/amended_analysis.py","v1")
    unchanged={name:inspect.getsource(getattr(new,name))==inspect.getsource(getattr(old,name)) for name in
               ("analyze_group","site_displacement_group","select_histories","correlation_parts","ratio_correlation","direction_parts")}
    assert all(unchanged.values())
    manifest=json.loads((RUN/"manifest.json").read_text());rows=[]
    for job in manifest["jobs"]:
        folder=RUN/"workers"/f"{job['slot']:04d}"
        row=dict(job,archive=str(folder/"trace.npz"),profile_archive=str(folder/"primary.npz"),fraction=job["mean_fraction"],case_id=f"slot_{job['slot']:04d}")
        case=new.load_case(row,RUN)
        with np.load(folder/"trace.npz") as z:
            flat=z["current_half_contrast_regions"]==0
            contrast=z["contrast"].copy();R=np.hypot(z["P_real"],z["P_imag"])
        assert np.all(contrast[flat]<=1e-12)
        assert not any(r["single_region_resolved"] for r,f in zip(case["scalar_rows"],flat) if f)
        assert not np.any(R[flat]>=.05)
        record=dict(slot=job["slot"],model=job["model"],kind=job["kind"],circumference=job["circumference"],fraction=job["mean_fraction"],
                    **case["audit"],flat_scalar_frames=int(flat.sum()),flat_scalar_direction_eligible=0,flat_scalar_site_eligible=0)
        rows.append(record)
    OUT.mkdir(parents=True,exist_ok=True);target=OUT/"corrected_analysis_source.py";shutil.copy2(SOURCE,target)
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,validated_cases=len(rows),profile_snapshots=sum(r["profile_frames"] for r in rows),
                original_failed_source_sha256=sha(RUN/"results/statistics_workflow/amended_analysis.py"),corrected_sha256=sha(target),
                corrected_source=str(target),regression_source=str(Path(__file__).resolve()),regression_sha256=sha(__file__),
                unchanged_readout_function_sources=unchanged,statistical_or_scalar_eligibility_change=False,tolerance_loosened=False,
                producer_manifest_sha256=sha(RUN/"manifest.json"),original_frozen_analysis_sha256=sha(RUN/"protocol/analysis_source.py"),
                cases_with_unobserved_saved_centers=[r["slot"] for r in rows if r["unobserved_flat_held_center_frames"]],
                unobserved_saved_center_frames=sum(r["unobserved_flat_held_center_frames"] for r in rows),
                maximum_identifiable_center_error=max((r["maximum_observable_center_error"] or 0) for r in rows),cases=rows)
    (OUT/"REGRESSION.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:result[k] for k in ("passed","validated_cases","profile_snapshots","cases_with_unobserved_saved_centers","unobserved_saved_center_frames","maximum_identifiable_center_error","corrected_sha256")},indent=2))

if __name__=="__main__":main()
