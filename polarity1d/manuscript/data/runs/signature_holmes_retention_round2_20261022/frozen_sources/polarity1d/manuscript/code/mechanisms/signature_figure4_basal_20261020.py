"""One fresh C0=1.5 Spring basal using the unchanged native preparation rule.

Freezing requires the completed E1/E2 native gate. The native worker is the
existing preparation driver, copied without edits; this module independently
qualifies its saved coupled-state history before E3 can consume the checkpoint.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shlex
import shutil

import numpy as np
from scipy.signal import find_peaks

import signature_figure4_campaign_20261020 as C


def prepare(out,gate):
    out=Path(out).resolve();gate=Path(gate).resolve();release=C.read(gate)
    if out.exists() or not out.is_relative_to(C.MS/"data/runs"):raise ValueError("Choose a fresh native basal run")
    if release.get("e1_e2_complete") is not True or not release.get("native_stage_gates"):
        raise ValueError("All E1/E2 native conditions must qualify before E3 basal preparation")
    for path,digest in release["native_stage_gates"].items():
        if C.sha(path)!=digest:raise ValueError("E1/E2 native completion gate changed")
    parent=C.read(C.SOURCE/"manifest.json")
    for relative,digest in parent["source_sha256"].items():
        if C.sha(C.SOURCE/relative)!=digest:raise ValueError("Retained native E0 source changed")
    preparation_source=C.MS/"data/runs/revision_ring_spring_20261019"
    preparation=C.read(preparation_source/"manifest.json");source_root=preparation_source/"frozen_sources/code"
    for relative,digest in preparation["source_sha256"].items():
        if C.sha(preparation_source/"frozen_sources"/relative)!=digest:raise ValueError("Retained native basal preparation source changed")
    paths=[source_root/"mechanisms"/name for name in ("revision_spring_ring_basals_20261019.py","spring_ring_model_20261019.py","spring_ring_closure_20261017.py")]
    paths.extend(source_root/"core/models"/name for name in ("spring.py","__init__.py"))
    for path in paths[1:]:
        platform="chemical" if path.name=="__init__.py" else "spring"
        relative=str(Path("frozen_sources")/platform/path.relative_to(preparation_source/"frozen_sources"))
        if parent["source_sha256"].get(relative)!=C.sha(path):raise ValueError("Spring native source differs from retained E0: "+str(path))
    for name in ("frozen_sources","preparations","logs","results","inputs"):(out/name).mkdir(parents=True)
    hashes={}
    for path in paths:
        relative=path.relative_to(preparation_source/"frozen_sources");target=out/"frozen_sources"/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);hashes[str(relative)]=C.sha(target)
    retained=C.MS/"data/runs/spring_ring_preparation_pilot_20261017/workers/qualified_basal.npz"
    shutil.copy2(retained,out/"inputs/retained_C2.npz")
    shutil.copy2(gate,out/"inputs/e1_e2_native_gate.json")
    qualifier=out/"frozen_sources"/Path(__file__).name;shutil.copy2(__file__,qualifier)
    helper=qualifier.with_name(Path(C.__file__).name);shutil.copy2(C.__file__,helper)
    manifest=dict(generated_at=C.now(),model="Spring mean-force global-extension circular analogue",jobs=[dict(slot=0,resting_circumference=1.5,cells=150,reused=False)],
        initial="Fresh native uniform concentrations; no remapped or interpolated basal",preparation_cue_duration=10.,basal_release=13000.,cue_peak=4.,cue_full_physical_support=.2,
        source_sha256=hashes,retained_input_sha256=C.sha(out/"inputs/retained_C2.npz"),retained_original=str(retained),
        retained_original_audit=str(retained.parents[1]/"results/audit.json"),
        retained_C2_input_is_unused_in_fresh_worker=True,solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.),
        qualification=dict(tail_window=1000.,scaled_range=1e-4,rhs_projected_to_1000=1e-4,minimum_contrast=.05,one_half_contrast_region=True,one_resolved_peak=True,physical_amount_relative_tolerance=1e-7),
        e1_e2_completion_gate=dict(path=str(out/"inputs/e1_e2_native_gate.json"),sha256=C.sha(out/"inputs/e1_e2_native_gate.json")),
        original_preparation_manifest=dict(path=str(preparation_source/"manifest.json"),sha256=C.sha(preparation_source/"manifest.json")),
        independent_qualifier_sources={str(qualifier):C.sha(qualifier),str(helper):C.sha(helper)},automatic_submission=False)
    C.write(out/"manifest.json",manifest);C.write(out/"freeze.json",dict(manifest_sha256=C.sha(out/"manifest.json")))
    worker=out/"worker.sh";driver=out/"frozen_sources/code/mechanisms/revision_spring_ring_basals_20261019.py"
    worker.write_text('#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec '+shlex.join([C.PYTHON,'-B',str(driver),'worker','--run',str(out),'--slot','0'])+'\n')
    return dict(run=str(out),manifest_sha256=C.sha(out/"manifest.json"),integrations=0,submitted=False)


def qualify(out):
    out=Path(out).resolve();manifest=C.read(out/"manifest.json")
    if C.sha(out/"manifest.json")!=C.read(out/"freeze.json")["manifest_sha256"]:raise ValueError("Frozen native basal manifest changed")
    for relative,digest in manifest["source_sha256"].items():
        if C.sha(out/"frozen_sources"/relative)!=digest:raise ValueError("Frozen native basal source changed")
    for path,digest in manifest["independent_qualifier_sources"].items():
        if C.sha(path)!=digest:raise ValueError("Frozen basal qualifier changed")
    status=C.read(out/"preparations/C1.5/status.json");folder=out/"preparations/C1.5"
    if status.get("status")!="qualified" or status.get("complete") is not True:raise ValueError("Native C1.5 preparation did not qualify; no retuning permitted")
    for path,digest in ((status["history_path"],status["history_sha256"]),(status["checkpoint"],status["checkpoint_sha256"])):
        if C.sha(path)!=digest:raise ValueError("Native basal scientific output changed")
    for name in list(C.sys.modules):
        if name in ("spring_ring_model_20261019","spring_ring_closure_20261017") or name=="models" or name.startswith("models."):
            C.sys.modules.pop(name)
    C.sys.path.insert(0,str(out/"frozen_sources/code/mechanisms"))
    spec=C.importlib.util.spec_from_file_location("qualified_f4_intermediate_spring",out/"frozen_sources/code/mechanisms/spring_ring_model_20261019.py")
    native=C.importlib.util.module_from_spec(spec);C.sys.modules[spec.name]=native;spec.loader.exec_module(native);model=native.Closure(1.5,150)
    with np.load(status["history_path"]) as data:t=data["time"].copy();states=data["states"].copy();positions=data["material_positions"].copy()
    with np.load(status["checkpoint"]) as data:final=data["state"].copy()
    np.testing.assert_array_equal(positions,model.material_positions);np.testing.assert_array_equal(states[0],model.initial());np.testing.assert_array_equal(states[-1],final)
    if states.shape!=(len(t),301) or t[0]!=0 or t[-1]!=13010. or not np.all(np.diff(t)>0):raise ValueError("Native fresh preparation history is incomplete")
    if not np.isfinite(states).all() or states[:,:-1].min() < -5e-10 or np.any(states[:,-1]<=0):raise ValueError("Native preparation positivity failed")
    scales=np.r_[np.full(300,2.),1.5];active=final[:150];contrast=float(np.ptp(active));mask=active>=active.min()+contrast*.5
    regions=int(np.sum(mask&~np.roll(mask,1))) if contrast>1e-12 else 0
    peaks,_=find_peaks(np.tile(active,3),prominence=.05);peaks=peaks[(peaks>=150)&(peaks<300)]-150
    pool=states[:,-1]*np.mean(states[:,:150]+states[:,150:-1],axis=1);drift=float(np.max(abs(pool/pool[0]-1)))
    tail_range=float(np.max(np.ptp(states[t>=t[-1]-1000],axis=0)/scales));rhs=float(np.max(abs(model.rhs(final,0.))/scales)*1000.)
    passed=contrast>=.05 and regions==1 and len(peaks)==1 and drift<=1e-7 and tail_range<=1e-4 and rhs<=1e-4
    report=dict(generated_utc=C.now(),passed=passed,native_manifest_sha256=C.sha(out/"manifest.json"),native_status_sha256=C.sha(folder/"status.json"),
        exact_native_uniform_initial_state=True,all_coupled_fields_qualified=True,contrast=contrast,half_contrast_regions=regions,resolved_peaks=len(peaks),
        full_physical_amount_relative_drift=drift,tail_1000_scaled_range=tail_range,rhs_projected_1000=rhs,
        checkpoint=status["checkpoint"],checkpoint_sha256=status["checkpoint_sha256"],native_kinetics_retuned=False)
    C.write(out/"results/independent_qualification_audit.json",report)
    if not passed:raise ValueError("Independent fresh C1.5 basal qualification failed")
    record=dict(model="spring_ring",size=1.5,cells=150,source=status["checkpoint"],input_sha256=status["checkpoint_sha256"],
        state_sha256=__import__('hashlib').sha256(np.asarray(final,dtype='<f8',order='C').tobytes()).hexdigest(),basal_circumference=float(final[-1]),qualified=True,
        qualification_status=str(folder/"status.json"),qualification_status_sha256=C.sha(folder/"status.json"),
        independent_qualification_audit=str(out/"results/independent_qualification_audit.json"),independent_qualification_audit_sha256=C.sha(out/"results/independent_qualification_audit.json"),
        native_kinetics_sources_identical_to_retained_E0=True)
    C.write(out/"results/qualified_basal_record.json",record);return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest="action",required=True)
    prep=sub.add_parser("prepare");prep.add_argument("--out",type=Path,required=True);prep.add_argument("--gate",type=Path,required=True)
    work=sub.add_parser("qualify");work.add_argument("--out",type=Path,required=True)
    args=parser.parse_args();result=prepare(args.out,args.gate) if args.action=="prepare" else qualify(args.out)
    print(__import__('json').dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
