"""Full retained-grid first-angular-moment extraction; no native integration.

Use the reviewed global_polarity_vector_20261017 moment definitions unchanged.
All original histories and scientific sources remain immutable. Plot/report
assembly is a later step after the common-active recovery and cue checks arrive.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name]="1"
sys.dont_write_bytecode=True
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
SOURCE_RUN=ROOT/"data/runs/chemical_ring_functional_v2_20261017"
DEFAULT=ROOT/"data/runs/vision_ring_moments_20261017"
DEFINITION=Path(__file__).with_name("global_polarity_vector_20261017.py")
EPOCHS=((1000,0,"cue withdrawal"),(4000,3000,"3000 cue-free units"),
        (9000,8000,"8000 cue-free units"),(14000,13000,"13000 cue-free units"))
GATES=(.05,.1,.2)


def now():return datetime.now(timezone.utc).isoformat()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(2**20),b""):h.update(block)
    return h.hexdigest()


def write_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+"\n")


def read_csv(path):
    with Path(path).open(newline="") as stream:return list(csv.DictReader(stream))


def write_csv(path,rows):
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def label(gate):return str(gate).replace(".","p")


def load_phase(path,offset):
    with np.load(path) as z:
        t=z["time"].copy()+offset;x=z["x"].copy()
        primary=z["states"][:,0].copy() if "states" in z else z["primary"].copy()
        legacy=z["harmonic_real"]+1j*z["harmonic_imag"] if "harmonic_real" in z else None
        schema="full native states" if "states" in z else "primary history plus complete native snapshots"
    return t,x,primary,legacy,schema


def main(run):
    began=time.monotonic()
    if (run/"audit.json").exists():raise RuntimeError("Completed extraction exists; preserve it")
    for directory in (run,run/"analysis_sources",run/"results"):
        directory.mkdir(parents=True,exist_ok=True)
    copied=run/"analysis_sources"/DEFINITION.name
    shutil.copy2(DEFINITION,copied)
    shutil.copy2(__file__,run/"analysis_sources"/Path(__file__).name)
    spec=importlib.util.spec_from_file_location("reviewed_ring_moment",copied)
    reviewed=importlib.util.module_from_spec(spec);spec.loader.exec_module(reviewed)
    sources={str(Path(__file__).resolve()):sha(__file__),str(DEFINITION):sha(DEFINITION)}
    for path in (SOURCE_RUN/"results/matrix.csv",SOURCE_RUN/"results/audit.json",SOURCE_RUN/"results/source_sha256.json"):
        sources[str(path)]=sha(path)
    matrix=read_csv(SOURCE_RUN/"results/matrix.csv")
    native_audit=json.loads((SOURCE_RUN/"results/audit.json").read_text())
    expected=json.loads((SOURCE_RUN/"results/source_sha256.json").read_text())
    assert native_audit["passed"] and len(matrix)==867
    assert sum(r["status"]=="complete" for r in matrix)==652
    analytic=reviewed.analytic_checks()
    write_json(run/"manifest.json",dict(created_at=now(),source_matrix=str(SOURCE_RUN/"results/matrix.csv"),
        source_matrix_sha256=sources[str(SOURCE_RUN/"results/matrix.csv")],definition_source=str(DEFINITION),
        definition_sha256=sha(DEFINITION),generator_sha256=sha(__file__),planned_cells=867,completed_native_conditions=652,
        readout_onset_times=[e[0] for e in EPOCHS],post_withdrawal_times=[e[1] for e in EPOCHS],
        primary_angle_gate=.1,sensitivity_angle_gates=[.05,.2],new_integrations=0,
        definitions="P=M1/M1(0); R=abs(P); q1=abs(M1)/M0; q2=abs(M2)/M0; M0=integral primary activity. Reviewed moment() unchanged.",
        scientific_limits="Activity moment, not cell motion, unique site identity, migration direction or deterministic persistence fit."))
    summaries=[];readouts=[];source_rows=[];trace_cases=[];offsets=[0]
    traces={key:[] for key in ("time","M0","M1_real","M1_imag","M2_real","M2_imag","P_x","P_y","R","q1","q2")}
    maximum_bg_relative_error=0.;maximum_legacy_error=0.;total_frames=0
    for matrix_index,row in enumerate(matrix):
        rec={key:row.get(key,"") for key in ("model","circumference","fraction","realized_fraction","short_route","long_route",
            "status","checkpoint_sha256","final_location_association","final_peak_count","final_single_region_spans_old_and_cue")}
        rec["matrix_row_index"]=matrix_index
        if row["status"]!="complete":summaries.append(rec);continue
        C=float(row["circumference"])
        phases=[]
        for phase,key,offset in (("cue","cue_path",0),("recovery","recovery_path",1000)):
            path=Path(row[key]);measured=sha(path);assert measured==expected[str(path)],path
            sources[str(path)]=measured
            phase_data=load_phase(path,offset);phases.append(phase_data)
            source_rows.append(dict(model=row["model"],circumference=C,fraction=float(row["fraction"]),phase=phase,
                path=str(path),sha256=measured,logical_bytes=path.stat().st_size,allocated_bytes=path.stat().st_blocks*512,
                source_schema=phase_data[4],samples=len(phase_data[0]),cells=len(phase_data[1]),phase_offset=offset))
        tc,x,pc,lc,_=phases[0];tr,xr,pr,lr,_=phases[1]
        np.testing.assert_array_equal(x,xr);np.testing.assert_array_equal(pc[-1],pr[0])
        t=np.r_[tc,tr[1:]];primary=np.concatenate((pc,pr[1:]),axis=0)
        assert t[0]==0 and t[-1]==14000 and np.all(np.diff(t)>0)
        M0=primary.sum(axis=1)*(C/len(x));assert np.all(M0>0)
        M1=reviewed.moment(primary,x,C);M2=reviewed.moment(primary,x,C,order=2)
        assert abs(M1[0])/M0[0]>1e-8
        P=M1/M1[0];R=abs(P);angles=np.angle(P)*180/np.pi
        q1=abs(M1)/M0;q2=abs(M2)/M0
        # Evaluate the defining cancellation directly on every native frame.
        bg=primary.min(axis=1)
        bg_M1=reviewed.moment(primary-bg[:,None],x,C)
        bg_error=float(np.max(abs(bg_M1-M1))/max(float(M0.max()),1.))
        assert bg_error<1e-12
        maximum_bg_relative_error=max(maximum_bg_relative_error,bg_error)
        legacy_error=""
        if lc is not None and lr is not None:
            legacy=np.r_[lc,lr[1:]];legacy_error=float(np.max(abs(abs(legacy)-q1)))
            assert legacy_error<1e-12
            maximum_legacy_error=max(maximum_legacy_error,legacy_error)
        index=len(trace_cases);trace_cases.append(dict(trace_index=index,matrix_row_index=matrix_index,
            model=row["model"],circumference=C,fraction=float(row["fraction"]),offset_start=offsets[-1],frames=len(t)))
        rec.update(trace_index=index,frames=len(t),max_recovery_sample_gap=float(np.diff(tr).max()),
            basal_M0=float(M0[0]),basal_M1_real=float(M1[0].real),basal_M1_imag=float(M1[0].imag),
            basal_first_moment_abs=float(abs(M1[0])),basal_q1=float(q1[0]),basal_q2=float(q2[0]),
            sampled_minimum_R=float(R.min()),sampled_time_at_minimum_R=float(t[np.argmin(R)]),
            sampled_minimum_R_cue_on=float(R[t<=1000].min()),sampled_minimum_R_cue_free=float(R[t>=1000].min()),
            minimum_primary=float(primary.min()),maximum_primary=float(primary.max()),
            uniform_background_cancellation_relative_error=bg_error,legacy_q1_max_error=legacy_error)
        for gate in GATES:
            eligible=R>=gate;tag=label(gate)
            rec[f"sampled_maximum_unmasked_abs_angle_deg_gate_{tag}"]=float(abs(angles[eligible]).max())
            rec[f"masked_frames_gate_{tag}"]=int((~eligible).sum())
        for onset,postoff,description in EPOCHS:
            ii=np.flatnonzero(np.isclose(t,onset,rtol=0,atol=1e-10));assert len(ii)==1
            i=int(ii[0]);out={key:rec[key] for key in ("model","circumference","fraction","realized_fraction","short_route","long_route")}
            out.update(time_from_onset=onset,time_after_withdrawal=postoff,readout_description=description,
                M0=float(M0[i]),M1_real=float(M1[i].real),M1_imag=float(M1[i].imag),M2_real=float(M2[i].real),M2_imag=float(M2[i].imag),
                P_x=float(P[i].real),P_y=float(P[i].imag),R=float(R[i]),q1=float(q1[i]),q2=float(q2[i]))
            for gate in GATES:
                tag=label(gate);angle=float(angles[i]) if R[i]>=gate else ""
                out[f"angle_deg_gate_{tag}"]=angle;out[f"angle_available_gate_{tag}"]=bool(R[i]>=gate)
                rec[f"angle_deg_t{onset}_gate_{tag}"]=angle
            for key in ("P_x","P_y","R","q1","q2","M0"):
                rec[f"{key}_t{onset}"]=out[key]
            readouts.append(out)
        for key,array in (("time",t),("M0",M0),("M1_real",M1.real),("M1_imag",M1.imag),("M2_real",M2.real),
                          ("M2_imag",M2.imag),("P_x",P.real),("P_y",P.imag),("R",R),("q1",q1),("q2",q2)):
            traces[key].append(np.asarray(array,dtype=np.float64))
        offsets.append(offsets[-1]+len(t));total_frames+=len(t);summaries.append(rec)
        if index%25==0:
            write_json(run/"status.json",dict(generated_at=now(),producer="running retained-data extraction",complete=False,
                       conditions_complete=index+1,total_conditions=652,elapsed_seconds=time.monotonic()-began))
    assert len(summaries)==867 and len(trace_cases)==652 and len(readouts)==2608
    traces={key:np.concatenate(value) for key,value in traces.items()}
    traces["offsets"]=np.asarray(offsets,dtype=np.int64)
    traces["matrix_row_index"]=np.asarray([r["matrix_row_index"] for r in trace_cases],dtype=np.int64)
    np.savez_compressed(run/"results/moment_traces.npz",**traces)
    for name,rows in (("matrix_summary",summaries),("readout_sensitivity",readouts),("trace_index",trace_cases),("sources",source_rows)):
        write_csv(run/"results"/(name+".csv"),rows)
    write_json(run/"source_sha256.json",sources)
    outputs={str(p.relative_to(run)):sha(p) for p in (run/"results").iterdir() if p.is_file()}
    allocated=sum(p.stat().st_blocks*512 for p in run.rglob("*") if p.is_file())
    audit=dict(generated_at=now(),passed=True,producer="completed retained-data extraction",new_integrations=0,
        planned_matrix_cells=867,complete_conditions=652,unavailable_cells=215,readout_rows=2608,sampled_frames=total_frames,
        exact_primary_phase_continuity=True,all_requested_readout_times_exactly_available=True,
        reviewed_definition_sha256=sha(copied),reviewed_analytic_checks=analytic,
        maximum_uniform_background_cancellation_error_relative_to_active_amount=maximum_bg_relative_error,
        maximum_legacy_q1_error=maximum_legacy_error,source_phase_files=len(source_rows),unique_source_phase_files=len(set(r["path"] for r in source_rows)),
        normalized_angle_primary_gate=.1,normalized_angle_sensitivity_gates=[.05,.2],
        scientific_limits="Sampled deterministic activity-vector summaries; no migration direction, unique protrusion identity or persistence fit.",
        elapsed_seconds=time.monotonic()-began,allocated_bytes_before_audit=allocated,outputs_sha256=outputs)
    write_json(run/"audit.json",audit)
    write_json(run/"status.json",dict(generated_at=now(),producer="finished",complete=True,audit=audit))
    print(json.dumps(audit,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--run",type=Path,default=DEFAULT)
    main(parser.parse_args().run.resolve())
