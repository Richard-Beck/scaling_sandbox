"""Immutable dense-distance stationary cues on a narrow native size range.

Uses the accepted native worker without changing kinetics, solver or readouts.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
import numpy as np
import pandas as pd
import revision_initial_arc_sizes_20261019 as arc

sizes=arc.prior
base=sizes.base
ROOT=arc.ROOT
DEFAULT=ROOT/"data/runs/revision_arc_narrow_grid_20261019"
PREP=ROOT/"data/runs/revision_arc_narrow_basals_20261019"
OUT=arc.OUT.parent/"narrow_arc_experiment"
CHEMICAL_D=(.5,1.25,2.5,5.,7.5,9.,10.,11.,12.,13.,14.)
SPRING_D=(.06508769377179444,.13017538754358887,.26035077508717774,.5207015501743555,.7,.85,.95,1.05,1.15,1.25)
sha,read,write,now=arc.sha,arc.read,arc.write,arc.now


def core_at_level(q,u,origin,level):
    x=(q-origin+.5)%1-.5;order=np.argsort(x);xx,uu=x[order],u[order]
    threshold=uu.min()+level*np.ptp(uu);inside=np.flatnonzero(uu>=threshold);a,b=inside[0],inside[-1]
    assert np.all(np.diff(inside)==1) and a>0 and b<len(xx)-1
    left=xx[a-1]+(threshold-uu[a-1])/(uu[a]-uu[a-1])*(xx[a]-xx[a-1])
    right=xx[b]+(threshold-uu[b])/(uu[b+1]-uu[b])*(xx[b+1]-xx[b])
    return float(left),float(right)


def geometry(left,right,d,w,C):
    amount=sum(max(0.,min(right,d+w/2+j*C)-max(left,d-w/2+j*C)) for j in (-1,0,1))
    fraction=float(np.clip(amount/w,0,1))
    return dict(left=left,right=right,width=right-left,edge_distance=d-right,
                centre="inside" if left<=d<=right else "outside",
                footprint="inside" if fraction>=1-1e-10 else "outside" if fraction<=1e-10 else "straddling",overlap=fraction,intersection=amount)


def candidates():
    result=[];sources={}
    for path in (arc.OUT/"supporting/data/stationary_arclength_summary.csv",sizes.OUT/"supporting/data/stationary_sizes_summary.csv"):
        sources[str(path)]=sha(path)
        for row in pd.read_csv(path).to_dict("records"):
            if row["status"]=="complete":
                row.update(candidate_table=str(path),candidate_condition_id=int(row["condition_id"]));result.append(row)
    path=ROOT/"data/runs/chemical_ring_functional_v2_20261017/results/matrix.csv";sources[str(path)]=sha(path)
    for row in pd.read_csv(path).to_dict("records"):
        if row["status"]=="complete" and row["circumference"] in (30.,35.,40.):
            result.append(dict(model=row["model"],size_parameter=row["circumference"],fraction=row["fraction"],cue_duration=1000.,
                               actual_cue_path=row["cue_path"],actual_cue_sha256=sha(row["cue_path"]),
                               actual_recovery_path=row["recovery_path"],actual_recovery_sha256=sha(row["recovery_path"]),candidate_table=str(path),candidate_condition_id=None))
    return result,sources


def basal_records():
    accepted=read(arc.ANGLE_RUN/"manifest.json");prepplan=read(PREP/"manifest.json");status=read(PREP/"results/preparation_summary.json")
    assert status["complete"] and status["fresh_qualified"]==2
    records=[]
    for model in base.MODELS:
        for size in ((2.5,2.75,3.) if model=="spring" else (30.,35.,40.)):
            if size in ((3.,) if model=="spring" else (30.,40.)):
                b=next(b for b in accepted["basals"] if b["model"]==model and b["size_parameter"]==size)
                original=arc.ANGLE_RUN/b["input"]
                assert sha(original)==b["input_sha256"]
                qualification=b["qualification"];scales=next(c["fixed_scales"] for c in accepted["cases"] if c["model"]==model and c["size_parameter"]==size)
            elif model=="spring":
                s=next(s for s in status["records"] if s["job"]["size_parameter"]==size);assert s["status"]=="qualified"
                original=Path(s["checkpoint"]);assert sha(original)==s["checkpoint_sha256"]
                qualification=dict(status="qualified",status_path=str(PREP/"workers"/f'{s["job"]["slot"]:03d}'/"status.json"),history_path=s["history_path"],history_sha256=s["history_sha256"])
                qualification["status_sha256"]=sha(qualification["status_path"]);scales=s["fixed_scales"]
            else:
                b=next(b for b in prepplan["exact_existing_chemical_preparations"] if b["model"]==model)
                original=PREP/b["input"];assert sha(original)==b["input_sha256"]
                qualification=dict(status="qualified",status_path=b["qualification_status_path"],status_sha256=b["qualification_status_sha256"],preparation_manifest=b["preparation_manifest"],preparation_manifest_sha256=b["preparation_manifest_sha256"])
                scales=b["fixed_scales"]
            records.append(dict(model=model,size_parameter=size,original=original,qualification=qualification,fixed_scales=scales))
    return records


def freeze(run,proof):
    if run.exists():raise FileExistsError("Fresh immutable dense-distance campaign required")
    assert proof and proof.exists(),"An independent basal qualification record is required"
    qualification=read(proof)
    assert qualification["passed"] and qualification["qualified_basals"]==5
    assert qualification["manifest_sha256"]==sha(PREP/"manifest.json")
    for folder in ("inputs","frozen_sources","workers","logs","results"):(run/folder).mkdir(parents=True)
    previous=read(arc.DEFAULT/"manifest.json");sources={}
    assert list(SPRING_D[:4])==previous["initial_arc_grids"]["spring"]
    for rel,digest in previous["source_sha256"].items():
        p=arc.DEFAULT/"frozen_sources"/rel;assert sha(p)==digest
        target=run/"frozen_sources"/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);sources[rel]=digest
    relative=Path(__file__).relative_to(ROOT);target=run/"frozen_sources"/relative;shutil.copy2(__file__,target);sources[str(relative)]=sha(target)
    for module in (arc,sizes,base):assert sha(module.__file__)==sources[str(Path(module.__file__).relative_to(ROOT))]
    available,table_sources=candidates();cases=[];jobs=[];contexts=[];basals=[];inputs={};rejected=[]
    for record in basal_records():
        model,size=record["model"],record["size_parameter"];spring=model=="spring"
        qualified=[b for b in qualification["records"] if b["model"]==model and b["size_parameter"]==size]
        if qualified:assert sha(record["original"])==qualified[0]["qualified_checkpoint_sha256"]
        path=run/"inputs"/f"{model}_C{size:g}_basal.npz";shutil.copy2(record["original"],path);inputs[str(path.relative_to(run))]=sha(path)
        with np.load(path) as z:state=z["state"].copy()
        n=(len(state)-1)//2 if spring else state.shape[1];C=float(state[-1]) if spring else size
        assert n==int(round(100*size if spring else 20*size))
        q=(np.arange(n)+.5)/n if spring else base.ring.RingGrid(C,n).x/C;q0=float(q[0]);u,v=(state[:n],state[n:-1]) if spring else (state[0],state[1])
        core=base.core_info(q,u,q0);M1=C*np.mean(u*np.exp(2j*np.pi*q));offset=float(np.degrees(np.angle(M1))-360*q0)
        basal=dict(model=model,size_parameter=size,cells=n,basal_circumference=C,input=str(path.relative_to(run)),input_sha256=sha(path),
                   original_path=str(record["original"]),original_sha256=sha(record["original"]),qualification=record["qualification"],
                   preparation_axis_fraction=q0,initial_M1_absolute_phase_deg=float(np.degrees(np.angle(M1))),initial_M1_offset_from_preparation_deg=offset)
        basals.append(basal)
        relevant=["code/mechanisms/spring_ring_model_20261019.py","code/mechanisms/spring_ring_closure_20261017.py","code/core/models/spring.py"] if spring else ["code/mechanisms/ring_distance_prototype_20261017.py",f"code/core/models/{model}.py"]+[p for p in sources if p.startswith("code/core/src/polarity1d/")]
        dynamics={rel:sources[rel] for rel in relevant};width=.2 if spring else 3.;peak=base.PEAKS[model]
        for index,d in enumerate(SPRING_D if spring else CHEMICAL_D):
            fraction=d/C;assert 0<fraction<=.5
            centre=(q0+fraction)%1;label=f"d_{index+1}"
            if spring:signal,dC=sizes.physical_raised_cosine(sizes.Closure(size,n),state,centre,peak,width);evaluation="continuous physical raised cosine, fixed material centre, no snapping"
            else:signal,evaluation=base.chemical_cue(base.ring.RingGrid(C,n),fraction,peak,width);dC=np.zeros(n)
            cuefile=run/"inputs"/f"{model}_C{size:g}_{label}_cue.npz"
            np.savez_compressed(cuefile,material_positions=q,initial_cue=signal,initial_cue_dC=dC,centre_fraction=centre,initial_state=state,prescribed_angle_radians=2*np.pi*fraction,prescribed_initial_arc=d)
            inputs[str(cuefile.relative_to(run))]=sha(cuefile)
            g={}
            for level,suffix in ((.5,"50"),(.1,"10")):
                edges=core_at_level(q,u,q0,level);a=geometry(edges[0]*C,edges[1]*C,d,width,C)
                g.update({f"initial_edge_left_{suffix}":a["left"],f"initial_edge_right_{suffix}":a["right"],f"initial_core_width_{suffix}":a["width"],
                          f"signed_centre_distance_from_positive_edge_{suffix}":a["edge_distance"],f"centre_relative_to_core_{suffix}":a["centre"],
                          f"footprint_relative_to_core_{suffix}":a["footprint"],f"fraction_cue_support_inside_core_{suffix}":a["overlap"]})
            patch=base.interval_weights(q,centre,width/(2*C));excess=C*np.mean(u-core["background"])
            context=dict(model=model,size_parameter=size,size_label=f'C0={size:g}' if spring else f'C={size:g}',circumference=C,resting_circumference=size if spring else None,cells=n,
                         arm="narrow_size_dense_initial_physical_arc",placement_label=label,angle_label=label,initial_arc_grid_index=index,common_initial_arc=d,physical_distance_at_onset=d,
                         fraction=fraction,requested_angle_deg=360*fraction,requested_angle_radians=2*np.pi*fraction,preparation_axis_fraction=q0,preparation_axis_absolute_deg=360*q0,
                         initial_M1_absolute_phase_deg=float(np.degrees(np.angle(M1))),initial_M1_offset_from_preparation_deg=offset,target_angle_in_P_coordinates_deg=360*fraction-offset,
                         cue_centre_fraction=centre,cue_centre_physical_at_onset=centre*C,full_physical_support=width,initial_angular_support_deg=360*width/C,requested_peak=peak,
                         initial_sampled_cue_peak=float(signal.max()),initial_cue_integral=float(C*signal.mean()),initial_cue_shape_sha256=hashlib.sha256(np.asarray(signal,dtype="<f8").tobytes()).hexdigest(),cue_evaluation=evaluation,
                         original_core_left_physical=core["left"]*C,original_core_right_physical=core["right"]*C,original_core_physical_width=core["width"]*C,
                         cue_support_core_intersection=g["fraction_cue_support_inside_core_50"]*width,overlap_fraction_cue_support=g["fraction_cue_support_inside_core_50"],
                         overlap_fraction_original_core=g["fraction_cue_support_inside_core_50"]*width/(core["width"]*C),initial_overlap=g["fraction_cue_support_inside_core_50"]>1e-12,
                         separate_challenger_interpretation_eligible=g["fraction_cue_support_inside_core_50"]<=1e-12,cue_support_contains_preparation_axis=d<=width/2,
                         initial_unit_peak_cue_excess_fraction_total=float(C*np.mean(signal/peak*(u-core["background"]))/excess),
                         initial_geometric_support_excess_fraction_total=float(C*((u-core["background"])@patch)/excess),
                         initial_cue_weighted_active=float(signal@u/signal.sum()),initial_cue_weighted_inactive=float(signal@v/signal.sum()),
                         input=str(path.relative_to(run)),input_sha256=sha(path),cue_input=str(cuefile.relative_to(run)),cue_input_sha256=sha(cuefile),fixed_scales=record["fixed_scales"],**g)
            contexts.append(context)
            for duration in ((.1,1.,10.,100.,1000.) if spring else (1.,10.,100.,300.,1000.)):
                case=dict(**context,condition_id=len(cases),cue_duration=duration,recovery=13000.,native_units_uncalibrated=True)
                options=[r for r in available if r["model"]==model and r["size_parameter"]==size and r["cue_duration"]==duration and abs(r["fraction"]-fraction)<1e-14]
                proof_row=None;seen=set()
                for candidate in options:
                    if candidate["actual_cue_path"] in seen:continue
                    seen.add(candidate["actual_cue_path"])
                    try:
                        proof_row=arc.reuse_proof(case,candidate,run,dynamics)
                        assert read(Path(proof_row["source_run"])/"manifest.json")["solver"]==previous["solver"]
                        proof_row.update(identical_solver_configuration=True,source_summary=candidate["candidate_table"],source_summary_condition_id=candidate["candidate_condition_id"])
                        break
                    except Exception as error:rejected.append(dict(condition_id=case["condition_id"],candidate=candidate["actual_cue_path"],reason=f"{type(error).__name__}: {error}"))
                if proof_row:case.update(execution="exact_archive_reuse",reuse=proof_row)
                else:case.update(execution="new_native_branch",worker_slot=len(jobs));jobs.append(case)
                cases.append(case)
    assert len(cases)==645 and len(contexts)==129 and len(basals)==12
    # Multiple separated cues must be available at the conservative 10% edge
    # simultaneously for every size, not only at one selected basal.
    external={}
    for model in base.MODELS:
        selected=[c for c in contexts if c["model"]==model]
        values=sorted({c["physical_distance_at_onset"] for c in selected})
        external[model]=[d for d in values if all(c["footprint_relative_to_core_10"]=="outside" for c in selected if c["physical_distance_at_onset"]==d)]
        assert len(external[model])>=3,(model,external[model])
    manifest=dict(generated_at=now(),experiment="Authorized narrow-size dense-initial-arc stationary-cue experiment",requested_conditions=len(cases),cases=cases,jobs=jobs,basals=basals,
                  source_sha256=sources,input_sha256=inputs,source_table_sha256=table_sources,rejected_reuse_candidates=rejected,
                  fixed_chemical_sizes=[30.,35.,40.],resting_Spring_sizes=[2.5,2.75,3.],initial_arc_grids=dict(chemical=list(CHEMICAL_D),spring=list(SPRING_D)),
                  independent_basal_qualification=str(proof),independent_basal_qualification_sha256=sha(proof),basal_stage_manifest=str(PREP/"manifest.json"),basal_stage_manifest_sha256=sha(PREP/"manifest.json"),
                  common_fully_outside_positions_at_10percent=external,overlap_cases_retained=True,excluded_for_overlap=0,
                  solver=previous["solver"],cue_full_physical_support=dict(chemical=3.,spring=.2),recovery=13000.,native_worker_implementation=previous["native_worker_implementation"],
                  angular_reference=previous["angular_reference"],Spring_coordinate_rule=previous["Spring_coordinate_rule"],
                  interpretation="Existing qualified native basals and unchanged kinetics; no size/kinetic retuning. Internal, straddling and fully outside initial footprints retained and classified independently at50/10% contrast. Activity-vector rates are not peak translation; weak R/near-antipodal projection require explicit context.")
    write(run/"manifest.json",manifest);pd.DataFrame(contexts).to_csv(run/"results/initial_context.csv",index=False)
    write(run/"freeze.json",dict(generated_at=now(),requested=len(cases),exact_reuse=len(cases)-len(jobs),new_native=len(jobs),qualified_basals=12,initial_placements=len(contexts),
                                initial_overlap_conditions=sum(c["initial_overlap"] for c in cases),sources=len(sources),inputs=len(inputs),manifest_sha256=sha(run/"manifest.json"),common_fully_outside_positions_at_10percent=external))
    collect(run);print(json.dumps(read(run/"freeze.json"),indent=2))


def collect(run):
    plan=read(run/"manifest.json");rows=[]
    for case in plan["cases"]:
        row={k:v for k,v in case.items() if k not in ("reuse",)}
        if case["execution"]=="exact_archive_reuse":row.update(status="complete",**{k:v for k,v in case["reuse"].items() if k.startswith("actual_")})
        else:
            path=run/"workers"/f'{case["worker_slot"]:03d}'/"status.json";status=read(path) if path.exists() else dict(status="not_started")
            row["status"]=status["status"]
            if status["status"]=="complete":
                for phase in status["phases"]:
                    prefix="actual_cue" if phase["kind"]=="challenge" else "actual_recovery";row[prefix+"_path"]=str(run/phase["path"]);row[prefix+"_sha256"]=phase["sha256"]
        rows.append(row)
    path=run/"results/condition_inventory.csv";tmp=path.with_suffix(".csv.tmp");pd.DataFrame(rows).to_csv(tmp,index=False);tmp.replace(path)
    counts=dict(generated_at=now(),requested=len(rows),complete=sum(r["status"]=="complete" for r in rows),exact_reuse=len(rows)-len(plan["jobs"]),new_native=len(plan["jobs"]),
                states=pd.Series([r["status"] for r in rows]).value_counts().to_dict())
    write(run/"results/collection.json",counts);print(json.dumps(counts,indent=2))


def submit(run):
    assert not(run/"launch.json").exists();plan=read(run/"manifest.json");assert sha(run/"manifest.json")==read(run/"freeze.json")["manifest_sha256"]
    for rel,digest in plan["source_sha256"].items():assert sha(run/"frozen_sources"/rel)==digest
    for rel,digest in plan["input_sha256"].items():assert sha(run/rel)==digest
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([arc.PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:05:00",f'--array=0-{len(plan["jobs"])-1}',"--job-name=arc-narrow-dense",
          "--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command]
    job=subprocess.check_output(args,text=True).strip();write(run/"launch.json",dict(generated_at=now(),array_job_id=job,sbatch_arguments=args,manifest_sha256=sha(run/"manifest.json")))
    print(json.dumps(read(run/"launch.json"),indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("action",choices=["freeze","submit","worker","collect"]);p.add_argument("--run",type=Path,default=DEFAULT);p.add_argument("--slot",type=int);p.add_argument("--qualification-proof",type=Path)
    a=p.parse_args();run=a.run.resolve()
    if a.action=="freeze":freeze(run,a.qualification_proof.resolve() if a.qualification_proof else None)
    elif a.action=="worker":sizes.worker(run,a.slot)
    else:globals()[a.action](run)
