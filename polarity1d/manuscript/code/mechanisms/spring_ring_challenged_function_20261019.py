"""Matched active-loss assay at exact Spring-ring cue-on t10 states.

The existing basal campaign remains immutable. Every challenged measurement
uses its basal core/taper and the same ongoing physical-support cue, including
the cue's circumference derivative in the analytic Jacobian.
"""
from pathlib import Path
import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
sys.dont_write_bytecode=True
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import spring_ring_function_20261019 as basal
from spring_ring_model_20261019 import Closure,physical_raised_cosine

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/"data/runs/revision_spring_challenged_function_20261019"
SOURCE=ROOT/"data/runs/revision_spring_distance_duration_20261019"
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/spring_challenged_function"
PYTHON=basal.PYTHON
read,write,sha,csvwrite=basal.read,basal.write,basal.sha,basal.csvwrite


def launch(run,source,topic):
    if run.exists():raise FileExistsError("Preserve prior native output; fresh run required")
    for name in ("inputs","frozen_sources","workers","logs","results"):(run/name).mkdir(parents=True)
    for name in ("vectors","data","metadata","viewers"):(topic/"supporting"/name).mkdir(parents=True,exist_ok=True)
    paths=[Path(__file__),Path(basal.__file__),Path(__file__).with_name("spring_ring_model_20261019.py"),
        Path(__file__).with_name("spring_ring_closure_20261017.py"),ROOT/"code/core/models/spring.py",ROOT/"code/core/models/__init__.py"]
    frozen={}
    for p in paths:
        rel=p.relative_to(ROOT);dest=run/"frozen_sources"/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest);frozen[str(rel)]=sha(dest)
    sm=read(source/"manifest.json");jobs=[j for j in sm["jobs"] if j["cue_duration"]==1000.]
    assert len(jobs)==14
    cases=[];external={str(source/"manifest.json"):sha(source/"manifest.json")}
    for job in jobs:
        status_path=source/"workers"/f'{job["slot"]:03d}'/"status.json";status=read(status_path)
        assert status["complete"] and status["status"]=="complete"
        phase=next(p for p in status["phases"] if p["kind"]=="challenge")
        archive=source/phase["path"];assert sha(archive)==phase["sha256"]
        source_basal=source/job["input"];assert sha(source_basal)==job["input_sha256"]
        with np.load(archive) as z:
            exact=np.flatnonzero(z["time"]==10.);assert len(exact)==1,"Missing exact native t10; no interpolation allowed"
            sample_index=int(exact[0]);state=z["states"][sample_index].copy()
        with np.load(source_basal) as z:reference=z["state"].copy()
        model=Closure(job["resting_circumference"],job["cells"]);geom=basal.geometry(model,reference)
        for epsilon in basal.EPS:basal.kick(model,state,geom["direction"],epsilon)
        target=run/"inputs"/f'case_{len(cases):02d}.npz'
        np.savez_compressed(target,state=state,basal_state=reference,source_time_native=10.)
        cases.append(dict(slot=len(cases),model="spring",stage="challenged_t10",cue_on=True,time_from_onset=10.,
            resting_circumference=job["resting_circumference"],cells=job["cells"],fraction=job["fraction"],
            cue_center_fraction=job["cue_center_fraction"],cue_peak=sm["cue_peak"],cue_full_physical_support=sm["cue_full_physical_support"],
            input=str(target.relative_to(run)),input_sha256=sha(target),source_native_job_slot=job["slot"],
            source_archive=str(archive),source_archive_sha256=sha(archive),source_sample_index=sample_index,source_exact_time=10.,
            source_basal=str(source_basal),source_basal_sha256=sha(source_basal)))
        for p in (status_path,archive,source_basal):external[str(p)]=sha(p)
    manifest=dict(generated_at=basal.now(),experiment="Spring ring matched common active-loss assay at exact cue-on t10",cases=cases,
        exclusions=[e for e in sm["exclusions"] if e.get("cue_duration")==1000.],
        frozen_source_sha256=frozen,external_source_sha256=external,topic=str(topic),
        primary_tau=basal.TAU,epsilons=basal.EPS,primary_epsilon=.005,
        sample_times=np.unique(np.r_[0.,np.logspace(-7,np.log10(basal.TAU),121),basal.TAU]).tolist(),
        tight=dict(rtol=1e-10,atol=1e-12,max_step=.01),tighter=dict(rtol=1e-12,atol=1e-14,max_step=.005),
        gate=dict(error_max=1e-5,linearity_relative=.05,linearity_error_allowance=3.,saturation_max=.2,prefix=True),
        cue_context="Physical support remains0.2 in evolving C; both sham/probe retain cue peak4 and material centre. Exact cue_dC enters analytic Jacobian.",
        disturbance="Frozen original basal half-contrast material core and cosine-tapered basal excess. Active loss transferred pointwise to inactive; initial C unchanged.",
        comparison="Matched counterpart of chemical common active-loss t10 assay; not the older inactive-withdrawal coefficient and not general stability.")
    write(run/"manifest.json",manifest)
    runner=run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command=shlex.join([PYTHON,str(runner),"worker","--run",str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    args=["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1","--mem=3G","--time=00:30:00","--array=0-13",
        "--job-name=Spring-challenged-function","--output="+str(run/"logs/%A_%a.out"),"--error="+str(run/"logs/%A_%a.err"),"--wrap","exec "+command]
    jid=subprocess.check_output(args,text=True).strip();write(run/"launch.json",dict(generated_at=basal.now(),array_job_id=jid,manifest_sha256=sha(run/"manifest.json"),command=args))
    write(run/"results/status.json",dict(generated_at=basal.now(),complete=False,producer="submitted",cases=14))
    print(json.dumps(dict(array_job_id=jid,run=str(run),cases=len(cases)),indent=2))


def verify(run):
    m=read(run/"manifest.json")
    for p,h in m["frozen_source_sha256"].items():assert sha(run/"frozen_sources"/p)==h
    for c in m["cases"]:assert sha(run/c["input"])==c["input_sha256"]
    return m


def integrate(model,state,case,times,settings):
    def signal(y):return physical_raised_cosine(model,y,case["cue_center_fraction"],case["cue_peak"],case["cue_full_physical_support"])
    solution=solve_ivp(lambda t,y:model.rhs(y,signal(y)[0]),(0.,float(times[-1])),state,method="BDF",t_eval=times,
        jac=lambda t,y:model.jacobian(y,*signal(y)),**settings)
    if not solution.success:raise RuntimeError(solution.message)
    states=solution.y.T;assert np.isfinite(states).all() and states[:,:-1].min()>=-5e-10
    amount=states[:,-1]*np.mean(states[:,:model.cells]+states[:,model.cells:-1],axis=1)
    error=float(abs(amount/amount[0]-1.).max());assert error<=1e-8
    return states,dict(nfev=solution.nfev,njev=solution.njev,nlu=solution.nlu,physical_amount_relative_drift=error)


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Native trajectory integration requires SLURM")
    m=verify(run);case=m["cases"][slot];folder=run/"workers"/f"{slot:02d}";folder.mkdir(parents=True,exist_ok=True)
    write(folder/"status.json",dict(generated_at=basal.now(),complete=False,producer="running",slurm_job_id=os.environ["SLURM_JOB_ID"]))
    with np.load(run/case["input"]) as z:state=z["state"].copy();reference=z["basal_state"].copy()
    model=Closure(case["resting_circumference"],case["cells"]);geom=basal.geometry(model,reference);times=np.asarray(m["sample_times"])
    arrays=dict(time=times,initial_state=state,basal_state=reference,weights_material=geom["weights_material"],direction=geom["direction"],material_positions=model.material_positions)
    numerical={}
    for accuracy in ("tight","tighter"):
        for label,epsilon in [("sham",None)]+[(f"eps{e:g}",e) for e in basal.EPS]:
            initial=state.copy() if epsilon is None else basal.kick(model,state,geom["direction"],epsilon)
            arrays[accuracy+"_"+label],numerical[accuracy+"_"+label]=integrate(model,initial,case,times,m[accuracy])
        np.testing.assert_array_equal(arrays[accuracy+"_sham"],arrays[accuracy+"_eps0"])
    archive=folder/"recovery.npz";np.savez_compressed(archive,**arrays)
    write(folder/"status.json",dict(generated_at=basal.now(),complete=True,producer="finished",slurm_job_id=os.environ["SLURM_JOB_ID"],
        archive=str(archive),archive_sha256=sha(archive),manifest_sha256=sha(run/"manifest.json"),numerical=numerical,
        full_state_zero_kick_identity=True,cue_and_circumference_derivative_used=True))


def stage_features(model,state,reference,geom,case):
    n=model.cells;u,v,C=model.fields(state);ub,vb,Cb=model.fields(reference)
    cue,cue_dC=physical_raised_cosine(model,state,case["cue_center_fraction"],case["cue_peak"],case["cue_full_physical_support"])
    basal_cue,_=physical_raised_cosine(model,reference,case["cue_center_fraction"],case["cue_peak"],case["cue_full_physical_support"])
    gain=float((cue*v).sum()/cue.sum());basal_gain=float((basal_cue*vb).sum()/basal_cue.sum())
    increment=model.rhs(state,cue)[:n]-model.rhs(state,0.)[:n]
    np.testing.assert_allclose(increment,cue*v,rtol=2e-12,atol=2e-12)
    peaks,_=find_peaks(np.tile(u,3),prominence=.05);peaks=peaks[(peaks>=n)&(peaks<2*n)]-n
    resolved=[int(i) for i in peaks if u[i]>=geom["threshold"] and geom["weights_material"][i]>0]
    inside=geom["weights_material"]>0;high=u>=geom["threshold"]
    old_peak=float(u[inside].max());core=float(C*((u-geom["background"])@geom["weights_material"]))
    old_reference=float(Cb*((ub-geom["background"])@geom["weights_material"]))
    return dict(physical_circumference_at_probe=C,basal_physical_circumference=Cb,
        cue_gain_native=gain,basal_cue_gain_native=basal_gain,gain_relative_to_same_basal= gain/basal_gain,
        native_increment_identity_max_error=float(abs(increment-cue*v).max()),cue_circumference_derivative_maxabs=float(abs(cue_dC).max()),
        initial_original_core_excess_ratio=core/old_reference,original_core_peak_above_initial_halfcontrast=old_peak>=geom["threshold"],
        original_region_resolved_peak_available=bool(resolved),original_region_resolved_peak_cells="/".join(map(str,resolved)),
        response_scope="original-region-associated peak present; identity untracked" if resolved else "regional conditional response; no resolved original-region peak",
        threshold_regions=int(np.sum(high&~np.roll(high,1))),
        initial_cue_core_overlap=bool(np.any((basal_cue>0)&inside)),
        current_cue_core_overlap=bool(np.any((cue>0)&inside)),
        core_width_native_at_probe=float(C*geom["weights_material"].sum()),
        initial_short_edge_gap=case["fraction"]*Cb-.5*geom["width_native"]-.5*case["cue_full_physical_support"])


def collect(run):
    m=verify(run);topic=Path(m["topic"]);rows=[];curves=[];features=[];audits=[];missing=[]
    reference_path=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/spring_function/supporting/data/spring_basal_recovery.csv"
    references=list(csv.DictReader(reference_path.open()))
    for case in m["cases"]:
        p=run/"workers"/f'{case["slot"]:02d}'/"status.json"
        if not p.exists() or not read(p).get("complete"):missing.append(case["slot"]);continue
        status=read(p);archive=Path(status["archive"]);assert sha(archive)==status["archive_sha256"]
        assert status["manifest_sha256"]==sha(run/"manifest.json")
        with np.load(archive) as z:arrays={k:z[k].copy() for k in z.files}
        with np.load(run/case["input"]) as z:
            np.testing.assert_array_equal(arrays["initial_state"],z["state"]);np.testing.assert_array_equal(arrays["basal_state"],z["basal_state"])
        model=Closure(case["resting_circumference"],case["cells"]);geom=basal.geometry(model,arrays["basal_state"])
        np.testing.assert_array_equal(arrays["weights_material"],geom["weights_material"]);np.testing.assert_array_equal(arrays["direction"],geom["direction"])
        for accuracy in ("tight","tighter"):
            np.testing.assert_array_equal(arrays[accuracy+"_sham"],arrays[accuracy+"_eps0"])
            for epsilon in basal.EPS:np.testing.assert_array_equal(arrays[f"{accuracy}_eps{epsilon:g}"][0],basal.kick(model,arrays["initial_state"],geom["direction"],epsilon))
        measured,history,gates=basal.evaluate(model,arrays,m["primary_tau"])
        feature=stage_features(model,arrays["initial_state"],arrays["basal_state"],geom,case)
        for row in measured:
            ref=next(r for r in references if float(r["resting_circumference"])==case["resting_circumference"] and float(r["epsilon"])==row["epsilon"])
            rows.append(dict(slot=case["slot"],model="spring",stage="challenged_t10",resting_circumference=case["resting_circumference"],
                fraction=case["fraction"],cue_on=True,archive=str(archive),archive_sha256=status["archive_sha256"],
                F_basal=float(ref["F"]),F_over_same_basal=row["F"]/float(ref["F"]),**row,**feature))
        for row in history:curves.append(dict(slot=case["slot"],resting_circumference=case["resting_circumference"],fraction=case["fraction"],**row))
        features.append(dict(slot=case["slot"],resting_circumference=case["resting_circumference"],fraction=case["fraction"],**feature))
        audits.append(dict(slot=case["slot"],source_exact_time=10.,cue_on=True,exact_saved_state=True,
            zero_kick_all_fields_identity=True,conservative_pointwise_transfer=True,initial_circumference_unchanged=True,
            cue_dC_used=True,**gates))
    if missing:print(json.dumps(dict(complete=False,missing=missing)));return
    assert len(rows)==42 and len(audits)==14
    data=topic/"supporting/data";meta=topic/"supporting/metadata"
    csvwrite(data/"spring_challenged_t10_recovery.csv",rows);csvwrite(data/"spring_challenged_t10_curves.csv",curves)
    csvwrite(data/"spring_challenged_t10_gain_geometry.csv",features)
    csvwrite(data/"unavailable_original_conditions.csv",m["exclusions"])
    result=dict(generated_at=basal.now(),complete=True,passed=all(r["passes"] for r in audits),
        common_clock_native=m["primary_tau"],matched_challenged_cases=14,new_native_trajectories=140,
        resolved_original_region_peak_absent=sum(not r["original_region_resolved_peak_available"] for r in features),
        fixed_basal_disturbance_and_material_core=True,cue_on_in_both_sham_and_probe=True,moving_C_cue_Jacobian_derivative=True,
        cases=audits,scope="Matched chemical common active-loss t10 counterpart; no inactive-pool withdrawal coefficients or off-epoch extensions.")
    write(meta/"verification.json",result);write(run/"results/verification.json",result)
    write(run/"results/status.json",dict(generated_at=basal.now(),complete=True,producer="finished",all_cases_qualified=result["passed"],cases=14))
    render(topic,rows)
    message="""# Spring ring challenged-state common response\n\nThe 14 eligible 1000-unit native branches supply exact saved states at cue time10. The opposing physical-support0.2 cue (peak4, fixed material centre) remains ON in both sham and probe. Its circumference derivative enters the analytic Jacobian as C evolves. No state interpolation is used.\n\nThe conservative activity-loss perturbation uses each size's original basal half-contrast material core and cosine-tapered basal excess direction, with epsilon0.001/0.005/0.01. Initial circumference is unchanged; all native chemical and circumference states evolve. This is the same active-loss assay as the chemical-ring t10 counterpart, not the older inactive-withdrawal coefficient.\n\n![Common challenged response](spring_challenged_t10_function.png)\n\nThe table preserves every empirical numerical, prefix-saturation and amplitude-linearity gate, plus initial/current cue overlap, actual sham/probe changes and availability of an original-region-associated resolved peak. Where no such peak is resolved, the coefficient is explicitly a regional conditional response. Peak presence concerns fixed geography and does not establish site ancestry. Positive deficit restoration is not actual original-region maintenance under the challenger.\n\nThe fixed native clock is0.01778279410038923. Basal and challenged fractions are compared only within the same resting size. The single excluded C01,f0.30 preparation geometry remains explicit. Native units and the mean-force/stretch/dilution/diffusion ring analogue remain distinct from the interval model and from cell movement.\n\n[Complete measurements](supporting/data/spring_challenged_t10_recovery.csv) · [Gain and geometry](supporting/data/spring_challenged_t10_gain_geometry.csv) · [Unavailable condition](supporting/data/unavailable_original_conditions.csv) · [Numerical verification](supporting/metadata/verification.json)\n"""
    (topic/"README.md").write_text(message)
    artifacts=[p for p in topic.rglob("*") if p.is_file() and p.name!="manifest.json"]
    write(meta/"manifest.json",dict(generated_at=basal.now(),run_manifest_sha256=sha(run/"manifest.json"),run=str(run),
        frozen_source_sha256=m["frozen_source_sha256"],external_source_sha256=m["external_source_sha256"],
        basal_comparison_table=str(reference_path),basal_comparison_sha256=sha(reference_path),
        collector_source_path=str(Path(__file__)),collector_source_sha256=sha(Path(__file__)),
        artifact_sha256={str(p.relative_to(topic)):sha(p) for p in artifacts},new_native_trajectories=140))
    print(json.dumps(result,indent=2))


def render(topic,rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(11.7,4.7));fig.subplots_adjust(left=.07,right=.985,top=.80,bottom=.30,wspace=.42)
    for C,color in ((1.,"#3976a8"),(2.,"#c96926"),(3.,"#2a9288")):
        group=sorted([r for r in rows if r["resting_circumference"]==C and r["epsilon"]==.005],key=lambda r:r["fraction"])
        for axis,key in zip(axes,("gain_relative_to_same_basal","F_over_same_basal","sham_core_change")):
            axis.plot([r["fraction"] for r in group],[r[key] for r in group],marker="o",color=color,lw=1.3,ms=4,label=f"C₀{C:g}")
            for r in group:
                if not r["original_region_resolved_peak_available"]:axis.scatter(r["fraction"],r[key],s=46,facecolor="white",edgecolor=color,zorder=4)
    for axis,title,label in zip(axes,("A  Input gain after10 cue units","B  Conditional restoration after10","C  Actual unperturbed regional change"),
        ("Gain at t10 / own basal gain","Fτ at t10 / own basal Fτ","Sham original-core active amount change")):
        axis.set(title=title,xlabel="Initial cue offset / circumference",ylabel=label);axis.grid(axis="y",color="#eeeeee",lw=.7)
    axes[0].legend(frameon=False,fontsize=8);axes[1].axhline(1,color="#aaaaaa",ls="--",lw=.7);axes[2].axhline(0,color="#aaaaaa",lw=.7)
    fig.suptitle("Spring ring | Matched cue-on common active-loss response at exact native t10",fontsize=12,y=.98)
    fig.text(.53,.045,"14 eligible native states; same frozen basal core/taper and activity→inactive transfer. Fixed native clock τ=0.0177827941; primary ε=0.5%.\n"
        "Cue remains on in sham/probe, with full circumference feedback and physical-support derivative. Open points: no resolved original-region peak.\n"
        "Restoring a sham-relative deficit does not establish actual incumbent maintenance; positive/negative sham changes remain explicit.",ha="center",fontsize=8,linespacing=1.4)
    for extension in ("png","pdf","svg"):
        p=(topic if extension=="png" else topic/"supporting/vectors")/f"spring_challenged_t10_function.{extension}"
        fig.savefig(p,dpi=200 if extension=="png" else None,bbox_inches="tight")
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=("launch","worker","collect"));parser.add_argument("--run",type=Path,default=RUN)
    parser.add_argument("--source",type=Path,default=SOURCE);parser.add_argument("--topic",type=Path,default=TOPIC);parser.add_argument("--slot",type=int)
    args=parser.parse_args()
    if args.action=="launch":launch(args.run.resolve(),args.source.resolve(),args.topic.resolve())
    elif args.action=="worker":worker(args.run.resolve(),args.slot)
    else:collect(args.run.resolve())


if __name__=="__main__":main()
