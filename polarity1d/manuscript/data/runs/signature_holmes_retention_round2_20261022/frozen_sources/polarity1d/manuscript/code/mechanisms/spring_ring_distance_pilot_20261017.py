"""Bounded Spring ring distance pilot from one exact retained basal state."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import shlex
import shutil
import signal
import subprocess
import sys
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import BDF
from scipy.signal import find_peaks
from spring_ring_closure_20261017 import Closure
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
PREP = ROOT / "data/runs/spring_ring_preparation_pilot_20261017"
DEFAULT = ROOT / "data/runs/spring_ring_distance_pilot_20261017"
FIG = ROOT / "figures/length_mechanisms_20261008/dev/spring_ring_distance_pilot_20261017"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False,default=json_native) + "\n")
    temporary.replace(path)


def json_native(item):
    if isinstance(item,np.generic):return item.item()
    raise TypeError(f"Cannot encode {type(item).__name__}")


def signal_values(closure,state,center,peak,width=.2):
    q = closure.material_positions
    distances = np.abs((q-q[center]+.5) % 1-.5) * float(state[-1])
    return np.where(distances < width/2, peak*.5*(1+np.cos(2*np.pi*distances/width)),0.)


def sample_times(end):
    values=np.r_[np.arange(0,10,.1),np.arange(10,100,1.),np.arange(100,end,5.),float(end)]
    return np.unique(values[values <= end])


def launch(run):
    if run.exists():
        raise ValueError("Use fresh output; retained evidence is immutable")
    prep=read(PREP / "workers/status.json")
    assert prep["status"]=="qualifying_single_site" and read(PREP / "results/audit.json")["passed"]
    checkpoint=PREP / prep["checkpoint"]
    assert sha(checkpoint)==prep["checkpoint_sha256"]
    for directory in ("inputs","frozen_sources","workers","logs","results"):
        (run / directory).mkdir(parents=True)
    shutil.copy2(checkpoint,run / "inputs/basal.npz")
    source_paths=[Path(__file__),Path(__file__).with_name("spring_ring_closure_20261017.py"),
        Path(__file__).with_name("length_dev_output_layout.py"),ROOT / "code/core/models/spring.py"]
    hashes={}
    for path in source_paths:
        relative=path.relative_to(ROOT)
        target=run / "frozen_sources" / relative
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target);hashes[str(relative)]=sha(target)
    prior=read(PREP / "manifest.json")
    for relative in ("code/mechanisms/spring_ring_closure_20261017.py","code/core/models/spring.py"):
        assert hashes[relative]==prior["source_sha256"][relative]
    closure=Closure(2.,200)
    with np.load(checkpoint) as z:
        initial=z["state"].copy()
    u=initial[:200];pole=int(np.argmax(u));threshold=float(u.min()+.5*np.ptp(u));core=u>=threshold
    jobs=[]
    for slot,offset in enumerate((60,80,100,-80)):
        center=(pole+offset)%200
        support=signal_values(closure,initial,center,4)>0
        overlap=bool(np.any(core&support))
        jobs.append(dict(slot=slot,offset_cells=offset,offset_fraction=offset/200,
            center_cell=center,shortest_initial_arc=min(abs(offset),200-abs(offset))*initial[-1]/200,
            overlap_with_original_core=overlap,eligible=not overlap,
            target_region_cells=np.flatnonzero(support).tolist(),reflection_control=offset<0))
    manifest=dict(generated_at=now(),model="Spring circular global-extension analogue",
        resting_circumference=2.,cells=200,initial_circumference=float(initial[-1]),
        initial_amount=4.,incumbent_cell=pole,original_core_cells=np.flatnonzero(core).tolist(),
        original_half_contrast_threshold=threshold,jobs=jobs,maintained_cue_duration=1000.,
        recovery_duration=13000.,cue_peak=4.,total_physical_support=.2,
        source_sha256=hashes,input_sha256=sha(run / "inputs/basal.npz"),
        preparation_checkpoint=str(checkpoint),preparation_checkpoint_sha256=sha(checkpoint),
        solver=dict(method="BDF",rtol=1e-9,atol=1e-11,max_step=10.,
            jacobian="Full dense numerical Jacobian; global mean-force/dilution dependencies included"),
        sample_rule="dt0.1 first10,dt1 to100,dt5 thereafter, restarted at withdrawal",
        soft_seconds=3300,pool_relative_tolerance=1e-7,allowed_negative_roundoff=5e-10,
        mechanics="Cdot=2(mean_q F(u)-.01(C-C0)); mean after nonlinear switch; dynamic circumference, native dilution",
        cue_geometry="Fixed material/angular centre; instantaneous physical width0.2; peak4; raised cosine",
        readouts="Raw activity in fixed original material core and initial cue target mask; full fields and circular components/peaks",
        exclusions="No parameter/duration optimization; no spontaneous-probability or migration claim; no native2D equivalence")
    write(run / "manifest.json",manifest)
    eligible=[str(job["slot"]) for job in jobs if job["eligible"]]
    if not eligible:
        raise RuntimeError("No nonoverlapping challenger positions")
    runner=run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command=shlex.join(["env","PYTHONDONTWRITEBYTECODE=1","OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1","MKL_NUM_THREADS=1",PYTHON,str(runner),"worker","--run",str(run)])
    command+=' --slot "$SLURM_ARRAY_TASK_ID"'
    job_id=subprocess.check_output(["sbatch","--parsable","--qos=small","--ntasks=1","--cpus-per-task=1",
        "--mem=3G","--time=01:00:00","--signal=B:USR1@60","--array="+','.join(eligible),
        "--job-name=spring-distance","--output="+str(run / "logs/%A_%a.out"),
        "--error="+str(run / "logs/%A_%a.err"),"--wrap","exec "+command],text=True).strip()
    write(run / "launch.json",dict(generated_at=now(),array_job_id=job_id,manifest_sha256=sha(run / "manifest.json")))
    print(json.dumps(dict(job_id=job_id,run=str(run),jobs=jobs),indent=2))


def integrate(closure,initial,center,peak,end,manifest,deadline):
    function=lambda t,y:closure.rhs(y,signal_values(closure,y,center,peak,manifest["total_physical_support"]))
    solver=BDF(function,0.,initial,end,rtol=manifest["solver"]["rtol"],atol=manifest["solver"]["atol"],max_step=10.)
    wanted=sample_times(end);times=[0.];states=[initial.copy()];index=1
    minimum=float(np.min(initial[:-1]));worst=0.;status="complete";message=""
    while solver.status=="running":
        if STOP or time.monotonic()>=deadline:
            status="timeout";message="Worker deadline or SLURM advance signal";break
        message=solver.step() or ""
        if solver.status=="failed":status="numerical_error";break
        if not np.isfinite(solver.y).all():status="nonfinite_state";break
        minimum=min(minimum,float(np.min(solver.y[:-1])))
        worst=max(worst,abs(closure.amount(solver.y)/manifest["initial_amount"]-1))
        interpolant=solver.dense_output()
        while index<len(wanted) and wanted[index]<=solver.t:
            times.append(float(wanted[index]));states.append(interpolant(wanted[index]).copy());index+=1
        if minimum < -manifest["allowed_negative_roundoff"] or worst > manifest["pool_relative_tolerance"]:
            status="numerical_error";message="Negative field or physical-amount drift";break
    if times[-1]!=solver.t:
        times.append(float(solver.t));states.append(solver.y.copy())
    return dict(time=np.asarray(times),states=np.asarray(states),status=status,message=message,
        minimum_concentration=minimum,max_amount_error=worst,nfev=int(solver.nfev),njev=int(solver.njev),nlu=int(solver.nlu))


def worker(run,slot):
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Trajectory integration requires SLURM")
    signal.signal(signal.SIGUSR1,lambda *args:globals().__setitem__("STOP",True))
    signal.signal(signal.SIGTERM,lambda *args:globals().__setitem__("STOP",True))
    manifest=read(run / "manifest.json");job=manifest["jobs"][slot]
    assert job["eligible"] and sha(run / "inputs/basal.npz")==manifest["input_sha256"]
    for relative,digest in manifest["source_sha256"].items():assert sha(run / "frozen_sources" / relative)==digest
    closure=Closure(2.,200)
    with np.load(run / "inputs/basal.npz") as z:state=z["state"].copy()
    directory=run / "workers" / f"{slot:02d}";directory.mkdir()
    record=dict(generated_at=now(),producer="running",complete=False,job=job,phases=[],slurm_job_id=os.environ["SLURM_JOB_ID"])
    write(directory / "status.json",record);deadline=time.monotonic()+manifest["soft_seconds"]
    try:
        for kind,peak,end in (("maintained_cue",4.,1000.),("recovery",0.,13000.)):
            record.update(generated_at=now(),status=kind);write(directory / "status.json",record)
            result=integrate(closure,state,job["center_cell"],peak,end,manifest,deadline)
            path=directory / (kind+".npz")
            np.savez_compressed(path,time=result["time"],states=result["states"],material_positions=closure.material_positions,
                center_cell=np.asarray(job["center_cell"]),cue_peak=np.asarray(peak),total_physical_support=np.asarray(.2))
            metadata={key:value for key,value in result.items() if key not in ("time","states")}
            metadata.update(path=str(path.relative_to(run)),sha256=sha(path),phase=kind)
            record["phases"].append(metadata);state=result["states"][-1].copy()
            if result["status"]!="complete":break
        record.update(generated_at=now(),producer="finished",complete=True,status=record["phases"][-1]["status"])
    except Exception as error:
        record.update(generated_at=now(),producer="finished",complete=True,status="error",error=f"{type(error).__name__}: {error}")
        write(directory / "status.json",record);raise
    write(directory / "status.json",record);print(json.dumps(record,indent=2))


def peak_readout(u,threshold):
    mask=u>=threshold
    regions=0 if not mask.any() else (1 if mask.all() else int(np.sum(mask & ~np.roll(mask,1))))
    peaks,properties=find_peaks(np.tile(u,3),prominence=0)
    keep=(peaks>=len(u))&(peaks<2*len(u))
    peaks=peaks[keep]-len(u);prominences=properties["prominences"][keep]
    resolved=peaks[prominences>=.05]
    return dict(frozen_half_contrast_regions=regions,resolved_peaks=len(resolved),
        peak_cells=resolved.tolist(),peak_values=u[resolved].tolist(),peak_prominences=prominences[prominences>=.05].tolist())


def load_job(run,job,manifest):
    record=read(run / "workers" / f'{job["slot"]:02d}' / "status.json")
    assert record["complete"] and record["status"]=="complete"
    time_values=[];states=[];offset=0.
    for phase in record["phases"]:
        assert sha(run / phase["path"])==phase["sha256"]
        with np.load(run / phase["path"]) as z:t,y=z["time"].copy(),z["states"].copy()
        take=slice(None) if not states else slice(1,None)
        time_values.extend(t[take]+offset);states.extend(y[take]);offset+=t[-1]
    return record,np.asarray(time_values),np.asarray(states)


def collect(run,figdir):
    import csv
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    manifest=read(run / "manifest.json")
    assert sha(run / "manifest.json")==read(run / "launch.json")["manifest_sha256"]
    for relative,digest in manifest["source_sha256"].items():assert sha(run / "frozen_sources" / relative)==digest
    with np.load(run / "inputs/basal.npz") as z:initial=z["state"].copy()
    assert sha(run / "inputs/basal.npz")==manifest["input_sha256"]==manifest["preparation_checkpoint_sha256"]
    core=np.asarray(manifest["original_core_cells"]);n=manifest["cells"]
    jobs=[];rows=[];summaries=[];source_hashes={"manifest.json":sha(run / "manifest.json"),"inputs/basal.npz":sha(run / "inputs/basal.npz")}
    for job in manifest["jobs"]:
        if not job["eligible"]:continue
        record,t,y=load_job(run,job,manifest)
        assert np.isfinite(y).all() and y[:,:-1].min()>=-manifest["allowed_negative_roundoff"]
        np.testing.assert_array_equal(y[0],initial)
        target=np.asarray(job["target_region_cells"])
        amount=y[:,-1]*np.mean(y[:,:n]+y[:,n:-1],axis=1)
        assert np.max(np.abs(amount/4-1))<1e-7
        old=np.mean(y[:,core],axis=1);new=np.mean(y[:,target],axis=1)
        features=[]
        for index,(tt,state) in enumerate(zip(t,y)):
            feature=peak_readout(state[:n],manifest["original_half_contrast_threshold"]);features.append(feature)
            rows.append(dict(slot=job["slot"],offset_fraction=job["offset_fraction"],time=tt,cue_on=tt<=1000,
                circumference=state[-1],old_core_mean=old[index],target_mean=new[index],
                old_center_active=state[manifest["incumbent_cell"]],target_center_active=state[job["center_cell"]],
                target_region_peak=float(np.max(state[target])),
                old_core_amount=state[-1]/n*np.sum(state[core]),target_amount=state[-1]/n*np.sum(state[target]),
                primary_peak=state[:n].max(),primary_minimum=state[:n].min(),
                frozen_half_contrast_regions=feature["frozen_half_contrast_regions"],resolved_peaks=feature["resolved_peaks"],
                physical_amount=amount[index],amount_relative_error=amount[index]/4-1))
        on=t<=1000;off=t>=1000
        two=np.array([f["resolved_peaks"]>=2 for f in features]) & on
        two_times=t[two]
        separated=np.array([f["frozen_half_contrast_regions"]>=2 for f in features]) & on
        separated_times=t[separated]
        summary=dict(offset_fraction=job["offset_fraction"],initial_shortest_arc=job["shortest_initial_arc"],
            original_core_cells=len(core),target_cells=len(target),original_overlap=job["overlap_with_original_core"],
            max_amount_relative_error=float(np.max(np.abs(amount/4-1))),
            max_cue_resolved_peaks=int(max(f["resolved_peaks"] for f,selected in zip(features,on) if selected)),
            sampled_two_resolved_peak_range=[float(two_times[0]),float(two_times[-1])] if len(two_times) else [],
            sampled_two_half_contrast_component_range=[float(separated_times[0]),float(separated_times[-1])] if len(separated_times) else [],
            first_target_mean_at_least_one=float(t[np.flatnonzero(new>=1)[0]]) if np.any(new>=1) else None,
            cue_end=next(rows[-len(t)+i] for i,tt in enumerate(t) if tt==1000),final=rows[-1],
            final_peak_features=features[-1],circumference_range=[float(y[:,-1].min()),float(y[:,-1].max())])
        translated=np.r_[np.roll(initial[:n],job["offset_cells"]),np.roll(initial[n:-1],job["offset_cells"]),initial[-1]]
        summary["final_profile_max_residual_vs_translated_initial"]=float(np.max(np.abs(y[-1,:-1]-translated[:-1])))
        summary["final_circumference_residual_vs_initial"]=float(y[-1,-1]-initial[-1])
        summaries.append(summary);jobs.append((job,t,y,old,new,features))
        source_hashes[f'workers/{job["slot"]:02d}/status.json']=sha(run / "workers" / f'{job["slot"]:02d}' / "status.json")
        for phase in record["phases"]:source_hashes[phase["path"]]=phase["sha256"]
    control=None
    if len(jobs)==4:
        positive=jobs[1];negative=jobs[3]
        np.testing.assert_array_equal(positive[1],negative[1])
        reflection=(-np.arange(n)+2*manifest["incumbent_cell"])%n
        transformed=np.c_[negative[2][:,reflection],negative[2][:,n+reflection],negative[2][:,-1]]
        delta=np.abs(positive[2]-transformed)
        control=dict(offsets=[.4,-.4],max_absolute_chemical_difference=float(delta[:,:-1].max()),
                     max_absolute_circumference_difference=float(delta[:,-1].max()))
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.titlesize":9,"axes.labelsize":8,
        "xtick.labelsize":7,"ytick.labelsize":7,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(8.0,7.5));grid=fig.add_gridspec(3,3,height_ratios=[1.2,1.,.72],left=.08,right=.94,bottom=.09,top=.83,wspace=.35,hspace=.68)
    shift=n//4-manifest["incumbent_cell"];q=(np.arange(n)+.5)/n
    primary_jobs=[entry for entry in jobs if not entry[0]["reflection_control"]]
    lower=min(entry[2][:,:n].min() for entry in primary_jobs);upper=max(entry[2][:,:n].max() for entry in primary_jobs)
    old_color,new_color="#235789","#b44436";heat_axes=[]
    for column,(job,t,y,old,new,features) in enumerate(primary_jobs):
        axis=fig.add_subplot(grid[0,column]);heat_axes.append(axis);mesh=axis.pcolormesh(q,t,np.roll(y[:,:n],shift,axis=1),shading="auto",cmap="viridis",vmin=lower,vmax=upper,rasterized=True)
        axis.set_yscale("symlog",linthresh=10);axis.set_ylim(0,t[-1]);axis.axhline(1000,color="white",ls=":",lw=1)
        axis.axvline(q[(job["center_cell"]+shift)%n],color="#ea8678",ls=":",lw=.9)
        axis.set(title=f'{chr(65+column)}   Initial d/C = {job["offset_fraction"]:.1f}\nd = {job["shortest_initial_arc"]:.3f}',xlabel="Material position around ring")
        if column==0:axis.set_ylabel("Time since cue onset")
        axis=fig.add_subplot(grid[1,column]);axis.plot(t,old,color=old_color,label="Original core");axis.plot(t,new,color=new_color,label="Cue target")
        axis.axvline(1000,color="#666666",ls=":",lw=1);axis.set_xscale("symlog",linthresh=10);axis.set_xlim(0,t[-1])
        axis.set_ylim(.65,1.7)
        axis.set(title=f'{chr(68+column)}   Fixed-region activity',xlabel="Time since cue onset")
        if column==0:axis.set_ylabel("Mean active Rac");axis.legend(frameon=False,fontsize=7)
    bounds=heat_axes[0].get_position();color_axis=fig.add_axes([.955,bounds.y0,.012,bounds.height]);fig.colorbar(mesh,cax=color_axis,label="Active Rac")
    axis=fig.add_subplot(grid[2,:])
    for (job,t,y,old,new,features),color in zip(primary_jobs,["#235789","#c29c51","#6c5480"]):
        axis.plot(t,y[:,-1],color=color,label=f'd/C = {job["offset_fraction"]:.1f}')
    axis.axvline(1000,color="#666666",ls=":",lw=1);axis.set_xscale("symlog",linthresh=10);axis.set_xlim(0,max(entry[1][-1] for entry in primary_jobs))
    axis.set(title="G   Circumference responds during the contest",xlabel="Time since cue onset",ylabel="Circumference C");axis.legend(frameon=False,ncol=3,fontsize=7)
    fig.suptitle("Spring circular global-extension analogue\nSeparation changes early site competition",fontsize=11,y=.98)
    files=[]
    for extension in ("png","pdf","svg"):
        path=dev_output_path(figdir,"spring_ring_distance_pilot."+extension);fig.savefig(path,dpi=250,facecolor="white",bbox_inches="tight");files.append(path)
    plt.close(fig)
    csv_path=dev_output_path(figdir,"fixed_region_timecourses.csv")
    with csv_path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    files.append(csv_path)
    result=dict(generated_at=now(),complete=True,initial_checkpoint_sha256=manifest["input_sha256"],jobs=summaries,reflection_control=control,
        scope="Bounded deterministic trajectories at fixed initial circumference, one cue peak/support/duration; no parameter optimization, noise, migration or universal threshold")
    write(run / "results/summary.json",result)
    path=dev_output_path(figdir,"summary.json");write(path,result);files.append(path)
    caption=("Spring circular global-extension analogue, first distance pilot. Every trial starts from the same exact qualified C0 = 2, N = 200 checkpoint "
        "with actual initial C = 2.04418073. Challenger material offsets 60, 80 and 100 cells give initial shortest arcs 0.6133, 0.8177 and 1.0221; "
        "none overlaps the original half-contrast core. Cue peak 4, total instantaneous physical support 0.2, maintained 1000 units then removed for 13000. "
        "Cue centre is fixed in material coordinates while width adapts to instantaneous circumference. A–C: raw active fields; "
        "white horizontal dotted line marks withdrawal, red vertical dotted line the target position. Time uses a linear scale through 10 then logarithmic spacing. "
        "Fields are cyclically shifted only for display, placing the original incumbent near material position 0.25; endpoints connect. "
        "The near cue bridges to the incumbent above the frozen initial half-contrast level by sampled t = 0.2, although two resolved maxima remain through t = 0.8. "
        "At offsets 0.4 and 0.5, separated competitor growth accompanies old-site decay: two half-contrast components remain through t = 6.6 and 6.9, "
        "and the old maximum extinguishes near its original position rather than merging spatially into the target. Two resolved maxima remain through t = 7.7 and 8.2. "
        "All target-region maxima exceed activity 1 at the first saved sample t = 0.1. Target means cross 1 within the shared sampled bracket (0.1, 0.2]; "
        "this does not establish identical activation kinetics. D–F: mean raw activity in frozen original-core and initial target masks; "
        "late differences in original-core activity include geometric overlap with the relocated broad cap and do not show a surviving second site. "
        "G: evolved circumference, which mediates retained global feedback and can differ after the same start. "
        "A reflected −0.4 control is archived separately. Native reaction, diffusion, force, feedback and dilution coefficients remain unchanged; "
        "periodic geometry and the mean-force circumference closure are explicit adaptations, not native 2D mechanics. "
        "The full dense numerical Jacobian includes global dependencies. No dose/duration optimization or spontaneous competitor probability is inferred. Generated "+now()+".")
    path=dev_output_path(figdir,"caption.txt");path.write_text(caption+"\n");files.append(path)
    write(dev_output_path(figdir,"manifest.json"),dict(generated_at=now(),artifact_root="../..",run=str(run),source_sha256=source_hashes,
        sha256={dev_artifact_key(figdir,p):sha(p) for p in files},renderer_sha256=sha(Path(__file__)),canonical_models_changed=False))
    print(json.dumps(result,indent=2,default=json_native))


def audit(run):
    """Read-only reconstruction independent of the worker's status decisions."""
    import csv
    manifest=read(run / "manifest.json")
    assert sha(run / "manifest.json")==read(run / "launch.json")["manifest_sha256"]
    for relative,digest in manifest["source_sha256"].items():assert sha(run / "frozen_sources" / relative)==digest
    assert sha(run / "inputs/basal.npz")==manifest["input_sha256"]==manifest["preparation_checkpoint_sha256"]
    with np.load(run / "inputs/basal.npz") as z:initial=z["state"].copy()
    n=manifest["cells"];closure=Closure(2.,n);scales=np.r_[np.full(2*n,2.),2.]
    report=[]
    for job in manifest["jobs"]:
        if not job["eligible"]:continue
        directory=run / "workers" / f'{job["slot"]:02d}';record=read(directory / "status.json")
        assert record["complete"] and record["status"]=="complete"
        previous=initial.copy();worst=0.
        for phase in record["phases"]:
            assert sha(run / phase["path"])==phase["sha256"]
            with np.load(run / phase["path"]) as z:
                t,y=z["time"].copy(),z["states"].copy();peak=float(z["cue_peak"])
            np.testing.assert_array_equal(y[0],previous)
            expected_duration=1000. if phase["phase"]=="maintained_cue" else 13000.
            np.testing.assert_array_equal(t,sample_times(expected_duration))
            assert peak==(4. if phase["phase"]=="maintained_cue" else 0.)
            assert np.isfinite(y).all() and np.min(y[:,:-1])>=-5e-10
            amounts=y[:,-1]*(y[:,:n]+y[:,n:-1]).sum(axis=1)/n
            worst=max(worst,float(np.max(np.abs(amounts/4-1))))
            assert worst<1e-7
            previous=y[-1].copy()
        # Reconstruct the final unforced field/C stationarity on frozen scales.
        final=y[t>=t[-1]-1000]
        scaled_range=float(np.max(np.ptp(final,axis=0)/scales))
        scaled_drift=float(np.max(np.abs(y[-1]-y[-2])/scales)*1000/(t[-1]-t[-2]))
        rhs=float(np.max(np.abs(closure.rhs(y[-1],np.zeros(n)))/scales)*1000)
        stderr=(run / "logs" / (read(run / "launch.json")["array_job_id"]+f'_{job["slot"]}.err')).read_text()
        report.append(dict(offset_fraction=job["offset_fraction"],exact_initial_and_withdrawal_continuity=True,
            all_samples_finite=True,correct_time_grid=True,max_sampled_amount_relative_error=worst,
            final_1000_scaled_full_field_and_C_range=scaled_range,final_1000_projected_drift=scaled_drift,
            final_rhs_scaled_to_1000=rhs,final_settling_pass=bool(max(scaled_range,scaled_drift,rhs)<1e-4),stderr=stderr))
    result=dict(generated_at=now(),passed=True,trajectories=report,
        scope="Exact saved states, immutable sources, cue phases, complete sample grids, physical amount and final field/C stationarity; no independent integration replay")
    write(run / "results/audit.json",result);print(json.dumps(result,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=["launch","worker","collect","audit"])
    parser.add_argument("--run",type=Path,default=DEFAULT);parser.add_argument("--slot",type=int);parser.add_argument("--fig",type=Path,default=FIG)
    args=parser.parse_args()
    if args.action=="launch":launch(args.run.resolve())
    elif args.action=="worker":worker(args.run.resolve(),args.slot)
    elif args.action=="audit":audit(args.run.resolve())
    else:collect(args.run.resolve(),args.fig.resolve())
