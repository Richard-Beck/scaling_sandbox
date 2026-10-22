"""Reviewed Spring circular global-extension analogue: SLURM preparation only."""
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
from spring_ring_closure_20261017 import Closure, native
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/spring_ring_preparation_pilot_20261017"
FIG = ROOT / "figures/length_mechanisms_20261008/dev/spring_ring_preparation_pilot_20261017"
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
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def launch(run):
    if run.exists():
        raise ValueError("Use a fresh run directory; no retained outputs are replaced")
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run / name).mkdir(parents=True)
    sources = [Path(__file__), Path(__file__).with_name("spring_ring_closure_20261017.py"),
               Path(__file__).with_name("length_dev_output_layout.py"),
               ROOT / "code/core/models/spring.py"]
    hashes = {}
    for path in sources:
        relative = path.relative_to(ROOT)
        target = run / "frozen_sources" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = sha(target)
    closure = Closure(resting_circumference=2.0, cells=200)
    initial = closure.initial()
    np.savez_compressed(run / "inputs/uniform.npz", state=initial,
                        material_positions=closure.material_positions)
    manifest = dict(generated_at=now(), model="Spring circular global-extension analogue",
        resting_circumference=2.0, cells=200, preparation_peak=4.0,
        total_physical_support=.2, preparation_durations=[10.0, 100.0, 1000.0],
        frozen_scales=dict(active=2.0, inactive=2.0, circumference=2.0),
        rule=dict(chunk=1000.0, tolerance=1e-4, consecutive_windows=2,
                  holdout=10000.0, holdout_tolerance=1e-3,
                  max_basal_per_attempt=30000.0, allowed_negative_roundoff=5e-10,
                  pool_relative_tolerance=1e-7, minimum_contrast=.05,
                  minimum_peak_prominence=.05),
        solver=dict(method="BDF", rtol=1e-9, atol=1e-11, max_step=10.0,
                    jacobian="Full dense numerical Jacobian; no sparsity mask. Mean-force and dilution couple all active cells globally."),
        soft_seconds=3300, source_sha256=hashes,
        input_sha256=sha(run / "inputs/uniform.npz"),
        new_assumptions=["Periodic one-dimensional signalling domain",
            "Cdot=2(mean_q F(u)-.01(C-C0)); mean after nonlinear switch",
            "Uniform material deformation and one global circumference",
            "Fresh uniform start removes native interval cosine; local cue selects site"],
        unchanged="Native local reaction, force switch, diffusion, deactivation and dilution coefficients",
        native_2d_equivalence=False, challenger_trials=0,
        reviewed_by="Root reviewed equation (3) and algebraic prototype before this bounded run")
    write(run / "manifest.json", manifest)
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", PYTHON, str(runner), "worker", "--run", str(run)])
    job = subprocess.check_output(["sbatch", "--parsable", "--qos=small", "--ntasks=1",
        "--cpus-per-task=1", "--mem=3G", "--time=01:00:00", "--signal=B:USR1@60",
        "--job-name=spring-ring-prep", "--output="+str(run / "logs/%j.out"),
        "--error="+str(run / "logs/%j.err"), "--wrap", "exec "+command], text=True).strip()
    write(run / "launch.json", dict(generated_at=now(), slurm_job_id=job,
                                     manifest_sha256=sha(run / "manifest.json")))
    print(json.dumps(dict(job_id=job, run=str(run)), indent=2))


def cue_values(closure, state, peak, width):
    circumference = float(state[-1])
    q = closure.material_positions
    distance = np.abs((q-q[0]+.5) % 1.0-.5) * circumference
    return np.where(distance < width / 2,
                    peak * .5 * (1+np.cos(2*np.pi*distance/width)), 0.0)


def integrate(closure, initial, peak, end, interval, manifest, deadline):
    signal_width = manifest["total_physical_support"]
    def function(t, state):
        return closure.rhs(state, cue_values(closure, state, peak, signal_width)
                           if peak else np.zeros(closure.cells))
    settings = manifest["solver"]
    solver = BDF(function, 0.0, initial, end, rtol=settings["rtol"],
                 atol=settings["atol"], max_step=settings["max_step"])
    wanted = np.arange(0, end+interval/2, interval)
    wanted = wanted[wanted <= end]
    if wanted[-1] != end:
        wanted = np.r_[wanted, end]
    times, states = [0.0], [initial.copy()]
    index, status, message = 1, "complete", ""
    reference_amount = manifest["resting_circumference"] * 2
    worst_pool, minimum = 0.0, float(np.min(initial[:-1]))
    while solver.status == "running":
        if STOP or time.monotonic() >= deadline:
            status, message = "timeout", "Worker soft deadline or SLURM advance signal"
            break
        message = solver.step() or ""
        if solver.status == "failed":
            status = "numerical_error"
            break
        minimum = min(minimum, float(np.min(solver.y[:-1])))
        worst_pool = max(worst_pool, abs(closure.amount(solver.y)/reference_amount-1))
        dense = solver.dense_output()
        while index < len(wanted) and wanted[index] <= solver.t:
            times.append(float(wanted[index]))
            states.append(dense(wanted[index]).copy())
            index += 1
        if minimum < -manifest["rule"]["allowed_negative_roundoff"] or worst_pool > manifest["rule"]["pool_relative_tolerance"]:
            status, message = "numerical_error", "Negative concentration or physical-amount drift"
            break
    if times[-1] != solver.t:
        times.append(float(solver.t))
        states.append(solver.y.copy())
    return dict(time=np.asarray(times), states=np.asarray(states), status=status, message=message,
        max_pool_relative_error=worst_pool, minimum_concentration=minimum,
        nfev=int(solver.nfev), njev=int(solver.njev), nlu=int(solver.nlu))


def diagnostics(times, states, manifest):
    n, rule = manifest["cells"], manifest["rule"]
    keep = times >= times[-1]-rule["chunk"]
    scales = np.r_[np.full(2*n, 2.0), manifest["resting_circumference"]]
    ranges = np.ptp(states[keep], axis=0)/scales
    drift = np.abs(states[-1]-states[-2])/(times[-1]-times[-2])*rule["chunk"]/scales
    return dict(steady=bool(max(ranges.max(), drift.max()) <= rule["tolerance"]),
        active_range=float(ranges[:n].max()), inactive_range=float(ranges[n:-1].max()),
        circumference_range=float(ranges[-1]), active_projected_drift=float(drift[:n].max()),
        inactive_projected_drift=float(drift[n:-1].max()), circumference_projected_drift=float(drift[-1]))


def architecture(state, n):
    u = state[:n]
    contrast = float(np.ptp(u))
    mask = u >= u.min()+.5*contrast
    regions = 0 if contrast < 1e-12 else (1 if mask.all() else int(np.sum(mask & ~np.roll(mask, 1))))
    peaks, properties = find_peaks(np.tile(u, 3), prominence=0)
    select = (peaks >= n) & (peaks < 2*n)
    prominences = properties["prominences"][select]
    return dict(contrast=contrast, connected_regions=regions,
        resolved_peaks=int(np.sum(prominences >= .05)),
        peak_cells=(peaks[select]-n).tolist(), prominences=prominences.tolist(),
        peak=float(u.max()), minimum=float(u.min()), circumference=float(state[-1]))


def worker(run):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("All trajectory integration must run via SLURM")
    signal.signal(signal.SIGUSR1, lambda *args: globals().__setitem__("STOP", True))
    signal.signal(signal.SIGTERM, lambda *args: globals().__setitem__("STOP", True))
    manifest = read(run / "manifest.json")
    for relative, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    assert sha(run / "inputs/uniform.npz") == manifest["input_sha256"]
    closure = Closure(manifest["resting_circumference"], manifest["cells"])
    with np.load(run / "inputs/uniform.npz") as z:
        initial = z["state"].copy()
    record = dict(generated_at=now(), producer="running", complete=False,
                  slurm_job_id=os.environ["SLURM_JOB_ID"], attempts=[])
    path = run / "workers/status.json"
    write(path, record)
    deadline, qualified = time.monotonic()+manifest["soft_seconds"], False
    try:
        for duration in manifest["preparation_durations"]:
            directory = run / "workers" / f"attempt_{duration:g}"
            directory.mkdir()
            attempt = dict(duration=duration, status="running", phases=[], basal_elapsed=0.0)
            record["attempts"].append(attempt)
            def phase(kind, state, peak, end, interval):
                record.update(generated_at=now(), status=kind, current_duration=duration)
                write(path, record)
                result = integrate(closure, state, peak, end, interval, manifest, deadline)
                archive = directory / f'{len(attempt["phases"]):03d}_{kind}.npz'
                np.savez_compressed(archive, time=result["time"], states=result["states"],
                    material_positions=closure.material_positions, cue_peak=np.asarray(peak),
                    total_physical_cue_support=np.asarray(manifest["total_physical_support"]))
                metadata = {key:value for key,value in result.items() if key not in ("time","states")}
                metadata.update(path=str(archive.relative_to(run)), sha256=sha(archive), kind=kind)
                attempt["phases"].append(metadata)
                write(path, record)
                return result
            prepared = phase("preparation", initial.copy(), manifest["preparation_peak"], duration, min(.1,duration/100))
            if prepared["status"] != "complete":
                attempt["status"] = prepared["status"]
                break
            state, consecutive = prepared["states"][-1].copy(), 0
            while attempt["basal_elapsed"] < manifest["rule"]["max_basal_per_attempt"]:
                chunk = phase("basal_window", state, 0.0, 1000.0, 5.0)
                state = chunk["states"][-1].copy()
                attempt["basal_elapsed"] += float(chunk["time"][-1])
                if chunk["status"] != "complete":
                    attempt["status"] = chunk["status"]
                    break
                steady = diagnostics(chunk["time"],chunk["states"],manifest)
                attempt["phases"][-1]["diagnostics"] = steady
                consecutive = consecutive+1 if steady["steady"] else 0
                if consecutive < 2:
                    continue
                if manifest["rule"]["max_basal_per_attempt"]-attempt["basal_elapsed"] < 10000:
                    break
                target = state.copy()
                hold = phase("holdout",state,0.0,10000.0,5.0)
                state = hold["states"][-1].copy()
                attempt["basal_elapsed"] += float(hold["time"][-1])
                if hold["status"] != "complete":
                    attempt["status"] = hold["status"]
                    break
                scales = np.r_[np.full(2*closure.cells,2.),closure.resting_circumference]
                excursion = float(np.max(np.abs(hold["states"]-target[None])/scales[None]))
                final_steady = diagnostics(hold["time"],hold["states"],manifest)
                shape = architecture(state,closure.cells)
                hold_check = dict(excursion=excursion,final_steady=final_steady,architecture=shape)
                attempt["phases"][-1]["diagnostics"] = hold_check
                if excursion <= .001 and final_steady["steady"]:
                    attempt["verified_settled"] = True
                    if shape["contrast"] < .05:
                        attempt["status"] = "homogeneous_or_low_contrast"
                    elif shape["connected_regions"] != 1 or shape["resolved_peaks"] != 1:
                        attempt["status"] = "multisite_or_ambiguous"
                    else:
                        qualified, attempt["status"] = True,"qualifying_single_site"
                        checkpoint = run / "workers/qualified_basal.npz"
                        np.savez_compressed(checkpoint,state=state,
                            material_positions=closure.material_positions,
                            resting_circumference=np.asarray(closure.resting_circumference))
                        record.update(checkpoint=str(checkpoint.relative_to(run)),checkpoint_sha256=sha(checkpoint))
                        write(checkpoint.with_suffix(".json"),dict(generated_at=now(),
                            origin="Fresh periodic uniform start; cue-selected incumbent",
                            model=manifest["model"],checkpoint_sha256=sha(checkpoint),
                            preparation_duration=duration,basal_elapsed=attempt["basal_elapsed"],
                            holdout=hold_check, source_sha256=manifest["source_sha256"],
                            native_2d_equivalence=False, challenger_trials=0))
                    break
                consecutive = 0
            if attempt["status"] == "running":
                attempt["status"] = "nonsettled"
            attempt["final_architecture"] = architecture(state,closure.cells)
            attempt["final_amount_relative_error"] = closure.amount(state)/(2*closure.resting_circumference)-1
            write(path,record)
            if qualified or attempt["status"] in ("timeout","numerical_error"):
                break
        record.update(generated_at=now(),producer="finished",complete=True,
                      status="qualifying_single_site" if qualified else record["attempts"][-1]["status"])
    except Exception as error:
        record.update(generated_at=now(),producer="finished",complete=True,status="error",
                      error=f"{type(error).__name__}: {error}")
        write(path,record)
        raise
    write(path,record)
    print(json.dumps(record,indent=2))


def audit(run):
    manifest,record = read(run / "manifest.json"),read(run / "workers/status.json")
    assert record["complete"]
    assert sha(run / "manifest.json") == read(run / "launch.json")["manifest_sha256"]
    for relative,digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    assert sha(run / "inputs/uniform.npz") == manifest["input_sha256"]
    closure = Closure(manifest["resting_circumference"],manifest["cells"])
    with np.load(run / "inputs/uniform.npz") as z:
        initial = z["state"].copy()
    np.testing.assert_array_equal(initial,closure.initial())
    worst,phase_count,hold_check = 0.,0,None
    for attempt in record["attempts"]:
        previous,windows = initial.copy(),[]
        for phase in attempt["phases"]:
            assert sha(run / phase["path"]) == phase["sha256"]
            with np.load(run / phase["path"]) as z:
                times,states,peak = z["time"].copy(),z["states"].copy(),float(z["cue_peak"])
            np.testing.assert_array_equal(states[0],previous)
            assert np.isfinite(states).all() and np.all(np.diff(times)>0)
            assert np.min(states[:,:-1]) >= -manifest["rule"]["allowed_negative_roundoff"]
            amounts = states[:,-1]*np.mean(states[:,:closure.cells]+states[:,closure.cells:-1],axis=1)
            worst = max(worst,float(np.max(np.abs(amounts/(2*closure.resting_circumference)-1))))
            assert worst <= manifest["rule"]["pool_relative_tolerance"]
            if phase["kind"] == "preparation":
                assert peak == 4 and times[-1] == attempt["duration"]
            else:
                assert peak == 0
                if phase["kind"] == "basal_window":
                    recomputed = diagnostics(times,states,manifest)
                    assert recomputed == phase["diagnostics"]
                    windows.append(recomputed["steady"])
                else:
                    assert windows[-2:] == [True,True]
                    scales = np.r_[np.full(2*closure.cells,2.),closure.resting_circumference]
                    excursion = float(np.max(np.abs(states-previous[None])/scales[None]))
                    hold_check = dict(excursion=excursion,final_steady=diagnostics(times,states,manifest),
                                      architecture=architecture(states[-1],closure.cells))
                    assert hold_check == phase["diagnostics"]
                    if attempt.get("verified_settled"):
                        assert excursion <= .001 and hold_check["final_steady"]["steady"]
            previous = states[-1].copy()
            phase_count += 1
        if attempt["status"] == "qualifying_single_site":
            assert hold_check["architecture"]["contrast"] >= .05
            assert hold_check["architecture"]["connected_regions"] == 1
            assert hold_check["architecture"]["resolved_peaks"] == 1
    if "checkpoint" in record:
        assert sha(run / record["checkpoint"]) == record["checkpoint_sha256"]
        with np.load(run / record["checkpoint"]) as z:
            np.testing.assert_array_equal(z["state"],previous)
    stationary_derivative = closure.rhs(previous,np.zeros(closure.cells))
    assert np.isfinite(stationary_derivative).all()
    derivative_scales = np.r_[np.full(2*closure.cells,2.),closure.resting_circumference]
    projected_rhs = float(np.max(np.abs(stationary_derivative)/derivative_scales)*1000)
    if record["status"] == "qualifying_single_site":
        assert projected_rhs < 1e-4
    stderr_path=run / "logs" / (read(run / "launch.json")["slurm_job_id"]+".err")
    stderr=stderr_path.read_text() if stderr_path.exists() else ""
    result = dict(generated_at=now(),passed=True,status=record["status"],
        phase_archives=phase_count,max_sampled_amount_relative_error=worst,holdout=hold_check,
        final_native_rhs_scaled_to_1000=projected_rhs,
        stderr=stderr,
        stderr_interpretation="A SciPy BDF internal-difference RuntimeWarning is retained verbatim; all stored states and final derivatives are finite and qualification checks pass. Its internal cause is not established by this audit." if stderr else "",
        scope="Saved-state identities, full-field/C qualification and physical amount; no mesh convergence or independent replay")
    write(run / "results/audit.json",result)
    print(json.dumps(result,indent=2))


def collect(run,figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    manifest,record = read(run / "manifest.json"),read(run / "workers/status.json")
    checked = read(run / "results/audit.json")
    assert checked["passed"]
    n = manifest["cells"]
    attempt = record["attempts"][-1]
    times,states,offset = [],[],0.
    sources = {"manifest.json":sha(run / "manifest.json"),
               "workers/status.json":sha(run / "workers/status.json"),
               "results/audit.json":sha(run / "results/audit.json")}
    for phase in attempt["phases"]:
        sources[phase["path"]] = sha(run / phase["path"])
        assert sources[phase["path"]] == phase["sha256"]
        with np.load(run / phase["path"]) as z:
            t,y = z["time"].copy(),z["states"].copy()
        take = slice(None) if not times else slice(1,None)
        times.extend(t[take]+offset);states.extend(y[take]);offset+=t[-1]
    times,states = np.asarray(times),np.asarray(states)
    q = (np.arange(n)+.5)/n
    with np.load(run / attempt["phases"][0]["path"]) as z:
        withdrawn = z["states"][-1].copy()
    display_shift = n//2-int(np.argmax(states[-1,:n]))
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,
        "axes.titlesize":9,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
        "pdf.fonttype":42,"svg.fonttype":"none"})
    figure,axes = plt.subplots(2,2,figsize=(7.1,5.0))
    figure.subplots_adjust(left=.10,right=.94,bottom=.12,top=.84,wspace=.44,hspace=.58)
    colors = ["#888888","#c29c51","#235789"]
    for state,label,color in zip([states[0],withdrawn,states[-1]],
        ["Uniform start","Cue withdrawal","Final basal"],colors):
        axes[0,0].plot(q,np.roll(state[:n],display_shift),label=label,color=color,lw=1.4)
    axes[0,0].set(title="A   Single active site",xlabel="Material position around ring",ylabel="Active Rac")
    axes[0,0].legend(frameon=False,fontsize=7)
    basal = times>=attempt["duration"]
    heat = axes[1,0].pcolormesh(q,times[basal]-attempt["duration"],
        np.roll(states[basal,:n],display_shift,axis=1),shading="auto",cmap="viridis",rasterized=True)
    figure.colorbar(heat,ax=axes[1,0],fraction=.045,pad=.03)
    axes[1,0].set(title="C   Basal continuation",xlabel="Material position around ring",ylabel="Time since cue withdrawal")
    axes[0,1].plot(times,states[:,-1]-manifest["resting_circumference"],color="#235789")
    axes[0,1].axvline(attempt["duration"],color="#c29c51",ls=":")
    axes[0,1].set(title="B   Extension remains dynamic",xlabel="Time since preparation onset",ylabel="C − C₀")
    right = axes[0,1].secondary_yaxis("right",functions=(lambda x:3+80*x,lambda x:(x-3)/80))
    right.set_ylabel("Global deactivation coefficient")
    axes[0,1].set_xscale("symlog",linthresh=10)
    axes[0,1].set_xlim(0,times[-1])
    amount = states[:,-1]*np.mean(states[:,:n]+states[:,n:-1],axis=1)
    axes[1,1].plot(times,(amount/(2*manifest["resting_circumference"])-1)*1e9,color="#235789")
    axes[1,1].axvline(attempt["duration"],color="#c29c51",ls=":")
    axes[1,1].set(title="D   Conserved physical amount",xlabel="Time since preparation onset",ylabel="Relative amount error (×10⁻⁹)")
    axes[1,1].set_xscale("symlog",linthresh=10)
    axes[1,1].set_xlim(0,times[-1])
    figure.suptitle("Spring: circular global-extension analogue\nSingle active site retained after cue withdrawal" if record["status"]=="qualifying_single_site" else "Spring: circular global-extension analogue\n"+record["status"].replace("_"," "),fontsize=11,y=.98)
    files=[]
    for extension in ("png","pdf","svg"):
        path=dev_output_path(figdir,"spring_ring_preparation_feasibility."+extension)
        figure.savefig(path,dpi=250,facecolor="white",bbox_inches="tight");files.append(path)
    plt.close(figure)
    caption=("Spring circular global-extension analogue, preparation feasibility only. New: periodic signalling, "
        "one uniformly deforming circumference, and Cdot=2(mean F(u)−0.01(C−C0)), averaging after the native nonlinear force switch. "
        "Native local reaction/diffusion/force functions, extension-dependent deactivation coefficients and dilution are retained. "
        "Fresh uniform u=.05645, v=2−u, C0=2; no native interval cosine or imported endpoint checkpoint. "
        f"A centred raised-cosine cue selects the site: peak 4, total instantaneous physical support 0.2, duration {attempt['duration']:g}. "
        "A: active profiles, cyclic shift for display only; position 0 joins 1. B: physical extension and corresponding global deactivation. "
        "C: complete active-field basal history. D: relative error in physical amount C×mean(u+v); no clipping or renormalization. "
        f"Basal continuation {attempt['basal_elapsed']:g} native units; qualification tests full active/inactive fields AND circumference "
        "on frozen scales 2, 2, C0, requiring two 1000-unit stable windows and a 10000-unit holdout. "
        "Full dense numerical Jacobian includes global mean-force coupling. This is not the native 2D continuum mechanics, "
        "a mesh-convergence study or a challenger-distance experiment. Generated "+now()+"; producer "+record["producer"]+".")
    caption_path=dev_output_path(figdir,"caption.txt");caption_path.write_text(caption+"\n");files.append(caption_path)
    summary=dict(generated_at=now(),status=record["status"],complete=record["complete"],
        attempts=[dict(duration=a["duration"],status=a["status"],basal_elapsed=a["basal_elapsed"]) for a in record["attempts"]],
        final_architecture=attempt.get("final_architecture"),max_amount_error=checked["max_sampled_amount_relative_error"],
        checkpoint=record.get("checkpoint"),display_only_cyclic_shift_cells=display_shift,
        native_2d_equivalence=False,challenger_trials=0)
    summary_path=dev_output_path(figdir,"summary.json");write(summary_path,summary);files.append(summary_path)
    write(dev_output_path(figdir,"manifest.json"),dict(generated_at=now(),artifact_root="../..",run=str(run),
        source_sha256=sources,sha256={dev_artifact_key(figdir,path):sha(path) for path in files},
        renderer_sha256=sha(Path(__file__)),canonical_models_changed=False))
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("action",choices=["launch","worker","audit","collect"])
    parser.add_argument("--run",type=Path,default=DEFAULT)
    parser.add_argument("--fig",type=Path,default=FIG)
    args=parser.parse_args()
    if args.action=="launch":launch(args.run.resolve())
    elif args.action=="worker":worker(args.run.resolve())
    elif args.action=="audit":audit(args.run.resolve())
    else:collect(args.run.resolve(),args.fig.resolve())
