"""Bounded SLURM-only preparation feasibility for three development rings."""
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
import scipy
from scipy.integrate import BDF
import ring_distance_prototype_20261017 as ring
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/ring_preparation_pilot_20261017"
FIG = ROOT / "figures/length_mechanisms_20261008/dev/ring_preparation_pilot_20261017"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
SCALES = {"wave_pinning": [2.2683, 2.2683], "goryachev": [2.75033, 2.75033],
          "holmes_model3": [7.5, 7.5, 3.1, 3.1]}
PEAKS = {"wave_pinning": 0.1, "goryachev": 0.1, "holmes_model3": 1.0}
STOP = False


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def init(run, cells=400):
    if run.exists():
        raise ValueError("Choose a fresh pilot directory; existing output is preserved")
    if "data/runs" not in str(run):
        raise ValueError("Development computation must be under data/runs")
    for name in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run / name).mkdir(parents=True)
    sources = [Path(__file__), Path(__file__).with_name("ring_distance_prototype_20261017.py"),
               Path(__file__).with_name("length_dev_output_layout.py")]
    for directory in (ROOT / "code/core/models", ROOT / "code/core/src/polarity1d"):
        sources += sorted(directory.glob("*.py"))
    hashes = {}
    for path in sources:
        relative = path.relative_to(ROOT)
        destination = run / "frozen_sources" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        hashes[str(relative)] = sha(destination)
    grid = ring.RingGrid(40, cells)
    jobs = []
    for slot, key in enumerate(ring.KEYS):
        model = ring.native_model(key)
        state = ring.native_rest(key, grid)
        assert np.all(np.ptp(state, axis=1) == 0)
        path = run / "inputs" / f"{key}.npz"
        np.savez_compressed(path, state=state, x=grid.x)
        means = ring.pool_amounts(model, grid, state) / grid.circumference
        expected = [SCALES[key][group[0]] for group in model.conserved_groups]
        np.testing.assert_allclose(means, expected, atol=1e-14)
        jobs.append(dict(slot=slot, model=key, input=str(path.relative_to(run)), input_sha256=sha(path),
                         state_names=list(model.state_names), diffusivities=list(model.diffusivities),
                         conserved_groups=[list(group) for group in model.conserved_groups],
                         fixed_field_scales=dict(zip(model.state_names, SCALES[key])),
                         prep_peak=PEAKS[key], native_pool_means=means.tolist()))
    manifest = dict(created_at=now(), status="Frozen before SLURM launch", jobs=jobs,
        circumference=40, cells=cells, prep_total_support=3, prep_center_cell=0,
        preparation_durations=[10.0, 100.0, 1000.0], basal_external_input=0.0,
        rule=dict(chunk=1000.0, tolerance=1e-4, consecutive_windows=2, holdout=10000.0,
                  holdout_tolerance=1e-3, max_total_basal_time_per_attempt=30000.0,
                  minimum_primary_contrast_native=0.05, minimum_resolved_peak_prominence_native=0.05,
                  allowed_negative_roundoff=5e-10, pool_relative_tolerance=1e-7),
        scales_definition="Each field uses its native conserved-group homogeneous mean, frozen before launch; no future-endpoint normalization.",
        single_site_definition="Verified settled basal state; primary contrast >=0.05 native (existing manuscript eligibility floor); exactly one periodic half-contrast region and one peak with prominence >=0.05 native. All peaks/prominences and translations remain recorded. These are pilot setup diagnostics, not challenger gates.",
        solver=dict(method="BDF", rtol=1e-9, atol=1e-11, max_step=10.0,
                    preparation_sampling="min(0.5,duration/100)", basal_sampling=5.0),
        time_budget=dict(slurm_seconds=3600, worker_soft_seconds=3300),
        native_start="Fresh native homogeneous state; Wave-pinning rounded native rest retained. No endpoint checkpoint imported.",
        cue_protocol_comparison="Ring preparation is a symmetric centred raised cosine of total width 3, first duration 10, peaks WP 0.1/Goryachev 0.1/Holmes 1. The current interval three-test response instead uses endpoint patches (relative 20% for the illustrative condition), including Goryachev K64-scaled peak 6.4. Native local reactions/diffusivities and conserved-pool means are unchanged; ring and interval preparation/response doses are not asserted equal.",
        source_sha256=hashes, integration_scope="Preparation only; no challenger runs or duration optimization.",
        environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, executable=PYTHON))
    write(run / "manifest.json", manifest)
    write(run / "status.json", dict(generated_at=now(), producer="not started", complete=False))
    return manifest


def launch(run, cells=400):
    manifest = init(run, cells)
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = "exec " + shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", PYTHON, str(runner), "worker", "--run", str(run)])
    command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G",
            "--time=01:00:00", "--signal=B:USR1@60", "--array=0-2", "--job-name=ring-prep",
            "--output=" + str(run / "logs/%A_%a.out"), "--error=" + str(run / "logs/%A_%a.err"),
            "--wrap", command]
    job_id = subprocess.check_output(args, text=True).strip()
    write(run / "launch.json", dict(generated_at=now(), array_job_id=job_id,
                                   command=args, manifest_sha256=sha(run / "manifest.json")))
    print(json.dumps(dict(run=str(run), array_job_id=job_id, jobs=len(manifest["jobs"])), indent=2))


def samples(end, interval):
    result = np.arange(0, end + interval / 2, interval)
    result = result[result <= end]
    return np.r_[result, end] if result[-1] != end else result


def integrate(model, grid, initial, cue, end, interval, manifest, deadline, initial_pools):
    """BDF accepted-step loop preserving partial full-field output on timeout."""
    species = len(model.state_names)
    tolerances = manifest["solver"]
    solver = BDF(lambda t, y: ring.native_rhs(model, grid, y.reshape(species, grid.cells), cue).ravel(),
        0, initial.ravel(), end, rtol=tolerances["rtol"], atol=tolerances["atol"],
        max_step=tolerances["max_step"], jac_sparsity=ring.periodic_sparsity(species, grid.cells))
    wanted = samples(end, interval)
    times, states = [0.0], [initial.copy()]
    index, status, message = 1, "complete", ""
    worst_pool, minimum = 0.0, float(np.min(initial))
    while solver.status == "running":
        if STOP or time.monotonic() >= deadline:
            status, message = "timeout", "Worker soft deadline or SLURM advance signal"
            break
        message = solver.step() or ""
        if solver.status == "failed":
            status = "numerical_error"
            break
        state = solver.y.reshape(species, grid.cells)
        minimum = min(minimum, float(np.min(state)))
        worst_pool = max(worst_pool, float(np.max(np.abs(ring.pool_amounts(model, grid, state) / initial_pools - 1))))
        dense = solver.dense_output()
        while index < len(wanted) and wanted[index] <= solver.t:
            times.append(float(wanted[index]))
            states.append(dense(wanted[index]).reshape(species, grid.cells).copy())
            index += 1
        if minimum < -manifest["rule"]["allowed_negative_roundoff"] or worst_pool > manifest["rule"]["pool_relative_tolerance"]:
            status, message = "numerical_error", "Negative concentration or conserved-pool drift"
            break
    if times[-1] != solver.t:
        times.append(float(solver.t))
        states.append(solver.y.reshape(species, grid.cells).copy())
    return dict(time=np.asarray(times), states=np.asarray(states), status=status, message=message,
                max_pool_relative_error=worst_pool, minimum_concentration=minimum,
                nfev=int(solver.nfev), njev=int(solver.njev), nlu=int(solver.nlu))


def diagnostics(time_values, states, scales, rule):
    keep = time_values >= time_values[-1] - rule["chunk"]
    values, time_values = states[keep], time_values[keep]
    variation = np.max(np.ptp(values, axis=0), axis=1) / scales
    drift = np.max(np.abs(values[-1] - values[-2]), axis=1) / np.diff(time_values[-2:])[0] * rule["chunk"] / scales
    return dict(steady=bool(np.max(np.r_[variation, drift]) <= rule["tolerance"]),
                relative_range=variation.tolist(), projected_relative_drift=drift.tolist())


def architecture(grid, state, basal_reference=None):
    primary = state[0]
    contrast = float(np.ptp(primary))
    if contrast < 1e-12:
        return dict(contrast=contrast, connected_regions=0, resolved_peaks=0, peak_cells=[], prominence=[], largest_peak_displacement=0)
    result = ring.spatial_readout(grid, primary, primary if basal_reference is None else basal_reference[0])
    result["contrast"] = contrast
    result["resolved_peaks"] = int(np.sum(np.asarray(result["prominence"]) >= 0.05))
    return result


def worker(run, slot):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("All pilot integration must run inside the reviewed SLURM job")
    global STOP
    signal.signal(signal.SIGUSR1, lambda signum, frame: globals().__setitem__("STOP", True))
    signal.signal(signal.SIGTERM, lambda signum, frame: globals().__setitem__("STOP", True))
    manifest = read(run / "manifest.json")
    for path, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / path) == digest, path
    job = manifest["jobs"][slot]
    assert sha(run / job["input"]) == job["input_sha256"]
    directory = run / "workers" / f"{slot:02d}"
    directory.mkdir()
    record = dict(generated_at=now(), producer="running", status="starting", job=job,
                  slurm_job_id=os.environ["SLURM_JOB_ID"], attempts=[])
    status_path = directory / "status.json"
    write(status_path, record)
    grid, model = ring.RingGrid(manifest["circumference"], manifest["cells"]), ring.native_model(job["model"])
    with np.load(run / job["input"]) as z:
        homogeneous = z["state"].copy()
    scales = np.array([job["fixed_field_scales"][name] for name in job["state_names"]])
    pools = ring.pool_amounts(model, grid, homogeneous)
    cue = ring.ring_cue(grid, 0, manifest["prep_total_support"], job["prep_peak"])
    zero = np.zeros(grid.cells)
    deadline = time.monotonic() + manifest["time_budget"]["worker_soft_seconds"]
    try:
        qualified = False
        for duration in manifest["preparation_durations"]:
            attempt_dir = directory / f"attempt_{duration:g}"
            attempt_dir.mkdir()
            attempt = dict(duration=duration, status="running", phases=[], basal_elapsed=0)
            record["attempts"].append(attempt)
            phase_index = 0

            def phase(kind, state, signal_values, end, interval):
                nonlocal phase_index
                record.update(generated_at=now(), status=kind, current_attempt_duration=duration)
                write(status_path, record)
                result = integrate(model, grid, state, signal_values, end, interval, manifest, deadline, pools)
                path = attempt_dir / f"{phase_index:03d}_{kind}.npz"
                np.savez_compressed(path, time=result["time"], states=result["states"], x=grid.x,
                                    state_names=np.asarray(job["state_names"]), stimulus=signal_values)
                metadata = {key: value for key, value in result.items() if key not in ("time", "states")}
                metadata.update(path=str(path.relative_to(run)), sha256=sha(path), phase=kind,
                                phase_index=phase_index, stored_samples=len(result["time"]))
                attempt["phases"].append(metadata)
                phase_index += 1
                write(status_path, record)
                return result

            prepared = phase("preparation", homogeneous.copy(), cue, duration, min(0.5, duration / 100))
            if prepared["status"] != "complete":
                attempt["status"] = prepared["status"]
                break
            state, consecutive = prepared["states"][-1].copy(), 0
            previous = state.copy()
            while attempt["basal_elapsed"] < manifest["rule"]["max_total_basal_time_per_attempt"]:
                chunk = phase("basal_window", state, zero, 1000, 5)
                state = chunk["states"][-1].copy()
                attempt["basal_elapsed"] += float(chunk["time"][-1])
                if chunk["status"] != "complete":
                    attempt["status"] = chunk["status"]
                    break
                metric = diagnostics(chunk["time"], chunk["states"], scales, manifest["rule"])
                metric["translation"] = ring.translation_compare(state, previous, scales)
                metric["architecture"] = architecture(grid, state)
                attempt["phases"][-1]["diagnostics"] = metric
                previous = state.copy()
                consecutive = consecutive + 1 if metric["steady"] else 0
                if consecutive < 2:
                    continue
                remaining = manifest["rule"]["max_total_basal_time_per_attempt"] - attempt["basal_elapsed"]
                if remaining < 10000:
                    break
                target = state.copy()
                hold = phase("holdout", state, zero, 10000, 5)
                state = hold["states"][-1].copy()
                attempt["basal_elapsed"] += float(hold["time"][-1])
                if hold["status"] != "complete":
                    attempt["status"] = hold["status"]
                    break
                excursion = float(np.max(np.abs(hold["states"] - target[None]) / scales[None, :, None]))
                final_steady = diagnostics(hold["time"], hold["states"], scales, manifest["rule"])
                shape = architecture(grid, state, target)
                hold_check = dict(excursion=excursion, final_steady=final_steady,
                                  translation=ring.translation_compare(state, target, scales), architecture=shape)
                attempt["phases"][-1]["diagnostics"] = hold_check
                if excursion <= 0.001 and final_steady["steady"]:
                    attempt["verified_settled"] = True
                    if shape["contrast"] < 0.05:
                        attempt["status"] = "homogeneous_or_low_contrast"
                    elif shape["connected_regions"] != 1 or shape["resolved_peaks"] != 1:
                        attempt["status"] = "multisite_or_ambiguous"
                    else:
                        attempt["status"] = "qualifying_single_site"
                        qualified = True
                        checkpoint = directory / "qualified_basal.npz"
                        np.savez_compressed(checkpoint, state=state, x=grid.x,
                            model=np.asarray(job["model"]), state_names=np.asarray(job["state_names"]),
                            circumference=np.asarray(grid.circumference), cells=np.asarray(grid.cells))
                        write(checkpoint.with_suffix(".json"), dict(generated_at=now(),
                            origin="fresh_periodic_preparation", boundary_condition="periodic",
                            checkpoint_sha256=sha(checkpoint), native_input_sha256=job["input_sha256"],
                            preparation_duration=duration, basal_elapsed=attempt["basal_elapsed"],
                            fixed_field_scales=job["fixed_field_scales"], holdout=hold_check,
                            source_sha256=manifest["source_sha256"],
                            scope="Preparation feasibility only; mesh adequacy and challenger outcomes remain unestablished."))
                        record["checkpoint"] = str(checkpoint.relative_to(run))
                        record["checkpoint_sha256"] = sha(checkpoint)
                    break
                consecutive, previous = 0, state.copy()
            if attempt["status"] == "running":
                motion = ring.translation_compare(state, prepared["states"][-1], scales)
                attempt["status"] = "drifting_or_nonsettled" if motion["best_integer_shift_cells"] else "nonsettled"
                attempt["translation_from_withdrawal"] = motion
            attempt["final_architecture"] = architecture(grid, state)
            attempt["final_pool_relative_error"] = (ring.pool_amounts(model, grid, state) / pools - 1).tolist()
            write(status_path, record)
            if qualified or attempt["status"] in ("timeout", "numerical_error"):
                break
        record.update(generated_at=now(), producer="finished", complete=True,
            status="qualifying_single_site" if qualified else record["attempts"][-1]["status"])
    except Exception as error:
        record.update(generated_at=now(), producer="finished", complete=True, status="error",
                      error=f"{type(error).__name__}: {error}")
        write(status_path, record)
        raise
    write(status_path, record)
    print(json.dumps(dict(model=job["model"], status=record["status"], attempts=len(record["attempts"])), indent=2))


def status(run):
    manifest = read(run / "manifest.json")
    records = []
    for job in manifest["jobs"]:
        path = run / "workers" / f'{job["slot"]:02d}' / "status.json"
        records.append(read(path) if path.exists() else dict(job=job, status="not_started", complete=False))
    result = dict(generated_at=now(), producer="finished" if all(row.get("complete") for row in records) else "running",
                  complete=all(row.get("complete") for row in records), workers=records)
    write(run / "status.json", result)
    return result


def audit(run):
    """Bounded independent read-only reconstruction of saved qualification."""
    manifest = read(run / "manifest.json")
    assert sha(run / "manifest.json") == read(run / "launch.json")["manifest_sha256"]
    for path, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / path) == digest, path
    records = []
    for job in manifest["jobs"]:
        assert sha(run / job["input"]) == job["input_sha256"]
        with np.load(run / job["input"]) as z:
            initial = z["state"].copy()
        assert np.all(np.ptp(initial, axis=1) == 0)
        directory = run / "workers" / f'{job["slot"]:02d}'
        record = read(directory / "status.json")
        assert record["complete"]
        scales = np.asarray([job["fixed_field_scales"][name] for name in job["state_names"]])
        groups = job["conserved_groups"]
        dx = manifest["circumference"] / manifest["cells"]
        initial_amounts = np.array([initial[group].sum() * dx for group in groups])
        maximum_pool_error, maximum_excursion = 0.0, 0.0
        windows, final_state, phase_count = [], None, 0
        for attempt in record["attempts"]:
            previous = initial.copy()
            for phase in attempt["phases"]:
                path = run / phase["path"]
                assert sha(path) == phase["sha256"]
                with np.load(path) as z:
                    times, states, stimulus = z["time"].copy(), z["states"].copy(), z["stimulus"].copy()
                assert states.shape[1:] == initial.shape and np.isfinite(states).all()
                np.testing.assert_array_equal(states[0], previous)
                assert np.all(np.diff(times) > 0)
                assert np.min(states) >= -manifest["rule"]["allowed_negative_roundoff"]
                amounts = np.column_stack([states[:, group, :].sum(axis=(1, 2)) * dx for group in groups])
                maximum_pool_error = max(maximum_pool_error, float(np.max(np.abs(amounts / initial_amounts - 1))))
                assert maximum_pool_error <= manifest["rule"]["pool_relative_tolerance"]
                if phase["phase"] == "preparation":
                    assert times[-1] == attempt["duration"] and stimulus.max() == job["prep_peak"]
                else:
                    np.testing.assert_array_equal(stimulus, np.zeros(manifest["cells"]))
                    if phase["phase"] == "basal_window":
                        ranges = np.ptp(states, axis=0).max(axis=1) / scales
                        drift = np.abs(states[-1] - states[-2]).max(axis=1) / (times[-1]-times[-2]) * 1000 / scales
                        windows.append(bool(np.max(np.r_[ranges, drift]) <= 1e-4))
                    elif phase["phase"] == "holdout":
                        assert windows[-2:] == [True, True]
                        maximum_excursion = max(maximum_excursion, float(np.max(np.abs(states - previous[None]) / scales[None,:,None])))
                        keep = times >= times[-1] - 1000
                        ranges = np.ptp(states[keep], axis=0).max(axis=1) / scales
                        drift = np.abs(states[-1] - states[-2]).max(axis=1) / (times[-1]-times[-2]) * 1000 / scales
                        assert maximum_excursion <= 1e-3 and np.max(np.r_[ranges, drift]) <= 1e-4
                previous, final_state = states[-1].copy(), states[-1].copy()
                phase_count += 1
        if record["status"] == "qualifying_single_site":
            path = run / record["checkpoint"]
            assert sha(path) == record["checkpoint_sha256"]
            with np.load(path) as z:
                np.testing.assert_array_equal(z["state"], final_state)
            assert read(path.with_suffix(".json"))["origin"] == "fresh_periodic_preparation"
            primary = final_state[0]
            contrast = float(np.ptp(primary))
            mask = primary >= primary.min() + 0.5 * contrast
            components = int(np.sum(mask & ~np.roll(mask, 1)))
            assert contrast >= 0.05 and components == 1
        records.append(dict(model=job["model"], status=record["status"], full_field_phase_archives=phase_count,
            native_homogeneous_start_verified=True, exact_phase_and_checkpoint_continuity=True,
            max_sampled_pool_relative_error=maximum_pool_error, max_scaled_holdout_excursion=maximum_excursion,
            all_qualification_criteria_recomputed=True))
    result = dict(generated_at=now(), passed=True, models=records,
        scope="Saved full-field preparation feasibility only; no independent trajectory replay, mesh claim or challenger outcome.")
    write(run / "results/independent_audit.json", result)
    print(json.dumps(result, indent=2))
    return result


def collect(run, figdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    current = status(run)
    manifest = read(run / "manifest.json")
    summaries = []
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":8, "axes.titlesize":9,
        "axes.labelsize":8, "xtick.labelsize":7, "ytick.labelsize":7, "pdf.fonttype":42, "svg.fonttype":"none"})
    figure, axes = plt.subplots(2, 3, figsize=(8.2, 5.4))
    figure.subplots_adjust(left=.08, right=.985, top=.85, bottom=.11, wspace=.42, hspace=.48)
    labels = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "holmes_model3":"Holmes 3"}
    source_hashes = {"manifest.json":sha(run / "manifest.json")}
    audit_path = run / "results/independent_audit.json"
    if audit_path.exists():
        assert read(audit_path)["passed"]
        source_hashes["results/independent_audit.json"] = sha(audit_path)
    for column, record in enumerate(current["workers"]):
        job = record["job"]
        status_path = run / "workers" / f'{job["slot"]:02d}' / "status.json"
        if status_path.exists():
            source_hashes[str(status_path.relative_to(run))] = sha(status_path)
        if not record.get("attempts"):
            axes[0,column].set_title(labels[job["model"]] + "\nnot started")
            continue
        attempt = record["attempts"][-1]
        phases = attempt["phases"]
        prep = next((item for item in phases if item["phase"]=="preparation"), None)
        if prep is None:
            continue
        with np.load(run / prep["path"]) as z:
            x, initial, withdrawn = z["x"].copy(), z["states"][0].copy(), z["states"][-1].copy()
        times, fields, offset = [], [], 0.0
        for item in phases:
            source_hashes[item["path"]] = sha(run / item["path"])
            assert source_hashes[item["path"]] == item["sha256"]
            if item["phase"] == "preparation":
                continue
            with np.load(run / item["path"]) as z:
                t, values = z["time"].copy(), z["states"][:,0].copy()
            times.extend((t if not times else t[1:]) + offset)
            fields.extend(values if not fields else values[1:])
            offset += float(t[-1])
        final = np.asarray(fields)[-1] if fields else withdrawn[0]
        display_shift = len(x)//2 - int(np.argmax(final))
        axes[0,column].plot(x, np.roll(initial[0],display_shift), color="#777777", ls=":", label="Native start")
        axes[0,column].plot(x, np.roll(withdrawn[0],display_shift), color="#c29c51", lw=1.2, label="Cue withdrawal")
        axes[0,column].plot(x, np.roll(final,display_shift), color="#235789", lw=1.5, label="Latest basal")
        axes[0,column].set(title=labels[job["model"]] + "\n" + record["status"].replace("_"," "),
                           xlabel="Ring position (native length)", xlim=(0,40), xticks=[0,20,40])
        if column==0:
            axes[0,column].set_ylabel("Primary activity\n(native concentration)")
        if fields:
            mesh=axes[1,column].pcolormesh(x, np.asarray(times), np.roll(np.asarray(fields),display_shift,axis=1), shading="auto", cmap="viridis", rasterized=True)
            figure.colorbar(mesh, ax=axes[1,column], fraction=.05, pad=.02)
        axes[1,column].set(xlabel="Ring position (native length)", xlim=(0,40), xticks=[0,20,40])
        if column==0:
            axes[1,column].set_ylabel("Time after cue withdrawal\n(model-native time)")
        summaries.append(dict(model=job["model"], status=record["status"], preparation_duration=attempt["duration"],
            basal_elapsed=attempt["basal_elapsed"], attempts=[dict(duration=item["duration"],status=item["status"]) for item in record["attempts"]],
            checkpoint=record.get("checkpoint"), display_only_cyclic_shift_cells=display_shift,
            final_architecture=attempt.get("final_architecture")))
    handles, labels_legend = axes[0,0].get_legend_handles_labels()
    figure.legend(handles,labels_legend,loc="upper center",bbox_to_anchor=(.52,.995),ncol=3,frameon=False,fontsize=8)
    files=[]
    for extension in ("png","pdf","svg"):
        path=dev_output_path(figdir, "ring_preparation_feasibility."+extension)
        figure.savefig(path,dpi=250,facecolor="white",bbox_inches="tight")
        files.append(path)
    plt.close(figure)
    caption=("Preparation feasibility only. Fresh native homogeneous starts and a symmetric raised-cosine preparation cue of total support 3, duration 10, and native peaks Wave-pinning 0.1/Goryachev 0.1/Holmes 1, followed by zero external basal stimulus on C = 40, N = "+str(manifest["cells"])+" periodic rings. Native reactions, diffusivities and pool means are unchanged; this is not dose matching to the current interval relative-20% response protocol, whose Goryachev peak is 6.4. Profiles show native start, cue withdrawal and final basal state; heatmaps retain full primary fields through basal continuation. Arrays are cyclically rotated only for display to centre the incumbent; positions 0 and 40 are connected, and saved checkpoints remain unrotated. Each model retains its native units and concentration scale. Qualification requires two stable 1000-unit full-field windows and a 10000-unit holdout on frozen conserved-pool scales. No challenger experiment or mesh-convergence claim is included. Generated "+now()+"; producer "+current["producer"]+".")
    path=dev_output_path(figdir,"caption.txt");path.write_text(caption+"\n");files.append(path)
    path=dev_output_path(figdir,"feasibility_summary.json");write(path,dict(generated_at=now(),producer=current["producer"],complete=current["complete"],models=summaries));files.append(path)
    write(dev_output_path(figdir,"manifest.json"),dict(generated_at=now(),producer=current["producer"],artifact_root="../..",
        run=str(run),source_sha256=source_hashes,sha256={dev_artifact_key(figdir,p):sha(p) for p in files},
        renderer=str(Path(__file__).relative_to(ROOT)),renderer_sha256=sha(Path(__file__)),
        canonical_models_changed=False,challenger_trials=0,mesh_adequacy_established=False))
    print(json.dumps(dict(complete=current["complete"], output=str(figdir),
        models=[{key: row[key] for key in ("model","status","preparation_duration","basal_elapsed")} for row in summaries]),indent=2))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("action",choices=["launch","worker","status","audit","collect"])
    parser.add_argument("--run",type=Path,default=DEFAULT)
    parser.add_argument("--slot",type=int)
    parser.add_argument("--cells",type=int,default=400)
    parser.add_argument("--fig",type=Path,default=FIG)
    args=parser.parse_args();run=args.run.resolve()
    if args.action=="launch":launch(run,args.cells)
    elif args.action=="worker":worker(run,args.slot)
    elif args.action=="status":
        result=status(run);print(json.dumps(dict(complete=result["complete"],models=[dict(model=r["job"]["model"],status=r["status"]) for r in result["workers"]]),indent=2))
    elif args.action=="audit":audit(run)
    else:collect(run,args.fig.resolve())


if __name__=="__main__":
    main()
