"""Bounded Holmes mesh and conserved active-Rac imbalance verification.

Integrations run only through SLURM. Native equations and archived experiments
are preserved; every campaign executes frozen inputs and source copies.
"""
from pathlib import Path
import argparse
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
from scipy.signal import find_peaks, peak_widths
import ring_preparation_pilot_20261017 as prep
import chemical_ring_distance_pilot_20261017 as distance
import ring_distance_prototype_20261017 as ring

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "data/runs/holmes_distance_robustness_20261017"
OLD = ROOT / "data/runs/chemical_ring_distance_pilot_20261017"
KEY = "holmes_model3"
SCALES = np.asarray([7.5, 7.5, 3.1, 3.1])
PYTHON = prep.PYTHON
read, write, sha, now = prep.read, prep.write, prep.sha, prep.now


def source_copy(run):
    paths = [Path(__file__)] + [Path(__file__).with_name(name) for name in (
        "ring_preparation_pilot_20261017.py", "chemical_ring_distance_pilot_20261017.py",
        "ring_distance_prototype_20261017.py", "length_dev_output_layout.py")]
    for folder in (ROOT / "code/core/models", ROOT / "code/core/src/polarity1d"):
        paths += sorted(folder.glob("*.py"))
    hashes = {}
    for path in paths:
        relative = path.relative_to(ROOT)
        destination = run / "frozen_sources" / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        hashes[str(relative)] = sha(destination)
    return hashes


def site_cores(grid, state, old_origin, cue_origin):
    """Freeze each terminal site's local half-prominence interval and taper."""
    u = state[0]
    tiled = np.tile(u, 3)
    peaks, properties = find_peaks(tiled, prominence=0.05)
    keep = (peaks >= grid.cells) & (peaks < 2 * grid.cells)
    peaks = peaks[keep]
    assert len(peaks) == 2, "The perturbation requires two resolved terminal sites"
    widths, levels, left, right = peak_widths(tiled, peaks, rel_height=0.5)
    records, smooth, amounts, cores = [], [], [], []
    background = float(u.min())
    for label, origin in (("old", old_origin), ("cue", cue_origin)):
        selected = int(np.argmin(grid.distance(grid.x[peaks % grid.cells], origin)))
        peak = int(peaks[selected] % grid.cells)
        center = float(grid.x[peak])
        a = float((left[selected] - peaks[selected]) * grid.dx)
        b = float((right[selected] - peaks[selected]) * grid.dx)
        core = distance.weights_about(grid, center, a, b)
        arc = grid.arc(grid.x, center)
        side = np.where(arc < 0, -a, b)
        weight = np.where((arc > a) & (arc < b), 0.5 * (1 + np.cos(np.pi * arc / side)), 0)
        weight = weight / (weight.sum() * grid.dx)
        amount = float(np.maximum(u - background, 0) @ core)
        records.append(dict(label=label, peak_cell=peak, center=center, left=a,
            right=b, half_prominence_width=float(widths[selected] * grid.dx),
            half_prominence_level=float(levels[selected]), above_background_amount=amount,
            normalized_weight_integral=float(weight.sum() * grid.dx)))
        smooth.append(weight); amounts.append(amount); cores.append(core)
    assert records[0]["peak_cell"] != records[1]["peak_cell"]
    epsilon = 0.01 * min(amounts)
    delta = epsilon * (smooth[1] - smooth[0])
    assert abs(delta.sum() * grid.dx) < 1e-13
    return dict(background=background, epsilon=epsilon, sites=records), np.asarray(smooth), np.asarray(cores), delta


def submit(run, action, name, array=None):
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = "exec " + shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", PYTHON, str(runner), action, "--run", str(run)])
    if array is not None:
        command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G",
        "--time=01:00:00", "--signal=B:USR1@60", "--job-name=" + name,
        "--output=" + str(run / f"logs/{action}_%A_%a.out"),
        "--error=" + str(run / f"logs/{action}_%A_%a.err")]
    if array is not None:
        args += ["--array=" + array]
    args += ["--wrap", command]
    job = subprocess.check_output(args, text=True).strip()
    write(run / f"launch_{action}.json", dict(generated_at=now(), job_id=job, command=args))
    return job


def launch(run):
    if run.exists():
        raise ValueError("Use a new run directory; frozen experiments are immutable")
    for name in ("inputs", "workers", "logs", "results"):
        (run / name).mkdir(parents=True)
    hashes = source_copy(run)
    # Initialize and freeze precisely one preparation task before any submission.
    preparation = run / "fine_preparation"
    pmanifest = prep.init(preparation, 1600)
    job = pmanifest["jobs"][2]
    job["slot"] = 0
    pmanifest.update(jobs=[job], preparation_durations=[10.0],
        integration_scope="Holmes-only N1600 fresh 10-unit preparation; no preparation-duration scan")
    write(preparation / "manifest.json", pmanifest)
    oldmanifest = read(OLD / "manifest.json")
    oldbase = oldmanifest["models"][KEY]
    for name, source in (("original_basal.npz", OLD / oldbase["checkpoint"]),
                         ("d16_terminal.npz", OLD / "workers/11/recovery.npz")):
        shutil.copy2(source, run / "inputs" / name)
    grid = ring.RingGrid(40, 800)
    with np.load(run / "inputs/d16_terminal.npz") as z:
        state = z["states"][-1].copy()
        assert z["time"][-1] == 13000
    core = oldbase["fixed_original_core"]
    origin = core["origin"]
    perturbation, smooth, cores, delta = site_cores(grid, state, origin, (origin + 16) % 40)
    for slot, sign in enumerate((0, -1, 1)):
        initial = state.copy()
        initial[0] += sign * delta
        assert initial.min() >= 0
        np.testing.assert_array_equal(initial[1:], state[1:])
        np.testing.assert_allclose(ring.pool_amounts(ring.native_model(KEY), grid, initial),
            ring.pool_amounts(ring.native_model(KEY), grid, state), rtol=0, atol=1e-12)
        np.savez_compressed(run / f"inputs/imbalance_{slot}.npz", state=initial,
            reference=state, delta=sign * delta, x=grid.x, smooth_weights=smooth,
            fractional_core_weights=cores, sign=sign)
    manifest = dict(created_at=now(), experiment="Bounded Holmes distance robustness",
        source_sha256=hashes, archived_experiment=str(OLD), old_manifest_sha256=sha(OLD / "manifest.json"),
        input_sha256={p.name: sha(p) for p in (run / "inputs").glob("*.npz")},
        state_names=oldbase["state_names"], fixed_field_scales=SCALES.tolist(),
        solver=oldmanifest["solver"], limits=oldmanifest["limits"],
        fine_protocol=dict(circumference=40, cells=1600, preparation_duration=10,
            distances=[12, 16], cue_width=3, cue_peak=1, cue_duration=1000, basal_duration=13000,
            same_qualified_basal="All fine challenges use one exact qualified N1600 checkpoint"),
        perturbation=dict(**perturbation, convention="positive transfers active Rac from old site core to cue site core",
            definition="1% of smaller fixed half-prominence-core above-minimum active Rac amount; normalized smooth raised-cosine tapers; inactive Rac and both Rho fields unchanged",
            duration=13000, signs=[0, -1, 1]),
        sampling=oldmanifest["sampling"],
        scope="Mesh verification at exactly two distances and finite-window imbalance robustness; no global stability, eigenvalue, or distance-law claim")
    write(run / "manifest.json", manifest)
    jobs = dict(preparation=submit(run, "prepare", "holmes-fine-prep"),
                imbalance=submit(run, "imbalance", "holmes-imbalance", "0-2"))
    write(preparation / "launch.json", dict(generated_at=now(),
        manifest_sha256=sha(preparation / "manifest.json"), job_id=jobs["preparation"],
        wrapper_launch_record=str(run / "launch_prepare.json")))
    write(run / "launch.json", dict(generated_at=now(), jobs=jobs, manifest_sha256=sha(run / "manifest.json")))
    print(json.dumps(jobs, indent=2))


def check_sources(run):
    manifest = read(run / "manifest.json")
    for path, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / path) == digest, path
    for name, digest in manifest["input_sha256"].items():
        assert sha(run / "inputs" / name) == digest
    return manifest


def prepare(run):
    assert os.environ.get("SLURM_JOB_ID")
    check_sources(run)
    prep.worker(run / "fine_preparation", 0)
    record = read(run / "fine_preparation/workers/00/status.json")
    if record["status"] != "qualifying_single_site":
        return
    audit = prep.audit(run / "fine_preparation")
    assert audit["passed"]
    checkpoint = run / "fine_preparation" / record["checkpoint"]
    shutil.copy2(checkpoint, run / "inputs/fine_basal.npz")
    write(run / "fine_basal.json", dict(generated_at=now(), checkpoint_sha256=sha(checkpoint),
        source_checkpoint=str(checkpoint), independent_preparation_audit=audit))
    print(json.dumps(dict(fine_challenge_job=submit(run, "fine", "holmes-fine-distance", "0-1")), indent=2))


def worker(run, action, slot):
    assert os.environ.get("SLURM_JOB_ID"), "Every integration requires SLURM"
    signal.signal(signal.SIGUSR1, lambda s, f: setattr(distance, "STOP", True))
    signal.signal(signal.SIGTERM, lambda s, f: setattr(distance, "STOP", True))
    manifest = check_sources(run)
    model = ring.native_model(KEY)
    directory = run / "workers" / f"{action}_{slot}"
    directory.mkdir()
    record = dict(generated_at=now(), producer="running", complete=False, status="starting",
        action=action, slot=slot, slurm_job_id=os.environ["SLURM_JOB_ID"], phases=[])
    if action == "fine":
        path = run / "inputs/fine_basal.npz"
        assert sha(path) == read(run / "fine_basal.json")["checkpoint_sha256"]
        with np.load(path) as z:
            original = z["state"].copy()
        grid = ring.RingGrid(40, 1600)
        d = (12, 16)[slot]
        core = distance.original_core(grid, original[0])
        center = (core["peak_cell"] + int(round(d / grid.dx))) % grid.cells
        cue = ring.ring_cue(grid, center, 3, 1)
        assert not np.any((core["weights"] > 0) & (cue > 0))
        state = original.copy()
        phases = [("cue", 1000., cue), ("recovery", 13000., np.zeros(grid.cells))]
        record.update(distance=d, input_sha256=sha(path), cue_center_cell=center,
                      original_peak_cell=core["peak_cell"])
    else:
        path = run / f"inputs/imbalance_{slot}.npz"
        with np.load(path) as z:
            state = z["state"].copy()
        with np.load(run / "inputs/original_basal.npz") as z:
            original = z["state"].copy()
        grid = ring.RingGrid(40, 800)
        core = distance.original_core(grid, original[0])
        center = (core["peak_cell"] + 320) % grid.cells
        phases = [("continuation", 13000., np.zeros(grid.cells))]
        record.update(sign=(0, -1, 1)[slot], input_sha256=sha(path), cue_center_cell=center)
    pools = ring.pool_amounts(model, grid, state)
    write(directory / "status.json", record)
    deadline = time.monotonic() + 3300
    try:
        for name, end, cue in phases:
            record.update(status=name, generated_at=now())
            write(directory / "status.json", record)
            result = distance.phase(model, grid, state, cue, end, manifest, deadline, pools)
            metrics = distance.readouts(grid, result["states"], result["time"], original, core, center, model)
            path = directory / f"{name}.npz"
            np.savez_compressed(path, time=result["time"], states=result["states"], x=grid.x,
                stimulus=cue, original_core_weights=core["weights"], **metrics)
            info = {key: value for key, value in result.items() if key not in ("time", "states")}
            info.update(path=str(path.relative_to(run)), sha256=sha(path), phase=name)
            record["phases"].append(info)
            state = result["states"][-1].copy()
            if result["status"] != "complete":
                record["status"] = result["status"]
                break
        else:
            record.update(status="complete", tail=distance.tail_diagnostic(result["time"], result["states"], SCALES))
        record.update(generated_at=now(), producer="finished", complete=True)
    except Exception as error:
        record.update(generated_at=now(), producer="finished", complete=True,
            status="error", error=f"{type(error).__name__}: {error}")
        write(directory / "status.json", record)
        raise
    write(directory / "status.json", record)
    print(json.dumps(record, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["launch", "prepare", "fine", "imbalance"])
    parser.add_argument("--run", type=Path, default=RUN)
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    if args.action == "launch":
        launch(args.run.resolve())
    elif args.action == "prepare":
        prepare(args.run.resolve())
    else:
        worker(args.run.resolve(), args.action, args.slot)


if __name__ == "__main__":
    main()
