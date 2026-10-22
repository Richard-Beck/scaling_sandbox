"""Two-size bridge pilot; joint size/separation change at fixed physical mesh."""
from pathlib import Path
import argparse
import copy
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
import ring_preparation_pilot_20261017 as prep
import chemical_ring_distance_pilot_20261017 as distance
import ring_distance_prototype_20261017 as ring

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / "data/runs/holmes_size_bridge_20261017"
REFERENCE = ROOT / "data/runs/holmes_distance_robustness_20261017"
PYTHON = prep.PYTHON
KEY = "holmes_model3"
SCALES = np.asarray([7.5, 7.5, 3.1, 3.1])
read, write, sha, now = prep.read, prep.write, prep.sha, prep.now


def freeze(run):
    paths = [Path(__file__)] + [Path(__file__).with_name(name) for name in (
        "ring_preparation_pilot_20261017.py", "chemical_ring_distance_pilot_20261017.py",
        "ring_distance_prototype_20261017.py", "length_dev_output_layout.py")]
    for folder in (ROOT / "code/core/models", ROOT / "code/core/src/polarity1d"):
        paths += sorted(folder.glob("*.py"))
    hashes = {}
    for source in paths:
        rel = source.relative_to(ROOT); dest = run / "frozen_sources" / rel
        dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, dest)
        hashes[str(rel)] = sha(dest)
    return hashes


def launch(run):
    if run.exists():
        raise ValueError("Choose a new run; frozen data are retained")
    for name in ("inputs/reference_C40", "workers", "logs", "results"):
        (run / name).mkdir(parents=True)
    hashes = freeze(run)
    template_path = REFERENCE / "fine_preparation/manifest.json"
    template = read(template_path)
    jobs = []
    for slot, circumference, cells in ((0, 30, 1200), (1, 50, 2000)):
        p = run / f"preparation_C{circumference}"
        for name in ("inputs", "workers", "logs", "results"):
            (p / name).mkdir(parents=True)
        (p / "frozen_sources").symlink_to(run / "frozen_sources", target_is_directory=True)
        grid = ring.RingGrid(circumference, cells)
        model = ring.native_model(KEY)
        state = ring.native_rest(KEY, grid)
        assert np.all(np.ptp(state, axis=1) == 0)
        path = p / f"inputs/{KEY}.npz"
        np.savez_compressed(path, state=state, x=grid.x)
        base = copy.deepcopy(template["jobs"][0])
        base.update(slot=0, input=str(path.relative_to(p)), input_sha256=sha(path),
            native_pool_means=(ring.pool_amounts(model, grid, state) / circumference).tolist())
        np.testing.assert_allclose(base["native_pool_means"], [7.5, 3.1], atol=1e-14)
        pmanifest = copy.deepcopy(template)
        pmanifest.update(created_at=now(), circumference=circumference, cells=cells,
            jobs=[base], source_sha256=hashes, preparation_durations=[10.0],
            integration_scope="Exactly one fresh Holmes preparation at this size; no alternative duration attempts")
        write(p / "manifest.json", pmanifest)
        write(p / "status.json", dict(generated_at=now(), complete=False, producer="not started"))
        jobs.append(dict(slot=slot, circumference=circumference, cells=cells, dx=.025,
            distance=.4*circumference, offset_cells=int(.4*cells), preparation=str(p.relative_to(run)),
            preparation_manifest_sha256=sha(p / "manifest.json")))
    copied = {}
    for name, source in (("basal.npz", REFERENCE / "inputs/fine_basal.npz"),
        ("cue.npz", REFERENCE / "workers/fine_1/cue.npz"),
        ("recovery.npz", REFERENCE / "workers/fine_1/recovery.npz"),
        ("status.json", REFERENCE / "workers/fine_1/status.json")):
        dest = run / "inputs/reference_C40" / name
        shutil.copy2(source, dest); copied[str(dest.relative_to(run))] = sha(dest)
    oldmanifest = read(REFERENCE / "manifest.json")
    manifest = dict(created_at=now(), experiment="Bounded joint size/separation bridge", jobs=jobs,
        source_sha256=hashes, reference_input_sha256=copied,
        reference_origin=str(REFERENCE), reference_run_manifest_sha256=sha(REFERENCE / "manifest.json"),
        preparation_protocol_template_sha256=sha(template_path),
        comparison=dict(circumferences=[30, 40, 50], cells=[1200, 1600, 2000], distances=[12, 16, 20],
            challenger_fraction=.4, physical_dx=.025, preparation_peak=1, preparation_duration=10,
            cue_width=3, cue_peak=1, cue_duration=1000, recovery_duration=13000),
        solver=oldmanifest["solver"], limits=oldmanifest["limits"], sampling=oldmanifest["sampling"],
        scope="Exactly two new preparations and two challenges; C40 reused. This joint size/separation change also changes material and incumbent width; no causal separation of these effects.")
    write(run / "manifest.json", manifest)
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = "exec " + shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", PYTHON, str(runner), "worker", "--run", str(run)])
    command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    args = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G",
        "--time=01:00:00", "--signal=B:USR1@60", "--array=0-1", "--job-name=holmes-size-bridge",
        "--output="+str(run/"logs/%A_%a.out"), "--error="+str(run/"logs/%A_%a.err"), "--wrap", command]
    job_id = subprocess.check_output(args, text=True).strip()
    write(run / "launch.json", dict(generated_at=now(), job_id=job_id,
        manifest_sha256=sha(run/"manifest.json"), command=args))
    for job in jobs:
        p = run / job["preparation"]
        write(p / "launch.json", dict(generated_at=now(), job_id=job_id, slot=job["slot"],
            manifest_sha256=sha(p/"manifest.json"), wrapper_record=str(run/"launch.json")))
    print(json.dumps(dict(job_id=job_id, workers=len(jobs), run=str(run)), indent=2))


def worker(run, slot):
    assert os.environ.get("SLURM_JOB_ID"), "All integration requires SLURM"
    manifest = read(run / "manifest.json")
    assert sha(run/"manifest.json") == read(run/"launch.json")["manifest_sha256"]
    for path, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / path) == digest
    job = manifest["jobs"][slot]
    directory = run / "workers" / f"{slot:02d}"
    directory.mkdir()
    status = dict(generated_at=now(), producer="running", complete=False, status="preparation", job=job, phases=[])
    write(directory/"status.json", status)
    p = run / job["preparation"]
    prep.worker(p, 0)
    pstatus = read(p/"workers/00/status.json")
    status["preparation_status"] = pstatus["status"]
    if pstatus["status"] != "qualifying_single_site":
        status.update(generated_at=now(), producer="finished", complete=True, status="preparation_failed")
        write(directory/"status.json", status)
        return
    audit = prep.audit(p)
    assert audit["passed"]
    source = p / pstatus["checkpoint"]
    checkpoint = directory / "qualified_basal.npz"
    shutil.copy2(source, checkpoint)
    assert sha(checkpoint) == pstatus["checkpoint_sha256"]
    with np.load(checkpoint) as z:
        original = z["state"].copy()
    grid = ring.RingGrid(job["circumference"], job["cells"])
    model = ring.native_model(KEY)
    core = distance.original_core(grid, original[0])
    center = (core["peak_cell"] + job["offset_cells"]) % grid.cells
    cue = ring.ring_cue(grid, center, 3, 1)
    status.update(checkpoint_sha256=sha(checkpoint), original_peak_cell=core["peak_cell"],
        cue_center_cell=center, original_core_width=core["width"])
    if np.any((core["weights"] > 0) & (cue > 0)):
        status.update(generated_at=now(), producer="finished", complete=True, status="geometric_overlap_exclusion")
        write(directory/"status.json", status)
        return
    signal.signal(signal.SIGUSR1, lambda s,f: setattr(distance,"STOP",True))
    signal.signal(signal.SIGTERM, lambda s,f: setattr(distance,"STOP",True))
    state = original.copy(); pools = ring.pool_amounts(model, grid, original)
    deadline = time.monotonic() + 3300
    try:
        for name,end,stimulus in (("cue",1000.,cue),("recovery",13000.,np.zeros(grid.cells))):
            status.update(generated_at=now(),status=name); write(directory/"status.json",status)
            result = distance.phase(model, grid, state, stimulus, end, manifest, deadline, pools)
            metrics = distance.readouts(grid, result["states"], result["time"], original, core, center, model)
            path = directory / f"{name}.npz"
            np.savez_compressed(path,time=result["time"],states=result["states"],x=grid.x,
                stimulus=stimulus,original_core_weights=core["weights"],**metrics)
            info = {key:value for key,value in result.items() if key not in ("time","states")}
            info.update(phase=name,path=str(path.relative_to(run)),sha256=sha(path))
            status["phases"].append(info); state=result["states"][-1].copy()
            if result["status"] != "complete":
                status["status"] = result["status"]
                break
        else:
            status.update(status="complete",tail=distance.tail_diagnostic(result["time"],result["states"],SCALES))
        status.update(generated_at=now(),producer="finished",complete=True)
    except Exception as error:
        status.update(generated_at=now(),producer="finished",complete=True,status="error",error=f"{type(error).__name__}: {error}")
        write(directory/"status.json",status)
        raise
    write(directory/"status.json",status)
    print(json.dumps(status,indent=2))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("action",choices=["launch","worker"])
    parser.add_argument("--run",type=Path,default=RUN); parser.add_argument("--slot",type=int)
    args=parser.parse_args()
    if args.action=="launch": launch(args.run.resolve())
    else: worker(args.run.resolve(),args.slot)


if __name__=="__main__": main()
