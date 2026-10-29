"""Minimal Spring per-cell profile replay for the R2 display rebuild.

The source arm retained only modes 0/1. This replays its fixed first replicate
at each outer size, with unchanged dynamics and random streams, and adds native
active profiles to the first 200 post-burn time units. No estimator is changed.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "polarity1d/manuscript").is_dir())
MS = ROOT / "polarity1d/manuscript"
SOURCE = MS / "data/runs/ring_r2_per_cell_20261027"
RUN = MS / "data/runs/ring_r2_display_profiles_20261029"
PY = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def now():
    return datetime.now(timezone.utc).isoformat()


def prepare():
    if (RUN / "manifest.json").exists():
        raise RuntimeError("Replay is already prepared; frozen inputs must remain unchanged")
    original = SOURCE / "frozen_sources/engine.py"
    code = original.read_text()
    # Add output writes only. All solver, hazard, event and native dynamics
    # expressions remain byte-for-byte identical to the source engine.
    changes = [
        ("    nm = 2\n", "    nm = 2\n    profiles = np.lib.format.open_memmap(folder / 'u_profiles.npy', mode='w+', dtype='float32', shape=(len(times), n))\n"),
        ('            saved["u"][next_sample] = coefficients[0]\n', '            saved["u"][next_sample] = coefficients[0]\n            profiles[next_sample] = fields[0]\n'),
        ('    np.save(folder / "time.npy", times)\n', '    profiles.flush()\n    np.save(folder / "time.npy", times)\n'),
    ]
    for old, new in changes:
        if code.count(old) != 1:
            raise RuntimeError("Output transformation no longer matches the frozen source")
        code = code.replace(old, new)
    for part in ("frozen_sources", "workers", "logs"):
        (RUN / part).mkdir(parents=True, exist_ok=True)
    for part in ("native", "inputs"):
        (RUN / part).symlink_to(SOURCE / part, target_is_directory=True)
    (RUN / "frozen_sources/engine.py").write_text(code)
    shutil.copy2(__file__, RUN / "frozen_sources/replay.py")
    source = read(SOURCE / "manifest.json")
    jobs = []
    for slot in (80, 116):
        original_job = source["jobs"][slot]
        assert original_job["model"] == "spring_ring"
        assert original_job["replicate"] == 0
        assert original_job["rate_level"] == "low" and original_job["amplitude_level"] == "high"
        job = dict(original_job)
        job.update(slot=len(jobs), replay_source_slot=slot, sample_duration=200., wall_seconds=1200.)
        jobs.append(job)
    manifest = dict(
        created_utc=now(), jobs=jobs, maximum_extension_level=0,
        basals=[b for b in source["basals"] if b["model"] == "spring_ring" and b["size"] in (1., 3.)],
        source_campaign=str(SOURCE.relative_to(ROOT)),
        source_manifest_sha256=sha(SOURCE / "manifest.json"),
        source_engine_sha256=sha(original),
        source_sha256={str(p.relative_to(RUN)): sha(p) for p in sorted((RUN / "frozen_sources").glob("*.py"))},
        output_only_changes=changes,
        selection="Fixed source replicate 0 at both outer sizes; mechanically first 200 post-burn native time units; no seed/window search",
        changes_to_jobs="Only local slot, output duration (12800 to 200), wall-clock limit and replay source reference; original .05 cadence retained",
        outputs="u_profiles.npy: frame by native cell-centre, float32 active concentration; original modes0/1, time, circumference and events also retained",
        scientific_status="Display-only replay; original arm remains authoritative for all persistence estimates",
    )
    write(RUN / "manifest.json", manifest)
    write(RUN / "freeze.json", dict(manifest_sha256=sha(RUN / "manifest.json")))
    (RUN / "worker.sh").write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec " +
        shlex.join([PY, "-B", str(RUN / "frozen_sources/replay.py"), "worker", "--slot"]) + ' "$SLURM_ARRAY_TASK_ID"\n')
    (RUN / "worker.sh").chmod(0o750)
    (RUN / "README.md").write_text(
        "# Minimal Spring per-cell display replay\n\n"
        "Fixed original replicate 0, native resting sizes 1 and 3, low event density and high amplitude. "
        "Native basal states, seeds, stochastic streams, BDF dynamics, hazard and tolerances are unchanged. "
        "Only direct native active-profile writes are added. The first 200 post-burn time units are retained at the original 0.05 cadence. "
        "Original campaign records are read only and remain authoritative for persistence estimates. "
        "Selection is fixed mechanically without consulting heading changes.\n\n"
        "Run `ring_r2_display_profiles_20261029.py prepare`, `submit`, then `validate`. "
        "The validation compares all shared events and native-cadence P/C with the original arm.\n")


def verify():
    manifest = read(RUN / "manifest.json")
    assert sha(RUN / "manifest.json") == read(RUN / "freeze.json")["manifest_sha256"]
    assert sha(SOURCE / "manifest.json") == manifest["source_manifest_sha256"]
    for name, digest in manifest["source_sha256"].items():
        assert sha(RUN / name) == digest
    return manifest


def worker(slot):
    verify()
    spec = importlib.util.spec_from_file_location("r2_display_engine", RUN / "frozen_sources/engine.py")
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    engine.RUN = RUN
    engine.worker(slot, 0)


def submit():
    verify()
    cmd = ["sbatch", "--parsable", "--partition=red", "--qos=small", "--cpus-per-task=1", "--mem=2G", "--time=00:30:00",
           "--job-name=R2_profiles", "--array=0-1", "--output=" + str(RUN / "logs/%A_%a.out"),
           "--error=" + str(RUN / "logs/%A_%a.err"), str(RUN / "worker.sh")]
    result = subprocess.run(cmd, text=True, capture_output=True)
    write(RUN / "submission.json", dict(generated_utc=now(), command=cmd, returncode=result.returncode,
                                          stdout=result.stdout, stderr=result.stderr, array_throttle=None))
    if result.returncode:
        raise RuntimeError(result.stderr)
    print(result.stdout.strip())


def validate():
    manifest = verify()
    records = []
    for job in manifest["jobs"]:
        folder = RUN / "workers" / f'{job["slot"]:03d}' / "segments/00"
        source = SOURCE / "workers" / f'{job["replay_source_slot"]:03d}' / "segments/00"
        status = read(folder / "status.json")
        assert status["complete"]
        times = np.load(folder / "time.npy")
        source_times = np.load(source / "heading_time.npy", mmap_mode="r")
        index = np.searchsorted(source_times, times)
        assert np.array_equal(source_times[index], times)
        u = np.load(folder / "u_modes.npy")
        P = u[:, 1].conj() / u[:, 0]
        source_P = np.load(source / "P.npy", mmap_mode="r")[index]
        C = np.load(folder / "circumference.npy")
        source_C = np.load(source / "circumference.npy", mmap_mode="r")[index]
        profiles = np.load(folder / "u_profiles.npy", mmap_mode="r")
        events = np.loadtxt(folder / "events.csv", delimiter=",", ndmin=2)
        source_events = np.loadtxt(source / "events.csv", delimiter=",", ndmin=2)
        source_events = source_events[source_events[:, 0] <= times[-1] + 1e-8]
        assert events.shape == source_events.shape
        assert np.allclose(events, source_events, rtol=0, atol=1e-7)
        assert np.allclose(P, source_P, rtol=2e-5, atol=2e-6)
        assert np.allclose(C, source_C, rtol=2e-7, atol=2e-8)
        n = profiles.shape[1]
        phase = np.exp(-1j * np.pi * np.arange(2) / n)
        reconstructed = np.fft.rfft(profiles, axis=1)[:, :2] * phase / n
        assert np.allclose(reconstructed, u, rtol=2e-5, atol=2e-7)
        records.append(dict(size=job["size"], seed=job["seed"], source_slot=job["replay_source_slot"],
            replay_slot=job["slot"], frames=len(times), native_cells=n,
            native_interval=[float(times[0]), float(times[-1])], events=len(events),
            max_event_abs_difference=float(np.max(abs(events-source_events))),
            max_P_abs_difference=float(np.max(abs(P-source_P))),
            max_C_abs_difference=float(np.max(abs(C-source_C))),
            max_profile_mode_abs_difference=float(np.max(abs(reconstructed-u))),
            numerical_checks=status["checks"],
            allocated_bytes=sum(p.stat().st_blocks*512 for p in folder.iterdir() if p.is_file()),
            output_sha256={p.name:sha(p) for p in folder.glob("*.npy")}))
    write(RUN / "validation.json", dict(generated_utc=now(), passed=True, records=records,
        criteria="All shared events match at absolute1e-7; P rtol2e-5/atol2e-6; C rtol2e-7/atol2e-8; full-profile modes agree with saved modes",
        original_records_read_only=True))
    print(json.dumps(read(RUN / "validation.json"), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "submit", "worker", "validate"))
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    if args.command == "worker":
        worker(args.slot)
    else:
        globals()[args.command]()
