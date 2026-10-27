"""Compact Poisson ring timing assay and reusable production evolution engine.

Native BDF evolution, exact event restarts, physical Gaussian pulses and
state-dependent perimeter intensity for the moving Spring ring. No profile
histories are written. The default CLI submits timing pilots only; the separate
production orchestrator freezes its own manifest and imports this engine.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy import sparse
from scipy.integrate import BDF
from scipy.optimize import brentq
from scipy.signal import fftconvolve

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "polarity1d/manuscript").is_dir())
MS = ROOT / "polarity1d/manuscript"
BASE = MS / "data/runs/signature_figure4_round3_compact_20261022"
RUN = MS / "data/runs/poisson_ring_modes_timing_20261027"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring_ring")
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes", spring_ring="Spring")
REFERENCE = dict(wave_pinning=.1, goryachev=6.4, holmes_model3=1., spring_ring=4.)
LEVELS = ("low", "medium", "high")
RATES = (.0001, .0003, .001)
AMPLITUDES = (.3, 1., 3.)


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def table(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def configuration(model, size, rate_index, amplitude_index, replicate=0):
    spring = model == "spring_ring"
    # Spring retains its uncalibrated native length/time scales.
    ell, tau = (.15, 5.) if spring else (1.5, 50.)
    rho = RATES[rate_index] * (100 if spring else 1)
    return dict(model=model, size=float(size), rate_level=LEVELS[rate_index],
                amplitude_level=LEVELS[amplitude_index], rho=rho,
                amplitude_mean=REFERENCE[model] * AMPLITUDES[amplitude_index],
                amplitude_reference_fraction=AMPLITUDES[amplitude_index],
                reference=REFERENCE[model], ell=ell, sigma=ell / math.sqrt(2), tau=tau,
                replicate=replicate, frame_interval=.05 if spring else .5,
                burn=100. if spring else 1000., sample_duration=200. if spring else 2000.,
                rtol=1e-8, atol=1e-10, max_step=.1 if spring else 1.,
                # Includes the mode at exactly the retained minimum wavelength.
                maximum_mode=int(math.floor(size / ell + 1e-9)))


def prepare():
    if RUN.exists():
        raise ValueError("Choose a fresh immutable campaign directory")
    for part in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (RUN / part).mkdir(parents=True)
    shutil.copytree(BASE / "frozen_sources", RUN / "native/frozen_sources")
    worker = RUN / "frozen_sources" / Path(__file__).name
    shutil.copy2(__file__, worker)
    basals = []
    for model in MODELS:
        for size in ((1., 2., 3.) if model == "spring_ring" else (20., 40., 60.)):
            p = BASE / f"inputs/{model}_C{size:g}.npz"
            dest = RUN / "inputs" / p.name
            shutil.copy2(p, dest)
            basals.append(dict(model=model, size=size, input=str(dest.relative_to(RUN)), sha256=sha(dest)))
    jobs = []
    for model in MODELS:
        sizes = (1., 2., 3.) if model == "spring_ring" else (20., 40., 60.)
        for size in sizes:
            job = configuration(model, size, 1, 1)
            jobs.append(job)
        # Largest-size high-rate/high-amplitude tail stresses event overhead.
        jobs.append(configuration(model, sizes[-1], 2, 2))
    for slot, job in enumerate(jobs):
        job.update(slot=slot, seed=202610270000 + slot)
    production = [configuration(model, size, ir, ia, rep)
                  for model in MODELS
                  for size in ((1., 2., 3.) if model == "spring_ring" else (20., 40., 60.))
                  for ir in range(3) for ia in range(3) for rep in range(4)]
    assert len(production) == 432
    sources = {str(p.relative_to(RUN)): sha(p) for p in sorted((RUN / "native").rglob("*.py"))}
    sources[str(worker.relative_to(RUN))] = sha(worker)
    manifest = dict(created_utc=now(), state="TIMING_PILOTS_ONLY", production_submitted=False,
                    timing_jobs=jobs, basals=basals, source_sha256=sources,
                    proposed_production_count=432, proposed_production_conditions=108,
                    production_repeats=4, production_horizon="Not selected or authorized for launch yet",
                    rate_levels_chemical=RATES, rate_levels_spring=[100 * r for r in RATES],
                    amplitude_reference_fractions=AMPLITUDES,
                    forcing="Homogeneous Poisson events per physical perimeter; exponential independent positive peak marks; uniform material centres; fixed physical Gaussian sigma=ell/sqrt(2) and exponential temporal decay tau",
                    spring_intensity="Exact time change: cumulative hazard derivative rho*C(t); unit-exponential hazard thresholds; no fixed basal-perimeter approximation",
                    spring_pulse_width="Recomputed at instantaneous physical circumference, including analytic cue derivative in the Jacobian",
                    tail_cutoff_tau=32, solver="Native BDF, exact burst restarts, no state clipping or mass adjustment",
                    retained="Per-frame complex64 u and v modes, float64 time and Spring circumference; no profiles or separate heading/coherence",
                    resources=dict(qos="small", cpus=1, memory="2G", walltime="00:30:00", array_throttle=None))
    write(RUN / "manifest.json", manifest)
    write(RUN / "freeze.json", dict(manifest_sha256=sha(RUN / "manifest.json")))
    script = RUN / "worker.sh"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nexec "
                      + shlex.join([PYTHON, "-B", str(worker), "worker", "--slot"]) + ' "$SLURM_ARRAY_TASK_ID"\n')
    script.chmod(0o750)
    (RUN / "README.md").write_text(
        "# Poisson ring Fourier timing pilot\n\n"
        "Only the 16 bounded timing jobs are submitted. The proposed 432 production histories are held for user timing review.\n\n"
        "Native BDF evolves the exact retained polarized ring preparations. Samples retain u_modes.npy and v_modes.npy, ordered frame then nonnegative mode including zero, plus time.npy and Spring circumference.npy. Fourier phases include the half-cell offset, so the heading vector is conj(u_modes[:,1])/u_modes[:,0]. No profile histories are saved.\n\n"
        "Chemical rings: C=20/40/60, correlation length 1.5, pulse SD=1.5/sqrt(2), temporal decay 50. Spring: native resting sizes 1/2/3, correlation length .15, pulse SD=.15/sqrt(2), decay 5; these scales are not biologically calibrated between models.\n\n"
        "Arrival densities are .0001/.0003/.001 in chemical native units and .01/.03/.1 in Spring native units. Mean exponential peak amplitudes are .3/1/3 times each model's established reference (.1, 6.4, 1, 4 respectively). Neither total input nor realized noise is normalized by circumference.\n\n"
        "Burn-in is deliberately short for resource timing (1000 chemical, 100 Spring), followed by 2000/200 sampled time units. These are timing benchmarks, not qualified stationary persistence estimates. Production burn-in, horizon and frame cadence remain to be selected after reviewing this pilot.\n\n"
        "Forcing pulses retain fixed physical width as Spring deforms. Its instantaneous arrival rate is rho*C(t). Pulses older than 32 decay times are omitted (relative amplitude below 1.27e-14). Native basal chemistry and all conserved pools remain active.\n")
    print(json.dumps(dict(run=str(RUN), timing_jobs=len(jobs), production_submitted=False)))


class Bursts:
    def __init__(self, job, cells):
        streams = np.random.SeedSequence(job["seed"]).spawn(3)
        self.wait = np.random.default_rng(streams[0])
        self.marks = np.random.default_rng(streams[1])
        self.positions = np.random.default_rng(streams[2])
        self.job = job
        self.q = (np.arange(cells) + .5) / cells
        self.events = []
        self.count = 0
        self.peak_sum = 0.
        self.peak_max = 0.

    def threshold(self):
        return float(self.wait.exponential())

    def append(self, t):
        peak = float(self.marks.exponential(self.job["amplitude_mean"]))
        center = float(self.positions.uniform())
        self.events.append((t, peak, center))
        if hasattr(self, "log_path"):
            with self.log_path.open("a") as event_file:
                event_file.write(f"{t:.17g},{peak:.17g},{center:.17g}\n")
        self.count += 1
        self.peak_sum += peak
        self.peak_max = max(peak, self.peak_max)

    def value(self, t, circumference):
        signal = np.zeros(len(self.q))
        derivative = np.zeros(len(self.q))
        tau, sigma = self.job["tau"], self.job["sigma"]
        # Mutation only removes terms already below the documented cutoff.
        for birth, peak, center in self.events:
            age = t - birth
            if age > 32 * tau:
                continue
            distance_fraction = (self.q - center + .5) % 1 - .5
            # Periodization remains smooth at the seam and opposite point.
            for image in (-1, 0, 1):
                d = distance_fraction + image
                pulse = peak * np.exp(-age / tau - .5 * (circumference * d / sigma)**2)
                signal += pulse
                derivative -= circumference * d*d / sigma**2 * pulse
        return signal, derivative


def native(platform):
    p = RUN / "native/frozen_sources/campaign.py"
    spec = importlib.util.spec_from_file_location("poisson_native_loader", p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.native_module(RUN / "native", platform)


def worker(slot, level=0):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Evolution requires SLURM")
    manifest = read(RUN / "manifest.json")
    if sha(RUN / "manifest.json") != read(RUN / "freeze.json")["manifest_sha256"]:
        raise ValueError("Manifest changed")
    for name, digest in manifest["source_sha256"].items():
        if sha(RUN / name) != digest:
            raise ValueError("Frozen source changed: " + name)
    job = dict(manifest.get("jobs", manifest.get("timing_jobs"))[slot])
    production = "jobs" in manifest
    if production:
        if not 0 <= level <= manifest["maximum_extension_level"]:
            raise ValueError("Undeclared extension level")
        job.update(extension_level=level, initial_sample_duration=job["sample_duration"],
                   sample_duration=job["sample_duration"] * 2**level)
    folder = RUN / "workers" / f"{slot:03d}"
    if production:
        folder = folder / "segments" / f"{level:02d}"
    folder.mkdir(parents=True, exist_ok=True)
    if (folder / "status.json").exists() and read(folder / "status.json").get("complete"):
        return
    started = time.monotonic()
    write(folder / "status.json", dict(complete=False, state="RUNNING", job=job, generated_utc=now()))
    try:
        evolve(job, folder, started)
    except BaseException as exc:
        write(folder / "status.json", dict(complete=False, state="FAILED", job=job,
                                           generated_utc=now(), error=repr(exc)))
        raise


def evolve(job, folder, started):
    basal = next(b for b in read(RUN / "manifest.json")["basals"]
                 if b["model"] == job["model"] and b["size"] == job["size"])
    p = RUN / basal["input"]
    if sha(p) != basal["sha256"]:
        raise ValueError("Basal changed")
    with np.load(p) as z:
        state = z["state"].copy()
    spring = job["model"] == "spring_ring"
    N = native("spring" if spring else "chemical")
    if spring:
        n = (len(state) - 1) // 2
        closure = N.Closure(job["size"], n)
        initial_pools = np.array([closure.amount(state)])
    else:
        nf, n = state.shape
        grid = N.ring.RingGrid(job["size"], n)
        model = N.ring.native_model(job["model"])
        initial_pools = N.ring.pool_amounts(model, grid, state)
        state = state.ravel()
    physical_dimension = len(state)
    forcing = Bursts(job, n)
    forcing.log_path = folder / "events.csv"
    y = np.r_[state, 0.]
    t0 = 0.
    level = job.get("extension_level", 0)
    if level:
        previous = folder.parent / f"{level-1:02d}"
        if not read(previous / "status.json").get("complete"):
            raise ValueError("Cannot extend an incomplete prefix")
        checkpoint = read(previous / "checkpoint.json")
        with np.load(previous / "checkpoint.npz") as z:
            y = z["state"].copy()
        if sha(previous / "checkpoint.npz") != checkpoint["state_sha256"]:
            raise ValueError("Checkpoint changed")
        t0 = checkpoint["native_time"]
        forcing.events = [tuple(e) for e in checkpoint["events"]]
        forcing.count = checkpoint["event_count"]
        forcing.peak_sum = checkpoint["peak_sum"]
        forcing.peak_max = checkpoint["peak_max"]
        for name in ("wait", "marks", "positions"):
            getattr(forcing, name).bit_generator.state = checkpoint["rng"][name]

    def rhs(t, y):
        physical = y[:-1]
        C = physical[-1] if spring else job["size"]
        cue, _ = forcing.value(t, C)
        if spring:
            delta = closure.rhs(physical, cue)
        else:
            delta = N.ring.native_rhs(model, grid, physical.reshape(nf, n), cue).ravel()
        return np.r_[delta, job["cell_rate"]]

    kwargs = dict(rtol=job["rtol"], atol=job["atol"], max_step=job["max_step"])
    if spring:
        def jac(t, y):
            cue, dC = forcing.value(t, y[-2])
            J = sparse.csc_matrix(closure.jacobian(y[:-1], cue, dC))
            bottom = sparse.csc_matrix((1, physical_dimension))
            return sparse.bmat([[J, None], [bottom, sparse.csc_matrix((1, 1))]], format="csc")
        kwargs["jac"] = jac
    else:
        kwargs["jac_sparsity"] = sparse.block_diag(
            (N.ring.periodic_sparsity(nf, n), sparse.csr_matrix((1, 1))), format="csc")

    end = job["burn"] + job["sample_duration"]
    first = t0 + job["frame_interval"] if level else job["burn"]
    times = first + np.arange(int(round((end-first) / job["frame_interval"])) + 1) * job["frame_interval"]
    nm = 2
    saved = {species: np.lib.format.open_memmap(folder / f"{species}_modes.npy", mode="w+", dtype="complex64", shape=(len(times), nm))
             for species in ("u",)}
    circles = np.lib.format.open_memmap(folder / "circumference.npy", mode="w+", dtype="float64", shape=(len(times),)) if spring else None
    phase = np.exp(-1j * np.pi * np.arange(nm) / n)
    next_sample = 0
    t = t0
    threshold = checkpoint["next_hazard_threshold"] if level else forcing.threshold()
    minimum = float(state[:-1].min() if spring else state.min())
    mass_error = 0.
    steps = nfev = njev = nlu = 0
    last_report = time.monotonic()
    sampling_start = None
    sample_events_start = None
    sample_hazard_start = None
    solver_start = time.monotonic()
    solver = BDF(rhs, t, y, end, **kwargs)
    while t < end - 1e-10:
        message = solver.step()
        if solver.status == "failed":
            raise RuntimeError(message)
        steps += 1
        dense = solver.dense_output()
        crossed = solver.y[-1] >= threshold
        stop = brentq(lambda s: float(dense(s)[-1] - threshold), solver.t_old, solver.t, xtol=1e-11) if crossed else solver.t
        while next_sample < len(times) and times[next_sample] <= stop + 1e-10:
            if sampling_start is None:
                sampling_start = time.monotonic()
                sample_events_start = forcing.count
                sample_hazard_start = float(dense(first)[-1])
            sample = dense(float(times[next_sample]))[:-1]
            fields = sample[:-1].reshape(2, n) if spring else sample.reshape(nf, n)
            coefficients = np.fft.rfft(fields[:2], axis=1)[:, :nm] * phase / n
            saved["u"][next_sample] = coefficients[0]
            
            if spring:
                circles[next_sample] = sample[-1]
                pools = np.array([closure.amount(sample)])
            else:
                pools = N.ring.pool_amounts(model, grid, fields)
            minimum = min(minimum, float(fields.min()))
            mass_error = max(mass_error, float(np.max(abs(pools / initial_pools - 1))))
            next_sample += 1
        t = float(stop)
        y = np.asarray(dense(t), float)
        minimum = min(minimum, float(y[:-2].min() if spring else y[:-1].min()))
        if not np.isfinite(y).all() or minimum < -5e-9 or mass_error > 1e-7:
            raise ValueError(f"Numerical guard failed: minimum={minimum}, pool={mass_error}")
        if crossed:
            nfev += solver.nfev
            njev += solver.njev
            nlu += solver.nlu
            forcing.append(t)
            # Prune only after accepted event times; BDF may evaluate times out of order.
            forcing.events = [e for e in forcing.events if t - e[0] < 32 * job["tau"]]
            threshold += forcing.threshold()
            if t < end - 1e-10:
                solver = BDF(rhs, t, y, end, **kwargs)
        if time.monotonic() - last_report > 30:
            for array in saved.values():
                array.flush()
            write(folder / "progress.json", dict(generated_utc=now(), time=t, target=end,
                  frames=next_sample, elapsed_seconds=time.monotonic()-solver_start, events=forcing.count))
            last_report = time.monotonic()
        if time.monotonic() - solver_start > job.get("wall_seconds", 1650):
            raise TimeoutError("Worker exceeded its declared runtime limit")
    completed = time.monotonic()
    if next_sample != len(times):
        raise ValueError("Missing recorded Fourier frames")
    for array in saved.values():
        array.flush()
    if circles is not None:
        circles.flush()
    np.save(folder / "time.npy", times)
    if "extension_level" in job:
        tmp = folder / "checkpoint.npz.tmp"
        with tmp.open("wb") as handle:
            np.savez_compressed(handle, state=y)
        tmp.replace(folder / "checkpoint.npz")
        write(folder / "checkpoint.json", dict(native_time=t, next_hazard_threshold=threshold,
              state_sha256=sha(folder / "checkpoint.npz"), events=forcing.events,
              event_count=forcing.count, peak_sum=forcing.peak_sum, peak_max=forcing.peak_max,
              rng={name: getattr(forcing, name).bit_generator.state for name in ("wait", "marks", "positions")}))
    io_completed = time.monotonic()
    sampling_seconds = io_completed - sampling_start
    burn_seconds = sampling_start - solver_start
    data_bytes = sum(p.stat().st_size for p in folder.glob("*.npy"))
    checks = dict(minimum_concentration=minimum, max_relative_pool_error=mass_error,
                  events_total=forcing.count, events_postburn=forcing.count-sample_events_start,
                  cumulative_hazard_postburn=float(y[-1]-sample_hazard_start),
                  empirical_mean_peak=forcing.peak_sum/max(1, forcing.count), maximum_peak=forcing.peak_max,
                  frames=next_sample, modes=nm, data_bytes=data_bytes,
                  integration_steps=steps, nfev=nfev+solver.nfev, njev=njev+solver.njev, nlu=nlu+solver.nlu)
    write(folder / "status.json", dict(complete=True, state="COMPLETE", generated_utc=now(), job=job,
          slurm_job_id=os.environ["SLURM_JOB_ID"], node=os.uname().nodename,
          setup_seconds=solver_start-started, burn_wall_seconds=burn_seconds,
          postburn_wall_seconds=sampling_seconds, postburn_duration=end-max(t0, job["burn"]),
          seconds_per_10000_postburn=sampling_seconds*10000/(end-max(t0, job["burn"])),
          elapsed_seconds=io_completed-started, checks=checks,
          caveat="Short benchmark; startup/short-transient cost and machine load may differ in long runs"))
    print(json.dumps(read(folder / "status.json")), flush=True)


def acf_raw(coefficients):
    """Uncentered complex-mode numerator; each lag uses its actual pair count."""
    x = np.asarray(coefficients, dtype=np.complex128)
    count = len(x)
    raw = fftconvolve(x, np.conj(x[::-1]), mode="full")[count-1:count+count//2]
    return raw.real / np.arange(count, count-len(raw), -1)


def half_rate(curve, cadence):
    if curve[1] <= .5:
        return None, "below_frame_resolution", math.log(2) / cadence
    hits = np.flatnonzero(curve[1:] <= .5) + 1
    if not len(hits):
        return None, "beyond_analysis_lag", math.log(2) / ((len(curve)-1)*cadence)
    j = int(hits[0])
    lag = cadence * ((j-1) + (curve[j-1]-.5)/(curve[j-1]-curve[j]))
    return math.log(2) / lag, "resolved", None


def collect():
    m = read(RUN / "manifest.json")
    rows = []
    incomplete = []
    for job in m["timing_jobs"]:
        folder = RUN / "workers" / f"{job['slot']:03d}"
        p = folder / "status.json"
        if not p.exists() or not read(p).get("complete"):
            incomplete.append(job["slot"])
            continue
        s = read(p)
        z = np.load(folder / "u_modes.npy", mmap_mode="r")
        crossing = []
        for n in (1, job["maximum_mode"]):
            raw = acf_raw(z[:, n])
            curve = raw/raw[0]
            gamma, status, bound = half_rate(curve, job["frame_interval"])
            crossing.append((gamma, status, bound))
        rows.append(dict(model=job["model"], size=job["size"], rate_level=job["rate_level"],
                         amplitude_level=job["amplitude_level"], rho=job["rho"],
                         amplitude_mean=job["amplitude_mean"], sample_duration=job["sample_duration"],
                         burn_wall_seconds=s["burn_wall_seconds"],
                         sample_wall_seconds=s["postburn_wall_seconds"],
                         minutes_per_10000=s["seconds_per_10000_postburn"]/60,
                         events_postburn=s["checks"]["events_postburn"],
                         data_mib=s["checks"]["data_bytes"]/2**20,
                         gamma1=crossing[0][0], mode1_status=crossing[0][1], mode1_bound=crossing[0][2],
                         gamma_shortest=crossing[1][0], shortest_status=crossing[1][1],
                         pool_error=s["checks"]["max_relative_pool_error"],
                         minimum=s["checks"]["minimum_concentration"], node=s["node"]))
    if rows:
        table(RUN / "results/timing.csv", rows)
    write(RUN / "results/status.json", dict(generated_utc=now(), complete=len(rows), expected=16,
          missing_or_failed=incomplete, production_submitted=False))
    print(json.dumps(dict(rows=rows, incomplete=incomplete, production_submitted=False), indent=2))


def submit():
    m = read(RUN / "manifest.json")
    if (RUN / "submission.json").exists():
        raise ValueError("Timing pilots already submitted")
    cmd = ["sbatch", "--parsable", "--job-name=poisson_timing", "--partition=red", "--qos=small",
           "--cpus-per-task=1", "--mem=2G", "--time=00:30:00", "--array=0-15",
           "--output="+str(RUN / "logs/%A_%a.out"), "--error="+str(RUN / "logs/%A_%a.err"), str(RUN / "worker.sh")]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    write(RUN / "submission.json", dict(created_utc=now(), job_id=result.stdout.strip(), command=cmd,
          scope="16 bounded timing pilots only", production_submitted=False))
    print(result.stdout.strip())


def checks():
    # Analytic autocorrelation check including a non-real complex mean.
    x = np.array([1+2j, 3-1j, -2+.5j, 4+3j, -.5-2j])
    got = acf_raw(x)
    wanted = np.array([np.real(np.mean(x[k:]*x[:len(x)-k].conj())) for k in range(len(got))])
    assert np.allclose(got, wanted, atol=1e-12)
    rate, status, _ = half_rate(np.array([1., .75, .5, .25]), 2.)
    assert status == "resolved" and np.isclose(rate, math.log(2)/4)
    job = configuration("wave_pinning", 60., 1, 1)
    job["seed"] = 20261027
    b = Bursts(job, 1200)
    b.events = [(0., 2., .3)]
    a, dC = b.value(3., 60.)
    finite = (b.value(3., 60.00001)[0]-b.value(3., 59.99999)[0])/.00002
    assert np.max(abs(finite-dC)) < 1e-9
    assert np.allclose(b.value(3.+job["tau"], 60.)[0], a/math.e)
    rng = np.random.default_rng(17)
    waits = rng.exponential(size=200000)
    marks = rng.exponential(job["amplitude_mean"], size=200000)
    assert abs(waits.mean()-1) < .01
    assert abs(marks.mean()/job["amplitude_mean"]-1) < .01
    result = dict(generated_utc=now(), autocorrelation_pair_count="pass", half_crossing="pass",
                  physical_pulse_derivative="pass", exponential_waits_and_marks="pass")
    print(json.dumps(result))
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=("prepare", "worker", "submit", "collect", "checks"))
    ap.add_argument("--slot", type=int)
    args = ap.parse_args()
    if args.command == "worker":
        worker(args.slot)
    else:
        globals()[args.command]()
