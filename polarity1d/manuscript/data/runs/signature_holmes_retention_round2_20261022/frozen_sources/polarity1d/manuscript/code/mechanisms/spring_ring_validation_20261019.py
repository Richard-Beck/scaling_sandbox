"""Independent bounded numerical validation of the declared Spring ring analogue.

Fresh evidence only. Integrations require SLURM; retained 17 October outputs
are read without mutation. No chemical clipping or amount renormalization.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from spring_ring_model_20261019 import Closure, physical_raised_cosine, CuePhase, run_protocol

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / "data/runs/spring_ring_validation_20261019"
PREP = ROOT / "data/runs/spring_ring_preparation_pilot_20261017"
OLD = ROOT / "data/runs/spring_ring_distance_pilot_20261017"
PYTHON = "/home/4473331/.conda/envs/scaling_sandbox/bin/python"


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
        raise ValueError("Fresh run directory required; retained evidence is immutable")
    for directory in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run / directory).mkdir(parents=True)
    paths = [Path(__file__), Path(__file__).with_name("spring_ring_model_20261019.py"),
             Path(__file__).with_name("spring_ring_closure_20261017.py"),
             Path(__file__).with_name("spring_ring_model_checks_20261019.py"),
             ROOT / "code/core/models/spring.py"]
    hashes = {}
    for source in paths:
        relative = source.relative_to(ROOT)
        target = run / "frozen_sources" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        hashes[str(relative)] = sha(target)
    assert read(PREP / "results/audit.json")["passed"]
    shutil.copy2(PREP / "workers/qualified_basal.npz", run / "inputs/basal.npz")
    jobs = []
    for name, n, fresh, method, rtol, shift, perturb, duration in (
        ("fresh_200", 200, True, "BDF", 1e-9, 0, False, 1000.),
        ("fresh_400", 400, True, "BDF", 1e-9, 0, False, 1000.),
        ("fresh_200_tight", 200, True, "BDF", 1e-11, 0, False, 1000.),
        ("retained_radau", 200, False, "Radau", 1e-9, 0, False, 1000.),
        ("retained_tight", 200, False, "BDF", 1e-11, 0, False, 1000.),
        ("rotation_73", 200, False, "BDF", 1e-9, 73, False, 1000.),
        ("conserved_perturbation", 200, False, "BDF", 1e-9, 0, True, 0.),
        ("finite_0p1", 200, False, "BDF", 1e-9, 0, False, .1),
        ("finite_10", 200, False, "BDF", 1e-9, 0, False, 10.),
    ):
        jobs.append(dict(slot=len(jobs), name=name, cells=n, fresh_preparation=fresh,
                         method=method, rtol=rtol, atol=rtol / 100, shift_cells=shift,
                         perturb=perturb, cue_duration=duration))
    manifest = dict(generated_at=now(), scope="C0=2 bounded independent numerical validation",
        model="Spring circular global-extension analogue", resting_circumference=2.,
        preparation_duration=10., basal_duration=13000., recovery_duration=13000.,
        cue_peak=4., cue_physical_support=.2, challenger_offset_fraction=.4,
        solver_max_step=10., jobs=jobs, source_sha256=hashes,
        input_sha256=sha(run / "inputs/basal.npz"),
        original_preparation_sha256=sha(PREP / "workers/qualified_basal.npz"),
        original_distance_manifest_sha256=sha(OLD / "manifest.json"),
        amount_relative_tolerance=1e-7, negative_tolerance=5e-10,
        tail_full_field_scaled_tolerance=1e-4, stationary_rhs_scaled_tolerance=1e-4,
        spatial_comparison_rule=dict(max_field_difference=.005,
            max_circumference_difference=.0002, max_vector_magnitude_difference=.005,
            reason="Bounded sensitivity thresholds on native fields; not a mathematical convergence proof"),
        perturbation=dict(seed=20261019, amplitude=.001,
            rule="Independent seeded bounded active increments; opposite inactive increments preserve each cell's total concentration exactly"),
        interpretations=["No native 2D mechanical equivalence", "No size law or noisy directional persistence",
                         "Fixed local kinetic coefficients; force averaging remains the declared new closure"])
    write(run / "manifest.json", manifest)
    runner = run / "frozen_sources" / Path(__file__).relative_to(ROOT)
    command = shlex.join([PYTHON, str(runner), "worker", "--run", str(run)])
    command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    job = subprocess.check_output(["sbatch", "--parsable", "--qos=small", "--ntasks=1",
        "--cpus-per-task=1", "--mem=3G", "--time=00:30:00", "--array=0-8",
        "--job-name=spring-validate", "--output=" + str(run / "logs/%A_%a.out"),
        "--error=" + str(run / "logs/%A_%a.err"), "--wrap", "exec " + command], text=True).strip()
    write(run / "launch.json", dict(generated_at=now(), array_job_id=job,
                                     manifest_sha256=sha(run / "manifest.json")))
    collect(run)
    print(json.dumps(dict(array_job_id=job, run=str(run), jobs=jobs), indent=2))


def sample_times(end):
    return np.unique(np.r_[np.arange(0, min(10., end), .1),
        np.arange(10., min(100., end), 1.), np.arange(100., end, 5.), end])


def features(state, n):
    u, v, circumference = state[:n], state[n:-1], float(state[-1])
    harmonic = np.mean(u * np.exp(2j * np.pi * np.arange(n) / n)) / np.mean(u)
    contrast = float(np.ptp(u))
    mask = u >= u.min() + .5 * contrast
    regions = 0 if contrast < 1e-12 else (1 if mask.all() else int(np.sum(mask & ~np.roll(mask, 1))))
    peaks, props = find_peaks(np.tile(u, 3), prominence=.05)
    peaks = peaks[(peaks >= n) & (peaks < 2*n)] - n
    return dict(circumference=circumference, active_minimum=float(u.min()),
        active_peak=float(u.max()), contrast=contrast, half_contrast_regions=regions,
        resolved_peaks=len(peaks), peak_cells=peaks.tolist(),
        vector_magnitude=float(abs(harmonic)),
        vector_direction_fraction=float(np.angle(harmonic) / (2*np.pi) % 1),
        amount=float(circumference * np.mean(u + v)))


def worker(run, slot):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Trajectory integration requires SLURM")
    manifest = read(run / "manifest.json")
    for relative, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    assert sha(run / "inputs/basal.npz") == manifest["input_sha256"]
    job = manifest["jobs"][slot]
    directory = run / "workers" / f"{slot:02d}"
    directory.mkdir()
    record = dict(generated_at=now(), producer="running", status="initializing",
                  complete=False, job=job, slurm_job_id=os.environ["SLURM_JOB_ID"], phases=[])
    record_path = directory / "status.json"
    write(record_path, record)
    closure = Closure(2., job["cells"])
    n = closure.cells
    if job["fresh_preparation"]:
        state = closure.initial()
    else:
        with np.load(run / "inputs/basal.npz") as z:
            state = z["state"].copy()
    if job["shift_cells"]:
        shift = job["shift_cells"]
        state = np.r_[np.roll(state[:n], shift), np.roll(state[n:-1], shift), state[-1]]
    if job["perturb"]:
        change = np.random.default_rng(20261019).uniform(-.001, .001, n)
        state[:n] += change
        state[n:-1] -= change
        np.savez_compressed(directory / "perturbed_initial.npz", state=state)
    started = time.monotonic()

    def evolve(kind, duration, peak, center):
        nonlocal state
        record.update(generated_at=now(), status=kind)
        write(record_path, record)
        def cue(y):
            return physical_raised_cosine(closure, y, center, peak, .2)
        def fun(t, y):
            return closure.rhs(y, cue(y)[0])
        def jac(t, y):
            return closure.jacobian(y, *cue(y))
        t_eval = sample_times(duration)
        solution = solve_ivp(fun, (0., duration), state, method=job["method"],
            rtol=job["rtol"], atol=job["atol"], max_step=10., t_eval=t_eval, jac=jac)
        states = solution.y.T
        amounts = states[:, -1] * np.mean(states[:, :n] + states[:, n:-1], axis=1)
        scales = np.r_[np.full(2*n, 2.), 2.]
        tail = solution.t >= max(0., duration - 1000.)
        tail_range = float(np.max(np.ptp(states[tail], axis=0) / scales))
        final_rhs = float(np.max(np.abs(fun(duration, states[-1])) / scales) * 1000.)
        archive = directory / (kind + ".npz")
        np.savez_compressed(archive, time=solution.t, states=states,
            material_positions=closure.material_positions, cue_peak=np.asarray(peak),
            center_fraction=np.asarray(center), total_physical_support=np.asarray(.2))
        phase = dict(kind=kind, status="complete" if solution.success else "numerical_error",
            message=solution.message, duration=duration, end_time=float(solution.t[-1]),
            path=str(archive.relative_to(run)), sha256=sha(archive), nfev=int(solution.nfev),
            njev=int(solution.njev), nlu=int(solution.nlu),
            max_sampled_amount_relative_error=float(np.max(np.abs(amounts / 4. - 1))),
            minimum_concentration=float(states[:, :-1].min()), all_finite=bool(np.isfinite(states).all()),
            final_1000_scaled_range=tail_range, final_rhs_scaled_to_1000=final_rhs,
            final=features(states[-1], n))
        record["phases"].append(phase)
        write(record_path, record)
        if not solution.success or len(solution.t) != len(t_eval):
            raise RuntimeError(solution.message)
        if not phase["all_finite"] or phase["minimum_concentration"] < -5e-10 or phase["max_sampled_amount_relative_error"] > 1e-7:
            raise RuntimeError("Finite/nonnegative/amount validation failed")
        state = states[-1].copy()

    try:
        center = float(closure.material_positions[job["shift_cells"] % n])
        if job["fresh_preparation"]:
            evolve("preparation", 10., 4., center)
            evolve("basal", 13000., 0., center)
        if job["cue_duration"]:
            evolve("challenge", job["cue_duration"], 4., (center + .4) % 1.)
            evolve("recovery", 13000., 0., (center + .4) % 1.)
        else:
            evolve("basal_perturbed", 13000., 0., center)
        final = record["phases"][-1]
        if job["name"] == "finite_0p1":
            with np.load(run / "inputs/basal.npz") as z:
                smoke_initial = z["state"].copy()
            smoke = run_protocol(closure,
                [CuePhase(.1, 4., center + .4, .2, .1), CuePhase(1., 0., center + .4, .2, .5)],
                initial_state=smoke_initial)
            _, yc = load_phase(run, record, "challenge")
            tr, yr = load_phase(run, record, "recovery")
            reference = np.vstack([smoke_initial, yc[-1],
                yr[np.flatnonzero(np.isclose(tr, .5))[0]],
                yr[np.flatnonzero(np.isclose(tr, 1.))[0]]])
            residual = float(np.max(np.abs(smoke["states"] - reference)))
            assert residual < 1e-8
            record["reusable_protocol_smoke"] = dict(passed=True,
                max_complete_state_difference=residual,
                exact_time=[float(t) for t in smoke["time"]],
                max_saved_amount_relative_error=smoke["max_saved_amount_relative_error"])
        record.update(status="complete", final=final["final"],
            tail_stationary=bool(final["final_1000_scaled_range"] <= 1e-4 and
                                 final["final_rhs_scaled_to_1000"] <= 1e-4))
    except Exception as error:
        record.update(status="error", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        record.update(generated_at=now(), producer="finished", complete=True,
                      elapsed_seconds=time.monotonic() - started)
        write(record_path, record)
    print(json.dumps(record, indent=2))


def load_phase(run, record, kind):
    phase = next(item for item in record["phases"] if item["kind"] == kind)
    path = run / phase["path"]
    assert sha(path) == phase["sha256"]
    with np.load(path) as z:
        return z["time"].copy(), z["states"].copy()


def periodic_resample(states, source_n, target_n):
    # All prepared poles lie at the first cell; compare positions relative to it.
    source = np.arange(source_n + 1) / source_n
    target = np.arange(target_n) / target_n
    fields = []
    for start in (0, source_n):
        block = states[:, start:start+source_n]
        closed = np.c_[block, block[:, 0]]
        fields.append(np.asarray([np.interp(target, source, row) for row in closed]))
    return np.c_[fields[0], fields[1], states[:, -1]]


def compare(a, b, n):
    difference = np.abs(a - b)
    return dict(max_chemical_field_difference=float(difference[:, :-1].max()),
        max_circumference_difference=float(difference[:, -1].max()),
        final_chemical_field_difference=float(difference[-1, :-1].max()),
        final_circumference_difference=float(difference[-1, -1]),
        max_vector_magnitude_difference=float(max(abs(features(x, n)["vector_magnitude"] -
            features(y, n)["vector_magnitude"]) for x, y in zip(a, b))))


def collect(run):
    manifest, launch = read(run / "manifest.json"), read(run / "launch.json")
    assert sha(run / "manifest.json") == launch["manifest_sha256"]
    for relative, digest in manifest["source_sha256"].items():
        assert sha(run / "frozen_sources" / relative) == digest
    records = []
    for job in manifest["jobs"]:
        path = run / "workers" / f'{job["slot"]:02d}' / "status.json"
        records.append(read(path) if path.exists() else dict(job=job, producer="pending", status="pending", complete=False))
    result = dict(generated_at=now(), producer="finished" if all(r["complete"] for r in records) else "running",
        status="pending", complete=all(r["complete"] for r in records),
        array_job_id=launch["array_job_id"], jobs=records, comparisons={},
        vector_definition="sum u_j exp(2πij/N) / sum u_j; direction is angle/(2π), arbitrary material origin at cell 0",
        limitations=manifest["interpretations"])
    by_name = {r["job"]["name"]: r for r in records if r.get("status") == "complete"}
    if "fresh_200" in by_name:
        for other, label in (("fresh_400", "grid_200_vs_400"), ("fresh_200_tight", "fresh_tolerance")):
            if other not in by_name:
                continue
            comparisons = {}
            for phase in ("preparation", "basal", "challenge", "recovery"):
                ta, ya = load_phase(run, by_name["fresh_200"], phase)
                tb, yb = load_phase(run, by_name[other], phase)
                np.testing.assert_array_equal(ta, tb)
                if other == "fresh_400":
                    yb = periodic_resample(yb, 400, 200)
                comparisons[phase] = compare(ya, yb, 200)
            result["comparisons"][label] = comparisons
    for name in ("retained_radau", "retained_tight", "rotation_73"):
        if name not in by_name:
            continue
        comparisons = {}
        for phase, oldphase in (("challenge", "maintained_cue"), ("recovery", "recovery")):
            t, y = load_phase(run, by_name[name], phase)
            if name == "rotation_73":
                y = np.c_[np.roll(y[:, :200], -73, axis=1), np.roll(y[:, 200:-1], -73, axis=1), y[:, -1]]
            with np.load(OLD / "workers/01" / (oldphase + ".npz")) as z:
                np.testing.assert_allclose(t, z["time"], rtol=0, atol=1e-12)
                reference = z["states"].copy()
            comparisons[phase] = compare(reference, y, 200)
        result["comparisons"][name + "_vs_retained"] = comparisons
    if "conserved_perturbation" in by_name:
        _, y = load_phase(run, by_name["conserved_perturbation"], "basal_perturbed")
        with np.load(run / "inputs/basal.npz") as z:
            reference = z["state"].copy()
        result["comparisons"]["conserved_perturbation"] = compare(reference[None], y[-1:], 200)
    if result["complete"]:
        errors = [r["job"]["name"] for r in records if r["status"] != "complete"]
        result["execution_errors"] = errors
        result["status"] = "execution_errors" if errors else "complete_requires_review"
        audit_path = run / "results/audit.json"
        if not errors and audit_path.exists():
            audit = read(audit_path)
            assert audit["manifest_sha256"] == sha(run / "manifest.json")
            result["audit"] = dict(path="results/audit.json", sha256=sha(audit_path),
                generated_at=audit["generated_at"], passed=audit["passed"])
            result["status"] = audit["status"]
    write(run / "results/summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("launch", "worker", "collect"))
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    if args.action == "launch":
        launch(args.run.resolve())
    elif args.action == "worker":
        worker(args.run.resolve(), args.slot)
    else:
        result = collect(args.run.resolve())
        print(json.dumps({k: v for k, v in result.items() if k != "jobs"}, indent=2))
