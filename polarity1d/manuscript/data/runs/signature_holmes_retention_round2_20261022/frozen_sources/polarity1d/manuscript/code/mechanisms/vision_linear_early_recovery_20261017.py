"""Documented computational amendment: retain only the required Spring prefix.

Equations, perturbations, clocks, tolerances and every original prefix sample
are unchanged. Compare the adjacent predeclared completed L0=.9 neighbor before
accepting the L0=.8 result. The original t=1 attempt remains preserved.
"""
from pathlib import Path
import argparse
import copy
import json
import os
import shlex
import shutil
import subprocess
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import solve_ivp
import all_five_incumbent_recovery_20261017 as assay

ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT/"data/runs/vision_linear_20261017/recovery"
DEFAULT = PARENT/"early_window_amendment"
TOPIC = ROOT/"figures/length_mechanisms_20261008/dev/vision_linear_20261017/supporting/metadata"
sha, read, write, csvwrite = assay.sha, assay.read, assay.write, assay.csvwrite


def initial_cases():
    parent = assay.verify(PARENT)
    cases = []
    for size in (.8, .9):
        old = next(c for c in parent["cases"] if c["model"] == "spring" and c["size_native"] == size)
        case = copy.deepcopy(old); case.update(slot=len(cases), parent_slot=old["slot"],
                                               amendment_role="required_target" if size == .8 else "predeclared_neighbor_check")
        cases.append(case)
    return parent, cases


def preflight():
    parent, cases = initial_cases(); horizon = parent["fixed_model_horizons"]["spring"]
    times = np.asarray(parent["continuous_samples"])
    prefix = times[times <= horizon]
    assert prefix[0] == 0 and prefix[-1] == horizon
    for case in cases:
        state, fields, x = assay.load_state(PARENT/case["input"], "spring")
        g = assay.initial_geometry(case, state, fields, x)
        for epsilon in assay.EPSILONS: assay.kick(state, "spring", g["direction"], epsilon)
        if case["size_native"] == .9:
            folder = PARENT/"workers"/f'{case["parent_slot"]:02d}'
            status = read(folder/"status.json"); assert status["complete"]
            assert sha(folder/"curves.npz") == status["curves_sha256"]
    return dict(generated_at=assay.now(), passed=True, new_integrations=0, sizes=[.8, .9],
                unchanged_selected_horizon=horizon, original_prefix_samples=prefix.tolist(),
                exact_original_prefix_retained=True,
                neighbor_comparison="Every sampled prefix |F_new-F_old| <= 3*(E_new+E_old); same predeclared empirical-sensitivity allowance",
                later_candidates="unmeasured for the required target, never inferred",
                original_parent_job=read(PARENT/"launch.json")["job_id"])


def initialize(run):
    if run.exists(): raise ValueError("Fresh amendment destination required")
    check = preflight(); parent, cases = initial_cases(); horizon = check["unchanged_selected_horizon"]
    for folder in ("inputs", "workers", "logs", "results", "frozen_sources"): (run/folder).mkdir(parents=True)
    sources = [Path(__file__), Path(assay.__file__), PARENT/"manifest.json", PARENT/"provenance.json",
               TOPIC/"EARLY_WINDOW_AMENDMENT.md"]
    for folder, pattern in (("code/core", "*.py"), ("code/core", "*.json"),
                            ("code/core/models", "*.py"), ("code/core/src/polarity1d", "*.py")):
        sources.extend(sorted((ROOT/folder).glob(pattern)))
    hashes = {}
    for source in sources:
        rel = source.relative_to(ROOT); target = run/"frozen_sources"/rel
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(source, target); hashes[str(rel)] = sha(target)
    for case in cases:
        source = PARENT/case["input"]; target = run/"inputs"/f'case_{case["slot"]:02d}.npz'
        shutil.copy2(source, target); assert sha(target) == case["checkpoint_sha256"]
        case["input"] = str(target.relative_to(run))
    manifest = copy.deepcopy(parent)
    manifest.update(generated_at=assay.now(), experiment="Same fixed Spring recovery window, without unused later dynamics",
                    cases=cases, continuous_samples=check["original_prefix_samples"],
                    continuous_candidates=[h for h in parent["continuous_candidates"] if h <= horizon],
                    unmeasured_later_candidates=[h for h in parent["continuous_candidates"] if h > horizon],
                    integration_end_native=horizon, frozen_source_sha256=hashes,
                    parent_run=str(PARENT), preflight=check,
                    stopping_rule="Two cases only; qualify target and compare adjacent completed neighbor, then stop")
    write(run/"manifest.json", manifest); write(run/"provenance.json", dict(generated_at=assay.now(), manifest_sha256=sha(run/"manifest.json")))
    return manifest


def integrate_prefix(case, state, manifest, reference=False):
    assert case["model"] == "spring"
    settings = manifest["continuous_tighter" if reference else "continuous_tight"]
    shape = state.shape
    fun = lambda t, y: assay.continuous_rhs(case, y.reshape(shape)).ravel()
    solution = solve_ivp(fun, (0., manifest["integration_end_native"]), state.ravel(), method="LSODA",
                         t_eval=np.asarray(manifest["continuous_samples"]), **settings)
    if not solution.success: raise RuntimeError(solution.message)
    states = solution.y.T.reshape((len(solution.t),)+shape)
    if not np.isfinite(states).all(): raise FloatingPointError("Nonfinite native prefix")
    if float(np.min(states[:, :-2])) < -max(50*settings["atol"], 1e-7): raise FloatingPointError("Materially negative native concentration")
    fields = dict(active=states[:, :-2:2], inactive=states[:, 1:-2:2], length=states[:, -1]-states[:, -2],
                  left_endpoint=states[:, -2], right_endpoint=states[:, -1])
    return dict(time=solution.t, fields=fields, states=states,
                numerical=dict(kind="tighter_same_window_prefix" if reference else "tight_same_window_prefix",
                               nfev=int(solution.nfev), njev=int(solution.njev), nlu=int(solution.nlu),
                               minimum_raw_chemical_value=float(np.min(states[:, :-2])), clipping_applied=False,
                               integration_end_native=manifest["integration_end_native"]))


def worker(run, slot):
    original = assay.integrate
    try:
        assay.integrate = integrate_prefix
        assay.worker(run, slot)
    finally:
        assay.integrate = original


def submit(run):
    assay.verify(run)
    if (run/"launch.json").exists(): raise ValueError("Already launched")
    frozen = run/"frozen_sources/code/mechanisms"/Path(__file__).name
    command = shlex.join([sys.executable, str(frozen), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    options = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G", "--time=00:20:00",
               "--job-name=vision-linear-prefix", "--array=0-1", "--chdir="+str(run),
               "--output="+str(run/"logs/prefix_%A_%a.out"), "--error="+str(run/"logs/prefix_%A_%a.err")]
    job = subprocess.check_output(options+["--wrap", command], text=True).strip().split(";")[0]
    write(run/"launch.json", dict(generated_at=assay.now(), job_id=job, command=command, cases=2))
    return job


def collect(run):
    manifest = assay.verify(run); rows = []; curves = []; comparison = []
    horizon = manifest["integration_end_native"]
    for case in manifest["cases"]:
        folder = run/"workers"/f'{case["slot"]:02d}'; status = read(folder/"status.json")
        assert status["complete"] and sha(folder/"curves.npz") == status["curves_sha256"]
        with np.load(folder/"curves.npz") as z:
            a, b = assay.evaluate_case(case, z, manifest); curves.extend(a); rows.extend(b)
            selected = next(r for r in b if r["horizon_native"] == horizon)
            assert selected["passes"], selected
            if case["size_native"] == .9:
                parent_folder = Path(manifest["parent_run"])/"workers"/f'{case["parent_slot"]:02d}'
                old_status = read(parent_folder/"status.json")
                assert sha(parent_folder/"curves.npz") == old_status["curves_sha256"]
                with np.load(parent_folder/"curves.npz") as old:
                    old_curves, _ = assay.evaluate_case(case, old, manifest)
                for epsilon in assay.EPSILONS[1:]:
                    for row in [r for r in a if r["epsilon"] == epsilon]:
                        prior = next(r for r in old_curves if r["epsilon"] == epsilon and r["time_native"] == row["time_native"])
                        difference = abs(row["restored_fraction"]-prior["restored_fraction"])
                        allowance = 3*(row["empirical_numerical_sensitivity"]+prior["empirical_numerical_sensitivity"])
                        comparison.append(dict(epsilon=epsilon, time_native=row["time_native"], absolute_F_difference=difference,
                                               declared_sensitivity_allowance=allowance, passes=difference <= allowance))
    assert all(r["passes"] for r in comparison), "Neighbor prefix differs beyond declared empirical sensitivity allowance"
    csvwrite(run/"results/candidate_gates.csv", rows); csvwrite(run/"results/short_curves.csv", curves)
    csvwrite(run/"results/neighbor_prefix_comparison.csv", comparison)
    write(run/"results/verification.json", dict(generated_at=assay.now(), passed=True, cases=2,
          exact_clock_and_prefix_preserved=True, target_and_neighbor_clock_gates_pass=True,
          neighbor_prefix_comparison_pass=True,
          max_neighbor_absolute_F_difference=max(r["absolute_F_difference"] for r in comparison),
          unmeasured_later_candidates=manifest["unmeasured_later_candidates"],
          source_sha256=sha(Path(__file__))))
    write(run/"results/status.json", dict(generated_at=assay.now(), complete=True, producer="finished", passed=True))
    return [r for r in rows if r["horizon_native"] == horizon]


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("action", choices=("preflight", "init", "submit", "worker", "collect"))
    parser.add_argument("--run", type=Path, default=DEFAULT); parser.add_argument("--slot", type=int)
    args = parser.parse_args(); run = args.run.resolve()
    if args.action == "preflight":
        result = preflight(); write(TOPIC/"EARLY_WINDOW_PREFLIGHT.json", result); print(json.dumps(result, indent=2))
    elif args.action == "init": initialize(run)
    elif args.action == "submit": print(submit(run))
    elif args.action == "worker": worker(run, args.slot)
    else: print(json.dumps(collect(run), indent=2))


if __name__ == "__main__": main()
