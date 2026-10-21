"""Fresh Figure 3 single-region basals; freeze-only launch payload by default.

Never reuses old endpoints, retunes kinetics or searches preparation duration.
The fixed native preparation is qualified before the separate long no-cue
control, and the challenge checkpoint remains the exact qualified pre-control
state. Full integrations are restricted to explicitly launched SLURM workers.
"""
from pathlib import Path
import argparse
import os
import shlex
import shutil
import sys
import time

import numpy as np
import signature_figure3_native_20261020 as native
from signature_figure3_geometry_20261020 import architecture

ROOT = native.ROOT
DEFAULT = ROOT/"data/runs/signature_figure3_basals_20261020"
read, write, sha, now = native.read, native.write, native.sha, native.now


def sources():
    paths = [Path(__file__)]+[Path(__file__).with_name(name) for name in (
        "signature_figure3_native_20261020.py", "signature_figure3_geometry_20261020.py",
        "ring_distance_prototype_20261017.py", "length_dev_output_layout.py",
        "spring_ring_model_20261019.py", "spring_ring_closure_20261017.py")]
    for directory in (ROOT/"code/core/models", ROOT/"code/core/src/polarity1d"):
        paths += sorted(directory.glob("*.py"))
    return paths


def freeze(run):
    if run.exists():
        raise FileExistsError("Choose a fresh directory; frozen sources and retained evidence are immutable")
    for directory in ("inputs", "frozen_sources", "workers", "logs", "results"):
        (run/directory).mkdir(parents=True)
    hashes = {}
    for path in sources():
        relative = path.relative_to(ROOT)
        target = run/"frozen_sources"/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[str(relative)] = sha(target)
    jobs, input_hashes = [], {}
    for model in native.MODELS:
        for size in native.SIZES[model]:
            adapter = native.Adapter(model, size)
            state = adapter.initial()
            path = run/"inputs"/f"{model}_C{size:g}_native_homogeneous.npz"
            np.savez_compressed(path, state=state, material_positions=adapter.q,
                                cells=adapter.cells, size_parameter=size, model=model)
            input_hashes[str(path.relative_to(run))] = sha(path)
            jobs.append(dict(slot=len(jobs), model=model, size_parameter=size, cells=adapter.cells,
                reference_density=100 if adapter.spring else 20,
                nominal_spacing=.01 if adapter.spring else .05,
                input=str(path.relative_to(run)), input_sha256=sha(path),
                preparation_peak=native.PREP_PEAKS[model], preparation_duration=10.,
                preparation_center_fraction=float(adapter.q[0]), cue_support=adapter.support,
                no_cue_control_horizon=5300. if adapter.spring else 53000.,
                no_cue_control_reason="Longest possible complete A2 protocol: 2*t_handover_cap + recovery",
                fixed_scales=adapter.scales.tolist(),
                native_initial_amounts=adapter.amounts(state).tolist()))
    manifest = dict(generated_at=now(), experiment="Figure 3 fresh three-size ring basals",
        jobs=jobs, requested_core_basals=12, existing_checkpoints_reused=0,
        native_kinetics_changed=False, source_sha256=hashes, input_sha256=input_hashes,
        solver=dict(method="BDF", rtol=1e-9, atol=1e-11, max_step=10.),
        rule=dict(window=1000., consecutive_windows=2, window_tolerance=1e-4,
            holdout=10000., holdout_tolerance=1e-3, maximum_qualification_release=30000.,
            final_window_tolerance=1e-4, minimum_contrast=.05,
            resolved_peak_prominence=.05, one_half_contrast_region=True, one_resolved_peak=True,
            minimum_concentration=-5e-10, conserved_amount_relative_tolerance=1e-7),
        long_no_cue_gate=dict(maximum_scaled_excursion=1e-3, final1000_scaled_range=1e-4,
                              exact_rhs_projected1000=1e-4, same_architecture=True),
        preparation_search=False,
        failure_policy="Report unqualified preparation; do not tune or replace checkpoint. A largest-size failure permits two-size model analysis.",
        spring_geometry="Existing mean-force global-extension circular analogue; sizes label resting circumference; physical circumference remains dynamic and is saved.",
        evidence_status="Prospective; no qualified basals or figure conclusions are asserted by freezing",
        time_budget=dict(worker_soft_seconds=3300, slurm_wall_seconds=3600))
    write(run/"manifest.json", manifest)
    write(run/"freeze.json", dict(generated_at=now(), manifest_sha256=sha(run/"manifest.json")))
    runner = run/"frozen_sources"/Path(__file__).relative_to(ROOT)
    command = "exec "+shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1",
        "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1", native.PYTHON, str(runner), "worker", "--run", str(run)])
    command += ' --slot "$SLURM_ARRAY_TASK_ID"'
    arguments = ["sbatch", "--parsable", "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=3G",
        "--time=01:00:00", "--signal=B:USR1@60", "--array=0-11", "--job-name=figure3-fresh-basals",
        "--output="+str(run/"logs/%A_%a.out"), "--error="+str(run/"logs/%A_%a.err"), "--wrap", command]
    write(run/"slurm_payload.json", dict(generated_at=now(), submitted=False,
        sbatch_arguments=arguments, manifest_sha256=sha(run/"manifest.json"),
        resource_reason="One native BDF solve per task; sparse chemical and existing dense Spring Jacobians; 3 GiB covers full-field preparation and long control histories.",
        no_concurrency_throttle=True))
    return preflight(run)


def verify(run):
    manifest = read(run/"manifest.json")
    assert sha(run/"manifest.json") == read(run/"freeze.json")["manifest_sha256"]
    for relative, digest in manifest["source_sha256"].items():
        assert sha(run/"frozen_sources"/relative) == digest, relative
    for relative, digest in manifest["input_sha256"].items():
        assert sha(run/relative) == digest, relative
    return manifest


def preflight(run):
    manifest = verify(run)
    assert len(manifest["jobs"]) == 12
    checks = []
    for job in manifest["jobs"]:
        adapter = native.Adapter(job["model"], job["size_parameter"], job["cells"])
        with np.load(run/job["input"]) as z:
            state = z["state"].copy()
            np.testing.assert_array_equal(z["material_positions"], adapter.q)
        np.testing.assert_array_equal(state, adapter.initial())
        fields = state[:-1].reshape(2, -1) if adapter.spring else state.reshape(-1, adapter.cells)
        assert np.all(np.ptp(fields, axis=1) == 0)
        assert job["cells"] == int(round(job["size_parameter"]*job["reference_density"]))
        rhs = adapter.rhs(state, job["preparation_peak"], job["preparation_center_fraction"])
        assert np.isfinite(rhs).all()
        if adapter.spring:
            d_amount = adapter.circumference(state)*np.mean(rhs[:adapter.cells]+rhs[adapter.cells:-1])
            d_amount += rhs[-1]*np.mean(state[:adapter.cells]+state[adapter.cells:-1])
            conservation_error = abs(float(d_amount))
        else:
            conservation_error = float(np.abs(adapter.amounts(rhs)).max())
        assert conservation_error < 1e-10
        checks.append(dict(model=job["model"], size=job["size_parameter"], cells=job["cells"],
                           fresh_homogeneous=True, rhs_conservation_absolute_error=conservation_error,
                           long_no_cue_horizon=job["no_cue_control_horizon"]))
    result = dict(generated_at=now(), passed=True, submitted=False, native_integrations=0,
                  source_files=len(manifest["source_sha256"]), input_hashes_verified=12,
                  manifest_sha256=sha(run/"manifest.json"), checks=checks,
                  qualification="Pending actual SLURM trajectories; this is algebraic/input validation")
    write(run/"results/preflight.json", result)
    return result


def worker(run, slot):
    native.require_slurm()
    manifest = verify(run)
    job = manifest["jobs"][slot]
    adapter = native.Adapter(job["model"], job["size_parameter"], job["cells"])
    directory = run/"workers"/f"{slot:03d}"
    directory.mkdir()
    record = dict(generated_at=now(), producer="running", complete=False, status="preparing",
                  job=job, slurm_job_id=os.environ["SLURM_JOB_ID"], phases=[],
                  fixed_preparation_duration=10., preparation_duration_search=False)
    path = directory/"status.json"
    write(path, record)
    deadline = time.monotonic()+manifest["time_budget"]["worker_soft_seconds"]
    with np.load(run/job["input"]) as z:
        state = z["state"].copy()

    def phase(label, initial, duration, peak=0., interval=5.):
        record.update(generated_at=now(), status=label)
        write(path, record)
        result = native.integrate(adapter, initial, duration, peak,
            job["preparation_center_fraction"], deadline=deadline, interval=interval)
        archive = directory/f'{len(record["phases"]):03d}_{label}.npz'
        metadata = native.save_phase(archive, result, adapter, peak, job["preparation_center_fraction"])
        metadata["kind"] = label
        metadata["path"] = str(archive.relative_to(run))
        record["phases"].append(metadata)
        write(path, record)
        if result["status"] != "complete":
            raise RuntimeError("Incomplete native phase: "+result["status"])
        return result

    try:
        result = phase("preparation", state, 10., job["preparation_peak"], .1)
        state, consecutive, elapsed, qualified = result["states"][-1].copy(), 0, 0., False
        while elapsed+12000 <= 30000 or (consecutive > 0 and elapsed+11000 <= 30000):
            result = phase("basal_window", state, 1000.)
            state = result["states"][-1].copy()
            elapsed += 1000.
            metric = native.window_metric(result["time"], result["states"], adapter.scales)
            record["phases"][-1]["diagnostic"] = metric
            consecutive = consecutive+1 if metric["passed"] else 0
            if consecutive < 2:
                continue
            initial_holdout = state.copy()
            holdout = phase("qualification_holdout", state, 10000.)
            elapsed += 10000.
            state = holdout["states"][-1].copy()
            excursion = float(np.max(np.abs(holdout["states"]-initial_holdout)/adapter.scales))
            tail = native.window_metric(holdout["time"], holdout["states"], adapter.scales)
            features = architecture(adapter.primary(state))
            passed = excursion <= 1e-3 and tail["passed"]
            record["qualification"] = dict(two_consecutive_windows=consecutive,
                scaled_holdout_excursion=excursion, final_window=tail, architecture=features,
                release_elapsed=elapsed, settled=passed)
            if passed:
                qualified = features["contrast"] >= .05 and features["half_contrast_regions"] == features["resolved_peaks"] == 1
                break
            consecutive = 0
        if not qualified:
            record.update(status="unqualified_single_region_basal", seed_status="fixed_preparation_unqualified",
                          final_architecture=architecture(adapter.primary(state)))
        else:
            checkpoint = directory/"qualified_pre_control_basal.npz"
            np.savez_compressed(checkpoint, state=state, material_positions=adapter.q,
                model=adapter.key, size_parameter=adapter.size, cells=adapter.cells,
                physical_circumference=adapter.circumference(state), fixed_scales=adapter.scales)
            record.update(checkpoint=str(checkpoint.relative_to(run)), checkpoint_sha256=sha(checkpoint),
                          origin="fresh_homogeneous_native_preparation", seed_status="qualified_fixed_duration_seed")
            control_initial = state.copy()
            control = phase("long_no_cue_control", state, job["no_cue_control_horizon"], interval=20.)
            excursion = float(np.max(np.abs(control["states"]-control_initial)/adapter.scales))
            tail = native.window_metric(control["time"], control["states"], adapter.scales)
            rhs = float(np.max(np.abs(adapter.rhs(control["states"][-1]))/adapter.scales)*1000.)
            features = architecture(adapter.primary(control["states"][-1]))
            passed = excursion <= 1e-3 and tail["passed"] and rhs <= 1e-4
            passed &= features["contrast"] >= .05 and features["half_contrast_regions"] == features["resolved_peaks"] == 1
            record.update(status="qualified" if passed else "long_no_cue_control_failed",
                long_no_cue=dict(duration=job["no_cue_control_horizon"], max_scaled_excursion=excursion,
                    final_window=tail, exact_rhs_projected1000=rhs, final_architecture=features, passed=bool(passed)),
                challenge_checkpoint_is_exact_pre_control_state=True)
    except Exception as error:
        record.update(status="error", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        record.update(generated_at=now(), producer="finished", complete=True)
        write(path, record)


def audit(run):
    manifest = verify(run)
    rows = []
    for job in manifest["jobs"]:
        path = run/"workers"/f'{job["slot"]:03d}'/"status.json"
        if not path.exists():
            rows.append(dict(job=job, complete=False, status="not_started")); continue
        record = read(path)
        if not record.get("complete"):
            rows.append(record); continue
        adapter = native.Adapter(job["model"], job["size_parameter"], job["cells"])
        with np.load(run/job["input"]) as z:
            previous = z["state"].copy()
        consecutive, maximum_error = 0, 0.
        recomputed_qualification, recomputed_long_control = False, False
        for phase in record["phases"]:
            archive = run/phase["path"]
            assert sha(archive) == phase["sha256"]
            with np.load(archive) as z:
                tt, yy = z["time"].copy(), z["states"].copy()
            np.testing.assert_array_equal(yy[0], previous)
            assert np.isfinite(yy).all() and np.all(np.diff(tt) > 0)
            assert min(adapter.concentration_minimum(y) for y in yy) >= -5e-10
            amounts = np.asarray([adapter.amounts(y) for y in yy])
            maximum_error = max(maximum_error, float(np.max(np.abs(amounts/np.asarray(job["native_initial_amounts"])-1))))
            assert maximum_error <= 1e-7
            if phase["kind"] == "basal_window":
                metric = native.window_metric(tt, yy, adapter.scales)
                consecutive = consecutive+1 if metric["passed"] else 0
            if phase["kind"] == "qualification_holdout":
                excursion = float(np.max(np.abs(yy-yy[0])/adapter.scales))
                tail = native.window_metric(tt, yy, adapter.scales)
                features = architecture(adapter.primary(yy[-1]))
                recomputed_qualification = consecutive >= 2 and excursion <= 1e-3 and tail["passed"] and features["contrast"] >= .05 and features["half_contrast_regions"] == features["resolved_peaks"] == 1
                if not recomputed_qualification:
                    consecutive = 0
            if phase["kind"] == "long_no_cue_control":
                assert tt[-1] == job["no_cue_control_horizon"]
                with np.load(run/record["checkpoint"]) as z:
                    np.testing.assert_array_equal(yy[0], z["state"])
                excursion = float(np.max(np.abs(yy-yy[0])/adapter.scales))
                tail = native.window_metric(tt, yy, adapter.scales)
                rhs = float(np.max(np.abs(adapter.rhs(yy[-1]))/adapter.scales)*1000)
                shape = architecture(adapter.primary(yy[-1]))
                recomputed_long_control = excursion <= 1e-3 and tail["passed"] and rhs <= 1e-4 and shape["contrast"] >= .05 and shape["half_contrast_regions"] == shape["resolved_peaks"] == 1
            previous = yy[-1].copy()
        if record["status"] == "qualified":
            assert sha(run/record["checkpoint"]) == record["checkpoint_sha256"]
            assert recomputed_qualification and recomputed_long_control
        record["independent_saved_field_audit"] = dict(qualification_recomputed=bool(recomputed_qualification),
            long_control_recomputed=bool(recomputed_long_control), max_saved_amount_error=maximum_error,
            hashes_and_full_native_phase_continuity_verified=True)
        rows.append(record)
    result = dict(generated_at=now(), complete=all(r.get("complete", False) for r in rows),
        qualified=sum(r["status"] == "qualified" for r in rows), requested=12, records=rows,
        native_trajectory_replay=False, evidence_scope="Independent recomputation from exact saved native fields")
    write(run/"results/basal_audit.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "preflight", "worker", "audit"))
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    result = worker(args.run.resolve(), args.slot) if args.action == "worker" else globals()[args.action](args.run.resolve())
    if result is not None:
        print(__import__("json").dumps(result if args.action != "audit" else
            {k: result[k] for k in ("generated_at", "complete", "qualified", "requested")}, indent=2))


if __name__ == "__main__":
    main()
