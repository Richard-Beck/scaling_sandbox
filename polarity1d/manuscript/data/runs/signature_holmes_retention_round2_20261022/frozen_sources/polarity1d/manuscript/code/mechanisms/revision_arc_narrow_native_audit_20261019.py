"""Terminal native validation of the narrow-size dense-initial-arc grid; no integration."""
from pathlib import Path
import json
import os
import subprocess
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import revision_arc_narrow_grid_20261019 as producer


def audit():
    run = producer.DEFAULT
    manifest = producer.read(run / "manifest.json")
    frozen = producer.read(run / "freeze.json")
    launch = producer.read(run / "launch.json")
    assert producer.sha(run / "manifest.json") == frozen["manifest_sha256"] == launch["manifest_sha256"]
    requested=len(manifest["cases"]); new_native=len(manifest["jobs"]); exact_reuses=requested-new_native
    inventory = pd.read_csv(run / "results/condition_inventory.csv")
    assert len(inventory) == requested and (inventory.status == "complete").all()
    for rel, digest in manifest["source_sha256"].items():
        assert producer.sha(run / "frozen_sources" / rel) == digest
    for rel, digest in manifest["input_sha256"].items():
        assert producer.sha(run / rel) == digest
    scheduler = subprocess.check_output(["sacct", "-n", "-X", "-j", launch["array_job_id"],
        "--format=JobIDRaw,JobID,State,ExitCode", "-P"], text=True)
    records = [line.split("|") for line in scheduler.splitlines() if line.strip()]
    assert len(records) == new_native, records
    assert all(r[2] == "COMPLETED" and r[3] == "0:0" for r in records), records
    cases, frames, full_phases, compact_phases = [], 0, 0, 0
    minimum, maximum_pool_error = float("inf"), 0.
    checked_source_paths = {}
    for case in manifest["cases"]:
        cid = case["condition_id"]
        row = inventory[inventory.condition_id == cid].iloc[0]
        spring, n, Cb = case["model"] == "spring", case["cells"], case["circumference"]
        with np.load(run / case["input"]) as z:
            initial = z["state"].copy()
        with np.load(run / case["cue_input"]) as z:
            q, signal = z["material_positions"].copy(), z["initial_cue"].copy()
        if spring:
            closure = producer.sizes.Closure(case["resting_circumference"], n)
            expected, expected_dC = producer.sizes.physical_raised_cosine(closure, initial, case["cue_centre_fraction"],
                case["requested_peak"], case["full_physical_support"])
            np.testing.assert_array_equal(expected, signal)
            with np.load(run / case["cue_input"]) as z:
                np.testing.assert_array_equal(expected_dC, z["initial_cue_dC"])
            amount0 = Cb * (initial[:n] + initial[n:-1]).mean()
        else:
            model = producer.base.ring.native_model(case["model"])
            grid = producer.base.ring.RingGrid(Cb, n)
            amount0 = producer.base.ring.pool_amounts(model, grid, initial)
        endpoints, case_frames, case_min, case_error = [], 0, float("inf"), 0.
        for phase, horizon in (("cue", case["cue_duration"]), ("recovery", case["recovery"])):
            path = Path(row["actual_" + phase + "_path"])
            assert producer.sha(path) == row["actual_" + phase + "_sha256"]
            origin = producer.base.source_run(path)
            old = producer.read(origin / "manifest.json")
            assert old["solver"] == manifest["solver"]
            relevant = ["code/mechanisms/spring_ring_model_20261019.py", "code/mechanisms/spring_ring_closure_20261017.py", "code/core/models/spring.py"] if spring else [
                "code/mechanisms/ring_distance_prototype_20261017.py", f'code/core/models/{case["model"]}.py'] + [
                    p for p in manifest["source_sha256"] if p.startswith("code/core/src/polarity1d/")]
            for rel in relevant:
                source = origin / "frozen_sources" / rel
                if str(source) not in checked_source_paths:
                    digest = producer.sha(source)
                    assert digest == manifest["source_sha256"][rel]
                    checked_source_paths[str(source)] = digest
            with np.load(path) as z:
                time = z["time"]
                assert time[0] == 0. and time[-1] == horizon and (np.diff(time) > 0).all()
                case_frames += len(time)
                if "states" in z:
                    states = z["states"].copy()
                    assert len(states) == len(time)
                    full_phases += 1
                else:
                    primary = z["primary"]
                    assert len(primary) == len(time) and np.isfinite(primary).all() and primary.min() >= -5e-10
                    assert z["snapshot_times"][0] == 0. and z["snapshot_times"][-1] == horizon
                    states = z["snapshot_states"].copy()
                    compact_phases += 1
                assert np.isfinite(states).all()
                fields = states[:, :-1] if spring else states
                local_min = float(fields.min())
                assert local_min >= -5e-10
                case_min = min(case_min, local_min)
                if spring:
                    assert (states[:, -1] > 0).all()
                    amount = states[:, -1] * (states[:, :n] + states[:, n:-1]).mean(axis=1)
                    np.testing.assert_array_equal(z["material_positions"], q)
                    centre_key = "cue_centre_fraction" if "cue_centre_fraction" in z else "cue_center_fraction"
                    assert float(z[centre_key]) == case["cue_centre_fraction"]
                    assert float(z["cue_peak"]) == (case["requested_peak"] if phase == "cue" else 0.)
                    if "full_physical_support" in z:
                        assert float(z["full_physical_support"]) == case["full_physical_support"]
                    if "resting_circumference" in z:
                        assert float(z["resting_circumference"]) == case["resting_circumference"]
                else:
                    np.testing.assert_array_equal(z["x"], grid.x)
                    np.testing.assert_array_equal(z["stimulus"], signal if phase == "cue" else np.zeros(n))
                    amount = np.stack([states[:, list(group), :].sum(axis=(1, 2)) * grid.dx for group in model.conserved_groups], axis=1)
                pool_error = float(np.max(np.abs(amount / amount0 - 1.)))
                assert pool_error < 1e-7
                case_error = max(case_error, pool_error)
                endpoints.append((states[0].copy(), states[-1].copy()))
        np.testing.assert_array_equal(endpoints[0][0], initial)
        np.testing.assert_array_equal(endpoints[0][1], endpoints[1][0])
        if case["execution"] == "new_native_branch":
            status = producer.read(run / "workers" / f'{case["worker_slot"]:03d}' / "status.json")
            assert status["status"] == "complete" and status["complete"] and status["condition_id"] == cid
            assert len(status["phases"]) == 2
            for item, phase in zip(status["phases"], ("cue", "recovery")):
                assert item["sha256"] == row["actual_" + phase + "_sha256"]
        frames += case_frames
        minimum, maximum_pool_error = min(minimum, case_min), max(maximum_pool_error, case_error)
        cases.append(dict(condition_id=cid, execution=case["execution"], frames=case_frames, native_minimum=case_min,
                          maximum_relative_pool_error=case_error, exact_basal_start=True, exact_withdrawal_restart=True,
                          actual_cue_sha256=row.actual_cue_sha256, actual_recovery_sha256=row.actual_recovery_sha256))
    result = dict(generated_at=producer.now(), status="passed", array_job_id=launch["array_job_id"],
        manifest_sha256=producer.sha(run / "manifest.json"), inventory_sha256=producer.sha(run / "results/condition_inventory.csv"),
        audit_script_sha256=producer.sha(__file__), requested=requested, complete=requested, exact_reused_pairs=exact_reuses, new_native_pairs=new_native,
        scheduler_completed_zero_exit=new_native, archive_hashes_verified=2*requested, frozen_sources_verified=len(manifest["source_sha256"]), frozen_inputs_verified=len(manifest["input_sha256"]),
        original_source_paths_verified=len(checked_source_paths), saved_primary_frames=frames, full_state_phases=full_phases,
        compact_primary_plus_full_snapshot_phases=compact_phases, minimum_native_field=minimum,
        maximum_relative_conserved_amount_error=maximum_pool_error, exact_basal_starts=requested, exact_cue_to_release_restarts=requested,
        overlap_conditions_retained=sum(c["initial_overlap"] for c in manifest["cases"]), excluded_for_overlap=0,
        scope="Every saved full state; compact histories checked for primary-field finiteness/nonnegativity and all available full snapshots. Full-field tail conclusions remain separate in derived tables.")
    pd.DataFrame(cases).to_csv(run / "results/native_terminal_inventory.csv", index=False)
    producer.write(run / "results/native_terminal_validation.json", result)
    (run / "results/scheduler_terminal.txt").write_text(scheduler)
    producer.write(run / "results/native_source_paths_verified.json", checked_source_paths)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    audit()
