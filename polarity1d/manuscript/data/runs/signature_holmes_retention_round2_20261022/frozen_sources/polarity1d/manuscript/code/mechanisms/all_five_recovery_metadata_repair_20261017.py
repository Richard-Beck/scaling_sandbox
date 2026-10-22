"""Recover Spring completion metadata from complete immutable pilot archives.

The original run's Spring workers saved all scientific arrays, then failed while
JSON-encoding a NumPy int32 solver counter. This performs no integration and
preserves original status/error records and frozen sources unchanged.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import sys
import numpy as np


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("run", type=Path)
    run = parser.parse_args().run.resolve()
    frozen = run/"frozen_sources/code/mechanisms/all_five_incumbent_recovery_20261017.py"
    specification = importlib.util.spec_from_file_location("frozen_incumbent_pilot", frozen)
    pilot = importlib.util.module_from_spec(specification); specification.loader.exec_module(pilot)
    manifest = pilot.verify(run)
    recovered = []
    for case in manifest["cases"]:
        if case["model"] != "spring": continue
        folder = run/"workers"/f'{case["slot"]:02d}'
        original = pilot.read(folder/"status.json")
        assert not original["complete"]
        errors = sorted((run/"logs").glob(f'*_{case["slot"]}.err'))
        assert len(errors) == 1
        error = errors[0].read_text()
        assert "Object of type int32 is not JSON serializable" in error
        assert 'write(folder/"status.json", result)' in error
        archive_path = folder/"curves.npz"; before = pilot.sha(archive_path)
        state, fields, x = pilot.load_state(run/case["input"], "spring")
        geometry = pilot.initial_geometry(case, state, fields, x)
        pools, differences = {}, {}
        maximum_initial_roundoff = 0.
        with np.load(archive_path) as z:
            trials = {}
            for prefix in ("tight", "tighter"):
                for tag in ("sham", "eps0", "eps0.001", "eps0.005", "eps0.01"):
                    name = prefix+"_"+tag
                    time = z[name+"_time"]
                    np.testing.assert_array_equal(time, manifest["continuous_samples"])
                    assert time[-1] == 1.
                    names = ("active", "inactive", "length", "left_endpoint", "right_endpoint")
                    f = {field: z[name+"_field_"+field].copy() for field in names}
                    assert all(np.isfinite(value).all() for value in f.values())
                    values = np.empty((len(time), len(state)))
                    values[:, :-2:2] = f["active"]; values[:, 1:-2:2] = f["inactive"]
                    values[:, -2] = f["left_endpoint"]; values[:, -1] = f["right_endpoint"]
                    epsilon = 0. if tag == "sham" else float(tag[3:])
                    expected = pilot.kick(state, "spring", geometry["direction"], epsilon)
                    initial_roundoff = float(np.max(abs(values[0]-expected)))
                    maximum_initial_roundoff = max(maximum_initial_roundoff, initial_roundoff)
                    # LSODA's polynomial output at t=0 can differ from its exact
                    # supplied restart input by one machine ulp.
                    np.testing.assert_allclose(values[0], expected, rtol=0, atol=2e-14)
                    trial = dict(time=time, fields=f, states=values)
                    core, whole = pilot.amounts(trial, case, geometry)
                    cr, wr = pilot.amount_rates(trial, case, geometry)
                    np.testing.assert_allclose(core, z[name+"_core_amount"], rtol=1e-13, atol=1e-14)
                    np.testing.assert_allclose(whole, z[name+"_whole_amount"], rtol=1e-13, atol=1e-14)
                    np.testing.assert_allclose(cr, z[name+"_core_rate"], rtol=1e-11, atol=1e-13)
                    np.testing.assert_allclose(wr, z[name+"_whole_rate"], rtol=1e-11, atol=1e-13)
                    pools[name] = pilot.pools(trial, case)
                    trials[name] = trial
            for prefix in ("tight", "tighter"):
                for field in trials[prefix+"_sham"]["fields"]:
                    np.testing.assert_array_equal(trials[prefix+"_sham"]["fields"][field], trials[prefix+"_eps0"]["fields"][field])
            for epsilon in pilot.EPSILONS:
                a, b = trials[f"tight_eps{epsilon:g}"], trials[f"tighter_eps{epsilon:g}"]
                differences[f"eps{epsilon:g}"] = {name: float(np.max(abs(a["fields"][name]-b["fields"][name]))) for name in a["fields"]}
        assert pilot.sha(archive_path) == before
        result = dict(generated_at=pilot.now(), producer="recovered_after_metadata_serialization_failure", complete=True,
            status="complete", case=case, curves_sha256=before, original_status_sha256=pilot.sha(folder/"status.json"),
            original_error_log_sha256=pilot.sha(errors[0]), original_error_log=str(errors[0].relative_to(run)),
            recovery_script_sha256=pilot.sha(Path(__file__)), new_integrations=0,
            geometry={k: v for k, v in geometry.items() if np.isscalar(v)},
            zero_kick_full_field_identity=True, tracked_pool_checks=pools,
            primary_reference_full_field_maxabs=differences,
            numerical={"completion_record": dict(kind="Full scientific archive independently reconstructed and verified; original solver counters unavailable")},
            independent_completion_evidence=dict(trajectories=10, exact_retained_restart_inputs_used=True,
                maximum_LSODA_initial_output_roundoff=maximum_initial_roundoff,
                all_prespecified_sample_times_present=True, all_trajectories_end_at_native_time=1.,
                moving_domain_physical_core_and_whole_amounts_recomputed=True,
                native_core_and_whole_rhs_rates_recomputed=True,
                original_status_and_failed_error_log_preserved=True,
                raw_archive_hash_unchanged=True))
        pilot.write(folder/"recovered_status.json", result)
        recovered.append(dict(slot=case["slot"], size=case["size_native"], archive_sha256=before,
                              maximum_pool_drift=max(v["max_relative_error"] for group in pools.values() for v in group)))
    pilot.write(run/"results/completion_metadata_recovery.json", dict(generated_at=pilot.now(),
        new_integrations=0, original_frozen_sources_unchanged=True,
        recovery_script_sha256=pilot.sha(Path(__file__)), cases=recovered))
    print(json.dumps(recovered, indent=2))


if __name__ == "__main__": main()
