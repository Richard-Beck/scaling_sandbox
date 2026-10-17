"""Read-only numerical/data audit of the common ring active-loss assay."""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import vision_ring_recovery_20261017 as producer
import ring_distance_prototype_20261017 as ring


def audit(run):
    manifest = producer.verify(run)
    table = list(csv.DictReader((run/"results/recovery.csv").open()))
    assert len(table) == 2106
    lookup = {(int(r["id"]), float(r["epsilon"])):r for r in table}
    checks = []; worst_pool = 0.; worst_rate = 0.; curve_count = 0
    for case in manifest["cases"]:
        folder = run/"workers"/f'{case["worker_slot"]:02d}'/f'case_{case["id"]:04d}'
        status = producer.read(folder/"status.json")
        assert status["complete"] and status["status"] == "complete"
        assert producer.sha(folder/"curves.npz") == status["curves_sha256"]
        state, basal, grid, native, geometry, cue = producer.load_case(run, case)
        with np.load(folder/"curves.npz") as z:
            np.testing.assert_array_equal(z["source_state"], state)
            np.testing.assert_array_equal(z["basal_state"], basal)
            np.testing.assert_array_equal(z["stimulus"], cue)
            np.testing.assert_array_equal(z["core_weights"], geometry["weights"])
            np.testing.assert_array_equal(z["direction"], geometry["direction"])
            times = z["time"]; assert times[0] == 0 and times[-1] == case["tau"] and np.all(np.diff(times) > 0)
            for precision in ("tight", "tighter"):
                for label, epsilon in [("sham", 0.)]+[(f"eps{e:g}", e) for e in producer.EPSILONS]:
                    name = precision+"_"+label
                    initial = z[name+"_initial_state"]; final = z[name+"_final_state"]
                    removed = epsilon*geometry["direction"]
                    np.testing.assert_array_equal(initial[0], state[0]-removed)
                    np.testing.assert_array_equal(initial[1], state[1]+removed)
                    np.testing.assert_array_equal(initial[2:], state[2:])
                    assert initial.min() >= 0 and final.min() >= -1e-7
                    np.testing.assert_allclose(initial[0]+initial[1], state[0]+state[1], atol=3e-14, rtol=0)
                    for index, full in ((0, initial), (-1, final)):
                        core = float(full[0]@geometry["weights"]); whole = float(full[0].sum()*grid.dx)
                        np.testing.assert_allclose(z[name+"_core_amount"][index], core, atol=1e-12, rtol=1e-13)
                        np.testing.assert_allclose(z[name+"_whole_amount"][index], whole, atol=1e-12, rtol=1e-13)
                        reaction = native.reaction(0., grid.x, full, cue)[0]
                        diffusion = native.diffusivities[0]*ring.periodic_laplacian(full[0], grid.dx)
                        rate = float((reaction+diffusion)@geometry["weights"])
                        error = abs(z[name+"_core_rate"][index]-rate)
                        worst_rate = max(worst_rate, float(error))
                        np.testing.assert_allclose(z[name+"_core_rate"][index], rate, atol=1e-11, rtol=1e-11)
                        np.testing.assert_allclose(z[name+"_pools"][index], ring.pool_amounts(native, grid, full), atol=1e-11, rtol=1e-13)
                    drift = float(np.max(abs(z[name+"_pools"]/z[name+"_pools"][0]-1)))
                    assert drift <= 1e-8; worst_pool = max(worst_pool, drift)
                    if label in ("sham", "eps0.005"):
                        np.testing.assert_array_equal(z[name+"_snapshot_states"][0], initial)
                        np.testing.assert_array_equal(z[name+"_snapshot_states"][-1], final)
                    curve_count += 1
                for key in ("initial_state", "final_state", "core_amount", "whole_amount", "core_rate", "whole_rate", "pools"):
                    np.testing.assert_array_equal(z[precision+"_sham_"+key], z[precision+"_eps0_"+key])
            stats = {}
            for epsilon in (.001, .005, .01):
                label = f"eps{epsilon:g}"
                sham = z["tight_sham_core_amount"]; probe = z["tight_"+label+"_core_amount"]
                deficit = sham-probe
                ref_deficit = z["tighter_sham_core_amount"]-z["tighter_"+label+"_core_amount"]
                F = (deficit[0]-deficit)/deficit[0]
                ref_F = (ref_deficit[0]-ref_deficit)/ref_deficit[0]
                expected_D0 = epsilon*float(geometry["direction"]@geometry["weights"])
                np.testing.assert_allclose(deficit[0], expected_D0, atol=1e-12, rtol=1e-10)
                zero = sham-z["tight_eps0_core_amount"]
                E = np.maximum.reduce([2*abs(F-ref_F), 2*abs(zero)/deficit[0],
                    64*np.finfo(float).eps*(abs(sham)+abs(probe)+abs(sham[0])+abs(probe[0]))/deficit[0]])
                row = lookup[(case["id"],epsilon)]
                for field,value in dict(F=F[-1], K=F[-1]/case["tau"], D0=deficit[0], Dtau=deficit[-1], E=E[-1],
                    sham_core_change=sham[-1]-sham[0], probe_core_change=probe[-1]-probe[0]).items():
                    np.testing.assert_allclose(float(row[field]), value, atol=3e-12, rtol=1e-10)
                stats[epsilon] = dict(F=F, E=E)
            mask = times > 0
            Emax = max(float(v["E"][mask].max()) for v in stats.values())
            saturation = max(float((abs(v["F"])+v["E"])[mask].max()) for v in stats.values())
            excess = -np.inf
            for ea,eb in ((.001,.005),(.001,.01),(.005,.01)):
                a,b = stats[ea],stats[eb]
                excess = max(excess,float((abs(a["F"]-b["F"])-.05*np.maximum(abs(a["F"]),abs(b["F"]))-3*(a["E"]+b["E"]))[mask].max()))
            passed = Emax <= 1e-5 and saturation <= .2 and excess <= 0
            assert all((lookup[(case["id"],e)]["passes"] == "True") == passed for e in (.001,.005,.01))
            checks.append(dict(id=case["id"], model=case["model"], stage=case["stage"],
                frozen_inputs_and_direction_exact=True, initial_kick_and_untouched_fields_exact=True,
                endpoint_native_amounts_rates_pools_pass=True, sham_zero_identity_pass=True,
                recomputed_prefix_Emax=Emax, recomputed_prefix_saturation=saturation,
                recomputed_linearity_excess=excess, passes=bool(passed),
                curves_path=str(folder/"curves.npz"), curves_sha256=status["curves_sha256"]))
    result = dict(generated_at=producer.now(), passed=all(c["passes"] for c in checks),
        cases=len(checks), complete_trajectory_trials=curve_count, amplitude_summary_rows=len(table),
        max_pool_relative_drift=worst_pool, max_endpoint_rate_allocation_discrepancy=worst_rate,
        full_grid_gate_failures=[c for c in checks if not c["passes"]], cases_checked=checks,
        auditor_sha256=producer.sha(Path(__file__)), producer_frozen_source_hashes=manifest["source_sha256"],
        interpretation="Read-only equation/data audit with native endpoint RHS recomputation; not an independent re-integration or certified error bound.")
    producer.write(run/"results/independent_data_audit.json", result)
    return {k:v for k,v in result.items() if k not in ("cases_checked", "producer_frozen_source_hashes")}


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--run", type=Path, default=producer.DEFAULT)
    args=p.parse_args(); print(json.dumps(audit(args.run.resolve()), indent=2))
