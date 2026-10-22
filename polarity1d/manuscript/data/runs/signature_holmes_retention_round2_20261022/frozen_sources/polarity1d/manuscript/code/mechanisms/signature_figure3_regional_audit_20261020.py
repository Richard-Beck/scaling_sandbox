"""Independent reconstruction of corrected Figure 3 A2/B saved native fields."""
from pathlib import Path
import argparse
import json

import numpy as np

import signature_figure3_native_20261020 as native
import signature_figure3_campaign_20261020 as legacy
import signature_figure3_regional_campaign_20261020 as campaign
import signature_figure3_regions_20261020 as regions
import signature_figure3_geometry_20261020 as geometry
from signature_figure3_analysis_20261020 import write_csv

read, write, sha, now = native.read, native.write, native.sha, native.now


def audit(run, selected_slots=None):
    plan = legacy.verify(run)
    for module in (native, regions, geometry, campaign):
        relative = str(Path(module.__file__).relative_to(native.ROOT))
        assert sha(Path(module.__file__)) == plan["source_sha256"][relative], "Analysis/source definition mismatch: "+relative
    jobs = plan["jobs"] if selected_slots is None else [plan["jobs"][s] for s in selected_slots]
    rows, checks = [], []
    for job in jobs:
        placement = plan["placements"][job["placement_id"]]
        row = dict(slot=job["slot"], arm=job["arm"], model=placement["model"], size_parameter=placement["size_parameter"],
            control=job.get("control"), mirrored=job.get("mirrored", False),
            placement_id=placement["placement_id"], initial_edge_gap_10=placement["edge_gap_10"],
            initial_classification_10=placement["classification_10"], initial_classification_50=placement["classification_50"],
            basal_circumference=placement["basal_circumference"], initial_10_percent_width=placement["initial_10_percent_contour"]["width"],
            cue_support=placement["support"], initial_center_distance=placement["center_distance"],
            initial_short_route=placement["short_route"], initial_long_route=placement["long_route"])
        statuspath = run/"workers"/f'{job["slot"]:04d}'/"status.json"
        if not statuspath.exists():
            rows.append(dict(**row, execution_status="not_started", scientific_status="unavailable")); continue
        record = read(statuspath)
        if record["status"] != "complete":
            rows.append(dict(**row, execution_status=record["status"], scientific_status="unavailable", error=record.get("error"))); continue
        placement = record["placement"]
        adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
        with np.load(run/placement["input"]) as z:
            initial = z["state"].copy()
        previous, elapsed = initial.copy(), 0.
        initial_amounts = adapter.amounts(initial)
        all_times, all_states, phase_records = [], [], []
        tail = None
        for phase in record["phases"]:
            archive = run/phase["path"]
            assert sha(archive) == phase["sha256"]
            with np.load(archive) as z:
                times, states = z["time"].copy(), z["states"].copy()
            np.testing.assert_array_equal(states[0], previous)
            assert np.isfinite(states).all() and np.all(np.diff(times) > 0)
            assert abs(phase["absolute_start"]-elapsed) <= 1e-8
            amount_error = float(np.max(abs(np.asarray([adapter.amounts(y) for y in states])/initial_amounts-1)))
            minimum = min(adapter.concentration_minimum(y) for y in states)
            assert amount_error <= 1e-7 and minimum >= -5e-10
            assert phase["max_amount_relative_error"] <= 1e-7 and phase["minimum_concentration"] >= -5e-10
            use = slice(None) if not all_times else slice(1, None)
            all_times.extend((elapsed+times[use]).tolist()); all_states.extend(states[use])
            elapsed += float(times[-1]); previous = states[-1].copy()
            phase_records.append(dict(kind=phase["kind"], amount_error=amount_error, minimum=minimum, full_native_frames=len(times)))
            if phase["kind"].startswith("cue_free"):
                tail = (times, states)
        assert abs(elapsed-record["completed_native_horizon"]) <= 1e-7
        assert sha(run/record["readouts"]) == record["readouts_sha256"]
        yy, tt = np.asarray(all_states), np.asarray(all_times)
        p = yy[:, :adapter.cells]
        vector = p @ np.exp(2j*np.pi*adapter.q)/p.sum(axis=1)
        reconstructed = [campaign.spatial(adapter, y, placement, detailed=False) for y in yy]
        with np.load(run/record["readouts"]) as z:
            np.testing.assert_array_equal(z["time"], tt)
            np.testing.assert_allclose(z["vector_real"], vector.real, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(z["vector_imag"], vector.imag, rtol=1e-12, atol=1e-14)
            np.testing.assert_array_equal(z["state_class"], [r["state_class"] for r in reconstructed])
            for name in ("incumbent_excess", "challenger_excess"):
                values = np.asarray([r.get(name, np.nan) for r in reconstructed])
                np.testing.assert_allclose(z[name], values, rtol=1e-12, atol=1e-12, equal_nan=True)
        row.update(execution_status="complete", scientific_status="verified_native_fields", completed_horizon=elapsed,
            status_path=str(statuspath), status_sha256=sha(statuspath), final_vector_magnitude=float(abs(vector[-1])),
            final_vector_heading=float(np.angle(vector[-1])) if abs(vector[-1]) > 1e-10 else None,
            minimum_vector_magnitude=float(abs(vector).min()), heading_coverage_005=float(np.mean(abs(vector) >= .05)))
        if tail is not None:
            finite = native.finite_tail(adapter, *tail)
            assert finite["both_late_windows_steady"] == record["finite_tail"]["both_late_windows_steady"]
            final_fate = campaign.fate(adapter, previous, placement)
            assert final_fate["winner"] == record["final_winner"]["winner"]
            row.update(winner=final_fate["winner"], final_topology=final_fate["state_class"], finite_tail_pass=finite["both_late_windows_steady"])
        if job["arm"] == "A2":
            source = read(Path(job["source_refinement_status"]))
            assert sha(Path(job["source_refinement_status"])) == job["source_refinement_status_sha256"]
            np.testing.assert_allclose(job["cue_duration"], job["multiplier"]*source["handover_time"], rtol=1e-13)
            row.update(multiplier=job["multiplier"], cue_duration=job["cue_duration"], own_handover_time=source["handover_time"])
        else:
            row.update(target_ratio=job["target_ratio"], seed_status=record["seed_status"], B_claim_eligible=record["B_distance_claim_eligible"])
            if record.get("withdrawal_time") is not None:
                archive = run/record["exact_rhs_archive"]
                assert sha(archive) == record["exact_rhs_sha256"]
                with np.load(archive) as z:
                    withdrawal = z["state"].copy()
                    np.testing.assert_array_equal(z["cue_free_rhs"], adapter.rhs(withdrawal))
                    np.testing.assert_array_equal(z["cue_on_rhs"], adapter.rhs(withdrawal, placement["cue_peak"], placement["cue_center_fraction"]))
                spatial = campaign.spatial(adapter, withdrawal, placement)
                gap = campaign.measured_gap(adapter, withdrawal, placement, .1)
                assert gap["available"] == record["measured_withdrawal_gap"]["available"]
                if gap["available"]:
                    np.testing.assert_allclose(gap["measured_gap"], record["measured_withdrawal_gap"]["measured_gap"], atol=1e-10)
                if record["seed_event"]["kind"] == "continuous_two_region_ratio":
                    assert spatial["two_region_available"]
                    difference = abs(spatial["challenger_excess"]-job["target_ratio"]*spatial["incumbent_excess"])/record["basal_complete_excess"]
                    assert difference <= 1e-7
                budget = regions.frozen_basin_exchange_budget(adapter, withdrawal, spatial)
                assert budget["available"] == record["instantaneous_cue_free_exchange"]["available"]
                if budget["available"]:
                    for calculated, saved in zip(budget["regions"], record["instantaneous_cue_free_exchange"]["regions"]):
                        for name in ("exact_rate", "chemical_rate", "diffusion_rate", "concentration_dilution_rate", "physical_weight_dilation_rate"):
                            np.testing.assert_allclose(calculated[name], saved[name], rtol=1e-12, atol=1e-12)
                if record["B_distance_claim_eligible"]:
                    assert gap["available"] and record["seed_event"]["continuous_root"]
                    assert record["seed_status"] == "ratio_reached_with_two_regions"
                    resolution = campaign.resolution_summary(adapter, *tail, placement)
                    assert resolution["primary"] == record["B_resolution"]["primary"]
                    row.update(resolution_status=resolution["primary"]["status"], resolution_lower=resolution["primary"]["lower"],
                        resolution_upper=resolution["primary"]["upper"], one_region_persists_through_cap=resolution["one_region_persists_through_cap"])
                row.update(withdrawal_time=record["withdrawal_time"], actual_ratio=record["actual_withdrawal_ratio"],
                    withdrawal_circumference=adapter.circumference(withdrawal), measured_gap_10=gap.get("measured_gap"),
                    measured_center_distance=gap.get("center_distance"), measured_long_center_distance=gap.get("long_center_distance"),
                    gap_status=gap["reason"], peak_center_distance=spatial.get("peak_center_distance"))
        checks.append(dict(slot=job["slot"], phases=phase_records, event_RHS_readout_and_finite_tail_reconstructed=True))
        rows.append(row)
    scope = "all" if selected_slots is None else "selected"
    result = dict(generated_at=now(), requested=len(jobs), verified=len(checks), complete=len(checks) == len(jobs),
        source_manifest_sha256=sha(run/"manifest.json"), records=checks,
        scope="Independent saved-native-field hashes/continuity/conservation/positivity, complete-region readouts, native exact RHS, event equality, global vectors and finite-tail gates; no trajectory replay")
    write(run/f"results/{scope}_independent_native_audit.json", result)
    write_csv(run/f"results/{scope}_evidence.csv", rows)
    if selected_slots is None and plan["arm"] == "B":
        reference = {(r["placement_id"], r.get("target_ratio")):r for r in rows if not r.get("control")}
        mirrors = []
        for row in rows:
            if row.get("control") != "mirror_placement":
                continue
            normal = reference.get((row["placement_id"], row.get("target_ratio")))
            if normal is None or row["execution_status"] != normal["execution_status"] or row["execution_status"] != "complete":
                mirrors.append(dict(slot=row["slot"], status="unavailable")); continue
            assert row["seed_status"] == normal["seed_status"]
            assert row.get("winner") == normal.get("winner")
            assert row.get("resolution_status") == normal.get("resolution_status")
            mirrors.append(dict(slot=row["slot"], status="verified", seed_status=row["seed_status"],
                winner=row.get("winner"), resolution_status=row.get("resolution_status"),
                withdrawal_time_difference=abs(row.get("withdrawal_time", 0.)-normal.get("withdrawal_time", 0.))))
        write(run/"results/mirror_comparison.json", dict(generated_at=now(), mirrors=mirrors,
            requested=sum(j.get("control") == "mirror_placement" for j in plan["jobs"]),
            verified=sum(r["status"] == "verified" for r in mirrors)))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--slots", type=str)
    args = parser.parse_args()
    slots = [int(v) for v in args.slots.split(",")] if args.slots else None
    result = audit(args.run.resolve(), slots)
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))


if __name__ == "__main__":
    main()
