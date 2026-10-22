"""Independent saved-native-field audit and complete Figure 3 evidence tables.

Runs are censused from their frozen job plan, never discovered by successful
output globbing. No missing, failed or unqualified case becomes a lower bound.
Scientific plots/translation claims remain gated on complete control evidence.
"""
from pathlib import Path
import argparse
import csv
import json

import numpy as np
import signature_figure3_native_20261020 as native
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_campaign_20261020 as campaign

read, write, sha, now = native.read, native.write, native.sha, native.now


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = sorted(set().union(*(r.keys() for r in rows))) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def phase_fields(run, phase):
    path = run/phase["path"]
    assert sha(path) == phase["sha256"], path
    with np.load(path) as z:
        return z["time"].copy(), z["states"].copy()


def audit(run, stage):
    plan = campaign.verify(run)
    stagepath = run/f"plans/{stage}.json"
    assert sha(stagepath) == read(run/f"plans/{stage}.sha256.json")["sha256"]
    jobs = read(stagepath)["jobs"]
    summaries, checks = [], []
    for job in jobs:
        placement = plan["placements"][job["placement_id"]]
        base = dict(stage=stage, slot=job["slot"], arm=job["arm"],
            model=placement["model"], size_parameter=placement["size_parameter"],
            basal_circumference=placement["basal_circumference"], cells=placement["cells"],
            placement_id=job["placement_id"], gap=placement["edge_gap_10"],
            gap_50=placement["edge_gap_50"], classification=placement["classification_10"],
            classification_50=placement["classification_50"],
            center_distance=placement["center_distance"], target_angle_radians=placement["target_angle_radians"],
            control=job.get("control", ""), target_ratio=job.get("target_ratio"),
            pulse_multiplier=job.get("multiplier"), cue_duration=job.get("cue_duration"))
        path = run/"workers"/stage/f'{job["slot"]:04d}'/"status.json"
        if not path.exists():
            summaries.append(dict(**base, execution_status="not_started", scientific_availability="unavailable")); continue
        record = read(path)
        if not record.get("complete"):
            summaries.append(dict(**base, execution_status=record["status"], scientific_availability="running")); continue
        if record["status"] != "complete":
            summaries.append(dict(**base, execution_status=record["status"], scientific_availability="unavailable_native_error")); continue
        adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
        with np.load(run/placement["input"]) as z:
            original = z["state"].copy()
        previous, total, amount_error, minimum = original.copy(), 0., 0., np.inf
        reference_amounts = adapter.amounts(original)
        tt, yy, phase_check = [], [], []
        tail_pair = None
        for phase in record["phases"]:
            t, y = phase_fields(run, phase)
            np.testing.assert_array_equal(y[0], previous)
            assert np.isfinite(y).all() and (len(t) == 1 or np.all(np.diff(t) > 0))
            assert y.shape[1:] == original.shape
            assert phase["status"] in ("complete", "event", "event_at_start")
            assert abs(phase["absolute_start"]-total) < 1e-8
            amounts = np.asarray([adapter.amounts(state) for state in y])
            amount_error = max(amount_error, float(abs(amounts/reference_amounts-1).max()))
            minimum = min(minimum, min(adapter.concentration_minimum(state) for state in y))
            assert amount_error <= 1e-7 and minimum >= -5e-10
            assert phase["max_amount_relative_error"] <= 1e-7 and phase["minimum_concentration"] >= -5e-10
            use = slice(None) if not tt else slice(1, None)
            tt.extend((total+t[use]).tolist()); yy.extend(y[use])
            total += float(t[-1])
            previous = y[-1].copy()
            phase_check.append(dict(kind=phase["kind"], full_native_fields=len(t),
                                    hash_verified=True, exact_continuity=True))
            if phase["kind"] in ("cue_free_tail", "cue_free_resolution_tail"):
                tail_pair = (t, y)
        np.testing.assert_allclose(total, record["completed_native_horizon"], atol=1e-7)
        t, y = np.asarray(tt), np.asarray(yy)
        assert sha(run/record["readouts"]) == record["readouts_sha256"]
        with np.load(run/record["readouts"]) as z:
            np.testing.assert_array_equal(z["time"], t)
            p = y[:, :adapter.cells]
            complex_vector = p @ np.exp(2j*np.pi*adapter.q)/p.sum(axis=1)
            np.testing.assert_allclose(z["vector_real"], complex_vector.real, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(z["vector_imag"], complex_vector.imag, rtol=1e-12, atol=1e-14)
        summary = dict(**base, execution_status="complete", scientific_availability="observed_native_run",
            completed_horizon=total, final_vector_magnitude=float(abs(complex_vector[-1])),
            final_vector_heading=float(np.angle(complex_vector[-1])) if abs(complex_vector[-1]) > 1e-10 else None,
            vector_coverage_005=float(np.mean(abs(complex_vector) >= .05)),
            vector_coverage_010=float(np.mean(abs(complex_vector) >= .1)),
            vector_coverage_020=float(np.mean(abs(complex_vector) >= .2)),
            minimum_vector_magnitude=float(abs(complex_vector).min()),
            minimum_concentration=minimum, maximum_amount_relative_error=amount_error,
            status_path=str(path), status_sha256=sha(path))
        if job["arm"] == "A1":
            event = record.get("handover_time")
            summary.update(handover_time=event, handover_lower_bound=record.get("handover_lower_bound"),
                handover_status=record["handover_status"], winner=record["final_winner"]["winner"])
            if event is not None:
                indices = np.flatnonzero(abs(t-event) < 1e-7)
                assert len(indices)
                value = campaign.scalar(adapter, y[int(indices[0])], adapter.primary(original), record["placement"])
                if event > 0:
                    assert abs(value["challenger_fraction"]-value["incumbent_fraction"]) <= 1e-6
                summary["event_root_excess_difference"] = value["challenger_fraction"]-value["incumbent_fraction"]
        elif job["arm"] == "A2":
            summary.update(winner=record["final_winner"]["winner"],
                           finite_tail_pass=record["finite_tail"]["both_late_windows_steady"])
        else:
            summary.update(seed_status=record["seed_status"], B_claim_eligible=record.get("B_distance_claim_eligible", False),
                withdrawal_time=record.get("withdrawal_time"), actual_ratio=record.get("actual_withdrawal_ratio"))
            measured = record["measured_withdrawal_gap"]
            summary["withdrawal_gap_available"] = measured["available"]
            summary["measured_withdrawal_gap"] = measured.get("measured_gap")
            if record.get("withdrawal_time") is not None:
                rhs_path = run/record["exact_rhs_archive"]
                assert sha(rhs_path) == record["exact_rhs_sha256"]
                with np.load(rhs_path) as z:
                    withdrawal = z["state"].copy()
                    np.testing.assert_array_equal(z["cue_free_rhs"], adapter.rhs(withdrawal))
                    np.testing.assert_array_equal(z["cue_on_rhs"], adapter.rhs(withdrawal, placement["cue_peak"], placement["cue_center_fraction"]))
                metric = campaign.current_regions(adapter, withdrawal, adapter.primary(original), placement)
                assert metric["available"] == measured["available"]
                if metric["available"]:
                    np.testing.assert_allclose(metric["measured_gap"], measured["measured_gap"], atol=1e-10)
                if summary["B_claim_eligible"]:
                    assert metric["available"] and np.isfinite(metric["measured_gap"])
                    assert record["signed_excess_ratio_gate"]["available"]
                    assert record["seed_status"] == "ratio_reached_with_two_regions"
                    resolution = record["B_resolution"]["resolution"]
                    summary.update(resolution_lower=resolution["lower"], resolution_upper=resolution["upper"],
                                   resolution_status=resolution["status"], winner=record["B_resolution"]["final_winner"]["winner"])
            if "finite_tail" in record:
                summary["finite_tail_pass"] = record["finite_tail"]["both_late_windows_steady"]
        if tail_pair is not None:
            independent_tail = native.finite_tail(adapter, *tail_pair)
            assert independent_tail["both_late_windows_steady"] == record["finite_tail"]["both_late_windows_steady"]
        checks.append(dict(slot=job["slot"], phases=phase_check, max_sampled_amount_error=amount_error,
            exact_event_and_RHS_checked=True, no_trajectory_reintegration=True))
        summaries.append(summary)
    complete = len(checks) == len(jobs)
    result = dict(generated_at=now(), stage=stage, requested=len(jobs), independently_checked=len(checks),
        complete=complete, passed=complete, records=checks,
        source_manifest_sha256=sha(run/"manifest.json"),
        audit_scope="Exact saved native-field hashes/phase continuity/amount/positivity, event-root equality, same-run vector, withdrawal RHS and finite-tail reconstruction; not an independent trajectory replay")
    write(run/f"results/{stage}_independent_native_audit.json", result)
    write_csv(run/f"results/{stage}_evidence.csv", summaries)
    return result, summaries


def A2_brackets(run):
    plan = campaign.verify(run)
    stage = read(run/"plans/A2.json")
    grouped = {}
    for job in stage["jobs"]:
        grouped.setdefault(job["placement_id"], []).append(job)
    rows = []
    for placement_id, jobs in grouped.items():
        placement = plan["placements"][placement_id]
        tested = []
        for job in sorted(jobs, key=lambda j: j["cue_duration"]):
            statuspath = run/"workers/A2"/f'{job["slot"]:04d}'/"status.json"
            if not statuspath.exists() or read(statuspath).get("status") != "complete":
                tested.append(dict(duration=job["cue_duration"], status="unavailable")); continue
            record = read(statuspath)
            tested.append(dict(duration=job["cue_duration"], status=record["final_winner"]["winner"],
                               finite_tail=record["finite_tail"]["both_late_windows_steady"]))
        wins = [i for i, value in enumerate(tested) if value["status"] == "challenger"]
        first = wins[0] if wins else None
        missing = any(v["status"] == "unavailable" for v in tested)
        row = dict(placement_id=placement_id, model=placement["model"], size_parameter=placement["size_parameter"],
            gap=placement["edge_gap_10"], classification=placement["classification_10"],
            all_tested_durations_and_fates=json.dumps(tested), native_trials_complete=not missing,
            threshold_scope="Tested-grid shortest observed switch, not a universal or monotonic commitment threshold")
        if missing:
            row.update(status="unavailable_incomplete_trials")
        elif first is None:
            row.update(status="no_observed_switch_at_tested_durations", lower=max(v["duration"] for v in tested), upper=None)
        else:
            previous_failures = [v["duration"] for v in tested[:first] if v["status"] != "challenger"]
            row.update(status="bracketed_shortest_tested_switch", lower=max(previous_failures) if previous_failures else 0.,
                upper=tested[first]["duration"],
                later_non_switch_after_first_switch=any(v["status"] != "challenger" for v in tested[first+1:]))
        rows.append(row)
    for unavailable in stage["unavailable"]:
        rows.append(dict(placement_id=unavailable["placement_id"], status=unavailable["status"],
                         native_trials_complete=False, unavailable_reason="Own A1 handover did not provide a positive observed duration"))
    write_csv(run/"results/A2_switching_brackets.csv", rows)
    return rows


def mirror_comparison(run, evidence):
    primary = {row["placement_id"]: row for row in evidence if row["arm"] == "A1" and not row.get("control")}
    rows = []
    for row in evidence:
        if row.get("control") != "mirror_placement":
            continue
        reference = primary[row["placement_id"]]
        complete = row["execution_status"] == reference["execution_status"] == "complete"
        result = dict(model=row["model"], size_parameter=row["size_parameter"], gap=row["gap"],
            mirror_slot=row["slot"], reference_slot=reference["slot"], complete=complete)
        if complete:
            t1, t2 = reference.get("handover_time"), row.get("handover_time")
            result.update(same_winner=reference["winner"] == row["winner"],
                reference_handover=t1, mirror_handover=t2,
                relative_handover_difference=abs(t1-t2)/max(t1, t2, 1.) if t1 is not None and t2 is not None else None,
                same_censoring=(t1 is None) == (t2 is None))
        rows.append(result)
    write_csv(run/"results/mirror_and_antipode_comparison.csv", rows)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("audit", "A2-brackets"))
    parser.add_argument("--run", type=Path, default=campaign.DEFAULT)
    parser.add_argument("--stage", default="primary")
    args = parser.parse_args()
    run = args.run.resolve()
    if args.action == "audit":
        result, evidence = audit(run, args.stage)
        if args.stage == "primary":
            mirror_comparison(run, evidence)
        print(json.dumps({k: v for k, v in result.items() if k != "records"}, indent=2))
    else:
        print(json.dumps(dict(rows=len(A2_brackets(run))), indent=2))


if __name__ == "__main__":
    main()
