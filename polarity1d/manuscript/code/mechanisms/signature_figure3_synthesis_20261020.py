"""Audited Figure 3 evidence tables and support-limited geometry remaps."""
from pathlib import Path
import argparse
import json

import numpy as np

import signature_figure3_native_20261020 as native
import signature_figure3_campaign_20261020 as legacy
import signature_figure3_geometry_20261020 as geometry
from signature_figure3_analysis_20261020 import write_csv

ROOT = native.ROOT
PRIMARY = ROOT/"data/runs/signature_figure3_campaign_20261020"
REGIONS = ROOT/"data/runs/signature_figure3_region_proposal_v2_20261020"
REFINED = ROOT/"data/runs/signature_figure3_handover_refinement_v2_20261020"
A2 = ROOT/"data/runs/signature_figure3_regional_A2_20261020"
B = ROOT/"data/runs/signature_figure3_regional_B_20261020"
DESTINATION = ROOT/"data/runs/signature_figure3_synthesis_20261020"
read, write, sha, now = native.read, native.write, native.sha, native.now


def A1_evidence(primary=PRIMARY, regional=REGIONS, refined=REFINED):
    plan = legacy.verify(primary)
    assert read(primary/"results/primary_independent_native_audit.json")["complete"]
    refinement = legacy.verify(refined)
    audit = read(refined/"results/independent_refinement_audit.json")
    assert audit["complete"] and audit["verified"] == len(refinement["jobs"])
    by_original = {refinement["jobs"][r["slot"]]["original_slot"]:r for r in audit["records"]}
    rows = []
    for job in plan["stages"]["primary"]:
        statuspath = primary/"workers/primary"/f'{job["slot"]:04d}'/"status.json"
        status = read(statuspath)
        placement = status["placement"]
        proposalpath = regional/"results"/f'primary_{job["slot"]:04d}_placement_{placement["placement_id"]:03d}_region_proposal.json'
        proposal = read(proposalpath)
        row = dict(slot=job["slot"], placement_id=placement["placement_id"], control=job.get("control"),
            model=placement["model"], size_parameter=placement["size_parameter"], basal_circumference=placement["basal_circumference"],
            gap=placement["edge_gap_10"], gap_50=placement["edge_gap_50"],
            classification=placement["classification_10"], classification_50=placement["classification_50"],
            initial_10_percent_width=placement["initial_10_percent_contour"]["width"], cue_support=placement["support"],
            center_distance=placement["center_distance"], target_angle_radians=placement["target_angle_radians"],
            short_route=placement["short_route"], long_route=placement["long_route"],
            basal_complete_excess=proposal["initial_complete_basal_excess"], final_topology=proposal["final_state_class"],
            original_native_horizon=status["completed_cue_horizon"],
            handover_time=None, handover_lower_bound=proposal["complete_region_handover_lower_bound"],
            handover_status=proposal["complete_region_handover_status"],
            native_status=str(statuspath), native_status_sha256=sha(statuspath),
            region_proposal=str(proposalpath), region_proposal_sha256=sha(proposalpath))
        if job["slot"] in by_original:
            record = by_original[job["slot"]]
            assert sha(Path(record["status_path"])) == record["status_sha256"]
            row.update(handover_time=record["handover_time"], handover_status="verified_continuous_complete_region_root",
                refinement_status=record["status_path"], refinement_status_sha256=record["status_sha256"])
        rows.append(row)
    assert len(rows) == 176
    return rows


def curve_value(rows, gap):
    """Piecewise-linear interpolation only between adjacent observed events.

    Exact censored placements retain their lower bound. No extrapolation or
    interpolation through a censored/topology-unavailable point is performed.
    """
    data = sorted((r for r in rows if r["classification"] == "external" and not r.get("control")), key=lambda r:r["gap"])
    if not data:
        return dict(available=False, status="unavailable_empty_curve")
    x = np.asarray([0. if abs(r["gap"]) < 1e-10 else r["gap"] for r in data])
    if gap < x[0]-1e-9 or gap > x[-1]+1e-9:
        return dict(available=False, status="outside_observed_external_gap_support", support=[float(x[0]), float(x[-1])])
    exact = np.flatnonzero(abs(x-gap) < 1e-9)
    if len(exact):
        r = data[int(exact[0])]
        if r.get("handover_time") is not None:
            return dict(available=True, status="observed_native_root", value=r["handover_time"], source_slots=[r["slot"]])
        return dict(available=False, status=r["handover_status"], lower_bound=r.get("handover_lower_bound"), source_slots=[r["slot"]])
    upper = int(np.searchsorted(x, gap)); lower = upper-1
    pair = [data[lower], data[upper]]
    if any(r.get("handover_time") is None for r in pair):
        return dict(available=False, status="interpolation_blocked_by_unavailable_or_censored_neighbor",
            neighbor_statuses=[r["handover_status"] for r in pair], source_slots=[r["slot"] for r in pair])
    fraction = float((gap-x[lower])/(x[upper]-x[lower]))
    return dict(available=True, status="piecewise_linear_between_adjacent_observed_roots",
        value=(1-fraction)*pair[0]["handover_time"]+fraction*pair[1]["handover_time"],
        source_slots=[r["slot"] for r in pair], interpolation_fraction=fraction)


def fixed_angle_remaps(rows):
    output = []
    for model in native.MODELS:
        by_size = {size:[r for r in rows if r["model"] == model and r["size_parameter"] == size and not r.get("control")]
                   for size in native.SIZES[model]}
        small = by_size[native.SIZES[model][0]]
        C0, w0, s0 = (small[0][k] for k in ("basal_circumference", "initial_10_percent_width", "cue_support"))
        for size, actual in by_size.items():
            C, width, support = (actual[0][k] for k in ("basal_circumference", "initial_10_percent_width", "cue_support"))
            for theta_degrees in (30., 45., 60., 75., 90., 105., 120., 135., 150., 165., 180.):
                theta = np.deg2rad(theta_degrees)
                g0, gap = theta*C0/(2*np.pi)-(w0+s0)/2, theta*C/(2*np.pi)-(width+support)/2
                baseline, remapped, measured = curve_value(small, g0), curve_value(small, gap), curve_value(actual, gap)
                row = dict(model=model, size_parameter=size, theta_degrees=theta_degrees, basal_circumference=C,
                    initial_10_percent_width=width, cue_support=support, baseline_gap=g0, measured_size_gap=gap,
                    small_same_angle_status=baseline["status"], geometry_only_status=remapped["status"],
                    measured_status=measured["status"], small_same_angle_time=baseline.get("value"),
                    geometry_only_time=remapped.get("value"), measured_time=measured.get("value"),
                    geometry_time_shift=None, matched_gap_size_time_shift=None,
                    all_components_available=baseline["available"] and remapped["available"] and measured["available"],
                    source_queries=json.dumps(dict(baseline=baseline, geometry_only=remapped, measured=measured)))
                if row["all_components_available"]:
                    row.update(geometry_time_shift=remapped["value"]-baseline["value"],
                        matched_gap_size_time_shift=measured["value"]-remapped["value"])
                    np.testing.assert_allclose(row["geometry_time_shift"]+row["matched_gap_size_time_shift"], measured["value"]-baseline["value"])
                output.append(row)
    return output


def A2_brackets(run=A2):
    plan = legacy.verify(run)
    audit = read(run/"results/all_independent_native_audit.json")
    assert audit["complete"] and audit["verified"] == len(plan["jobs"])
    grouped = {}
    for job in plan["jobs"]:
        grouped.setdefault(job["placement_id"], []).append(job)
    rows = []
    for placement_id, jobs in grouped.items():
        placement = plan["placements"][placement_id]
        tested = []
        for job in sorted(jobs, key=lambda j:j["cue_duration"]):
            statuspath = run/"workers"/f'{job["slot"]:04d}'/"status.json"
            status = read(statuspath)
            assert status["status"] == "complete"
            tested.append(dict(slot=job["slot"], multiplier=job["multiplier"], duration=job["cue_duration"],
                winner=status["final_winner"]["winner"], topology=status["final_winner"]["state_class"],
                finite_tail_pass=status["finite_tail"]["both_late_windows_steady"]))
        wins = [i for i,r in enumerate(tested) if r["winner"] == "challenger"]
        first = wins[0] if wins else None
        unresolved_shorter = first is not None and any(r["winner"] == "unresolved" for r in tested[:first])
        lower = max((r["duration"] for r in (tested if first is None else tested[:first]) if r["winner"] == "incumbent"), default=0.)
        rows.append(dict(placement_id=placement_id, model=placement["model"], size_parameter=placement["size_parameter"],
            gap=placement["edge_gap_10"], classification=placement["classification_10"],
            handover_time=jobs[0]["handover_time"], lower_confirmed_incumbent_duration=lower,
            shortest_observed_switch=None if first is None else tested[first]["duration"],
            status="no_observed_switch_at_tested_durations" if first is None else "shortest_tested_switch_with_unresolved_shorter_trials" if unresolved_shorter else "bracketed_shortest_tested_switch",
            later_non_switch_after_first_switch=False if first is None else any(r["winner"] != "challenger" for r in tested[first+1:]),
            tested_trials=json.dumps(tested), native_trials=5,
            threshold_scope="Finite tested grid; neither monotonic commitment nor full native-field stability is assumed"))
    rows.extend(dict(**r, native_trials=0) for r in plan["unavailable"])
    assert len(rows) == 143
    return rows


def control_handovers(A1):
    run = ROOT/"data/runs/signature_figure3_control_handover_refinement_20261020"
    source = ROOT/"data/runs/signature_figure3_A1_controls_20261020"
    plan, originals = legacy.verify(run), legacy.verify(source)
    audit = read(run/"results/independent_refinement_audit.json")
    assert audit["complete"]
    by_original = {plan["jobs"][r["slot"]]["original_slot"]:r for r in audit["records"]}
    output = []
    for job in originals["stages"]["controls"]:
        placement = originals["placements"][job["placement_id"]]
        reference = next(r for r in A1 if not r.get("control") and r["model"] == placement["model"] and
            r["size_parameter"] == 37.8 and abs(r["gap"]-placement["matched_gap"]) < 1e-9)
        record = by_original.get(job["slot"])
        current = None if record is None else record["handover_time"]
        baseline = reference["handover_time"]
        output.append(dict(model=placement["model"], size_parameter=37.8, gap=placement["matched_gap"],
            variant=job["control"], cells=placement["cells"], baseline_handover=baseline,
            control_handover=current, same_event_availability=(current is None) == (baseline is None),
            relative_time_difference=None if current is None or baseline is None else abs(current-baseline)/baseline))
    return output


def threshold_handovers(A1):
    output = []
    for row in A1:
        proposal = read(Path(row["region_proposal"]))
        trace = Path(proposal["traces"])
        assert sha(trace) == proposal["traces_sha256"]
        with np.load(trace) as z:
            times, two = z["time"].copy(), z["two_region_available"].astype(bool)
            for fraction, key in ((.25, "dominance_025"), (.5, "dominance_050"), (.75, "dominance_075")):
                value = z[key]
                eligible = two[:-1] & two[1:] & np.isfinite(value[:-1]) & np.isfinite(value[1:])
                found = np.flatnonzero(eligible & (value[:-1] < 0) & (value[1:] >= 0))
                i = int(found[0])+1 if len(found) else None
                output.append(dict(slot=row["slot"], placement_id=row["placement_id"], model=row["model"],
                    size_parameter=row["size_parameter"], gap=row["gap"], classification=row["classification"], control=row["control"],
                    contour_fraction=fraction, lower=None if i is None else float(times[i-1]),
                    upper=None if i is None else float(times[i]),
                    status="observed_saved_two_region_bracket" if i is not None else "unavailable_no_continuous_sampled_two_region_crossing",
                    convention="Integrate original excess only above each region's own local-peak contrast contour; watershed basins and common global floor unchanged"))
    return output


def B_evidence(run=B):
    plan = legacy.verify(run)
    auditpath = run/"results/all_independent_native_audit.json"
    audit = read(auditpath)
    assert audit["complete"] and audit["verified"] == len(plan["jobs"])
    rows = []
    for job in plan["jobs"]:
        placement = plan["placements"][job["placement_id"]]
        statuspath = run/"workers"/f'{job["slot"]:04d}'/"status.json"
        record = read(statuspath)
        assert record["status"] == "complete"
        row = dict(slot=job["slot"], placement_id=job["placement_id"], model=placement["model"],
            size_parameter=placement["size_parameter"], target_ratio=job["target_ratio"], control=job.get("control"),
            initial_gap=placement["edge_gap_10"], seed_status=record["seed_status"],
            primary_claim_eligible=record["B_distance_claim_eligible"],
            seed_lower_bound=record.get("seed_lower_bound"), native_status=str(statuspath), native_status_sha256=sha(statuspath))
        if record.get("withdrawal_time") is None:
            rows.append(row); continue
        assert sha(run/record["exact_rhs_archive"]) == record["exact_rhs_sha256"]
        with np.load(run/record["exact_rhs_archive"]) as z:
            state = z["state"].copy()
        C = float(state[-1]) if placement["model"] == "spring" else placement["size_parameter"]
        row.update(withdrawal_time=record["withdrawal_time"], withdrawal_circumference=C,
            actual_ratio=record["actual_withdrawal_ratio"], withdrawal_topology=record["withdrawal_spatial"]["state_class"],
            gap_status=record["measured_withdrawal_gap"]["reason"], winner=record["final_winner"]["winner"],
            finite_tail_pass=record["finite_tail"]["both_late_windows_steady"])
        gaps = record["measured_withdrawal_gap"]
        if gaps["available"]:
            descriptors = gaps["region_descriptors"]
            inc, chal = descriptors[gaps["incumbent_index"]], descriptors[gaps["challenger_index"]]
            delta = gaps["center_distance"]
            np.testing.assert_allclose(delta-gaps["measured_gap"], (inc["width"]+chal["width"])/2, atol=1e-10)
            signed_angle = 2*np.pi*float(geometry.arc(chal["midpoint"], inc["midpoint"], C))/C
            row.update(measured_gap=gaps["measured_gap"], measured_center_distance=delta,
                measured_long_center_distance=gaps["long_center_distance"],
                measured_long_edge_gap=gaps["long_center_distance"]-(inc["width"]+chal["width"])/2,
                incumbent_10_percent_width=inc["width"], challenger_10_percent_width=chal["width"],
                potential_angular_step=abs(signed_angle), potential_signed_angular_step=signed_angle,
                angle_definition="Between actual global10%-contour region midpoints at withdrawal; potential challenger win step, not assumed realized heading change")
        raw = record["B_resolution"] if row["primary_claim_eligible"] else record["B_resolution"]["raw_resolution"]
        resolution = raw["primary"]
        row.update(resolution_status=resolution["status"], resolution_lower=resolution["lower"], resolution_upper=resolution["upper"],
            one_region_persists_through_cap=raw["one_region_persists_through_cap"],
            primary_resolution_definition="First single resolved peak; contour coalescence retained separately")
        for fraction in (.1, .25, .5, .75):
            threshold = record["gap_threshold_sensitivity"][str(fraction)]
            row[f"gap_{fraction}"] = threshold.get("measured_gap")
            row[f"gap_status_{fraction}"] = threshold["reason"]
        if row["primary_claim_eligible"]:
            angle = row["potential_angular_step"]
            lo, hi = resolution["lower"], resolution["upper"]
            row.update(potential_angular_rate_lower=0. if hi is None else angle/hi,
                potential_angular_rate_upper=None if lo == 0 else angle/lo,
                angular_rate_status="upper_bound_from_resolution_lower_bound" if hi is None else "interval_from_native_saved_resolution_bracket")
        assert sha(run/record["readouts"]) == record["readouts_sha256"]
        with np.load(run/record["readouts"]) as z:
            times = z["time"]
            vector = z["vector_real"]+1j*z["vector_imag"]
            withdraw_index = int(np.argmin(abs(times-record["withdrawal_time"])))
            assert abs(times[withdraw_index]-record["withdrawal_time"]) < 1e-7
            endpoint = len(times)-1 if resolution["upper"] is None else int(np.argmin(abs(times-record["withdrawal_time"]-resolution["upper"])))
            if resolution["upper"] is not None:
                assert abs(times[endpoint]-record["withdrawal_time"]-resolution["upper"]) < 1e-7
            def turn(first, second):
                return float(geometry.arc(np.angle(second), np.angle(first), 2*np.pi)) if min(abs(first), abs(second)) >= .05 else None
            row.update(vector_magnitude_basal=float(abs(vector[0])), vector_magnitude_withdrawal=float(abs(vector[withdraw_index])),
                vector_magnitude_resolution_or_cap=float(abs(vector[endpoint])),
                vector_turn_basal_to_resolution_or_cap=turn(vector[0], vector[endpoint]),
                vector_turn_withdrawal_to_resolution_or_cap=turn(vector[withdraw_index], vector[endpoint]),
                vector_turn_basal_to_final=turn(vector[0], vector[-1]),
                vector_turn_definition="Signed shortest endpoint heading change; magnitude>=0.05 at both endpoints. No unwrapping through cancellation; basal and withdrawal references retained separately.")
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--destination", type=Path, default=DESTINATION)
    parser.add_argument("--include-A2", action="store_true")
    parser.add_argument("--include-B", action="store_true")
    args = parser.parse_args()
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    A1 = A1_evidence()
    write_csv(destination/"A1_complete_region_evidence.csv", A1)
    write_csv(destination/"A1_fixed_angle_geometry_remap.csv", fixed_angle_remaps(A1))
    write_csv(destination/"A1_control_handover_convergence.csv", control_handovers(A1))
    write_csv(destination/"A1_threshold_handover_sensitivity.csv", threshold_handovers(A1))
    if args.include_A2:
        write_csv(destination/"A2_all_switching_brackets.csv", A2_brackets())
    if args.include_B:
        write_csv(destination/"B_complete_region_evidence_and_geometry.csv", B_evidence())
        write_csv(destination/"B_control_evidence_and_geometry.csv", B_evidence(ROOT/"data/runs/signature_figure3_regional_B_controls_20261020"))
    write(destination/"generation.json", dict(generated_at=now(), A1_cases=len(A1),
        source_refinement_audit_sha256=sha(REFINED/"results/independent_refinement_audit.json"),
        source_native_primary_audit_sha256=sha(PRIMARY/"results/primary_independent_native_audit.json"),
        native_integrations=0, source_code_sha256=sha(Path(__file__)), A2_included=args.include_A2, B_included=args.include_B,
        B_audit_sha256=sha(B/"results/all_independent_native_audit.json") if args.include_B else None,
        stage="Audited native evidence and support-limited measured-geometry remaps; unavailable/censored/topology outcomes retained"))
    print(json.dumps(dict(destination=str(destination), A1_cases=len(A1))))


if __name__ == "__main__":
    main()
