"""Read-only regional proposal audit against retained complete A1 fields.

This does not alter the approved native trajectories or calibrate A2/B. Saved
sample brackets expose whether a continuous two-region handover exists; a
single-region reorientation remains a separate measured outcome.
"""
from pathlib import Path
import argparse
import csv
import json

import numpy as np
import signature_figure3_native_20261020 as native
import signature_figure3_regions_20261020 as regions
import signature_figure3_campaign_20261020 as campaign
from signature_figure3_analysis_20261020 import write_csv

read, write, sha, now = native.read, native.write, native.sha, native.now


def analyze_case(run, job, plan, destination, stage="primary"):
    statuspath = run/"workers"/stage/f'{job["slot"]:04d}'/"status.json"
    record = read(statuspath)
    # The worker records the actual mirrored cue axis. Reusing the manifest's
    # unmirrored placement would silently mislabel symmetry-control regions.
    placement = record["placement"]
    adapter = native.Adapter(placement["model"], placement["size_parameter"], placement["cells"])
    if record["status"] != "complete":
        return dict(placement_id=placement["placement_id"], available=False, reason=record["status"])
    old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
    cue = placement["cue_center_fraction"]
    arrays = {key: [] for key in ("time", "regional_dominance", "allocated_region_dominance", "whole_field_dominance",
        "incumbent_excess", "challenger_excess", "two_region_available", "resolved_peaks", "full_excess_amount",
        "vector_magnitude", "vector_heading", "state_class", "dominance_025", "dominance_050", "dominance_075")}
    baseline_excess, event_records = None, []
    previous_d, previous_time, last_state_class = None, None, None
    for phase_index, phase in enumerate(record["phases"]):
        archive = run/phase["path"]
        assert sha(archive) == phase["sha256"]
        with np.load(archive) as z:
            tt, yy = z["time"].copy(), z["states"].copy()
        for index in range(0 if phase_index == 0 else 1, len(tt)):
            time_value = phase["absolute_start"]+float(tt[index])
            state = yy[index]
            C = adapter.circumference(state)
            value = regions.readout(adapter.primary(state), C, old, cue)
            if baseline_excess is None:
                baseline_excess = value["full_excess_amount"]
            two = value["two_region_available"]
            d = value["regional_dominance"]
            if two and previous_d is None and d >= 0:
                event_records.append(dict(event="two_region_birth_already_challenger_dominant",
                    time=time_value, dominance=d,
                    note="No continuous handover bracket may be inferred across an unavailable one-region frame"))
            if two and previous_d is not None and previous_d < 0 <= d:
                event_records.append(dict(event="continuous_two_region_handover_saved_bracket", lower=previous_time,
                    upper=time_value, lower_dominance=previous_d, upper_dominance=d))
            if last_state_class != value["state_class"]:
                event_records.append(dict(event="spatial_state_transition", time=time_value,
                    previous=last_state_class, current=value["state_class"],
                    two_region_handover=False))
            previous_d, previous_time = (d if two else None), time_value
            last_state_class = value["state_class"]
            vector = value["global_vector"]
            scalar_values = dict(time=time_value, regional_dominance=d,
                allocated_region_dominance=value.get("allocated_region_dominance"),
                whole_field_dominance=value["whole_field_diagnostic"].get("dominance"),
                incumbent_excess=value.get("incumbent_excess"), challenger_excess=value.get("challenger_excess"),
                two_region_available=two, resolved_peaks=value["resolved_peaks"],
                full_excess_amount=value["full_excess_amount"], vector_magnitude=vector["magnitude"],
                vector_heading=vector["heading_radians"], state_class=value["state_class"],
                dominance_025=value.get("contour_sensitivity", {}).get("0.25", {}).get("dominance"),
                dominance_050=value.get("contour_sensitivity", {}).get("0.5", {}).get("dominance"),
                dominance_075=value.get("contour_sensitivity", {}).get("0.75", {}).get("dominance"))
            for key, item in scalar_values.items():
                arrays[key].append(item if item is not None else np.nan)
    destination.mkdir(parents=True, exist_ok=True)
    stem = f'{stage}_{job["slot"]:04d}_placement_{placement["placement_id"]:03d}'
    archive = destination/f'{stem}_region_traces.npz'
    np.savez_compressed(archive, **{k: np.asarray(v) for k, v in arrays.items()})
    events = [e for e in event_records if e["event"] == "continuous_two_region_handover_saved_bracket"]
    dominant_births = [e for e in event_records if e["event"] == "two_region_birth_already_challenger_dominant"]
    if events:
        handover_status = "observed_saved_native_two_region_bracket"
    elif dominant_births:
        handover_status = "unavailable_continuous_crossing_at_two_region_birth"
    elif arrays["two_region_available"][-1]:
        handover_status = "right_censored_at_completed_native_horizon"
    else:
        handover_status = "unavailable_two_region_topology_at_final_state"
    result = dict(generated_at=now(), slot=job["slot"], control=job.get("control"),
        placement_id=placement["placement_id"], model=placement["model"],
        size_parameter=placement["size_parameter"], gap=placement["edge_gap_10"],
        classification=placement["classification_10"], source_status=str(statuspath), source_status_sha256=sha(statuspath),
        old_fixed_footprint_handover=record.get("handover_time"),
        complete_region_handover_first_bracket=events[0] if events else None,
        complete_region_handover_events=len(events), two_region_frames=int(np.sum(arrays["two_region_available"])),
        complete_region_handover_status=handover_status,
        complete_region_handover_lower_bound=float(arrays["time"][-1]) if handover_status == "right_censored_at_completed_native_horizon" else None,
        sampled_frames=len(arrays["time"]), final_state_class=arrays["state_class"][-1],
        final_allocated_dominance=arrays["allocated_region_dominance"][-1],
        final_whole_field_dominance=arrays["whole_field_dominance"][-1],
        initial_complete_basal_excess=baseline_excess, traces=str(archive), traces_sha256=sha(archive),
        all_state_transitions=event_records,
        scientific_status="Proposed readout on unchanged retained native fields; no dependent calibration is authorized by this diagnostic")
    # NaN is meaningful in NPZ arrays, while JSON uses explicit null statuses.
    for key in ("final_allocated_dominance", "final_whole_field_dominance"):
        if not np.isfinite(result[key]):
            result[key] = None
    write(destination/f'{stem}_region_proposal.json', result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, default=campaign.DEFAULT)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--include-mirrors", action="store_true")
    parser.add_argument("--placement-id", type=int)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--stage", default="primary")
    args = parser.parse_args()
    run = args.run.resolve()
    plan = campaign.verify(run)
    candidates = [j for j in plan["stages"][args.stage] if j["arm"] == "A1" and (args.include_mirrors or not j.get("control"))]
    if args.placement_id is not None:
        candidates = [j for j in candidates if j["placement_id"] == args.placement_id]
    elif not args.all:
        # One antipodal, external trajectory at each of all twelve basals.
        candidates = [j for j in candidates if plan["placements"][j["placement_id"]]["requested_role"] == "antipode"]
    destination = args.destination.resolve() if args.destination else run/"results/complete_region_proposal"
    rows = []
    for job in candidates:
        result = analyze_case(run, job, plan, destination, args.stage)
        rows.append({k: v for k, v in result.items() if k != "all_state_transitions"})
        print(json.dumps({k: result[k] for k in ("placement_id", "model", "size_parameter", "complete_region_handover_events", "final_state_class")}), flush=True)
    csv_rows = [{k: json.dumps(v) if isinstance(v, dict) else v for k, v in row.items()} for row in rows]
    name = "all_native_and_mirror_summary.csv" if args.all and args.include_mirrors else "all_case_summary.csv" if args.all else "proposal_sample_summary.csv"
    write_csv(destination/name, csv_rows)
    write(destination/"last_reanalysis.json", dict(generated_at=now(), cases=len(rows),
        native_integrations=0, source_manifest_sha256=sha(run/"manifest.json"),
        region_source_sha256=sha(Path(regions.__file__)), reanalysis_source_sha256=sha(Path(__file__)),
        scientific_status="Definition proposal pending review; existing frozen simulation sources unchanged"))


if __name__ == "__main__":
    main()
