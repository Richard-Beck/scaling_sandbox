"""Direct native global-heading checks for proposed Figure 3 translations.

R=|M1|/M0 is instantaneous coherence. It is not Figure 4's basal-M1-normalized
vector magnitude. Every direction-quality gap breaks angle unwrapping.
"""
from pathlib import Path
import argparse
import json

import numpy as np

import signature_figure3_native_20261020 as native
import signature_figure3_campaign_20261020 as legacy
import signature_figure3_geometry_20261020 as geometry
import signature_figure3_synthesis_20261020 as synthesis
from signature_figure3_analysis_20261020 import write_csv


def direction_intervals(vector, floor):
    """Usable endpoint directions with a conservative linear-field chord gate.

    Opposite usable endpoints can straddle cancellation between saved frames.
    Their straight complex-moment chord must also stay above the gate before
    any interpolation/unwrapping across that interval is accepted.
    """
    vector = np.asarray(vector)
    R = abs(vector)
    eligible = R > floor if floor == 1e-10 else R >= floor
    delta = np.diff(vector)
    denominator = abs(delta)**2
    projection = np.divide(-np.real(vector[:-1]*np.conj(delta)), denominator,
                           out=np.zeros_like(denominator), where=denominator > 0)
    minimum = abs(vector[:-1]+np.clip(projection, 0., 1.)*delta)
    chord = minimum > floor if floor == 1e-10 else minimum >= floor
    return eligible, eligible[:-1] & eligible[1:] & chord, minimum


def target_bracket(times, vector, target, sign, floor):
    """Saved-field heading crossing, with cancellation intervals left open."""
    times, vector = np.asarray(times), np.asarray(vector)
    eligible, intervals, minimum = direction_intervals(vector, floor)
    if target <= 1e-10 and eligible[0]:
        return dict(status="already_at_target_at_start", lower=0., upper=0., prior_heading_gap=False)
    starts = np.flatnonzero(eligible & ~np.r_[False, intervals])
    stops = np.flatnonzero(eligible & ~np.r_[intervals, False])+1
    phi0 = float(np.angle(vector[0]))
    for start, stop in zip(starts, stops):
        progress = sign*np.unwrap(np.angle(vector[start:stop]*np.exp(-1j*phi0)))
        if start > 0 and progress[0] >= target:
            return dict(status="target_attained_when_heading_returns_after_gap", lower=None, upper=None,
                first_available_attainment_time=float(times[start]),
                unobserved_direction_gap_last_sample=float(times[start-1]), prior_heading_gap=True,
                note="No crossing time is inferred across a coherence/direction-quality gap")
        found = np.flatnonzero((progress[:-1] < target) & (progress[1:] >= target))
        if len(found):
            index = int(start+found[0])+1
            return dict(status="observed_saved_valid_segment_crossing", lower=float(times[index-1]), upper=float(times[index]),
                prior_heading_gap=bool((~intervals[:index]).any()),
                first_global_crossing_proven=not bool((~intervals[:index]).any()))
    if (~intervals).any():
        return dict(status="no_observed_target_crossing_with_heading_gaps", lower=None, upper=None,
            completed_horizon=float(times[-1]), prior_heading_gap=True,
            note="The complete native horizon is retained; undefined/gated direction intervals prevent a global timing lower bound")
    return dict(status="right_censored_saved_heading_at_completed_horizon", lower=float(times[-1]), upper=None,
                prior_heading_gap=False, scope="No crossing among complete retained native field samples")


def analyze(destination):
    A1 = synthesis.A1_evidence()
    rows, targets = [], []
    for row in A1:
        status = native.read(Path(row["native_status"]))
        placement = status["placement"]
        archive = synthesis.PRIMARY/status["readouts"]
        assert native.sha(archive) == status["readouts_sha256"]
        with np.load(archive) as z:
            times = z["time"].copy()
            vector = z["vector_real"]+1j*z["vector_imag"]
        phi0 = float(np.angle(vector[0]))
        sign = -1. if status["job"].get("mirrored") else 1.
        theta = sign*float(geometry.arc(2*np.pi*placement["cue_center_fraction"], phi0, 2*np.pi))
        if abs(placement["center_distance"]) <= 1e-9:
            theta = 0.
        if abs(abs(theta)-np.pi) < 1e-9:
            theta = np.pi
        if abs(theta) < 1e-10:
            theta = 0.
        assert theta >= 0
        result = dict(slot=row["slot"], placement_id=row["placement_id"], model=row["model"],
            size_parameter=row["size_parameter"], classification=row["classification"], control=row["control"],
            gap=row["gap"], cue_axis_target_angle=theta, native_handover_time=row["handover_time"],
            handover_status=row["handover_status"], angle_at_handover=None, coherence_at_handover=None,
            angle_at_native_cap=sign*float(geometry.arc(np.angle(vector[-1]), phi0, 2*np.pi)) if abs(vector[-1]) > 1e-10 else None,
            native_cap=float(times[-1]), minimum_sampled_coherence=float(abs(vector).min()),
            final_coherence=float(abs(vector[-1])), native_vector_archive=str(archive), native_vector_archive_sha256=native.sha(archive))
        if row["handover_time"] is not None:
            refined = native.read(Path(row["refinement_status"]))
            v = refined["event_readout"]["global_vector"]
            R = v["magnitude"]
            result.update(coherence_at_handover=R,
                angle_at_handover=sign*float(geometry.arc(v["heading_radians"], phi0, 2*np.pi)) if v["direction_available"] else None,
                handover_angle_fraction_of_requested=None if not v["direction_available"] or theta <= 1e-10 else sign*float(geometry.arc(v["heading_radians"], phi0, 2*np.pi))/theta)
        for floor in (1e-10, .05, .1, .2):
            eligible, intervals, minimum = direction_intervals(vector, floor)
            suffix = "direction" if floor == 1e-10 else str(floor)
            result["saved_frame_coverage_"+suffix] = float(np.mean(eligible))
            result["linear_chord_guarded_interval_time_coverage_"+suffix] = float(np.sum(np.diff(times)*intervals)/times[-1])
            result["heading_gap_intervals_"+suffix] = int(np.sum(~intervals))
            for fraction in (1., .5, .9, .99):
                value = target_bracket(times, vector, fraction*theta, sign, floor)
                targets.append(dict(slot=row["slot"], placement_id=row["placement_id"], model=row["model"],
                    size_parameter=row["size_parameter"], classification=row["classification"], control=row["control"],
                    gap=row["gap"], coherence_floor=floor, target_fraction=fraction,
                    target_angle=fraction*theta, **value))
        rows.append(result)
    destination.mkdir(parents=True, exist_ok=True)
    write_csv(destination/"A1_direct_heading_check.csv", rows)
    write_csv(destination/"A1_heading_target_time_brackets.csv", targets)
    native.write(destination/"heading_audit.json", dict(generated_at=native.now(), native_cases=len(rows),
        target_gate_combinations=len(targets), native_integrations=0, source_code_sha256=native.sha(Path(__file__)),
        coherence_definition="R=abs(first circular activity moment)/total instantaneous active amount; different from basal-M1-normalized Figure4 P",
        crossing_scope="Retained native field sample brackets only; cancellation breaks all unwrapping, including any intervening straight complex-moment chord below the quality floor. This conservative interpolation guard is not an asserted native dense trajectory. Full-angle target primary;50/90/99% sensitivities retained.",
        raw_direction_floor=1e-10, reporter_quality_floor=.05, other_quality_floors=[.1,.2]))
    return dict(cases=len(rows), target_gate_combinations=len(targets))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--destination", type=Path, default=synthesis.DESTINATION)
    args = parser.parse_args()
    print(json.dumps(analyze(args.destination.resolve())))


if __name__ == "__main__":
    main()
