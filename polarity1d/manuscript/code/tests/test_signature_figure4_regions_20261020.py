"""Ancestry, physical coordinates and censoring of imaging-style region tracks."""
from pathlib import Path
import sys
import json

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"mechanisms"))
from signature_figure4_regions_20261020 import _complete_peak_basins, _dominant, components, frame_regions, imaging_variants, survival_curve, track_profiles
from signature_figure4_regions_campaign_20261020 import analyze_case, sha


def bump(n=100, center=0., width=.07, amplitude=1.):
    q=(np.arange(n)+.5)/n
    distance=(q-center+.5) % 1-.5
    return .1+amplitude*np.exp(-.5*(distance/width)**2)


def test_periodic_edge_is_one_complete_region():
    mask=np.array([True, True, False, False, True])
    result=components(mask)
    assert len(result)==1
    assert set(result[0])=={0,1,4}
    profile=bump()
    regions=frame_regions(profile,10.,1.)["regions"]
    assert len(regions)==1
    assert abs(np.angle(np.exp(1j*regions[0]["heading"])))<1e-12


def test_weak_profile_has_no_qualified_regions_even_with_relative_contour_peak():
    result=frame_regions(bump(amplitude=.01),10.,1.)
    assert result["low_contrast"] and result["regions"]==[]
    with pytest.raises(ValueError):
        frame_regions(bump(),10.,float("nan"))


def test_full_component_amount_and_width_scale_with_actual_circumference():
    a=frame_regions(bump(center=.3),10.,1.)["regions"][0]
    b=frame_regions(bump(center=.3),20.,1.)["regions"][0]
    assert b["width"]==pytest.approx(2*a["width"])
    assert b["excess"]==pytest.approx(2*a["excess"])
    assert b["width_fraction"]==a["width_fraction"]


def test_translating_region_keeps_identity_across_the_ring_seam():
    profiles=np.array([bump(center=c) for c in [.98, .99, .0, .01, .02]])
    result=track_profiles(np.arange(5.),profiles,10.,basal_contrast=1.)
    assert len(result["tracks"])==1 and not result["births"]
    assert result["tracks"][0]["left_censored"]
    assert result["tracks"][0]["right_censored"]
    assert result["survival"]["observed"]["incident_tracks"]==0


def test_newborn_takeover_needs_two_observed_frames_above_three_to_one():
    old=bump(center=.2,amplitude=1.)-.1
    new=bump(center=.7,amplitude=1.)-.1
    profiles=np.array([.1+old,.1+old+new,.1+.1*old+new,.1+.1*old+new])
    result=track_profiles(np.arange(4.),profiles,10.,basal_contrast=1.)
    assert len(result["births"])==1
    assert result["births"][0]["opposite_hemisphere_birth"]
    assert len(result["takeovers"])==1
    assert result["takeovers"][0]["time"]==3.
    assert result["takeovers"][0]["first_qualifying_frame_time"]==2.
    assert result["takeovers"][0]["dominance"]>.5


def test_initial_tracks_are_not_counted_as_new_births_or_incident_lifetimes():
    profiles=np.array([bump(center=.2),bump(center=.2),np.full(100,.1)])
    result=track_profiles(np.arange(3.),profiles,10.,basal_contrast=1.)
    assert result["summary"]["observed_births"]==0
    assert result["tracks"][0]["end_reason"]=="low_contrast_unresolved"
    assert result["tracks"][0]["right_censored"]
    assert result["survival"]["observed"]["excluded_unknown_birth"]==1


def test_birth_rate_uses_physical_perimeter_time_exposure():
    first=bump(center=.2)-.1;second=bump(center=.7)-.1
    profiles=np.array([.1+first,.1+first+second,.1+first+second])
    a=track_profiles(np.arange(3.),profiles,[10.,10.,10.],basal_contrast=1.)
    b=track_profiles(np.arange(3.),profiles,[10.,20.,30.],basal_contrast=1.)
    assert a["summary"]["observed_births"]==b["summary"]["observed_births"]==1
    assert b["summary"]["perimeter_time_exposure"]==40.
    assert a["summary"]["births_per_time_per_perimeter"]==pytest.approx(2*b["summary"]["births_per_time_per_perimeter"])


def test_merging_regions_does_not_invent_a_newborn_or_a_takeover():
    old=bump(center=.3,width=.06)-.1
    new=bump(center=.7,width=.06)-.1
    bridge=bump(center=.5,width=.17,amplitude=.9)-.1
    profiles=np.array([.1+old+new,.1+old+new+bridge,.1+old+new+bridge])
    result=track_profiles(np.arange(3.),profiles,10.,basal_contrast=1.)
    assert [r["region_count"] for r in result["frames"]]==[2,1,1]
    assert not result["births"] and not result["takeovers"]
    assert any(event["kind"]=="merge" for event in result["overlap_events"])
    assert any(track["end_reason"]=="merger_identity_loss" and track["right_censored"]
               for track in result["tracks"])


def test_broad_incumbent_below_global_contour_keeps_its_complete_amount():
    C=30.;x=(np.arange(600)+.5)*C/600
    def triangle(center,radius,height):
        distance=np.abs((x-center+C/2) % C-C/2)
        return height*np.maximum(1-distance/radius,0.)
    old=triangle(6.,3.,1.)
    new=triangle(21.,.4,1.2)
    profiles=np.array([.2+old,.2+.45*old+new,.2+.45*old+new])
    result=track_profiles(np.arange(3.),profiles,C,basal_contrast=1.)
    assert [f["region_count"] for f in result["frames"]]==[1,1,1]
    assert result["frames"][1]["latent_peak_basins"]==1
    assert result["frames"][1]["regions"][0]["excess"]==pytest.approx(.48,abs=1e-12)
    assert not result["takeovers"]
    # The previous implementation invented E_old=0 at contour disappearance
    # and returned D=1. Actual complete amounts are1.35 versus0.48.
    assert (.48-1.35)/(.48+1.35)<0


def test_contour_disappearance_without_resolved_native_peak_is_not_zero_amount():
    old=bump(center=.2,amplitude=1.)-.1
    new=bump(center=.7,amplitude=1.)-.1
    profiles=np.array([.1+old,.1+new,.1+new,.1+new])
    result=track_profiles(np.arange(4.),profiles,10.,basal_contrast=1.)
    assert result["births"] and not result["takeovers"]
    assert any(t["takeover_unscorable_reason"] is not None for t in result["tracks"])


def test_pairwise_win_with_a_larger_third_region_does_not_preempt_later_global_takeover():
    old=bump(center=.1,width=.035)-.1
    new=bump(center=.43,width=.035)-.1
    third=bump(center=.76,width=.035)-.1
    profiles=np.array([.1+old,.1+old+.8*new,
        .1+.1*old+.65*new+1.2*third,.1+.1*old+.65*new+1.2*third,
        .1+.1*old+1.1*new+.8*third,.1+.1*old+1.1*new+.8*third])
    result=track_profiles(np.arange(6.),profiles,20.,basal_contrast=1.)
    newborn=result["births"][0]["track_id"]
    assert next(e for e in result["pairwise_takeovers"] if e["newborn_track_id"]==newborn)["time"]==3.
    actual=next(e for e in result["takeovers"] if e["newborn_track_id"]==newborn)
    assert actual["time"]==5. and actual["first_qualifying_frame_time"]==4.
    assert actual["globally_dominant_required"]


def test_native_peak_basins_partition_all_excess_and_active_flat_saddles_are_ambiguous():
    profile=np.array([0.,2.,1.,1.,3.,0.,0.,0.])
    atoms=_complete_peak_basins(profile,8.)
    assert sum(a["excess"] for a in atoms)==pytest.approx(profile.sum())
    assert len(atoms)==2 and all(not a["amount_available"] for a in atoms)
    assert _dominant([dict(excess=5.,amount_available=False),
                      dict(excess=1.,amount_available=True)]) is None


def test_complete_basin_partition_is_invariant_to_periodic_grid_rotation():
    p=bump(center=.17,width=.09)+.6*(bump(center=.69,width=.04)-.1)
    a=sorted(r["excess"] for r in _complete_peak_basins(p,17.))
    b=sorted(r["excess"] for r in _complete_peak_basins(np.roll(p,37),17.))
    assert a==pytest.approx(b,abs=1e-12)


def test_kaplan_meier_retains_right_censoring_and_observation_interval_sensitivity():
    tracks=[dict(left_censored=False,right_censored=False,duration=2.,duration_lower=1.,duration_upper=3.),
            dict(left_censored=False,right_censored=True,duration=4.,duration_lower=4.,duration_upper=5.),
            dict(left_censored=True,right_censored=True,duration=10.)]
    curve=survival_curve(tracks)
    assert curve["incident_tracks"]==2
    assert curve["curve"][1]["survival"]==.5
    assert curve["curve"][-1]["survival"]==.5
    assert survival_curve(tracks,event_clock="lower")["curve"][1]["time"]==1.
    assert survival_curve(tracks,event_clock="upper")["curve"][1]["time"]==3.


def test_optical_smoothing_preserves_periodic_mean_and_uses_physical_length():
    p=bump(center=.03)
    raw=frame_regions(p,10.,1.)
    blurred=frame_regions(p,10.,1.,smoothing=.5)
    assert blurred["smoothed_mean"]==pytest.approx(raw["smoothed_mean"],abs=1e-14)
    assert blurred["contrast"]<raw["contrast"]


def test_imaging_variants_keep_analysis_end_and_process_basal_with_same_optics():
    p=bump();profiles=np.tile(p,(6,1))
    variants=list(imaging_variants(np.arange(6.),profiles,10.,basal_contrast=float(np.ptp(p)),
        kernel_length=.5,basal_profile=p,basal_circumference=10.,thresholds=(.5,),downsampling=(4,)))
    assert len(variants)==2
    assert [f["time"] for f in variants[0]["result"]["frames"]]==[0.,4.,5.]
    assert variants[1]["observed_basal_contrast"]<variants[0]["observed_basal_contrast"]
    assert variants[1]["result"]["summary"]["low_contrast_time_fraction"]==0


def test_case_adapter_hashes_exact_basals_and_emits_serializable_common_imaging_view(tmp_path):
    native=tmp_path/"native"
    native.mkdir()
    p=bump()
    t=np.arange(0.,401.,40.)
    np.savez(native/"profiles.npz",time=t,primary=np.tile(p,(len(t),1)),
             circumference=10.,material_positions=(np.arange(len(p))+.5)/len(p))
    np.savez(native/"basal.npz",primary=p,circumference=10.)
    for name in ("scalar.json","checkpoint.json","status.json"):
        (native/name).write_text("{}\n")
    case=dict(id="test",model="wave_pinning",native_cells=len(p),
              profile_cadence=40.,native_horizon=400.,actualC0=10.,
              basal_contrast=float(np.ptp(p)))
    for key,name in (("profile","profiles.npz"),("scalar","scalar.json"),
                     ("basal_profile","basal.npz"),("basal_checkpoint","checkpoint.json"),
                     ("source_status","status.json")):
        case[key+"_path"]=str(native/name)
        case[key+"_sha256"]=sha(native/name)
    record=analyze_case(case,tmp_path/"outputs",inventory_sha256="test_inventory")
    persisted=json.loads((tmp_path/"outputs"/"result.json").read_text())
    assert persisted["passed"] and record["primary"]["downsample_factor"]==2
    assert record["primary"]["smoothed"] and len(record["sensitivities"])==24
    assert record["actual_analysis_window"]==[200.,400.]
    with pytest.raises(FileExistsError):
        analyze_case(case,tmp_path/"outputs",inventory_sha256="test_inventory")
