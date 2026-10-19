"""Direct-field readouts for the separate common-initial-arc stationary arm.

No native integrations or figure rendering; all exports replace atomically.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import fcntl
import io
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import revision_initial_arc_sizes_20261019 as producer
import revision_standard_angles_duration_20261019 as base
import revision_standard_angles_report_20261019 as prior
import revision_stationary_characterization_20261019 as rich

RUN=producer.DEFAULT
OUT=base.TOPIC.parent/"standard_arclength_sizes"
SOURCE={}
TERMINAL=False


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def record(path):
    path=Path(path);SOURCE[str(path)]=sha(path);return path
def dest(name,group="data"):
    parent=OUT/"supporting"/group;parent.mkdir(parents=True,exist_ok=True);return parent/name
def jsonwrite(name,value):
    target=dest(name,"metadata");temporary=target.with_suffix(target.suffix+".tmp");temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");temporary.replace(target)
def csvwrite(name,frame):
    target=dest(name);temporary=target.with_suffix(target.suffix+".tmp");frame.to_csv(temporary,index=False);temporary.replace(target)


def protocol_key(case):
    """Identify executed protocols, not labels for equivalent coordinates."""
    fields=[case["model"],case["input_sha256"],float(case["cue_duration"]).hex(),
            float(case["recovery"]).hex(),float(case["full_physical_support"]).hex(),
            float(case["requested_peak"]).hex()]
    if case["model"]=="spring":
        fields += [float(case["resting_circumference"]).hex(),float(case["cue_centre_fraction"]).hex()]
    else:
        fields += [case["initial_cue_shape_sha256"]]
    return hashlib.sha256(json.dumps(fields,separators=(",",":")).encode()).hexdigest()


def reconcile_arms(manifest,frame,previous):
    oldmanifest=json.loads(record(base.ROOT/"data/runs/revision_standard_angles_sizes_20261019/manifest.json").read_text())
    oldcases=oldmanifest["cases"];oldkeys={c["condition_id"]:protocol_key(c) for c in oldcases}
    newkeys={c["condition_id"]:protocol_key(c) for c in manifest["cases"]}
    common=set(oldkeys.values())&set(newkeys.values())
    assert len(set(oldkeys.values()))==len(set(newkeys.values()))==240
    assert len(common)==125 and len(set(oldkeys.values())|set(newkeys.values()))==355
    allrows=[]
    for arm,table,keys in (("fixed_angle",previous,oldkeys),("common_initial_physical_arc",frame,newkeys)):
        for row in table.to_dict("records"):
            row.update(experiment_arm=arm,unique_protocol_id=keys[int(row["condition_id"])],arm_row_id=f'{arm}:{int(row["condition_id"])}')
            allrows.append(row)
    union=pd.DataFrame(allrows);assert len(union)==480
    csvwrite("stationary_two_arms_all480.csv",union)
    identities=[]
    for key,group in union.groupby("unique_protocol_id",sort=False):
        identities.append(dict(unique_protocol_id=key,model=group.iloc[0]["model"],size_parameter=group.iloc[0]["size_parameter"],
            cue_duration=group.iloc[0]["cue_duration"],arm_rows=len(group),arm_row_ids=";".join(group.arm_row_id),
            complete_arm_rows=int((group.status=="complete").sum())))
    csvwrite("stationary_two_arms_unique355.csv",pd.DataFrame(identities))
    pairs=[];reference_fractions=(.0625,.125,.25,.5)
    indices=sorted({int(c["initial_arc_grid_index"]) for c in manifest["cases"]});assert len(indices)==4
    fractions_by_index=dict(zip(indices,reference_fractions))
    for case in manifest["cases"]:
        # Pair the design levels that coincide at chemical C20 / Spring C0=1.
        # Exact shared protocols at other levels have separate unique IDs above.
        arc=frame[frame.condition_id==case["condition_id"]].iloc[0]
        fraction=fractions_by_index[int(case["initial_arc_grid_index"])]
        candidates=previous[(previous.model==case["model"]) & (previous.size_parameter==case["size_parameter"]) &
                            (previous.cue_duration==case["cue_duration"]) & np.isclose(previous.fraction,fraction,rtol=0,atol=1e-15)]
        assert len(candidates)==1
        angle=candidates.iloc[0]
        pair=dict(model=case["model"],size_parameter=case["size_parameter"],cue_duration=case["cue_duration"],
                  reference_level_index=case["initial_arc_grid_index"],arc_condition_id=case["condition_id"],
                  angle_condition_id=int(angle.condition_id),same_executed_protocol=newkeys[case["condition_id"]]==oldkeys[int(angle.condition_id)])
        for label,row in (("angle_arm",angle),("arc_arm",arc)):
            for key in ("status","physical_distance_at_onset","requested_angle_deg","target_angle_in_P_coordinates_deg","fraction","circumference","original_core_physical_width",
                        "cue_support_core_intersection","overlap_fraction_cue_support","overlap_fraction_original_core","initial_overlap",
                        "initial_cue_weighted_active","initial_cue_weighted_inactive","R_withdrawal","R_final","angle_withdrawal_raw_deg",
                        "angle_final_raw_deg","angle_withdrawal_deg","angle_final_deg","fraction_requested_turn_final","target_error_final_deg",
                        "M0_ratio_withdrawal","M0_ratio_final","q1_withdrawal","q1_final","fixed_region_excess_withdrawal","fixed_region_excess_final",
                        "full_field_tail_gate","full_field_tail_status","actual_cue_path","actual_cue_sha256","actual_recovery_path","actual_recovery_sha256"):
                pair[label+"_"+key]=row.get(key)
        pairs.append(pair)
    csvwrite("stationary_two_arms_reference_level_pairs.csv",pd.DataFrame(pairs))
    return dict(arm_rows=480,unique_protocol_conditions=355,shared_protocol_conditions=125,reference_level_pairs=240,
                same_protocol_reference_level_pairs=sum(p["same_executed_protocol"] for p in pairs))


def write_readme(generated,complete):
    text=f"""# Common-initial-arc stationary-cue arm

Generated: {generated}. Complete: {complete}/240 requested conditions. No overlap is excluded.

Chemical rings use C = 20/30/40 and initial arc offsets 1.25/2.5/5/10. Spring uses resting C₀ = 1/2/3 and the same four physical initial offsets defined by the exact qualified C₀=1 basal circumference times 1/16, 1/8, 1/4 and 1/2. Cue centres remain stationary in material coordinates; their physical offset subsequently scales with Spring's instantaneous circumference. Native supports, peaks, durations, basals and 13,000-unit cue-free continuation retain the first arm's definitions. Distance and angle remain linked at each fixed circumference.

- [Current 240 requested conditions and available readouts](supporting/data/stationary_arclength_summary.csv)
- [Initial context, overlap and exact cue shapes](supporting/data/stationary_arclength_initial_context.csv)
- [Initial native profiles and cue samples](supporting/data/stationary_arclength_initial_profiles_and_cues.csv)
- [Saved regional first-entry and retained-entry brackets](supporting/data/stationary_arclength_crossings.csv)
- [Final native primary profiles](supporting/data/stationary_arclength_final_profiles.csv)
- [Both arms: 480 rows](supporting/data/stationary_two_arms_all480.csv) and [355 unique protocols](supporting/data/stationary_two_arms_unique355.csv)
- [Paired design levels](supporting/data/stationary_two_arms_reference_level_pairs.csv): levels coincide at chemical C20 / Spring C0=1, then compare fixed angle with fixed initial arc across size. Shared protocol IDs separately capture all 125 exact reused conditions, including matches at different design-level indices.
- [Latest collection status](supporting/metadata/latest_status.json), [validation](supporting/metadata/validation.json) and [source hashes](supporting/metadata/source_sha256.json)

Raw M₁ components use fixed physical ring axes. P = M₁/M₁(0) aligns and scales the initial vector once; R = |P|, physical M₀ is total active amount, and q₁ = |M₁|/M₀. Moment quadrature is instantaneous C times the native-cell mean. Endpoint angles are masked at R < 0.1; raw phases remain diagnostics. Compare alignment using wrapped target_error_final_deg and the actual target_angle_in_P_coordinates_deg, which includes the initial moment's measured offset from the preparation axis. Fraction of requested turn uses the raw signed endpoint angle only when R ≥ 0.1; its branch discontinuity near ±180° makes it a raw diagnostic rather than a universal progress or success measure. The fixed original material region subtracts its transported initial background amount and does not track a site's ancestry. Event times are saved-sample brackets, separate for cue-on and release; sparse duration levels do not provide refined minimum cue durations. Full-field tail flags are empirical finite-window excursions; retained primary-only archives keep the full-field flag unavailable.

All compact files replace atomically; latest_status.json is written last and hashes that generation's files. Missing native branches remain pending/error rows. These exports contain no new integration or figure layout. The old angle-arm files are unchanged.
"""
    target=OUT/"README.md";temporary=target.with_suffix(".md.tmp");temporary.write_text(text);temporary.replace(target)


def collect():
    manifest=json.loads(record(RUN/"manifest.json").read_text())
    assert sha(RUN/"manifest.json")==json.loads(record(RUN/"freeze.json").read_text())["manifest_sha256"]
    inventory_source=RUN/"results/condition_inventory.csv";inventory_bytes=inventory_source.read_bytes()
    inventory_snapshot=dest("source_condition_inventory.csv");temporary=inventory_snapshot.with_suffix(".csv.tmp")
    temporary.write_bytes(inventory_bytes);temporary.replace(inventory_snapshot);record(inventory_snapshot)
    inventory=pd.read_csv(io.BytesIO(inventory_bytes))
    assert len(inventory)==240
    for rel,digest in manifest["source_sha256"].items():assert sha(RUN/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(RUN/rel)==digest
    rows=[];events=[];traces={};profiles=[];contexts=[];audits=[];basal_context={};context_cues=set()
    old=pd.read_csv(record(base.TOPIC.parent/"standard_angles_sizes/supporting/data/stationary_sizes_summary.csv"))
    oldevents=pd.read_csv(record(base.TOPIC.parent/"standard_angles_sizes/supporting/data/stationary_sizes_crossings.csv"))
    oldevent_lookup={(int(e["condition_id"]),e["threshold"],e["phase"],e["event"]):e for e in oldevents.to_dict("records")}
    oldtrace=np.load(record(base.TOPIC.parent/"standard_angles_sizes/supporting/data/stationary_sizes_traces.npz"))
    maximum_reuse_replay_error=0.;maximum_shared_trace_error=0.;maximum_shared_geometry_error=0.;maximum_shared_event_error=0.;shared_replays=0;shared_events=0
    for case in manifest["cases"]:
        cid=case["condition_id"];status=inventory[inventory.condition_id==cid].iloc[0]
        row={k:v for k,v in case.items() if k not in ("reuse","fixed_scales")};row["status"]=status["status"]
        row["source_kind"]=case["execution"]
        if case["input"] not in basal_context:
            with np.load(RUN/case["input"]) as z:initial_context=z["state"].copy()
            with np.load(RUN/case["cue_input"]) as z:initial_q=z["material_positions"].copy()
            initial_u=initial_context[:case["cells"]] if case["model"]=="spring" else initial_context[0]
            initial_M0=case["circumference"]*initial_u.mean()
            initial_M1=case["circumference"]*np.mean(initial_u*np.exp(2j*np.pi*initial_q))
            basal_context[case["input"]]=dict(M0_initial=float(initial_M0),M1_initial_real=float(initial_M1.real),M1_initial_imag=float(initial_M1.imag),
                M1_initial_magnitude=float(abs(initial_M1)),q1_initial=float(abs(initial_M1)/initial_M0),initial_active_contrast=float(np.ptp(initial_u)))
        row.update(basal_context[case["input"]])
        if case["cue_input"] not in context_cues:
            with np.load(RUN/case["input"]) as z:initial=z["state"].copy()
            with np.load(RUN/case["cue_input"]) as z:q=z["material_positions"].copy();signal=z["initial_cue"].copy()
            u0=initial[:case["cells"]] if case["model"]=="spring" else initial[0]
            core0=base.core_info(q,u0,case["preparation_axis_fraction"])
            aq=(360*(q-case["preparation_axis_fraction"])+120)%360-120
            for j in range(case["cells"]):
                contexts.append(dict(model=case["model"],size_parameter=case["size_parameter"],size_label=case["size_label"],placement_label=case["placement_label"],
                    material_position=q[j],normal_angle_from_preparation_deg=aq[j],basal_active=u0[j],basal_active_normalized_contrast=(u0[j]-u0.min())/np.ptp(u0),
                    initial_cue=signal[j],initial_cue_unit_peak=signal[j]/case["requested_peak"],fixed_original_core_cell_material_weight=core0["weights"][j]))
            context_cues.add(case["cue_input"])
        if status["status"]!="complete":rows.append(row);continue
        spring=case["model"]=="spring";Cb=case["circumference"];n=case["cells"];T=case["cue_duration"]
        with np.load(RUN/case["input"]) as z:initial=z["state"].copy()
        with np.load(RUN/case["cue_input"]) as z:q=z["material_positions"].copy();signal=z["initial_cue"].copy()
        parts=[];complete_states=[]
        for label in ("cue","recovery"):
            path=record(status["actual_"+label+"_path"]);digest=status["actual_"+label+"_sha256"];assert sha(path)==digest
            row["actual_"+label+"_path"]=str(path);row["actual_"+label+"_sha256"]=digest
            part=prior.load_phase(path,spring,Cb);np.testing.assert_array_equal(part["q"],q);parts.append(part)
            with np.load(path) as z:
                if spring:
                    key="cue_centre_fraction" if "cue_centre_fraction" in z else "cue_center_fraction"
                    assert float(z[key])==case["cue_centre_fraction"] and float(z["cue_peak"])==(case["requested_peak"] if label=="cue" else 0.)
                else:np.testing.assert_array_equal(z["stimulus"],signal if label=="cue" else np.zeros(n))
                if part["full"] is not None:complete_states.append(part["full"])
                else:
                    assert z["snapshot_times"][0]==0. and z["snapshot_times"][-1]==(T if label=="cue" else 13000.)
                    complete_states.append(z["snapshot_states"].copy())
        np.testing.assert_array_equal(complete_states[0][0],initial)
        np.testing.assert_array_equal(complete_states[0][-1],complete_states[1][0])
        assert parts[0]["time"][-1]==T and parts[1]["time"][-1]==13000.
        t=np.r_[parts[0]["time"],parts[1]["time"][1:]+T];joff=len(parts[0]["time"])-1
        u=np.vstack((parts[0]["u"],parts[1]["u"][1:]));C=np.r_[parts[0]["C"],parts[1]["C"][1:]]
        u0=initial[:n] if spring else initial[0];M0=C*u.mean(axis=1);M1=C*(u@np.exp(2j*np.pi*q))/n;M2=C*(u@np.exp(4j*np.pi*q))/n
        assert (M0>0).all();P=M1/M1[0];R=abs(P);angle=prior.heading(P);q1=abs(M1)/M0;q2=abs(M2)/M0;M0ratio=M0/M0[0]
        core=base.core_info(q,u0,case["preparation_axis_fraction"]);w=core["weights"];A=C*(u@w)
        bg0=Cb*core["background"]*w.sum();E0=A[0]-bg0;E=(A-bg0)/E0
        eligible=min(R[joff],R[-1])>=.1;target=case["target_angle_in_P_coordinates_deg"]
        row.update(M0_initial=float(M0[0]),M1_initial_real=float(M1[0].real),M1_initial_imag=float(M1[0].imag),
            M1_withdrawal_real=float(M1[joff].real),M1_withdrawal_imag=float(M1[joff].imag),M1_final_real=float(M1[-1].real),M1_final_imag=float(M1[-1].imag),
            P_withdrawal_real=float(P[joff].real),P_withdrawal_imag=float(P[joff].imag),P_final_real=float(P[-1].real),P_final_imag=float(P[-1].imag),
            M0_withdrawal=float(M0[joff]),M0_final=float(M0[-1]),M0_ratio_withdrawal=float(M0ratio[joff]),M0_ratio_final=float(M0ratio[-1]),
            M0_ratio_sampled_minimum=float(M0ratio.min()),M0_ratio_sampled_maximum=float(M0ratio.max()),
            M0_ratio_sampled_minimum_time=float(t[M0ratio.argmin()]),M0_ratio_sampled_maximum_time=float(t[M0ratio.argmax()]),
            M0_ratio_cue_on_sampled_minimum=float(M0ratio[:joff+1].min()),M0_ratio_cue_on_sampled_maximum=float(M0ratio[:joff+1].max()),
            R_withdrawal=float(R[joff]),R_final=float(R[-1]),R_sampled_minimum=float(R.min()),R_sampled_maximum=float(R.max()),
            R_sampled_minimum_time=float(t[R.argmin()]),R_cue_on_sampled_minimum=float(R[:joff+1].min()),R_release_sampled_minimum=float(R[joff:].min()),
            M0_ratio_at_R_sampled_minimum=float(M0ratio[R.argmin()]),q2_at_R_sampled_minimum=float(q2[R.argmin()]),
            q1_withdrawal=float(q1[joff]),q1_final=float(q1[-1]),q2_withdrawal=float(q2[joff]),q2_final=float(q2[-1]),
            angle_withdrawal_raw_deg=float(angle[joff]),angle_final_raw_deg=float(angle[-1]),
            angle_withdrawal_deg=float(angle[joff]) if R[joff]>=.1 else None,angle_final_deg=float(angle[-1]) if R[-1]>=.1 else None,
            release_signed_turn_deg=float(rich.wrap_error(angle[-1],angle[joff])) if eligible else None,
            release_target_error_change_deg=float(abs(rich.wrap_error(angle[-1],target))-abs(rich.wrap_error(angle[joff],target))) if eligible else None,
            target_error_final_deg=float(rich.wrap_error(angle[-1],target)) if R[-1]>=.1 else None,
            fraction_requested_turn_final=float(angle[-1]/case["requested_angle_deg"]) if R[-1]>=.1 else None,
            fixed_region_excess_initial=1.,fixed_region_excess_withdrawal=float(E[joff]),fixed_region_excess_final=float(E[-1]),
            fixed_region_excess_cue_on_sampled_minimum=float(E[:joff+1].min()),fixed_region_excess_cue_on_sampled_minimum_time=float(t[:joff+1][E[:joff+1].argmin()]),
            fixed_region_excess_release_sampled_minimum=float(E[joff:].min()),fixed_region_excess_release_sampled_minimum_time=float(t[joff:][E[joff:].argmin()]),
            fixed_region_excess_sampled_minimum=float(E.min()),fixed_region_excess_sampled_minimum_time=float(t[E.argmin()]),
            fixed_region_release_net_change=float(E[-1]-E[joff]),fixed_region_rebound_from_sampled_minimum=float(E[-1]-E.min()),
            fixed_region_physical_amount_withdrawal=float(A[joff]),fixed_region_physical_amount_final=float(A[-1]),
            fixed_region_current_background_excess_final=float((C[-1]*((u[-1]-core["background"])@w))/E0),
            physical_C_withdrawal=float(C[joff]),physical_C_final=float(C[-1]),physical_C_min=float(C.min()),physical_C_max=float(C.max()))
        final=complete_states[1][-1];factor=np.exp(-2j*np.pi*np.fft.fftfreq(n)*float(angle[-1])/360*n)
        final_u=final[:n] if spring else final[0];shifted_u=np.fft.ifft(np.fft.fft(u0)*factor).real
        primary_error=float(abs(final_u-shifted_u).max()/np.ptp(u0));peaks,regions=prior.periodic_counts(final_u)
        if spring:
            shifted=np.r_[shifted_u,np.fft.ifft(np.fft.fft(initial[n:-1])*factor).real,initial[-1]]
            full_error=float(abs(final[:-1]-shifted[:-1]).max())
        else:shifted=np.fft.ifft(np.fft.fft(initial,axis=1)*factor[None,:],axis=1).real;full_error=float(abs(final-shifted).max())
        scales=case["fixed_scales"] if spring else np.asarray(case["fixed_scales"])[:,None]
        fulltail=prior.tail_error(parts[1],scales)
        primarytail=float(np.ptp(parts[1]["u"][parts[1]["time"]>=12000],axis=0).max()/np.ptp(u0))
        pool_error=0.
        for states in complete_states:
            if spring:
                amount=states[:,-1]*(states[:,:n]+states[:,n:-1]).mean(axis=1);initialpool=Cb*(initial[:n]+initial[n:-1]).mean()
                pool_error=max(pool_error,float(abs(amount/initialpool-1).max()))
            else:
                model=base.ring.native_model(case["model"])
                for group in model.conserved_groups:
                    initialpool=Cb*initial[list(group)].sum()/n;amount=Cb*states[:,list(group)].sum(axis=(1,2))/n
                    pool_error=max(pool_error,float(abs(amount/initialpool-1).max()))
        assert pool_error<1e-7
        row.update(final_resolved_peaks=peaks,final_halfcontrast_regions=regions,
            final_primary_rotated_basal_relative_contrast_error=primary_error,final_full_rotated_basal_absolute_error=full_error,
            final_full_field_tail_excursion=fulltail,full_field_tail_gate=None if fulltail is None else bool(fulltail<=1e-5),
            full_field_tail_status="unavailable_from_retained_primary_only_history" if fulltail is None else ("passes_finite_tail_gate" if fulltail<=1e-5 else "fails_finite_tail_gate"),
            final_primary_tail_relative_contrast_excursion=primarytail,conserved_pool_max_relative_error=pool_error)
        if case.get("source_angle_arm_condition_id") is not None:
            previous=old[old.condition_id==case["source_angle_arm_condition_id"]].iloc[0]
            for phase in ("cue","recovery"):assert row["actual_"+phase+"_sha256"]==previous["actual_"+phase+"_sha256"]
            metrics=[k for k in previous.index if k in row and k not in case and isinstance(row[k],(float,int,np.number)) and not isinstance(row[k],bool)]
            replay=max((abs(row[k]-previous[k]) for k in metrics if not pd.isna(previous[k]) and row[k] is not None),default=0.)
            assert replay<1e-10;maximum_reuse_replay_error=max(maximum_reuse_replay_error,replay);shared_replays+=1
            assert row["full_field_tail_status"]==previous["full_field_tail_status"]
            geometry=["physical_distance_at_onset","fraction","requested_angle_deg","target_angle_in_P_coordinates_deg","circumference",
                      "full_physical_support","requested_peak","initial_sampled_cue_peak","initial_cue_integral","initial_angular_support_deg",
                      "original_core_physical_width","cue_support_core_intersection","overlap_fraction_cue_support","overlap_fraction_original_core",
                      "initial_cue_weighted_active","initial_cue_weighted_inactive"]
            geometry_error=max(abs(row[k]-previous[k]) for k in geometry)
            assert geometry_error<1e-10;maximum_shared_geometry_error=max(maximum_shared_geometry_error,geometry_error)
            assert bool(row["initial_overlap"])==bool(previous["initial_overlap"])
        rows.append(row)
        for threshold in (.9,.5,.1):
            for phase,keep in (("cue_on",t<=T),("release",t>=T),("entire_history",np.ones(len(t),bool))):
                tt=t[keep];ee=E[keep];mask=ee<=threshold;cross=np.flatnonzero(mask)
                meta=dict(condition_id=cid,model=case["model"],size_parameter=case["size_parameter"],size_label=case["size_label"],placement_label=case["placement_label"],
                    requested_angle_deg=case["requested_angle_deg"],cue_duration=T,initial_overlap=case["initial_overlap"],phase=phase,threshold=threshold,
                    phase_start=float(tt[0]),phase_horizon=float(tt[-1]),sampled_minimum=float(ee.min()),sampled_minimum_time=float(tt[ee.argmin()]),
                    rises_above_threshold_after_first_saved_entry=bool(np.any(~mask[int(cross[0]):])) if len(cross) else False,fraction_at_phase_end=float(ee[-1]))
                for event,bracket in (("first_saved_entry",rich.first_bracket(tt,mask)),("entry_retained_through_all_later_saved_phase_frames",rich.retained_bracket(tt,mask))):
                    entry=dict(**meta,event=event,**bracket,phase_local_lower=bracket["lower"]-tt[0],
                               phase_local_upper=None if bracket["upper"] is None else bracket["upper"]-tt[0])
                    events.append(entry)
                    if case.get("source_angle_arm_condition_id") is not None:
                        before=oldevent_lookup[case["source_angle_arm_condition_id"],threshold,phase,event]
                        assert entry["status"]==before["status"]
                        for key in ("lower","upper","phase_start","phase_horizon","sampled_minimum","sampled_minimum_time","fraction_at_phase_end"):
                            if entry[key] is None:assert pd.isna(before[key])
                            else:
                                error=abs(entry[key]-before[key]);assert error<1e-10
                                maximum_shared_event_error=max(maximum_shared_event_error,error)
                        assert entry["rises_above_threshold_after_first_saved_entry"]==before["rises_above_threshold_after_first_saved_entry"]
                        shared_events+=1
        for key,value in dict(time=t,phase=np.where(t<=T,0,1),M0=M0,M0_ratio=M0ratio,M1_real=M1.real,M1_imag=M1.imag,P_real=P.real,P_imag=P.imag,R=R,
            angle_raw_deg=angle,angle_masked_deg=np.where(R>=.1,angle,np.nan),q1=q1,q2=q2,fixed_region_transported_background_excess=E,
            fixed_region_physical_active_amount=A,C=C,prescribed_physical_offset=case["fraction"]*C).items():
            traces[f"c{cid}_{key}"]=value
            if case.get("source_angle_arm_condition_id") is not None:
                before=oldtrace[f'c{case["source_angle_arm_condition_id"]}_{key}']
                if key=="phase":np.testing.assert_array_equal(value,before)
                else:
                    np.testing.assert_allclose(value,before,rtol=0,atol=1e-10,equal_nan=True)
                    error=float(np.nanmax(abs(value-before)))
                    maximum_shared_trace_error=max(maximum_shared_trace_error,error)
        aq=(360*(q-case["preparation_axis_fraction"])+120)%360-120
        for j in range(n):profiles.append(dict(condition_id=cid,model=case["model"],size_parameter=case["size_parameter"],placement_label=case["placement_label"],cue_duration=T,
            material_position=q[j],normal_angle_from_preparation_deg=aq[j],final_active=final_u[j],rotated_basal_active=shifted_u[j]))
        audits.append(dict(condition_id=cid,pool_error=pool_error,exact_basal_restart=True,exact_phase_restart=True,actual_hashes_and_cue_profiles_match=True,
            primary_rotated_basal_error=primary_error,full_field_tail_status=row["full_field_tail_status"]))
    oldtrace.close()
    assert shared_replays==125 and len(context_cues)==48
    frame=pd.DataFrame(rows);csvwrite("stationary_arclength_summary.csv",frame)
    csvwrite("stationary_arclength_crossings.csv",pd.DataFrame(events));csvwrite("stationary_arclength_final_profiles.csv",pd.DataFrame(profiles))
    csvwrite("stationary_arclength_initial_profiles_and_cues.csv",pd.DataFrame(contexts));csvwrite("stationary_arclength_initial_context.csv",pd.read_csv(record(RUN/"results/initial_context.csv")))
    tracepath=dest("stationary_arclength_traces.npz");temporary=tracepath.with_name(tracepath.stem+".tmp.npz");np.savez_compressed(temporary,**traces);temporary.replace(tracepath)
    reconciliation=reconcile_arms(manifest,frame,old)
    generated=datetime.now(timezone.utc).isoformat();complete=int((frame.status=="complete").sum());allocated=sum(p.stat().st_blocks*512 for p in RUN.rglob("*") if p.is_file())
    validation=dict(generated_at=generated,passed_for_completed_subset=True,all_requested_complete=complete==240,requested=240,complete=complete,
        mutable_inventory_source=str(inventory_source),inventory_snapshot="supporting/data/source_condition_inventory.csv",inventory_snapshot_sha256=sha(inventory_snapshot),
        exact_reuses=240-len(manifest["jobs"]),new_native=len(manifest["jobs"]),initial_overlap_conditions=int(frame.initial_overlap.sum()),excluded_for_overlap=0,
        maximum_angle_arm_shared_scalar_replay_error=maximum_reuse_replay_error,maximum_angle_arm_shared_trace_replay_error=maximum_shared_trace_error,
        maximum_angle_arm_shared_geometry_replay_error=maximum_shared_geometry_error,maximum_angle_arm_shared_regional_event_replay_error=maximum_shared_event_error,
        shared_angle_arm_histories_replayed=shared_replays,shared_angle_arm_regional_events_replayed=shared_events,arm_reconciliation=reconciliation,
        maximum_saved_complete_state_pool_error=max((a["pool_error"] for a in audits),default=0.),
        all_frozen_source_and_input_hashes_match=True,all_completed_actual_archive_hashes_and_native_cue_profiles_match=True,records=audits,
        field_tail_counts=frame.full_field_tail_status.value_counts().to_dict(),threshold_records=len(events),raw_campaign_allocated_bytes=allocated,
        raw_campaign_large_manifest_required=bool(allocated>=2*1024**3))
    jsonwrite("validation.json",validation)
    for module in (producer,base,prior,rich):record(module.__file__)
    record(__file__);jsonwrite("source_sha256.json",SOURCE)
    write_readme(generated,complete)
    artifacts={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob("*") if p.is_file() and p.name!="latest_status.json" and not p.name.endswith((".tmp",".lock"))}
    producer_state="finished" if complete==240 else ("terminal_native_outputs_incomplete" if TERMINAL else ("waiting_for_native_jobs" if not (RUN/"launch.json").exists() else "native_campaign_in_progress"))
    jsonwrite("latest_status.json",dict(generated_at=generated,producer=producer_state,
        requested=240,complete=complete,pending_or_error=240-complete,status_counts=frame.status.value_counts().to_dict(),terminal_collection=TERMINAL,
        field_tail_counts=validation["field_tail_counts"],excluded_for_overlap=0,all_requested_complete=complete==240,
        compact_outputs_are_atomic=True,arm_reconciliation=reconciliation,artifact_sha256=artifacts))
    print(json.dumps({k:v for k,v in validation.items() if k!="records"},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--run",type=Path,default=RUN);parser.add_argument("--out",type=Path,default=OUT);parser.add_argument("--terminal",action="store_true");args=parser.parse_args()
    RUN=args.run.resolve();OUT=args.out.resolve();TERMINAL=args.terminal
    assert OUT!=base.TOPIC.parent/"standard_angles_sizes"
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/".collector.lock").open("w") as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
        collect()
