"""Authoritative direct-field observables for the immutable 240-condition grid."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import revision_standard_angles_sizes_20261019 as producer
import revision_standard_angles_duration_20261019 as base
import revision_standard_angles_report_20261019 as prior
import revision_stationary_characterization_20261019 as rich

RUN=producer.DEFAULT
OUT=producer.OUT
SOURCE={}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def record(path):
    path=Path(path);SOURCE[str(path)]=sha(path);return path
def dest(name,group="data"):
    parent=OUT/"supporting"/group;parent.mkdir(parents=True,exist_ok=True);return parent/name
def jsonwrite(name,value):
    target=dest(name,"metadata");temporary=target.with_suffix(target.suffix+".tmp");temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n");temporary.replace(target)
def csvwrite(name,frame):
    target=dest(name);temporary=target.with_suffix(target.suffix+".tmp");frame.to_csv(temporary,index=False);temporary.replace(target)


def collect():
    manifest=json.loads(record(RUN/"manifest.json").read_text());inventory=pd.read_csv(record(RUN/"results/condition_inventory.csv"))
    assert len(inventory)==240
    for rel,digest in manifest["source_sha256"].items():assert sha(RUN/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(RUN/rel)==digest
    rows=[];events=[];traces={};profiles=[];contexts=[];audits=[];basal_context={}
    old=pd.read_csv(record(base.TOPIC.parent/"stationary_characterization/supporting/data/stationary_whole_grid_summary.csv"))
    maximum_reuse_replay_error=0.
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
        if case.get("source_condition_id") is not None:
            previous=old[old.condition_id==case["source_condition_id"]].iloc[0]
            metrics=["M0_ratio_withdrawal","M0_ratio_final","R_withdrawal","R_final","fixed_region_excess_withdrawal","fixed_region_excess_final","R_sampled_minimum"]
            replay=max(abs(row[k]-previous[k]) for k in metrics);assert replay<1e-10;maximum_reuse_replay_error=max(maximum_reuse_replay_error,replay)
        rows.append(row)
        for threshold in (.9,.5,.1):
            for phase,keep in (("cue_on",t<=T),("release",t>=T),("entire_history",np.ones(len(t),bool))):
                tt=t[keep];ee=E[keep];mask=ee<=threshold;cross=np.flatnonzero(mask)
                meta=dict(condition_id=cid,model=case["model"],size_parameter=case["size_parameter"],size_label=case["size_label"],angle_label=case["angle_label"],
                    requested_angle_deg=case["requested_angle_deg"],cue_duration=T,initial_overlap=case["initial_overlap"],phase=phase,threshold=threshold,
                    phase_start=float(tt[0]),phase_horizon=float(tt[-1]),sampled_minimum=float(ee.min()),sampled_minimum_time=float(tt[ee.argmin()]),
                    rises_above_threshold_after_first_saved_entry=bool(np.any(~mask[int(cross[0]):])) if len(cross) else False,fraction_at_phase_end=float(ee[-1]))
                events.append(dict(**meta,event="first_saved_entry",**rich.first_bracket(tt,mask)))
                events.append(dict(**meta,event="entry_retained_through_all_later_saved_phase_frames",**rich.retained_bracket(tt,mask)))
        for key,value in dict(time=t,phase=np.where(t<=T,0,1),M0=M0,M0_ratio=M0ratio,M1_real=M1.real,M1_imag=M1.imag,P_real=P.real,P_imag=P.imag,R=R,
            angle_raw_deg=angle,angle_masked_deg=np.where(R>=.1,angle,np.nan),q1=q1,q2=q2,fixed_region_transported_background_excess=E,
            fixed_region_physical_active_amount=A,C=C,prescribed_physical_offset=case["fraction"]*C).items():traces[f"c{cid}_{key}"]=value
        aq=(360*(q-case["preparation_axis_fraction"])+120)%360-120
        for j in range(n):profiles.append(dict(condition_id=cid,model=case["model"],size_parameter=case["size_parameter"],angle_label=case["angle_label"],cue_duration=T,
            material_position=q[j],normal_angle_from_preparation_deg=aq[j],final_active=final_u[j],rotated_basal_active=shifted_u[j]))
        if T==(.1 if spring else 1.):
            for j in range(n):contexts.append(dict(model=case["model"],size_parameter=case["size_parameter"],size_label=case["size_label"],angle_label=case["angle_label"],
                material_position=q[j],normal_angle_from_preparation_deg=aq[j],basal_active=u0[j],basal_active_normalized_contrast=(u0[j]-u0.min())/np.ptp(u0),
                initial_cue=signal[j],initial_cue_unit_peak=signal[j]/case["requested_peak"],fixed_original_core_cell_material_weight=w[j]))
        audits.append(dict(condition_id=cid,pool_error=pool_error,exact_basal_restart=True,exact_phase_restart=True,actual_hashes_and_cue_profiles_match=True,
            primary_rotated_basal_error=primary_error,full_field_tail_status=row["full_field_tail_status"]))
    frame=pd.DataFrame(rows);csvwrite("stationary_sizes_summary.csv",frame)
    csvwrite("stationary_sizes_crossings.csv",pd.DataFrame(events));csvwrite("stationary_sizes_final_profiles.csv",pd.DataFrame(profiles))
    csvwrite("stationary_sizes_initial_profiles_and_cues.csv",pd.DataFrame(contexts));csvwrite("stationary_sizes_initial_context.csv",pd.read_csv(RUN/"results/initial_context.csv"))
    tracepath=dest("stationary_sizes_traces.npz");temporary=tracepath.with_name(tracepath.stem+".tmp.npz");np.savez_compressed(temporary,**traces);temporary.replace(tracepath)
    generated=datetime.now(timezone.utc).isoformat();complete=int((frame.status=="complete").sum());allocated=sum(p.stat().st_blocks*512 for p in RUN.rglob("*") if p.is_file())
    validation=dict(generated_at=generated,passed_for_completed_subset=True,all_requested_complete=complete==240,requested=240,complete=complete,
        exact_reuses=240-len(manifest["jobs"]),new_native=len(manifest["jobs"]),initial_overlap_conditions=int(frame.initial_overlap.sum()),excluded_for_overlap=0,
        maximum_original80_scalar_replay_error=maximum_reuse_replay_error,maximum_saved_complete_state_pool_error=max((a["pool_error"] for a in audits),default=0.),
        all_frozen_source_and_input_hashes_match=True,all_completed_actual_archive_hashes_and_native_cue_profiles_match=True,records=audits,
        field_tail_counts=frame.full_field_tail_status.value_counts().to_dict(),threshold_records=len(events),raw_campaign_allocated_bytes=allocated,
        raw_campaign_large_manifest_required=bool(allocated>=2*1024**3))
    jsonwrite("validation.json",validation)
    for module in (producer,base,prior,rich):record(module.__file__)
    record(__file__);jsonwrite("source_sha256.json",SOURCE)
    (OUT/"README.md").write_text(f"""# Three-size stationary-cue experiment

Generated: {generated}. Complete: {complete}/240 requested conditions. All 125 initial cue–core overlaps remain included; no new protocol or overlap exclusion is introduced.

The exact qualified chemical basals use C = 20/30/40; Spring uses resting C₀ = 1/2/3 with instantaneous physical C(t). Each basal is tested at prepared-axis cue-centre offsets π, π/2, π/4 and π/8 and five native durations, followed by 13,000 cue-free units. Continuous centres are never snapped. Physical supports and native peaks/clocks remain unchanged, so angular footprints vary with size. Native material totals and basal architecture also change with size; these are not isolated causal geometry contrasts.

- [Authoritative rich 240-condition summary](supporting/data/stationary_sizes_summary.csv)
- [Actual initial profiles and cue shapes](supporting/data/stationary_sizes_initial_profiles_and_cues.csv)
- [Input anchor, physical/angular width and overlap diagnostics](supporting/data/stationary_sizes_initial_context.csv)
- [All final primary profiles and shifted-basal comparisons](supporting/data/stationary_sizes_final_profiles.csv)
- [Regional threshold brackets and rebounds](supporting/data/stationary_sizes_crossings.csv)
- [Native replay and source validation](supporting/metadata/validation.json)

Raw M₁ components use fixed ring reference axes. P = M₁/M₁(0) rotates those axes once by the initial phase and rescales magnitude; P(0) = 1 + 0i. Physical M₀ is total active amount, R = |P| and q₁ = |M₁|/M₀. Both normalizations and instantaneous C are retained. Angles and angular changes are masked at endpoint R < 0.1; raw phases remain separate diagnostics. The fixed original material-region excess subtracts transported initial background amount; it measures regional activity rather than site ancestry. Saved-sample extrema and threshold brackets do not establish refined minimum cue duration.

The accompanying response renderer lives in `response/` and consumes these atomic compact tables. Earlier 80-condition checkpoints and all raw archives remain preserved.
""")
    print(json.dumps({k:v for k,v in validation.items() if k!="records"},indent=2))


if __name__=="__main__":collect()
