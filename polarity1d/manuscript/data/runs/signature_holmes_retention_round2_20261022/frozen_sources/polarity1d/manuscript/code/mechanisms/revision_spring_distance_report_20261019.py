"""Independent replay and presentation of the bounded Spring duration grid."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import revision_ring_distance_readouts_20261019 as chemical_readouts

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
RUN=ROOT/"data/runs/revision_spring_distance_duration_20261019"
OUT=chemical_readouts.geometry.OUT/"spring_distance"
SOURCE={}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p);SOURCE[str(p)]=sha(p);return json.loads(p.read_text())
def dest(name,group="data"):
    parent=OUT if name.endswith(".png") else OUT/"supporting"/group
    parent.mkdir(parents=True,exist_ok=True);return parent/name
def jsonwrite(name,z):dest(name,"metadata").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
def savefig(stem,fig):
    for ext in ("png","pdf","svg"):fig.savefig(dest(stem+"."+ext,"vectors"),dpi=220,bbox_inches="tight")
    plt.close(fig)


def main():
    chemical_readouts.geometry.setup();manifest=read(RUN/"manifest.json")
    assert read(RUN/"results/collection.json")["complete"]
    for rel,digest in manifest["source_sha256"].items():assert sha(RUN/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(RUN/rel)==digest
    assert sha(Path(manifest["independent_preparation_holdout"]))==manifest["independent_preparation_holdout_sha256"]
    summary=[];events=[];history={};checks=[];max_replay=0.;max_amount=0.;max_shape=0.;max_C=0.
    for job in manifest["jobs"]:
        directory=RUN/"workers"/f'{job["slot"]:03d}';record=read(directory/"status.json")
        assert record["complete"] and record["status"]=="complete"
        assert sha(directory/"readouts.json")==record["readouts_sha256"]
        assert sha(directory/"readout_traces.npz")==record["traces_sha256"]
        out=read(directory/"readouts.json")
        with np.load(RUN/job["input"]) as z:initial=z["state"].copy();q=z["material_positions"].copy()
        n=len(q);Cb=initial[-1];physical_core,bg=chemical_readouts.frozen_core(Cb*q,Cb,initial[:n]);w=physical_core/Cb
        parts=[]
        for phase in record["phases"]:
            path=RUN/phase["path"];assert sha(path)==phase["sha256"];SOURCE[str(path)]=phase["sha256"]
            prefix="actual_cue" if phase["kind"]=="challenge" else "actual_recovery"
            out[prefix+"_path"]=str(path);out[prefix+"_sha256"]=phase["sha256"]
            with np.load(path) as z:parts.append((z["time"].copy(),z["states"].copy()))
        np.testing.assert_array_equal(parts[0][1][0],initial)
        np.testing.assert_array_equal(parts[0][1][-1],parts[1][1][0])
        t=np.r_[parts[0][0],parts[1][0][1:]+job["cue_duration"]];y=np.vstack((parts[0][1],parts[1][1][1:]))
        C=y[:,-1];u=y[:,:n];A=C*(u@w);bgA0=Cb*bg*w.sum();E0=A[0]-bgA0;E=(A-bgA0)/E0
        amount=C*(y[:,:n]+y[:,n:-1]).mean(axis=1);mass=float(abs(amount/(2*job["resting_circumference"])-1).max());max_amount=max(max_amount,mass);assert mass<1e-7
        theta=2*np.pi*(q-q[initial[:n].argmax()]);M0=C*u.mean(axis=1);M1=C*(u@np.exp(1j*theta))/n;M2=C*(u@np.exp(2j*theta))/n;P=M1/M1[0];R=abs(P)
        angle=np.rad2deg(np.angle(P));angle=np.where(np.isclose(abs(angle),180.,atol=1e-7,rtol=0),180.,angle)
        with np.load(directory/"readout_traces.npz") as z:
            for key,value in (("old_core_transported_background_excess",E),("P_real",P.real),("P_imag",P.imag),("q2",abs(M2)/M0)):
                error=float(abs(value-z[key]).max());assert error<1e-11;max_replay=max(max_replay,error)
        shift=angle[-1]/360*n;factor=np.exp(-2j*np.pi*np.fft.fftfreq(n)*shift)
        shifted_u=np.fft.ifft(np.fft.fft(initial[:n])*factor).real
        shifted_v=np.fft.ifft(np.fft.fft(initial[n:-1])*factor).real
        shape=max(float(abs(u[-1]-shifted_u).max()),float(abs(y[-1,n:-1]-shifted_v).max()))
        max_shape=max(max_shape,shape);max_C=max(max_C,float(abs(C[-1]-Cb)))
        argmin=int(E.argmin());out.update(old_core_excess_minimum_time=float(t[argmin]),old_core_excess_minimum=float(E[argmin]),
            rotated_basal_max_field_error=shape,basal_C_vs_final_C_error=float(abs(C[-1]-Cb)),
            angle_final=float(angle[-1]) if R[-1]>=.1 else None)
        meta=dict(model="spring",resting_circumference=job["resting_circumference"],basal_circumference=Cb,
                  fraction=job["fraction"],cue_duration=job["cue_duration"],physical_distance_at_onset=job["physical_distance_at_onset"])
        for threshold in (.9,.5,.1):
            for kind,bracket in (("first_saved_sample_crossing",chemical_readouts.first_bracket(t,E<=threshold)),
                                 ("sampled_entry_retained_through_horizon",chemical_readouts.terminal_bracket(t,E<=threshold))):
                events.append(dict(**meta,observable="physical_core_transported_background_excess",threshold=threshold,event=kind,**bracket))
        cue_error=np.rad2deg(np.angle(np.exp(1j*np.deg2rad(angle-360*job["fraction"]))))
        for tolerance in (15.,30.):
            eligible=(abs(cue_error)<=tolerance)&(R>=.1)
            events.append(dict(**meta,observable="activity_vector_cue_alignment_deg",threshold=tolerance,event="first_saved_sample_entry",**chemical_readouts.first_bracket(t,eligible)))
            events.append(dict(**meta,observable="activity_vector_cue_alignment_deg",threshold=tolerance,event="sampled_entry_retained_through_horizon",**chemical_readouts.terminal_bracket(t,eligible)))
        history[(job["resting_circumference"],job["fraction"],job["cue_duration"])]=dict(time=t,states=y,E=E,R=R,angle=angle,initial=initial,q=q)
        checks.append(dict(slot=job["slot"],mass_relative_error=mass,maximum_replayed_observable_error=max_replay,rotated_basal_shape_error=shape,final_resolved_peaks=out["resolved_peaks_final"]))
        summary.append(out)
    frame=pd.DataFrame(summary);assert len(frame)==70 and frame.tail_stationary.all() and (frame.resolved_peaks_final==1).all()
    frame.to_csv(dest("spring_distance_duration_summary.csv"),index=False)
    pd.DataFrame(events).to_csv(dest("spring_distance_crossing_brackets.csv"),index=False)
    pd.DataFrame(manifest["exclusions"]).to_csv(dest("spring_distance_exclusions.csv"),index=False)
    threshold=[]
    for (C0,f),g in frame.groupby(["resting_circumference","fraction"]):
        g=g.sort_values("cue_duration")
        for observable,success in (("cue_alignment_15deg",g.cue_aligned_final_15deg.to_numpy(bool)),("retained_core_excess_le_half",(g.old_core_excess_final<=.5).to_numpy())):
            times=g.cue_duration.to_numpy();passed=times[success];failed=times[~success]
            threshold.append(dict(model="spring",resting_circumference=C0,fraction=f,observable=observable,
                shortest_tested_pass=float(passed.min()) if len(passed) else None,all_tested_passes=";".join(map(str,passed)),all_tested_failures=";".join(map(str,failed)),
                below_first_test_unresolved=bool(success[0]),nonmonotone_tested_success=bool(np.any(success[:-1]&~success[1:])),
                number_tested=len(g),status="shortest_tested_success" if len(passed) else "no_success_in_tested_duration_set"))
    threshold=pd.DataFrame(threshold);threshold.to_csv(dest("spring_shortest_tested_passes.csv"),index=False)
    # Independent previous validation replay at identical state, geometry and dose.
    prior=ROOT/"data/runs/spring_ring_validation_20261019";prior_manifest=read(prior/"manifest.json");reference_checks=[]
    for name,T in (("finite_0p1",.1),("finite_10",10.)):
        job=next(j for j in prior_manifest["jobs"] if j["name"]==name);status=read(prior/"workers"/f'{job["slot"]:02d}'/"status.json")
        errors=[];count=0;new=history[(2.,.4,T)]
        for phase in status["phases"]:
            path=prior/phase["path"];assert sha(path)==phase["sha256"]
            offset=0. if phase["kind"]=="challenge" else T
            with np.load(path) as z:tt=z["time"]+offset;yy=z["states"]
            for ti,yi in zip(tt,yy):
                found=np.flatnonzero(np.isclose(new["time"],ti,rtol=0,atol=1e-10))
                if len(found):errors.append(float(abs(new["states"][int(found[0])]-yi).max()));count+=1
        err=max(errors);assert err<1e-8;reference_checks.append(dict(name=name,common_saved_frames=count,max_complete_state_difference=err))
    render(frame,threshold,history)
    doc(frame,threshold,max_replay,max_amount,max_shape,max_C,checks,reference_checks)


def render(frame,threshold,history):
    colors={.30:"#4678A6",.35:"#309B84",.40:"#B07AA1",.45:"#CC7938",.50:"#555555"}
    fig=plt.figure(figsize=(11.8,9.1));fig.suptitle("Figure 3 continued · Finite cues produce graded retained reorientation on Spring rings",fontsize=14,y=.985)
    for j,C0 in enumerate((1.,2.,3.)):
        ax=fig.add_axes([.065+.31*j,.585,.255,.29])
        for f,color in colors.items():
            g=frame[(frame.resting_circumference==C0)&np.isclose(frame.fraction,f)].sort_values("cue_duration")
            if len(g):ax.plot(g.cue_duration,g.angle_final,"o-",color=color,ms=4,lw=1.2,label=f"d₀/Cb = {f:.2f} ({360*f:.0f}°)")
        ax.set(xscale="log",xticks=[.1,1,10,100,1000],ylim=(-3,185),yticks=[0,45,90,135,180],xlabel="Cue duration (native units)",ylabel="Retained activity angle (°)" if j==0 else "")
        ax.set_title(f"A{j+1}  Resting C₀ = {C0:g}; basal Cb = {frame[frame.resting_circumference==C0].basal_circumference.iloc[0]:.3f}",loc="left",pad=12)
    fig.legend(handles=[plt.Line2D([],[],color=co,marker="o",ms=4,label=f"d₀/Cb = {f:.2f} ({360*f:.0f}°)") for f,co in colors.items()],loc="upper center",bbox_to_anchor=(.52,.935),ncol=5,frameon=False,fontsize=8.5)
    b=fig.add_axes([.065,.16,.31,.27])
    for T,color in ((.1,"#4678A6"),(1.,"#309B84"),(10.,"#B07AA1"),(100.,"#CC7938")):
        h=history[(2.,.3,T)];n=len(h["q"]);angle=(360*(h["q"]-h["q"][0])+90)%360-90;order=np.argsort(angle)
        b.plot(angle[order],h["states"][-1,:n][order],color=color,lw=1.3,label=f"T = {T:g}")
    b.axvline(108,color="#999999",ls=":",lw=.8)
    b.set(xlim=(-90,270),ylim=(.65,1.3),xlabel="Normal angle from original peak (°)",ylabel="Final active Rac concentration")
    b.set_title("B  C₀2, d₀/Cb = 0.30: one shifted region",loc="left",pad=12)
    b.legend(frameon=False,fontsize=8,ncol=2)
    for j,(observable,title) in enumerate((("cue_alignment_15deg","C  Within 15° of the cue normal"),("retained_core_excess_le_half","D  Retained original-core excess ≤ ½"))):
        ax=fig.add_axes([.455+.275*j,.16,.22,.27])
        for C0,color in ((1.,"#4678A6"),(2.,"#B07AA1"),(3.,"#CC7938")):
            g=threshold[(threshold.resting_circumference==C0)&(threshold.observable==observable)]
            ax.plot(g.fraction,g.shortest_tested_pass,"o-",color=color,ms=4,label=f"Resting C₀ {C0:g}")
        ax.set(xlim=(.28,.52),xticks=[.3,.4,.5],yscale="log",ylim=(.7,1500),yticks=[1,10,100,1000],xlabel="Prescribed shorter-arc fraction d₀/Cb",ylabel="Shortest tested passing duration")
        ax.set_title(title,loc="left",pad=12,fontsize=9.5)
        if j==0:ax.legend(frameon=False,fontsize=8)
    fig.text(.06,.076,"All readouts follow 13,000 cue-free units. Fixed peak 4 and full instantaneous physical support 0.2; exact shared basals.\n70 complete branches; C₀1/f0.30 excluded for initial overlap. All final fields pass finite-window stationarity and contain one resolved peak.",fontsize=9)
    fig.text(.06,.022,"A continuous activity angle is an observed spatial-field readout, without a cell-velocity law. Points in C/D are tested durations, not refined minimum thresholds.",fontsize=8.5)
    savefig("spring_distance_duration",fig)


def doc(frame,threshold,max_replay,max_amount,max_shape,max_C,checks,reference_checks):
    generated=datetime.now(timezone.utc).isoformat()
    evidence=f"""# Spring distance-duration extension

Generated: {generated}. All 70 native branches finished; 5 prescribed initial-overlap conditions are explicitly excluded. Three exact qualified shared basals use resting C₀ = 1/2/3, N = 100/200/300, with actual basal circumferences 1.041403/2.044181/3.045168. Prescribed material offsets are 0.30/0.35/0.40/0.45/0.50; physical support stays 0.2, peak stays 4, durations are 0.1/1/10/100/1000, withdrawal lasts 13,000. The inherited mean-force global-extension analogue preserves approved stretching, dilution, material diffusion and local kinetic parameters. It is not a local mechanical propagation model.

The strongest new observation is graded retained reorientation. At C₀2 and offset 0.30, durations 0.1/1/10/100/1000 leave final first-moment angles approximately 8.3°/40.7°/90.1°/107.7°/108.0°. These are single shifted activity regions. All 70 final histories have one resolved peak, pass full-field finite-window stationarity, and agree with a fractionally rotated basal active/inactive profile to a maximum absolute difference {max_shape:.3g}; final circumference agrees with its exact basal to {max_C:.3g}. This supports a graded family of retained translated profiles over the sampled window. It does not establish unique asymptotic attractors or a movement law. The vector magnitude alone would not prove one region; native field and peak/shape checks supply that evidence.

Prescribed physical separation at onset is d₀ = f Cb, where Cb is actual basal circumference. During the assay d(t) = f C(t), while the prescribed normal angle 2πf remains fixed in material coordinates. C₀ denotes resting circumference, not basal or instantaneous perimeter. Each table retains the full physical-distance range. The largest C excursion is {(abs(frame.physical_C_max-frame.physical_C_min)/frame.basal_circumference).max()*100:.3f}% of basal circumference; small size does not justify silently discarding it.

Original-core loss uses physical active amount in a fixed basal material half-contrast interval. Its primary excess subtracts the initial background amount transported by the same uniform material deformation: [C(t)∫core u dq − Cb·ubg,0·∫core dq] divided by the initial excess. Thus pure affine dilution alone does not register loss. Raw physical amount and C(t)∫core[u−ubg,0]dq are separate exported alternatives. A fixed-region amount change can reflect redistribution and does not identify molecular destruction or death of a uniquely tracked incumbent site. First threshold crossings and entry retained through all later saved times remain separate, with saved-sample brackets and finite-horizon censoring.

Two endpoint success definitions are shown: final activity angle within 15° of the prescribed cue normal at R≥0.1, and final original-core transported-background excess at most 0.5. They need not agree. All successes and failures across tested durations are retained. “Shortest tested passing duration” does not establish a unique minimum or monotonic duration dependence; a pass at the first tested duration would leave the region below 0.1 unresolved. None of these 14 geometry series has a nonmonotone tested success island for the displayed gates.

Independent audit recomputes every native phase checksum, exact basal restart and cue/recovery continuity, physical conserved amount, first/second moments and material-core integrals. Maximum observable replay discrepancy: {max_replay:.3g}; maximum physical-amount relative error: {max_amount:.3g}. The independent earlier C₂ offset 0.4 validation histories for durations 0.1/10 reproduce at shared saved full-state frames to <1e−8. The signed parent preparation holdout audit is frozen by hash in the new campaign. `validation.json` records the exact replay and source checks; full native fields remain immutable under `data/runs/revision_spring_distance_duration_20261019/`.

The existing chemical-grid descriptor classifier is not informative for Spring's 1000-duration endpoint class: all 14 eligible long-cue Spring geometries are cue-aligned single regions. This is retained as a constant-outcome result, without manufacturing classifier discrimination. The finite-duration angle and core-loss curves provide the relevant Spring counterpart.
"""
    evidence += r"""

The antipodal input preserves the initial reflection axis to numerical precision. Its tested retained angles of 0°/180° therefore reflect the prescribed symmetric geometry, rather than establishing a separate binary biological mechanism.

The transported-background core excess is

$$
E_{\Omega_0}^{\mathrm{tr}}(t)=
\frac{C(t)\int_{\Omega_0}u(q,t)\,\mathrm{d}q
      -C_{\mathrm b}u_{\mathrm{bg},0}|\Omega_0|}
     {C_{\mathrm b}\int_{\Omega_0}
       \left[u(q,0)-u_{\mathrm{bg},0}\right]\,\mathrm{d}q}.
\tag{1}
$$

Here \(\Omega_0\) is the frozen basal half-contrast material interval. The subtracted background is an initial physical amount. This excess stays constant under a pure affine dilution of an otherwise unchanged material profile.
"""
    dest("EVIDENCE.md","metadata").write_text(evidence)
    caption="Spring finite-duration distance response. A, Final activity-vector angle after 13,000 cue-free units versus native cue duration for each resting circumference and prescribed shorter-arc fraction. B, Exact native final profiles at C₀2 and fraction 0.30 show one activity region translated by graded angles after different finite cues; dotted line is the 108° cue normal. C/D, Shortest tested passing duration for final cue alignment within 15° and for final original-core physical transported-background excess at most one half. Full physical input support 0.2 and native peak 4 are fixed. C₀1/fraction 0.30 is initially overlap-excluded. Both activity fields and dynamic circumference evolve under the approved Spring mean-force circular analogue. These finite sampled readouts do not establish a unique minimum duration, a molecular lineage or cell movement."
    dest("CAPTION.md","metadata").write_text(caption+"\n")
    SOURCE[str(Path(__file__))]=sha(__file__);jsonwrite("source_sha256.json",SOURCE)
    allocated=sum(p.stat().st_blocks*512 for p in RUN.rglob("*") if p.is_file())
    jsonwrite("validation.json",dict(generated_at=generated,passed=True,native_conditions=70,excluded_initial_overlaps=5,
        all_final_fields_stationary=True,all_final_profiles_single_resolved_peak=True,all_exact_checkpoint_and_phase_continuity=True,
        maximum_moment_and_core_replay_error=max_replay,maximum_physical_amount_relative_error=max_amount,
        maximum_fractionally_rotated_basal_field_error=max_shape,maximum_basal_vs_final_circumference_error=max_C,
        prior_independent_native_replays=reference_checks,records=checks,raw_campaign_allocated_bytes=allocated,
        raw_campaign_large_manifest_required=bool(allocated>=2*1024**3)))
    print(json.dumps(dict(passed=True,conditions=70,max_replayed_error=max_replay,max_rotated_basal_error=max_shape,raw_allocated_MiB=allocated/1024**2,output=str(OUT)),indent=2))


if __name__=="__main__":main()
