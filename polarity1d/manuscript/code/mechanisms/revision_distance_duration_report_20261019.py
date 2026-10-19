"""Audit the chemical extension and compare finite duration outcomes on four rings."""
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
import revision_ring_distance_readouts_20261019 as retained
import ring_distance_prototype_20261017 as ring

ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
RUN=ROOT/"data/runs/revision_chemical_distance_duration_20261019"
SPRING=ROOT/"data/runs/revision_spring_distance_duration_20261019"
OUT=retained.OUT
SOURCE={}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p);SOURCE[str(p)]=sha(p);return json.loads(p.read_text())
def dest(name,group="data"):
    parent=OUT if name.endswith(".png") else OUT/"supporting"/group
    parent.mkdir(parents=True,exist_ok=True);return parent/name
def jsonwrite(name,z):dest(name,"metadata").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
def savefig(name,fig):
    for ext in ("png","pdf","svg"):fig.savefig(dest(name+"."+ext,"vectors"),dpi=220,bbox_inches="tight")
    plt.close(fig)


def audit_chemical(manifest):
    checks=[];max_prefix=0.;max_replay=0.;max_pool=0.
    for job in manifest["jobs"]:
        directory=RUN/"workers"/f'{job["slot"]:03d}';status=read(directory/"status.json")
        assert status["status"]=="complete" and status["complete"]
        assert sha(directory/"readouts.json")==status["readouts_sha256"]
        assert sha(directory/"readout_traces.npz")==status["traces_sha256"]
        out=read(directory/"readouts.json");parts=[]
        with np.load(RUN/job["input"]) as z:initial=z["state"].copy()
        for phase in status["phases"]:
            path=RUN/phase["path"];assert sha(path)==phase["sha256"];SOURCE[str(path)]=phase["sha256"]
            with np.load(path) as z:parts.append((z["time"].copy(),z["states"].copy(),z["x"].copy()))
        np.testing.assert_array_equal(parts[0][1][0],initial);np.testing.assert_array_equal(parts[0][1][-1],parts[1][1][0])
        t=np.r_[parts[0][0],parts[1][0][1:]+job["duration"]];y=np.concatenate((parts[0][1],parts[1][1][1:]));x=parts[0][2]
        model=ring.native_model(job["model"]);dx=40/len(x);pool_errors=[]
        for group in model.conserved_groups:
            amount=dx*y[:,list(group),:].sum(axis=(1,2));pool_errors.append(float(abs(amount/amount[0]-1).max()))
        pool_error=max(pool_errors);assert pool_error<1e-7;max_pool=max(max_pool,pool_error)
        w,bg=retained.frozen_core(x,40.,initial[0]);A=y[:,0]@w;E=(A-bg*w.sum())/(A[0]-bg*w.sum())
        theta=2*np.pi*(x-x[initial[0].argmax()])/40.;M1=dx*(y[:,0]@np.exp(1j*theta));P=M1/M1[0]
        with np.load(directory/"readout_traces.npz") as z:
            error=max(float(abs(E-z["old_core_excess"]).max()),float(abs(P.real-z["P_real"]).max()),float(abs(P.imag-z["P_imag"]).max()))
            assert error<1e-10;max_replay=max(max_replay,error)
        assert sha(Path(job["source_cue"]))==job["source_cue_sha256"]
        with np.load(job["source_cue"]) as z:
            tt=z["time"];pp=z["states"][:,0] if "states" in z else z["primary"]
            differences=[]
            for ti,pi in zip(parts[0][0],parts[0][1][:,0]):
                match=np.flatnonzero(np.isclose(tt,ti,rtol=0,atol=1e-10))
                if len(match):differences.append(float(abs(pi-pp[int(match[0])]).max()))
            prefix=max(differences)/max(1.,float(pp.max()));assert prefix<1e-7;max_prefix=max(max_prefix,prefix)
        checks.append(dict(model=job["model"],fraction=job["fraction"],duration=job["duration"],
            phase_hashes_match=True,exact_native_restart=True,relative_pool_error=pool_error,
            maximum_observable_replay_error=error,source_1000_prefix_relative_error=prefix,
            final_field_tail_stationary=bool(out["tail_stationary"])))
    return dict(records=checks,maximum_pool_error=max_pool,maximum_observable_replay_error=max_replay,maximum_retained_prefix_scaled_error=max_prefix)


def merge():
    manifest=read(RUN/"manifest.json");assert read(RUN/"results/collection.json")["complete"]
    for rel,digest in manifest["source_sha256"].items():assert sha(RUN/"frozen_sources"/rel)==digest
    for rel,digest in manifest["input_sha256"].items():assert sha(RUN/rel)==digest
    audit=audit_chemical(manifest)
    new=pd.read_csv(RUN/"results/distance_duration_readouts.csv");new["source_kind"]="new_chemical_duration_branch"
    old=pd.read_csv(OUT/"supporting/data/distance_response_summary.csv")
    epochs=pd.read_csv(OUT/"supporting/data/distance_response_readouts.csv")
    references=[]
    for job in manifest["reused_references"]:
        r=old[(old.model==job["model"])&(old.circumference==40)&np.isclose(old.fraction,job["fraction"])].iloc[0]
        off=epochs[(epochs.model==job["model"])&(epochs.circumference==40)&np.isclose(epochs.fraction,job["fraction"])&(epochs.time==1000)].iloc[0]
        angle=r.final_angle_deg;error=r.final_cue_angle_error_deg
        references.append(dict(**r.to_dict(),status="complete",source_kind="exact_retained_1000_reference",
            actual_cue_path=job["source_cue"],actual_cue_sha256=job["source_cue_sha256"],
            actual_recovery_path=job["source_recovery"],actual_recovery_sha256=job["source_recovery_sha256"],
            R_final=r.final_R,R_at_withdrawal=off.R,angle_final=angle,angle_at_withdrawal=off.angle_deg,
            cue_angle_error_final=error,cue_aligned_final_15deg=bool(pd.notna(error) and abs(error)<=15 and r.final_R>=.1),
            tail_stationary=None,field_tail_status="full_field_tail_not_remeasured_from_compact_reference",
            cue_on_loss_0p5_status=r.half_loss_status,cue_on_loss_0p5_lower=r.half_loss_lower,cue_on_loss_0p5_upper=r.half_loss_upper))
    chemical=pd.concat((new,pd.DataFrame(references)),ignore_index=True)
    chemical["field_tail_status"]=chemical.apply(lambda r:r.field_tail_status if pd.notna(r.field_tail_status) else ("stationary_at_declared_finite_gate" if r.tail_stationary else "fails_declared_field_tail_gate"),axis=1)
    chemical.to_csv(dest("chemical_distance_duration_with_reference.csv"),index=False)
    spring=pd.read_csv(retained.geometry.OUT/"spring_distance/supporting/data/spring_distance_duration_summary.csv")
    spring["source_kind"]="new_Spring_duration_branch";spring["model"]="spring";spring["circumference"]=spring.resting_circumference
    spring["field_tail_status"]="stationary_at_declared_finite_gate"
    spring["physical_distance"]=spring.physical_distance_at_onset;spring["angular_separation_deg"]=360*spring.fraction
    combined=pd.concat((chemical,spring),ignore_index=True)
    combined["size_parameter"]=np.where(combined.model=="spring","resting circumference C0","fixed circumference C")
    exclusions=pd.DataFrame(read(SPRING/"manifest.json")["exclusions"]);exclusions["model"]="spring";exclusions["source_kind"]="initial_overlap_exclusion"
    exclusions["circumference"]=exclusions.resting_circumference
    exclusions["size_parameter"]="resting circumference C0"
    all_requested=pd.concat((combined,exclusions),ignore_index=True)
    all_requested.to_csv(dest("distance_duration_all_models.csv"),index=False)
    assert len(chemical)==75 and len(combined)==145 and len(all_requested)==150
    minimum=[]
    for (model,C,f),g in combined.groupby(["model","circumference","fraction"]):
        g=g.sort_values("cue_duration");times=g.cue_duration.to_numpy();R=g.R_final.to_numpy();error=g.cue_angle_error_final.to_numpy();E=g.old_core_excess_final.to_numpy()
        specs=[("cue_alignment_deg",tol,(abs(error)<=tol)&(R>=.1)) for tol in (10.,15.,20.)]
        specs += [("retained_original_core_excess_fraction",threshold,E<=threshold) for threshold in (.4,.5,.6)]
        for observable,threshold,success in specs:
            passed=times[success];failed=times[~success];first=float(passed.min()) if len(passed) else None
            prior=failed[failed<first] if first is not None else np.array([])
            minimum.append(dict(model=model,circumference=C,fraction=f,observable=observable,threshold=threshold,
                shortest_tested_pass=first,adjacent_tested_failure=float(prior.max()) if len(prior) else None,
                status="shortest_tested_success" if len(passed) else "no_success_in_tested_duration_set",
                first_positive_duration_already_passes=bool(success[0]),nonmonotone_tested_success=bool(np.any(success[:-1]&~success[1:])),
                all_tested_passes=";".join(map(str,passed)),all_tested_failures=";".join(map(str,failed)),number_tested=len(g)))
    minimum=pd.DataFrame(minimum);minimum.to_csv(dest("distance_duration_shortest_tested_sensitivity.csv"),index=False)
    render(combined,minimum)
    docs(combined,minimum,audit)


def render(frame,minimum):
    selected=[("wave_pinning",40,.2,.5),("goryachev",40,.1,.5),("holmes_model3",40,.25,.5),("spring",2,.3,.5)]
    colors=("#4678A6","#CB7737")
    fig=plt.figure(figsize=(11.8,10.5));fig.suptitle("Figure 3 continued · Distance and finite cue duration alter retained activity orientation",fontsize=14,y=.985)
    for i,(model,C,near,far) in enumerate(selected):
        left=.075+.48*(i%2);bottom=.675-.27*(i//2)
        ax=fig.add_axes([left,bottom,.37,.19])
        for f,color,label in ((near,colors[0],"Near"),(far,colors[1],"Far")):
            g=frame[(frame.model==model)&(frame.circumference==C)&np.isclose(frame.fraction,f)].sort_values("cue_duration")
            coordinate="d₀/Cb" if model=="spring" else "d/C"
            ax.plot(g.cue_duration,g.angle_final,"o-",color=color,ms=4,lw=1.2,label=f"{label}: {coordinate} = {f:.2f}, target {360*f:.0f}°")
            failed=g[g.field_tail_status=="fails_declared_field_tail_gate"]
            ax.scatter(failed.cue_duration,failed.angle_final,s=31,facecolor="white",edgecolor=color,linewidth=1.3,zorder=5)
            ax.axhline(360*f,color=color,ls=":",lw=.7,alpha=.7)
            if model=="holmes_model3" and label=="Far":ax.text(.03,.72,"Far T = 1000: R = 0.040; angle masked",transform=ax.transAxes,fontsize=8,color=color)
        ax.set(xscale="log",xlim=(.07 if model=="spring" else .7,1500),ylim=(-7,191),yticks=[0,45,90,135,180],xlabel="Cue duration (native units)",ylabel="Final activity angle (°)")
        name="Spring" if model=="spring" else retained.geometry.NAMES[model]
        ax.set_title(f"A{i+1}  {name} · "+("resting C₀ = 2; basal Cb = 2.044" if model=="spring" else "fixed C = 40"),loc="left",pad=8)
        ax.legend(frameon=True,facecolor="white",edgecolor="none",framealpha=.95,fontsize=8,loc="upper left",bbox_to_anchor=(.03,1.005))
    b=fig.add_axes([.075,.09,.37,.205]);manifest=read(SPRING/"manifest.json")
    for T,color in ((.1,"#4678A6"),(1.,"#309B84"),(10.,"#B07AA1"),(100.,"#CC7938")):
        job=next(j for j in manifest["jobs"] if j["resting_circumference"]==2. and np.isclose(j["fraction"],.3) and j["cue_duration"]==T)
        directory=SPRING/"workers"/f'{job["slot"]:03d}'
        with np.load(directory/"recovery.npz") as z:y=z["states"][-1];q=z["material_positions"];n=len(q)
        angle=(360*(q-q[0])+90)%360-90;order=np.argsort(angle)
        b.plot(angle[order],y[:n][order],color=color,lw=1.3,label=f"T = {T:g}")
    b.axvline(108,color="#888888",ls=":",lw=.7);b.set(xlim=(-90,270),ylim=(.65,1.3),xlabel="Normal angle from original peak (°)",ylabel="Final active Rac concentration")
    b.set_title("B  Spring near cue: one translated activity region",loc="left",pad=10);b.legend(frameon=False,fontsize=8,ncol=2)
    ax=fig.add_axes([.565,.073,.395,.23]);ax.axis("off");cells=[]
    for model,C,near,far in selected:
        for f,label in ((near,"near"),(far,"far")):
            g=minimum[(minimum.model==model)&(minimum.circumference==C)&np.isclose(minimum.fraction,f)]
            a=g[(g.observable=="cue_alignment_deg")&(g.threshold==15)].iloc[0]
            e=g[(g.observable=="retained_original_core_excess_fraction")&np.isclose(g.threshold,.5)].iloc[0]
            long=frame[(frame.model==model)&(frame.circumference==C)&np.isclose(frame.fraction,f)&(frame.cue_duration==1000)].iloc[0]
            short=lambda v:"—" if pd.isna(v) else f"{v:g}"
            name="Spring" if model=="spring" else retained.geometry.NAMES[model]
            cells.append([name+" "+label,short(a.shortest_tested_pass),short(e.shortest_tested_pass),f"{long.R_final:.3f}"])
    table=ax.table(cellText=cells,colLabels=["Same cases", "Cue ±15°", "Core ≤½", "R at T1000"],cellLoc="center",loc="center",colWidths=[.40,.20,.20,.20])
    table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.25)
    for (row,col),cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd");cell.set_linewidth(.5)
        if row==0:cell.set_facecolor("#eeeeee")
    ax.set_title("C  Shortest tested passes depend on the endpoint",loc="left",pad=12,fontsize=9.5)
    fig.text(.065,.925,"All points follow 13,000 cue-free units. Curves connect tested durations; dotted lines mark prescribed cue normals.",fontsize=9.5)
    fig.text(.065,.903,"Open circles: failed separate field-tail gate. The full-field tail was not remeasured for reused chemical T = 1000 references.",fontsize=8.5)
    fig.text(.065,.026,"Native clocks are not calibrated across models. —: no pass in the tested set. Core loss is a fixed-region activity readout, not site ancestry.\nChemical full support 3, native peaks 0.1/6.4/1; Spring full instantaneous physical support 0.2, peak 4. Full maps, tail flags and threshold sensitivities remain linked.",fontsize=8.5)
    savefig("revision_figure3_distance_duration",fig)


def docs(frame,minimum,audit):
    generated=datetime.now(timezone.utc).isoformat()
    evidence=f"""# Four-model finite-duration distance response

Generated: {generated}. The targeted extension contains 60 fresh chemical branches and 15 exact retained 1000-duration references, plus 70 fresh Spring branches and 5 excluded initial-overlap conditions. The combined requested table has 150 rows. Chemical tests are fixed C = 40, T = 1/10/100/300/1000; Spring tests use resting C₀1/2/3 and T = 0.1/1/10/100/1000. These native clocks are not a cross-model rate calibration. The original complete 652-history grid retains systematic size-distance evidence at cue duration 1000.

Selected near/far placements are declared in the second main sheet: Wave-pinning f = 0.20/0.50, Goryachev f = 0.10/0.50, Holmes f = 0.25/0.50 at C = 40; Spring f = 0.30/0.50 at resting C₀ = 2. Full five-placement duration grids remain in `distance_duration_all_models.csv` and `chemical_distance_duration_with_reference.csv`. The selected curves show distinct response routes, not a universal ordering by physical distance.

Wave-pinning and Spring can retain graded intermediate angles after finite near cues. Goryachev's cue-associated retained response needs a longer tested cue at the farther placement, consistent with the independently observed distance-dependent cue-on core-depletion times. Holmes can shift/merge near sites or retain distinct far regions while the global activity vector loses magnitude; it has no within-15° cue-normal endpoint in this bounded duration grid. Region activity, vector orientation and magnitude, and native peak organization are distinct readouts.

Shortest tested passes are exported separately for angle tolerances 10°/15°/20° at R≥0.1 and retained original-core excess fractions 0.4/0.5/0.6. Every tested successful island and failure remains explicit. At Spring C₀ = 2/f = 0.30, T = 10 leaves angle 90.12° for target 108°: it fails 15° but passes 20°, so the shortest tested angle pass changes 100→10 when the endpoint changes. No finite curve supplies a unique minimum pulse duration, monotonic success law, or exclusion of untested short successful islands. First saved crossing and entry retained through all later saved frames remain distinct statistics.

All 60 fresh chemical branches reproduce the exact initial native state and cue/recovery restart. Independently recomputed conserved-pool relative error is at most {audit['maximum_pool_error']:.3g}; core/vector replay error is at most {audit['maximum_observable_replay_error']:.3g}. Their cue-on primary profiles reproduce common saved times in the retained 1000 histories to maximum scaled discrepancy {audit['maximum_retained_prefix_scaled_error']:.3g}. Three far Holmes T = 300 branches fail the separate full-field final-tail gate; all completed finite readouts remain visible with that flag. For the compact retained 1000 references, no full-field tail qualification is invented where it was not remeasured. Stationarity is distinct from an activity-angle success.

Spring physical geometry is audited separately in `../spring_distance/supporting/metadata/validation.json`: physical amount includes C(t), and the fixed-material-core excess subtracts the transported initial background amount so pure affine dilution is not mistaken for loss. The first moment used here is the physical M₁ divided by initial M₁, not the API's M₁/M₀ normalization; both R and q₁ remain exported. Full native profile and fractional-rotation checks support the single translated Spring regions.
"""
    evidence += "\n## Shared Spring preparations and function assays\n\nThe [exact shared preparations](../../../../../../../../data/runs/revision_ring_spring_20261019/results/preparation_summary.json) and their independent holdout audit provide the common initial states for deterministic and stochastic experiments. The [basal function report](../../../../spring_function/README.md) links organization, native inactive-weighted gain and the common active-loss recovery assay at all three Spring sizes. The [challenged t10 function report](../../../../spring_challenged_function/README.md) links all 14 eligible exact challenge states and their regional recovery/gain measurements. Ten challenged regions lack a resolved old-associated peak, so those results remain explicitly fixed-region responses. They are not a partition of replacement time or an identified-site stability test. The [interactive distance-response explorer](../../../../supporting/viewers/distance_response.html) exposes the full curves and declared endpoint choices.\n"
    evidence += "\nThe [Spring descriptor coverage record](../../../spring_distance/supporting/metadata/descriptor_coverage.json) closes the standard-long-cue counterpart: all 14 eligible geometries have one cue-associated peak, so categorical descriptor ranking is uninformative. Continuous angle/core measurements remain available. Historical exact chemical interval/ring maps do not imply a nonuniform Spring interval-reflection identity under the different endpoint-force versus mean-force closures.\n"
    dest("DURATION_EVIDENCE.md","metadata").write_text(evidence)
    caption="Distance-duration response on four supported ring models. A1–A4, Final activity angle after 13,000 cue-free units for selected near/far prescribed cue centres; dotted lines mark the corresponding cue-normal angles. Exact qualified basals, native kinetics and fixed physical support are held within each model. Chemical sizes are C = 40; Spring uses resting C₀ = 2 and actual basal Cb = 2.044. Tested chemical durations are 1/10/100/300/1000, Spring 0.1/1/10/100/1000. At far Holmes T = 1000, R=0.040 masks direction rather than defining a zero angle. B, Native Spring near-cue final profiles show one full-strength region retained at graded translated positions, validated against a rotated exact basal profile. C, Shortest tested passing durations for two explicitly different endpoints: final angle within 15° of the cue normal at R ≥ 0.1, and final original-core excess no more than one half. A dash means no success in the tested duration set. Near Wave-pinning passes the angle gate at 1000 while retaining more than half the original-core excess, illustrating endpoint dependence. These are finite-window observations, not a unique minimum-duration search or migration prediction."
    caption += " Near denotes the selected nearest eligible prescribed placement at the displayed size; it is not a cross-model matched physical distance. Open symbols mark measured failure of the separate full-field final-tail criterion; compact retained chemical T = 1000 references have no invented tail flag. Lines guide the eye between tested points, with intermediate durations unmeasured. Antipodal inputs preserve the initial reflection axis, so 0°/180° far-cue outcomes do not establish a separate binary biological mechanism."
    caption += " The Holmes cue-associated single peak on the first sheet is at C = 30; these duration series use C = 40."
    dest("DURATION_CAPTION.md","metadata").write_text(caption+"\n")
    content=dict(assets=[dict(path=str((OUT/"revision_figure3_distance_duration.png").relative_to(REPO)),caption="Compact four-model duration-response continuation with native Spring profile evidence.",page_size=[11.8,11.7],page_title="Figure 3 continued. Distance and finite cue duration",minimum_printed_width_inches=10.)],legend=[caption],
        narrative=["Finite near cues can leave one full-strength activity region at a graded intermediate angular position; retained angle need not be a binary old/cue outcome.","Distance changes regional depletion, retained angle and vector magnitude differently across native models. The shortest tested successful duration depends on the explicitly chosen endpoint."],
        sources=[dict(label="Four-model duration evidence and limits",path=str(dest("DURATION_EVIDENCE.md","metadata").relative_to(REPO))),dict(label="All150requested duration conditions",path=str(dest("distance_duration_all_models.csv").relative_to(REPO))),dict(label="Endpoint sensitivity and all tested success/failure sets",path=str(dest("distance_duration_shortest_tested_sensitivity.csv").relative_to(REPO))),dict(label="Spring native field and geometry audit",path=str((retained.geometry.OUT/"spring_distance/supporting/metadata/validation.json").relative_to(REPO)))])
    content["supplementary"]=[
        dict(title="Complete fixed-duration distance response in all 652 native chemical rings",path=str((OUT/"distance_response_all_sizes.png").relative_to(REPO))),
        dict(title="All three Spring sizes and five finite-duration placements",path=str((retained.geometry.OUT/"spring_distance/spring_distance_duration.png").relative_to(REPO))),
        dict(title="All 150 requested duration conditions and actual native archive provenance",path=str((OUT/"supporting/data/distance_duration_all_models.csv").relative_to(REPO))),
        dict(title="Chemical 75-condition duration grid including exact retained references",path=str((OUT/"supporting/data/chemical_distance_duration_with_reference.csv").relative_to(REPO))),
        dict(title="Endpoint threshold sensitivity and complete tested success/failure sets",path=str((OUT/"supporting/data/distance_duration_shortest_tested_sensitivity.csv").relative_to(REPO))),
        dict(title="All original-core and activity-alignment crossing brackets",path=str((OUT/"supporting/data/distance_crossing_brackets.csv").relative_to(REPO))),
        dict(title="Spring transported-background core and alignment crossing brackets",path=str((retained.geometry.OUT/"spring_distance/supporting/data/spring_distance_crossing_brackets.csv").relative_to(REPO))),
        dict(title="Spring basal organization, gain and common active-loss recovery",path=str((retained.geometry.OUT.parent/"spring_function/README.md").relative_to(REPO))),
        dict(title="Spring challenged-state conditional regional function measurements",path=str((retained.geometry.OUT.parent/"spring_challenged_function/README.md").relative_to(REPO))),
        dict(title="Interactive distance-response explorer",path=str((retained.geometry.OUT.parent/"supporting/viewers/distance_response.html").relative_to(REPO)))]
    content["sources"].extend([
        dict(label="Exact shared Spring preparations and qualification",path=str((ROOT/"data/runs/revision_ring_spring_20261019/results/preparation_summary.json").relative_to(REPO))),
        dict(label="Independent Spring preparation windows and holdout",path=str((ROOT/"data/runs/revision_ring_spring_20261019/results/root_independent_preparation_holdout.json").relative_to(REPO))),
        dict(label="Spring basal function verification",path=str((retained.geometry.OUT.parent/"spring_function/supporting/metadata/verification.json").relative_to(REPO))),
        dict(label="Spring challenged-state function verification",path=str((retained.geometry.OUT.parent/"spring_challenged_function/supporting/metadata/verification.json").relative_to(REPO))),
        dict(label="Declared Spring ring model and approved geometric carryover",path=str((ROOT/"SPRING_RING_STATUS_20261019.md").relative_to(REPO)))])
    content["sources"].append(dict(label="Spring categorical invariance, continuous descriptors and interval-correspondence scope",path=str((retained.geometry.OUT/"spring_distance/supporting/metadata/descriptor_coverage.json").relative_to(REPO))))
    jsonwrite("content_continued.json",content)
    SOURCE[str(Path(__file__))]=sha(__file__);jsonwrite("duration_source_sha256.json",SOURCE)
    allocated=sum(p.stat().st_blocks*512 for p in RUN.rglob("*") if p.is_file())
    jsonwrite("duration_validation.json",dict(generated_at=generated,passed=True,native_chemical_branches=60,retained_chemical_references=15,
        native_Spring_branches=70,excluded_Spring_conditions=5,combined_requested_conditions=150,
        raw_chemical_allocated_bytes=allocated,large_manifest_required=bool(allocated>=2*1024**3),
        source_sha256={str(RUN/"manifest.json"):sha(RUN/"manifest.json")},**audit))
    print(json.dumps(dict(passed=True,chemical_branches=60,combined_requested=150,raw_chemical_MiB=allocated/1024**2,max_replay=audit["maximum_observable_replay_error"],output=str(OUT)),indent=2))


if __name__=="__main__":merge()
