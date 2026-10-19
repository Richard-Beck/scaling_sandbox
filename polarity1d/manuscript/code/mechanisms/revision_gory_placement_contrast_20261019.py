"""A controlled placement contrast within the existing stationary-cue grid."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import revision_standard_angles_sizes_20261019 as campaign
import revision_standard_angles_report_20261019 as report
import revision_standard_angles_duration_20261019 as native

ROOT=native.ROOT
REPO=ROOT.parents[1]
OUT=campaign.OUT/"placement_contrast"
SOURCE={}
ON="#327F8B"
FREE="#A74C84"
GREY="#666666"


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source(path):path=Path(path);SOURCE[str(path)]=sha(path);return path
def dest(name,group="data"):
    directory=OUT if name.endswith(".png") else OUT/"supporting"/group
    directory.mkdir(parents=True,exist_ok=True);return directory/name
def write(name,z):dest(name,"metadata").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
def savefig(name,fig):
    for ext in ("png","pdf","svg"):fig.savefig(dest(name+"."+ext,"vectors"),dpi=230,bbox_inches="tight")
    plt.close(fig)


def main():
    table=pd.read_csv(source(campaign.OUT/"supporting/data/stationary_sizes_summary.csv"))
    selected=table[(table.model=="goryachev")&(table.size_parameter==40)&(table.cue_duration==300)&table.angle_label.isin(["pi_4","pi_2"])].sort_values("requested_angle_deg")
    assert len(selected)==2 and selected.status.eq("complete").all()
    records=[];profiles=[];traces=[];histories={};initials=[];cues=[];maximum_replay=0.
    with np.load(source(campaign.OUT/"supporting/data/stationary_sizes_traces.npz")) as z:reference={k:z[k].copy() for k in z.files if k.startswith(("c108_","c113_"))}
    for row in selected.itertuples():
        phases=[]
        for label in ("cue","recovery"):
            path=source(getattr(row,"actual_"+label+"_path"));assert sha(path)==getattr(row,"actual_"+label+"_sha256")
            with np.load(path) as z:
                phases.append((z["time"].copy(),z["states"].copy(),z["x"].copy()))
                if label=="cue":signal=z["stimulus"].copy()
                else:assert not np.count_nonzero(z["stimulus"])
        np.testing.assert_array_equal(phases[0][1][-1],phases[1][1][0]);assert phases[0][0][-1]==300. and phases[1][0][-1]==13000.
        initial=phases[0][1][0];initials.append(initial);cues.append(signal)
        x=phases[0][2];q=x/40;n=len(q);joff=len(phases[0][0])-1
        t=np.r_[phases[0][0],phases[1][0][1:]+300];y=np.concatenate((phases[0][1],phases[1][1][1:]))
        u=y[:,0];M0=40*u.mean(axis=1);M1=40*(u@np.exp(2j*np.pi*q))/n;P=M1/M1[0];R=abs(P);angle=report.heading(P)
        core=native.core_info(q,initial[0],row.preparation_axis_fraction);w=core["weights"];A=40*(u@w);B0=40*core["background"]*w.sum();E=(A-B0)/(A[0]-B0)
        for key,value in (("time",t),("M0",M0),("P_real",P.real),("P_imag",P.imag),("R",R),("fixed_region_transported_background_excess",E)):
            error=float(abs(value-reference[f"c{row.condition_id}_{key}"]).max());assert error<1e-10;maximum_replay=max(maximum_replay,error)
        M0ratio=M0/M0[0];release=phases[1][0]
        weighted_inactive=float(signal@initial[1]/signal.sum())
        initial_transfer=float(40*np.mean(signal*initial[1]))
        model=native.ring.native_model("goryachev");grid=native.ring.RingGrid(40,n)
        increment=native.ring.native_rhs(model,grid,initial,signal)-native.ring.native_rhs(model,grid,initial,np.zeros(n))
        assert abs(40*increment[0].mean()-initial_transfer)<1e-10
        withdrawal_peaks,_=report.periodic_counts(u[joff]);final_peaks,_=report.periodic_counts(u[-1])
        context=dict(condition_id=int(row.condition_id),requested_angle_deg=float(row.requested_angle_deg),cue_duration=300.,circumference=40.,
            full_physical_support=float(row.full_physical_support),initial_angular_support_deg=float(row.initial_angular_support_deg),
            requested_peak=float(row.requested_peak),sampled_peak=float(signal.max()),physical_cue_integral=float(40*signal.mean()),
            cue_centre_fraction=float(row.cue_centre_fraction),physical_cue_centre=float(row.cue_centre_physical_at_onset),
            cue_support_core_overlap=float(row.overlap_fraction_cue_support),unit_peak_cue_weighted_initial_excess_fraction=float(row.initial_unit_peak_cue_excess_fraction_total),
            initial_input_weighted_inactive_concentration=weighted_inactive,initial_input_induced_physical_active_transfer_rate=initial_transfer,
            withdrawal_resolved_peaks=withdrawal_peaks,final_resolved_peaks=final_peaks,
            angle_withdrawal_deg=float(angle[joff]),angle_final_deg=float(angle[-1]),R_withdrawal=float(R[joff]),R_final=float(R[-1]),
            M0_ratio_withdrawal=float(M0ratio[joff]),M0_ratio_final=float(M0ratio[-1]),
            M0_ratio_release_minimum=float(M0ratio[joff:].min()),M0_ratio_release_maximum=float(M0ratio[joff:].max()),
            original_region_excess_withdrawal=float(E[joff]),original_region_excess_final=float(E[-1]),
            original_region_release_minimum=float(E[joff:].min()),original_region_release_rebound=float(E[-1]-E[joff:].min()),
            full_field_tail_status=row.full_field_tail_status,full_field_tail_excursion=float(row.final_full_field_tail_excursion),
            cue_path=row.actual_cue_path,cue_sha256=row.actual_cue_sha256,recovery_path=row.actual_recovery_path,recovery_sha256=row.actual_recovery_sha256)
        records.append(context)
        positions=(360*(q-row.preparation_axis_fraction)+90)%360-90
        physical=(40*(q-row.preparation_axis_fraction)+10)%40-10
        for state,index in (("initial",0),("withdrawal",joff),("retained",len(t)-1)):
            for j in range(n):profiles.append(dict(condition_id=row.condition_id,requested_angle_deg=row.requested_angle_deg,phase=state,time=float(t[index]),
                material_position=q[j],normal_angle_from_preparation_deg=positions[j],physical_arc_from_preparation=physical[j],active=u[index,j],
                active_change_from_initial=u[index,j]-u[0,j],initial_cue=signal[j]))
        for j in range(len(t)):traces.append(dict(condition_id=row.condition_id,requested_angle_deg=row.requested_angle_deg,time=t[j],
            phase="cue_on" if j<=joff else "cue_free",time_since_withdrawal=t[j]-300,M0=M0[j],M0_ratio=M0ratio[j],M1_real=M1[j].real,M1_imag=M1[j].imag,
            P_real=P[j].real,P_imag=P[j].imag,R=R[j],angle_raw_deg=angle[j],angle_deg=angle[j] if R[j]>=.1 else None,
            original_region_transported_background_excess=E[j],original_region_physical_active_amount=A[j]))
        histories[row.requested_angle_deg]=dict(time=t,states=y,q=q,positions=positions,joff=joff,M0ratio=M0ratio,P=P,R=R,E=E,initial_cue=signal)
    np.testing.assert_array_equal(initials[0],initials[1]);np.testing.assert_array_equal(np.roll(cues[0],100),cues[1])
    assert all(r["requested_peak"]==6.4 and r["sampled_peak"]==6.4 and r["full_physical_support"]==3. and r["cue_duration"]==300 for r in records)
    assert records[0]["physical_cue_integral"]==records[1]["physical_cue_integral"]
    assert all(r["cue_support_core_overlap"]==0. for r in records)
    assert all(r["full_field_tail_status"]=="passes_finite_tail_gate" for r in records)
    summary=pd.DataFrame(records);summary.to_csv(dest("placement_contrast_summary.csv"),index=False)
    pd.DataFrame(profiles).to_csv(dest("placement_contrast_profiles.csv"),index=False);pd.DataFrame(traces).to_csv(dest("placement_contrast_traces.csv"),index=False)
    render(summary,histories);docs(summary,maximum_replay)


def render(summary,histories):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.labelsize":9.5,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(10.7,10.4));grid=fig.add_gridspec(3,2,left=.075,right=.975,bottom=.14,top=.84,hspace=.64,wspace=.27,height_ratios=[1.0,1.15,.85])
    fig.suptitle("Same basal state and prescribed cue · Placement changes retention",fontsize=14.5,y=.985)
    fig.text(.055,.943,"Goryachev C = 40; stationary cue T = 300; identical full support 3, native coefficient peak 6.4 and spatial integral 9.6.",fontsize=10)
    fig.text(.055,.917,"The two withdrawal angles are similar (29.6° / 33.5°); the complete fields and complex activity vectors differ.",fontsize=10)
    handles=[plt.Line2D([],[],color=GREY,ls="--",label="Initial"),plt.Line2D([],[],color=ON,label="At withdrawal / cue on"),plt.Line2D([],[],color=FREE,label="Retained / cue free")]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.52,.895),ncol=3,frameon=False,fontsize=9)
    for column,target in enumerate((45.,90.)):
        h=histories[target];r=summary[summary.requested_angle_deg==target].iloc[0];joff=h["joff"]
        ax=fig.add_subplot(grid[0,column]);order=np.argsort(h["positions"])
        ax.axvspan(target-13.5,target+13.5,color="#eeeeee",zorder=0)
        ax.plot(h["positions"][order],h["states"][0,0,order],color=GREY,ls="--",lw=1)
        ax.plot(h["positions"][order],h["states"][joff,0,order],color=ON,lw=1.4)
        ax.plot(h["positions"][order],h["states"][-1,0,order],color=FREE,lw=1.4)
        ax.set(xlim=(-90,270),ylim=(-1,69),xticks=[0,45,90,180,270],xlabel="Normal angle from prepared axis (°)",ylabel="Native active concentration")
        ax.set_title(f"A{column+1}  Cue at {target:g}° · actual local profiles",loc="left",pad=10,fontsize=10.5)
        ax.spines[["top","right"]].set_visible(False)
        ax=fig.add_subplot(grid[1,column]);circle=np.linspace(0,np.pi/2,100)
        ax.plot(np.cos(circle),np.sin(circle),color="#bbbbbb",ls=":",lw=.8)
        ax.axhline(0,color="#dddddd",lw=.7);ax.axvline(0,color="#dddddd",lw=.7)
        P=h["P"]
        ax.plot(P[:joff+1].real,P[:joff+1].imag,color=ON,lw=1.5)
        ax.plot(P[joff:].real,P[joff:].imag,color=FREE,lw=1.5)
        for part,color in ((P[:joff+1],ON),(P[joff:],FREE)):
            path=np.r_[0,np.cumsum(abs(np.diff(part)))];mid=np.argmin(abs(path-path[-1]*.55));end=min(len(part)-1,mid+max(1,len(part)//30))
            ax.annotate("",xy=(part[end].real,part[end].imag),xytext=(part[mid].real,part[mid].imag),arrowprops=dict(arrowstyle="->",color=color,lw=1.25))
        ax.plot(P[0].real,P[0].imag,"o",color=GREY,ms=4)
        ax.plot(P[joff].real,P[joff].imag,"o",mfc="white",mec=ON,ms=6,mew=1.4,zorder=5)
        ax.plot(P[-1].real,P[-1].imag,"D",color=FREE,ms=5,zorder=6)
        ax.plot([0,np.cos(np.deg2rad(target))],[0,np.sin(np.deg2rad(target))],color="#cccccc",ls="--",lw=.7)
        ax.set(xlim=(-.07,1.07),ylim=(-.07,1.07),xticks=[0,.5,1],yticks=[0,.5,1],xlabel="Re P",ylabel="Im P")
        ax.set_aspect("equal",adjustable="box");ax.spines[["top","right"]].set_visible(False)
        ax.set_title(f"B{column+1}  {r.angle_withdrawal_deg:.1f}° → {r.angle_final_deg:.2f}°;  R {r.R_withdrawal:.3f} → 1",loc="left",pad=8,fontsize=10)
        ax=fig.add_subplot(grid[2,column]);tt=h["time"][joff:]-300
        ax.plot(tt,h["E"][joff:],color="#333333",lw=1.5,label="Original-region excess")
        ax.plot(tt,h["R"][joff:],color="#4477AA",lw=1.2,label="R")
        ax.plot(tt,h["M0ratio"][joff:],color="#9b8250",ls="--",lw=1.1,label="M₀ / initial M₀")
        ax.set_xscale("symlog",linthresh=2.)
        ax.set(xlim=(0,13000),ylim=(-.025,1.06),xlabel="Cue-free elapsed time (native units)",ylabel="Relative activity / amount")
        ax.set_xticks([0,100,1000,13000]);ax.set_xticklabels(["0","100","1000","13000"]);ax.spines[["top","right"]].set_visible(False)
        ax.set_title(f"C{column+1}  Original-region excess {r.original_region_excess_withdrawal:.3f} → {r.original_region_excess_final:.3f}",loc="left",pad=7,fontsize=10)
        if column==0:ax.legend(frameon=False,fontsize=7.7,loc="center left",bbox_to_anchor=(.01,.55))
    fig.text(.055,.075,"A: full-ring profiles; grey band marks input support. B: equal scales; circle = withdrawal, diamond = +13,000. C: time is linear near zero, logarithmic later.",fontsize=8.6)
    fig.text(.055,.049,"Only cue placement changes; withdrawal states differ. Both initial supports miss the half-contrast core, with different small basal tails.",fontsize=8.6)
    fig.text(.055,.023,"Regional activity does not track site lineage. This selected contrast does not identify a general fate predictor, basin boundary or unique mechanism.",fontsize=8.4)
    savefig("gory_same_basal_placement",fig)


def docs(summary,error):
    generated=datetime.now(timezone.utc).isoformat()
    caption="Controlled placement contrast in the same qualified Goryachev C = 40 basal state. The stationary native activation-rate coefficient has identical full physical support 3, peak 6.4, spatial integral 9.6 and duration 300; the cue centre alone changes from 45° to 90°. This matches prescribed input-coefficient exposure, rather than activated biochemical amount: the native channel is stimulus × local inactive concentration. Both supports are disjoint from the original basal half-contrast core, with measured nonzero basal-tail weighting retained. A, Actual initial, withdrawal and retained primary fields on one full turn use a common native concentration scale. B, Complex P = M₁(t)/M₁(0) trajectories use identical real/imaginary scales; teal is cue on, magenta is cue free, the open circle is withdrawal and the diamond is the endpoint after 13,000 cue-free units. Similar withdrawal headings (29.59° versus 33.48°) conceal different vector magnitudes (.924 versus .721) and different local organizations. During release, the 45° condition continues to 44.10°, whereas the 90° condition returns to −.028°. C, Physical excess activity in the fixed original basal material region, R = |P| and physical total active amount M₀/M₀(0) are plotted through cue-free relaxation. The time axis is symlog: linear within two native units of withdrawal and logarithmic thereafter. Both endpoint fields pass the independent finite-tail gate. Similar phase at withdrawal does not imply matched full states or matched P, and this selected pair does not identify a general fate predictor or unique competition mechanism. Activity orientation and local redistribution do not establish site lineage or cell movement."
    for a,b in (("from45° to90°","from 45° to 90°"),("after13,000","after 13,000"),("versus33.48°","versus 33.48°"),("versus.721","versus .721"),("the45°","the 45°"),("the90°","the 90°"),("to44.10°","to 44.10°"),("to−.028°","to −.028°")):caption=caption.replace(a,b)
    dest("CAPTION.md","metadata").write_text(caption+"\n")
    a,b=summary.iloc[0],summary.iloc[1]
    evidence=f"""# Same-basal placement contrast

Generated: {generated}. The selected pair comes from the already completed240-condition stationary-cue experiment, without new simulations. Both cases use Goryachev C = 40 and cue duration300; centres45° and90° are conditions113 and108. Selection is explicit: their withdrawal angles are similar, while their cue-free fates differ markedly. This illustrates information lost by reducing the state to its phase at withdrawal. It is a selected teaching contrast, not a matched-state experiment or a statistical predictor.

The two full native initial states are exactly equal. The actual cue arrays are exact cyclic translations by 100 of 800 cells, so coefficient peak, support, integral and shape are equal. Both have analytic/sample coefficient peak 6.4, physical full support 3, angular full support 27° and physical input-coefficient spatial integral 9.6. Cue duration 300 and cue-free horizon 13,000 are identical, with exact cue/release state continuity. This matches prescribed input-coefficient exposure rather than activated biochemical amount, because the native activation channel is stimulus × local inactive concentration. Initial cue-weighted inactive concentrations are {a.initial_input_weighted_inactive_concentration:.8f} and {b.initial_input_weighted_inactive_concentration:.8f}; the corresponding integrated initial active-transfer rates are {a.initial_input_induced_physical_active_transfer_rate:.8f} and {b.initial_input_induced_physical_active_transfer_rate:.8f}. These native RHS increments are checked directly. The cue centre is defined from the prepared axis; the measured initial M₁ phase offset is retained in the authoritative 240 table. The two inputs have zero geometric half-contrast-core intersection. Their unit-peak cue-weighted basal excess fractions are {a.unit_peak_cue_weighted_initial_excess_fraction:.6g} and {b.unit_peak_cue_weighted_initial_excess_fraction:.6g}; both are small but not declared identical or zero.

At withdrawal, the 45° condition has angle {a.angle_withdrawal_deg:.6f}°, R {a.R_withdrawal:.6f}, total active amount ratio {a.M0_ratio_withdrawal:.6f} and original-region excess {a.original_region_excess_withdrawal:.6f}. The 90° condition has angle {b.angle_withdrawal_deg:.6f}°, R {b.R_withdrawal:.6f}, amount ratio {b.M0_ratio_withdrawal:.6f} and regional excess {b.original_region_excess_withdrawal:.6f}. Both withdrawal profiles have two prominence-qualified peaks; both final profiles have one. These are different physical states and different complex P values, despite similar headings.

After 13,000 cue-free native time units, retained angles are {a.angle_final_deg:.6f}° and {b.angle_final_deg:.6f}°, both with R ≈ 1 and total activity returned to basal. Original-region excess is {a.original_region_excess_final:.8g} versus {b.original_region_excess_final:.8g}. In the 90° case, regional excess rebounds by {b.original_region_release_rebound:.6f} from its cue-free sampled minimum; the 45° case instead continues losing regional excess. During release the total active amount ratios span [{a.M0_ratio_release_minimum:.6f}, {a.M0_ratio_release_maximum:.6f}] and [{b.M0_ratio_release_minimum:.6f}, {b.M0_ratio_release_maximum:.6f}]. These changes in total activity cannot by themselves represent the distinct local distributions and resultant-vector trajectories.

Raw physical M₁ components use fixed ring reference axes. Dividing by the complex initial M₁ rescales and rotates axes once, giving P(0) = 1 + 0i. P is not a complex notation for cell velocity. M₀ is the physical primary active amount; it differs from conserved total protein including inactive fields. The regional excess uses the frozen basal half-contrast material interval and transported initial background amount, matching the full240 experiment. Its change can reflect redistribution or local chemistry and does not identify ancestry, molecular destruction or a unique competition mechanism.

Native phase hashes and exact phase continuity are checked. Directly reconstructed M₀, P, R and regional traces match the authoritative 240 traces to maximum discrepancy {error:.3g}. Both full fields pass their independently measured final 1000-unit tail gate. Panel C uses a symlog time axis, linear within two native units of withdrawal and logarithmic thereafter; no time averaging is applied. Finite-window retention does not prove asymptotic attraction, and the pair alone does not locate a basin boundary or establish how untested placements behave.
"""
    replacements={"completed240":"completed 240","duration300":"duration 300","centres45° and90°":"centres 45° and 90°","conditions113 and108":"conditions 113 and 108", "by100 of800":"by 100 of 800","peak6.4":"peak 6.4","support3":"support 3","support27°":"support 27°","integral9.6":"integral 9.6","horizon13,000":"horizon 13,000","authoritative240":"authoritative 240","are{":"are {","The90°":"The 90°","the90°":"the 90°","the45°":"the 45°","After13,000":"After 13,000","full240":"full 240","final1000":"final 1000","maximum{":"maximum {"}
    for x,y in replacements.items():evidence=evidence.replace(x,y)
    dest("EVIDENCE.md","metadata").write_text(evidence)
    (OUT/"README.md").write_text("""# Controlled same-basal placement contrast

- [Review sheet](gory_same_basal_placement.png) · [PDF](supporting/vectors/gory_same_basal_placement.pdf)
- [Exact endpoint/input comparison](supporting/data/placement_contrast_summary.csv)
- [Actual native profiles](supporting/data/placement_contrast_profiles.csv)
- [Reconstructed physical/vector/regional trajectories](supporting/data/placement_contrast_traces.csv)
- [Caption](supporting/metadata/CAPTION.md) · [Selection rationale and evidence](supporting/metadata/EVIDENCE.md)

These two existing Goryachev C40 conditions share the exact basal state and prescribed input-coefficient profile/duration. Similar withdrawal angles lead to further rotation at 45° versus return at 90°. The complete activity vector and native local fields distinguish the states. The pair is selected to teach this distinction and does not identify a universal mechanism or endpoint predictor.
""")
    write("validation.json",dict(generated_at=generated,passed=True,conditions=2,new_integrations=0,condition_ids=[113,108],
        exact_same_full_basal_state=True,exact_cue_rotation_cells=100,identical_sampled_peak_support_integral_duration=True,
        both_zero_initial_core_overlap=True,exact_phase_continuity=True,both_full_field_tail_gates_pass=True,maximum_authoritative_trace_replay_error=error))
    recordpaths=[Path(__file__),Path(native.__file__),Path(report.__file__),Path(campaign.__file__),ROOT/"code/core/models/goryachev.py",ROOT/"code/mechanisms/ring_distance_prototype_20261017.py"]
    for path in recordpaths:source(path)
    write("source_sha256.json",SOURCE)
    write("content_checkpoint.json",dict(title="Controlled cue-placement contrast",role="Selected existing pair within the stationary experiment; no new protocol or final layout.",
        assets=[dict(path=str((OUT/"gory_same_basal_placement.png").relative_to(REPO)),caption="Same basal and prescribed coefficient profile, different placement: similar withdrawal angles can precede distinct retained states.")],legend=[caption],
        sources=[dict(label="Controlled native inputs and endpoints",path=str(dest("placement_contrast_summary.csv").relative_to(REPO))),dict(label="Selection rationale, vector/region definitions and limits",path=str(dest("EVIDENCE.md","metadata").relative_to(REPO)))]))
    print(json.dumps(dict(passed=True,conditions=2,new_integrations=0,maximum_trace_replay_error=error,output=str(OUT)),indent=2))


if __name__=="__main__":main()
