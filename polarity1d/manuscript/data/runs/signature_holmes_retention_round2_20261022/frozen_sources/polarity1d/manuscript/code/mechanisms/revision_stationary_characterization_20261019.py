"""Characterize cue-on response versus cue-free retention in the same 80 histories."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import revision_standard_angles_duration_20261019 as native
import revision_standard_angles_report_20261019 as prior

ROOT=native.ROOT
REPO=ROOT.parents[1]
RUN=native.DEFAULT
OUT=native.TOPIC.parent/"stationary_characterization"
SOURCE={}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p):
    p=Path(p);SOURCE[str(p)]=sha(p);return p
def dest(name,group="data"):
    folder=OUT if name.endswith(".png") else OUT/"supporting"/group
    folder.mkdir(parents=True,exist_ok=True);return folder/name
def write(name,z):dest(name,"metadata").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
def savefig(stem,fig):
    for ext in ("png","pdf","svg"):fig.savefig(dest(stem+"."+ext,"vectors"),dpi=220,bbox_inches="tight")
    plt.close(fig)


def wrap_error(angle,target):
    return np.rad2deg(np.angle(np.exp(1j*np.deg2rad(angle-target))))


def first_bracket(t,mask):
    j=np.flatnonzero(mask)
    if not len(j):return dict(status="not_observed_by_phase_horizon",lower=float(t[-1]),upper=None)
    k=int(j[0]);return dict(status="already_true_at_phase_start" if k==0 else "saved_sample_bracket",lower=float(t[k-1]) if k else float(t[0]),upper=float(t[k]))


def retained_bracket(t,mask):
    if not mask[-1]:return dict(status="not_retained_at_phase_horizon",lower=float(t[-1]),upper=None)
    lastfalse=np.flatnonzero(~mask)
    k=int(lastfalse[-1]+1) if len(lastfalse) else 0
    return dict(status="true_at_all_saved_phase_frames" if k==0 else "saved_sample_bracket",lower=float(t[k-1]) if k else float(t[0]),upper=float(t[k]))


def analyze():
    original_path=native.TOPIC/"supporting/data/standard_angles_summary.csv"
    original=pd.read_csv(record(original_path));manifest=json.loads(record(RUN/"manifest.json").read_text())
    with np.load(record(native.TOPIC/"supporting/data/standard_angles_traces.npz")) as z:reference={k:z[k].copy() for k in z.files}
    rows=[];events=[];alltraces={};checks=[];max_error=0.
    for case in manifest["cases"]:
        cid=case["condition_id"];row=original[original.condition_id==cid].iloc[0].to_dict();spring=case["model"]=="spring"
        Cb=case["circumference"]
        assert sha(RUN/case["input"])==case["input_sha256"]
        with np.load(RUN/case["input"]) as z:initial=z["state"].copy()
        parts=[]
        for label in ("cue","recovery"):
            path=record(row["actual_"+label+"_path"]);assert sha(path)==row["actual_"+label+"_sha256"]
            parts.append(prior.load_phase(path,spring,Cb))
        q=parts[0]["q"];n=len(q);u0=initial[:n] if spring else initial[0]
        assert np.array_equal(parts[0]["u"][0],u0) and np.array_equal(parts[0]["u"][-1],parts[1]["u"][0])
        T=case["cue_duration"];joff=len(parts[0]["time"])-1
        t=np.r_[parts[0]["time"],parts[1]["time"][1:]+T]
        u=np.vstack((parts[0]["u"],parts[1]["u"][1:]));C=np.r_[parts[0]["C"],parts[1]["C"][1:]]
        M0=C*u.mean(axis=1);M1=C*(u@np.exp(2j*np.pi*q))/n;M2=C*(u@np.exp(4j*np.pi*q))/n
        P=M1/M1[0];R=abs(P);angle=prior.heading(P);q1=abs(M1)/M0;q2=abs(M2)/M0
        core=native.core_info(q,u0,case["preparation_axis_fraction"]);w=core["weights"]
        core_amount=C*(u@w);background0=Cb*core["background"]*w.sum();initial_excess=core_amount[0]-background0
        E=(core_amount-background0)/initial_excess;M0ratio=M0/M0[0]
        for key,value in (("time",t),("P_real",P.real),("P_imag",P.imag),("R",R),("q1",q1),("q2",q2),("fixed_region_transported_background_excess",E),("C",C)):
            error=float(abs(value-reference[f"c{cid}_{key}"]).max());assert error<1e-10;max_error=max(max_error,error)
        row.update(M0_initial=float(M0[0]),M1_initial_real=float(M1[0].real),M1_initial_imag=float(M1[0].imag),
            M0_withdrawal=float(M0[joff]),M0_final=float(M0[-1]),M0_ratio_withdrawal=float(M0ratio[joff]),M0_ratio_final=float(M0ratio[-1]),
            M0_ratio_sampled_minimum=float(M0ratio.min()),M0_ratio_sampled_maximum=float(M0ratio.max()),
            M0_ratio_sampled_minimum_time=float(t[M0ratio.argmin()]),M0_ratio_sampled_maximum_time=float(t[M0ratio.argmax()]),
            M0_ratio_cue_on_sampled_minimum=float(M0ratio[:joff+1].min()),M0_ratio_cue_on_sampled_maximum=float(M0ratio[:joff+1].max()),
            R_sampled_minimum=float(R.min()),R_sampled_minimum_time=float(t[R.argmin()]),R_cue_on_sampled_minimum=float(R[:joff+1].min()),
            M0_ratio_at_R_sampled_minimum=float(M0ratio[R.argmin()]),q2_at_R_sampled_minimum=float(q2[R.argmin()]),
            R_release_sampled_minimum=float(R[joff:].min()),R_sampled_maximum=float(R.max()),
            q1_withdrawal=float(q1[joff]),q1_final=float(q1[-1]),q2_withdrawal=float(q2[joff]),q2_final=float(q2[-1]),
            angle_withdrawal_raw_deg=float(angle[joff]),angle_final_raw_deg=float(angle[-1]),
            angle_withdrawal_deg=float(angle[joff]) if R[joff]>=.1 else None,
            angle_final_deg=float(angle[-1]) if R[-1]>=.1 else None,
            release_signed_turn_deg=float(wrap_error(angle[-1],angle[joff])) if min(R[joff],R[-1])>=.1 else None,
            release_target_error_change_deg=float(abs(wrap_error(angle[-1],case["target_angle_in_P_coordinates_deg"]))-abs(wrap_error(angle[joff],case["target_angle_in_P_coordinates_deg"]))) if min(R[joff],R[-1])>=.1 else None)
        row.update(fixed_region_excess_initial=1.,fixed_region_excess_withdrawal=float(E[joff]),fixed_region_excess_final=float(E[-1]),
            fixed_region_excess_cue_on_sampled_minimum=float(E[:joff+1].min()),
            fixed_region_excess_cue_on_sampled_minimum_time=float(t[:joff+1][E[:joff+1].argmin()]),
            fixed_region_excess_release_sampled_minimum=float(E[joff:].min()),
            fixed_region_excess_release_sampled_minimum_time=float(t[joff:][E[joff:].argmin()]),
            fixed_region_excess_sampled_minimum=float(E.min()),fixed_region_excess_sampled_minimum_time=float(t[E.argmin()]),
            fixed_region_release_net_change=float(E[-1]-E[joff]),fixed_region_rebound_from_sampled_minimum=float(E[-1]-E.min()),
            fixed_region_physical_amount_withdrawal=float(core_amount[joff]),fixed_region_physical_amount_final=float(core_amount[-1]),
            physical_C_withdrawal=float(C[joff]),physical_C_final=float(C[-1]),physical_C_min=float(C.min()),physical_C_max=float(C.max()))
        for threshold in (.9,.5,.1):
            for phase,keep in (("cue_on",t<=T),("release",t>=T),("entire_history",np.ones(len(t),bool))):
                tt=t[keep];ee=E[keep];mask=ee<=threshold;cross=np.flatnonzero(mask)
                meta=dict(condition_id=cid,model=case["model"],angle_label=case["angle_label"],requested_angle_deg=case["requested_angle_deg"],cue_duration=T,
                    initial_overlap=case["initial_overlap"],phase=phase,threshold=threshold,phase_start=float(tt[0]),phase_horizon=float(tt[-1]),
                    sampled_minimum=float(ee.min()),sampled_minimum_time=float(tt[ee.argmin()]),
                    rises_above_threshold_after_first_saved_entry=bool(np.any(~mask[int(cross[0]):])) if len(cross) else False,
                    fraction_at_phase_end=float(ee[-1]))
                events.append(dict(**meta,event="first_saved_entry",**first_bracket(tt,mask)))
                events.append(dict(**meta,event="entry_retained_through_all_later_saved_phase_frames",**retained_bracket(tt,mask)))
        rows.append(row)
        traces=dict(time=t,phase=np.where(t<=T,0,1),M0=M0,M0_ratio=M0ratio,M1_real=M1.real,M1_imag=M1.imag,
            P_real=P.real,P_imag=P.imag,R=R,angle_raw_deg=angle,angle_masked_deg=np.where(R>=.1,angle,np.nan),q1=q1,q2=q2,
            fixed_region_transported_background_excess=E,fixed_region_physical_active_amount=core_amount,C=C,
            prescribed_physical_offset=case["fraction"]*C)
        alltraces.update({f"c{cid}_{key}":value for key,value in traces.items()})
        checks.append(dict(condition_id=cid,actual_phase_hashes_match=True,exact_primary_restart_and_continuity=True,maximum_prior_trace_replay_error=max_error))
    frame=pd.DataFrame(rows);frame.to_csv(dest("stationary_whole_grid_summary.csv"),index=False)
    pd.DataFrame(events).to_csv(dest("stationary_fixed_region_crossings.csv"),index=False)
    np.savez_compressed(dest("stationary_whole_grid_traces.npz"),**alltraces)
    compact=["condition_id","model","angle_label","cue_duration","initial_overlap","angle_withdrawal_deg","angle_final_deg","R_withdrawal","R_final","R_sampled_minimum",
        "M0_ratio_withdrawal","M0_ratio_final","M0_ratio_sampled_minimum","M0_ratio_sampled_maximum","fixed_region_excess_withdrawal","fixed_region_excess_final",
        "fixed_region_excess_cue_on_sampled_minimum","fixed_region_rebound_from_sampled_minimum","release_signed_turn_deg","full_field_tail_status"]
    frame[compact].to_csv(dest("stationary_endpoint_pairs.csv"),index=False)
    render(frame);docs(frame,max_error,checks)


def setup():plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.labelsize":10,"pdf.fonttype":42,"svg.fonttype":"none"})


def phase_curves(ax,g,off,final,color):
    ax.plot(g.cue_duration,g[off],"o--",color=color,lw=1,ms=4,mfc="white",mew=1)
    ax.plot(g.cue_duration,g[final],"o-",color=color,lw=1.35,ms=4)
    for _,r in g.iterrows():
        if pd.notna(r[off]) and pd.notna(r[final]):ax.plot([r.cue_duration]*2,[r[off],r[final]],color=color,lw=.65,alpha=.65)


def common(ax,model,ylabel):
    ax.set(xscale="log",xlim=(.07 if model=="spring" else .7,1500),xlabel="Cue duration (native units)",ylabel=ylabel)
    ax.spines[["top","right"]].set_visible(False);ax.tick_params(labelsize=9)


def legend(fig,top):
    anglehandles=[plt.Line2D([],[],color=prior.COLORS[label],lw=1.5,label=prior.LABELS[label]) for label,_ in native.ANGLES]
    phasehandles=[plt.Line2D([],[],color="#555555",marker="o",mfc="white",ls="--",ms=4,label="At withdrawal"),plt.Line2D([],[],color="#555555",marker="o",ls="-",ms=4,label="After cue-free retention")]
    fig.legend(handles=anglehandles,loc="upper center",bbox_to_anchor=(.52,top),ncol=4,frameon=False,fontsize=10)
    fig.legend(handles=phasehandles,loc="upper center",bbox_to_anchor=(.52,top-.045),ncol=2,frameon=False,fontsize=9)


def render(frame):
    setup();fig,axes=plt.subplots(2,2,figsize=(10.8,8.7))
    fig.subplots_adjust(left=.08,right=.975,bottom=.15,top=.78,hspace=.4,wspace=.27)
    fig.suptitle("Cue-on deflection and cue-free retention are different readouts",fontsize=15,y=.985)
    fig.text(.055,.939,"The same 80 stationary-cue histories: open markers at withdrawal; filled markers after 13,000 cue-free units.",fontsize=10)
    legend(fig,.91)
    for i,model in enumerate(native.MODELS):
        ax=axes.flat[i]
        for label,_ in native.ANGLES:
            g=frame[(frame.model==model)&(frame.angle_label==label)].sort_values("cue_duration")
            phase_curves(ax,g,"angle_withdrawal_deg","angle_final_deg",prior.COLORS[label])
            ax.axhline(g.target_angle_in_P_coordinates_deg.iloc[0],color=prior.COLORS[label],ls=":",lw=.6,alpha=.55)
            failed=g[g.full_field_tail_status=="fails_finite_tail_gate"]
            ax.scatter(failed.cue_duration,failed.angle_final_deg,marker="s",facecolor="none",edgecolor=prior.COLORS[label],s=60,lw=.8,zorder=6)
        common(ax,model,"Activity angle relative to initial M₁ (°)")
        ax.set(ylim=(-9,190),yticks=[0,22.5,45,90,135,180]);ax.set_title(f"{chr(65+i)}  {prior.panel_title(model)}",loc="left",fontsize=10.5,pad=10)
        if model=="goryachev":ax.annotate("33.5° → 0°",xy=(300,-.03),xytext=(35,57),fontsize=8.5,arrowprops=dict(arrowstyle="-",lw=.6,color="#777777"))
        if model=="holmes_model3":ax.text(.03,.78,"180°, T1000: both angles masked",transform=ax.transAxes,fontsize=8.3,color="#555555")
    fig.text(.055,.081,"Thin vertical lines pair the same duration. Dotted lines mark exact prescribed cue normals; connecting curves guide the eye between tested points.",fontsize=8.7)
    fig.text(.055,.054,"Displayed angle is masked at R < 0.1. Square outline: retained full-field tail gate failed. Four compact reuses have no full-field tail qualification.",fontsize=8.7)
    fig.text(.055,.027,"All 35 initial overlaps are retained. Native widths, peaks and clocks remain model-specific; this activity-vector readout does not specify cell motion.",fontsize=8.7)
    savefig("stationary_withdrawal_retention",fig)
    fig,axes=plt.subplots(4,2,figsize=(10.8,12.7))
    fig.subplots_adjust(left=.08,right=.975,bottom=.1,top=.84,hspace=.62,wspace=.27)
    fig.suptitle("Magnitude and total activity during the same stationary-cue experiment",fontsize=14.5,y=.985)
    fig.text(.055,.95,"Both quantities use physical amount, including Spring's changing circumference. Endpoint pairs are shown against sampled trajectory ranges.",fontsize=9.2)
    legend(fig,.927)
    for i,model in enumerate(native.MODELS):
        for j,(off,final,low,high,ylabel,title) in enumerate((
            ("R_withdrawal","R_final","R_sampled_minimum","R_sampled_maximum","R = |M₁/M₁(initial)|","Relative first-moment magnitude"),
            ("M0_ratio_withdrawal","M0_ratio_final","M0_ratio_sampled_minimum","M0_ratio_sampled_maximum","M₀/M₀(initial)","Total active amount"))):
            ax=axes[i,j]
            for label,_ in native.ANGLES:
                g=frame[(frame.model==model)&(frame.angle_label==label)].sort_values("cue_duration");color=prior.COLORS[label]
                ax.fill_between(g.cue_duration,g[low],g[high],color=color,alpha=.07,lw=0)
                phase_curves(ax,g,off,final,color)
            ax.axhline(1,color="#aaaaaa",ls=":",lw=.7);common(ax,model,ylabel)
            ax.set_title(f"{NAMES[model]} · {title}",loc="left",fontsize=9.8,pad=7)
    fig.text(.055,.048,"Shading spans sampled minima and maxima over cue-on plus release; it is not uncertainty. Filled final points are distinct from cue-on activity growth.",fontsize=8.5)
    fig.text(.055,.024,"Activity amount and first-moment magnitude are separate quantities: opposing regions can reduce the resultant without depleting total activity.",fontsize=8.5)
    savefig("stationary_amount_magnitude",fig)


NAMES=prior.NAMES


def docs(frame,error,checks):
    generated=datetime.now(timezone.utc).isoformat();events=pd.read_csv(dest("stationary_fixed_region_crossings.csv"))
    caption="Cue withdrawal versus retained activity orientation in the same standardized-angle stationary-cue experiment. A–D, The open point is measured at cue withdrawal and the filled point after 13,000 cue-free native time units; thin vertical segments connect each pair, dashed curves connect withdrawal points and solid curves connect retained points. Four shared colours denote prepared-axis cue-centre offsets 180°/90°/45°/22.5°; dotted lines mark the exact targets in P coordinates. Every one of the 80 conditions remains, including 35 initial cue–core overlaps, with the exact same basals, fixed physical cue supports and native dose/time protocols as the standardized-angle checkpoint. The displayed activity angle is arg P for physical P=M₁(t)/M₁(0), masked at R < 0.1. In the paired magnitude/amount support sheet, R is distinguished from physical total activity M₀(t)/M₀(0); translucent envelopes span sampled minima/maxima over the full history and are not error bands. The square outline flags the independently measured full-field tail-gate failure at Holmes 180° / T300; four compact reuses have no full-field qualification. Cue-on deflection can relax back toward the original state or progress further toward the cue after withdrawal. No refined minimum-duration claim, ancestry inference or movement law is assigned."
    dest("CAPTION.md","metadata").write_text(caption+"\n")
    examples=frame[((frame.model=="goryachev")&(frame.angle_label.isin(["pi_2","pi_4"]))&(frame.cue_duration==300))|((frame.model=="wave_pinning")&(frame.angle_label=="pi_2")&(frame.cue_duration==10))|((frame.model=="spring")&(frame.angle_label=="pi_2")&(frame.cue_duration==1))]
    exampletext="\n".join(f"- {NAMES[r.model]}, target {r.requested_angle_deg:g}°, duration {r.cue_duration:g}: withdrawal {r.angle_withdrawal_deg:.3f}° → retained {r.angle_final_deg:.3f}°; R {r.R_withdrawal:.3f} → {r.R_final:.3f}." for r in examples.itertuples())
    evidence=f"""# Cue-on response and cue-free retention in the same 80 histories

Generated: {generated}. This is a complementary analysis of the existing 80 standardized-angle histories. No new integrations, sizes, cue protocols or moving cues are introduced. The exact existing initial overlaps, physical circumference histories, continuous target phases, native archive provenance and independent field-tail flags remain attached to every row.

The primary sheet directly pairs withdrawal and retained angles at each tested duration. The paired magnitude/amount sheet shows why activity-vector magnitude and total active amount must remain distinct. The supplied tables contain sampled extrema, fixed-region depletion/rebound and threshold brackets so these quantities need not be inferred from a single endpoint.

## Observed examples

{exampletext}

These examples show that cue-on angular deflection is not a retained-angle surrogate. Post-withdrawal evolution may return toward the original direction or continue toward the input normal. A first moment alone does not identify site lineage: the existing native final-profile checks show all 60 non-antipodal endpoints are one translated full-strength basal region, while Holmes antipodal cancellation requires its separate local-field evidence.

## Measured quantities

The raw physical first moment M₁ has real and imaginary components along fixed ring reference axes. Normalizing by the complex initial M₁ both rescales the vector and rotates its reference axes once by the initial M₁ phase, placing P(0) at 1 + 0i. P therefore measures orientation relative to the initial resultant, with no subsequent rotating reference frame. Its magnitude depends on both total activity and spatial alignment. M₀ is physical total active amount; the native API's M₁/M₀ differs from P = M₁/M₁(0). We export physical M₀ and M₁ directly from every native active field, plus both normalizations. Spring uses instantaneous C(t), and its prescribed material cue offset remains fixed while its physical offset changes as fraction × C(t).

For every history the compact summary exports withdrawal/final angles and R, physical M₀ ratios at withdrawal/final and sampled min/max, sampled minimum R, cue-on and release fixed-original-region minima, final regional value, net release change and rebound from the observed minimum. “Minimum” here always means minimum over saved samples, not a proven continuous-time extremum. Raw angles are preserved separately; plotted angles are masked below R = 0.1.

The fixed-original-region observable is the physical active amount in the frozen basal half-contrast material interval, minus its transported initial background amount, divided by initial excess. Pure affine dilution leaves it unchanged. It can change through profile translation, redistribution or local chemistry, so depletion is not equivalent to destruction or death of an independently tracked incumbent. Overlapping inputs cannot be relabelled as a separate challenger. Rebound is quantified as final excess minus its sampled minimum and as final minus withdrawal; threshold recrossing is explicit.

The crossing table has {len(events)} records: 80 conditions × 3 levels (0.9/0.5/0.1 of initial regional excess) × 3 phases (cue-on/release/full) × 2 events (first saved entry or entry retained through all later saved phase frames). Times use cue-onset origin. Each first/retained entry has a saved-sample bracket or explicit finite-horizon censoring, and the table reports whether the region rises above that threshold after its first sampled entry. Saved sampling can miss shorter excursions; these brackets do not establish a minimum effective cue duration.

The independent existing whole-grid trace replay agrees to {error:.3g}. Actual native phase hashes and exact primary restart/continuity are checked again. Prior full-field tail flags remain 75 passes / 1 failure / 4 unavailable; activity-angle behaviour is not used to overwrite them. Native physical input full widths stay 3 for chemical models and 0.2 for Spring. Shared centre angles do not match physical or angular footprint, dose or native clocks across models.

## Outputs

- `stationary_endpoint_pairs.csv`: compact paired observables across all 80 conditions.
- `stationary_whole_grid_summary.csv`: all paired and transient observables, exact native provenance, anchor/overlap/geometry and field-tail flags.
- `stationary_whole_grid_traces.npz`: keys `c{{condition_id}}_{{observable}}`, including direct physical M0, M1_real, M1_imag, P_real, P_imag, R, q1, q2, C, fixed-region amount/excess and masked/raw angles.
- `stationary_fixed_region_crossings.csv`: all censoring-aware saved-sample brackets and recrossing flags.

The selected local-field/vector sketch produced by the other figure worker is a separate view of these same histories. Existing main/prototype art remains preserved; this checkpoint does not choose the final Figure 3 composition.
"""
    # Keep prose readable after inserting compact variable names above.
    for a,b in (("same 80","same 80"),("all80","all 80"),("all 60","all 60"),("80conditions","80 conditions"),("35initial","35 initial"),("C(t)","C(t)")):
        evidence=evidence.replace(a,b)
    dest("EVIDENCE.md","metadata").write_text(evidence)
    with dest("EVIDENCE.md","metadata").open("a") as f:
        f.write(r"""

$$
M_n(t)=C(t)\int_0^1u(q,t)e^{2\pi i n q}\,\mathrm{d}q,
\quad P(t)=\frac{M_1(t)}{M_1(0)},
\quad R(t)=|P(t)|,
\quad q_1(t)=\frac{|M_1(t)|}{M_0(t)}.
\tag{1}
$$

$$
E_{\Omega_0}^{\mathrm{tr}}(t)=
\frac{C(t)\int_{\Omega_0}u(q,t)\,\mathrm{d}q-C_{\mathrm b}u_{\mathrm{bg},0}|\Omega_0|}
     {C_{\mathrm b}\int_{\Omega_0}[u(q,0)-u_{\mathrm{bg},0}]\,\mathrm{d}q}.
\tag{2}
$$
""")
    (OUT/"README.md").write_text("""# Stationary-cue response and retention

This characterizes the same 80 standardized-angle histories without new integration or a new protocol.

- [Paired withdrawal and retained orientation](stationary_withdrawal_retention.png) · [PDF](supporting/vectors/stationary_withdrawal_retention.pdf)
- [Activity-vector magnitude and total active amount](stationary_amount_magnitude.png)
- [Compact paired endpoint table](supporting/data/stationary_endpoint_pairs.csv)
- [Full summary with exact source paths and flags](supporting/data/stationary_whole_grid_summary.csv)
- [Fixed-original-region crossings and rebounds](supporting/data/stationary_fixed_region_crossings.csv)
- [Caption](supporting/metadata/CAPTION.md) · [Evidence and definitions](supporting/metadata/EVIDENCE.md) · [Validation](supporting/metadata/validation.json)

Cue-on deflection can disappear or progress during cue-free relaxation. The complex activity vector, its magnitude and total active amount provide complementary measurements. Fixed-region depletion and rebound quantify regional activity without assigning site ancestry. All initial overlaps remain included; no minimum-effective-duration claim is introduced.
""")
    write("validation.json",dict(generated_at=generated,passed=True,conditions=80,new_integrations=0,initial_overlap_conditions=int(frame.initial_overlap.sum()),
        maximum_prior_trace_replay_error=error,all_actual_phase_hashes_match=True,all_exact_primary_restarts=True,threshold_records=len(events),
        inherited_full_field_tail_counts=frame.full_field_tail_status.value_counts().to_dict(),records=checks))
    content=dict(title="Stationary-cue response versus retained activity",role="Complementary whole-grid characterization of the same 80conditions; composition remains open.",
        assets=[dict(path=str((OUT/"stationary_withdrawal_retention.png").relative_to(REPO)),caption="Direct withdrawal/retention pairs at every requested angle and duration.")],legend=[caption],
        sources=[dict(label="All 80 paired observables and provenance",path=str(dest("stationary_whole_grid_summary.csv").relative_to(REPO))),dict(label="Regional threshold crossings and rebound",path=str(dest("stationary_fixed_region_crossings.csv").relative_to(REPO))),dict(label="Definitions and evidence",path=str(dest("EVIDENCE.md","metadata").relative_to(REPO)))],
        supplementary=[dict(title="Whole-grid physical activity amount and first-moment magnitude",path=str((OUT/"stationary_amount_magnitude.png").relative_to(REPO)))])
    write("content_checkpoint.json",content);record(__file__);record(Path(prior.__file__));record(Path(native.__file__));write("source_sha256.json",SOURCE)
    print(json.dumps(dict(passed=True,conditions=80,new_integrations=0,maximum_prior_replay_error=error,threshold_records=len(events),output=str(OUT)),indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--render-only",action="store_true");args=parser.parse_args()
    if args.render_only:render(pd.read_csv(dest("stationary_whole_grid_summary.csv")))
    else:analyze()
