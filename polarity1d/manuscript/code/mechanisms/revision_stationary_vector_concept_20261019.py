"""Separate same-80 stationary-cue complex-vector concept sketch; no integrations."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

REPO = Path(__file__).resolve().parents[4]
SOURCE = REPO / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/standard_angles"
OUT = SOURCE / "vector_concept"
DATA = OUT / "supporting/data"
META = OUT / "supporting/metadata"
VECTORS = OUT / "supporting/vectors"
CASES = [("wave_pinning", "pi", 300., "Wave-pinning", "Retained reversal"),
         ("goryachev", "pi_2", 300., "Goryachev", "Cue excursion and return"),
         ("holmes_model3", "pi", 1000., "Holmes 3", "Two peaks, weak resultant"),
         ("spring", "pi_2", 10., "Spring ring", "Retained graded turn")]
ON = "#cf7a26"
OFF = "#367fb0"
MINIMUM = "#975ba5"
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    SOURCES[str(path)] = sha(path)
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def load_phase(path, spring, C):
    SOURCES[str(path)] = sha(path)
    with np.load(path) as z:
        t = z["time"].copy()
        state = z["states"].copy() if "states" in z else None
        if spring:
            n = (state.shape[1] - 1) // 2
            u, circumference, q = state[:, :n], state[:, -1], z["material_positions"].copy()
        else:
            u = state[:, 0] if state is not None else z["primary"].copy()
            q = z["x"].copy() / C
            circumference = np.full(len(t), C)
    return dict(time=t, u=u, C=circumference, q=q)


def load_cases():
    rows = read_csv(SOURCE / "supporting/data/standard_angles_summary.csv")
    assert len(rows) == 80
    SOURCES[str(SOURCE / "supporting/data/standard_angles_traces.npz")] = sha(SOURCE / "supporting/data/standard_angles_traces.npz")
    arrays = np.load(SOURCE / "supporting/data/standard_angles_traces.npz")
    result = []
    for model, angle, duration, name, outcome in CASES:
        row, = [r for r in rows if r["model"] == model and r["angle_label"] == angle and float(r["cue_duration"]) == duration]
        phases = [load_phase(row["actual_"+key+"_path"], model == "spring", float(row["circumference"])) for key in ("cue", "recovery")]
        for key in ("cue", "recovery"):
            assert SOURCES[row["actual_"+key+"_path"]] == row["actual_"+key+"_sha256"]
        a, b = phases
        np.testing.assert_array_equal(a["q"], b["q"])
        np.testing.assert_array_equal(a["u"][-1], b["u"][0])
        assert a["time"][-1] == duration and b["time"][-1] == 13000.
        q = a["q"]
        t = np.r_[a["time"], b["time"][1:] + duration]
        C = np.r_[a["C"], b["C"][1:]]
        u = np.vstack((a["u"], b["u"][1:]))
        M0 = C * u.mean(axis=1)
        M1 = C * (u @ np.exp(2j*np.pi*q)) / len(q)
        P = M1/M1[0]
        R = abs(P)
        cid = row["condition_id"]
        np.testing.assert_array_equal(t, arrays["c"+cid+"_time"])
        reference = arrays["c"+cid+"_P_real"] + 1j*arrays["c"+cid+"_P_imag"]
        discrepancy = float(abs(P-reference).max())
        assert discrepancy < 1e-12
        w = len(a["time"]) - 1
        minimum = int(R.argmin())
        initial_phase = np.rad2deg(np.angle(M1[0]))
        x = (360*q-initial_phase+90)%360-90
        order = np.argsort(x)
        result.append(dict(row=row, model=model, name=name, outcome=outcome, duration=duration,
            t=t, u=u, C=C, q=q, M0=M0, M1=M1, P=P, R=R, withdrawal=w, minimum=minimum,
            normal=x, order=order, compact_P_max_error=discrepancy))
    arrays.close()
    return result


def save(fig, name):
    paths = []
    for ext in ("png", "pdf", "svg"):
        path = (OUT if ext == "png" else VECTORS) / (name+"."+ext)
        fig.savefig(path, dpi=220, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def concept_sheet(cases):
    fig, axes = plt.subplots(3, 4, figsize=(14.2, 11.3), gridspec_kw=dict(height_ratios=(1.05, 1.15, 1)))
    maxP = max(float(abs(c["P"]).max()) for c in cases)
    limit = max(1.4, 1.12*maxP)
    for j, c in enumerate(cases):
        row=c["row"];u=c["u"];t=c["t"];w=c["withdrawal"];m=c["minimum"]
        ax=axes[0,j]
        target=float(row["target_angle_in_P_coordinates_deg"])
        width=float(row["initial_angular_support_deg"])
        ax.axvspan(target-width/2, target+width/2, color=ON, alpha=.1)
        for idx,color,style,label in ((0,".2",":","Initial"),(m,MINIMUM,"-.","Sampled min R"),(w,ON,"--","Withdrawal"),(-1,OFF,"-","Final")):
            ax.plot(c["normal"][c["order"]],u[idx,c["order"]],color=color,ls=style,lw=1.35,label=label)
        ax.set_xlim(-90,270);ax.set_ylim(bottom=0);ax.set_xticks((0,90,180,270))
        ax.set_xlabel("Outward-normal angle from initial M₁ (°)", fontsize=8)
        ax.set_ylabel("Native active concentration u", fontsize=8)
        ax.set_title(f"A{j+1}  {c['name']} · {float(row['requested_angle_deg']):g}°, T={c['duration']:g}\n{c['outcome']}",loc="left",fontsize=9.4,fontweight="bold")
        ax.text(.02,.97,f"Input peak {float(row['requested_peak']):g}; initial width {width:.2f}°\nmin-R profile t={t[m]:g}",transform=ax.transAxes,va="top",fontsize=6.8)
        ax=axes[1,j]
        ax.add_patch(Circle((0,0),1,fill=False,color=".78",ls=":",lw=.8))
        ax.add_patch(Circle((0,0),.1,facecolor=".92",edgecolor=".8",lw=.5))
        ax.axhline(0,color=".88",lw=.6);ax.axvline(0,color=".88",lw=.6)
        p=c["P"]
        ax.plot(p[:w+1].real,p[:w+1].imag,color=ON,lw=1.5)
        ax.plot(p[w:].real,p[w:].imag,color=OFF,lw=1.5)
        for idx,marker,color in ((0,"o",".2"),(w,"s",ON),(-1,"D",OFF)):
            ax.scatter(p[idx].real,p[idx].imag,s=31,marker=marker,facecolors="white",edgecolors=color,zorder=4,lw=1.3)
        ax.scatter(p[m].real,p[m].imag,marker="x",color=MINIMUM,s=26,zorder=5,lw=1)
        for segment,color in ((p[:w+1],ON),(p[w:],OFF)):
            # An arrow on each nontrivial phase, chosen by sampled path length.
            length=np.r_[0,np.cumsum(abs(np.diff(segment)))];total=length[-1]
            if total>1e-6:
                k=max(1,int(np.searchsorted(length,.55*total)))
                k=min(k,len(segment)-1)
                before=max(0,k-2)
                ax.annotate("",xy=(segment[k].real,segment[k].imag),xytext=(segment[before].real,segment[before].imag),arrowprops=dict(arrowstyle="->",color=color,lw=1.1),zorder=3)
        ax.set_aspect("equal",adjustable="box");ax.set_xlim(-limit,limit);ax.set_ylim(-limit,limit)
        ax.set_xlabel("Re P",fontsize=8);ax.set_ylabel("Im P",fontsize=8)
        ax.set_title(f"B{j+1}  Physical activity-vector trajectory",loc="left",fontsize=8.7,fontweight="bold")
        finaltheta="masked (R<.1)" if c["R"][-1]<.1 else f"{np.rad2deg(np.angle(p[-1])):.2f}°"
        if np.isclose(abs(np.rad2deg(np.angle(p[-1]))),180):finaltheta="180°"
        ax.text(.03,.03,f"Final heading: {finaltheta}",transform=ax.transAxes,fontsize=7)
        ax=axes[2,j]
        ax.axvspan(0,c["duration"],color=ON,alpha=.09)
        ax.axvline(c["duration"],color=ON,lw=.8,ls=":")
        ax.plot(t,c["R"],color="#6650a5",lw=1.4,label="R = |P|")
        ax.plot(t,c["M0"]/c["M0"][0],color="#398768",lw=1.2,ls="--",label="M₀ / M₀(0)")
        ax.axhline(.1,color=".65",lw=.6,ls=":")
        for y,color in ((c["R"],"#6650a5"),(c["M0"]/c["M0"][0],"#398768")):
            ax.scatter([t[w],t[-1]],[y[w],y[-1]],s=14,facecolors="white",edgecolors=color,zorder=3)
        ax.set_xscale("symlog",linthresh=c["duration"]/10)
        ax.set_xlim(0,t[-1]);ax.set_ylim(0,max(1.4,maxP*1.07))
        ax.set_xlabel(f"Native time (linear 0–{c['duration']/10:g}, then log)",fontsize=7.6)
        ax.set_ylabel("Initial-normalized physical quantity",fontsize=8)
        ax.set_title(f"C{j+1}  Resultant versus total activity",loc="left",fontsize=8.7,fontweight="bold")
        if j==0:ax.legend(frameon=False,fontsize=7,loc="lower left")
    fig.suptitle("Stationary cues · spatial activity builds a complex vector",fontsize=16,y=.995)
    fig.text(.5,.951,r"$M_1(t)=\sum_j u_j(t)\,\Delta s_j(t)\,e^{i\theta_j},\qquad \Delta s_j(t)=C(t)/N,\qquad P(t)=M_1(t)/M_1(0)$",ha="center",fontsize=15)
    fig.text(.5,.918,"Each contribution is active concentration × physical segment length × outward normal; opposite contributions can cancel.  M₀ = ΣuⱼΔsⱼ.",ha="center",fontsize=9)
    fig.text(.5,.890,"Four actual conditions from the same 80-case experiment · cue fixed in place, then 13,000 cue-free native units · native amplitudes and clocks differ",ha="center",fontsize=8.5)
    handles=[Line2D([],[],color=color,ls=style,label=label) for color,style,label in ((".2",":","Initial profile"),(MINIMUM,"-.","Sampled minimum-R profile"),(ON,"--","Withdrawal profile"),(OFF,"-","Final profile"),(ON,"-","Cue on"),(OFF,"-","Cue released"))]
    handles += [Line2D([],[],color=color,marker=marker,ls="none",mfc="white",label=label) for marker,color,label in (("o",".2","Initial P"),("s",ON,"Withdrawal P"),("D",OFF,"Final P"))]
    fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.5,-.008),ncol=5,fontsize=7.5,frameon=False)
    fig.text(.5,.041,"P planes have equal real/imaginary scales and a unit-circle reference. Shaded R<.1 disk: heading is not interpreted. Profile shading: initial cue footprint.",ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.08,1,.905),h_pad=2.0,w_pad=2)
    return save(fig,"stationary_vector_concept")


def native_changes(cases):
    fig,axes=plt.subplots(1,4,figsize=(14,5.4))
    for j,c in enumerate(cases):
        ax=axes[j];delta=c["u"]-c["u"][0]
        vmax=float(abs(delta).max())
        mesh=ax.pcolormesh(c["normal"][c["order"]],c["t"],delta[:,c["order"]],cmap="RdBu_r",vmin=-vmax,vmax=vmax,shading="nearest",rasterized=True)
        ax.set_yscale("symlog",linthresh=c["duration"]/10);ax.set_ylim(c["t"][-1],0)
        ax.axhline(c["duration"],color="black",ls=":",lw=.9)
        ax.set_xlim(-90,270);ax.set_xticks((0,90,180,270));ax.set_xlabel("Normal angle from initial M₁ (°)",fontsize=8)
        ax.set_ylabel("Native time (linear near 0; then log)",fontsize=8)
        ax.set_title(f"{c['name']}\n{c['outcome']}",loc="left",fontsize=9.3)
        cb=fig.colorbar(mesh,ax=ax,orientation="horizontal",fraction=.05,pad=.19)
        cb.set_label("Δu from initial (native concentration)",fontsize=7.3);cb.ax.tick_params(labelsize=7)
    fig.suptitle("Actual spatial changes through cue and release",fontsize=15)
    fig.text(.5,.912,"Full native histories of the same four cases · dashed line: cue withdrawal · signed concentration changes, with no peak normalization",ha="center",fontsize=8.5)
    fig.tight_layout(rect=(0,.02,1,.87),w_pad=2)
    return save(fig,"stationary_vector_native_changes")


def records(cases):
    rows=[];traces={};profiles=[]
    for c in cases:
        w=c["withdrawal"];m=c["minimum"];p=c["P"];ratio=c["M0"]/c["M0"][0]
        row=dict(model=c["model"],condition_id=c["row"]["condition_id"],cue_angle_degrees=float(c["row"]["requested_angle_deg"]),cue_duration=c["duration"],
            initial_M0=float(c["M0"][0]),initial_M1_real=float(c["M1"][0].real),initial_M1_imag=float(c["M1"][0].imag),
            minimum_sampled_R=float(c["R"][m]),minimum_R_native_time=float(c["t"][m]),M0_ratio_at_minimum_R=float(ratio[m]),
            withdrawal_R=float(c["R"][w]),withdrawal_P_real=float(p[w].real),withdrawal_P_imag=float(p[w].imag),withdrawal_M0_ratio=float(ratio[w]),
            final_R=float(c["R"][-1]),final_P_real=float(p[-1].real),final_P_imag=float(p[-1].imag),final_M0_ratio=float(ratio[-1]),
            final_angle_degrees=float(np.rad2deg(np.angle(p[-1]))),final_heading_defined=bool(c["R"][-1]>=.1),
            circumference_min=float(c["C"].min()),circumference_max=float(c["C"].max()),final_resolved_peaks=int(c["row"]["final_resolved_peaks"]),
            compact_P_maximum_difference=c["compact_P_max_error"],native_cue_archive=c["row"]["actual_cue_path"],native_recovery_archive=c["row"]["actual_recovery_path"])
        rows.append(row)
        for k,a in (("time",c["t"]),("M0",c["M0"]),("M1_real",c["M1"].real),("M1_imag",c["M1"].imag),("P_real",p.real),("P_imag",p.imag),("C",c["C"])):
            traces[c["model"]+"_"+k]=a
        for label,index in (("initial",0),("sampled_minimum_R",m),("withdrawal",w),("final",len(c["t"])-1)):
            for cell in range(len(c["q"])):
                profiles.append(dict(model=c["model"],snapshot=label,time=float(c["t"][index]),material_position=float(c["q"][cell]),normal_degrees=float(c["normal"][cell]),physical_ds=float(c["C"][index]/len(c["q"])),native_active=float(c["u"][index,cell])))
    for filename,items in (("selected_case_measurements.csv",rows),("selected_native_profiles.csv",profiles)):
        with (DATA/filename).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(items[0]));writer.writeheader();writer.writerows(items)
    np.savez_compressed(DATA/"selected_physical_moment_traces.npz",**traces)
    return rows


def write_docs(rows,artifacts):
    caption="Stationary-cue complex-vector concept sketch from four selected conditions in the completed 80-condition experiment. A1–A4 show actual native active-concentration profiles at the initial time, the sampled minimum R, cue withdrawal and the final 13,000-unit cue-free endpoint. The minimum-R profile is taken from the exact same stored native frame as its moment. Shading marks the initial cue footprint. B1–B4 show measured P=M₁(t)/M₁(0) in the complex plane, with orange cue-on and blue post-withdrawal paths; circle/square/diamond mark initial/withdrawal/final values. Dividing by the initial complex M₁ rotates the initial direction to the positive real axis and scales by its initial magnitude; it does not divide by the current total active amount. Thus R can exceed 1, as during the Spring cue, or become small through cancellation of opposed contributions while M₀ remains substantial, as in Holmes. All planes have the same equal real/imaginary scales, a unit-circle reference and an R<.1 disk where heading is not interpreted. C1–C4 retain full native-time traces of R and physical total active amount M₀/M₀(0), with explicitly labelled logarithmic compression beyond the early linear interval. Physical M₁ and M₀ use C(t)/N quadrature, including Spring's changing circumference. Native concentrations, input peaks, widths and clocks differ across models. This review sheet does not replace existing Figure 3 art or characterize moving cues."
    caption += " Spring's sampled minimum R differs from 1 only at numerical roundoff, so that snapshot is bookkeeping rather than a weak-resultant episode. For its 90° input, 35.6% of the initial cue support overlaps the original half-contrast core; that case does not begin with two separate activity sites."
    (META/"CAPTION.md").write_text("# Stationary-cue concept caption\n\n"+caption+"\n")
    text=r"""# Spatial activity and the physical complex vector

This separate reviewable sketch uses the same four suggested cases from the completed 80-condition stationary-cue checkpoint. It adds no native integration and leaves all previous art, content and outline files intact.

- [Complex-vector concept sheet](stationary_vector_concept.png) · [PDF](supporting/vectors/stationary_vector_concept.pdf) · [SVG](supporting/vectors/stationary_vector_concept.svg)
- [Actual signed spatial changes](stationary_vector_native_changes.png)
- [Selected measurements](supporting/data/selected_case_measurements.csv) · [Native profile snapshots](supporting/data/selected_native_profiles.csv)
- [Caption](supporting/metadata/CAPTION.md) · [Source and validation record](supporting/metadata/validation.json)

For outward normal angle θⱼ and native active concentration uⱼ, the physical ring quadrature is

$$
M_1(t)=\sum_j u_j(t)\,\Delta s_j(t)e^{i\theta_j},\qquad
\Delta s_j(t)=\frac{C(t)}{N},\qquad
P(t)=\frac{M_1(t)}{M_1(0)}.
\tag{1}
$$

$$
M_0(t)=\sum_j u_j(t)\,\Delta s_j(t),\qquad R(t)=|P(t)|.
\tag{2}
$$

Each complex contribution combines activity, physical segment length and outward-normal direction. Opposed contributions can cancel despite substantial total activity. P is normalized by the initial complex physical moment: division rotates the initial direction to the positive real axis and scales by the initial magnitude, giving P(0)=1. It does not normalize by the total active amount at each frame, so R may exceed 1. The physical M₁/M₀ shape-normalized helper is not substituted. Spring uses its actual circumference at every sample. No profile is divided by an independent history maximum.

## What the actual four cases show

- Wave-pinning, 180°, T=300: the measured real-axis path changes sign near the origin and retains a reversed final full-strength vector. The minimum sampled magnitude is about .00262 at t=205; this sampling does not assert an exactly sampled zero.
- Goryachev, 90°, T=300: the cue produces a substantial excursion to 33.48° with R=.721 at withdrawal, followed by return to −.0278° and R≈1. Its raw active peak changes from about 64.95 to 37.31 at withdrawal. The near-original-axis endpoint does not imply an absent response.
- Holmes 3, 180°, T=1000: the final two-peak profile has R=.03993 while physical total active amount is .99746 of initial. Direction is masked, and loss of the resultant does not mean loss of nearly all activity.
- Spring, 90°, T=10: the cue produces 76.09° with R=1.281 at withdrawal and retains 72.44° with R≈1 after release. The final native profile is one shifted full-strength activity region. Circumference is included in both physical moments. Its sampled minimum R differs from 1 only at numerical roundoff, rather than marking a weak-resultant episode. About 35.6% of the initial cue support overlaps the original half-contrast core, so the case does not start with two separate activity sites.

The antipodal chemical paths have imaginary components at numerical roundoff; equal complex-plane scales prevent numerical asymmetry from being magnified. The R<.1 disk concerns heading interpretation, not validity of the measured weak complex vector. The displayed profile time at minimum R is the exact stored native frame, with no nearest-profile substitution. The signed Δu history companion exposes spatial changes in native concentration units; its colour limits are explicitly numerical rather than a peak normalization.

The physical total active amount M₀ is distinct from total conserved active+inactive pool. Near-constant M₀ with weak M₁ supports cancellation or redistribution of activity rather than global activity disappearance; it does not assign a particular site lineage or functional competition. Shared cue-centre angles do not match native cue dose, footprint or timescale. A finite 13,000-unit withdrawal does not prove asymptotic stability. Activity orientation is not a measured movement vector.

## Source and validation

The selected native challenge/recovery archives are identified and hashed in the measurement table and source JSON. Physical P is recomputed directly from every selected native profile and matches the retained standard-angle trace. All four suggested cases were retained after inspecting their actual response. The complete 80-condition experiment remains in the parent standard_angles directory.
"""
    (OUT/"README.md").write_text(text)
    validation=dict(generated_at=datetime.now(timezone.utc).isoformat(),new_native_integrations=0,selected_cases=rows,
        cases_retained_as_requested=True,physical_quadrature="C(t)/N",M1_normalization="initial complex physical M1",profiles_peak_normalized=False,
        sources=[dict(path=p,sha256=h) for p,h in SOURCES.items()],source_script=str(Path(__file__)),source_script_sha256=sha(__file__),
        artifacts=[dict(path=str(p),sha256=sha(p)) for p in artifacts],visual_inspection="pending",independent_review="pending")
    (META/"validation.json").write_text(json.dumps(validation,indent=2))


def main():
    for folder in (OUT,DATA,META,VECTORS):folder.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.spines.top":False,"axes.spines.right":False,
        "axes.linewidth":.65,"xtick.labelsize":7,"ytick.labelsize":7,"pdf.fonttype":42,"svg.fonttype":"none"})
    cases=load_cases();rows=records(cases);artifacts=concept_sheet(cases);artifacts+=native_changes(cases);write_docs(rows,artifacts)
    print(json.dumps(dict(output=str(OUT),cases=4,artifacts=len(artifacts))))


if __name__=="__main__":main()
