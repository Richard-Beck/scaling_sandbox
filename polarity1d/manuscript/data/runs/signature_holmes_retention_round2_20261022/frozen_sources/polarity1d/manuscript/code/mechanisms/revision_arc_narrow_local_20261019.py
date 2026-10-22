"""One clean external native size contrast from the complete stationary experiment."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT=Path(__file__).resolve().parents[2]
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
OUT=TOPIC/"local_size_contrast"


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    for leaf in ("data","metadata","vectors"):(OUT/"supporting"/leaf).mkdir(parents=True,exist_ok=True)
    p=TOPIC/"supporting/data/narrow_arc_summary.csv";q=TOPIC/"cue_on_core/supporting/data/cue_on_core_native_rates.csv"
    summary=pd.read_csv(p);rates=pd.read_csv(q,float_precision="round_trip");rows=summary[summary.condition_id.isin([39,149])].sort_values("size_parameter")
    assert len(rows)==2 and (rows.model=="wave_pinning").all() and (rows.physical_distance_at_onset==11).all() and (rows.cue_duration==1000).all()
    assert (rows.footprint_relative_to_core_10=="outside").all() and (rows.full_field_tail_status=="passes_finite_tail_gate").all()
    data=[];export=[];source={str(p):sha(p),str(q):sha(q)}
    for r in rows.to_dict("records"):
        archive=Path(r["actual_cue_path"]);assert sha(archive)==r["actual_cue_sha256"];source[str(archive)]=sha(archive)
        with np.load(archive) as z:
            u=z["states"][:,0].copy();time=z["time"].copy();position=z["x"].copy()
        relative=(position-r["preparation_axis_fraction"]*r["circumference"]+r["circumference"]/2)%r["circumference"]-r["circumference"]/2
        order=np.argsort(relative);relative=relative[order];u=u[:,order]
        native=rates[rates.condition_id==r["condition_id"]].sort_values("frame")
        np.testing.assert_array_equal(time,native.cue_time)
        data.append(dict(row=r,u=u,time=time,position=relative,native=native))
        for j,t in enumerate(time):
            export.extend(dict(condition_id=r["condition_id"],size_parameter=r["size_parameter"],cue_time=float(t),initial_physical_arc=float(s),native_active=float(v)) for s,v in zip(relative,u[j]))
    pd.DataFrame(export).to_csv(OUT/"supporting/data/native_cue_on_fields.csv",index=False)
    rows.to_csv(OUT/"supporting/data/selected_conditions.csv",index=False)
    rates[rates.condition_id.isin(rows.condition_id)].to_csv(OUT/"supporting/data/complex_P_and_native_rate.csv",index=False)
    plt.rcParams.update({"font.size":10,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(14,6.6));grid=fig.add_gridspec(1,3,left=.065,right=.97,top=.76,bottom=.23,wspace=.40,width_ratios=[1,1,1.08])
    axes=[fig.add_subplot(grid[0,j]) for j in range(3)];cmap=plt.get_cmap("viridis");norm=Normalize(0,max(d["u"].max() for d in data));colors=["#c98a20","#564d86"]
    fig.text(.04,.958,"The same external physical cue, two native incumbent sizes",fontsize=18,weight="bold")
    fig.text(.04,.907,"Wave-pinning · d₀ = 11 · full support 3 · coefficient peak 0.1 · held cue T = 1000 · native kinetics and exact basals retained",fontsize=10.5)
    for ax,d,color in zip(axes[:2],data,colors):
        r=d["row"];image=ax.pcolormesh(d["position"],d["time"],d["u"],cmap=cmap,norm=norm,shading="nearest",rasterized=True)
        ax.set_yscale("symlog",linthresh=.1);ax.set_ylim(0,1000);ax.set_yticks([0,.1,1,10,100,1000],["0",".1","1","10","100","1000"])
        ax.set_xlim(-20,20)
        if r["circumference"]<40:
            ax.axvspan(-20,-r["circumference"]/2,color="#e7eaec",zorder=0)
            ax.axvspan(r["circumference"]/2,20,color="#e7eaec",zorder=0)
        ax.axvline(11,color="white",ls=":",lw=1.1)
        ax.axvline(r["initial_edge_left_10"],color="#f0b875",ls="--",lw=.8);ax.axvline(r["initial_edge_right_10"],color="#f0b875",ls="--",lw=.8)
        ax.plot([9.5,12.5],[0,0],color="#cf852b",lw=5,solid_capstyle="butt",clip_on=False)
        ax.set_title(f"C={r['size_parameter']:g}; actual target {r['target_angle_in_P_coordinates_deg']:.0f}°",loc="left",fontsize=11,color=color)
        ax.set_xlabel("Physical arc from initial preparation axis",fontsize=9)
    axes[0].set_ylabel("Cue-on time (native units)",fontsize=10)
    for d,color in zip(data,colors):
        r,n=d["row"],d["native"]
        axes[2].plot(n.P_real,n.P_imag,color=color,lw=1.8,label=f"C={r['size_parameter']:g}")
        axes[2].scatter(n.P_real.iloc[-1],n.P_imag.iloc[-1],s=38,color=color,marker="s",zorder=4)
        target=np.deg2rad(r["target_angle_in_P_coordinates_deg"])
        axes[2].plot([0,np.cos(target)],[0,np.sin(target)],color=color,lw=.7,ls=":")
        for wanted in [10,100,300]:
            j=int(np.argmin(abs(n.cue_time.to_numpy()-wanted)));a=n.iloc[j]
            axes[2].scatter(a.P_real,a.P_imag,s=16,color=color,zorder=4)
            axes[2].annotate(f"{a.cue_time:.0f}",(a.P_real,a.P_imag),xytext=(3,4 if color==colors[0] else -10),textcoords="offset points",fontsize=7.3,color=color)
    angle=np.linspace(0,2*np.pi,300);axes[2].plot(np.cos(angle),np.sin(angle),color="#aeb9bf",ls="--",lw=.7)
    axes[2].axhline(0,color="#ccd3d7",lw=.6);axes[2].axvline(0,color="#ccd3d7",lw=.6);axes[2].scatter([1],[0],s=24,color="#253a44",zorder=5)
    axes[2].text(.77,-.16,"initial P = 1",fontsize=8);axes[2].set_xlim(-1.15,1.25);axes[2].set_ylim(-1.2,1.2);axes[2].set_aspect("equal",adjustable="box")
    axes[2].set_xlabel("Re P",fontsize=10);axes[2].set_ylabel("Im P",fontsize=10);axes[2].set_title("The actual global complex P path",loc="left",fontsize=11)
    axes[2].legend(loc="upper right",frameon=False,fontsize=9)
    for ax in axes:ax.spines[["top","right"]].set_visible(False);ax.tick_params(labelsize=8.5)
    cax=fig.add_axes((.26,.15,.26,.018));bar=fig.colorbar(image,cax=cax,orientation="horizontal");bar.set_label("Native active concentration (shared scale)",fontsize=8.5,labelpad=2);bar.ax.tick_params(labelsize=8)
    fig.text(.04,.068,"Every saved native cue-on frame is shown. Dashed: initial 10%-contrast edges. White/orange: cue centre/support. Grey padding is outside the shorter domain.",fontsize=8.8)
    fig.text(.04,.042,"Time is linear to 0.1 then log. Complex axes have equal scale; dashed circle is R = 1; squares mark withdrawal and numbers give actual native times.",fontsize=9)
    fig.text(.04,.016,"Physical cue arc is matched; native target angle and incumbent architecture differ. These are activity patterns and resultants, not site ancestry or cell movement.",fontsize=9)
    for ext in ("png","pdf","svg"):
        path=OUT/"external_local_size_contrast.png" if ext=="png" else OUT/"supporting/vectors"/("external_local_size_contrast."+ext)
        fig.savefig(path,dpi=180,facecolor="white")
    plt.close(fig)
    stamp=datetime.now(timezone.utc).isoformat()
    (OUT/"README.md").write_text(f"""# One native external-cue size contrast

Generated {stamp}. [PNG](external_local_size_contrast.png) · [PDF](supporting/vectors/external_local_size_contrast.pdf) · [SVG](supporting/vectors/external_local_size_contrast.svg) · [every native field sample](supporting/data/native_cue_on_fields.csv) · [complexP and exact native rates](supporting/data/complex_P_and_native_rate.csv) · [all case metadata](supporting/data/selected_conditions.csv).

This selected illustration fixes model(Wave-pinning), initial physical cue distance11, native fullsupport3, coefficientpeak0.1 and duration1000, and changes native basal circumference30→40. Both whole cue footprints are outside the initial10%-contrast regions, both branches retain full native fields, neither cue-on vector reachesR<0.1, and both later finite field tails pass. These are new645-grid conditionIDs39/149; they are selected for a clean geometric size contrast, not retained success. Targetangles are132°/99°. Initial firstmoments and incumbent shapes change together, so this is not a causal partition of circumference mapping versus incumbent architecture.

The two native kymographs retain every saved cue-on frame, actual cell locations, common zero-based native concentration limits and an explicitly compressed native time axis(linear0–0.1, logarithmic afterward). Initial10%-contrast edges are operational geometric readouts. The orange bar marks full cue support; periodic domain coordinates are relative to the declared initial preparation axis. The complexP plane uses exact normalized physical moments in initial-aligned Cartesian axes with equal scales. Squares are withdrawal; labelled times are actual nearby saved times, not interpolated field states. The R=1 reference circle distinguishes direction from moment magnitude. No site lineage, peak speed or cell motion is inferred. Full distance/size/rate coverage and all other durations remain in the comprehensive supplement.
""")
    result=dict(generated_at=stamp,passed=True,condition_ids=[39,149],native_field_samples=len(export),new_integrations=0,source_sha256=source,
                renderer_sha256=sha(__file__),artifact_sha256={str(x.relative_to(OUT)):sha(x) for x in OUT.rglob("*") if x.is_file() and x.name!="verification.json"})
    (OUT/"supporting/metadata/verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(OUT/"external_local_size_contrast.png")


if __name__=="__main__":main()
