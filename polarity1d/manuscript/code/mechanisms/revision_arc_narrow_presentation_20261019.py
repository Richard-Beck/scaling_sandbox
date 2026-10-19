"""Present the narrow-range native cue-on experiment from atomic readouts only."""
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
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
NAMES=dict(wave_pinning="Wave-pinning",goryachev="Goryachev",holmes_model3="Holmes",spring="Spring")


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clean(values):return [float(v) if np.isfinite(v) else None for v in values]


def main():
    data=OUT/"supporting/data";metadata=OUT/"supporting/metadata";vectors=OUT/"supporting/vectors";vectors.mkdir(parents=True,exist_ok=True)
    status=json.loads((metadata/"latest_status.json").read_text())
    paths=[data/"narrow_arc_summary.csv",data/"narrow_arc_traces.npz",data/"narrow_arc_cue_on_RHS_rates.csv"]
    source={str(p):sha(p) for p in paths}
    for p in paths:assert source[str(p)]==status["artifact_sha256"][str(p.relative_to(OUT))]
    table=pd.read_csv(paths[0]);rhs=pd.read_csv(paths[2]);cases=table[table.cue_duration==1000.].copy();assert len(table)==645 and len(cases)==129
    curves=[]
    with np.load(paths[1]) as traces:
        for r in cases.to_dict("records"):
            cid=int(r["condition_id"]);complete=r["status"]=="complete"
            item=dict(condition_id=cid,model=r["model"],size=float(r["size_parameter"]),status=r["status"],
                      initial_arc=float(r["physical_distance_at_onset"]),target_deg=float(r["target_angle_in_P_coordinates_deg"]),
                      initial_C=float(r["circumference"]),full_physical_support=float(r["full_physical_support"]),
                      centre50=r["centre_relative_to_core_50"],centre10=r["centre_relative_to_core_10"],
                      footprint50=r["footprint_relative_to_core_50"],footprint10=r["footprint_relative_to_core_10"],
                      edge50=float(r["initial_edge_right_50"]),edge10=float(r["initial_edge_right_10"]),
                      signed_edge50=float(r["signed_centre_distance_from_positive_edge_50"]),signed_edge10=float(r["signed_centre_distance_from_positive_edge_10"]),
                      near_support_edge_gap10=float(r["signed_cue_near_edge_gap_from_positive_core_edge_10"]),
                      source_kind=r["source_kind"])
            if complete:
                mask=traces[f"c{cid}_phase"]==0
                phase=traces[f"c{cid}_angle_masked_deg"][mask].copy();phase[np.isclose(abs(phase),180,atol=1e-6,rtol=0)]=180.
                item.update(time=clean(traces[f"c{cid}_time"][mask]),heading=clean(phase),R=clean(traces[f"c{cid}_R"][mask]),
                            closure=clean(traces[f"c{cid}_angular_target_error_closed"][mask]),C=clean(traces[f"c{cid}_C"][mask]),
                            prescribed_arc=clean(traces[f"c{cid}_prescribed_physical_offset"][mask]),
                            M0_ratio=clean(traces[f"c{cid}_M0_ratio"][mask]),regional_excess=clean(traces[f"c{cid}_fixed_region_transported_background_excess"][mask]),
                            rate_left=clean(traces[f"c{cid}_cue_rate_left"]),rate_right=clean(traces[f"c{cid}_cue_rate_right"]),
                            rate=clean(traces[f"c{cid}_cue_sampled_angular_secant"]),
                            weak_frames=int(r["cue_weak_saved_frames"]),minimum_R=float(r["R_cue_on_sampled_minimum"]),
                            full_state_coverage=r["cue_archive_full_state_coverage"],field_tail=r["full_field_tail_status"],
                            actual_cue_path=r["actual_cue_path"],actual_cue_sha256=r["actual_cue_sha256"])
                anchors=rhs[rhs.condition_id==cid].to_dict("records")
                item["RHS_anchors"]=[{key:(None if pd.isna(value) else value) for key,value in a.items()} for a in anchors]
            else:
                item.update(time=[],heading=[],R=[],closure=[],C=[],prescribed_arc=[],M0_ratio=[],regional_excess=[],rate_left=[],rate_right=[],rate=[],RHS_anchors=[])
            curves.append(item)
    p=data/"narrow_cue_on129_interactive.json";temp=p.with_suffix(".tmp.json");temp.write_text(json.dumps(curves,separators=(",",":"),allow_nan=False));temp.replace(p)
    render(curves,status)
    report=dict(generated_at=datetime.now(timezone.utc).isoformat(),requested_conditions=645,available_conditions=status["complete"],
                requested_longest_cue_histories=129,available_longest_cue_histories=sum(c["status"]=="complete" for c in curves),
                status="complete_native_readout_coverage" if status["complete"]==645 else "explicit_partial_snapshot",source_sha256=source,
                collector_generation=status["generated_at"],new_integrations=0,continuous_native_RHS_scan=False,
                method="Saved physical activity-vector curves and actual-interval angular secants. Weak R<0.1 is a gap; no unwrapping. Initial spatial classes remain readout definitions, not pathway or ancestry assignments.",
                displayed_static_scope="Only shared distances whose full footprints clear the initial10%-contrast region at all three sizes. Other regimes remain in full interactive data.",
                renderer_sha256=sha(__file__),interactive_sha256=sha(p))
    (metadata/"presentation.json").write_text(json.dumps(report,indent=2)+"\n")
    # The collector may see the preceding JSON generation before this presenter
    # replaces it. Synchronize this display-only artifact hash after publishing;
    # scientific CSV/NPZ hashes and collector-generation time remain unchanged.
    status["artifact_sha256"]["supporting/data/narrow_cue_on129_interactive.json"]=sha(p)
    status["presentation_metadata_updated_at"]=report["generated_at"]
    q=metadata/"latest_status.tmp.json";q.write_text(json.dumps(status,indent=2)+"\n");q.replace(metadata/"latest_status.json")
    print(json.dumps(report,indent=2))


def render(curves,status):
    plt.rcParams.update({"font.size":10,"pdf.fonttype":42,"svg.fonttype":"none"})
    for model,name in NAMES.items():
        selected=[c for c in curves if c["model"]==model];sizes=sorted({c["size"] for c in selected})
        offsets=sorted(d for d in {c["initial_arc"] for c in selected} if all(c["footprint10"]=="outside" for c in selected if c["initial_arc"]==d))
        colours={d:plt.get_cmap("viridis")(.08+.84*j/max(1,len(offsets)-1)) for j,d in enumerate(offsets)}
        fig,axes=plt.subplots(3,3,figsize=(14,11.5));fig.subplots_adjust(left=.075,right=.985,top=.79,bottom=.145,hspace=.34,wspace=.24)
        available=sum(c["status"]=="complete" for c in selected)
        fig.text(.04,.955,f"{name} · Cue-on response beyond the initial 10% edge",fontsize=19,weight="bold")
        fig.text(.04,.920,f"Narrow-range experiment · T = 1000 native units · {available}/{len(selected)} longest-cue histories available in this snapshot",fontsize=11)
        fig.legend(handles=[Line2D([0],[0],color=colours[d],lw=1.8,label=f"d={d:.5g}") for d in offsets],loc="upper left",bbox_to_anchor=(.04,.90),ncol=5,frameon=False,fontsize=10)
        model_rates=[v for c in selected if c["status"]=="complete" and c["initial_arc"] in offsets for v in c["rate"] if v is not None]
        rate_bounds=(min(-.001,min(model_rates,default=0)*1.1),max(.001,max(model_rates,default=0)*1.1))
        for col,size in enumerate(sizes):
            panel=[c for c in selected if c["size"]==size and c["initial_arc"] in offsets]
            complete=[c for c in panel if c["status"]=="complete"]
            for c in complete:
                color=colours[c["initial_arc"]]
                axes[0,col].plot(c["time"],c["heading"],color=color,lw=1.5)
                axes[0,col].axhline(c["target_deg"],color=color,ls=":",lw=.6,alpha=.65)
                axes[1,col].plot(c["time"],c["R"],color=color,lw=1.5)
                segments=[[(left,rate),(right,rate)] for left,right,rate in zip(c["rate_left"],c["rate_right"],c["rate"]) if rate is not None]
                if segments:
                    axes[2,col].add_collection(LineCollection(segments,colors=[color],linewidths=.8))
                    axes[2,col].autoscale_view(scalex=False,scaley=True)
            C=next(c["initial_C"] for c in selected if c["size"]==size)
            title=f"C={size:g}" if model!="spring" else f"C₀={size:g}; initial C={C:.5f}"
            axes[0,col].set_title(title+f"\n{len(complete)}/{len(panel)} outside histories available",loc="left",fontsize=11)
            if not complete:
                for ax in axes[:,col]:ax.text(.5,.5,"Native history pending",ha="center",va="center",transform=ax.transAxes,color="#697781")
            axes[0,col].set_ylim(-4,184);axes[1,col].axhline(.1,color="#999999",ls="--",lw=.7);axes[1,col].set_ylim(0,2)
            for ax in axes[:,col]:
                ax.set_xscale("symlog",linthresh=.1);ax.set_xlim(0,1000);ax.set_xticks([0,.1,1,10,100,1000],["0",".1","1","10","100","1000"]);ax.minorticks_off();ax.spines[["top","right"]].set_visible(False)
            axes[2,col].set_yscale("symlog",linthresh=.001);axes[2,col].axhline(0,color="#bbbbbb",lw=.6)
            axes[2,col].set_ylim(*rate_bounds)
            axes[2,col].set_xlabel("Time since cue onset (native units)")
        axes[0,0].set_ylabel("Activity-vector heading (degrees)");axes[1,0].set_ylabel("R = |P|");axes[2,0].set_ylabel("Sampled angular secant\n(degrees / native unit)")
        fig.text(.04,.092,"Shared distances only: full supports clear all three initial 10%-contrast regions. Dotted heading lines mark actual cue targets.",fontsize=10)
        fig.text(.04,.068,"Rate segments span their actual saved intervals, not instantaneous native-RHS values. Weak R < 0.1 creates gaps; no phase is unwrapped across them.",fontsize=9.5)
        fig.text(.04,.044,"Time is linear 0–0.1, logarithmic thereafter. Rate y axes are linear near zero and logarithmic in magnitude beyond ±0.001; scales are not cross-model clocks.",fontsize=9.2)
        fig.text(.04,.022,"A global angular rate includes projection and cancellation; it is not local-peak travel. Near-antipodal targets require the accompanying R and native profiles.",fontsize=9.5)
        for ext in ("png","pdf","svg"):
            path=OUT/("narrow_external_"+model+".png") if ext=="png" else OUT/"supporting/vectors"/("narrow_external_"+model+"."+ext)
            fig.savefig(path,dpi=180,facecolor="white")
        plt.close(fig)
        print(OUT/("narrow_external_"+model+".png"),flush=True)


if __name__=="__main__":main()
