"""Systematic distance/size view of existing exact native RHS rate anchors."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
OUT=TOPIC/"rate_map"
NAMES=dict(wave_pinning="Wave-pinning",goryachev="Goryachev",holmes_model3="Holmes",spring="Spring")
MARKERS=dict(inside="o",straddling="D",outside="^")
COLORS=["#c98a20","#3c8894","#564d86"]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def draw(rows,time,tag):
    plt.rcParams.update({"font.size":11,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig,axes=plt.subplots(2,2,figsize=(13,10.5))
    fig.subplots_adjust(left=.08,right=.98,top=.83,bottom=.17,hspace=.48,wspace=.27)
    fig.text(.045,.965,"Global reorientation rate across physical cue distance and size",fontsize=18,weight="bold")
    fig.text(.045,.930,f"Exact native right-hand side at cue-on time t={time:g} · 129 prescribed placements · all initial footprint classes retained",fontsize=11)
    handles=[Line2D([0],[0],marker=m,color="#53646e",lw=0,markersize=6,label="Initial full support "+r) for r,m in MARKERS.items()]
    fig.legend(handles=handles,loc="upper left",bbox_to_anchor=(.04,.912),ncol=3,frameon=False,fontsize=10)
    for ax,(model,name) in zip(axes.flat,NAMES.items()):
        part=rows[rows.model==model];legend=[]
        for color,size in zip(COLORS,sorted(part.size_parameter.unique())):
            r=part[part.size_parameter==size].sort_values("physical_distance_at_onset")
            x=r.physical_distance_at_onset.to_numpy();y=r[f"omega_t{int(time)}"].to_numpy()
            ax.plot(x,y,color=color,lw=1.0,alpha=.85)
            for regime,marker in MARKERS.items():
                keep=r.footprint_relative_to_core_10==regime
                ax.scatter(x[keep],y[keep],color=color,marker=marker,s=35,edgecolors="white",linewidths=.45,zorder=3)
            label=("C₀=" if model=="spring" else "C=")+f"{size:g}"
            legend.append(Line2D([0],[0],color=color,lw=1.5,label=label))
        ax.legend(handles=legend,fontsize=9,frameon=False,loc="upper left")
        ax.set_title(name,loc="left",weight="bold",fontsize=13)
        ax.set_xlabel("Prescribed initial physical centre arc d₀",fontsize=10)
        ax.set_ylabel("Angular rate (degrees / native time)",fontsize=10)
        ax.axhline(0,color="#b9c2c7",lw=.7);ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(labelsize=9)
        if model=="spring":ax.set_xticks([0,.25,.5,.75,1,1.25])
    fig.text(.045,.112,"Markers classify the complete initial cue support against the interpolated 10%-contrast edge; centre membership is separately reported.",fontsize=9.5)
    fig.text(.045,.088,"All 129 rates are observed with R ≥ 0.1 at this instant. Lines guide the eye between tested positions; no continuous-distance fit is asserted.",fontsize=9.5)
    fig.text(.045,.064,"The rate is activity-vector phase change, not local-peak speed. Initial projection, current moment magnitude and native organization all contribute.",fontsize=9.4)
    fig.text(.045,.040,"Same native coefficients/supports/basals within each model. Holding d₀ fixed changes angular target and basal architecture; model clocks are uncalibrated.",fontsize=9.1)
    for ext in ("png","pdf","svg"):
        path=OUT/(tag+".png") if ext=="png" else OUT/"supporting/vectors"/(tag+"."+ext)
        fig.savefig(path,dpi=180,facecolor="white")
    plt.close(fig)


def main():
    for leaf in ("data","metadata","vectors"):(OUT/"supporting"/leaf).mkdir(parents=True,exist_ok=True)
    p=TOPIC/"supporting/data/narrow_arc_summary.csv";q=TOPIC/"supporting/data/narrow_arc_cue_on_RHS_rates.csv"
    table=pd.read_csv(p);rhs=pd.read_csv(q)
    selected=table[table.cue_duration==1000].copy();assert len(selected)==129 and (selected.status=="complete").all()
    for t in (0,10):
        rates=rhs[(rhs.requested_time==t)&rhs.condition_id.isin(selected.condition_id)].set_index("condition_id")
        assert len(rates)==129 and rates.index.is_unique and rates.angular_rate_available.all()
        assert rates.RHS_available.all() and (rates.R>=.1).all()
        selected[f"omega_t{t}"]=selected.condition_id.map(rates.native_RHS_angular_rate_deg_per_native_time)
        selected[f"R_t{t}"]=selected.condition_id.map(rates.R)
    selected.to_csv(OUT/"supporting/data/all129_native_rate_anchors.csv",index=False)
    draw(selected,10,"cue_on_rate_distance_size")
    draw(selected,0,"onset_rate_distance_size")
    stamp=datetime.now(timezone.utc).isoformat()
    text=f"""# Distance, size and global cue-on angular rate

Generated {stamp}. [Time10 PNG](cue_on_rate_distance_size.png) · [PDF](supporting/vectors/cue_on_rate_distance_size.pdf) · [SVG](supporting/vectors/cue_on_rate_distance_size.svg) · [onset context](onset_rate_distance_size.png) · [all129 exact case values and source archives](supporting/data/all129_native_rate_anchors.csv).

These are all129 prescribed placements in the narrow-range645-condition experiment, using eachT1000branch. Time10 and onset were already exported native-RHS anchors; this renderer performs no integration or additional RHS evaluation. Every displayed rate has retained full-state/snapshot support and R≥0.1 at that instant. Missing or weak values would not be assigned a phase, but none occurs at these two displayed times. Later weak intervals remain in complete trajectory sources.

The main time10 map is a systematic observation at a declared native clock, not a claim that it summarizes the whole trajectory. The onset map is context for finite-width input projection and initial moment effects. Native-rate histories from twelve common fully external placements and every saved-interval secant across all645branches remain linked in the full supplement. Native clocks are uncalibrated across models. Holding initial physical arc fixed changes target angle and basal architecture, so the size comparison does not isolate a causal geometric or incumbent-remodeling speed.

Circles, diamonds and triangles refer to complete initial supports inside, straddling or outside the10%-contrast region; centre membership and50% sensitivity are separately exported. These are spatial classes, not causal pathway identities. Lines connect tested points as visual guides without a fitted distance law. All source archive identities are retained in the table. All native kinetics, cue coefficients/supports, solver and qualified basals remain unchanged.
"""
    (OUT/"README.md").write_text(text)
    result=dict(generated_at=stamp,passed=True,conditions=129,observed_times=[0,10],usable_rates=258,new_integrations=0,new_RHS_evaluations=0,
                source_sha256={str(p):sha(p),str(q):sha(q)},renderer_sha256=sha(__file__),
                artifact_sha256={str(x.relative_to(OUT)):sha(x) for x in OUT.rglob("*") if x.is_file() and x.name!="verification.json"})
    (OUT/"supporting/metadata/verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(OUT/"cue_on_rate_distance_size.png")


if __name__=="__main__":main()
