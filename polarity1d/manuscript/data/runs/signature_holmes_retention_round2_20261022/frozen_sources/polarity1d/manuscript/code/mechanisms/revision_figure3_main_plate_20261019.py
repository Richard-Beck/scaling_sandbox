"""Final page composition from the three accepted current-experiment views.

Reads compact plotted-data CSVs only. No native integrations or new readouts.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[key]="1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
TOPIC=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
OUT=TOPIC/"main_plate"
NAMES=dict(wave_pinning="Wave-pinning",goryachev="Goryachev",holmes_model3="Holmes",spring="Spring")
COLORS=["#c98a20","#3c8894","#564d86"]
MARKERS=dict(inside="o",straddling="D",outside="^")


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(rel,sources):
    p=TOPIC/rel;sources[str(p)]=sha(p);return pd.read_csv(p,float_precision="round_trip")


def main():
    for leaf in ("data","metadata","vectors"):(OUT/"supporting"/leaf).mkdir(parents=True,exist_ok=True)
    sources={}
    cells=read("vector_construction/supporting/data/native_cell_contributions.csv",sources)
    sectors=read("vector_construction/supporting/data/exact_sector_sums.csv",sources)
    states=read("vector_construction/supporting/data/two_state_measurements.csv",sources)
    anchors=read("rate_map/supporting/data/all129_native_rate_anchors.csv",sources)
    fields=read("local_size_contrast/supporting/data/native_cue_on_fields.csv",sources)
    paths=read("local_size_contrast/supporting/data/complex_P_and_native_rate.csv",sources)
    cases=read("local_size_contrast/supporting/data/selected_conditions.csv",sources)
    assert len(states)==2 and len(anchors)==129 and set(fields.condition_id)=={39,149}
    plt.rcParams.update({"font.size":8,"axes.labelsize":8,"axes.titlesize":9,"xtick.labelsize":7.5,"ytick.labelsize":7.5,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(7.2,9.6));manifest={"A":{},"B":{},"C":{}}
    fig.text(.02,.986,"A  Local activity adds to a global vector",weight="bold",fontsize=10)
    norm=Normalize(0,float(cells.active.max()));cmap=plt.get_cmap("viridis");scale=.83
    for j,label in enumerate(["initial","cue_withdrawal"]):
        ax=fig.add_axes((.025+.27*j,.785,.25,.185));d=cells[cells.state==label].sort_values("cell");s=sectors[sectors.state==label].sort_values("sector");m=states[states.state==label].iloc[0]
        alpha=np.deg2rad(d.normal_angle_deg.to_numpy());n=len(d);half=np.pi/n
        start=np.c_[np.cos(alpha-half),np.sin(alpha-half)];end=np.c_[np.cos(alpha+half),np.sin(alpha+half)]
        ax.add_collection(LineCollection(np.stack([start,end],axis=1),array=d.active.to_numpy(),norm=norm,cmap=cmap,linewidth=7.8,capstyle="butt",zorder=3))
        for r in s.to_dict("records"):
            value=complex(r["normalized_contribution_real"],r["normalized_contribution_imag"]);tail=1.12*np.exp(2j*np.pi*r["sector"]/12);v=scale*value;length=abs(v)
            ax.arrow(tail.real,tail.imag,v.real,v.imag,width=min(.011,.18*length),head_width=min(.065,.55*length),head_length=min(.085,.45*length),length_includes_head=True,color="#327aa5",zorder=4)
        P=complex(m.P_real,m.P_imag);v=scale*P
        ax.arrow(0,0,v.real,v.imag,width=.021,head_width=.10,head_length=.14,length_includes_head=True,color="#b93e50",zorder=6)
        ax.scatter([0],[0],s=4,color="#b93e50",zorder=7)
        target=np.linspace(np.deg2rad(150),np.deg2rad(186),80);ax.plot(1.09*np.cos(target),1.09*np.sin(target),color="#cf852b",lw=1.6)
        ax.set_xlim(-1.52,1.52);ax.set_ylim(-1.52,1.52);ax.set_aspect("equal");ax.axis("off")
        centre=.15+.27*j
        fig.text(centre,.778,"Initial" if j==0 else "Cue withdrawal",ha="center",fontsize=8.5,weight="bold")
        fig.text(centre,.760,rf"$R={m.R:.3f}$",ha="center",fontsize=8.5,color="#b93e50")
        fig.text(centre,.743,rf"$M_0/M_0(0)={m.M0_ratio:.3f}$",ha="center",fontsize=8,color="#2d7864")
        manifest["A"][label]=dict(native_cells=n,P=[m.P_real,m.P_imag],R=m.R,M0_ratio=m.M0_ratio,sector_vector_sum=[float(s.normalized_contribution_real.sum()),float(s.normalized_contribution_imag.sum())],arrow_scale=scale,equal_xy=True)
    fig.text(.60,.966,"Holmes · C=30\nd₀=14 · T=300",fontsize=8.5,linespacing=1.5)
    fig.text(.60,.905,r"$M_0=C\langle u\rangle$",fontsize=10)
    fig.text(.60,.871,r"$M_1=C\langle u e^{2\pi iq}\rangle$",fontsize=10)
    fig.text(.60,.837,r"$P=M_1/M_1(0),\quad R=|P|$",fontsize=9.5)
    fig.text(.60,.810,"Blue: sector contributions\nRed: resultant P",fontsize=8,linespacing=1.4)
    cax=fig.add_axes((.63,.775,.28,.009));cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal",ticks=[0,1,2]);cb.ax.tick_params(labelsize=7,length=2,pad=1);cb.set_label("Native active concentration",fontsize=7.5,labelpad=2)
    fig.text(.02,.710,"B  Cue-on angular rate at t=10",fontsize=10,weight="bold")
    # One shared colour and footprint key for all four distance–size axes.
    fig.text(.07,.691,"C / C₀:",fontsize=8)
    for k,(color,label) in enumerate(zip(COLORS,["30 / 2.5","35 / 2.75","40 / 3"])):
        x=.22+.24*k;fig.add_artist(Line2D([x,x+.035],[.694,.694],transform=fig.transFigure,color=color,lw=1.4));fig.text(x+.045,.691,label,fontsize=8)
    for k,(label,marker) in enumerate(MARKERS.items()):
        x=.22+.24*k;fig.add_artist(Line2D([x+.015],[.670],transform=fig.transFigure,marker=marker,color="#52646f",lw=0,markersize=4));fig.text(x+.045,.667,label,fontsize=8)
    fig.text(.07,.667,"Full support:",fontsize=8)
    axesB=[]
    for i,(model,name) in enumerate(NAMES.items()):
        row,col=divmod(i,2);ax=fig.add_axes((.09+.48*col,.48-.195*row,.37,.145));axesB.append(ax)
        part=anchors[anchors.model==model];coordinates=[]
        for color,size in zip(COLORS,sorted(part.size_parameter.unique())):
            d=part[part.size_parameter==size].sort_values("physical_distance_at_onset");x=d.physical_distance_at_onset.to_numpy();y=d.omega_t10.to_numpy();ax.plot(x,y,color=color,lw=.9)
            for regime,marker in MARKERS.items():
                keep=d.footprint_relative_to_core_10==regime;ax.scatter(x[keep],y[keep],marker=marker,color=color,s=10,edgecolors="white",linewidths=.20,zorder=3)
            coordinates.extend(d[["condition_id","physical_distance_at_onset","omega_t10","size_parameter","footprint_relative_to_core_10"]].to_dict("records"))
        ax.set_title(name,loc="left",fontsize=9,weight="bold",pad=3);ax.axhline(0,color="#bcc5ca",lw=.45)
        ax.spines[["top","right"]].set_visible(False);ax.tick_params(length=2,pad=2)
        if model=="spring":ax.set_xticks([0,.5,1,1.25],["0","0.5","1","1.25"])
        else:ax.set_xticks([0,5,10,14])
        if row==1:ax.set_xlabel("Initial physical centre arc d₀",fontsize=8,labelpad=2)
        if col==0:ax.set_ylabel("ω (° / native time)",fontsize=8,labelpad=2)
        manifest["B"][model]=coordinates
    fig.text(.02,.234,"C  Wave-pinning: the same external physical cue",fontsize=10,weight="bold")
    normC=Normalize(0,float(fields.native_active.max()));axesC=[]
    for j,cid in enumerate([39,149]):
        ax=fig.add_axes((.09+.28*j,.059,.23,.145));axesC.append(ax)
        d=fields[fields.condition_id==cid];m=d.pivot(index="cue_time",columns="initial_physical_arc",values="native_active").sort_index();r=cases[cases.condition_id==cid].iloc[0]
        image=ax.pcolormesh(m.columns.to_numpy(),m.index.to_numpy(),m.to_numpy(),shading="nearest",cmap="viridis",norm=normC,rasterized=True)
        ax.set_yscale("symlog",linthresh=.1);ax.set_ylim(0,1000);ax.set_yticks([0,.1,1,10,100,1000],["0",".1","1","10","100","1000"]);ax.minorticks_off();ax.set_xlim(-20,20);ax.set_xticks([-20,0,20])
        if r.size_parameter<40:ax.axvspan(-20,-15,color="#e7eaec",zorder=0);ax.axvspan(15,20,color="#e7eaec",zorder=0)
        for edge in [r.initial_edge_left_10,r.initial_edge_right_10]:ax.axvline(edge,color="#f1b76f",ls="--",lw=.5)
        ax.axvline(11,color="white",ls=":",lw=.6);ax.plot([9.5,12.5],[0,0],color="#cf852b",lw=2.5,clip_on=False)
        ax.set_title(f"C={r.size_parameter:g}",color=COLORS[0 if j==0 else 2],loc="left",fontsize=8.5,pad=3)
        ax.set_xlabel("Initial physical arc",fontsize=7.8,labelpad=2);ax.spines[["top","right"]].set_visible(False);ax.tick_params(length=2,pad=2,labelsize=7)
        if j==0:ax.set_ylabel("Cue-on time",fontsize=8,labelpad=2)
        manifest["C"][str(cid)]=dict(shape=list(m.shape),time_min=float(m.index.min()),time_max=float(m.index.max()),physical_axis_limits=[-20,20],actual_C=r.circumference,field_samples=int(m.size),native_field_array_sha256=hashlib.sha256(np.ascontiguousarray(m.to_numpy()).tobytes()).hexdigest())
    cax=fig.add_axes((.614,.065,.008,.135));cb=fig.colorbar(image,cax=cax,ticks=[0,.5,1]);cb.ax.tick_params(labelsize=7,length=2,pad=1);cb.set_label("u",fontsize=8,labelpad=2)
    ax=fig.add_axes((.718,.044,.25,.171));axesC.append(ax)
    for cid,color in [(39,COLORS[0]),(149,COLORS[2])]:
        d=paths[paths.condition_id==cid].sort_values("frame");r=cases[cases.condition_id==cid].iloc[0]
        ax.plot(d.P_real,d.P_imag,color=color,lw=1.1);ax.scatter(d.P_real.iloc[-1],d.P_imag.iloc[-1],s=12,color=color,marker="s",zorder=4)
        theta=np.deg2rad(r.target_angle_in_P_coordinates_deg);ax.plot([0,np.cos(theta)],[0,np.sin(theta)],color=color,lw=.5,ls=":")
        for wanted in [100,300]:
            point=d.loc[(d.cue_time-wanted).abs().idxmin()];ax.scatter(point.P_real,point.P_imag,s=5,color=color);ax.annotate(f"{point.cue_time:.0f}",(point.P_real,point.P_imag),xytext=(2,3 if cid==39 else -7),textcoords="offset points",color=color,fontsize=6.8)
        manifest["C"][str(cid)]["complex_P_frames"]=len(d)
        manifest["C"][str(cid)]["complex_P_array_sha256"]=hashlib.sha256(np.ascontiguousarray(d[["P_real","P_imag"]].to_numpy()).tobytes()).hexdigest()
    angles=np.linspace(0,2*np.pi,300);ax.plot(np.cos(angles),np.sin(angles),ls="--",lw=.5,color="#b9c3c9");ax.axhline(0,color="#c8d1d6",lw=.45);ax.axvline(0,color="#c8d1d6",lw=.45);ax.scatter([1],[0],s=8,color="#213c49",zorder=4)
    ax.set_xlim(-1.2,1.2);ax.set_ylim(-1.2,1.2);ax.set_aspect("equal");ax.set_xticks([-1,0,1]);ax.set_yticks([-1,0,1]);ax.set_xlabel("Re P",fontsize=8,labelpad=2);ax.set_ylabel("Im P",fontsize=8,labelpad=1);ax.set_title("P path",fontsize=8.5,loc="left",pad=3);ax.spines[["top","right"]].set_visible(False);ax.tick_params(length=2,pad=2,labelsize=7)
    fig.canvas.draw()
    for axis in [fig.axes[0],fig.axes[1],ax]:
        p0=axis.transData.transform([0,0]);sx=np.linalg.norm(axis.transData.transform([1,0])-p0);sy=np.linalg.norm(axis.transData.transform([0,1])-p0);assert abs(sx-sy)<1e-10
    for ext in ("png","pdf","svg"):
        path=OUT/"figure3_main.png" if ext=="png" else OUT/"supporting/vectors"/("figure3_main."+ext)
        fig.savefig(path,dpi=300,facecolor="white")
    plt.close(fig)
    pd.DataFrame(anchors).to_csv(OUT/"supporting/data/displayed_rate_points.csv",index=False)
    stamp=datetime.now(timezone.utc).isoformat()
    report=dict(generated_at=stamp,passed=True,figure_inches=[7.2,9.6],purpose="Final page composition from three existing accepted data views; no new science.",new_integrations=0,new_RHS_evaluations=0,
                renderer_sha256=sha(__file__),source_sha256=sources,artist_mapping=manifest,equal_xy_axes=3,nominal_font_points=dict(axis=8,tick=7.5,panel_label=10,time_annotation=6.8),
                caption_word_count=len((OUT/"CAPTION.md").read_text().split()) if (OUT/"CAPTION.md").exists() else None,
                artifact_sha256={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob("*") if p.is_file() and p.name!="verification.json"})
    (OUT/"supporting/metadata/verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(OUT/"figure3_main.png")


if __name__=="__main__":main()
