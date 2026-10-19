"""Display the qualified narrow-range basals and every frozen cue placement."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,os
for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[k]="1"
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/"data/runs/revision_arc_narrow_grid_20261019"
OUT=ROOT/"figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment"
MODELS=[("wave_pinning","Wave-pinning"),("goryachev","Goryachev"),("holmes_model3","Holmes"),("spring","Spring")]
COLORS=dict(inside="#287ba0",straddling="#c58a25",outside="#a84b5c")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    plan=json.loads((RUN/"manifest.json").read_text());assert sha(RUN/"manifest.json")==json.loads((RUN/"freeze.json").read_text())["manifest_sha256"]
    plt.rcParams.update({"font.size":9.5,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig,axes=plt.subplots(4,3,figsize=(14.5,18),sharey=True);sources={str(RUN/"manifest.json"):sha(RUN/"manifest.json")}
    for ri,(model,name) in enumerate(MODELS):
        basals=sorted([b for b in plan["basals"] if b["model"]==model],key=lambda b:b["size_parameter"])
        half=max(b["basal_circumference"] for b in basals)/2
        for ci,b in enumerate(basals):
            cases=[c for c in plan["cases"] if c["model"]==model and c["size_parameter"]==b["size_parameter"] and c["cue_duration"]==1000.]
            cases.sort(key=lambda c:c["physical_distance_at_onset"]);c=cases[0];C=b["basal_circumference"];ax=axes[ri,ci]
            path=RUN/b["input"];assert sha(path)==b["input_sha256"];sources[str(path)]=sha(path)
            with np.load(path) as z:state=z["state"].copy()
            with np.load(RUN/c["cue_input"]) as z:q=z["material_positions"].copy()
            u=state[:b["cells"]] if model=="spring" else state[0]
            x=((q-b["preparation_axis_fraction"]+.5)%1-.5)*C;order=np.argsort(x)
            for lo,hi in ((-half,-C/2),(C/2,half)):
                if hi>lo+1e-12:ax.axvspan(lo,hi,color="#f1f3f4",zorder=0)
            ax.axvspan(c["initial_edge_left_10"],c["initial_edge_right_10"],color="#f3f5f6",zorder=0)
            ax.axvspan(c["initial_edge_left_50"],c["initial_edge_right_50"],color="#e4edf3",zorder=0)
            ax.plot(x[order],((u-u.min())/np.ptp(u))[order],color="#172f3c",lw=1.6)
            for suffix,style,color in (("50","--","#6a7c88"),("10",":","#8998a2")):
                for side in ("left","right"):ax.axvline(c[f"initial_edge_{side}_{suffix}"],ls=style,color=color,lw=.8)
            for j,case in enumerate(cases):
                d=case["physical_distance_at_onset"];w=case["full_physical_support"];y=-.15-.12*j;color=COLORS[case["footprint_relative_to_core_10"]]
                for shift in (-C,0,C):
                    lo=max(-C/2,d-w/2+shift);hi=min(C/2,d+w/2+shift)
                    if hi>lo:ax.plot([lo,hi],[y,y],color=color,lw=4.2,solid_capstyle="butt")
                ax.plot(d,y,'o',color=color,ms=3.7,mec="white",mew=.45)
                ax.text(-half+.04*half,y,f'd={d:.4g}',color=color,fontsize=8,ha="left",va="center",bbox=dict(fc="white",ec="none",alpha=.88,pad=.3))
            size=f'C₀={b["size_parameter"]:g}; Cᵦ={C:.5f}' if model=="spring" else f'C={C:g}'
            ax.set_title(f'{name} · {size}\n10% edge = {c["initial_edge_right_10"]:.5g}',loc="left",fontsize=10.8)
            ax.set_xlim(-half,half);ax.set_ylim(-1.49,1.1);ax.set_yticks([0,.5,1]);ax.spines[["top","right"]].set_visible(False);ax.tick_params(labelsize=8.5)
            if ci==0:ax.set_ylabel("Activity above minimum / basal contrast\nFull cue supports below",fontsize=9.5)
            if ri==3:ax.set_xlabel("Signed initial physical arc from prepared axis",fontsize=9.3)
    fig.subplots_adjust(left=.075,right=.985,top=.865,bottom=.13,wspace=.20,hspace=.34)
    fig.text(.075,.97,"Qualified narrow-range design: internal, edge and outside cues",fontsize=19,weight="bold")
    fig.text(.075,.944,"All 129 frozen placements on the actual 12 basals; 645 cue-duration conditions. No parameters or basals are retuned.",fontsize=11)
    fig.text(.075,.923,"Bars: full physical support. Dots: cue centres. Colours classify the whole footprint against the initial 10%-contrast region.",fontsize=10)
    for j,(name,color) in enumerate(COLORS.items()):fig.text(.075+.22*j,.076,name.capitalize()+" support",color=color,fontsize=11,weight="bold")
    fig.text(.075,.054,"Black profile: [u − min(u)] / [max(u) − min(u)]. Dashed edges: 50%; dotted edges: 10%. These are operational activity thresholds.",fontsize=9.7)
    fig.text(.075,.034,"Common physical x scales within model rows; grey outer areas lie beyond smaller rings. Supports wrap periodically near the antipode.",fontsize=9.7)
    fig.text(.075,.014,"Spring distances are initial physical arcs of stationary material centres; subsequent physical arcs vary with C(t). No site ancestry or cell-motion claim.",fontsize=9.7)
    vectors=OUT/"supporting/vectors";vectors.mkdir(parents=True,exist_ok=True)
    outputs=[]
    for ext in ("png","pdf","svg"):
        path=OUT/"narrow_initial_design.png" if ext=="png" else vectors/("narrow_initial_design."+ext)
        fig.savefig(path,dpi=180,facecolor="white");outputs.append(path)
    plt.close(fig)
    source=OUT/"supporting/metadata/narrow_initial_design.json"
    source.write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),native_integration="none",qualified_basals=12,placements=129,conditions=645,
        source_sha256=sources,renderer_sha256=sha(__file__),artifact_sha256={str(p):sha(p) for p in outputs}),indent=2)+"\n")
    print(OUT/"narrow_initial_design.png")


if __name__=="__main__":main()
