"""Retained-data Figure 2 components; no integrations or new main-figure variant.

The common active-incumbent recovery assay and global polarity component are
owned separately. Legacy inactive-withdrawal coefficients are never imported.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import re
import sys

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.patches import Polygon, Patch, Rectangle
from matplotlib.lines import Line2D
import length_dev_output_layout as layout

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
TOPIC = ROOT / "figures/length_mechanisms_20261008/dev/expanded_ring_components_20261017"
RUN = ROOT / "data/runs/chemical_ring_functional_v2_20261017"
PAIR_CSV = ROOT / "figures/length_mechanisms_20261008/dev/chemical_ring_functional_20261017/supporting/data/ring_functional_ring_linear_mirror_comparison.csv"
MODELS = ("wave_pinning", "goryachev", "holmes_model3")
NAMES = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "holmes_model3":"Holmes 3"}
LABELS = ("Single peak near old; identity unresolved", "Cue-associated peak", "Distinct old + cue peaks",
          "One region spans old + cue", "Shifted / merged; unresolved", "No resolved peak", "Ambiguous multiple peaks")
COLORS = ("#4678a3", "#d77c35", "#8b65a9", "#69a7a0", "#bcb069", "#ececec", "#6e777c")
EXAMPLES = (("1",30,12,"Cue-associated single peak"), ("2",40,12,"Merged / bridged region"),
            ("3",40,16,"Two distinct associated peaks"))
ARTIFACTS = []
SOURCES = {}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b""):
            digest.update(chunk)
    return digest.hexdigest()


def source(path, expected=None):
    path = Path(path).resolve()
    if str(path) not in SOURCES:
        SOURCES[str(path)] = sha(path)
    if expected is not None:
        assert SOURCES[str(path)] == expected, path
    return path


def read_csv(path):
    with source(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(row, key, default=np.nan):
    try:return float(row[key]) if row[key] not in ("",None) else default
    except (KeyError,ValueError,TypeError):return default


def out(name):
    return layout.dev_output_path(TOPIC,name)


def csvwrite(name, rows):
    path = out(name)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    ARTIFACTS.append(path)


def export(fig, stem, caption):
    for extension in (".png",".pdf",".svg"):
        path=out(stem+extension)
        fig.savefig(path,dpi=200 if extension==".png" else None,bbox_inches="tight")
        ARTIFACTS.append(path)
    plt.close(fig)
    path=out(stem+"_caption.txt");path.write_text(caption+"\n");ARTIFACTS.append(path)


def example_row(matrix, circumference, distance):
    matches=[r for r in matrix if r["model"]=="holmes_model3" and number(r,"circumference")==circumference
             and np.isclose(number(r,"short_route"),distance)]
    assert len(matches)==1 and matches[0]["status"]=="complete"
    return matches[0]


def baseline(pairs):
    # Figure 1B now uses the native finite-volume cue-weighted derivative.
    # Compare its endpoint support3 with the reflected full-ring support6;
    # retain the source point and actual full-support3 ring readouts separately.
    for row in pairs:
        L,C=number(row,"linear_L"),number(row,"ring_C")
        with np.load(row["linear_checkpoint"]) as z:
            linear=z["state"].copy();lx=z["x"].copy()
        with np.load(row["ring_checkpoint"]) as z:
            ring=z["state"].copy();rx=z["x"].copy()
        divisor=7.5 if row["model"]=="holmes_model3" else 1.
        linear_distance=L-lx
        linear_cue=np.where(linear_distance<=3,.5*(1+np.cos(np.pi*np.minimum(linear_distance/3,1))),0.)
        linear_weights=np.full(len(lx),L/len(lx))
        center=(int(np.argmax(ring[0]))+len(rx)//2)%len(rx)
        ring_distance=np.abs((rx-rx[center]+C/2)%C-C/2)
        ring_cue=np.where(ring_distance<=3,.5*(1+np.cos(np.pi*np.minimum(ring_distance/3,1))),0.)
        ring_weights=np.full(len(rx),C/len(rx))
        linear_gain=float((linear_cue*linear[1]/divisor)@linear_weights/(linear_cue@linear_weights))
        ring_gain=float((ring_cue*ring[1]/divisor)@ring_weights/(ring_cue@ring_weights))
        row.update(linear_end_cue_support=3,ring_matched_cue_full_support=6,
                   linear_matched_cue_weighted_gain=linear_gain,ring_matched_cue_weighted_gain=ring_gain,
                   matched_cue_weighted_gain_ring_over_linear=ring_gain/linear_gain,
                   actual_ring_challenge_full_support=3)
    metrics=(
        ("Full half-contrast width",2.,"ring_half_width_over_linear"),
        ("Full excess-active core amount",2.,"ring_half_excess_over_linear"),
        ("Absolute primary peak",1.,"peak_ring_over_linear"),
        ("Matched cue-weighted input gain",1.,"matched_cue_weighted_gain_ring_over_linear"),
        ("Opposing-point gain (support)",1.,"point_gain_ring_over_linear"),
    )
    table=[]; summaries=[]
    for label,expected,key in metrics:
        row=[label,f"{expected:g}"]
        for model in MODELS:
            values=np.array([number(r,key)*expected for r in pairs if r["model"]==model])
            row.append(f"{values.min():.6f}–{values.max():.6f}")
            summaries.append(dict(model=model,observable=label,expected_ring_over_line=expected,
                                  measured_ratio_min=float(values.min()),measured_ratio_max=float(values.max()),
                                  maximum_relative_deviation=float(np.max(abs(values/expected-1)))))
        table.append(row)
    fig,ax=plt.subplots(figsize=(13.5,4.5));ax.axis("off")
    headers=["Measured ring / interval ratio","Mirror expectation","Wave-pinning (4 pairs)","Goryachev (4 pairs)","Holmes 3 (5 pairs)"]
    tab=ax.table(cellText=table,colLabels=headers,cellLoc="center",colWidths=(.27,.14,.197,.197,.197),bbox=(0,.38,1,.51))
    tab.auto_set_font_size(False);tab.set_fontsize(10)
    for (row,col),cell in tab.get_celld().items():
        cell.set_edgecolor("#d4d9dd");cell.set_linewidth(.6)
        if row==0:cell.set_facecolor("#e8eef2");cell.set_text_props(weight="bold")
        elif col==0:cell.set_text_props(ha="left")
    ax.set_title(r"Native baseline correspondence: ring circumference $C=2L$",fontsize=16,pad=18)
    ax.text(0,.27,"WP/G: L = 5, 15, 25, 35  ↔  C = 10, 30, 50, 70\nHolmes: L = 10, 15, 20, 25, 30  ↔  C = 20, 30, 40, 50, 60",fontsize=10,va="top",transform=ax.transAxes)
    ax.text(.57,.27,"Matched cue: endpoint support 3 ↔ ring full support 6.\nPoint gain remains a distinct supporting readout.",fontsize=10,va="top",transform=ax.transAxes)
    ax.text(0,.06,"Independently prepared native basals; all 26 checkpoint hashes verified. Ring h = 0.05; interval uses 201 cells.\nA support-3 endpoint cue reflects to full ring support 6; the retained ring challenges use support 3. Finite cues are not mirror-equivalent.",fontsize=10,va="top",transform=ax.transAxes)
    fig.tight_layout()
    export(fig,"expanded_ring_baseline_correspondence","Compact measured C=2L correspondence for all 13 retained native ring/interval pairs, without interpolation or substitution of checkpoints. Full ring width and excess-active core amount are approximately twice their interval counterparts; absolute peak and native cue susceptibility are approximately equal. The matched cue-weighted gain applies the same Figure1B finite-volume average of the native cue derivative to endpoint support3 and its reflected full-ring support6. The opposing-point gain is retained as a distinct supporting readout. Every entry reports the minimum–maximum ratio across all listed matches. These are independently prepared states on corresponding basal branches, not byte-identical reflected full states and not a continuum-convergence result. The ring has mesh spacing0.05; the 201-cell interval mesh spacing varies withL. The actual map/outcome ring challenges retain full support3 and are not finite-cue-equivalent to the interval support3 protocol. No integration is needed to evaluate these basal cue derivatives. The new common active-incumbent recovery baseline is being measured separately; legacy inactive-withdrawal coefficients are not displayed or substituted.")
    csvwrite("expanded_ring_baseline_ratio_summary.csv",summaries)
    # Deliberately exclude legacy inactive-withdrawal coefficients from this
    # new component's data table; original source CSV remains immutable.
    legacy={"linear_conditional_repair","ring_basal_conditional_repair","repair_ring_over_linear"}
    csvwrite("expanded_ring_baseline_pair_measurements.csv",[{k:v for k,v in r.items() if k not in legacy} for r in pairs])
    return summaries


def category(text):
    if text.startswith("one old-core-associated") or text=="one peak associated with old core":return 0
    if text=="one peak associated with cue patch":return 1
    if text=="distinct peaks associated with old core and cue patch":return 2
    if text.startswith("one connected region spanning"):return 3
    if text.startswith("one shifted"):return 4
    if text.startswith("no peak"):return 5
    return 6


def map_cells(ax, data, model, sizes, key, norm, cmap, categorical=False, highlights=False):
    for iy,circumference in enumerate(sizes):
        rows=sorted([r for r in data if r["model"]==model and number(r,"circumference")==circumference],key=lambda r:number(r,"fraction"))
        centers=np.array([number(r,"realized_fraction",number(r,"fraction")) for r in rows])
        edges=np.r_[centers[0]-.0125,(centers[:-1]+centers[1:])/2,centers[-1]+.0125]
        for j,row in enumerate(rows):
            complete=row["status"]=="complete"
            preparation=row["status"]=="unqualified_preparation"
            value=category(row[key]) if complete and categorical else number(row,key)
            face=cmap(norm(value)) if complete else "#777777" if preparation else "#dedede"
            ax.add_patch(Rectangle((edges[j],iy-.5),edges[j+1]-edges[j],1,facecolor=face,edgecolor="white",lw=.3))
            if not complete:
                ax.plot(centers[j],iy,"_" if preparation else "x",color="white" if preparation else "#888888",ms=3,mew=.65)
            if highlights and model=="holmes_model3" and complete:
                for label,C,d,_ in EXAMPLES:
                    if C==circumference and np.isclose(number(row,"short_route"),d):
                        ax.text(centers[j],iy,label,ha="center",va="center",fontsize=8,fontweight="bold",
                                bbox=dict(boxstyle="circle,pad=.18",facecolor="white",edgecolor="black",lw=1.2))
    ax.set(xlim=(.0875,.5125),ylim=(-.5,len(sizes)-.5),xlabel="Realized shorter-route fraction",ylabel="Circumference (native units)")
    ax.set_xticks((.1,.2,.3,.4,.5));ax.set_yticks(range(len(sizes)),[f"{s:g}" for s in sizes]);ax.tick_params(axis="y",labelsize=8)
    ax.set_title(NAMES[model],fontsize=12,pad=9)


def gain_maps(matrix,sizes):
    references={model:number(next(r for r in matrix if r["model"]==model and number(r,"circumference")==30 and number(r,"fraction")==.5),"initial_cue_weighted_gain") for model in MODELS}
    data=[]
    for row in matrix:
        rec={key:row.get(key,"") for key in ("model","circumference","fraction","realized_fraction","short_route","long_route","status","initial_cue_weighted_gain")}
        rec["own_model_C30_antipodal_reference"]=references[row["model"]]
        rec["gain_over_own_C30_antipodal"]=number(row,"initial_cue_weighted_gain")/references[row["model"]] if row["status"]=="complete" else ""
        data.append(rec)
    values=[number(r,"gain_over_own_C30_antipodal") for r in data if r["status"]=="complete"]
    assert min(values)>0 and max(values)<=3
    norm=Normalize(0,3);cmap=plt.get_cmap("viridis")
    fig,axes=plt.subplots(1,3,figsize=(13.5,6.7))
    for ax,model in zip(axes,MODELS):map_cells(ax,data,model,sizes,"gain_over_own_C30_antipodal",norm,cmap)
    fig.suptitle("Initial challenged-site gain / each model's own C30 antipodal gain",fontsize=15,y=.99)
    fig.subplots_adjust(left=.055,right=.91,bottom=.27,top=.88,wspace=.25)
    cbax=fig.add_axes((.93,.27,.016,.53));cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbax);cb.set_ticks((0,1,2,3));cb.set_label("Dimensionless onset gain")
    fig.text(.055,.095,"One shared 0–3 scale retains nearly flat location responses: largest within-size range 0.0392% (WP), 0.0891% (G), 0.00208% (Holmes).\nShorter route = fraction × C; complementary route = C − shorter route. Every placement at a given model/C uses its exact same complete basal state.",fontsize=9.5)
    fig.legend(handles=[Line2D([],[],marker="x",ls="",color="#888888",label="198 cue/core-overlap exclusions"),
                        Patch(facecolor="#777777",label="17 unavailable Holmes C80 preparations")],loc="lower center",ncol=2,frameon=False,fontsize=9,bbox_to_anchor=(.5,.005))
    export(fig,"expanded_ring_initial_gain_maps","All 867 prospective size/location cells. The native exact cue-RHS derivative is cue-weighted inactive substrate in Wave-pinning/Goryachev and cue-weighted inactive Rac divided by 7.5 in Holmes. Values are divided by each model's own C30 antipodal (d15) value and displayed on the same fixed 0–3 color scale, without clipping; all 652 completed gains lie inside this range. This within-model normalization does not cross-calibrate native time or cue units. The onset derivative is not the total basal activation propensity. Realized shorter-route fractions retain center rounding; both physical routes are exported. Gray crosses preserve 198 cue/core-overlap exclusions; dark gray preserves the 17 unavailable Holmes C80 preparations. Location dependence is nearly absent and is not magnified by separate narrow scales. The new common recovery map is not part of this retained-data component.")
    csvwrite("expanded_ring_initial_gain_map.csv",data)
    return references,float(min(values)),float(max(values))


def outcome_maps(matrix,sizes):
    fig,axes=plt.subplots(1,3,figsize=(13.5,7.2))
    for ax,model in zip(axes,MODELS):map_cells(ax,matrix,model,sizes,"final_location_association",Normalize(-.5,6.5),ListedColormap(COLORS),categorical=True,highlights=True)
    fig.suptitle("Complete endpoint context at t = 14000, with the three Holmes histories highlighted",fontsize=14.5,y=.99)
    fig.subplots_adjust(left=.055,right=.99,bottom=.34,top=.88,wspace=.25)
    handles=[Patch(facecolor=color,label=label) for color,label in zip(COLORS,LABELS)]
    handles.extend([Line2D([],[],marker="x",ls="",color="#888888",label="Cue/core overlap excluded"),Patch(facecolor="#777777",label="Holmes C80 unavailable")])
    fig.legend(handles=handles,loc="lower center",ncol=3,frameon=False,fontsize=8.5,bbox_to_anchor=(.5,.035))
    fig.text(.055,.205,"① Holmes C30 d12: cue-associated single peak     ② Holmes C40 d12: merged / bridged region     ③ Holmes C40 d16: two distinct sites\nCircled cells use exact native h = 0.05 histories. Endpoint labels are descriptive geographic associations, not stable-branch or tracked-fate classifications.",fontsize=9.5)
    export(fig,"expanded_ring_outcome_context","The complete retained 867-cell endpoint context for Wave-pinning, Goryachev and Holmes at 14000 native units, with all preparation/geometry exclusions retained. The three highlighted exact Holmes cells are 1: C30 d12 (fraction 0.4), 2: C40 d12 (fraction 0.3), and 3: C40 d16 (fraction 0.4), corresponding to the linked native kymograph component. These illustrative conditions were selected using already observed activity patterns, not an outcome-blind validation set. The C30-to-C40 comparison changes circumference, material, architecture and complementary route; the two C40 conditions share identical complete basal Rac/Rho fields and change cue location. Geographic endpoint associations alone do not establish site identity or stability; the full native histories provide the selected examples' temporal context.")
    csvwrite("expanded_ring_outcome_map.csv",[{key:row.get(key,"") for key in ("model","circumference","fraction","realized_fraction","short_route","long_route","status","final_peak_count","final_location_association","final_single_region_spans_old_and_cue","checkpoint_sha256")} for row in matrix])


def native_histories(matrix, audited_sources):
    records=[]; data={}; histories=[]
    for label,C,d,role in EXAMPLES:
        row=example_row(matrix,C,d);assert number(row,"dx")==.05 and number(row,"cells")==C/.05
        phases=[]
        for phase,offset in (("cue",0),("recovery",1000)):
            path=source(row[phase+"_path"],audited_sources[row[phase+"_path"]])
            with np.load(path) as z:
                assert "states" in z.files and z["states"].shape[1]==4
                state=z["states"].copy();time=z["time"].copy()+offset;x=z["x"].copy();stimulus=z["stimulus"].copy()
            phases.append((time,x,state,stimulus))
            records.append(dict(example=label,role=role,circumference=C,short_route=d,long_route=C-d,
                                fraction=number(row,"fraction"),dx=.05,phase=phase,phase_offset=offset,
                                source_path=str(path),source_sha256=sha(path),checkpoint_sha256=row["checkpoint_sha256"],
                                native_fields=4,cells=len(x),source_samples=len(time),display_start=0,display_end=4000,
                                cue_full_support=3,cue_peak=1,cue_duration=1000,final_location_association=row["final_location_association"]))
        np.testing.assert_array_equal(phases[0][1],phases[1][1]);np.testing.assert_array_equal(phases[0][2][-1],phases[1][2][0])
        time=np.r_[phases[0][0],phases[1][0][1:]];primary=np.r_[phases[0][2][:,0],phases[1][2][1:,0]]
        assert np.all(np.diff(time)>0) and time[0]==0 and time[-1]==14000
        choose=time<=4000
        # Put the exact original peak coordinate at zero and cut the periodic
        # display at -10. Column permutation preserves every native sample.
        status=json.loads(source(row["status_path"],audited_sources[row["status_path"]]).read_text())
        origin=phases[0][1][status["core"]["peak_cell"]]
        display_x=(phases[0][1]-origin+10)%C-10
        order=np.argsort(display_x);display_x=display_x[order]
        display_primary=primary[choose][:,order]
        data[f"example{label}_time"]=time[choose];data[f"example{label}_x"]=display_x;data[f"example{label}_primary"]=display_primary
        histories.append((label,C,d,role,time[choose],display_x,display_primary,primary.max(),phases[0][2][0],phases[0][3]))
    # Location control really shares every native field, not just primary Rac.
    np.testing.assert_array_equal(histories[1][-2],histories[2][-2])
    np.testing.assert_array_equal(np.sort(histories[1][-1]),np.sort(histories[2][-1]))
    assert records[0]["checkpoint_sha256"]!=records[2]["checkpoint_sha256"]==records[4]["checkpoint_sha256"]
    path=out("expanded_ring_Holmes_native_display_samples.npz");np.savez_compressed(path,**data);ARTIFACTS.append(path)
    csvwrite("expanded_ring_Holmes_example_sources.csv",records)
    vmax=max(h[7] for h in histories);assert all(h[6].min()>=0 for h in histories)
    fig,axes=plt.subplots(1,3,figsize=(13.5,6.2),sharex=True,sharey=True)
    for ax,(label,C,d,role,time,x,primary,_,_,_) in zip(axes,histories):
        ax.set_facecolor("#bdbdbd")
        # A periodic duplicate closes the cell intersected by the plot cut;
        # the interior arrays are unchanged and no new field is interpolated.
        closed_x=np.r_[x,x[0]+C];closed_primary=np.c_[primary,primary[:,0]]
        mesh=ax.pcolormesh(closed_x,time,closed_primary,cmap="magma",vmin=0,vmax=vmax,shading="nearest",rasterized=True)
        ax.hlines(1000,-10,C-10,color="#83d6e2",ls="--",lw=.9)
        ax.axvline(d,color="#83d6e2",ls=":",lw=.7)
        if C==30:
            ax.axvline(20,color="black",lw=.8)
            ax.text(25,2100,"outside C30 ring",ha="center",va="center",rotation=90,color="#4c4c4c",fontsize=10)
        ax.set(xlim=(-10,30),ylim=(0,4000),xlabel="Position relative to original peak (native units)",
               xticks=(-10,0,10,20,30),yticks=(0,1000,2000,3000,4000))
        ax.set_title(f"{label}  Holmes C{C}, d{d}\n{role}",fontsize=11,pad=10)
    axes[0].set_ylabel("Time from cue onset (native units)")
    fig.subplots_adjust(left=.06,right=.925,bottom=.28,top=.80,wspace=.15)
    cbax=fig.add_axes((.942,.26,.014,.47));cb=fig.colorbar(mesh,cax=cbax);cb.set_label("Native active Rac (shared unclipped scale)")
    fig.suptitle("Three exact native Holmes histories: cue on 0–1000, then cue free",fontsize=15,y=.97)
    fig.text(.06,.08,"1 → 2 changes circumference: complementary route 18 → 28, with changed material and incumbent architecture.\n2 → 3 relocates the cue on the identical complete C40 Rac/Rho state. Full source histories extend to t = 14000; this display is linear t = 0–4000.\nPeriodic cut at −10 keeps the original site intact at 0; −10 and C−10 are identical endpoints. All three use h = 0.05, cue peak 1 and support 3.",fontsize=9.5)
    export(fig,"expanded_ring_Holmes_three_histories","Three exact native h=0.05 Holmes Rac histories selected from the highlighted complete source map: 1 C30 d12 has a cue-associated surviving single peak; 2 C40 d12 forms one merged/bridged region spanning original and challenged coordinates; 3 C40 d16 retains two distinct associated sites. The shared absolute active-Rac scale starts at zero and its maximum is the largest value in all three full 14000-unit histories, so nothing in the linear 0–4000 display is clipped. The shared physical coordinate is -10..30 relative to the original peak at0. Periodic column reindexing keeps that original site intact; -10 and C-10 are the same periodic endpoint. C30 occupies -10..20, with gray20..30 outside its domain, not zero activity. Cues remain at12/16. A repeated boundary column closes the periodic image without interpolation. Cyan horizontal lines mark cue withdrawal at1000. The C30/C40 d12 comparison changes size, conserved material, incumbent architecture and complementary route18/28 together. The two C40 histories share exactly identical complete basal Rac/Rho states and cue profile values while cue location differs. Different recovery sampling (50 versus5 native units in the later source histories) limits temporal precision; no event-time interpolation is used. The saved display arrays preserve exact native samples and fields, with no regridding. These selected histories supply temporal context beyond endpoint geographic flags; they do not establish universal stability or a pure distance effect.")
    return records,vmax


def handoff(references,gainrange,vmax):
    path=TOPIC/"COMPONENTS.md"
    text=f"""# Expanded ring Figure 2 components

Generated {datetime.now(timezone.utc).isoformat()}. Retained-data component package for one forthcoming updated checkpoint. No native integration, new scientific condition or main-figure variant was produced here. The [revision plan]({os.path.relpath(REPO/'current manuscript figs/figure_revision_plan_20261017.md',TOPIC).replace(' ','%20')}) and [ring evidence design]({os.path.relpath(REPO/'current manuscript figs/ring_figure2_evidence_design_20261017.md',TOPIC).replace(' ','%20')}) define the intended sequence.

1. [Baseline correspondence](expanded_ring_baseline_correspondence.png) ([PDF](supporting/vectors/expanded_ring_baseline_correspondence.pdf), [SVG](supporting/vectors/expanded_ring_baseline_correspondence.svg)) summarizes all 13 exact native size matches at C = 2L. Full width and excess-active core amount double; peak and native input gain agree closely. Its matched cue-weighted gain uses the same Figure1B finite-volume average and a geometrically reflected input: linear endpoint support3 versus ring full support6. Point gains remain distinct supporting readouts. The [13-pair table](supporting/data/expanded_ring_baseline_pair_measurements.csv) retains both definitions, raw native values, declared cue widths, mesh sizes and all26checkpoint hashes; the original source CSV is unchanged. These independently prepared basals are not byte-identical reflected states. The Figure1 L19/L39 highlights have no exact C38/C78 counterpart in the retained ring grid.
2. [Initial gain maps](expanded_ring_initial_gain_maps.png) ([PDF](supporting/vectors/expanded_ring_initial_gain_maps.pdf), [SVG](supporting/vectors/expanded_ring_initial_gain_maps.svg)) retain all 867 planned cells, including 198 overlap exclusions and 17 unavailable Holmes C80 preparations. Each model is divided by its own exact C30 antipodal gain; the shared fixed 0–3 scale contains the full range {gainrange[0]:.6g}–{gainrange[1]:.6g}. Tiny location effects remain visually tiny. [Map data](supporting/data/expanded_ring_initial_gain_map.csv) retain both realized physical routes and every status. Reference gains are WP {references['wave_pinning']:.9g}, Goryachev {references['goryachev']:.9g}, Holmes {references['holmes_model3']:.9g} in native units.
3. [Complete outcome context](expanded_ring_outcome_context.png) ([PDF](supporting/vectors/expanded_ring_outcome_context.pdf), [SVG](supporting/vectors/expanded_ring_outcome_context.svg)) retains all three models and marks the three selected Holmes source cells. Labels describe endpoint geography at14000; they do not establish tracked fate or stable branches. [Map data](supporting/data/expanded_ring_outcome_map.csv) retain all planned cells and flags.
4. [Three Holmes histories](expanded_ring_Holmes_three_histories.png) ([PDF](supporting/vectors/expanded_ring_Holmes_three_histories.pdf), [SVG](supporting/vectors/expanded_ring_Holmes_three_histories.svg)) use exact native h=0.05 histories: 1 C30 d12, cue-associated single peak; 2 C40 d12, merger/bridge; 3 C40 d16, two distinct sites. All use the same unclipped native active-Rac range 0–{vmax:.9g} and linear0–4000 time. Periodic column reindexing gives common x−10..30 with the original peak intact at0, cues12/16 and identical periodic endpoints−10/C−10; C30 occupies−10..20 with gray20..30 outside. Comparison1→2 changes size, material, incumbent architecture and complementary route; comparison2→3 relocates the same cue within an exactly identical complete C40 Rac/Rho starting state. [Source table](supporting/data/expanded_ring_Holmes_example_sources.csv) and [exact native display arrays](supporting/data/expanded_ring_Holmes_native_display_samples.npz) preserve the three cases without replacing or modifying raw sources.

The baseline now compares **the same cue-weighted susceptibility definition as Figure1B**, with geometrically matched endpoint support3 and reflected ring full support6. The original point gains remain linked in the table as distinct supporting readouts. The maps and finite ring histories retain their original full-support3 challenge protocol; their weighted gain is therefore not the baseline's full-support6 average, even though remote substrate is nearly uniform and the numerical differences are small. Finite-cue dose/coverage and response thresholds are not mirror-equivalent.

**The common active-incumbent recovery baseline and challenged-state recovery maps are pending separate native measurements.** This package deliberately imports no legacy inactive-withdrawal coefficients into new Figure2C. Those values cannot fill the new assay. Basal recovery has no future cue-location variable; any location map must use a declared challenged snapshot under the finalized assay. The global activity-polarity vector component is owned separately and is not duplicated here.

All native source files are referenced in place. The [component manifest](supporting/metadata/expanded_ring_components.manifest.json) records sources, artifact hashes, native continuity checks and display definitions. The previous checkpoint and all prior report pixels remain provenance. Reproduce with:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/expanded_ring_components_20261017.py
```

These components are ready for assembly into one updated checkpoint after the new common-assay results arrive. The combined figure still requires its fresh actual scientific and visual critique before the user's requested pause.
"""
    path.write_text(text);ARTIFACTS.append(path)


def main():
    source(__file__)
    pairs=read_csv(PAIR_CSV);matrix=read_csv(RUN/"results/matrix.csv")
    audit=json.loads(source(RUN/"results/audit.json").read_text());assert audit["passed"] and len(matrix)==867
    assert len(pairs)==13 and all(number(r,"ring_C")==2*number(r,"linear_L") for r in pairs)
    for row in pairs:
        source(row["linear_checkpoint"],row["linear_checkpoint_sha256"]);source(row["ring_checkpoint"],row["ring_checkpoint_sha256"])
    for name in ("figure_revision_plan_20261017.md","ring_figure2_evidence_design_20261017.md"):
        source(REPO/"current manuscript figs"/name)
    source(Path(__file__).with_name("all_five_linear_design_evidence_20261017.py"))
    audited_sources=json.loads(source(RUN/"results/source_sha256.json").read_text())
    sizes=sorted({number(r,"circumference") for r in matrix})
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                         "pdf.fonttype":42,"svg.fonttype":"none"})
    summaries=baseline(pairs)
    references,low,high=gain_maps(matrix,sizes)
    outcome_maps(matrix,sizes)
    records,vmax=native_histories(matrix,audited_sources)
    handoff(references,(low,high),vmax)
    links=0
    for target in re.findall(r'\]\(([^)]+)\)',(TOPIC/"COMPONENTS.md").read_text()):
        if not target.endswith("expanded_ring_components.manifest.json"):
            assert (TOPIC/target.replace('%20',' ')).exists(),target
        links+=1
    counts={key:sum(r["status"]==key for r in matrix) for key in sorted({r["status"] for r in matrix})}
    record=dict(generated_at=datetime.now(timezone.utc).isoformat(),artifact_root="../..",
                sha256={str(p.relative_to(TOPIC)):sha(p) for p in ARTIFACTS},source_sha256=SOURCES,
                audit=dict(passed=True,baseline_native_pairs=13,baseline_checkpoints_verified=26,planned_cells=867,
                           coverage=counts,shared_gain_scale=[0,3],normalized_gain_range=[low,high],references=references,
                           native_Holmes_phase_files=6,native_Holmes_field_count=4,exact_phase_continuity=True,
                           exact_complete_C40_basal_identity=True,identical_C40_cue_value_multiset=True,
                           native_mesh_spacing=.05,display_time_window=[0,4000],linear_time_axis=True,
                           shared_unclipped_active_Rac_scale=[0,vmax],local_links_checked=links,
                           raw_sources_modified=False,new_integrations=0,legacy_repair_coefficients_imported=False,
                           common_active_recovery_measurements="pending separate package",global_vector_component="owned separately"))
    path=out("expanded_ring_components.manifest.json");path.write_text(json.dumps(record,indent=2)+"\n")
    print(json.dumps({"topic":str(TOPIC),"artifacts":len(ARTIFACTS),"sources":len(SOURCES),"audit":record["audit"]},indent=2))


if __name__=="__main__":main()
