"""One endpoint-vector-grid alternative from the already validated 240 histories."""
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
from matplotlib.patches import Rectangle,Circle
from matplotlib.colors import LinearSegmentedColormap,Normalize
import revision_standard_angles_sizes_20261019 as campaign

ROOT=campaign.ROOT
REPO=ROOT.parents[1]
OUT=campaign.OUT/"vector_grid"
CSV=campaign.OUT/"supporting/data/stationary_sizes_summary.csv"
TRACES=campaign.OUT/"supporting/data/stationary_sizes_traces.npz"
EXPECTED_CSV="761bf294ad1056620724d2fc8ddc9ca9c2bcd10beed9e8dd2fb3071152bde0d8"
GLYPH_SCALE=.28
ROW_PITCH=.7
R_GATE=.1
MODELS=("wave_pinning","goryachev","holmes_model3","spring")
NAMES={"wave_pinning":"Wave-pinning","goryachev":"Goryachev","holmes_model3":"Holmes","spring":"Spring"}
ANGLES=("pi","pi_2","pi_4","pi_8")
SOURCE={}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source(path):path=Path(path);SOURCE[str(path)]=sha(path);return path
def dest(name,group="data"):
    parent=OUT if name.endswith(".png") else OUT/"supporting"/group
    parent.mkdir(parents=True,exist_ok=True);return parent/name
def write(name,z):dest(name,"metadata").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")


def main():
    assert sha(CSV)==EXPECTED_CSV,"The authoritative table changed; review provenance before rendering"
    frame=pd.read_csv(source(CSV));assert len(frame)==240 and frame.status.eq("complete").all()
    values=[];error=0.;maximum_heading_roundoff=0.
    with np.load(source(TRACES)) as z:
        for row in frame.itertuples():
            cid=row.condition_id;P=complex(z[f"c{cid}_P_real"][-1],z[f"c{cid}_P_imag"][-1]);R=abs(P)
            error=max(error,abs(R-row.R_final));assert abs(R-row.R_final)<1e-12
            weak=R<R_GATE
            if not weak:
                delta=float(np.rad2deg(np.angle(np.exp(1j*(np.angle(P)-np.deg2rad(row.angle_final_deg))))))
                # The accepted table canonicalizes phases within 1e-6° of ±180°.
                assert abs(delta)<1.000001e-6
                maximum_heading_roundoff=max(maximum_heading_roundoff,abs(delta))
            model_row=MODELS.index(row.model);sizes=(1.,2.,3.) if row.model=="spring" else (20.,30.,40.)
            size_column=sizes.index(row.size_parameter);angle_row=ANGLES.index(row.angle_label)
            times=(.1,1.,10.,100.,1000.) if row.model=="spring" else (1.,10.,100.,300.,1000.)
            duration_column=times.index(row.cue_duration)
            fraction=None if weak else float(row.fraction_requested_turn_final)
            values.append(dict(condition_id=cid,model=row.model,model_row=model_row,size_parameter=row.size_parameter,size_label=row.size_label,
                basal_circumference=row.circumference,resting_circumference=row.resting_circumference if pd.notna(row.resting_circumference) else None,
                size_column=size_column,angle_label=row.angle_label,requested_angle_deg=row.requested_angle_deg,angle_row=angle_row,cue_duration=row.cue_duration,
                duration_column=duration_column,origin_x=float(duration_column),origin_y=float(ROW_PITCH*(3-angle_row)),
                P_real=float(P.real),P_imag=float(P.imag),R=R,display_angle_deg=None if weak else float(row.angle_final_deg),
                retained_fraction_requested_turn=fraction,background_fraction_clipped=None if weak else float(np.clip(fraction,0,1)),
                glyph_kind="unoriented_magnitude_circle" if weak else "directed_vector",
                glyph_scale=GLYPH_SCALE,arrow_delta_x=None if weak else float(GLYPH_SCALE*P.real),arrow_delta_y=None if weak else float(GLYPH_SCALE*P.imag),
                magnitude_circle_radius=float(GLYPH_SCALE*R) if weak else None,
                initial_guide_delta_x=GLYPH_SCALE,initial_guide_delta_y=0.,
                full_field_tail_status=row.full_field_tail_status,field_tail_failed=row.full_field_tail_status=="fails_finite_tail_gate",
                field_tail_unavailable=row.full_field_tail_status=="unavailable_from_retained_primary_only_history",
                initial_overlap=bool(row.initial_overlap),overlap_fraction_cue_support=row.overlap_fraction_cue_support,
                initial_angular_support_deg=row.initial_angular_support_deg,full_physical_support=row.full_physical_support,
                target_angle_in_P_coordinates_deg=row.target_angle_in_P_coordinates_deg,
                actual_cue_path=row.actual_cue_path,actual_cue_sha256=row.actual_cue_sha256,
                actual_recovery_path=row.actual_recovery_path,actual_recovery_sha256=row.actual_recovery_sha256))
    glyphs=pd.DataFrame(values);glyphs.to_csv(dest("endpoint_vector_glyphs.csv"),index=False)
    render(glyphs);docs(glyphs,error,maximum_heading_roundoff)


def render(glyphs):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.labelsize":9,"pdf.fonttype":42,"svg.fonttype":"none"})
    cmap=LinearSegmentedColormap.from_list("retained_fraction",["#ffffff","#d7e9e7","#8dbdb7"]);norm=Normalize(0,1)
    fig,axes=plt.subplots(4,3,figsize=(12.1,12.1))
    fig.subplots_adjust(left=.067,right=.982,bottom=.105,top=.815,wspace=.23,hspace=.35)
    fig.suptitle("One stationary-cue experiment · All 240 retained activity vectors",fontsize=16,y=.987)
    fig.text(.04,.952,"Arrows show actual final P = M₁/M₁(initial), after 13,000 cue-free units. Pale rightward guides mark the initial vector P = 1.",fontsize=10)
    fig.text(.04,.929,"Each panel: four cue-centre angles × five tested durations. Columns are duration categories; native clocks are not calibrated across models.",fontsize=9.5)
    fig.text(.04,.905,"Weak R < 0.1: an unoriented circle has radius proportional to R; its measured value is printed. □ failed field-tail gate; △ unavailable full-field tail.",fontsize=9.2)
    cax=fig.add_axes([.20,.867,.60,.012]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation="horizontal",ticks=[0,.25,.5,.75,1])
    cb.ax.tick_params(labelsize=8,pad=2);cb.set_label("Retained fraction of requested turn (colour clipped to 0–1; heading shown only for R ≥ 0.1)",fontsize=8.8,labelpad=3)
    for model in MODELS:
        mi=MODELS.index(model);sizes=(1.,2.,3.) if model=="spring" else (20.,30.,40.)
        for si,size in enumerate(sizes):
            ax=axes[mi,si];panel=glyphs[(glyphs.model==model)&(glyphs.size_parameter==size)]
            assert len(panel)==20
            for row in panel.itertuples():
                x,y=row.origin_x,row.origin_y
                face="#f0f0f0" if row.glyph_kind=="unoriented_magnitude_circle" else cmap(norm(row.background_fraction_clipped))
                ax.add_patch(Rectangle((x-.48,y-.34),.96,.68,facecolor=face,edgecolor="#dddddd",lw=.4,zorder=0))
                ax.plot([x,x+GLYPH_SCALE],[y,y],color="#a8a8a8",lw=.65,zorder=1)
                ax.plot(x,y,".",color="#aaaaaa",ms=1.7,zorder=2)
                if row.glyph_kind=="directed_vector":
                    ax.quiver([x],[y],[row.arrow_delta_x],[row.arrow_delta_y],angles="xy",scale_units="xy",scale=1.,units="xy",
                        width=.014,headwidth=3.5,headlength=3.5,headaxislength=3.3,minlength=0.,pivot="tail",color="#202020",zorder=4)
                else:
                    ax.add_patch(Circle((x,y),row.magnitude_circle_radius,facecolor="#555555",edgecolor="#333333",lw=.5,zorder=4))
                    ax.text(x,y-.21,f"R={row.R:.4f}",ha="center",va="center",fontsize=6.7,color="#333333")
                if row.field_tail_failed:ax.plot(x+.33,y+.22,marker="s",mfc="none",mec="#a55545",mew=.95,ms=4,zorder=6)
                if row.field_tail_unavailable:ax.plot(x+.33,y+.22,marker="^",mfc="none",mec="#777777",mew=.8,ms=4,zorder=6)
            labels=[]
            for label in ANGLES:
                row=panel[panel.angle_label==label].iloc[0]
                labels.append(f"{row.requested_angle_deg:g}°"+("*" if row.initial_overlap else ""))
            times=(.1,1.,10.,100.,1000.) if model=="spring" else (1.,10.,100.,300.,1000.)
            ax.set(xlim=(-.5,4.5),ylim=(-.35,2.45),xticks=range(5),yticks=[ROW_PITCH*y for y in [3,2,1,0]])
            ax.set_xticklabels([f"{t:g}" for t in times],fontsize=8);ax.set_yticklabels(labels,fontsize=8)
            ax.tick_params(length=0,pad=3);ax.set_aspect("equal",adjustable="box")
            for spine in ax.spines.values():spine.set_visible(False)
            title=f"{NAMES[model]} · C = {size:g}" if model!="spring" else f"Spring · C₀ = {size:g}; Cb = {panel.basal_circumference.iloc[0]:.3f}"
            ax.set_title(title,loc="left",fontsize=9.8,pad=6)
            if si==0:ax.set_ylabel("Prescribed cue-centre angle",fontsize=8.7,labelpad=5)
    fig.text(.5,.067,"Tested cue duration (native units; categories)",ha="center",fontsize=10)
    fig.text(.04,.042,"* Cue support intersects the initial half-contrast core: all 125 overlap cases are retained. Physical support stays 3 for chemical models and 0.2 for Spring.",fontsize=9)
    fig.text(.04,.019,"Arrows and circles use one global magnitude scale and equal data-coordinate scales. This is an endpoint-vector view, with no success threshold or minimum-duration claim.",fontsize=8.8)
    for ext in ("png","pdf","svg"):fig.savefig(dest("stationary_endpoint_vector_grid."+ext,"vectors"),dpi=230,bbox_inches="tight")
    plt.close(fig)


def docs(glyphs,error,maximum_heading_roundoff):
    generated=datetime.now(timezone.utc).isoformat();weak=glyphs.glyph_kind.eq("unoriented_magnitude_circle")
    caption="Final activity-vector grid for the complete stationary-cue experiment. Four model rows and three native-size columns contain all 240 requested conditions; each panel has four prescribed cue-centre-angle rows and five tested native-duration columns. Duration positions are categories rather than uniformly elapsed time. Each directed glyph is the measured final normalized complex P = M₁(t)/M₁(0), drawn from the cell origin with one global display factor of 0.28 and equal x/y coordinate scales. Its angle and length therefore encode phase and R = |P| directly. A pale rightward segment marks initial P = 1. When R < 0.1 the panel draws no heading: an unoriented disk has radius 0.28R and prints R. Background colour encodes retained angle divided by the requested angle for strong vectors, clipped only for the displayed 0–1 colour range; raw fractions remain exported. A corner square marks a failed independent full-field tail gate, and a triangle marks unavailable full-field tail qualification. Asterisks on angle rows mark initial cue-support/half-contrast-core intersection; all 125 overlap cases remain. Chemical sizes are fixed C = 20/30/40; Spring uses resting C₀ = 1/2/3 with actual basal Cb shown. All endpoints follow 13,000 cue-free native time units. Fixed physical widths, native peaks and clocks remain model-specific, so this is not an isolated causal geometry comparison. No binary success definition, refined minimum duration, site ancestry or cell movement is inferred. This is one review candidate from the same experiment; the existing response atlas remains unchanged."
    dest("CAPTION.md","metadata").write_text(caption+"\n")
    manifest=dict(generated_at=generated,role="One alternative visual candidate; no additional native analysis or protocol",requested=240,displayed=240,
        source_csv=str(CSV),source_csv_sha256=EXPECTED_CSV,source_traces=str(TRACES),source_traces_sha256=sha(TRACES),
        renderer=str(Path(__file__)),renderer_sha256=sha(__file__),glyph_table=str(dest("endpoint_vector_glyphs.csv")),glyph_table_sha256=sha(dest("endpoint_vector_glyphs.csv")),
        data_mapping=dict(panel_rows=list(MODELS),panel_size_columns=dict(chemical=[20,30,40],spring_resting=[1,2,3]),
            angle_rows=[180,90,45,22.5],duration_categories=dict(chemical=[1,10,100,300,1000],spring=[.1,1,10,100,1000]),
            origin_x="duration category index 0…4",origin_y="0.7 times (3 minus angle-row index 0…3)"),
        glyph_mapping=dict(display_factor=GLYPH_SCALE,strong_R_minimum=R_GATE,arrow_delta="0.28 times actual final (P_real,P_imag)",
            weak_direction="not displayed",weak_disk_radius="0.28 times measured R",initial_guide=[GLYPH_SCALE,0.],background="phase/requested turn, colour clipped to 0…1"),
        quiver=dict(angles="xy",scale_units="xy",scale=1.,units="xy",width=.014,headwidth=3.5,headlength=3.5,headaxislength=3.3,minlength=0.,pivot="tail"),
        axes=dict(aspect="equal",x_limits=[-.5,4.5],y_limits=[-.35,2.45]),
        weak_R_cells=int(weak.sum()),field_tail_failures=int(glyphs.field_tail_failed.sum()),field_tail_unavailable=int(glyphs.field_tail_unavailable.sum()),
        initial_overlap_conditions=int(glyphs.initial_overlap.sum()),maximum_R_source_replay_error=float(error),maximum_heading_difference_due_to_table_180degree_canonicalization=maximum_heading_roundoff,
        all_directed_glyphs_use_actual_P=True,all_weak_directions_masked=True,previous_atlas_preserved=True)
    write("renderer_manifest.json",manifest)
    (OUT/"README.md").write_text("""# Endpoint-vector grid review candidate

One systematic view of the same complete 240-condition stationary-cue experiment. This is a visual alternative for discussion, not a settled manuscript layout.

- [Vector grid](stationary_endpoint_vector_grid.png) · [PDF](supporting/vectors/stationary_endpoint_vector_grid.pdf) · [SVG](supporting/vectors/stationary_endpoint_vector_grid.svg)
- [Exact 240 glyph values and display flags](supporting/data/endpoint_vector_glyphs.csv)
- [Caption](supporting/metadata/CAPTION.md)
- [Source and renderer mapping manifest](supporting/metadata/renderer_manifest.json)

Arrow direction and length are the actual final P and R. Weak-vector cells show measured magnitude without assigning a heading. Independent field-tail flags and initial overlap remain separate. The data and existing atlases are unchanged; no new simulations or diagnostic sweep was performed.
""")
    write("content_checkpoint.json",dict(title="All 240 endpoint activity vectors",role="Single alternative visual candidate from the same stationary experiment",assets=[dict(path=str((OUT/"stationary_endpoint_vector_grid.png").relative_to(REPO)),caption="Actual final P direction and magnitude across all requested sizes, placements and durations.")],legend=[caption],
        sources=[dict(label="Exact displayed glyph values",path=str(dest("endpoint_vector_glyphs.csv").relative_to(REPO))),dict(label="Frozen source and rendering map",path=str(dest("renderer_manifest.json","metadata").relative_to(REPO)))]))
    source(__file__);write("source_sha256.json",SOURCE)
    print(json.dumps(dict(passed=True,displayed=240,weak_R_cells=int(weak.sum()),field_tail_failures=int(glyphs.field_tail_failed.sum()),field_tail_unavailable=int(glyphs.field_tail_unavailable.sum()),
        overlap_conditions=int(glyphs.initial_overlap.sum()),maximum_R_replay_error=error,output=str(OUT)),indent=2))


if __name__=="__main__":main()
