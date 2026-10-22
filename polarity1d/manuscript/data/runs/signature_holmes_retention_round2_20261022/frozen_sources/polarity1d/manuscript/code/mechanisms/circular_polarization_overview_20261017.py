"""Five-model circular-preparation overview from retained audited profiles only."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/length_mechanisms_20261008/dev/circular_polarization_overview_20261017"
CHEM = ROOT / "data/runs/ring_preparation_pilot_N800_20261017"
SPRING = ROOT / "data/runs/spring_ring_preparation_pilot_20261017"
DEBELLY = ROOT / "data/runs/debelly_planar_loop_preparation_20261017"
SOURCES = {}
PHASES = ("native_start", "cue_withdrawal", "cue_free_endpoint")
PHASE_STYLE = {"native_start":("#777777", ":", "Native start"),
               "cue_withdrawal":("#BD9138", "-", "Cue withdrawal"),
               "cue_free_endpoint":("#235789", "-", "Cue-free endpoint")}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path, expected=None):
    path = Path(path)
    digest = sha(path)
    if expected is not None:
        assert digest == expected, path
    SOURCES[str(path.relative_to(ROOT))] = digest
    return path


def read(path):
    return json.loads(source(path).read_text())


def profile(x, values, circumference, origin, archive, key):
    x, values = np.asarray(x), np.asarray(values)
    assert x.shape == values.shape and np.isfinite(values).all()
    arc = (x-origin+circumference/2) % circumference-circumference/2
    order = np.argsort(arc)
    return dict(arc=arc[order], values=values[order], circumference=float(circumference),
        physical_origin=float(origin), archive=str(Path(archive).relative_to(ROOT)), key=key,
        native_cell_index=order)


def chemical_profiles():
    manifest = read(CHEM / "manifest.json")
    audit = read(CHEM / "results/independent_audit.json")
    assert audit["passed"] and manifest["cells"] == 800
    output = {}
    labels = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "holmes_model3":"Holmes 3"}
    for job in manifest["jobs"]:
        record = read(CHEM / "workers" / f'{job["slot"]:02d}' / "status.json")
        assert record["status"] == "qualifying_single_site"
        prep = next(item for item in record["attempts"][0]["phases"] if item["phase"] == "preparation")
        prep_path = source(CHEM / prep["path"], prep["sha256"])
        with np.load(prep_path) as z:
            assert z["time"][-1] == 10
            x, start, withdrawal = z["x"].copy(), z["states"][0,0].copy(), z["states"][-1,0].copy()
        final_path = source(CHEM / record["checkpoint"],record["checkpoint_sha256"])
        final_meta = read(final_path.with_suffix(".json"))
        with np.load(final_path) as z:
            final = z["state"][0].copy()
            np.testing.assert_array_equal(z["x"],x)
        origin = float(x[np.argmax(final)])
        model = job["model"]
        output[model] = dict(label=labels[model], qualified=True, status="qualified single site",
            geometry="Periodic conservative cell-centred reaction–diffusion ring; native local reactions, diffusivities and pool means unchanged.",
            cells=800, primary="Rac" if model=="holmes_model3" else "u",
            concentration_unit="µM" if model=="holmes_model3" else "native conc.",
            length_unit="µm" if model=="holmes_model3" else "native length",
            arc_convention="Fixed cell-centred ring arclength",
            cue=dict(peak=job["prep_peak"], total_support=3, duration=10),
            basal_time=final_meta["basal_elapsed"], holdout=10000,
            profiles={
                "native_start":profile(x,start,40,origin,prep_path,"states[0,0]"),
                "cue_withdrawal":profile(x,withdrawal,40,origin,prep_path,"states[-1,0]"),
                "cue_free_endpoint":profile(x,final,40,origin,final_path,"state[0]"),
            })
    return output


def spring_profiles():
    manifest = read(SPRING / "manifest.json")
    audit = read(SPRING / "results/audit.json")
    record = read(SPRING / "workers/status.json")
    assert audit["passed"] and record["status"] == "qualifying_single_site"
    input_path = source(SPRING / "inputs/uniform.npz",manifest["input_sha256"])
    prep = next(item for item in record["attempts"][0]["phases"] if item["kind"]=="preparation")
    prep_path = source(SPRING / prep["path"],prep["sha256"])
    final_path = source(SPRING / record["checkpoint"],record["checkpoint_sha256"])
    read(final_path.with_suffix(".json"))
    with np.load(input_path) as z:
        start,q = z["state"].copy(),z["material_positions"].copy()
    with np.load(prep_path) as z:
        withdrawal=z["states"][-1].copy()
        assert z["time"][-1]==10
        np.testing.assert_array_equal(z["states"][0],start)
        np.testing.assert_array_equal(z["material_positions"],q)
    with np.load(final_path) as z:
        final=z["state"].copy()
        np.testing.assert_array_equal(z["material_positions"],q)
    count=len(q)
    incumbent=int(np.argmax(final[:count]))
    profiles={}
    for phase,state,path,key in (("native_start",start,input_path,"state[:200]"),
                                  ("cue_withdrawal",withdrawal,prep_path,"states[-1,:200]"),
                                  ("cue_free_endpoint",final,final_path,"state[:200]")):
        circumference=float(state[-1])
        profiles[phase]=profile(q*circumference,state[:count],circumference,
            float(q[incumbent]*circumference),path,key)
    return dict(label="Spring analogue",qualified=True,status="qualified single site",cells=count,
        primary="Active Rac",concentration_unit="native conc.",length_unit="native length",
        arc_convention="Instantaneous physical arc q C(t); each phase retains its own circumference",
        geometry="Circular global-extension analogue: periodic signalling, uniform material deformation and circumference driven by the spatial mean after the native nonlinear force switch. Native coefficients are retained; native 2D mechanics equivalence is not claimed.",
        cue=dict(peak=4,total_physical_support=0.2,duration=10),basal_time=13000,holdout=10000,
        scipy_warning_retained_in_source_audit=True,profiles=profiles)


def debelly_profiles():
    status_path=DEBELLY/"status.json"
    record=read(status_path)
    original_audit_path=DEBELLY/"independent_audit.json"
    continuation=DEBELLY/"continuation"
    audit_path=continuation/"independent_audit.json"
    if not record.get("complete") or not audit_path.exists():
        raise RuntimeError("deBelly prospective verification is not ready; no five-panel overview has been fabricated.")
    cfg=read(DEBELLY/"inputs.json")
    original_audit=read(original_audit_path)
    audit=read(audit_path)
    verification_status=read(continuation/"status.json")
    assert verification_status["complete"] and verification_status["audited"]
    assert sha(original_audit_path)==audit["original_20k_audit_sha256"]
    assert sha(status_path)==audit["original_status_sha256"]
    identity=read(continuation/"restart_identity_before_evolution.json")
    assert identity["all_exact"] and identity["before_first_evolution"]
    initial_path=source(continuation/"initial_state.npz",identity["start_checkpoint_sha256"])
    with np.load(source(DEBELLY/"half_loop_checkpoint.npz")) as original,np.load(initial_path) as restart:
        for key in identity["max_absolute_errors"]:
            np.testing.assert_array_equal(original[key],restart[key])
    read(DEBELLY/"manifest.json")
    read(continuation/"manifest.json")
    pulse_path=source(DEBELLY/"phases/pulse.npz")
    half_path=source(continuation/"half_loop_checkpoint.npz")
    full_path=source(continuation/"full_loop_checkpoint.npz")
    with np.load(pulse_path) as z:
        x,start,withdrawal=z["x"].copy(),z["field_rac"][0].copy(),z["field_rac"][-1].copy()
        assert z["time"][-1]==cfg["cue_duration"]
    with np.load(half_path) as z:
        final_half=z["rac"].copy()
    with np.load(full_path) as z:
        s,final=z["s"].copy(),z["rac"].copy()
    C=float(cfg["circumference"])
    np.testing.assert_array_equal(s,np.r_[x,C-x[::-1]])
    np.testing.assert_array_equal(final,np.r_[final_half,final_half[::-1]])
    profiles={
        "native_start":profile(s,np.r_[start,start[::-1]],C,0,pulse_path,"field_rac[0] mirrored evenly"),
        "cue_withdrawal":profile(s,np.r_[withdrawal,withdrawal[::-1]],C,0,pulse_path,"field_rac[-1] mirrored evenly"),
        "cue_free_endpoint":profile(s,final,C,0,full_path,"rac"),
    }
    qualified=bool(audit.get("success",False))
    return dict(label="deBelly point-cut loop",qualified=qualified,
        status="qualified after prospective continuation" if qualified else "qualification unresolved after continuation",cells=len(s),
        primary="Rac",concentration_unit="native conc.",length_unit="µm",
        arc_convention="Fixed reference-circle arclength centred on the point cut; not reconstructed nonlinear deformed-contour arclength",
        geometry=cfg["mechanics"]+"; "+cfg["gap_rule"]+". The two mirrored flanks belong to one site around the point cut. This is not a full nonlinear 2D contour simulation.",
        cue=dict(peak=cfg["cue_peak"],half_support=cfg["cue_width_half"],duration=cfg["cue_duration"]),
        basal_time=30000,holdout=10000,prospective_holdout=[20000,30000],
        original_20000_cap_success=original_audit["success"],
        original_20000_cap_preserved_unresolved=True,
        full_native_state_stationary=audit["full_native_state_stationary"],
        cue_free_single_incumbent=audit["cue_free_single_incumbent"],profiles=profiles)


def main():
    source(Path(__file__))
    source(Path(__file__).with_name("length_dev_output_layout.py"))
    models=chemical_profiles()
    models["spring"]=spring_profiles()
    models["debelly"]=debelly_profiles()
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.titlesize":10,
        "axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
        "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig=plt.figure(figsize=(7.8,5.5))
    grid=fig.add_gridspec(2,6,left=.085,right=.975,bottom=.10,top=.86,wspace=1.2,hspace=.66)
    positions=[grid[0,0:2],grid[0,2:4],grid[0,4:6],grid[1,1:3],grid[1,3:5]]
    axes=[fig.add_subplot(position) for position in positions]
    rows=[]
    for ax,(key,model),letter in zip(axes,models.items(),"ABCDE"):
        for phase in PHASES:
            p=model["profiles"][phase]
            color,style,label=PHASE_STYLE[phase]
            ax.plot(p["arc"],p["values"],color=color,ls=style,lw=1.5 if phase!="cue_withdrawal" else 1.25,label=label)
            for index in range(len(p["arc"])):
                rows.append(dict(model=key,phase=phase,native_cell_index=int(p["native_cell_index"][index]),
                    arc_from_incumbent_native=float(p["arc"][index]),primary_value_native=float(p["values"][index]),
                    primary=model["primary"],concentration_unit=model["concentration_unit"],length_unit=model["length_unit"],
                    arc_convention=model["arc_convention"],
                    phase_circumference_native=p["circumference"],physical_origin=p["physical_origin"],
                    source_archive=p["archive"],source_array=p["key"]))
        C=model["profiles"]["cue_free_endpoint"]["circumference"]
        extent=4 if key=="goryachev" else max(p["circumference"] for p in model["profiles"].values())/2
        coordinate_label="Reference arc from cut" if key=="debelly" else "Arc from incumbent"
        ax.set(xlim=(-extent,extent),xticks=[-extent,0,extent],
            xlabel=f'{coordinate_label}\n({model["length_unit"]})',
            ylabel=f'{model["primary"]} ({model["concentration_unit"]})')
        status="Qualified" if model["qualified"] else "Not qualified"
        ax.set_title(model["label"],pad=20)
        ax.text(.5,1.02,f'{status}; C = {C:.4g}',transform=ax.transAxes,ha="center",va="bottom",fontsize=7.3,
            color="#333333" if model["qualified"] else "#A33A22")
        if key=="goryachev":
            ax.text(.97,.93,"Site view",transform=ax.transAxes,ha="right",va="top",fontsize=7,color="#555555")
        ax.grid(axis="y",alpha=.12,lw=.6)
        box=ax.get_position();fig.text(box.x0-.055,box.y1+.052,letter,fontweight="bold",fontsize=11)
    handles=[Line2D([],[],color=style[0],ls=style[1],lw=1.6,label=style[2]) for style in PHASE_STYLE.values()]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.53,.99),ncol=3,frameon=False,fontsize=8)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();clipped=[]
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box=artist.get_window_extent(renderer)
            if box.width and box.height and not (box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.x1 and box.y1<=fig.bbox.y1):
                clipped.append(artist.get_text())
    assert not clipped,clipped
    files=[]
    for extension in ("png","pdf","svg"):
        path=dev_output_path(OUT,"circular_polarization_overview."+extension)
        fig.savefig(path,dpi=300,facecolor="white");files.append(path)
    plt.close(fig)
    path=dev_output_path(OUT,"native_primary_profiles.csv")
    with path.open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    files.append(path)
    count=sum(model["qualified"] for model in models.values())
    caption=(f"Provisional circular-polarization overview: {count} of five model constructions meet their specified cue-free qualification checks. "
        "Raw primary profiles show each native homogeneous start, cue withdrawal and final cue-free endpoint; amplitudes are not normalized across models. "
        "A–C: Wave-pinning, Goryachev and Holmes 3 retain native local reactions, diffusivities and conserved-pool means on conservative periodic rings (C = 40, N = 800). Goryachev is shown within ±4 native length units to resolve its site. "
        "D: Spring uses a declared circular global-extension analogue with periodic signalling and mean native nonlinear force; its circumference evolves, and each profile uses that phase’s actual physical arc. Native 2D mechanics equivalence is not claimed. "
        "E: deBelly uses a declared tangential point-cut planar-loop reduction on a fixed reference circle (C = 62 µm, two mirrored 31-µm half-arcs); the cue-created zero-width cut persists during basal recovery. Its two flanks are one connected front, not two fronts; full nonlinear contour mechanics are not simulated. "
        "All profiles are centred only for display; saved states remain unchanged. Chemical/Spring preparation lasts 10 native units (total cue support 3/0.2; peaks WP 0.1, Goryachev 0.1, Holmes 1, Spring 4); deBelly uses a 60-s cue of peak 3.6 and 6.2-µm support on each flank. "
        "These geometry-specific cues are not asserted dose-matched to the interval experiments or each other. Qualifications include native full-field stationarity and a 10000-unit cue-free holdout on model-specific frozen scales. "
        "deBelly qualifies through a separate prospective 20000–30000-s verification from its exact 20000-s state; the original capped audit remains unresolved and unchanged. "
        "This is first-milestone preparation evidence, not a selected main figure, a universal geometry result, length-scaling law, challenger/takeover test or migration-persistence claim.")
    path=dev_output_path(OUT,"caption.txt");path.write_text(caption+"\n");files.append(path)
    summary={key:{name:value for name,value in model.items() if name!="profiles"} for key,model in models.items()}
    path=dev_output_path(OUT,"qualification_and_geometry.json")
    path.write_text(json.dumps(dict(generated_at=datetime.now(timezone.utc).isoformat(),qualified=count,total=5,
        models=summary,source_array_checks=True,all_profile_values_raw=True,text_within_page_bounds=True,
        visual_inspection="Pending rendered PNG inspection."),indent=2)+"\n");files.append(path)
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer="finished",artifact_root="../..",
        status="Provisional first-milestone overview",new_integrations=0,challenger_trials=0,
        source_sha256=SOURCES,sha256={dev_artifact_key(OUT,path):sha(path) for path in files},
        renderer_sha256=sha(Path(__file__)),
        reproduction="PYTHONDONTWRITEBYTECODE=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python "+str(Path(__file__)))
    dev_output_path(OUT,"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(dict(output=str(OUT),qualified=count,total=5,models={key:model["status"] for key,model in models.items()}),indent=2))


if __name__=="__main__":
    main()
