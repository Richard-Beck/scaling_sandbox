"""Read-only physical cyclic alignment of N400/N800 periodic basal states."""
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
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
RUNS = {400: ROOT / "data/runs/ring_preparation_pilot_20261017",
        800: ROOT / "data/runs/ring_preparation_pilot_N800_20261017"}
OUT = ROOT / "figures/length_mechanisms_20261008/dev/ring_preparation_mesh_20261017"
SOURCES = {}


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


def write_csv(name, rows):
    path = dev_output_path(OUT, name)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def width(x, values):
    level = values.min() + 0.5 * np.ptp(values)
    inside = np.flatnonzero(values >= level)
    left, right = int(inside[0]), int(inside[-1])
    assert left > 0 and right < len(x) - 1
    assert np.all(np.diff(inside) == 1)
    xleft = x[left - 1] + (level-values[left - 1])/(values[left]-values[left - 1])*(x[left]-x[left - 1])
    xright = x[right] + (level-values[right])/(values[right + 1]-values[right])*(x[right + 1]-x[right])
    return float(xright-xleft)


def main():
    source(Path(__file__))
    source(Path(__file__).with_name("length_dev_output_layout.py"))
    data = {}
    for cells, run in RUNS.items():
        audit = read(run / "results/independent_audit.json")
        assert audit["passed"]
        manifest = read(run / "manifest.json")
        assert manifest["cells"] == cells and manifest["circumference"] == 40
        for job in manifest["jobs"]:
            record = read(run / "workers" / f'{job["slot"]:02d}' / "status.json")
            assert record["status"] == "qualifying_single_site" and len(record["attempts"]) == 1
            assert record["attempts"][0]["duration"] == 10
            path = source(run / record["checkpoint"], record["checkpoint_sha256"])
            with np.load(path) as z:
                state, x = z["state"].copy(), z["x"].copy()
            origin = float(x[np.argmax(state[0])])
            coordinate = (x - origin + 20) % 40 - 20
            order = np.argsort(coordinate)
            data[(job["model"], cells)] = dict(x=coordinate[order], state=state[:,order],
                origin=origin, width=width(coordinate[order],state[0,order]),
                peak=float(state[0].max()), job=job, manifest=manifest)
    names = {"wave_pinning":"Wave-pinning", "goryachev":"Goryachev", "holmes_model3":"Holmes 3"}
    summary, fields, plots = [], [], {}
    for key in names:
        coarse, fine = data[(key,400)], data[(key,800)]
        assert coarse["job"]["prep_peak"] == fine["job"]["prep_peak"]
        assert coarse["job"]["diffusivities"] == fine["job"]["diffusivities"]
        assert coarse["job"]["fixed_field_scales"] == fine["job"]["fixed_field_scales"]
        common = coarse["x"]
        interpolated = np.stack([np.interp(common, fine["x"], values, period=40) for values in fine["state"]])
        scales = np.asarray([coarse["job"]["fixed_field_scales"][name] for name in coarse["job"]["state_names"]])
        difference = coarse["state"] - interpolated
        native_max = np.max(np.abs(difference),axis=1)
        normalized_max = np.max(np.abs(difference)/scales[:,None],axis=1)
        summary.append(dict(model=key, N400_peak=coarse["peak"], N800_peak=fine["peak"],
            peak_change_fraction=fine["peak"]/coarse["peak"]-1,
            N400_half_contrast_width=coarse["width"], N800_half_contrast_width=fine["width"],
            width_change_fraction=fine["width"]/coarse["width"]-1,
            N400_measured_peak_position=coarse["origin"], N800_measured_peak_position=fine["origin"],
            fine_physical_alignment_shift=coarse["origin"]-fine["origin"],
            max_all_field_difference_over_conserved_pool_scale=float(normalized_max.max()),
            max_primary_difference_native=float(native_max[0]),
            max_primary_difference_over_N800_peak=float(native_max[0]/fine["peak"])))
        for field_index,name in enumerate(coarse["job"]["state_names"]):
            for index, position in enumerate(common):
                fields.append(dict(model=key, field=name, arc_from_incumbent_peak_native=float(position),
                    N400_native=float(coarse["state"][field_index,index]),
                    N800_at_common_position_native=float(interpolated[field_index,index]),
                    N400_minus_N800_native=float(difference[field_index,index]),
                    fixed_conserved_pool_scale=float(scales[field_index])))
        plots[key] = dict(coarse=coarse,fine=fine,position=common,
                          maximum_field_error=np.max(np.abs(difference)/scales[:,None],axis=0))
    files = [write_csv("mesh_summary.csv",summary),write_csv("aligned_native_fields.csv",fields)]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,"axes.titlesize":10,
        "axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7,
        "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    fig,axes=plt.subplots(2,3,figsize=(7.8,4.9))
    fig.subplots_adjust(left=.085,right=.98,top=.88,bottom=.12,wspace=.42,hspace=.47)
    for col,(key,label) in enumerate(names.items()):
        plot=plots[key];c,f=plot["coarse"],plot["fine"]
        axes[0,col].plot(c["x"],c["state"][0],color="#0072B2",lw=1.7,label="N = 400")
        axes[0,col].plot(f["x"],f["state"][0],color="#D55E00",lw=1.3,ls="--",label="N = 800")
        axes[0,col].set_title(label)
        axes[1,col].plot(plot["position"],plot["maximum_field_error"],color="#333333",lw=1.3)
        extent=4 if key=="goryachev" else 20
        for row in range(2):
            axes[row,col].set(xlim=(-extent,extent),xticks=[-extent,0,extent])
            axes[row,col].grid(axis="y",alpha=.12,lw=.6)
        axes[1,col].set_xlabel("Arc from incumbent peak\n(native length)")
        if col==0:
            axes[0,col].set_ylabel("Primary activity\n(native concentration)")
            axes[1,col].set_ylabel("Largest field difference\n/ fixed pool scale")
        if key=="goryachev":
            axes[0,col].text(.98,.94,"Site view",transform=axes[0,col].transAxes,ha="right",va="top",fontsize=7,color="#555555")
        axes[1,col].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.53,.99),ncol=2,frameon=False,fontsize=8)
    for extension in ("png","pdf","svg"):
        path=dev_output_path(OUT,"ring_preparation_mesh_comparison."+extension)
        fig.savefig(path,dpi=300,facecolor="white")
        files.append(path)
    plt.close(fig)
    caption=("Preparation mesh sensitivity at C = 40. Every N = 400 and N = 800 run starts from a fresh native homogeneous state, uses the same total-support-3 centred preparation cue for 10 native units (peaks Wave-pinning 0.1, Goryachev 0.1, Holmes 1), and passes the fixed full-field basal-window and holdout checks. Top: final native primary profiles after physical cyclic alignment of the measured primary peak; the Goryachev site is shown over ±4 native length units. Bottom: the largest absolute difference over all native fields, each divided by its conserved-group mean, on the N = 400 physical coordinates after periodic interpolation of the N = 800 representation. The stored finite-volume cell averages are represented at their native cell centres; no checkpoint is altered or amplitude fitted. Peak/width changes and full aligned fields are exported. This supports preparation feasibility on both meshes; a two-grid comparison does not establish continuum convergence or challenger-response accuracy. Ring cue geometry and preparation doses are not asserted equal to the current interval response protocol.")
    path=dev_output_path(OUT,"caption.txt");path.write_text(caption+"\n");files.append(path)
    verification=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,
        qualifying_models_on_each_grid=3,new_integrations_in_renderer=0,challenger_trials=0,
        alignment="Physical arc coordinate relative to each measured primary-peak cell; common circumference40. The resulting origin shifts are recorded explicitly; no amplitude normalization or state mutation.",
        qualification="Both mesh-specific independent saved-field audits pass; unchanged fixed scales and native kinetic/diffusion parameters.",
        interpretation="Observed N400-to-N800 sensitivity only, not a universal grid-adequacy pass/fail rule.",
        visual_inspection="Pending rendered PNG inspection.",summary=summary)
    path=dev_output_path(OUT,"verification.json");path.write_text(json.dumps(verification,indent=2)+"\n");files.append(path)
    manifest=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer="finished",artifact_root="../..",
        source_sha256=SOURCES,sha256={dev_artifact_key(OUT,path):sha(path) for path in files},
        renderer_sha256=sha(Path(__file__)),reproduction="PYTHONDONTWRITEBYTECODE=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python "+str(Path(__file__)))
    dev_output_path(OUT,"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(dict(output=str(OUT),summary=summary),indent=2))


if __name__=="__main__":
    main()
