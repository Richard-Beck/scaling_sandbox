"""Map accepted fixed-initial-arc cues against operational native basal edges.

Read-only compact-source derivation; no native integration or changed protocol.
"""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter, defaultdict
import csv
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GROUP = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
SOURCE = GROUP / "standard_arclength_sizes/supporting/data"
OUT = GROUP / "arc_review"
DATA = OUT / "supporting/data"
META = OUT / "supporting/metadata"
VECTORS = OUT / "supporting/vectors"
MODELS = [("wave_pinning", "Wave-pinning"), ("goryachev", "Goryachev"),
          ("holmes_model3", "Holmes"), ("spring", "Spring")]
COLORS = ["#236f92", "#348b68", "#cf942c", "#a4486a"]


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def edge(x, u, level):
    """Interpolate the connected threshold region containing the declared axis."""
    y = (u - u.min()) / np.ptp(u)
    k = int(np.argmin(np.abs(x)))
    assert abs(x[k]) < 1e-10 and y[k] >= level
    left = right = k
    while left > 0 and y[left - 1] >= level: left -= 1
    while right + 1 < len(x) and y[right + 1] >= level: right += 1
    assert 0 < left <= right < len(x) - 1
    xl = x[left - 1] + (level - y[left - 1]) * (x[left] - x[left - 1]) / (y[left] - y[left - 1])
    xr = x[right] + (level - y[right]) * (x[right + 1] - x[right]) / (y[right + 1] - y[right])
    return float(xl), float(xr)


def overlap(core, d, width, circumference):
    lo, hi = core
    amount = sum(max(0., min(hi, d + width / 2 + j * circumference)
                     - max(lo, d - width / 2 + j * circumference)) for j in (-1, 0, 1))
    fraction = float(np.clip(amount / width, 0, 1))
    status = "inside" if fraction >= 1 - 1e-10 else "outside" if fraction <= 1e-10 else "straddling"
    return amount, fraction, status


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    for p in (DATA, META, VECTORS): p.mkdir(parents=True, exist_ok=True)
    context_path = SOURCE / "stationary_arclength_initial_context.csv"
    profile_path = SOURCE / "stationary_arclength_initial_profiles_and_cues.csv"
    context = list(csv.DictReader(context_path.open()))
    raw = defaultdict(list)
    for r in csv.DictReader(profile_path.open()):
        raw[(r["model"], float(r["size_parameter"]), r["placement_label"])].append(r)
    assert len(context) == len(raw) == 48
    basals, placements, profiles = {}, [], []
    max_edge_error = max_overlap_error = 0.
    for r in context:
        model, size, label = r["model"], float(r["size_parameter"]), r["placement_label"]
        c, d, w = map(float, (r["circumference"], r["physical_distance_at_onset"], r["full_physical_support"]))
        p = raw[(model, size, label)]
        angle = np.array([float(v["normal_angle_from_preparation_deg"]) for v in p])
        x = ((angle + 180) % 360 - 180) * c / 360
        order = np.argsort(x); x = x[order]
        u = np.array([float(v["basal_active"]) for v in p])[order]
        key = model, size
        if key not in basals:
            edges = {level: edge(x, u, level) for level in (.1, .5)}
            basals[key] = dict(x=x, u=u, edges=edges, context=r)
            for xx, uu in zip(x, u):
                profiles.append(dict(model=model, size_parameter=size, circumference=c,
                                     signed_physical_position_from_preparation=float(xx),
                                     basal_active=float(uu), contrast_normalized_active=float((uu-u.min())/np.ptp(u))))
        else:
            assert np.array_equal(x, basals[key]["x"]) and np.array_equal(u, basals[key]["u"])
        b = basals[key]
        max_edge_error = max(max_edge_error, abs(b["edges"][.5][0] - float(r["original_core_left_physical"])),
                             abs(b["edges"][.5][1] - float(r["original_core_right_physical"])))
        out = dict(model=model, size_parameter=size, size_label=r["size_label"], circumference=c,
                   placement_label=label, initial_arc_grid_index=int(r["initial_arc_grid_index"]),
                   initial_physical_centre_offset=d, actual_target_angle_in_P_coordinates_deg=float(r["target_angle_in_P_coordinates_deg"]),
                   full_physical_cue_support=w, initial_angular_cue_support_deg=float(r["initial_angular_support_deg"]),
                   native_cue_peak=float(r["requested_peak"]), sampled_cue_peak=float(r["initial_sampled_cue_peak"]),
                   prescribed_input_coefficient_integral=float(r["initial_cue_integral"]),
                   initial_unit_peak_cue_excess_fraction_total=float(r["initial_unit_peak_cue_excess_fraction_total"]),
                   initial_cue_weighted_active=float(r["initial_cue_weighted_active"]),
                   initial_cue_weighted_inactive=float(r["initial_cue_weighted_inactive"]),
                   basal_input=r["input"], basal_input_sha256=r["input_sha256"], cue_input=r["cue_input"], cue_input_sha256=r["cue_input_sha256"])
        for level, suffix in ((.5, "50"), (.1, "10")):
            lo, hi = b["edges"][level]
            amount, fraction, status = overlap((lo, hi), d, w, c)
            out.update({f"initial_edge_left_{suffix}": lo, f"initial_edge_right_{suffix}": hi,
                        f"initial_core_width_{suffix}": hi-lo,
                        f"signed_centre_distance_from_positive_edge_{suffix}": d-hi,
                        f"centre_relative_to_core_{suffix}": "inside" if lo-1e-10 <= d <= hi+1e-10 else "outside",
                        f"footprint_relative_to_core_{suffix}": status,
                        f"cue_support_intersection_{suffix}": amount,
                        f"fraction_cue_support_inside_core_{suffix}": fraction})
        max_overlap_error = max(max_overlap_error, abs(out["fraction_cue_support_inside_core_50"] - float(r["overlap_fraction_cue_support"])))
        placements.append(out)
    assert len(basals) == 12 and max_edge_error < 1e-10 and max_overlap_error < 1e-10
    placements.sort(key=lambda r: ([m for m, _ in MODELS].index(r["model"]), r["size_parameter"], r["initial_arc_grid_index"]))
    write_csv(DATA / "initial_placement48.csv", placements)
    write_csv(DATA / "initial_basal12_profiles.csv", profiles)
    baseline = []
    for (model, size), b in basals.items():
        row = dict(model=model, size_parameter=size, size_label=b["context"]["size_label"], circumference=float(b["context"]["circumference"]),
                   native_active_min=float(b["u"].min()), native_active_max=float(b["u"].max()), native_active_contrast=float(np.ptp(b["u"])))
        for level, suffix in ((.5,"50"),(.1,"10")):
            row.update({f"edge_left_{suffix}":b["edges"][level][0], f"edge_right_{suffix}":b["edges"][level][1],
                        f"width_{suffix}":b["edges"][level][1]-b["edges"][level][0]})
        baseline.append(row)
    write_csv(DATA / "initial_basal12_edges.csv", baseline)
    plt.rcParams.update({"font.size":10, "axes.titlesize":11, "svg.fonttype":"none", "pdf.fonttype":42})
    fig, axes = plt.subplots(4, 3, figsize=(14.5, 14), sharey=True)
    for ri, (model, title) in enumerate(MODELS):
        sizes = [1,2,3] if model == "spring" else [20,30,40]
        half = max(float(basals[(model, s)]["context"]["circumference"]) for s in sizes) / 2
        for ci, size in enumerate(sizes):
            ax = axes[ri, ci]; b = basals[(model,size)]
            c = float(b["context"]["circumference"])
            for a, z in ((-half, -c/2), (c/2, half)):
                if z > a+1e-12: ax.axvspan(a,z,color="#f1f3f5",zorder=0)
            ax.axvspan(*b["edges"][.5], color="#e9f0f5", zorder=0)
            for lv, style, color in ((.5,"--","#657785"),(.1,":","#9aa5ad")):
                for xx in b["edges"][lv]: ax.axvline(xx, ls=style, color=color, lw=1, zorder=1)
                ax.axhline(lv, color=color, ls=style, lw=.7, alpha=.55)
            ax.plot(b["x"], (b["u"]-b["u"].min())/np.ptp(b["u"]), color="#182d38", lw=1.8, zorder=2)
            local = [r for r in placements if r["model"] == model and r["size_parameter"] == size]
            for j, r in enumerate(local):
                yy = -.18-j*.17; d=r["initial_physical_centre_offset"]; w=r["full_physical_cue_support"]
                for shift in (-c,0,c):
                    a=max(-c/2,d-w/2+shift); z=min(c/2,d+w/2+shift)
                    if z>a: ax.plot([a,z],[yy,yy],color=COLORS[j],lw=5,solid_capstyle="butt",zorder=4)
                ax.plot([d],[yy],marker="o",ms=4,mec="white",mew=.6,color=COLORS[j],zorder=5)
                label=f'd={d:.4g} · {r["footprint_relative_to_core_50"]}'
                ax.text(-half+.035*half, yy, label, ha="left",va="center",fontsize=8.1,color=COLORS[j],
                        bbox=dict(fc="white",ec="none",alpha=.88,pad=.3), zorder=6)
            w50=b["edges"][.5][1]-b["edges"][.5][0]; w10=b["edges"][.1][1]-b["edges"][.1][0]
            size_label=f'C₀={size}; Cᵦ={c:.4f}' if model=="spring" else f'C={size}'
            ax.set_title(f'{title} · {size_label}\nw₅₀={w50:.4g}; w₁₀={w10:.4g}',loc="left",pad=6)
            ax.set_xlim(-half,half); ax.set_ylim(-.83,1.1); ax.set_yticks([0,.5,1]); ax.set_yticklabels(["0",".5","1"])
            ax.spines[["top","right"]].set_visible(False); ax.spines["bottom"].set_color("#97a6ad")
            ax.tick_params(labelsize=9)
            if ci==0: ax.set_ylabel("Activity above minimum / basal contrast\nCue supports below",fontsize=10)
            if ri==3: ax.set_xlabel("Signed initial physical arc from prepared axis",fontsize=10)
    fig.subplots_adjust(left=.075,right=.985,top=.86,bottom=.13,wspace=.19,hspace=.47)
    fig.suptitle("Where the existing cues fall on the initial incumbent",x=.075,y=.975,ha="left",fontsize=20,fontweight="bold")
    fig.text(.075,.938,"Same four initial physical offsets at every size; native basals and cue widths remain unchanged.",fontsize=12)
    fig.text(.075,.915,"Black: [u − min(u)] / [max(u) − min(u)]. Dashed edges: 50% contrast; dotted: 10%.",fontsize=10.5)
    fig.text(.075,.067,"Coloured bars are the full prescribed cue support; dots are continuous cue centres. Bar labels classify support at 50% contrast.",fontsize=10)
    fig.text(.075,.047,"Physical x scale is shared across sizes within each model row. Grey outer areas are beyond the smaller ring; antipodal support wraps.",fontsize=10)
    fig.text(.075,.027,"These operational edges are readouts, not biological boundaries. Initial cue overlap is retained; it does not establish a distinct challenger or site ancestry.",fontsize=10)
    fig.savefig(OUT / "initial_basal_placement_map.png",dpi=180,facecolor="white")
    for suffix in ("pdf","svg"): fig.savefig(VECTORS / ("initial_basal_placement_map."+suffix),facecolor="white")
    plt.close(fig)
    generated=datetime.now(timezone.utc).isoformat()
    metrics=dict(generated_at=generated,status="passed",native_work="none; accepted compact native-profile/context tables only",
                 basals=12,placements=48,condition_rows=240,
                 edge_definition="connected initial region containing the declared preparation axis, at min(u)+level*(max(u)-min(u)); linear interpolation between native cell centres",
                 signed_distance_definition="prescribed positive initial centre arc minus positive initial operational edge; negative is inside",
                 sensitivity="10% of initial contrast, separate from the main 50% criterion; not a statistical confidence interval",
                 footprint_definition="length of full physical cue support inside periodic operational region, divided by full cue support; inside/straddling/outside",
                 maximum_halfcontrast_edge_replay_error=max_edge_error, maximum_halfcontrast_overlap_replay_error=max_overlap_error,
                 halfcontrast_centre_counts=dict(Counter(r["centre_relative_to_core_50"] for r in placements)),
                 halfcontrast_footprint_counts=dict(Counter(r["footprint_relative_to_core_50"] for r in placements)),
                 tenpercent_centre_counts=dict(Counter(r["centre_relative_to_core_10"] for r in placements)),
                 tenpercent_footprint_counts=dict(Counter(r["footprint_relative_to_core_10"] for r in placements)),
                 centre_class_changes=sum(r["centre_relative_to_core_50"]!=r["centre_relative_to_core_10"] for r in placements),
                 footprint_class_changes=sum(r["footprint_relative_to_core_50"]!=r["footprint_relative_to_core_10"] for r in placements),
                 sources={str(p):sha(p) for p in (context_path,profile_path,Path(__file__))})
    (META / "initial_placement_map_validation.json").write_text(json.dumps(metrics,indent=2)+"\n")
    artifacts=[DATA / "initial_placement48.csv", DATA / "initial_basal12_edges.csv",DATA / "initial_basal12_profiles.csv",
               OUT / "initial_basal_placement_map.png", VECTORS / "initial_basal_placement_map.pdf",VECTORS / "initial_basal_placement_map.svg",
               META / "initial_placement_map_validation.json"]
    (META / "initial_placement_map_artifacts.json").write_text(json.dumps(dict(generated_at=generated,artifacts={str(p):sha(p) for p in artifacts}),indent=2)+"\n")
    print(json.dumps(metrics,indent=2)); print(OUT / "initial_basal_placement_map.png")


if __name__ == "__main__": main()
