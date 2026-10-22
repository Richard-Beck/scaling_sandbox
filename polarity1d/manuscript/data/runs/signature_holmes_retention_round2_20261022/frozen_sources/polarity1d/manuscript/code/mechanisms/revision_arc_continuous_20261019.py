"""Forty-eight saved longest-cue trajectories; no integration or RHS scan."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
GROUP = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
OUT = GROUP / "arc_review"
DATA = GROUP / "standard_arclength_sizes/supporting/data"
COLORS = ("#236f92", "#348b68", "#cf942c", "#a4486a")
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes", spring="Spring")


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dest(name, kind="data"):
    path = OUT / name if name.endswith(".png") else OUT / "supporting" / kind / name
    path.parent.mkdir(parents=True, exist_ok=True); return path
def wrapped(x): return np.degrees(np.angle(np.exp(1j * np.radians(x))))
def clean(values): return [None if not np.isfinite(x) else float(x) for x in values]


def main():
    # Keep the first delivered review at the same relative-link depth.
    snapshot = OUT / "index.pre_continuous.snapshot.html"
    if not snapshot.exists():
        shutil.copy2(OUT / "index.html", snapshot)
        for name in ("source_sha256", "delivery_check"):
            original = OUT / (name + ".json")
            if original.exists(): shutil.copy2(original, OUT / (name + ".pre_continuous.snapshot.json"))
    csvpath, tracepath = DATA / "stationary_arclength_summary.csv", DATA / "stationary_arclength_traces.npz"
    assert sha(csvpath) == "9b00954df3a3ff50f65de37ea9f3d82968bab10196b8cdc10a5d054d5329e278"
    assert sha(tracepath) == "bbac92256b738f4b84f95440c716cabfa583f06ee07c52a7cdcafd9eaff5714f"
    frame = pd.read_csv(csvpath); frame = frame[frame.cue_duration == 1000].copy(); assert len(frame) == 48
    curves, points, intervals, metadata = [], [], [], []
    with np.load(tracepath) as z:
        for row in frame.itertuples():
            cid = row.condition_id; keep = z[f"c{cid}_phase"] == 0
            t = z[f"c{cid}_time"][keep]; P = z[f"c{cid}_P_real"][keep] + 1j * z[f"c{cid}_P_imag"][keep]
            physical_C = z[f"c{cid}_C"][keep]
            physical_offset = physical_C * row.physical_distance_at_onset / row.circumference
            R = abs(P); raw = np.degrees(np.angle(P)); strong = R >= .1
            heading = raw.copy(); heading[np.isclose(abs(heading), 180., rtol=0, atol=1e-6)] = 180.
            heading[~strong] = np.nan
            target = row.target_angle_in_P_coordinates_deg
            closure = 1 - abs(wrapped(raw - target)) / abs(wrapped(-target)); closure[~strong] = np.nan
            antipodal = abs(abs(wrapped(target)) - 180.) < 1e-6
            delta = np.degrees(np.angle(P[1:] * np.conj(P[:-1])))
            available = strong[1:] & strong[:-1] & (not antipodal)
            secant = np.where(available, delta / np.diff(t), np.nan)
            with np.load(row.actual_cue_path) as archive:
                fullstate = "states" in archive.files
            meta = dict(condition_id=int(cid), model=row.model, size=float(row.size_parameter),
                initial_arc=float(row.physical_distance_at_onset), target_deg=float(target), overlap=float(row.overlap_fraction_cue_support),
                actual_basal_C=float(row.circumference), full_field_tail_status=row.full_field_tail_status,
                cue_archive_has_full_states=bool(fullstate), antipodal_axis=bool(antipodal),
                angular_rate_withheld=bool(antipodal), weak_saved_frames=int((~strong).sum()),
                saved_frames=len(t), minimum_R=float(R.min()), minimum_saved_interval=float(np.diff(t).min()),
                maximum_saved_interval=float(np.diff(t).max()), maximum_saved_angular_increment_deg=float(abs(delta).max()),
                actual_cue_path=row.actual_cue_path, actual_cue_sha256=row.actual_cue_sha256)
            metadata.append(meta)
            for j in range(len(t)):
                points.append(dict(condition_id=cid, time=float(t[j]), P_real=float(P[j].real), P_imag=float(P[j].imag),
                    R=float(R[j]), raw_phase_deg=float(raw[j]), display_heading_deg=float(heading[j]),
                    angular_error_closed=float(closure[j]), heading_available=bool(strong[j]),
                    physical_C=float(physical_C[j]), prescribed_physical_arc=float(physical_offset[j])))
            for j in range(len(secant)):
                intervals.append(dict(condition_id=cid, t_left=float(t[j]), t_right=float(t[j + 1]), delta_time=float(t[j + 1] - t[j]),
                    wrapped_phase_increment_deg=float(delta[j]), R_left=float(R[j]), R_right=float(R[j + 1]),
                    sampled_angular_secant_deg_per_native_time=float(secant[j]), rate_available=bool(available[j]),
                    unavailable_reason="antipodal_axis_rate_withheld" if antipodal else ("weak_endpoint" if not available[j] else "")))
            curves.append(dict(**meta, time=clean(t), heading=clean(heading), R=clean(R), closure=clean(closure),
                physical_C=clean(physical_C), prescribed_physical_arc=clean(physical_offset),
                raw_phase=clean(raw), rate_left=clean(t[:-1]), rate_right=clean(t[1:]), sampled_rate=clean(secant)))
    assert len(curves) == 48 and sum(c["antipodal_axis"] for c in curves) == 4
    pd.DataFrame(metadata).to_csv(dest("cue_on48_conditions.csv"), index=False)
    pd.DataFrame(points).to_csv(dest("cue_on48_samples.csv"), index=False)
    pd.DataFrame(intervals).to_csv(dest("cue_on48_sampled_secants.csv"), index=False)
    dest("cue_on48_interactive.json").write_text(json.dumps(curves, separators=(",", ":"), allow_nan=False))
    render(curves)
    report = dict(generated_at=datetime.now(timezone.utc).isoformat(), histories=48, models=4, native_cue_duration=1000,
        saved_cue_on_frames=len(points), intervals=len(intervals), antipodal_rate_withheld=4,
        histories_with_weak_saved_frames=sum(c["weak_saved_frames"] > 0 for c in curves),
        weak_saved_frames=sum(c["weak_saved_frames"] for c in curves), full_state_histories=sum(c["cue_archive_has_full_states"] for c in curves),
        retained_primary_only_histories=sum(not c["cue_archive_has_full_states"] for c in curves),
        angular_secants_available=int(sum(np.isfinite(i["sampled_angular_secant_deg_per_native_time"]) for i in intervals)),
        maximum_nonantipodal_saved_phase_increment_deg=max(c["maximum_saved_angular_increment_deg"] for c in curves if not c["antipodal_axis"]),
        method="Actual physicalP/time caches only. Sampled secant=arg(P_right*conj(P_left))/actual_delta_t, available only for strong endpoints and non-antipodal protocols. No nativeRHS scan, integration or inferred inactive field.",
        heading_convention="Raw phase retained; values within1e-6degree of either180degree branch displayed as+180. WeakR<.1 is a gap. No phase unwrapping across a weak interval.",
        time_display="Native time, linear0...0.1 and logarithmic thereafter; duration1000 in each model's native clock.",
        snapshot=str(snapshot), snapshot_sha256=sha(snapshot),
        source_sha256={str(csvpath):sha(csvpath), str(tracepath):sha(tracepath), str(Path(__file__)):sha(__file__)})
    dest("cue_on48_validation.json", "metadata").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))


def render(curves):
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10.5, "pdf.fonttype":42, "svg.fonttype":"none"})
    for model, name in NAMES.items():
        selected = [c for c in curves if c["model"] == model]
        sizes = (1.,2.,3.) if model == "spring" else (20.,30.,40.)
        arcs = sorted(set(c["initial_arc"] for c in selected)); assert len(arcs) == 4
        fig, axes = plt.subplots(2,3,figsize=(13.6,9.3))
        fig.subplots_adjust(left=.065,right=.985,bottom=.18,top=.80,wspace=.23,hspace=.33)
        fig.text(.035,.95,f"{name} · Saved cue-on activity-vector histories",fontsize=18,weight="bold")
        fig.text(.035,.912,"Same four initial physical cue offsets, three native sizes. The cue remains on throughout 0–1000; withdrawal is not shown.",fontsize=11)
        fig.legend(handles=[Line2D([0],[0],color=color,lw=1.8,label=f"d = {d:.5g}") for color,d in zip(COLORS,arcs)],
                   loc="upper left",bbox_to_anchor=(.035,.895),ncol=4,frameon=False,fontsize=11)
        maximum_R = max(max(c["R"]) for c in selected)
        for col,size in enumerate(sizes):
            panel = sorted([c for c in selected if c["size"] == size],key=lambda c:c["initial_arc"])
            for c,color in zip(panel,COLORS):
                axes[0,col].plot(c["time"],c["heading"],color=color,lw=1.7)
                axes[0,col].axhline(c["target_deg"],color=color,ls=":",lw=.7,alpha=.75)
                axes[1,col].plot(c["time"],c["R"],color=color,lw=1.7)
            title=f"C = {size:g}" if model != "spring" else f"C₀ = {size:g}; basal C = {panel[0]['actual_basal_C']:.4f}"
            targets=" / ".join(f"{c['target_deg']:.2f}".rstrip("0").rstrip(".") for c in panel)
            axes[0,col].set_title(title+"\nTargets: "+targets+"°",fontsize=10.5,loc="left")
            axes[0,col].set_ylim(-4,184)
            axes[1,col].set_ylim(0,max(1.15,maximum_R*1.06))
            axes[1,col].axhline(1,color="#bbbbbb",ls=":",lw=.75)
            axes[1,col].axhline(.1,color="#999999",ls="--",lw=.7)
            for ax in axes[:,col]:
                ax.set_xscale("symlog",linthresh=.1);ax.set_xlim(0,1000)
                ax.set_xticks([0,.1,1,10,100,1000],["0","0.1","1","10","100","1000"])
                ax.minorticks_off();ax.spines[["top","right"]].set_visible(False)
            axes[1,col].set_xlabel("Time since cue onset (native units)")
        axes[0,0].set_ylabel("Activity-vector heading (degrees)")
        axes[1,0].set_ylabel("R = |M₁ / M₁(initial)|")
        fig.text(.035,.087,"Dotted heading lines mark the actual targets. R < 0.1 produces heading gaps; no phase is joined across cancellation.",fontsize=10)
        fig.text(.035,.062,"The smallest-size antipodal heading switch is symmetry-constrained and passes through a weak resultant; it is not a measured finite rapid rotation.",fontsize=9.8)
        fig.text(.035,.037,"Native clock: linear 0–0.1, logarithmic thereafter. All initial overlaps retained; exact cue-footprint overlap and sample intervals are exported.",fontsize=9.8)
        for ext in ("png","pdf","svg"):fig.savefig(dest("cue_on_"+model+"."+ext,"vectors"),dpi=210)
        plt.close(fig)
        print(str(dest("cue_on_"+model+".png")),flush=True)


if __name__=="__main__":main()
