"""Timestamped exploratory broad-strength review from collected native data only.

No integration, protocol change or main-figure replacement. Each invocation saves
the exact collected CSV/status snapshot used by its exports.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import io
import json
import os
import time

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[4]
MS = ROOT / "polarity1d/manuscript"
RUN = MS / "data/runs/revision_noise_strength_broad_20261019"
OUT = MS / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/broad_strength_campaign_20261019/review"
MODELS = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes 3", spring_ring="Spring ring")
GRID = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, .1, .3, 1.])
COLORS = ["#2676a6", "#28866b", "#c77722"]


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def savejson(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def snapshot(run, out):
    partial = run / "results/partial"
    for _ in range(30):
        before = (partial / "status.json").read_bytes()
        csvpath = partial / "condition_window_summaries.csv"
        stat = csvpath.stat()
        data = csvpath.read_bytes()
        after = (partial / "status.json").read_bytes()
        # Collector writes the CSV before its corresponding generation receipt.
        if before == after and stat.st_mtime_ns <= (partial / "status.json").stat().st_mtime_ns and stat.st_mtime_ns == csvpath.stat().st_mtime_ns:
            break
        time.sleep(1)
    else:
        raise RuntimeError("Could not acquire a consistent atomic collector snapshot")
    (out / "supporting/data/condition_window_snapshot.csv").write_bytes(data)
    (out / "supporting/metadata/generation_snapshot.json").write_bytes(after)
    return pd.read_csv(io.BytesIO(data)), json.loads(after), hashlib.sha256(data).hexdigest()


def marker(ax, x, y, color, available, target, size=5):
    if not np.isfinite(y) or available == 0:
        return
    symbol = "s" if target == 12 else "o"
    ax.plot(x, y, marker=symbol, markersize=size, markeredgecolor=color,
            markerfacecolor=color if available == target else "white", linestyle="none", zorder=5)


def main_sheet(df, status):
    fig, axes = plt.subplots(4, 3, figsize=(13.8, 14.6))
    fig.subplots_adjust(left=.075, right=.98, bottom=.12, top=.865, hspace=.55, wspace=.28)
    fig.suptitle("Existing-law strength sweep: exploratory organization readouts", x=.075, ha="left", y=.982, fontsize=19)
    fig.text(.075, .954, "Positive input mean and fluctuation amplitude scale together; each model keeps its native law and three original sizes.", fontsize=11)
    obs = status["data_observed_at"].replace("T", " ")[:19]
    fig.text(.075, .934, f"Native status observed {obs} UTC: {status['producer']} · {status['new_completed_outputs']}/256 new outputs · {status['failures']} failures", fontsize=10, color="#4d5b64")
    fig.text(.075, .915, "New levels: two-seed target (exploratory). Original 0.001/0.01: twelve seeds. Whiskers show observed seed range, not confidence intervals.", fontsize=10)
    phasekey = [Line2D([], [], color="#32434a", lw=1.8, label="solid: raw regions / noisy R / amount A"),
                Line2D([], [], color="#32434a", lw=1.8, ls="--", label="dashed: qualified regions / mean-control R / coherence Q")]
    fig.legend(handles=phasekey, loc="upper left", bbox_to_anchor=(.068, .902), ncol=2, frameon=False, fontsize=10, handlelength=3.3)
    rowdata = df[(df.window == "full")]
    for ri, (model, name) in enumerate(MODELS.items()):
        sizes = [1., 2., 3.] if model == "spring_ring" else [20., 40., 60.]
        modeldata = rowdata[rowdata.model == model]
        for ci, ax in enumerate(axes[ri]):
            ax.set_xscale("log"); ax.set_xlim(7e-5, 1.45)
            ax.set_xticks([1e-4, 1e-3, 1e-2, .1, 1.])
            ax.set_xticklabels([r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", "1"])
            ax.grid(axis="y", color="#e6ebed", lw=.7); ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(labelsize=9)
            ax.set_xlabel("Input-strength factor f", fontsize=10)
        for si, (size, color) in enumerate(zip(sizes, COLORS)):
            d = modeldata[(modeldata["size"] == size) & (modeldata.kind == "noise")].set_index("fraction").reindex(GRID)
            control = modeldata[(modeldata["size"] == size) & (modeldata.kind == "mean")].set_index("fraction").reindex(GRID)
            specifications = [
                [(d, "multiregion_fraction", "-"), (d, "qualified_multiregion_fraction", "--")],
                [(d, "mean_R", "-"), (control, "mean_R", "--")],
                [(d, "mean_A", "-"), (d, "mean_Q", "--")],
            ]
            for ci, specs in enumerate(specifications):
                ax = axes[ri, ci]
                for data, metric, ls in specs:
                    y = data[metric].to_numpy(float)
                    ax.plot(GRID, y, color=color, ls=ls, lw=1.7 if ls == "-" else 1.35, alpha=1 if ls == "-" else .8)
                    if data is d:
                        lo = data[metric + "_seed_min"].to_numpy(float)
                        hi = data[metric + "_seed_max"].to_numpy(float)
                        ax.vlines(GRID, lo, hi, color=color, lw=.8, alpha=.45)
                        for f, value, n, target in zip(GRID, y, data.available_cases, data.target_cases):
                            marker(ax, f, value, color, n, target, size=4.3 if ls == "-" else 3.)
            missing = d.available_cases.fillna(0).to_numpy() == 0
            axes[ri, 0].plot(GRID[missing], np.full(np.sum(missing), -.115 - .065 * si),
                            transform=axes[ri, 0].get_xaxis_transform(), linestyle="none", marker="o", markersize=4,
                            markerfacecolor="white", markeredgecolor=color, clip_on=False)
        axes[ri, 0].set_ylim(-.025, 1.025)
        axes[ri, 0].set_ylabel("Fraction of sampled times", fontsize=10)
        axes[ri, 0].set_title(f"{chr(65+ri)}1  {name}: multiple regions", loc="left", fontsize=11)
        axes[ri, 1].set_title(f"{chr(65+ri)}2  First-moment magnitude", loc="left", fontsize=11)
        axes[ri, 1].set_ylabel(r"Mean $R=|P|$", fontsize=10)
        axes[ri, 2].set_title(f"{chr(65+ri)}3  Active amount and coherence", loc="left", fontsize=11)
        axes[ri, 2].set_ylabel("Mean normalized quantity", fontsize=10)
        for ci in (1, 2):
            ax = axes[ri, ci]; ax.axhline(1, color="#b4bdc3", lw=.7, zorder=0)
            ax.set_ylim(bottom=0)
        labels = [("C₀" if model == "spring_ring" else "C") + f"={size:g}" for size in sizes]
        axes[ri, 1].legend([Line2D([], [], color=c, lw=2) for c in COLORS], labels, loc="best", frameon=False, fontsize=8, ncol=3, columnspacing=.9)
    fig.text(.075, .077, "Regions: connected arcs above current minimum + half current contrast. Dashed occupancy also requires contrast ≥10% of own basal contrast.", fontsize=9.7)
    fig.text(.075, .059, r"$A=M_0/M_0(0)$; $Q=q_1/q_1(0)$; $R=AQ$ at each frame. Mean $R$ is the mean of the paired product, not the product of separate means.", fontsize=9.7)
    fig.text(.075, .041, "Filled circles: complete n=2; squares: n=12. Open data markers: partial n=1; open markers below occupancy panels: no noisy path available.", fontsize=9.7)
    fig.text(.075, .023, "Full windows: chemical 0–32000 (cadence 5); Spring 0–3200 (cadence 0.5). Connecting segments guide the eye. Native units are not calibrated across models.", fontsize=9.3)
    return fig


def profile_sheet(run, out):
    manifest = json.loads((run / "manifest.json").read_text())
    selected = []
    fields = []
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 6.9))
    fig.subplots_adjust(left=.075, right=.975, bottom=.25, top=.76, wspace=.3)
    fig.suptitle("What multiple threshold regions look like in native profiles", x=.075, ha="left", y=.975, fontsize=18)
    fig.text(.075, .926, "C=20, f=1; first prespecified pilot path per model. Snapshot conditioning is diagnostic, not representative sampling.", fontsize=11)
    snapshot_time = None
    for i, model in enumerate(("wave_pinning", "goryachev")):
        job, = [j for j in manifest["jobs"] if j["model"] == model and j["size"] == 20 and j["fraction"] == 1 and j["kind"] == "noise" and j["replicate"] == 0 and j["variant"] == "standard"]
        control, = [j for j in manifest["jobs"] if j["model"] == model and j["size"] == 20 and j["fraction"] == 1 and j["kind"] == "mean" and j["variant"] == "standard"]
        p = run / "workers" / f"{job['slot']:04d}"
        cp = run / "workers" / f"{control['slot']:04d}"
        with np.load(p / "primary.npz") as z:
            t, u, x = z["time"], z["primary"], z["x"]
        with np.load(p / "trace.npz") as z:
            scalar_t = z["time"]; ix = np.searchsorted(scalar_t, t)
            np.testing.assert_array_equal(scalar_t[ix], t)
            R = abs(z["P_real"] + 1j*z["P_imag"])[ix]
            A = (z["M0"]/z["M0"][0])[ix]; Q = (z["q1"]/z["q1"][0])[ix]
            regions = z["current_half_contrast_regions"][ix]
        if snapshot_time is None:
            eligible = (regions > 1) & (np.ptp(u, axis=1) >= .1*np.ptp(u[0])) & (R >= .75)
            if not eligible.any():
                raise RuntimeError("Prespecified WP profile qualification has no sampled frame")
            k = np.flatnonzero(eligible)[np.argmin(abs(t[eligible] - t[-1]/2))]
            snapshot_time = float(t[k])
        else:
            k, = np.flatnonzero(t == snapshot_time)
        with np.load(cp / "primary.npz") as z:
            ck, = np.flatnonzero(z["time"] == snapshot_time)
            mean = z["primary"][ck]
        threshold = float(u[k].min() + .5*np.ptp(u[k]))
        ax = axes[i]
        xx = np.r_[x[-1]-20, x, x[0]+20]
        for field, color, label, ls in ((u[0], "#879198", "Exact basal", "-"), (u[k], "#bc6332", "Noisy snapshot", "-"), (mean, "#263b49", "Matched mean input", "--")):
            ax.plot(xx, np.r_[field[-1], field, field[0]], color=color, lw=1.6, label=label, ls=ls)
        ax.axhline(threshold, color="#bc6332", lw=.9, ls=":", label="Current half-contrast")
        ax.set(xlim=(0, 20), ylim=(0, None), xlabel="Native physical arc coordinate", ylabel="Native active concentration")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(f"{MODELS[model]} · native t={snapshot_time:g}", loc="left", fontsize=12, y=1.115)
        ax.text(0, 1.035, f"{regions[k]} threshold regions; R={R[k]:.3f}, A={A[k]:.3f}, Q={Q[k]:.3f}", transform=ax.transAxes, fontsize=10)
        ax.legend(frameon=False, fontsize=9, loc="upper right")
        if model == "goryachev":
            inset = ax.inset_axes([.44, .29, .52, .34])
            inset.plot(x, u[k], color="#bc6332", lw=1.1)
            inset.plot(x, mean, color="#263b49", lw=1., ls="--")
            inset.axhline(threshold, color="#bc6332", lw=.7, ls=":")
            inset.set_xlim(0,20);inset.tick_params(labelsize=7)
            inset.set_title("Native-scale enlargement", fontsize=8)
        selected.append(dict(model=model,size=20,fraction=1,seed=job["seed"],slot=job["slot"],snapshot_time=snapshot_time,
                             selection="WP: qualifying profile nearest midwindow; Gory: same native time on first prespecified path",
                             regions=int(regions[k]),R=float(R[k]),A=float(A[k]),Q=float(Q[k]),threshold=threshold,
                             trace=str(p/"trace.npz"),trace_sha256=sha(p/"trace.npz"),profile=str(p/"primary.npz"),profile_sha256=sha(p/"primary.npz"),
                             mean_profile=str(cp/"primary.npz"),mean_profile_sha256=sha(cp/"primary.npz"),basal=str(run/job["input"]),basal_sha256=job["input_sha256"]))
        fields.extend(dict(model=model,x=float(a),initial=float(b),noisy=float(c),mean_control=float(d),time=snapshot_time) for a,b,c,d in zip(x,u[0],u[k],mean))
    fig.text(.075, .163, "WP selection: among its saved frames with >1 current half-contrast region, contrast ≥10% of basal, and R≥0.75, choose nearest the window midpoint.", fontsize=10)
    fig.text(.075, .137, "Gory uses the same native time on its own independent path. The two models have separate native concentration scales; neither profile is peak-normalized.", fontsize=10)
    fig.text(.075, .111, "The matched mean control receives the declared ensemble mean at every location. A low first moment or many threshold regions alone does not identify a mechanism.", fontsize=10)
    fig.text(.075, .085, "Threshold arcs are morphology descriptors: they do not establish equal-strength peaks, functional competition, site lineage, or an event class.", fontsize=10)
    pd.DataFrame(selected).to_csv(out / "supporting/data/profile_selection.csv", index=False)
    pd.DataFrame(fields).to_csv(out / "supporting/data/native_profile_cells.csv", index=False)
    return fig, selected


def export(fig, out, stem):
    fig.savefig(out / f"{stem}.png", dpi=170)
    for ext in ("pdf", "svg"):
        fig.savefig(out / "supporting/vectors" / f"{stem}.{ext}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, default=RUN); p.add_argument("--out", type=Path, default=OUT)
    a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    for sub in ("data", "metadata", "vectors"):
        (a.out / "supporting" / sub).mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans", "svg.fonttype":"none", "pdf.fonttype":42})
    df, status, csvsha = snapshot(a.run, a.out)
    fig = main_sheet(df, status); export(fig, a.out, "broad_strength_exploratory")
    profiles, selection = profile_sheet(a.run, a.out); export(profiles, a.out, "broad_strength_native_profiles")
    with PdfPages(a.out / "broad_strength_review.pdf") as pdf:
        pdf.savefig(fig); pdf.savefig(profiles)
    plt.close("all")
    outputs = {str(q.relative_to(a.out)):sha(q) for q in sorted(a.out.rglob("*")) if q.is_file() and q.name != "RENDER_RECEIPT.json"}
    savejson(a.out / "supporting/metadata/RENDER_RECEIPT.json", dict(generated_at=now(),native_observed_at=status["data_observed_at"],
              source=str(Path(__file__).resolve()),source_sha256=sha(__file__),manifest_sha256=sha(a.run/"manifest.json"),
              collector_snapshot_sha256=csvsha,new_completed_outputs=status["new_completed_outputs"],scientific_cases=status["scientific_cases_available"],
              native_producer_last_observed=status["producer"],main_figure_changed=False,outputs=outputs))
    print(json.dumps(dict(out=str(a.out),generated_at=now(),completed=status["new_completed_outputs"],profile_times=[r["snapshot_time"] for r in selection]),indent=2))


if __name__ == "__main__":
    main()
