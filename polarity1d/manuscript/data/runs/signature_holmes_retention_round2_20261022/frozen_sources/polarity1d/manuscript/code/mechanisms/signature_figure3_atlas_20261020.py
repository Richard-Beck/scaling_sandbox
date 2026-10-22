"""Complete retained native histories for every Figure 3 A1/A2/B/control case."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shlex
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[4]
MANUSCRIPT = REPO/"polarity1d/manuscript"
RUNS = MANUSCRIPT/"data/runs"
OUT = MANUSCRIPT/"figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance"
DEFAULT = RUNS/"signature_figure3_native_atlas_20261020"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def freeze(run, out):
    if run.exists():
        raise FileExistsError("Fresh immutable native-atlas run required")
    for folder in ("logs", "workers", "results"):
        (run/folder).mkdir(parents=True)
    shutil.copy2(__file__, run/"frozen_renderer.py")
    cases = []
    definitions = [("A1", "signature_figure3_campaign_20261020", "primary"),
        ("A2", "signature_figure3_regional_A2_20261020", None),
        ("B", "signature_figure3_regional_B_20261020", None),
        ("A1mesh", "signature_figure3_A1_controls_20261020", "controls"),
        ("Bmesh", "signature_figure3_regional_B_controls_20261020", None)]
    for stage, dirname, substage in definitions:
        source = RUNS/dirname
        manifest = read(source/"manifest.json")
        jobs = manifest["stages"][substage] if substage else manifest["jobs"]
        audit = source/f"results/{substage}_independent_native_audit.json" if substage else source/"results/all_independent_native_audit.json"
        evidence = read(audit)
        assert evidence["complete"] and evidence.get("verified", evidence.get("independently_checked")) == len(jobs)
        for job in jobs:
            statuspath = source/"workers"/substage/f'{job["slot"]:04d}'/"status.json" if substage else source/"workers"/f'{job["slot"]:04d}'/"status.json"
            status = read(statuspath)
            assert status["status"] == "complete"
            placement = status["placement"]
            case = dict(case_id=f'{stage}_{job["slot"]:04d}', stage=stage, source_run=str(source),
                source_manifest_sha256=sha(source/"manifest.json"), source_audit_sha256=sha(audit),
                source_status=str(statuspath), source_status_sha256=sha(statuspath),
                model=placement["model"], size_parameter=placement["size_parameter"],
                initial_gap=placement["edge_gap_10"], classification=placement["classification_10"],
                target_ratio=job.get("target_ratio"), multiplier=job.get("multiplier"), control=job.get("control"))
            if substage:
                proposals = (RUNS/"signature_figure3_region_proposal_v2_20261020/results" if stage == "A1" else source/"results/subcell_regions")
                proposalpath = proposals/f'{substage}_{job["slot"]:04d}_placement_{placement["placement_id"]:03d}_region_proposal.json'
                proposal = read(proposalpath)
                case.update(region_trace=proposal["traces"], region_trace_sha256=proposal["traces_sha256"])
            cases.append(case)
    chunks = [list(range(i, min(i+20, len(cases)))) for i in range(0, len(cases), 20)]
    plan = dict(cases=cases, chunks=chunks, output=str(out), renderer_sha256=sha(run/"frozen_renderer.py"),
        purpose="All retained native active-field histories/profiles, complete-region amounts and raw activity-vector coherence/direction; full native-species archives linked", native_integrations=0)
    write(run/"manifest.json", plan)
    command = "exec "+shlex.join(["env", "PYTHONDONTWRITEBYTECODE=1", "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1",
        "MKL_NUM_THREADS=1", "/home/4473331/.conda/envs/scaling_sandbox/bin/python", str(run/"frozen_renderer.py"), "worker", "--run", str(run)])+' --slot "$SLURM_ARRAY_TASK_ID"'
    write(run/"slurm_payload.json", dict(jobs=len(chunks), cases=len(cases), submitted=False, sbatch_arguments=["sbatch", "--parsable",
        "--qos=small", "--ntasks=1", "--cpus-per-task=1", "--mem=4G", "--time=01:00:00",
        f"--array=0-{len(chunks)-1}", "--job-name=figure3-native-atlas", "--output="+str(run/"logs/%A_%a.out"),
        "--error="+str(run/"logs/%A_%a.err"), "--wrap", command]))
    return dict(cases=len(cases), jobs=len(chunks))


def render(case, output):
    source = Path(case["source_run"])
    statuspath = Path(case["source_status"])
    assert sha(statuspath) == case["source_status_sha256"]
    record = read(statuspath)
    placement = record["placement"]
    cells = placement["cells"]
    tt, yy = [], []
    for phase in record["phases"]:
        path = source/phase["path"]
        assert sha(path) == phase["sha256"]
        with np.load(path) as z:
            t, y = z["time"].copy(), z["states"][:, :cells].copy()
        use = slice(None) if not tt else slice(1, None)
        tt.extend((phase["absolute_start"]+t[use]).tolist()); yy.extend(y[use])
    t, y = np.asarray(tt), np.asarray(yy)
    old = placement["initial_10_percent_contour"]["origin"]/placement["basal_circumference"]
    q = ((np.arange(cells)+.5)/cells-old+.25) % 1
    order = np.argsort(q)
    q, y = q[order], y[:, order]
    trace = Path(case["region_trace"]) if "region_trace" in case else source/record["readouts"]
    assert sha(trace) == case.get("region_trace_sha256", record.get("readouts_sha256"))
    with np.load(trace) as z:
        rt = z["time"].copy()
        incumbent, challenger = z["incumbent_excess"].copy(), z["challenger_excess"].copy()
    vpath = source/record["readouts"]
    assert sha(vpath) == record["readouts_sha256"]
    with np.load(vpath) as z:
        vt = z["time"].copy()
        vector = z["vector_real"]+1j*z["vector_imag"]
    plt.rcParams.update({"font.size":10, "pdf.fonttype":42})
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 7.7))
    ax = axes[0, 0]
    image = ax.pcolormesh(q, t, y, shading="nearest", cmap="viridis", rasterized=True)
    ax.set_yscale("symlog", linthresh=1. if case["model"] == "spring" else 10.)
    ax.invert_yaxis()
    ax.set_xlabel("Native material arc fraction (rotated for display)")
    ax.set_ylabel("Native time, full retained horizon")
    figure.colorbar(image, ax=ax, label="Native active signal")
    withdrawal = record.get("withdrawal_time")
    if withdrawal is not None:
        ax.axhline(withdrawal, color="white", ls="--", lw=1.)
    ax = axes[0, 1]
    wanted = [0., t[-1]] if withdrawal is None else [0., withdrawal, withdrawal+(t[-1]-withdrawal)/2, t[-1]]
    used = set()
    for point in wanted:
        index = int(np.argmin(abs(t-point)))
        if index in used: continue
        used.add(index)
        ax.plot(q, y[index], label=f"t={t[index]:.5g}", lw=1.3)
    ax.set_xlabel("Native material arc fraction (rotated)")
    ax.set_ylabel("Native active Rac" if case["model"] == "holmes_model3" else "Native active concentration")
    ax.legend(fontsize=8, ncol=2)
    ax = axes[1, 0]
    ax.plot(rt, incumbent, color="#21689a", label="Incumbent-associated complete region")
    ax.plot(rt, challenger, color="#e58025", label="Challenger-associated complete region")
    ax.set_xscale("symlog", linthresh=1. if case["model"] == "spring" else 10.)
    ax.set_xlabel("Native time"); ax.set_ylabel("Active excess (concentration × length)")
    ax.legend(fontsize=8)
    if withdrawal is not None: ax.axvline(withdrawal, color="#555555", ls="--", lw=.8)
    ax = axes[1, 1]
    ax.plot(vt, abs(vector), color="#1a7466", label="Coherence R")
    ax.axhline(.05, color="#aaaaaa", ls=":", lw=.8)
    ax.set_ylabel("R = |M₁| / M₀", color="#1a7466")
    ax.set_xscale("symlog", linthresh=1. if case["model"] == "spring" else 10.)
    ax.set_xlabel("Native time")
    angle = ax.twinx()
    eligible = abs(vector) >= .05
    angle.scatter(vt[eligible], np.angle(vector[eligible]), s=2, color="#74418a", alpha=.6)
    angle.set_ylabel("Wrapped heading (rad), R≥0.05", color="#74418a")
    angle.set_ylim(-np.pi, np.pi)
    for axis in axes.flat:
        axis.spines[["top", "right"]].set_visible(False)
    title = f'{case["stage"]} · {case["model"]} · size={case["size_parameter"]:g} · initial g={case["initial_gap"]:.4g} · {case["classification"]}'
    if case.get("target_ratio") is not None: title += f' · withdrawal ratio={case["target_ratio"]:g}'
    if case.get("multiplier") is not None: title += f' · pulse={case["multiplier"]:g}× own handover'
    if case.get("control"): title += ' · '+case["control"]
    figure.suptitle(title, fontsize=12)
    figure.text(.5, .013, "Every retained native active-field sample is shown. Dashed line: cue withdrawal. Region identity is spatial, not molecular ancestry.\n"
        "NaN amount intervals remain gaps. Wrapped headings are points; no unwrapping through cancellation. Full native-species archives are linked in the reader.",
        ha="center", fontsize=8)
    figure.tight_layout(rect=(0, .06, 1, .94))
    target = output/"supporting/viewers/native_histories"/(case["case_id"]+".png")
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=130, facecolor="white")
    plt.close(figure)
    return dict(case_id=case["case_id"], png=str(target), png_sha256=sha(target), native_frames=len(t), native_cells=cells,
                source_status_sha256=case["source_status_sha256"], full_horizon=float(t[-1]))


def worker(run, slot):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Complete native atlas rendering requires SLURM")
    plan = read(run/"manifest.json")
    assert sha(Path(__file__)) == plan["renderer_sha256"]
    output = Path(plan["output"])
    records = [render(plan["cases"][i], output) for i in plan["chunks"][slot]]
    write(run/"workers"/f"{slot:04d}.json", dict(complete=True, records=records))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "worker"))
    parser.add_argument("--run", type=Path, default=DEFAULT)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--slot", type=int)
    args = parser.parse_args()
    result = freeze(args.run.resolve(), args.output.resolve()) if args.action == "freeze" else worker(args.run.resolve(), args.slot)
    if result is not None: print(json.dumps(result))


if __name__ == "__main__":
    main()
