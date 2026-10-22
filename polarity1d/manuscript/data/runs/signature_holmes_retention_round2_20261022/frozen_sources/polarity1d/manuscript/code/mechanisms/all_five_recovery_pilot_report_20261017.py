"""Render the completed bounded recovery pilot without further integration."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import shutil
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/"data/runs/all_five_incumbent_recovery_pilot_20261017"
OUT = ROOT/"figures/length_mechanisms_20261008/dev/all_five_recovery_pilot_20261017"
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", debelly="deBelly", holmes_model3="Holmes 3", spring="Spring")
COLORS = dict(short="#0072b2", middle="#7b4990", long="#d55e00")
FILES = []


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def csvread(path): return list(csv.DictReader(Path(path).open()))


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        path = OUT/(stem+".png") if ext == "png" else OUT/"supporting/vectors"/(stem+"."+ext)
        fig.savefig(path, dpi=220 if ext == "png" else None)
        FILES.append(path)
    plt.close(fig)


def render(clocks, curves):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 5, figsize=(13.8, 6.7), sharey="row")
    fig.subplots_adjust(left=.065, right=.985, top=.805, bottom=.18, wspace=.27, hspace=.75)
    for col, clock in enumerate(clocks["models"]):
        model = clock["model"]; tau = clock["selected_horizon_native"]
        top, bottom = axes[:, col]
        assert tau is not None
        selected = [r for r in curves if r["model"] == model and float(r["epsilon"]) == .005]
        native_sizes = []
        for label, color in COLORS.items():
            rows = [r for r in selected if r["label"] == label]
            times = np.array([float(r["time_native"]) for r in rows])
            F = np.array([float(r["restored_fraction"]) for r in rows])
            E = np.array([float(r["empirical_numerical_sensitivity"]) for r in rows])
            keep = times <= tau*(1+1e-12)
            top.plot(times[keep]/tau, F[keep], color=color, lw=1.5)
            top.fill_between(times[keep]/tau, F[keep]-E[keep], F[keep]+E[keep], color=color, alpha=.2)
            endpoint = int(np.argmin(abs(times-tau)))
            top.scatter([1], [F[endpoint]], color=color, s=17, zorder=4)
            positive = times > 0
            bottom.plot(times[positive], 1-F[positive], color=color, lw=1.5)
            bottom.scatter([tau], [1-F[endpoint]], color=color, s=17, zorder=4)
            native_sizes.append(float(rows[0]["size_native"]))
        top.set(xlim=(0, 1.05), ylim=(-.005, .21), xticks=[0, .5, 1.], yticks=[0, .1, .2], xlabel="Time / own selected window")
        top.axhline(.2, color="#999999", lw=.7, ls="--")
        top.set_title(f'{"ABCDE"[col]}  {NAMES[model]}', loc="left", pad=44, fontsize=10)
        top.text(.5, 1.055, f'Sizes {" / ".join(f"{size:g}" for size in native_sizes)}\nτ = {tau:.6g} native time',
                 transform=top.transAxes, ha="center", va="bottom", fontsize=8.2)
        bottom.set(xscale="log", ylim=(0, 1.04), yticks=[0, .2, .4, .6, .8, 1.], xlabel="Time after kick (native)")
        bottom.axhline(.8, color="#999999", lw=.7, ls="--")
        bottom.axvline(tau, color="#555555", lw=.75, ls=":")
        bottom.set_title(f'{"FGHIJ"[col]}  Complete short history', loc="left", fontsize=9.5, pad=7)
        for ax in (top, bottom): ax.grid(axis="y", color="#eeeeee", lw=.6)
    axes[0, 0].set_ylabel("Induced core deficit restored, F")
    axes[1, 0].set_ylabel("Induced core deficit remaining, D(t) / D(0)")
    fig.suptitle("A common active-site deactivation gives a resolved early recovery readout in all five linear models", x=.52, y=.985, fontsize=12)
    fig.legend(handles=[Line2D([], [], color=color, lw=1.7, label=label.capitalize()+" eligible size") for label, color in COLORS.items()],
               loc="upper center", bbox_to_anchor=(.52, .945), ncol=3, frameon=False, fontsize=9)
    fig.text(.525, .045, "0.5% basal excess × smooth incumbent taper removed; explicit inactive pools receive it where present.\n"
             "deBelly changes Rac only, retaining all mechanical/controller state. A–E use different declared native windows; clocks are not calibrated across models.\n"
             "Positive conditional recovery concerns the paired deficit, not generic stability or resistance to a challenger. Native grid retained; Goryachev L39 core has four centers.",
             ha="center", va="center", fontsize=8.2, linespacing=1.5)
    save(fig, "all_five_recovery_pilot_curves")


def main():
    for category in ("vectors", "data", "metadata", "viewers"):
        (OUT/"supporting"/category).mkdir(parents=True, exist_ok=True)
    clocks = read(RUN/"results/clock_decision.json")
    verify = read(RUN/"results/verification.json")
    assert read(RUN/"results/status.json")["complete"]
    curves = csvread(RUN/"results/short_curves.csv")
    gates = csvread(RUN/"results/candidate_gates.csv")
    render(clocks, curves)
    selected = []
    for clock in clocks["models"]:
        middle = clock["primary_results"][1]["primary_K"]
        for row in clock["primary_results"]:
            detail = next(r for r in curves if r["model"] == row["model"] and float(r["size_native"]) == row["size_native"]
                          and float(r["epsilon"]) == .005 and abs(float(r["time_native"])-row["horizon_native"]) < 1e-14)
            selected.append({**row, "K_relative_to_middle": row["primary_K"]/middle,
                             **{key: float(detail[key]) for key in ("initial_core_damage", "core_deficit", "whole_active_deficit",
                                  "outside_core_deficit", "sham_core_change", "probe_core_change")}})
    table = OUT/"supporting/data/all_five_recovery_pilot_selected.csv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0])); writer.writeheader(); writer.writerows(selected)
    FILES.append(table)
    for filename in ("clock_decision.json", "verification.json", "completion_metadata_recovery.json"):
        path = OUT/"supporting/metadata"/("all_five_recovery_pilot_"+filename)
        shutil.copy2(RUN/"results"/filename, path); FILES.append(path)
    for filename in ("short_curves.csv", "candidate_gates.csv"):
        path = OUT/"supporting/data"/("all_five_recovery_pilot_"+filename)
        shutil.copy2(RUN/"results"/filename, path); FILES.append(path)
    model_lines = []
    for clock in clocks["models"]:
        rows = clock["primary_results"]; ks = [r["primary_K"] for r in rows]
        model_lines.append(f'| {NAMES[clock["model"]]} | '+" / ".join(f'{r["size_native"]:g}' for r in rows)+f' | {clock["selected_horizon_native"]:.8g} | '+" / ".join(f'{value:.7g}' for value in ks)+f' | {ks[2]/ks[1]:.6g} |')
    failures = []
    for clock in clocks["models"]:
        group = [r for r in gates if r["model"] == clock["model"]]
        failures.append(f'| {NAMES[clock["model"]]} | '+" | ".join(str(sum(r[key] == "False" for r in group)) for key in ("error_pass", "linearity_pass", "saturation_pass"))+" |")
    report = f'''# All-five linear active-site recovery pilot

**Result.** The same localized activity-loss question gives a numerically resolved, amplitude-consistent early response in all five retained linear models. Goryachev shows the largest middle-to-long increase. The other models change more modestly over those pilot pairs; broadening alone does not imply the same functional response.

![All-five native recovery pilot](all_five_recovery_pilot_curves.png)

**Message conveyed.** A small incumbent-targeted activity reduction can be compared across models without pretending deBelly has an inactive Rac pool. The resulting deficit shrinks at every selected clock in all 15 native cases. This supplies a candidate all-five Figure 1 C readout. It is not an all-length result or a universal stability law.

**Why this figure.** The early traces display the actual quantity summarized by the coefficient on one shared restored-fraction scale. The complete short histories show how the deficit continues to shrink beyond the selected early window and why later clocks exceed the 20% range. Their logarithmic time axes do not establish accelerating recovery. Both panels retain short, middle and long cases. Error bands are shown but are smaller than line width. Different selected windows do not calibrate native rates across models; compare size effects within a model.

## Selected finite-window responses

The coefficient is

\\[
K_{{\\tau}}=\\frac{{D(0)-D(\\tau)}}{{\\tau D(0)}}.
\\tag{{1}}
\\]

Here D is the difference between separately evolved sham and probe active amounts inside the frozen original half-contrast core. Its native inverse-time units differ across models. Each model uses one window for all three sizes and all amplitudes. The final column compares long versus middle **at the same within-model clock**, not models against each other.

| Model | Short / middle / long nominal size | Selected τ, native time | Kτ: short / middle / long | Long / middle Kτ |
|---|---|---:|---|---:|
{chr(10).join(model_lines)}

For WP and Goryachev the middle/long pair is 19 → 39, matching Figure 1 A. Other pilot pairs are deBelly 27 → 39, Holmes 28 → 50 and Spring L₀ 1.6 → 2.6; they must not be relabelled as Figure 1 A's highlighted pairs. The full all-97 C campaign has **not** been launched. All-five initial cue gains are separately available from existing states in the [design evidence](../all_five_linear_function_20261017/DESIGN.md).

At the selected clocks, all primary changes exceed three times the empirical sensitivity estimate. Maximum prefix sensitivity is below 2.33 × 10⁻⁷, against the fixed 10⁻⁵ limit. All selected prefix saturation and amplitude-linearity checks pass. Actual probe core changes are positive in all 15 cases, while sham changes are zero to roundoff except the small deBelly L18 basal drift (4.17 × 10⁻⁷ active amount over 5 s). In later challenged-state tests this need not remain true: a positive conditional coefficient alone cannot establish net maintenance during competition.

## Intervention and interpretation

The [frozen design](../../../../data/runs/all_five_incumbent_recovery_pilot_20261017/frozen_sources/figures/length_mechanisms_20261008/dev/all_five_linear_function_20261017/DESIGN.md) defines a raised-cosine taper anchored to the original basal half-contrast core, multiplied by its excess-active profile. Epsilon = 0.005 attenuates local basal excess by at most 0.5%; sensitivities are 0.001 and 0.01. Epsilon is not an identical integrated fractional core damage across models. WP/G/H/Spring transfer the removed amount pointwise to their explicit corresponding inactive pool. deBelly removes Rac activity alone; Rho, MCA and all mechanical and endpoint-controller fields start unchanged. Native coupled dynamics then evolve with zero extra cue, without clamps, field freezes or later material return.

The perturbation is sampled at native nodes or finite-volume centers. An edge finite-volume cell can straddle the interpolated core edge, so the full-domain and outside-core differences are exported separately. Local recovery includes reaction and transport, not only molecular reactivation. Spring's core stays fixed in material coordinates, with sham and probe each integrated using its own evolving physical length. deBelly's native half-second stagger and 0.98 endpoint relaxation remain unchanged. Its dense reference tests arithmetic of this same discrete update; it does not establish temporal convergence.

## Complete clock and failure record

Continuous candidate clocks are exactly logspace(−6, 0, 25); deBelly candidates are 0.5, 1, 2, 5 and 10 s. All candidates are sampled, with a denser continuous history and every native deBelly half-step. Every sampled prefix through a candidate must satisfy empirical sensitivity E ≤ 10⁻⁵, pairwise 5% amplitude agreement with a 3E allowance, and |F| + E ≤ 0.2. No minimum effect is required. The largest passing clock is selected once per model.

Counts below are **case × candidate** exclusions over the complete prescribed list; overlapping failures would count in more than one column. All exclusions here arise from exceeding the early-response range, rather than numerical or linearity failures. They are not failures of native polarization or eventual repair.

| Model | Error exclusions | Linearity exclusions | Saturation exclusions |
|---|---:|---:|---:|
{chr(10).join(failures)}

The selected clocks are pilot decisions for later root review. A ring port must check the extended size range and challenged states before treating these clocks as portable. No per-size or separation-specific clock adjustment is allowed.

## Verification, provenance and limits

All 15 exact source checkpoint hashes, the 97-case selection table, nonnegative kicks, conservative compensation where defined and untouched other initial fields were checked before freezing. Zero-kick executions reproduce the same sham fields exactly. Both Holmes pools, WP/G/Spring active-plus-inactive material and deBelly bound-plus-unbound MCA were tracked. Maximum pool drift was {verify["maximum_pool_relative_drift"]:.3g}. All deBelly banded and independent dense solves passed the scaled linear residual criterion (maximum {verify["maximum_debelly_scaled_linear_residual"]:.3g}); all 11 field/controller differences are retained.

Native meshes are unchanged. Goryachev L39's interpolated core contains only four sample centers; its response is a result of that retained discretization, not a mesh-convergence result. No mesh campaign was added to this pilot.

SLURM array 21125088 completed the scientific integrations for all 15 cases. The three Spring workers saved complete scientific archives and then failed while serializing NumPy-int32 LSODA counters. A separate [hashed metadata-recovery audit](../../../../data/runs/all_five_incumbent_recovery_pilot_20261017/results/completion_metadata_recovery.json) reconstructed all ten trajectories per affected case, recomputed moving-domain amounts/native RHS rates and pool checks, verified unchanged archive hashes and independent zero-kick identities, and preserved the original failed statuses/error logs. No integration was repeated. Those original solver counters are unavailable; no values were invented. LSODA's stored t = 0 polynomial output differed from the exact supplied restart by at most one machine ulp (2.23 × 10⁻¹⁶), included in the roundoff allowance.

- [Selected 15-case table](supporting/data/all_five_recovery_pilot_selected.csv), [all short curves](supporting/data/all_five_recovery_pilot_short_curves.csv), and [every candidate gate](supporting/data/all_five_recovery_pilot_candidate_gates.csv).
- [Clock decision](supporting/metadata/all_five_recovery_pilot_clock_decision.json), [full verification](supporting/metadata/all_five_recovery_pilot_verification.json), and [metadata recovery](supporting/metadata/all_five_recovery_pilot_completion_metadata_recovery.json).
- [Frozen run manifest and sources](../../../../data/runs/all_five_incumbent_recovery_pilot_20261017/manifest.json), [original job launch](../../../../data/runs/all_five_incumbent_recovery_pilot_20261017/launch.json), and [native full-field archives](../../../../data/runs/all_five_incumbent_recovery_pilot_20261017/workers/).
- [PDF](supporting/vectors/all_five_recovery_pilot_curves.pdf) and [SVG](supporting/vectors/all_five_recovery_pilot_curves.svg).

**Stopping point:** the bounded pilot is complete. Root review and fresh actual-result critique precede any full grid or updated main-text coefficient claim. No migration model or deBelly ring extension was run.
'''
    path = OUT/"RESULTS.md"; path.write_text(report); FILES.append(path)
    sources = [Path(__file__), ROOT/"code/mechanisms/all_five_incumbent_recovery_20261017.py",
               ROOT/"code/mechanisms/all_five_recovery_metadata_repair_20261017.py",
               RUN/"manifest.json", RUN/"provenance.json", RUN/"results/clock_decision.json",
               RUN/"results/verification.json", RUN/"results/short_curves.csv", RUN/"results/candidate_gates.csv"]
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), artifact_root="../..",
                    source_sha256={str(p.relative_to(ROOT.parents[1])): sha(p) for p in sources},
                    sha256={str(p.relative_to(OUT)): sha(p) for p in FILES},
                    integration_scope="Completed bounded15-case pilot only; no new integration from report generation")
    (OUT/"supporting/metadata/all_five_recovery_pilot_report.manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(dict(report=str(OUT/"RESULTS.md"), cases=len(selected), models=len(clocks["models"]), exports=len(FILES)), indent=2))


if __name__ == "__main__": main()
