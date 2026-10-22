"""Linear checkpoint evidence report from completed, audited native campaigns."""
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
RUN = ROOT/"data/runs/vision_linear_20261017"
OUT = ROOT/"figures/length_mechanisms_20261008/dev/vision_linear_20261017"
GAIN = ROOT/"figures/length_mechanisms_20261008/dev/all_five_linear_function_20261017/supporting/data/native_gain_design_evidence.csv"
FOLDS = ROOT/"figures/length_mechanisms_20261008/dev/linear_then_ring_checkpoint_20261017/supporting/data/all_five_incumbent_size_folds.csv"
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
NAMES = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "deBelly", "Holmes 3", "Spring")))
COLORS = dict(zip(MODELS, ("#007d73", "#893f82", "#d55e00", "#0072b2", "#8c6d31")))
REFERENCES = dict(zip(MODELS, (19., 19., 19., 20., 1.)))
PAIRS = dict(zip(MODELS, ((19., 39.), (19., 39.), (19., 39.), (20., 40.), (1., 2.))))
FILES = []


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path): return json.loads(Path(path).read_text())
def rows(path): return list(csv.DictReader(Path(path).open()))
def flag(value): return value is True or value == "True"


def copy_file(source, relative):
    target = OUT/relative; target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target); FILES.append(target)
    return target


def csvwrite(path, data):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    FILES.append(path)


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        target = OUT/(stem+".png") if ext == "png" else OUT/"supporting/vectors"/(stem+"."+ext)
        fig.savefig(target, dpi=220 if ext == "png" else None, bbox_inches="tight")
        FILES.append(target)
    plt.close(fig)


def normalize(gains, recoveries):
    output = []
    for model in MODELS:
        reference_gain = next(r for r in gains if r["model"] == model and float(r["size_native"]) == REFERENCES[model])
        reference_repair = next(r for r in recoveries if r["model"] == model and float(r["size_native"]) == REFERENCES[model])
        assert flag(reference_repair["coefficient_available"])
        for row in [r for r in recoveries if r["model"] == model]:
            gain = next(r for r in gains if r["model"] == model and float(r["size_native"]) == float(row["size_native"]))
            output.append(dict(model=model, size_native=float(row["size_native"]), physical_length=float(row["physical_length"]),
                reference_size_native=REFERENCES[model], size_relative_to_reference=float(row["size_native"])/REFERENCES[model],
                cue_weighted_gain_native=float(gain["cue_weighted_gain_native"]),
                cue_gain_relative_to_reference=float(gain["cue_weighted_gain_native"])/float(reference_gain["cue_weighted_gain_native"]),
                horizon_native=float(row["horizon_native"]), restored_fraction=float(row["restored_fraction"]),
                conditional_K_per_native_time=float(row["conditional_K_per_native_time"]), empirical_E=float(row["empirical_E"]),
                repair_relative_to_reference=float(row["conditional_K_per_native_time"])/float(reference_repair["conditional_K_per_native_time"]),
                coefficient_available=flag(row["coefficient_available"]), model_clock_qualified=flag(row["model_clock_qualified"]),
                checkpoint_sha256=row["checkpoint_sha256"], recovery_archive_sha256=row["archive_sha256"]))
    assert len(output) == 97
    return output


def function_figure(folds, normalized):
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 4.15))
    fig.subplots_adjust(left=.135, right=.985, bottom=.24, top=.79, wspace=.55)
    ax = axes[0]
    for field, marker, label, delta in (("excess_active_fold", "o", "Excess active amount", -.17),
            ("width_fold", "s", "Half-contrast width", 0.), ("peak_concentration_fold", "^", "Peak concentration", .17)):
        ax.scatter([float(r[field]) for r in folds], np.arange(5)+delta, marker=marker, s=26,
                   facecolor="#333333" if marker == "o" else "white", edgecolor="#333333", lw=.8, label=label)
    ax.set(yticks=range(5), yticklabels=[NAMES[r["model"]] for r in folds], ylim=(4.6, -.6), xlim=(.8, 3.2),
           xlabel="Larger / smaller value", title="A  Established incumbent organization")
    ax.axvline(1, ls="--", color="#999999", lw=.7)
    ax.legend(loc="lower left", bbox_to_anchor=(-.12, 1.07), frameon=False, fontsize=7.7)
    for axis, key, title in ((axes[1], "cue_gain_relative_to_reference", "B  Initial fixed-support cue gain"),
                             (axes[2], "repair_relative_to_reference", "C  Common local-core recovery")):
        all_values = []
        for model in MODELS:
            group = sorted([r for r in normalized if r["model"] == model], key=lambda r:r["size_native"])
            valid = [r for r in group if key != "repair_relative_to_reference" or r["coefficient_available"]]
            x = [r["size_relative_to_reference"] for r in valid]; y = [r[key] for r in valid]
            all_values += y
            axis.plot(x, y, color=COLORS[model], lw=1.3, marker="o", ms=2.8, label=NAMES[model])
            for size in PAIRS[model]:
                item = next(r for r in valid if r["size_native"] == size)
                axis.scatter([item["size_relative_to_reference"]], [item[key]], facecolors="white", edgecolors=COLORS[model], s=25, lw=.8, zorder=4)
        axis.set_yscale("log", base=2)
        low = 2.**np.floor(np.log2(min(all_values))); high = 2.**np.ceil(np.log2(max(all_values)))
        ticks = 2.**np.arange(np.floor(np.log2(low)), np.ceil(np.log2(high))+1)
        axis.set(ylim=(low*.88, high*1.12), yticks=ticks,
                 yticklabels=[f"{v:g}" if v >= 1 else f"1/{1/v:g}" for v in ticks], xlim=(0., 2.75), xticks=[0., 1., 2.],
                 xlabel="Native size / own reference", ylabel="Relative to own reference (log₂ scale)", title=title)
        axis.axhline(1, ls="--", lw=.7, color="#999999"); axis.grid(axis="y", color="#eeeeee", lw=.6)
    fig.legend(handles=[Line2D([], [], color=COLORS[m], lw=1.4, label=NAMES[m]) for m in MODELS],
               loc="upper center", bbox_to_anchor=(.67, .99), ncol=3, frameon=False, fontsize=8)
    fig.text(.55, .045, "All five established linear models; 97 exact qualified basals. References: L19, Holmes L20, Spring L₀1.\n"
             "A size pairs: 19→39; Holmes 20→40; Spring L₀1→2. B and C have separate shared logarithmic axes.\n"
             "C uses one fixed native window per model; rates are not calibrated across models. Positive conditional repair is not generic stability.",
             ha="center", va="center", fontsize=8.1, linespacing=1.5)
    save(fig, "vision_linear_all_five_function")


def finite_detail(finite, samples):
    primary = [r for r in finite if r["status"] == "complete" and flag(r["primary"])]
    fig, axes = plt.subplots(2, 5, figsize=(13.8, 6.1), sharey="row")
    fig.subplots_adjust(left=.065, right=.99, top=.80, bottom=.13, wspace=.25, hspace=.55)
    for col, model in enumerate(MODELS):
        group = sorted([r for r in primary if r["model"] == model], key=lambda r:float(r["size_native"]))
        x = [float(r["size_native"]) for r in group]
        top, bottom = axes[:, col]
        for key, color, label in (("withdrawal_native_orientation", "#b24a33", "At cue withdrawal"),
                                  ("post13000_native_orientation", "#684c8b", "13,000 after withdrawal")):
            top.plot(x, [float(r[key]) for r in group], color=color, marker="o", ms=2.5, lw=1.2, label=label)
        for site, color, label in (("old", "#325c9b", "Original core"), ("challenger", "#b74b32", "Mirrored challenger core")):
            bottom.plot(x, [float(r[f"post13000_{site}_excess_ratio"]) for r in group], color=color,
                        marker="o", ms=2.5, lw=1.2, label=label)
        top.set_title(f'{NAMES[model]}\nn = {len(group)}; cue {50 if model == "spring" else 1000} native units', fontsize=9.5, pad=10)
        top.set(ylim=(-1.02, 1.02), yticks=[-1, 0, 1]); top.axhline(0, color="#aaaaaa", lw=.6)
        bottom.set(ylim=(-.03, 1.07), yticks=[0, .5, 1.], xlabel="Resting L₀" if model == "spring" else "Interval length")
        for axis in (top, bottom): axis.grid(axis="y", color="#eeeeee", lw=.6)
    axes[0, 0].set_ylabel("Native activity polarity\n(signed first moment)")
    axes[1, 0].set_ylabel("Core excess after 13,000 /\ninitial incumbent excess")
    axes[0, 0].legend(frameon=False, fontsize=7, loc="lower left")
    axes[1, 0].legend(frameon=False, fontsize=7, loc="center left")
    fig.suptitle("Finite competition distinguishes cue-on activity from retained local organization", y=.985, fontsize=12)
    fig.text(.53, .91, "Fixed endpoint support 3 (Spring 0.2), unchanged native peaks, all 97 qualified sizes; 13 preparation-ineligible sizes retained in data.\n"
             "Shared scales within each row. First moments do not uniquely identify a site, protrusion or cell movement.", ha="center", fontsize=8.5, linespacing=1.5)
    save(fig, "vision_linear_finite_retained_detail")


def main():
    for sub in ("vectors", "data", "metadata", "viewers"): (OUT/"supporting"/sub).mkdir(parents=True, exist_ok=True)
    recovery = rows(RUN/"recovery/results/recovery_summary.csv")
    gain = rows(GAIN); finite = rows(RUN/"finite/results/metrics.csv"); samples = rows(RUN/"finite/results/samples.csv")
    assert len(recovery) == 97 and len(gain) == 97
    assert read(RUN/"recovery/results/status.json")["complete"] and read(RUN/"finite/results/status.json")["complete"]
    normalized = normalize(gain, recovery)
    function_figure(rows(FOLDS), normalized); finite_detail(finite, samples)
    csvwrite(OUT/"supporting/data/vision_linear_normalized_function.csv", normalized)
    copy_file(RUN/"recovery/results/recovery_summary.csv", "supporting/data/vision_linear_recovery_summary.csv")
    copy_file(RUN/"finite/results/metrics.csv", "supporting/data/vision_linear_finite_metrics.csv")
    copy_file(RUN/"finite/results/samples.csv", "supporting/data/vision_linear_finite_samples.csv")
    copy_file(GAIN, "supporting/data/vision_linear_gain_all97.csv")
    copy_file(FOLDS, "supporting/data/vision_linear_incumbent_size_pairs.csv")
    function_lines = []
    finite_lines = []
    for model in MODELS:
        lo, hi = [next(r for r in normalized if r["model"] == model and r["size_native"] == size) for size in PAIRS[model]]
        function_lines.append(f'| {NAMES[model]} | {lo["size_native"]:g} → {hi["size_native"]:g} | {hi["cue_gain_relative_to_reference"]:.6g} | {hi["repair_relative_to_reference"]:.6g} | {lo["horizon_native"]:.8g} |')
        group = [r for r in finite if r["model"] == model and r["status"] == "complete" and flag(r["primary"])]
        finite_lines.append(f'| {NAMES[model]} | {len(group)} | '+" / ".join(str(sum(flag(r[f"post10000_{key}"]) for r in group)) for key in ("orientation_success", "full_polarity_success", "field_tail_pass"))+" | "+" / ".join(str(sum(flag(r[f"post13000_{key}"]) for r in group)) for key in ("orientation_success", "full_polarity_success", "field_tail_pass"))+" |")
    rverify = read(RUN/"recovery/results/verification.json")
    all_clocks = rverify["all_fixed_clocks_qualified"]
    assert all_clocks, "Report fixed-clock failures before rendering an unqualified common trend"
    text = f'''# All-five linear checkpoint evidence: organization, response and finite competition

**Message.** Enlargement increases incumbent material in every established model, while its consequences differ. The complete gain and common local-recovery series show strong size dependence in Goryachev and smaller or differently structured responses in the other models. A fixed local finite cue then exposes replacement, incumbent return and persistent two-region activity. These finite outcomes are displayed alongside the transient assays rather than inferred from them.

## Organization and local function across every qualified size

![All-five established linear organization and complete functional series](vision_linear_all_five_function.png)

**Why this figure.** The accepted amount/width/peak comparison anchors the functional question. Every qualified size enters both response curves. Separate shared log₂ axes retain the full short-cell range and keep nearly unchanged responses near one; no model receives its own stretched scale. Open markers identify the exact organization size pairs. All values are linear; ring coefficients do not enter this figure.

| Model | Exact size pair | Initial gain fold | Common local-repair fold | Fixed recovery τ, native units |
|---|---|---:|---:|---:|
{chr(10).join(function_lines)}

Initial gain is the extra primary-activity production per unit native activation input, averaged over the fixed endpoint cue. It is inactive substrate in WP/G/Spring, inactive Rac divided by 7.5 in Holmes, and 0.04/(1+Rho²) in deBelly. Native channels and clocks are not calibrated between models.

The common perturbation removes epsilon = 0.005 times a smooth, fixed basal-incumbent excess profile, with 0.001/0.01 amplitude checks. Explicit inactive pools receive that material where present; deBelly changes Rac alone while MCA/mechanics/controller states start unchanged. Every native field then evolves. Within the frozen original core,

\\[
D(t)=A_{{\\mathrm{{sham}}}}(t)-A_{{\\mathrm{{probe}}}}(t),
\\qquad K_{{\\tau}}=\\frac{{D(0)-D(\\tau)}}{{\\tau D(0)}}.
\\tag{{1}}
\\]

All 97 measurements qualify at their fixed model clocks under the reviewed empirical numerical, amplitude-linearity and unsaturated-prefix gates. Actual sham/probe changes and rates, core/whole/outside deficits and every gate remain in the [recovery summary](supporting/data/vision_linear_recovery_summary.csv) and [complete audit](../../../../data/runs/vision_linear_20261017/recovery/results/verification.json). The positive coefficient concerns this activity-loss direction, not general stability or resistance to an actual challenger. In deBelly it does not measure the difficulty of rewriting MCA/mechanical organization, which the kick leaves initially intact.

## Substantive finite competition and withdrawal

![All-five finite competition: global activity bias and retained local organization](vision_linear_finite_retained_detail.png)

**Why this figure.** Native first moments give a bounded signed activity summary during and after a substantial cue. The second row shows the local regional amounts that a first moment can obscure: near-zero polarity can accompany two active regions. All five models and 97 qualified sizes remain visible. Each row has a common scale. The [full-range cue-on/core comparison](vision_linear_finite_competition.png) retains deBelly's large direct cue response.

The right-end cosine has fixed support 3 (Spring 0.2) and unchanged model-native peaks 0.1/6.4/3.6/1/4 for WP/G/deBelly/Holmes/Spring. Primary duration is 1000 native units, Spring 50, followed by 13,000 cue-free units. These are standardized within-model exposures, not universal commitment or physical-time calibration. The prior fixed-support reference strength clocks are 73.5, 377.5, 335.79, 2360 and 0.31014 respectively: even Holmes' reference does not need to commit within 1000. Half/double-duration checks at the prespecified 15 short/middle/long anchor states are retained; amplitudes and durations were fixed before inspecting the new histories.

| Model | Qualified sizes | 10,000 off: orientation / strength / field flag | 13,000 off: orientation / strength / field flag |
|---|---:|---|---|
{chr(10).join(finite_lines)}

The new full native histories reproduce every existing primary 10,000-unit scalar and field flag. Continuing to 13,000 leaves those counts unchanged, while the earlier 3000-unit window can differ: deBelly has 12 orientation, 5 strength and 0 field passes then, compared with 12/11/12 at 8000–13,000. A rightward scalar bias is therefore distinct from original strength and full-field stationarity. These finite-window flags do not prove attractors.

At the largest intervals, WP L39 and Holmes L50 can retain two threshold-defined active regions after withdrawal; Goryachev L39 instead restores its original incumbent despite challenger activity during the cue. All deBelly/Spring primary cases end right-oriented. deBelly L18's mirrored-core excess is about 7.99 initial-incumbent units during the cue and 1.013 after withdrawal. Its nonconserved Rac source and slower native fields make cue-on amplitude a different observation from retained organization. The complete nonmonotone/null size responses and unreached crossing clocks are preserved in the [finite table](supporting/data/vision_linear_finite_metrics.csv), [declared time samples](supporting/data/vision_linear_finite_samples.csv) and [finite report](supporting/metadata/FINITE_RESULTS.md).

## Measurements motivated by these results

| Theoretical observation | Corresponding measurement question |
|---|---|
| Incumbent material grows through different peak/width changes. | How do reporter peak, active-region extent and integrated signal change with cell size? Molecular amount requires calibration. |
| Equal native local input can produce different early activation. | Under a fixed measured local stimulus, how fast does a competing activity region initially grow? |
| A specified local activity deficit can recover without rewriting all slower fields. | After brief local suppression, measure local loss/rebound alongside the incumbent's other organizational features. |
| Cue-on competition and post-cue survival differ; two regions can cancel a global bias. | Track incumbent/challenger activity, extinction/merger or coexistence, and global activity-polarity magnitude/direction during and after a finite cue. |

WGD-associated enlarged-cell persistence motivates these measurements. Neither a first moment nor a local recovery coefficient is a cell velocity law, and these simulations do not establish that size alone causes the migration observations.

## Evidence, numerical scope and completion

The [supplied mechanism PDF](../../length_mechanisms.pdf) already establishes incumbent enlargement and selected finite collapse/retained-reversal examples. The previous canonical duration suite, identified in the [long-cue inventory](supporting/metadata/LONG_CUE_INVENTORY.md), already contains the primary scalar branches used here. This campaign adds the complete common active-loss recovery series, independently audited fixed-exposure native histories/regions and 13,000-unit comparison. It does not claim every finite fate as a new discovery or import ring results into linear evidence.

All 97 qualified exact checkpoints are retained, plus 13 preparation-ineligible rows in finite coverage. C reuses 15 exact reviewed pilot archives and adds 82 missing cases. Every main measurement uses its completed original full-span archive. The slow Spring L₀0.8 original finished before any substitution or cancellation; its later candidates are measured. A separately frozen, root-authorized early-window run at L₀0.8 and its adjacent L₀0.9 validated the effect of changing the adaptive solver's final time. It passes the unchanged sensitivity allowance and remains additional validation, unused by the main collection. Only that separate run has unmeasured later candidates. See the [amendment definition](supporting/metadata/EARLY_WINDOW_AMENDMENT.md), [fresh audit and original-versus-amendment comparisons](supporting/metadata/EARLY_WINDOW_REVIEW.md) and [amendment results](../../../../data/runs/vision_linear_20261017/recovery/early_window_amendment/results/verification.json).

Goryachev L39's original core has four native cell centers; native solver/reference checks do not establish spatial convergence. Core membership/background stay frozen, raw and signed excess are separate, and outside-core transport is not relabelled whole-activity restoration. Spring uses each evolving physical length; deBelly keeps independent endpoints, native 0.5-unit stepping and all 11 fields/controllers. Native full state is continuous at withdrawal; cue-dependent derived rates can jump appropriately.

- [Design and approved protocol](DESIGN.md), [fresh design review](supporting/metadata/DESIGN_REVIEW.md), and [original long-cue inventory](supporting/metadata/LONG_CUE_INVENTORY.md).
- [Recovery inputs/producer](../../../../data/runs/vision_linear_20261017/recovery/manifest.json), [all recovery curves/gates](../../../../data/runs/vision_linear_20261017/recovery/results/), and [finite inputs/producer](../../../../data/runs/vision_linear_20261017/finite/manifest.json).
- [Finite archive and scalar audit](../../../../data/runs/vision_linear_20261017/finite/results/verification.json); every one of the 97 primary and 30 anchor branches is complete, with 13 unavailable preparations visible.
- [Compact normalized functional table](supporting/data/vision_linear_normalized_function.csv) and [gain source](supporting/data/vision_linear_gain_all97.csv). PNG figures have PDF/SVG counterparts under supporting/vectors.

The [fresh actual-result critique](supporting/metadata/FRESH_RESULT_CRITIQUE.md) is the independent review of these actual outputs. This evidence is supplied to the checkpoint assembler; no further size/duration search is launched by report generation.
'''
    path = OUT/"RESULTS.md"; path.write_text(text); FILES.append(path)
    source_files = [Path(__file__), GAIN, FOLDS, RUN/"recovery/results/recovery_summary.csv", RUN/"recovery/results/verification.json",
                    RUN/"finite/results/metrics.csv", RUN/"finite/results/samples.csv", RUN/"finite/results/verification.json",
                    OUT/"supporting/metadata/FINITE_RENDER_MANIFEST.json"]
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), artifact_root="../..",
        source_sha256={str(p.relative_to(ROOT.parents[1])): sha(p) for p in source_files},
        sha256={str(p.relative_to(OUT)): sha(p) for p in FILES}, new_integrations=0)
    (OUT/"supporting/metadata/vision_linear_report.manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(dict(report=str(path), recovery_cases=len(recovery), finite_rows=len(finite), exports=len(FILES)), indent=2))


if __name__ == "__main__": main()
