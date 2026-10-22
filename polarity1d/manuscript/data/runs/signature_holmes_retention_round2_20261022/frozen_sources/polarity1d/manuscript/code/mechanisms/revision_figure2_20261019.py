"""Reanalyse retained deterministic evidence for revision-program Figure 2.

No integration and no writes to retained evidence. Geometry/stochastic work is
owned by separate revision programs. Regenerate with the existing environment.
"""
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
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
DEV = ROOT / "figures/length_mechanisms_20261008/dev"
OUT = DEV / "revision_program_20261019/figure2"
DATA = OUT / "supporting/data"
META = OUT / "supporting/metadata"
VECTORS = OUT / "supporting/vectors"
FINITE = ROOT / "data/runs/vision_linear_20261017/finite"
RECOVERY = ROOT / "data/runs/vision_linear_20261017/recovery/results"
CONTROLS = ROOT / "data/runs/three_test_mechanisms_20261006"
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
NAMES = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "deBelly", "Holmes 3", "Spring")))
COLORS = dict(zip(MODELS, ("#3976a8", "#c96926", "#2a9288", "#8560a4", "#a34f7f")))
MARKERS = dict(zip(MODELS, ("o", "s", "^", "D", "v")))
REF = dict(zip(MODELS, (19., 19., 19., 20., 1.)))
PAIRS = dict(zip(MODELS, ((19., 39.), (19., 39.), (19., 39.), (20., 40.), (1., 2.))))
FIELDS = dict(zip(MODELS, ("active", "active", "rac", "rac_active", "active")))
FATE_COLORS = {"challenger only": "#347f96", "both cores": "#ac7d24", "incumbent only": "#9b465a", "neither core": "#999999"}
SOURCES = {}
ARTIFACTS = []
CHECKS = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def source(path, expected=None):
    path = Path(path).resolve()
    digest = sha(path)
    if expected is not None:
        assert digest == expected, path
    SOURCES[str(path.relative_to(REPO))] = digest
    return path


def rows(path):
    with source(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def read(path):
    return json.loads(source(path).read_text())


def flag(value):
    return value is True or value == "True"


def export(name, records):
    path = DATA / name
    fields = list(dict.fromkeys(k for r in records for k in r))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    ARTIFACTS.append(path)


def write(path, text):
    path.write_text(text)
    ARTIFACTS.append(path)


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        path = (OUT if ext == "png" else VECTORS) / f"{stem}.{ext}"
        fig.savefig(path, dpi=210 if ext == "png" else None, bbox_inches="tight")
        ARTIFACTS.append(path)
    plt.close(fig)


def close(actual, expected, rtol=3e-10, atol=3e-12):
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)


def function_data():
    folds = rows(DEV / "integrated_framework_candidate_20261017/supporting/data/all_five_incumbent_size_folds.csv")
    profiles = rows(DEV / "integrated_framework_candidate_20261017/supporting/data/all_five_native_basal_profiles.csv")
    gains = rows(DEV / "all_five_linear_function_20261017/supporting/data/native_gain_design_evidence.csv")
    gain_identity = read(DEV / "all_five_linear_function_20261017/supporting/metadata/native_gain_design_evidence.json")
    assert gain_identity["retained_checkpoint_hashes_checked"] == 97
    assert gain_identity["maximum_native_increment_identity_error"] < 1e-10
    recovery = rows(RECOVERY / "recovery_summary.csv")
    verify = read(RECOVERY / "verification.json")
    assert verify["passed"] and verify["completed_cases"] == 97
    assert verify["all_fixed_clocks_qualified"] and verify["zero_kick_all_fields_identity"]
    assert len(gains) == len(recovery) == 97 and len(folds) == 5
    normalized, summary = [], []
    for model in MODELS:
        gr = [r for r in gains if r["model"] == model]
        rr = [r for r in recovery if r["model"] == model]
        g0 = next(r for r in gr if float(r["size_native"]) == REF[model])
        r0 = next(r for r in rr if float(r["size_native"]) == REF[model])
        for g in gr:
            r = next(r for r in rr if r["size_native"] == g["size_native"])
            assert flag(r["coefficient_available"])
            close(float(r["restored_fraction"]), float(r["horizon_native"]) * float(r["conditional_K_per_native_time"]), rtol=3e-14, atol=3e-17)
            checkpoint = source(ROOT / g["checkpoint"], g["checkpoint_sha256"])
            assert g["checkpoint_sha256"] == r["checkpoint_sha256"]
            # Recompute the cue-weighted native channel from complete basal fields.
            with np.load(checkpoint) as z:
                fields = {k[6:]: z[k][-1] for k in z.files if k.startswith("field_")}
                x = z["x"]
                physical = float(g["physical_length"])
                if model == "spring":
                    x = np.linspace(0, physical, len(fields["active"]))
                    weights = np.full(len(x), physical/(len(x)-1)); weights[[0, -1]] *= .5
                else:
                    weights = np.full(len(x), physical/len(x))
                coefficient = .04/(1+fields["rho"]**2) if model == "debelly" else fields["inactive"] if model == "spring" else z["state"][1]/(7.5 if model == "holmes_model3" else 1.)
                width = .2 if model == "spring" else 3.
                distance = physical-x
                cue = np.where(distance <= width, .5*(1+np.cos(np.pi*np.minimum(distance/width, 1))), 0.)
                computed = float((cue*coefficient) @ weights / (cue @ weights))
                close(computed, float(g["cue_weighted_gain_native"]), rtol=1e-12, atol=1e-13)
            normalized.append(dict(model=model, size_native=float(g["size_native"]),
                size_over_reference=float(g["size_native"])/REF[model],
                gain_native=computed, gain_over_reference=computed/float(g0["cue_weighted_gain_native"]),
                restored_fraction=float(r["restored_fraction"]),
                repair_over_reference=float(r["restored_fraction"])/float(r0["restored_fraction"]),
                clock_native=float(r["horizon_native"]), checkpoint=g["checkpoint"], checkpoint_sha256=g["checkpoint_sha256"]))
        pair = [next(r for r in normalized if r["model"] == model and r["size_native"] == size) for size in PAIRS[model]]
        f = next(r for r in folds if r["model"] == model)
        for key in ("excess_active", "width", "peak_concentration"):
            close(float(f[key+"_fold"]), float(f[key+"_large"])/float(f[key+"_small"]))
        summary.append(dict(model=model, size_small=PAIRS[model][0], size_large=PAIRS[model][1],
            amount_fold=float(f["excess_active_fold"]), width_fold=float(f["width_fold"]),
            peak_fold=float(f["peak_concentration_fold"]),
            gain_fold=pair[1]["gain_over_reference"]/pair[0]["gain_over_reference"],
            repair_fold=pair[1]["repair_over_reference"]/pair[0]["repair_over_reference"]))
    CHECKS["basal_checkpoints_hash_verified_and_gains_recomputed"] = 97
    CHECKS["retained_native_increment_identity_max_error"] = gain_identity["maximum_native_increment_identity_error"]
    CHECKS["qualified_recovery_identity_and_full_campaign_audit"] = 97
    export("organization_size_pairs.csv", summary)
    export("native_basal_profiles.csv", profiles)
    export("all97_gain_recovery.csv", normalized)
    return folds, profiles, normalized, summary


def fate(old, challenger, threshold=1.):
    return ("both cores" if old >= threshold and challenger >= threshold else
            "incumbent only" if old >= threshold else
            "challenger only" if challenger >= threshold else "neither core")


def finite_data():
    original = rows(FINITE / "results/metrics.csv")
    verify = read(FINITE / "results/verification.json")
    assert verify["passed"] and verify["complete_primary"] == 97 and verify["complete_sensitivity"] == 30
    result, sensitivity, displays, sample_rows = [], [], {}, []
    for row in original:
        if row["status"] != "complete":
            continue
        archive = source(FINITE / row["archive"], row["archive_sha256"])
        with np.load(archive) as z:
            field = FIELDS[row["model"]]
            cue, post = z["cue_field_"+field], z["post_field_"+field]
            np.testing.assert_array_equal(cue[-1], post[0])
            b = float(row["basal_background"])
            scale = float(row["initial_halfcontrast_threshold"])-b
            physical = float(row["physical_length"])
            q = z["material_coordinate"] if row["model"] == "spring" else z["x"]/float(row["size_native"])
            edge = float(row["core_width_native"])/physical
            old = q <= edge; new = q >= 1-edge
            cue_peak = (cue[:, new].max(axis=1)-b)/scale
            post_old = (post[:, old].max(axis=1)-b)/scale
            post_new = (post[:, new].max(axis=1)-b)/scale
            close(cue_peak, z["cue_regional_challenger_peak_threshold_ratio"])
            close(post_old, z["post_regional_old_peak_threshold_ratio"])
            close(post_new, z["post_regional_challenger_peak_threshold_ratio"])
            close(post_old[-1], float(row["post13000_old_peak_threshold_ratio"]))
            close(post_new[-1], float(row["post13000_challenger_peak_threshold_ratio"]))
            # Independent regional amount calculation, including evolving Spring length.
            for phase, u in (("cue", cue), ("post", post)):
                length = z[phase+"_field_length"] if row["model"] == "spring" else physical
                weights = z["core_weights_material"]
                amount = (u-b) @ weights[::-1] * length / float(row["initial_incumbent_excess"])
                close(amount, z[phase+"_regional_challenger_excess_ratio"])
            r = dict(model=row["model"], size_native=float(row["size_native"]),
                cue_duration_native=float(row["cue_duration_native"]), primary=flag(row["primary"]),
                cue_core_overlap=flag(row["cue_core_overlap"]), anchor=flag(row["anchorflag"]),
                cue_challenger_max_halfcontrast_ratio=float(cue_peak.max()),
                cue_challenger_end_halfcontrast_ratio=float(cue_peak[-1]),
                post13000_old_halfcontrast_ratio=float(post_old[-1]),
                post13000_challenger_halfcontrast_ratio=float(post_new[-1]),
                post13000_challenger_excess_over_initial_incumbent=float(z["post_regional_challenger_excess_ratio"][-1]),
                post13000_old_excess_over_initial_incumbent=float(z["post_regional_old_excess_ratio"][-1]),
                final_core_geography=fate(post_old[-1], post_new[-1]),
                above_halfcontrast_during_cue=bool(cue_peak.max() >= 1.),
                orientation_success=flag(row["post13000_orientation_success"]),
                field_tail_pass=flag(row["post13000_field_tail_pass"]),
                connected_regions=int(round(float(row["post13000_sites_above_initial_halfcontrast"]))),
                archive=str(archive.relative_to(REPO)), archive_sha256=row["archive_sha256"])
            result.append(r)
            for threshold in (.75, 1., 1.25):
                sensitivity.append(dict(model=r["model"], size_native=r["size_native"], primary=r["primary"],
                    cue_duration_native=r["cue_duration_native"], halfcontrast_ratio_threshold=threshold,
                    cue_activation=bool(cue_peak.max() >= threshold),
                    final_core_geography=fate(post_old[-1], post_new[-1], threshold)))
            selected = ((row["model"] == "wave_pinning" and float(row["size_native"]) == 39. and flag(row["primary"])) or
                        (row["model"] == "debelly" and float(row["size_native"]) == 39. and float(row["cue_duration_native"]) in (500., 1000.)))
            largest = {"wave_pinning":39., "goryachev":39., "debelly":39., "holmes_model3":50., "spring":2.6}
            if selected or (flag(row["primary"]) and float(row["size_native"]) == largest[row["model"]]):
                key = f'{row["model"]}_{r["size_native"]:g}_{r["cue_duration_native"]:g}'
                T = r["cue_duration_native"]
                t = np.r_[z["cue_time"], T+z["post_time"][1:]]
                u = np.concatenate((cue, post[1:]))
                displays[key] = dict(time=t, field=u, fraction=q, row=r,
                    old_peak=np.r_[(cue[:, old].max(axis=1)-b)/scale, post_old[1:]],
                    new_peak=np.r_[cue_peak, post_new[1:]])
                for phase in ("cue", "post"):
                    for i, time in enumerate(z[phase+"_time"]):
                        sample_rows.append(dict(example=key, phase=phase, time_native=float(time),
                            old_peak_halfcontrast_ratio=float(z[phase+"_regional_old_peak_threshold_ratio"][i]),
                            challenger_peak_halfcontrast_ratio=float(z[phase+"_regional_challenger_peak_threshold_ratio"][i]),
                            old_excess_ratio=float(z[phase+"_regional_old_excess_ratio"][i]),
                            challenger_excess_ratio=float(z[phase+"_regional_challenger_excess_ratio"][i])))
    assert len(result) == 127 and sum(r["primary"] for r in result) == 97
    export("all127_native_activation_and_fate.csv", result)
    export("threshold_sensitivity.csv", sensitivity)
    export("selected_native_regional_histories.csv", sample_rows)
    export("preparation_unavailable.csv", [r for r in original if r["status"] != "complete"])
    CHECKS["full_finite_archive_sha256_and_native_peak_amount_recalculation"] = len(result)
    CHECKS["native_withdrawal_identity"] = len(result)
    return result, sensitivity, displays


def control_data():
    manifest = read(CONTROLS / "manifest.json")
    assert read(CONTROLS / "results/status.json")["complete"]
    source(CONTROLS / "results/cell_evidence.json")
    table = rows(CONTROLS / "results/trajectory_metrics.csv")
    # No legacy deBelly trajectory is admitted: corrected evidence comes below.
    selected = [r for r in table if r["model"] == "spring" and r["kind"] in ("feedback_off_during_pulse", "matched_gradient_peak")]
    selected += [r for r in table if r["model"] in ("wave_pinning", "holmes_model3") and r["geometry"] == "fixed_physical" and
                 float(r["length"]) in ({39.} if r["model"] == "wave_pinning" else {28., 50.}) and
                 float(r["duration"]) == 5120. and r["kind"] in ("native", "mobility_4x", "mesh_401")]
    controls = []
    for row in selected:
        slot = int(row["slot"])
        meta = read(CONTROLS / f"trajectories/trajectory_{slot:04d}.json")
        archive = source(CONTROLS / f"trajectories/trajectory_{slot:04d}.npz", meta["archive_sha256"])
        source(ROOT / meta["job"]["case"]["checkpoint"], meta["checkpoint_sha256"])
        with np.load(archive) as z:
            tail = z["recovery_time"] >= 9500-1e-8
            ori, contrast = z["recovery_orientation"][tail], z["recovery_contrast"][tail]
            success = bool(ori.min() > .05 and contrast.min() >= .05 and np.ptp(ori) < .005)
            assert success == meta["metrics"]["orientation_success"]
            initial_magnitude = abs(float(meta["job"]["case"]["initial_orientation"]))
            strength_error = float(np.max(abs(ori/initial_magnitude-1.)))
            strength_success = bool(success and strength_error <= .01)
            assert strength_success == meta["metrics"]["full_polarity_success"]
            close(ori[-1], meta["metrics"]["final_orientation"])
        controls.append(dict(model=row["model"], geometry=row["geometry"], size_native=float(row["length"]),
            cue_duration_native=float(row["duration"]), intervention=row["kind"],
            orientation_success=success, final_orientation=float(ori[-1]),
            full_polarity_success=strength_success, maximum_relative_strength_error=strength_error,
            initial_scalar_magnitude=initial_magnitude,
            field_tail_pass=meta["metrics"]["field_tail_pass"], native_recovery_horizon=10000,
            archive=str(archive.relative_to(REPO)), archive_sha256=meta["archive_sha256"]))
    assert len(controls) == 22
    CHECKS["selected_evolved_control_gates_recomputed_from_native_histories"] = len(controls)
    wp = {r["intervention"]:r for r in controls if r["model"] == "wave_pinning"}
    assert not wp["native"]["orientation_success"] and wp["mobility_4x"]["orientation_success"]
    assert not wp["native"]["full_polarity_success"] and not wp["mobility_4x"]["full_polarity_success"]
    CHECKS["WP_mobility_endpoint_scope"] = dict(
        intervention_scope="Fourfold inactive diffusion throughout cue and recovery",
        same_fixed_support_pulse_native=5120,
        native_orientation=wp["native"]["final_orientation"],
        mobility4x_orientation=wp["mobility_4x"]["final_orientation"],
        native_orientation_gate_pass=False, mobility4x_orientation_gate_pass=True,
        native_strength_gate_pass=False, mobility4x_strength_gate_pass=False)
    export("selected_evolved_interventions.csv", controls)
    bolus = rows(ROOT / "data/runs/inactive_pool_figure_exports_20261008/response_vs_length.csv")
    source(ROOT / "data/runs/inactive_pool_figure_exports_20261008/README.md")
    selected_bolus = [r for r in bolus if float(r["duration"]) == .1 and float(r["epsilon"]) == .005]
    assert len(selected_bolus) == 36
    export("separate_inactive_bolus.csv", selected_bolus)
    ablation_root = ROOT / "data/runs/length_mechanisms_20261008/narrative_20261009/resting_ablation"
    ablation = rows(ablation_root / "ablation_endpoint_metrics.csv")
    validation = read(ablation_root / "ablation_validation.json")
    assert validation["passed"] and len(validation["checks"]) == 6
    assert all(r["verified_steady_state"] and r["exact_restart_verified"] and r["native_trace_hash_verified"] for r in validation["checks"])
    for condition in ("native", "ablation"):
        values = [abs(float(r["left_boundary_displacement"])-float(r["right_boundary_displacement"])) for r in ablation if r["condition"] == condition]
        close(max(values)-min(values), validation["displacement_spread_um"][condition])
    CHECKS["corrected_resting_ablation_validation"] = validation
    export("corrected_debelly_resting_ablation.csv", ablation)
    return controls, selected_bolus, ablation


def reversal_data():
    root = ROOT / "data/runs/debelly_corrected_three_test_20261008/results"
    validation = read(root / "validation.json")
    assert validation["passed"] and validation["independent_audit_passed"]
    source(root / "independent_audit.json", validation["independent_audit_sha256"])
    endpoints = rows(root / "endpoints.csv")
    assert len(endpoints) == 660
    selected = [r for r in endpoints if r["model"] in ("wave_pinning", "goryachev") and
                r["condition"] == "fixed_physical" and r["definition"] == "full_polarity"]
    export("lead_pair_fixed_support_strength_thresholds.csv", selected)
    # Retain all five and every cue class as scope/context, including nonmonotone deBelly.
    export("all_five_three_cue_strength_thresholds.csv", [r for r in endpoints if r["definition"] == "full_polarity"])
    CHECKS["canonical_combined_endpoint_table_rows"] = len(endpoints)
    CHECKS["lead_pair_successful_duration_strictly_increasing"] = {}
    for model in ("wave_pinning", "goryachev"):
        group = sorted([r for r in selected if r["model"] == model and r["status"] == "success"], key=lambda r:float(r["length"]))
        assert np.all(np.diff([float(r["endpoint"]) for r in group]) > 0)
        CHECKS["lead_pair_successful_duration_strictly_increasing"][model] = dict(
            successful_lengths=len(group), minimum_length=float(group[0]["length"]), maximum_length=float(group[-1]["length"]),
            first_duration=float(group[0]["endpoint"]), last_duration=float(group[-1]["endpoint"]),
            censored_lengths=sum(r["status"] == "censored" for r in selected if r["model"] == model))
    return selected


def claim_figure(endpoints, normalized, finite, displays):
    fig = plt.figure(figsize=(11.8, 9.4))
    gs = fig.add_gridspec(2, 3, height_ratios=(1., 1.04), left=.075, right=.985, top=.865, bottom=.20, wspace=.40, hspace=.65)
    ax = fig.add_subplot(gs[0, 0])
    for model in ("wave_pinning", "goryachev"):
        group = sorted([r for r in endpoints if r["model"] == model], key=lambda r:float(r["length"]))
        ref = float(next(r for r in group if float(r["length"]) == 19.)["endpoint"])
        valid = [r for r in group if r["status"] == "success"]
        ax.plot([float(r["length"])/19. for r in valid], [float(r["endpoint"])/ref for r in valid], color=COLORS[model], lw=1.6)
        for r in valid:
            ax.scatter(float(r["length"])/19., float(r["endpoint"])/ref, s=25, facecolor=COLORS[model] if flag(r["selected_field_tail_pass"]) else "white", edgecolor=COLORS[model], linewidth=.9, zorder=4)
        for r in group:
            if r["status"] == "censored":
                ax.scatter(float(r["length"])/19., float(r["maximum_tested_duration"])/ref, marker="^", s=32, facecolor="white", edgecolor=COLORS[model], lw=1)
    ax.set_yscale("log")
    ax.set(xlabel="Interval length / 19", ylabel="Shortest strength-passing duration\n/ own L19 duration", xlim=(0, 2.15),
        title="A  Similar qualitative reversal resistance")
    ax.grid(axis="y", color="#eeeeee", lw=.7)
    ax.text(.03, .98, "Open circles: field flag fails\nOpen triangles: tested cap", transform=ax.transAxes, va="top", fontsize=7.4)
    for col, key, title in ((1, "gain_over_reference", "B  Different challenger responsiveness"), (2, "repair_over_reference", "C  Different incumbent restoration")):
        ax = fig.add_subplot(gs[0, col])
        for model in MODELS:
            group = sorted([r for r in normalized if r["model"] == model], key=lambda r:r["size_native"])
            lead = model in ("wave_pinning", "goryachev")
            ax.plot([r["size_over_reference"] for r in group], [r[key] for r in group], color=COLORS[model],
                lw=1.6 if lead else .85, alpha=1 if lead else .35, marker="o", ms=2.8 if lead else 2, zorder=3 if lead else 1)
            if lead:
                pair = [r for r in group if r["size_native"] in PAIRS[model]]
                ax.scatter([r["size_over_reference"] for r in pair], [r[key] for r in pair], s=32, facecolor="white", edgecolor=COLORS[model], zorder=5, lw=1)
                last = pair[-1]
                ax.annotate(f'{last[key]:.3g}×', (last["size_over_reference"], last[key]), xytext=(4, 5 if model == "wave_pinning" else -12),
                    textcoords="offset points", color=COLORS[model], fontsize=8.3)
        ax.set(xlabel="Native size / own reference", ylabel="Value / own reference", xlim=(0, 2.7), xticks=(0, 1, 2), title=title)
        ax.set_yscale("log", base=2)
        ax.axhline(1, color="#aaaaaa", lw=.7, ls="--")
        ax.grid(axis="y", color="#eeeeee", lw=.7)
    for col, model, title in ((0, "wave_pinning", "D  A robust challenger; two retained regions"), (1, "goryachev", "E  A weaker challenger; incumbent return")):
        ax = fig.add_subplot(gs[1, col])
        d = displays[f"{model}_39_1000"]
        keep = d["time"] <= 4000
        im = ax.pcolormesh(d["fraction"], d["time"][keep]/1000., d["field"][keep], shading="auto", cmap="magma", vmin=0, rasterized=True)
        ax.axhline(1, color="#26b9bd", lw=1.2)
        ax.set(xlabel="Fractional interval position", ylabel="Time / cue duration", title=title+f"\n{NAMES[model]} L39; cue 1000", ylim=(4, 0))
        if model == "goryachev":
            i = int(np.argmin(abs(d["time"]-1000.)))
            b, scale = float(d["field"][0].min()), .5*float(np.ptp(d["field"][0]))
            ax.text(.97, .88, f'At withdrawal\nChallenger peak {b+scale*d["new_peak"][i]:.2f}\nIncumbent peak {b+scale*d["old_peak"][i]:.2f}',
                transform=ax.transAxes, ha="right", va="top", fontsize=7.7, color="white")
        cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.025); cb.ax.tick_params(labelsize=7)
        cb.set_label("Native active concentration", fontsize=7.8)
    ax = fig.add_subplot(gs[1, 2])
    for r in finite:
        ax.scatter(r["cue_challenger_max_halfcontrast_ratio"], r["post13000_challenger_excess_over_initial_incumbent"],
            s=26 if r["primary"] else 15, marker=MARKERS[r["model"]], facecolors=FATE_COLORS[r["final_core_geography"]] if r["primary"] else "none",
            edgecolors=FATE_COLORS[r["final_core_geography"]], linewidth=.65, alpha=.75)
    ax.set_xscale("log")
    ax.axvline(1, color="#aaaaaa", lw=.8, ls="--")
    ax.set(xlabel="Maximum challenger during cue\n/ initial half-contrast above background", ylabel="Challenger excess after 13,000 off\n/ initial incumbent core excess", title="F  Activation alone does not establish fate",
        ylim=(-.065, 1.12), xlim=(.075, 55.))
    ax.legend(handles=[Patch(facecolor=FATE_COLORS[k], label=k.capitalize()) for k in ("challenger only", "both cores", "incumbent only")], frameon=False, loc="upper left", fontsize=7.4)
    ax.text(.98, .07, "97 primary exposures\n+ 30 duration controls", transform=ax.transAxes, ha="right", fontsize=7.4)
    fig.legend(handles=[Line2D([], [], color=COLORS[m], lw=1.7 if m in ("wave_pinning", "goryachev") else 1.,
        alpha=1 if m in ("wave_pinning", "goryachev") else .4, label=NAMES[m]) for m in MODELS],
        ncol=5, loc="upper center", bbox_to_anchor=(.53, .954), frameon=False, fontsize=8.2)
    fig.suptitle("Figure 2 | Similar reversal trends, distinct challenge responses", fontsize=13., y=.995)
    fig.text(.53, .055, "A: fixed-support scalar-strength thresholds; similar means increasing over observed successes, not equal slopes or calibrated times.\n"
        "B–C: all 97 basals; WP/Goryachev highlighted. C is a specified local active-loss response, not general stability.\n"
        "D–E: exact native fields, model-specific zero-based scales; cyan = withdrawal. F: final fixed-core geography; filled = primary, open = duration control.",
        ha="center", fontsize=8.0, linespacing=1.5)
    save(fig, "figure2_mechanisms")


def measurement_context(folds, normalized, finite, displays):
    fig = plt.figure(figsize=(12.3, 9.1))
    gs = fig.add_gridspec(2, 3, height_ratios=(1., 1.08), left=.075, right=.985, top=.88, bottom=.22, wspace=.39, hspace=.61)
    ax = fig.add_subplot(gs[0, 0])
    for field, marker, label, delta in (("excess_active_fold", "o", "Amount", -.16), ("width_fold", "s", "Width", 0), ("peak_concentration_fold", "^", "Peak", .16)):
        ax.scatter([float(r[field]) for r in folds], np.arange(5)+delta, marker=marker, s=27,
            color="#333333" if marker == "o" else "white", edgecolor="#333333", lw=.8, label=label)
    ax.set(yticks=range(5), yticklabels=[NAMES[r["model"]] for r in folds], ylim=(4.5, -.5), xlim=(.8, 3.18),
        xlabel="Larger / smaller value", title="A  More incumbent material; different organization")
    ax.axvline(1, color="#aaaaaa", lw=.7, ls="--")
    ax.legend(frameon=False, fontsize=7.8, loc="upper left", bbox_to_anchor=(-.12, 1.26), ncol=3)
    for col, key, title in ((1, "gain_over_reference", "B  Initial opposing-input gain"), (2, "repair_over_reference", "C  Local deficit restoration")):
        axis = fig.add_subplot(gs[0, col])
        for model in MODELS:
            group = sorted([r for r in normalized if r["model"] == model], key=lambda r:r["size_native"])
            axis.plot([r["size_over_reference"] for r in group], [r[key] for r in group], color=COLORS[model], lw=1.4, marker="o", ms=2.6)
            pair = [r for r in group if r["size_native"] in PAIRS[model]]
            axis.scatter([r["size_over_reference"] for r in pair], [r[key] for r in pair], s=26, facecolor="white", edgecolor=COLORS[model], zorder=4, lw=.8)
        axis.set(xlabel="Native size / own reference", ylabel="Value / own reference", xlim=(0, 2.7), xticks=(0, 1, 2), title=title)
        axis.set_yscale("log", base=2)
        axis.axhline(1, color="#aaaaaa", lw=.7, ls="--")
        axis.grid(axis="y", color="#eeeeee", lw=.7)
    fig.legend([Line2D([], [], color=COLORS[m], lw=1.6, label=NAMES[m]) for m in MODELS],
        [NAMES[m] for m in MODELS], ncol=5, loc="upper center", bbox_to_anchor=(.61, .962), frameon=False, fontsize=8.2)
    ax = fig.add_subplot(gs[1, 0])
    for r in finite:
        ax.scatter(r["cue_challenger_max_halfcontrast_ratio"], r["post13000_challenger_excess_over_initial_incumbent"],
            s=27 if r["primary"] else 15, marker=MARKERS[r["model"]],
            facecolors=FATE_COLORS[r["final_core_geography"]] if r["primary"] else "none",
            edgecolors=FATE_COLORS[r["final_core_geography"]], linewidth=.65, alpha=.75, zorder=3)
    ax.set_xscale("log")
    ax.set(xlabel="Maximum challenger peak during cue\n/ initial half-contrast above background",
        ylabel="Challenger core excess after 13,000 off\n/ initial incumbent core excess",
        title="D  Strong activation has different later fates", ylim=(-.065, 1.12), xlim=(.075, 55.))
    ax.axvline(1, color="#aaaaaa", lw=.8, ls="--")
    ax.legend(handles=[Patch(facecolor=FATE_COLORS[k], label=k.capitalize()) for k in ("challenger only", "both cores", "incumbent only")],
        frameon=False, loc="upper left", fontsize=7.5)
    ax.text(.98, .07, "Filled: 97 primary exposures\nOpen: 30 duration controls", transform=ax.transAxes, ha="right", fontsize=7.6)
    ax = fig.add_subplot(gs[1, 1])
    d = displays["wave_pinning_39_1000"]
    keep = d["time"] <= 4000
    image = ax.pcolormesh(d["fraction"], d["time"][keep]/1000, d["field"][keep], shading="auto", cmap="magma", vmin=0, rasterized=True)
    ax.axhline(1, color="#26b9bd", lw=1.2)
    ax.set(xlabel="Fractional interval position", ylabel="Time / cue duration", title="E  Two regions survive withdrawal", ylim=(4, 0))
    ax.text(.5, -.23, "Wave-pinning L39; 1000-unit cue", transform=ax.transAxes, ha="center", fontsize=8.2)
    cb = fig.colorbar(image, ax=ax, fraction=.045, pad=.025)
    cb.ax.tick_params(labelsize=7); cb.set_label("Native active concentration", fontsize=8)
    ax = fig.add_subplot(gs[1, 2])
    for T, ls in ((500, "--"), (1000, "-")):
        d = displays[f"debelly_39_{T}"]
        for key, color, name in (("old_peak", "#325c9b", "Incumbent"), ("new_peak", "#b74b32", "Challenger")):
            ax.plot(d["time"], d[key], ls=ls, lw=1.35, color=color, label=f"{name}; cue {T}")
        ax.axvline(T, color="#555555", lw=.7, ls=ls)
    ax.axhline(1, color="#aaaaaa", lw=.7, ls=":")
    ax.set(xlabel="Time from cue onset (native seconds)", ylabel="Core peak / initial half-contrast\nabove background", title="F  Withdrawal reveals incumbent return", ylim=(-.2, 23), xlim=(0, 14000))
    ax.set_xscale("symlog", linthresh=1000)
    ax.set_xticks([0, 500, 1000, 3000, 13000], ["0", "500", "1000", "3000", "13,000"])
    ax.tick_params(axis="x", labelsize=7)
    ax.legend(frameon=False, fontsize=7.1, loc="upper right")
    ax.text(.5, -.23, "Corrected deBelly L39; same complete basal state", transform=ax.transAxes, ha="center", fontsize=8.2)
    fig.suptitle("Figure 2 supporting evidence | All-five organization, response and post-withdrawal ambiguity", fontsize=12.4, y=.995)
    fig.text(.53, .055, "A: near-doubling pairs L19→39; Holmes 20→40; Spring L₀1→2. B–C: all 97 qualified basals; native units and clocks differ.\n"
        "D: core geography is threshold-defined, not a scalar replacement gate or tracked site identity. E: cyan marks withdrawal.\n"
        "F: dashed 500-unit cue returns to the incumbent; solid 1000-unit cue retains the challenger. Each has 13,000 cue-free units.",
        ha="center", fontsize=8.1, linespacing=1.5)
    save(fig, "figure2_measurement_context")


def supporting_native(profiles, displays):
    fig, axes = plt.subplots(2, 5, figsize=(13.4, 6.6))
    fig.subplots_adjust(left=.055, right=.985, top=.86, bottom=.20, wspace=.46, hspace=.7)
    for j, model in enumerate(MODELS):
        for size, style in zip(PAIRS[model], ("--", "-")):
            group = [r for r in profiles if r["model"] == model and float(r["nominal_size"]) == size]
            axes[0, j].plot([float(r["physical_position"]) for r in group], [float(r["active_concentration"]) for r in group],
                color=COLORS[model], ls=style, lw=1.5, label=f"Size {size:g}")
        axes[0, j].set(title=NAMES[model], xlabel="Physical native position", ylabel="Native active concentration")
        axes[0, j].legend(frameon=False, fontsize=7)
        d = next(d for d in displays.values() if d["row"]["model"] == model and d["row"]["primary"])
        T = d["row"]["cue_duration_native"]; keep = d["time"] <= 4*T
        im = axes[1, j].pcolormesh(d["fraction"], d["time"][keep]/T, d["field"][keep], shading="auto", cmap="magma", vmin=0, rasterized=True)
        axes[1, j].axhline(1, color="#26b9bd", lw=1.)
        axes[1, j].set(xlabel="Native position fraction", ylabel="Time / cue duration", ylim=(4, 0), title=f'Size {d["row"]["size_native"]:g}; cue {T:g}')
        cb = fig.colorbar(im, ax=axes[1, j], fraction=.05, pad=.025); cb.ax.tick_params(labelsize=6.5)
    fig.suptitle("Figure 2 supporting evidence | Native basal organization and selected complete-state histories", fontsize=12, y=.98)
    fig.text(.52, .055, "Top: exact near-doubling basal profiles, physical native coordinates. Bottom: prespecified largest qualified size in every model;\n"
        "zero-based, model-specific native activity scales; cyan = withdrawal. Spring position is material fraction. Full post-cue horizon is 13,000.\n"
        "The fixed-support corrected deBelly history is distinct from the prescribed relative-20% mechanism illustration in the retained source report.", ha="center", fontsize=8.2, linespacing=1.5)
    save(fig, "figure2_native_support")


def supporting_causal(controls, bolus, ablation):
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.9))
    fig.subplots_adjust(left=.1, right=.98, top=.86, bottom=.22, hspace=.62, wspace=.32)
    ax = axes[0, 0]
    labels = [("wave_pinning", 39.), ("holmes_model3", 28.), ("holmes_model3", 50.)]
    for i, (model, size) in enumerate(labels):
        for kind, offset, color, label in (("native", -.16, "#555555", "Native"), ("mobility_4x", .16, "#3c8b96", "4× inactive diffusion")):
            r = next(r for r in controls if r["model"] == model and r["size_native"] == size and r["intervention"] == kind)
            ax.bar(i+offset, r["final_orientation"], width=.29, color=color, label=label if i==0 else None)
        mesh = [r for r in controls if r["model"] == model and r["size_native"] == size and r["intervention"] == "mesh_401"]
        if mesh:
            ax.scatter(i-.16, mesh[0]["final_orientation"], marker="x", color="#d09328", s=45, zorder=5, label="401-cell control" if i==0 else None)
    ax.axhline(.05, color="#999999", ls="--", lw=.8)
    ax.set(xticks=range(3), xticklabels=["WP L39", "Holmes L28", "Holmes L50"], ylabel="Final native activity orientation", title="A  Mobility rescue is conditional", ylim=(-.17, .105))
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left")
    ax.text(.5, -.23, "Same 5120-unit fixed-support pulse; 10,000 off", transform=ax.transAxes, ha="center", fontsize=8)
    ax = axes[0, 1]
    geometries = ("fixed_physical", "fixed_relative", "fixed_slope")
    for i, geometry in enumerate(geometries):
        group = [r for r in controls if r["model"] == "spring" and r["geometry"] == geometry and r["intervention"] == "feedback_off_during_pulse"]
        for r in group:
            ax.scatter(r["size_native"], i, s=145, marker="s", facecolor="#3c8b96" if r["orientation_success"] else "#9b465a", edgecolor="white", linewidth=.7)
            ax.text(r["size_native"], i+.23, f'{r["cue_duration_native"]:.3g}', ha="center", fontsize=7)
    ax.set(yticks=range(3), yticklabels=["Fixed support 0.2", "Relative support 20%", "Constant slope"], ylim=(2.5, -.5),
        xlim=(.35, 2.85), xticks=(.6, 1.6, 2.2, 2.6), xlabel="Spring resting length", title="B  Mechanical feedback matters by cue geometry")
    ax.legend(handles=[Patch(facecolor="#3c8b96", label="Still passes"), Patch(facecolor="#9b465a", label="Loses native success")], frameon=False, loc="lower center", ncol=2, fontsize=7.6)
    ax.text(.5, -.23, "Only extension-dependent deactivation removed during pulse", transform=ax.transAxes, ha="center", fontsize=8)
    ax = axes[1, 0]
    for model in ("wave_pinning", "goryachev"):
        group = sorted([r for r in bolus if r["model"] == model], key=lambda r:float(r["length"]))
        ax.plot([float(r["length"]) for r in group], [float(r["reaction_gain_per_damage"]) for r in group], color=COLORS[model], marker="o", ms=3, lw=1.4, label=NAMES[model])
    ax.set(xlabel="Native interval length", ylabel="Immediate excess net reaction\n/ active material lost (native time⁻¹)", title="C  A separate substrate disturbance tests restoration")
    ax.legend(frameon=False, fontsize=8)
    ax.text(.5, -.25, "0.5% whole inactive-profile removal for 0.1; exact return", transform=ax.transAxes, ha="center", fontsize=8)
    ax = axes[1, 1]
    for condition, style, label in (("native", "-", "Native"), ("ablation", "--", "Tension→Rho increment removed")):
        group = sorted([r for r in ablation if r["condition"] == condition], key=lambda r:float(r["length_um"]))
        ax.plot([float(r["length_um"]) for r in group], [abs(float(r["left_boundary_displacement"])-float(r["right_boundary_displacement"])) for r in group],
            color="#2a9288", ls=style, marker="o", lw=1.5, label=label)
    ax.set(xlabel="Corrected deBelly interval length (µm)", ylabel="Settled boundary displacement (µm)", title="D  Resting mechanics has a causal tension contribution")
    ax.legend(frameon=False, fontsize=7.4)
    ax.text(.5, -.25, "No countercue; identical corrected starting states; qualified holdout", transform=ax.transAxes, ha="center", fontsize=8)
    fig.suptitle("Figure 2 supporting evidence | Selected interventions distinguish specific contributions", fontsize=12.2, y=.98)
    fig.text(.54, .045, "A–B are fixed-dose evolved controls, not new intervention-specific duration thresholds. B labels give each native successful pulse duration;\n"
        "native feedback returns for recovery. C is a different perturbation from the common active-loss assay. D concerns resting organization.\n"
        "These interventions do not assign additive causal fractions to polarity replacement or establish a shared physical clock between models.", ha="center", fontsize=8.2, linespacing=1.5)
    save(fig, "figure2_selected_interventions")


def discrimination(finite, sensitivity, controls):
    summary = []
    for subset_name, subset in (("primary", [r for r in finite if r["primary"]]), ("all_exposures", finite),
                                ("nonoverlapping_primary", [r for r in finite if r["primary"] and not r["cue_core_overlap"]])):
        for level in (.75, 1., 1.25):
            activated = [r for r in subset if r["cue_challenger_max_halfcontrast_ratio"] >= level]
            challenger_only = [r for r in subset if fate(r["post13000_old_halfcontrast_ratio"], r["post13000_challenger_halfcontrast_ratio"], level) == "challenger only"]
            true_positive = sum(r in challenger_only for r in activated)
            summary.append(dict(subset=subset_name, threshold_halfcontrast_ratio=level, exposures=len(subset),
                cue_activated=len(activated), challenger_only_final=len(challenger_only),
                activated_and_challenger_only=true_positive,
                activated_but_other_final=len(activated)-true_positive,
                challenger_only_without_cue_activation=len(challenger_only)-true_positive,
                positive_predictive_fraction=true_positive/len(activated),
                comment="Descriptive deterministic condition counts; sizes and doses are not independent biological replicates."))
    export("candidate_signature_discrimination.csv", summary)
    duration = []
    for model in MODELS:
        sizes = sorted({r["size_native"] for r in finite if r["model"] == model and r["anchor"]})
        for size in sizes:
            group = sorted([r for r in finite if r["model"] == model and r["size_native"] == size], key=lambda r:r["cue_duration_native"])
            assert len(group) == 3
            duration.append(dict(model=model, size_native=size, durations="/".join(f'{r["cue_duration_native"]:g}' for r in group),
                core_geographies=" / ".join(r["final_core_geography"] for r in group),
                fate_invariant_across_half_primary_double=len({r["final_core_geography"] for r in group}) == 1,
                orientation_successes="/".join(str(r["orientation_success"]) for r in group)))
    export("matched_state_duration_contrasts.csv", duration)
    mobility = [r for r in controls if r["intervention"] == "mobility_4x"]
    spring = [r for r in controls if r["intervention"] == "feedback_off_during_pulse"]
    CHECKS["discrimination_summary"] = summary
    CHECKS["duration_anchor_states"] = len(duration)
    CHECKS["duration_invariant_core_geographies"] = sum(r["fate_invariant_across_half_primary_double"] for r in duration)
    CHECKS["selected_mobility_orientation_gate_passes"] = sum(r["orientation_success"] for r in mobility)
    CHECKS["spring_feedback_loss_by_geometry"] = {g:sum(not r["orientation_success"] for r in spring if r["geometry"] == g) for g in ("fixed_physical", "fixed_relative", "fixed_slope")}
    return summary, duration


def documents(summary, normalized, finite, discrimination_rows, duration):
    matrix = [
        dict(explanation="Limited challenger initiation", passive="Opposite-region inactive concentration, or an inhibitory field; conditional on known input channel.",
            stimulation="Cue-weighted incremental activation and initial challenger growth under fixed physical support and peak.",
            active_loss="No unique prediction for incumbent deficit restoration.", intervention="No cross-model selective substrate-rescue experiment in this package; bolus removal is a separate incumbent disturbance.",
            observed_discrimination="Goryachev near-doubling gain 0.654; WP/Holmes/Spring near 1. Lower gain co-occurs with greater restoration and is not an independent causal partition.",
            limits="Concentration is not amount. Growing ramp peak or support changes initial production independently of this mechanism."),
        dict(explanation="Stronger local restoration of the incumbent", passive="Peak/width/amount distinguish organization but do not measure restoring response.",
            stimulation="Incumbent loss and rebound must be tracked separately from challenger activation.",
            active_loss="Sham-referenced small active deficit restoration at one fixed within-model clock; all native fields evolve.",
            intervention="Goryachev/WP whole-inactive-profile bolus with exact pointwise return probes a different disturbance; Goryachev has much larger size dependence in immediate reaction per damage.",
            observed_discrimination="Common restoration fold: Goryachev 3.44, Spring 1.61, WP/Holmes/deBelly about 1.01–1.10. All Spring primary exposures still replace, so larger restoration alone does not predict fate.",
            limits="Perturbation direction and clock matter. Positive deficit restoration is not generic stability, actual incumbent maintenance or an exhaustive cause of resistance."),
        dict(explanation="Transport and spatial redistribution limit replacement", passive="Broad incumbent with similar peak/gain is compatible with a redistribution burden; not uniquely diagnostic.",
            stimulation="Challenger activation with prolonged old-region retention or co-retention; native spatial histories distinguish this from weak initiation.",
            active_loss="Near-flat local restoration can coexist with strong size dependence of replacement; it does not exclude transport.",
            intervention="At WP L39, fourfold inactive diffusion throughout cue and recovery changes retained orientation from 0.04146 (fail) to 0.06215 (pass) after the same fixed-support 5120-unit pulse. Both fail panel A's stronger full-polarity-strength gate. Holmes L28/L50 remain outside the orientation gate at that dose.",
            observed_discrimination="Robust activation does not discriminate challenger-only retention from both-core retention. The WP control supports a transport contribution to retained orientation at a fixed dose; Holmes nonrescue does not establish absence of transport effects.",
            limits="Transport, total material and initial organization covary with size. No intervention-specific threshold curve or reduction of panel A's strength threshold has been measured; no purely transport-isolated explanation is claimed."),
        dict(explanation="Mechanical or slower-field organization changes response", passive="Extension, endpoint attachment and local flow can complement activity; activity alone may conceal these states.",
            stimulation="Track cue-on overshoot, post-withdrawal return and slower native fields; corrected deBelly39 changes fate between 500 and 1000 units from the same full state.",
            active_loss="Rac-only deBelly kick initially preserves MCA/mechanics; common active-loss restoration cannot quantify their rewriting difficulty.",
            intervention="Spring pulse-only extension-feedback removal loses all four fixed-slope successes but only one of three in each patch geometry. Corrected deBelly tension→Rho removal changes resting displacement spread.",
            observed_discrimination="Spring causal effect depends on stimulus geometry; corrected deBelly resting ablation is not a causal reversal-time result. A scalar high activity peak is insufficient.",
            limits="The Spring and deBelly mechanisms are different. Evolved intervention, passive temporal association and fixed-state budget are kept at their own inference levels.")]
    export("mechanism_observable_matrix.csv", matrix)
    headers = "| Explanation | Passive observation | Stimulation-dependent signature | Small active-loss assay | Evolved intervention / causal scope |\n|---|---|---|---|---|\n"
    matrix_md = headers+"\n".join("| "+" | ".join(r[k] for k in ("explanation", "passive", "stimulation", "active_loss", "intervention"))+" |" for r in matrix)
    write(OUT / "MECHANISM_OBSERVABLE_MATRIX.md", "# Mechanism × observable matrix\n\n"+matrix_md+"\n\n"+"\n\n".join("**"+r["explanation"]+".** "+r["observed_discrimination"]+" "+r["limits"] for r in matrix)+"\n")
    p = next(r for r in discrimination_rows if r["subset"] == "primary" and r["threshold_halfcontrast_ratio"] == 1.)
    a = next(r for r in discrimination_rows if r["subset"] == "all_exposures" and r["threshold_halfcontrast_ratio"] == 1.)
    folds_text = "\n".join(f'| {NAMES[r["model"]]} | {r["size_small"]:g}→{r["size_large"]:g} | {r["amount_fold"]:.4g} | {r["width_fold"]:.4g} | {r["peak_fold"]:.4g} | {r["gain_fold"]:.4g} | {r["repair_fold"]:.4g} |' for r in summary)
    caption = [
        "Similar size-dependent reversal behavior can arise through distinct mechanisms, which predict different observable responses to polarity challenges. Wave-pinning and Goryachev provide the lead comparison: both require longer fixed-support stimulation at increasing size over their observed scalar-strength successes, while their gain and local-restoration responses differ. Similarity is qualitative rather than equality of slopes or calibrated physical durations.",
        "A, Shortest tested fixed-support pulse duration passing the original scalar-strength replacement gate, divided by each model's own L19 duration (WP 73.5; Goryachev 377.5 native units). Both curves strictly increase over observed successes: WP L5–31, 5.906–4540; Goryachev L3–39, 33.726–3540. Open triangles retain WP L33–39 censoring at the 5120-unit tested cap; no time is assigned beyond that cap. Open circles retain scalar successes failing the separate native-field tail flag (WP L21–31). Preparations are unavailable at WP L1/L3 and Goryachev L1. Strength requires reversal of the signed scalar polarity with magnitude within 1% of its initial value, in addition to the original orientation, contrast and final-window persistence conditions; this does not prove complete field equilibration.",
        "B, Calculated cue-weighted incremental primary-activity production per native opposing input, for all 97 qualified basals. The fixed endpoint raised-cosine support is 3, or 0.2 in Spring. Gain and size divide by each model's own reference (L19, Holmes L20, Spring L₀1). Goryachev's near-doubling gain is 0.654 of the reference; lower available substrate supports weaker initial activation under the same local input. This calculation is conditional on the native input channel.",
        "C, Sham-relative deficit restoration after 0.5% of a smooth frozen-incumbent excess profile is deactivated; removed material enters an explicit inactive pool where present. Corrected deBelly changes Rac alone and initially preserves MCA/mechanics. All native fields evolve. The repaired fraction at one fixed native clock per model is divided by its own reference. Clocks are WP/Holmes 0.177827941, Goryachev 0.005623413, deBelly 5 and Spring 0.017782794; 0.1%/1% amplitude and numerical checks qualify all 97 cases. This measures response to a declared deficit, not general stability.",
        "D–E, Exact native WP/Goryachev L39 fields during a 1000-unit fixed-support cue and the next 3000 units; cyan marks withdrawal. WP develops a robust challenger and retains two threshold-defined regions through the complete 13,000-unit off-window. Goryachev develops a smaller but nonzero challenger, below the original half-contrast, then restores its incumbent. Native concentration scales are model-specific and start at zero; these are different within-model native doses, not calibrated equal cross-model inputs.",
        f"F, Native histories reanalysed for every primary exposure (97 filled points) and every prespecified half/double-duration anchor (30 open points). Horizontal coordinate is the largest challenger-core peak during the cue, above frozen basal background and divided by initial half-contrast. Vertical coordinate is challenger-core excess after 13,000 cue-free units divided by initial incumbent excess. Colour gives final original/mirrored-core geography: challenger only, both cores or incumbent only, using the same half-contrast threshold; it is not a tracked site identity or scalar replacement gate. Model symbols are WP circle, Goryachev square, deBelly triangle-up, Holmes diamond and Spring triangle-down. Among {p['cue_activated']} primary exposures reaching the cue-on threshold, {p['activated_but_other_final']} retain both cores; a strong challenger is therefore insufficient to establish replacement. Across all 127 exposures, {a['activated_but_other_final']} threshold-positive cues have another final geography, and {a['challenger_only_without_cue_activation']} challenger-only outcomes cross the threshold only after withdrawal. This quantifies observational ambiguity rather than identifying a unique molecular mechanism.",
        "Cue peaks are WP 0.1, Goryachev 6.4, corrected deBelly 3.6, Holmes 1 and Spring 4; primary duration is 1000 native units, or 50 in Spring. Native times, channels and sizes are not calibrated across models. Thirteen preparation-ineligible sizes remain in the source table. The supporting WP mobility control passes retained orientation at a fixed 5120-unit pulse, while both native and fourfold-diffusion conditions fail A's stronger full-polarity-strength gate; fourfold diffusion acts throughout cue and recovery. Selected interventions and complete native profiles/histories are supporting evidence; no new integrations were run."
    ]
    write(OUT / "CAPTION.md", "# Figure 2. Similar reversal trends, distinct challenge responses\n\n"+"\n\n".join(caption)+"\n")
    evidence = f'''# Figure 2 evidence, methods and limits

The user's central claim is: **Similar size-dependent reversal behavior can arise through distinct mechanisms, which predict different observable responses to polarity challenges.** The figure begins with two increasing fixed-support duration curves, then tests mechanism-linked responses. Wave-pinning and Goryachev share the qualitative size trend over successful scalar-strength cases, not a common exponent, slope or calibrated physical clock. WP L5–31 increases 5.906→4540 and is censored at L33–39; Goryachev L3–39 increases 33.726→3540. WP L21–31 scalar successes fail the separate native-field tail flag and remain visible. Corrected deBelly is nonmonotone; Spring changes trend with cue geometry; the claim is not extended to every model and stimulus. All five × three cue classes are preserved in [the full threshold context](supporting/data/all_five_three_cue_strength_thresholds.csv).

The lead mechanisms differ in measured challenge responses. Goryachev combines lower available substrate and stronger local restoration; WP retains nearly constant initial gain and recovery near its reference while its larger incumbent broadens and redistribution/competition contributes. The positive causal evidence is conditional: fourfold inactive mobility changes WP's retained orientation from failing to passing at a fixed dose, while both conditions fail panel A's stronger full-polarity-strength gate; common deficit restoration directly measures one restoring response; other evolved controls show cue-dependent mechanical contributions. The ambiguity analysis in F is a negative discrimination result, not a molecular mechanism classifier.

## Exact size contrasts

| Model | Pair | Excess amount | Width | Peak | Initial gain | Common repaired fraction |
|---|---|---:|---:|---:|---:|---:|
{folds_text}

## Candidate signatures across size and stimulus conditions

All 97 primary fixed-support trajectories and all 30 prespecified duration controls are read from their original hash-checked native fields. For a fixed original/mirrored core, let the basal background be the minimum initial concentration and let half-contrast be half its initial range. A peak ratio of one is the initial half-height threshold. Native fields independently reproduce every archived peak and regional challenger amount, including Spring's evolving physical length. Final core geography uses the two peaks after 13,000 units; full connected-region counts and scalar orientation/field gates remain separate columns.

For the separate small active-loss assay, let the original-core amount deficit between paired sham and probe be $D(t)$. The plotted quantity is the restored fraction at a fixed model-native clock, normalized by its own reference value:

\\[
D(t)=A_{{\\mathrm{{sham}}}}(t)-A_{{\\mathrm{{probe}}}}(t),
\\qquad F_{{\\tau}}=\\frac{{D(0)-D(\\tau)}}{{D(0)}},
\\qquad K_{{\\tau}}=\\frac{{F_{{\\tau}}}}{{\\tau}}.
\\tag{{1}}
\\]

The fixed within-model clock makes the size folds of $F_{{\\tau}}$ and $K_{{\\tau}}$ identical; clocks and rates are not ranked across models.

At that threshold, {p['cue_activated']}/{p['exposures']} primary exposures reach a strong cue-on challenger; {p['activated_but_other_final']} of those retain both original and challenger cores. One primary Goryachev condition ends challenger-only although its cue-on challenger remains below threshold, exposing delayed growth after withdrawal. Including the duration controls yields {a['cue_activated']} activated exposures, of which {a['activated_but_other_final']} do not end challenger-only: co-retention and corrected deBelly39's short-pulse incumbent return are distinct. These are deterministic condition counts; correlated sizes and doses are not biological replicates and no sampling confidence interval is implied.

The [discrimination table](supporting/data/candidate_signature_discrimination.csv) repeats this test at 0.75, 1 and 1.25 times the initial half-contrast and separately removes cue/core-overlap cases. The [matched-state duration table](supporting/data/matched_state_duration_contrasts.csv) retains all 15 short/middle/long anchors. {sum(r['fate_invariant_across_half_primary_double'] for r in duration)}/15 core geographies are unchanged over half/primary/double exposure; corrected deBelly39 changes from incumbent-only at 500 to challenger-only at 1000/2000. Goryachev19 can grow above threshold only after withdrawal at 500 units. Continuous small cue responses are never relabelled as no activation.

The gain and recovery series are not independent causal axes: in Goryachev both change with size and the same conserved pool. Positive deficit restoration also occurs in Spring, whose primary long exposures all end challenger-only. Thus organization alone, gain alone, positive restoration alone and a strong challenger alone each leave observational equivalences. The [mechanism × observable matrix](MECHANISM_OBSERVABLE_MATRIX.md) records positive support, failures of discrimination and the intervention needed to test the remaining ambiguity.

## Selected causal evidence, preserved separately

The [supporting intervention figure](figure2_selected_interventions.png) reads retained evolved controls, without new threshold searches. Scalar orientation and strength gates are independently recomputed from the final 500 units of the 10,000-unit recovery arrays. At Wave-pinning L39, the same 5120-unit fixed-support pulse yields retained orientation 0.04146 natively (orientation gate fails) and 0.06215 with fourfold inactive diffusion (orientation gate passes). The intervention remains active throughout both cue and recovery. Both conditions fail panel A's stronger full-polarity-strength gate: the initial scalar magnitude is 0.43234 and its 1% acceptance band is 0.42802–0.43666. Both pass the separate native-field tail flag. Thus the control supports a transport contribution to retained orientation at this fixed dose; it does not establish a reduction in the strength-duration threshold plotted in A. A 401-cell control remains outside the orientation gate. Holmes L28/L50 also remain outside that gate under fourfold inactive Rac/Rho diffusion throughout cue and recovery at the same tested duration. Their nonrescue does not establish absence of transport contributions or the impossibility of a longer successful pulse.

Spring removes only extension-dependent global deactivation during each original successful native pulse, restoring native feedback during recovery. All four fixed-slope successes are lost, whereas only one of three is lost under either patch geometry. The fixed-slope input peak normally grows with length; matching its peak to the reference separately defeats three of four responses. This is a cue-dependent causal effect, not a general mechanical explanation inferred from model identity.

The separate whole-inactive-profile bolus removes 0.5% for 0.1 native units then returns the exact original profile pointwise. Its immediate reaction per core active deficit increases from 0.97309 to 15.28539 across Goryachev L9→39, versus 0.61734→0.79019 in WP. This global substrate disturbance and the main local active-loss assay have distinct amplitudes, clocks, directions and recovery tails; their coefficients are never merged. Native kinetics remain unchanged, so this tests restoration after disturbance rather than an isolated autocatalytic eigenvalue.

Corrected deBelly's tension-to-Rho resting ablation starts from identical corrected checkpoints and has qualified settling plus a 10,000-second holdout. The displacement spread across L23/L31/L39 decreases from 1.968526 to 0.976774 µm (50.38%), while Rac contrast remains 0.205–0.210. L31/L39 normalized orientations fall below 0.05. This is a causal contribution to resting mechanical organization, not a measured contribution to reversal duration. No deprecated centered-endpoint deBelly control from the original three-test mechanism campaign enters this analysis.

## Reproduction and inference limits

Run `/home/4473331/.conda/envs/scaling_sandbox/bin/python polarity1d/manuscript/code/mechanisms/revision_figure2_20261019.py` from the repository root. This is a single-core render/reanalysis of existing deterministic data, without integrations or campaign changes. All source and artifact hashes appear in [manifest.json](supporting/metadata/manifest.json); [verification.json](supporting/metadata/verification.json) records native recomputation and existing completed scientific audits. The source tables preserve preparation unavailability, cue/core overlap, full native field flags and discrete region counts.

Core geography concerns fixed spatial regions, not site ancestry, unique protrusions or a migration law. Native Goryachev L39 has four sample centres in its original core; existing solver checks do not establish spatial convergence. The original fixed-support cues, common recovery disturbance and selected intervention doses are not experimental calibrations across models. The faint but positive Goryachev challenger is a measured response, even when below a half-contrast gate. Threshold sensitivity and full size curves prevent a chosen cutoff or nearly doubled pair from replacing the full result. Figure 3 owns ring placement/geometry; no ring-derived number is substituted here.
'''
    write(OUT / "EVIDENCE_METHODS_LIMITS.md", evidence)
    title = "Figure 2. Similar reversal trends, distinct challenge responses"
    base = str(OUT.relative_to(REPO))
    content = dict(title=title, role="Similar size-dependent reversal behavior can arise through distinct mechanisms, which predict different observable responses to polarity challenges.",
        assets=[dict(path=base+"/figure2_mechanisms.png", caption=title, page_size=[11.8, 9.4], page_title=title, minimum_printed_width_inches=10)],
        legend=caption,
        narrative=["Wave-pinning and Goryachev both show increasing fixed-support scalar-strength duration over observed successes, but their magnitude, eligibility and censoring differ. This is qualitative similarity, not a shared scaling law.",
            "Incumbent material increases in every model, but concentration, width, initial responsiveness and local restoration change differently.",
            "Substrate availability and a declared local deficit assay jointly support a coupled activation/restoration explanation in Goryachev; their correlation does not partition reversal time.",
            f"A strong challenger is insufficient: {p['activated_but_other_final']} of {p['cue_activated']} activated primary exposures retain both original and challenger cores. Corrected deBelly39 exposes additional incumbent return after a shorter exposure.",
            "Transport and mechanical interventions show causal contributions under specified cues. The WP fixed-dose mobility control passes retained orientation while both conditions fail A's stronger strength gate; fourfold inactive diffusion acts throughout cue and recovery. Cue geometry changes intervention effects, so no single activity signature uniquely identifies an explanation."],
        measurements=["Measure incumbent peak, physical extent and integrated excess signal before stimulation, keeping concentration-like intensity distinct from calibrated amount.",
            "Measure initial challenger activation under characterized fixed local input, and sham-relative incumbent deficit restoration under a specified small suppression.",
            "Track original and challenger regions throughout stimulation and withdrawal; pair their local activity with substrate, attachment, flow or extension measurements when testing the corresponding mechanism."],
        limitations=["Native clocks and channels differ; the assays are not cross-model biological calibrations.",
            "Positive restoration, strong activation and final fixed-core geography each leave observational equivalences; no unique molecular model, general stability or velocity law is inferred.",
            "Deterministic sizes and doses are conditions rather than independent experimental replicates. Finite-window patterns do not prove attractors."],
        sources=[dict(label="Evidence, methods and limits", path=base+"/EVIDENCE_METHODS_LIMITS.md"),
            dict(label="Mechanism × observable matrix", path=base+"/MECHANISM_OBSERVABLE_MATRIX.md"),
            dict(label="All 127 native activation/fate observations", path=base+"/supporting/data/all127_native_activation_and_fate.csv"),
            dict(label="Threshold and dose discrimination", path=base+"/supporting/data/candidate_signature_discrimination.csv")],
        supplementary=[dict(title="All-five response context and post-withdrawal ambiguity", assets=[dict(path=base+"/figure2_measurement_context.png")],
            legend=["All-five organization/gain/recovery, all 127 native cue/fate observations and exact corrected deBelly39 duration contrast: strong challenger activation followed by incumbent return after 500 units, versus challenger retention after 1000, with 13,000 cue-free units in both. Symlog time is linear through 1000 then logarithmic."], narrative=[], measurements=[], limitations=[], sources=[]),
            dict(title="Selected evolved interventions", assets=[dict(path=base+"/figure2_selected_interventions.png")],
            legend=["WP L39 retained-orientation control: fourfold inactive diffusion throughout cue and recovery changes final orientation from 0.04146 (fail) to 0.06215 (pass) after the same 5120-unit pulse. Both conditions fail main panel A's stronger full-polarity-strength gate. The sheet also retains cue-specific Spring feedback removal, separate Goryachev/WP inactive-bolus restoration and corrected deBelly resting ablation. Exact protocols and limits are in the linked evidence file."],
            narrative=[], measurements=[], limitations=["The WP orientation-gate result does not establish a change in panel A's strength threshold. Fixed-dose controls do not locate new duration thresholds; resting ablation is separate from reversal causation."], sources=[]),
            dict(title="Native profiles and illustrative histories", assets=[dict(path=base+"/figure2_native_support.png")], legend=["Exact basal profiles at the stated near-doubling pairs and prespecified largest qualified size in each model."], narrative=[], measurements=[], limitations=[], sources=[])])
    write(META / "content.json", json.dumps(content, indent=2)+"\n")
    write(OUT / "README.md", "# Revision-program Figure 2\n\n**Similar size-dependent reversal behavior can arise through distinct mechanisms, which predict different observable responses to polarity challenges.**\n\n[Main figure](figure2_mechanisms.png) · [PDF](supporting/vectors/figure2_mechanisms.pdf) · [SVG](supporting/vectors/figure2_mechanisms.svg)\n\n[Caption](CAPTION.md) · [Evidence, methods and limits](EVIDENCE_METHODS_LIMITS.md) · [Mechanism × observable matrix](MECHANISM_OBSERVABLE_MATRIX.md)\n\n[All-five measurement context](figure2_measurement_context.png) · [Selected interventions](figure2_selected_interventions.png) · [Native profiles and histories](figure2_native_support.png)\n\nSource tables and numerical audits are under supporting/. All producers are finished; this package creates no simulations.\n")


def main():
    for folder in (OUT, DATA, META, VECTORS, OUT / "supporting/viewers"):
        folder.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":8.4, "axes.titlesize":9., "axes.titlepad":11,
        "axes.spines.top":False, "axes.spines.right":False, "legend.fontsize":8., "pdf.fonttype":42})
    source(Path(__file__))
    folds, profiles, normalized, pair_summary = function_data()
    finite, sensitivity, displays = finite_data()
    controls, bolus, ablation = control_data()
    discrimination_rows, duration = discrimination(finite, sensitivity, controls)
    endpoints = reversal_data()
    claim_figure(endpoints, normalized, finite, displays)
    measurement_context(folds, normalized, finite, displays)
    supporting_native(profiles, displays)
    supporting_causal(controls, bolus, ablation)
    documents(pair_summary, normalized, finite, discrimination_rows, duration)
    verification = dict(generated_at=datetime.now(timezone.utc).isoformat(), passed=True,
        no_new_integrations=True, source_producer_state="all retained campaigns finished", checks=CHECKS)
    write(META / "verification.json", json.dumps(verification, indent=2)+"\n")
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), artifact_root="../..", source_sha256=SOURCES,
        sha256={str(p.relative_to(OUT)):sha(p) for p in ARTIFACTS}, new_integrations=0)
    (META / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(dict(output=str(OUT), native_finite_branches=len(finite), basal_gains=len(normalized), evolved_controls=len(controls),
        checks=str(META / "verification.json"), exports=len(ARTIFACTS)), indent=2))


if __name__ == "__main__":
    main()
