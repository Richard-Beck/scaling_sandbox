"""Render-only assembly for the three-figure vision-aligned checkpoint.

Native campaigns are owned separately. This module never integrates a model or
substitutes an older assay for a missing measurement. Final rendering is added
only after the completed input contracts are supplied and checked.
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
from matplotlib.colors import Normalize, ListedColormap
from matplotlib.patches import Patch
import length_dev_output_layout as layout
import expanded_ring_components_20261017 as ring_components

ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT.parents[1]
DEV = ROOT / "figures/length_mechanisms_20261008/dev"
OUT = DEV / "vision_aligned_checkpoint_20261017"
MODELS = ("wave_pinning", "goryachev", "debelly", "holmes_model3", "spring")
RING_MODELS = ("wave_pinning", "goryachev", "holmes_model3")
NAMES = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "deBelly", "Holmes 3", "Spring")))
COLORS = dict(zip(MODELS, ("#3976a8", "#c96926", "#2a9288", "#8560a4", "#a34f7f")))
REFS = {"wave_pinning":19., "goryachev":19., "debelly":19., "holmes_model3":20., "spring":1.}


class Evidence:
    """Record exact source/output identities without modifying any source."""

    def __init__(self):
        self.sources = {}
        self.artifacts = []

    @staticmethod
    def sha(path):
        digest = hashlib.sha256()
        with Path(path).open("rb") as stream:
            for block in iter(lambda: stream.read(2**20), b""):
                digest.update(block)
        return digest.hexdigest()

    def source(self, path, expected=None):
        path = Path(path).resolve()
        digest = self.sha(path)
        if expected is not None:
            assert digest == expected, path
        self.sources[str(path.relative_to(REPO))] = digest
        return path

    def rows(self, path, required=()):
        with self.source(path).open(newline="") as stream:
            reader = csv.DictReader(stream)
            assert set(required).issubset(reader.fieldnames or ()), (path, required)
            return list(reader)

    def filtered_rows(self, path, predicate, required=()):
        """Stream large ensemble exports; retained frames are exact samples."""
        with self.source(path).open(newline="") as stream:
            reader = csv.DictReader(stream)
            assert set(required).issubset(reader.fieldnames or ()), (path, required)
            return [row for row in reader if predicate(row)]

    def json(self, path):
        return json.loads(self.source(path).read_text())

    def npz(self, path, expected=None):
        with np.load(self.source(path, expected)) as archive:
            return {key: archive[key].copy() for key in archive.files}

    def save(self, fig, stem):
        assert stem.startswith("vision_aligned_")
        for extension in (".png", ".pdf", ".svg"):
            path = layout.dev_output_path(OUT, stem + extension)
            fig.savefig(path, dpi=200 if extension == ".png" else None, bbox_inches="tight")
            self.artifacts.append(path)
        plt.close(fig)

    def export(self, name, data):
        path = layout.dev_output_path(OUT, name)
        keys = list(dict.fromkeys(key for row in data for key in row))
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys)
            writer.writeheader()
            writer.writerows(data)
        self.artifacts.append(path)
        return path

    def export_npz(self, name, data):
        path = layout.dev_output_path(OUT,name)
        np.savez_compressed(path,**data)
        self.artifacts.append(path)
        return path

    def manifest(self, checks):
        path = layout.dev_output_path(OUT, "vision_aligned_checkpoint.manifest.json")
        record = dict(generated_at=datetime.now(timezone.utc).isoformat(),
                      artifact_root="../..", source_sha256=self.sources,
                      sha256={str(p.relative_to(OUT)):self.sha(p) for p in self.artifacts},
                      checks=checks)
        path.write_text(json.dumps(record, indent=2) + "\n")
        return path


def linear_organization(evidence):
    topic = DEV / "integrated_framework_candidate_20261017/supporting/data"
    folds = evidence.rows(topic / "all_five_incumbent_size_folds.csv",
                          ("model", "excess_active_fold", "width_fold", "peak_concentration_fold"))
    profiles = evidence.rows(topic / "all_five_native_basal_profiles.csv",
                             ("model", "nominal_size", "physical_position", "active_concentration"))
    assert len(folds) == 5 and {r["model"] for r in folds} == set(MODELS)
    return folds, profiles


def linear_gain(evidence):
    path = DEV / "all_five_linear_function_20261017/supporting/data/native_gain_design_evidence.csv"
    raw = evidence.rows(path, ("model", "size_native", "physical_length", "cue_weighted_gain_native"))
    assert len(raw) == 97
    result = []
    for row in raw:
        ref, = [r for r in raw if r["model"] == row["model"]
                and float(r["size_native"]) == REFS[row["model"]]]
        result.append(dict(**row,
                           reference_nominal_size=REFS[row["model"]],
                           reference_gain_native=float(ref["cue_weighted_gain_native"]),
                           relative_nominal_size=float(row["size_native"]) / REFS[row["model"]],
                           gain_over_own_reference=float(row["cue_weighted_gain_native"])
                               / float(ref["cue_weighted_gain_native"])))
    return result


def recovery_rows(evidence, path, expected_cases):
    """Adapt the native assay collector's clock-decision structure.

    Gate failures remain explicit rows with their measured values/status. The
    renderer never upgrades a failed measurement to a qualified coefficient.
    """
    record = evidence.json(path)
    assert record["complete"], path
    rows = []
    for block in record["models"]:
        tau = float(block["selected_horizon_native"])
        for result in block["primary_results"]:
            F = float(result["primary_F"])
            K = float(result["primary_K"])
            np.testing.assert_allclose(F, tau * K, rtol=3e-14, atol=3e-17)
            rows.append(dict(model=result["model"],
                             size_native=float(result["size_native"]),
                             horizon_native=tau, restored_fraction=F,
                             conditional_K_per_native_time=K,
                             empirical_E=float(result["primary_E"]),
                             qualified=bool(result["passes"])))
    assert len(rows) == expected_cases
    assert len({(r["model"],r["size_native"]) for r in rows}) == len(rows)
    return record, rows


def common_linear_recovery(evidence):
    path = ROOT / "data/runs/vision_linear_20261017/recovery/results/recovery_summary.csv"
    status = evidence.json(path.parent / "status.json")
    audit = evidence.json(path.parent / "verification.json")
    assert status["complete"] and status["all_fixed_clocks_qualified"]
    assert audit["passed"] and audit["completed_cases"] == 97
    assert audit["frozen_inputs_and_sources_verified"] and audit["zero_kick_all_fields_identity"]
    raw = evidence.rows(path, ("model","size_native","horizon_native",
                               "restored_fraction","conditional_K_per_native_time",
                               "coefficient_available"))
    assert len(raw) == 97
    rows = []
    for row in raw:
        result = dict(row)
        result["qualified"] = row["coefficient_available"] == "True"
        result["restored_fraction"] = float(row["restored_fraction"])
        result["size_native"] = float(row["size_native"])
        np.testing.assert_allclose(result["restored_fraction"],
                                   float(row["horizon_native"])
                                   * float(row["conditional_K_per_native_time"]),
                                   rtol=3e-14,atol=3e-17)
        rows.append(result)
    assert len({(r["model"],r["size_native"]) for r in rows}) == 97
    return rows


def linear_finite(evidence):
    root = ROOT / "data/runs/vision_linear_20261017/finite"
    rows = evidence.rows(root / "results/metrics.csv")
    checks = evidence.json(root / "results/verification.json")
    assert checks["passed"] and checks["complete_primary"] == 97
    assert checks["complete_sensitivity"] == 30 and checks["ineligible_visible"] == 13
    primary = [r for r in rows if r["primary"] == "True" and r["status"] == "complete"]
    assert len(primary) == 97
    display,selected = {},[]
    fields = {"wave_pinning":"active","goryachev":"active","debelly":"rac",
              "holmes_model3":"rac_active","spring":"active"}
    for model in MODELS:
        row = max([r for r in primary if r["model"] == model],
                  key=lambda r:float(r["size_native"]))
        archive = evidence.npz(root / row["archive"],row["archive_sha256"])
        field = fields[model];T = float(row["cue_duration_native"])
        cue,post = archive["cue_field_"+field],archive["post_field_"+field]
        np.testing.assert_array_equal(cue[-1],post[0])
        t = np.r_[archive["cue_time"],T+archive["post_time"][1:]]
        u = np.r_[cue,post[1:]]
        assert t[0] == 0 and t[-1] == T + 13000 and np.all(np.diff(t)>0)
        keep = t <= 4*T
        display[model] = dict(time_native=t[keep],time_over_duration=t[keep]/T,
                              fractional_position=archive["x"]/float(row["size_native"]),
                              primary=u[keep],full_history_min=float(u.min()),
                              full_history_max=float(u.max()),cue_duration_native=T,
                              size_native=float(row["size_native"]))
        if model == "spring":
            length = np.r_[archive["cue_field_length"],archive["post_field_length"][1:]]
            display[model]["physical_length"] = length[keep]
        selected.append(dict(**row,source_archive_path=str(root / row["archive"]),
                             selection_rule="largest eligible native size in every model",
                             display_start_native=0,display_end_native=4*T,
                             coordinate="initial material fraction" if model=="spring"
                             else "fractional nominal native interval",
                             field=field))
    return rows,primary,display,selected


def matched_native_profiles(evidence):
    """Actual C30 and L15 profiles, with interval reflection for display only.

    Cell-center positions are retained. Ring positions are cyclically permuted
    to [-L,L); interval samples are reflected about their geometric boundary.
    This is a plotting operation, never a replacement/preparation of a state.
    """
    path = DEV / "expanded_ring_components_20261017/supporting/data/expanded_ring_baseline_pair_measurements.csv"
    pairs = evidence.rows(path, ("model", "linear_L", "ring_C", "linear_checkpoint", "ring_checkpoint"))
    assert len(pairs) == 13
    for pair in pairs:
        for geometry in ("linear","ring"):
            evidence.source(pair[geometry+"_checkpoint"],pair[geometry+"_checkpoint_sha256"])
    profiles = []
    for model in RING_MODELS:
        pair, = [r for r in pairs if r["model"] == model and float(r["linear_L"]) == 15.]
        linear = evidence.npz(pair["linear_checkpoint"], pair["linear_checkpoint_sha256"])
        ring = evidence.npz(pair["ring_checkpoint"], pair["ring_checkpoint_sha256"])
        L, C = float(pair["linear_L"]), float(pair["ring_C"])
        assert C == 2 * L and np.argmax(linear["state"][0]) == 0
        ring_x = (ring["x"] + C / 2) % C - C / 2
        order = np.argsort(ring_x)
        for x, a in zip(ring_x[order], ring["state"][0, order]):
            profiles.append(dict(model=model, geometry="native ring", position=float(x),
                                 active_concentration=float(a), linear_L=L, ring_C=C))
        for sign in (-1, 1):
            for x, a in zip(sign * linear["x"], linear["state"][0]):
                profiles.append(dict(model=model, geometry="interval reflected for display",
                                     position=float(x), active_concentration=float(a), linear_L=L, ring_C=C))
    return pairs, profiles


def retained_ring_data(evidence):
    """Read exact retained finite-challenge context; no new assay coefficients."""
    topic = DEV / "expanded_ring_components_20261017/supporting/data"
    gains = evidence.rows(topic / "expanded_ring_initial_gain_map.csv")
    outcomes = evidence.rows(topic / "expanded_ring_outcome_map.csv")
    display = evidence.npz(topic / "expanded_ring_Holmes_native_display_samples.npz")
    vectors = evidence.rows(DEV / "global_polarity_vector_20261017/supporting/data/global_polarity_traces.csv")
    evidence.source(Path(ring_components.__file__))
    assert len(gains) == len(outcomes) == 867
    return gains, outcomes, display, vectors


def common_ring_recovery(evidence):
    root = ROOT / "data/runs/vision_ring_recovery_20261017/results"
    raw = evidence.rows(root / "recovery.csv")
    primary = [r for r in raw if float(r["epsilon"]) == .005]
    basals = {(r["model"],float(r["circumference"])):r for r in primary if r["stage"] == "basal"}
    assert len(basals) == 50 and len(primary) == 702
    matrix = evidence.rows(root / "challenged_full_matrix.csv")
    assert len(matrix) == 867
    result = []
    for row in matrix:
        row = dict(row)
        if row["status"] == "complete":
            assert float(row["epsilon"]) == .005 and row["passes"] == "True"
            basal = basals[(row["model"],float(row["circumference"]))]
            assert row["tau"] == basal["tau"]
            F, F0, tau = float(row["F"]),float(basal["F"]),float(row["tau"])
            np.testing.assert_allclose(F, tau * float(row["K"]), rtol=3e-14)
            row.update(basal_F=F0, basal_K=float(basal["K"]),
                       recovery_over_basal=F/F0)
            assert 0 < F/F0 < 1 and float(row["sham_core_change"]) < 0
        result.append(row)
    pairs = evidence.rows(root / "ring_linear_common_recovery_pairs.csv")
    assert len(pairs) == 13
    for pair in pairs:
        assert float(pair["ring_circumference"]) == 2 * float(pair["linear_length"])
        assert pair["linear_coefficient_available"] == pair["ring_prefix_gates_pass"] == "True"
        evidence.source(pair["linear_source_csv"],pair["linear_source_csv_sha256"])
    return list(basals.values()),result,pairs


def noise_tables(evidence, folder, prefix="vision_noise_"):
    """Stream only plotted quantities, retaining their native 5-unit cadence."""
    folder = Path(folder)
    audit = evidence.json(folder.parent/"metadata"/(prefix+"analysis_audit.json"))
    assert audit["cases"]==243 and audit["groups"]==45
    assert audit["maxlag_native_time"]==8000 and audit["bootstrap_draws"]==1000
    assert not audit["exponential_fit"] and not audit["burn_in_selected"]
    evidence.source(audit["generator_path"],audit["generator_sha256"])
    seeds = evidence.rows(folder / (prefix+"seed_summary.csv"))
    noise = [r for r in seeds if r["condition"] == "noise"]
    assert len(seeds) == 243 and len(noise) == 216
    counts = {}
    for row in noise:
        key = (row["model"],float(row["circumference"]),float(row["fraction"]))
        counts[key] = counts.get(key,0) + 1
    assert len(counts) == 18 and set(counts.values()) == {12}
    assert set(counts) == {(m,C,f) for m in RING_MODELS
                          for C in (20.,40.,60.) for f in (.001,.01)}
    assert sum(r["condition"]=="mean" for r in seeds) == 18
    assert sum(r["condition"]=="zero" for r in seeds) == 9
    traces = evidence.filtered_rows(
        folder / (prefix+"ensemble_traces.csv"),
        lambda r:r["metric"]=="R" or
            (r["metric"]=="initial_memory" and r["condition"]=="noise"),
        ("model","circumference","condition","fraction","time","metric",
         "n_seeds","mean","ci_lower","ci_upper"))
    assert max(float(r["time"]) for r in traces) == 32000
    lags = evidence.filtered_rows(folder / (prefix+"lag_correlations.csv"),
                                  lambda r:r["window"]=="full")
    persistence = evidence.filtered_rows(folder / (prefix+"persistence.csv"),
                                         lambda r:r["window"]=="full")
    physical = evidence.filtered_rows(
        folder / (prefix+"site_displacement.csv"),
        lambda r:r["window"]=="full" and float(r["lag"])==8000)
    selected = evidence.rows(folder / (prefix+"representative_histories.csv"))
    prescribed = evidence.rows(ROOT / "data/runs/vision_noise_main_20261017/results/statistics_workflow/all_prescribed_cases.csv")
    assert len(prescribed) == 243
    evidence.source(audit["index_path"],audit["index_sha256"])
    for name in ("seed_summary.csv","ensemble_traces.csv","lag_correlations.csv",
                 "persistence.csv","site_displacement.csv","representative_histories.csv"):
        path = folder/(prefix+name)
        digest = evidence.sources[str(path.resolve().relative_to(REPO))]
        assert digest == audit["outputs"][str(path.relative_to(folder.parent.parent))]
    assert max(float(r["restricted_horizon"]) for r in persistence) == 8000
    return dict(seeds=seeds,traces=traces,lags=lags,persistence=persistence,
                physical=physical,selected=selected,prescribed=prescribed,
                display_time_step_native=5)


def noise_series(rows, model, C, f, condition, metric=None):
    block = [r for r in rows if r["model"]==model
             and float(r["circumference"])==C
             and float(r["fraction"])==f and r["condition"]==condition
             and (metric is None or r.get("metric")==metric)]
    time_key = "time" if metric is not None else "lag"
    return sorted(block,key=lambda r:float(r[time_key]))


def plot_noise_series(ax, rows, xkey, ykey, color, style, ci=True, lw=1.1):
    assert rows
    x = np.array([float(r[xkey])/1000 for r in rows])
    y = np.array([float(r[ykey]) for r in rows])
    ax.plot(x,y,color=color,ls=style,lw=lw)
    if ci:
        lo = np.array([float(r["ci_lower"]) for r in rows])
        hi = np.array([float(r["ci_upper"]) for r in rows])
        keep = np.isfinite(lo) & np.isfinite(hi)
        ax.fill_between(x,lo,hi,where=keep,color=color,alpha=.09,lw=0)


def noise_history_examples(evidence, tables):
    choices = (("wave_pinning","noise"),("goryachev","noise"),
               ("goryachev","mean"),("holmes_model3","noise"))
    histories,selected = [],[]
    for model,condition in choices:
        if condition=="noise":
            row, = [r for r in tables["selected"] if r["model"]==model
                    and float(r["circumference"])==20 and float(r["fraction"])==.01
                    and r["role"]=="within_condition_descriptor_medoid"]
            path,expected = row["profile_archive"],row["profile_archive_sha256"]
            case_id = row["case_id"]
            selection = "condition-descriptor medoid among 12 declared realizations"
        else:
            row, = [r for r in tables["prescribed"] if r["model"]==model
                    and float(r["circumference"])==20 and float(r["fraction"])==.01
                    and r["condition"]=="mean"]
            path,expected = row["profile_archive"],row["profile_sha256"]
            case_id = row.get("case_id",Path(path).parent.name)
            selection = "deterministic equal-mean control paired to the Goryachev condition"
        raw = evidence.npz(path,expected)
        assert float(raw["circumference"])==20 and raw["time"][0]==0 and raw["time"][-1]==32000
        origin = raw["x"][int(np.argmax(raw["primary"][0]))]
        x = (raw["x"]-origin+5)%20-5
        order = np.argsort(x)
        histories.append(dict(model=model,condition=condition,x=x[order],
                              time=raw["time"],primary=raw["primary"][:,order]))
        selected.append(dict(model=model,circumference=20,fraction=.01,
                             condition=condition,case_id=case_id,selection=selection,
                             source_profile=str(path),source_profile_sha256=expected,
                             periodic_cut_native=-5,initial_peak_coordinate=0))
    return histories,selected


def axes_style():
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":8,
                         "axes.labelsize":8, "xtick.labelsize":7, "ytick.labelsize":7,
                         "axes.spines.top":False, "axes.spines.right":False,
                         "pdf.fonttype":42, "svg.fonttype":"none"})


def draw_organization(ax, folds):
    for key, marker, label, offset in (
            ("excess_active_fold", "o", "Excess active amount", -.18),
            ("width_fold", "s", "Core width", 0),
            ("peak_concentration_fold", "^", "Peak concentration", .18)):
        ax.scatter([float(r[key]) for r in folds], np.arange(5) + offset,
                   marker=marker, s=23, facecolor="#333333" if marker == "o" else "white",
                   edgecolor="#333333", lw=.8, label=label)
    labels = [NAMES[r["model"]] + "\n" + ("L₀ 1→2" if r["model"] == "spring"
              else f'{float(r["nominal_size_small"]):g}→{float(r["nominal_size_large"]):g}')
              for r in folds]
    ax.set(yticks=np.arange(5), yticklabels=labels, xlim=(.8,3.2),
           xticks=[1,2,3], ylim=(4.6,-.6), xlabel="Larger / smaller value")
    ax.tick_params(axis="y", length=0, labelsize=7)
    ax.axvline(1, color="#aaaaaa", ls="--", lw=.7)
    ax.set_title("A  Incumbent organization", loc="left", fontsize=9)


def draw_gain(ax, gains):
    for model, marker in zip(MODELS, ("o","s","D","^","v")):
        rows = [r for r in gains if r["model"] == model]
        ax.plot([r["relative_nominal_size"] for r in rows],
                [r["gain_over_own_reference"] for r in rows],
                color=COLORS[model], marker=marker, ms=2.6, lw=1, label=NAMES[model],
                mec="white", mew=.2)
    ax.axhline(1, color="#bdbdbd", ls="--", lw=.7)
    ax.set_yscale("log", base=2)
    ax.set(xlim=(.1,2.75), ylim=(.5,8), xticks=[.5,1,1.5,2,2.5],
           yticks=[.5,1,2,4,8], yticklabels=["½","1","2","4","8"],
           xlabel="Nominal size / own reference", ylabel="Gain / own reference (log₂)")
    ax.set_title("B  Initial opposing input gain", loc="left", fontsize=9)


def draw_recovery_folds(ax, rows):
    """Draw qualified within-model folds; callers supply explicit statuses.

    Canonical fields: model,size_native,restored_fraction,qualified. Clocks and
    raw fractions remain in the evidence CSV/caption. Unqualified cases appear
    as gaps; the exact own-reference case must qualify, rather than silently
    substituting a convenient reference.
    """
    all_values = []
    for model, marker in zip(MODELS, ("o","s","D","^","v")):
        block = sorted([r for r in rows if r["model"] == model],
                       key=lambda r:float(r["size_native"]))
        ref, = [r for r in block if float(r["size_native"]) == REFS[model]]
        assert ref["qualified"] and float(ref["restored_fraction"]) > 0
        values = [float(r["restored_fraction"]) / float(ref["restored_fraction"])
                  if r["qualified"] else np.nan for r in block]
        all_values.extend(v for v in values if np.isfinite(v))
        ax.plot([float(r["size_native"]) / REFS[model] for r in block], values,
                color=COLORS[model], marker=marker, ms=2.6, lw=1,
                mec="white", mew=.2, label=NAMES[model])
    ax.axhline(1, color="#bdbdbd", ls="--", lw=.7)
    ax.set(xlim=(.1,2.75), ylim=(min(0,np.floor(min(all_values))),
                               max(4,np.ceil(max(all_values)))),
           xticks=[.5,1,1.5,2,2.5],
           xlabel="Nominal size / own reference", ylabel="Restored fraction / own reference")
    ax.set_title("C  Local deficit restoration", loc="left", fontsize=9)


def draw_matched_profile(ax, profiles, model):
    for geometry, style, color, label in (
            ("interval reflected for display", "--", "#333333", "Reflected interval"),
            ("native ring", "-", COLORS[model], "Native ring")):
        rows = sorted([r for r in profiles if r["model"] == model and r["geometry"] == geometry],
                      key=lambda r:r["position"])
        ax.plot([r["position"] for r in rows], [r["active_concentration"] for r in rows],
                ls=style, color=color, lw=1.1, label=label)
    ax.set(xlim=(-15,15), ylim=(0,None), xticks=[-15,0,15])
    ax.set_title(NAMES[model], fontsize=8)
    ax.set_xlabel("Boundary coordinate (native)", fontsize=7)
    ax.set_ylabel("Native active concentration", fontsize=7)
    ax.tick_params(labelsize=6.6)


def draw_ring_map(ax, rows, model, key, norm, cmap, categorical=False):
    sizes = sorted({float(r["circumference"]) for r in rows})
    ring_components.map_cells(ax, rows, model, sizes, key, norm, cmap,
                              categorical=categorical, highlights=True)
    ax.set_title(NAMES[model], fontsize=8, pad=6)
    ax.set_xlabel("Shorter-route fraction d/C", fontsize=7)
    ax.set_ylabel("Sampled circumference C", fontsize=7)
    ax.tick_params(labelsize=6.3)


def draw_ring_history(ax, display, index, C, distance, vmax):
    x = display[f"example{index}_x"]
    p = display[f"example{index}_primary"]
    t = display[f"example{index}_time"]
    mesh = ax.pcolormesh(np.r_[x,x[0]+C], t, np.c_[p,p[:,0]],
                         cmap="magma", vmin=0, vmax=vmax, shading="nearest",
                         rasterized=True)
    ax.axvline(0, color="white", ls=":", lw=.6)
    ax.axvline(distance, color="#83d6e2", ls=":", lw=.6)
    ax.hlines(1000, -10, C-10, color="#83d6e2", ls="--", lw=.8)
    if C == 30:
        ax.axvspan(20,30,color="#d0d0d0")
        ax.text(25,2400,"Outside C30",rotation=90,ha="center",va="center",
                fontsize=6.5,color="#555555")
    ax.set(xlim=(-10,30), ylim=(4000,0), xticks=[-10,0,10,20,30],
           yticks=[0,1000,2000,4000], xlabel="Position from original peak",
           ylabel="Native time")
    labels = {1:"Cue-associated single peak",2:"Merged / bridged region",3:"Two distinct sites"}
    ax.set_title(f"E{index}  Holmes C{C} d{distance}\n{labels[index]}", fontsize=8, pad=7)
    ax.tick_params(labelsize=6.5)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    return mesh


def vector_trace(rows, C, distance):
    block = [r for r in rows if r["model"] == "holmes_model3"
             and float(r["circumference"]) == C
             and np.isclose(float(r["fraction"]),distance/C)]
    assert block
    result = {key:np.array([float(r[key]) for r in block])
              for key in ("time","vector_x","vector_y","relative_magnitude")}
    assert np.all(np.diff(result["time"]) > 0)
    return result


def draw_linked_vectors(axis_vector, axis_magnitude, rows):
    theta = np.linspace(0,2*np.pi,361)
    axis_vector.plot(np.cos(theta),np.sin(theta),color="#d1d1d1",lw=.7)
    axis_vector.axhline(0,color="#dddddd",lw=.5)
    axis_vector.axvline(0,color="#dddddd",lw=.5)
    for index,C,distance,color,style in (
            (1,30,12,ring_components.COLORS[1],"--"),
            (2,40,12,ring_components.COLORS[3],"-"),
            (3,40,16,ring_components.COLORS[2],"-")):
        z = vector_trace(rows,C,distance)
        keep = z["time"] <= 4000
        axis_vector.plot(z["vector_x"][keep],z["vector_y"][keep],color=color,ls=style,lw=1.2)
        for moment,marker in ((1000,"s"),(4000,"^")):
            j, = np.flatnonzero(z["time"] == moment)
            axis_vector.plot(z["vector_x"][j],z["vector_y"][j],marker=marker,ms=4,
                             color=color,mec="white",mew=.3)
        end, = np.flatnonzero(z["time"] == 4000)
        axis_vector.text(z["vector_x"][end]+(.08 if index==3 else -.12),
                         z["vector_y"][end]+.06,str(index),color=color,
                         fontsize=7,fontweight="bold")
        axis_magnitude.plot(z["time"][keep],z["relative_magnitude"][keep],
                            color=color,ls=style,lw=1.2,
                            label=f"{index}: cue {360*distance/C:g}°")
    axis_vector.plot(1,0,"o",mfc="white",mec="black",ms=4,zorder=8)
    axis_vector.set(xlim=(-1.15,1.15),ylim=(-1.15,1.15),xticks=[-1,0,1],
                    yticks=[-1,0,1],xlabel="Re P",ylabel="Im P")
    axis_vector.set_aspect("equal",adjustable="box")
    axis_vector.set_title("F  Activity-vector space (not cell position)",loc="left",fontsize=8,pad=7)
    axis_vector.tick_params(labelsize=6.5)
    z = vector_trace(rows,40,20)
    keep = z["time"] <= 4000
    axis_magnitude.plot(z["time"][keep],z["relative_magnitude"][keep],color="#707070",
                        ls=":",lw=1.2,label="Cancellation: cue 180°")
    axis_magnitude.axvline(1000,color="#888888",ls="--",lw=.6)
    axis_magnitude.axhline(1,color="#bdbdbd",ls=":",lw=.6)
    axis_magnitude.set(xlim=(0,4000),ylim=(0,1.15),xticks=[0,1000,2000,4000],
                       yticks=[0,.5,1],xlabel="Native time from cue onset",
                       ylabel="|P| / own initial magnitude")
    axis_magnitude.set_title("G  Oriented activity magnitude",loc="left",fontsize=8,pad=7)
    axis_magnitude.tick_params(labelsize=6.5)
    axis_magnitude.xaxis.label.set_size(7)
    axis_magnitude.yaxis.label.set_size(7)
    axis_magnitude.legend(loc="upper right",bbox_to_anchor=(1,.28),frameon=False,
                          fontsize=6.1,ncol=2,handlelength=1.5,columnspacing=.7)


def figure2_response(evidence, profiles, gainmap, recoverymap, recovery_pairs):
    fig = plt.figure(figsize=(9.4,10.6))
    fig.text(.51,.985,"Figure 2 · Ring organization and local response",ha="center",
             va="top",fontsize=11,weight="bold")
    fig.text(.51,.955,"Continued review sheet 1 of 2 · supported chemical rings",
             ha="center",fontsize=8)
    fig.text(.06,.925,"A  An actual matched native interval and ring",fontsize=9)
    for left,model in zip((.07,.385,.70),RING_MODELS):
        draw_matched_profile(fig.add_axes([left,.80,.235,.095]),profiles,model)
    fig.legend(handles=[
        Line2D([],[],color="#333333",ls="--",label="Interval L15, reflected for display"),
        Line2D([],[],color="#777777",ls="-",label="Independently prepared ring C30")],
        loc="center",bbox_to_anchor=(.51,.76),ncol=2,frameon=False,fontsize=7)
    largest = max(abs(float(r["ring_over_linear_K"])-1) for r in recovery_pairs)
    fig.text(.06,.728,
             f"All 13 C = 2L pairs: ring width/amount ≈ 2×; peak/matched gain ≈ 1×; recovery agrees within {100*largest:.2f}%.\n"
             "Independent native meshes/states. Matched gain: endpoint support 3 ↔ full ring support 6; challenges below use full support 3.",
             va="top",fontsize=7,linespacing=1.5)
    fig.text(.06,.678,"B  Initial challenged-site gain / own C30 antipodal gain",fontsize=9)
    for left,model in zip((.07,.385,.70),RING_MODELS):
        draw_ring_map(fig.add_axes([left,.438,.235,.21]),gainmap,model,
                      "gain_over_own_C30_antipodal",Normalize(0,3),plt.get_cmap("viridis"))
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,3),cmap="viridis"),
                      cax=fig.add_axes([.95,.438,.012,.21]))
    cb.set_ticks([0,1,2,3]);cb.ax.tick_params(labelsize=6.5)
    cb.set_label("Gain / own C30 antipodal gain",fontsize=7)
    fig.text(.06,.381,
             "One 0–3 scale: initial gain varies by <0.09% across placements at fixed C.\n"
             "Numbers link to sustained responses on sheet 2. Rows are discrete sampled circumferences.",
             va="top",fontsize=7,linespacing=1.5)
    fig.text(.06,.323,"C  Common active-loss response at t = 10 / the same size's basal response",fontsize=9)
    for left,model in zip((.07,.385,.70),RING_MODELS):
        draw_ring_map(fig.add_axes([left,.083,.235,.21]),recoverymap,model,
                      "recovery_over_basal",Normalize(0,1.1),plt.get_cmap("viridis"))
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1.1),cmap="viridis"),
                      cax=fig.add_axes([.95,.083,.012,.21]))
    cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6.5)
    cb.set_label("F at challenged t10 / own basal F",fontsize=7)
    fig.text(.06,.017,
             "Same original core/damage direction; native τ = 0.1778 (WP/Holmes), 0.005623 (Goryachev). Positive conditional F coexists with actual loss in every sham core.\n"
             "867 planned cells: 652 complete; × = 198 cue/core overlaps; dark gray = 17 unavailable Holmes C80 basals. Native input/time units are not cross-calibrated.",
             va="top",fontsize=7,linespacing=1.5)
    evidence.save(fig,"vision_aligned_figure2_response")


def figure2_competition(evidence, outcomes, display, vectors):
    fig = plt.figure(figsize=(9.4,10.6))
    fig.text(.51,.985,"Figure 2 · Sustained competition and global activity polarity",
             ha="center",va="top",fontsize=11,weight="bold")
    fig.text(.51,.954,"Continued review sheet 2 of 2 · full support 3; 1000 on, 13,000 off; native peaks WP 0.1 / G 6.4 / H 1",
             ha="center",fontsize=8)
    for left,model in zip((.07,.385,.70),RING_MODELS):
        draw_ring_map(fig.add_axes([left,.715,.235,.19]),outcomes,model,
                      "final_location_association",Normalize(-.5,6.5),
                      ListedColormap(ring_components.COLORS),categorical=True)
    fig.text(.06,.93,"D  Complete local-pattern context at t = 14,000",fontsize=9)
    handles = [Patch(facecolor=ring_components.COLORS[j],label=ring_components.LABELS[j])
               for j in range(5)]
    handles += [Line2D([],[],marker="x",ls="",color="#888888",label="Cue/core overlap excluded"),
                Patch(facecolor="#777777",label="Holmes C80 unavailable")]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.51,.672),ncol=3,
               frameon=False,fontsize=6.4,handlelength=1.3,columnspacing=1,labelspacing=.5)
    vmax = max(float(display[f"example{i}_primary"].max()) for i in (1,2,3))
    for left,index,C,distance in zip((.07,.385,.70),(1,2,3),(30,40,40),(12,12,16)):
        mesh = draw_ring_history(fig.add_axes([left,.355,.235,.225]),
                                  display,index,C,distance,vmax)
    cb = fig.colorbar(mesh,cax=fig.add_axes([.95,.355,.012,.225]))
    cb.ax.tick_params(labelsize=6.5);cb.set_label("Native active Rac",fontsize=7)
    fig.text(.06,.303,
             "1→2 changes size/material/architecture/routes and cue angle (144°→108°). Only 2→3 moves the cue in the identical full basal state (108°→144°).\n"
             "Holmes peak 1; cue off at 1000. Native histories show 0–4000; map endpoints are at 14,000. Gray lies outside C30; old peak remains intact at the periodic cut.",
             fontsize=7,linespacing=1.4)
    f = fig.add_axes([.10,.065,.35,.17]);g = fig.add_axes([.59,.065,.35,.17])
    draw_linked_vectors(f,g,vectors)
    fig.text(.50,.013,
             "P = M₁(t)/M₁(0): initial (1,0); magnitude may exceed 1. Squares: cue off; triangles: t = 4000.\n"
             "Near-zero resultant need not mean inactivity: antipodal Holmes retains two sites and nearly unchanged total Rac. See the linked full-grid/vector evidence.",
             ha="center",fontsize=7,linespacing=1.4)
    evidence.save(fig,"vision_aligned_figure2_competition")


def figure1(evidence, folds, gains, recovery, finite, primary, display):
    fig = plt.figure(figsize=(9.8,11.4))
    fig.text(.52,.992,"Figure 1 · Enlargement, local response and sustained competition",
             ha="center",va="top",fontsize=11,weight="bold")
    a = fig.add_axes([.135,.685,.205,.18]);draw_organization(a,folds)
    b = fig.add_axes([.45,.685,.205,.18]);draw_gain(b,gains)
    c = fig.add_axes([.765,.685,.205,.18]);draw_recovery_folds(c,recovery)
    for ax in (a,b,c):ax.tick_params(labelsize=6.8)
    fig.legend(*a.get_legend_handles_labels(),loc="upper left",
               bbox_to_anchor=(.12,.962),frameon=False,fontsize=6.8,labelspacing=.35)
    fig.legend(*b.get_legend_handles_labels(),loc="upper center",
               bbox_to_anchor=(.72,.962),frameon=False,ncol=3,fontsize=6.8,
               columnspacing=.9,handlelength=1.5)
    fig.text(.06,.631,
             "All five established linear models. B/C: own references L19, Holmes L20, Spring resting L₀ = 1; all 97 eligible sizes retained.\n"
             "C uses the common active-loss assay; native τ: WP/H 0.1778, G 0.005623, deBelly 5, Spring 0.01778. Compare size folds within each model.",
             va="top",fontsize=7,linespacing=1.5)
    fig.text(.06,.568,"D  After the same prolonged local cue and 13,000 cue-free units: every eligible size",fontsize=9)
    old_color,new_color = "#4678a3","#d77c35"
    final_values = [float(r[key]) for r in primary
                    for key in ("post13000_old_excess_ratio","post13000_challenger_excess_ratio")]
    lim = (min(-.04,min(final_values)-.025),max(1.10,max(final_values)+.04))
    for index,(left,model) in enumerate(zip((.08,.265,.45,.635,.82),MODELS),1):
        ax = fig.add_axes([left,.38,.14,.14])
        rows = sorted([r for r in primary if r["model"] == model],
                      key=lambda r:float(r["size_native"]))
        for key,color,label in (("post13000_old_excess_ratio",old_color,"Original core"),
                                ("post13000_challenger_excess_ratio",new_color,"Mirrored core")):
            x = [float(r["size_native"])/REFS[model] for r in rows]
            y = [float(r[key]) for r in rows]
            ax.plot(x,y,color=color,marker="o",ms=2.5,lw=1,label=label)
            overlap = [(xx,yy) for xx,yy,r in zip(x,y,rows) if r["cue_core_overlap"]=="True"]
            if overlap:
                ax.scatter(*np.asarray(overlap).T,s=22,facecolor="white",edgecolor=color,lw=.9,zorder=5)
            ax.plot(x[-1],y[-1],"s",mfc="none",mec="black",ms=5,mew=.7)
        n_bad = sum(r["model"]==model and r["status"]=="preparation_ineligible" for r in finite)
        ax.set(xlim=(.1,2.75),ylim=lim,xticks=[.5,1,2,2.5],yticks=[0,.5,1],
               xlabel="Size / own reference")
        ax.set_title(f"{NAMES[model]}\n{len(rows)} eligible; {n_bad} ineligible",fontsize=7.1,pad=7)
        ax.axhline(0,color="#d0d0d0",lw=.6);ax.axhline(1,color="#dddddd",lw=.6,ls=":")
        ax.tick_params(labelsize=6.4);ax.xaxis.label.set_size(6.7)
        ax.annotate(f"E{index}",(x[-1],max(float(rows[-1]["post13000_old_excess_ratio"]),
                                        float(rows[-1]["post13000_challenger_excess_ratio"]))),
                    xytext=(-2,7),textcoords="offset points",ha="right",fontsize=6.5)
        if index==1:
            ax.set_ylabel("Core excess / initial\nincumbent excess",fontsize=7)
            handles,labels = ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc="center",bbox_to_anchor=(.79,.55),ncol=2,
               frameon=False,fontsize=7,handlelength=1.5)
    fig.text(.06,.314,
             "One common final-stage y scale; no negative excess is clipped. Fixed original and mirrored regions are not tracked site identities.\n"
             "Squares select the prespecified largest eligible size shown below; open circles flag cue/core overlap. Ineligible preparations receive no response.",
             va="top",fontsize=7,linespacing=1.5)
    fig.text(.06,.265,"E  Exact native histories at the largest eligible size in every model",fontsize=9)
    for index,(left,model) in enumerate(zip((.08,.265,.45,.635,.82),MODELS),1):
        ax = fig.add_axes([left,.065,.14,.15]);data = display[model]
        mesh = ax.pcolormesh(data["fractional_position"],data["time_over_duration"],
                             data["primary"],cmap="magma",vmin=0,
                             vmax=data["full_history_max"],shading="nearest",rasterized=True)
        ax.axhline(1,color="#83d6e2",ls="--",lw=.8)
        ax.set(xlim=(0,1),ylim=(4,0),xticks=[0,.5,1],yticks=[0,1,2,4])
        ax.tick_params(labelsize=6.4)
        size = data["size_native"];T = data["cue_duration_native"]
        ax.set_title(f"E{index}  {NAMES[model]}\n"
                     f"{'L₀' if model=='spring' else 'L'} = {size:g}; cue T = {T:g}",fontsize=7.1,pad=7)
        if index==1:ax.set_ylabel("Time / native cue duration\n(clocks differ by model)",fontsize=7)
        if model=="goryachev":
            row = next(r for r in primary if r["model"]==model
                       and float(r["size_native"])==size)
            ax.text(.56,2.6,
                    f'Cue-end peaks:\nold {float(row["withdrawal_old_peak"]):.2f}\n'
                    f'new {float(row["withdrawal_challenger_peak"]):.2f}',
                    color="#dddddd",fontsize=6.0,ha="center",va="center",linespacing=1.4)
        cb = fig.colorbar(mesh,cax=fig.add_axes([left,.015,.14,.006]),orientation="horizontal")
        cb.set_ticks([0,data["full_history_max"]])
        cb.set_ticklabels(["0",f'{data["full_history_max"]:.3g}'])
        cb.ax.tick_params(labelsize=5.8,length=2)
    fig.text(.52,.037,"Fractional native coordinate: incumbent at 0, challenge at 1; Spring uses its material coordinate",
             ha="center",fontsize=7)
    fig.text(.52,-.022,
             "Zero-based native concentration bars differ by model. Cyan marks withdrawal; displayed histories cover 0–4T; D uses the full 13,000-unit follow-up.\n"
             "Fixed endpoint support 3 (Spring 0.2); native peaks WP 0.1 / G 6.4 / deBelly 3.6 / H 1 / Spring 4. Half/double durations remain linked sensitivities.",
             ha="center",fontsize=7,linespacing=1.5)
    evidence.save(fig,"vision_aligned_figure1_linear")


def figure3(evidence, tables, histories):
    size_colors = {20.:"#3976a8",40.:"#c96926",60.:"#8560a4"}
    styles = {.001:"--",.01:"-"}
    fig = plt.figure(figsize=(9.8,12.6))
    fig.text(.51,.987,"Figure 3 · Replicated activity-polarity memory under local fluctuations",
             ha="center",va="top",fontsize=11,weight="bold")
    fig.text(.06,.951,
             "12 independent realizations per model/C/strength; C = 20, 40, 60; duration 32,000 native units. Added input is positive, with both declared means retained.\n"
             "Fixed local statistics: μ / native cue peak = 0.001 or 0.01; Gaussian spatial parameter 1.5; normalized OU knots τ = 50, spacing 5. Native units differ by model.",
             va="top",fontsize=7.2,linespacing=1.5)
    handles = [Line2D([],[],color=c,label=f"C{C:g}") for C,c in size_colors.items()]
    handles += [Line2D([],[],color="#444444",ls="--",label="f = 0.001"),
                Line2D([],[],color="#444444",ls="-",label="f = 0.01"),
                Line2D([],[],color="#444444",ls=":",label="High-f mean control"),
                Line2D([],[],color="#999999",ls="-.",label="Zero-input control")]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.51,.916),ncol=4,
               frameon=False,fontsize=7,columnspacing=1.1,handlelength=1.8)
    fig.text(.06,.857,"A  Initial-axis memory ⟨Re P(t)⟩: both strengths; bands are 95% seed-bootstrap intervals",fontsize=9)
    memory_values = [float(r[k]) for r in tables["traces"]
                     if r["metric"]=="initial_memory" and r["condition"]=="noise"
                     for k in ("mean","ci_lower","ci_upper") if np.isfinite(float(r[k]))]
    memory_lim = (min(-1.1,np.floor(min(memory_values)*10)/10),
                  max(1.2,np.ceil(max(memory_values)*10)/10))
    R_values = [float(r[k]) for r in tables["traces"] if r["metric"]=="R"
                for k in ("mean","ci_lower","ci_upper") if np.isfinite(float(r[k]))]
    Rmax = max(1.2,np.ceil(max(R_values)*10)/10)
    for left,model in zip((.08,.39,.70),RING_MODELS):
        a = fig.add_axes([left,.70,.235,.12])
        for C,color in size_colors.items():
            for f,style in styles.items():
                rows = noise_series(tables["traces"],model,C,f,"noise","initial_memory")
                plot_noise_series(a,rows,"time","mean",color,style)
        a.set(xlim=(0,32),ylim=memory_lim,xticks=[0,16,32],yticks=[-1,0,1],
              xlabel="Native time (×10³)")
        a.set_title(NAMES[model],fontsize=8,pad=7);a.axhline(0,color="#dddddd",lw=.6)
        if model=="wave_pinning":a.set_ylabel("Mean initial-axis memory",fontsize=7)
        a.tick_params(labelsize=6.6);a.xaxis.label.set_size(7)
    fig.text(.06,.647,"B  High strength f = 0.01: mean magnitude ⟨|P|⟩ with deterministic-mean and zero-input controls",fontsize=9)
    for left,model in zip((.08,.39,.70),RING_MODELS):
        b = fig.add_axes([left,.49,.235,.12])
        for C,color in size_colors.items():
            for condition,style,lw in (("noise","-",1.1),("mean",":",1),("zero","-.",.65)):
                f = 0 if condition=="zero" else .01
                rows = noise_series(tables["traces"],model,C,f,condition,"R")
                plot_noise_series(b,rows,"time","mean",
                                  "#999999" if condition=="zero" else color,style,
                                  ci=condition=="noise",lw=lw)
        b.set(xlim=(0,32),ylim=(0,Rmax),xticks=[0,16,32],xlabel="Native time (×10³)")
        if model=="wave_pinning":b.set_ylabel("Mean relative magnitude",fontsize=7)
        b.tick_params(labelsize=6.6);b.xaxis.label.set_size(7)
    fig.text(.06,.437,"C  Full-window vector lag correlation: both strengths; magnitude contributes to this unmasked observable",fontsize=9)
    for left,model in zip((.08,.39,.70),RING_MODELS):
        c = fig.add_axes([left,.29,.235,.11])
        for C,color in size_colors.items():
            for f,style in styles.items():
                rows = noise_series(tables["lags"],model,C,f,"noise")
                finite = [float(r["correlation"]) for r in rows if np.isfinite(float(r["correlation"]))]
                assert all(abs(v)<=1+1e-6 for v in finite)
                plot_noise_series(c,rows,"lag","correlation",color,style)
        c.set(xlim=(0,8),ylim=(-1.05,1.05),xticks=[0,4,8],yticks=[-1,0,1],
              xlabel="Native lag (×10³)")
        c.axhline(0,color="#dddddd",lw=.6)
        if model=="wave_pinning":c.set_ylabel("Vector lag correlation",fontsize=7)
        c.tick_params(labelsize=6.6);c.xaxis.label.set_size(7)
    fig.text(.06,.227,"D  Exact C20, f = 0.01 local histories: condition medoids and the equal-mean Goryachev control",fontsize=9)
    max_by_model = {model:max(float(r["primary"].max()) for r in histories if r["model"]==model)
                    for model in RING_MODELS}
    for index,(left,data) in enumerate(zip((.07,.30,.53,.76),histories),1):
        ax = fig.add_axes([left,.025,.18,.15]);u=data["primary"];x=data["x"]
        mesh = ax.pcolormesh(np.r_[x,x[0]+20],data["time"]/1000,np.c_[u,u[:,0]],
                             cmap="magma",vmin=0,vmax=max_by_model[data["model"]],
                             shading="nearest",rasterized=True)
        ax.axvline(0,color="white",ls=":",lw=.5)
        ax.set(xlim=(-5,15),ylim=(32,0),xticks=[-5,0,5,15],yticks=[0,16,32])
        ax.tick_params(labelsize=6.3)
        title = f'D{index}  {NAMES[data["model"]]}\n' + ("Equal-mean control" if data["condition"]=="mean" else "Noisy-input medoid")
        ax.set_title(title,fontsize=7.1,pad=7)
        if index==1:ax.set_ylabel("Native time (×10³)",fontsize=7)
        cb = fig.colorbar(mesh,cax=fig.add_axes([left,-.018,.18,.006]),orientation="horizontal")
        cb.set_ticks([0,max_by_model[data["model"]]])
        cb.set_ticklabels(["0",f'{max_by_model[data["model"]]:.3g}'])
        cb.ax.tick_params(labelsize=5.8,length=2)
    fig.text(.51,-.005,"Boundary arclength from the initial peak (native units); periodic cut at −5 keeps that peak intact",
             ha="center",fontsize=7)
    fig.text(.51,-.054,
             "Zero-based native concentration bars; both Goryachev histories share a scale. Local histories are individual selected realizations, not ensemble profiles.\n"
             "Read correlations with magnitude, gated-direction coverage and physical peak displacement in linked evidence; no fitted decay time or cell velocity is assigned.",
             ha="center",fontsize=7,linespacing=1.5)
    evidence.save(fig,"vision_aligned_figure3_noise")


def render_ring_component():
    axes_style()
    evidence = Evidence()
    evidence.source(__file__)
    pairs,profiles = matched_native_profiles(evidence)
    gains,outcomes,display,vectors = retained_ring_data(evidence)
    basals,recovery,recovery_pairs = common_ring_recovery(evidence)
    figure2_response(evidence,profiles,gains,recovery,recovery_pairs)
    figure2_competition(evidence,outcomes,display,vectors)
    evidence.export("vision_aligned_matched_native_profiles.csv",profiles)
    evidence.export("vision_aligned_ring_conditional_recovery.csv",recovery)
    path = evidence.manifest(dict(checkpoint_complete=False,
                                   completed_actual_component="Figure 2",
                                   ring_basal_states=len(basals),
                                   ring_planned_cells=len(recovery),
                                   new_native_integrations=0))
    print(path)


def render_deterministic_component():
    axes_style()
    evidence = Evidence();evidence.source(__file__)
    folds,basal_profiles = linear_organization(evidence)
    gains = linear_gain(evidence)
    recovery = common_linear_recovery(evidence)
    finite,primary,display,selected = linear_finite(evidence)
    pairs,ring_profiles = matched_native_profiles(evidence)
    gainmap,outcomes,ring_display,vectors = retained_ring_data(evidence)
    basals,recoverymap,recovery_pairs = common_ring_recovery(evidence)
    figure1(evidence,folds,gains,recovery,finite,primary,display)
    figure2_response(evidence,ring_profiles,gainmap,recoverymap,recovery_pairs)
    figure2_competition(evidence,outcomes,ring_display,vectors)
    evidence.export("vision_aligned_linear_gain.csv",gains)
    evidence.export("vision_aligned_linear_recovery.csv",recovery)
    evidence.export("vision_aligned_linear_longest_examples.csv",selected)
    arrays = {model+"_"+key:np.asarray(value)
              for model,block in display.items() for key,value in block.items()}
    evidence.export_npz("vision_aligned_linear_native_display.npz",arrays)
    evidence.export("vision_aligned_matched_native_profiles.csv",ring_profiles)
    evidence.export("vision_aligned_ring_conditional_recovery.csv",recoverymap)
    path = evidence.manifest(dict(checkpoint_complete=False,
                                   completed_actual_component="Figures 1 and 2",
                                   linear_gain_cases=len(gains),
                                   linear_common_recovery_cases=len(recovery),
                                   linear_finite_primary_cases=len(primary),
                                   linear_finite_sensitivities=30,
                                   linear_ineligible_preparations=13,
                                   ring_basal_states=len(basals),
                                   ring_planned_cells=len(recoverymap),
                                   new_native_integrations=0))
    print(path)


def render_full_checkpoint():
    axes_style()
    evidence = Evidence();evidence.source(__file__)
    noise = noise_tables(evidence,DEV/"vision_noise_20261017/supporting/data")
    histories,noise_selected = noise_history_examples(evidence,noise)
    folds,basal_profiles = linear_organization(evidence)
    gains = linear_gain(evidence)
    recovery = common_linear_recovery(evidence)
    finite,primary,display,selected = linear_finite(evidence)
    pairs,ring_profiles = matched_native_profiles(evidence)
    gainmap,outcomes,ring_display,vectors = retained_ring_data(evidence)
    basals,recoverymap,recovery_pairs = common_ring_recovery(evidence)
    figure1(evidence,folds,gains,recovery,finite,primary,display)
    figure2_response(evidence,ring_profiles,gainmap,recoverymap,recovery_pairs)
    figure2_competition(evidence,outcomes,ring_display,vectors)
    figure3(evidence,noise,histories)
    evidence.export("vision_aligned_linear_gain.csv",gains)
    evidence.export("vision_aligned_linear_recovery.csv",recovery)
    evidence.export("vision_aligned_linear_longest_examples.csv",selected)
    arrays = {model+"_"+key:np.asarray(value)
              for model,block in display.items() for key,value in block.items()}
    evidence.export_npz("vision_aligned_linear_native_display.npz",arrays)
    evidence.export("vision_aligned_matched_native_profiles.csv",ring_profiles)
    evidence.export("vision_aligned_ring_conditional_recovery.csv",recoverymap)
    evidence.export("vision_aligned_noise_examples.csv",noise_selected)
    path = evidence.manifest(dict(checkpoint_complete=True,
                                   working_figures=3,figure2_review_sheets=2,
                                   linear_gain_cases=len(gains),
                                   linear_common_recovery_cases=len(recovery),
                                   linear_finite_primary_cases=len(primary),
                                   linear_finite_sensitivities=30,
                                   linear_ineligible_preparations=13,
                                   ring_basal_states=len(basals),
                                   ring_planned_cells=len(recoverymap),
                                   noise_cases=243,noise_realizations=216,
                                   noise_realizations_per_group=12,
                                   noise_scalar_plot_cadence=5,
                                   new_native_integrations=0))
    print(path)


if __name__ == "__main__":
    if sys.argv[1:] == ["--figure2-only"]:
        render_ring_component()
    elif sys.argv[1:] == ["--deterministic-only"]:
        render_deterministic_component()
    elif sys.argv[1:] == ["--full"]:
        render_full_checkpoint()
    else:
        raise SystemExit("Final rendering requires the complete, independently checked campaign adapters.")
