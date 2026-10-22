"""Supporting four-panel candidate from existing audited chemistry evidence.

No integration, canonical-data write, or report mutation is performed. This is
a provisional supporting draft, not an agreed main-text figure.

Run with the existing scaling_sandbox Python, with one BLAS/OpenMP thread.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import subprocess
import sys

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"
sys.dont_write_bytecode = True

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "code/core"), str(ROOT / "code/core/src")]
from models.goryachev import Goryachev
from models.wave_pinning import WavePinning

OUT = ROOT / "figures/length_mechanisms_20261008/dev/chemistry_main_text_candidate_20261017"
POOL = ROOT / "data/runs/inactive_pool_figure_exports_20261008"
RAW = ROOT / "data/runs/inactive_pool_restoring_response_20261008"
COMPACT = ROOT / "data/runs/length_mechanisms_20261008/comparative_revision"
INITIAL = ROOT / "data/runs/initial_polarity_length_corrected_20261008"
REPORT = ROOT / "figures/length_mechanisms_20261008"
MODELS = {
    "goryachev": ("Goryachev", Goryachev(), "#893F82"),
    "wave_pinning": ("Wave-pinning", WavePinning(), "#007D73"),
}
LENGTH_COLORS = {9: "#0072B2", 21: "#8B3F82", 39: "#D55E00"}
STEM = "chemistry_main_text_candidate"
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def register(path, expected=None):
    path = Path(path)
    actual = sha(path)
    if expected is not None:
        assert actual == expected, f"Source hash mismatch: {path}"
    SOURCES[str(path.relative_to(ROOT))] = actual
    return path


def read(path):
    return json.loads(register(path).read_text())


def rows(path):
    with register(path).open(newline="") as file:
        return list(csv.DictReader(file))


def load(path):
    with np.load(register(path), allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def write_json(name, value):
    path = dev_output_path(OUT, name)
    path.write_text(json.dumps(value, indent=2) + "\n")
    return path


def write_csv(name, records):
    path = dev_output_path(OUT, name)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path


def reaction(model, u, v, stimulus=0.0):
    u, v, stimulus = np.broadcast_arrays(u, v, stimulus)
    result = model.reaction(
        0.0, np.zeros(u.size), np.stack([u.ravel(), v.ravel()]), stimulus.ravel()
    )[0]
    return result.reshape(u.shape)


def direct_reaction(name, u, v):
    if name == "goryachev":
        return (0.33 / 101 * u * u + 0.67 / 101 * u) * v - 0.01733 * u
    return (0.067 + u * u / (1.0 + u * u)) * v - u


def check_report():
    final = read(REPORT / "revision_20261014_draft_alignment/equation_typesetting/final_audit.json")
    assert final["passed"]
    pdf = register(REPORT / "length_mechanisms.pdf", final["pdf_sha256"])
    text = subprocess.run(["pdftotext", str(pdf), "-"], check=True,
                          capture_output=True, text=True).stdout
    pages = text.split("\f")
    required = {
        3: "Inactive depletion raises the barrier to a new Goryachev pole",
        4: "Long Goryachev cells resist damage and rebound more strongly",
    }
    for page, title in required.items():
        assert title in " ".join(pages[page - 1].split())
    assets = REPORT / "revision_20261014/visual_palette/assets"
    checked_artifacts = 0
    for stem in ("chemical_hurdles", "incumbent_restoration"):
        manifest = read(assets / f"{stem}.manifest.json")
        assert manifest["passed"] and manifest["new_integrations"] == 0
        for name, digest in manifest["artifacts"].items():
            register(assets / name, digest)
            checked_artifacts += 1
        for name, digest in manifest["sources"].items():
            register(ROOT / name, digest)
    return dict(pdf_sha256=sha(pdf), pages=list(required),
                current_report_figure_artifacts_verified=checked_artifacts)


def audit_and_measure():
    report_identity = check_report()
    pool_provenance = read(POOL / "provenance.json")
    for name, digest in pool_provenance["inputs"].items():
        register(POOL / name, digest)
    for name, digest in pool_provenance["sources"].items():
        register(ROOT / name, digest)
    for name, digest in pool_provenance["models"].items():
        register(ROOT / f"code/core/models/{name}.py", digest)
    comparative_provenance = read(COMPACT / "provenance.json")
    for name, digest in comparative_provenance["sources"].items():
        register(ROOT / name, digest)
    assert read(COMPACT / "independent_observable_audit.json")["passed"]
    assert read(RAW / "results/independent_audit.json")["passed"]
    initial_provenance = read(INITIAL / "report_provenance.json")
    metrics_path = register(INITIAL / "steady_state_metrics.csv",
                            initial_provenance["summary_csv_sha256"])
    metrics = rows(metrics_path)
    assert len(metrics) == len(initial_provenance["checkpoints"]) == 97
    # Check the whole current checkpoint catalog, then independently compare
    # the 36 plotted chemistry cases to their exact saved state arrays below.
    for item in initial_provenance["checkpoints"]:
        register(ROOT / item["path"], item["sha256"])
    steady = load(POOL / "evidence_D0p1.npz")
    compact = load(COMPACT / "compact_evidence.npz")
    response_rows = rows(POOL / "response_vs_length.csv")
    manifest = read(RAW / "manifest.json")
    selected_jobs = [job for job in manifest["jobs"] if job["duration"] == 0.1]
    assert len(selected_jobs) == 36
    records, profiles, cases = [], [], []
    max_source_identity_error = 0.0
    for name, (_, model, _) in MODELS.items():
        # The compact maps are an additional independent cross-check of the
        # native basal reaction identity used for every plotted initial state.
        U, V = np.meshgrid(compact[f"{name}_reaction_u"], compact[f"{name}_reaction_v"])
        native = reaction(model, U, V)
        np.testing.assert_allclose(native, direct_reaction(name, U, V), rtol=2e-13, atol=2e-14)
        np.testing.assert_allclose(native, compact[f"{name}_basal_net_reaction"], rtol=2e-13, atol=2e-14)
        rr = sorted([row for row in response_rows if row["model"] == name
                     and float(row["duration"]) == 0.1
                     and float(row["epsilon"]) == 0.005], key=lambda row: float(row["length"]))
        assert len(rr) == 18
        baseline_v = next(float(row["opposite_inactive"]) for row in metrics
                          if row["model"] == name and float(row["domain_length"]) == 9)
        baseline_gain = float(rr[0]["reaction_gain_per_damage"])
        for row in rr:
            length = float(row["length"])
            L = int(length)
            metric = next(row for row in metrics if row["model"] == name
                          and float(row["domain_length"]) == length)
            job = next(job for job in selected_jobs if job["case"]["model"] == name
                       and job["case"]["length"] == length)
            case = job["case"]
            assert case["checkpoint"] == metric["checkpoint"]
            assert case["checkpoint_sha256"] == metric["checkpoint_sha256"]
            cp = load(register(ROOT / metric["checkpoint"], metric["checkpoint_sha256"]))
            state = cp["state"]
            np.testing.assert_array_equal(state, steady[f"{name}_L{L}_steady"])
            np.testing.assert_allclose(state[1, -1], float(metric["opposite_inactive"]), rtol=2e-15)
            base = RAW / row["archive"]
            pulse_meta, sham_meta = read(base.with_suffix(".json")), read(base.with_name("epsilon_0.json"))
            assert pulse_meta["job"] == sham_meta["job"] == job
            assert pulse_meta["epsilon"] == 0.005 and sham_meta["epsilon"] == 0
            z = load(register(base, pulse_meta["archive_sha256"]))
            sham = load(register(base.with_name("epsilon_0.npz"), sham_meta["archive_sha256"]))
            for archive, epsilon in ((z, 0.005), (sham, 0.0)):
                np.testing.assert_array_equal(archive["originalstate"], state)
                np.testing.assert_allclose(archive["removed_inactive"], epsilon * state[1], rtol=0, atol=1e-13)
                np.testing.assert_allclose(archive["pulse_active"][0], state[0], rtol=0, atol=1e-12)
                np.testing.assert_allclose(archive["pulse_inactive"][0], (1 - epsilon) * state[1], rtol=0, atol=1e-12)
                returned = np.stack([archive["pulse_active"][-1],
                                     archive["pulse_inactive"][-1] + epsilon * state[1]])
                np.testing.assert_array_equal(archive["returnedstate"], returned)
                np.testing.assert_allclose(archive["recovery_active"][0], returned[0], rtol=0, atol=1e-10)
                np.testing.assert_allclose(archive["recovery_inactive"][0], returned[1], rtol=0, atol=1e-10)
                assert archive["pulse_time"][-1] == 0.1 and archive["recovery_time"][0] == 0
            x = z["x"]
            dx = length / len(x)
            np.testing.assert_allclose(x, (np.arange(len(x)) + 0.5) * dx, atol=1e-14)
            width = float(metric["old_pole_width"])
            np.testing.assert_allclose(width, job["width"], rtol=2e-14)
            weights = np.clip(width - np.arange(len(x)) * dx, 0.0, dx)
            u0, v0 = state
            # Equal local cue: the additive initial rate per unit local s is v.
            # Verify at every grid point for several s, including native 0.1.
            for strength in (0.01, 0.1, 1.0):
                delta = reaction(model, u0, v0, strength) - reaction(model, u0, v0)
                error = float(np.max(np.abs(delta - strength * v0)))
                max_source_identity_error = max(max_source_identity_error, error)
                np.testing.assert_allclose(delta, strength * v0, rtol=2e-12, atol=2e-14)
            pert_u, pert_v = z["recovery_active"][0], z["recovery_inactive"][0]
            sham_u, sham_v = sham["recovery_active"][0], sham["recovery_inactive"][0]
            # Independent finite-volume integrals at the exact return state.
            damage = float((sham_u - pert_u) @ weights)
            excess_reaction = float((direct_reaction(name, pert_u, pert_v)
                                     - direct_reaction(name, sham_u, sham_v)) @ weights)
            gain = excess_reaction / damage
            # CSV was written by subtracting already-integrated core amounts;
            # the cancellation order can differ in the final few digits.
            np.testing.assert_allclose(damage, float(row["active_damage_at_release"]), rtol=2e-9, atol=2e-12)
            np.testing.assert_allclose(excess_reaction, float(row["excess_net_activation_at_release"]), rtol=2e-9, atol=2e-12)
            np.testing.assert_allclose(gain, float(row["reaction_gain_per_damage"]), rtol=2e-9, atol=2e-9)
            np.testing.assert_allclose(float(row["excess_net_activation_at_release"])
                                       / float(row["active_damage_at_release"]),
                                       float(row["reaction_gain_per_damage"]), rtol=2e-15)
            if L in LENGTH_COLORS:
                prefix = f"{name}_L{L}_"
                np.testing.assert_array_equal(compact[prefix + "recovery_time"], z["recovery_time"])
                exported_gain = (z["recovery_core_reaction"] - sham["recovery_core_reaction"]) / float(row["active_damage_at_release"])
                np.testing.assert_allclose(compact[prefix + "reaction_per_initial_damage"], exported_gain, rtol=2e-13, atol=2e-13)
                for i, (xi, active, inactive) in enumerate(zip(x, u0, v0)):
                    profiles.append(dict(model=name, length=length, cell_index=i,
                                         position_native=xi, active_native=active,
                                         inactive_native=inactive, checkpoint=metric["checkpoint"]))
            cue_raw = float(metric["opposite_inactive"])
            gain_raw = float(row["reaction_gain_per_damage"])
            records.append(dict(model=name, length_native=length,
                                opposite_inactive_native=cue_raw,
                                initial_added_activation_per_unit_local_cue_native=cue_raw,
                                initial_added_activation_at_local_cue_0p1_native=0.1 * cue_raw,
                                initial_activation_relative_to_L9=cue_raw / baseline_v,
                                activation_L9_reference_native=baseline_v,
                                reaction_gain_per_damage_native=gain_raw,
                                restoring_gain_relative_to_L9=gain_raw / baseline_gain,
                                restoration_L9_reference_native=baseline_gain,
                                fixed_core_width_native=width,
                                active_core_damage_native=float(row["active_damage_at_release"]),
                                excess_core_reaction_native=float(row["excess_net_activation_at_release"]),
                                fractional_core_damage=float(row["fractional_damage"]),
                                inactive_removal_fraction=0.005, hold_duration_native=0.1,
                                checkpoint=metric["checkpoint"], checkpoint_sha256=metric["checkpoint_sha256"],
                                perturbation_archive=str(base.relative_to(ROOT))))
            cases.append(dict(model=name, length=length, source_slot=job["slot"],
                              checkpoint_and_export_array_identity=True, pointwise_removal_return=True,
                              fixed_core_damage_and_reaction_recomputed=True))
    assert len(records) == 36 and len(profiles) == 6 * 201
    for name in MODELS:
        first = next(row for row in records if row["model"] == name and row["length_native"] == 9)
        assert first["initial_activation_relative_to_L9"] == first["restoring_gain_relative_to_L9"] == 1
    audit = dict(passed=True, new_integrations=0, current_report_identity=report_identity,
                 verified_catalog_checkpoint_hashes=97, exact_chemistry_checkpoints=36,
                 raw_archives_rechecked=72, initial_source_identity="f(u,v;s) - f(u,v;0) = s v",
                 initial_source_identity_max_absolute_error=max_source_identity_error,
                 native_basal_maps_and_closed_form_identity=True,
                 normalization="Each model divided by its own exact L=9 value; raw values retained in CSV.",
                 cases=cases)
    return records, profiles, audit


def render(records, profiles):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "lines.linewidth": 1.7,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })
    fig = plt.figure(figsize=(7.35, 5.65))
    grid = fig.add_gridspec(2, 2, left=0.102, right=0.90, top=0.875,
                           bottom=0.10, wspace=0.60, hspace=0.63)
    axes = [fig.add_subplot(grid[i, j]) for i in range(2) for j in range(2)]
    twins = []
    for col, (name, (label, _, _)) in enumerate(MODELS.items()):
        ax = axes[col]
        twin = ax.twinx()
        twin.spines["right"].set_visible(True)
        twins.append(twin)
        for length, color in LENGTH_COLORS.items():
            selected = [row for row in profiles if row["model"] == name and row["length"] == length]
            x = [row["position_native"] for row in selected]
            ax.plot(x, [row["active_native"] for row in selected], color=color)
            twin.plot(x, [row["inactive_native"] for row in selected], color=color,
                      ls=(0, (4, 2)), lw=1.5)
        ax.set(title=f"{label}: basal profiles", xlabel="Position, x (native units)",
               ylabel="Active, u (native conc.)", xlim=(0, 40), ylim=(0, 125 if col == 0 else 1.25))
        twin.set(ylabel="Inactive, v (native conc.)", ylim=(0, 0.30 if col == 0 else 2.0))
        ax.set_xticks([0, 10, 20, 30, 40])
        ax.yaxis.set_major_locator(MaxNLocator(4))
        twin.yaxis.set_major_locator(MaxNLocator(4))
        ax.grid(axis="y", alpha=0.12, lw=0.6)
    legend = [Line2D([], [], color=color, label=f"L = {length}")
              for length, color in LENGTH_COLORS.items()]
    legend += [Line2D([], [], color="#333333", label="Active u"),
               Line2D([], [], color="#333333", ls=(0, (4, 2)), label="Inactive v")]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.985),
               ncol=5, frameon=False, handlelength=2, columnspacing=1.3)
    for name, (label, _, color) in MODELS.items():
        selected = [row for row in records if row["model"] == name]
        x = [row["length_native"] for row in selected]
        for ax, key in ((axes[2], "initial_activation_relative_to_L9"),
                        (axes[3], "restoring_gain_relative_to_L9")):
            ax.plot(x, [row[key] for row in selected], color=color, marker="o",
                    ms=3.0, mew=0, label=label)
        axes[3].annotate(f"{selected[-1]['restoring_gain_relative_to_L9']:.1f}×",
                         (39, selected[-1]["restoring_gain_relative_to_L9"]),
                         xytext=(-2, 7),
                         textcoords="offset points", ha="right", color=color, fontsize=8.5)
    axes[2].set(title="Initial opposing-cue activation", ylabel="Added activation / local cue\n(relative to L = 9)",
                ylim=(0, 1.12), yticks=[0, 0.25, 0.5, 0.75, 1])
    axes[2].legend(frameon=False, loc="lower left", handlelength=1.8)
    axes[3].set(title="Immediate incumbent restoration", ylabel="Excess reaction / damage\n(relative to L = 9)",
                ylim=(0, 18), yticks=[0, 5, 10, 15])
    for ax in axes[2:]:
        ax.set(xlabel="Cell length, L (native units)", xlim=(8, 40), xticks=[9, 21, 39])
        ax.axhline(1, color="#aaaaaa", ls=":", lw=0.8, zorder=0)
        ax.grid(axis="y", alpha=0.12, lw=0.6)
    # Place panel letters in figure coordinates, outside the plot and above
    # the y-axis tick region, so twin labels cannot obscure them.
    for letter, ax in zip("ABCD", axes):
        box = ax.get_position()
        fig.text(box.x0 - 0.065, box.y1 + 0.027, letter, fontsize=12, fontweight="bold")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    clipped = []
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box = artist.get_window_extent(renderer)
            if box.width > 0 and box.height > 0 and not (box.x0 >= bounds.x0 and box.x1 <= bounds.x1
                                                        and box.y0 >= bounds.y0 and box.y1 <= bounds.y1):
                clipped.append(artist.get_text())
    assert not clipped, f"Text exceeds page bounds: {clipped}"
    paths = []
    for extension in ("png", "pdf", "svg"):
        path = dev_output_path(OUT, f"{STEM}.{extension}")
        fig.savefig(path, dpi=320, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths, dict(text_within_page_bounds=True, panel_count=4,
                       native_profile_axes_have_separate_scales=True,
                       size_inches=[7.35, 5.65], png_dpi=320,
                       visual_inspection="Pending manual inspection of rendered PNG.")


def main():
    register(Path(__file__))
    register(Path(__file__).with_name("length_dev_output_layout.py"))
    records, profiles, audit = audit_and_measure()
    outputs = [write_csv("length_response_measurements.csv", records),
               write_csv("basal_profiles.csv", profiles)]
    paths, layout = render(records, profiles)
    outputs.extend(paths)
    caption = (
        "Provisional supporting draft: two local chemical effects that could help retain polarity. "
        "A–B, Exact basal active u (solid, left axes) and inactive v (dashed, right axes) profiles "
        "at L = 9, 21 and 39, in model-native units with separately scaled concentration axes. "
        "C, Initial added activation per unit equal local opposing cue at the opposite-end cell; "
        "the exact source identity f(u,v;s) − f(u,v;0) = s v makes this the initial local inactive "
        "concentration. D, Immediate excess basal reaction per active-core damage after a 0.5% "
        "whole-inactive-profile removal for 0.1 native time units and exact pointwise return. "
        "Both measurements in D subtract a matched sham within the fixed original half-contrast "
        "core; reaction excludes diffusion. C–D show all 18 audited lengths from 9 to 39, normalized "
        "to each model’s own L = 9 value. At L = 39, Goryachev reaches 0.333× initial activation "
        "and 15.7× immediate restoring gain; Wave-pinning reaches 0.99990× and 1.28×. These local "
        "diagnostics leave challenger formation, survival and takeover unresolved. They do not "
        "measure late recovery, global stability, a causal partition of reversal time, migration "
        "persistence or a calibrated biological WGD effect. No new simulations; this is not a "
        "selected main-text figure."
    )
    caption_path = dev_output_path(OUT, "caption.txt")
    caption_path.write_text(caption + "\n")
    outputs.append(caption_path)
    endpoint_summary = {}
    for name in MODELS:
        last = next(row for row in records if row["model"] == name and row["length_native"] == 39)
        endpoint_summary[name] = {key: last[key] for key in (
            "initial_activation_relative_to_L9", "restoring_gain_relative_to_L9",
            "opposite_inactive_native", "reaction_gain_per_damage_native")}
    outputs.append(write_json("evidence_gaps.json", dict(
        status="Supporting candidate; not an agreed main-text figure.",
        supported="Two local chemical diagnostics that could contribute to polarity retention.",
        evidence_gap="The figure does not determine whether an alternative front can form, survive cue withdrawal, or take over the incumbent.",
        additional_limits=[
            "The cue response is pointwise and initial; fixed physical or relative cue-support integrals are not plotted.",
            "The incumbent perturbation creates no competing front and removes the entire inactive profile.",
            "The restoring gain is immediate, excludes diffusion and depends on this specified damage direction and core.",
            "The diagnostic is not a global stability eigenvalue or a late recovery rate.",
            "The coupled effects do not quantify additive causal contributions to reversal time.",
            "Native 1D polarity readouts do not simulate migration persistence or calibrate biological WGD fold changes.",
        ], endpoint_summary=endpoint_summary)))
    audit["layout"] = layout
    audit_path = write_json("verification.json", audit)
    outputs.append(audit_path)
    outputs.append(write_json("methods.json", dict(
        status="Provisional supporting draft", integration="none",
        original_perturbation=read(POOL / "provenance.json")["definition"],
        source_figures="Current length_mechanisms.pdf Figures 3–4, hash-verified against final equation-typesetting audit.",
        profiles="Exact saved original states; 201 finite-volume cell centres. No interpolation or fitted profiles.",
        activation="Exact initial added active reaction per unit local stimulus at the opposite-end cell; equals v at that cell.",
        restoration="(perturbed minus sham) core-integrated net reaction at return / (sham minus perturbed) core-integrated active amount at return.",
        core="Each length's original boundary-connected half-contrast core, kept fixed; partial final cell included.",
        normalization="Each model's C and D curves divided by their own raw L=9 reference; both raw denominators and values are exported.",
        comparison="Normalized trends compare length sensitivity within each model; raw native rate scales differ between models.",
        all_checkpoints="97 current catalog hashes verified; 36 exact chemistry checkpoints compared to frozen exports and original perturbation/sham arrays.")))
    timestamp = datetime.now(timezone.utc).isoformat()
    manifest = dict(
        generated_at=timestamp, producer="finished", status="Provisional supporting candidate",
        artifact_root="../..", integration="none", new_integrations=0,
        renderer=str(Path(__file__).relative_to(ROOT)), renderer_sha256=sha(Path(__file__)),
        caption=caption, source_sha256=SOURCES,
        sha256={dev_artifact_key(OUT, path): sha(path) for path in outputs},
        current_report_sha256=audit["current_report_identity"]["pdf_sha256"],
        reproduction="PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__)),
    )
    write_json("manifest.json", manifest)
    print(json.dumps(dict(output=str(OUT), passed=True, source_hashes=len(SOURCES),
                          endpoint_summary=endpoint_summary), indent=2))


if __name__ == "__main__":
    main()
