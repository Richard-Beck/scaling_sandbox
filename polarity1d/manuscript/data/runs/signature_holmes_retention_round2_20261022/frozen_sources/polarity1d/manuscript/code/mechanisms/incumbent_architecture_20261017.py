"""Provisional basal active-region architecture draft; no new integration."""
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
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/length_mechanisms_20261008/dev/incumbent_architecture_20261017"
INITIAL = ROOT / "data/runs/initial_polarity_length_corrected_20261008"
POOL = ROOT / "data/runs/inactive_pool_figure_exports_20261008"
STEM = "incumbent_architecture"
MODELS = {"wave_pinning": ("Wave-pinning", "#007D73"), "goryachev": ("Goryachev", "#893F82")}
LENGTH_COLORS = {9: "#0072B2", 21: "#8B3F82", 39: "#D55E00"}
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path, expected=None):
    path = Path(path)
    digest = sha(path)
    if expected is not None:
        assert digest == expected, f"Source hash mismatch: {path}"
    SOURCES[str(path.relative_to(ROOT))] = digest
    return path


def read(path):
    return json.loads(source(path).read_text())


def write_json(name, data):
    path = dev_output_path(OUT, name)
    path.write_text(json.dumps(data, indent=2) + "\n")
    return path


def write_csv(name, rows):
    path = dev_output_path(OUT, name)
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def measure():
    provenance = read(INITIAL / "report_provenance.json")
    with source(INITIAL / "steady_state_metrics.csv", provenance["summary_csv_sha256"]).open(newline="") as file:
        raw = list(csv.DictReader(file))
    raw = [row for row in raw if row["model"] in MODELS]
    assert len(raw) == 41
    records, profiles, checks = [], [], []
    pool_provenance = read(POOL / "provenance.json")
    with np.load(source(POOL / "evidence_D0p1.npz", pool_provenance["inputs"]["evidence_D0p1.npz"])) as archive:
        snapshots = {f"{model}_L{L}_steady": archive[f"{model}_L{L}_steady"].copy()
                     for model in MODELS for L in LENGTH_COLORS}
    for model in MODELS:
        rows = sorted([row for row in raw if row["model"] == model], key=lambda row: float(row["domain_length"]))
        reference = next(row for row in rows if float(row["domain_length"]) == 9)
        width9, peak9 = float(reference["old_pole_width"]), float(reference["old_pole_peak_active"])
        assert len(rows) == (20 if model == "wave_pinning" else 21)
        for row in rows:
            length = float(row["domain_length"])
            width, peak = float(row["old_pole_width"]), float(row["old_pole_peak_active"])
            records.append(dict(model=model, length_native=length, length_relative_to_L9=length / 9,
                incumbent_width_native=width, width_reference_L9_native=width9,
                width_relative_to_L9=width / width9, incumbent_peak_native=peak,
                peak_reference_L9_native=peak9, peak_relative_to_L9=peak / peak9,
                incumbent_width_fraction_of_length=width / length,
                incumbent_active_amount_native=float(row["old_pole_total_active"]),
                incumbent_active_amount_above_background_native=float(row["old_pole_excess_active"]),
                basal_background_active_native=float(row["background_active"]),
                checkpoint=row["checkpoint"], checkpoint_sha256=row["checkpoint_sha256"]))
            if int(length) not in LENGTH_COLORS:
                continue
            with np.load(source(ROOT / row["checkpoint"], row["checkpoint_sha256"])) as cp:
                state, x = cp["state"].copy(), cp["x"].copy()
            np.testing.assert_array_equal(state, snapshots[f"{model}_L{int(length)}_steady"])
            np.testing.assert_allclose(x, (np.arange(len(x)) + 0.5) * length / len(x), atol=1e-14)
            u = state[0]
            np.testing.assert_allclose(u.max(), peak, rtol=2e-15)
            level = u.min() + 0.5 * np.ptp(u)
            index = int(np.flatnonzero(u < level)[0])
            measured_width = x[index - 1] + (level - u[index - 1]) / (u[index] - u[index - 1]) * (x[index] - x[index - 1])
            np.testing.assert_allclose(measured_width, width, rtol=2e-14)
            for i in range(len(x)):
                profiles.append(dict(model=model, length_native=length, cell_index=i,
                    position_native=float(x[i]), active_native=float(u[i]),
                    inactive_native=float(state[1, i])))
            checks.append(dict(model=model, length=length, exact_checkpoint_and_export_identity=True,
                native_peak_and_half_contrast_width_recomputed=True,
                whole_pool_density_at_basal=float(state.sum() / len(x))))
    for row in records:
        if row["length_native"] == 9:
            assert row["length_relative_to_L9"] == row["width_relative_to_L9"] == row["peak_relative_to_L9"] == 1
    return records, profiles, checks


def render(records, profiles):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 9, "axes.titlesize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(7.35, 5.65))
    fig.subplots_adjust(left=0.105, right=0.975, top=0.875, bottom=0.105, wspace=0.36, hspace=0.63)
    for col, (model, (label, _)) in enumerate(MODELS.items()):
        ax = axes[0, col]
        for length, color in LENGTH_COLORS.items():
            selected = [row for row in profiles if row["model"] == model and row["length_native"] == length]
            ax.plot([row["position_native"] for row in selected], [row["active_native"] for row in selected], color=color, lw=1.8)
        ax.set(title=f"{label}: basal active profiles", ylabel="Active u (native concentration)",
               xlabel="Position, x (native units)", xlim=(0, 40 if col == 0 else 4),
               ylim=(0, 1.25 if col == 0 else 125),
               xticks=[0, 10, 20, 30, 40] if col == 0 else [0, 1, 2, 3, 4],
               yticks=[0, 0.4, 0.8, 1.2] if col == 0 else [0, 40, 80, 120])
        if col == 1:
            ax.text(0.98, 0.93, "Incumbent-end view", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, color="#555555")
        ax.grid(axis="y", alpha=0.12, lw=0.6)
    for model, (label, color) in MODELS.items():
        selected = [row for row in records if row["model"] == model]
        x = [row["length_relative_to_L9"] for row in selected]
        for ax, key in ((axes[1, 0], "width_relative_to_L9"), (axes[1, 1], "peak_relative_to_L9")):
            ax.plot(x, [row[key] for row in selected], color=color, marker="o", ms=3.0, mew=0, lw=1.8, label=label)
    for ax in axes[1]:
        ax.plot([0, 4.55], [0, 4.55], color="#999999", ls="--", lw=1.0, zorder=0, label="Proportional to length")
        ax.axhline(1, color="#999999", ls=":", lw=1.0, zorder=0, label="Constant")
        ax.axvline(1, color="#dddddd", lw=0.7, zorder=0)
        ax.set(xlim=(0.2, 4.55), ylim=(0, 4.7), xticks=[1, 2, 3, 4], yticks=[0, 1, 2, 3, 4],
               xlabel="Length relative to L = 9")
        ax.grid(axis="y", alpha=0.12, lw=0.6)
    axes[1, 0].set(title="Width of the incumbent active region", ylabel="Width relative to L = 9")
    axes[1, 1].set(title="Peak active concentration", ylabel="Peak relative to L = 9")
    handles, labels = axes[1, 0].get_legend_handles_labels()
    axes[1, 0].legend(handles[:2], labels[:2], loc="upper left", frameon=False, fontsize=8)
    handles, labels = axes[1, 1].get_legend_handles_labels()
    axes[1, 1].legend(handles[2:], labels[2:], loc="upper left", frameon=False, fontsize=7.5)
    fig.legend(handles=[Line2D([], [], color=color, lw=1.8, label=f"L = {length}")
                        for length, color in LENGTH_COLORS.items()],
               loc="upper center", bbox_to_anchor=(0.54, 0.985), ncol=3, frameon=False, fontsize=8.5)
    for ax, letter in zip(axes.flat, "ABCD"):
        box = ax.get_position()
        fig.text(box.x0 - 0.073, box.y1 + 0.027, letter, fontsize=12, fontweight="bold")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    clipped = []
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box = artist.get_window_extent(renderer)
            if box.width and box.height and not (box.x0 >= 0 and box.x1 <= fig.bbox.x1 and box.y0 >= 0 and box.y1 <= fig.bbox.y1):
                clipped.append(artist.get_text())
    assert not clipped, clipped
    files = []
    for extension in ("png", "pdf", "svg"):
        path = dev_output_path(OUT, f"{STEM}.{extension}")
        fig.savefig(path, dpi=320, facecolor="white")
        files.append(path)
    plt.close(fig)
    return files


def main():
    source(Path(__file__))
    source(Path(__file__).with_name("length_dev_output_layout.py"))
    records, profiles, checks = measure()
    files = render(records, profiles)
    files += [write_csv("incumbent_architecture_measurements.csv", records),
              write_csv("basal_active_profiles.csv", profiles)]
    caption = (
        "Discussion draft: added material is organized differently in two basal polarity models. "
        "A–B, Exact active profiles at L = 9, 21 and 39 in native units; Wave-pinning shows the full "
        "domains, while Goryachev shows x = 0–4 to resolve its narrow incumbent site. C–D, Incumbent "
        "width and peak concentration versus length, each divided by its own model’s L = 9 value. "
        "Width is the first interpolated boundary-connected crossing of the profile minimum plus "
        "half its contrast. All eligible lengths are included: Wave-pinning L = 5–39 (20 lengths), "
        "Goryachev L = 3–39 (21). Dashed and dotted guides indicate proportional-to-length and "
        "constant responses; they are not fitted laws. For a 4.33× increase in length from 9 to 39, "
        "Wave-pinning width rises 4.23× while its peak rises 1.02×; Goryachev width rises 1.12× while "
        "its peak rises 4.23×. Goryachev width first decreases to a minimum near L = 13, then rises; "
        "the shortest eligible Wave-pinning peak is 0.88× its L = 9 value. The conserved mean u + v "
        "is fixed (2.2683 in Wave-pinning; 2.75033 in Goryachev), so total available material is "
        "proportional to length by construction. The result concerns its organization. Active-region width "
        "does not measure a physical protrusion, force or migration persistence. No new integrations."
    )
    path = dev_output_path(OUT, "caption.txt")
    path.write_text(caption + "\n")
    files.append(path)
    summary = {model: {"shortest_length": min(row["length_native"] for row in records if row["model"] == model),
        "L39": next(row for row in records if row["model"] == model and row["length_native"] == 39)} for model in MODELS}
    files.append(write_json("verification.json", dict(passed=True, new_integrations=0,
        source_measurements_sha256_verified=True, all_eligible_model_metrics_included=41,
        six_profile_checkpoints_rechecked=checks,
        raw_measurement_denominators_and_normalized_values_exported=True,
        reference="Each model's exact L=9 width and peak; length divided by 9.",
        guides="Reference slopes only; no regressions or scaling-law fit.",
        text_within_page_bounds=True, visual_inspection="Pending final PNG inspection.", summary=summary)))
    files.append(write_json("evidence_gaps.json", dict(
        status="Provisional discussion draft of site organization only.",
        message="Wave-pinning mainly widens its active region; Goryachev mainly raises active concentration.",
        exceptions="The shortest eligible cells differ from the long-cell pattern; Goryachev width is mildly nonmonotone, with a minimum near L13.",
        limits=["Half-contrast biochemical width is not physical protrusion width.",
                "Peak and width are descriptive; they do not alone establish whole-pool allocation or a causal mechanism.",
                "More total available material follows from extending a fixed-density mass-conserved setup.",
                "This figure does not establish forces, migration persistence or biological WGD fold calibration."],
        why_this_display="Native profiles show what changes spatially; common reference-normalized width and peak curves test widening and concentration changes across every eligible sampled length.")))
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer="finished",
        artifact_root="../..", status="Provisional discussion draft", new_integrations=0,
        source_sha256=SOURCES, renderer_sha256=sha(Path(__file__)), caption=caption,
        sha256={dev_artifact_key(OUT, file): sha(file) for file in files},
        reproduction="PYTHONDONTWRITEBYTECODE=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__)))
    write_json("manifest.json", manifest)
    print(json.dumps(dict(output=str(OUT), passed=True, source_hashes=len(SOURCES), summary=summary), indent=2))


if __name__ == "__main__":
    main()
