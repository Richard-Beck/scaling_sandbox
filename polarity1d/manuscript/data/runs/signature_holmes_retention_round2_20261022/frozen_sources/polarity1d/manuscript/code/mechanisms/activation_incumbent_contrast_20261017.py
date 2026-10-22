"""Small discussion draft: challenged-end buildup and incumbent coexistence.

Read six retained trajectories only; no integrations or canonical writes.
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
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "figures/length_mechanisms_20261008"
OUT = REPORT / "dev/activation_incumbent_contrast_20261017"
GOR = ROOT / "data/runs/length_revision_20261009"
BUDGET = ROOT / "data/runs/debelly_mca_speed_decomposition_20261010"
ENDPOINTS = REPORT / "revision_20261014_draft_alignment/endpoint_native_cue_series.csv"
STEM = "activation_incumbent_contrast"
COLORS = ["#0072B2", "#8B3F82", "#D55E00"]
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


def series():
    traces, checks, records = {}, [], []
    manifest = read(GOR / "manifest.json")
    for slot, length in enumerate((9, 21, 39)):
        path = GOR / f"workers/{slot:02d}.npz"
        meta = read(path.with_suffix(".json"))
        job = meta["job"]
        case = job["case"]
        assert job == manifest["jobs"][slot]
        assert case["model"] == "goryachev" and case["length"] == length
        assert case["condition"] == "fixed_relative" and case["response_amplitude"] == 6.4
        assert case["width_fraction_at_reference"] == 0.2 and job["duration"] > 180
        cp = source(ROOT / case["checkpoint"], case["checkpoint_sha256"])
        with np.load(cp) as checkpoint:
            basal = checkpoint["state"].copy()
        with np.load(source(path, meta["archive_sha256"])) as z:
            np.testing.assert_array_equal(z["originalstate"], basal)
            mask = z["pulse_time"] <= 180
            time = z["pulse_time"][mask].copy()
            active = z["pulse_active"][mask][:, [0, -1]].copy()
            np.testing.assert_array_equal(z["pulse_active"], z["pulse_field_active"])
            np.testing.assert_allclose(active[0], basal[0, [0, -1]], rtol=2e-13, atol=1e-15)
            x = z["pulse_x"].copy()
        traces[("goryachev", length)] = dict(time=time, old=active[:, 0], new=active[:, 1],
            basal=float(basal[0, 0]), duration=job["duration"], time_unit="native units",
            x=x[[0, -1]], cue_peak=6.4, native_field_tail="not assessed for this maintained-cue trajectory")
        checks.append(dict(model="goryachev", length=length, archive_sha256_verified=True,
                           exact_original_checkpoint=True, native_endpoint_cells=[0, 200],
                           common_window_before_withdrawal=True))
    budget_manifest = read(BUDGET / "manifest.json")
    assert budget_manifest["geometry"] == "fixed_relative"
    with source(ENDPOINTS).open(newline="") as file:
        endpoint_rows = list(csv.DictReader(file))
    for length in (23, 31, 39):
        slot = budget_manifest["lengths_um"].index(length)
        path = BUDGET / f"workers/{slot:03d}.npz"
        meta = read(path.with_suffix(".json"))
        assert meta["passed"] and meta["parameters"]["center_boundary_motion"] is False
        check = meta["checks"][0]
        assert check["length_um"] == length
        cp = source(ROOT / check["checkpoint"], meta["source_sha256"][check["checkpoint"]])
        original_meta = read((ROOT / check["archive"]).with_suffix(".json"))
        case = original_meta["job"]["case"]
        assert case["model"] == "debelly" and case["length"] == length
        assert case["condition"] == "fixed_relative" and case["amplitude"] == 3.6
        assert case["width_fraction_at_reference"] == 0.2
        assert original_meta["archive_sha256"] == meta["source_sha256"][check["archive"]]
        with np.load(source(path, meta["npz_sha256"])) as z:
            time = z[f"L{length}_time"].copy()
            active = z[f"L{length}_rac"][:, [0, -1]].copy()
            x = z[f"L{length}_x"].copy()
        exported = [row for row in endpoint_rows if int(row["length_um"]) == length]
        np.testing.assert_array_equal(time, [float(row["elapsed_time_s"]) for row in exported])
        np.testing.assert_array_equal(active[:, 0], [float(row["old_Rac"]) for row in exported])
        np.testing.assert_array_equal(active[:, 1], [float(row["new_Rac"]) for row in exported])
        assert time[-1] == check["cue_duration_s"] > 180
        assert all(float(row["cue_duration_s"]) == time[-1] for row in exported)
        with np.load(cp) as z:
            np.testing.assert_array_equal(active[0], z["state_rac"][[0, -1]])
        mask = time <= 180
        traces[("debelly", length)] = dict(time=time[mask], old=active[mask, 0], new=active[mask, 1],
            basal=float(active[0, 0]), duration=check["cue_duration_s"], time_unit="s",
            x=x[[0, -1]], cue_peak=3.6, native_field_tail=check["native_field_tail_pass"])
        checks.append(dict(model="debelly", length=length, archive_sha256_verified=True,
                           corrected_independent_endpoints=True, exact_original_checkpoint=True,
                           full_endpoint_export_exactly_matches_native_arrays=True,
                           common_window_before_withdrawal=True,
                           native_field_tail_pass=check["native_field_tail_pass"]))
    for (model, length), trace in traces.items():
        assert len(trace["time"]) == 361
        np.testing.assert_array_equal(trace["time"], np.arange(0, 180.5, 0.5))
        relative = trace["old"] / trace["basal"]
        assert relative[0] == 1.0
        trace["relative"] = relative
        for i, time in enumerate(trace["time"]):
            records.append(dict(model=model, length_native=length, elapsed_time=float(time),
                time_unit=trace["time_unit"], challenged_activity_native=float(trace["new"][i]),
                incumbent_activity_native=float(trace["old"][i]),
                incumbent_basal_native=trace["basal"], incumbent_fraction_of_basal=float(relative[i]),
                challenged_cell_index=200, incumbent_cell_index=0,
                incumbent_reference_position_native=float(trace["x"][0]),
                challenged_reference_position_native=float(trace["x"][1]),
                relative_cue_support=0.2, cue_peak_native=trace["cue_peak"],
                retained_cue_duration=trace["duration"], native_field_tail_pass=trace["native_field_tail"]))
    return traces, checks, records


def render(traces):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(7.35, 5.6))
    grid = fig.add_gridspec(2, 2, left=0.105, right=0.975, bottom=0.10,
                           top=0.81, hspace=0.29, wspace=0.40)
    axes = np.array([[fig.add_subplot(grid[i, j]) for j in range(2)] for i in range(2)])
    for col, (model, lengths, label, protocol, unit) in enumerate([
        ("goryachev", (9, 21, 39), "Goryachev", "Maintained cue; peak 6.4", "native units"),
        ("debelly", (23, 31, 39), "Corrected deBelly", "Selected finite pulses; peak 3.6", "s"),
    ]):
        for length, color in zip(lengths, COLORS):
            trace = traces[(model, length)]
            axes[0, col].plot(trace["time"], trace["new"], color=color, lw=1.8)
            axes[1, col].plot(trace["time"], trace["relative"], color=color, lw=1.8)
        box = axes[0, col].get_position()
        center = (box.x0 + box.x1) / 2
        fig.text(center, 0.966, label, ha="center", va="top", fontsize=11, fontweight="bold")
        fig.text(center, 0.927, protocol, ha="center", va="top", fontsize=8.5)
        handles = [Line2D([], [], color=color, lw=1.8, label=f"L = {length}")
                   for length, color in zip(lengths, COLORS)]
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(center, 0.895),
                   frameon=False, ncol=3, fontsize=8, handlelength=1.6, columnspacing=1.0)
        species = "u" if col == 0 else "Rac"
        axes[0, col].set(ylabel=f"Challenged-end {species}\n(native concentration)",
                         ylim=(0, 22.5 if col == 0 else 5), yticks=[0, 10, 20] if col == 0 else [0, 2, 4])
        axes[1, col].set(ylabel=f"Incumbent {species}\n(fraction of basal)", ylim=(0, 1.075), yticks=[0, 0.5, 1],
                         xlabel=f"Time after cue onset ({unit})")
        axes[1, col].axhline(1, color="#999999", ls=":", lw=0.8, zorder=0)
        for row in range(2):
            axes[row, col].set(xlim=(0, 180), xticks=[0, 60, 120, 180])
            axes[row, col].grid(axis="y", alpha=0.13, lw=0.6)
        axes[0, col].tick_params(labelbottom=False)
    # Keep the inset wholly above all main-trajectory values in its horizontal
    # footprint; the main direct histories must remain continuously visible.
    inset = axes[0, 0].inset_axes([0.08, 0.68, 0.35, 0.24])
    for length in (9, 21, 39):
        trace = traces[("goryachev", length)]
        footprint = (trace["time"] >= 0.08 * 180) & (trace["time"] <= 0.43 * 180)
        assert np.max(trace["new"][footprint]) < 0.60 * 22.5
    for length, color in zip((9, 21, 39), COLORS):
        trace = traces[("goryachev", length)]
        mask = trace["time"] <= 10
        inset.plot(trace["time"][mask], trace["new"][mask], color=color, lw=1.2)
    inset.set(xlim=(0, 10), ylim=(0, 3), xticks=[0, 10], yticks=[0, 3])
    inset.set_title("First 10 native units", fontsize=7.4, pad=3)
    inset.tick_params(labelsize=7, length=2, pad=1)
    for ax, letter in zip(axes.flat, "ABCD"):
        box = ax.get_position()
        fig.text(box.x0 - 0.074, box.y1 + 0.012, letter, fontsize=12, fontweight="bold")
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
    traces, checks, records = series()
    files = render(traces)
    path = dev_output_path(OUT, "source_series.csv")
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    files.append(path)
    caption = (
        "Discussion draft: challenged-end buildup and persistent incumbent activity are distinct observations. "
        "A–B, Raw primary activity at the challenged right endpoint-adjacent cell. C–D, Activity at the "
        "incumbent left endpoint-adjacent cell divided by its own exact basal value; the normalized axes "
        "are shared. Columns show Goryachev (u) and corrected independent-endpoint deBelly (Rac). "
        "Colors indicate short, intermediate and long cells, with actual lengths given in each column. "
        "All curves show the common first 0–180 cue-on interval within each model, before any withdrawal. "
        "Goryachev uses maintained relative-20% cues of native peak 6.4; deBelly uses selected finite "
        "relative-20% pulses of peak 3.6. Native time and seconds, concentrations and doses are not "
        "calibrated across models. Goryachev responds immediately (inset), but longer cells show much "
        "less substantial challenger buildup during this window and retain more incumbent activity. "
        "deBelly develops strong local challenged-end Rac while incumbent Rac persists, especially at "
        "L = 39. No front-establishment threshold is imposed. These trajectories do not prove migration "
        "persistence, failed takeover or isolated causation. Final relaxed outcomes are outside the "
        "display: the selected deBelly L = 23 and 39 trajectories fail their final native-field-tail "
        "checks; L = 31 passes. No new integrations."
    )
    caption_path = dev_output_path(OUT, "caption.txt")
    caption_path.write_text(caption + "\n")
    files.append(caption_path)
    summary = []
    for (model, length), trace in traces.items():
        summary.append(dict(model=model, length=length, time_180_unit=trace["time_unit"],
                            challenged_activity_180=float(trace["new"][-1]),
                            incumbent_fraction_180=float(trace["relative"][-1]),
                            cue_duration=trace["duration"], native_field_tail_pass=trace["native_field_tail"]))
    verification = dict(passed=True, new_integrations=0, trajectory_count=6, plotted_samples=2166,
        common_window="0 to 180, at the native 0.5-unit or 0.5-s saved times; all before withdrawal",
        normalization="Old endpoint activity / that endpoint's own exact basal activity; all time-zero fractions equal one.",
        cross_model_limits="Model-native concentration, time and cue peak have no shared biological calibration.",
        visual_inspection="Pending rendered PNG inspection.", text_within_page_bounds=True,
        checks=checks, endpoint_summary=summary)
    files.append(write_json("verification.json", verification))
    files.append(write_json("evidence_gaps.json", dict(
        status="Small discussion draft; not an integrated or selected main-text figure.",
        message="Length-dependent challenged-end buildup differs from persistence of the existing incumbent.",
        limits=["No universal activity threshold or dynamic establishment gate is introduced.",
                "The Goryachev inset exposes its immediate small response; it is not nonresponsive.",
                "Strong local Rac activation does not establish a surviving spatial front or takeover.",
                "A common cue-on window removes withdrawal timing from the displayed interval, but does not match model kinetics, doses or length ranges.",
                "The deBelly L23 and L39 final native-field-tail failures remain separate from their observed early cue-on response.",
                "No migration persistence or isolated causal mechanism is established."],
        why_this_display="Direct paired trajectories keep both challenger and incumbent visible at the same elapsed times; raw challenger values avoid normalizing away differences, while basal-normalized incumbent values expose retention without conflating different initial amplitudes.")))
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer="finished",
        status="Small provisional discussion draft", artifact_root="../..", new_integrations=0,
        renderer=str(Path(__file__).relative_to(ROOT)), renderer_sha256=sha(Path(__file__)),
        source_sha256=SOURCES, sha256={dev_artifact_key(OUT, file): sha(file) for file in files},
        caption=caption, reproduction="PYTHONDONTWRITEBYTECODE=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__)))
    write_json("manifest.json", manifest)
    print(json.dumps(dict(output=str(OUT), passed=True, source_hashes=len(SOURCES), endpoint_summary=summary), indent=2))


if __name__ == "__main__":
    main()
