"""Paired stationary-cue review pages from the two accepted arms; no native work."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
GROUP = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
OUT = GROUP / "paired_stationary_review"
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring")
NAMES = dict(wave_pinning="Wave-pinning", goryachev="Goryachev", holmes_model3="Holmes", spring="Spring")
BASENAMES = dict(wave_pinning="wave_pinning_paired_vectors", goryachev="goryachev_paired_vectors", holmes_model3="holmes_paired_vectors", spring="spring_paired_vectors")
ARMS = ("fixed_angle", "common_initial_physical_arc")
SCALE, PITCH, R_GATE = .24, .95, .1
ORANGE, BLACK = "#df7b25", "#111111"
EXPECTED = {
    "standard_angles_sizes/supporting/data/stationary_sizes_summary.csv": "761bf294ad1056620724d2fc8ddc9ca9c2bcd10beed9e8dd2fb3071152bde0d8",
    "standard_angles_sizes/supporting/data/stationary_sizes_traces.npz": "cf07adf328d05325342d07d8468a8823454a9e4efb628ca407f5c06992249d7b",
    "standard_arclength_sizes/supporting/data/stationary_arclength_summary.csv": "9b00954df3a3ff50f65de37ea9f3d82968bab10196b8cdc10a5d054d5329e278",
    "standard_arclength_sizes/supporting/data/stationary_arclength_traces.npz": "bbac92256b738f4b84f95440c716cabfa583f06ee07c52a7cdcafd9eaff5714f",
    "standard_arclength_sizes/supporting/data/stationary_two_arms_all480.csv": "b07dd2b161cb4361510c90db0f182ce1e5e8aab8ea8f0072c0b3959cc6de6629",
}
SOURCES = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dest(name, group="metadata"):
    path = OUT / name if name.endswith(".png") or name == "README.md" else OUT / "supporting" / group / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write(name, value):
    dest(name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def wrap(degrees):
    return float(np.rad2deg(np.angle(np.exp(1j * np.deg2rad(degrees)))))


def degree_label(value):
    return f"{round(value):d}°" if abs(value - round(value)) < 1e-6 else f"{value:.2f}°"


def load_values():
    for relative, digest in EXPECTED.items():
        path = GROUP / relative
        assert sha(path) == digest, f"Accepted source changed: {path}"
        SOURCES[str(path)] = digest
    union = pd.read_csv(GROUP / "standard_arclength_sizes/supporting/data/stationary_two_arms_all480.csv")
    assert len(union) == 480 and union.unique_protocol_id.nunique() == 355
    protocol_ids = {(r.experiment_arm, int(r.condition_id)): r.unique_protocol_id for r in union.itertuples()}
    shared = set(union.groupby("unique_protocol_id").size().loc[lambda x: x == 2].index)
    records, maximum_error = [], 0.
    for arm, folder, prefix in ((ARMS[0], "standard_angles_sizes", "stationary_sizes"),
                                (ARMS[1], "standard_arclength_sizes", "stationary_arclength")):
        frame = pd.read_csv(GROUP / folder / "supporting/data" / (prefix + "_summary.csv"))
        assert len(frame) == 240 and frame.status.eq("complete").all()
        with np.load(GROUP / folder / "supporting/data" / (prefix + "_traces.npz")) as traces:
            for row in frame.itertuples():
                cid = int(row.condition_id)
                phase, time = traces[f"c{cid}_phase"], traces[f"c{cid}_time"]
                off = int(np.flatnonzero(phase == 0)[-1])
                assert time[off] == row.cue_duration
                P = traces[f"c{cid}_P_real"] + 1j * traces[f"c{cid}_P_imag"]
                withdrawal, final = complex(P[off]), complex(P[-1])
                maximum_error = max(maximum_error, abs(abs(withdrawal) - row.R_withdrawal), abs(abs(final) - row.R_final))
                assert abs(abs(withdrawal) - row.R_withdrawal) < 1e-12 and abs(abs(final) - row.R_final) < 1e-12
                for value, angle in ((withdrawal, row.angle_withdrawal_deg), (final, row.angle_final_deg)):
                    if abs(value) >= R_GATE:
                        assert abs(wrap(np.rad2deg(np.angle(value)) - angle)) < 1.000001e-6
                target = float(row.target_angle_in_P_coordinates_deg)
                initial_error = abs(wrap(-target))
                assert initial_error > 0.
                final_error = None if abs(final) < R_GATE else abs(wrap(np.rad2deg(np.angle(final)) - target))
                closed = None if final_error is None else float(1. - final_error / initial_error)
                placement_row = ("pi", "pi_2", "pi_4", "pi_8").index(row.angle_label) if arm == ARMS[0] else 3 - int(row.initial_arc_grid_index)
                times = (.1, 1., 10., 100., 1000.) if row.model == "spring" else (1., 10., 100., 300., 1000.)
                column = times.index(row.cue_duration)
                uid = protocol_ids[(arm, cid)]
                value = dict(arm=arm, arm_row=ARMS.index(arm), condition_id=cid, unique_protocol_id=uid,
                    shared_with_other_arm=uid in shared, model=row.model, size_parameter=float(row.size_parameter),
                    size_label=row.size_label, basal_circumference=float(row.circumference),
                    resting_circumference=float(row.resting_circumference) if pd.notna(row.resting_circumference) else None,
                    placement_row=placement_row, cue_duration=float(row.cue_duration), duration_column=column,
                    origin_x=float(column), origin_y=float(PITCH * (3 - placement_row)), glyph_scale=SCALE,
                    requested_initial_arc=float(row.physical_distance_at_onset), requested_angle_deg=float(row.requested_angle_deg),
                    actual_target_angle_in_P_coordinates_deg=target, initial_angular_target_distance_deg=initial_error,
                    final_angular_target_distance_deg=final_error, fraction_initial_angular_error_closed=closed,
                    display_background_clipped=None if closed is None else float(np.clip(closed, 0., 1.)),
                    colour_clipped_below_zero=bool(closed is not None and closed < 0.),
                    colour_clipped_above_one=bool(closed is not None and closed > 1.),
                    initial_overlap=bool(row.initial_overlap), overlap_fraction_cue_support=float(row.overlap_fraction_cue_support),
                    overlap_fraction_original_core=float(row.overlap_fraction_original_core),
                    full_physical_support=float(row.full_physical_support), initial_angular_support_deg=float(row.initial_angular_support_deg),
                    requested_peak=float(row.requested_peak), initial_sampled_cue_peak=float(row.initial_sampled_cue_peak),
                    initial_input_coefficient_integral=float(row.initial_cue_integral),
                    full_field_tail_status=row.full_field_tail_status,
                    field_tail_failed=row.full_field_tail_status == "fails_finite_tail_gate",
                    field_tail_unavailable=row.full_field_tail_status == "unavailable_from_retained_primary_only_history",
                    physical_C_withdrawal=float(row.physical_C_withdrawal), physical_C_final=float(row.physical_C_final),
                    actual_cue_path=row.actual_cue_path, actual_cue_sha256=row.actual_cue_sha256,
                    actual_recovery_path=row.actual_recovery_path, actual_recovery_sha256=row.actual_recovery_sha256)
                for label, vector in (("withdrawal", withdrawal), ("final", final)):
                    strong = abs(vector) >= R_GATE
                    value.update({label + "_P_real": float(vector.real), label + "_P_imag": float(vector.imag),
                        label + "_R": float(abs(vector)), label + "_raw_phase_deg": float(np.rad2deg(np.angle(vector))),
                        label + "_heading_displayed": bool(strong),
                        label + "_arrow_delta_x": float(SCALE * vector.real) if strong else None,
                        label + "_arrow_delta_y": float(SCALE * vector.imag) if strong else None,
                        label + "_magnitude_circle_radius": float(SCALE * abs(vector)) if not strong else None})
                records.append(value)
    glyphs = pd.DataFrame(records)
    assert len(glyphs) == 480 and glyphs.unique_protocol_id.nunique() == 355
    assert glyphs.initial_overlap.sum() == 280
    assert glyphs.field_tail_failed.sum() == 2 and glyphs.field_tail_unavailable.sum() == 17
    # All tips fit the same tile bounds; no magnitude clipping or cell-dependent scaling.
    assert max(glyphs.withdrawal_R.max(), glyphs.final_R.max()) * SCALE < min(.48, PITCH / 2)
    glyphs.to_csv(dest("paired_vector_glyphs.csv", "data"), index=False)
    return glyphs, maximum_error


def draw_vector(ax, row, label, color, width, zorder):
    x, y = row.origin_x, row.origin_y
    if getattr(row, label + "_heading_displayed"):
        dx, dy = getattr(row, label + "_arrow_delta_x"), getattr(row, label + "_arrow_delta_y")
        ax.quiver([x], [y], [dx], [dy], angles="xy", scale_units="xy", scale=1., units="xy", width=width,
                  headwidth=3.2, headlength=3.8, headaxislength=3.5, minlength=0., pivot="tail", color=color, zorder=zorder)
        if label == "withdrawal":
            ax.plot(x + dx, y + dy, marker="o", mfc="none", mec=color, mew=.7, ms=3., zorder=zorder + 2)
    else:
        ax.add_patch(Circle((x, y), getattr(row, label + "_magnitude_circle_radius"),
            facecolor="none" if label == "withdrawal" else color, edgecolor=color, lw=.9, zorder=zorder))


def render(glyphs):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42, "svg.fonttype": "none"})
    cmap = LinearSegmentedColormap.from_list("angular_error_closed", ["#ffffff", "#d9e8e6", "#8dbdb7"])
    norm = Normalize(0, 1)
    transforms = []
    for model in MODELS:
        fig = plt.figure(figsize=(14.2, 10.5))
        fig.text(.035, .968, f"{NAMES[model]} · Stationary cue at fixed angle or fixed initial arc", fontsize=17, weight="bold")
        fig.text(.035, .936, "P = M₁ / M₁(initial): direction and length show the activity vector. Grey rightward guide: initial P = 1.", fontsize=11.5)
        handles = [Line2D([0], [0], color=ORANGE, lw=1.2, marker="o", mfc="none", ms=4, label="At cue withdrawal"),
                   Line2D([0], [0], color=BLACK, lw=1.8, label="After 13,000 cue-free units")]
        fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.028, .923), ncol=2, frameon=False, fontsize=11.5, handlelength=2.)
        cax = fig.add_axes([.50, .891, .43, .014])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", ticks=[0, .5, 1])
        cb.ax.tick_params(labelsize=9, pad=2, length=2)
        cb.set_label("Fraction of initial angular target error closed (colour clipped to 0–1)", fontsize=10, labelpad=4)
        sizes = (1., 2., 3.) if model == "spring" else (20., 30., 40.)
        for ai, arm in enumerate(ARMS):
            fig.text(.035, .821 if ai == 0 else .463, "A · Fixed cue-centre angle" if ai == 0 else "B · Fixed initial cue-centre arc", fontsize=13, weight="bold")
            for si, size in enumerate(sizes):
                ax = fig.add_axes([.105 + .31 * si, .535 if ai == 0 else .177, .245, .252])
                panel = glyphs[(glyphs.model == model) & (glyphs.arm == arm) & (glyphs.size_parameter == size)]
                assert len(panel) == 20
                for row in panel.itertuples():
                    x, y = row.origin_x, row.origin_y
                    face = cmap(norm(row.display_background_clipped)) if row.final_heading_displayed else "#ededed"
                    ax.add_patch(Rectangle((x - .48, y - PITCH / 2 + .01), .96, PITCH - .02,
                        facecolor=face, edgecolor="#d8d8d8", lw=.5, zorder=0))
                    ax.plot([x, x + SCALE], [y, y], color="#aaaaaa", lw=.65, zorder=1)
                    ax.plot(x, y, ".", color="#999999", ms=2., zorder=2)
                    draw_vector(ax, row, "withdrawal", ORANGE, .008, 3)
                    draw_vector(ax, row, "final", BLACK, .014, 4)
                    if not row.final_heading_displayed:
                        ax.text(x, y - .30, f"Rᶠ={row.final_R:.4f}", ha="center", va="center", fontsize=7.3)
                    if row.field_tail_failed:
                        ax.plot(x + .35, y + .34, marker="s", mfc="none", mec="#a55545", mew=1.1, ms=5., zorder=8)
                    if row.field_tail_unavailable:
                        ax.plot(x + .35, y + .34, marker="^", mfc="none", mec="#777777", mew=.85, ms=4.7, zorder=8)
                labels = []
                for placement in range(4):
                    row = panel[panel.placement_row == placement].iloc[0]
                    star = "*" if row.initial_overlap else ""
                    target = degree_label(row.actual_target_angle_in_P_coordinates_deg)
                    labels.append(target + star if ai == 0 else f"d={row.requested_initial_arc:.4g}\n{target}{star}")
                times = (.1, 1., 10., 100., 1000.) if model == "spring" else (1., 10., 100., 300., 1000.)
                ax.set(xlim=(-.5, 4.5), ylim=(-PITCH / 2, 3 * PITCH + PITCH / 2),
                       xticks=range(5), yticks=[PITCH * y for y in (3, 2, 1, 0)])
                ax.set_xticklabels([f"{t:g}" for t in times], fontsize=10)
                ax.set_yticklabels(labels, fontsize=10, linespacing=1.05)
                ax.tick_params(length=0, pad=4)
                ax.set_aspect("equal", adjustable="box")
                for spine in ax.spines.values():
                    spine.set_visible(False)
                title = f"C = {size:g}" if model != "spring" else f"C₀ = {size:g}; basal C = {panel.basal_circumference.iloc[0]:.4f}"
                ax.set_title(title, fontsize=11.5, pad=7, loc="left")
                ax.set_xlabel("Tested cue duration (native units)", fontsize=9.8, labelpad=5)
                transforms.append((model, ai, si, ax))
        width = .2 if model == "spring" else 3.
        peak = dict(wave_pinning=.1, goryachev=6.4, holmes_model3=1., spring=4.)[model]
        widths = [glyphs[(glyphs.model == model) & (glyphs.size_parameter == size)].initial_angular_support_deg.iloc[0] for size in sizes]
        fig.text(.035, .103, f"Native full physical support = {width:g}; prescribed peak = {peak:g}. Initial angular supports (left → right): " + ", ".join(f"{x:.2f}°" for x in widths) + ".", fontsize=10)
        fig.text(.035, .080, "* Cue footprint intersects the original half-contrast core; every overlap is retained. Spring centres stay fixed in material coordinates.", fontsize=10)
        fig.text(.035, .057, "Weak R < 0.1: magnitude-only circles, with no assigned heading. □ failed finite full-field tail check; △ unavailable full-field tail check.", fontsize=10)
        fig.text(.035, .034, "One vector scale throughout: Δglyph = 0.24 P; equal x/y scales. Duration columns are categories. Both arms: 480 rows, 355 unique histories.", fontsize=10)
        fig.text(.035, .012, "Target-relative colour uses wrapped angular distances; it is neither a binary success gate nor a minimum-duration estimate. Native clocks differ between models.", fontsize=9.5)
        fig.canvas.draw()
        for m, a, s, ax in transforms[-6:]:
            point, px, py = ax.transData.transform([[0., 0.], [1., 0.], [0., 1.]])
            assert abs(np.linalg.norm(px - point) - np.linalg.norm(py - point)) < 1e-10
        for extension in ("png", "pdf", "svg"):
            fig.savefig(dest(BASENAMES[model] + "." + extension, "vectors"), dpi=220)
        plt.close(fig)


def documents(glyphs, error):
    generated = datetime.now(timezone.utc).isoformat()
    caption = """Paired stationary-cue experiment, with one review page per model. Each page has two experiment arms and three native sizes. The top row fixes the prescribed cue-centre angle; the bottom row fixes the initial physical arc from the declared preparation axis. Every panel contains all four placements and five tested native durations. Bottom-row labels give initial physical arc d and the actual target angle in normalized P coordinates at that size. The arms coincide at the smallest size; 480 arm-labelled rows represent 355 unique native histories, including 125 exact reused conditions.

Orange vectors show withdrawal P and black vectors show retained P after 13,000 cue-free units. Both are actual physical first moments normalized by the same condition's initial complex moment: P = M₁(t)/M₁(0), with M₁ = C(t) mean[u exp(2πiq)]. Normalization rescales and rotates the reference axes once so P(initial) = 1; it does not use a moving coordinate frame. A pale rightward guide has length one in P units. Every plotted component uses the same display factor 0.24 and equal x/y data-coordinate scales. Weak R = |P| < 0.1 is displayed as an unoriented circle of radius 0.24R, without assigning a heading.

Background colour reports one minus the final wrapped angular distance to the actual target divided by the initial wrapped angular distance to that target. Thus a smaller requested turn is not itself interpreted as weaker retention. Only the displayed colour is clipped to 0–1; raw angular distances, signed P components and the unclipped value are exported. Weak-final-vector cells have no angular colour. This target-proximity measure is not a binary success gate or a refined minimum cue duration.

Asterisks mark initial cue-footprint/original-half-contrast-core overlap; all 280 arm-labelled overlap conditions remain. A square marks a failed finite full-field tail check and a triangle marks unavailable full-field tail qualification, separately from vector magnitude. Fixed physical cue support and prescribed peak retain native values; initial angular footprint consequently changes across size. Spring uses resting C₀ and actual basal C, a fixed material cue centre and a fixed instantaneous physical support. Its prescribed physical arc subsequently scales with C(t). These figures characterize local activity and its resultant, without a cell movement law, site ancestry or an isolated causal geometry attribution. Duration columns are categories; native clocks are not calibrated across models. The old angle-arm figures and review packet are unchanged.
"""
    dest("CAPTION.md").write_text(caption)
    manifest = dict(generated_at=generated, role="Paired-arm review graphics; no native integrations or new outcome thresholds",
        displayed_arm_rows=480, unique_protocol_histories=355, exact_shared_conditions=125,
        models=list(MODELS), pages=[BASENAMES[m] + ".png" for m in MODELS],
        source_sha256=SOURCES, renderer=str(Path(__file__)), renderer_sha256=sha(__file__),
        glyph_table=str(dest("paired_vector_glyphs.csv", "data")), glyph_table_sha256=sha(dest("paired_vector_glyphs.csv", "data")),
        mapping=dict(arm_rows=list(ARMS), size_columns=dict(chemical=[20, 30, 40], spring_resting=[1, 2, 3]),
            placement_rows="Descending target angles in angle arm; descending physical arcs in fixed-arc arm", duration_columns="Native tested durations as categorical indices0..4",
            origin_x="duration column index", origin_y="0.95*(3-placement_row)", global_P_display_scale=SCALE,
            axes_aspect="equal", axis_limits=dict(x=[-.5, 4.5], y=[-PITCH / 2, 3 * PITCH + PITCH / 2]),
            background="1-abs(wrap(actual_final_phase-actual_target))/abs(wrap(-actual_target)); finalR>=.1 only; displaycolour clipped0..1",
            withdrawal_style="thin orange vector with open tip", final_style="black vector", weak_style="unoriented circle of radius0.24R"),
        quiver=dict(angles="xy", scale_units="xy", scale=1., units="xy", withdrawal_width=.008, final_width=.014,
            headwidth=3.2, headlength=3.8, headaxislength=3.5, minlength=0., pivot="tail"),
        withdrawal_weak=int((~glyphs.withdrawal_heading_displayed).sum()), final_weak=int((~glyphs.final_heading_displayed).sum()),
        field_tail_failures=int(glyphs.field_tail_failed.sum()), field_tail_unavailable=int(glyphs.field_tail_unavailable.sum()),
        initial_overlap_rows=int(glyphs.initial_overlap.sum()), maximum_R_replay_error=float(error),
        colour_raw_below_zero=int(glyphs.colour_clipped_below_zero.sum()), colour_raw_above_one=int(glyphs.colour_clipped_above_one.sum()),
        old_angle_outputs_preserved=True)
    write("renderer_manifest.json", manifest)
    links = "\n".join(f"- [{NAMES[m]}]({BASENAMES[m]}.png) · [PDF](supporting/vectors/{BASENAMES[m]}.pdf) · [SVG](supporting/vectors/{BASENAMES[m]}.svg)" for m in MODELS)
    dest("README.md").write_text("# Paired stationary-cue review\n\nTwo arms of the same stationary-cue experiment, shown in four model pages. These are new review graphics; the old angle-only art and packet are preserved.\n\n" + links + "\n\n- [Exact plotted vector, target, colour and flag values](supporting/data/paired_vector_glyphs.csv)\n- [Caption and definitions](supporting/metadata/CAPTION.md)\n- [Source and renderer manifest](supporting/metadata/renderer_manifest.json)\n- [Source hashes](supporting/metadata/source_sha256.json)\n\nEvery cell pairs cue withdrawal (orange) with 13,000-unit cue-free retention (black). The global vector scale is shared across every model, size and arm. Target-relative backgrounds use wrapped angular distances. All overlaps, weak-vector magnitude displays and independent field-tail statuses remain visible. No native simulation, threshold search or new mechanism attribution is introduced.\n")
    SOURCES[str(Path(__file__))] = sha(__file__)
    write("source_sha256.json", SOURCES)
    artifacts = {str(path.relative_to(OUT)): sha(path) for path in OUT.rglob("*") if path.is_file() and path.name != "DELIVERY_CHECK.json"}
    write("DELIVERY_CHECK.json", dict(generated_at=generated, source_sha256=SOURCES, artifact_sha256=artifacts,
        png_pages=4, vector_pages=8, glyph_rows=480, unique_histories=355, all_sources_verified=True,
        all_pixels_use_equal_xy_coordinates=True, all_vector_lengths_unclipped=True))
    print(json.dumps(dict(passed=True, pages=4, displayed=480, unique_histories=355, output=str(OUT)), indent=2))


if __name__ == "__main__":
    values, error = load_values()
    render(values)
    documents(values, error)
