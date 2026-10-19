"""Three compact native-history contrasts from the accepted fixed-initial-arc arm."""
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import revision_standard_angles_report_20261019 as prior
import ring_distance_prototype_20261017 as ring

ROOT = Path(__file__).resolve().parents[2]
GROUP = ROOT / "figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3"
OUT = GROUP / "arc_illustrative_histories"
DATA = GROUP / "standard_arclength_sizes/supporting/data"
CARDS = (
    ("holmes_fixed_arc_size", "Holmes · Size changes cue-on reorientation", (132, 172), "C = 20", "C = 40"),
    ("wave_pinning_arc_near_far", "Wave-pinning · Placement changes cue-on reorientation", (43, 58), "d = 1.25", "d = 10"),
    ("goryachev_arc_exception", "Goryachev · Near and far cues follow different time courses", (113, 118), "d = 5", "d = 10"),
)
COLORS = ("#236f92", "#a4486a")
SOURCE = {}


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def record(path):
    path = Path(path); SOURCE[str(path)] = sha(path); return path
def dest(name, group="data"):
    path = OUT / name if name.endswith(".png") or name == "README.md" else OUT / "supporting" / group / name
    path.parent.mkdir(parents=True, exist_ok=True); return path
def write(name, value): dest(name, "metadata").write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
def wrap(x): return np.degrees(np.angle(np.exp(1j * np.radians(x))))


def load():
    tablepath = record(DATA / "stationary_arclength_summary.csv")
    assert sha(tablepath) == "9b00954df3a3ff50f65de37ea9f3d82968bab10196b8cdc10a5d054d5329e278"
    frame = pd.read_csv(tablepath); assert len(frame) == 240 and frame.status.eq("complete").all()
    glyphpath = record(GROUP / "paired_stationary_review/supporting/data/paired_vector_glyphs.csv")
    glyphs = pd.read_csv(glyphpath); a = glyphs[glyphs.arm == "common_initial_physical_arc"].copy()
    a["withdrawal_error_closed"] = np.where(a.withdrawal_heading_displayed,
        1 - np.abs(wrap(a.withdrawal_raw_phase_deg - a.actual_target_angle_in_P_coordinates_deg)) / a.initial_angular_target_distance_deg, np.nan)
    a["retained_error_closed"] = a.fraction_initial_angular_error_closed
    a.to_csv(dest("arc_target_closure240.csv"), index=False)
    adjacent = []
    for model, part in a.groupby("model"):
        for metric in ("withdrawal_error_closed", "retained_error_closed"):
            for label, keys, order in (("increasing_size", ["requested_initial_arc", "cue_duration"], "size_parameter"),
                                       ("increasing_initial_arc", ["size_parameter", "cue_duration"], "requested_initial_arc")):
                for _, group in part.groupby(keys):
                    group = group.sort_values(order)
                    for j in range(len(group) - 1):
                        left, right = group.iloc[j], group.iloc[j + 1]
                        adjacent.append(dict(model=model, metric=metric, comparison=label,
                            condition_id_first=int(left.condition_id), condition_id_second=int(right.condition_id),
                            size_first=float(left.size_parameter), size_second=float(right.size_parameter),
                            initial_arc_first=float(left.requested_initial_arc), initial_arc_second=float(right.requested_initial_arc),
                            cue_duration=float(left.cue_duration), closure_first=float(left[metric]), closure_second=float(right[metric]),
                            closure_change=float(right[metric] - left[metric])))
    pd.DataFrame(adjacent).to_csv(dest("arc_adjacent_comparisons680.csv"), index=False); assert len(adjacent) == 680
    traces = np.load(record(DATA / "stationary_arclength_traces.npz"))
    cases, profiles, curves, endpoints, rates = {}, [], [], [], []
    for _, _, ids, _, _ in CARDS:
        for cid in ids:
            row = frame[frame.condition_id == cid].iloc[0]
            assert row.full_field_tail_status == "passes_finite_tail_gate"
            parts = []
            for phase in ("cue", "recovery"):
                path = record(row["actual_" + phase + "_path"]); assert sha(path) == row["actual_" + phase + "_sha256"]
                origin = next(parent for parent in path.parents if (parent / "manifest.json").exists())
                kernels = [Path(ring.__file__), ROOT / "code/core/models" / (row.model + ".py")]
                kernels += list((ROOT / "code/core/src/polarity1d").glob("*.py"))
                for kernel in kernels:
                    assert sha(kernel) == sha(origin / "frozen_sources" / kernel.relative_to(ROOT))
                    record(kernel)
                with np.load(path) as z:
                    assert "states" in z, "Selected examples require full native fields"
                    parts.append(dict(time=z["time"].copy(), states=z["states"].copy(), x=z["x"].copy(), stimulus=z["stimulus"].copy()))
            np.testing.assert_array_equal(parts[0]["states"][-1], parts[1]["states"][0])
            C, T = row.circumference, row.cue_duration
            t = np.r_[parts[0]["time"], T + parts[1]["time"][1:]]
            u = np.vstack([parts[0]["states"][:, 0], parts[1]["states"][1:, 0]])
            q = parts[0]["x"] / C
            M1 = C * np.mean(u * np.exp(2j * np.pi * q), axis=1); P = M1 / M1[0]
            np.testing.assert_allclose(P.real, traces[f"c{cid}_P_real"], rtol=0, atol=2e-14)
            np.testing.assert_allclose(P.imag, traces[f"c{cid}_P_imag"], rtol=0, atol=2e-14)
            signed = ((q - row.preparation_axis_fraction + .5) % 1 - .5) * C
            order = np.argsort(signed); x = signed[order]
            off = len(parts[0]["time"]) - 1
            angle = np.where(abs(P) >= .1, wrap(np.degrees(np.angle(P))), np.nan)
            model, grid = ring.native_model(row.model), ring.RingGrid(C, len(q))
            weights = np.exp(2j * np.pi * q)
            omega = []
            for j, state in enumerate(parts[0]["states"]):
                derivative = ring.native_rhs(model, grid, state, parts[0]["stimulus"])
                dM1 = C * np.mean(derivative[0] * weights)
                strong = abs(P[j]) >= .1
                value = float(np.degrees(np.imag(np.conj(M1[j]) * dM1) / abs(M1[j])**2)) if strong else np.nan
                omega.append(value)
                rates.append(dict(condition_id=cid, time_from_cue_onset=float(parts[0]["time"][j]),
                    actual_target_deg=row.target_angle_in_P_coordinates_deg, angle_deg=float(angle[j]), R=float(abs(P[j])),
                    native_RHS_angular_rate_deg_per_native_time=value,
                    native_RHS_rate_per_actual_target_per_native_time=float(value / abs(wrap(row.target_angle_in_P_coordinates_deg))),
                    heading_and_rate_available=bool(strong)))
            omega = np.asarray(omega)
            profile = dict(initial=u[0, order], withdrawal=u[off, order], final=u[-1, order])
            for phase, values in profile.items():
                for j in range(len(x)):
                    profiles.append(dict(condition_id=cid, model=row.model, size=row.size_parameter, phase=phase,
                        signed_physical_position=float(x[j]), native_active=float(values[j]), cue_coefficient=float(parts[0]["stimulus"][order[j]])))
            for j in range(len(t)):
                curves.append(dict(condition_id=cid, time_from_onset=float(t[j]), time_from_withdrawal=float(t[j] - T),
                    phase="cue_on" if j <= off else "cue_free", P_real=float(P[j].real), P_imag=float(P[j].imag),
                    R=float(abs(P[j])), angle_deg=float(angle[j]), M0=float(C * u[j].mean())))
            counts = {key: prior.periodic_counts(value) for key, value in profile.items()}
            metric = a[a.condition_id == cid].iloc[0]
            endpoint = dict(condition_id=cid, model=row.model, size=row.size_parameter, cue_duration=T,
                initial_arc=row.physical_distance_at_onset, actual_target_deg=row.target_angle_in_P_coordinates_deg,
                initial_core_width=row.original_core_physical_width, overlap_fraction_cue_support=row.overlap_fraction_cue_support,
                withdrawal_angle_deg=row.angle_withdrawal_deg, final_angle_deg=row.angle_final_deg,
                withdrawal_R=row.R_withdrawal, final_R=row.R_final,
                cue_on_angular_secant_deg_per_native_time=float(angle[off] / T),
                cue_on_target_normalized_secant_per_native_time=float(angle[off] / abs(wrap(row.target_angle_in_P_coordinates_deg)) / T),
                maximum_saved_native_RHS_angular_rate=float(np.nanmax(omega)),
                minimum_saved_native_RHS_angular_rate=float(np.nanmin(omega)),
                withdrawal_error_closed=metric.withdrawal_error_closed, retained_error_closed=metric.retained_error_closed,
                initial_peaks=counts["initial"][0], withdrawal_peaks=counts["withdrawal"][0], final_peaks=counts["final"][0],
                full_field_tail_excursion=row.final_full_field_tail_excursion, full_field_tail_status=row.full_field_tail_status,
                cue_peak=float(parts[0]["stimulus"].max()), input_coefficient_integral=float(C * parts[0]["stimulus"].mean()))
            endpoints.append(endpoint)
            cases[cid] = dict(row=row, x=x, profile=profile, time=t - T, cue_time=parts[0]["time"],
                omega=omega, P=P, angle=angle, off=off, endpoint=endpoint)
    traces.close()
    pd.DataFrame(profiles).to_csv(dest("selected_native_profiles.csv"), index=False)
    pd.DataFrame(curves).to_csv(dest("selected_vector_histories.csv"), index=False)
    pd.DataFrame(endpoints).to_csv(dest("selected_conditions.csv"), index=False)
    pd.DataFrame(rates).to_csv(dest("selected_cue_on_native_RHS_rates.csv"), index=False)
    record(ring.__file__)
    for name in ("wave_pinning", "holmes_model3", "goryachev"):
        record(ROOT / "code/core/models" / (name + ".py"))
    return cases


def render(cases):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "pdf.fonttype": 42, "svg.fonttype": "none"})
    for filename, title, ids, first, second in CARDS:
        pair = [cases[cid] for cid in ids]
        assert abs(pair[0]["endpoint"]["input_coefficient_integral"] - pair[1]["endpoint"]["input_coefficient_integral"]) < 1e-12
        fig = plt.figure(figsize=(12.2, 10.6))
        fig.text(.045, .958, title, fontsize=18, weight="bold")
        T = pair[0]["row"].cue_duration
        fig.text(.045, .923, f"Cue-on response for {T:g} native units. Rate is measured from the native RHS at each saved state; no reintegration.", fontsize=11.3)
        angle_ax = fig.add_axes([.085, .620, .385, .245]); rate_ax = fig.add_axes([.565, .620, .385, .245])
        R_ax = fig.add_axes([.085, .420, .865, .115])
        for j, (case, label) in enumerate(zip(pair, (first, second))):
            row, tt = case["row"], case["cue_time"]
            angle_ax.plot(tt, case["angle"][:case["off"] + 1], color=COLORS[j], lw=1.8, label=label)
            angle_ax.axhline(row.target_angle_in_P_coordinates_deg, color=COLORS[j], ls=":", lw=.9)
            rate_ax.plot(tt, case["omega"], color=COLORS[j], lw=1.8, label=label)
            R_ax.plot(tt, abs(case["P"][:case["off"] + 1]), color=COLORS[j], lw=1.8, label=label)
        ticks = sorted(set([0., 1., 10., 100., float(T)]))
        for ax in (angle_ax, rate_ax, R_ax):
            ax.set_xscale("symlog", linthresh=1.)
            ax.set_xlim(0, T)
            ax.set_xticks(ticks, [f"{t:g}" for t in ticks]); ax.minorticks_off()
            ax.set_xlabel("Time since cue onset (native units)", fontsize=10)
            ax.spines[["top", "right"]].set_visible(False)
        angle_ax.set_ylabel("Activity-vector angle (degrees)")
        angle_ax.set_title("Cue-on direction; dotted lines mark targets", fontsize=11.5, loc="left")
        rate_ax.set_ylabel("Angular rate (degrees / native unit)")
        rate_ax.set_title("Instantaneous cue-on angular rate", fontsize=11.5, loc="left")
        rate_ax.axhline(0, color="#bbbbbb", lw=.6)
        R_ax.set_ylabel("R = |P|")
        R_ax.set_title("Cue-on vector magnitude", fontsize=10.5, loc="left")
        R_ax.axhline(1., color="#aaaaaa", ls=":", lw=.8)
        angle_ax.legend(loc="best", frameon=False, fontsize=10)
        axes = [fig.add_axes([.085 + .48 * j, .190, .385, .140]) for j in range(2)]
        ymax = max(np.max(v) for c in pair for v in c["profile"].values()) * 1.08
        halfspan = max(c["row"].circumference for c in pair) / 2
        for j, (ax, case, label) in enumerate(zip(axes, pair, (first, second))):
            row = case["row"]; d = row.physical_distance_at_onset
            ax.axvspan(row.original_core_left_physical, row.original_core_right_physical, color="#999999", alpha=.12, zorder=0)
            ax.axvspan(d - 1.5, d + 1.5, color="#d78637", alpha=.17, zorder=0)
            for phase, color, ls in (("initial", "#777777", "--"), ("withdrawal", "#d77b29", "-"), ("final", "#111111", "-")):
                ax.plot(case["x"], case["profile"][phase], color=color, ls=ls, lw=1.65)
            ax.set(xlim=(-halfspan, halfspan), ylim=(0, ymax), xlabel="Physical arc from original axis", ylabel="Native active")
            ax.set_title(f"C = {row.size_parameter:g} · d = {d:g}; target {row.target_angle_in_P_coordinates_deg:g}°", loc="left", color=COLORS[j], fontsize=10.8, weight="bold")
            if row.circumference / 2 < halfspan:
                for left, right in ((-halfspan, -row.circumference / 2), (row.circumference / 2, halfspan)):
                    ax.axvspan(left, right, color="#eeeeee", zorder=-1)
                ax.text(-.75 * halfspan, .8 * ymax, "outside\nC = 20", ha="center", fontsize=7.8, color="#666666")
                ax.text(.75 * halfspan, .8 * ymax, "outside\nC = 20", ha="center", fontsize=7.8, color="#666666")
            fig.text(.085 + .48 * j, .112,
                f"Core width {row.original_core_physical_width:.3f}; cue overlap {100 * row.overlap_fraction_cue_support:.1f}%.\n"
                f"Withdrawal {row.angle_withdrawal_deg:.2f}°; after 13,000 cue-free units {row.angle_final_deg:.2f}°.", fontsize=9.7)
        fig.legend(handles=[Line2D([0], [0], color="#777777", ls="--", label="Initial"),
            Line2D([0], [0], color="#d77b29", label="At withdrawal"), Line2D([0], [0], color="#111111", label="After 13,000 cue-free units"),
            Patch(facecolor="#dddddd", label="Original half-contrast core"), Patch(facecolor="#f3d8b8", label="Cue footprint")],
            loc="upper left", bbox_to_anchor=(.055, .370), ncol=5, frameon=False, fontsize=8.6, handlelength=1.5, columnspacing=.9)
        for ax in axes:
            ax.spines[["top", "right"]].set_visible(False)
        peak = pair[0]["endpoint"]["cue_peak"]
        fig.text(.045, .072, f"Native cue: full physical support 3, prescribed peak {peak:g}; matched input-coefficient integral {pair[0]['endpoint']['input_coefficient_integral']:g}, not activated amount.", fontsize=9.8)
        fig.text(.045, .049, "Cue-on time axes: linear from 0–1, logarithmic thereafter. Rates are angular activity-vector rates; they are not peak travel or cell velocity.", fontsize=9.8)
        fig.text(.045, .026, "Weak R < 0.1 would mask heading and rate independently; no weak interval occurs here. All selected histories pass the finite full-field tail check.", fontsize=9.8)
        for extension in ("png", "pdf", "svg"):
            fig.savefig(dest(filename + "." + extension, "vectors"), dpi=210)
        plt.close(fig)
        print(str(OUT / (filename + ".png")), flush=True)


def docs():
    record(__file__)
    write("source_sha256.json", SOURCE)
    write("selection_and_validation.json", dict(generated_at=datetime.now(timezone.utc).isoformat(),
        selection=[dict(basename=name, title=title, fixed_arc_condition_ids=list(ids)) for name, title, ids, _, _ in CARDS],
        selected_conditions=6, selected_native_archives=12, full_fields_required=True, all_selected_full_tail_pass=True,
        rationale="Cue-on native angular trajectories and instantaneous rates for a same-arc size contrast, a near/far placement contrast and Goryachev's distinct time courses. Retained endpoints are secondary profile context. No unique causal mechanism inferred.",
        instantaneous_rate="degrees(Im(conj(M1)*native_RHS_dM1)/abs(M1)^2), measured only whenR>=.1; no differentiation across a weak-vector gap",
        angular_and_target_normalized_rates_exported=True,
        branch_safe_metric="1-abs(wrap(phase-target))/abs(wrap(-target)); independently at withdrawal and finite retained endpoint, withR>=.1",
        all_grid_readout_rows=240, adjacent_comparison_records=680, comparison_records_are_not_independent_replicates=True,
        no_new_integration=True, previous_paired_pdf_unchanged=True))
    dest("CAPTION.md", "metadata").write_text(r"""# Cue-on reorientation illustrations

Three selected contrasts show the activity-vector heading, instantaneous angular rate and magnitude during the prescribed stationary cue. Local native profiles at onset, withdrawal and after the finite 13,000-unit continuation supply secondary spatial context. The latter endpoints are not rate measurements. Each panel retains the same native physical-distance and concentration scale within its contrast. The fixed-arc size comparison also changes the qualified basal organization and target normal angle.

The physical first moment and its instantaneous phase rate are

\[
M_1(t)=C(t)\,\operatorname{mean}_q\!\left[u(q,t)e^{2\pi i q}\right],
\qquad
\omega(t)=\frac{\operatorname{Im}\!\left[\overline{M_1(t)}\,\dot M_1(t)\right]}{|M_1(t)|^2}.
\tag{1}
\]

The displayed rates convert radians to degrees and use each exact saved native full state with the unchanged native right-hand side. The cue-on value at time zero is the right derivative with the cue applied. No integration, phase interpolation across a weak interval or finite-pulse endpoint substitution is used. Heading and rate are independently unavailable when the normalized vector magnitude is below 0.1; none of these selected cue-on frames is weak. All selected targets are non-antipodal. Native time axes are linear from zero to one and logarithmic thereafter.

Absolute angular rate and angular rate divided by the actual requested target angle are different observables. Both are exported, together with explicitly named onset-to-withdrawal secants and actual native time samples. In particular, the selected Goryachev farther cue has a larger absolute angular secant but a smaller fraction of its larger requested turn per native unit. The curves do not establish peak transport, cell velocity or a unique causal role of geometry. Full-field tail qualification applies to the secondary finite retained endpoint and passes for all six selected histories. Equal prescribed input-coefficient exposure does not imply equal activated biochemical amount. Core and cue-footprint overlap describe fixed spatial regions rather than site ancestry.
""")
    links = "\n".join(f"- [{title}]({name}.png) · [PDF](supporting/vectors/{name}.pdf)" for name, title, _, _, _ in CARDS)
    dest("README.md").write_text("# Fixed-arc cue-on rate illustrations\n\nThree selected review contrasts foreground actual cue-on heading, instantaneous native angular rate and vector magnitude. Initial, withdrawal and retained profiles supply secondary context. The four-page paired-grid review remains unchanged.\n\n" + links + "\n\n- [Exact selected condition values and named secants](supporting/data/selected_conditions.csv)\n- [Cue-on rates from native full-state right-hand sides](supporting/data/selected_cue_on_native_RHS_rates.csv), including absolute and actual-target-normalized rates\n- [Actual native profiles](supporting/data/selected_native_profiles.csv) and [physical-vector histories](supporting/data/selected_vector_histories.csv)\n- [Caption, equation and limitations](supporting/metadata/CAPTION.md)\n- [All 240 branch-safe cue-withdrawal/retained values](supporting/data/arc_target_closure240.csv) and [680 adjacent descriptive comparisons](supporting/data/arc_adjacent_comparisons680.csv)\n- [Selection and validation](supporting/metadata/selection_and_validation.json) · [Source hashes](supporting/metadata/source_sha256.json)\n\nThe qualitative size and placement trends are not pointwise universal. Angular rate and target-relative rate differ, and finite retained outcomes do not measure cue-on rate. The cards show native local activity and its first-moment resultant, without interpreting fixed regions as site ancestry or prescribing cell movement. Common initial physical cue offsets across size also change target normal angles and basal organization. Every selected condition retains the same native kinetics, cue channel, width, amplitude and finite release horizon; input-coefficient exposure is not equal activated biochemical amount.\n")
    artifacts = {str(path.relative_to(OUT)): sha(path) for path in OUT.rglob("*") if path.is_file() and path.name != "DELIVERY_CHECK.json"}
    write("DELIVERY_CHECK.json", dict(source_sha256=SOURCE, artifact_sha256=artifacts))


if __name__ == "__main__":
    cases = load(); render(cases); docs()
