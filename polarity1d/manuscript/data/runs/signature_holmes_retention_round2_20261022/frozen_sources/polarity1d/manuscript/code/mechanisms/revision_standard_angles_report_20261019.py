"""Replay the standardized-angle checkpoint, preserving every initial overlap."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import argparse
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import revision_standard_angles_duration_20261019 as native

ROOT = native.ROOT
REPO = ROOT.parents[1]
RUN = native.DEFAULT
OUT = native.TOPIC
NAMES = {"wave_pinning": "Wave-pinning", "goryachev": "Goryachev", "holmes_model3": "Holmes", "spring": "Spring"}
COLORS = {"pi": "#444444", "pi_2": "#4477AA", "pi_4": "#228877", "pi_8": "#CC7744"}
LABELS = {"pi": "π = 180°", "pi_2": "π/2 = 90°", "pi_4": "π/4 = 45°", "pi_8": "π/8 = 22.5°"}
SOURCE = {}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def record(path):
    path = Path(path)
    SOURCE[str(path)] = sha(path)
    return path


def dest(name, group="data"):
    folder = OUT if name.endswith(".png") else OUT / "supporting" / group
    folder.mkdir(parents=True, exist_ok=True)
    return folder / name


def write_json(name, value):
    dest(name, "metadata").write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def savefig(stem, fig):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(dest(stem + "." + ext, "vectors"), dpi=220, bbox_inches="tight")
    plt.close(fig)


def load_phase(path, spring, Cb):
    record(path)
    with np.load(path) as z:
        t = z["time"].copy()
        full = z["states"].copy() if "states" in z else None
        if spring:
            n = (full.shape[1] - 1) // 2
            u = full[:, :n]
            C = full[:, -1]
            q = z["material_positions"].copy()
        else:
            u = full[:, 0] if full is not None else z["primary"].copy()
            C = np.full(len(t), Cb)
            q = z["x"].copy() / Cb
    return dict(time=t, u=u, C=C, full=full, q=q)


def heading(P):
    angle = np.rad2deg(np.angle(P))
    return np.where(np.isclose(np.abs(angle), 180, rtol=0, atol=1e-6), 180., angle)


def periodic_counts(profile):
    n = len(profile)
    contrast = np.ptp(profile)
    if contrast <= 1e-12:
        return 0, 0
    peaks, _ = find_peaks(np.tile(profile, 3), prominence=.05 * contrast)
    count = np.count_nonzero((peaks >= n) & (peaks < 2*n))
    inside = profile >= profile.min() + .5 * contrast
    regions = np.count_nonzero(inside & ~np.roll(inside, 1))
    return int(count), int(regions)


def bracket(t, mask):
    found = np.flatnonzero(mask)
    if not len(found):
        return dict(status="right_censored", lower=float(t[-1]), upper=None)
    j = int(found[0])
    return dict(status="already_true_at_start" if j == 0 else "sample_bracket", lower=float(t[j-1]) if j else 0., upper=float(t[j]))


def tail_error(phase, scales, window=1000.):
    keep = phase["time"] >= phase["time"][-1] - window
    if phase["full"] is None:
        return None
    y = phase["full"][keep]
    # Empirical complete-field excursion normalized by frozen species scales.
    return float(np.max(np.ptp(y, axis=0) / np.asarray(scales)))


def main():
    manifest = json.loads(record(RUN / "manifest.json").read_text())
    inventory = pd.read_csv(record(RUN / "results/condition_inventory.csv"))
    assert len(inventory) == 80 and (inventory.status == "complete").all()
    for rel, digest in manifest["source_sha256"].items():
        assert sha(RUN / "frozen_sources" / rel) == digest
    for rel, digest in manifest["input_sha256"].items():
        assert sha(RUN / rel) == digest
    rows, events, profile_rows, context_rows, audits = [], [], [], [], []
    trajectories = {}
    for case in manifest["cases"]:
        cid = case["condition_id"]
        found = inventory[inventory.condition_id == cid].iloc[0]
        spring = case["model"] == "spring"
        with np.load(RUN / case["input"]) as z:
            initial = z["state"].copy()
        with np.load(RUN / case["cue_input"]) as z:
            q = z["material_positions"].copy()
            signal = z["initial_cue"].copy()
        n = len(q)
        Cb = case["circumference"]
        u0 = initial[:n] if spring else initial[0]
        core = native.core_info(q, u0, case["preparation_axis_fraction"])
        w = core["weights"]
        background_amount = Cb * core["background"] * w.sum()
        initial_excess = Cb * (u0 @ w) - background_amount
        parts = []
        for label in ("cue", "recovery"):
            path = Path(found["actual_" + label + "_path"])
            assert sha(path) == found["actual_" + label + "_sha256"]
            parts.append(load_phase(path, spring, Cb))
            np.testing.assert_array_equal(parts[-1]["q"],q)
            with np.load(path) as z:
                if spring:
                    key="cue_centre_fraction" if "cue_centre_fraction" in z else "cue_center_fraction"
                    assert float(z[key])==case["cue_centre_fraction"]
                    assert float(z["cue_peak"])==(case["requested_peak"] if label=="cue" else 0.)
                else:
                    np.testing.assert_array_equal(z["stimulus"],signal if label=="cue" else np.zeros(n))
        np.testing.assert_array_equal(native.full_state(found.actual_cue_path, 0.), initial)
        np.testing.assert_array_equal(native.full_state(found.actual_cue_path, case["cue_duration"]), native.full_state(found.actual_recovery_path, 0.))
        assert parts[0]["time"][-1] == case["cue_duration"] and parts[1]["time"][-1] == 13000.
        t = np.r_[parts[0]["time"], parts[1]["time"][1:] + case["cue_duration"]]
        u = np.vstack((parts[0]["u"], parts[1]["u"][1:]))
        C = np.r_[parts[0]["C"], parts[1]["C"][1:]]
        M0 = C * u.mean(axis=1)
        M1 = C * (u @ np.exp(2j*np.pi*q)) / n
        M2 = C * (u @ np.exp(4j*np.pi*q)) / n
        P = M1 / M1[0]
        R = np.abs(P)
        angle = heading(P)
        target = case["target_angle_in_P_coordinates_deg"]
        error = np.rad2deg(np.angle(np.exp(1j*np.deg2rad(angle-target))))
        amount = C * (u @ w)
        E = (amount - background_amount) / initial_excess
        current_bg_E = C * ((u-core["background"]) @ w) / initial_excess
        phaseindex = len(parts[0]["time"]) - 1
        final = native.full_state(found.actual_recovery_path, 13000.)
        final_u = final[:n] if spring else final[0]
        assert np.array_equal(final_u, u[-1])
        peak_count, region_count = periodic_counts(final_u)
        shift = float(angle[-1]) / 360 * n
        factor = np.exp(-2j*np.pi*np.fft.fftfreq(n)*shift)
        if spring:
            shiftbase = np.r_[np.fft.ifft(np.fft.fft(initial[:n])*factor).real,
                              np.fft.ifft(np.fft.fft(initial[n:-1])*factor).real, initial[-1]]
            rotated_error = np.abs(final-shiftbase)
            shape_error = float(np.max(rotated_error[:-1]))
            relative_shape = float(np.max(rotated_error[:-1] / np.asarray(case["fixed_scales"])[:-1]))
        else:
            shiftbase = np.fft.ifft(np.fft.fft(initial, axis=1)*factor[None, :], axis=1).real
            rotated_error = np.abs(final-shiftbase)
            shape_error = float(rotated_error.max())
            relative_shape = float(np.max(rotated_error / np.asarray(case["fixed_scales"])[:, None]))
        primary_shift = np.fft.ifft(np.fft.fft(u0)*factor).real
        primary_shift_rel = float(np.max(np.abs(final_u-primary_shift))/np.ptp(u0))
        scales = case["fixed_scales"] if spring else np.asarray(case["fixed_scales"])[:, None]
        full_tail = tail_error(parts[1], scales)
        primary_tail = float(np.max(np.ptp(parts[1]["u"][parts[1]["time"] >= 12000], axis=0)) / np.ptp(u0))
        pool_error = 0.
        for part in parts:
            if part["full"] is None:
                continue
            if spring:
                pool = part["C"] * (part["full"][:, :n]+part["full"][:, n:-1]).mean(axis=1)
                pool_error = max(pool_error, float(np.max(np.abs(pool/(Cb*(initial[:n]+initial[n:-1]).mean())-1))))
            else:
                model = native.ring.native_model(case["model"])
                for group in model.conserved_groups:
                    initial_pool = Cb * initial[list(group)].sum()/n
                    pool = Cb * part["full"][:, list(group)].sum(axis=(1, 2))/n
                    pool_error = max(pool_error, float(np.max(np.abs(pool/initial_pool-1))))
        assert pool_error < 1e-7 and np.isfinite(P).all()
        meta = {k:v for k,v in case.items() if k not in ("reuse", "fixed_scales")}
        row = dict(**meta, status="complete", source_kind=case["execution"],
            **{key:found[key] for key in found.index if key.startswith("actual_")},
            angle_final_raw_deg=float(angle[-1]), angle_final_deg=float(angle[-1]) if R[-1]>=.1 else None,
            R_final=float(R[-1]), q1_final=float(abs(M1[-1])/M0[-1]), q2_final=float(abs(M2[-1])/M0[-1]),
            target_error_final_deg=float(error[-1]) if R[-1]>=.1 else None,
            fraction_requested_turn_final=float(angle[-1]/case["requested_angle_deg"]) if R[-1]>=.1 else None,
            angle_withdrawal_deg=float(angle[phaseindex]) if R[phaseindex]>=.1 else None,
            R_withdrawal=float(R[phaseindex]), fixed_region_excess_withdrawal=float(E[phaseindex]),
            fixed_region_excess_final=float(E[-1]), fixed_region_excess_minimum=float(E.min()),
            fixed_region_excess_minimum_time=float(t[E.argmin()]),
            fixed_region_physical_active_amount_final=float(amount[-1]),
            fixed_region_current_background_excess_final=float(current_bg_E[-1]),
            final_resolved_peaks=peak_count, final_halfcontrast_regions=region_count,
            final_primary_rotated_basal_relative_contrast_error=primary_shift_rel,
            final_full_rotated_basal_absolute_error=shape_error,
            final_full_rotated_basal_scale_error=relative_shape,
            final_full_field_tail_excursion=full_tail,
            full_field_tail_gate=None if full_tail is None else bool(full_tail<=1e-5),
            full_field_tail_status="unavailable_from_retained_primary_only_history" if full_tail is None else ("passes_finite_tail_gate" if full_tail<=1e-5 else "fails_finite_tail_gate"),
            final_primary_tail_relative_contrast_excursion=primary_tail,
            physical_C_min=float(C.min()), physical_C_max=float(C.max()), physical_C_final=float(C[-1]),
            conserved_pool_max_relative_error=pool_error)
        rows.append(row)
        for threshold in (.9, .5, .1):
            for phase, keep in (("cue_on", t<=case["cue_duration"]), ("all_saved", np.ones(len(t), bool))):
                events.append(dict(condition_id=cid, model=case["model"], angle_label=case["angle_label"], cue_duration=case["cue_duration"],
                    phase=phase, observable="fixed_original_material_region_transported_background_excess", threshold=threshold,
                    **bracket(t[keep], E[keep]<=threshold), rebound_after_first=bool(np.any(E[np.flatnonzero(E<=threshold)[0]:]>threshold)) if np.any(E<=threshold) else False))
        aq = (360*(q-case["preparation_axis_fraction"])+120)%360-120
        if case["cue_duration"] == (0.1 if spring else 1.):
            for j in range(n):
                context_rows.append(dict(model=case["model"], angle_label=case["angle_label"], material_position=q[j], normal_angle_from_preparation_deg=aq[j],
                    basal_active=u0[j], basal_active_normalized_contrast=(u0[j]-u0.min())/np.ptp(u0), initial_cue=signal[j], initial_cue_unit_peak=signal[j]/case["requested_peak"],
                    fixed_original_core_cell_material_weight=w[j]))
        for j in range(n):
            profile_rows.append(dict(condition_id=cid, model=case["model"], angle_label=case["angle_label"], cue_duration=case["cue_duration"],
                material_position=q[j], normal_angle_from_preparation_deg=aq[j], final_active=final_u[j], rotated_basal_active=primary_shift[j]))
        trajectories[str(cid)] = dict(time=t, P_real=P.real, P_imag=P.imag, R=R, angle_deg=angle,
                                     q1=np.abs(M1)/M0, q2=np.abs(M2)/M0, fixed_region_transported_background_excess=E,
                                     fixed_region_physical_active_amount=amount, C=C)
        audits.append(dict(condition_id=cid, exact_basal_restart=True, exact_phase_restart=True, phase_hashes_match=True,
                           pool_relative_error=pool_error, primary_shape_error=primary_shift_rel, full_field_tail_excursion=full_tail))
    frame = pd.DataFrame(rows)
    frame.to_csv(dest("standard_angles_summary.csv"), index=False)
    pd.DataFrame(events).to_csv(dest("standard_angles_fixed_region_crossings.csv"), index=False)
    profiles = pd.DataFrame(profile_rows)
    profiles.to_csv(dest("standard_angles_final_profiles.csv"), index=False)
    contexts = pd.DataFrame(context_rows)
    contexts.to_csv(dest("standard_angles_initial_profiles_and_cues.csv"), index=False)
    pd.read_csv(RUN/"results/initial_context.csv").to_csv(dest("standard_angles_initial_context.csv"), index=False)
    np.savez_compressed(dest("standard_angles_traces.npz"), **{f"c{cid}_{key}": value for cid,z in trajectories.items() for key,value in z.items()})
    render(frame, contexts)
    docs(frame, audits)


def panel_title(model):
    return NAMES[model] + (" · C₀ = 2, basal Cb = 2.044" if model=="spring" else " · C = 40")


def curve_sheet(frame, observable, stem, title, ylabel, limits, targets=False):
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.4))
    fig.subplots_adjust(left=.085, right=.975, bottom=.16, top=.81, hspace=.43, wspace=.28)
    fig.suptitle(title, fontsize=15, y=.985)
    fig.text(.06,.938,"Exact prepared basals; all four requested cue-centre angles, including initial overlaps. Readout after 13,000 cue-free units.",fontsize=10)
    handles = [plt.Line2D([],[], color=COLORS[label], marker="o", ms=4, label=LABELS[label]) for label,_ in native.ANGLES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.52,.91), ncol=4, frameon=False, fontsize=10)
    for i, model in enumerate(native.MODELS):
        ax = axes.flat[i]
        for label, _ in native.ANGLES:
            g = frame[(frame.model==model)&(frame.angle_label==label)].sort_values("cue_duration")
            ax.plot(g.cue_duration,g[observable],"o-",color=COLORS[label],ms=4,lw=1.3)
            failed = g[g.full_field_tail_status=="fails_finite_tail_gate"]
            ax.scatter(failed.cue_duration,failed[observable],facecolor="white",edgecolor=COLORS[label],s=30,linewidth=1.4,zorder=5)
            if targets:
                ax.axhline(g.target_angle_in_P_coordinates_deg.iloc[0],ls=":",lw=.8,color=COLORS[label],alpha=.7)
        weak = frame[(frame.model==model)&(frame.R_final<.1)]
        if len(weak) and targets:
            ax.text(.035,.07,"Direction masked at R < 0.1",transform=ax.transAxes,fontsize=8.5,color="#555555")
        if not targets:
            msg="Non-antipodal curves overlap at R ≈ 1" if model=="holmes_model3" else "All four curves overlap at R ≈ 1"
            ax.text(.04,.1,msg,transform=ax.transAxes,fontsize=8.5,color="#555555")
        support = frame[frame.model==model].initial_angular_support_deg.iloc[0]
        ax.set_title(f"{chr(65+i)}  {panel_title(model)}",loc="left",fontsize=10.8,pad=9)
        ax.text(.02,1.015,f"Full cue width: {support:.2f}°",transform=ax.transAxes,fontsize=8.5,color="#555555")
        ax.set(xscale="log",xlim=(.07 if model=="spring" else .7,1500),ylim=limits,xlabel="Cue duration (native units)",ylabel=ylabel)
        if targets:ax.set_yticks([0,22.5,45,90,135,180])
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(labelsize=9)
    fig.text(.06,.085,"Shared centre angles retain native physical cue widths and peaks; neither cue footprint nor native clocks are matched across models.",fontsize=9)
    fig.text(.06,.058,"Lines connect tested points. Open markers fail the separate full-field tail gate; reused primary-only histories have no invented full-field flag.",fontsize=8.5)
    fig.text(.06,.032,"Initial support–core overlap occurs at: WP 45°/22.5°; Holmes 45°/22.5°; Spring 90°/45°/22.5°. No overlap conditions were excluded.",fontsize=8.5)
    savefig(stem,fig)


def render(frame,contexts):
    plt.rcParams.update({"font.size":9,"axes.labelsize":10,"font.family":"DejaVu Sans","svg.fonttype":"none","pdf.fonttype":42})
    curve_sheet(frame,"angle_final_deg","standard_angles_orientation","Standardized cue angles · Retained activity orientation", "Final activity angle (°)",(-8,190),True)
    upper=max(1.12,float(frame.R_final.max())*1.1)
    curve_sheet(frame,"R_final","standard_angles_magnitude","Standardized cue angles · Relative activity-moment magnitude", "R = |M₁(final)/M₁(initial)|",(-.02,upper))
    fig,axes=plt.subplots(2,2,figsize=(10.8,8.5))
    fig.subplots_adjust(left=.08,right=.98,bottom=.14,top=.83,hspace=.5,wspace=.25)
    fig.suptitle("Initial activity profiles and continuous cue footprints",fontsize=15,y=.98)
    fig.text(.055,.935,"Grey: basal active profile normalized by its own contrast. Colours: cue shape divided by requested peak; dashed line: basal half contrast.",fontsize=9.2)
    handles=[plt.Line2D([],[],color=COLORS[label],label=LABELS[label]) for label,_ in native.ANGLES]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.52,.91),ncol=4,frameon=False)
    for i,model in enumerate(native.MODELS):
        ax=axes.flat[i];z=contexts[(contexts.model==model)&(contexts.angle_label=="pi")].sort_values("normal_angle_from_preparation_deg")
        ax.fill_between(z.normal_angle_from_preparation_deg,0,z.basal_active_normalized_contrast,color="#dddddd")
        ax.plot(z.normal_angle_from_preparation_deg,z.basal_active_normalized_contrast,color="#777777",lw=1.1)
        ax.axhline(.5,color="#aaaaaa",ls="--",lw=.7)
        for label,_ in native.ANGLES:
            g=contexts[(contexts.model==model)&(contexts.angle_label==label)].sort_values("normal_angle_from_preparation_deg")
            ax.plot(g.normal_angle_from_preparation_deg,g.initial_cue_unit_peak,color=COLORS[label],lw=1.2)
        ax.set(xlim=(-120,240),ylim=(-.02,1.07),xlabel="Normal angle from prepared axis (°)",ylabel="Normalized basal profile / cue shape",xticks=[-90,0,45,90,180])
        ax.set_xticklabels(["−90","0","45","90","180"])
        ax.set_title(f"{chr(65+i)}  {panel_title(model)}",loc="left",pad=10,fontsize=10.5)
        ax.spines[["top","right"]].set_visible(False)
    fig.text(.055,.067,"Physical full support: chemical 3 (27° at C = 40); Spring 0.2 (35.22° at its basal circumference). Native peaks: WP 0.1, Goryachev 6.4, Holmes 1, Spring 4.",fontsize=8.7)
    fig.text(.055,.035,"All overlap is retained. Overlapping inputs perturb the original activity region; fixed-region depletion cannot identify a separate challenger or site ancestry.",fontsize=8.7)
    savefig("standard_angles_initial_context",fig)


def docs(frame,audits):
    generated=stamp()
    compact_cols=["model","angle_label","requested_angle_deg","cue_duration","angle_final_deg","R_final","target_error_final_deg","fraction_requested_turn_final","initial_overlap","overlap_fraction_cue_support","full_field_tail_status","final_resolved_peaks","final_primary_rotated_basal_relative_contrast_error"]
    frame[compact_cols].to_csv(dest("standard_angles_endpoints.csv"),index=False)
    caption="Standardized cue-centre angles on four native ring models. A–D, Final physical first-moment orientation after the cue and 13,000 cue-free native time units, for prescribed prepared-axis offsets π, π/2, π/4 and π/8. All 80 requested conditions are retained, including 35 initial overlap cases; 48 trajectories are new native branches and 32 are exact archive reuses with frozen state, cue, dynamics and timing proofs. Chemical basals use fixed C = 40 and durations 1/10/100/300/1000; Spring uses resting C₀ = 2, actual basal Cb = 2.04418 and durations 0.1/1/10/100/1000. The cue centres use the declared preparation axis without a fresh argmax; plotted P = M₁(t)/M₁(0) uses the measured initial physical first-moment phase and stores the exact target offset. Spring π/8 remains at its continuous half-cell centre. Dotted lines mark prescribed target normals in P coordinates; lines between measurements guide the eye. Direction is masked when R = |P| < 0.1. Open markers flag failure of the independent full-field final-1000-unit excursion gate; reused primary-only histories carry an explicit unavailable full-field flag. All native physical cue widths and peaks are preserved: chemical full support 3 (27°), Spring instantaneous physical full support 0.2 (35.22° initially), with peaks 0.1/6.4/1/4 for WP/Goryachev/Holmes/Spring. Shared centre angles do not match cue footprint, dose or native timescales across models. Initial support–core intersections and cue-weighted initial activity are exported; overlapping cases describe local-region reorientation, not competition between initially separate sites. Activity orientation is not a movement law. Antipodal inputs preserve a reflection axis; their axial outcomes do not establish a distinct binary mechanism."
    dest("CAPTION.md","metadata").write_text(caption+"\n")
    width=frame.groupby("model").initial_angular_support_deg.first()
    overlap=frame.drop_duplicates(["model","angle_label"])[["model","angle_label","overlap_fraction_cue_support","initial_unit_peak_cue_excess_fraction_total"]]
    overlap_table = " | ".join(overlap.columns)+"\n"+" | ".join(["---"]*len(overlap.columns))+"\n"+"\n".join(" | ".join(str(x) for x in r) for r in overlap.itertuples(index=False,name=None))
    lines=["# Standardized-angle checkpoint", "", f"Generated: {generated}.", "", "The user requested the same sensible centre angles before judging which geometry claim should lead Figure 3. This separate checkpoint preserves all current main and earlier prototype figures. It does not choose a final manuscript layout.", "", "The frozen 80-condition experiment uses the exact current basals and native protocols, with 32 proved exact reuses and 48 new SLURM branches (array 21187059). All 80 complete; all 35 overlap rows remain. The primary readouts are retained angle and magnitude; no arbitrary alignment threshold or refined minimum duration is promoted.", "", "## Inputs and overlap", "", "Cue centres are the declared prepared-axis coordinate plus the requested fraction, never a plateau argmax. Chemical preparation axis has absolute normal 0.225°; Spring 0.9°. Measured initial physical M₁ phases and their small offsets are recorded per condition. Spring π/8 uses continuous material centre 0.065, without snapping to a cell or renormalizing its slightly reduced sampled maximum. The analytic peak remains 4. Both the sampled shape hash and actual physical cue integral are frozen.", "", "The support–core overlap fraction divides the physical intersection between the continuous full cue-support interval and linearly interpolated basal half-contrast core by the full physical cue support. This geometric measure is distinct from cue-weighted basal excess activity divided by total basal excess. Thus Goryachev can have no half-contrast-core intersection while still receiving input on nonzero basal tails.", "", overlap_table, "", "Native angular full widths are " + "; ".join(f"{NAMES[m]} {w:.5f}°" for m,w in width.items()) + ". Shared centre angles are not matched footprint or dose controls.", "", "## Readouts and limits", "", "The physical moments use instantaneous circumference, Mₙ = C(t) mean[u exp(2πinq)]. P = M₁/M₁(0), R = |P|, qₙ = |Mₙ|/M₀. The angle is reported only for R ≥ 0.1, with raw angles also retained. Target error is circular; the fractional requested turn divides the retained angle by the requested prepared-axis offset. This exposes why a common ±15° success tolerance is especially permissive at 22.5°. No shortest-duration claim is inferred here.", "", "Fixed-region summaries use the original basal material half-contrast interval and physical amount, subtracting its transported initial background amount before normalization. They describe regional activity, potentially changed by translation, redistribution or chemistry; they do not distinguish a unique incumbent lineage, destruction or separate challenger when the initial cue overlaps the same region. Crossing brackets are saved-sample brackets, with censoring and rebound flags.", "", "Full native endpoint profiles and their Fourier-shifted basal counterparts remain exported. A full-field final-tail gate is separately reported as the largest species-scale-normalized excursion over the last 1000 cue-free units (≤10⁻⁵ passes); compact reused primary-only trajectories have no invented full-field qualification. Native pool conservation and exact state/phase restarts are replayed independently of worker summaries. A finite 13,000-unit withdrawal does not prove asymptotic attraction. Native clocks are not calibrated across models; the global mean-force Spring closure does not supply a local mechanical propagation time.", "", "## Files", "", "- `standard_angles_summary.csv`: all 80 conditions, exact archive paths/hashes, anchor/overlap diagnostics, continuous endpoints and field flags.", "- `standard_angles_endpoints.csv`: compact review table.", "- `standard_angles_initial_context.csv`: 16 frozen input/overlap records.", "- `standard_angles_initial_profiles_and_cues.csv`: actual sampled basal/cue fields.", "- `standard_angles_final_profiles.csv`: all final primary fields and rotated basal comparisons.", "- `standard_angles_fixed_region_crossings.csv`: threshold brackets with censoring and rebound flags.", "- `standard_angles_traces.npz`: native sampled moment, geometry and fixed-region summaries.", "", "The separate supporting magnitude and initial-input-context sheets retain PNG/PDF/SVG exports. Existing geometry and selected near/far duration sheets remain unchanged."]
    dest("EVIDENCE.md","metadata").write_text("\n".join(lines)+"\n")
    with dest("EVIDENCE.md","metadata").open("a") as f:
        nonantipodal=frame[frame.angle_label!="pi"]
        f.write(f"\nAll {len(nonantipodal)} non-antipodal final primary profiles contain one resolved peak and match a fractional rotation of their own exact basal primary profile to at most {nonantipodal.final_primary_rotated_basal_relative_contrast_error.max():.4g} of its initial contrast. Their R values range from {nonantipodal.R_final.min():.15g} to {nonantipodal.R_final.max():.15g}. These native profile checks distinguish a translated full-strength region from vector cancellation.\n")
        f.write(r"""

The physical moments and fixed-region excess are

$$
M_n(t)=C(t)\int_0^1 u(q,t)e^{2\pi i n q}\,\mathrm{d}q,
\qquad P(t)=\frac{M_1(t)}{M_1(0)},
\qquad R(t)=|P(t)|.
\tag{1}
$$

$$
E_{\Omega_0}^{\mathrm{tr}}(t)=
\frac{C(t)\int_{\Omega_0}u(q,t)\,\mathrm{d}q
      -C_{\mathrm b}u_{\mathrm{bg},0}|\Omega_0|}
     {C_{\mathrm b}\int_{\Omega_0}
       [u(q,0)-u_{\mathrm{bg},0}]\,\mathrm{d}q}.
\tag{2}
$$

Here \(\Omega_0\) is the fixed basal half-contrast material interval; the subtracted background is its transported initial physical amount. Pure affine dilution of an otherwise unchanged material profile preserves this excess.
""")
    allocated=sum(p.stat().st_blocks*512 for p in RUN.rglob("*") if p.is_file())
    validation=dict(generated_at=generated,passed=True,requested=80,complete=80,new_native=48,exact_reuse=32,initial_overlap_conditions=int(frame.initial_overlap.sum()),excluded_for_overlap=0,
        manifest_sha256=sha(RUN/"manifest.json"),all_frozen_sources_and_inputs_hash_verified=True,all_actual_native_archive_hashes_verified=True,
        all_exact_initial_and_phase_restarts=True,maximum_pool_relative_error=float(frame.conserved_pool_max_relative_error.max()),
        field_tail_counts=frame.full_field_tail_status.value_counts().to_dict(),records=audits,
        raw_campaign_allocated_bytes=allocated,raw_campaign_large_manifest_required=bool(allocated>=2*1024**3))
    write_json("validation.json",validation)
    content=dict(title="Standardized prescribed-angle duration checkpoint",role="Separate user-requested review checkpoint; existing main/prototypes preserved.",
        assets=[dict(path=str((OUT/"standard_angles_orientation.png").relative_to(REPO)),caption="Four-model comparison at π, π/2, π/4 and π/8, retaining all initial overlaps.")],
        legend=[caption],sources=[dict(label="Complete 80-condition continuous endpoints and provenance",path=str(dest("standard_angles_summary.csv").relative_to(REPO))),dict(label="Input geometry and overlap diagnostics",path=str(dest("standard_angles_initial_context.csv").relative_to(REPO))),dict(label="Evidence, definitions and limits",path=str(dest("EVIDENCE.md","metadata").relative_to(REPO)))],
        supplementary=[dict(title="Relative activity-moment magnitude",path=str((OUT/"standard_angles_magnitude.png").relative_to(REPO))),dict(title="Initial native profile and cue footprints",path=str((OUT/"standard_angles_initial_context.png").relative_to(REPO)))])
    write_json("content_checkpoint.json",content)
    relative_run=os.path.relpath(RUN,OUT)
    (OUT/"README.md").write_text(f"""# Standardized angle–duration checkpoint

This is the user-requested stationary-cue experiment at π, π/2, π/4 and π/8. It preserves the previous Figure 3 sheets while the experiment and interpretation are reviewed.

- [Orientation response](standard_angles_orientation.png) · [PDF](supporting/vectors/standard_angles_orientation.pdf)
- [Activity-moment magnitude](standard_angles_magnitude.png)
- [Initial basal profiles and cue footprints](standard_angles_initial_context.png)
- [Compact endpoint table](supporting/data/standard_angles_endpoints.csv)
- [Complete endpoints, native provenance and field flags](supporting/data/standard_angles_summary.csv)
- [Input anchor and overlap diagnostics](supporting/data/standard_angles_initial_context.csv)
- [Caption](supporting/metadata/CAPTION.md) · [Evidence and definitions](supporting/metadata/EVIDENCE.md) · [Validation](supporting/metadata/validation.json)

All 80 conditions are complete: 48 new native branches and 32 proved exact archive reuses. All 35 initial-overlap conditions remain. Chemical tests use fixed C = 40 and native durations 1/10/100/300/1000; Spring uses resting C₀ = 2 and native durations 0.1/1/10/100/1000. Every endpoint follows 13,000 cue-free units. Cue-centre angles are shared; native physical widths, peaks and clocks remain model-specific.

All 60 non-antipodal final primary profiles have one resolved peak and match a rotated copy of their exact basal profile to within 2.21×10⁻¹² of its contrast. Their graded retained angles describe one translated full-strength activity region. Antipodal Holmes instead loses first-moment magnitude: at T = 1000, R ≈ 0.040 and direction is masked. Native profiles and magnitude distinguish these outcomes.

Raw angle and R are primary. Target error and fractional requested turn remain exported; no refined minimum-duration or arbitrary small-angle alignment claim is made. Overlapping inputs perturb the original activity region, so fixed-region loss does not identify a separate challenger or site lineage. Activity orientation is not measured cell movement.

The [frozen native manifest]({relative_run}/manifest.json) and [collection counts]({relative_run}/results/collection.json) record the full experiment. SLURM array: 21187059. New native allocated storage: {allocated/1024**2:.1f} MiB.

To rebuild the three plotted sheets from retained compact tables, run `revision_standard_angles_report_20261019.py --render-only` in the existing environment. The default command independently replays the native archives before generating tables, plots and validation; that replay requires the preserved native data checkout.
""")
    record(__file__)
    write_json("source_sha256.json",SOURCE)
    print(json.dumps(dict(passed=True,complete=80,overlap_conditions=int(frame.initial_overlap.sum()),maximum_pool_error=float(frame.conserved_pool_max_relative_error.max()),field_tail_counts=validation["field_tail_counts"],raw_allocated_MiB=allocated/1024**2,output=str(OUT)),indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--render-only",action="store_true",help="Rebuild art from existing compact tables without changing native validation or tables")
    args=parser.parse_args()
    if args.render_only:
        render(pd.read_csv(dest("standard_angles_summary.csv")),pd.read_csv(dest("standard_angles_initial_profiles_and_cues.csv")))
    else:
        main()
