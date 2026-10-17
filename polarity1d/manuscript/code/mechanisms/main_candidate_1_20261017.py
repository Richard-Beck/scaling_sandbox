"""Six-panel chemical-response candidate, using retained evidence only.

The previous draft's audit checks provenance and exact states. This renderer
additionally reconstructs reaction and diffusion in the original fixed core;
it never calls an integration routine or writes to retained simulation data.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
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
from matplotlib.ticker import MaxNLocator

import chemistry_main_text_candidate_20261017 as retained
from length_dev_output_layout import dev_output_path, dev_artifact_key

ROOT = retained.ROOT
OUT = ROOT / "figures/length_mechanisms_20261008/dev/main_candidate_1_20261017"
STEM = "main_candidate_1"
COLORS = retained.LENGTH_COLORS
LABELS = {"goryachev": "Goryachev", "wave_pinning": "Wave-pinning"}
REVIEW_ROOT = ROOT.parents[1] / "current manuscript figs/reviews"


def write_json(name, obj):
    path = dev_output_path(OUT, name)
    path.write_text(json.dumps(obj, indent=2) + "\n")
    return path


def write_csv(name, records):
    path = dev_output_path(OUT, name)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path


def local_rates(name, u, v, length):
    """Original cell-centred, no-flux FV operator; no integration."""
    dx = length / u.shape[-1]
    lap = np.empty_like(u)
    lap[..., 0] = (u[..., 1] - u[..., 0]) / dx**2
    lap[..., -1] = (u[..., -2] - u[..., -1]) / dx**2
    lap[..., 1:-1] = (u[..., :-2] - 2 * u[..., 1:-1] + u[..., 2:]) / dx**2
    reaction = retained.direct_reaction(name, u, v)
    diffusion = (0.0025 if name == "goryachev" else 0.1) * lap
    return reaction, diffusion


def activation_factor(name, u):
    if name == "goryachev":
        return (0.33 * u * u + 0.67 * u) / 101
    return 0.067 + u * u / (1 + u * u)


def collect():
    records, profiles, audit = retained.audit_and_measure()
    summary = retained.rows(retained.POOL / "response_vs_length.csv")
    review_csv = REVIEW_ROOT / "incumbent_r2_return_state_decomposition.csv"
    with review_csv.open(newline="") as stream:
        independent = list(csv.DictReader(stream))
    retained.SOURCES[os.path.relpath(review_csv, ROOT)] = retained.sha(review_csv)
    review_md = REVIEW_ROOT / "incumbent_r2.md"
    retained.SOURCES[os.path.relpath(review_md, ROOT)] = retained.sha(review_md)
    r3_csv = REVIEW_ROOT / "incumbent_r3_capture_factorization.csv"
    with r3_csv.open(newline="") as stream:
        r3_factors = list(csv.DictReader(stream))
    for path in (r3_csv, REVIEW_ROOT / "incumbent_r3.md"):
        retained.SOURCES[os.path.relpath(path, ROOT)] = retained.sha(path)
    curves, checks, maxima = [], [], dict(rate_error=0.0, deficit_error=0.0,
                                         allocation_pointwise_error=0.0,
                                         returned_global_mass_error=0.0)
    for row in records:
        name, length = row["model"], row["length_native"]
        original = next(x for x in summary if x["model"] == name
                        and float(x["length"]) == length
                        and float(x["duration"]) == 0.1
                        and float(x["epsilon"]) == 0.005)
        path = ROOT / row["perturbation_archive"]
        z = retained.load(path)
        sham = retained.load(path.with_name("epsilon_0.npz"))
        width = row["fixed_core_width_native"]
        dx = length / len(z["x"])
        weights = np.clip(width - np.arange(len(z["x"])) * dx, 0.0, dx)
        r, d = local_rates(name, z["returnedstate"][0], z["returnedstate"][1], length)
        rs, ds = local_rates(name, sham["returnedstate"][0], sham["returnedstate"][1], length)
        damage = float((sham["returnedstate"][0] - z["returnedstate"][0]) @ weights)
        chemical = float((r - rs) @ weights)
        diffusive = float((d - ds) @ weights)
        total = chemical + diffusive
        up, vp = z["returnedstate"]
        us, vs = sham["returnedstate"]
        c = 0.01733 if name == "goryachev" else 1.0
        active_field = (activation_factor(name, up) - activation_factor(name, us)) * vs - c * (up - us)
        inactive_surplus = activation_factor(name, up) * (vp - vs)
        residual = float(np.max(np.abs(active_field + inactive_surplus - (r - rs))))
        np.testing.assert_allclose(active_field + inactive_surplus, r - rs, rtol=3e-8, atol=4e-14)
        active_gain = float(active_field @ weights) / damage
        substrate_gain = float(inactive_surplus @ weights) / damage
        np.testing.assert_allclose(active_gain + substrate_gain + diffusive / damage,
                                   total / damage, rtol=3e-8, atol=3e-9)
        global_v_surplus = float((vp - vs).sum() * dx)
        global_u_loss = float((us - up).sum() * dx)
        mass_error = abs(global_v_surplus - global_u_loss)
        np.testing.assert_allclose(global_v_surplus, global_u_loss, rtol=2e-7, atol=3e-11)
        maxima["allocation_pointwise_error"] = max(maxima["allocation_pointwise_error"], residual)
        maxima["returned_global_mass_error"] = max(maxima["returned_global_mass_error"], mass_error)
        review = next(r for r in independent if r["model"] == name and float(r["length_native"]) == length)
        for key, value in (("active_state_term_per_damage", active_gain),
                           ("inactive_state_term_per_damage", substrate_gain),
                           ("diffusive_gain_per_damage", diffusive / damage),
                           ("net_gain_per_damage", total / damage)):
            np.testing.assert_allclose(value, float(review[key]), rtol=3e-7, atol=3e-8)
        for key, value in (("excess_net_activation_at_release", chemical),
                           ("excess_total_material_rate_at_release", total),
                           ("reaction_gain_per_damage", chemical / damage),
                           ("total_gain_per_damage", total / damage)):
            np.testing.assert_allclose(value, float(original[key]), rtol=2e-9, atol=3e-11)
        row.update(excess_core_diffusion_native=diffusive,
                   excess_core_net_rate_native=total,
                   diffusion_gain_per_damage_native=diffusive / damage,
                   net_gain_per_damage_native=total / damage,
                   chemical_fraction_of_immediate_net_gain=chemical / total,
                   sustained_90_percent_recovery_time_native=float(original["ninety_percent_damage_recovery_time"]),
                   active_field_at_sham_inactive_per_damage_native=active_gain,
                   inactive_surplus_at_perturbed_active_per_damage_native=substrate_gain,
                   active_field_plus_native_diffusion_per_damage_native=active_gain + diffusive / damage,
                   global_inactive_surplus_at_return_native=global_v_surplus,
                   global_active_loss_at_return_native=global_u_loss,
                   global_inactive_surplus_per_core_damage=global_v_surplus / damage,
                   returned_global_mass_error_native=mass_error,
                   allocation_max_pointwise_error_native=residual)
        if int(length) in COLORS:
            for phase in ("pulse", "recovery"):
                u, v = z[f"{phase}_active"], z[f"{phase}_inactive"]
                us, vs = sham[f"{phase}_active"], sham[f"{phase}_inactive"]
                np.testing.assert_array_equal(z[f"{phase}_time"], sham[f"{phase}_time"])
                rp, dp = local_rates(name, u, v, length)
                r0, d0 = local_rates(name, us, vs, length)
                excess_r = (rp - r0) @ weights
                excess_d = (dp - d0) @ weights
                deficit = (us - u) @ weights
                sham_core = us @ weights
                recorded_deficit = sham[f"{phase}_core_active"] - z[f"{phase}_core_active"]
                recorded_total = z[f"{phase}_core_total_rate"] - sham[f"{phase}_core_total_rate"]
                np.testing.assert_allclose(deficit, recorded_deficit, rtol=2e-7, atol=4e-12)
                np.testing.assert_allclose(excess_r + excess_d, recorded_total, rtol=2e-7, atol=4e-12)
                maxima["deficit_error"] = max(maxima["deficit_error"], float(np.max(np.abs(deficit - recorded_deficit))))
                maxima["rate_error"] = max(maxima["rate_error"], float(np.max(np.abs(excess_r + excess_d - recorded_total))))
                times = z[f"{phase}_time"] + (-0.1 if phase == "pulse" else 0)
                for i, time in enumerate(times):
                    curves.append(dict(model=name, length_native=length, phase=phase,
                                       time_from_return_native=float(time),
                                       active_core_deficit_native=float(deficit[i]),
                                       active_core_loss_percent=float(100 * deficit[i] / sham_core[i]),
                                       active_deficit_over_damage_at_return=float(deficit[i] / damage),
                                       excess_core_reaction_native=float(excess_r[i]),
                                       excess_core_diffusion_native=float(excess_d[i]),
                                       excess_core_net_rate_native=float(excess_r[i] + excess_d[i])))
        checks.append(dict(model=name, length_native=length,
                           immediate_reaction_and_diffusion_recomputed=True,
                           original_summary_agrees=True,
                           ordered_active_substrate_split_reconstructed=True,
                           returned_global_inactive_surplus_matches_active_loss=True,
                           independent_R2_decomposition_agrees=True))
    for name in LABELS:
        rr = [r for r in records if r["model"] == name]
        baseline = rr[0]["net_gain_per_damage_native"]
        for r in rr:
            r["net_gain_relative_to_L9"] = r["net_gain_per_damage_native"] / baseline
        if name == "goryachev":
            assert all(r["active_field_at_sham_inactive_per_damage_native"] < 0 for r in rr)
            assert all(r["active_field_plus_native_diffusion_per_damage_native"] < 0 for r in rr)
    audit["new_candidate_reanalysis"] = dict(
        new_integrations=0, exact_return_state_budgets=36,
        full_selected_matched_histories=6, full_selected_phases=12,
        independent_core_budget_max_absolute_errors=maxima,
        ordered_split="Active change first at sham inactive field; inactive change second at perturbed active field. Cross term belongs to inactive contribution.",
        decomposition_is="Exact fixed-state algebraic allocation; no separately evolved pathway intervention.",
        checks=checks)
    r3_product_error = 0.0
    for factor in r3_factors:
        record = next(r for r in records if r["model"] == factor["model"]
                      and r["length_native"] == float(factor["length_native"]))
        product = float(factor["Abar_delta_v_native"]) * float(factor["surplus_exposure_per_damage"])
        r3_product_error = max(r3_product_error, abs(product - float(factor["substrate_allocation_per_damage_native"])))
        # Imported diagnostics use independently subtracted very small deficits;
        # use the existing budget's cancellation tolerance rather than bit identity.
        np.testing.assert_allclose(product, float(factor["substrate_allocation_per_damage_native"]), rtol=2e-9, atol=3e-10)
        np.testing.assert_allclose(product, record["inactive_surplus_at_perturbed_active_per_damage_native"], rtol=3e-7, atol=3e-8)
    audit["R3_imported_factorization"] = dict(rows=len(r3_factors), source=str(r3_csv),
                                               imported_products_match_existing_substrate_allocation=True,
                                               max_absolute_imported_product_error=r3_product_error,
                                               new_numerical_reanalysis=False)
    return records, profiles, curves, audit, r3_factors


def render(records, profiles, curves):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.titlesize": 10, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": .7,
                         "pdf.fonttype": 42, "svg.fonttype": "none",
                         "lines.linewidth": 1.65})
    fig = plt.figure(figsize=(7.6, 8.9))
    positions = [(0.108, .735, .302, .170), (.602, .735, .302, .170),
                 (.108, .365, .302, .175), (.602, .365, .302, .175),
                 (.108, .102, .302, .155), (.602, .102, .302, .155)]
    axes = [fig.add_axes(pos) for pos in positions]
    profile_handles = [Line2D([], [], color=c, label=f"L = {L}") for L, c in COLORS.items()]
    profile_handles += [Line2D([], [], color="#333333", label="Active u"),
                        Line2D([], [], color="#333333", ls=(0, (4, 2)), label="Inactive v")]
    fig.legend(handles=profile_handles,
               loc="upper center", bbox_to_anchor=(.50, .998), ncol=5,
               frameon=False, columnspacing=1.3, handlelength=2.0)
    for col, (name, label) in enumerate(LABELS.items()):
        ax = axes[col]
        twin = ax.twinx()
        twin.spines["right"].set_visible(True)
        for length, color in COLORS.items():
            rows = [r for r in profiles if r["model"] == name and r["length"] == length]
            x = np.array([r["position_native"] for r in rows])
            ax.plot(x, [r["active_native"] for r in rows], color=color)
            twin.plot(x, [r["inactive_native"] for r in rows], color=color, ls=(0, (4, 2)), lw=1.4)
        ax.set(xlabel="Position, x (native units)", ylabel="Active u (native conc.)",
               xlim=(0, 3 if col == 0 else 40),
               ylim=(0, 125 if col == 0 else 1.25))
        twin.set(ylabel="Inactive v (native conc.)", ylim=(0, .30 if col == 0 else 2))
        ax.set_xticks([0, 1, 2, 3] if col == 0 else [0, 20, 40])
        ax.yaxis.set_major_locator(MaxNLocator(4))
        twin.yaxis.set_major_locator(MaxNLocator(4))
        fig.text(positions[col][0] + positions[col][2] / 2, .941, label,
                 ha="center", fontsize=12, fontweight="bold")
        ax.set_title("Basal incumbent" + (": end view" if col == 0 else ": full cell"), pad=8)
        rr = [r for r in records if r["model"] == name]
        first, last = rr[0], rr[-1]
        vtext = (f"Opposite-end v: {first['opposite_inactive_native']:.3f} → "
                 f"{last['opposite_inactive_native']:.3f} ({last['initial_activation_relative_to_L9']:.3f}×)" if col == 0 else
                 "Opposite-end v: 1.808 → 1.808 (1.000×)")
        law = (r"$A(u)=\frac{0.33u^2+0.67u}{101},\quad c=0.01733$" if col == 0 else
               r"$A(u)=0.067+\frac{u^2}{1+u^2},\quad c=1$")
        fig.text(positions[col][0] + positions[col][2] / 2, .644, law,
                 fontsize=9, ha="center")
        fig.text(positions[col][0] + positions[col][2] / 2, .615, vtext,
                 fontsize=8.2, ha="center")
        ax = axes[2 + col]
        for length, color in COLORS.items():
            pulse = [r for r in curves if r["model"] == name and r["length_native"] == length
                     and r["phase"] == "pulse"]
            recovery = [r for r in curves if r["model"] == name and r["length_native"] == length
                        and r["phase"] == "recovery" and r["time_from_return_native"] <= 1]
            series = pulse + recovery[1:]
            ax.plot([r["time_from_return_native"] for r in series],
                    [r["active_core_loss_percent"] for r in series], color=color)
        ax.axvspan(-.1, 0, color="#e8e8e8", zorder=0)
        ax.axvline(0, color="#777777", ls=":", lw=.9)
        ax.axhline(0, color="#aaaaaa", lw=.6, zorder=0)
        ax.set(xlim=(-.12, 1), ylim=(-.00004, .00110) if col == 0 else (-.001, .053),
               xlabel="Time from return (native units)", ylabel="Core active loss (%)",
               xticks=[-.1, 0, .5, 1])
        ax.set_xticklabels(["−0.1", "0", "0.5", "1"])
        ax.get_xticklabels()[0].set_horizontalalignment("right")
        ax.get_xticklabels()[1].set_horizontalalignment("left")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_title("Damage and early recovery", pad=8)
        ax = axes[4 + col]
        length_grid = np.array([r["length_native"] for r in rr])
        active = np.array([r["active_field_at_sham_inactive_per_damage_native"] for r in rr])
        substrate = np.array([r["inactive_surplus_at_perturbed_active_per_damage_native"] for r in rr])
        diffusion = np.array([r["diffusion_gain_per_damage_native"] for r in rr])
        total = np.array([r["net_gain_per_damage_native"] for r in rr])
        ax.plot(length_grid, substrate, color="#999999", ls=(0, (5, 2)), lw=1.8, zorder=2)
        ax.plot(length_grid, active, color="#555555", ls=":", lw=1.7, zorder=2)
        ax.plot(length_grid, diffusion, color="#aaaaaa", ls="-.", lw=1.2, zorder=2)
        ax.plot(length_grid, total, color="#222222", lw=0, marker="o", ms=3.8,
                mfc="white", mec="#222222", mew=.65, zorder=3)
        for L, color in COLORS.items():
            i = int(np.flatnonzero(length_grid == L)[0])
            ax.plot(L, total[i], marker="o", color=color, ms=3.5, zorder=4)
        ax.axhline(0, color="#bbbbbb", lw=.65, zorder=0)
        ax.set(xlabel="Cell length, L (native units)", ylabel="Allocated rate / damage\n(1 / native time)",
               xlim=(7, 41), xticks=[9, 21, 39], ylim=(-.6, 17) if col == 0 else (-.025, .88))
        ax.set_title("Returned-state rate allocation", pad=8)
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.text(.03, .93 if col == 0 else .48, f"Net L39 / L9: {last['net_gain_relative_to_L9']:.2f}×", transform=ax.transAxes,
                ha="left", va="top", fontsize=8.0)
        if col == 0:
            ax.text(.03, .71, "Surplus supplies rebound\nActive-field term < 0\n−0.0252 to −0.0235", transform=ax.transAxes,
                    ha="left", va="top", fontsize=7.7, color="#444444")
    fig.text(.50, .676, r"$f(u,v)=A(u)v-cu$; equal local opposing cue: $\Delta f_{\mathrm{opp}}=s\,v_{\mathrm{opp}}$",
             ha="center", fontsize=9)
    fig.text(.50, .586, "Remove 0.5% of whole inactive profile at t = −0.1; return it pointwise at t = 0",
             ha="center", fontsize=8.3)
    fig.text(.50, .301, "After exact return: global inactive surplus balances global active loss",
             ha="center", fontsize=8.1)
    fig.legend(handles=[Line2D([], [], color="#999999", ls=(0, (5, 2)), lw=1.8,
                               label="Inactive surplus at perturbed u"),
                        Line2D([], [], color="#555555", ls=":", lw=1.7,
                               label="Active-field change at sham v"),
                        Line2D([], [], color="#aaaaaa", ls="-.", lw=1.2, label="Diffusion"),
                        Line2D([], [], color="#222222", lw=0, marker="o", mfc="white", ms=3.8,
                               label="Net original-core rate")],
               loc="lower center", bbox_to_anchor=(.50, .002), frameon=False,
               ncol=2, handlelength=2.6, columnspacing=2, fontsize=7.5)
    for letter, ax in zip("ABCDEF", axes):
        pos = ax.get_position()
        fig.text(pos.x0 - .069, pos.y1 + .019, letter, fontsize=12, fontweight="bold")
        ax.grid(axis="y", lw=.5, alpha=.14)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    clipped = []
    for text in fig.findobj(matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if box.width and box.height and not (box.x0 >= bounds.x0 and box.x1 <= bounds.x1
                                                 and box.y0 >= bounds.y0 and box.y1 <= bounds.y1):
                clipped.append(text.get_text())
    assert not clipped, clipped
    outputs = []
    for extension in ("png", "pdf", "svg"):
        path = dev_output_path(OUT, f"{STEM}.{extension}")
        fig.savefig(path, dpi=320, facecolor="white")
        outputs.append(path)
    plt.close(fig)
    return outputs, dict(panel_count=6, text_within_page_bounds=True,
                         model_columns=True, length_colors_consistent=True,
                         size_inches=[7.6, 8.9], png_dpi=320,
                         manual_visual_review="Passed: revised PNG inspected. Native activation laws, protocol, global mass-balance statement and allocation legend are separated from axes. Goryachev terms retain a linear axis with direct sign/magnitude annotation; net markers leave the dominant substrate curve visible. Length colors remain consistent.")


CAPTION = """Candidate 1, revision 2. Incumbent substrate capture and remote substrate availability in two conserved chemical models. Columns retain Goryachev and Wave-pinning; colors identify L = 9, 21 and 39 throughout. A–B: exact basal active u (solid, left axes) and inactive v (dashed, right axes), with the native basal law f(u,v) = A(u)v − cu. Goryachev is explicitly magnified to the incumbent end, x = 0–3; Wave-pinning shows the full physical domains. The opposite-end inactive values below compare L9 to L39 using the last native cell of each full state. Initial added reaction from an equal local opposing cue is exactly s times that local inactive concentration; its L39/L9 ratio is 0.333 in Goryachev and 0.99990 in Wave-pinning. Goryachev's activation factor grows strongly with its much taller incumbent; Wave-pinning mainly broadens while its sampled active plateau and local conversion factor change little. Fixed mean conserved-pool concentration means total material grows with length by construction; local parameters are unchanged.

C–D: active material lost from the fixed original boundary-connected half-contrast core, relative to its matched sham, as a percentage of that sham's core active material. The original inactive profile is reduced by 0.5% pointwise at t = −0.1 and held in an external nonreacting reservoir; both chemical fields continue evolving normally. At t = 0, exactly the removed original profile is added to the evolved inactive field, without resetting either field. Gray shading marks the hold. The curves show the original logarithmically sampled histories on linear axes through the final retained sample at or before t = 1; full 5000-unit recovery histories are exported. Active concentration is continuous at removal and return. Within-model spatial core and edge-cell weights are held fixed throughout each perturbation and sham. Fractional damage at return falls from 0.001002% to 0.000488% in Goryachev and from 0.04915% to 0.04821% in Wave-pinning; the damage itself is not matched across lengths or models.

E–F: exact returned-state allocation of the immediate original-core excess rate, divided by that core's active loss at return, at all 18 lengths. First change active from sham to perturbed while holding inactive at sham: [A(up) − A(us)]vs − c(up − us). Then change inactive from sham to perturbed at the perturbed active field: A(up)(vp − vs). The cross term belongs to the second contribution. Each expression is integrated over the fixed original core; native active-diffusion difference completes the actual net rate. This is a fixed-state algebraic allocation, not independently evolved pathway ablations. All rate axes remain linear. In Goryachev the active-field contribution is negative at every length, −0.0252 to −0.0235 per native time; the inactive-surplus term rises from 0.997 to 15.309 and supplies the positive chemical rebound under this ordering. In Wave-pinning both chemical terms are positive: active-field 0.097 to 0.187 and inactive-surplus 0.520 to 0.603. Net core repair per damage rises 15.64-fold in Goryachev and 1.25-fold in Wave-pinning. Colored markers identify the same selected lengths as A–D.

From L9 to L39, Goryachev's substrate contribution grows 15.35-fold through a 16.81-fold increase in its inactive-surplus-weighted native conversion coefficient and a 0.913-fold change in core surplus per core damage. The corresponding Wave-pinning factorization is 1.159-fold = 1.093-fold coefficient × 1.060-fold exposure. These are properties of the actual coupled returned states.

The returned-state inactive surplus is endogenous: globally it balances active loss relative to sham after the exact originally removed material has been returned. There is no external oversupply. The native laws and observed allocation support a substrate-capture interpretation under this disturbance, while reduced opposite-end inactive concentration limits initial cue-driven activation there. They do not independently measure a general capture capacity or intrinsic active-field stability, and no simultaneously competing front is created. Absolute native rates are not calibrated between models. Sources: current length_mechanisms.pdf Figures 3–4 and retained native archives, with the ordered split checked against fresh review R2. No new integrations; this is a revised provisional main-text candidate."""


MATH_METHODS = r"""# Caption and precise local measurements

""" + CAPTION + r"""

## Definitions

Let the superscripts \(p\) and \(s\) denote perturbed and matched-sham fields. Time zero is the pointwise inactive-profile return. The spatial region \(\Omega_0\) is each length's original boundary-connected half-contrast incumbent core, held fixed during both histories. All integrals below use the native finite-volume cell widths and a fractional weight for the core's last cell. The unchanged native laws are

$$
\begin{aligned}
f(u,v)&=A(u)v-cu,\\[3pt]
A_{\mathrm G}(u)&=\frac{0.33u^2+0.67u}{101},\qquad c_{\mathrm G}=0.01733,\\[3pt]
A_{\mathrm{WP}}(u)&=0.067+\frac{u^2}{1+u^2},\qquad c_{\mathrm{WP}}=1.
\end{aligned}
\tag{1}
$$

The initial added chemical source from an equal local opposing cue is

$$
\Delta f_{\mathrm{opp}}
=f(u_{\mathrm{opp}},v_{\mathrm{opp}};s)-f(u_{\mathrm{opp}},v_{\mathrm{opp}};0)
=s\,v_{\mathrm{opp}}.
\tag{2}
$$

This source identity does not include the basal reaction or diffusion. The original-core active deficit, and the plotted fractional loss, are

$$
\begin{aligned}
\mathcal D(t)&=\int_{\Omega_0}\!\left[u^{s}(x,t)-u^{p}(x,t)\right]\,\mathrm dx,\\[3pt]
\ell(t)&=100\,\frac{\mathcal D(t)}{\displaystyle\int_{\Omega_0}\!u^{s}(x,t)\,\mathrm dx}.
\end{aligned}
\tag{3}
$$

At the actual returned state, define \(\Delta u=u^p-u^s\) and \(\Delta v=v^p-v^s\). The exact ordered chemical split is

$$
\begin{aligned}
\Delta f
={}&\underbrace{\left[A(u^p)-A(u^s)\right]v^s-c\,\Delta u}_{\text{active-field change at sham inactive field}}\\
&+\underbrace{A(u^p)\,\Delta v}_{\text{inactive-surplus capture at perturbed active field}}.
\end{aligned}
\tag{4}
$$

The active change is evaluated first, and the inactive change second. This ordering assigns their cross term to the inactive contribution; reversing the order reallocates it. The allocated immediate rates, all evaluated at \(t=0^+\), are

$$
\begin{aligned}
\mathcal E_u&=\int_{\Omega_0}\!\left(\left[A(u^p)-A(u^s)\right]v^s-c\,\Delta u\right)\,\mathrm dx,\\[3pt]
\mathcal E_v&=\int_{\Omega_0}\!A(u^p)\,\Delta v\,\mathrm dx,\\[3pt]
\mathcal E_D&=\int_{\Omega_0}\!D_u\left[\frac{\partial^2u^{p}}{\partial x^2}-\frac{\partial^2u^{s}}{\partial x^2}\right]_{t=0^+}\,\mathrm dx.
\end{aligned}
\tag{5}
$$

The original cell-centred no-flux finite-volume Laplacian is used for the diffusion term, with native active diffusivities \(D_u=0.0025\) in Goryachev and \(D_u=0.1\) in Wave-pinning. The plotted gains are

$$
g_u=\frac{\mathcal E_u}{\mathcal D(0)},
\quad
g_v=\frac{\mathcal E_v}{\mathcal D(0)},
\quad
g_D=\frac{\mathcal E_D}{\mathcal D(0)},
\qquad
g_{\mathrm{net}}=g_u+g_v+g_D
=-\frac{\mathcal D'(0^+)}{\mathcal D(0)}.
\tag{6}
$$

Positive net gain describes instantaneous repair of this measured original-core deficit. The individual terms are an algebraic allocation at coupled states; they are not separately evolved interventions. The plotted percentage loss uses the time-dependent matched-sham core amount in its denominator, whereas gains use the deficit at return.

After exact return, global conservation gives

$$
\int_0^L\!\left(v^p-v^s\right)\,\mathrm dx
=\int_0^L\!\left(u^s-u^p\right)\,\mathrm dx.
\tag{7}
$$

The surplus is relative to sham and is balanced by globally lost active material. The reservoir is empty after return, and its original contents are exactly those withheld before the hold. The core deficit is only part of that global active loss; it need not equal the global inactive surplus.

## Bounded reanalysis and limits

All 36 exact return-state chemical, diffusion and total rates were reconstructed directly from native fields and checked against the archived summary. The ordered split was recomputed independently, compared to review R2's numerical table, and checked for pointwise closure and returned-state global conservation at all lengths. Both phases of all six displayed perturbed/sham history pairs were checked against recorded original-core observables. No solver or integration routine was called. Full basal profiles and all retained recovery samples through 5000 native time units are exported; no new convergence claim is made for the original 201-cell meshes.

The two models have separate native units. The same fractional inactive removal produces different active damage and different field perturbations. Their nonlinear activation laws explain how inactive substrate is converted at the actual states; this experiment does not independently vary architecture or measure a general capture capacity. In particular, Goryachev's small negative active-field term rules out interpreting this rebound as direct active-field restoration at unchanged inactive concentration under the stated allocation. Successful competition against a rival, global stability and migration persistence are separate outcomes.
"""


def main():
    retained.register(Path(__file__))
    retained.register(Path(retained.__file__))
    retained.register(Path(__file__).with_name("length_dev_output_layout.py"))
    records, profiles, curves, audit, r3_factors = collect()
    outputs = [write_csv("length_measurements.csv", records),
               write_csv("basal_profiles.csv", profiles),
               write_csv("core_damage_and_recovery.csv", curves),
               write_csv("R3_capture_factorization.csv", r3_factors)]
    rendered, layout = render(records, profiles, curves)
    outputs += rendered
    audit["layout"] = layout
    outputs.append(write_json("verification.json", audit))
    caption = dev_output_path(OUT, "caption.txt")
    caption.write_text(CAPTION + "\n")
    outputs.append(caption)
    math_caption = dev_output_path(OUT, "caption.md")
    r3_link = os.path.relpath(REVIEW_ROOT / "incumbent_r3.md", math_caption.parent)
    math_caption.write_text(MATH_METHODS + "\n## R3 factorization evidence\n\n"
                            + "The native coefficient $A(u)=\\partial f/\\partial v$ is the local response of activation to extra inactive concentration at fixed active concentration. "
                            + "The surplus-weighted coefficient averages $A(u^p)$ using the actual positive core inactive surplus $\\Delta v$ as weights; "
                            + "its product with core inactive surplus divided by core active deficit gives the plotted substrate allocation. "
                            + f"[Fresh review R3, equation (1)]({r3_link}) gives the exact definition; "
                            + "[the imported factorization data](../data/R3_capture_factorization.csv) preserve its four endpoint rows.\n")
    outputs.append(math_caption)
    outputs.append(write_json("methods.json", dict(
        new_integrations=0, source_figures=[3, 4],
        protocol=retained.read(retained.POOL / "provenance.json")["definition"],
        core="Original boundary-connected half-contrast core; its interpolated boundary is fixed, with a fractional edge-cell weight.",
        damage="Integral over original core of (sham active minus perturbed active) at return.",
        excess_reaction="Integral over original core of f(perturbed returned state) minus f(sham returned state).",
        chemical_allocation="Active first at sham inactive: [A(up)-A(us)]vs-c(up-us); inactive second at perturbed active: A(up)(vp-vs). Cross term belongs to inactive contribution. Each integral divided by core active loss at return.",
        excess_diffusion="Integral over original core of Du times (native no-flux FV Laplacian of perturbed active minus sham active).",
        immediate_net_gain="(Excess reaction plus excess diffusion) / damage at return. Positive means instantaneous repair of this measured core deficit.",
        returned_global_conservation="Global inactive surplus relative to sham equals global active loss after returning exactly the originally withheld profile. No external oversupply.",
        trajectory_deficit="100 times integral(sham active minus perturbed active) / integral(sham active), at each saved sample.",
        displayed_time_window="−0.1 to the final saved native sample at or before 1 after return; linear axes. All retained 5000-unit recovery samples exported.",
        profile_view="Goryachev x=0–3 incumbent-end view; Wave-pinning full domains. Full basal profiles are exported.",
        units="Every length/time/concentration/rate is model-native; cross-model absolute rates are not calibrated.",
        causal_limit="Ordered chemical allocation describes coupled returned states; architecture and substrate availability are not independently manipulated, and no rival competition or stability spectrum is evolved.")))
    outputs.append(write_json("review_response.json", dict(
        revision=2, review="current manuscript figs/reviews/incumbent_r2.md",
        previous_version="../archive_main_candidate_1_v1_20261017/",
        interpretation="The previous net repair diagnostic is real, but in Goryachev its positive chemical part is allocated to endogenous inactive-surplus capture, while the active-field-at-sham-inactive contribution is mildly negative.",
        resulting_change="Replace reaction-versus-net verification overplot with ordered active-field, inactive-surplus and diffusion contributions plus net. State native A(u) laws and returned-state global conservation. Frame message as conditional native substrate capture and remote availability.",
        R3_review="current manuscript figs/reviews/incumbent_r3.md",
        R3_response="R3 accepts the six-panel candidate. Add only its supplied coefficient-times-surplus fold factorization to caption/discussion and link the four imported endpoint rows. Clarify that Wave-pinning's sampled concentration changes little; do not infer deep saturation. No figure row or panel substitution.",
        allocation_limit="Exact algebra at actual coupled fields; not evolved pathway ablations or an independently measured general capture capacity.",
        new_integrations=0)))
    endpoints = {name: [r for r in records if r["model"] == name and r["length_native"] in (9, 39)] for name in LABELS}
    outputs.append(write_json("key_results.json", endpoints))
    manifest = dict(generated_at=datetime.now(timezone.utc).isoformat(), producer="finished",
                    status="Provisional main-text candidate 1, revision 2; not final story selection", revision=2, artifact_root="../..",
                    previous_version="../archive_main_candidate_1_v1_20261017/",
                    new_integrations=0, renderer=str(Path(__file__).relative_to(ROOT)),
                    source_sha256=retained.SOURCES,
                    artifacts={dev_artifact_key(OUT, path): retained.sha(path) for path in outputs},
                    reproduction="PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python " + str(Path(__file__)))
    write_json("manifest.json", manifest)
    print(json.dumps(dict(output=str(OUT), passed=True, new_integrations=0,
                          sources=len(retained.SOURCES), net_ratios={name: rows[-1]["net_gain_relative_to_L9"] for name, rows in endpoints.items()}), indent=2))


if __name__ == "__main__":
    main()
