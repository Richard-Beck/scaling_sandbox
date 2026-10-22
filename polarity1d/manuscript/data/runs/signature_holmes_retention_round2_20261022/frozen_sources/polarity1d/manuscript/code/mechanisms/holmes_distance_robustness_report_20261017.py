"""Read-only reconstruction and compact figures for bounded Holmes checks."""
from pathlib import Path
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
from scipy.signal import find_peaks, peak_widths
import holmes_distance_robustness_20261017 as campaign

ROOT, RUN, OLD = campaign.ROOT, campaign.RUN, campaign.OLD
OUT = ROOT / "figures/length_mechanisms_20261008/dev/holmes_distance_robustness_20261017"
SCALES = campaign.SCALES
NAMES = ["Active Rac", "Inactive Rac", "Active Rho", "Inactive Rho"]
COLORS = ["#343434", "#4575b4", "#d76e2e"]


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save(fig, name):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"{name}.{suffix}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def peaks(grid, state, origin, label):
    u = state[0]; tiled = np.tile(u, 3)
    cells, props = find_peaks(tiled, prominence=0.05)
    keep = (cells >= grid.cells) & (cells < 2 * grid.cells)
    cells = cells[keep]
    width, level, left, right = peak_widths(tiled, cells, rel_height=0.5)
    rows = []
    for i, cell in enumerate(cells):
        rows.append(dict(case=label, cell=int(cell % grid.cells),
            position=float(grid.arc(grid.x[cell % grid.cells], origin)),
            height=float(u[cell % grid.cells]), prominence=float(props["prominences"][keep][i]),
            half_prominence_width=float(width[i] * grid.dx),
            left_valley=float(tiled[props["left_bases"][keep][i]]),
            right_valley=float(tiled[props["right_bases"][keep][i]])))
    return rows


def load_phase(path, expected_initial=None, end=None):
    with np.load(path) as z:
        t, states = z["time"].copy(), z["states"].copy()
        stimulus = z["stimulus"].copy()
    assert np.isfinite(states).all() and states.min() >= -5e-10
    assert np.all(np.diff(t) > 0)
    if end is not None:
        assert t[-1] == end
    if expected_initial is not None:
        np.testing.assert_array_equal(states[0], expected_initial)
    grid = campaign.ring.RingGrid(40, states.shape[-1])
    pools = np.column_stack([states[:, group, :].sum(axis=(1, 2)) * grid.dx for group in ([0, 1], [2, 3])])
    pool_error = float(np.max(np.abs(pools / pools[0] - 1)))
    assert pool_error < 1e-7
    return t, states, stimulus, pool_error


def imbalance():
    manifest = campaign.check_sources(RUN)
    base = RUN / "inputs"
    with np.load(base / "imbalance_0.npz") as z:
        reference = z["reference"].copy(); weights = z["fractional_core_weights"].copy()
    grid = campaign.ring.RingGrid(40, 800)
    x = grid.arc(grid.x, .025 + 8) + 8; order = np.argsort(x)
    bg = manifest["perturbation"]["background"]
    records, summaries, traces, peakrows, profrows = [], [], [], [], []
    maximum_mass = 0.0
    for slot, sign in enumerate((0, -1, 1)):
        status = campaign.read(RUN / f"workers/imbalance_{slot}/status.json")
        assert status["status"] == "complete"
        path = base / f"imbalance_{slot}.npz"
        with np.load(path) as z:
            initial = z["state"].copy(); delta = z["delta"].copy()
            smooth = z["smooth_weights"].copy()
        np.testing.assert_array_equal(initial[1:], reference[1:])
        np.testing.assert_array_equal(initial[0], reference[0] + delta)
        assert initial.min() >= 0
        mass = float(delta.sum() * grid.dx)
        maximum_mass = max(maximum_mass, abs(mass))
        assert abs(mass) < 1e-13
        np.testing.assert_allclose(smooth.sum(axis=1) * grid.dx, [1, 1], rtol=0, atol=1e-14)
        phasepath = RUN / f"workers/imbalance_{slot}/continuation.npz"
        assert campaign.sha(phasepath) == status["phases"][0]["sha256"]
        t, states, signal, pool_error = load_phase(phasepath, initial, 13000)
        np.testing.assert_array_equal(signal, np.zeros(grid.cells))
        amounts = states[:, 0] @ weights.T
        excess = (states[:, 0] - bg) @ weights.T
        imbalance = (excess[:, 1] - excess[:, 0]) / excess.sum(axis=1)
        records.append(dict(sign=sign, time=t, states=states, delta=delta, imbalance=imbalance))
        for i, ti in enumerate(t):
            traces.append(dict(sign=sign, continuation_time=ti,
                old_fixed_core_amount=amounts[i, 0], cue_fixed_core_amount=amounts[i, 1],
                old_fixed_core_excess=excess[i, 0], cue_fixed_core_excess=excess[i, 1],
                signed_fixed_core_excess_imbalance=imbalance[i]))
        peakrows += peaks(grid, states[-1], .025, f"sign{sign:+d}_terminal")
        summaries.append(dict(sign=sign, initial_delta_integral=mass,
            moved_amount=float(np.maximum(delta, 0).sum() * grid.dx),
            maximum_initial_delta=float(np.max(np.abs(delta))),
            initial_minimum=float(initial.min()), pool_relative_error=pool_error,
            initial_imbalance=float(imbalance[0]), final_imbalance=float(imbalance[-1]),
            final_resolved_peaks=len(peaks(grid, states[-1], .025, ""))))
    unperturbed = records[0]
    for record, summary in zip(records, summaries):
        differences = record["states"][-1] - unperturbed["states"][-1]
        summary["final_max_scaled_full_field_difference_from_control"] = float(np.max(np.abs(differences) / SCALES[:, None]))
        summary["initial_imbalance_difference_from_control"] = float(record["imbalance"][0] - unperturbed["imbalance"][0])
        summary["final_imbalance_difference_from_control"] = float(record["imbalance"][-1] - unperturbed["imbalance"][-1])
        for j in range(grid.cells):
            for k, name in enumerate(NAMES):
                profrows.append(dict(sign=record["sign"], arc_from_original_peak=x[j], field=name,
                    initial=record["states"][0, k, j], terminal=record["states"][-1, k, j],
                    terminal_difference_from_control=differences[k, j]))
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.5), constrained_layout=True)
    for record, color in zip(records, COLORS):
        label = {0: "Unperturbed", -1: "Cue → old", 1: "Old → cue"}[record["sign"]]
        axes[0, 0].plot(x[order], record["delta"][order], color=color, label=label)
        axes[0, 1].plot(x[order], record["states"][-1, 0, order], color=color, label=label)
        axes[1, 0].plot(record["time"], record["imbalance"] - unperturbed["imbalance"], color=color, label=label)
        axes[1, 1].plot(x[order], (record["states"][-1, 0] - unperturbed["states"][-1, 0])[order], color=color)
    axes[0, 0].set(title="A  Initial active-Rac redistribution", ylabel="Δ active Rac")
    axes[0, 1].set(title="B  Two-site endpoint after 13,000 units", ylabel="Active Rac")
    axes[1, 0].set(title="C  Fixed-core imbalance relative to control", xlabel="Continuation time", ylabel="Δ signed excess imbalance")
    axes[1, 0].set_xscale("symlog", linthresh=1)
    axes[1, 0].set_xticks([0, 1, 10, 100, 1000, 13000], ["0", "1", "10", "100", "1000", "13000"])
    axes[1, 0].set_xlabel("Continuation time (linear to 1, log thereafter)")
    axes[1, 1].set(title="D  Terminal active-Rac differences", ylabel="Active Rac − control")
    axes[0, 0].legend(frameon=False, fontsize=9)
    for ax in (axes[0, 0], axes[0, 1], axes[1, 1]):
        ax.set_xlabel("Arc from original incumbent peak")
        ax.axvline(0, color=".7", lw=.7, ls=":"); ax.axvline(16, color=".7", lw=.7, ls=":")
    for ax in axes.flat:
        ax.spines[["top", "right"]].set_visible(False)
    save(fig, "imbalance_continuation")
    write_csv(OUT / "supporting/data/imbalance_summary.csv", summaries)
    write_csv(OUT / "supporting/data/imbalance_traces.csv", traces)
    write_csv(OUT / "supporting/data/imbalance_terminal_peaks.csv", peakrows)
    write_csv(OUT / "supporting/data/imbalance_full_field_profiles.csv", profrows)
    campaign.write(OUT / "supporting/metadata/imbalance_verification.json", dict(generated_at=campaign.now(), passed=True,
        exact_phase_start_and_native_fields_verified=True, maximum_abs_initial_active_Rac_mass_change=maximum_mass,
        all_fields_positive=True, native_pool_means=[7.5, 3.1], summaries=summaries))
    text = f"""# Holmes two-site pattern after conserved active-Rac imbalance

Generated {campaign.now()}; producer finished. All three bounded SLURM continuations are complete.

The two resolved Rac sites remain after a further 13,000 cue-free time units in the unperturbed continuation and both signs of the specified 1% redistribution. At the endpoint, about {abs(summaries[1]['final_imbalance_difference_from_control']/summaries[1]['initial_imbalance_difference_from_control']):.3%} of the initial fixed-core excess imbalance difference from the unperturbed control remains. Most reduction occurs during the first 100 units; the absolute residual then increases slightly from approximately 2.10 × 10⁻⁴ at 100 to 2.11 × 10⁻⁴ at 13,000. Thus the recorded trajectory shows rapid initial relaxation followed by a small, slowly increasing residual, rather than sustained damping. This is finite-window robustness to this perturbation, not an eigenvalue, proof of branch stability, or exclusion of slower competition.

![Conserved Rac imbalance and continuation](../../imbalance_continuation.png)

The figure shows the exact initial active-Rac change, raw terminal active-Rac profiles, the signed fixed-core excess imbalance relative to the concurrently evolved unperturbed control, and the small remaining profile difference. Panel C uses a time axis linear through 1 and logarithmic thereafter, exposing the rapid initial change and the later window together. The two site positions and entire profiles remain available; fixed original-core activity is a spatial readout, not a survival criterion. The control itself continues to change slightly.

All cases start from the exact N800 d16 endpoint after its original 13,000-unit basal window. The old and cue terminal local maxima are at −0.10 and +16.00 relative to the original incumbent. Each fixed site core is its interpolated half-prominence interval; the widths are 8.52355 and 7.92087. The fixed background is the starting profile minimum, {bg:.12g}. The smaller core's above-background active-Rac amount is 8.16912328; 1% is {manifest['perturbation']['epsilon']:.12g}. This transferred amount is only 0.0272304% of the conserved total Rac pool of 300; the experiment does not transfer 1% of total Rac. Smooth raised-cosine weights taper to zero at each core boundary, and each has unit integrated weight. The same amount is removed from one site and added to the other in active Rac. All other fields are identical at the initial instant; fields remain nonnegative. Maximum absolute integrated active-Rac change is {maximum_mass:.3g}. Native conserved Rac and Rho pools, kinetics, diffusivities and solver remain unchanged.

The signed excess imbalance is (cue excess − old excess)/(cue excess + old excess), integrating each concentration minus the fixed starting background over its fixed fractional-cell half-prominence core. The cores differ in width, so this is a defined local comparison, not a symmetry coordinate or conserved-pool partition. Global pool conservation is checked separately on every saved full state. The ± labels denote equal and opposite initial redistributions, not biological replicates.

The initial and terminal imbalances are {summaries[0]['initial_imbalance']:.9f} → {summaries[0]['final_imbalance']:.9f} (control), {summaries[1]['initial_imbalance']:.9f} → {summaries[1]['final_imbalance']:.9f} (cue to old), and {summaries[2]['initial_imbalance']:.9f} → {summaries[2]['final_imbalance']:.9f} (old to cue). Maximum sampled pool relative error is {max(r['pool_relative_error'] for r in summaries):.3g}. The remaining scaled full-field differences from control are {summaries[1]['final_max_scaled_full_field_difference_from_control']:.6g} and {summaries[2]['final_max_scaled_full_field_difference_from_control']:.6g}, normalized by frozen native field scales 7.5, 7.5, 3.1, 3.1.

This test challenges the near-symmetric two-site tail in a controlled direction while preserving the global Rac amount. It cannot establish asymptotic stability, arbitrary-perturbation robustness, a general separation law, or biological migration persistence. It does not change or reinterpret the native interval model.

Data: [exact perturbation manifest](../../../../../../data/runs/holmes_distance_robustness_20261017/manifest.json), [frozen inputs and exact deltas](../../../../../../data/runs/holmes_distance_robustness_20261017/inputs), [full trajectory archives](../../../../../../data/runs/holmes_distance_robustness_20261017/workers), [summary](../data/imbalance_summary.csv), [fixed-core traces](../data/imbalance_traces.csv), [terminal peaks/valleys/widths](../data/imbalance_terminal_peaks.csv), [all four initial/terminal fields](../data/imbalance_full_field_profiles.csv), and [verification](imbalance_verification.json). Array 21083244 used three one-CPU/3-GB QoS-small workers without a concurrency throttle. The unchanged native BDF tolerances are rtol 1e−9, atol 1e−11, maximum step 10; sampling is 0.1 through 10, 1 through 100, and 5 thereafter.
"""
    (OUT / "supporting/metadata/IMBALANCE_RESULTS.md").write_text(text)
    return summaries


def mesh():
    if not all((RUN / f"workers/fine_{slot}/status.json").exists() and campaign.read(RUN / f"workers/fine_{slot}/status.json")["complete"] for slot in (0, 1)):
        return None
    with np.load(RUN / "inputs/fine_basal.npz") as z:
        fine_basal = z["state"].copy()
    with np.load(RUN / "inputs/original_basal.npz") as z:
        coarse_basal = z["state"].copy()
    grid = campaign.ring.RingGrid(40, 800)
    x = grid.arc(grid.x, grid.x[0]); order = np.argsort(x)
    snapshot_rows, error_rows, peakrows, comparison = [], [], [], {}
    transient_profiles = None
    basename = campaign.sha(RUN / "inputs/fine_basal.npz")
    for slot, d, oldslot in ((0, 12, 10), (1, 16, 11)):
        status = campaign.read(RUN / f"workers/fine_{slot}/status.json")
        assert status["status"] == "complete" and status["input_sha256"] == basename
        previous = {800: coarse_basal, 1600: fine_basal}
        saved = {800: {"basal": coarse_basal}, 1600: {"basal": fine_basal}}
        for phase, end in (("cue", 1000), ("recovery", 13000)):
            phases = []
            for n, path in ((800, OLD / f"workers/{oldslot:02d}/{phase}.npz"),
                            (1600, RUN / f"workers/fine_{slot}/{phase}.npz")):
                t, y, signal, err = load_phase(path, previous[n], end)
                expected_cue = campaign.ring.ring_cue(campaign.ring.RingGrid(40, n), int(d / (40 / n)), 3, 1) if phase == "cue" else np.zeros(n)
                np.testing.assert_array_equal(signal, expected_cue)
                phases.append((t, y, err)); previous[n] = y[-1].copy()
                saved[n][phase] = y[-1].copy()
                if n == 1600:
                    info = next(p for p in status["phases"] if p["phase"] == phase)
                    assert campaign.sha(path) == info["sha256"]
            np.testing.assert_array_equal(phases[0][0], phases[1][0])
            # Incumbent-centred coordinates coincide at every second fine cell.
            difference = phases[1][1][:, :, ::2] - phases[0][1]
            for field, name in enumerate(NAMES):
                ti, cell = np.unravel_index(np.argmax(np.abs(difference[:, field])), difference[:, field].shape)
                error_rows.append(dict(distance=d, phase=phase, field=name,
                    max_abs_difference_all_sampled_times=float(np.max(np.abs(difference[:, field]))),
                    max_scaled_difference_all_sampled_times=float(np.max(np.abs(difference[:, field])) / SCALES[field]),
                    time_of_max_difference=float(phases[0][0][ti]),
                    arc_of_max_difference=float(x[cell]),
                    terminal_max_abs_difference=float(np.max(np.abs(difference[-1, field]))),
                    coarse_pool_relative_error=phases[0][2], fine_pool_relative_error=phases[1][2]))
                if d == 12 and phase == "recovery" and field == 2:
                    transient_profiles = (float(phases[0][0][ti]), phases[0][1][ti].copy(), phases[1][1][ti].copy())
        comparison[d] = saved
        for n in (800, 1600):
            g = campaign.ring.RingGrid(40, n)
            xx = g.arc(g.x, g.x[0])
            for phase, state in saved[n].items():
                peakrows += peaks(g, state, g.x[0], f"d{d}_N{n}_{phase}")
                for cell in range(n):
                    row = dict(distance=d, cells=n, phase=phase, arc_from_incumbent_peak=xx[cell])
                    row.update({name: state[k, cell] for k, name in enumerate(NAMES)})
                    snapshot_rows.append(row)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.4), constrained_layout=True)
    for row, d in enumerate((12, 16)):
        for col, phase in enumerate(("cue", "recovery")):
            ax = axes[row, col]
            for n, style, color in ((800, "--", "#c46a3a"), (1600, "-", "#244c75")):
                g = campaign.ring.RingGrid(40, n); xx = g.arc(g.x, g.x[0] + d/2) + d/2; sort = np.argsort(xx)
                ax.plot(xx[sort], comparison[d][n][phase][0, sort], ls=style, color=color, lw=1.7, label=f"N = {n}")
            ax.axvline(0, color=".65", ls=":", lw=.7); ax.axvline(d, color=".65", ls=":", lw=.7)
            ax.set(title=f"d = {d}: " + ("cue withdrawal" if phase == "cue" else "13,000 units after withdrawal"),
                ylabel="Active Rac", xlabel="Arc from original incumbent peak")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False)
    save(fig, "challenge_mesh_profiles")
    fig, axes = plt.subplots(2, 4, figsize=(13, 5.5), constrained_layout=True)
    for row, d in enumerate((12, 16)):
        for field in range(4):
            ax = axes[row, field]
            for n, style, color in ((800, "--", "#c46a3a"), (1600, "-", "#244c75")):
                g = campaign.ring.RingGrid(40, n); xx = g.arc(g.x, g.x[0] + d/2) + d/2; sort = np.argsort(xx)
                ax.plot(xx[sort], comparison[d][n]["recovery"][field, sort], ls=style, color=color, lw=1.4, label=f"N{n}")
            ax.set(title=f"d={d}, {NAMES[field]}", xlabel="Incumbent-centred arc")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False)
    save(fig, "challenge_mesh_all_fields")
    ti, coarse_transient, fine_transient = transient_profiles
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3), constrained_layout=True)
    for ax, field in zip(axes, (0, 2)):
        for n, state, style, color in ((800, coarse_transient, "--", "#c46a3a"), (1600, fine_transient, "-", "#244c75")):
            g = campaign.ring.RingGrid(40, n); xx = g.arc(g.x, g.x[0] + 6) + 6; sort = np.argsort(xx)
            ax.plot(xx[sort], state[field, sort], ls=style, color=color, label=f"N = {n}")
        ax.set(title=f"d12, {NAMES[field]}, {ti:g} units after withdrawal", xlabel="Arc from original incumbent peak", ylabel=NAMES[field])
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False)
    save(fig, "challenge_mesh_largest_transient")
    write_csv(OUT / "supporting/data/mesh_field_errors.csv", error_rows)
    write_csv(OUT / "supporting/data/mesh_peak_metrics.csv", peakrows)
    write_csv(OUT / "supporting/data/mesh_full_field_profiles.csv", snapshot_rows)
    maxerr = max(row["max_scaled_difference_all_sampled_times"] for row in error_rows)
    finalerr = max(row["terminal_max_abs_difference"] / SCALES[NAMES.index(row["field"])] for row in error_rows if row["phase"] == "recovery")
    campaign.write(OUT / "supporting/metadata/mesh_verification.json", dict(generated_at=campaign.now(), passed=True,
        shared_fine_checkpoint_sha256=basename, exact_phase_continuity=True, exact_cues_and_sampling=True,
        maximum_scaled_full_field_difference_all_saved_samples=maxerr, maximum_scaled_terminal_difference=finalerr,
        spatial_comparison="Every second fine cell matches coarse incumbent-centred coordinate; no optimized alignment, remapping, fitted shift, or normalization of concentrations", errors=error_rows))
    final = [p for p in peakrows if p["case"].endswith("recovery")]
    table = "\n".join(f"| {p['case']} | {p['position']:.3f} | {p['height']:.7f} | {p['half_prominence_width']:.5f} | {p['left_valley']:.7f}, {p['right_valley']:.7f} |" for p in final)
    text = f"""# Holmes distance contrast persists on the finer mesh

Generated {campaign.now()}; producer finished. Exactly two fresh N1600 challenge trajectories were integrated.

The d12 endpoint remains a single Rac cap between the original incumbent and cue, while d16 retains two separated Rac caps through the recorded 13,000-unit cue-free window. N800 and N1600 agree in these full spatial profiles, rather than only in original-core means or component counts. The largest scaled field difference across all saved cue/recovery samples is {maxerr:.6g}; the largest terminal scaled difference is {finalerr:.6g}, using frozen field scales 7.5, 7.5, 3.1, 3.1.

![Challenge mesh check](../../challenge_mesh_profiles.png)

These four raw active-Rac panels compare cue withdrawal and the terminal cue-free endpoint for both tested distances. Dashed orange is N800 and solid blue N1600. Dotted lines mark the original incumbent and cue centres. The circumference is 40 in both meshes. Width, concentration, peak positions and valleys are visible; no threshold replaces the profiles.

![All four terminal fields](../../challenge_mesh_all_fields.png)

The largest sampled mesh residual occurs during the d12 merger, at {ti:g} time units after withdrawal, in active Rho. Its absolute magnitude is {max(r['max_abs_difference_all_sampled_times'] for r in error_rows if r['distance']==12 and r['phase']=='recovery' and r['field']=='Active Rho'):.6g}. The raw profiles below expose that transient difference; the endpoint agreement alone does not establish identical merger timing or uniform transient accuracy.

![Largest sampled transient mesh discrepancy](../../challenge_mesh_largest_transient.png)

| Case | Peak position | Height | Half-prominence width | Left, right valleys |
|---|---:|---:|---:|---|
{table}

The N1600 basal state was prepared afresh from the native homogeneous state, with the same peak 1, total-width 3 raised-cosine cue for exactly 10 time units, then the unchanged full-field settling rule and additional 10,000-unit holdout. It qualified after 13,000 cue-free units. Both fine challenges restart the exact same saved qualified checkpoint `{basename}`; coarse and fine states are separately prepared native ring states, not interpolated initial conditions. Native Rac and Rho pool means remain 7.5 and 3.1. The fine preparation saved-field audit verifies exact continuity, positive fields and maximum pool relative error 4.44e−16. Its holdout excursion is 3.58e−16 on frozen scales.

The challenger is the same total-width 3, peak 1 symmetric raised cosine for 1000 units followed by 13,000 basal units. Only location changes within each mesh: d12 and d16 are exactly 480 and 640 cells at N1600, relative to the measured incumbent peak. Native reactions, diffusivities, conserved groups, BDF tolerances and sample times are unchanged. Comparisons use physical arc relative to each mesh's original incumbent centre; every second fine cell shares the coarse coordinate. No fitted translation or profile normalization is used. [All-time all-field residuals](../data/mesh_field_errors.csv) include every retained time sample. This finite-difference check supports resolution of the selected contrast; two meshes do not establish an asymptotic convergence order.

This is a bounded mesh challenge to one within-model observation. The finite d16 two-site window does not establish a stable branch or rule out slower competition. These two distances do not establish a transition location, a monotone distance law, minimum cue durations, or biological WGD/migration calibration.

Data: [full profile snapshots](../data/mesh_full_field_profiles.csv), [all peak/valley/width values](../data/mesh_peak_metrics.csv), [all-field errors](../data/mesh_field_errors.csv), [verification](mesh_verification.json), [fine preparation audit](../../../../../../data/runs/holmes_distance_robustness_20261017/fine_preparation/results/independent_audit.json), [fine challenge archives](../../../../../../data/runs/holmes_distance_robustness_20261017/workers), [frozen campaign](../../../../../../data/runs/holmes_distance_robustness_20261017/manifest.json). Preparation job 21083243 and challenge array 21087361 used QoS small, one CPU and 3 GB per worker, without an array throttle. The completed preparation wrapper omitted its nested launch record, so its automatic audit/handoff stopped after successful integration. The missing record was explicitly reconstructed from the original wrapper submission and unchanged manifest, the saved-field audit passed, then the two challenge workers were submitted. No frozen source or trajectory was changed.
"""
    (OUT / "supporting/metadata/MESH_RESULTS.md").write_text(text)
    return dict(maximum_scaled_full_field_difference=maxerr, terminal_scaled_difference=finalerr)


def main():
    for name in ("supporting/data", "supporting/metadata"):
        (OUT / name).mkdir(parents=True, exist_ok=True)
    perturb = imbalance(); fine = mesh()
    (OUT / "supporting/metadata/RESULTS.md").write_text("# Bounded Holmes distance verification\n\n" +
        f"Generated {campaign.now()}.\n\n" +
        "[Conserved active-Rac imbalance result](IMBALANCE_RESULTS.md) includes the completed three-worker continuation and exact initial redistribution.\n\n" +
        ("[N800/N1600 challenge mesh result](MESH_RESULTS.md) includes the completed fresh fine preparation and exactly d12/d16.\n" if fine else "The fresh N1600 challenges are still running; the fine-mesh result is pending.\n") +
        ("\n[Fresh review and response](REVIEW_RESPONSE.md).\n" if (OUT / "supporting/metadata/REVIEW_RESPONSE.md").exists() else ""))
    artifacts = [p for p in OUT.rglob("*") if p.is_file() and p.name != "manifest.json"]
    campaign.write(OUT / "supporting/metadata/manifest.json", dict(generated_at=campaign.now(),
        producer="finished" if fine else "fine mesh pending", source_sha256={str(Path(__file__).relative_to(ROOT)): campaign.sha(Path(__file__))},
        artifacts={str(p.relative_to(OUT)): campaign.sha(p) for p in artifacts},
        run_manifest_sha256=campaign.sha(RUN / "manifest.json"),
        reproduction=f"{campaign.PYTHON} {Path(__file__)}"))
    print(json.dumps(dict(output=str(OUT), imbalance=perturb, fine=fine), indent=2))


if __name__ == "__main__":
    main()
