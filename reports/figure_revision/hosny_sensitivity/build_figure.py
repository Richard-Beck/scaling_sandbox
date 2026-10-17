"""Render the standalone Hosny sensitivity review figure from analysis outputs."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports/figure_revision/figures/figure3_review/hosny_sensitivity"
SURFACE = OUT / "hosny_surface.csv"
SURFACE_ZOOM = OUT / "hosny_surface_zoom.csv"
SURFACE_ZOOM_META = OUT / "hosny_surface_zoom.json"
PROFILES = OUT / "hosny_drift_profiles.csv"
LEAVEOUT = OUT / "hosny_leaveout_comparison.csv"
ANALYSIS = OUT / "hosny_analysis.json"
FIT_PATH = ROOT / "wgd_langevin/gaussian_hierarchy/results/fit_05.json"
PNG = OUT / "hosny_sensitivity.png"

NAVY = "#19394d"
FIXED = "#168b87"
PROFILED = "#d87529"
LEVEL_COLORS = {1: "#25855a", 4: "#c48216", 9: "#b3444b"}
FIXED_LEVEL_COLORS = {1: "#9dcdb0", 4: "#ddc580", 9: "#dca1a5"}
STUDIES = ["Xuan2018", "Mallin2023", "Jemaa2017", "Jemaa2023", "Pooled"]
ANALYSES = ["full5_allstudies", "full5_remaining4_only", "leaveHosnyOut_refit"]
ANALYSIS_LABELS = ["5-study fit", "Same fit, pool 4", "Refit on 4"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.titleweight": "bold", "axes.spines.top": False,
                     "axes.spines.right": False})


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value in (None, ""):
        return float("nan")
    return float(value)


def is_true(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def surface_panel(ax, local_records, full_records, score_prefix: str, title: str, caption: str,
                  fit_tau_ratio: float, fit_sigma_ratio: float):
    tau = np.array(sorted({number(r, "tau_ratio") for r in local_records}))
    sigma = np.array(sorted({number(r, "sigma_ratio") for r in local_records}))
    xgrid, ygrid = np.meshgrid(tau, sigma, indexing="ij")
    methods = [
        (f"fixed_drift_delta_{score_prefix}", "fixed", "-"),
        (f"optimized_drift_delta_{score_prefix}", "refit", "--")]
    for key, method, style in methods:
        values = np.full((len(tau), len(sigma)), np.nan)
        lookup = {(number(r, "tau_ratio"), number(r, "sigma_ratio")): number(r, key)
                  for r in local_records}
        for i, tx in enumerate(tau):
            for j, sy in enumerate(sigma):
                values[i, j] = lookup[(tx, sy)]
        colors = [FIXED_LEVEL_COLORS[v] if method == "fixed" else LEVEL_COLORS[v]
                  for v in (1, 4, 9)]
        ax.contour(xgrid, ygrid, values, levels=[1, 4, 9], colors=colors,
                   linestyles=style, linewidths=2.45 if method == "fixed" else 1.2,
                   corner_mask=True, zorder=2 if method == "fixed" else 3)
    ax.scatter([fit_tau_ratio], [fit_sigma_ratio], marker="*", s=105, facecolor="#111111",
               edgecolor="white", linewidth=.65, zorder=6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("WGD+ / WGD− persistence (τ) fold")
    ax.set_ylabel("WGD+ / WGD− noise (σ) fold")
    ax.set_title(title, loc="left", color=NAVY, pad=5)
    ax.grid(which="major", color="#dce4e8", linewidth=.55, alpha=.8)
    ax.text(.02, .03, caption, transform=ax.transAxes, va="bottom", ha="left",
            fontsize=7, color="#526570",
            bbox=dict(boxstyle="round,pad=.25", facecolor="white", edgecolor="none", alpha=.88))
    ax.set_xlim(float(tau[0]), float(tau[-1]))
    ax.set_ylim(float(sigma[0]), float(sigma[-1]))

    # Keep the full requested surface in view while enlarging the fitted valley.
    ix = inset_axes(ax, width="28%", height="28%", loc="upper right", borderpad=.8)
    full_tau = np.array(sorted({number(r, "tau_ratio") for r in full_records}))
    full_sigma = np.array(sorted({number(r, "sigma_ratio") for r in full_records}))
    full_xgrid, full_ygrid = np.meshgrid(full_tau, full_sigma, indexing="ij")
    for key, method, style in methods:
        values = np.full((len(full_tau), len(full_sigma)), np.nan)
        lookup = {(number(r, "tau_ratio"), number(r, "sigma_ratio")): number(r, key)
                  for r in full_records}
        for i, tx in enumerate(full_tau):
            for j, sy in enumerate(full_sigma):
                values[i, j] = lookup[(tx, sy)]
        colors = [FIXED_LEVEL_COLORS[v] if method == "fixed" else LEVEL_COLORS[v]
                  for v in (1, 4, 9)]
        ix.contour(full_xgrid, full_ygrid, values, levels=[1, 4, 9], colors=colors,
                   linestyles=style, linewidths=.9 if method == "fixed" else .55,
                   corner_mask=True, zorder=2 if method == "fixed" else 3)
    ix.set_xscale("log"); ix.set_yscale("log")
    ix.set_xlim(float(full_tau[0]), float(full_tau[-1])); ix.set_ylim(float(full_sigma[0]), float(full_sigma[-1]))
    ix.scatter([1.0], [1.0], marker="o", s=20, facecolor="white", edgecolor="#34434b",
               linewidth=.7, zorder=5)
    ix.scatter([fit_tau_ratio], [fit_sigma_ratio], marker="*", s=38, facecolor="#111111",
               edgecolor="white", linewidth=.45, zorder=6)
    ix.set_title("Full", fontsize=6, pad=1)
    ix.tick_params(axis="both", which="both", labelsize=5, length=2)
    ix.grid(which="major", color="white", linewidth=.35, alpha=.7)
    ix.set_facecolor("white")
    for spine in ix.spines.values():
        spine.set_edgecolor("#89979e"); spine.set_linewidth(.55)
    ax.text(.02, .97, "Fixed/refitted drift contours nearly overlap", transform=ax.transAxes, va="top", ha="left",
            fontsize=6.6, color="#526570")


def drift_profile(ax, rows):
    profiles = [
        ("mu_reference", "Hosny WGD− drift", "#4169a1", "o"),
        ("mu_wgd", "Hosny WGD+ drift", "#a34e7b", "s"),
        ("drift_contrast", "WGD+ − WGD− drift", "#37834d", "^"),
    ]
    for name, label, color, marker in profiles:
        group = [r for r in rows if r["profile_parameter"] == name
                 and is_true(r.get("selected_start", False))
                 and is_true(r.get("success", False))]
        group.sort(key=lambda r: number(r, "fixed_value"))
        x = np.array([number(r, "fixed_value") for r in group])
        y = np.array([number(r, "delta_from_best_unrestricted") for r in group])
        ok = np.isfinite(x) & np.isfinite(y)
        ax.plot(x[ok], y[ok], marker=marker, ms=3.4, lw=1.2, color=color,
                label=label)
    for threshold, color in LEVEL_COLORS.items():
        ax.axhline(threshold, color=color, ls=":", lw=.9)
    ax.set_xlabel("Fixed signed drift during profile (µm/h)")
    ax.set_ylabel("Increase in total fit loss")
    ax.set_title("C  Drift profiles", loc="left", color=NAVY, pad=5)
    ax.grid(alpha=.22)
    ax.legend(frameon=False, fontsize=6.5, ncol=1, loc="upper center")


def tradeoff(ax, rows, hosny_tau_fold: float, hosny_sigma_fold: float):
    data = [r for r in rows if number(r, "delta_from_best_unrestricted") <= 9
            and is_true(r.get("selected_start", False))
            and is_true(r.get("success", False))]
    delta = np.array([number(r, "delta_from_best_unrestricted") for r in data])
    tau = np.array([number(r, "hosny_tau_fold") for r in data])
    sigma = np.array([number(r, "hosny_sigma_fold") for r in data])
    valid = np.isfinite(delta) & np.isfinite(tau) & np.isfinite(sigma)
    cmap = plt.get_cmap("viridis")
    points = ax.scatter(tau[valid], sigma[valid], c=delta[valid], cmap=cmap,
                        vmin=0, vmax=9, marker="o", s=38, alpha=.86,
                        edgecolor="white", linewidth=.45, zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.scatter([hosny_tau_fold], [hosny_sigma_fold], marker="o", s=96,
               facecolor="none", edgecolor="#111111", linewidth=1.3, zorder=5)
    ax.annotate("Selected fit", xy=(hosny_tau_fold, hosny_sigma_fold),
                xytext=(8, 10), textcoords="offset points", fontsize=7,
                color="#263943", arrowprops=dict(arrowstyle="-", lw=.7, color="#526570"))
    xlo, xhi = min(float(tau[valid].min()), hosny_tau_fold), max(float(tau[valid].max()), hosny_tau_fold)
    ylo, yhi = min(float(sigma[valid].min()), hosny_sigma_fold), max(float(sigma[valid].max()), hosny_sigma_fold)
    ax.set_xlim(xlo / 1.12, xhi * 1.12)
    ax.set_ylim(ylo / 1.15, yhi * 1.15)
    ax.set_xlabel("Hosny WGD+ / WGD− persistence fold")
    ax.set_ylabel("Hosny WGD+ / WGD− noise fold")
    ax.set_title("D  Hosny τ/noise profile trade-offs", loc="left", color=NAVY, pad=5)
    ax.grid(which="major", alpha=.2)
    cb = plt.colorbar(points, ax=ax, fraction=.05, pad=.03, ticks=[0, 1, 4, 9])
    cb.set_label("Increase in total fit objective (Δ)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.text(.02, .02, "Selected 1D profiles; all points use the same circle",
            transform=ax.transAxes, va="bottom", fontsize=6.8, color="#526570")


def leaveout_heatmap(ax, rows, field: str, letter: str, title: str,
                     full_fit_study_values: dict[str, tuple[float, float]]):
    selected = {r["analysis"]: r for r in rows if r["study"] == "Pooled"}
    matrix = []
    for study in STUDIES:
        vals = []
        for analysis in ANALYSES:
            if study == "Pooled":
                row = selected.get(analysis)
                col = "population_" + field
            else:
                row = next((r for r in rows if r["analysis"] == analysis and r["study"] == study), None)
                col = field
                if analysis == "full5_allstudies" and row is None:
                    value_pair = full_fit_study_values.get(study)
                    vals.append(value_pair[0 if field == "tau_fold" else 1]
                                if value_pair is not None else float("nan"))
                    continue
            vals.append(number(row, col) if row else float("nan"))
        matrix.append(vals)
    mat = np.asarray(matrix, float)
    logfold = np.log2(mat)
    finite = np.abs(logfold[np.isfinite(logfold)])
    limit = max(float(np.max(finite)) if finite.size else 1., .15)
    im = ax.imshow(logfold, cmap="coolwarm", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                   aspect="auto")
    ax.set_xticks(range(3), ANALYSIS_LABELS, rotation=0)
    ax.set_yticks(range(len(STUDIES)), STUDIES)
    ax.set_title(f"{letter}  {title}", loc="left", color=NAVY, pad=5)
    for i in range(len(STUDIES)):
        for j in range(3):
            value = mat[i, j]
            txt = f"{value:.3g}" if np.isfinite(value) else "—"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                    color="#192e3a" if abs(logfold[i, j]) < limit * .72 else "white")
    ax.set_xlabel("Fit / scoring comparison")
    ax.tick_params(axis="x", labelsize=6.5)
    ax.tick_params(axis="y", labelsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.03)
    cb.set_label("log₂ fold contrast", fontsize=6.8)
    cb.ax.tick_params(labelsize=6)


def main():
    required = [SURFACE, SURFACE_ZOOM, SURFACE_ZOOM_META, PROFILES, LEAVEOUT, ANALYSIS, FIT_PATH]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Wait for the analysis outputs before rendering: " + ", ".join(missing))
    surfaces = read_csv(SURFACE)
    zoom_surfaces = read_csv(SURFACE_ZOOM)
    zoom_meta = json.loads(SURFACE_ZOOM_META.read_text())
    profiles = read_csv(PROFILES)
    leaveout = read_csv(LEAVEOUT)
    meta = json.loads(ANALYSIS.read_text())
    if meta.get("smoke_run", True):
        raise ValueError("Refusing to render a smoke run as the review figure")
    if not meta.get("surface", {}).get("all_drift_optimizers_successful", False):
        raise ValueError("Analysis metadata does not confirm convergence of every surface drift refit")
    zoom_tau = {number(r, "tau_ratio") for r in zoom_surfaces}
    zoom_sigma = {number(r, "sigma_ratio") for r in zoom_surfaces}
    zoom_coordinates = {(number(r, "tau_ratio"), number(r, "sigma_ratio")) for r in zoom_surfaces}
    expected_zoom_points = len(zoom_tau) * len(zoom_sigma)
    if (len(zoom_surfaces) != expected_zoom_points or len(zoom_coordinates) != expected_zoom_points
            or any(not is_true(r.get("success", False)) or r.get("status") != "valid"
                   for r in zoom_surfaces)):
        raise ValueError("Focused surface must be a complete rectangular grid of converged points")
    grid_meta = zoom_meta.get("grid", {})
    if (int(grid_meta.get("tau_points", len(zoom_tau))) != len(zoom_tau)
            or int(grid_meta.get("sigma_points", len(zoom_sigma))) != len(zoom_sigma)
            or int(grid_meta.get("evaluated_points", len(zoom_surfaces))) != len(zoom_surfaces)
            or ("shape" in grid_meta and list(grid_meta["shape"]) != [len(zoom_tau), len(zoom_sigma)])):
        raise ValueError("Focused surface row count or axes do not match its metadata")
    if not all(is_true(r.get("selected_start", False)) and is_true(r.get("success", False))
               for r in profiles):
        profiles = [r for r in profiles if is_true(r.get("selected_start", False))
                    and is_true(r.get("success", False))]
    required_profiles = {"mu_reference", "mu_wgd", "drift_contrast"}
    if required_profiles != {r.get("profile_parameter") for r in profiles}:
        raise ValueError("Each drift profile must have selected converged points")
    needed_analyses = set(ANALYSES)
    if not needed_analyses.issubset({r.get("analysis") for r in leaveout}):
        raise ValueError("Leave-Hosny-out CSV must include all three fit/scoring states")

    # Independently reconstruct per-study and pooled contrasts from fit_05.
    fit_record = json.loads(FIT_PATH.read_text())
    fitted = np.exp(np.asarray(fit_record["z"], dtype=float)[:20].reshape(5, 2, 2))
    expected_tau = fitted[:, 1, 0] / fitted[:, 0, 0]
    expected_sigma = fitted[:, 1, 1] / fitted[:, 0, 1]
    expected_population_tau = float(np.exp(np.log(expected_tau).mean()))
    expected_population_sigma = float(np.exp(np.log(expected_sigma).mean()))
    expected_population_tau4 = float(np.exp(np.log(expected_tau[:4]).mean()))
    expected_population_sigma4 = float(np.exp(np.log(expected_sigma[:4]).mean()))
    full_fit_study_values = {study: (float(expected_tau[index]), float(expected_sigma[index]))
                             for index, study in enumerate(STUDIES[:-1])}
    for index, study in enumerate(STUDIES[:-1]):
        row = next((r for r in leaveout if r.get("analysis") == "full5_remaining4_only"
                    and r.get("study") == study), None)
        if (row is None or not np.isclose(number(row, "tau_fold"), expected_tau[index], rtol=1e-8)
                or not np.isclose(number(row, "sigma_fold"), expected_sigma[index], rtol=1e-8)):
            raise ValueError(f"Full-fit fold in leaveout table does not match fit_05: {study}")
    pooled = next(r for r in leaveout if r.get("analysis") == "full5_allstudies"
                  and r.get("study") == "Pooled")
    if (not np.isclose(number(pooled, "population_tau_fold"), expected_population_tau, rtol=1e-8)
            or not np.isclose(number(pooled, "population_sigma_fold"), expected_population_sigma, rtol=1e-8)):
        raise ValueError("Pooled full-fit fold in leaveout table does not match fit_05 log effects")
    pooled4 = next(r for r in leaveout if r.get("analysis") == "full5_remaining4_only"
                   and r.get("study") == "Pooled")
    if (not np.isclose(number(pooled4, "population_tau_fold"), expected_population_tau4, rtol=1e-8)
            or not np.isclose(number(pooled4, "population_sigma_fold"), expected_population_sigma4, rtol=1e-8)):
        raise ValueError("Four-study pooled fold from the full fit does not match fit_05")
    grid = next((r for r in surfaces
                 if np.isclose(number(r, "tau_ratio"), expected_population_tau, rtol=1e-10)
                 and np.isclose(number(r, "sigma_ratio"), expected_population_sigma, rtol=1e-10)), None)
    if grid is None or not (abs(number(grid, "fixed_drift_delta_data")) < 1e-7
                            and abs(number(grid, "fixed_drift_delta_total")) < 1e-7):
        raise ValueError("Relative surface does not contain the selected fit as a zero-delta reference")
    zoom_grid = next((r for r in zoom_surfaces
                      if np.isclose(number(r, "tau_ratio"), expected_population_tau, rtol=1e-10)
                      and np.isclose(number(r, "sigma_ratio"), expected_population_sigma, rtol=1e-10)), None)
    if zoom_grid is None or not (abs(number(zoom_grid, "fixed_drift_delta_data")) < 1e-7
                                 and abs(number(zoom_grid, "fixed_drift_delta_total")) < 1e-7):
        raise ValueError("Focused surface does not contain the selected fit at zero fixed-drift delta")

    fig = plt.figure(figsize=(15.6, 10.2), constrained_layout=False)
    gs = fig.add_gridspec(2, 3, left=.065, right=.985, top=.835, bottom=.075,
                          wspace=.37, hspace=.46, width_ratios=[1.05, 1.05, 1.0])
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[0, 2]); axD = fig.add_subplot(gs[1, 0])
    axE = fig.add_subplot(gs[1, 1]); axF = fig.add_subplot(gs[1, 2])
    fit_tau_ratio = number(meta["surface"], "fit_tau_ratio")
    fit_sigma_ratio = number(meta["surface"], "fit_sigma_ratio")
    surface_panel(axA, zoom_surfaces, surfaces, "data", "A  Data misfit (WGD− τ/σ fixed)", "Contour levels: Δ 1, 4, 9",
                  fit_tau_ratio, fit_sigma_ratio)
    surface_panel(axB, zoom_surfaces, surfaces, "total", "B  Total loss (WGD− τ/σ fixed)", "Contour levels: Δ 1, 4, 9",
                  fit_tau_ratio, fit_sigma_ratio)
    drift_profile(axC, profiles)
    hosny_tau_fold = float(fitted[4, 1, 0] / fitted[4, 0, 0])
    hosny_sigma_fold = float(fitted[4, 1, 1] / fitted[4, 0, 1])
    tradeoff(axD, profiles, hosny_tau_fold, hosny_sigma_fold)
    leaveout_heatmap(axE, leaveout, "tau_fold", "E", "Persistence after leaving Hosny out",
                     full_fit_study_values)
    leaveout_heatmap(axF, leaveout, "sigma_fold", "F", "Noise after leaving Hosny out",
                     full_fit_study_values)

    method_handles = [Line2D([0], [0], color="#8c9aa1", ls="-", lw=2.4, label="Hosny drifts fixed"),
                      Line2D([0], [0], color="#42545e", ls="--", lw=1.2, label="Both Hosny drifts refit")]
    level_handles = [Line2D([0], [0], color=LEVEL_COLORS[v], lw=1.5, label=f"Δ {v}") for v in (1, 4, 9)]
    fig.legend(handles=method_handles + level_handles, frameon=False, ncol=5,
               loc="upper center", bbox_to_anchor=(.52, .886), fontsize=7.5,
               handlelength=2.0, columnspacing=1.35)
    fig.suptitle("How much does the Hosny drift assumption change the fitted contrast?",
                 x=.065, y=.972, ha="left", color=NAVY, weight="bold", fontsize=15)
    fig.text(.065, .913,
             "Primary Gaussian OU hierarchy · profile panels reoptimize remaining parameters · Δ values are descriptive fit sensitivity",
             ha="left", color="#4d626e", fontsize=8.3)
    fig.text(.065, .018,
             "Objective is arbitrarily weighted and has no calibrated statistical or biological units; it is not a likelihood. "
             "Profile ranges are unions of sampled 1D fits, not joint uncertainty regions. Hosny contributes medians to fitting.",
             ha="left", color="#526570", fontsize=7.1)

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    summary = meta.get("summary_results", {})
    profile_ranges = summary.get("sampled_profile_parameter_ranges", {})
    range1 = profile_ranges.get("1", {})
    range9 = profile_ranges.get("9", {})
    pooled_comparisons = summary.get("pooled_fold_comparisons", {})
    full5 = pooled_comparisons.get("full5_allstudies", {})
    comp4 = pooled_comparisons.get("full5_remaining4_only", {})
    refit4 = pooled_comparisons.get("leaveHosnyOut_refit", {})
    result_text = (
        f" Sampled Hosny τ/σ folds range from {range1['hosny_tau_fold'][0]:.1f}–{range1['hosny_tau_fold'][1]:.1f}/"
        f"{range1['hosny_sigma_fold'][0]:.3f}–{range1['hosny_sigma_fold'][1]:.3f} at Δ≤1 to "
        f"{range9['hosny_tau_fold'][0]:.1f}–{range9['hosny_tau_fold'][1]:.1f}/"
        f"{range9['hosny_sigma_fold'][0]:.3f}–{range9['hosny_sigma_fold'][1]:.3f} at Δ≤9. "
        f"Pooled τ/σ folds are {full5['population_tau_fold']:.2f}/{full5['population_sigma_fold']:.3f} "
        f"(five studies), {comp4['population_tau_fold']:.2f}/{comp4['population_sigma_fold']:.3f} "
        f"(same fit, four-study composition), and {refit4['population_tau_fold']:.2f}/{refit4['population_sigma_fold']:.3f} "
        "(four-study refit)."
    )
    caption = (
        "Standalone Hosny drift sensitivity analysis for the primary Gaussian OU hierarchy. "
        "Panels A–B show the local relative persistence/noise surface with WGD− persistence/noise fixed; Hosny signed drifts are either fixed at the selected fit or both reoptimized at every point. Broad-grid context appears in each inset. "
        "Lines show changes in scaled data misfit (A) or regularized total loss (B); contour color identifies Δ 1, 4, and 9, and solid pale versus dashed saturated lines identify fixed versus refitted drifts. The two drift treatments nearly overlap. "
        "Panel C profiles each signed drift and their contrast while reoptimizing the remaining parameters. Panel D shows Hosny-specific persistence/noise folds from selected sampled one-dimensional profiles within total Δ 9; all points use the same circle and color represents increase in total fit objective. A hollow circle marks the selected fit. "
        "Panels E–F compare study-specific folds and their pooled fold from the five-study fit, the same full-fit parameters summarized across the remaining four studies, and a refit excluding Hosny. "
        "The objective uses heuristic residual scales and regularization, is not calibrated to meaningful statistical or biological units, and is not a likelihood. Δ values are descriptive sensitivity guides, not confidence intervals; unions of sampled one-dimensional profiles do not form a joint region."
        " Hosny path predictions use stable Rice calculations; a normal asymptotic fallback applies only in the extreme high drift-to-step-noise tail, as recorded in analysis provenance."
        + result_text
    )
    (OUT / "CAPTION.md").write_text(caption + "\n")
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_script": str(Path(__file__).relative_to(ROOT)),
        "figure": str(PNG.relative_to(ROOT)),
        "inputs_sha256": {str(p.relative_to(ROOT)): digest(p) for p in required},
        "analysis_selected_fit": meta.get("selected_fit"),
        "panels": ["focused fixed/refitted drift loss surfaces with broad-grid insets", "drift profiles", "profile tradeoffs", "leave-Hosny-out fold comparison"],
        "interpretation": "Heuristic regularized fit sensitivity; not a likelihood, confidence interval, probability, or joint parameter envelope.",
        "numerical_note": "Hosny path summaries use stable Rice calculations; the normal asymptotic fallback is limited to the extreme high drift-to-step-noise tail and is documented in the analysis provenance.",
        "validation": {"png_only": True, "all_input_files_present": True,
                       "surface_rows": len(surfaces), "zoom_surface_rows": len(zoom_surfaces),
                       "profile_rows": len(profiles),
                       "leaveout_rows": len(leaveout), "full_fit_folds_recomputed_from_fit_05": True,
                       "selected_fit_is_zero_delta_grid_reference": True,
                       "selected_fit_is_zero_delta_zoom_reference": True,
                       "zoom_grid_metadata": zoom_meta,
                       "finite_surface_rows": sum(
                           np.isfinite(number(r, "fixed_drift_delta_total")) and
                           np.isfinite(number(r, "optimized_drift_delta_total")) for r in surfaces)},
    }
    provenance["validation"]["finite_surface_rows"] = int(provenance["validation"]["finite_surface_rows"])
    provenance["hosny_fit_fold"] = {"tau_fold": hosny_tau_fold, "sigma_fold": hosny_sigma_fold}
    provenance["summary_results"] = meta.get("summary_results", {})
    (OUT / "figure_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({"figure": str(PNG), "panels": 6,
                      "surface_rows": len(surfaces), "profile_rows": len(profiles),
                      "leaveout_rows": len(leaveout)}, indent=2))


if __name__ == "__main__":
    main()
