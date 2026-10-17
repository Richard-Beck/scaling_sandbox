"""Build conditional persistence/noise loss surfaces for revised Figure 3."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUT = HERE / "figures" / "figure3_review"
FIT_DIR = ROOT / "wgd_langevin/gaussian_hierarchy"
RESULTS = FIT_DIR / "results"
sys.path.insert(0, str(FIT_DIR))

spec = importlib.util.spec_from_file_location("gaussian_hierarchy_fit", FIT_DIR / "fit.py")
fit = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = fit
spec.loader.exec_module(fit)

NAVY = "#19394d"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "pdf.fonttype": 42, "axes.titleweight": "semibold"})


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    OUT.mkdir(exist_ok=True)
    candidates = [json.loads(p.read_text()) for p in sorted(RESULTS.glob("fit_*.json"))]
    candidates = [r for r in candidates if r["scenario"] == "primary" and r["success"]]
    if not candidates:
        raise RuntimeError("No successful primary Gaussian-hierarchy fits found")
    best = min(candidates, key=lambda r: r["total_score"])
    z0 = np.asarray(best["z"], dtype=float)
    pars0 = z0[:20].reshape(5, 2, 2)
    centers = np.exp(pars0[:, :, :])
    group_names = [(0, "WGD−", "Reference"), (1, "WGD+", "WGD")]
    ngrid = 49
    # Cover a broad local range, while clipping the persistence grid to the
    # common domain allowed for every study-specific value in the swept group.
    group_specs = {}
    for group, label, condition in group_names:
        base_logs = pars0[:, group, :].copy()
        tau_center = float(np.exp(base_logs[:, 0].mean()))
        sigma_center = float(np.exp(base_logs[:, 1].mean()))
        tau_lo = max(-2.0, float(np.log(1e-6) - base_logs[:, 0].min()))
        tau_hi = min(2.0, float(np.log(12. - 1e-8) - base_logs[:, 0].max()))
        sigma_lo = max(-2.0, float(np.log(1e-4) - base_logs[:, 1].min()))
        sigma_hi = min(2.0, float(np.log(1e8) - base_logs[:, 1].max()))
        group_specs[group] = dict(label=label, condition=condition, base_logs=base_logs,
                                  center_tau=tau_center, center_sigma=sigma_center,
                                  tau_log_bounds=(tau_lo, tau_hi),
                                  sigma_log_bounds=(sigma_lo, sigma_hi))
    surfaces = {}
    rows = []
    for group, label, condition in group_names:
        spec_g = group_specs[group]
        base_logs = spec_g["base_logs"]
        tau_center = spec_g["center_tau"]
        sigma_center = spec_g["center_sigma"]
        tau_shifts = np.unique(np.r_[np.linspace(*spec_g["tau_log_bounds"], ngrid), 0.0])
        sigma_shifts = np.unique(np.r_[np.linspace(*spec_g["sigma_log_bounds"], ngrid), 0.0])
        tau_values = tau_center * np.exp(tau_shifts)
        sigma_values = sigma_center * np.exp(sigma_shifts)
        grid_shape = (len(tau_values), len(sigma_values))
        data_delta = np.full(grid_shape, np.nan)
        total_delta = np.full(grid_shape, np.nan)
        data_center = float(best["data_score"])
        total_center = float(best["total_score"])
        for ix, tau in enumerate(tau_values):
            dtau = np.log(tau / tau_center)
            for iy, sigma in enumerate(sigma_values):
                dsigma = np.log(sigma / sigma_center)
                valid_tau = spec_g["tau_log_bounds"][0] <= dtau <= spec_g["tau_log_bounds"][1]
                valid_sigma = spec_g["sigma_log_bounds"][0] <= dsigma <= spec_g["sigma_log_bounds"][1]
                if not (valid_tau and valid_sigma):
                    rows.append({"condition": label, "tau_h": tau,
                                 "sigma_um_h_minus_3_over_2": sigma,
                                 "delta_data_score": "", "delta_regularized_objective": "",
                                 "evaluated_data_score": "", "regularized_objective": "",
                                 "tau_log_shift": dtau, "sigma_log_shift": dsigma,
                                 "status": "masked_parameter_bound"})
                    continue
                z = z0.copy()
                modified = pars0.copy()
                modified[:, group, 0] += dtau
                modified[:, group, 1] += dsigma
                z[:20] = modified.ravel()
                dr, residual_rows = fit.data_residuals(fit.analytic(z))
                data_mask = np.asarray([
                    row["study"] != "Hosny_baseline" or row.get("quantile") == 1
                    for row in residual_rows
                ], dtype=bool)
                dr = dr[data_mask]
                data_score = float(dr @ dr)
                pool = fit.pooling_residuals(z, fit.SCENARIOS["primary"])
                prior = fit.persistence_residuals(z, "primary")
                total_score = float(data_score + pool @ pool + prior @ prior)
                data_delta[ix, iy] = data_score - data_center
                total_delta[ix, iy] = total_score - total_center
                rows.append({"condition": label, "tau_h": tau, "sigma_um_h_minus_3_over_2": sigma,
                             "delta_data_score": data_delta[ix, iy],
                             "delta_regularized_objective": total_delta[ix, iy],
                             "evaluated_data_score": data_score,
                             "regularized_objective": total_score,
                             "tau_log_shift": dtau, "sigma_log_shift": dsigma,
                             "status": "valid"})
        surfaces[group] = dict(label=label, condition=condition, tau=tau_values,
                               sigma=sigma_values, data=data_delta, total=total_delta,
                               center_tau=tau_center, center_sigma=sigma_center,
                               log_tau_bounds=spec_g["tau_log_bounds"],
                               log_sigma_bounds=spec_g["sigma_log_bounds"], grid_shape=grid_shape)

    all_data = np.array([r["delta_data_score"] for r in rows if r["status"] == "valid"], dtype=float)
    all_total = np.array([r["delta_regularized_objective"] for r in rows if r["status"] == "valid"], dtype=float)
    # Verify the stored optimizer point independently, using the exact fitting
    # mask for Hosny (median rows only).
    reference_resid, reference_rows = fit.data_residuals(fit.analytic(z0))
    fit_mask = np.asarray([r["study"] != "Hosny_baseline" or r.get("quantile") == 1
                           for r in reference_rows], dtype=bool)
    recomputed_data = float(reference_resid[fit_mask] @ reference_resid[fit_mask])
    reference_pool = fit.pooling_residuals(z0, fit.SCENARIOS["primary"])
    reference_prior = fit.persistence_residuals(z0, "primary")
    recomputed_total = float(recomputed_data + reference_pool @ reference_pool + reference_prior @ reference_prior)
    if not np.isclose(recomputed_data, best["data_score"], rtol=1e-9, atol=1e-9):
        raise AssertionError(f"Data-score validation failed: {recomputed_data} vs {best['data_score']}")
    if not np.isclose(recomputed_total, best["total_score"], rtol=1e-9, atol=1e-9):
        raise AssertionError(f"Objective validation failed: {recomputed_total} vs {best['total_score']}")
    if not np.isfinite(all_data).all() or not np.isfinite(all_total).all():
        raise AssertionError("Non-finite value in the evaluated surface grid")
    for group in (0, 1):
        s = surfaces[group]
        cx = int(np.flatnonzero(s["tau"] == s["center_tau"])[0])
        cy = int(np.flatnonzero(s["sigma"] == s["center_sigma"])[0])
        if not (np.isclose(s["data"][cx, cy], 0., atol=1e-9) and
                np.isclose(s["total"][cx, cy], 0., atol=1e-9)):
            raise AssertionError(f"Fitted group-mean marker is not the zero-delta reference for group {group}")
        for ix, tau in enumerate(s["tau"]):
            for iy, sigma in enumerate(s["sigma"]):
                if np.isfinite(s["data"][ix, iy]):
                    dtau = np.log(tau / s["center_tau"])
                    dsigma = np.log(sigma / s["center_sigma"])
                    if not (s["log_tau_bounds"][0] <= dtau <= s["log_tau_bounds"][1]):
                        raise AssertionError("Evaluated persistence exceeds parameter bounds")
                    if not (s["log_sigma_bounds"][0] <= dsigma <= s["log_sigma_bounds"][1]):
                        raise AssertionError("Evaluated noise exceeds parameter bounds")
    # Signed asinh display keeps the fitted valley visible beside distant,
    # poor-fitting parameter combinations. Raw deltas remain in the CSV.
    display_scale = 0.25
    transform = lambda value: np.arcsinh(value / display_scale)
    data_display = transform(all_data)
    total_display = transform(all_total)
    data_lim = float(np.nanmax(np.abs(data_display))) or 1.
    total_lim = float(np.nanmax(np.abs(total_display))) or 1.
    fig = plt.figure(figsize=(13.1, 9.0), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, left=.06, right=.88, bottom=.07, top=.87,
                          wspace=.16, hspace=.20)
    panels = [(0, "data", "WGD−: data misfit"), (1, "data", "WGD+: data misfit"),
              (0, "total", "WGD−: regularized objective"),
              (1, "total", "WGD+: regularized objective")]
    for idx, (group, which, title) in enumerate(panels):
        ax = fig.add_subplot(gs[idx // 2, idx % 2], projection="3d")
        s = surfaces[group]
        X, Y = np.meshgrid(np.log10(s["tau"]), np.log10(s["sigma"]), indexing="ij")
        raw_Z = s[which]
        Z = transform(raw_Z)
        lim = data_lim if which == "data" else total_lim
        norm = Normalize(-lim, lim)
        ax.computed_zorder = False
        ax.plot_surface(X, Y, Z, facecolors=cm.coolwarm(norm(Z)), rstride=1,
                        cstride=1, linewidth=0, antialiased=True, shade=False,
                        alpha=.96)
        ax.scatter([np.log10(s["center_tau"])], [np.log10(s["center_sigma"])], [0], marker="*",
                   s=140, c="#111111", edgecolors="white", linewidths=1.2,
                   depthshade=False, zorder=100, label=f'{s["condition"]} group mean')
        from matplotlib.ticker import FuncFormatter
        from matplotlib.ticker import MaxNLocator, NullLocator
        ax.xaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{10**val:.2g}"))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda val, pos: f"{10**val:.2g}"))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())
        ax.set_xlim(np.log10(s["tau"][0]), np.log10(s["tau"][-1]))
        ax.set_ylim(np.log10(s["sigma"][0]), np.log10(s["sigma"][-1]))
        ax.set_xlabel("Persistence τ (h)", labelpad=4)
        ax.set_ylabel("Noise σ (µm h⁻³ᐟ²)", labelpad=6)
        ax.set_zlabel("asinh(Δ data score / 0.25)" if which == "data" else "asinh(Δ objective / 0.25)", labelpad=3)
        ax.set_title(title, loc="left", pad=3, color=NAVY, fontsize=10)
        ax.view_init(elev=27, azim=-132)
        ax.tick_params(axis="both", labelsize=7, pad=1)
        ax.tick_params(axis="z", labelsize=7, pad=1)
        ax.legend(loc="upper left", frameon=False, fontsize=7)
        ax.set_zlim(float(np.nanmin(Z)), float(np.nanmax(Z)))
        ax.set_box_aspect((1.2, 1.0, .70))
    fig.suptitle("Persistence–noise cost surfaces around the selected primary fit",
                 x=.06, y=.975, ha="left", color=NAVY, weight="bold", fontsize=14)
    fig.text(.06, .935,
             "Each panel shifts one condition’s shared log means; study deviations and the other condition are fixed. "
             "Stars mark fitted geometric means. Vertical height and color use signed asinh scaling; raw deltas are in CSV.",
             fontsize=8.5, color="#384c59")
    # Distinct shared scales for the two score definitions.
    for which, lim, x in [("data", data_lim, .925), ("total", total_lim, .965)]:
        sm = plt.cm.ScalarMappable(norm=Normalize(-lim, lim), cmap="coolwarm")
        cax = fig.add_axes([x, .25, .012, .50])
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label("asinh(Δ data / 0.25)" if which == "data" else "asinh(Δ objective / 0.25)",
                     fontsize=7, labelpad=6)
        cb.ax.tick_params(labelsize=7)
    fig.savefig(OUT / "figure3_persistence_noise_surface.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # A paired heatmap makes the narrow minima and condition-specific markers
    # easier to compare than perspective alone.
    fig, axs = plt.subplots(2, 2, figsize=(12.0, 8.3), constrained_layout=False)
    fig.subplots_adjust(left=.09, right=.90, bottom=.10, top=.88, hspace=.30, wspace=.23)
    heatpanels = [(0, "data", "WGD−: data misfit"), (1, "data", "WGD+: data misfit"),
                  (0, "total", "WGD−: regularized objective"),
                  (1, "total", "WGD+: regularized objective")]
    for idx, (group, which, panel_title) in enumerate(heatpanels):
        ax = axs.flat[idx]
        s = surfaces[group]
        Z = transform(s[which])
        lim = data_lim if which == "data" else total_lim
        image = ax.pcolormesh(s["tau"], s["sigma"], Z.T, shading="auto", cmap="coolwarm",
                              norm=Normalize(-lim, lim), rasterized=True)
        ax.scatter(s["center_tau"], s["center_sigma"], marker="*", s=150,
                   c="#111111", edgecolors="white", linewidths=1.2, zorder=5)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Persistence τ (h)"); ax.set_ylabel("Noise σ (µm h⁻³ᐟ²)")
        ax.set_title(panel_title, loc="left", color=NAVY, fontsize=10, pad=6)
        ax.grid(which="major", color="white", alpha=.35, linewidth=.6)
        cb = fig.colorbar(image, ax=ax, fraction=.045, pad=.03)
        cb.set_label("asinh(Δ score / 0.25)", fontsize=7, labelpad=5)
        cb.ax.tick_params(labelsize=7)
    fig.suptitle("Conditional persistence–noise loss around the fitted group means",
                 x=.09, y=.97, ha="left", color=NAVY, weight="bold", fontsize=14)
    fig.text(.09, .925,
             "Signed data and regularized-objective changes; star = fitted geometric mean. "
             "Color uses asinh(Δ score / 0.25); raw values are in the companion CSV.",
             fontsize=8.5, color="#384c59")
    fig.savefig(OUT / "figure3_persistence_noise_surface_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    with (OUT / "figure3_persistence_noise_surface.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_fit": min((p for p in RESULTS.glob("fit_*.json")
                             if json.loads(p.read_text())["scenario"] == "primary"
                             and json.loads(p.read_text())["success"]),
                            key=lambda p: json.loads(p.read_text())["total_score"]).name,
        "selected_fit_sha256": sha(min((p for p in RESULTS.glob("fit_*.json")
                                        if json.loads(p.read_text())["scenario"] == "primary"
                                        and json.loads(p.read_text())["success"]),
                                       key=lambda p: json.loads(p.read_text())["total_score"])),
        "selected_total_score": best["total_score"],
        "selected_data_score": best["data_score"],
        "validation": {"recomputed_data_score": recomputed_data,
                       "stored_data_score": best["data_score"],
                       "recomputed_total_score": recomputed_total,
                       "stored_total_score": best["total_score"],
                       "fit_scores_match": True,
                       "all_evaluated_grid_values_finite": bool(np.isfinite(all_data).all() and np.isfinite(all_total).all()),
                       "evaluated_grid_points": int(len(all_data)),
                       "masked_grid_points": int(sum(r["status"] != "valid" for r in rows))},
        "definition": "For each condition, sweep its study-specific log(tau) and log(sigma) values by common offsets. Keep within-condition study deviations, the other condition, and Hosny signed drift fixed. Recompute the existing analytic data residual score and the full primary regularized objective.",
        "relative_value": "Each panel is score minus the selected-fit score; negative data deltas are retained. This is conditional sensitivity, not a profile fit or a likelihood surface.",
        "display_transform": "Signed asinh(delta / 0.25) for height and color; all raw score deltas are retained in the CSV.",
        "group_means": {surfaces[g]["label"]: {"tau_h": surfaces[g]["center_tau"],
                            "sigma_um_h_minus_3_over_2": surfaces[g]["center_sigma"]}
                        for g in (0, 1)},
        "limits": {surfaces[g]["label"]: {"tau_log_shift": surfaces[g]["log_tau_bounds"],
                    "sigma_log_shift": surfaces[g]["log_sigma_bounds"]} for g in (0, 1)},
        "physical_axes_by_condition": {surfaces[g]["label"]: {
                    "tau_h": [float(surfaces[g]["tau"][0]), float(surfaces[g]["tau"][-1])],
                    "sigma_um_h_minus_3_over_2": [float(surfaces[g]["sigma"][0]), float(surfaces[g]["sigma"][-1])]}
                    for g in (0, 1)},
        "inputs": {str(FIT_DIR / "fit.py"): sha(FIT_DIR / "fit.py"),
                   str(FIT_DIR / "model.py"): sha(FIT_DIR / "model.py")},
        "grid_shape_by_condition": {surfaces[g]["label"]: list(surfaces[g]["grid_shape"])
                                     for g in (0, 1)},
        "panels": ["WGD− and WGD+ data-score deltas", "WGD− and WGD+ regularized-objective deltas"]
    }
    (OUT / "figure3_persistence_noise_surface.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"fit": metadata["selected_fit"], "data_delta_range": [float(all_data.min()), float(all_data.max())],
                      "objective_delta_range": [float(all_total.min()), float(all_total.max())],
                      "rows": len(rows), "output": str(OUT / "figure3_persistence_noise_surface.png")}, indent=2))


if __name__ == "__main__":
    main()
