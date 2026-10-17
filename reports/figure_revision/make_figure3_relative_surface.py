"""Generate the relative WGD+ persistence/noise loss surface for Figure 3."""
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
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "figures" / "figure3_review"
FIT_DIR = ROOT / "wgd_langevin/gaussian_hierarchy"
RESULTS = FIT_DIR / "results"
sys.path.insert(0, str(FIT_DIR))
spec = importlib.util.spec_from_file_location("gaussian_hierarchy_fit", FIT_DIR / "fit.py")
fit = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = fit
spec.loader.exec_module(fit)

NAVY = "#19394d"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.titleweight": "semibold"})


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    candidates = [json.loads(p.read_text()) for p in sorted(RESULTS.glob("fit_*.json"))]
    candidates = [r for r in candidates if r["scenario"] == "primary" and r["success"]]
    best = min(candidates, key=lambda r: r["total_score"])
    fit_path = RESULTS / "fit_05.json"
    if best.get("data_score") != 11.473550794223446 or best.get("total_score") != 33.06008722062528:
        raise AssertionError("Selected optimum differs from the documented fit_05 scores")
    z0 = np.asarray(best["z"], dtype=float)
    pars0 = z0[:20].reshape(5, 2, 2)
    group_geo = np.exp(pars0.mean(axis=0))
    tau_minus, sigma_minus = group_geo[0]
    tau_plus, sigma_plus = group_geo[1]
    tau_fit_fold, sigma_fit_fold = tau_plus / tau_minus, sigma_plus / sigma_minus

    # A common shift of all WGD+ study log parameters preserves fitted study
    # deviations. Intersect study-level bounds before constructing the grid.
    plus = pars0[:, 1, :]
    tau_model_bounds = (float(np.log(1e-6) - plus[:, 0].min()),
                        float(np.log(12. - 1e-8) - plus[:, 0].max()))
    sigma_model_bounds = (float(np.log(1e-4) - plus[:, 1].min()),
                          float(np.log(1e8) - plus[:, 1].max()))
    # Focus on the interval needed to show both requested landmarks, clipping
    # only to actual study-level model parameter bounds.
    tau_shift_bounds = (max(float(np.log(.5 / tau_fit_fold)), tau_model_bounds[0]),
                        min(float(np.log(30. / tau_fit_fold)), tau_model_bounds[1]))
    sigma_shift_bounds = (max(float(np.log(.05 / sigma_fit_fold)), sigma_model_bounds[0]),
                          min(float(np.log(2. / sigma_fit_fold)), sigma_model_bounds[1]))
    tau_folds = np.unique(np.r_[np.geomspace(tau_fit_fold * np.exp(tau_shift_bounds[0]),
                                               tau_fit_fold * np.exp(tau_shift_bounds[1]), 61),
                                  1., tau_fit_fold])
    sigma_folds = np.unique(np.r_[np.geomspace(sigma_fit_fold * np.exp(sigma_shift_bounds[0]),
                                                 sigma_fit_fold * np.exp(sigma_shift_bounds[1]), 61),
                                    1., sigma_fit_fold])

    # The optimizer uses Hosny medians only. Keep its exact residual mask and
    # regularization terms at every point on this conditional surface.
    rows = []
    data_delta = np.full((len(tau_folds), len(sigma_folds)), np.nan)
    total_delta = np.full_like(data_delta, np.nan)
    for ix, tau_fold in enumerate(tau_folds):
        dtau = float(np.log(tau_fold / tau_fit_fold))
        for iy, sigma_fold in enumerate(sigma_folds):
            dsigma = float(np.log(sigma_fold / sigma_fit_fold))
            valid = (tau_shift_bounds[0] <= dtau <= tau_shift_bounds[1]
                     and sigma_shift_bounds[0] <= dsigma <= sigma_shift_bounds[1])
            row = {"tau_WGDplus_over_WGDminus": tau_fold,
                   "sigma_WGDplus_over_WGDminus": sigma_fold,
                   "delta_data_misfit": "", "delta_regularized_total_loss": "",
                   "data_misfit": "", "regularized_total_loss": "",
                   "tau_log_shift_from_fit": dtau, "sigma_log_shift_from_fit": dsigma,
                   "status": "valid" if valid else "masked_parameter_bound"}
            if valid:
                z = z0.copy()
                modified = pars0.copy()
                modified[:, 1, 0] += dtau
                modified[:, 1, 1] += dsigma
                z[:20] = modified.ravel()
                residual, residual_rows = fit.data_residuals(fit.analytic(z))
                mask = np.asarray([r["study"] != "Hosny_baseline" or r.get("quantile") == 1
                                   for r in residual_rows], dtype=bool)
                data = float(residual[mask] @ residual[mask])
                pooling = fit.pooling_residuals(z, fit.SCENARIOS["primary"])
                prior = fit.persistence_residuals(z, "primary")
                total = float(data + pooling @ pooling + prior @ prior)
                data_delta[ix, iy] = data - best["data_score"]
                total_delta[ix, iy] = total - best["total_score"]
                row.update(delta_data_misfit=data_delta[ix, iy],
                           delta_regularized_total_loss=total_delta[ix, iy],
                           data_misfit=data, regularized_total_loss=total)
            rows.append(row)

    # Directly verify the stored optimum and both requested reference points.
    ref_resid, ref_rows = fit.data_residuals(fit.analytic(z0))
    ref_mask = np.asarray([r["study"] != "Hosny_baseline" or r.get("quantile") == 1
                           for r in ref_rows], dtype=bool)
    checked_data = float(ref_resid[ref_mask] @ ref_resid[ref_mask])
    pool = fit.pooling_residuals(z0, fit.SCENARIOS["primary"])
    prior = fit.persistence_residuals(z0, "primary")
    checked_total = float(checked_data + pool @ pool + prior @ prior)
    if not np.isclose(checked_data, 11.473550794223446, rtol=0, atol=1e-10):
        raise AssertionError(f"Stored data score mismatch: {checked_data}")
    if not np.isclose(checked_total, 33.06008722062528, rtol=0, atol=1e-9):
        raise AssertionError(f"Stored total loss mismatch: {checked_total}")
    if not (np.isfinite(data_delta[np.isfinite(data_delta)]).all()
            and np.isfinite(total_delta[np.isfinite(total_delta)]).all()):
        raise AssertionError("Non-finite evaluated loss on grid")
    ix = int(np.flatnonzero(tau_folds == tau_fit_fold)[0])
    iy = int(np.flatnonzero(sigma_folds == sigma_fit_fold)[0])
    if not (np.isclose(data_delta[ix, iy], 0, atol=1e-9)
            and np.isclose(total_delta[ix, iy], 0, atol=1e-9)):
        raise AssertionError("Fitted WGD+ point is not present as the zero-delta reference")
    ex = int(np.flatnonzero(tau_folds == 1.0)[0])
    ey = int(np.flatnonzero(sigma_folds == 1.0)[0])
    if not (np.isfinite(data_delta[ex, ey]) and np.isfinite(total_delta[ex, ey])):
        raise AssertionError("Equal-group-means point is outside true parameter bounds")
    for row in rows:
        if row["status"] == "valid":
            assert tau_shift_bounds[0] <= row["tau_log_shift_from_fit"] <= tau_shift_bounds[1]
            assert sigma_shift_bounds[0] <= row["sigma_log_shift_from_fit"] <= sigma_shift_bounds[1]

    valid_data = data_delta[np.isfinite(data_delta)]
    valid_total = total_delta[np.isfinite(total_delta)]
    scale = .25
    transform = lambda v: np.arcsinh(v / scale)
    panels = [(data_delta, "Change in data misfit"),
              (total_delta, "Change in regularized total loss")]
    fig, axs = plt.subplots(1, 2, figsize=(12.0, 5.0), constrained_layout=False)
    fig.subplots_adjust(left=.08, right=.90, bottom=.15, top=.78, wspace=.30)
    for ax, (raw, title), values in zip(axs, panels, [valid_data, valid_total]):
        lim = float(np.max(np.abs(transform(values))))
        image = ax.pcolormesh(tau_folds, sigma_folds, transform(raw).T,
                              shading="auto", cmap="coolwarm", norm=Normalize(-lim, lim))
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.scatter([1], [1], marker="o", s=62, facecolor="white", edgecolor="#222", linewidth=1.1,
                   zorder=5, label="Equal group means")
        ax.scatter([tau_fit_fold], [sigma_fit_fold], marker="*", s=170, facecolor="#111",
                   edgecolor="white", linewidth=1.0, zorder=6,
                   label="Selected fitted WGD+ means")
        ax.set_xlabel("WGD+ / WGD− geometric-mean persistence")
        ax.set_ylabel("WGD+ / WGD− geometric-mean noise")
        ax.set_title(title, loc="left", color=NAVY, fontsize=10, pad=8)
        ax.grid(which="major", color="white", alpha=.4, linewidth=.6)
        cb = fig.colorbar(image, ax=ax, fraction=.046, pad=.04)
        ticks = np.linspace(-lim, lim, 5)
        cb.set_ticks(ticks)
        cb.set_ticklabels([f"{scale*np.sinh(t):.3g}" for t in ticks])
        cb.ax.tick_params(labelsize=7)
        ax.legend(loc="best", frameon=True, fontsize=7)
    fig.suptitle("Relative persistence–noise loss surface with WGD− fixed",
                 x=.08, y=.94, ha="left", color=NAVY, weight="bold", fontsize=14)
    fig.text(.08, .865,
             "Each point shifts WGD+ study log parameters together; fitted study deviations, WGD−, and Hosny drifts stay fixed. "
             "The circle marks equal geometric group means; the star marks fit_05. "
             "Color uses signed asinh(Δ / 0.25); colorbar ticks show raw score changes.", fontsize=8.2, color="#384c59")
    fig.savefig(OUT / "figure3_relative_persistence_noise_surface.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with (OUT / "figure3_relative_persistence_noise_surface.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_fit": fit_path.name, "selected_fit_sha256": sha(fit_path),
        "selected_data_misfit": best["data_score"], "selected_regularized_total_loss": best["total_score"],
        "group_geometric_means": {"WGD_minus": {"tau_h": float(tau_minus), "sigma_um_h_minus_3_over_2": float(sigma_minus)},
                                  "WGD_plus": {"tau_h": float(tau_plus), "sigma_um_h_minus_3_over_2": float(sigma_plus)}},
        "relative_markers": {"equal_group_means": [1.0, 1.0],
                              "selected_fitted_WGDplus_fold": [float(tau_fit_fold), float(sigma_fit_fold)]},
        "validation": {"recomputed_data_misfit": checked_data,
                       "recomputed_regularized_total_loss": checked_total,
                       "stored_scores_match": True,
                       "all_evaluated_grid_values_finite": True,
                       "evaluated_grid_points": int(np.isfinite(data_delta).sum()),
                       "masked_grid_points": int(sum(r["status"] != "valid" for r in rows)),
                       "equal_group_means_in_grid": bool(1.0 in tau_folds and 1.0 in sigma_folds),
                       "equal_group_means_evaluated_finite": bool(np.isfinite(data_delta[ex, ey]) and np.isfinite(total_delta[ex, ey])),
                       "selected_fit_in_grid": True,
                       "selected_fit_evaluated_finite": bool(np.isfinite(data_delta[ix, iy]) and np.isfinite(total_delta[ix, iy]))},
        "definition": "WGD− parameters are fixed. WGD+ study-specific log persistence/noise parameters are shifted together, preserving fitted deviations; other parameters and Hosny drifts are frozen. Data score uses the exact fitting residual mask, including Hosny medians only. Total loss adds the existing primary pooling and persistence regularizers.",
        "relative_loss": "Both scores are reported as changes from the selected fit. Equal group means at (1,1) are descriptive group-mean equality, not study-by-study parameter equality.",
        "display": "Heatmap color is signed asinh(delta / 0.25); colorbar tick labels report raw score changes.",
        "input_sha256": {str(FIT_DIR / "fit.py"): sha(FIT_DIR / "fit.py"),
                         str(FIT_DIR / "model.py"): sha(FIT_DIR / "model.py")},
        "grid": {"shape": [int(len(tau_folds)), int(len(sigma_folds))],
                 "tau_ratio_bounds": [float(tau_folds[0]), float(tau_folds[-1])],
                 "sigma_ratio_bounds": [float(sigma_folds[0]), float(sigma_folds[-1])]}
    }
    (OUT / "figure3_relative_persistence_noise_surface.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"fit": fit_path.name, "relative_markers": metadata["relative_markers"],
                      "evaluated_points": metadata["validation"]["evaluated_grid_points"],
                      "output": str(OUT / "figure3_relative_persistence_noise_surface.png")}, indent=2))


if __name__ == "__main__":
    main()
