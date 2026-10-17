"""Refine the local Hosny drift surface around the selected fit."""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CAMPAIGN = ROOT / "wgd_langevin/refined_ou_20261003"
OUT = CAMPAIGN / "results/hosny_sensitivity"
sys.path.insert(0, str(HERE))
import analyze as analysis  # imports the locally stabilized Rice calculations


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def contour_summary(rows: list[dict], labels: dict[str, str]) -> dict:
    summary = {}
    for label, field in labels.items():
        valid = [r for r in rows if r.get("status") == "valid" and r.get("success") in (True, "True")]
        values = np.asarray([float(r[field]) for r in valid], dtype=float)
        x = np.asarray([float(r["tau_ratio"]) for r in valid])
        y = np.asarray([float(r["sigma_ratio"]) for r in valid])
        summary[label] = {}
        for threshold in (1., 4., 9.):
            inside = values <= threshold
            summary[label][str(int(threshold))] = {
                "sampled_point_count": int(inside.sum()),
                "tau_ratio_bounds": [float(x[inside].min()), float(x[inside].max())] if inside.any() else None,
                "sigma_ratio_bounds": [float(y[inside].min()), float(y[inside].max())] if inside.any() else None}
    return summary


def main() -> None:
    fit_path = analysis.RESULTS / analysis.FIT_NAME
    base_obj = json.loads(fit_path.read_text())
    base = np.asarray(base_obj["z"], dtype=float)
    base_scores = analysis.scores(base)
    if not (np.isclose(base_scores["data_score"], analysis.FIT_DATA_SCORE, rtol=0, atol=1e-10)
            and np.isclose(base_scores["total_score"], analysis.FIT_TOTAL_SCORE, rtol=0, atol=1e-9)):
        raise AssertionError("The base fit scores no longer reproduce selected primary fit")

    pars = base[:20].reshape(5, 2, 2)
    geo = np.exp(pars.mean(axis=0))
    fit_tau = float(geo[1, 0] / geo[0, 0])
    fit_sigma = float(geo[1, 1] / geo[0, 1])
    plus = pars[:, 1, :]
    tau_bounds = (float(np.log(1e-6) - plus[:, 0].min()),
                  float(np.log(12. - 1e-8) - plus[:, 0].max()))
    sigma_bounds = (float(np.log(1e-4) - plus[:, 1].min()),
                    float(np.log(1e8) - plus[:, 1].max()))
    # Reserve one coordinate on each axis for the exact selected-fit value.
    tau_axis = np.unique(np.r_[np.geomspace(fit_tau / 1.45, fit_tau * 1.45, 60), fit_tau])
    sigma_axis = np.unique(np.r_[np.geomspace(fit_sigma / 1.45, fit_sigma * 1.45, 60), fit_sigma])
    rows = []
    for tau_ratio in tau_axis:
        dt = float(np.log(tau_ratio / fit_tau))
        if not tau_bounds[0] <= dt <= tau_bounds[1]:
            raise AssertionError(f"Local tau coordinate exceeds parameter bounds: {tau_ratio}")
        for sigma_ratio in sigma_axis:
            ds = float(np.log(sigma_ratio / fit_sigma))
            if not sigma_bounds[0] <= ds <= sigma_bounds[1]:
                raise AssertionError(f"Local sigma coordinate exceeds parameter bounds: {sigma_ratio}")
            z_fixed = base.copy()
            shifted = pars.copy()
            shifted[:, 1, 0] += dt
            shifted[:, 1, 1] += ds
            z_fixed[:20] = shifted.ravel()
            fixed = analysis.scores(z_fixed)
            z_opt, sol = analysis.run_drift_pair(z_fixed)
            optimized = analysis.scores(z_opt)
            in_bounds = bool(np.all(z_opt >= analysis.fit.LOW - 1e-10)
                             and np.all(z_opt <= analysis.fit.HIGH + 1e-10))
            finite = bool(np.isfinite(list(fixed.values())).all()
                          and np.isfinite(list(optimized.values())).all())
            success = bool(sol.success and in_bounds and finite)
            rows.append({
                "tau_ratio": float(tau_ratio), "sigma_ratio": float(sigma_ratio),
                "fixed_drift_reference": float(z_fixed[20]),
                "fixed_drift_wgd": float(z_fixed[21]),
                "fixed_drift_delta_data": fixed["data_score"] - base_scores["data_score"],
                "fixed_drift_delta_total": fixed["total_score"] - base_scores["total_score"],
                "optimized_drift_reference": float(z_opt[20]),
                "optimized_drift_wgd": float(z_opt[21]),
                "optimized_drift_delta": float(z_opt[21] - z_opt[20]),
                "optimized_drift_delta_data": optimized["data_score"] - base_scores["data_score"],
                "optimized_drift_delta_total": optimized["total_score"] - base_scores["total_score"],
                "z": json.dumps(z_opt.tolist()),
                "success": success, "nfev": int(sol.nfev),
                "status": "valid" if success else "optimizer_failed_or_invalid"})

    if not all(r["success"] for r in rows):
        failed = sum(not r["success"] for r in rows)
        raise RuntimeError(f"{failed} local surface points did not converge inside bounds")
    values = np.asarray([[float(r[k]) for k in (
        "fixed_drift_delta_data", "fixed_drift_delta_total",
        "optimized_drift_delta_data", "optimized_drift_delta_total")] for r in rows])
    if not np.isfinite(values).all():
        raise AssertionError("Non-finite loss in local surface")
    max_nonworsening = max(r["optimized_drift_delta_total"] - r["fixed_drift_delta_total"] for r in rows)
    if max_nonworsening > 1e-7:
        raise AssertionError("Drift reoptimization worsened the regularized total objective")

    OUT.mkdir(parents=True, exist_ok=True)
    zoom_csv = OUT / "hosny_surface_zoom.csv"
    with zoom_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    full_rows = read_rows(OUT / "hosny_surface.csv")
    contour_labels = {
        "fixed_data": "fixed_drift_delta_data",
        "fixed_total": "fixed_drift_delta_total",
        "optimized_data": "optimized_drift_delta_data",
        "optimized_total": "optimized_drift_delta_total"}
    meta = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selected_fit": analysis.FIT_NAME,
        "fit_sha256": sha(fit_path),
        "analysis_script_sha256": sha(HERE / "analyze.py"),
        "refinement_script_sha256": sha(Path(__file__).resolve()),
        "grid": {"tau_ratio_bounds": [float(tau_axis[0]), float(tau_axis[-1])],
                 "sigma_ratio_bounds": [float(sigma_axis[0]), float(sigma_axis[-1])],
                 "tau_points": int(len(tau_axis)), "sigma_points": int(len(sigma_axis)),
                 "shape": [int(len(tau_axis)), int(len(sigma_axis))],
                 "evaluated_points": len(rows),
                 "fit_ratio": [fit_tau, fit_sigma]},
        "validation": {"all_points_converged": all(r["success"] for r in rows),
                       "all_values_finite": bool(np.isfinite(values).all()),
                       "all_parameters_within_model_bounds": True,
                       "max_optimized_minus_fixed_total_delta": float(max_nonworsening),
                       "nonworsening_tolerance": 1e-7,
                       "selected_fit_scores_recomputed": base_scores},
        "sampled_contour_bounds_by_score_and_threshold": {
            "full_grid": contour_summary(full_rows, contour_labels),
            "local_refinement": contour_summary(rows, contour_labels)},
        "interpretation": "The local surface refines the fixed/refit drift Δ1/4/9 contours around selected primary fit. Contour bounds are sampled-grid extents, not continuous intervals.",
        "output_sha256": {zoom_csv.name: sha(zoom_csv)}}
    (OUT / "hosny_surface_zoom.json").write_text(json.dumps(meta, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"grid_points": len(rows), "all_converged": True,
                      "max_nonworsening": max_nonworsening,
                      "output": str(zoom_csv)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
