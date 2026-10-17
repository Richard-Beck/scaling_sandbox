"""Hosny drift and leave-one-study-out sensitivity analyses for Figure 3."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import i0e, i1e
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[4]
CAMPAIGN = ROOT / "wgd_langevin/refined_ou_20261003"
FIT_DIR = CAMPAIGN / "scripts"
RESULTS = CAMPAIGN / "results"
OUT = CAMPAIGN / "results/hosny_sensitivity"
sys.path.insert(0, str(FIT_DIR))
spec = importlib.util.spec_from_file_location("gaussian_hierarchy_fit", FIT_DIR / "fit.py")
fit = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = fit
spec.loader.exec_module(fit)
from fast_residuals import compile_residuals
fit.data_residuals = compile_residuals(fit)

_rice_mean_original = fit.rice.mean
_rice_ppf_original = fit.rice.ppf


def rice_mean_stable(b, loc=0, scale=1):
    """Rice mean using exponentially scaled Bessel functions at finite b."""
    b_arr = np.asarray(b, dtype=float)
    out = np.empty_like(b_arr)
    large = b_arr > 1e4
    # For very large noncentrality, Rice converges to a normal radial
    # displacement; its mean is nu + sigma^2/(2 nu) to first order.
    out[large] = scale * (b_arr[large] + .5 / b_arr[large]) + loc
    x = b_arr[~large] ** 2 / 4.
    out[~large] = (scale * np.sqrt(np.pi / 2.) *
                   ((1. + 2. * x) * i0e(x) + 2. * x * i1e(x)) + loc)
    return float(out) if out.ndim == 0 else out


def rice_ppf_stable(q, b, loc=0, scale=1):
    """Use a large-noncentrality normal approximation only where needed."""
    b_arr = np.asarray(b, dtype=float)
    if np.all(b_arr > 1e4):
        return scale * (b_arr + norm.ppf(q)) + loc
    if np.any(b_arr > 1e4):
        exact = _rice_ppf_original(q, b, loc=loc, scale=scale)
        approx = scale * (b_arr + norm.ppf(q)) + loc
        return np.where(b_arr > 1e4, approx, exact)
    return _rice_ppf_original(q, b, loc=loc, scale=scale)


fit.rice.mean = rice_mean_stable
fit.rice.ppf = rice_ppf_stable

_selected = min(((p, json.loads(p.read_text())) for p in RESULTS.glob("fit_*.json") if json.loads(p.read_text()).get("success") and json.loads(p.read_text()).get("scenario") == "primary"), key=lambda x: x[1]["total_score"])
FIT_NAME = _selected[0].name
FIT_DATA_SCORE = _selected[1]["data_score"]
FIT_TOTAL_SCORE = _selected[1]["total_score"]
FIT_MASK = lambda row: row["study"] != "Hosny_baseline" or row.get("quantile") == 1


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scores(z: np.ndarray) -> dict:
    data_r, rows = fit.data_residuals(fit.analytic(z))
    mask = np.asarray([FIT_MASK(row) for row in rows], dtype=bool)
    pool = fit.pooling_residuals(z, fit.SCENARIOS["primary"])
    prior = fit.persistence_residuals(z, "primary")
    data_score = float(data_r[mask] @ data_r[mask])
    pool_score = float(pool @ pool)
    prior_score = float(prior @ prior)
    return {"data_score": data_score, "pooling_penalty": pool_score,
            "persistence_penalty": prior_score,
            "total_score": data_score + pool_score + prior_score}


def hosny_residual(z: np.ndarray) -> np.ndarray:
    # Only Hosny medians depend on these drift coordinates. Direct extraction
    # is algebraically identical to filtering data_residuals, and avoids
    # repeatedly iterating 510 unrelated population MSD targets.
    pred = fit.analytic(z)["Hosny_baseline"]
    return np.asarray([(pred[o["group"]]["quantiles"][o["metric"]][1] - o["obs"]) / o["width"]
                       for o in fit.FIT_OBS["Hosny_baseline"] if o.get("quantile") == 1])


def run_drift_pair(z: np.ndarray, initial: np.ndarray | None = None) -> tuple[np.ndarray, object]:
    zstart = z.copy()
    if initial is not None:
        zstart[20:22] = initial
    objective = lambda drift: hosny_residual(np.r_[zstart[:20], drift])
    sol = least_squares(objective, zstart[20:22],
                        bounds=(fit.LOW[20:22], fit.HIGH[20:22]),
                        max_nfev=600, ftol=1e-10, xtol=1e-10, gtol=1e-8,
                        # Close over fixed trajectory parameters; only two
                        # signed Hosny drift coordinates move here.
                        )
    zout = z.copy(); zout[20:22] = sol.x
    return zout, sol


def surface(smoke: bool, base: np.ndarray, base_scores: dict, out_rows: list) -> dict:
    pars = base[:20].reshape(5, 2, 2)
    geo = np.exp(pars.mean(axis=0))
    fit_tau_ratio = float(geo[1, 0] / geo[0, 0])
    fit_sigma_ratio = float(geo[1, 1] / geo[0, 1])
    plus = pars[:, 1, :]
    tau_shift = (max(np.log(.5 / fit_tau_ratio), float(np.log(1e-6) - plus[:, 0].min())),
                 min(np.log(30. / fit_tau_ratio), float(np.log(12. - 1e-8) - plus[:, 0].max())))
    sigma_shift = (max(np.log(.05 / fit_sigma_ratio), float(np.log(1e-4) - plus[:, 1].min())),
                   min(np.log(2. / fit_sigma_ratio), float(np.log(1e8) - plus[:, 1].max())))
    n = 9 if smoke else 61
    taus = np.unique(np.r_[np.geomspace(fit_tau_ratio * np.exp(tau_shift[0]),
                                        fit_tau_ratio * np.exp(tau_shift[1]), n), 1., fit_tau_ratio])
    sigmas = np.unique(np.r_[np.geomspace(fit_sigma_ratio * np.exp(sigma_shift[0]),
                                          fit_sigma_ratio * np.exp(sigma_shift[1]), n), 1., fit_sigma_ratio])
    total = len(taus) * len(sigmas)
    for ix, tau_ratio in enumerate(taus):
        for iy, sigma_ratio in enumerate(sigmas):
            dt, ds = np.log(tau_ratio / fit_tau_ratio), np.log(sigma_ratio / fit_sigma_ratio)
            z = base.copy()
            shifted = pars.copy(); shifted[:, 1, 0] += dt; shifted[:, 1, 1] += ds
            z[:20] = shifted.ravel()
            fixed = scores(z)
            # Optimize drift residuals against all Hosny median rows at these
            # fixed tau/noise values. Regularizers do not involve drift.
            z_optimized, sol = run_drift_pair(z)
            optimized = scores(z_optimized)
            out_rows.append({
                "tau_ratio": float(tau_ratio), "sigma_ratio": float(sigma_ratio),
                "fixed_drift_reference": float(z[20]), "fixed_drift_wgd": float(z[21]),
                "fixed_drift_delta_data": fixed["data_score"] - base_scores["data_score"],
                "fixed_drift_delta_total": fixed["total_score"] - base_scores["total_score"],
                "optimized_drift_reference": float(z_optimized[20]),
                "optimized_drift_wgd": float(z_optimized[21]),
                "optimized_drift_delta": float(z_optimized[21] - z_optimized[20]),
                "optimized_drift_delta_data": optimized["data_score"] - base_scores["data_score"],
                "optimized_drift_delta_total": optimized["total_score"] - base_scores["total_score"],
                "z": json.dumps(z_optimized.tolist()),
                "success": bool(sol.success), "nfev": int(sol.nfev),
                "status": "valid" if sol.success else "optimizer_failed"})
    return {"fit_tau_ratio": fit_tau_ratio, "fit_sigma_ratio": fit_sigma_ratio,
            "tau_ratio_bounds": [float(taus[0]), float(taus[-1])],
            "sigma_ratio_bounds": [float(sigmas[0]), float(sigmas[-1])],
            "grid_points": total}


def packed_start(base: np.ndarray, kind: str, value: float, warm: np.ndarray | None) -> tuple[np.ndarray, np.ndarray, np.ndarray, callable]:
    seed = base if warm is None else warm
    if kind == "mu_reference":
        # Hold mu_ref=value; optimize all 20 log parameters and mu_WGD.
        q0 = np.r_[seed[:20], seed[21]]
        lo = np.r_[fit.LOW[:20], fit.LOW[21]]; hi = np.r_[fit.HIGH[:20], fit.HIGH[21]]
        def unpack(q):
            z = seed.copy(); z[:20] = q[:20]; z[20] = value; z[21] = q[20]; return z
    elif kind == "mu_wgd":
        q0 = np.r_[seed[:20], seed[20]]
        lo = np.r_[fit.LOW[:20], fit.LOW[20]]; hi = np.r_[fit.HIGH[:20], fit.HIGH[20]]
        def unpack(q):
            z = seed.copy(); z[:20] = q[:20]; z[20] = q[20]; z[21] = value; return z
    else:
        # delta=mu_WGD-mu_ref; optimize c=(mu_ref+mu_WGD)/2.
        c_lo = max(-200. - value / 2., -200. + value / 2.)
        c_hi = min(200. - value / 2., 200. + value / 2.)
        c0 = float((seed[20] + seed[21]) / 2.)
        q0 = np.r_[seed[:20], np.clip(c0, c_lo + 1e-10, c_hi - 1e-10)]
        lo = np.r_[fit.LOW[:20], c_lo]; hi = np.r_[fit.HIGH[:20], c_hi]
        def unpack(q):
            z = seed.copy(); z[:20] = q[:20]
            z[20], z[21] = q[20] - value / 2., q[20] + value / 2.
            return z
    return q0, lo, hi, unpack


def profile_residual(q: np.ndarray, unpack) -> np.ndarray:
    return fit.residual(unpack(q), "primary")


def profile_point(kind: str, value: float, base: np.ndarray, warm: np.ndarray | None,
                  rng: np.random.Generator, starts_n: int) -> tuple[list[dict], np.ndarray]:
    seed_vectors = [base.copy()]
    if warm is not None:
        seed_vectors.append(warm.copy())
    while len(seed_vectors) < starts_n:
        pert = base.copy()
        pert[:20] += rng.normal(0, .4, 20)
        pert[20:22] += rng.normal(0, 1.0, 2)
        seed_vectors.append(np.clip(pert, fit.LOW + 1e-8, fit.HIGH - 1e-8))
    records = []
    solutions = []
    for start_id, seed in enumerate(seed_vectors):
        q0, lo, hi, unpack = packed_start(base, kind, value, seed)
        q0 = np.clip(q0, lo + 1e-10, hi - 1e-10)
        sol = least_squares(profile_residual, q0, args=(unpack,), bounds=(lo, hi),
                            max_nfev=1200, ftol=1e-10, xtol=1e-10, gtol=1e-8)
        z = unpack(sol.x)
        sc = scores(z)
        pars = np.exp(z[:20].reshape(5, 2, 2))
        effects = pars[:, 1, :] / pars[:, 0, :]
        pooled = fit.hyperparameters(z)
        records.append({"profile_parameter": kind, "fixed_value": float(value),
                        "data_score": sc["data_score"], "pooling_penalty": sc["pooling_penalty"],
                        "persistence_penalty": sc["persistence_penalty"], "total_score": sc["total_score"],
                        "mu_reference": float(z[20]), "mu_wgd": float(z[21]),
                        "drift_contrast": float(z[21] - z[20]),
                        "hosny_tau_reference": float(pars[4, 0, 0]),
                        "hosny_tau_wgd": float(pars[4, 1, 0]),
                        "hosny_sigma_reference": float(pars[4, 0, 1]),
                        "hosny_sigma_wgd": float(pars[4, 1, 1]),
                        "tau_pooled_fold": float(pooled["tau"]["population_wgd_fold"]),
                        "sigma_pooled_fold": float(pooled["sigma"]["population_wgd_fold"]),
                        "hosny_tau_fold": float(effects[4, 0]),
                        "hosny_sigma_fold": float(effects[4, 1]),
                        "z": json.dumps(z.tolist()),
                        "success": bool(sol.success), "nfev": int(sol.nfev),
                        "optimality": float(sol.optimality), "start_id": start_id,
                        "selected_start": False})
        solutions.append(z)
    best_score = min(r["total_score"] for r in records)
    candidates = [i for i, r in enumerate(records) if r["success"]]
    if not candidates:
        raise RuntimeError(f"All multistart profile fits failed for {kind}={value}")
    best_idx = min(candidates, key=lambda i: records[i]["total_score"])
    for i, r in enumerate(records):
        r["selected_start"] = i == best_idx
    return records, solutions[best_idx]


def run_profiles(smoke: bool, base: np.ndarray, base_total: float) -> tuple[list[dict], dict]:
    rng = np.random.default_rng(202609290)
    profiles = []
    warm_by_kind = {k: None for k in ("mu_reference", "mu_wgd", "drift_contrast")}
    fit_mu0, fit_mu1 = map(float, base[20:22])
    # Broad initial points; extend outwards only if loss <9 at either end.
    initial_specs = {
        "mu_reference": np.linspace(-30., 30., 7),
        "mu_wgd": np.linspace(-30., 30., 7),
        "drift_contrast": np.linspace(-60., 60., 7),
    }
    initial_specs["mu_reference"] = np.unique(np.r_[initial_specs["mu_reference"], fit_mu0])
    initial_specs["mu_wgd"] = np.unique(np.r_[initial_specs["mu_wgd"], fit_mu1])
    initial_specs["drift_contrast"] = np.unique(np.r_[initial_specs["drift_contrast"], fit_mu1 - fit_mu0])
    starts_n = 2 if smoke else 3
    range_meta = {}
    for kind, pilot_values in initial_specs.items():
        # Coarse pilot locates the sampled region under Δ9. Then refine that
        # bracket densely, retaining zero and the unrestricted fit coordinate.
        if smoke:
            values = pilot_values
        else:
            pilot_scores = []
            warm = None
            for value in pilot_values:
                one, warm = profile_point(kind, float(value), base, warm, rng, 1)
                pilot_scores.append(one[0])
            low_i = [i for i, r in enumerate(pilot_scores)
                     if r["total_score"] - base_total <= 9.]
            if not low_i:
                raise RuntimeError(f"Pilot profile has no sampled Δ9 region for {kind}")
            lo_i, hi_i = max(0, min(low_i) - 1), min(len(pilot_values) - 1, max(low_i) + 1)
            dense = np.linspace(pilot_values[lo_i], pilot_values[hi_i], 31)
            # Keep the zero and unrestricted-fit landmarks exact, whether or
            # not they fall in the bracket found by the coarse pilot.
            landmark = {"mu_reference": [0., fit_mu0], "mu_wgd": [0., fit_mu1],
                        "drift_contrast": [0., fit_mu1 - fit_mu0]}[kind]
            values = np.unique(np.r_[dense, landmark])
        all_records = []
        warm = None
        for value in values:
            records, warm = profile_point(kind, float(value), base, warm, rng, starts_n)
            all_records.extend(records)
        scored = [r for r in all_records if r["selected_start"]]
        profbest = min(r["total_score"] for r in scored)
        low = [float(r["fixed_value"]) for r in scored
               if r["total_score"] - base_total <= 9.]
        for r in all_records:
            r["delta_from_best_unrestricted"] = r["total_score"] - base_total
        profiles.extend(all_records)
        range_meta[kind] = {"sampled_values": [float(x) for x in values],
                            "sampled_delta9_bracket": [min(low), max(low)] if low else None,
                            "best_profile_score": profbest,
                            "minimum_delta_vs_unrestricted": profbest - base_total,
                            "sampled_profile_sensitivity_only": True,
                            "descriptive_delta_thresholds": [1., 4., 9.]}
    return profiles, range_meta


def no_hosny_residual(q: np.ndarray, base: np.ndarray) -> np.ndarray:
    z = base.copy(); z[:16] = q
    r, rows = fit.data_residuals(fit.analytic(z))
    mask = np.asarray([row["study"] != "Hosny_baseline" for row in rows], dtype=bool)
    logs = z[:16].reshape(4, 2, 2)[:, :, :]
    baseline = logs[:, 0, :]; effect = logs[:, 1, :] - baseline
    pooling = np.r_[((baseline - baseline.mean(axis=0)) / fit.BASE_SD).ravel(),
                    ((effect - effect.mean(axis=0)) / fit.EFFECT_SD).ravel()]
    prior = ((logs[:, :, 0] - np.log(fit.PRIOR_MINUTES["primary"] / 60.)).ravel()
             / fit.PRIOR_LOG_SD)
    return np.r_[r[mask], pooling, prior]


def leaveout(smoke: bool, base: np.ndarray, full_scores: dict) -> tuple[list[dict], dict]:
    rows = []
    full_hyper = fit.hyperparameters(base)
    fullpars = np.exp(base[:20].reshape(5, 2, 2))
    studies = fit.STUDIES[:4]
    # Full fit, all studies: pooled hyperparameter estimate and full objective.
    rows.append({"analysis": "full5_allstudies", "study": "Pooled",
                 "population_tau_fold": full_hyper["tau"]["population_wgd_fold"],
                 "population_sigma_fold": full_hyper["sigma"]["population_wgd_fold"],
                 **full_scores, "z": json.dumps(base.tolist())})
    # Same full-fit parameters, scored/pooled on only the four non-Hosny studies.
    remaining4 = np.asarray(base[:16]).reshape(4, 2, 2)
    r4, obs4 = fit.data_residuals(fit.analytic(base))
    m4 = np.asarray([r["study"] != "Hosny_baseline" for r in obs4], dtype=bool)
    baseline4 = remaining4[:, 0, :]; effect4 = remaining4[:, 1, :] - baseline4
    pool4 = np.r_[((baseline4 - baseline4.mean(axis=0)) / fit.BASE_SD).ravel(),
                  ((effect4 - effect4.mean(axis=0)) / fit.EFFECT_SD).ravel()]
    prior4 = ((remaining4[:, :, 0] - np.log(10. / 60.)).ravel() / fit.PRIOR_LOG_SD)
    score4_data = float(r4[m4] @ r4[m4]); score4_pool = float(pool4 @ pool4)
    score4_prior = float(prior4 @ prior4)
    score4 = {"data_score": score4_data, "pooling_penalty": score4_pool,
              "persistence_penalty": score4_prior,
              "total_score": score4_data + score4_pool + score4_prior}
    pooled_tau4, pooled_sigma4 = np.exp(effect4.mean(axis=0))
    for i, study in enumerate(studies):
        rows.append({"analysis": "full5_remaining4_only", "study": study,
                     "tau_fold": float(fullpars[i, 1, 0] / fullpars[i, 0, 0]),
                     "sigma_fold": float(fullpars[i, 1, 1] / fullpars[i, 0, 1]),
                     "population_tau_fold": float(pooled_tau4),
                     "population_sigma_fold": float(pooled_sigma4),
                     **score4, "z": json.dumps(base[:16].tolist())})
    rows.append({"analysis": "full5_remaining4_only", "study": "Pooled",
                 "population_tau_fold": float(pooled_tau4),
                 "population_sigma_fold": float(pooled_sigma4), **score4,
                 "z": json.dumps(base[:16].tolist())})
    nstarts = 2 if smoke else 8
    optima = []
    for start_id in range(nstarts):
        q0 = base[:16].copy()
        if start_id:
            rng = np.random.default_rng(43000 + start_id)
            q0 += rng.normal(0, .7, 16)
        lo, hi = fit.LOW[:16], fit.HIGH[:16]
        sol = least_squares(no_hosny_residual, np.clip(q0, lo + 1e-10, hi - 1e-10),
                            args=(base,), bounds=(lo, hi), max_nfev=1200,
                            ftol=1e-10, xtol=1e-10, gtol=1e-8)
        z = base.copy(); z[:16] = sol.x
        pred_resid, obs_rows = fit.data_residuals(fit.analytic(z))
        mask = np.asarray([r["study"] != "Hosny_baseline" for r in obs_rows], dtype=bool)
        p = z[:16].reshape(4, 2, 2)
        baseline = p[:, 0, :]; effect = p[:, 1, :] - baseline
        pool_r = np.r_[((baseline - baseline.mean(axis=0)) / fit.BASE_SD).ravel(),
                       ((effect - effect.mean(axis=0)) / fit.EFFECT_SD).ravel()]
        prior_r = ((p[:, :, 0] - np.log(10. / 60.)).ravel() / fit.PRIOR_LOG_SD)
        ds = float(pred_resid[mask] @ pred_resid[mask]); ps = float(pool_r @ pool_r)
        prs = float(prior_r @ prior_r)
        optima.append({"z": z, "sol": sol, "scores": {"data_score": ds,
                       "pooling_penalty": ps, "persistence_penalty": prs,
                       "total_score": ds + ps + prs}})
    successful = [o for o in optima if o["sol"].success]
    if not successful:
        raise RuntimeError("All leave-Hosny-out refit starts failed")
    best = min(successful, key=lambda o: o["scores"]["total_score"])
    selected_start = min((i for i, o in enumerate(optima) if o["sol"].success),
                         key=lambda i: optima[i]["scores"]["total_score"])
    zbest = best["z"]; pbest = np.exp(zbest[:16].reshape(4, 2, 2))
    base4 = np.log(pbest)
    eff = base4[:, 1, :] - base4[:, 0, :]
    poptau, popsigma = np.exp(eff.mean(axis=0))
    for i, study in enumerate(studies):
        rows.append({"analysis": "leaveHosnyOut_refit", "study": study,
                     "tau_fold": float(pbest[i, 1, 0] / pbest[i, 0, 0]),
                     "sigma_fold": float(pbest[i, 1, 1] / pbest[i, 0, 1]),
                     "population_tau_fold": float(poptau), "population_sigma_fold": float(popsigma),
                     **best["scores"], "success": bool(best["sol"].success),
                     "nfev": int(best["sol"].nfev), "selected_start": selected_start,
                     "z": json.dumps(zbest[:16].tolist())})
    rows.append({"analysis": "leaveHosnyOut_refit", "study": "Pooled",
                 "population_tau_fold": float(poptau), "population_sigma_fold": float(popsigma),
                 **best["scores"], "success": bool(best["sol"].success),
                 "nfev": int(best["sol"].nfev), "selected_start": selected_start,
                 "z": json.dumps(zbest[:16].tolist())})
    summary = {"starts": [{"start_id": i, "total_score": o["scores"]["total_score"],
                           "success": bool(o["sol"].success), "nfev": int(o["sol"].nfev),
                           "optimality": float(o["sol"].optimality)} for i, o in enumerate(optima)],
               "best_total_score": best["scores"]["total_score"],
               "best_success": bool(best["sol"].success),
               "full5_scores": full_scores,
               "full5_same4_scores": score4,
               "refitted4_scores": best["scores"]}
    return rows, summary


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row.keys()))
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys); writer.writeheader(); writer.writerows(rows)


def summarize_results(profiles: list[dict], leave_rows: list[dict]) -> dict:
    selected = [r for r in profiles if r["selected_start"] and r["success"]]
    thresholds = {}
    for threshold in (1., 4., 9.):
        accepted = [r for r in selected if r["delta_from_best_unrestricted"] <= threshold]
        thresholds[str(int(threshold))] = {key: [float(min(r[key] for r in accepted)),
                                                   float(max(r[key] for r in accepted))]
            for key in ("hosny_tau_fold", "hosny_sigma_fold", "mu_reference", "mu_wgd", "drift_contrast")}
        thresholds[str(int(threshold))]["sampled_profile_fit_count"] = len(accepted)
    zeros = {}
    for kind in ("mu_reference", "mu_wgd", "drift_contrast"):
        row = next(r for r in selected if r["profile_parameter"] == kind
                   and np.isclose(r["fixed_value"], 0., atol=1e-12))
        zeros[kind] = {"delta_objective_vs_unrestricted": float(row["delta_from_best_unrestricted"]),
                       "total_score": float(row["total_score"])}
    pooled = {}
    for label in ("full5_allstudies", "full5_remaining4_only", "leaveHosnyOut_refit"):
        row = next(r for r in leave_rows if r["analysis"] == label and r["study"] == "Pooled")
        pooled[label] = {k: float(row[k]) for k in
                         ("population_tau_fold", "population_sigma_fold", "data_score",
                          "pooling_penalty", "persistence_penalty", "total_score")}
    near_rows = [r for r in selected if r["delta_from_best_unrestricted"] <= 9.]
    active = []
    for row in near_rows:
        z = np.asarray(json.loads(row["z"]), dtype=float)
        for i, (value, lo, hi) in enumerate(zip(z, fit.LOW, fit.HIGH)):
            if value - lo <= 1e-6:
                active.append({"profile_parameter": row["profile_parameter"], "fixed_value": float(row["fixed_value"]),
                               "coordinate": i, "bound": "lower"})
            if hi - value <= 1e-6:
                active.append({"profile_parameter": row["profile_parameter"], "fixed_value": float(row["fixed_value"]),
                               "coordinate": i, "bound": "upper"})
    return {"sampled_profile_parameter_ranges": thresholds,
            "zero_fixed_drift_profile_delta": zeros,
            "pooled_fold_comparisons": pooled,
            "active_model_bounds_among_sampled_delta9_profile_fits": active,
            "interpretation": "The unions summarize sampled one-dimensional profile fits only; they are not a joint uncertainty region or full parameter limits."}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true", help="small run to benchmark the optimizer")
    args = parser.parse_args()
    start_time = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    fit_path = RESULTS / FIT_NAME
    fit_obj = json.loads(fit_path.read_text())
    if fit_obj["total_score"] != FIT_TOTAL_SCORE or fit_obj["data_score"] != FIT_DATA_SCORE:
        raise AssertionError("selected primary fit provenance score mismatch")
    base = np.asarray(fit_obj["z"], dtype=float)
    base_scores = scores(base)
    if not (np.isclose(base_scores["total_score"], FIT_TOTAL_SCORE, rtol=0, atol=1e-9)
            and np.isclose(base_scores["data_score"], FIT_DATA_SCORE, rtol=0, atol=1e-10)):
        raise AssertionError(f"Recomputed full objective does not match selected primary fit: {base_scores}")
    direct_hosny=hosny_residual(base)
    full_r,full_rows=fit.data_residuals(fit.analytic(base))
    full_hosny=full_r[np.asarray([r['study']=='Hosny_baseline' and FIT_MASK(r) for r in full_rows])]
    if not np.allclose(direct_hosny,full_hosny,rtol=0,atol=1e-12):
        raise AssertionError('Direct Hosny drift objective differs from full residual filtering')
    test_b = np.asarray([0., .1, 1., 10.])
    rice_validation = {"test_noncentralities": test_b.tolist(),
                       "maximum_scaled_bessel_vs_scipy_mean_abs_difference": float(np.max(
                           np.abs(rice_mean_stable(test_b) - _rice_mean_original(test_b))))}
    if rice_validation["maximum_scaled_bessel_vs_scipy_mean_abs_difference"] > 1e-8:
        raise AssertionError(f"Stable Rice mean differs from SciPy for ordinary noncentralities: {rice_validation}")
    # Independent unrestricted optimizer check around the stored multistart
    # winner. Keep successful starts and require the refit not to worsen it.
    central_trials = []
    central_nstarts = 2 if args.smoke else 3
    for start_id in range(central_nstarts):
        q0 = base.copy()
        if start_id:
            rng = np.random.default_rng(52000 + start_id)
            q0[:20] += rng.normal(0, .1, 20)
            q0[20:22] += rng.normal(0, .2, 2)
        sol = least_squares(fit.residual, np.clip(q0, fit.LOW + 1e-10, fit.HIGH - 1e-10),
                            args=("primary",), bounds=(fit.LOW, fit.HIGH), max_nfev=1200,
                            ftol=1e-10, xtol=1e-10, gtol=1e-8)
        trial_score = float(sol.fun @ sol.fun)
        central_trials.append({"start_id": start_id, "score": trial_score,
                               "success": bool(sol.success), "nfev": int(sol.nfev),
                               "optimality": float(sol.optimality), "z": sol.x.tolist()})
    successful_central = [x for x in central_trials if x["success"]]
    if not successful_central:
        raise RuntimeError("All unrestricted central reoptimization starts failed")
    central_best = min(successful_central, key=lambda r: r["score"])
    if central_best["score"] > FIT_TOTAL_SCORE + 1e-7:
        raise AssertionError("Central unrestricted refit worsened the stored selected primary fit objective")
    print(json.dumps({"phase":"broad_surface", "elapsed_seconds":time.time()-start_time}),flush=True)
    surfaces = []
    surface_meta = surface(args.smoke, base, base_scores, surfaces)
    # Drift-only minimization must never worsen the fixed-drift score.
    if any(r["optimized_drift_delta_total"] > r["fixed_drift_delta_total"] + 1e-7
           for r in surfaces if r["success"]):
        raise AssertionError("Reoptimized drift worsened the total objective on the relative surface")
    print(json.dumps({"phase":"drift_profiles", "elapsed_seconds":time.time()-start_time}),flush=True)
    profiles, profile_meta = run_profiles(args.smoke, base, central_best["score"])
    print(json.dumps({"phase":"leave_hosny_out", "elapsed_seconds":time.time()-start_time}),flush=True)
    leave_rows, leave_meta = leaveout(args.smoke, base, base_scores)
    write_csv(OUT / "hosny_surface.csv", surfaces)
    write_csv(OUT / "hosny_drift_profiles.csv", profiles)
    write_csv(OUT / "hosny_leaveout_comparison.csv", leave_rows)
    # Explicitly audit every stored optimized parameter against model bounds.
    surface_bounds_ok = all(np.all(np.asarray(json.loads(r["z"])) >= fit.LOW - 1e-10)
                            and np.all(np.asarray(json.loads(r["z"])) <= fit.HIGH + 1e-10)
                            for r in surfaces)
    profile_selected = [r for r in profiles if r["selected_start"]]
    profile_bounds_ok = all(np.all(np.asarray(json.loads(r["z"])) >= fit.LOW - 1e-10)
                            and np.all(np.asarray(json.loads(r["z"])) <= fit.HIGH + 1e-10)
                            for r in profile_selected)
    leave_selected_z = np.asarray(json.loads(next(r["z"] for r in leave_rows
                                  if r["analysis"] == "leaveHosnyOut_refit" and r["study"] == "Pooled")))
    leave_bounds_ok = bool(np.all(leave_selected_z >= fit.LOW[:16] - 1e-10)
                           and np.all(leave_selected_z <= fit.HIGH[:16] + 1e-10))
    if not (surface_bounds_ok and profile_bounds_ok and leave_bounds_ok):
        raise AssertionError("One or more optimized outputs exceed their model parameter bounds")
    profile_start_spread = {}
    for kind in {r["profile_parameter"] for r in profiles}:
        ranges = []
        for value in {r["fixed_value"] for r in profiles if r["profile_parameter"] == kind}:
            same = [r["total_score"] for r in profiles if r["profile_parameter"] == kind
                    and r["fixed_value"] == value and r["success"]]
            if len(same) > 1:
                ranges.append(max(same) - min(same))
        profile_start_spread[kind] = max(ranges, default=0.)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "fit": FIT_NAME, "fit_sha256": sha(fit_path),
        "fit_scores": base_scores,
        "fixed_fit_score_validation": {"stored_data_score": FIT_DATA_SCORE,
             "stored_total_score": FIT_TOTAL_SCORE,
             "recomputed_data_score": base_scores["data_score"],
             "recomputed_total_score": base_scores["total_score"],
             "match": True},
        "rice_numeric_validation": rice_validation,
        "central_unrestricted_multistart": {"trials": central_trials,
              "selected_score": central_best["score"],
              "selected_z": central_best["z"],
              "not_worse_than_selected_primary_fit": central_best["score"] <= FIT_TOTAL_SCORE + 1e-7},
        "surface": {**surface_meta,
             "all_drift_optimizers_successful": all(r["success"] for r in surfaces),
             "max_optimized_minus_fixed_total_delta": max(r["optimized_drift_delta_total"] - r["fixed_drift_delta_total"] for r in surfaces),
             "all_fitted_parameters_within_model_bounds": surface_bounds_ok},
        "profile_optimization_audit": {"selected_profiles_successful": all(r["success"] for r in profile_selected),
             "all_selected_parameters_within_model_bounds": profile_bounds_ok,
             "max_multistart_score_spread_by_profile": profile_start_spread,
             "converged_profile_start_count": sum(r["success"] for r in profiles),
             "profile_start_count": len(profiles)},
        "leaveout_optimization_audit": {"selected_fit_successful": leave_meta["best_success"],
             "all_parameters_within_model_bounds": leave_bounds_ok},
        "profiles": profile_meta,
        "leave_hosny_out": leave_meta,
        "summary_results": summarize_results(profiles, leave_rows),
        "semantics": {"profile_ranges": "Sampled profile sensitivity from bounded multistart fits; not joint confidence intervals or full uncertainty limits.",
                      "delta_thresholds": "Changes in total objective of 1, 4, and 9 are descriptive thresholds selected before profiling.",
                      "excluded_hosny": "Hosny is omitted from data residuals, pooling centers, pooling residuals, and persistence penalties in the leave-Hosny-out refit."},
        "inputs_sha256": {str(FIT_DIR / "fit.py"): sha(FIT_DIR / "fit.py"),
                          str(FIT_DIR / "model.py"): sha(FIT_DIR / "model.py")},
        "smoke_run": bool(args.smoke), "elapsed_seconds": time.time() - start_time,
        "grid_counts": {"surface": len(surfaces), "profile_rows_all_starts": len(profiles),
                        "leaveout_rows": len(leave_rows)}
    }
    metadata["full_objective_equivalence"] = json.loads((OUT / "full_objective_equivalence.json").read_text())
    metadata["compiled_residual_source_sha256"] = sha(Path(__file__).with_name("fast_residuals.py"))
    (OUT / "hosny_analysis.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"smoke": args.smoke, "elapsed_seconds": metadata["elapsed_seconds"],
                      "surface_points": len(surfaces), "profile_rows": len(profiles),
                      "leaveout_rows": len(leave_rows), "out": str(OUT)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
