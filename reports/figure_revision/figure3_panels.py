"""Reusable data and panel helpers for revised Figure 3.

The module consumes the source observations and validated analytic fit outputs;
it does not perform or change any optimization.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import LogLocator, FuncFormatter

ROOT = Path(__file__).resolve().parents[2]
FIT_DIR = ROOT / "wgd_langevin" / "gaussian_hierarchy"
RESULTS_DIR = FIT_DIR / "results"
FIGURE_DIR = ROOT / "reports" / "figure_revision" / "figures" / "figure3_review"
OBSERVATIONS_CSV = ROOT / "wgd_langevin" / "data" / "observations.csv"
HOSNY_DATA = ROOT / "wgd_langevin" / "hosny" / "data.py"
UPDATED_OBSERVATIONS_CSV = FIGURE_DIR / "figure3_observations_updated.csv"
UPDATED_OBSERVATIONS_PROVENANCE = FIGURE_DIR / "figure3_observations_updated_provenance.json"
RELATIVE_SURFACE_CSV = FIGURE_DIR / "figure3_relative_persistence_noise_surface.csv"

sys.path.insert(0, str(FIT_DIR))
spec = importlib.util.spec_from_file_location("figure3_gaussian_fit", FIT_DIR / "fit.py")
fit = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = fit
spec.loader.exec_module(fit)

STUDY_LABELS = {
    "Xuan2018": "Xuan 2018",
    "Mallin2023": "Mallin 2023",
    "Jemaa2017": "Jemaa 2017",
    "Jemaa2023": "Jemaa 2023",
    "Hosny_baseline": "Hosny",
}
STUDY_ORDER = list(STUDY_LABELS)
METRIC_ORDER = {
    "Xuan2018": ["speed", "msd"],
    "Mallin2023": ["path_length", "directness"],
    "Jemaa2017": ["speed"],
    "Jemaa2023": ["speed"],
    "Hosny_baseline": ["path_length", "directness", "parallel"],
}
METRIC_LABELS = {"speed": "speed", "msd": "MSD", "path_length": "path length",
                 "directness": "directness", "parallel": "signed migration index"}
NAVY = "#19394d"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fit_score(record: dict, z: np.ndarray, scenario: str) -> dict:
    resid, rows = fit.data_residuals(fit.analytic(z))
    mask = np.asarray([r["study"] != "Hosny_baseline" or r.get("quantile") == 1 for r in rows])
    data = float(resid[mask] @ resid[mask])
    pool = fit.pooling_residuals(z, fit.SCENARIOS[scenario])
    prior = fit.persistence_residuals(z, scenario)
    components = {"data_score": data, "pooling_penalty": float(pool @ pool),
                  "persistence_prior_penalty": float(prior @ prior)}
    components["total_score"] = sum(components.values())
    for key, value in components.items():
        if not np.isclose(value, float(record[key]), rtol=0, atol=1e-8):
            raise ValueError(f"Stored {scenario} fit {key} does not validate: {value} != {record[key]}")
    return components


def load_fit_scenarios(root: Path | None = None) -> dict[str, dict]:
    """Load each best successful scenario fit and validate scores and rows."""
    results = (root / "wgd_langevin/gaussian_hierarchy/results") if root else RESULTS_DIR
    candidates = []
    for path in sorted(results.glob("fit_*.json")):
        record = json.loads(path.read_text())
        if record.get("success") and record.get("scenario") in fit.SCENARIOS:
            candidates.append((path, record))
    best = {}
    for scenario in fit.SCENARIOS:
        available = [(p, r) for p, r in candidates if r["scenario"] == scenario]
        if not available:
            raise FileNotFoundError(f"No successful fit found for scenario {scenario}")
        path, record = min(available, key=lambda pair: pair[1]["total_score"])
        z = np.asarray(record["z"], dtype=float)
        _fit_score(record, z, scenario)
        resid, rows = fit.data_residuals(fit.analytic(z))
        mask = np.asarray([r["study"] != "Hosny_baseline" or r.get("quantile") == 1 for r in rows])
        expected = [r for r, keep in zip(rows, mask) if keep]
        stored = record.get("rows", [])
        if len(expected) != len(stored):
            raise ValueError(f"Stored residual-row count for {path.name} does not match source observations")
        for a, b in zip(expected, stored):
            if (a["study"], a["group"], a["metric"], a.get("quantile")) != (
                    b["study"], b["group"], b["metric"], b.get("quantile")):
                raise ValueError(f"Stored row mapping mismatch in {path.name}")
            if not np.isclose(a["obs"], b["obs"], rtol=0, atol=1e-12):
                raise ValueError(f"Stored observation mismatch in {path.name}: {a} != {b}")
            if not np.isclose(a["predicted"], b["predicted"], rtol=1e-9, atol=1e-9):
                raise ValueError(f"Stored prediction mismatch in {path.name}: {a} != {b}")
        best[scenario] = {"record": record, "path": path, "z": z,
                          "scores": _fit_score(record, z, scenario)}
    return best


def _source_observation_maps() -> tuple[dict, dict]:
    fit_obs = {}
    full_obs = {}
    for study, observations in fit.FIT_OBS.items():
        for row in observations:
            key = (study, int(row["group"]), row["metric"])
            if study == "Hosny_baseline":
                quantile = int(row["quantile"])
                full_obs.setdefault(key, {})[quantile] = row
                if quantile == 1:
                    fit_obs[key] = row
            else:
                fit_obs[key] = row
                full_obs.setdefault(key, {})[1] = row
    return fit_obs, full_obs


def build_observation_rows(primary_fit: dict | None = None, root: Path | None = None) -> tuple[list[dict], dict]:
    """Build the nine point comparisons from source observations and fit_05.

    The observed fold is the WGD+ observed central value divided by WGD−.
    Hosny uses each group's median, as in the objective; its group quartiles
    are carried in separate columns and never converted to a ratio interval.
    Source tolerance bounds are retained separately, not depicted as errors.
    """
    if root:
        candidate_dir = root / "wgd_langevin/gaussian_hierarchy/results"
        fit_records = []
        for p in sorted(candidate_dir.glob("fit_*.json")):
            r = json.loads(p.read_text())
            if r.get("scenario") == "primary" and r.get("success"):
                fit_records.append((p, r))
        if not fit_records:
            raise FileNotFoundError("No successful primary Gaussian-hierarchy fit")
        path, record = min(fit_records, key=lambda pair: pair[1]["total_score"])
    elif primary_fit is None:
        path = RESULTS_DIR / "fit_05.json"
        record = json.loads(path.read_text())
    else:
        path = primary_fit.get("path") if isinstance(primary_fit, dict) else None
        record = primary_fit["record"] if "record" in primary_fit else primary_fit
    if record.get("scenario") != "primary" or not record.get("success"):
        raise ValueError("Observed/predicted panel requires a successful primary fit")
    z = np.asarray(record["z"], dtype=float)
    _fit_score(record, z, "primary")
    predicted, prediction_rows = fit.data_residuals(fit.analytic(z))
    # Source observations come from observations.csv for four non-Hosny
    # studies and data.py for Hosny, loaded through the model's canonical reader.
    selected_obs, all_obs = _source_observation_maps()
    stored_rows = {(r["study"], int(r["group"]), r["metric"], r.get("quantile")): r
                   for r in record["rows"]}
    pred_rows = {(r["study"], int(r["group"]), r["metric"], r.get("quantile")): r
                 for r in prediction_rows
                 if r["study"] != "Hosny_baseline" or r.get("quantile") == 1}
    rows = []
    for study in STUDY_ORDER:
        for metric in METRIC_ORDER[study]:
            k0, k1 = (study, 0, metric), (study, 1, metric)
            o0, o1 = selected_obs[k0], selected_obs[k1]
            q0, q1 = all_obs[k0], all_obs[k1]
            fit_key0 = (study, 0, metric, 1 if study == "Hosny_baseline" else None)
            fit_key1 = (study, 1, metric, 1 if study == "Hosny_baseline" else None)
            source_row0 = pred_rows[fit_key0]
            source_row1 = pred_rows[fit_key1]
            stored0, stored1 = stored_rows[fit_key0], stored_rows[fit_key1]
            if not (np.isclose(o0["obs"], stored0["obs"], rtol=0, atol=1e-12)
                    and np.isclose(o1["obs"], stored1["obs"], rtol=0, atol=1e-12)):
                raise ValueError(f"Primary fit and source data disagree for {study}/{metric}")
            pred0, pred1 = float(source_row0["predicted"]), float(source_row1["predicted"])
            obs0, obs1 = float(o0["obs"]), float(o1["obs"])
            row = {
                "study": study, "study_label": STUDY_LABELS[study], "metric": metric,
                "metric_label": METRIC_LABELS[metric],
                "observed_WGDminus": obs0, "observed_WGDplus": obs1,
                "observed_fold_WGDplus_over_WGDminus": obs1 / obs0,
                "predicted_WGDminus": pred0, "predicted_WGDplus": pred1,
                "predicted_fold_WGDplus_over_WGDminus": pred1 / pred0,
                "log2_predicted_over_observed_fold": float(np.log2((pred1 / pred0) / (obs1 / obs0))),
                "fit_convention": "group median" if study == "Hosny_baseline" else "source fit-role estimate",
                "Hosny_fit_quantile": 1 if study == "Hosny_baseline" else "",
                "source_fit_residual_reference": float(stored0["residual"]),
                "source_fit_residual_WGDplus": float(stored1["residual"]),
                "reference_source_low": o0.get("low", ""), "reference_source_high": o0.get("high", ""),
                "WGDplus_source_low": o1.get("low", ""), "WGDplus_source_high": o1.get("high", ""),
                "reference_source_width": o0.get("width", ""), "WGDplus_source_width": o1.get("width", ""),
                "reference_q1": q0.get(0, {}).get("obs", ""),
                "reference_q3": q0.get(2, {}).get("obs", ""),
                "WGDplus_q1": q1.get(0, {}).get("obs", ""),
                "WGDplus_q3": q1.get(2, {}).get("obs", ""),
                "ratio_error_bars": False,
                "source_mapping_reference": f"{study}:group0:{metric}",
                "source_mapping_WGDplus": f"{study}:group1:{metric}",
            }
            rows.append(row)
    if len(rows) != 9:
        raise AssertionError(f"Expected nine study/metric comparisons, got {len(rows)}")
    # The restored Jemaa speed rows must be present and derive from source
    # fit-role observations, not any prior figure CSV.
    for study in ("Jemaa2017", "Jemaa2023"):
        if not any(r["study"] == study and r["metric"] == "speed" for r in rows):
            raise AssertionError(f"Missing observed-vs-predicted speed fold for {study}")
    provenance = {
        "selected_fit": path.name if 'path' in locals() and path else "fit_05.json",
        "selected_scenario": "primary", "selected_total_score": float(record["total_score"]),
        "selected_data_score": float(record["data_score"]),
        "fit_validation": _fit_score(record, z, "primary"),
        "row_count": len(rows),
        "row_mappings": [{"row": i, "study": r["study"], "metric": r["metric"],
                           "source_mapping_reference": r["source_mapping_reference"],
                           "source_mapping_WGDplus": r["source_mapping_WGDplus"],
                           "Hosny_fit_quantile": r["Hosny_fit_quantile"]}
                          for i, r in enumerate(rows)],
        "observation_semantics": "Non-Hosny values use observations.csv rows whose role is fit. Hosny values use source quartile arrays and the median only for fitting. Parallel is the signed directional-displacement/path-length index. Quartiles and source tolerances are exported separately; no fold interval is inferred or plotted.",
        "source_sha256": {str(p.relative_to(ROOT)): _sha(p) for p in
                          (path, FIT_DIR / "fit.py", FIT_DIR / "model.py", OBSERVATIONS_CSV,
                           HOSNY_DATA) if p is not None},
    }
    return rows, provenance


def write_observation_table(root: Path | None = None) -> tuple[Path, Path]:
    scenarios = load_fit_scenarios(root=root)
    rows, provenance = build_observation_rows(primary_fit=scenarios["primary"], root=root)
    provenance["all_scenario_fits"] = {
        name: {"fit": value["path"].name, "sha256": _sha(value["path"]),
              "data_score": float(value["record"]["data_score"]),
              "total_score": float(value["record"]["total_score"])}
        for name, value in scenarios.items()}
    target_csv = (root / "reports/figure_revision/figures/figure3_review/figure3_observations_updated.csv") if root else UPDATED_OBSERVATIONS_CSV
    target_meta = target_csv.with_name("figure3_observations_updated_provenance.json")
    target_csv.parent.mkdir(parents=True, exist_ok=True)
    with target_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    provenance["output_sha256"] = {target_csv.name: _sha(target_csv)}
    target_meta.write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return target_csv, target_meta


def draw_observed_predicted(ax, rows: list[dict], title: str | None = None,
                            observed_color: str = "#5b7380",
                            predicted_color: str = "#d36c32"):
    """Draw all nine study/metric fold pairs on one shared logarithmic axis."""
    lookup = {(r["study"], r["metric"]): r for r in rows}
    order = [(study, metric) for study in STUDY_ORDER for metric in METRIC_ORDER[study]]
    y_positions = np.arange(len(order), dtype=float)
    for idx, (study, metric) in enumerate(order):
        r = lookup[(study, metric)]
        yobs, ypred = idx + .12, idx - .12
        observed = float(r["observed_fold_WGDplus_over_WGDminus"])
        predicted = float(r["predicted_fold_WGDplus_over_WGDminus"])
        ax.plot([observed, predicted], [yobs, ypred], color="#aebbc1", lw=.8, zorder=1)
        ax.scatter([observed], [yobs], marker="o", s=28, facecolor=observed_color,
                   edgecolor="white", linewidth=.45, zorder=3)
        ax.scatter([predicted], [ypred], marker="D", s=25, facecolor=predicted_color,
                   edgecolor="white", linewidth=.45, zorder=4)
    labels = [f"{STUDY_LABELS[s]} · {METRIC_LABELS[m]}" for s, m in order]
    ax.set_yticks(y_positions, labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.axvline(1., color="#33454e", linestyle="--", linewidth=.85, zorder=0)
    ax.set_xlabel("WGD+ / WGD− fold")
    if title:
        ax.set_title(title, loc="left", color=NAVY)
    ax.grid(axis="x", which="major", color="#dce4e8", linewidth=.6)
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1., 2., 5.)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.tick_params(axis="y", labelsize=7.6)
    ax.tick_params(axis="x", labelsize=7.2)
    ax.set_ylim(len(order) - .4, -.6)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=observed_color,
                              markeredgecolor="white", label="Observed fold", markersize=6),
                       Line2D([0], [0], marker="D", color="none", markerfacecolor=predicted_color,
                              markeredgecolor="white", label="Predicted fold", markersize=5.5),
                       Line2D([0], [0], color="#33454e", linestyle="--", label="Equal groups")],
              frameon=False, ncol=3, fontsize=7, loc="lower right")
    return {f"{s}:{m}": float(y) for y, (s, m) in zip(y_positions, order)}


def draw_relative_surface(ax, csv_path: Path | None = None,
                          score_column: str = "delta_data_misfit",
                          norm: Normalize | None = None, add_colorbar: bool = False,
                          title: str | None = None, cmap: str = "RdBu_r"):
    """Draw the saved WGD−-relative score surface; performs no refit."""
    path = csv_path or RELATIVE_SURFACE_CSV
    with Path(path).open(newline="") as f:
        rows = list(csv.DictReader(f))
    rows = [r for r in rows if r.get("status") == "valid" and r.get(score_column, "") != ""]
    tau = np.asarray(sorted({float(r["tau_WGDplus_over_WGDminus"]) for r in rows}), dtype=float)
    sigma = np.asarray(sorted({float(r["sigma_WGDplus_over_WGDminus"]) for r in rows}), dtype=float)
    if len(rows) != len(tau) * len(sigma):
        raise ValueError(f"Saved relative surface is not a complete rectangular grid: {path}")
    grid = np.full((len(tau), len(sigma)), np.nan)
    ix = {v: i for i, v in enumerate(tau)}; iy = {v: i for i, v in enumerate(sigma)}
    for r in rows:
        grid[ix[float(r["tau_WGDplus_over_WGDminus"])],
             iy[float(r["sigma_WGDplus_over_WGDminus"])]] = float(r[score_column])
    if not np.isfinite(grid).all():
        raise ValueError(f"Relative surface contains non-finite displayed scores: {score_column}")
    if norm is None:
        lim = float(np.max(np.abs(grid)))
        norm = Normalize(-lim, lim)
    mesh = ax.pcolormesh(tau, sigma, grid.T, shading="auto", cmap=cmap, norm=norm, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.scatter([1.], [1.], marker="o", s=37, facecolor="white", edgecolor="#242f35",
               linewidth=.8, zorder=5)
    ax.scatter([15.158309242829992], [.1381073764467221], marker="*", s=105,
               facecolor="#111111", edgecolor="white", linewidth=.65, zorder=6)
    ax.set_xlabel("WGD+ / WGD− persistence fold")
    ax.set_ylabel("WGD+ / WGD− noise fold")
    if title:
        ax.set_title(title, loc="left", color=NAVY)
    ax.grid(which="major", color="white", linewidth=.55, alpha=.52)
    if add_colorbar:
        cb = ax.figure.colorbar(mesh, ax=ax, fraction=.045, pad=.035)
        cb.set_label(score_column.replace("_", " "))
    return mesh


def main() -> None:
    """Export current observation folds from the validated primary fit."""
    fits = load_fit_scenarios()
    primary = fits["primary"]
    if primary["path"].name != "fit_05.json":
        raise ValueError(f"Expected fit_05 primary optimum, found {primary['path'].name}")
    csv_path, provenance_path = write_observation_table()
    print(json.dumps({"scenario_fits": {k: v["path"].name for k, v in fits.items()},
                      "observed_predicted_rows": 9,
                      "csv": str(csv_path), "provenance": str(provenance_path)}, indent=2))


if __name__ == "__main__":
    main()
