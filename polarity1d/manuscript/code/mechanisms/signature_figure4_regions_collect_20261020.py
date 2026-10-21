"""Whole-seed A5 aggregation; native movies and per-seed tracks stay immutable.

Survival is measured for incident observed contour components. Unresolved
amounts and identity loss remain visible rather than becoming zero takeovers.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import warnings

import numpy as np

from signature_figure4_regions_campaign_20261020 import sha, save


def quantiles(values):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanquantile(values, [.025, .975], axis=0)


def average(values, weights):
    values = np.asarray(values, float)
    valid = np.isfinite(values)
    denominator = weights @ valid.astype(float)
    return np.divide(weights @ np.where(valid, values, 0.), denominator,
        out=np.full(len(weights), np.nan), where=denominator > 0)


def estimands(result):
    s = result["summary"]
    births = result["births"]
    tracks = [t for t in result["tracks"] if not t["left_censored"]]
    angles = [b["distance_from_dominant_radians"] for b in births
              if b["distance_from_dominant_radians"] is not None]
    physical = [b["distance_from_dominant_physical"] for b in births
                if b["distance_from_dominant_physical"] is not None]
    unscorable = sum(t.get("takeover_unscorable_reason") is not None
                     or t.get("previous_dominant_id") is None for t in tracks)
    answer = {name: s[name] for name in (
        "region_count_mean", "low_contrast_time_fraction", "mean_region_width",
        "mean_region_width_fraction", "births_per_time_per_perimeter")}
    answer.update(
        multiple_front_fraction=sum(v for k, v in s["region_count_time_distribution"].items() if int(k) >= 2),
        single_front_fraction=s["region_count_time_distribution"].get("1", 0.),
        observed_births=s["observed_births"],
        observed_takeovers=s["takeovers"],
        observed_takeovers_per_birth=s["takeovers"]/len(births) if births else None,
        observed_pairwise_takeovers=s["pairwise_takeovers"],
        observed_pairwise_takeovers_per_birth=s["pairwise_takeovers"]/len(births) if births else None,
        unscorable_takeover_tracks=unscorable,
        ever_unscorable_takeover_fraction=unscorable/len(tracks) if tracks else None,
        opposite_birth_fraction=np.mean([b["opposite_hemisphere_birth"] for b in births
            if b["distance_from_dominant_radians"] is not None]) if angles else None,
        mean_birth_angle=np.mean(angles) if angles else None,
        mean_birth_physical_distance=np.mean(physical) if physical else None)
    for turning in result["turning"]:
        answer[f"mean_absolute_turn_{turning['frame_interval']}_frames"] = turning["mean_absolute_angle"]
    return {k: None if v is None else float(v) for k, v in answer.items()}


def survival_arrays(results, *, clock, horizon, landmarks=()):
    samples = []
    for result in results:
        tracks = [t for t in result["tracks"] if not t["left_censored"]]
        samples.append([(float(t["duration"] if t["right_censored"] or clock == "observed"
                              else t["duration_"+clock]), not t["right_censored"]) for t in tracks])
    grid = np.unique(np.r_[0., horizon, [v for v in landmarks if 0 <= v <= horizon],
        [min(t, horizon) for sample in samples for t, event in sample]])
    risk = np.asarray([[sum(t >= at for t, event in sample) for at in grid] for sample in samples], float)
    events = np.asarray([[sum(event and t == at for t, event in sample) for at in grid] for sample in samples], float)
    censor = np.asarray([[sum(not event and t == at for t, event in sample) for at in grid] for sample in samples], float)
    return grid, risk, events, censor


def weighted_survival(risk, events, weights):
    r, d = weights @ risk, weights @ events
    factor = 1.-np.divide(d, r, out=np.zeros_like(d), where=r > 0)
    survival = np.cumprod(factor, axis=1)
    # A zero-survival tail is identified. A positive tail after the last
    # censored component is not; do not extrapolate a flat KM tail as data.
    survival[(r <= 0) & (survival > 0)] = np.nan
    return survival


def rmst_bounds(grid, curves, horizon):
    """Integrate identified survival and retain the unidentified-tail bounds."""
    curves = np.atleast_2d(curves)
    lower = np.zeros(len(curves)); upper = np.zeros(len(curves))
    for j, curve in enumerate(curves):
        if not np.isfinite(curve[0]):
            lower[j], upper[j] = np.nan, np.nan
            continue
        valid = np.flatnonzero(np.isfinite(curve))
        last = int(valid[-1])
        area = float(np.sum(curve[:last]*np.diff(grid[:last+1])))
        lower[j] = area
        upper[j] = area+curve[last]*max(0., horizon-grid[last])
    return lower, upper


def histogram_rows(values, bins, weights, context):
    hist = np.asarray([np.histogram(v, bins=bins)[0] for v in values], float)
    total = hist.sum(axis=1)
    numerator, denominator = weights @ hist, weights @ total
    boot = np.divide(numerator, denominator[:, None], out=np.full_like(numerator, np.nan),
                     where=denominator[:, None] > 0)
    point = hist.sum(axis=0)/hist.sum() if hist.sum() else np.full(len(bins)-1, np.nan)
    ci = quantiles(boot)
    return [dict(context, bin_lower=float(a), bin_upper=float(b), probability=float(point[i]),
        ci_lower=float(ci[0, i]), ci_upper=float(ci[1, i]), observed_items=int(hist[:, i].sum()),
        observed_total=int(hist.sum())) for i, (a, b) in enumerate(zip(bins[:-1], bins[1:]))]


def aggregate(results, *, metadata, draws=2000, seed=2026102005):
    n = len(results)
    if not n:
        raise ValueError("An A5 condition needs complete independent seed measurements")
    deterministic = metadata.get("kind") in ("mean", "zero")
    if deterministic and n != 1:
        raise ValueError("A deterministic matched control cannot become independent seed replicates")
    weights = (np.ones((1, n)) if deterministic else
               np.random.default_rng(seed).multinomial(n, np.ones(n)/n, size=draws))
    replication = dict(independent_seeds=0 if deterministic else n,
        deterministic_controls=n if deterministic else 0,
        bootstrap_resamples=0 if deterministic else draws,
        uncertainty_unit="deterministic matched control; no bootstrap" if deterministic else "whole independent seed")
    metrics, distributions, survival = [], [], []
    samples = [estimands(r) for r in results]
    for name in samples[0]:
        values = np.asarray([np.nan if s[name] is None else s[name] for s in samples], float)
        boot = average(values, weights)
        ci = quantiles(boot)
        valid = np.isfinite(values)
        metrics.append(dict(metadata, statistic=name, estimate=float(values[valid].mean()) if valid.any() else None,
            ci_lower=float(ci[0]), ci_upper=float(ci[1]), **replication,
            seeds_with_estimand=int(valid.sum()), finite_bootstrap_fraction=float(np.isfinite(boot).mean()),
            limited_two_seed_ensemble=not deterministic and n == 2))
    counts = sorted({int(k) for r in results for k in r["summary"]["region_count_time_distribution"]})
    for count in counts:
        values = np.asarray([r["summary"]["region_count_time_distribution"].get(str(count), 0.) for r in results])
        ci = quantiles(average(values, weights))
        distributions.append(dict(metadata, statistic="simultaneous_region_count", value=count,
            probability=float(values.mean()), ci_lower=float(ci[0]), ci_upper=float(ci[1])))
    birth_angles = [[b["distance_from_dominant_radians"] for b in r["births"]
        if b["distance_from_dominant_radians"] is not None] for r in results]
    distributions += histogram_rows(birth_angles, np.linspace(0., np.pi, 19), weights,
        dict(metadata, statistic="birth_angle_from_previous_dominant"))
    for interval in (1, 2, 4, 8):
        values = [next(t["angles_radians"] for t in r["turning"] if t["frame_interval"] == interval) for r in results]
        distributions += histogram_rows(values, np.linspace(-np.pi, np.pi, 37), weights,
            dict(metadata, statistic="dominant_heading_turn", frame_interval=interval))
    horizon = 1600. if metadata["model"] in ("spring", "spring_ring") else 16000.
    # Keep the RMST comparison window fixed across experiments, while the
    # survival curves retain the full observed E4 lifetime range.
    survival_horizon = .5*metadata["native_horizon"] if metadata.get("native_horizon") else horizon
    if survival_horizon < horizon:
        raise ValueError("The fixed RMST window exceeds the observed lifetime window")
    for clock in ("observed", "lower", "upper"):
        grid, risk, event, censor = survival_arrays(results, clock=clock, horizon=survival_horizon,
                                                   landmarks=(horizon,))
        point = weighted_survival(risk, event, np.ones((1, n)))[0]
        boot = weighted_survival(risk, event, weights)
        ci = quantiles(boot)
        for i, time in enumerate(grid):
            survival.append(dict(metadata, lifetime_clock=clock, time=float(time),
                survival=float(point[i]), ci_lower=float(ci[0, i]), ci_upper=float(ci[1, i]),
                at_risk=int(risk[:, i].sum()), deaths=int(event[:, i].sum()), censored=int(censor[:, i].sum()),
                finite_bootstrap_fraction=float(np.isfinite(boot[:, i]).mean())))
        common = grid <= horizon
        lo, hi = rmst_bounds(grid[common], point[common], horizon)
        blo, bhi = rmst_bounds(grid[common], boot[:, common], horizon)
        metrics.append(dict(metadata, statistic="restricted_mean_lifetime_"+clock,
            estimate=float(lo[0]) if np.isfinite(lo[0]) and np.isclose(lo[0], hi[0]) else None,
            lower_bound=float(lo[0]), upper_bound=float(hi[0]),
            ci_lower=float(quantiles(blo)[0]), ci_upper=float(quantiles(bhi)[1]),
            integration_horizon=horizon, **dict(replication, uncertainty_unit=replication["uncertainty_unit"]+
                "; incident observed tracks; censored tail not extrapolated")))
    return metrics, distributions, survival


def write_csv(path, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(dict.fromkeys(k for row in rows for k in row))
    temporary = path.with_name(path.name+".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if isinstance(v, float) and not np.isfinite(v) else v) for k, v in row.items()})
    temporary.replace(path)


class CsvSink:
    """Stream condition groups, rather than retaining millions of rows in RAM."""
    context = ("model C_label f kind sigma experiment native_horizon threshold downsample_factor smoothed "
               "observation_cadence primary native_profile_cadence observation_proxy lifetimes_cadence_limited").split()
    fields = dict(
        metrics=("statistic estimate ci_lower ci_upper independent_seeds deterministic_controls bootstrap_resamples seeds_with_estimand "
                 "finite_bootstrap_fraction uncertainty_unit limited_two_seed_ensemble lower_bound "
                 "upper_bound integration_horizon").split(),
        distributions=("statistic value probability ci_lower ci_upper bin_lower bin_upper "
                       "observed_items observed_total frame_interval").split(),
        survival=("lifetime_clock time survival ci_lower ci_upper at_risk deaths censored "
                  "finite_bootstrap_fraction").split())

    def __init__(self, path, kind):
        self.path = Path(path)
        self.temporary = self.path.with_name(self.path.name+".tmp")
        self.stream = self.temporary.open("w", newline="")
        self.writer = csv.DictWriter(self.stream, fieldnames=self.context+self.fields[kind])
        self.writer.writeheader()
        self.rows = 0

    def extend(self, rows):
        for row in rows:
            self.writer.writerow({k: ("" if isinstance(v, float) and not np.isfinite(v) else v)
                                  for k, v in row.items()})
            self.rows += 1

    def finish(self):
        self.stream.close()
        self.temporary.replace(self.path)
        return dict(path=str(self.path.resolve()), sha256=sha(self.path), rows=self.rows)


def collect(inventory_path, workers, output, *, expected_inventory_sha256, worker_sources):
    inventory_path, workers, output = Path(inventory_path), Path(workers), Path(output)
    if sha(inventory_path) != expected_inventory_sha256:
        raise ValueError("Canonical region inventory identity failed")
    if (output/"complete.json").exists():
        raise FileExistsError("An audited collection is immutable")
    inventory = json.loads(inventory_path.read_text())
    cases = inventory["cases"]
    groups, receipts = defaultdict(list), []
    for case in cases:
        path = workers/case["id"]/"result.json"
        record = json.loads(path.read_text())
        if not record["passed"] or record["case"] != case or record["inventory_sha256"] != expected_inventory_sha256:
            raise ValueError("Case receipt differs from its canonical native input")
        if (record["source_sha256"] != worker_sources["signature_figure4_regions_campaign_20261020.py"]
                or record["tracking_source_sha256"] != worker_sources["signature_figure4_regions_20261020.py"]):
            raise ValueError("Region measurement was not produced by the declared frozen sources")
        receipts.append(dict(id=case["id"], path=str(path.resolve()), sha256=sha(path)))
        key = (case["model"], case["C_label"], case["f"], case["kind"], case.get("sigma", 1.),
               inventory.get("analysis_arm", case.get("experiment", "E0")), case.get("native_horizon"))
        groups[key].append(record)
    output.mkdir(parents=True, exist_ok=True)
    sinks = {kind: CsvSink(output/(kind+".csv"), kind) for kind in ("metrics", "distributions", "survival")}
    for group_index, (key, records) in enumerate(sorted(groups.items())):
        model, C, f, kind, sigma, experiment, native_horizon = key
        if len({(r["raw_profile_cadence"], r["common_imaging_cadence"]) for r in records}) != 1:
            raise ValueError("An A5 condition cannot mix native observation cadences")
        if kind == "noise" and any("seed" in r["case"] for r in records):
            seeds = [r["case"].get("seed") for r in records]
            if None in seeds or len(set(seeds)) != len(seeds):
                raise ValueError("Observation companions are not independent noisy seed replicates")
        variants = {(v["threshold"], v["downsample_factor"], v["smoothed"]) for r in records for v in r["sensitivities"]}
        variants.add((.5, records[0]["primary"]["downsample_factor"], True))
        for threshold, factor, smoothed in sorted(variants):
            results = []
            for record in records:
                options = record["sensitivities"]+[record["primary"]]
                selected = [v for v in options if (v["threshold"], v["downsample_factor"], v["smoothed"]) == (threshold, factor, smoothed)]
                if not selected:
                    raise ValueError("A requested sensitivity is missing from an independent seed")
                variant = selected[0]
                path = Path(variant["path"])
                if sha(path) != variant["sha256"]:
                    raise ValueError("Region track artifact identity failed")
                with gzip.open(path, "rt") as stream:
                    results.append(json.load(stream)["result"])
            cadence = records[0]["raw_profile_cadence"]*factor
            metadata = dict(model=model, C_label=C, f=f, kind=kind, sigma=sigma, experiment=experiment,
                native_horizon=native_horizon,
                threshold=threshold, downsample_factor=factor, smoothed=smoothed,
                observation_cadence=cadence, primary=threshold == .5 and smoothed
                    and cadence == records[0]["common_imaging_cadence"],
                native_profile_cadence=records[0]["raw_profile_cadence"],
                observation_proxy=True, lifetimes_cadence_limited=True)
            if metadata["kind"] != key[3]:
                raise ValueError("A sensitivity changed the canonical condition identity")
            a, b, c = aggregate(results, metadata=metadata, seed=2026102005+group_index)
            for sink_kind, rows in (("metrics", a), ("distributions", b), ("survival", c)):
                sinks[sink_kind].extend(rows)
    files = {kind: sink.finish() for kind, sink in sinks.items()}
    receipt_path = output/"native_receipts.json"
    save(receipt_path, receipts)
    save(output/"complete.json", dict(generated_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        inventory_sha256=expected_inventory_sha256, cases=len(cases), groups=len(groups),
        source_sha256=sha(__file__), outputs=files, bootstrap_resamples=2000,
        condition_identity_checked_against_inventory=True,
        group_kind_counts={case_kind: sum(key[3] == case_kind for key in groups)
                           for case_kind in sorted({key[3] for key in groups})},
        seed_receipts=dict(path=str(receipt_path.resolve()), sha256=sha(receipt_path)),
        scope="A5 observation statistics only; no migration inference or A6/H2 confirmation"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--inventory-sha256", required=True)
    parser.add_argument("--workers", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker-freeze", required=True, type=Path)
    parser.add_argument("--worker-freeze-sha256", required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Whole-campaign tracking/bootstrap aggregation requires SLURM")
    if sha(args.worker_freeze) != args.worker_freeze_sha256:
        raise ValueError("The tracking source freeze identity failed")
    collect(args.inventory, args.workers, args.output, expected_inventory_sha256=args.inventory_sha256,
            worker_sources=json.loads(args.worker_freeze.read_text())["source_sha256"])


if __name__ == "__main__":
    main()
