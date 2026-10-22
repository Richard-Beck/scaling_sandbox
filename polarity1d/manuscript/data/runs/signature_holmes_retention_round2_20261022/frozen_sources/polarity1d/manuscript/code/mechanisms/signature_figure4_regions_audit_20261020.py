"""Independent native-frame/qualification/budget audit of A5 imaging views.

This auditor does not import the tracking implementation. It checks native
layout, optical qualification, periodic counts, physical widths and amount
budgets. It does not claim an independent reconstruction of every latent
identity or takeover; those have separate analytic and peer-review evidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

from signature_figure4_regions_campaign_20261020 import save, sha


def frame_checks(primary, circumference, basal_contrast, frames, *, threshold=.5):
    p = np.asarray(primary, float)
    C = np.broadcast_to(np.asarray(circumference, float), (len(p),))
    floor, contrast = p.min(axis=1), np.ptp(p, axis=1)
    qualified = (contrast > 1e-12) & (contrast >= .1*basal_contrast)
    mask = (p >= floor[:, None]+threshold*contrast[:, None]) & qualified[:, None]
    count = np.sum(mask & ~np.roll(mask, 1, axis=1), axis=1)
    count[mask.all(axis=1)] = 1
    width = mask.mean(axis=1)*C
    total = (p.mean(axis=1)-floor)*C
    if len(frames) != len(p):
        raise ValueError("Archived imaging frames differ from native cadence")
    for i, frame in enumerate(frames):
        if (frame["region_count"] != count[i]
                or frame["low_contrast"] != (contrast[i] < .1*basal_contrast)
                or frame["qualified"] != bool(qualified[i])
                or not np.isclose(frame["contrast"], contrast[i], rtol=1e-12, atol=1e-12)
                or not np.isclose(frame["circumference"], C[i], rtol=1e-12, atol=1e-12)):
            raise ValueError("Native frame count/contrast/geometry differs from A5")
        regions = frame["regions"]
        if not np.isclose(sum(r["width"] for r in regions), width[i], rtol=1e-12, atol=1e-12):
            raise ValueError("Physical region widths do not cover the native contour")
        amount = sum(r["excess"] for r in regions)
        if qualified[i] and frame["latent_peak_basins"] == 0:
            if not np.isclose(amount, total[i], rtol=1e-10, atol=1e-10):
                raise ValueError("Complete region amounts discard native activity")
        if amount < -1e-10 or amount > total[i]+1e-10*max(1., total[i]):
            raise ValueError("Region amounts exceed the native complete-excess budget")
        if any(not r["amount_available"] for r in regions) and frame["dominant_track_id"] is not None:
            raise ValueError("An ambiguous amount silently supplied a dominant identity")
    return dict(frames=len(frames), qualified_frames=int(qualified.sum()),
                native_max_region_count=int(count.max()), complete_amount_budgets_checked=True)


def takeover_checks(result):
    times = [f["time"] for f in result["frames"]]
    by_time = {f["time"]: f for f in result["frames"]}
    for event in result["takeovers"]:
        old, new = event["previous_dominant_excess"], event["newborn_excess"]
        if not old > 0 or not new > 0 or not (new-old)/(new+old) > .5:
            raise ValueError("An observed takeover lacks two positive measured amounts")
        first, last = event["first_qualifying_frame_time"], event["time"]
        identity = event["newborn_track_id"]
        if (times.index(last)-times.index(first) != 1
                or any(by_time[t]["dominant_track_id"] != identity for t in (first, last))
                or event.get("globally_dominant_required") is not True):
            raise ValueError("Pairwise dominance was mistaken for a globally dominant observed newborn")


def audit_case(case, workers, inventory_sha256):
    measurement_path = Path(workers)/case["id"]/"result.json"
    record = json.loads(measurement_path.read_text())
    if not record["passed"] or record["case"] != case or record["inventory_sha256"] != inventory_sha256:
        raise ValueError("Measurement receipt disagrees with canonical input")
    expected_views = {(threshold, factor, optical) for threshold in (.3, .5, .7)
                      for factor in (1, 2, 4, 8) for optical in (False, True)}
    actual_views = [(v["threshold"], v["downsample_factor"], v["smoothed"])
                    for v in record["sensitivities"]]
    common_cadence = 8. if case["model"] == "spring_ring" else 80.
    if (len(actual_views) != 24 or set(actual_views) != expected_views
            or record["primary"]["threshold"] != .5 or not record["primary"]["smoothed"]
            or record["primary"]["downsample_factor"]*case["profile_cadence"] != common_cadence):
        raise ValueError("Incomplete sensitivity grid or incompatible primary imaging cadence")
    for key in ("profile", "scalar", "basal_profile", "basal_checkpoint", "source_status"):
        if record["native_input_sha256"][key] != case[key+"_sha256"]:
            raise ValueError("Measurement input hash differs from canonical native inventory")
    with np.load(case["basal_checkpoint_path"], allow_pickle=False) as z:
        state = z["state"]
        native_basal = state[:case["native_cells"]] if case["model"] == "spring_ring" else state[0]
    with np.load(case["basal_profile_path"], allow_pickle=False) as z:
        basal = z["primary"].copy(); C0 = float(z["circumference"])
    if not np.array_equal(native_basal, basal):
        raise ValueError("Compact basal is not the native active-state block")
    with np.load(case["profile_path"], allow_pickle=False) as z:
        time, p = z["time"].copy(), z["primary"].copy()
        C = np.broadcast_to(z["circumference"], time.shape).copy()
    selected = np.flatnonzero(time >= case["native_horizon"]/2-1e-10)
    factor = record["primary"]["downsample_factor"]
    ix = selected[np.unique(np.r_[np.arange(0, len(selected), factor), len(selected)-1])]
    kernel = .15 if case["model"] == "spring_ring" else 1.5
    optical_basal = float(np.ptp(gaussian_filter1d(basal, kernel*len(basal)/C0, mode="wrap")))
    if not np.isclose(record["primary"]["observed_basal_contrast"], optical_basal, rtol=1e-12, atol=1e-12):
        raise ValueError("Optical basal qualification differs from exact native reference")
    optical = np.array([gaussian_filter1d(row, kernel*len(row)/length, mode="wrap")
                        for row, length in zip(p[ix], C[ix])])
    primary_path = Path(record["primary"]["path"])
    if sha(primary_path) != record["primary"]["sha256"]:
        raise ValueError("Primary A5 imaging artifact identity failed")
    with gzip.open(primary_path, "rt") as stream:
        result = json.load(stream)["result"]
    if not np.array_equal(time[ix], [f["time"] for f in result["frames"]]):
        raise ValueError("Observed frame times differ from native profiles")
    checks = frame_checks(optical, C[ix], optical_basal, result["frames"])
    dt = np.diff(time[ix]); counts = np.array([f["region_count"] for f in result["frames"]])
    expected_mean = float(np.sum(counts[:-1]*dt)/dt.sum())
    if not np.isclose(result["summary"]["region_count_mean"], expected_mean, rtol=1e-12, atol=1e-12):
        raise ValueError("Region-count summary differs from native physical-time weighting")
    takeover_checks(result)
    views = {False: p[selected], True: np.array([
        gaussian_filter1d(row, kernel*len(row)/length, mode="wrap")
        for row, length in zip(p[selected], C[selected])])}
    primary_key = (.5, factor, True)
    views_checked, all_frames = 1, len(result["frames"])
    seen = {primary_key}
    for variant in record["sensitivities"]:
        key = (variant["threshold"], variant["downsample_factor"], variant["smoothed"])
        if key in seen:
            continue
        seen.add(key)
        variant_path = Path(variant["path"])
        if sha(variant_path) != variant["sha256"]:
            raise ValueError("Sensitivity tracking artifact identity failed")
        with gzip.open(variant_path, "rt") as stream:
            alternative = json.load(stream)["result"]
        at = np.unique(np.r_[np.arange(0, len(selected), key[1]), len(selected)-1])
        if not np.array_equal(time[selected[at]], [f["time"] for f in alternative["frames"]]):
            raise ValueError("Sensitivity observation times differ from native profiles")
        contrast_reference = optical_basal if key[2] else float(np.ptp(basal))
        if not np.isclose(variant["observed_basal_contrast"], contrast_reference, rtol=1e-12, atol=1e-12):
            raise ValueError("Sensitivity optical/raw basal reference differs from native")
        frame_checks(views[key[2]][at], C[selected[at]], contrast_reference,
                     alternative["frames"], threshold=key[0])
        takeover_checks(alternative)
        views_checked += 1
        all_frames += len(alternative["frames"])
    return dict(case_id=case["id"], passed=True, measurement_sha256=sha(measurement_path),
                primary_sha256=sha(primary_path), views_checked=views_checked,
                all_view_native_frames=all_frames, **checks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--inventory-sha256", required=True)
    parser.add_argument("--workers", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--shards", type=int, default=12)
    parser.add_argument("--shard", type=int, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Native campaign audit requires SLURM")
    if sha(args.inventory) != args.inventory_sha256:
        raise ValueError("Native inventory identity failed")
    cases = json.loads(args.inventory.read_text())["cases"]
    records = [audit_case(case, args.workers, args.inventory_sha256)
               for i, case in enumerate(cases) if i % args.shards == args.shard]
    save(args.output/f"shard_{args.shard:03d}.json", dict(passed=True,
        generated_utc=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        inventory_sha256=args.inventory_sha256, shard=args.shard, shards=args.shards,
        cases=len(records), records=records,
        scope="Native block layout, exact optical/raw baseline, every declared contour/cadence/optics frame qualification/count/geometry/complete-excess budget, strict observed-global event gates, and primary summary weighting. Latent identities and takeover completeness rely on separate analytic/peer-review evidence."))


if __name__ == "__main__":
    main()
