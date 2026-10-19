"""Bounded, read-only inventory of the existing four-model noisy-ring regimes.

No native integration, representative reselection, event classifier or artwork.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
DEV = ROOT / "polarity1d/manuscript/figures/length_mechanisms_20261008/dev"
FIG = DEV / "revision_program_20261019/figure4"
OLD = DEV / "vision_noise_20261017"
OUT = FIG / "regime_review_20261019"
DATA = OUT / "supporting/data"
META = OUT / "supporting/metadata"
MODELS = ("wave_pinning", "goryachev", "holmes_model3", "spring_ring")


def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def export(name, records):
    path = DATA / (name + ".csv")
    keys = list(dict.fromkeys(k for r in records for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(records)
    return path


def runs(mask, time):
    """Consecutive positive *sample* runs; no continuous-state duration claim."""
    edges = np.diff(np.r_[False, mask, False].astype(int))
    starts, stops = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1) - 1
    return [dict(start_index=int(a), last_index=int(b), first_positive_time=float(time[a]),
                 last_positive_time=float(time[b]), observed_span=float(time[b] - time[a]),
                 positive_samples=int(b - a + 1), left_window_censored=bool(a == 0),
                 right_window_censored=bool(b == len(time) - 1),
                 preceding_negative_time=float(time[a - 1]) if a else np.nan,
                 following_negative_time=float(time[b + 1]) if b + 1 < len(time) else np.nan)
            for a, b in zip(starts, stops)]


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)
    chemical = rows(FIG / "supporting/data/native_sources.csv")
    spring = rows(FIG / "supporting/data/spring_noise_seed_summary.csv")
    reps = rows(OLD / "supporting/data/vision_noise_representative_histories.csv")
    reps = [r for r in reps if r["role"] == "within_condition_descriptor_medoid"]
    reps += rows(FIG / "supporting/data/spring_noise_representative_histories.csv")
    repkeys = {(r["model"], int(r.get("slot", r.get("case_id", "slot_0").split("_")[-1]))) for r in reps}
    index = []
    for r in chemical:
        index.append(dict(model=r["model"], size=float(r["circumference"]), fraction=float(r["fraction"]),
                          condition=r["condition"], seed=int(r["seed"]), slot=int(r["slot"]),
                          archive=r["path"], expected_sha=r["sha256"]))
    for r in spring:
        index.append(dict(model=r["model"], size=float(r["resting_circumference"]), fraction=float(r["fraction"]),
                          condition=r["condition"], seed=int(r["seed"]), slot=int(r["slot"]),
                          archive=r["archive"], expected_sha=r["archive_sha256"]))
    assert len(index) == 324
    cases, episodes, medoids, examples = [], [], [], []
    cache = {}
    profiles_checked = 0
    for base in index:
        path = Path(base["archive"])
        status = json.loads((path.parent / "status.json").read_text())
        assert sha(path) == base["expected_sha"] == status["output_sha256"]["trace.npz"]
        with np.load(path) as z:
            t = z["time"].copy()
            reg = z["current_half_contrast_regions"].copy()
            peaks = z["peaks_prominence_ge_005"].copy()
            contrast = z["contrast"].copy()
            P = z["P_real"] + 1j * z["P_imag"]
            R, q1, M0 = abs(P), z["q1"].copy(), z["M0"].copy()
        dt = .5 if base["model"] == "spring_ring" else 5.
        assert len(t) == 6401 and np.all(np.diff(t) == dt)
        profile_path = path.parent / "primary.npz"
        profile_sha = sha(profile_path)
        assert profile_sha == status["output_sha256"]["primary.npz"]
        with np.load(profile_path) as z:
            pt, u = z["time"].copy(), z["primary"].copy()
        ix = np.searchsorted(t, pt)
        np.testing.assert_array_equal(t[ix], pt)
        pc = np.ptp(u, axis=1)
        mask = u >= u.min(axis=1)[:, None] + pc[:, None] / 2
        count = np.sum(mask & ~np.roll(mask, 1, axis=1), axis=1)
        count[mask.all(axis=1)] = 1
        count[pc <= 1e-12] = 0
        np.testing.assert_array_equal(count, reg[ix])
        profiles_checked += len(pt)
        ismedoid = (base["model"], base["slot"]) in repkeys
        common = {k: v for k, v in base.items() if k != "expected_sha"}
        common.update(scalar_sha256=base["expected_sha"], profile_archive=str(profile_path),
                      profile_sha256=profile_sha, existing_condition_medoid=ismedoid)
        masks = dict(multiregion=reg > 1,
                     resolved_multiregion=(reg > 1) & (contrast >= .05) & (peaks > 1),
                     weak_resultant=R < .1)
        bykind = {name: runs(m, t) for name, m in masks.items()}
        case = dict(**common, scalar_cadence=dt, profile_cadence=float(pt[1]-pt[0]),
                    observations=len(t), saved_profiles=len(pt), total_time=float(t[-1]),
                    multiregion_fraction=float(np.mean(reg > 1)),
                    resolved_multiregion_fraction=float(np.mean(masks["resolved_multiregion"])),
                    multiple_prominent_peaks_fraction=float(np.mean(peaks > 1)),
                    weak_R_fraction=float(np.mean(R < .1)), low_contrast_fraction=float(np.mean(contrast < .05)),
                    max_current_regions=int(reg.max()), max_prominent_peaks=int(peaks.max()),
                    mean_R=float(R.mean()), final_R=float(R[-1]),
                    mean_active_amount_ratio=float(np.mean(M0/M0[0])),
                    mean_q1_ratio=float(np.mean(q1/q1[0])), initial_q1=float(q1[0]))
        for kind, rr in bykind.items():
            spans = [r["observed_span"] for r in rr]
            case[kind+"_sample_runs"] = len(rr)
            case[kind+"_longest_observed_span"] = max(spans, default=0.)
            case[kind+"_median_observed_span"] = float(np.median(spans)) if spans else np.nan
            for j, r in enumerate(rr):
                a, b = r["start_index"], r["last_index"]
                episodes.append(dict(model=base["model"], size=base["size"], fraction=base["fraction"],
                    condition=base["condition"], seed=base["seed"], slot=base["slot"],
                    kind=kind, run_index=j, **r, max_regions=int(reg[a:b+1].max()),
                    max_prominent_peaks=int(peaks[a:b+1].max()),
                    minimum_contrast=float(contrast[a:b+1].min()), minimum_R=float(R[a:b+1].min())))
        cases.append(case)
        if ismedoid:
            medoids.append(case)
            for kind in ("multiregion", "resolved_multiregion"):
                if not bykind[kind]:
                    continue
                longest = max(bykind[kind], key=lambda r: (r["observed_span"], -r["start_index"]))
                a, b = longest["start_index"], longest["last_index"]
                profile_ids = np.flatnonzero((pt >= t[a]) & (pt <= t[b]))
                if not len(profile_ids):
                    continue
                middle = (t[a]+t[b])/2
                j = int(profile_ids[np.argmin(abs(pt[profile_ids]-middle))])
                k = int(ix[j])
                examples.append(dict(**common, kind=kind,
                    selection="Longest consecutive-positive sampled run within existing condition medoid; midpoint-nearest saved field",
                    **longest, profile_time=float(pt[j]), regions=int(reg[k]), prominent_peaks=int(peaks[k]),
                    contrast=float(contrast[k]), R=float(R[k]), q1=float(q1[k]),
                    active_amount_ratio=float(M0[k]/M0[0]), native_field_min=float(u[j].min()),
                    native_field_max=float(u[j].max()), profiles_in_observed_span=len(profile_ids)))
        cache[(base["model"], base["slot"])] = case
    # Existing inferential fractions must reproduce exactly before adding descriptors.
    oldseed = rows(FIG / "supporting/data/seed_observation_summary.csv") + spring
    maximum_existing_difference = 0.
    for r in oldseed:
        c = cache[(r["model"], int(r["slot"]))]
        for metric in ("multiregion_fraction", "weak_R_fraction", "mean_R", "final_R"):
            maximum_existing_difference = max(maximum_existing_difference, abs(float(r[metric])-c[metric]))
    assert maximum_existing_difference < 1e-12
    oldobs = rows(FIG / "supporting/data/observation_statistics.csv") + rows(FIG / "supporting/data/spring_noise_observation_statistics.csv")
    rng = np.random.default_rng(2026101991)
    groups = []
    for model in MODELS:
        sizes = (1., 2., 3.) if model == "spring_ring" else (20., 40., 60.)
        for size in sizes:
            for kind, fraction in (("noise", .001), ("noise", .01), ("mean", .001), ("mean", .01), ("basal", 0.)):
                cc = [r for r in cases if r["model"] == model and r["size"] == size and r["condition"] == kind and r["fraction"] == fraction]
                assert len(cc) == (12 if kind == "noise" else 1)
                summary = dict(model=model, size=size, fraction=fraction, condition=kind, cases=len(cc),
                    scalar_cadence=cc[0]["scalar_cadence"], profile_cadence=cc[0]["profile_cadence"],
                    observations_per_case=6401, total_time_per_case=cc[0]["total_time"],
                    seeds_with_multiregion=sum(r["multiregion_fraction"] > 0 for r in cc),
                    seeds_with_resolved_multiregion=sum(r["resolved_multiregion_fraction"] > 0 for r in cc),
                    seeds_with_multiple_peaks=sum(r["multiple_prominent_peaks_fraction"] > 0 for r in cc),
                    maximum_current_regions=max(r["max_current_regions"] for r in cc),
                    maximum_prominent_peaks=max(r["max_prominent_peaks"] for r in cc))
                weights = rng.multinomial(12, np.full(12, 1/12), 2000)/12 if kind == "noise" else None
                for metric in ("multiregion_fraction", "resolved_multiregion_fraction", "multiple_prominent_peaks_fraction", "weak_R_fraction", "low_contrast_fraction", "mean_R", "final_R"):
                    v = np.array([r[metric] for r in cc])
                    summary[metric] = float(v.mean())
                    lo, hi = np.quantile(weights @ v, [.025, .975]) if weights is not None else (np.nan, np.nan)
                    # Preserve previously accepted intervals for unchanged primary fractions.
                    if kind == "noise" and metric in ("multiregion_fraction", "weak_R_fraction", "mean_R", "final_R"):
                        original = [r for r in oldobs if r["model"] == model and float(r.get("circumference", r.get("resting_circumference"))) == size and float(r["fraction"]) == fraction and r["condition"] == "noise" and r["metric"] == metric]
                        assert len(original) == 1
                        np.testing.assert_allclose(float(original[0]["value"]), v.mean(), atol=1e-12, rtol=0)
                        lo, hi = float(original[0]["ci_lower"]), float(original[0]["ci_upper"])
                    summary[metric+"_ci_lower"], summary[metric+"_ci_upper"] = float(lo), float(hi)
                for name in ("multiregion", "resolved_multiregion", "weak_resultant"):
                    longest = [r[name+"_longest_observed_span"] for r in cc]
                    summary[name+"_longest_span_min"] = min(longest)
                    summary[name+"_longest_span_median"] = float(np.median(longest))
                    summary[name+"_longest_span_max"] = max(longest)
                    summary[name+"_total_sample_runs"] = sum(r[name+"_sample_runs"] for r in cc)
                groups.append(summary)
    paths = [export("case_inventory", cases), export("condition_inventory", groups),
             export("sampled_runs", episodes), export("existing_medoid_inventory", medoids),
             export("native_medoid_examples", examples)]
    inputs = [FIG / "supporting/data/native_sources.csv", FIG / "supporting/data/spring_noise_seed_summary.csv",
              FIG / "supporting/data/seed_observation_summary.csv", FIG / "supporting/data/observation_statistics.csv",
              FIG / "supporting/data/spring_noise_observation_statistics.csv",
              OLD / "supporting/data/vision_noise_representative_histories.csv",
              FIG / "supporting/data/spring_noise_representative_histories.csv"]
    old_art = json.loads((FIG / "supporting/metadata/figure4_with_spring_render_validation.json").read_text())["output_hashes"]
    main_art = [r for r in old_art if Path(r["path"]).stem in ("figure4_with_spring_dynamics", "figure4_with_spring_geometry_patterns")]
    assert len(main_art) == 6
    for r in main_art:
        assert sha(r["path"]) == r["sha256"]
    audit = dict(generated_at=datetime.now(timezone.utc).isoformat(), source=str(Path(__file__).relative_to(ROOT)),
                 source_sha256=sha(__file__), cases=len(cases), noisy_cases=sum(r["condition"]=="noise" for r in cases),
                 existing_medoids=len(medoids), profiles_region_replayed=profiles_checked,
                 maximum_existing_summary_difference=maximum_existing_difference,
                 bootstrap_new_descriptor_draws=2000, bootstrap_new_descriptor_seed=2026101991,
                 bootstrap_unit="whole independent seed; deterministic controls receive no bootstrap interval",
                 primary_fraction_intervals="Unchanged accepted intervals copied after point-estimate replay",
                 run_definition="Consecutive positive scalar samples, span=last sample minus first; no uninterrupted-duration or mechanistic event claim",
                 source_inputs=[dict(path=str(p.relative_to(ROOT)), sha256=sha(p)) for p in inputs],
                 current_main_artifacts_unchanged=main_art,
                 artifacts=[dict(path=str(p.relative_to(ROOT)), sha256=sha(p)) for p in paths])
    (META / "validation.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(json.dumps({k:v for k,v in audit.items() if k != "artifacts"}, indent=2))


if __name__ == "__main__":
    main()
