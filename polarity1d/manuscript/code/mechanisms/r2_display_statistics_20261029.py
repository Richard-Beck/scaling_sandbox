"""Read-only extraction of the three event cases in the R2 display rebuild.

No evolution, event reconstruction, condition selection or parameter changes.
The frozen campaign estimators and the preceding R2 whole-seed conventions
are reused. Run with the existing scaling_sandbox Python, under SLURM for
the full native-cadence record reads and FFTs.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import figure4_round3_estimators_20261022 as H
import ring_r2_campaign_20261027 as C

DATA = C.OUT / "supporting/data"
META = C.OUT / "supporting/metadata"
CASES = (("goryachev_perimeter", "goryachev", "per_perimeter"),
         ("spring_perimeter", "spring_ring", "per_perimeter"),
         ("spring_cell", "spring_ring", "per_cell"))
INPUTS = {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def source(path):
    path = Path(path)
    name = str(path.relative_to(C.ROOT))
    if name not in INPUTS:
        INPUTS[name] = dict(sha256=sha(path), bytes=path.stat().st_size)
    return path


def read(path):
    return C.read(source(path))


def load(path):
    return np.load(source(path), mmap_mode="r")


def cohort(run, model, size):
    manifest = read(run / "manifest.json")
    match = [c for c in manifest["conditions"]
             if c["model"] == model and c["size"] == size
             and c["rate_level"] == "low" and c["amplitude_level"] == "high"]
    assert len(match) == 1
    condition = match[0]
    expected = read(run / "expected_levels.json")
    levels = {int(expected[str(slot)]) for slot in condition["slots"]}
    assert len(levels) == 1
    level = levels.pop()
    for slot in condition["slots"]:
        for k in range(level + 1):
            status = read(C.segment(run, slot, k) / "status.json")
            assert status["complete"], (run, slot, k)
            if run == C.RUN:
                assert status.get("reduced_record")
    return manifest, condition, level


def record(run, job, level):
    parts = [C.segment(run, job["slot"], k) for k in range(level + 1)]
    modes = [load(p / "u_modes.npy") for p in parts]
    mode1 = np.concatenate([a[:, 1].astype(complex) for a in modes])
    if run == C.RUN:
        P = np.concatenate([load(p / "P.npy") for p in parts])
        t = np.concatenate([load(p / "heading_time.npy") for p in parts])
    else:
        P = np.concatenate([(a[:, 1].conj() / a[:, 0]).astype("complex64") for a in modes])
        t = np.concatenate([load(p / "time.npy") for p in parts])
    circumference = (np.concatenate([load(p / "circumference.npy") for p in parts])
                     if job["model"] == "spring_ring" else np.full(len(P), job["size"]))
    assert len(t) == len(P) == len(circumference)
    assert np.allclose(np.diff(t), job["frame_interval"], rtol=0, atol=1e-7)
    return mode1, P, t, circumference


def alpha_row(case_id, model, arm, groups, engine, production):
    sizes = (1., 2., 3.) if model == "spring_ring" else (20., 40., 60.)
    ss = [groups[s] for s in sizes]
    x = np.log([s[2].mean() for s in ss])
    a = -(x - x.mean()) / np.sum((x - x.mean())**2)
    point = [engine.half_rate(s[0].mean(0) / s[0][:, 0].mean(), s[1]) for s in ss]
    low, high = [], []
    for g, kind, bound in point:
        low.append(math.log(g if g else bound) if kind != "beyond_analysis_lag" else -math.inf)
        high.append(math.log(g if g else bound) if kind != "below_frame_resolution" else math.inf)
    lower = sum(ai * (l if ai > 0 else h) for ai, l, h in zip(a, low, high) if abs(ai) > 1e-14)
    upper = sum(ai * (h if ai > 0 else l) for ai, l, h in zip(a, low, high) if abs(ai) > 1e-14)
    resolved = all(k == "resolved" for _, k, _ in point)
    value = float(sum(ai * math.log(g) for ai, (g, _, _) in zip(a, point))) if resolved else None
    counts, _ = production.bootstrap_counts()
    counts = counts.astype(int)
    lookup = {tuple(c): ix for ix, c in enumerate(counts)}
    cached = []
    for raw, cadence, circles in ss:
        curves = counts @ raw
        cached.append([engine.half_rate(rr / rr[0], cadence) for rr in curves])
    # Preserve the preceding R2 bootstrap seed and draw order exactly.
    rng = np.random.default_rng(H.stable_seed(("R2alpha", arm, model, "low", "high")))
    bootstrap = []
    for draw in range(2000):
        ys, cs, rl, rh, kinds = [], [], [], [], []
        for si, (raw, cadence, circles) in enumerate(ss):
            c = rng.multinomial(4, np.full(4, .25))
            g, kind, bound = cached[si][lookup[tuple(c)]]
            ys.append(math.log(g if g else bound))
            cs.append(math.log(c @ circles / 4))
            kinds.append(kind)
            rl.append(math.log(g if g else bound) if kind != "beyond_analysis_lag" else -math.inf)
            rh.append(math.log(g if g else bound) if kind != "below_frame_resolution" else math.inf)
        xx = np.asarray(cs); xx -= xx.mean(); coef = -xx / (xx @ xx)
        dl = sum(ai * (l if ai > 0 else h) for ai, l, h in zip(coef, rl, rh) if abs(ai) > 1e-14)
        dh = sum(ai * (h if ai > 0 else l) for ai, l, h in zip(coef, rl, rh) if abs(ai) > 1e-14)
        bootstrap.append(dict(case_id=case_id, draw=draw, alpha_lower=dl, alpha_upper=dh,
                              status="resolved" if all(k == "resolved" for k in kinds) else "bounded"))
    status = ("resolved" if resolved else "lower_bound" if math.isfinite(lower) and not math.isfinite(upper)
              else "upper_bound" if math.isfinite(upper) and not math.isfinite(lower)
              else "unbounded" if not math.isfinite(lower) and not math.isfinite(upper) else "interval_bound")
    row = dict(case_id=case_id, model=model, arm=arm, rate_level="low", amplitude_level="high",
               status=status, alpha=value, alpha_lower=lower, alpha_upper=upper,
               alpha_ci_low=float(np.quantile([b["alpha_lower"] for b in bootstrap], .025, method="inverted_cdf")),
               alpha_ci_high=float(np.quantile([b["alpha_upper"] for b in bootstrap], .975, method="inverted_cdf")),
               bootstrap_resolved_fraction=sum(b["status"] == "resolved" for b in bootstrap) / 2000,
               rate_statuses=";".join(k for _, k, _ in point),
               mean_circumferences=";".join(f"{math.exp(xx):.12g}" for xx in x),
               analysis_lags=";".join(str((s[0].shape[1] - 1) * s[1]) for s in ss),
               crossing_lags=";".join(str(math.log(2) / g if g else None) for g, _, _ in point),
               seeds_per_size=4)
    return row, bootstrap


def main():
    DATA.mkdir(parents=True, exist_ok=True); META.mkdir(parents=True, exist_ok=True)
    production = C.legacy(); engine = production.engine()
    source(C.SOURCE / "frozen_sources/production.py")
    source(C.SOURCE / "frozen_sources/engine.py")
    source(Path(H.__file__))
    heading, sufficient, rates, alphas, bootstraps, selections = [], [], [], [], [], []
    raw_cache = {}
    for case_id, model, arm in CASES:
        sizes = (1., 2., 3.) if model == "spring_ring" else (20., 40., 60.)
        descriptors = {}
        for size in sizes:
            run = C.SOURCE if arm == "per_perimeter" or size == sizes[1] else C.RUN
            manifest, condition, level = cohort(run, model, size)
            job = manifest["jobs"][condition["slots"][0]]
            n = sum(len(load(C.segment(run, job["slot"], k) /
                             ("heading_time.npy" if run == C.RUN else "time.npy"))) for k in range(level + 1))
            descriptors[size] = run, manifest, condition, level, n, job
        cadence = descriptors[sizes[0]][5]["frame_interval"]
        maxlag = min(descriptors[size][4] // 2 * cadence for size in (sizes[0], sizes[-1]))
        lags = np.unique(np.rint(np.r_[np.linspace(cadence, maxlag, 161),
                                           np.geomspace(cadence, maxlag, 100)] / cadence).astype(int)) * cadence
        groups = {}
        for size in sizes:
            run, manifest, condition, level, _, job = descriptors[size]
            key = str(run), model, size
            moments, circles, raw = [], [], []
            for slot in condition["slots"]:
                jj = manifest["jobs"][slot]
                mode1, P, t, cc = record(run, jj, level)
                circles.append(float(cc.mean()))
                if key not in raw_cache:
                    raw.append(engine.acf_raw(mode1))
                selections.append(dict(case_id=case_id, model=model, arm=arm, size=size,
                                       slot=slot, replicate=jj["replicate"], seed=jj["seed"],
                                       run=str(run.relative_to(C.ROOT)), extension_level=level,
                                       record_start=float(t[0]), record_end=float(t[-1]),
                                       frame_interval=jj["frame_interval"],
                                       mode_stride=jj.get("mode_stride", 1), rho=jj["rho"],
                                       amplitude_mean=jj["amplitude_mean"], cell_rate=jj.get("cell_rate")))
                if size not in (sizes[0], sizes[-1]):
                    continue
                valid = (abs(P) >= .1) & np.isfinite(P)
                m = H.seed_moments(P, valid, cadence, lags)
                moments.append(m)
                zero_n = int(valid.sum())
                sufficient.append(dict(case_id=case_id, size=size, replicate=jj["replicate"], seed=jj["seed"],
                                       lag=0., cosine_sum=zero_n, eligible_pairs=zero_n, possible_pairs=len(P)))
                for i, lag in enumerate(lags):
                    sufficient.append(dict(case_id=case_id, size=size, replicate=jj["replicate"], seed=jj["seed"],
                                           lag=lag, cosine_sum=m[0, i], eligible_pairs=int(m[1, i]),
                                           possible_pairs=int(m[2, i])))
            if key not in raw_cache:
                mode_cadence = cadence * job.get("mode_stride", 1)
                raw_cache[key] = np.asarray(raw), mode_cadence, np.asarray(circles)
            groups[size] = raw_cache[key]
            rr, mode_cadence, circle_array = groups[size]
            g, kind, bound = engine.half_rate(rr.mean(0) / rr[:, 0].mean(), mode_cadence)
            lo, hi, los, his, probability = production.bootstrap_rates(rr, mode_cadence)
            rates.append(dict(case_id=case_id, model=model, arm=arm, size=size,
                              mean_physical_C=float(circle_array.mean()), gamma=g, status=kind, bound=bound,
                              gamma_ci_low=lo, gamma_ci_high=hi, ci_low_status=los, ci_high_status=his,
                              bootstrap_resolved_probability=probability, analysis_lag=(rr.shape[1]-1)*mode_cadence,
                              crossing_half_lag=math.log(2)/g if g else None,
                              mode_cadence=mode_cadence, extension_level=level, seeds=4))
            if moments:
                summary = H.summarize(moments, circles, ("R2display", case_id, size), lags, lags)
                zerorows = [r for r in sufficient if r["case_id"] == case_id and r["size"] == size and r["lag"] == 0]
                zero_n = sum(r["eligible_pairs"] for r in zerorows)
                zero_possible = sum(r["possible_pairs"] for r in zerorows)
                heading.append(dict(case_id=case_id, model=model, arm=arm, size=size, lag=0., correlation=1.,
                                    ci_low=1., ci_high=1., coverage=zero_n/zero_possible,
                                    eligible_pairs=zero_n, possible_pairs=zero_possible, seeds=4))
                for i, lag in enumerate(lags):
                    low, high = H.percentile(summary["Cboot"][:, i])
                    heading.append(dict(case_id=case_id, model=model, arm=arm, size=size, lag=float(lag),
                                        correlation=summary["C"][i], ci_low=low, ci_high=high,
                                        coverage=summary["coverage"][i], eligible_pairs=int(summary["eligible_pairs"][i]),
                                        possible_pairs=int(summary["possible_pairs"][i]), seeds=4))
        a, bs = alpha_row(case_id, model, arm, groups, engine, production)
        alphas.append(a); bootstraps.extend(bs)
        print(case_id, a["status"], a["alpha"], a["alpha_ci_low"], a["alpha_ci_high"], flush=True)
    for name, rows in (("model_heading.csv", heading), ("model_heading_seed_moments.csv.gz", sufficient),
                       ("model_alpha.csv", alphas), ("model_alpha_bootstrap.csv.gz", bootstraps),
                       ("model_rates.csv", rates), ("model_records.csv", selections)):
        pd.DataFrame(rows).to_csv(DATA/name, index=False)
    C.write(META / "display_statistics.json", dict(
        generated_utc=C.now(), extraction_code=str(Path(__file__).relative_to(C.ROOT)),
        extraction_code_sha256=sha(__file__), slurm_job_id=os.environ.get("SLURM_JOB_ID"),
        conditions="Low rate, high amplitude; exactly the requested three model/supply cases",
        source_records_read_only=True, inputs=INPUTS,
        heading="Original H.seed_moments/summarize; mean cos heading difference; endpoint |P|>=0.1; four seeds; 2000 whole-seed draws; zero-lag correlation identically one",
        heading_lags="Union of161 linear and100 logarithmic lags rounded to native cadence; common full half-record horizon for the outer sizes in each case",
        alpha="Unchanged frozen mode1 raw-pair ACF and first-half-height rate; negative three-size OLS log-rate/log-mean-physical-C slope; unchanged2000 whole-seed bootstrap draw seed/order",
        shared_spring_center="Per-cell arm uses original size2 per-perimeter trajectory, with its state-dependent hazard; outer per-cell trajectories use the fixed size2 basal-circumference rate",
        output_sha256={p.name: sha(p) for p in DATA.glob("model_*.csv*")}))


if __name__ == "__main__":
    main()
