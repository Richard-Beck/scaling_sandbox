"""Independent read-only analysis of the prospective noisy-ring experiment.

The producer is separately frozen. This module neither integrates a model nor
selects stochastic inputs. Correlations retain weak vectors and adaptation;
uncertainty resamples complete independent realizations, never time origins.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

GATES = (.05, .1, .2)
MODELS = ("wave_pinning", "goryachev", "holmes_model3")
MODEL_LABELS = dict(zip(MODELS, ("Wave-pinning", "Goryachev", "Holmes 3")))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    return path


def periodic_moments(primary, x, circumference):
    primary = np.asarray(primary, float)
    x = np.asarray(x, float)
    dx = circumference / len(x)
    assert primary.shape[-1] == len(x)
    np.testing.assert_allclose(np.diff(x), dx, atol=2e-12, rtol=2e-12)
    theta = 2 * np.pi * x / circumference
    return (dx * primary.sum(axis=-1),
            dx * (primary @ np.exp(1j * theta)),
            dx * (primary @ np.exp(2j * theta)))


def positive_lag_products(a, b=None, maxlag=None):
    """sum a[t+lag] conj(b[t]), with no circular wrap in time."""
    a = np.asarray(a)
    b = a if b is None else np.asarray(b)
    assert a.ndim == b.ndim == 1 and len(a) == len(b)
    n = len(a)
    maxlag = n - 1 if maxlag is None else int(maxlag)
    assert 0 <= maxlag < n
    size = 1 << (2 * n - 1).bit_length()
    return np.fft.ifft(np.fft.fft(a, size) * np.conj(np.fft.fft(b, size)))[:maxlag + 1]


def correlation_parts(vector, maxlag):
    """Same estimator, with positive tail sums and resolved endpoint products.

    Prefix subtraction destroys tiny later energies after a decaying transient.
    FFT absolute error also dominates tiny tail products. Use direct products
    when its conservative FFT roundoff estimate cannot resolve normalized C at
    1e-8, or the FFT product itself is near that error estimate. No C is clipped.
    """
    vector=np.asarray(vector,complex)
    count=len(vector)-np.arange(maxlag+1)
    power=np.abs(vector)**2
    prefix=np.cumsum(power)
    suffix=np.cumsum(power[::-1])[::-1]
    origin_sum=prefix[len(vector)-np.arange(maxlag+1)-1]
    later_sum=suffix[np.arange(maxlag+1)]
    fft_sum=positive_lag_products(vector,maxlag=maxlag).real
    product=fft_sum.copy()
    scale=float(np.sum(power))
    fft_error=64*np.finfo(float).eps*max(1.,np.log2(max(2,2*len(vector))))*scale
    energy_scale=np.sqrt(origin_sum)*np.sqrt(later_sum)
    fallback=(fft_error>1e-8*energy_scale)|(np.abs(fft_sum)<1000*fft_error)
    differences=np.zeros(maxlag+1)
    normalized_error=np.full(maxlag+1,np.nan)
    for lag in range(maxlag+1):
        if fallback[lag]:
            earlier=vector[:len(vector)-lag];later=vector[lag:]
            direct=float(np.vdot(earlier,later).real)
            differences[lag]=abs(direct-fft_sum[lag]);product[lag]=direct
            dot_error=(len(earlier)+8)*np.finfo(float).eps*float(np.sum(np.abs(earlier)*np.abs(later)))
            if energy_scale[lag]>0:normalized_error[lag]=dot_error/energy_scale[lag]
        elif energy_scale[lag]>0:
            normalized_error[lag]=fft_error/energy_scale[lag]
    subtraction_sum=(np.r_[0.,prefix][-1]-np.r_[0.,prefix][np.arange(maxlag+1)])
    audit=dict(direct_fallback=fallback,fft_numerator_difference=differences/count,
               tail_energy_subtraction_difference=np.abs(later_sum-subtraction_sum)/count,
               normalized_numerator_error_estimate=normalized_error,
               fft_error_estimate=fft_error/count)
    return product/count,origin_sum/count,later_sum/count,count,audit


def ratio_correlation(numerator, energy_origin, energy_later):
    denom = np.sqrt(energy_origin) * np.sqrt(energy_later)
    result = np.full(np.broadcast_shapes(numerator.shape, denom.shape), np.nan)
    np.divide(numerator, denom, out=result, where=denom > 0)
    finite = np.isfinite(result)
    assert not finite.any() or np.max(np.abs(result[finite])) <= 1 + 3e-12
    return result


def direction_parts(vector, maxlag, gate):
    magnitude = np.abs(vector)
    valid = magnitude >= gate
    unit = np.zeros(len(vector), complex)
    np.divide(vector, magnitude, out=unit, where=valid)
    numerator = positive_lag_products(unit, maxlag=maxlag).real
    eligible = positive_lag_products(valid.astype(float), maxlag=maxlag).real
    eligible = np.rint(eligible).astype(int)
    assert np.all(eligible >= 0)
    count = len(vector) - np.arange(maxlag + 1)
    return numerator, eligible, count


def percentile_bounds(values):
    """Undefined replication intervals remain NaN, with no all-NaN warnings."""
    values=np.asarray(values)
    lower=np.full(values.shape[1:],np.nan);upper=lower.copy()
    if len(values)<2:
        return lower,upper
    usable=np.isfinite(values).any(axis=0)
    if usable.any():
        lower[usable],upper[usable]=np.nanpercentile(values[:,usable],(2.5,97.5),axis=0)
    return lower,upper


def bootstrap_weights(n, draws, rng):
    if n < 2 or draws == 0:
        return np.empty((0, n))
    counts = np.zeros((draws, n), float)
    samples = rng.integers(0, n, (draws, n))
    for i in range(draws):
        counts[i] = np.bincount(samples[i], minlength=n) / n
    return counts


def first_crossing(lags, values, level=1 / math.e):
    below = np.flatnonzero(values <= level)
    if not len(below):
        return np.nan
    j = int(below[0])
    if j == 0 or not np.isfinite(values[j - 1]):
        return float(lags[j])
    # A descriptive piecewise-linear crossing on the observed curve, no fit.
    return float(np.interp(level, values[j:j + 1].tolist() + [values[j - 1]],
                           lags[j:j + 1].tolist() + [lags[j - 1]]))


def connected_regions(mask):
    mask = np.asarray(mask, bool)
    if not mask.any():
        return 0
    if mask.all():
        return 1
    return int(np.count_nonzero(mask & ~np.roll(mask, 1)))


def circular_components(mask):
    mask=np.asarray(mask,bool)
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(len(mask))]
    starts=np.flatnonzero(mask & ~np.roll(mask,1))
    components=[]
    for start in starts:
        cells=[];index=int(start)
        while mask[index]:
            cells.append(index);index=(index+1)%len(mask)
        components.append(np.asarray(cells,int))
    return components


def profile_diagnostics(primary, x, circumference, times, relative_magnitude=None,
                        physical_jump_scale=1.5):
    """Independently reproduce the producer's half-contrast native diagnostics."""
    rows=[];previous_position=None;unwrapped=0.;n=len(x);dx=circumference/n
    frozen_level=float(primary[0].min()+.5*np.ptp(primary[0]))
    for j,(time,u) in enumerate(zip(times,primary)):
        contrast=float(np.ptp(u));base=float(np.min(u))
        components=circular_components(u>=base+.5*contrast) if contrast>1e-12 else []
        frozen_regions=connected_regions(u>=frozen_level)
        peaks,properties=find_peaks(np.tile(u,3),prominence=0.)
        keep=(peaks>=n)&(peaks<2*n)
        peak_count=int(np.count_nonzero(properties["prominences"][keep]>=.05))
        records=[]
        for component in components:
            excess=u[component]-base
            mass=float(excess.sum()*dx)
            first=np.sum(excess*np.exp(2j*np.pi*x[component]/circumference))
            center=float(np.angle(first)%(2*np.pi)*circumference/(2*np.pi))
            records.append((mass,center))
        if records:
            position=max(records,key=lambda item:item[0])[1]
        elif previous_position is not None:
            position=previous_position
        else:
            position=float(x[int(np.argmax(u))])
        jump=0. if previous_position is None else (position-previous_position+circumference/2)%circumference-circumference/2
        unwrapped+=jump
        row=dict(time=float(time),contrast=contrast,baseline=base,peak=float(u.max()),
                 dominant_position=position,dominant_angle=2*np.pi*position/circumference,
                 cumulative_dominant_displacement=unwrapped,sampled_dominant_step=jump,
                 sampled_large_dominant_step=bool(abs(jump)>physical_jump_scale),
                 prominent_peak_count=peak_count,regions_frozen_half_contrast=frozen_regions,
                 regions_half_contrast=len(components),
                 single_region_resolved=bool(len(components)==1 and contrast>=.05 and peak_count<=1),
                 flat_profile=bool(contrast<=1e-12))
        if relative_magnitude is not None:
            row["relative_magnitude"]=float(relative_magnitude[j])
        rows.append(row);previous_position=position
    return rows


def selfcheck():
    rng = np.random.default_rng(918207)
    z = rng.normal(size=31) + 1j * rng.normal(size=31)
    maxlag = 9
    parts = correlation_parts(z, maxlag)
    errors = []
    for lag in range(maxlag + 1):
        a, b = z[lag:], z[:len(z) - lag]
        direct = (np.mean((a * b.conj()).real), np.mean(abs(b)**2), np.mean(abs(a)**2))
        errors.extend(abs(float(parts[k][lag]) - direct[k]) for k in range(3))
    assert max(errors) < 3e-15
    same = ratio_correlation(*correlation_parts(np.ones(30, complex), 10)[:3])
    np.testing.assert_allclose(same, 1, atol=3e-15)
    phase_rate = .17
    z = np.exp(1j * phase_rate * np.arange(40))
    np.testing.assert_allclose(ratio_correlation(*correlation_parts(z, 12)[:3]),
                               np.cos(phase_rate * np.arange(13)), atol=3e-15)
    weak = z.copy()
    weak[::4] *= .01
    for gate in GATES:
        numerator, eligible, available = direction_parts(weak, 9, gate)
        for lag in range(10):
            a, b = weak[lag:], weak[:len(weak) - lag]
            mask = (abs(a) >= gate) & (abs(b) >= gate)
            assert eligible[lag] == mask.sum()
            direct = np.sum(((a[mask] / abs(a[mask])) * (b[mask] / abs(b[mask])).conj()).real)
            assert abs(numerator[lag] - direct) < 1e-13
    assert connected_regions([1, 1, 0, 0, 1]) == 1
    assert connected_regions([1, 0, 1, 0]) == 2
    assert connected_regions([1, 1, 1]) == 1
    assert connected_regions([0, 0]) == 0
    n = 1000
    C = 10.
    x = (np.arange(n) + .5) * C / n
    theta = 2 * np.pi * x / C
    primary = np.stack([3 + np.cos(theta), 3 + np.cos(theta - .7), 3 + np.cos(2 * theta)])
    m0, m1, m2 = periodic_moments(primary, x, C)
    np.testing.assert_allclose(m0, 30, atol=1e-12)
    np.testing.assert_allclose(m1[:2], C / 2 * np.exp(1j * np.array([0., .7])), atol=1e-12)
    np.testing.assert_allclose(m1[-1], 0, atol=1e-12)
    np.testing.assert_allclose(m2[-1], C / 2, atol=1e-12)
    profiles = profile_diagnostics(primary, x, C, np.arange(3))
    assert [r["prominent_peak_count"] for r in profiles] == [1, 1, 2]
    assert first_crossing(np.arange(4), np.array([1., .8, .3, .2])) > 1
    return dict(maximum_direct_correlation_error=max(errors), harmonic_checks=True,
                periodic_peak_checks=True, magnitude_gate_pair_counts=True,
                rotating_vector_exact_correlation=True,
                uncertainty_unit="independent complete realization")


def load_case(row, campaign):
    """Producer adapter; accepts canonical trace/profile keys and aliases."""
    path = Path(row.get("archive", row.get("output", row.get("path", ""))))
    if not path.is_absolute():
        path = campaign / path
    with np.load(path, allow_pickle=False) as z:
        names = set(z.files)
        def read(*keys):
            for key in keys:
                if key in names:
                    return np.asarray(z[key]).copy()
            raise KeyError(f"None of {keys} in {path}: {sorted(names)}")
        time = read("time", "scalar_time", "trace_time")
        profile_path = Path(row.get("profile_archive", "")) if row.get("profile_archive") else None
        if profile_path is not None:
            if not profile_path.is_absolute():
                profile_path = campaign / profile_path
            with np.load(profile_path, allow_pickle=False) as profiles:
                x = np.asarray(profiles["x"]).copy()
                primary = np.asarray(profiles["primary"]).copy()
                profile_time = np.asarray(profiles["time"]).copy()
        else:
            x = read("x")
            primary = read("primary", "profiles")
            profile_time = read("profile_time", "profiles_time")
        M0 = read("M0", "m0", "total_active")
        stored_P = read("P_real", "polarity_x")+1j*read("P_imag", "polarity_y") if "P_real" in names or "polarity_x" in names else None
        stored_descriptors={key:read(key) for key in ("current_half_contrast_regions","frozen_threshold_regions","peaks_prominence_ge_005","dominant_region_position","contrast") if key in names}
        if "M1" in names:
            M1 = read("M1")
        else:
            M1 = read("M1_real", "m1_real", "harmonic_real") + 1j * read("M1_imag", "m1_imag", "harmonic_imag")
        if "M2" in names:
            M2 = read("M2")
        else:
            M2 = read("M2_real", "m2_real") + 1j * read("M2_imag", "m2_imag")
    C = float(row["circumference"])
    assert len(time) == len(M0) == len(M1) == len(M2)
    assert np.all(np.isfinite(primary)) and np.all(np.isfinite(M0))
    assert np.all(np.isfinite(M1)) and np.all(np.isfinite(M2))
    assert time[0] == 0 and np.all(np.diff(time) > 0)
    np.testing.assert_allclose(np.diff(time), np.diff(time)[0], rtol=1e-11, atol=1e-9)
    assert np.all(M0 > 0) and abs(M1[0]) / M0[0] > 1e-8
    indices = np.searchsorted(time, profile_time)
    assert np.all(indices < len(time))
    np.testing.assert_array_equal(time[indices], profile_time)
    pm0, pm1, pm2 = periodic_moments(primary, x, C)
    # Saved scalar moments must independently reproduce retained profiles.
    scale = max(1., float(np.max(M0)))
    moment_error = max(float(np.max(abs(pm0 - M0[indices]))),
                       float(np.max(abs(pm1 - M1[indices]))),
                       float(np.max(abs(pm2 - M2[indices]))))
    assert moment_error <= 3e-11 * scale, (path, moment_error)
    P = M1 / M1[0]
    if stored_P is not None:
        np.testing.assert_allclose(P,stored_P,atol=3e-14,rtol=3e-14)
    if row.get("input"):
        input_path=Path(row["input"])
        if not input_path.is_absolute():
            input_path=campaign/input_path
        with np.load(input_path) as basal:
            np.testing.assert_array_equal(primary[0],basal["state"][0])
    key = dict(model=row["model"], circumference=C, condition=row["condition"],
               fraction=float(row.get("fraction", row.get("strength_fraction", 0))),
               seed=str(row.get("seed", "")), case_id=row.get("case_id", path.stem))
    profile_rows = profile_diagnostics(primary, x, C, profile_time, abs(P[indices]))
    for native,independent in (("current_half_contrast_regions","regions_half_contrast"),
                               ("frozen_threshold_regions","regions_frozen_half_contrast"),
                               ("peaks_prominence_ge_005","prominent_peak_count")):
        if native in stored_descriptors:
            np.testing.assert_array_equal(stored_descriptors[native][indices],[r[independent] for r in profile_rows])
    centroid_validation = dict(available=False)
    if "dominant_region_position" in stored_descriptors:
        native_position = stored_descriptors["dominant_region_position"]
        independent_position = np.asarray([r["dominant_position"] for r in profile_rows])
        delta=(native_position[indices]-independent_position+C/2)%C-C/2
        resolved = np.asarray([r["regions_half_contrast"] > 0 for r in profile_rows])
        profile_flat = np.asarray([r["flat_profile"] for r in profile_rows])
        assert np.array_equal(~resolved, profile_flat)
        assert np.all(np.isfinite(native_position)) and np.all(np.isfinite(independent_position))
        # A current centroid exists only for a nonflat current region. Preserve
        # the original tolerance on every such independently observed profile.
        nonflat_error = float(np.max(abs(delta[resolved]), initial=0.))
        assert nonflat_error < 3e-10, (path, nonflat_error)
        scalar_flat = stored_descriptors["contrast"] <= 1e-12
        assert np.array_equal(scalar_flat, stored_descriptors["current_half_contrast_regions"] == 0)
        assert not scalar_flat[0] and not profile_flat[0]
        # Each cadence independently carries its preceding label when flat.
        # Labels may differ across cadences: they are not current site positions.
        scalar_cache_error = float(np.max(abs(native_position[1:][scalar_flat[1:]] - native_position[:-1][scalar_flat[1:]]), initial=0.))
        profile_cache_error = float(np.max(abs(independent_position[1:][profile_flat[1:]] - independent_position[:-1][profile_flat[1:]]), initial=0.))
        assert scalar_cache_error == 0. and profile_cache_error == 0.
        scalar_site_eligible = ((stored_descriptors["current_half_contrast_regions"] == 1)
                                & (stored_descriptors["peaks_prominence_ge_005"] <= 1)
                                & (stored_descriptors["contrast"] >= .05))
        assert not np.any(scalar_site_eligible & scalar_flat)
        centroid_validation = dict(available=True, nonflat_profile_frames=int(resolved.sum()),
            flat_profile_frames=int(profile_flat.sum()), scalar_flat_frames=int(scalar_flat.sum()),
            maximum_nonflat_position_error=nonflat_error,
            strict_nonflat_tolerance=3e-10,
            maximum_flat_cached_label_difference=float(np.max(abs(delta[profile_flat]), initial=0.)),
            scalar_cache_carry_error=scalar_cache_error, profile_cache_carry_error=profile_cache_error,
            flat_site_eligible_frames=0,
            interpretation="Flat cached labels carry last observations at separate 5/40-unit cadences; no current identifiable site is assigned")
    scalar_rows=[]
    if all(key in stored_descriptors for key in ("current_half_contrast_regions","peaks_prominence_ge_005","dominant_region_position","contrast")):
        for j,t in enumerate(time):
            scalar_rows.append(dict(time=float(t),dominant_position=float(stored_descriptors["dominant_region_position"][j]),
                                    relative_magnitude=float(abs(P[j])),
                                    single_region_resolved=bool(stored_descriptors["current_half_contrast_regions"][j]==1 and
                                        stored_descriptors["peaks_prominence_ge_005"][j]<=1 and stored_descriptors["contrast"][j]>=.05)))
    return dict(**key, time=time, P=P, R=abs(P), M0=M0,scalar_rows=scalar_rows,
                total_active_ratio=M0 / M0[0], q1=abs(M1) / M0, q2=abs(M2) / M0,
                x=x, profile_time=profile_time,
                profile_rows=profile_rows, archive=path,
                profile_archive=profile_path if profile_path is not None else path,
                audit=dict(path=str(path), sha256=sha(path),
                           profile_path=str(profile_path) if profile_path is not None else str(path),
                           profile_sha256=sha(profile_path) if profile_path is not None else sha(path), frames=len(time),
                           profile_frames=len(profile_time), moment_error=moment_error,
                           primary_minimum=float(primary.min()),
                           centroid_validation=centroid_validation))


def analyze_group(cases, maxlag_time, bootstrap_draws, rng):
    n = len(cases)
    time = cases[0]["time"]
    for case in cases[1:]:
        np.testing.assert_array_equal(case["time"], time)
    horizon = float(time[-1])
    dt = float(time[1] - time[0])
    maxlag = int(round(maxlag_time / dt))
    assert abs(maxlag * dt - maxlag_time) < 1e-8 and maxlag * dt <= horizon / 4 + 1e-8
    identity = {key: cases[0][key] for key in ("model", "circumference", "condition", "fraction")}
    identity["n_seeds"] = n
    weights = bootstrap_weights(n, bootstrap_draws, rng)
    ensemble_traces = []
    for metric in ("initial_memory", "vector_y", "R", "total_active_ratio", "q1", "q2",
                   "weak_fraction_gate0p05", "weak_fraction_gate0p1", "weak_fraction_gate0p2"):
        if metric == "initial_memory":
            values = np.stack([case["P"].real for case in cases])
        elif metric == "vector_y":
            values = np.stack([case["P"].imag for case in cases])
        elif metric.startswith("weak_fraction"):
            gate = float(metric.split("gate")[1].replace("p", "."))
            values = np.stack([case["R"] < gate for case in cases]).astype(float)
        else:
            values = np.stack([case[metric] for case in cases])
        if len(weights):
            lower, upper = percentile_bounds(weights @ values)
        else:
            lower = upper = np.full(len(time), np.nan)
        for i, t in enumerate(time):
            ensemble_traces.append(dict(**identity, time=float(t), metric=metric,
                                       mean=float(values[:, i].mean()), ci_lower=float(lower[i]), ci_upper=float(upper[i])))
    correlations, directions, seed_lags, persistence = [], [], [], []
    lags = np.arange(maxlag + 1) * dt
    windows = {"full":(0., horizon), "first_half":(0., horizon / 2), "last_half":(horizon / 2, horizon)}
    for window, (start, end) in windows.items():
        mask = (time >= start) & (time <= end)
        assert mask.sum() >= maxlag + 2
        parts = [correlation_parts(case["P"][mask], maxlag) for case in cases]
        a, e0, e1 = (np.stack([part[j] for part in parts]) for j in range(3))
        corr = ratio_correlation(a.mean(axis=0), e0.mean(axis=0), e1.mean(axis=0))
        seed_corr = ratio_correlation(a, e0, e1)
        if len(weights):
            boot_corr = ratio_correlation(weights @ a, weights @ e0, weights @ e1)
            lower, upper = percentile_bounds(boot_corr)
            integral_boot = np.trapezoid(boot_corr, lags, axis=1)
            integral_lower, integral_upper = np.percentile(integral_boot, (2.5, 97.5))
        else:
            lower = upper = np.full(len(lags), np.nan)
            integral_lower = integral_upper = np.nan
        for j, lag in enumerate(lags):
            correlations.append(dict(**identity, window=window, window_start=start, window_end=end,
                                     lag=float(lag), correlation=float(corr[j]), ci_lower=float(lower[j]), ci_upper=float(upper[j]),
                                     endpoint_energy_origin=float(e0[:,j].mean()),endpoint_energy_later=float(e1[:,j].mean()),
                                     endpoint_rms_origin=float(np.sqrt(e0[:,j].mean())),endpoint_rms_later=float(np.sqrt(e1[:,j].mean())),
                                     status=("undefined_zero_endpoint_energy" if not np.isfinite(corr[j]) else
                                             "weak_endpoint_energy_below_R0p05" if min(np.sqrt(e0[:,j].mean()),np.sqrt(e1[:,j].mean()))<.05 else "available"),
                                     direct_product_fallback_seed_fraction=float(np.mean([part[4]["direct_fallback"][j] for part in parts])),
                                     maximum_seed_fft_direct_numerator_difference=float(max(part[4]["fft_numerator_difference"][j] for part in parts)),
                                     maximum_seed_tail_energy_subtraction_difference=float(max(part[4]["tail_energy_subtraction_difference"][j] for part in parts)),
                                     restricted_horizon=maxlag_time, time_origins_per_seed=int(parts[0][3][j])))
            for k, case in enumerate(cases):
                seed_lags.append(dict(**{key:case[key] for key in ("model","circumference","condition","fraction","seed","case_id")},
                                      window=window, lag=float(lag), cross_product=float(a[k,j]),
                                      energy_origin=float(e0[k,j]), energy_later=float(e1[k,j]),
                                      correlation=float(seed_corr[k,j]),
                                      status="available" if np.isfinite(seed_corr[k,j]) else "undefined_zero_endpoint_energy",
                                      time_origins=int(parts[k][3][j]),
                                      direct_product_fallback=bool(parts[k][4]["direct_fallback"][j]),
                                      fft_direct_numerator_difference=float(parts[k][4]["fft_numerator_difference"][j]),
                                      tail_energy_subtraction_difference=float(parts[k][4]["tail_energy_subtraction_difference"][j]),
                                      normalized_numerator_error_estimate=float(parts[k][4]["normalized_numerator_error_estimate"][j])))
        persistence.append(dict(**identity, window=window, window_start=start, window_end=end,
                                restricted_horizon=maxlag_time, restricted_integral=float(np.trapezoid(corr, lags)),
                                ci_lower=float(integral_lower), ci_upper=float(integral_upper),
                                descriptive_first_1_over_e_crossing=first_crossing(lags, corr),
                                correlation_at_maxlag=float(corr[-1]),
                                crossed_1_over_e=bool(np.nanmin(corr) <= 1 / math.e),
                                normalized_correlation_status=("weak_endpoint_energy_below_R0p05" if
                                    np.min(np.sqrt(np.minimum(e0.mean(axis=0),e1.mean(axis=0))))<.05 else "available"),
                                minimum_endpoint_rms=float(np.min(np.sqrt(np.minimum(e0.mean(axis=0),e1.mean(axis=0))))),
                                window_sample_weak_fraction_R0p05=float(np.mean([case["R"][mask]<.05 for case in cases])),
                                exponential_fit=False, negative_correlations_retained=True))
        for gate in GATES:
            dparts = [direction_parts(case["P"][mask], maxlag, gate) for case in cases]
            dn = np.stack([part[0] for part in dparts])
            eligible = np.stack([part[1] for part in dparts])
            available = np.stack([part[2] for part in dparts])
            avg_eligible = eligible.mean(axis=0)
            directional = np.divide(dn.mean(axis=0), avg_eligible,
                                    out=np.full(len(lags), np.nan), where=avg_eligible > 0)
            fraction = avg_eligible / available.mean(axis=0)
            if len(weights):
                be = weights @ eligible
                bd = np.divide(weights @ dn, be, out=np.full_like(be,np.nan), where=be > 0)
                dl, du = percentile_bounds(bd)
                bootstrap_defined=np.isfinite(bd).mean(axis=0)
                fl, fu = percentile_bounds((weights @ eligible) / (weights @ available))
            else:
                dl=du=fl=fu=np.full(len(lags),np.nan)
                bootstrap_defined=np.full(len(lags),np.nan)
            for j, lag in enumerate(lags):
                directions.append(dict(**identity, window=window, gate=gate, lag=float(lag),
                                       direction_correlation=float(directional[j]), ci_lower=float(dl[j]),ci_upper=float(du[j]),
                                       eligible_pair_fraction=float(fraction[j]), eligible_ci_lower=float(fl[j]),eligible_ci_upper=float(fu[j]),
                                       bootstrap_defined_fraction=float(bootstrap_defined[j]),
                                       status="available" if avg_eligible[j] > 0 else "unavailable_no_eligible_pairs"))
    return ensemble_traces, correlations, directions, seed_lags, persistence


def site_displacement_group(cases,maxlag_time,bootstrap_draws,rng,cadence="profile"):
    """Endpoint-conditioned physical displacement, never an ancestry assignment."""
    time_key="time" if cadence=="scalar" else "profile_time"
    row_key="scalar_rows" if cadence=="scalar" else "profile_rows"
    times=cases[0][time_key]
    for case in cases[1:]:
        np.testing.assert_array_equal(case[time_key],times)
    dt=float(times[1]-times[0]);horizon=float(times[-1])
    maxlag=int(np.floor(maxlag_time/dt+1e-10))
    C=cases[0]["circumference"]
    identity={k:cases[0][k] for k in ("model","circumference","condition","fraction")}
    identity["n_seeds"]=len(cases)
    identity["site_cadence_native"]=dt
    identity["site_descriptor_source"]=cadence
    weights=bootstrap_weights(len(cases),bootstrap_draws,rng)
    rows=[]
    for window,(start,end) in {"full":(0.,horizon),"first_half":(0.,horizon/2),"last_half":(horizon/2,horizon)}.items():
        mask=(times>=start)&(times<=end)
        squared=np.zeros((len(cases),maxlag+1));cosine=squared.copy();eligible=squared.copy();available=squared.copy()
        for k,case in enumerate(cases):
            positions=np.asarray([r["dominant_position"] for r in case[row_key]])[mask]
            valid=np.asarray([r["single_region_resolved"] and r["relative_magnitude"]>=.1 for r in case[row_key]])[mask]
            for lag in range(maxlag+1):
                a=positions[lag:];b=positions[:len(positions)-lag]
                good=valid[lag:]&valid[:len(valid)-lag]
                displacement=(a-b+C/2)%C-C/2
                squared[k,lag]=np.sum(displacement[good]**2)
                cosine[k,lag]=np.sum(np.cos(2*np.pi*displacement[good]/C))
                eligible[k,lag]=good.sum();available[k,lag]=len(good)
        denominator=eligible.sum(axis=0)
        for metric,values in (("shortest_physical_displacement_squared",squared),("angular_cosine",cosine)):
            result=np.divide(values.sum(axis=0),denominator,out=np.full(maxlag+1,np.nan),where=denominator>0)
            if len(weights):
                be=weights@eligible
                boot=np.divide(weights@values,be,out=np.full_like(be,np.nan),where=be>0)
                lower,upper=percentile_bounds(boot)
                bootstrap_defined=np.isfinite(boot).mean(axis=0)
            else:
                lower=upper=np.full(maxlag+1,np.nan)
                bootstrap_defined=np.full(maxlag+1,np.nan)
            for j in range(maxlag+1):
                rows.append(dict(**identity,window=window,lag=float(j*dt),metric=metric,value=float(result[j]),
                                 ci_lower=float(lower[j]),ci_upper=float(upper[j]),
                                 eligible_pair_fraction=float(denominator[j]/available[:,j].sum()),
                                 status="available" if denominator[j]>0 else "unavailable_no_eligible_pairs",
                                 bootstrap_defined_fraction=float(bootstrap_defined[j]),
                                 same_site_ancestry_assumed=False))
    return rows


def select_histories(cases):
    """One descriptive medoid per noise group plus explicit cancellation control.

    Distance uses within-group median/IQR-standardized frozen descriptors.
    Ties use case_id. No preference for a size trend or a reorientation.
    """
    selected = []
    groups = {}
    for case in cases:
        if case["condition"] != "noise":
            continue
        key = tuple(case[k] for k in ("model","circumference","fraction"))
        groups.setdefault(key, []).append(case)
    for key, group in sorted(groups.items()):
        group = sorted(group,key=lambda c:c["case_id"])
        values=[]
        for case in group:
            times=case["time"]
            profiles=case["profile_rows"]
            values.append([float(case["P"][-1].real),float(case["R"].min()),
                           float(np.trapezoid(case["R"],times)/times[-1]),
                           float(np.mean([r["regions_half_contrast"]>1 for r in profiles]))])
        values=np.asarray(values)
        median=np.median(values,axis=0)
        scale=np.subtract(*np.percentile(values,(75,25),axis=0))
        scale=np.where(scale>1e-12,scale,1.)
        distance=np.sum(abs((values-median)/scale),axis=1)
        index=int(np.argmin(distance))
        case=group[index]
        selected.append(dict(model=key[0],circumference=key[1],fraction=key[2],
                             case_id=case["case_id"],seed=case["seed"],archive=str(case["archive"]),
                             archive_sha256=case["audit"]["sha256"],profile_archive=str(case["profile_archive"]),
                             profile_archive_sha256=case["audit"]["profile_sha256"],role="within_condition_descriptor_medoid",
                             medoid_distance=float(distance[index]),selection_n=len(group)))
    if groups:
        all_noise=[c for c in cases if c["condition"]=="noise"]
        case=min(all_noise,key=lambda c:(float(c["R"].min()),c["case_id"]))
        selected.append(dict(model=case["model"],circumference=case["circumference"],fraction=case["fraction"],
                             case_id=case["case_id"],seed=case["seed"],archive=str(case["archive"]),
                             archive_sha256=case["audit"]["sha256"],profile_archive=str(case["profile_archive"]),
                             profile_archive_sha256=case["audit"]["profile_sha256"],role="lowest_sampled_resultant_diagnostic",
                             medoid_distance=np.nan,selection_n=len(all_noise)))
    return selected


def analyze(index, output, maxlag, bootstrap_draws=1000, bootstrap_seed=920731, prefix="vision_noise_", site_cadence="profile"):
    index=Path(index).resolve();output=Path(output).resolve()
    manifest_audit={}
    if index.suffix.lower() == ".json":
        manifest=json.loads(index.read_text())
        for relative,digest in manifest["source_sha256"].items():
            assert sha(index.parent/"frozen_sources"/relative)==digest,relative
        assert sha(index.parent/"configuration.json")==manifest["configuration_sha256"]
        for path,digest in manifest.get("external_sha256",{}).items():
            assert sha(path)==digest,path
        manifest_audit=dict(frozen_source_count=len(manifest["source_sha256"]),
                            external_count=len(manifest.get("external_sha256",{})),
                            configuration_hash_verified=True,worker_output_hashes_verified=0,
                            native_horizon=manifest["horizon"],statuses=[])
        jobs=manifest["jobs"]
        rows=[]
        for job in jobs:
            folder=index.parent/"workers"/f"{int(job['slot']):04d}"
            status_path=folder/"status.json"
            assert status_path.exists(), f"Missing worker status: {status_path}"
            status=json.loads(status_path.read_text())
            state=status.get("status",status.get("state",""))
            assert status["manifest_sha256"]==sha(index)
            assert status.get("producer")=="terminal" and status.get("complete") and state=="complete",(folder,state)
            assert status["job"]==job
            assert status["native_time"]==manifest["horizon"]
            assert sha(index.parent/job["input"])==job["input_sha256"]
            for relative,digest in status["output_sha256"].items():
                assert sha(folder/relative)==digest,folder/relative
                manifest_audit["worker_output_hashes_verified"]+=1
            manifest_audit["statuses"].append(dict(slot=job["slot"],status_path=str(status_path),
                                                status_sha256=sha(status_path),elapsed_seconds=status["elapsed_seconds"],
                                                minimum_concentration=status["minimum_concentration"],
                                                max_pool_relative_error=status["max_pool_relative_error"],
                                                accepted_steps=status["accepted_steps"],nfev=status["nfev"]))
            rows.append(dict(job, status=state, model=job["model"],
                             condition=job.get("condition", {"basal":"zero"}.get(job["kind"],job["kind"])),
                             fraction=job.get("mean_fraction",job.get("fraction",0)),
                             archive=str(folder/"trace.npz"),profile_archive=str(folder/"primary.npz"),
                             case_id=f"slot_{int(job['slot']):04d}"))
    else:
        with index.open(newline="") as stream:
            rows=list(csv.DictReader(stream))
    assert rows, "No campaign cases"
    invalid=[r for r in rows if r.get("status","complete")!="complete"]
    assert not invalid, f"Incomplete or failed campaign rows must be resolved before final analysis: {invalid[:3]}"
    cases=[load_case(row,index.parent) for row in rows]
    keys=[(c["model"],c["circumference"],c["condition"],c["fraction"],c["seed"]) for c in cases]
    assert len(keys)==len(set(keys)), "Duplicate realization identities"
    groups={}
    for case in cases:
        group=tuple(case[k] for k in ("model","circumference","condition","fraction"))
        groups.setdefault(group,[]).append(case)
    accumulated=[[] for _ in range(5)]
    site_rows=[]
    rng=np.random.default_rng(bootstrap_seed)
    for key,group in sorted(groups.items()):
        group=sorted(group,key=lambda c:c["case_id"])
        result=analyze_group(group,maxlag,bootstrap_draws,rng)
        site_rows.extend(site_displacement_group(group,maxlag,bootstrap_draws,rng,site_cadence))
        for store,rows_here in zip(accumulated,result):
            store.extend(rows_here)
    data=output/"supporting/data"
    names=("vision_noise_ensemble_traces.csv","vision_noise_lag_correlations.csv",
           "vision_noise_direction_correlations.csv","vision_noise_seed_lag_parts.csv","vision_noise_persistence.csv")
    files=[write_csv(data/(prefix+name.removeprefix("vision_noise_")),rows_here) for name,rows_here in zip(names,accumulated)]
    files.append(write_csv(data/(prefix+"site_displacement.csv"),site_rows))
    del accumulated,result,site_rows
    seed_summary=[];profile_rows=[];seed_traces=[]
    for case in cases:
        ident={k:case[k] for k in ("model","circumference","condition","fraction","seed","case_id")}
        profiles=case["profile_rows"];times=case["time"]
        seed_summary.append(dict(**ident,archive=str(case["archive"]),archive_sha256=case["audit"]["sha256"],
                                 final_memory=float(case["P"][-1].real),minimum_R=float(case["R"].min()),
                                 time_mean_R=float(np.trapezoid(case["R"],times)/times[-1]),final_R=float(case["R"][-1]),
                                 final_total_active_ratio=float(case["total_active_ratio"][-1]),
                                 multi_region_sample_fraction=float(np.mean([r["regions_half_contrast"]>1 for r in profiles])),
                                 large_dominant_step_count=sum(r["sampled_large_dominant_step"] for r in profiles),
                                 cumulative_dominant_displacement=profiles[-1]["cumulative_dominant_displacement"],
                                 final_dominant_angle=profiles[-1]["dominant_angle"]))
        profile_rows.extend(dict(**ident,**row) for row in profiles)
        for j,t in enumerate(times):
            seed_traces.append(dict(**ident,time=float(t),vector_x=float(case["P"][j].real),
                                    vector_y=float(case["P"][j].imag),R=float(case["R"][j]),
                                    total_active=float(case["M0"][j]),total_active_ratio=float(case["total_active_ratio"][j]),
                                    q1=float(case["q1"][j]),q2=float(case["q2"][j])))
    files.extend([write_csv(data/(prefix+"seed_summary.csv"),seed_summary),
                  write_csv(data/(prefix+"profile_diagnostics.csv"),profile_rows),
                  write_csv(data/(prefix+"seed_traces.csv"),seed_traces),
                  write_csv(data/(prefix+"representative_histories.csv"),select_histories(cases))])
    metadata=output/"supporting/metadata";metadata.mkdir(parents=True,exist_ok=True)
    audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),analysis_only=True,
               generator_path=str(Path(__file__).resolve()),generator_sha256=sha(__file__),
               index_path=str(index),index_sha256=sha(index),cases=len(cases),groups=len(groups),
               maxlag_native_time=maxlag,bootstrap_draws=bootstrap_draws,bootstrap_seed=bootstrap_seed,
               uncertainty_unit="independent complete realization; deterministic/one-seed groups have no confidence interval",
               primary_correlation="uncentered raw vector; ratio of ensemble-averaged within-realization cross-products and endpoint energies",
               signed_integral=True,exponential_fit=False,burn_in_selected=False,site_descriptor_cadence=site_cadence,
               numerical_amendment="positive suffix energies and direct product fallback preserve the original estimator under decaying weak tails; no normalized-correlation clipping",
               centroid_validation_amendment="Strict original nonflat tolerance; independently verified flat cached labels at scalar/profile cadences; no current site assigned to flat frames; statistics unchanged",
               manifest_audit=manifest_audit,selfcheck=selfcheck(),sources=[case["audit"] for case in cases],
               outputs={str(path.relative_to(output)):sha(path) for path in files})
    path=metadata/(prefix+"analysis_audit.json")
    path.write_text(json.dumps(audit,indent=2,allow_nan=False)+"\n")
    return audit


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--maxlag",type=float)
    parser.add_argument("--bootstrap",type=int,default=1000)
    parser.add_argument("--selfcheck",action="store_true")
    parser.add_argument("--prefix",default="vision_noise_")
    parser.add_argument("--site-cadence",choices=("profile","scalar"),default="profile")
    args=parser.parse_args()
    if args.selfcheck:
        print(json.dumps(selfcheck(),indent=2));return
    assert args.index and args.output and args.maxlag is not None
    audit=analyze(args.index,args.output,args.maxlag,args.bootstrap,prefix=args.prefix,site_cadence=args.site_cadence)
    print(json.dumps({k:audit[k] for k in ("cases","groups","maxlag_native_time","bootstrap_draws")},indent=2))


if __name__=="__main__":
    main()
