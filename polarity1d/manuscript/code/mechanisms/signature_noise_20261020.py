"""Stationary nonperiodic native-input noise for the gap-hold assay.

Forcing is an explicit assay convention: weak positive activation input during
inter-pulse gaps. Chemical/mechanical parameters are unchanged. Controls use
the identical constant mean or zero input; final recovery is noise free.
"""
from __future__ import annotations

import hashlib
import math
import numpy as np
from scipy.signal import lfilter


class LinearGapNoise:
    """Separable spatial/temporal OU knots with declared positive input map.

    The spatial covariance is exp(-|x-y|/ell) on the actual nonperiodic cell
    grid. Temporal knot covariance is exp(-|t-s|/tau). Gaussian interpolation
    is divided by its analytic SD to preserve the declared marginal mean and
    variance between knots. It is a declared interpolated-knot process, not a
    claim of exact continuous OU paths. No realized-field mean normalization.
    """
    event_times = ()

    def __init__(self, length, cells, end_time, seed, mean,
                 ell=1.9, tau=10., knot=2.5, sigma=1., kind="noise",
                 transform="lognormal"):
        values = (length, end_time, mean, ell, tau, knot, sigma)
        if (not all(math.isfinite(v) for v in values)
                or min(length, end_time, ell, tau, knot) <= 0
                or mean < 0 or sigma < 0 or cells < 2
                or kind not in ("noise", "mean", "zero")
                or transform not in ("lognormal", "bounded_tanh")):
            raise ValueError("Invalid linear gap-noise parameters")
        self.length, self.cells, self.end_time = float(length), int(cells), float(end_time)
        self.seed, self.mean, self.ell = int(seed), float(mean), float(ell)
        self.tau, self.knot, self.sigma, self.kind = float(tau), float(knot), float(sigma), kind
        self.transform = transform
        self.temporal_rho = float(np.exp(-knot / tau))
        self.spatial_rho = float(np.exp(-(length / cells) / ell))
        self.knots = None
        if kind == "noise":
            count = int(np.ceil(end_time / knot)) + 1
            generator = np.random.Generator(np.random.PCG64(self.seed))
            white = generator.standard_normal((count, cells))
            # Unit marginal variance at the first spatial position as well as
            # every subsequent one; a stationary spatial AR, without wrapping.
            spatial_scale = np.sqrt(-np.expm1(-2 * length / cells / ell))
            white[:, 0] /= spatial_scale
            spatial = lfilter([spatial_scale], [1., -self.spatial_rho], white, axis=1)
            temporal_scale = np.sqrt(-np.expm1(-2 * knot / tau))
            spatial[0] /= temporal_scale
            self.knots = lfilter([temporal_scale], [1., -self.temporal_rho], spatial, axis=0)
            if not np.isfinite(self.knots).all():
                raise ValueError("Nonfinite stochastic forcing realization")

    def gaussian(self, time):
        if self.knots is None:
            raise ValueError("Deterministic controls have no Gaussian realization")
        if not -1e-9 <= time <= self.end_time + 1e-9:
            raise ValueError("Requested time is outside frozen noise realization")
        time = min(max(float(time), 0.), self.end_time)
        index = min(int(time / self.knot), len(self.knots) - 2)
        phase = time / self.knot - index
        variance = (1-phase)**2 + phase**2 + 2*phase*(1-phase)*self.temporal_rho
        return ((1-phase)*self.knots[index] + phase*self.knots[index+1]) / np.sqrt(variance)

    def __call__(self, time, x, length):
        x = np.asarray(x, dtype=float)
        if len(x) != self.cells or not np.isclose(length, self.length, rtol=0, atol=1e-9):
            raise ValueError("Noise realization must use its frozen nominal cell grid")
        if self.kind != "noise":
            return np.full_like(x, self.mean if self.kind == "mean" else 0.)
        return self.map_gaussian(self.gaussian(time))

    def map_gaussian(self, values):
        values=np.asarray(values)
        if self.transform == "bounded_tanh":
            return self.mean * (1. + np.tanh(self.sigma * values))
        return self.mean * np.exp(self.sigma * values - .5*self.sigma**2)

    def audit(self):
        if self.transform == "bounded_tanh":
            nodes, weights=np.polynomial.hermite.hermgauss(80)
            variance=self.mean**2 * float(np.sum(weights*np.tanh(self.sigma*np.sqrt(2.)*nodes)**2)/np.sqrt(np.pi))
        else:
            variance=self.mean**2 * np.expm1(self.sigma**2)
        result = dict(
            kind=self.kind, seed=self.seed, generator="numpy.Generator(PCG64)",
            cells=self.cells, nominal_length=self.length, end_time=self.end_time,
            ell=self.ell, tau=self.tau, knot=self.knot, sigma=self.sigma,
            spatial_covariance="exp(-abs(x-y)/ell), nonperiodic nominal cell grid",
            temporal_covariance="Exact stationary OU at knots; variance-normalized linear interpolation between knots",
            theoretical_input_mean=0. if self.kind == "zero" else self.mean,
            theoretical_input_variance=variance if self.kind == "noise" else 0.,
            input_transform=self.transform,
            input_formula=("mean*(1+tanh(sigma*Z))" if self.transform == "bounded_tanh"
                           else "mean*exp(sigma*Z-sigma^2/2)"),
            theoretical_hard_input_range=([0.,2*self.mean] if self.transform == "bounded_tanh" else [0.,None]),
            transformed_covariance="Nonlinear image of the declared Gaussian covariance; not itself OU or exponentially correlated.",
            application="Inter-pulse gaps only; no stochastic forcing during final native recovery",
            native_channel="Additional activation input; kinetics and mechanical coefficients unchanged",
        )
        if self.knots is not None:
            result.update(
                gaussian_knots_sha256=hashlib.sha256(
                    np.asarray(self.knots, dtype="<f8", order="C").tobytes()).hexdigest(),
                realized_gaussian_mean=float(self.knots.mean()),
                realized_gaussian_variance=float(self.knots.var()),
            )
            # Stream transformed knot moments without another full realization.
            total=squares=spatial=temporal=0.;count=spatial_count=temporal_count=0
            previous=None
            for start in range(0,len(self.knots),4096):
                values=self.map_gaussian(self.knots[start:start+4096])
                total+=float(values.sum());squares+=float(np.sum(values**2));count+=values.size
                spatial+=float(np.sum(values[:,:-1]*values[:,1:]));spatial_count+=values[:,:-1].size
                temporal+=float(np.sum(values[:-1]*values[1:]));temporal_count+=values[:-1].size
                if previous is not None:
                    temporal+=float(np.sum(previous*values[0]));temporal_count+=self.cells
                previous=values[-1].copy()
            mean=total/count;measured_variance=squares/count-mean**2
            result["transformed_knot_statistics"]=dict(mean=mean,variance=measured_variance,
                adjacent_spatial_correlation=(spatial/spatial_count-mean**2)/measured_variance if measured_variance>0 else None,
                adjacent_temporal_correlation=(temporal/temporal_count-mean**2)/measured_variance if measured_variance>0 else None,
                scope="Whole latent knot field; applied-gap statistics are audited separately.")
        return result


class ShiftedInput:
    """Evaluate a frozen global-time input inside a restarted native phase."""
    event_times = ()

    def __init__(self, stimulus, offset):
        self.stimulus, self.offset = stimulus, float(offset)

    def __call__(self, time, x, length):
        return self.stimulus(time + self.offset, x, length)
