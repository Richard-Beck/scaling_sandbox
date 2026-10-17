# Numerical amendment to the frozen vector-correlation implementation

The main producer, all frozen native sources, configuration, input/noise law,
realizations and original prospective analysis remain unchanged. This amendment
changes numerical evaluation of the already frozen correlation formula; it does
not change the scientific estimator, time-origin windows, lag grid, bootstrap,
selection rules or weak-resultant criteria.

## Reproduced failure

The completed high-mean Goryachev C20 control spans many orders of magnitude in
first-moment amplitude. In the original implementation, later endpoint energy
was computed by subtracting a large prefix sum from the total. Positive tiny
tail energies were lost to floating-point cancellation; FFT endpoint products
also became inaccurate relative to the tiny products. Directly applying the
original correlation to its complete trace reproduces the assertion failure,
including nonfinite normalized values. Those values are an analysis defect,
not evidence of exceptional persistence.

## Separate corrected source

The original `../../protocol/analysis_source.py` is preserved byte-for-byte.
The separate [corrected source](corrected_analysis_source.py) evaluates positive
later energies by suffix accumulation, computes the two square roots before
multiplying them, and uses direct endpoint dot products when the conservative
FFT roundoff estimate cannot resolve the normalized product or the numerator
itself lies near that estimate. No normalized correlation is clipped into a
permitted interval. Zero endpoint energies remain undefined, and weak endpoint
RMS values retain their declared 0.05 status. Seed-level output records the
fallback, numerator sensitivity and tail-subtraction sensitivity.

The [regression record](REGRESSION.json) ties both source hashes to the unchanged
producer manifest and configuration. It compares every one of 1601 prescribed
lags in all three windows for all 27 completed mean/zero controls directly
against endpoint dot products and positive energy sums: 81 windows in total.
The largest actual absolute correlation discrepancy is 1.281 × 10⁻¹².

Additional fixtures cover a decaying transient plus tiny oscillation, an impulse
plus a 10⁻¹²⁰ tail, uniformly tiny and large resolved vectors, and an
underflowed-energy tail. All pass; genuinely unrepresentable endpoint energies
remain undefined. The original analytic rotating-vector, periodic moment,
peak-count and magnitude-gate checks also pass.

This record and source are separate from the producer's frozen tree. A fresh
review must inspect the correction before the corrected implementation is used
for final ensemble statistics. The reviewed pilot artifacts are unchanged.
