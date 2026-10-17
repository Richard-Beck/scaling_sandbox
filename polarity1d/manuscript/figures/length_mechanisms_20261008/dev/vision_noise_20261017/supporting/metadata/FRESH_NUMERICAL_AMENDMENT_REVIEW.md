# Fresh review: correlation numerical amendment

Reviewed 17 October 2026 at 20:36 UTC. This is a bounded, independent preflight
of the analysis-only amendment, **not a review of the final replicated result**.
The final report, confidence intervals, all actual plots and selected histories
remain to be assessed after complete production and analysis.

**Disposition: no essential blocker to using the separately versioned corrected
implementation.** The change repairs evaluation of the declared estimator. It
does not change its scientific definition, independent-seed uncertainty unit,
windows, lag grid, magnitude gates, representative selection or simulation.

## Verified failure and correction

The actual completed Goryachev C20 uniform-mean control at $f=0.01$ (worker
0223) reproduces the failure. Its final relative moment is
$R=3.24\times10^{-17}$. Across the full-window 1601-lag grid, the frozen
implementation loses 1110 later endpoint energies to subtractive cancellation;
these give nonfinite normalized correlations. Its largest finite absolute
correlation is 2.453, and its own bound assertion fails. These are numerical
failures, rather than unusual biological persistence.

I inspected the complete difference from the frozen source. Positive suffix
sums replace subtraction of nearly equal prefix energies; separate endpoint
square roots avoid premature underflow of their product. Direct endpoint dot
products replace the FFT numerator when its conservative roundoff estimate
cannot resolve the normalized result, or the numerator is near that estimate.
The remaining changes export diagnostics and identify this amendment. There
is no correlation clipping, mean subtraction, changed weighting, discarded
weak vector or substitute trajectory.

The preflight independently sums real-coordinate endpoint products and positive
energies at **all 1601 lags in all three windows of all 27 deterministic
controls**. The maximum actual correlation difference is
$1.28009\times10^{-12}$. The corrected full-window extreme control uses 1279
direct-product fallbacks, agrees within $7.90\times10^{-13}$ and retains finite
positive endpoint energies at every lag. These tiny energies remain subject
to the declared weak-resultant interpretation; numerical availability does not
make a vanished activity axis meaningful.

The five supplied extreme fixtures pass independent recomputation. The maximum
fixture difference is $1.20049\times10^{-12}$. The underflowed tail remains
undefined at 1600 nonzero lags. An additional exactly-zero vector remains
undefined at all 1601 lags, and a direction gate with no qualifying endpoint
pairs retains zero eligible pairs rather than creating a heading. The analytic
rotating-vector, periodic moment/peak and magnitude-gate self-checks pass.
Undefined energy statuses and unavailable direction statuses remain explicit
in the source's table-writing paths; the amendment does not change them.

## Provenance and scope

Independent hashes match the regression record and producer manifest:

- Original frozen analysis: `607de8cde2676aa0fff731ce4e9adc87750396594a3049112b89196cf1ad9a19`.
- Separate corrected source: `304c00f1f99174463faed8a984b647975e0462dd5986f62bcfc7ce1fcb7e8350`.
- Producer manifest: `3faa28c0553cf1838f99e6e4bb74c3a081f2a20704939edb413e5a7a6cba8c81`.
- Configuration: `9e8965322eab1e0d9bdc303de2f1bf645f380e29fdb402a21e2542a93d244ba1`.

All 22 frozen native executable sources, 11 frozen protocol copies, nine copied
basals, 18 original input/provenance references and 216 output hashes from the
27 completed controls pass. All control scalar grids reach 32000 in five-unit
increments. Each control status names the unchanged manifest. The frozen
`protocol/analysis_source.py` remains byte-identical to its original plan-freeze
hash. This check neither modified original sources/archives nor launched an
integration. It created only small review artifacts, below the large-output
manifest threshold.

The [independent checks](FRESH_NUMERICAL_AMENDMENT_CHECKS.json) and
[review script](review_numerical_amendment.py) preserve the bounded evidence.
The [amendment and regression](../../../../../../data/runs/vision_noise_main_20261017/results/analysis_numerical_amendment/AMENDMENT.md)
preserve the separate corrected source and original-source hashes. The final
analysis must identify the corrected source explicitly and preserve the raw
weak-energy correlation and diagnostic columns. The final actual-result review
will still examine all outcomes, independent-seed intervals, mean controls,
direction eligibility, physical versus angular motion, native profiles and
selection rationale. The high-mean Goryachev C20 control belongs to a distinct
forced-state regime, so this numerical correction does not authorize treating
it as a persistence comparison around an intact original pole.
