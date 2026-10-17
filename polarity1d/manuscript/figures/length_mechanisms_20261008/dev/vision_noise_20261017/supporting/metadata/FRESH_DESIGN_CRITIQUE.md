# Fresh design critique: noisy-ring polarity persistence

Reviewed 17 October 2026, before stochastic integration. Scope: repository and
manuscript instructions, `MANUSCRIPT_VISION.md`, the revised [design](../../DESIGN.md),
native reaction channels and retained C60 qualification metadata. No integrations,
source changes or job launches were performed during this review.

**Assessment:** the bounded 12-noise + 12-mean + six-zero pilot is scientifically
defensible. The revised mean-channel wording and prospective statistical
appendix resolve the essential conceptual/statistical concerns below. The spatial
stationarity/kernel-width wording has also been corrected. Accept this bounded
pilot subject to frozen-config/preflight checks of cadence, prominence and
undefined denominators; no new assay,
noise-strength search or model selection is needed. This pilot directly tests
activity organization under declared fluctuations
and keeps all three supported ring models. It does not establish migration
persistence, noise robustness in general, or a causal effect of size alone.
The main ensemble must remain conditional on review of actual pilot evidence.

## Essential corrections and frozen analysis requirements

1. **Uniform-input argument: corrected during review.** A spatially uniform activation *rate*
   has no first input harmonic, but its active-material transfer is proportional
   to patterned inactive substrate. Native RHS inspection gives added transfer
   $Iv$ in Wave-pinning/Goryachev and $I R_{\mathrm i}/7.5$ in Holmes. Thus uniform
   $I=\mu$ generally contributes directly to the active first moment:

   $$
   \left.\frac{\mathrm dM_1}{\mathrm dt}\right|_{\rm added}
   =\mu V_1
   \quad\text{or}\quad
   \frac{\mu R_{{\rm i},1}}{7.5}. \tag{1}
   $$

   The revised design correctly removes its original “zero directly” statement.
   The exact same-channel mean controls are essential. Their nonlinear effect cannot be
   removed by subtracting a deterministic vector trace. Holmes must retain all
   four fields, its factor-two kinetics and Rac activation 3.5; this is the
   existing manuscript regime, not a uniquely published default.

2. **Keep the precise interpolation claim.** The selected deterministic variance
   normalization preserves standard Gaussian marginals and hence ensemble input
   mean $\mu$ at every time. It does not produce stationary continuous OU
   covariance. The revised wording correctly declares phase-dependent temporal
   covariance and distinguishes stationary knots/spatial statistics from
   interpolation in time. Its explicit kernel standard deviation
   $\ell/\sqrt{2}$ produces the intended Gaussian covariance width $\ell$ before
   finite-ring periodization. Export actual covariance and check
   regenerated input mean, variance and local correlations. A shared fine OU
   path/subsampled coarse path is an appropriate representation check; it changes
   the interpolant and is distinct from tightening the ODE solver on one fixed
   forcing path.

3. **Lag estimator and windows: resolved in the prospective appendix.** Its
   real cross-product numerator, both endpoint-energy denominators, equal seed
   weighting, time-origin windows, maximum lag and signed integral quadrature are
   appropriate. Preserve the output/lag cadence in the frozen configuration and
   mark a zero-energy denominator undefined. Use the same declared windows/lags
   across sizes within a model.
   The raw uncentered vector correlation includes retained mean polarity; its
   nonzero plateau and finite-window integral are not a stationary fluctuation
   correlation time. First/last-half differences diagnose adaptation rather than
   justify choosing a favorable burn-in. Resample independent seeds for
   uncertainty, retaining each whole trajectory. One pilot seed supports no size
   inference; 12 ensemble seeds may still leave rare switches or decay unresolved.
   Flat curves require an unresolved finite-window statement, not fitted infinity.

4. **Expose geometric and forcing differences alongside persistence.** Fixed
   local covariance gives more noisy boundary at larger circumference; it does
   not equalize integrated forcing or its uniform/first angular modes. The new
   appendix appropriately exports those diagnostics without renormalizing
   realizations. A translated
   unchanged profile has smaller angular displacement for the same physical
   displacement on a larger ring. Therefore angular-vector memory alone cannot
   establish slower local drift or stronger incumbent retention. Keep physical
   displacement, profile alignment, magnitude, total activity, second harmonic
   and multisite intervals adjacent to the vector result. Frozen local-peak
   prominence/region conventions must reject numerical background ripples while
   preserving raw profiles; a change of dominant maximum alone is not replacement.

5. **Keep comparative presentation honest.** Preserve both strengths and all
   sizes/models, including weak, destroyed, nonmonotone and unavailable cases.
   Use common dimensionless memory/magnitude axes and common limits within each
   model's native-time comparison; native cue fractions and times do not calibrate
   biological doses or rank model persistence. Choose kymographs by a declared
   seedwise rule and link their scalar traces. Experimental implications concern
   time-resolved boundary activity/asymmetry, local-site displacement and
   competition; connecting those measurements to protrusion or velocity remains
   a separate empirical question.

## Optional ideas, not pilot prerequisites

A supplementary centered vector covariance could distinguish fluctuations
around a retained mean direction, if adaptation is sufficiently resolved. It
must accompany the raw memory result rather than replace it. More noise laws,
correlation scales, amplitudes, optimized preparations or a velocity law are not
needed for this bounded test.

**Anticipated user critique:** “Does the comparative result survive seeing the
mean controls, every seed, the local activity patterns and physical drift on
honest shared scales?” The design can answer that question. A polished layout or
a favorable example cannot answer it before the actual measurements exist.
