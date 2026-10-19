# Standardized-angle checkpoint

Generated: 2026-10-19T16:17:26.586384+00:00.

The user requested the same sensible centre angles before judging which geometry claim should lead Figure 3. This separate checkpoint preserves all current main and earlier prototype figures. It does not choose a final manuscript layout.

The frozen 80-condition experiment uses the exact current basals and native protocols, with 32 proved exact reuses and 48 new SLURM branches (array 21187059). All 80 complete; all 35 overlap rows remain. The primary readouts are retained angle and magnitude; no arbitrary alignment threshold or refined minimum duration is promoted.

## Inputs and overlap

Cue centres are the declared prepared-axis coordinate plus the requested fraction, never a plateau argmax. Chemical preparation axis has absolute normal 0.225°; Spring 0.9°. Measured initial physical M₁ phases and their small offsets are recorded per condition. Spring π/8 uses continuous material centre 0.065, without snapping to a cell or renormalizing its slightly reduced sampled maximum. The analytic peak remains 4. Both the sampled shape hash and actual physical cue integral are frozen.

The support–core overlap fraction divides the physical intersection between the continuous full cue-support interval and linearly interpolated basal half-contrast core by the full physical cue support. This geometric measure is distinct from cue-weighted basal excess activity divided by total basal excess. Thus Goryachev can have no half-contrast-core intersection while still receiving input on nonzero basal tails.

model | angle_label | overlap_fraction_cue_support | initial_unit_peak_cue_excess_fraction_total
--- | --- | --- | ---
wave_pinning | pi | 0.0 | 1.6236355547901732e-13
wave_pinning | pi_2 | 0.0 | 0.0001063098307381067
wave_pinning | pi_4 | 0.8800107233974135 | 0.10306995862036293
wave_pinning | pi_8 | 1.0 | 0.12268918672229838
goryachev | pi | 0.0 | 6.207441400554662e-22
goryachev | pi_2 | 0.0 | 7.276558244160445e-11
goryachev | pi_4 | 0.0 | 2.704968797518808e-05
goryachev | pi_8 | 0.0 | 0.015966450483554658
holmes_model3 | pi | 0.0 | 3.3599819351312725e-12
holmes_model3 | pi_2 | 0.0 | 0.003782576394155654
holmes_model3 | pi_4 | 1.0 | 0.0904719362048865
holmes_model3 | pi_8 | 1.0 | 0.09062144338311513
spring | pi | 0.0 | 8.885337725048679e-06
spring | pi_2 | 0.3560143098954516 | 0.04129704469016072
spring | pi_4 | 0.9999999999999999 | 0.1000413294106167
spring | pi_8 | 1.0 | 0.1034620100112941

Native angular full widths are Goryachev 27.00000°; Holmes 27.00000°; Spring 35.22193°; Wave-pinning 27.00000°. Shared centre angles are not matched footprint or dose controls.

## Readouts and limits

The physical moments use instantaneous circumference, Mₙ = C(t) mean[u exp(2πinq)]. P = M₁/M₁(0), R = |P|, qₙ = |Mₙ|/M₀. The angle is reported only for R ≥ 0.1, with raw angles also retained. Target error is circular; the fractional requested turn divides the retained angle by the requested prepared-axis offset. This exposes why a common ±15° success tolerance is especially permissive at 22.5°. No shortest-duration claim is inferred here.

Fixed-region summaries use the original basal material half-contrast interval and physical amount, subtracting its transported initial background amount before normalization. They describe regional activity, potentially changed by translation, redistribution or chemistry; they do not distinguish a unique incumbent lineage, destruction or separate challenger when the initial cue overlaps the same region. Crossing brackets are saved-sample brackets, with censoring and rebound flags.

Full native endpoint profiles and their Fourier-shifted basal counterparts remain exported. A full-field final-tail gate is separately reported as the largest species-scale-normalized excursion over the last 1000 cue-free units (≤10⁻⁵ passes); compact reused primary-only trajectories have no invented full-field qualification. Native pool conservation and exact state/phase restarts are replayed independently of worker summaries. A finite 13,000-unit withdrawal does not prove asymptotic attraction. Native clocks are not calibrated across models; the global mean-force Spring closure does not supply a local mechanical propagation time.

## Files

- `standard_angles_summary.csv`: all 80 conditions, exact archive paths/hashes, anchor/overlap diagnostics, continuous endpoints and field flags.
- `standard_angles_endpoints.csv`: compact review table.
- `standard_angles_initial_context.csv`: 16 frozen input/overlap records.
- `standard_angles_initial_profiles_and_cues.csv`: actual sampled basal/cue fields.
- `standard_angles_final_profiles.csv`: all final primary fields and rotated basal comparisons.
- `standard_angles_fixed_region_crossings.csv`: threshold brackets with censoring and rebound flags.
- `standard_angles_traces.npz`: native sampled moment, geometry and fixed-region summaries.

The separate supporting magnitude and initial-input-context sheets retain PNG/PDF/SVG exports. Existing geometry and selected near/far duration sheets remain unchanged.

All 60 non-antipodal final primary profiles contain one resolved peak and match a fractional rotation of their own exact basal primary profile to at most 2.205e-12 of its initial contrast. Their R values range from 0.999999999999995 to 1.00000000000025. These native profile checks distinguish a translated full-strength region from vector cancellation.


The physical moments and fixed-region excess are

$$
M_n(t)=C(t)\int_0^1 u(q,t)e^{2\pi i n q}\,\mathrm{d}q,
\qquad P(t)=\frac{M_1(t)}{M_1(0)},
\qquad R(t)=|P(t)|.
\tag{1}
$$

$$
E_{\Omega_0}^{\mathrm{tr}}(t)=
\frac{C(t)\int_{\Omega_0}u(q,t)\,\mathrm{d}q
      -C_{\mathrm b}u_{\mathrm{bg},0}|\Omega_0|}
     {C_{\mathrm b}\int_{\Omega_0}
       [u(q,0)-u_{\mathrm{bg},0}]\,\mathrm{d}q}.
\tag{2}
$$

Here \(\Omega_0\) is the fixed basal half-contrast material interval; the subtracted background is its transported initial physical amount. Pure affine dilution of an otherwise unchanged material profile preserves this excess.
