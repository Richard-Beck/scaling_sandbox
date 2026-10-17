# Main polarity figure

A. Recorded Wave-pinning preparation at length 15: nominal left cosine patch,
withdrawal and verified basal endpoint. The relative-20% width is 3 here.
All three response tests use the exact nominal relative-20% preparation checkpoint
at each model/length, with no reconditioning or reflection.

B. Two finite right cosine patches have fixed physical width 3 (Spring native
0.2) or relative width 20%. Their peak stays independent of length. Goryachev's
single fixed multiplier is selected by an audited bounded pilot, shared across
both patches and all lengths; every other model retains its nominal peak.
The third challenge is the existing finite g*x countergradient, with constant
physical slope and length-dependent peak. All pulses withdraw to basal input.

C. Spearman rank correlation between nominal length (native resting L0 for
Spring) and shortest passing tested cue duration, separately for each scalar
gate. The same inclusion rule applies to every model and challenge: include
every audited scalar success, irrespective of its separate native-field tail
flag, with at least four distinct successful lengths and nonconstant durations.
Tied durations receive average ranks. No significance or effect-size cutoff is
applied; weak, zero and negative correlations receive colors on the same fixed
-1 to +1 scale. Constant durations have undefined correlation; fewer than four
successful lengths is unavailable. Neither is assigned an artificial zero.
All 30 current cells are available. In-cell n counts the included scalar
successes. Native-field flags remain separate, visible as open markers in
Supplement 1 and as pass/fail/unknown counts in data/derived/length_correlations.csv.
Censored and preparation-ineligible lengths are excluded, never assigned a
duration. These success-conditional summaries can omit difficult lengths and
summarize monotonic association rather than curve shape or a scaling law.
Duration excludes preparation and the 10,000-unit recovery. The 1% strength
gate compares scalar polarity, not full-field convergence. Coefficients are
independently checked using pairwise average ranks and a centered dot product.
The older field-stable-only OLS table is retained in corrected_length_trends.csv
for provenance and is not displayed in panel C.
Statistical definition: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html

# Supplement 1

Five model rows by exactly three response-test columns. Orientation requires all
201 final-500 basal samples to exceed 0.05, with contrast at least 0.05 and
orientation range below 0.005. Strength additionally requires polarity within
1% of reversed initial scalar magnitude. Open success markers fail the separate
0.1% native-field range/projected-drift flag. Triangles mark explicit censor
limits; ineligible preparation is separate. Each endpoint is the shortest tested
passing pulse; every sampled fail-to-pass boundary was refined within 1%.
Duration excludes preparation and 10,000-unit recovery. Exact unchanged trials
can be imported with their immutable native arrays; altered stimuli or starting
checkpoints are integrated afresh. All arrays and imported identities are audited.

## Corrected deBelly boundary method

The deBelly row replaces all centered-model results with fresh relative-20% establishment/release and reversal simulations using independent local Rac/tension endpoint commands (center_boundary_motion=False). Reflection symmetry is retained; inactive endpoints remain prescribed constraints, not stress-free migration. The unchanged four models retain their audited native trials. Old centered-model figures are deprecated and archived.
