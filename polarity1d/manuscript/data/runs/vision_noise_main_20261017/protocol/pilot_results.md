# Noisy-ring pilot: actual response regimes and numerical feasibility

**All 30 prescribed pilot cases completed.** This is one noisy realization at
each model/size/strength, plus its deterministic mean and zero-input controls;
it is not a replicated persistence comparison. Native time/input units are
uncalibrated across models. Main-ensemble launch remains a separate decision.
Generated 2026-10-17T20:27:03.328669+00:00.

## Message and decision supported

The frozen input produces a useful range without changing its parameters to
obtain a preferred size trend. Wave-pinning changes weakly over this horizon;
Holmes exhibits axis changes while retaining a strong resultant; Goryachev at
C20 and the higher mean reaches weak/multiregion states. Crucially, **the
Goryachev C20 deterministic high-mean control itself loses polarity**. That
condition belongs to a different forced-state regime; its noisy response is
not a clean comparison of directional memory around the original polarized
branch. C60 retains a substantial resultant under the same local mean.

These are actual one-seed regimes for assessing the predeclared ensemble, not
proof of a size trend. Both strengths and every control remain visible. The
local noise mean and physical correlation length/time are held fixed within
model; total input and its whole-ring angular modes vary with circumference.

## Input and readouts

The native positive activation input is

$$
I(s,t)=fA_m\exp[Z(s,t)-1/2],\qquad
f\in\{0.001,0.01\},\qquad A_m\in\{0.1,6.4,1\}.\tag{1}
$$

The three $A_m$ values belong to Wave-pinning, Goryachev and Holmes respectively.
They are native cue references, not equal biological doses. The Gaussian field
has unit ensemble pointwise variance, periodized Gaussian spatial covariance
with physical Gaussian width 1.5 and exact OU knots with correlation time 50.
Knots are separated by 5; linear interpolation is divided by its analytic
standard deviation. This is a declared interpolated stochastic law, not an
exact continuous OU path. The ensemble input mean is $fA_m$; no realization is
renormalized, demeaned or clipped. Uniform mean-input controls can contribute
to the active first moment through multiplication by the nonuniform inactive
field; a uniform input harmonic alone does not imply zero stimulus moment.

The primary state vector is $P(t)=M_1(t)/M_1(0)$, initially $(1,0)$.
Initial-axis memory $\operatorname{Re}P$ is separate from the uncentered
lag correlation. Magnitude, total activity, both harmonic ratios and native
profiles accompany it. Correlations use every 5-unit lag through 1000 and are
computed separately over full/first/last-half windows without selecting burn-in.
The signed restricted integral retains negative correlations and has no forced
exponential fit. A one-seed pilot has **no replication confidence interval**.
The [frozen statistical definitions](../vision_noise_20261017/DESIGN.md) give
all equations, segmentation rules and interpretation limits.

## Initial memory and magnitude, including controls

![Every pilot memory trace](vision_noise_pilot_initial_memory.png)

Solid lines are individual noisy histories, dashed lines their deterministic
means and dotted lines zero input. Blue/orange denote C20/C60. The common axis
keeps weak changes visibly small. Initial memory can change because magnitude
changes, direction changes, or both; it is not automatically a turning rate.

![Every pilot magnitude trace](vision_noise_pilot_magnitude.png)

The lower high-input Goryachev C20 mean curve establishes the distinct
forced-state regime. A small first moment alone does not identify homogeneous
activity: the local profiles and second harmonic remain necessary.

## Lag memory and direction are kept distinct

![Vector and direction correlations](vision_noise_pilot_correlations.png)

Solid curves are full-window vector correlations; dashed curves are
magnitude-gated direction cosines. Direction results require both endpoint
magnitudes to exceed 0.1. Their usable-pair fractions and 0.05/0.2 sensitivities
are retained in the [direction table](supporting/data/vision_noise_pilot_direction_correlations.csv).
The deterministic mean is dotted where endpoint RMS magnitudes exceed 0.05.
Its weak Goryachev C20 values are withheld from this display, explicitly flagged
and still retained in the raw table. Normalizing a nearly vanished vector must
not produce a claim of persistent heading.

Most conditions have not resolved strong decorrelation by lag 1000. The
restricted integral is an observation-window descriptor, not an extrapolated
persistence time. The [full/half-window table](supporting/data/vision_noise_pilot_persistence.csv)
exposes adaptation and weak endpoint energies; it should be read beside these
raw curves rather than used as a replacement for them.

| Model | C | f | Noisy final Re P | Noisy min R | Noisy final R | Mean final R | Noisy C(1000) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Wave-pinning | 20 | 0.001 | 0.9992 | 0.9914 | 0.9995 | 1 | 0.9994 |
| Wave-pinning | 60 | 0.001 | 1.0006 | 0.9972 | 1.0006 | 1 | 1.0000 |
| Wave-pinning | 20 | 0.01 | 1.0022 | 0.9568 | 1.0027 | 0.9999 | 0.9698 |
| Wave-pinning | 60 | 0.01 | 0.9960 | 0.9659 | 1.0029 | 1.001 | 0.9948 |
| Goryachev | 20 | 0.001 | 0.9702 | 0.8884 | 0.9703 | 0.9562 | 0.9995 |
| Goryachev | 60 | 0.001 | 0.9857 | 0.9722 | 0.9857 | 0.9853 | 1.0000 |
| Goryachev | 20 | 0.01 | 0.2384 | 0.0072 | 0.3140 | 3.007e-13 | 0.2284 |
| Goryachev | 60 | 0.01 | 0.8103 | 0.7751 | 0.8107 | 0.8361 | 0.9980 |
| Holmes 3 | 20 | 0.001 | 0.9958 | 0.9954 | 0.9990 | 1.001 | 0.9983 |
| Holmes 3 | 60 | 0.001 | 0.9982 | 0.9965 | 1.0005 | 1.001 | 0.9996 |
| Holmes 3 | 20 | 0.01 | 0.8422 | 0.8614 | 1.0048 | 1.013 | 0.8814 |
| Holmes 3 | 60 | 0.01 | 0.9658 | 0.9616 | 1.0105 | 1.013 | 0.9906 |

## Pair coverage and physical versus angular motion

![Coverage and physical displacement](vision_noise_pilot_coverage_geometry.png)

The single-region endpoint diagnostic additionally requires one current
half-contrast region, at most one peak with native prominence at least 0.05,
contrast at least 0.05 and R at least 0.1 at both endpoints. Missing eligible
pairs remain unavailable. In particular, high-input Goryachev has zero eligible
single-region pairs at lag 1000, and low-input coverage is below 6% there; this
pilot supplies no resolved Goryachev physical-drift comparison under that
criterion. The native histories remain the appropriate evidence about its
local pattern dynamics.

The high-input Wave-pinning example illustrates the geometric distinction:
C20/C60 direction cosines at lag 1000 are 0.9699/0.9949, whereas eligible physical
RMS displacements are 0.7884/0.9853 native length units (pair coverage
98.7%/90.7%). Here the larger ring has stronger angular memory while its
observed endpoint physical displacement is larger. This one-seed comparison
is not a population trend; it shows why angular memory must not be relabelled
slower physical site motion. All panels keep the same physical-displacement
axis; the small Goryachev values have very low eligibility, rather than proving
exceptionally stable same-site motion.

## Every noisy native history

![All twelve noisy histories](vision_noise_pilot_all_noise_histories.png)

All twelve stochastic histories are shown, so no pilot example is selected to
support enlargement increasing persistence. Each row uses one native activity
colour scale across both sizes and strengths; the periodic cut keeps the
original incumbent at position zero. The same angular axis spans different
physical distances at C20 and C60. Quantitative
[endpoint-conditioned physical displacement](supporting/data/vision_noise_pilot_site_displacement.csv)
therefore accompanies angular correlations. These 20-unit profile-cadence endpoints do
not prove continued site identity. A dominant-site change alone is not called
replacement. The [profile table](supporting/data/vision_noise_pilot_profile_diagnostics.csv)
retains current/frozen half-contrast region counts, prominent peaks and tracking
ambiguities.

The lowest noisy Goryachev C20 high-input resultant occurs at $t=3770$:
$R=0.00721$, while contrast is 3.281, $q_2=0.1970$ and two current half-contrast
regions with three prominent peaks remain. Its deterministic mean control ends
with contrast approximately 2.95×10⁻¹² and R approximately 3.01×10⁻¹³. Thus the
noisy low first moment is not equivalent to homogeneous activity. These exact
scalar events are retained in the corresponding worker trace archives.

The medoid candidate table is [retained](supporting/data/vision_noise_pilot_representative_histories.csv)
for the later replicated figure. With one seed per condition its within-group
medoids are simply those individual histories, not representative population
examples. A separate lowest-resultant history is labelled a diagnostic.

## Actual input statistics and numerical accounting

The [twelve input audits](supporting/data/vision_noise_pilot_input_audits.csv) preserve
realized local means/second moments, spatial and temporal covariance checks,
and the zeroth/first whole-ring input modes alongside their finite-ring
expectations. Finite-realization deviations are reported, not corrected by
rescaling. Whole-ring mode amplitudes must not be mistaken for equal forcing
at different circumferences merely because local covariance is fixed.

The independent read-only analysis verified 22
frozen executable source hashes, 12 original
input/provenance hashes and 240
worker-output hashes. All 30 initial primary profiles equal their exact copied
basals. All scalar and profile grids reach 4000 at their prescribed 5/20-unit
cadences. Every retained profile independently reproduces the raw zeroth/first/
second moments and native region/peak counts. Maximum raw-moment discrepancy
is 2.27e-13; maximum conserved-pool
relative error is 2.61e-12.
The maximum zero-control final-memory deviation is 2.22e-16.
Worker elapsed time ranges from 0.79
to 119.59 seconds.

The [analysis audit](supporting/metadata/vision_noise_pilot_analysis_audit.json) records
complete exact sources and outputs. The analysis also passes direct versus FFT
lag-product checks, an analytic rotating-vector correlation, periodic harmonic
and peak tests, and magnitude-gate pair-count tests. The [complete numerical comparison](../../../../data/runs/vision_noise_numerical_audit_20261017/NUMERICAL_COMPARISON.json)
contains six independently verified C60 high-input checks. Tightening ODE
tolerances for the identical forcing changes $P$ by at most
$4.53\times10^{-8}$, $1.40\times10^{-9}$ and $7.70\times10^{-8}$
in Wave-pinning, Goryachev and Holmes respectively. The coupled 5-to-2.5-unit
interpolation comparison changes maximum $P$ by 0.00546, 0.00176 and 0.01807;
maximum raw-profile changes are 2.31%, 0.219% and 12.74% of basal contrast, and
maximum dominant-region position differences are 0.0566, 0.0156 and 0.1866
native length units. Native profiles are therefore not numerically identical.
This comparison changes the declared interpolated driving path and measures
sensitivity to that construction; it is not additional ODE error or proof of
continuous-OU convergence.

Independent direct recomputation of full-window vector correlations finds
maximum coarse/fine curve differences 0.000110, 0.0000176 and 0.000563. Their
restricted integrals change 997.469→997.534, 998.324→998.323 and
992.855→992.469 native time units. Thus these global pilot summaries are less
sensitive than individual local profiles over this check. Prominent-peak counts
can change under the interpolation refinement, especially in Goryachev; this
reinforces the explicit coverage limits on its local-site displacement summary.
A different uncoupled random realization would not be a convergence check.

## Why these displays and what could be measured

Memory and magnitude traces distinguish maintained directional asymmetry from
loss of its resultant; mean controls expose a change of the forced state.
Vector/direction curves prevent amplitude fluctuations being relabelled as
turning. The native histories explain whether vector changes accompany a moving
single site, multiple regions or loss of local contrast. The same ingredients
could be obtained from a whole-boundary activity reporter, with physical site
positions retained alongside the angular vector. This proposes an observable
polarity-memory measurement; it supplies no cell velocity law or calibrated
migration persistence estimate.

## Reviewed pilot and bounded main continuation

The [fresh actual-result critique](supporting/metadata/FRESH_RESULT_CRITIQUE.md)
finds no remaining essential blocker after inspecting all five actual plots,
controls, raw archives and complete numerical comparisons. The parent has
reviewed this evidence and authorized the bounded main ensemble. No extra
amplitude has been selected from these outcomes, and no missing or failed
realization is silently censored.

The main design retains both declared strengths and explicitly labels the
high-mean Goryachev C20 regime. It uses 12 independent seeds per
model/size/strength at C20/40/60 over 32000 native units, with a shared maximum
lag of 8000. The longer prospective horizon addresses unresolved weak-pilot
decay without selecting a favourable size trend; it is not yet an ensemble
result. Main site diagnostics use the scalar cadence of 5, independently
checked at every saved profile snapshot; main saved profiles use cadence 40.
Pilot physical-displacement diagnostics remain at their declared cadence of 20.
This report records the completed pilot, not the outcome of that continuation.
