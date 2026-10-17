# Shared KD mechanisms: PACC/non-PACC parameter comparison

Review figure: [hosny_shared_parameter_comparison.png](results/hosny_shared_parameter_comparison.png).
This presentation compares the four **baseline** parameters, with shared HIF/RhoA
mechanism subsets across PACC and non-PACC populations. It does not force equal
parameter values, equal proportional effects, or equal absolute effects.

## Architecture and representative model

Matching population-specific fits from the exhaustive sweep are paired by HIF and
RhoA subsets. Their arithmetic mean loss is exactly the ten-source-data joint
objective, because each population contributes five equally weighted source
summaries. Continuous parameters remain separate between populations; therefore
this subset-only constraint requires no new unrestricted continuous optimization.
All 225 nonempty subset pairs and the joint no-change reference are retained.

The representative is HIF-kappa / RhoA-kappa: the best shared model at the minimum
nonzero tested complexity. This is an explicit parsimony/presentation choice, not
a formal model selection rule or claim of full adequacy. Both populations share
the same allowed KD mechanisms but have separate baseline and KD parameter values.

Complexity k counts extra KD parameters **per population**. Across both populations
there are eight baseline parameters plus 2k KD parameters. k=2 therefore means
12 fitted scalar values in total, not ten. Shared subsets are discrete structural
restrictions, not shared scalar effect sizes.

## Reading the figure

A/B show the combined frontier and its mechanism subsets. C compares the four
representative baseline parameter values directly. D shows their PACC/non-PACC
ratios for every nonzero shared-frontier model without averaging parameters.
E checks each ratio *within the representative model*: fix that baseline ratio,
then jointly refit all remaining baseline and KD parameters in both populations.
All four E panels use the same error scale. These are descriptive loss profiles,
not confidence intervals. No arbitrary acceptance threshold is used.

For a fixed memory, rate or concentration ratio, the corresponding baseline
parameters are linked during optimization. A fixed kick-amplitude ratio is handled
by an exact joint quadratic minimization against both populations' path targets.
The amplitude remains shared across baseline and KDs within each population.
The concentrated noise parameters of both KDs are reoptimized in every profile.

The ratio grid spans 0.03–30 and explicitly includes 1. Two differential-evolution
seeds and a local start are evaluated per point using the existing simulation bank.
No new biological data, heterogeneity, or observation-model freedom was introduced.

The profiles are conditional on the representative mechanism and supplied-data
assumptions. Cross-frontier agreement alone does not establish identification:
parameters can compensate within each model. The unresolved directional-index
normalization, approximate quantiles and lack of replicate-level observations
remain limitations. No likelihood or parameter probability is inferred.

## Validation and reproduction

    sbatch wgd_langevin/hosny/poisson_velocity/shared_mechanisms/run.slurm

Uses the existing scaling_sandbox environment, four CPUs, 6 GB, red/small QoS.
Profiles are checkpointed in results/profiles.json. Do not reuse them after
changing model/source/loss assumptions. The figure is the only image export (PNG).

The constrained-amplitude objective is checked against direct evaluation of the
observation loss, and representative scores against the established fitter.
Every profile enforces its requested ratio to numerical tolerance. Independent
simulations with two separate seeds (32,768 cells/state each) verify the four
equal-ratio predictions against the representative model; reproduce using
`validate_equal.py`. They check numerical predictions, not new experimental data.

Outputs include all shared-model scores, frontier parameter values/ratios, complete
profiled parameter sets, representative predictions and equality validation.
Source/input hashes and modeling definitions are in metadata.json. The original
population-specific exhaustive architecture and its figure are retained separately.
