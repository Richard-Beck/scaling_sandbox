# Exhaustive mechanism comparison: findings

The clearest economical PACC explanation is **directionality reduction after both
HIF and RhoA KD**. At the same two-parameter budget as the original hypothesis,
its independently simulated error is 0.0273 versus 0.0512 for HIF directionality
reduction plus RhoA memory reduction. The next-best two-parameter PACC mechanism
has error 0.0422, a considerably larger gap than the measured Monte Carlo repeat
variation. This supports changing the mechanistic interpretation within this model;
it does not establish a unique cellular cause.

The non-PACC result is less specific. Its nominal best two-parameter model changes
HIF memory and RhoA kick amplitude, with error 0.0440 versus 0.0492 for the original
hypothesis. But the four combinations of HIF memory/amplitude and RhoA
memory/amplitude have errors 0.04403–0.04429. These are near-ties relative to the
approximately 0.0006 maximum observed repeat discrepancy. It would be misleading
to single out one of those four mechanisms as clearly identified.

## Best fit at each complexity

| KD-specific parameters | non-PACC error | PACC error |
|---|---:|---:|
| 0 | 0.05594 | 0.06115 |
| 2 | 0.04403 | 0.02735 |
| 3 | 0.03572 | 0.02521 |
| 4 | 0.03220 | 0.02123 |
| 5 | 0.02793 | 0.01588 |
| 6 | 0.02421 | 0.01202 |
| 7 | 0.02455 | 0.01091 |
| 8 | 0.02462 | 0.01068 |

Every model has four additional jointly fitted baseline parameters. The exact
winning subsets appear in panel B and `best_at_complexity.csv`. The independently
evaluated non-PACC frontier ends at six KD parameters; its tiny apparent advantage
over seven/eight should not be treated as a statistically resolved preference.
The optimizer's simulation-bank loss respects every subset/superset inequality;
the small reversals are in independent simulation performance. PACC improvements
become small near seven/eight parameters. No formal elbow or preferred global
complexity has been selected.

## Illustrated equal-complexity models

For PACCs, the best two-parameter fit shares tau, kick rate and kick amplitude
across baseline and both KDs. Kappa is approximately 0.512 at baseline, 0.099 after
HIF KD and 0.013 after RhoA KD. This suppresses the excess forward movement left
by the original RhoA-memory-only hypothesis. Other distributional discrepancies
remain visible in panel C.

For non-PACCs, the illustrated nominal winner reduces tau from about 0.040 to
0.029 h after HIF KD and reduces kick amplitude from about 256 to 205 um/h after
RhoA KD. Its close alternatives prevent a unique interpretation of these changes.

The short fitted memory times and large velocity-jump amplitudes are not direct
measurements of instantaneous cell speeds. Their product can preserve sampled
motion while individual values vary substantially below the 30-minute cadence.
These uncertainties persist despite enumerating all sharing structures.

## Validation

452 joint models were fitted: 225 nonempty HIF/RhoA subset pairs plus a joint
no-change null for each PACC class. All four parameters can increase or decrease
when selected; unselected parameters equal the jointly fitted baseline. The
original decrease-constrained hypothesis is an additional marked reference.

All simulator checks passed. The accelerated objective agrees with the original
objective across every sharing structure. Every subset/superset loss inequality
was checked. Every fitted model was independently evaluated using 32,768 cells
per state; all best-at-complexity models received a second independent evaluation.
Detailed original-versus-best-two comparisons use 65,536 cells/state.

Fresh simulation is a numerical validation of the fitted predictions, not an
experimental train/test split. The source contains approximate aggregate quartiles,
not raw trajectories or replicate-specific distributions. No p-values or parameter
confidence intervals are inferred from the Pareto calculation.

All best-at-complexity representatives had a converged selected optimization run.
Among other models, 24 retained a lower local-search candidate despite a non-success
termination flag, and one retained an embedded smaller-model candidate. All three
DE searches and their outcomes remain recorded. Exhaustive subset coverage should
not be confused with a proof of globally optimal continuous parameters.

The initial full sweep and validation took about 169 seconds wall time and 620
CPU-seconds (four allocated CPUs), excluding figure generation. A checkpointed
reproduction/figure run completed as SLURM job 20513998. The original sweep job
20513965 completed all inference outputs; its scheduled render ran before the
plotting file had been created, and rendering was subsequently completed and
verified through the successful checkpointed run.
