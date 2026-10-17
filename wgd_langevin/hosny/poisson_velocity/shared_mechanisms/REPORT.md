# Parameter contrasts under shared KD mechanisms

The representative shared model changes directionality after both HIF and RhoA
KD. It is the best shared explanation at the simplest tested nonzero complexity.
It provides a readable reference for parameter comparisons, while the frontier
and ratio profiles expose how conditional those comparisons are.

## Baseline point estimates and cross-model sensitivity

| Parameter | Representative non-PACC | Representative PACC | Representative PACC/N ratio | Range across nonzero shared frontier |
|---|---:|---:|---:|---:|
| Memory tau (h) | 0.020 | 0.040 | 2.00 | 2.00–20.96 |
| Kick rate lambda (/h) | 0.752 | 0.545 | 0.725 | 0.696–0.809 |
| Kick magnitude scale A (um/h) | 428.4 | 447.3 | 1.044 | 0.100–1.044 |
| Directionality concentration kappa | 0.230 | 0.512 | 2.233 | 1.810–2.246 |

The concentration and rate contrasts remain in the same direction with similar
magnitudes across the frontier. Memory stays above one in these optimized point
estimates, but its magnitude varies greatly. Kick amplitude changes direction
across models. None of those observations substitutes for within-model profiling.

## Within-model ratio profiles

Fixing a baseline PACC/non-PACC ratio to one, and jointly refitting everything else
under the representative HIF-kappa/RhoA-kappa architecture, gives:

| Equal-parameter constraint | Increase in bank loss | Increase in independent loss, two seeds |
|---|---:|---:|
| Same memory | 0.00010 | 0.00001 to 0.00006 |
| Same kick rate | 0.00143 | 0.00127 to 0.00155 |
| Same kick amplitude | 0.00005 | -0.00025 to -0.00019 |
| Same directionality | 0.00489 | 0.00438 to 0.00438 |

The small negative independent amplitude penalties reflect simulation-bank /
optimization differences; they are not evidence that adding a constraint must
improve an exactly optimized model.

Within these assumptions, greater PACC directionality is the clearest parameter
contrast. A lower PACC kick rate is a weaker but consistent feature. A memory
difference is **not identified merely because all frontier point estimates show
one**: equal memory fits essentially as well after compensating changes. Kick
amplitude is also weakly identified. These are comparative descriptive statements,
not significance claims or confidence bounds.

The broad memory/amplitude profiles express a physical sampling degeneracy: when
kicks decay well within the 30-minute observation interval, decay time and jump
magnitude can compensate while preserving sampled displacement. The huge jump
amplitudes in the table are not measurements of instantaneous cell speeds.

For reference, representative HIF/RhoA kappa values are approximately 0.032/0 in
non-PACCs and 0.099/0.013 in PACCs. The presentation constraint shares which
parameters change; it does not impose equal KD response magnitudes or signs.

## Presentation choices

The main figure focuses on baseline parameter comparisons, model choice and
identification, rather than repeating the track/observation panels. Representative
predictions for all five trajectory metrics are retained in CSV, and the earlier
exhaustive overview remains available. The representative is a parsimony anchor,
not a declaration that more flexible models are unnecessary or that all summaries
are reproduced adequately.

No parameter averaging, model probabilities, confidence thresholds or additional
literature observations were introduced. Conclusions remain conditional on the
unresolved directional-index denominator and the stationary homogeneous model.
