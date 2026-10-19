# Figure 2: spatial-resolution qualification of the response options

Independent review, 2026-10-19T17:11:25.458420+00:00. This note preserves the already hashed main independent review and examines only the selected native response lineage and existing cross-checks. No integration or grid campaign was run.

The proposed qualitative argument is supported for the declared native model realizations: the same native activation-coefficient channel produces different size-dependent finite activity responses, and a declared small active-to-inactive kick produces different finite original-region deficits. The exact sampled profiles, coefficients, regions, clocks and native mesh must accompany that claim. These are conditional model-response options, rather than spatially converged quantitative predictions or calibrated optical protocols.

## What is and is not checked

The selected L19/L39 activation and recovery states use 201 cells. Physical cell spacing therefore changes from 19/201 = 0.09453 to 39/201 = 0.19403. Goryachev's interpolated original half-contrast core widths are 0.71419 and 0.79923, containing eight and four native sample centres. The native width ratio 1.119 and early repaired-fraction ratio 3.439 are point estimates on those discretizations, with spatial error not quantified.

Independent replay verifies the actual saved activation fields and all three paired recovery amplitudes, conservative kicks, zero-kick identity, tighter-solver comparisons, wider observation windows and measured deficit crossover. Those checks bound errors in the saved native response and its temporal integration. They do not refine the spatial mesh. The source explicitly says native meshes remain unchanged and no mesh campaign was added: [all_five_recovery_pilot_report_20261017.py](../../../../../../code/mechanisms/all_five_recovery_pilot_report_20261017.py) (line 177) and [RESULTS.md](../../../vision_linear_20261017/RESULTS.md) (line 68). No controlled spatial refinement for these exact L19/L39 activation/recovery experiments was found in the inspected lineage and report/code indices.

The roughly 0.2% difference in sampled cue peaks is much smaller than the roughly 35% observed excitation contrast. This rules out a simple proportional cue-sampling explanation; it is not a spatial-error bound or a proof excluding nonlinear sensitivity. The instantaneous native source is linear in the coefficient and inactive substrate, but the finite fields also include reaction and transport. Matched excitation-amplitude sensitivity is unavailable in this selected set. Consequently, a cue-peak correction cannot establish continuum convergence or isolate the causal contribution of substrate depletion.

The recovery plot retains raw initial deficit and its evolution: its whole-core loss is only about 0.22% for WP and 0.285% for Goryachev. Goryachev's restored-fraction fold weakens from 3.439 at the common early clock to 2.040 at 0.1 and 1.281 at 1. These finite native observations address the overly narrow early-window concern without establishing reporter precision, optical timing, spatial convergence or follow-up beyond 1.

## Existing corroboration with a different scope

The already available [ring_linear_common_recovery_pairs.csv](../../../../../../data/runs/vision_ring_recovery_20261017/results/ring_linear_common_recovery_pairs.csv) and [ring_linear_common_recovery_pairs_audit.json](../../../../../../data/runs/vision_ring_recovery_20261017/results/ring_linear_common_recovery_pairs_audit.json) compare the same active-loss definition and native clock between independently prepared interval states and ring states at C = 2L. The four Goryachev pairs are:

| Interval L / ring C | Relative ring-minus-interval K difference |
|---|---:|
| 5 / 10 | +0.186% |
| 15 / 30 | -0.175% |
| 25 / 50 | -0.686% |
| 35 / 70 | -1.016% |

These pairs provide corroboration of the qualitative restoration trend on different native meshes. They use independent preparations and different ring/interval domains, omit the selected L19/L39 pair, and compare the early coefficient rather than the full 1-unit recovery response. They are not a controlled spatial-refinement experiment. The existing [vision_ring_report_20261017.py](../../../../../../code/mechanisms/vision_ring_report_20261017.py) (line 382) explicitly makes this distinction; the largest Goryachev difference also exceeds its paired temporal numerical error.

The exact mirrored-operator fixture is a separate algebraic/numerical identity. It does not turn independently prepared, unmatched meshes into a convergence series.

## Recommended claim boundary

Keep the original all-five organization material and these finite native responses available for author discussion, with actual requested/sampled input, native support, frozen-region convention and four-centre limitation visible. The data motivate measurements that may distinguish the stated mechanisms under a defined perturbation. Do not present precise native width, early recovery or size folds as resolution-independent PDE predictions, a unique mechanism identifier or an experimentally calibrated effect size. A stronger spatial claim remains unqualified by the inspected evidence; this note does not authorize or request another campaign.
