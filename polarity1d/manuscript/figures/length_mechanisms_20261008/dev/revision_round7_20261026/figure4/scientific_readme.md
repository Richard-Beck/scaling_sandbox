# Unexpected

The common weak input f = 0.005 lies near the C20 relocation transition in Goryachev but well below it in Wave-pinning and Holmes. Weak input therefore has very different correlation scales across facets. Strong-input activity is not always a single relocating domain; the Wave-pinning examples show a patch field. At strong input Spring correlation curves nearly coincide. At weak input Goryachev C40 and C60 remain near 1 through the measured window, so their relative rates are unresolved. The strong-input readout is orientation of whole-ring activity imbalance, not persistence of an established polarity domain. This distinction limits the polarity-to-migration interpretation.

## Evidence and definitions

No new simulations are used. Baseline curves and α confidence intervals preserve the round6 estimates. Model ACFs are recomputed from the same complete native compact histories at dense early lags, with the original |P| ≥ 0.1 endpoint rule and 2000 seed-bootstrap resamples; their original 96–384 seeds are preserved. Baseline fields and seed-0 Wave-pinning examples are re-exported as compact CSVs. C uses the original |P| ≥ 0.1 endpoint eligibility mask, while Figure5D intentionally has none. Native f is a fraction of the model reference input, not a common cellular dose.

## Baselines and model readouts

The three synthetic baselines have α = 2, 0 and −1 by construction. Physical drift is fixed across size in the first case. Jump rate is respectively fixed or proportional to C in the other two. P = M1/M0 is the normalized complex first moment of activity; heading is arg(P). The native local exponent is estimated in the original adjacent-size neighborhoods at a shared lag with eligible endpoint correlation and coverage. Missing α estimates are not censored values and receive no marker.

## Reproduction and provenance

Run `round7_persistence_render_20261026.py --figure4` in the existing scaling_sandbox Python environment. `supporting/metadata/provenance.json` hashes the unchanged source tables and exemplar archives. Data-generation timestamps are stored there. No dropped round6 panels or supplement are generated.

## Prospective cellular assays and predicted signatures

In matched WGD− and WGD+ cells, image full-perimeter polarity activity under matched fluctuating input and measure heading autocorrelation together with first-mode coherence. Classify a coherent single leading region separately from coexisting regions or a low-coherence patch field. The models predict longer coherent heading memory in larger WGD+ cells at weak input; activity-imbalance orientation in a patch field is a separate outcome.

Estimate heading-decorrelation rate against measured perimeter at weak and strong input. Near-quadratic persistence with perimeter is the unchanged-physical-drift baseline (Holmes/Spring-like); a larger weak-input increase paired with smaller physical boundary drift would favour Wave-pinning-like transport resistance. The size advantage should narrow as input strengthens and can become small, as in Spring. These are prospective model-class signatures, not a calibrated prediction of cell-velocity persistence.
