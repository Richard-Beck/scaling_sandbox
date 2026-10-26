# Unexpected

The common weak input f = 0.005 lies near the C20 relocation transition in Goryachev but well below it in Wave-pinning and Holmes. Weak input therefore has very different correlation scales across facets. Strong-input activity is not always a single relocating domain; the Wave-pinning examples show a patch field. This figure describes heading persistence and leaves the drift/jump interpretation to Figure5.

## Evidence and definitions

No new simulations are used. Baseline curves and α confidence intervals preserve the round6 estimates. Model ACFs are recomputed from the same complete native compact histories at dense early lags, with the original |P| ≥ 0.1 endpoint rule and 2000 seed-bootstrap resamples; their original 96–384 seeds are preserved. Baseline fields and seed-0 Wave-pinning examples are re-exported as compact CSVs. C uses the original |P| ≥ 0.1 endpoint eligibility mask, while Figure5D intentionally has none. Native f is a fraction of the model reference input, not a common cellular dose.

## Baselines and model readouts

The three synthetic baselines have α = 2, 0 and −1 by construction. Physical drift is fixed across size in the first case. Jump rate is respectively fixed or proportional to C in the other two. P = M1/M0 is the normalized complex first moment of activity; heading is arg(P). The native local exponent is estimated in the original adjacent-size neighborhoods at a shared lag with eligible endpoint correlation and coverage. Missing α estimates are not censored values and receive no marker.

## Reproduction and provenance

Run `round7_persistence_render_20261026.py --figure4` in the existing scaling_sandbox Python environment. `supporting/metadata/provenance.json` hashes the unchanged source tables and exemplar archives. Data-generation timestamps are stored there. No dropped round6 panels or supplement are generated.
