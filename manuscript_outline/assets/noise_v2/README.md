# Noisy-ring interpretation supplement, outline V2

This render-only supplement uses the completed 243-case campaign: 216 noisy
realizations and 27 deterministic controls. It adds no native integration and
does not modify the original data, estimator, or canonical figures.

- `noise_memory_geometry_v2` is the proposed continued main sheet: actual
  restricted correlation areas, direction, final magnitude and matched physical
  displacement/coverage at both strengths for all three models.
- `noise_adaptation_v2` retains full, first-half and last-half diagnostics at the
  same lag horizon, without a selected burn-in.
- `noise_native_examples_v2` shows exact high-input C20/C60 within-condition
  medoids. The compact NPZ retains every saved native primary profile and the
  full scalar traces used for these six examples. The selection JSON records
  all candidate distances, relative native paths and original source hashes.
- `noise_interpretation_v2.content.json` supplies the detailed captions,
  narrative, measurement suggestions and limits for outline integration.

Rebuild the three PNG/PDF sheets from the retained, committed inputs:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python -B scripts/render_noise_interpretation_v2.py
```

The default reads only this folder's compact CSV/CSV.gz/NPZ/JSON files and does
not access the original workers or bulk analysis tables. The same script in a
relocated checkout resolves these inputs relative to its own repository root.
The environment needs NumPy and Matplotlib; this repository already has them in
the `scaling_sandbox` Conda environment.

The optional `--audit-native` mode requires the original local campaign. It
reexports the compact inputs and checks every native trace hash and final
magnitude, independently recomputes the lag-8000 direction/site values in all
54 noisy condition/window combinations, and reintegrates all 54 restricted
areas from every retained source lag. The completed one-time results are in
`noise_reconciliation_v2.json`; the default also rechecks the 54 integrals and
compares every plotted direction, site and magnitude source row with its
retained export. Stored 95% whole-seed bootstrap intervals retain their source
values and hashes; a different bootstrap random stream is not substituted.

The scientific claim is conditional activity memory under the declared input
law. Angular geometry, forced-state changes and single-site eligibility limit
its interpretation as physical site motion. In particular, Goryachev's sparse
physical comparison is unavailable for a population curve and is never drawn
as zero drift. The data do not establish a cell velocity law, causal mechanism
for migration persistence, or an extrapolated stationary persistence time.
