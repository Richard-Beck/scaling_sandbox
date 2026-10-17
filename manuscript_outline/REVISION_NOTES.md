# Second outline: scientific integration and validation

The initial edition was committed and pushed immediately after integration as
`6199cdc33fea355133e4ede149131455dc5ccb91`. Its 31-page PDF, input document and
build records are preserved. The subsequent independent
[coherence review](../reports/manuscript_coherence_review_v1_20261017.md) identified
substantive interpretive and presentation repairs. This second edition makes
those repairs using completed results; no new integrations or fitted parameters
were required.

The hypothesis remains that WGD-associated enlargement can alter polarity
organization and replacement in ways relevant to migration persistence. Length
is a between-cell size comparison under fixed local kinetics. It does not stand
for elongation during a turn. Deterministic replacement, noisy activity memory,
observed migration summaries and conditional OU velocity memory remain separate
evidence. Their connection is an experimental hypothesis.

## Repairs made

1. **Mechanisms and measurements.** Two readable continued Figure 2 sheets join
   each model's mechanism, evidence and corresponding proposed measurement.
   The detailed legend distinguishes evolved mobility/deactivation interventions,
   corrected deBelly resting-state ablation, endpoint budgets and descriptive
   clocks. The common gain/recovery assays are useful matched observables, not
   a causal partition of all three cue-duration curves. Corrected deBelly's
   nonmonotone full curve precedes selected rising-range explanations. Separate
   Goryachev perturbations retain their own clocks and recovery behaviour.
2. **Physical interpretation of noise.** Figure 4 now includes both input
   strengths, signed finite-window areas, direction-only memory, final magnitude,
   physical single-site displacement, exact eligible-pair coverage and all
   full/first/last windows. Matched C20/C60 native examples use declared medoids
   and shared within-model color scales. Larger Wave-pinning and Holmes rings
   retain more angular memory at high input. Holmes does not establish smaller
   physical displacement; its site coverage decreases. Goryachev's small-ring
   forced-state change and unavailable population drift remain explicit.
   Undecayed finite-window areas are not assigned fitted persistence times.
3. **Migration panel B.** The author-requested observed points/model lines and
   local curve-shape display are implemented. Points sample digitized raster
   curves; the original experimental lag grid is unknown. Matched log-lag
   windows expose the observed/model curvature discrepancy and different
   crossovers, approximately 2.32 versus 1.11 hours. Numerical fitting inputs,
   parameters and other panels remain unchanged. The detail sheet preserves
   legibility; slopes are descriptive new calculations, not author-reported
   exponents or independent observations.
4. **Reproducibility and numbering.** A portable wrapper reconstructs only the
   original checkpoint's 64 declared inputs from the first Git publication.
   Original replay reproduces all four supplied PNGs byte-for-byte. The numbered
   variant changes only title numerals, with pixel checks excluding every
   scientific panel. The new noise renderer uses compact retained inputs by
   default; optional native auditing is separate. The first edition's original
   renderer is unchanged.

The approved corrected five-model, three-challenge Spearman comparison remains
the foundation. Full curves, censoring, preparation ineligibility and separate
native-field flags remain visible. The major SUM-159 WGD+/− motility-assay
section is reserved in Figure 6; no experimental results are invented. Proposed
mechanism measurements remain alongside their theoretical motivation rather
than substituting for the new motility results.

## Rebuild and audit scope

Use the existing `scaling_sandbox` environment. From the repository root:

```sh
python -B scripts/render_mechanism_crosswalk_v2.py
python -B manuscript_outline/prepare_revision.py
python -B scripts/build_manuscript_outline_v2.py --input manuscript_outline/outline_v2.json --output manuscript_outline/outline_v2.pdf
python -B scripts/verify_manuscript_outline.py manuscript_outline/outline_v2.json --require-tracked --report manuscript_outline/outline_v2.validation.json
```

The retained PNGs make outline regeneration independent of simulation archives.
The second-edition builder extends the frozen first renderer with explicit
foldout page sizes for dense figures and verifies at least ten inches of printed
image width, preserving the crosswalk's minimum nine-point type. It changes no
scientific content or first-edition file.
For source-level figure replay, see
[portable mechanisms](PORTABLE_MECHANISM_RENDER.md),
[migration preservation checks](assets/migration_v2/validation.json) and
[noise provenance](assets/noise_v2/noise_interpretation_v2.manifest.json).
Default noise replay is `python -B scripts/render_noise_interpretation_v2.py`;
`--audit-native` separately requires the local full campaign archive.
These checks validate source identity, calculations and packaging; they do not
establish biological validity or a unique polarity-to-velocity mapping.

Build manifests record exact current inputs. The first edition's recorded
hashes refer to its published Git snapshot: later updates by the separate
checkpoint producer can change live documentation. Verify historical inputs
against `6199cdc`, rather than rewriting the first PDF or its original manifest.

The first publication passed its complete GitHub Pages build and deployment
([run 39069424354](https://github.com/Richard-Beck/scaling_sandbox/actions/runs/39069424354)).
The second publication's final audit records its separate checks and deployment
status. All large native outputs remain local with measured storage and purpose
in [DATA_MANIFEST.md](../DATA_MANIFEST.md). Reviewed source and compact evidence
are committed; the preexisting local deletion set is preserved as workspace
state and is not silently committed or reversed.
