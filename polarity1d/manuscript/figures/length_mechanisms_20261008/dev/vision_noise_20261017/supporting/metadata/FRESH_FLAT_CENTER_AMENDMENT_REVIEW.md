# Fresh review: flat-profile center validation

Reviewed 17 October 2026 at 20:53 UTC. **The narrow adapter amendment is
acceptable for the analysis-only retry.** This does not review or accept final
replicated scientific results. The earlier numerical-amendment review and
all original frozen sources remain unchanged.

The operative reviewed source is
`data/runs/vision_noise_main_20261017/results/statistics_centroid_amendment_20261017/analysis.py`,
SHA-256 `61622d8d21102397d38e9dc4050cc83852a60169c0e08fdedff494e1ab415dc7`.
I inspected its complete difference from the already reviewed numerical
correction. Only `load_case` and the top-level audit text in `analyze` change;
all correlation, uncertainty, directional, site-displacement and selection
functions are identical at the Python syntax-tree level.

An independent actual-data check isolates the earlier assertion to worker
0223, the Goryachev C20 higher-mean control. Its native scalar descriptors
become flat at $t=4145$, and the first saved flat profile is at $t=4160$.
There is then no current half-contrast component or measured site center.
The five-unit producer and forty-unit profile analysis carry different last
nonflat labels, separated by $3.00279\times10^{-5}$ native length units.
Requiring those historical labels to agree would compare different previous
observation times, rather than the same current identifiable site.

Across **all 243 raw scalar traces and saved primary-profile grids**, I verified
that the flat flags exactly coincide with zero current components, and that
the labels carry their preceding value during scalar flat intervals. There
are 5572 scalar flat samples and 697 saved flat profiles, all in worker 0223.
No flat sample qualifies for the unchanged physical endpoint assay: its
contrast must be at least 0.05 and it must have exactly one current region.
The original direction gate still uses the moment magnitude, independently of
the cached local-center label. No cached flat center enters a local displacement.

The change preserves the original $3\times10^{-10}$ equality tolerance on every
profile with an actual component, validates flat flags/counts/carry separately,
and exports excluded-frame counts and the flat-label discrepancy. It changes
neither eligibility nor any value used in statistical or representative-history
calculations. The source's regression covers four declared cases, including
the failing control; **it is not an all-243 center regression**. The retry still
requires the strict observable-center comparison on all 243 cases before
statistics can be exported. Final acceptance will verify that complete audit.

The [independent checks](FRESH_FLAT_CENTER_CHECKS.json) record the source hash,
all-243 raw flag/eligibility checks and the bounded scope. Failed job 21126773
and its original analysis source remain preserved. No integration, noise-path
change, threshold adjustment, tolerance relaxation or trajectory omission is
needed for this correction.
