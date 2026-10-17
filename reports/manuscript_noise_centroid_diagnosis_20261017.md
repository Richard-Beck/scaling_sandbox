# Read-only diagnosis: flat-state centroid validation

Checked at approximately 20:49 UTC on 17 October 2026. No native archive,
producer, frozen analysis or statistics source was modified by this diagnosis.

SLURM statistics job 21126773 failed with exit 1:0 after 54 s. Its failed
assertion compares `dominant_region_position` from the producer's scalar traces
with the independently reconstructed position from saved profiles, with a
3 × 10⁻¹⁰ native-length tolerance. The producer and all 243 native integrations
remain complete and pass the independent execution/noise-regeneration audit.

The failing case is slot 0223: Goryachev, circumference 20, deterministic
uniform mean input 0.064 (fraction 0.01 of the native cue peak). Direct
recomputation using the failed analysis module finds:

- All nonflat saved profiles agree **exactly** in dominant-region position.
- 697 saved profiles fail, beginning at t = 4160. Every failing frame has
  contrast at most 10⁻¹² and zero current half-contrast regions.
- The native cached position is 0.024942652872420003; the independently
  reconstructed cached position is 0.02491262497562861, giving a circular
  difference of 3.0027896791295916 × 10⁻⁵.

Both implementations deliberately carry the last valid position after no
current region is resolved. The native producer updates this cache at scalar
cadence 5; the independent reconstruction sees saved profiles at cadence 40.
Their last resolved observations can therefore differ. The cached value is
not a current physical activity-site centroid in a numerically flat state.

The bounded correction is to require the original strict centroid agreement
on frames with current nonflat regions, and separately validate finite/cache
behavior at each declared cadence in flat intervals. Preserve those intervals
as flat/unavailable for site interpretation. Do not broaden the tolerance for
resolved profiles, assign a heading or site in the flat state, change the
correlation estimator, or reintegrate the verified native ensemble.

Any corrected statistics implementation should remain separate from the
failed frozen workflow source, preserve that failure record and demonstrate
the resolved/nonflat identity and expected cache behavior before a replacement
statistics submission. The direct scalar-cadence site diagnostic already
requires contrast at least 0.05 and one resolved region, so these flat labels
are ineligible for its physical displacement comparison.
