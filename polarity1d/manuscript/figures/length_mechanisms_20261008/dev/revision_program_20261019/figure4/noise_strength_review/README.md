# Existing weak/strong histories and native snapshots

Author-requested discussion artifact,19 October2026. These four pages expose
existing input regimes and local patterns without replacing either main
Figure4 sheet. No native integration, new forcing or representative reselection
was performed.

[Combined four-page PDF](noise_strength_native_comparison.pdf) ·
[All-size prevalence and sampled-run inventory](../regime_review_20261019/README.md) ·
[Exact sources and snapshot selection](supporting/data/snapshot_sources_and_selection.csv) ·
[Validation and hashes](supporting/metadata/validation.json)

## Pages

| Model and size | Native snapshot time | Weak / strong medoid seeds | Observed strong profile |
|---|---:|---|---|
| [Wave-pinning C20](noise_strength_wave_pinning.png) |16000|202610180000 /202610180018|One half-contrast region; the retained strong medoid has no multiregion sample.|
| [Holmes C60](noise_strength_holmes_model3.png) |9560|202610180195 /202610180205|Two half-contrast regions, three prominent peaks; R=0.795.|
| [Spring C0=3](noise_strength_spring_ring.png) |204|202610190455 /202610190471|Two half-contrast regions, two prominent peaks; R=0.159. The two plot edges join around the ring.|
| [Goryachev C20](noise_strength_goryachev.png) |17000|202610180080 /202610180086|Three half-contrast regions, three prominent peaks; R=0.103. Its uniform mean-input control has lost polarity.|

Each PNG has matching individual PDF/SVG exports under supporting/vectors/.
Every history is an unchanged within-condition descriptor medoid selected
from12 independent seeds. Weak and strong histories use different realizations;
they do not share a matched forcing path. The starting active fields are
checked to be exactly identical within each model/size comparison.

## Caption and interpretation

**A–B:** Full saved native active-field histories at the weaker f=0.001 and
stronger f=0.01 inputs, shown with one shared zero-based concentration scale
per page. Time runs downward. Dashed horizontal lines mark the exact saved
native time used in both profiles; each column states the actual time and seed.
The horizontal coordinate is material boundary position in cycles, not physical
distance traveled. The left and right edges are the same periodic location.
Chemical histories have801 saved fields at cadence40 over32000 native time
units; Spring has801 fields at cadence4 over3200. Differences in native units
and forcing statistics preclude calibrated cross-model time or amplitude ranks.

**C:** Actual initial, weak-input, strong-input and uniform stronger-mean-control
profiles. The initial curve is at t=0; the other three use the identical marked
snapshot time. Curves retain their native amplitudes without per-history peak
normalization. The Goryachev inset is explicitly a native-scale vertical
enlargement of the same curves; the full profiles and kymographs retain their
shared scale. Legend R values use the physical first moment normalized by
its own initial complex value:

\[
M_0(t)=\frac{C(t)}{N}\sum_{j=1}^{N}u_j(t),
\qquad
M_1(t)=\frac{C(t)}{N}\sum_{j=1}^{N}u_j(t)e^{2\pi i q_j},
\qquad
P(t)=\frac{M_1(t)}{M_1(0)},\quad R(t)=|P(t)|.
\tag{1}
\]

Spring includes its actual instantaneous physical circumference C(t).
R is relative first-moment magnitude, not total active amount or proof of
coherent profile translation. A low resultant may accompany multiple activity
regions while activity remains locally substantial. No movement law or site
ancestry is inferred.

The two declared inputs scale **positive mean activation and fluctuation
amplitude together**. With a unit-variance Gaussian field and fixed lognormal
shape parameter1, the input is f times the native reference amplitude times
exp(Z−1/2). Thus a tenfold increase in f raises both ensemble mean and standard
deviation tenfold. It is not a zero-mean variance sweep. Native reference
amplitudes are WP0.1, Goryachev6.4, Holmes1 and Spring4. Chemical correlation
scales are ell1.5/tau50; Spring has its separately declared ell0.15/tau5
transported-material law. The precise finite-knot interpolation conventions
remain those of the accepted campaigns.

Goryachev C20 requires the visible control qualification: under the stronger
uniform mean alone, the first moment collapses and the selected profile is
nearly homogeneous. The stronger noisy history therefore describes a
forced-pattern regime, rather than clean wandering or breakup of a pole that
would remain intact under the mean input. The other three selected mean
controls retain a single half-contrast region. No-extra-input controls and
all other sizes/strengths remain accessible through the inventory.

## Transparent snapshot selection

The histories were fixed by the existing medoid selection. Within the strong
medoid, Holmes, Spring and Goryachev snapshots are the saved field nearest
the midpoint of its longest consecutive scalar-sample run having multiple
current half-contrast regions, contrast at least0.05 and more than one peak
of prominence at least0.05. These are inherited native identification thresholds,
not a new functional competition classifier. **Prominence0.05 is an absolute
model-native concentration threshold**, not5% of a profile's range. It does
not establish equivalent secondary-site strength across models.

| Strong medoid | Qualified positive-sample run | Observed span | Chosen saved time |
|---|---|---:|---:|
| Holmes C60 |9385–9740|355|9560|
| Spring C0=3 |130–280.5|150.5|204|
| Goryachev C20 |16670–17300|630|17000|

WP uses the middle saved frame at16000 because its C20 strong medoid has no
multiregion sample. Weak and mean-control profiles use the **exact same saved
time**, without interpolation. These snapshots are selected diagnostics inside
unchanged medoids, not estimates of a typical marginal state or an example
chosen from the most extreme seed. Consecutive positive samples do not prove
uninterrupted multiregion organization between samples. Scalar cadence is5
for chemicals and0.5 for Spring; all selected field times coincide exactly
with scalar times.

The [all-size inventory](../regime_review_20261019/README.md) gives the needed
population context: no sampled multiregion occupancy at f=0.001; at f=0.01,
largest-size Holmes and Spring means are17.15% and14.12%, respectively, with
whole-seed uncertainty and heterogeneous episode spans. The selected Holmes
medoid itself has1.78% multiregion occupancy, illustrating why the picture is
not a substitute for group coverage. Multiple regions and prominent peaks
remain descriptive morphology; they do not establish competition, switching
events, uninterrupted lineage or migration.

## Reproduction and compact data

[Renderer](../../../../../../code/mechanisms/revision_figure4_strength_review_20261019.py)
reads only the audited native archives and existing medoid inventory.
The source/selection CSV records all12 source cases, exact times, scalar and
profile cadences, native circumferences, physical moments and archive hashes.
[Native snapshot profiles](supporting/data/native_snapshot_profiles.csv) retain
every plotted field sample and its physical quadrature weight; the
[compact arrays](supporting/data/native_comparison_arrays.npz) retain the
histories, profiles, material positions, C(t), P/R and counts used by the pages.

The automatic audit verifies all12 trace/profile source hashes, exact common
initial fields, exact saved times and physical P replay. It also checks the
six original main PNG/PDF/SVG hashes remain unchanged. The four pages were
visually inspected after rendering. Individual and combined exports, compact
data and renderer hashes are recorded in the validation JSON. Generation uses
one CPU, approximately14seconds and0.5GiB peak memory; no heavy computation or
new simulation was needed.
