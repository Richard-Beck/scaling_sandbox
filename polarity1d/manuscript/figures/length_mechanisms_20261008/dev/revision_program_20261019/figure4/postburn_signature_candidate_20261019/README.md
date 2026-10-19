# Figure 4 signature candidate: second-half observations

[Memory versus circumference](postburn_memory_coverage.png) ·
[Temporal direction curves](postburn_temporal_memory.png) ·
[Paired temporal coverage](postburn_temporal_coverage.png) ·
[Two-page temporal PDF](supporting/vectors/postburn_temporal_curves_with_coverage.pdf) ·
[Plot data](supporting/data/postburn_memory_coverage_points.csv) ·
[Generation/source receipt](supporting/metadata/MEMORY_RENDER_RECEIPT.json)

[Short native histories with linked profiles](postburn_local_examples.png) ·
[All108 model/size/strength trajectory examples,12 pages](supporting/vectors/all_conditions_second_half_native_atlas.pdf) ·
[All-strength temporal/window supplement,9 pages](supporting/vectors/all_models_sizes_strengths_temporal_windows.pdf) ·
[Organization, amount/coherence and widths,9 pages](supporting/vectors/all_models_sizes_strengths_organization_windows.pdf) ·
[Late-quarter sensitivity,2 pages](supporting/vectors/late_quarter_sensitivity.pdf) ·
[Native selection/provenance](supporting/data/all_condition_native_atlas_selection.csv)

This is a separate discussion candidate, generated on 19 October 2026. It does
not replace the accepted main sheets or the comprehensive
[broad-strength report](../broad_strength_campaign_20261019/review_complete_scientific/report.html).
All native production is terminal. The underlying complete-scientific data were
generated at 20:02:33 UTC; each new export records its own generation time.

The first candidate asks how activity-heading memory varies with circumference
across the input regimes already sampled. We discard the first half of each
trajectory and use chemicals t=16000–32000 and Spring t=1600–3200. This is an
explicit observed-window choice, not proof of equilibrium or stationarity.
Quarter/half-window checks and instantaneous spatial-width summaries were
completed in read-only SLURM job21200796. The spatial/window supplement passed
an independent source and circular-geometry audit. The local spatial examples
below use stored native fields; no additional native integration was performed.

The compact narrative being developed has two parts: circumference and input
regime jointly modify activity-heading retention; loss of heading memory can
accompany distinct local organizations. The current memory/coverage sheet and
short native histories are candidate main evidence. Detailed organization,
region-width, controls and observed-window suites remain linked supplements.
No competition or migration mechanism is assigned from these descriptors.

The upper row reports the pooled cosine of the heading difference at native
lag8000 (chemicals) or800 (Spring), using pairs whose two endpoint magnitudes
R are at least0.1. The lower row retains the number of usable pairs divided by
all available pairs. A direction estimate is unavailable below100 eligible pairs
or1% coverage. The selected60 points all meet this reporting gate, including
the lower-coverage Goryachev points. These are conditional activity-heading
correlations; they are not cell-migration autocorrelations or estimates of
uninterrupted site identity.

Five common strengths are shown:0.001,0.003,0.01,0.03,0.1. The original0.001/0.01
conditions have twelve seeds (filled squares); newer strengths have two (open
circles). Thin whiskers are the range of individually supported seed estimates,
not confidence intervals. Group estimates pool eligible cosine sums and counts,
so their estimand differs from an unweighted average of seed correlations.
Connecting lines guide the eye between three sampled sizes. Spring's size axis
is the actual physical basal circumference; physical activity moments continue
to use the evolving circumference C(t).

The supplemental temporal sheets show all six stored lag checkpoints at every
model/size for these same five strengths. Read direction and coverage together;
the fixed-lag size comparison alone cannot identify a decorrelation timescale.
Native times and forcing correlations differ between the chemical and Spring
models, and are not experimentally calibrated. Increasing f scales positive
input mean and fluctuation amplitude together.

The useful provisional contrast is that loss of directional memory can coexist
with high direction-pair coverage in WP/Holmes/Spring, whereas high-input
Goryachev increasingly loses usable direction with size. Local organization,
active amount and coherence must distinguish the underlying profile changes.
There is no assumption that multiple regions compete or that a polarity heading
defines cell movement. The link to migration is a testable activity-persistence
analogy; it is not a velocity law or a prediction of migration persistence.

Spring's original selected upper-strength numerical receipt remains failed for
sampled multiregion-run q90. No duration or event-rate inference is made here.
The supported scope and the unchanged failed evidence remain explicit in the
comprehensive report. Following separate root acceptance of the ten-metric scope
and independent exact-manifest review, replication was submitted as21200762 at
20:38:41 UTC. The present plots retain the completed576-case pilot snapshot;
they do not include subsequent replication cases.

The renderer is
[revision_figure4_postburn_signature_candidate_20261019.py](../../../../../../code/mechanisms/revision_figure4_postburn_signature_candidate_20261019.py).
It imports no native solver and performs no native integration. Every plotted
memory/coverage value comes directly from the independently audited
[all-window memory table](supporting/data/source_memory_all_windows.csv), whose
source SHA is fixed in the renderer and receipt. PNG/PDF/SVG exports and compact
plot tables are retained together under this folder.

The short local-history sheet uses the largest original size at f=0.1, showing
the fixed window T/2–5T/8 and snapshot9T/16. Its initial, noisy and deterministic
mean-control profiles retain their native amplitudes. Fields are sampled every40
native units for chemicals and4 for Spring. The Goryachev enlargement is on a
native concentration scale; a count of small threshold regions in that profile
must be read with its5%-of-basal contrast. Spring uses material coordinate q,
without treating stretch as translation. The actual C(t) remains in all physical
moments. Snapshots were fixed by time, with no search for extra peaks.

Condition medoids minimize the sum of pairwise L1 distances across seven
second-half descriptors (R,A,Q,contrast,qualified multiple-region fraction,
weak-R fraction,low-contrast fraction), each divided by its within-condition
range, using scale1 for a range below1e-12. Ties within1e-12 select the smallest
seed. With two paths both tie: this does not establish typicality. The twelve-page
atlas displays one such unchanged path at every model/size/strength, with a
shared zero-based native activity scale across the nine strengths on each page.
Per-panel percentages describe the selected path, not ensemble occupancy.
Every seed's scalar summaries and all native archive identities remain in the
comprehensive report. The full score table and exact main-field/profile arrays
are retained under supporting/data.

The nine-page temporal supplement covers all108 noisy condition groups, all
nine strengths, the full/first-half/second-half windows and all six lag
checkpoints. It separates conditional heading correlation, usable-pair coverage
and uncentered first-moment correlation. The latter uses all pairs and lag-energy
normalization; it includes magnitude/forced-state adaptation and remains a
descriptive supplemental quantity outside the narrowed ten-metric numerical
acceptance. It is not a centered stationary-process autocorrelation.

Sources and output identities are recorded separately in
[LOCAL_RENDER_RECEIPT.json](supporting/metadata/LOCAL_RENDER_RECEIPT.json) and
[FULL_TEMPORAL_RENDER_RECEIPT.json](supporting/metadata/FULL_TEMPORAL_RENDER_RECEIPT.json).

The organization supplement retains raw and contrast-qualified multiple-region
occupancy, low contrast, activity magnitude, active amount/coherence and matched
ensemble-mean controls for every strength, size and full/half window. Conditional
region widths join periodic seam-crossing arcs and use each frame's physical
circumference. Largest widths average qualified frames; second-largest widths
average qualified multiple-region frames; condition summaries weight contributing
seed descriptors equally. Qualification coverage and seed contributors accompany
the width tables. Native cell widths and profile cadence limit spatial and
temporal resolution. Late-quarter differences use paired seeds and explicitly
retain finite-window sensitivity rather than certifying stationarity.
