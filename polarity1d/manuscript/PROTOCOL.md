# Corrected polarity establishment and matched reversal

Updated 2026-10-03 for the user's authorized clean rerun. The active manuscript
uses five models: Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring. LEGI and
Otsuji are excluded from reversal. Historical linear-conditioned curves and
clocks are deprecated and supply no active manuscript results.

## 1. Evidence and current experimental chain

Completed preparation: independent native rest -> finite left raised-cosine
patch -> restore basal input -> verified release. New matched reversal:
qualifying selected full checkpoint -> mirrored finite right patch -> basal
recovery. Use nominal amplitude factor 1 and all three current width conditions
at the models' existing lengths. Never substitute old linear-conditioned states.

[data/establishment_basal](data/establishment_basal/) retains corrected
preparation endpoints, native checkpoints and audit evidence. The original
campaign is `polarity1d/four_assay/pilots/polarity_establishment_ref15_relative20_basal_20261003/`.
It completed 3,234 combinations without execution errors: 2,476 polarized and
758 unpolarized. All 462 LEGI cases return to basal homogeneity, supporting its
exclusion. Otsuji's reversal exclusion is the user's explicit scope decision.
Neither enters the new length-response summary.

The main figure uses the corrected chain and newly generated matched reversal
measurements. Until those measurements pass collection/audits, any displayed
reversal schematic is labeled proposed; old curves cannot fill missing results.

## 2. Preparation inputs and matching cue

The extra patch is `A * (1 + cos(pi * d / w)) / 2` inside support
`0 <= d <= w`, zero outside. Left distance is x; right distance is L - x.
The native peak A is 0.1 for Wave-pinning/Goryachev, 1 for Holmes 3, 3.6 for
deBelly and 4 for Spring. Extra input returns to zero at pulse cutoff;
intrinsic basal reactions remain active.

| Condition | Non-Spring models | Spring |
|---|---|---|
| `fixed_10pct` | Width 1.5, reference 15 | Width 0.1, native resting reference 1 |
| `fixed_20pct` | Width 3, reference 15 | Width 0.2, native resting reference 1 |
| `relative_20pct` | 20% of native cue-domain length | 20% of instantaneous physical length |

An oversized fixed support is evaluated on the available domain without
rescaling peak or width. There is no old half-domain restriction or equal-integral
amplitude condition. Spring's reference 1 is not calibrated microns. deBelly
uses its native nominal cue coordinates and retains displacement separately.

Selected preparation checkpoints come from the bounded duration search:
shortest tested left pulse within 1% of the strongest observed verified left
endpoint, with left strength and contrast >= 0.05. Release requires two passing
1,000-unit windows at fixed-scale range/projected drift <= 0.0001 and a 10,000-unit
holdout with excursion <= 0.001 and a passing final window; cap 60,000.
Each original duration starts independently from native rest. Zero selected left
cue duration means no polarizing cue was needed. Keep its actual interpretation.

## 3. Exact restart, eligibility and response

Verify selected checkpoint hashes, nominal factor 1, width identity and qualified
left state before submitting trials. Preparation-ineligible cases are not reversal
failures. Reproducibility checks may continue a copy under no cue and compare all
native fields with the recorded fixed scales; every response still starts from
the original selected checkpoint. No reflection/reconditioning/state replacement.

Every tested right-pulse duration independently restarts the same full state.
After cutoff recover for 10,000 native units. Both scalar gates evaluate all
201 original native samples in the final 500 units, with a separate field flag:

| Measure | Tail criterion | Interpretation |
|---|---|---|
| Orientation reversal | Orientation > 0.05 and contrast >= 0.05 throughout; orientation range < 0.005 | Stable observable right orientation over the measured window |
| Scalar strength replacement | Orientation gate plus every orientation within 1% of absolute initial left orientation | Comparable reversed scalar first-moment strength |
| Field-tail flag | Each monitored field pointwise range and projected 500-unit drift <= 0.1% frozen scale | Independent native-field stability diagnostic |

The field flag remains visible and participates in any explicitly stable-field
fit mask. Scalar success alone does not establish full-profile convergence or
replacement. Zero right-pulse control should not reverse a qualifying left state;
a passing zero control is a preparation-instability flag, not a reversal time zero.

## 4. Bounded nonmonotonic duration search and audit

Freeze the new model grids: zero plus inherited dense numerical duration
coordinates plus 17 logarithmic samples over these current preparation ranges.
Reuse no old outcomes. The inherited coordinates alone have preparation-dependent
bounds and are insufficient for the new protocol.

| Model | Dense coordinates | Added log range |
|---|---|---|
| Wave-pinning/Goryachev | 291 samples, 10–5,120 | 0.05–5,120 |
| Holmes 3 | 341 samples, 0.25–5,120 | 0.05–5,120 |
| deBelly | 58 samples, 10–1,000 | 0.5–2,000 |
| Spring | 111 samples, 0.02–5 | 0.002–50 |

Refine every sampled fail-to-pass boundary for both scalar gates with bounded
local sampling toward 1% duration brackets. If the first positive duration passes,
use bounded halvings and retain a lower-range unresolved flag when appropriate.
Do not bisect the full range assuming monotonic success. Keep all sampled success
intervals, nonmonotonic flags and unresolved brackets. A refined endpoint remains
the shortest **tested** passing duration: sampling can miss narrow intervals.
No-success cases are censored within the reported maximum; missing execution is
an error. Both gates share identical tested trials.

Collection audits verify coverage, all tail samples, selected minimum duration,
full-strength implication of orientation success, field flags, original checkpoint
hashes and frozen scientific sources. Save full final states/native fields and
record preparation provenance. Endpoint clocks exclude preparation and recovery.

## 5. Supplements and personal diagnostic boundary

Corrected establishment heatmaps retain selected left cue duration, total
establishment time, final polarity and apparent endpoint count. Total establishment
time is cue duration plus persistent scalar-polarity arrival within 1% (floor
0.000001), excluding verification; profile clusters are not exact attractors.
The detailed reversal supplement is rebuilt from new matched results only.

[figures/diagnostic_reversal_fields/index.html](figures/diagnostic_reversal_fields/index.html)
is the historical personal field explorer. Its old width references,
linear-conditioned inputs, LEGI zero-total-input artifact and Otsuji ambiguity
are explicitly diagnostic. Keep only its needed dependency/provenance closure.
It supplies no manuscript curve or fit. [ASSUMPTIONS.md](ASSUMPTIONS.md) records
native model conventions and practical interpretation limits.
