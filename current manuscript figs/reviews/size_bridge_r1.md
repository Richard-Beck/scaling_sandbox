# Independent review: bounded Holmes size/separation pilot

17 October 2026. Fresh-context critique of the completed
[three-size result](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/holmes_size_bridge_20261017/supporting/metadata/RESULTS.md).
No integrations were run and no authored result-package files were edited.
The reviewed export manifest was generated at 15:52:55 UTC; its producer is finished.

## Decision

**This is a credible, bounded size-response observation, with a more specific
message than a change in cap count: the challenger replaces the incumbent at
C30, while the two sites remain distinct at C40 and C50.** All three challenges
use the same local cue and fixed pool densities. Consequently, this is a useful
bridge from the fixed-circumference separation experiment to an explicit
enlargement path. It is also a joint intervention in circumference, material,
incumbent width, absolute separation and both periodic routes. The result does
not identify which of those changes causes the difference.

**C30 is not the smaller-ring version of the fixed-C40, d12 merger.** Its old
cap disappears near its original position; the surviving cap remains near the
cue and broadens. The fixed-C40, d12 case instead fills the intervening valley
and ends with a cap between the two original sites. This distinction must
survive any condensation into a main figure.

The preparation, accounting and retained fields adequately support the scoped
pilot. No additional simulation is necessary before discussing these results
with the user. They strengthen the claim that a prescribed size change can
alter spatial competition, but they do not resolve the chemical mechanism or
the connection to migration persistence.

## What the actual histories show

The [early spatial figure](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/holmes_size_bridge_20261017/size_bridge_spatial.png)
now extends through 4000 s from cue onset, so it includes the C30 transition.
The [full-history companion](../../polarity1d/manuscript/figures/length_mechanisms_20261008/dev/holmes_size_bridge_20261017/size_bridge_full_history.png)
retains the complete 1000-s cue and 13,000-s recovery. Both are necessary to
interpret a terminal one-cap count correctly.

At C30, two distinct maxima are present at withdrawal and for most of the next
2000 s. The old maximum remains near −0.075 µm while declining: its Rac is
1.4166 µM at 1750 s off and 1.2941 µM at 1950 s off. The cue-local maximum
remains near 12 µm. By 2000 s off, Rac at the original position has fallen to
0.3886 µM and the sole resolved maximum is near 11.9 µm. The terminal cap is
centered at 11.95 µm, with the original basal width and plateau height to the
reported precision. This is loss of the old site followed by a cue-centered
one-cap pattern, rather than continuous motion of the old cap to the cue or
merger into an intersite cap.

The prominence-0.05 old maximum is last detected at 1985 s off and absent at
the next 1990-s sample. That threshold crossing is a diagnostic timestamp;
the spatial history, cap positions and raw decline supply the response-mode
interpretation. The gap rises during the subsequent expansion of the surviving
cap, so late gap filling alone must not be mistaken for the earlier C40 merger.

At C40 and C50, both old and challenged sites remain distinct throughout the
recorded recovery. Terminal old-site peaks are approximately 1.4847 µM; the
challenged peaks are approximately 1.4845 and 1.4842 µM. The old cap becomes
narrower than its corresponding basal cap. Thus incumbent retention concerns
the persistence of a localized region, not preservation of its initial width
or active amount.

Initial challenger activation is nearly equal along the size path. Integrating
over the same physical width-3 cue region and dividing by its width gives
mean Rac at 1 s of 0.657726, 0.657752 and 0.657779 µM for C30, C40 and C50.
At withdrawal the corresponding means are 1.96151, 1.96250 and 1.96321 µM.
The different final modes therefore arise after substantial challenger
activation in every case; this is not evidence that a larger cell simply
requires a longer stimulus to activate the challenger.

## The useful physical clue, and its limit

The retained profiles expose an architectural change that is more informative
than a generic travel-distance account:

| Circumference / separation (µm) | Basal cap width (µm) | Old-cap width at withdrawal (µm) | Terminal mode |
|---|---:|---:|---|
| 30 / 12 | 12.380 | 3.941 | Old-site loss; one cue-aligned cap |
| 40 / 16 | 16.538 | 8.226 | Two retained caps |
| 50 / 20 | 20.696 | 12.453 | Two retained caps |

Widths are the reported local half-prominence widths; the basal single-cap
definition coincides with its half-contrast width. The plateau remains near
1.484 µM as the basal domain broadens with circumference. Making a challenger
of comparable local amplitude leaves a much narrower residual old region at
C30. At C40/C50 the terminal challenger widths are similar, approximately
7.92/7.87 µm, while the terminal old-region widths grow from 8.52 to 12.73 µm.
This motivates a hypothesis about how enlarged domains allocate activity
between sites. It does not demonstrate a minimum sustainable cap width,
a material-capacity threshold or an independently identified inhibitory signal.

The total active Rac amount is also informative but is not itself conserved.
The terminal values are approximately 24.8259, 33.0450 and 41.3485 µM·µm,
compared with basal values of 24.8259, 33.1293 and 41.4326. Nearly the same
activity within each size is reorganized into one broad region or two narrower
regions. Conservation applies to active-plus-inactive Rac and the corresponding
Rho pool, not to the active Rac component alone.

The already completed C40, d12 experiment prevents a separation-only reading:
at the same physical d12, C30 loses its incumbent whereas C40 merges the two
regions. Moving the cue from d12 to d16 at fixed C40 then changes merger to
two-site retention. These existing comparisons are complementary. They do not
form a full causal partition: changing C at fixed d still changes the other
periodic route, material and the freshly prepared incumbent.

The specified size path is scientifically legitimate as a model of enlargement
at fixed concentrations and local kinetic parameters. Calling it a joint
intervention does not invalidate the size-associated outcome; it limits the
explanation of that outcome. Its total Rac/Rho amounts grow from 225/93 through
300/124 to 375/155 µM·µm, its two routes grow from 12/18 to 16/24 to 20/30 µm,
and the fixed 3-µm cue occupies 10%, 7.5% and 6% of the circumference.

## Qualification and numerical support

The two new basals start afresh from native homogeneous fields and qualify
under the declared unchanged full-field windows and 10,000-s holdout. I
independently checked their saved phase hashes, uniform starts and pool means,
exact phase chains, qualifying consecutive windows, holdout excursions,
checkpoint continuity and single-cap endpoints. The maximum scaled holdout
excursions are approximately 7.16 × 10⁻¹⁶ and 8.60 × 10⁻¹⁶ for C30/C50.
No failing preparation is silently treated as a response outcome.

All six saved challenge phases have the expected exact initial arrays and
sampling, prescribed cues, positive finite concentrations and conserved pool
amounts to floating-point precision. The integrated local cue rate is 1.5
µM·µm/s at each size. All three responses pass the stated late full-field
stationarity diagnostics. C40/C50 still exhibit small changes, with maximal
scaled late-window metrics of about 2.71 × 10⁻⁵ and 1.41 × 10⁻⁵, respectively,
against the 10⁻⁴ diagnostic tolerance. Passing that finite criterion is not a
stability calculation.

The same physical mesh spacing, 0.025 µm, is used throughout, and the fixed-C40
response has its separate finer-grid and small-imbalance evidence. The new
C30/C50 outcomes do not themselves have a two-grid challenge comparison or
matched signed-imbalance continuations. Those controls should not be borrowed
as if directly performed at the new sizes. This does not prevent presenting
the present bounded pilot; it limits any stronger robustness claim or precise
transition-time claim.

I verified 24 frozen-source hashes, four frozen C40 reference inputs, all 15
current export artifacts, the report-source and run-manifest identities, and
all 14 local result-document links. These are accounting and saved-field
checks, not independent integration replay or asymptotic analysis.

## What to add to the main candidate

Keep the fixed-C40 comparison primary when the main question is the effect of
challenger location: it changes cue position while holding the incumbent and
material fixed. The new three-size result should initially remain a clearly
linked supporting result and a direct sentence explaining its replacement
versus retention contrast. C50 adds confirmation along the chosen path, but
does not add a third distinct mechanism and need not occupy another main panel.

If the selected main figure must visibly include a size intervention, a
reasonable compact alternative is three columns: **C30/d12, C40/d12,
C40/d16**, each with the spatial history and basal/withdrawal/final profiles.
That would show replacement, merger and two-site retention while preserving
the controlled fixed-C distance contrast. Use a common 0–4000-s early clock,
label the 13,000-s-off terminal profiles, and retain C50 and complete histories
in support. This is a proposed layout for selection with the user, not a
recommendation to automatically enlarge the current main figure.

Neither layout demonstrates increased migration persistence. Keeping an old
Rac site despite a challenger is a possible component of directional memory,
but two retained caps can also represent competing directions. The model does
not specify how their activity produces whole-cell motion. Observed migration
summaries, conditional OU persistence estimates, interval reversal readouts
and these periodic spatial patterns remain distinct quantities.

## Anticipated user critique

**Anticipated user critique — hypothetical, not an actual user quotation:**
“This does connect the test to an actual increase in size. But the small cell
has not merged its fronts: the old one dies and the new one stays near the
cue. Do not put those two outcomes under the same one-cap label. The useful
point may be that a comparable challenger leaves a much narrower incumbent
in the small cell. Show me that physical organization. You have changed
material, the incumbent and the routes together, so do not tell me you have
identified separation as the cause. And retaining two Rac caps still does not
tell me why the enlarged cell migrates in one direction for longer. Keep the
main figure readable and tell me exactly which step of that argument the
data support.”

The present result answers the first request with actual spatial histories
and a measurable architectural comparison. It supplies a candidate bridge from
enlargement to altered competition. It leaves the chemical cause and behavioral
link open; another duration-versus-size plot would not resolve either gap.

## Review scope

Read the root, manuscript and development-output instructions; the current
report status and four-assay README; the working figure context, manuscript
workflow/protocol/assumptions and three-test guidance; the migration figure's
source caption and OU results; both current main candidates and the preceding
distance/robustness reviews. Inspected all three current size figures, the
result document, source, audits, manifests, CSV exports and complete native
cue/recovery arrays, together with the new saved preparation phases. The
review evaluates the finished pilot and presentation, not the selection of a
final manuscript figure.
