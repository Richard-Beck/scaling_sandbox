# Fresh-context review of candidate 1, revision 2

Reviewed 17 October 2026 against the actual six-panel PNG, candidate discussion,
technical caption, exported measurements, native reaction definitions and retained
returned-state arrays. Repository and manuscript guidance, current comparison
status and migration-figure context were read before interpretation. No native
integration or candidate-file edit was performed.

**Decision: retain revision 2 as a useful provisional main-text candidate. A panel
substitution is not essential.** The profiles, native laws, actual damage/recovery
traces and returned-state allocation now form a coherent mechanism. The earlier
intrinsic-restoration interpretation has been corrected explicitly. One focused
caption/discussion addition should make the remaining inference concrete: identify
how much of the length-dependent substrate response comes from a larger native
conversion coefficient, rather than greater exposure to inactive surplus. The
retained states answer that question without another experiment or diagnostic row.

**Decisive point: state the property that grows.**

The native coefficient \(A(u)=\partial f/\partial v\) is the instantaneous local
activation response per increment in inactive concentration at fixed active
concentration. This is the physical meaning of “capture” here: conversion of the
shared inactive pool into active material. It is not a separate binding process or
a measured force. The current final row displays a product of this coefficient
with the particular returned inactive surplus. Calling the product stronger
capture is defensible for this disturbance, but leaves the user's central question
partly implicit.

The exact retained-state factorization is

$$
g_v
=\frac{\displaystyle\int_{\Omega_0} A(u^p)\,\Delta v\,\mathrm dx}{\mathcal D(0)}
=\underbrace{\frac{\displaystyle\int_{\Omega_0} A(u^p)\,\Delta v\,\mathrm dx}
{\displaystyle\int_{\Omega_0}\Delta v\,\mathrm dx}}_{\overline A_{\Delta v}}
\underbrace{\frac{\displaystyle\int_{\Omega_0}\Delta v\,\mathrm dx}{\mathcal D(0)}}_{\text{core surplus per core active deficit}}.
\tag{1}
$$

Here the core, deficit, perturbed fields and matched sham are exactly those in the
candidate. The first factor has units of inverse native time; the second is
dimensionless. The inactive surplus in these cores is positive. This is an exact
description of the actual coupled state, not an independently evolved intervention.

| Model | Length | Surplus-weighted conversion coefficient, \(\overline A_{\Delta v}\) | Core surplus / core deficit | Product, \(g_v\) |
| --- | ---: | ---: | ---: | ---: |
| Goryachev | 9 | 1.87702 | 0.531411 | 0.997471 |
| Goryachev | 39 | 31.55190 | 0.485196 | 15.30885 |
| Wave-pinning | 9 | 0.560248 | 0.928784 | 0.520350 |
| Wave-pinning | 39 | 0.612353 | 0.984652 | 0.602955 |

Thus Goryachev's substrate contribution grows **15.35-fold** through a
**16.81-fold** increase in its surplus-weighted conversion coefficient, while
surplus per core deficit actually falls to **0.913 times** its initial value.
Wave-pinning's coefficient grows **1.093-fold**, exposure **1.060-fold**, and their
product **1.159-fold**. These are substrate-term ratios; the already reported net
rate ratios remain 15.64 and 1.25.

The four full-precision endpoint rows are retained in the
[factorization CSV](incumbent_r3_capture_factorization.csv).

The coefficient increase also exists in the undisturbed starting states, so this
is not solely an artefact of surplus weighting. The spatial mean of native
\(A(u)\) over each original core rises from **1.83350 to 31.54567** in Goryachev,
versus **0.551040 to 0.608498** in Wave-pinning. The larger Goryachev peak therefore
changes the local conversion response strongly. Wave-pinning adds chiefly spatial
extent and active amount while preserving a similar plateau concentration.
Broadening still increases its integrated activation response; a nearly unchanged
local coefficient does not make its larger incumbent functionally irrelevant.

One concise addition giving the coefficient and exposure fold changes is enough
in the candidate discussion and caption. The existing laws and profiles already
show the state change, and E–F show its consequence under the specified challenge.
Equation (1) and the numerical table can remain linked technical evidence. I do
not recommend replacing a panel merely to display this factorization.

**Decisive boundary: describe the observed state response, not saturation or
competition as an isolated cause.**

The sign of recovery after returning substrate and the global mass identity are
expected consequences of this protocol. They are checks, not the principal new
finding. The informative result is the strong length dependence above, coupled
with the lower remote inactive concentration in Goryachev. That combination has
a clear possible consequence: stronger incumbent conversion of available inactive
material alongside a smaller initial added source from an equal local remote cue.

The current scope statements correctly stop before successful rival suppression,
takeover, general stability or migration persistence. Preserve them. No rival is
present in the removal/return assay, and lower initial cue-driven activation does
not determine the rival's subsequent fate. The two chemical models remain
candidate physical explanations motivated by WGD-associated migration persistence,
not a demonstrated bridge to observed migration or conditional OU persistence.

The phrase “Wave-pinning's factor saturates” is mathematically correct but should
not become “its sampled incumbent is already fully saturated.” At its observed
peaks, \(A\) is approximately 0.615–0.626, well below its asymptote of 1.067.
Its coefficient changes little chiefly because the observed plateau concentration
changes little. The model laws explain the consequences of the sampled states;
the comparison does not independently isolate saturation as the cause of their
different morphologies.

**Optional polish; no additional experiment is required.**

- The exact chemical split is useful evidence, but its ordering need not dominate
  the main narrative. I checked both natural orderings across all 36 returned
  states. Goryachev's integrated active-field term remains negative in both;
  Wave-pinning's remains positive. The cross term is at most 0.0020% of the
  Goryachev substrate contribution and 0.0431% of the Wave-pinning contribution.
  Keep the formal ordering in technical methods and state the biological result
  directly in the main discussion. The sign conclusion is not fragile here.
- The PNG is legible and its column/color logic works. Four allocation symbols,
  two nearly invisible Goryachev terms, a conservation sentence and repeated
  qualifications give accounting more prominence than it needs. Small terms can
  stay as they are for this reviewable candidate; diffusion and ordering detail
  can move out of the main figure face in a later publication edit. Adding another
  panel or row would worsen the balance.
- Keep “immediate” attached to the 15.64-fold rate result. It is not a 15.64-fold
  reduction in full recovery time, a finite-damage tolerance, or an eigenvalue of
  the full spatial system. The actual traces appropriately show the intervention's
  unequal and very small damage.
- Preserve the explicit Goryachev end-view label and the distinction between
  concentration and total amount. Do not make its small width change a central
  mechanistic result on the retained length-dependent 201-cell spacing.
- In the technical caption after equation (5), “the second line” should refer to
  the diffusion line, which is the third line.

**Anticipated user critique — explicitly hypothetical, not actual quotations.**

- “What property of the larger incumbent grew? The new wording still makes me
  work through an accounting identity to find the answer.” The conversion-factor
  and exposure numbers above answer this directly without another figure.
- “Of course activation resumes after putting substrate back. Why is this an
  interesting size effect?” The approximately 17-fold coefficient increase,
  despite slightly less core surplus per damage, is the nontrivial distinction.
- “The broad Wave-pinning site contains more active material too. Do not imply
  it has gained nothing because its peak is unchanged.” Distinguish local
  conversion sensitivity from integrated extent and total amount.
- “This could favor the old site, but you have not shown a competing site losing.”
  Correct: the proposed consequence is mechanistically motivated, while the
  measured outcomes remain conditional rebound and initial remote activation.

**Evidence and reproducibility of this review.**

Numbers were reconstructed directly from `supporting/data/length_measurements.csv`
and its linked native `epsilon_0.005.npz` and `epsilon_0.npz` archives under
`polarity1d/manuscript/data/runs/inactive_pool_restoring_response_20261008/`.
The native 201-cell quadrature and fractional edge weight used by the candidate
were retained. Equation (1), basal core means and reversed-order cross terms were
evaluated algebraically from those arrays using the native laws in
`code/core/models/goryachev.py` and `wave_pinning.py`. No solver or integration
routine was called, and original evidence was not modified.
