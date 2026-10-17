# Prospective supported-ring design critique

Reviewed 17 October 2026 against the repository instructions,
`MANUSCRIPT_VISION.md`, the three-figure revision plan, the
[prospective plan](../../PLAN.md), the
[accepted linear pilot](../../../all_five_recovery_pilot_20261017/RESULTS.md),
and the retained ring producer and moment-extraction sources. This review ran
no integrations, changed no native equations or retained data, and generated
no figure variant.

**Verdict: accept the bounded 15-state ring port preflight.** The plan addresses
the intended scientific questions without replacing the all-five linear figure
or the separate replicated-noise figure. The 702-state recovery grid is
conditional on actual algebraic and clock results; this prospective acceptance
does not certify those results. No broad campaign redesign is needed.

**Actual port addendum.** The subsequently completed 15-state pilot is accepted
at the reviewed clocks, clearing the clock condition for the full recovery
grid. The [algebraic record](RECOVERY_ALGEBRAIC_PREFLIGHT.json) matches the
producer's even-reflection, core/taper/kick and periodic-mode fixtures; maximum
mirrored-RHS disagreement is 2.06 × 10⁻¹⁵. This review independently verified all
32 frozen-source hashes and all 702 copied case-input hashes. The separate
[input audit](RECOVERY_INPUT_AUDIT.json) records 652 exact t = 10 complete-state
comparisons and all prescribed nonnegative kicks.

Independently recomputing the gates from all 15 saved pilot NPZs reproduces all
45 amplitude responses: maximum prefix E is 2.85962 × 10⁻⁷, maximum |F| + E is
0.19121824, and maximum amplitude-linearity excess is −1.34607 × 10⁻⁹. All 15
archive hashes and the pilot CSV hash verify. Saved initial kicks exactly
match the frozen basal direction, compensation and untouched other fields;
saved sham/zero initial and final native fields are identical. The source also
checks complete sham/zero histories during execution. There is no remaining
essential recovery-grid launch blocker; the full grid still requires its own
reported statewise gates and completion audit.

## Conditions before the full recovery grid

1. **Port the intervention, geometry and cue continuation together.** Preserve
   the original basal background, excess profile, two half-contrast core radii,
   compact cosine taper and fractional core weights when probing a challenged
   snapshot. At the same complete t = 10 state, sham, independent zero kick and
   each active-to-inactive transfer must evolve with the identical retained
   cue. All remaining native fields evolve. Verify positive initial core
   deficit, unchanged other fields, conservative corresponding pools, periodic
   seam handling, and zero-kick identity. An unavailable nonnegative kick is a
   declared unavailable amplitude/state; clipping, substituting a current-state
   taper or changing damage by separation would change the assay.
2. **Use the actual fixed-clock prefix gates.** Check WP/Holmes at
   0.1778279410038923 and Goryachev at 0.005623413251903491 native time, across
   all three positive amplitudes and every sampled prefix. Reproduce the
   reviewed empirical sensitivity, pairwise 5% amplitude agreement with the
   3(Ea + Eb) allowance, and |F| + E ≤ 0.2 rule. Tighter integration of the same
   mesh supplies tolerance sensitivity, not spatial convergence. The 15-state
   selection covers three basal sizes and two challenged extremes per model;
   it is an efficient portability check rather than a proof for intermediate
   states. Apply and export these gates throughout the full grid too. A failed
   state must not receive a qualified coefficient through a private clock or
   omitted sensitivity amplitude. Coordinate any revised clock once per model
   with the linear analysis before making a common-clock comparison.
3. **Separate reflection identity from independent-state agreement.** On a
   matching 2N-cell ring, even reflection must reproduce the native Neumann
   reaction/diffusion operator, mirrored taper and integrated core weights.
   Whole-incumbent width and amount double, whereas the local coefficient and
   peak agree. The actual independent basals retain their distinct native
   meshes and preparation histories; quantify their observed differences
   rather than describing the algebraic construction as their numerical
   validation. The existing Goryachev linear L39 core has four sample centers,
   so tighter time integration cannot establish a continuum size effect.
4. **Keep all planned cells and failure reasons.** The 867-cell inventory is
   652 complete challenges, 198 overlap exclusions and 17 unavailable Holmes
   C80 cases. Fifty qualified basals plus those 652 snapshots give 702 probes,
   not 867 valid recovery values. Distinguish an unavailable preparation,
   excluded cue geometry, unavailable kick, execution failure and failed
   numerical/early-response gate. Missing execution is not a biological
   outcome. Record measured storage and source/input hashes at completion.

The common readout remains a signed local-core deficit response:

\[
D(t)=A_{\mathrm{sham,core}}(t)-A_{\mathrm{probe,core}}(t),
\qquad
K_{\tau}=\frac{1-D(\tau)/D(0)}{\tau}.
\tag{1}
\]

A positive value can coexist with loss of actual core activity under challenge.
Export and interpret sham/probe changes and rates, as well as whole-ring and
outside-core deficits. Recovery includes reaction and transport and does not
isolate intrinsic stability or establish a causal role for incumbent width.

## Finite-cue and vector evidence

The nine declared anchors adequately cover the illustrated response classes
and all three models for a bounded duration check. Eighteen new responses at
500 and 2000 cue units, reusing 1000, can show whether those selected outcomes
depend on duration; they cannot certify every grid cell as duration-independent.
Likewise, the complete-grid 3000/8000/13000 cue-free readouts assess dependence
on the declared observation window, not convergence to an asymptotic endpoint.
Preserve duration-sensitive outcomes instead of changing the protocol to make
the existing map invariant. State peak, support, duration and recovery window
next to each comparison. These are standardized native protocols within each
model, without a cross-model biochemical or time calibration.

Use exact complete restart fields only: primary-only intermediate movies do
not supply hidden inactive/Rho fields at 500. A 2000-unit cue may continue the
exact complete retained 1000-unit state. Check phase continuity, positivity and
native pools; compare any duplicated primary prefixes on their exact sampled
times. All recovery readout times are measured from the corresponding cue
withdrawal, including the duration controls.

**Resolved before duration launch.** The initial duration source compared a
replayed T = 500 prefix only at complete native states 0, 1 and 10, leaving its
duplicated primary history near withdrawal unchecked. The corrected source
compares the primary field at all exact shared samples through 500, including
500, and records sample count, final time and raw/scaled disagreement. It
retains the 0/1/10 complete-state checks. This review verified the corrected
source, all 26 frozen duration-source hashes and its manifest provenance. The
complete t = 1000 cue continuation and T = 500 restart selection are correct.
The compact source supplies primary activity, rather than hidden native fields,
at 500. No remaining essential design blocker prevents the declared bounded
duration execution; the actual prefix, positivity and pool results remain to
be checked.

The read-only global moment plan is appropriate: retain P = M1(t)/M1(0), its
magnitude, the masked angle, total activity and mass-normalized first and
second moments. The extractor checks exact requested readout availability,
phase continuity, source hashes and uniform-background cancellation. Minimum
magnitude and angle excursions are **sampled** summaries. Low first moment can
mean opposed active sites, so the Holmes cancellation example and native local
histories are essential. These deterministic vectors supply activity-pattern
evidence; Figure 3 alone owns the replicated noisy persistence analysis.

## Assembly safeguards and scientific message

The existing `ring_figure2_evidence_design_20261017.py` exports
`linear_conditional_repair` and `ring_basal_conditional_repair` from the earlier
inactive-withdrawal/return tables. Those columns cannot populate the new Figure
2C baseline or heatmap. Use the new common active-deactivation evidence and
retain the older assay only as explicitly identified support. Basal recovery
has no future cue-location variable; a separation map must say **t = 10,
cue retained**, and preferably display the change from its own basal value.

The baseline gain correspondence must compare the same local observable.
Point-local gain and finite-width cue-weighted gain can both be reported but
should not be treated as identical measurements. A linear endpoint support of
3 reflects to full ring support 6. The retained ring challenge has full support
3, so the C = 2L basal construction supplies no finite-response equivalence.
Shared scales must retain genuinely flat separation dependence and permit
negative recovery coefficients if observed.

The useful outcome is model-specific coupling, or dissociation, between
incumbent organization, initial challenger gain, conditional local recovery and
the later spatial competition. The Holmes C40 d12/d16 comparison relocates the
cue in the identical complete basal; C30 d12 versus C40 d12 also changes total
material, architecture and the complementary route. Endpoint location labels
do not prove site lineage; merger, return and replacement claims need their
highlighted histories. Experimental proposals can pair local activity-loss
responses and challenger establishment with global activity-polarity magnitude
and reorientation, plus local protrusion histories. No cell velocity law or
inference of measured migration persistence is supplied here.

**Anticipated user critique (hypothetical, not a quotation):** the next
checkpoint should put gain/recovery and the declared long-cue outcomes together
in a readable ring argument, show the common-assay baseline rather than the old
repair map, and preserve the full model-dependent and nonmonotone evidence.
Another transient-only or example-only sheet would leave the stated objective
unfinished. The current bounded completion plan can meet that objective once
its actual gates and outcomes have been audited.
