# Report work: current status

Updated 2026-09-29 from the user's status clarification. This note supersedes
earlier claims that the four-test comparison is scientifically complete or final.

## Complete for now

- **Assay explanation:** [version 2 PDF](four_test_suite/publication_figures/assay_explanation_v2.pdf), [PNG](four_test_suite/publication_figures/assay_explanation_v2.png), and [notes/provenance](four_test_suite/publication_figures/ASSAY_EXPLANATION_V2_NOTES.md). This is the approved keeper; the original explanation draft was rejected.
- **Figure 3:** [canonical integrated PNG](figure_revision/figures/figure3.png), [caption](figure_revision/figures/figure3_review/figure3_updated_caption.md), and [handoff](figure_revision/figures/figure3_review/FINAL_HANDOFF.md). The existing Figure 3 PDF is stale; use the PNG.
- **Hosny sensitivity supplement:** [PNG](figure_revision/figures/figure3_review/hosny_sensitivity/hosny_sensitivity.png) and [caption](figure_revision/figures/figure3_review/hosny_sensitivity/CAPTION.md).

Preserve these completed outputs, sources and provenance. No report PDF has been
rebuilt to integrate the current standalone figures.

## Four-assay suite: ongoing diagnosis and consolidation

The current seven-model comparison is [revision 5 PDF](four_test_suite/revision5/results/four_test_all_models.pdf)
([PNG](four_test_suite/revision5/results/four_test_all_models.png)). Its current
individual panels are [competition](four_test_suite/revision5/results/competition.pdf),
[gradient reversal](four_test_suite/revision5/results/gradient_reversal.pdf),
[finite countercue](four_test_suite/revision5/results/finite_countercue.pdf), and
[periodic amplitude/phase](four_test_suite/revision5/results/periodic_harmonics.pdf).
See [revision-5 protocol notes](four_test_suite/revision5/RESULTS_README.md) and
[task handoff](four_test_suite/TASK_HANDOFF.md) for numerical coverage and history.

Revision 5 is a completed numerical snapshot, **not a final consolidated suite**.
Diagnosing and interpreting its length responses remains active work. The spring
competition investigation triggered a cross-panel audit of whether preparations
are steady when the initial cue is released or changed. The current audit covers
competition and gradient reversal across all seven models; it does not certify
every preparation in all four assays.

The [preparation audit](four_test_suite/preparation_audit/README.md) now reports
308 cases: 117 steady, 173 not steady and flagged for later, and 18 unverified
under the spring pause. See the [model/test summary](four_test_suite/preparation_audit/results/REPORT.md)
and [rerun flags](four_test_suite/preparation_audit/results/flagged_for_later_rerun.csv).
These findings qualify the existing length curves; many experiments are expected
to rerun following audit review. This status note does not launch reruns.

The [spring steady-cue sweep](four_test_suite/publication_figures/spring_competition/steady_cue_full/README.md)
and its [diagnostic figure](four_test_suite/publication_figures/spring_competition/steady_cue_full/results/spring_competition_steady_full.pdf)
are retained as diagnostic evidence, not a validated length-scaling result.
Spring competition remains paused with **NO CLEAR LENGTH CONTRAST**; follow
[the suite's user decisions](four_test_suite/AGENTS.md) before any further runs.

## Outstanding scientific work

First review the preparation audit, rerun the affected experiments as authorized,
and consolidate a stable four-assay figure suite. Then explore, explain and audit
the length-dependent effects **or lack of effects**, and communicate those findings
with well-designed figures. Do not treat successful execution, endpoint coverage
or numerical validation alone as completion of that scientific interpretation.

## Superseding campaign guidance — 2026-10-01

Revision 6 completed all preparations and response jobs on September 29, with
both PNGs generated; see the final section of four_test_suite/TASK_HANDOFF.md.
The user has now requested standardized modal cue/preparation/timing protocols
and regeneration of both figures. This work is in four_test_suite/revision7/;
read its README.md and check launch.json/results/status.json for current state.
The earlier statements above identifying revision 5 as current and reruns as
still awaiting authorization are historical. Revision-6 outputs are preserved.
Spring competition remains paused; scientific review remains outstanding.

## User decision: revision 6 restored as baseline — 2026-10-01

The user rejected revision 7 because numerical stimulus strengths were changed
without being requested. Revision 6 results are the baseline. Revision 7 is
superseded and must not be used as the accepted comparison, even though its
producer completed. The user now requests discussion of possible low-impact
standardizations only; no corrected campaign or additional reruns are authorized.
Keep revision-6 outputs and selected stimulus strengths intact. Spring
competition remains paused. This decision supersedes earlier revision-7 guidance.

## Revision 8 authorized — 2026-10-01

The user subsequently authorized a full four-test rerun from revision 6 with
three timing/sampling standardizations and the successful spring cue/countercue
preparation. See revision8/README.md: shared 1000-unit steady confirmation,
10000-unit countercue recovery, common response sampling, spring left Gaussian
conditioning then autonomous settling. All existing stimulus strengths/shapes,
other preparation histories and periodic starting states are retained. Requested
outputs are the suite and settling-time PDFs, with PNG previews. Revision 7
remains rejected; revision 6 is preserved as baseline. Spring competition stays
paused. Check revision8/launch.json and results/status.json for current progress.

## Revision 8 completed — 2026-10-01

All 441 preparations and 4096 response jobs completed with exit code 0; both
requested PDFs and PNG previews are in four_test_suite/revision8/results/
(relative to reports; omit four_test_suite/ when reading the suite handoff).
Final collection has complete duration coverage and no execution errors. Only
one category differs from revision 6: Otsuji competition at 25 micrometres is
now unsettled. Spring countercue succeeds at all 21 lengths with unchanged
thresholds. See revision8/FINAL_HANDOFF.md and comparison_to_revision6.json
for details. Revision 6 remains the preserved baseline, revision 7 is rejected,
and scientific review remains outstanding.

## Revision-8 visualization update — 2026-10-01

User requested mode-colored outcome ranges, length regressions, a merged plot
and model-by-metric slope heatmap. Nineteen Goryachev competition successes
(18 at zero, one within the first 2.5-unit sample) are classified as censored
with mode resolved_at_withdrawal; original statuses/endpoints remain in raw
records, original CSV columns and results/before_mode_revision/. This is not
conventional right-censoring and does not indicate failure to resolve. Updated
exports have 498 success, 161 censored, 48 unreliable, 35 ineligible, 7 unsettled
and 21 paused cases. All 66 eligible series have OLS lines (>3 valid points);
periods and periodic outputs are fitted separately. Combined plot:
revision8/results/length_response_and_settling.png; slope heatmap:
revision8/results/length_response_gradients.png. Sources and provenance are in
revision8_visualization/render.py and results/visualization_provenance.json.
Use this renderer for current exports; frozen build_report.py restores the
original classification. No simulations were rerun or modified.
