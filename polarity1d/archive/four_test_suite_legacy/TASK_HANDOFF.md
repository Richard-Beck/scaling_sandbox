# Seven-model comparison: task handoff

Revision-5 output generation completed 2026-09-29. **The suite remains provisional:** length-response diagnosis, preparation-audit review and later reruns are ongoing work. Read [current project status](../CURRENT_STATUS.md) and [the preparation audit](preparation_audit/README.md). Revision 5 is the current numerical snapshot; earlier revisions are retained for comparison.

## Current user decisions — 2026-09-29

- **KEEPER:** [assay explanation v2](publication_figures/assay_explanation_v2.pdf)
  ([PNG](publication_figures/assay_explanation_v2.png)). The user approved this
  cue-kymograph/time-coloured-profile figure. Preserve its sources and provenance.
  The first flowchart-based assay draft was rejected.
- **PAUSED: spring competition — NO CLEAR LENGTH CONTRAST.** Do not resume,
  extend, refine, retry or attempt to rescue the workflow or explain its old U-shape
  without an explicit new user request. Retain original and steady-cue results as
  diagnostic evidence, not a validated length effect. See
  [workflow status](publication_figures/spring_competition/WORKFLOW_STATUS.json)
  and [agent instructions](AGENTS.md).
- **Initial-cue audit completed:** [14-combination report](preparation_audit/results/REPORT.md),
  [all length diagnostics](preparation_audit/results/case_diagnostics.csv), and
  [cases flagged for later reruns](preparation_audit/results/flagged_for_later_rerun.csv).
  Across 308 original preparation cases, 117 pass the operational stationarity
  checks, 173 do not, and 18 paused spring-competition cases lack archived full
  preparation fields and remain unverified. The user confirmed flags should apply
  to preparations **not** at steady state. No corrected endpoint reruns were
  launched, and spring competition was not simulated further. See
  [criteria and scope](preparation_audit/README.md). Rerun flags are not authorization
  to start campaigns; the spring pause still applies.

## Where to find the results

All paths below are relative to `reports/four_test_suite/` in `~/projects/scaling_sandbox`.

| Output | Location |
|---|---|
| Main seven-model comparison | [PDF](revision5/results/four_test_all_models.pdf) · [PNG](revision5/results/four_test_all_models.png) |
| Finite-countercue figure | [PDF](revision5/results/finite_countercue.pdf) · [PNG](revision5/results/finite_countercue.png) |
| Periodic amplitude and phase | [PDF](revision5/results/periodic_harmonics.pdf) |
| Numerical results | [Endpoints](revision5/results/endpoints.csv) · [all countercue trials](revision5/results/finite_duration_tests.csv) · [harmonic fits](revision5/results/periodic_harmonics.csv) |
| Protocol details and qualifications | [Revision 5 notes](revision5/RESULTS_README.md) |

The comparison covers wave-pinning, Otsuji, Goryachev, LEGI, deBelly, Holmes Model 3, and spring. Its four assays are peak competition, gradient reversal, finite countercue, and periodic direction following; the periodic assay occupies two columns (amplitude and phase).

Case JSONs are in [revision5/data/](revision5/data/), and new duration-chunk results in [revision5/trial_data/](revision5/trial_data/). The [case manifest](revision5/manifest.json), [trial manifest](revision5/trial_manifest.json), [collection map](revision5/collection_manifest.json), and [frozen provenance](revision5/frozen_provenance.json) identify source records, reused results, and executable hashes. Some diagnostics remain in earlier revisions and are referenced by those records.

## What changed, and why

1. **Initial screen:** built the seven-model/four-test comparison on sparse length and duration grids. Poor endpoint coverage prompted a model/protocol postmortem.
2. **Revision 2 — model and stimulation repairs:** restored deBelly's established reflection-symmetric baseline and activation defaults. Spring was calibrated against its standalone positive control using native drive amplitude 4. Holmes required correcting missing factors of two in the inhibition kinetics, then selecting Rac activation 3.5 from a bounded polarization search; this is a selected regime, not a published default. Regimes were fixed across lengths. Search records: [deBelly](revision2/DEBELLY_SEARCH.md), [spring](revision2/SPRING_SEARCH.md), [Holmes](revision2/HOLMES_SEARCH.md).
3. **Revision 3 — undo periodic measurement drift:** the initial comparison had substituted per-switch latency without agreement. Restored the previous edge-orientation harmonic amplitude and signed phase lag, excluding the first cycle (fit over P through 6P). The accepted gradient-reversal protocol was retained.
4. **Revision 4 — restore length resolution:** replaced sparse sampling with 1–39 µm every 2 µm for the original four models and deBelly; Holmes uses 10–50 µm every 2 µm, and spring L0=0.6–2.6 every 0.1 native unit. Existing off-grid points were retained. This changed sampled domain lengths, not the 201-cell simulation mesh.
5. **Revision 5 — restore duration resolution:** deBelly's nine-point countercue grid was causing artificial steps. Restored the historical 50-point sweep from 10–1,000 s, retaining old points (58 total). Refined the other grids to 291 candidates for each original model, 341 for Holmes, and 111 for spring. Added 31,802 simulations without changing cues or scoring. deBelly thresholds changed from 250 to 232.245 s at 21 µm and from 800 to 636.327 s at 39 µm, matching the earlier Rac thresholds at those lengths.

## Validation and remaining qualifications

Review caught a Holmes dispatch error that would have used the old adapter; corrected routing and representative old/new-duration smoke tests passed before production. SLURM campaigns used frozen sources and per-job hash checks. Final audits verified all 770 cases, 616 unchanged non-countercue records, 1,565 unchanged old trials, and all 33,367 exported countercue trials. No execution failures remained; 21 preparation-ineligible cases explicitly skipped new duration tests.

Thresholds are the shortest **tested** passing durations. Seven wave-pinning/Holmes cases have nonmonotonic pass/fail outcomes near the settling boundary; these are documented and retained without smoothing. Historical deBelly Figure 2 also used different preparation timing and x50 criteria, so finer sampling does not imply complete protocol equivalence. Model-native time units and spring's uncalibrated length units remain explicit.


## Revision 6 steady-state rerun campaign — authorized 2026-09-29

The user requested and authorized settling-time measurements, full-suite reruns from qualifying states, a new suite figure and a settling companion. Clarifications: **keep spring competition paused**; **keep existing periodic preparation**. Both gradient reversal and finite countercue are included independently. See [campaign methods](revision6/README.md), [SLURM launch/stage state](revision6/launch.json), and [generated results/status](revision6/results/README.md). Revision-6 figures are preliminary until the producer and all duration chunks have finished; do not substitute revision-5 endpoints into missing cells. Each preparation phase has a 10,000-native-unit cap and a full-state settling test. Finite countercue checks its cue-conditioned and autonomous preparation phases. Failed preparations disqualify only the affected model/test/length; execution errors are separate. Exact full mechanical-state deBelly restarts were validated. The full dense revision-5 duration grids and existing response scoring are retained. The assay explanation v2 remains a keeper and was not modified.

Launch confirmed: preparation array **20642503** finished all 441 cases (401 qualified; 34 steady but ineligible; 6 unsettled). Response array **20642946** contains 4097 jobs, including 3523 finite duration chunks; progress job **20644318** and final figure job **20644319** use the corrected explicit report service. Frozen simulation inputs were not changed. Six unsettled cases: deBelly finite-countercue no-cue preparation at L=11,12,13,15,17,18. The settling companion now contains all preparation measurements; the suite figure fills as complete responses/sweeps arrive. Check launch.json and results/status.json for current completion state.

Final collection completed **2026-09-29 21:37 UTC**: all 441 preparation jobs and 4097 response jobs completed with exit code 0, full duration coverage verified, no missing outputs or execution errors. `revision6/launch.json` stage is `finished`; `results/status.json` has `complete: true`. Both new figures are complete, with 518 successful assay records, 142 censored, 48 unreliable periodic fits, 35 ineligible, 6 unsettled and 21 paused spring-competition records. The extra ineligible case beyond the 34 preparation exclusions lost the required two peaks during the retained competition withdrawal taper.

## Revision 7 protocol standardization — authorized 2026-10-01

The user requested a modal ruleset across models for cue shape, application
history, conditioning, settling and observation times; disagreements in times
use the longest existing value. The request supersedes the revision-6 decision
that periodic preparation remain unchanged. New protocols and frozen inputs are
in `revision7/README.md`; revision-6 outputs remain preserved. Spring competition
remains paused. Revision 7 uses identical numerical amplitudes and relative cue
geometry, common left-gradient conditioning (including spring and periodic),
1000-unit settling windows, 10000-unit countercue recovery, and the union of all
existing dense countercue duration grids. See revision7/launch.json and
revision7/results/status.json for execution and completeness, not stale earlier
status paragraphs. Native model units and activation coupling remain distinct.

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
