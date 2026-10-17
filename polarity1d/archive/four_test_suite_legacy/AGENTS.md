# Four-test suite: current user decisions

Read `TASK_HANDOFF.md` before interpreting results or starting work here.

- **Spring competition investigation is PAUSED by the user (2026-09-29).** Its
  current conclusion is **NO CLEAR LENGTH CONTRAST**. Do not resume, extend,
  refine, retry, smooth, or attempt to rescue the spring competition workflow,
  including the steady-cue sweep, tolerance studies, feedback interventions,
  or an explanation of its old U-shaped curve. Resume only if the user explicitly
  requests it. Preserve existing data as diagnostics; do not present them as a
  validated length effect. This restriction covers workflow entry points in
  `publication_figures/spring_competition/` and equivalent new implementations.
- The user authorized the revision-6 steady-state preparation and endpoint
  campaign on 2026-09-29, including BOTH gradient reversal and finite countercue.
  See `revision6/README.md` and `revision6/launch.json`. The user explicitly kept
  spring competition paused and periodic preparation unchanged. All other
  authorized assays may be rerun from qualifying full-state checkpoints.
  Preserve frozen executable inputs during this campaign; diagnose execution
  errors separately from scientific failure to settle within the 10,000-unit cap.
- `publication_figures/assay_explanation_v2.{pdf,png,svg}` is a user-approved
  **KEEPER**. Preserve its sources, reconstructed fields and provenance. Do not
  replace or redesign it without a user request. The original `assay_explanation.*`
  draft was rejected. No flowcharts; minimize explanation inside plots.
- The prior audit flags remain historical diagnostics. New reruns are authorized
  by the explicit revision-6 request, not by those flags. Spring competition
  remains excluded.

- For revision-6 report regeneration use
  `revision6_report_service/service.py report`. Legacy adapters alter sys.path;
  a bare `import build_report` in a new driver can import revision 3. The active
  progress/final figure jobs use explicit revision-6 module loading. See
  `revision6/report_orchestration_repair.json`; scientific inputs remain frozen.

- On 2026-10-01 the user authorized the revision-7 standardized modal protocol
  campaign. This supersedes the revision-6 periodic-preparation exception.
  Read revision7/README.md; preserve its frozen inputs while running. Spring
  competition remains paused. Retain revision-6 artifacts for comparison.

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
