> Active manuscript work moved to [../manuscript/](../manuscript/AGENTS.md). Do not extend this deprecated suite or promote its old reversal figures. The notes below preserve historical rules and frozen keeper protections.

# Active normalized two-assay suite

Read README.md and results/status.json before interpreting or modifying results.
The user authorized normalized cues at all existing lengths, and revision 11
uses fully settled conditioning -> release profiles and persistent 1% arrival clocks. Competition/periodic tracking
remain deferred. The figure shows shared conditioning once, separate responses,
and no cumulative fits.

- Gradient reversal and finite countercue share the same linear cue shape and
  amplitude within each model. Preserve this normalization and retained native
  initial states, eligibility/response gates, horizons and duration grids. Settling
  uses the user-authorized continuation/holdout rules and persistent 1% arrivals.
  Every model conditions directly from its native initial state. Do not add
  spontaneous pre-settling or reflection to either active assay.
- `baseline/` remains the immutable accepted v8 reference. Validate paired
  readout preservation and the conditioning -> release checkpoint graph. All
  revision-11 response inputs intentionally change; do not demand exact old readouts.
- Competition and periodic tracking are inactive. Do not restore either unless
  requested. Spring competition additionally remains paused for lack of a clear
  length contrast; do not retry, refine or rescue it.
- The archived publication_figures/assay_explanation_v2 is an approved KEEPER;
  preserve its sources, reconstructed fields and provenance.
- Runtime must not import reports, historical archives or prior revisions.
- Use the existing environment in environment.json and SLURM for full reruns.
- Do not modify executable inputs while a frozen campaign is running.
- Numerical completion does not establish the length-response interpretation.
