# Active normalized two-assay suite

Read README.md and results/status.json before interpreting or modifying results.
The user authorized normalized cues at all existing lengths, and revision 10
removes Spring gradient pre-settling/reflection. Competition/periodic tracking
remain deferred. The figure shows shared conditioning once, separate responses,
and no cumulative fits.

- Gradient reversal and finite countercue share the same linear cue shape and
  amplitude within each model. Preserve this normalization and retained native
  initial states, timing, settling/eligibility gates, horizons and duration grids.
  Every model conditions directly from its native initial state. Do not add
  spontaneous pre-settling or reflection to either active assay.
- `baseline/` remains the immutable accepted v8 reference. Validate paired
  readout preservation and exact unchanged-protocol controls, not numerical
  equality of intentionally changed finite-cue experiments.
- Competition and periodic tracking are inactive. Do not restore either unless
  requested. Spring competition additionally remains paused for lack of a clear
  length contrast; do not retry, refine or rescue it.
- The archived publication_figures/assay_explanation_v2 is an approved KEEPER;
  preserve its sources, reconstructed fields and provenance.
- Runtime must not import reports, historical archives or prior revisions.
- Use the existing environment in environment.json and SLURM for full reruns.
- Do not modify executable inputs while a frozen campaign is running.
- Numerical completion does not establish the length-response interpretation.
