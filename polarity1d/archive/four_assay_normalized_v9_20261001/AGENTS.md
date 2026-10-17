# Active normalized two-assay suite

Read README.md and results/status.json before interpreting or modifying results.
The user authorized revision 9 on 2026-10-01: apply the successful normalized-cue
pilot at all existing sampled lengths; defer competition and periodic tracking.

- Gradient reversal and finite countercue share the same linear cue shape and
  amplitude within each model. Preserve this normalization and retained native
  preparation, timing, settling/eligibility gates, horizons and duration grids.
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
