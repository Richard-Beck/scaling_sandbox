# Revision 8: bounded protocol standardizations from revision 6

Authorized 2026-10-01: rerun the four-test suite and generate the suite and cue
settling PDFs, using consistent recovery, settling confirmation and response
sampling, plus the tested spring cue/countercue preparation. Revision 6 remains
preserved as the baseline; rejected revision 7 is not an input.

## Changes applied

1. All existing steady preparation phases use a 1000-native-unit confirmation
   window, 100-unit detection intervals, and the existing 10000-unit cap. The
   0.1% pointwise change/projected drift test covers every evolving species and
   applicable mechanics. The polarity eligibility window is also 1000 units.
   deBelly retains its additional Rac/MCA half-mass position checks.
2. All finite countercues recover without stimulus for 10000 native units.
   The existing final 500-unit response success window and thresholds remain.
3. Response samples: 2.5 native units for competition, gradient reversal and
   countercue recovery; finite-pulse samples min(2.5, duration/20), providing at
   least 20 intervals for brief pulses; periodic samples P/80 for every model.
   Model-specific solvers and tolerances are retained. Endpoint time resolution
   is now 2.5 units. Preparations keep 5-unit output and 100-unit settling checks.
4. Spring finite countercue now directly applies a left-edge Gaussian cue at
   existing amplitude 4, width parameter 0.2, to its default initial state until
   full-state settling. It then settles without cue before an independently
   restarted right Gaussian challenge. No spontaneous preconditioning or
   reflection in this assay. This extends the successful L0=1 pilot to all lengths.

Spring's projected cue drift uses the solver's exact nonlinear activation term;
revision 6 used cue*v as an approximation. This correction makes the full-state
settling check consistent with the equations actually integrated, without
changing model dynamics or stimulus strength.

## Explicitly retained

All stimulus strengths and shapes from revision 6; original model coefficients
and initial states; revision-5 lengths and each model's dense duration grid;
competition's existing widths and model-specific taper/abrupt withdrawal;
20-unit gradient crossfade and 20000-unit competition/reversal observation;
response readouts, thresholds and dwell times; native model units. Holmes keeps
localized finite training and its stronger linear challenge. Periodic starting
states remain the revision-3 preparations (including Holmes fixed training and
spring fixed spontaneous preparation), and P=100/400, six cycles and scoring
P through 6P remain unchanged. The longer confirmation window applies to
existing steady-preparation phases, not to legacy periodic initialization.
Spring gradient preparation retains spontaneous settling/reflection before
left-gradient settling. Spring competition remains explicitly paused.

Remaining differences are cue geometry/strength and coupling, native units,
model initial states/coefficient regimes, periodic initialization, competition
withdrawal, and checks needed for extra species/mechanics. No additional
harmonization of these protocols is implied.

## Execution and outputs

Existing scaling_sandbox Conda environment. SLURM QoS small, one core and 3 GB
per array slot; no concurrency throttle. Sources, model code, this protocol,
manifest and duration grids are frozen and hash-checked. Report modules load
explicitly to avoid legacy import shadowing. No revision-6 endpoints are reused.

`launch.json` records jobs and execution stage. `results/status.json` records
completeness; missing execution records are distinct from scientific censoring.
Figures update as data arrive, and at final collection. Each PNG/PDF shows its
generation timestamp. HTML uses manual refresh only.

Requested PDFs: `results/four_test_all_models.pdf`, `results/time_to_settle.pdf`.
PNG/SVG previews, endpoints, settling diagnostics and provenance are also saved.
Completed computation does not establish scientific interpretation or approval.
